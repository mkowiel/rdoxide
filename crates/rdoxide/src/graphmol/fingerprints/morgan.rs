//! Morgan fingerprint generator (upstream `GraphMol/Fingerprints/MorganGenerator.{h,cpp}`
//! and `MorganFingerprints::getConnectivityInvariants` from `FingerprintUtil.cpp`).

use std::cmp::Ordering;
use std::collections::HashSet;

use super::{AdditionalOutput, FingerprintArguments, hash_combine, hash_pair, hash_vect};
use crate::error::{RDKitError, Result};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::chirality::get_use_legacy_stereo_perception;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// `MorganArguments` (the Morgan-specific part).
#[derive(Clone, Debug)]
pub struct MorganArguments {
    pub radius: u32,
    pub only_nonzero_invariants: bool,
    pub include_redundant_environments: bool,
}

impl MorganArguments {
    pub fn info_string(&self) -> String {
        format!(
            "MorganArguments onlyNonzeroInvariants={} radius={}",
            u8::from(self.only_nonzero_invariants),
            self.radius
        )
    }
}

/// `MorganAtomEnv`
#[derive(Clone, Debug)]
pub struct MorganAtomEnv {
    pub code: u64,
    pub atom_id: u32,
    pub layer: u32,
}

impl MorganAtomEnv {
    pub(super) fn update_additional_output(
        &self,
        mol: &ROMol,
        ao: &mut AdditionalOutput,
        bit_id: u64,
    ) {
        if let Some(m) = &mut ao.bit_info_map {
            m.entry(bit_id)
                .or_default()
                .push((self.atom_id, self.layer));
        }
        if let Some(c) = &mut ao.atom_counts {
            c[self.atom_id as usize] += 1;
        }
        if let Some(a) = &mut ao.atom_to_bits {
            a[self.atom_id as usize].push(bit_id);
        }
        if let Some(apb) = &mut ao.atoms_per_bit {
            let mut involved = vec![self.atom_id as i32];
            if self.layer > 0 {
                let n = mol.num_atoms() as usize;
                let dm = crate::graphmol::matrices::get_distance_mat(mol);
                for i in 0..n {
                    let d = dm[self.atom_id as usize * n + i] + 0.1;
                    if (d as u32) <= self.layer && i as u32 != self.atom_id {
                        involved.push(i as i32);
                    }
                }
            }
            apb.entry(bit_id).or_default().push(involved);
        }
    }
}

/// `MorganFingerprints::getConnectivityInvariants`
pub fn get_connectivity_invariants(mol: &ROMol, include_ring_membership: bool) -> Result<Vec<u32>> {
    let pt = PeriodicTable::get();
    let mut invars = Vec::with_capacity(mol.num_atoms() as usize);
    for i in 0..mol.num_atoms() {
        let atom = mol.atom(i);
        let mut comps = vec![
            atom.atomic_num(),
            mol.total_degree(i)?,
            mol.total_num_hs(i, true)?,
            atom.formal_charge() as u32,
        ];
        let delta_mass = (atom.mass() - pt.get_atomic_weight(atom.atomic_num())?) as i32;
        comps.push(delta_mass as u32);
        if include_ring_membership && mol.ring_info().num_atom_rings(i)? != 0 {
            comps.push(1);
        }
        invars.push(hash_vect(&comps));
    }
    Ok(invars)
}

/// `MorganBondInvGenerator::getBondInvariants` (with the new stereo
/// perception the molecule needs its CIP labels)
pub fn morgan_bond_invariants(mol: &ROMol, use_bond_types: bool, use_chirality: bool) -> Vec<u32> {
    (0..mol.num_bonds())
        .map(|i| {
            let bond = mol.bond(i);
            let mut inv: i32 = 1;
            if use_bond_types {
                if !use_chirality
                    || bond.bond_type() != BondType::Double
                    || bond.stereo() == BondStereo::StereoNone
                {
                    inv = bond.bond_type() as i32;
                } else {
                    let mut stereo = bond.stereo();
                    // with the new stereo perception the E/Z labels come from
                    // the CIP labels; other stereo is used as it is
                    if !get_use_legacy_stereo_perception() {
                        match bond.get_prop(cp::_CIPCode).map(|v| v.to_rd_string()) {
                            Some(c) if c == "E" => stereo = BondStereo::StereoE,
                            Some(c) if c == "Z" => stereo = BondStereo::StereoZ,
                            _ => {}
                        }
                    }
                    inv = 100 + 10 * bond.bond_type() as i32 + stereo as i32;
                }
            }
            inv as u32
        })
        .collect()
}

/// A `boost::dynamic_bitset<>` over the bonds.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct BondSet(Vec<u64>);

impl BondSet {
    fn new(n: u32) -> Self {
        BondSet(vec![0; n.div_ceil(64) as usize])
    }
    fn set(&mut self, i: u32) {
        self.0[(i / 64) as usize] |= 1 << (i % 64);
    }
    fn or_assign(&mut self, o: &Self) {
        for (a, b) in self.0.iter_mut().zip(&o.0) {
            *a |= b;
        }
    }
    /// `operator<` of equally-sized dynamic_bitsets: compares from the highest block
    fn cmp_bits(&self, o: &Self) -> Ordering {
        if self.0.is_empty() {
            return Ordering::Equal;
        }
        for (a, b) in self.0.iter().rev().zip(o.0.iter().rev()) {
            match a.cmp(b) {
                Ordering::Equal => continue,
                other => return other,
            }
        }
        Ordering::Equal
    }
}

/// `MorganEnvGenerator::getEnvironments`
pub fn get_environments(
    mol: &ROMol,
    margs: &MorganArguments,
    args: &FingerprintArguments,
    from_atoms: Option<&[u32]>,
    atom_invariants: &[u32],
    bond_invariants: &[u32],
) -> Result<Vec<MorganAtomEnv>> {
    let n_atoms = mol.num_atoms();
    let n_bonds = mol.num_bonds();
    if atom_invariants.len() < n_atoms as usize {
        return Err(RDKitError::precondition("bad atom invariants size"));
    }
    if bond_invariants.len() < n_bonds as usize {
        return Err(RDKitError::precondition("bad bond invariants size"));
    }
    let mut result = Vec::with_capacity(((margs.radius + 1) * n_atoms) as usize);

    // OutputType values; Morgan codes always fit in 32 bits
    let mut current: Vec<u64> = atom_invariants.iter().map(|&v| u64::from(v)).collect();
    let mut next_layer = vec![0u64; n_atoms as usize];

    let mut include_atoms = vec![from_atoms.is_none(); n_atoms as usize];
    if let Some(fa) = from_atoms {
        for &i in fa {
            *include_atoms
                .get_mut(i as usize)
                .ok_or_else(|| RDKitError::precondition("bad fromAtoms index"))? = true;
        }
    }
    let mut chiral_atoms = vec![false; n_atoms as usize];
    let mut neighborhoods: HashSet<BondSet> = HashSet::new();
    let mut atom_neighborhoods = vec![BondSet::new(n_bonds); n_atoms as usize];
    let mut round_neighborhoods = atom_neighborhoods.clone();
    let mut dead_atoms = vec![false; n_atoms as usize];

    let atom_order: Vec<u32> = if margs.only_nonzero_invariants {
        let mut ordering: Vec<(i32, u32)> = (0..n_atoms)
            .map(|i| (i32::from(current[i as usize] == 0), i))
            .collect();
        ordering.sort();
        ordering.into_iter().map(|(_, i)| i).collect()
    } else {
        (0..n_atoms).collect()
    };

    for i in 0..n_atoms {
        if include_atoms[i as usize] && (!margs.only_nonzero_invariants || current[i as usize] != 0)
        {
            result.push(MorganAtomEnv {
                code: current[i as usize],
                atom_id: i,
                layer: 0,
            });
        }
    }

    let mut nbr_invs: Vec<(i32, u32)> = Vec::with_capacity(8);
    for layer in 0..margs.radius {
        // (neighborhood, invariant, atom index)
        let mut all_this_round: Vec<(BondSet, u64, u32)> = Vec::new();
        for &atom_idx in &atom_order {
            let ai = atom_idx as usize;
            if dead_atoms[ai] {
                continue;
            }
            if mol.atom_neighbors(atom_idx).next().is_none() {
                dead_atoms[ai] = true;
                continue;
            }
            nbr_invs.clear();
            for bond in mol.atom_bonds(atom_idx) {
                round_neighborhoods[ai].set(bond.idx());
                let o = bond.other_atom_idx(atom_idx)?;
                let nb = atom_neighborhoods[o as usize].clone();
                round_neighborhoods[ai].or_assign(&nb);
                nbr_invs.push((
                    bond_invariants[bond.idx() as usize] as i32,
                    current[o as usize] as u32,
                ));
            }
            nbr_invs.sort();
            let mut invar: u32 = layer;
            hash_combine(&mut invar, current[ai] as u32);
            let atom = mol.atom(atom_idx);
            let mut looks_chiral = atom.chiral_tag() != ChiralType::Unspecified;
            for (k, &(bt, inv)) in nbr_invs.iter().enumerate() {
                hash_combine(&mut invar, hash_pair(bt as u32, inv));
                if args.include_chirality && looks_chiral && !chiral_atoms[ai] {
                    if bt != BondType::Single as i32 || (k > 0 && inv == nbr_invs[k - 1].1) {
                        looks_chiral = false;
                    }
                }
            }
            if args.include_chirality && looks_chiral {
                chiral_atoms[ai] = true;
                let cip = atom
                    .get_prop(cp::_CIPCode)
                    .map(|v| v.to_rd_string())
                    .unwrap_or_default();
                let v = match cip.as_str() {
                    "R" => 3,
                    "S" => 2,
                    _ => 1,
                };
                hash_combine(&mut invar, v);
            }
            next_layer[ai] = u64::from(invar);
            all_this_round.push((round_neighborhoods[ai].clone(), u64::from(invar), atom_idx));
        }

        all_this_round.sort_by(|a, b| a.0.cmp_bits(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        for (nbhd, inv, idx) in all_this_round {
            if margs.include_redundant_environments || !neighborhoods.contains(&nbhd) {
                if (!margs.only_nonzero_invariants || atom_invariants[idx as usize] != 0)
                    && include_atoms[idx as usize]
                {
                    result.push(MorganAtomEnv {
                        code: inv,
                        atom_id: idx,
                        layer: layer + 1,
                    });
                    neighborhoods.insert(nbhd);
                }
            } else {
                dead_atoms[idx as usize] = true;
            }
        }
        std::mem::swap(&mut current, &mut next_layer);
        next_layer.iter_mut().for_each(|v| *v = 0);
        atom_neighborhoods = round_neighborhoods.clone();
    }
    Ok(result)
}

/// `defaultFeatureSmarts`: Donor, Acceptor, Aromatic, Halogen, Basic, Acidic
/// (adapted from Gobbi and Poppinger, Biotech. Bioeng. 61 47-54 (1998)).
pub const DEFAULT_FEATURE_SMARTS: [&str; 6] = [
    "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
    "[a]",
    "[F,Cl,Br,I]",
    "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
    "[$([C,S](=[O,S,P])-[O;H1,-1])]",
];

fn default_feature_matchers() -> Result<&'static [ROMol]> {
    static MATCHERS: std::sync::OnceLock<Vec<ROMol>> = std::sync::OnceLock::new();
    if let Some(m) = MATCHERS.get() {
        return Ok(m);
    }
    let params = crate::graphmol::smiles_parse::smarts::SmartsParserParams::default();
    let mut v = Vec::new();
    for sma in DEFAULT_FEATURE_SMARTS {
        let m = crate::graphmol::smiles_parse::smarts::mol_from_smarts(sma, &params)?
            .ok_or_else(|| RDKitError::invariant("bad smarts"))?;
        v.push(m);
    }
    Ok(MATCHERS.get_or_init(|| v))
}

/// `MorganFingerprints::getFeatureInvariants`
pub fn get_feature_invariants(mol: &ROMol, patterns: Option<&[ROMol]>) -> Result<Vec<u32>> {
    let queries = match patterns {
        Some(p) => p,
        None => default_feature_matchers()?,
    };
    let mut invars = vec![0u32; mol.num_atoms() as usize];
    let ps = crate::graphmol::substruct::SubstructMatchParameters::default();
    for (i, q) in queries.iter().enumerate() {
        let mask = 1u32.wrapping_shl(i as u32);
        for m in crate::graphmol::substruct::substruct_match(mol, q, &ps)? {
            for (_, a) in m {
                invars[a as usize] |= mask;
            }
        }
    }
    Ok(invars)
}
