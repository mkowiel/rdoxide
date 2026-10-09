//! `HetAtomTautomerv2` / `HetAtomProtomerv2` (upstream `TautomerHashv2`):
//! finds the mobile-hydrogen conjugated systems of a molecule, makes their
//! bonds aromatic and counts the hydrogens and charges they hold.

use std::collections::VecDeque;
use std::sync::OnceLock;

use super::{
    Cx, EXCLUDE_FROM_TAUTOMERISM_PROP, add_cx_extensions, clean_up, ensure_rings, in_ring, smiles,
};
use crate::error::Result;
use crate::graphmol::ChiralType;
use crate::graphmol::atom::HybridizationType;
use crate::graphmol::bond::{Bond, BondStereo, BondType};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::cx_write::cx_fields;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;

/// Bond is involved in a carboxyl, amide, ... (`bondFlagCarboxyl`).
const BOND_FLAG_CARBOXYL: u64 = 1;
/// Bond carries the exclusion property (`bondFlagProperty`).
const BOND_FLAG_PROPERTY: u64 = 2;

/// Functional groups whose bonds never take part in a tautomeric system; the
/// first two atoms of each pattern are the bond that gets flagged.
const EXCLUDED_GROUPS: [&str; 12] = [
    // one side of the "amide", with an ugly exclusion for amidine
    "[C;!$(C-C(=[NH])-[NH2])]-[C;!$(C(-C)(=[NH])-[NH2])](=[O,N,S])-[O,N,S]",
    "[A;!$(C=[O,N,S])]-[O,N,S]-C=[O,N,S]", // the other side
    "[OH0,SH0]-C=[O,N,S]",                 // "esters" and "carboxyls"
    "C-[N;$([N+]=C(N)(N)),$(N-C(N)=N)]",   // guanidine like
    "[C]-[c](:[o,n,s]):[o,n,s]",           // a limited version for aromatic systems
    "*[SD4](=*)=*",                        // sulfates, sulfonyls, etc.
    "*=[SD4;$(S(=*)=*)](*)*",              // the other side
    "*-[H0]=,#[C,N]=,#*",                  // isocyanates, azides, et al
    "*[#7+]-[O-]",                         // nitro groups and aromatic n-oxides
    "[O,S;H]-P=O",                         // phosphoric acid, etc.
    "[#6]-P=O",                            // phosphoric acid, etc.
    "[#6]-N=[SD4]=*",                      // don't know what this one is called
];

fn excluded_group_queries() -> Result<&'static [ROMol]> {
    static QUERIES: OnceLock<Vec<ROMol>> = OnceLock::new();
    if let Some(q) = QUERIES.get() {
        return Ok(q);
    }
    let params = SmartsParserParams::default();
    let parsed = EXCLUDED_GROUPS
        .iter()
        .map(|sma| {
            mol_from_smarts(sma, &params)?
                .ok_or_else(|| crate::RDKitError::invariant("bad MolHash SMARTS"))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(QUERIES.get_or_init(|| parsed))
}

/// `details::getBondFlags`
fn bond_flags(mol: &ROMol) -> Result<Vec<u64>> {
    let mut flags: Vec<u64> = mol
        .bonds()
        .iter()
        .map(|b| {
            if b.has_prop(EXCLUDE_FROM_TAUTOMERISM_PROP) {
                BOND_FLAG_PROPERTY
            } else {
                0
            }
        })
        .collect();
    for query in excluded_group_queries()? {
        for found in substruct_match(mol, query, &SubstructMatchParameters::default())? {
            let bond = mol
                .get_bond_between_atoms(found[0].1 as u32, found[1].1 as u32)?
                .expect("the first two query atoms are bonded");
            flags[bond.idx() as usize] |= BOND_FLAG_CARBOXYL;
        }
    }
    Ok(flags)
}

/// The molecule plus the per-atom/bond exclusion flags and the start bonds,
/// with the predicates of the upstream helper functions.
struct Systems<'a> {
    mol: &'a ROMol,
    atom_flags: Vec<bool>,
    bond_flags: Vec<u64>,
    start_bonds: Vec<bool>,
    conjugated_neighbors: Vec<u32>,
}

fn is_unsaturated_bond(b: &Bond) -> bool {
    b.is_aromatic()
        || matches!(
            b.bond_type(),
            BondType::Aromatic | BondType::Double | BondType::Triple
        )
}

fn is_hetero(mol: &ROMol, idx: u32) -> bool {
    let n = mol.atom(idx).atomic_num();
    n != 6 && n > 1
}

impl<'a> Systems<'a> {
    fn new(mol: &'a ROMol) -> Result<Self> {
        let atom_flags = mol
            .atoms()
            .iter()
            .map(|a| a.has_prop(EXCLUDE_FROM_TAUTOMERISM_PROP))
            .collect();
        let mut s = Systems {
            mol,
            atom_flags,
            bond_flags: bond_flags(mol)?,
            start_bonds: vec![false; mol.num_bonds() as usize],
            conjugated_neighbors: Vec::new(),
        };
        s.start_bonds = (0..mol.num_bonds())
            .map(|b| s.is_possible_starting_bond(b))
            .collect();
        s.conjugated_neighbors = (0..mol.num_atoms())
            .map(|a| s.num_conjugated_neighbors(a))
            .collect();
        Ok(s)
    }

    fn num_hs(&self, atom: u32) -> u32 {
        self.mol.total_num_hs(atom, false).unwrap_or(0)
    }

    /// `queryAtomUnsaturated`
    fn is_unsaturated(&self, atom: u32) -> bool {
        match (
            self.mol.total_degree(atom),
            self.mol.atom(atom).total_valence(),
        ) {
            (Ok(degree), Ok(valence)) => degree < valence,
            _ => false,
        }
    }

    /// Candidate atoms have implicit Hs or are unsaturated. Being aromatic is
    /// not enough: in `Cn1cccc1` the N is aromatic but not tautomeric.
    fn is_candidate(&self, atom: u32) -> bool {
        !self.atom_flags[atom as usize] && (self.num_hs(atom) > 0 || self.is_unsaturated(atom))
    }

    /// A bond between a candidate heteroatom and an unsaturated atom.
    fn is_possible_starting_bond(&self, bond: u32) -> bool {
        let b = self.mol.bond(bond);
        let (beg, end) = (b.begin_atom_idx(), b.end_atom_idx());
        if self.bond_flags[bond as usize] != 0
            || self.atom_flags[beg as usize]
            || self.atom_flags[end as usize]
        {
            return false;
        }
        let hetero_beg = is_hetero(self.mol, beg) && self.is_candidate(beg);
        let hetero_end = is_hetero(self.mol, end) && self.is_candidate(end);
        // no check for aromaticity here, see `is_candidate`
        (hetero_beg && self.is_unsaturated(end)) || (hetero_end && self.is_unsaturated(beg))
    }

    fn has_start_bond(&self, atom: u32) -> bool {
        self.mol
            .atom_bond_indices(atom)
            .any(|b| self.start_bonds[b as usize])
    }

    /// Neighbors of `atom` with a conjugated or starting bond.
    fn num_conjugated_neighbors(&self, atom: u32) -> u32 {
        self.mol
            .atom_neighbors(atom)
            .filter(|&n| {
                self.mol
                    .atom_bonds(n)
                    .any(|b| b.is_conjugated() || self.start_bonds[b.idx() as usize])
            })
            .count() as u32
    }

    /// Prevents "overreach" with enamines: the first bond of `CNC=C` is kept
    /// out of the system, giving `[CH3]-[N]:[C]:[C]`.
    fn is_overreach(&self, atom: u32, oatom: u32, bond: u32, nbr_bond: u32) -> bool {
        (self.start_bonds[bond as usize] || self.has_start_bond(atom))
            && !self.start_bonds[nbr_bond as usize]
            && is_hetero(self.mol, atom)
            && !is_unsaturated_bond(self.mol.bond(nbr_bond))
            && self.conjugated_neighbors[oatom as usize] < 2
    }

    /// Does not extend to a stereocenter over a single non-conjugated bond,
    /// unless the source atom has a non-aromatic unsaturated bond.
    fn should_skip_chiral_neighbor(&self, source: u32, target: u32, bond: u32) -> bool {
        if self.mol.atom(target).hybridization() != HybridizationType::SP3
            || self.num_hs(target) != 1
        {
            return false;
        }
        let connecting = self.mol.bond(bond);
        if is_unsaturated_bond(connecting) || connecting.is_conjugated() {
            return false;
        }
        !self
            .mol
            .atom_bonds(source)
            .any(|b| !b.is_aromatic() && is_unsaturated_bond(b))
    }

    /// Stops an aromatic carbon without H from extending into an exocyclic
    /// C=C, which would destroy its E/Z stereo.
    fn leaves_aromatic_for_vinyl(&self, atom: u32, other: u32) -> bool {
        let (a, o) = (self.mol.atom(atom), self.mol.atom(other));
        let (aromatic, plain) = match (a.is_aromatic(), o.is_aromatic()) {
            (true, false) => (atom, other),
            (false, true) => (other, atom),
            _ => return false,
        };
        if self.num_hs(aromatic) != 0
            || self.mol.atom(aromatic).atomic_num() != 6
            || self.mol.atom(plain).atomic_num() != 6
        {
            return false;
        }
        self.mol.atom_bonds(plain).any(|b| {
            b.bond_type() == BondType::Double
                && !b.is_aromatic()
                && b.other_atom_idx(plain)
                    .is_ok_and(|e| self.mol.atom(e).atomic_num() == 6)
        })
    }

    /// `skipNeighborBond`
    fn skip_neighbor_bond(&self, atom: u32, other: u32, nbr_bond: u32) -> bool {
        let b = self.mol.bond(nbr_bond);
        if self.bond_flags[nbr_bond as usize] != 0
            || self.atom_flags[other as usize]
            || self.atom_flags[atom as usize]
        {
            return true;
        }
        if !b.is_aromatic() && self.leaves_aromatic_for_vinyl(atom, other) {
            return true;
        }
        (!self.is_candidate(other) && !self.has_start_bond(other))
            || (!is_unsaturated_bond(b) && !b.is_conjugated() && !self.has_start_bond(atom))
    }

    /// Both directions of [`skip_neighbor_bond`](Self::skip_neighbor_bond)
    /// and the other reasons not to follow `nbr_bond` from `atom`.
    fn skip_bond(&self, atom: u32, oatom: u32, from_bond: u32, nbr_bond: u32) -> bool {
        self.is_overreach(atom, oatom, from_bond, nbr_bond)
            || self.should_skip_chiral_neighbor(atom, oatom, nbr_bond)
            || (self.skip_neighbor_bond(atom, oatom, nbr_bond)
                && self.skip_neighbor_bond(oatom, atom, nbr_bond))
    }

    /// Grows the conjugated system of the start bond `start` and returns its
    /// bonds, or `None` when it has no mobile hydrogen.
    fn grow_system(
        &self,
        start: u32,
        considered: &mut [bool],
        to_modify: &mut [bool],
    ) -> Option<()> {
        let mol = self.mol;
        let nb = mol.num_bonds() as usize;
        let (mut system, mut system_atoms, mut visited_atoms) = (
            vec![false; nb],
            vec![false; mol.num_atoms() as usize],
            vec![false; mol.num_atoms() as usize],
        );
        let (mut donor_cs, mut hetero_hs) = (0u32, 0u32);
        let mut queue: VecDeque<u32> = VecDeque::new();
        let sb = mol.bond(start);
        for atm in [sb.begin_atom_idx(), sb.end_atom_idx()] {
            if mol.atom(atm).atomic_num() == 6 {
                donor_cs += u32::from(self.num_hs(atm) > 0);
            } else if is_hetero(mol, atm) {
                hetero_hs += self.num_hs(atm);
            }
            for nbr_bond in mol.atom_bond_indices(atm).collect::<Vec<_>>() {
                if nbr_bond == start || considered[nbr_bond as usize] {
                    continue;
                }
                let oatom = mol.bond(nbr_bond).other_atom_idx(atm).ok()?;
                if self.skip_bond(atm, oatom, start, nbr_bond) {
                    continue;
                }
                if !queue.contains(&nbr_bond) {
                    queue.push_back(nbr_bond);
                }
                considered[start as usize] = true;
                system[start as usize] = true;
                system_atoms[sb.begin_atom_idx() as usize] = true;
                system_atoms[sb.end_atom_idx() as usize] = true;
            }
        }
        while let Some(bnd) = queue.pop_front() {
            if considered[bnd as usize] {
                continue;
            }
            considered[bnd as usize] = true;
            system[bnd as usize] = true;
            let b = mol.bond(bnd);
            system_atoms[b.begin_atom_idx() as usize] = true;
            system_atoms[b.end_atom_idx() as usize] = true;
            for atm in [b.begin_atom_idx(), b.end_atom_idx()] {
                if visited_atoms[atm as usize] {
                    continue;
                }
                if mol.atom(atm).atomic_num() == 6 {
                    if self.num_hs(atm) > 0 {
                        donor_cs += 1;
                        visited_atoms[atm as usize] = true;
                    }
                } else if mol.atom(atm).atomic_num() > 1 {
                    hetero_hs += self.num_hs(atm);
                    visited_atoms[atm as usize] = true;
                }
                for nbr_bnd in mol.atom_bond_indices(atm).collect::<Vec<_>>() {
                    if considered[nbr_bnd as usize] || queue.contains(&nbr_bnd) {
                        continue;
                    }
                    let oatom = mol.bond(nbr_bnd).other_atom_idx(atm).ok()?;
                    if !self.skip_bond(atm, oatom, bnd, nbr_bnd) {
                        queue.push_back(nbr_bnd);
                    }
                }
            }
        }
        // at least two bonds and one active H
        if system.iter().filter(|&&s| s).count() > 1 && (hetero_hs > 0 || donor_cs > 0) {
            // the traversal may miss ring closure bonds: every bond between
            // two atoms of the system belongs to it
            for atom in (0..mol.num_atoms()).filter(|&a| system_atoms[a as usize]) {
                for b in mol.atom_bonds(atom) {
                    if b.other_atom_idx(atom)
                        .is_ok_and(|o| system_atoms[o as usize])
                    {
                        to_modify[b.idx() as usize] = true;
                    }
                }
            }
            Some(())
        } else {
            // forget the bonds marked as considered for this system
            for (i, &in_system) in system.iter().enumerate() {
                if in_system {
                    considered[i] = false;
                }
            }
            None
        }
    }

    /// Two tautomers must hash alike when only one of them has the bond next
    /// to a mobile H conjugated: extend the system over bonds whose far
    /// neighbors already hold a start or system bond.
    fn extend_to_neighbors(&self, to_modify: &mut [bool]) {
        let mol = self.mol;
        for start in (0..mol.num_bonds()).filter(|&b| self.start_bonds[b as usize]) {
            let sb = mol.bond(start);
            for atm in [sb.begin_atom_idx(), sb.end_atom_idx()] {
                for nbr_bond in mol.atom_bond_indices(atm).collect::<Vec<_>>() {
                    if nbr_bond == start || to_modify[nbr_bond as usize] {
                        continue;
                    }
                    let Ok(oatom) = mol.bond(nbr_bond).other_atom_idx(atm) else {
                        continue;
                    };
                    if self.num_hs(oatom) == 0 {
                        continue;
                    }
                    // Do not pull a potential stereocenter into the system:
                    // chain atoms are judged structurally (SP3, one H), ring
                    // atoms only when chirality is annotated, since ring SP3+1H
                    // atoms often do take part in ring tautomerism.
                    if mol.atom(oatom).hybridization() == HybridizationType::SP3
                        && self.num_hs(oatom) == 1
                        && (!in_ring(mol, oatom)
                            || mol.atom(oatom).chiral_tag() != ChiralType::Unspecified)
                    {
                        continue;
                    }
                    let near_system = mol.atom_neighbors(oatom).filter(|&n| n != atm).any(|n| {
                        mol.atom_bond_indices(n)
                            .any(|b| to_modify[b as usize] || self.start_bonds[b as usize])
                    });
                    if near_system {
                        to_modify[nbr_bond as usize] = true;
                    }
                }
            }
        }
    }
}

/// `TautomerHashv2` (`proto`: the protomer hash)
pub(super) fn tautomer_hash_v2(mol: &mut RWMol, proto: bool, cx: Cx) -> Result<String> {
    ensure_rings(mol)?;
    let nb = mol.num_bonds() as usize;
    let (mut to_modify, atoms_to_modify) = {
        let systems = Systems::new(mol)?;
        let (mut considered, mut to_modify) = (vec![false; nb], vec![false; nb]);
        for b in 0..nb {
            if to_modify[b] || considered[b] || !systems.start_bonds[b] {
                continue;
            }
            systems.grow_system(b as u32, &mut considered, &mut to_modify);
        }
        systems.extend_to_neighbors(&mut to_modify);
        let mut atoms = vec![false; mol.num_atoms() as usize];
        for b in (0..nb).filter(|&b| to_modify[b]) {
            atoms[mol.bond(b as u32).begin_atom_idx() as usize] = true;
            atoms[mol.bond(b as u32).end_atom_idx() as usize] = true;
        }
        (to_modify, atoms)
    };

    for b in (0..nb).filter(|&b| to_modify[b]) {
        // E/Z stereo of exocyclic double bonds is kept (as DOUBLE bonds) so
        // that geometric isomers, e.g. of hydrazones, stay distinct
        let bond = mol.bond(b as u32);
        let exocyclic_stereo = bond.stereo() != BondStereo::StereoNone
            && in_ring(mol, bond.begin_atom_idx()) != in_ring(mol, bond.end_atom_idx());
        let bond = mol.bond_mut(b as u32);
        if exocyclic_stereo {
            bond.set_is_aromatic(false);
        } else {
            bond.set_bond_type(BondType::Aromatic);
            bond.set_is_aromatic(true);
            bond.set_stereo(BondStereo::StereoNone)?;
        }
    }
    to_modify.clear();

    let (mut hcount, mut charge) = (0u32, 0i32);
    for atom in (0..mol.num_atoms()).filter(|&a| atoms_to_modify[a as usize]) {
        charge += mol.atom(atom).formal_charge();
        hcount += mol.total_num_hs(atom, false)?;
        let a = mol.atom_mut(atom);
        a.set_is_aromatic(false);
        a.set_formal_charge(0);
        a.set_no_implicit(true);
        a.set_num_explicit_hs(0);
    }
    // upstream tests `bitset::empty()`, i.e. size zero, here
    if mol.num_bonds() > 0 || mol.num_atoms() > 0 {
        clean_up(mol)?;
    }
    let mut result = smiles(mol, true, cx)?;
    result += &if proto {
        format!("_{}", hcount as i32 - charge)
    } else {
        format!("_{hcount}_{charge}")
    };
    add_cx_extensions(mol, &mut result, cx, cx_fields::RADICALS)?;
    Ok(result)
}
