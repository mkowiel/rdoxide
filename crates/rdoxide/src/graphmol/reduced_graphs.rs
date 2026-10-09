//! Extended reduced graphs and ErG fingerprints (upstream
//! `GraphMol/ReducedGraphs/ReducedGraphs.{h,cpp}`).
//!
//! Stiefl et al., "ErG: 2D Pharmacophore Descriptions for Scaffold Hopping",
//! J. Chem. Inf. Model. 46, 208 (2006). A molecule is reduced to its
//! pharmacophoric features (donor, acceptor, positive, negative,
//! hydrophobic) plus a dummy atom at every ring centroid; the fingerprint
//! counts feature pairs by type and topological distance.

use std::sync::OnceLock;

use crate::error::{RDKitError, Result};
use crate::graphmol::atom::{Atom, HybridizationType};
use crate::graphmol::bond::BondType;
use crate::graphmol::matrices::get_distance_mat;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::{RDProps, RDValue};

/// The atom property holding the feature types (`Vec<i32>`) of an atom of
/// the reduced graph.
pub const ERG_ATOM_TYPES_PROP: &str = "_ErGAtomTypes";

/// Number of pharmacophore features, `nFeatures`.
pub const N_FEATURES: usize = 5;

/// `smartsPatterns`: donor, acceptor, positive, negative and hydrophobic
/// (which must stay last: its index doubles as the "aliphatic ring" type).
const FEATURE_SMARTS: [&str; N_FEATURES] = [
    "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O;H0;v2]),$([O,S;v1;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o;+0;!$([o]:n);!$([o]:c:n)])]",
    "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
    "[$([C,S](=[O,S,P])-[O;H1,-1])]",
    "[$([C;D3,D4](-[CH3])-[CH3]),$([S;D2](-C)-C)]",
];

/// For each feature type, whether each atom of the molecule has it.
pub type AtomTypes = Vec<Vec<bool>>;

fn feature_patterns() -> Result<&'static [ROMol]> {
    static PATTERNS: OnceLock<Vec<ROMol>> = OnceLock::new();
    if let Some(p) = PATTERNS.get() {
        return Ok(p);
    }
    let params = SmartsParserParams::default();
    let parsed = FEATURE_SMARTS
        .iter()
        .map(|sma| {
            mol_from_smarts(sma, &params)?.ok_or_else(|| RDKitError::invariant("bad smarts"))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(PATTERNS.get_or_init(|| parsed))
}

/// `getErGAtomTypes`: which atoms match which pharmacophore pattern.
pub fn get_erg_atom_types(mol: &ROMol) -> Result<AtomTypes> {
    let params = SubstructMatchParameters::default();
    feature_patterns()?
        .iter()
        .map(|pattern| {
            let mut has_type = vec![false; mol.num_atoms() as usize];
            for m in substruct_match(mol, pattern, &params)? {
                has_type[m[0].1 as usize] = true;
            }
            Ok(has_type)
        })
        .collect()
}

/// The feature types of an atom of a reduced graph.
fn erg_types_of(atom: &Atom) -> Result<&[i32]> {
    match atom.get_prop(ERG_ATOM_TYPES_PROP) {
        Some(RDValue::VecInt(v)) => Ok(v),
        _ => Err(RDKitError::Key(ERG_ATOM_TYPES_PROP.into())),
    }
}

/// `generateMolExtendedReducedGraph`: the molecule with a dummy atom bonded
/// to every ring of fewer than eight atoms, and without ring atoms of degree
/// two that carry no feature. Every atom has its feature types in
/// [`ERG_ATOM_TYPES_PROP`]; a dummy atom has the "aromatic" type (one past
/// the last feature type) or the "aliphatic" one (the last feature type).
pub fn generate_mol_extended_reduced_graph(
    mol: &ROMol,
    atom_types: Option<&[Vec<bool>]>,
) -> Result<RWMol> {
    let computed;
    let atom_types = match atom_types {
        Some(t) => t,
        None => {
            computed = get_erg_atom_types(mol)?;
            &computed
        }
    };
    if !mol.ring_info().is_initialized() {
        return Err(RDKitError::precondition("RingInfo not initialized"));
    }
    let aliphatic_flag = atom_types.len() as i32 - 1;
    let aromatic_flag = atom_types.len() as i32;

    let mut res = mol.clone();
    for i in 0..mol.num_atoms() {
        let types: Vec<i32> = (0..atom_types.len())
            .filter(|&t| atom_types[t][i as usize])
            .map(|t| t as i32)
            .collect();
        res.atom_mut(i)
            .set_prop(ERG_ATOM_TYPES_PROP, types, false)?;
    }

    // a dummy atom at the centroid of every small ring
    for ring in mol.ring_info().atom_rings() {
        if ring.len() >= 8 {
            continue;
        }
        let centroid = res.add_atom(Atom::new(0), false);
        let (mut n_aromatic, mut n_sp2) = (0, 0);
        for &idx in ring {
            res.add_bond(idx as u32, centroid, BondType::Single)?;
            let atom = mol.atom(idx as u32);
            if atom.is_aromatic() {
                n_aromatic += 1;
            } else if atom.hybridization() == HybridizationType::SP2 {
                n_sp2 += 1;
            }
        }
        let flag = if n_aromatic >= 2 || n_sp2 >= (ring.len() / 2) as i32 {
            aromatic_flag
        } else {
            aliphatic_flag
        };
        res.atom_mut(centroid)
            .set_prop(ERG_ATOM_TYPES_PROP, vec![flag], false)?;
    }

    // ring atoms of degree two without features carry no information
    res.begin_batch_edit()?;
    for i in (0..mol.num_atoms()).rev() {
        if mol.ring_info().num_atom_rings(i)? != 0
            && mol.atom_degree(i) == 2
            && erg_types_of(res.atom(i))?.is_empty()
        {
            res.remove_atom(i)?;
        }
    }
    res.commit_batch_edit()?;
    Ok(res)
}

/// `generateErGFingerprintForReducedGraph`: counts of feature pairs by type
/// pair and topological distance in `min_path..=max_path`, with
/// `fuzz_increment` added to the neighboring distance bins.
pub fn generate_erg_fingerprint_for_reduced_graph(
    mol: &ROMol,
    atom_types: Option<&[Vec<bool>]>,
    fuzz_increment: f64,
    min_path: u32,
    max_path: u32,
) -> Result<Vec<f64>> {
    if max_path <= min_path {
        return Err(RDKitError::precondition("maxPath<=minPath"));
    }
    // (the pair table has one extra type for ring centroids)
    let n_types = atom_types.map_or(N_FEATURES, <[_]>::len) + 1;
    let n_bins = (max_path - min_path + 1) as usize;
    let size = n_types * (n_types + 1) / 2 * n_bins;
    let mut res = vec![0.0; size];

    let dist_mat = get_distance_mat(mol);
    let n_atoms = mol.num_atoms() as usize;
    let types: Vec<&[i32]> = mol
        .atoms()
        .iter()
        .map(erg_types_of)
        .collect::<Result<_>>()?;

    for i in 0..n_atoms {
        if types[i].is_empty() {
            continue;
        }
        for j in i + 1..n_atoms {
            if types[j].is_empty() {
                continue;
            }
            let dist = dist_mat[i * n_atoms + j] as i64;
            if dist < i64::from(min_path) || dist > i64::from(max_path) {
                continue;
            }
            for &ti in types[i] {
                for &tj in types[j] {
                    let (lo, hi) = (ti.min(tj) as usize, ti.max(tj) as usize);
                    let block = if lo > 0 {
                        lo * n_types - lo * (lo + 1) / 2
                    } else {
                        0
                    } + hi;
                    let bin = block * n_bins + (dist - i64::from(min_path)) as usize;
                    res[bin] += 1.0;
                    if dist > i64::from(min_path) && fuzz_increment > 0.0 {
                        res[bin - 1] += fuzz_increment;
                    }
                    if dist < i64::from(max_path) && fuzz_increment > 0.0 && bin + 1 < size {
                        res[bin + 1] += fuzz_increment;
                    }
                }
            }
        }
    }
    Ok(res)
}

/// `getErGFingerprint`: the ErG fingerprint of a molecule.
pub fn get_erg_fingerprint(
    mol: &ROMol,
    atom_types: Option<&[Vec<bool>]>,
    fuzz_increment: f64,
    min_path: u32,
    max_path: u32,
) -> Result<Vec<f64>> {
    let reduced = generate_mol_extended_reduced_graph(mol, atom_types)?;
    generate_erg_fingerprint_for_reduced_graph(
        &reduced,
        atom_types,
        fuzz_increment,
        min_path,
        max_path,
    )
}
