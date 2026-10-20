//! Atom feature vectors (upstream `GraphMol/Descriptors/AtomFeat.cpp`): the
//! one-hot encoded symbol, degree, hybridization, valence, charge, ring sizes
//! and hydrogen count of an atom, optionally with its chirality.
//!
//! Copyright (C) 2020 Guillaume GODIN and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::graphmol::atom::{ChiralType, HybridizationType, ValenceType};
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::romol::ROMol;

pub const ATOM_FEAT_VERSION: &str = "1.0.0";

const SYMBOLS: [&str; 13] = [
    "B", "C", "N", "O", "S", "F", "Si", "P", "Cl", "Br", "I", "H", "*",
];
const HYBRIDIZATIONS: [HybridizationType; 5] = [
    HybridizationType::SP,
    HybridizationType::SP2,
    HybridizationType::SP3,
    HybridizationType::SP3D,
    HybridizationType::SP3D2,
];
const CHIRAL_TAGS: [ChiralType; 3] = [
    ChiralType::TetrahedralCW,
    ChiralType::TetrahedralCCW,
    ChiralType::Other,
];

/// `AtomFeatVect`: 49 values, 52 with `add_chiral`.
pub fn atom_feat_vect(mol: &mut ROMol, atom_id: u32, add_chiral: bool) -> Result<Vec<f64>> {
    mol.get_atom(atom_id)?;
    if !mol.ring_info().is_sssr_or_better() {
        find_sssr(mol, false, false)?;
    }
    let atom = mol.atom(atom_id);
    let mut feats = Vec::with_capacity(52);
    let mut one_hot = |hit: bool| feats.push(f64::from(u8::from(hit)));

    // symbol, with a trailing "unknown" slot
    let symbol = atom.symbol();
    for s in SYMBOLS {
        one_hot(s == symbol);
    }
    one_hot(!SYMBOLS.contains(&symbol.as_str()));
    let degree = mol.atom_degree(atom_id);
    for i in 0..7 {
        one_hot(degree == i);
    }
    for h in HYBRIDIZATIONS {
        one_hot(atom.hybridization() == h);
    }
    let implicit = atom.valence(ValenceType::Implicit)?;
    for i in 0..7 {
        one_hot(implicit == i);
    }
    for charge in -1..2 {
        one_hot(atom.formal_charge() == charge);
    }
    for size in 3..9 {
        one_hot(mol.ring_info().is_atom_in_ring_of_size(atom_id, size)?);
    }
    one_hot(atom.is_aromatic());
    let total_hs = mol.total_num_hs(atom_id, false)?;
    for i in 0..5 {
        one_hot(total_hs == i);
    }
    feats.push(1.0 / f64::from(mol.num_atoms()));
    if add_chiral {
        for tag in CHIRAL_TAGS {
            feats.push(f64::from(u8::from(atom.chiral_tag() == tag)));
        }
    }
    Ok(feats)
}
