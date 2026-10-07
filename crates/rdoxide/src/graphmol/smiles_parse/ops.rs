//! Port of RDKit `Code/GraphMol/SmilesParse/SmilesParseOps.cpp` (the parts
//! used by the SMILES parser) and `Canon::chiralAtomNeedsTagInversion`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;

use crate::graphmol::atom::{Atom, ChiralType, ValenceType};
use crate::graphmol::bond::{Bond, BondDir, BondType};
use crate::graphmol::chirality::{get_allow_nontetrahedral_chirality, nontetrahedral};
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};

/// A `SmilesParseException` message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesParseError(pub String);

/// `detail::printSyntaxErrorMessage`
pub fn print_syntax_error_message(
    input: &str,
    err_message: &str,
    bad_token_position: usize,
    input_type: &str,
) {
    const ERROR_SIZE: usize = 41;
    const PREFIX_SIZE: usize = ERROR_SIZE / 2;
    let truncate = |pos: usize| -> &str {
        let n = input.len();
        if pos >= PREFIX_SIZE && pos + PREFIX_SIZE < n {
            &input[pos - PREFIX_SIZE..pos - PREFIX_SIZE + ERROR_SIZE]
        } else if pos >= PREFIX_SIZE {
            &input[(pos - PREFIX_SIZE).min(n)..]
        } else {
            &input[..n.min(ERROR_SIZE)]
        }
    };
    let num_dashes = if bad_token_position >= PREFIX_SIZE {
        PREFIX_SIZE
    } else {
        bad_token_position.saturating_sub(1)
    };
    crate::rd_log_error!("{input_type} Parse Error: {err_message} while parsing: {input}");
    crate::rd_log_error!(
        "{input_type} Parse Error: check for mistakes around position {bad_token_position}:"
    );
    crate::rd_log_error!("{}", truncate(bad_token_position.saturating_sub(1)));
    crate::rd_log_error!("{}^", "~".repeat(num_dashes));
}

/// `GetUnspecifiedBondType`
pub fn get_unspecified_bond_type(a1: &Atom, a2: &Atom) -> BondType {
    if a1.is_aromatic() && a2.is_aromatic() {
        BondType::Aromatic
    } else {
        BondType::Single
    }
}

/// `CheckRingClosureBranchStatus`
pub fn check_ring_closure_branch_status(mol: &mut RWMol, atom: u32) {
    let deg = mol.atom_degree(atom);
    let a = mol.atom(atom);
    if atom != mol.num_atoms() - 1
        && (deg == 1 || (deg == 2 && atom != 0) || (deg == 3 && atom == 0))
        && matches!(
            a.chiral_tag(),
            ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
        )
    {
        mol.atom_mut(atom).invert_chirality();
    }
}

/// Dihedral-direction swap used when closing rings (`swapBondDirIfNeeded`).
fn swap_bond_dir_if_needed(bond1: &mut Bond, bond2: &Bond) {
    if bond1.bond_dir() == BondDir::None && bond2.bond_dir() != BondDir::None {
        bond1.set_bond_dir(bond2.bond_dir());
        if bond1.begin != bond2.begin {
            match bond1.bond_dir() {
                BondDir::EndDownRight => bond1.set_bond_dir(BondDir::EndUpRight),
                BondDir::EndUpRight => bond1.set_bond_dir(BondDir::EndDownRight),
                _ => {}
            }
        }
    }
}

fn could_be_ring_closure(v: i32) -> bool {
    (0..100_000).contains(&v)
}

fn ring_closures(mol: &RWMol, atom: u32) -> Vec<i32> {
    match mol.atom(atom).props.get(cp::_RingClosures) {
        Some(RDValue::VecInt(v)) => v.clone(),
        _ => Vec::new(),
    }
}

/// `CloseMolRings`: turn the pairs of ring-closure partial bonds into real
/// bonds. Ring numbers are processed in increasing order (std::map order).
pub fn close_mol_rings(
    mol: &mut RWMol,
    ring_atoms: &mut BTreeMap<i32, Vec<u32>>,
    ring_bonds: &mut BTreeMap<i32, Vec<Bond>>,
    tolerate_partials: bool,
) -> Result<(), SmilesParseError> {
    let marks: Vec<i32> = ring_atoms.keys().copied().collect();
    for mark in marks {
        if !could_be_ring_closure(mark) {
            continue;
        }
        let atoms = ring_atoms.get(&mark).cloned().unwrap_or_default();
        let mut removed: Vec<u32> = Vec::new();
        let mut i = 0;
        while i < atoms.len() {
            let atom1 = atoms[i];
            i += 1;
            if !tolerate_partials && i == atoms.len() {
                return Err(SmilesParseError("unclosed ring".into()));
            } else if i < atoms.len() && atoms[i] == atom1 {
                return Err(SmilesParseError(format!(
                    "duplicated ring closure {mark} bonds atom {atom1} to itself"
                )));
            } else if i < atoms.len() && mol.bond_idx_between(atom1, atoms[i]).is_some() {
                return Err(SmilesParseError(format!(
                    "ring closure {mark} duplicates bond between atom {atom1} and atom {}",
                    atoms[i]
                )));
            } else if i < atoms.len() {
                let atom2 = atoms[i];
                i += 1;
                let bonds = ring_bonds.get_mut(&mark).expect("Missing bond bookmark");
                let mut bond1 = bonds.remove(0);
                let bond2 = bonds.remove(0);
                if let Some(v) = bond2.props.get("_cxsmilesBondIdx").cloned() {
                    bond1.props.set_val("_cxsmilesBondIdx", v);
                }
                let mut matched;
                if !bond1.has_prop(cp::_unspecifiedOrder) {
                    matched = bond1;
                    match matched.bond_type() {
                        BondType::DativeL => {
                            matched.begin = atom2;
                            matched.end = atom1;
                            matched.set_bond_type(BondType::Dative);
                        }
                        BondType::DativeR => {
                            matched.end = atom2;
                            matched.set_bond_type(BondType::Dative);
                        }
                        _ => matched.end = atom2,
                    }
                    swap_bond_dir_if_needed(&mut matched, &bond2);
                } else {
                    matched = bond2;
                    match matched.bond_type() {
                        BondType::DativeL => {
                            matched.begin = atom1;
                            matched.end = atom2;
                            matched.set_bond_type(BondType::Dative);
                        }
                        BondType::DativeR => {
                            matched.end = atom1;
                            matched.set_bond_type(BondType::Dative);
                        }
                        _ => matched.end = atom1,
                    }
                    swap_bond_dir_if_needed(&mut matched, &bond1);
                }
                if matched.bond_type() == BondType::Unspecified && !matched.has_query() {
                    let t = get_unspecified_bond_type(mol.atom(atom1), mol.atom(atom2));
                    matched.set_bond_type(t);
                }
                if matched.bond_type() == BondType::Aromatic {
                    matched.set_is_aromatic(true);
                }
                let n = mol
                    .add_bond_object(matched)
                    .map_err(|e| SmilesParseError(e.to_string()))?;
                let bond_idx = n as i32 - 1;
                for a in [atom1, atom2] {
                    let mut closures = ring_closures(mol, a);
                    if let Some(p) = closures.iter().position(|&c| c == -(mark + 1)) {
                        closures[p] = bond_idx;
                    }
                    mol.atom_mut(a).props.set_val(cp::_RingClosures, closures);
                }
                removed.push(atom1);
                removed.push(atom2);
            }
        }
        if let Some(list) = ring_atoms.get_mut(&mark) {
            for a in removed {
                if let Some(p) = list.iter().position(|&x| x == a) {
                    list.remove(p);
                }
            }
            if list.is_empty() {
                ring_atoms.remove(&mark);
            }
        }
    }
    Ok(())
}

/// `checkChiralPermutation` limits
fn permutation_limit(tag: ChiralType) -> Option<i32> {
    match tag {
        ChiralType::Tetrahedral | ChiralType::Allene => Some(2),
        ChiralType::SquarePlanar => Some(3),
        ChiralType::Octahedral => Some(30),
        ChiralType::TrigonalBipyramidal => Some(20),
        _ => None,
    }
}

/// `CheckChiralitySpecifications`
pub fn check_chirality_specifications(
    mol: &mut RWMol,
    strict: bool,
) -> Result<(), SmilesParseError> {
    for i in 0..mol.num_atoms() {
        let tag = mol.atom(i).chiral_tag();
        if tag <= ChiralType::Other {
            continue;
        }
        let Some(limit) = permutation_limit(tag) else {
            continue;
        };
        let Some(perm) = mol
            .atom(i)
            .props
            .get(cp::_chiralPermutation)
            .and_then(RDValue::as_int)
        else {
            continue;
        };
        if perm < 0 || perm > limit {
            let error = format!("Invalid chiral specification on atom {i}");
            crate::rd_log_warning!("{error}");
            if strict {
                return Err(SmilesParseError(error));
            }
        }
        if tag == ChiralType::Tetrahedral {
            let a = mol.atom_mut(i);
            if perm == 0 || perm == 1 {
                a.set_chiral_tag(ChiralType::TetrahedralCCW);
                a.clear_prop(cp::_chiralPermutation);
            } else if perm == 2 {
                a.set_chiral_tag(ChiralType::TetrahedralCW);
                a.clear_prop(cp::_chiralPermutation);
            }
        }
    }
    Ok(())
}

/// `SetUnspecifiedBondTypes`
pub fn set_unspecified_bond_types(mol: &mut RWMol) {
    for i in 0..mol.num_bonds() {
        if mol.bond(i).has_prop(cp::_unspecifiedOrder) {
            let (b, e) = (mol.bond(i).begin, mol.bond(i).end);
            let t = get_unspecified_bond_type(mol.atom(b), mol.atom(e));
            let bond = mol.bond_mut(i);
            bond.set_bond_type(t);
            bond.set_is_aromatic(t == BondType::Aromatic);
        }
    }
}

/// `GetBondOrdering`: the bond order used to interpret a chiral tag written
/// in a SMILES string. Returns the number of ring closures.
pub fn get_bond_ordering(mol: &RWMol, atom: u32) -> (Vec<i32>, usize) {
    let closures = ring_closures(mol, atom);
    let mut neighbors: Vec<(u32, i64)> = vec![(atom, -1)];
    for (nbr, bidx) in mol.adj[atom as usize].iter().copied() {
        if !closures.contains(&(bidx as i32)) {
            neighbors.push((nbr, i64::from(bidx)));
        }
    }
    neighbors.sort();
    let mut ordering = Vec::new();
    for (n, b) in neighbors {
        if n == atom && b == -1 {
            ordering.extend_from_slice(&closures);
        } else {
            ordering.push(b as i32);
        }
    }
    (ordering, closures.len())
}

/// `Canon::details::isUnsaturated`
pub fn is_unsaturated(mol: &RWMol, atom: u32) -> bool {
    mol.atom_bonds(atom)
        .any(|b| b.bond_type_as_double().unwrap_or(0.0) > 1.0)
}

/// `Canon::details::atomHasFourthValence`
pub fn atom_has_fourth_valence(atom: &Atom) -> bool {
    atom.num_explicit_hs() == 1
        || (!atom.needs_update_property_cache()
            && atom.valence(ValenceType::Implicit).unwrap_or(0) == 1)
}

/// `Canon::chiralAtomNeedsTagInversion`
pub fn chiral_atom_needs_tag_inversion(
    mol: &RWMol,
    atom: u32,
    is_atom_first: bool,
    num_closures: usize,
) -> bool {
    let a = mol.atom(atom);
    mol.atom_degree(atom) == 3
        && ((is_atom_first && a.num_explicit_hs() == 1)
            || (!atom_has_fourth_valence(a) && num_closures == 1 && !is_unsaturated(mol, atom)))
}

/// `AdjustAtomChiralityFlags`
pub fn adjust_atom_chirality_flags(mol: &mut RWMol) {
    for i in 0..mol.num_atoms() {
        let tag = mol.atom(i).chiral_tag();
        match tag {
            ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW => {
                let (ordering, num_closures) = get_bond_ordering(mol, i);
                let mut n_swaps = mol.perturbation_order(i, &ordering).unwrap_or(0);
                if chiral_atom_needs_tag_inversion(
                    mol,
                    i,
                    mol.atom(i).has_prop(cp::_SmilesStart),
                    num_closures,
                ) {
                    n_swaps += 1;
                }
                if n_swaps % 2 != 0 {
                    mol.atom_mut(i).invert_chirality();
                }
            }
            ChiralType::SquarePlanar | ChiralType::TrigonalBipyramidal | ChiralType::Octahedral => {
                let (mut bonds, _) = get_bond_ordering(mol, i);
                let ref_max = nontetrahedral::get_max_nbors(tag) as usize;
                if bonds.len() < ref_max {
                    let missing = ref_max - bonds.len();
                    let pos = if mol.atom(i).has_prop(cp::_SmilesStart) {
                        0
                    } else {
                        1.min(bonds.len())
                    };
                    for _ in 0..missing {
                        bonds.insert(pos, -1);
                    }
                }
                let perm = nontetrahedral::get_chiral_permutation(mol, i, &bonds, true);
                mol.atom_mut(i).props.set_val(cp::_chiralPermutation, perm);
            }
            _ => {}
        }
    }
}

/// `MolOps::getAttachmentPointLabelNumber`
pub fn get_attachment_point_label_number(mol: &RWMol, atom: u32) -> u32 {
    let a = mol.atom(atom);
    if a.atomic_num() != 0 || mol.atom_degree(atom) != 1 {
        return 0;
    }
    let Some(RDValue::String(label)) = a.props.get(cp::atomLabel) else {
        return 0;
    };
    let Some(suffix) = label.strip_prefix("_AP") else {
        return 0;
    };
    if suffix.is_empty() || !suffix.bytes().all(|c| c.is_ascii_digit()) {
        return 0;
    }
    suffix.parse().unwrap_or(0)
}

/// `CleanupAfterParsing`
pub fn cleanup_after_parsing(mol: &mut RWMol) {
    for i in 0..mol.num_atoms() {
        let ap = get_attachment_point_label_number(mol, i);
        let a = mol.atom_mut(i);
        a.clear_prop(cp::_RingClosures);
        a.clear_prop(cp::_SmilesStart);
        if ap == 1 || ap == 2 {
            a.props.set_val(cp::_fromAttachPoint, ap);
        }
    }
    for i in 0..mol.num_bonds() {
        let b = mol.bond_mut(i);
        b.clear_prop(cp::_unspecifiedOrder);
        b.clear_prop("_cxsmilesBondIdx");
    }
    if !get_allow_nontetrahedral_chirality() {
        let mut need_warn = false;
        for i in 0..mol.num_atoms() {
            let a = mol.atom_mut(i);
            if a.has_prop(cp::_chiralPermutation) {
                need_warn = true;
                a.clear_prop(cp::_chiralPermutation);
            }
            if a.chiral_tag() > ChiralType::Other {
                need_warn = true;
                a.set_chiral_tag(ChiralType::Unspecified);
            }
        }
        if need_warn {
            crate::rd_log_warning!(
                "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
            );
        }
    }
}
