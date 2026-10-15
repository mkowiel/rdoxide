//! Labels and helpers shared by the R-group decomposition (upstream
//! `RGroupDecomposition/RGroupUtils.{h,cpp}`).

use std::collections::{BTreeMap, BTreeSet};

use super::params::rgroup_labelling as labelling;
use crate::graphmol::atom::Atom;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_warning};

/// Temporary R-label of a core atom.
pub const RLABEL: &str = "tempRlabel";
/// How the R-label was found, a [`Labelling`] as `i32`.
pub const RLABEL_TYPE: &str = "tempRlabelType";
/// Index of an atom in the core, kept on the atoms of the matching molecule.
pub const RLABEL_CORE_INDEX: &str = "rLabelCoreIndex";
/// R-labels of a side chain dummy atom (`Vec<i32>`).
pub const SIDECHAIN_RLABELS: &str = "sideChainRlabels";
/// Set on an R group molecule once it has been relabelled.
pub const DONE: &str = "RLABEL_PROCESSED";
/// Set on dummy atoms that came from the input molecule.
pub const RGROUP_INPUT_DUMMY: &str = "_rgroupInputDummy";
/// Set on the dummies added to unlabelled core atoms.
pub const UNLABELED_CORE_ATTACHMENT: &str = "unlabeledCoreAttachment";

/// Label of the core used before a core was found (`EMPTY_CORE_LABEL`).
pub const EMPTY_CORE_LABEL: i32 = -100000;

/// Where an R-label of a core came from (upstream `Labelling`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Labelling {
    RGroupLabels = 0,
    IsotopeLabels = 1,
    AtomMapLabels = 2,
    IndexLabels = 3,
    DummyLabels = 4,
    InternalLabels = 5,
}

impl Labelling {
    pub fn from_i32(v: i32) -> Option<Labelling> {
        Some(match v {
            0 => Labelling::RGroupLabels,
            1 => Labelling::IsotopeLabels,
            2 => Labelling::AtomMapLabels,
            3 => Labelling::IndexLabels,
            4 => Labelling::DummyLabels,
            5 => Labelling::InternalLabels,
            _ => return None,
        })
    }

    /// `labellingToString`
    pub fn name(self) -> &'static str {
        match self {
            Labelling::RGroupLabels => "RGroupLabels",
            Labelling::IsotopeLabels => "IsotopeLabels",
            Labelling::AtomMapLabels => "AtomMapLabels",
            Labelling::IndexLabels => "IndexLabels",
            Labelling::DummyLabels => "DummyLabels",
            Labelling::InternalLabels => "InternalLabels",
        }
    }
}

pub(crate) fn int_prop(atom: &Atom, key: &str) -> Option<i32> {
    atom.get_prop(key).and_then(RDValue::as_int)
}

pub(crate) fn set_int_prop(atom: &mut Atom, key: &str, value: i32) {
    atom.dict_mut().set_val(key, value);
}

/// `getRlabels`: the atoms of `mol` carrying an R-label, by label.
pub fn get_rlabels(mol: &RWMol) -> Result<BTreeMap<i32, u32>> {
    let mut atoms = BTreeMap::new();
    for (idx, atom) in mol.atoms().iter().enumerate() {
        if let Some(rlabel) = int_prop(atom, RLABEL) {
            if atoms.insert(rlabel, idx as u32).is_some() {
                return Err(RDKitError::invariant("Duplicate labels in rgroup core!"));
            }
        }
    }
    Ok(atoms)
}

/// `clearInputLabels`
pub fn clear_input_labels(atom: &mut Atom) {
    let _ = atom.set_atom_map_num(0, false);
    atom.clear_prop(cp::_MolFileRLabel);
}

/// `setLabel`: gives `atom` the R-label `label` unless it is 0. Duplicated
/// labels are replaced by `max_label` when `relabel` is set, else rejected.
pub fn set_label(
    atom: &mut Atom,
    mut label: i32,
    labels: &mut BTreeSet<i32>,
    max_label: &mut i32,
    relabel: bool,
    ty: Labelling,
) -> Result<bool> {
    match ty {
        Labelling::IsotopeLabels => atom.set_isotope(0),
        Labelling::AtomMapLabels => {
            let _ = atom.set_atom_map_num(0, false);
        }
        Labelling::RGroupLabels if atom.has_prop(cp::_MolFileRLabel) => {
            atom.clear_prop(cp::_MolFileRLabel);
            atom.set_isotope(0);
        }
        _ => {}
    }
    if label == 0 {
        return Ok(false);
    }
    if labels.contains(&label) {
        if !relabel {
            return Err(RDKitError::Value(format!(
                "Duplicate label in input, current type is:{}",
                ty.name()
            )));
        }
        if ty == Labelling::InternalLabels {
            rd_log_warning!("Relabelling existing label");
        }
        label = *max_label;
    }
    set_int_prop(atom, RLABEL, label);
    set_int_prop(atom, RLABEL_TYPE, ty as i32);
    labels.insert(label);
    *max_label = (*max_label).max(label + 1);
    Ok(true)
}

/// `isUserRLabel`
pub fn is_user_rlabel(atom: &Atom) -> bool {
    atom.has_prop(RLABEL)
        && int_prop(atom, RLABEL_TYPE)
            .and_then(Labelling::from_i32)
            .is_some_and(|t| t != Labelling::IndexLabels)
}

/// `isDummyRGroupAttachment`
pub fn is_dummy_rgroup_attachment(mol: &ROMol, idx: u32) -> bool {
    let atom = mol.atom(idx);
    if atom.atomic_num() != 0 || mol.atom_degree(idx) != 1 {
        return false;
    }
    if is_user_rlabel(atom) {
        return true;
    }
    atom.get_prop(UNLABELED_CORE_ATTACHMENT)
        .and_then(RDValue::as_bool)
        .unwrap_or(false)
}

/// `isAtomWithMultipleNeighborsOrNotDummyRGroupAttachment`
pub fn is_atom_with_multiple_neighbors_or_not_dummy_rgroup_attachment(
    mol: &ROMol,
    idx: u32,
) -> bool {
    mol.atom_degree(idx) > 1 || !is_dummy_rgroup_attachment(mol, idx)
}

/// `isAnyAtomWithMultipleNeighborsOrNotUserRLabel`
pub fn is_any_atom_with_multiple_neighbors_or_not_user_rlabel(mol: &ROMol, idx: u32) -> bool {
    if mol.atom(idx).atomic_num() != 0 {
        return false;
    }
    is_atom_with_multiple_neighbors_or_not_dummy_rgroup_attachment(mol, idx)
}

/// `relabelMappedDummies`: relabels the dummy atoms carrying an R-group
/// mapping so that renderers show them as `R#`.
pub fn relabel_mapped_dummies(mol: &mut ROMol, input_labels: u32, output_labels: u32) {
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom_mut(idx);
        if atom.atomic_num() != 0 {
            continue;
        }
        let mut atom_map_num = 0u32;
        if input_labels & labelling::ATOM_MAP != 0 {
            atom_map_num = atom.atom_map_num().unsigned_abs();
        }
        if atom_map_num == 0 && input_labels & labelling::ISOTOPE != 0 {
            atom_map_num = atom.isotope();
        }
        if atom_map_num == 0 && input_labels & labelling::MDL_RGROUP != 0 {
            atom_map_num = atom
                .get_prop(cp::_MolFileRLabel)
                .and_then(RDValue::as_uint)
                .unwrap_or(0);
        }
        if atom_map_num == 0 {
            continue;
        }
        let r_label = format!("R{atom_map_num}");
        let _ = atom.set_atom_map_num(
            if output_labels & labelling::ATOM_MAP != 0 {
                atom_map_num as i32
            } else {
                0
            },
            false,
        );
        atom.set_isotope(if output_labels & labelling::ISOTOPE != 0 {
            atom_map_num
        } else {
            0
        });
        if output_labels & labelling::MDL_RGROUP != 0 {
            atom.dict_mut().set_val(cp::_MolFileRLabel, atom_map_num);
            atom.dict_mut().set_val(cp::dummyLabel, r_label);
        } else {
            atom.clear_prop(cp::_MolFileRLabel);
            atom.clear_prop(cp::dummyLabel);
        }
    }
}
