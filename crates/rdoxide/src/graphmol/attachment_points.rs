//! Explicit attachment points: dummy atoms marked as the places where a
//! molecule is attached to the rest of a structure, and their MDL `ATTCHPT`
//! form, an annotation on the atom they are attached to
//! (`expandAttachmentPoints`, `collapseAttachmentPoints`,
//! `isMarkedAttachmentPoint` and friends of `MolOps.cpp`).
//!
//! Ported from RDKit `Code/GraphMol/MolOps.cpp` (2026.09.1).
//!
//! Copyright (C) 2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::atom::Atom;
use super::bond::{BondDir, BondType};
use super::mol_ops::add_hs::set_terminal_atom_coords;
use super::queries::make_atom_null_query;
use super::romol::{ROMol, RWMol};
use super::smiles_parse::ops::get_attachment_point_label_number;
use crate::Result;
use crate::rd_log_warning;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};

/// `MolOps::attachmentPointLabelPrefix`: the prefix of numbered explicit
/// attachment-point atom labels.
pub const ATTACHMENT_POINT_LABEL_PREFIX: &str = "_AP";

fn int_prop(mol: &ROMol, atom: u32, key: &str) -> Option<i32> {
    mol.atom(atom).get_prop(key).and_then(RDValue::as_int)
}

/// `MolOps::isMarkedAttachmentPoint`: a degree-one dummy atom with the
/// `_fromAttachPoint` property or a valid `_AP<n>` label. This is about the
/// identity of the atom only, not about whether it can be collapsed.
pub fn is_marked_attachment_point(mol: &ROMol, atom: u32) -> bool {
    let a = mol.atom(atom);
    a.atomic_num() == 0
        && mol.atom_degree(atom) == 1
        && (a.has_prop(cp::_fromAttachPoint) || get_attachment_point_label_number(mol, atom) != 0)
}

/// `details::isAttachmentPoint`: whether a dummy atom could be collapsed.
pub fn is_attachment_point(mol: &ROMol, atom: u32, marked_only: bool) -> bool {
    let a = mol.atom(atom);
    if a.atomic_num() != 0 || mol.atom_degree(atom) != 1 {
        return false;
    }
    if marked_only && !is_marked_attachment_point(mol, atom) {
        return false;
    }
    // we know that the atom is degree 1
    let bond = mol.atom_bonds(atom).next().expect("a degree one atom");
    if !matches!(bond.bond_type(), BondType::Single | BondType::Unspecified)
        || bond.bond_dir() != BondDir::None
    {
        return false;
    }
    match a.query() {
        // a * from SMARTS
        Some(q) if !q.negate && q.description() == "AtomNull" => true,
        // a * from CXSMILES
        Some(q) => q.negate && q.description() == "AtomAtomicNum" && q.val == 1,
        None => true,
    }
}

/// `details::addExplicitAttachmentPoint`: adds a dummy atom marked as the
/// attachment point `val` to an atom and returns its index.
fn add_explicit_attachment_point(
    mol: &mut RWMol,
    atom_idx: u32,
    val: u32,
    add_as_query: bool,
    add_coords: bool,
) -> Result<u32> {
    let mut dummy = Atom::new(0);
    if add_as_query {
        dummy.set_query(Some(make_atom_null_query()));
    }
    dummy.props.set_val(cp::_fromAttachPoint, val);
    let idx = mol.add_atom(dummy, false);
    mol.add_bond(atom_idx, idx, BondType::Single)?;
    mol.update_atom_property_cache(idx, false)?;
    if add_coords && mol.num_conformers() > 0 {
        set_terminal_atom_coords(mol, idx, atom_idx)?;
    }
    Ok(idx)
}

/// `MolOps::expandAttachmentPoints`: attachment points encoded as `ATTCHPT`
/// properties are added to the graph as dummy atoms.
pub fn expand_attachment_points(
    mol: &mut RWMol,
    add_as_queries: bool,
    add_coords: bool,
) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        let Some(value) = int_prop(mol, idx, cp::molAttachPoint) else {
            continue;
        };
        let targets: Vec<u32> = [
            (1, value == 1 || value == -1),
            (2, value == 2 || value == -1),
        ]
        .into_iter()
        .filter_map(|(t, wanted)| wanted.then_some(t))
        .collect();
        if targets.is_empty() {
            rd_log_warning!(
                "Invalid value for molAttachPoint: {value} on atom {idx}. Not expanding this attachment point."
            );
            continue;
        }
        for tval in targets {
            mol.atom_mut(idx).clear_prop(cp::molAttachPoint);
            add_explicit_attachment_point(mol, idx, tval, add_as_queries, add_coords)?;
        }
    }
    Ok(())
}

/// `MolOps::collapseAttachmentPoints`: dummy atoms are removed and replaced
/// with attachment point annotations on the atoms they were attached to.
pub fn collapse_attachment_points(mol: &mut RWMol, marked_only: bool) -> Result<()> {
    let mut removed_any = false;
    let mut attach_labels = vec![0i32; mol.num_atoms() as usize];

    for idx in 0..mol.num_atoms() {
        if !is_attachment_point(mol, idx, marked_only) {
            continue;
        }
        let native = int_prop(mol, idx, cp::_fromAttachPoint);
        let mut value = native.unwrap_or(0);
        if marked_only && !(0..=2).contains(&value) {
            rd_log_warning!(
                "Invalid value for _fromAttachPoint: {value} on atom {idx}. Not collapsing this atom"
            );
            continue;
        }
        if marked_only && native.is_none() {
            // _AP<n> labels identify explicit attachment points, but n is
            // not an MDL ATTCHPT position. Treat a label-only atom as
            // position 1.
            value = 1;
        }
        if !marked_only && value == 0 {
            value = 1;
        }
        let bond = mol.atom_bonds(idx).next().expect("a degree one atom");
        if !matches!(bond.bond_type(), BondType::Single | BondType::Unspecified)
            || bond.bond_dir() != BondDir::None
        {
            continue;
        }
        let other = bond.other_atom_idx(idx)?;
        if attach_labels[other as usize] != 0 {
            if attach_labels[other as usize] != -1 {
                value = -1;
            } else {
                rd_log_warning!(
                    "More than two attachment points on atom {other}. Attachment point {idx} will not be collapsed."
                );
                continue;
            }
        }
        if !removed_any {
            mol.begin_batch_edit()?;
            removed_any = true;
        }
        attach_labels[other as usize] = value;
        mol.remove_atom(idx)?;
    }
    // set the attachment point labels
    for (idx, &label) in attach_labels.iter().enumerate() {
        if label != 0 {
            mol.atom_mut(idx as u32)
                .props
                .set_val(cp::molAttachPoint, label);
        }
    }
    if removed_any {
        mol.commit_batch_edit()?;
    }
    Ok(())
}
