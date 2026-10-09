//! `SubstanceGroup`: collections of atoms and bonds with properties, as in
//! the MDL formats (polymers, data fields, mixtures, ...).
//!
//! Ported from RDKit `Code/GraphMol/SubstanceGroup.{h,cpp}` (2026.09.1).
//!
//! Copyright (C) 2018-2020 Greg Landrum and T5 Informatics GmbH.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Upstream groups keep a pointer to the molecule that owns them; here the
//! molecule owns its groups, and the operations that need to look at the
//! molecule (index validation, bookmarks, bond classification) take it as an
//! argument.

use std::collections::BTreeMap;

use crate::geometry::Point3D;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::{Dict, RDProps, RDValue};
use crate::{RDKitError, Result};

const SOURCE: &str = "Code/GraphMol/SubstanceGroup.cpp";

/// Whether a group bond crosses the group boundary (`SubstanceGroup::BondType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SGroupBondType {
    /// external/crossing bond
    XBond,
    /// internal/contained bond
    CBond,
}

/// A bracket: three points (the third is for 3D-aware formats).
pub type Bracket = [Point3D; 3];

/// A SAP line (`SubstanceGroup::AttachPoint`); `lv_idx` is -1 when not set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttachPoint {
    pub a_idx: u32,
    pub lv_idx: i32,
    pub id: String,
}

/// V3000 `CSTATE`; the vector is ignored when comparing.
#[derive(Debug, Clone)]
pub struct CState {
    pub bond_idx: u32,
    pub vector: Point3D,
}

impl PartialEq for CState {
    fn eq(&self, other: &Self) -> bool {
        self.bond_idx == other.bond_idx
    }
}

/// `SubstanceGroup`
#[derive(Debug, Clone)]
pub struct SubstanceGroup {
    props: Dict,
    is_valid: bool,
    atoms: Vec<u32>,
    parent_atoms: Vec<u32>,
    bonds: Vec<u32>,
    brackets: Vec<Bracket>,
    cstates: Vec<CState>,
    attach_points: Vec<AttachPoint>,
}

impl RDProps for SubstanceGroup {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

/// Equality ignores brackets and cstates, which involve coordinates (and, as
/// upstream, the properties).
impl PartialEq for SubstanceGroup {
    fn eq(&self, other: &Self) -> bool {
        self.atoms == other.atoms
            && self.parent_atoms == other.parent_atoms
            && self.bonds == other.bonds
            && self.attach_points == other.attach_points
    }
}

fn remove_element(container: &mut Vec<u32>, element: u32) {
    if let Some(pos) = container.iter().position(|&e| e == element) {
        container.remove(pos);
    }
}

fn check_atom_range(mol: &ROMol, atoms: &[u32]) -> Result<()> {
    match atoms.iter().find(|&&a| a >= mol.num_atoms()) {
        Some(a) => Err(RDKitError::Value(format!("Atom index {a} is out of range"))),
        None => Ok(()),
    }
}

fn group_error(message: impl Into<String>) -> RDKitError {
    // `SubstanceGroupException` is a `std::runtime_error`
    RDKitError::Runtime(message.into())
}

impl SubstanceGroup {
    /// `SubstanceGroup(mol, type)`: the type is required, other properties
    /// depend on it.
    pub fn new(sgroup_type: &str) -> Self {
        let mut props = Dict::default();
        props.set_val("TYPE", sgroup_type.to_string());
        SubstanceGroup {
            props,
            is_valid: true,
            atoms: Vec::new(),
            parent_atoms: Vec::new(),
            bonds: Vec::new(),
            brackets: Vec::new(),
            cstates: Vec::new(),
            attach_points: Vec::new(),
        }
    }

    /// invalid groups must be ignored
    pub fn is_valid(&self) -> bool {
        self.is_valid
    }
    pub fn set_is_valid(&mut self, valid: bool) {
        self.is_valid = valid;
    }

    pub fn atoms(&self) -> &[u32] {
        &self.atoms
    }
    pub fn parent_atoms(&self) -> &[u32] {
        &self.parent_atoms
    }
    pub fn bonds(&self) -> &[u32] {
        &self.bonds
    }
    pub fn brackets(&self) -> &[Bracket] {
        &self.brackets
    }
    pub fn brackets_mut(&mut self) -> &mut Vec<Bracket> {
        &mut self.brackets
    }
    pub fn cstates(&self) -> &[CState] {
        &self.cstates
    }
    pub fn cstates_mut(&mut self) -> &mut Vec<CState> {
        &mut self.cstates
    }
    pub fn attach_points(&self) -> &[AttachPoint] {
        &self.attach_points
    }
    pub fn attach_points_mut(&mut self) -> &mut Vec<AttachPoint> {
        &mut self.attach_points
    }

    pub fn clear_brackets(&mut self) {
        self.brackets.clear();
    }
    pub fn clear_cstates(&mut self) {
        self.cstates.clear();
    }
    pub fn clear_attach_points(&mut self) {
        self.attach_points.clear();
    }

    /// `setAtoms`: the indices must exist in `mol`
    pub fn set_atoms(&mut self, mol: &ROMol, atoms: Vec<u32>) -> Result<()> {
        check_atom_range(mol, &atoms)?;
        self.atoms = atoms;
        Ok(())
    }

    /// `setParentAtoms`
    pub fn set_parent_atoms(&mut self, mol: &ROMol, atoms: Vec<u32>) -> Result<()> {
        check_atom_range(mol, &atoms)?;
        self.parent_atoms = atoms;
        Ok(())
    }

    /// `setBonds`
    pub fn set_bonds(&mut self, mol: &ROMol, bonds: Vec<u32>) -> Result<()> {
        if let Some(b) = bonds.iter().find(|&&b| b >= mol.num_bonds()) {
            return Err(RDKitError::Value(format!("Bond index {b} is out of range")));
        }
        self.bonds = bonds;
        Ok(())
    }

    pub fn add_atom_with_idx(&mut self, mol: &ROMol, idx: u32) -> Result<()> {
        if idx >= mol.num_atoms() {
            return Err(RDKitError::Value("Atom index out of range".into()));
        }
        self.atoms.push(idx);
        Ok(())
    }

    pub fn add_atom_with_bookmark(&mut self, mol: &ROMol, mark: i32) -> Result<()> {
        self.atoms.push(mol.get_atom_with_bookmark(mark)?);
        Ok(())
    }

    pub fn add_parent_atom_with_idx(&mut self, mol: &ROMol, idx: u32) -> Result<()> {
        if idx >= mol.num_atoms() {
            return Err(RDKitError::Value("Atom index out of range".into()));
        }
        if !self.atoms.contains(&idx) {
            return Err(group_error(format!(
                "Atom {idx} is not a member of current SubstanceGroup"
            )));
        }
        self.parent_atoms.push(idx);
        Ok(())
    }

    pub fn add_parent_atom_with_bookmark(&mut self, mol: &ROMol, mark: i32) -> Result<()> {
        let idx = mol.get_atom_with_bookmark(mark)?;
        if !self.atoms.contains(&idx) {
            return Err(group_error(format!(
                "Atom with bookmark {mark} is not a member of current SubstanceGroup "
            )));
        }
        self.parent_atoms.push(idx);
        Ok(())
    }

    pub fn add_bond_with_idx(&mut self, mol: &ROMol, idx: u32) -> Result<()> {
        if idx >= mol.num_bonds() {
            return Err(RDKitError::Value("Bond index out of range".into()));
        }
        self.bonds.push(idx);
        Ok(())
    }

    pub fn add_bond_with_bookmark(&mut self, mol: &ROMol, mark: i32) -> Result<()> {
        self.bonds.push(mol.get_bond_with_bookmark(mark)?);
        Ok(())
    }

    // These can leave attachment points and CStates in an invalid state.
    pub fn remove_atom_with_idx(&mut self, idx: u32) {
        remove_element(&mut self.atoms, idx);
    }
    pub fn remove_parent_atom_with_idx(&mut self, idx: u32) {
        remove_element(&mut self.parent_atoms, idx);
    }
    pub fn remove_bond_with_idx(&mut self, idx: u32) {
        remove_element(&mut self.bonds, idx);
    }

    pub fn add_bracket(&mut self, bracket: Bracket) {
        self.brackets.push(bracket);
    }

    pub fn add_cstate(&mut self, mol: &ROMol, bond_idx: u32, vector: Point3D) -> Result<()> {
        if self.bonds.is_empty() {
            return Err(RDKitError::precondition_at(
                "no bonds",
                "!d_bonds.empty()",
                SOURCE,
                181,
            ));
        }
        if self.bond_type(mol, bond_idx)? != SGroupBondType::XBond {
            return Err(group_error(format!(
                "Bond with index {bond_idx} is not an XBOND for current SubstanceGroup"
            )));
        }
        self.cstates.push(CState { bond_idx, vector });
        Ok(())
    }

    pub fn add_attach_point(&mut self, a_idx: u32, lv_idx: i32, id: &str) {
        self.attach_points.push(AttachPoint {
            a_idx,
            lv_idx,
            id: id.to_string(),
        });
    }

    /// Is the bond a crossing (XBOND) or contained (CBOND) bond of the group?
    pub fn bond_type(&self, mol: &ROMol, bond_idx: u32) -> Result<SGroupBondType> {
        if !self.bonds.contains(&bond_idx) {
            return Err(RDKitError::precondition_at(
                "bond is not part of the SubstanceGroup",
                "std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end()",
                SOURCE,
                202,
            ));
        }
        let bond = mol.get_bond(bond_idx)?;
        let begin_in = self.atoms.contains(&bond.begin_atom_idx());
        let end_in = self.atoms.contains(&bond.end_atom_idx());
        match (begin_in, end_in) {
            (true, true) => Ok(SGroupBondType::CBond),
            (true, false) | (false, true) => Ok(SGroupBondType::XBond),
            (false, false) => Err(group_error(format!(
                "Neither beginning nor ending atoms of bond {bond_idx} is in this SubstanceGroup."
            ))),
        }
    }

    /// Adjusts the atom indices for the removal of atom `atom_idx` from the
    /// molecule: higher indices move down by one. Fails if the group
    /// includes the atom. Returns whether anything changed.
    pub fn adjust_to_removed_atom(&mut self, atom_idx: u32) -> Result<bool> {
        let contains = || {
            group_error("adjustToRemovedAtom() called on SubstanceGroup which contains the atom")
        };
        let mut changed = false;
        for aid in self.atoms.iter_mut().chain(self.parent_atoms.iter_mut()) {
            if *aid == atom_idx {
                return Err(contains());
            }
            if *aid > atom_idx {
                changed = true;
                *aid -= 1;
            }
        }
        for ap in &mut self.attach_points {
            if ap.a_idx == atom_idx || ap.lv_idx == atom_idx as i32 {
                return Err(contains());
            }
            if ap.a_idx > atom_idx {
                changed = true;
                ap.a_idx -= 1;
            }
            if ap.lv_idx > atom_idx as i32 {
                changed = true;
                ap.lv_idx -= 1;
            }
        }
        Ok(changed)
    }

    /// Is the atom part of the definition of this group?
    pub fn includes_atom(&self, atom_idx: u32) -> bool {
        self.atoms.contains(&atom_idx)
            || self.parent_atoms.contains(&atom_idx)
            || self
                .attach_points
                .iter()
                .any(|ap| ap.a_idx == atom_idx || ap.lv_idx == atom_idx as i32)
    }

    /// The bond version of [`Self::adjust_to_removed_atom`].
    pub fn adjust_to_removed_bond(&mut self, bond_idx: u32) -> Result<bool> {
        let contains = || {
            group_error("adjustToRemovedBond() called on SubstanceGroup which contains the bond")
        };
        let mut changed = false;
        for bid in &mut self.bonds {
            if *bid == bond_idx {
                return Err(contains());
            }
            if *bid > bond_idx {
                changed = true;
                *bid -= 1;
            }
        }
        for cs in &mut self.cstates {
            if cs.bond_idx == bond_idx {
                return Err(contains());
            }
            if cs.bond_idx > bond_idx {
                changed = true;
                cs.bond_idx -= 1;
            }
        }
        Ok(changed)
    }

    /// Is the bond part of the definition of this group?
    pub fn includes_bond(&self, bond_idx: u32) -> bool {
        self.bonds.contains(&bond_idx) || self.cstates.iter().any(|cs| cs.bond_idx == bond_idx)
    }

    /// The `PARENT` property (an index into the `index` properties).
    fn parent(&self) -> Option<u32> {
        self.props.get("PARENT").and_then(RDValue::as_uint)
    }

    fn index_prop(&self) -> Option<u32> {
        self.props.get("index").and_then(RDValue::as_uint)
    }

    /// The upstream `operator<<`: `<index> <TYPE> Brk: .. CSt: .. AtPt: ..`
    /// followed by the member lists.
    pub fn describe(&self, index_in_mol: usize) -> String {
        let mut s = format!(
            "{index_in_mol} {}",
            self.props
                .get("TYPE")
                .map(RDValue::to_rd_string)
                .unwrap_or_default()
        );
        if !self.brackets.is_empty() {
            s += &format!(" Brk: {}", self.brackets.len());
        }
        if !self.cstates.is_empty() {
            s += &format!(" CSt: {}", self.cstates.len());
        }
        if !self.attach_points.is_empty() {
            s += &format!(" AtPt: {}", self.attach_points.len());
        }
        for (label, list) in [
            ("Atoms", &self.atoms),
            ("PAtoms", &self.parent_atoms),
            ("Bonds", &self.bonds),
        ] {
            if !list.is_empty() {
                s += &format!(" {label}: {{ ");
                for i in list {
                    s += &format!("{i} ");
                }
                s += "}";
            }
        }
        s
    }
}

/// `SubstanceGroupChecks`
pub mod checks {
    use super::*;

    pub const SGROUP_TYPES: [&str; 15] = [
        // polymer sgroups:
        "SRU", "MON", "COP", "CRO", "GRA", "MOD", "MER", "ANY", // formulations/mixtures:
        "COM", "MIX", "FOR", // other
        "SUP", "MUL", "DAT", "GEN",
    ];
    pub const SGROUP_SUBTYPES: [&str; 3] = ["ALT", "RAN", "BLO"];
    pub const SGROUP_CONNECT_TYPES: [&str; 3] = ["HH", "HT", "EU"];
    pub const SGROUP_CLASSES: [&str; 18] = [
        "AA",
        "dAA",
        "DNA",
        "RNA",
        "SUGAR",
        "BASE",
        "PHOSPHATE",
        "LINKER",
        "CHEM",
        "LGRP",
        "MODAA",
        "MODdAA",
        "MODDNA",
        "MODRNA",
        "XLINKAA",
        "XLINKdAA",
        "XLINKDNA",
        "XLINKRNA",
    ];

    pub fn is_valid_type(t: &str) -> bool {
        SGROUP_TYPES.contains(&t)
    }
    pub fn is_valid_sub_type(t: &str) -> bool {
        SGROUP_SUBTYPES.contains(&t)
    }
    pub fn is_valid_connect_type(t: &str) -> bool {
        SGROUP_CONNECT_TYPES.contains(&t)
    }
    pub fn is_valid_class(c: &str) -> bool {
        SGROUP_CLASSES.contains(&c)
    }

    /// No group of `mol` uses `id` as its `ID` property.
    pub fn is_substance_group_id_free(mol: &ROMol, id: u32) -> bool {
        !mol.substance_groups()
            .iter()
            .any(|sg| sg.props.get("ID").and_then(RDValue::as_uint) == Some(id))
    }
}

/// `addSubstanceGroup`: adds the group to the molecule and returns its index.
pub fn add_substance_group(mol: &mut ROMol, sgroup: SubstanceGroup) -> u32 {
    let groups = mol.substance_groups_mut();
    groups.push(sgroup);
    groups.len() as u32 - 1
}

/// Which kind of molecule component is being removed.
trait Removed {
    fn included_in(sg: &SubstanceGroup, idx: u32) -> bool;
    fn adjust(sg: &mut SubstanceGroup, idx: u32) -> Result<bool>;
}

struct RemovedAtom;
struct RemovedBond;

impl Removed for RemovedAtom {
    fn included_in(sg: &SubstanceGroup, idx: u32) -> bool {
        sg.includes_atom(idx)
    }
    fn adjust(sg: &mut SubstanceGroup, idx: u32) -> Result<bool> {
        sg.adjust_to_removed_atom(idx)
    }
}

impl Removed for RemovedBond {
    fn included_in(sg: &SubstanceGroup, idx: u32) -> bool {
        sg.includes_bond(idx)
    }
    fn adjust(sg: &mut SubstanceGroup, idx: u32) -> Result<bool> {
        sg.adjust_to_removed_bond(idx)
    }
}

/// Is `idx` or one of its ancestors (through `PARENT`) being removed?
fn removed_parent_in_hierarchy(
    idx: usize,
    groups: &[SubstanceGroup],
    to_remove: &[bool],
    index_lookup: &BTreeMap<u32, usize>,
) -> bool {
    if to_remove[idx] {
        return true;
    }
    groups[idx]
        .parent()
        .and_then(|parent| index_lookup.get(&parent))
        .is_some_and(|&p| removed_parent_in_hierarchy(p, groups, to_remove, index_lookup))
}

/// Drops the groups that include component `idx` (and their descendants) and
/// renumbers the other groups for its removal.
fn remove_substance_groups_referencing<R: Removed>(mol: &mut RWMol, idx: u32) -> Result<()> {
    let groups = std::mem::take(mol.substance_groups_mut());
    if groups.is_empty() {
        return Ok(());
    }
    let to_remove: Vec<bool> = groups.iter().map(|sg| R::included_in(sg, idx)).collect();
    let parents_present = groups.iter().any(|sg| sg.has_prop("PARENT"));
    // lookup from the "index" property to the position in the original list
    let mut index_lookup = BTreeMap::new();
    if parents_present && to_remove.iter().any(|&r| r) {
        for (i, sg) in groups.iter().enumerate() {
            if let Some(index) = sg.index_prop() {
                index_lookup.insert(index, i);
            }
        }
    }
    let keep: Vec<bool> = groups
        .iter()
        .enumerate()
        .map(|(i, sg)| {
            if to_remove[i] {
                return false;
            }
            let Some(parent) = sg.parent().filter(|_| parents_present) else {
                return !sg.has_prop("PARENT") || !parents_present;
            };
            match index_lookup.get(&parent) {
                // our parent isn't around, so it isn't being removed
                None => true,
                Some(&p) => {
                    !to_remove[p]
                        && !removed_parent_in_hierarchy(p, &groups, &to_remove, &index_lookup)
                }
            }
        })
        .collect();
    let mut kept = Vec::with_capacity(groups.len());
    for (mut sg, keep) in groups.into_iter().zip(keep) {
        if keep {
            R::adjust(&mut sg, idx)?;
            kept.push(sg);
        }
    }
    *mol.substance_groups_mut() = kept;
    Ok(())
}

/// `removeSubstanceGroupsReferencingAtom`: also renumbers the remaining
/// groups for the removal of the atom. (Deleting groups that contain the atom
/// is the safe choice, since it is unknown whether keeping them is fine.)
pub fn remove_substance_groups_referencing_atom(mol: &mut RWMol, idx: u32) -> Result<()> {
    remove_substance_groups_referencing::<RemovedAtom>(mol, idx)
}

/// `removeSubstanceGroupsReferencingBond`
pub fn remove_substance_groups_referencing_bond(mol: &mut RWMol, idx: u32) -> Result<()> {
    remove_substance_groups_referencing::<RemovedBond>(mol, idx)
}

/// Renumbers the groups of `mol` after `other` was inserted behind
/// `orig_num_atoms` atoms and `orig_num_bonds` bonds (`insertSubstanceGroups`).
pub(crate) fn insert_substance_groups(
    mol: &mut RWMol,
    other: &ROMol,
    orig_num_atoms: u32,
    orig_num_bonds: u32,
) -> Result<()> {
    for sgroup in other.substance_groups() {
        let mut sgroup = sgroup.clone();
        let shift = |v: &[u32], by: u32| v.iter().map(|i| i + by).collect::<Vec<_>>();
        let atoms = shift(&sgroup.atoms, orig_num_atoms);
        sgroup.set_atoms(mol, atoms)?;
        let bonds = shift(&sgroup.bonds, orig_num_bonds);
        sgroup.set_bonds(mol, bonds)?;
        let parent_atoms = shift(&sgroup.parent_atoms, orig_num_atoms);
        sgroup.set_parent_atoms(mol, parent_atoms)?;
        for cstate in &mut sgroup.cstates {
            cstate.bond_idx += orig_num_bonds;
        }
        for sap in &mut sgroup.attach_points {
            sap.a_idx += orig_num_atoms;
            if sap.lv_idx != -1 {
                sap.lv_idx += orig_num_atoms as i32;
            }
        }
        add_substance_group(mol, sgroup);
    }
    Ok(())
}
