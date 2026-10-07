//! Ported from RDKit `Code/GraphMol/StereoGroup.{h,cpp}`.
//!
//! Copyright (C) 2018-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Upstream stores `Atom*`/`Bond*`; here groups store atom/bond indices which
//! the molecule renumbers when atoms or bonds are removed.

use crate::{RDKitError, Result};

/// `StereoGroupType`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StereoGroupType {
    Absolute = 0,
    Or = 1,
    And = 2,
}

impl StereoGroupType {
    pub fn name(self) -> &'static str {
        match self {
            StereoGroupType::Absolute => "STEREO_ABSOLUTE",
            StereoGroupType::Or => "STEREO_OR",
            StereoGroupType::And => "STEREO_AND",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    pub(crate) group_type: StereoGroupType,
    pub(crate) atoms: Vec<u32>,
    pub(crate) bonds: Vec<u32>,
    pub(crate) read_id: u32,
    pub(crate) write_id: u32,
}

fn check_for_dupes(v: &[u32], what: &str) -> Result<()> {
    for (i, x) in v.iter().enumerate() {
        if v[..i].contains(x) {
            return Err(RDKitError::Value(format!(
                "Duplicate {what} in StereoGroup"
            )));
        }
    }
    Ok(())
}

impl StereoGroup {
    pub fn new(
        group_type: StereoGroupType,
        atoms: Vec<u32>,
        bonds: Vec<u32>,
        read_id: u32,
    ) -> Result<Self> {
        check_for_dupes(&atoms, "atom")?;
        check_for_dupes(&bonds, "bond")?;
        Ok(StereoGroup {
            group_type,
            atoms,
            bonds,
            read_id,
            write_id: 0,
        })
    }
    pub fn group_type(&self) -> StereoGroupType {
        self.group_type
    }
    pub fn atoms(&self) -> &[u32] {
        &self.atoms
    }
    pub fn bonds(&self) -> &[u32] {
        &self.bonds
    }
    pub fn read_id(&self) -> u32 {
        self.read_id
    }
    pub fn write_id(&self) -> u32 {
        self.write_id
    }
    pub fn set_write_id(&mut self, id: u32) {
        self.write_id = id;
    }
}

/// `removeAtomFromGroups` + renumbering of higher atom indices.
pub(crate) fn remove_atom_from_groups(atom: u32, groups: &mut Vec<StereoGroup>) {
    for g in groups.iter_mut() {
        if let Some(pos) = g.atoms.iter().position(|&a| a == atom) {
            g.atoms.remove(pos);
        }
    }
    groups.retain(|g| !(g.atoms.is_empty() && g.bonds.is_empty()));
}

/// `removeBondFromGroups`
pub(crate) fn remove_bond_from_groups(bond: u32, groups: &mut Vec<StereoGroup>) {
    for g in groups.iter_mut() {
        if let Some(pos) = g.bonds.iter().position(|&b| b == bond) {
            g.bonds.remove(pos);
        }
    }
    groups.retain(|g| !(g.atoms.is_empty() && g.bonds.is_empty()));
}

/// `removeGroupsWithAtom`
pub fn remove_groups_with_atom(atom: u32, groups: &mut Vec<StereoGroup>) {
    groups.retain(|g| !g.atoms.contains(&atom));
}

/// `removeGroupsWithBond`
pub fn remove_groups_with_bond(bond: u32, groups: &mut Vec<StereoGroup>) {
    groups.retain(|g| !g.bonds.contains(&bond));
}

/// `assignStereoGroupIds`: give every group without a write id a unique one.
pub fn assign_stereo_group_ids(groups: &mut [StereoGroup]) {
    if groups.is_empty() {
        return;
    }
    let mut and_ids: Vec<bool> = Vec::new();
    let mut or_ids: Vec<bool> = Vec::new();
    for g in groups.iter_mut() {
        let ids = match g.group_type {
            StereoGroupType::And => &mut and_ids,
            StereoGroupType::Or => &mut or_ids,
            StereoGroupType::Absolute => continue,
        };
        let id = g.write_id as usize;
        if id == 0 {
            continue;
        }
        if id >= ids.len() {
            ids.resize(id + 1, false);
        }
        if ids[id] {
            crate::rd_log_warning!(
                "StereoGroup ID {} is used by more than one group, and will be reassined",
                id
            );
            g.write_id = 0;
        } else {
            ids[id] = true;
        }
    }
    let (mut next_and, mut next_or) = (0usize, 0usize);
    for g in groups.iter_mut() {
        let (ids, next) = match g.group_type {
            StereoGroupType::And => (&and_ids, &mut next_and),
            StereoGroupType::Or => (&or_ids, &mut next_or),
            StereoGroupType::Absolute => continue,
        };
        if g.write_id == 0 {
            *next += 1;
            while *next < ids.len() && ids[*next] {
                *next += 1;
            }
            g.write_id = *next as u32;
        }
    }
}
