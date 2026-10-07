//! Ported from RDKit `Code/GraphMol/Conformer.{h,cpp}`.
//!
//! Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::geometry::Point3D;
use crate::rdgeneral::{Dict, RDProps};
use crate::{RDKitError, Result};

/// A set of atomic coordinates (`RDKit::Conformer`).
#[derive(Debug, Clone, PartialEq)]
pub struct Conformer {
    pub(crate) is_3d: bool,
    pub(crate) id: u32,
    pub(crate) has_owner: bool,
    pub(crate) positions: Vec<Point3D>,
    pub(crate) props: Dict,
}

impl Default for Conformer {
    fn default() -> Self {
        Conformer::new(0)
    }
}

impl RDProps for Conformer {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

impl Conformer {
    pub fn new(num_atoms: usize) -> Self {
        Conformer {
            is_3d: true,
            id: 0,
            has_owner: false,
            positions: vec![Point3D::default(); num_atoms],
            props: Dict::new(),
        }
    }

    /// Copy that is not owned by any molecule (copy-constructor semantics).
    pub fn copy_detached(&self) -> Conformer {
        let mut c = self.clone();
        c.has_owner = false;
        c
    }

    pub fn num_atoms(&self) -> u32 {
        self.positions.len() as u32
    }
    pub fn positions(&self) -> &[Point3D] {
        &self.positions
    }
    pub fn positions_mut(&mut self) -> &mut Vec<Point3D> {
        &mut self.positions
    }
    pub fn id(&self) -> u32 {
        self.id
    }
    pub fn set_id(&mut self, id: u32) {
        self.id = id;
    }
    pub fn is_3d(&self) -> bool {
        self.is_3d
    }
    pub fn set_3d(&mut self, v: bool) {
        self.is_3d = v;
    }
    pub fn has_owning_mol(&self) -> bool {
        self.has_owner
    }

    /// `getAtomPos` (range-checked like upstream's `URANGE_CHECK`).
    pub fn atom_pos(&self, idx: u32) -> Result<Point3D> {
        self.positions
            .get(idx as usize)
            .copied()
            .ok_or_else(|| range_error(idx, self.num_atoms()))
    }

    pub fn atom_pos_mut(&mut self, idx: u32) -> Result<&mut Point3D> {
        let n = self.num_atoms();
        self.positions
            .get_mut(idx as usize)
            .ok_or_else(|| range_error(idx, n))
    }

    /// `setAtomPos`: grows the conformer when needed.
    pub fn set_atom_pos(&mut self, idx: u32, p: Point3D) -> Result<()> {
        if idx == u32::MAX {
            return Err(RDKitError::Value("atom index overflow".into()));
        }
        if idx as usize >= self.positions.len() {
            self.positions.resize(idx as usize + 1, Point3D::default());
        }
        self.positions[idx as usize] = p;
        Ok(())
    }

    pub fn resize(&mut self, n: usize) {
        self.positions.resize(n, Point3D::default());
    }
}

/// `URANGE_CHECK` failure message.
pub(crate) fn range_error(idx: u32, max: u32) -> RDKitError {
    RDKitError::range("atomId", idx, max)
}

/// `hasNonZeroZCoords`
pub fn has_non_zero_z_coords(conf: &Conformer) -> bool {
    conf.positions.iter().any(|p| p.z.abs() > 1e-3)
}
