//! Port of RDKit `Code/GraphMol/MolBundle.h`: groups of related molecules.
//! Boost-serialization pickles are not supported.

use std::sync::Arc;

use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{Dict, RDProps};
use crate::{RDKitError, Result};

/// `MolBundle`; with `fixed_size` it is a `FixedMolSizeMolBundle`, whose
/// molecules all have the same numbers of atoms and bonds.
#[derive(Debug, Clone, Default)]
pub struct MolBundle {
    props: Dict,
    mols: Vec<Arc<ROMol>>,
    fixed_size: bool,
}

impl RDProps for MolBundle {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

impl MolBundle {
    pub fn new() -> Self {
        Self::default()
    }

    /// `FixedMolSizeMolBundle`
    pub fn new_fixed_size() -> Self {
        MolBundle {
            fixed_size: true,
            ..Default::default()
        }
    }

    pub fn mols(&self) -> &[Arc<ROMol>] {
        &self.mols
    }

    /// `addMol`: returns the new size.
    pub fn add_mol(&mut self, mol: Arc<ROMol>) -> Result<usize> {
        if self.fixed_size
            && let Some(first) = self.mols.first()
        {
            if mol.num_atoms() != first.num_atoms() {
                return Err(RDKitError::Value(
                    "all molecules in a bundle must have the same number of atoms".into(),
                ));
            }
            if mol.num_bonds() != first.num_bonds() {
                return Err(RDKitError::Value(
                    "all molecules in a bundle must have the same number of bonds".into(),
                ));
            }
        }
        self.mols.push(mol);
        Ok(self.mols.len())
    }

    pub fn size(&self) -> usize {
        self.mols.len()
    }

    pub fn is_empty(&self) -> bool {
        self.mols.is_empty()
    }

    /// `getMol`
    pub fn mol(&self, idx: usize) -> Result<&Arc<ROMol>> {
        self.mols
            .get(idx)
            .ok_or_else(|| RDKitError::Index(format!("{idx}")))
    }
}
