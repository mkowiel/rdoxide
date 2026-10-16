//! `rdkit.Chem.rdchem.MolBundle` (upstream `GraphMol/Wrap/MolBundle.cpp`): a
//! group of related molecules. Not supported: pickles, the substructure
//! search methods and the properties.

use std::sync::Arc;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use rdoxide::graphmol::mol_bundle::MolBundle;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// A class for storing groups of related molecules.
#[pyclass(name = "MolBundle", module = "rdkit.Chem.rdchem")]
#[derive(Default)]
pub struct PyMolBundle {
    pub bundle: MolBundle,
}

impl PyMolBundle {
    pub fn wrap(bundle: MolBundle) -> Self {
        PyMolBundle { bundle }
    }

    fn get(&self, idx: i64) -> PyResult<PyMol> {
        let mol = usize::try_from(idx)
            .ok()
            .and_then(|i| self.bundle.mol(i).ok())
            .ok_or_else(|| PyIndexError::new_err(idx.to_string()))?;
        Ok(PyMol {
            mol: mol.as_ref().clone(),
        })
    }
}

#[pymethods]
impl PyMolBundle {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn __getitem__(&self, idx: i64) -> PyResult<PyMol> {
        self.get(idx)
    }

    fn __len__(&self) -> usize {
        self.bundle.size()
    }

    /// Adds a molecule (a copy of it is stored); returns the new size.
    #[pyo3(name = "AddMol")]
    fn add_mol(&mut self, nmol: PyRef<'_, PyMol>) -> PyResult<usize> {
        self.bundle.add_mol(Arc::new(nmol.mol.clone())).py()
    }

    #[pyo3(name = "GetMol")]
    fn get_mol(&self, idx: i64) -> PyResult<PyMol> {
        self.get(idx)
    }

    #[pyo3(name = "Size")]
    fn size(&self) -> usize {
        self.bundle.size()
    }
}
