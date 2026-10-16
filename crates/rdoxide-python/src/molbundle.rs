//! `rdoxide.Chem.rdchem.MolBundle` (upstream `GraphMol/Wrap/MolBundle.cpp`): a
//! group of related molecules. Not supported: pickles and the properties.

use std::sync::Arc;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::mol_bundle::MolBundle;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::substruct::{
    Molecules, get_substruct_match_impl, get_substruct_matches_impl, has_substruct_match_impl,
};

/// A class for storing groups of related molecules.
#[pyclass(name = "MolBundle", module = "rdoxide.Chem.rdchem", subclass)]
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

    /// Queries whether or not any molecule in the bundle contains a
    /// particular substructure (the query is a molecule or a MolBundle), or
    /// with SubstructMatchParameters.
    #[pyo3(name = "HasSubstructMatch", signature = (*args, **kwargs))]
    fn has_substruct_match<'py>(
        slf: &Bound<'py, Self>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<bool> {
        has_substruct_match_impl(&Molecules::Bundle(slf.clone()), args, kwargs)
    }

    /// Returns the indices of the atoms from the first molecule in the bundle
    /// that matches a substructure query (a molecule or a MolBundle).
    #[pyo3(name = "GetSubstructMatch", signature = (*args, **kwargs))]
    fn get_substruct_match<'py>(
        slf: &Bound<'py, Self>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        get_substruct_match_impl(&Molecules::Bundle(slf.clone()), args, kwargs)
    }

    /// Returns a tuple of all the indices of the atoms from the first
    /// molecule in the bundle that matches a substructure query (a molecule
    /// or a MolBundle).
    #[pyo3(name = "GetSubstructMatches", signature = (*args, **kwargs))]
    fn get_substruct_matches<'py>(
        slf: &Bound<'py, Self>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        get_substruct_matches_impl(&Molecules::Bundle(slf.clone()), args, kwargs)
    }
}

/// A MolBundle whose molecules all have the same numbers of atoms and bonds.
#[pyclass(
    name = "FixedMolSizeMolBundle",
    module = "rdoxide.Chem.rdchem",
    extends = PyMolBundle
)]
pub struct PyFixedMolSizeMolBundle;

#[pymethods]
impl PyFixedMolSizeMolBundle {
    #[new]
    fn new() -> PyClassInitializer<Self> {
        PyClassInitializer::from(PyMolBundle::wrap(MolBundle::new_fixed_size()))
            .add_subclass(PyFixedMolSizeMolBundle)
    }
}

/// Whether MolBundles can be serialized (pickled); rdoxide does not support
/// Boost serialization.
#[pyfunction]
#[pyo3(name = "MolBundleCanSerialize")]
pub fn mol_bundle_can_serialize() -> bool {
    false
}
