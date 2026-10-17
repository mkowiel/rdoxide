//! `MolToCMLBlock` and `MolToCMLFile` (upstream `Wrap/rdmolfiles.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::cml_writer::{mol_to_cml_block, mol_to_cml_file};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Writes a CML block for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - confId: (optional) selects which conformation to output
///     - kekulize: (optional) triggers kekulization of the molecule before it's written
#[pyfunction]
#[pyo3(name = "MolToCMLBlock", signature = (mol, confId = -1, kekulize = true))]
#[allow(non_snake_case)]
fn py_mol_to_cml_block(
    mol: PyRef<'_, PyMol>,
    confId: i32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] kekulize: bool,
) -> PyResult<String> {
    mol_to_cml_block(&mol.mol, confId, kekulize).py()
}

/// Writes a CML file for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - filename: the file to write to
///     - confId: (optional) selects which conformation to output
///     - kekulize: (optional) triggers kekulization of the molecule before it's written
#[pyfunction]
#[pyo3(name = "MolToCMLFile", signature = (mol, filename, confId = -1, kekulize = true))]
#[allow(non_snake_case)]
fn py_mol_to_cml_file(
    mol: PyRef<'_, PyMol>,
    filename: &str,
    confId: i32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] kekulize: bool,
) -> PyResult<()> {
    mol_to_cml_file(&mol.mol, filename, confId, kekulize).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_to_cml_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_cml_file, m)?)
}
