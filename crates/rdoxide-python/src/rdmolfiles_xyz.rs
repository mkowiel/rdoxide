//! `MolFromXYZBlock`, `MolFromXYZFile`, `MolToXYZBlock` and `MolToXYZFile`
//! (upstream `Wrap/rdmolfiles.cpp`).

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::xyz::{
    mol_from_xyz_block, mol_from_xyz_file, mol_to_xyz_block, mol_to_xyz_file,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// The upstream wrappers log `FileParseException`s as warnings and swallow
/// every other error, returning `None`.
fn finish(res: rdoxide::Result<Option<rdoxide::graphmol::RWMol>>) -> Option<PyMol> {
    match res {
        Ok(mol) => mol.map(|mol| PyMol { mol }),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            None
        }
        Err(_) => None,
    }
}

/// Construct a molecule from an XYZ file.
///
///   ARGUMENTS:
///
///     - xyzFileName: name of the file to read
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromXYZFile", signature = (xyzFileName))]
#[allow(non_snake_case)]
fn py_mol_from_xyz_file(py: Python<'_>, xyzFileName: &str) -> Option<PyMol> {
    finish(py.detach(|| mol_from_xyz_file(xyzFileName)))
}

/// Construct a molecule from an XYZ string.
///
///   ARGUMENTS:
///
///     - xyzBlock: the XYZ data to read
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromXYZBlock", signature = (xyzBlock))]
#[allow(non_snake_case)]
fn py_mol_from_xyz_block(py: Python<'_>, xyzBlock: &Bound<'_, PyAny>) -> PyResult<Option<PyMol>> {
    // pyObjectToString accepts str or bytes
    let text = match xyzBlock.extract::<String>() {
        Ok(s) => s,
        Err(_) => String::from_utf8_lossy(&xyzBlock.extract::<Vec<u8>>()?).into_owned(),
    };
    Ok(finish(py.detach(|| mol_from_xyz_block(&text))))
}

/// Returns a XYZ block for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - precision: precision of the coordinates
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToXYZBlock", signature = (mol, confId = -1, precision = 6))]
#[allow(non_snake_case)]
fn py_mol_to_xyz_block(mol: PyRef<'_, PyMol>, confId: i32, precision: usize) -> PyResult<String> {
    mol_to_xyz_block(&mol.mol, confId, precision).py()
}

/// Writes a XYZ file for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - filename: the file to write to
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - precision: precision of the coordinates
#[pyfunction]
#[pyo3(name = "MolToXYZFile", signature = (mol, filename, confId = -1, precision = 6))]
#[allow(non_snake_case)]
fn py_mol_to_xyz_file(
    mol: PyRef<'_, PyMol>,
    filename: &str,
    confId: i32,
    precision: usize,
) -> PyResult<()> {
    mol_to_xyz_file(&mol.mol, filename, confId, precision).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_from_xyz_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_from_xyz_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_xyz_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_xyz_file, m)?)?;
    Ok(())
}
