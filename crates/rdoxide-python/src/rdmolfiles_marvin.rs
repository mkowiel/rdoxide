//! `MolFromMrvFile` and `MolFromMrvBlock` (upstream `Wrap/rdmolfiles.cpp`).

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::marvin::{
    MrvParserParams, mol_from_mrv_block, mol_from_mrv_file,
};

use crate::errors::to_py;
use crate::mol::PyMol;

/// A missing file is an `IOError`; other `FileParseException`s are logged as
/// warnings and every other error gives `None`.
fn finish(res: rdoxide::Result<Option<rdoxide::graphmol::RWMol>>) -> PyResult<Option<PyMol>> {
    match res {
        Ok(mol) => Ok(mol.map(|mol| PyMol { mol })),
        Err(e @ RDKitError::BadFile(_)) => Err(to_py(e)),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Ok(None)
        }
        Err(_) => Ok(None),
    }
}

/// Construct a molecule from a Marvin (Mrv) file.
///
///   ARGUMENTS:
///
///     - fileName: name of the file to read
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to true.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromMrvFile", signature = (molFileName, sanitize = true, removeHs = true))]
#[allow(non_snake_case)]
fn py_mol_from_mrv_file(
    py: Python<'_>,
    molFileName: &str,
    sanitize: bool,
    removeHs: bool,
) -> PyResult<Option<PyMol>> {
    let params = MrvParserParams {
        sanitize,
        remove_hs: removeHs,
    };
    finish(py.detach(|| mol_from_mrv_file(molFileName, &params)))
}

/// Construct a molecule from a Marvin (mrv) block.
///
///   ARGUMENTS:
///
///     - molBlock: string containing the Marvin block
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromMrvBlock", signature = (mrvBlock, sanitize = true, removeHs = true))]
#[allow(non_snake_case)]
fn py_mol_from_mrv_block(
    py: Python<'_>,
    mrvBlock: &Bound<'_, PyAny>,
    sanitize: bool,
    removeHs: bool,
) -> PyResult<Option<PyMol>> {
    // pyObjectToString accepts str or bytes
    let text = match mrvBlock.extract::<String>() {
        Ok(s) => s,
        Err(_) => String::from_utf8_lossy(&mrvBlock.extract::<Vec<u8>>()?).into_owned(),
    };
    let params = MrvParserParams {
        sanitize,
        remove_hs: removeHs,
    };
    finish(py.detach(|| mol_from_mrv_block(&text, &params).map(Some)))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_from_mrv_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_from_mrv_block, m)?)?;
    Ok(())
}
