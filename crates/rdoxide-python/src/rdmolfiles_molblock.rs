//! `MolFromMolBlock` / `MolFromMolFile` (upstream `Wrap/rdmolfiles.cpp`).

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::{
    MolFileParserParams, mol_from_mol_block_bytes, mol_from_mol_file,
};

use crate::errors::to_py;
use crate::mol::PyMol;

/// `pyObjectToString`: accepts str or bytes.
fn py_object_to_bytes(obj: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(s) = obj.extract::<String>() {
        return Ok(s.into_bytes());
    }
    obj.extract()
}

/// Upstream logs `FileParseException`s as warnings and swallows everything
/// else; only `BadFileException` reaches Python (as `OSError`).
fn finish(res: rdoxide::Result<Option<rdoxide::graphmol::RWMol>>) -> PyResult<Option<PyMol>> {
    match res {
        Ok(Some(mol)) => Ok(Some(PyMol { mol })),
        Ok(None) => Ok(None),
        Err(e @ RDKitError::BadFile(_)) => Err(to_py(e)),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Ok(None)
        }
        Err(_) => Ok(None),
    }
}

fn params(sanitize: bool, remove_hs: bool, strict_parsing: bool) -> MolFileParserParams {
    MolFileParserParams {
        sanitize,
        remove_hs,
        strict_parsing,
        ..MolFileParserParams::default()
    }
}

/// Construct a molecule from a Mol block.
///
///   ARGUMENTS:
///
///     - molBlock: string containing the Mol block
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///     - strictParsing: (optional) if this is false, the parser is more lax about.
///       correctness of the content.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromMolBlock", signature = (molBlock, sanitize = true, removeHs = true, strictParsing = true))]
#[allow(non_snake_case)]
fn mol_from_mol_block(
    py: Python<'_>,
    molBlock: &Bound<'_, PyAny>,
    sanitize: bool,
    removeHs: bool,
    strictParsing: bool,
) -> PyResult<Option<PyMol>> {
    let data = py_object_to_bytes(molBlock)?;
    let p = params(sanitize, removeHs, strictParsing);
    finish(py.detach(|| mol_from_mol_block_bytes(&data, &p)))
}

/// Construct a molecule from a Mol file.
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
///     - strictParsing: (optional) if this is false, the parser is more lax about.
///       correctness of the content.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromMolFile", signature = (molFileName, sanitize = true, removeHs = true, strictParsing = true))]
#[allow(non_snake_case)]
fn mol_from_mol_file_py(
    py: Python<'_>,
    molFileName: String,
    sanitize: bool,
    removeHs: bool,
    strictParsing: bool,
) -> PyResult<Option<PyMol>> {
    let p = params(sanitize, removeHs, strictParsing);
    finish(py.detach(|| mol_from_mol_file(&molFileName, &p)))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(mol_from_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(mol_from_mol_file_py, m)?)?;
    Ok(())
}
