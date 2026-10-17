//! `MolFromTPLFile`, `MolFromTPLBlock`, `MolToTPLFile` and `MolToTPLBlock`
//! (upstream `Wrap/rdmolfiles.cpp`).

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::tpl::{
    TplParserParams, mol_from_tpl_block, mol_from_tpl_file, mol_to_tpl_file, mol_to_tpl_text,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Construct a molecule from a TPL file.
///
///   ARGUMENTS:
///
///     - fileName: name of the file to read
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - skipFirstConf: (optional) skips reading the first conformer.
///       Defaults to False.
///       This should be set to True when reading TPLs written by
///       the CombiCode.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromTPLFile", signature = (fileName, sanitize = true, skipFirstConf = false))]
#[allow(non_snake_case)]
fn py_mol_from_tpl_file(
    py: Python<'_>,
    #[pyo3(from_py_with = crate::pyconv::path_like)] fileName: String,
    #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] skipFirstConf: bool,
) -> PyResult<Option<PyMol>> {
    let params = TplParserParams {
        sanitize,
        skip_first_conf: skipFirstConf,
    };
    match py.detach(|| mol_from_tpl_file(&fileName, &params)) {
        Ok(mol) => Ok(mol.map(|mol| PyMol { mol })),
        // a missing file is an IOError, any other failure gives None
        Err(e @ RDKitError::BadFile(_)) => Err(crate::errors::to_py(e)),
        Err(_) => Ok(None),
    }
}

/// Construct a molecule from a TPL block.
///
///   ARGUMENTS:
///
///     - tplBlock: string containing the TPL block
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - skipFirstConf: (optional) skips reading the first conformer.
///       Defaults to False.
///       This should be set to True when reading TPLs written by
///       the CombiCode.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromTPLBlock", signature = (tplBlock, sanitize = true, skipFirstConf = false))]
#[allow(non_snake_case)]
fn py_mol_from_tpl_block(
    py: Python<'_>,
    #[pyo3(from_py_with = crate::pyconv::text_string)] tplBlock: String,
    #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] skipFirstConf: bool,
) -> Option<PyMol> {
    let params = TplParserParams {
        sanitize,
        skip_first_conf: skipFirstConf,
    };
    py.detach(|| mol_from_tpl_block(&tplBlock, &params))
        .ok()
        .flatten()
        .map(|mol| PyMol { mol })
}

/// Writes a molecule to a TPL file.
///
///   ARGUMENTS:
///
///     - mol: the molecule
///     - fileName: name of the file to write
///     - partialChargeProp: name of the property to use for partial charges
///       Defaults to '_GasteigerCharge'.
///     - writeFirstConfTwice: Defaults to False.
///       This should be set to True when writing TPLs to be read by
///       the CombiCode.
#[pyfunction]
#[pyo3(name = "MolToTPLFile", signature = (mol, fileName, partialChargeProp = "_GasteigerCharge", writeFirstConfTwice = false))]
#[allow(non_snake_case)]
fn py_mol_to_tpl_file(
    mol: PyRef<'_, PyMol>,
    #[pyo3(from_py_with = crate::pyconv::path_like)] fileName: String,
    partialChargeProp: &str,
    #[pyo3(from_py_with = crate::pyconv::truthy)] writeFirstConfTwice: bool,
) -> PyResult<()> {
    mol_to_tpl_file(&mol.mol, &fileName, partialChargeProp, writeFirstConfTwice).py()
}

/// Returns the Tpl block for a molecule.
///
///   ARGUMENTS:
///
///     - mol: the molecule
///     - partialChargeProp: name of the property to use for partial charges
///       Defaults to '_GasteigerCharge'.
///     - writeFirstConfTwice: Defaults to False.
///       This should be set to True when writing TPLs to be read by
///       the CombiCode.
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToTPLBlock", signature = (mol, partialChargeProp = "_GasteigerCharge", writeFirstConfTwice = false))]
#[allow(non_snake_case)]
fn py_mol_to_tpl_block(
    mol: PyRef<'_, PyMol>,
    partialChargeProp: &str,
    #[pyo3(from_py_with = crate::pyconv::truthy)] writeFirstConfTwice: bool,
) -> String {
    mol_to_tpl_text(&mol.mol, partialChargeProp, writeFirstConfTwice)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_from_tpl_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_from_tpl_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_tpl_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_tpl_block, m)?)?;
    Ok(())
}
