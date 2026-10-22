//! `rdoxide.Chem.rdinchi`: InChI support (upstream
//! `External/INCHI-API/Wrap/pyInchi.cpp`).
//!
//! The full upstream surface works (the `/z` polymer layer excepted, see
//! docs/INCHI_PORT.md).

use pyo3::prelude::*;
use rdoxide::external::inchi_api::{
    ExtraInchiReturnValues, get_inchi_version, inchi_to_inchi_key, inchi_to_mol,
    mol_block_to_inchi, mol_to_inchi, mol_to_inchi_key,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

fn opt(options: &str) -> Option<&str> {
    (!options.is_empty()).then_some(options)
}

/// return a ROMol for a InChI string
///   Returns:
///     a tuple with:
///       - the molecule
///       - the return code from the InChI conversion
///       - a string with any messages from the InChI conversion
///       - a string with any log messages from the InChI conversion
#[pyfunction]
#[pyo3(name = "InchiToMol", signature = (inchi, sanitize = true, removeHs = true))]
#[allow(non_snake_case)]
fn inchi_to_mol_py(
    py: Python<'_>,
    inchi: &str,
    sanitize: bool,
    removeHs: bool,
) -> PyResult<(Option<PyMol>, i32, String, String)> {
    let mut rv = ExtraInchiReturnValues::default();
    let mol = py
        .detach(|| inchi_to_mol(inchi, &mut rv, sanitize, removeHs))
        .py()?;
    Ok((
        mol.map(|mol| PyMol { mol }),
        rv.return_code,
        rv.message,
        rv.log,
    ))
}

/// return the InChI for a ROMol molecule.
///
///   Arguments:
///     - mol: the molecule to use.
///     - options: the InChI generation options.
///       Options should be prefixed with either a - or a /
///       Available options are explained in the InChI technical FAQ:
///       https://www.inchi-trust.org/technical-faq/#15.14
///       and the User Guide:
///       https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
///   Returns:
///     a tuple with:
///       - the InChI
///       - the return code from the InChI conversion
///       - a string with any messages from the InChI conversion
///       - a string with any log messages from the InChI conversion
///       - a string with the InChI AuxInfo
#[pyfunction]
#[pyo3(name = "MolToInchi", signature = (mol, options = String::new()))]
fn mol_to_inchi_py(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    options: String,
) -> PyResult<(String, i32, String, String, String)> {
    let mut rv = ExtraInchiReturnValues::default();
    let mol = &mol.mol;
    let inchi = py
        .detach(|| mol_to_inchi(mol, &mut rv, opt(&options)))
        .py()?;
    Ok((inchi, rv.return_code, rv.message, rv.log, rv.aux_info))
}

/// return the InChI for a ROMol molecule.
///
///   Arguments:
///     - molblock: the mol block to use.
///     - options: the InChI generation options.
///       Options should be prefixed with either a - or a /
///       Available options are explained in the InChI technical FAQ:
///       https://www.inchi-trust.org/technical-faq/#15.14
///       and the User Guide:
///       https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
///   Returns:
///     a tuple with:
///       - the InChI
///       - the return code from the InChI conversion
///       - a string with any messages from the InChI conversion
///       - a string with any log messages from the InChI conversion
///       - a string with the InChI AuxInfo
#[pyfunction]
#[pyo3(name = "MolBlockToInchi", signature = (molblock, options = String::new()))]
fn mol_block_to_inchi_py(
    py: Python<'_>,
    molblock: &str,
    options: String,
) -> PyResult<(String, i32, String, String, String)> {
    let mut rv = ExtraInchiReturnValues::default();
    let inchi = py
        .detach(|| mol_block_to_inchi(molblock, &mut rv, opt(&options)))
        .py()?;
    Ok((inchi, rv.return_code, rv.message, rv.log, rv.aux_info))
}

/// return the InChI key for an InChI string
#[pyfunction]
#[pyo3(name = "InchiToInchiKey", signature = (inchi))]
fn inchi_to_inchi_key_py(inchi: &str) -> String {
    inchi_to_inchi_key(inchi)
}

/// return the InChI key for a ROMol molecule.
///
///   Arguments:
///     - mol: the molecule to use.
///     - options: the InChI generation options.
///       Options should be prefixed with either a - or a /
///       Available options are explained in the InChI technical FAQ:
///       https://www.inchi-trust.org/technical-faq/#15.14
///       and the User Guide available from:
///       https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
///   Returns: the InChI key
#[pyfunction]
#[pyo3(name = "MolToInchiKey", signature = (mol, options = String::new()))]
fn mol_to_inchi_key_py(py: Python<'_>, mol: PyRef<'_, PyMol>, options: String) -> PyResult<String> {
    let mol = &mol.mol;
    py.detach(|| mol_to_inchi_key(mol, opt(&options))).py()
}

/// returns the version of the InChI software being used
#[pyfunction]
#[pyo3(name = "GetInchiVersion")]
fn get_inchi_version_py() -> &'static str {
    get_inchi_version()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(inchi_to_mol_py, m)?)?;
    m.add_function(wrap_pyfunction!(mol_to_inchi_py, m)?)?;
    m.add_function(wrap_pyfunction!(mol_block_to_inchi_py, m)?)?;
    m.add_function(wrap_pyfunction!(inchi_to_inchi_key_py, m)?)?;
    m.add_function(wrap_pyfunction!(mol_to_inchi_key_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_inchi_version_py, m)?)?;
    Ok(())
}
