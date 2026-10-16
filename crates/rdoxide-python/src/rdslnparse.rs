//! `rdoxide.Chem.rdSLNParse`: Sybyl line notation (upstream
//! `GraphMol/SLNParse/Wrap/rdSLNParse.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::sln_parse::{sln_query_to_mol, sln_to_mol};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Construct a molecule from an SLN string.
///
///     ARGUMENTS:
///
///     - SLN: the SLN string
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///   RETURNS:
///
///     a Mol object, None on failure.
///
///   NOTE: the SLN should not contain query information or properties. To build a
///     query from SLN, use MolFromQuerySLN.
#[pyfunction]
#[pyo3(name = "MolFromSLN", signature = (SLN, sanitize = true, debugParser = false))]
#[allow(non_snake_case)]
fn mol_from_sln(
    py: Python<'_>,
    SLN: &str,
    sanitize: bool,
    debugParser: bool,
) -> PyResult<Option<PyMol>> {
    let res = py
        .detach(|| sln_to_mol(SLN, sanitize, i32::from(debugParser)))
        .py()?;
    Ok(res.map(|mol| PyMol { mol }))
}

/// Construct a query molecule from an SLN string.
///
///   ARGUMENTS:
///
///     - SLN: the SLN string
///
///     - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached
///       heavy atoms. Defaults to False.
///
///   RETURNS:
///
///     a Mol object suitable for using in substructure queries, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromQuerySLN", signature = (SLN, mergeHs = true, debugParser = false))]
#[allow(non_snake_case)]
fn mol_from_query_sln(
    py: Python<'_>,
    SLN: &str,
    mergeHs: bool,
    debugParser: bool,
) -> PyResult<Option<PyMol>> {
    let res = py
        .detach(|| sln_query_to_mol(SLN, mergeHs, i32::from(debugParser)))
        .py()?;
    Ok(res.map(|mol| PyMol { mol }))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing classes and functions for working with Sybyl line notation (SLN).",
    )?;
    m.add_function(wrap_pyfunction!(mol_from_sln, m)?)?;
    m.add_function(wrap_pyfunction!(mol_from_query_sln, m)?)?;
    Ok(())
}
