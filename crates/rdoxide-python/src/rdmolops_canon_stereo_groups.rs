//! `StereoGroupAbsOptions` and `CanonicalizeStereoGroups` of
//! `rdoxide.Chem.rdmolops` (upstream `GraphMol/Wrap/MolOps.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::graphmol::smiles_parse::canonicalize_stereo_groups::{
    StereoGroupAbsOptions, canonicalize_stereo_groups,
};

use crate::enums::add_enum;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

const ABS_OPTIONS: &str = "StereoGroupAbsOptions";

/// Rationalize Enhanced Stereo indications to a canonical form
///
///   ARGUMENTS:
///
///     - molecule: the molecule to update
///     - StereoGroupAbsOptions outputAbsoluteGroups: controls output of abs groups:
///       one of: OnlyIncludeWhenOtherGroupsExist, NeverInclude, AlwaysInclude
///     - maxStereoGroups: maximum number of OR or AND stereo groups to process
///       (default is 12)
#[pyfunction]
#[pyo3(name = "CanonicalizeStereoGroups", signature = (mol, outputAbsoluteGroups = None, maxStereoGroups = 12))]
#[allow(non_snake_case)]
fn py_canonicalize_stereo_groups(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    outputAbsoluteGroups: Option<i64>,
    maxStereoGroups: u32,
) -> PyResult<PyMol> {
    let options = match outputAbsoluteGroups.unwrap_or(0) {
        0 => StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist,
        1 => StereoGroupAbsOptions::NeverInclude,
        2 => StereoGroupAbsOptions::AlwaysInclude,
        _ => return Err(PyValueError::new_err("invalid StereoGroupAbsOptions value")),
    };
    let m = &mol.mol;
    let res = py
        .detach(|| canonicalize_stereo_groups(m, options, maxStereoGroups))
        .py()?;
    Ok(PyMol { mol: res })
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum(
        m,
        ABS_OPTIONS,
        &[
            ("OnlyIncludeWhenOtherGroupsExist", 0),
            ("NeverInclude", 1),
            ("AlwaysInclude", 2),
        ],
        None,
    )?;
    m.add_function(wrap_pyfunction!(py_canonicalize_stereo_groups, m)?)
}
