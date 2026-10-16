//! `rdoxide.Chem.rdPartialCharges` (upstream `PartialCharges/Wrap/rdPartialCharges.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::partial_charges::gasteiger::{self, DEFAULT_ITERATIONS};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// Compute Gasteiger partial charges for molecule
///
///  The charges are computed using an iterative procedure presented in
///
///  Ref : J.Gasteiger, M. Marseli, Iterative Equalization of Oribital Electronegatiity
///  A Rapid Access to Atomic Charges, Tetrahedron Vol 36 p3219 1980
///
///  The computed charges are stored on each atom are stored a computed property ( under the name
///  _GasteigerCharge). In addition, each atom also stored the total charge for the implicit hydrogens
///  on the atom (under the property name _GasteigerHCharge)
///
///  ARGUMENTS:
///
///     - mol : the molecule of interrest
///     - nIter : number of iteration (defaults to 12)
///     - throwOnParamFailure : toggles whether or not an exception should be raised if parameters
///       for an atom cannot be found.  If this is false (the default), all parameters for unknown
///       atoms will be set to zero.  This has the effect of removing that atom from the iteration.
///
#[pyfunction]
#[pyo3(
    name = "ComputeGasteigerCharges",
    signature = (mol, nIter = DEFAULT_ITERATIONS, throwOnParamFailure = None)
)]
#[allow(non_snake_case)]
fn compute_gasteiger_charges(
    mol: &Bound<'_, PyMol>,
    nIter: i32,
    throwOnParamFailure: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let throw = bool_arg(throwOnParamFailure, false)?;
    gasteiger::compute_gasteiger_charges(&mut mol.borrow_mut().mol, nIter, throw).py()?;
    Ok(())
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to set partial charges - currently Gasteiger Charges",
    )?;
    m.add_function(wrap_pyfunction!(compute_gasteiger_charges, m)?)
}
