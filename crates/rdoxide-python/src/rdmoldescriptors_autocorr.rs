//! `rdkit.Chem.rdMolDescriptors.MQNs_` and `CalcAUTOCORR2D` (upstream
//! `Descriptors/MQN.cpp`, `AUTOCORR2D.cpp` and their wrappers in
//! `Wrap/rdMolDescriptors.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::descriptors::autocorr2d::{AUTOCORR2D_VERSION, autocorr2d};
use rdoxide::graphmol::descriptors::mqn::calc_mqns;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// returns the MQN descriptors for a molecule
#[pyfunction]
#[pyo3(name = "MQNs_", signature = (mol, force = None))]
fn py_mqns(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    force: Option<&Bound<'_, PyAny>>,
) -> PyResult<Vec<u32>> {
    // `force` is accepted for API compatibility; upstream ignores it
    bool_arg(force, false)?;
    let m = &mol.borrow().mol;
    py.detach(|| calc_mqns(m)).py()
}

/// Returns 2D Autocorrelation descriptors vector
#[pyfunction]
#[pyo3(name = "CalcAUTOCORR2D", signature = (mol, CustomAtomProperty = ""))]
#[allow(non_snake_case)]
fn py_autocorr2d(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    CustomAtomProperty: &str,
) -> PyResult<Vec<f64>> {
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    py.detach(|| autocorr2d(core, CustomAtomProperty)).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("_CalcAUTOCORR2D_version", AUTOCORR2D_VERSION)?;
    m.add_function(wrap_pyfunction!(py_mqns, m)?)?;
    m.add_function(wrap_pyfunction!(py_autocorr2d, m)?)
}
