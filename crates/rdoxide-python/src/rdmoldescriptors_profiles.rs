//! `CalcAUTOCORR3D`, `CalcRDF`, `CalcMORSE`, `CalcCoulombMat` and `CalcEEMcharges` and `CalcWHIM` of `rdkit.Chem.rdMolDescriptors`
//! (upstream `Descriptors/AUTOCORR3D.cpp`, `RDF.cpp`, `MORSE.cpp`,
//! `CoulombMat.cpp`, `EEM.cpp`, `WHIM.cpp` and their
//! wrappers in `Wrap/rdMolDescriptors.cpp`).

use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::descriptors::autocorr3d::{AUTOCORR3D_VERSION, autocorr3d};
use rdoxide::graphmol::descriptors::coulomb_mat::{COULOMB_MAT_VERSION, coulomb_mat};
use rdoxide::graphmol::descriptors::eem::{EEM_VERSION, eem};
use rdoxide::graphmol::descriptors::morse::{MORSE_VERSION, morse};
use rdoxide::graphmol::descriptors::rdf::{RDF_VERSION, rdf};
use rdoxide::graphmol::descriptors::whim::{WHIM_VERSION, whim};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors_3d::require_conformer;
use crate::vects::PyVectD;

/// Returns 3D Autocorrelation descriptors vector
#[pyfunction]
#[pyo3(name = "CalcAUTOCORR3D", signature = (mol, confId = -1, CustomAtomProperty = ""))]
#[allow(non_snake_case)]
fn calc_autocorr3d(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    confId: i32,
    CustomAtomProperty: &str,
) -> PyResult<Vec<f64>> {
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    require_conformer(core, "Code/GraphMol/Descriptors/AUTOCORR3D.cpp", 233)?;
    py.detach(|| autocorr3d(core, confId, CustomAtomProperty))
        .py()
}

/// Returns radial distribution fonction descriptors (RDF)
#[pyfunction]
#[pyo3(name = "CalcRDF", signature = (mol, confId = -1, CustomAtomProperty = ""))]
#[allow(non_snake_case)]
fn calc_rdf(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    confId: i32,
    CustomAtomProperty: &str,
) -> PyResult<Vec<f64>> {
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    require_conformer(core, "Code/GraphMol/Descriptors/RDF.cpp", 192)?;
    py.detach(|| rdf(core, confId, CustomAtomProperty)).py()
}

/// Returns Molecule Representation of Structures based on Electron diffraction descriptors
#[pyfunction]
#[pyo3(name = "CalcMORSE", signature = (mol, confId = -1, CustomAtomProperty = ""))]
#[allow(non_snake_case)]
fn calc_morse(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    confId: i32,
    CustomAtomProperty: &str,
) -> PyResult<Vec<f64>> {
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    require_conformer(core, "Code/GraphMol/Descriptors/MORSE.cpp", 172)?;
    py.detach(|| morse(core, confId, CustomAtomProperty)).py()
}

/// Returns severals Coulomb randomized matrices
#[pyfunction]
#[pyo3(name = "CalcCoulombMat", signature = (mol, confId = -1))]
#[allow(non_snake_case)]
fn calc_coulomb_mat<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    confId: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    let m = mol.borrow();
    let core = &m.mol;
    require_conformer(core, "Code/GraphMol/Descriptors/CoulombMat.cpp", 40)?;
    let rows = py.detach(|| coulomb_mat(core, confId)).py()?;
    let rows = rows
        .into_iter()
        .map(|r| PyVectD::new_py(py, r))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, rows)
}

/// Returns EEM atomic partial charges
#[pyfunction]
#[pyo3(name = "CalcEEMcharges", signature = (mol, confId = -1))]
#[allow(non_snake_case)]
fn calc_eem_charges(py: Python<'_>, mol: &Bound<'_, PyMol>, confId: i32) -> PyResult<Vec<f64>> {
    let m = mol.borrow();
    let core = &m.mol;
    require_conformer(core, "Code/GraphMol/Descriptors/EEM.cpp", 235)?;
    py.detach(|| eem(core, confId)).py()
}

/// Returns the WHIM descriptors vector
#[pyfunction]
#[pyo3(name = "CalcWHIM", signature = (mol, confId = -1, thresh = 0.001, CustomAtomProperty = ""))]
#[allow(non_snake_case)]
fn calc_whim(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    confId: i32,
    thresh: f64,
    CustomAtomProperty: &str,
) -> PyResult<Vec<f64>> {
    let m = mol.borrow();
    let core = &m.mol;
    require_conformer(core, "Code/GraphMol/Descriptors/WHIM.cpp", 382)?;
    py.detach(|| whim(core, confId, thresh, CustomAtomProperty))
        .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("_CalcAUTOCORR3D_version", AUTOCORR3D_VERSION)?;
    m.add("_CalcRDF_version", RDF_VERSION)?;
    m.add("_CalcMORSE_version", MORSE_VERSION)?;
    m.add("_CalcCoulombMat_version", COULOMB_MAT_VERSION)?;
    m.add_function(wrap_pyfunction!(calc_coulomb_mat, m)?)?;
    m.add("_CalcEMMcharges_version", EEM_VERSION)?;
    m.add_function(wrap_pyfunction!(calc_eem_charges, m)?)?;
    m.add("_CalcWHIM_version", WHIM_VERSION)?;
    m.add_function(wrap_pyfunction!(calc_whim, m)?)?;
    m.add_function(wrap_pyfunction!(calc_autocorr3d, m)?)?;
    m.add_function(wrap_pyfunction!(calc_rdf, m)?)?;
    m.add_function(wrap_pyfunction!(calc_morse, m)?)
}
