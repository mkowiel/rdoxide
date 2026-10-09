//! `CalcAUTOCORR3D`, `CalcRDF` and `CalcMORSE` of `rdkit.Chem.rdMolDescriptors`
//! (upstream `Descriptors/AUTOCORR3D.cpp`, `RDF.cpp`, `MORSE.cpp` and their
//! wrappers in `Wrap/rdMolDescriptors.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::descriptors::autocorr3d::{AUTOCORR3D_VERSION, autocorr3d};
use rdoxide::graphmol::descriptors::morse::{MORSE_VERSION, morse};
use rdoxide::graphmol::descriptors::rdf::{RDF_VERSION, rdf};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors_3d::require_conformer;

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

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("_CalcAUTOCORR3D_version", AUTOCORR3D_VERSION)?;
    m.add("_CalcRDF_version", RDF_VERSION)?;
    m.add("_CalcMORSE_version", MORSE_VERSION)?;
    m.add_function(wrap_pyfunction!(calc_autocorr3d, m)?)?;
    m.add_function(wrap_pyfunction!(calc_rdf, m)?)?;
    m.add_function(wrap_pyfunction!(calc_morse, m)?)
}
