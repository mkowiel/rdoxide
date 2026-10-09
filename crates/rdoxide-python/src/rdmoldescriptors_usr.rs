//! The USR/USRCAT functions of `rdkit.Chem.rdMolDescriptors` (upstream
//! `Descriptors/USRDescriptor.cpp` and its wrappers in
//! `Wrap/rdMolDescriptors.cpp`).

use pyo3::exceptions::{PyIndexError, PyOverflowError, PyValueError};
use pyo3::prelude::*;
use rdoxide::geometry::Point3D;
use rdoxide::graphmol::descriptors::usr::{
    USR_BLOCK, calc_usr_distributions, calc_usr_distributions_from_points,
    calc_usr_from_distributions, calc_usr_score, usr, usrcat,
};

use crate::errors::{Check, IntoPyResult, Site, check_failed, expect_list, no_converter};
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;

const USR_FILE: &str = "Code/GraphMol/Descriptors/USRDescriptor.cpp";

fn require_usr_input(mol: &rdoxide::graphmol::romol::ROMol) -> PyResult<()> {
    if mol.num_conformers() == 0 {
        return Err(PyValueError::new_err("no conformers"));
    }
    if mol.num_atoms() < 3 {
        return Err(PyValueError::new_err("too few atoms (minimum three)"));
    }
    Ok(())
}

/// `python::extract<unsigned int>` of every element of `group`.
fn unsigned_atoms(group: &Bound<'_, PyAny>) -> PyResult<Vec<u32>> {
    group
        .extract::<Vec<Bound<'_, PyAny>>>()?
        .iter()
        .map(|n| match n.extract::<i64>() {
            Ok(v) if v < 0 => Err(PyOverflowError::new_err(
                "can't convert negative value to unsigned int",
            )),
            _ => n.extract::<u32>(),
        })
        .collect()
}

/// The `Point3D` instances of a Python sequence.
fn points_of(seq: &Bound<'_, PyAny>) -> PyResult<Vec<Point3D>> {
    seq.extract::<Vec<Bound<'_, PyAny>>>()?
        .iter()
        .map(|item| match item.extract::<PyRef<'_, PyPoint3D>>() {
            Ok(p) => Ok(p.p),
            Err(_) => Err(no_converter("RDGeom::Point3D", item)),
        })
        .collect()
}

fn new_points(py: Python<'_>, pts: &[Point3D]) -> PyResult<Vec<Py<PyPoint3D>>> {
    pts.iter().map(|&p| Py::new(py, PyPoint3D { p })).collect()
}

/// Returns a USR descriptor for one conformer of a molecule
#[pyfunction]
#[pyo3(name = "GetUSR", signature = (mol, confId = -1))]
#[allow(non_snake_case)]
fn get_usr(py: Python<'_>, mol: &Bound<'_, PyMol>, confId: i32) -> PyResult<Vec<f64>> {
    let m = mol.borrow();
    let core = &m.mol;
    require_usr_input(core)?;
    py.detach(|| usr(core, confId)).py()
}

/// Returns the four USR distance distributions for a set of coordinates
#[pyfunction]
#[pyo3(name = "GetUSRDistributions", signature = (coords, points = None))]
fn get_usr_distributions(
    py: Python<'_>,
    coords: &Bound<'_, PyAny>,
    points: Option<&Bound<'_, PyAny>>,
) -> PyResult<Vec<Vec<f64>>> {
    let coords = points_of(coords)?;
    if coords.is_empty() {
        return Err(PyValueError::new_err("no coordinates"));
    }
    let (dist, pts) = calc_usr_distributions(&coords).py()?;
    if let Some(out) = points.filter(|p| !p.is_none()) {
        let list = expect_list(out)?;
        for p in new_points(py, &pts)? {
            list.append(p)?;
        }
    }
    Ok(dist.into())
}

/// Returns the USR distance distributions for a set of coordinates and points
#[pyfunction]
#[pyo3(name = "GetUSRDistributionsFromPoints", signature = (coords, points))]
fn get_usr_distributions_from_points(
    coords: &Bound<'_, PyAny>,
    points: &Bound<'_, PyAny>,
) -> PyResult<Vec<Vec<f64>>> {
    let coords = points_of(coords)?;
    if coords.is_empty() {
        return Err(PyValueError::new_err("no coordinates"));
    }
    let points = points_of(points)?;
    if points.is_empty() {
        return Err(PyValueError::new_err("no points"));
    }
    Ok(calc_usr_distributions_from_points(&coords, &points))
}

/// Returns the USR descriptor from a set of distance distributions
#[pyfunction]
#[pyo3(name = "GetUSRFromDistributions", signature = (distances))]
fn get_usr_from_distributions(distances: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    let rows = distances.extract::<Vec<Bound<'_, PyAny>>>()?;
    if rows.is_empty() {
        return Err(PyValueError::new_err("no distances"));
    }
    let dist = rows
        .iter()
        .map(|row| {
            let d = row.extract::<Vec<f64>>()?;
            if d.is_empty() {
                return Err(PyValueError::new_err("distances missing"));
            }
            Ok(d)
        })
        .collect::<PyResult<Vec<_>>>()?;
    if dist.len() != 4 {
        // the 12-value output buffer upstream only fits four distributions
        return Err(check_failed(
            Site {
                check: Check::Precondition,
                file: USR_FILE,
                line: 300,
                expr: "descriptor.size() == 3 * dist.size()",
            },
            "descriptor must have 3 times more elements than dist",
        ));
    }
    Ok(calc_usr_from_distributions(&dist))
}

/// Returns the USR score for two USR or USRCAT descriptors
#[pyfunction]
#[pyo3(name = "GetUSRScore", signature = (descriptor1, descriptor2, weights = None))]
fn get_usr_score(
    descriptor1: &Bound<'_, PyAny>,
    descriptor2: &Bound<'_, PyAny>,
    weights: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let d1 = descriptor1.extract::<Vec<f64>>()?;
    let d2 = descriptor2.extract::<Vec<f64>>()?;
    if d1.len() != d2.len() {
        return Err(PyValueError::new_err(
            "descriptors must have the same length",
        ));
    }
    let num_weights = d1.len() / USR_BLOCK;
    let given = match weights {
        Some(w) => w.extract::<Vec<f64>>()?,
        None => Vec::new(),
    };
    let w = match given.len() {
        0 => vec![1.0; num_weights],
        n if n == num_weights => given,
        _ => return Err(PyValueError::new_err("number of weights is not correct")),
    };
    calc_usr_score(&d1, &d2, &w).py()
}

/// Returns a USRCAT descriptor for one conformer of a molecule
#[pyfunction]
#[pyo3(name = "GetUSRCAT", signature = (mol, atomSelections = None, confId = -1))]
#[allow(non_snake_case)]
fn get_usrcat(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    atomSelections: Option<&Bound<'_, PyAny>>,
    confId: i32,
) -> PyResult<Vec<f64>> {
    let m = mol.borrow();
    let core = &m.mol;
    require_usr_input(core)?;
    // selections hold 1-based atom numbers
    let selections = match atomSelections.filter(|s| !s.is_none()) {
        None => None,
        Some(sel) => {
            let list = expect_list(sel)?;
            if list.is_empty() {
                return Err(PyValueError::new_err("empty atom selections"));
            }
            let ids = list
                .iter()
                .map(|group| {
                    unsigned_atoms(&group)?
                        .into_iter()
                        .map(|n| {
                            n.checked_sub(1)
                                .ok_or_else(|| PyIndexError::new_err("atom selections are 1-based"))
                        })
                        .collect::<PyResult<Vec<u32>>>()
                })
                .collect::<PyResult<Vec<_>>>()?;
            Some(ids)
        }
    };
    py.detach(|| usrcat(core, selections.as_deref(), confId))
        .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_usr, m)?)?;
    m.add_function(wrap_pyfunction!(get_usr_distributions, m)?)?;
    m.add_function(wrap_pyfunction!(get_usr_distributions_from_points, m)?)?;
    m.add_function(wrap_pyfunction!(get_usr_from_distributions, m)?)?;
    m.add_function(wrap_pyfunction!(get_usr_score, m)?)?;
    m.add_function(wrap_pyfunction!(get_usrcat, m)?)
}
