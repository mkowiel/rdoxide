//! Subgraph / path enumeration functions of `rdoxide.Chem.rdmolops`
//! (upstream `GraphMol/Wrap/MolOps.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::subgraphs;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::vects::{PyListVectI, PyVectI};

/// Finds all subgraphs of a particular length in a molecule
#[pyfunction]
#[pyo3(name = "FindAllSubgraphsOfLengthN", signature = (mol, length, useHs = false, rootedAtAtom = -1))]
#[allow(non_snake_case)]
fn find_all_subgraphs_of_length_n(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    length: u32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useHs: bool,
    rootedAtAtom: i32,
) -> PyResult<Py<PyAny>> {
    let r = subgraphs::find_all_subgraphs_of_length_n(&mol.mol, length, useHs, rootedAtAtom, None)
        .py()?;
    PyListVectI::new_py(py, r)
}

/// Finds all subgraphs of a particular length in a molecule
#[pyfunction]
#[pyo3(name = "FindAllSubgraphsOfLengthMToN", signature = (mol, min, max, useHs = false, rootedAtAtom = -1))]
#[allow(non_snake_case)]
fn find_all_subgraphs_of_length_m_to_n<'py>(
    py: Python<'py>,
    mol: PyRef<'_, PyMol>,
    min: u32,
    max: u32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useHs: bool,
    rootedAtAtom: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    if min > max {
        return Err(PyValueError::new_err("lowerLen > upperLen"));
    }
    let m = subgraphs::find_all_subgraphs_of_lengths_m_to_n(
        &mol.mol,
        min,
        max,
        useHs,
        rootedAtAtom,
        None,
    )
    .py()?;
    let mut out = Vec::new();
    for l in min..=max {
        let paths = m.get(&(l as i32)).cloned().unwrap_or_default();
        let items = paths
            .into_iter()
            .map(|p| PyTuple::new(py, p))
            .collect::<PyResult<Vec<_>>>()?;
        out.push(PyList::new(py, items)?);
    }
    PyTuple::new(py, out)
}

/// Finds unique subgraphs of a particular length in a molecule
#[pyfunction]
#[pyo3(name = "FindUniqueSubgraphsOfLengthN", signature = (mol, length, useHs = false, useBO = true, rootedAtAtom = -1))]
#[allow(non_snake_case)]
fn find_unique_subgraphs_of_length_n(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    length: u32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useHs: bool,
    useBO: bool,
    rootedAtAtom: i32,
) -> PyResult<Py<PyAny>> {
    let r = subgraphs::find_unique_subgraphs_of_length_n(
        &mol.mol,
        length,
        useHs,
        useBO,
        rootedAtAtom,
        None,
    )
    .py()?;
    PyListVectI::new_py(py, r)
}

/// Finds all paths of a particular length in a molecule
#[pyfunction]
#[pyo3(name = "FindAllPathsOfLengthN", signature = (mol, length, useBonds = true, useHs = false, rootedAtAtom = -1, onlyShortestPaths = false))]
#[allow(non_snake_case)]
fn find_all_paths_of_length_n(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    length: u32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useBonds: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useHs: bool,
    rootedAtAtom: i32,
    onlyShortestPaths: bool,
) -> PyResult<Py<PyAny>> {
    let r = subgraphs::find_all_paths_of_length_n(
        &mol.mol,
        length,
        useBonds,
        useHs,
        rootedAtAtom,
        onlyShortestPaths,
        None,
    )
    .py()?;
    PyListVectI::new_py(py, r)
}

/// Finds the bonds within a certain radius of an atom in a molecule
#[pyfunction]
#[pyo3(name = "FindAtomEnvironmentOfRadiusN", signature = (mol, radius, rootedAtAtom, useHs = false, enforceSize = true, atomMap = None))]
#[allow(non_snake_case)]
fn find_atom_environment_of_radius_n(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    radius: u32,
    rootedAtAtom: u32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useHs: bool,
    enforceSize: bool,
    atomMap: Option<&Bound<'_, PyAny>>,
) -> PyResult<Py<PyAny>> {
    let (path, amap) = subgraphs::find_atom_environment_of_radius_n(
        &mol.mol,
        radius,
        rootedAtAtom,
        useHs,
        enforceSize,
    )
    .py()?;
    if let Some(d) = atomMap
        && !d.is_none()
    {
        let d = d.cast::<PyDict>()?;
        d.clear();
        for (k, v) in amap {
            d.set_item(k, v)?;
        }
    }
    PyVectI::new_py(py, path)
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_all_subgraphs_of_length_n, m)?)?;
    m.add_function(wrap_pyfunction!(find_all_subgraphs_of_length_m_to_n, m)?)?;
    m.add_function(wrap_pyfunction!(find_unique_subgraphs_of_length_n, m)?)?;
    m.add_function(wrap_pyfunction!(find_all_paths_of_length_n, m)?)?;
    m.add_function(wrap_pyfunction!(find_atom_environment_of_radius_n, m)?)?;
    Ok(())
}
