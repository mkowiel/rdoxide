//! Matrix functions of `rdoxide.Chem.rdmolops` (upstream `GraphMol/Wrap/MolOps.cpp`):
//! `GetDistanceMatrix`, `Get3DDistanceMatrix` and `GetAdjacencyMatrix`.
//!
//! Like upstream, the matrices are cached on the molecule as computed
//! properties (`force` recomputes, `prefix` namespaces the cache key).
//! The `_Paths` companion matrices upstream also stores are not kept.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::graphmol::matrices;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// A row-major `n x n` buffer as a NumPy array of the given dtype.
fn square_array<'py>(
    py: Python<'py>,
    bytes: Vec<u8>,
    n: usize,
    dtype: &str,
) -> PyResult<Bound<'py, PyAny>> {
    // raw bytes in a bytearray give a writable array without Python lists
    let buffer = pyo3::types::PyByteArray::new(py, &bytes);
    let np = py.import("numpy")?;
    np.call_method1("frombuffer", (buffer, dtype))?
        .call_method1("reshape", ((n, n),))
}

fn f64_bytes(values: &[f64]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_ne_bytes()).collect()
}

/// Returns the molecule's topological distance matrix.
///
/// ARGUMENTS:
///
///   - mol: the molecule to use
///
///   - useBO: (optional) toggles use of bond orders in calculating the distance matrix
///     Default value is 0.
///
///   - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
///     matrix (to return a "Barysz" distance matrix).
///     Default value is 0.
///
///   - force: (optional) forces the calculation to proceed, even if there is a cached value.
///     Default value is 0.
///
///   - prefix: (optional, internal use) sets the prefix used in the property cache
///     Default value is .
///
///  RETURNS:
///
///    a Numeric array of floats with the distance matrix
#[pyfunction]
#[pyo3(name = "GetDistanceMatrix", signature = (mol, useBO = None, useAtomWts = None, force = None, prefix = ""))]
#[allow(non_snake_case)]
fn get_distance_matrix<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    useBO: Option<&Bound<'py, PyAny>>,
    useAtomWts: Option<&Bound<'py, PyAny>>,
    force: Option<&Bound<'py, PyAny>>,
    prefix: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let use_bo = bool_arg(useBO, false)?;
    let use_atom_wts = bool_arg(useAtomWts, false)?;
    let force = bool_arg(force, false)?;
    let (n, d) = {
        let mut m = mol.borrow_mut();
        let d = matrices::get_distance_mat_cached(&mut m.mol, use_bo, use_atom_wts, force, prefix)
            .py()?;
        (m.mol.num_atoms() as usize, d)
    };
    square_array(py, f64_bytes(&d), n, "float64")
}

/// Returns the molecule's 3D distance matrix.
///
/// ARGUMENTS:
///
///   - mol: the molecule to use
///
///   - confId: (optional) chooses the conformer Id to use
///     Default value is -1.
///
///   - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
///     matrix (to return a "Barysz" distance matrix).
///     Default value is 0.
///
///   - force: (optional) forces the calculation to proceed, even if there is a cached value.
///     Default value is 0.
///
///   - prefix: (optional, internal use) sets the prefix used in the property cache
///     Default value is .
///
///  RETURNS:
///
///    a Numeric array of floats with the distance matrix
#[pyfunction]
#[pyo3(name = "Get3DDistanceMatrix", signature = (mol, confId = -1, useAtomWts = None, force = None, prefix = ""))]
#[allow(non_snake_case)]
fn get_3d_distance_matrix<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    confId: i32,
    useAtomWts: Option<&Bound<'py, PyAny>>,
    force: Option<&Bound<'py, PyAny>>,
    prefix: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let use_atom_wts = bool_arg(useAtomWts, false)?;
    let force = bool_arg(force, false)?;
    let (n, d) = {
        let mut m = mol.borrow_mut();
        let d =
            matrices::get_3d_distance_mat_cached(&mut m.mol, confId, use_atom_wts, force, prefix)
                .map_err(|_| PyValueError::new_err("Bad Conformer Id"))?;
        (m.mol.num_atoms() as usize, d)
    };
    square_array(py, f64_bytes(&d), n, "float64")
}

/// Returns the molecule's adjacency matrix.
///
/// ARGUMENTS:
///
///   - mol: the molecule to use
///
///   - useBO: (optional) toggles use of bond orders in the matrix. Default value is 0.
///
///   - emptyVal: (optional) sets the elements of the matrix that have no bonds.
///     Default value is 0.
///
///   - force: (optional) forces the calculation to proceed, even if there is a cached value.
///     Default value is 0.
///
///   - prefix: (optional, internal use) sets the prefix used in the property cache
///     Default value is .
///
///  RETURNS:
///
///    a Numeric array of floats containing the adjacency matrix
#[pyfunction]
#[pyo3(name = "GetAdjacencyMatrix", signature = (mol, useBO = None, emptyVal = 0, force = None, prefix = ""))]
#[allow(non_snake_case)]
fn get_adjacency_matrix<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    useBO: Option<&Bound<'py, PyAny>>,
    emptyVal: i32,
    force: Option<&Bound<'py, PyAny>>,
    prefix: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let use_bo = bool_arg(useBO, false)?;
    let force = bool_arg(force, false)?;
    let (n, a) = {
        let mut m = mol.borrow_mut();
        let a = matrices::get_adjacency_matrix_cached(&mut m.mol, use_bo, emptyVal, force, prefix)
            .py()?;
        (m.mol.num_atoms() as usize, a)
    };
    if use_bo {
        square_array(py, f64_bytes(&a), n, "float64")
    } else {
        // `std::round` of each entry, as an int array
        let ints: Vec<u8> = a
            .iter()
            .flat_map(|v| (v.round() as i32).to_ne_bytes())
            .collect();
        square_array(py, ints, n, "intc")
    }
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_distance_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(get_3d_distance_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(get_adjacency_matrix, m)?)?;
    Ok(())
}
