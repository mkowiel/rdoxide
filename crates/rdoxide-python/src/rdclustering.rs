//! `rdoxide.ML.Cluster.Clustering` (upstream `ML/Cluster/Murtagh/Clustering.cpp`):
//! Murtagh's hierarchical clustering of points or of precomputed distances.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::ml::cluster::murtagh::{Dendrogram, Linkage, cluster_points, hierarchical_cluster};

use crate::errors::IntoPyResult;
use crate::pyconv::float64_flat;

/// A float64 C-contiguous copy of the numpy array `data`, flattened, its shape
/// and whether `data` itself was such an array. Anything that is not a numpy array is rejected like upstream.
fn float_array<'py>(
    py: Python<'py>,
    data: &Bound<'py, PyAny>,
    ndim: usize,
) -> PyResult<(Vec<f64>, Vec<usize>, bool)> {
    let np = py.import("numpy")?;
    if !data.is_instance(&np.getattr("ndarray")?)? {
        return Err(PyValueError::new_err("PyArray_Type expected as input"));
    }
    let arr = np.call_method1("ascontiguousarray", (data, "float64"))?;
    let shape: Vec<usize> = arr.getattr("shape")?.extract()?;
    if shape.len() != ndim {
        return Err(PyValueError::new_err(format!(
            "object of wrong depth for desired array: expected {ndim} dimensions"
        )));
    }
    // when no conversion was needed the clustering works on the caller's array
    let is_argument = arr.is(data);
    Ok((float64_flat(&arr)?, shape, is_argument))
}

fn linkage(option: i64) -> PyResult<Linkage> {
    Linkage::from_option(option)
        .ok_or_else(|| PyValueError::new_err(format!("unknown clustering option {option}")))
}

/// The `(ia, ib, crit)` numpy arrays of a clustering history.
fn history<'py>(py: Python<'py>, d: Dendrogram) -> PyResult<Bound<'py, PyTuple>> {
    let np = py.import("numpy")?;
    let array = |values: Bound<'py, PyAny>, dtype: &str| -> PyResult<Bound<'py, PyAny>> {
        np.call_method1("array", (values, dtype))
    };
    PyTuple::new(
        py,
        [
            array(d.ia.into_pyobject(py)?.into_any(), "l")?,
            array(d.ib.into_pyobject(py)?.into_any(), "l")?,
            array(d.crit.into_pyobject(py)?.into_any(), "d")?,
        ],
    )
}

/// TODO: provide docstring
#[pyfunction]
#[pyo3(name = "MurtaghCluster", signature = (data, nPts, sz, option))]
#[allow(non_snake_case)]
fn murtagh_cluster<'py>(
    py: Python<'py>,
    data: &Bound<'py, PyAny>,
    nPts: usize,
    sz: usize,
    option: i64,
) -> PyResult<Bound<'py, PyTuple>> {
    let (points, shape, _) = float_array(py, data, 2)?;
    if nPts > shape[0] || sz > shape[1] {
        return Err(PyValueError::new_err(
            "nPts and sz must not exceed the shape of the data",
        ));
    }
    let method = linkage(option)?;
    // rows are `shape[1]` values wide, however many columns are used
    let packed: Vec<f64> = (0..nPts)
        .flat_map(|i| points[i * shape[1]..i * shape[1] + sz].iter().copied())
        .collect();
    // the pairwise distances dominate for big inputs: use all cores
    let threads = if nPts * nPts * sz.max(1) > 1 << 20 {
        0
    } else {
        1
    };
    let d = py
        .detach(|| cluster_points(&packed, nPts, sz, method, threads))
        .py()?;
    history(py, d)
}

/// TODO: provide docstring
#[pyfunction]
#[pyo3(name = "MurtaghDistCluster", signature = (data, nPts, option))]
#[allow(non_snake_case)]
fn murtagh_dist_cluster<'py>(
    py: Python<'py>,
    data: &Bound<'py, PyAny>,
    nPts: usize,
    option: i64,
) -> PyResult<Bound<'py, PyTuple>> {
    let (mut dists, _, is_argument) = float_array(py, data, 1)?;
    let method = linkage(option)?;
    let d = py
        .detach(|| hierarchical_cluster(nPts, method, &mut dists))
        .py()?;
    if is_argument {
        // the clustering works in place: the caller's array is left holding
        // the updated dissimilarities, as with upstream
        let np = py.import("numpy")?;
        np.call_method1("copyto", (data, np.call_method1("array", (dists,))?))?;
    }
    history(py, d)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(murtagh_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(murtagh_dist_cluster, m)?)
}
