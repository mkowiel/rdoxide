//! `rdoxide.DataManip.Metric.rdMetricMatrixCalc` (upstream
//! `DataManip/MetricMatrixCalc/Wrap/rdMetricMatrixCalc.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use rdoxide::RDKitError;
use rdoxide::data_manip::metric_matrix::{euclidean_distance_matrix, tanimoto_matrix};

use crate::datastructs::{PyExplicitBitVect, PySparseBitVect};
use crate::errors::IntoPyResult;
use crate::pyconv::{float64_array, float64_values};

/// Threads for a matrix of `work` elementary steps: big matrices use every
/// core, small ones are not worth the thread start-up.
fn threads_for(work: usize) -> i32 {
    if work > 1 << 20 { 0 } else { 1 }
}

/// A numeric table as row-major values with its shape: a 2-D numpy array of
/// any integer or floating point type, or a sequence of sequences of numbers.
fn descriptor_table(
    py: Python<'_>,
    descrip_mat: &Bound<'_, PyAny>,
) -> PyResult<(Vec<f64>, usize, usize)> {
    let np = py.import("numpy")?;
    if descrip_mat.is_instance(&np.getattr("ndarray")?)? {
        let kind: String = descrip_mat.getattr("dtype")?.getattr("kind")?.extract()?;
        if !matches!(kind.as_str(), "i" | "u" | "f") {
            return Err(PyValueError::new_err(
                "The array has to be of type int, float, or double for GetEuclideanDistMat",
            ));
        }
        let (values, shape) = float64_values(py, descrip_mat)?;
        let &[n_rows, n_cols] = shape.as_slice() else {
            return Err(PyValueError::new_err(
                "GetEuclideanDistMat needs a two-dimensional array",
            ));
        };
        if n_rows == 0 || n_cols == 0 {
            return Err(RDKitError::invariant("empty descriptor array")).py();
        }
        return Ok((values, n_rows, n_cols));
    }
    // a list or tuple of sequences of numbers
    let n_rows = descrip_mat.len()?;
    if n_rows == 0 {
        return Err(RDKitError::invariant("Empty list passed in")).py();
    }
    let mut values = Vec::new();
    let mut n_cols = 0;
    for i in 0..n_rows {
        let row = descrip_mat.get_item(i)?;
        let row: Vec<f64> = (0..row.len()?)
            .map(|k| row.get_item(k)?.extract::<f64>())
            .collect::<PyResult<_>>()?;
        if i == 0 {
            n_cols = row.len();
        } else if row.len() != n_cols {
            return Err(PyValueError::new_err(
                "All subsequences must be the same length",
            ));
        }
        values.extend(row);
    }
    Ok((values, n_rows, n_cols))
}

/// Compute the distance matrix from a descriptor matrix using the Euclidean distance metric
///
///   ARGUMENTS:
///
///     descripMat - A python object of any one of the following types
///                    1. A numeric array of dimensions n by m where n is the number of items in the data set
///                        and m is the number of descriptors
///                    2. A list of Numeric Vectors (or 1D arrays), each entry in the list corresponds
///                        to descriptor vector for one item
///                    3. A list (or tuple) of lists (or tuples) of values, where the values can be extracted to
///                        double.
///
///   RETURNS:
///     A numeric one-dimensional array containing the lower triangle elements of the symmetric distance matrix
#[pyfunction]
#[pyo3(name = "GetEuclideanDistMat", signature = (descripMat))]
#[allow(non_snake_case)]
fn get_euclidean_dist_mat<'py>(
    py: Python<'py>,
    descripMat: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    let (values, n_rows, n_cols) = descriptor_table(py, descripMat)?;
    let threads = threads_for(n_rows * n_rows * n_cols);
    let matrix = py
        .detach(|| euclidean_distance_matrix(&values, n_rows, n_cols, threads))
        .py()?;
    float64_array(py, &matrix)
}

/// The Tanimoto matrix of a sequence of explicit or sparse bit vectors.
fn tanimoto_of_sequence<'py>(
    py: Python<'py>,
    bit_vect_list: &Bound<'py, PyAny>,
    as_distance: bool,
) -> PyResult<Bound<'py, PyAny>> {
    let n = bit_vect_list.len()?;
    if n <= 1 {
        return Err(RDKitError::invariant("at least two bit vectors are needed")).py();
    }
    let first = bit_vect_list.get_item(0)?;
    let threads = threads_for(n * n);
    let matrix = if first.extract::<PyRef<'_, PyExplicitBitVect>>().is_ok() {
        let guards = (0..n)
            .map(|i| {
                bit_vect_list
                    .get_item(i)?
                    .extract::<PyRef<'_, PyExplicitBitVect>>()
                    .map_err(PyErr::from)
            })
            .collect::<PyResult<Vec<_>>>()?;
        let items: Vec<_> = guards.iter().map(|g| &g.bv).collect();
        py.detach(|| tanimoto_matrix(&items, as_distance, threads))
            .py()?
    } else if first.extract::<PyRef<'_, PySparseBitVect>>().is_ok() {
        let guards = (0..n)
            .map(|i| {
                bit_vect_list
                    .get_item(i)?
                    .extract::<PyRef<'_, PySparseBitVect>>()
                    .map_err(PyErr::from)
            })
            .collect::<PyResult<Vec<_>>>()?;
        let items: Vec<_> = guards.iter().map(|g| &g.bv).collect();
        py.detach(|| tanimoto_matrix(&items, as_distance, threads))
            .py()?
    } else {
        return Err(PyValueError::new_err(
            "GetTanimotoDistMat can only take a sequence of ExplicitBitVects or SparseBitvects",
        ));
    };
    float64_array(py, &matrix)
}

/// Compute the distance matrix from a list of BitVects using the Tanimoto distance metric
///
///   ARGUMENTS:
///
///     bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors,
///                   needs to be expanded to support a list of SparseBitVects
///
///   RETURNS:
///     A numeric 1 dimensional array containing the lower triangle elements of the
///     symmetric distance matrix
#[pyfunction]
#[pyo3(name = "GetTanimotoDistMat", signature = (bitVectList))]
#[allow(non_snake_case)]
fn get_tanimoto_dist_mat<'py>(
    py: Python<'py>,
    bitVectList: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    tanimoto_of_sequence(py, bitVectList, true)
}

/// Compute the similarity matrix from a list of BitVects
///
///   ARGUMENTS:
///
///     bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors,
///                   needs to be expanded to support a list of SparseBitVects
///
///   RETURNS:
///     A numeric 1 dimensional array containing the lower triangle elements of the symmetric similarity matrix
#[pyfunction]
#[pyo3(name = "GetTanimotoSimMat", signature = (bitVectList))]
#[allow(non_snake_case)]
fn get_tanimoto_sim_mat<'py>(
    py: Python<'py>,
    bitVectList: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    tanimoto_of_sequence(py, bitVectList, false)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing the calculator for metric matrix calculation, \ne.g. similarity and distance matrices",
    )?;
    m.add_function(wrap_pyfunction!(get_euclidean_dist_mat, m)?)?;
    m.add_function(wrap_pyfunction!(get_tanimoto_dist_mat, m)?)?;
    m.add_function(wrap_pyfunction!(get_tanimoto_sim_mat, m)?)?;
    Ok(())
}
