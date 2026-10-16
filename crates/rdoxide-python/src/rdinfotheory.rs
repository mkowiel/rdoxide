//! `rdoxide.ML.InfoTheory.rdInfoTheory` (upstream `ML/InfoTheory/Wrap/rdInfoTheory.cpp`):
//! `InfoEntropy`, `InfoGain` and `ChiSquare` of numpy count arrays. The bit
//! ranker and the correlation-matrix generator are not ported yet.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::ml::info_theory::{chi_square, info_entropy, info_entropy_gain};

use crate::errors::to_py;

/// The data of a numpy array in one of the types upstream instantiates its
/// templates for (anything else is converted to doubles, as upstream does).
enum Counts {
    F32(Vec<f32>),
    F64(Vec<f64>),
    I32(Vec<i32>),
    I64(Vec<i64>),
}

/// Runs `$body` with `$data` bound to the slice of whichever type `$counts` holds.
macro_rules! with_counts {
    ($counts:expr, |$data:ident| $body:expr) => {
        match &$counts {
            Counts::F32($data) => $body,
            Counts::F64($data) => $body,
            Counts::I32($data) => $body,
            Counts::I64($data) => $body,
        }
    };
}

fn conversion_error() -> PyErr {
    PyValueError::new_err("could not convert argument to a numeric array")
}

/// `nativeOrDoubleType` + `PyArray_ContiguousFromObject`: the contiguous data
/// of `obj` with `min_dim..=max_dim` dimensions, and its shape.
fn load_counts(
    py: Python<'_>,
    obj: &Bound<'_, PyAny>,
    min_dim: usize,
    max_dim: usize,
) -> PyResult<(Counts, Vec<usize>)> {
    let np = py.import("numpy")?;
    if !obj.is_instance(&np.getattr("ndarray")?)? {
        return Err(PyValueError::new_err("Expecting a Numeric array object"));
    }
    let dtype = obj.getattr("dtype")?;
    // NPY_FLOAT, NPY_DOUBLE, NPY_INT, NPY_LONG and NPY_LONGLONG are read as they are
    let native = match dtype.getattr("num")?.extract::<i32>()? {
        11 => Some("f32"),
        12 => Some("f64"),
        5 | 7 | 9 => Some(if dtype.getattr("itemsize")?.extract::<usize>()? == 4 {
            "i32"
        } else {
            "i64"
        }),
        _ => None,
    };
    let target = match native {
        Some(_) => dtype,
        None => {
            let double = np.getattr("float64")?;
            if !np
                .call_method1("can_cast", (&dtype, &double, "safe"))?
                .is_truthy()?
            {
                return Err(conversion_error());
            }
            double
        }
    };
    let arr = np
        .call_method1("ascontiguousarray", (obj, target))
        .map_err(|_| conversion_error())?;
    let shape: Vec<usize> = arr.getattr("shape")?.extract()?;
    if shape.len() < min_dim || shape.len() > max_dim {
        return Err(conversion_error());
    }
    let flat = arr.call_method0("ravel")?.call_method0("tolist")?;
    let counts = match native {
        Some("f32") => Counts::F32(flat.extract::<Vec<f32>>()?),
        Some("i32") => Counts::I32(flat.extract()?),
        Some("i64") => Counts::I64(flat.extract()?),
        _ => Counts::F64(flat.extract()?),
    };
    Ok((counts, shape))
}

/// calculates the informational entropy of the values in an array
///
///   ARGUMENTS:
///
///     - resMat: pointer to a long int array containing the data
///     - dim: long int containing the length of the _tPtr_ array.
///
///   RETURNS:
///
///     a double
#[pyfunction]
#[pyo3(name = "InfoEntropy", signature = (resArr))]
#[allow(non_snake_case)]
fn info_entropy_py(py: Python<'_>, resArr: &Bound<'_, PyAny>) -> PyResult<f64> {
    let (counts, shape) = load_counts(py, resArr, 1, 1)?;
    // we are expecting a 1 dimensional array
    if shape[0] == 0 {
        return Err(to_py(RDKitError::invariant("")));
    }
    Ok(with_counts!(counts, |data| info_entropy(data)))
}

/// Calculates the information gain for a variable
///
///    ARGUMENTS:
///
///      - varMat: a Numeric Array object
///        varMat is a Numeric array with the number of possible occurrences
///          of each result for reach possible value of the given variable.
///
///        So, for a variable which adopts 4 possible values and a result which
///          has 3 possible values, varMat would be 4x3
///
///    RETURNS:
///
///      - a Python float object
///
///    NOTES
///
///      - this is a dropin replacement for _PyInfoGain()_ in entropy.py
#[pyfunction]
#[pyo3(name = "InfoGain", signature = (resArr))]
#[allow(non_snake_case)]
fn info_gain_py(py: Python<'_>, resArr: &Bound<'_, PyAny>) -> PyResult<f64> {
    let (counts, shape) = load_counts(py, resArr, 2, 2)?;
    Ok(with_counts!(counts, |data| info_entropy_gain(
        data, shape[0], shape[1]
    )))
}

/// Calculates the chi squared value for a variable
///
///    ARGUMENTS:
///
///      - varMat: a Numeric Array object
///        varMat is a Numeric array with the number of possible occurrences
///          of each result for reach possible value of the given variable.
///
///        So, for a variable which adopts 4 possible values and a result which
///          has 3 possible values, varMat would be 4x3
///
///    RETURNS:
///
///      - a Python float object
#[pyfunction]
#[pyo3(name = "ChiSquare", signature = (resArr))]
#[allow(non_snake_case)]
fn chi_square_py(py: Python<'_>, resArr: &Bound<'_, PyAny>) -> PyResult<f64> {
    let (counts, shape) = load_counts(py, resArr, 2, 2)?;
    Ok(with_counts!(counts, |data| chi_square(
        data, shape[0], shape[1]
    )))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing bunch of functions for information metrics and a ranker to rank bits",
    )?;
    m.add_function(wrap_pyfunction!(info_entropy_py, m)?)?;
    m.add_function(wrap_pyfunction!(info_gain_py, m)?)?;
    m.add_function(wrap_pyfunction!(chi_square_py, m)?)
}
