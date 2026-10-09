//! Conversions of Python arguments shared by the bindings (upstream
//! `RDBoost/Wrap.h`).

use pyo3::exceptions::{PyAttributeError, PyOverflowError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyInt;

/// Boost's `unsigned int` converter: an int, or the errors Boost raises for
/// negative numbers and for objects that are not integers.
pub fn extract_unsigned(obj: &Bound<'_, PyAny>) -> PyResult<u32> {
    if obj.is_instance_of::<PyInt>() {
        if let Ok(v) = obj.extract::<u32>() {
            return Ok(v);
        }
        return Err(PyOverflowError::new_err(
            if obj.extract::<i64>().is_ok_and(|v| v < 0) {
                "can't convert negative value to unsigned int"
            } else {
                "Python int too large to convert to C long"
            },
        ));
    }
    Err(crate::errors::no_converter("unsigned int", obj))
}

/// `pythonObjectToVect<unsigned>(obj, maxV)`: the elements of a sequence, or
/// `None` when `obj` is missing or falsy (None, 0, an empty sequence). With
/// `max_v`, elements that are not below it are rejected.
pub fn python_object_to_vect(
    obj: Option<&Bound<'_, PyAny>>,
    max_v: Option<u32>,
) -> PyResult<Option<Vec<u32>>> {
    let Some(obj) = obj else { return Ok(None) };
    if !obj.is_truthy()? {
        return Ok(None);
    }
    if !obj.hasattr("__iter__")? {
        return Err(PyAttributeError::new_err(format!(
            "'{}' object has no attribute '__iter__'",
            obj.get_type().name()?
        )));
    }
    let mut res = Vec::new();
    for item in obj.try_iter()? {
        let v = extract_unsigned(&item?)?;
        if max_v.is_some_and(|max| v >= max) {
            return Err(PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        res.push(v);
    }
    Ok(Some(res))
}

/// A NumPy `float64` array holding `values`, built from raw bytes so that
/// large results do not pass through Python lists.
pub fn float64_array<'py>(py: Python<'py>, values: &[f64]) -> PyResult<Bound<'py, PyAny>> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_ne_bytes()).collect();
    // a bytearray gives a writable array without another copy
    let buffer = pyo3::types::PyByteArray::new(py, &bytes);
    let np = py.import("numpy")?;
    np.call_method1("frombuffer", (buffer, np.getattr("float64")?))
}

/// The values of `array` converted to `float64`, in C order, and its shape.
pub fn float64_values(
    py: Python<'_>,
    array: &Bound<'_, PyAny>,
) -> PyResult<(Vec<f64>, Vec<usize>)> {
    let np = py.import("numpy")?;
    let contiguous = np.call_method1("ascontiguousarray", (array, "float64"))?;
    let shape: Vec<usize> = contiguous.getattr("shape")?.extract()?;
    let bytes: Vec<u8> = contiguous.call_method0("tobytes")?.extract()?;
    let values = bytes
        .chunks_exact(8)
        .map(|c| f64::from_ne_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]))
        .collect();
    Ok((values, shape))
}
