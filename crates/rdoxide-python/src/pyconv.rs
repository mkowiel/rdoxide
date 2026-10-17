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

/// Sets the declared attribute `name` of `obj` through its data descriptor on
/// the class, and rejects every other name as upstream's `safeSetattr` does.
///
/// A `__setattr__` of a `#[pyclass]` cannot delegate to `object.__setattr__`:
/// before Python 3.13 that call refuses types whose `tp_setattro` is not
/// the inherited one ("can't apply this __setattr__").
pub fn set_declared_attribute(
    obj: &Bound<'_, PyAny>,
    name: &str,
    value: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let descriptor = obj
        .get_type()
        .getattr(name)
        .map_err(|_| PyAttributeError::new_err(format!("Cannot set unknown attribute '{name}'")))?;
    if !descriptor.hasattr("__set__")? {
        return Err(PyAttributeError::new_err(format!(
            "attribute '{name}' of '{}' object is read-only",
            obj.get_type().name()?
        )));
    }
    descriptor.call_method1("__set__", (obj, value))?;
    Ok(())
}

/// `from_py_with` converter for `bool` arguments: boost::python (and so the
/// upstream API) accepts any truthy/falsy object, such as `0` and `1`.
pub fn truthy(obj: &Bound<'_, PyAny>) -> PyResult<bool> {
    obj.is_truthy()
}

/// `from_py_with` converter for file names: a `str` or any `os.PathLike`.
pub fn path_like(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    obj.py()
        .import("os")?
        .call_method1("fspath", (obj,))?
        .extract::<String>()
}

/// `from_py_with` converter for optional file names (see [`path_like`]).
pub fn path_like_opt(obj: &Bound<'_, PyAny>) -> PyResult<Option<String>> {
    if obj.is_none() {
        Ok(None)
    } else {
        path_like(obj).map(Some)
    }
}

/// `from_py_with` converter for text that may be given as `str` or `bytes`.
pub fn text_bytes(obj: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(b) = obj.extract::<Vec<u8>>() {
        if obj.is_instance_of::<pyo3::types::PyBytes>() {
            return Ok(b);
        }
    }
    Ok(obj.extract::<String>()?.into_bytes())
}

/// `from_py_with` converter for text given as `str` or `bytes`
/// (`pyObjectToString`).
pub fn text_string(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    Ok(String::from_utf8_lossy(&text_bytes(obj)?).into_owned())
}
