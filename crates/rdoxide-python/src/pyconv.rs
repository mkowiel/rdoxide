//! Conversions of Python arguments shared by the bindings (upstream
//! `RDBoost/Wrap.h`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

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
    let mut res = Vec::new();
    for item in obj.try_iter()? {
        let v: u32 = item?.extract()?;
        if max_v.is_some_and(|max| v >= max) {
            return Err(PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        res.push(v);
    }
    Ok(Some(res))
}
