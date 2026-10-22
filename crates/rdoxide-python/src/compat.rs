//! Types of `rdoxide._compat`: tuples that also equal lists and exceptions
//! that derive from the exception classes of both RDKit wrapper flavours.

use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::PyTuple;

fn compat_attr<'py>(
    py: Python<'py>,
    cell: &'static PyOnceLock<Py<PyAny>>,
    name: &str,
) -> PyResult<Bound<'py, PyAny>> {
    Ok(cell
        .get_or_try_init(py, || {
            Ok::<_, PyErr>(py.import("rdoxide._compat")?.getattr(name)?.unbind())
        })?
        .bind(py)
        .clone())
}

/// A `TupleList` (a tuple that equals the same-valued list) with these items.
pub fn tuple_list<'py>(
    py: Python<'py>,
    items: Bound<'py, PyTuple>,
) -> PyResult<Bound<'py, PyTuple>> {
    static CLS: PyOnceLock<Py<PyAny>> = PyOnceLock::new();
    Ok(compat_attr(py, &CLS, "TupleList")?
        .call1((items,))?
        .cast_into::<PyTuple>()?)
}

macro_rules! dual_error {
    ($fn:ident, $name:literal) => {
        pub fn $fn(py: Python<'_>, msg: impl Into<String>) -> PyErr {
            static CLS: PyOnceLock<Py<PyAny>> = PyOnceLock::new();
            match compat_attr(py, &CLS, $name) {
                Ok(cls) => match cls.call1((msg.into(),)) {
                    Ok(e) => PyErr::from_value(e),
                    Err(e) => e,
                },
                Err(e) => e,
            }
        }
    };
}

dual_error!(overflow_type_error, "OverflowTypeError");
dual_error!(value_type_error, "ValueTypeError");
dual_error!(runtime_type_error, "RuntimeTypeError");
