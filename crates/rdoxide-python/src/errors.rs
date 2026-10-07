//! Translation of rdoxide errors into the Python exceptions RDKit raises
//! (see upstream `Code/RDBoost/Wrap/RDBase.cpp` exception translators).

use pyo3::exceptions::{PyIndexError, PyKeyError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use rdoxide::RDKitError;

pub fn to_py(err: RDKitError) -> PyErr {
    match err {
        RDKitError::Invariant { .. } => PyRuntimeError::new_err(err.to_string()),
        RDKitError::Value(m) => PyValueError::new_err(m),
        RDKitError::Key(m) => PyKeyError::new_err(m),
        RDKitError::Index(m) => PyIndexError::new_err(m),
        RDKitError::Sanitize(e) => crate::rdchem::sanitize_error(&e),
    }
}

/// `result.map_err(to_py)` shorthand.
pub trait IntoPyResult<T> {
    fn py(self) -> PyResult<T>;
}

impl<T> IntoPyResult<T> for Result<T, RDKitError> {
    fn py(self) -> PyResult<T> {
        self.map_err(to_py)
    }
}
