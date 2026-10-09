//! Translation of rdoxide errors into the Python exceptions RDKit raises
//! (see upstream `Code/RDBoost/Wrap/RDBase.cpp` exception translators).

use pyo3::exceptions::{PyIOError, PyIndexError, PyKeyError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use rdoxide::RDKitError;

pub fn to_py(err: RDKitError) -> PyErr {
    match err {
        RDKitError::Invariant { .. } => PyRuntimeError::new_err(err.to_string()),
        RDKitError::Value(m) => PyValueError::new_err(m),
        RDKitError::Conformer(_) => PyValueError::new_err("Bad Conformer Id"),
        RDKitError::Key(m) => PyKeyError::new_err(m),
        RDKitError::Index(m) => PyIndexError::new_err(m),
        RDKitError::Runtime(m) | RDKitError::FileParse(m) => PyRuntimeError::new_err(m),
        RDKitError::BadFile(m) => PyIOError::new_err(m),
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

/// The kind of upstream check that failed.
#[derive(Clone, Copy)]
pub enum Check {
    Precondition,
    Invariant,
}

/// A failed `PRECONDITION`/`CHECK_INVARIANT` of the C++ sources. The error
/// text names the C++ source line, so each call site passes its own.
#[derive(Clone, Copy)]
pub struct Site {
    pub check: Check,
    pub file: &'static str,
    pub line: u32,
    pub expr: &'static str,
}

/// The `RuntimeError` of a failed upstream check; also logged like upstream
/// does.
pub fn check_failed(site: Site, message: &str) -> PyErr {
    let kind = match site.check {
        Check::Precondition => "Pre-condition Violation",
        Check::Invariant => "Invariant Violation",
    };
    let body = format!(
        "{kind}\n{message}\nViolation occurred on line {} in file {}\nFailed Expression: {}\n",
        site.line, site.file, site.expr
    );
    rdoxide::rd_log_error!("\n\n****\n{body}****\n");
    PyRuntimeError::new_err(format!(
        "{}\tRDKIT: {}\n\tBOOST: \n",
        body.replace('\n', "\n\t").trim_end_matches('\t'),
        rdoxide::rdgeneral::RDKIT_VERSION
    ))
}
