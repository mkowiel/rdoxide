//! Translation of rdoxide errors into the Python exceptions RDKit raises
//! (see upstream `Code/RDBoost/Wrap/RDBase.cpp` exception translators).

use pyo3::exceptions::{
    PyIOError, PyIndexError, PyKeyError, PyNotImplementedError, PyRuntimeError, PyTypeError,
    PyValueError,
};
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyList};
use rdoxide::RDKitError;

pub fn to_py(err: RDKitError) -> PyErr {
    match err {
        RDKitError::Invariant { .. } => PyRuntimeError::new_err(err.to_string()),
        RDKitError::Value(m) => PyValueError::new_err(m),
        RDKitError::Conformer(_) => PyValueError::new_err("Bad Conformer Id"),
        RDKitError::Key(m) => PyKeyError::new_err(m),
        RDKitError::Index(m) => PyIndexError::new_err(m),
        RDKitError::Runtime(m) | RDKitError::FileParse(m) | RDKitError::UnhandledFeature(m) => {
            PyRuntimeError::new_err(m)
        }
        RDKitError::BadFile(m) => PyIOError::new_err(m),
        RDKitError::ChemicalReaction(m) => {
            PyValueError::new_err(format!("ChemicalParserException: {m}"))
        }
        RDKitError::ChemicalReactionParser(m) => {
            PyValueError::new_err(format!("ChemicalReactionParserException: {m}"))
        }
        RDKitError::Sanitize(e) => crate::rdchem::sanitize_error(&e),
        RDKitError::NotImplemented(m) => PyNotImplementedError::new_err(m),
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

fn type_name(obj: &Bound<'_, PyAny>) -> String {
    obj.get_type()
        .name()
        .map_or_else(|_| "?".into(), |n| n.to_string())
}

/// `python::extract<python::list>(obj)` of the Boost wrappers: the list, or
/// the `TypeError` Boost raises for other objects.
pub fn expect_list<'a, 'py>(obj: &'a Bound<'py, PyAny>) -> PyResult<&'a Bound<'py, PyList>> {
    obj.cast::<PyList>().map_err(|_| {
        PyTypeError::new_err(format!(
            "Expecting an object of type list; got an object of type {} instead",
            type_name(obj)
        ))
    })
}

/// `python::extract<python::dict>(obj)` of the Boost wrappers: the dict, or
/// the `TypeError` Boost raises for other objects.
pub fn expect_dict<'a, 'py>(obj: &'a Bound<'py, PyAny>) -> PyResult<&'a Bound<'py, PyDict>> {
    obj.cast::<PyDict>().map_err(|_| {
        PyTypeError::new_err(format!(
            "Expecting an object of type dict; got an object of type {} instead",
            type_name(obj)
        ))
    })
}

/// The `TypeError` Boost raises when no converter produces the C++ type
/// `cpp_type` from `obj`.
pub fn no_converter(cpp_type: &str, obj: &Bound<'_, PyAny>) -> PyErr {
    PyTypeError::new_err(format!(
        "No registered converter was able to produce a C++ rvalue of type {cpp_type} from this Python object of type {}",
        type_name(obj)
    ))
}

/// `Boost.Python.ArgumentError`, a `TypeError` subclass.
fn argument_error_type(py: Python<'_>) -> PyResult<&Bound<'_, pyo3::types::PyType>> {
    static TYPE: PyOnceLock<Py<pyo3::types::PyType>> = PyOnceLock::new();
    TYPE.get_or_try_init(py, || {
        PyErr::new_type(
            py,
            c"Boost.Python.ArgumentError",
            None,
            Some(&py.get_type::<PyTypeError>()),
            None,
        )
    })
    .map(|t| t.bind(py))
}

/// The `Boost.Python.ArgumentError` raised when no overload of `func`
/// matches the call. `got` is the Python-side signature (`Mol, Mol, tuple`),
/// `expected` the C++ one.
pub fn argument_error(py: Python<'_>, func: &str, got: &[String], expected: &str) -> PyErr {
    argument_error_in(py, "None", func, got, expected)
}

/// [`argument_error`] for a method; `scope` is the Python class name.
pub fn argument_error_in(
    py: Python<'_>,
    scope: &str,
    func: &str,
    got: &[String],
    expected: &str,
) -> PyErr {
    let msg = format!(
        "Python argument types in\n    {scope}.{func}({})\ndid not match C++ signature:\n    {expected}",
        got.join(", ")
    );
    match argument_error_type(py) {
        Ok(t) => PyErr::from_type(t.clone(), msg),
        Err(e) => e,
    }
}

/// Name of the Python type of `obj` as Boost prints it in signatures.
pub fn boost_type_name(obj: &Bound<'_, PyAny>) -> String {
    let n = type_name(obj);
    if obj.is_instance_of::<pyo3::types::PyString>() {
        "str".into()
    } else {
        n
    }
}
