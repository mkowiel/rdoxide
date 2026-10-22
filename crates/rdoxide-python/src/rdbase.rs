//! `rdoxide.rdBase` — port of upstream `Code/RDBoost/Wrap/RDBase.cpp`.

use std::sync::Arc;

use pyo3::prelude::*;
use rdoxide::rdgeneral::log::{self, CaptureLog, Level, LogSink, LogStateSetter};

/// Sends each line to `sys.stderr` (upstream `PyErrStream`).
struct PyStderrSink;
impl LogSink for PyStderrSink {
    fn write(&self, _level: Level, text: &str) {
        Python::attach(|py| {
            if let Ok(stderr) = py.import("sys").and_then(|s| s.getattr("stderr")) {
                for line in text.lines() {
                    let _ = stderr.call_method1("write", (format!("{line}\n"),));
                }
            }
        });
    }
}

/// Sends each line to `logging.getLogger("rdoxide").<level>` (upstream `PyLogStream`).
struct PyLoggerSink;
impl LogSink for PyLoggerSink {
    fn write(&self, level: Level, text: &str) {
        let method = match level {
            Level::Debug => "debug",
            Level::Info => "info",
            Level::Warning => "warning",
            Level::Error => "error",
        };
        Python::attach(|py| {
            let logger = py
                .import("logging")
                .and_then(|l| l.call_method1("getLogger", ("rdoxide",)));
            match logger {
                Ok(logger) => {
                    for line in text.lines() {
                        if let Err(e) = logger.call_method1(method, (line,)) {
                            e.print(py);
                        }
                    }
                }
                Err(e) => e.print(py),
            }
        });
    }
}

#[pyfunction]
#[pyo3(name = "LogToCppStreams")]
/// Initialize RDKit logs with C++ streams
fn log_to_cpp_streams() {
    log::init_logs();
}

#[pyfunction]
#[pyo3(name = "LogToPythonLogger")]
/// Initialize RDKit logs with Python's logging module
fn log_to_python_logger() {
    log::set_destinations(|_| Arc::new(PyLoggerSink));
}

#[pyfunction]
#[pyo3(name = "LogToPythonStderr")]
/// Initialize RDKit logs with Python's stderr stream
fn log_to_python_stderr() {
    log::set_destinations(|_| Arc::new(PyStderrSink));
}

#[pyfunction]
#[pyo3(name = "WrapLogs")]
/// Tee the RDKit logs to Python's stderr stream
fn wrap_logs() {
    log::set_tees(|_| Arc::new(PyStderrSink));
}

#[pyfunction]
#[pyo3(name = "EnableLog")]
fn enable_log(spec: &str) {
    log::enable_logs(spec);
}

#[pyfunction]
#[pyo3(name = "DisableLog")]
fn disable_log(spec: &str) {
    log::disable_logs(spec);
}

#[pyfunction]
#[pyo3(name = "LogStatus")]
fn log_status() -> String {
    log::log_status()
}

#[pyfunction]
#[pyo3(name = "LogDebugMsg")]
/// Log a message to the RDKit debug logs
fn log_debug_msg(py: Python<'_>, msg: &str) {
    py.detach(|| log::log(Level::Debug, msg));
}

#[pyfunction]
#[pyo3(name = "LogInfoMsg")]
/// Log a message to the RDKit info logs
fn log_info_msg(py: Python<'_>, msg: &str) {
    py.detach(|| log::log(Level::Info, msg));
}

#[pyfunction]
#[pyo3(name = "LogWarningMsg")]
/// Log a message to the RDKit warning logs
fn log_warning_msg(py: Python<'_>, msg: &str) {
    py.detach(|| log::log(Level::Warning, msg));
}

#[pyfunction]
#[pyo3(name = "LogErrorMsg")]
/// Log a message to the RDKit error logs
fn log_error_msg(py: Python<'_>, msg: &str) {
    py.detach(|| log::log(Level::Error, msg));
}

#[pyfunction]
#[pyo3(name = "LogMessage")]
/// Log a message to any rdApp.* log
fn log_message(py: Python<'_>, spec: &str, msg: &str) {
    if let Some(level) = Level::from_spec(spec) {
        py.detach(|| log::log(level, msg));
    }
}

#[pyfunction]
#[pyo3(name = "AttachFileToLog", signature = (spec, filename, delay = 100))]
/// Causes the log to write to a file
fn attach_file_to_log(spec: &str, filename: &str, delay: i32) {
    // no-op upstream as well
    let _ = (spec, filename, delay);
}

#[pyfunction]
#[pyo3(name = "SeedRandomNumberGenerator")]
/// Provides a seed to the standard C random number generator
/// This does not affect pure Python code, but is relevant to some of the RDKit C++ components.
fn seed_random_number_generator(seed: u32) {
    rdoxide::rdgeneral::random::seed(seed);
}

#[pyfunction]
#[pyo3(name = "_version")]
/// Deprecated, use the constant rdkitVersion instead
fn version() -> &'static str {
    rdoxide::rdgeneral::RDKIT_VERSION
}

/// Temporarily block logs from outputting while this instance is in scope.
#[pyclass(name = "BlockLogs", module = "rdoxide.rdBase", unsendable)]
struct BlockLogs {
    setter: Option<LogStateSetter>,
}

#[pymethods]
impl BlockLogs {
    #[new]
    fn new() -> Self {
        BlockLogs {
            setter: Some(LogStateSetter::new()),
        }
    }
    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }
    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        &mut self,
        exc_type: Bound<'_, PyAny>,
        exc_val: Bound<'_, PyAny>,
        traceback: Bound<'_, PyAny>,
    ) {
        let _ = (exc_type, exc_val, traceback);
        self.setter = None;
    }
}

/// Captures messages from rdErrorLog while this instance is in scope.
/// Can be used as a context manager. The ``messages`` property is
/// accessible both inside the context and after it exits.
/// Nesting is supported: inner captures shadow outer ones.
#[pyclass(name = "CaptureErrorLog", module = "rdoxide.rdBase", unsendable)]
struct PyCaptureErrorLog {
    capture: CaptureLog,
}

#[pymethods]
impl PyCaptureErrorLog {
    #[new]
    fn new() -> Self {
        PyCaptureErrorLog {
            capture: CaptureLog::error(),
        }
    }
    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }
    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        &mut self,
        exc_type: Bound<'_, PyAny>,
        exc_val: Bound<'_, PyAny>,
        traceback: Bound<'_, PyAny>,
    ) {
        let _ = (exc_type, exc_val, traceback);
        self.capture.finish();
    }
    /// Messages captured from rdErrorLog.
    #[getter]
    fn messages(&self) -> String {
        self.capture.messages()
    }
}

/// Returns the maximum number of threads any rdoxide operation uses.
///
/// It comes from the RDOXIDE_NUM_THREADS environment variable (a number, or
/// "all"; a negative number means all but that many) and defaults to all
/// hardware threads. Functions with a numThreads argument use the maximum
/// when it is 0 and never exceed it.
#[pyfunction]
#[pyo3(name = "GetMaxThreads")]
fn get_max_threads() -> usize {
    rdoxide::rdgeneral::threads::max_threads()
}

/// Sets the maximum number of threads any rdoxide operation uses: a positive
/// number is that many threads, 0 all hardware threads and a negative number
/// all but that many.
#[pyfunction]
#[pyo3(name = "SetMaxThreads", signature = (numThreads))]
#[allow(non_snake_case)]
fn set_max_threads(numThreads: i32) {
    rdoxide::rdgeneral::threads::set_max_threads(numThreads);
}

/// Wraps a Python file object (upstream `rdBase.streambuf`); suppliers and
/// writers only need its `read`/`write`/`flush`, which are forwarded.
#[pyclass(name = "streambuf", module = "rdoxide.rdBase")]
#[derive(Debug)]
pub struct PyStreambuf {
    file: Py<PyAny>,
}

#[pymethods]
impl PyStreambuf {
    #[new]
    #[pyo3(signature = (python_file_obj, buffer_size = 0))]
    fn new(python_file_obj: Py<PyAny>, buffer_size: usize) -> Self {
        let _ = buffer_size;
        Self {
            file: python_file_obj,
        }
    }

    fn __getattr__(&self, py: Python<'_>, name: &str) -> PyResult<Py<PyAny>> {
        Ok(self.file.bind(py).getattr(name)?.unbind())
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyStreambuf>()?;
    m.add(
        "__doc__",
        "Module containing basic definitions for wrapped C++ code\n\n",
    )?;
    log::init_logs();
    crate::vects::add_classes(m)?;
    m.add("rdkitVersion", rdoxide::rdgeneral::RDKIT_VERSION)?;
    m.add("boostVersion", "")?;
    m.add(
        "rdkitBuild",
        // upstream: "<system>|<processor>|<compiler>|<build type>"; tests read the compiler field
        format!(
            "{}|{}|rustc|Release",
            std::env::consts::OS,
            std::env::consts::ARCH
        ),
    )?;
    m.add("rdoxideVersion", rdoxide::rdgeneral::RDOXIDE_VERSION)?;
    // PyO3 behaves like nanobind (TypeError on bad arguments, builtin
    // functions do not bind as methods), so upstream Python code must take
    // its nanobind code paths.
    m.add("_wrapperType", "nanobind")?;
    m.add("_serializationEnabled", true)?;
    m.add("_iostreamsEnabled", false)?;
    m.add("_multithreadedEnabled", true)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(get_max_threads, m)?)?;
    m.add_function(wrap_pyfunction!(set_max_threads, m)?)?;
    m.add_function(wrap_pyfunction!(log_to_cpp_streams, m)?)?;
    m.add_function(wrap_pyfunction!(log_to_python_logger, m)?)?;
    m.add_function(wrap_pyfunction!(log_to_python_stderr, m)?)?;
    m.add_function(wrap_pyfunction!(wrap_logs, m)?)?;
    m.add_function(wrap_pyfunction!(enable_log, m)?)?;
    m.add_function(wrap_pyfunction!(disable_log, m)?)?;
    m.add_function(wrap_pyfunction!(log_status, m)?)?;
    m.add_function(wrap_pyfunction!(log_debug_msg, m)?)?;
    m.add_function(wrap_pyfunction!(log_info_msg, m)?)?;
    m.add_function(wrap_pyfunction!(log_warning_msg, m)?)?;
    m.add_function(wrap_pyfunction!(log_error_msg, m)?)?;
    m.add_function(wrap_pyfunction!(log_message, m)?)?;
    m.add_function(wrap_pyfunction!(attach_file_to_log, m)?)?;
    m.add_function(wrap_pyfunction!(seed_random_number_generator, m)?)?;
    m.add_class::<BlockLogs>()?;
    m.add_class::<PyCaptureErrorLog>()?;
    Ok(())
}
