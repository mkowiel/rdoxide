//! Output target shared by the text writers (`SDWriter`, `SmilesWriter`): a
//! named file, stdout (`"-"`) or a Python file-like object.

use std::io::Write;

use pyo3::exceptions::{PyIOError, PyRuntimeError};
use pyo3::prelude::*;

enum Target {
    File(std::io::BufWriter<std::fs::File>),
    Stdout,
    Py(Py<PyAny>),
}

pub struct WriterSink(Target);

/// The kind of upstream check that fails when a closed writer is used.
#[derive(Clone, Copy)]
pub enum Check {
    Precondition,
    Invariant,
}

/// Where upstream checks `dp_ostream`: the error text names the C++ source
/// line, so each call site passes its own.
#[derive(Clone, Copy)]
pub struct Site {
    pub check: Check,
    pub file: &'static str,
    pub line: u32,
}

/// The error of using a closed writer; also logged like upstream does.
pub fn no_stream(site: Site) -> PyErr {
    let kind = match site.check {
        Check::Precondition => "Pre-condition Violation",
        Check::Invariant => "Invariant Violation",
    };
    let body = format!(
        "{kind}\nno output stream\nViolation occurred on line {} in file {}\nFailed Expression: dp_ostream\n",
        site.line, site.file
    );
    rdoxide::rd_log_error!("\n\n****\n{body}****\n");
    PyRuntimeError::new_err(format!(
        "{}\tRDKIT: {}\n\tBOOST: \n",
        body.replace('\n', "\n\t").trim_end_matches('\t'),
        rdoxide::rdgeneral::RDKIT_VERSION
    ))
}

impl WriterSink {
    /// A string names the output file (`"-"` is stdout); anything else is
    /// taken as a file-like object with a `write` method.
    pub fn open(target: &Bound<'_, PyAny>) -> PyResult<Self> {
        let Ok(name) = target.extract::<String>() else {
            return Ok(WriterSink(Target::Py(target.clone().unbind())));
        };
        if name == "-" {
            return Ok(WriterSink(Target::Stdout));
        }
        let f = std::fs::File::create(&name)
            .map_err(|_| PyIOError::new_err(format!("File error: Bad output file {name}")))?;
        Ok(WriterSink(Target::File(std::io::BufWriter::new(f))))
    }

    pub fn write(&mut self, py: Python<'_>, text: &str) -> PyResult<()> {
        match &mut self.0 {
            Target::File(f) => f.write_all(text.as_bytes()),
            Target::Stdout => std::io::stdout().write_all(text.as_bytes()),
            Target::Py(o) => return o.call_method1(py, "write", (text,)).map(|_| ()),
        }
        .map_err(|e| PyIOError::new_err(e.to_string()))
    }

    pub fn flush(&mut self, py: Python<'_>) -> PyResult<()> {
        match &mut self.0 {
            Target::File(f) => {
                let _ = f.flush();
            }
            Target::Stdout => {
                let _ = std::io::stdout().flush();
            }
            Target::Py(o) => {
                if o.bind(py).hasattr("flush")? {
                    o.call_method0(py, "flush")?;
                }
            }
        }
        Ok(())
    }
}

impl Drop for WriterSink {
    fn drop(&mut self) {
        if let Target::File(f) = &mut self.0 {
            let _ = f.flush();
        }
    }
}

/// Shared by the writers' `close`: flushes and drops an open sink.
pub fn close_sink(py: Python<'_>, sink: &mut Option<WriterSink>) -> PyResult<()> {
    if let Some(s) = sink.as_mut() {
        s.flush(py)?;
    }
    *sink = None;
    Ok(())
}
