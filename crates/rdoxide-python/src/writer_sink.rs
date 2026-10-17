//! Output target shared by the text writers (`SDWriter`, `SmilesWriter`): a
//! named file, stdout (`"-"`) or a Python file-like object.

use std::io::Write;

use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;

use crate::errors::Site;
use crate::errors::check_failed;

enum Target {
    File(std::io::BufWriter<std::fs::File>),
    Stdout,
    Py(Py<PyAny>),
}

pub struct WriterSink(Target);

/// The error of using a closed writer.
pub fn no_stream(site: Site) -> PyErr {
    check_failed(site, "no output stream")
}

impl WriterSink {
    /// A string (or path) names the output file (`"-"` is stdout); anything else is
    /// taken as a file-like object with a `write` method.
    pub fn open(target: &Bound<'_, PyAny>) -> PyResult<Self> {
        let name = match target.extract::<String>() {
            Ok(name) => name,
            Err(_) if target.hasattr("__fspath__")? => crate::pyconv::path_like(target)?,
            Err(_) => {
                // upstream's streambuf insists on a text-mode file object
                let text_base = target.py().import("io")?.getattr("TextIOBase")?;
                if !target.is_instance(&text_base)? {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "Need a text mode file object like StringIO or a file opened with mode 't'",
                    ));
                }
                return Ok(WriterSink(Target::Py(target.clone().unbind())));
            }
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
