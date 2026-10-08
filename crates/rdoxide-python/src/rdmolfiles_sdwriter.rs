//! `SDWriter` (upstream `Wrap/SDWriter.cpp`): writes to a named file,
//! stdout (`"-"`) or a Python file-like object.

use std::io::Write;

use pyo3::exceptions::{PyIOError, PyRuntimeError};
use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::sd_writer::mol_to_sd_text;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

enum Sink {
    File(std::io::BufWriter<std::fs::File>),
    Stdout,
    Py(Py<PyAny>),
}

/// A class for writing molecules to SD files.
///
///   Usage examples:
///
///     1) writing to a named file:
///
///        >>> writer = SDWriter('out.sdf')
///        >>> for mol in list_of_mols:
///        ...    writer.write(mol)
///
///     2) writing to a file-like object:
///
///        >>> import gzip
///        >>> outf=gzip.open('out.sdf.gz','wt+')
///        >>> writer = SDWriter(outf)
///        >>> for mol in list_of_mols:
///        ...   writer.write(mol)
///        >>> writer.close()
///        >>> outf.close()
///
///   By default all non-private molecular properties are written to the SD file.
///   This can be changed using the SetProps method:
///
///        >>> writer = SDWriter('out.sdf')
///        >>> writer.SetProps(['prop1','prop2'])
#[pyclass(name = "SDWriter", module = "rdkit.Chem.rdmolfiles", unsendable)]
pub struct PySDWriter {
    sink: Option<Sink>,
    props: Vec<String>,
    molid: u32,
    kekulize: bool,
    force_v3000: bool,
}

fn no_stream() -> PyErr {
    PyRuntimeError::new_err("Pre-condition Violation\n\tno output stream")
}

impl PySDWriter {
    fn emit(&mut self, py: Python<'_>, text: &str) -> PyResult<()> {
        match self.sink.as_mut().ok_or_else(no_stream)? {
            Sink::File(f) => f
                .write_all(text.as_bytes())
                .map_err(|e| PyIOError::new_err(e.to_string())),
            Sink::Stdout => std::io::stdout()
                .write_all(text.as_bytes())
                .map_err(|e| PyIOError::new_err(e.to_string())),
            Sink::Py(o) => o.call_method1(py, "write", (text,)).map(|_| ()),
        }
    }
}

#[pymethods]
impl PySDWriter {
    /// Constructor.
    ///
    ///    If a string argument is provided, it will be treated as the name of the output file.
    ///    If a file-like object is provided, output will be sent there.
    #[new]
    #[allow(non_snake_case)]
    fn new(fileName: &Bound<'_, PyAny>) -> PyResult<Self> {
        let sink = if let Ok(name) = fileName.extract::<String>() {
            if name == "-" {
                Sink::Stdout
            } else {
                let f = std::fs::File::create(&name).map_err(|_| {
                    PyIOError::new_err(format!("File error: Bad output file {name}"))
                })?;
                Sink::File(std::io::BufWriter::new(f))
            }
        } else {
            Sink::Py(fileName.clone().unbind())
        };
        Ok(PySDWriter {
            sink: Some(sink),
            props: Vec::new(),
            molid: 0,
            kekulize: true,
            force_v3000: false,
        })
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        &mut self,
        py: Python<'_>,
        exc_type: &Bound<'_, PyAny>,
        exc_val: &Bound<'_, PyAny>,
        traceback: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        let _ = (exc_type, exc_val, traceback);
        self.close(py)?;
        Ok(false)
    }

    /// Sets the properties to be written to the output file
    ///
    ///   ARGUMENTS:
    ///
    ///     - props: a list or tuple of property names
    #[pyo3(name = "SetProps", signature = (props))]
    fn set_props(&mut self, props: Vec<String>) {
        if self.molid > 0 {
            rdoxide::rd_log_warning!(
                "WARNING: Setting property list after a few molecules have been written\n"
            );
        }
        self.props = props;
    }

    /// Writes a molecule to the output file.
    ///
    ///   ARGUMENTS:
    ///
    ///     - mol: the Mol to be written
    ///     - confId: (optional) ID of the conformation to write
    #[pyo3(signature = (mol, confId = -1))]
    #[allow(non_snake_case)]
    fn write(&mut self, py: Python<'_>, mol: PyRef<'_, PyMol>, confId: i32) -> PyResult<()> {
        if self.sink.is_none() {
            return Err(no_stream());
        }
        let text = mol_to_sd_text(
            &mol.mol,
            confId,
            self.kekulize,
            self.force_v3000,
            self.molid as i32,
            Some(&self.props),
        )
        .py()?;
        self.emit(py, &text)?;
        self.molid += 1;
        Ok(())
    }

    /// Flushes the output file (forces the disk file to be updated).
    fn flush(&mut self, py: Python<'_>) -> PyResult<()> {
        match self.sink.as_mut().ok_or_else(no_stream)? {
            Sink::File(f) => {
                let _ = f.flush();
            }
            Sink::Stdout => {
                let _ = std::io::stdout().flush();
            }
            Sink::Py(o) => {
                if o.bind(py).hasattr("flush")? {
                    o.call_method0(py, "flush")?;
                }
            }
        }
        Ok(())
    }

    /// Flushes the output file and closes it. The Writer cannot be used after this.
    fn close(&mut self, py: Python<'_>) -> PyResult<()> {
        if self.sink.is_some() {
            self.flush(py)?;
        }
        self.sink = None;
        Ok(())
    }

    /// Returns the number of molecules written so far.
    #[pyo3(name = "NumMols")]
    fn num_mols(&self) -> u32 {
        self.molid
    }

    /// Sets whether or not V3000 mol file writing is being forced.
    #[pyo3(name = "SetForceV3000", signature = (val))]
    fn set_force_v3000(&mut self, val: bool) {
        self.force_v3000 = val;
    }

    /// Returns whether or not V3000 mol file writing is being forced.
    #[pyo3(name = "GetForceV3000")]
    fn get_force_v3000(&self) -> bool {
        self.force_v3000
    }

    /// Sets whether or not molecules are kekulized on writing.
    #[pyo3(name = "SetKekulize", signature = (val))]
    fn set_kekulize(&mut self, val: bool) {
        self.kekulize = val;
    }

    /// Returns whether or not molecules are kekulized on writing.
    #[pyo3(name = "GetKekulize")]
    fn get_kekulize(&self) -> bool {
        self.kekulize
    }

    /// returns the SD text for a molecule
    #[staticmethod]
    #[pyo3(name = "GetText", signature = (mol, confId = -1, kekulize = true, force_v3000 = false, molid = -1))]
    #[allow(non_snake_case)]
    fn get_text(
        mol: PyRef<'_, PyMol>,
        confId: i32,
        kekulize: bool,
        force_v3000: bool,
        molid: i32,
    ) -> PyResult<String> {
        mol_to_sd_text(&mol.mol, confId, kekulize, force_v3000, molid, None).py()
    }
}

impl Drop for PySDWriter {
    fn drop(&mut self) {
        if let Some(Sink::File(f)) = self.sink.as_mut() {
            let _ = f.flush();
        }
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySDWriter>()
}
