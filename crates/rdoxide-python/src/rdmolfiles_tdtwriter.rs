//! `TDTWriter` (upstream `Wrap/TDTWriter.cpp`): writes to a named file,
//! stdout (`"-"`) or a Python file-like object.

use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::tdt_writer::TdtWriter;

use crate::errors::IntoPyResult;
use crate::errors::{Check, Site};
use crate::mol::PyMol;
use crate::writer_sink::{WriterSink, close_sink, no_stream};

const WRITE: Site = Site {
    check: Check::Invariant,
    file: "Code/GraphMol/FileParsers/TDTWriter.cpp",
    line: 87,
    expr: "dp_ostream",
};
const FLUSH: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/FileParsers/MolWriters.h",
    line: 89,
    expr: "dp_ostream",
};

/// A class for writing molecules to TDT files.
#[pyclass(name = "TDTWriter", module = "rdoxide.Chem.rdmolfiles", unsendable)]
pub struct PyTdtWriter {
    sink: Option<WriterSink>,
    writer: TdtWriter,
}

#[pymethods]
impl PyTdtWriter {
    /// Constructor.
    ///
    ///    If a string argument is provided, it will be treated as the name of the output file.
    ///    If a file-like object is provided, output will be sent there.
    #[new]
    #[pyo3(signature = (fileObj))]
    #[allow(non_snake_case)]
    fn new(fileObj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(PyTdtWriter {
            sink: Some(WriterSink::open(fileObj)?),
            writer: TdtWriter::new(),
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
        self.writer.set_props(props);
    }

    /// Writes a molecule to the output file.
    ///
    ///   ARGUMENTS:
    ///
    ///     - mol: the Mol to be written
    ///     - confId: (optional) ID of the conformation to write
    #[pyo3(signature = (mol, confId = -1))]
    #[allow(non_snake_case)]
    fn write(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>, confId: i32) -> PyResult<()> {
        let sink = self.sink.as_mut().ok_or_else(|| no_stream(WRITE))?;
        let text = {
            let mut m = mol.borrow_mut();
            self.writer.format_mol(&mut m.mol, confId).py()?
        };
        sink.write(py, &text)
    }

    /// Flushes the output file (forces the disk file to be updated).
    fn flush(&mut self, py: Python<'_>) -> PyResult<()> {
        self.sink
            .as_mut()
            .ok_or_else(|| no_stream(FLUSH))?
            .flush(py)
    }

    /// Flushes the output file and closes it. The Writer cannot be used after this.
    fn close(&mut self, py: Python<'_>) -> PyResult<()> {
        close_sink(py, &mut self.sink)
    }

    /// Returns the number of molecules written so far.
    #[pyo3(name = "NumMols")]
    fn num_mols(&self) -> u32 {
        self.writer.num_mols()
    }

    /// causes 2D conformations to be written (default is 3D conformations)
    #[pyo3(name = "SetWrite2D", signature = (state = true))]
    fn set_write_2d(&mut self, state: bool) {
        self.writer.set_write_2d(state);
    }

    #[pyo3(name = "GetWrite2D")]
    fn get_write_2d(&self) -> bool {
        self.writer.write_2d()
    }

    /// causes names to be written to the output file as NAME records
    #[pyo3(name = "SetWriteNames", signature = (state = true))]
    fn set_write_names(&mut self, state: bool) {
        self.writer.set_write_names(state);
    }

    #[pyo3(name = "GetWriteNames")]
    fn get_write_names(&self) -> bool {
        self.writer.write_names()
    }

    /// sets the number of digits to be written for coordinates
    #[pyo3(name = "SetNumDigits", signature = (numDigits))]
    #[allow(non_snake_case)]
    fn set_num_digits(&mut self, numDigits: usize) {
        self.writer.set_num_digits(numDigits);
    }

    #[pyo3(name = "GetNumDigits")]
    fn get_num_digits(&self) -> usize {
        self.writer.num_digits()
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTdtWriter>()
}
