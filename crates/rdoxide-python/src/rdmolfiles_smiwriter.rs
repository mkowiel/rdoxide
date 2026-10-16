//! `SmilesWriter` (upstream `Wrap/SmilesWriter.cpp`): writes to a named file,
//! stdout (`"-"`) or a Python file-like object.

use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::smiles_writer::SmilesWriter;

use crate::errors::IntoPyResult;
use crate::errors::{Check, Site};
use crate::mol::PyMol;
use crate::writer_sink::{WriterSink, close_sink, no_stream};

const FILE: &str = "Code/GraphMol/FileParsers/SmilesWriter.cpp";
const WRITE: Site = Site {
    check: Check::Invariant,
    file: FILE,
    line: 116,
    expr: "dp_ostream",
};
const FLUSH: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/FileParsers/MolWriters.h",
    line: 89,
    expr: "dp_ostream",
};

/// A class for writing molecules to text files.
#[pyclass(name = "SmilesWriter", module = "rdoxide.Chem.rdmolfiles", unsendable)]
pub struct PySmilesWriter {
    sink: Option<WriterSink>,
    writer: SmilesWriter,
}

#[pymethods]
impl PySmilesWriter {
    /// Constructor.
    ///
    ///    ARGUMENTS:
    ///
    ///      - fileName: name of the output file. ('-' to write to stdout)
    ///      - delimiter: (optional) delimiter to be used to separate entries on each line.
    ///      - nameHeader: (optional) text to use for the name column in the header line.
    ///                    If this is blank, names will not be included in the output.
    ///      - includeHeader: (optional) toggles inclusion of a header line in the output file.
    ///      - isomericSmiles: (optional) toggles output of isomeric smiles (includes stereochem information).
    ///      - kekuleSmiles: (optional) toggles output of kekule smiles (no aromatic bonds for molecules that have been kekulized).
    #[new]
    #[pyo3(signature = (fileObj, delimiter = " ", nameHeader = "Name", includeHeader = true, isomericSmiles = true, kekuleSmiles = false))]
    #[allow(non_snake_case)]
    fn new(
        fileObj: &Bound<'_, PyAny>,
        delimiter: &str,
        nameHeader: &str,
        includeHeader: bool,
        isomericSmiles: bool,
        kekuleSmiles: bool,
    ) -> PyResult<Self> {
        Ok(PySmilesWriter {
            sink: Some(WriterSink::open(fileObj)?),
            writer: SmilesWriter::new(
                delimiter,
                nameHeader,
                includeHeader,
                isomericSmiles,
                kekuleSmiles,
            ),
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
    ///     - confId: (optional) ignored
    #[pyo3(signature = (mol, confId = -1))]
    #[allow(non_snake_case)]
    fn write(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>, confId: i32) -> PyResult<()> {
        let _ = confId;
        let sink = self.sink.as_mut().ok_or_else(|| no_stream(WRITE))?;
        let text = {
            let mut m = mol.borrow_mut();
            self.writer.format_mol(&mut m.mol).py()?
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
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySmilesWriter>()
}
