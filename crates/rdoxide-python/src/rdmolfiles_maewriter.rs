//! `MaeWriter` (upstream `Wrap/MaeWriter.cpp`): writes to a named file or a
//! text-mode Python file-like object.

use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::mae_writer::{MaeWriter, mol_to_mae_block};

use crate::errors::IntoPyResult;
use crate::errors::{Check, Site};
use crate::mol::PyMol;
use crate::writer_sink::{WriterSink, close_sink, no_stream};

const FILE: &str = "Code/GraphMol/FileParsers/MaeWriter.cpp";
const FLUSH: Site = Site {
    check: Check::Precondition,
    file: FILE,
    line: 565,
    expr: "dp_ostream",
};
const WRITE: Site = Site {
    check: Check::Precondition,
    file: FILE,
    line: 589,
    expr: "dp_ostream",
};

/// An experimental class for writing molecules to Maestro files.
///
///   Usage examples:
///
///     1) writing to a named file:
///
///        >>> writer = MaeWriter('out.mae')
///        >>> for mol in list_of_mols:
///        ...    writer.write(mol)
///
///     2) writing to a file-like object:
///
///        >>> import gzip
///        >>> outf=gzip.open('out.mae.gz','wt+')
///        >>> writer = MaeWriter(outf)
///        >>> for mol in list_of_mols:
///        ...   writer.write(mol)
///        >>> writer.close()
///        >>> outf.close()
///
///   By default all non-private molecule, atom and bond properties are written
///   to the Maestro file. This can be changed using the SetProps method:
///
///        >>> writer = MaeWriter('out.mae')
///        >>> writer.SetProps(['prop1','prop2'])
///
///   Properties that are specified, but are not present will be ignored.
///
///   Kekulization is mandatory, as the Maestro format does not have
///   the concept of an aromatic bond
///
///   As this is an experimental writer, many features are not supported yet,
///   e.g. chirality and bond stereo labels, stereo groups, substance groups,
///   isotopes, or even dummy atoms. Note that these aren't supported by
///   MaeMolSupplier either.
#[pyclass(name = "MaeWriter", module = "rdoxide.Chem.rdmolfiles", unsendable)]
pub struct PyMaeWriter {
    sink: Option<WriterSink>,
    writer: MaeWriter,
}

#[pymethods]
impl PyMaeWriter {
    #[new]
    #[pyo3(signature = (fileobj))]
    fn new(fileobj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(PyMaeWriter {
            sink: Some(WriterSink::open(fileobj)?),
            writer: MaeWriter::new(),
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

    /// Sets the atom and mol properties to be written to the output file
    ///
    ///   ARGUMENTS:
    ///
    ///     - props: a list or tuple of atom and mol property names
    #[pyo3(name = "SetProps", signature = (props_list))]
    fn set_props(&mut self, props_list: Vec<String>) {
        self.writer.set_props(props_list);
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
        let mol = mol.borrow().mol.clone();
        let writer = &mut self.writer;
        let text = py.detach(|| writer.write(&mol, confId)).py()?;
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

    /// returns the Maestro ct block text for a molecule
    #[staticmethod]
    #[pyo3(name = "GetText", signature = (mol, confId = -1, props_list = Vec::new()))]
    #[allow(non_snake_case)]
    fn get_text(
        py: Python<'_>,
        mol: &Bound<'_, PyMol>,
        confId: i32,
        props_list: Vec<String>,
    ) -> PyResult<String> {
        let mol = mol.borrow().mol.clone();
        py.detach(|| mol_to_mae_block(&mol, confId, &props_list))
            .py()
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMaeWriter>()
}
