//! `SDWriter` (upstream `Wrap/SDWriter.cpp`): writes to a named file,
//! stdout (`"-"`) or a Python file-like object.

use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::sd_writer::mol_to_sd_text;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::writer_sink::{Check, Site, WriterSink, close_sink, no_stream};

const WRITE: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/FileParsers/SDWriter.cpp",
    line: 171,
};
const FLUSH: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/FileParsers/MolWriters.h",
    line: 169,
};

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
    sink: Option<WriterSink>,
    props: Vec<String>,
    molid: u32,
    kekulize: bool,
    force_v3000: bool,
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
        let sink = WriterSink::open(fileName)?;
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
            return Err(no_stream(WRITE));
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
        self.sink
            .as_mut()
            .ok_or_else(|| no_stream(WRITE))?
            .write(py, &text)?;
        self.molid += 1;
        Ok(())
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

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySDWriter>()
}
