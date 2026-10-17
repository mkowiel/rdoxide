//! `MultithreadedSDMolSupplier` and `MultithreadedSmilesMolSupplier`
//! (upstream `Wrap/Multithreaded{SD,Smiles}MolSupplier.cpp`).

use pyo3::exceptions::PyStopIteration;
use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::mol_file_parser::MolFileParserParams;
use rdoxide::graphmol::file_parsers::multithreaded_supplier::{
    MultithreadedMolSupplier, SupplierParams,
};
use rdoxide::graphmol::file_parsers::smiles_mol_supplier::SmilesMolSupplierParams;

use crate::mol::PyMol;
use crate::rdmolfiles_smisupplier::translate;

/// The shared `__next__` of the two classes. Like upstream, a request after
/// the last record (possible once, for an empty input) yields `None`.
fn next_item(py: Python<'_>, supplier: &mut MultithreadedMolSupplier) -> PyResult<Py<PyAny>> {
    // (blocks while the next record is being parsed: let other threads run)
    let at_end = py.detach(|| supplier.at_end());
    if at_end {
        return Err(PyStopIteration::new_err("End of supplier hit"));
    }
    match py.detach(|| supplier.next_mol()).flatten() {
        Some(mol) => Ok(Py::new(py, PyMol { mol })?.into_any()),
        None => Ok(py.None()),
    }
}

/// A class which concurrently supplies molecules from a text file.
///   Please note that this class is still a bit experimental and the API may
///   change in future releases.
///
///   Usage examples:
///
///     1) Lazy evaluation: the molecules might not be constructed until we ask for them:
///
///        >>> suppl = MultithreadedSDMolSupplier('in.sdf')
///        >>> for mol in suppl:
///        ...    if(mol):
///        ...      mol.GetNumAtoms()
///
///     2) Lazy evaluation 2:
///
///        >>> suppl = MultithreadedSDMolSupplier('in.sdf')
///        >>> while (!suppl.atEnd()):
///        ...    mol = next(mol)
///        ...    if(mol):
///        ...      mol.GetNumAtoms()
///
///   The molecules come back in the order of the file.
#[pyclass(
    name = "MultithreadedSDMolSupplier",
    module = "rdoxide.Chem.rdmolfiles",
    unsendable
)]
pub struct PyMultithreadedSDMolSupplier {
    supplier: MultithreadedMolSupplier,
}

#[pymethods]
impl PyMultithreadedSDMolSupplier {
    /// Constructor
    ///
    /// ARGUMENTS:
    ///
    ///   - fileName: name of the file to be read
    ///
    ///   - sanitize: (optional) toggles sanitization of molecules as they are read.
    ///     Defaults to true.
    ///
    ///   - removeHs: (optional) removes Hs. Defaults to true.
    ///
    ///   - strictParsing: (optional) allows strict or lax parsing. Defaults to true.
    ///
    ///   - numWriterThreads: (optional) number of writer threads. Defaults to 1.
    ///
    ///   - sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.
    ///
    ///   - sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.
    #[new]
    #[pyo3(signature = (fileName, sanitize = true, removeHs = true, strictParsing = true, numWriterThreads = 1, sizeInputQueue = 5, sizeOutputQueue = 5))]
    #[allow(non_snake_case)]
    fn new(
        fileName: &str,
        #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
        #[pyo3(from_py_with = crate::pyconv::truthy)] removeHs: bool,
        #[pyo3(from_py_with = crate::pyconv::truthy)] strictParsing: bool,
        numWriterThreads: i32,
        sizeInputQueue: usize,
        sizeOutputQueue: usize,
    ) -> PyResult<Self> {
        let parse_params = MolFileParserParams {
            sanitize,
            remove_hs: removeHs,
            strict_parsing: strictParsing,
            ..MolFileParserParams::default()
        };
        let params = SupplierParams {
            num_writer_threads: numWriterThreads,
            size_input_queue: sizeInputQueue,
            size_output_queue: sizeOutputQueue,
        };
        MultithreadedMolSupplier::from_sd_file(fileName, parse_params, params)
            .map(|supplier| PyMultithreadedSDMolSupplier { supplier })
            .map_err(translate)
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        &mut self,
        exc_type: &Bound<'_, PyAny>,
        exc_val: &Bound<'_, PyAny>,
        traceback: &Bound<'_, PyAny>,
    ) -> bool {
        let _ = (exc_type, exc_val, traceback);
        self.supplier.close();
        false
    }

    /// Returns the next molecule in the file. Raises _StopIteration_ on EOF.
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        next_item(py, &mut self.supplier)
    }

    /// Returns true if we have read all records else false.
    #[pyo3(name = "atEnd")]
    fn at_end(&mut self, py: Python<'_>) -> bool {
        py.detach(|| self.supplier.at_end())
    }

    /// Returns the record id for the last extracted item.
    #[pyo3(name = "GetLastRecordId")]
    fn get_last_record_id(&self) -> u32 {
        self.supplier.last_record_id()
    }

    /// Returns the text for the last extracted item.
    #[pyo3(name = "GetLastItemText")]
    fn get_last_item_text(&self) -> String {
        self.supplier.last_item_text()
    }

    /// returns whether or not any property lists that are present will be processed when reading molecules
    #[pyo3(name = "GetProcessPropertyLists")]
    fn get_process_property_lists(&self) -> bool {
        self.supplier.process_property_lists()
    }

    /// sets whether or not any property lists that are present will be processed when reading molecules
    #[pyo3(name = "SetProcessPropertyLists", signature = (val))]
    fn set_process_property_lists(&self, val: bool) {
        self.supplier.set_process_property_lists(val);
    }
}

/// A class which concurrently supplies molecules from a text file.
///   Please note that this class is still a bit experimental and the API may
///   change in future releases.
///
///   Usage examples:
///
///     1) Lazy evaluation: the molecules might not be constructed until we ask for them:
///
///        >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
///        >>> for mol in suppl:
///        ...    if(mol):
///        ...      mol.GetNumAtoms()
///
///     2) Lazy evaluation 2:
///
///        >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
///        >>> while (!suppl.atEnd()):
///        ...    mol = next(mol)
///        ...    if(mol):
///        ...      mol.GetNumAtoms()
///
///   The molecules come back in the order of the file.
#[pyclass(
    name = "MultithreadedSmilesMolSupplier",
    module = "rdoxide.Chem.rdmolfiles",
    unsendable
)]
pub struct PyMultithreadedSmilesMolSupplier {
    supplier: MultithreadedMolSupplier,
}

#[pymethods]
impl PyMultithreadedSmilesMolSupplier {
    /// Constructor
    ///
    /// ARGUMENTS:
    ///
    ///   - fileName: name of the file to be read
    ///
    ///   - delimiter: (optional) text delimiter (a string).  Defauts to ' \t'.
    ///
    ///   - smilesColumn: (optional) index of the column containing the SMILES
    ///     data.  Defaults to 0.
    ///
    ///   - nameColumn: (optional) index of the column containing molecule names.
    ///     Defaults to 1.
    ///
    ///   - titleLine: (optional) set this toggle if the file contains a title line.
    ///     Defaults to true.
    ///
    ///   - sanitize: (optional) toggles sanitization of molecules as they are read.
    ///     Defaults to true.
    ///
    ///   - numWriterThreads: (optional) number of writer threads. Defaults to 1.
    ///
    ///   - sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.
    ///
    ///   - sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.
    #[new]
    #[pyo3(signature = (fileName, delimiter = " \t", smilesColumn = 0, nameColumn = 1, titleLine = true, sanitize = true, numWriterThreads = 1, sizeInputQueue = 5, sizeOutputQueue = 5))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn new(
        fileName: &str,
        delimiter: &str,
        smilesColumn: i32,
        nameColumn: i32,
        #[pyo3(from_py_with = crate::pyconv::truthy)] titleLine: bool,
        #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
        numWriterThreads: i32,
        sizeInputQueue: usize,
        sizeOutputQueue: usize,
    ) -> PyResult<Self> {
        let mut parse_params = SmilesMolSupplierParams {
            delimiter: delimiter.to_string(),
            smiles_column: smilesColumn,
            name_column: nameColumn,
            title_line: titleLine,
            ..SmilesMolSupplierParams::default()
        };
        parse_params.parse_parameters.sanitize = sanitize;
        let params = SupplierParams {
            num_writer_threads: numWriterThreads,
            size_input_queue: sizeInputQueue,
            size_output_queue: sizeOutputQueue,
        };
        MultithreadedMolSupplier::from_smiles_file(fileName, parse_params, params)
            .map(|supplier| PyMultithreadedSmilesMolSupplier { supplier })
            .map_err(translate)
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        &mut self,
        exc_type: &Bound<'_, PyAny>,
        exc_val: &Bound<'_, PyAny>,
        traceback: &Bound<'_, PyAny>,
    ) -> bool {
        let _ = (exc_type, exc_val, traceback);
        self.supplier.close();
        false
    }

    /// Returns the next molecule in the file. Raises _StopIteration_ on EOF.
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        next_item(py, &mut self.supplier)
    }

    /// Returns true if we have read all records else false.
    #[pyo3(name = "atEnd")]
    fn at_end(&mut self, py: Python<'_>) -> bool {
        py.detach(|| self.supplier.at_end())
    }

    /// Returns the record id for the last extracted item.
    #[pyo3(name = "GetLastRecordId")]
    fn get_last_record_id(&self) -> u32 {
        self.supplier.last_record_id()
    }

    /// Returns the text for the last extracted item.
    #[pyo3(name = "GetLastItemText")]
    fn get_last_item_text(&self) -> String {
        self.supplier.last_item_text()
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMultithreadedSDMolSupplier>()?;
    m.add_class::<PyMultithreadedSmilesMolSupplier>()
}
