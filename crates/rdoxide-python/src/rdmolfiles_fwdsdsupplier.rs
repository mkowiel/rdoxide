//! `ForwardSDMolSupplier` (upstream `Wrap/ForwardSDMolSupplier.cpp` and
//! `Wrap/MolSupplier.h`): SD records from a binary file object or a file name.

use pyo3::exceptions::{PyIOError, PyStopIteration, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::file_parsers::sd_mol_supplier::ForwardSDMolSupplier;

use crate::errors::{Check, Site, check_failed, to_py};
use crate::mol::PyMol;

const NO_STREAM: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp",
    line: 332,
    expr: "dp_inStream",
};

/// Where the records come from: read once, on first use.
enum Source {
    /// A binary file-like object that is read to its end when first needed.
    Pending(Py<PyAny>),
    Open(ForwardSDMolSupplier),
    Closed,
}

/// A class which supplies molecules from file-like object containing SD data.
///
///   Usage examples:
///
///     1) Lazy evaluation: the molecules are not constructed until we ask for them:
///
///        >>> suppl = ForwardSDMolSupplier(file('in.sdf'))
///        >>> for mol in suppl:
///        ...    if mol is not None: mol.GetNumAtoms()
///
///     2) we can also read from compressed files:
///
///        >>> import gzip
///        >>> suppl = ForwardSDMolSupplier(gzip.open('in.sdf.gz'))
///        >>> for mol in suppl:
///        ...   if mol is not None: print mol.GetNumAtoms()
///
///   Properties in the SD file are used to set properties on each molecule.
///   The properties are accessible using the mol.GetProp(propName) method.
#[pyclass(
    name = "ForwardSDMolSupplier",
    module = "rdoxide.Chem.rdmolfiles",
    unsendable
)]
pub struct PyForwardSDMolSupplier {
    source: Source,
    params: MolFileParserParams,
    process_property_lists: bool,
}

/// Everything a binary file object yields (an object without `read` is an
/// empty stream, as with upstream's stream buffer).
fn read_all(file: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    if !file.hasattr("read")? {
        return Ok(Vec::new());
    }
    let data = file.call_method0("read")?;
    Ok(data.cast::<PyBytes>()?.as_bytes().to_vec())
}

impl PyForwardSDMolSupplier {
    /// The supplier, reading a pending file object on first use.
    fn supplier(&mut self, py: Python<'_>) -> PyResult<&mut ForwardSDMolSupplier> {
        if let Source::Pending(file) = &self.source {
            let data = read_all(file.bind(py))?;
            let mut suppl = ForwardSDMolSupplier::from_bytes(data, self.params.clone());
            suppl.set_process_property_lists(self.process_property_lists);
            self.source = Source::Open(suppl);
        }
        match &mut self.source {
            Source::Open(s) => Ok(s),
            _ => Err(check_failed(NO_STREAM, "no stream")),
        }
    }
}

#[pymethods]
impl PyForwardSDMolSupplier {
    #[new]
    #[pyo3(signature = (fileobj, sanitize = true, removeHs = true, strictParsing = true))]
    #[allow(non_snake_case)]
    fn new(
        py: Python<'_>,
        fileobj: &Bound<'_, PyAny>,
        #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
        #[pyo3(from_py_with = crate::pyconv::truthy)] removeHs: bool,
        #[pyo3(from_py_with = crate::pyconv::truthy)] strictParsing: bool,
    ) -> PyResult<Self> {
        let params = MolFileParserParams {
            sanitize,
            remove_hs: removeHs,
            strict_parsing: strictParsing,
            ..MolFileParserParams::default()
        };
        let source = if let Ok(name) = fileobj.extract::<String>() {
            Source::Open(
                ForwardSDMolSupplier::from_file(&name, params.clone()).map_err(|e| match e {
                    RDKitError::BadFile(m) => PyIOError::new_err(format!("File error: {m}")),
                    e => to_py(e),
                })?,
            )
        } else {
            let text_base = py.import("io")?.getattr("TextIOBase")?;
            if fileobj.is_instance(&text_base)? {
                return Err(PyValueError::new_err(
                    "Need a binary mode file object like BytesIO or a file opened with mode 'b'",
                ));
            }
            Source::Pending(fileobj.clone().unbind())
        };
        Ok(PyForwardSDMolSupplier {
            source,
            params,
            process_property_lists: true,
        })
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> {
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
        self.source = Source::Closed;
        false
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
    /// (returns a Python object: PyO3 would turn `Option::None` into
    /// StopIteration, but upstream yields `None` for unparsable records)
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let stop = || PyStopIteration::new_err("End of supplier hit");
        let suppl = self.supplier(py)?;
        if suppl.at_end().map_err(to_py)? {
            return Err(stop());
        }
        let mol = match py.detach(|| suppl.next_mol()) {
            Ok(m) => m.map(|mol| PyMol { mol }),
            Err(e @ RDKitError::FileParse(_)) => return Err(to_py(e)),
            Err(_) => None,
        };
        if suppl.at_end().map_err(to_py)? && suppl.eof_hit_on_read() {
            return Err(stop());
        }
        Ok(mol.into_pyobject(py)?.into_any().unbind())
    }

    /// Returns whether or not we have hit EOF.
    #[pyo3(name = "atEnd")]
    fn at_end(&mut self, py: Python<'_>) -> PyResult<bool> {
        self.supplier(py)?.at_end().map_err(to_py)
    }

    /// Returns whether or EOF was hit while parsing the previous entry.
    #[pyo3(name = "GetEOFHitOnRead")]
    fn get_eof_hit_on_read(&mut self, py: Python<'_>) -> PyResult<bool> {
        Ok(self.supplier(py)?.eof_hit_on_read())
    }

    /// returns whether or not any property lists that are present will be processed when reading molecules
    #[pyo3(name = "GetProcessPropertyLists")]
    fn get_process_property_lists(&self) -> bool {
        match &self.source {
            Source::Open(s) => s.process_property_lists(),
            _ => self.process_property_lists,
        }
    }

    /// sets whether or not any property lists that are present will be processed when reading molecules
    #[pyo3(name = "SetProcessPropertyLists", signature = (val))]
    fn set_process_property_lists(&mut self, val: bool) {
        self.process_property_lists = val;
        if let Source::Open(s) = &mut self.source {
            s.set_process_property_lists(val);
        }
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyForwardSDMolSupplier>()
}
