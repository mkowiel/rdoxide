//! `MaeMolSupplier` (upstream `Wrap/MaeMolSupplier.cpp` and `Wrap/MolSupplier.h`).

use pyo3::exceptions::{PyIndexError, PyStopIteration, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::mae_mol_supplier::{MaeMolSupplier, MaeMolSupplierParams};

use crate::errors::to_py;
use crate::mol::PyMol;
use crate::rdmolfiles_sdsupplier::translate;

const SOURCE: &str = "Code/GraphMol/FileParsers/MaeMolSupplier.cpp";

/// A class which supplies molecules from file-like object containing Maestro data.
///
///   Usage examples:
///
///     1) Lazy evaluation: the molecules are not constructed until we ask for them:
///
///        >>> suppl = MaeMolSupplier(file('in.mae'))
///        >>> for mol in suppl:
///        ...    if mol is not None: mol.GetNumAtoms()
///
///     2) we can also read from compressed files:
///
///        >>> import gzip
///        >>> suppl = MaeMolSupplier(gzip.open('in.maegz'))
///        >>> for mol in suppl:
///        ...   if mol is not None: print mol.GetNumAtoms()
///
///   Properties in the Maestro file are used to set properties on each molecule.
///   The properties are accessible using the mol.GetProp(propName) method.
#[pyclass(
    name = "MaeMolSupplier",
    module = "rdoxide.Chem.rdmolfiles",
    unsendable
)]
pub struct PyMaeMolSupplier {
    inner: Option<MaeMolSupplier>,
}

fn params(sanitize: bool, remove_hs: bool) -> MaeMolSupplierParams {
    MaeMolSupplierParams {
        sanitize,
        remove_hs,
    }
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

impl PyMaeMolSupplier {
    fn supplier(&mut self) -> PyResult<&mut MaeMolSupplier> {
        self.inner.as_mut().ok_or_else(|| {
            to_py(RDKitError::precondition_at(
                "no stream",
                "dp_sInStream != nullptr",
                SOURCE,
                595,
            ))
        })
    }
}

#[pymethods]
impl PyMaeMolSupplier {
    #[new]
    #[pyo3(signature = (fileobj = None, sanitize = true, removeHs = true))]
    #[allow(non_snake_case)]
    fn new(
        py: Python<'_>,
        fileobj: Option<&Bound<'_, PyAny>>,
        #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
        #[pyo3(from_py_with = crate::pyconv::truthy)] removeHs: bool,
    ) -> PyResult<Self> {
        let params = params(sanitize, removeHs);
        let Some(fileobj) = fileobj else {
            return Ok(PyMaeMolSupplier { inner: None });
        };
        let inner = if let Ok(name) = fileobj.extract::<String>() {
            MaeMolSupplier::from_file(&name, params).map_err(translate)?
        } else {
            let text_base = py.import("io")?.getattr("TextIOBase")?;
            if fileobj.is_instance(&text_base)? {
                return Err(PyValueError::new_err(
                    "Need a binary mode file object like BytesIO or a file opened with mode 'b'",
                ));
            }
            MaeMolSupplier::from_bytes(read_all(fileobj)?, params).map_err(translate)?
        };
        Ok(PyMaeMolSupplier { inner: Some(inner) })
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
        self.inner = None;
        false
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let suppl = self.supplier()?;
        if suppl.at_end() {
            return Err(PyStopIteration::new_err("End of supplier hit"));
        }
        let mol = match py.detach(|| suppl.next_mol()) {
            Ok(m) => Some(PyMol { mol: m }),
            Err(e @ RDKitError::FileParse(_)) => return Err(translate(e)),
            Err(_) => None,
        };
        Ok(mol.into_pyobject(py)?.into_any().unbind())
    }

    fn __getitem__(&mut self, py: Python<'_>, idx: i64) -> PyResult<Option<PyMol>> {
        let mut idx = idx;
        let suppl = self.supplier()?;
        if idx < 0 {
            idx += i64::from(suppl.length().map_err(translate)?);
            if idx < 0 {
                return Err(PyIndexError::new_err("invalid index"));
            }
        }
        match py.detach(|| suppl.get(idx as u32)) {
            Ok(mol) => Ok(Some(PyMol { mol })),
            Err(_) => {
                if suppl.at_end() {
                    Err(PyIndexError::new_err("invalid index"))
                } else {
                    Ok(None)
                }
            }
        }
    }

    /// Resets our position in the file to the beginning.
    fn reset(&mut self) -> PyResult<()> {
        self.supplier()?.reset().map_err(translate)
    }

    fn __len__(&mut self, py: Python<'_>) -> PyResult<usize> {
        let suppl = self.supplier()?;
        py.detach(|| suppl.length())
            .map(|n| n as usize)
            .map_err(translate)
    }

    /// Sets the text to be parsed
    #[pyo3(name = "SetData", signature = (data, sanitize = true, removeHs = true))]
    #[allow(non_snake_case)]
    fn set_data(
        &mut self,
        #[pyo3(from_py_with = crate::pyconv::text_bytes)] data: Vec<u8>,
        #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
        #[pyo3(from_py_with = crate::pyconv::truthy)] removeHs: bool,
    ) -> PyResult<()> {
        self.inner =
            Some(MaeMolSupplier::from_bytes(data, params(sanitize, removeHs)).map_err(translate)?);
        Ok(())
    }

    /// Returns whether or not we have hit EOF.
    #[pyo3(name = "atEnd")]
    fn at_end(&mut self) -> PyResult<bool> {
        Ok(self.supplier()?.at_end())
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMaeMolSupplier>()
}
