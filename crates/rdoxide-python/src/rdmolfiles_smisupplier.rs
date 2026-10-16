//! `SmilesMolSupplier` (upstream `Wrap/SmilesMolSupplier.cpp` and
//! `Wrap/MolSupplier.h`).

use pyo3::exceptions::{PyIOError, PyIndexError, PyRuntimeError, PyStopIteration};
use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::smiles_mol_supplier::{
    SmilesMolSupplier, SmilesMolSupplierParams,
};

use crate::errors::to_py;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// The translators registered in rdmolfiles.cpp.
pub(crate) fn translate(e: RDKitError) -> PyErr {
    match e {
        RDKitError::BadFile(m) => PyIOError::new_err(format!("File error: {m}")),
        RDKitError::FileParse(m) => PyRuntimeError::new_err(format!("File parsing error: {m}")),
        e => to_py(e),
    }
}

/// The wrapper's arguments as `SmilesMolSupplierParams` (the v1 API only
/// carries `sanitize` into the SMILES parser parameters).
fn params(
    delimiter: &str,
    smiles_column: i32,
    name_column: i32,
    title_line: Option<&Bound<'_, PyAny>>,
    sanitize: Option<&Bound<'_, PyAny>>,
) -> PyResult<SmilesMolSupplierParams> {
    let mut p = SmilesMolSupplierParams {
        delimiter: delimiter.to_string(),
        smiles_column,
        name_column,
        title_line: bool_arg(title_line, true)?,
        ..SmilesMolSupplierParams::default()
    };
    p.parse_parameters.sanitize = bool_arg(sanitize, true)?;
    Ok(p)
}

/// A class which supplies molecules from a text file.
///
///   Usage examples:
///
///     1) Lazy evaluation: the molecules are not constructed until we ask for them:
///
///        >>> suppl = SmilesMolSupplier('in.smi')
///        >>> for mol in suppl:
///        ...    mol.GetNumAtoms()
///
///     2) Lazy evaluation 2:
///
///        >>> suppl = SmilesMolSupplier('in.smi')
///        >>> mol1 = next(suppl)
///        >>> mol2 = next(suppl)
///        >>> suppl.reset()
///        >>> mol3 = next(suppl)
///        # mol3 and mol1 are the same:
///        >>> MolToSmiles(mol3)==MolToSmiles(mol1)
///
///     3) Random Access:  all molecules are constructed as soon as we ask for the
///        length:
///
///        >>> suppl = SmilesMolSupplier('in.smi')
///        >>> nMols = len(suppl)
///        >>> for i in range(nMols):
///        ...   suppl[i].GetNumAtoms()
///
///   If the input file has a title line and more than two columns (smiles and id), the
///   additional columns will be used to set properties on each molecule.  The properties
///   are accessible using the mol.GetProp(propName) method.
#[pyclass(
    name = "SmilesMolSupplier",
    module = "rdoxide.Chem.rdmolfiles",
    unsendable
)]
pub struct PySmilesMolSupplier {
    inner: SmilesMolSupplier,
}

#[pymethods]
impl PySmilesMolSupplier {
    /// Constructor
    ///
    ///   ARGUMENTS:
    ///
    ///     - fileName: name of the file to be read
    ///
    ///     - delimiter: (optional) text delimiter (a string).  Defauts to ' '.
    ///
    ///     - smilesColumn: (optional) index of the column containing the SMILES
    ///       data.  Defaults to 0.
    ///
    ///     - nameColumn: (optional) index of the column containing molecule names.
    ///       Defaults to 1.
    ///
    ///     - titleLine: (optional) set this toggle if the file contains a title line.
    ///       Defaults to 1.
    ///
    ///     - sanitize: (optional) toggles sanitization of molecules as they are read.
    ///       Defaults to 1.
    #[new]
    #[pyo3(signature = (data = None, delimiter = " ", smilesColumn = 0, nameColumn = 1, titleLine = None, sanitize = None))]
    #[allow(non_snake_case)]
    fn new(
        data: Option<String>,
        delimiter: &str,
        smilesColumn: i32,
        nameColumn: i32,
        titleLine: Option<&Bound<'_, PyAny>>,
        sanitize: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let inner = match data {
            None => SmilesMolSupplier::default(),
            Some(f) => {
                let p = params(delimiter, smilesColumn, nameColumn, titleLine, sanitize)?;
                SmilesMolSupplier::from_file(&f, p).map_err(translate)?
            }
        };
        Ok(PySmilesMolSupplier { inner })
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
        self.inner.close();
        false
    }

    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyRefMut<'_, Self>> {
        slf.inner.reset().map_err(translate)?;
        Ok(slf)
    }

    /// Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
    /// (returns a Python object: PyO3 would turn `Option::None` into
    /// StopIteration, but upstream yields `None` for unparsable records)
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        if self.inner.at_end().map_err(translate)? {
            return Err(PyStopIteration::new_err("End of supplier hit"));
        }
        let inner = &mut self.inner;
        let m = match py.detach(|| inner.next_mol()) {
            Ok(m) => m.map(|mol| PyMol { mol }),
            Err(e @ RDKitError::FileParse(_)) => return Err(translate(e)),
            Err(_) => None,
        };
        Ok(m.into_pyobject(py)?.into_any().unbind())
    }

    fn __getitem__(&mut self, py: Python<'_>, idx: i64) -> PyResult<Option<PyMol>> {
        let mut idx = idx;
        if idx < 0 {
            idx += self.inner.length().map_err(translate)? as i64;
            if idx < 0 {
                return Err(PyIndexError::new_err("invalid index"));
            }
        }
        let inner = &mut self.inner;
        match py.detach(|| inner.get(idx as usize)) {
            Ok(m) => Ok(m.map(|mol| PyMol { mol })),
            Err(_) => {
                if self.inner.at_end().unwrap_or(false) {
                    Err(PyIndexError::new_err("invalid index"))
                } else {
                    Ok(None)
                }
            }
        }
    }

    /// Resets our position in the file to the beginning.
    fn reset(&mut self) -> PyResult<()> {
        self.inner.reset().map_err(translate)
    }

    fn __len__(&mut self, py: Python<'_>) -> PyResult<usize> {
        let inner = &mut self.inner;
        py.detach(|| inner.length()).map_err(translate)
    }

    /// Sets the text to be parsed
    #[pyo3(name = "SetData", signature = (data, delimiter = " ", smilesColumn = 0, nameColumn = 1, titleLine = None, sanitize = None))]
    #[allow(non_snake_case)]
    fn set_data(
        &mut self,
        data: &str,
        delimiter: &str,
        smilesColumn: i32,
        nameColumn: i32,
        titleLine: Option<&Bound<'_, PyAny>>,
        sanitize: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let p = params(delimiter, smilesColumn, nameColumn, titleLine, sanitize)?;
        self.inner.set_data(data.as_bytes(), p);
        Ok(())
    }

    /// returns the text for an item
    #[pyo3(name = "GetItemText", signature = (index))]
    fn get_item_text(&mut self, index: u32) -> PyResult<String> {
        self.inner.get_item_text(index as usize).map_err(translate)
    }
}

/// Creates a supplier over the SMILES text in `text`
#[pyfunction]
#[pyo3(name = "SmilesMolSupplierFromText", signature = (text, delimiter = " ", smilesColumn = 0, nameColumn = 1, titleLine = None, sanitize = None))]
#[allow(non_snake_case)]
fn smiles_supplier_from_text(
    text: &str,
    delimiter: &str,
    smilesColumn: i32,
    nameColumn: i32,
    titleLine: Option<&Bound<'_, PyAny>>,
    sanitize: Option<&Bound<'_, PyAny>>,
) -> PyResult<PySmilesMolSupplier> {
    let p = params(delimiter, smilesColumn, nameColumn, titleLine, sanitize)?;
    let mut inner = SmilesMolSupplier::default();
    inner.set_data(text.as_bytes(), p);
    Ok(PySmilesMolSupplier { inner })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySmilesMolSupplier>()?;
    m.add_function(wrap_pyfunction!(smiles_supplier_from_text, m)?)?;
    Ok(())
}
