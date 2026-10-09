//! `SDMolSupplier` (upstream `Wrap/SDMolSupplier.cpp` and `Wrap/MolSupplier.h`).

use pyo3::exceptions::{PyIOError, PyIndexError, PyRuntimeError, PyStopIteration};
use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::file_parsers::sd_mol_supplier::SDMolSupplier;

use crate::errors::to_py;
use crate::mol::PyMol;

/// The translators registered in rdmolfiles.cpp.
fn translate(e: RDKitError) -> PyErr {
    match e {
        RDKitError::BadFile(m) => PyIOError::new_err(format!("File error: {m}")),
        RDKitError::FileParse(m) => PyRuntimeError::new_err(format!("File parsing error: {m}")),
        e => to_py(e),
    }
}

fn params(sanitize: bool, remove_hs: bool, strict_parsing: bool) -> MolFileParserParams {
    MolFileParserParams {
        sanitize,
        remove_hs,
        strict_parsing,
        ..MolFileParserParams::default()
    }
}

/// A class which supplies molecules from an SD file.
///
///   Usage examples:
///
///     1) Lazy evaluation: the molecules are not constructed until we ask for them:
///
///        >>> suppl = SDMolSupplier('in.sdf')
///        >>> for mol in suppl:
///        ...    mol.GetNumAtoms()
///
///     2) Lazy evaluation 2:
///
///        >>> suppl = SDMolSupplier('in.sdf')
///        >>> mol1 = next(suppl)
///        >>> mol2 = next(suppl)
///        >>> suppl.reset()
///        >>> mol3 = next(suppl)
///        # mol3 and mol1 are the same:
///        >>> MolToSmiles(mol3)==MolToSmiles(mol1)
///
///     3) Random Access:
///
///        >>> suppl = SDMolSupplier('in.sdf')
///        >>> mol1 = suppl[0]
///        >>> mol2 = suppl[1]
///        # NOTE: this will generate an IndexError if the supplier doesn't have that many
///        molecules.
///
///     4) Random Access 2:  looping over all molecules
///
///        >>> suppl = SDMolSupplier('in.sdf')
///        >>> nMols = len(suppl)
///        >>> for i in range(nMols):
///        ...   suppl[i].GetNumAtoms()
///
///   Properties in the SD file are used to set properties on each molecule.
///   The properties are accessible using the mol.GetProp(propName) method.
#[pyclass(
    name = "SDMolSupplier",
    module = "rdkit.Chem.rdmolfiles",
    subclass,
    unsendable
)]
pub struct PySDMolSupplier {
    inner: SDMolSupplier,
}

#[pymethods]
impl PySDMolSupplier {
    #[new]
    #[pyo3(signature = (fileName = None, sanitize = true, removeHs = true, strictParsing = true))]
    #[allow(non_snake_case)]
    fn new(
        fileName: Option<String>,
        sanitize: bool,
        removeHs: bool,
        strictParsing: bool,
    ) -> PyResult<Self> {
        let inner = match fileName {
            None => SDMolSupplier::default(),
            Some(f) => SDMolSupplier::from_file(&f, params(sanitize, removeHs, strictParsing))
                .map_err(translate)?,
        };
        Ok(PySDMolSupplier { inner })
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
    #[pyo3(name = "SetData", signature = (data, sanitize = true, removeHs = true, strictParsing = true))]
    #[allow(non_snake_case)]
    fn set_data(&mut self, data: &str, sanitize: bool, removeHs: bool, strictParsing: bool) {
        self.inner
            .set_data(data.as_bytes(), params(sanitize, removeHs, strictParsing));
    }

    /// Sets the locations of mol beginnings in the input stream. Be *very* careful with this method.
    #[pyo3(name = "_SetStreamIndices", signature = (locs))]
    fn set_stream_indices(&mut self, locs: Vec<i64>) -> PyResult<()> {
        let locs = locs.into_iter().map(|l| l.max(0) as usize).collect();
        self.inner.set_stream_indices(locs).map_err(translate)
    }

    /// returns the text for an item
    #[pyo3(name = "GetItemText", signature = (index))]
    fn get_item_text(&mut self, index: u32) -> PyResult<String> {
        let t = self
            .inner
            .get_item_text(index as usize)
            .map_err(translate)?;
        Ok(String::from_utf8_lossy(&t).into_owned())
    }

    /// Returns whether or not we have hit EOF.
    #[pyo3(name = "atEnd")]
    fn at_end(&self) -> PyResult<bool> {
        self.inner.at_end().map_err(translate)
    }

    /// returns whether or not any property lists that are present will be processed when reading molecules
    #[pyo3(name = "GetProcessPropertyLists")]
    fn get_process_property_lists(&self) -> bool {
        self.inner.process_property_lists
    }

    /// sets whether or not any property lists that are present will be processed when reading molecules
    #[pyo3(name = "SetProcessPropertyLists", signature = (val))]
    fn set_process_property_lists(&mut self, val: bool) {
        self.inner.process_property_lists = val;
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySDMolSupplier>()?;
    Ok(())
}
