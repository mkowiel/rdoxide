//! `TDTMolSupplier` (upstream `Wrap/TDTMolSupplier.cpp`).

use pyo3::exceptions::{PyIndexError, PyStopIteration};
use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::tdt_mol_supplier::{TdtMolSupplier, TdtMolSupplierParams};

use crate::mol::PyMol;
use crate::rdmolfiles_smisupplier::translate;

fn params(
    name_record: &str,
    conf_id_2d: i32,
    conf_id_3d: i32,
    sanitize: bool,
) -> TdtMolSupplierParams {
    let mut p = TdtMolSupplierParams {
        name_record: name_record.to_string(),
        conf_id_2d,
        conf_id_3d,
        ..TdtMolSupplierParams::default()
    };
    p.parse_parameters.sanitize = sanitize;
    p
}

/// A class which supplies molecules from a TDT file.
///
///   Usage examples:
///
///     1) Lazy evaluation: the molecules are not constructed until we ask for them:
///
///        >>> suppl = TDTMolSupplier('in.smi')
///        >>> for mol in suppl:
///        ...    mol.GetNumAtoms()
///
///     2) Lazy evaluation 2:
///
///        >>> suppl = TDTMolSupplier('in.smi')
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
///        >>> suppl = TDTMolSupplier('in.smi')
///        >>> nMols = len(suppl)
///        >>> for i in range(nMols):
///        ...   suppl[i].GetNumAtoms()
///
///   Properties in the file are used to set properties on each molecule.
///   The properties are accessible using the mol.GetProp(propName) method.
#[pyclass(
    name = "TDTMolSupplier",
    module = "rdoxide.Chem.rdmolfiles",
    unsendable
)]
pub struct PyTdtMolSupplier {
    inner: TdtMolSupplier,
}

#[pymethods]
impl PyTdtMolSupplier {
    #[new]
    #[pyo3(signature = (fileName = None, nameRecord = "", confId2D = -1, confId3D = -1, sanitize = true))]
    #[allow(non_snake_case)]
    fn new(
        fileName: Option<&str>,
        nameRecord: &str,
        confId2D: i32,
        confId3D: i32,
        sanitize: bool,
    ) -> PyResult<Self> {
        let inner = match fileName {
            None => TdtMolSupplier::new(),
            Some(f) => {
                TdtMolSupplier::from_file(f, params(nameRecord, confId2D, confId3D, sanitize))
                    .map_err(translate)?
            }
        };
        Ok(PyTdtMolSupplier { inner })
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
    #[pyo3(name = "SetData", signature = (data, nameRecord = "", confId2D = -1, confId3D = -1, sanitize = true))]
    #[allow(non_snake_case)]
    fn set_data(
        &mut self,
        data: &str,
        nameRecord: &str,
        confId2D: i32,
        confId3D: i32,
        sanitize: bool,
    ) {
        self.inner.set_data(
            data.as_bytes(),
            params(nameRecord, confId2D, confId3D, sanitize),
        );
    }

    /// returns the text for an item
    #[pyo3(name = "GetItemText", signature = (index))]
    fn get_item_text(&mut self, index: u32) -> PyResult<String> {
        self.inner.get_item_text(index as usize).map_err(translate)
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTdtMolSupplier>()
}
