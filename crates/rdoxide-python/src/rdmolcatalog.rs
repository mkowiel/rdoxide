//! `rdoxide.Chem.rdMolCatalog` (upstream `MolCatalog/Wrap/rdMolCatalog.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple};
use rdoxide::graphmol::mol_catalog::{MolCatalog, MolCatalogEntry, MolCatalogParams};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// A molecule catalog entry.
#[pyclass(
    name = "MolCatalogEntry",
    module = "rdoxide.Chem.rdMolCatalog",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyMolCatalogEntry {
    entry: MolCatalogEntry,
}

#[pymethods]
impl PyMolCatalogEntry {
    #[new]
    #[pyo3(signature = (pickle = None))]
    fn py_new(pickle: Option<&[u8]>) -> PyResult<Self> {
        Ok(Self {
            entry: match pickle {
                Some(p) => MolCatalogEntry::from_pickle(p).py()?,
                None => MolCatalogEntry::default(),
            },
        })
    }

    #[pyo3(name = "Serialize")]
    fn serialize<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, &self.entry.serialize().py()?))
    }

    #[pyo3(name = "GetDescription")]
    fn get_description(&self) -> String {
        rdoxide::catalogs::CatalogEntry::description(&self.entry)
    }

    #[pyo3(name = "SetDescription")]
    fn set_description(&mut self, val: &str) {
        self.entry.set_description(val);
    }

    #[pyo3(name = "GetMol")]
    fn get_mol(&self) -> PyResult<PyMol> {
        let mol = self
            .entry
            .mol()
            .ok_or_else(|| PyValueError::new_err("the entry has no molecule"))?;
        Ok(PyMol { mol: mol.clone() })
    }

    #[pyo3(name = "SetMol")]
    fn set_mol(&mut self, mol: &Bound<'_, PyMol>) {
        self.entry.set_mol(mol.borrow().mol.clone());
    }

    #[pyo3(name = "GetOrder")]
    fn get_order(&self) -> u32 {
        rdoxide::catalogs::CatalogEntry::order(&self.entry) as u32
    }

    #[pyo3(name = "SetOrder")]
    fn set_order(&mut self, order: u32) {
        self.entry.set_order(order);
    }

    fn __reduce__<'py>(
        slf: &Bound<'py, Self>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyTuple>)> {
        let py = slf.py();
        let pkl = PyBytes::new(py, &slf.borrow().entry.serialize().py()?);
        Ok((py.get_type::<Self>().into_any(), PyTuple::new(py, [pkl])?))
    }
}

/// A hierarchical catalog of molecules.
#[pyclass(name = "MolCatalog", module = "rdoxide.Chem.rdMolCatalog")]
pub struct PyMolCatalog {
    cat: MolCatalog,
}

#[pymethods]
impl PyMolCatalog {
    #[new]
    fn py_new(pickle: &[u8]) -> PyResult<Self> {
        Ok(Self {
            cat: MolCatalog::from_pickle(pickle).py()?,
        })
    }

    #[pyo3(name = "GetNumEntries")]
    fn get_num_entries(&self) -> u32 {
        self.cat.num_entries()
    }

    #[pyo3(name = "GetFPLength")]
    fn get_fp_length(&self) -> u32 {
        self.cat.fp_length()
    }

    #[pyo3(name = "Serialize")]
    fn serialize<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, &self.cat.serialize().py()?))
    }

    #[pyo3(name = "GetBitDescription")]
    fn get_bit_description(&self, idx: u32) -> PyResult<String> {
        let entry = self.cat.entry_with_bit_id(idx).py()?;
        entry
            .map(rdoxide::catalogs::CatalogEntry::description)
            .ok_or_else(|| PyValueError::new_err("no entry has this bit id"))
    }

    #[pyo3(name = "GetBitEntryId")]
    fn get_bit_entry_id(&self, idx: u32) -> PyResult<i32> {
        self.cat.id_of_entry_with_bit_id(idx).py()
    }

    #[pyo3(name = "GetEntryBitId")]
    fn get_entry_bit_id(&self, idx: u32) -> PyResult<i32> {
        Ok(rdoxide::catalogs::CatalogEntry::bit_id(
            self.cat.entry_with_idx(idx).py()?,
        ))
    }

    #[pyo3(name = "GetEntryDescription")]
    fn get_entry_description(&self, idx: u32) -> PyResult<String> {
        Ok(rdoxide::catalogs::CatalogEntry::description(
            self.cat.entry_with_idx(idx).py()?,
        ))
    }

    #[pyo3(name = "GetEntryDownIds")]
    fn get_entry_down_ids(&self, idx: u32) -> PyResult<Vec<u32>> {
        Ok(self.cat.down_entry_list(idx).py()?.to_vec())
    }

    #[pyo3(name = "AddEntry")]
    fn add_entry(&mut self, entry: &Bound<'_, PyMolCatalogEntry>) -> u32 {
        self.cat.add_entry(entry.borrow().entry.clone(), true)
    }

    #[pyo3(name = "AddEdge")]
    fn add_edge(&mut self, id1: u32, id2: u32) -> PyResult<()> {
        self.cat.add_edge(id1, id2).py()
    }

    fn __reduce__<'py>(
        slf: &Bound<'py, Self>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyTuple>)> {
        let py = slf.py();
        let pkl = PyBytes::new(py, &slf.borrow().cat.serialize().py()?);
        Ok((py.get_type::<Self>().into_any(), PyTuple::new(py, [pkl])?))
    }
}

#[pyfunction]
#[pyo3(name = "CreateMolCatalog")]
fn create_mol_catalog() -> PyMolCatalog {
    PyMolCatalog {
        cat: MolCatalog::new(MolCatalogParams::default()),
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMolCatalog>()?;
    m.add_class::<PyMolCatalogEntry>()?;
    m.add_function(wrap_pyfunction!(create_mol_catalog, m)?)?;
    Ok(())
}
