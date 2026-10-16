//! `rdoxide.Chem.rdfragcatalog` (upstream `FragCatalog/Wrap`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use rdoxide::catalogs::CatalogEntry;
use rdoxide::graphmol::frag_catalog::{
    FragCatGenerator, FragCatParams, FragCatalog, FragFPGenerator,
};

use crate::datastructs::PyExplicitBitVect;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Parameters of a fragment catalog.
#[pyclass(
    name = "FragCatParams",
    module = "rdoxide.Chem.rdfragcatalog",
    from_py_object
)]
#[derive(Clone)]
pub struct PyFragCatParams {
    params: FragCatParams,
}

#[pymethods]
impl PyFragCatParams {
    #[new]
    #[pyo3(signature = (lLen, uLen, fgroupFilename, tol = 1e-8))]
    #[allow(non_snake_case)]
    fn py_new(lLen: u32, uLen: u32, fgroupFilename: &str, tol: f64) -> PyResult<Self> {
        Ok(Self {
            params: FragCatParams::new(lLen, uLen, fgroupFilename, tol).py()?,
        })
    }

    #[pyo3(name = "GetTypeString")]
    fn get_type_string(&self) -> String {
        rdoxide::catalogs::CatalogParams::type_str(&self.params).to_string()
    }

    #[pyo3(name = "GetUpperFragLength")]
    fn get_upper_frag_length(&self) -> u32 {
        self.params.upper_frag_length()
    }

    #[pyo3(name = "GetLowerFragLength")]
    fn get_lower_frag_length(&self) -> u32 {
        self.params.lower_frag_length()
    }

    #[pyo3(name = "GetTolerance")]
    fn get_tolerance(&self) -> f64 {
        self.params.tolerance()
    }

    #[pyo3(name = "GetNumFuncGroups")]
    fn get_num_func_groups(&self) -> u32 {
        self.params.num_func_groups()
    }

    #[pyo3(name = "GetFuncGroup")]
    fn get_func_group(&self, fid: u32) -> PyResult<PyMol> {
        Ok(PyMol {
            mol: self.params.func_group(fid).py()?.clone(),
        })
    }

    #[pyo3(name = "Serialize")]
    fn serialize(&self) -> PyResult<String> {
        self.params.serialize().py()
    }
}

/// A fragment catalog.
#[pyclass(name = "FragCatalog", module = "rdoxide.Chem.rdfragcatalog")]
pub struct PyFragCatalog {
    cat: FragCatalog,
}

impl PyFragCatalog {
    fn bit_entry(&self, idx: u32) -> PyResult<&rdoxide::graphmol::frag_catalog::FragCatalogEntry> {
        self.cat
            .entry_with_bit_id(idx)
            .py()?
            .ok_or_else(|| PyValueError::new_err("no entry has this bit id"))
    }
}

fn func_group_ids(entry: &rdoxide::graphmol::frag_catalog::FragCatalogEntry) -> Vec<i32> {
    entry.func_group_map().values().flatten().copied().collect()
}

#[pymethods]
impl PyFragCatalog {
    /// `FragCatalog(params)` or `FragCatalog(pickle)`
    #[new]
    fn py_new(arg: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(p) = arg.extract::<PyFragCatParams>() {
            return Ok(Self {
                cat: FragCatalog::new(p.params),
            });
        }
        let pkl: Vec<u8> = arg.extract()?;
        Ok(Self {
            cat: FragCatalog::from_pickle(&pkl).py()?,
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

    #[pyo3(name = "GetCatalogParams")]
    fn get_catalog_params(&self) -> PyResult<PyFragCatParams> {
        self.cat
            .catalog_params()
            .cloned()
            .map(|params| PyFragCatParams { params })
            .ok_or_else(|| PyValueError::new_err("the catalog has no parameters"))
    }

    #[pyo3(name = "Serialize")]
    fn serialize<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, &self.cat.serialize().py()?))
    }

    #[pyo3(name = "GetBitDescription")]
    fn get_bit_description(&self, idx: u32) -> PyResult<String> {
        Ok(self.bit_entry(idx)?.description())
    }

    #[pyo3(name = "GetBitOrder")]
    fn get_bit_order(&self, idx: u32) -> PyResult<i32> {
        Ok(self.bit_entry(idx)?.order())
    }

    #[pyo3(name = "GetBitFuncGroupIds")]
    fn get_bit_func_group_ids(&self, idx: u32) -> PyResult<Vec<i32>> {
        Ok(func_group_ids(self.bit_entry(idx)?))
    }

    #[pyo3(name = "GetBitEntryId")]
    fn get_bit_entry_id(&self, idx: u32) -> PyResult<i32> {
        self.cat.id_of_entry_with_bit_id(idx).py()
    }

    #[pyo3(name = "GetEntryBitId")]
    fn get_entry_bit_id(&self, idx: u32) -> PyResult<i32> {
        Ok(self.cat.entry_with_idx(idx).py()?.bit_id())
    }

    #[pyo3(name = "GetEntryDescription")]
    fn get_entry_description(&self, idx: u32) -> PyResult<String> {
        Ok(self.cat.entry_with_idx(idx).py()?.description())
    }

    #[pyo3(name = "GetEntryOrder")]
    fn get_entry_order(&self, idx: u32) -> PyResult<i32> {
        Ok(self.cat.entry_with_idx(idx).py()?.order())
    }

    #[pyo3(name = "GetEntryFuncGroupIds")]
    fn get_entry_func_group_ids(&self, idx: u32) -> PyResult<Vec<i32>> {
        Ok(func_group_ids(self.cat.entry_with_idx(idx).py()?))
    }

    #[pyo3(name = "GetEntryDownIds")]
    fn get_entry_down_ids(&self, idx: u32) -> PyResult<Vec<u32>> {
        Ok(self.cat.down_entry_list(idx).py()?.to_vec())
    }

    #[pyo3(name = "GetBitDiscrims")]
    fn get_bit_discrims(&self, idx: u32) -> PyResult<Vec<f64>> {
        let (a, b, c) = self.bit_entry(idx)?.discrims().py()?;
        Ok(vec![f64::from(a), b as f64, c as f64])
    }

    fn __reduce__<'py>(
        slf: &Bound<'py, Self>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, pyo3::types::PyTuple>)> {
        let py = slf.py();
        let pkl = PyBytes::new(py, &slf.borrow().cat.serialize().py()?);
        Ok((
            py.get_type::<Self>().into_any(),
            pyo3::types::PyTuple::new(py, [pkl])?,
        ))
    }
}

/// Generates the fragments of molecules into a catalog.
#[pyclass(name = "FragCatGenerator", module = "rdoxide.Chem.rdfragcatalog")]
pub struct PyFragCatGenerator;

#[pymethods]
impl PyFragCatGenerator {
    #[new]
    fn py_new() -> Self {
        PyFragCatGenerator
    }

    #[pyo3(name = "AddFragsFromMol")]
    fn add_frags_from_mol(
        &self,
        mol: &Bound<'_, PyMol>,
        fcat: &Bound<'_, PyFragCatalog>,
    ) -> PyResult<u32> {
        FragCatGenerator
            .add_frags_from_mol(&mol.borrow().mol, &mut fcat.borrow_mut().cat)
            .py()
    }
}

/// Generates fingerprints from a fragment catalog.
#[pyclass(name = "FragFPGenerator", module = "rdoxide.Chem.rdfragcatalog")]
pub struct PyFragFPGenerator;

#[pymethods]
impl PyFragFPGenerator {
    #[new]
    fn py_new() -> Self {
        PyFragFPGenerator
    }

    #[pyo3(name = "GetFPForMol")]
    fn get_fp_for_mol(
        &self,
        py: Python<'_>,
        mol: &Bound<'_, PyMol>,
        fcat: &Bound<'_, PyFragCatalog>,
    ) -> PyResult<PyExplicitBitVect> {
        let mol = &mol.borrow().mol;
        let fcat_ref = fcat.borrow();
        let cat = &fcat_ref.cat;
        let fp = py.detach(|| FragFPGenerator.fp_for_mol(mol, cat)).py()?;
        Ok(PyExplicitBitVect::wrap(fp))
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFragCatParams>()?;
    m.add_class::<PyFragCatalog>()?;
    m.add_class::<PyFragCatGenerator>()?;
    m.add_class::<PyFragFPGenerator>()?;
    Ok(())
}
