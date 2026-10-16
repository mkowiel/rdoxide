//! `rdkit.Chem.rdDeprotect` (upstream `GraphMol/Deprotect/Wrap/`).

use pyo3::prelude::*;
use rdoxide::graphmol::deprotect::{self, DeprotectData};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// DeprotectData class, contains a single deprotection reaction and information
///
///  deprotectdata.deprotection_class - functional group being protected
///  deprotectdata.reaction_smarts - reaction smarts used for deprotection
///  deprotectdata.abbreviation - common abbreviation for the protecting group
///  deprotectdata.full_name - full name for the protecting group
#[pyclass(
    name = "DeprotectData",
    module = "rdkit.Chem.rdDeprotect",
    skip_from_py_object
)]
#[derive(Clone)]
struct PyDeprotectData {
    data: DeprotectData,
}

#[pymethods]
impl PyDeprotectData {
    /// Construct a new DeprotectData instance.
    ///   >>> reaction_class = "amine"
    ///   >>> reaction_smarts = "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]"
    ///   >>> abbreviation = "Boc"
    ///   >>> full_name = "tert-butyloxycarbonyl"
    ///   >>> data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
    ///   >>> assert data.isValid()
    #[new]
    #[pyo3(signature = (deprotection_class, reaction_smarts, abbreviation, full_name))]
    fn new(
        deprotection_class: &str,
        reaction_smarts: &str,
        abbreviation: &str,
        full_name: &str,
    ) -> PyResult<Self> {
        let data = DeprotectData::new(
            deprotection_class,
            reaction_smarts,
            abbreviation,
            full_name,
            "",
        )
        .py()?;
        Ok(PyDeprotectData { data })
    }

    #[getter]
    fn deprotection_class(&self) -> &str {
        &self.data.deprotection_class
    }

    #[getter]
    fn full_name(&self) -> &str {
        &self.data.full_name
    }

    #[getter]
    fn abbreviation(&self) -> &str {
        &self.data.abbreviation
    }

    #[getter]
    fn reaction_smarts(&self) -> &str {
        &self.data.reaction_smarts
    }

    #[getter]
    fn example(&self) -> &str {
        &self.data.example
    }

    /// Returns True if the DeprotectData has a valid reaction
    #[pyo3(name = "isValid")]
    fn is_valid(&self) -> bool {
        self.data.is_valid()
    }
}

/// The deprotections of an optional Python iterable (`None`: the defaults).
fn deprotections(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Vec<DeprotectData>> {
    match obj {
        None => Ok(deprotect::get_deprotections().to_vec()),
        Some(o) if o.is_none() => Ok(deprotect::get_deprotections().to_vec()),
        Some(o) => o
            .try_iter()?
            .map(|item| Ok(item?.extract::<PyRef<'_, PyDeprotectData>>()?.data.clone()))
            .collect(),
    }
}

/// Return the default list of deprotections
#[pyfunction]
#[pyo3(name = "GetDeprotections")]
fn py_get_deprotections() -> Vec<PyDeprotectData> {
    deprotect::get_deprotections()
        .iter()
        .map(|d| PyDeprotectData { data: d.clone() })
        .collect()
}

/// Return the deprotected version of the molecule.
#[pyfunction]
#[pyo3(name = "Deprotect", signature = (mol, deprotections = None))]
fn py_deprotect(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    deprotections: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    let list = self::deprotections(deprotections)?;
    let source = &mol.mol;
    let res = py.detach(|| deprotect::deprotect(source, &list)).py()?;
    Ok(PyMol { mol: res })
}

/// Deprotects the molecule in place.
#[pyfunction]
#[pyo3(name = "DeprotectInPlace", signature = (mol, deprotections = None))]
fn py_deprotect_in_place(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    deprotections: Option<&Bound<'_, PyAny>>,
) -> PyResult<bool> {
    let list = self::deprotections(deprotections)?;
    let mol = &mut mol.mol;
    py.detach(|| deprotect::deprotect_in_place(mol, &list)).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDeprotectData>()?;
    m.add_function(wrap_pyfunction!(py_get_deprotections, m)?)?;
    m.add_function(wrap_pyfunction!(py_deprotect, m)?)?;
    m.add_function(wrap_pyfunction!(py_deprotect_in_place, m)?)?;
    Ok(())
}
