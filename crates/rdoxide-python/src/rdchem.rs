//! `rdkit.Chem.rdchem` — port of upstream `Code/GraphMol/Wrap/rdchem.cpp`
//! and `Code/GraphMol/Wrap/Table.cpp`.

use pyo3::prelude::*;
use rdoxide::graphmol::PeriodicTable;
use rdoxide::graphmol::periodic_table::{AtomicData, Element};

use crate::errors::IntoPyResult;

/// Accepts an atomic number or an element symbol, like upstream's overloads.
#[derive(FromPyObject)]
enum ElementArg {
    Num(u32),
    Symbol(String),
}

fn lookup(elem: &ElementArg) -> PyResult<&'static AtomicData> {
    let pt = PeriodicTable::get();
    match elem {
        ElementArg::Num(n) => pt.resolve(Element::Num(*n)),
        ElementArg::Symbol(s) => pt.resolve(Element::Symbol(s)),
    }
    .py()
}

/// A class which stores information from the Periodic Table.
///
///   It is not possible to create a PeriodicTable object directly from Python,
///   use GetPeriodicTable() to get the global table.
///
///   The PeriodicTable object can be queried for a variety of properties:
///
///     - GetAtomicWeight
///
///     - GetAtomicNumber
///
///     - GetElementSymbol
///
///     - GetElementName
///
///     - GetRvdw (van der Waals radius)
///
///     - GetRCovalent (covalent radius)
///
///     - GetDefaultValence
///
///     - GetValenceList
///
///     - GetNOuterElecs (number of valence electrons)
///
///     - GetMostCommonIsotope
///
///     - GetMostCommonIsotopeMass
///
///     - GetRb0
///
///     - GetAbundanceForIsotope
///
///     - GetMassForIsotope
///
///   When it makes sense, these can be queried using either an atomic number (integer)
///   or an atomic symbol (string)
#[pyclass(name = "PeriodicTable", module = "rdkit.Chem.rdchem", frozen)]
#[derive(Debug)]
pub struct PyPeriodicTable;

#[pymethods]
impl PyPeriodicTable {
    #[pyo3(name = "GetAtomicWeight")]
    fn get_atomic_weight(&self, elem: ElementArg) -> PyResult<f64> {
        Ok(lookup(&elem)?.mass)
    }
    #[pyo3(name = "GetAtomicNumber")]
    fn get_atomic_number(&self, symbol: &str) -> PyResult<u32> {
        PeriodicTable::get().get_atomic_number(symbol).py()
    }
    #[pyo3(name = "GetElementSymbol")]
    fn get_element_symbol(&self, atomic_number: u32) -> PyResult<&'static str> {
        PeriodicTable::get().get_element_symbol(atomic_number).py()
    }
    #[pyo3(name = "GetElementName")]
    fn get_element_name(&self, atomic_number: u32) -> PyResult<&'static str> {
        PeriodicTable::get().get_element_name(atomic_number).py()
    }
    #[pyo3(name = "GetRvdw")]
    fn get_rvdw(&self, elem: ElementArg) -> PyResult<f64> {
        Ok(lookup(&elem)?.r_vdw)
    }
    #[pyo3(name = "GetRcovalent")]
    fn get_rcovalent(&self, elem: ElementArg) -> PyResult<f64> {
        Ok(lookup(&elem)?.r_cov)
    }
    #[pyo3(name = "GetRb0")]
    fn get_rb0(&self, elem: ElementArg) -> PyResult<f64> {
        Ok(lookup(&elem)?.r_b0)
    }
    #[pyo3(name = "GetDefaultValence")]
    fn get_default_valence(&self, elem: ElementArg) -> PyResult<i32> {
        Ok(lookup(&elem)?.valences[0])
    }
    #[pyo3(name = "GetValenceList")]
    fn get_valence_list(&self, elem: ElementArg) -> PyResult<Vec<i32>> {
        Ok(lookup(&elem)?.valences.clone())
    }
    #[pyo3(name = "GetNOuterElecs")]
    fn get_n_outer_elecs(&self, elem: ElementArg) -> PyResult<i32> {
        Ok(lookup(&elem)?.n_outer_elecs)
    }
    #[pyo3(name = "GetMostCommonIsotope")]
    fn get_most_common_isotope(&self, elem: ElementArg) -> PyResult<i32> {
        Ok(lookup(&elem)?.common_isotope)
    }
    #[pyo3(name = "GetMostCommonIsotopeMass")]
    fn get_most_common_isotope_mass(&self, elem: ElementArg) -> PyResult<f64> {
        Ok(lookup(&elem)?.common_isotope_mass)
    }
    #[pyo3(name = "GetAbundanceForIsotope")]
    fn get_abundance_for_isotope(&self, elem: ElementArg, isotope: u32) -> PyResult<f64> {
        Ok(lookup(&elem)?.isotopes.get(&isotope).map_or(0.0, |x| x.1))
    }
    #[pyo3(name = "GetMassForIsotope")]
    fn get_mass_for_isotope(&self, elem: ElementArg, isotope: u32) -> PyResult<f64> {
        Ok(lookup(&elem)?.isotopes.get(&isotope).map_or(0.0, |x| x.0))
    }
    #[pyo3(name = "GetMaxAtomicNumber")]
    fn get_max_atomic_number(&self) -> u32 {
        PeriodicTable::get().get_max_atomic_number()
    }
    #[pyo3(name = "GetRow")]
    fn get_row(&self, elem: ElementArg) -> PyResult<u32> {
        Ok(lookup(&elem)?.row)
    }
    #[pyo3(name = "MoreElectroNegative")]
    fn more_electronegative(&self, anum1: u32, anum2: u32) -> PyResult<bool> {
        PeriodicTable::get().more_electronegative(anum1, anum2).py()
    }
}

#[pyfunction]
#[pyo3(name = "GetPeriodicTable")]
/// Returns the application's PeriodicTable instance.
fn get_periodic_table(py: Python<'_>) -> PyResult<Py<PyPeriodicTable>> {
    static TABLE: pyo3::sync::PyOnceLock<Py<PyPeriodicTable>> = pyo3::sync::PyOnceLock::new();
    TABLE
        .get_or_try_init(py, || Py::new(py, PyPeriodicTable))
        .map(|t| t.clone_ref(py))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing the core chemistry functionality of the RDKit",
    )?;
    m.add_class::<PyPeriodicTable>()?;
    m.add_function(wrap_pyfunction!(get_periodic_table, m)?)?;
    Ok(())
}
