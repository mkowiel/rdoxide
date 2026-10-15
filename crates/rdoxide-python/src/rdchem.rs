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

pyo3::create_exception!(
    rdkit.Chem.rdchem,
    MolSanitizeException,
    pyo3::exceptions::PyValueError
);
pyo3::create_exception!(
    rdkit.Chem.rdchem,
    AtomSanitizeException,
    MolSanitizeException
);
pyo3::create_exception!(
    rdkit.Chem.rdchem,
    AtomValenceException,
    AtomSanitizeException
);
pyo3::create_exception!(
    rdkit.Chem.rdchem,
    AtomKekulizeException,
    AtomSanitizeException
);
pyo3::create_exception!(rdkit.Chem.rdchem, KekulizeException, MolSanitizeException);

/// exception arising from sanitization (the `cause` of the Python exceptions)
#[pyclass(
    name = "_cppMolSanitizeException",
    module = "rdkit.Chem.rdchem",
    frozen
)]
#[derive(Debug)]
pub struct SanitizeCause {
    pub kind: &'static str,
    pub message: String,
    pub atom_idx: Option<u32>,
    pub atom_indices: Vec<u32>,
}

#[pymethods]
impl SanitizeCause {
    #[pyo3(name = "Message")]
    fn message(&self) -> String {
        self.message.clone()
    }
    #[pyo3(name = "GetType")]
    fn get_type(&self) -> &'static str {
        self.kind
    }
    #[pyo3(name = "GetAtomIdx")]
    fn get_atom_idx(&self) -> PyResult<u32> {
        self.atom_idx
            .ok_or_else(|| pyo3::exceptions::PyAttributeError::new_err("GetAtomIdx"))
    }
    #[pyo3(name = "GetAtomIndices")]
    fn get_atom_indices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        pyo3::types::PyTuple::new(py, &self.atom_indices)
    }
}

/// Build the Python exception for a sanitization failure.
pub fn sanitize_error(err: &rdoxide::SanitizeError) -> PyErr {
    use rdoxide::SanitizeError as S;
    Python::attach(|py| {
        let (kind, message, atom_idx, atom_indices) = match err {
            S::Generic { message } => ("MolSanitizeException", message.clone(), None, vec![]),
            S::AtomValence { atom_idx, message } => (
                "AtomValenceException",
                message.clone(),
                Some(*atom_idx),
                vec![],
            ),
            S::AtomKekulize { atom_idx, message } => (
                "AtomKekulizeException",
                message.clone(),
                Some(*atom_idx),
                vec![],
            ),
            S::Kekulize {
                atom_indices,
                message,
            } => (
                "KekulizeException",
                message.clone(),
                None,
                atom_indices.clone(),
            ),
        };
        let exc_type = match err {
            S::Generic { .. } => py.get_type::<MolSanitizeException>(),
            S::AtomValence { .. } => py.get_type::<AtomValenceException>(),
            S::AtomKekulize { .. } => py.get_type::<AtomKekulizeException>(),
            S::Kekulize { .. } => py.get_type::<KekulizeException>(),
        };
        let build = || -> PyResult<PyErr> {
            let cause = Py::new(
                py,
                SanitizeCause {
                    kind,
                    message: message.clone(),
                    atom_idx,
                    atom_indices,
                },
            )?;
            let inst = exc_type.call1((message.clone(),))?;
            inst.setattr("cause", &cause)?;
            // upstream sets `cause` on the exception type as well
            exc_type.setattr("cause", &cause)?;
            Ok(PyErr::from_value(inst))
        };
        build().unwrap_or_else(|e| e)
    })
}

#[pyfunction]
fn tossit() -> PyResult<()> {
    Err(pyo3::exceptions::PyIndexError::new_err("1"))
}

/// Get the current global mol pickler options.
#[pyfunction]
#[pyo3(name = "GetDefaultPickleProperties")]
fn get_default_pickle_properties() -> u32 {
    rdoxide::graphmol::pickle::get_default_pickle_properties()
}

/// Set the current global mol pickler options.
#[pyfunction]
#[pyo3(name = "SetDefaultPickleProperties")]
fn set_default_pickle_properties(arg1: u32) {
    rdoxide::graphmol::pickle::set_default_pickle_properties(arg1)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing the core chemistry functionality of the RDKit",
    )?;
    m.add_class::<PyPeriodicTable>()?;
    crate::substruct::init(m)?;
    crate::rdresonance::init(m)?;
    m.add_function(wrap_pyfunction!(get_periodic_table, m)?)?;

    crate::rdchem_stereo_info::add_chem(m)?;
    use crate::enums::{add_enum, add_enum_exported};
    use crate::enums_impl::*;
    add_enum(
        m,
        "HybridizationType",
        &PyHybridizationType::members(),
        Some("Hybridization type"),
    )?;
    add_enum_exported(m, "ChiralType", &PyChiralType::members(), None)?;
    add_enum_exported(m, "ValenceType", &PyValenceType::members(), None)?;
    add_enum(
        m,
        "AtomMonomerType",
        &crate::monomer_info::PyAtomMonomerType::members(),
        None,
    )?;
    add_enum(m, "BondType", &PyBondType::members(), None)?;
    add_enum(m, "BondDir", &PyBondDir::members(), None)?;
    add_enum(m, "BondStereo", &PyBondStereo::members(), None)?;
    add_enum_exported(m, "StereoGroupType", &PyStereoGroupType::members(), None)?;
    add_enum_exported(
        m,
        "CompositeQueryType",
        &[
            ("COMPOSITE_AND", 0),
            ("COMPOSITE_OR", 1),
            ("COMPOSITE_XOR", 2),
        ],
        None,
    )?;
    add_enum_exported(
        m,
        "ResonanceFlags",
        &[
            ("ALLOW_INCOMPLETE_OCTETS", 1),
            ("ALLOW_CHARGE_SEPARATION", 2),
            ("KEKULE_ALL", 4),
            ("UNCONSTRAINED_CATIONS", 8),
            ("UNCONSTRAINED_ANIONS", 16),
        ],
        None,
    )?;
    add_enum_exported(
        m,
        "PropertyPickleOptions",
        &[
            ("NoProps", 0),
            ("MolProps", 0x1),
            ("AtomProps", 0x2),
            ("BondProps", 0x4),
            ("QueryAtomData", 0x2),
            ("PrivateProps", 0x10),
            ("ComputedProps", 0x20),
            ("AllProps", 0xFFFF),
            ("CoordsAsDouble", 0x10000),
            ("NoConformers", 0x20000),
        ],
        None,
    )?;
    m.add_function(wrap_pyfunction!(get_default_pickle_properties, m)?)?;
    m.add_function(wrap_pyfunction!(set_default_pickle_properties, m)?)?;

    let py = m.py();
    m.add_class::<SanitizeCause>()?;
    m.add(
        "MolSanitizeException",
        py.get_type::<MolSanitizeException>(),
    )?;
    m.add(
        "AtomSanitizeException",
        py.get_type::<AtomSanitizeException>(),
    )?;
    m.add(
        "AtomValenceException",
        py.get_type::<AtomValenceException>(),
    )?;
    m.add(
        "AtomKekulizeException",
        py.get_type::<AtomKekulizeException>(),
    )?;
    m.add("KekulizeException", py.get_type::<KekulizeException>())?;

    m.add_class::<crate::atom::PyAtom>()?;
    m.add_class::<crate::atom::PyQueryAtom>()?;
    crate::rdchem_atom_funcs::add(m)?;
    crate::monomer_info::init(m)?;
    m.add_class::<crate::bond::PyBond>()?;
    m.add_class::<crate::bond::PyQueryBond>()?;
    m.add_class::<crate::conformer::PyConformer>()?;
    m.add_class::<crate::mol::PyMol>()?;
    m.add_class::<crate::mol::PyRWMol>()?;
    m.add_class::<crate::mol::PyEditableMol>()?;
    m.add_class::<crate::mol::AtomSeq>()?;
    m.add_class::<crate::mol::BondSeq>()?;
    m.add_class::<crate::ringinfo::PyRingInfo>()?;
    m.add_class::<crate::stereo_group::PyStereoGroup>()?;
    m.add_function(wrap_pyfunction!(
        crate::stereo_group::create_stereo_group,
        m
    )?)?;
    crate::substance_group::add(m)?;
    m.add_function(wrap_pyfunction!(tossit, m)?)?;
    Ok(())
}
