//! `rdoxide.Chem.rdFMCS` (upstream `FMCS/Wrap/rdFMCS.cpp`): the maximum common
//! substructure of a set of molecules.

use std::sync::{Arc, Mutex};

use pyo3::exceptions::{PyAttributeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::fmcs::{
    self, AtomComparator, AtomTyper, BondComparator, BondTyper, Graph, MCSAtomCompareParameters,
    MCSBondCompareParameters, MCSParameters, MCSProgressData, MCSResult, RingComparator,
};
use rdoxide::{RDKitError, Result as CoreResult};

use crate::enums::{add_enum, to_py};
use crate::mol::PyMol;
use crate::pyconv::set_declared_attribute;
use crate::rdmolfiles_write::Args;

const ATOM_COMPARE: &str = "AtomCompare";
const BOND_COMPARE: &str = "BondCompare";
const RING_COMPARE: &str = "RingCompare";

fn atom_comparator_value(c: AtomComparator) -> i64 {
    c as i64
}

fn atom_comparator_from(v: i64) -> PyResult<AtomComparator> {
    Ok(match v {
        0 => AtomComparator::Any,
        1 => AtomComparator::Elements,
        2 => AtomComparator::Isotopes,
        3 => AtomComparator::AnyHeavyAtom,
        _ => return Err(PyValueError::new_err("Unknown AtomComparator")),
    })
}

fn bond_comparator_from(v: i64) -> PyResult<BondComparator> {
    Ok(match v {
        0 => BondComparator::Any,
        1 => BondComparator::Order,
        2 => BondComparator::OrderExact,
        _ => return Err(PyValueError::new_err("Unknown BondComparator")),
    })
}

fn ring_comparator_from(v: i64) -> PyResult<RingComparator> {
    Ok(match v {
        0 => RingComparator::IgnoreRingFusion,
        1 => RingComparator::PermissiveRingFusion,
        2 => RingComparator::StrictRingFusion,
        _ => return Err(PyValueError::new_err("Unknown RingComparator")),
    })
}

// ------------------------------------------------------------------ parameter classes

macro_rules! bool_props {
    ($class:ident, $inner:ident, $( ($py:ident, $get:ident, $set:ident, $field:ident, $doc:literal) ),* $(,)?) => {
        #[pymethods]
        #[allow(non_snake_case)]
        impl $class {
            $(
                #[doc = $doc]
                #[getter($py)]
                fn $get(&self) -> bool {
                    self.p.$field
                }
                #[setter($py)]
                fn $set(&mut self, value: bool) {
                    self.p.$field = value;
                }
            )*
        }
    };
}

/// Parameters controlling how atom-atom matching is done
#[pyclass(
    name = "MCSAtomCompareParameters",
    module = "rdoxide.Chem.rdFMCS",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyAtomCompareParameters {
    p: MCSAtomCompareParameters,
}

bool_props!(
    PyAtomCompareParameters,
    p,
    (
        MatchValences,
        get_match_valences,
        set_match_valences,
        match_valences,
        "include atom valences in the match"
    ),
    (
        MatchChiralTag,
        get_match_chiral_tag,
        set_match_chiral_tag,
        match_chiral_tag,
        "include atom chirality in the match"
    ),
    (
        MatchFormalCharge,
        get_match_formal_charge,
        set_match_formal_charge,
        match_formal_charge,
        "include formal charge in the match"
    ),
    (
        RingMatchesRingOnly,
        get_ring_matches_ring_only,
        set_ring_matches_ring_only,
        ring_matches_ring_only,
        "ring atoms are only allowed to match other ring atoms"
    ),
    (
        CompleteRingsOnly,
        get_complete_rings_only,
        set_complete_rings_only,
        complete_rings_only,
        "results cannot include lone ring atoms"
    ),
    (
        MatchIsotope,
        get_match_isotope,
        set_match_isotope,
        match_isotope,
        "use isotope atom queries in MCSResults"
    ),
);

#[pymethods]
#[allow(non_snake_case)]
impl PyAtomCompareParameters {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Require atoms to be within this many angstroms in 3D
    #[getter(MaxDistance)]
    fn get_max_distance(&self) -> f64 {
        self.p.max_distance
    }
    #[setter(MaxDistance)]
    fn set_max_distance(&mut self, value: f64) {
        self.p.max_distance = value;
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        set_declared_attribute(slf.as_any(), name, value)
    }
}

/// Parameters controlling how bond-bond matching is done
#[pyclass(
    name = "MCSBondCompareParameters",
    module = "rdoxide.Chem.rdFMCS",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyBondCompareParameters {
    p: MCSBondCompareParameters,
}

bool_props!(
    PyBondCompareParameters,
    p,
    (
        RingMatchesRingOnly,
        get_ring_matches_ring_only,
        set_ring_matches_ring_only,
        ring_matches_ring_only,
        "ring bonds are only allowed to match other ring bonds"
    ),
    (
        CompleteRingsOnly,
        get_complete_rings_only,
        set_complete_rings_only,
        complete_rings_only,
        "results cannot include partial rings"
    ),
    (
        MatchFusedRings,
        get_match_fused_rings,
        set_match_fused_rings,
        match_fused_rings,
        "enforce check on ring fusion, i.e. alpha-methylnaphthalene\nwon't match beta-methylnaphtalene, but decalin\nwill match cyclodecane unless MatchFusedRingsStrict is True"
    ),
    (
        MatchFusedRingsStrict,
        get_match_fused_rings_strict,
        set_match_fused_rings_strict,
        match_fused_rings_strict,
        "only enforced if MatchFusedRings is True; the ring fusion\nmust be the same in both query and target, i.e. decalin\nwon't match cyclodecane"
    ),
    (
        MatchStereo,
        get_match_stereo,
        set_match_stereo,
        match_stereo,
        "include bond stereo in the comparison"
    ),
);

#[pymethods]
impl PyBondCompareParameters {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        set_declared_attribute(slf.as_any(), name, value)
    }
}

/// Information about the MCS progress
#[pyclass(
    name = "MCSProgressData",
    module = "rdoxide.Chem.rdFMCS",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyProgressData {
    d: MCSProgressData,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyProgressData {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// number of atoms in MCS
    #[getter(numAtoms)]
    fn num_atoms(&self) -> u32 {
        self.d.num_atoms
    }

    /// number of bonds in MCS
    #[getter(numBonds)]
    fn num_bonds(&self) -> u32 {
        self.d.num_bonds
    }

    /// number of processed seeds
    #[getter(seedProcessed)]
    fn seed_processed(&self) -> u32 {
        self.d.seed_processed
    }
}

// ------------------------------------------------------------------ base classes of the callbacks

/// Base class. Subclass and override
/// MCSAtomCompare.__call__() to define custom
/// atom compare functions, then set MCSParameters.AtomTyper
/// to an instance of the subclass
#[pyclass(name = "MCSAtomCompare", module = "rdoxide.Chem.rdFMCS", subclass)]
pub struct PyAtomCompare;

#[pymethods]
#[allow(non_snake_case)]
impl PyAtomCompare {
    #[new]
    #[pyo3(signature = (*_args, **_kwargs))]
    fn new(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
        PyAtomCompare
    }

    /// Return True if both atoms are, or are not, in a ring
    #[pyo3(signature = (parameters, mol1, atom1, mol2, atom2))]
    fn CheckAtomRingMatch(
        &self,
        parameters: &PyAtomCompareParameters,
        mol1: &PyMol,
        atom1: u32,
        mol2: &PyMol,
        atom2: u32,
    ) -> bool {
        fmcs::check_atom_ring_match(&parameters.p, &mol1.mol, atom1, &mol2.mol, atom2)
    }

    /// Return True if both atoms have the same formal charge
    #[pyo3(signature = (parameters, mol1, atom1, mol2, atom2))]
    fn CheckAtomCharge(
        &self,
        parameters: &PyAtomCompareParameters,
        mol1: &PyMol,
        atom1: u32,
        mol2: &PyMol,
        atom2: u32,
    ) -> bool {
        let _ = parameters;
        fmcs::check_atom_charge(&mol1.mol, atom1, &mol2.mol, atom2)
    }

    /// Return True if both atoms have, or have not, a chiral tag
    #[pyo3(signature = (parameters, mol1, atom1, mol2, atom2))]
    fn CheckAtomChirality(
        &self,
        parameters: &PyAtomCompareParameters,
        mol1: &PyMol,
        atom1: u32,
        mol2: &PyMol,
        atom2: u32,
    ) -> bool {
        let _ = parameters;
        fmcs::check_atom_chirality(&mol1.mol, atom1, &mol2.mol, atom2)
    }

    /// override to implement custom atom comparison
    #[pyo3(signature = (parameters, mol1, atom1, mol2, atom2))]
    fn __call__(
        &self,
        parameters: &PyAtomCompareParameters,
        mol1: &PyMol,
        atom1: u32,
        mol2: &PyMol,
        atom2: u32,
    ) -> PyResult<bool> {
        let _ = (parameters, mol1, atom1, mol2, atom2);
        Err(PyAttributeError::new_err(
            "The __call__() method must be overridden in the rdFMCS.MCSAtomCompare subclass",
        ))
    }
}

/// Base class. Subclass and override
/// MCSBondCompare.__call__() to define custom
/// bond compare functions, then set MCSParameters.BondTyper
/// to an instance of the subclass
#[pyclass(name = "MCSBondCompare", module = "rdoxide.Chem.rdFMCS", subclass)]
pub struct PyBondCompare;

#[pymethods]
#[allow(non_snake_case)]
impl PyBondCompare {
    #[new]
    #[pyo3(signature = (*_args, **_kwargs))]
    fn new(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
        PyBondCompare
    }

    /// Return True if both bonds have, or have not, a stereo descriptor
    #[pyo3(signature = (parameters, mol1, bond1, mol2, bond2))]
    fn CheckBondStereo(
        &self,
        parameters: &PyBondCompareParameters,
        mol1: &PyMol,
        bond1: u32,
        mol2: &PyMol,
        bond2: u32,
    ) -> bool {
        let _ = parameters;
        fmcs::check_bond_stereo(&mol1.mol, bond1, &mol2.mol, bond2)
    }

    /// Return True if both bonds are, or are not, part of a ring
    #[pyo3(signature = (parameters, mol1, bond1, mol2, bond2))]
    fn CheckBondRingMatch(
        &self,
        parameters: &PyBondCompareParameters,
        mol1: &PyMol,
        bond1: u32,
        mol2: &PyMol,
        bond2: u32,
    ) -> bool {
        fmcs::check_bond_ring_match(&parameters.p, &mol1.mol, bond1, &mol2.mol, bond2)
    }

    /// override to implement custom bond comparison
    #[pyo3(signature = (parameters, mol1, bond1, mol2, bond2))]
    fn __call__(
        &self,
        parameters: &PyBondCompareParameters,
        mol1: &PyMol,
        bond1: u32,
        mol2: &PyMol,
        bond2: u32,
    ) -> PyResult<bool> {
        let _ = (parameters, mol1, bond1, mol2, bond2);
        Err(PyAttributeError::new_err(
            "The __call__() method must be overridden in the rdFMCS.MCSBondCompare subclass",
        ))
    }
}

macro_rules! callback_base {
    ($class:ident, $name:literal, $doc:literal, $subclass:literal, ($($arg:ident: $ty:ty),*)) => {
        #[doc = $doc]
        #[pyclass(name = $name, module = "rdoxide.Chem.rdFMCS", subclass)]
        pub struct $class;

        #[pymethods]
        impl $class {
            #[new]
            #[pyo3(signature = (*_args, **_kwargs))]
            fn new(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
                $class
            }

            #[pyo3(signature = (*_args, **_kwargs))]
            fn __call__(
                &self,
                _args: &Bound<'_, PyTuple>,
                _kwargs: Option<&Bound<'_, PyDict>>,
            ) -> PyResult<bool> {
                Err(PyAttributeError::new_err(concat!(
                    "The __call__() method must be overridden in the rdFMCS.",
                    $subclass,
                    " subclass"
                )))
            }
        }
    };
}

callback_base!(
    PyProgress,
    "MCSProgress",
    "Base class. Subclass and override\nMCSProgress.__call__()\nto define a custom callback function",
    "MCSProgress",
    ()
);
callback_base!(
    PyFinalMatchCheck,
    "MCSFinalMatchCheck",
    "Base class. Subclass and override\nMCSFinalMatchCheck.__call__()\nto define a custom boolean callback function.\nReturning True will cause the growing seed to be accepted,\nFalse to be rejected",
    "MCSFinalMatchCheck",
    ()
);
callback_base!(
    PyAcceptance,
    "MCSAcceptance",
    "Base class. Subclass and override\nMCSAcceptance.__call__()\nto define a custom boolean callback function.\nReturning True will cause the MCS candidate to be accepted,\nFalse to be rejected",
    "MCSAcceptance",
    ()
);

// ------------------------------------------------------------------ python callbacks

/// The Python objects a parameter set was configured with.
#[derive(Debug, Default)]
struct PyCallbacks {
    atom_comp: Option<Py<PyAny>>,
    bond_comp: Option<Py<PyAny>>,
    progress: Option<Py<PyAny>>,
    final_check: Option<Py<PyAny>>,
    accept: Option<Py<PyAny>>,
}

impl PyCallbacks {
    fn clone_ref(&self, py: Python<'_>) -> PyCallbacks {
        let c = |o: &Option<Py<PyAny>>| o.as_ref().map(|o| o.clone_ref(py));
        PyCallbacks {
            atom_comp: c(&self.atom_comp),
            bond_comp: c(&self.bond_comp),
            progress: c(&self.progress),
            final_check: c(&self.final_check),
            accept: c(&self.accept),
        }
    }
}

/// The Python object of a molecule handed to the search, so callbacks see the
/// user's own `Mol` instead of a copy.
#[derive(Default)]
struct MolTable {
    entries: Vec<(usize, Py<PyAny>)>,
}

impl MolTable {
    fn py_mol(&self, py: Python<'_>, mol: &ROMol) -> PyResult<Py<PyAny>> {
        let addr = mol as *const ROMol as usize;
        if let Some((_, o)) = self.entries.iter().find(|(a, _)| *a == addr) {
            return Ok(o.clone_ref(py));
        }
        // a molecule the search had to copy (no ring information)
        Ok(Py::new(
            py,
            PyMol {
                mol: ROMol::copy_from(mol, false, -1),
            },
        )?
        .into_any())
    }
}

/// Carries a Python exception out of a callback: the search itself only knows
/// `RDKitError`.
type ErrorSlot = Arc<Mutex<Option<PyErr>>>;

fn callback_error(slot: &ErrorSlot, e: PyErr) -> RDKitError {
    *slot.lock().unwrap_or_else(|p| p.into_inner()) = Some(e);
    RDKitError::Runtime("exception in a Python callback".into())
}

/// Everything the Python-calling closures share.
#[derive(Clone)]
struct CallbackContext {
    mols: Arc<MolTable>,
    callbacks: Arc<PyCallbacks>,
    error: ErrorSlot,
    /// snapshot of the parameters handed to callbacks that take them
    base: Arc<MCSParameters>,
}

impl CallbackContext {
    fn call<R>(&self, f: impl FnOnce(Python<'_>) -> PyResult<R>) -> CoreResult<R> {
        Python::attach(f).map_err(|e| callback_error(&self.error, e))
    }

    fn py_parameters(
        &self,
        py: Python<'_>,
        params: &MCSParameters,
    ) -> PyResult<Py<PyMCSParameters>> {
        let mut p = PyMCSParameters::from_core(py, params)?;
        p.cb = self.callbacks.clone_ref(py);
        Py::new(py, p)
    }
}

fn atom_idx_match_tuple<'py>(
    py: Python<'py>,
    c1: &[u32],
    c2: &[u32],
    query: &Graph,
    target: &Graph,
) -> PyResult<Bound<'py, PyTuple>> {
    let pairs: Vec<(i64, i64)> = (0..query.num_vertices())
        .map(|i| (i64::from(query.atom(c1[i])), i64::from(target.atom(c2[i]))))
        .collect();
    PyTuple::new(py, pairs)
}

fn bond_idx_match_tuple<'py>(
    py: Python<'py>,
    c1: &[u32],
    c2: &[u32],
    mol1: &ROMol,
    query: &Graph,
    mol2: &ROMol,
    target: &Graph,
) -> PyResult<Bound<'py, PyTuple>> {
    let mut pairs: Vec<(i64, i64)> = Vec::with_capacity(query.num_edges());
    for &(s, t, _) in query.edges() {
        let qb = mol1
            .bond_idx_between(query.atom(c1[s as usize]), query.atom(c1[t as usize]))
            .ok_or_else(|| PyValueError::new_err("query bond not found"))?;
        let tb = mol2
            .bond_idx_between(target.atom(c2[s as usize]), target.atom(c2[t as usize]))
            .ok_or_else(|| PyValueError::new_err("target bond not found"))?;
        pairs.push((i64::from(qb), i64::from(tb)));
    }
    PyTuple::new(py, pairs)
}

fn pairs_tuple<'py>(py: Python<'py>, pairs: &[(i32, i32)]) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, pairs.iter().copied())
}

// ------------------------------------------------------------------ MCSParameters

/// Parameters controlling how the MCS is constructed
#[pyclass(
    name = "MCSParameters",
    module = "rdoxide.Chem.rdFMCS",
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PyMCSParameters {
    base: MCSParameters,
    atom_cp: Py<PyAtomCompareParameters>,
    bond_cp: Py<PyBondCompareParameters>,
    cb: PyCallbacks,
}

impl PyMCSParameters {
    fn from_core(py: Python<'_>, params: &MCSParameters) -> PyResult<Self> {
        Ok(PyMCSParameters {
            base: params.clone(),
            atom_cp: Py::new(
                py,
                PyAtomCompareParameters {
                    p: params.atom_compare_parameters.clone(),
                },
            )?,
            bond_cp: Py::new(
                py,
                PyBondCompareParameters {
                    p: params.bond_compare_parameters.clone(),
                },
            )?,
            cb: PyCallbacks::default(),
        })
    }

    /// The parameters of the search, with closures calling the configured
    /// Python objects.
    fn to_core(&self, py: Python<'_>, mols: Arc<MolTable>, error: ErrorSlot) -> MCSParameters {
        let mut p = self.base.clone();
        p.atom_compare_parameters = self.atom_cp.borrow(py).p.clone();
        p.bond_compare_parameters = self.bond_cp.borrow(py).p.clone();
        let ctx = CallbackContext {
            mols,
            callbacks: Arc::new(self.cb.clone_ref(py)),
            error,
            base: Arc::new(p.clone()),
        };
        if ctx.callbacks.atom_comp.is_some() {
            let c = ctx.clone();
            p.atom_typer = AtomTyper::Custom(Arc::new(move |params, mol1, a1, mol2, a2| {
                c.call(|py| {
                    let obj = c.callbacks.atom_comp.as_ref().expect("atom comparator");
                    let ps = Py::new(py, PyAtomCompareParameters { p: params.clone() })?;
                    let m1 = c.mols.py_mol(py, mol1)?;
                    let m2 = c.mols.py_mol(py, mol2)?;
                    obj.bind(py)
                        .call_method1("__call__", (ps, m1, a1, m2, a2))?
                        .is_truthy()
                })
            }));
        }
        if ctx.callbacks.bond_comp.is_some() {
            let c = ctx.clone();
            p.bond_typer = BondTyper::Custom(Arc::new(move |params, mol1, b1, mol2, b2| {
                c.call(|py| {
                    let obj = c.callbacks.bond_comp.as_ref().expect("bond comparator");
                    let ps = Py::new(py, PyBondCompareParameters { p: params.clone() })?;
                    let m1 = c.mols.py_mol(py, mol1)?;
                    let m2 = c.mols.py_mol(py, mol2)?;
                    obj.bind(py)
                        .call_method1("__call__", (ps, m1, b1, m2, b2))?
                        .is_truthy()
                })
            }));
        }
        if ctx.callbacks.progress.is_some() {
            let c = ctx.clone();
            p.progress_callback = Some(Arc::new(move |stat, params| {
                c.call(|py| {
                    let obj = c.callbacks.progress.as_ref().expect("progress callback");
                    let pd = Py::new(py, PyProgressData { d: *stat })?;
                    let ps = c.py_parameters(py, params)?;
                    obj.bind(py).call_method1("__call__", (pd, ps))?.is_truthy()
                })
            }));
        }
        if ctx.callbacks.final_check.is_some() {
            let c = ctx.clone();
            p.final_match_checker = Some(Arc::new(
                move |c1, c2, mol1, query, mol2, target, params| {
                    c.call(|py| {
                        let obj = c
                            .callbacks
                            .final_check
                            .as_ref()
                            .expect("final match checker");
                        let ps = c.py_parameters(py, params)?;
                        let atoms = atom_idx_match_tuple(py, c1, c2, query, target)?;
                        let bonds = bond_idx_match_tuple(py, c1, c2, mol1, query, mol2, target)?;
                        let m1 = c.mols.py_mol(py, mol1)?;
                        let m2 = c.mols.py_mol(py, mol2)?;
                        obj.bind(py)
                            .call_method1("__call__", (m1, m2, atoms, bonds, ps))?
                            .is_truthy()
                    })
                },
            ));
        }
        if ctx.callbacks.accept.is_some() {
            let c = ctx.clone();
            p.should_accept_mcs = Some(Arc::new(
                move |query, target, atom_pairs, bond_pairs, params| {
                    c.call(|py| {
                        let obj = c.callbacks.accept.as_ref().expect("acceptance callback");
                        let ps = c.py_parameters(py, params)?;
                        let atoms = pairs_tuple(py, atom_pairs)?;
                        let bonds = pairs_tuple(py, bond_pairs)?;
                        let q = c.mols.py_mol(py, query)?;
                        let t = c.mols.py_mol(py, target)?;
                        obj.bind(py)
                            .call_method1("__call__", (q, t, atoms, bonds, ps))?
                            .is_truthy()
                    })
                },
            ));
        }
        let _ = &ctx.base;
        p
    }
}

/// The class of `obj` must define a callable `__call__` itself (upstream looks
/// in the type's own dictionary).
fn require_own_call(obj: &Bound<'_, PyAny>, message: &str) -> PyResult<()> {
    let dict = obj.get_type().getattr("__dict__")?;
    match dict.call_method1("get", ("__call__",)) {
        Ok(f) if !f.is_none() && f.is_callable() => Ok(()),
        _ => Err(PyTypeError::new_err(message.to_string())),
    }
}

fn is_enum_instance(py: Python<'_>, obj: &Bound<'_, PyAny>, name: &str) -> PyResult<Option<i64>> {
    // enum members are ints of the registered enum class
    let cls = to_py(py, name, 0)?;
    if obj.is_instance(&cls.bind(py).get_type())? {
        return Ok(Some(obj.extract::<i64>()?));
    }
    Ok(None)
}

macro_rules! scalar_props {
    ($( ($py:ident, $get:ident, $set:ident, $field:ident, $ty:ty, $doc:literal) ),* $(,)?) => {
        #[pymethods]
        #[allow(non_snake_case)]
        impl PyMCSParameters {
            $(
                #[doc = $doc]
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.base.$field.clone()
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.base.$field = value;
                }
            )*
        }
    };
}

scalar_props!(
    (
        MaximizeBonds,
        get_maximize_bonds,
        set_maximize_bonds,
        maximize_bonds,
        bool,
        "toggles maximizing the number of bonds (instead of the number of atoms)"
    ),
    (
        Threshold,
        get_threshold,
        set_threshold,
        threshold,
        f64,
        "fraction of the dataset that must contain the MCS"
    ),
    (
        Timeout,
        get_timeout,
        set_timeout,
        timeout,
        u32,
        "timeout (in seconds) for the calculation"
    ),
    (
        Verbose,
        get_verbose,
        set_verbose,
        verbose,
        bool,
        "toggles verbose mode"
    ),
    (
        InitialSeed,
        get_initial_seed,
        set_initial_seed,
        initial_seed,
        String,
        "SMILES string to be used as the seed of the MCS"
    ),
    (
        StoreAll,
        get_store_all,
        set_store_all,
        store_all,
        bool,
        "toggles storage of degenerate MCSs"
    ),
);

#[pymethods]
#[allow(non_snake_case)]
impl PyMCSParameters {
    #[new]
    fn new(py: Python<'_>) -> PyResult<Self> {
        Self::from_core(py, &MCSParameters::default())
    }

    /// parameters for comparing atoms
    #[getter(AtomCompareParameters)]
    fn get_atom_compare_parameters(&self, py: Python<'_>) -> Py<PyAtomCompareParameters> {
        self.atom_cp.clone_ref(py)
    }
    #[setter(AtomCompareParameters)]
    fn set_atom_compare_parameters(
        &mut self,
        py: Python<'_>,
        value: PyRef<'_, PyAtomCompareParameters>,
    ) {
        self.atom_cp.borrow_mut(py).p = value.p.clone();
    }

    /// parameters for comparing bonds
    #[getter(BondCompareParameters)]
    fn get_bond_compare_parameters(&self, py: Python<'_>) -> Py<PyBondCompareParameters> {
        self.bond_cp.clone_ref(py)
    }
    #[setter(BondCompareParameters)]
    fn set_bond_compare_parameters(
        &mut self,
        py: Python<'_>,
        value: PyRef<'_, PyBondCompareParameters>,
    ) {
        self.bond_cp.borrow_mut(py).p = value.p.clone();
    }

    /// atom typer to be used. Must be one of the
    /// members of the rdFMCS.AtomCompare class or
    /// an instance of a user-defined subclass of
    /// rdFMCS.MCSAtomCompare
    #[getter(AtomTyper)]
    fn get_atom_typer(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        if let Some(obj) = &self.cb.atom_comp {
            return Ok(obj.clone_ref(py));
        }
        let c = match self.base.atom_typer {
            AtomTyper::Any => AtomComparator::Any,
            AtomTyper::Elements => AtomComparator::Elements,
            AtomTyper::Isotopes => AtomComparator::Isotopes,
            AtomTyper::AnyHeavyAtom => AtomComparator::AnyHeavyAtom,
            AtomTyper::Custom(_) => return Err(PyTypeError::new_err("Unknown AtomTyper")),
        };
        to_py(py, ATOM_COMPARE, atom_comparator_value(c))
    }
    #[setter(AtomTyper)]
    fn set_atom_typer(&mut self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Some(v) = is_enum_instance(py, value, ATOM_COMPARE)? {
            self.base.set_atom_typer_from_enum(atom_comparator_from(v)?);
            self.cb.atom_comp = None;
            return Ok(());
        }
        require_own_call(
            value,
            "AtomTyper must be an AtomCompare enum value or an instance of a rdFMCS.MCSAtomCompare subclass with a callable __call__() method",
        )?;
        self.cb.atom_comp = Some(value.clone().unbind());
        Ok(())
    }

    /// bond typer to be used. Must be one of the
    /// members of the rdFMCS.BondCompare class or
    /// an instance of a user-defined subclass of
    /// rdFMCS.MCSBondCompare
    #[getter(BondTyper)]
    fn get_bond_typer(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        if let Some(obj) = &self.cb.bond_comp {
            return Ok(obj.clone_ref(py));
        }
        let c = match self.base.bond_typer {
            BondTyper::Any => BondComparator::Any,
            BondTyper::Order => BondComparator::Order,
            BondTyper::OrderExact => BondComparator::OrderExact,
            BondTyper::Custom(_) => return Err(PyTypeError::new_err("Unknown BondTyper")),
        };
        to_py(py, BOND_COMPARE, c as i64)
    }
    #[setter(BondTyper)]
    fn set_bond_typer(&mut self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Some(v) = is_enum_instance(py, value, BOND_COMPARE)? {
            self.base.set_bond_typer_from_enum(bond_comparator_from(v)?);
            self.cb.bond_comp = None;
            return Ok(());
        }
        require_own_call(
            value,
            "BondTyper must be a BondCompare enum value or an instance of a rdFMCS.MCSBondCompare subclass with a callable __call__() method",
        )?;
        self.cb.bond_comp = Some(value.clone().unbind());
        Ok(())
    }

    /// progress callback class. Must be a
    /// user-defined subclass of rdFMCS.MCSProgress
    #[getter(ProgressCallback)]
    fn get_progress_callback(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.cb.progress.as_ref().map(|o| o.clone_ref(py))
    }
    #[setter(ProgressCallback)]
    fn set_progress_callback(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        require_own_call(
            value,
            "The __call__() method must be overridden in the rdFMCS.MCSProgress subclass",
        )?;
        self.cb.progress = Some(value.clone().unbind());
        Ok(())
    }

    /// seed final match checker callback class. Must be a
    /// user-defined subclass of rdFMCS.MCSFinalMatchCheck
    #[getter(FinalMatchChecker)]
    fn get_final_match_checker(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.cb.final_check.as_ref().map(|o| o.clone_ref(py))
    }
    #[setter(FinalMatchChecker)]
    fn set_final_match_checker(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        require_own_call(
            value,
            "The __call__() method must be overridden in the rdFMCS.MCSFinalMatchCheck subclass",
        )?;
        self.cb.final_check = Some(value.clone().unbind());
        Ok(())
    }

    /// MCS acceptance callback class. Must be a
    /// user-defined subclass of rdFMCS.MCSAcceptance
    #[getter(ShouldAcceptMCS)]
    fn get_should_accept_mcs(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.cb.accept.as_ref().map(|o| o.clone_ref(py))
    }
    #[setter(ShouldAcceptMCS)]
    fn set_should_accept_mcs(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        require_own_call(
            value,
            "The __call__() method must be overridden in the rdFMCS.MCSAcceptance subclass",
        )?;
        self.cb.accept = Some(value.clone().unbind());
        Ok(())
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        set_declared_attribute(slf.as_any(), name, value)
    }
}

// ------------------------------------------------------------------ results and FindMCS

/// used to return MCS results
#[pyclass(name = "MCSResult", module = "rdoxide.Chem.rdFMCS")]
pub struct PyMCSResult {
    r: MCSResult,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyMCSResult {
    /// number of atoms in MCS
    #[getter(numAtoms)]
    fn num_atoms(&self) -> u32 {
        self.r.num_atoms
    }

    /// number of bonds in MCS
    #[getter(numBonds)]
    fn num_bonds(&self) -> u32 {
        self.r.num_bonds
    }

    /// query molecule for the MCS
    #[getter(queryMol)]
    fn query_mol(&self) -> Option<PyMol> {
        self.r.query_mol.as_ref().map(|m| PyMol {
            mol: ROMol::copy_from(m, false, -1),
        })
    }

    /// SMARTS string for the MCS
    #[getter(smartsString)]
    fn smarts_string(&self) -> String {
        self.r.smarts_string.clone()
    }

    /// if True, the MCS calculation did not finish
    #[getter(canceled)]
    fn canceled(&self) -> bool {
        self.r.canceled
    }

    /// Dictionary collecting all degenerate (SMARTS, queryMol) pairs
    /// (empty if MCSParameters.StoreAll is False)
    #[getter(degenerateSmartsQueryMolDict)]
    fn degenerate_smarts_query_mol_dict<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for (smarts, mol) in &self.r.degenerate_smarts_query_mol_dict {
            d.set_item(
                smarts,
                PyMol {
                    mol: ROMol::copy_from(mol, false, -1),
                },
            )?;
        }
        Ok(d)
    }
}

/// The molecule objects of a list and their borrows.
type BorrowedMols<'py> = (Vec<Bound<'py, PyMol>>, Vec<PyRef<'py, PyMol>>);

/// The molecules of a list, borrowed for the duration of a search.
fn borrow_mols<'py>(mols: &Bound<'py, PyAny>) -> PyResult<BorrowedMols<'py>> {
    let mut objects = Vec::new();
    let mut refs = Vec::new();
    for item in mols.try_iter()? {
        let item = item?;
        if item.is_none() {
            return Err(PyValueError::new_err("molecule is None"));
        }
        let b = item.cast_into::<PyMol>()?;
        refs.push(b.borrow());
        objects.push(b);
    }
    Ok((objects, refs))
}

fn run_find_mcs(
    py: Python<'_>,
    mols: &Bound<'_, PyAny>,
    build: impl FnOnce(Python<'_>, Arc<MolTable>, ErrorSlot) -> PyResult<MCSParameters>,
) -> PyResult<PyMCSResult> {
    let (objects, refs) = borrow_mols(mols)?;
    let table = Arc::new(MolTable {
        entries: refs
            .iter()
            .zip(&objects)
            .map(|(r, o)| {
                (
                    &r.mol as *const ROMol as usize,
                    o.clone().into_any().unbind(),
                )
            })
            .collect(),
    });
    let error: ErrorSlot = Arc::new(Mutex::new(None));
    let params = build(py, table, Arc::clone(&error))?;
    let core: Vec<&ROMol> = refs.iter().map(|r| &r.mol).collect();
    let res = py.detach(|| fmcs::find_mcs(&core, Some(&params)));
    match res {
        Ok(r) => Ok(PyMCSResult { r }),
        Err(e) => match error.lock().unwrap_or_else(|p| p.into_inner()).take() {
            Some(py_err) => Err(py_err),
            None => Err(crate::errors::to_py(e)),
        },
    }
}

/// Find the MCS for a set of molecules
///
/// FindMCS(mols, maximizeBonds=True, threshold=1.0, timeout=3600,
///         verbose=False, matchValences=False, ringMatchesRingOnly=False,
///         completeRingsOnly=False, matchChiralTag=False,
///         atomCompare=AtomCompare.CompareElements,
///         bondCompare=BondCompare.CompareOrder,
///         ringCompare=RingCompare.IgnoreRingFusion, seedSmarts="")
/// FindMCS(mols, parameters)
#[pyfunction]
#[pyo3(name = "FindMCS", signature = (*args, **kwargs))]
fn py_find_mcs(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyMCSResult> {
    let a = Args {
        name: "FindMCS",
        args,
        kwargs,
    };
    let mols = a
        .get(0, "mols")?
        .ok_or_else(|| PyTypeError::new_err("FindMCS() missing argument 'mols'"))?;
    if let Some(params) = a.get(1, "parameters")?
        && let Ok(params) = params.cast::<PyMCSParameters>()
    {
        let params = params.borrow();
        return run_find_mcs(py, &mols, |py, table, error| {
            Ok(params.to_core(py, table, error))
        });
    }
    let enum_arg = |pos: usize, key: &str, default: i64| -> PyResult<i64> {
        a.get(pos, key)?.map_or(Ok(default), |v| v.extract())
    };
    let mut p = MCSParameters {
        maximize_bonds: a.flag(1, "maximizeBonds", true)?,
        threshold: a.get(2, "threshold")?.map_or(Ok(1.0), |v| v.extract())?,
        timeout: a.get(3, "timeout")?.map_or(Ok(3600), |v| v.extract())?,
        verbose: a.flag(4, "verbose", false)?,
        initial_seed: a
            .get(12, "seedSmarts")?
            .map_or(Ok(String::new()), |v| v.extract())?,
        ..Default::default()
    };
    p.atom_compare_parameters.match_valences = a.flag(5, "matchValences", false)?;
    let ring_matches_ring_only = a.flag(6, "ringMatchesRingOnly", false)?;
    let complete_rings_only = a.flag(7, "completeRingsOnly", false)?;
    p.atom_compare_parameters.match_chiral_tag = a.flag(8, "matchChiralTag", false)?;
    p.atom_compare_parameters.ring_matches_ring_only = ring_matches_ring_only;
    p.set_atom_typer_from_enum(atom_comparator_from(enum_arg(9, "atomCompare", 1)?)?);
    p.set_bond_typer_from_enum(bond_comparator_from(enum_arg(10, "bondCompare", 1)?)?);
    let ring_comp = ring_comparator_from(enum_arg(11, "ringCompare", 0)?)?;
    p.bond_compare_parameters.ring_matches_ring_only = ring_matches_ring_only;
    p.bond_compare_parameters.complete_rings_only = complete_rings_only;
    p.bond_compare_parameters.match_fused_rings = ring_comp != RingComparator::IgnoreRingFusion;
    p.bond_compare_parameters.match_fused_rings_strict =
        ring_comp == RingComparator::StrictRingFusion;
    run_find_mcs(py, &mols, |_, _, _| Ok(p))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing a C++ implementation of the FMCS algorithm",
    )?;
    add_enum(
        m,
        ATOM_COMPARE,
        &[
            ("CompareAny", 0),
            ("CompareElements", 1),
            ("CompareIsotopes", 2),
            ("CompareAnyHeavyAtom", 3),
        ],
        None,
    )?;
    add_enum(
        m,
        BOND_COMPARE,
        &[
            ("CompareAny", 0),
            ("CompareOrder", 1),
            ("CompareOrderExact", 2),
        ],
        None,
    )?;
    add_enum(
        m,
        RING_COMPARE,
        &[
            ("IgnoreRingFusion", 0),
            ("PermissiveRingFusion", 1),
            ("StrictRingFusion", 2),
        ],
        None,
    )?;
    m.add_class::<PyMCSResult>()?;
    m.add_class::<PyMCSParameters>()?;
    m.add_class::<PyAtomCompareParameters>()?;
    m.add_class::<PyBondCompareParameters>()?;
    m.add_class::<PyProgressData>()?;
    m.add_class::<PyAtomCompare>()?;
    m.add_class::<PyBondCompare>()?;
    m.add_class::<PyProgress>()?;
    m.add_class::<PyFinalMatchCheck>()?;
    m.add_class::<PyAcceptance>()?;
    m.add_function(wrap_pyfunction!(py_find_mcs, m)?)?;
    Ok(())
}
