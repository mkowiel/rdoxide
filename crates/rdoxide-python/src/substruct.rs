//! Substructure matching on `Mol` (`HasSubstructMatch`, `GetSubstructMatch`,
//! `GetSubstructMatches`) and `SubstructMatchParameters` — port of the
//! corresponding parts of upstream `Code/GraphMol/Wrap/Mol.cpp` and
//! `substructmethods.h`.

use std::sync::{Arc, Mutex};

use pyo3::exceptions::{PyRuntimeError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::substruct::{
    AtomCoordsMatcher, FinalMatchCheck, MatchCheck, MatchVect, SubstructMatchParameters,
    substruct_match,
};
use rdoxide::graphmol::substruct_exists::substruct_has_match;

use crate::atom::{AtomHandle, PyAtom};
use crate::bond::{BondHandle, PyBond};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::molbundle::PyMolBundle;

/// Parameters controlling substructure matching
#[pyclass(
    name = "SubstructMatchParameters",
    module = "rdoxide.Chem.rdchem",
    skip_from_py_object
)]
#[derive(Clone)]
#[allow(non_snake_case)]
pub struct PySubstructMatchParameters {
    /// Use chirality in determining whether or not atoms/bonds match
    #[pyo3(get, set)]
    useChirality: bool,
    /// take enhanced stereochemistry into account while doing the match. This only has an effect if useChirality is also True.
    #[pyo3(get, set)]
    useEnhancedStereo: bool,
    /// aromatic and conjugated bonds match each other
    #[pyo3(get, set)]
    aromaticMatchesConjugated: bool,
    /// aromatic and single or double bonds match each other
    #[pyo3(get, set)]
    aromaticMatchesSingleOrDouble: bool,
    /// use generic groups (=homology groups) as a post-filtering step (if any are present in the molecule)
    #[pyo3(get, set)]
    useGenericMatchers: bool,
    /// Consider query-query matches, not just simple matches
    #[pyo3(get, set)]
    useQueryQueryMatches: bool,
    /// Allow recursive queries
    #[pyo3(get, set)]
    recursionPossible: bool,
    /// uniquify (by atom index) match results
    #[pyo3(get, set)]
    uniquify: bool,
    /// maximum number of matches to return
    #[pyo3(get, set)]
    maxMatches: u32,
    /// maximum number of recursive matches to find
    #[pyo3(get, set)]
    maxRecursiveMatches: u32,
    /// number of threads to use when multi-threading is possible.0 selects the number of concurrent threads supported by thehardware. negative values are added to the number of concurrentthreads supported by the hardware.
    #[pyo3(get, set)]
    numThreads: i32,
    /// bond properties that must be equivalent in order to match.
    #[pyo3(get, set)]
    bondProperties: Vec<String>,
    /// atom properties that must be equivalent in order to match.
    #[pyo3(get, set)]
    atomProperties: Vec<String>,
    /// If set, query atoms and bonds with specified stereochemistry will match atoms and bonds with unspecified stereochemistry.
    #[pyo3(get, set)]
    specifiedStereoQueryMatchesUnspecified: bool,
    /// if set, only the extraAtomCheck will be used to determine whether or not atoms match
    #[pyo3(get, set)]
    extraAtomCheckOverridesDefaultCheck: bool,
    /// if set, only the extraBondCheck will be used to determine whether or not bonds match
    #[pyo3(get, set)]
    extraBondCheckOverridesDefaultCheck: bool,
    extra_atom_check: Option<Arc<Callback>>,
    extra_bond_check: Option<Arc<Callback>>,
    extra_final_check: Option<Arc<Callback>>,
}

#[pymethods]
impl PySubstructMatchParameters {
    #[new]
    fn py_new() -> Self {
        Self::from_core(&SubstructMatchParameters::default())
    }

    /// allows you to provide a function that will be called
    /// for each atom pair that matches during substructure searching,
    /// after all other comparisons have passed.
    /// The function should return true or false indicating whether or not
    /// that atom-match should be accepted. An AtomCoordsMatcher can be
    /// given instead of a function.
    #[pyo3(name = "setExtraAtomCheckFunc", signature = (func))]
    fn set_extra_atom_check_func(&mut self, func: &Bound<'_, PyAny>) {
        self.extra_atom_check = Some(Arc::new(Callback::new(func)));
    }

    /// allows you to provide a function that will be called
    /// for each bond pair that matches during substructure searching,
    /// after all other comparisons have passed.
    /// The function should return true or false indicating whether or not
    /// that bond-match should be accepted.
    #[pyo3(name = "setExtraBondCheckFunc", signature = (func))]
    fn set_extra_bond_check_func(&mut self, func: &Bound<'_, PyAny>) {
        self.extra_bond_check = Some(Arc::new(Callback::new(func)));
    }

    /// allows you to provide a function that will be called to validate
    /// a final match. The function is called with the molecule and the
    /// list of matched atom indices and should return true or false
    /// indicating whether or not that match should be accepted.
    #[pyo3(name = "setExtraFinalCheck", signature = (func))]
    fn set_extra_final_check(&mut self, func: &Bound<'_, PyAny>) {
        self.extra_final_check = Some(Arc::new(Callback::new(func)));
    }
}

/// A Python object registered as an extra check: a function, or an
/// `AtomCoordsMatcher` (evaluated without the GIL).
#[derive(Debug)]
enum Callback {
    Function(Py<PyAny>),
    Coords(Py<PyAtomCoordsMatcher>),
}

impl Callback {
    fn new(func: &Bound<'_, PyAny>) -> Self {
        match func.cast::<PyAtomCoordsMatcher>() {
            Ok(m) => Callback::Coords(m.clone().unbind()),
            Err(_) => Callback::Function(func.clone().unbind()),
        }
    }
}

/// The first error raised by a Python callback; matching stops accepting
/// pairs once it is set and the error is re-raised after the search.
type ErrorSlot = Arc<Mutex<Option<PyErr>>>;

/// Calls `f` with the GIL held, remembering the first error and answering
/// `false` for failed calls.
fn guarded(errors: &ErrorSlot, f: impl FnOnce(Python<'_>) -> PyResult<bool>) -> bool {
    Python::attach(|py| match f(py) {
        Ok(accepted) => accepted,
        Err(e) => {
            errors
                .lock()
                .unwrap_or_else(|p| p.into_inner())
                .get_or_insert(e);
            false
        }
    })
}

/// Allows using atom coordinates as part of substructure matching
#[pyclass(
    name = "AtomCoordsMatcher",
    module = "rdoxide.Chem.rdchem",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct PyAtomCoordsMatcher {
    /// reference conformer ID
    #[pyo3(get, set)]
    refConfId: i32,
    /// query conformer ID
    #[pyo3(get, set)]
    queryConfId: i32,
    /// squared distance tolerance
    #[pyo3(get, set)]
    tol2: f64,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyAtomCoordsMatcher {
    /// constructor taking reference and query conformer IDs and a distance tolerance
    #[new]
    #[pyo3(signature = (refConfId = -1, queryConfId = -1, tol = 1e-4))]
    fn py_new(refConfId: i32, queryConfId: i32, tol: f64) -> Self {
        PyAtomCoordsMatcher {
            refConfId,
            queryConfId,
            tol2: tol * tol,
        }
    }

    fn __call__(&self, py: Python<'_>, queryAtom: &PyAtom, targetAtom: &PyAtom) -> PyResult<bool> {
        let (qm, qi) = atom_location(queryAtom)?;
        let (tm, ti) = atom_location(targetAtom)?;
        let (qm, tm) = (qm.borrow(py), tm.borrow(py));
        Ok(self.core().matches(&qm.mol, qi, &tm.mol, ti))
    }
}

impl PyAtomCoordsMatcher {
    fn core(&self) -> AtomCoordsMatcher {
        AtomCoordsMatcher {
            ref_conf_id: self.refConfId,
            query_conf_id: self.queryConfId,
            tol2: self.tol2,
        }
    }
}

fn atom_location(atom: &PyAtom) -> PyResult<(&Py<PyMol>, u32)> {
    match &atom.handle {
        AtomHandle::InMol { mol, idx } => Ok((mol, *idx)),
        AtomHandle::Owned(_) => Err(PyRuntimeError::new_err("atom has no owning molecule")),
    }
}

impl PySubstructMatchParameters {
    /// The Python object for the given core parameters (callbacks unset).
    pub(crate) fn from_core(d: &SubstructMatchParameters) -> Self {
        PySubstructMatchParameters {
            useChirality: d.use_chirality,
            useEnhancedStereo: d.use_enhanced_stereo,
            aromaticMatchesConjugated: d.aromatic_matches_conjugated,
            aromaticMatchesSingleOrDouble: d.aromatic_matches_single_or_double,
            useGenericMatchers: d.use_generic_matchers,
            useQueryQueryMatches: d.use_query_query_matches,
            recursionPossible: d.recursion_possible,
            uniquify: d.uniquify,
            maxMatches: d.max_matches,
            maxRecursiveMatches: d.max_recursive_matches,
            numThreads: d.num_threads,
            bondProperties: d.bond_properties.clone(),
            atomProperties: d.atom_properties.clone(),
            specifiedStereoQueryMatchesUnspecified: d.specified_stereo_query_matches_unspecified,
            extraAtomCheckOverridesDefaultCheck: d.extra_atom_check_overrides_default,
            extraBondCheckOverridesDefaultCheck: d.extra_bond_check_overrides_default,
            extra_atom_check: None,
            extra_bond_check: None,
            extra_final_check: None,
        }
    }

    /// The core parameters for matching `query` against `mol`; the Python
    /// callbacks see atoms and bonds of those two molecules.
    fn to_core(
        &self,
        py: Python<'_>,
        mol: &Bound<'_, PyMol>,
        query: &Bound<'_, PyMol>,
        errors: &ErrorSlot,
    ) -> SubstructMatchParameters {
        let pair = Pair {
            mol: mol.clone().unbind(),
            query: query.clone().unbind(),
            errors: errors.clone(),
        };
        SubstructMatchParameters {
            extra_atom_check: self
                .extra_atom_check
                .as_deref()
                .map(|cb| pair.pair_check(py, cb, PairKind::Atom)),
            extra_atom_check_overrides_default: self.extraAtomCheckOverridesDefaultCheck,
            extra_bond_check: self
                .extra_bond_check
                .as_deref()
                .map(|cb| pair.pair_check(py, cb, PairKind::Bond)),
            extra_bond_check_overrides_default: self.extraBondCheckOverridesDefaultCheck,
            extra_final_check: self
                .extra_final_check
                .as_deref()
                .map(|cb| pair.final_check(py, cb)),
            ..self.base_params()
        }
    }

    /// The parameters without Python callbacks.
    pub(crate) fn base_params(&self) -> SubstructMatchParameters {
        SubstructMatchParameters {
            use_chirality: self.useChirality,
            use_enhanced_stereo: self.useEnhancedStereo,
            aromatic_matches_conjugated: self.aromaticMatchesConjugated,
            aromatic_matches_single_or_double: self.aromaticMatchesSingleOrDouble,
            use_generic_matchers: self.useGenericMatchers,
            use_query_query_matches: self.useQueryQueryMatches,
            recursion_possible: self.recursionPossible,
            uniquify: self.uniquify,
            max_matches: self.maxMatches,
            max_recursive_matches: self.maxRecursiveMatches,
            num_threads: self.numThreads,
            bond_properties: self.bondProperties.clone(),
            atom_properties: self.atomProperties.clone(),
            specified_stereo_query_matches_unspecified: self.specifiedStereoQueryMatchesUnspecified,
            ..SubstructMatchParameters::default()
        }
    }
}

/// The molecules a search runs on, as seen by Python callbacks.
struct Pair {
    mol: Py<PyMol>,
    query: Py<PyMol>,
    errors: ErrorSlot,
}

#[derive(Clone, Copy)]
enum PairKind {
    Atom,
    Bond,
}

impl Pair {
    /// `extraAtomCheck` / `extraBondCheck`: `func(queryItem, molItem)`.
    fn pair_check(&self, py: Python<'_>, cb: &Callback, kind: PairKind) -> MatchCheck {
        let func = match cb {
            Callback::Function(f) => f.clone_ref(py),
            Callback::Coords(m) => return m.borrow(py).core().into(),
        };
        let (mol, query, errors) = (
            self.mol.clone_ref(py),
            self.query.clone_ref(py),
            self.errors.clone(),
        );
        MatchCheck::new(move |_, qi, _, mi| {
            guarded(&errors, |py| {
                let item = |owner: &Py<PyMol>, idx: u32| -> PyResult<Py<PyAny>> {
                    let owner = owner.clone_ref(py);
                    Ok(match kind {
                        PairKind::Atom => Py::new(
                            py,
                            PyAtom {
                                handle: AtomHandle::InMol { mol: owner, idx },
                            },
                        )?
                        .into_any(),
                        PairKind::Bond => Py::new(
                            py,
                            PyBond {
                                handle: BondHandle::InMol { mol: owner, idx },
                            },
                        )?
                        .into_any(),
                    })
                };
                func.bind(py)
                    .call1((item(&query, qi)?, item(&mol, mi)?))?
                    .is_truthy()
            })
        })
    }

    /// `extraFinalCheck`: `func(mol, [matched atom indices])`.
    fn final_check(&self, py: Python<'_>, cb: &Callback) -> FinalMatchCheck {
        let Callback::Function(func) = cb else {
            return FinalMatchCheck::new(|_, _| true);
        };
        let (func, mol, errors) = (
            func.clone_ref(py),
            self.mol.clone_ref(py),
            self.errors.clone(),
        );
        FinalMatchCheck::new(move |_, aids| {
            guarded(&errors, |py| {
                func.bind(py)
                    .call1((mol.clone_ref(py), aids.to_vec()))?
                    .is_truthy()
            })
        })
    }
}

struct Args<'a, 'py> {
    name: &'static str,
    args: &'a Bound<'py, PyTuple>,
    kwargs: Option<&'a Bound<'py, PyDict>>,
}

impl<'py> Args<'_, 'py> {
    fn get(&self, pos: usize, key: &str) -> PyResult<Option<Bound<'py, PyAny>>> {
        let kw = match self.kwargs {
            Some(k) => k.get_item(key)?,
            None => None,
        };
        let p = if pos < self.args.len() {
            Some(self.args.get_item(pos)?)
        } else {
            None
        };
        match (p, kw) {
            (Some(_), Some(_)) => Err(PyTypeError::new_err(format!(
                "{}() got multiple values for argument '{key}'",
                self.name
            ))),
            (a, b) => Ok(a.or(b)),
        }
    }
    fn flag(&self, pos: usize, key: &str, default: bool) -> PyResult<bool> {
        self.get(pos, key)?.map_or(Ok(default), |v| v.is_truthy())
    }
    fn query(&self) -> PyResult<Molecules<'py>> {
        let q = self.get(0, "query")?.ok_or_else(|| {
            PyTypeError::new_err(format!("{}() missing required argument 'query'", self.name))
        })?;
        if let Ok(b) = q.cast::<PyMolBundle>() {
            return Ok(Molecules::Bundle(b.clone()));
        }
        Ok(Molecules::Mol(q.cast::<PyMol>()?.clone()))
    }
    fn params(&self) -> Option<Bound<'py, PySubstructMatchParameters>> {
        self.get(1, "params")
            .ok()
            .flatten()
            .and_then(|p| p.cast_into::<PySubstructMatchParameters>().ok())
    }
}

/// What a method needs to run a search: either the user's parameters
/// (possibly with Python callbacks) or plain core parameters.
#[derive(Clone)]
enum Params<'py> {
    User(Bound<'py, PySubstructMatchParameters>),
    Core(SubstructMatchParameters),
}

/// A molecule or a bundle of molecules, as the subject or query of a search.
pub(crate) enum Molecules<'py> {
    Mol(Bound<'py, PyMol>),
    Bundle(Bound<'py, PyMolBundle>),
}

impl<'py> Molecules<'py> {
    fn mols(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyMol>>> {
        match self {
            Molecules::Mol(m) => Ok(vec![m.clone()]),
            Molecules::Bundle(b) => b
                .borrow()
                .bundle
                .mols()
                .iter()
                .map(|m| {
                    Bound::new(
                        py,
                        PyMol {
                            mol: m.as_ref().clone(),
                        },
                    )
                })
                .collect(),
        }
    }
}

/// `SubstructMatch` for molecules and bundles: the first subject molecule and
/// query molecule combination that matches gives the result.
fn run_any(
    py: Python<'_>,
    subject: &Molecules<'_>,
    query: &Molecules<'_>,
    params: Params<'_>,
) -> PyResult<Vec<MatchVect>> {
    let queries = query.mols(py)?;
    for mol in subject.mols(py)? {
        for q in &queries {
            let res = run(py, &mol, q, params.clone())?;
            if !res.is_empty() {
                return Ok(res);
            }
        }
    }
    Ok(Vec::new())
}

fn run(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    query: &Bound<'_, PyMol>,
    params: Params<'_>,
) -> PyResult<Vec<MatchVect>> {
    run_with(py, mol, query, params, substruct_match)
}

/// Whether `query` matches `mol` (the existence-only search).
fn run_has(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    query: &Bound<'_, PyMol>,
    params: Params<'_>,
) -> PyResult<bool> {
    run_with(py, mol, query, params, substruct_has_match)
}

fn run_with<T: Send>(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    query: &Bound<'_, PyMol>,
    params: Params<'_>,
    search: impl Fn(&ROMol, &ROMol, &SubstructMatchParameters) -> rdoxide::Result<T> + Send + Sync,
) -> PyResult<T> {
    let errors = ErrorSlot::default();
    let ps = match params {
        Params::User(p) => p.borrow().to_core(py, mol, query, &errors),
        Params::Core(ps) => ps,
    };
    let result = {
        let (mol, query) = (mol.borrow(), query.borrow());
        let (m, q) = (&mol.mol, &query.mol);
        py.detach(|| search(m, q, &ps))
    };
    if let Some(e) = errors.lock().unwrap_or_else(|p| p.into_inner()).take() {
        return Err(e);
    }
    result.py()
}

fn run_any_has(
    py: Python<'_>,
    subject: &Molecules<'_>,
    query: &Molecules<'_>,
    params: Params<'_>,
) -> PyResult<bool> {
    let queries = query.mols(py)?;
    for mol in subject.mols(py)? {
        for q in &queries {
            if run_has(py, &mol, q, params.clone())? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn match_tuple<'py>(py: Python<'py>, m: &MatchVect) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, m.iter().map(|&(_, a)| a))
}

pub(crate) fn has_substruct_match_impl<'py>(
    subject: &Molecules<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<bool> {
    let py = args.py();
    let a = Args {
        name: "HasSubstructMatch",
        args,
        kwargs,
    };
    let query = a.query()?;
    let ps = match a.params() {
        Some(p) => Params::User(p),
        None => Params::Core(SubstructMatchParameters {
            recursion_possible: a.flag(1, "recursionPossible", true)?,
            use_chirality: a.flag(2, "useChirality", false)?,
            use_query_query_matches: a.flag(3, "useQueryQueryMatches", false)?,
            max_matches: 1,
            ..SubstructMatchParameters::default()
        }),
    };
    run_any_has(py, subject, &query, ps)
}

pub(crate) fn get_substruct_match_impl<'py>(
    subject: &Molecules<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let py = args.py();
    let a = Args {
        name: "GetSubstructMatch",
        args,
        kwargs,
    };
    let query = a.query()?;
    let ps = match a.params() {
        Some(p) => Params::User(p),
        None => Params::Core(SubstructMatchParameters {
            use_chirality: a.flag(1, "useChirality", false)?,
            use_query_query_matches: a.flag(2, "useQueryQueryMatches", false)?,
            max_matches: 1,
            ..SubstructMatchParameters::default()
        }),
    };
    let res = run_any(py, subject, &query, ps)?;
    match res.first() {
        Some(m) => match_tuple(py, m),
        None => Ok(PyTuple::empty(py)),
    }
}

pub(crate) fn get_substruct_matches_impl<'py>(
    subject: &Molecules<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let py = args.py();
    let a = Args {
        name: "GetSubstructMatches",
        args,
        kwargs,
    };
    let query = a.query()?;
    let ps = match a.params() {
        Some(p) => Params::User(p),
        None => Params::Core(SubstructMatchParameters {
            uniquify: a.flag(1, "uniquify", true)?,
            use_chirality: a.flag(2, "useChirality", false)?,
            use_query_query_matches: a.flag(3, "useQueryQueryMatches", false)?,
            max_matches: match a.get(4, "maxMatches")? {
                Some(v) => v.extract()?,
                None => 1000,
            },
            ..SubstructMatchParameters::default()
        }),
    };
    let res = run_any(py, subject, &query, ps)?;
    let items = res
        .iter()
        .map(|m| match_tuple(py, m))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, items)
}

#[pymethods]
impl PyMol {
    /// Queries whether or not the molecule contains a particular substructure.
    ///
    ///   ARGUMENTS:
    ///     - query: a Molecule
    ///
    ///     - recursionPossible: (optional)
    ///
    ///     - useChirality: enables the use of stereochemistry in the matching
    ///
    ///     - useQueryQueryMatches: use query-query matching logic
    ///
    ///   RETURNS: True or False
    ///
    /// Alternatively: HasSubstructMatch(query, params) with SubstructMatchParameters.
    #[pyo3(name = "HasSubstructMatch", signature = (*args, **kwargs))]
    fn has_substruct_match<'py>(
        slf: &Bound<'py, Self>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<bool> {
        has_substruct_match_impl(&Molecules::Mol(slf.clone()), args, kwargs)
    }

    /// Returns the indices of the molecule's atoms that match a substructure query.
    ///
    ///   ARGUMENTS:
    ///     - query: a Molecule
    ///
    ///     - useChirality: enables the use of stereochemistry in the matching
    ///
    ///     - useQueryQueryMatches: use query-query matching logic
    ///
    ///   RETURNS: a tuple of integers
    ///
    ///   NOTES:
    ///      - only a single match is returned
    ///      - the ordering of the indices corresponds to the atom ordering
    ///          in the query. For example, the first index is for the atom in
    ///          this molecule that matches the first atom in the query.
    ///
    /// Alternatively: GetSubstructMatch(query, params) with SubstructMatchParameters.
    #[pyo3(name = "GetSubstructMatch", signature = (*args, **kwargs))]
    fn get_substruct_match<'py>(
        slf: &Bound<'py, Self>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        get_substruct_match_impl(&Molecules::Mol(slf.clone()), args, kwargs)
    }

    /// Returns tuples of the indices of the molecule's atoms that match a substructure query.
    ///
    ///   ARGUMENTS:
    ///     - query: a Molecule.
    ///     - uniquify: (optional) determines whether or not the matches are uniquified.
    ///                 Defaults to 1.
    ///
    ///     - useChirality: enables the use of stereochemistry in the matching
    ///
    ///     - useQueryQueryMatches: use query-query matching logic
    ///
    ///     - maxMatches: The maximum number of matches that will be returned.
    ///                   In high-symmetry cases with medium-sized molecules, it is
    ///                   very easy to end up with a combinatorial explosion in the
    ///                   number of possible matches. This argument prevents that from
    ///                   having unintended consequences
    ///
    ///   RETURNS: a tuple of tuples of integers
    ///
    ///   NOTE:
    ///      - the ordering of the indices corresponds to the atom ordering
    ///          in the query. For example, the first index is for the atom in
    ///          this molecule that matches the first atom in the query.
    ///
    /// Alternatively: GetSubstructMatches(query, params) with SubstructMatchParameters.
    #[pyo3(name = "GetSubstructMatches", signature = (*args, **kwargs))]
    fn get_substruct_matches<'py>(
        slf: &Bound<'py, Self>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        get_substruct_matches_impl(&Molecules::Mol(slf.clone()), args, kwargs)
    }
}

/// This function is included to speed substructure queries from databases,
/// it's probably not of general interest.
///
///   ARGUMENTS:
///     - pkl: a Molecule pickle
///
///     - query: a Molecule
///
///     - recursionPossible: (optional)
///
///     - useChirality: (optional)
///
///     - useQueryQueryMatches: use query-query matching logic
///
///   RETURNS: True or False
#[pyfunction]
#[pyo3(name = "_HasSubstructMatchStr", signature = (pkl, query, recursionPossible = true, useChirality = false, useQueryQueryMatches = false))]
#[allow(non_snake_case)]
fn has_substruct_match_str(
    py: Python<'_>,
    #[pyo3(from_py_with = crate::pyconv::text_bytes)] pkl: Vec<u8>,
    query: PyRef<'_, PyMol>,
    #[pyo3(from_py_with = crate::pyconv::truthy)] recursionPossible: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useChirality: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] useQueryQueryMatches: bool,
) -> PyResult<bool> {
    let ps = SubstructMatchParameters {
        recursion_possible: recursionPossible,
        use_chirality: useChirality,
        use_query_query_matches: useQueryQueryMatches,
        max_matches: 1,
        ..SubstructMatchParameters::default()
    };
    let query = &query.mol;
    py.detach(|| {
        let mol = rdoxide::graphmol::pickle::mol_from_pickle(&pkl, 0)
            .map_err(|_| pyo3::exceptions::PyValueError::new_err("Null Molecule"))?;
        substruct_match(&mol, query, &ps)
            .map(|m| !m.is_empty())
            .map_err(crate::errors::to_py)
    })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(has_substruct_match_str, m)?)?;
    m.add_class::<PySubstructMatchParameters>()?;
    m.add_class::<PyAtomCoordsMatcher>()?;
    Ok(())
}
