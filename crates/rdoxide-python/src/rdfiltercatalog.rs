//! `rdkit.Chem.rdfiltercatalog` (upstream `FilterCatalog/nbWrap/FilterCatalog.cpp`):
//! structural-alert filters, with matchers that may be written in Python.

use std::any::Any;
use std::sync::{Arc, Mutex};

use pyo3::exceptions::{PyIndexError, PyNotImplementedError, PyRuntimeError, PyTypeError};
use pyo3::prelude::*;
use pyo3::pyclass_init::PyClassInitializer;
use pyo3::types::{PyBytes, PyDict, PyList, PyString, PyTuple};
use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::filter_catalog::catalogs;
use rdoxide::graphmol::filter_catalog::{
    And, ExclusionList, FilterCatalog, FilterCatalogEntry, FilterCatalogParams,
    FilterHierarchyMatcher, FilterMatch, FilterMatcher, HierarchyNode, Not, Or, SharedMatcher,
    SmartsMatcher, filter_catalog_can_serialize, functional_group_hierarchy, run_filter_catalog,
};
use rdoxide::graphmol::functional_group_hierarchy::flattened_functional_group_hierarchy;
use rdoxide::{RDKitError, Result as CoreResult};

use crate::enums::add_enum;
use crate::errors::to_py;
use crate::mol::PyMol;

const MODULE: &str = "rdkit.Chem.rdfiltercatalog";

/// A Python exception raised inside a Python-defined matcher; the matcher
/// interface only carries `RDKitError`, so the exception waits here until
/// the Python-facing function that started the match returns it.
static PENDING: Mutex<Option<PyErr>> = Mutex::new(None);

fn stash(e: PyErr) -> RDKitError {
    *PENDING.lock().unwrap_or_else(|p| p.into_inner()) = Some(e);
    RDKitError::Runtime("exception in a Python filter matcher".into())
}

fn map_err(e: RDKitError) -> PyErr {
    match PENDING.lock().unwrap_or_else(|p| p.into_inner()).take() {
        Some(pe) => pe,
        None => to_py(e),
    }
}

// ------------------------------------------------------------------ Python-defined matchers

/// A matcher implemented by a Python subclass of `FilterMatcher`.
struct PyMatcher {
    obj: Py<PyAny>,
}

impl std::fmt::Debug for PyMatcher {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("PyMatcher")
    }
}

fn py_mol(py: Python<'_>, mol: &ROMol) -> PyResult<Py<PyAny>> {
    Ok(Py::new(
        py,
        PyMol {
            mol: ROMol::copy_from(mol, false, -1),
        },
    )?
    .into_any())
}

impl FilterMatcher for PyMatcher {
    fn is_valid(&self) -> bool {
        Python::attach(|py| {
            self.obj
                .bind(py)
                .call_method0("IsValid")
                .and_then(|r| r.is_truthy())
                .unwrap_or_else(|e| {
                    stash(e);
                    false
                })
        })
    }

    fn name(&self) -> String {
        Python::attach(|py| {
            self.obj
                .bind(py)
                .call_method0("GetName")
                .and_then(|r| r.extract::<String>())
                .unwrap_or_else(|e| {
                    stash(e);
                    String::new()
                })
        })
    }

    fn get_matches(&self, mol: &ROMol, matches: &mut Vec<FilterMatch>) -> CoreResult<bool> {
        Python::attach(|py| {
            let mut run = || -> PyResult<bool> {
                let list = PyList::empty(py);
                let matched = self
                    .obj
                    .bind(py)
                    .call_method1("GetMatches", (py_mol(py, mol)?, &list))?
                    .is_truthy()?;
                if matched {
                    for item in list.iter() {
                        matches.push(item.extract::<PyRef<'_, PyFilterMatch>>()?.m.clone());
                    }
                }
                Ok(matched)
            };
            run().map_err(stash)
        })
    }

    fn has_match(&self, mol: &ROMol) -> CoreResult<bool> {
        Python::attach(|py| {
            let run = || -> PyResult<bool> {
                self.obj
                    .bind(py)
                    .call_method1("HasMatch", (py_mol(py, mol)?,))?
                    .is_truthy()
            };
            run().map_err(stash)
        })
    }

    fn copy(&self) -> SharedMatcher {
        Python::attach(|py| {
            Arc::new(PyMatcher {
                obj: self.obj.clone_ref(py),
            })
        })
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

// ------------------------------------------------------------------ matcher classes

/// Base class for matching molecules to filters.
///
///  A FilterMatcherBase supplies the following API
///  - IsValid() returns True if the matcher is valid for use, False otherwise
///  - HasMatch(mol) returns True if the molecule matches the filter
///  - GetMatches(mol) -> [FilterMatch, FilterMatch] returns all the FilterMatch data
///        that matches the molecule
///
/// print( FilterMatcherBase ) will print user-friendly information about the filter
/// Note that a FilterMatcherBase can be combined from many FilterMatcherBases
/// This is why GetMatches can return multiple FilterMatcherBases.
#[pyclass(
    name = "FilterMatcher",
    module = "rdkit.Chem.rdfiltercatalog",
    subclass
)]
pub struct PyFilterMatcher {
    name: String,
}

/// The native matcher behind a Python object, or `None` for a plain
/// (possibly Python-subclassed) `FilterMatcher`.
fn native_matcher(obj: &Bound<'_, PyAny>) -> Option<SharedMatcher> {
    if let Ok(m) = obj.cast::<PySmartsMatcher>() {
        return Some(Arc::new(m.borrow().m.clone()));
    }
    if let Ok(m) = obj.cast::<PyExclusionList>() {
        return Some(Arc::new(m.borrow().m.clone()));
    }
    if let Ok(m) = obj.cast::<PyFilterHierarchyMatcher>() {
        let node = m
            .borrow()
            .node
            .read()
            .unwrap_or_else(|e| e.into_inner())
            .clone();
        return Some(Arc::new(node));
    }
    if let Ok(m) = obj.cast::<PyAnd>() {
        return Some(Arc::new(m.borrow().m.clone()));
    }
    if let Ok(m) = obj.cast::<PyOr>() {
        return Some(Arc::new(m.borrow().m.clone()));
    }
    if let Ok(m) = obj.cast::<PyNot>() {
        return Some(Arc::new(m.borrow().m.clone()));
    }
    None
}

/// Any matcher object as a native matcher (a copy for the native classes, a
/// callback wrapper for Python subclasses).
fn to_matcher(obj: &Bound<'_, PyAny>) -> PyResult<SharedMatcher> {
    if let Some(m) = native_matcher(obj) {
        return Ok(m);
    }
    if obj.is_instance_of::<PyFilterMatcher>() {
        return Ok(Arc::new(PyMatcher {
            obj: obj.clone().unbind(),
        }));
    }
    Err(PyTypeError::new_err("expected a FilterMatcher"))
}

fn named(name: &str) -> PyFilterMatcher {
    PyFilterMatcher {
        name: name.to_string(),
    }
}

/// The Python object of a native matcher.
fn matcher_to_py(py: Python<'_>, m: &SharedMatcher) -> PyResult<Py<PyAny>> {
    let any = m.as_any();
    macro_rules! wrap {
        ($core:ty, $py:ident, $make:expr) => {
            if let Some(c) = any.downcast_ref::<$core>() {
                let sub = $make(c);
                return Ok(Py::new(
                    py,
                    PyClassInitializer::from(named(&m.name())).add_subclass(sub),
                )?
                .into_any());
            }
        };
    }
    wrap!(SmartsMatcher, PySmartsMatcher, |c: &SmartsMatcher| {
        PySmartsMatcher { m: c.clone() }
    });
    wrap!(ExclusionList, PyExclusionList, |c: &ExclusionList| {
        PyExclusionList { m: c.clone() }
    });
    wrap!(And, PyAnd, |c: &And| PyAnd { m: c.clone() });
    wrap!(Or, PyOr, |c: &Or| PyOr { m: c.clone() });
    wrap!(Not, PyNot, |c: &Not| PyNot { m: c.clone() });
    wrap!(
        FilterHierarchyMatcher,
        PyFilterHierarchyMatcher,
        |c: &FilterHierarchyMatcher| {
            PyFilterHierarchyMatcher {
                node: Arc::new(std::sync::RwLock::new(c.clone())),
            }
        }
    );
    if let Some(p) = any.downcast_ref::<PyMatcher>() {
        return Ok(p.obj.clone_ref(py));
    }
    Err(PyTypeError::new_err("unknown matcher type"))
}

fn matches_to_list<'py>(
    py: Python<'py>,
    matches: Vec<FilterMatch>,
) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty(py);
    for m in matches {
        list.append(Py::new(py, PyFilterMatch { m })?)?;
    }
    Ok(list)
}

fn not_overridden(what: &str) -> PyErr {
    PyNotImplementedError::new_err(format!(
        "{what} must be implemented in a subclass of FilterMatcher"
    ))
}

#[pymethods]
#[allow(non_snake_case)]
impl PyFilterMatcher {
    #[new]
    #[pyo3(signature = (*args, **_kwargs))]
    fn new(args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
        // a Python subclass may take any arguments of its own
        let name = args
            .get_item(0)
            .ok()
            .and_then(|a| a.extract::<String>().ok())
            .unwrap_or_default();
        PyFilterMatcher { name }
    }

    /// Initialize the matcher with its name (for Python subclasses).
    #[pyo3(signature = (*args, **_kwargs))]
    fn __init__(&mut self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) {
        if let Some(n) = args
            .get_item(0)
            .ok()
            .and_then(|n| n.extract::<String>().ok())
        {
            self.name = n;
        }
    }

    /// Return True if the filter matcher is valid, False otherwise
    fn IsValid(slf: &Bound<'_, Self>) -> PyResult<bool> {
        match native_matcher(slf.as_any()) {
            Some(m) => Ok(m.is_valid()),
            None => Err(not_overridden("IsValid")),
        }
    }

    /// Returns True if mol matches the filter
    fn HasMatch(slf: &Bound<'_, Self>, mol: &Bound<'_, PyMol>) -> PyResult<bool> {
        match native_matcher(slf.as_any()) {
            Some(m) => m.has_match(&mol.borrow().mol).map_err(map_err),
            None => Err(not_overridden("HasMatch")),
        }
    }

    /// Returns the list of matching subfilters mol matches any filter
    fn GetMatches<'py>(
        slf: &Bound<'py, Self>,
        mol: &Bound<'py, PyMol>,
    ) -> PyResult<Bound<'py, PyList>> {
        let py = slf.py();
        let Some(m) = native_matcher(slf.as_any()) else {
            return Err(not_overridden("GetMatches"));
        };
        let mut matches = Vec::new();
        if m.get_matches(&mol.borrow().mol, &mut matches)
            .map_err(map_err)?
        {
            return matches_to_list(py, matches);
        }
        Ok(PyList::empty(py))
    }

    fn GetName(slf: &Bound<'_, Self>) -> String {
        match native_matcher(slf.as_any()) {
            Some(m) => m.name(),
            None => slf.borrow().name.clone(),
        }
    }

    fn __str__(slf: &Bound<'_, Self>) -> String {
        Self::GetName(slf)
    }
}

/// Smarts Matcher Filter
///  basic constructors:
///    SmartsMatcher( name, smarts_pattern, minCount=1, maxCount=UINT_MAX )
///    SmartsMatcher( name, molecule, minCount=1, maxCount=UINT_MAX )
///
///   note: If the supplied smarts pattern is not valid, the IsValid() function will
///    return False
#[pyclass(name = "SmartsMatcher", module = "rdkit.Chem.rdfiltercatalog", extends = PyFilterMatcher)]
pub struct PySmartsMatcher {
    m: SmartsMatcher,
}

#[pymethods]
#[allow(non_snake_case)]
impl PySmartsMatcher {
    #[new]
    #[pyo3(signature = (*args))]
    fn new(args: &Bound<'_, PyTuple>) -> PyResult<PyClassInitializer<Self>> {
        let min = |pos: usize| -> PyResult<u32> {
            args.get_item(pos).ok().map_or(Ok(1), |v| v.extract())
        };
        let max = |pos: usize| -> PyResult<u32> {
            args.get_item(pos)
                .ok()
                .map_or(Ok(u32::MAX), |v| v.extract())
        };
        let (m, name) = match args.len() {
            0 => {
                return Err(PyTypeError::new_err(
                    "SmartsMatcher() needs a name or a molecule",
                ));
            }
            1 => {
                let first = args.get_item(0)?;
                if let Ok(mol) = first.cast::<PyMol>() {
                    let pattern = Arc::new(ROMol::copy_from(&mol.borrow().mol, false, -1));
                    (
                        SmartsMatcher::from_pattern(
                            rdoxide::graphmol::filter_catalog::SMARTS_MATCH_NAME_DEFAULT,
                            pattern,
                            1,
                            u32::MAX,
                        ),
                        rdoxide::graphmol::filter_catalog::SMARTS_MATCH_NAME_DEFAULT.to_string(),
                    )
                } else {
                    let name: String = first.extract()?;
                    (SmartsMatcher::with_name(&name), name)
                }
            }
            _ => {
                let name: String = args.get_item(0)?.extract()?;
                let second = args.get_item(1)?;
                let m = if let Ok(mol) = second.cast::<PyMol>() {
                    let pattern = Arc::new(ROMol::copy_from(&mol.borrow().mol, false, -1));
                    SmartsMatcher::from_pattern(&name, pattern, min(2)?, max(3)?)
                } else {
                    let smarts: String = second.extract()?;
                    SmartsMatcher::from_smarts(&name, &smarts, min(2)?, max(3)?)
                };
                (m, name)
            }
        };
        Ok(PyClassInitializer::from(named(&name)).add_subclass(PySmartsMatcher { m }))
    }

    /// Set the pattern molecule (or the smarts pattern) for the SmartsMatcher
    fn SetPattern(&mut self, pat: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Ok(mol) = pat.cast::<PyMol>() {
            self.m
                .set_pattern(Arc::new(ROMol::copy_from(&mol.borrow().mol, false, -1)));
        } else {
            self.m.set_pattern_smarts(&pat.extract::<String>()?);
        }
        Ok(())
    }

    fn GetPattern(&self) -> Option<PyMol> {
        self.m.pattern().map(|p| PyMol {
            mol: ROMol::copy_from(p, false, -1),
        })
    }

    /// Get the minimum times pattern must appear for the filter to match
    fn GetMinCount(&self) -> u32 {
        self.m.min_count()
    }

    /// Set the minimum times pattern must appear to match
    fn SetMinCount(&mut self, count: u32) {
        self.m.set_min_count(count);
    }

    /// Get the maximum times pattern can appear for the filter to match
    fn GetMaxCount(&self) -> u32 {
        self.m.max_count()
    }

    /// Set the maximum times pattern can appear for the filter to match
    fn SetMaxCount(&mut self, count: u32) {
        self.m.set_max_count(count);
    }
}

#[pyclass(name = "ExclusionList", module = "rdkit.Chem.rdfiltercatalog", extends = PyFilterMatcher)]
pub struct PyExclusionList {
    m: ExclusionList,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyExclusionList {
    #[new]
    fn new() -> PyClassInitializer<Self> {
        PyClassInitializer::from(named("Not any of")).add_subclass(PyExclusionList {
            m: ExclusionList::new(),
        })
    }

    /// Set a list of FilterMatcherBases that should not appear in a molecule
    fn SetExclusionPatterns(slf: &Bound<'_, Self>, list: &Bound<'_, PyAny>) -> PyResult<()> {
        // the list may contain the matcher itself, so convert before borrowing mutably
        let mut patterns = Vec::new();
        for item in list.try_iter()? {
            patterns.push(to_matcher(&item?)?.copy());
        }
        slf.borrow_mut().m.set_exclusion_patterns(patterns);
        Ok(())
    }

    /// Add a FilterMatcherBase that should not appear in a molecule
    fn AddPattern(slf: &Bound<'_, Self>, base: &Bound<'_, PyAny>) -> PyResult<()> {
        let m = to_matcher(base)?;
        slf.borrow_mut().m.add_pattern(m.as_ref()).map_err(map_err)
    }
}

macro_rules! binary_op {
    ($class:ident, $pyname:literal, $core:ident, $default:literal) => {
        #[pyclass(name = $pyname, module = "rdkit.Chem.rdfiltercatalog.FilterMatchOps", extends = PyFilterMatcher)]
        pub struct $class {
            m: $core,
        }

        #[pymethods]
        impl $class {
            #[new]
            fn new(
                arg1: &Bound<'_, PyAny>,
                arg2: &Bound<'_, PyAny>,
            ) -> PyResult<PyClassInitializer<Self>> {
                let m = $core::new(to_matcher(arg1)?.copy(), to_matcher(arg2)?.copy());
                Ok(PyClassInitializer::from(named($default)).add_subclass($class { m }))
            }
        }
    };
}

binary_op!(PyAnd, "And", And, "And");
binary_op!(PyOr, "Or", Or, "Or");

#[pyclass(name = "Not", module = "rdkit.Chem.rdfiltercatalog.FilterMatchOps", extends = PyFilterMatcher)]
pub struct PyNot {
    m: Not,
}

#[pymethods]
impl PyNot {
    #[new]
    fn new(arg1: &Bound<'_, PyAny>) -> PyResult<PyClassInitializer<Self>> {
        Ok(PyClassInitializer::from(named("Not")).add_subclass(PyNot {
            m: Not::new(to_matcher(arg1)?.copy()),
        }))
    }
}

/// Hierarchical Filter
///  basic constructors:
///    FilterHierarchyMatcher( matcher )
///    where can be any FilterMatcherBase (SmartsMatcher, etc)
///  FilterHierarchyMatcher's have children and can form matching
///   trees.  When GetFilterMatches is called, the most specific (
///   i.e. lowest node in a branch) is returned.
#[pyclass(name = "FilterHierarchyMatcher", module = "rdkit.Chem.rdfiltercatalog", extends = PyFilterMatcher)]
pub struct PyFilterHierarchyMatcher {
    node: HierarchyNode,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyFilterHierarchyMatcher {
    #[new]
    #[pyo3(signature = (matcher = None))]
    fn new(matcher: Option<&Bound<'_, PyAny>>) -> PyResult<PyClassInitializer<Self>> {
        let node = match matcher {
            Some(m) => FilterHierarchyMatcher::with_matcher(to_matcher(m)?.as_ref()),
            None => FilterHierarchyMatcher::new(),
        };
        Ok(
            PyClassInitializer::from(named("Unnamed FilterMatcherBase")).add_subclass(
                PyFilterHierarchyMatcher {
                    node: Arc::new(std::sync::RwLock::new(node)),
                },
            ),
        )
    }

    /// Set the filtermatcher pattern for this node. An empty node is
    /// considered a root node and passes along the matches to the children.
    fn SetPattern(&self, matcher: &Bound<'_, PyAny>) -> PyResult<()> {
        let m = to_matcher(matcher)?;
        self.node
            .write()
            .unwrap_or_else(|e| e.into_inner())
            .set_pattern(m.as_ref())
            .map_err(map_err)
    }

    /// Add a child node to this hierarchy.
    fn AddChild(
        &self,
        py: Python<'_>,
        hierarchy: &Bound<'_, PyFilterHierarchyMatcher>,
    ) -> PyResult<Py<PyAny>> {
        let child = hierarchy
            .borrow()
            .node
            .read()
            .unwrap_or_else(|e| e.into_inner())
            .clone();
        let added = self
            .node
            .write()
            .unwrap_or_else(|e| e.into_inner())
            .add_child(&child)
            .map_err(map_err)?;
        Ok(Py::new(
            py,
            PyClassInitializer::from(named("Unnamed FilterMatcherBase"))
                .add_subclass(PyFilterHierarchyMatcher { node: added }),
        )?
        .into_any())
    }
}

/// Object that holds the result of running FilterMatcherBase::GetMatches
///
///  - filterMatch holds the FilterMatchBase that triggered the match
///  - atomPairs holds the [ (query_atom_idx, target_atom_idx) ] pairs for the matches.
///
/// Note that some matches may not have atom pairs (especially matches that
/// use FilterMatchOps.Not)
#[pyclass(
    name = "FilterMatch",
    module = "rdkit.Chem.rdfiltercatalog",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyFilterMatch {
    m: FilterMatch,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyFilterMatch {
    #[new]
    #[pyo3(signature = (filter, atomPairs))]
    fn new(filter: &Bound<'_, PyAny>, atomPairs: Vec<(i32, i32)>) -> PyResult<Self> {
        Ok(PyFilterMatch {
            m: FilterMatch {
                filter_match: to_matcher(filter)?,
                atom_pairs: atomPairs,
            },
        })
    }

    #[getter(filterMatch)]
    fn filter_match(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        matcher_to_py(py, &self.m.filter_match)
    }

    #[getter(atomPairs)]
    fn atom_pairs(&self) -> Vec<(i32, i32)> {
        self.m.atom_pairs.clone()
    }
}

// ------------------------------------------------------------------ entries and catalogs

/// FilterCatalogEntry
/// A filter catalog entry is an entry in a filter catalog.
/// Each filter is named and is used to flag a molecule usually for some
/// undesirable property.
#[pyclass(
    name = "FilterCatalogEntry",
    module = "rdkit.Chem.rdfiltercatalog",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyFilterCatalogEntry {
    e: Arc<FilterCatalogEntry>,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyFilterCatalogEntry {
    #[new]
    #[pyo3(signature = (name = None, matcher = None))]
    fn new(name: Option<&str>, matcher: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let e = match (name, matcher) {
            (Some(name), Some(m)) => {
                FilterCatalogEntry::with_matcher(name, Some(to_matcher(m)?.copy()))
            }
            (None, None) => FilterCatalogEntry::new(),
            _ => return Err(PyTypeError::new_err("FilterCatalogEntry(name, matcher)")),
        };
        Ok(PyFilterCatalogEntry { e: Arc::new(e) })
    }

    fn IsValid(&self) -> bool {
        self.e.is_valid()
    }

    /// Get the description of the catalog entry
    fn GetDescription(&self) -> String {
        self.e.description()
    }

    /// Set the description of the catalog entry
    fn SetDescription(&self, description: &str) {
        self.e.set_description(description);
    }

    /// Retrieve the list of filters that match the molecule
    fn GetFilterMatches<'py>(
        &self,
        py: Python<'py>,
        mol: &Bound<'py, PyMol>,
    ) -> PyResult<Bound<'py, PyList>> {
        let mut matches = Vec::new();
        if self
            .e
            .filter_matches(&mol.borrow().mol, &mut matches)
            .map_err(map_err)?
        {
            return matches_to_list(py, matches);
        }
        Ok(PyList::empty(py))
    }

    /// Returns True if the catalog entry contains filters that match the molecule
    fn HasFilterMatch(&self, mol: &Bound<'_, PyMol>) -> PyResult<bool> {
        self.e.has_filter_match(&mol.borrow().mol).map_err(map_err)
    }

    fn Serialize(&self) -> PyResult<()> {
        Err(PyRuntimeError::new_err(
            "Boost SERIALIZATION is not enabled",
        ))
    }

    fn GetPropList(&self) -> Vec<String> {
        self.e.prop_list()
    }

    fn SetProp(&self, key: &str, val: &str) {
        self.e.set_prop(key, val);
    }

    fn GetProp(&self, key: &str) -> PyResult<String> {
        self.e.prop(key).map_err(map_err)
    }

    fn ClearProp(&self, key: &str) {
        self.e.clear_prop(key);
    }
}

/// Parameters of a FilterCatalog
#[pyclass(
    name = "FilterCatalogParams",
    module = "rdkit.Chem.rdfiltercatalog",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyFilterCatalogParams {
    p: FilterCatalogParams,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyFilterCatalogParams {
    #[new]
    #[pyo3(signature = (catalogs = None))]
    fn new(catalogs: Option<u32>) -> Self {
        PyFilterCatalogParams {
            p: catalogs.map_or_else(FilterCatalogParams::new, FilterCatalogParams::with_catalogs),
        }
    }

    fn AddCatalog(&mut self, catalogs: u32) -> bool {
        self.p.add_catalog(catalogs)
    }
}

#[pyclass(name = "FilterCatalog", module = "rdkit.Chem.rdfiltercatalog")]
pub struct PyFilterCatalog {
    c: FilterCatalog,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyFilterCatalog {
    #[new]
    #[pyo3(signature = (arg = None))]
    fn new(arg: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let c = match arg {
            None => FilterCatalog::new(),
            Some(a) => {
                if let Ok(p) = a.cast::<PyFilterCatalogParams>() {
                    FilterCatalog::from_params(&p.borrow().p).map_err(map_err)?
                } else if a.is_instance_of::<PyBytes>() || a.is_instance_of::<PyString>() {
                    return Err(PyRuntimeError::new_err(
                        "Boost SERIALIZATION is not enabled",
                    ));
                } else {
                    FilterCatalog::from_catalogs(a.extract::<u32>()?).map_err(map_err)?
                }
            }
        };
        Ok(PyFilterCatalog { c })
    }

    fn Serialize(&self) -> PyResult<()> {
        Err(PyRuntimeError::new_err(
            "Pickling of FilterCatalog instances is not enabled",
        ))
    }

    fn __getstate__(&self) -> PyResult<()> {
        self.Serialize()
    }

    /// Add a FilterCatalogEntry to the catalog
    fn AddEntry(&mut self, entry: &Bound<'_, PyFilterCatalogEntry>) {
        self.c
            .add_entry(Arc::new(entry.borrow().e.as_ref().clone()));
    }

    /// Remove the given entry (or the entry at an index) from the catalog
    fn RemoveEntry(&mut self, obj: &Bound<'_, PyAny>) -> PyResult<bool> {
        if let Ok(idx) = obj.extract::<usize>() {
            if obj.is_instance_of::<pyo3::types::PyInt>() {
                return Ok(self.c.remove_entry_at(idx));
            }
        }
        let entry = obj.cast::<PyFilterCatalogEntry>()?;
        Ok(self.c.remove_entry(&entry.borrow().e))
    }

    /// Returns the number of entries in the catalog
    fn GetNumEntries(&self) -> usize {
        self.c.num_entries()
    }

    /// Return the FilterCatalogEntry at the specified index
    fn GetEntryWithIdx(&self, idx: usize) -> PyResult<PyFilterCatalogEntry> {
        self.GetEntry(idx)
    }

    fn GetEntry(&self, idx: usize) -> PyResult<PyFilterCatalogEntry> {
        match self.c.entry(idx) {
            Ok(e) => Ok(PyFilterCatalogEntry { e: Arc::clone(e) }),
            Err(_) => Err(PyIndexError::new_err(idx.to_string())),
        }
    }

    /// Returns True if the catalog has an entry that matches mol
    fn HasMatch(&self, mol: &Bound<'_, PyMol>) -> PyResult<bool> {
        self.c.has_match(&mol.borrow().mol).map_err(map_err)
    }

    /// Return the first catalog entry that matches mol
    fn GetFirstMatch(&self, mol: &Bound<'_, PyMol>) -> PyResult<Option<PyFilterCatalogEntry>> {
        Ok(self
            .c
            .first_match(&mol.borrow().mol)
            .map_err(map_err)?
            .map(|e| PyFilterCatalogEntry { e }))
    }

    /// Return all catalog entries that match mol
    fn GetMatches(&self, mol: &Bound<'_, PyMol>) -> PyResult<Vec<PyFilterCatalogEntry>> {
        Ok(self
            .c
            .matches(&mol.borrow().mol)
            .map_err(map_err)?
            .into_iter()
            .map(|e| PyFilterCatalogEntry { e })
            .collect())
    }

    /// Return every matching filter from all catalog entries that match mol
    fn GetFilterMatches<'py>(
        &self,
        py: Python<'py>,
        mol: &Bound<'py, PyMol>,
    ) -> PyResult<Bound<'py, PyList>> {
        let matches = self.c.filter_matches(&mol.borrow().mol).map_err(map_err)?;
        matches_to_list(py, matches)
    }
}

/// Returns True if the FilterCatalog is serializable (requires boost serialization)
#[pyfunction]
#[pyo3(name = "FilterCatalogCanSerialize")]
fn py_can_serialize() -> bool {
    filter_catalog_can_serialize()
}

/// Run the filter catalog on the input list of smiles strings.
/// Use numThreads=0 to use all available processors.
/// Returns a vector of vectors. For each input smiles, a vector of
/// FilterCatalogEntry objects are returned for each matched filter. If a
/// molecule matches no filter, the vector will be empty. If a smiles string
/// can't be parsed, a 'Bad smiles' entry is returned.
#[pyfunction]
#[pyo3(name = "RunFilterCatalog", signature = (filterCatalog, smiles, numThreads = 1))]
#[allow(non_snake_case)]
fn py_run_filter_catalog(
    py: Python<'_>,
    filterCatalog: &Bound<'_, PyFilterCatalog>,
    smiles: Vec<String>,
    numThreads: i32,
) -> PyResult<Vec<Vec<PyFilterCatalogEntry>>> {
    let catalog = &filterCatalog.borrow().c;
    let results = py.detach(|| run_filter_catalog(catalog, &smiles, numThreads));
    results
        .into_iter()
        .map(|r| {
            r.map(|v| v.into_iter().map(|e| PyFilterCatalogEntry { e }).collect())
                .map_err(map_err)
        })
        .collect()
}

/// Returns the functional group hierarchy filter catalog
#[pyfunction]
#[pyo3(name = "GetFunctionalGroupHierarchy")]
fn py_functional_group_hierarchy() -> PyFilterCatalog {
    PyFilterCatalog {
        c: functional_group_hierarchy().clone(),
    }
}

/// Returns the flattened functional group hierarchy as a dictionary
///  of name:ROMOL_SPTR substructure items
#[pyfunction]
#[pyo3(name = "GetFlattenedFunctionalGroupHierarchy", signature = (normalized = false))]
fn py_flattened_hierarchy<'py>(py: Python<'py>, normalized: bool) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    for (name, mol) in flattened_functional_group_hierarchy(normalized) {
        d.set_item(
            name,
            PyMol {
                mol: ROMol::copy_from(mol, false, -1),
            },
        )?;
    }
    Ok(d)
}

/// Create an integer pair (tuple) for use in MatchTypeVect
#[pyfunction]
#[pyo3(name = "IntPair")]
fn py_int_pair(a: i32, b: i32) -> (i32, i32) {
    (a, b)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add("MatchTypeVect", py.import("builtins")?.getattr("list")?)?;
    add_enum(
        m,
        "FilterCatalogs",
        &[
            ("PAINS_A", i64::from(catalogs::PAINS_A)),
            ("PAINS_B", i64::from(catalogs::PAINS_B)),
            ("PAINS_C", i64::from(catalogs::PAINS_C)),
            ("PAINS", i64::from(catalogs::PAINS)),
            ("BRENK", i64::from(catalogs::BRENK)),
            ("NIH", i64::from(catalogs::NIH)),
            ("ZINC", i64::from(catalogs::ZINC)),
            ("CHEMBL_Glaxo", i64::from(catalogs::CHEMBL_GLAXO)),
            ("CHEMBL_Dundee", i64::from(catalogs::CHEMBL_DUNDEE)),
            ("CHEMBL_BMS", i64::from(catalogs::CHEMBL_BMS)),
            ("CHEMBL_SureChEMBL", i64::from(catalogs::CHEMBL_SURECHEMBL)),
            ("CHEMBL_MLSMR", i64::from(catalogs::CHEMBL_MLSMR)),
            (
                "CHEMBL_Inpharmatica",
                i64::from(catalogs::CHEMBL_INPHARMATICA),
            ),
            ("CHEMBL_LINT", i64::from(catalogs::CHEMBL_LINT)),
            ("CHEMBL", i64::from(catalogs::CHEMBL)),
            ("ALL", i64::from(catalogs::ALL)),
        ],
        None,
    )?;
    m.add_class::<PyFilterCatalogParams>()?;
    // FilterCatalogParams.FilterCatalogs, as in upstream
    let params_cls = m.getattr("FilterCatalogParams")?;
    params_cls.setattr("FilterCatalogs", m.getattr("FilterCatalogs")?)?;
    m.delattr("FilterCatalogs")?;
    m.getattr("__all__")?
        .call_method1("remove", ("FilterCatalogs",))?;
    m.add_class::<PyFilterMatch>()?;
    m.add_class::<PyFilterMatcher>()?;
    m.add("FilterMatcherBase", m.getattr("FilterMatcher")?)?;
    m.add("PythonFilterMatcher", m.getattr("FilterMatcher")?)?;
    m.add_class::<PySmartsMatcher>()?;
    m.add_class::<PyExclusionList>()?;
    m.add_class::<PyFilterHierarchyMatcher>()?;
    m.add_class::<PyFilterCatalogEntry>()?;
    m.add_class::<PyFilterCatalog>()?;
    m.add_function(wrap_pyfunction!(py_can_serialize, m)?)?;
    m.add_function(wrap_pyfunction!(py_run_filter_catalog, m)?)?;
    m.add_function(wrap_pyfunction!(py_functional_group_hierarchy, m)?)?;
    m.add_function(wrap_pyfunction!(py_flattened_hierarchy, m)?)?;
    m.add_function(wrap_pyfunction!(py_int_pair, m)?)?;

    let ops = PyModule::new(py, "FilterMatchOps")?;
    ops.setattr("__name__", format!("{MODULE}.FilterMatchOps"))?;
    ops.add_class::<PyAnd>()?;
    ops.add_class::<PyOr>()?;
    ops.add_class::<PyNot>()?;
    m.add("FilterMatchOps", &ops)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(format!("{MODULE}.FilterMatchOps"), &ops)?;
    Ok(())
}
