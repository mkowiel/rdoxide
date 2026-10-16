//! `rdkit.Chem.rdSubstructLibrary` (upstream
//! `GraphMol/SubstructLibrary/Wrap/SubstructLibraryWrap.cpp`).
//!
//! Libraries pickle to this port's own (text) format.

use std::collections::HashMap;

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::pyclass_init::PyClassInitializer;
use pyo3::types::{PyBool, PyBytes, PyDict, PyTuple};
use rdoxide::graphmol::substruct::SubstructMatchParameters;
use rdoxide::graphmol::substruct_library::{
    KeyFromPropHolder, MolHolder, PatternHolder, Shared, SubstructLibrary, shared,
};

use crate::datastructs::{PyExplicitBitVect, bytes_arg};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::substruct::PySubstructMatchParameters;

fn read<T>(holder: &Shared<T>) -> std::sync::RwLockReadGuard<'_, T> {
    holder
        .read()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

fn write<T>(holder: &Shared<T>) -> std::sync::RwLockWriteGuard<'_, T> {
    holder
        .write()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

// ------------------------------------------------------------ molecule holders

/// Base class for holding molecules used in the Substructure Library.
/// Instantiations of this class are passed into the SubstructureLibrary.
/// The API is quite simple:
///   AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule
///   GetMol(idx) -> return the molecule at index idx
#[pyclass(
    name = "MolHolderBase",
    module = "rdkit.Chem.rdSubstructLibrary",
    subclass
)]
pub struct PyMolHolderBase {
    inner: Shared<MolHolder>,
}

#[pymethods]
impl PyMolHolderBase {
    fn __len__(&self) -> usize {
        read(&self.inner).len() as usize
    }

    /// Adds molecule to the molecule holder
    #[pyo3(name = "AddMol", signature = (m))]
    fn add_mol(&self, m: PyRef<'_, PyMol>) -> PyResult<u32> {
        write(&self.inner).add_mol(&m.mol).py()
    }

    /// Returns a particular molecule in the molecule holder
    ///
    ///   ARGUMENTS:
    ///     - idx: which molecule to return
    ///
    ///   NOTE: molecule indices start at 0
    #[pyo3(name = "GetMol", signature = (arg1))]
    fn get_mol(&self, arg1: u32) -> PyResult<Option<PyMol>> {
        molecule(&self.inner, arg1)
    }
}

fn molecule(holder: &Shared<MolHolder>, idx: u32) -> PyResult<Option<PyMol>> {
    Ok(read(holder)
        .get_mol(idx)
        .py()?
        .map(|m| PyMol { mol: (*m).clone() }))
}

macro_rules! mol_holder_class {
    ($py:ident, $name:literal, $variant:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyclass(name = $name, module = "rdkit.Chem.rdSubstructLibrary", extends = PyMolHolderBase)]
        pub struct $py;

        #[pymethods]
        impl $py {
            #[new]
            fn new() -> PyClassInitializer<Self> {
                PyClassInitializer::from(PyMolHolderBase {
                    inner: shared(MolHolder::$variant(Vec::new())),
                })
                .add_subclass($py)
            }
        }
    };
}

mol_holder_class!(
    PyMolHolder,
    "MolHolder",
    Mols,
    "Holds raw in-memory molecules\n  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule\n  GetMol(idx,sanitize=True) -> return the molecule at index idx\n"
);
mol_holder_class!(
    PyCachedMolHolder,
    "CachedMolHolder",
    Binary,
    "Holds molecules in their binary representation.\nThis allows more molecules to be held in memory at a time\n  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule\n\n  AddBinary(data) -> adds a picked molecule molecule to the molecule holder, returns index of molecule\n                     The data is stored as-is, no checking is done for validity.\n  GetMol(idx) -> return the molecule at index idx\n"
);
mol_holder_class!(
    PyCachedSmilesMolHolder,
    "CachedSmilesMolHolder",
    Smiles,
    "Holds molecules as smiles string\nThis allows more molecules to be held in memory at a time\n  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule\n\n  AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule\n                       The smiles is stored as-is, no checking is done for validity.\n  GetMol(idx) -> return the molecule at index idx\n"
);
mol_holder_class!(
    PyCachedTrustedSmilesMolHolder,
    "CachedTrustedSmilesMolHolder",
    TrustedSmiles,
    "Holds molecules as trusted smiles string\nThis allows more molecules to be held in memory at a time and avoids RDKit sanitization\noverhead.\nSee: http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html\n  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule\n\n  AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule\n                       The smiles is stored as-is, no checking is done for validity.\n  GetMol(idx,s) -> return the molecule at index idx, \n              note, only light sanitization is done here, for instance\n              the molecules RingInfo is not initialized\n"
);

#[pymethods]
impl PyCachedMolHolder {
    /// Add a pickled binary molecule to the holder, returns index of molecule
    ///
    /// The data is stored as-is, no checking is done for validity.
    #[pyo3(name = "AddBinary", signature = (pickle))]
    fn add_binary(slf: PyRef<'_, Self>, pickle: &Bound<'_, PyAny>) -> PyResult<u32> {
        let base = slf.as_super();
        let data = bytes_arg(pickle)?;
        write(&base.inner).add_binary(data).py()
    }
}

macro_rules! add_smiles {
    ($($py:ident),*) => {$(
        #[pymethods]
        impl $py {
            /// Add a smiles to the dataset, returns index of molecule
            ///
            /// The smiles is stored as-is, no checking is done for validity.
            #[pyo3(name = "AddSmiles", signature = (smiles))]
            fn add_smiles(slf: PyRef<'_, Self>, smiles: &str) -> PyResult<u32> {
                write(&slf.as_super().inner).add_smiles(smiles).py()
            }
        }
    )*};
}
add_smiles!(PyCachedSmilesMolHolder, PyCachedTrustedSmilesMolHolder);

/// The Python object (of the matching subclass) that shares `inner`.
fn wrap_mol_holder(py: Python<'_>, inner: Shared<MolHolder>) -> PyResult<Py<PyAny>> {
    let kind = match &*read(&inner) {
        MolHolder::Mols(_) => 0,
        MolHolder::Binary(_) => 1,
        MolHolder::Smiles(_) => 2,
        MolHolder::TrustedSmiles(_) => 3,
    };
    let base = PyClassInitializer::from(PyMolHolderBase { inner });
    Ok(match kind {
        0 => Py::new(py, base.add_subclass(PyMolHolder))?.into_any(),
        1 => Py::new(py, base.add_subclass(PyCachedMolHolder))?.into_any(),
        2 => Py::new(py, base.add_subclass(PyCachedSmilesMolHolder))?.into_any(),
        _ => Py::new(py, base.add_subclass(PyCachedTrustedSmilesMolHolder))?.into_any(),
    })
}

// ------------------------------------------------------------ pattern holders

/// Base class for the fingerprints used to rule out impossible matches
#[pyclass(
    name = "FPHolderBase",
    module = "rdkit.Chem.rdSubstructLibrary",
    subclass
)]
pub struct PyFPHolderBase {
    inner: Shared<PatternHolder>,
}

#[pymethods]
impl PyFPHolderBase {
    fn __len__(&self) -> usize {
        read(&self.inner).len() as usize
    }

    /// Adds a molecule to the fingerprinter
    #[pyo3(name = "AddMol", signature = (m))]
    fn add_mol(&self, m: PyRef<'_, PyMol>) -> PyResult<u32> {
        write(&self.inner).add_mol(&m.mol).py()
    }

    /// Adds a raw bit vector to the fingerprinter
    #[pyo3(name = "AddFingerprint", signature = (v))]
    fn add_fingerprint(&self, v: PyRef<'_, PyExplicitBitVect>) -> u32 {
        write(&self.inner).add_fingerprint(v.bv.clone())
    }

    /// Return the bit vector at the specified index
    #[pyo3(name = "GetFingerprint", signature = (idx))]
    fn get_fingerprint(&self, idx: u32) -> PyResult<PyExplicitBitVect> {
        Ok(PyExplicitBitVect::wrap(
            read(&self.inner).get_fingerprint(idx).py()?.clone(),
        ))
    }

    /// Returns True if the given query bit vector is contained in the fingerprint at the specified index
    #[pyo3(name = "PassesFilter", signature = (idx, query))]
    fn passes_filter(&self, idx: u32, query: PyRef<'_, PyExplicitBitVect>) -> PyResult<bool> {
        read(&self.inner).passes_filter(idx, &query.bv).py()
    }

    /// Make the search fingerprint for a query molecule
    #[pyo3(name = "MakeFingerprint", signature = (mol))]
    fn make_fingerprint(&self, mol: PyRef<'_, PyMol>) -> PyResult<PyExplicitBitVect> {
        Ok(PyExplicitBitVect::wrap(
            read(&self.inner).make_fingerprint(&mol.mol).py()?,
        ))
    }
}

/// Holds fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.
#[pyclass(name = "PatternHolder", module = "rdkit.Chem.rdSubstructLibrary", extends = PyFPHolderBase, subclass)]
pub struct PyPatternHolder;

#[pymethods]
impl PyPatternHolder {
    #[new]
    #[pyo3(signature = (numBits = PatternHolder::DEFAULT_NUM_BITS))]
    #[allow(non_snake_case)]
    fn new(numBits: u32) -> PyClassInitializer<Self> {
        PyClassInitializer::from(PyFPHolderBase {
            inner: shared(PatternHolder::new(numBits)),
        })
        .add_subclass(PyPatternHolder)
    }
}

/// Holds tautomeric fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.
/// These fingerprints are designed to be used with TautomerQueries.
#[pyclass(name = "TautomerPatternHolder", module = "rdkit.Chem.rdSubstructLibrary", extends = PyPatternHolder)]
pub struct PyTautomerPatternHolder;

#[pymethods]
impl PyTautomerPatternHolder {
    #[new]
    #[pyo3(signature = (numBits = PatternHolder::DEFAULT_NUM_BITS))]
    #[allow(non_snake_case)]
    fn new(numBits: u32) -> PyClassInitializer<Self> {
        PyClassInitializer::from(PyFPHolderBase {
            inner: shared(PatternHolder::new_tautomeric(numBits)),
        })
        .add_subclass(PyPatternHolder)
        .add_subclass(PyTautomerPatternHolder)
    }
}

fn wrap_fp_holder(py: Python<'_>, inner: Shared<PatternHolder>) -> PyResult<Py<PyAny>> {
    let tautomeric = read(&inner).is_tautomeric();
    let base = PyClassInitializer::from(PyFPHolderBase { inner }).add_subclass(PyPatternHolder);
    Ok(if tautomeric {
        Py::new(py, base.add_subclass(PyTautomerPatternHolder))?.into_any()
    } else {
        Py::new(py, base)?.into_any()
    })
}

// ------------------------------------------------------------------ key holders

/// Holds keys to return external references to the molecules in the molholder.
/// By default use the _Name property but can be overridden to be any property
#[pyclass(
    name = "KeyHolderBase",
    module = "rdkit.Chem.rdSubstructLibrary",
    subclass
)]
pub struct PyKeyHolderBase {
    inner: Shared<KeyFromPropHolder>,
}

#[pymethods]
impl PyKeyHolderBase {
    fn __len__(&self) -> usize {
        read(&self.inner).len() as usize
    }

    /// Add a key to the key holder, getting it from the molecule
    #[pyo3(name = "AddMol", signature = (m))]
    fn add_mol(&self, m: PyRef<'_, PyMol>) -> u32 {
        write(&self.inner).add_mol(&m.mol)
    }

    /// Add a key to the key holder, must be done in the same order as the molecules
    #[pyo3(name = "AddKey", signature = (arg1))]
    fn add_key(&self, arg1: String) -> u32 {
        write(&self.inner).add_key(arg1)
    }

    /// Return the key at the specified index
    #[pyo3(name = "GetKey", signature = (arg1))]
    fn get_key(&self, arg1: u32) -> PyResult<String> {
        Ok(read(&self.inner).get_key(arg1).py()?.to_string())
    }

    /// Returns the keys for the given indices
    #[pyo3(name = "GetKeys", signature = (indices))]
    fn get_keys(&self, indices: Vec<u32>) -> PyResult<Vec<String>> {
        read(&self.inner).get_keys(&indices).py()
    }
}

/// Holds keys to return external references to the molecules in the molholder.
/// By default use the _Name property but can be overridden to be any property
#[pyclass(name = "KeyFromPropHolder", module = "rdkit.Chem.rdSubstructLibrary", extends = PyKeyHolderBase)]
pub struct PyKeyFromPropHolder;

#[pymethods]
impl PyKeyFromPropHolder {
    #[new]
    #[pyo3(signature = (propname = "_Name"))]
    fn new(propname: &str) -> PyClassInitializer<Self> {
        PyClassInitializer::from(PyKeyHolderBase {
            inner: shared(KeyFromPropHolder::new(propname)),
        })
        .add_subclass(PyKeyFromPropHolder)
    }

    /// Return the property used to retrieve the keys
    #[pyo3(name = "GetPropName")]
    fn get_prop_name(slf: PyRef<'_, Self>) -> String {
        read(&slf.as_super().inner).prop_name().to_string()
    }
}

fn wrap_key_holder(py: Python<'_>, inner: Shared<KeyFromPropHolder>) -> PyResult<Py<PyAny>> {
    let init =
        PyClassInitializer::from(PyKeyHolderBase { inner }).add_subclass(PyKeyFromPropHolder);
    Ok(Py::new(py, init)?.into_any())
}

// ------------------------------------------------------------ the library

/// Positional and keyword arguments of the `GetMatches` family, resolved to
/// upstream's overloads: `([startIdx, endIdx,] recursionPossible=True,
/// useChirality=True, useQueryQueryMatches=False, numThreads=-1[, maxResults])`
/// or `([startIdx, endIdx,] parameters, numThreads=-1[, maxResults])`.
struct SearchArgs {
    range: Option<(u32, u32)>,
    params: SubstructMatchParameters,
    num_threads: i32,
    max_results: i32,
}

fn is_plain_int(o: &Bound<'_, PyAny>) -> bool {
    o.is_instance_of::<pyo3::types::PyInt>() && !o.is_instance_of::<PyBool>()
}

fn parse_search_args(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
    default_max_results: Option<i32>,
) -> PyResult<SearchArgs> {
    let mut positional: Vec<Bound<'_, PyAny>> = args.iter().collect();
    let mut named: HashMap<String, Bound<'_, PyAny>> = HashMap::new();
    if let Some(kw) = kwargs {
        for (k, v) in kw.iter() {
            named.insert(k.extract()?, v);
        }
    }
    // the optional index range comes first
    let mut range = None;
    if positional.len() >= 2 && is_plain_int(&positional[0]) && is_plain_int(&positional[1]) {
        range = Some((positional[0].extract()?, positional[1].extract()?));
        positional.drain(..2);
    } else if named.contains_key("startIdx") || named.contains_key("endIdx") {
        match (named.remove("startIdx"), named.remove("endIdx")) {
            (Some(s), Some(e)) => range = Some((s.extract()?, e.extract()?)),
            _ => {
                return Err(PyTypeError::new_err(
                    "startIdx and endIdx must be given together",
                ));
            }
        }
    }
    let with_parameters = named.contains_key("parameters")
        || positional
            .first()
            .is_some_and(|p| p.extract::<PyRef<'_, PySubstructMatchParameters>>().is_ok());
    let mut names: Vec<&str> = if with_parameters {
        vec!["parameters", "numThreads"]
    } else {
        vec![
            "recursionPossible",
            "useChirality",
            "useQueryQueryMatches",
            "numThreads",
        ]
    };
    if default_max_results.is_some() {
        names.push("maxResults");
    }
    if positional.len() > names.len() {
        return Err(PyTypeError::new_err("too many positional arguments"));
    }
    for (name, value) in names.iter().zip(positional) {
        if named.insert((*name).to_string(), value).is_some() {
            return Err(PyTypeError::new_err(format!(
                "multiple values for argument '{name}'"
            )));
        }
    }
    let mut take = |name: &str| named.remove(name);
    let params = if with_parameters {
        match take("parameters") {
            Some(p) => p
                .extract::<PyRef<'_, PySubstructMatchParameters>>()?
                .base_params(),
            None => return Err(PyTypeError::new_err("missing argument 'parameters'")),
        }
    } else {
        let mut p = SubstructMatchParameters {
            recursion_possible: true,
            use_chirality: true,
            use_query_query_matches: false,
            ..Default::default()
        };
        if let Some(v) = take("recursionPossible") {
            p.recursion_possible = v.extract()?;
        }
        if let Some(v) = take("useChirality") {
            p.use_chirality = v.extract()?;
        }
        if let Some(v) = take("useQueryQueryMatches") {
            p.use_query_query_matches = v.extract()?;
        }
        p
    };
    let num_threads = take("numThreads").map_or(Ok(-1), |v| v.extract())?;
    let max_results = match (take("maxResults"), default_max_results) {
        (Some(v), _) => v.extract()?,
        (None, d) => d.unwrap_or(-1),
    };
    if let Some(extra) = named.keys().next() {
        return Err(PyTypeError::new_err(format!(
            "unexpected argument '{extra}'"
        )));
    }
    Ok(SearchArgs {
        range,
        params,
        num_threads,
        max_results,
    })
}

/// A query of the `GetMatches` family.
enum Query {
    Mol(Box<rdoxide::graphmol::ROMol>),
    Tautomer(rdoxide::graphmol::tautomer_query::TautomerQuery),
    Bundle(rdoxide::graphmol::mol_bundle::MolBundle),
    Extended(rdoxide::graphmol::generalized_substruct::ExtendedQueryMol),
}

fn query_of(query: &Bound<'_, PyAny>) -> PyResult<Query> {
    if let Ok(m) = query.extract::<PyRef<'_, PyMol>>() {
        return Ok(Query::Mol(Box::new(m.mol.clone())));
    }
    if let Ok(t) = query.extract::<PyRef<'_, crate::rdtautomerquery::PyTautomerQuery>>() {
        return Ok(Query::Tautomer(t.query.clone()));
    }
    if let Ok(b) = query.extract::<PyRef<'_, crate::molbundle::PyMolBundle>>() {
        return Ok(Query::Bundle(b.bundle.clone()));
    }
    if let Ok(x) = query.extract::<PyRef<'_, crate::rdgeneralizedsubstruct::PyExtendedQueryMol>>() {
        return Ok(Query::Extended(x.query.clone()));
    }
    Err(PyTypeError::new_err(
        "query must be a Mol, a TautomerQuery, a MolBundle or an ExtendedQueryMol",
    ))
}

/// SubstructLibrary: This provides a simple API for substructure searching large datasets
/// The SubstructLibrary takes full advantage of available threads during the search operation.
#[pyclass(name = "SubstructLibrary", module = "rdkit.Chem.rdSubstructLibrary")]
pub struct PySubstructLibrary {
    lib: SubstructLibrary,
}

impl PySubstructLibrary {
    fn range(&self, range: Option<(u32, u32)>) -> (u32, u32) {
        range.unwrap_or((0, self.lib.len()))
    }
}

#[pymethods]
impl PySubstructLibrary {
    /// `SubstructLibrary()`, `SubstructLibrary(molecules[, fingerprints][, keys])` or
    /// `SubstructLibrary(pickle)`
    #[new]
    #[pyo3(signature = (molecules = None, fingerprints = None, keys = None))]
    fn new(
        molecules: Option<&Bound<'_, PyAny>>,
        fingerprints: Option<&Bound<'_, PyAny>>,
        keys: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let Some(molecules) = molecules else {
            return Ok(PySubstructLibrary {
                lib: SubstructLibrary::default(),
            });
        };
        if molecules.is_instance_of::<pyo3::types::PyString>()
            || molecules.is_instance_of::<PyBytes>()
        {
            let text = String::from_utf8_lossy(&bytes_arg(molecules)?).into_owned();
            return Ok(PySubstructLibrary {
                lib: SubstructLibrary::from_text(&text).py()?,
            });
        }
        let mols = molecules
            .extract::<PyRef<'_, PyMolHolderBase>>()?
            .inner
            .clone();
        let (mut fps, mut key_holder) = (None, None);
        // the optional holders may come in either order, or be None
        for holder in [fingerprints, keys].into_iter().flatten() {
            if holder.is_none() {
                continue;
            }
            if let Ok(f) = holder.extract::<PyRef<'_, PyFPHolderBase>>() {
                fps = Some(f.inner.clone());
            } else if let Ok(k) = holder.extract::<PyRef<'_, PyKeyHolderBase>>() {
                key_holder = Some(k.inner.clone());
            } else {
                return Err(PyTypeError::new_err(
                    "expected a fingerprint holder or a key holder",
                ));
            }
        }
        Ok(PySubstructLibrary {
            lib: SubstructLibrary::new(mols, fps, key_holder),
        })
    }

    /// Returns the molecule holder of the library
    #[pyo3(name = "GetMolHolder")]
    fn get_mol_holder(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        wrap_mol_holder(py, self.lib.mol_holder().clone())
    }

    /// Returns the fingerprint holder of the library (None if it has none)
    #[pyo3(name = "GetFpHolder")]
    fn get_fp_holder(&self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        self.lib
            .fp_holder()
            .map(|h| wrap_fp_holder(py, h.clone()))
            .transpose()
    }

    /// Returns the key holder of the library (None if it has none)
    #[pyo3(name = "GetKeyHolder")]
    fn get_key_holder(&self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        self.lib
            .key_holder()
            .map(|h| wrap_key_holder(py, h.clone()))
            .transpose()
    }

    /// Adds a molecule to the substruct library
    #[pyo3(name = "AddMol", signature = (mol))]
    fn add_mol(&self, mol: PyRef<'_, PyMol>) -> PyResult<u32> {
        self.lib.add_mol(&mol.mol).py()
    }

    /// Get the matching indices for the query
    ///
    ///   ARGUMENTS:
    ///     - query: the query molecule
    ///     - startIdx, endIdx: optional range of molecules to search
    ///     - recursionPossible, useChirality, useQueryQueryMatches: matching options
    ///       (or a SubstructMatchParameters instance as `parameters`)
    ///     - numThreads: number of threads, -1 uses all cores [default -1]
    ///     - maxResults: maximum number of results, -1 returns all [default 1000]
    #[pyo3(name = "GetMatches", signature = (query, *args, **kwargs))]
    fn get_matches(
        &self,
        py: Python<'_>,
        query: &Bound<'_, PyAny>,
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Vec<u32>> {
        let a = parse_search_args(args, kwargs, Some(1000))?;
        let (start, end) = self.range(a.range);
        let lib = self.lib.clone();
        match query_of(query)? {
            Query::Mol(query) => py.detach(|| {
                lib.get_matches(&query, start, end, &a.params, a.num_threads, a.max_results)
            }),
            Query::Tautomer(query) => py.detach(|| {
                lib.get_tautomer_matches(
                    &query,
                    start,
                    end,
                    &a.params,
                    a.num_threads,
                    a.max_results,
                )
            }),
            Query::Bundle(query) => py.detach(|| {
                lib.get_bundle_matches(&query, start, end, &a.params, a.num_threads, a.max_results)
            }),
            Query::Extended(query) => py.detach(|| {
                lib.get_xq_matches(&query, start, end, &a.params, a.num_threads, a.max_results)
            }),
        }
        .py()
    }

    /// Return the number of matches for the query
    #[pyo3(name = "CountMatches", signature = (query, *args, **kwargs))]
    fn count_matches(
        &self,
        py: Python<'_>,
        query: &Bound<'_, PyAny>,
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<u32> {
        let a = parse_search_args(args, kwargs, None)?;
        let (start, end) = self.range(a.range);
        let lib = self.lib.clone();
        match query_of(query)? {
            Query::Mol(query) => {
                py.detach(|| lib.count_matches(&query, start, end, &a.params, a.num_threads))
            }
            Query::Tautomer(query) => py.detach(|| {
                lib.count_tautomer_matches(&query, start, end, &a.params, a.num_threads)
            }),
            Query::Bundle(query) => {
                py.detach(|| lib.count_bundle_matches(&query, start, end, &a.params, a.num_threads))
            }
            Query::Extended(query) => {
                py.detach(|| lib.count_xq_matches(&query, start, end, &a.params, a.num_threads))
            }
        }
        .py()
    }

    /// Return True if any molecule in the library matches the query
    #[pyo3(name = "HasMatch", signature = (query, *args, **kwargs))]
    fn has_match(
        &self,
        py: Python<'_>,
        query: &Bound<'_, PyAny>,
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<bool> {
        let a = parse_search_args(args, kwargs, None)?;
        let (start, end) = self.range(a.range);
        let lib = self.lib.clone();
        match query_of(query)? {
            Query::Mol(query) => {
                py.detach(|| lib.has_match(&query, start, end, &a.params, a.num_threads))
            }
            Query::Tautomer(query) => {
                py.detach(|| lib.has_tautomer_match(&query, start, end, &a.params, a.num_threads))
            }
            Query::Bundle(query) => {
                py.detach(|| lib.has_bundle_match(&query, start, end, &a.params, a.num_threads))
            }
            Query::Extended(query) => {
                py.detach(|| lib.has_xq_match(&query, start, end, &a.params, a.num_threads))
            }
        }
        .py()
    }

    /// Returns a particular molecule in the molecule holder
    ///
    ///   ARGUMENTS:
    ///     - idx: which molecule to return
    ///
    ///   NOTE: molecule indices start at 0
    #[pyo3(name = "GetMol", signature = (idx))]
    fn get_mol(&self, idx: u32) -> PyResult<Option<PyMol>> {
        molecule(self.lib.mol_holder(), idx)
    }

    /// Sets the search order for the library
    ///
    ///   ARGUMENTS:
    ///     - order: sequence of molecule indices (None clears the order)
    ///
    ///   NOTE: molecule indices start at 0
    #[pyo3(name = "SetSearchOrder", signature = (seq))]
    fn set_search_order(&mut self, seq: Option<Vec<u32>>) -> PyResult<()> {
        self.lib.set_search_order(seq.unwrap_or_default()).py()
    }

    /// Returns the search order for the library
    ///
    ///   NOTE: molecule indices start at 0
    #[pyo3(name = "GetSearchOrder")]
    fn get_search_order<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.lib.search_order())
    }

    fn __len__(&self) -> usize {
        self.lib.len() as usize
    }

    /// Serialize a substructure library to a python text stream.
    /// The stream can be a file in text mode or an io.StringIO type object
    #[pyo3(name = "ToStream", signature = (stream))]
    fn to_stream(&self, stream: &Bound<'_, PyAny>) -> PyResult<()> {
        stream.call_method1("write", (self.lib.to_text().py()?,))?;
        Ok(())
    }

    /// Deserialize a substructure library from a python bytes stream.
    /// The stream is a binary file or an io.BytesIO type object
    #[pyo3(name = "InitFromStream", signature = (stream))]
    fn init_from_stream(&mut self, stream: &Bound<'_, PyAny>) -> PyResult<()> {
        let data = bytes_arg(&stream.call_method0("read")?)?;
        self.lib = SubstructLibrary::from_text(&String::from_utf8_lossy(&data)).py()?;
        Ok(())
    }

    /// Returns the serialized (pickled) library
    #[pyo3(name = "Serialize")]
    fn serialize<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, self.lib.to_text().py()?.as_bytes()))
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let data = slf.borrow().serialize(py)?;
        PyTuple::new(
            py,
            [
                slf.get_type().into_any(),
                PyTuple::new(py, [data])?.into_any(),
            ],
        )
    }
}

/// Returns True if the SubstructLibrary is serializable
#[pyfunction]
#[pyo3(name = "SubstructLibraryCanSerialize")]
fn substruct_library_can_serialize() -> bool {
    true
}

/// Add pattern fingerprints to the given library, use numThreads=-1 to use all available cores
#[pyfunction]
#[pyo3(name = "AddPatterns", signature = (sslib, patterns = None, numThreads = None))]
#[allow(non_snake_case)]
fn add_patterns(
    py: Python<'_>,
    sslib: &Bound<'_, PySubstructLibrary>,
    patterns: Option<&Bound<'_, PyAny>>,
    numThreads: Option<i32>,
) -> PyResult<()> {
    // AddPatterns(sslib, numThreads=1) or AddPatterns(sslib, patterns, numThreads=1)
    let (holder, threads) = match patterns {
        Some(p) if p.is_instance_of::<PyFPHolderBase>() => (
            Some(p.extract::<PyRef<'_, PyFPHolderBase>>()?.inner.clone()),
            numThreads.unwrap_or(1),
        ),
        Some(p) if p.is_none() => {
            return Err(PyValueError::new_err("Null PatternHolder"));
        }
        Some(p) => (None, p.extract::<i32>()?),
        None => (None, numThreads.unwrap_or(1)),
    };
    let mut lib = sslib.borrow().lib.clone();
    let lib = py
        .detach(|| lib.add_patterns(holder, threads).map(|()| lib))
        .py()?;
    sslib.borrow_mut().lib = lib;
    Ok(())
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing a class for substructure searching large sets of molecules.",
    )?;
    m.add_class::<PyMolHolderBase>()?;
    m.add_class::<PyMolHolder>()?;
    m.add_class::<PyCachedMolHolder>()?;
    m.add_class::<PyCachedSmilesMolHolder>()?;
    m.add_class::<PyCachedTrustedSmilesMolHolder>()?;
    m.add_class::<PyFPHolderBase>()?;
    m.add_class::<PyPatternHolder>()?;
    m.add_class::<PyTautomerPatternHolder>()?;
    m.add_class::<PyKeyHolderBase>()?;
    m.add_class::<PyKeyFromPropHolder>()?;
    m.add_class::<PySubstructLibrary>()?;
    m.add_function(wrap_pyfunction!(substruct_library_can_serialize, m)?)?;
    m.add_function(wrap_pyfunction!(add_patterns, m)?)?;
    Ok(())
}
