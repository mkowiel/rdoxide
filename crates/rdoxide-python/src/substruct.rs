//! Substructure matching on `Mol` (`HasSubstructMatch`, `GetSubstructMatch`,
//! `GetSubstructMatches`) and `SubstructMatchParameters` — port of the
//! corresponding parts of upstream `Code/GraphMol/Wrap/Mol.cpp` and
//! `substructmethods.h`.

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Parameters controlling substructure matching
#[pyclass(
    name = "SubstructMatchParameters",
    module = "rdkit.Chem.rdchem",
    get_all,
    set_all,
    skip_from_py_object
)]
#[derive(Clone)]
#[allow(non_snake_case)]
pub struct PySubstructMatchParameters {
    /// Use chirality in determining whether or not atoms/bonds match
    useChirality: bool,
    /// take enhanced stereochemistry into account while doing the match. This only has an effect if useChirality is also True.
    useEnhancedStereo: bool,
    /// aromatic and conjugated bonds match each other
    aromaticMatchesConjugated: bool,
    /// aromatic and single or double bonds match each other
    aromaticMatchesSingleOrDouble: bool,
    /// use generic groups (=homology groups) as a post-filtering step (if any are present in the molecule)
    useGenericMatchers: bool,
    /// Consider query-query matches, not just simple matches
    useQueryQueryMatches: bool,
    /// Allow recursive queries
    recursionPossible: bool,
    /// uniquify (by atom index) match results
    uniquify: bool,
    /// maximum number of matches to return
    maxMatches: u32,
    /// maximum number of recursive matches to find
    maxRecursiveMatches: u32,
    /// number of threads to use when multi-threading is possible.0 selects the number of concurrent threads supported by thehardware. negative values are added to the number of concurrentthreads supported by the hardware.
    numThreads: i32,
    /// bond properties that must be equivalent in order to match.
    bondProperties: Vec<String>,
    /// atom properties that must be equivalent in order to match.
    atomProperties: Vec<String>,
    /// If set, query atoms and bonds with specified stereochemistry will match atoms and bonds with unspecified stereochemistry.
    specifiedStereoQueryMatchesUnspecified: bool,
}

#[pymethods]
impl PySubstructMatchParameters {
    #[new]
    fn py_new() -> Self {
        let d = SubstructMatchParameters::default();
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
            bondProperties: d.bond_properties,
            atomProperties: d.atom_properties,
            specifiedStereoQueryMatchesUnspecified: d.specified_stereo_query_matches_unspecified,
        }
    }
}

impl PySubstructMatchParameters {
    pub(crate) fn to_core(&self) -> SubstructMatchParameters {
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
        }
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
    fn query(&self) -> PyResult<Bound<'py, PyMol>> {
        let q = self.get(0, "query")?.ok_or_else(|| {
            PyTypeError::new_err(format!("{}() missing required argument 'query'", self.name))
        })?;
        Ok(q.cast::<PyMol>()?.clone())
    }
    fn params(&self) -> Option<SubstructMatchParameters> {
        self.get(1, "params").ok().flatten().and_then(|p| {
            p.cast::<PySubstructMatchParameters>()
                .ok()
                .map(|p| p.borrow().to_core())
        })
    }
}

fn run(
    py: Python<'_>,
    mol: &PyMol,
    query: &Bound<'_, PyMol>,
    ps: &SubstructMatchParameters,
) -> PyResult<Vec<MatchVect>> {
    let q = query.borrow();
    let (m, q) = (&mol.mol, &q.mol);
    py.detach(|| substruct_match(m, q, ps)).py()
}

fn match_tuple<'py>(py: Python<'py>, m: &MatchVect) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, m.iter().map(|&(_, a)| a))
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
    fn has_substruct_match(
        &self,
        py: Python<'_>,
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<bool> {
        let a = Args {
            name: "HasSubstructMatch",
            args,
            kwargs,
        };
        let query = a.query()?;
        let ps = match a.params() {
            Some(p) => p,
            None => SubstructMatchParameters {
                recursion_possible: a.flag(1, "recursionPossible", true)?,
                use_chirality: a.flag(2, "useChirality", false)?,
                use_query_query_matches: a.flag(3, "useQueryQueryMatches", false)?,
                max_matches: 1,
                ..SubstructMatchParameters::default()
            },
        };
        Ok(!run(py, self, &query, &ps)?.is_empty())
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
        &self,
        py: Python<'py>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let a = Args {
            name: "GetSubstructMatch",
            args,
            kwargs,
        };
        let query = a.query()?;
        let ps = match a.params() {
            Some(p) => p,
            None => SubstructMatchParameters {
                use_chirality: a.flag(1, "useChirality", false)?,
                use_query_query_matches: a.flag(2, "useQueryQueryMatches", false)?,
                max_matches: 1,
                ..SubstructMatchParameters::default()
            },
        };
        let res = run(py, self, &query, &ps)?;
        match res.first() {
            Some(m) => match_tuple(py, m),
            None => Ok(PyTuple::empty(py)),
        }
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
        &self,
        py: Python<'py>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let a = Args {
            name: "GetSubstructMatches",
            args,
            kwargs,
        };
        let query = a.query()?;
        let ps = match a.params() {
            Some(p) => p,
            None => SubstructMatchParameters {
                uniquify: a.flag(1, "uniquify", true)?,
                use_chirality: a.flag(2, "useChirality", false)?,
                use_query_query_matches: a.flag(3, "useQueryQueryMatches", false)?,
                max_matches: match a.get(4, "maxMatches")? {
                    Some(v) => v.extract()?,
                    None => 1000,
                },
                ..SubstructMatchParameters::default()
            },
        };
        let res = run(py, self, &query, &ps)?;
        let items = res
            .iter()
            .map(|m| match_tuple(py, m))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySubstructMatchParameters>()?;
    Ok(())
}
