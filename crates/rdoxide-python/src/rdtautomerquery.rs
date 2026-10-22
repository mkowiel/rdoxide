//! `rdoxide.Chem.rdTautomerQuery` (upstream `TautomerQuery/Wrap/rdTautomerQuery.cpp`):
//! substructure queries that match all tautomers of a molecule.
//!
//! Python callbacks in `SubstructMatchParameters` are not used by these
//! searches, and the query cannot be serialized (upstream needs
//! boost::serialization for that).

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyTuple};
use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::generalized_substruct::{
    tautomer_query_from_binary, tautomer_query_to_binary,
};
use rdoxide::graphmol::substruct::{MatchVect, SubstructMatchParameters};
use rdoxide::graphmol::tautomer_query::TautomerQuery;

use crate::datastructs::PyExplicitBitVect;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::substruct::PySubstructMatchParameters;

/// The tautomer query class: a query that matches all the tautomeric forms
/// of a molecule.
#[pyclass(
    name = "TautomerQuery",
    module = "rdoxide.Chem.rdTautomerQuery",
    skip_from_py_object
)]
pub struct PyTautomerQuery {
    pub(crate) query: TautomerQuery,
}

fn py_mol(py: Python<'_>, mol: &ROMol) -> PyResult<Py<PyMol>> {
    Py::new(
        py,
        PyMol {
            mol: ROMol::copy_from(mol, false, -1),
        },
    )
}

fn match_tuple<'py>(py: Python<'py>, m: &MatchVect) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, m.iter().map(|&(_, a)| a))
}

/// Positional (after the target) or keyword argument.
fn arg<'py>(
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    pos: usize,
    key: &str,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    if let Some(v) = args.iter().nth(pos) {
        return Ok(Some(v));
    }
    kwargs.map_or(Ok(None), |k| k.get_item(key))
}

fn flag(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
    pos: usize,
    key: &str,
    default: bool,
) -> PyResult<bool> {
    arg(args, kwargs, pos, key)?.map_or(Ok(default), |v| v.is_truthy())
}

/// The search parameters: a `SubstructMatchParameters` as first extra
/// argument (or `params=`), else the individual flags.
fn parameters(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
    with_uniquify: bool,
    max_matches: u32,
) -> PyResult<SubstructMatchParameters> {
    if let Some(p) = arg(args, kwargs, 0, "params")? {
        if let Ok(p) = p.cast::<PySubstructMatchParameters>() {
            return Ok(p.borrow().base_params());
        }
    }
    let offset = usize::from(with_uniquify);
    let mut ps = SubstructMatchParameters {
        use_chirality: flag(args, kwargs, offset, "useChirality", false)?,
        use_query_query_matches: flag(args, kwargs, offset + 1, "useQueryQueryMatches", false)?,
        max_matches,
        ..SubstructMatchParameters::default()
    };
    if with_uniquify {
        ps.uniquify = flag(args, kwargs, 0, "uniquify", true)?;
        if let Some(m) = arg(args, kwargs, 3, "maxMatches")? {
            ps.max_matches = m.extract()?;
        }
    }
    Ok(ps)
}

#[pymethods]
impl PyTautomerQuery {
    #[new]
    #[pyo3(signature = (mol, tautomerTransformFile = None))]
    #[allow(non_snake_case)]
    fn py_new(mol: &Bound<'_, PyAny>, tautomerTransformFile: Option<&str>) -> PyResult<Self> {
        if let Ok(pkl) = mol.cast::<pyo3::types::PyBytes>() {
            return Ok(PyTautomerQuery {
                query: tautomer_query_from_binary(pkl.as_bytes()).py()?,
            });
        }
        let mol = mol.cast::<PyMol>()?.borrow();
        Ok(PyTautomerQuery {
            query: TautomerQuery::from_mol(&mol.mol, tautomerTransformFile.unwrap_or("")).py()?,
        })
    }

    /// Returns True if the target contains a match of any tautomer of the query.
    /// Optional arguments: recursionPossible, useChirality, useQueryQueryMatches,
    /// or a SubstructMatchParameters instance.
    #[pyo3(name = "IsSubstructOf", signature = (target, *args, **kwargs))]
    fn is_substruct_of(
        &self,
        py: Python<'_>,
        target: &Bound<'_, PyMol>,
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<bool> {
        let mut ps = parameters_for_single(args, kwargs, true)?;
        ps.max_matches = 1;
        let mol = target.borrow().mol.clone();
        let query = &self.query;
        py.detach(|| query.is_substruct_of(&mol, &ps)).py()
    }

    /// Returns the indices of the target atoms matching the query (first match only).
    #[pyo3(name = "GetSubstructMatch", signature = (target, *args, **kwargs))]
    fn get_substruct_match<'py>(
        &self,
        py: Python<'py>,
        target: &Bound<'py, PyMol>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let mut ps = parameters_for_single(args, kwargs, false)?;
        ps.max_matches = 1;
        let mol = target.borrow().mol.clone();
        let query = &self.query;
        let (matches, _) = py.detach(|| query.substruct_of(&mol, &ps)).py()?;
        match matches.first() {
            Some(m) => match_tuple(py, m),
            None => Ok(PyTuple::empty(py)),
        }
    }

    /// Returns tuples of the indices of the target atoms matching the query.
    /// Optional arguments: uniquify, useChirality, useQueryQueryMatches,
    /// maxMatches, or a SubstructMatchParameters instance.
    #[pyo3(name = "GetSubstructMatches", signature = (target, *args, **kwargs))]
    fn get_substruct_matches<'py>(
        &self,
        py: Python<'py>,
        target: &Bound<'py, PyMol>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let ps = parameters(args, kwargs, true, 1000)?;
        let mol = target.borrow().mol.clone();
        let query = &self.query;
        let (matches, _) = py.detach(|| query.substruct_of(&mol, &ps)).py()?;
        let items = matches
            .iter()
            .map(|m| match_tuple(py, m))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    /// Like GetSubstructMatches, but returns (match, matching tautomer) pairs.
    #[pyo3(name = "GetSubstructMatchesWithTautomers", signature = (target, *args, **kwargs))]
    fn get_substruct_matches_with_tautomers<'py>(
        &self,
        py: Python<'py>,
        target: &Bound<'py, PyMol>,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let ps = parameters(args, kwargs, true, 1000)?;
        let mol = target.borrow().mol.clone();
        let query = &self.query;
        let (matches, tautomers) = py.detach(|| query.substruct_of(&mol, &ps)).py()?;
        let items = matches
            .iter()
            .zip(&tautomers)
            .map(|(m, t)| {
                PyTuple::new(
                    py,
                    [
                        match_tuple(py, m)?.into_any(),
                        py_mol(py, t)?.into_bound(py).into_any(),
                    ],
                )
            })
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    /// Returns the pattern fingerprint of the template molecule.
    #[pyo3(name = "PatternFingerprintTemplate", signature = (fingerprintSize = 2048))]
    #[allow(non_snake_case)]
    fn pattern_fingerprint_template(&self, fingerprintSize: u32) -> PyResult<PyExplicitBitVect> {
        Ok(PyExplicitBitVect::wrap(
            self.query
                .pattern_fingerprint_template(fingerprintSize)
                .py()?,
        ))
    }

    /// Returns the template molecule (a copy).
    #[pyo3(name = "GetTemplateMolecule")]
    fn template_molecule(&self, py: Python<'_>) -> PyResult<Py<PyMol>> {
        py_mol(py, self.query.template_molecule())
    }

    /// Returns the indices of the atoms modified by the tautomer transforms.
    #[pyo3(name = "GetModifiedAtoms")]
    fn modified_atoms(&self) -> Vec<u32> {
        self.query.modified_atoms().to_vec()
    }

    /// Returns the indices of the bonds modified by the tautomer transforms.
    #[pyo3(name = "GetModifiedBonds")]
    fn modified_bonds(&self) -> Vec<u32> {
        self.query.modified_bonds().to_vec()
    }

    /// Returns the tautomers of the query molecule.
    #[pyo3(name = "GetTautomers")]
    fn tautomers<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let items = self
            .query
            .tautomers()
            .iter()
            .map(|t| py_mol(py, t))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    /// Returns the binary form (a pickle) of the TautomerQuery.
    #[pyo3(name = "ToBinary")]
    fn to_binary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyBytes>> {
        let data = py.detach(|| tautomer_query_to_binary(&self.query)).py()?;
        Ok(pyo3::types::PyBytes::new(py, &data))
    }

    /// `(type(self), (self.ToBinary(),))`
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let args = PyTuple::new(py, [slf.call_method0("ToBinary")?])?;
        PyTuple::new(py, [slf.get_type().into_any(), args.into_any()])
    }
}

/// Parameters of the single-result searches: a `SubstructMatchParameters`, or
/// `recursionPossible` (when `with_recursion`), `useChirality` and
/// `useQueryQueryMatches`.
fn parameters_for_single(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
    with_recursion: bool,
) -> PyResult<SubstructMatchParameters> {
    if let Some(p) = arg(args, kwargs, 0, "params")? {
        if let Ok(p) = p.cast::<PySubstructMatchParameters>() {
            return Ok(p.borrow().base_params());
        }
    }
    let offset = usize::from(with_recursion);
    let mut ps = SubstructMatchParameters {
        use_chirality: flag(args, kwargs, offset, "useChirality", false)?,
        use_query_query_matches: flag(args, kwargs, offset + 1, "useQueryQueryMatches", false)?,
        ..SubstructMatchParameters::default()
    };
    if with_recursion {
        ps.recursion_possible = flag(args, kwargs, 0, "recursionPossible", true)?;
    }
    Ok(ps)
}

/// Returns the pattern fingerprint of a target, for screening before a tautomer search.
#[pyfunction]
#[pyo3(name = "PatternFingerprintTautomerTarget", signature = (target, fingerprintSize = 2048))]
#[allow(non_snake_case)]
fn pattern_fingerprint_tautomer_target(
    target: &Bound<'_, PyAny>,
    fingerprintSize: u32,
) -> PyResult<PyExplicitBitVect> {
    let mol = target
        .cast::<PyMol>()
        .map_err(|_| PyTypeError::new_err("target must be a Mol"))?
        .borrow();
    Ok(PyExplicitBitVect::wrap(
        TautomerQuery::pattern_fingerprint_target(&mol.mol, fingerprintSize).py()?,
    ))
}

/// Returns True if the TautomerQuery is serializable (requires that the RDKit was built with boost::serialization)
#[pyfunction]
#[pyo3(name = "TautomerQueryCanSerialize")]
fn tautomer_query_can_serialize() -> bool {
    true
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTautomerQuery>()?;
    m.add_function(wrap_pyfunction!(pattern_fingerprint_tautomer_target, m)?)?;
    m.add_function(wrap_pyfunction!(tautomer_query_can_serialize, m)?)?;
    Ok(())
}
