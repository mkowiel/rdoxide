//! `rdkit.Chem.rdGeneralizedSubstruct` (upstream
//! `GeneralizedSubstruct/Wrap/rdGeneralizedSubstruct.cpp`): generalized
//! substructure searching with extended query molecules.
//!
//! The binary and JSON forms are rdoxide's own (see
//! `rdoxide::graphmol::generalized_substruct`), not upstream's.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple};
use rdoxide::graphmol::generalized_substruct::{
    ExtendedQueryMol, create_extended_query_mol, has_substruct_match_xq,
    pattern_fingerprint_target_mol, substruct_match_xq,
};
use rdoxide::graphmol::substruct::{MatchVect, SubstructMatchParameters};

use crate::datastructs::PyExplicitBitVect;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmolops_adjust_query::PyAdjustQueryParameters;
use crate::substruct::PySubstructMatchParameters;

/// Extended query molecule for use in generalized substructure searching.
#[pyclass(
    name = "ExtendedQueryMol",
    module = "rdkit.Chem.rdGeneralizedSubstruct",
    skip_from_py_object
)]
pub struct PyExtendedQueryMol {
    pub(crate) query: ExtendedQueryMol,
}

#[pymethods]
impl PyExtendedQueryMol {
    /// constructor from either a binary string (from ToBinary()) or a JSON string.
    #[new]
    #[pyo3(signature = (text, isJSON = false))]
    #[allow(non_snake_case)]
    fn py_new(text: &Bound<'_, PyAny>, isJSON: bool) -> PyResult<Self> {
        let data: Vec<u8> = match text.extract::<Vec<u8>>() {
            Ok(b) => b,
            Err(_) => text.extract::<String>()?.into_bytes(),
        };
        Ok(Self {
            query: ExtendedQueryMol::from_text(&data, isJSON).py()?,
        })
    }

    #[pyo3(name = "InitFromBinary")]
    fn init_from_binary(&mut self, pkl: &[u8]) -> PyResult<()> {
        self.query = ExtendedQueryMol::from_binary(pkl).py()?;
        Ok(())
    }

    #[pyo3(name = "InitFromJSON")]
    fn init_from_json(&mut self, text: &str) -> PyResult<()> {
        self.query = ExtendedQueryMol::from_json(text).py()?;
        Ok(())
    }

    #[pyo3(name = "ToBinary")]
    fn to_binary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let data = py.detach(|| self.query.to_binary()).py()?;
        Ok(PyBytes::new(py, &data))
    }

    #[pyo3(name = "ToJSON")]
    fn to_json(&self) -> PyResult<String> {
        self.query.to_json().py()
    }

    #[pyo3(name = "PatternFingerprintQuery", signature = (fingerprintSize = 2048))]
    #[allow(non_snake_case)]
    fn pattern_fingerprint_query(&self, fingerprintSize: u32) -> PyResult<PyExplicitBitVect> {
        Ok(PyExplicitBitVect::wrap(
            self.query.pattern_fingerprint_query(fingerprintSize).py()?,
        ))
    }
}

fn search_params(
    params: Option<&Bound<'_, PySubstructMatchParameters>>,
) -> SubstructMatchParameters {
    params.map_or_else(SubstructMatchParameters::default, |p| {
        p.borrow().base_params()
    })
}

fn match_tuple<'py>(py: Python<'py>, m: &MatchVect) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, m.iter().map(|&(_, a)| a))
}

/// determines whether or not a molecule is a match to a generalized substructure query
#[pyfunction]
#[pyo3(name = "MolHasSubstructMatch", signature = (mol, query, params = None))]
fn mol_has_substruct_match(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    query: &Bound<'_, PyExtendedQueryMol>,
    params: Option<&Bound<'_, PySubstructMatchParameters>>,
) -> PyResult<bool> {
    let ps = search_params(params);
    let (mol, query) = (&mol.borrow().mol, &query.borrow().query);
    py.detach(|| has_substruct_match_xq(mol, query, &ps)).py()
}

/// returns first match (if any) of a molecule to a generalized substructure query
#[pyfunction]
#[pyo3(name = "MolGetSubstructMatch", signature = (mol, query, params = None))]
fn mol_get_substruct_match<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    query: &Bound<'py, PyExtendedQueryMol>,
    params: Option<&Bound<'py, PySubstructMatchParameters>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let mut ps = search_params(params);
    ps.max_matches = 1;
    let (mol, query) = (&mol.borrow().mol, &query.borrow().query);
    let matches = py.detach(|| substruct_match_xq(mol, query, &ps)).py()?;
    match matches.first() {
        Some(m) => match_tuple(py, m),
        None => PyTuple::new(py, Vec::<i32>::new()),
    }
}

/// returns all matches (if any) of a molecule to a generalized substructure query
#[pyfunction]
#[pyo3(name = "MolGetSubstructMatches", signature = (mol, query, params = None))]
fn mol_get_substruct_matches<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    query: &Bound<'py, PyExtendedQueryMol>,
    params: Option<&Bound<'py, PySubstructMatchParameters>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let ps = search_params(params);
    let (mol, query) = (&mol.borrow().mol, &query.borrow().query);
    let matches = py.detach(|| substruct_match_xq(mol, query, &ps)).py()?;
    let items = matches
        .iter()
        .map(|m| match_tuple(py, m))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, items)
}

/// Creates a pattern fingerprint for a target molecule that is compatible with an extended query
#[pyfunction]
#[pyo3(name = "PatternFingerprintTarget", signature = (target, fingerprintSize = 2048))]
#[allow(non_snake_case)]
fn pattern_fingerprint_target(
    target: &Bound<'_, PyMol>,
    fingerprintSize: u32,
) -> PyResult<PyExplicitBitVect> {
    Ok(PyExplicitBitVect::wrap(
        pattern_fingerprint_target_mol(&target.borrow().mol, fingerprintSize).py()?,
    ))
}

/// Creates an ExtendedQueryMol from the input molecule
///
/// This takes a query molecule and, conceptually, performs the following steps to
/// produce an ExtendedQueryMol:
///
///   1. Enumerates features like Link Nodes and SRUs
///   2. Converts everything into TautomerQueries
///   3. Runs adjustQueryProperties()
///
/// Each step is optional
#[pyfunction]
#[pyo3(
    name = "CreateExtendedQueryMol",
    signature = (mol, doEnumeration = true, doTautomers = true, adjustQueryProperties = false, adjustQueryParameters = None)
)]
#[allow(non_snake_case)]
fn create_extended_query(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    doEnumeration: bool,
    doTautomers: bool,
    adjustQueryProperties: bool,
    adjustQueryParameters: Option<&Bound<'_, PyAdjustQueryParameters>>,
) -> PyResult<PyExtendedQueryMol> {
    let params = adjustQueryParameters.map_or_else(Default::default, |p| p.borrow().params.clone());
    let mol = &mol.borrow().mol;
    let query = py
        .detach(|| {
            create_extended_query_mol(
                mol,
                doEnumeration,
                doTautomers,
                adjustQueryProperties,
                &params,
            )
        })
        .py()?;
    Ok(PyExtendedQueryMol { query })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions for generalized substructure searching",
    )?;
    m.add_class::<PyExtendedQueryMol>()?;
    m.add_function(wrap_pyfunction!(mol_has_substruct_match, m)?)?;
    m.add_function(wrap_pyfunction!(mol_get_substruct_match, m)?)?;
    m.add_function(wrap_pyfunction!(mol_get_substruct_matches, m)?)?;
    m.add_function(wrap_pyfunction!(pattern_fingerprint_target, m)?)?;
    m.add_function(wrap_pyfunction!(create_extended_query, m)?)?;
    Ok(())
}
