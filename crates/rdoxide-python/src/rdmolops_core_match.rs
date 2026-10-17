//! `GetMostSubstitutedCoreMatch` and `SortMatchesByDegreeOfCoreSubstitution`
//! of `rdoxide.Chem.rdmolops` (upstream `GraphMol/Wrap/MolOps.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::substruct::{
    MatchVect, get_most_substituted_core_match, sort_matches_by_degree_of_core_substitution,
};

use crate::mol::PyMol;

/// `seqOfSeqsToMatchVectTypeVect`: the (query index, molecule index) pairs of
/// each match, given as the molecule indices in query order.
fn matches_from(matches: &Bound<'_, PyAny>) -> PyResult<Vec<MatchVect>> {
    let mut res = Vec::new();
    for item in matches.try_iter()? {
        let atoms: Vec<i32> = item?.extract()?;
        if atoms.is_empty() {
            return Err(PyValueError::new_err("tuples in matches must not be empty"));
        }
        res.push(
            atoms
                .into_iter()
                .enumerate()
                .map(|(i, a)| (i as i32, a))
                .collect(),
        );
    }
    if res.is_empty() {
        return Err(PyValueError::new_err("matches must not be empty"));
    }
    Ok(res)
}

fn match_to_tuple<'py>(py: Python<'py>, m: &MatchVect) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, m.iter().map(|&(_, a)| a))
}

/// Postprocesses the results of a mol.GetSubstructMatches(core) call
/// where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups).
/// It returns the match with the largest number of non-hydrogen matches to
/// the terminal dummy atoms.
///
///   ARGUMENTS:
///
///     - mol: the molecule GetSubstructMatches was run on
///
///     - core: the molecule used as a substructure query
///
///     - matches: the result returned by GetSubstructMatches
///
///   RETURNS: the tuple where terminal dummy atoms in the core match the largest
///            number of non-hydrogen atoms in mol
#[pyfunction]
#[pyo3(name = "GetMostSubstitutedCoreMatch", signature = (mol, core, matches))]
fn py_get_most_substituted_core_match<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    core: PyRef<'py, PyMol>,
    matches: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyTuple>> {
    let matches = matches_from(matches)?;
    match_to_tuple(
        py,
        get_most_substituted_core_match(&mol.mol, &core.mol, &matches),
    )
}

/// Postprocesses the results of a mol.GetSubstructMatches(core) call
/// where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups).
/// It returns a copy of matches sorted by decreasing number of non-hydrogen matches
/// to the terminal dummy atoms.
///
///   ARGUMENTS:
///
///     - mol: the molecule GetSubstructMatches was run on
///
///     - core: the molecule used as a substructure query
///
///     - matches: the result returned by GetSubstructMatches
///
///   RETURNS: a copy of matches sorted by decreasing number of non-hydrogen matches
///            to the terminal dummy atoms
#[pyfunction]
#[pyo3(name = "SortMatchesByDegreeOfCoreSubstitution", signature = (mol, core, matches))]
fn py_sort_matches_by_degree_of_core_substitution<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    core: PyRef<'py, PyMol>,
    matches: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyTuple>> {
    let matches = matches_from(matches)?;
    let sorted = sort_matches_by_degree_of_core_substitution(&mol.mol, &core.mol, &matches);
    let tuples = sorted
        .iter()
        .map(|m| match_to_tuple(py, m))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, tuples)
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_get_most_substituted_core_match, m)?)?;
    m.add_function(wrap_pyfunction!(
        py_sort_matches_by_degree_of_core_substitution,
        m
    )?)?;
    Ok(())
}
