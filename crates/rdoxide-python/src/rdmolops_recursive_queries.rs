//! `AddRecursiveQuery`, `MolAddRecursiveQueries` and `ParseMolQueryDefFile`
//! of `rdoxide.Chem.rdmolops` (upstream `GraphMol/Wrap/MolOps.cpp`).

use std::collections::BTreeMap;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use rdoxide::graphmol::chem_transforms::recursive_queries::{
    add_recursive_queries, add_recursive_query, parse_query_def_file, parse_query_def_text,
};
use rdoxide::graphmol::romol::ROMol;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Modifies a molecule by adding a recursive query to one of its atoms
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///
///     - query: the molecule to be used as the recursive query (this will be copied)
///
///     - atomIdx: the atom to modify
///
///     - preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved
///
///   RETURNS: None
#[pyfunction]
#[pyo3(name = "AddRecursiveQuery", signature = (mol, query, atomIdx, preserveExistingQuery = true))]
#[allow(non_snake_case)]
fn py_add_recursive_query(
    mut mol: PyRefMut<'_, PyMol>,
    query: PyRef<'_, PyMol>,
    atomIdx: u32,
    #[pyo3(from_py_with = crate::pyconv::truthy)] preserveExistingQuery: bool,
) -> PyResult<()> {
    add_recursive_query(&mut mol.mol, &query.mol, atomIdx, preserveExistingQuery).py()
}

/// Adds named recursive queries to atoms
#[pyfunction]
#[pyo3(name = "MolAddRecursiveQueries", signature = (mol, queries, propName))]
#[allow(non_snake_case)]
fn py_mol_add_recursive_queries(
    mut mol: PyRefMut<'_, PyMol>,
    queries: &Bound<'_, PyDict>,
    propName: &str,
) -> PyResult<()> {
    let mut replacements: BTreeMap<String, ROMol> = BTreeMap::new();
    for (k, v) in queries.iter() {
        let m: PyRef<'_, PyMol> = v.extract()?;
        replacements.insert(k.extract()?, m.mol.clone());
    }
    add_recursive_queries(&mut mol.mol, &replacements, propName).py()?;
    Ok(())
}

/// reads query definitions from a simply formatted file
#[pyfunction]
#[pyo3(name = "ParseMolQueryDefFile", signature = (fileobj, standardize = true, delimiter = "\t", comment = "//", nameColumn = 0, smartsColumn = 1))]
#[allow(non_snake_case)]
fn py_parse_mol_query_def_file<'py>(
    py: Python<'py>,
    fileobj: &Bound<'py, PyAny>,
    #[pyo3(from_py_with = crate::pyconv::truthy)] standardize: bool,
    delimiter: &str,
    comment: &str,
    nameColumn: usize,
    smartsColumn: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let defs = if let Ok(name) = fileobj.extract::<String>() {
        parse_query_def_file(
            &name,
            standardize,
            delimiter,
            comment,
            nameColumn,
            smartsColumn,
        )
        .py()?
    } else {
        let text = fileobj.call_method0("read")?;
        let text = match text.extract::<String>() {
            Ok(t) => t,
            Err(_) => String::from_utf8_lossy(&text.extract::<Vec<u8>>()?).into_owned(),
        };
        parse_query_def_text(
            &text,
            standardize,
            delimiter,
            comment,
            nameColumn,
            smartsColumn,
        )
    };
    let res = PyDict::new(py);
    for (name, mol) in defs {
        res.set_item(name, Py::new(py, PyMol { mol })?)?;
    }
    Ok(res)
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_add_recursive_query, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_add_recursive_queries, m)?)?;
    m.add_function(wrap_pyfunction!(py_parse_mol_query_def_file, m)?)
}
