//! `rdkit.Chem.rdReducedGraphs` (upstream
//! `GraphMol/ReducedGraphs/Wrap/rdReducedGraphs.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::graphmol::reduced_graphs::{
    generate_erg_fingerprint_for_reduced_graph, generate_mol_extended_reduced_graph,
    get_erg_fingerprint,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::pyconv::float64_array;

/// Custom atom types are not supported, as upstream.
fn reject_atom_types(atom_types: Option<&Bound<'_, PyAny>>) -> PyResult<()> {
    match atom_types {
        Some(t) if t.is_truthy()? => Err(PyValueError::new_err(
            "specification of atom types not yet supported",
        )),
        _ => Ok(()),
    }
}

/// Returns the reduced graph for a molecule
#[pyfunction]
#[pyo3(name = "GenerateMolExtendedReducedGraph", signature = (mol, atomTypes = None))]
#[allow(non_snake_case)]
fn generate_mol_extended_reduced_graph_py(
    mol: PyRef<'_, PyMol>,
    atomTypes: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    reject_atom_types(atomTypes)?;
    let reduced = generate_mol_extended_reduced_graph(&mol.mol, None).py()?;
    Ok(PyMol { mol: reduced })
}

/// Returns the ErG fingerprint vector for a reduced graph
#[pyfunction]
#[pyo3(name = "GenerateErGFingerprintForReducedGraph", signature = (mol, atomTypes = None, fuzzIncrement = 0.3, minPath = 1, maxPath = 15))]
#[allow(non_snake_case)]
fn generate_erg_fingerprint_for_reduced_graph_py<'py>(
    py: Python<'py>,
    mol: PyRef<'_, PyMol>,
    atomTypes: Option<&Bound<'_, PyAny>>,
    fuzzIncrement: f64,
    minPath: u32,
    maxPath: u32,
) -> PyResult<Bound<'py, PyAny>> {
    reject_atom_types(atomTypes)?;
    let fp =
        generate_erg_fingerprint_for_reduced_graph(&mol.mol, None, fuzzIncrement, minPath, maxPath)
            .py()?;
    float64_array(py, &fp)
}

/// Returns the ErG fingerprint vector for a molecule
#[pyfunction]
#[pyo3(name = "GetErGFingerprint", signature = (mol, atomTypes = None, fuzzIncrement = 0.3, minPath = 1, maxPath = 15))]
#[allow(non_snake_case)]
fn get_erg_fingerprint_py<'py>(
    py: Python<'py>,
    mol: PyRef<'_, PyMol>,
    atomTypes: Option<&Bound<'_, PyAny>>,
    fuzzIncrement: f64,
    minPath: u32,
    maxPath: u32,
) -> PyResult<Bound<'py, PyAny>> {
    reject_atom_types(atomTypes)?;
    let fp = get_erg_fingerprint(&mol.mol, None, fuzzIncrement, minPath, maxPath).py()?;
    float64_array(py, &fp)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to generate and work with reduced graphs",
    )?;
    m.add_function(wrap_pyfunction!(generate_mol_extended_reduced_graph_py, m)?)?;
    m.add_function(wrap_pyfunction!(
        generate_erg_fingerprint_for_reduced_graph_py,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(get_erg_fingerprint_py, m)?)
}
