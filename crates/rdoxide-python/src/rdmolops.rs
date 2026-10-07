//! `rdkit.Chem.rdmolops` — port of upstream `Code/GraphMol/Wrap/MolOps.cpp`
//! (ring finding so far).

use pyo3::prelude::*;
use rdoxide::graphmol::find_rings::{self, SymmetrizeSSSRAlgorithm};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

fn algorithm_from_py(obj: Option<&Bound<'_, PyAny>>) -> PyResult<SymmetrizeSSSRAlgorithm> {
    Ok(match obj.map(|o| o.extract::<i64>()).transpose()? {
        None | Some(0) => SymmetrizeSSSRAlgorithm::Default,
        Some(1) => SymmetrizeSSSRAlgorithm::Legacy,
        Some(2) => SymmetrizeSSSRAlgorithm::Rdl,
        Some(v) => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "invalid SymmetrizeSSSRAlgorithm {v}"
            )));
        }
    })
}

/// Get the smallest set of simple rings for a molecule.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use.
///     - includeDativeBonds: whether or not dative bonds should be included in the ring finding.
///     - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring finding.
///
///   RETURNS: a sequence of sequences containing the rings found as atom ids
///          The length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.
#[pyfunction]
#[pyo3(name = "GetSSSR", signature = (mol, includeDativeBonds = false, includeHydrogenBonds = false))]
#[allow(non_snake_case)]
fn get_sssr(
    mut mol: PyRefMut<'_, PyMol>,
    includeDativeBonds: bool,
    includeHydrogenBonds: bool,
) -> PyResult<Vec<Vec<i32>>> {
    find_rings::find_sssr(&mut mol.mol, includeDativeBonds, includeHydrogenBonds).py()
}

/// Get a symmetrized SSSR for a molecule.
///
///   The symmetrized SSSR is at least as large as the SSSR for a molecule.
///   In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be
///   a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).
#[pyfunction]
#[pyo3(name = "GetSymmSSSR", signature = (mol, includeDativeBonds = false, includeHydrogenBonds = false, algorithm = None, recalcSSSR = true))]
#[allow(non_snake_case)]
fn get_symm_sssr(
    mut mol: PyRefMut<'_, PyMol>,
    includeDativeBonds: bool,
    includeHydrogenBonds: bool,
    algorithm: Option<&Bound<'_, PyAny>>,
    recalcSSSR: bool,
) -> PyResult<Vec<Vec<i32>>> {
    let alg = algorithm_from_py(algorithm)?;
    find_rings::symmetrize_sssr(
        &mut mol.mol,
        alg,
        recalcSSSR,
        includeDativeBonds,
        includeHydrogenBonds,
    )
    .py()
}

/// Does a non-SSSR ring finding for a molecule.
#[pyfunction]
#[pyo3(name = "FastFindRings")]
fn fast_find_rings(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    find_rings::fast_find_rings(&mut mol.mol).py()
}

/// Generate Unique Ring Families.
#[pyfunction]
#[pyo3(name = "FindRingFamilies", signature = (mol, includeDativeBonds = false, includeHydrogenBonds = false))]
#[allow(non_snake_case)]
fn find_ring_families(
    mut mol: PyRefMut<'_, PyMol>,
    includeDativeBonds: bool,
    includeHydrogenBonds: bool,
) -> PyResult<()> {
    find_rings::find_ring_families(&mut mol.mol, includeDativeBonds, includeHydrogenBonds).py()
}

/// sets usage of the legacy symmetric SSSR code during sanitization
#[pyfunction]
#[pyo3(name = "SetUseLegacyRingFinding")]
fn set_use_legacy_ring_finding(val: bool) {
    find_rings::set_use_legacy_ring_finding(val);
}

/// returns whether or not the legacy symmetric SSSR code is being used during sanitization
#[pyfunction]
#[pyo3(name = "GetUseLegacyRingFinding")]
fn get_use_legacy_ring_finding() -> bool {
    find_rings::get_use_legacy_ring_finding()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing RDKit functionality for manipulating molecules.",
    )?;
    crate::enums::add_enum(
        m,
        "SymmetrizeSSSRAlgorithm",
        &[("DEFAULT", 0), ("LEGACY", 1), ("RDL", 2)],
        None,
    )?;
    m.add_function(wrap_pyfunction!(get_sssr, m)?)?;
    m.add_function(wrap_pyfunction!(get_symm_sssr, m)?)?;
    m.add_function(wrap_pyfunction!(fast_find_rings, m)?)?;
    m.add_function(wrap_pyfunction!(find_ring_families, m)?)?;
    m.add_function(wrap_pyfunction!(set_use_legacy_ring_finding, m)?)?;
    m.add_function(wrap_pyfunction!(get_use_legacy_ring_finding, m)?)?;
    Ok(())
}
