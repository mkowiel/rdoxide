//! `rdoxide.Chem.rdmolops` — port of upstream `Code/GraphMol/Wrap/MolOps.cpp`
//! (ring finding so far).

use pyo3::prelude::*;
use rdoxide::graphmol::aromaticity::{self, AromaticityModel};
use rdoxide::graphmol::conjug_hybrid;
use rdoxide::graphmol::find_rings::{self, SymmetrizeSSSRAlgorithm};
use rdoxide::graphmol::kekulize;
use rdoxide::graphmol::mol_ops::{self, sanitize_flags as sf};

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

/// Kekulizes the molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the
///       molecule will be marked non-aromatic following the kekulization.
///       Default value is False.
///     - canonical: (optional) whether to canonicalize the kekulization.
///
///   NOTES:
///
///     - The molecule is modified in place.
///     - this does not modify query bonds which have bond type queries (like those
///       which come from SMARTS) or rings containing them.
///     - even if clearAromaticFlags is False the BondType for all modified
///       aromatic bonds will be changed from AROMATIC to SINGLE or DOUBLE
///       Kekulization.
#[pyfunction]
#[pyo3(name = "Kekulize", signature = (mol, clearAromaticFlags = false, canonical = true))]
#[allow(non_snake_case)]
fn py_kekulize(
    mut mol: PyRefMut<'_, PyMol>,
    clearAromaticFlags: bool,
    canonical: bool,
) -> PyResult<()> {
    kekulize::kekulize(
        &mut mol.mol,
        clearAromaticFlags,
        canonical,
        kekulize::DEFAULT_MAX_BACKTRACKS,
    )
    .py()
}

/// Kekulizes the molecule if possible. Otherwise the molecule is not modified
#[pyfunction]
#[pyo3(name = "KekulizeIfPossible", signature = (mol, clearAromaticFlags = false, canonical = true))]
#[allow(non_snake_case)]
fn py_kekulize_if_possible(
    mut mol: PyRefMut<'_, PyMol>,
    clearAromaticFlags: bool,
    canonical: bool,
) -> PyResult<()> {
    kekulize::kekulize_if_possible(
        &mut mol.mol,
        clearAromaticFlags,
        canonical,
        kekulize::DEFAULT_MAX_BACKTRACKS,
    )
    .py()
    .map(|_| ())
}

/// Kekulize, check valencies, set aromaticity, conjugation and hybridization
///
///     - The molecule is modified in place.
///
///     - If sanitization fails, an exception will be thrown unless catchErrors is set
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///     - sanitizeOps: (optional) sanitization operations to be carried out
///       these should be constructed by or'ing together the
///       operations in rdoxide.Chem.SanitizeFlags
///     - catchErrors: (optional) if provided, instead of raising an exception
///       when sanitization fails (the default behavior), the
///       operations that failed (as defined in rdoxide.Chem.SanitizeFlags)
///       is returned. Zero is returned on success.
#[pyfunction]
#[pyo3(name = "SanitizeMol", signature = (mol, sanitizeOps = None, catchErrors = false))]
#[allow(non_snake_case)]
fn py_sanitize_mol(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    sanitizeOps: Option<u64>,
    catchErrors: bool,
) -> PyResult<Py<PyAny>> {
    let ops = sanitizeOps.map_or(sf::SANITIZE_ALL, |v| v as u32);
    let failed = match mol_ops::sanitize_mol(&mut mol.mol, ops) {
        Ok(()) => 0,
        Err((flag, e)) => {
            if !catchErrors {
                return Err(crate::errors::to_py(e));
            }
            flag
        }
    };
    crate::enums::to_py(py, "SanitizeFlags", i64::from(failed))
}

/// Does aromaticity perception
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - model: the model to use
#[pyfunction]
#[pyo3(name = "SetAromaticity", signature = (mol, model = 0))]
fn py_set_aromaticity(mut mol: PyRefMut<'_, PyMol>, model: i64) -> PyResult<()> {
    let model = match model {
        0x0 => AromaticityModel::Default,
        0x1 => AromaticityModel::Rdkit,
        0x2 => AromaticityModel::Simple,
        0x4 => AromaticityModel::Mdl,
        0x8 => AromaticityModel::Mmff94,
        _ => AromaticityModel::Custom,
    };
    aromaticity::set_aromaticity(&mut mol.mol, model)
        .py()
        .map(|_| ())
}

/// finds conjugated bonds
#[pyfunction]
#[pyo3(name = "SetConjugation")]
fn py_set_conjugation(mut mol: PyRefMut<'_, PyMol>) {
    conjug_hybrid::set_conjugation(&mut mol.mol);
}

/// Assigns hybridization states to atoms
#[pyfunction]
#[pyo3(name = "SetHybridization")]
fn py_set_hybridization(mut mol: PyRefMut<'_, PyMol>) {
    conjug_hybrid::set_hybridization(&mut mol.mol);
}

/// Assigns radical counts to atoms
#[pyfunction]
#[pyo3(name = "AssignRadicals")]
fn py_assign_radicals(mut mol: PyRefMut<'_, PyMol>) {
    mol_ops::assign_radicals(&mut mol.mol);
}

/// cleans up certain common bad functionalities in the molecule
#[pyfunction]
#[pyo3(name = "Cleanup")]
fn py_cleanup(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    mol_ops::clean_up(&mut mol.mol).py()
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
    crate::enums::add_enum_exported(
        m,
        "SanitizeFlags",
        &[
            ("SANITIZE_NONE", 0x0),
            ("SANITIZE_CLEANUP", 0x1),
            ("SANITIZE_PROPERTIES", 0x2),
            ("SANITIZE_SYMMRINGS", 0x4),
            ("SANITIZE_KEKULIZE", 0x8),
            ("SANITIZE_FINDRADICALS", 0x10),
            ("SANITIZE_SETAROMATICITY", 0x20),
            ("SANITIZE_SETCONJUGATION", 0x40),
            ("SANITIZE_SETHYBRIDIZATION", 0x80),
            ("SANITIZE_CLEANUPCHIRALITY", 0x100),
            ("SANITIZE_ADJUSTHS", 0x200),
            ("SANITIZE_CLEANUP_ORGANOMETALLICS", 0x400),
            ("SANITIZE_CLEANUPATROPISOMERS", 0x800),
            ("SANITIZE_ALL", 0xFFF_FFFF),
        ],
        None,
    )?;
    crate::enums::add_enum_exported(
        m,
        "AromaticityModel",
        &[
            ("AROMATICITY_DEFAULT", 0x0),
            ("AROMATICITY_RDKIT", 0x1),
            ("AROMATICITY_SIMPLE", 0x2),
            ("AROMATICITY_MDL", 0x4),
            ("AROMATICITY_MMFF94", 0x8),
            ("AROMATICITY_CUSTOM", 0xFFF_FFFF),
        ],
        None,
    )?;
    crate::rdmolops_hs::init(m)?;
    m.add_function(wrap_pyfunction!(py_kekulize, m)?)?;
    m.add_function(wrap_pyfunction!(py_kekulize_if_possible, m)?)?;
    m.add_function(wrap_pyfunction!(py_sanitize_mol, m)?)?;
    crate::rdmolops_subgraphs::add(m)?;
    crate::rdmolops_fingerprints::add(m)?;
    crate::rdmolops_matrices::add(m)?;
    crate::rdmolops_stereo_groups::add(m)?;
    crate::rdmolops_adjust_query::add(m)?;
    crate::rdmolops_stereo_misc::add(m)?;
    crate::rdchem_stereo_info::add_ops(m)?;
    crate::rdmolops_structure::add(m)?;
    crate::rdmolops_subset::add(m)?;
    crate::rdmolops_fragments::add(m)?;
    crate::rdmolops_transforms::add(m)?;
    crate::rdmolops_wedging::add(m)?;
    crate::rdmolops_zip::add(m)?;
    m.add_function(wrap_pyfunction!(py_set_aromaticity, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_conjugation, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_hybridization, m)?)?;
    m.add_function(wrap_pyfunction!(py_assign_radicals, m)?)?;
    m.add_function(wrap_pyfunction!(py_cleanup, m)?)?;
    m.add_function(wrap_pyfunction!(get_sssr, m)?)?;
    m.add_function(wrap_pyfunction!(get_symm_sssr, m)?)?;
    m.add_function(wrap_pyfunction!(fast_find_rings, m)?)?;
    m.add_function(wrap_pyfunction!(find_ring_families, m)?)?;
    m.add_function(wrap_pyfunction!(set_use_legacy_ring_finding, m)?)?;
    m.add_function(wrap_pyfunction!(get_use_legacy_ring_finding, m)?)?;
    Ok(())
}
