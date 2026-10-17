//! `MolFromMol2File` and `MolFromMol2Block` (upstream `Wrap/rdmolfiles.cpp`).

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::mol2::{
    Mol2ParserParams, mol_from_mol2_block, mol_from_mol2_file,
};

use crate::mol::PyMol;

/// Construct a molecule from a Tripos Mol2 file.
///
///   NOTE:
///     The parser expects the atom-typing scheme used by Corina.
///     Atom types from Tripos' dbtranslate are less supported.
///     Other atom typing schemes are unlikely to work.
///
///   ARGUMENTS:
///
///     - mol2FileName: name of the file to read
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to true.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///     - cleanupSubstructures: (optional) toggles standardizing some
///       substructures found in mol2 files.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromMol2File", signature = (mol2FileName, sanitize = true, removeHs = true, cleanupSubstructures = true))]
#[allow(non_snake_case)]
fn py_mol_from_mol2_file(
    py: Python<'_>,
    #[pyo3(from_py_with = crate::pyconv::path_like)] mol2FileName: String,
    #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] removeHs: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] cleanupSubstructures: bool,
) -> PyResult<Option<PyMol>> {
    let params = Mol2ParserParams {
        sanitize,
        remove_hs: removeHs,
        cleanup_substructures: cleanupSubstructures,
    };
    match py.detach(|| mol_from_mol2_file(&mol2FileName, &params)) {
        Ok(mol) => Ok(mol.map(|mol| PyMol { mol })),
        // a missing file is an IOError, any other failure gives None
        Err(e @ RDKitError::BadFile(_)) => Err(crate::errors::to_py(e)),
        Err(_) => Ok(None),
    }
}

/// Construct a molecule from a Tripos Mol2 block.
///
///   NOTE:
///     The parser expects the atom-typing scheme used by Corina.
///     Atom types from Tripos' dbtranslate are less supported.
///     Other atom typing schemes are unlikely to work.
///
///   ARGUMENTS:
///
///     - mol2Block: string containing the Mol2 block
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///     - cleanupSubstructures: (optional) toggles standardizing some
///       substructures found in mol2 files.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromMol2Block", signature = (mol2Block, sanitize = true, removeHs = true, cleanupSubstructures = true))]
#[allow(non_snake_case)]
fn py_mol_from_mol2_block(
    py: Python<'_>,
    mol2Block: &str,
    #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] removeHs: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] cleanupSubstructures: bool,
) -> Option<PyMol> {
    let params = Mol2ParserParams {
        sanitize,
        remove_hs: removeHs,
        cleanup_substructures: cleanupSubstructures,
    };
    py.detach(|| mol_from_mol2_block(mol2Block, &params))
        .ok()
        .flatten()
        .map(|mol| PyMol { mol })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_from_mol2_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_from_mol2_block, m)?)
}
