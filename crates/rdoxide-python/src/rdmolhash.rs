//! `rdoxide.Chem.rdMolHash` (upstream `GraphMol/MolHash/Wrap/rdMolHash.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::graphmol::mol_hash::{EXCLUDE_FROM_TAUTOMERISM_PROP, HashFunction, mol_hash};

use crate::enums::add_enum;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

const FUNCTIONS: [(&str, HashFunction); 19] = [
    ("AnonymousGraph", HashFunction::AnonymousGraph),
    ("ElementGraph", HashFunction::ElementGraph),
    ("CanonicalSmiles", HashFunction::CanonicalSmiles),
    ("MurckoScaffold", HashFunction::MurckoScaffold),
    ("ExtendedMurcko", HashFunction::ExtendedMurcko),
    ("MolFormula", HashFunction::MolFormula),
    ("AtomBondCounts", HashFunction::AtomBondCounts),
    ("DegreeVector", HashFunction::DegreeVector),
    ("Mesomer", HashFunction::Mesomer),
    ("HetAtomTautomer", HashFunction::HetAtomTautomer),
    ("HetAtomProtomer", HashFunction::HetAtomProtomer),
    ("RedoxPair", HashFunction::RedoxPair),
    ("Regioisomer", HashFunction::Regioisomer),
    ("NetCharge", HashFunction::NetCharge),
    ("SmallWorldIndexBR", HashFunction::SmallWorldIndexBR),
    ("SmallWorldIndexBRL", HashFunction::SmallWorldIndexBRL),
    (
        "ArthorSubstructureOrder",
        HashFunction::ArthorSubstructureOrder,
    ),
    ("HetAtomTautomerv2", HashFunction::HetAtomTautomerv2),
    ("HetAtomProtomerv2", HashFunction::HetAtomProtomerv2),
];

/// Generate a hash for a molecule. The func argument determines which hash is generated.
#[pyfunction]
#[pyo3(name = "MolHash", signature = (mol, func, useCxSmiles = false, cxFlagsToSkip = 0))]
#[allow(non_snake_case)]
fn mol_hash_py(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    func: i64,
    useCxSmiles: bool,
    cxFlagsToSkip: u32,
) -> PyResult<String> {
    let func = FUNCTIONS
        .iter()
        .map(|(_, f)| *f)
        .find(|f| *f as i64 == func)
        .ok_or_else(|| PyValueError::new_err("func must be a HashFunction value"))?;
    let mut copy = mol.mol.clone();
    py.detach(|| mol_hash(&mut copy, func, useCxSmiles, cxFlagsToSkip))
        .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to generate hashes for molecules",
    )?;
    let members: Vec<(&str, i64)> = FUNCTIONS.iter().map(|(n, f)| (*n, *f as i64)).collect();
    add_enum(m, "HashFunction", &members, None)?;
    m.add_function(wrap_pyfunction!(mol_hash_py, m)?)?;
    m.add("excludeFromTautomerismProp", EXCLUDE_FROM_TAUTOMERISM_PROP)
}
