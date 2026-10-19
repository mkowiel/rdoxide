//! `rdoxide.Chem.rdmolfiles` — port of upstream `Code/GraphMol/Wrap/rdmolfiles.cpp`
//! (SMILES input so far).

use std::collections::BTreeMap;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use rdoxide::graphmol::smiles_parse::{self, SmilesParserParams};

use crate::mol::PyMol;

/// Parameters controlling SMILES Parsing
#[pyclass(
    name = "SmilesParserParams",
    module = "rdoxide.Chem.rdmolfiles",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PySmilesParserParams {
    /// controls the amount of debugging information produced
    #[pyo3(get, set, name = "debugParse")]
    debug_parse: i32,
    /// controls whether or not the molecule name is also parsed
    #[pyo3(get, set, name = "parseName")]
    parse_name: bool,
    /// controls whether or not the CXSMILES extensions are parsed
    #[pyo3(get, set, name = "allowCXSMILES")]
    allow_cx_smiles: bool,
    /// controls whether or not problems in CXSMILES parsing causes molecule parsing to fail
    #[pyo3(get, set, name = "strictCXSMILES")]
    strict_cx_smiles: bool,
    /// controls whether or not the molecule is sanitized before being returned
    #[pyo3(get, set)]
    sanitize: bool,
    /// controls whether or not Hs are removed before the molecule is returned
    #[pyo3(get, set, name = "removeHs")]
    remove_hs: bool,
}

#[pymethods]
impl PySmilesParserParams {
    #[new]
    fn py_new() -> Self {
        let d = SmilesParserParams::default();
        PySmilesParserParams {
            debug_parse: 0,
            parse_name: d.parse_name,
            allow_cx_smiles: d.allow_cx_smiles,
            strict_cx_smiles: d.strict_cx_smiles,
            sanitize: d.sanitize,
            remove_hs: d.remove_hs,
        }
    }
}

impl PySmilesParserParams {
    fn to_core(&self) -> SmilesParserParams {
        SmilesParserParams {
            sanitize: self.sanitize,
            allow_cx_smiles: self.allow_cx_smiles,
            strict_cx_smiles: self.strict_cx_smiles,
            parse_name: self.parse_name,
            remove_hs: self.remove_hs,
            debug_parse: self.debug_parse != 0,
            ..SmilesParserParams::default()
        }
    }
}

/// `pyObjectToString`: accepts str or bytes.
fn py_object_to_string(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = obj.extract::<String>() {
        return Ok(s);
    }
    let b: Vec<u8> = obj.extract()?;
    Ok(String::from_utf8_lossy(&b).into_owned())
}

fn parse(py: Python<'_>, smiles: &str, params: &SmilesParserParams) -> Option<PyMol> {
    // upstream catches every exception and returns None
    match py.detach(|| smiles_parse::mol_from_smiles(smiles, params)) {
        Ok(Some(mol)) => Some(PyMol { mol }),
        _ => None,
    }
}

/// Construct a molecule from a SMILES string.
///
///   ARGUMENTS:
///
///     - SMILES: the smiles string
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - replacements: (optional) a dictionary of replacement strings (see below)
///       Defaults to {}.
///
///   RETURNS:
///
///     a Mol object, None on failure.
///
/// Alternatively: MolFromSmiles(SMILES, params) with a SmilesParserParams object.
#[pyfunction]
#[pyo3(name = "MolFromSmiles", signature = (SMILES, sanitize = None, replacements = None, params = None))]
#[allow(non_snake_case)]
fn mol_from_smiles(
    py: Python<'_>,
    SMILES: &Bound<'_, PyAny>,
    sanitize: Option<&Bound<'_, PyAny>>,
    replacements: Option<&Bound<'_, PyDict>>,
    params: Option<PyRef<'_, PySmilesParserParams>>,
) -> PyResult<Option<PyMol>> {
    let smiles = py_object_to_string(SMILES)?;
    // MolFromSmiles(SMILES, params) overload
    if let Some(s) = sanitize {
        if let Ok(p) = s.cast::<PySmilesParserParams>() {
            return Ok(parse(py, &smiles, &p.borrow().to_core()));
        }
    }
    if let Some(p) = params {
        return Ok(parse(py, &smiles, &p.to_core()));
    }
    let sanitize = match sanitize {
        Some(s) => s.is_truthy()?,
        None => true,
    };
    let mut repl = BTreeMap::new();
    if let Some(d) = replacements {
        for (k, v) in d.iter() {
            repl.insert(k.extract::<String>()?, v.extract::<String>()?);
        }
    }
    let p = SmilesParserParams {
        sanitize,
        remove_hs: sanitize,
        replacements: repl,
        ..SmilesParserParams::default()
    };
    Ok(parse(py, &smiles, &p))
}

/// Returns the canonical atom ranking for each atom of a molecule fragment.
/// If breakTies is False, this returns the symmetry class for each atom.  The symmetry
/// class is used by the canonicalization routines to type each atom based on the whole
/// chemistry of the molecular graph.  Any atom with the same rank (symmetry class) is
/// indistinguishable.
#[pyfunction]
#[pyo3(name = "CanonicalRankAtoms", signature = (mol, breakTies = true, includeChirality = true, includeIsotopes = true, includeAtomMaps = true, includeChiralPresence = false))]
#[allow(non_snake_case)]
fn canonical_rank_atoms(
    mut mol: PyRefMut<'_, PyMol>,
    breakTies: bool,
    includeChirality: bool,
    includeIsotopes: bool,
    includeAtomMaps: bool,
    includeChiralPresence: bool,
) -> PyResult<Vec<u32>> {
    let opts = rdoxide::graphmol::new_canon::RankOptions {
        break_ties: breakTies,
        include_chirality: includeChirality,
        include_isotopes: includeIsotopes,
        include_atom_maps: includeAtomMaps,
        include_chiral_presence: includeChiralPresence,
        ..Default::default()
    };
    rdoxide::graphmol::new_canon::rank_mol_atoms(&mut mol.mol, opts).map_err(crate::errors::to_py)
}

/// Returns the canonical atom ranking for each atom of a molecule fragment
///   See help(CanonicalRankAtoms) for more information.
#[pyfunction]
#[pyo3(name = "CanonicalRankAtomsInFragment", signature = (mol, atomsToUse, bondsToUse = None, atomSymbols = None, breakTies = true, includeChirality = true, includeIsotopes = true, includeAtomMaps = true, includeChiralPresence = false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn canonical_rank_atoms_in_fragment(
    mut mol: PyRefMut<'_, PyMol>,
    atomsToUse: Vec<i64>,
    bondsToUse: Option<Vec<i64>>,
    atomSymbols: Option<Vec<String>>,
    breakTies: bool,
    includeChirality: bool,
    includeIsotopes: bool,
    includeAtomMaps: bool,
    includeChiralPresence: bool,
) -> PyResult<Vec<i32>> {
    let na = mol.mol.num_atoms() as usize;
    let nb = mol.mol.num_bonds() as usize;
    if atomsToUse.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "atomsToUse must not be empty",
        ));
    }
    if atomSymbols.as_ref().is_some_and(|a| a.len() != na) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "length of atom symbol list != number of atoms",
        ));
    }
    let mut atoms = vec![false; na];
    for a in atomsToUse {
        if a < 0 || a as usize >= na {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        atoms[a as usize] = true;
    }
    let mut bonds = vec![false; nb];
    for b in bondsToUse.unwrap_or_default() {
        if b < 0 || b as usize >= nb {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        bonds[b as usize] = true;
    }
    let _ = includeAtomMaps;
    // NB: upstream passes includeChiralPresence in the includeAtomMaps slot
    let opts = rdoxide::graphmol::new_canon::FragmentRankOptions {
        atom_symbols: atomSymbols.as_deref(),
        break_ties: breakTies,
        include_chirality: includeChirality,
        include_isotopes: includeIsotopes,
        include_atom_maps: includeChiralPresence,
        include_chiral_presence: false,
        ..Default::default()
    };
    let ranks =
        rdoxide::graphmol::new_canon::rank_fragment_atoms(&mut mol.mol, &atoms, &bonds, &opts)
            .map_err(crate::errors::to_py)?;
    Ok(ranks
        .iter()
        .enumerate()
        .map(|(i, &r)| if atoms[i] { r as i32 } else { -1 })
        .collect())
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing RDKit functionality for working with molecular file formats.",
    )?;
    m.add_class::<PySmilesParserParams>()?;
    crate::rdmolfiles_write::init(m)?;
    crate::rdmolfiles_smarts::init(m)?;
    crate::rdmolfiles_molblock::init(m)?;
    crate::rdmolfiles_scsr::init(m)?;
    crate::rdmolfiles_xyz::init(m)?;
    crate::rdmolfiles_cml::init(m)?;
    crate::rdmolfiles_mol2::init(m)?;
    crate::rdmolfiles_svg::init(m)?;
    crate::rdmolfiles_tpl::init(m)?;
    crate::rdmolfiles_proplists::init(m)?;
    crate::rdmolfiles_marvin::init(m)?;
    crate::rdmolfiles_pdb::init(m)?;
    crate::rdmolfiles_sequence::init(m)?;
    crate::rdmolfiles_molwrite::init(m)?;
    crate::rdmolfiles_sdsupplier::init(m)?;
    crate::rdmolfiles_fwdsdsupplier::init(m)?;
    crate::rdmolfiles_smisupplier::init(m)?;
    crate::rdmolfiles_mtsupplier::init(m)?;
    crate::rdmolfiles_sdwriter::init(m)?;
    crate::rdmolfiles_smiwriter::init(m)?;
    crate::rdmolfiles_tdtsupplier::init(m)?;
    crate::rdmolfiles_tdtwriter::init(m)?;
    m.add_function(wrap_pyfunction!(mol_from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(canonical_rank_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(canonical_rank_atoms_in_fragment, m)?)?;
    Ok(())
}
