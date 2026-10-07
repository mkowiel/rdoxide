//! `rdkit.Chem.rdmolfiles` — port of upstream `Code/GraphMol/Wrap/rdmolfiles.cpp`
//! (SMILES input so far).

use std::collections::BTreeMap;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use rdoxide::graphmol::smiles_parse::{self, SmilesParserParams};

use crate::mol::PyMol;

/// Parameters controlling SMILES Parsing
#[pyclass(
    name = "SmilesParserParams",
    module = "rdkit.Chem.rdmolfiles",
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

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing RDKit functionality for working with molecular file formats.",
    )?;
    m.add_class::<PySmilesParserParams>()?;
    m.add_function(wrap_pyfunction!(mol_from_smiles, m)?)?;
    Ok(())
}
