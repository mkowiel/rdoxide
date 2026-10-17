//! `rdoxide.Chem.rdmolfiles`: SMARTS parsing — port of the corresponding parts
//! of upstream `Code/GraphMol/Wrap/rdmolfiles.cpp`.

use std::collections::BTreeMap;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::smiles_parse::smarts::{
    SmartsParserParams, atom_from_smarts, bond_from_smarts, mol_from_smarts,
};
use rdoxide::graphmol::smiles_parse::{atom_from_smiles, bond_from_smiles};

use crate::atom::PyAtom;
use crate::bond::PyBond;
use crate::mol::PyMol;

/// Parameters controlling SMARTS Parsing
#[pyclass(
    name = "SmartsParserParams",
    module = "rdoxide.Chem.rdmolfiles",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PySmartsParserParams {
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
    /// toggles merging H atoms in the SMARTS into neighboring atoms
    #[pyo3(get, set, name = "mergeHs")]
    merge_hs: bool,
}

#[pymethods]
impl PySmartsParserParams {
    #[new]
    fn py_new() -> Self {
        let d = SmartsParserParams::default();
        PySmartsParserParams {
            debug_parse: 0,
            parse_name: d.parse_name,
            allow_cx_smiles: d.allow_cx_smiles,
            strict_cx_smiles: d.strict_cx_smiles,
            merge_hs: d.merge_hs,
        }
    }
}

impl PySmartsParserParams {
    fn to_core(&self) -> SmartsParserParams {
        SmartsParserParams {
            debug_parse: self.debug_parse != 0,
            parse_name: self.parse_name,
            allow_cx_smiles: self.allow_cx_smiles,
            strict_cx_smiles: self.strict_cx_smiles,
            merge_hs: self.merge_hs,
            ..SmartsParserParams::default()
        }
    }
}

fn py_object_to_string(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = obj.extract::<String>() {
        return Ok(s);
    }
    let b: Vec<u8> = obj.extract()?;
    Ok(String::from_utf8_lossy(&b).into_owned())
}

fn parse(py: Python<'_>, smarts: &str, params: &SmartsParserParams) -> Option<PyMol> {
    // upstream catches every exception and returns None
    match py.detach(|| mol_from_smarts(smarts, params)) {
        Ok(Some(mol)) => Some(PyMol { mol }),
        _ => None,
    }
}

/// Construct a molecule from a SMARTS string.
///
///   ARGUMENTS:
///
///     - SMARTS: the smarts string
///
///     - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached
///       atoms.  So, for example, 'C[H]' becomes '[C;!H0]'.
///       Defaults to 0.
///
///     - replacements: (optional) a dictionary of replacement strings (see below)
///       Defaults to {}. See the documentation for MolFromSmiles for an explanation.
///
///   RETURNS:
///
///     a Mol object, None on failure.
///
/// Alternatively: MolFromSmarts(SMARTS, params) with a SmartsParserParams object.
#[pyfunction]
#[pyo3(name = "MolFromSmarts", signature = (*args, **kwargs))]
fn py_mol_from_smarts(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Option<PyMol>> {
    let get = |pos: usize, key: &str| -> PyResult<Option<Bound<'_, PyAny>>> {
        let kw = match kwargs {
            Some(k) => k.get_item(key)?,
            None => None,
        };
        let p = if pos < args.len() {
            Some(args.get_item(pos)?)
        } else {
            None
        };
        match (p, kw) {
            (Some(_), Some(_)) => Err(PyTypeError::new_err(format!(
                "MolFromSmarts() got multiple values for argument '{key}'"
            ))),
            (a, b) => Ok(a.or(b)),
        }
    };
    let smarts_obj = get(0, "SMARTS")?.ok_or_else(|| {
        PyTypeError::new_err("MolFromSmarts() missing required argument 'SMARTS'")
    })?;
    let smarts = py_object_to_string(&smarts_obj)?;
    if let Some(p) = get(1, "params")? {
        if let Ok(p) = p.cast::<PySmartsParserParams>() {
            let params = p.borrow().to_core();
            return Ok(parse(py, &smarts, &params));
        }
    }
    let merge_hs = match get(1, "mergeHs")? {
        Some(v) => v.is_truthy()?,
        None => false,
    };
    let mut replacements = BTreeMap::new();
    if let Some(d) = get(2, "replacements")? {
        let d = d.cast::<PyDict>()?;
        for (k, v) in d.iter() {
            replacements.insert(k.extract::<String>()?, v.extract::<String>()?);
        }
    }
    let params = SmartsParserParams {
        merge_hs,
        replacements,
        ..SmartsParserParams::default()
    };
    Ok(parse(py, &smarts, &params))
}

/// Construct an atom from a SMILES string
#[pyfunction]
#[pyo3(name = "AtomFromSmiles")]
#[allow(non_snake_case)]
fn py_atom_from_smiles(py: Python<'_>, SMILES: &str) -> PyResult<Option<Py<PyAtom>>> {
    atom_from_smiles(SMILES)
        .map(|a| PyAtom::owned(py, a))
        .transpose()
}

/// Construct a bond from a SMILES string
#[pyfunction]
#[pyo3(name = "BondFromSmiles")]
#[allow(non_snake_case)]
fn py_bond_from_smiles(py: Python<'_>, SMILES: &str) -> PyResult<Option<Py<PyBond>>> {
    bond_from_smiles(SMILES)
        .map(|b| PyBond::owned(py, b))
        .transpose()
}

/// Construct an atom from a SMARTS string
#[pyfunction]
#[pyo3(name = "AtomFromSmarts")]
#[allow(non_snake_case)]
fn py_atom_from_smarts(py: Python<'_>, SMARTS: &str) -> PyResult<Option<Py<PyAtom>>> {
    atom_from_smarts(SMARTS)
        .map(|a| PyAtom::owned(py, a))
        .transpose()
}

/// Construct a bond from a SMARTS string
#[pyfunction]
#[pyo3(name = "BondFromSmarts")]
#[allow(non_snake_case)]
fn py_bond_from_smarts(py: Python<'_>, SMILES: &str) -> PyResult<Option<Py<PyBond>>> {
    bond_from_smarts(SMILES)
        .map(|b| PyBond::owned(py, b))
        .transpose()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySmartsParserParams>()?;
    m.add_function(wrap_pyfunction!(py_mol_from_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_atom_from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_bond_from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_atom_from_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_bond_from_smarts, m)?)?;
    Ok(())
}
