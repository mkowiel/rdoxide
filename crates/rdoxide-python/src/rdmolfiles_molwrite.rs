//! Mol block writers (upstream `Wrap/rdmolfiles.cpp`): `MolWriterParams`,
//! `MolToMolBlock`, `MolToV3KMolBlock`, `MolToV2KMolBlock`, `MolToMolFile`,
//! `MolToV3KMolFile`.

use pyo3::exceptions::{PyIOError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::file_parsers::mol_file_writer::{
    MolWriterParams, mol_to_mol_block, mol_to_v2k_mol_block, mol_to_v3k_mol_block,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmolfiles_write::Args;

/// Parameters controlling Mol writing
#[pyclass(
    name = "MolWriterParams",
    module = "rdkit.Chem.rdmolfiles",
    get_all,
    set_all,
    skip_from_py_object
)]
#[derive(Clone)]
#[allow(non_snake_case)]
pub struct PyMolWriterParams {
    /// toggles inclusion of stereochemistry information (default=True)
    includeStereo: bool,
    /// triggers kekulization of the molecule before it is written (default=True)
    kekulize: bool,
    /// force generation a V3000 mol block (happens automatically with more than 999 atoms or bonds)(default=False)
    forceV3000: bool,
    /// precision of coordinates (only available in V3000)(default=false)
    precision: u32,
}

#[pymethods]
impl PyMolWriterParams {
    #[new]
    fn new() -> Self {
        let d = MolWriterParams::default();
        PyMolWriterParams {
            includeStereo: d.include_stereo,
            kekulize: d.kekulize,
            forceV3000: d.force_v3000,
            precision: d.precision,
        }
    }
}

impl PyMolWriterParams {
    fn to_core(&self) -> MolWriterParams {
        MolWriterParams {
            include_stereo: self.includeStereo,
            kekulize: self.kekulize,
            force_v3000: self.forceV3000,
            precision: self.precision,
        }
    }
}

#[derive(Clone, Copy)]
enum Kind {
    Default,
    V3K,
    V2K,
}

fn run(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    p: MolWriterParams,
    conf_id: i32,
    kind: Kind,
) -> PyResult<String> {
    let mol = mol.cast::<PyMol>()?;
    let m = mol.borrow();
    let core = &m.mol;
    py.detach(|| match kind {
        Kind::Default => mol_to_mol_block(core, &p, conf_id),
        Kind::V3K => mol_to_v3k_mol_block(core, &p, conf_id),
        Kind::V2K => mol_to_v2k_mol_block(core, &p, conf_id),
    })
    .py()
}

/// Shared argument handling: `(mol, params, confId)` or
/// `(mol, includeStereo, confId, kekulize[, forceV3000])`.
fn dispatch(
    py: Python<'_>,
    name: &'static str,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
    file_arg: bool,
    kind: Kind,
) -> PyResult<(String, Option<String>)> {
    let a = Args { name, args, kwargs };
    let mol = a
        .get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err(format!("{name}() missing required argument 'mol'")))?;
    let off = usize::from(file_arg);
    let file = if file_arg {
        Some(
            a.get(1, "filename")?
                .ok_or_else(|| {
                    PyTypeError::new_err(format!("{name}() missing required argument 'filename'"))
                })?
                .extract::<String>()?,
        )
    } else {
        None
    };
    if let Some(p) = a.get(1 + off, "params")?
        && let Ok(p) = p.cast::<PyMolWriterParams>()
    {
        let params = p.borrow().to_core();
        let conf_id = a.int(2 + off, "confId", -1)?;
        return Ok((run(py, &mol, params, conf_id, kind)?, file));
    }
    let params = MolWriterParams {
        include_stereo: a.flag(1 + off, "includeStereo", true)?,
        kekulize: a.flag(3 + off, "kekulize", true)?,
        force_v3000: if matches!(kind, Kind::Default) {
            a.flag(4 + off, "forceV3000", false)?
        } else {
            false
        },
        ..MolWriterParams::default()
    };
    let conf_id = a.int(2 + off, "confId", -1)?;
    Ok((run(py, &mol, params, conf_id, kind)?, file))
}

fn write_file(name: &str, text: &str) -> PyResult<()> {
    std::fs::write(name, text)
        .map_err(|_| PyIOError::new_err(format!("File error: Bad output file {name}")))
}

/// Returns a Mol block for a molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule
///     - includeStereo: (optional) toggles inclusion of stereochemical
///       information in the output
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - kekulize: (optional) triggers kekulization of the molecule before it's written,
///       as suggested by the MDL spec.
///     - forceV3000 (optional) force generation a V3000 mol block (happens automatically with
///       more than 999 atoms or bonds)
///
///   RETURNS:
///
///     a string
///
/// Alternatively: MolToMolBlock(mol, params, confId=-1) with a MolWriterParams object.
#[pyfunction]
#[pyo3(name = "MolToMolBlock", signature = (*args, **kwargs))]
fn py_mol_to_mol_block(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    Ok(dispatch(py, "MolToMolBlock", args, kwargs, false, Kind::Default)?.0)
}

/// Returns a V3000 Mol block for a molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule
///     - includeStereo: (optional) toggles inclusion of stereochemical
///       information in the output
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - kekulize: (optional) triggers kekulization of the molecule before it's written,
///       as suggested by the MDL spec.
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToV3KMolBlock", signature = (*args, **kwargs))]
fn py_mol_to_v3k_mol_block(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    Ok(dispatch(py, "MolToV3KMolBlock", args, kwargs, false, Kind::V3K)?.0)
}

/// Returns a V2000 Mol block for a molecule
///    ARGUMENTS:
///
///      - mol: the molecule
///      - params: the MolWriterParams
///      - confId: (optional) selects which conformation to output (-1 = default)
///
///    RETURNS:
///
///      a string
///
///    NOTE: this function throws a ValueError if the molecule has more than 999 atoms, bonds, or SGroups
#[pyfunction]
#[pyo3(name = "MolToV2KMolBlock", signature = (mol, params = None, confId = -1))]
#[allow(non_snake_case)]
fn py_mol_to_v2k_mol_block(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    params: Option<PyRef<'_, PyMolWriterParams>>,
    confId: i32,
) -> PyResult<String> {
    let p = params.map_or_else(MolWriterParams::default, |p| p.to_core());
    run(py, mol, p, confId, Kind::V2K)
}

/// Writes a Mol file for a molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule
///     - filename: the file to write to
///     - includeStereo: (optional) toggles inclusion of stereochemical
///       information in the output
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - kekulize: (optional) triggers kekulization of the molecule before it's written,
///       as suggested by the MDL spec.
///     - forceV3000 (optional) force generation a V3000 mol block (happens automatically with
///       more than 999 atoms or bonds)
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToMolFile", signature = (*args, **kwargs))]
fn py_mol_to_mol_file(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let (text, file) = dispatch(py, "MolToMolFile", args, kwargs, true, Kind::Default)?;
    write_file(&file.unwrap_or_default(), &text)
}

/// Writes a V3000 Mol file for a molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule
///     - filename: the file to write to
///     - includeStereo: (optional) toggles inclusion of stereochemical
///       information in the output
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - kekulize: (optional) triggers kekulization of the molecule before it's written,
///       as suggested by the MDL spec.
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToV3KMolFile", signature = (*args, **kwargs))]
fn py_mol_to_v3k_mol_file(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let (text, file) = dispatch(py, "MolToV3KMolFile", args, kwargs, true, Kind::V3K)?;
    write_file(&file.unwrap_or_default(), &text)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMolWriterParams>()?;
    m.add_function(wrap_pyfunction!(py_mol_to_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_v3k_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_v2k_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_mol_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_v3k_mol_file, m)?)?;
    Ok(())
}
