//! `rdoxide.Avalon.pyAvalonTools`: Avalon fingerprints, canonical SMILES and 2D coordinates (upstream
//! `External/AvalonTools/nbWrap/pyAvalonTools.cpp`). The structure checker (`CheckMolecule`, ...) is
//! not available yet.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rdoxide::external::avalon_tools::{
    DEFAULT_BIT_FLAGS, get_avalon_count_fp, get_avalon_count_fp_from_data, get_avalon_fp,
    get_avalon_fp_as_words, get_avalon_fp_as_words_from_data, get_avalon_fp_from_data,
    get_canon_smiles, get_canon_smiles_from_data, set_2d_coords, set_2d_coords_from_data,
};

use crate::datastructs::PyExplicitBitVect;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

type Args<'a, 'py> = &'a Bound<'py, pyo3::types::PyTuple>;
type Kwargs<'a, 'py> = Option<&'a Bound<'py, pyo3::types::PyDict>>;

/// the argument called `name` at position `pos` of the overload, or its keyword, or `default`
fn arg<'py, T: for<'a> FromPyObject<'a, 'py>>(
    args: Args<'_, 'py>,
    kwargs: Kwargs<'_, 'py>,
    pos: usize,
    name: &str,
    default: Option<T>,
) -> PyResult<T>
where
    for<'a> <T as FromPyObject<'a, 'py>>::Error: Into<PyErr>,
{
    let v = match args.get_item(pos) {
        Ok(v) => Some(v),
        Err(_) => match kwargs {
            Some(d) => d.get_item(name)?,
            None => None,
        },
    };
    match v {
        Some(v) => v.extract::<T>().map_err(Into::into),
        None => default.ok_or_else(|| {
            pyo3::exceptions::PyTypeError::new_err(format!("missing required argument '{name}'"))
        }),
    }
}

fn no_molecule() -> PyErr {
    PyRuntimeError::new_err("could not build a molecule from the data")
}

/// returns an Avalon fingerprint for a molecule or for molecule data (SMILES or MDL mol block)
///
/// The calling conventions are `GetAvalonFP(mol, nBits=512, isQuery=False, resetVect=False,
/// bitFlags=avalonSimilarityBits)` and `GetAvalonFP(molData, isSmiles, nBits=512, ...)`.
#[pyfunction]
#[pyo3(name = "GetAvalonFP", signature = (mol, *args, **kwargs))]
fn get_avalon_fp_py(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    args: &Bound<'_, pyo3::types::PyTuple>,
    kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
) -> PyResult<PyExplicitBitVect> {
    if let Ok(data) = mol.extract::<String>() {
        let is_smiles: bool = arg(args, kwargs, 0, "isSmiles", None)?;
        let n_bits: u32 = arg(args, kwargs, 1, "nBits", Some(512))?;
        let is_query: bool = arg(args, kwargs, 2, "isQuery", Some(false))?;
        let bit_flags: i32 = arg(args, kwargs, 4, "bitFlags", Some(DEFAULT_BIT_FLAGS))?;
        let res = py
            .detach(|| get_avalon_fp_from_data(&data, is_smiles, n_bits, is_query, bit_flags))
            .py()?
            .ok_or_else(no_molecule)?;
        return Ok(PyExplicitBitVect::wrap(res));
    }
    let mol: PyRef<'_, PyMol> = mol.extract()?;
    let n_bits: u32 = arg(args, kwargs, 0, "nBits", Some(512))?;
    let is_query: bool = arg(args, kwargs, 1, "isQuery", Some(false))?;
    let bit_flags: i32 = arg(args, kwargs, 3, "bitFlags", Some(DEFAULT_BIT_FLAGS))?;
    let m = &mol.mol;
    let res = py
        .detach(|| get_avalon_fp(m, n_bits, is_query, bit_flags))
        .py()?;
    Ok(PyExplicitBitVect::wrap(res))
}

/// returns the Avalon fingerprint of a molecule (or molecule data) as a list of 32 bit words
#[pyfunction]
#[pyo3(name = "GetAvalonFPAsWords", signature = (mol, *args, **kwargs))]
fn get_avalon_fp_as_words_py(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    args: &Bound<'_, pyo3::types::PyTuple>,
    kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
) -> PyResult<Vec<u32>> {
    if let Ok(data) = mol.extract::<String>() {
        let is_smiles: bool = arg(args, kwargs, 0, "isSmiles", None)?;
        let n_bits: u32 = arg(args, kwargs, 1, "nBits", Some(512))?;
        let is_query: bool = arg(args, kwargs, 2, "isQuery", Some(false))?;
        let bit_flags: i32 = arg(args, kwargs, 4, "bitFlags", Some(DEFAULT_BIT_FLAGS))?;
        return py
            .detach(|| {
                get_avalon_fp_as_words_from_data(&data, is_smiles, n_bits, is_query, bit_flags)
            })
            .py()?
            .ok_or_else(no_molecule);
    }
    let mol: PyRef<'_, PyMol> = mol.extract()?;
    let n_bits: u32 = arg(args, kwargs, 0, "nBits", Some(512))?;
    let is_query: bool = arg(args, kwargs, 1, "isQuery", Some(false))?;
    let bit_flags: i32 = arg(args, kwargs, 3, "bitFlags", Some(DEFAULT_BIT_FLAGS))?;
    let m = &mol.mol;
    py.detach(|| get_avalon_fp_as_words(m, n_bits, is_query, bit_flags))
        .py()
}

/// returns the Avalon count fingerprint of a molecule (or molecule data)
#[pyfunction]
#[pyo3(name = "GetAvalonCountFP", signature = (mol, *args, **kwargs))]
fn get_avalon_count_fp_py(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    args: &Bound<'_, pyo3::types::PyTuple>,
    kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
) -> PyResult<Vec<u32>> {
    if let Ok(data) = mol.extract::<String>() {
        let is_smiles: bool = arg(args, kwargs, 0, "isSmiles", None)?;
        let n_bits: u32 = arg(args, kwargs, 1, "nBits", Some(512))?;
        let is_query: bool = arg(args, kwargs, 2, "isQuery", Some(false))?;
        let bit_flags: i32 = arg(args, kwargs, 3, "bitFlags", Some(DEFAULT_BIT_FLAGS))?;
        return py
            .detach(|| get_avalon_count_fp_from_data(&data, is_smiles, n_bits, is_query, bit_flags))
            .py()?
            .ok_or_else(no_molecule);
    }
    let mol: PyRef<'_, PyMol> = mol.extract()?;
    let n_bits: u32 = arg(args, kwargs, 0, "nBits", Some(512))?;
    let is_query: bool = arg(args, kwargs, 1, "isQuery", Some(false))?;
    let bit_flags: i32 = arg(args, kwargs, 2, "bitFlags", Some(DEFAULT_BIT_FLAGS))?;
    let m = &mol.mol;
    py.detach(|| get_avalon_count_fp(m, n_bits, is_query, bit_flags))
        .py()
}

/// returns the canonical SMILES of a molecule, or of molecule data (`GetCanonSmiles(molData, isSmiles, flags=-1)`)
#[pyfunction]
#[pyo3(name = "GetCanonSmiles", signature = (mol, *args, **kwargs))]
fn get_canon_smiles_py(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    args: &Bound<'_, pyo3::types::PyTuple>,
    kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
) -> PyResult<String> {
    if let Ok(data) = mol.extract::<String>() {
        let is_smiles: bool = arg(args, kwargs, 0, "isSmiles", None)?;
        let flags: i32 = arg(args, kwargs, 1, "flags", Some(-1))?;
        return Ok(py.detach(|| get_canon_smiles_from_data(&data, is_smiles, flags)));
    }
    let mol: PyRef<'_, PyMol> = mol.extract()?;
    let flags: i32 = arg(args, kwargs, 0, "flags", Some(-1))?;
    let m = &mol.mol;
    py.detach(|| get_canon_smiles(m, flags)).py()
}

/// generates 2D coordinates with the Avalon layout: `Generate2DCoords(mol, clearConfs=True)` returns the
/// conformer id, `Generate2DCoords(molData, isSmiles)` returns an MDL mol block
#[pyfunction]
#[pyo3(name = "Generate2DCoords", signature = (mol, *args, **kwargs))]
fn generate_2d_coords_py(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    args: &Bound<'_, pyo3::types::PyTuple>,
    kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
) -> PyResult<Py<PyAny>> {
    if let Ok(data) = mol.extract::<String>() {
        let is_smiles: bool = arg(args, kwargs, 0, "isSmiles", None)?;
        let mb = py.detach(|| set_2d_coords_from_data(&data, is_smiles));
        return Ok(mb.into_pyobject(py)?.into_any().unbind());
    }
    let mut mol: PyRefMut<'_, PyMol> = mol.extract()?;
    let clear: bool = arg(args, kwargs, 0, "clearConfs", Some(true))?;
    let id = set_2d_coords(&mut mol.mol, clear).py()?;
    Ok(id.into_pyobject(py)?.into_any().unbind())
}

/// the parts of AvalonTools that are not ported raise an `ImportError`, so that code guarding
/// `import MolKey` with `except ImportError` keeps treating AvalonTools as unavailable (dunder
/// names still raise `AttributeError` for the import machinery and introspection)
const MISSING_ATTRIBUTES: &std::ffi::CStr = c"
def __getattr__(name):
    if name.startswith('__'):
        raise AttributeError(name)
    raise ImportError(f'pyAvalonTools.{name} is not available in rdoxide yet')
";

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let globals = pyo3::types::PyDict::new(py);
    py.run(MISSING_ATTRIBUTES, Some(&globals), None)?;
    m.setattr("__getattr__", globals.get_item("__getattr__")?)?;
    m.add_function(wrap_pyfunction!(get_avalon_fp_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_avalon_fp_as_words_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_avalon_count_fp_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_canon_smiles_py, m)?)?;
    m.add_function(wrap_pyfunction!(generate_2d_coords_py, m)?)?;
    m.add("avalonSSSBits", 0x007FFF)?;
    m.add("avalonSimilarityBits", 0xF07FFF)?;
    Ok(())
}
