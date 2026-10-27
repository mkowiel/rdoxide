//! `rdoxide.Avalon.pyAvalonTools`: Avalon fingerprints (upstream
//! `External/AvalonTools/Wrap/pyAvalonTools.cpp`). Functions that need the Avalon structure
//! check, SMILES reader or 2D layout are not available yet.

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use rdoxide::external::avalon_tools::{
    DEFAULT_BIT_FLAGS, get_avalon_count_fp, get_avalon_fp, get_avalon_fp_as_words,
    get_avalon_fp_from_molblock,
};

use crate::datastructs::PyExplicitBitVect;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// returns an Avalon fingerprint for a molecule (or a MOL block when `isSmiles` is False)
///
/// The calling conventions are `GetAvalonFP(mol, nBits=512, isQuery=False, resetVect=False,
/// bitFlags=15761407)` and `GetAvalonFP(molBlock, isSmiles=False, nBits=512, ...)`.
#[pyfunction]
#[pyo3(name = "GetAvalonFP", signature = (mol, *args, **kwargs))]
fn get_avalon_fp_py(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    args: &Bound<'_, pyo3::types::PyTuple>,
    kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
) -> PyResult<PyExplicitBitVect> {
    let kw = |name: &str| -> PyResult<Option<Bound<'_, PyAny>>> {
        match kwargs {
            Some(d) => d.get_item(name),
            None => Ok(None),
        }
    };
    if let Ok(data) = mol.extract::<String>() {
        // (data, isSmiles, nBits, isQuery, resetVect, bitFlags)
        let pos = |i: usize| args.get_item(i).ok();
        let is_smiles: bool = match pos(0).or(kw("isSmiles")?) {
            Some(v) => v.extract()?,
            None => true,
        };
        let n_bits: u32 = match pos(1).or(kw("nBits")?) {
            Some(v) => v.extract()?,
            None => 512,
        };
        let is_query: bool = match pos(2).or(kw("isQuery")?) {
            Some(v) => v.extract()?,
            None => false,
        };
        let bit_flags: i32 = match pos(4).or(kw("bitFlags")?) {
            Some(v) => v.extract()?,
            None => DEFAULT_BIT_FLAGS,
        };
        if is_smiles {
            return Err(PyNotImplementedError::new_err(
                "SMILES input needs the Avalon SMILES reader and layout, which are not ported yet",
            ));
        }
        let res = py
            .detach(|| get_avalon_fp_from_molblock(&data, n_bits, is_query, bit_flags))
            .py()?;
        return Ok(PyExplicitBitVect::wrap(res));
    }
    let mol: PyRef<'_, PyMol> = mol.extract()?;
    let pos = |i: usize| args.get_item(i).ok();
    let n_bits: u32 = match pos(0).or(kw("nBits")?) {
        Some(v) => v.extract()?,
        None => 512,
    };
    let is_query: bool = match pos(1).or(kw("isQuery")?) {
        Some(v) => v.extract()?,
        None => false,
    };
    let bit_flags: i32 = match pos(3).or(kw("bitFlags")?) {
        Some(v) => v.extract()?,
        None => DEFAULT_BIT_FLAGS,
    };
    let m = &mol.mol;
    let res = py
        .detach(|| get_avalon_fp(m, n_bits, is_query, bit_flags))
        .py()?;
    Ok(PyExplicitBitVect::wrap(res))
}

/// returns the Avalon fingerprint of a molecule as a list of 32 bit words
#[pyfunction]
#[pyo3(name = "GetAvalonFPAsWords", signature = (mol, nBits = 512, isQuery = false, resetVect = false, bitFlags = DEFAULT_BIT_FLAGS))]
#[allow(non_snake_case)]
fn get_avalon_fp_as_words_py(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    nBits: u32,
    isQuery: bool,
    resetVect: bool,
    bitFlags: i32,
) -> PyResult<Vec<u32>> {
    let _ = resetVect;
    let m = &mol.mol;
    py.detach(|| get_avalon_fp_as_words(m, nBits, isQuery, bitFlags))
        .py()
}

/// returns the Avalon count fingerprint of a molecule
#[pyfunction]
#[pyo3(name = "GetAvalonCountFP", signature = (mol, nBits = 512, isQuery = false, bitFlags = DEFAULT_BIT_FLAGS))]
#[allow(non_snake_case)]
fn get_avalon_count_fp_py(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    nBits: u32,
    isQuery: bool,
    bitFlags: i32,
) -> PyResult<Vec<u32>> {
    let m = &mol.mol;
    py.detach(|| get_avalon_count_fp(m, nBits, isQuery, bitFlags))
        .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_avalon_fp_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_avalon_fp_as_words_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_avalon_count_fp_py, m)?)?;
    Ok(())
}
