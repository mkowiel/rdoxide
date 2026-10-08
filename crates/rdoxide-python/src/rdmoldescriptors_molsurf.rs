//! `rdkit.Chem.rdMolDescriptors`: Labute ASA, TPSA and the VSA binnings
//! (upstream `Descriptors/MolSurf.cpp` and its wrappers in
//! `Wrap/rdMolDescriptors.cpp`).

use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use rdoxide::graphmol::descriptors::mol_surf::{self as surf, LABUTE_ASA_VERSION, TPSA_VERSION};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;
use crate::vects::PyVectD;

/// `CalcLabuteASA(mol, includeHs=True, force=False)`
#[pyfunction]
#[pyo3(name = "CalcLabuteASA", signature = (mol, includeHs = None, force = None))]
#[allow(non_snake_case)]
fn calc_labute_asa(
    mol: &Bound<'_, PyMol>,
    includeHs: Option<&Bound<'_, PyAny>>,
    force: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let (include_hs, force) = (bool_arg(includeHs, true)?, bool_arg(force, false)?);
    surf::calc_labute_asa(&mut mol.borrow_mut().mol, include_hs, force).py()
}

/// `_CalcLabuteASAContribs(mol, includeHs=True, force=False)`: the atom
/// contributions (a `_vectd`) and the contribution of the implicit Hs.
#[pyfunction]
#[pyo3(name = "_CalcLabuteASAContribs", signature = (mol, includeHs = None, force = None))]
#[allow(non_snake_case)]
fn calc_labute_asa_contribs<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    includeHs: Option<&Bound<'py, PyAny>>,
    force: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (include_hs, force) = (bool_arg(includeHs, true)?, bool_arg(force, false)?);
    let c = surf::get_labute_atom_contribs(&mut mol.borrow_mut().mol, include_hs, force).py()?;
    PyTuple::new(
        py,
        [
            PyVectD::new_py(py, c.atoms)?,
            c.h.into_pyobject(py)?.into_any().unbind(),
        ],
    )
}

/// `CalcTPSA(mol, force=False, includeSandP=False)`
#[pyfunction]
#[pyo3(name = "CalcTPSA", signature = (mol, force = None, includeSandP = None))]
#[allow(non_snake_case)]
fn calc_tpsa(
    mol: &Bound<'_, PyMol>,
    force: Option<&Bound<'_, PyAny>>,
    includeSandP: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let (force, s_and_p) = (bool_arg(force, false)?, bool_arg(includeSandP, false)?);
    surf::calc_tpsa(&mut mol.borrow_mut().mol, force, s_and_p).py()
}

/// `_CalcTPSAContribs(mol, force=False, includeSandP=False)`: tuple of the atom contributions
#[pyfunction]
#[pyo3(name = "_CalcTPSAContribs", signature = (mol, force = None, includeSandP = None))]
#[allow(non_snake_case)]
fn calc_tpsa_contribs<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    force: Option<&Bound<'py, PyAny>>,
    includeSandP: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (force, s_and_p) = (bool_arg(force, false)?, bool_arg(includeSandP, false)?);
    let (contribs, _) =
        surf::get_tpsa_atom_contribs(&mut mol.borrow_mut().mol, force, s_and_p).py()?;
    PyTuple::new(py, contribs)
}

/// The bins of `SlogP_VSA_`/`SMR_VSA_`: `None` (use the defaults) for a missing
/// or empty sequence, as `if (bins) { if (len(bins)) ... }` upstream.
fn optional_bins(bins: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<f64>>> {
    match bins {
        Some(b) if b.is_truthy()? => Ok(Some(all_bins(b)?).filter(|v| !v.is_empty())),
        _ => Ok(None),
    }
}

/// `len(bins)` doubles (`TypeError` for objects without a length).
fn all_bins(bins: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    (0..bins.len()?)
        .map(|i| bins.get_item(i)?.extract::<f64>())
        .collect()
}

fn to_list<'py>(py: Python<'py>, values: Vec<f64>) -> PyResult<Bound<'py, PyList>> {
    PyList::new(py, values)
}

/// `SlogP_VSA_(mol, bins=[], force=False)`
#[pyfunction]
#[pyo3(name = "SlogP_VSA_", signature = (mol, bins = None, force = None))]
fn slogp_vsa<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    bins: Option<&Bound<'py, PyAny>>,
    force: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyList>> {
    let (bins, force) = (optional_bins(bins)?, bool_arg(force, false)?);
    to_list(
        py,
        surf::calc_slogp_vsa(&mut mol.borrow_mut().mol, bins.as_deref(), force).py()?,
    )
}

/// `SMR_VSA_(mol, bins=[], force=False)`
#[pyfunction]
#[pyo3(name = "SMR_VSA_", signature = (mol, bins = None, force = None))]
fn smr_vsa<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    bins: Option<&Bound<'py, PyAny>>,
    force: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyList>> {
    let (bins, force) = (optional_bins(bins)?, bool_arg(force, false)?);
    to_list(
        py,
        surf::calc_smr_vsa(&mut mol.borrow_mut().mol, bins.as_deref(), force).py()?,
    )
}

/// `PEOE_VSA_(mol, bins=[], force=False)`
#[pyfunction]
#[pyo3(name = "PEOE_VSA_", signature = (mol, bins = None, force = None))]
fn peoe_vsa<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    bins: Option<&Bound<'py, PyAny>>,
    force: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyList>> {
    let (bins, force) = (optional_bins(bins)?, bool_arg(force, false)?);
    to_list(
        py,
        surf::calc_peoe_vsa(&mut mol.borrow_mut().mol, bins.as_deref(), force).py()?,
    )
}

/// `CustomProp_VSA_(mol, customPropName, bins, force=False)`
#[pyfunction]
#[pyo3(name = "CustomProp_VSA_", signature = (mol, customPropName, bins, force = None))]
#[allow(non_snake_case)]
fn custom_prop_vsa<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    customPropName: &str,
    bins: &Bound<'py, PyAny>,
    force: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyList>> {
    let (bins, force) = (all_bins(bins)?, bool_arg(force, false)?);
    to_list(
        py,
        surf::calc_custom_prop_vsa(&mut mol.borrow_mut().mol, customPropName, &bins, force).py()?,
    )
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_labute_asa, m)?)?;
    m.add_function(wrap_pyfunction!(calc_labute_asa_contribs, m)?)?;
    m.add_function(wrap_pyfunction!(calc_tpsa, m)?)?;
    m.add_function(wrap_pyfunction!(calc_tpsa_contribs, m)?)?;
    m.add_function(wrap_pyfunction!(slogp_vsa, m)?)?;
    m.add_function(wrap_pyfunction!(smr_vsa, m)?)?;
    m.add_function(wrap_pyfunction!(peoe_vsa, m)?)?;
    m.add_function(wrap_pyfunction!(custom_prop_vsa, m)?)?;
    m.add("_CalcLabuteASA_version", LABUTE_ASA_VERSION)?;
    m.add("_CalcTPSA_version", TPSA_VERSION)?;
    Ok(())
}
