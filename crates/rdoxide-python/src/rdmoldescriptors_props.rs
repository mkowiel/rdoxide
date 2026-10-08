//! `rdkit.Chem.rdMolDescriptors`: molecular weights, formula and atom counts
//! (upstream `Descriptors/MolDescriptors.cpp` and its wrappers in
//! `Wrap/rdMolDescriptors.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::mol_props;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// `_CalcMolWt(mol, onlyHeavy=False)`: average molecular weight.
#[pyfunction]
#[pyo3(name = "_CalcMolWt", signature = (mol, onlyHeavy = None))]
#[allow(non_snake_case)]
fn calc_mol_wt(mol: &Bound<'_, PyMol>, onlyHeavy: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
    let only_heavy = bool_arg(onlyHeavy, false)?;
    mol_props::get_avg_mol_wt(&mol.borrow().mol, only_heavy).py()
}

/// `CalcExactMolWt(mol, onlyHeavy=False)`: monoisotopic molecular weight.
#[pyfunction]
#[pyo3(name = "CalcExactMolWt", signature = (mol, onlyHeavy = None))]
#[allow(non_snake_case)]
fn calc_exact_mol_wt(
    mol: &Bound<'_, PyMol>,
    onlyHeavy: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let only_heavy = bool_arg(onlyHeavy, false)?;
    mol_props::get_exact_mol_wt(&mol.borrow().mol, only_heavy).py()
}

/// `CalcMolFormula(mol, separateIsotopes=False, abbreviateHIsotopes=True)`
#[pyfunction]
#[pyo3(
    name = "CalcMolFormula",
    signature = (mol, separateIsotopes = None, abbreviateHIsotopes = None)
)]
#[allow(non_snake_case)]
fn calc_mol_formula(
    mol: &Bound<'_, PyMol>,
    separateIsotopes: Option<&Bound<'_, PyAny>>,
    abbreviateHIsotopes: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let separate = bool_arg(separateIsotopes, false)?;
    let abbreviate = bool_arg(abbreviateHIsotopes, true)?;
    mol_props::get_mol_formula(&mol.borrow().mol, separate, abbreviate).py()
}

/// `CalcNumHeavyAtoms(mol)`
#[pyfunction]
#[pyo3(name = "CalcNumHeavyAtoms", signature = (mol))]
fn calc_num_heavy_atoms(mol: &Bound<'_, PyMol>) -> u32 {
    mol.borrow().mol.num_heavy_atoms()
}

/// `CalcNumAtoms(mol)`: all atoms, implicit Hs included.
#[pyfunction]
#[pyo3(name = "CalcNumAtoms", signature = (mol))]
fn calc_num_atoms(mol: &Bound<'_, PyMol>) -> PyResult<u32> {
    mol.borrow().mol.num_atoms_with_hs().py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_mol_wt, m)?)?;
    m.add_function(wrap_pyfunction!(calc_exact_mol_wt, m)?)?;
    m.add_function(wrap_pyfunction!(calc_mol_formula, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_heavy_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_atoms, m)?)?;
    m.add("_CalcMolWt_version", "1.0.0")?;
    m.add("_CalcExactMolWt_version", "1.0.0")?;
    m.add("_CalcMolFormula_version", "1.3.0")?;
    m.add("_CalcNumHeavyAtoms_version", "1.0.0")?;
    m.add("_CalcNumAtoms_version", "1.0.0")?;
    Ok(())
}
