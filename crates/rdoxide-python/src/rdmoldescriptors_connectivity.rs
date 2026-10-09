//! `rdkit.Chem.rdMolDescriptors`: Chi, Hall-Kier alpha, Kappa and Phi
//! (upstream `Descriptors/ConnectivityDescriptors.cpp` and its wrappers in
//! `Wrap/rdMolDescriptors.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::graphmol::descriptors::connectivity::{
    self as conn, CHI_VERSION, HALL_KIER_ALPHA_VERSION, KAPPA_VERSION, PHI_VERSION,
};

use crate::errors::{IntoPyResult, expect_list};
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// Defines `Calc<name>(mol, force=False)` over a `&mut ROMol` (the Chi values
/// are cached on the molecule).
macro_rules! chi_fn {
    ($fname:ident, $pyname:literal, $core:path) => {
        #[pyfunction]
        #[pyo3(name = $pyname, signature = (mol, force = None))]
        fn $fname(mol: &Bound<'_, PyMol>, force: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
            let force = bool_arg(force, false)?;
            $core(&mut mol.borrow_mut().mol, force).py()
        }
    };
}

chi_fn!(calc_chi0v, "CalcChi0v", conn::calc_chi0v);
chi_fn!(calc_chi1v, "CalcChi1v", conn::calc_chi1v);
chi_fn!(calc_chi0n, "CalcChi0n", conn::calc_chi0n);
chi_fn!(calc_chi1n, "CalcChi1n", conn::calc_chi1n);

/// Defines `Calc<name>(mol, force=False)` for the path-based Chi of fixed order.
macro_rules! chi_order_fn {
    ($fname:ident, $pyname:literal, $core:path, $order:literal) => {
        #[pyfunction]
        #[pyo3(name = $pyname, signature = (mol, force = None))]
        fn $fname(mol: &Bound<'_, PyMol>, force: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
            let force = bool_arg(force, false)?;
            $core(&mut mol.borrow_mut().mol, $order, force).py()
        }
    };
}

chi_order_fn!(calc_chi2v, "CalcChi2v", conn::calc_chi_nv, 2);
chi_order_fn!(calc_chi3v, "CalcChi3v", conn::calc_chi_nv, 3);
chi_order_fn!(calc_chi4v, "CalcChi4v", conn::calc_chi_nv, 4);
chi_order_fn!(calc_chi2n, "CalcChi2n", conn::calc_chi_nn, 2);
chi_order_fn!(calc_chi3n, "CalcChi3n", conn::calc_chi_nn, 3);
chi_order_fn!(calc_chi4n, "CalcChi4n", conn::calc_chi_nn, 4);

/// `CalcChiNv(mol, n, force=False)`: Chi of order `n` using valence deltas.
#[pyfunction]
#[pyo3(name = "CalcChiNv", signature = (mol, n, force = None))]
fn calc_chi_nv(mol: &Bound<'_, PyMol>, n: u32, force: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
    let force = bool_arg(force, false)?;
    conn::calc_chi_nv(&mut mol.borrow_mut().mol, n, force).py()
}

/// `CalcChiNn(mol, n, force=False)`: Chi of order `n` using the number of outer electrons.
#[pyfunction]
#[pyo3(name = "CalcChiNn", signature = (mol, n, force = None))]
fn calc_chi_nn(mol: &Bound<'_, PyMol>, n: u32, force: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
    let force = bool_arg(force, false)?;
    conn::calc_chi_nn(&mut mol.borrow_mut().mol, n, force).py()
}

/// `CalcHallKierAlpha(mol, atomContribs=None)`: with a list of one entry per
/// atom, the per-atom alphas are stored into it.
#[pyfunction]
#[pyo3(name = "CalcHallKierAlpha", signature = (mol, atomContribs = None))]
#[allow(non_snake_case)]
fn calc_hall_kier_alpha(
    mol: &Bound<'_, PyMol>,
    atomContribs: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let mol = mol.borrow();
    let n_atoms = mol.mol.num_atoms() as usize;
    // `None` means "do not report the contributions"
    let list = match atomContribs {
        Some(obj) if !obj.is_none() => Some(expect_list(obj)?),
        _ => None,
    };
    if let Some(list) = list {
        if list.len() != n_atoms {
            return Err(PyValueError::new_err(
                "length of atomContribs list != number of atoms",
            ));
        }
        let mut contribs = vec![0.0; n_atoms];
        let res = conn::calc_hall_kier_alpha(&mol.mol, Some(&mut contribs)).py()?;
        for (i, c) in contribs.into_iter().enumerate() {
            list.set_item(i, c)?;
        }
        return Ok(res);
    }
    conn::calc_hall_kier_alpha(&mol.mol, None).py()
}

/// Defines `Calc<name>(mol)` for the descriptors that only read the molecule.
macro_rules! mol_fn {
    ($fname:ident, $pyname:literal, $core:path) => {
        #[pyfunction]
        #[pyo3(name = $pyname, signature = (mol))]
        fn $fname(mol: &Bound<'_, PyMol>) -> PyResult<f64> {
            $core(&mol.borrow().mol).py()
        }
    };
}

mol_fn!(calc_kappa1, "CalcKappa1", conn::calc_kappa1);
mol_fn!(calc_kappa2, "CalcKappa2", conn::calc_kappa2);
mol_fn!(calc_kappa3, "CalcKappa3", conn::calc_kappa3);
mol_fn!(calc_phi, "CalcPhi", conn::calc_phi);

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_chi0v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi1v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi2v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi3v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi4v, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi0n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi1n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi2n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi3n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi4n, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_nv, m)?)?;
    m.add_function(wrap_pyfunction!(calc_chi_nn, m)?)?;
    m.add_function(wrap_pyfunction!(calc_hall_kier_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(calc_kappa1, m)?)?;
    m.add_function(wrap_pyfunction!(calc_kappa2, m)?)?;
    m.add_function(wrap_pyfunction!(calc_kappa3, m)?)?;
    m.add_function(wrap_pyfunction!(calc_phi, m)?)?;
    for name in [
        "Chi0v", "Chi1v", "Chi2v", "Chi3v", "Chi4v", "ChiNv", "Chi0n", "Chi1n", "Chi2n", "Chi3n",
        "Chi4n", "ChiNn",
    ] {
        m.add(format!("_Calc{name}_version"), CHI_VERSION)?;
    }
    m.add("_CalcHallKierAlpha_version", HALL_KIER_ALPHA_VERSION)?;
    for name in ["Kappa1", "Kappa2", "Kappa3"] {
        m.add(format!("_Calc{name}_version"), KAPPA_VERSION)?;
    }
    m.add("_CalcPhi_version", PHI_VERSION)?;
    Ok(())
}
