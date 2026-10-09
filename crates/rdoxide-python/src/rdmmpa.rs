//! `rdkit.Chem.rdMMPA` (upstream `GraphMol/MMPA/Wrap/rdMMPA.cpp`).

use pyo3::exceptions::{PyAttributeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::mmpa::{DEFAULT_PATTERN, MmpaPair, fragment_mol, fragment_mol_at_bonds};
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmolfiles_write::Args;

fn smiles_of(mol: &ROMol) -> PyResult<String> {
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        ..SmilesWriteParams::default()
    };
    mol_to_smiles(&mut mol.clone(), &params).py()
}

/// A tuple of (core, side chains) pairs, as molecules or SMILES strings.
fn pairs_to_py<'py>(
    py: Python<'py>,
    pairs: Vec<MmpaPair>,
    as_mols: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    let mut out = Vec::with_capacity(pairs.len());
    for (core, side_chains) in pairs {
        let item = if as_mols {
            let core = match core {
                Some(mol) => Py::new(py, PyMol { mol })?.into_any(),
                None => py.None(),
            };
            let side_chains = Py::new(py, PyMol { mol: side_chains })?.into_any();
            PyTuple::new(py, [core, side_chains])?
        } else {
            let core = core
                .as_ref()
                .map(smiles_of)
                .transpose()?
                .unwrap_or_default();
            PyTuple::new(py, [core, smiles_of(&side_chains)?])?
        };
        out.push(item);
    }
    PyTuple::new(py, out)
}

fn count(args: &Args<'_, '_>, pos: usize, key: &str, default: usize) -> PyResult<usize> {
    args.get(pos, key)?.map_or(Ok(default), |v| v.extract())
}

fn bonds_from(obj: &Bound<'_, PyAny>) -> PyResult<Vec<u32>> {
    if obj.try_iter().is_err() {
        return Err(PyAttributeError::new_err(format!(
            "'{}' object has no attribute '__iter__'",
            obj.get_type().name()?
        )));
    }
    obj.extract()
}

fn pattern(args: &Args<'_, '_>, pos: usize) -> PyResult<String> {
    args.get(pos, "pattern")?
        .map_or_else(|| Ok(DEFAULT_PATTERN.to_string()), |v| v.extract())
}

/// Does the fragmentation necessary for an MMPA analysis
///
/// Overloads (as in the C++ wrapper):
///   FragmentMol(mol, maxCuts=3, maxCutBonds=20, pattern=..., resultsAsMols=True)
///   FragmentMol(mol, minCuts, maxCuts, maxCutBonds, pattern=..., resultsAsMols=True)
///   FragmentMol(mol, bondsToCut, minCuts=1, maxCuts=3, resultsAsMols=True)
#[pyfunction]
#[pyo3(name = "FragmentMol", signature = (*args, **kwargs))]
fn py_fragment_mol<'py>(
    py: Python<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let a = Args {
        name: "FragmentMol",
        args,
        kwargs,
    };
    let mol = a
        .get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err("FragmentMol() missing required argument 'mol'"))?;
    let mol = mol.cast::<PyMol>()?.borrow().mol.clone();
    // Boost.Python tries the overloads last-registered first and the
    // `bondsToCut` one takes any object, so a second positional argument
    // always selects it (and a non-sequence then fails inside).
    let bonds_to_cut: Option<Vec<u32>> = if args.len() == 2 || (3..=5).contains(&args.len()) {
        Some(bonds_from(&args.get_item(1)?)?)
    } else {
        a.get(1, "bondsToCut")?
            .filter(|_| args.len() < 2)
            .map(|v| bonds_from(&v))
            .transpose()?
    };

    if let Some(bonds) = bonds_to_cut {
        if bonds.is_empty() {
            return Err(PyValueError::new_err("bondsToCut must be non-empty"));
        }
        let min_cuts = count(&a, 2, "minCuts", 1)?;
        let max_cuts = count(&a, 3, "maxCuts", 3)?;
        let as_mols = a.flag(4, "resultsAsMols", true)?;
        let res = py
            .detach(|| fragment_mol_at_bonds(&mol, &bonds, min_cuts, max_cuts))
            .py()?;
        return pairs_to_py(py, res, as_mols);
    }
    let by_min_cuts = a.get(1, "minCuts")?.is_some() || args.len() > 5;
    let (min_cuts, max_cuts, max_cut_bonds, pat, as_mols) = if by_min_cuts {
        let min_cuts = a
            .get(1, "minCuts")?
            .ok_or_else(|| PyTypeError::new_err("FragmentMol() missing argument 'minCuts'"))?
            .extract()?;
        (
            min_cuts,
            a.get(2, "maxCuts")?
                .ok_or_else(|| PyTypeError::new_err("FragmentMol() missing argument 'maxCuts'"))?
                .extract()?,
            a.get(3, "maxCutBonds")?
                .ok_or_else(|| {
                    PyTypeError::new_err("FragmentMol() missing argument 'maxCutBonds'")
                })?
                .extract()?,
            pattern(&a, 4)?,
            a.flag(5, "resultsAsMols", true)?,
        )
    } else {
        (
            1,
            count(&a, 1, "maxCuts", 3)?,
            count(&a, 2, "maxCutBonds", 20)?,
            pattern(&a, 3)?,
            a.flag(4, "resultsAsMols", true)?,
        )
    };
    let res = py
        .detach(|| fragment_mol(&mol, min_cuts, max_cuts, max_cut_bonds, &pat))
        .py()?
        .unwrap_or_default();
    pairs_to_py(py, res, as_mols)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing a C++ implementation of code for doing MMPA",
    )?;
    m.add_function(wrap_pyfunction!(py_fragment_mol, m)?)
}
