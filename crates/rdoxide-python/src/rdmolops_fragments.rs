//! Bond fragmentation of `rdkit.Chem.rdmolops` (upstream
//! `GraphMol/Wrap/MolOps.cpp`): `FragmentOnBonds`, `FragmentOnSomeBonds` and
//! `FragmentOnBRICSBonds`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use rdoxide::graphmol::bond::BondType;
use rdoxide::graphmol::chem_transforms::{
    fragment_on_bonds, fragment_on_brics_bonds, fragment_on_some_bonds,
};

use crate::enums_impl::PyBondType;
use crate::errors::{Check, IntoPyResult, Site, check_failed};
use crate::mol::PyMol;
use crate::pyconv::python_object_to_vect;
use crate::rdmoldescriptors::bool_arg;

const MOL_FRAGMENTER: &str = "Code/GraphMol/ChemTransforms/MolFragmenter.cpp";

/// The failed `PRECONDITION` of `fragmentOnBonds` / `fragmentOnSomeBonds`
/// that `site_line` (the last line of the upstream statement) belongs to.
fn precondition(line: u32, expr: &'static str, message: &str) -> PyErr {
    check_failed(
        Site {
            check: Check::Precondition,
            file: MOL_FRAGMENTER,
            line,
            expr,
        },
        message,
    )
}

/// The upstream preconditions on the per-bond arguments. `lines` are the
/// statement-end lines of the dummy-label and bond-type checks.
fn check_per_bond_args(
    num_bonds: usize,
    dummy_labels: Option<&Vec<(u32, u32)>>,
    bond_types: Option<&Vec<BondType>>,
    lines: (u32, u32),
) -> PyResult<()> {
    if dummy_labels.is_some_and(|d| d.len() < num_bonds) {
        return Err(precondition(
            lines.0,
            "(!dummyLabels || dummyLabels->size() >= bondIndices.size())",
            "bad dummyLabel vector",
        ));
    }
    if bond_types.is_some_and(|b| b.len() != num_bonds) {
        return Err(precondition(
            lines.1,
            "(!bondTypes || bondTypes->size() == bondIndices.size())",
            "bad bondType vector",
        ));
    }
    Ok(())
}

fn bond_indices(mol: &PyMol, obj: &Bound<'_, PyAny>) -> PyResult<Vec<u32>> {
    python_object_to_vect(Some(obj), Some(mol.mol.num_bonds()))?
        .ok_or_else(|| PyValueError::new_err("empty bond indices"))
}

/// `[(v1, v2), ...]` from a sequence of pairs; `None` unless truthy.
fn dummy_labels(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<(u32, u32)>>> {
    let Some(obj) = obj.filter(|o| o.is_truthy().unwrap_or(false)) else {
        return Ok(None);
    };
    (0..obj.len()?)
        .map(|i| {
            let pair = obj.get_item(i)?;
            Ok((pair.get_item(0)?.extract()?, pair.get_item(1)?.extract()?))
        })
        .collect::<PyResult<Vec<_>>>()
        .map(Some)
}

fn bond_types(obj: Option<&Bound<'_, PyAny>>, num_bonds: usize) -> PyResult<Option<Vec<BondType>>> {
    let Some(obj) = obj.filter(|o| o.is_truthy().unwrap_or(false)) else {
        return Ok(None);
    };
    if obj.len()? != num_bonds {
        return Err(PyValueError::new_err("bondTypes shorter than bondIndices"));
    }
    (0..num_bonds)
        .map(|i| PyBondType::extract(&obj.get_item(i)?))
        .collect::<PyResult<Vec<_>>>()
        .map(Some)
}

/// Return a new molecule with all specified bonds broken
///
///   ARGUMENTS:
///
///       - mol            - the molecule to be modified
///       - bondIndices    - indices of the bonds to be broken
///       - addDummies  - toggles addition of dummy atoms to indicate where
///         bonds were broken
///       - dummyLabels - used to provide the labels to be used for the dummies.
///         the first element in each pair is the label for the dummy
///         that replaces the bond's beginAtom, the second is for the
///         dummy that replaces the bond's endAtom. If not provided, the
///         dummies are labeled with atom indices.
///       - bondTypes - used to provide the bond type to use between the
///         fragments and the dummy atoms. If not provided, defaults to single.
///       - cutsPerAtom - used to return the number of cuts made at each atom.
///
///   RETURNS:
///       a new Mol with the modifications
#[pyfunction]
#[pyo3(
    name = "FragmentOnBonds",
    signature = (mol, bondIndices, addDummies = None, dummyLabels = None, bondTypes = None, cutsPerAtom = None)
)]
#[allow(non_snake_case)]
fn py_fragment_on_bonds(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    bondIndices: &Bound<'_, PyAny>,
    addDummies: Option<&Bound<'_, PyAny>>,
    dummyLabels: Option<&Bound<'_, PyAny>>,
    bondTypes: Option<&Bound<'_, PyAny>>,
    cutsPerAtom: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    let indices = bond_indices(&mol, bondIndices)?;
    let add_dummies = bool_arg(addDummies, true)?;
    let labels = dummy_labels(dummyLabels)?;
    let types = bond_types(bondTypes, indices.len())?;
    // a non-empty list receives the number of cuts per atom
    let cuts_list = cutsPerAtom
        .and_then(|c| c.cast::<PyList>().ok())
        .filter(|l| !l.is_empty());
    let mut cuts = match cuts_list {
        Some(list) => {
            if list.len() < mol.mol.num_atoms() as usize {
                return Err(PyValueError::new_err(
                    "cutsPerAtom shorter than the number of atoms",
                ));
            }
            Some(vec![0u32; list.len()])
        }
        None => None,
    };
    check_per_bond_args(indices.len(), labels.as_ref(), types.as_ref(), (438, 440))?;
    if cuts
        .as_ref()
        .is_some_and(|c| c.len() != mol.mol.num_atoms() as usize)
    {
        return Err(precondition(
            442,
            "(!nCutsPerAtom || nCutsPerAtom->size() == mol.getNumAtoms())",
            "bad nCutsPerAtom vector",
        ));
    }
    let core = &mol.mol;
    let res = py
        .detach(|| {
            fragment_on_bonds(
                core,
                &indices,
                add_dummies,
                labels.as_deref(),
                types.as_deref(),
                cuts.as_deref_mut(),
            )
        })
        .py()?;
    if let (Some(list), Some(cuts)) = (cuts_list, cuts) {
        for (i, c) in cuts.iter().take(core.num_atoms() as usize).enumerate() {
            list.set_item(i, *c)?;
        }
    }
    Ok(PyMol { mol: res })
}

/// fragment on some bonds
#[pyfunction]
#[pyo3(
    name = "FragmentOnSomeBonds",
    signature = (mol, bondIndices, numToBreak = 1, addDummies = None, dummyLabels = None, bondTypes = None, returnCutsPerAtom = None)
)]
#[allow(non_snake_case)]
fn py_fragment_on_some_bonds<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    bondIndices: &Bound<'py, PyAny>,
    numToBreak: u32,
    addDummies: Option<&Bound<'py, PyAny>>,
    dummyLabels: Option<&Bound<'py, PyAny>>,
    bondTypes: Option<&Bound<'py, PyAny>>,
    returnCutsPerAtom: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let indices = bond_indices(&mol, bondIndices)?;
    let add_dummies = bool_arg(addDummies, true)?;
    let labels = dummy_labels(dummyLabels)?;
    let types = bond_types(bondTypes, indices.len())?;
    let want_cuts = bool_arg(returnCutsPerAtom, false)?;
    check_per_bond_args(indices.len(), labels.as_ref(), types.as_ref(), (317, 319))?;
    let core = &mol.mol;
    let (mols, cuts) = py
        .detach(|| {
            fragment_on_some_bonds(
                core,
                &indices,
                numToBreak,
                add_dummies,
                labels.as_deref(),
                types.as_deref(),
                want_cuts,
            )
        })
        .py()?;
    let frags = PyTuple::new(
        py,
        mols.into_iter()
            .map(|mol| Py::new(py, PyMol { mol }))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    if want_cuts {
        let cuts = PyTuple::new(
            py,
            cuts.into_iter()
                .map(|c| PyTuple::new(py, c))
                .collect::<PyResult<Vec<_>>>()?,
        )?;
        PyTuple::new(py, [frags.into_any(), cuts.into_any()])
    } else {
        Ok(frags)
    }
}

/// Return a new molecule with all BRICS bonds broken
#[pyfunction]
#[pyo3(name = "FragmentOnBRICSBonds", signature = (mol))]
fn py_fragment_on_brics_bonds(py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
    let core = &mol.mol;
    let res = py.detach(|| fragment_on_brics_bonds(core)).py()?;
    Ok(PyMol { mol: res })
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_fragment_on_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(py_fragment_on_some_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(py_fragment_on_brics_bonds, m)?)?;
    Ok(())
}
