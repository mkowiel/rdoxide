//! Module-level atom helpers of `rdoxide.Chem.rdchem` (upstream
//! `GraphMol/Wrap/Atom.cpp`): MDL R labels, aliases and values, the
//! supplemental SMILES label and `GetNumPiElectrons`.

use pyo3::prelude::*;
use rdoxide::graphmol::atom::{
    get_atom_alias, get_atom_rlabel, get_atom_value, get_supplemental_smiles_label,
    num_pi_electrons, set_atom_alias, set_atom_rlabel, set_atom_value,
    set_supplemental_smiles_label,
};

use crate::atom::PyAtom;
use crate::errors::IntoPyResult;

/// Returns the atom's MDL AtomRLabel (this is an integer from 0 to 99)
#[pyfunction]
#[pyo3(name = "GetAtomRLabel", signature = (atom))]
fn py_get_atom_rlabel(py: Python<'_>, atom: PyRef<'_, PyAtom>) -> PyResult<i32> {
    atom.with(py, get_atom_rlabel)
}

/// Sets the atom's MDL RLabel (this is an integer from 0 to 99).
/// Setting to 0 clears the rlabel.
#[pyfunction]
#[pyo3(name = "SetAtomRLabel", signature = (atom, rlabel))]
fn py_set_atom_rlabel(py: Python<'_>, mut atom: PyRefMut<'_, PyAtom>, rlabel: i32) -> PyResult<()> {
    atom.with_mut(py, |a| set_atom_rlabel(a, rlabel))?.py()
}

/// Returns the atom's MDL alias text
#[pyfunction]
#[pyo3(name = "GetAtomAlias", signature = (atom))]
fn py_get_atom_alias(py: Python<'_>, atom: PyRef<'_, PyAtom>) -> PyResult<String> {
    atom.with(py, get_atom_alias)
}

/// Sets the atom's MDL alias text.
/// Setting to an empty string clears the alias.
#[pyfunction]
#[pyo3(name = "SetAtomAlias", signature = (atom, rlabel))]
fn py_set_atom_alias(py: Python<'_>, mut atom: PyRefMut<'_, PyAtom>, rlabel: &str) -> PyResult<()> {
    atom.with_mut(py, |a| set_atom_alias(a, rlabel))
}

/// Returns the atom's MDL alias text
#[pyfunction]
#[pyo3(name = "GetAtomValue", signature = (atom))]
fn py_get_atom_value(py: Python<'_>, atom: PyRef<'_, PyAtom>) -> PyResult<String> {
    atom.with(py, get_atom_value)
}

/// Sets the atom's MDL alias text.
/// Setting to an empty string clears the alias.
#[pyfunction]
#[pyo3(name = "SetAtomValue", signature = (atom, rlabel))]
fn py_set_atom_value(py: Python<'_>, mut atom: PyRefMut<'_, PyAtom>, rlabel: &str) -> PyResult<()> {
    atom.with_mut(py, |a| set_atom_value(a, rlabel))
}

/// Gets the supplemental smiles label on an atom, returns an empty string if not present.
#[pyfunction]
#[pyo3(name = "GetSupplementalSmilesLabel", signature = (atom))]
fn py_get_supplemental_smiles_label(py: Python<'_>, atom: PyRef<'_, PyAtom>) -> PyResult<String> {
    atom.with(py, get_supplemental_smiles_label)
}

/// Sets a supplemental label on an atom that is written to the smiles string.
///
/// >>> m = Chem.MolFromSmiles("C")
/// >>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')
/// >>> Chem.MolToSmiles(m)
/// 'C<xxx>'
#[pyfunction]
#[pyo3(name = "SetSupplementalSmilesLabel", signature = (atom, label))]
fn py_set_supplemental_smiles_label(
    py: Python<'_>,
    mut atom: PyRefMut<'_, PyAtom>,
    label: &str,
) -> PyResult<()> {
    atom.with_mut(py, |a| set_supplemental_smiles_label(a, label))
}

/// Returns the number of electrons an atom is using for pi bonding
#[pyfunction]
#[pyo3(name = "GetNumPiElectrons", signature = (atom))]
fn py_get_num_pi_electrons(atom: PyRef<'_, PyAtom>) -> PyResult<u32> {
    atom.with_owner(|mol, idx| num_pi_electrons(mol, idx).py())
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_get_atom_rlabel, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_atom_rlabel, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_atom_alias, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_atom_alias, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_atom_value, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_atom_value, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_supplemental_smiles_label, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_supplemental_smiles_label, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_num_pi_electrons, m)?)?;
    Ok(())
}
