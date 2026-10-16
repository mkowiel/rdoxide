//! Stereo-perception helpers and small atom/molecule predicates of
//! `rdoxide.Chem.rdmolops` (upstream `GraphMol/Wrap/MolOps.cpp`):
//! `AssignChiralTypesFromBondDirs`, `AssignStereochemistryFrom3D`,
//! `AssignAtomChiralTags*`, `DetectBondStereo(c|C)hemistry`,
//! `SetBondStereoFromDirections`, `SetTerminalAtomCoords`, `Cleanup*`,
//! `NeedsHs`, `CountAtomElec` and `AtomHasConjugatedBond`.

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::aromaticity::count_atom_elec;
use rdoxide::graphmol::chirality::bond_dirs::set_bond_stereo_from_directions;
use rdoxide::graphmol::chirality::from_coords::{
    assign_chiral_types_from_3d, assign_chiral_types_from_bond_dirs,
    assign_chiral_types_from_mol_parity, assign_stereochemistry_from_3d,
    detect_bond_stereochemistry,
};
use rdoxide::graphmol::conjug_hybrid::atom_has_conjugated_bond;
use rdoxide::graphmol::mol_ops::add_hs::{needs_hs, set_terminal_atom_coords};
use rdoxide::graphmol::mol_ops::{
    clean_up_organometallics, cleanup_atropisomers, cleanup_chirality,
};

use crate::atom::{AtomHandle, PyAtom};
use crate::conformer::{ConfHandle, PyConformer};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Uses bond directions to assign ChiralTypes to a molecule's atoms.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - confId: (optional) the conformation to use
///     - replaceExistingTags: (optional) replace any existing information about stereochemistry
#[pyfunction]
#[pyo3(name = "AssignChiralTypesFromBondDirs", signature = (mol, confId = -1, replaceExistingTags = true))]
#[allow(non_snake_case)]
fn py_assign_chiral_types_from_bond_dirs(
    mut mol: PyRefMut<'_, PyMol>,
    confId: i32,
    replaceExistingTags: bool,
) -> PyResult<()> {
    assign_chiral_types_from_bond_dirs(&mut mol.mol, confId, replaceExistingTags).py()
}

/// Uses a conformer (should be 3D) to assign ChiralTypes to a molecule's atoms
///         and stereo flags to its bonds
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - confId: (optional) the conformation to use
///     - replaceExistingTags: (optional) replace any existing information about stereochemistry
#[pyfunction]
#[pyo3(name = "AssignStereochemistryFrom3D", signature = (mol, confId = -1, replaceExistingTags = true))]
#[allow(non_snake_case)]
fn py_assign_stereochemistry_from_3d(
    mut mol: PyRefMut<'_, PyMol>,
    confId: i32,
    replaceExistingTags: bool,
) -> PyResult<()> {
    assign_stereochemistry_from_3d(&mut mol.mol, confId, replaceExistingTags).py()
}

/// Sets the chiral tags on a molecule's atoms based on
///   a 3D conformation.
///   NOTE that this does not check to see if atoms are chiral centers (i.e. all
///   substituents are different), it merely sets the chiral type flags based on the
///   coordinates and atom ordering. Use AssignStereochemistryFrom3D() if you
///   want chiral flags only on actual stereocenters.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - confId: the conformer id to use, -1 for the default
///     - replaceExistingTags: if True, existing stereochemistry information will be cleared
///     before running the calculation.
#[pyfunction]
#[pyo3(name = "AssignAtomChiralTagsFromStructure", signature = (mol, confId = -1, replaceExistingTags = true))]
#[allow(non_snake_case)]
fn py_assign_atom_chiral_tags_from_structure(
    mut mol: PyRefMut<'_, PyMol>,
    confId: i32,
    replaceExistingTags: bool,
) -> PyResult<()> {
    assign_chiral_types_from_3d(&mut mol.mol, confId, replaceExistingTags).py()
}

/// Sets the chiral tags on a molecule's atoms based on
///   the molParity atom property.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - replaceExistingTags: if True, existing stereochemistry information will be cleared
///     before running the calculation.
#[pyfunction]
#[pyo3(name = "AssignAtomChiralTagsFromMolParity", signature = (mol, replaceExistingTags = true))]
#[allow(non_snake_case)]
fn py_assign_atom_chiral_tags_from_mol_parity(
    mut mol: PyRefMut<'_, PyMol>,
    replaceExistingTags: bool,
) -> PyResult<()> {
    assign_chiral_types_from_mol_parity(&mut mol.mol, replaceExistingTags).py()
}

/// Assign stereochemistry to bonds based on coordinates and a conformer.
///         DEPRECATED
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///     - conformer: Conformer providing the coordinates
#[pyfunction]
#[pyo3(name = "DetectBondStereoChemistry", signature = (mol, conformer))]
fn py_detect_bond_stereo_chemistry(
    mol: &Bound<'_, PyMol>,
    conformer: &Bound<'_, PyConformer>,
) -> PyResult<()> {
    let id = match &conformer.borrow().handle {
        ConfHandle::InMol { mol: owner, id } if owner.is(mol) => *id,
        ConfHandle::Owned(_) => {
            return Err(crate::errors::to_py(RDKitError::precondition_at(
                "no owner",
                "dp_mol",
                "Code/GraphMol/Conformer.h",
                92,
            )));
        }
        ConfHandle::InMol { .. } => {
            return Err(crate::errors::to_py(RDKitError::precondition_at(
                "conformer does not belong to molecule",
                "&(conf->getOwningMol()) == &mol",
                "Code/GraphMol/FileParsers/MolFileStereochem.cpp",
                216,
            )));
        }
    };
    let mut mol = mol.borrow_mut();
    let conf_id = i32::try_from(id).unwrap_or(-1);
    detect_bond_stereochemistry(&mut mol.mol, conf_id).py()
}

/// DEPRECATED
///
///     - mol: the molecule to be modified
///     - confId: Conformer to use for the coordinates
#[pyfunction]
#[pyo3(name = "DetectBondStereochemistry", signature = (mol, confId = -1))]
#[allow(non_snake_case)]
fn py_detect_bond_stereochemistry(mut mol: PyRefMut<'_, PyMol>, confId: i32) -> PyResult<()> {
    detect_bond_stereochemistry(&mut mol.mol, confId).py()
}

/// Uses the directions of neighboring bonds to set cis/trans stereo on double bonds.
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
#[pyfunction]
#[pyo3(name = "SetBondStereoFromDirections", signature = (mol))]
fn py_set_bond_stereo_from_directions(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    set_bond_stereo_from_directions(&mut mol.mol).py()
}

/// Sets Cartesian coordinates for a terminal atom.
///
///   Useful for growing a molecule (e.g. by adding Hs) when the geometry is known.
///   NOTE: this sets the appropriate coordinates in all of the molecule's conformers
///   ARGUMENTS:
///
///     - mol: the molecule the atoms belong to.
///     - idx: index of the terminal atom whose coordinates are set.
///     - mol: index of the bonded neighbor atom.
///
///   RETURNS: Nothing
#[pyfunction]
#[pyo3(name = "SetTerminalAtomCoords", signature = (mol, idx, otherIdx))]
#[allow(non_snake_case)]
fn py_set_terminal_atom_coords(
    mut mol: PyRefMut<'_, PyMol>,
    idx: u32,
    otherIdx: u32,
) -> PyResult<()> {
    set_terminal_atom_coords(&mut mol.mol, idx, otherIdx).py()
}

/// Removes chirality information from atoms that cannot be stereocenters
#[pyfunction]
#[pyo3(name = "CleanupChirality", signature = (mol))]
fn py_cleanup_chirality(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    cleanup_chirality(&mut mol.mol).py()
}

/// Removes atropisomer stereo from bonds that cannot be atropisomers
#[pyfunction]
#[pyo3(name = "CleanupAtropisomers", signature = (mol))]
fn py_cleanup_atropisomers(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    cleanup_atropisomers(&mut mol.mol).py()
}

/// Cleans up organometallic structures (dative bonds, charges) in place
#[pyfunction]
#[pyo3(name = "CleanupOrganometallics", signature = (mol))]
fn py_cleanup_organometallics(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    clean_up_organometallics(&mut mol.mol).py()
}

/// returns whether or not the molecule needs to have Hs added
#[pyfunction]
#[pyo3(name = "NeedsHs", signature = (mol))]
fn py_needs_hs(mol: PyRef<'_, PyMol>) -> PyResult<bool> {
    needs_hs(&mol.mol).py()
}

/// returns the number of electrons available on an atom to donate for aromaticity
#[pyfunction]
#[pyo3(name = "CountAtomElec", signature = (atom))]
fn py_count_atom_elec(atom: PyRef<'_, PyAtom>) -> PyResult<i32> {
    if matches!(atom.handle, AtomHandle::Owned(_)) {
        // upstream asks the atom for its implicit Hs before it needs the owner
        return Err(crate::errors::to_py(RDKitError::precondition_at(
            "getNumImplicitHs() called without preceding call to calcImplicitValence()",
            "d_implicitValence > -1",
            "Code/GraphMol/Atom.cpp",
            307,
        )));
    }
    atom.with_owner(|mol, idx| Ok(count_atom_elec(mol, idx)))
}

/// returns whether or not the atom is involved in a conjugated bond
#[pyfunction]
#[pyo3(name = "AtomHasConjugatedBond", signature = (atom))]
fn py_atom_has_conjugated_bond(atom: PyRef<'_, PyAtom>) -> PyResult<bool> {
    atom.with_owner(|mol, idx| Ok(atom_has_conjugated_bond(mol, idx)))
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_assign_chiral_types_from_bond_dirs, m)?)?;
    m.add_function(wrap_pyfunction!(py_assign_stereochemistry_from_3d, m)?)?;
    m.add_function(wrap_pyfunction!(
        py_assign_atom_chiral_tags_from_structure,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        py_assign_atom_chiral_tags_from_mol_parity,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(py_detect_bond_stereo_chemistry, m)?)?;
    m.add_function(wrap_pyfunction!(py_detect_bond_stereochemistry, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_bond_stereo_from_directions, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_terminal_atom_coords, m)?)?;
    m.add_function(wrap_pyfunction!(py_cleanup_chirality, m)?)?;
    m.add_function(wrap_pyfunction!(py_cleanup_atropisomers, m)?)?;
    m.add_function(wrap_pyfunction!(py_cleanup_organometallics, m)?)?;
    m.add_function(wrap_pyfunction!(py_needs_hs, m)?)?;
    m.add_function(wrap_pyfunction!(py_count_atom_elec, m)?)?;
    m.add_function(wrap_pyfunction!(py_atom_has_conjugated_bond, m)?)?;
    Ok(())
}
