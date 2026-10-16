//! Attachment-point helpers of `rdoxide.Chem.rdmolops` (upstream
//! `GraphMol/Wrap/MolOps.cpp`): `ExpandAttachmentPoints`,
//! `CollapseAttachmentPoints`, `GetAttachmentPointLabelNumber`,
//! `IsMarkedAttachmentPoint` and `ATTACHMENT_POINT_LABEL_PREFIX`.

use pyo3::prelude::*;
use rdoxide::graphmol::attachment_points::{
    ATTACHMENT_POINT_LABEL_PREFIX, collapse_attachment_points, expand_attachment_points,
    is_marked_attachment_point,
};
use rdoxide::graphmol::smiles_parse::ops::get_attachment_point_label_number;

use crate::atom::PyAtom;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// attachment points encoded as attachPt properties are added to the graph as dummy atoms
///
///   Arguments:
///    - mol: molecule to be modified
///    - addAsQueries: if true, the dummy atoms will be added as null queries
///         (i.e. they will match any atom in a substructure search)
///    - addCoords: if true and the molecule has one or more conformers,
///         positions for the attachment points will be added to the conformer(s)
#[pyfunction]
#[pyo3(name = "ExpandAttachmentPoints", signature = (mol, addAsQueries = true, addCoords = true))]
#[allow(non_snake_case)]
fn py_expand_attachment_points(
    mut mol: PyRefMut<'_, PyMol>,
    addAsQueries: bool,
    addCoords: bool,
) -> PyResult<()> {
    expand_attachment_points(&mut mol.mol, addAsQueries, addCoords).py()
}

/// dummy atoms in the graph are removed and replaced with attachment point annotations on the attached atoms
///
///   Arguments:
///    - mol: molecule to be modified
///    - markedOnly: if true, only dummy atoms with the _fromAttachPoint property
///      or a valid _AP<n> atom label will be collapsed. The label suffix is not
///      interpreted as an MDL attachment-point position.
///
///   In order for a dummy atom to be considered for collapsing it must have:
///    - degree 1 with a single or unspecified bond
///    - the bond to it can not be wedged
///    - either no query or be an AtomNullQuery
#[pyfunction]
#[pyo3(name = "CollapseAttachmentPoints", signature = (mol, markedOnly = true))]
#[allow(non_snake_case)]
fn py_collapse_attachment_points(mut mol: PyRefMut<'_, PyMol>, markedOnly: bool) -> PyResult<()> {
    collapse_attachment_points(&mut mol.mol, markedOnly).py()
}

/// returns the positive integer from a valid _AP<n> attachment-point label
///
///   The atom must be a degree-one dummy atom. Returns 0 if it does not have a
///   valid numbered attachment-point label. The returned number is a label
///   identifier, not an MDL ATTCHPT position.
///
///   Arguments:
///    - atom: the atom to inspect
#[pyfunction]
#[pyo3(name = "GetAttachmentPointLabelNumber", signature = (atom))]
fn py_get_attachment_point_label_number(atom: PyRef<'_, PyAtom>) -> PyResult<u32> {
    atom.with_owner(|mol, idx| Ok(get_attachment_point_label_number(mol, idx)))
}

/// returns whether an atom is a marked explicit attachment point
///
///   A marked attachment point is a degree-one dummy atom with the
///   _fromAttachPoint property or a valid _AP<n> atom label, where n is a
///   positive decimal integer. This checks attachment-point identity only, not
///   whether the atom can currently be collapsed. In particular, an attachment
///   point connected by a wedged bond is still considered marked.
///
///   Arguments:
///    - atom: the atom to inspect
#[pyfunction]
#[pyo3(name = "IsMarkedAttachmentPoint", signature = (atom))]
fn py_is_marked_attachment_point(atom: PyRef<'_, PyAtom>) -> PyResult<bool> {
    atom.with_owner(|mol, idx| Ok(is_marked_attachment_point(mol, idx)))
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_expand_attachment_points, m)?)?;
    m.add_function(wrap_pyfunction!(py_collapse_attachment_points, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_attachment_point_label_number, m)?)?;
    m.add_function(wrap_pyfunction!(py_is_marked_attachment_point, m)?)?;
    m.add(
        "ATTACHMENT_POINT_LABEL_PREFIX",
        ATTACHMENT_POINT_LABEL_PREFIX,
    )?;
    Ok(())
}
