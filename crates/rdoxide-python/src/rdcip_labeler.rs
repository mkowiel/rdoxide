//! `rdoxide.Chem.rdCIPLabeler` (upstream `GraphMol/CIPLabeler/Wrap/rdCIPLabeler.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::cip_labeler::assign_cip_labels_selected;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::pyconv::python_object_to_vect;

/// New implementation of Stereo assignment using a true CIP ranking.
/// On return:  The molecule to contains CIP flags
/// Errors:  when maxRecursiveIterations is exceeded, throws a MaxIterationsExceeded error
/// ARGUMENTS:
///
///  - mol: the molecule
///  - atomsToLabel: (optional) list of atoms to label
///  - bondsToLabel: (optional) list of bonds to label
///  - maxRecursiveIterations: (optional) protects against pseudo-infinite
/// recursion for highly symmetrical structures. The limit is shared by the preliminary and full labeling passes.
///  A value of 1,250,000 takes about 1 second. Most structures require fewer than 10,000 recursive comparisons.
///  A peptide with MW~3000 took about 100 comparisons, and a 20,000 MW protein took about 600 comparisons
/// (0 = default - no limit)
#[pyfunction]
#[pyo3(name = "AssignCIPLabels", signature = (mol, atomsToLabel = None, bondsToLabel = None, maxRecursiveIterations = 0))]
#[allow(non_snake_case)]
fn py_assign_cip_labels(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    atomsToLabel: Option<&Bound<'_, PyAny>>,
    bondsToLabel: Option<&Bound<'_, PyAny>>,
    maxRecursiveIterations: u32,
) -> PyResult<()> {
    let num_atoms = mol.mol.num_atoms();
    let num_bonds = mol.mol.num_bonds();
    let bitset = |obj: Option<&Bound<'_, PyAny>>, size: u32| -> PyResult<Vec<bool>> {
        let mut bits = vec![false; size as usize];
        for idx in python_object_to_vect(obj, Some(size))?.unwrap_or_default() {
            bits[idx as usize] = true;
        }
        Ok(bits)
    };
    let mut atoms = bitset(atomsToLabel, num_atoms)?;
    let mut bonds = bitset(bondsToLabel, num_bonds)?;
    if !atomsToLabel.map_or(Ok(false), |o| o.is_truthy())?
        && !bondsToLabel.map_or(Ok(false), |o| o.is_truthy())?
    {
        atoms.fill(true);
        bonds.fill(true);
    }
    let mol = &mut mol.mol;
    py.detach(|| assign_cip_labels_selected(mol, &atoms, &bonds, maxRecursiveIterations))
        .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing a function to assign stereochemical labels based on an accurate CIP rules implementation. This algoritm is a port of https://github.com/SiMolecule/centres, which was originally written by John Mayfield. The original algorithm is described in:\n\nHanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, A., Redkin, D.\nAlgorithmic Analysis of Cahn--Ingold--Prelog Rules of Stereochemistry:\nProposals for Revised Rules and a Guide for Machine Implementation.\nJ. Chem. Inf. Model. 2018, 58, 1755-1765.\n",
    )?;
    m.add_function(wrap_pyfunction!(py_assign_cip_labels, m)?)?;
    Ok(())
}
