//! Substructure transforms and stereo-bond marking of `rdkit.Chem.rdmolops`
//! (upstream `GraphMol/Wrap/MolOps.cpp`): `DeleteSubstructs` and
//! `FindPotentialStereoBonds`.

use pyo3::prelude::*;
use rdoxide::graphmol::chem_transforms::delete_substructs;
use rdoxide::graphmol::chirality::legacy::find_potential_stereo_bonds;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// Removes atoms matching a substructure query from a molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///
///     - query: the molecule to be used as a substructure query
///
///     - onlyFrags: (optional) if this toggle is set, atoms will only be removed if
///       the entire fragment in which they are found is matched by the query.
///       See below for examples.
///       Default value is 0 (remove the atoms whether or not the entire fragment matches)
///
///     - useChirality: (optional) match the substructure query using chirality
///
///   RETURNS: a new molecule with the substructure removed
///
///   NOTES:
///
///     - The original molecule is *not* modified.
///
///   EXAMPLES:
///
///    The following examples substitute SMILES/SMARTS strings for molecules, you'd have
///    to actually use molecules:
///
///     - DeleteSubstructs('CCOC','OC') -> 'CC'
///
///     - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'
///
///     - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'
///
///     - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'
#[pyfunction]
#[pyo3(name = "DeleteSubstructs", signature = (mol, query, onlyFrags = None, useChirality = None))]
#[allow(non_snake_case)]
fn py_delete_substructs(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    query: PyRef<'_, PyMol>,
    onlyFrags: Option<&Bound<'_, PyAny>>,
    useChirality: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    let only_frags = bool_arg(onlyFrags, false)?;
    let use_chirality = bool_arg(useChirality, false)?;
    let (mol, query) = (&mol.mol, &query.mol);
    let res = py
        .detach(|| delete_substructs(mol, query, only_frags, use_chirality))
        .py()?;
    Ok(PyMol { mol: res })
}

/// Find bonds that could be cis/trans in a molecule and mark them as Bond.STEREOANY.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - cleanIt: (optional) if this option is set to true, any previous marking of _CIPCode
///                on the bond is cleared - otherwise it is left untouched
#[pyfunction]
#[pyo3(name = "FindPotentialStereoBonds", signature = (mol, cleanIt = None))]
#[allow(non_snake_case)]
fn py_find_potential_stereo_bonds(
    mut mol: PyRefMut<'_, PyMol>,
    cleanIt: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let clean_it = bool_arg(cleanIt, false)?;
    find_potential_stereo_bonds(&mut mol.mol, clean_it).py()
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_delete_substructs, m)?)?;
    m.add_function(wrap_pyfunction!(py_find_potential_stereo_bonds, m)?)?;
    Ok(())
}
