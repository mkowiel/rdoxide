//! Structure-level functions of `rdoxide.Chem.rdmolops` (upstream
//! `GraphMol/Wrap/MolOps.cpp`): `GetMolFrags`, `CombineMols`, `RenumberAtoms`,
//! `GetFormalCharge`, `GetShortestPath` and `RemoveStereochemistry`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::geometry::Point3D;
use rdoxide::graphmol::chem_transforms::combine_mols;
use rdoxide::graphmol::matrices::get_shortest_path;
use rdoxide::graphmol::mol_ops::frags::get_mol_frags_with_query;
use rdoxide::graphmol::mol_ops::haptic::{dative_bonds_to_haptic, haptic_bonds_to_dative};
use rdoxide::graphmol::mol_ops::{
    get_formal_charge, get_mol_frags_as_mols, get_mol_frags_atoms, get_mol_frags_mapping,
    remove_stereochemistry, renumber_atoms,
};

use crate::errors::{
    Check, IntoPyResult, Site, argument_error, boost_type_name, check_failed, expect_list,
};
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;
use crate::pyconv::python_object_to_vect;
use crate::rdmoldescriptors::bool_arg;

/// Finds the disconnected fragments from a molecule.
///
///   For example, for the molecule 'CC(=O)[O-].[NH3+]C' GetMolFrags() returns
///   ((0, 1, 2, 3), (4, 5))
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - asMols: (optional) if this is provided and true, the fragments
///       will be returned as molecules instead of atom ids.
///     - sanitizeFrags: (optional) if this is provided and true, the fragments
///       molecules will be sanitized before returning them.
///     - frags: (optional, defaults to None) if asMols is true and this is provided
///        as an empty list, the result will be mol.GetNumAtoms() long on return and
///        will contain the fragment assignment for each Atom
///     - fragsMolAtomMapping: (optional, defaults to None) if asMols is true and this
///       is provided as an empty list, the result will be numFrags long on
///       return, and each entry will contain the indices of the Atoms in that fragment:
///       [(0, 1, 2, 3), (4, 5)]
///
///   RETURNS: a tuple of tuples with IDs for the atoms in each fragment
///            or a tuple of molecules.
#[pyfunction]
#[pyo3(name = "GetMolFrags", signature = (mol, asMols = None, sanitizeFrags = None, frags = None, fragsMolAtomMapping = None))]
#[allow(non_snake_case)]
fn py_get_mol_frags<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    asMols: Option<&Bound<'py, PyAny>>,
    sanitizeFrags: Option<&Bound<'py, PyAny>>,
    frags: Option<&Bound<'py, PyAny>>,
    fragsMolAtomMapping: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let as_mols = bool_arg(asMols, false)?;
    let sanitize = bool_arg(sanitizeFrags, true)?;
    let core = &mol.mol;
    if !as_mols {
        let atoms = get_mol_frags_atoms(core);
        return PyTuple::new(
            py,
            atoms
                .into_iter()
                .map(|f| PyTuple::new(py, f))
                .collect::<PyResult<Vec<_>>>()?,
        );
    }
    let frags = frags.filter(|f| !f.is_none());
    let mapping_out = fragsMolAtomMapping.filter(|f| !f.is_none());
    let (mols, mapping) = py
        .detach(|| get_mol_frags_as_mols(core, sanitize, true))
        .py()?;
    if let Some(out) = frags {
        let list = expect_list(out)?;
        let (_, assignment) = get_mol_frags_mapping(core);
        for f in assignment {
            list.append(f)?;
        }
    }
    if let Some(out) = mapping_out {
        let list = expect_list(out)?;
        for atoms in mapping {
            list.append(PyTuple::new(py, atoms)?)?;
        }
    }
    PyTuple::new(
        py,
        mols.into_iter()
            .map(|m| Py::new(py, PyMol { mol: m }))
            .collect::<PyResult<Vec<_>>>()?,
    )
}

/// Combine the atoms from two molecules to produce a third
#[pyfunction]
#[pyo3(name = "CombineMols", signature = (mol1, mol2, offset = None))]
fn py_combine_mols(
    py: Python<'_>,
    mol1: PyRef<'_, PyMol>,
    mol2: PyRef<'_, PyMol>,
    offset: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    let offset = match offset {
        None => Point3D::new(0.0, 0.0, 0.0),
        Some(o) => {
            o.extract::<PyRef<'_, PyPoint3D>>()
                .map_err(|_| {
                    argument_error(
                        py,
                        "CombineMols",
                        &["Mol".into(), "Mol".into(), boost_type_name(o)],
                        "CombineMols(RDKit::ROMol mol1, RDKit::ROMol mol2, RDGeom::Point3D offset=<rdoxide.Geometry.rdGeometry.Point3D object at 0x0>)",
                    )
                })?
                .p
        }
    };
    let (m1, m2) = (&mol1.mol, &mol2.mol);
    let res = py.detach(|| combine_mols(m1, m2, offset)).py()?;
    Ok(PyMol { mol: res })
}

/// One way of showing haptic bonds (such as cyclopentadiene to
/// iron in ferrocene) is to use a dummy atom with a dative bond to the
/// iron atom with the bond labelled with the atoms involved in the
/// organic end of the bond.  Another way is to have explicit dative
/// bonds from the atoms of the haptic group to the metal atom.  This
/// function converts the former representation to the latter.
///
/// ARGUMENTS:
///
///   - mol: the molecule to use
///
/// RETURNS:
///   a modified copy of the molecule
#[pyfunction]
#[pyo3(name = "HapticBondsToDative", signature = (mol))]
fn py_haptic_bonds_to_dative(py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
    let mut res = mol.mol.clone();
    py.detach(|| haptic_bonds_to_dative(&mut res)).py()?;
    Ok(PyMol { mol: res })
}

/// Does the reverse of hapticBondsToDative.  If there are multiple
/// contiguous atoms attached by dative bonds to an atom (probably a metal
/// atom), the dative bonds will be replaced by a dummy atom in their
/// centre attached to the (metal) atom by a dative bond, which is
/// labelled with ENDPTS of the atoms that had the original dative bonds.
///
/// ARGUMENTS:
///
///   - mol: the molecule to use
///
/// RETURNS:
///   a modified copy of the molecule
#[pyfunction]
#[pyo3(name = "DativeBondsToHaptic", signature = (mol))]
fn py_dative_bonds_to_haptic(py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
    let mut res = mol.mol.clone();
    py.detach(|| dative_bonds_to_haptic(&mut res)).py()?;
    Ok(PyMol { mol: res })
}

/// The PDB residue name (`pdb_residue_name`) or chain id of an atom, or "".
fn pdb_field(atom: &rdoxide::graphmol::atom::Atom, residue: bool) -> String {
    match atom.monomer_info() {
        Some(info) if info.is_pdb_residue() => {
            if residue {
                info.residue_name.clone()
            } else {
                info.chain_id.clone()
            }
        }
        _ => String::new(),
    }
}

fn split_mol_by_pdb<'py>(
    py: Python<'py>,
    mol: &PyMol,
    white_list: Option<Vec<String>>,
    negate_list: bool,
    residue: bool,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let frags = get_mol_frags_with_query(
        &mol.mol,
        |_, a| pdb_field(a, residue),
        false,
        white_list.as_deref(),
        negate_list,
    )
    .py()?;
    let res = pyo3::types::PyDict::new(py);
    for (key, mol) in frags {
        res.set_item(key, Py::new(py, PyMol { mol })?)?;
    }
    Ok(res)
}

/// Splits a molecule into pieces based on PDB residue information.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - whiteList: only residues in this list will be returned
///     - negateList: if set, negates the white list inclusion logic
///
///   RETURNS: a dictionary keyed by residue name with molecules as the values
#[pyfunction]
#[pyo3(name = "SplitMolByPDBResidues", signature = (mol, whiteList = None, negateList = false))]
#[allow(non_snake_case)]
fn py_split_mol_by_pdb_residues<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    whiteList: Option<Vec<String>>,
    #[pyo3(from_py_with = crate::pyconv::truthy)] negateList: bool,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    split_mol_by_pdb(py, &mol, whiteList, negateList, true)
}

/// Splits a molecule into pieces based on PDB chain information.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - whiteList: only residues in this list will be returned
///     - negateList: if set, negates the white list inclusion logic
///
///   RETURNS: a dictionary keyed by chain id with molecules as the values
#[pyfunction]
#[pyo3(name = "SplitMolByPDBChainId", signature = (mol, whiteList = None, negateList = false))]
#[allow(non_snake_case)]
fn py_split_mol_by_pdb_chain_id<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    whiteList: Option<Vec<String>>,
    #[pyo3(from_py_with = crate::pyconv::truthy)] negateList: bool,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    split_mol_by_pdb(py, &mol, whiteList, negateList, false)
}

/// Returns a copy of a molecule with renumbered atoms
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///
///     - newOrder: the new ordering the atoms (should be numAtoms long)
///       for example: if newOrder is [3,2,0,1], then atom 3 in the original
///       molecule will be atom 0 in the new one
#[pyfunction]
#[pyo3(name = "RenumberAtoms", signature = (mol, newOrder))]
#[allow(non_snake_case)]
fn py_renumber_atoms(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    newOrder: &Bound<'_, PyAny>,
) -> PyResult<PyMol> {
    let core = &mol.mol;
    let n = core.num_atoms();
    if newOrder.len()? < n as usize {
        return Err(PyValueError::new_err(
            "atomCounts shorter than the number of atoms",
        ));
    }
    let order = python_object_to_vect(Some(newOrder), Some(n))?
        .ok_or_else(|| PyValueError::new_err("newOrder argument must be non-empty"))?;
    if order.len() != n as usize {
        return Err(check_failed(
            Site {
                check: Check::Precondition,
                file: "Code/GraphMol/Renumber.cpp",
                line: 22,
                expr: "newOrder.size() == nAts",
            },
            "bad newOrder size",
        ));
    }
    let res = py.detach(|| renumber_atoms(core, &order)).py()?;
    Ok(PyMol { mol: res })
}

/// Returns the formal charge for the molecule.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
#[pyfunction]
#[pyo3(name = "GetFormalCharge", signature = (mol))]
fn py_get_formal_charge(mol: PyRef<'_, PyMol>) -> i32 {
    get_formal_charge(&mol.mol)
}

/// Find the shortest path between two atoms using the Bellman-Ford algorithm.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - idx1: index of the first atom
///     - idx2: index of the second atom
#[pyfunction]
#[pyo3(name = "GetShortestPath", signature = (mol, aid1, aid2))]
fn py_get_shortest_path<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    aid1: i32,
    aid2: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    let n = mol.mol.num_atoms() as i32;
    if aid1 < 0 || aid1 >= n || aid2 < 0 || aid2 >= n {
        return Err(PyValueError::new_err("bad atom index"));
    }
    if aid1 == aid2 {
        return Err(check_failed(
            Site {
                check: Check::Invariant,
                file: "Code/GraphMol/Matrices.cpp",
                line: 339,
                expr: "aid1 != aid2",
            },
            "",
        ));
    }
    let path = get_shortest_path(&mol.mol, aid1 as u32, aid2 as u32).py()?;
    PyTuple::new(py, path)
}

/// Removes all stereochemistry info from the molecule.
#[pyfunction]
#[pyo3(name = "RemoveStereochemistry", signature = (mol))]
fn py_remove_stereochemistry(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    remove_stereochemistry(&mut mol.mol).py()
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_get_mol_frags, m)?)?;
    m.add_function(wrap_pyfunction!(py_combine_mols, m)?)?;
    m.add_function(wrap_pyfunction!(py_renumber_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(py_split_mol_by_pdb_residues, m)?)?;
    m.add_function(wrap_pyfunction!(py_split_mol_by_pdb_chain_id, m)?)?;
    m.add_function(wrap_pyfunction!(py_haptic_bonds_to_dative, m)?)?;
    m.add_function(wrap_pyfunction!(py_dative_bonds_to_haptic, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_formal_charge, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_shortest_path, m)?)?;
    m.add_function(wrap_pyfunction!(py_remove_stereochemistry, m)?)?;
    Ok(())
}
