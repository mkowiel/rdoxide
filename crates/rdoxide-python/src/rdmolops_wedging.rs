//! Bond wedging of `rdkit.Chem.rdmolops` (upstream `GraphMol/Wrap/MolOps.cpp`):
//! `BondWedgingParameters`, `WedgeMolBonds`, `WedgeBond`,
//! `ReapplyMolBlockWedging` and `RemoveNonExplicit3DChirality`.

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::Conformer;
use rdoxide::graphmol::chirality::remove_non_explicit_3d_chirality;
use rdoxide::graphmol::wedge_bonds::{
    BondWedgingParameters, reapply_mol_block_wedging, wedge_bond, wedge_mol_bonds,
};

use crate::bond::PyBond;
use crate::conformer::{ConfHandle, PyConformer};
use crate::errors::{IntoPyResult, to_py};
use crate::mol::PyMol;

/// Parameters controlling how bond wedging is done.
#[pyclass(
    name = "BondWedgingParameters",
    module = "rdkit.Chem.rdmolops",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyBondWedgingParameters {
    /// If this is enabled then two bonds will be wedged at chiral
    /// centers subject to the following constraints:
    ///   1. ring bonds will not be wedged
    ///   2. bonds to chiral centers will not be wedged
    ///   3. bonds separated by more than 120 degrees will not be
    ///       wedged
    #[pyo3(get, set, name = "wedgeTwoBondsIfPossible")]
    wedge_two_bonds_if_possible: bool,
}

#[pymethods]
impl PyBondWedgingParameters {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }
}

impl From<Option<PyRef<'_, PyBondWedgingParameters>>> for PyBondWedgingParameters {
    fn from(params: Option<PyRef<'_, PyBondWedgingParameters>>) -> Self {
        params.map(|p| p.clone()).unwrap_or_default()
    }
}

/// Set the wedging on single bonds in a molecule.
///    The wedging scheme used is that from Mol files.
///
///   ARGUMENTS:
///
///     - molecule: the molecule to update
///     - conformer: the conformer to use to determine wedge direction
#[pyfunction]
#[pyo3(name = "WedgeMolBonds", signature = (mol, conformer, params = None))]
fn py_wedge_mol_bonds(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    conformer: &Bound<'_, PyConformer>,
    params: Option<PyRef<'_, PyBondWedgingParameters>>,
) -> PyResult<()> {
    // upstream checks the owner while looking for atropisomer bonds
    match &conformer.borrow().handle {
        ConfHandle::InMol { mol: owner, .. } if owner.is(mol) => {}
        ConfHandle::InMol { .. } => {
            return Err(to_py(RDKitError::precondition_at(
                "conformer does not belong to molecule",
                "conf == nullptr || &(conf->getOwningMol()) == &mol",
                "Code/GraphMol/Atropisomers.cpp",
                1231,
            )));
        }
        ConfHandle::Owned(_) => {
            return Err(to_py(RDKitError::precondition_at(
                "no owner",
                "dp_mol",
                "Code/GraphMol/Conformer.h",
                92,
            )));
        }
    }
    // copied first: the conformer belongs to `mol`, which is borrowed mutably below
    let conf = conformer.borrow().with(py, Conformer::clone)?;
    let params = PyBondWedgingParameters::from(params);
    let params = BondWedgingParameters {
        wedge_two_bonds_if_possible: params.wedge_two_bonds_if_possible,
    };
    wedge_mol_bonds(&mut mol.borrow_mut().mol, Some(&conf), &params).py()
}

/// Set the wedging on an individual bond from a molecule.
///    The wedging scheme used is that from Mol files.
///   ARGUMENTS:
///     - bond: the bond to update
///     - atom ID: the atom from which to do the wedging
///     - conformer: the conformer to use to determine wedge direction
#[pyfunction]
#[pyo3(name = "WedgeBond", signature = (bond, fromAtomIdx, conf))]
#[allow(non_snake_case)]
fn py_wedge_bond(
    py: Python<'_>,
    bond: &Bound<'_, PyBond>,
    fromAtomIdx: u32,
    conf: &Bound<'_, PyConformer>,
) -> PyResult<()> {
    let same_mol = match (&bond.borrow().handle, &conf.borrow().handle) {
        (
            crate::bond::BondHandle::InMol { mol: bond_mol, .. },
            ConfHandle::InMol { mol: conf_mol, .. },
        ) => bond_mol.is(conf_mol),
        _ => false,
    };
    if !same_mol {
        return Err(to_py(RDKitError::precondition_at(
            "bond and conformer do not belong to same molecule",
            "&conf->getOwningMol() == &bond->getOwningMol()",
            "Code/GraphMol/WedgeBonds.cpp",
            544,
        )));
    }
    let conf = conf.borrow().with(py, Conformer::clone)?;
    bond.borrow()
        .with_mol_mut(py, |mol, idx| wedge_bond(mol, idx, fromAtomIdx, &conf))?
        .py()
}

/// Set the wedging to that which was read from the original
///      MolBlock, over-riding anything that was originally there.
///
///           ARGUMENTS:
///
///             - molecule: the molecule to update
///             - allBondTypes: reapply the wedging also on bonds other
///               than single and aromatic ones
///             - verify: if true, the function will check that the wedges are only
///               applied in sensible places (i.e.single bonds connected to chiral
///               centers or atropisomeric bonds)
#[pyfunction]
#[pyo3(name = "ReapplyMolBlockWedging", signature = (mol, allBondTypes = true, verify = false))]
#[allow(non_snake_case)]
fn py_reapply_mol_block_wedging(
    mut mol: PyRefMut<'_, PyMol>,
    allBondTypes: bool,
    verify: bool,
) -> PyResult<()> {
    reapply_mol_block_wedging(&mut mol.mol, allBondTypes, verify).py()
}

/// Remove chiral markings that were derived from a 3D mol but were not
///         explicity marked in the mol block. (wedge bond or CFG indication
///
///           ARGUMENTS:
///
///             - molecule: the molecule to update
#[pyfunction]
#[pyo3(name = "RemoveNonExplicit3DChirality", signature = (mol))]
fn py_remove_non_explicit_3d_chirality(mut mol: PyRefMut<'_, PyMol>) {
    remove_non_explicit_3d_chirality(&mut mol.mol);
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBondWedgingParameters>()?;
    m.add_function(wrap_pyfunction!(py_wedge_mol_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(py_wedge_bond, m)?)?;
    m.add_function(wrap_pyfunction!(py_reapply_mol_block_wedging, m)?)?;
    m.add_function(wrap_pyfunction!(py_remove_non_explicit_3d_chirality, m)?)?;
    Ok(())
}
