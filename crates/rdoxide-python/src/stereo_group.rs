//! `rdchem.StereoGroup` — port of upstream `Code/GraphMol/Wrap/StereoGroup.cpp`.

use pyo3::prelude::*;
use rdoxide::graphmol::stereo_group::{StereoGroup, StereoGroupType};

use crate::atom::PyAtom;
use crate::bond::PyBond;
use crate::enums_impl::PyStereoGroupType;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// A collection of atoms with a defined stereochemical relationship.
///
/// Used to help represent a sample with unknown stereochemistry, or that is a mix
/// of diastereomers.
#[pyclass(name = "StereoGroup", module = "rdoxide.Chem.rdchem")]
#[derive(Debug)]
pub struct PyStereoGroup {
    pub group: StereoGroup,
    pub mol: Option<Py<PyMol>>,
}

#[pymethods]
impl PyStereoGroup {
    /// Returns the StereoGroupType.
    #[pyo3(name = "GetGroupType")]
    fn get_group_type(&self) -> PyStereoGroupType {
        PyStereoGroupType(self.group.group_type())
    }
    /// access the atoms in the StereoGroup.
    #[pyo3(name = "GetAtoms")]
    fn get_atoms(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAtom>>> {
        let Some(mol) = &self.mol else {
            return Ok(Vec::new());
        };
        self.group
            .atoms()
            .iter()
            .map(|&a| PyAtom::in_mol(py, mol.clone_ref(py), a))
            .collect()
    }
    /// access the bonds in the StereoGroup.
    #[pyo3(name = "GetBonds")]
    fn get_bonds(&self, py: Python<'_>) -> PyResult<Vec<Py<PyBond>>> {
        let Some(mol) = &self.mol else {
            return Ok(Vec::new());
        };
        self.group
            .bonds()
            .iter()
            .map(|&b| PyBond::in_mol(py, mol.clone_ref(py), b))
            .collect()
    }
    /// return the StereoGroup's original ID.
    /// Note that the ID only makes sense for AND/OR groups.
    #[pyo3(name = "GetReadId")]
    fn get_read_id(&self) -> u32 {
        self.group.read_id()
    }
    /// return the StereoGroup's ID that will be exported.
    /// Note that the ID only makes sense for AND/OR groups.
    #[pyo3(name = "GetWriteId")]
    fn get_write_id(&self) -> u32 {
        self.group.write_id()
    }
    /// return the StereoGroup's ID that will be exported.
    #[pyo3(name = "SetWriteId")]
    fn set_write_id(&mut self, id: u32) {
        self.group.set_write_id(id);
    }
}

/// creates a StereoGroup associated with a molecule from a list of atom Ids
#[pyfunction]
#[pyo3(name = "CreateStereoGroup", signature = (stereoGroupType, mol, atomIds = Vec::new(), bondIds = Vec::new(), readId = 0))]
#[allow(non_snake_case)]
pub fn create_stereo_group(
    stereoGroupType: &Bound<'_, PyAny>,
    mol: &Bound<'_, PyMol>,
    atomIds: Vec<u32>,
    bondIds: Vec<u32>,
    readId: u32,
) -> PyResult<PyStereoGroup> {
    let t: StereoGroupType = PyStereoGroupType::extract(stereoGroupType)?;
    {
        let m = mol.borrow();
        for &a in &atomIds {
            m.mol.get_atom(a).py()?;
        }
        for &b in &bondIds {
            m.mol.get_bond(b).py()?;
        }
    }
    let group = StereoGroup::new(t, atomIds, bondIds, readId).py()?;
    Ok(PyStereoGroup {
        group,
        mol: Some(mol.clone().unbind()),
    })
}
