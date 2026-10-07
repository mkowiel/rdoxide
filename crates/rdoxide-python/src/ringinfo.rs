//! `rdchem.RingInfo` — port of upstream `Code/GraphMol/Wrap/RingInfo.cpp`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::{FindRingType, RingInfo};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// contains information about a molecule's rings
#[pyclass(name = "RingInfo", module = "rdkit.Chem.rdchem")]
#[derive(Debug)]
pub struct PyRingInfo {
    pub mol: Py<PyMol>,
}

impl PyRingInfo {
    fn with<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&RingInfo) -> rdoxide::Result<R>,
    ) -> PyResult<R> {
        f(self.mol.bind(py).borrow().mol.ring_info()).py()
    }
    fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut RingInfo) -> rdoxide::Result<R>,
    ) -> PyResult<R> {
        f(self.mol.bind(py).borrow_mut().mol.ring_info_mut()).py()
    }
}

fn tuples<'py>(py: Python<'py>, rings: &[Vec<i32>]) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(
        py,
        rings.iter().map(|r| PyTuple::new(py, r).expect("tuple")),
    )
}

#[pymethods]
impl PyRingInfo {
    #[pyo3(name = "IsAtomInRingOfSize")]
    fn is_atom_in_ring_of_size(&self, py: Python<'_>, idx: u32, size: u32) -> PyResult<bool> {
        self.with(py, |r| r.is_atom_in_ring_of_size(idx, size))
    }
    #[pyo3(name = "MinAtomRingSize")]
    fn min_atom_ring_size(&self, py: Python<'_>, idx: u32) -> PyResult<u32> {
        self.with(py, |r| r.min_atom_ring_size(idx))
    }
    #[pyo3(name = "AreAtomsInSameRing")]
    fn are_atoms_in_same_ring(&self, py: Python<'_>, idx1: u32, idx2: u32) -> PyResult<bool> {
        self.with(py, |r| r.are_atoms_in_same_ring(idx1, idx2))
    }
    #[pyo3(name = "AreAtomsInSameRingOfSize")]
    fn are_atoms_in_same_ring_of_size(
        &self,
        py: Python<'_>,
        idx1: u32,
        idx2: u32,
        size: u32,
    ) -> PyResult<bool> {
        self.with(py, |r| r.are_atoms_in_same_ring_of_size(idx1, idx2, size))
    }
    #[pyo3(name = "IsBondInRingOfSize")]
    fn is_bond_in_ring_of_size(&self, py: Python<'_>, idx: u32, size: u32) -> PyResult<bool> {
        self.with(py, |r| r.is_bond_in_ring_of_size(idx, size))
    }
    #[pyo3(name = "MinBondRingSize")]
    fn min_bond_ring_size(&self, py: Python<'_>, idx: u32) -> PyResult<u32> {
        self.with(py, |r| r.min_bond_ring_size(idx))
    }
    #[pyo3(name = "AreBondsInSameRing")]
    fn are_bonds_in_same_ring(&self, py: Python<'_>, idx1: u32, idx2: u32) -> PyResult<bool> {
        self.with(py, |r| r.are_bonds_in_same_ring(idx1, idx2))
    }
    #[pyo3(name = "AreBondsInSameRingOfSize")]
    fn are_bonds_in_same_ring_of_size(
        &self,
        py: Python<'_>,
        idx1: u32,
        idx2: u32,
        size: u32,
    ) -> PyResult<bool> {
        self.with(py, |r| r.are_bonds_in_same_ring_of_size(idx1, idx2, size))
    }
    #[pyo3(name = "NumAtomRings")]
    fn num_atom_rings(&self, py: Python<'_>, idx: u32) -> PyResult<u32> {
        self.with(py, |r| r.num_atom_rings(idx))
    }
    #[pyo3(name = "NumBondRings")]
    fn num_bond_rings(&self, py: Python<'_>, idx: u32) -> PyResult<u32> {
        self.with(py, |r| r.num_bond_rings(idx))
    }
    #[pyo3(name = "NumRings")]
    fn num_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, RingInfo::num_rings)
    }
    #[pyo3(name = "IsRingFused")]
    #[allow(non_snake_case)]
    fn is_ring_fused(&self, py: Python<'_>, ringIdx: u32) -> PyResult<bool> {
        self.with_mut(py, |r| r.is_ring_fused(ringIdx))
    }
    #[pyo3(name = "AreRingsFused")]
    #[allow(non_snake_case)]
    fn are_rings_fused(&self, py: Python<'_>, ring1Idx: u32, ring2Idx: u32) -> PyResult<bool> {
        self.with_mut(py, |r| r.are_rings_fused(ring1Idx, ring2Idx))
    }
    #[pyo3(name = "NumFusedBonds")]
    #[allow(non_snake_case)]
    fn num_fused_bonds(&self, py: Python<'_>, ringIdx: u32) -> PyResult<u32> {
        self.with_mut(py, |r| r.num_fused_bonds(ringIdx))
    }
    #[pyo3(name = "AtomRings")]
    fn atom_rings<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let rings = self.with(py, |r| Ok(r.atom_rings().to_vec()))?;
        tuples(py, &rings)
    }
    #[pyo3(name = "BondRings")]
    fn bond_rings<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let rings = self.with(py, |r| Ok(r.bond_rings().to_vec()))?;
        tuples(py, &rings)
    }
    #[pyo3(name = "AtomMembers")]
    fn atom_members<'py>(&self, py: Python<'py>, idx: u32) -> PyResult<Bound<'py, PyTuple>> {
        let m = self.with(py, |r| r.atom_members(idx).map(<[i32]>::to_vec))?;
        PyTuple::new(py, m)
    }
    #[pyo3(name = "BondMembers")]
    fn bond_members<'py>(&self, py: Python<'py>, idx: u32) -> PyResult<Bound<'py, PyTuple>> {
        let m = self.with(py, |r| r.bond_members(idx).map(<[i32]>::to_vec))?;
        PyTuple::new(py, m)
    }
    #[pyo3(name = "AtomRingSizes")]
    fn atom_ring_sizes<'py>(&self, py: Python<'py>, idx: u32) -> PyResult<Bound<'py, PyTuple>> {
        let m = self.with(py, |r| r.atom_ring_sizes(idx))?;
        PyTuple::new(py, m)
    }
    #[pyo3(name = "BondRingSizes")]
    fn bond_ring_sizes<'py>(&self, py: Python<'py>, idx: u32) -> PyResult<Bound<'py, PyTuple>> {
        let m = self.with(py, |r| r.bond_ring_sizes(idx))?;
        PyTuple::new(py, m)
    }
    #[pyo3(name = "NumRingFamilies")]
    fn num_ring_families(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, RingInfo::num_ring_families)
    }
    #[pyo3(name = "AtomRingFamilies")]
    fn atom_ring_families<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let rings = self.with(py, |r| Ok(r.atom_ring_families().to_vec()))?;
        tuples(py, &rings)
    }
    #[pyo3(name = "BondRingFamilies")]
    fn bond_ring_families<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let rings = self.with(py, |r| Ok(r.bond_ring_families().to_vec()))?;
        tuples(py, &rings)
    }
    #[pyo3(name = "AreRingFamiliesInitialized")]
    fn are_ring_families_initialized(&self) -> bool {
        false
    }
    /// Adds a ring to the set. Be very careful with this operation.
    #[pyo3(name = "AddRing")]
    #[allow(non_snake_case)]
    fn add_ring(&self, py: Python<'_>, atomIds: Vec<i32>, bondIds: Vec<i32>) -> PyResult<()> {
        if atomIds.len() != bondIds.len() {
            return Err(PyValueError::new_err("list sizes must match"));
        }
        self.with_mut(py, |r| {
            if !r.is_initialized() {
                r.initialize(FindRingType::OtherOrUnknown);
            }
            r.add_ring(atomIds, bondIds).map(|_| ())
        })
    }
}
