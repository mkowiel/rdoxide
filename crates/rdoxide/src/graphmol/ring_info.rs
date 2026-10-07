//! Ported from RDKit `Code/GraphMol/RingInfo.{h,cpp}`.
//!
//! Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::{RDKitError, Result};

/// `FIND_RING_TYPE`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum FindRingType {
    Fast,
    Sssr,
    SymmSssr,
    #[default]
    OtherOrUnknown,
}

/// Ring membership information (`RDKit::RingInfo`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RingInfo {
    initialized: bool,
    find_type: FindRingType,
    atom_members: Vec<Vec<i32>>,
    bond_members: Vec<Vec<i32>>,
    atom_rings: Vec<Vec<i32>>,
    bond_rings: Vec<Vec<i32>>,
    atom_ring_families: Vec<Vec<i32>>,
    bond_ring_families: Vec<Vec<i32>>,
    fused_rings: Vec<Vec<bool>>,
    num_fused_bonds: Vec<u32>,
}

fn not_init() -> RDKitError {
    RDKitError::precondition("RingInfo not initialized")
}

impl RingInfo {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_initialized(&self) -> bool {
        self.initialized
    }
    pub fn initialize(&mut self, ring_type: FindRingType) {
        self.initialized = true;
        self.find_type = ring_type;
    }
    pub fn ring_type(&self) -> FindRingType {
        self.find_type
    }
    pub fn is_find_fast_or_better(&self) -> bool {
        self.initialized
            && matches!(
                self.find_type,
                FindRingType::Fast | FindRingType::Sssr | FindRingType::SymmSssr
            )
    }
    pub fn is_sssr_or_better(&self) -> bool {
        self.initialized && matches!(self.find_type, FindRingType::Sssr | FindRingType::SymmSssr)
    }
    pub fn is_symm_sssr(&self) -> bool {
        self.initialized && self.find_type == FindRingType::SymmSssr
    }

    /// `reset(resetRingFamilies)`
    pub fn reset(&mut self, reset_ring_families: bool) {
        if !self.initialized {
            return;
        }
        self.initialized = false;
        self.find_type = FindRingType::OtherOrUnknown;
        self.atom_members.clear();
        self.bond_members.clear();
        self.atom_rings.clear();
        self.bond_rings.clear();
        if reset_ring_families {
            self.reset_ring_families();
        }
    }

    pub fn reset_ring_families(&mut self) {
        self.atom_ring_families.clear();
        self.bond_ring_families.clear();
    }

    pub fn preallocate(&mut self, num_atoms: usize, num_bonds: usize) {
        self.atom_members.resize(num_atoms, Vec::new());
        self.bond_members.resize(num_bonds, Vec::new());
    }

    fn check(&self) -> Result<()> {
        if self.initialized {
            Ok(())
        } else {
            Err(not_init())
        }
    }

    /// `addRing`: returns the new number of rings.
    pub fn add_ring(&mut self, atom_indices: Vec<i32>, bond_indices: Vec<i32>) -> Result<u32> {
        self.check()?;
        if atom_indices.len() != bond_indices.len() {
            return Err(RDKitError::precondition("length mismatch"));
        }
        let ring_idx = self.atom_rings.len() as i32;
        for &i in &atom_indices {
            if i as usize >= self.atom_members.len() {
                self.atom_members.resize(i as usize + 1, Vec::new());
            }
            self.atom_members[i as usize].push(ring_idx);
        }
        for &i in &bond_indices {
            if i as usize >= self.bond_members.len() {
                self.bond_members.resize(i as usize + 1, Vec::new());
            }
            self.bond_members[i as usize].push(ring_idx);
        }
        self.atom_rings.push(atom_indices);
        self.bond_rings.push(bond_indices);
        Ok(self.atom_rings.len() as u32)
    }

    pub fn num_rings(&self) -> Result<u32> {
        self.check()?;
        Ok(self.atom_rings.len() as u32)
    }

    pub fn atom_rings(&self) -> &[Vec<i32>] {
        &self.atom_rings
    }
    pub fn bond_rings(&self) -> &[Vec<i32>] {
        &self.bond_rings
    }

    pub fn atom_members(&self, idx: u32) -> Result<&[i32]> {
        self.check()?;
        Ok(self
            .atom_members
            .get(idx as usize)
            .map_or(&[][..], Vec::as_slice))
    }
    pub fn bond_members(&self, idx: u32) -> Result<&[i32]> {
        self.check()?;
        Ok(self
            .bond_members
            .get(idx as usize)
            .map_or(&[][..], Vec::as_slice))
    }

    pub fn num_atom_rings(&self, idx: u32) -> Result<u32> {
        Ok(self.atom_members(idx)?.len() as u32)
    }
    pub fn num_bond_rings(&self, idx: u32) -> Result<u32> {
        Ok(self.bond_members(idx)?.len() as u32)
    }

    pub fn atom_ring_sizes(&self, idx: u32) -> Result<Vec<i32>> {
        Ok(self
            .atom_members(idx)?
            .iter()
            .map(|&r| self.atom_rings[r as usize].len() as i32)
            .collect())
    }
    pub fn bond_ring_sizes(&self, idx: u32) -> Result<Vec<i32>> {
        Ok(self
            .bond_members(idx)?
            .iter()
            .map(|&r| self.bond_rings[r as usize].len() as i32)
            .collect())
    }

    pub fn is_atom_in_ring_of_size(&self, idx: u32, size: u32) -> Result<bool> {
        Ok(self
            .atom_members(idx)?
            .iter()
            .any(|&r| self.atom_rings[r as usize].len() == size as usize))
    }
    pub fn is_bond_in_ring_of_size(&self, idx: u32, size: u32) -> Result<bool> {
        Ok(self
            .bond_members(idx)?
            .iter()
            .any(|&r| self.bond_rings[r as usize].len() == size as usize))
    }

    /// `minAtomRingSize`: 0 when not in a ring. Ties resolve to the first ring.
    pub fn min_atom_ring_size(&self, idx: u32) -> Result<u32> {
        Ok(self
            .atom_members(idx)?
            .iter()
            .map(|&r| self.atom_rings[r as usize].len() as u32)
            .min()
            .unwrap_or(0))
    }
    pub fn min_bond_ring_size(&self, idx: u32) -> Result<u32> {
        Ok(self
            .bond_members(idx)?
            .iter()
            .map(|&r| self.bond_rings[r as usize].len() as u32)
            .min()
            .unwrap_or(0))
    }

    fn in_same_ring_of_size(
        members: &[Vec<i32>],
        rings: &[Vec<i32>],
        i1: u32,
        i2: u32,
        size: u32,
    ) -> bool {
        let (Some(m1), Some(m2)) = (members.get(i1 as usize), members.get(i2 as usize)) else {
            return false;
        };
        // merge walk over the (sorted) ring-index lists, like upstream
        let (mut a, mut b) = (0, 0);
        while a < m1.len() && b < m2.len() {
            if m1[a] < m2[b] {
                a += 1;
            } else if m1[a] > m2[b] {
                b += 1;
            } else if size == 0 || rings[m1[a] as usize].len() == size as usize {
                return true;
            } else {
                a += 1;
                b += 1;
            }
        }
        false
    }

    pub fn are_atoms_in_same_ring_of_size(&self, i1: u32, i2: u32, size: u32) -> Result<bool> {
        self.check()?;
        Ok(Self::in_same_ring_of_size(
            &self.atom_members,
            &self.atom_rings,
            i1,
            i2,
            size,
        ))
    }
    pub fn are_atoms_in_same_ring(&self, i1: u32, i2: u32) -> Result<bool> {
        self.are_atoms_in_same_ring_of_size(i1, i2, 0)
    }
    pub fn are_bonds_in_same_ring_of_size(&self, i1: u32, i2: u32, size: u32) -> Result<bool> {
        self.check()?;
        Ok(Self::in_same_ring_of_size(
            &self.bond_members,
            &self.bond_rings,
            i1,
            i2,
            size,
        ))
    }
    pub fn are_bonds_in_same_ring(&self, i1: u32, i2: u32) -> Result<bool> {
        self.are_bonds_in_same_ring_of_size(i1, i2, 0)
    }

    fn init_fused_rings(&mut self) {
        if self.fused_rings.len() == self.bond_rings.len() {
            return;
        }
        let n = self.bond_rings.len();
        self.fused_rings = vec![vec![false; n]; n];
        for ring_indices in &self.bond_members {
            for i in 0..ring_indices.len().saturating_sub(1) {
                for j in i + 1..ring_indices.len() {
                    let (r1, r2) = (ring_indices[i] as usize, ring_indices[j] as usize);
                    self.fused_rings[r1][r2] = true;
                    self.fused_rings[r2][r1] = true;
                }
            }
        }
    }

    fn ring_bounds(&self, ring: u32, what: &str) -> Result<()> {
        if (ring as usize) < self.fused_rings.len() {
            Ok(())
        } else {
            Err(RDKitError::precondition(format!("{what} out of bounds")))
        }
    }

    pub fn is_ring_fused(&mut self, ring: u32) -> Result<bool> {
        self.init_fused_rings();
        self.ring_bounds(ring, "ringIdx")?;
        Ok(self.fused_rings[ring as usize].iter().any(|&b| b))
    }

    pub fn are_rings_fused(&mut self, r1: u32, r2: u32) -> Result<bool> {
        self.init_fused_rings();
        self.ring_bounds(r1, "ring1Idx")?;
        self.ring_bounds(r2, "ring2Idx")?;
        Ok(self.fused_rings[r1 as usize]
            .get(r2 as usize)
            .copied()
            .unwrap_or(false))
    }

    pub fn num_fused_bonds(&mut self, ring: u32) -> Result<u32> {
        if ring as usize >= self.bond_rings.len() {
            return Err(RDKitError::precondition("ringIdx out of bounds"));
        }
        if self.num_fused_bonds.len() != self.bond_rings.len() {
            self.num_fused_bonds = self
                .bond_rings
                .iter()
                .map(|r| {
                    r.iter()
                        .filter(|&&b| {
                            self.bond_members
                                .get(b as usize)
                                .is_some_and(|m| m.len() > 1)
                        })
                        .count() as u32
                })
                .collect();
        }
        Ok(self.num_fused_bonds[ring as usize])
    }

    pub fn num_fused_ring_neighbors(&mut self, ring: u32) -> Result<u32> {
        Ok(self.fused_ring_neighbors(ring)?.len() as u32)
    }

    pub fn fused_ring_neighbors(&mut self, ring: u32) -> Result<Vec<u32>> {
        self.init_fused_rings();
        self.ring_bounds(ring, "ringIdx")?;
        Ok(self.fused_rings[ring as usize]
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| i as u32)
            .collect())
    }

    pub fn add_ring_family(
        &mut self,
        atom_indices: Vec<i32>,
        bond_indices: Vec<i32>,
    ) -> Result<u32> {
        self.check()?;
        self.atom_ring_families.push(atom_indices);
        self.bond_ring_families.push(bond_indices);
        Ok(self.atom_ring_families.len() as u32)
    }
    pub fn num_ring_families(&self) -> Result<u32> {
        self.check()?;
        Ok(self.atom_ring_families.len() as u32)
    }
    pub fn atom_ring_families(&self) -> &[Vec<i32>] {
        &self.atom_ring_families
    }
    pub fn bond_ring_families(&self) -> &[Vec<i32>] {
        &self.bond_ring_families
    }
}

/// `RingUtils::normalizeRing`: rotate so the smallest index is first and
/// orient so that the second element is smaller than the last.
pub fn normalize_ring(ring: &mut [i32]) {
    if ring.is_empty() {
        return;
    }
    let min_pos = ring
        .iter()
        .enumerate()
        .min_by_key(|(_, v)| **v)
        .map(|(i, _)| i)
        .unwrap_or(0);
    ring.rotate_left(min_pos);
    if ring.len() > 2 && ring[ring.len() - 1] < ring[1] {
        ring[1..].reverse();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_ring_info() {
        let mut ri = RingInfo::new();
        assert!(ri.num_rings().is_err());
        ri.initialize(FindRingType::Sssr);
        // two fused rings sharing atoms 2,3 / bond 2
        ri.add_ring(vec![0, 1, 2, 3], vec![0, 1, 2, 3]).unwrap();
        ri.add_ring(vec![2, 3, 4, 5, 6], vec![2, 4, 5, 6, 7])
            .unwrap();
        assert_eq!(ri.num_rings().unwrap(), 2);
        assert_eq!(ri.num_atom_rings(2).unwrap(), 2);
        assert_eq!(ri.atom_ring_sizes(3).unwrap(), vec![4, 5]);
        assert_eq!(ri.min_atom_ring_size(3).unwrap(), 4);
        assert_eq!(ri.min_atom_ring_size(42).unwrap(), 0);
        assert!(ri.are_atoms_in_same_ring(0, 3).unwrap());
        assert!(!ri.are_atoms_in_same_ring(0, 5).unwrap());
        assert!(ri.are_atoms_in_same_ring_of_size(2, 3, 5).unwrap());
        assert!(ri.is_ring_fused(0).unwrap());
        assert!(ri.are_rings_fused(0, 1).unwrap());
        assert_eq!(ri.num_fused_bonds(1).unwrap(), 1);
        assert_eq!(ri.fused_ring_neighbors(0).unwrap(), vec![1]);
        ri.reset(true);
        assert!(!ri.is_initialized());
    }

    #[test]
    fn normalize() {
        let mut r = vec![3, 4, 1, 2];
        normalize_ring(&mut r);
        assert_eq!(r, vec![1, 2, 3, 4]);
        let mut r = vec![3, 2, 1, 4];
        normalize_ring(&mut r);
        assert_eq!(r, vec![1, 2, 3, 4]);
    }
}
