//! Ring and element helpers of the MMFF atom typer (upstream
//! `ForceFieldHelpers/MMFF/AtomTyper.cpp`).

use std::collections::{BTreeMap, BTreeSet};

use crate::error::Result;
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::ROMol;

/// Periodic table row of an element, 0 for hydrogen.
pub fn periodic_table_row(atomic_num: u32) -> u32 {
    match atomic_num {
        3..=10 => 1,
        11..=18 => 2,
        19..=36 => 3,
        37..=54 => 4,
        _ => 0,
    }
}

/// Periodic table row of an element counted from 1 for helium; hydrogen has
/// row 0 and transition metals have their row multiplied by 10.
pub fn periodic_table_row_hl(atomic_num: u32) -> u32 {
    let row = match atomic_num {
        2 => 1,
        3..=10 => 2,
        11..=18 => 3,
        19..=36 => 4,
        37..=54 => 5,
        _ => 0,
    };
    if matches!(atomic_num, 21..=30 | 39..=48) {
        row * 10
    } else {
        row
    }
}

/// Whether all the bonds along the ring (except the closing one) are aromatic.
pub fn is_ring_aromatic(mol: &ROMol, ring: &[i32]) -> Result<bool> {
    for pair in ring.windows(2) {
        let aromatic = mol
            .get_bond_between_atoms(pair[0] as u32, pair[1] as u32)?
            .is_some_and(|b| b.bond_type() == BondType::Aromatic);
        if !aromatic {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether the atom is in an aromatic ring of the given size.
pub fn is_atom_in_aromatic_ring_of_size(mol: &ROMol, idx: u32, size: u32) -> Result<bool> {
    if !mol.atom(idx).is_aromatic() {
        return Ok(false);
    }
    for ring in mol.ring_info().atom_rings() {
        if ring.len() == size as usize
            && ring.contains(&(idx as i32))
            && is_ring_aromatic(mol, ring)?
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether the atom is a nitrogen with an oxygen of degree 1 attached.
pub fn is_atom_n_oxide(mol: &ROMol, idx: u32) -> Result<bool> {
    if mol.atom(idx).atomic_num() != 7 || mol.total_degree(idx)? < 3 {
        return Ok(false);
    }
    for nbr in mol.atom_neighbors(idx) {
        if mol.atom(nbr).atomic_num() == 8 && mol.total_degree(nbr)? == 1 {
            return Ok(true);
        }
    }
    Ok(false)
}

fn bonded(mol: &ROMol, a: u32, b: u32) -> bool {
    mol.bond_idx_between(a, b).is_some()
}

/// Neighbors of `idx` other than `skip`.
fn neighbors_except(mol: &ROMol, idx: u32, skip: u32) -> BTreeSet<u32> {
    mol.atom_neighbors(idx).filter(|&n| n != skip).collect()
}

/// 3 or 4 if the angle idx1-idx2-idx3 is in a ring of that size, else 0.
pub fn is_angle_in_ring_of_size_3_or_4(mol: &ROMol, idx1: u32, idx2: u32, idx3: u32) -> u32 {
    if !(bonded(mol, idx1, idx2) && bonded(mol, idx2, idx3)) {
        return 0;
    }
    if bonded(mol, idx3, idx1) {
        return 3;
    }
    let s1 = neighbors_except(mol, idx1, idx2);
    let s2 = neighbors_except(mol, idx3, idx2);
    if s1.intersection(&s2).next().is_some() {
        4
    } else {
        0
    }
}

/// 4 or 5 if the torsion idx1-idx2-idx3-idx4 is in a ring of that size, else 0.
pub fn is_torsion_in_ring_of_size_4_or_5(
    mol: &ROMol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
) -> u32 {
    if !(bonded(mol, idx1, idx2) && bonded(mol, idx2, idx3) && bonded(mol, idx3, idx4)) {
        return 0;
    }
    if bonded(mol, idx4, idx1) {
        return 4;
    }
    let s1 = neighbors_except(mol, idx1, idx2);
    let s2 = neighbors_except(mol, idx4, idx3);
    if s1.intersection(&s2).next().is_some() {
        5
    } else {
        0
    }
}

/// Whether all the atoms are in one ring of the given size.
pub fn are_atoms_in_same_ring_of_size(mol: &ROMol, ring_size: u32, atoms: &[u32]) -> bool {
    mol.ring_info()
        .atom_rings()
        .iter()
        .filter(|r| r.len() == ring_size as usize)
        .any(|r| atoms.iter().all(|&a| r.contains(&(a as i32))))
}

/// Whether both atoms are aromatic and in one aromatic ring.
pub fn are_atoms_in_same_aromatic_ring(mol: &ROMol, idx1: u32, idx2: u32) -> Result<bool> {
    if !(mol.atom(idx1).is_aromatic() && mol.atom(idx2).is_aromatic()) {
        return Ok(false);
    }
    for ring in mol.ring_info().atom_rings() {
        if ring.contains(&(idx1 as i32))
            && ring.contains(&(idx2 as i32))
            && is_ring_aromatic(mol, ring)?
        {
            return Ok(true);
        }
    }
    Ok(false)
}

const IS_AROMATIC_BIT: u32 = 1 << 31;

#[derive(Debug, Default)]
struct RingMembership {
    in_aromatic_ring: bool,
    /// ring indices, with `IS_AROMATIC_BIT` set for aromatic rings
    rings: BTreeSet<u32>,
}

/// The rings each atom belongs to, by ring size (upstream
/// `RingMembershipSize`).
#[derive(Debug)]
pub struct RingMembershipSize {
    by_size: BTreeMap<usize, BTreeMap<u32, RingMembership>>,
}

impl RingMembershipSize {
    pub fn new(mol: &ROMol) -> Result<Self> {
        let mut by_size: BTreeMap<usize, BTreeMap<u32, RingMembership>> = BTreeMap::new();
        for (ring_idx, ring) in mol.ring_info().atom_rings().iter().enumerate() {
            let aromatic = is_ring_aromatic(mol, ring)?;
            let mut flagged = ring_idx as u32;
            if aromatic {
                flagged |= IS_AROMATIC_BIT;
            }
            let members = by_size.entry(ring.len()).or_default();
            for &atom in ring {
                let m = members.entry(atom as u32).or_default();
                m.rings.insert(flagged);
                if aromatic {
                    m.in_aromatic_ring = true;
                }
            }
        }
        Ok(RingMembershipSize { by_size })
    }

    pub fn is_atom_in_aromatic_ring_of_size(&self, idx: u32, size: usize) -> bool {
        self.by_size
            .get(&size)
            .and_then(|m| m.get(&idx))
            .is_some_and(|m| m.in_aromatic_ring)
    }

    pub fn are_atoms_in_same_aromatic_ring(&self, idx1: u32, idx2: u32) -> bool {
        self.by_size
            .values()
            .any(|members| match (members.get(&idx1), members.get(&idx2)) {
                (Some(a), Some(b)) => a
                    .rings
                    .intersection(&b.rings)
                    .any(|r| r & IS_AROMATIC_BIT != 0),
                _ => false,
            })
    }

    pub fn are_atoms_in_same_ring_of_size(&self, size: usize, atoms: &[u32]) -> bool {
        let Some(members) = self.by_size.get(&size) else {
            return false;
        };
        let Some(first) = atoms.first().and_then(|a| members.get(a)) else {
            return false;
        };
        let mut common = first.rings.clone();
        let mut same = false;
        for atom in &atoms[1..] {
            if common.is_empty() {
                break;
            }
            same = false;
            let Some(m) = members.get(atom) else { break };
            common = common.intersection(&m.rings).copied().collect();
            same = !common.is_empty();
        }
        same
    }
}
