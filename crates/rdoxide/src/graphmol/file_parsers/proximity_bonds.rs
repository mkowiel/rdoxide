//! Bond perception from 3D coordinates ("connect the dots").
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/ProximityBonds.{h,cpp}` (2026.09.1);
//! the PDB residue handling is not ported because atoms carry no monomer
//! info yet.
//!
//! Copyright (C) 2013 Greg Landrum and other RDKit contributors (originally by
//! Roger Sayle). Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::graphmol::bond::BondType;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::RWMol;

/// do not bond hydrogens to each other
pub const CTD_IGNORE_H_H_CONTACTS: u32 = 0x1;
/// remove the H-H bond of a multivalent hydrogen when it also has a non-H neighbor
pub const CTD_QUICKREMOVE_H_H_CONTACTS: u32 = 0x2;

const EXTDIST: f64 = 0.45;
const MINDIST2: f64 = 0.16;
const MAXDIST: f64 = 5.45;
const MAXDIST2: f64 = 29.7025;

const HASHSIZE: usize = 1024;
const HASHMASK: i32 = 1023;
const HASHX: i32 = 571;
const HASHY: i32 = 127;
const HASHZ: i32 = 3;

/// Coordinates and radius are single precision as upstream.
#[derive(Clone, Copy, Default)]
struct ProximityEntry {
    x: f32,
    y: f32,
    z: f32,
    r: f32,
    atm: u32,
    hash: i32,
    next: Option<usize>,
    elem: u32,
}

fn is_bonded(p: &ProximityEntry, q: &ProximityEntry, flags: u32) -> bool {
    if flags & CTD_IGNORE_H_H_CONTACTS != 0 && p.elem == 1 && q.elem == 1 {
        return false;
    }
    let dx = f64::from(p.x) - f64::from(q.x);
    let mut dist2 = dx * dx;
    if dist2 > MAXDIST2 {
        return false;
    }
    let dy = f64::from(p.y) - f64::from(q.y);
    dist2 += dy * dy;
    if dist2 > MAXDIST2 {
        return false;
    }
    let dz = f64::from(p.z) - f64::from(q.z);
    dist2 += dz * dz;
    if !(MINDIST2..=MAXDIST2).contains(&dist2) {
        return false;
    }
    let radius = f64::from(p.r) + f64::from(q.r) + EXTDIST;
    dist2 <= radius * radius
}

/// `ConnectTheDots`: adds single bonds between atoms that are closer than the
/// sum of their covalent radii (plus a tolerance) in the default conformer.
/// Does nothing without a conformer.
pub fn connect_the_dots(mol: &mut RWMol, flags: u32) -> Result<()> {
    if mol.num_conformers() == 0 {
        return Ok(());
    }
    let table = PeriodicTable::get();
    let count = mol.num_atoms() as usize;
    let conf = mol.conformer(-1)?.clone();
    let mut hash_table = [None::<usize>; HASHSIZE];
    let mut entries: Vec<ProximityEntry> = Vec::with_capacity(count);

    for i in 0..count {
        let elem = mol.atom(i as u32).atomic_num();
        let p = conf.atom_pos(i as u32)?;
        let mut entry = ProximityEntry {
            x: p.x as f32,
            y: p.y as f32,
            z: p.z as f32,
            r: table.get_rcovalent(elem)? as f32,
            atm: i as u32,
            elem,
            ..ProximityEntry::default()
        };
        let hash = HASHX * (p.x / MAXDIST) as i32
            + HASHY * (p.y / MAXDIST) as i32
            + HASHZ * (p.z / MAXDIST) as i32;
        for dx in [-HASHX, 0, HASHX] {
            for dy in [-HASHY, 0, HASHY] {
                for dz in [-HASHZ, 0, HASHZ] {
                    let probe = hash + dx + dy + dz;
                    let mut list = hash_table[(probe & HASHMASK) as usize];
                    while let Some(j) = list {
                        let other = entries[j];
                        if other.hash == probe
                            && is_bonded(&entry, &other, flags)
                            && mol.get_bond_between_atoms(entry.atm, other.atm)?.is_none()
                        {
                            mol.add_bond(entry.atm, other.atm, BondType::Single)?;
                        }
                        list = other.next;
                    }
                }
            }
        }
        let slot = (hash & HASHMASK) as usize;
        entry.next = hash_table[slot];
        hash_table[slot] = Some(i);
        entry.hash = hash;
        entries.push(entry);
    }
    cleanup_multivalent_hydrogens(mol, flags)
}

/// ConnectTheDots can attach a hydrogen to several atoms: keep only the
/// closest neighbor.
fn cleanup_multivalent_hydrogens(mol: &mut RWMol, flags: u32) -> Result<()> {
    let conf = mol.conformer(-1)?.clone();
    for i in 0..mol.num_atoms() {
        if mol.atom(i).atomic_num() != 1 || mol.atom_bonds(i).count() <= 1 {
            continue;
        }
        if flags & CTD_QUICKREMOVE_H_H_CONTACTS != 0 {
            // an H neighbor and a non-H neighbor: remove the bond to the H
            let mut bond_to_h = None;
            let mut bond_to_non_h = None;
            for bond in mol.atom_bonds(i) {
                let other = bond.other_atom_idx(i)?;
                if mol.atom(other).atomic_num() == 1 {
                    bond_to_h = Some((bond.begin_atom_idx(), bond.end_atom_idx()));
                } else {
                    bond_to_non_h = Some(bond.idx());
                }
            }
            if let (Some((a, b)), Some(_)) = (bond_to_h, bond_to_non_h) {
                mol.remove_bond(a, b)?;
            }
        }
        if mol.atom_bonds(i).count() == 1 {
            continue;
        }
        // cut all but the shortest bond
        let p = conf.atom_pos(i)?;
        let neighbors: Vec<u32> = mol.atom_neighbors(i).collect();
        let mut best = 10000.0f32;
        let mut best_idx = mol.num_atoms() + 1;
        for &n in &neighbors {
            let d = (p - conf.atom_pos(n)?).length() as f32;
            if d < best {
                best = d;
                best_idx = n;
            }
        }
        mol.begin_batch_edit()?;
        for &n in &neighbors {
            if n == best_idx {
                if let Some(bond) = mol.get_bond_between_atoms(i, n)?.map(|b| b.idx()) {
                    mol.bond_mut(bond).set_bond_type(BondType::Single);
                }
            } else {
                mol.remove_bond(i, n)?;
            }
        }
        mol.commit_batch_edit()?;
    }
    Ok(())
}
