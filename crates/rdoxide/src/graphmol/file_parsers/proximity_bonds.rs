//! Bond perception from 3D coordinates ("connect the dots").
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/ProximityBonds.{h,cpp}` (2026.09.1).
//!
//! Copyright (C) 2013 Greg Landrum and other RDKit contributors (originally by
//! Roger Sayle). Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::BondType;
use crate::graphmol::monomer_info::{AtomMonomerInfo, AtomMonomerType};
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

/// The monomer info of `atom` if it is a PDB residue.
fn pdb_info(atom: &Atom) -> Option<&AtomMonomerInfo> {
    atom.monomer_info()
        .filter(|i| i.monomer_type == AtomMonomerType::PdbResidue)
}

fn insertion_code(i: &AtomMonomerInfo) -> &str {
    i.pdb.as_ref().map_or("", |p| p.insertion_code.as_str())
}

fn is_hetero(i: &AtomMonomerInfo) -> bool {
    i.pdb.as_ref().is_some_and(|p| p.is_hetero_atom)
}

/// `SamePDBResidue`
pub fn same_pdb_residue(p: &AtomMonomerInfo, q: &AtomMonomerInfo) -> bool {
    p.residue_number == q.residue_number
        && p.residue_name == q.residue_name
        && p.chain_id == q.chain_id
        && insertion_code(p) == insertion_code(q)
}

/// `IsBlacklistedAtom`: everything except non-metals and metalloids
/// (B-O, Si-S, Ge-Se, Sb-Te).
fn is_blacklisted_atom(atom: &Atom) -> bool {
    let elem = atom.atomic_num();
    !((5..=8).contains(&elem)
        || (14..=16).contains(&elem)
        || (32..=34).contains(&elem)
        || (51..=52).contains(&elem))
}

/// `IsBlacklistedPair`: bonds between different PDB residues that involve a
/// metal, noble gas or halogen (or a water) are not wanted.
pub fn is_blacklisted_pair(beg: &Atom, end: &Atom) -> bool {
    let (Some(bi), Some(ei)) = (pdb_info(beg), pdb_info(end)) else {
        return false;
    };
    !same_pdb_residue(bi, ei)
        && (is_blacklisted_atom(beg)
            || is_blacklisted_atom(end)
            || bi.residue_name == "HOH"
            || ei.residue_name == "HOH")
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
                            && !is_blacklisted_pair(mol.atom(entry.atm), mol.atom(other.atm))
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
        let atom_info = mol.atom(i).monomer_info();
        for &n in &neighbors {
            let d = (p - conf.atom_pos(n)?).length() as f32;
            let same_residue = match (atom_info, mol.atom(n).monomer_info()) {
                (Some(a), Some(b)) => a.residue_number == b.residue_number,
                _ => true,
            };
            if d < best && same_residue {
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

/// `StandardPDBDoubleBond(rescode, atm1, atm2)`: the C=O style double bonds of
/// the standard residues. Atom names are the four-character PDB names; the
/// order of the two atoms does not matter.
fn standard_double_bond(residue: &str, a1: &str, a2: &str) -> bool {
    let pair = |x: &str, y: &str| (a1 == x && a2 == y) || (a1 == y && a2 == x);
    let c_o = || pair(" C  ", " O  ");
    let res = residue.get(..3).unwrap_or(residue);
    match res {
        "ACE" | "ALA" | "CYS" | "GLY" | "ILE" | "LEU" | "LYS" | "MET" | "PRO" | "SER" | "THR"
        | "VAL" => c_o(),
        "ARG" => c_o() || pair(" CZ ", " NH2"),
        "ASN" | "ASP" => c_o() || pair(" CG ", " OD1"),
        "GLN" | "GLU" => c_o() || pair(" CD ", " OE1"),
        "HIS" => c_o() || pair(" CD2", " CG ") || pair(" CE1", " ND1"),
        "PHE" | "TYR" => {
            c_o() || pair(" CD1", " CG ") || pair(" CD2", " CE2") || pair(" CE1", " CZ ")
        }
        "TRP" => {
            c_o()
                || pair(" CD1", " CG ")
                || pair(" CD2", " CE2")
                || pair(" CE3", " CZ3")
                || pair(" CH2", " CZ2")
        }
        "  A" | " DA" => {
            pair(" C6 ", " N1 ")
                || pair(" C2 ", " N3 ")
                || pair(" C4 ", " C5 ")
                || pair(" C8 ", " N7 ")
                || pair(" OP1", " P  ")
        }
        "  G" | " DG" => {
            pair(" C6 ", " O6 ")
                || pair(" C2 ", " N3 ")
                || pair(" C4 ", " C5 ")
                || pair(" C8 ", " N7 ")
                || pair(" OP1", " P  ")
        }
        "  C" | " DC" => {
            pair(" C2 ", " O2 ")
                || pair(" C4 ", " N3 ")
                || pair(" C5 ", " C6 ")
                || pair(" OP1", " P  ")
        }
        "  T" | " DT" | "  U" | " DU" => {
            pair(" C2 ", " O2 ")
                || pair(" C4 ", " O4 ")
                || pair(" C5 ", " C6 ")
                || pair(" OP1", " P  ")
        }
        _ => false,
    }
}

/// `StandardPDBDoubleBond(mol, beg, end)`
fn is_standard_double_bond(mol: &RWMol, beg: u32, end: u32) -> bool {
    let (Some(bi), Some(ei)) = (pdb_info(mol.atom(beg)), pdb_info(mol.atom(end))) else {
        return false;
    };
    if !same_pdb_residue(bi, ei) || is_hetero(bi) || is_hetero(ei) {
        return false;
    }
    if !standard_double_bond(&bi.residue_name, &bi.name, &ei.name) {
        return false;
    }
    // neither end may already carry a double bond
    [beg, end]
        .iter()
        .all(|&a| mol.atom_bonds(a).all(|b| b.bond_type() != BondType::Double))
}

/// `StandardPDBResidueBondOrders`: assigns double bonds in the standard amino
/// acid and nucleotide residues.
pub fn standard_pdb_residue_bond_orders(mol: &mut RWMol) {
    for idx in 0..mol.num_bonds() {
        let b = mol.bond(idx);
        if b.bond_type() != BondType::Single {
            continue;
        }
        let (beg, end) = (b.begin_atom_idx(), b.end_atom_idx());
        if is_standard_double_bond(mol, beg, end) {
            mol.bond_mut(idx).set_bond_type(BondType::Double);
        }
    }
}
