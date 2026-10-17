//! Conversion between haptic bonds (a dummy atom with a dative bond to the
//! metal and `_MolFileBondEndPts`) and plain dative bonds from each atom of
//! the ligand (`hapticBondsToDative`/`dativeBondsToHaptic` of `MolOps.cpp`,
//! RDKit 2026.09.1).
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, VecDeque};

use crate::Result;
use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::queries::make_atom_null_query;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::common_properties as cp;

/// `details::hapticBondEndpoints`: the (zero-based) atom indices of the
/// `_MolFileBondEndPts` property, whose first number is the count and the
/// others count from 1.
pub fn haptic_bond_endpoints(bond: &Bond) -> Vec<i32> {
    let Some(RDValue::String(endpts)) = bond.get_prop(cp::_MolFileBondEndPts) else {
        return Vec::new();
    };
    let Some(inner) = endpts.strip_prefix('(').and_then(|s| s.strip_suffix(')')) else {
        return Vec::new();
    };
    inner
        .split(' ')
        .filter(|t| !t.is_empty())
        .skip(1)
        .filter_map(|t| t.parse::<f64>().ok())
        .map(|v| (v - 1.0) as i32)
        .collect()
}

/// `MolOps::hapticBondsToDative`
pub fn haptic_bonds_to_dative(mol: &mut RWMol) -> Result<()> {
    mol.begin_batch_edit()?;
    for bi in 0..mol.num_bonds() {
        let bond = mol.bond(bi);
        if bond.bond_type() != BondType::Dative {
            continue;
        }
        let oats = haptic_bond_endpoints(bond);
        if oats.is_empty() {
            continue;
        }
        let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let (dummy, metal) = if mol.atom(begin).atomic_num() == 0 {
            (begin, end)
        } else if mol.atom(end).atomic_num() == 0 {
            (end, begin)
        } else {
            continue;
        };
        for oat in oats {
            if oat >= 0 && (oat as u32) < mol.num_atoms() {
                mol.add_bond(oat as u32, metal, BondType::Dative)?;
            }
        }
        mol.remove_atom(dummy)?;
    }
    mol.commit_batch_edit()
}

/// The sets of more than one contiguous atoms of `all_ats`.
fn contiguous_atoms(mol: &ROMol, all_ats: &[u32]) -> Vec<Vec<u32>> {
    let mut contig: Vec<Vec<u32>> = Vec::new();
    let mut done = vec![false; mol.num_atoms() as usize];
    let mut in_all = vec![false; mol.num_atoms() as usize];
    for &a in all_ats {
        in_all[a as usize] = true;
    }
    for &start in all_ats {
        if done[start as usize] {
            continue;
        }
        contig.push(Vec::new());
        let mut to_do: VecDeque<u32> = VecDeque::from([start]);
        while let Some(next) = to_do.pop_front() {
            if !done[next as usize] {
                done[next as usize] = true;
                contig.last_mut().expect("just pushed").push(next);
            }
            for nbr in mol.atom_neighbors(next) {
                if in_all[nbr as usize] && !done[nbr as usize] {
                    to_do.push_back(nbr);
                }
            }
        }
        if contig.last().is_some_and(|c| c.len() < 2) {
            contig.pop();
        }
    }
    contig
}

/// Adds a dummy atom centred on `haptic_atoms`, with a dative bond from it to
/// the metal.
fn add_haptic_bond(mol: &mut RWMol, metal_idx: u32, haptic_atoms: &[u32]) -> Result<()> {
    let mut dummy = Atom::new(0);
    dummy.set_query(Some(make_atom_null_query()));
    let dummy_idx = mol.add_atom(dummy, true);
    for conf in mol.conformers_mut() {
        let mut pos = Point3D::default();
        for &ha in haptic_atoms {
            pos += conf.atom_pos(ha)?;
        }
        pos /= haptic_atoms.len() as f64;
        conf.set_atom_pos(dummy_idx, pos)?;
    }
    let bond_idx = mol.add_bond(dummy_idx, metal_idx, BondType::Dative)? - 1;
    // the first number is the count, the rest count from 1
    let numbers: Vec<String> = haptic_atoms.iter().map(|ha| (ha + 1).to_string()).collect();
    let endpts = format!("({} {})", haptic_atoms.len(), numbers.join(" "));
    let bond = mol.bond_mut(bond_idx);
    bond.props.set_val(cp::_MolFileBondEndPts, endpts);
    bond.props
        .set_val(cp::_MolFileBondAttach, "ALL".to_string());
    Ok(())
}

/// `MolOps::dativeBondsToHaptic`
pub fn dative_bonds_to_haptic(mol: &mut RWMol) -> Result<()> {
    // the atoms with a dative bond to each (assumed metal) end atom
    let mut dative_atoms: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for b in mol.bonds() {
        if b.bond_type() == BondType::Dative {
            dative_atoms
                .entry(b.end_atom_idx())
                .or_default()
                .push(b.begin_atom_idx());
        }
    }
    mol.begin_batch_edit()?;
    for (metal, atoms) in &dative_atoms {
        // each contiguous set is the end points of a haptic bond to the metal
        for ca in contiguous_atoms(mol, atoms) {
            add_haptic_bond(mol, *metal, &ca)?;
            for cat in ca {
                mol.remove_bond(*metal, cat)?;
            }
        }
    }
    mol.commit_batch_edit()
}
