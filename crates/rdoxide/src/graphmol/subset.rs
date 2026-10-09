//! Extraction of sub-molecules (`copyMolSubset`): copies selected atoms and
//! bonds together with their conformer coordinates and the substance and
//! stereo groups that
//! touch the selection.
//!
//! Ported from RDKit `Code/GraphMol/Subset.cpp` (2026.09.1).
//!
//! Copyright (C) 2025 Hussein Faara, Brian Kelley and other RDKit Contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!

use std::collections::BTreeMap;

use super::bond::{Bond, BondStereo};
use super::conformer::Conformer;
use super::mol_ops::{sanitize_flags::SANITIZE_ALL, sanitize_mol};
use super::query_ops::{query_atom_from_atom, query_bond_from_bond};
use super::romol::{ROMol, RWMol};
use super::stereo_group::StereoGroup;
use super::substance_group::add_substance_group;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::dict::COMPUTED_PROP_NAME;
use crate::{RDKitError, Result};

/// `SubsetMethod`: what the indices of a path stand for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SubsetMethod {
    /// the given atoms and the bonds between them
    #[default]
    BondsBetweenAtoms,
    /// the given bonds and their atoms
    Bonds,
}

/// `SubsetOptions`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubsetOptions {
    /// sanitize the subset
    pub sanitize: bool,
    /// clear all computed properties of the subset
    pub clear_computed_props: bool,
    /// return the subset with query atoms and bonds
    pub copy_as_query: bool,
    /// copy the coordinates of the conformers
    pub copy_coordinates: bool,
    /// kept for API compatibility; upstream copies every conformer
    pub conformer_idx: u32,
    /// how a path is interpreted; ignored when atoms and bonds are given
    pub method: SubsetMethod,
}

impl Default for SubsetOptions {
    fn default() -> Self {
        SubsetOptions {
            sanitize: false,
            clear_computed_props: false,
            copy_as_query: false,
            copy_coordinates: true,
            conformer_idx: u32::MAX,
            method: SubsetMethod::default(),
        }
    }
}

/// `SubsetInfo`: what was selected and where it ended up.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SubsetInfo {
    pub selected_atoms: Vec<bool>,
    pub selected_bonds: Vec<bool>,
    /// original atom index -> index in the subset
    pub atom_mapping: BTreeMap<u32, u32>,
    /// original bond index -> index in the subset
    pub bond_mapping: BTreeMap<u32, u32>,
}

/// The neighbor of `dbl_atom` that is neither `stereo_atom` nor the other end
/// of the double bond (`getOtherAtomIdx`).
fn other_stereo_neighbor(mol: &ROMol, bond: u32, dbl_atom: u32, stereo_atom: u32) -> Option<u32> {
    let b = mol.bond(bond);
    let (mut ref_atom, mut other) = (b.begin_atom_idx(), b.end_atom_idx());
    if other == dbl_atom {
        std::mem::swap(&mut ref_atom, &mut other);
    }
    mol.atom_neighbors(ref_atom)
        .find(|&n| n != stereo_atom && n != other)
}

/// `handleBondStereo`: translates the stereo atoms of a copied double bond
/// into the subset, moving to the other neighbor where one is missing, and
/// turns E/Z into CIS/TRANS since the CIP ranks may have changed.
fn handle_bond_stereo(
    extracted: &mut Bond,
    ref_mol: &ROMol,
    ref_bond: u32,
    atom_mapping: &BTreeMap<u32, u32>,
) -> Result<()> {
    if extracted.stereo_atoms().len() != 2 {
        return Ok(());
    }
    let rb = ref_mol.bond(ref_bond);
    let atoms = extracted.stereo_atoms().to_vec();
    let mut needs_swap = false;
    let mut resolve = |stereo_atom: i32, end_atom: u32| -> Option<u32> {
        if let Some(&mapped) = atom_mapping.get(&(stereo_atom as u32)) {
            return Some(mapped);
        }
        if ref_mol.atom_degree(end_atom) < 3 {
            return None;
        }
        let other = other_stereo_neighbor(ref_mol, ref_bond, end_atom, stereo_atom as u32)?;
        let mapped = atom_mapping.get(&other).copied()?;
        needs_swap = !needs_swap;
        Some(mapped)
    };
    let first = resolve(atoms[0], rb.begin_atom_idx());
    let second = first.and_then(|_| resolve(atoms[1], rb.end_atom_idx()));
    let (Some(first), Some(second)) = (first, second) else {
        extracted.stereo_atoms_mut().clear();
        extracted.set_stereo(BondStereo::StereoNone)?;
        return Ok(());
    };
    *extracted.stereo_atoms_mut() = vec![first as i32, second as i32];
    let stereo = rb.stereo();
    let flipped = match (needs_swap, stereo) {
        (false, BondStereo::StereoZ) => Some(BondStereo::StereoCis),
        (false, BondStereo::StereoE) => Some(BondStereo::StereoTrans),
        (true, BondStereo::StereoZ | BondStereo::StereoCis) => Some(BondStereo::StereoTrans),
        (true, BondStereo::StereoE | BondStereo::StereoTrans) => Some(BondStereo::StereoCis),
        _ => None,
    };
    if let Some(stereo) = flipped {
        extracted.set_stereo(stereo)?;
    }
    Ok(())
}

fn copy_selected_atoms_and_bonds(
    extracted: &mut RWMol,
    mol: &ROMol,
    info: &mut SubsetInfo,
    options: &SubsetOptions,
) -> Result<()> {
    for atom in mol.atoms() {
        if !info.selected_atoms[atom.idx() as usize] {
            continue;
        }
        let mut copy = if options.copy_as_query {
            query_atom_from_atom(atom)?
        } else {
            atom.copy_detached()
        };
        RDProps::clear_computed_props(&mut copy);
        info.atom_mapping
            .insert(atom.idx(), extracted.add_atom(copy, false));
    }
    for bond in mol.bonds() {
        if !info.selected_bonds[bond.idx() as usize] {
            continue;
        }
        let (Some(&begin), Some(&end)) = (
            info.atom_mapping.get(&bond.begin_atom_idx()),
            info.atom_mapping.get(&bond.end_atom_idx()),
        ) else {
            return Err(RDKitError::Value(
                "copyMolSubset: subset bonds contain atoms not contained in subset atoms".into(),
            ));
        };
        let mut copy = if options.copy_as_query {
            query_bond_from_bond(bond)
        } else {
            bond.copy_detached()
        };
        handle_bond_stereo(&mut copy, mol, bond.idx(), &info.atom_mapping)?;
        copy.set_begin_atom_idx_unchecked(begin);
        copy.set_end_atom_idx_unchecked(end);
        let num_bonds = extracted.add_bond_object(copy)?;
        info.bond_mapping.insert(bond.idx(), num_bonds - 1);
    }
    Ok(())
}

/// Substance groups survive when everything they refer to was selected.
fn copy_selected_substance_groups(
    extracted: &mut RWMol,
    mol: &ROMol,
    info: &SubsetInfo,
) -> Result<()> {
    let all_selected =
        |members: &[u32], selected: &[bool]| members.iter().all(|&m| selected[m as usize]);
    let map = |members: &[u32], mapping: &BTreeMap<u32, u32>| -> Result<Vec<u32>> {
        members
            .iter()
            .map(|m| {
                mapping
                    .get(m)
                    .copied()
                    .ok_or_else(|| RDKitError::invariant("unmapped SubstanceGroup member"))
            })
            .collect()
    };
    for group in mol.substance_groups() {
        if !(all_selected(group.atoms(), &info.selected_atoms)
            && all_selected(group.bonds(), &info.selected_bonds)
            && all_selected(group.parent_atoms(), &info.selected_atoms))
        {
            continue;
        }
        let mut copy = group.clone();
        copy.set_atoms(extracted, map(group.atoms(), &info.atom_mapping)?)?;
        copy.set_parent_atoms(extracted, map(group.parent_atoms(), &info.atom_mapping)?)?;
        copy.set_bonds(extracted, map(group.bonds(), &info.bond_mapping)?)?;
        add_substance_group(extracted, copy);
    }
    Ok(())
}

/// Stereo groups survive when they touch the selection; their members that
/// were not selected are dropped.
fn copy_selected_stereo_groups(
    extracted: &mut RWMol,
    mol: &ROMol,
    info: &SubsetInfo,
) -> Result<()> {
    let touches = |members: &[u32], selected: &[bool]| {
        members.is_empty() || members.iter().any(|&m| selected[m as usize])
    };
    let mut groups = Vec::new();
    for group in mol.stereo_groups() {
        if !(touches(group.atoms(), &info.selected_atoms)
            && touches(group.bonds(), &info.selected_bonds))
        {
            continue;
        }
        let atoms = group
            .atoms()
            .iter()
            .filter_map(|a| info.atom_mapping.get(a).copied())
            .collect();
        let bonds = group
            .bonds()
            .iter()
            .filter_map(|b| info.bond_mapping.get(b).copied())
            .collect();
        let mut copy = StereoGroup::new(group.group_type(), atoms, bonds, group.read_id())?;
        copy.set_write_id(group.write_id());
        groups.push(copy);
    }
    extracted.set_stereo_groups(groups)
}

fn copy_coords(extracted: &mut RWMol, mol: &ROMol, info: &SubsetInfo) -> Result<()> {
    for conf in mol.conformers() {
        let mut copy = Conformer::new(extracted.num_atoms() as usize);
        copy.is_3d = conf.is_3d();
        for (&from, &to) in &info.atom_mapping {
            copy.positions[to as usize] = conf.positions()[from as usize];
        }
        copy.set_id(conf.id());
        extracted.add_conformer(copy, false)?;
    }
    Ok(())
}

/// The molecule-level computed properties of `src` (`copyComputedProps`).
fn copy_computed_props(src: &ROMol, dst: &mut RWMol) {
    dst.update_props(&src.props, false);
    for key in dst.get_prop_list(true, false) {
        if key != COMPUTED_PROP_NAME {
            dst.clear_prop(&key);
        }
    }
}

fn copy_selection(mol: &ROMol, info: &mut SubsetInfo, options: &SubsetOptions) -> Result<RWMol> {
    let mut extracted = RWMol::new();
    copy_selected_atoms_and_bonds(&mut extracted, mol, info, options)?;
    copy_selected_substance_groups(&mut extracted, mol, info)?;
    copy_selected_stereo_groups(&mut extracted, mol, info)?;
    if options.copy_coordinates {
        copy_coords(&mut extracted, mol, info)?;
    }
    if options.sanitize {
        sanitize_mol(&mut extracted, SANITIZE_ALL).map_err(|(_, e)| e)?;
    }
    if options.clear_computed_props {
        extracted.clear_computed_props(true);
    } else {
        copy_computed_props(mol, &mut extracted);
    }
    Ok(extracted)
}

fn select_indices(indices: &[u32], count: u32) -> Result<Vec<bool>> {
    let mut selected = vec![false; count as usize];
    for &idx in indices {
        match selected.get_mut(idx as usize) {
            Some(slot) => *slot = true,
            None => return Err(RDKitError::Index((idx as i32).to_string())),
        }
    }
    Ok(selected)
}

/// `copyMolSubset(mol, atoms, bonds, subsetInfo, options)`: exactly the given
/// atoms and bonds (`options.method` is ignored). Selecting everything gives a
/// plain copy and an empty `SubsetInfo`.
pub fn copy_mol_subset_atoms_bonds(
    mol: &ROMol,
    atoms: &[u32],
    bonds: &[u32],
    options: &SubsetOptions,
) -> Result<(RWMol, SubsetInfo)> {
    if atoms.len() == mol.num_atoms() as usize && bonds.len() == mol.num_bonds() as usize {
        return Ok((ROMol::copy_from(mol, false, -1), SubsetInfo::default()));
    }
    let mut info = SubsetInfo {
        selected_atoms: select_indices(atoms, mol.num_atoms())?,
        selected_bonds: select_indices(bonds, mol.num_bonds())?,
        ..SubsetInfo::default()
    };
    let res = copy_selection(mol, &mut info, options)?;
    Ok((res, info))
}

/// `copyMolSubset(mol, path, subsetInfo, options)`: the atoms (or bonds,
/// depending on `options.method`) in `path`; out-of-range indices are ignored.
pub fn copy_mol_subset(
    mol: &ROMol,
    path: &[u32],
    options: &SubsetOptions,
) -> Result<(RWMol, SubsetInfo)> {
    let (num_atoms, num_bonds) = (mol.num_atoms(), mol.num_bonds());
    let mut info = SubsetInfo {
        selected_atoms: vec![false; num_atoms as usize],
        selected_bonds: vec![false; num_bonds as usize],
        ..SubsetInfo::default()
    };
    match options.method {
        SubsetMethod::BondsBetweenAtoms => {
            for &idx in path.iter().filter(|&&i| i < num_atoms) {
                info.selected_atoms[idx as usize] = true;
            }
            for bond in mol.bonds() {
                info.selected_bonds[bond.idx() as usize] = info.selected_atoms
                    [bond.begin_atom_idx() as usize]
                    && info.selected_atoms[bond.end_atom_idx() as usize];
            }
        }
        SubsetMethod::Bonds => {
            for &idx in path.iter().filter(|&&i| i < num_bonds) {
                info.selected_bonds[idx as usize] = true;
                let bond = mol.bond(idx);
                info.selected_atoms[bond.begin_atom_idx() as usize] = true;
                info.selected_atoms[bond.end_atom_idx() as usize] = true;
            }
        }
    }
    let res = copy_selection(mol, &mut info, options)?;
    Ok((res, info))
}
