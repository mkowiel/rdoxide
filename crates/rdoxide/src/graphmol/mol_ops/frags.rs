//! Fragment handling: `MolOps::getMolFrags` and the `copyMolSubset` path it
//! uses for small, simple fragments.
//!
//! Ported from RDKit `Code/GraphMol/MolOps.cpp` and `Code/GraphMol/Subset.cpp`
//! (2026.09.1).
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! SubstanceGroups are not ported yet; their checks are no-ops.

use std::collections::BTreeMap;

use super::{sanitize_flags::SANITIZE_ALL, sanitize_mol};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::BondStereo;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::StereoGroup;
use crate::rdgeneral::RDProps;
use crate::{RDKitError, Result};

/// `getMolFrags(mol, mapping)`: the connected component of every atom
/// (`boost::connected_components`: components are numbered in order of
/// their lowest atom index). Returns the number of fragments.
pub fn get_mol_frags_mapping(mol: &ROMol) -> (u32, Vec<i32>) {
    let n = mol.num_atoms() as usize;
    let mut mapping = vec![-1i32; n];
    let mut n_frags = 0;
    let mut stack = Vec::new();
    for start in 0..n {
        if mapping[start] >= 0 {
            continue;
        }
        mapping[start] = n_frags;
        stack.push(start as u32);
        while let Some(a) = stack.pop() {
            for nbr in mol.atom_neighbors(a) {
                if mapping[nbr as usize] < 0 {
                    mapping[nbr as usize] = n_frags;
                    stack.push(nbr);
                }
            }
        }
        n_frags += 1;
    }
    (n_frags as u32, mapping)
}

/// `getMolFrags(mol, frags)`: atom indices of each fragment.
pub fn get_mol_frags_atoms(mol: &ROMol) -> Vec<Vec<i32>> {
    let (n, mapping) = get_mol_frags_mapping(mol);
    let mut frags = vec![Vec::new(); n as usize];
    for (i, &m) in mapping.iter().enumerate() {
        frags[m as usize].push(i as i32);
    }
    frags
}

/// `getOtherAtomIdx` (Subset.cpp)
fn subset_other_atom_idx(mol: &ROMol, bond: u32, dbl_atom: u32, stereo_atom: u32) -> Option<u32> {
    let b = mol.bond(bond);
    let (mut ref_atom, mut other) = (b.begin_atom_idx(), b.end_atom_idx());
    if other == dbl_atom {
        std::mem::swap(&mut ref_atom, &mut other);
    }
    mol.atom_neighbors(ref_atom)
        .find(|&n| n != stereo_atom && n != other)
}

/// `handleBondStereo` (Subset.cpp)
fn handle_bond_stereo(
    extracted: &mut crate::graphmol::bond::Bond,
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
    let mut resolve = |sa: i32, end_atom: u32| -> Option<u32> {
        if let Some(&m) = atom_mapping.get(&(sa as u32)) {
            return Some(m);
        }
        if ref_mol.atom_degree(end_atom) < 3 {
            return None;
        }
        let o = subset_other_atom_idx(ref_mol, ref_bond, end_atom, sa as u32)?;
        let m = atom_mapping.get(&o).copied()?;
        needs_swap = !needs_swap;
        Some(m)
    };
    let m1 = resolve(atoms[0], rb.begin_atom_idx());
    let m2 = m1.and_then(|_| resolve(atoms[1], rb.end_atom_idx()));
    let (Some(m1), Some(m2)) = (m1, m2) else {
        extracted.stereo_atoms_mut().clear();
        extracted.set_stereo(BondStereo::StereoNone)?;
        return Ok(());
    };
    *extracted.stereo_atoms_mut() = vec![m1 as i32, m2 as i32];
    let rs = rb.stereo();
    if !needs_swap {
        if rs == BondStereo::StereoZ {
            extracted.set_stereo(BondStereo::StereoCis)?;
        } else if rs == BondStereo::StereoE {
            extracted.set_stereo(BondStereo::StereoTrans)?;
        }
    } else if matches!(rs, BondStereo::StereoZ | BondStereo::StereoCis) {
        extracted.set_stereo(BondStereo::StereoTrans)?;
    } else if matches!(rs, BondStereo::StereoE | BondStereo::StereoTrans) {
        extracted.set_stereo(BondStereo::StereoCis)?;
    }
    Ok(())
}

/// `copyMolSubset(mol, atoms, opts)` with `SubsetMethod::BONDS_BETWEEN_ATOMS`,
/// `clearComputedProps = true` and no sanitization.
fn copy_mol_subset_bonds_between_atoms(
    mol: &ROMol,
    atoms: &[i32],
    copy_coordinates: bool,
) -> Result<RWMol> {
    let na = mol.num_atoms() as usize;
    let mut selected_atoms = vec![false; na];
    for &a in atoms {
        if (a as usize) < na {
            selected_atoms[a as usize] = true;
        }
    }
    let selected_bonds: Vec<bool> = mol
        .bonds()
        .iter()
        .map(|b| {
            selected_atoms[b.begin_atom_idx() as usize] && selected_atoms[b.end_atom_idx() as usize]
        })
        .collect();
    // copyMolSubset(mol, atoms, bonds, ...) returns a plain copy for a full selection
    if selected_atoms.iter().all(|&x| x) && selected_bonds.iter().all(|&x| x) {
        // the full-copy shortcut does not clear computed props
        return Ok(ROMol::copy_from(mol, false, -1));
    }
    let mut res = RWMol::new();
    let mut atom_mapping: BTreeMap<u32, u32> = BTreeMap::new();
    let mut bond_mapping: BTreeMap<u32, u32> = BTreeMap::new();
    for a in mol.atoms() {
        if !selected_atoms[a.idx() as usize] {
            continue;
        }
        let mut c = a.copy_detached();
        RDProps::clear_computed_props(&mut c);
        atom_mapping.insert(a.idx(), res.add_atom(c, false));
    }
    for b in mol.bonds() {
        if !selected_bonds[b.idx() as usize] {
            continue;
        }
        let mut c = b.copy_detached();
        handle_bond_stereo(&mut c, mol, b.idx(), &atom_mapping)?;
        c.set_begin_atom_idx_unchecked(atom_mapping[&b.begin_atom_idx()]);
        c.set_end_atom_idx_unchecked(atom_mapping[&b.end_atom_idx()]);
        let n = res.add_bond_object(c)?;
        bond_mapping.insert(b.idx(), n - 1);
    }
    // TODO(SubstanceGroups): copySelectedSubstanceGroups
    let mut sgs = Vec::new();
    for sg in mol.stereo_groups() {
        let sel = |objs: &[u32], s: &[bool]| objs.is_empty() || objs.iter().any(|&o| s[o as usize]);
        if !(sel(sg.atoms(), &selected_atoms) && sel(sg.bonds(), &selected_bonds)) {
            continue;
        }
        let atoms: Vec<u32> = sg
            .atoms()
            .iter()
            .filter_map(|a| atom_mapping.get(a).copied())
            .collect();
        let bonds: Vec<u32> = sg
            .bonds()
            .iter()
            .filter_map(|b| bond_mapping.get(b).copied())
            .collect();
        let mut nsg = StereoGroup::new(sg.group_type(), atoms, bonds, sg.read_id())?;
        nsg.write_id = sg.write_id;
        sgs.push(nsg);
    }
    res.set_stereo_groups(sgs)?;
    if copy_coordinates {
        for conf in mol.conformers() {
            let mut c = Conformer::new(res.num_atoms() as usize);
            c.is_3d = conf.is_3d();
            for (&from, &to) in &atom_mapping {
                c.positions[to as usize] = conf.positions()[from as usize];
            }
            c.set_id(conf.id());
            res.add_conformer(c, false)?;
        }
    }
    res.clear_computed_props(true);
    Ok(res)
}

/// `getMolFrags(mol, sanitizeFrags, frags, fragsMolAtomMapping, copyConformers)`
/// with `asMols=true`: returns the fragments and their atom mappings.
pub fn get_mol_frags_as_mols(
    mol: &ROMol,
    sanitize_frags: bool,
    copy_conformers: bool,
) -> Result<(Vec<RWMol>, Vec<Vec<i32>>)> {
    let (n_frags, frags) = get_mol_frags_mapping(mol);
    let mut res = Vec::new();
    let mut mapping = Vec::new();
    if n_frags == 1 {
        res.push(ROMol::copy_from(mol, false, -1));
        mapping.push((0..mol.num_atoms() as i32).collect());
    } else {
        for i in 0..n_frags as i32 {
            let comp: Vec<i32> = (0..mol.num_atoms() as i32)
                .filter(|&a| frags[a as usize] == i)
                .collect();
            let in_frag = |a: u32| frags[a as usize] == i;
            let challenging = || {
                for &idx in &comp {
                    let t = mol.atom(idx as u32).chiral_tag();
                    if t != ChiralType::Unspecified && t != ChiralType::Other {
                        return true;
                    }
                    for b in mol.atom_bonds(idx as u32) {
                        if in_frag(b.other_atom_idx(idx as u32).unwrap_or(u32::MAX))
                            && !matches!(b.stereo(), BondStereo::StereoNone | BondStereo::StereoAny)
                        {
                            return true;
                        }
                    }
                }
                // TODO(SubstanceGroups)
                for sg in mol.stereo_groups() {
                    if sg.atoms().iter().any(|&a| in_frag(a)) {
                        return true;
                    }
                    if sg.bonds().iter().any(|&b| {
                        let bd = mol.bond(b);
                        in_frag(bd.begin_atom_idx()) && in_frag(bd.end_atom_idx())
                    }) {
                        return true;
                    }
                }
                false
            };
            if comp.len() == 1 || (n_frags > 3 && !challenging()) {
                res.push(copy_mol_subset_bonds_between_atoms(
                    mol,
                    &comp,
                    copy_conformers,
                )?);
            } else {
                let mut frag = ROMol::copy_from(mol, false, -1);
                frag.begin_batch_edit()?;
                for idx in 0..mol.num_atoms() {
                    if !in_frag(idx) {
                        frag.remove_atom(idx)?;
                    }
                }
                frag.commit_batch_edit()?;
                res.push(frag);
            }
            mapping.push(comp);
        }
    }
    if !copy_conformers {
        for f in &mut res {
            f.clear_conformers();
        }
    }
    if sanitize_frags {
        for f in &mut res {
            sanitize_mol(f, SANITIZE_ALL).map_err(|(_, e)| e)?;
        }
    }
    if res.is_empty() && mol.num_atoms() > 0 {
        return Err(RDKitError::invariant("no fragments found"));
    }
    Ok((res, mapping))
}
