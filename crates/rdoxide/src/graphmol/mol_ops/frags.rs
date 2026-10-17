//! Fragment handling: `MolOps::getMolFrags`.
//!
//! Ported from RDKit `Code/GraphMol/MolOps.cpp` (2026.09.1).
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!

use super::{sanitize_flags::SANITIZE_ALL, sanitize_mol};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::BondStereo;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::subset::{SubsetOptions, copy_mol_subset};
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
                if mol.substance_groups().iter().any(|sg| {
                    sg.atoms()
                        .iter()
                        .chain(sg.parent_atoms())
                        .any(|&a| in_frag(a))
                }) {
                    return true;
                }
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
                let atoms: Vec<u32> = comp.iter().map(|&a| a as u32).collect();
                let options = SubsetOptions {
                    clear_computed_props: true,
                    copy_coordinates: copy_conformers,
                    ..SubsetOptions::default()
                };
                res.push(copy_mol_subset(mol, &atoms, &options)?.0);
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

/// `MolOps::getMolFragsWithQuery`: groups the atoms of `mol` by the value of
/// `query` (optionally restricted to, or excluding, the values of
/// `white_list`) and returns one molecule, with the bonds between its atoms
/// and their coordinates, for each value.
pub fn get_mol_frags_with_query<K: Ord + Clone>(
    mol: &ROMol,
    query: impl Fn(&ROMol, &crate::graphmol::atom::Atom) -> K,
    sanitize_frags: bool,
    white_list: Option<&[K]>,
    negate_list: bool,
) -> Result<std::collections::BTreeMap<K, RWMol>> {
    use crate::graphmol::conformer::Conformer;
    use std::collections::BTreeMap;

    let n = mol.num_atoms();
    let mut assignments: Vec<Option<K>> = vec![None; n as usize];
    let mut ids: Vec<i32> = vec![-1; n as usize];
    let mut res: BTreeMap<K, RWMol> = BTreeMap::new();
    for i in 0..n {
        let atom = mol.atom(i);
        let key = query(mol, atom);
        if let Some(list) = white_list {
            let found = list.contains(&key);
            if found == negate_list {
                continue;
            }
        }
        let frag = res.entry(key.clone()).or_default();
        ids[i as usize] = frag.add_atom(atom.clone(), false) as i32;
        assignments[i as usize] = Some(key.clone());
        // bonds to the atoms already in the same fragment
        for nbr in mol.atom_neighbors(i) {
            if nbr < i && assignments[nbr as usize].as_ref() == Some(&key) {
                let Some(bidx) = mol.bond_idx_between(i, nbr) else {
                    continue;
                };
                let mut bond = mol.bond(bidx).clone();
                let (b, e) = (bond.begin_atom_idx(), bond.end_atom_idx());
                bond.set_begin_atom_idx_unchecked(ids[b as usize] as u32);
                bond.set_end_atom_idx_unchecked(ids[e as usize] as u32);
                frag.add_bond_object(bond)?;
            }
        }
    }
    // conformers
    for conf in mol.conformers() {
        for frag in res.values_mut() {
            let mut c = Conformer::new(frag.num_atoms() as usize);
            c.set_id(conf.id());
            c.set_3d(conf.is_3d());
            frag.add_conformer(c, false)?;
        }
        for i in 0..n {
            let Some(key) = &assignments[i as usize] else {
                continue;
            };
            let frag = res.get_mut(key).expect("fragment exists");
            frag.conformer_mut(conf.id() as i32)?
                .set_atom_pos(ids[i as usize] as u32, conf.atom_pos(i)?)?;
        }
    }
    if sanitize_frags {
        for frag in res.values_mut() {
            sanitize_mol(frag, SANITIZE_ALL).map_err(|(_, e)| e)?;
        }
    }
    Ok(res)
}
