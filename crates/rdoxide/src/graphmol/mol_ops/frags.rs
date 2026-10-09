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
