//! Substructure-based molecule transforms (upstream
//! `GraphMol/ChemTransforms/ChemTransforms.cpp`): `deleteSubstructs`,
//! `replaceSubstructs`, `replaceSidechains`, `MurckoDecompose` and
//! `combineMols`; `replaceCore` lives in [`replace_core`] and the bond
//! fragmentation (`MolFragmenter.cpp`) in [`mol_fragmenter`].
//!
//! Ported from RDKit 2026.09.1.
//!
//! Copyright (C) 2006-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

pub mod mol_fragmenter;
pub mod replace_core;

use std::collections::{BTreeMap, BTreeSet};

use crate::geometry::Point3D;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::BondType;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::matrices::get_distance_mat_cached;
use crate::graphmol::mol_ops::frags::get_mol_frags_atoms;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::StereoGroup;
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_warning};

pub use mol_fragmenter::{
    FragmenterBondType, fragment_on_bond_patterns, fragment_on_bonds, fragment_on_brics_bonds,
    fragment_on_some_bonds,
};
pub use replace_core::{replace_core, replace_core_with_match};

/// `details::updateSubMolConfs`: rebuilds the conformers of `res` from those
/// of `mol`, leaving out the removed atoms.
pub(crate) fn update_sub_mol_confs(mol: &ROMol, res: &mut RWMol, removed: &[bool]) -> Result<()> {
    res.clear_conformers();
    for conf in mol.conformers() {
        let mut new_conf = Conformer::new(res.num_atoms() as usize);
        new_conf.set_id(conf.id());
        new_conf.set_3d(conf.is_3d());
        let kept = conf
            .positions()
            .iter()
            .zip(removed)
            .filter(|(_, removed)| !**removed)
            .map(|(p, _)| *p);
        for (i, p) in kept.enumerate() {
            new_conf.set_atom_pos(i as u32, p)?;
        }
        res.add_conformer(new_conf, false)?;
    }
    Ok(())
}

/// `details::copyStereoGroups`: the stereo groups of `mol` restricted to the
/// atoms (`old index -> index in new_mol`) that made it into `new_mol`.
pub(crate) fn copy_stereo_groups(
    atom_map: &BTreeMap<u32, u32>,
    mol: &ROMol,
    new_mol: &mut RWMol,
) -> Result<()> {
    if mol.stereo_groups().is_empty() {
        return Ok(());
    }
    let mut groups = Vec::new();
    for group in mol.stereo_groups() {
        let atoms: Vec<u32> = group
            .atoms()
            .iter()
            .filter_map(|a| atom_map.get(a).copied())
            .collect();
        let bonds: Vec<u32> = group
            .bonds()
            .iter()
            .filter_map(|&b| {
                let bond = mol.bond(b);
                let first = atom_map.get(&bond.begin_atom_idx())?;
                let second = atom_map.get(&bond.end_atom_idx())?;
                new_mol.bond_idx_between(*first, *second)
            })
            .collect();
        if !atoms.is_empty() {
            groups.push(StereoGroup::new(
                group.group_type(),
                atoms,
                bonds,
                group.read_id(),
            )?);
        }
    }
    new_mol.set_stereo_groups(groups)
}

/// The first match of `query` in `mol`, as `SubstructMatch(mol, query,
/// matchVect, recursionPossible, useChirality)` finds it.
fn first_match(mol: &ROMol, query: &ROMol, use_chirality: bool) -> Result<Option<MatchVect>> {
    let params = SubstructMatchParameters {
        use_chirality,
        max_matches: 1,
        ..SubstructMatchParameters::default()
    };
    Ok(substruct_match(mol, query, &params)?.into_iter().next())
}

/// `deleteSubstructs(mol, query, onlyFrags, useChirality)`: a copy of `mol`
/// without the atoms matched by `query`. With `only_frags` an atom is removed
/// only when a whole fragment is matched.
pub fn delete_substructs(
    mol: &ROMol,
    query: &ROMol,
    only_frags: bool,
    use_chirality: bool,
) -> Result<RWMol> {
    let mut res = mol.clone();
    let params = SubstructMatchParameters {
        use_chirality,
        ..SubstructMatchParameters::default()
    };
    let matches = substruct_match(&res, query, &params)?;
    if matches.is_empty() {
        return Ok(res);
    }
    let matches: Vec<Vec<i32>> = matches
        .into_iter()
        .map(|m| m.into_iter().map(|(_, atom)| atom).collect())
        .collect();

    let mut to_delete: BTreeSet<i32> = BTreeSet::new();
    if only_frags {
        for mut frag in get_mol_frags_atoms(&res) {
            frag.sort_unstable();
            let whole_frag_matched = matches.iter().find(|m| {
                let mut sorted = (*m).clone();
                sorted.sort_unstable();
                sorted == frag
            });
            if let Some(m) = whole_frag_matched {
                to_delete.extend(m);
            }
        }
    } else {
        for m in &matches {
            to_delete.extend(m);
        }
    }
    if to_delete.is_empty() {
        return Ok(res);
    }

    res.begin_batch_edit()?;
    let mut removed = vec![false; mol.num_atoms() as usize];
    for &idx in &to_delete {
        removed[idx as usize] = true;
        res.remove_atom(idx as u32)?;
    }
    res.commit_batch_edit()?;
    update_sub_mol_confs(mol, &mut res, &removed)?;
    res.clear_computed_props(true);
    // update our properties, but allow unhappiness:
    res.update_property_cache(false)?;
    Ok(res)
}

/// `combineMols(mol1, mol2, offset)`: the atoms of both molecules in one
/// molecule. The coordinates of `mol2` are shifted by `offset` in every
/// conformer that both molecules have.
pub fn combine_mols(mol1: &ROMol, mol2: &ROMol, offset: Point3D) -> Result<RWMol> {
    let mut res = mol1.clone();
    let n_atoms1 = res.num_atoms();
    res.insert_mol(mol2)?;

    if mol1.num_conformers() > 0 && mol2.num_conformers() > 0 {
        if mol1.num_conformers() != mol2.num_conformers() {
            rd_log_warning!("combineMols: molecules have unequal numbers of conformers");
        }
        for conf1 in res.conformers_mut() {
            match mol2.conformer(conf1.id() as i32) {
                Ok(conf2) => {
                    for (i, &p) in conf2.positions().iter().enumerate() {
                        conf1.set_atom_pos(i as u32 + n_atoms1, p + offset)?;
                    }
                }
                Err(_) => {
                    rd_log_warning!("combineMols: conformer id {} not found in mol2", conf1.id())
                }
            }
        }
    }
    res.clear_computed_props(true);
    Ok(res)
}

/// `replaceSubstructs(mol, query, replacement, replaceAll,
/// replacementConnectionPoint, useChirality)`: copies of `mol` in which the
/// atoms matched by `query` are replaced by `replacement`, bonded through its
/// atom `replacement_connection_point` to what the first matched atom was
/// bonded to. With `replace_all` there is a single result with every match
/// replaced, otherwise one result per match.
///
/// As upstream, only the `replace_all` results get their conformers updated
/// and their cached properties cleared.
pub fn replace_substructs(
    mol: &ROMol,
    query: &ROMol,
    replacement: &ROMol,
    replace_all: bool,
    replacement_connection_point: u32,
    use_chirality: bool,
) -> Result<Vec<RWMol>> {
    if replacement_connection_point >= replacement.num_atoms() {
        return Err(RDKitError::precondition("bad replacementConnectionPoint"));
    }
    let params = SubstructMatchParameters {
        use_chirality,
        ..SubstructMatchParameters::default()
    };
    let matches = substruct_match(mol, query, &params)?;
    // nothing matched: a copy of the starting molecule
    if matches.is_empty() {
        let mut copy = mol.clone();
        copy.clear_computed_props(false);
        return Ok(vec![copy]);
    }

    let mut res: Vec<RWMol> = Vec::new();
    let mut removed = vec![false; mol.num_atoms() as usize];
    let mut del_list: BTreeSet<u32> = BTreeSet::new();
    for fgp_match in &matches {
        let matched: Vec<u32> = fgp_match.iter().map(|&(_, a)| a as u32).collect();
        let sorted: BTreeSet<u32> = matched.iter().copied().collect();

        if !replace_all || res.is_empty() {
            res.push(mol.clone());
        }
        let new_mol = res.last_mut().expect("a result was just pushed");
        // the replacement's atoms are numbered from here on
        let num_orig_atoms = new_mol.num_atoms();
        new_mol.insert_mol(replacement)?;

        // duplicate the bonds from the first matched atom to the remainder
        // of the molecule onto the replacement's connection point
        let nbrs: Vec<u32> = new_mol.atom_neighbors(matched[0]).collect();
        for nbr in nbrs {
            if !sorted.contains(&nbr) {
                let bond_idx = new_mol
                    .bond_idx_between(matched[0], nbr)
                    .ok_or_else(|| RDKitError::invariant("required bond not found"))?;
                let bond_type = new_mol.bond(bond_idx).bond_type();
                new_mol.add_bond(
                    num_orig_atoms + replacement_connection_point,
                    nbr,
                    bond_type,
                )?;
            }
        }

        if replace_all {
            del_list.extend(&sorted);
        } else {
            new_mol.begin_batch_edit()?;
            for &idx in &sorted {
                removed[idx as usize] = true;
                new_mol.remove_atom(idx)?;
            }
            new_mol.commit_batch_edit()?;
        }
    }

    if !del_list.is_empty() {
        if replace_all {
            let new_mol = &mut res[0];
            new_mol.begin_batch_edit()?;
            for &idx in &del_list {
                removed[idx as usize] = true;
                new_mol.remove_atom(idx)?;
            }
            new_mol.commit_batch_edit()?;
        }
        // clear conformers and computed props and do basic updates on the
        // resulting molecules, but allow unhappiness
        for new_mol in &mut res {
            update_sub_mol_confs(mol, new_mol, &removed)?;
            new_mol.clear_computed_props(true);
            new_mol.update_property_cache(false)?;
        }
    }
    Ok(res)
}

/// `replaceSidechains(mol, coreQuery, useChirality)`: the core matched by
/// `core_query` with its sidechains replaced by dummy atoms, labelled (as
/// isotopes) from 1. `None` when the core is not found.
pub fn replace_sidechains(
    mol: &ROMol,
    core_query: &ROMol,
    use_chirality: bool,
) -> Result<Option<RWMol>> {
    let Some(match_v) = first_match(mol, core_query, use_chirality)? else {
        return Ok(None);
    };
    if match_v.is_empty() {
        return Ok(None);
    }
    let n_atoms = mol.num_atoms() as usize;
    let mut matching = vec![false; n_atoms];
    for &(_, a) in &match_v {
        matching[a as usize] = true;
    }

    let mut new_mol = mol.clone();
    let mut keep = vec![false; n_atoms];
    let mut dummies: Vec<u32> = Vec::new();
    for &(q, a) in &match_v {
        let (q, a) = (q as u32, a as u32);
        keep[a as usize] = true;
        // an atom of higher degree than its core atom is an attachment point
        if new_mol.atom_degree(a) > core_query.atom_degree(q) {
            let nbrs: Vec<u32> = new_mol.atom_neighbors(a).collect();
            for nbr in nbrs {
                if matching[nbr as usize] {
                    continue;
                }
                // this neighbor isn't in the match: make it a dummy atom
                keep[nbr as usize] = true;
                dummies.push(nbr);
                if let Some(bond_idx) = new_mol.bond_idx_between(a, nbr) {
                    let bond = new_mol.bond_mut(bond_idx);
                    bond.set_is_aromatic(false);
                    bond.set_bond_type(BondType::Single);
                }
                let atom = new_mol.atom_mut(nbr);
                atom.set_atomic_num(0);
                atom.set_num_explicit_hs(0);
                atom.set_is_aromatic(false);
                atom.set_isotope(dummies.len() as u32);
            }
        }
    }

    let mut removed = vec![false; n_atoms];
    new_mol.begin_batch_edit()?;
    for idx in 0..n_atoms {
        if !keep[idx] {
            new_mol.remove_atom(idx as u32)?;
            removed[idx] = true;
        }
    }
    // bonds between newly added dummies go as well
    for (i, &d1) in dummies.iter().enumerate() {
        for &d2 in &dummies[i + 1..] {
            if new_mol.bond_idx_between(d1, d2).is_some() {
                new_mol.remove_bond(d1, d2)?;
            }
        }
    }
    new_mol.commit_batch_edit()?;

    update_sub_mol_confs(mol, &mut new_mol, &removed)?;
    // clear computed props and do basic updates on the resulting molecule,
    // but allow unhappiness
    new_mol.clear_computed_props(true);
    new_mol.update_property_cache(false)?;
    Ok(Some(new_mol))
}

/// `MurckoDecompose(mol)`: the scaffold of `mol` (ring systems and the
/// linkers between them, plus exocyclic double bonds). Like upstream this
/// caches the topological distance matrix on `mol`.
pub fn murcko_decompose(mol: &mut ROMol) -> Result<RWMol> {
    let mut res = mol.clone();
    let n_atoms = res.num_atoms() as usize;
    if n_atoms == 0 {
        return Ok(res);
    }

    // start by getting the shortest paths matrix:
    get_distance_mat_cached(mol, false, false, true, "")?;
    let Some(RDValue::VecInt(path_mat)) = mol.get_prop(cp::DistanceMatrix_Paths) else {
        return Err(RDKitError::invariant("no shortest paths matrix"));
    };

    let mut keep = vec![false; n_atoms];
    for (i, k) in keep.iter_mut().enumerate() {
        *k = res.ring_info().num_atom_rings(i as u32)? > 0;
    }
    // mark the atoms along the shortest path between each pair of rings as
    // keepers
    let rings = res.ring_info().atom_rings();
    for (i, ring_i) in rings.iter().enumerate() {
        for ring_j in &rings[i + 1..] {
            let mut atom_i = ring_i[0];
            let atom_j = ring_j[0];
            while atom_i != atom_j {
                keep[atom_i as usize] = true;
                atom_i = path_mat[atom_j as usize * n_atoms + atom_i as usize];
                // the disconnected case
                if atom_i < 0 {
                    break;
                }
            }
        }
    }

    let mut removed = vec![false; n_atoms];
    res.begin_batch_edit()?;
    for i in 0..n_atoms as u32 {
        if keep[i as usize] {
            continue;
        }
        let mut remove_it = true;
        // check if the atom has a neighboring keeper:
        let nbrs: Vec<u32> = res.atom_neighbors(i).collect();
        for nbr in nbrs {
            if !keep[nbr as usize] {
                continue;
            }
            let bond_type = res
                .bond_idx_between(i, nbr)
                .map(|b| res.bond(b).bond_type());
            let nbr_atom = res.atom_mut(nbr);
            if bond_type == Some(BondType::Double) {
                remove_it = false;
                break;
            } else if nbr_atom.is_aromatic() && nbr_atom.atomic_num() != 6 {
                // fix aromatic heteroatoms:
                nbr_atom.set_num_explicit_hs(1);
            } else if nbr_atom.is_aromatic()
                && nbr_atom.atomic_num() == 6
                && nbr_atom.formal_charge() == 1
            {
                // fix aromatic carbocations
                nbr_atom.set_num_explicit_hs(1);
            } else if nbr_atom.no_implicit() || nbr_atom.chiral_tag() != ChiralType::Unspecified {
                nbr_atom.set_no_implicit(false);
                nbr_atom.set_num_explicit_hs(0);
                nbr_atom.set_chiral_tag(ChiralType::Unspecified);
            }
        }
        if remove_it {
            res.remove_atom(i)?;
            removed[i as usize] = true;
        }
    }
    res.commit_batch_edit()?;

    update_sub_mol_confs(mol, &mut res, &removed)?;
    res.clear_computed_props(true);
    Ok(res)
}
