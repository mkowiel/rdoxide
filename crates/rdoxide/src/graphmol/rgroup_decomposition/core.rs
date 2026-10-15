//! A core of the decomposition and its matching against molecules (upstream
//! `RGroupDecomposition/RGroupCore.{h,cpp}`).

use std::collections::{BTreeMap, BTreeSet};

use super::params::RGroupDecompositionParameters;
use super::params::rgroup_labels::ISOTOPE_LABELS;
use super::utils::{
    Labelling, RLABEL, RLABEL_CORE_INDEX, RLABEL_TYPE, UNLABELED_CORE_ATTACHMENT, int_prop,
    is_atom_with_multiple_neighbors_or_not_dummy_rgroup_attachment, is_dummy_rgroup_attachment,
    is_user_rlabel, set_int_prop,
};
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::BondStereo;
use crate::graphmol::chem_transforms::{copy_stereo_groups, update_sub_mol_confs};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_ops::sanitize_flags::{SANITIZE_CLEANUP, SANITIZE_SYMMRINGS};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::substruct::{FinalCheck, MatchVect, SubstructMatchParameters, bond_compat};
use crate::rdgeneral::RDProps;
use crate::{RDKitError, Result};

/// All the sequences taking one element of each list in `v`; with
/// `allow_duplicates` off an element is used at most once per sequence.
fn cartesian_product(v: &[Vec<i32>], allow_duplicates: bool) -> Vec<Vec<i32>> {
    let mut s: Vec<Vec<i32>> = vec![Vec::new()];
    for u in v {
        let mut r = Vec::new();
        for x in &s {
            for &y in u {
                if allow_duplicates || !x.contains(&y) {
                    let mut next = x.clone();
                    next.push(y);
                    r.push(next);
                }
            }
        }
        if r.is_empty() {
            return r;
        }
        s = r;
    }
    s
}

/// `RCore`
#[derive(Debug, Clone, Default)]
pub struct RCore {
    pub core: RWMol,
    /// the core without the terminal R-label dummies, matched against molecules
    pub matching_mol: RWMol,
    pub labelled_core: RWMol,
    core_atoms_with_user_labels: Vec<bool>,
    pub number_user_rgroups: usize,
    /// terminal user-labelled dummy atom of the core -> its neighbor
    terminal_rgroup_atom_to_neighbor: BTreeMap<i32, i32>,
}

impl RCore {
    pub fn new(core: RWMol) -> Result<Self> {
        let mut rcore = RCore {
            core,
            ..Default::default()
        };
        rcore.init()?;
        Ok(rcore)
    }

    /// `init`
    pub fn init(&mut self) -> Result<()> {
        self.find_indices_with_rlabel();
        self.number_user_rgroups = self
            .core_atoms_with_user_labels
            .iter()
            .filter(|&&b| b)
            .count();
        self.build_matching_mol()
    }

    fn find_indices_with_rlabel(&mut self) {
        self.core_atoms_with_user_labels = self
            .core
            .atoms()
            .iter()
            .map(|a| int_prop(a, RLABEL).is_some_and(|l| l > 0))
            .collect();
    }

    pub fn is_core_atom_user_labelled(&self, idx: i32) -> bool {
        self.core_atoms_with_user_labels
            .get(idx as usize)
            .copied()
            .unwrap_or(false)
    }

    pub fn is_terminal_rgroup_with_user_label(&self, idx: i32) -> bool {
        self.terminal_rgroup_atom_to_neighbor.contains_key(&idx)
    }

    /// `buildMatchingMol`
    fn build_matching_mol(&mut self) -> Result<()> {
        let mut mol = self.core.clone();
        self.terminal_rgroup_atom_to_neighbor.clear();
        mol.begin_batch_edit()?;
        for idx in 0..mol.num_atoms() {
            set_int_prop(mol.atom_mut(idx), RLABEL_CORE_INDEX, idx as i32);
            if mol.atom(idx).atomic_num() == 0
                && mol.atom_degree(idx) == 1
                && is_dummy_rgroup_attachment(&mol, idx)
            {
                if let Some(q) = mol.atom(idx).query() {
                    let description = q.description();
                    if description != "AtomNull" && description != "AtomAtomicNum" {
                        let ty = int_prop(mol.atom(idx), RLABEL_TYPE);
                        let is_isotope_label = description == "AtomIsotope"
                            && ty == Some(Labelling::IsotopeLabels as i32);
                        debug_assert_eq!(Labelling::IsotopeLabels as u32, ISOTOPE_LABELS);
                        if !is_isotope_label {
                            continue;
                        }
                    }
                }
                let neighbor = mol
                    .atom_neighbors(idx)
                    .next()
                    .ok_or_else(|| RDKitError::invariant("terminal dummy without neighbor"))?;
                self.terminal_rgroup_atom_to_neighbor
                    .entry(idx as i32)
                    .or_insert(neighbor as i32);
                mol.remove_atom(idx)?;
            }
        }
        mol.commit_batch_edit()?;
        self.matching_mol = mol;
        Ok(())
    }

    fn matching_index_to_core_index(&self, matching_index: i32) -> Result<i32> {
        int_prop(
            self.matching_mol.atom(matching_index as u32),
            RLABEL_CORE_INDEX,
        )
        .ok_or_else(|| RDKitError::invariant("Matched atom missing core index"))
    }

    /// `replaceCoreAtom`
    fn replace_core_atom(
        mol: &mut RWMol,
        atom_idx: u32,
        other: &ROMol,
        other_idx: u32,
    ) -> Result<()> {
        let other_atom = other.atom(other_idx);
        let atomic_number = other_atom.atomic_num();
        let was_dummy = mol.atom(atom_idx).atomic_num() == 0;
        if was_dummy || mol.atom(atom_idx).has_query() {
            if mol.atom(atom_idx).has_query() {
                let new_atom = Atom::new(atomic_number);
                mol.replace_atom(atom_idx, &new_atom, true)?;
            } else {
                mol.atom_mut(atom_idx).set_atomic_num(atomic_number as i32);
            }
        }
        {
            let target = mol.atom_mut(atom_idx);
            target.set_is_aromatic(other_atom.is_aromatic());
            target.set_formal_charge(other_atom.formal_charge());
        }
        if was_dummy {
            let num_hs = other
                .atom_neighbors(other_idx)
                .filter(|&n| other.atom(n).atomic_num() == 1)
                .count() as u32;
            let total_hs = other.total_num_hs(other_idx, false)?;
            let target = mol.atom_mut(atom_idx);
            target.set_no_implicit(true);
            target.set_num_explicit_hs(num_hs + total_hs);
            mol.update_atom_property_cache(atom_idx, false)?;
        }
        Ok(())
    }

    /// `replaceCoreAtomsWithMolMatches`: the core with its dummy and query
    /// atoms and bonds replaced by those of `mol` they matched.
    pub fn replace_core_atoms_with_mol_matches(
        &self,
        mol: &ROMol,
        match_: &MatchVect,
    ) -> Result<RWMol> {
        let mut replaced = self.core.clone();
        for &(core_idx, mol_idx) in match_ {
            if is_atom_with_multiple_neighbors_or_not_dummy_rgroup_attachment(
                &replaced,
                core_idx as u32,
            ) {
                Self::replace_core_atom(&mut replaced, core_idx as u32, mol, mol_idx as u32)?;
            }
        }
        let lookup: BTreeMap<i32, i32> = match_.iter().copied().collect();
        for bond_idx in 0..replaced.num_bonds() {
            if !replaced.bond(bond_idx).has_query() {
                continue;
            }
            let (begin, end) = {
                let b = replaced.bond(bond_idx);
                (b.begin_atom_idx(), b.end_atom_idx())
            };
            let mol_begin = lookup.get(&(begin as i32)).copied().unwrap_or(0);
            let mol_end = lookup.get(&(end as i32)).copied().unwrap_or(0);
            match mol.bond_idx_between(mol_begin as u32, mol_end as u32) {
                None => {
                    if replaced.atom(begin).atomic_num() != 0
                        && replaced.atom(end).atomic_num() != 0
                    {
                        return Err(RDKitError::invariant(
                            "Failed to find core bond in molecule",
                        ));
                    }
                }
                Some(mol_bond) => {
                    let mol_bond = mol.bond(mol_bond);
                    let mut new_bond = crate::graphmol::bond::Bond::new(mol_bond.bond_type());
                    new_bond.set_is_aromatic(mol_bond.is_aromatic());
                    replaced.replace_bond(bond_idx, &new_bond, true, false)?;
                }
            }
        }
        if mol.num_conformers() > 0 {
            if replaced.num_conformers() == 0 {
                let conf = Conformer::new(replaced.num_atoms() as usize);
                replaced.add_conformer(conf, true)?;
            }
            let positions: Vec<_> = {
                let mol_conf = mol.conformer(-1)?;
                match_
                    .iter()
                    .map(|&(c, m)| (c as u32, mol_conf.atom_pos(m as u32)))
                    .collect::<Vec<_>>()
            };
            let replaced_conf = replaced.conformer_mut(-1)?;
            for (core_idx, pos) in positions {
                replaced_conf.set_atom_pos(core_idx, pos?)?;
            }
        } else {
            replaced.clear_conformers();
        }
        Ok(replaced)
    }

    /// `extractCoreFromMolMatch`: the part of `mol` matched by the core, with
    /// dummy atoms for the attachment points of the R groups.
    pub fn extract_core_from_mol_match(
        &self,
        mol: &ROMol,
        match_: &MatchVect,
        params: &RGroupDecompositionParameters,
    ) -> Result<RWMol> {
        let mut extracted = mol.clone();
        let mut keep = vec![false; mol.num_atoms() as usize];
        let mut new_bonds = Vec::new();
        // new dummy atom -> the molecule atom it stands for
        let mut dummy_atom_map: Vec<(u32, u32)> = Vec::new();
        // molecule atom -> its atom in `extracted`
        let mut mol_atom_map: BTreeMap<u32, u32> = BTreeMap::new();

        for &(core_idx, target_idx) in match_ {
            let (core_idx, target_idx) = (core_idx as u32, target_idx as u32);
            let query_atom = self.core.atom(core_idx);
            if let Some(label) = int_prop(query_atom, RLABEL) {
                set_int_prop(extracted.atom_mut(target_idx), RLABEL, label);
            }
            if let Some(ty) = int_prop(query_atom, RLABEL_TYPE) {
                set_int_prop(extracted.atom_mut(target_idx), RLABEL_TYPE, ty);
            }
            if query_atom.atomic_num() == 0
                && query_atom.has_prop(RLABEL)
                && self.core.atom_degree(core_idx) == 1
            {
                continue;
            }
            keep[target_idx as usize] = true;
            mol_atom_map.insert(target_idx, target_idx);
            let mut neighbor_number: i32 = -1;
            let mut is_chiral = matches!(
                extracted.atom(target_idx).chiral_tag(),
                ChiralType::TetrahedralCCW | ChiralType::TetrahedralCW
            );
            if is_chiral && !params.substructmatch_params.use_chirality {
                extracted
                    .atom_mut(target_idx)
                    .set_chiral_tag(ChiralType::Unspecified);
                is_chiral = false;
            }
            let target_neighbors: Vec<u32> = extracted.atom_neighbors(target_idx).collect();
            let mut hydrogens_to_add = Vec::new();
            for target_neighbor in target_neighbors {
                neighbor_number += 1;
                let query_neighbor = match_.iter().find(|&&(q, m)| {
                    m as u32 == target_neighbor
                        && self.core.bond_idx_between(core_idx, q as u32).is_some()
                });
                let Some(&(query_neighbor_idx, _)) = query_neighbor else {
                    if extracted.atom(target_neighbor).atomic_num() == 1 {
                        hydrogens_to_add.push(target_neighbor);
                        mol_atom_map.insert(target_neighbor, target_neighbor);
                    } else if is_chiral {
                        is_chiral = false;
                        extracted
                            .atom_mut(target_idx)
                            .set_chiral_tag(ChiralType::Unspecified);
                    }
                    continue;
                };
                let query_neighbor_atom = self.core.atom(query_neighbor_idx as u32);
                if !(query_neighbor_atom.atomic_num() == 0
                    && query_neighbor_atom.has_prop(RLABEL)
                    && self.core.atom_degree(query_neighbor_idx as u32) == 1)
                {
                    continue;
                }
                // upstream copies the `Atom` part only: the query is not kept
                let mut new_dummy = query_neighbor_atom.copy_detached();
                new_dummy.set_query(None);
                new_dummy.clear_computed_props();
                let new_dummy_idx = extracted.add_atom(new_dummy, false);
                dummy_atom_map.push((new_dummy_idx, target_neighbor));
                if new_dummy_idx as usize >= keep.len() {
                    keep.resize(new_dummy_idx as usize + 1, false);
                }
                keep[new_dummy_idx as usize] = true;
                let connecting_idx = extracted
                    .bond_idx_between(target_idx, target_neighbor)
                    .ok_or_else(|| RDKitError::invariant("connecting bond not found"))?;
                let mut connecting = extracted.bond(connecting_idx).copy_detached();
                if connecting.stereo() > BondStereo::StereoAny {
                    connecting.set_stereo(BondStereo::StereoAny)?;
                }
                if connecting.begin_atom_idx() == target_neighbor {
                    connecting.set_begin_atom_idx_unchecked(new_dummy_idx);
                } else {
                    connecting.set_end_atom_idx_unchecked(new_dummy_idx);
                }
                new_bonds.push(connecting);

                let other_bonds: Vec<u32> = extracted
                    .atom_bond_indices(target_idx)
                    .filter(|&b| b != connecting_idx)
                    .collect();
                for b in other_bonds {
                    let bond = extracted.bond_mut(b);
                    if bond.stereo() > BondStereo::StereoAny {
                        for stereo_atom in bond.stereo_atoms_mut() {
                            if *stereo_atom == target_neighbor as i32 {
                                *stereo_atom = new_dummy_idx as i32;
                            }
                        }
                    }
                }
                if is_chiral {
                    let switch_it = match extracted.atom_degree(target_idx) {
                        4 => neighbor_number % 2 == 0,
                        3 => neighbor_number == 1,
                        _ => false,
                    };
                    neighbor_number -= 1;
                    if switch_it {
                        extracted.atom_mut(target_idx).invert_chirality();
                    }
                }
            }
            for h in hydrogens_to_add {
                keep[h as usize] = true;
            }
        }
        for bond in new_bonds {
            extracted.add_bond_object(bond)?;
        }

        extracted.begin_batch_edit()?;
        let mut removed = vec![false; extracted.num_atoms() as usize];
        for idx in 0..extracted.num_atoms() {
            if !keep.get(idx as usize).copied().unwrap_or(false) {
                extracted.remove_atom(idx)?;
                removed[idx as usize] = true;
            }
        }
        extracted.commit_batch_edit()?;
        // index of an atom of the unreduced molecule in `extracted`
        let mut new_index = vec![u32::MAX; removed.len()];
        let mut next = 0;
        for (old, &gone) in removed.iter().enumerate() {
            if !gone {
                new_index[old] = next;
                next += 1;
            }
        }

        update_sub_mol_confs(mol, &mut extracted, &removed)?;
        for conf in mol.conformers() {
            let id = conf.id() as i32;
            for &(dummy, target) in &dummy_atom_map {
                let pos = conf.atom_pos(target)?;
                extracted
                    .conformer_mut(id)?
                    .set_atom_pos(new_index[dummy as usize], pos)?;
            }
        }

        if mol.num_conformers() == 0 && self.core.num_conformers() > 0 {
            let mut mol_copy = mol.clone();
            for conf in self.core.conformers() {
                let mut new_conf = Conformer::new(mol.num_atoms() as usize);
                new_conf.set_id(conf.id());
                new_conf.set_3d(conf.is_3d());
                for &(core_idx, target_idx) in match_ {
                    new_conf.set_atom_pos(target_idx as u32, conf.atom_pos(core_idx as u32)?)?;
                }
                mol_copy.add_conformer(new_conf, false)?;
            }
            update_sub_mol_confs(&mol_copy, &mut extracted, &removed)?;
            mol_copy.clear_conformers();
            for idx in 0..extracted.num_atoms() {
                if !is_user_rlabel(extracted.atom(idx)) {
                    continue;
                }
                let Some(label) = int_prop(extracted.atom(idx), RLABEL) else {
                    continue;
                };
                let core_atom = self
                    .core
                    .atoms()
                    .iter()
                    .position(|a| int_prop(a, RLABEL) == Some(label));
                if let Some(core_atom) = core_atom {
                    for (i, conf) in self.core.conformers().iter().enumerate() {
                        let pos = conf.atom_pos(core_atom as u32)?;
                        extracted.conformer_mut(i as i32)?.set_atom_pos(idx, pos)?;
                    }
                }
            }
        }

        let atom_map: BTreeMap<u32, u32> = mol_atom_map
            .iter()
            .map(|(&old, &same)| (old, new_index[same as usize]))
            .filter(|&(_, new)| new != u32::MAX)
            .collect();
        copy_stereo_groups(&atom_map, mol, &mut extracted)?;
        extracted.clear_computed_props(true);
        extracted.update_property_cache(false)?;
        // failing to sanitize is fine, the core only needs its ring info
        let _ = sanitize_mol(&mut extracted, SANITIZE_SYMMRINGS | SANITIZE_CLEANUP);
        Ok(extracted)
    }

    /// `checkAllBondsToRGroupPresent`
    pub fn check_all_bonds_to_rgroup_present(
        &self,
        mol: &ROMol,
        target_rgroup_idx: i32,
        target_to_core_indices: &[Vec<i32>],
    ) -> Result<bool> {
        let mut core_neighbors = BTreeSet::new();
        for nbr in mol.atom_neighbors(target_rgroup_idx as u32) {
            if mol.atom(nbr).atomic_num() < 1 {
                continue;
            }
            let core_atoms = &target_to_core_indices[nbr as usize];
            if core_atoms.is_empty() {
                continue;
            }
            if core_atoms.len() != 1 {
                return Err(RDKitError::invariant(
                    "Target atom neighboring R group should have exactly one match in core",
                ));
            }
            let core_idx = core_atoms[0];
            let core_atom = self.core.atom(core_idx as u32);
            if !(core_atom.atomic_num() == 0 && is_user_rlabel(core_atom)) {
                core_neighbors.insert(core_idx);
            }
        }
        if core_neighbors.is_empty() {
            return Err(RDKitError::invariant(
                "Unable to find target atom(s) matching core for attachment point",
            ));
        }
        if core_neighbors.len() == 1 {
            return Ok(true);
        }
        let core_rgroup_indices = &target_to_core_indices[target_rgroup_idx as usize];
        if core_rgroup_indices.len() != core_neighbors.len() {
            return Ok(false);
        }
        let mut bond_indices = BTreeSet::new();
        for &neighbor in &core_neighbors {
            for &rgroup in core_rgroup_indices {
                if let Some(bond) = self.core.bond_idx_between(neighbor as u32, rgroup as u32) {
                    bond_indices.insert(bond);
                }
            }
        }
        Ok(bond_indices.len() == core_neighbors.len())
    }

    /// `matchTerminalUserRGroups`: extends a match of the matching molecule
    /// (core without the terminal R-label dummies) by the ways the dummies
    /// can be matched to neighbors in `target`.
    pub fn match_terminal_user_rgroups(
        &self,
        target: &RWMol,
        match_: &MatchVect,
        sss_params: &SubstructMatchParameters,
    ) -> Result<Vec<MatchVect>> {
        let match_: MatchVect = match_
            .iter()
            .map(|&(q, t)| Ok((self.matching_index_to_core_index(q)?, t)))
            .collect::<Result<_>>()?;
        let mut all_mappings = Vec::new();
        if self.terminal_rgroup_atom_to_neighbor.is_empty() {
            all_mappings.push(match_);
            return Ok(all_mappings);
        }
        let match_map: BTreeMap<i32, i32> = match_.iter().copied().collect();
        let mut mapped_target = vec![false; target.num_atoms() as usize];
        for &(_, t) in &match_ {
            mapped_target[t as usize] = true;
        }

        let mut missing_dummies: Vec<i32> = Vec::new();
        let mut available_mappings_for_dummy: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
        let mut symmetric_hydrogens = vec![false; target.num_atoms() as usize];
        let mut target_atom_bonds_to_core_counts: BTreeMap<i32, i32> = BTreeMap::new();

        for core_idx in 0..self.core.num_atoms() {
            let core_atom = self.core.atom(core_idx);
            if core_atom.atomic_num() == 1
                || self.is_terminal_rgroup_with_user_label(core_idx as i32)
            {
                continue;
            }
            let mut dummy_indexes: Vec<i32> = self
                .core
                .atom_neighbors(core_idx)
                .filter(|&n| self.is_terminal_rgroup_with_user_label(n as i32))
                .map(|n| n as i32)
                .collect();
            if dummy_indexes.is_empty() {
                continue;
            }
            let label_of = |idx: i32| int_prop(self.core.atom(idx as u32), RLABEL).unwrap_or(0);
            // user labels (positive) ascending, then the other labels descending
            dummy_indexes.sort_by(|&a, &b| {
                let (l1, l2) = (label_of(a), label_of(b));
                let less = |x: i32, y: i32| if x > 0 && y > 0 { x < y } else { x > y };
                if less(l1, l2) {
                    std::cmp::Ordering::Less
                } else if less(l2, l1) {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Equal
                }
            });

            let mut neighbor_dummy_lists: Vec<Vec<i32>> = Vec::new();
            for &dummy_idx in &dummy_indexes {
                let neighbor_idx = self.terminal_rgroup_atom_to_neighbor[&dummy_idx];
                let core_bond = self
                    .core
                    .bond_idx_between(dummy_idx as u32, neighbor_idx as u32)
                    .ok_or_else(|| RDKitError::invariant("core bond not found"))?;
                let target_idx = match_map.get(&neighbor_idx).copied().unwrap_or(0);
                let mut available = Vec::new();
                for nbr in target.atom_neighbors(target_idx as u32) {
                    if mapped_target[nbr as usize] {
                        continue;
                    }
                    let target_bond = target
                        .bond_idx_between(target_idx as u32, nbr)
                        .ok_or_else(|| RDKitError::invariant("target bond not found"))?;
                    if bond_compat(&self.core, core_bond, target, target_bond, sss_params)? {
                        available.push(nbr as i32);
                        *target_atom_bonds_to_core_counts
                            .entry(nbr as i32)
                            .or_insert(0) += 1;
                    }
                }
                if available.len() > 1
                    && available
                        .iter()
                        .all(|&i| target.atom(i as u32).atomic_num() == 1)
                {
                    let single = available
                        .iter()
                        .copied()
                        .find(|&i| !symmetric_hydrogens[i as usize])
                        .unwrap_or(available[0]);
                    available = vec![single];
                    symmetric_hydrogens[single as usize] = true;
                }
                neighbor_dummy_lists.push(available);
            }

            // every dummy may also stay unmapped, marked by a negative number
            let mut start = 0;
            let with_unmapped: Vec<Vec<i32>> = neighbor_dummy_lists
                .iter()
                .map(|v| {
                    let mut v = v.clone();
                    start -= 1;
                    v.push(start);
                    v
                })
                .collect();
            let best = best_dummy_mapping(cartesian_product(&with_unmapped, false));
            for (i, &dummy_idx) in dummy_indexes.iter().enumerate() {
                if best[i] < 0 {
                    let dummy = self.core.atom(dummy_idx as u32);
                    if dummy.has_prop(UNLABELED_CORE_ATTACHMENT)
                        || (is_user_rlabel(dummy)
                            && (core_atom.atomic_num() == 0 || core_atom.has_query()))
                    {
                        missing_dummies.push(dummy_idx);
                    } else {
                        return Ok(all_mappings);
                    }
                } else {
                    available_mappings_for_dummy.insert(dummy_idx, neighbor_dummy_lists[i].clone());
                }
            }
        }
        if available_mappings_for_dummy.is_empty() {
            all_mappings.push(match_);
            return Ok(all_mappings);
        }

        let dummies_with_mapping: Vec<i32> = available_mappings_for_dummy.keys().copied().collect();
        let available_lists: Vec<Vec<i32>> =
            available_mappings_for_dummy.values().cloned().collect();
        let allow_duplicates = target_atom_bonds_to_core_counts.values().any(|&c| c > 1);
        let all_available = cartesian_product(&available_lists, allow_duplicates);
        if all_available.is_empty() {
            all_mappings.push(match_);
            return Ok(all_mappings);
        }
        let size = all_available[0].len() + match_.len();

        // when dummies are missing the core is checked without them
        let has_missing = !missing_dummies.is_empty();
        let mut check_core: Option<RWMol> = None;
        let mut core_to_check: BTreeMap<i32, u32> = BTreeMap::new();
        if has_missing {
            let mut check = self.core.clone();
            for idx in 0..check.num_atoms() {
                set_int_prop(check.atom_mut(idx), "__core_index__", idx as i32);
            }
            missing_dummies.sort_by(|a, b| b.cmp(a));
            for &index in &missing_dummies {
                let neighbor = check
                    .atom_neighbors(index as u32)
                    .next()
                    .ok_or_else(|| RDKitError::invariant("dummy without neighbor"))?;
                check.remove_atom(index as u32)?;
                // removal shifted the indices above `index`
                let neighbor = if neighbor > index as u32 {
                    neighbor - 1
                } else {
                    neighbor
                };
                check.update_atom_property_cache(neighbor, false)?;
            }
            for idx in 0..check.num_atoms() {
                let core_index = int_prop(check.atom(idx), "__core_index__")
                    .ok_or_else(|| RDKitError::invariant("missing core index"))?;
                core_to_check.insert(core_index, idx);
            }
            check_core = Some(check);
        }

        let mut query_indices = vec![0u32; size];
        let mut target_indices = vec![0u32; size];
        for (position, &(q, t)) in match_.iter().enumerate() {
            query_indices[position] = if has_missing {
                core_to_check.get(&q).copied().unwrap_or(0)
            } else {
                q as u32
            };
            target_indices[position] = t as u32;
        }
        let query_matching_mol: &RWMol = check_core.as_ref().unwrap_or(&self.core);
        let mut final_check = FinalCheck::new(query_matching_mol, target, sss_params);
        let mut target_bonds_present = vec![false; target.num_bonds() as usize];

        for dummy_mapping in &all_available {
            if match_.len() + dummy_mapping.len() != size {
                return Err(RDKitError::invariant("Size error in dummy mapping"));
            }
            let mut duplicate_bonds = false;
            target_bonds_present.iter_mut().for_each(|b| *b = false);
            for (i, &mapped_target_idx) in dummy_mapping.iter().enumerate() {
                let position = match_.len() + i;
                let dummy = dummies_with_mapping[i];
                query_indices[position] = if has_missing {
                    core_to_check.get(&dummy).copied().unwrap_or(0)
                } else {
                    dummy as u32
                };
                target_indices[position] = mapped_target_idx as u32;
                if allow_duplicates {
                    let neighbor_idx = self.terminal_rgroup_atom_to_neighbor[&dummy];
                    let target_neighbor = match_map.get(&neighbor_idx).copied().unwrap_or(0);
                    let bond = target
                        .bond_idx_between(mapped_target_idx as u32, target_neighbor as u32)
                        .ok_or_else(|| RDKitError::invariant("Matching target bond not found"))?;
                    if target_bonds_present[bond as usize] {
                        duplicate_bonds = true;
                        break;
                    }
                    target_bonds_present[bond as usize] = true;
                }
            }
            if !duplicate_bonds && final_check.check(&query_indices, &target_indices)? {
                let mut with_dummy = match_.clone();
                for (i, &t) in dummy_mapping.iter().enumerate() {
                    with_dummy.push((dummies_with_mapping[i], t));
                }
                all_mappings.push(with_dummy);
            }
        }
        Ok(all_mappings)
    }
}

/// The preferred assignment of the dummies of a core atom: the one leaving
/// the fewest unmapped (negative entries); among assignments leaving equally
/// many unmapped, one that maps earlier dummies is preferred. The first
/// candidate wins ties.
fn best_dummy_mapping(candidates: Vec<Vec<i32>>) -> Vec<i32> {
    let unmapped = |c: &[i32]| c.iter().filter(|&&v| v < 0).count();
    let better = |a: &[i32], b: &[i32]| -> bool {
        let (ua, ub) = (unmapped(a), unmapped(b));
        if ua != ub {
            return ua < ub;
        }
        if ua == 0 {
            return false;
        }
        for (&v1, &v2) in a.iter().zip(b) {
            if (v1 > 0 && v2 > 0) || (v1 < 0 && v2 < 0) {
                continue;
            }
            if v1 > 0 && v2 < 0 {
                return true;
            }
            if v1 < 0 && v2 > 0 {
                break;
            }
        }
        false
    };
    let mut iter = candidates.into_iter();
    let mut best = iter.next().unwrap_or_default();
    for c in iter {
        if better(&c, &best) {
            best = c;
        }
    }
    best
}
