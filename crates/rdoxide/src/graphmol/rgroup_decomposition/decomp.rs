//! The R-group decomposition driver (upstream
//! `RGroupDecomposition/RGroupDecomp.{h,cpp}`).

use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, RwLock};
use std::time::Instant;

use super::data::{RData, RGroupData, RGroupMatch};
use super::decomp_data::{RGroupDecompData, RGroupDecompositionProcessResult};
use super::params::{RGroupDecompositionParameters, rgroup_matching};
use super::score::check_for_timeout;
use super::utils::{RGROUP_INPUT_DUMMY, RLABEL, SIDECHAIN_RLABELS, int_prop};
use crate::graphmol::chem_transforms::replace_core_with_match;
use crate::graphmol::mol_ops::add_hs::{
    AddHsParameters, RemoveHsParameters, add_hs, remove_hs, set_terminal_atom_coords,
};
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::log::LogStateSetter;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_debug};

const TARGET_ATOM_IDX: &str = "__rgdTargetAtomIdx";
const TARGET_BOND_IDX: &str = "__rgdTargetBondIdx";
const MAX_PERMUTATIONS: usize = 100_000;

/// One decomposed molecule: R group label -> molecule.
pub type RGroupRow = BTreeMap<String, ROMol>;
pub type RGroupRows = Vec<RGroupRow>;
/// One R group over all the molecules.
pub type RGroupColumn = Vec<ROMol>;
pub type RGroupColumns = BTreeMap<String, RGroupColumn>;

/// `UsedLabelMap`: for every final label whether the molecule being output
/// used it and whether the user defined it.
#[derive(Debug, Clone)]
struct UsedLabelMap {
    map: BTreeMap<i32, (bool, bool)>,
}

impl UsedLabelMap {
    fn new(mapping: &BTreeMap<i32, i32>) -> Self {
        UsedLabelMap {
            map: mapping.iter().map(|(&k, &v)| (v, (false, k > 0))).collect(),
        }
    }

    fn has(&self, label: i32) -> bool {
        self.map.contains_key(&label)
    }

    fn is_used(&self, label: i32) -> bool {
        self.map[&label].0
    }

    fn set_used(&mut self, label: i32) {
        self.map.entry(label).or_default().0 = true;
    }

    fn is_user_defined(&self, label: i32) -> bool {
        self.map[&label].1
    }
}

/// `RGroupDecomposition`
#[derive(Debug)]
pub struct RGroupDecomposition {
    data: RGroupDecompData,
}

impl RGroupDecomposition {
    pub fn new(core: &ROMol, params: RGroupDecompositionParameters) -> Result<Self> {
        Self::new_multi(&[core], params)
    }

    pub fn new_multi(cores: &[&ROMol], params: RGroupDecompositionParameters) -> Result<Self> {
        Ok(RGroupDecomposition {
            data: RGroupDecompData::new(cores, params)?,
        })
    }

    pub fn params(&self) -> &RGroupDecompositionParameters {
        &self.data.params
    }

    fn label_atom_bond_indices(mol: &mut RWMol) {
        for idx in 0..mol.num_atoms() {
            mol.atom_mut(idx)
                .dict_mut()
                .set_val(TARGET_ATOM_IDX, idx as i32);
        }
        for idx in 0..mol.num_bonds() {
            mol.bond_mut(idx)
                .dict_mut()
                .set_val(TARGET_BOND_IDX, idx as i32);
        }
    }

    /// `setTargetAtomBondIndices`: records on `mol` which atoms and bonds of
    /// the target it came from, for highlighting.
    fn set_target_atom_bond_indices(&self, mol: &mut ROMol, include_bonds_to_rlabels: bool) {
        let remove_hs_post_match = self.data.params.remove_hydrogens_post_match;
        let mut atom_indices = vec![-1i32; mol.num_atoms() as usize];
        let mut bond_indices = vec![-1i32; mol.num_bonds() as usize];
        let is_hydrogen = RGroupData::is_mol_hydrogen(mol);
        for idx in 0..mol.num_atoms() {
            let atom = mol.atom_mut(idx);
            let Some(target_idx) = int_prop(atom, TARGET_ATOM_IDX) else {
                continue;
            };
            atom.clear_prop(TARGET_ATOM_IDX);
            if (atom.atomic_num() == 1 && remove_hs_post_match && !is_hydrogen)
                || (atom.atomic_num() == 0 && !atom.has_prop(RGROUP_INPUT_DUMMY))
            {
                continue;
            }
            atom_indices[idx as usize] = target_idx;
        }
        for idx in 0..mol.num_bonds() {
            let (begin, end) = {
                let b = mol.bond(idx);
                (b.begin_atom_idx(), b.end_atom_idx())
            };
            let Some(RDValue::Int(target_idx)) = mol.bond(idx).get_prop(TARGET_BOND_IDX).cloned()
            else {
                continue;
            };
            mol.bond_mut(idx).clear_prop(TARGET_BOND_IDX);
            let (begin_atom, end_atom) = (mol.atom(begin), mol.atom(end));
            if (begin_atom.atomic_num() == 1 || end_atom.atomic_num() == 1)
                && remove_hs_post_match
                && !is_hydrogen
            {
                continue;
            }
            if !include_bonds_to_rlabels
                && ((begin_atom.atomic_num() == 0 && !begin_atom.has_prop(RGROUP_INPUT_DUMMY))
                    || (end_atom.atomic_num() == 0 && !end_atom.has_prop(RGROUP_INPUT_DUMMY)))
            {
                continue;
            }
            bond_indices[idx as usize] = target_idx;
        }
        // trailing unset entries are dropped, then the unset ones
        let highlight_atoms: Vec<i32> = trim_unset(atom_indices);
        let highlight_bonds: Vec<i32> = trim_unset(bond_indices);
        mol.dict_mut()
            .set_val(cp::_rgroupTargetAtoms, highlight_atoms);
        mol.dict_mut()
            .set_val(cp::_rgroupTargetBonds, highlight_bonds);
    }

    /// `getMatchingCoreIdx`: the index of the core matching `mol` (-1 if
    /// none) and the unique matches of the heavy atoms of `mol`.
    pub fn matching_core_idx(&self, mol: &ROMol) -> Result<(i32, Vec<MatchVect>)> {
        let mut rwmol = mol.clone();
        let (core_idx, matches) = self.matching_core_internal(&mut rwmol)?;
        let num_atoms = mol.num_atoms() as i32;
        let mut unique = BTreeSet::new();
        for m in matches {
            let mut heavy: MatchVect = m.into_iter().filter(|&(_, t)| t < num_atoms).collect();
            heavy.sort();
            unique.insert(heavy);
        }
        Ok((core_idx, unique.into_iter().collect()))
    }

    /// `getMatchingCoreInternal`: adds hydrogens to `mol`, finds the core
    /// matching it best and the matches on that core.
    fn matching_core_internal(&self, mol: &mut RWMol) -> Result<(i32, Vec<MatchVect>)> {
        add_hs(
            mol,
            &AddHsParameters {
                explicit_only: false,
                add_coords: true,
                ..Default::default()
            },
            None,
        )?;
        let params = &self.data.params;
        let mut sss_params = params.substructmatch_params.clone();
        sss_params.uniquify = false;
        sss_params.recursion_possible = true;

        let mut core_idx = -1;
        let mut best_core = None;
        let mut global_min_heavy_nbrs = -1i32;
        let mut filtered: Vec<MatchVect> = Vec::new();
        for (&key, core) in &self.data.cores {
            if params.do_tautomers {
                return Err(RDKitError::Runtime(
                    "doTautomers is not supported: TautomerQuery is not ported".into(),
                ));
            }
            let base_matches = substruct_match(mol, &core.matching_mol, &sss_params)?;
            let mut tmatches: Vec<MatchVect> = Vec::new();
            for base in &base_matches {
                tmatches.extend(core.match_terminal_user_rgroups(mol, base, &sss_params)?);
            }
            if tmatches.is_empty() {
                continue;
            }
            let mut heavy_nbrs = vec![0i32; tmatches.len()];
            for (i, mv) in tmatches.iter().enumerate() {
                let mut passes_filter = params.only_match_at_rgroups;
                let mut target_to_core: Vec<Vec<i32>> = vec![Vec::new(); mol.num_atoms() as usize];
                for &(q, t) in mv {
                    target_to_core[t as usize].push(q);
                }
                for &(q, t) in mv {
                    if !core.is_core_atom_user_labelled(q) {
                        for nbr in mol.atom_neighbors(t as u32) {
                            if mol.atom(nbr).atomic_num() != 1
                                && target_to_core[nbr as usize].is_empty()
                            {
                                if params.only_match_at_rgroups {
                                    passes_filter = false;
                                    break;
                                }
                                heavy_nbrs[i] += 1;
                            }
                        }
                    } else if core.is_terminal_rgroup_with_user_label(q)
                        && params.only_match_at_rgroups
                        && !core.check_all_bonds_to_rgroup_present(mol, t, &target_to_core)?
                    {
                        passes_filter = false;
                    }
                    if !passes_filter && params.only_match_at_rgroups {
                        break;
                    }
                }
                if passes_filter {
                    filtered.push(mv.clone());
                }
            }
            if !params.only_match_at_rgroups {
                let min_heavy_nbrs = *heavy_nbrs.iter().min().expect("non-empty matches");
                let better = match best_core {
                    None => true,
                    Some(best_key) => {
                        let best: &super::core::RCore = &self.data.cores[&best_key];
                        min_heavy_nbrs < global_min_heavy_nbrs
                            && !substruct_match(&core.core, &best.core, &sss_params)?.is_empty()
                    }
                };
                if better {
                    filtered = tmatches
                        .into_iter()
                        .zip(&heavy_nbrs)
                        .filter(|&(_, &h)| h <= min_heavy_nbrs)
                        .map(|(m, _)| m)
                        .collect();
                    global_min_heavy_nbrs = min_heavy_nbrs;
                    best_core = Some(key);
                    core_idx = key;
                    if global_min_heavy_nbrs == 0 {
                        break;
                    }
                }
            } else if !filtered.is_empty() {
                best_core = Some(key);
                core_idx = key;
                break;
            }
        }
        if best_core.is_none() {
            filtered.clear();
        }
        Ok((core_idx, filtered))
    }

    /// `add`: adds a molecule to the decomposition. Returns its index, -1 if
    /// no core matches, or -2 if it has no side chains.
    pub fn add(&mut self, inmol: &ROMol) -> Result<i32> {
        let mut mol = inmol.clone();
        let (core_idx, mut tmatches) = self.matching_core_internal(&mut mol)?;
        if core_idx == -1 {
            rd_log_debug!("No core matches");
            return Ok(-1);
        }
        if self.data.params.include_target_mol_in_results {
            Self::label_atom_bond_indices(&mut mol);
        }
        if tmatches.len() > 1 {
            let strategy = self.data.params.matching_strategy;
            if strategy == rgroup_matching::NO_SYMMETRIZATION
                || (self.data.matches.is_empty() && strategy == rgroup_matching::GREEDY)
            {
                tmatches.truncate(1);
            }
        }

        for idx in 0..mol.num_atoms() {
            let atom = mol.atom_mut(idx);
            if atom.atomic_num() == 0 {
                atom.dict_mut().set_val(RGROUP_INPUT_DUMMY, true);
                atom.set_isotope(0);
                atom.set_atom_map_num(0, false)?;
                atom.clear_prop(cp::_MolFileRLabel);
                atom.dict_mut().set_val(cp::dummyLabel, "*");
            }
        }

        let non_unique_matches = if self.data.params.allow_multiple_cores_in_same_mol {
            split_non_unique_matches(&tmatches, mol.num_atoms())
        } else {
            vec![tmatches]
        };

        let mut potential_matches: Vec<RGroupMatch> = Vec::new();
        for split_match in &non_unique_matches {
            for tmatch in split_match {
                if let Some(index) =
                    self.add_match(inmol, &mol, core_idx, tmatch, &mut potential_matches)?
                {
                    return Ok(index);
                }
            }
            if potential_matches.is_empty() {
                rd_log_debug!("No attachment points in side chains");
                return Ok(-2);
            }
            if self.data.params.matching_strategy != rgroup_matching::GA {
                let n = self.data.matches[self.data.previous_match_size..]
                    .iter()
                    .fold(1usize, |n, m| n.saturating_mul(m.len()));
                if n.saturating_mul(potential_matches.len()) > MAX_PERMUTATIONS {
                    self.data.process(self.data.prune_permutations, false)?;
                }
            }
            self.data
                .matches
                .push(std::mem::take(&mut potential_matches));
        }
        if !self.data.matches.is_empty() {
            let params = &self.data.params;
            let strategy = params.matching_strategy;
            if strategy & rgroup_matching::GREEDY != 0
                || (strategy & rgroup_matching::GREEDY_CHUNKS != 0
                    && params.chunk_size != 0
                    && self.data.matches.len() % params.chunk_size as usize == 0)
            {
                self.data.process(self.data.prune_permutations, false)?;
            }
        }
        Ok(self.data.matches.len() as i32 - 1)
    }

    /// The body of the loop over the matches in `add`: decomposes `mol` for
    /// one match of the core, appending the result to `potential_matches`.
    /// `Some(index)` when a new core was found and the molecule was added
    /// again.
    fn add_match(
        &mut self,
        inmol: &ROMol,
        mol: &RWMol,
        core_idx: i32,
        tmatch: &MatchVect,
        potential_matches: &mut Vec<RGroupMatch>,
    ) -> Result<Option<i32>> {
        let rcore = &self.data.cores[&core_idx];
        let core_copy = rcore.replace_core_atoms_with_mol_matches(mol, tmatch)?;
        let Some(t_mol) = replace_core_with_match(mol, &core_copy, tmatch, false, true, false)?
        else {
            return Ok(None);
        };
        let mut matched: BTreeMap<i32, RGroupData> = BTreeMap::new();
        let (fragments, _) = get_mol_frags_as_mols(&t_mol, false, true)?;
        let num_fragments = fragments.len();
        let mut new_labels = Vec::new();
        for (i, mut new_mol) in fragments.into_iter().enumerate() {
            let mut rlabels_on_side_chain: Vec<i32> = Vec::new();
            new_mol.dict_mut().set_val("core", core_idx);
            new_mol
                .dict_mut()
                .set_val("idx", self.data.matches.len() as i32);
            new_mol.dict_mut().set_val("frag_idx", i as i32);
            for atom_idx in 0..new_mol.num_atoms() {
                if new_mol.atom(atom_idx).atomic_num() != 0 {
                    continue;
                }
                if new_mol.atom(atom_idx).has_prop(RGROUP_INPUT_DUMMY) {
                    new_mol.atom_mut(atom_idx).clear_prop(RGROUP_INPUT_DUMMY);
                    continue;
                }
                let core_atom_index = new_mol.atom(atom_idx).isotope();
                let core_atom = rcore.core.atom(core_atom_index);
                let Some(rlabel) = int_prop(core_atom, RLABEL) else {
                    continue;
                };
                let atom = new_mol.atom_mut(atom_idx);
                let mut on_atom: Vec<i32> = match atom.get_prop(SIDECHAIN_RLABELS) {
                    Some(RDValue::VecInt(v)) => v.clone(),
                    _ => Vec::new(),
                };
                on_atom.push(rlabel);
                atom.dict_mut().set_val(SIDECHAIN_RLABELS, on_atom);
                new_labels.push(rlabel); // keep track of all labels used
                rlabels_on_side_chain.push(rlabel);
                let first_bond = new_mol.atom_bond_indices(atom_idx).next();
                if let Some(bond_idx) = first_bond {
                    let connecting = new_mol.bond_mut(bond_idx);
                    if connecting.stereo() > crate::graphmol::bond::BondStereo::StereoAny {
                        connecting.set_stereo(crate::graphmol::bond::BondStereo::StereoAny)?;
                    }
                }
            }
            if self.data.params.include_target_mol_in_results {
                self.set_target_atom_bond_indices(&mut new_mol, true);
            }
            if !rlabels_on_side_chain.is_empty() {
                let new_mol = Arc::new(new_mol);
                for &rlabel in &rlabels_on_side_chain {
                    matched
                        .entry(rlabel)
                        .or_default()
                        .add(&new_mol, &rlabels_on_side_chain)?;
                }
            } else if num_fragments == 1 {
                // need to make a new core
                let mut new_core = mol.clone();
                for &(core_atom, mol_atom) in tmatch {
                    let rlabel = int_prop(rcore.core.atom(core_atom as u32), RLABEL);
                    let atom = new_core.atom_mut(mol_atom as u32);
                    if let Some(rlabel) = rlabel {
                        atom.dict_mut().set_val(RLABEL, rlabel);
                    }
                    atom.dict_mut().set_val("keep", true);
                }
                new_core.begin_batch_edit()?;
                for idx in 0..new_core.num_atoms() {
                    if !new_core.atom(idx).has_prop("keep") {
                        new_core.remove_atom(idx)?;
                    }
                }
                new_core.commit_batch_edit()?;
                if new_core.num_atoms() > 0 {
                    let smiles =
                        mol_to_smiles(&mut new_core.clone(), &SmilesWriteParams::default())?;
                    if !self.data.new_cores.contains_key(&smiles) {
                        let key = self.data.new_core_label;
                        self.data.new_core_label -= 1;
                        self.data.new_cores.insert(smiles, key);
                        self.data
                            .cores
                            .insert(key, super::core::RCore::new(new_core)?);
                        return self.add(inmol).map(Some);
                    }
                }
            }
        }
        self.data.labels.extend(new_labels);

        if !matched.is_empty() {
            let rgroups: BTreeMap<i32, RData> = matched
                .into_iter()
                .map(|(label, data)| (label, Arc::new(RwLock::new(data))))
                .collect();
            let rcore = &self.data.cores[&core_idx];
            let num_user_groups = rgroups
                .iter()
                .filter(|(label, rg)| **label > 0 && !rg.read().expect("rgroup lock").is_hydrogen)
                .count();
            let number_missing = rcore
                .number_user_rgroups
                .checked_sub(num_user_groups)
                .ok_or_else(|| RDKitError::invariant("Data error in missing user rgroup count"))?;
            let mut extracted =
                rcore.extract_core_from_mol_match(mol, tmatch, &self.data.params)?;
            if self.data.params.include_target_mol_in_results {
                self.set_target_atom_bond_indices(&mut extracted, false);
            }
            let mut potential = RGroupMatch::new(core_idx, number_missing, rgroups, extracted);
            if self.data.params.include_target_mol_in_results {
                potential.set_target_molecule_for_highlights(mol.clone());
            }
            potential_matches.push(potential);
        }
        Ok(None)
    }

    /// `process`
    pub fn process(&mut self) -> bool {
        self.process_and_score().success
    }

    /// `processAndScore`: failures are reported as an unsuccessful result.
    pub fn process_and_score(&mut self) -> RGroupDecompositionProcessResult {
        self.data
            .process(self.data.prune_permutations, true)
            .unwrap_or(RGroupDecompositionProcessResult {
                success: false,
                score: -1.0,
            })
    }

    /// `getRGroupLabels`
    pub fn rgroup_labels(&self) -> Result<Vec<String>> {
        Ok(self.rgroups_as_columns()?.into_keys().collect())
    }

    /// `outputCoreMolecule`: the matched core with the unused R group
    /// attachment points removed.
    fn output_core_molecule(&self, matched: &RGroupMatch, used: &UsedLabelMap) -> Result<ROMol> {
        let mut core = matched.matched_core.lock().expect("core lock").clone();
        let params = &self.data.params;
        const RETAINED: &str = "__rgdRetainedLabel";
        for idx in (0..core.num_atoms()).rev() {
            if core.atom(idx).atomic_num() != 0 {
                continue;
            }
            let label = self.data.get_rlabel(core.atom(idx))?;
            if core.atom_neighbors(idx).next().is_none() {
                continue;
            }
            let is_user_defined = used.has(label) && used.is_user_defined(label);
            let is_used = used.has(label) && used.is_used(label);
            if !is_used && (!is_user_defined || params.remove_all_hydrogen_rgroups_and_labels) {
                core.atom_mut(idx).set_atomic_num(1);
                core.update_atom_property_cache(idx, false)?;
            } else {
                // marked, as removing hydrogens renumbers the atoms
                core.atom_mut(idx).dict_mut().set_val(RETAINED, label);
            }
        }
        if params.remove_hydrogens_post_match {
            let _blocker = LogStateSetter::new();
            remove_hs(&mut core, &RemoveHsParameters::default(), false)?;
            core.update_property_cache(false)?;
        }
        for idx in 0..core.num_atoms() {
            let Some(label) = int_prop(core.atom(idx), RETAINED) else {
                continue;
            };
            core.atom_mut(idx).clear_prop(RETAINED);
            if core.num_conformers() == 0 || (used.has(label) && used.is_user_defined(label)) {
                continue;
            }
            let old_label = self
                .data
                .final_rlabel_mapping
                .iter()
                .find(|&(_, &v)| v == label)
                .map(|(&k, _)| k);
            if old_label.is_some_and(|old| matched.rgroups.contains_key(&old)) {
                let neighbor = core.atom_neighbors(idx).next();
                if let Some(neighbor) = neighbor {
                    set_terminal_atom_coords(&mut core, idx, neighbor)?;
                }
            }
        }
        if !core.ring_info().is_initialized() {
            crate::graphmol::find_rings::symmetrize_sssr(
                &mut core,
                crate::graphmol::find_rings::SymmetrizeSSSRAlgorithm::Default,
                true,
                false,
                false,
            )?;
        }
        Ok(core)
    }

    /// `getRGroupsAsRows`
    pub fn rgroups_as_rows(&self) -> Result<RGroupRows> {
        let permutation = self.data.current_best_permutation()?;
        let used_label_map = UsedLabelMap::new(&self.data.final_rlabel_mapping);
        let mut rows = Vec::new();
        for it in &permutation {
            let mut seen = used_label_map.clone();
            let mut row = RGroupRow::new();
            if self.data.params.include_target_mol_in_results {
                if let Some(target) =
                    it.target_molecule_for_highlights(self.data.params.remove_hydrogens_post_match)?
                {
                    row.entry(RGroupData::MOL_LABEL.to_string())
                        .or_insert(target);
                }
            }
            for (label, rgroup) in &it.rgroups {
                let real_label = self.real_label(*label)?;
                seen.set_used(real_label);
                row.entry(RGroupData::rgroup_label(real_label))
                    .or_insert_with(|| rgroup.read().expect("rgroup lock").combined_mol.clone());
            }
            row.entry(RGroupData::CORE_LABEL.to_string())
                .or_insert(self.output_core_molecule(it, &seen)?);
            rows.push(row);
        }
        Ok(rows)
    }

    fn real_label(&self, label: i32) -> Result<i32> {
        self.data
            .final_rlabel_mapping
            .get(&label)
            .copied()
            .ok_or_else(|| {
                RDKitError::invariant("unprocessed rlabel, please call process() first.")
            })
    }

    /// `getRGroupsAsColumns`
    pub fn rgroups_as_columns(&self) -> Result<RGroupColumns> {
        let permutation = self.data.current_best_permutation()?;
        let mut groups = RGroupColumns::new();
        let mut with_real_mol: BTreeSet<String> = BTreeSet::new();
        with_real_mol.insert(RGroupData::CORE_LABEL.to_string());
        if self.data.params.include_target_mol_in_results {
            with_real_mol.insert(RGroupData::MOL_LABEL.to_string());
        }
        let used_label_map = UsedLabelMap::new(&self.data.final_rlabel_mapping);
        for it in &permutation {
            let mut seen = used_label_map.clone();
            if self.data.params.include_target_mol_in_results {
                if let Some(target) =
                    it.target_molecule_for_highlights(self.data.params.remove_hydrogens_post_match)?
                {
                    groups
                        .entry(RGroupData::MOL_LABEL.to_string())
                        .or_default()
                        .push(target);
                }
            }
            for (label, rgroup) in &it.rgroups {
                let real_label = self.real_label(*label)?;
                let rg = rgroup.read().expect("rgroup lock");
                if !rg.combined_mol.has_prop(super::utils::DONE) {
                    return Err(RDKitError::invariant("Not done! Call process()"));
                }
                if seen.is_used(real_label) {
                    return Err(RDKitError::invariant(
                        "R group label appears multiple times!",
                    ));
                }
                seen.set_used(real_label);
                let name = RGroupData::rgroup_label(real_label);
                groups
                    .entry(name.clone())
                    .or_default()
                    .push(rg.combined_mol.clone());
                with_real_mol.insert(name);
            }
            groups
                .entry(RGroupData::CORE_LABEL.to_string())
                .or_default()
                .push(self.output_core_molecule(it, &seen)?);
            for &real_label in self.data.final_rlabel_mapping.values() {
                if !seen.is_used(real_label) {
                    groups
                        .entry(RGroupData::rgroup_label(real_label))
                        .or_default()
                        .push(RWMol::new());
                }
            }
        }
        groups.retain(|name, _| with_real_mol.contains(name));
        Ok(groups)
    }
}

/// Drops the unset (-1) entries.
fn trim_unset(indices: Vec<i32>) -> Vec<i32> {
    indices.into_iter().filter(|&i| i != -1).collect()
}

/// `splitNonUniqueMatches`: groups the matches covering the same atoms.
fn split_non_unique_matches(tmatches: &[MatchVect], n_atoms: u32) -> Vec<Vec<MatchVect>> {
    let mut out: Vec<Vec<MatchVect>> = Vec::new();
    let mut atom_sets: Vec<Vec<bool>> = Vec::new();
    for m in tmatches {
        let mut set = vec![false; n_atoms as usize];
        for &(_, t) in m {
            set[t as usize] = true;
        }
        match atom_sets.iter().position(|s| *s == set) {
            None => {
                atom_sets.push(set);
                out.push(vec![m.clone()]);
            }
            Some(i) => out[i].push(m.clone()),
        }
    }
    out
}

/// `RGroupDecompose`: decomposes `mols` on `cores`, returning the rows and
/// the indices of the molecules no core matched.
pub fn rgroup_decompose_rows(
    cores: &[&ROMol],
    mols: &[&ROMol],
    options: RGroupDecompositionParameters,
) -> Result<(RGroupRows, Vec<u32>)> {
    let mut decomp = RGroupDecomposition::new_multi(cores, options)?;
    let unmatched = decompose(&mut decomp, mols)?;
    Ok((decomp.rgroups_as_rows()?, unmatched))
}

/// `RGroupDecompose` with the result as columns.
pub fn rgroup_decompose_columns(
    cores: &[&ROMol],
    mols: &[&ROMol],
    options: RGroupDecompositionParameters,
) -> Result<(RGroupColumns, Vec<u32>)> {
    let mut decomp = RGroupDecomposition::new_multi(cores, options)?;
    let unmatched = decompose(&mut decomp, mols)?;
    Ok((decomp.rgroups_as_columns()?, unmatched))
}

fn decompose(decomp: &mut RGroupDecomposition, mols: &[&ROMol]) -> Result<Vec<u32>> {
    let t0 = Instant::now();
    let mut unmatched = Vec::new();
    for (i, mol) in mols.iter().enumerate() {
        if decomp.add(mol)? == -1 {
            unmatched.push(i as u32);
        }
        check_for_timeout(t0, decomp.params().timeout)?;
    }
    decomp.process();
    Ok(unmatched)
}

/// Quiet the unused import lint until the matching code needs it.
#[allow(dead_code)]
fn _unused(_: &SubstructMatchParameters) {}
