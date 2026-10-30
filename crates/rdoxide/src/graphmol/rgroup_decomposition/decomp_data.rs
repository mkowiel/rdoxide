//! The state of a decomposition: cores, matches, scoring and relabelling
//! (upstream `RGroupDecomposition/RGroupDecompData.{h,cpp}`).

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use super::core::RCore;
use super::data::{RGroupData, RGroupMatch};
use super::fingerprint_score::{FingerprintVarianceScoreData, fingerprint_variance_score};
use super::ga::RGroupGa;
use super::params::{
    RGroupDecompositionParameters, rgroup_labelling, rgroup_matching, rgroup_score,
};
use super::score::{CartesianProduct, RGroupScorer, check_for_timeout};
use super::utils::{
    DONE, EMPTY_CORE_LABEL, RGROUP_INPUT_DUMMY, RLABEL, RLABEL_CORE_INDEX, RLABEL_TYPE,
    SIDECHAIN_RLABELS, UNLABELED_CORE_ATTACHMENT, get_rlabels, int_prop,
    is_any_atom_with_multiple_neighbors_or_not_user_rlabel,
};
use crate::graphmol::atom::{Atom, set_atom_rlabel};
use crate::graphmol::bond::BondType;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs, set_terminal_atom_coords};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::log::LogStateSetter;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `RGroupDecompositionProcessResult`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RGroupDecompositionProcessResult {
    pub success: bool,
    pub score: f64,
}

/// Labels handed out while relabelling.
#[derive(Debug, Default)]
pub(crate) struct UsedLabels {
    labels_used: BTreeSet<i32>,
}

impl UsedLabels {
    fn add(&mut self, rlabel: i32) -> bool {
        self.labels_used.insert(rlabel)
    }

    /// The smallest unused positive label.
    fn next(&mut self) -> i32 {
        let mut i = 1;
        while self.labels_used.contains(&i) {
            i += 1;
        }
        self.labels_used.insert(i);
        i
    }
}

/// `RGroupDecompData`
#[derive(Debug)]
pub struct RGroupDecompData {
    pub cores: BTreeMap<i32, RCore>,
    /// new "cores" found along the way, by SMILES
    pub new_cores: BTreeMap<String, i32>,
    pub new_core_label: i32,
    pub previous_match_size: usize,
    pub prune_permutations: bool,
    pub params: RGroupDecompositionParameters,
    pub matches: Vec<Vec<RGroupMatch>>,
    pub labels: BTreeSet<i32>,
    pub permutation: Vec<usize>,
    prune_length: usize,
    pruned_fingerprint_variance_score_data: FingerprintVarianceScoreData,
    pub final_rlabel_mapping: BTreeMap<i32, i32>,
    scorer: RGroupScorer,
}

impl RGroupDecompData {
    pub fn new(input_cores: &[&ROMol], params: RGroupDecompositionParameters) -> Result<Self> {
        let mut data = RGroupDecompData {
            cores: BTreeMap::new(),
            new_cores: BTreeMap::new(),
            new_core_label: EMPTY_CORE_LABEL,
            previous_match_size: 0,
            prune_permutations: true,
            params,
            matches: Vec::new(),
            labels: BTreeSet::new(),
            permutation: Vec::new(),
            prune_length: 0,
            pruned_fingerprint_variance_score_data: FingerprintVarianceScoreData::default(),
            final_rlabel_mapping: BTreeMap::new(),
            scorer: RGroupScorer::default(),
        };
        for core in input_cores {
            data.add_input_core(core)?;
        }
        data.prepare_cores()?;
        Ok(data)
    }

    fn add_input_core(&mut self, input_core: &ROMol) -> Result<()> {
        if self.params.do_enumeration {
            let bundle = crate::graphmol::mol_enumerator::enumerate_default(input_core, 0)?;
            if !bundle.is_empty() {
                for core in bundle.mols() {
                    self.add_core(core)?;
                }
                return Ok(());
            }
        }
        self.add_core(input_core)
    }

    fn add_core(&mut self, input_core: &ROMol) -> Result<()> {
        let mut core = input_core.clone();
        if self.params.allow_multiple_rgroups_on_unlabelled && !self.params.only_match_at_rgroups {
            self.params
                .add_dummy_atoms_to_unlabelled_core_atoms(&mut core)?;
        }
        let key = self.cores.len() as i32;
        self.cores.insert(key, RCore::new(core)?);
        Ok(())
    }

    fn prepare_cores(&mut self) -> Result<()> {
        let keys: Vec<i32> = self.cores.keys().copied().collect();
        for key in keys {
            let align_core = (key != 0).then(|| self.cores[&0].core.clone());
            let core = self.cores.get_mut(&key).expect("core key");
            if !self
                .params
                .prepare_core(&mut core.core, align_core.as_ref())?
            {
                return Err(RDKitError::invariant("Could not prepare at least one core"));
            }
            core.init()?;
            core.labelled_core = core.core.clone();
        }
        Ok(())
    }

    /// `setRlabel`
    pub fn set_rlabel(&self, atom: &mut Atom, rlabel: i32) -> Result<()> {
        if rlabel <= 0 {
            return Err(RDKitError::precondition("RLabels must be >0"));
        }
        let labelling = self.params.rgroup_labelling;
        if labelling & rgroup_labelling::ATOM_MAP != 0 {
            atom.set_atom_map_num(rlabel, false)?;
        }
        if labelling & rgroup_labelling::MDL_RGROUP != 0 {
            atom.dict_mut()
                .set_val(cp::dummyLabel, format!("R{rlabel}"));
            set_atom_rlabel(atom, rlabel)?;
        } else {
            atom.clear_prop(cp::dummyLabel);
        }
        if labelling & rgroup_labelling::ISOTOPE != 0 {
            atom.set_isotope(rlabel as u32);
        }
        Ok(())
    }

    /// `getRlabel`
    pub fn get_rlabel(&self, atom: &Atom) -> Result<i32> {
        let labelling = self.params.rgroup_labelling;
        if labelling & rgroup_labelling::ATOM_MAP != 0 {
            return Ok(atom.atom_map_num());
        }
        if labelling & rgroup_labelling::ISOTOPE != 0 {
            return Ok(atom.isotope() as i32);
        }
        if labelling & rgroup_labelling::MDL_RGROUP != 0 {
            if let Some(label) = atom.get_prop(cp::_MolFileRLabel).and_then(RDValue::as_uint) {
                return Ok(label as i32);
            }
        }
        Err(RDKitError::invariant("no valid r label found"))
    }

    fn score_from_pruned_data(&mut self, permutation: &[usize], reset: bool) -> Result<f64> {
        if self.params.score_method != rgroup_score::FINGERPRINT_VARIANCE {
            return Err(RDKitError::precondition(
                "Scoring method is not fingerprint variance!",
            ));
        }
        if permutation.len() < self.prune_length {
            return Err(RDKitError::precondition("Illegal permutation prune length"));
        }
        if (permutation.len() as f64) < self.prune_length as f64 * 1.5 {
            for pos in self.prune_length..permutation.len() {
                self.pruned_fingerprint_variance_score_data
                    .add_variance_data(pos, permutation[pos], &self.matches, &self.labels)?;
            }
            let score = self
                .pruned_fingerprint_variance_score_data
                .fingerprint_variance_group_score()?;
            if reset {
                for pos in self.prune_length..permutation.len() {
                    self.pruned_fingerprint_variance_score_data
                        .remove_variance_data(pos, permutation[pos], &self.matches, &self.labels)?;
                }
            } else {
                self.prune_length = permutation.len();
            }
            Ok(score)
        } else if reset {
            fingerprint_variance_score(permutation, &self.matches, &self.labels, None)
        } else {
            self.pruned_fingerprint_variance_score_data.clear();
            self.prune_length = permutation.len();
            fingerprint_variance_score(
                permutation,
                &self.matches,
                &self.labels,
                Some(&mut self.pruned_fingerprint_variance_score_data),
            )
        }
    }

    /// `prune`: keeps only the current best match of each molecule.
    fn prune(&mut self) -> Result<()> {
        if self.permutation.len() > self.matches.len() {
            return Err(RDKitError::precondition(
                "permutation.size() should be <= matches.size()",
            ));
        }
        let offset = self.matches.len() - self.permutation.len();
        for (mol_idx, &p) in self.permutation.iter().enumerate() {
            let keep = self.matches[mol_idx + offset]
                .get(p)
                .cloned()
                .ok_or_else(|| RDKitError::Index("permutation index out of range".into()))?;
            self.matches[mol_idx + offset] = vec![keep];
        }
        self.permutation = vec![0; self.permutation.len()];
        if self.params.score_method == rgroup_score::FINGERPRINT_VARIANCE
            && self.params.matching_strategy != rgroup_matching::GA
        {
            let permutation = self.permutation.clone();
            self.score_from_pruned_data(&permutation, false)?;
        }
        Ok(())
    }

    /// `GetCurrentBestPermutation`: the best match of every molecule, without
    /// the R groups that are hydrogen in all of them.
    pub fn current_best_permutation(&self) -> Result<Vec<RGroupMatch>> {
        let remove_all_hydrogen_rgroups = self.params.remove_all_hydrogen_rgroups
            || self.params.remove_all_hydrogen_rgroups_and_labels;
        let is_pruned = self.permutation.len() < self.matches.len();
        let mut results = Vec::with_capacity(self.matches.len());
        for (i, candidates) in self.matches.iter().enumerate() {
            let pi = if is_pruned {
                0
            } else {
                *self
                    .permutation
                    .get(i)
                    .ok_or_else(|| RDKitError::Index("vector::_M_range_check".into()))?
            };
            results.push(
                candidates
                    .get(pi)
                    .cloned()
                    .ok_or_else(|| RDKitError::Index("vector::_M_range_check".into()))?,
            );
        }

        // label -> cores having it
        let mut label_cores: BTreeMap<i32, BTreeSet<i32>> = BTreeMap::new();
        let mut cores_visited = BTreeSet::new();
        for position in &results {
            if cores_visited.insert(position.core_idx) {
                if let Some(core) = self.cores.get(&position.core_idx) {
                    for label in get_rlabels(&core.core)?.keys() {
                        label_cores
                            .entry(*label)
                            .or_default()
                            .insert(position.core_idx);
                    }
                }
            }
        }

        let mut labels_to_erase = BTreeSet::new();
        for &label in &self.labels {
            if label > 0 && !remove_all_hydrogen_rgroups {
                continue;
            }
            let all_h = results.iter().all(|position| {
                let has_core = label_cores
                    .get(&label)
                    .is_some_and(|cores| cores.contains(&position.core_idx));
                !(has_core
                    && position
                        .rgroups
                        .get(&label)
                        .is_some_and(|rg| !rg.read().expect("rgroup lock").is_hydrogen))
            });
            if all_h {
                labels_to_erase.insert(label);
                for position in &mut results {
                    position.rgroups.remove(&label);
                }
            }
        }

        for position in &results {
            let mut core = position.matched_core.lock().expect("core lock");
            for idx in 0..core.num_atoms() {
                let atom = core.atom(idx);
                let Some(atom_label) = int_prop(atom, RLABEL) else {
                    continue;
                };
                if atom.atomic_num() != 0 {
                    continue;
                }
                if atom_label > 0 && !self.params.remove_all_hydrogen_rgroups_and_labels {
                    continue;
                }
                if labels_to_erase.contains(&atom_label) {
                    let atom = core.atom_mut(idx);
                    atom.set_atomic_num(1);
                    atom.clear_prop(RLABEL);
                    atom.clear_prop(RLABEL_TYPE);
                    atom.clear_prop(UNLABELED_CORE_ATTACHMENT);
                    core.update_atom_property_cache(idx, false)?;
                }
            }
        }
        Ok(results)
    }

    /// `addCoreUserLabels`
    fn add_core_user_labels(core: &RWMol, user_labels: &mut BTreeSet<i32>) -> Result<()> {
        for label in get_rlabels(core)?.keys() {
            if *label > 0 {
                user_labels.insert(*label);
            }
        }
        Ok(())
    }

    /// `addAtoms`: attaches the new atoms to the atoms they are paired with.
    fn add_atoms(mol: &mut RWMol, atoms_to_add: Vec<(u32, Atom)>) -> Result<()> {
        for (anchor, atom) in atoms_to_add {
            let new_idx = mol.add_atom(atom, false);
            mol.add_bond(anchor, new_idx, BondType::Single)?;
            if mol.num_conformers() > 0 {
                set_terminal_atom_coords(mol, new_idx, anchor)?;
            }
        }
        Ok(())
    }

    /// `replaceHydrogenCoreDummy`: when the R group at `current_label` is
    /// hydrogen, a hydrogen neighbor of `atom` becomes the dummy instead.
    fn replace_hydrogen_core_dummy(
        &self,
        matched: &RGroupMatch,
        core: &mut RWMol,
        atom_idx: u32,
        current_label: i32,
        rlabel: i32,
    ) -> Result<bool> {
        let Some(group) = matched.rgroups.get(&current_label) else {
            return Ok(false);
        };
        if !group.read().expect("rgroup lock").is_hydrogen {
            return Ok(false);
        }
        let neighbors: Vec<u32> = core.atom_neighbors(atom_idx).collect();
        for n in neighbors {
            if core.atom(n).atomic_num() == 1 {
                let neighbor = core.atom_mut(n);
                neighbor.set_atomic_num(0);
                self.set_rlabel(neighbor, rlabel)?;
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Adds a dummy atom to `atom_idx` unless a hydrogen could be turned into
    /// it (when matching) or the match has no R group at `label`.
    fn attach_dummy(
        &self,
        core: &mut RWMol,
        atom_idx: u32,
        label: i32,
        rlabel: i32,
        matched: Option<&RGroupMatch>,
        atoms_to_add: &mut Vec<(u32, Atom)>,
    ) -> Result<()> {
        let mut add_new = true;
        if let Some(matched) = matched {
            add_new = !self.replace_hydrogen_core_dummy(matched, core, atom_idx, label, rlabel)?;
            if add_new && !matched.rgroups.contains_key(&label) {
                add_new = false;
            }
        }
        if add_new {
            let mut new_atom = Atom::new(0);
            self.set_rlabel(&mut new_atom, rlabel)?;
            atoms_to_add.push((atom_idx, new_atom));
        }
        Ok(())
    }

    /// `relabelCore`
    fn relabel_core(
        &self,
        core: &mut RWMol,
        mappings: &mut BTreeMap<i32, i32>,
        used_labels: &mut UsedLabels,
        index_labels: &BTreeSet<i32>,
        extra_atom_rlabels: &BTreeMap<i32, Vec<i32>>,
        matched: Option<&RGroupMatch>,
    ) -> Result<()> {
        let atoms = get_rlabels(core)?;
        let mut atoms_to_add: Vec<(u32, Atom)> = Vec::new();
        for (&user_label, &atom_idx) in &atoms {
            if user_label < 0 {
                continue;
            }
            mappings.insert(user_label, user_label);
            used_labels.add(user_label);
            if core.atom(atom_idx).atomic_num() == 0 && core.atom_degree(atom_idx) == 1 {
                let mut atom = core.atom(atom_idx).clone();
                self.set_rlabel(&mut atom, user_label)?;
                *core.atom_mut(atom_idx) = atom;
            } else {
                self.attach_dummy(
                    core,
                    atom_idx,
                    user_label,
                    user_label,
                    matched,
                    &mut atoms_to_add,
                )?;
            }
        }

        for &new_label in index_labels {
            let Some(&atom_idx) = atoms.get(&new_label) else {
                continue;
            };
            let rlabel = match mappings.get(&new_label) {
                Some(&r) => r,
                None => {
                    let r = used_labels.next();
                    mappings.insert(new_label, r);
                    r
                }
            };
            if core.atom(atom_idx).atomic_num() == 0
                && !is_any_atom_with_multiple_neighbors_or_not_user_rlabel(core, atom_idx)
            {
                let mut atom = core.atom(atom_idx).clone();
                self.set_rlabel(&mut atom, rlabel)?;
                *core.atom_mut(atom_idx) = atom;
            } else {
                self.attach_dummy(
                    core,
                    atom_idx,
                    new_label,
                    rlabel,
                    matched,
                    &mut atoms_to_add,
                )?;
            }
        }

        for (label, extra) in extra_atom_rlabels {
            let Some(&atom_idx) = atoms.get(label) else {
                continue; // label not used in the rgroup
            };
            for _ in 0..extra.len() {
                let rlabel = used_labels.next();
                if core.atom(atom_idx).atomic_num() <= 1 {
                    return Err(RDKitError::invariant(
                        "Multiple attachments to a dummy (or hydrogen) is weird.",
                    ));
                }
                let mut new_atom = Atom::new(0);
                self.set_rlabel(&mut new_atom, rlabel)?;
                atoms_to_add.push((atom_idx, new_atom));
            }
        }

        Self::add_atoms(core, atoms_to_add)?;
        for &atom_idx in atoms.values() {
            let atom = core.atom_mut(atom_idx);
            atom.clear_prop(RLABEL);
            atom.clear_prop(RLABEL_TYPE);
        }
        core.update_property_cache(false)
    }

    /// `relabelRGroup`
    fn relabel_rgroup(&self, rgroup: &mut RGroupData, mappings: &BTreeMap<i32, i32>) -> Result<()> {
        if rgroup.combined_mol.has_prop(DONE) {
            rgroup.labelled = true;
            return Ok(());
        }
        let mol = &mut rgroup.combined_mol;
        mol.dict_mut().set_val(DONE, true);
        let mut atoms_to_add: Vec<(u32, Atom)> = Vec::new();
        let mut core_index_to_atomic_num: BTreeMap<i32, u32> = BTreeMap::new();
        for idx in 0..mol.num_atoms() {
            if let Some(RDValue::VecInt(rlabels)) =
                mol.atom(idx).get_prop(SIDECHAIN_RLABELS).cloned()
            {
                mol.atom_mut(idx).set_isotope(0);
                for rlabel in rlabels {
                    let label = *mappings
                        .get(&rlabel)
                        .ok_or_else(|| RDKitError::invariant("Unprocessed mapping"))?;
                    let atom = mol.atom_mut(idx);
                    if atom.atomic_num() == 0 {
                        if !atom.has_prop(RGROUP_INPUT_DUMMY) {
                            self.set_rlabel(atom, label)?;
                        }
                    } else if atom.has_prop(RLABEL_CORE_INDEX) {
                        atom.set_atomic_num(0);
                        self.set_rlabel(atom, label)?;
                    } else {
                        let mut new_atom = Atom::new(0);
                        self.set_rlabel(&mut new_atom, label)?;
                        atoms_to_add.push((idx, new_atom));
                    }
                }
            }
            if let Some(core_index) = int_prop(mol.atom(idx), RLABEL_CORE_INDEX) {
                let atom = mol.atom_mut(idx);
                core_index_to_atomic_num.insert(core_index, atom.atomic_num());
                atom.set_atomic_num(0);
            }
        }
        Self::add_atoms(mol, atoms_to_add)?;
        if self.params.remove_hydrogens_post_match {
            let _blocker = LogStateSetter::new();
            remove_hs(mol, &RemoveHsParameters::default(), false)?;
        }
        mol.update_property_cache(false)?;
        rgroup.labelled = true;
        let mol = &mut rgroup.combined_mol;
        for idx in 0..mol.num_atoms() {
            if let Some(core_index) = int_prop(mol.atom(idx), RLABEL_CORE_INDEX) {
                let atomic_num = core_index_to_atomic_num
                    .get(&core_index)
                    .copied()
                    .unwrap_or(0);
                let atom = mol.atom_mut(idx);
                atom.set_atomic_num(atomic_num as i32);
                atom.set_no_implicit(true);
                atom.clear_prop(RLABEL_CORE_INDEX);
            }
            mol.atom_mut(idx).clear_prop(SIDECHAIN_RLABELS);
        }
        Ok(())
    }

    /// `relabel`: gives the R groups their final labels.
    fn relabel(&mut self) -> Result<()> {
        let best = self.current_best_permutation()?;
        let mut user_labels = BTreeSet::new();
        let mut index_labels = BTreeSet::new();
        let mut extra_atom_rlabels: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
        for it in &best {
            for (&label, rgroup) in &it.rgroups {
                if label > 0 {
                    user_labels.insert(label);
                }
                if label < 0 && !self.params.only_match_at_rgroups {
                    index_labels.insert(label);
                }
                for (rlabel, count) in rgroup.read().expect("rgroup lock").num_bonds_to_rlabels() {
                    if count > 1 {
                        extra_atom_rlabels.insert(rlabel, vec![0; (count - 1) as usize]);
                    }
                }
            }
        }
        for core in self.cores.values_mut() {
            core.labelled_core = core.core.clone();
            Self::add_core_user_labels(&core.labelled_core, &mut user_labels)?;
        }

        let mut mapping = BTreeMap::new();
        let mut used_labels = UsedLabels::default();
        for &label in &user_labels {
            used_labels.add(label);
        }
        let mut labelled_cores = BTreeMap::new();
        for (&key, core) in &self.cores {
            let mut labelled = core.labelled_core.clone();
            self.relabel_core(
                &mut labelled,
                &mut mapping,
                &mut used_labels,
                &index_labels,
                &extra_atom_rlabels,
                None,
            )?;
            labelled_cores.insert(key, labelled);
        }
        for (key, labelled) in labelled_cores {
            self.cores.get_mut(&key).expect("core key").labelled_core = labelled;
        }

        for it in &best {
            for rgroup in it.rgroups.values() {
                let mut rg = rgroup.write().expect("rgroup lock");
                self.relabel_rgroup(&mut rg, &mapping)?;
            }
            let mut matched_core = it.matched_core.lock().expect("core lock");
            self.relabel_core(
                &mut matched_core,
                &mut mapping,
                &mut used_labels,
                &index_labels,
                &extra_atom_rlabels,
                Some(it),
            )?;
        }

        let unique_mapped: BTreeSet<i32> = mapping.values().copied().collect();
        if mapping.len() != unique_mapped.len() {
            return Err(RDKitError::invariant(
                "Error in uniqueness of final RLabel mapping",
            ));
        }
        if unique_mapped.len() != user_labels.len() + index_labels.len() {
            return Err(RDKitError::invariant("Error in final RMapping size"));
        }
        self.final_rlabel_mapping = mapping;
        Ok(())
    }

    /// `score`
    fn score(&mut self, permutation: &[usize], pruned: bool) -> Result<f64> {
        match self.params.score_method {
            rgroup_score::MATCH => {
                self.scorer
                    .match_score(permutation, &self.matches, &self.labels)
            }
            rgroup_score::FINGERPRINT_VARIANCE => {
                if pruned {
                    self.score_from_pruned_data(permutation, true)
                } else {
                    fingerprint_variance_score(permutation, &self.matches, &self.labels, None)
                }
            }
            _ => Ok(f64::NAN),
        }
    }

    /// `process`: scores the permutations of the R group assignments and keeps
    /// the best; with `finalize` the R groups also get their final labels.
    pub fn process(
        &mut self,
        prune_matches: bool,
        finalize: bool,
    ) -> Result<RGroupDecompositionProcessResult> {
        if self.matches.is_empty() {
            return Ok(RGroupDecompositionProcessResult {
                success: false,
                score: -1.0,
            });
        }
        let t0 = Instant::now();
        self.scorer.start_processing();

        if self.params.matching_strategy == rgroup_matching::GA {
            let mut ga = RGroupGa::new(
                &self.matches,
                &self.labels,
                &self.params,
                (self.params.timeout >= 0.0).then_some(t0),
            );
            if u64::from(ga.number_permutations()) < 100 * ga.population_size() as u64 {
                self.params.matching_strategy = rgroup_matching::EXHAUSTIVE;
            } else {
                let result = if self.params.ga_number_runs > 1 {
                    let results = ga.run_batch()?;
                    // the first of the best scores wins
                    results.into_iter().reduce(|best, r| {
                        if best.scorer.best_score() < r.scorer.best_score() {
                            r
                        } else {
                            best
                        }
                    })
                } else {
                    Some(ga.run(1)?)
                };
                if let Some(result) = result {
                    self.scorer = result.scorer;
                }
            }
        }

        let mut iterator = None;
        if self.params.matching_strategy != rgroup_matching::GA {
            let mut offset = 0;
            if prune_matches && self.params.score_method != rgroup_score::FINGERPRINT_VARIANCE {
                offset = self.previous_match_size;
            }
            self.previous_match_size = self.matches.len();
            let sizes: Vec<usize> = self.matches[offset..].iter().map(Vec::len).collect();
            self.permutation = vec![0; sizes.len()];
            let mut cartesian = CartesianProduct::new(&sizes);
            let mut count = 0usize;
            while cartesian.next() {
                if count > cartesian.max_permutations {
                    return Err(RDKitError::Value("next() did not finish".into()));
                }
                let permutation = cartesian.permutation.clone();
                let new_score = self.score(&permutation, true)?;
                // heuristic to overcome floating point comparison issues
                if (new_score - self.scorer.best_score()).abs() < 1e-6 {
                    self.scorer.push_tie_to_store(&permutation);
                } else if new_score > self.scorer.best_score() {
                    self.scorer.set_best_permutation(&permutation, new_score);
                    self.scorer.clear_tie_store();
                    self.scorer.push_tie_to_store(&permutation);
                }
                count += 1;
            }
            iterator = Some(cartesian);
        }

        if self.scorer.tie_store_size() > 1 {
            self.scorer.break_ties(
                &self.matches,
                &self.labels,
                iterator.as_ref(),
                t0,
                self.params.timeout,
            )?;
            self.scorer.clear_tie_store();
        } else {
            check_for_timeout(t0, self.params.timeout)?;
        }
        self.permutation = self.scorer.best_permutation().to_vec();
        if prune_matches || finalize {
            self.prune()?;
        }
        if finalize {
            self.relabel()?;
        }
        Ok(RGroupDecompositionProcessResult {
            success: true,
            score: self.scorer.best_score(),
        })
    }
}
