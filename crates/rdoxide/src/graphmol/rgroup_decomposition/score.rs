//! Scoring of the permutations of R group assignments (upstream
//! `RGroupDecomposition/RGroupScore.{h,cpp}`).

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::time::Instant;

use super::data::RGroupMatch;
use crate::{RDKitError, Result};

/// `checkForTimeout`: an error once `timeout` seconds have passed since `t0`;
/// timeouts of 0 or less never expire.
pub fn check_for_timeout(t0: Instant, timeout: f64) -> Result<()> {
    if timeout > 0.0 && t0.elapsed().as_secs_f64() >= timeout {
        return Err(RDKitError::Runtime("operation timed out".into()));
    }
    Ok(())
}

/// `CartesianProduct`: iterates over all the permutations that pick one
/// element of each row.
#[derive(Debug, Clone)]
pub struct CartesianProduct {
    pub permutation: Vec<usize>,
    sizes: Vec<usize>,
    bases: VecDeque<usize>,
    pub max_permutations: usize,
    permutation_count: usize,
}

impl CartesianProduct {
    pub fn new(sizes: &[usize]) -> Self {
        let mut bases = VecDeque::new();
        let mut max_permutations = 1usize;
        for &size in sizes {
            bases.push_front(max_permutations);
            max_permutations = max_permutations.wrapping_mul(size);
        }
        CartesianProduct {
            permutation: vec![0; sizes.len()],
            sizes: sizes.to_vec(),
            bases,
            max_permutations,
            permutation_count: 0,
        }
    }

    /// Moves to the next permutation; `false` after the last one.
    pub fn next(&mut self) -> bool {
        self.permutation_count += 1;
        if self.permutation_count == 1 {
            return true;
        }
        self.increment(0)
    }

    pub fn value_of(&self, p: &[usize]) -> usize {
        p.iter()
            .enumerate()
            .map(|(i, &v)| self.bases[i].wrapping_mul(v))
            .fold(0usize, usize::wrapping_add)
    }

    fn increment(&mut self, row: usize) -> bool {
        if self.permutation_count > self.max_permutations || row >= self.sizes.len() {
            return false;
        }
        self.permutation[row] += 1;
        if self.permutation[row] > self.sizes[row] - 1 {
            self.permutation[row] = 0;
            return self.increment(row + 1);
        }
        true
    }
}

#[derive(Debug, Clone, Default)]
struct RLabelData {
    num_rgroups: usize,
    match_set_vect: Vec<BTreeMap<String, u32>>,
    linker_match_set: BTreeMap<BTreeSet<i32>, usize>,
}

/// `RGroupScorer::State`
#[derive(Debug, Clone, Default)]
struct State {
    n: usize,
    num_added_rgroups: usize,
    num_matched_user_rgroups: usize,
    heavy_count_per_label: BTreeMap<i32, i32>,
    label_data_map: BTreeMap<i32, RLabelData>,
    permutation: Vec<usize>,
}

impl State {
    /// `computeTieBreakingCriteria`
    fn compute_tie_breaking_criteria(
        &mut self,
        matches: &[Vec<RGroupMatch>],
        ordered_labels: &[i32],
        heavy_counts: &mut [i32],
    ) {
        let offset = matches.len() - self.permutation.len();
        for (i, &label) in ordered_labels.iter().enumerate() {
            for (m, &p) in self.permutation.iter().enumerate() {
                let found = matches[m + offset][p].rgroups.get(&label);
                if found.is_some_and(|rg| !rg.read().expect("rgroup lock").is_hydrogen) {
                    if label < 0 && heavy_counts[i] == 0 {
                        self.num_added_rgroups += 1;
                    } else if label > 0 {
                        self.num_matched_user_rgroups += 1;
                    }
                    heavy_counts[i] += 1;
                }
            }
        }
    }
}

/// `RGroupScorer`: scores permutations, remembers the best ones and breaks
/// ties between them.
#[derive(Debug, Clone, Default)]
pub struct RGroupScorer {
    best_score: f64,
    current: State,
    initial: State,
    saved: State,
    store: VecDeque<State>,
}

impl RGroupScorer {
    pub fn best_score(&self) -> f64 {
        self.best_score
    }

    pub fn best_permutation(&self) -> &[usize] {
        &self.saved.permutation
    }

    pub fn tie_store_size(&self) -> usize {
        self.store.len()
    }

    pub fn clear_tie_store(&mut self) {
        self.store.clear();
    }

    pub fn push_tie_to_store(&mut self, permutation: &[usize]) {
        self.current.permutation = permutation.to_vec();
        self.store.push_back(self.current.clone());
    }

    pub fn set_best_permutation(&mut self, permutation: &[usize], score: f64) {
        self.best_score = score;
        self.current.permutation = permutation.to_vec();
        self.saved = self.current.clone();
    }

    /// `startProcessing`
    pub fn start_processing(&mut self) {
        self.initial = self.saved.clone();
        self.best_score = -f64::MAX;
        self.clear_tie_store();
    }

    /// `matchScore`: the more R groups of the molecules are equal, the higher
    /// the score.
    pub fn match_score(
        &mut self,
        permutation: &[usize],
        matches: &[Vec<RGroupMatch>],
        labels: &BTreeSet<i32>,
    ) -> Result<f64> {
        if permutation.len() > matches.len() {
            return Err(RDKitError::precondition(
                "permutation.size() should be <= matches.size()",
            ));
        }
        let offset = matches.len() - permutation.len();
        self.current = self.initial.clone();
        let mut num_rgroups: BTreeMap<i32, usize> = BTreeMap::new();
        for (m, &p) in permutation.iter().enumerate() {
            for label in matches[m + offset][p].rgroups.keys() {
                let count = num_rgroups.entry(*label).or_insert(0);
                *count += 1;
                self.current.n = self.current.n.max(*count);
            }
        }
        let n = self.current.n;
        let mut score = 0.0;
        for &l in labels {
            let label_data = self.current.label_data_map.entry(l).or_default();
            for (m, &p) in permutation.iter().enumerate() {
                let Some(rg) = matches[m + offset][p].rgroups.get(&l) else {
                    continue;
                };
                let rg = rg.read().expect("rgroup lock");
                label_data.num_rgroups += 1;
                if rg.is_linker {
                    *label_data
                        .linker_match_set
                        .entry(rg.attachments.clone())
                        .or_insert(0) += 1;
                }
                for (i, smiles) in rg.smiles_vect.iter().enumerate() {
                    if i == label_data.match_set_vect.len() {
                        label_data.match_set_vect.push(BTreeMap::new());
                    }
                    *label_data.match_set_vect[i]
                        .entry(smiles.clone())
                        .or_insert(0) += 1;
                }
            }

            let mut temp_score = 0.0;
            let num_sets = label_data.match_set_vect.len();
            let num_with_label = label_data.num_rgroups;
            for match_set in &mut label_data.match_set_vect {
                if n > num_with_label {
                    match_set.insert(String::new(), (n - num_with_label) as u32);
                }
                let mut counts: Vec<u32> = match_set.values().copied().collect();
                counts.sort_by(|a, b| b.cmp(a));
                for (i, &count) in counts.iter().enumerate() {
                    let lscore = f64::from(count) / ((i + 1) * matches.len()) as f64;
                    temp_score += lscore * lscore;
                }
                temp_score /= num_sets as f64;
            }

            let mut max_linker_matches = 0usize;
            for (attachments, &count) in &label_data.linker_match_set {
                if (attachments.len() > 1 || count > 1) && attachments.len() > max_linker_matches {
                    max_linker_matches = attachments.len().max(count);
                }
            }
            let mut increment = 1.0;
            let mut linker_increment = 1.0;
            if max_linker_matches > 0 {
                linker_increment = max_linker_matches as f64 / matches.len() as f64;
            } else {
                increment = temp_score;
            }
            score += increment * linker_increment;
        }
        Ok(score)
    }

    /// `breakTies`: among the permutations of the tie store, prefers those
    /// matching more user R groups, adding fewer new ones, having more heavy
    /// R groups, and finally the last in iteration order.
    pub fn break_ties(
        &mut self,
        matches: &[Vec<RGroupMatch>],
        labels: &BTreeSet<i32>,
        iterator: Option<&CartesianProduct>,
        t0: Instant,
        timeout: f64,
    ) -> Result<()> {
        let mut max_perm_value = 0usize;
        self.current = self.saved.clone();
        self.current.num_added_rgroups = labels.len();
        // user labels ascending, then the others by decreasing label
        let mut ordered_labels: Vec<i32> = labels.iter().copied().filter(|&l| l >= 0).collect();
        ordered_labels.extend(labels.iter().rev().copied().filter(|&l| l < 0));
        let mut heavy_iter = self.current.heavy_count_per_label.iter();
        let mut largest_heavy_counts: Vec<i32> = Vec::with_capacity(labels.len());
        for &label in &ordered_labels {
            let mut count = 0;
            if let Some((&first, &second)) = heavy_iter.next() {
                if (first > 0) == (label > 0) {
                    count = second;
                }
            }
            largest_heavy_counts.push(count);
        }
        let initial_heavy_counts = largest_heavy_counts.clone();
        while let Some(mut state) = self.store.pop_front() {
            let mut heavy_counts = initial_heavy_counts.clone();
            state.compute_tie_breaking_criteria(matches, &ordered_labels, &mut heavy_counts);
            let perm_value = iterator.map_or(max_perm_value, |it| it.value_of(&state.permutation));
            let take = if state.num_matched_user_rgroups > self.current.num_matched_user_rgroups {
                true
            } else if state.num_matched_user_rgroups == self.current.num_matched_user_rgroups {
                if state.num_added_rgroups < self.current.num_added_rgroups {
                    true
                } else if state.num_added_rgroups == self.current.num_added_rgroups {
                    heavy_counts > largest_heavy_counts
                        || (heavy_counts == largest_heavy_counts && perm_value > max_perm_value)
                } else {
                    false
                }
            } else {
                false
            };
            if take {
                self.current = state;
                largest_heavy_counts = heavy_counts;
                max_perm_value = perm_value;
            }
            check_for_timeout(t0, timeout)?;
        }
        self.current.heavy_count_per_label.clear();
        for (&label, &count) in ordered_labels.iter().zip(&largest_heavy_counts) {
            self.current.heavy_count_per_label.insert(label, count);
        }
        self.saved = self.current.clone();
        Ok(())
    }
}
