//! `SearchResults.cpp`

use std::collections::HashMap;

use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

#[derive(Default)]
pub struct SearchResults {
    hit_molecules: Vec<ROMol>,
    mol_names: HashMap<String, usize>,
    max_num_results: u64,
    timed_out: bool,
    cancelled: bool,
    best_hit_found: Option<ROMol>,
}

fn name_of(m: &ROMol) -> String {
    m.get_prop(cp::_Name)
        .map(|v| v.to_rd_string())
        .unwrap_or_default()
}

fn sim_of(m: &ROMol) -> Option<f64> {
    m.get_prop("Similarity").and_then(|v| v.as_double())
}

impl SearchResults {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with(
        mols: Vec<ROMol>,
        best_hit: Option<ROMol>,
        max_num_res: u64,
        timed_out: bool,
        cancelled: bool,
    ) -> Self {
        SearchResults {
            hit_molecules: mols,
            mol_names: HashMap::new(),
            max_num_results: max_num_res,
            timed_out,
            cancelled,
            best_hit_found: best_hit,
        }
    }
    pub fn max_num_results(&self) -> u64 {
        self.max_num_results
    }
    pub fn hit_molecules(&self) -> &[ROMol] {
        &self.hit_molecules
    }
    pub fn best_hit(&self) -> Option<&ROMol> {
        self.best_hit_found.as_ref()
    }
    pub fn timed_out(&self) -> bool {
        self.timed_out
    }
    pub fn cancelled(&self) -> bool {
        self.cancelled
    }

    /// Merges `other` into this, keeping only molecules with unique names;
    /// of two of the same name the one with the higher "Similarity" is kept.
    pub fn merge_results(&mut self, other: &mut SearchResults) {
        self.max_num_results += other.max_num_results;
        self.timed_out |= other.timed_out;
        self.cancelled |= other.cancelled;
        if self.mol_names.is_empty() {
            for (i, m) in self.hit_molecules.iter().enumerate() {
                self.mol_names.insert(name_of(m), i);
            }
        }
        self.hit_molecules.reserve(other.hit_molecules.len());
        for mol in std::mem::take(&mut other.hit_molecules) {
            let name = name_of(&mol);
            match self.mol_names.get(&name) {
                None => {
                    self.mol_names.insert(name, self.hit_molecules.len());
                    self.hit_molecules.push(mol);
                }
                Some(&i) => {
                    if let (Some(cur), Some(new)) = (sim_of(&self.hit_molecules[i]), sim_of(&mol))
                        && new > cur
                    {
                        self.hit_molecules[i] = mol;
                    }
                }
            }
        }
        match (&self.best_hit_found, &other.best_hit_found) {
            (Some(a), Some(b)) => {
                if sim_of(a).unwrap_or(0.0) < sim_of(b).unwrap_or(0.0) {
                    self.best_hit_found = other.best_hit_found.clone();
                }
            }
            (None, Some(b)) => self.best_hit_found = Some(b.clone()),
            _ => {}
        }
        other.max_num_results = 0;
    }
}

impl Clone for SearchResults {
    fn clone(&self) -> Self {
        SearchResults {
            hit_molecules: self.hit_molecules.clone(),
            mol_names: HashMap::new(),
            max_num_results: self.max_num_results,
            timed_out: self.timed_out,
            cancelled: self.cancelled,
            best_hit_found: self.best_hit_found.clone(),
        }
    }
}
