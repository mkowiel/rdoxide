//! The maximum common substructure search (upstream
//! `FMCS/MaximumCommonSubgraph.{h,cpp}`): seeds are grown bond by bond in
//! the first (smallest) molecule and matched against all the others.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use std::time::Instant;

use super::cache::{DuplicatedSeedCache, SubstructureCache, seed_hash_key};
use super::graph::{Graph, NOT_SET};
use super::params::{
    AtomTyper, BondTyper, MCSParameters, MCSProgressData, MCSResult, RingComparator,
};
use super::seed::{MolFragment, Seed, TargetMatch};
use super::target::{FinalCheck, MatchTable, Target, substruct_match_custom_table};
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::Bond;
use crate::graphmol::queries::{
    CompositeQueryType, Query, QueryKind, make_atom_in_n_rings_query, make_atom_in_ring_query,
    make_atom_isotope_query, make_atom_min_ring_size_query, make_atom_num_query,
    make_bond_is_in_ring_query, make_bond_order_equals_query,
};
use crate::graphmol::query_ops::{atom_expand_query, bond_expand_query};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::graphmol::{BondStereo, BondType, ROMol};
use crate::{RDKitError, Result};

/// The current best result (`MaximumCommonSubgraph::MCS`).
#[derive(Debug, Clone, Default)]
struct Mcs<'a> {
    atoms: Vec<u32>,
    bonds: Vec<u32>,
    query_molecule: Option<&'a ROMol>,
    targets: Arc<Vec<Target<'a>>>,
}

/// `MaximumCommonSubgraph`
#[derive(Debug)]
pub struct MaximumCommonSubgraph<'a> {
    start: Instant,
    stat: MCSProgressData,
    parameters: MCSParameters,
    /// Whether matches are verified by `final_match_check` (disabled while
    /// the initial seeds are made).
    final_check_enabled: bool,
    /// minimum number of matches
    threshold_count: u32,
    molecules: Vec<&'a ROMol>,
    query_atom_labels: Vec<u32>,
    query_bond_labels: Vec<u32>,
    hash_cache: SubstructureCache,
    query_atom_match_table: MatchTable,
    query_bond_match_table: MatchTable,
    duplicate_cache: DuplicatedSeedCache,
    query_molecule: Option<&'a ROMol>,
    query_molecule_matched_bonds: usize,
    query_molecule_matched_atoms: u32,
    query_molecule_single_matched_atom: Option<u32>,
    targets: Arc<Vec<Target<'a>>>,
    /// sorted by the number of bonds, the biggest first
    seeds: Vec<Seed>,
    mcs_idx: Mcs<'a>,
    degenerate_mcs_map: BTreeMap<Vec<u32>, Mcs<'a>>,
}

struct LabelDefinition {
    item_index: usize,
    value: u32,
}

fn query_is_atom_in_ring(mol: &ROMol, atom: u32) -> Result<bool> {
    Ok(mol.ring_info().num_atom_rings(atom)? > 0)
}

impl<'a> MaximumCommonSubgraph<'a> {
    pub fn new(params: &MCSParameters) -> Self {
        let mut parameters = params.clone();
        if parameters.atom_compare_parameters.match_chiral_tag {
            parameters.bond_compare_parameters.match_stereo = true;
        }
        MaximumCommonSubgraph {
            start: Instant::now(),
            stat: MCSProgressData::default(),
            parameters,
            final_check_enabled: true,
            threshold_count: 0,
            molecules: Vec::new(),
            query_atom_labels: Vec::new(),
            query_bond_labels: Vec::new(),
            hash_cache: SubstructureCache::default(),
            query_atom_match_table: MatchTable::default(),
            query_bond_match_table: MatchTable::default(),
            duplicate_cache: DuplicatedSeedCache::default(),
            query_molecule: None,
            query_molecule_matched_bonds: 0,
            query_molecule_matched_atoms: 0,
            query_molecule_single_matched_atom: None,
            targets: Arc::new(Vec::new()),
            seeds: Vec::new(),
            mcs_idx: Mcs::default(),
            degenerate_mcs_map: BTreeMap::new(),
        }
    }

    pub fn parameters(&self) -> &MCSParameters {
        &self.parameters
    }

    pub fn parameters_mut(&mut self) -> &mut MCSParameters {
        &mut self.parameters
    }

    /// The molecule seeds are grown in.
    pub fn query_molecule(&self) -> &'a ROMol {
        self.query_molecule.expect("query molecule is set")
    }

    pub fn max_number_bonds(&self) -> usize {
        self.mcs_idx.bonds.len()
    }

    pub fn max_number_atoms(&self) -> usize {
        self.mcs_idx.atoms.len()
    }

    fn final_check(&self) -> Option<FinalCheck<'_>> {
        self.final_check_enabled.then_some(FinalCheck {
            params: &self.parameters,
        })
    }

    fn clear(&mut self) {
        self.targets = Arc::new(Vec::new());
        self.molecules.clear();
        self.start = Instant::now();
    }

    /// `init`: the labels of the query atoms and bonds and the match tables
    /// of the query against itself and against every target.
    fn init(&mut self, start_idx: usize) -> Result<()> {
        let query = self.molecules[start_idx];
        self.query_molecule = Some(query);
        self.query_atom_labels.clear();
        self.query_bond_labels.clear();
        self.query_atom_match_table.clear();
        self.query_bond_match_table.clear();
        self.duplicate_cache.clear();
        self.hash_cache.clear();

        let p = &self.parameters;
        let n_atoms = query.num_atoms() as usize;
        self.query_atom_match_table.resize(n_atoms, n_atoms);
        for aj in 0..n_atoms {
            for ai in 0..n_atoms {
                let m = p.atom_typer.compare(
                    &p.atom_compare_parameters,
                    query,
                    ai as u32,
                    query,
                    aj as u32,
                )?;
                self.query_atom_match_table.set(ai, aj, m);
            }
        }
        let n_bonds = query.num_bonds() as usize;
        self.query_bond_match_table.resize(n_bonds, n_bonds);
        for aj in 0..n_bonds {
            for ai in 0..n_bonds {
                let m = p.bond_typer.compare(
                    &p.bond_compare_parameters,
                    query,
                    ai as u32,
                    query,
                    aj as u32,
                )?;
                self.query_bond_match_table.set(ai, aj, m);
            }
        }

        let mut current_label_value = 1u32;
        let mut labels: Vec<LabelDefinition> = Vec::new();
        self.query_atom_labels = vec![NOT_SET; n_atoms];
        for ai in 0..n_atoms {
            let atom = query.atom(ai as u32);
            let label = match &p.atom_typer {
                AtomTyper::Any => 1,
                AtomTyper::Elements => {
                    atom.atomic_num()
                        | if p.atom_compare_parameters.match_valences {
                            atom.total_valence()? >> 8
                        } else {
                            0
                        }
                }
                AtomTyper::Isotopes => {
                    atom.atomic_num()
                        | (atom.isotope() >> 8)
                        | if p.atom_compare_parameters.match_valences {
                            atom.total_valence()? >> 16
                        } else {
                            0
                        }
                }
                // AnyHeavyAtom and user defined functors
                typer => {
                    let mut found = None;
                    for label in &labels {
                        if typer.compare(
                            &p.atom_compare_parameters,
                            query,
                            label.item_index as u32,
                            query,
                            ai as u32,
                        )? {
                            found = Some(label.value);
                            break;
                        }
                    }
                    match found {
                        Some(v) => v,
                        None => {
                            current_label_value += 1;
                            labels.push(LabelDefinition {
                                item_index: ai,
                                value: current_label_value,
                            });
                            current_label_value
                        }
                    }
                }
            };
            self.query_atom_labels[ai] = label;
        }

        labels.clear();
        current_label_value = 1;
        self.query_bond_labels = vec![NOT_SET; n_bonds];
        for aj in 0..n_bonds {
            let bond = query.bond(aj as u32);
            let label = match &p.bond_typer {
                BondTyper::Any => 1,
                BondTyper::OrderExact => bond.bond_type() as u32 + 1,
                BondTyper::Order => {
                    use BondType::*;
                    let order = match bond.bond_type() {
                        Aromatic | OneAndAHalf => Single,
                        TwoAndAHalf => Double,
                        ThreeAndAHalf => Triple,
                        FourAndAHalf => Quadruple,
                        FiveAndAHalf => Quintuple,
                        other => other,
                    };
                    order as u32 + 1
                }
                typer => {
                    let mut found = None;
                    for label in &labels {
                        if typer.compare(
                            &p.bond_compare_parameters,
                            query,
                            label.item_index as u32,
                            query,
                            aj as u32,
                        )? {
                            found = Some(label.value);
                            break;
                        }
                    }
                    match found {
                        Some(v) => v,
                        None => {
                            current_label_value += 1;
                            labels.push(LabelDefinition {
                                item_index: aj,
                                value: current_label_value,
                            });
                            current_label_value
                        }
                    }
                }
            };
            self.query_bond_labels[aj] = label;
        }

        // upstream starts the targets at the second molecule
        let mut targets = Vec::with_capacity(self.molecules.len().saturating_sub(1));
        for &mol in &self.molecules[1..] {
            let mut topology = Graph::new();
            for a in 0..mol.num_atoms() {
                topology.add_atom(a);
            }
            for b in mol.bonds() {
                topology.add_bond(b.idx(), b.begin_atom_idx(), b.end_atom_idx());
            }
            let (nq, nt) = (n_atoms, mol.num_atoms() as usize);
            let mut atom_match_table = MatchTable::default();
            atom_match_table.resize(nq, nt);
            for aj in 0..nt {
                for ai in 0..nq {
                    let m = p.atom_typer.compare(
                        &p.atom_compare_parameters,
                        query,
                        ai as u32,
                        mol,
                        aj as u32,
                    )?;
                    atom_match_table.set(ai, aj, m);
                }
            }
            let (nq, nt) = (n_bonds, mol.num_bonds() as usize);
            let mut bond_match_table = MatchTable::default();
            bond_match_table.resize(nq, nt);
            for aj in 0..nt {
                for ai in 0..nq {
                    let m = p.bond_typer.compare(
                        &p.bond_compare_parameters,
                        query,
                        ai as u32,
                        mol,
                        aj as u32,
                    )?;
                    bond_match_table.set(ai, aj, m);
                }
            }
            targets.push(Target {
                molecule: mol,
                topology,
                atom_match_table,
                bond_match_table,
            });
        }
        self.targets = Arc::new(targets);
        Ok(())
    }

    /// Adds `seed` to the seed set, keeping it sorted by the number of bonds
    /// (after the seeds with at least as many).
    fn add_seed(&mut self, seed: Seed) {
        let pos = self
            .seeds
            .iter()
            .position(|s| s.num_bonds() < seed.num_bonds())
            .unwrap_or(self.seeds.len());
        self.seeds.insert(pos, seed);
    }

    /// `makeInitialSeeds`
    fn make_initial_seeds(&mut self) -> Result<()> {
        let query = self.query_molecule();
        let mut excluded_bonds = vec![false; query.num_bonds() as usize];
        self.seeds.clear();
        self.query_molecule_matched_bonds = 0;
        self.query_molecule_matched_atoms = 0;
        self.query_molecule_single_matched_atom = None;
        let n_targets = self.targets.len();
        if !self.parameters.initial_seed.is_empty() {
            // user defined seed
            let initial_seed_molecule = mol_from_smarts(
                &self.parameters.initial_seed,
                &SmartsParserParams::default(),
            )?
            .ok_or_else(|| RDKitError::Value("invalid InitialSeed SMARTS".into()))?;
            let matching_substructs = substruct_match(
                query,
                &initial_seed_molecule,
                &SubstructMatchParameters::default(),
            )?;
            for ms in &matching_substructs {
                let mut seed = Seed::new();
                seed.set_store_all_degenerate_mcs(self.parameters.store_all);
                seed.excluded_bonds = excluded_bonds.clone();
                seed.match_result = vec![TargetMatch::default(); n_targets];
                let mut initial_seed_to_query_atom: BTreeMap<u32, u32> = BTreeMap::new();
                for &(sai, qai) in ms {
                    seed.add_atom(qai as u32);
                    initial_seed_to_query_atom.insert(sai as u32, qai as u32);
                }
                for &(sai, _) in ms {
                    for bond_idx in initial_seed_molecule.atom_bond_indices(sai as u32) {
                        let initial_bond = initial_seed_molecule.bond(bond_idx);
                        let qai1 = initial_seed_to_query_atom[&initial_bond.begin_atom_idx()];
                        let qai2 = initial_seed_to_query_atom[&initial_bond.end_atom_idx()];
                        let b = query
                            .bond_idx_between(qai1, qai2)
                            .ok_or_else(|| RDKitError::invariant("bond must not be NULL"))?;
                        if !seed.excluded_bonds[b as usize] {
                            seed.add_bond(query, b)?;
                            seed.excluded_bonds[b as usize] = true;
                        }
                    }
                }
                seed.compute_remaining_size(query);
                if self.check_if_match_and_append(&mut seed)? {
                    self.query_molecule_matched_bonds = seed.num_bonds();
                }
            }
            if self.seeds.is_empty() {
                crate::rd_log_warning!(
                    "The provided InitialSeed is not an MCS and will be ignored"
                );
            }
        }
        if self.seeds.is_empty() {
            // a seed from each query bond
            for bond in 0..query.num_bonds() {
                let b = query.bond(bond);
                let mut seed = Seed::new();
                seed.set_store_all_degenerate_mcs(self.parameters.store_all);
                seed.match_result = vec![TargetMatch::default(); n_targets];
                seed.add_atom(b.begin_atom_idx());
                seed.add_atom(b.end_atom_idx());
                seed.excluded_bonds = excluded_bonds.clone(); // all bonds from first to current
                seed.add_bond(query, bond)?;
                excluded_bonds[bond as usize] = true;
                seed.compute_remaining_size(query);
                if self.check_if_match_and_append(&mut seed)? {
                    self.query_molecule_matched_bonds += 1;
                } else {
                    for s in &mut self.seeds {
                        s.excluded_bonds[bond as usize] = true;
                    }
                }
            }
        }

        let nq = query.num_atoms();
        let targets = Arc::clone(&self.targets);
        for i in 0..nq {
            let query_mol_atom_in_ring = query_is_atom_in_ring(query, i)?;
            let mut matched = 0;
            let mut cand_single: Option<u32> = None;
            for tag in targets.iter() {
                let nt = tag.molecule.num_atoms();
                for aj in 0..nt {
                    if tag.atom_match_table.at(i as usize, aj as usize) {
                        let target_atom_in_ring = query_is_atom_in_ring(tag.molecule, aj)?;
                        matched += 1;
                        if !(self.parameters.bond_compare_parameters.complete_rings_only
                            && (query_mol_atom_in_ring || target_atom_in_ring))
                        {
                            let should_accept = match &self.parameters.should_accept_mcs {
                                None => true,
                                Some(accept) => accept(
                                    query,
                                    tag.molecule,
                                    &[(i as i32, aj as i32)],
                                    &[],
                                    &self.parameters,
                                )?,
                            };
                            if should_accept {
                                cand_single = Some(i);
                            }
                        }
                        break;
                    }
                }
            }
            if matched > 0 && matched >= self.threshold_count {
                self.query_molecule_matched_atoms += 1;
                if let Some(cand) = cand_single {
                    self.query_molecule_single_matched_atom =
                        Some(match self.query_molecule_single_matched_atom {
                            None => cand,
                            Some(cur) => {
                                // std::max(cand, cur, comp)
                                let (a, b) = (query.atom(cand), query.atom(cur));
                                let less = |a: &Atom, b: &Atom| {
                                    let (da, db) =
                                        (query.atom_degree(a.idx()), query.atom_degree(b.idx()));
                                    if da != db {
                                        da < db
                                    } else if a.formal_charge() != b.formal_charge() {
                                        a.formal_charge() < b.formal_charge()
                                    } else if a.atomic_num() != b.atomic_num() {
                                        a.atomic_num() < b.atomic_num()
                                    } else {
                                        a.idx() < b.idx()
                                    }
                                };
                                if less(a, b) { cur } else { cand }
                            }
                        });
                }
            }
        }
        Ok(())
    }

    /// `growSeeds`: returns false if canceled.
    fn grow_seeds(&mut self) -> Result<bool> {
        let mut mcs_found = false;
        let mut canceled = false;
        while !self.seeds.is_empty() {
            if self.max_number_bonds() == self.query_molecule_matched_bonds {
                break; // MCS == query
            }
            // The seed being grown is taken out of the set, so new seeds can
            // be added; they all go in front of it.
            let mut si = self.seeds.remove(0);
            let before = self.seeds.len();
            si.grow(self)?;
            let added = self.seeds.len() - before;
            self.seeds.insert(added, si);

            let fs = &self.seeds[0];
            let (fs_atoms, fs_bonds) = (fs.num_atoms(), fs.num_bonds());
            let (max_atoms, max_bonds) = (self.max_number_atoms(), self.max_number_bonds());
            let mut possible_mcs = if !self.parameters.maximize_bonds {
                fs_atoms > max_atoms || (fs_atoms == max_atoms && fs_bonds > max_bonds)
            } else {
                fs_bonds > max_bonds || (fs_bonds == max_bonds && fs_atoms > max_atoms)
            };
            let is_degenerate_mcs = fs_bonds == max_bonds && fs_atoms == max_atoms;
            if !possible_mcs && self.parameters.store_all {
                possible_mcs = is_degenerate_mcs;
            }
            if possible_mcs && self.parameters.bond_compare_parameters.complete_rings_only {
                possible_mcs = check_if_rings_are_closed(
                    self.query_molecule(),
                    &self.seeds[0].fragment,
                    self.parameters.atom_compare_parameters.complete_rings_only,
                )?;
            }
            if possible_mcs {
                possible_mcs = check_if_should_accept_mcs(
                    &self.seeds[0].fragment,
                    self.query_molecule(),
                    &self.targets,
                    &self.parameters,
                )?;
            }
            if possible_mcs {
                mcs_found = true;
                self.mcs_idx.atoms = self.seeds[0].fragment.atoms.clone();
                self.mcs_idx.bonds = self.seeds[0].fragment.bonds.clone();
                if self.parameters.store_all {
                    if !is_degenerate_mcs {
                        self.degenerate_mcs_map.clear();
                    }
                    let mut key = self.mcs_idx.bonds.clone();
                    key.sort_unstable();
                    let mut value = self.mcs_idx.clone();
                    value.query_molecule = self.query_molecule;
                    value.targets = Arc::clone(&self.targets);
                    self.degenerate_mcs_map.insert(key, value);
                }
            }
            if self.seeds[added].growing_stage == NOT_SET {
                self.seeds.remove(added); // finished
            }
            self.stat.num_atoms = self.max_number_atoms() as u32;
            self.stat.num_bonds = self.max_number_bonds() as u32;
            if !self.progress()? {
                canceled = true;
                break;
            }
        }
        if mcs_found {
            // postponed copy of the current set of molecules
            self.mcs_idx.query_molecule = self.query_molecule;
            self.mcs_idx.targets = Arc::clone(&self.targets);
        }
        Ok(!canceled)
    }

    /// The progress callback; false stops the search.
    fn progress(&self) -> Result<bool> {
        match &self.parameters.progress_callback {
            Some(cb) => cb(&self.stat, &self.parameters),
            None => {
                let elapsed = self.start.elapsed().as_micros() as u64;
                Ok(self.parameters.timeout == 0
                    || elapsed <= u64::from(self.parameters.timeout) * 1_000_000)
            }
        }
    }

    /// `createSeedFromMCS`: the MCS as a seed of target `new_query_target`.
    fn create_seed_from_mcs(&self, new_query_target: usize, new_seed: &mut Seed) -> Result<bool> {
        let mcs_query = self.mcs_idx.query_molecule.expect("MCS query molecule");
        let mut mcs = Seed::new();
        mcs.set_store_all_degenerate_mcs(self.parameters.store_all);
        mcs.excluded_bonds = vec![false; mcs_query.num_bonds() as usize];
        let mut mcs_atom_idx_map = vec![0u32; mcs_query.num_atoms() as usize];
        for &atom in &self.mcs_idx.atoms {
            mcs_atom_idx_map[atom as usize] = mcs.add_atom(atom);
        }
        for &bond in &self.mcs_idx.bonds {
            mcs.add_bond(mcs_query, bond)?;
        }
        let new_query = &self.mcs_idx.targets[new_query_target];
        let Some(matched) = substruct_match_custom_table(
            &new_query.topology,
            new_query.molecule,
            &mcs.topology,
            mcs_query,
            &new_query.atom_match_table,
            &new_query.bond_match_table,
            self.final_check(),
        )?
        else {
            return Ok(false);
        };
        let mut target_atoms = vec![NOT_SET; mcs.num_atoms()];
        new_seed.excluded_bonds = vec![false; new_query.molecule.num_bonds() as usize];
        for &(ai, tj) in &matched {
            target_atoms[ai as usize] = new_query.topology.atom(tj);
            new_seed.add_atom(new_query.topology.atom(tj));
        }
        for &bond in &self.mcs_idx.bonds {
            let b = mcs_query.bond(bond);
            let i = mcs_atom_idx_map[b.begin_atom_idx() as usize];
            let j = mcs_atom_idx_map[b.end_atom_idx() as usize];
            let tb = new_query
                .molecule
                .bond_idx_between(target_atoms[i as usize], target_atoms[j as usize])
                .ok_or_else(|| RDKitError::invariant("tb most not be NULL"))?;
            new_seed.add_bond(new_query.molecule, tb)?;
        }
        new_seed.compute_remaining_size(new_query.molecule);
        Ok(true)
    }

    /// `find`: the MCS of `mols` (all with initialized ring information).
    pub fn find(&mut self, src_mols: &[&'a ROMol]) -> Result<MCSResult> {
        self.clear();
        let mut res = MCSResult::default();
        if src_mols.len() < 2 {
            return Err(RDKitError::Runtime(
                "FMCS. Invalid argument. mols.size() must be at least 2".into(),
            ));
        }
        if self.parameters.threshold > 1.0 {
            return Err(RDKitError::Runtime(
                "FMCS. Invalid argument. Parameter Threshold must be 1.0 or less.".into(),
            ));
        }
        let n = src_mols.len() as i32;
        self.threshold_count = ((n - 1)
            .min(1.max((f64::from(n) * self.parameters.threshold).ceil() as i32 - 1)))
            as u32;

        if self.parameters.atom_compare_parameters.complete_rings_only {
            self.parameters.bond_compare_parameters.complete_rings_only = true;
        }
        if self.parameters.bond_compare_parameters.complete_rings_only {
            self.parameters
                .bond_compare_parameters
                .ring_matches_ring_only = true;
        }
        if self.parameters.atom_compare_parameters.complete_rings_only {
            self.parameters
                .atom_compare_parameters
                .ring_matches_ring_only = true;
        }

        self.molecules = src_mols.to_vec();
        self.molecules.sort_by_key(|m| m.num_bonds()); // stable
        let mut start_idx = 0usize;
        let end_idx = self.molecules.len() - self.threshold_count as usize;
        while start_idx < end_idx && self.molecules[start_idx].num_atoms() == 0 {
            start_idx += 1;
        }
        let mut are_seeds_empty = false;
        let mut i = start_idx;
        while i < end_idx && !are_seeds_empty && !res.canceled {
            self.init(start_idx)?;
            if self.targets.is_empty() {
                break;
            }
            self.final_check_enabled = false;
            let seeds = self.make_initial_seeds();
            self.final_check_enabled = true;
            seeds?;

            are_seeds_empty = self.seeds.is_empty();
            res.canceled = !(are_seeds_empty || self.grow_seeds()?);
            if self.max_number_bonds() == 0 {
                self.mcs_idx = Mcs::default();
                self.make_initial_seeds()?; // check all possible initial seeds
                if !are_seeds_empty {
                    let query = self.query_molecule();
                    let fs = &self.seeds[0];
                    let single_ring_bond = fs.fragment.bonds.len() == 1
                        && query.ring_info().num_bond_rings(fs.fragment.bonds[0])? > 0;
                    if (self.max_number_bonds() == 1
                        || !(self.parameters.bond_compare_parameters.complete_rings_only
                            && single_ring_bond))
                        && check_if_should_accept_mcs(
                            &fs.fragment,
                            query,
                            &self.targets,
                            &self.parameters,
                        )?
                    {
                        self.mcs_idx.query_molecule = Some(query);
                        self.mcs_idx.targets = Arc::clone(&self.targets);
                        self.mcs_idx.atoms = fs.fragment.atoms.clone();
                        self.mcs_idx.bonds = fs.fragment.bonds.clone();
                    }
                }
                if self.mcs_idx.query_molecule.is_none()
                    && let Some(atom) = self.query_molecule_single_matched_atom
                {
                    self.mcs_idx.query_molecule = self.query_molecule;
                    self.mcs_idx.targets = Arc::clone(&self.targets);
                    self.mcs_idx.atoms = vec![atom];
                    self.mcs_idx.bonds = Vec::new();
                }
            } else if i + 1 < end_idx {
                let mut seed = Seed::new();
                if self.create_seed_from_mcs(i, &mut seed)? {
                    // the MCS is matched with the new query
                    self.seeds.push(seed);
                }
                // change the query molecule for threshold < 1
                self.molecules.swap(start_idx, i + 1);
            }
            i += 1;
        }

        res.num_atoms = self.max_number_atoms() as u32;
        if res.num_atoms == 0 && self.query_molecule_single_matched_atom.is_some() {
            res.num_atoms = 1;
        }
        res.num_bonds = self.max_number_bonds() as u32;

        if res.num_bonds > 0 || self.query_molecule_single_matched_atom.is_some() {
            if !self.parameters.store_all {
                let (smarts, mol) = self.generate_result_smarts_and_query_mol(&self.mcs_idx)?;
                res.smarts_string = smarts;
                res.query_mol = Some(mol);
            } else {
                for mcs in self.degenerate_mcs_map.values() {
                    let (smarts, mol) = self.generate_result_smarts_and_query_mol(mcs)?;
                    res.degenerate_smarts_query_mol_dict
                        .entry(smarts)
                        .or_insert(mol);
                }
            }
        }
        self.clear();
        Ok(res)
    }

    /// `checkIfMatchAndAppend`: matches the seed against the targets and, if
    /// it matches, adds it to the seed set.
    pub(crate) fn check_if_match_and_append(&mut self, seed: &mut Seed) -> Result<bool> {
        let mut found_in_cache = false;
        let mut cache_key = 0u64;
        let mut cache_entry: Option<usize> = None;

        if let Some(value) = self.duplicate_cache.find(&seed.dup_cache_key) {
            if !value {
                return Ok(false); // mismatched
            }
            found_in_cache = true;
        }
        let found_in_dup_cache = found_in_cache;
        if !found_in_cache {
            cache_key = seed_hash_key(seed, &self.query_atom_labels, &self.query_bond_labels);
            cache_entry = self.hash_cache.find(cache_key);
            if let Some(entry) = cache_entry {
                // possibly found: check for a hash collision
                let query = self.query_molecule();
                for g in self.hash_cache.entry(entry) {
                    if g.num_vertices() != seed.num_atoms() || g.num_edges() != seed.num_bonds() {
                        continue;
                    }
                    found_in_cache = substruct_match_custom_table(
                        g,
                        query,
                        &seed.topology,
                        query,
                        &self.query_atom_match_table,
                        &self.query_bond_match_table,
                        self.final_check(),
                    )?
                    .is_some();
                    if !found_in_cache {
                        break;
                    }
                }
            }
        }
        let found = found_in_cache || self.match_seed(seed)?;

        if found {
            // store the new seed
            self.add_seed(seed.clone());
            if !found_in_dup_cache && seed.num_bonds() >= 3 {
                self.duplicate_cache.add(&seed.dup_cache_key, true);
            }
            if !found_in_cache {
                self.hash_cache.add(seed, cache_key, cache_entry);
            }
        } else if seed.num_bonds() > 3 {
            // cache mismatched duplicates too
            self.duplicate_cache.add(&seed.dup_cache_key, false);
        }
        Ok(found)
    }

    /// `match`: whether the seed is contained in enough targets.
    pub(crate) fn match_seed(&mut self, seed: &mut Seed) -> Result<bool> {
        let targets = Arc::clone(&self.targets);
        let query = self.query_molecule();
        let max_miss = targets.len() as u32 - self.threshold_count;
        let mut missing = 0;
        let mut passed = 0;
        for (itarget, tag) in targets.iter().enumerate() {
            let mut target_matched = false;
            if !seed.match_result.is_empty() && !seed.match_result[itarget].empty {
                target_matched = self.match_incremental_fast(seed, itarget)?;
            }
            if !target_matched {
                // slow full match
                let matched = substruct_match_custom_table(
                    &tag.topology,
                    tag.molecule,
                    &seed.topology,
                    query,
                    &tag.atom_match_table,
                    &tag.bond_match_table,
                    self.final_check(),
                )?;
                target_matched = matched.is_some();
                match matched {
                    Some(m) => {
                        if seed.match_result.is_empty() {
                            seed.match_result = vec![TargetMatch::default(); targets.len()];
                        }
                        let mut tm = std::mem::take(&mut seed.match_result[itarget]);
                        tm.init(seed, &m, query, tag.molecule);
                        seed.match_result[itarget] = tm;
                    }
                    None => {
                        if !seed.match_result.is_empty() {
                            seed.match_result[itarget].clear();
                        }
                    }
                }
            }
            if target_matched {
                passed += 1;
                if passed >= self.threshold_count {
                    break; // it's enough
                }
            } else {
                missing += 1;
                if missing > max_miss {
                    break;
                }
            }
        }
        Ok(missing <= max_miss)
    }

    /// `matchIncrementalFast`: extends the stored match of the parent seed
    /// by the bonds added since.
    fn match_incremental_fast(&self, seed: &mut Seed, itarget: usize) -> Result<bool> {
        let target = &self.targets[itarget];
        let tmol = target.molecule;
        let mut m = std::mem::take(&mut seed.match_result[itarget]);
        let result = self.extend_match(seed, &mut m, target, tmol);
        seed.match_result[itarget] = m;
        result
    }

    fn extend_match(
        &self,
        seed: &Seed,
        m: &mut TargetMatch,
        target: &Target<'a>,
        tmol: &ROMol,
    ) -> Result<bool> {
        if m.empty {
            return Ok(false);
        }
        let mut matched = false;
        for new_bond_seed_idx in m.matched_bond_size..seed.num_bonds() {
            matched = false;
            let mut atom_added = false;
            let new_bond_query_idx = seed.fragment.bonds[new_bond_seed_idx];
            let (begin, end) = seed.fragment.bond_ends[new_bond_seed_idx];
            let i = seed.fragment.seed_atom_idx_map[&begin] as usize;
            let j = seed.fragment.seed_atom_idx_map[&end] as usize;
            let (source_seed_idx, other_seed_idx) = if i >= m.matched_atom_size {
                (j, i)
            } else {
                (i, j)
            };
            let other_query_idx = seed.fragment.atoms[other_seed_idx];
            let source_query_idx = seed.fragment.atoms[source_seed_idx];
            let source_target_idx = m.target_atom_idx[source_query_idx as usize];
            let mut tb: Option<u32> = None;
            let mut other_target_idx = NOT_SET;

            if other_seed_idx < m.matched_atom_size {
                other_target_idx = m.target_atom_idx[other_query_idx as usize];
                tb = tmol.bond_idx_between(source_target_idx, other_target_idx);
                if let Some(tbi) = tb
                    && !m.visited_target_bonds[tbi as usize]
                {
                    matched = target
                        .bond_match_table
                        .at(new_bond_query_idx as usize, tbi as usize);
                }
            } else {
                for (nbr, tbi) in tmol
                    .atom_neighbors(source_target_idx)
                    .zip(tmol.atom_bond_indices(source_target_idx))
                {
                    tb = Some(tbi);
                    if m.visited_target_bonds[tbi as usize] {
                        continue;
                    }
                    other_target_idx = nbr;
                    if m.visited_target_atoms[other_target_idx as usize] {
                        continue;
                    }
                    matched = target
                        .atom_match_table
                        .at(other_query_idx as usize, other_target_idx as usize)
                        && target
                            .bond_match_table
                            .at(new_bond_query_idx as usize, tbi as usize);
                    if matched {
                        atom_added = true;
                        break;
                    }
                }
            }

            if matched {
                // update the match history
                let tbi = tb.expect("matched bond");
                if atom_added {
                    m.matched_atom_size += 1;
                    m.target_atom_idx[other_query_idx as usize] = other_target_idx;
                    m.visited_target_atoms[other_target_idx as usize] = true;
                }
                m.matched_bond_size += 1;
                m.target_bond_idx[new_bond_query_idx as usize] = tbi;
                m.visited_target_bonds[tbi as usize] = true;
            } else {
                m.clear();
                return Ok(false);
            }
        }

        if m.matched_atom_size != seed.num_atoms() || m.matched_bond_size != seed.num_bonds() {
            // number of unique items
            m.clear();
            return Ok(false);
        }
        if matched && self.final_check_enabled {
            let c1: Vec<u32> = (0..seed.num_atoms() as u32).collect();
            let c2: Vec<u32> = (0..seed.num_atoms() as u32)
                .map(|si| m.target_atom_idx[seed.topology.atom(si) as usize])
                .collect();
            matched = super::final_match::final_match_check(
                &c1,
                &c2,
                self.query_molecule(),
                &seed.topology,
                tmol,
                &target.topology,
                &self.parameters,
            )?; // check CHIRALITY
            if !matched {
                m.clear();
            }
        }
        Ok(matched)
    }

    /// `generateResultSMARTSAndQueryMol`
    fn generate_result_smarts_and_query_mol(&self, mcs: &Mcs<'a>) -> Result<(String, Arc<ROMol>)> {
        let mcs_query = mcs.query_molecule.expect("MCS query molecule");
        let query = self.query_molecule();
        let mut seed = Seed::new();
        seed.set_store_all_degenerate_mcs(self.parameters.store_all);
        seed.excluded_bonds = vec![false; mcs_query.num_bonds() as usize];
        let mut atom_idx_map = vec![0u32; mcs_query.num_atoms() as usize];
        // key is a unique bond type / atomic number
        let mut bond_match_set: Vec<BTreeSet<BondType>> = vec![BTreeSet::new(); mcs.bonds.len()];
        let mut bond_stereo_set: Vec<BTreeMap<BondType, BondStereo>> =
            vec![BTreeMap::new(); mcs.bonds.len()];
        let mut atom_match_set: Vec<BTreeMap<u32, (u32, ChiralType)>> =
            vec![BTreeMap::new(); mcs.atoms.len()];
        for &atom in &mcs.atoms {
            atom_idx_map[atom as usize] = seed.num_atoms() as u32;
            seed.add_atom(atom);
        }
        for &bond in &mcs.bonds {
            seed.add_bond(mcs_query, bond)?;
        }

        if !mcs.bonds.is_empty() {
            for tag in mcs.targets.iter() {
                // there is no bond match info
                let Some(matched) = substruct_match_custom_table(
                    &tag.topology,
                    tag.molecule,
                    &seed.topology,
                    query,
                    &tag.atom_match_table,
                    &tag.bond_match_table,
                    self.final_check(),
                )?
                else {
                    continue;
                };
                let mut target_atoms = vec![NOT_SET; seed.num_atoms()];
                for &(ai, tj) in &matched {
                    let target_atom = tag.topology.atom(tj);
                    target_atoms[ai as usize] = target_atom;
                    let ta = tag.molecule.atom(target_atom);
                    if ta.atomic_num()
                        != mcs_query
                            .atom(seed.fragment.atoms[ai as usize])
                            .atomic_num()
                    {
                        atom_match_set[ai as usize]
                            .insert(ta.atomic_num(), (ta.atomic_num(), ta.chiral_tag()));
                    }
                }
                for (bi, &bond) in mcs.bonds.iter().enumerate() {
                    let b = mcs_query.bond(bond);
                    let i = atom_idx_map[b.begin_atom_idx() as usize];
                    let j = atom_idx_map[b.end_atom_idx() as usize];
                    let (ti, tj) = (target_atoms[i as usize], target_atoms[j as usize]);
                    if let Some(tbi) = tag.molecule.bond_idx_between(ti, tj) {
                        let tb = tag.molecule.bond(tbi);
                        if b.bond_type() != tb.bond_type() {
                            bond_match_set[bi].insert(tb.bond_type());
                            bond_stereo_set[bi].insert(tb.bond_type(), tb.stereo());
                        }
                    }
                }
            }
        }

        let mut mol = ROMol::new();
        let ri = mcs_query.ring_info();
        let need_atom_ring_queries = self
            .parameters
            .atom_compare_parameters
            .ring_matches_ring_only
            || self
                .parameters
                .bond_compare_parameters
                .match_fused_rings_strict;
        let mut mcs_ring_is_complete: Vec<bool> = Vec::new();
        if need_atom_ring_queries {
            mcs_ring_is_complete = vec![true; ri.num_rings()? as usize];
            let mut query_bond_in_mcs = vec![false; mcs_query.num_bonds() as usize];
            for &bond in &mcs.bonds {
                query_bond_in_mcs[bond as usize] = true;
            }
            for (ring_idx, bond_ring) in ri.bond_rings().iter().enumerate() {
                if bond_ring.iter().any(|&b| !query_bond_in_mcs[b as usize]) {
                    mcs_ring_is_complete[ring_idx] = false;
                }
            }
        }
        for (ai, &atom_idx) in mcs.atoms.iter().enumerate() {
            let atom = mcs_query.atom(atom_idx);
            let num_atom_rings = ri.num_atom_rings(atom_idx)?;
            let mut a = Atom::new(0);
            if matches!(self.parameters.atom_typer, AtomTyper::Isotopes)
                || self.parameters.atom_compare_parameters.match_isotope
            {
                // '[0*]-[0*]-[13*]' for CC[13NH2]
                a.set_query(Some(make_atom_isotope_query(atom.isotope() as i32)));
            } else {
                a.set_query(Some(make_atom_num_query(atom.atomic_num() as i32)));
                for &(atomic_num, chiral_tag) in atom_match_set[ai].values() {
                    atom_expand_query(
                        &mut a,
                        make_atom_num_query(atomic_num as i32),
                        CompositeQueryType::Or,
                        true,
                    )?;
                    if self.parameters.atom_compare_parameters.match_chiral_tag
                        && matches!(
                            chiral_tag,
                            ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
                        )
                    {
                        a.set_chiral_tag(chiral_tag);
                    }
                }
            }
            if need_atom_ring_queries {
                let rings_of_atom = ri.atom_members(atom_idx)?;
                let num_complete_rings = rings_of_atom
                    .iter()
                    .filter(|&&r| mcs_ring_is_complete[r as usize])
                    .count();
                if self
                    .parameters
                    .atom_compare_parameters
                    .ring_matches_ring_only
                    && num_complete_rings == 0
                {
                    let mut q = make_atom_in_ring_query();
                    q.negate = num_atom_rings == 0;
                    atom_expand_query(&mut a, q, CompositeQueryType::And, true)?;
                } else if self
                    .parameters
                    .bond_compare_parameters
                    .match_fused_rings_strict
                    && num_atom_rings == 1
                    && num_complete_rings == 1
                {
                    let ring_size = ri.atom_rings()[rings_of_atom[0] as usize].len();
                    let mut q = Query::new(QueryKind::Or);
                    q.description = "AtomOr".to_string();
                    q.children
                        .push(make_atom_min_ring_size_query(ring_size as i32));
                    let mut q2 = make_atom_in_n_rings_query(1);
                    q2.negate = true;
                    q.children.push(q2);
                    atom_expand_query(&mut a, q, CompositeQueryType::And, true)?;
                }
            }
            mol.add_atom(a, true);
        }
        for (bi, &bond_idx) in mcs.bonds.iter().enumerate() {
            let bond = mcs_query.bond(bond_idx);
            let mut b = Bond::new(BondType::Unspecified);
            b.set_begin_atom_idx_unchecked(atom_idx_map[bond.begin_atom_idx() as usize]);
            b.set_end_atom_idx_unchecked(atom_idx_map[bond.end_atom_idx() as usize]);
            b.set_query(Some(make_bond_order_equals_query(bond.bond_type())));
            for &bond_type in &bond_match_set[bi] {
                bond_expand_query(
                    &mut b,
                    make_bond_order_equals_query(bond_type),
                    CompositeQueryType::Or,
                    true,
                )?;
                if self.parameters.bond_compare_parameters.match_stereo
                    && bond_stereo_set[bi][&bond_type] > BondStereo::StereoAny
                {
                    b.set_stereo(bond_stereo_set[bi][&bond_type])?;
                }
            }
            if self
                .parameters
                .bond_compare_parameters
                .ring_matches_ring_only
                || self
                    .parameters
                    .bond_compare_parameters
                    .match_fused_rings_strict
            {
                let num_bond_rings = ri.num_bond_rings(bond_idx)?;
                let mut q = make_bond_is_in_ring_query();
                q.negate = num_bond_rings == 0;
                bond_expand_query(&mut b, q, CompositeQueryType::And, true)?;
            }
            mol.add_bond_object(b)?;
        }
        let smarts = mol_to_smarts(&mut mol, &SmilesWriteParams::default())?;
        Ok((smarts, Arc::new(mol)))
    }
}

/// `checkIfRingsAreClosed`
fn check_if_rings_are_closed(
    om: &ROMol,
    fs: &MolFragment,
    no_lone_ring_atoms: bool,
) -> Result<bool> {
    if fs.bonds.is_empty() && fs.atoms.is_empty() {
        return Ok(true);
    }
    let ri = om.ring_info();
    let n_rings = ri.num_rings()? as usize;
    if n_rings == 0 {
        return Ok(true);
    }
    let mut mcs_bonds = vec![false; om.num_bonds() as usize];
    let mut mcs_non_fused_rings = vec![false; n_rings];
    let mut mcs_fused_rings = vec![false; n_rings];
    for &bi in &fs.bonds {
        mcs_bonds[bi as usize] = true;
        if ri.num_bond_rings(bi)? == 1 {
            mcs_non_fused_rings[ri.bond_members(bi)?[0] as usize] = true;
        }
    }
    for ring_idx in 0..n_rings {
        if !mcs_non_fused_rings[ring_idx] {
            continue;
        }
        for &bi in &ri.bond_rings()[ring_idx] {
            let mut keep_bond = false;
            for &member_of in ri.bond_members(bi as u32)? {
                if member_of as usize == ring_idx {
                    keep_bond = true;
                } else if mcs_non_fused_rings[member_of as usize] {
                    keep_bond = false;
                    break;
                }
            }
            if keep_bond && !mcs_bonds[bi as usize] {
                return Ok(false);
            }
        }
    }
    if no_lone_ring_atoms {
        for &ai in &fs.atoms {
            let ring_indices = ri.atom_members(ai)?;
            if !ring_indices.is_empty()
                && !ring_indices
                    .iter()
                    .any(|&r| mcs_non_fused_rings[r as usize])
            {
                return Ok(false);
            }
        }
    }
    if !mcs_non_fused_rings.iter().any(|&b| b) {
        for &bi in &fs.bonds {
            if ri.num_bond_rings(bi)? > 1 {
                for &ring_idx in ri.bond_members(bi)? {
                    mcs_fused_rings[ring_idx as usize] = true;
                }
            }
        }
    }
    if mcs_fused_rings.iter().any(|&b| b) {
        for ring_idx in 0..n_rings {
            if !mcs_fused_rings[ring_idx] {
                continue;
            }
            if ri.bond_rings()[ring_idx]
                .iter()
                .all(|&bi| mcs_bonds[bi as usize])
            {
                return Ok(true);
            }
        }
        return Ok(false);
    }
    Ok(true)
}

/// `checkIfShouldAcceptMCS`
fn check_if_should_accept_mcs(
    f: &MolFragment,
    query: &ROMol,
    targets: &[Target<'_>],
    p: &MCSParameters,
) -> Result<bool> {
    let Some(accept) = &p.should_accept_mcs else {
        return Ok(true);
    };
    if f.bonds.is_empty() {
        return Ok(true);
    }
    let mut seed = Seed::new(); // result MCS
    seed.excluded_bonds = vec![false; query.num_bonds() as usize];
    for &atom in &f.atoms {
        seed.add_atom(atom);
    }
    for &bond in &f.bonds {
        seed.add_bond(query, bond)?;
    }
    for target in targets {
        let Some(matched) = substruct_match_custom_table(
            &target.topology,
            target.molecule,
            &seed.topology,
            query,
            &target.atom_match_table,
            &target.bond_match_table,
            Some(FinalCheck { params: p }),
        )?
        else {
            continue;
        };
        let mut atom_pair_match: Vec<(i32, i32)> = Vec::with_capacity(matched.len());
        let mut query_to_target_atom_idx = vec![NOT_SET; query.num_atoms() as usize];
        for &(m1, m2) in &matched {
            let query_atom_idx = seed.topology.atom(m1);
            let target_atom_idx = target.topology.atom(m2);
            atom_pair_match.push((query_atom_idx as i32, target_atom_idx as i32));
            query_to_target_atom_idx[query_atom_idx as usize] = target_atom_idx;
        }
        let mut bond_pair_match: Vec<(i32, i32)> = Vec::with_capacity(seed.topology.num_edges());
        for &(s, t, _) in seed.topology.edges() {
            let query_begin = seed.topology.atom(s);
            let query_end = seed.topology.atom(t);
            let query_bond = query
                .bond_idx_between(query_begin, query_end)
                .ok_or_else(|| RDKitError::invariant("query bond not found"))?;
            let target_begin = query_to_target_atom_idx[query_begin as usize];
            let target_end = query_to_target_atom_idx[query_end as usize];
            if target_begin == NOT_SET || target_end == NOT_SET {
                return Err(RDKitError::invariant("unmatched MCS atom"));
            }
            let target_bond = target
                .molecule
                .bond_idx_between(target_begin, target_end)
                .ok_or_else(|| RDKitError::invariant("target bond not found"))?;
            bond_pair_match.push((query_bond as i32, target_bond as i32));
        }
        if !accept(
            query,
            target.molecule,
            &atom_pair_match,
            &bond_pair_match,
            p,
        )? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// `findMCS(mols, params)`. Molecules without ring information are treated
/// as having no rings, like upstream.
pub fn find_mcs(mols: &[&ROMol], params: Option<&MCSParameters>) -> Result<MCSResult> {
    use crate::graphmol::ring_info::FindRingType;
    let default_params = MCSParameters::default();
    let params = params.unwrap_or(&default_params);
    // molecules without ring info get an (empty) initialized one on a copy
    let owned: Vec<Option<ROMol>> = mols
        .iter()
        .map(|m| {
            (!m.ring_info().is_initialized()).then(|| {
                let mut copy = ROMol::copy_from(m, false, -1);
                copy.ring_info_mut()
                    .initialize(FindRingType::OtherOrUnknown);
                copy
            })
        })
        .collect();
    let refs: Vec<&ROMol> = mols
        .iter()
        .zip(&owned)
        .map(|(&m, o)| o.as_ref().unwrap_or(m))
        .collect();
    let mut fmcs = MaximumCommonSubgraph::new(params);
    fmcs.find(&refs)
}

/// `findMCS(mols, maximizeBonds, threshold, timeout, ...)`
#[allow(clippy::too_many_arguments)]
pub fn find_mcs_with_options(
    mols: &[&ROMol],
    maximize_bonds: bool,
    threshold: f64,
    timeout: u32,
    verbose: bool,
    match_valences: bool,
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
    match_chiral_tag: bool,
    atom_comp: super::params::AtomComparator,
    bond_comp: super::params::BondComparator,
    ring_comp: RingComparator,
) -> Result<MCSResult> {
    let mut ps = MCSParameters {
        maximize_bonds,
        threshold,
        timeout,
        verbose,
        ..Default::default()
    };
    ps.set_atom_typer_from_enum(atom_comp);
    ps.atom_compare_parameters.match_valences = match_valences;
    ps.atom_compare_parameters.match_chiral_tag = match_chiral_tag;
    ps.atom_compare_parameters.ring_matches_ring_only = ring_matches_ring_only;
    ps.set_bond_typer_from_enum(bond_comp);
    ps.bond_compare_parameters.ring_matches_ring_only = ring_matches_ring_only;
    ps.bond_compare_parameters.complete_rings_only = complete_rings_only;
    ps.bond_compare_parameters.match_fused_rings = ring_comp != RingComparator::IgnoreRingFusion;
    ps.bond_compare_parameters.match_fused_rings_strict =
        ring_comp == RingComparator::StrictRingFusion;
    find_mcs(mols, Some(&ps))
}

/// `findMCS_P`: parameters as JSON.
pub fn find_mcs_json(mols: &[&ROMol], params_json: &str) -> Result<MCSResult> {
    let mut p = MCSParameters::default();
    super::params::parse_mcs_parameters_json(params_json, &mut p)?;
    find_mcs(mols, Some(&p))
}
