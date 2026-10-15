//! Seeds: connected fragments of the query molecule that are grown bond by
//! bond (upstream `FMCS/Seed.{h,cpp}`, `FMCS/TargetMatch.{h,cpp}`).

use std::collections::{BTreeMap, BTreeSet};

use super::cache::{Composition2N, DupCacheKey};
use super::graph::{Graph, NOT_SET};
use super::mcs::MaximumCommonSubgraph;
use crate::graphmol::ROMol;
use crate::{RDKitError, Result};

/// A match between a seed and a target: pairs of (seed vertex, target vertex).
pub type MatchV = Vec<(u32, u32)>;

/// `MolFragment`: a fragment of a molecule (all indices refer to it).
#[derive(Debug, Clone, Default)]
pub struct MolFragment {
    pub atoms: Vec<u32>,
    pub bonds: Vec<u32>,
    /// (begin atom, end atom) of each bond in `bonds`
    pub bond_ends: Vec<(u32, u32)>,
    /// molecule atom index -> position in `atoms`
    pub seed_atom_idx_map: BTreeMap<u32, u32>,
}

/// `NewBond`: an outgoing bond scheduled to be added to a seed.
#[derive(Debug, Clone)]
pub struct NewBond {
    pub bond_idx: u32,
    /// the atom at the far end of the bond
    pub new_atom_idx: u32,
    /// position in the seed of the far atom when it is already there (ring
    /// closure), `NOT_SET` otherwise
    pub end_atom_idx: u32,
}

/// `TargetMatch`: the match of a seed against one target, kept to extend it
/// incrementally when the seed grows.
#[derive(Debug, Clone)]
pub struct TargetMatch {
    pub empty: bool,
    pub matched_atom_size: usize,
    pub matched_bond_size: usize,
    pub target_atom_idx: Vec<u32>,
    pub target_bond_idx: Vec<u32>,
    pub visited_target_bonds: Vec<bool>,
    pub visited_target_atoms: Vec<bool>,
}

impl Default for TargetMatch {
    fn default() -> Self {
        TargetMatch {
            empty: true,
            matched_atom_size: 0,
            matched_bond_size: 0,
            target_atom_idx: Vec::new(),
            target_bond_idx: Vec::new(),
            visited_target_bonds: Vec::new(),
            visited_target_atoms: Vec::new(),
        }
    }
}

impl TargetMatch {
    pub fn clear(&mut self) {
        self.empty = true;
        self.target_atom_idx.clear();
        self.target_bond_idx.clear();
        self.visited_target_bonds.clear();
        self.visited_target_atoms.clear();
    }

    /// `TargetMatch::init`
    pub fn init(&mut self, seed: &Seed, matched: &MatchV, query: &ROMol, target: &ROMol) {
        self.target_atom_idx = vec![NOT_SET; query.num_atoms() as usize];
        self.target_bond_idx = vec![NOT_SET; query.num_bonds() as usize];
        self.visited_target_bonds = vec![false; target.num_bonds() as usize];
        self.visited_target_atoms = vec![false; target.num_atoms() as usize];
        self.matched_atom_size = matched.len();
        for &(seed_atom, target_atom) in matched {
            self.target_atom_idx[seed.fragment.atoms[seed_atom as usize] as usize] = target_atom;
            self.visited_target_atoms[target_atom as usize] = true;
        }
        self.matched_bond_size = 0;
        for (&bond, &(i, j)) in seed.fragment.bonds.iter().zip(&seed.fragment.bond_ends) {
            let ti = self.target_atom_idx[i as usize];
            let tj = self.target_atom_idx[j as usize];
            if let Some(tb) = target.bond_idx_between(ti, tj) {
                self.matched_bond_size += 1;
                self.target_bond_idx[bond as usize] = tb;
                self.visited_target_bonds[tb as usize] = true;
            }
        }
        self.empty = false;
    }
}

/// `Seed`
#[derive(Debug, Clone)]
pub struct Seed {
    /// all directly connected outgoing bonds, for multistage growing
    new_bonds: Vec<NewBond>,
    store_all_degenerate_mcs: bool,
    /// 0 new seed; `NOT_SET` finished; n > 0 in progress
    pub growing_stage: u32,
    pub fragment: MolFragment,
    /// topology with references to the source molecule
    pub topology: Graph,
    pub excluded_bonds: Vec<bool>,
    /// in this subgraph, for the performance of future growing
    pub last_added_atoms_begin_idx: usize,
    pub remaining_bonds: u32,
    pub remaining_atoms: u32,
    pub dup_cache_key: DupCacheKey,
    /// for each target
    pub match_result: Vec<TargetMatch>,
}

impl Default for Seed {
    fn default() -> Self {
        Seed {
            new_bonds: Vec::new(),
            store_all_degenerate_mcs: false,
            growing_stage: 0,
            fragment: MolFragment::default(),
            topology: Graph::new(),
            excluded_bonds: Vec::new(),
            last_added_atoms_begin_idx: 0,
            remaining_bonds: 0,
            remaining_atoms: 0,
            dup_cache_key: DupCacheKey::default(),
            match_result: Vec::new(),
        }
    }
}

impl Seed {
    pub fn new() -> Self {
        Seed::default()
    }

    pub fn num_atoms(&self) -> usize {
        self.fragment.atoms.len()
    }

    pub fn num_bonds(&self) -> usize {
        self.fragment.bonds.len()
    }

    pub fn set_store_all_degenerate_mcs(&mut self, value: bool) {
        self.store_all_degenerate_mcs = value;
    }

    /// `createFromParent`: the part of the parent that a child inherits.
    pub fn create_from_parent(&self) -> Seed {
        Seed {
            new_bonds: Vec::new(),
            store_all_degenerate_mcs: self.store_all_degenerate_mcs,
            growing_stage: 0,
            fragment: self.fragment.clone(),
            topology: self.topology.clone(),
            excluded_bonds: self.excluded_bonds.clone(),
            last_added_atoms_begin_idx: self.num_atoms(),
            remaining_bonds: self.remaining_bonds,
            remaining_atoms: self.remaining_atoms,
            dup_cache_key: self.dup_cache_key.clone(),
            match_result: Vec::new(),
        }
    }

    pub fn can_grow_bigger_than(&self, max_bonds: usize, max_atoms: usize) -> bool {
        let bonds = self.remaining_bonds as usize + self.num_bonds();
        let atoms = self.remaining_atoms as usize + self.num_atoms();
        bonds > max_bonds
            || (bonds == max_bonds
                && (atoms > max_atoms || (self.store_all_degenerate_mcs && atoms == max_atoms)))
    }

    /// `addAtom`: adds a molecule atom, returns its position in the seed.
    pub fn add_atom(&mut self, atom: u32) -> u32 {
        let i = self.fragment.atoms.len() as u32;
        self.fragment.atoms.push(atom);
        self.fragment.seed_atom_idx_map.insert(atom, i);
        self.topology.add_atom(atom);
        self.dup_cache_key.add_atom(atom);
        i
    }

    /// `addBond`: adds a molecule bond whose atoms are already in the seed.
    pub fn add_bond(&mut self, mol: &ROMol, bond: u32) -> Result<usize> {
        if self.excluded_bonds[bond as usize] {
            return Err(RDKitError::invariant("bond is already excluded"));
        }
        self.excluded_bonds[bond as usize] = true;
        let b = mol.bond(bond);
        let (begin, end) = (b.begin_atom_idx(), b.end_atom_idx());
        self.fragment.bonds.push(bond);
        self.fragment.bond_ends.push((begin, end));
        let i = self.fragment.seed_atom_idx_map[&begin];
        let j = self.fragment.seed_atom_idx_map[&end];
        self.topology.add_bond(bond, i, j);
        self.dup_cache_key.add_bond(bond);
        Ok(self.num_bonds())
    }

    /// `computeRemainingSize`: counts what could still be added to the seed.
    pub fn compute_remaining_size(&mut self, qmol: &ROMol) {
        self.remaining_bonds = 0;
        self.remaining_atoms = 0;
        let mut end_atom_stack: Vec<u32> = Vec::new();
        let mut visited_bonds = self.excluded_bonds.clone();
        let mut visited_atoms = vec![false; qmol.num_atoms() as usize];
        for &atom in &self.fragment.atoms {
            visited_atoms[atom as usize] = true;
        }
        // 1. direct neighbours of the atoms added last
        for seed_atom_idx in self.last_added_atoms_begin_idx..self.num_atoms() {
            let atom = self.fragment.atoms[seed_atom_idx];
            for (other, bond) in adjacent(qmol, atom) {
                if !visited_bonds[bond as usize] {
                    self.remaining_bonds += 1;
                    visited_bonds[bond as usize] = true;
                    if !visited_atoms[other as usize] {
                        self.remaining_atoms += 1;
                        visited_atoms[other as usize] = true;
                        end_atom_stack.push(other);
                    }
                }
            }
        }
        // 2. go deep
        while let Some(ai) = end_atom_stack.pop() {
            for (other, bond) in adjacent(qmol, ai) {
                if !visited_bonds[bond as usize] {
                    self.remaining_bonds += 1;
                    visited_bonds[bond as usize] = true;
                    if !visited_atoms[other as usize] {
                        self.remaining_atoms += 1;
                        visited_atoms[other as usize] = true;
                        end_atom_stack.push(other);
                    }
                }
            }
        }
    }

    fn add_new_bond_from_atom(&mut self, qmol: &ROMol, src_atom: u32, bond: u32) {
        let end_atom = qmol.bond(bond).other_atom_idx(src_atom).expect("bond atom");
        let end_atom_idx = self
            .fragment
            .atoms
            .iter()
            .position(|&a| a == end_atom)
            .map_or(NOT_SET, |p| p as u32);
        self.new_bonds.push(NewBond {
            bond_idx: bond,
            new_atom_idx: end_atom,
            end_atom_idx,
        });
    }

    /// `fillNewBonds`
    fn fill_new_bonds(
        &mut self,
        qmol: &ROMol,
        mut mcs: Option<&mut MaximumCommonSubgraph<'_>>,
    ) -> Result<()> {
        let mut excluded_bonds = self.excluded_bonds.clone();
        let ri = qmol.ring_info();
        // all atoms added on previous growing only
        for src_atom_idx in self.last_added_atoms_begin_idx..self.num_atoms() {
            let atom = self.fragment.atoms[src_atom_idx];
            for (_, bi) in adjacent(qmol, atom) {
                // already in the seed or in NewBonds from another atom in a RING
                if excluded_bonds[bi as usize] {
                    continue;
                }
                excluded_bonds[bi as usize] = true;
                if let Some(mcs) = mcs.as_deref_mut() {
                    if mcs.parameters().bond_compare_parameters.complete_rings_only
                        && ri.num_bond_rings(bi)? == 1
                        && !self
                            .can_add_all_non_fused_ring_bonds_connected_to_bond(atom, bi, mcs)?
                    {
                        continue;
                    }
                }
                self.add_new_bond_from_atom(qmol, atom, bi);
            }
        }
        Ok(())
    }

    /// `canAddAllNonFusedRingBondsConnectedToBond`
    fn can_add_all_non_fused_ring_bonds_connected_to_bond(
        &self,
        src_atom: u32,
        bond: u32,
        mcs: &mut MaximumCommonSubgraph<'_>,
    ) -> Result<bool> {
        let mol = mcs.query_molecule();
        let ri = mol.ring_info();
        let ring = ri.bond_members(bond)?[0] as usize;
        let bond_ring = &ri.bond_rings()[ring];
        let mut non_fused_ring_bond_indices: BTreeSet<u32> = BTreeSet::new();
        let mut connected_atom_indices = vec![false; mol.num_atoms() as usize];
        let mut seed = self.create_from_parent();
        for &bi in bond_ring {
            if bi as u32 != bond && ri.num_bond_rings(bi as u32)? == 1 {
                non_fused_ring_bond_indices.insert(bi as u32);
            }
        }
        let mut curr_atom = src_atom;
        let mut curr_bond = Some(bond);
        let excluded_bonds = seed.excluded_bonds.clone();
        while let Some(cb) = curr_bond {
            if !seed.excluded_bonds[cb as usize] {
                let b = mol.bond(cb);
                connected_atom_indices[b.begin_atom_idx() as usize] = true;
                connected_atom_indices[b.end_atom_idx() as usize] = true;
                seed.add_new_bond_from_atom(mol, curr_atom, cb);
            }
            curr_bond = None;
            for &cand in &non_fused_ring_bond_indices {
                let cand_bond = mol.bond(cand);
                if connected_atom_indices[cand_bond.begin_atom_idx() as usize] {
                    curr_atom = cand_bond.begin_atom_idx();
                } else if connected_atom_indices[cand_bond.end_atom_idx() as usize] {
                    curr_atom = cand_bond.end_atom_idx();
                } else {
                    continue;
                }
                non_fused_ring_bond_indices.remove(&cand);
                curr_bond = Some(cand);
                break;
            }
        }
        if seed.new_bonds.is_empty() {
            return Ok(false);
        }
        let new_bonds = std::mem::take(&mut seed.new_bonds);
        add_new_bonds_to_seed(&new_bonds, mol, &mut seed)?;
        seed.new_bonds = new_bonds;
        seed.match_result = self.match_result.clone();
        seed.excluded_bonds = excluded_bonds;
        let params = mcs.parameters_mut();
        let orig_fused = params.bond_compare_parameters.match_fused_rings;
        let orig_strict = params.bond_compare_parameters.match_fused_rings_strict;
        params.bond_compare_parameters.match_fused_rings = false;
        params.bond_compare_parameters.match_fused_rings_strict = false;
        let res = mcs.match_seed(&mut seed);
        let params = mcs.parameters_mut();
        params.bond_compare_parameters.match_fused_rings = orig_fused;
        params.bond_compare_parameters.match_fused_rings_strict = orig_strict;
        res
    }

    /// `grow`: adds the seeds obtained by adding bonds to this seed to the
    /// seed set of `mcs`, in several stages (each call does one).
    pub fn grow(&mut self, mcs: &mut MaximumCommonSubgraph<'_>) -> Result<()> {
        if !self.can_grow_bigger_than(mcs.max_number_bonds(), mcs.max_number_atoms()) {
            self.growing_stage = NOT_SET; // finished
            return Ok(());
        }
        let qmol = mcs.query_molecule();
        if self.growing_stage == 0 {
            // 0. fill out the list of all directly connected outgoing bonds
            self.fill_new_bonds(qmol, Some(&mut *mcs))?;
            if self.new_bonds.is_empty() {
                self.growing_stage = NOT_SET; // finished
                return Ok(());
            }
            // 1. check and add the biggest child seed with all outgoing bonds
            let mut seed = self.create_from_parent();
            add_new_bonds_to_seed(&self.new_bonds, qmol, &mut seed)?;
            if !seed.can_grow_bigger_than(mcs.max_number_bonds(), mcs.max_number_atoms()) {
                self.growing_stage = NOT_SET;
                return Ok(()); // all children are too small
            }
            seed.match_result = self.match_result.clone();
            // this seed + all extern bonds is a part of the MCS
            let all_matched = mcs.check_if_match_and_append(&mut seed)?;
            self.growing_stage = 1;
            if all_matched && self.new_bonds.len() > 1 {
                return Ok(()); // grow deep first; postpone the next steps
            }
        }
        // 2. check and add all 2^N-1-1 other possible seeds
        if self.new_bonds.len() == 1 {
            self.growing_stage = NOT_SET;
            return Ok(()); // everything has been done
        }
        // OPTIMISATION: check each individual bond first; if this seed plus
        // the individual bond is no MCS, drop the bond from the growing
        let mut num_erased_new_bonds = 0;
        for i in 0..self.new_bonds.len() {
            let new_bond = self.new_bonds[i].clone();
            let mut seed = self.create_from_parent();
            if new_bond.end_atom_idx == NOT_SET {
                seed.add_atom(new_bond.new_atom_idx);
            }
            seed.add_bond(qmol, new_bond.bond_idx)?;
            seed.compute_remaining_size(qmol);
            if seed.can_grow_bigger_than(mcs.max_number_bonds(), mcs.max_number_atoms()) {
                if !self.match_result.is_empty() {
                    seed.match_result = self.match_result.clone();
                }
                if !mcs.check_if_match_and_append(&mut seed)? {
                    // exclude this new bond from growing this seed
                    self.new_bonds[i].bond_idx = NOT_SET;
                    num_erased_new_bonds += 1;
                }
            }
        }
        if num_erased_new_bonds > 0 {
            self.new_bonds.retain(|b| b.bond_idx != NOT_SET);
        }
        // add all other from 2^k-1 possible seeds, k = number of new bonds
        if self.new_bonds.len() > 1 {
            if self.new_bonds.len() >= 64 {
                return Err(RDKitError::Runtime(
                    "Number of new external bonds of a seed must be below 64".into(),
                ));
            }
            let max_composition_value = Composition2N::compute_2n(self.new_bonds.len() as u32) - 1;
            let mut composition = Composition2N::new(max_composition_value, max_composition_value);
            while composition.generate_next() {
                // exclude already processed individual external bonds
                if composition.is_2_power() {
                    continue;
                }
                // exclude the already processed combination of all bonds
                if num_erased_new_bonds == 0 && composition.bit_set() == max_composition_value {
                    continue;
                }
                let mut seed = self.create_from_parent();
                let mut new_atoms_set = vec![false; qmol.num_atoms() as usize];
                for (i, new_bond) in self.new_bonds.iter().enumerate() {
                    if composition.is_set(i) {
                        if new_bond.end_atom_idx == NOT_SET
                            && !new_atoms_set[new_bond.new_atom_idx as usize]
                        {
                            seed.add_atom(new_bond.new_atom_idx);
                            new_atoms_set[new_bond.new_atom_idx as usize] = true;
                        }
                        seed.add_bond(qmol, new_bond.bond_idx)?;
                    }
                }
                seed.compute_remaining_size(qmol);
                if seed.can_grow_bigger_than(mcs.max_number_bonds(), mcs.max_number_atoms()) {
                    seed.match_result = self.match_result.clone();
                    mcs.check_if_match_and_append(&mut seed)?;
                }
            }
        }
        self.growing_stage = NOT_SET; // finished
        Ok(())
    }
}

/// (neighbor atom, bond) pairs of an atom, in the order of the molecule.
fn adjacent(mol: &ROMol, atom: u32) -> Vec<(u32, u32)> {
    mol.atom_neighbors(atom)
        .zip(mol.atom_bond_indices(atom))
        .collect()
}

/// `addNewBondsToSeed`: adds all of `new_bonds` to `seed`, which must be a
/// child of the seed that produced them. Returns the atoms that were added.
fn add_new_bonds_to_seed(
    new_bonds: &[NewBond],
    qmol: &ROMol,
    seed: &mut Seed,
) -> Result<Vec<bool>> {
    let mut new_atoms_set = vec![false; qmol.num_atoms() as usize];
    for new_bond in new_bonds {
        if new_bond.end_atom_idx == NOT_SET && !new_atoms_set[new_bond.new_atom_idx as usize] {
            // (a new bond may close a ring with an atom added just now)
            seed.add_atom(new_bond.new_atom_idx);
            new_atoms_set[new_bond.new_atom_idx as usize] = true;
        }
        seed.add_bond(qmol, new_bond.bond_idx)?;
    }
    seed.remaining_bonds = seed.remaining_bonds.wrapping_sub(new_bonds.len() as u32); // added ALL
    seed.remaining_atoms = seed
        .remaining_atoms
        .wrapping_sub(new_atoms_set.iter().filter(|&&b| b).count() as u32);
    Ok(new_atoms_set)
}
