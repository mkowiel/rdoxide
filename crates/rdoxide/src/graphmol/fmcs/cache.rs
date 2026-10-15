//! The caches and helpers of the seed growing (upstream
//! `FMCS/DuplicatedSeedCache.h`, `FMCS/SubstructureCache.h` and
//! `FMCS/Composition2N.h`).

use std::collections::HashMap;

use super::graph::Graph;
use super::seed::Seed;

/// Key of the duplicated seed cache: the sorted atom and bond indices of the
/// seed's fragment.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct DupCacheKey {
    atoms: Vec<u32>,
    bonds: Vec<u32>,
}

impl DupCacheKey {
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn add_atom(&mut self, i: u32) {
        let pos = self.atoms.partition_point(|&x| x < i);
        self.atoms.insert(pos, i);
    }

    pub fn add_bond(&mut self, i: u32) {
        let pos = self.bonds.partition_point(|&x| x < i);
        self.bonds.insert(pos, i);
    }
}

/// `DuplicatedSeedCache`: whether a seed with these atoms and bonds matched.
#[derive(Debug, Default)]
pub struct DuplicatedSeedCache {
    index: HashMap<DupCacheKey, bool>,
    max_atoms: usize,
}

impl DuplicatedSeedCache {
    pub fn clear(&mut self) {
        self.index.clear();
        self.max_atoms = 0;
    }

    pub fn find(&self, key: &DupCacheKey) -> Option<bool> {
        if key.num_atoms() > self.max_atoms {
            return None; // fast check against the biggest key in the cache
        }
        self.index.get(key).copied()
    }

    pub fn add(&mut self, key: &DupCacheKey, found: bool) {
        self.max_atoms = self.max_atoms.max(key.num_atoms());
        self.index.entry(key.clone()).or_insert(found);
    }
}

/// Hash of a seed (`SubstructureCache::HashKey`): a Morgan-like code over the
/// atom and bond labels of the query molecule. Arithmetic follows the LP64
/// build of upstream (`unsigned long` is 64 bit, `unsigned int` 32 bit).
pub fn seed_hash_key(seed: &Seed, query_atom_labels: &[u32], query_bond_labels: &[u32]) -> u64 {
    let nv = seed.num_atoms();
    let ne = seed.num_bonds();
    let n_iterations = ne.min(5);
    let mut curr: Vec<u64> = (0..nv)
        .map(|i| u64::from(query_atom_labels[seed.fragment.atoms[i] as usize]))
        .collect();
    let mut prev = vec![0u64; nv];
    for _ in 0..n_iterations {
        prev.copy_from_slice(&curr);
        for seed_bond_idx in 0..ne {
            let bond = seed.fragment.bonds[seed_bond_idx];
            let order = query_bond_labels[bond as usize];
            let (begin, end) = seed.fragment.bond_ends[seed_bond_idx];
            let atom1 = seed.fragment.seed_atom_idx_map[&begin] as usize;
            let atom2 = seed.fragment.seed_atom_idx_map[&end] as usize;
            let v1 = prev[atom1] as u32;
            let v2 = prev[atom2] as u32;
            let term = |v: u32| {
                u64::from(
                    v.wrapping_mul(v)
                        .wrapping_add(v.wrapping_add(23).wrapping_mul(order.wrapping_add(1721))),
                )
            };
            curr[atom1] = curr[atom1].wrapping_add(term(v2));
            curr[atom2] = curr[atom2].wrapping_add(term(v1));
        }
    }
    let mut result = 0u64;
    for code in curr {
        result = result.wrapping_add(code.wrapping_mul(code.wrapping_add(6849)).wrapping_add(29));
    }
    result
}

/// `SubstructureCache`: seeds that are known to be contained in the query,
/// stored by hash; a hash is not unique, so each entry holds several graphs.
#[derive(Debug, Default)]
pub struct SubstructureCache {
    value_storage: Vec<Vec<Graph>>,
    numeric_index: HashMap<u64, usize>,
}

impl SubstructureCache {
    pub fn clear(&mut self) {
        self.value_storage.clear();
        self.numeric_index.clear();
    }

    /// The index of the entry of `key`, if any.
    pub fn find(&self, key: u64) -> Option<usize> {
        self.numeric_index.get(&key).copied()
    }

    pub fn entry(&self, index: usize) -> &[Graph] {
        &self.value_storage[index]
    }

    /// Stores the topology of `seed` in the entry `entry` (a new one when
    /// `None`).
    pub fn add(&mut self, seed: &Seed, key: u64, entry: Option<usize>) {
        let index = match entry {
            Some(i) => i,
            None => {
                self.value_storage.push(Vec::new());
                self.value_storage.len() - 1
            }
        };
        self.value_storage[index].push(seed.topology.clone());
        // like upstream, the index points at the last storage entry
        let last = self.value_storage.len() - 1;
        self.numeric_index.entry(key).or_insert(last);
    }
}

/// `Composition2N`: generator of the 2^N-1 possible bit combinations, the
/// biggest seeds first.
#[derive(Debug)]
pub struct Composition2N {
    bits: u64,
    inverse_bits: u64,
    max_value: u64,
    value_mask: u64,
}

impl Composition2N {
    pub fn new(max_value: u64, value_mask: u64) -> Self {
        Composition2N {
            bits: 0,
            inverse_bits: 0,
            max_value,
            value_mask,
        }
    }

    pub fn compute_2n(power: u32) -> u64 {
        1u64 << power
    }

    pub fn bit_set(&self) -> u64 {
        self.inverse_bits
    }

    pub fn generate_next(&mut self) -> bool {
        self.bits += 1;
        if self.bits <= self.max_value {
            self.inverse_bits = (!self.bits).wrapping_add(1) & self.value_mask;
            true
        } else {
            false
        }
    }

    /// exactly one bit is set
    pub fn is_2_power(&self) -> bool {
        self.bit_set().count_ones() == 1
    }

    pub fn is_set(&self, bit: usize) -> bool {
        self.bit_set() & (1u64 << bit) != 0
    }
}
