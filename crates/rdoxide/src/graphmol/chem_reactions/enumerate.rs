//! Reaction enumeration: strategies for walking through the combinations of
//! building blocks, and the `EnumerateLibrary` engine built on them.
//!
//! Ported from RDKit `Code/GraphMol/ChemReactions/Enumerate/` (2026.09.1).
//! Serialization through boost archives is not provided, which corresponds to
//! an RDKit built without `RDK_USE_BOOST_SERIALIZATION`.
//!
//! Copyright (C) 2015 Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeSet;
use std::fmt::Write as _;

use super::{ChemicalReaction, runner};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::random::MinStdRand;
use crate::{RDKitError, Result, rd_log_info};

/// `EnumerationTypes::BBS`: the building blocks of each reactant.
pub type Bbs = Vec<Vec<ROMol>>;
/// `EnumerationTypes::RGROUPS`: one index into each set of building blocks.
pub type RGroups = Vec<u64>;

/// `EnumerationStrategyBase::EnumerationOverflow`
pub const ENUMERATION_OVERFLOW: u64 = u64::MAX;

/// `computeNumProducts`: the product of the sizes, or the overflow marker.
pub fn compute_num_products(sizes: &[u64]) -> u64 {
    let mut n: u64 = 1;
    for &s in sizes {
        if s != 0 && u64::MAX / s < n {
            return ENUMERATION_OVERFLOW;
        }
        n *= s;
    }
    n
}

/// `getSizesFromBBs`
pub fn sizes_from_bbs(bbs: &Bbs) -> RGroups {
    bbs.iter().map(|b| b.len() as u64).collect()
}

/// `getReactantsFromRGroups`
pub fn reactants_from_rgroups(bbs: &Bbs, rgroups: &[u64]) -> Result<Vec<ROMol>> {
    if bbs.len() != rgroups.len() {
        return Err(RDKitError::precondition(
            "BBS and RGROUPS must have the same # reactants",
        ));
    }
    bbs.iter()
        .zip(rgroups)
        .map(|(b, &i)| {
            b.get(i as usize)
                .cloned()
                .ok_or_else(|| RDKitError::precondition("RGROUPS index out of range"))
        })
        .collect()
}

/// A bounded random index generator: `boost::uniform_int_distribution<>(0, n - 1)`.
fn draw(rng: &mut MinStdRand, size: u64) -> u64 {
    rng.uniform_int(0, (size as i64 - 1) as i32) as u64
}

#[derive(Debug, Clone, Default)]
struct RandomAllBbs {
    offset: usize,
    max_offset: usize,
}

#[derive(Debug, Clone, Default)]
struct EvenPairs {
    used_count: Vec<i64>,
    var_used: Vec<Vec<u64>>,
    pair_used: Vec<Vec<u64>>,
    pair_counts: Vec<Vec<u64>>,
    selected: BTreeSet<u64>,
    seed: u64,
    m: u64,
    a: u64,
    b: u64,
    nslack: u64,
    min_nslack: u64,
    rejected_period: u64,
    rejected_unique: u64,
    rejected_slack_condition: u64,
    rejected_bb_sampling_condition: u64,
}

#[derive(Debug, Clone)]
enum Kind {
    Cartesian,
    Random(MinStdRand),
    RandomAllBbs(MinStdRand, RandomAllBbs),
    EvenPairs(Box<EvenPairs>),
}

/// The enumeration strategies (`CartesianProductStrategy`,
/// `RandomSampleStrategy`, `RandomSampleAllBBsStrategy` and
/// `EvenSamplePairsStrategy`) behind the `EnumerationStrategyBase` interface.
#[derive(Debug, Clone)]
pub struct EnumerationStrategy {
    permutation: RGroups,
    sizes: RGroups,
    num_permutations: u64,
    processed: u64,
    kind: Kind,
}

impl EnumerationStrategy {
    fn with_kind(kind: Kind) -> Self {
        EnumerationStrategy {
            permutation: Vec::new(),
            sizes: Vec::new(),
            num_permutations: 0,
            processed: 0,
            kind,
        }
    }

    pub fn cartesian_product() -> Self {
        Self::with_kind(Kind::Cartesian)
    }

    pub fn random_sample() -> Self {
        Self::with_kind(Kind::Random(MinStdRand::new(1)))
    }

    pub fn random_sample_all_bbs() -> Self {
        Self::with_kind(Kind::RandomAllBbs(
            MinStdRand::new(1),
            RandomAllBbs::default(),
        ))
    }

    pub fn even_sample_pairs() -> Self {
        Self::with_kind(Kind::EvenPairs(Box::default()))
    }

    /// `type()`
    pub fn type_name(&self) -> &'static str {
        match self.kind {
            Kind::Cartesian => "CartesianProductStrategy",
            Kind::Random(_) => "RandomSampleStrategy",
            Kind::RandomAllBbs(..) => "RandomSampleAllBBsStrategy",
            Kind::EvenPairs(_) => "EvenSamplePairsStrategy",
        }
    }

    /// `initialize(rxn, bbs)`: only the number of building blocks of each
    /// reactant matters to the strategies.
    pub fn initialize(&mut self, sizes: &[u64]) -> Result<()> {
        self.sizes = sizes.to_vec();
        self.permutation = vec![0; sizes.len()];
        self.num_permutations = compute_num_products(sizes);
        self.processed = 0;
        match &mut self.kind {
            Kind::Cartesian => {}
            Kind::Random(_) => {}
            Kind::RandomAllBbs(_, st) => {
                st.offset = 0;
                st.max_offset = sizes.iter().copied().max().unwrap_or(0) as usize;
            }
            Kind::EvenPairs(st) => {
                if self.num_permutations == ENUMERATION_OVERFLOW {
                    return Err(RDKitError::precondition(
                        "Cannot represent all permutations for the even sampler",
                    ));
                }
                st.initialize(sizes, self.num_permutations);
            }
        }
        Ok(())
    }

    /// `operator bool`: are there more permutations left?
    pub fn has_next(&self) -> bool {
        match self.kind {
            Kind::Cartesian => {
                self.num_permutations == ENUMERATION_OVERFLOW
                    || self.processed < self.num_permutations
            }
            _ => true,
        }
    }

    /// `getPosition()`
    pub fn position(&self) -> &[u64] {
        &self.permutation
    }

    /// `getNumPermutations()`
    pub fn num_permutations(&self) -> u64 {
        self.num_permutations
    }

    /// `getPermutationIdx()`
    pub fn permutation_idx(&self) -> u64 {
        self.processed
    }

    /// `next()`: advance and return the new position.
    pub fn advance(&mut self) -> Result<&[u64]> {
        match &mut self.kind {
            Kind::Cartesian => {
                if self.processed != 0 {
                    increment_cartesian(
                        &mut self.permutation,
                        &self.sizes,
                        self.num_permutations,
                        self.processed,
                    );
                }
                self.processed += 1;
            }
            Kind::Random(rng) => {
                for (p, &s) in self.permutation.iter_mut().zip(&self.sizes) {
                    *p = draw(rng, s);
                }
                self.processed += 1;
            }
            Kind::RandomAllBbs(rng, st) => {
                if st.offset >= st.max_offset {
                    for (p, &s) in self.permutation.iter_mut().zip(&self.sizes) {
                        *p = draw(rng, s);
                    }
                    st.offset = 0;
                } else {
                    for (p, &s) in self.permutation.iter_mut().zip(&self.sizes) {
                        *p = (*p + 1) % s;
                    }
                    st.offset += 1;
                }
                self.processed += 1;
            }
            Kind::EvenPairs(st) => {
                st.nslack = 0;
                let mut found = None;
                'outer: while self.processed < self.num_permutations {
                    for _ in 0..st.m {
                        st.seed = (st.seed * st.a + st.b) % st.m;
                        if st.seed > self.num_permutations {
                            st.rejected_period += 1;
                            continue;
                        } else if st.selected.contains(&st.seed) {
                            st.rejected_unique += 1;
                            continue;
                        }
                        let seed = st.seed;
                        if st.try_add(seed, &self.sizes, &mut self.permutation) {
                            self.processed += 1;
                            found = Some(seed);
                            break 'outer;
                        }
                    }
                    // loosen heuristic
                    st.nslack += 1;
                    st.min_nslack += 1;
                }
                let Some(seed) = found else {
                    return Err(RDKitError::Runtime("Ran out of molecules".into()));
                };
                decode(seed, &self.sizes, &mut self.permutation);
            }
        }
        Ok(&self.permutation)
    }

    /// `skip(count)`
    pub fn skip(&mut self, count: u64) -> Result<()> {
        for _ in 0..count {
            self.advance()?;
        }
        Ok(())
    }

    /// `EvenSamplePairsStrategy::stats()`; empty for the other strategies.
    pub fn stats(&self) -> Option<String> {
        match &self.kind {
            Kind::EvenPairs(st) => Some(st.stats(&self.sizes, self.processed)),
            _ => None,
        }
    }
}

/// `CartesianProductStrategy::next(rowToIncrement)`, iteratively.
fn increment_cartesian(perm: &mut [u64], sizes: &[u64], num_permutations: u64, processed: u64) {
    let has_next = num_permutations == ENUMERATION_OVERFLOW || processed < num_permutations;
    if !has_next {
        return;
    }
    for row in 0..perm.len() {
        perm[row] += 1;
        if perm[row] < sizes[row] {
            return;
        }
        perm[row] = 0;
    }
}

/// `EvenSamplePairsStrategy::decode`: unpack an integer into an RGroup selection.
fn decode(mut seed: u64, sizes: &[u64], perm: &mut [u64]) {
    for j in (0..sizes.len()).rev() {
        perm[j] = seed % sizes[j];
        seed /= sizes[j];
    }
}

impl EvenPairs {
    fn initialize(&mut self, sizes: &[u64], num_permutations: u64) {
        let npos = sizes.len();
        self.used_count = vec![0; npos];
        self.var_used = sizes.iter().map(|&s| vec![0; s as usize]).collect();
        let nmonomers: usize = sizes.iter().map(|&s| s as usize).sum();
        self.pair_used = vec![vec![0; nmonomers]; nmonomers];
        self.pair_counts = vec![vec![0; npos]; npos];
        // find the modulus of the random number generator
        self.m = 1;
        while self.m < num_permutations {
            self.m *= 2;
        }
        self.a = 5;
        self.b = 7;
        self.seed = 0;
        self.nslack = 0;
        self.rejected_period = 0;
        self.rejected_unique = 0;
        self.rejected_slack_condition = 0;
        self.rejected_bb_sampling_condition = 0;
        self.selected.clear();
    }

    /// Tries to add the encoded position; rejected unless the building
    /// blocks and building block pairs stay evenly sampled.
    fn try_add(&mut self, iseed: u64, rgroups: &[u64], digits: &mut [u64]) -> bool {
        decode(iseed, rgroups, digits);
        let n = rgroups.len();
        let mut islack: u64 = 0;
        for i in 0..n {
            let used = self.var_used[i][digits[i] as usize];
            if used != 0 {
                islack += used;
            }
            if islack > self.nslack {
                self.rejected_slack_condition += 1;
                return false;
            }
        }

        islack = 0;
        let mut ioffset = 0usize;
        for i in 0..n {
            let mut joffset = 0usize;
            for j in 0..n {
                if j == i {
                    continue;
                }
                let ii = digits[i] as usize + ioffset;
                let jj = digits[j] as usize + joffset;
                if self.pair_used[ii][jj] > 0 {
                    let numer = self.pair_used[ii][jj] as f64;
                    let denom = ((rgroups[i] as f64) * (rgroups[j] as f64)).sqrt();
                    islack = (numer / denom) as i32 as u64;
                }
                joffset += rgroups[j] as usize;
            }
            ioffset += rgroups[i] as usize;
        }
        if islack > self.nslack {
            self.rejected_bb_sampling_condition += 1;
            return false;
        }

        // keep track of building block usage
        for i in 0..n {
            let d = digits[i] as usize;
            if self.var_used[i][d] == 0 {
                self.used_count[i] += 1;
            }
            self.var_used[i][d] += 1;
            if self.used_count[i] == rgroups[i] as i64 {
                // complete variable scan => initialize
                if self.nslack > self.min_nslack && rgroups[i] > 1 {
                    self.nslack = self.min_nslack;
                }
                self.used_count[i] = 0;
                for j in 0..rgroups[i] as usize {
                    self.var_used[i][j] = self.var_used[i][j].wrapping_sub(1);
                    if self.var_used[i][j] > 0 {
                        self.used_count[i] += 1;
                    }
                }
            }
        }

        // keep track of building block pair usage
        let mut ioffset = 0usize;
        for i in 0..n {
            let mut joffset = 0usize;
            for j in 0..n {
                if j != i {
                    let ii = digits[i] as usize + ioffset;
                    let jj = digits[j] as usize + joffset;
                    if self.pair_used[ii][jj] == 0 {
                        self.pair_counts[i][j] += 1;
                    }
                    self.pair_used[ii][jj] += 1;
                    if self.pair_counts[i][j] >= rgroups[i] * rgroups[j] {
                        // all pairs visited
                        if self.nslack > self.min_nslack && (rgroups[i] > 1 || rgroups[j] > 1) {
                            self.nslack = self.min_nslack;
                        }
                        self.pair_counts[i][j] = 0;
                        for ii in 0..rgroups[i] as usize {
                            for jj in 0..rgroups[j] as usize {
                                let cell = &mut self.pair_used[ioffset + ii][joffset + jj];
                                *cell = cell.wrapping_sub(1);
                                if *cell > 0 {
                                    self.pair_counts[i][j] += 1;
                                }
                            }
                        }
                    }
                }
                joffset += rgroups[j] as usize;
            }
            ioffset += rgroups[i] as usize;
        }

        self.selected.insert(iseed);
        true
    }

    fn stats(&self, nvars: &[u64], processed: u64) -> String {
        let npos = nvars.len();
        let mut ss = String::from("#BEGIN# BBSTAT\n");
        for i in 0..npos {
            if nvars[i] == 1 {
                continue;
            }
            let maxcount = self.var_used[i].iter().copied().max().unwrap_or(0);
            let _ = write!(
                ss,
                "{}\t{}\t{:6.2}",
                i + 1,
                nvars[i],
                processed as f64 / nvars[i] as f64
            );
            for l in 0..=maxcount {
                let mut n = 0;
                for j in 0..nvars[i] as usize {
                    if self.var_used[i][j] == l {
                        n += 1;
                    }
                    if n > 0 {
                        let _ = write!(ss, "\t{l}|{n}");
                    }
                }
                ss.push('\n');
            }
        }
        ss.push_str("#END# BBSTAT\n#BEGIN# PAIRSTAT\n");
        let mut ioffset = 0usize;
        for i in 0..npos {
            if nvars[i] != 1 {
                let mut joffset = 0usize;
                for j in 0..npos {
                    if nvars[j] != 1 && j > i {
                        let mut maxcount = 0;
                        for ii in 0..nvars[i] as usize {
                            for jj in 0..nvars[j] as usize {
                                maxcount = maxcount.max(self.pair_used[ii + ioffset][jj + joffset]);
                            }
                        }
                        let _ = write!(
                            ss,
                            "{}\t{}\t{}\t{}\t{:6.2}",
                            i + 1,
                            j + 1,
                            nvars[i],
                            nvars[j],
                            processed as f64 / (nvars[i] * nvars[j]) as f64
                        );
                        for l in 0..=maxcount {
                            let mut n = 0;
                            for ii in 0..nvars[i] as usize {
                                for jj in 0..nvars[j] as usize {
                                    if l == self.pair_used[ii + ioffset][jj + joffset] {
                                        n += 1;
                                    }
                                }
                            }
                            if n > 0 {
                                let _ = write!(ss, "\t{l}|{n}");
                            }
                        }
                        ss.push('\n');
                    }
                    joffset += nvars[j] as usize;
                }
            }
            ioffset += nvars[i] as usize;
        }
        ss.push_str("#END# PAIRSTAT\n");
        let _ = writeln!(ss, "Rejected Period: {}", self.rejected_period);
        let _ = writeln!(ss, "Rejected (dupes): {}", self.rejected_unique);
        let _ = writeln!(
            ss,
            "Rejected Slack Conditions: {}",
            self.rejected_slack_condition
        );
        let _ = writeln!(
            ss,
            "Rejected Pair Sampling: {}",
            self.rejected_bb_sampling_condition
        );
        ss
    }
}

// ---------------------------------------------------------------- library

/// `EnumerationParams`
#[derive(Debug, Clone, Copy)]
pub struct EnumerationParams {
    /// how many times the reactant template can match a reagent
    pub reagent_max_match_count: i32,
    /// require that all products of the reagent plus the product templates
    /// pass sanitization
    pub sane_partial_products: bool,
}

impl Default for EnumerationParams {
    fn default() -> Self {
        EnumerationParams {
            reagent_max_match_count: i32::MAX,
            sane_partial_products: false,
        }
    }
}

fn count_matches(bb: &ROMol, query: &ROMol, max_matches: u32) -> Result<usize> {
    let ps = SubstructMatchParameters {
        uniquify: true,
        use_chirality: true,
        use_query_query_matches: false,
        recursion_possible: true,
        max_matches: max_matches + 1,
        num_threads: 1,
        ..SubstructMatchParameters::default()
    };
    Ok(substruct_match(bb, query, &ps)?.len())
}

/// `removeNonmatchingReagents`
pub fn remove_nonmatching_reagents(
    rxn: &ChemicalReaction,
    bbs: &Bbs,
    params: &EnumerationParams,
) -> Result<Bbs> {
    if bbs.len() > rxn.num_reactant_templates() {
        return Err(RDKitError::precondition(
            "Number of Reagents not compatible with reaction templates",
        ));
    }
    let mut result: Bbs = vec![Vec::new(); bbs.len()];
    for (reactant_idx, reagents) in bbs.iter().enumerate() {
        let mut removed = 0;
        let max_matches = if params.reagent_max_match_count == i32::MAX {
            0
        } else {
            params.reagent_max_match_count as u32
        };
        let template = &rxn.reactants()[reactant_idx];
        for mol in reagents {
            let matches = count_matches(mol, template, max_matches)?;
            let mut remove = matches == 0 || matches > params.reagent_max_match_count as usize;
            if !remove && params.sane_partial_products {
                // see if any sane products come out of this reagent
                for partial in runner::run_reactant(rxn, mol, reactant_idx, 1)? {
                    let sane = partial
                        .into_iter()
                        .filter(|p| sanitize_mol(&mut p.clone(), SANITIZE_ALL).is_ok())
                        .count();
                    if sane == 0 {
                        // if any product template has no sane products, we bail
                        remove = true;
                        break;
                    }
                }
            }
            if remove {
                removed += 1;
            } else {
                result[reactant_idx].push(mol.clone());
            }
        }
        if removed > 0 {
            rd_log_info!("Removed {removed} non matching reagents at template {reactant_idx}\n");
        }
    }
    Ok(result)
}

/// `EnumerateLibrary` (and `EnumerateLibraryBase`): enumerates the products
/// of a reaction over sets of building blocks.
#[derive(Debug, Clone, Default)]
pub struct EnumerateLibrary {
    rxn: ChemicalReaction,
    enumerator: Option<EnumerationStrategy>,
    initial_enumerator: Option<EnumerationStrategy>,
    bbs: Bbs,
}

impl EnumerateLibrary {
    /// An empty library, as made by the default constructor.
    pub fn empty() -> Self {
        Self::default()
    }

    /// `EnumerateLibrary(rxn, reagents, [enumerator,] params)`; the default
    /// strategy is the Cartesian product.
    pub fn new(
        rxn: &ChemicalReaction,
        reagents: &Bbs,
        enumerator: Option<&EnumerationStrategy>,
        params: &EnumerationParams,
    ) -> Result<Self> {
        let mut rxn = rxn.clone();
        rxn.init_reactant_matchers(false);
        let bbs = remove_nonmatching_reagents(&rxn, reagents, params)?;
        let mut enumerator = enumerator
            .cloned()
            .unwrap_or_else(EnumerationStrategy::cartesian_product);
        enumerator.initialize(&sizes_from_bbs(&bbs))?;
        Ok(EnumerateLibrary {
            rxn,
            initial_enumerator: Some(enumerator.clone()),
            enumerator: Some(enumerator),
            bbs,
        })
    }

    fn strategy(&self) -> Result<&EnumerationStrategy> {
        self.enumerator
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("Null Enumerator"))
    }

    /// `operator bool`
    pub fn has_next(&self) -> Result<bool> {
        Ok(self
            .enumerator
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("Null enumeration strategy"))?
            .has_next())
    }

    /// `reset()`
    pub fn reset(&mut self) {
        if let Some(initial) = &self.initial_enumerator {
            self.enumerator = Some(initial.clone());
        }
    }

    /// `resetState()`
    pub fn reset_state(&mut self) -> Result<()> {
        let initial = self
            .initial_enumerator
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("Unset initial enumerator"))?;
        self.enumerator = Some(initial.clone());
        Ok(())
    }

    pub fn reaction(&self) -> &ChemicalReaction {
        &self.rxn
    }

    /// `getEnumerator()`
    pub fn enumerator(&self) -> Result<&EnumerationStrategy> {
        self.strategy()
    }

    /// `getPosition()`
    pub fn position(&self) -> Result<&[u64]> {
        Ok(self.strategy()?.position())
    }

    /// `getReagents()`
    pub fn reagents(&self) -> &Bbs {
        &self.bbs
    }

    /// `next()`: the products of the next combination of building blocks.
    pub fn next_products(&mut self) -> Result<Vec<Vec<ROMol>>> {
        if !self.has_next()? {
            return Err(RDKitError::precondition("No more enumerations"));
        }
        let indices = self
            .enumerator
            .as_mut()
            .ok_or_else(|| RDKitError::precondition("Null enumeration strategy"))?
            .advance()?
            .to_vec();
        let reactants = reactants_from_rgroups(&self.bbs, &indices)?;
        runner::run_reactants(&self.rxn, &reactants, 1000, 1)
    }

    /// `nextSmiles()`
    pub fn next_smiles(&mut self) -> Result<Vec<Vec<String>>> {
        let products = self.next_products()?;
        products
            .into_iter()
            .map(|set| {
                set.into_iter()
                    .map(|mut m| mol_to_smiles(&mut m, &SmilesWriteParams::default()))
                    .collect()
            })
            .collect()
    }
}

/// `EnumerateLibraryCanSerialize`: boost archives are not provided.
pub fn enumerate_library_can_serialize() -> bool {
    false
}
