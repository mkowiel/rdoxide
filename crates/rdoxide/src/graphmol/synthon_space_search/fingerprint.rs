//! `SynthonSpaceFingerprintSearcher.cpp`

use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

use super::details;
use super::hit_set::SynthonSpaceHitSet;
use super::params::{SynthonSpaceSearchParams, ThreadMode};
use super::search_results::SearchResults;
use super::searcher::{
    SearchResultCallback, SearcherBase, SynthonSearcher, check_possible_hits, default_quick_verify,
    default_verify_hit, search, search_with_callback,
};
use super::space::SynthonSpace;
use super::synthon_set::SynthonSet;
use crate::datastructs::bit_ops::{Metric, similarity};
use crate::datastructs::bit_vect::{BitVect, ExplicitBitVect};
use crate::graphmol::fingerprints::{FingerprintFuncArguments, FingerprintGenerator};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::{RDKitError, Result};

/// `FingerprintGenerator::getFingerprint(mol)`
pub fn get_fingerprint(fp_gen: &FingerprintGenerator, mol: &ROMol) -> Result<ExplicitBitVect> {
    fp_gen.fingerprint(
        mol,
        &mut FingerprintFuncArguments {
            from_atoms: None,
            ignore_atoms: None,
            conf_id: -1,
            additional_output: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
        },
    )
}

fn tanimoto(a: &ExplicitBitVect, b: &ExplicitBitVect) -> f64 {
    similarity(Metric::Tanimoto, a, b).unwrap_or(0.0)
}

pub struct FingerprintSearcher<'a> {
    base: SearcherBase<'a>,
    fp_gen: &'a FingerprintGenerator,
    query_fp: ExplicitBitVect,
    /// keyed by the address of the fragment
    frag_fps: HashMap<usize, Arc<ExplicitBitVect>>,
}

fn find_synthon_search_start(
    num_frag_set_bits: u32,
    similarity_cutoff: f64,
    set_num: usize,
    reaction: &SynthonSet,
    pool: &[super::synthon::Synthon],
) -> usize {
    let min_bits = (similarity_cutoff * f64::from(num_frag_set_bits)) as u32;
    let bits_at = |n: usize| {
        let o = reaction
            .fingerprint_ordered_synthon_num(set_num, n)
            .unwrap();
        pool[reaction.synthons()[set_num][o].1]
            .fp()
            .map_or(0, |f| f.num_on_bits())
    };
    let (mut first, mut last) = (0, reaction.synthons()[set_num].len());
    while first < last {
        let mid = first + (last - first) / 2;
        if bits_at(mid) < min_bits {
            first = mid + 1;
        } else {
            last = mid;
        }
    }
    first
}

fn get_hit_synthons(
    frag_fps: &[&ExplicitBitVect],
    similarity_cutoff: f64,
    reaction: &SynthonSet,
    pool: &[super::synthon::Synthon],
    synthon_set_order: &[usize],
) -> Vec<Vec<usize>> {
    let n_sets = reaction.synthons().len();
    let mut synthons_to_use: Vec<Vec<bool>> = reaction
        .synthons()
        .iter()
        .map(|s| vec![false; s.len()])
        .collect();
    let mut frag_sims: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_sets];
    for i in 0..synthon_set_order.len() {
        let so = synthon_set_order[i];
        let max_bits = if similarity_cutoff <= 0.0 {
            frag_fps[i].num_bits()
        } else {
            1 + (f64::from(frag_fps[i].num_on_bits()) / similarity_cutoff) as u32
        };
        let start = find_synthon_search_start(
            frag_fps[i].num_on_bits(),
            similarity_cutoff,
            so,
            reaction,
            pool,
        );
        let synthons = &reaction.synthons()[so];
        let mut frag_matched = false;
        for j in start..synthons.len() {
            // Search in the sorted order, stopping when the number of bits
            // goes above max_bits.
            let synthon_num = reaction.fingerprint_ordered_synthon_num(so, j).unwrap();
            let sfp = pool[synthons[synthon_num].1].fp().unwrap();
            if sfp.num_on_bits() > max_bits {
                break;
            }
            let sim = tanimoto(frag_fps[i], sfp);
            if sim >= similarity_cutoff {
                synthons_to_use[so][synthon_num] = true;
                frag_sims[so].push((synthon_num, sim));
                frag_matched = true;
            }
        }
        if !frag_matched {
            return Vec::new();
        }
    }
    details::expand_bit_set(&mut synthons_to_use);
    let mut ret = details::bit_sets_to_vectors(&synthons_to_use);
    // Order the synthons in descending order of their similarity to the
    // corresponding fragment FP.
    for i in 0..frag_fps.len() {
        if frag_sims[i].is_empty() {
            // filled in by expandBitSet: use all the synthons
            frag_sims[i] = (0..synthons_to_use[i].len()).map(|j| (j, 0.0)).collect();
        } else {
            libcxx_sort_by(&mut frag_sims[i], |a, b| a.1 > b.1);
        }
        ret[i] = frag_sims[i].iter().map(|fs| fs.0).collect();
    }
    ret
}

impl<'a> FingerprintSearcher<'a> {
    pub fn new(
        query: &ROMol,
        fp_gen: &'a FingerprintGenerator,
        params: &SynthonSpaceSearchParams,
        space: &'a SynthonSpace,
        check_type: bool,
    ) -> Result<Self> {
        if check_type && space.has_fingerprints() && fp_gen.info_string() != space.fp_type() {
            return Err(RDKitError::Runtime(format!(
                "The search fingerprints must match those in the database.  You are searching with {} vs {} in the database.",
                fp_gen.info_string(),
                space.fp_type()
            )));
        }
        Ok(FingerprintSearcher {
            base: SearcherBase::new(query, params, space)?,
            fp_gen,
            query_fp: get_fingerprint(fp_gen, query)?,
            frag_fps: HashMap::new(),
        })
    }
}

impl<'a> SynthonSearcher<'a> for FingerprintSearcher<'a> {
    fn base(&self) -> &SearcherBase<'a> {
        &self.base
    }

    fn extra_search_setup(
        &mut self,
        frag_sets: &mut Vec<Vec<Arc<ROMol>>>,
        _end_time: Option<Instant>,
    ) -> Result<bool> {
        let map = details::map_frags_by_smiles(frag_sets)?;
        for (_smi, members) in map {
            let first = &frag_sets[members[0].0][members[0].1];
            let fp = Arc::new(get_fingerprint(self.fp_gen, first)?);
            for &(i, j) in &members {
                self.frag_fps
                    .insert(Arc::as_ptr(&frag_sets[i][j]) as usize, fp.clone());
            }
        }
        Ok(true)
    }

    fn search_frag_set(
        &self,
        frag_set: &[Arc<ROMol>],
        reaction: &'a SynthonSet,
    ) -> Result<Vec<SynthonSpaceHitSet<'a>>> {
        let mut results = Vec::new();
        // It can't be a hit if there are more fragments than synthon sets.
        if frag_set.len() > reaction.synthons().len() {
            return Ok(results);
        }
        let pool = self.base.space.pool();
        let frag_fps: Vec<&ExplicitBitVect> = frag_set
            .iter()
            .map(|f| &*self.frag_fps[&(Arc::as_ptr(f) as usize)])
            .collect();
        let params = &self.base.params;
        for order in details::perm_m_from_n(frag_set.len(), reaction.synthons().len()) {
            // For fingerprints the isotope numbers are ignored.
            let these = get_hit_synthons(
                &frag_fps,
                params.similarity_cutoff - params.frag_similarity_adjuster,
                reaction,
                pool,
                &order,
            );
            if !these.is_empty() {
                let hs = SynthonSpaceHitSet::new_fp(reaction, pool, &these, frag_set);
                if hs.num_hits > 0 {
                    results.push(hs);
                }
            }
        }
        Ok(results)
    }

    fn quick_verify(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> bool {
        if !default_quick_verify(&self.base, hitset, synth_nums) {
            return false;
        }
        let p = &self.base.params;
        self.approx_similarity(hitset, synth_nums)
            >= p.similarity_cutoff - p.approx_similarity_adjuster
    }

    fn approx_similarity(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> f64 {
        let fp = hitset.fp.as_ref().expect("fingerprint hit set");
        let mut full = fp.synthon_fps[0][synth_nums[0]].clone();
        for i in 1..synth_nums.len() {
            full = full.or(fp.synthon_fps[i][synth_nums[i]]);
        }
        if let Some(a) = fp.add_fp {
            full = full.or(a);
        }
        // the subtract FP has already had its bits flipped
        if let Some(s) = fp.subtract_fp {
            full = full.and(s);
        }
        tanimoto(&full, &self.query_fp)
    }

    fn verify_hit(&self, hit: &mut ROMol, rxn_id: &str, synth_names: &[&str]) -> bool {
        if !default_verify_hit(&self.base, hit) {
            return false;
        }
        let Ok(fp) = get_fingerprint(self.fp_gen, hit) else {
            return false;
        };
        let sim = tanimoto(&fp, &self.query_fp);
        let name = details::build_product_name(rxn_id, synth_names);
        if sim > self.base.best_similarity_so_far() {
            let _ = hit.set_prop(cp::_Name, name.as_str(), false);
            self.base.update_best_hit_so_far(hit, sim);
        }
        if sim >= self.base.params.similarity_cutoff {
            let _ = hit.set_prop("Similarity", sim, false);
            let _ = hit.set_prop(cp::_Name, name.as_str(), false);
            return true;
        }
        false
    }
}

impl SynthonSpace {
    fn ensure_fingerprints(&mut self, fp_gen: &FingerprintGenerator) -> Result<()> {
        if self.has_fingerprints() && fp_gen.info_string() != self.fp_type() {
            return Err(RDKitError::Runtime(format!(
                "The search fingerprints must match those in the database.  You are searching with {} vs {} in the database.",
                fp_gen.info_string(),
                self.fp_type()
            )));
        }
        self.build_synthon_fingerprints(fp_gen)
    }

    /// `fingerprintSearch`. Builds the synthon fingerprints first if the
    /// space doesn't have any yet.
    pub fn fingerprint_search(
        &mut self,
        query: &ROMol,
        fp_gen: &FingerprintGenerator,
        params: &SynthonSpaceSearchParams,
    ) -> Result<SearchResults> {
        check_query(query)?;
        self.ensure_fingerprints(fp_gen)?;
        let mut s = FingerprintSearcher::new(query, fp_gen, params, self, false)?;
        search(&mut s, ThreadMode::ThreadFragments)
    }

    /// `fingerprintSearch` returning the hits to a callback in chunks.
    pub fn fingerprint_search_incremental(
        &mut self,
        query: &ROMol,
        fp_gen: &FingerprintGenerator,
        cb: SearchResultCallback<'_>,
        params: &SynthonSpaceSearchParams,
    ) -> Result<()> {
        check_query(query)?;
        self.ensure_fingerprints(fp_gen)?;
        let mut s = FingerprintSearcher::new(query, fp_gen, params, self, false)?;
        search_with_callback(&mut s, cb, ThreadMode::ThreadFragments)
    }

    /// `fingerprintSearch` on the possible hits between two lines of the
    /// possible hits file.
    pub fn fingerprint_search_lines(
        &self,
        query: &ROMol,
        fp_gen: &FingerprintGenerator,
        params: &SynthonSpaceSearchParams,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SearchResults> {
        check_query(query)?;
        check_lines(start_line, finish_line)?;
        let s = FingerprintSearcher::new(query, fp_gen, params, self, false)?;
        check_possible_hits(&s, start_line, finish_line)
    }
}

pub(super) fn check_query(query: &ROMol) -> Result<()> {
    if query.num_atoms() == 0 {
        return Err(RDKitError::precondition("Search query must contain atoms."));
    }
    Ok(())
}

pub(super) fn check_lines(start_line: u64, finish_line: u64) -> Result<()> {
    if start_line >= finish_line {
        return Err(RDKitError::precondition(
            "Search startLine must be less than finishLine.",
        ));
    }
    Ok(())
}
