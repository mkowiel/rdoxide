//! `SynthonSpaceRascalSearcher.cpp`

use std::sync::Arc;
use std::time::Instant;

use super::details;
use super::fingerprint::{check_lines, check_query};
use super::hit_set::SynthonSpaceHitSet;
use super::params::{SynthonSpaceSearchParams, ThreadMode};
use super::search_results::SearchResults;
use super::searcher::{
    SearchResultCallback, SearcherBase, SynthonSearcher, check_possible_hits, default_quick_verify,
    default_verify_hit, search, search_with_callback,
};
use super::space::SynthonSpace;
use super::synthon::Synthon;
use super::synthon_set::SynthonSet;
use crate::Result;
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::rascal_mces::{RascalOptions, rascal_mces};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;

pub struct RascalSearcher<'a> {
    base: SearcherBase<'a>,
    options: RascalOptions,
    frag_options: RascalOptions,
}

fn find_synthon_search_start(
    frag_atoms_and_bonds: u32,
    cutoff: f64,
    set_num: usize,
    reaction: &SynthonSet,
    pool: &[Synthon],
) -> usize {
    let min = (cutoff * f64::from(frag_atoms_and_bonds)) as u32;
    let (mut first, mut last) = (0, reaction.synthons()[set_num].len());
    while first < last {
        let mid = first + (last - first) / 2;
        let o = reaction.rascal_ordered_synthon_num(set_num, mid).unwrap();
        let sm = pool[reaction.synthons()[set_num][o].1].search_mol();
        if sm.num_atoms() + sm.num_bonds() < min {
            first = mid + 1;
        } else {
            last = mid;
        }
    }
    first
}

fn get_hit_synthons(
    frag_set: &[Arc<ROMol>],
    options: &RascalOptions,
    reaction: &SynthonSet,
    pool: &[Synthon],
    synthon_set_order: &[usize],
) -> Result<Vec<Vec<usize>>> {
    // Match fragments against the smallest synthon sets first.
    let mut frag_orders: Vec<(usize, usize)> = synthon_set_order
        .iter()
        .enumerate()
        .map(|(i, &so)| (i, reaction.synthons()[so].len()))
        .collect();
    libcxx_sort_by(&mut frag_orders, |a, b| a.1 < b.1);
    let mut synthons_to_use: Vec<Vec<bool>> = reaction
        .synthons()
        .iter()
        .map(|s| vec![false; s.len()])
        .collect();
    for i in 0..synthon_set_order.len() {
        let frag_num = frag_orders[i].0;
        let so = synthon_set_order[frag_num];
        let frag = &frag_set[frag_num];
        let frag_atoms_and_bonds = frag.num_atoms() + frag.num_bonds();
        let max = 1 + (f64::from(frag_atoms_and_bonds) / options.similarity_threshold) as u32;
        let start = find_synthon_search_start(
            frag_atoms_and_bonds,
            options.similarity_threshold,
            so,
            reaction,
            pool,
        );
        let synthons = &reaction.synthons()[so];
        let mut frag_matched = false;
        for j in start..synthons.len() {
            let synthon_num = reaction.rascal_ordered_synthon_num(so, j).unwrap();
            let sm = pool[synthons[synthon_num].1].search_mol();
            if sm.num_atoms() + sm.num_bonds() > max {
                break;
            }
            if !rascal_mces(frag, sm, options)?.is_empty() {
                synthons_to_use[so][synthon_num] = true;
                frag_matched = true;
            }
        }
        if !frag_matched {
            return Ok(Vec::new());
        }
    }
    details::expand_bit_set(&mut synthons_to_use);
    let mut ret = details::bit_sets_to_vectors(&synthons_to_use);
    for (i, r) in ret.iter_mut().enumerate() {
        let si = &reaction.synthons()[i];
        libcxx_sort_by(r, |&a, &b| {
            pool[si[a].1].orig_mol().num_atoms() < pool[si[b].1].orig_mol().num_atoms()
        });
    }
    Ok(ret)
}

impl<'a> RascalSearcher<'a> {
    pub fn new(
        query: &ROMol,
        options: &RascalOptions,
        params: &SynthonSpaceSearchParams,
        space: &'a SynthonSpace,
    ) -> Result<Self> {
        let mut frag_options = options.clone();
        frag_options.similarity_threshold -= params.frag_similarity_adjuster;
        Ok(RascalSearcher {
            base: SearcherBase::new(query, params, space)?,
            options: options.clone(),
            frag_options,
        })
    }
}

impl<'a> SynthonSearcher<'a> for RascalSearcher<'a> {
    fn base(&self) -> &SearcherBase<'a> {
        &self.base
    }

    fn extra_search_setup(
        &mut self,
        frag_sets: &mut Vec<Vec<Arc<ROMol>>>,
        end_time: Option<Instant>,
    ) -> Result<bool> {
        let mut num_done = 100;
        for frag_set in frag_sets.iter_mut() {
            for frag in frag_set.iter_mut() {
                if !frag.ring_info().is_initialized() {
                    let f = Arc::get_mut(frag).expect("unshared fragment");
                    symmetrize_sssr(f, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
                }
                num_done -= 1;
                if num_done == 0 {
                    num_done = 100;
                    if details::check_time_out(end_time) {
                        return Ok(false);
                    }
                }
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
        if frag_set.len() > reaction.synthons().len() {
            return Ok(results);
        }
        let pool = self.base.space.pool();
        for order in details::perm_m_from_n(frag_set.len(), reaction.synthons().len()) {
            // Rascal ignores isotope numbers, which makes things easier.
            let these = get_hit_synthons(frag_set, &self.frag_options, reaction, pool, &order)?;
            if !these.is_empty() {
                let hs = SynthonSpaceHitSet::new(reaction, pool, &these, frag_set);
                if hs.num_hits > 0 {
                    results.push(hs);
                }
            }
        }
        Ok(results)
    }

    fn quick_verify(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> bool {
        default_quick_verify(&self.base, hitset, synth_nums)
            && self.approx_similarity(hitset, synth_nums) >= self.options.similarity_threshold
    }

    fn approx_similarity(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> f64 {
        // If the query is an exact substructure of the product, that's an
        // upper bound on the Johnson similarity.
        let qbit = (self.base.query.num_atoms() + self.base.query.num_bonds()) as i32;
        let (mut num_atoms, mut num_bonds, mut num_conns) = (0i32, 0i32, 0i32);
        for (i, &n) in synth_nums.iter().enumerate() {
            let s = hitset.synthons_to_use[i][n].1;
            num_atoms += s.num_heavy_atoms() as i32;
            num_bonds += s.orig_mol().num_bonds() as i32 - s.num_dummies() as i32;
            num_conns += s.num_dummies() as i32;
        }
        // Pairs of dummies will form a bond in the final result.
        num_conns /= 2;
        f64::from(qbit) / f64::from(num_atoms + num_bonds + num_conns)
    }

    fn verify_hit(&self, hit: &mut ROMol, rxn_id: &str, synth_names: &[&str]) -> bool {
        if !default_verify_hit(&self.base, hit) {
            return false;
        }
        let Ok(res) = rascal_mces(hit, &self.base.query, &self.options) else {
            return false;
        };
        // Rascal reports all matches that proceed to full MCES elucidation,
        // even if the final similarity ends up below the threshold.
        let Some(first) = res.first() else {
            return false;
        };
        let sim = first.similarity();
        let name = details::build_product_name(rxn_id, synth_names);
        if sim > self.base.best_similarity_so_far() {
            let _ = hit.set_prop(cp::_Name, name.as_str(), false);
            self.base.update_best_hit_so_far(hit, sim);
        }
        if sim >= self.options.similarity_threshold {
            let _ = hit.set_prop("Similarity", sim, false);
            let _ = hit.set_prop(cp::_Name, name.as_str(), false);
            return true;
        }
        false
    }
}

impl SynthonSpace {
    /// `rascalSearch`
    pub fn rascal_search(
        &self,
        query: &ROMol,
        options: &RascalOptions,
        params: &SynthonSpaceSearchParams,
    ) -> Result<SearchResults> {
        check_query(query)?;
        let mut s = RascalSearcher::new(query, options, params, self)?;
        search(&mut s, ThreadMode::ThreadFragments)
    }

    /// `rascalSearch` returning the hits to a callback in chunks.
    pub fn rascal_search_incremental(
        &self,
        query: &ROMol,
        options: &RascalOptions,
        cb: SearchResultCallback<'_>,
        params: &SynthonSpaceSearchParams,
    ) -> Result<()> {
        check_query(query)?;
        let mut s = RascalSearcher::new(query, options, params, self)?;
        search_with_callback(&mut s, cb, ThreadMode::ThreadFragments)
    }

    /// `rascalSearch` on the possible hits between two lines of the possible
    /// hits file.
    pub fn rascal_search_lines(
        &self,
        query: &ROMol,
        options: &RascalOptions,
        params: &SynthonSpaceSearchParams,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SearchResults> {
        check_query(query)?;
        let _ = check_lines;
        let s = RascalSearcher::new(query, options, params, self)?;
        check_possible_hits(&s, start_line, finish_line)
    }
}
