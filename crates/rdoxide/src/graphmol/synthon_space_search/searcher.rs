//! `SynthonSpaceSearcher.cpp`: the search driver shared by the searchers.

use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use super::details::{self, Stepper};
use super::hit_set::SynthonSpaceHitSet;
use super::params::{FragSetUniquifyMode, SynthonSpaceSearchParams, ThreadMode};
use super::search_results::SearchResults;
use super::space::SynthonSpace;
use super::synthon_set::SynthonSet;
use crate::graphmol::chem_transforms::mol_zip::{MolzipLabel, MolzipParams, molzip_single};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::random::Mt19937;
use crate::rdgeneral::threads::{num_threads_to_use, par_map};
use crate::{RDKitError, Result};

/// The synthon combinations to try: a hit set and an index into each of its
/// synthon lists.
pub type ToTry<'h, 'a> = Vec<(&'h SynthonSpaceHitSet<'a>, Vec<usize>)>;
/// Receives chunks of hits; returns true to stop the search.
pub type SearchResultCallback<'f> = &'f mut dyn FnMut(&mut Vec<ROMol>) -> Result<bool>;

/// State common to all the searchers.
pub struct SearcherBase<'a> {
    pub space: &'a SynthonSpace,
    pub query: ROMol,
    pub params: SynthonSpaceSearchParams,
    pub frag_set_uniquify_mode: FragSetUniquifyMode,
    rand_gen: Mutex<Option<Mt19937>>,
    best: Mutex<(f64, Option<ROMol>)>,
}

impl<'a> SearcherBase<'a> {
    pub fn new(
        query: &ROMol,
        params: &SynthonSpaceSearchParams,
        space: &'a SynthonSpace,
    ) -> Result<Self> {
        if params.random_sample && params.max_hits == -1 {
            return Err(RDKitError::Runtime(
                "Random sample is incompatible with maxHits of -1.".into(),
            ));
        }
        let rand_gen = if params.random_sample {
            let seed = if params.random_seed == -1 {
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0, |d| d.as_nanos() as u32)
            } else {
                params.random_seed as u32
            };
            Some(Mt19937::new(seed))
        } else {
            None
        };
        // For the fragmentation, it is often useful to be able to keep track
        // of the original indices.
        let mut query = query.clone();
        for i in 0..query.num_atoms() {
            query.atom_mut(i).set_prop("ORIG_IDX", i, false)?;
        }
        for i in 0..query.num_bonds() {
            query.bond_mut(i).set_prop("ORIG_IDX", i, false)?;
        }
        Ok(SearcherBase {
            space,
            query,
            params: params.clone(),
            frag_set_uniquify_mode: FragSetUniquifyMode::BySmiles,
            rand_gen: Mutex::new(rand_gen),
            best: Mutex::new((0.0, None)),
        })
    }

    /// `updateBestHitSoFar`
    pub fn update_best_hit_so_far(&self, poss_best: &ROMol, sim: f64) {
        let mut best = self.best.lock().unwrap();
        if sim > best.0 {
            best.0 = sim;
            let mut m = poss_best.clone();
            let _ = m.set_prop("Similarity", sim, false);
            best.1 = Some(m);
        }
    }
    pub fn best_similarity_so_far(&self) -> f64 {
        self.best.lock().unwrap().0
    }
    fn take_best_hit(&self) -> Option<ROMol> {
        self.best.lock().unwrap().1.take()
    }

    fn shuffle<T>(&self, v: &mut [T]) {
        let mut g = self.rand_gen.lock().unwrap();
        let g = g.as_mut().expect("random generator");
        for i in (1..v.len()).rev() {
            let j = g.uniform_u32(0, i as u32) as usize;
            v.swap(i, j);
        }
    }
}

pub trait SynthonSearcher<'a>: Sync {
    fn base(&self) -> &SearcherBase<'a>;

    /// Searches the fragment set against the synthon set.
    fn search_frag_set(
        &self,
        frag_set: &[Arc<ROMol>],
        reaction: &'a SynthonSet,
    ) -> Result<Vec<SynthonSpaceHitSet<'a>>>;

    /// An approximate similarity between the hit and query, used to order the
    /// possible hits.
    fn approx_similarity(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> f64;

    /// Checks that can be done without building the molecule.
    fn quick_verify(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> bool {
        default_quick_verify(self.base(), hitset, synth_nums)
    }

    /// Checks that the molecule is definitely a hit.
    fn verify_hit(&self, mol: &mut ROMol, _rxn_id: &str, _synth_names: &[&str]) -> bool {
        default_verify_hit(self.base(), mol)
    }

    fn build_hit(
        &self,
        hitset: &SynthonSpaceHitSet<'a>,
        synth_nums: &[usize],
        synth_names: &mut Vec<&'a str>,
    ) -> Option<ROMol> {
        let mut synths = Vec::with_capacity(synth_nums.len());
        for (i, &n) in synth_nums.iter().enumerate() {
            let (name, synthon) = hitset.synthons_to_use[i][n];
            synths.push(synthon.orig_mol());
            synth_names[i] = name;
        }
        details::build_product(&synths)
    }

    fn num_query_fragments_required(&self) -> u32 {
        self.base().space.max_num_synthons() as u32
    }

    fn extra_search_setup(
        &mut self,
        _frag_sets: &mut Vec<Vec<Arc<ROMol>>>,
        _end_time: Option<Instant>,
    ) -> Result<bool> {
        Ok(true)
    }

    fn process_to_try_set(
        &self,
        to_try: &mut ToTry<'_, 'a>,
        end_time: Option<Instant>,
        results: &mut Vec<ROMol>,
        num_hits_found: &AtomicI64,
        num_poss_hits_written: &mut u64,
    ) -> Result<()> {
        default_process_to_try_set(
            self,
            to_try,
            end_time,
            results,
            num_hits_found,
            num_poss_hits_written,
        )
    }
}

pub fn default_quick_verify(
    base: &SearcherBase<'_>,
    hitset: &SynthonSpaceHitSet<'_>,
    synth_nums: &[usize],
) -> bool {
    let p = &base.params;
    let synthon = |i: usize| hitset.synthons_to_use[i][synth_nums[i]].1;
    if p.min_hit_heavy_atoms != 0 || p.max_hit_heavy_atoms > -1 {
        let n: u32 = (0..synth_nums.len())
            .map(|i| synthon(i).num_heavy_atoms())
            .sum();
        if n < p.min_hit_heavy_atoms
            || (p.max_hit_heavy_atoms > -1 && i64::from(n) > i64::from(p.max_hit_heavy_atoms))
        {
            return false;
        }
    }
    if p.min_hit_mol_wt > 0.0 || p.max_hit_mol_wt > -1.0 {
        let wt: f64 = (0..synth_nums.len()).map(|i| synthon(i).mol_wt()).sum();
        if wt < p.min_hit_mol_wt || (p.max_hit_mol_wt > 0.0 && wt > p.max_hit_mol_wt) {
            return false;
        }
    }
    if p.min_hit_chiral_atoms != 0 || p.max_hit_chiral_atoms > -1 {
        let n: u32 = (0..synth_nums.len())
            .map(|i| synthon(i).num_chiral_atoms())
            .sum();
        let exc: u32 = (0..synth_nums.len())
            .map(|i| synthon(i).num_chiral_atoms_exc_dummies())
            .sum();
        // n is the upper bound on the number of chiral atoms in the hit,
        // exc the lower bound.
        if n < p.min_hit_chiral_atoms
            || (p.max_hit_chiral_atoms > -1 && i64::from(exc) > i64::from(p.max_hit_chiral_atoms))
        {
            return false;
        }
    }
    true
}

pub fn default_verify_hit(base: &SearcherBase<'_>, mol: &mut ROMol) -> bool {
    let p = &base.params;
    if p.min_hit_chiral_atoms != 0 || p.max_hit_chiral_atoms != 0 {
        let Ok((n, _)) = details::count_chiral_atoms(mol) else {
            return false;
        };
        if n < p.min_hit_chiral_atoms
            || (p.max_hit_chiral_atoms > -1 && i64::from(n) > i64::from(p.max_hit_chiral_atoms))
        {
            return false;
        }
    }
    true
}

/// `buildAndVerifyHit`
pub fn build_and_verify_hit<'a, S: SynthonSearcher<'a> + ?Sized>(
    s: &S,
    hitset: &SynthonSpaceHitSet<'a>,
    synth_nums: &[usize],
) -> Option<ROMol> {
    if !s.quick_verify(hitset, synth_nums) {
        return None;
    }
    let mut synth_names: Vec<&'a str> = vec![""; synth_nums.len()];
    let mut prod = s.build_hit(hitset, synth_nums, &mut synth_names)?;
    // Do a final check of the whole thing. The fragments can match the
    // synthons while the final product doesn't match (an aromatic ring
    // formed by joining the synthons, for example).
    if !s.verify_hit(&mut prod, hitset.reaction.id(), &synth_names) {
        return None;
    }
    Some(prod)
}

fn sort_hits(hits: &mut Vec<ROMol>) {
    if hits.first().is_some_and(|h| h.has_prop("Similarity")) {
        let mut keys: Vec<(usize, f64, u32)> = hits
            .iter()
            .enumerate()
            .map(|(i, h)| {
                (
                    i,
                    h.get_prop("Similarity")
                        .and_then(|v| v.as_double())
                        .unwrap_or(0.0),
                    h.num_atoms(),
                )
            })
            .collect();
        libcxx_sort_by(
            &mut keys,
            |l, r| {
                if l.1 == r.1 { l.2 < r.2 } else { l.1 > r.1 }
            },
        );
        let mut old: Vec<Option<ROMol>> = std::mem::take(hits).into_iter().map(Some).collect();
        *hits = keys.iter().map(|k| old[k.0].take().unwrap()).collect();
    }
}

fn was_enough_possible_hits(params: &SynthonSpaceSearchParams, written: u64) -> bool {
    params.write_possible_hits_and_stop && written >= params.max_possible_hits_to_write
}

fn write_possible_hits(
    file: &str,
    to_try: &ToTry<'_, '_>,
    max_to_write: u64,
    written: &mut u64,
) -> Result<()> {
    use std::io::Write;
    if *written >= max_to_write {
        return Ok(());
    }
    let mut out = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(file)
        .map_err(|e| RDKitError::Runtime(format!("Couldn't open {file}: {e}")))?;
    let mut text = String::new();
    for (hitset, nums) in to_try {
        let n = nums.len();
        for i in 0..n {
            text.push_str(hitset.synthons_to_use[i][nums[i]].1.smiles());
            if n > 1 && i < n - 1 {
                text.push('.');
            }
        }
        text.push(' ');
        text.push_str(&details::build_product_name_hs(hitset, nums));
        text.push('\n');
        *written += 1;
        if *written == max_to_write {
            break;
        }
    }
    out.write_all(text.as_bytes())
        .map_err(|e| RDKitError::Runtime(format!("Couldn't write {file}: {e}")))
}

/// `makeHitsFromToTry`
pub fn make_hits_from_to_try<'a, S: SynthonSearcher<'a> + ?Sized>(
    s: &S,
    to_try: &ToTry<'_, 'a>,
    end_time: Option<Instant>,
    results: &mut Vec<ROMol>,
    num_hits_found: &AtomicI64,
    num_poss_hits_written: &mut u64,
) -> Result<()> {
    let params = &s.base().params;
    if !params.possible_hits_file.is_empty() {
        write_possible_hits(
            &params.possible_hits_file,
            to_try,
            params.max_possible_hits_to_write,
            num_poss_hits_written,
        )?;
    }
    results.clear();
    if params.write_possible_hits_and_stop {
        return Ok(());
    }
    let max_hits = params.max_hits;
    let hit_start = params.hit_start;
    let next = AtomicUsize::new(0);
    let found: Mutex<Vec<(usize, ROMol)>> = Mutex::new(Vec::new());
    let work = || {
        let mut num_tries = 1;
        loop {
            let this_try = next.fetch_add(1, Ordering::Relaxed);
            if this_try >= to_try.len() {
                break;
            }
            if let Some(prod) = build_and_verify_hit(s, to_try[this_try].0, &to_try[this_try].1) {
                found.lock().unwrap().push((this_try, prod));
                let n = num_hits_found.fetch_add(1, Ordering::Relaxed) + 1;
                if max_hits != -1 && n >= max_hits + hit_start {
                    break;
                }
            }
            // Don't check the time every go, as it's quite expensive.
            num_tries -= 1;
            if num_tries == 0 {
                num_tries = 100;
                if details::check_time_out(end_time) {
                    break;
                }
            }
        }
    };
    let n_threads = num_threads_to_use(params.num_threads);
    if n_threads > 1 {
        std::thread::scope(|scope| {
            for _ in 0..n_threads {
                scope.spawn(work);
            }
        });
    } else {
        work();
    }
    let mut found = found.into_inner().unwrap();
    found.sort_by_key(|f| f.0);
    results.extend(found.into_iter().map(|f| f.1));
    Ok(())
}

fn sort_to_try_by_approx_similarity<'a, S: SynthonSearcher<'a> + ?Sized>(
    s: &S,
    to_try: &mut ToTry<'_, 'a>,
) {
    let mut tmp: Vec<(usize, f64)> = to_try
        .iter()
        .enumerate()
        .map(|(i, (hs, nums))| (i, s.approx_similarity(hs, nums)))
        .collect();
    libcxx_sort_by(&mut tmp, |l, r| l.1 > r.1);
    if let (Some(first), Some(last)) = (tmp.first(), tmp.last()) {
        crate::rd_log_info!(
            "Best approx similarity : {} and worst : {}",
            first.1,
            last.1
        );
    }
    let new_to_try: ToTry<'_, 'a> = tmp.iter().map(|p| to_try[p.0].clone()).collect();
    *to_try = new_to_try;
}

pub fn default_process_to_try_set<'a, S: SynthonSearcher<'a> + ?Sized>(
    s: &S,
    to_try: &mut ToTry<'_, 'a>,
    end_time: Option<Instant>,
    results: &mut Vec<ROMol>,
    num_hits_found: &AtomicI64,
    num_poss_hits_written: &mut u64,
) -> Result<()> {
    // Different fragmentations can produce overlapping synthon lists in the
    // same reaction, so remove the duplicates.
    details::sort_and_uniquify_to_try(to_try);
    if s.base().params.random_sample {
        s.base().shuffle(to_try);
    } else {
        sort_to_try_by_approx_similarity(s, to_try);
    }
    make_hits_from_to_try(
        s,
        to_try,
        end_time,
        results,
        num_hits_found,
        num_poss_hits_written,
    )
}

fn sort_hit_sets(hitsets: &mut [SynthonSpaceHitSet<'_>]) {
    hitsets.sort_by(|h1, h2| {
        if h1.reaction.id() == h2.reaction.id() {
            h1.num_hits.cmp(&h2.num_hits)
        } else {
            h1.reaction.id().cmp(h2.reaction.id())
        }
    });
}

/// `buildAllHits`
fn build_all_hits<'a, S: SynthonSearcher<'a> + ?Sized>(
    s: &S,
    hitsets: &[SynthonSpaceHitSet<'a>],
    end_time: Option<Instant>,
    timed_out: &mut bool,
    results: &mut Vec<ROMol>,
) -> Result<()> {
    let params = &s.base().params;
    if hitsets.is_empty() {
        crate::rd_log_warning!("No possible synthon matches.");
        return Ok(());
    }
    if !params.possible_hits_file.is_empty() {
        // The file will be appended to, so clear it out at the start.
        std::fs::write(&params.possible_hits_file, "").ok();
    }
    let mut to_try: ToTry<'_, 'a> = Vec::new();
    let mut enough_hits = false;
    let num_hits_found = AtomicI64::new(0);
    let max_hits = params.max_hits;
    let hit_start = params.hit_start;
    let mut num_poss_hits_written = 0u64;
    'outer: for hitset in hitsets {
        let sizes: Vec<usize> = hitset.synthons_to_use.iter().map(Vec::len).collect();
        let mut stepper = Stepper::new(&sizes);
        while !stepper.at_end() {
            to_try.push((hitset, stepper.curr_state.clone()));
            if to_try.len() as u64 == params.to_try_chunk_size {
                let mut part = Vec::new();
                s.process_to_try_set(
                    &mut to_try,
                    end_time,
                    &mut part,
                    &num_hits_found,
                    &mut num_poss_hits_written,
                )?;
                results.append(&mut part);
                enough_hits = (max_hits != -1
                    && num_hits_found.load(Ordering::Relaxed) >= max_hits + hit_start)
                    || was_enough_possible_hits(params, num_poss_hits_written);
                *timed_out = details::check_time_out(end_time);
                to_try.clear();
                if enough_hits || *timed_out {
                    break 'outer;
                }
            }
            stepper.step();
        }
    }
    // Do any remaining.
    if !enough_hits && !*timed_out && !to_try.is_empty() {
        let mut part = Vec::new();
        s.process_to_try_set(
            &mut to_try,
            end_time,
            &mut part,
            &num_hits_found,
            &mut num_poss_hits_written,
        )?;
        results.append(&mut part);
    }

    sort_hits(results);
    // The multi-threaded versions might produce more hits than requested.
    if max_hits != -1 && results.len() as i64 > max_hits + hit_start {
        results.truncate((max_hits + hit_start) as usize);
    }
    // Now get rid of any hits before hit_start.
    if hit_start != 0 {
        if (hit_start as usize) < results.len() {
            results.drain(..hit_start as usize);
        } else {
            results.clear();
        }
    }
    Ok(())
}

fn build_hits<'a, S: SynthonSearcher<'a> + ?Sized>(
    s: &S,
    hitsets: &mut [SynthonSpaceHitSet<'a>],
    end_time: Option<Instant>,
    timed_out: &mut bool,
    results: &mut Vec<ROMol>,
) -> Result<()> {
    if hitsets.is_empty() {
        return Ok(());
    }
    if s.base().params.random_sample {
        s.base().shuffle(hitsets);
    } else {
        sort_hit_sets(hitsets);
    }
    build_all_hits(s, hitsets, end_time, timed_out, results)
}

fn search_reaction<'a, S: SynthonSearcher<'a>>(
    s: &S,
    reaction: &'a SynthonSet,
    end_time: Option<Instant>,
    num_threads: i32,
    fragments: &[Vec<Arc<ROMol>>],
) -> Result<Vec<SynthonSpaceHitSet<'a>>> {
    let timed_out = AtomicBool::new(false);
    let count = AtomicUsize::new(0);
    let idx: Vec<usize> = (0..fragments.len()).collect();
    let per_frag: Vec<Result<Vec<SynthonSpaceHitSet<'a>>>> = par_map(&idx, num_threads, |&i| {
        if timed_out.load(Ordering::Relaxed) {
            return Ok(Vec::new());
        }
        if count.fetch_add(1, Ordering::Relaxed) % 100 == 0 && details::check_time_out(end_time) {
            timed_out.store(true, Ordering::Relaxed);
            return Ok(Vec::new());
        }
        s.search_frag_set(&fragments[i], reaction)
    });
    let mut hits = Vec::new();
    for r in per_frag {
        hits.extend(r?);
    }
    Ok(hits)
}

fn do_the_search<'a, S: SynthonSearcher<'a>>(
    s: &S,
    frag_sets: &[Vec<Arc<ROMol>>],
    end_time: Option<Instant>,
    timed_out: &mut bool,
    tot_hits: &mut u64,
    thread_mode: ThreadMode,
) -> Result<Vec<SynthonSpaceHitSet<'a>>> {
    let base = s.base();
    let space: &'a SynthonSpace = base.space;
    let reaction_hits: Vec<Vec<SynthonSpaceHitSet<'a>>> = match thread_mode {
        ThreadMode::ThreadFragments => {
            // Do the reactions one at a time, parallelising the fragSet
            // searches.
            let mut v = Vec::new();
            for rxn in space.reactions() {
                v.push(search_reaction(
                    s,
                    rxn,
                    end_time,
                    base.params.num_threads,
                    frag_sets,
                )?);
            }
            v
        }
        ThreadMode::ThreadReactions => {
            // Parallelise at the reaction level.
            let rxns: Vec<&'a SynthonSet> = space.reactions().iter().collect();
            let stop = AtomicBool::new(false);
            let res: Vec<Result<Vec<SynthonSpaceHitSet<'a>>>> =
                par_map(&rxns, base.params.num_threads, |rxn| {
                    if stop.load(Ordering::Relaxed) || details::check_time_out(end_time) {
                        stop.store(true, Ordering::Relaxed);
                        return Ok(Vec::new());
                    }
                    search_reaction(s, rxn, end_time, 1, frag_sets)
                });
            res.into_iter().collect::<Result<Vec<_>>>()?
        }
    };
    let mut all_hits = Vec::new();
    *tot_hits = 0;
    for rh in reaction_hits {
        *tot_hits += rh.iter().map(|h| h.num_hits as u64).sum::<u64>();
        all_hits.extend(rh);
    }
    *timed_out = details::check_time_out(end_time);
    Ok(all_hits)
}

fn assemble_hit_sets<'a, S: SynthonSearcher<'a>>(
    s: &mut S,
    end_time: Option<Instant>,
    timed_out: &mut bool,
    tot_hits: &mut u64,
    thread_mode: ThreadMode,
) -> Result<Vec<SynthonSpaceHitSet<'a>>> {
    let (query, max_frag_sets, num_threads, mode) = {
        let b = s.base();
        (
            b.query.clone(),
            b.params.max_num_frag_sets,
            b.params.num_threads,
            b.frag_set_uniquify_mode,
        )
    };
    let n_required = s.num_query_fragments_required();
    let mut fragments = details::split_molecule(
        &query,
        n_required,
        max_frag_sets,
        end_time,
        num_threads,
        mode,
        timed_out,
    )?;
    if *timed_out {
        return Ok(Vec::new());
    }
    if !s.extra_search_setup(&mut fragments, end_time)? {
        return Ok(Vec::new());
    }
    do_the_search(&*s, &fragments, end_time, timed_out, tot_hits, thread_mode)
}

/// `SynthonSpaceSearcher::search(threadMode)`
pub fn search<'a, S: SynthonSearcher<'a>>(
    s: &mut S,
    thread_mode: ThreadMode,
) -> Result<SearchResults> {
    let mut results = Vec::new();
    let end_time = (s.base().params.time_out > 0)
        .then(|| Instant::now() + Duration::from_secs(s.base().params.time_out));
    let mut timed_out = false;
    let mut tot_hits = 0;
    let mut all_hits = assemble_hit_sets(s, end_time, &mut timed_out, &mut tot_hits, thread_mode)?;
    if !timed_out && s.base().params.build_hits {
        build_hits(&*s, &mut all_hits, end_time, &mut timed_out, &mut results)?;
    }
    Ok(SearchResults::with(
        results,
        s.base().take_best_hit(),
        tot_hits,
        timed_out,
        false,
    ))
}

/// `SynthonSpaceSearcher::search(cb, threadMode)`
pub fn search_with_callback<'a, S: SynthonSearcher<'a>>(
    s: &mut S,
    cb: SearchResultCallback<'_>,
    thread_mode: ThreadMode,
) -> Result<()> {
    let mut timed_out = false;
    let mut tot_hits = 0;
    let mut all_hits = assemble_hit_sets(s, None, &mut timed_out, &mut tot_hits, thread_mode)?;
    sort_hit_sets(&mut all_hits);
    let params = s.base().params.clone();
    let mut to_try: ToTry<'_, 'a> = Vec::new();
    let num_hits_found = AtomicI64::new(0);
    let mut num_poss_hits_found = 0u64;
    let mut stop = false;
    let enough =
        |n: &AtomicI64| params.max_hits != -1 && n.load(Ordering::Relaxed) >= params.max_hits;
    'outer: for hitset in &all_hits {
        let sizes: Vec<usize> = hitset.synthons_to_use.iter().map(Vec::len).collect();
        let mut stepper = Stepper::new(&sizes);
        while !stepper.at_end() {
            to_try.push((hitset, stepper.curr_state.clone()));
            if to_try.len() as u64 == params.to_try_chunk_size {
                let mut part = Vec::new();
                s.process_to_try_set(
                    &mut to_try,
                    None,
                    &mut part,
                    &num_hits_found,
                    &mut num_poss_hits_found,
                )?;
                stop = cb(&mut part)?;
                to_try.clear();
                if stop || enough(&num_hits_found) {
                    break 'outer;
                }
            }
            stepper.step();
        }
    }
    // Do any remaining.
    if !stop && !enough(&num_hits_found) && !to_try.is_empty() {
        let mut part = Vec::new();
        s.process_to_try_set(
            &mut to_try,
            None,
            &mut part,
            &num_hits_found,
            &mut num_poss_hits_found,
        )?;
        cb(&mut part)?;
    }
    Ok(())
}

fn read_poss_hits_lines(file: &str, start: u64, finish: u64) -> Vec<(String, String)> {
    let Ok(text) = std::fs::read_to_string(file) else {
        crate::rd_log_error!("Possible hits file {file} not found, so no results generated.");
        return Vec::new();
    };
    let mut lines = text.split('\n');
    for _ in 0..start {
        if lines.next().is_none() {
            return Vec::new();
        }
    }
    let mut ret = Vec::new();
    let mut i = start;
    while i < finish {
        let Some(line) = lines.next() else { break };
        i += 1;
        if line.is_empty() {
            continue;
        }
        let line = line.strip_suffix('\r').unwrap_or(line);
        let mut bits = line.split(['\t', ' ']).filter(|b| !b.is_empty());
        let Some(smiles) = bits.next() else { continue };
        ret.push((smiles.to_string(), bits.next().unwrap_or("").to_string()));
    }
    ret
}

/// `SynthonSpaceSearcher::checkPossibleHits`: builds the possible hits in
/// `possibleHitsFile` between the given lines and returns those that match.
pub fn check_possible_hits<'a, S: SynthonSearcher<'a>>(
    s: &S,
    start_line: u64,
    finish_line: u64,
) -> Result<SearchResults> {
    let params = &s.base().params;
    if params.possible_hits_file.is_empty() {
        crate::rd_log_warning!("No possible hits file given, so no results generated.");
        return Ok(SearchResults::new());
    }
    let check_lines = read_poss_hits_lines(&params.possible_hits_file, start_line, finish_line);
    if check_lines.is_empty() {
        return Ok(SearchResults::new());
    }
    let mz = MolzipParams {
        label: MolzipLabel::Isotope,
        align_coordinates: false,
        ..Default::default()
    };
    let result_mols: Vec<Option<ROMol>> =
        par_map(&check_lines, params.num_threads, |(smiles, name)| {
            let mol = mol_from_smiles(smiles, &SmilesParserParams::default())
                .ok()
                .flatten()?;
            let prod = (|| -> Result<ROMol> {
                let mut zip = molzip_single(&mol, &mz)?
                    .ok_or_else(|| RDKitError::Value("molzip failed".into()))?;
                sanitize_mol(&mut zip, SANITIZE_ALL).map_err(|(_, e)| e)?;
                Ok(zip)
            })();
            let mut prod = match prod {
                Ok(p) => p,
                Err(e) => {
                    println!(
                        "Failed to zip {} for {name} because \n{e}",
                        details::smiles_of(&mol)
                    );
                    return None;
                }
            };
            let _ = prod.set_prop(cp::_Name, name.as_str(), false);
            let mut synth_names: Vec<&str> = name.split(';').collect();
            let rxn_id = synth_names.pop().unwrap_or("InvalidName");
            s.verify_hit(&mut prod, rxn_id, &synth_names)
                .then_some(prod)
        });
    let mut result_mols: Vec<ROMol> = result_mols.into_iter().flatten().collect();
    sort_hits(&mut result_mols);
    Ok(SearchResults::with(
        result_mols,
        s.base().take_best_hit(),
        check_lines.len() as u64,
        false,
        false,
    ))
}
