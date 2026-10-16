//! `SynthonSpaceSubstructureSearcher.cpp`

use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

use super::details::{self, ConnectorBits};
use super::hit_set::SynthonSpaceHitSet;
use super::params::{SynthonSpaceSearchParams, ThreadMode};
use super::search_results::SearchResults;
use super::searcher::{
    SearchResultCallback, SearcherBase, SynthonSearcher, check_possible_hits, search,
    search_with_callback,
};
use super::space::SynthonSpace;
use super::synthon::PATT_FP_NUM_BITS;
use super::synthon_set::SynthonSet;
use crate::Result;
use crate::datastructs::bit_ops::all_probe_bits_match;
use crate::datastructs::bit_vect::{BitVect, ExplicitBitVect};
use crate::geometry::Point3D;
use crate::graphmol::chem_transforms::combine_mols;
use crate::graphmol::fingerprints::pattern::pattern_fingerprint_mol;
use crate::graphmol::generalized_substruct::ExtendedQueryMol;
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::queries::{CompositeQueryType, make_atom_isotope_query, make_atom_type_query};
use crate::graphmol::query_ops::atom_expand_query;
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;

/// The per-fragment information used to screen synthons.
struct FragInfo {
    patt_fp: ExplicitBitVect,
    conn_regs: Vec<ROMol>,
    conn_reg_smis: Vec<String>,
    conn_reg_fps: Vec<ExplicitBitVect>,
}

pub struct SubstructureSearcher<'a> {
    base: SearcherBase<'a>,
    match_params: SubstructMatchParameters,
    /// keyed by the address of the fragment
    frag_info: HashMap<usize, Arc<FragInfo>>,
}

fn key(frag: &Arc<ROMol>) -> usize {
    Arc::as_ptr(frag) as usize
}

fn has_match(mol: &ROMol, query: &ROMol, ps: &SubstructMatchParameters) -> bool {
    substruct_match(mol, query, ps).is_ok_and(|m| !m.is_empty())
}

fn merge_fragments_2(
    i: usize,
    one: usize,
    two: usize,
    others: &[usize],
    frag_set: &[Arc<ROMol>],
) -> Result<Vec<Arc<ROMol>>> {
    let mut ret = vec![Arc::new((*frag_set[i]).clone())];
    ret.push(Arc::new(combine_mols(
        &frag_set[one],
        &frag_set[two],
        Point3D::default(),
    )?));
    for &o in others {
        if o != one && o != two {
            ret.push(Arc::new((*frag_set[o]).clone()));
        }
    }
    Ok(ret)
}

fn merge_fragments_4(
    i: usize,
    one: usize,
    two: usize,
    three: usize,
    four: usize,
    frag_set: &[Arc<ROMol>],
) -> Result<Vec<Arc<ROMol>>> {
    Ok(vec![
        Arc::new((*frag_set[i]).clone()),
        Arc::new(combine_mols(
            &frag_set[one],
            &frag_set[two],
            Point3D::default(),
        )?),
        Arc::new(combine_mols(
            &frag_set[three],
            &frag_set[four],
            Point3D::default(),
        )?),
    ])
}

fn merge_ring_forming_frags(frag_set: &[Arc<ROMol>]) -> Result<Vec<Vec<Arc<ROMol>>>> {
    let mut out = Vec::new();
    // There can be 1 or 2 ring closures in a reaction, so either 2 of the
    // synthons will be bidentate (1 ring closure), or 1 synthon will be
    // 4-dentate.
    let conn_patts = details::get_connector_patterns(frag_set);
    for i in 0..conn_patts.len() {
        if conn_patts[i].count_ones() > 1 {
            let others: Vec<usize> = (0..conn_patts.len()).filter(|&j| j != i).collect();
            match others.len() {
                2 => out.push(merge_fragments_2(
                    i, others[0], others[1], &others, frag_set,
                )?),
                3 => {
                    out.push(merge_fragments_2(
                        i, others[0], others[1], &others, frag_set,
                    )?);
                    out.push(merge_fragments_2(
                        i, others[0], others[2], &others, frag_set,
                    )?);
                    out.push(merge_fragments_2(
                        i, others[1], others[2], &others, frag_set,
                    )?);
                }
                4 => {
                    out.push(merge_fragments_4(
                        i, others[0], others[1], others[2], others[3], frag_set,
                    )?);
                    out.push(merge_fragments_4(
                        i, others[0], others[2], others[1], others[3], frag_set,
                    )?);
                    out.push(merge_fragments_4(
                        i, others[1], others[2], others[0], others[3], frag_set,
                    )?);
                }
                _ => {}
            }
        }
    }
    Ok(out)
}

fn screen_synthons_with_fps(
    patt_fps: &[&ExplicitBitVect],
    reaction: &SynthonSet,
    pool: &[super::synthon::Synthon],
    synthon_order: &[usize],
) -> Vec<Vec<bool>> {
    let mut passed: Vec<Vec<bool>> = reaction
        .synthons()
        .iter()
        .map(|s| vec![false; s.len()])
        .collect();
    let mut frags_matched = vec![false; synthon_order.len()];
    for (i, &so) in synthon_order.iter().enumerate() {
        let set = &reaction.synthons()[so];
        for (j, (_, idx)) in set.iter().enumerate() {
            let sfp = pool[*idx].patt_fp();
            if patt_fps[i].num_on_bits() <= sfp.num_on_bits()
                && all_probe_bits_match(patt_fps[i], sfp).unwrap_or(false)
            {
                passed[so][j] = true;
                frags_matched[i] = true;
            }
        }
        // If nothing matched this fragment, the whole thing's a bust.
        if !frags_matched[i] {
            break;
        }
    }
    if !frags_matched.iter().all(|&b| b) {
        for p in &mut passed {
            p.iter_mut().for_each(|b| *b = false);
        }
    }
    passed
}

fn find_synthon_search_start(
    frag_num_atoms: u32,
    set_num: usize,
    reaction: &SynthonSet,
    pool: &[super::synthon::Synthon],
) -> usize {
    let atoms_at = |n: usize| {
        let o = reaction
            .substructure_ordered_synthon_num(set_num, n)
            .unwrap();
        pool[reaction.synthons()[set_num][o].1]
            .search_mol()
            .num_atoms()
    };
    let mut first = 0;
    if frag_num_atoms <= atoms_at(first) {
        return first;
    }
    // binary_search_leftmost
    let mut last = reaction.synthons()[set_num].len();
    while first < last {
        let mid = first + (last - first) / 2;
        if atoms_at(mid) < frag_num_atoms {
            first = mid + 1;
        } else {
            last = mid;
        }
    }
    first
}

fn get_hit_synthons(
    mol_frags: &[ROMol],
    passed_screens: &[Vec<bool>],
    reaction: &SynthonSet,
    pool: &[super::synthon::Synthon],
    synthon_set_order: &[usize],
    ps: &SubstructMatchParameters,
) -> Vec<Vec<usize>> {
    // Match the fragments against the synthon sets, the smallest set first,
    // because if a fragment doesn't match anything the whole thing's a bust.
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
    for (frag_num, _) in &frag_orders {
        let frag_num = *frag_num;
        let so = synthon_set_order[frag_num];
        let set = &reaction.synthons()[so];
        let passed = &passed_screens[so];
        let mut frag_matched = false;
        let start = find_synthon_search_start(mol_frags[frag_num].num_atoms(), so, reaction, pool);
        for j in start..set.len() {
            let synthon_num = reaction.substructure_ordered_synthon_num(so, j).unwrap();
            if passed[synthon_num]
                && has_match(
                    pool[set[synthon_num].1].search_mol(),
                    &mol_frags[frag_num],
                    ps,
                )
            {
                synthons_to_use[so][synthon_num] = true;
                frag_matched = true;
            }
        }
        if !frag_matched {
            return Vec::new();
        }
    }
    details::expand_bit_set(&mut synthons_to_use);
    let mut ret = details::bit_sets_to_vectors(&synthons_to_use);
    // sort the selected synthons into ascending order of number of atoms
    for (i, r) in ret.iter_mut().enumerate() {
        let synthons_i = &reaction.synthons()[i];
        libcxx_sort_by(r, |&a, &b| {
            pool[synthons_i[a].1].orig_mol().num_atoms()
                < pool[synthons_i[b].1].orig_mol().num_atoms()
        });
    }
    ret
}

impl<'a> SubstructureSearcher<'a> {
    pub fn new(
        query: &ROMol,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
        space: &'a SynthonSpace,
    ) -> Result<Self> {
        let mut ps = match_params.clone();
        // only the existence of a match matters when screening
        ps.max_matches = 1;
        Ok(SubstructureSearcher {
            base: SearcherBase::new(query, params, space)?,
            match_params: ps,
            frag_info: HashMap::new(),
        })
    }

    fn info(&self, frag: &Arc<ROMol>) -> &FragInfo {
        &self.frag_info[&key(frag)]
    }

    fn check_connector_regions(&self, reaction: &SynthonSet, frag_set: &[Arc<ROMol>]) -> bool {
        let rxn_regs = reaction.connector_regions();
        let rxn_smis = reaction.connector_region_smiles();
        let rxn_fps = reaction.conn_reg_fps();
        for frag in frag_set {
            let info = self.info(frag);
            if info.conn_regs.is_empty() {
                continue;
            }
            let mut found = false;
            'regs: for j in 0..info.conn_regs.len() {
                for k in 0..rxn_regs.len() {
                    if rxn_smis[k] == info.conn_reg_smis[j] {
                        found = true;
                        break 'regs;
                    }
                    if info.conn_reg_fps[j].num_on_bits() <= rxn_fps[k].num_on_bits()
                        && all_probe_bits_match(&info.conn_reg_fps[j], &rxn_fps[k]).unwrap_or(false)
                        && has_match(
                            &rxn_regs[k],
                            &info.conn_regs[j],
                            &SubstructMatchParameters {
                                max_matches: 1,
                                ..Default::default()
                            },
                        )
                    {
                        found = true;
                        break 'regs;
                    }
                }
            }
            if !found {
                return false;
            }
        }
        true
    }
}

impl<'a> SynthonSearcher<'a> for SubstructureSearcher<'a> {
    fn base(&self) -> &SearcherBase<'a> {
        &self.base
    }

    fn num_query_fragments_required(&self) -> u32 {
        let max_synthons = self.base.space.max_num_synthons() as u32;
        let max_conns = self.base.space.max_num_connectors();
        max_synthons.max(max_conns + 1)
    }

    fn extra_search_setup(
        &mut self,
        frag_sets: &mut Vec<Vec<Arc<ROMol>>>,
        end_time: Option<Instant>,
    ) -> Result<bool> {
        if self.base.space.has_ring_former() {
            // If there is a reaction that has a pair of synthons that form a
            // ring, extra merged fragments are required.
            let mut extra = Vec::new();
            for frag_set in frag_sets.iter() {
                // (for a fragSet of size 2 both fragments will be bidentate
                // so might match, or at least one of them isn't)
                if frag_set.len() > 2 {
                    extra.extend(merge_ring_forming_frags(frag_set)?);
                }
            }
            frag_sets.extend(extra);
        }
        let map = details::map_frags_by_smiles(frag_sets)?;
        if details::check_time_out(end_time) {
            return Ok(false);
        }
        let patt_fp_size = self.base.space.pattern_fp_size()?;
        let mut said_something = false;
        for (_smi, members) in map {
            // Query atoms may define the environment of the fragment that a
            // potentially matching synthon may not have, so they need to be
            // made more generic.
            for &(i, j) in &members {
                let frag = Arc::get_mut(&mut frag_sets[i][j]).expect("unshared fragment");
                if details::remove_query_atoms(frag)? && !said_something {
                    said_something = true;
                    crate::rd_log_info!("Complex queries can be slow.");
                }
            }
            let first = &frag_sets[members[0].0][members[0].1];
            let patt_fp = pattern_fingerprint_mol(first, patt_fp_size, false)?;
            let mut info = FragInfo {
                patt_fp,
                conn_regs: Vec::new(),
                conn_reg_smis: Vec::new(),
                conn_reg_fps: Vec::new(),
            };
            if let Some(conn_region) = details::build_conn_region(first)? {
                let (regs, _) = get_mol_frags_as_mols(&conn_region, false, true)?;
                for cr in regs {
                    info.conn_reg_fps
                        .push(pattern_fingerprint_mol(&cr, PATT_FP_NUM_BITS, false)?);
                    info.conn_reg_smis.push(details::smiles_of(&cr));
                    info.conn_regs.push(cr);
                }
            }
            let info = Arc::new(info);
            for &(i, j) in &members {
                self.frag_info.insert(key(&frag_sets[i][j]), info.clone());
            }
        }
        // Now apply the largest fragment heuristic to the fragments.
        for fs in frag_sets.iter_mut() {
            let mut fps: Vec<(usize, u32)> = fs
                .iter()
                .enumerate()
                .map(|(i, f)| (i, self.info(f).patt_fp.num_on_bits()))
                .collect();
            libcxx_sort_by(&mut fps, |a, b| a.1 > b.1);
            let new_frags: Vec<Arc<ROMol>> = fps.iter().map(|p| fs[p.0].clone()).collect();
            *fs = new_frags;
        }
        Ok(true)
    }

    fn search_frag_set(
        &self,
        frag_set: &[Arc<ROMol>],
        reaction: &'a SynthonSet,
    ) -> Result<Vec<SynthonSpaceHitSet<'a>>> {
        let mut results = Vec::new();
        let pool = self.base.space.pool();
        let patt_fps: Vec<&ExplicitBitVect> =
            frag_set.iter().map(|f| &self.info(f).patt_fp).collect();
        let conns: ConnectorBits = details::get_connector_pattern(frag_set);
        // It can't be a hit if the number of fragments is more than 1 plus
        // the number of bonds being formed in the reaction.
        if frag_set.len() > reaction.connectors().count_ones() as usize + 1 {
            return Ok(results);
        }
        // Check that all the frags have a connector region that matches
        // something in this reaction set.
        if !self.check_connector_regions(reaction, frag_set) {
            return Ok(results);
        }
        let mut frag_set_cp: Vec<ROMol> = frag_set.iter().map(|f| (**f).clone()).collect();

        // All the possible permutations of connector numbers compatible with
        // the number of synthon sets in this reaction.
        let conn_combs =
            details::get_connector_permutations(&frag_set_cp, conns, reaction.connectors());

        // Select only the synthons that have fingerprints that are a
        // superset of the fragment fingerprints; try all synthon orders.
        let synthon_set_orders =
            details::perm_m_from_n(frag_set_cp.len(), reaction.synthons().len());
        for sso in &synthon_set_orders {
            let passed_screens = screen_synthons_with_fps(&patt_fps, reaction, pool, sso);
            if passed_screens.iter().all(|s| !s.iter().any(|&b| b)) {
                continue;
            }
            // Find all synthons that match the fragments with each
            // connector combination.
            for conn_comb in &conn_combs {
                for (i, assignments) in conn_comb.iter().enumerate() {
                    for &(atom_idx, isotope) in assignments {
                        let atom = frag_set_cp[i].atom_mut(atom_idx);
                        atom.set_isotope(isotope);
                        if atom.has_query() {
                            atom.set_query(Some(make_atom_type_query(0, false)));
                            atom_expand_query(
                                atom,
                                make_atom_isotope_query(isotope as i32),
                                CompositeQueryType::And,
                                true,
                            )?;
                        }
                    }
                }
                let these = get_hit_synthons(
                    &frag_set_cp,
                    &passed_screens,
                    reaction,
                    pool,
                    sso,
                    &self.match_params,
                );
                if !these.is_empty() {
                    let hs = SynthonSpaceHitSet::new(reaction, pool, &these, frag_set);
                    if hs.num_hits > 0 {
                        results.push(hs);
                    }
                }
            }
        }
        Ok(results)
    }

    fn approx_similarity(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> f64 {
        let qbit = f64::from(self.base.query.num_atoms());
        let num_atoms: u32 = synth_nums
            .iter()
            .enumerate()
            .map(|(i, &n)| hitset.synthons_to_use[i][n].1.num_heavy_atoms())
            .sum();
        qbit / f64::from(num_atoms)
    }

    fn verify_hit(&self, hit: &mut ROMol, rxn_id: &str, synth_names: &[&str]) -> bool {
        if !super::searcher::default_verify_hit(&self.base, hit) {
            return false;
        }
        if has_match(hit, &self.base.query, &self.match_params) {
            let name = details::build_product_name(rxn_id, synth_names);
            let _ = hit.set_prop(cp::_Name, name, false);
            return true;
        }
        false
    }
}

impl SynthonSpace {
    pub fn reactions(&self) -> &[SynthonSet] {
        &self.reactions
    }

    /// `substructureSearch`
    pub fn substructure_search(
        &self,
        query: &ROMol,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
    ) -> Result<SearchResults> {
        if query.num_atoms() == 0 {
            return Err(crate::RDKitError::precondition(
                "Search query must contain atoms.",
            ));
        }
        let mut s = SubstructureSearcher::new(query, match_params, params, self)?;
        search(&mut s, ThreadMode::ThreadReactions)
    }

    /// `substructureSearch` returning the hits to a callback in chunks.
    pub fn substructure_search_incremental(
        &self,
        query: &ROMol,
        cb: SearchResultCallback<'_>,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
    ) -> Result<()> {
        if query.num_atoms() == 0 {
            return Err(crate::RDKitError::precondition(
                "Search query must contain atoms.",
            ));
        }
        let mut s = SubstructureSearcher::new(query, match_params, params, self)?;
        search_with_callback(&mut s, cb, ThreadMode::ThreadReactions)
    }

    /// `substructureSearch` on the possible hits between two lines of the
    /// possible hits file.
    pub fn substructure_search_lines(
        &self,
        query: &ROMol,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SearchResults> {
        if query.num_atoms() == 0 {
            return Err(crate::RDKitError::precondition(
                "Search query must contain atoms.",
            ));
        }
        if start_line >= finish_line {
            return Err(crate::RDKitError::precondition(
                "Search startLine must be less than finishLine.",
            ));
        }
        let s = SubstructureSearcher::new(query, match_params, params, self)?;
        check_possible_hits(&s, start_line, finish_line)
    }
}

impl SynthonSpace {
    /// `substructureSearch` with an extended query (a molecule, bundle,
    /// tautomer query or bundle of tautomer queries).
    pub fn substructure_search_xq(
        &self,
        query: &ExtendedQueryMol,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
    ) -> Result<SearchResults> {
        self.extended_search(query, match_params, params, None)
    }

    /// `substructureSearch` with an extended query on the possible hits
    /// between two lines of the possible hits file.
    pub fn substructure_search_xq_lines(
        &self,
        query: &ExtendedQueryMol,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SearchResults> {
        if start_line >= finish_line {
            return Err(crate::RDKitError::precondition(
                "Search startLine must be less than finishLine.",
            ));
        }
        self.extended_search(query, match_params, params, Some((start_line, finish_line)))
    }

    fn search_one(
        &self,
        mol: &ROMol,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
        lines: Option<(u64, u64)>,
    ) -> Result<SearchResults> {
        match lines {
            None => self.substructure_search(mol, match_params, params),
            Some((s, f)) => self.substructure_search_lines(mol, match_params, params, s, f),
        }
    }

    fn search_many<'m>(
        &self,
        mols: impl Iterator<Item = &'m ROMol>,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
        lines: Option<(u64, u64)>,
    ) -> Result<SearchResults> {
        let mut results = SearchResults::new();
        let mut tmp_params = params.clone();
        for m in mols {
            let mut these = self.search_one(m, match_params, &tmp_params, lines)?;
            if tmp_params.max_hits != -1 {
                tmp_params.max_hits -= these.hit_molecules().len() as i64;
            }
            results.merge_results(&mut these);
        }
        Ok(results)
    }

    fn extended_search(
        &self,
        query: &ExtendedQueryMol,
        match_params: &SubstructMatchParameters,
        params: &SynthonSpaceSearchParams,
        lines: Option<(u64, u64)>,
    ) -> Result<SearchResults> {
        match query {
            ExtendedQueryMol::Mol(m) => self.search_one(m, match_params, params, lines),
            ExtendedQueryMol::Bundle(b) => {
                self.search_many(b.mols().iter().map(|m| &**m), match_params, params, lines)
            }
            ExtendedQueryMol::Tautomer(tq) => self.search_many(
                tq.tautomers().iter().map(|m| &**m),
                match_params,
                params,
                lines,
            ),
            ExtendedQueryMol::TautomerBundle(tqs) => {
                let mut results = SearchResults::new();
                let mut tmp_params = params.clone();
                for tq in tqs {
                    let mut these = self.search_many(
                        tq.tautomers().iter().map(|m| &**m),
                        match_params,
                        &tmp_params,
                        lines,
                    )?;
                    if tmp_params.max_hits != -1 {
                        tmp_params.max_hits -= these.hit_molecules().len() as i64;
                    }
                    results.merge_results(&mut these);
                }
                Ok(results)
            }
        }
    }
}
