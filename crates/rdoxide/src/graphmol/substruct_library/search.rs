//! The multithreaded search of the library (upstream `internalGetMatches`).
//!
//! Worker `t` of `n` screens the library positions `start + t, start + t + n,
//! ...`. With a `maxResults` limit every worker is allotted a share of it and
//! the workers that stopped early afterwards catch up with the most advanced
//! one, so that the final (sorted and trimmed) result does not depend on the
//! number of threads.

use std::collections::HashMap;
use std::sync::Arc;

use super::holders::{MolHolder, PatternHolder};
use crate::datastructs::ExplicitBitVect;
use crate::error::{RDKitError, Result};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::generalized_substruct::{ExtendedQueryMol, has_substruct_match_xq};
use crate::graphmol::queries::QueryKind;
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::SubstructMatchParameters;
use crate::graphmol::substruct_exists::substruct_has_match;
use crate::graphmol::tautomer_query::TautomerQuery;
use crate::rdgeneral::threads::{num_threads_to_use, par_map};

/// True if the query depends on ring information that a trusted-SMILES
/// molecule does not carry (`query_needs_rings`).
pub(super) fn query_needs_rings(query: &ROMol) -> bool {
    let atoms = query.atoms().iter().any(|a| {
        a.query().is_some_and(|q| {
            q.describe(0).contains("Ring")
                || (q.description == "RecursiveStructure"
                    && matches!(&q.kind, QueryKind::Recursive { mol, .. } if query_needs_rings(&mol.0)))
        })
    });
    atoms
        || query
            .bonds()
            .iter()
            .any(|b| b.query().is_some_and(|q| q.describe(0).contains("Ring")))
}

/// `query_needs_rings` for an extended query.
pub(super) fn xq_needs_rings(xq: &ExtendedQueryMol) -> bool {
    match xq {
        ExtendedQueryMol::Mol(m) => query_needs_rings(m),
        ExtendedQueryMol::Bundle(b) => b.mols().iter().any(|m| query_needs_rings(m)),
        ExtendedQueryMol::Tautomer(tq) => query_needs_rings(tq.template_molecule()),
        ExtendedQueryMol::TautomerBundle(v) => {
            v.iter().any(|tq| query_needs_rings(tq.template_molecule()))
        }
    }
}

/// The kinds of query the library can be searched with.
#[derive(Clone, Copy)]
pub(super) enum QueryRef<'a> {
    Mol(&'a ROMol),
    Tautomer(&'a TautomerQuery),
    Extended(&'a ExtendedQueryMol),
}

impl QueryRef<'_> {
    pub(super) fn needs_rings(&self) -> bool {
        match self {
            QueryRef::Mol(q) => query_needs_rings(q),
            QueryRef::Tautomer(tq) => query_needs_rings(tq.template_molecule()),
            QueryRef::Extended(x) => xq_needs_rings(x),
        }
    }
}

/// Everything a worker needs; shared by reference between the workers.
pub(super) struct Search<'a> {
    pub mols: &'a MolHolder,
    pub fps: Option<&'a PatternHolder>,
    pub query: QueryRef<'a>,
    /// Library positions to leave out (molecules an earlier query of a
    /// bundle already matched).
    pub skip: Option<&'a [bool]>,
    /// The query's pattern fingerprint (`None` when there is nothing to filter
    /// with).
    pub query_bits: Option<ExplicitBitVect>,
    pub params: SubstructMatchParameters,
    pub needs_rings: bool,
    pub search_order: &'a [u32],
}

/// What a worker found and the position where it stopped.
struct Outcome {
    hits: Vec<u32>,
    end: u32,
}

impl Search<'_> {
    fn passes_filter(&self, idx: u32) -> Result<bool> {
        match (self.fps, &self.query_bits) {
            (Some(fps), Some(bits)) => fps.passes_filter(idx, bits),
            _ => Ok(true),
        }
    }

    fn matches(&self, idx: u32) -> Result<bool> {
        if self.skip.is_some_and(|s| s[idx as usize]) || !self.passes_filter(idx)? {
            return Ok(false);
        }
        let Some(mol) = self.mols.get_mol(idx)? else {
            return Ok(false);
        };
        let mol: Arc<ROMol> = if self.needs_rings && !mol.ring_info().is_symm_sssr() {
            let mut m = (*mol).clone();
            symmetrize_sssr(&mut m, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
            Arc::new(m)
        } else {
            mol
        };
        Ok(match self.query {
            QueryRef::Tautomer(tq) => !tq.substruct_of(&mol, &self.params)?.0.is_empty(),
            QueryRef::Mol(q) => substruct_has_match(&mol, q, &self.params)?,
            QueryRef::Extended(x) => has_substruct_match_xq(&mol, x, &self.params)?,
        })
    }

    /// `SubSearcher`: screens positions `start, start + step, ...` below `end`;
    /// stops after `max_results` hits when that is positive, reporting the
    /// position of the last hit as the end.
    fn sub_search(&self, start: u32, mut end: u32, step: u32, max_results: i32) -> Result<Outcome> {
        if !self.search_order.is_empty() && (self.search_order.len() as u32) < end {
            return Err(RDKitError::precondition("bad searchOrder data"));
        }
        let mut hits = Vec::new();
        let mut idx = start;
        while idx < end {
            let sidx = self.search_order.get(idx as usize).copied().unwrap_or(idx);
            if self.matches(sidx)? {
                hits.push(sidx);
                if max_results > 0 && hits.len() as i32 == max_results {
                    end = idx;
                    break;
                }
            }
            idx += step;
        }
        Ok(Outcome { hits, end })
    }
}

/// Positions of the library entries in the search order.
fn positions(order: &[u32]) -> HashMap<u32, u32> {
    order
        .iter()
        .enumerate()
        .map(|(p, &v)| (v, p as u32))
        .collect()
}

/// `internalGetMatches` for a plain query: the indices of the matches (in
/// search order for one thread, sorted by search position otherwise).
pub(super) fn run(
    search: &Search,
    start_idx: u32,
    end_idx: u32,
    num_threads: i32,
    max_results: i32,
) -> Result<Vec<u32>> {
    let size = search.mols.len();
    if start_idx >= size {
        return Err(RDKitError::precondition("startIdx out of bounds"));
    }
    let order = search.search_order;
    if !order.is_empty() && start_idx as usize >= order.len() {
        return Err(RDKitError::precondition("startIdx out of bounds"));
    }
    if end_idx <= start_idx {
        return Err(RDKitError::precondition("endIdx > startIdx"));
    }
    // do not do any work if no results were requested
    if max_results == 0 {
        return Ok(Vec::new());
    }
    let mut end_idx = end_idx.min(size);
    if !order.is_empty() {
        end_idx = end_idx.min(order.len() as u32);
    }
    let n = (num_threads_to_use(num_threads) as u32).min(end_idx);
    if n <= 1 {
        return Ok(search.sub_search(start_idx, end_idx, 1, max_results)?.hits);
    }

    // split maxResults over the workers, the excess going to the first ones
    let share = if max_results > 0 {
        max_results / n as i32
    } else {
        max_results
    };
    let excess = if max_results > 0 {
        (max_results % n as i32) as u32
    } else {
        0
    };
    let workers: Vec<u32> = (0..n).collect();
    let mut outcomes = par_map(&workers, n as i32, |&t| {
        let limit = share + i32::from(t < excess);
        search.sub_search(start_idx + t, end_idx, n, limit)
    })
    .into_iter()
    .collect::<Result<Vec<_>>>()?;

    if max_results > 0 {
        // let the workers that stopped early screen as far as the most
        // productive one did
        let max_end = outcomes.iter().map(|o| o.end).max().unwrap_or(end_idx);
        let behind: Vec<(usize, u32)> = outcomes
            .iter()
            .enumerate()
            .filter(|(_, o)| o.end < max_end)
            .map(|(t, o)| (t, o.end + n))
            .collect();
        let extra = par_map(&behind, n as i32, |&(_, from)| {
            search.sub_search(from, max_end, n, -1)
        });
        for ((t, _), more) in behind.iter().zip(extra) {
            outcomes[*t].hits.extend(more?.hits);
        }
    }

    let mut idxs: Vec<u32> = outcomes.into_iter().flat_map(|o| o.hits).collect();
    // sorting makes the result independent of the number of threads
    let to_position = (!order.is_empty()).then(|| positions(order));
    if let Some(pos) = &to_position {
        idxs.iter_mut().for_each(|v| *v = pos[v]);
    }
    idxs.sort_unstable();
    if max_results > 0 {
        idxs.truncate(max_results as usize);
    }
    if to_position.is_some() {
        idxs.iter_mut().for_each(|v| *v = order[*v as usize]);
    }
    Ok(idxs)
}
