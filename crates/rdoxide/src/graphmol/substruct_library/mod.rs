//! Multithreaded substructure search of a molecule library (upstream
//! `GraphMol/SubstructLibrary`).
//!
//! A [`SubstructLibrary`] couples a [`MolHolder`] with an optional
//! [`PatternHolder`] (pattern fingerprints that rule out impossible matches
//! cheaply) and an optional [`KeyFromPropHolder`] (user keys). Searches are
//! spread over worker threads and return the same indices for any thread
//! count.
//!
//! Not ported: boost-serialization pickles; libraries
//! serialize to this crate's own format ([`SubstructLibrary::to_bytes`]).

mod holders;
mod search;
mod serialization;

use std::sync::{Arc, PoisonError, RwLock, RwLockReadGuard};

use search::QueryRef;

pub use holders::{KeyFromPropHolder, MolHolder, PatternHolder};

use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::error::{RDKitError, Result};
use crate::graphmol::generalized_substruct::ExtendedQueryMol;
use crate::graphmol::mol_bundle::MolBundle;
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::SubstructMatchParameters;
use crate::graphmol::tautomer_query::TautomerQuery;
use crate::rd_log_warning;
use crate::rdgeneral::threads::par_map;

/// A holder shared between a library and the code that built it, so that
/// molecules added through either are seen by both (as with upstream's
/// `shared_ptr` holders).
pub type Shared<T> = Arc<RwLock<T>>;

pub fn shared<T>(value: T) -> Shared<T> {
    Arc::new(RwLock::new(value))
}

fn read<T>(holder: &Shared<T>) -> RwLockReadGuard<'_, T> {
    holder.read().unwrap_or_else(PoisonError::into_inner)
}

/// `SubstructLibrary`
#[derive(Debug, Clone)]
pub struct SubstructLibrary {
    mols: Shared<MolHolder>,
    fps: Option<Shared<PatternHolder>>,
    keys: Option<Shared<KeyFromPropHolder>>,
    search_order: Vec<u32>,
}

impl Default for SubstructLibrary {
    fn default() -> Self {
        Self::new(shared(MolHolder::default()), None, None)
    }
}

impl SubstructLibrary {
    pub fn new(
        mols: Shared<MolHolder>,
        fps: Option<Shared<PatternHolder>>,
        keys: Option<Shared<KeyFromPropHolder>>,
    ) -> Self {
        SubstructLibrary {
            mols,
            fps,
            keys,
            search_order: Vec::new(),
        }
    }

    pub fn mol_holder(&self) -> &Shared<MolHolder> {
        &self.mols
    }

    pub fn fp_holder(&self) -> Option<&Shared<PatternHolder>> {
        self.fps.as_ref()
    }

    pub fn key_holder(&self) -> Option<&Shared<KeyFromPropHolder>> {
        self.keys.as_ref()
    }

    /// Adds a molecule (and its fingerprint and key) and returns its index.
    pub fn add_mol(&self, mol: &ROMol) -> Result<u32> {
        let idx = self
            .mols
            .write()
            .unwrap_or_else(PoisonError::into_inner)
            .add_mol(mol)?;
        if let Some(fps) = &self.fps {
            let fp_idx = fps
                .write()
                .unwrap_or_else(PoisonError::into_inner)
                .add_mol(mol)?;
            if fp_idx != idx {
                return Err(RDKitError::invariant(
                    "#mols different than #fingerprints in SubstructLibrary",
                ));
            }
        }
        if let Some(keys) = &self.keys {
            let key_idx = keys
                .write()
                .unwrap_or_else(PoisonError::into_inner)
                .add_mol(mol);
            if key_idx != idx {
                return Err(RDKitError::invariant(
                    "#mols different than #keys in SubstructLibrary",
                ));
            }
        }
        Ok(idx)
    }

    /// The molecule at `idx` (`IndexError` when out of range, `None` for a
    /// stored SMILES that does not parse).
    pub fn get_mol(&self, idx: u32) -> Result<Option<Arc<ROMol>>> {
        read(&self.mols).get_mol(idx)
    }

    pub fn len(&self) -> u32 {
        read(&self.mols).len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Sets the order in which molecules are screened (an `IndexError` if an
    /// index is out of range); an empty order restores the natural one.
    pub fn set_search_order(&mut self, order: Vec<u32>) -> Result<()> {
        let size = self.len();
        if let Some(&bad) = order.iter().find(|&&i| i >= size) {
            return Err(RDKitError::Index(bad.to_string()));
        }
        self.search_order = order;
        Ok(())
    }

    pub fn search_order(&self) -> &[u32] {
        &self.search_order
    }

    /// `getMatches`: the indices of the molecules between `start_idx` and
    /// `end_idx` (exclusive) that contain `query`, at most `max_results` of
    /// them when that is non-negative (-1: all). Uses `num_threads` workers
    /// (<= 0: all cores minus `|num_threads|`).
    pub fn get_matches(
        &self,
        query: &ROMol,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
        max_results: i32,
    ) -> Result<Vec<u32>> {
        self.search(
            QueryRef::Mol(query),
            None,
            start_idx,
            end_idx,
            params,
            num_threads,
            max_results,
        )
    }

    /// `getMatches` for a tautomer query. The pattern fingerprints are only
    /// used when they are tautomeric ones.
    pub fn get_tautomer_matches(
        &self,
        query: &TautomerQuery,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
        max_results: i32,
    ) -> Result<Vec<u32>> {
        self.search(
            QueryRef::Tautomer(query),
            None,
            start_idx,
            end_idx,
            params,
            num_threads,
            max_results,
        )
    }

    /// `countMatches` for a tautomer query.
    pub fn count_tautomer_matches(
        &self,
        query: &TautomerQuery,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<u32> {
        Ok(self
            .get_tautomer_matches(query, start_idx, end_idx, params, num_threads, -1)?
            .len() as u32)
    }

    /// `hasMatch` for a tautomer query.
    pub fn has_tautomer_match(
        &self,
        query: &TautomerQuery,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<bool> {
        Ok(!self
            .get_tautomer_matches(query, start_idx, end_idx, params, num_threads, 1)?
            .is_empty())
    }

    /// The pattern fingerprint of the query to screen the library with
    /// (`Bits`), and the fingerprints to screen against; `None` when there
    /// is nothing to filter with.
    fn query_bits<'a>(
        fps: Option<&'a PatternHolder>,
        query: QueryRef,
    ) -> Result<(Option<&'a PatternHolder>, Option<ExplicitBitVect>)> {
        let Some(f) = fps else {
            return Ok((None, None));
        };
        let not_tautomeric = || {
            rd_log_warning!(
                "Pattern fingerprints for tautomersearch aren't tautomer fingerprints, ignoring..."
            );
            Ok((None, None))
        };
        // the bits every member of a bundle has
        let common = |fp_of: &mut dyn FnMut(usize) -> Result<ExplicitBitVect>, n: usize| {
            let mut bits: Option<ExplicitBitVect> = None;
            for i in 0..n {
                let fp = fp_of(i)?;
                bits = Some(match bits {
                    None => fp,
                    Some(b) => b.and(&fp),
                });
            }
            Ok::<_, RDKitError>(bits)
        };
        match query {
            QueryRef::Mol(q) => Ok((Some(f), Some(f.make_fingerprint(q)?))),
            QueryRef::Tautomer(tq) if f.is_tautomeric() => Ok((
                Some(f),
                Some(tq.pattern_fingerprint_template(f.num_bits())?),
            )),
            QueryRef::Tautomer(_) => not_tautomeric(),
            QueryRef::Extended(ExtendedQueryMol::Mol(q)) => {
                Ok((Some(f), Some(f.make_fingerprint(q)?)))
            }
            QueryRef::Extended(ExtendedQueryMol::Bundle(b)) => Ok((
                Some(f),
                common(&mut |i| f.make_fingerprint(&b.mols()[i]), b.size())?,
            )),
            QueryRef::Extended(ExtendedQueryMol::Tautomer(tq)) if f.is_tautomeric() => Ok((
                Some(f),
                Some(tq.pattern_fingerprint_template(f.num_bits())?),
            )),
            QueryRef::Extended(ExtendedQueryMol::TautomerBundle(v)) if f.is_tautomeric() => Ok((
                Some(f),
                common(
                    &mut |i| v[i].pattern_fingerprint_template(f.num_bits()),
                    v.len(),
                )?,
            )),
            QueryRef::Extended(_) => not_tautomeric(),
        }
    }

    fn search(
        &self,
        query: QueryRef,
        skip: Option<&[bool]>,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
        max_results: i32,
    ) -> Result<Vec<u32>> {
        let mols = read(&self.mols);
        let fps = self.fps.as_ref().map(read);
        let (fps_ref, query_bits) = Self::query_bits(fps.as_deref(), query)?;
        let search = search::Search {
            mols: &mols,
            fps: fps_ref,
            query,
            skip,
            query_bits,
            // only whether there is a match matters
            params: SubstructMatchParameters {
                max_matches: 1,
                ..params.clone()
            },
            needs_rings: query.needs_rings(),
            search_order: &self.search_order,
        };
        search::run(&search, start_idx, end_idx, num_threads, max_results)
    }

    /// `getMatches` for an extended query.
    pub fn get_xq_matches(
        &self,
        query: &ExtendedQueryMol,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
        max_results: i32,
    ) -> Result<Vec<u32>> {
        self.search(
            QueryRef::Extended(query),
            None,
            start_idx,
            end_idx,
            params,
            num_threads,
            max_results,
        )
    }

    /// `countMatches` for an extended query.
    pub fn count_xq_matches(
        &self,
        query: &ExtendedQueryMol,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<u32> {
        Ok(self
            .get_xq_matches(query, start_idx, end_idx, params, num_threads, -1)?
            .len() as u32)
    }

    /// `hasMatch` for an extended query.
    pub fn has_xq_match(
        &self,
        query: &ExtendedQueryMol,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<bool> {
        Ok(!self
            .get_xq_matches(query, start_idx, end_idx, params, num_threads, 1)?
            .is_empty())
    }

    /// `getMatches` for a bundle: the molecules matching any of its members,
    /// each reported once, member by member (`molbundleGetMatches`).
    pub fn get_bundle_matches(
        &self,
        query: &MolBundle,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
        max_results: i32,
    ) -> Result<Vec<u32>> {
        let mut found = vec![false; self.len() as usize];
        let mut res: Vec<u32> = Vec::new();
        for member in query.mols() {
            let remaining = if max_results > 0 {
                max_results - res.len() as i32
            } else {
                max_results
            };
            let hits = self.search(
                QueryRef::Mol(member),
                Some(&found),
                start_idx,
                end_idx,
                params,
                num_threads,
                remaining,
            )?;
            hits.iter().for_each(|&i| found[i as usize] = true);
            res.extend(hits);
            if max_results > 0 && res.len() as i32 >= max_results {
                break;
            }
        }
        Ok(res)
    }

    /// `countMatches` for a bundle.
    pub fn count_bundle_matches(
        &self,
        query: &MolBundle,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<u32> {
        Ok(self
            .get_bundle_matches(query, start_idx, end_idx, params, num_threads, -1)?
            .len() as u32)
    }

    /// `hasMatch` for a bundle.
    pub fn has_bundle_match(
        &self,
        query: &MolBundle,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<bool> {
        Ok(!self
            .get_bundle_matches(query, start_idx, end_idx, params, num_threads, 1)?
            .is_empty())
    }

    /// `countMatches`
    pub fn count_matches(
        &self,
        query: &ROMol,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<u32> {
        Ok(self
            .get_matches(query, start_idx, end_idx, params, num_threads, -1)?
            .len() as u32)
    }

    /// `hasMatch`
    pub fn has_match(
        &self,
        query: &ROMol,
        start_idx: u32,
        end_idx: u32,
        params: &SubstructMatchParameters,
        num_threads: i32,
    ) -> Result<bool> {
        Ok(!self
            .get_matches(query, start_idx, end_idx, params, num_threads, 1)?
            .is_empty())
    }

    /// `addPatterns`: fingerprints all molecules of a library that has none
    /// yet, in parallel, with `patterns` (default: a [`PatternHolder`]). The
    /// holder is shared with the caller; fingerprints it held are replaced.
    pub fn add_patterns(
        &mut self,
        patterns: Option<Shared<PatternHolder>>,
        num_threads: i32,
    ) -> Result<()> {
        if self.fps.is_some() {
            return Err(RDKitError::precondition(
                "Substruct library already has fingerprints",
            ));
        }
        let patterns = patterns.unwrap_or_else(|| shared(PatternHolder::default()));
        let fps = {
            let maker = read(&patterns);
            let mols = read(&self.mols);
            let ids: Vec<u32> = (0..mols.len()).collect();
            par_map(&ids, num_threads, |&i| match mols.get_mol(i)? {
                Some(mol) => maker.make_fingerprint(&mol),
                // an entry that cannot be built gets an empty fingerprint
                None => maker.make_fingerprint(&ROMol::new()),
            })
            .into_iter()
            .collect::<Result<Vec<_>>>()?
        };
        *patterns
            .write()
            .unwrap_or_else(PoisonError::into_inner)
            .fingerprints_mut() = fps;
        self.fps = Some(patterns);
        Ok(())
    }
}
