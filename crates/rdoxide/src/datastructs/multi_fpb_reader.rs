//! Searching several FPB files at once (upstream
//! `DataStructs/MultiFPBReader.{h,cpp}`). The readers are searched in
//! parallel when more than one thread is requested.

use std::sync::{Arc, RwLock};

use super::fpb_reader::FpbReader;
use crate::error::{RDKitError, Result};
use crate::rdgeneral::threads::par_map;

/// An [`FpbReader`] that can be shared, e.g. between Python objects and a
/// [`MultiFpbReader`].
pub type SharedFpbReader = Arc<RwLock<FpbReader>>;

/// A `(similarity, fingerprint index, reader index)` search hit.
pub type MultiNeighbor = (f64, u32, u32);

/// A set of [`FpbReader`]s searched together.
#[derive(Debug, Default)]
pub struct MultiFpbReader {
    readers: Vec<SharedFpbReader>,
    initialized: bool,
    init_on_search: bool,
}

fn read_lock(reader: &SharedFpbReader) -> std::sync::RwLockReadGuard<'_, FpbReader> {
    reader.read().unwrap_or_else(|e| e.into_inner())
}

fn write_lock(reader: &SharedFpbReader) -> std::sync::RwLockWriteGuard<'_, FpbReader> {
    reader.write().unwrap_or_else(|e| e.into_inner())
}

impl MultiFpbReader {
    /// With `init_on_search` the child readers are only initialized by the
    /// first search, which is useful with large files.
    pub fn new(init_on_search: bool) -> Self {
        MultiFpbReader {
            readers: Vec::new(),
            initialized: false,
            init_on_search,
        }
    }

    /// Initializes every child reader and checks that their fingerprints
    /// have the same length.
    pub fn init(&mut self) -> Result<()> {
        let mut n_bits = 0;
        for reader in &self.readers {
            write_lock(reader).init()?;
            let this = read_lock(reader).n_bits()?;
            if n_bits == 0 {
                n_bits = this;
            } else if this != n_bits {
                return Err(RDKitError::Value(
                    "bit lengths of child readers don't match".into(),
                ));
            }
        }
        self.initialized = true;
        Ok(())
    }

    /// Number of readers.
    pub fn len(&self) -> usize {
        self.readers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.readers.is_empty()
    }

    /// Number of bits in the fingerprints.
    pub fn n_bits(&self) -> Result<u32> {
        if self.readers.is_empty() {
            return Err(RDKitError::precondition("no readers"));
        }
        if !self.initialized {
            return Err(RDKitError::precondition("not initialized"));
        }
        read_lock(&self.readers[0]).n_bits()
    }

    pub fn reader(&self, which: usize) -> Result<&SharedFpbReader> {
        self.readers
            .get(which)
            .ok_or_else(|| RDKitError::range("which", which, self.readers.len()))
    }

    /// Adds a reader (initializing it if this reader already is) and returns
    /// the number of readers.
    pub fn add_reader(&mut self, reader: SharedFpbReader) -> Result<usize> {
        if self.initialized {
            write_lock(&reader).init()?;
        }
        self.readers.push(reader);
        Ok(self.readers.len())
    }

    /// Runs `search` on every reader, on up to `num_threads` threads (zero
    /// or negative: all hardware threads minus `|num_threads|`), and returns
    /// the per-reader results in reader order.
    fn search_all<T: Send>(
        &self,
        num_threads: i32,
        search: impl Fn(&FpbReader) -> Result<Vec<T>> + Sync,
    ) -> Result<Vec<Vec<T>>> {
        if !(self.initialized || self.init_on_search) {
            return Err(RDKitError::precondition("not initialized"));
        }
        let per_reader = par_map(&self.readers, num_threads, |reader| {
            if self.init_on_search {
                write_lock(reader).init()?;
            }
            search(&read_lock(reader))
        });
        per_reader.into_iter().collect()
    }

    fn merge_neighbors(per_reader: Vec<Vec<(f64, u32)>>) -> Vec<MultiNeighbor> {
        let mut res: Vec<MultiNeighbor> = per_reader
            .into_iter()
            .zip(0u32..)
            .flat_map(|(hits, reader)| hits.into_iter().map(move |(sim, idx)| (sim, idx, reader)))
            .collect();
        // most similar first, then by reader and index
        res.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.2.cmp(&b.2)).then(a.1.cmp(&b.1)));
        res
    }

    /// Tanimoto neighbors over all readers, most similar first.
    pub fn get_tanimoto_neighbors(
        &self,
        query: &[u8],
        threshold: f64,
        num_threads: i32,
    ) -> Result<Vec<MultiNeighbor>> {
        let hits = self.search_all(num_threads, |r| {
            r.get_tanimoto_neighbors(query, threshold, true)
        })?;
        Ok(Self::merge_neighbors(hits))
    }

    /// Tversky neighbors over all readers, most similar first.
    pub fn get_tversky_neighbors(
        &self,
        query: &[u8],
        ca: f64,
        cb: f64,
        threshold: f64,
        num_threads: i32,
    ) -> Result<Vec<MultiNeighbor>> {
        let hits = self.search_all(num_threads, |r| {
            r.get_tversky_neighbors(query, ca, cb, threshold, true)
        })?;
        Ok(Self::merge_neighbors(hits))
    }

    /// `(fingerprint index, reader index)` of every fingerprint containing
    /// all bits of `query`, ordered by index and then reader.
    pub fn get_containing_neighbors(
        &self,
        query: &[u8],
        num_threads: i32,
    ) -> Result<Vec<(u32, u32)>> {
        let hits = self.search_all(num_threads, |r| r.get_containing_neighbors(query))?;
        let mut res: Vec<(u32, u32)> = hits
            .into_iter()
            .zip(0u32..)
            .flat_map(|(idxs, reader)| idxs.into_iter().map(move |idx| (idx, reader)))
            .collect();
        res.sort_unstable();
        Ok(res)
    }
}
