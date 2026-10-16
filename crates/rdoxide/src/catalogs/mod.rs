//! Port of RDKit `Code/Catalogs`: the catalog framework (`CatalogEntry`,
//! `CatalogParams`, `HierarchCatalog`) and its binary pickle format.

use std::collections::BTreeMap;

use crate::{RDKitError, Result};

pub const VERSION_MAJOR: i32 = 1;
pub const VERSION_MINOR: i32 = 0;
pub const VERSION_PATCH: i32 = 0;
pub const ENDIAN_ID: i32 = 0xDEADBEEF_u32 as i32;

/// Writer of the little-endian stream format of `RDKit::streamWrite`.
#[derive(Debug, Default)]
pub struct StreamWriter(pub Vec<u8>);

impl StreamWriter {
    pub fn i32(&mut self, v: i32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    pub fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    pub fn f64(&mut self, v: f64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    pub fn bytes(&mut self, v: &[u8]) {
        self.0.extend_from_slice(v);
    }
    /// A string as its length (int) followed by the characters.
    pub fn string(&mut self, s: &str) {
        self.i32(s.len() as i32);
        self.bytes(s.as_bytes());
    }
}

/// Reader of the stream format.
#[derive(Debug)]
pub struct StreamReader<'a> {
    pub data: &'a [u8],
    pub pos: usize,
}

impl<'a> StreamReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        StreamReader { data, pos: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.data.len())
            .ok_or_else(|| RDKitError::Value("Bad pickle format: truncated data".into()))?;
        let s = &self.data[self.pos..end];
        self.pos = end;
        Ok(s)
    }
    pub fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }
    pub fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }
    pub fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(
            self.take(8)?.try_into().expect("8 bytes"),
        ))
    }
    pub fn string(&mut self) -> Result<String> {
        let n = self.i32()?;
        let n = usize::try_from(n)
            .map_err(|_| RDKitError::Value("Bad pickle format: negative length".into()))?;
        Ok(String::from_utf8_lossy(self.take(n)?).into_owned())
    }
}

/// `RDCatalog::CatalogParams`
pub trait CatalogParams: Sized + Clone + Default {
    fn type_str(&self) -> &str;
    fn to_stream(&self, w: &mut StreamWriter);
    fn init_from_stream(&mut self, r: &mut StreamReader) -> Result<()>;
}

/// `RDCatalog::CatalogEntry`, with the order the hierarchical catalog files
/// it under.
pub trait CatalogEntry: Sized {
    type Order: Ord + Copy;
    fn bit_id(&self) -> i32;
    fn set_bit_id(&mut self, id: i32);
    fn description(&self) -> String;
    fn order(&self) -> Self::Order;
    fn to_stream(&self, w: &mut StreamWriter) -> Result<()>;
    fn from_stream(r: &mut StreamReader) -> Result<Self>;
}

/// `RDCatalog::HierarchCatalog`: entries arranged in a directed graph.
///
/// *Indices* refer to the entries of the underlying structure, *bit ids* to
/// the "interesting" entries (those that were added with `update_fp_length`).
#[derive(Debug, Clone)]
pub struct HierarchCatalog<E: CatalogEntry, P: CatalogParams> {
    fp_length: u32,
    params: Option<P>,
    entries: Vec<E>,
    /// outgoing edges, in the order they were added
    down: Vec<Vec<u32>>,
    order_map: BTreeMap<E::Order, Vec<u32>>,
}

impl<E: CatalogEntry, P: CatalogParams> Default for HierarchCatalog<E, P> {
    fn default() -> Self {
        HierarchCatalog {
            fp_length: 0,
            params: None,
            entries: Vec::new(),
            down: Vec::new(),
            order_map: BTreeMap::new(),
        }
    }
}

fn range_error(idx: u32) -> RDKitError {
    RDKitError::Index(idx.to_string())
}

impl<E: CatalogEntry, P: CatalogParams> HierarchCatalog<E, P> {
    pub fn new(params: P) -> Self {
        HierarchCatalog {
            params: Some(params),
            ..Default::default()
        }
    }

    pub fn fp_length(&self) -> u32 {
        self.fp_length
    }

    pub fn set_fp_length(&mut self, val: u32) {
        self.fp_length = val;
    }

    pub fn catalog_params(&self) -> Option<&P> {
        self.params.as_ref()
    }

    /// `setCatalogParams`: only once.
    pub fn set_catalog_params(&mut self, params: P) -> Result<()> {
        if self.params.is_some() {
            return Err(RDKitError::precondition(
                "A parameter object already exists on the catalog",
            ));
        }
        self.params = Some(params);
        Ok(())
    }

    pub fn num_entries(&self) -> u32 {
        self.entries.len() as u32
    }

    /// `addEntry`: the index of the new entry. With `update_fp_length` the
    /// entry gets the next bit id.
    pub fn add_entry(&mut self, mut entry: E, update_fp_length: bool) -> u32 {
        if update_fp_length {
            entry.set_bit_id(self.fp_length as i32);
            self.fp_length += 1;
        }
        let eid = self.entries.len() as u32;
        self.order_map.entry(entry.order()).or_default().push(eid);
        self.entries.push(entry);
        self.down.push(Vec::new());
        eid
    }

    /// `addEdge`: from `id1` to `id2` (a no-op when it exists already).
    pub fn add_edge(&mut self, id1: u32, id2: u32) -> Result<()> {
        let n = self.num_entries();
        if id1 >= n {
            return Err(range_error(id1));
        }
        if id2 >= n {
            return Err(range_error(id2));
        }
        if !self.down[id1 as usize].contains(&id2) {
            self.down[id1 as usize].push(id2);
        }
        Ok(())
    }

    pub fn entry_with_idx(&self, idx: u32) -> Result<&E> {
        self.entries
            .get(idx as usize)
            .ok_or_else(|| range_error(idx))
    }

    pub fn entry_with_idx_mut(&mut self, idx: u32) -> Result<&mut E> {
        self.entries
            .get_mut(idx as usize)
            .ok_or_else(|| range_error(idx))
    }

    /// `getIdOfEntryWithBitId`: -1 when there is none.
    pub fn id_of_entry_with_bit_id(&self, bit_id: u32) -> Result<i32> {
        if bit_id >= self.fp_length {
            return Err(range_error(bit_id));
        }
        Ok((bit_id as usize..self.entries.len())
            .find(|&i| self.entries[i].bit_id() == bit_id as i32)
            .map_or(-1, |i| i as i32))
    }

    pub fn entry_with_bit_id(&self, bit_id: u32) -> Result<Option<&E>> {
        let id = self.id_of_entry_with_bit_id(bit_id)?;
        Ok((id >= 0).then(|| &self.entries[id as usize]))
    }

    /// `getDownEntryList`
    pub fn down_entry_list(&self, idx: u32) -> Result<&[u32]> {
        self.down
            .get(idx as usize)
            .map(Vec::as_slice)
            .ok_or_else(|| range_error(idx))
    }

    /// `getEntriesOfOrder`
    pub fn entries_of_order(&self, order: E::Order) -> Result<&[u32]> {
        self.order_map
            .get(&order)
            .map(Vec::as_slice)
            .ok_or_else(|| {
                RDKitError::invariant(
                    " catalog does not contain any entries of the order specified",
                )
            })
    }

    /// `toStream`/`Serialize`
    pub fn serialize(&self) -> Result<Vec<u8>> {
        let params = self
            .params
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("NULL parameter object"))?;
        let mut w = StreamWriter::default();
        w.i32(ENDIAN_ID);
        w.i32(VERSION_MAJOR);
        w.i32(VERSION_MINOR);
        w.i32(VERSION_PATCH);
        w.u32(self.fp_length);
        w.u32(self.num_entries());
        params.to_stream(&mut w);
        for e in &self.entries {
            e.to_stream(&mut w)?;
        }
        for children in &self.down {
            w.u32(children.len() as u32);
            children.iter().for_each(|&c| w.i32(c as i32));
        }
        Ok(w.0)
    }

    /// `initFromString`
    pub fn from_pickle(data: &[u8]) -> Result<Self> {
        let mut r = StreamReader::new(data);
        for _ in 0..4 {
            r.i32()?;
        }
        let mut cat = Self {
            fp_length: r.u32()?,
            ..Self::default()
        };
        let n = r.u32()?;
        let mut params = P::default();
        params.init_from_stream(&mut r)?;
        cat.params = Some(params);
        for _ in 0..n {
            let e = E::from_stream(&mut r)?;
            cat.add_entry(e, false);
        }
        for i in 0..n {
            let nbrs = r.u32()?;
            for _ in 0..nbrs {
                let c = r.i32()?;
                cat.add_edge(i, c as u32)?;
            }
        }
        Ok(cat)
    }
}
