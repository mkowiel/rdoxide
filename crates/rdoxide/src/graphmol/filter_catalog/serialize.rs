//! Binary form of filter catalogs, entries and matchers.
//!
//! RDKit uses Boost serialization archives; rdoxide uses its own format (a
//! tag, then the fields, with query molecules stored as RDKit pickles), so
//! the strings are not exchangeable with RDKit.

use std::sync::{Arc, RwLock};

use super::catalog::{FilterCatalog, FilterCatalogParams};
use super::entry::FilterCatalogEntry;
use super::matchers::{
    And, ExclusionList, FilterHierarchyMatcher, FilterMatcher, Not, Or, SharedMatcher,
    SmartsMatcher,
};
use crate::graphmol::pickle::{flags, mol_from_pickle_stream, pickle_mol};
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

const MAGIC_CATALOG: &[u8; 4] = b"RDXC";
const MAGIC_ENTRY: &[u8; 4] = b"RDXE";

const T_AND: u8 = 1;
const T_OR: u8 = 2;
const T_NOT: u8 = 3;
const T_SMARTS: u8 = 4;
const T_EXCLUSION: u8 = 5;
const T_HIERARCHY: u8 = 6;

fn bad() -> RDKitError {
    RDKitError::Value("invalid filter catalog data".into())
}

struct Writer(Vec<u8>);

impl Writer {
    fn u8(&mut self, v: u8) {
        self.0.push(v);
    }
    fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn bytes(&mut self, b: &[u8]) {
        self.u32(b.len() as u32);
        self.0.extend_from_slice(b);
    }
    fn str(&mut self, s: &str) {
        self.bytes(s.as_bytes());
    }
    fn mol(&mut self, mol: &ROMol) -> Result<()> {
        let pkl = pickle_mol(mol, flags::ALL_PROPS)?;
        self.bytes(&pkl);
        Ok(())
    }
    fn opt_matcher(&mut self, m: Option<&SharedMatcher>) -> Result<()> {
        match m {
            None => self.u8(0),
            Some(m) => {
                self.u8(1);
                self.matcher(m.as_ref())?;
            }
        }
        Ok(())
    }
    fn matcher(&mut self, m: &dyn FilterMatcher) -> Result<()> {
        let any = m.as_any();
        if let Some(a) = any.downcast_ref::<And>() {
            self.u8(T_AND);
            self.opt_matcher(a.arg1.as_ref())?;
            self.opt_matcher(a.arg2.as_ref())?;
        } else if let Some(a) = any.downcast_ref::<Or>() {
            self.u8(T_OR);
            self.opt_matcher(a.arg1.as_ref())?;
            self.opt_matcher(a.arg2.as_ref())?;
        } else if let Some(a) = any.downcast_ref::<Not>() {
            self.u8(T_NOT);
            self.opt_matcher(a.arg1.as_ref())?;
        } else if let Some(s) = any.downcast_ref::<SmartsMatcher>() {
            self.u8(T_SMARTS);
            self.str(&s.name());
            self.u32(s.min_count());
            self.u32(s.max_count());
            match s.pattern() {
                None => self.u8(0),
                Some(p) => {
                    self.u8(1);
                    self.mol(p)?;
                }
            }
        } else if let Some(e) = any.downcast_ref::<ExclusionList>() {
            self.u8(T_EXCLUSION);
            self.u32(e.exclusion_patterns().len() as u32);
            for p in e.exclusion_patterns() {
                self.matcher(p.as_ref())?;
            }
        } else if let Some(h) = any.downcast_ref::<FilterHierarchyMatcher>() {
            self.hierarchy(h)?;
        } else {
            return Err(RDKitError::Value(
                "cannot serialize this filter matcher".into(),
            ));
        }
        Ok(())
    }
    fn hierarchy(&mut self, h: &FilterHierarchyMatcher) -> Result<()> {
        self.u8(T_HIERARCHY);
        self.opt_matcher(h.matcher())?;
        self.u32(h.children().len() as u32);
        for c in h.children() {
            self.hierarchy(&c.read().unwrap_or_else(|e| e.into_inner()))?;
        }
        Ok(())
    }
    fn entry(&mut self, e: &FilterCatalogEntry) -> Result<()> {
        self.opt_matcher(e.matcher())?;
        let keys = e.prop_list();
        self.u32(keys.len() as u32);
        for k in keys {
            self.str(&k);
            self.str(&e.prop(&k)?);
        }
        Ok(())
    }
}

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let end = self.pos.checked_add(n).ok_or_else(bad)?;
        let b = self.data.get(self.pos..end).ok_or_else(bad)?;
        self.pos = end;
        Ok(b)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }
    fn bytes(&mut self) -> Result<&[u8]> {
        let n = self.u32()? as usize;
        self.take(n)
    }
    fn string(&mut self) -> Result<String> {
        String::from_utf8(self.bytes()?.to_vec()).map_err(|_| bad())
    }
    fn mol(&mut self) -> Result<ROMol> {
        let data = self.bytes()?;
        let mut p = 0;
        mol_from_pickle_stream(data, &mut p, flags::ALL_PROPS)
    }
    fn opt_matcher(&mut self) -> Result<Option<SharedMatcher>> {
        match self.u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.matcher()?)),
            _ => Err(bad()),
        }
    }
    fn matcher(&mut self) -> Result<SharedMatcher> {
        Ok(match self.u8()? {
            T_AND => Arc::new(And {
                arg1: self.opt_matcher()?,
                arg2: self.opt_matcher()?,
            }),
            T_OR => Arc::new(Or {
                arg1: self.opt_matcher()?,
                arg2: self.opt_matcher()?,
            }),
            T_NOT => Arc::new(Not {
                arg1: self.opt_matcher()?,
            }),
            T_SMARTS => {
                let name = self.string()?;
                let (min, max) = (self.u32()?, self.u32()?);
                let mut m = SmartsMatcher::with_name(&name);
                m.set_min_count(min);
                m.set_max_count(max);
                if self.u8()? == 1 {
                    m.set_pattern(Arc::new(self.mol()?));
                }
                Arc::new(m)
            }
            T_EXCLUSION => {
                let n = self.u32()?;
                let pats = (0..n).map(|_| self.matcher()).collect::<Result<Vec<_>>>()?;
                Arc::new(ExclusionList::from_patterns(pats))
            }
            T_HIERARCHY => Arc::new(self.hierarchy_body()?),
            _ => return Err(bad()),
        })
    }
    /// The body of a hierarchy node whose tag was already read.
    fn hierarchy_body(&mut self) -> Result<FilterHierarchyMatcher> {
        let matcher = self.opt_matcher()?;
        let n = self.u32()?;
        let mut children = Vec::new();
        for _ in 0..n {
            if self.u8()? != T_HIERARCHY {
                return Err(bad());
            }
            children.push(Arc::new(RwLock::new(self.hierarchy_body()?)));
        }
        Ok(FilterHierarchyMatcher::from_parts(matcher, children))
    }
    fn entry(&mut self) -> Result<FilterCatalogEntry> {
        let matcher = self.opt_matcher()?;
        let entry = FilterCatalogEntry::with_matcher("", matcher);
        entry.clear_prop("description");
        let n = self.u32()?;
        for _ in 0..n {
            let k = self.string()?;
            let v = self.string()?;
            entry.set_prop(&k, &v);
        }
        Ok(entry)
    }
}

/// `FilterCatalogEntry::Serialize`
pub fn entry_to_binary(entry: &FilterCatalogEntry) -> Result<Vec<u8>> {
    let mut w = Writer(MAGIC_ENTRY.to_vec());
    w.entry(entry)?;
    Ok(w.0)
}

/// `FilterCatalogEntry(pickle)`
pub fn entry_from_binary(data: &[u8]) -> Result<FilterCatalogEntry> {
    if data.len() < 4 || &data[..4] != MAGIC_ENTRY {
        return Err(bad());
    }
    Reader { data, pos: 4 }.entry()
}

/// `FilterCatalog::Serialize`
pub fn catalog_to_binary(catalog: &FilterCatalog) -> Result<Vec<u8>> {
    let mut w = Writer(MAGIC_CATALOG.to_vec());
    match catalog.catalog_params() {
        None => w.u8(0),
        Some(p) => {
            w.u8(1);
            w.u32(p.catalogs().len() as u32);
            p.catalogs().iter().for_each(|&c| w.u32(c));
        }
    }
    w.u32(catalog.num_entries() as u32);
    for e in catalog.entries() {
        w.entry(e)?;
    }
    Ok(w.0)
}

/// `FilterCatalog(pickle)`
pub fn catalog_from_binary(data: &[u8]) -> Result<FilterCatalog> {
    if data.len() < 4 || &data[..4] != MAGIC_CATALOG {
        return Err(bad());
    }
    let mut r = Reader { data, pos: 4 };
    let params = match r.u8()? {
        0 => None,
        1 => {
            let n = r.u32()?;
            let flags = (0..n).map(|_| r.u32()).collect::<Result<Vec<_>>>()?;
            Some(FilterCatalogParams::from_flags(flags))
        }
        _ => return Err(bad()),
    };
    let n = r.u32()?;
    let entries = (0..n)
        .map(|_| r.entry().map(Arc::new))
        .collect::<Result<Vec<_>>>()?;
    Ok(FilterCatalog::from_parts(params, entries))
}
