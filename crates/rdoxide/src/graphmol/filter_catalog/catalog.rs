//! `FilterCatalog` and its parameters (upstream `FilterCatalog/FilterCatalog.{h,cpp}`,
//! `Filters.{h,cpp}`, `FilterCatalogRunner.cpp`, `FunctionalGroupHierarchy.cpp`).
//! Serialization (Boost) is not supported, as in an upstream build without it.

use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock};

use super::entry::FilterCatalogEntry;
use super::matchers::{
    FilterHierarchyMatcher, FilterMatch, FilterMatcher, HierarchyNode, SmartsMatcher,
};
use crate::graphmol::ROMol;
use crate::graphmol::functional_group_hierarchy::functional_group_definitions;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

/// `FilterCatalogParams::FilterCatalogs` bit flags.
pub mod catalogs {
    pub const PAINS_A: u32 = 1 << 1;
    pub const PAINS_B: u32 = 1 << 2;
    pub const PAINS_C: u32 = 1 << 3;
    pub const PAINS: u32 = PAINS_A | PAINS_B | PAINS_C;
    pub const BRENK: u32 = 1 << 4;
    pub const NIH: u32 = 1 << 5;
    pub const ZINC: u32 = 1 << 6;
    pub const CHEMBL_GLAXO: u32 = 1 << 7;
    pub const CHEMBL_DUNDEE: u32 = 1 << 8;
    pub const CHEMBL_BMS: u32 = 1 << 9;
    pub const CHEMBL_SURECHEMBL: u32 = 1 << 10;
    pub const CHEMBL_MLSMR: u32 = 1 << 11;
    pub const CHEMBL_INPHARMATICA: u32 = 1 << 12;
    pub const CHEMBL_LINT: u32 = 1 << 13;
    pub const CHEMBL: u32 = CHEMBL_GLAXO
        | CHEMBL_DUNDEE
        | CHEMBL_BMS
        | CHEMBL_SURECHEMBL
        | CHEMBL_MLSMR
        | CHEMBL_INPHARMATICA
        | CHEMBL_LINT;
    pub const ALL: u32 = PAINS | BRENK | NIH | ZINC | CHEMBL;
}

/// One built-in filter (`FilterData_t`).
#[derive(Debug)]
struct FilterData {
    name: String,
    smarts: String,
    max: u32,
}

/// A built-in catalog: its properties and filters.
#[derive(Debug)]
struct BuiltinCatalog {
    props: Vec<(String, String)>,
    filters: Vec<FilterData>,
}

fn builtin_catalogs() -> &'static BTreeMap<&'static str, BuiltinCatalog> {
    static DATA: OnceLock<BTreeMap<&'static str, BuiltinCatalog>> = OnceLock::new();
    DATA.get_or_init(|| {
        let text = include_str!("../../../data/filtercatalog.txt");
        let mut map: BTreeMap<&'static str, BuiltinCatalog> = BTreeMap::new();
        let mut current: Option<&'static str> = None;
        for line in text.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            if let Some(name) = line.strip_prefix('@') {
                current = Some(name);
                map.insert(
                    name,
                    BuiltinCatalog {
                        props: Vec::new(),
                        filters: Vec::new(),
                    },
                );
                continue;
            }
            let cat = map.get_mut(current.expect("catalog header")).expect("catalog");
            let f: Vec<&str> = line.split('\t').collect();
            match f[0] {
                "P" => cat.props.push((f[1].to_string(), f[2].to_string())),
                "E" => cat.filters.push(FilterData {
                    name: f[1].to_string(),
                    smarts: f[2].to_string(),
                    max: f[3].parse().expect("max count"),
                }),
                other => unreachable!("bad filter catalog line type {other}"),
            }
        }
        map
    })
}

fn catalog_key(flag: u32) -> Option<&'static str> {
    Some(match flag {
        catalogs::BRENK => "BRENK",
        catalogs::NIH => "NIH",
        catalogs::PAINS_A => "PAINS_A",
        catalogs::PAINS_B => "PAINS_B",
        catalogs::PAINS_C => "PAINS_C",
        catalogs::ZINC => "ZINC",
        catalogs::CHEMBL_GLAXO => "CHEMBL_Glaxo",
        catalogs::CHEMBL_DUNDEE => "CHEMBL_Dundee",
        catalogs::CHEMBL_BMS => "CHEMBL_BMS",
        catalogs::CHEMBL_SURECHEMBL => "CHEMBL_SureChEMBL",
        catalogs::CHEMBL_MLSMR => "CHEMBL_MLSMR",
        catalogs::CHEMBL_INPHARMATICA => "CHEMBL_Inpharmatica",
        catalogs::CHEMBL_LINT => "CHEMBL_LINT",
        _ => return None,
    })
}

/// `FilterCatalogParams`
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FilterCatalogParams {
    catalogs: Vec<u32>,
}

impl FilterCatalogParams {
    pub fn new() -> Self {
        FilterCatalogParams::default()
    }

    pub fn with_catalogs(catalogs: u32) -> Self {
        let mut p = FilterCatalogParams::new();
        p.add_catalog(catalogs);
        p
    }

    /// Adds every known catalog among the bits of `catalog`; true if any was
    /// added.
    pub fn add_catalog(&mut self, catalog: u32) -> bool {
        let mut added = false;
        for i in 0..u32::BITS {
            let bit = catalog & (1 << i) & catalogs::ALL;
            if bit != 0 && catalog_key(bit).is_some() {
                self.catalogs.push(bit);
                added = true;
            }
        }
        added
    }

    pub fn catalogs(&self) -> &[u32] {
        &self.catalogs
    }

    /// `fillCatalog`
    pub fn fill_catalog(&self, catalog: &mut FilterCatalog) -> Result<()> {
        let data = builtin_catalogs();
        for &flag in &self.catalogs {
            let key = catalog_key(flag).expect("known catalog");
            let builtin = &data[key];
            for f in &builtin.filters {
                let entry = make_filter_catalog_entry(f, &builtin.props).ok_or_else(|| {
                    let name = builtin
                        .props
                        .iter()
                        .find(|(k, _)| k == "FilterSet")
                        .map_or("Unnamed internal catalog", |(_, v)| v.as_str());
                    RDKitError::Value(format!("Bad entry in built-in filter catalog: {name}"))
                })?;
                catalog.add_entry(Arc::new(entry));
            }
        }
        Ok(())
    }
}

/// `MakeFilterCatalogEntry`; `None` when the SMARTS does not parse.
fn make_filter_catalog_entry(
    data: &FilterData,
    props: &[(String, String)],
) -> Option<FilterCatalogEntry> {
    let matcher = SmartsMatcher::from_smarts(
        &data.name,
        &data.smarts,
        if data.max > 0 { data.max + 1 } else { 1 },
        u32::MAX,
    );
    let entry = FilterCatalogEntry::with_matcher(&data.name, Some(Arc::new(matcher)));
    if !entry.is_valid() {
        return None;
    }
    for (k, v) in props {
        entry.set_prop(k, v);
    }
    Some(entry)
}

/// The number of built-in entries of a catalog flag (`GetNumEntries`).
pub fn num_builtin_entries(catalog: u32) -> usize {
    catalog_key(catalog).map_or(0, |k| builtin_catalogs()[k].filters.len())
}

/// `FilterCatalog`
#[derive(Debug, Clone, Default)]
pub struct FilterCatalog {
    params: Option<FilterCatalogParams>,
    entries: Vec<Arc<FilterCatalogEntry>>,
}

impl FilterCatalog {
    pub fn new() -> Self {
        FilterCatalog::default()
    }

    pub fn from_params(params: &FilterCatalogParams) -> Result<Self> {
        let mut c = FilterCatalog::new();
        c.set_catalog_params(params)?;
        Ok(c)
    }

    pub fn from_catalogs(catalogs: u32) -> Result<Self> {
        Self::from_params(&FilterCatalogParams::with_catalogs(catalogs))
    }

    /// `setCatalogParams`: replaces the entries by those of the parameters.
    pub fn set_catalog_params(&mut self, params: &FilterCatalogParams) -> Result<()> {
        self.entries.clear();
        self.params = Some(params.clone());
        params.fill_catalog(self)
    }

    pub fn catalog_params(&self) -> Option<&FilterCatalogParams> {
        self.params.as_ref()
    }

    pub fn add_entry(&mut self, entry: Arc<FilterCatalogEntry>) -> usize {
        self.entries.push(entry);
        self.entries.len() - 1
    }

    pub fn remove_entry_at(&mut self, idx: usize) -> bool {
        if idx < self.entries.len() {
            self.entries.remove(idx);
            return true;
        }
        false
    }

    pub fn remove_entry(&mut self, entry: &Arc<FilterCatalogEntry>) -> bool {
        match self.entries.iter().position(|e| Arc::ptr_eq(e, entry)) {
            Some(i) => {
                self.entries.remove(i);
                true
            }
            None => false,
        }
    }

    pub fn num_entries(&self) -> usize {
        self.entries.len()
    }

    pub fn entry_with_idx(&self, idx: usize) -> Option<&Arc<FilterCatalogEntry>> {
        self.entries.get(idx)
    }

    /// `getEntry`: an out of range index is an `IndexError`.
    pub fn entry(&self, idx: usize) -> Result<&Arc<FilterCatalogEntry>> {
        self.entries.get(idx).ok_or(RDKitError::Index(idx.to_string()))
    }

    /// `getIdxForEntry`: `u32::MAX` when the entry is not in the catalog.
    pub fn idx_for_entry(&self, entry: &Arc<FilterCatalogEntry>) -> u32 {
        self.entries
            .iter()
            .position(|e| Arc::ptr_eq(e, entry))
            .map_or(u32::MAX, |i| i as u32)
    }

    pub fn has_match(&self, mol: &ROMol) -> Result<bool> {
        Ok(self.first_match(mol)?.is_some())
    }

    pub fn first_match(&self, mol: &ROMol) -> Result<Option<Arc<FilterCatalogEntry>>> {
        for e in &self.entries {
            if e.has_filter_match(mol)? {
                return Ok(Some(Arc::clone(e)));
            }
        }
        Ok(None)
    }

    pub fn matches(&self, mol: &ROMol) -> Result<Vec<Arc<FilterCatalogEntry>>> {
        let mut result = Vec::new();
        for e in &self.entries {
            if e.has_filter_match(mol)? {
                result.push(Arc::clone(e));
            }
        }
        Ok(result)
    }

    pub fn filter_matches(&self, mol: &ROMol) -> Result<Vec<FilterMatch>> {
        let mut result = Vec::new();
        for e in &self.entries {
            e.filter_matches(mol, &mut result)?;
        }
        Ok(result)
    }
}

/// Whether catalogs can be serialized (upstream: Boost serialization built in).
pub fn filter_catalog_can_serialize() -> bool {
    false
}

/// `RunFilterCatalog`: the entries matching each SMILES; an unparsable SMILES
/// gives the single entry "no valid RDKit molecule". The result does not
/// depend on the thread count.
pub fn run_filter_catalog(
    catalog: &FilterCatalog,
    smiles: &[String],
    num_threads: i32,
) -> Vec<Result<Vec<Arc<FilterCatalogEntry>>>> {
    let bad_smiles = Arc::new(FilterCatalogEntry::with_matcher("no valid RDKit molecule", None));
    par_map(smiles, num_threads, |smi| {
        match mol_from_smiles(smi, &SmilesParserParams::default()) {
            Ok(Some(mol)) => catalog.matches(&mol),
            _ => Ok(vec![Arc::clone(&bad_smiles)]),
        }
    })
}

/// `GetFunctionalGroupHierarchy`: one entry per top-level functional group,
/// each a [`FilterHierarchyMatcher`] tree.
pub fn functional_group_hierarchy() -> &'static FilterCatalog {
    static CATALOG: OnceLock<FilterCatalog> = OnceLock::new();
    CATALOG.get_or_init(|| {
        let mut catalog = FilterCatalog::new();
        let mut toplevel: Vec<HierarchyNode> = Vec::new();
        let mut stack: Vec<HierarchyNode> = Vec::new();
        for (level, name, pattern) in functional_group_definitions() {
            let node = FilterHierarchyMatcher::with_matcher(&SmartsMatcher::from_pattern(
                name,
                Arc::new(pattern),
                1,
                u32::MAX,
            ));
            if level == 0 {
                let top = Arc::new(std::sync::RwLock::new(node));
                toplevel.push(Arc::clone(&top));
                stack.clear();
                stack.push(top);
            } else {
                let child = stack[level - 1]
                    .write()
                    .unwrap_or_else(|e| e.into_inner())
                    .add_child(&node)
                    .expect("valid functional group");
                stack.truncate(level);
                stack.push(child);
            }
        }
        for top in toplevel {
            let snapshot = top.read().unwrap_or_else(|e| e.into_inner()).clone();
            let name = snapshot.name();
            catalog.add_entry(Arc::new(FilterCatalogEntry::with_matcher(
                &name,
                Some(Arc::new(snapshot)),
            )));
        }
        catalog
    })
}
