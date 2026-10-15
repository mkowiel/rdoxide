//! `FilterCatalogEntry` (upstream `FilterCatalog/FilterCatalogEntry.{h,cpp}`):
//! a named matcher with string properties.

use std::sync::RwLock;

use super::matchers::{FilterMatch, SharedMatcher};
use crate::graphmol::ROMol;
use crate::rdgeneral::Dict;
use crate::{RDKitError, Result};

const DESCRIPTION: &str = "description";

#[derive(Debug, Default)]
pub struct FilterCatalogEntry {
    matcher: Option<SharedMatcher>,
    props: RwLock<Dict>,
}

impl Clone for FilterCatalogEntry {
    fn clone(&self) -> Self {
        FilterCatalogEntry {
            matcher: self.matcher.clone(),
            props: RwLock::new(self.props_read().clone()),
        }
    }
}

impl FilterCatalogEntry {
    /// An entry without matcher or description.
    pub fn new() -> Self {
        FilterCatalogEntry::default()
    }

    pub fn with_matcher(name: &str, matcher: Option<SharedMatcher>) -> Self {
        let entry = FilterCatalogEntry {
            matcher,
            props: RwLock::new(Dict::new()),
        };
        entry.set_description(name);
        entry
    }

    fn props_read(&self) -> std::sync::RwLockReadGuard<'_, Dict> {
        self.props.read().unwrap_or_else(|e| e.into_inner())
    }

    fn props_write(&self) -> std::sync::RwLockWriteGuard<'_, Dict> {
        self.props.write().unwrap_or_else(|e| e.into_inner())
    }

    pub fn is_valid(&self) -> bool {
        self.matcher.as_ref().is_some_and(|m| m.is_valid())
    }

    pub fn matcher(&self) -> Option<&SharedMatcher> {
        self.matcher.as_ref()
    }

    pub fn description(&self) -> String {
        self.props_read()
            .get(DESCRIPTION)
            .map(|v| v.to_rd_string())
            .unwrap_or_default()
    }

    pub fn set_description(&self, description: &str) {
        self.props_write().set_val(DESCRIPTION, description);
    }

    pub fn prop_list(&self) -> Vec<String> {
        self.props_read().keys().map(str::to_string).collect()
    }

    pub fn set_prop(&self, key: &str, val: &str) {
        self.props_write().set_val(key, val);
    }

    /// `getProp<std::string>`: a missing key is a `KeyError`.
    pub fn prop(&self, key: &str) -> Result<String> {
        self.props_read()
            .get(key)
            .map(|v| v.to_rd_string())
            .ok_or_else(|| RDKitError::Key(key.to_string()))
    }

    pub fn has_prop(&self, key: &str) -> bool {
        self.props_read().has_val(key)
    }

    pub fn clear_prop(&self, key: &str) {
        self.props_write().clear_val(key);
    }

    /// Appends the matches; false when the entry is invalid or does not match.
    pub fn filter_matches(&self, mol: &ROMol, matches: &mut Vec<FilterMatch>) -> Result<bool> {
        match &self.matcher {
            Some(m) if m.is_valid() => m.get_matches(mol, matches),
            _ => Ok(false),
        }
    }

    pub fn has_filter_match(&self, mol: &ROMol) -> Result<bool> {
        match &self.matcher {
            Some(m) if m.is_valid() => m.has_match(mol),
            _ => Ok(false),
        }
    }
}
