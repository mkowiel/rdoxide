//! Filter matchers (upstream `FilterCatalog/FilterMatcherBase.h`,
//! `FilterMatchers.{h,cpp}`): objects deciding whether a molecule matches a
//! filter, and combinators of them.

use std::any::Any;
use std::fmt::Debug;
use std::sync::{Arc, RwLock};

use crate::graphmol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::{RDKitError, Result};

pub const DEFAULT_FILTERMATCHERBASE_NAME: &str = "Unnamed FilterMatcherBase";
pub const SMARTS_MATCH_NAME_DEFAULT: &str = "Unnamed SmartsMatcher";

/// A matcher shared between catalogs, entries and combinators.
pub type SharedMatcher = Arc<dyn FilterMatcher>;

/// `FilterMatcherBase`
pub trait FilterMatcher: Send + Sync + Debug {
    fn is_valid(&self) -> bool;

    fn name(&self) -> String;

    /// Appends the matches to `matches`; returns whether the molecule matches.
    fn get_matches(&self, mol: &ROMol, matches: &mut Vec<FilterMatch>) -> Result<bool>;

    fn has_match(&self, mol: &ROMol) -> Result<bool>;

    /// A copy of the matcher (sub-matchers are shared).
    fn copy(&self) -> SharedMatcher;

    fn as_any(&self) -> &dyn Any;
}

/// `FilterMatch`: the matcher that triggered and the (query atom, target
/// atom) pairs of the match.
#[derive(Debug, Clone)]
pub struct FilterMatch {
    pub filter_match: SharedMatcher,
    pub atom_pairs: MatchVect,
}

impl PartialEq for FilterMatch {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::addr_eq(
            Arc::as_ptr(&self.filter_match),
            Arc::as_ptr(&other.filter_match),
        ) && self.atom_pairs == other.atom_pairs
    }
}

fn arg_name(arg: &Option<SharedMatcher>) -> String {
    arg.as_ref()
        .map_or_else(|| "<nullmatcher>".to_string(), |a| a.name())
}

fn arg_valid(arg: &Option<SharedMatcher>) -> bool {
    arg.as_ref().is_some_and(|a| a.is_valid())
}

fn required<'a>(arg: &'a Option<SharedMatcher>, what: &str) -> Result<&'a SharedMatcher> {
    arg.as_ref()
        .ok_or_else(|| RDKitError::precondition(what.to_string()))
}

/// `FilterMatchOps::And`
#[derive(Debug, Clone, Default)]
pub struct And {
    pub arg1: Option<SharedMatcher>,
    pub arg2: Option<SharedMatcher>,
}

impl And {
    pub fn new(arg1: SharedMatcher, arg2: SharedMatcher) -> Self {
        And {
            arg1: Some(arg1),
            arg2: Some(arg2),
        }
    }
}

impl FilterMatcher for And {
    fn is_valid(&self) -> bool {
        arg_valid(&self.arg1) && arg_valid(&self.arg2)
    }

    fn name(&self) -> String {
        format!("({} And {})", arg_name(&self.arg1), arg_name(&self.arg2))
    }

    fn has_match(&self, mol: &ROMol) -> Result<bool> {
        let msg = "FilterMatchOps::And is not valid, null arg1 or arg2";
        if !self.is_valid() {
            return Err(RDKitError::precondition(msg));
        }
        Ok(required(&self.arg1, msg)?.has_match(mol)?
            && required(&self.arg2, msg)?.has_match(mol)?)
    }

    fn get_matches(&self, mol: &ROMol, match_vect: &mut Vec<FilterMatch>) -> Result<bool> {
        let msg = "FilterMatchOps::And is not valid, null arg1 or arg2";
        if !self.is_valid() {
            return Err(RDKitError::precondition(msg));
        }
        let mut matches = Vec::new();
        if required(&self.arg1, msg)?.get_matches(mol, &mut matches)?
            && required(&self.arg2, msg)?.get_matches(mol, &mut matches)?
        {
            *match_vect = matches;
            return Ok(true);
        }
        Ok(false)
    }

    fn copy(&self) -> SharedMatcher {
        Arc::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// `FilterMatchOps::Or`
#[derive(Debug, Clone, Default)]
pub struct Or {
    pub arg1: Option<SharedMatcher>,
    pub arg2: Option<SharedMatcher>,
}

impl Or {
    pub fn new(arg1: SharedMatcher, arg2: SharedMatcher) -> Self {
        Or {
            arg1: Some(arg1),
            arg2: Some(arg2),
        }
    }
}

impl FilterMatcher for Or {
    fn is_valid(&self) -> bool {
        arg_valid(&self.arg1) && arg_valid(&self.arg2)
    }

    fn name(&self) -> String {
        format!("({} Or {})", arg_name(&self.arg1), arg_name(&self.arg2))
    }

    fn has_match(&self, mol: &ROMol) -> Result<bool> {
        let msg = "Or is not valid, null arg1 or arg2";
        if !self.is_valid() {
            return Err(RDKitError::precondition(msg));
        }
        Ok(required(&self.arg1, msg)?.has_match(mol)?
            || required(&self.arg2, msg)?.has_match(mol)?)
    }

    fn get_matches(&self, mol: &ROMol, match_vect: &mut Vec<FilterMatch>) -> Result<bool> {
        let msg = "FilterMatchOps::Or is not valid, null arg1 or arg2";
        if !self.is_valid() {
            return Err(RDKitError::precondition(msg));
        }
        let res1 = required(&self.arg1, msg)?.get_matches(mol, match_vect)?;
        let res2 = required(&self.arg2, msg)?.get_matches(mol, match_vect)?;
        Ok(res1 || res2)
    }

    fn copy(&self) -> SharedMatcher {
        Arc::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// `FilterMatchOps::Not`
#[derive(Debug, Clone, Default)]
pub struct Not {
    pub arg1: Option<SharedMatcher>,
}

impl Not {
    pub fn new(arg1: SharedMatcher) -> Self {
        Not { arg1: Some(arg1) }
    }
}

impl FilterMatcher for Not {
    fn is_valid(&self) -> bool {
        arg_valid(&self.arg1)
    }

    fn name(&self) -> String {
        format!("(Not {})", arg_name(&self.arg1))
    }

    fn has_match(&self, mol: &ROMol) -> Result<bool> {
        let arg = required(&self.arg1, "FilterMatchOps::Not: arg1 is null")?;
        if !arg.is_valid() {
            return Err(RDKitError::precondition(
                "FilterMatchOps::Not: arg1 is null",
            ));
        }
        Ok(!arg.has_match(mol)?)
    }

    fn get_matches(&self, mol: &ROMol, _: &mut Vec<FilterMatch>) -> Result<bool> {
        let arg = required(&self.arg1, "FilterMatchOps::Not: arg1 is null")?;
        if !arg.is_valid() {
            return Err(RDKitError::precondition(
                "FilterMatchOps::Not: arg1 is null",
            ));
        }
        let mut match_vect = Vec::new();
        Ok(!arg.get_matches(mol, &mut match_vect)?)
    }

    fn copy(&self) -> SharedMatcher {
        Arc::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// `SmartsMatcher`: a pattern that must occur between `min_count` and
/// `max_count` times.
#[derive(Debug, Clone)]
pub struct SmartsMatcher {
    name: String,
    pattern: Option<Arc<ROMol>>,
    min_count: u32,
    max_count: u32,
}

fn parse_smarts(smarts: &str) -> Option<Arc<ROMol>> {
    let params = SmartsParserParams {
        merge_hs: true,
        ..SmartsParserParams::default()
    };
    mol_from_smarts(smarts, &params)
        .ok()
        .flatten()
        .map(Arc::new)
}

impl SmartsMatcher {
    /// A matcher without a pattern (invalid until one is set).
    pub fn with_name(name: &str) -> Self {
        SmartsMatcher {
            name: name.to_string(),
            pattern: None,
            min_count: 0,
            max_count: u32::MAX,
        }
    }

    pub fn from_pattern(name: &str, pattern: Arc<ROMol>, min_count: u32, max_count: u32) -> Self {
        SmartsMatcher {
            name: name.to_string(),
            pattern: Some(pattern),
            min_count,
            max_count,
        }
    }

    /// The pattern comes from SMARTS; an unparsable SMARTS leaves the matcher
    /// invalid.
    pub fn from_smarts(name: &str, smarts: &str, min_count: u32, max_count: u32) -> Self {
        SmartsMatcher {
            name: name.to_string(),
            pattern: parse_smarts(smarts),
            min_count,
            max_count,
        }
    }

    pub fn pattern(&self) -> Option<&Arc<ROMol>> {
        self.pattern.as_ref()
    }

    pub fn set_pattern(&mut self, pattern: Arc<ROMol>) {
        self.pattern = Some(pattern);
    }

    pub fn set_pattern_smarts(&mut self, smarts: &str) {
        self.pattern = parse_smarts(smarts);
    }

    pub fn min_count(&self) -> u32 {
        self.min_count
    }

    pub fn set_min_count(&mut self, v: u32) {
        self.min_count = v;
    }

    pub fn max_count(&self) -> u32 {
        self.max_count
    }

    pub fn set_max_count(&mut self, v: u32) {
        self.max_count = v;
    }

    fn pattern_or_err(&self) -> Result<&Arc<ROMol>> {
        self.pattern
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("bad on pattern"))
    }

    fn count_ok(&self, count: u32) -> bool {
        count >= self.min_count && (self.max_count == u32::MAX || count <= self.max_count)
    }

    fn single_match_only(&self) -> bool {
        self.min_count == 1 && self.max_count == u32::MAX
    }
}

impl FilterMatcher for SmartsMatcher {
    fn is_valid(&self) -> bool {
        self.pattern.is_some()
    }

    fn name(&self) -> String {
        self.name.clone()
    }

    fn get_matches(&self, mol: &ROMol, match_vect: &mut Vec<FilterMatch>) -> Result<bool> {
        let pattern = self.pattern_or_err()?;
        if self.single_match_only() {
            let params = SubstructMatchParameters {
                max_matches: 1,
                ..Default::default()
            };
            let matches = substruct_match(mol, pattern, &params)?;
            if let Some(first) = matches.into_iter().next() {
                match_vect.push(FilterMatch {
                    filter_match: self.copy(),
                    atom_pairs: first,
                });
                return Ok(true);
            }
            return Ok(false);
        }
        let matches = substruct_match(mol, pattern, &SubstructMatchParameters::default())?;
        let on_pat_exists = self.count_ok(matches.len() as u32);
        if on_pat_exists {
            let clone = self.copy();
            for m in matches {
                match_vect.push(FilterMatch {
                    filter_match: Arc::clone(&clone),
                    atom_pairs: m,
                });
            }
        }
        Ok(on_pat_exists)
    }

    fn has_match(&self, mol: &ROMol) -> Result<bool> {
        let pattern = self.pattern_or_err()?;
        if self.single_match_only() {
            let params = SubstructMatchParameters {
                max_matches: 1,
                ..Default::default()
            };
            return Ok(!substruct_match(mol, pattern, &params)?.is_empty());
        }
        let matches = substruct_match(mol, pattern, &SubstructMatchParameters::default())?;
        Ok(self.count_ok(matches.len() as u32))
    }

    fn copy(&self) -> SharedMatcher {
        Arc::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// `ExclusionList`: matches molecules that match none of the patterns.
#[derive(Debug, Clone, Default)]
pub struct ExclusionList {
    off_patterns: Vec<SharedMatcher>,
}

impl ExclusionList {
    pub fn new() -> Self {
        ExclusionList::default()
    }

    pub fn from_patterns(patterns: Vec<SharedMatcher>) -> Self {
        ExclusionList {
            off_patterns: patterns,
        }
    }

    pub fn add_pattern(&mut self, base: &dyn FilterMatcher) -> Result<()> {
        if !base.is_valid() {
            return Err(RDKitError::precondition("Invalid FilterMatcherBase"));
        }
        self.off_patterns.push(base.copy());
        Ok(())
    }

    pub fn exclusion_patterns(&self) -> &[SharedMatcher] {
        &self.off_patterns
    }

    pub fn set_exclusion_patterns(&mut self, patterns: Vec<SharedMatcher>) {
        self.off_patterns = patterns;
    }

    fn check(&self, mol: &ROMol) -> Result<bool> {
        if !self.is_valid() {
            return Err(RDKitError::precondition(
                "ExclusionList: one of the exclusion pattens is invalid",
            ));
        }
        let mut result = true;
        for p in &self.off_patterns {
            if !result {
                break;
            }
            result &= !p.has_match(mol)?;
        }
        Ok(result)
    }
}

impl FilterMatcher for ExclusionList {
    fn is_valid(&self) -> bool {
        self.off_patterns.iter().all(|p| p.is_valid())
    }

    fn name(&self) -> String {
        let mut res = "(Not any of".to_string();
        for p in &self.off_patterns {
            res.push(' ');
            res.push_str(&p.name());
        }
        res.push(')');
        res
    }

    fn get_matches(&self, mol: &ROMol, _: &mut Vec<FilterMatch>) -> Result<bool> {
        self.check(mol)
    }

    fn has_match(&self, mol: &ROMol) -> Result<bool> {
        self.check(mol)
    }

    fn copy(&self) -> SharedMatcher {
        Arc::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// A node of a [`FilterHierarchyMatcher`] tree, shared between copies of the
/// parent like upstream's `shared_ptr` children.
pub type HierarchyNode = Arc<RwLock<FilterHierarchyMatcher>>;

/// `FilterHierarchyMatcher`: a matcher with child matchers; when a molecule
/// matches, the matches of the most specific children are returned.
#[derive(Debug, Clone, Default)]
pub struct FilterHierarchyMatcher {
    children: Vec<HierarchyNode>,
    matcher: Option<SharedMatcher>,
}

impl FilterHierarchyMatcher {
    /// A root node without a matcher.
    pub fn new() -> Self {
        FilterHierarchyMatcher::default()
    }

    pub fn with_matcher(matcher: &dyn FilterMatcher) -> Self {
        FilterHierarchyMatcher {
            children: Vec::new(),
            matcher: Some(matcher.copy()),
        }
    }

    pub fn matcher(&self) -> Option<&SharedMatcher> {
        self.matcher.as_ref()
    }

    pub fn children(&self) -> &[HierarchyNode] {
        &self.children
    }

    /// Rebuilds a node from its parts (deserialization).
    pub fn from_parts(matcher: Option<SharedMatcher>, children: Vec<HierarchyNode>) -> Self {
        FilterHierarchyMatcher { children, matcher }
    }

    pub fn set_pattern(&mut self, matcher: &dyn FilterMatcher) -> Result<()> {
        if !matcher.is_valid() {
            return Err(RDKitError::precondition(
                "Adding invalid patterns is not allowed.",
            ));
        }
        self.matcher = Some(matcher.copy());
        Ok(())
    }

    /// Adds a copy of `hierarchy` as a child and returns the stored child.
    pub fn add_child(&mut self, hierarchy: &FilterHierarchyMatcher) -> Result<HierarchyNode> {
        if !hierarchy.matcher.as_ref().is_some_and(|m| m.is_valid()) {
            return Err(RDKitError::precondition(
                "Only one root node is allowed in a FilterHierarchyMatcher",
            ));
        }
        let child = Arc::new(RwLock::new(hierarchy.clone()));
        self.children.push(Arc::clone(&child));
        Ok(child)
    }

    fn root_matcher(&self) -> Result<&SharedMatcher> {
        self.matcher
            .as_ref()
            .ok_or_else(|| RDKitError::precondition("FilterHierarchyMatcher has no matcher"))
    }
}

impl FilterMatcher for FilterHierarchyMatcher {
    fn is_valid(&self) -> bool {
        self.matcher.as_ref().is_some_and(|m| m.is_valid())
    }

    fn name(&self) -> String {
        match &self.matcher {
            Some(m) => m.name(),
            None => "FilterMatcherHierarchy root".to_string(),
        }
    }

    fn get_matches(&self, mol: &ROMol, m: &mut Vec<FilterMatch>) -> Result<bool> {
        let mut temp = Vec::new();
        let result = self.root_matcher()?.get_matches(mol, &mut temp)?;
        if result {
            let mut children = Vec::new();
            for child in &self.children {
                child
                    .read()
                    .unwrap_or_else(|e| e.into_inner())
                    .get_matches(mol, &mut children)?;
            }
            if children.is_empty() {
                m.extend(temp);
            } else {
                m.extend(children);
            }
        }
        Ok(result)
    }

    fn has_match(&self, mol: &ROMol) -> Result<bool> {
        let mut temp = Vec::new();
        self.get_matches(mol, &mut temp)
    }

    fn copy(&self) -> SharedMatcher {
        Arc::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}
