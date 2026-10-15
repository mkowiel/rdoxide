//! Substructure matching: VF2 (`vf2.hpp`), `SubstructMatch`,
//! `SubstructMatchCount`, atom/bond compatibility and recursive queries.
//!
//! Ported from RDKit `Code/GraphMol/Substruct/` (2026.09.1).
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//! The VF2 state machine is an RDKit modification of vflib-2.0 by P. Foggia.
//!
//! Upstream fills the sets of recursive (`$(...)`) queries in place while
//! matching and clears them afterwards; here those results live in a
//! [`RecursionSets`] side table keyed by query node for the duration of the
//! search, which is observably equivalent.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt;
use std::sync::Arc;

use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::queries::{
    AtomQuery, BondQuery, QueryKind, RecursionSets, atom_queries_match, has_recursive_query,
};
use crate::graphmol::romol::ROMol;
use crate::graphmol::stereo_group::StereoGroupType;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::ptree::PTree;
use crate::rdgeneral::utils::count_swaps_to_interconvert;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// A match: (query atom index, molecule atom index) pairs ordered by query
/// atom (`MatchVectType`).
pub type MatchVect = Vec<(i32, i32)>;

/// A callback comparing one query atom/bond with one molecule atom/bond:
/// `(query, query index, molecule, molecule index)` (`extraAtomCheck` /
/// `extraBondCheck`). It may be called from several threads.
#[derive(Clone)]
pub struct MatchCheck(Arc<MatchCheckFn>);

type MatchCheckFn = dyn Fn(&ROMol, u32, &ROMol, u32) -> bool + Send + Sync;

impl MatchCheck {
    pub fn new(f: impl Fn(&ROMol, u32, &ROMol, u32) -> bool + Send + Sync + 'static) -> Self {
        MatchCheck(Arc::new(f))
    }

    fn call(&self, query: &ROMol, qi: u32, mol: &ROMol, mi: u32) -> bool {
        (self.0)(query, qi, mol, mi)
    }
}

/// A callback validating a finished match: `(molecule, matched atom indices
/// in query atom order)` (`extraFinalCheck`).
#[derive(Clone)]
pub struct FinalMatchCheck(Arc<FinalMatchCheckFn>);

type FinalMatchCheckFn = dyn Fn(&ROMol, &[u32]) -> bool + Send + Sync;

impl FinalMatchCheck {
    pub fn new(f: impl Fn(&ROMol, &[u32]) -> bool + Send + Sync + 'static) -> Self {
        FinalMatchCheck(Arc::new(f))
    }
}

/// `AtomCoordsMatchFunctor`: atoms match when their positions in the
/// reference (molecule) and query conformers are within a tolerance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomCoordsMatcher {
    pub ref_conf_id: i32,
    pub query_conf_id: i32,
    /// squared distance tolerance
    pub tol2: f64,
}

impl AtomCoordsMatcher {
    pub fn new(ref_conf_id: i32, query_conf_id: i32, tol: f64) -> Self {
        AtomCoordsMatcher {
            ref_conf_id,
            query_conf_id,
            tol2: tol * tol,
        }
    }

    /// Molecules without conformers (or without the requested one) never match.
    pub fn matches(&self, query: &ROMol, qi: u32, mol: &ROMol, mi: u32) -> bool {
        if query.num_conformers() == 0 || mol.num_conformers() == 0 {
            return false;
        }
        let pos = |m: &ROMol, conf: i32, idx: u32| m.conformer(conf)?.atom_pos(idx);
        match (
            pos(query, self.query_conf_id, qi),
            pos(mol, self.ref_conf_id, mi),
        ) {
            (Ok(q), Ok(t)) => (q - t).length_sq() <= self.tol2,
            _ => false,
        }
    }
}

impl From<AtomCoordsMatcher> for MatchCheck {
    fn from(m: AtomCoordsMatcher) -> Self {
        MatchCheck::new(move |q, qi, mol, mi| m.matches(q, qi, mol, mi))
    }
}

macro_rules! callback_traits {
    ($($t:ty),*) => {$(
        impl fmt::Debug for $t {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(stringify!($t))
            }
        }
        /// Callbacks are equal when they are the same closure.
        impl PartialEq for $t {
            fn eq(&self, other: &Self) -> bool {
                Arc::ptr_eq(&self.0, &other.0)
            }
        }
        impl Eq for $t {}
    )*};
}
callback_traits!(MatchCheck, FinalMatchCheck);

/// `SubstructMatchParameters`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubstructMatchParameters {
    pub use_chirality: bool,
    pub use_enhanced_stereo: bool,
    pub aromatic_matches_conjugated: bool,
    pub use_query_query_matches: bool,
    pub use_generic_matchers: bool,
    pub recursion_possible: bool,
    pub uniquify: bool,
    pub max_matches: u32,
    pub num_threads: i32,
    pub atom_properties: Vec<String>,
    pub bond_properties: Vec<String>,
    pub max_recursive_matches: u32,
    pub specified_stereo_query_matches_unspecified: bool,
    pub aromatic_matches_single_or_double: bool,
    /// Called for atom pairs that passed every other comparison.
    pub extra_atom_check: Option<MatchCheck>,
    /// Use only `extra_atom_check` to decide whether atoms match.
    pub extra_atom_check_overrides_default: bool,
    /// Called for bond pairs that passed every other comparison.
    pub extra_bond_check: Option<MatchCheck>,
    /// Use only `extra_bond_check` to decide whether bonds match.
    pub extra_bond_check_overrides_default: bool,
    /// Called with each finished match.
    pub extra_final_check: Option<FinalMatchCheck>,
}

impl Default for SubstructMatchParameters {
    fn default() -> Self {
        SubstructMatchParameters {
            use_chirality: false,
            use_enhanced_stereo: false,
            aromatic_matches_conjugated: false,
            use_query_query_matches: false,
            use_generic_matchers: false,
            recursion_possible: true,
            uniquify: true,
            max_matches: 1000,
            num_threads: 1,
            atom_properties: Vec::new(),
            bond_properties: Vec::new(),
            max_recursive_matches: 1000,
            specified_stereo_query_matches_unspecified: false,
            aromatic_matches_single_or_double: false,
            extra_atom_check: None,
            extra_atom_check_overrides_default: false,
            extra_bond_check: None,
            extra_bond_check_overrides_default: false,
            extra_final_check: None,
        }
    }
}

const NULL_NODE: u32 = u32::MAX;

// ================================================================== VF2

/// What VF2 needs from a graph: adjacency lists of (neighbor, edge id) in the
/// order of the upstream `boost::adjacency_list`.
pub(crate) trait Vf2Graph {
    fn num_nodes(&self) -> u32;
    fn adjacency(&self, node: u32) -> &[(u32, u32)];

    fn degree(&self, node: u32) -> u32 {
        self.adjacency(node).len() as u32
    }

    /// The first edge from `a` to `b` (`boost::edge`).
    fn edge_between(&self, a: u32, b: u32) -> Option<u32> {
        self.adjacency(a)
            .iter()
            .find(|&&(n, _)| n == b)
            .map(|&(_, e)| e)
    }
}

impl Vf2Graph for ROMol {
    fn num_nodes(&self) -> u32 {
        self.num_atoms()
    }

    fn adjacency(&self, node: u32) -> &[(u32, u32)] {
        &self.adj[node as usize]
    }
}

/// The VF2 sub-state (`boost::detail::VF2SubState`, no node sorting).
struct Vf2<'a, G, VC, EC, MC> {
    g1: &'a G,
    g2: &'a G,
    vc: VC,
    ec: EC,
    mc: MC,
    n1: u32,
    n2: u32,
    core_len: u32,
    t1_len: u32,
    t2_len: u32,
    core_1: Vec<u32>,
    core_2: Vec<u32>,
    term_1: Vec<u32>,
    term_2: Vec<u32>,
}

/// `Pair` with the "VF2 Plus" neighbor iterator: (g2 vertex, position).
struct Pair {
    n1: u32,
    n2: u32,
    nbr: Option<(u32, usize)>,
}

impl<G: Vf2Graph, VC, EC, MC> Vf2<'_, G, VC, EC, MC>
where
    VC: FnMut(u32, u32) -> Result<bool>,
    EC: FnMut(u32, u32) -> Result<bool>,
    MC: FnMut(&[u32], &[u32]) -> Result<bool>,
{
    fn is_goal(&self) -> bool {
        self.core_len == self.n1
    }

    fn is_dead(&self) -> bool {
        self.n1 > self.n2 || self.t1_len > self.t2_len
    }

    fn next_pair(&self, pair: &mut Pair) -> bool {
        if pair.n1 == NULL_NODE {
            pair.n1 = 0;
        }
        if pair.n2 == NULL_NODE {
            pair.n2 = 0;
        } else {
            pair.n2 += 1;
        }
        if self.t1_len > self.core_len && self.t2_len > self.core_len {
            while pair.n1 < self.n1
                && (self.core_1[pair.n1 as usize] != NULL_NODE
                    || self.term_1[pair.n1 as usize] == 0)
            {
                pair.n1 += 1;
                pair.n2 = 0;
            }
            if pair.nbr.is_none() {
                let mapped = self
                    .g1
                    .adjacency(pair.n1)
                    .iter()
                    .map(|&(n, _)| n)
                    .find(|&n| self.core_1[n as usize] != NULL_NODE)
                    .expect("terminal node has a mapped neighbor");
                pair.nbr = Some((self.core_1[mapped as usize], 0));
            }
        } else {
            while pair.n1 < self.n1 && self.core_1[pair.n1 as usize] != NULL_NODE {
                pair.n1 += 1;
                pair.n2 = 0;
            }
        }
        if let Some((v, ref mut pos)) = pair.nbr {
            let adj = self.g2.adjacency(v);
            while *pos < adj.len() && self.core_2[adj[*pos].0 as usize] != NULL_NODE {
                *pos += 1;
            }
            if *pos < adj.len() {
                pair.n2 = adj[*pos].0;
                *pos += 1;
            } else {
                pair.n2 = self.n2;
            }
        } else if self.t1_len > self.core_len && self.t2_len > self.core_len {
            while pair.n2 < self.n2
                && (self.core_2[pair.n2 as usize] != NULL_NODE
                    || self.term_2[pair.n2 as usize] == 0)
            {
                pair.n2 += 1;
            }
        } else {
            while pair.n2 < self.n2 && self.core_2[pair.n2 as usize] != NULL_NODE {
                pair.n2 += 1;
            }
        }
        pair.n1 < self.n1 && pair.n2 < self.n2
    }

    fn is_feasible_pair(&mut self, node1: u32, node2: u32) -> Result<bool> {
        if self.g1.degree(node1) > self.g2.degree(node2) {
            return Ok(false);
        }
        if !(self.vc)(node1, node2)? {
            return Ok(false);
        }
        for &(other1, e1) in self.g1.adjacency(node1) {
            let other2 = self.core_1[other1 as usize];
            if other2 != NULL_NODE {
                match self.g2.edge_between(node2, other2) {
                    Some(e2) => {
                        if !(self.ec)(e1, e2)? {
                            return Ok(false);
                        }
                    }
                    None => return Ok(false),
                }
            }
        }
        Ok(true)
    }

    fn add_pair(&mut self, node1: u32, node2: u32) {
        self.core_len += 1;
        let cl = self.core_len;
        if self.term_1[node1 as usize] == 0 {
            self.term_1[node1 as usize] = cl;
            self.t1_len += 1;
        }
        if self.term_2[node2 as usize] == 0 {
            self.term_2[node2 as usize] = cl;
            self.t2_len += 1;
        }
        self.core_1[node1 as usize] = node2;
        self.core_2[node2 as usize] = node1;
        for &(o, _) in self.g1.adjacency(node1) {
            if self.term_1[o as usize] == 0 {
                self.term_1[o as usize] = cl;
                self.t1_len += 1;
            }
        }
        for &(o, _) in self.g2.adjacency(node2) {
            if self.term_2[o as usize] == 0 {
                self.term_2[o as usize] = cl;
                self.t2_len += 1;
            }
        }
    }

    fn back_track(&mut self, node1: u32, node2: u32) {
        let cl = self.core_len;
        if self.term_1[node1 as usize] == cl {
            self.term_1[node1 as usize] = 0;
            self.t1_len -= 1;
        }
        for &(o, _) in self.g1.adjacency(node1) {
            if self.term_1[o as usize] == cl {
                self.term_1[o as usize] = 0;
                self.t1_len -= 1;
            }
        }
        if self.term_2[node2 as usize] == cl {
            self.term_2[node2 as usize] = 0;
            self.t2_len -= 1;
        }
        for &(o, _) in self.g2.adjacency(node2) {
            if self.term_2[o as usize] == cl {
                self.term_2[o as usize] = 0;
                self.t2_len -= 1;
            }
        }
        self.core_1[node1 as usize] = NULL_NODE;
        self.core_2[node2 as usize] = NULL_NODE;
        self.core_len -= 1;
    }

    fn core_set(&self) -> (Vec<u32>, Vec<u32>) {
        let mut c1 = Vec::with_capacity(self.core_len as usize);
        let mut c2 = Vec::with_capacity(self.core_len as usize);
        for i in 0..self.n1 {
            if self.core_1[i as usize] != NULL_NODE {
                c1.push(i);
                c2.push(self.core_1[i as usize]);
            }
        }
        (c1, c2)
    }

    /// `MatchAll`: returns true when the caller must stop.
    fn match_all(
        &mut self,
        res: &mut Vec<Vec<(u32, u32)>>,
        lim: u32,
        count_only: bool,
    ) -> Result<bool> {
        if self.is_goal() {
            let (c1, c2) = self.core_set();
            if (self.mc)(&c1, &c2)? {
                if count_only {
                    res.push(Vec::new());
                } else {
                    res.push(c1.into_iter().zip(c2).collect());
                }
                return Ok(lim != 0 && res.len() >= lim as usize);
            }
        }
        if self.is_dead() {
            return Ok(false);
        }
        let mut pair = Pair {
            n1: NULL_NODE,
            n2: NULL_NODE,
            nbr: None,
        };
        while self.next_pair(&mut pair) {
            if self.is_feasible_pair(pair.n1, pair.n2)? {
                self.add_pair(pair.n1, pair.n2);
                if self.match_all(res, lim, count_only)? {
                    return Ok(true);
                }
                self.back_track(pair.n1, pair.n2);
            }
        }
        Ok(false)
    }
}

/// `boost::vf2_all(g1, g2, vc, ec, mc, F, max_results)`
pub(crate) fn vf2_all<G, VC, EC, MC>(
    g1: &G,
    g2: &G,
    vc: VC,
    ec: EC,
    mc: MC,
    max_results: u32,
    count_only: bool,
) -> Result<Vec<Vec<(u32, u32)>>>
where
    G: Vf2Graph,
    VC: FnMut(u32, u32) -> Result<bool>,
    EC: FnMut(u32, u32) -> Result<bool>,
    MC: FnMut(&[u32], &[u32]) -> Result<bool>,
{
    let (n1, n2) = (g1.num_nodes(), g2.num_nodes());
    let mut s = Vf2 {
        g1,
        g2,
        vc,
        ec,
        mc,
        n1,
        n2,
        core_len: 0,
        t1_len: 0,
        t2_len: 0,
        core_1: vec![NULL_NODE; n1 as usize],
        core_2: vec![NULL_NODE; n2 as usize],
        term_1: vec![0; n1 as usize],
        term_2: vec![0; n2 as usize],
    };
    let mut res = Vec::new();
    s.match_all(&mut res, max_results, count_only)?;
    Ok(res)
}

// ================================================================== compatibility

/// `Atom::Match` / `QueryAtom::Match` of query atom `qi` against `mi`.
fn atom_match(query: &ROMol, qi: u32, mol: &ROMol, mi: u32, rec: &RecursionSets) -> Result<bool> {
    let qa = query.atom(qi);
    match qa.query() {
        Some(q) => q.matches_with(mol, mi, Some(rec)),
        None => Ok(qa.matches(mol.atom(mi))),
    }
}

/// `Bond::Match` / `QueryBond::Match`
fn bond_match(query: &ROMol, qb: u32, mol: &ROMol, mb: u32, rec: &RecursionSets) -> Result<bool> {
    let b = query.bond(qb);
    match b.query() {
        Some(q) => q.matches_with(mol, mb, Some(rec)),
        None => Ok(b.matches(mol.bond(mb))),
    }
}

/// `queriesMatch` (QueryBond.cpp)
fn bond_queries_match(q1: &BondQuery, q2: &BondQuery) -> bool {
    const EQ: [&str; 6] = [
        "BondRingSize",
        "BondMinRingSize",
        "BondOrder",
        "BondDir",
        "BondInRing",
        "BondInNRings",
    ];
    let (d1, d2) = (q1.description.as_str(), q2.description.as_str());
    if d1 == "BondNull" || d2 == "BondNull" {
        true
    } else if d1 == "BondOr" {
        q1.children.iter().any(|c1| {
            if d2 == "BondOr" {
                q2.children.iter().any(|c2| bond_queries_match(c1, c2))
            } else {
                bond_queries_match(c1, q2)
            }
        })
    } else if d1 == "BondAnd" {
        q1.children.iter().all(|c1| {
            if d2 == "BondAnd" {
                q2.children.iter().any(|c2| bond_queries_match(c1, c2))
            } else {
                bond_queries_match(c1, q2)
            }
        })
    } else if d2 == "BondOr" {
        q2.children.iter().any(|c2| bond_queries_match(q1, c2))
    } else if d2 == "BondAnd" {
        // note: upstream inverts the child results here
        !q2.children.iter().any(|c2| bond_queries_match(q1, c2))
    } else if EQ.contains(&d1) {
        if q1.negate == q2.negate {
            q1.val == q2.val
        } else {
            q1.val != q2.val
        }
    } else {
        false
    }
}

/// `propertyCompat`
fn property_compat(p1: &impl RDProps, p2: &impl RDProps, properties: &[String]) -> bool {
    let as_str = |v: &RDValue| v.to_rd_string();
    for prop in properties {
        match (p1.get_prop(prop), p2.get_prop(prop)) {
            (Some(a), Some(b)) => {
                if as_str(a) != as_str(b) {
                    return false;
                }
            }
            (None, None) => {}
            _ => return false,
        }
    }
    true
}

struct Ctx<'a> {
    query: &'a ROMol,
    mol: &'a ROMol,
    ps: &'a SubstructMatchParameters,
    rec: &'a RecursionSets,
}

fn has_chiral_label(t: ChiralType) -> bool {
    matches!(t, ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW)
}

impl Ctx<'_> {
    /// `atomCompat`
    fn atom_compat(&self, qi: u32, mi: u32) -> Result<bool> {
        let (qa, ma) = (self.query.atom(qi), self.mol.atom(mi));
        let extra = self.ps.extra_atom_check.as_ref();
        if let (true, Some(check)) = (self.ps.extra_atom_check_overrides_default, extra) {
            return Ok(check.call(self.query, qi, self.mol, mi));
        }
        let res = if self.ps.use_query_query_matches && qa.has_query() && ma.has_query() {
            match (qa.query(), ma.query()) {
                (Some(q1), Some(q2)) => atom_queries_match(q1, q2),
                _ => false,
            }
        } else {
            atom_match(self.query, qi, self.mol, mi, self.rec)?
        };
        if !res {
            return Ok(false);
        }
        if !self.ps.atom_properties.is_empty() && !property_compat(qa, ma, &self.ps.atom_properties)
        {
            return Ok(false);
        }
        if extra.is_some_and(|check| !check.call(self.query, qi, self.mol, mi)) {
            return Ok(false);
        }
        Ok(true)
    }

    /// `AtomLabelFunctor`
    fn atom_label(&self, qi: u32, mi: u32) -> Result<bool> {
        if self.ps.use_chirality
            && has_chiral_label(self.query.atom(qi).chiral_tag())
            && !self.ps.specified_stereo_query_matches_unspecified
            && !has_chiral_label(self.mol.atom(mi).chiral_tag())
        {
            return Ok(false);
        }
        self.atom_compat(qi, mi)
    }

    /// `bondCompat`
    fn bond_compat(&self, qb: u32, mb: u32) -> Result<bool> {
        let (b1, b2) = (self.query.bond(qb), self.mol.bond(mb));
        let extra = self.ps.extra_bond_check.as_ref();
        if let (true, Some(check)) = (self.ps.extra_bond_check_overrides_default, extra) {
            return Ok(check.call(self.query, qb, self.mol, mb));
        }
        let conj_sd = |b: &crate::graphmol::bond::Bond| {
            b.is_conjugated() && matches!(b.bond_type(), BondType::Single | BondType::Double)
        };
        let sd = |b: &crate::graphmol::bond::Bond| {
            matches!(b.bond_type(), BondType::Single | BondType::Double)
        };
        let (t1, t2) = (b1.bond_type(), b2.bond_type());
        let res = if self.ps.use_query_query_matches && b1.has_query() && b2.has_query() {
            match (b1.query(), b2.query()) {
                (Some(q1), Some(q2)) => bond_queries_match(q1, q2),
                _ => false,
            }
        } else if self.ps.aromatic_matches_conjugated
            && !b1.has_query()
            && !b2.has_query()
            && ((t1 == BondType::Aromatic && t2 == BondType::Aromatic)
                || (t1 == BondType::Aromatic && conj_sd(b2))
                || (t2 == BondType::Aromatic && conj_sd(b1)))
        {
            true
        } else if self.ps.aromatic_matches_single_or_double
            && !b1.has_query()
            && !b2.has_query()
            && ((t1 == BondType::Aromatic && t2 == BondType::Aromatic)
                || (t1 == BondType::Aromatic && sd(b2))
                || (t2 == BondType::Aromatic && sd(b1)))
        {
            true
        } else {
            bond_match(self.query, qb, self.mol, mb, self.rec)?
        };
        if !res {
            return Ok(false);
        }
        if t1 == BondType::Dative && t2 == BondType::Dative {
            // the query atom is matched against the molecule atom
            if !atom_match(
                self.query,
                b1.begin_atom_idx(),
                self.mol,
                b2.begin_atom_idx(),
                self.rec,
            )? || !atom_match(
                self.query,
                b1.end_atom_idx(),
                self.mol,
                b2.end_atom_idx(),
                self.rec,
            )? {
                return Ok(false);
            }
        }
        if !self.ps.bond_properties.is_empty() && !property_compat(b1, b2, &self.ps.bond_properties)
        {
            return Ok(false);
        }
        if extra.is_some_and(|check| !check.call(self.query, qb, self.mol, mb)) {
            return Ok(false);
        }
        Ok(true)
    }

    /// `BondLabelFunctor`
    fn bond_label(&self, qb: u32, mb: u32) -> Result<bool> {
        if self.ps.use_chirality {
            let q = self.query.bond(qb);
            if q.bond_type() == BondType::Double && q.stereo() > BondStereo::StereoAny {
                let m = self.mol.bond(mb);
                if m.bond_type() == BondType::Double
                    && !self.ps.specified_stereo_query_matches_unspecified
                    && m.stereo() <= BondStereo::StereoAny
                {
                    return Ok(false);
                }
            }
        }
        self.bond_compat(qb, mb)
    }
}

/// `Chirality::translateEZLabelToCisTrans`
fn ez_to_cis_trans(s: BondStereo) -> BondStereo {
    match s {
        BondStereo::StereoE => BondStereo::StereoTrans,
        BondStereo::StereoZ => BondStereo::StereoCis,
        other => other,
    }
}

/// `MolMatchFinalCheckFunctor`
struct FinalCheck<'a> {
    query: &'a ROMol,
    mol: &'a ROMol,
    ps: &'a SubstructMatchParameters,
    mol_stereo_groups: BTreeMap<u32, usize>,
    matches_seen: HashSet<Vec<bool>>,
}

impl<'a> FinalCheck<'a> {
    fn new(query: &'a ROMol, mol: &'a ROMol, ps: &'a SubstructMatchParameters) -> Self {
        let mut mol_stereo_groups = BTreeMap::new();
        if ps.use_enhanced_stereo {
            for (gi, sg) in mol.stereo_groups().iter().enumerate() {
                if sg.group_type() == StereoGroupType::Absolute {
                    continue;
                }
                for &a in sg.atoms() {
                    mol_stereo_groups.insert(a, gi);
                }
            }
        }
        FinalCheck {
            query,
            mol,
            ps,
            mol_stereo_groups,
            matches_seen: HashSet::new(),
        }
    }

    fn check(&mut self, q_c: &[u32], m_c: &[u32]) -> Result<bool> {
        let nq = self.query.num_atoms() as usize;
        if let Some(FinalMatchCheck(f)) = &self.ps.extra_final_check {
            if !f(self.mol, &m_c[..nq]) {
                return Ok(false);
            }
        }
        let mut seen_key = Vec::new();
        if self.ps.uniquify {
            seen_key = vec![false; self.mol.num_atoms() as usize];
            for &m in &m_c[..nq] {
                seen_key[m as usize] = true;
            }
            if self.matches_seen.contains(&seen_key) {
                return Ok(false);
            }
        }
        if !self.ps.use_chirality {
            if self.ps.uniquify {
                self.matches_seen.insert(seen_key);
            }
            return Ok(true);
        }
        let mut matches: BTreeMap<u32, bool> = BTreeMap::new();
        for i in 0..nq {
            let q_at = q_c[i];
            let qa = self.query.atom(q_at);
            let q_deg = self.query.atom_degree(q_at);
            if q_deg < 3 || !has_chiral_label(qa.chiral_tag()) {
                continue;
            }
            let m_at = m_c[i];
            let ma = self.mol.atom(m_at);
            if !has_chiral_label(ma.chiral_tag()) {
                if self.ps.specified_stereo_query_matches_unspecified {
                    continue;
                }
                return Ok(false);
            }
            if q_deg > self.mol.atom_degree(m_at) {
                return Ok(false);
            }
            let mut q_order: Vec<i32> = Vec::new();
            let mut m_order: Vec<i32> = Vec::new();
            for j in 0..nq {
                let qb = self.query.bond_idx_between(q_at, q_c[j]);
                let mb = self.mol.bond_idx_between(m_at, m_c[j]);
                if let (Some(qb), Some(mb)) = (qb, mb) {
                    m_order.push(mb as i32);
                    q_order.push(qb as i32);
                    if m_order.len() as u32 == q_deg {
                        break;
                    }
                }
            }
            if q_order.len() as u32 != q_deg {
                return Err(RDKitError::invariant("missing matches"));
            }
            let q_perm = self.query.perturbation_order(q_at, &q_order)?;
            let unmatched = self.mol.atom_degree(m_at) as usize - m_order.len();
            m_order.extend(std::iter::repeat_n(-1, unmatched));
            let mo_order: Vec<i32> = self
                .mol
                .atom_bond_indices(m_at)
                .map(|b| {
                    if m_order.contains(&(b as i32)) {
                        b as i32
                    } else {
                        -1
                    }
                })
                .collect();
            let m_perm = count_swaps_to_interconvert(&mo_order, &m_order)? as i32;
            let require_match = q_perm % 2 == m_perm % 2;
            let labels_match = qa.chiral_tag() == ma.chiral_tag();
            let match_ok = require_match == labels_match;
            if self.mol_stereo_groups.contains_key(&m_at) {
                matches.insert(m_at, match_ok);
            } else if !match_ok {
                return Ok(false);
            }
        }
        let mut q_to_mol: BTreeMap<u32, u32> = BTreeMap::new();
        for j in 0..nq {
            q_to_mol.insert(q_c[j], m_c[j]);
        }
        if self.ps.use_enhanced_stereo && !self.enhanced_stereo_ok(&q_to_mol, &matches) {
            return Ok(false);
        }
        for qb in self.query.bonds() {
            if qb.bond_type() != BondType::Double || qb.stereo() <= BondStereo::StereoAny {
                continue;
            }
            if qb.stereo_atoms().len() != 2 {
                continue;
            }
            let m_begin = q_to_mol[&qb.begin_atom_idx()];
            let m_end = q_to_mol[&qb.end_atom_idx()];
            let mb = self
                .mol
                .bond_idx_between(m_begin, m_end)
                .ok_or_else(|| RDKitError::invariant("Matching bond not found"))?;
            let mb = self.mol.bond(mb);
            if mb.bond_type() != BondType::Double {
                continue;
            }
            if !self.ps.specified_stereo_query_matches_unspecified
                && mb.stereo() <= BondStereo::StereoAny
            {
                return Ok(false);
            }
            if mb.stereo_atoms().len() != 2 {
                continue;
            }
            let qsa = qb.stereo_atoms();
            let msa = mb.stereo_atoms();
            let qm = |a: i32| q_to_mol.get(&(a as u32)).map_or(u32::MAX, |&v| v);
            let (e1, e2) = if m_begin == mb.begin_atom_idx() {
                (
                    (qm(qsa[0]) == msa[0] as u32) as u32,
                    (qm(qsa[1]) == msa[1] as u32) as u32,
                )
            } else {
                (
                    (qm(qsa[0]) == msa[1] as u32) as u32,
                    (qm(qsa[1]) == msa[0] as u32) as u32,
                )
            };
            let total = e1 + e2;
            let (ms, qs) = (ez_to_cis_trans(mb.stereo()), ez_to_cis_trans(qb.stereo()));
            if ms == qs && total == 1 {
                return Ok(false);
            }
            if ms != qs && total != 1 {
                return Ok(false);
            }
        }
        if self.ps.uniquify {
            self.matches_seen.insert(seen_key);
        }
        Ok(true)
    }

    /// `enhancedStereoIsOK`
    fn enhanced_stereo_ok(
        &self,
        q_to_mol: &BTreeMap<u32, u32>,
        matches: &BTreeMap<u32, bool>,
    ) -> bool {
        let mut mol_atoms_to_query_groups: BTreeMap<u32, usize> = BTreeMap::new();
        for (qgi, sg) in self.query.stereo_groups().iter().enumerate() {
            if sg.group_type() == StereoGroupType::Absolute {
                continue;
            }
            let is_and = sg.group_type() == StereoGroupType::And;
            for &a in sg.atoms() {
                let m = q_to_mol.get(&a).copied().unwrap_or(0);
                match self.mol_stereo_groups.get(&m) {
                    None => return false,
                    Some(&mg) => {
                        if is_and
                            && self.mol.stereo_groups()[mg].group_type() != StereoGroupType::And
                        {
                            return false;
                        }
                    }
                }
                mol_atoms_to_query_groups.insert(m, qgi);
            }
        }
        for sg in self.mol.stereo_groups() {
            if sg.group_type() == StereoGroupType::Absolute {
                continue;
            }
            let mut does_match = false;
            let mut seen = false;
            let mut q_group: Option<usize> = None;
            for &a in sg.atoms() {
                let Some(&this_does_match) = matches.get(&a) else {
                    continue;
                };
                let this_q_group = mol_atoms_to_query_groups.get(&a).copied();
                if !seen {
                    does_match = this_does_match;
                    q_group = this_q_group;
                    seen = true;
                } else if does_match != this_does_match || this_q_group != q_group {
                    return false;
                }
            }
        }
        true
    }
}

// ================================================================== recursion

/// `detail::MatchSubqueries`: fills the result sets of every recursive query
/// below `query` (reached through `ptr`, the node's address in its owner).
fn match_subqueries(
    mol: &ROMol,
    query: &AtomQuery,
    ps: &SubstructMatchParameters,
    subquery_map: &mut BTreeMap<u32, BTreeSet<i32>>,
    rec: &mut RecursionSets,
) -> Result<()> {
    if let QueryKind::Recursive {
        mol: qmol,
        serial_number,
        ..
    } = &query.kind
    {
        let key = query.node_key();
        let mut set = BTreeSet::new();
        let mut done = false;
        if *serial_number != 0 {
            if let Some(s) = subquery_map.get(serial_number) {
                set = s.clone();
                done = true;
            }
        }
        if !done {
            let starts = recursive_matcher(mol, &qmol.0, subquery_map, ps, rec)?;
            set.extend(starts);
            if *serial_number != 0 {
                subquery_map.insert(*serial_number, set.clone());
            }
        }
        rec.insert(key, set);
    }
    for c in &query.children {
        match_subqueries(mol, c, ps, subquery_map, rec)?;
    }
    Ok(())
}

/// `detail::RecursiveMatcher`: the molecule atoms matched by the first (or
/// `_queryRootAtom`) query atom.
fn recursive_matcher(
    mol: &ROMol,
    query: &ROMol,
    subquery_map: &mut BTreeMap<u32, BTreeSet<i32>>,
    ps: &SubstructMatchParameters,
    rec: &mut RecursionSets,
) -> Result<Vec<i32>> {
    let lps = SubstructMatchParameters {
        max_matches: ps.max_recursive_matches.max(ps.max_matches),
        uniquify: false,
        ..ps.clone()
    };
    for i in 0..query.num_atoms() {
        if let Some(q) = query.atom(i).query() {
            match_subqueries(mol, q, &lps, subquery_map, rec)?;
        }
    }
    let pms = run_vf2(mol, query, &lps, rec, false)?;
    let root = query.get_prop(cp::_queryRootAtom).and_then(RDValue::as_int);
    let mut matches = Vec::new();
    for pairs in pms {
        match root {
            None => {
                if let Some(&(_, m)) = pairs.first() {
                    matches.push(m as i32);
                }
            }
            Some(r) => match pairs.iter().find(|&&(q, _)| q == r as u32) {
                Some(&(_, m)) => matches.push(m as i32),
                None => crate::rd_log_error!("no match found for queryRootAtom"),
            },
        }
        if matches.len() == lps.max_matches as usize {
            break;
        }
    }
    Ok(matches)
}

/// The atoms of `mol` that match the query atom `qa`, which is what
/// iterating a `QueryAtomIterator` (`Mol.GetAtomsMatchingQuery`) yields.
/// Recursive queries are evaluated against `mol` first.
pub fn atoms_matching_query(mol: &ROMol, qa: &Atom) -> Result<Vec<u32>> {
    let Some(query) = qa.query() else {
        return Ok((0..mol.num_atoms())
            .filter(|&i| qa.matches(mol.atom(i)))
            .collect());
    };
    let mut rec = RecursionSets::new();
    if has_recursive_query(Some(query)) {
        let ps = SubstructMatchParameters::default();
        match_subqueries(mol, query, &ps, &mut BTreeMap::new(), &mut rec)?;
    }
    let mut matching = Vec::new();
    for idx in 0..mol.num_atoms() {
        if query.matches_with(mol, idx, Some(&rec))? {
            matching.push(idx);
        }
    }
    Ok(matching)
}

fn run_vf2(
    mol: &ROMol,
    query: &ROMol,
    ps: &SubstructMatchParameters,
    rec: &RecursionSets,
    count_only: bool,
) -> Result<Vec<Vec<(u32, u32)>>> {
    let ctx = Ctx {
        query,
        mol,
        ps,
        rec,
    };
    let mut fc = FinalCheck::new(query, mol, ps);
    vf2_all(
        query,
        mol,
        |i, j| ctx.atom_label(i, j),
        |i, j| ctx.bond_label(i, j),
        |c1, c2| fc.check(c1, c2),
        ps.max_matches,
        count_only,
    )
}

fn prepare_recursion(
    mol: &ROMol,
    query: &ROMol,
    ps: &SubstructMatchParameters,
) -> Result<RecursionSets> {
    let mut rec = RecursionSets::new();
    if ps.recursion_possible {
        let mut subquery_map = BTreeMap::new();
        for i in 0..query.num_atoms() {
            if let Some(q) = query.atom(i).query() {
                match_subqueries(mol, q, ps, &mut subquery_map, &mut rec)?;
            }
        }
    }
    Ok(rec)
}

fn check_params(ps: &SubstructMatchParameters) -> Result<()> {
    if ps.use_generic_matchers {
        return Err(RDKitError::Value(
            "useGenericMatchers is not supported yet".into(),
        ));
    }
    Ok(())
}

/// `SubstructMatch(mol, query, params)`
pub fn substruct_match(
    mol: &ROMol,
    query: &ROMol,
    ps: &SubstructMatchParameters,
) -> Result<Vec<MatchVect>> {
    check_params(ps)?;
    let (mn, qn) = (mol.num_atoms(), query.num_atoms());
    if mn == 0 || qn == 0 || qn > mn {
        return Ok(Vec::new());
    }
    let rec = prepare_recursion(mol, query, ps)?;
    let pms = run_vf2(mol, query, ps, &rec, false)?;
    let mut matches = Vec::with_capacity(pms.len());
    let mut match_vect: MatchVect = vec![(0, 0); qn as usize];
    for pairs in pms {
        for (q, m) in pairs {
            match_vect[q as usize] = (q as i32, m as i32);
        }
        matches.push(match_vect.clone());
    }
    Ok(matches)
}

/// `SubstructMatchCount(mol, query, params)`
pub fn substruct_match_count(
    mol: &ROMol,
    query: &ROMol,
    ps: &SubstructMatchParameters,
) -> Result<u32> {
    check_params(ps)?;
    if mol.num_atoms() == 0 || query.num_atoms() == 0 {
        return Ok(0);
    }
    let rec = prepare_recursion(mol, query, ps)?;
    Ok(run_vf2(mol, query, ps, &rec, true)?.len() as u32)
}

/// `removeDuplicates(matches, nAtoms)`
pub fn remove_duplicates(matches: &mut Vec<MatchVect>, n_atoms: u32) {
    let mut seen: HashSet<Vec<bool>> = HashSet::new();
    matches.retain(|m| {
        let mut v = vec![false; n_atoms as usize];
        for &(_, a) in m {
            v[a as usize] = true;
        }
        seen.insert(v)
    });
}

impl SubstructMatchParameters {
    /// `substructMatchParamsToJSON`: the serializable options as JSON (callbacks
    /// and property lists are not included). Unlike upstream, the text is compact.
    pub fn to_json(&self) -> String {
        let mut pt = PTree::new();
        pt.put("useChirality", self.use_chirality);
        pt.put("useEnhancedStereo", self.use_enhanced_stereo);
        pt.put(
            "aromaticMatchesConjugated",
            self.aromatic_matches_conjugated,
        );
        pt.put("useQueryQueryMatches", self.use_query_query_matches);
        pt.put("recursionPossible", self.recursion_possible);
        pt.put("uniquify", self.uniquify);
        pt.put("maxMatches", self.max_matches);
        pt.put("maxRecursiveMatches", self.max_recursive_matches);
        pt.put("numThreads", self.num_threads);
        pt.put(
            "specifiedStereoQueryMatchesUnspecified",
            self.specified_stereo_query_matches_unspecified,
        );
        pt.put(
            "aromaticMatchesSingleOrDouble",
            self.aromatic_matches_single_or_double,
        );
        pt.to_json()
    }

    /// `updateSubstructMatchParamsFromJSON`: options missing from `json` keep
    /// their current value; an empty string changes nothing.
    pub fn update_from_json(&mut self, json: &str) -> crate::Result<()> {
        if json.is_empty() {
            return Ok(());
        }
        let pt = PTree::from_json(json)?;
        self.use_chirality = pt.get_bool("useChirality", self.use_chirality);
        self.use_enhanced_stereo = pt.get_bool("useEnhancedStereo", self.use_enhanced_stereo);
        self.aromatic_matches_conjugated = pt.get_bool(
            "aromaticMatchesConjugated",
            self.aromatic_matches_conjugated,
        );
        self.use_query_query_matches =
            pt.get_bool("useQueryQueryMatches", self.use_query_query_matches);
        self.recursion_possible = pt.get_bool("recursionPossible", self.recursion_possible);
        self.uniquify = pt.get_bool("uniquify", self.uniquify);
        self.max_matches = pt.get_u32("maxMatches", self.max_matches);
        self.max_recursive_matches = pt.get_u32("maxRecursiveMatches", self.max_recursive_matches);
        self.num_threads = pt.get_i32("numThreads", self.num_threads);
        self.specified_stereo_query_matches_unspecified = pt.get_bool(
            "specifiedStereoQueryMatchesUnspecified",
            self.specified_stereo_query_matches_unspecified,
        );
        self.aromatic_matches_single_or_double = pt.get_bool(
            "aromaticMatchesSingleOrDouble",
            self.aromatic_matches_single_or_double,
        );
        Ok(())
    }
}
