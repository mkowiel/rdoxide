//! Atom and bond queries.
//!
//! Ported from RDKit `Code/Query/*.h` (Query, EqualityQuery, Less/Greater
//! queries, RangeQuery, SetQuery, And/Or/XOrQuery, NullQueryAlgebra) and
//! `Code/GraphMol/QueryOps.{h,cpp}`, `QueryAtom.{h,cpp}`, `QueryBond.{h,cpp}`
//! (2026.09.1).
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Upstream's class hierarchy (`Query<int, Atom const*, true>` and
//! subclasses, with `std::function` data functions) becomes one [`Query`]
//! struct with a [`QueryKind`] and a data-function enum ([`AtomFn`] /
//! [`BondFn`]). Descriptions, type labels and `getFullDescription` text are
//! kept identical, since SMARTS writing and `DescribeQuery` depend on them.

use std::collections::BTreeSet;
use std::fmt::Debug;
use std::sync::Arc;

use crate::graphmol::atom::{ChiralType, ValenceType};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::query_ops::query_is_atom_bridgehead;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `Queries::CompositeQueryType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompositeQueryType {
    And,
    Or,
    XOr,
    OrSequence,
}

/// `queryCmp(v1, v2, tol)`
fn query_cmp(v1: i32, v2: i32, tol: i32) -> i32 {
    let diff = v1.wrapping_sub(v2);
    if diff <= tol {
        if diff >= -tol { 0 } else { -1 }
    } else {
        1
    }
}

fn query_cmp_f(v1: f64, v2: f64, tol: f64) -> i32 {
    let diff = v1 - v2;
    if diff <= tol {
        if diff >= -tol { 0 } else { -1 }
    } else {
        1
    }
}

/// The value of a `HasPropWithValue` query.
#[derive(Debug, Clone, PartialEq)]
pub enum PropValue {
    Int(i32),
    UInt(u32),
    Bool(bool),
    Float(f32),
    Double(f64),
    String(String),
}

/// The query molecule of a recursive (`$(...)`) query; compared by identity.
#[derive(Debug, Clone)]
pub struct QueryMol(pub Arc<ROMol>);

impl PartialEq for QueryMol {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

/// What a query node does (upstream: the concrete `Query` subclass).
#[derive(Debug, Clone, PartialEq)]
pub enum QueryKind {
    /// base `Query`: `matchFunc(dataFunc(x))`, or `bool(dataFunc(x))`
    Base {
        always_true: bool,
    },
    Equality,
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    Range {
        lower: i32,
        upper: i32,
        lower_open: bool,
        upper_open: bool,
    },
    Set(BTreeSet<i32>),
    And,
    Or,
    XOr,
    /// `AtomRingQuery`: a negative value means "in any ring"
    AtomRing,
    /// `RecursiveStructureQuery`
    Recursive {
        mol: QueryMol,
        serial_number: u32,
        set: BTreeSet<i32>,
        initialized: bool,
    },
    /// `HasPropQuery`
    HasProp(String),
    /// `HasPropWithValueQuery`
    HasPropWithValue {
        prop: String,
        val: PropValue,
        tol: f64,
    },
}

/// A query data function (`setDataFunc`): extracts an `int` from the target.
pub trait QueryTarget: Debug + Clone + PartialEq + Send + Sync + 'static {
    fn eval(&self, mol: &ROMol, idx: u32) -> Result<i32>;
    fn props(mol: &ROMol, idx: u32) -> &crate::rdgeneral::Dict;
}

/// A query on atoms or bonds (`Queries::Query<int, T const*, true>`).
#[derive(Debug, Clone, PartialEq)]
pub struct Query<F: QueryTarget> {
    pub kind: QueryKind,
    pub description: String,
    pub type_label: String,
    pub negate: bool,
    pub val: i32,
    pub tol: i32,
    pub children: Vec<Query<F>>,
    pub data_func: Option<F>,
}

pub type AtomQuery = Query<AtomFn>;

/// Match results of recursive queries, keyed by [`Query::node_key`].
pub type RecursionSets = std::collections::HashMap<usize, BTreeSet<i32>>;
pub type BondQuery = Query<BondFn>;

impl<F: QueryTarget> Query<F> {
    pub fn new(kind: QueryKind) -> Self {
        Query {
            kind,
            description: String::new(),
            type_label: String::new(),
            negate: false,
            val: 0,
            tol: 0,
            children: Vec::new(),
            data_func: None,
        }
    }

    fn with(kind: QueryKind, val: i32, func: F, descr: &str) -> Self {
        let mut q = Query::new(kind);
        q.val = val;
        q.data_func = Some(func);
        q.description = descr.to_string();
        q
    }

    pub fn description(&self) -> &str {
        &self.description
    }

    fn data(&self, mol: &ROMol, idx: u32) -> Result<i32> {
        match &self.data_func {
            Some(f) => f.eval(mol, idx),
            None => Err(RDKitError::precondition("no data function")),
        }
    }

    /// `Match(what)`
    pub fn matches(&self, mol: &ROMol, idx: u32) -> Result<bool> {
        self.matches_with(mol, idx, None)
    }

    /// `Match(what)` where the sets of recursive (`$(...)`) queries come
    /// from `recursion` (keyed by the address of the query node) when given;
    /// upstream fills those sets in place during substructure matching.
    pub fn matches_with(
        &self,
        mol: &ROMol,
        idx: u32,
        recursion: Option<&RecursionSets>,
    ) -> Result<bool> {
        let cmp_val = |q: &Self| -> Result<i32> { Ok(query_cmp(q.val, q.data(mol, idx)?, q.tol)) };
        let res = match &self.kind {
            QueryKind::Base { always_true } => {
                let v = self.data(mol, idx)?;
                if *always_true { true } else { v != 0 }
            }
            QueryKind::Equality => cmp_val(self)? == 0,
            QueryKind::Less => cmp_val(self)? < 0,
            QueryKind::LessEqual => cmp_val(self)? <= 0,
            QueryKind::Greater => cmp_val(self)? > 0,
            QueryKind::GreaterEqual => cmp_val(self)? >= 0,
            QueryKind::Range {
                lower,
                upper,
                lower_open,
                upper_open,
            } => {
                let v = self.data(mol, idx)?;
                let l = query_cmp(*lower, v, self.tol);
                let u = query_cmp(*upper, v, self.tol);
                let lr = if *lower_open { l < 0 } else { l <= 0 };
                let ur = if *upper_open { u > 0 } else { u >= 0 };
                lr && ur
            }
            QueryKind::Set(set) => set.contains(&self.data(mol, idx)?),
            QueryKind::Recursive { set, .. } => {
                let v = self.data(mol, idx)?;
                match recursion {
                    Some(r) => r.get(&self.node_key()).is_some_and(|s| s.contains(&v)),
                    None => set.contains(&v),
                }
            }
            QueryKind::And => {
                let mut r = true;
                for c in &self.children {
                    if !c.matches_with(mol, idx, recursion)? {
                        r = false;
                        break;
                    }
                }
                r
            }
            QueryKind::Or => {
                let mut r = false;
                for c in &self.children {
                    if c.matches_with(mol, idx, recursion)? {
                        r = true;
                        break;
                    }
                }
                r
            }
            QueryKind::XOr => {
                let mut r = false;
                for c in &self.children {
                    if c.matches_with(mol, idx, recursion)? {
                        if r {
                            r = false;
                            break;
                        }
                        r = true;
                    }
                }
                r
            }
            QueryKind::AtomRing => {
                let v = self.data(mol, idx)?;
                if self.val < 0 {
                    v != 0
                } else {
                    query_cmp(v, self.val, self.tol) == 0
                }
            }
            QueryKind::HasProp(name) => F::props(mol, idx).has_val(name),
            QueryKind::HasPropWithValue { prop, val, tol } => match F::props(mol, idx).get(prop) {
                None => false,
                Some(v) => prop_value_matches(v, val, *tol),
            },
        };
        Ok(res != self.negate)
    }

    /// Identity of this query node (used to key recursive-query results).
    pub fn node_key(&self) -> usize {
        self as *const Self as usize
    }

    /// `getFullDescription()`
    pub fn full_description(&self) -> String {
        let d = &self.description;
        match &self.kind {
            QueryKind::Equality
            | QueryKind::AtomRing
            | QueryKind::HasProp(_)
            | QueryKind::HasPropWithValue { .. } => {
                format!(
                    "{d} {} {} val",
                    self.val,
                    if self.negate { "!=" } else { "=" }
                )
            }
            QueryKind::Less => format!(
                "{d} {}{}",
                self.val,
                if self.negate { " ! < " } else { " < " }
            ),
            QueryKind::LessEqual => format!(
                "{d} {}{}",
                self.val,
                if self.negate { " ! <= " } else { " <= " }
            ),
            QueryKind::Greater => format!(
                "{d} {}{}",
                self.val,
                if self.negate { " ! > " } else { " > " }
            ),
            QueryKind::GreaterEqual => format!(
                "{d} {}{}",
                self.val,
                if self.negate { " ! >= " } else { " >= " }
            ),
            QueryKind::Range { lower, upper, .. } => {
                format!(
                    "{d}{} {lower} val {upper}",
                    if self.negate { " ! " } else { "" }
                )
            }
            QueryKind::Set(set) | QueryKind::Recursive { set, .. } => {
                let mut s = format!("{d} val{}", if self.negate { " not in " } else { " in (" });
                for v in set {
                    s += &format!("{v}, ");
                }
                s.push(')');
                s
            }
            _ => {
                if self.negate {
                    format!("not {d}")
                } else {
                    d.clone()
                }
            }
        }
    }

    /// `detail::qhelper` (used by `describeQuery`)
    pub fn describe(&self, depth: usize) -> String {
        let mut res = "  ".repeat(depth);
        res += &self.full_description();
        res.push('\n');
        for c in &self.children {
            res += &c.describe(depth + 1);
        }
        res
    }
}

fn prop_value_matches(v: &RDValue, target: &PropValue, tol: f64) -> bool {
    match target {
        PropValue::String(s) => match v {
            RDValue::String(x) => x == s,
            _ => false,
        },
        PropValue::Int(t) => v
            .as_int()
            .is_some_and(|x| query_cmp(x, *t, tol as i32) == 0),
        PropValue::UInt(t) => v
            .as_uint()
            .is_some_and(|x| query_cmp(x as i32, *t as i32, tol as i32) == 0),
        PropValue::Bool(t) => v.as_bool().is_some_and(|x| x == *t),
        PropValue::Float(t) => v
            .as_double()
            .is_some_and(|x| query_cmp_f(x, *t as f64, tol) == 0),
        PropValue::Double(t) => v.as_double().is_some_and(|x| query_cmp_f(x, *t, tol) == 0),
    }
}

// ------------------------------------------------------------------ atom data functions

/// Atom data functions (`queryAtom*` in QueryOps.h).
#[derive(Debug, Clone, PartialEq)]
pub enum AtomFn {
    Aromatic,
    Aliphatic,
    ExplicitDegree,
    TotalDegree,
    NonHydrogenDegree,
    HeavyAtomDegree,
    HCount,
    ImplicitHCount,
    HasImplicitH,
    ImplicitValence,
    ExplicitValence,
    TotalValence,
    Unsaturated,
    Num,
    Type,
    Mass,
    Isotope,
    FormalCharge,
    NegativeFormalCharge,
    Hybridization,
    NumRadicalElectrons,
    HasChiralTag,
    MissingChiralTag,
    HasHeteroatomNbrs,
    NumHeteroatomNbrs,
    HasAliphaticHeteroatomNbrs,
    NumAliphaticHeteroatomNbrs,
    InNRings,
    InRing,
    HasRingBond,
    IsBridgehead,
    MinRingSize,
    RingBondCount,
    RingMembership,
    /// `queryAtomIsInRingOfSize(at, tgt)`
    RingSize(i32),
    /// `queryAtomIsInRingOfSize(at, lower, upper, lowerOpen, upperOpen)`
    RingSizeRange {
        lower: i32,
        upper: i32,
        lower_open: bool,
        upper_open: bool,
    },
    /// `RecursiveStructureQuery::getAtIdx`
    AtIdx,
    /// `nullDataFun`
    Null,
}

/// `makeAtomType`
pub fn make_atom_type(atomic_num: i32, aromatic: bool) -> i32 {
    atomic_num + 1000 * aromatic as i32
}

/// `getAtomTypeIsAromatic`
pub fn get_atom_type_is_aromatic(val: i32) -> bool {
    val > 1000
}

/// `getAtomTypeAtomicNum`
pub fn get_atom_type_atomic_num(val: i32) -> i32 {
    if val > 1000 { val - 1000 } else { val }
}

/// `massIntegerConversionFactor`
pub const MASS_INTEGER_CONVERSION_FACTOR: i32 = 1000;

/// `queryAtomIsInRingOfSize(at, lower, upper, lowerOpen, upperOpen)`
fn atom_ring_size_in_range(
    mol: &ROMol,
    idx: u32,
    lower: i32,
    upper: i32,
    lo: bool,
    uo: bool,
) -> Result<i32> {
    for rs in mol.ring_info().atom_ring_sizes(idx)? {
        if (rs > lower || (rs == lower && !lo))
            && (upper < 0 || (rs < upper || (rs == upper && !uo)))
        {
            return Ok(rs);
        }
    }
    Ok(if lower > -1 {
        -1
    } else if upper > -1 {
        i32::MAX
    } else {
        0
    })
}

impl QueryTarget for AtomFn {
    fn eval(&self, mol: &ROMol, idx: u32) -> Result<i32> {
        let at = mol.atom(idx);
        let ri = mol.ring_info();
        let nbr_count = |f: &dyn Fn(&crate::graphmol::atom::Atom) -> bool| {
            mol.atom_neighbors(idx).filter(|&n| f(mol.atom(n))).count() as i32
        };
        let hetero = |n: &crate::graphmol::atom::Atom| n.atomic_num() != 6 && n.atomic_num() != 1;
        let aliph_hetero = |n: &crate::graphmol::atom::Atom| !n.is_aromatic() && hetero(n);
        let ring_bonds = || {
            mol.atom_bond_indices(idx)
                .filter(|&b| ri.num_bond_rings(b).unwrap_or(0) > 0)
                .count() as i32
        };
        Ok(match self {
            AtomFn::Aromatic => at.is_aromatic() as i32,
            AtomFn::Aliphatic => !at.is_aromatic() as i32,
            AtomFn::ExplicitDegree => mol.atom_degree(idx) as i32,
            AtomFn::TotalDegree => mol.total_degree(idx)? as i32,
            AtomFn::NonHydrogenDegree => nbr_count(&|n| n.atomic_num() != 1 || n.isotope() > 1),
            AtomFn::HeavyAtomDegree => nbr_count(&|n| n.atomic_num() > 1),
            AtomFn::HCount => mol.total_num_hs(idx, true)? as i32,
            AtomFn::ImplicitHCount => mol.total_num_hs(idx, false)? as i32,
            AtomFn::HasImplicitH => (mol.total_num_hs(idx, false)? > 0) as i32,
            AtomFn::ImplicitValence => at.valence(ValenceType::Implicit)? as i32,
            AtomFn::ExplicitValence => {
                at.valence(ValenceType::Explicit)? as i32 - at.num_explicit_hs() as i32
            }
            AtomFn::TotalValence => at.total_valence()? as i32,
            AtomFn::Unsaturated => (mol.total_degree(idx)? < at.total_valence()?) as i32,
            AtomFn::Num => at.atomic_num() as i32,
            AtomFn::Type => make_atom_type(at.atomic_num() as i32, at.is_aromatic()),
            AtomFn::Mass => (MASS_INTEGER_CONVERSION_FACTOR as f64 * at.mass()).round() as i32,
            AtomFn::Isotope => at.isotope() as i32,
            AtomFn::FormalCharge => at.formal_charge(),
            AtomFn::NegativeFormalCharge => -at.formal_charge(),
            AtomFn::Hybridization => at.hybridization() as i32,
            AtomFn::NumRadicalElectrons => at.num_radical_electrons() as i32,
            AtomFn::HasChiralTag => (at.chiral_tag() != ChiralType::Unspecified) as i32,
            AtomFn::MissingChiralTag => {
                (at.chiral_tag() == ChiralType::Unspecified && at.has_prop(cp::_ChiralityPossible))
                    as i32
            }
            AtomFn::HasHeteroatomNbrs => (nbr_count(&hetero) > 0) as i32,
            AtomFn::NumHeteroatomNbrs => nbr_count(&hetero),
            AtomFn::HasAliphaticHeteroatomNbrs => (nbr_count(&aliph_hetero) > 0) as i32,
            AtomFn::NumAliphaticHeteroatomNbrs => nbr_count(&aliph_hetero),
            AtomFn::InNRings | AtomFn::RingMembership => ri.num_atom_rings(idx)? as i32,
            AtomFn::InRing => (ri.num_atom_rings(idx)? != 0) as i32,
            AtomFn::HasRingBond => (ring_bonds() > 0) as i32,
            AtomFn::IsBridgehead => query_is_atom_bridgehead(mol, idx) as i32,
            AtomFn::MinRingSize => ri.min_atom_ring_size(idx)? as i32,
            AtomFn::RingBondCount => ring_bonds(),
            AtomFn::RingSize(tgt) => {
                if *tgt >= 0 && ri.is_atom_in_ring_of_size(idx, *tgt as u32)? {
                    *tgt
                } else {
                    0
                }
            }
            AtomFn::RingSizeRange {
                lower,
                upper,
                lower_open,
                upper_open,
            } => atom_ring_size_in_range(mol, idx, *lower, *upper, *lower_open, *upper_open)?,
            AtomFn::AtIdx => idx as i32,
            AtomFn::Null => 1,
        })
    }

    fn props(mol: &ROMol, idx: u32) -> &crate::rdgeneral::Dict {
        mol.atom(idx).dict()
    }
}

// ------------------------------------------------------------------ bond data functions

/// Bond data functions (`queryBond*` in QueryOps.h).
#[derive(Debug, Clone, PartialEq)]
pub enum BondFn {
    Order,
    IsSingleOrAromatic,
    IsDoubleOrAromatic,
    IsSingleOrDouble,
    IsSingleOrDoubleOrAromatic,
    Dir,
    InNRings,
    HasStereo,
    InRing,
    MinRingSize,
    /// `queryBondIsInRingOfSize<tgt>`
    RingSize(i32),
    /// `nullDataFun`
    Null,
}

impl QueryTarget for BondFn {
    fn eval(&self, mol: &ROMol, idx: u32) -> Result<i32> {
        let b = mol.bond(idx);
        let bt = b.bond_type();
        let ri = mol.ring_info();
        Ok(match self {
            BondFn::Order => bt as i32,
            BondFn::IsSingleOrAromatic => {
                matches!(bt, BondType::Single | BondType::Aromatic) as i32
            }
            BondFn::IsDoubleOrAromatic => {
                matches!(bt, BondType::Double | BondType::Aromatic) as i32
            }
            BondFn::IsSingleOrDouble => matches!(bt, BondType::Single | BondType::Double) as i32,
            BondFn::IsSingleOrDoubleOrAromatic => {
                matches!(bt, BondType::Single | BondType::Double | BondType::Aromatic) as i32
            }
            BondFn::Dir => b.bond_dir() as i32,
            BondFn::InNRings => ri.num_bond_rings(idx)? as i32,
            BondFn::HasStereo => (b.stereo() > BondStereo::StereoNone) as i32,
            BondFn::InRing => (ri.num_bond_rings(idx)? != 0) as i32,
            BondFn::MinRingSize => ri.min_bond_ring_size(idx)? as i32,
            BondFn::RingSize(tgt) => {
                if ri.is_bond_in_ring_of_size(idx, *tgt as u32)? {
                    *tgt
                } else {
                    0
                }
            }
            BondFn::Null => 1,
        })
    }

    fn props(mol: &ROMol, idx: u32) -> &crate::rdgeneral::Dict {
        mol.bond(idx).dict()
    }
}

// ------------------------------------------------------------------ constructors

/// `makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, func, descr)`
pub fn make_atom_simple_query(what: i32, func: AtomFn, descr: &str) -> AtomQuery {
    Query::with(QueryKind::Equality, what, func, descr)
}

/// `makeAtomRangeQuery`
pub fn make_atom_range_query(
    lower: i32,
    upper: i32,
    lower_open: bool,
    upper_open: bool,
    func: AtomFn,
    descr: &str,
) -> AtomQuery {
    let mut q = Query::new(QueryKind::Range {
        lower,
        upper,
        lower_open,
        upper_open,
    });
    q.data_func = Some(func);
    q.description = descr.to_string();
    q
}

macro_rules! atom_queries {
    ($($(#[$m:meta])* $name:ident($($arg:ident)?) => $val:expr, $func:expr, $descr:literal;)*) => {
        $(
            $(#[$m])*
            pub fn $name($($arg: i32)?) -> AtomQuery {
                make_atom_simple_query($val, $func, $descr)
            }
        )*
    };
}

atom_queries! {
    /// `makeAtomNumQuery`
    make_atom_num_query(what) => what, AtomFn::Num, "AtomAtomicNum";
    /// `makeAtomImplicitValenceQuery`
    make_atom_implicit_valence_query(what) => what, AtomFn::ImplicitValence, "AtomImplicitValence";
    /// `makeAtomExplicitValenceQuery`
    make_atom_explicit_valence_query(what) => what, AtomFn::ExplicitValence, "AtomExplicitValence";
    /// `makeAtomTotalValenceQuery`
    make_atom_total_valence_query(what) => what, AtomFn::TotalValence, "AtomTotalValence";
    /// `makeAtomExplicitDegreeQuery`
    make_atom_explicit_degree_query(what) => what, AtomFn::ExplicitDegree, "AtomExplicitDegree";
    /// `makeAtomTotalDegreeQuery`
    make_atom_total_degree_query(what) => what, AtomFn::TotalDegree, "AtomTotalDegree";
    /// `makeAtomHeavyAtomDegreeQuery`
    make_atom_heavy_atom_degree_query(what) => what, AtomFn::HeavyAtomDegree, "AtomHeavyAtomDegree";
    /// `makeAtomHCountQuery`
    make_atom_h_count_query(what) => what, AtomFn::HCount, "AtomHCount";
    /// `makeAtomImplicitHCountQuery`
    make_atom_implicit_h_count_query(what) => what, AtomFn::ImplicitHCount, "AtomImplicitHCount";
    /// `makeAtomHasImplicitHQuery`
    make_atom_has_implicit_h_query() => 1, AtomFn::HasImplicitH, "AtomHasImplicitH";
    /// `makeAtomAromaticQuery`
    make_atom_aromatic_query() => 1, AtomFn::Aromatic, "AtomIsAromatic";
    /// `makeAtomAliphaticQuery`
    make_atom_aliphatic_query() => 1, AtomFn::Aliphatic, "AtomIsAliphatic";
    /// `makeAtomUnsaturatedQuery`
    make_atom_unsaturated_query() => 1, AtomFn::Unsaturated, "AtomUnsaturated";
    /// `makeAtomMassQuery`
    make_atom_mass_query(what) => MASS_INTEGER_CONVERSION_FACTOR * what, AtomFn::Mass, "AtomMass";
    /// `makeAtomIsotopeQuery`
    make_atom_isotope_query(what) => what, AtomFn::Isotope, "AtomIsotope";
    /// `makeAtomFormalChargeQuery`
    make_atom_formal_charge_query(what) => what, AtomFn::FormalCharge, "AtomFormalCharge";
    /// `makeAtomNegativeFormalChargeQuery`
    make_atom_negative_formal_charge_query(what) => what, AtomFn::NegativeFormalCharge, "AtomNegativeFormalCharge";
    /// `makeAtomHybridizationQuery`
    make_atom_hybridization_query(what) => what, AtomFn::Hybridization, "AtomHybridization";
    /// `makeAtomNumRadicalElectronsQuery`
    make_atom_num_radical_electrons_query(what) => what, AtomFn::NumRadicalElectrons, "AtomNumRadicalElectrons";
    /// `makeAtomHasChiralTagQuery`
    make_atom_has_chiral_tag_query() => 1, AtomFn::HasChiralTag, "AtomHasChiralTag";
    /// `makeAtomMissingChiralTagQuery`
    make_atom_missing_chiral_tag_query() => 1, AtomFn::MissingChiralTag, "AtomMissingChiralTag";
    /// `makeAtomInRingQuery`
    make_atom_in_ring_query() => 1, AtomFn::InRing, "AtomInRing";
    /// `makeAtomInNRingsQuery`
    make_atom_in_n_rings_query(what) => what, AtomFn::InNRings, "AtomInNRings";
    /// `makeAtomMinRingSizeQuery`
    make_atom_min_ring_size_query(what) => what, AtomFn::MinRingSize, "AtomMinRingSize";
    /// `makeAtomHasRingBondQuery`
    make_atom_has_ring_bond_query() => 1, AtomFn::HasRingBond, "AtomHasRingBond";
    /// `makeAtomNumHeteroatomNbrsQuery`
    make_atom_num_heteroatom_nbrs_query(what) => what, AtomFn::NumHeteroatomNbrs, "AtomNumHeteroatomNeighbors";
    /// `makeAtomHasHeteroatomNbrsQuery`
    make_atom_has_heteroatom_nbrs_query() => 1, AtomFn::HasHeteroatomNbrs, "AtomHasHeteroatomNeighbors";
    /// `makeAtomNumAliphaticHeteroatomNbrsQuery`
    make_atom_num_aliphatic_heteroatom_nbrs_query(what) => what, AtomFn::NumAliphaticHeteroatomNbrs, "AtomNumAliphaticHeteroatomNeighbors";
    /// `makeAtomHasAliphaticHeteroatomNbrsQuery`
    make_atom_has_aliphatic_heteroatom_nbrs_query() => 1, AtomFn::HasAliphaticHeteroatomNbrs, "AtomHasAliphaticHeteroatomNeighbors";
    /// `makeAtomNonHydrogenDegreeQuery`
    make_atom_non_hydrogen_degree_query(what) => what, AtomFn::NonHydrogenDegree, "AtomNonHydrogenDegree";
    /// `makeAtomIsBridgeheadQuery`
    make_atom_is_bridgehead_query() => 1, AtomFn::IsBridgehead, "AtomIsBridgehead";
}

/// `makeAtomTypeQuery`
pub fn make_atom_type_query(num: i32, aromatic: bool) -> AtomQuery {
    make_atom_simple_query(make_atom_type(num, aromatic), AtomFn::Type, "AtomType")
}

/// `makeAtomRingBondCountQuery` (an `AtomRingQuery`)
pub fn make_atom_ring_bond_count_query(what: i32) -> AtomQuery {
    Query::with(
        QueryKind::AtomRing,
        what,
        AtomFn::RingBondCount,
        "AtomRingBondCount",
    )
}

/// `AtomRingQuery(v)` ("AtomInNRings" with `queryAtomRingMembership`)
pub fn make_atom_ring_query(what: i32) -> AtomQuery {
    Query::with(
        QueryKind::AtomRing,
        what,
        AtomFn::RingMembership,
        "AtomInNRings",
    )
}

/// `makeAtomInRingOfSizeQuery(tgt)`
pub fn make_atom_in_ring_of_size_query(tgt: i32) -> AtomQuery {
    make_atom_simple_query(tgt, AtomFn::RingSize(tgt), "AtomRingSize")
}

/// `makeAtomInRingOfSizeQuery(lower, upper, lowerOpen, upperOpen)`
pub fn make_atom_in_ring_of_size_range_query(
    lower: i32,
    upper: i32,
    lower_open: bool,
    upper_open: bool,
) -> AtomQuery {
    make_atom_range_query(
        lower,
        upper,
        lower_open,
        upper_open,
        AtomFn::RingSizeRange {
            lower,
            upper,
            lower_open,
            upper_open,
        },
        "range_AtomRingSize",
    )
}

/// `makeAtomNullQuery`
pub fn make_atom_null_query() -> AtomQuery {
    let mut q = Query::new(QueryKind::Base { always_true: true });
    q.data_func = Some(AtomFn::Null);
    q.description = "AtomNull".into();
    q
}

fn atom_or(children: Vec<AtomQuery>, negate: bool, label: &str) -> AtomQuery {
    let mut q = Query::new(QueryKind::Or);
    q.description = "AtomOr".into();
    q.negate = negate;
    q.children = children;
    q.type_label = label.into();
    q
}

/// `makeQAtomQuery`
pub fn make_q_atom_query() -> AtomQuery {
    atom_or(
        vec![make_atom_num_query(6), make_atom_num_query(1)],
        true,
        "Q",
    )
}

/// `makeQHAtomQuery`
pub fn make_qh_atom_query() -> AtomQuery {
    let mut q = make_atom_num_query(6);
    q.negate = true;
    q.type_label = "QH".into();
    q
}

/// `makeAAtomQuery`
pub fn make_a_atom_query() -> AtomQuery {
    let mut q = make_atom_num_query(1);
    q.negate = true;
    q.type_label = "A".into();
    q
}

/// `makeAHAtomQuery`
pub fn make_ah_atom_query() -> AtomQuery {
    let mut q = make_atom_null_query();
    q.type_label = "AH".into();
    q
}

/// `makeXAtomQuery`
pub fn make_x_atom_query() -> AtomQuery {
    atom_or(
        [9, 17, 35, 53, 85]
            .into_iter()
            .map(make_atom_num_query)
            .collect(),
        false,
        "X",
    )
}

/// `makeXHAtomQuery`
pub fn make_xh_atom_query() -> AtomQuery {
    let mut q = make_x_atom_query();
    q.children.push(make_atom_num_query(1));
    q.type_label = "XH".into();
    q
}

/// `makeMHAtomQuery`
pub fn make_mh_atom_query() -> AtomQuery {
    const NON_METALS: [i32; 22] = [
        0, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86,
    ];
    atom_or(
        NON_METALS.into_iter().map(make_atom_num_query).collect(),
        true,
        "MH",
    )
}

/// `makeMAtomQuery`
pub fn make_m_atom_query() -> AtomQuery {
    let mut q = make_mh_atom_query();
    q.children.push(make_atom_num_query(1));
    q.type_label = "M".into();
    q
}

/// `convertComplexNameToQuery`
pub fn complex_name_to_query(symb: &str) -> Result<AtomQuery> {
    Ok(match symb {
        "Q" => make_q_atom_query(),
        "QH" => make_qh_atom_query(),
        "A" => make_a_atom_query(),
        "AH" => make_ah_atom_query(),
        "X" => make_x_atom_query(),
        "XH" => make_xh_atom_query(),
        "M" => make_m_atom_query(),
        "MH" => make_mh_atom_query(),
        _ => return Err(RDKitError::invariant("bad complex query symbol")),
    })
}

/// `RecursiveStructureQuery(query, serialNumber)`
pub fn make_recursive_structure_query(mol: ROMol, serial_number: u32) -> AtomQuery {
    let mut q = Query::new(QueryKind::Recursive {
        mol: QueryMol(Arc::new(mol)),
        serial_number,
        set: BTreeSet::new(),
        initialized: false,
    });
    q.data_func = Some(AtomFn::AtIdx);
    q.description = "RecursiveStructure".into();
    q
}

/// `makeHasPropQuery<T>(property)`
pub fn make_has_prop_query<F: QueryTarget>(property: &str) -> Query<F> {
    let mut q = Query::new(QueryKind::HasProp(property.to_string()));
    q.description = "HasProp".into();
    q
}

/// `makePropQuery<T>(propname, val, tolerance)`
pub fn make_prop_query<F: QueryTarget>(property: &str, val: PropValue, tol: f64) -> Query<F> {
    let mut q = Query::new(QueryKind::HasPropWithValue {
        prop: property.to_string(),
        val,
        tol,
    });
    q.description = "HasPropWithValue".into();
    q
}

fn bond_equals(val: i32, func: BondFn, descr: &str, label: &str) -> BondQuery {
    let mut q = Query::with(QueryKind::Equality, val, func, descr);
    q.type_label = label.into();
    q
}

/// `makeBondOrderEqualsQuery`
pub fn make_bond_order_equals_query(what: BondType) -> BondQuery {
    bond_equals(what as i32, BondFn::Order, "BondOrder", "BondOrder")
}

/// `makeSingleOrAromaticBondQuery`
pub fn make_single_or_aromatic_bond_query() -> BondQuery {
    bond_equals(
        1,
        BondFn::IsSingleOrAromatic,
        "SingleOrAromaticBond",
        "BondOrder",
    )
}

/// `makeDoubleOrAromaticBondQuery`
pub fn make_double_or_aromatic_bond_query() -> BondQuery {
    bond_equals(
        1,
        BondFn::IsDoubleOrAromatic,
        "DoubleOrAromaticBond",
        "BondOrder",
    )
}

/// `makeSingleOrDoubleBondQuery`
pub fn make_single_or_double_bond_query() -> BondQuery {
    bond_equals(
        1,
        BondFn::IsSingleOrDouble,
        "SingleOrDoubleBond",
        "BondOrder",
    )
}

/// `makeSingleOrDoubleOrAromaticBondQuery`
pub fn make_single_or_double_or_aromatic_bond_query() -> BondQuery {
    bond_equals(
        1,
        BondFn::IsSingleOrDoubleOrAromatic,
        "SingleOrDoubleOrAromaticBond",
        "BondOrder",
    )
}

/// `makeBondDirEqualsQuery`
pub fn make_bond_dir_equals_query(what: BondDir) -> BondQuery {
    bond_equals(what as i32, BondFn::Dir, "BondDir", "")
}

/// `makeBondHasStereoQuery`
pub fn make_bond_has_stereo_query() -> BondQuery {
    bond_equals(1, BondFn::HasStereo, "BondStereo", "")
}

/// `makeBondIsInRingQuery`
pub fn make_bond_is_in_ring_query() -> BondQuery {
    bond_equals(1, BondFn::InRing, "BondInRing", "")
}

/// `makeBondInNRingsQuery`
pub fn make_bond_in_n_rings_query(what: i32) -> BondQuery {
    bond_equals(what, BondFn::InNRings, "BondInNRings", "")
}

/// `makeBondMinRingSizeQuery`
pub fn make_bond_min_ring_size_query(what: i32) -> BondQuery {
    bond_equals(what, BondFn::MinRingSize, "BondMinRingSize", "")
}

/// `makeBondInRingOfSizeQuery` (`RANGE_CHECK(3, tgt, 20)`)
pub fn make_bond_in_ring_of_size_query(tgt: i32) -> Result<BondQuery> {
    if !(3..=20).contains(&tgt) {
        return Err(RDKitError::range("tgt", tgt, 20));
    }
    Ok(bond_equals(tgt, BondFn::RingSize(tgt), "BondRingSize", ""))
}

/// `makeBondNullQuery`
pub fn make_bond_null_query() -> BondQuery {
    let mut q = Query::new(QueryKind::Base { always_true: true });
    q.data_func = Some(BondFn::Null);
    q.description = "BondNull".into();
    q
}

// ------------------------------------------------------------------ expandQuery

/// `mergeNullQueries` (NullQueryAlgebra.h). Returns the merged query.
fn merge_null_queries<F: QueryTarget>(
    mut ret: Query<F>,
    ret_null: bool,
    mut other: Query<F>,
    other_null: bool,
    how: CompositeQueryType,
) -> Query<F> {
    fn null_first<F: QueryTarget>(
        ret: Query<F>,
        other: Query<F>,
        how: CompositeQueryType,
    ) -> Query<F> {
        let negated = ret.negate;
        match how {
            CompositeQueryType::And => {
                if !negated {
                    other
                } else {
                    ret
                }
            }
            CompositeQueryType::Or => {
                if negated {
                    other
                } else {
                    ret
                }
            }
            CompositeQueryType::XOr => {
                let mut r = other;
                if !negated {
                    r.negate = !r.negate;
                }
                r
            }
            CompositeQueryType::OrSequence => ret,
        }
    }
    if ret_null && other_null {
        let (nq, no) = (ret.negate, other.negate);
        match how {
            CompositeQueryType::And => {
                if !nq && no {
                    ret.negate = true;
                }
            }
            CompositeQueryType::Or => {
                if nq && !no {
                    ret.negate = false;
                }
            }
            CompositeQueryType::XOr => {
                if !nq && !no {
                    ret.negate = true;
                } else if (nq as u8 + no as u8) == 1 {
                    ret.negate = false;
                }
            }
            CompositeQueryType::OrSequence => {}
        }
        ret
    } else if ret_null {
        null_first(ret, other, how)
    } else {
        std::mem::swap(&mut ret, &mut other);
        null_first(ret, other, how)
    }
}

/// `QueryAtom::expandQuery` / `QueryBond::expandQuery`: combines `current`
/// with `what`. `prefix` is "Atom" or "Bond".
pub fn expand_query<F: QueryTarget>(
    current: Query<F>,
    what: Query<F>,
    how: CompositeQueryType,
    maintain_order: bool,
    prefix: &str,
) -> Result<Query<F>> {
    let null_descr = format!("{prefix}Null");
    let this_null = current.description == null_descr;
    let other_null = what.description == null_descr;
    if this_null || other_null {
        return Ok(merge_null_queries(
            current, this_null, what, other_null, how,
        ));
    }
    let (kind, descr) = match how {
        CompositeQueryType::And => (QueryKind::And, "And"),
        CompositeQueryType::Or => (QueryKind::Or, "Or"),
        CompositeQueryType::XOr => (QueryKind::XOr, "Xor"),
        CompositeQueryType::OrSequence => {
            return Err(RDKitError::invariant("unrecognized combination query"));
        }
    };
    let mut q = Query::new(kind);
    q.description = format!("{prefix}{descr}");
    if maintain_order {
        q.children = vec![current, what];
    } else {
        q.children = vec![what, current];
    }
    Ok(q)
}

// ------------------------------------------------------------------ query inspection

/// `isComplexQuery(const Bond*)`
pub fn is_complex_bond_query(q: Option<&BondQuery>) -> bool {
    let Some(q) = q else {
        return false;
    };
    if q.negate {
        return true;
    }
    match q.description.as_str() {
        "BondOrder" | "SingleOrAromaticBond" => false,
        "BondAnd" | "BondXor" => true,
        "BondOr" if q.children.len() == 2 => {
            for c in &q.children {
                if c.description != "BondOrder" || c.negate {
                    return true;
                }
                if c.val != BondType::Single as i32 && c.val != BondType::Aromatic as i32 {
                    return true;
                }
            }
            false
        }
        _ => true,
    }
}

fn complex_query_helper(q: &AtomQuery, has_at_num: &mut bool) -> bool {
    if q.negate {
        return true;
    }
    match q.description.as_str() {
        "AtomAtomicNum" | "AtomType" => {
            *has_at_num = true;
            false
        }
        "AtomOr" | "AtomXor" => true,
        "AtomAnd" => q
            .children
            .iter()
            .any(|c| complex_query_helper(c, has_at_num)),
        _ => false,
    }
}

/// `isComplexQuery(const Atom*)`
pub fn is_complex_atom_query(q: Option<&AtomQuery>) -> bool {
    let Some(q) = q else {
        return false;
    };
    if q.negate {
        return true;
    }
    match q.description.as_str() {
        "AtomNull" | "AtomAtomicNum" | "AtomType" => false,
        "AtomOr" | "AtomXor" => true,
        "AtomAnd" => {
            let mut has_at_num = false;
            if complex_query_helper(q, &mut has_at_num) {
                return true;
            }
            !has_at_num
        }
        _ => true,
    }
}

fn atom_list_helper(q: &AtomQuery, ignore_negation: bool) -> bool {
    if !ignore_negation && q.negate {
        return false;
    }
    match q.description.as_str() {
        "AtomAtomicNum" | "AtomType" => true,
        "AtomOr" => q
            .children
            .iter()
            .all(|c| atom_list_helper(c, ignore_negation)),
        _ => false,
    }
}

/// `isAtomListQuery`
pub fn is_atom_list_query(q: Option<&AtomQuery>, atomic_num: u32) -> bool {
    let Some(q) = q else {
        return false;
    };
    if q.description == "AtomOr" {
        q.children.iter().all(|c| atom_list_helper(c, false))
    } else if q.negate && atom_list_helper(q, true) {
        true
    } else {
        q.description == "AtomAtomicNum" && q.val != atomic_num as i32
    }
}

/// `getAtomListQueryVals`
pub fn get_atom_list_query_vals(q: &AtomQuery, vals: &mut Vec<i32>) -> Result<()> {
    let fix = |v: i32| if v >= 1000 { v - 1000 } else { v };
    match q.description.as_str() {
        "AtomOr" => {
            for c in &q.children {
                let d = c.description.as_str();
                if c.negate || !matches!(d, "AtomOr" | "AtomAtomicNum" | "AtomType") {
                    return Err(RDKitError::Value("bad query type1".into()));
                }
                match d {
                    "AtomOr" => get_atom_list_query_vals(c, vals)?,
                    "AtomAtomicNum" => vals.push(c.val),
                    _ => vals.push(fix(c.val)),
                }
            }
        }
        "AtomAtomicNum" => vals.push(q.val),
        "AtomType" => vals.push(fix(q.val)),
        _ => return Err(RDKitError::invariant("bad query type")),
    }
    Ok(())
}

/// `isAtomAromatic(const Atom*)` for an atom of `mol`
pub fn is_atom_aromatic(mol: &ROMol, idx: u32) -> bool {
    let a = mol.atom(idx);
    let Some(q) = a.query() else {
        return crate::graphmol::atom::is_aromatic_atom(mol, idx);
    };
    match q.description.as_str() {
        "AtomAtomicNum" => a.is_aromatic(),
        "AtomIsAromatic" => !q.negate,
        "AtomIsAliphatic" => q.negate,
        "AtomType" => get_atom_type_is_aromatic(q.val) != q.negate,
        "AtomAnd" => {
            let c = &q.children;
            if !c.is_empty() && c[0].description == "AtomAtomicNum" {
                if q.negate {
                    false
                } else if c.len() > 1 && c[1].description == "AtomIsAliphatic" {
                    false
                } else {
                    c.len() > 1 && c[1].description == "AtomIsAromatic"
                }
            } else {
                false
            }
        }
        _ => false,
    }
}

/// `isAtomDummy`
pub fn is_atom_dummy(mol: &ROMol, idx: u32) -> bool {
    let a = mol.atom(idx);
    match a.query() {
        None => a.atomic_num() == 0,
        Some(q) => !q.negate && q.description == "AtomNull",
    }
}

const BOND_ORDER_QUERY_FUNCTIONS: [&str; 5] = [
    "BondOrder",
    "SingleOrAromaticBond",
    "DoubleOrAromaticBond",
    "SingleOrDoubleBond",
    "SingleOrDoubleOrAromaticBond",
];

/// `QueryOps::hasBondTypeQuery`
pub fn has_bond_type_query(q: &BondQuery) -> bool {
    if q.type_label == "BondOrder"
        || (q.type_label.is_empty() && BOND_ORDER_QUERY_FUNCTIONS.contains(&q.description.as_str()))
    {
        return true;
    }
    q.children.iter().any(has_bond_type_query)
}

fn complex_bond_type_helper(q: &BondQuery, mut seen_bond_order: bool) -> bool {
    let is_bond_order = q.description == "BondOrder";
    if BOND_ORDER_QUERY_FUNCTIONS.contains(&q.description.as_str())
        && (seen_bond_order || !is_bond_order || q.negate)
    {
        return true;
    }
    for c in &q.children {
        if complex_bond_type_helper(c, seen_bond_order | is_bond_order) {
            return true;
        }
        if c.description == "BondOrder" {
            seen_bond_order = true;
        }
    }
    false
}

/// `QueryOps::hasComplexBondTypeQuery`
pub fn has_complex_bond_type_query(q: &BondQuery) -> bool {
    complex_bond_type_helper(q, false)
}

/// `detail::hasRecursiveQuery`
fn has_recursive(q: &AtomQuery, check_initialized: bool) -> bool {
    if let QueryKind::Recursive { initialized, .. } = &q.kind {
        if !check_initialized || !initialized {
            return true;
        }
    }
    q.children
        .iter()
        .any(|c| has_recursive(c, check_initialized))
}

/// `hasRecursiveQuery(atom)`
pub fn has_recursive_query(q: Option<&AtomQuery>) -> bool {
    q.is_some_and(|q| has_recursive(q, false))
}

/// `hasUninitializedRecursiveQuery(atom)`
pub fn has_uninitialized_recursive_query(q: Option<&AtomQuery>) -> bool {
    q.is_some_and(|q| has_recursive(q, true))
}

/// `QueryOps::completeMolQueries`
pub fn complete_mol_queries(mol: &mut ROMol, magic_val: u32) -> Result<()> {
    fn complete(q: &mut AtomQuery, mol: &ROMol, idx: u32, magic: u32) -> Result<()> {
        let is_eq = matches!(
            q.kind,
            QueryKind::Equality
                | QueryKind::Less
                | QueryKind::LessEqual
                | QueryKind::Greater
                | QueryKind::GreaterEqual
                | QueryKind::AtomRing
                | QueryKind::HasProp(_)
                | QueryKind::HasPropWithValue { .. }
        );
        if is_eq && q.val as u32 == magic {
            q.val = q.data(mol, idx)?;
        }
        for c in &mut q.children {
            complete(c, mol, idx, magic)?;
        }
        Ok(())
    }
    for i in 0..mol.num_atoms() {
        if let Some(mut q) = mol.atom_mut(i).query.take() {
            let r = complete(&mut q, mol, i, magic_val);
            mol.atom_mut(i).query = Some(q);
            r?;
        }
    }
    Ok(())
}

/// The `queriesMatch` equality descriptions (QueryAtom.cpp).
const EQUALITY_QUERIES: [&str; 20] = [
    "AtomType",
    "AtomRingBondCount",
    "AtomRingSize",
    "AtomMinRingSize",
    "AtomImplicitValence",
    "AtomExplicitValence",
    "AtomTotalValence",
    "AtomAtomicNum",
    "AtomExplicitDegree",
    "AtomTotalDegree",
    "AtomHCount",
    "AtomIsAromatic",
    "AtomIsAliphatic",
    "AtomUnsaturated",
    "AtomMass",
    "AtomFormalCharge",
    "AtomNegativeFormalCharge",
    "AtomHybridization",
    "AtomInRing",
    "AtomInNRings",
];

/// `queriesMatch` (QueryAtom.cpp): whether two atom queries can match each other.
pub fn atom_queries_match(q1: &AtomQuery, q2: &AtomQuery) -> bool {
    let local = |a: &AtomQuery, b: &AtomQuery| {
        if a.negate == b.negate {
            a.val == b.val
        } else {
            a.val != b.val
        }
    };
    let (d1, d2) = (q1.description.as_str(), q2.description.as_str());
    if d1 == "AtomNull" || d2 == "AtomNull" {
        true
    } else if d1 == "AtomOr" {
        q1.children.iter().any(|c1| {
            if d2 == "AtomOr" {
                q2.children.iter().any(|c2| atom_queries_match(c1, c2))
            } else {
                atom_queries_match(c1, q2)
            }
        })
    } else if d1 == "AtomAnd" {
        q1.children.iter().all(|c1| {
            if d2 == "AtomAnd" {
                q2.children.iter().any(|c2| atom_queries_match(c1, c2))
            } else {
                atom_queries_match(c1, q2)
            }
        })
    } else if d2 == "AtomOr" {
        q2.children.iter().any(|c2| atom_queries_match(q1, c2))
    } else if d2 == "AtomAnd" {
        q2.children.iter().all(|c2| atom_queries_match(q1, c2))
    } else if d1 == d2 && EQUALITY_QUERIES.contains(&d1) {
        local(q1, q2)
    } else {
        false
    }
}
