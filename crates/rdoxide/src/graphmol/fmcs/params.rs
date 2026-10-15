//! Parameters, comparators and results of the FMCS algorithm (upstream
//! `FMCS/FMCS.h` and the comparator part of `FMCS/FMCS.cpp`).

use std::sync::Arc;

use super::graph::Graph;
use crate::Result;
use crate::graphmol::{BondStereo, BondType, ROMol, atom::ChiralType};

/// `AtomComparator`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomComparator {
    Any,
    Elements,
    Isotopes,
    AnyHeavyAtom,
}

/// `BondComparator`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondComparator {
    Any,
    Order,
    OrderExact,
}

/// `RingComparator`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingComparator {
    IgnoreRingFusion,
    PermissiveRingFusion,
    StrictRingFusion,
}

/// `MCSAtomCompareParameters`
#[derive(Debug, Clone, PartialEq)]
pub struct MCSAtomCompareParameters {
    pub match_valences: bool,
    pub match_chiral_tag: bool,
    pub match_formal_charge: bool,
    pub ring_matches_ring_only: bool,
    pub complete_rings_only: bool,
    pub match_isotope: bool,
    pub max_distance: f64,
}

impl Default for MCSAtomCompareParameters {
    fn default() -> Self {
        MCSAtomCompareParameters {
            match_valences: false,
            match_chiral_tag: false,
            match_formal_charge: false,
            ring_matches_ring_only: false,
            complete_rings_only: false,
            match_isotope: false,
            max_distance: -1.0,
        }
    }
}

/// `MCSBondCompareParameters`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MCSBondCompareParameters {
    pub ring_matches_ring_only: bool,
    pub complete_rings_only: bool,
    pub match_fused_rings: bool,
    pub match_fused_rings_strict: bool,
    pub match_stereo: bool,
}

/// A user-defined atom comparison (`MCSAtomCompareFunction`); the user data
/// of upstream is captured by the closure.
pub type AtomCompareFn =
    Arc<dyn Fn(&MCSAtomCompareParameters, &ROMol, u32, &ROMol, u32) -> Result<bool> + Send + Sync>;
/// A user-defined bond comparison (`MCSBondCompareFunction`).
pub type BondCompareFn =
    Arc<dyn Fn(&MCSBondCompareParameters, &ROMol, u32, &ROMol, u32) -> Result<bool> + Send + Sync>;
/// `MCSAcceptanceFunction`: query, target, (query, target) atom pairs and
/// bond pairs of a candidate MCS.
pub type AcceptanceFn = Arc<
    dyn Fn(&ROMol, &ROMol, &[(i32, i32)], &[(i32, i32)], &MCSParameters) -> Result<bool>
        + Send
        + Sync,
>;
/// `MCSFinalMatchCheckFunction`: the match as parallel arrays of query and
/// target graph vertices, then the query molecule and graph and the target
/// molecule and graph.
pub type FinalMatchCheckFn = Arc<
    dyn Fn(&[u32], &[u32], &ROMol, &Graph, &ROMol, &Graph, &MCSParameters) -> Result<bool>
        + Send
        + Sync,
>;
/// `MCSProgressCallback`: return false to interrupt the search.
pub type ProgressFn = Arc<dyn Fn(&MCSProgressData, &MCSParameters) -> Result<bool> + Send + Sync>;

/// The atom comparison of a search; the predefined ones are recognized by
/// the algorithm (they allow cheaper atom labels).
#[derive(Clone)]
pub enum AtomTyper {
    Any,
    Elements,
    Isotopes,
    AnyHeavyAtom,
    Custom(AtomCompareFn),
}

impl std::fmt::Debug for AtomTyper {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AtomTyper::Any => f.write_str("Any"),
            AtomTyper::Elements => f.write_str("Elements"),
            AtomTyper::Isotopes => f.write_str("Isotopes"),
            AtomTyper::AnyHeavyAtom => f.write_str("AnyHeavyAtom"),
            AtomTyper::Custom(_) => f.write_str("Custom"),
        }
    }
}

impl AtomTyper {
    pub fn compare(
        &self,
        p: &MCSAtomCompareParameters,
        mol1: &ROMol,
        atom1: u32,
        mol2: &ROMol,
        atom2: u32,
    ) -> Result<bool> {
        match self {
            AtomTyper::Any => mcs_atom_compare_any(p, mol1, atom1, mol2, atom2),
            AtomTyper::Elements => mcs_atom_compare_elements(p, mol1, atom1, mol2, atom2),
            AtomTyper::Isotopes => mcs_atom_compare_isotopes(p, mol1, atom1, mol2, atom2),
            AtomTyper::AnyHeavyAtom => mcs_atom_compare_any_heavy_atom(p, mol1, atom1, mol2, atom2),
            AtomTyper::Custom(f) => f(p, mol1, atom1, mol2, atom2),
        }
    }
}

impl From<AtomComparator> for AtomTyper {
    fn from(c: AtomComparator) -> Self {
        match c {
            AtomComparator::Any => AtomTyper::Any,
            AtomComparator::Elements => AtomTyper::Elements,
            AtomComparator::Isotopes => AtomTyper::Isotopes,
            AtomComparator::AnyHeavyAtom => AtomTyper::AnyHeavyAtom,
        }
    }
}

/// The bond comparison of a search.
#[derive(Clone)]
pub enum BondTyper {
    Any,
    Order,
    OrderExact,
    Custom(BondCompareFn),
}

impl std::fmt::Debug for BondTyper {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BondTyper::Any => f.write_str("Any"),
            BondTyper::Order => f.write_str("Order"),
            BondTyper::OrderExact => f.write_str("OrderExact"),
            BondTyper::Custom(_) => f.write_str("Custom"),
        }
    }
}

impl BondTyper {
    pub fn compare(
        &self,
        p: &MCSBondCompareParameters,
        mol1: &ROMol,
        bond1: u32,
        mol2: &ROMol,
        bond2: u32,
    ) -> Result<bool> {
        match self {
            BondTyper::Any => mcs_bond_compare_any(p, mol1, bond1, mol2, bond2),
            BondTyper::Order => mcs_bond_compare_order(p, mol1, bond1, mol2, bond2, true),
            BondTyper::OrderExact => mcs_bond_compare_order(p, mol1, bond1, mol2, bond2, false),
            BondTyper::Custom(f) => f(p, mol1, bond1, mol2, bond2),
        }
    }
}

impl From<BondComparator> for BondTyper {
    fn from(c: BondComparator) -> Self {
        match c {
            BondComparator::Any => BondTyper::Any,
            BondComparator::Order => BondTyper::Order,
            BondComparator::OrderExact => BondTyper::OrderExact,
        }
    }
}

/// `MCSProgressData`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MCSProgressData {
    pub num_atoms: u32,
    pub num_bonds: u32,
    pub seed_processed: u32,
}

/// `MCSParameters`
#[derive(Clone)]
pub struct MCSParameters {
    pub store_all: bool,
    pub maximize_bonds: bool,
    /// fraction of the molecules that must contain the MCS (1.0: all)
    pub threshold: f64,
    /// in seconds, 0: none
    pub timeout: u32,
    pub verbose: bool,
    pub atom_compare_parameters: MCSAtomCompareParameters,
    pub bond_compare_parameters: MCSBondCompareParameters,
    pub atom_typer: AtomTyper,
    pub bond_typer: BondTyper,
    /// Replaces the default (timeout) check when set; return false to stop.
    pub progress_callback: Option<ProgressFn>,
    /// Accepts or rejects a growing MCS candidate.
    pub final_match_checker: Option<FinalMatchCheckFn>,
    /// Accepts or rejects a fully grown MCS candidate.
    pub should_accept_mcs: Option<AcceptanceFn>,
    /// SMARTS of a user-defined seed, or empty
    pub initial_seed: String,
}

impl Default for MCSParameters {
    fn default() -> Self {
        MCSParameters {
            store_all: false,
            maximize_bonds: true,
            threshold: 1.0,
            timeout: 0,
            verbose: false,
            atom_compare_parameters: MCSAtomCompareParameters::default(),
            bond_compare_parameters: MCSBondCompareParameters::default(),
            atom_typer: AtomTyper::Elements,
            bond_typer: BondTyper::Order,
            progress_callback: None,
            final_match_checker: None,
            should_accept_mcs: None,
            initial_seed: String::new(),
        }
    }
}

impl std::fmt::Debug for MCSParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MCSParameters")
            .field("store_all", &self.store_all)
            .field("maximize_bonds", &self.maximize_bonds)
            .field("threshold", &self.threshold)
            .field("timeout", &self.timeout)
            .field("atom_typer", &self.atom_typer)
            .field("bond_typer", &self.bond_typer)
            .finish_non_exhaustive()
    }
}

impl MCSParameters {
    /// `setMCSAtomTyperFromEnum`
    pub fn set_atom_typer_from_enum(&mut self, c: AtomComparator) {
        self.atom_typer = c.into();
    }

    /// `setMCSAtomTyperFromConstChar`: unknown names are ignored.
    pub fn set_atom_typer_from_str(&mut self, name: &str) {
        let c = match name {
            "Any" => AtomComparator::Any,
            "Elements" => AtomComparator::Elements,
            "Isotopes" => AtomComparator::Isotopes,
            "AnyHeavy" => AtomComparator::AnyHeavyAtom,
            _ => return,
        };
        self.set_atom_typer_from_enum(c);
    }

    /// `setMCSBondTyperFromEnum`
    pub fn set_bond_typer_from_enum(&mut self, c: BondComparator) {
        self.bond_typer = c.into();
    }

    /// `setMCSBondTyperFromConstChar`: unknown names are ignored.
    pub fn set_bond_typer_from_str(&mut self, name: &str) {
        let c = match name {
            "Any" => BondComparator::Any,
            "Order" => BondComparator::Order,
            "OrderExact" => BondComparator::OrderExact,
            _ => return,
        };
        self.set_bond_typer_from_enum(c);
    }
}

/// `MCSResult`
#[derive(Debug, Clone, Default)]
pub struct MCSResult {
    pub num_atoms: u32,
    pub num_bonds: u32,
    pub smarts_string: String,
    /// interrupted by the timeout or the progress callback; the MCS found so
    /// far is still valid
    pub canceled: bool,
    pub query_mol: Option<Arc<ROMol>>,
    /// all degenerate MCSs (SMARTS -> query) when `store_all` is set
    pub degenerate_smarts_query_mol_dict: std::collections::BTreeMap<String, Arc<ROMol>>,
}

impl MCSResult {
    pub fn is_completed(&self) -> bool {
        !self.canceled
    }
}

// ------------------------------------------------------------------ atom comparison

/// `checkAtomRingMatch`
pub fn check_atom_ring_match(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
) -> bool {
    if p.ring_matches_ring_only {
        let in_ring1 = mol1.ring_info().num_atom_rings(atom1).unwrap_or(0) > 0;
        let in_ring2 = mol2.ring_info().num_atom_rings(atom2).unwrap_or(0) > 0;
        in_ring1 == in_ring2
    } else {
        true
    }
}

/// `checkAtomCharge`
pub fn check_atom_charge(mol1: &ROMol, atom1: u32, mol2: &ROMol, atom2: u32) -> bool {
    mol1.atom(atom1).formal_charge() == mol2.atom(atom2).formal_charge()
}

fn is_tetrahedral(tag: ChiralType) -> bool {
    matches!(tag, ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW)
}

/// `checkAtomChirality`
pub fn check_atom_chirality(mol1: &ROMol, atom1: u32, mol2: &ROMol, atom2: u32) -> bool {
    if is_tetrahedral(mol1.atom(atom1).chiral_tag()) {
        return is_tetrahedral(mol2.atom(atom2).chiral_tag());
    }
    true
}

/// `checkAtomDistance`: the molecules need a conformer.
pub fn check_atom_distance(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
) -> Result<bool> {
    let pos1 = mol1.conformer(-1)?.atom_pos(atom1)?;
    let pos2 = mol2.conformer(-1)?.atom_pos(atom2)?;
    Ok((pos1 - pos2).length() <= p.max_distance)
}

/// The checks shared by all atom comparators once the element/isotope
/// criterion has passed.
fn atom_common_checks(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
) -> Result<bool> {
    if p.match_chiral_tag && !check_atom_chirality(mol1, atom1, mol2, atom2) {
        return Ok(false);
    }
    if p.match_formal_charge && !check_atom_charge(mol1, atom1, mol2, atom2) {
        return Ok(false);
    }
    if p.max_distance > 0.0 && !check_atom_distance(p, mol1, atom1, mol2, atom2)? {
        return Ok(false);
    }
    if p.ring_matches_ring_only {
        return Ok(check_atom_ring_match(p, mol1, atom1, mol2, atom2));
    }
    Ok(true)
}

/// `MCSAtomCompareAny`
pub fn mcs_atom_compare_any(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
) -> Result<bool> {
    atom_common_checks(p, mol1, atom1, mol2, atom2)
}

/// `MCSAtomCompareElements`
pub fn mcs_atom_compare_elements(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
) -> Result<bool> {
    let (a1, a2) = (mol1.atom(atom1), mol2.atom(atom2));
    if a1.atomic_num() != a2.atomic_num() {
        return Ok(false);
    }
    if p.match_valences && a1.total_valence()? != a2.total_valence()? {
        return Ok(false);
    }
    atom_common_checks(p, mol1, atom1, mol2, atom2)
}

/// `MCSAtomCompareIsotopes`
pub fn mcs_atom_compare_isotopes(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
) -> Result<bool> {
    if mol1.atom(atom1).isotope() != mol2.atom(atom2).isotope() {
        return Ok(false);
    }
    atom_common_checks(p, mol1, atom1, mol2, atom2)
}

/// `MCSAtomCompareAnyHeavyAtom`
pub fn mcs_atom_compare_any_heavy_atom(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
) -> Result<bool> {
    let (n1, n2) = (mol1.atom(atom1).atomic_num(), mol2.atom(atom2).atomic_num());
    if n1 == n2 || (n1 > 1 && n2 > 1) {
        return mcs_atom_compare_any(p, mol1, atom1, mol2, atom2);
    }
    Ok(false)
}

// ------------------------------------------------------------------ bond comparison

/// `BondMatchOrderMatrix::isEqual`
fn bond_types_equal(t1: BondType, t2: BondType, ignore_aromatization: bool) -> bool {
    use BondType::*;
    if t1 == t2 || t1 == Unspecified || t2 == Unspecified || t1 == Zero || t2 == Zero {
        return true;
    }
    if !ignore_aromatization {
        return false;
    }
    let pairs = [
        (Single, Aromatic),
        (Single, OneAndAHalf),
        (Double, TwoAndAHalf),
        (Triple, ThreeAndAHalf),
        (Quadruple, FourAndAHalf),
        (Quintuple, FiveAndAHalf),
    ];
    pairs
        .iter()
        .any(|&(a, b)| (t1 == a && t2 == b) || (t1 == b && t2 == a))
}

/// `checkBondStereo`
pub fn check_bond_stereo(mol1: &ROMol, bond1: u32, mol2: &ROMol, bond2: u32) -> bool {
    let (b1, b2) = (mol1.bond(bond1), mol2.bond(bond2));
    let (bs1, bs2) = (b1.stereo(), b2.stereo());
    if b1.bond_type() == BondType::Double && b2.bond_type() == BondType::Double {
        if bs1 > BondStereo::StereoAny && bs2 <= BondStereo::StereoAny {
            return false;
        }
        return bs1 == bs2;
    }
    true
}

/// `havePairOfCompatibleRings`
pub fn have_pair_of_compatible_rings(mol1: &ROMol, bond1: u32, mol2: &ROMol, bond2: u32) -> bool {
    let (ri1, ri2) = (mol1.ring_info(), mol2.ring_info());
    let (bond_rings1, bond_rings2) = (ri1.bond_rings(), ri2.bond_rings());
    let members = |ri: &crate::graphmol::ring_info::RingInfo, b| {
        ri.bond_members(b).map(<[i32]>::to_vec).unwrap_or_default()
    };
    for ring_idx1 in members(ri1, bond1) {
        let ring1 = &bond_rings1[ring_idx1 as usize];
        let fused1 = ring_is_fused(ri1, ring_idx1 as u32);
        for ring_idx2 in members(ri2, bond2) {
            let ring2 = &bond_rings2[ring_idx2 as usize];
            if ring1.len() == ring2.len() {
                return true;
            }
            if fused1 && ring2.len() > ring1.len() {
                return true;
            }
            if ring_is_fused(ri2, ring_idx2 as u32) && ring1.len() > ring2.len() {
                return true;
            }
        }
    }
    false
}

/// `RingInfo::isRingFused` without the lazily built cache: a ring is fused
/// when one of its bonds belongs to another ring too.
pub(crate) fn ring_is_fused(ri: &crate::graphmol::ring_info::RingInfo, ring: u32) -> bool {
    ri.bond_rings()[ring as usize]
        .iter()
        .any(|&b| ri.num_bond_rings(b as u32).unwrap_or(0) > 1)
}

/// `checkBondRingMatch`
pub fn check_bond_ring_match(
    p: &MCSBondCompareParameters,
    mol1: &ROMol,
    bond1: u32,
    mol2: &ROMol,
    bond2: u32,
) -> bool {
    let in_ring1 = mol1.ring_info().num_bond_rings(bond1).unwrap_or(0) > 0;
    let in_ring2 = mol2.ring_info().num_bond_rings(bond2).unwrap_or(0) > 0;
    let mut res = in_ring1 == in_ring2;
    if p.complete_rings_only && in_ring1 && in_ring2 {
        res = have_pair_of_compatible_rings(mol1, bond1, mol2, bond2);
    }
    res
}

/// `MCSBondCompareAny`
pub fn mcs_bond_compare_any(
    p: &MCSBondCompareParameters,
    mol1: &ROMol,
    bond1: u32,
    mol2: &ROMol,
    bond2: u32,
) -> Result<bool> {
    if p.match_stereo && !check_bond_stereo(mol1, bond1, mol2, bond2) {
        return Ok(false);
    }
    if p.ring_matches_ring_only {
        return Ok(check_bond_ring_match(p, mol1, bond1, mol2, bond2));
    }
    Ok(true)
}

/// `MCSBondCompareOrder` (`ignore_aromatization`) and
/// `MCSBondCompareOrderExact`
pub fn mcs_bond_compare_order(
    p: &MCSBondCompareParameters,
    mol1: &ROMol,
    bond1: u32,
    mol2: &ROMol,
    bond2: u32,
    ignore_aromatization: bool,
) -> Result<bool> {
    let (t1, t2) = (mol1.bond(bond1).bond_type(), mol2.bond(bond2).bond_type());
    if !bond_types_equal(t1, t2, ignore_aromatization) {
        return Ok(false);
    }
    if p.match_stereo && !check_bond_stereo(mol1, bond1, mol2, bond2) {
        return Ok(false);
    }
    if p.ring_matches_ring_only {
        return Ok(check_bond_ring_match(p, mol1, bond1, mol2, bond2));
    }
    Ok(true)
}

/// Parses `MCSParameters` from JSON (`parseMCSParametersJSON`).
pub fn parse_mcs_parameters_json(json: &str, params: &mut MCSParameters) -> Result<()> {
    use crate::rdgeneral::ptree::PTree;
    if json.is_empty() {
        return Ok(());
    }
    let pt = PTree::from_json(json)?;
    let boolean = |key: &str, default: bool| pt.get_bool(key, default);
    let double = |key: &str, default: f64| {
        pt.get_str(key)
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(default)
    };
    let p = params;
    p.maximize_bonds = boolean("MaximizeBonds", p.maximize_bonds);
    p.threshold = double("Threshold", p.threshold);
    p.timeout = pt.get_u32("Timeout", p.timeout);
    let (a, b) = (
        &mut p.atom_compare_parameters,
        &mut p.bond_compare_parameters,
    );
    a.match_valences = boolean("MatchValences", a.match_valences);
    a.match_chiral_tag = boolean("MatchChiralTag", a.match_chiral_tag);
    a.match_formal_charge = boolean("MatchFormalCharge", a.match_formal_charge);
    a.ring_matches_ring_only = boolean("RingMatchesRingOnly", a.ring_matches_ring_only);
    a.max_distance = double("MaxDistance", a.max_distance);
    b.ring_matches_ring_only = boolean("RingMatchesRingOnly", b.ring_matches_ring_only);
    a.ring_matches_ring_only = boolean("AtomRingMatchesRingOnly", a.ring_matches_ring_only);
    b.ring_matches_ring_only = boolean("BondRingMatchesRingOnly", b.ring_matches_ring_only);
    b.complete_rings_only = boolean("CompleteRingsOnly", b.complete_rings_only);
    a.complete_rings_only = boolean("AtomCompleteRingsOnly", a.complete_rings_only);
    b.complete_rings_only = boolean("BondCompleteRingsOnly", b.complete_rings_only);
    b.match_fused_rings = boolean("MatchFusedRings", b.match_fused_rings);
    b.match_fused_rings_strict = boolean("MatchFusedRingsStrict", b.match_fused_rings_strict);
    b.match_stereo = boolean("MatchStereo", b.match_stereo);
    p.store_all = boolean("StoreAll", p.store_all);
    p.set_atom_typer_from_str(pt.get_str("AtomCompare").unwrap_or("def"));
    p.set_bond_typer_from_str(pt.get_str("BondCompare").unwrap_or("def"));
    p.initial_seed = pt.get_str("InitialSeed").unwrap_or("").to_string();
    Ok(())
}
