//! `MolOps::adjustQueryProperties`: turns a molecule into a more useful
//! substructure query (degree and ring queries, generic atoms and bonds,
//! MDL five-ring aromaticity, ...).
//!
//! Ported from RDKit `Code/GraphMol/AdjustQuery.cpp` (2026.09.1), including
//! the JSON parameter format (`parseAdjustQueryParametersFromJSON`).
//!
//! Copyright (C) 2015-2020 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondStereo, BondType};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::mol_ops::sanitize_flags::{SANITIZE_SETAROMATICITY, SANITIZE_SYMMRINGS};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::queries::{
    AtomQuery, BondQuery, CompositeQueryType, Query, QueryKind, make_atom_explicit_degree_query,
    make_atom_heavy_atom_degree_query, make_atom_in_n_rings_query, make_atom_in_ring_query,
    make_atom_null_query, make_bond_null_query, make_bond_order_equals_query,
    make_single_or_aromatic_bond_query, make_single_or_double_or_aromatic_bond_query,
};
use crate::graphmol::query_ops::{atom_expand_query, query_atom_from_atom};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::ptree::PTree;
use crate::{RDKitError, Result, rd_log_warning};

/// `MolOps::AdjustQueryWhichFlags`: which atoms/bonds an adjustment skips.
pub mod adjust_flags {
    pub const ADJUST_IGNORENONE: u32 = 0x0;
    pub const ADJUST_IGNORECHAINS: u32 = 0x1;
    pub const ADJUST_IGNOREDUMMIES: u32 = 0x2;
    pub const ADJUST_IGNORERINGS: u32 = 0x4;
    pub const ADJUST_IGNORENONDUMMIES: u32 = 0x8;
    pub const ADJUST_IGNOREMAPPED: u32 = 0x10;
    pub const ADJUST_IGNOREALL: u32 = 0xFFF_FFFF;
}
use adjust_flags::*;

const CONJUGATED_OR_AROMATIC: &str = "_conjugatedOrAromatic";

/// `MolOps::AdjustQueryParameters`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdjustQueryParameters {
    /// add degree queries
    pub adjust_degree: bool,
    pub adjust_degree_flags: u32,
    /// add ring-count queries
    pub adjust_ring_count: bool,
    pub adjust_ring_count_flags: u32,
    /// convert dummy atoms without isotope labels to any-atom queries
    pub make_dummies_queries: bool,
    /// perceive and set aromaticity
    pub aromatize_if_possible: bool,
    /// convert bonds to generic queries (any bonds)
    pub make_bonds_generic: bool,
    pub make_bonds_generic_flags: u32,
    /// convert atoms to generic queries (any atoms)
    pub make_atoms_generic: bool,
    pub make_atoms_generic_flags: u32,
    /// adjust the heavy-atom degree instead of the overall degree
    pub adjust_heavy_degree: bool,
    pub adjust_heavy_degree_flags: u32,
    /// add ring-chain queries
    pub adjust_ring_chain: bool,
    pub adjust_ring_chain_flags: u32,
    /// remove stereochemistry from double bonds that do not have `molStereoCare`
    pub use_stereo_care_for_bonds: bool,
    /// bond queries in conjugated five-rings become SINGLE|DOUBLE|AROMATIC
    pub adjust_conjugated_five_rings: bool,
    /// the 5-ring aromaticity behavior of the (former) MDL software
    pub set_mdl_five_ring_aromaticity: bool,
    /// single bonds from aromatic or conjugated atoms to degree-one
    /// neighbors become SINGLE|AROMATIC
    pub adjust_single_bonds_to_degree_one_neighbors: bool,
    /// non-ring single bonds between aromatic or conjugated atoms become
    /// SINGLE|AROMATIC
    pub adjust_single_bonds_between_aromatic_atoms: bool,
}

impl Default for AdjustQueryParameters {
    fn default() -> Self {
        AdjustQueryParameters {
            adjust_degree: true,
            adjust_degree_flags: ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS,
            adjust_ring_count: false,
            adjust_ring_count_flags: ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS,
            make_dummies_queries: true,
            aromatize_if_possible: true,
            make_bonds_generic: false,
            make_bonds_generic_flags: ADJUST_IGNORENONE,
            make_atoms_generic: false,
            make_atoms_generic_flags: ADJUST_IGNORENONE,
            adjust_heavy_degree: false,
            adjust_heavy_degree_flags: ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS,
            adjust_ring_chain: false,
            adjust_ring_chain_flags: ADJUST_IGNORENONE,
            use_stereo_care_for_bonds: false,
            adjust_conjugated_five_rings: false,
            set_mdl_five_ring_aromaticity: false,
            adjust_single_bonds_to_degree_one_neighbors: false,
            adjust_single_bonds_between_aromatic_atoms: false,
        }
    }
}

impl AdjustQueryParameters {
    /// `AdjustQueryParameters::noAdjustments`: everything disabled
    pub fn no_adjustments() -> Self {
        AdjustQueryParameters {
            adjust_degree: false,
            make_dummies_queries: false,
            aromatize_if_possible: false,
            ..Default::default()
        }
    }
}

fn is_mapped(atom: &Atom) -> bool {
    atom.has_prop(cp::molAtomMapNumber)
}

/// Does an adjustment with these `flags` apply to an atom with the given
/// ring membership, (original) atomic number and mapping?
fn applies(flags: u32, in_rings: bool, atomic_num: u32, mapped: bool) -> bool {
    !((flags & ADJUST_IGNORECHAINS != 0 && !in_rings)
        || (flags & ADJUST_IGNORERINGS != 0 && in_rings)
        || (flags & ADJUST_IGNOREDUMMIES != 0 && atomic_num == 0)
        || (flags & ADJUST_IGNORENONDUMMIES != 0 && atomic_num != 0)
        || (flags & ADJUST_IGNOREMAPPED != 0 && mapped))
}

fn ensure_symm_sssr(mol: &mut RWMol) -> Result<()> {
    if !mol.ring_info().is_symm_sssr() {
        symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    }
    Ok(())
}

fn query_bond(query: BondQuery) -> Bond {
    let mut bond = Bond::new(BondType::Unspecified);
    bond.set_query(Some(query));
    bond
}

fn query_atom(query: AtomQuery) -> Atom {
    let mut atom = Atom::new(0);
    atom.set_query(Some(query));
    atom
}

/// `BOND_OR_QUERY` of two bond orders
fn order_or_aromatic(order: BondType) -> Bond {
    let mut or = Query::new(QueryKind::Or);
    or.children = vec![
        make_bond_order_equals_query(order),
        make_bond_order_equals_query(BondType::Aromatic),
    ];
    or.description = "BondOr".into();
    query_bond(or)
}

/// Five-rings with at least three conjugated bonds match both aromatic and
/// aliphatic rings: their plain bonds get a SINGLE|DOUBLE|AROMATIC query.
fn adjust_conjugated_five_rings(mol: &mut RWMol) -> Result<()> {
    ensure_symm_sssr(mol)?;
    let rings: Vec<Vec<i32>> = mol.ring_info().bond_rings().to_vec();
    let replacement = query_bond(make_single_or_double_or_aromatic_bond_query());
    for ring in rings.iter().filter(|r| r.len() == 5) {
        let conjugated = ring
            .iter()
            .filter(|&&bi| mol.bond(bi as u32).is_conjugated())
            .count();
        if conjugated < 3 {
            continue;
        }
        for &bi in ring {
            let bi = bi as u32;
            let (begin, end) = (mol.bond(bi).begin_atom_idx(), mol.bond(bi).end_atom_idx());
            for atom in [begin, end] {
                mol.atom_mut(atom)
                    .set_prop(CONJUGATED_OR_AROMATIC, 1i32, true)?;
            }
            if matches!(
                mol.bond(bi).bond_type(),
                BondType::Single | BondType::Double | BondType::Aromatic
            ) {
                if mol.bond(bi).has_query() {
                    rd_log_warning!(
                        "adjustConjugatedFiveRings: replacing a bond that already has a query"
                    );
                }
                mol.replace_bond(bi, &replacement, false, false)?;
            }
        }
    }
    Ok(())
}

fn is_aromatic_or_conjugated(atom: &Atom) -> bool {
    atom.is_aromatic() || atom.has_prop(CONJUGATED_OR_AROMATIC)
}

/// Single bonds from aromatic atoms may match aromatic bonds: either to
/// degree-one neighbors or between two aromatic atoms.
fn adjust_single_bonds_from_aromatic_atoms(
    mol: &mut RWMol,
    to_degree_one_neighbors: bool,
    between_aromatic_atoms: bool,
) -> Result<()> {
    if !to_degree_one_neighbors && !between_aromatic_atoms {
        return Ok(());
    }
    let replacement = query_bond(make_single_or_aromatic_bond_query());
    ensure_symm_sssr(mol)?;
    for bi in 0..mol.num_bonds() {
        let bond = mol.bond(bi);
        if bond.has_query() || bond.bond_type() != BondType::Single {
            continue;
        }
        let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let begin_aromatic = is_aromatic_or_conjugated(mol.atom(begin));
        let end_aromatic = is_aromatic_or_conjugated(mol.atom(end));
        let replace = if to_degree_one_neighbors && begin_aromatic != end_aromatic {
            (begin_aromatic && mol.atom_degree(end) == 1)
                || (end_aromatic && mol.atom_degree(begin) == 1)
        } else {
            between_aromatic_atoms && begin_aromatic && end_aromatic
        };
        if replace {
            mol.replace_bond(bi, &replacement, false, false)?;
        }
    }
    Ok(())
}

/// Aromatic five-rings of carbons and one "A" atom from a CTAB match both
/// aromatic and aliphatic rings: `A1:c:c:c:c:1` becomes
/// `A1-,:c=,:c-,:c=,:c-,:1`.
fn set_mdl_aromaticity(mol: &mut RWMol) -> Result<()> {
    ensure_symm_sssr(mol)?;
    let rings: Vec<Vec<i32>> = mol.ring_info().atom_rings().to_vec();
    for ring in rings.iter().filter(|r| r.len() == 5) {
        let ring: Vec<u32> = ring.iter().map(|&a| a as u32).collect();
        if let Some(dummy) = mdl_five_ring_dummy(mol, &ring)? {
            replace_mdl_ring_bonds(mol, &ring, ring[dummy])?;
        }
    }
    Ok(())
}

/// The position of the single "A" atom of a ring that qualifies for MDL
/// aromaticity (fully aromatic, only C and "A", no query bonds).
fn mdl_five_ring_dummy(mol: &RWMol, ring: &[u32]) -> Result<Option<usize>> {
    let mut dummy = None;
    for (i, &ai) in ring.iter().enumerate() {
        let atom = mol.atom(ai);
        if !atom.is_aromatic() {
            return Ok(None);
        }
        if atom.atomic_num() == 0 && atom.query().is_some_and(|q| q.type_label == "A") {
            if dummy.replace(i).is_some() {
                return Ok(None); // a second dummy
            }
        } else if atom.atomic_num() != 6 {
            return Ok(None);
        }
        let previous = ring[(i + ring.len() - 1) % ring.len()];
        let bond = mol
            .get_bond_between_atoms(ai, previous)?
            .ok_or_else(|| RDKitError::invariant("expected bond not found"))?;
        if bond.has_query() {
            return Ok(None);
        }
    }
    Ok(dummy)
}

fn replace_between(mol: &mut RWMol, a: u32, b: u32, with: &Bond) -> Result<()> {
    let idx = mol
        .get_bond_between_atoms(a, b)?
        .map(Bond::idx)
        .ok_or_else(|| RDKitError::invariant("expected bond not found"))?;
    mol.replace_bond(idx, with, false, false)
}

/// Three layers: the dummy (0), its two neighbors (1) and the two atoms that
/// are not connected to it (2).
fn replace_mdl_ring_bonds(mol: &mut RWMol, ring: &[u32], l0: u32) -> Result<()> {
    let (mut l1, mut l2) = (Vec::new(), Vec::new());
    for &ai in ring.iter().filter(|&&ai| ai != l0) {
        if mol.get_bond_between_atoms(ai, l0)?.is_some() {
            l1.push(ai);
        } else {
            l2.push(ai);
        }
    }
    if l1.len() != 2 || l2.len() != 2 {
        return Err(RDKitError::invariant("bad layer size"));
    }
    let single_aromatic = order_or_aromatic(BondType::Single);
    let double_aromatic = order_or_aromatic(BondType::Double);
    for &ai in &l1 {
        replace_between(mol, ai, l0, &single_aromatic)?;
        let partner = if mol.get_bond_between_atoms(ai, l2[0])?.is_some() {
            l2[0]
        } else {
            l2[1]
        };
        replace_between(mol, ai, partner, &double_aromatic)?;
    }
    replace_between(mol, l2[0], l2[1], &single_aromatic)
}

/// Makes sure atom `idx` carries a query (copying the atom's properties) and
/// returns it.
fn ensure_query_atom(mol: &mut RWMol, idx: u32) -> Result<&mut Atom> {
    if !mol.atom(idx).has_query() {
        let replacement = query_atom_from_atom(mol.atom(idx))?;
        mol.replace_atom(idx, &replacement, true)?;
    }
    Ok(mol.atom_mut(idx))
}

/// Everything the per-atom adjustments need to know about atom `idx` before
/// any of them modifies it.
struct AtomContext {
    in_rings: u32,
    atomic_num: u32,
}

impl AtomContext {
    fn applies(&self, mol: &RWMol, idx: u32, flags: u32) -> bool {
        applies(
            flags,
            self.in_rings > 0,
            self.atomic_num,
            is_mapped(mol.atom(idx)),
        )
    }
}

fn adjust_atoms(mol: &mut RWMol, params: &AdjustQueryParameters, orig_nums: &[u32]) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        // the atom's properties may be unavailable once it is a query atom
        let ctx = AtomContext {
            in_rings: mol.ring_info().num_atom_rings(idx)?,
            atomic_num: orig_nums[idx as usize],
        };
        if params.make_dummies_queries
            && ctx.atomic_num == 0
            && !mol.atom(idx).has_query()
            && mol.atom(idx).isotope() == 0
        {
            mol.replace_atom(idx, &query_atom(make_atom_null_query()), true)?;
        }
        if params.adjust_degree && ctx.applies(mol, idx, params.adjust_degree_flags) {
            let degree = mol.atom_degree(idx) as i32;
            expand(mol, idx, make_atom_explicit_degree_query(degree))?;
        }
        if params.adjust_heavy_degree && ctx.applies(mol, idx, params.adjust_heavy_degree_flags) {
            let heavy = mol.total_degree(idx)? as i32 - mol.total_num_hs(idx, true)? as i32;
            expand(mol, idx, make_atom_heavy_atom_degree_query(heavy))?;
        }
        if params.adjust_ring_count && ctx.applies(mol, idx, params.adjust_ring_count_flags) {
            expand(mol, idx, make_atom_in_n_rings_query(ctx.in_rings as i32))?;
        }
        if params.adjust_ring_chain && ctx.applies(mol, idx, params.adjust_ring_chain_flags) {
            let mut in_ring = make_atom_in_ring_query();
            in_ring.negate = ctx.in_rings == 0;
            expand(mol, idx, in_ring)?;
        }
    }
    Ok(())
}

fn expand(mol: &mut RWMol, idx: u32, query: AtomQuery) -> Result<()> {
    let atom = ensure_query_atom(mol, idx)?;
    atom_expand_query(atom, query, CompositeQueryType::And, true)
}

/// `parseWhichString`: the flags named in a `|`-separated string.
fn parse_which_string(txt: &str) -> Result<u32> {
    let mut res = ADJUST_IGNORENONE;
    for token in txt.split('|').filter(|t| !t.is_empty()) {
        res |= match token {
            "IGNORENONE" => ADJUST_IGNORENONE,
            "IGNORERINGS" => ADJUST_IGNORERINGS,
            "IGNORECHAINS" => ADJUST_IGNORECHAINS,
            "IGNOREDUMMIES" => ADJUST_IGNOREDUMMIES,
            "IGNORENONDUMMIES" => ADJUST_IGNORENONDUMMIES,
            "IGNOREALL" => ADJUST_IGNOREALL,
            _ => {
                return Err(RDKitError::Value(format!(
                    "Unknown flag value: '{token}'. Flags ignored."
                )));
            }
        };
    }
    Ok(res)
}

/// `parseAdjustQueryParametersFromJSON`: updates the parameters from a JSON
/// object. Keys that are absent leave the parameter as it was.
pub fn parse_adjust_query_parameters_from_json(
    params: &mut AdjustQueryParameters,
    json: &str,
) -> Result<()> {
    if json.is_empty() {
        return Err(RDKitError::precondition("empty JSON provided"));
    }
    let pt = PTree::from_json(json)?;
    macro_rules! bools {
        ($($key:literal => $field:ident),* $(,)?) => {
            $( params.$field = pt.get_bool($key, params.$field); )*
        };
    }
    bools! {
        "adjustDegree" => adjust_degree,
        "adjustRingCount" => adjust_ring_count,
        "makeDummiesQueries" => make_dummies_queries,
        "aromatizeIfPossible" => aromatize_if_possible,
        "makeBondsGeneric" => make_bonds_generic,
        "makeAtomsGeneric" => make_atoms_generic,
        "adjustHeavyDegree" => adjust_heavy_degree,
        "adjustRingChain" => adjust_ring_chain,
        "useStereoCareForBonds" => use_stereo_care_for_bonds,
        "adjustConjugatedFiveRings" => adjust_conjugated_five_rings,
        "setMDLFiveRingAromaticity" => set_mdl_five_ring_aromaticity,
        "adjustSingleBondsToDegreeOneNeighbors" => adjust_single_bonds_to_degree_one_neighbors,
        "adjustSingleBondsBetweenAromaticAtoms" => adjust_single_bonds_between_aromatic_atoms,
    }
    macro_rules! flags {
        ($($key:literal => $field:ident),* $(,)?) => {
            $(
                let which = pt.get_str($key).unwrap_or("").to_uppercase();
                if !which.is_empty() {
                    params.$field = parse_which_string(&which)?;
                }
            )*
        };
    }
    flags! {
        "adjustDegreeFlags" => adjust_degree_flags,
        "adjustHeavyDegreeFlags" => adjust_heavy_degree_flags,
        "adjustRingCountFlags" => adjust_ring_count_flags,
        "makeBondsGenericFlags" => make_bonds_generic_flags,
        "makeAtomsGenericFlags" => make_atoms_generic_flags,
        "adjustRingChainFlags" => adjust_ring_chain_flags,
    }
    Ok(())
}

/// `MolOps::adjustQueryProperties(mol, params)` in place.
pub fn adjust_query_properties_in_place(
    mol: &mut RWMol,
    params: &AdjustQueryParameters,
) -> Result<()> {
    if params.aromatize_if_possible {
        sanitize_mol(mol, SANITIZE_SYMMRINGS | SANITIZE_SETAROMATICITY).map_err(|(_, e)| e)?;
    } else {
        ensure_symm_sssr(mol)?;
    }
    let orig_nums: Vec<u32> = mol.atoms().iter().map(Atom::atomic_num).collect();

    if params.make_atoms_generic {
        for idx in 0..mol.num_atoms() {
            let in_rings = mol.ring_info().num_atom_rings(idx)? > 0;
            let flags = params.make_atoms_generic_flags;
            if !((flags & ADJUST_IGNORECHAINS != 0 && !in_rings)
                || (flags & ADJUST_IGNORERINGS != 0 && in_rings)
                || (flags & ADJUST_IGNOREMAPPED != 0 && is_mapped(mol.atom(idx))))
            {
                mol.replace_atom(idx, &query_atom(make_atom_null_query()), true)?;
            }
        }
    }
    if params.make_bonds_generic {
        for idx in 0..mol.num_bonds() {
            let in_rings = mol.ring_info().num_bond_rings(idx)? > 0;
            let flags = params.make_bonds_generic_flags;
            if !((flags & ADJUST_IGNORECHAINS != 0 && !in_rings)
                || (flags & ADJUST_IGNORERINGS != 0 && in_rings))
            {
                let generic = query_bond(make_bond_null_query());
                mol.replace_bond(idx, &generic, true, false)?;
            }
        }
    }
    adjust_atoms(mol, params, &orig_nums)?;
    if params.use_stereo_care_for_bonds {
        for bi in 0..mol.num_bonds() {
            let bond = mol.bond_mut(bi);
            let stereo_care = bond
                .props
                .get(cp::molStereoCare)
                .and_then(crate::rdgeneral::RDValue::as_int)
                .is_some_and(|v| v != 0);
            if bond.bond_type() == BondType::Double
                && bond.stereo() > BondStereo::StereoAny
                && !stereo_care
            {
                bond.set_stereo(BondStereo::StereoNone)?;
            }
        }
    }
    if params.set_mdl_five_ring_aromaticity {
        set_mdl_aromaticity(mol)?;
    }
    if params.adjust_conjugated_five_rings {
        adjust_conjugated_five_rings(mol)?;
    }
    adjust_single_bonds_from_aromatic_atoms(
        mol,
        params.adjust_single_bonds_to_degree_one_neighbors,
        params.adjust_single_bonds_between_aromatic_atoms,
    )?;
    if params.set_mdl_five_ring_aromaticity {
        // (upstream only cleans up in this case)
        for idx in 0..mol.num_atoms() {
            mol.atom_mut(idx).clear_prop(CONJUGATED_OR_AROMATIC);
        }
    }
    Ok(())
}

/// `MolOps::adjustQueryProperties(mol, params)`: a modified copy.
pub fn adjust_query_properties(mol: &ROMol, params: &AdjustQueryParameters) -> Result<RWMol> {
    let mut res = mol.clone();
    adjust_query_properties_in_place(&mut res, params)?;
    Ok(res)
}
