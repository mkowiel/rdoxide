//! Ported from RDKit `ChemReactions/ReactionRunner.{h,cpp}`: running a
//! reaction on reactant molecules to build the products.
//!
//! The entry points are [`run_reactants`] (all combinations of template
//! matches), [`run_reactant`] (one reactant against one template) and
//! [`run_reactant_in_place`] (single-reactant, single-product reactions that
//! only modify the reactant). Products of different match combinations are
//! independent, so they are built on up to `num_threads` threads.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::ChemicalReaction;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::chirality::bond_dirs::{
    find_stereo_atoms, has_stereo_bond_dir, neighboring_directed_bond,
    set_double_bond_neighbor_directions,
};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::StereoGroup;
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::threads::par_map;
use crate::rdgeneral::utils::count_swaps_to_interconvert;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_warning};

/// Marks product atoms that were dummies in the product template.
const WAS_DUMMY: &str = "was_dummy";
/// Marks product double bonds whose stereo is unknown in the template.
const UNKNOWN_STEREO_RXN_BOND: &str = "_UnknownStereoRxnBond";

fn reaction_error(message: &str) -> RDKitError {
    RDKitError::ChemicalReaction(message.into())
}

/// `ReactantProductAtomMapping`
struct AtomMapping {
    mapped_atoms: Vec<bool>,
    skipped_atoms: Vec<bool>,
    /// reactant atom -> the product atoms it maps to
    react_prod: BTreeMap<u32, Vec<u32>>,
    prod_react: BTreeMap<u32, u32>,
    /// pairs of (sorted) atom-map numbers that are bonded in the template
    template_bonds: BTreeSet<(i32, i32)>,
}

impl AtomMapping {
    /// `reactProdAtomMap[idx]` (which inserts an empty list when missing)
    fn products_of(&mut self, reactant_atom: u32) -> Vec<u32> {
        self.react_prod.entry(reactant_atom).or_default().clone()
    }
}

// ---------------------------------------------------------------- matching

/// `getReactantMatchesToTemplate`: all matches of the template, minus those
/// that touch protected atoms.
fn reactant_matches_to_template(
    reactant: &ROMol,
    templ: &ROMol,
    max_matches: u32,
    ssparams: &SubstructMatchParameters,
) -> Result<Vec<MatchVect>> {
    let mut params = ssparams.clone();
    params.uniquify = false;
    params.max_matches = max_matches;
    let matches = substruct_match(reactant, templ, &params)?;
    Ok(matches
        .into_iter()
        .filter(|m| {
            !m.iter()
                .any(|&(_, ri)| reactant.atom(ri as u32).has_prop(cp::_protected))
        })
        .collect())
}

/// `getReactantMatches`: the matches of every (or one) reactant against its
/// template; `None` when a reactant has no match.
fn reactant_matches(
    reactants: &[&ROMol],
    rxn: &ChemicalReaction,
    max_matches: u32,
    match_single: Option<usize>,
) -> Result<Option<Vec<Vec<MatchVect>>>> {
    if reactants.len() != rxn.num_reactant_templates() {
        return Err(RDKitError::precondition("reactant size mismatch"));
    }
    let mut by_reactant = vec![Vec::new(); reactants.len()];
    for (i, templ) in rxn.reactants().iter().enumerate() {
        if match_single.is_none() || match_single == Some(i) {
            let matches = reactant_matches_to_template(
                reactants[i],
                templ,
                max_matches,
                rxn.substruct_params(),
            )?;
            if matches.is_empty() {
                return Ok(None);
            }
            by_reactant[i] = matches;
        }
    }
    Ok(Some(by_reactant))
}

/// `recurseOverReactantCombinations`; returns false when the product limit
/// stopped the enumeration.
fn recurse_over_combinations(
    by_reactant: &[Vec<MatchVect>],
    out: &mut Vec<Vec<MatchVect>>,
    level: usize,
    combination: &[MatchVect],
    max_products: u32,
) -> bool {
    let limit_hit =
        |out: &Vec<Vec<MatchVect>>| max_products > 0 && out.len() >= max_products as usize;
    if limit_hit(out) {
        return false;
    }
    let mut keep_going = true;
    for m in &by_reactant[level] {
        let mut prod = combination.to_vec();
        prod[level] = m.clone();
        if level == by_reactant.len() - 1 {
            if limit_hit(out) {
                keep_going = false;
                break;
            }
            out.push(prod);
        } else {
            keep_going =
                recurse_over_combinations(by_reactant, out, level + 1, &prod, max_products);
        }
    }
    keep_going
}

/// `generateReactantCombinations`
fn reactant_combinations(by_reactant: &[Vec<MatchVect>], max_products: u32) -> Vec<Vec<MatchVect>> {
    let mut out = Vec::new();
    let start = vec![MatchVect::new(); by_reactant.len()];
    if !recurse_over_combinations(by_reactant, &mut out, 0, &start, max_products) {
        rd_log_warning!("Maximum product count hit {max_products}, stopping reaction early...\n");
    }
    out
}

// ------------------------------------------------------- product building

/// `updatePropsFromImplicitProps`: applies the charge / H count / isotope
/// queries recorded on a template atom to `atom`; true if `atom` changed.
fn update_props_from_implicit_props(template: &Atom, atom: &mut Atom) -> bool {
    let mut res = false;
    if let Some(val) = template
        .get_prop(cp::_QueryFormalCharge)
        .and_then(RDValue::as_int)
        && val != atom.formal_charge()
    {
        atom.set_formal_charge(val);
        res = true;
    }
    if let Some(val) = template
        .get_prop(cp::_QueryHCount)
        .and_then(RDValue::as_uint)
        && (!atom.no_implicit() || atom.num_explicit_hs() != val)
    {
        atom.set_num_explicit_hs(val);
        atom.set_no_implicit(true); // github #1544
        res = true;
    }
    if let Some(val) = template
        .get_prop(cp::_QueryIsotope)
        .and_then(RDValue::as_uint)
        && val != atom.isotope()
    {
        atom.set_isotope(val);
        res = true;
    }
    res
}

/// `updateImplicitAtomProperties`: copies the properties the product
/// template leaves unspecified from the reactant atom.
fn update_implicit_atom_properties(prod: &mut Atom, react: &Atom) {
    if prod.atomic_num() != react.atomic_num() {
        return;
    }
    if !prod.has_prop(cp::_QueryFormalCharge) {
        prod.set_formal_charge(react.formal_charge());
    }
    if !prod.has_prop(cp::_QueryIsotope) {
        prod.set_isotope(react.isotope());
    }
    if !prod.has_prop(cp::_ReactionDegreeChanged) && !prod.has_prop(cp::_QueryHCount) {
        prod.set_num_explicit_hs(react.num_explicit_hs());
        prod.set_no_implicit(react.no_implicit());
    }
}

fn map_number(atom: &Atom) -> Option<i32> {
    atom.get_prop(cp::molAtomMapNumber)
        .and_then(RDValue::as_int)
}

/// `convertTemplateToMol`: a plain molecule from a product template, with
/// the template's atom-map numbers moved to bookmarks (and `old_mapno`).
fn convert_template_to_mol(templ: &ROMol) -> Result<RWMol> {
    let mut res = RWMol::new();
    for oatom in templ.atoms() {
        let mut new_atom = oatom.copy_detached();
        // upstream copies the Atom part only: queries are dropped
        new_atom.set_query(None);
        let map_num = map_number(&new_atom);
        if map_num.is_some() {
            new_atom.clear_prop(cp::molAtomMapNumber);
        }
        new_atom.set_chiral_tag(ChiralType::Unspecified);
        if oatom
            .get_prop(cp::molInversionFlag)
            .and_then(RDValue::as_int)
            == Some(4)
        {
            new_atom.set_chiral_tag(oatom.chiral_tag());
        }
        let template_copy = new_atom.clone();
        update_props_from_implicit_props(&template_copy, &mut new_atom);
        let idx = res.add_atom(new_atom, false);
        if let Some(n) = map_num {
            res.set_atom_bookmark(idx, n);
            res.atom_mut(idx).dict_mut().set_val(cp::reactionMapNum, n);
        }
    }
    for old in templ.bonds() {
        let n = res.add_bond(old.begin_atom_idx(), old.end_atom_idx(), old.bond_type())?;
        let new_idx = n - 1;
        let (begin_aromatic, end_aromatic) = (
            res.atom(old.begin_atom_idx()).is_aromatic(),
            res.atom(old.end_atom_idx()).is_aromatic(),
        );
        let new_b = res.bond_mut(new_idx);
        new_b.set_bond_dir(old.bond_dir());
        if let Some(q) = old.query() {
            let descr = q.description();
            if descr == "BondOr" && old.bond_type() == BondType::Single {
                if begin_aromatic && end_aromatic {
                    new_b.set_bond_type(BondType::Aromatic);
                    new_b.set_is_aromatic(true);
                } else {
                    new_b.set_bond_type(BondType::Single);
                    new_b.set_is_aromatic(false);
                }
            } else if descr == "BondNull" {
                new_b.dict_mut().set_val(cp::NullBond, 1i32);
            }
        }
        if old.bond_type() == BondType::Double {
            let (start, end) = (old.begin_atom_idx(), old.end_atom_idx());
            if templ.atom_degree(start) > 1
                && templ.atom_degree(end) > 1
                && (neighboring_directed_bond(templ, start).is_none()
                    || neighboring_directed_bond(templ, end).is_none())
            {
                res.bond_mut(new_idx)
                    .dict_mut()
                    .set_val(UNKNOWN_STEREO_RXN_BOND, 1i32);
            }
        }
        res.bond_mut(new_idx).update_props(old.dict(), true);
    }
    Ok(res)
}

/// `getAtomMappingsReactantProduct`
fn atom_mappings(
    match_: &MatchVect,
    reactant_template: &ROMol,
    product: &RWMol,
    num_react_atoms: u32,
) -> Result<AtomMapping> {
    let mut mapping = AtomMapping {
        mapped_atoms: vec![false; num_react_atoms as usize],
        skipped_atoms: vec![false; num_react_atoms as usize],
        react_prod: BTreeMap::new(),
        prod_react: BTreeMap::new(),
        template_bonds: BTreeSet::new(),
    };
    for bond in reactant_template.bonds() {
        let a = reactant_template.atom(bond.begin_atom_idx()).atom_map_num();
        let b = reactant_template.atom(bond.end_atom_idx()).atom_map_num();
        mapping.template_bonds.insert((a.min(b), a.max(b)));
    }
    for &(ti, ri) in match_ {
        let ri = ri as u32;
        let template_atom = reactant_template.atom(ti as u32);
        match map_number(template_atom) {
            Some(n) if product.has_atom_bookmark(n) => {
                for &p in product.get_all_atoms_with_bookmark(n)? {
                    mapping.react_prod.entry(ri).or_default().push(p);
                    mapping.mapped_atoms[ri as usize] = true;
                    if p >= product.num_atoms() {
                        return Err(RDKitError::invariant("yikes!"));
                    }
                    mapping.prod_react.insert(p, ri);
                }
            }
            _ => mapping.skipped_atoms[ri as usize] = true,
        }
    }
    Ok(mapping)
}

/// The product atom bonded to `atom` among `candidates`
/// (`reactProdMapAnchorIdx`).
fn anchor_idx(product: &ROMol, atom: u32, candidates: &[u32]) -> Result<u32> {
    if candidates.len() == 1 {
        return Ok(candidates[0]);
    }
    candidates
        .iter()
        .copied()
        .find(|&c| product.bond_idx_between(atom, c).is_some())
        .ok_or_else(|| RDKitError::invariant("match not found"))
}

/// A stereo bond end: the stereo ("anchor") neighbor and the other,
/// non-anchor neighbor if there is one (`StereoBondEndCap`).
struct StereoBondEndCap {
    anchor: u32,
    non_anchor: Option<u32>,
}

impl StereoBondEndCap {
    fn new(mol: &ROMol, atom: u32, other_double_bond_atom: u32, anchor: u32) -> Result<Self> {
        if mol.total_degree(atom)? > 3 {
            return Err(RDKitError::precondition(
                "Stereo Bond extremes must have less than four neighbors",
            ));
        }
        let non_anchor = mol
            .atom_neighbors(atom)
            .find(|&n| n != other_double_bond_atom && n != anchor);
        Ok(StereoBondEndCap { anchor, non_anchor })
    }

    /// The product atoms that can stand for this end's anchor, and whether
    /// the stereo label flips because the non-anchor neighbor is used.
    fn product_anchor_candidates(&self, mapping: &AtomMapping) -> (Vec<u32>, bool) {
        if let Some(c) = mapping.react_prod.get(&self.anchor) {
            return (c.clone(), false);
        }
        if let Some(c) = self.non_anchor.and_then(|n| mapping.react_prod.get(&n)) {
            return (c.clone(), true);
        }
        (Vec::new(), false)
    }
}

/// `forwardReactantBondStereo`
fn forward_reactant_bond_stereo(
    mapping: &AtomMapping,
    product: &mut RWMol,
    p_bond: u32,
    reactant: &ROMol,
    r_bond: u32,
) -> Result<()> {
    let rb = reactant.bond(r_bond);
    if rb.stereo() <= BondStereo::StereoAny {
        return Err(RDKitError::precondition(
            "bond in reactant must have defined stereo",
        ));
    }
    let (r_start, r_end) = (rb.begin_atom_idx(), rb.end_atom_idx());
    let stereo_atoms = find_stereo_atoms(reactant, r_bond)?;
    if stereo_atoms.len() != 2 {
        rd_log_warning!(
            "WARNING: neither stereo atoms nor CIP codes found for double bond. Stereochemistry info will not be propagated to product.\n"
        );
        product
            .bond_mut(p_bond)
            .set_stereo(BondStereo::StereoNone)?;
        return Ok(());
    }
    let mut start = StereoBondEndCap::new(reactant, r_start, r_end, stereo_atoms[0] as u32)?;
    let mut end = StereoBondEndCap::new(reactant, r_end, r_start, stereo_atoms[1] as u32)?;
    let p_begin = product.bond(p_bond).begin_atom_idx();
    let p_end = product.bond(p_bond).end_atom_idx();
    let begin_react = mapping.prod_react.get(&p_begin).copied().unwrap_or(0);
    if begin_react == r_end {
        std::mem::swap(&mut start, &mut end);
    } else if begin_react != r_start {
        return Err(RDKitError::Runtime(
            "Reactant and Product bond ends do not match".into(),
        ));
    }
    let (start_candidates, start_swap) = start.product_anchor_candidates(mapping);
    let (end_candidates, end_swap) = end.product_anchor_candidates(mapping);
    if start_candidates.is_empty() || end_candidates.is_empty() {
        return Ok(());
    }
    let p_start_anchor = anchor_idx(product, p_begin, &start_candidates)?;
    let p_end_anchor = anchor_idx(product, p_end, &end_candidates)?;
    if product.bond_idx_between(p_begin, p_start_anchor).is_none()
        || product.bond_idx_between(p_end, p_end_anchor).is_none()
    {
        rd_log_warning!(
            "stereo atoms in input cannot be mapped to output (atoms are no longer bonded)\n"
        );
        return Ok(());
    }
    product.set_bond_stereo_atoms(p_bond, p_start_anchor, p_end_anchor)?;
    let flip = start_swap != end_swap;
    let cis = matches!(rb.stereo(), BondStereo::StereoCis | BondStereo::StereoZ);
    let stereo = if cis != flip {
        BondStereo::StereoCis
    } else {
        BondStereo::StereoTrans
    };
    product.bond_mut(p_bond).set_stereo(stereo)
}

/// `translateProductStereoBondDirections`
fn translate_product_stereo_bond_directions(
    product: &mut RWMol,
    p_bond: u32,
    start: u32,
    end: u32,
) -> Result<()> {
    let (pb_begin, pb_end) = (
        product.bond(p_bond).begin_atom_idx(),
        product.bond(p_bond).end_atom_idx(),
    );
    let (sb, eb) = (product.bond(start), product.bond(end));
    let p_start_anchor = sb.other_atom_idx(pb_begin)?;
    let p_end_anchor = eb.other_atom_idx(pb_end)?;
    let mut same_dir = sb.bond_dir() == eb.bond_dir();
    if sb.begin_atom_idx() == pb_begin {
        same_dir = !same_dir;
    }
    if eb.begin_atom_idx() != pb_end {
        same_dir = !same_dir;
    }
    product.set_bond_stereo_atoms(p_bond, p_start_anchor, p_end_anchor)?;
    product.bond_mut(p_bond).set_stereo(if same_dir {
        BondStereo::StereoTrans
    } else {
        BondStereo::StereoCis
    })
}

/// `updateStereoBonds`: double-bond stereo of the product, from the
/// template's bond directions or the mapped reactant bond.
fn update_stereo_bonds(product: &mut RWMol, reactant: &ROMol, mapping: &AtomMapping) -> Result<()> {
    for b in 0..product.num_bonds() {
        if product.bond(b).bond_type() != BondType::Double {
            continue;
        }
        if product.bond(b).has_prop(UNKNOWN_STEREO_RXN_BOND) {
            let bond = product.bond_mut(b);
            bond.set_stereo(BondStereo::StereoNone)?;
            bond.clear_prop(UNKNOWN_STEREO_RXN_BOND);
            continue;
        }
        let (begin, end) = (
            product.bond(b).begin_atom_idx(),
            product.bond(b).end_atom_idx(),
        );
        let start_dir = neighboring_directed_bond(product, begin);
        let end_dir = neighboring_directed_bond(product, end);
        if let (Some(s), Some(e)) = (start_dir, end_dir) {
            if !has_stereo_bond_dir(product.bond(s)) || !has_stereo_bond_dir(product.bond(e)) {
                return Err(RDKitError::precondition(
                    "Both neighboring bonds must have bond directions",
                ));
            }
            translate_product_stereo_bond_directions(product, b, s, e)?;
            continue;
        }
        let (Some(&rb_begin), Some(&rb_end)) =
            (mapping.prod_react.get(&begin), mapping.prod_react.get(&end))
        else {
            continue;
        };
        let Some(r_bond) = reactant.bond_idx_between(rb_begin, rb_end) else {
            continue;
        };
        let rb = reactant.bond(r_bond);
        if rb.bond_type() == BondType::Double {
            if rb.stereo() > BondStereo::StereoAny {
                forward_reactant_bond_stereo(mapping, product, b, reactant, r_bond)?;
            } else if rb.stereo() == BondStereo::StereoAny {
                product.bond_mut(b).set_stereo(BondStereo::StereoAny)?;
            }
        }
    }
    Ok(())
}

/// `setReactantBondPropertiesToProduct`: bonds the template leaves open
/// ("null" or query bonds) take their type, query and properties from the
/// mapped reactant bond.
fn set_reactant_bond_properties_to_product(
    product: &mut RWMol,
    reactant: &ROMol,
    mapping: &AtomMapping,
) -> Result<()> {
    for bidx in 0..product.num_bonds() {
        let (pb, pe) = (
            product.bond(bidx).begin_atom_idx(),
            product.bond(bidx).end_atom_idx(),
        );
        let (Some(&rb), Some(&re)) = (mapping.prod_react.get(&pb), mapping.prod_react.get(&pe))
        else {
            continue;
        };
        let Some(r_bond_idx) = reactant.bond_idx_between(rb, re) else {
            continue;
        };
        let r_bond = reactant.bond(r_bond_idx);
        let p_bond = product.bond(bidx);
        if !p_bond.has_prop(cp::NullBond)
            && !p_bond.has_prop(cp::_MolFileBondQuery)
            && !r_bond.has_query()
        {
            continue;
        }
        if let Some(q) = r_bond.query() {
            let mut replacement = Bond::new(r_bond.bond_type());
            replacement.set_query(Some(q.clone()));
            product.replace_bond(bidx, &replacement, false, true)?;
        } else {
            product.bond_mut(bidx).set_bond_type(r_bond.bond_type());
        }
        let p_bond = product.bond_mut(bidx);
        if r_bond.bond_type() == BondType::Double && r_bond.bond_dir() == BondDir::EitherDouble {
            p_bond.set_bond_dir(BondDir::EitherDouble);
        }
        p_bond.set_is_aromatic(r_bond.is_aromatic());
        p_bond.update_props(r_bond.dict(), false);
        p_bond.clear_prop(cp::NullBond);
    }
    Ok(())
}

/// `checkProductChirality`
fn check_product_chirality(reactant_chirality: ChiralType, product_atom: &mut Atom) {
    let flag = product_atom
        .get_prop(cp::molInversionFlag)
        .and_then(RDValue::as_int)
        .unwrap_or(0);
    match flag {
        0 | 2 => product_atom.set_chiral_tag(reactant_chirality),
        1 => {
            if reactant_chirality != ChiralType::TetrahedralCW
                && reactant_chirality != ChiralType::TetrahedralCCW
            {
                rd_log_warning!("unsupported chiral type on reactant atom ignored\n");
            } else {
                product_atom.set_chiral_tag(reactant_chirality);
                product_atom.invert_chirality();
            }
        }
        3 => product_atom.set_chiral_tag(ChiralType::Unspecified),
        4 => {}
        _ => rd_log_warning!(
            "unrecognized chiral inversion/retention flag on product atom ignored\n"
        ),
    }
}

/// `setReactantAtomPropertiesToProduct`
fn set_reactant_atom_properties_to_product(
    product_atom: &mut Atom,
    reactant_atom: &Atom,
    set_implicit_properties: bool,
    reactant_id: u32,
) {
    if product_atom.atomic_num() == 0 || product_atom.has_prop(cp::_MolFileAtomQuery) {
        product_atom.set_atomic_num(reactant_atom.atomic_num() as i32);
        product_atom.set_is_aromatic(reactant_atom.is_aromatic());
        if !set_implicit_properties {
            product_atom.set_isotope(reactant_atom.isotope());
        }
        product_atom.clear_prop(cp::dummyLabel);
        product_atom.clear_prop(cp::_MolFileRLabel);
        product_atom.dict_mut().set_val(WAS_DUMMY, true);
    } else {
        product_atom.clear_prop(WAS_DUMMY);
    }
    product_atom.dict_mut().set_val(
        cp::reactantAtomIdx,
        RDValue::UnsignedInt(reactant_atom.idx()),
    );
    product_atom
        .dict_mut()
        .set_val(cp::reactantIdx, RDValue::UnsignedInt(reactant_id));
    if set_implicit_properties {
        update_implicit_atom_properties(product_atom, reactant_atom);
    }
    if is_specified_chirality(reactant_atom.chiral_tag())
        && product_atom.has_prop(cp::molInversionFlag)
    {
        check_product_chirality(reactant_atom.chiral_tag(), product_atom);
    }
    if let Some(info) = reactant_atom.monomer_info() {
        product_atom.set_monomer_info(Some(info.clone()));
    }
}

fn is_specified_chirality(tag: ChiralType) -> bool {
    tag != ChiralType::Unspecified && tag != ChiralType::Other
}

/// `addBondToProduct`: the new bond's index.
fn add_bond_to_product(orig: &Bond, product: &mut RWMol, begin: u32, end: u32) -> Result<u32> {
    let n = match orig.query() {
        None => product.add_bond(begin, end, orig.bond_type())?,
        Some(q) => {
            let mut bond = Bond::new(orig.bond_type());
            bond.set_begin_atom_idx_unchecked(begin);
            bond.set_end_atom_idx_unchecked(end);
            bond.set_query(Some(q.clone()));
            product.add_bond_object(bond)?
        }
    };
    Ok(n - 1)
}

/// `addMissingProductBonds`
fn add_missing_product_bonds(
    orig: &Bond,
    product: &mut RWMol,
    mapping: &mut AtomMapping,
) -> Result<()> {
    let begin = mapping.products_of(orig.begin_atom_idx());
    let end = mapping.products_of(orig.end_atom_idx());
    if begin.len() != end.len() {
        return Err(RDKitError::invariant(
            "Different number of start-end points for product bonds.",
        ));
    }
    for (b, e) in begin.into_iter().zip(end) {
        add_bond_to_product(orig, product, b, e)?;
    }
    Ok(())
}

/// `addMissingProductAtom`: copies an unmapped reactant atom (and its bond
/// to the already present neighbor) into the product.
fn add_missing_product_atom(
    react_atom: u32,
    react_neighbor: u32,
    prod_neighbor: u32,
    product: &mut RWMol,
    reactant: &ROMol,
    mapping: &mut AtomMapping,
    reactant_id: u32,
) -> Result<()> {
    let mut new_atom = reactant.atom(react_atom).copy_detached();
    new_atom
        .dict_mut()
        .set_val(cp::reactantAtomIdx, RDValue::UnsignedInt(react_atom));
    new_atom
        .dict_mut()
        .set_val(cp::reactantIdx, RDValue::UnsignedInt(reactant_id));
    let product_idx = product.add_atom(new_atom, false);
    mapping
        .react_prod
        .entry(react_atom)
        .or_default()
        .push(product_idx);
    mapping.prod_react.insert(product_idx, react_atom);
    let orig = reactant
        .bond_idx_between(react_neighbor, react_atom)
        .map(|b| reactant.bond(b))
        .ok_or_else(|| RDKitError::invariant("missing reactant bond"))?;
    let (mut begin, mut end) = (product_idx, prod_neighbor);
    if orig.begin_atom_idx() == react_neighbor {
        std::mem::swap(&mut begin, &mut end);
    }
    let b = add_bond_to_product(orig, product, begin, end)?;
    let prod_bond = product.bond_mut(b);
    if orig.bond_type() == BondType::Double && orig.bond_dir() == BondDir::EitherDouble {
        prod_bond.set_bond_dir(BondDir::EitherDouble);
    }
    prod_bond.update_props(orig.dict(), true);
    Ok(())
}

/// `addReactantNeighborsToProduct`: walks out from a mapped reactant atom
/// and copies everything not covered by the template into the product.
fn add_reactant_neighbors_to_product(
    reactant: &ROMol,
    reactant_atom: u32,
    product: &mut RWMol,
    visited: &mut [bool],
    chiral_atoms_to_check: &mut Vec<u32>,
    mapping: &mut AtomMapping,
    reactant_id: u32,
) -> Result<()> {
    let mut stack: VecDeque<u32> = VecDeque::new();
    stack.push_back(reactant_atom);
    while let Some(current) = stack.pop_front() {
        let current_products = mapping.react_prod.get(&current).cloned().ok_or_else(|| {
            RDKitError::invariant("reactant atom on traversal stack not present in product.")
        })?;
        visited[current as usize] = true;
        for nbr in reactant.atom_neighbors(current).collect::<Vec<_>>() {
            if visited[nbr as usize] || mapping.skipped_atoms[nbr as usize] {
                continue;
            }
            let orig_bond = || -> Result<&Bond> {
                reactant
                    .bond_idx_between(current, nbr)
                    .map(|b| reactant.bond(b))
                    .ok_or_else(|| RDKitError::invariant("missing reactant bond"))
            };
            if mapping.mapped_atoms[nbr as usize] {
                if !mapping.mapped_atoms[current as usize] {
                    if !mapping.react_prod.contains_key(&nbr) {
                        return Err(RDKitError::invariant(
                            "reactant atom not present in product.",
                        ));
                    }
                    add_missing_product_bonds(orig_bond()?, product, mapping)?;
                } else {
                    let prod_begin = mapping.react_prod[&current][0];
                    let prod_end = mapping.react_prod[&nbr][0];
                    if product.bond_idx_between(prod_begin, prod_end).is_none() {
                        let a1 = product
                            .atom(prod_begin)
                            .get_prop(cp::reactionMapNum)
                            .and_then(RDValue::as_int);
                        let a2 = product
                            .atom(prod_end)
                            .get_prop(cp::reactionMapNum)
                            .and_then(RDValue::as_int);
                        let (Some(a1), Some(a2)) = (a1, a2) else {
                            return Err(RDKitError::invariant("atoms should be mapped in product"));
                        };
                        if !mapping.template_bonds.contains(&(a1.min(a2), a1.max(a2))) {
                            add_missing_product_bonds(orig_bond()?, product, mapping)?;
                        }
                    }
                }
            } else if mapping.react_prod.contains_key(&nbr) {
                add_missing_product_bonds(orig_bond()?, product, mapping)?;
            } else {
                for &i in &current_products {
                    add_missing_product_atom(
                        nbr,
                        current,
                        i,
                        product,
                        reactant,
                        mapping,
                        reactant_id,
                    )?;
                }
                stack.push_back(nbr);
                if reactant.atom(nbr).chiral_tag() != ChiralType::Unspecified {
                    chiral_atoms_to_check.push(nbr);
                }
            }
        }
    }
    Ok(())
}

/// `checkAndCorrectChiralityOfMatchingAtomsInProduct`: re-derives the
/// parity of mapped chiral atoms whose neighbors changed.
fn check_and_correct_chirality_of_matching_atoms(
    reactant: &ROMol,
    reactant_atom_idx: u32,
    product: &mut RWMol,
    mapping: &mut AtomMapping,
) -> Result<()> {
    let reactant_atom = reactant.atom(reactant_atom_idx);
    for product_atom_idx in mapping.products_of(reactant_atom_idx) {
        let inversion_flag = product
            .atom(product_atom_idx)
            .get_prop(cp::molInversionFlag)
            .and_then(RDValue::as_int)
            .unwrap_or(0);
        if !is_specified_chirality(reactant_atom.chiral_tag()) || inversion_flag > 2 {
            continue;
        }
        let (rd, pd) = (
            reactant.atom_degree(reactant_atom_idx),
            product.atom_degree(product_atom_idx),
        );
        if rd < 3 || pd < 3 || rd.abs_diff(pd) > 1 {
            continue;
        }
        let mut n_unknown = 0;
        let mut r_order: Vec<i32> = reactant
            .atom_bond_indices(reactant_atom_idx)
            .map(|b| b as i32)
            .collect();
        let mut p_order: Vec<i32> = Vec::new();
        for nbr in product.atom_neighbors(product_atom_idx).collect::<Vec<_>>() {
            let rbond = mapping
                .prod_react
                .get(&nbr)
                .and_then(|&r| reactant.bond_idx_between(reactant_atom_idx, r));
            match rbond {
                None => {
                    n_unknown += 1;
                    if n_unknown > 1 {
                        break;
                    }
                    p_order.push(-1);
                }
                Some(b) => p_order.push(b as i32),
            }
        }
        if n_unknown == 1 {
            if rd == pd {
                let unmatched = reactant
                    .atom_bond_indices(reactant_atom_idx)
                    .map(|b| b as i32)
                    .find(|b| !p_order.contains(b));
                if let (Some(unmatched), Some(pos)) =
                    (unmatched, p_order.iter().position(|&x| x == -1))
                {
                    p_order[pos] = unmatched;
                }
                n_unknown = 0;
                if p_order.contains(&-1) {
                    return Err(RDKitError::invariant("extra unmapped atom"));
                }
            } else if pd > rd {
                if let Some(pos) = p_order.iter().position(|&x| x == -1) {
                    p_order.remove(pos);
                }
                n_unknown = 0;
                if p_order.contains(&-1) {
                    return Err(RDKitError::invariant("extra unmapped atom"));
                }
            }
        }
        if n_unknown == 0 {
            if rd > pd {
                // drop reactant bonds that are not in the product
                let mut i = 0;
                while i < r_order.len() && r_order.len() > p_order.len() {
                    if p_order.contains(&r_order[i]) {
                        i += 1;
                    } else {
                        r_order.remove(i);
                    }
                }
            }
            let atom = product.atom_mut(product_atom_idx);
            atom.set_chiral_tag(reactant_atom.chiral_tag());
            let n_swaps = count_swaps_to_interconvert(&r_order, &p_order)?;
            let mut invert = n_swaps % 2 == 1;
            if inversion_flag == 1 {
                invert = !invert;
            }
            if invert {
                atom.invert_chirality();
            }
        }
    }
    Ok(())
}

/// `checkAndCorrectChiralityOfProduct`: chiral atoms copied in from the
/// reactant (not part of the template).
fn check_and_correct_chirality_of_product(
    reactant: &ROMol,
    chiral_atoms: &[u32],
    product: &mut RWMol,
    mapping: &mut AtomMapping,
) -> Result<()> {
    for &ra in chiral_atoms {
        let react_atom = reactant.atom(ra);
        if react_atom.chiral_tag() == ChiralType::Unspecified {
            return Err(RDKitError::invariant("missing atom chirality."));
        }
        let react_degree = reactant.atom_degree(ra);
        for (i, product_atom_idx) in mapping.products_of(ra).into_iter().enumerate() {
            if react_atom.chiral_tag() != product.atom(product_atom_idx).chiral_tag() {
                return Err(RDKitError::invariant("invalid product chirality."));
            }
            if react_degree != product.atom_degree(product_atom_idx) {
                product
                    .atom_mut(product_atom_idx)
                    .set_chiral_tag(ChiralType::Unspecified);
            } else if matches!(
                react_atom.chiral_tag(),
                ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
            ) {
                let mut new_order: Vec<i32> = Vec::new();
                for rb in reactant.atom_bond_indices(ra).collect::<Vec<_>>() {
                    let other = reactant.bond(rb).other_atom_idx(ra)?;
                    let prods = mapping
                        .react_prod
                        .get(&other)
                        .ok_or_else(|| RDKitError::invariant("other atom from bond not mapped."))?;
                    let neighbor = *prods
                        .get(i)
                        .ok_or_else(|| RDKitError::invariant("other atom from bond not mapped."))?;
                    let pbond = product
                        .bond_idx_between(product_atom_idx, neighbor)
                        .ok_or_else(|| {
                            RDKitError::invariant("no matching bond found in product")
                        })?;
                    new_order.push(pbond as i32);
                }
                let n_swaps = product.perturbation_order(product_atom_idx, &new_order)?;
                if n_swaps % 2 == 1 {
                    product.atom_mut(product_atom_idx).invert_chirality();
                }
            }
        }
    }
    Ok(())
}

/// `copyEnhancedStereoGroups`: carries the reactant's stereo groups over to
/// the product atoms that kept their stereo.
fn copy_enhanced_stereo_groups(
    reactant: &ROMol,
    product: &mut RWMol,
    mapping: &AtomMapping,
) -> Result<()> {
    let mut new_groups: Vec<StereoGroup> = Vec::new();
    for sg in reactant.stereo_groups() {
        let mut atoms = Vec::new();
        for &ra in sg.atoms() {
            let Some(prods) = mapping.react_prod.get(&ra) else {
                continue;
            };
            for &pa in prods {
                let atom = product.atom(pa);
                if atom.chiral_tag() == ChiralType::Unspecified {
                    continue;
                }
                if atom
                    .get_prop(cp::molInversionFlag)
                    .and_then(RDValue::as_int)
                    == Some(4)
                {
                    continue;
                }
                atoms.push(pa);
            }
        }
        if !atoms.is_empty() {
            new_groups.push(StereoGroup::new(
                sg.group_type(),
                atoms,
                Vec::new(),
                sg.read_id(),
            )?);
        }
    }
    if !new_groups.is_empty() {
        new_groups.extend(product.stereo_groups().iter().cloned());
        product.set_stereo_groups(new_groups)?;
    }
    Ok(())
}

/// `generateProductConformers`
fn generate_product_conformers(
    product_conf: &mut Conformer,
    reactant: &ROMol,
    mapping: &AtomMapping,
) -> Result<()> {
    let Some(react_conf) = reactant.conformers().first() else {
        return Ok(());
    };
    if react_conf.is_3d() {
        product_conf.set_3d(true);
    }
    for (&ra, prods) in &mapping.react_prod {
        if prods.len() > 1 {
            rd_log_warning!(
                "reactant atom match more than one product atom, coordinates need to be revised\n"
            );
        }
        for &pa in prods {
            product_conf.set_atom_pos(pa, react_conf.atom_pos(ra)?)?;
        }
    }
    Ok(())
}

/// `addReactantAtomsAndBonds`: merges one matched reactant into the product.
#[allow(clippy::too_many_arguments)]
fn add_reactant_atoms_and_bonds(
    rxn: &ChemicalReaction,
    product: &mut RWMol,
    reactant: &ROMol,
    match_: &MatchVect,
    reactant_template: &ROMol,
    product_conf: Option<&mut Conformer>,
    reactant_id: u32,
) -> Result<()> {
    let mut mapping = atom_mappings(match_, reactant_template, product, reactant.num_atoms())?;
    let mut visited = vec![false; reactant.num_atoms() as usize];
    set_reactant_bond_properties_to_product(product, reactant, &mapping)?;

    let mut chiral_atoms_to_check: Vec<u32> = Vec::new();
    for &(_, reactant_atom_idx) in match_ {
        let reactant_atom_idx = reactant_atom_idx as u32;
        if !mapping.mapped_atoms[reactant_atom_idx as usize] {
            continue;
        }
        let reactant_atom = reactant.atom(reactant_atom_idx);
        for product_atom_idx in mapping.products_of(reactant_atom_idx) {
            set_reactant_atom_properties_to_product(
                product.atom_mut(product_atom_idx),
                reactant_atom,
                rxn.implicit_properties_flag(),
                reactant_id,
            );
            if let Some(q) = reactant_atom.query() {
                let mut new_atom = product.atom(product_atom_idx).copy_detached();
                new_atom.set_query(Some(q.clone()));
                product.replace_atom(product_atom_idx, &new_atom, true)?;
            }
        }
        add_reactant_neighbors_to_product(
            reactant,
            reactant_atom_idx,
            product,
            &mut visited,
            &mut chiral_atoms_to_check,
            &mut mapping,
            reactant_id,
        )?;
        check_and_correct_chirality_of_matching_atoms(
            reactant,
            reactant_atom_idx,
            product,
            &mut mapping,
        )?;
    }
    check_and_correct_chirality_of_product(
        reactant,
        &chiral_atoms_to_check,
        product,
        &mut mapping,
    )?;
    update_stereo_bonds(product, reactant, &mapping)?;
    copy_enhanced_stereo_groups(reactant, product, &mapping)?;
    if let Some(conf) = product_conf {
        conf.resize(product.num_atoms() as usize);
        generate_product_conformers(conf, reactant, &mapping)?;
    }
    Ok(())
}

/// `copyTemplateStereoGroupsToMol`: stereo groups of the product template
/// (all atoms mapped) override those inherited from the reactants.
fn copy_template_stereo_groups_to_mol(template: &ROMol, product: &mut RWMol) -> Result<()> {
    if template.stereo_groups().is_empty() {
        return Ok(());
    }
    let mut in_template_groups = vec![false; product.num_atoms() as usize];
    let mut new_groups: Vec<StereoGroup> = Vec::new();
    for sg in template.stereo_groups() {
        let mut keep = true;
        let mut atoms = Vec::new();
        for &ta in sg.atoms() {
            let map_num = template.atom(ta).atom_map_num();
            if map_num == 0 {
                keep = false;
                break;
            }
            for pa in 0..product.num_atoms() {
                if product
                    .atom(pa)
                    .get_prop(cp::reactionMapNum)
                    .and_then(RDValue::as_int)
                    == Some(map_num)
                {
                    atoms.push(pa);
                    in_template_groups[pa as usize] = true;
                }
            }
        }
        if keep && !atoms.is_empty() {
            new_groups.push(StereoGroup::new(
                sg.group_type(),
                atoms,
                Vec::new(),
                sg.read_id(),
            )?);
        }
    }
    if !new_groups.is_empty() {
        for sg in product.stereo_groups().to_vec() {
            let overlapping = sg
                .atoms()
                .iter()
                .filter(|&&a| in_template_groups[a as usize])
                .count();
            if overlapping == 0 {
                new_groups.push(sg);
            } else if overlapping < sg.atoms().len() {
                let atoms = sg
                    .atoms()
                    .iter()
                    .copied()
                    .filter(|&a| !in_template_groups[a as usize])
                    .collect();
                new_groups.push(StereoGroup::new(
                    sg.group_type(),
                    atoms,
                    Vec::new(),
                    sg.read_id(),
                )?);
            }
        }
        product.set_stereo_groups(new_groups)?;
    }
    Ok(())
}

/// `generateOneProductSet`: one product per product template for one
/// combination of reactant matches.
fn generate_one_product_set(
    rxn: &ChemicalReaction,
    reactants: &[&ROMol],
    matches: &[MatchVect],
) -> Result<Vec<ROMol>> {
    if reactants.len() != matches.len() {
        return Err(RDKitError::precondition("vector size mismatch"));
    }
    let do_confs = reactants.iter().any(|r| r.num_conformers() > 0);
    let do_bond_dirs = reactants.iter().any(|r| {
        r.bonds()
            .iter()
            .any(|b| b.bond_type() == BondType::Single && b.bond_dir() != BondDir::None)
    });
    let mut res = Vec::with_capacity(rxn.num_product_templates());
    for ptempl in rxn.products() {
        let mut product = convert_template_to_mol(ptempl)?;
        let mut conf = do_confs.then(|| {
            let mut c = Conformer::new(0);
            c.set_3d(false);
            c
        });
        for (reactant_id, rtempl) in rxn.reactants().iter().enumerate() {
            add_reactant_atoms_and_bonds(
                rxn,
                &mut product,
                reactants[reactant_id],
                &matches[reactant_id],
                rtempl,
                conf.as_mut(),
                reactant_id as u32,
            )?;
        }
        if let Some(c) = conf {
            product.add_conformer(c, true)?;
        }
        if do_bond_dirs {
            set_double_bond_neighbor_directions(&mut product, None)?;
        }
        if !ptempl.stereo_groups().is_empty() {
            copy_template_stereo_groups_to_mol(ptempl, &mut product)?;
        }
        product.update_property_cache(false)?;
        res.push(product);
    }
    Ok(res)
}

// ---------------------------------------------------------- entry points

/// `run_Reactants`: runs the reaction on `reactants` (one per reactant
/// template) and returns the products of every combination of template
/// matches, at most `max_products` of them (`0`: no limit). Combinations are
/// built on up to `num_threads` threads (`<= 0`: all but `|n|` hardware
/// threads); the result does not depend on the thread count.
pub fn run_reactants(
    rxn: &ChemicalReaction,
    reactants: &[ROMol],
    max_products: u32,
    num_threads: i32,
) -> Result<Vec<Vec<ROMol>>> {
    if !rxn.is_initialized() {
        return Err(reaction_error(
            "initMatchers() must be called before runReactants()",
        ));
    }
    if reactants.len() != rxn.num_reactant_templates() {
        return Err(reaction_error(
            "Number of reactants provided does not match number of reactant templates.",
        ));
    }
    if rxn.num_product_templates() == 0 {
        return Ok(Vec::new());
    }
    let refs: Vec<&ROMol> = reactants.iter().collect();
    let Some(by_reactant) = reactant_matches(&refs, rxn, max_products, None)? else {
        return Ok(Vec::new());
    };
    let combinations = reactant_combinations(&by_reactant, max_products);
    par_map(&combinations, num_threads, |combination| {
        generate_one_product_set(rxn, &refs, combination)
    })
    .into_iter()
    .collect()
}

/// `run_Reactant`: runs one molecule against the reactant template at
/// `reactant_idx` of a multi-reactant reaction (the other reactants are
/// empty placeholders); one product set per match.
pub fn run_reactant(
    rxn: &ChemicalReaction,
    reactant: &ROMol,
    reactant_idx: usize,
    num_threads: i32,
) -> Result<Vec<Vec<ROMol>>> {
    if !rxn.is_initialized() {
        return Err(reaction_error(
            "initMatchers() must be called before runReactants()",
        ));
    }
    if rxn.num_product_templates() == 0 {
        return Ok(Vec::new());
    }
    if reactant_idx >= rxn.num_reactant_templates() {
        return Err(RDKitError::precondition("reactantIdx out of bounds"));
    }
    let placeholder = ROMol::new();
    let refs: Vec<&ROMol> = (0..rxn.num_reactant_templates())
        .map(|i| {
            if i == reactant_idx {
                reactant
            } else {
                &placeholder
            }
        })
        .collect();
    let Some(mut by_reactant) = reactant_matches(&refs, rxn, 1000, Some(reactant_idx))? else {
        return Ok(Vec::new());
    };
    let matches = std::mem::take(&mut by_reactant[reactant_idx]);
    let per_product: Vec<Vec<MatchVect>> = matches
        .into_iter()
        .map(|m| {
            let mut all = vec![MatchVect::new(); rxn.num_reactant_templates()];
            all[reactant_idx] = m;
            all
        })
        .collect();
    par_map(&per_product, num_threads, |m| {
        generate_one_product_set(rxn, &refs, m)
    })
    .into_iter()
    .collect()
}

/// Map number -> atom index of the mapped atoms of a template.
fn mapped_atoms_of(templ: &ROMol) -> BTreeMap<i32, u32> {
    templ
        .atoms()
        .iter()
        .enumerate()
        .filter_map(|(i, a)| map_number(a).filter(|&n| n != 0).map(|n| (n, i as u32)))
        .collect()
}

/// `updateAtomsModifiedByReaction`
fn update_atoms_modified_by_reaction(
    reactant: &mut RWMol,
    reactant_template: &ROMol,
    product_template: &ROMol,
    product_atom_map: &BTreeMap<i32, u32>,
    reactant_product_map: &BTreeMap<i32, u32>,
    match_: &MatchVect,
) -> Result<bool> {
    let mut modified = false;
    for (&map_num, &r_template_idx) in reactant_product_map {
        let r_atom = reactant_template.atom(r_template_idx);
        let p_atom = product_template.atom(product_atom_map[&map_num]);
        let atom_idx = match_[r_template_idx as usize].1 as u32;
        if r_atom.atomic_num() != p_atom.atomic_num()
            && (p_atom.atomic_num() != 0 || !p_atom.has_query())
        {
            reactant
                .atom_mut(atom_idx)
                .set_atomic_num(p_atom.atomic_num() as i32);
            modified = true;
        }
        if update_props_from_implicit_props(p_atom, reactant.atom_mut(atom_idx)) {
            modified = true;
        }
        let Some(flag) = p_atom
            .get_prop(cp::molInversionFlag)
            .and_then(RDValue::as_int)
        else {
            continue;
        };
        let tag = reactant.atom(atom_idx).chiral_tag();
        match flag {
            0 | 2 => {}
            1 => {
                if is_specified_chirality(tag) {
                    reactant.atom_mut(atom_idx).invert_chirality();
                    modified = true;
                }
            }
            3 => {
                reactant
                    .atom_mut(atom_idx)
                    .set_chiral_tag(ChiralType::Unspecified);
                modified = true;
            }
            4 => {
                reactant
                    .atom_mut(atom_idx)
                    .set_chiral_tag(p_atom.chiral_tag());
                modified = true;
                let p_idx = product_atom_map[&map_num];
                let p_order: Vec<i32> = product_template
                    .atom_neighbors(p_idx)
                    .map(|n| product_template.atom(n).atom_map_num())
                    .filter(|&n| n != 0)
                    .collect();
                let a_order: Vec<i32> = reactant
                    .atom_neighbors(atom_idx)
                    .filter_map(|n| {
                        match_
                            .iter()
                            .find(|&&(_, ri)| ri as u32 == n)
                            .map(|&(ti, _)| reactant_template.atom(ti as u32).atom_map_num())
                    })
                    .filter(|&n| n != 0)
                    .collect();
                if p_order.len() == a_order.len()
                    && count_swaps_to_interconvert(&a_order, &p_order)? % 2 == 1
                {
                    reactant.atom_mut(atom_idx).invert_chirality();
                }
            }
            _ => rd_log_warning!(
                "unrecognized chiral inversion/retention flag on product atom ignored\n"
            ),
        }
    }
    Ok(modified)
}

/// `updateBondsModifiedByReaction`
fn update_bonds_modified_by_reaction(
    reactant: &mut RWMol,
    reactant_template: &ROMol,
    product_template: &ROMol,
    product_atom_map: &BTreeMap<i32, u32>,
    reactant_product_map: &BTreeMap<i32, u32>,
    match_: &MatchVect,
) -> Result<bool> {
    let mut modified = false;
    let reactant_atom = |template_idx: u32| match_[template_idx as usize].1 as u32;
    for (&map_num, &r_idx) in reactant_product_map {
        let p_idx = product_atom_map[&map_num];
        for nbr in product_template.atom_neighbors(p_idx) {
            let nbr_map = product_template.atom(nbr).atom_map_num();
            let Some(&nbr_r_idx) = (nbr_map != 0)
                .then(|| reactant_product_map.get(&nbr_map))
                .flatten()
            else {
                continue;
            };
            let p_bond = product_template
                .bond_idx_between(p_idx, nbr)
                .map(|b| product_template.bond(b))
                .ok_or_else(|| {
                    RDKitError::invariant("missing bond between known neighbors in product")
                })?;
            match reactant_template.bond_idx_between(r_idx, nbr_r_idx) {
                Some(rb) => {
                    let r_bond = reactant_template.bond(rb);
                    if p_bond.bond_type() != BondType::Unspecified
                        && p_bond.bond_type() != r_bond.bond_type()
                    {
                        let bond = reactant
                            .bond_idx_between(
                                reactant_atom(r_bond.begin_atom_idx()),
                                reactant_atom(r_bond.end_atom_idx()),
                            )
                            .ok_or_else(|| {
                                RDKitError::invariant(
                                    "missing bond between known neighbors in reactant",
                                )
                            })?;
                        reactant.bond_mut(bond).set_bond_type(p_bond.bond_type());
                        modified = true;
                    }
                }
                None => {
                    let bond =
                        reactant.bond_idx_between(reactant_atom(r_idx), reactant_atom(nbr_r_idx));
                    match bond {
                        None => {
                            let beg = reactant_atom(
                                reactant_product_map[&product_template
                                    .atom(p_bond.begin_atom_idx())
                                    .atom_map_num()],
                            );
                            let end = reactant_atom(
                                reactant_product_map
                                    [&product_template.atom(p_bond.end_atom_idx()).atom_map_num()],
                            );
                            add_bond_to_product(p_bond, reactant, beg, end)?;
                            modified = true;
                        }
                        Some(b) if reactant.bond(b).bond_type() != p_bond.bond_type() => {
                            reactant.bond_mut(b).set_bond_type(p_bond.bond_type());
                            modified = true;
                        }
                        Some(_) => {}
                    }
                }
            }
        }
        for nbr in reactant_template.atom_neighbors(r_idx) {
            let nbr_map = reactant_template.atom(nbr).atom_map_num();
            if nbr_map != 0
                && let Some(&p_nbr) = product_atom_map.get(&nbr_map)
                && product_template.bond_idx_between(p_idx, p_nbr).is_none()
            {
                reactant.remove_bond(reactant_atom(r_idx), reactant_atom(nbr))?;
                modified = true;
            }
        }
    }
    Ok(modified)
}

/// `run_Reactant` (in place): applies a single-reactant, single-product
/// reaction that adds no atoms directly to `reactant`; with
/// `remove_unmatched_atoms` atoms not reachable from the reacting atoms are
/// dropped. Returns whether the molecule was modified.
pub fn run_reactant_in_place(
    rxn: &ChemicalReaction,
    reactant: &mut RWMol,
    remove_unmatched_atoms: bool,
) -> Result<bool> {
    if rxn.num_reactant_templates() != 1 {
        return Err(RDKitError::precondition("only one reactant supported"));
    }
    if rxn.num_product_templates() != 1 {
        return Err(RDKitError::precondition("only one product supported"));
    }
    if !rxn.is_initialized() {
        return Err(reaction_error(
            "initMatchers() must be called before runReactants()",
        ));
    }
    let reactant_template = &rxn.reactants()[0];
    let product_template = &rxn.products()[0];
    let product_atom_map = mapped_atoms_of(product_template);
    let reactant_product_map: BTreeMap<i32, u32> = mapped_atoms_of(reactant_template)
        .into_iter()
        .filter(|(n, _)| product_atom_map.contains_key(n))
        .collect();
    if product_template
        .atoms()
        .iter()
        .any(|a| a.atom_map_num() == 0 || !reactant_product_map.contains_key(&a.atom_map_num()))
    {
        return Err(reaction_error(
            "single component reactions which add atoms in the product are not supported",
        ));
    }
    let matches =
        reactant_matches_to_template(reactant, reactant_template, 1, rxn.substruct_params())?;
    let Some(match_) = matches.first() else {
        return Ok(false);
    };

    let mut atoms_to_remove = vec![false; reactant.num_atoms() as usize];
    // template atoms with no counterpart in the product
    for (i, atom) in reactant_template.atoms().iter().enumerate() {
        let n = atom.atom_map_num();
        if n == 0 || !reactant_product_map.contains_key(&n) {
            atoms_to_remove[match_[i].1 as usize] = true;
        }
    }
    if remove_unmatched_atoms {
        traverse_to_find_atoms_to_remove(reactant, reactant_template, &mut atoms_to_remove, match_);
    }
    let mut modified = false;
    reactant.begin_batch_edit()?;
    modified |= update_atoms_modified_by_reaction(
        reactant,
        reactant_template,
        product_template,
        &product_atom_map,
        &reactant_product_map,
        match_,
    )?;
    modified |= update_bonds_modified_by_reaction(
        reactant,
        reactant_template,
        product_template,
        &product_atom_map,
        &reactant_product_map,
        match_,
    )?;
    if atoms_to_remove.iter().any(|&r| r) {
        modified = true;
        for (i, remove) in atoms_to_remove.iter().enumerate() {
            if *remove {
                reactant.remove_atom(i as u32)?;
            }
        }
    }
    reactant.commit_batch_edit()?;
    Ok(modified)
}

/// `traverseToFindAtomsToRemove`: also removes the reactant atoms that are
/// only reachable through atoms being removed.
fn traverse_to_find_atoms_to_remove(
    reactant: &ROMol,
    templ: &ROMol,
    atoms: &mut [bool],
    match_: &MatchVect,
) {
    let mut to_remove: Vec<bool> = atoms.iter().map(|&a| !a).collect();
    for &(_, ri) in match_ {
        to_remove[ri as usize] = false;
    }
    for &(ti, ri) in match_ {
        let mut to_consider: VecDeque<u32> = VecDeque::new();
        if templ.atom(ti as u32).atom_map_num() != 0 && !atoms[ri as usize] {
            to_consider.push_back(ri as u32);
        }
        while let Some(atom) = to_consider.pop_back() {
            to_remove[atom as usize] = false;
            for nbr in reactant.atom_neighbors(atom) {
                if to_remove[nbr as usize] {
                    to_consider.push_front(nbr);
                }
            }
        }
    }
    for (a, r) in atoms.iter_mut().zip(to_remove) {
        *a |= r;
    }
}

/// `reduceProductToSideChains`: strips the scaffold atoms from a product,
/// leaving the side chains attached to dummy atoms (`add_dummy_atoms`) or
/// annotated with `_rgroupAtomMaps` / `_rgroupBonds`.
pub fn reduce_product_to_side_chains(product: &ROMol, add_dummy_atoms: bool) -> Result<ROMol> {
    let mut mol = product.clone();
    let num_atoms = mol.num_atoms();
    let bookmark_of = |mol: &ROMol, atom: u32| -> i32 {
        mol.atom_bookmarks
            .iter()
            .find(|(_, atoms)| atoms.contains(&atom))
            .map_or(-1, |(&mark, _)| mark)
    };
    let mut atoms_to_remove = Vec::new();
    for scaffold in (0..num_atoms).rev() {
        let a = mol.atom(scaffold);
        if !(a.has_prop(cp::reactionMapNum) || !a.has_prop(cp::reactantAtomIdx)) {
            continue;
        }
        // (atom, bond type, atom map number or -1)
        let mut bonds_to_product: Vec<(u32, BondType, i32)> = Vec::new();
        for nbr in mol.atom_neighbors(scaffold).collect::<Vec<_>>() {
            let n = mol.atom(nbr);
            if !n.has_prop(cp::reactionMapNum) && n.has_prop(cp::reactantAtomIdx) {
                let bond_type = mol
                    .bond_idx_between(scaffold, nbr)
                    .map(|b| mol.bond(b).bond_type())
                    .ok_or_else(|| RDKitError::invariant("missing bond"))?;
                let mapno = if n.has_prop(WAS_DUMMY) {
                    n.get_prop(cp::reactionMapNum)
                        .and_then(RDValue::as_int)
                        .unwrap_or(-1)
                } else {
                    -1
                };
                bonds_to_product.push((nbr, bond_type, mapno));
            }
        }
        let mapno = if bonds_to_product.is_empty() {
            -1
        } else {
            bookmark_of(&mol, scaffold)
        };
        atoms_to_remove.push(scaffold);
        if bonds_to_product.is_empty() {
            continue;
        }
        if add_dummy_atoms {
            let idx = mol.add_atom(Atom::new(0), true);
            for &(r_atom, bond_type, bmap) in &bonds_to_product {
                mol.add_bond(idx, r_atom, bond_type)?;
                let atommapno = if bmap == -1 { mapno } else { bmap };
                if atommapno != 0 {
                    mol.atom_mut(idx)
                        .dict_mut()
                        .set_val(cp::molAtomMapNumber, atommapno);
                }
            }
        } else {
            for &(r_atom, bond_type, bmap) in &bonds_to_product {
                let atommapno = if bmap == -1 { mapno } else { bmap };
                if mapno != -1 {
                    let atom = mol.atom_mut(r_atom);
                    let mut rgroups = match atom.get_prop(cp::_rgroupAtomMaps) {
                        Some(RDValue::VecInt(v)) => v.clone(),
                        _ => Vec::new(),
                    };
                    let mut bonds = match atom.get_prop("_rgroupBonds") {
                        Some(RDValue::VecInt(v)) => v.clone(),
                        _ => Vec::new(),
                    };
                    rgroups.push(atommapno);
                    bonds.push(bond_type as i32);
                    atom.dict_mut().set_val(cp::_rgroupAtomMaps, rgroups);
                    atom.dict_mut().set_val("_rgroupBonds", bonds);
                }
            }
        }
    }
    mol.begin_batch_edit()?;
    for ai in atoms_to_remove {
        mol.remove_atom(ai)?;
    }
    mol.commit_batch_edit()?;
    Ok(mol)
}
