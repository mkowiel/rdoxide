//! SMARTS writer: `MolToSmarts`, `MolFragmentToSmarts`, `GetAtomSmarts`,
//! `GetBondSmarts`.
//!
//! Ported from RDKit `Code/GraphMol/SmilesParse/SmartsWrite.{h,cpp}`
//! (2026.09.1), including the CXSMARTS writers.
//!
//! Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::cx_write::{cx_fields, get_cx_extensions};
use super::write::{SmilesWriteParams, get_atom_smiles, in_organic_subset};
use crate::graphmol::atom::{Atom, ChiralType, HybridizationType};
use crate::graphmol::bond::{Bond, BondDir, BondType};
use crate::graphmol::canon::{AtomColor, FragmentOptions, MolStackElem, canonicalize_fragment};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::{AtomQuery, BondQuery, MASS_INTEGER_CONVERSION_FACTOR, QueryKind};
use crate::graphmol::ring_info::FindRingType;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

const HAS_AND: u32 = 0x1;
const HAS_LOWAND: u32 = 0x2;
const HAS_OR: u32 = 0x4;
const HAS_RECURSION: u32 = 0x8;

const QATOM_HAS_STEREO_SET: &str = "_qatomHasStereoSet";

/// The value of an `ATOM_EQUALS_QUERY` (the `dynamic_cast` succeeds for
/// equality queries and their subclasses).
fn eq_val<F: crate::graphmol::queries::QueryTarget>(
    q: &crate::graphmol::queries::Query<F>,
) -> Option<i32> {
    match q.kind {
        QueryKind::Equality
        | QueryKind::Less
        | QueryKind::LessEqual
        | QueryKind::Greater
        | QueryKind::GreaterEqual
        | QueryKind::AtomRing
        | QueryKind::HasProp(_)
        | QueryKind::HasPropWithValue { .. } => Some(q.val),
        _ => None,
    }
}

/// `std::string::find(s) > 0 && < length()`: present, but not at position 0.
fn contains_not_at_start(haystack: &str, needle: &str) -> bool {
    haystack.find(needle).is_some_and(|p| p > 0)
}

/// `_combineChildSmarts`
fn combine_child_smarts(
    cs1: &str,
    features1: u32,
    cs2: &str,
    features2: u32,
    descrip: &str,
    features: &mut u32,
) -> Result<String> {
    let mut res = String::new();
    if contains_not_at_start(descrip, "Or") {
        if (features1 & HAS_LOWAND != 0 && features1 & HAS_OR != 0)
            || (features2 & HAS_LOWAND != 0 && features2 & HAS_OR != 0)
        {
            return Err(RDKitError::Value(
                "This is a non-smartable query - OR above and below AND in the binary tree".into(),
            ));
        }
        res += cs1;
        if !(cs1.is_empty() || cs2.is_empty()) {
            res.push(',');
        }
        res += cs2;
        *features |= HAS_OR;
    } else if contains_not_at_start(descrip, "And") {
        let symb = if features1 & HAS_OR != 0 || features2 & HAS_OR != 0 {
            *features |= HAS_LOWAND;
            ";"
        } else {
            *features |= HAS_AND;
            "&"
        };
        res += cs1;
        if !(cs1.is_empty() || cs2.is_empty()) {
            res += symb;
        }
        res += cs2;
    } else {
        return Err(RDKitError::Value(format!(
            "Don't know how to combine using {descrip}"
        )));
    }
    *features |= features1;
    *features |= features2;
    Ok(res)
}

/// Context for writing one atom: the atom and whether its molecule asked
/// for isomeric output (`_doIsoSmiles`).
struct AtomCtx<'a> {
    atom: &'a mut Atom,
    iso: bool,
}

impl AtomCtx<'_> {
    fn has_symbol(&self) -> bool {
        self.atom.has_prop(cp::smilesSymbol)
    }

    /// The chirality suffix appended by `getAtomSmartsSimple` (sets
    /// `_qatomHasStereoSet` so it is written once per atom).
    fn chirality_suffix(&mut self, need_paren: &mut bool) -> String {
        if self.iso
            && self.atom.chiral_tag() != ChiralType::Unspecified
            && !self.atom.has_prop(QATOM_HAS_STEREO_SET)
            && !self.atom.has_prop(cp::_brokenChirality)
        {
            self.atom.props.set_val(QATOM_HAS_STEREO_SET, 1i32);
            match self.atom.chiral_tag() {
                ChiralType::TetrahedralCW => {
                    *need_paren = true;
                    return "@@".into();
                }
                ChiralType::TetrahedralCCW => {
                    *need_paren = true;
                    return "@".into();
                }
                _ => {}
            }
        }
        String::new()
    }
}

/// `getAtomSmartsSimple`
fn atom_smarts_simple(
    ctx: &mut AtomCtx<'_>,
    query: &AtomQuery,
    need_paren: &mut bool,
    check_for_symbol: bool,
) -> Result<String> {
    let equery = eq_val(query);
    let mut descrip = query.description.as_str();
    #[derive(PartialEq)]
    enum Mods {
        None,
        Range,
        Less,
        Greater,
    }
    let mut mods = Mods::None;
    if let Some(d) = descrip.strip_prefix("range_") {
        mods = Mods::Range;
        descrip = d;
    } else if let Some(d) = descrip.strip_prefix("less_") {
        mods = Mods::Less;
        descrip = d;
    } else if let Some(d) = descrip.strip_prefix("greater_") {
        mods = Mods::Greater;
        descrip = d;
    }
    let mut res = String::new();
    let mut has_val = false;
    let mut simple = |s: &str, val: bool, paren: bool, res: &mut String, has_val: &mut bool| {
        res.push_str(s);
        *has_val = val;
        *need_paren = paren;
    };
    match descrip {
        "AtomImplicitHCount" => simple("h", true, true, &mut res, &mut has_val),
        "AtomHasImplicitH" => simple("h", false, true, &mut res, &mut has_val),
        "AtomTotalValence" => simple("v", true, true, &mut res, &mut has_val),
        "AtomAtomicNum" => {
            if !ctx.has_symbol() {
                simple("#", true, true, &mut res, &mut has_val);
            }
        }
        "AtomExplicitDegree" => simple("D", true, true, &mut res, &mut has_val),
        "AtomNonHydrogenDegree" => simple("d", true, true, &mut res, &mut has_val),
        "AtomTotalDegree" => simple("X", true, true, &mut res, &mut has_val),
        "AtomHasRingBond" => simple("x", false, true, &mut res, &mut has_val),
        "AtomHCount" => simple("H", true, true, &mut res, &mut has_val),
        "AtomIsAliphatic" => simple("A", false, false, &mut res, &mut has_val),
        "AtomIsAromatic" => simple("a", false, false, &mut res, &mut has_val),
        "AtomNull" => simple("*", false, false, &mut res, &mut has_val),
        "AtomInRing" => simple("R", false, true, &mut res, &mut has_val),
        "AtomMinRingSize" => simple("r", true, true, &mut res, &mut has_val),
        "AtomRingSize" => simple("k", true, true, &mut res, &mut has_val),
        "AtomInNRings" => {
            let hv = mods == Mods::None && equery.is_some_and(|v| v >= 0);
            simple("R", hv, true, &mut res, &mut has_val);
        }
        "AtomHasHeteroatomNeighbors" => simple("z", false, true, &mut res, &mut has_val),
        "AtomNumHeteroatomNeighbors" => simple("z", true, true, &mut res, &mut has_val),
        "AtomHasAliphaticHeteroatomNeighbors" => simple("Z", false, true, &mut res, &mut has_val),
        "AtomNumAliphaticHeteroatomNeighbors" => simple("Z", true, true, &mut res, &mut has_val),
        "AtomFormalCharge" | "AtomNegativeFormalCharge" => {
            let val = equery.unwrap_or(0);
            let neg_first = descrip == "AtomFormalCharge";
            res.push(if (val < 0) == neg_first { '-' } else { '+' });
            if mods == Mods::None && val.abs() != 1 {
                res += &val.abs().to_string();
            }
            *need_paren = true;
        }
        "AtomHybridization" if equery.is_some() => {
            res.push('^');
            let v = equery.unwrap_or(0);
            let h = [
                (HybridizationType::S, "0"),
                (HybridizationType::SP, "1"),
                (HybridizationType::SP2, "2"),
                (HybridizationType::SP3, "3"),
                (HybridizationType::SP3D, "4"),
                (HybridizationType::SP3D2, "5"),
            ];
            if let Some((_, s)) = h.iter().find(|(t, _)| *t as i32 == v) {
                res += s;
            }
            *need_paren = true;
        }
        "AtomMass" if equery.is_some() => {
            res += &format!("{}*", equery.unwrap_or(0) / MASS_INTEGER_CONVERSION_FACTOR);
            *need_paren = true;
        }
        "AtomIsotope" if equery.is_some() => {
            res += &format!("{}*", equery.unwrap_or(0));
            *need_paren = true;
        }
        "AtomRingBondCount" => simple("x", true, true, &mut res, &mut has_val),
        "AtomUnsaturated" => simple("$(*=,:,#*)", false, true, &mut res, &mut has_val),
        "AtomType" if equery.is_some() => {
            let v = equery.unwrap_or(0);
            let (at_num, aromatic) = if v > 1000 {
                (v - 1000, true)
            } else {
                (v, false)
            };
            if !check_for_symbol || !ctx.has_symbol() {
                let mut symbol = PeriodicTable::get()
                    .get_element_symbol(at_num as u32)?
                    .to_string();
                if aromatic {
                    let mut b = symbol.into_bytes();
                    b[0] = b[0].to_ascii_lowercase();
                    symbol = String::from_utf8(b).unwrap_or_default();
                }
                res += &symbol;
                if !in_organic_subset(at_num) {
                    *need_paren = true;
                }
            } else {
                res.push(if aromatic { 'a' } else { 'A' });
            }
        }
        _ => {
            crate::rd_log_warning!("Cannot write SMARTS for query type : {descrip}. Ignoring it.");
            res.push('*');
        }
    }
    if mods != Mods::None {
        res.push('{');
        match mods {
            Mods::Less => res += &format!("{}-", equery.unwrap_or(0)),
            Mods::Range => {
                let QueryKind::Range { lower, upper, .. } = query.kind else {
                    return Err(RDKitError::invariant(
                        "query could not be converted to range query",
                    ));
                };
                res += &format!("{lower}-{upper}");
            }
            Mods::Greater => res += &format!("-{}", equery.unwrap_or(0)),
            Mods::None => {}
        }
        res.push('}');
    } else if has_val {
        res += &equery.unwrap_or(0).to_string();
    }
    res += &ctx.chirality_suffix(need_paren);
    Ok(res)
}

/// `getRecursiveStructureQuerySmarts`
fn recursive_structure_query_smarts(
    query: &AtomQuery,
    params: &SmilesWriteParams,
) -> Result<String> {
    let QueryKind::Recursive { mol, .. } = &query.kind else {
        return Err(RDKitError::precondition("bad query"));
    };
    let (s, _, _) = mol_to_smarts_impl(&mol.0, params, None, None)?;
    let mut res = format!("$({s})");
    if query.negate {
        res = format!("!{res}");
    }
    Ok(res)
}

/// `_recurseGetSmarts`
fn recurse_get_smarts(
    ctx: &mut AtomCtx<'_>,
    node: &AtomQuery,
    negate: bool,
    features: &mut u32,
    params: &SmilesWriteParams,
) -> Result<String> {
    let mut descrip = node.description.clone();
    let mut child1_features = 0u32;
    if node.children.len() < 2 {
        return Err(RDKitError::invariant("Not enough children on the query"));
    }
    let child1 = &node.children[0];
    let dsc1 = child1.description.as_str();
    let mut need_paren = false;
    let csmarts1 = if dsc1 == "RecursiveStructure" {
        *features |= HAS_RECURSION;
        recursive_structure_query_smarts(child1, params)?
    } else if dsc1 != "AtomOr" && dsc1 != "AtomAnd" {
        let s = atom_smarts_simple(ctx, child1, &mut need_paren, descrip == "AtomAnd")?;
        if negate ^ child1.negate {
            format!("!{s}")
        } else {
            s
        }
    } else {
        recurse_get_smarts(
            ctx,
            child1,
            negate ^ child1.negate,
            &mut child1_features,
            params,
        )?
    };
    if negate {
        if descrip == "AtomOr" {
            descrip = "AtomAnd".into();
        } else if descrip == "AtomAnd" {
            descrip = "AtomOr".into();
        }
    }
    let mut res = csmarts1;
    for child2 in &node.children[1..] {
        let mut child2_features = 0u32;
        let dsc2 = child2.description.as_str();
        let csmarts2 = if dsc2 == "RecursiveStructure" {
            *features |= HAS_RECURSION;
            recursive_structure_query_smarts(child2, params)?
        } else if dsc2 != "AtomOr" && dsc2 != "AtomAnd" {
            let s = atom_smarts_simple(ctx, child2, &mut need_paren, false)?;
            if negate ^ child2.negate {
                format!("!{s}")
            } else {
                s
            }
        } else {
            recurse_get_smarts(
                ctx,
                child2,
                negate ^ child2.negate,
                &mut child2_features,
                params,
            )?
        };
        res = combine_child_smarts(
            &res,
            child1_features,
            &csmarts2,
            child2_features,
            &descrip,
            features,
        )?;
    }
    Ok(res)
}

/// `getBasicBondRepr`
fn basic_bond_repr(
    typ: BondType,
    dir: BondDir,
    reverse_dative: bool,
    params: &SmilesWriteParams,
) -> String {
    let dir_or = |default: &str| -> String {
        if params.do_isomeric_smiles {
            match dir {
                BondDir::EndDownRight => return "\\".into(),
                BondDir::EndUpRight => return "/".into(),
                _ => {}
            }
        }
        default.into()
    };
    match typ {
        BondType::Single => dir_or("-"),
        BondType::Double => "=".into(),
        BondType::Triple => "#".into(),
        BondType::Quadruple => "$".into(),
        BondType::Aromatic => dir_or(":"),
        BondType::Dative => {
            if params.include_dative_bonds {
                if reverse_dative {
                    "<-".into()
                } else {
                    "->".into()
                }
            } else {
                "-".into()
            }
        }
        BondType::Zero => "~".into(),
        _ => String::new(),
    }
}

/// `getBondSmartsSimple`
fn bond_smarts_simple(
    bond: &Bond,
    q: &BondQuery,
    atom_to_left: i32,
    params: &SmilesWriteParams,
) -> Result<String> {
    let equery = eq_val(q);
    Ok(match q.description.as_str() {
        "BondNull" => "~".into(),
        "BondInRing" => "@".into(),
        "SingleOrAromaticBond" => match bond.bond_dir() {
            BondDir::EndDownRight => "\\".into(),
            BondDir::EndUpRight => "/".into(),
            _ => String::new(),
        },
        "SingleOrDoubleBond" => "-,=".into(),
        "DoubleOrAromaticBond" => "=,:".into(),
        "SingleOrDoubleOrAromaticBond" => "-,=,:".into(),
        "BondDir" if equery.is_some() => {
            let v = equery.unwrap_or(0);
            if v == BondDir::EndDownRight as i32 {
                "\\".into()
            } else if v == BondDir::EndUpRight as i32 {
                "/".into()
            } else {
                return Err(RDKitError::Value(
                    "Can't write smarts for this bond dir type".into(),
                ));
            }
        }
        "BondOrder" if equery.is_some() => {
            let reverse = atom_to_left >= 0 && bond.begin_atom_idx() != atom_to_left as u32;
            let bt = BondType::from_u8(equery.unwrap_or(0) as u8).unwrap_or(BondType::Unspecified);
            basic_bond_repr(bt, bond.bond_dir(), reverse, params)
        }
        d => {
            return Err(RDKitError::Value(format!(
                "Can't write smarts for this query bond type: {d}"
            )));
        }
    })
}

/// `_recurseBondSmarts`
fn recurse_bond_smarts(
    bond: &Bond,
    node: &BondQuery,
    negate: bool,
    atom_to_left: i32,
    features: &mut u32,
    params: &SmilesWriteParams,
) -> Result<String> {
    let mut descrip = node.description.clone();
    if node.children.len() != 2 {
        return Err(RDKitError::invariant("Too many children on the query"));
    }
    let (child1, child2) = (&node.children[0], &node.children[1]);
    let (mut f1, mut f2) = (0u32, 0u32);
    let is_comp = |d: &str| d == "BondOr" || d == "BondAnd";
    let mut csmarts1;
    let mut csmarts2 = String::new();
    if !is_comp(&child1.description) {
        csmarts1 = bond_smarts_simple(bond, child1, atom_to_left, params)?;
        if negate ^ child1.negate {
            csmarts1 = format!("!{csmarts1}");
        }
    } else {
        csmarts1 = recurse_bond_smarts(
            bond,
            child1,
            negate ^ child1.negate,
            atom_to_left,
            &mut f1,
            params,
        )?;
    }
    if !is_comp(&child2.description) {
        csmarts2 = bond_smarts_simple(bond, child2, atom_to_left, params)?;
        if negate ^ child2.negate {
            csmarts2 = format!("!{csmarts2}");
        }
    } else {
        // note: upstream assigns this to csmarts1
        csmarts1 = recurse_bond_smarts(
            bond,
            child2,
            negate ^ child2.negate,
            atom_to_left,
            &mut f2,
            params,
        )?;
    }
    if negate {
        if descrip == "BondOr" {
            descrip = "BondAnd".into();
        } else if descrip == "BondAnd" {
            descrip = "BondOr".into();
        }
    }
    combine_child_smarts(&csmarts1, f1, &csmarts2, f2, &descrip, features)
}

/// `getNonQueryAtomSmarts`
fn non_query_atom_smarts(ctx: &mut AtomCtx<'_>) -> Result<String> {
    let atom = &*ctx.atom;
    let mut res = String::from("[");
    if atom.isotope() != 0 {
        res += &atom.isotope().to_string();
    }
    if let Some(RDValue::String(s)) = atom.get_prop(cp::smilesSymbol) {
        res += s;
    } else if in_organic_subset(atom.atomic_num() as i32) {
        res += &format!("#{}", atom.atomic_num());
    } else {
        res += &atom.symbol();
    }
    let mut added_chirality = false;
    if ctx.iso
        && atom.chiral_tag() != ChiralType::Unspecified
        && !atom.has_prop(QATOM_HAS_STEREO_SET)
        && !atom.has_prop(cp::_brokenChirality)
    {
        let tag = atom.chiral_tag();
        ctx.atom.props.set_val(QATOM_HAS_STEREO_SET, 1i32);
        match tag {
            ChiralType::TetrahedralCW => {
                res += "@@";
                added_chirality = true;
            }
            ChiralType::TetrahedralCCW => {
                res += "@";
                added_chirality = true;
            }
            _ => {}
        }
    }
    let atom = &*ctx.atom;
    if added_chirality && atom.num_explicit_hs() == 1 {
        res.push('H');
    }
    let chg = atom.formal_charge();
    if chg != 0 {
        match chg {
            -1 => res.push('-'),
            1 => res.push('+'),
            c if c < 0 => res += &c.to_string(),
            c => res += &format!("+{c}"),
        }
    }
    if let Some(m) = atom
        .get_prop(cp::molAtomMapNumber)
        .and_then(RDValue::as_int)
    {
        res += &format!(":{m}");
    }
    res.push(']');
    Ok(res)
}

/// `getNonQueryBondSmarts`
fn non_query_bond_smarts(bond: &Bond, atom_to_left: i32, params: &SmilesWriteParams) -> String {
    if bond.is_aromatic() {
        if params.do_isomeric_smiles {
            match bond.bond_dir() {
                BondDir::EndDownRight => return "\\".into(),
                BondDir::EndUpRight => return "/".into(),
                _ => {}
            }
        }
        ":".into()
    } else {
        let reverse = atom_to_left >= 0 && bond.begin_atom_idx() != atom_to_left as u32;
        basic_bond_repr(bond.bond_type(), bond.bond_dir(), reverse, params)
    }
}

/// `SmartsWrite::GetAtomSmarts` for the atom at `idx` of `mol`.
pub fn get_atom_smarts(mol: &mut ROMol, idx: u32, params: &SmilesWriteParams) -> Result<String> {
    let iso = mol.has_prop(cp::_doIsoSmiles);
    if mol
        .atom(idx)
        .query()
        .is_some_and(|q| q.description.is_empty())
    {
        return get_atom_smiles(mol, idx, &SmilesWriteParams::default());
    }
    let atom = mol.atom_mut(idx);
    get_atom_smarts_for(atom, iso, params)
}

/// `SmartsWrite::GetAtomSmarts` for a free-standing atom (no owning
/// molecule).
pub fn get_atom_smarts_for(
    atom: &mut Atom,
    iso: bool,
    params: &SmilesWriteParams,
) -> Result<String> {
    let mut ctx = AtomCtx { atom, iso };
    let Some(query) = ctx.atom.query().cloned() else {
        return non_query_atom_smarts(&mut ctx);
    };
    let mut need_paren = false;
    let mut features = 0u32;
    let mut res;
    let descrip = query.description.as_str();
    if descrip == "AtomOr" || descrip == "AtomAnd" {
        need_paren = true;
        res = recurse_get_smarts(&mut ctx, &query, query.negate, &mut features, params)?;
        if res.len() == 1 {
            need_paren = false;
        }
    } else if descrip == "RecursiveStructure" {
        res = recursive_structure_query_smarts(&query, params)?;
        need_paren = true;
    } else {
        res = atom_smarts_simple(&mut ctx, &query, &mut need_paren, true)?;
        if query.negate {
            res = format!("!{res}");
            need_paren = true;
        }
    }
    if let Some(m) = ctx.atom.get_prop(cp::molAtomMapNumber) {
        need_paren = true;
        res = format!("{res}:{}", m.to_rd_string());
    }
    if let Some(RDValue::String(symbol)) = ctx.atom.get_prop(cp::smilesSymbol) {
        need_paren = true;
        res = if res.is_empty() {
            symbol.clone()
        } else {
            format!("{symbol};{res}")
        };
    }
    if need_paren {
        res = format!("[{res}]");
    }
    Ok(res)
}

/// `SmartsWrite::GetBondSmarts`
pub fn get_bond_smarts(
    bond: &Bond,
    params: &SmilesWriteParams,
    atom_to_left: i32,
) -> Result<String> {
    let Some(query) = bond.query() else {
        return Ok(non_query_bond_smarts(bond, atom_to_left, params));
    };
    let mut features = 0u32;
    let d = query.description.as_str();
    if d == "BondAnd" || d == "BondOr" {
        recurse_bond_smarts(
            bond,
            query,
            query.negate,
            atom_to_left,
            &mut features,
            params,
        )
    } else {
        let mut res = if query.negate {
            "!".to_string()
        } else {
            String::new()
        };
        res += &bond_smarts_simple(bond, query, atom_to_left, params)?;
        Ok(res)
    }
}

/// `FragmentSmartsConstruct`
#[allow(clippy::too_many_arguments)]
fn fragment_smarts_construct(
    mol: &mut ROMol,
    atom_idx: u32,
    colors: &mut [AtomColor],
    ranks: &[u32],
    params: &SmilesWriteParams,
    atom_ordering: &mut Vec<u32>,
    bond_ordering: &mut Vec<u32>,
    atoms_in_play: &[bool],
    bonds_in_play: Option<&[bool]>,
) -> Result<String> {
    mol.ring_info_mut().reset(true);
    mol.ring_info_mut().initialize(FindRingType::SymmSssr);
    for i in 0..mol.num_atoms() {
        mol.update_atom_property_cache(i, false)?;
    }
    let mut stack = Vec::with_capacity((mol.num_atoms() + mol.num_bonds()) as usize);
    canonicalize_fragment(
        mol,
        atom_idx,
        colors,
        ranks,
        &mut stack,
        FragmentOptions {
            atoms_in_play: Some(atoms_in_play),
            bonds_in_play,
            bond_symbols: None,
            do_isomeric_smiles: params.do_isomeric_smiles,
            do_random: false,
            do_chiral_inversions: true,
        },
    )?;
    mol.ring_info_mut().reset(true);
    let mut res = String::new();
    for elem in stack {
        match elem {
            MolStackElem::Atom(a) => {
                res += &get_atom_smarts(mol, a, params)?;
                atom_ordering.push(a);
            }
            MolStackElem::Bond(b, left) => {
                res += &get_bond_smarts(mol.bond(b), params, left)?;
                bond_ordering.push(b);
            }
            MolStackElem::Ring(n) => {
                if n < 10 {
                    res += &n.to_string();
                } else {
                    res += &format!("%{n}");
                }
            }
            MolStackElem::BranchOpen(_) => res.push('('),
            MolStackElem::BranchClose(_) => res.push(')'),
        }
    }
    Ok(res)
}

/// `molToSmarts`: returns the SMARTS and the atom/bond output orders.
fn mol_to_smarts_impl(
    inmol: &ROMol,
    params: &SmilesWriteParams,
    colors: Option<Vec<AtomColor>>,
    bonds_in_play: Option<&[bool]>,
) -> Result<(String, Vec<u32>, Vec<u32>)> {
    let n = inmol.num_atoms() as usize;
    if n == 0 {
        return Ok((String::new(), Vec::new(), Vec::new()));
    }
    if params.rooted_at_atom >= n as i32 {
        return Err(RDKitError::precondition("bad atom index"));
    }
    let mut colors = colors.unwrap_or_else(|| vec![AtomColor::White; n]);
    let atoms_in_play: Vec<bool> = colors.iter().map(|c| *c == AtomColor::White).collect();
    let mut mol = ROMol::copy_from(inmol, false, -1);
    let ranks: Vec<u32> = (0..n as u32).collect();
    if params.do_isomeric_smiles {
        mol.set_prop(cp::_doIsoSmiles, 1i32, false)?;
    }
    let mut atom_ordering = Vec::new();
    let mut bond_ordering = Vec::new();
    let mut res = String::new();
    while colors.contains(&AtomColor::White) {
        let next = if params.rooted_at_atom > -1
            && colors[params.rooted_at_atom as usize] == AtomColor::White
        {
            params.rooted_at_atom as u32
        } else {
            let mut next_rank = n as u32 + 1;
            let mut next = 0u32;
            for i in 0..n {
                if colors[i] != AtomColor::White {
                    continue;
                }
                let t = mol.atom(i as u32).chiral_tag();
                if t != ChiralType::TetrahedralCCW && t != ChiralType::TetrahedralCW {
                    next = i as u32;
                    break;
                }
                if ranks[i] < next_rank {
                    next_rank = ranks[i];
                    next = i as u32;
                }
            }
            next
        };
        res += &fragment_smarts_construct(
            &mut mol,
            next,
            &mut colors,
            &ranks,
            params,
            &mut atom_ordering,
            &mut bond_ordering,
            &atoms_in_play,
            bonds_in_play,
        )?;
        if colors.contains(&AtomColor::White) {
            res.push('.');
        }
    }
    Ok((res, atom_ordering, bond_ordering))
}

/// `MolToSmarts(mol, params)`; sets `_smilesAtomOutputOrder` /
/// `_smilesBondOutputOrder` on `mol`.
pub fn mol_to_smarts(mol: &mut ROMol, params: &SmilesWriteParams) -> Result<String> {
    if mol.num_atoms() == 0 {
        return Ok(String::new());
    }
    let (res, ao, bo) = mol_to_smarts_impl(mol, params, None, None)?;
    mol.set_prop(cp::_smilesAtomOutputOrder, ao, true)?;
    mol.set_prop(cp::_smilesBondOutputOrder, bo, true)?;
    Ok(res)
}

/// `MolFragmentToSmarts(mol, params, atomsToUse, bondsToUse)`
pub fn mol_fragment_to_smarts(
    mol: &mut ROMol,
    params: &SmilesWriteParams,
    atoms_to_use: &[i32],
    bonds_to_use: Option<&[i32]>,
) -> Result<String> {
    if atoms_to_use.is_empty() {
        return Err(RDKitError::precondition("no atoms provided"));
    }
    if bonds_to_use.is_some_and(|b| b.is_empty()) {
        return Err(RDKitError::precondition("no bonds provided"));
    }
    let n = mol.num_atoms() as usize;
    if n == 0 {
        return Ok(String::new());
    }
    let bonds_in_play = match bonds_to_use {
        Some(bs) => {
            let mut v = vec![false; mol.num_bonds() as usize];
            for &b in bs {
                *v.get_mut(b as usize)
                    .ok_or_else(|| RDKitError::range("idx", b, mol.num_bonds()))? = true;
            }
            Some(v)
        }
        None => None,
    };
    let mut colors = vec![AtomColor::Black; n];
    for &a in atoms_to_use {
        *colors
            .get_mut(a as usize)
            .ok_or_else(|| RDKitError::range("idx", a, n))? = AtomColor::White;
    }
    let ps = SmilesWriteParams {
        rooted_at_atom: -1,
        ..params.clone()
    };
    let (res, ao, bo) = mol_to_smarts_impl(mol, &ps, Some(colors), bonds_in_play.as_deref())?;
    mol.set_prop(cp::_smilesAtomOutputOrder, ao, true)?;
    mol.set_prop(cp::_smilesBondOutputOrder, bo, true)?;
    Ok(res)
}

/// Appends the CX extensions of `mol` (which must have just been written, so
/// that its output order is known) to a SMARTS.
fn with_cx_extensions(mol: &mut ROMol, mut smarts: String) -> Result<String> {
    if !smarts.is_empty() {
        let extensions = get_cx_extensions(mol, cx_fields::ALL)?;
        if !extensions.is_empty() {
            smarts.push(' ');
            smarts += &extensions;
        }
    }
    Ok(smarts)
}

/// `MolToCXSmarts`: the SMARTS of a molecule followed by its CX extensions.
/// Dative bonds are written as the single bonds they are in SMARTS.
pub fn mol_to_cx_smarts(mol: &mut ROMol, params: &SmilesWriteParams) -> Result<String> {
    let ps = SmilesWriteParams {
        include_dative_bonds: false,
        ..params.clone()
    };
    let smarts = mol_to_smarts(mol, &ps)?;
    with_cx_extensions(mol, smarts)
}

/// `MolFragmentToCXSmarts`
pub fn mol_fragment_to_cx_smarts(
    mol: &mut ROMol,
    params: &SmilesWriteParams,
    atoms_to_use: &[i32],
    bonds_to_use: Option<&[i32]>,
) -> Result<String> {
    let smarts = mol_fragment_to_smarts(mol, params, atoms_to_use, bonds_to_use)?;
    with_cx_extensions(mol, smarts)
}
