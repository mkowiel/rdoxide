//! CDXML parser.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/CDXMLParser.{h,cpp}` (the
//! build without the ChemDraw library: only the XML flavour is read, binary
//! CDX is not supported and molecules cannot be written).
//!
//! Copyright (c) 2022 Brian P Kelley. Licensed under the BSD 3-Clause license
//! (see LICENSES/RDKit-BSD-3-Clause.txt).

mod xml;

use std::collections::{BTreeMap, HashMap};

use self::xml::Node;
use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::chem_transforms::mol_zip::{MolzipLabel, MolzipParams, molzip_single};
use crate::graphmol::chirality::clear_single_bond_dir_flags;
use crate::graphmol::chirality::from_coords::{
    assign_chiral_types_from_bond_dirs, detect_bond_stereochemistry,
};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::find_rings::fast_find_rings;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::query_hs::{has_query_hs, merge_query_hs};
use crate::graphmol::mol_ops::sanitize_flags::{SANITIZE_ALL, SANITIZE_CLEANUP};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::queries::CompositeQueryType;
use crate::graphmol::queries::{
    AtomFn, AtomQuery, Query, QueryKind, make_a_atom_query, make_atom_explicit_degree_query,
    make_atom_implicit_h_count_query, make_atom_null_query, make_atom_num_query,
    make_atom_range_query, make_atom_ring_bond_count_query, make_atom_unsaturated_query,
    make_bond_is_in_ring_query, make_bond_null_query, make_bond_order_equals_query,
    make_double_or_aromatic_bond_query, make_m_atom_query, make_mh_atom_query, make_q_atom_query,
    make_single_or_aromatic_bond_query, make_single_or_double_bond_query, make_x_atom_query,
};
use crate::graphmol::query_ops::{
    atom_expand_query, bond_expand_query, query_atom_from_atom, replace_atom_with_query_atom,
};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

const NEEDS_FUSE: &str = "CDXML_NEEDS_FUSE";
const CDXML_FRAG_ID: &str = "CDXML_FRAG_ID";
const FUSE_LABEL: &str = "CDXML_NODE_ID";
const CDX_SCHEME_ID: &str = "CDX_SCHEME_ID";
const CDX_STEP_ID: &str = "CDX_STEP_ID";
const CDX_REAGENT_ID: &str = "CDX_REAGENT_ID";
const CDX_PRODUCT_ID: &str = "CDX_PRODUCT_ID";
const CDX_AGENT_ID: &str = "CDX_AGENT_ID";
const CDX_ATOM_POS: &str = "CDX_ATOM_POS";
const CDX_ATOM_ID: &str = "_CDX_ATOM_ID";
const CDX_BOND_ID: &str = "_CDX_BOND_ID";
const CDX_BOND_ORDERING: &str = "CDX_BOND_ORDERING";
const CDXML_FREE_SITES_PROP: &str = "_cdxmlFreeSites";
const CDXML_RING_BOND_COUNT_AS_DRAWN_PROP: &str = "_cdxmlRingBondCountAsDrawn";
const CDXML_LINK_NODE_MIN_REP_PROP: &str = "_cdxmlLinkNodeMinRep";
const CDXML_LINK_NODE_MAX_REP_PROP: &str = "_cdxmlLinkNodeMaxRep";
const CDXML_VARIABLE_ATTACHMENT_ENDPOINTS_PROP: &str = "_cdxmlVariableAttachmentEndpoints";
const RDKIT_DEPICT_BONDLENGTH: f64 = 1.5;

/// `CDXMLFormat`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CdxmlFormat {
    Cdxml = 0,
    Cdx = 1,
    Auto = 2,
}

/// `CDXMLParserParams`
#[derive(Debug, Clone)]
pub struct CdxmlParserParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub format: CdxmlFormat,
    pub parse_queries: bool,
    pub strict_query_parsing: bool,
}

impl Default for CdxmlParserParams {
    fn default() -> Self {
        CdxmlParserParams {
            sanitize: true,
            remove_hs: true,
            format: CdxmlFormat::Auto,
            parse_queries: false,
            strict_query_parsing: false,
        }
    }
}

/// `hasChemDrawCDXSupport`
pub fn has_chem_draw_cdx_support() -> bool {
    false
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BondQueryType {
    None,
    Any,
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Topology {
    None,
    Ring,
    Chain,
}

struct BondInfo {
    bond_id: i32,
    start: i32,
    end: i32,
    order: BondType,
    display: String,
    query_type: BondQueryType,
    topology: Topology,
}

#[derive(Default)]
struct SchemeInfo {
    scheme_id: i32,
    step_id: i32,
    products: Vec<u32>,
    reactants: Vec<u32>,
    above_arrow: Vec<u32>,
    below_arrow: Vec<u32>,
    atom_map: Vec<u32>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Unsaturation {
    None,
    MustBeAbsent,
    MustBePresent,
}

fn invariant(msg: &str) -> RDKitError {
    RDKitError::invariant(msg)
}

// ----- number parsing with the semantics of std::stoi / istream extraction -----

/// `std::stoi`: skips leading whitespace, reads an optionally signed integer
/// prefix and ignores whatever follows.
fn stoi(s: &str) -> Option<i32> {
    let t = s.trim_start();
    let bytes = t.as_bytes();
    let mut end = 0;
    if end < bytes.len() && (bytes[end] == b'-' || bytes[end] == b'+') {
        end += 1;
    }
    let digits = end;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end == digits {
        return None;
    }
    t[..end].parse().ok()
}

/// The `ptree::get<int>` conversion: the whole text (surrounding blanks
/// aside) must be an integer.
fn ptree_int(s: &str) -> Option<i32> {
    s.trim().parse().ok()
}

/// Repeated `istream >> int`: stops at the first token that does not start
/// with an integer; an integer prefix of a token counts.
fn to_vec_i32(s: &str) -> Vec<i32> {
    let mut out = Vec::new();
    let mut rest = s;
    loop {
        let t = rest.trim_start();
        if t.is_empty() {
            break;
        }
        let Some(v) = stoi(t) else { break };
        out.push(v);
        let b = t.as_bytes();
        let mut end = 0;
        if b[end] == b'-' || b[end] == b'+' {
            end += 1;
        }
        while end < b.len() && b[end].is_ascii_digit() {
            end += 1;
        }
        rest = &t[end..];
        // the remainder of a token such as "1.5" is not an integer
        if rest.chars().next().is_some_and(|c| !c.is_whitespace()) {
            break;
        }
    }
    out
}

fn to_vec_u32(s: &str) -> Vec<u32> {
    to_vec_i32(s).into_iter().map(|v| v as u32).collect()
}

fn to_vec_f64(s: &str) -> Vec<f64> {
    let mut out = Vec::new();
    for tok in s.split_whitespace() {
        match tok.parse::<f64>() {
            Ok(v) => out.push(v),
            Err(_) => break,
        }
    }
    out
}

// ----- molecule helpers -----

fn get_fuse_label(mol: &ROMol, idx: u32) -> u32 {
    match mol.atom(idx).get_prop(FUSE_LABEL) {
        Some(RDValue::UnsignedInt(u)) => *u,
        Some(RDValue::Int(i)) => *i as u32,
        _ => 0,
    }
}

fn set_fuse_label(atom: &mut Atom, idx: u32) -> Result<()> {
    if idx != 0 {
        atom.set_prop(FUSE_LABEL, idx, false)
    } else {
        atom.clear_prop(FUSE_LABEL);
        Ok(())
    }
}

fn prop_u32(atom: &Atom, key: &str) -> Option<u32> {
    match atom.get_prop(key)? {
        RDValue::UnsignedInt(u) => Some(*u),
        RDValue::Int(i) => Some(*i as u32),
        _ => None,
    }
}

fn prop_i32(atom: &Atom, key: &str) -> Option<i32> {
    match atom.get_prop(key)? {
        RDValue::Int(i) => Some(*i),
        RDValue::UnsignedInt(u) => Some(*u as i32),
        _ => None,
    }
}

/// `addquery`: gives the atom at `idx` the query `qry`, an optional label
/// and `noImplicit`.
fn add_query(mol: &mut RWMol, qry: AtomQuery, symbol: &str, idx: u32) -> Result<()> {
    let mut qa = query_atom_from_atom(mol.atom(idx))?;
    qa.set_query(Some(qry));
    qa.set_no_implicit(true);
    mol.replace_atom(idx, &qa, false)?;
    if !symbol.is_empty() {
        mol.atom_mut(idx).set_prop(cp::atomLabel, symbol, false)?;
    }
    Ok(())
}

fn to_query_atom(mol: &mut RWMol, idx: u32) -> Result<()> {
    replace_atom_with_query_atom(mol, idx)?;
    mol.atom_mut(idx).set_no_implicit(true);
    Ok(())
}

fn expand(mol: &mut RWMol, idx: u32, q: AtomQuery) -> Result<()> {
    atom_expand_query(mol.atom_mut(idx), q, CompositeQueryType::And, true)
}

fn is_r_query_label(label: &str) -> bool {
    label.starts_with('R')
}

fn is_any_non_hydrogen_query_label(label: &str) -> bool {
    label == "A" || label == "a"
}

fn is_wildcard_query_label(label: &str) -> bool {
    label == "*" || is_r_query_label(label)
}

#[allow(clippy::too_many_arguments)]
fn apply_atom_query_restrictions(
    mol: &mut RWMol,
    idx: u32,
    params: &CdxmlParserParams,
    restrict_implicit_hydrogens: bool,
    ring_bond_count: i32,
    ring_bond_count_at_least: bool,
    substituent_count: i32,
    max_substituent_count: i32,
    unsaturation: Unsaturation,
) -> Result<()> {
    if !params.parse_queries {
        return Ok(());
    }
    if !restrict_implicit_hydrogens
        && ring_bond_count < 0
        && !ring_bond_count_at_least
        && substituent_count < 0
        && max_substituent_count < 0
        && unsaturation == Unsaturation::None
    {
        return Ok(());
    }
    to_query_atom(mol, idx)?;
    if restrict_implicit_hydrogens {
        expand(mol, idx, make_atom_implicit_h_count_query(0))?;
    }
    if ring_bond_count >= 0 {
        expand(mol, idx, make_atom_ring_bond_count_query(ring_bond_count))?;
    } else if ring_bond_count_at_least {
        expand(
            mol,
            idx,
            Query::with(
                QueryKind::LessEqual,
                4,
                AtomFn::RingBondCount,
                "less_AtomRingBondCount",
            ),
        )?;
    }
    if substituent_count >= 0 {
        expand(mol, idx, make_atom_explicit_degree_query(substituent_count))?;
    }
    if max_substituent_count >= 0 {
        expand(
            mol,
            idx,
            make_atom_range_query(
                0,
                max_substituent_count,
                false,
                false,
                AtomFn::ExplicitDegree,
                "range_AtomExplicitDegree",
            ),
        )?;
    }
    if unsaturation != Unsaturation::None {
        let mut q = make_atom_unsaturated_query();
        if unsaturation == Unsaturation::MustBeAbsent {
            q.negate = true;
        }
        expand(mol, idx, q)?;
    }
    Ok(())
}

fn apply_deferred_free_sites(mol: &mut RWMol, params: &CdxmlParserParams) -> Result<()> {
    if !params.parse_queries {
        return Ok(());
    }
    for idx in 0..mol.num_atoms() {
        let Some(free_sites) = prop_i32(mol.atom(idx), CDXML_FREE_SITES_PROP) else {
            continue;
        };
        let min_degree = mol.atom_degree(idx) as i32;
        to_query_atom(mol, idx)?;
        expand(
            mol,
            idx,
            make_atom_range_query(
                min_degree,
                min_degree + free_sites,
                false,
                false,
                AtomFn::ExplicitDegree,
                "range_AtomExplicitDegree",
            ),
        )?;
        mol.atom_mut(idx).clear_prop(CDXML_FREE_SITES_PROP);
    }
    Ok(())
}

fn apply_deferred_ring_bond_count_as_drawn(
    mol: &mut RWMol,
    params: &CdxmlParserParams,
) -> Result<()> {
    if !params.parse_queries {
        return Ok(());
    }
    if !(0..mol.num_atoms()).any(|i| mol.atom(i).has_prop(CDXML_RING_BOND_COUNT_AS_DRAWN_PROP)) {
        return Ok(());
    }
    fast_find_rings(mol)?;
    for idx in 0..mol.num_atoms() {
        if !mol.atom(idx).has_prop(CDXML_RING_BOND_COUNT_AS_DRAWN_PROP) {
            continue;
        }
        let mut ring_bond_count = 0;
        let bond_idxs: Vec<u32> = mol.atom_bonds(idx).map(Bond::idx).collect();
        for b in bond_idxs {
            if mol.ring_info().num_bond_rings(b)? > 0 {
                ring_bond_count += 1;
            }
        }
        to_query_atom(mol, idx)?;
        expand(mol, idx, make_atom_ring_bond_count_query(ring_bond_count))?;
        mol.atom_mut(idx)
            .clear_prop(CDXML_RING_BOND_COUNT_AS_DRAWN_PROP);
    }
    Ok(())
}

fn apply_deferred_link_nodes(mol: &mut RWMol) -> Result<()> {
    let mut link_nodes = String::new();
    for idx in 0..mol.num_atoms() {
        let (Some(min_rep), Some(max_rep)) = (
            prop_i32(mol.atom(idx), CDXML_LINK_NODE_MIN_REP_PROP),
            prop_i32(mol.atom(idx), CDXML_LINK_NODE_MAX_REP_PROP),
        ) else {
            continue;
        };
        mol.atom_mut(idx).clear_prop(CDXML_LINK_NODE_MIN_REP_PROP);
        mol.atom_mut(idx).clear_prop(CDXML_LINK_NODE_MAX_REP_PROP);
        if min_rep <= 0 || max_rep < min_rep {
            rd_log_warning!("Invalid LinkNode repeat range on atom {}", idx);
            continue;
        }
        let mut nbrs: Vec<u32> = mol.atom_neighbors(idx).collect();
        if nbrs.len() != 2 {
            rd_log_warning!(
                "Only LinkNodes with exactly two neighbors are supported on atom {}",
                idx
            );
            continue;
        }
        if nbrs[0] > nbrs[1] {
            nbrs.swap(0, 1);
        }
        if !link_nodes.is_empty() {
            link_nodes.push('|');
        }
        link_nodes += &format!(
            "{} {} 2 {} {} {} {}",
            min_rep,
            max_rep,
            idx + 1,
            nbrs[0] + 1,
            idx + 1,
            nbrs[1] + 1
        );
    }
    if !link_nodes.is_empty() {
        if let Some(RDValue::String(existing)) = mol.get_prop(cp::molFileLinkNodes)
            && !existing.is_empty()
        {
            link_nodes = format!("{existing}|{link_nodes}");
        }
        mol.set_prop(cp::molFileLinkNodes, link_nodes, false)?;
    }
    Ok(())
}

fn apply_deferred_variable_attachments(mol: &mut RWMol) -> Result<()> {
    let mut atom_id_to_idx: HashMap<u32, u32> = HashMap::new();
    for idx in 0..mol.num_atoms() {
        if let Some(id) = prop_u32(mol.atom(idx), CDX_ATOM_ID) {
            atom_id_to_idx.insert(id, idx);
        }
    }
    for idx in 0..mol.num_atoms() {
        let Some(RDValue::VecUnsignedInt(ids)) = mol
            .atom(idx)
            .get_prop(CDXML_VARIABLE_ATTACHMENT_ENDPOINTS_PROP)
            .cloned()
        else {
            continue;
        };
        mol.atom_mut(idx)
            .clear_prop(CDXML_VARIABLE_ATTACHMENT_ENDPOINTS_PROP);
        if ids.is_empty() {
            continue;
        }
        if mol.atom_degree(idx) != 1 {
            rd_log_warning!(
                "Only attachment-point nodes with a single substituent bond are supported on atom {}",
                idx
            );
            continue;
        }
        let bond_idx = mol.atom_bonds(idx).next().expect("one bond").idx();
        let mut end_points = format!("({}", ids.len());
        let mut missing = false;
        for id in &ids {
            match atom_id_to_idx.get(id) {
                Some(i) => end_points += &format!(" {}", i + 1),
                None => {
                    rd_log_warning!("Attachment endpoint {} not found in molecule", id);
                    missing = true;
                    break;
                }
            }
        }
        if missing {
            continue;
        }
        end_points.push(')');
        let bond = mol.bond_mut(bond_idx);
        bond.set_prop(cp::_MolFileBondEndPts, end_points, false)?;
        bond.set_prop(cp::_MolFileBondAttach, "ANY", false)?;
    }
    Ok(())
}

// ----- fragment fusing -----

#[derive(Default)]
struct FragmentReplacement {
    replacement_atom: Option<u32>,
    fragment_atoms: Vec<u32>,
}

impl FragmentReplacement {
    fn replace(&self, mol: &mut RWMol) -> Result<bool> {
        let Some(rep) = self.replacement_atom else {
            return Ok(true);
        };
        let Some(RDValue::VecInt(bond_ordering)) =
            mol.atom(rep).get_prop(CDX_BOND_ORDERING).cloned()
        else {
            return Ok(true);
        };
        let replacement_bonds: Vec<u32> = mol.atom_bonds(rep).map(Bond::idx).collect();
        for b in replacement_bonds {
            let cdx_bond_id = match mol.bond(b).get_prop(CDX_BOND_ID) {
                Some(RDValue::Int(i)) => Some(*i),
                Some(RDValue::UnsignedInt(u)) => Some(*u as i32),
                _ => None,
            };
            let Some(bond_id) = cdx_bond_id else {
                rd_log_warning!("bond missing internal CDX BOND id, can't attach fragment at bond");
                return Ok(false);
            };
            let Some(pos) = bond_ordering.iter().position(|&x| x == bond_id) else {
                return Ok(false);
            };
            if pos >= self.fragment_atoms.len() {
                rd_log_warning!(
                    "bond ordering and number of atoms in fragment mismatch, can't attach fragment at bond:{}",
                    bond_id
                );
                return Ok(false);
            }
            let xatom = self.fragment_atoms[pos];
            let (begin, end, bond_type) = {
                let bond = mol.bond(b);
                (bond.begin_atom_idx(), bond.end_atom_idx(), bond.bond_type())
            };
            let xbonds: Vec<(u32, u32)> = mol
                .atom_bonds(xatom)
                .map(|xb| (xb.idx(), xb.other_atom_idx(xatom).expect("other atom")))
                .collect();
            for (_, other) in xbonds {
                if begin == rep {
                    mol.add_bond(other, end, bond_type)?;
                } else {
                    mol.add_bond(begin, other, bond_type)?;
                }
            }
        }
        mol.remove_atom(rep)?;
        for &atom in &self.fragment_atoms {
            mol.remove_atom(atom)?;
        }
        Ok(true)
    }
}

fn replace_fragments(mol: &mut RWMol) -> Result<bool> {
    let mut replacements: BTreeMap<i32, FragmentReplacement> = BTreeMap::new();
    for idx in 0..mol.num_atoms() {
        let label = get_fuse_label(mol, idx) as i32;
        if label == 0 {
            continue;
        }
        let frag = replacements.entry(label).or_default();
        if mol.atom(idx).has_prop(CDX_BOND_ORDERING) {
            frag.replacement_atom = Some(idx);
        } else {
            frag.fragment_atoms.push(idx);
        }
    }
    let mut success = true;
    mol.begin_batch_edit()?;
    for replacement in replacements.values() {
        if !replacement.replace(mol)? {
            success = false;
            break;
        }
    }
    mol.commit_batch_edit()?;
    Ok(success)
}

fn clear_internal_cdx_props(mols: &mut [RWMol]) {
    for mol in mols {
        for idx in 0..mol.num_atoms() {
            let atom = mol.atom_mut(idx);
            atom.clear_prop(CDX_ATOM_ID);
            atom.clear_prop(CDX_BOND_ORDERING);
        }
        for idx in 0..mol.num_bonds() {
            mol.bond_mut(idx).clear_prop(CDX_BOND_ID);
        }
    }
}

fn scale_bonds(mol: &ROMol, conf: &mut Conformer, target: f64, bond_length: f64) {
    let avg = if bond_length < 0.0 {
        // no bond length for the document: scale the average bond length
        let mut sum = 0.0;
        for b in 0..mol.num_bonds() {
            let bond = mol.bond(b);
            let pos = conf.positions();
            sum +=
                (pos[bond.begin_atom_idx() as usize] - pos[bond.end_atom_idx() as usize]).length();
        }
        sum / f64::from(mol.num_bonds())
    } else {
        bond_length
    };
    if avg > 0.0 {
        let scale = target / avg;
        for p in conf.positions_mut() {
            p.x *= scale;
            p.y *= scale;
            p.z *= scale;
        }
    }
}

struct StereoGroupInfo {
    sgroup: i32,
    group_type: StereoGroupType,
    atoms: Vec<u32>,
}

fn stereo_group_order(t: StereoGroupType) -> u8 {
    match t {
        StereoGroupType::Absolute => 0,
        StereoGroupType::Or => 1,
        StereoGroupType::And => 2,
    }
}

/// Document-wide state of one parse.
struct State {
    /// CDXML atom id -> atom index in its fragment's molecule
    ids: HashMap<u32, u32>,
    missing_frag_id: i32,
}

fn attr_int(node: &Node, key: &str, default: i32) -> i32 {
    node.attr(key).and_then(ptree_int).unwrap_or(default)
}

fn bad_attribute(frag_id: i32, node: &str, attr: &str, value: &str) -> bool {
    rd_log_error!(
        "Failed to parse XML fragment {} node: {} attribute: {}: {}",
        frag_id,
        node,
        attr,
        value
    );
    false
}

/// `parse_fragment`: adds the atoms and bonds of a `<fragment>` to `mol`.
/// Returns whether the fragment is usable.
fn parse_fragment(
    mol: &mut RWMol,
    frag: &Node,
    state: &mut State,
    params: &CdxmlParserParams,
    external_attachment: i32,
) -> Result<bool> {
    let mut frag_id = attr_int(frag, "id", -1);
    if frag_id == -1 {
        rd_log_warning!("Invalid or missing fragment id from CDXML fragment, assigning new one...");
        frag_id = state.missing_frag_id;
        state.missing_frag_id -= 1;
    }
    mol.set_prop(CDXML_FRAG_ID, frag_id, false)?;
    let mut atom_id: i32 = -1;
    let mut bonds: Vec<BondInfo> = Vec::new();
    let mut sgroups: BTreeMap<(i32, u8), StereoGroupInfo> = BTreeMap::new();
    let mut skip_fragment = false;

    for node in &frag.children {
        if node.name == "n" {
            let mut elemno: i32 = 6;
            let mut num_hydrogens: i32 = 0;
            let mut charge: i32 = 0;
            let mut atommap: i32 = 0;
            let mut mergeparent: i32 = -1;
            let mut rgroup_num: i32 = -1;
            let mut isotope: i32 = 0;
            let mut sgroup: i32 = -1;
            let mut explicit_hs = false;
            let mut grouptype = StereoGroupType::Absolute;
            let mut query_label = String::new();
            let mut bond_ordering: Vec<i32> = Vec::new();
            let mut elementlist: Vec<i32> = Vec::new();
            let mut atom_coords: Vec<f64> = Vec::new();
            let mut nodetype = String::new();
            let mut generic_nickname = String::new();
            let mut negative_elementlist = false;
            let mut restrict_implicit_hydrogens = false;
            let mut ring_bond_count: i32 = -1;
            let mut ring_bond_count_at_least = false;
            let mut ring_bond_count_as_drawn = false;
            let mut substituent_count: i32 = -1;
            let mut max_substituent_count: i32 = -1;
            let mut free_sites: i32 = -1;
            let mut link_count_low: i32 = -1;
            let mut link_count_high: i32 = -1;
            let mut variable_attachment_ids: Vec<u32> = Vec::new();
            let mut restrict_rxn_change = false;
            let mut rxn_stereo: i32 = 0;
            let mut unsaturation = Unsaturation::None;

            if node.attrs.is_empty() {
                return Err(RDKitError::FileParse("No such node (<xmlattr>)".into()));
            }
            for (key, value) in &node.attrs {
                macro_rules! int {
                    () => {
                        match stoi(value) {
                            Some(v) => v,
                            None => return Ok(bad_attribute(frag_id, &node.name, key, value)),
                        }
                    };
                }
                match key.as_str() {
                    "id" => {
                        atom_id = int!();
                        if state.ids.contains_key(&(atom_id as u32)) {
                            rd_log_error!(
                                "Warning, duplicated atom id {} skipping fragment",
                                atom_id
                            );
                            skip_fragment = true;
                            break;
                        }
                    }
                    "Element" => elemno = int!(),
                    "NumHydrogens" => {
                        num_hydrogens = int!();
                        explicit_hs = true;
                    }
                    "Charge" => charge = int!(),
                    "Isotope" => isotope = int!(),
                    "NodeType" => {
                        nodetype = value.clone();
                        match nodetype.as_str() {
                            "Nickname" | "Fragment" => {
                                elemno = 0;
                                atommap = atom_id;
                            }
                            "ExternalConnectionPoint" => {
                                if external_attachment <= 0 {
                                    rd_log_error!(
                                        "External Connection Point is not set skipping fragment"
                                    );
                                    skip_fragment = true;
                                    break;
                                }
                                elemno = 0;
                                atommap = external_attachment;
                                mergeparent = external_attachment;
                            }
                            "GenericNickname" => {
                                'outer: for tnode in &node.children {
                                    if tnode.name == "t" {
                                        for snode in &tnode.children {
                                            if snode.name == "s" {
                                                let s = snode.data.as_str();
                                                if !s.is_empty() {
                                                    if is_r_query_label(s) {
                                                        if s.len() > 1 {
                                                            rgroup_num = match stoi(&s[1..]) {
                                                                Some(v) => v,
                                                                None => {
                                                                    return Ok(bad_attribute(
                                                                        frag_id, &node.name, key,
                                                                        value,
                                                                    ));
                                                                }
                                                            };
                                                        }
                                                        elemno = 0;
                                                        query_label = if generic_nickname.is_empty()
                                                        {
                                                            "R".to_string()
                                                        } else {
                                                            generic_nickname.clone()
                                                        };
                                                    } else if is_any_non_hydrogen_query_label(s)
                                                        || matches!(s, "Q" | "X" | "M" | "MH" | "*")
                                                    {
                                                        query_label = s.to_string();
                                                        elemno = 0;
                                                    }
                                                }
                                                break;
                                            }
                                        }
                                        break 'outer;
                                    }
                                }
                            }
                            "ElementList" => query_label = "ElementList".to_string(),
                            "LinkNode" => {
                                link_count_low = 1;
                                link_count_high = 1;
                            }
                            "VariableAttachment" | "MultiAttachment" if params.parse_queries => {
                                elemno = 0;
                            }
                            _ => {}
                        }
                    }
                    "GenericNickname" => generic_nickname = value.clone(),
                    "LinkCountLow" => link_count_low = int!(),
                    "LinkCountHigh" => link_count_high = int!(),
                    "Attachments"
                        if (nodetype == "VariableAttachment" || nodetype == "MultiAttachment")
                            && params.parse_queries =>
                    {
                        variable_attachment_ids = to_vec_u32(value);
                    }
                    "BondOrdering" => bond_ordering = to_vec_i32(value),
                    "ElementList" => {
                        let mut text = value.as_str();
                        if let Some(rest) = text.strip_prefix("NOT ") {
                            negative_elementlist = true;
                            text = rest;
                        }
                        elementlist = to_vec_i32(text);
                    }
                    "ImplicitHydrogens" => {
                        restrict_implicit_hydrogens = value == "yes" || value == "true";
                    }
                    "RingBondCount" => match value.as_str() {
                        "NoRingBonds" => ring_bond_count = 0,
                        "SimpleRing" => ring_bond_count = 2,
                        "Fusion" => ring_bond_count = 3,
                        "SpiroOrHigher" => ring_bond_count_at_least = true,
                        "AsDrawn" => ring_bond_count_as_drawn = true,
                        "Unspecified" => {}
                        other => rd_log_warning!(
                            "Unhandled RingBondCount query value {} ignoring",
                            other
                        ),
                    },
                    "UnsaturatedBonds" => match value.as_str() {
                        "MustBePresent" => unsaturation = Unsaturation::MustBePresent,
                        "MustBeAbsent" => unsaturation = Unsaturation::MustBeAbsent,
                        "Unspecified" => {}
                        other => rd_log_warning!(
                            "Unhandled UnsaturatedBonds query value {} ignoring",
                            other
                        ),
                    },
                    "SubstituentsExactly" => substituent_count = int!(),
                    "SubstituentsUpTo" => max_substituent_count = int!(),
                    "FreeSites" => free_sites = int!(),
                    "RxnChange" => {
                        restrict_rxn_change = value == "yes" || value == "true" || value == "1";
                    }
                    "RxnStereo" => match value.as_str() {
                        "Inversion" => rxn_stereo = 1,
                        "Retention" => rxn_stereo = 2,
                        "Unspecified" => {}
                        other => {
                            rd_log_warning!("Unhandled RxnStereo query value {} ignoring", other)
                        }
                    },
                    "p" => atom_coords = to_vec_f64(value),
                    "EnhancedStereoGroupNum" => sgroup = int!(),
                    "EnhancedStereoType" => match value.as_str() {
                        "And" => grouptype = StereoGroupType::And,
                        "Or" => grouptype = StereoGroupType::Or,
                        "Absolute" => grouptype = StereoGroupType::Absolute,
                        other => {
                            rd_log_warning!("Unhandled enhanced stereo type {} ignoring", other);
                        }
                    },
                    _ => {}
                }
            }
            if nodetype == "GenericNickname"
                && query_label.is_empty()
                && !generic_nickname.is_empty()
            {
                if is_r_query_label(&generic_nickname) {
                    if generic_nickname.len() > 1 {
                        rgroup_num = stoi(&generic_nickname[1..])
                            .ok_or_else(|| RDKitError::invariant("stoi: no conversion"))?;
                    }
                    query_label = generic_nickname.clone();
                    elemno = 0;
                } else if is_any_non_hydrogen_query_label(&generic_nickname)
                    || matches!(generic_nickname.as_str(), "Q" | "X" | "M" | "MH" | "*")
                {
                    query_label = generic_nickname.clone();
                    elemno = 0;
                }
            }

            // add the atom
            if atom_id == -1 {
                return Err(invariant("Uninitialized atom id in cdxml."));
            }
            let mut rd_atom = Atom::new(elemno as u32);
            rd_atom.set_formal_charge(charge);
            rd_atom.set_num_explicit_hs(num_hydrogens as u32);
            rd_atom.set_no_implicit(explicit_hs);
            rd_atom.set_isotope(isotope as u32);
            if rgroup_num >= 0 {
                rd_atom.set_atom_map_num(rgroup_num, false)?;
            }
            set_fuse_label(&mut rd_atom, atommap as u32)?;
            if !bond_ordering.is_empty() {
                rd_atom.set_prop(CDX_BOND_ORDERING, bond_ordering, false)?;
            }
            if mergeparent > 0 {
                rd_atom.set_prop("MergeParent", mergeparent, false)?;
            }
            rd_atom.set_prop(CDX_ATOM_POS, atom_coords, false)?;
            rd_atom.set_prop(CDX_ATOM_ID, atom_id as u32, false)?;
            let idx = mol.add_atom(rd_atom, true);
            if !query_label.is_empty() {
                if is_wildcard_query_label(&query_label) {
                    add_query(mol, make_atom_null_query(), &query_label, idx)?;
                } else if is_any_non_hydrogen_query_label(&query_label) {
                    add_query(mol, make_a_atom_query(), &query_label, idx)?;
                } else if query_label == "Q" {
                    add_query(mol, make_q_atom_query(), &query_label, idx)?;
                } else if query_label == "X" {
                    add_query(mol, make_x_atom_query(), &query_label, idx)?;
                } else if query_label == "M" {
                    add_query(mol, make_m_atom_query(), &query_label, idx)?;
                } else if query_label == "MH" {
                    add_query(mol, make_mh_atom_query(), &query_label, idx)?;
                } else if query_label == "ElementList" {
                    if elementlist.is_empty() {
                        rd_log_warning!("ElementList is empty, ignoring...");
                    } else {
                        let mut q = Query::new(QueryKind::Or);
                        q.description = "AtomOr".to_string();
                        for at_num in &elementlist {
                            q.children.push(make_atom_num_query(*at_num));
                        }
                        q.negate = negative_elementlist;
                        add_query(mol, q, &query_label, idx)?;
                        mol.atom_mut(idx).set_atomic_num(elementlist[0]);
                    }
                } else {
                    mol.atom_mut(idx)
                        .set_prop(cp::atomLabel, query_label.as_str(), false)?;
                }
            }
            apply_atom_query_restrictions(
                mol,
                idx,
                params,
                restrict_implicit_hydrogens,
                ring_bond_count,
                ring_bond_count_at_least,
                substituent_count,
                max_substituent_count,
                unsaturation,
            )?;
            {
                let atom = mol.atom_mut(idx);
                if ring_bond_count_as_drawn {
                    atom.set_prop(CDXML_RING_BOND_COUNT_AS_DRAWN_PROP, 1i32, false)?;
                }
                if free_sites >= 0 {
                    atom.set_prop(CDXML_FREE_SITES_PROP, free_sites, false)?;
                }
                if link_count_low > 0 && link_count_high >= link_count_low {
                    atom.set_prop(CDXML_LINK_NODE_MIN_REP_PROP, link_count_low, false)?;
                    atom.set_prop(CDXML_LINK_NODE_MAX_REP_PROP, link_count_high, false)?;
                }
                if !variable_attachment_ids.is_empty() {
                    atom.set_prop(
                        CDXML_VARIABLE_ATTACHMENT_ENDPOINTS_PROP,
                        variable_attachment_ids,
                        false,
                    )?;
                }
                if params.parse_queries && restrict_rxn_change {
                    atom.set_prop(cp::molRxnExactChange, 1i32, false)?;
                }
                if params.parse_queries && rxn_stereo != 0 {
                    atom.set_prop(cp::molInversionFlag, rxn_stereo, false)?;
                }
            }
            if sgroup != -1 {
                let entry = sgroups
                    .entry((sgroup, stereo_group_order(grouptype)))
                    .or_insert_with(|| StereoGroupInfo {
                        sgroup,
                        group_type: grouptype,
                        atoms: Vec::new(),
                    });
                entry.atoms.push(idx);
            }
            state.ids.insert(atom_id as u32, idx);
            if nodetype == "Nickname" || nodetype == "Fragment" {
                for fragment in &node.children {
                    if fragment.name == "fragment" {
                        if !parse_fragment(mol, fragment, state, params, atom_id)? {
                            skip_fragment = true;
                            break;
                        }
                        mol.set_prop(NEEDS_FUSE, true, false)?;
                        // parse_fragment set the id of the nested fragment
                        mol.set_prop(CDXML_FRAG_ID, frag_id, false)?;
                    }
                }
            }
        } else if node.name == "b" {
            let mut bond_id: i32 = -1;
            let mut start_atom: i32 = -1;
            let mut end_atom: i32 = -1;
            let mut order = BondType::Single;
            let mut display = String::new();
            let mut query_type = BondQueryType::None;
            let mut topology = Topology::None;
            if node.attrs.is_empty() {
                return Err(RDKitError::FileParse("No such node (<xmlattr>)".into()));
            }
            for (key, value) in &node.attrs {
                macro_rules! int {
                    () => {
                        match stoi(value) {
                            Some(v) => v,
                            None => return Ok(bad_attribute(frag_id, &node.name, key, value)),
                        }
                    };
                }
                match key.as_str() {
                    "id" => bond_id = int!(),
                    "B" => start_atom = int!(),
                    "E" => end_atom = int!(),
                    "Order" => {
                        let tokens: Vec<&str> = value.split_whitespace().collect();
                        let t = |i: usize| tokens[i];
                        let parsed = if tokens.len() == 1 && value == "1.5" {
                            order = BondType::Aromatic;
                            true
                        } else if tokens.len() == 1 && value == "any" {
                            order = BondType::Unspecified;
                            query_type = BondQueryType::Any;
                            true
                        } else if tokens.len() == 1 && value == "dative" {
                            order = BondType::Dative;
                            true
                        } else if tokens.len() == 1 && (value == "hydrogen" || value == "Hydrogen")
                        {
                            order = BondType::Hydrogen;
                            true
                        } else if tokens.len() == 2
                            && ((t(0) == "1" && t(1) == "2") || (t(0) == "2" && t(1) == "1"))
                        {
                            order = BondType::Single;
                            query_type = BondQueryType::SingleOrDouble;
                            true
                        } else if tokens.len() == 2
                            && ((t(0) == "1" && (t(1) == "1.5" || t(1) == "4"))
                                || (t(1) == "1" && (t(0) == "1.5" || t(0) == "4")))
                        {
                            order = BondType::Single;
                            query_type = BondQueryType::SingleOrAromatic;
                            true
                        } else if tokens.len() == 2
                            && ((t(0) == "2" && (t(1) == "1.5" || t(1) == "4"))
                                || (t(1) == "2" && (t(0) == "1.5" || t(0) == "4")))
                        {
                            order = BondType::Double;
                            query_type = BondQueryType::DoubleOrAromatic;
                            true
                        } else {
                            match stoi(value) {
                                Some(1) => {
                                    order = BondType::Single;
                                    true
                                }
                                Some(2) => {
                                    order = BondType::Double;
                                    true
                                }
                                Some(3) => {
                                    order = BondType::Triple;
                                    true
                                }
                                Some(4) => {
                                    order = BondType::Quadruple;
                                    true
                                }
                                _ => false,
                            }
                        };
                        if !parsed {
                            return Ok(bad_attribute(frag_id, &node.name, key, value));
                        }
                    }
                    "Topology" => match value.as_str() {
                        "Ring" => topology = Topology::Ring,
                        "Chain" => topology = Topology::Chain,
                        "RingOrChain" | "Unspecified" => {}
                        _ => return Ok(bad_attribute(frag_id, &node.name, key, value)),
                    },
                    "Display" => display = value.clone(),
                    _ => {}
                }
            }
            let info = BondInfo {
                bond_id,
                start: start_atom,
                end: end_atom,
                order,
                display,
                query_type,
                topology,
            };
            // BondInfo::validate
            let valid = match (
                state.ids.get(&(info.start as u32)),
                state.ids.get(&(info.end as u32)),
            ) {
                (Some(&s), Some(&e)) => s < mol.num_atoms() && e < mol.num_atoms() && s != e,
                _ => false,
            };
            if !valid {
                rd_log_error!("Bad bond in CDXML skipping fragment {}...", frag_id);
                skip_fragment = true;
                break;
            }
            bonds.push(info);
        }
    }

    if !skip_fragment {
        for bond in &mut bonds {
            let mut swap = false;
            if bond.display == "WedgeEnd" {
                swap = true;
                bond.display = "WedgeBegin".to_string();
            }
            if bond.display == "WedgedHashEnd" {
                swap = true;
                bond.display = "WedgedHashBegin".to_string();
            }
            let mut start_idx = state.ids[&(bond.start as u32)];
            let mut end_idx = state.ids[&(bond.end as u32)];
            if swap {
                std::mem::swap(&mut start_idx, &mut end_idx);
            }
            let mut qb: Option<Bond> = None;
            match bond.query_type {
                BondQueryType::Any => {
                    let mut b = Bond::new(bond.order);
                    b.set_query(Some(make_bond_null_query()));
                    qb = Some(b);
                }
                BondQueryType::SingleOrDouble => {
                    let mut b = Bond::new(bond.order);
                    b.set_query(Some(make_single_or_double_bond_query()));
                    qb = Some(b);
                }
                BondQueryType::SingleOrAromatic => {
                    let mut b = Bond::new(bond.order);
                    b.set_query(Some(make_single_or_aromatic_bond_query()));
                    qb = Some(b);
                }
                BondQueryType::DoubleOrAromatic => {
                    let mut b = Bond::new(bond.order);
                    b.set_query(Some(make_double_or_aromatic_bond_query()));
                    qb = Some(b);
                }
                BondQueryType::None => {}
            }
            if bond.topology != Topology::None {
                // `QueryBond(bondType)` starts with a bond order query
                let b = qb.get_or_insert_with(|| {
                    let mut b = Bond::new(bond.order);
                    b.set_query(Some(make_bond_order_equals_query(bond.order)));
                    b
                });
                let mut topology_query = make_bond_is_in_ring_query();
                if bond.topology == Topology::Chain {
                    topology_query.negate = true;
                }
                bond_expand_query(b, topology_query, CompositeQueryType::And, true)?;
            }
            let bond_idx = if let Some(mut b) = qb {
                b.set_begin_atom_idx_unchecked(start_idx);
                b.set_end_atom_idx_unchecked(end_idx);
                mol.add_bond_object(b)? - 1
            } else {
                mol.add_bond(start_idx, end_idx, bond.order)? - 1
            };
            if bond.order == BondType::Aromatic {
                mol.bond_mut(bond_idx).set_is_aromatic(true);
                mol.atom_mut(state.ids[&(bond.end as u32)])
                    .set_is_aromatic(true);
                mol.atom_mut(state.ids[&(bond.start as u32)])
                    .set_is_aromatic(true);
            }
            let bnd = mol.bond_mut(bond_idx);
            bnd.set_prop(CDX_BOND_ID, bond.bond_id, false)?;
            match bond.display.as_str() {
                "WedgeBegin" => {
                    bnd.set_bond_dir(BondDir::BeginWedge);
                    bnd.set_prop(cp::_MolFileBondCfg, 1i32, false)?;
                }
                "WedgedHashBegin" => {
                    bnd.set_bond_dir(BondDir::BeginDash);
                    bnd.set_prop(cp::_MolFileBondCfg, 3i32, false)?;
                }
                "Wavy" => match bond.order {
                    BondType::Single => {
                        bnd.set_bond_dir(BondDir::Unknown);
                        bnd.set_prop(cp::_MolFileBondCfg, 2i32, false)?;
                    }
                    BondType::Double => {
                        bnd.set_bond_dir(BondDir::EitherDouble);
                        bnd.set_stereo(BondStereo::StereoAny)?;
                    }
                    _ => rd_log_warning!(
                        "ignoring Wavy bond set on a non double bond id: {}",
                        bond.bond_id
                    ),
                },
                _ => {}
            }
        }
        apply_deferred_ring_bond_count_as_drawn(mol, params)?;
        apply_deferred_free_sites(mol, params)?;
        apply_deferred_link_nodes(mol)?;
        apply_deferred_variable_attachments(mol)?;
    }

    if !sgroups.is_empty() {
        let mut groups = Vec::new();
        for info in sgroups.values() {
            let g_id = if info.group_type != StereoGroupType::Absolute && info.sgroup > 0 {
                info.sgroup as u32
            } else {
                0
            };
            groups.push(StereoGroup::new(
                info.group_type,
                info.atoms.clone(),
                Vec::new(),
                g_id,
            )?);
        }
        mol.set_stereo_groups(groups)?;
    }
    Ok(!skip_fragment)
}

#[allow(clippy::too_many_arguments)]
fn set_reaction_data(
    kind: &str,
    prop: &str,
    scheme: &SchemeInfo,
    frag_ids: &[u32],
    fragments: &HashMap<u32, usize>,
    grouped_fragments: &BTreeMap<u32, Vec<i32>>,
    mols: &mut [RWMol],
) -> Result<()> {
    let mut reagent_idx = 0u32;
    for idx in frag_ids {
        let Some(group) = grouped_fragments.get(idx) else {
            rd_log_warning!(
                "CDXMLParser: Schema {} step {} {} reaction fragment {} not found in document.",
                scheme.scheme_id,
                scheme.step_id,
                kind,
                idx
            );
            continue;
        };
        for reaction_fragment_id in group {
            let Some(&mol_idx) = fragments.get(&(*reaction_fragment_id as u32)) else {
                rd_log_warning!(
                    "CDXMLParser: Schema {} step {} {} fragment {} not found in document.",
                    scheme.scheme_id,
                    scheme.step_id,
                    kind,
                    idx
                );
                continue;
            };
            let mol = &mut mols[mol_idx];
            mol.set_prop(CDX_SCHEME_ID, scheme.scheme_id, false)?;
            mol.set_prop(CDX_STEP_ID, scheme.step_id, false)?;
            mol.set_prop(prop, reagent_idx, false)?;
        }
        reagent_idx += 1;
    }
    Ok(())
}

/// Everything collected while walking the document.
struct Doc {
    mols: Vec<RWMol>,
    fragment_lookup: HashMap<u32, usize>,
    grouped_fragments: BTreeMap<u32, Vec<i32>>,
    schemes: Vec<SchemeInfo>,
}

/// `visit_children`: fragments may be nested in groups, so the walk is recursive.
fn visit_children(
    node: &Node,
    state: &mut State,
    doc: &mut Doc,
    bond_length: f64,
    params: &CdxmlParserParams,
    mut group_id: i32,
) -> Result<()> {
    let zip_params = MolzipParams {
        label: MolzipLabel::AtomProperty,
        atom_property: FUSE_LABEL.to_string(),
        enforce_valence_rules: false,
        ..Default::default()
    };
    for frag in &node.children {
        if frag.name == "fragment" {
            let mut mol = RWMol::new();
            if !parse_fragment(&mut mol, frag, state, params, -1)? {
                continue;
            }
            let frag_id = match mol.get_prop(CDXML_FRAG_ID) {
                Some(RDValue::Int(i)) => *i as u32,
                _ => 0,
            };
            doc.fragment_lookup.insert(frag_id, doc.mols.len());
            if group_id != -1 {
                doc.grouped_fragments
                    .entry(group_id as u32)
                    .or_default()
                    .push(frag_id as i32);
            } else {
                doc.grouped_fragments
                    .entry(frag_id)
                    .or_default()
                    .push(frag_id as i32);
            }
            if mol.has_prop(NEEDS_FUSE) {
                mol.clear_prop(NEEDS_FUSE);
                let fused = if replace_fragments(&mut mol)? {
                    molzip_single(&mol, &zip_params)
                } else {
                    rd_log_warning!("Failed replacement of fragment skipping... {}", frag_id);
                    continue;
                };
                let mut fused = match fused {
                    Ok(Some(f)) => f,
                    _ => {
                        rd_log_warning!("Failed fusion of fragment skipping... {}", frag_id);
                        continue;
                    }
                };
                fused.set_prop(CDXML_FRAG_ID, frag_id as i32, false)?;
                doc.mols.push(fused);
            } else {
                doc.mols.push(mol);
            }
            let res = doc.mols.last_mut().expect("molecule");
            let mut conf = Conformer::new(res.num_atoms() as usize);
            conf.set_3d(false);
            let mut has_conf = false;
            for idx in 0..res.num_atoms() {
                let mut p = Point3D::new(0.0, 0.0, 0.0);
                if let Some(RDValue::VecDouble(coord)) = res.atom(idx).get_prop(CDX_ATOM_POS) {
                    has_conf = true;
                    if coord.len() == 2 {
                        p.x = coord[0];
                        // CDXML uses an inverted coordinate system
                        p.y = -coord[1];
                        p.z = 0.0;
                    }
                }
                conf.set_atom_pos(idx, p)?;
                res.atom_mut(idx).clear_prop(CDX_ATOM_POS);
            }
            if has_conf {
                scale_bonds(res, &mut conf, RDKIT_DEPICT_BONDLENGTH, bond_length);
                let conf_idx = res.add_conformer(conf, true)?;
                assign_chiral_types_from_bond_dirs(res, conf_idx as i32, true)?;
            }
            // the wedging is no longer needed once atom stereo is perceived
            clear_single_bond_dir_flags(res, false);
            if params.sanitize {
                let sanitized: Result<()> = (|| {
                    if params.remove_hs {
                        // bond stereo must be found before the Hs go away
                        sanitize_mol(res, SANITIZE_CLEANUP).map_err(|(_, e)| e)?;
                        detect_bond_stereochemistry(res, -1)?;
                        if params.parse_queries && has_query_hs(res)?.0 {
                            merge_query_hs(res, false, false)?;
                        }
                        remove_hs(res, &RemoveHsParameters::default(), true)?;
                    } else {
                        sanitize_mol(res, SANITIZE_ALL).map_err(|(_, e)| e)?;
                        detect_bond_stereochemistry(res, -1)?;
                    }
                    Ok(())
                })();
                if sanitized.is_err() {
                    rd_log_warning!(
                        "CDXMLParser: failed sanitizing skipping fragment {}",
                        frag_id
                    );
                    doc.mols.pop();
                    continue;
                }
                assign_stereochemistry(res, true, true, true)?;
            } else {
                detect_bond_stereochemistry(res, -1)?;
            }
        } else if frag.name == "scheme" {
            let scheme_id = attr_int(frag, "id", -1);
            for step in &frag.children {
                if step.name != "step" {
                    continue;
                }
                if step.attrs.is_empty() {
                    return Err(RDKitError::FileParse("No such node (<xmlattr>)".into()));
                }
                let mut scheme = SchemeInfo {
                    scheme_id,
                    step_id: attr_int(step, "id", -1),
                    ..Default::default()
                };
                for (key, value) in &step.attrs {
                    match key.as_str() {
                        "ReactionStepProducts" => scheme.products = to_vec_u32(value),
                        "ReactionStepReactants" => scheme.reactants = to_vec_u32(value),
                        "ReactionStepObjectsAboveArrow" => scheme.above_arrow = to_vec_u32(value),
                        "ReactionStepObjectsBelowArrow" => scheme.below_arrow = to_vec_u32(value),
                        "ReactionStepAtomMap" => scheme.atom_map = to_vec_u32(value),
                        _ => {}
                    }
                }
                doc.schemes.push(scheme);
            }
        } else {
            if frag.name == "group" {
                group_id = frag
                    .attr("id")
                    .and_then(ptree_int)
                    .ok_or_else(|| RDKitError::FileParse("No such node (<xmlattr>.id)".into()))?;
            }
            visit_children(frag, state, doc, bond_length, params, group_id)?;
        }
    }
    Ok(())
}

/// `MolsFromCDXML`: the molecules of a CDXML document.
pub fn mols_from_cdxml(cdxml: &str, params: &CdxmlParserParams) -> Result<Vec<RWMol>> {
    if params.format == CdxmlFormat::Cdx {
        return Err(RDKitError::FileParse(
            "Full ChemDraw support is not enabled, cannot parse CDX files".into(),
        ));
    }
    let roots = xml::parse(cdxml)?;
    let mut state = State {
        ids: HashMap::new(),
        missing_frag_id: -1,
    };
    let mut doc = Doc {
        mols: Vec::new(),
        fragment_lookup: HashMap::new(),
        grouped_fragments: BTreeMap::new(),
        schemes: Vec::new(),
    };
    for root in &roots {
        if root.name != "CDXML" {
            continue;
        }
        let bond_length = root
            .attr("BondLength")
            .and_then(|v| v.trim().parse::<f64>().ok())
            .ok_or_else(|| RDKitError::FileParse("No such node (<xmlattr>.BondLength)".into()))?;
        for node in &root.children {
            if node.name == "page" {
                visit_children(node, &mut state, &mut doc, bond_length, params, -1)?;
            }
        }
    }
    if !doc.schemes.is_empty() {
        let mut fragments: HashMap<u32, usize> = HashMap::new();
        let mut atoms: HashMap<u32, (usize, u32)> = HashMap::new();
        for (mol_idx, mol) in doc.mols.iter().enumerate() {
            let idx = match mol.get_prop(CDXML_FRAG_ID) {
                Some(RDValue::Int(i)) => *i as u32,
                _ => 0,
            };
            fragments.insert(idx, mol_idx);
            for a in 0..mol.num_atoms() {
                if let Some(id) = prop_u32(mol.atom(a), CDX_ATOM_ID) {
                    atoms.insert(id, (mol_idx, a));
                }
            }
        }
        let schemes = std::mem::take(&mut doc.schemes);
        for scheme in &schemes {
            set_reaction_data(
                "ReactionStepReactants",
                CDX_REAGENT_ID,
                scheme,
                &scheme.reactants,
                &fragments,
                &doc.grouped_fragments,
                &mut doc.mols,
            )?;
            set_reaction_data(
                "ReactionStepProducts",
                CDX_PRODUCT_ID,
                scheme,
                &scheme.products,
                &fragments,
                &doc.grouped_fragments,
                &mut doc.mols,
            )?;
            let mut agents = scheme.above_arrow.clone();
            agents.extend_from_slice(&scheme.below_arrow);
            set_reaction_data(
                "ReactionStepAgents",
                CDX_AGENT_ID,
                scheme,
                &agents,
                &fragments,
                &doc.grouped_fragments,
                &mut doc.mols,
            )?;
            if scheme.atom_map.len() % 2 != 0 {
                rd_log_warning!(
                    "CDXMLParser: Schema {} step {} ReactionStepAtomMap has odd number of entries, skipping schema...",
                    scheme.scheme_id,
                    scheme.step_id
                );
                continue;
            }
            for i in 0..scheme.atom_map.len() / 2 {
                for id in [scheme.atom_map[i * 2], scheme.atom_map[i * 2 + 1]] {
                    match atoms.get(&id) {
                        Some(&(m, a)) => {
                            doc.mols[m]
                                .atom_mut(a)
                                .set_atom_map_num(i as i32 + 1, false)?;
                        }
                        None => rd_log_warning!(
                            "CDXMLParser: Schema {} step {} ReactionStepAtomMap cannot find atom with node id {} skipping schema...",
                            scheme.scheme_id,
                            scheme.step_id,
                            id
                        ),
                    }
                }
            }
        }
    }
    clear_internal_cdx_props(&mut doc.mols);
    Ok(doc.mols)
}

/// `MolsFromCDXMLFile`
pub fn mols_from_cdxml_file(file_name: &str, params: &CdxmlParserParams) -> Result<Vec<RWMol>> {
    let text = std::fs::read(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    mols_from_cdxml(&String::from_utf8_lossy(&text), params)
}

/// `MolToCDXMLBlock` (writing needs the ChemDraw library).
pub fn mol_to_cdxml_block(_mol: &RWMol, _format: CdxmlFormat) -> Result<String> {
    Err(RDKitError::FileParse(
        "RDKit build withoutChemDraw writing support. ".into(),
    ))
}
