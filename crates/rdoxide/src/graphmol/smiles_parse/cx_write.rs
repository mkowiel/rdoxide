//! CXSMILES output: the `|...|` extension block written after a SMILES
//! (upstream `SmilesParse/CXSmilesOps.cpp`, the `SmilesWrite` half, and
//! `MolToCXSmiles` of `SmilesWrite.cpp`).
//!
//! Not ported: the link-node (`LN:`) block, which needs `MolEnumerator`.

use std::collections::BTreeMap;

use super::write::{SmilesWriteParams, mol_to_smiles_detail};
use crate::graphmol::atropisomers::{
    get_atropisomer_atoms_and_bonds, wedge_bonds_from_atropisomers,
};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO;
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::file_parsers::mol_file_writer::get_mol_file_bond_stereo_info;
use crate::graphmol::kekulize::{DEFAULT_MAX_BACKTRACKS, kekulize};
use crate::graphmol::mol_ops::cleanup_stereo_groups;
use crate::graphmol::romol::ROMol;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType, assign_stereo_group_ids};
use crate::graphmol::wedge_bonds::{
    WedgeBonds, WedgeInfo, pick_bonds_to_wedge, reapply_mol_block_wedging,
};
use crate::rd_log_warning;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::strings::format_g;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `SmilesWrite::CXSmilesFields`: which extensions to write.
pub mod cx_fields {
    pub const NONE: u32 = 0;
    pub const ATOM_LABELS: u32 = 1 << 0;
    pub const MOLFILE_VALUES: u32 = 1 << 1;
    pub const COORDS: u32 = 1 << 2;
    pub const RADICALS: u32 = 1 << 3;
    pub const ATOM_PROPS: u32 = 1 << 4;
    pub const LINKNODES: u32 = 1 << 5;
    pub const ENHANCEDSTEREO: u32 = 1 << 6;
    pub const SGROUPS: u32 = 1 << 7;
    pub const POLYMER: u32 = 1 << 8;
    pub const BOND_CFG: u32 = 1 << 9;
    pub const BOND_ATROPISOMER: u32 = 1 << 10;
    pub const COORDINATE_BONDS: u32 = 1 << 11;
    pub const HYDROGEN_BONDS: u32 = 1 << 12;
    pub const ZERO_BONDS: u32 = 1 << 13;
    pub const ALL: u32 = 0x7fff_ffff;
    pub const ALL_BUT_COORDS: u32 = ALL ^ COORDS;
}

/// `RestoreBondDirOption`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RestoreBondDirOption {
    /// Restore the wedging of the original mol block.
    True,
    /// Clear all bond direction information.
    Clear,
}

const PSEUDOATOMS: [&str; 2] = ["Pol", "Mod"];
/// The prefix of attachment point labels (`MolOps::attachmentPointLabelPrefix`).
const ATTACHMENT_POINT_LABEL_PREFIX: &str = "_AP";

// https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html#escaping
const SGROUP_ALLOWED: &str = "><\"!@#$%()[]./\\?-+*^_~= ";
const ATOM_PROP_ALLOWED: &str = "><\"!@#$%()[]./\\?-+*^_~= ";
const LABEL_ALLOWED: &str = "><\"!@#%()[]./\\?-+*^_~=,: ";

/// `sgroupTypemap` (CXSMILES code -> molfile type)
const SGROUP_TYPEMAP: [(&str, &str); 15] = [
    ("n", "SRU"),
    ("mon", "MON"),
    ("mer", "MER"),
    ("co", "COP"),
    ("xl", "CRO"),
    ("mod", "MOD"),
    ("mix", "MIX"),
    ("f", "FOR"),
    ("any", "ANY"),
    ("gen", "GEN"),
    ("c", "COM"),
    ("grf", "GRA"),
    ("alt", "COP"),
    ("ran", "COP"),
    ("blk", "COP"),
];

/// `quote_string`: characters that are not alphanumeric or allowed are
/// written as `&#code;`.
fn quote_string(txt: &str, allowed: &str) -> String {
    let mut res = String::new();
    for &c in txt.as_bytes() {
        if !allowed.is_empty() && !c.is_ascii_alphanumeric() && !allowed.as_bytes().contains(&c) {
            // `static_cast<unsigned>(char)`: bytes above 127 sign-extend
            res += &format!("&#{};", i32::from(c as i8) as u32);
        } else {
            res.push(c as char);
        }
    }
    res
}

fn string_prop(p: &impl RDProps, key: &str) -> Option<String> {
    p.get_prop(key).map(RDValue::to_rd_string)
}

fn uint_prop(p: &impl RDProps, key: &str) -> Option<u32> {
    p.get_prop(key).and_then(RDValue::as_uint)
}

fn uint_vec_prop(p: &impl RDProps, key: &str) -> Option<Vec<u32>> {
    match p.get_prop(key)? {
        RDValue::VecUnsignedInt(v) => Some(v.clone()),
        RDValue::VecInt(v) => Some(v.iter().map(|&i| i as u32).collect()),
        _ => None,
    }
}

/// Position of every atom (or bond) in an output order.
fn reverse_order(order: &[u32], n: u32) -> Vec<u32> {
    let mut rev = vec![0u32; n as usize];
    for (i, &o) in order.iter().enumerate() {
        rev[o as usize] = i as u32;
    }
    rev
}

/// `std::find(order, v) - order.begin()`
fn position(order: &[u32], v: u32) -> u32 {
    order.iter().position(|&o| o == v).unwrap_or(order.len()) as u32
}

fn join<T: ToString>(items: &[T], sep: &str) -> String {
    items.iter().map(T::to_string).collect::<Vec<_>>().join(sep)
}

// ------------------------------------------------------------ enhanced stereo

/// `Atropisomers::getAllAtomIdsForStereoGroup`: the atoms of the group plus
/// the atoms next to a group bond that carry a wedge or hash.
fn stereo_group_atom_ids(mol: &ROMol, group: &StereoGroup, wedge: &WedgeBonds) -> Vec<u32> {
    let mut ids: Vec<u32> = group.atoms().to_vec();
    for &bond in group.bonds() {
        let b = mol.bond(bond);
        for atom in [b.begin_atom_idx(), b.end_atom_idx()] {
            for other in mol.atom_bonds(atom) {
                if other.idx() == bond {
                    continue;
                }
                let marked = matches!(other.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash)
                    || matches!(wedge.get(&other.idx()), Some(WedgeInfo::Atropisomer(_)));
                if marked && !ids.contains(&atom) {
                    ids.push(atom);
                }
            }
        }
    }
    ids
}

fn get_enhanced_stereo_block(mol: &ROMol, atom_order: &[u32], wedge: &WedgeBonds) -> String {
    if mol.stereo_groups().is_empty() {
        return String::new();
    }
    let rev = reverse_order(atom_order, mol.num_atoms());
    let mut sorted: Vec<(StereoGroup, Vec<u32>)> = mol
        .stereo_groups()
        .iter()
        .filter_map(|sg| {
            let mut idxs: Vec<u32> = stereo_group_atom_ids(mol, sg, wedge)
                .into_iter()
                .map(|a| rev[a as usize])
                .collect();
            idxs.sort_unstable();
            (!idxs.is_empty()).then(|| (sg.clone(), idxs))
        })
        .collect();
    // by type, write id and atom indices
    sorted.sort_by(|(a, ia), (b, ib)| {
        (a.group_type() as u8, a.write_id(), ia).cmp(&(b.group_type() as u8, b.write_id(), ib))
    });
    let mut groups: Vec<StereoGroup> = sorted.iter().map(|(g, _)| g.clone()).collect();
    assign_stereo_group_ids(&mut groups);
    let mut res = String::new();
    for (g, (_, atoms)) in groups.iter().zip(&sorted) {
        match g.group_type() {
            StereoGroupType::Absolute => res += "a:",
            StereoGroupType::Or => res += &format!("o{}:", g.write_id()),
            StereoGroupType::And => res += &format!("&{}:", g.write_id()),
        }
        for a in atoms {
            res += &format!("{a},");
        }
    }
    if res.ends_with(',') {
        res.pop();
    }
    res
}

// ------------------------------------------------------------ substance groups

fn sgroup_hierarchy_block(mol: &mut ROMol) -> String {
    if mol.substance_groups().is_empty() {
        return String::new();
    }
    let mut order: BTreeMap<u32, u32> = BTreeMap::new();
    let mut parent_present = false;
    let count = mol.substance_groups().len();
    for i in 0..count {
        let sg = &mut mol.substance_groups_mut()[i];
        if sg.has_prop("_cxsmilesOutputIndex") {
            let idx = uint_prop(sg, "index").unwrap_or(i as u32);
            order.insert(idx, uint_prop(sg, "_cxsmilesOutputIndex").unwrap_or(0));
            sg.clear_prop("_cxsmilesOutputIndex");
        }
        parent_present |= sg.has_prop("PARENT");
    }
    if !parent_present {
        return String::new();
    }
    let mut accum: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (i, sg) in mol.substance_groups().iter().enumerate() {
        if let Some(parent) = uint_prop(sg, "PARENT")
            && order.contains_key(&parent)
        {
            let idx = uint_prop(sg, "index").unwrap_or(i as u32);
            if order.contains_key(&idx) {
                accum.entry(order[&parent]).or_default().push(order[&idx]);
            }
        }
    }
    if accum.is_empty() {
        return String::new();
    }
    let mut res = String::from("SgH:");
    for (parent, children) in &accum {
        res += &format!("{parent}:{},", join(children, "."));
    }
    res.trim_end_matches(',').to_string()
}

fn sgroup_polymer_block(mol: &mut ROMol, atom_order: &[u32], bond_order: &[u32]) -> String {
    if mol.substance_groups().is_empty() {
        return String::new();
    }
    let mut output_index = uint_prop(mol, "_cxsmilesOutputIndex").unwrap_or(0);
    let rev_atoms = reverse_order(atom_order, mol.num_atoms());
    // the first code (alphabetically) of every type
    let mut reverse_typemap: BTreeMap<&str, &str> = BTreeMap::new();
    let mut sorted_map = SGROUP_TYPEMAP;
    sorted_map.sort();
    for (code, typ) in sorted_map {
        reverse_typemap.entry(typ).or_insert(code);
    }
    let mut res = String::new();
    for sg in mol.substance_groups_mut() {
        let Some(typ) = string_prop(sg, "TYPE") else {
            continue;
        };
        let Some(&code) = reverse_typemap.get(typ.as_str()) else {
            continue;
        };
        sg.set_prop("_cxsmilesOutputIndex", output_index, false)
            .ok();
        output_index += 1;

        res += "Sg:";
        match (typ.as_str(), string_prop(sg, "SUBTYPE").as_deref()) {
            ("COP", Some("ALT")) => res += "alt",
            ("COP", Some("RAN")) => res += "ran",
            ("COP", Some("BLO")) => res += "blk",
            _ => res += code,
        }
        res.push(':');
        for &a in sg.atoms() {
            res += &format!("{},", rev_atoms[a as usize]);
        }
        // overwrites the last "," (or the ":" when there are no atoms)
        res.pop();
        res.push(':');
        if let Some(label) = string_prop(sg, "LABEL") {
            res += &quote_string(&label, SGROUP_ALLOWED);
        }
        res.push(':');
        if let Some(connect) = string_prop(sg, "CONNECT") {
            res += &connect.to_lowercase();
        }
        res.push(':');
        if let Some(head) = uint_vec_prop(sg, "XBHEAD").filter(|v| v.len() > 1) {
            res += &join(
                &head
                    .iter()
                    .map(|&v| bond_order[v as usize])
                    .collect::<Vec<_>>(),
                ",",
            );
        }
        res.push(':');
        if let Some(tail) = uint_vec_prop(sg, "XBCORR").filter(|v| v.len() > 2) {
            let bonds: Vec<u32> = tail
                .iter()
                .skip(1)
                .step_by(2)
                .map(|&v| bond_order[v as usize])
                .collect();
            res += &join(&bonds, ",");
        }
        res += ":,";
    }
    while res.ends_with(',') {
        res.pop();
    }
    mol.set_prop("_cxsmilesOutputIndex", output_index, false)
        .ok();
    res
}

fn sgroup_data_block(mol: &mut ROMol, atom_order: &[u32]) -> String {
    if mol.substance_groups().is_empty() {
        return String::new();
    }
    let mut output_index = uint_prop(mol, "_cxsmilesOutputIndex").unwrap_or(0);
    let rev = reverse_order(atom_order, mol.num_atoms());
    let mut res = String::new();
    for sg in mol.substance_groups_mut() {
        if string_prop(sg, "TYPE").as_deref() != Some("DAT") {
            continue;
        }
        sg.set_prop("_cxsmilesOutputIndex", output_index, false)
            .ok();
        output_index += 1;
        res += "SgD:";
        // the atom order is not canonicalized: it may carry meaning
        for &a in sg.atoms() {
            res += &format!("{},", rev[a as usize]);
        }
        res.pop();
        res.push(':');
        if let Some(name) = string_prop(sg, "FIELDNAME").filter(|s| !s.is_empty()) {
            res += &quote_string(&name, SGROUP_ALLOWED);
        }
        res.push(':');
        if let Some(RDValue::VecString(fields)) = sg.get_prop("DATAFIELDS")
            && !fields.is_empty()
        {
            res += &fields
                .iter()
                .map(|f| quote_string(f, SGROUP_ALLOWED))
                .collect::<Vec<_>>()
                .join(",");
        }
        res.push(':');
        if let Some(op) = string_prop(sg, "QUERYOP").filter(|s| !s.is_empty()) {
            res += &op;
        }
        res.push(':');
        if let Some(info) = string_prop(sg, "FIELDINFO").filter(|s| !s.is_empty()) {
            res += &quote_string(&info, SGROUP_ALLOWED);
        }
        res.push(':');
        if let Some(tag) = string_prop(sg, "FIELDTAG").filter(|s| !s.is_empty()) {
            res += &quote_string(&tag, SGROUP_ALLOWED);
        }
        // the coordinates are not written
        res += ":,";
    }
    if res.ends_with(',') {
        res.pop();
    }
    mol.set_prop("_cxsmilesOutputIndex", output_index, false)
        .ok();
    res
}

// ------------------------------------------------------------ atom blocks

fn atom_label_block(mol: &ROMol, atom_order: &[u32]) -> String {
    let mut labels = Vec::with_capacity(atom_order.len());
    for &idx in atom_order {
        let atom = mol.atom(idx);
        let dummy_label = string_prop(atom, cp::dummyLabel);
        let attach = atom
            .get_prop(cp::_fromAttachPoint)
            .and_then(RDValue::as_int);
        let label = if let Some(l) = string_prop(atom, cp::_QueryAtomGenericLabel) {
            quote_string(&format!("{l}_p"), LABEL_ALLOWED)
        } else if atom.atomic_num() == 0
            && let Some(l) = dummy_label.filter(|l| PSEUDOATOMS.contains(&l.as_str()))
        {
            quote_string(&format!("{l}_p"), LABEL_ALLOWED)
        } else if atom.atomic_num() == 0
            && let Some(v @ (1 | 2)) = attach
        {
            quote_string(
                &format!("{ATTACHMENT_POINT_LABEL_PREFIX}{v}"),
                LABEL_ALLOWED,
            )
        } else if let Some(l) = string_prop(atom, cp::atomLabel) {
            quote_string(&l, LABEL_ALLOWED)
        } else {
            String::new()
        };
        labels.push(label);
    }
    let res = labels.join(";");
    // nothing found: no block
    if res.chars().all(|c| c == ';') {
        String::new()
    } else {
        res
    }
}

fn value_block(mol: &ROMol, atom_order: &[u32], prop: &str) -> String {
    atom_order
        .iter()
        .map(|&idx| {
            string_prop(mol.atom(idx), prop)
                .map_or_else(String::new, |l| quote_string(&l, ATOM_PROP_ALLOWED))
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn radical_block(mol: &ROMol, atom_order: &[u32]) -> String {
    let mut rads: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (i, &idx) in atom_order.iter().enumerate() {
        let n = mol.atom(idx).num_radical_electrons();
        if n > 0 {
            rads.entry(n).or_default().push(i);
        }
    }
    let mut res = String::new();
    for (n, atoms) in rads {
        match n {
            1 => res += "^1:",
            2 => res += "^2:",
            3 => res += "^5:",
            _ => rd_log_warning!("unsupported number of radical electrons {n}"),
        }
        for a in atoms {
            res += &format!("{a},");
        }
    }
    res
}

fn zero_small_vals(v: f64) -> f64 {
    if v.abs() < 1e-4 { 0.0 } else { v }
}

fn coords_block(mol: &ROMol, atom_order: &[u32]) -> Result<String> {
    let conf = mol.conformer(-1)?;
    let mut parts = Vec::with_capacity(atom_order.len());
    for &idx in atom_order {
        let pt = conf.atom_pos(idx)?;
        let mut s = format!(
            "{},{},",
            format_g(zero_small_vals(pt.x), 6),
            format_g(zero_small_vals(pt.y), 6)
        );
        if conf.is_3d() {
            let z = format_g(zero_small_vals(pt.z), 6);
            if z != "0" {
                s += &z;
            }
        }
        parts.push(s);
    }
    Ok(parts.join(";"))
}

fn atom_props_block(mol: &ROMol, atom_order: &[u32]) -> String {
    let skip = [
        cp::atomLabel,
        cp::molFileValue,
        cp::molParity,
        cp::molAtomMapNumber,
        cp::molStereoCare,
        cp::molRxnExactChange,
        cp::molInversionFlag,
    ];
    let mut res = String::new();
    for (which, &idx) in atom_order.iter().enumerate() {
        let atom = mol.atom(idx);
        let is_attachment = atom.atomic_num() == 0 && atom.has_prop(cp::_fromAttachPoint);
        for name in atom.get_prop_list(false, false) {
            if skip.contains(&name.as_str()) {
                continue;
            }
            let value = string_prop(atom, &name).unwrap_or_default();
            if name == "dummyLabel"
                && (is_attachment || value == "*" || PSEUDOATOMS.contains(&value.as_str()))
            {
                // a pseudoatom or attachment point
                continue;
            }
            if res.is_empty() {
                res += "atomProp";
            }
            res += &format!(
                ":{which}.{}.{}",
                quote_string(&name, ATOM_PROP_ALLOWED),
                quote_string(&value, ATOM_PROP_ALLOWED)
            );
        }
    }
    res
}

// ------------------------------------------------------------ bond blocks

/// Does the wedge of `bond` (at `first_atom`) flip when the atoms are
/// reordered for the SMILES? Only needed for atropisomers without coordinates.
fn atropisomer_wedge_flips(
    mol: &ROMol,
    nbr_bond: u32,
    first_atom: u32,
    atom_order: &[u32],
) -> Result<bool> {
    let ends = get_atropisomer_atoms_and_bonds(mol, nbr_bond)
        .ok_or_else(|| RDKitError::Value("Internal error - should not occur".into()))?;
    let nb = mol.bond(nbr_bond);
    let mut swaps = 0;
    if position(atom_order, nb.begin_atom_idx()) > position(atom_order, nb.end_atom_idx()) {
        swaps += 1;
    }
    for (atom, bonds) in &ends {
        if *atom == first_atom {
            // swapping atoms on the side of the wedge does not change it
            continue;
        }
        if bonds.len() == 2 {
            let other = |b: u32| mol.bond(b).other_atom_idx(*atom);
            let (a, b) = (other(bonds[0])?, other(bonds[1])?);
            if position(atom_order, a) > position(atom_order, b) {
                swaps += 1;
            }
        }
    }
    Ok(swaps % 2 == 1)
}

fn bond_config_block(
    mol: &ROMol,
    atom_order: &[u32],
    bond_order: &[u32],
    coords_included: bool,
    wedge: &WedgeBonds,
    conf: Option<&Conformer>,
    atropisomer_only: bool,
) -> Result<String> {
    let mut parts: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    for (i, &idx) in bond_order.iter().enumerate() {
        let bond = mol.bond(idx);
        let mut wedge_start = bond.begin_atom_idx();
        if !crate::graphmol::canon::can_have_direction(bond.bond_type()) {
            continue;
        }
        // favor the wedge state of the bond
        let mut bd = match bond.bond_dir() {
            d @ (BondDir::BeginDash | BondDir::BeginWedge | BondDir::Unknown) => d,
            _ => BondDir::None,
        };
        if atropisomer_only && bd == BondDir::None {
            continue;
        }
        // is it a wedge next to an atropisomeric bond?
        let mut is_atropisomer = false;
        if matches!(bd, BondDir::BeginDash | BondDir::BeginWedge) {
            let first_atom = bond.begin_atom_idx();
            for nbr in mol.atom_bonds(first_atom) {
                if nbr.idx() == idx {
                    continue;
                }
                if matches!(
                    nbr.stereo(),
                    BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW
                ) {
                    is_atropisomer = true;
                    if !coords_included
                        && atropisomer_wedge_flips(mol, nbr.idx(), first_atom, atom_order)?
                    {
                        bd = if bd == BondDir::BeginWedge {
                            BondDir::BeginDash
                        } else {
                            BondDir::BeginWedge
                        };
                    }
                    break;
                }
            }
        }
        if atropisomer_only {
            if !is_atropisomer {
                continue;
            }
        } else {
            // a wedging caused by a chiral atom
            if bd == BondDir::None
                && let Some(cfg) = uint_prop(bond, cp::_MolFileBondCfg)
            {
                bd = match cfg {
                    1 => BondDir::BeginWedge,
                    2 => BondDir::Unknown,
                    3 => BondDir::BeginDash,
                    _ => BondDir::None,
                };
            }
            if bd == BondDir::None && coords_included {
                let (code, reverse) = get_mol_file_bond_stereo_info(mol, idx, wedge, conf)?;
                bd = match code {
                    1 => BondDir::BeginWedge,
                    3 => BondDir::Unknown,
                    6 => BondDir::BeginDash,
                    _ => BondDir::None,
                };
                if reverse {
                    wedge_start = bond.end_atom_idx();
                }
            }
        }
        let beg_order = position(atom_order, wedge_start);
        let kind = match bd {
            BondDir::Unknown => "w",
            // up and down wedges only with coordinates or for atropisomers
            BondDir::BeginWedge if coords_included || is_atropisomer => "wU",
            BondDir::BeginDash if coords_included || is_atropisomer => "wD",
            _ => continue,
        };
        parts
            .entry(kind)
            .or_default()
            .push(format!("{beg_order}.{i}"));
    }
    Ok(parts
        .into_iter()
        .map(|(k, v)| format!("{k}:{}", v.join(",")))
        .collect::<Vec<_>>()
        .join(","))
}

/// The `C:` (coordinate bonds) and `H:` (hydrogen bonds) blocks.
fn bond_type_block(
    mol: &ROMol,
    bond_type: BondType,
    symbol: &str,
    atom_order: &[u32],
    bond_order: &[u32],
) -> String {
    let items: Vec<String> = bond_order
        .iter()
        .enumerate()
        .filter(|&(_, &idx)| mol.bond(idx).bond_type() == bond_type)
        .map(|(i, &idx)| {
            format!(
                "{}.{i}",
                position(atom_order, mol.bond(idx).begin_atom_idx())
            )
        })
        .collect();
    if items.is_empty() {
        String::new()
    } else {
        format!("{symbol}:{}", items.join(","))
    }
}

fn zero_bonds_block(mol: &ROMol, bond_order: &[u32]) -> String {
    let items: Vec<String> = bond_order
        .iter()
        .enumerate()
        .filter(|&(_, &idx)| mol.bond(idx).bond_type() == BondType::Zero)
        .map(|(i, _)| i.to_string())
        .collect();
    if items.is_empty() {
        String::new()
    } else {
        format!("Z:{}", items.join(","))
    }
}

/// The `c:`, `t:` and `ctu:` blocks of ring double bonds.
fn ringbond_cistrans_block(mol: &ROMol, atom_order: &[u32], bond_order: &[u32]) -> String {
    let rinfo = mol.ring_info();
    if !rinfo.is_initialized() {
        return String::new();
    }
    let (mut c, mut t, mut ctu) = (String::new(), String::new(), String::new());
    let push = |s: &mut String, prefix: &str, label: &str| {
        s.push_str(if s.is_empty() { prefix } else { "," });
        s.push_str(label);
    };
    for (i, &idx) in bond_order.iter().enumerate() {
        if rinfo.num_bond_rings(idx).unwrap_or(0) == 0
            || rinfo.min_bond_ring_size(idx).unwrap_or(0) < MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO
        {
            continue;
        }
        let bond = mol.bond(idx);
        if !matches!(bond.bond_type(), BondType::Double | BondType::Aromatic) {
            continue;
        }
        let stereo = bond.stereo();
        if !matches!(
            stereo,
            BondStereo::StereoAny | BondStereo::StereoCis | BondStereo::StereoTrans
        ) {
            continue;
        }
        let label = i.to_string();
        if stereo == BondStereo::StereoAny {
            push(&mut ctu, "ctu:", &label);
            continue;
        }
        let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let stereo_atoms = bond.stereo_atoms();
        let mut need_swap = false;
        // upstream indexes `atomOrder` with a comparison result here
        // (`atomOrder[nbr->getIdx() < o1]`); the same is done on purpose
        for (atom, other_end, sa) in [(beg, end, stereo_atoms[0]), (end, beg, stereo_atoms[1])] {
            if mol.atom_degree(atom) > 2 {
                let o1 = atom_order[sa as usize];
                for nbr in mol.atom_neighbors(atom) {
                    if nbr == other_end || nbr == sa as u32 {
                        continue;
                    }
                    if atom_order[usize::from(nbr < o1)] != 0 {
                        need_swap = !need_swap;
                    }
                }
            }
        }
        if stereo == BondStereo::StereoCis || need_swap {
            push(&mut c, "c:", &label);
        } else {
            push(&mut t, "t:", &label);
        }
    }
    c + &t + &ctu
}

fn append_to_extension(addition: &str, base: &mut String) {
    if !addition.is_empty() {
        if base.len() > 1 {
            base.push(',');
        }
        base.push_str(addition);
    }
}

/// `checkCXFeatures`
fn check_cx_features(mol: &ROMol) {
    if mol.has_prop(cp::molFileLinkNodes) {
        rd_log_warning!("CX Extensions: mol has link nodes which are not currently supported");
    }
    if mol
        .substance_groups()
        .iter()
        .any(|sg| sg.has_prop("PARENT"))
    {
        rd_log_warning!("CX Extensions: Substance group hierarchy is not always preserved.");
    }
}

/// `getCXExtensions(mol, flags)`: the `|...|` block for a molecule that was
/// just written to SMILES (`_smilesAtomOutputOrder` and
/// `_smilesBondOutputOrder` must be set); empty if there is nothing to write.
pub fn get_cx_extensions(mol: &mut ROMol, flags: u32) -> Result<String> {
    use cx_fields as f;
    check_cx_features(mol);
    let (Some(atom_order), Some(bond_order)) = (
        uint_vec_prop(mol, cp::_smilesAtomOutputOrder),
        uint_vec_prop(mol, cp::_smilesBondOutputOrder),
    ) else {
        return Err(RDKitError::Value(
            "Input molecule does not have the required smiles ordering properties set".into(),
        ));
    };
    let mut res = String::from("|");

    let need_labels = atom_order.iter().any(|&i| {
        let a = mol.atom(i);
        a.has_prop(cp::atomLabel)
            || a.has_prop(cp::_QueryAtomGenericLabel)
            || a.has_prop(cp::dummyLabel)
            || a.has_prop(cp::_fromAttachPoint)
    });
    let need_values = atom_order
        .iter()
        .any(|&i| mol.atom(i).has_prop(cp::molFileValue));
    if flags & f::COORDS != 0 && mol.num_conformers() > 0 {
        res += &format!("({})", coords_block(mol, &atom_order)?);
    }
    if flags & f::ATOM_LABELS != 0 && need_labels {
        let labels = atom_label_block(mol, &atom_order);
        if !labels.is_empty() {
            if res.len() > 1 {
                res.push(',');
            }
            res += &format!("${labels}$");
        }
    }
    if flags & f::MOLFILE_VALUES != 0 && need_values {
        if res.len() > 1 {
            res.push(',');
        }
        res += &format!("$_AV:{}$", value_block(mol, &atom_order, cp::molFileValue));
    }
    let radicals = radical_block(mol, &atom_order);
    if flags & f::RADICALS != 0 && !radicals.is_empty() {
        if res.len() > 1 {
            res.push(',');
        }
        res += &radicals;
        if res.ends_with(',') {
            res.pop();
        }
    }
    if flags & f::ATOM_PROPS != 0 {
        append_to_extension(&atom_props_block(mol, &atom_order), &mut res);
    }

    let conf: Option<Conformer> = (mol.num_conformers() > 0 && flags & f::COORDS != 0)
        .then(|| mol.conformer(-1).cloned())
        .transpose()?;
    let mut wedge = WedgeBonds::new();
    if flags & f::BOND_CFG != 0 {
        wedge = pick_bonds_to_wedge(mol, conf.as_ref())?;
        let coords_included = flags & f::COORDS != 0 && mol.num_conformers() > 0;
        let block = bond_config_block(
            mol,
            &atom_order,
            &bond_order,
            coords_included,
            &wedge,
            conf.as_ref(),
            false,
        )?;
        append_to_extension(&block, &mut res);
        append_to_extension(
            &ringbond_cistrans_block(mol, &atom_order, &bond_order),
            &mut res,
        );
    } else if flags & f::BOND_ATROPISOMER != 0 {
        // BOND_CFG includes the atropisomer wedging
        wedge_bonds_from_atropisomers(mol, conf.as_ref(), &mut wedge)?;
        let block = bond_config_block(
            mol,
            &atom_order,
            &bond_order,
            conf.is_some(),
            &wedge,
            conf.as_ref(),
            true,
        )?;
        append_to_extension(&block, &mut res);
    }
    if flags & f::COORDINATE_BONDS != 0 {
        append_to_extension(
            &bond_type_block(mol, BondType::Dative, "C", &atom_order, &bond_order),
            &mut res,
        );
    }
    if flags & f::HYDROGEN_BONDS != 0 {
        append_to_extension(
            &bond_type_block(mol, BondType::Hydrogen, "H", &atom_order, &bond_order),
            &mut res,
        );
    }
    if flags & f::ZERO_BONDS != 0 {
        append_to_extension(&zero_bonds_block(mol, &bond_order), &mut res);
    }
    if flags & f::ENHANCEDSTEREO != 0 {
        append_to_extension(
            &get_enhanced_stereo_block(mol, &atom_order, &wedge),
            &mut res,
        );
    }
    if flags & f::SGROUPS != 0 {
        let block = sgroup_data_block(mol, &atom_order);
        append_to_extension(&block, &mut res);
    }
    if flags & f::POLYMER != 0 {
        let block = sgroup_polymer_block(mol, &atom_order, &bond_order);
        append_to_extension(&block, &mut res);
    }
    if flags & (f::SGROUPS | f::POLYMER) != 0 {
        let block = sgroup_hierarchy_block(mol);
        append_to_extension(&block, &mut res);
    }
    mol.clear_prop("_cxsmilesOutputIndex");
    Ok(if res.len() > 1 {
        res + "|"
    } else {
        String::new()
    })
}

/// `MolToCXSmiles`: canonical SMILES followed by the extension block.
pub fn mol_to_cx_smiles(
    mol: &ROMol,
    params: &SmilesWriteParams,
    mut flags: u32,
    restore_bond_dirs: RestoreBondDirOption,
) -> Result<String> {
    let mut trwmol = mol.clone();
    let include_stereo_groups = flags & cx_fields::ENHANCEDSTEREO != 0;
    let mut params = params.clone();
    // kekulize here, so that the wedges are not put on bonds that would
    // otherwise be aromatic
    if params.do_kekule {
        kekulize(&mut trwmol, false, false, DEFAULT_MAX_BACKTRACKS)?;
        params.do_kekule = false;
    }
    let mut res = mol_to_smiles_detail(&mut trwmol, &params, true, include_stereo_groups)?;
    if res.is_empty() {
        return Ok(res);
    }
    match restore_bond_dirs {
        RestoreBondDirOption::True => reapply_mol_block_wedging(&mut trwmol, true, false)?,
        RestoreBondDirOption::Clear => {
            for idx in 0..trwmol.num_bonds() {
                if !crate::graphmol::canon::can_have_direction(trwmol.bond(idx).bond_type()) {
                    continue;
                }
                // keep wiggly bonds for the bond config block
                let cfg = uint_prop(trwmol.bond(idx), cp::_MolFileBondCfg);
                let bond = trwmol.bond_mut(idx);
                if cfg == Some(2) {
                    bond.set_bond_dir(BondDir::Unknown);
                } else {
                    if bond.bond_dir() != BondDir::None {
                        bond.set_bond_dir(BondDir::None);
                    }
                    bond.clear_prop(cp::_MolFileBondCfg);
                }
            }
        }
    }
    if !params.do_isomeric_smiles {
        flags &= !(cx_fields::ENHANCEDSTEREO | cx_fields::BOND_CFG);
    }
    if params.clean_stereo {
        if trwmol.needs_update_property_cache() {
            trwmol.update_property_cache(false)?;
        }
        assign_stereochemistry(&mut trwmol, true, false, false)?;
        cleanup_stereo_groups(&mut trwmol)?;
    }
    let ext = get_cx_extensions(&mut trwmol, flags)?;
    if !ext.is_empty() {
        res += " ";
        res += &ext;
    }
    Ok(res)
}
