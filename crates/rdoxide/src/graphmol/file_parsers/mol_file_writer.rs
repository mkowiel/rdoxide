//! MDL Mol block writer (`MolToMolBlock`, `MolToV3KMolBlock`,
//! `MolToV2KMolBlock`) and bond wedging (`Chirality::pickBondsToWedge`).
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MolFileWriter.cpp`,
//! `Code/GraphMol/WedgeBonds.cpp`, `Code/GraphMol/Chirality.cpp`
//! (`GetMolFileBondStereoInfo` and helpers) and `StereoGroup.cpp`
//! (`assignStereoGroupIds`), 2026.09.1.
//!
//! Not ported yet (reported as errors): SubstanceGroups (including the ones
//! the writer creates for SMARTS queries, generic groups and zero-order
//! bonds), wedging of atropisomers, and 2D coordinate generation for
//! molecules without a conformer when `includeStereo` is set.
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::f64::consts::PI;

use crate::geometry::Point3D;
use crate::graphmol::atom::{ChiralType, ValenceType};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::get_use_legacy_stereo_perception;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::{
    AtomQuery, BondQuery, get_atom_list_query_vals, is_atom_list_query,
};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::in_organic_subset;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

const MAX_V2000_COORD: f64 = 100000.;
const MIN_V2000_COORD: f64 = -10000.;

/// `MolWriterParams`
#[derive(Debug, Clone, PartialEq)]
pub struct MolWriterParams {
    pub include_stereo: bool,
    pub kekulize: bool,
    pub force_v3000: bool,
    pub precision: u32,
}

impl Default for MolWriterParams {
    fn default() -> Self {
        MolWriterParams {
            include_stereo: true,
            kekulize: true,
            force_v3000: false,
            precision: 6,
        }
    }
}

fn not_ported(what: &str) -> RDKitError {
    RDKitError::Runtime(format!("{what} is not supported by rdoxide yet"))
}

fn int_prop(props: &impl RDProps, k: &str) -> Option<i32> {
    props.get_prop(k).and_then(RDValue::as_int)
}

fn uint_prop(props: &impl RDProps, k: &str) -> Option<u32> {
    props.get_prop(k).and_then(RDValue::as_uint)
}

fn str_prop(props: &impl RDProps, k: &str) -> Option<String> {
    match props.get_prop(k) {
        Some(RDValue::String(s)) => Some(s.clone()),
        Some(v) => Some(v.to_rd_string()),
        None => None,
    }
}

// ------------------------------------------------------------------ wedging

/// `Chirality::WedgeInfoBase`
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WedgeInfo {
    /// `WedgeInfoChiral`: wedge from this atom
    Chiral(u32),
    /// `WedgeInfoAtropisomer`
    Atropisomer(BondDir),
}

pub type WedgeBonds = BTreeMap<u32, WedgeInfo>;

fn is_tetra(t: ChiralType) -> bool {
    matches!(t, ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW)
}

/// `detail::countChiralNbrs`
fn count_chiral_nbrs(mol: &ROMol, no_nbrs: i32) -> (bool, Vec<i32>) {
    let mut n = vec![no_nbrs; mol.num_atoms() as usize];
    for b in mol.bonds() {
        if matches!(
            b.bond_dir(),
            BondDir::BeginWedge | BondDir::BeginDash | BondDir::Unknown
        ) {
            if is_tetra(mol.atom(b.begin_atom_idx()).chiral_tag()) {
                n[b.begin_atom_idx() as usize] = no_nbrs + 1;
            } else if is_tetra(mol.atom(b.end_atom_idx()).chiral_tag()) {
                n[b.end_atom_idx() as usize] = no_nbrs + 1;
            }
        }
    }
    let mut chi_nbrs = false;
    for at in mol.atoms() {
        let i = at.idx() as usize;
        if n[i] > no_nbrs || !is_tetra(at.chiral_tag()) {
            continue;
        }
        n[i] = 0;
        chi_nbrs = true;
        for nb in mol.atom_neighbors(at.idx()) {
            let na = mol.atom(nb);
            if na.atomic_num() == 1 {
                n[i] -= 10;
                continue;
            }
            if is_tetra(na.chiral_tag()) {
                n[i] -= 1;
            }
        }
    }
    (chi_nbrs, n)
}

/// `getDoubleBondPresence`
fn double_bond_presence(mol: &ROMol, atom: u32) -> (i32, i32, i32) {
    let (mut d, mut known, mut any) = (0, 0, 0);
    for b in mol.atom_bonds(atom) {
        if b.bond_type() == BondType::Double {
            d += 1;
            if b.stereo() == BondStereo::StereoAny {
                any += 1;
            } else if b.stereo() as u8 > BondStereo::StereoAny as u8 {
                known += 1;
            }
        }
    }
    (d, known, any)
}

/// `detail::pickBondToWedge`
fn pick_bond_to_wedge(
    mol: &mut ROMol,
    atom: u32,
    n_chiral_nbrs: &[i32],
    wedge_bonds: &WedgeBonds,
    no_nbrs: i32,
) -> Result<Option<u32>> {
    if !mol.ring_info().is_sssr_or_better() {
        find_sssr(mol, false, false)?;
    }
    let mut scores: Vec<(i32, u32)> = Vec::new();
    for bid in mol.atom_bond_indices(atom).collect::<Vec<_>>() {
        let b = mol.bond(bid);
        if b.bond_type() != BondType::Single || wedge_bonds.contains_key(&bid) {
            continue;
        }
        let o = b.other_atom_idx(atom)?;
        let oa = mol.atom(o);
        if oa.atomic_num() == 1 {
            scores.push((-1000000, bid));
            continue;
        }
        let mut score = oa.atomic_num() as i32
            + 100 * mol.atom_degree(o) as i32
            + 1000 * i32::from(oa.chiral_tag() != ChiralType::Unspecified);
        if n_chiral_nbrs[o as usize] < no_nbrs {
            score -= 100000 * n_chiral_nbrs[o as usize];
        }
        score += 10000 * mol.ring_info().num_atom_rings(o)? as i32;
        score += 20000 * mol.ring_info().num_bond_rings(bid)? as i32;
        let (d, known, any) = double_bond_presence(mol, o);
        score += 11000 * d + 12000 * known + 23000 * any;
        if oa.has_prop(cp::_fromAttachPoint) {
            score += 500000;
        }
        scores.push((score, bid));
    }
    // std::min_element: the first minimum
    Ok(scores.iter().min_by_key(|p| *p).map(|p| p.1))
}

/// `Chirality::pickBondsToWedge(mol, nullptr, conf)`
pub fn pick_bonds_to_wedge(mol: &mut ROMol, conf: Option<&Conformer>) -> Result<WedgeBonds> {
    const NO_NBRS: i32 = 100;
    let mut indices: Vec<u32> = (0..mol.num_atoms()).collect();
    let (chi_nbrs, n) = count_chiral_nbrs(mol, NO_NBRS);
    if chi_nbrs {
        libcxx_sort_by(&mut indices, |a, b| n[*a as usize] < n[*b as usize]);
    }
    let mut wedge = WedgeBonds::new();
    for idx in indices {
        if n[idx as usize] > NO_NBRS {
            continue;
        }
        if !is_tetra(mol.atom(idx).chiral_tag()) {
            break;
        }
        if let Some(b) = pick_bond_to_wedge(mol, idx, &n, &wedge, NO_NBRS)? {
            wedge.insert(b, WedgeInfo::Chiral(idx));
        }
    }
    wedge_bonds_from_atropisomers(mol, conf)?;
    Ok(wedge)
}

/// `Atropisomers::wedgeBondsFromAtropisomers`: only the candidate check is
/// ported.
fn wedge_bonds_from_atropisomers(mol: &mut ROMol, _conf: Option<&Conformer>) -> Result<()> {
    if !mol.ring_info().is_sssr_or_better() {
        find_sssr(mol, false, false)?;
    }
    for b in mol.bonds() {
        if b.bond_type() != BondType::Single
            || !matches!(
                b.stereo(),
                BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW
            )
        {
            continue;
        }
        let (bd, ed) = (
            mol.total_degree(b.begin_atom_idx())?,
            mol.total_degree(b.end_atom_idx())?,
        );
        if (2..=3).contains(&bd) && (2..=3).contains(&ed) {
            return Err(not_ported("wedging atropisomers"));
        }
    }
    Ok(())
}

/// `detail::determineBondWedgeState(bond, fromAtomIdx, conf)`
fn determine_bond_wedge_state_from(
    mol: &ROMol,
    bond: u32,
    from_atom: u32,
    conf: Option<&Conformer>,
) -> Result<BondDir> {
    let b = mol.bond(bond);
    if b.bond_type() != BondType::Single {
        return Err(RDKitError::precondition("bad bond order for wedging"));
    }
    let res = b.bond_dir();
    let Some(conf) = conf else {
        return Ok(res);
    };
    let (atom, bond_atom) = if b.begin_atom_idx() == from_atom {
        (b.begin_atom_idx(), b.end_atom_idx())
    } else {
        (b.end_atom_idx(), b.begin_atom_idx())
    };
    let chiral = mol.atom(atom).chiral_tag();
    if !is_tetra(chiral) {
        return Err(RDKitError::invariant("bad chiral type for wedging"));
    }
    let mut center = conf.atom_pos(atom)?;
    let mut tmp = conf.atom_pos(bond_atom)?;
    center.z = 0.0;
    tmp.z = 0.0;
    let Ok(ref_vect) = center.direction_vector(&tmp) else {
        return Ok(res);
    };
    let mut idxs: Vec<i32> = vec![bond as i32];
    let mut angles: Vec<f64> = vec![0.0];
    for nb in mol.atom_bond_indices(atom) {
        if nb == bond {
            continue;
        }
        let o = mol.bond(nb).other_atom_idx(atom)?;
        let mut p = conf.atom_pos(o)?;
        p.z = 0.0;
        let Ok(v) = center.direction_vector(&p) else {
            return Ok(res);
        };
        let mut angle = ref_vect.signed_angle_to(&v);
        if angle < 0.0 {
            angle += 2. * PI;
        }
        let pos = angles
            .iter()
            .position(|&a| a.partial_cmp(&angle) != Some(std::cmp::Ordering::Less))
            .unwrap_or(angles.len());
        angles.insert(pos, angle);
        idxs.insert(pos, nb as i32);
    }
    let mut n_swaps = mol.perturbation_order(atom, &idxs)?;
    if angles.len() == 3 {
        let angle_tol = PI * 1.9 / 180.;
        if angles[2] - angles[1] >= PI - angle_tol {
            n_swaps += 1;
        }
    }
    Ok(if chiral == ChiralType::TetrahedralCCW {
        if n_swaps % 2 == 1 {
            BondDir::BeginDash
        } else {
            BondDir::BeginWedge
        }
    } else if n_swaps % 2 == 1 {
        BondDir::BeginWedge
    } else {
        BondDir::BeginDash
    })
}

/// `detail::determineBondWedgeState(bond, wedgeBonds, conf)`
pub fn determine_bond_wedge_state(
    mol: &ROMol,
    bond: u32,
    wedge: &WedgeBonds,
    conf: Option<&Conformer>,
) -> Result<BondDir> {
    match wedge.get(&bond) {
        None => Ok(mol.bond(bond).bond_dir()),
        Some(WedgeInfo::Atropisomer(d)) => Ok(*d),
        Some(WedgeInfo::Chiral(a)) => determine_bond_wedge_state_from(mol, bond, *a, conf),
    }
}

/// `detail::isBondPotentialStereoBond`
fn is_bond_potential_stereo_bond(mol: &ROMol, bond: u32) -> Result<bool> {
    let b = mol.bond(bond);
    if b.bond_type() != BondType::Double {
        return Ok(false);
    }
    let (ba, ea) = (b.begin_atom_idx(), b.end_atom_idx());
    let (bd, ed) = (mol.total_degree(ba)?, mol.total_degree(ea)?);
    if bd > 1
        && bd < 4
        && ed > 1
        && ed < 4
        && mol.total_num_hs(ba, true)? < 2
        && mol.total_num_hs(ea, true)? < 2
    {
        for ring in mol.ring_info().bond_rings() {
            if ring.len() < 8 && ring.contains(&(bond as i32)) {
                return Ok(false);
            }
        }
        Ok(true)
    } else {
        Ok(false)
    }
}

/// `Chirality::canBeStereoBond`
fn can_be_stereo_bond(mol: &ROMol, bond: u32) -> bool {
    let b = mol.bond(bond);
    if !matches!(b.bond_type(), BondType::Double | BondType::Aromatic) {
        return false;
    }
    let rank_key = if get_use_legacy_stereo_perception() {
        cp::_CIPRank
    } else {
        cp::_ChiralAtomRank
    };
    for atom in [b.begin_atom_idx(), b.end_atom_idx()] {
        let mut ranks: Vec<i32> = Vec::new();
        for nb in mol.atom_bond_indices(atom) {
            if nb == bond {
                continue;
            }
            let nbb = mol.bond(nb);
            if nbb.bond_type() != BondType::Single {
                continue;
            }
            if matches!(nbb.bond_dir(), BondDir::EndUpRight | BondDir::EndDownRight) {
                return false;
            }
            if nbb.bond_dir() == BondDir::Unknown && nbb.begin_atom_idx() == atom {
                return false;
            }
            let other = if nbb.begin_atom_idx() == atom {
                nbb.end_atom_idx()
            } else {
                nbb.begin_atom_idx()
            };
            let rank = mol
                .atom(other)
                .get_prop(rank_key)
                .and_then(RDValue::as_int)
                .unwrap_or(-1);
            if rank >= 0 {
                if ranks.contains(&rank) {
                    return false;
                }
                ranks.push(rank);
            }
        }
    }
    true
}

/// `Chirality::shouldBeACrossedBond`
fn should_be_a_crossed_bond(mol: &ROMol, bond: u32) -> Result<bool> {
    let b = mol.bond(bond);
    let (ba, ea) = (b.begin_atom_idx(), b.end_atom_idx());
    if b.stereo() == BondStereo::StereoAny {
        for atom in [ba, ea] {
            for nb in mol.atom_bonds(atom) {
                if nb.bond_dir() == BondDir::Unknown && nb.begin_atom_idx() == atom {
                    return Ok(false);
                }
            }
        }
        return Ok(true);
    }
    if b.stereo() != BondStereo::StereoNone {
        return Ok(false);
    }
    if !is_bond_potential_stereo_bond(mol, bond)? {
        return Ok(false);
    }
    if b.bond_dir() == BondDir::EitherDouble {
        return Ok(true);
    }
    let deg = |a: u32| mol.atom_degree(a);
    let excess = |a: u32| -> Result<i64> {
        Ok(mol.atom(a).total_valence()? as i64 - mol.total_degree(a)? as i64)
    };
    if deg(ba) > 1
        && deg(ea) > 1
        && excess(ba)? == 1
        && excess(ea)? == 1
        && can_be_stereo_bond(mol, bond)
    {
        return Ok(true);
    }
    Ok(false)
}

/// `BondGetDirCode`
fn bond_get_dir_code(dir: BondDir) -> i32 {
    match dir {
        BondDir::None => 0,
        BondDir::BeginWedge => 1,
        BondDir::BeginDash => 6,
        BondDir::Unknown => 4,
        BondDir::EitherDouble => 3,
        _ => 0,
    }
}

/// `GetMolFileBondStereoInfo`: `(dirCode, reverse)`
fn get_mol_file_bond_stereo_info(
    mol: &ROMol,
    bond: u32,
    wedge: &WedgeBonds,
    conf: Option<&Conformer>,
) -> Result<(i32, bool)> {
    let b = mol.bond(bond);
    let mut reverse = false;
    let mut dir = BondDir::None;
    if crate::graphmol::atropisomers::can_have_direction(b.bond_type()) {
        // canHaveDirection(const Bond&): SINGLE or AROMATIC; determining a
        // wedge state needs a single bond
        dir = if b.bond_type() == BondType::Single {
            determine_bond_wedge_state(mol, bond, wedge, conf)?
        } else {
            match wedge.get(&bond) {
                Some(_) => determine_bond_wedge_state(mol, bond, wedge, conf)?,
                None => b.bond_dir(),
            }
        };
        if matches!(
            dir,
            BondDir::BeginDash | BondDir::BeginWedge | BondDir::Unknown
        ) && let Some(WedgeInfo::Chiral(a)) = wedge.get(&bond)
            && *a != b.begin_atom_idx()
        {
            reverse = true;
        }
    } else if b.bond_type() == BondType::Double && should_be_a_crossed_bond(mol, bond)? {
        dir = BondDir::EitherDouble;
    }
    Ok((bond_get_dir_code(dir), reverse))
}

// ------------------------------------------------------------------ queries

fn query_bond_topology(q: &BondQuery) -> i32 {
    let mut q = q;
    if q.description == "BondAnd" && !q.negate && q.children.len() == 2 {
        let (c1, c2) = (&q.children[0], &q.children[1]);
        if (c1.description == "BondInRing") != (c2.description == "BondInRing") {
            q = if c1.description == "BondInRing" {
                c1
            } else {
                c2
            };
        }
    }
    if q.description == "BondInRing" {
        if q.negate { 2 } else { 1 }
    } else {
        0
    }
}

fn query_bond_symbol(q: &BondQuery) -> i32 {
    if q.description == "BondOrder" || query_bond_topology(q) != 0 {
        return 0;
    }
    let mut q = q;
    if q.description == "BondAnd" && !q.negate && q.children.len() == 2 {
        if q.children[1].description == "BondInRing" {
            q = &q.children[0];
        } else if q.children[0].description == "BondInRing" {
            q = &q.children[1];
        }
    }
    let d = q.description.as_str();
    if d == "BondOr" && !q.negate {
        if q.children.len() == 2 {
            let (c1, c2) = (&q.children[0], &q.children[1]);
            if c1.description == "BondOrder"
                && !c1.negate
                && c2.description == "BondOrder"
                && !c2.negate
            {
                let (mut t1, mut t2) = (c1.val, c2.val);
                if t1 > t2 {
                    std::mem::swap(&mut t1, &mut t2);
                }
                let (s, d2, a) = (
                    BondType::Single as i32,
                    BondType::Double as i32,
                    BondType::Aromatic as i32,
                );
                if t1 == s && t2 == d2 {
                    return 5;
                } else if t1 == s && t2 == a {
                    return 6;
                } else if t1 == d2 && t2 == a {
                    return 7;
                }
            }
        }
        8
    } else if d == "SingleOrAromaticBond" && !q.negate {
        6
    } else if d == "SingleOrDoubleBond" && !q.negate {
        5
    } else if d == "DoubleOrAromaticBond" && !q.negate {
        7
    } else {
        8
    }
}

fn is_atom_r_group(mol: &ROMol, idx: u32) -> bool {
    let a = mol.atom(idx);
    a.atomic_num() == 0 && a.has_prop(cp::_MolFileRLabel)
}

/// `hasComplexQuery`
fn has_complex_query(mol: &ROMol, idx: u32) -> bool {
    let a = mol.atom(idx);
    let Some(q) = a.query() else {
        return false;
    };
    if q.description == "AtomAtomicNum" && q.val == a.atomic_num() as i32 {
        return false;
    }
    if q.description == "AtomAnd"
        && q.children
            .first()
            .is_some_and(|c| c.description == "AtomAtomicNum")
    {
        return false;
    }
    true
}

fn is_list_query(mol: &ROMol, idx: u32) -> bool {
    let a = mol.atom(idx);
    is_atom_list_query(a.query(), a.atomic_num())
}

/// Whether `moveAdditionalPropertiesToSGroups` would add SubstanceGroups.
fn needs_additional_sgroups(mol: &ROMol) -> bool {
    let redundant = |q: &AtomQuery| {
        q.description == "AtomAnd"
            && q.children.len() == 2
            && q.children[0].description == "AtomAtomicNum"
            && !q.children[0].negate
            && !q.children[1].negate
            && matches!(
                q.children[1].description.as_str(),
                "AtomIsotope" | "AtomFormalCharge"
            )
    };
    for a in mol.atoms() {
        if a.has_prop(cp::_QueryAtomGenericLabel) {
            return true;
        }
        if let Some(q) = a.query() {
            if a.has_prop(cp::MRV_SMA) {
                return true;
            }
            if !is_list_query(mol, a.idx())
                && q.description != "AtomNull"
                && !q.description.starts_with("AtomType")
                && !q.description.starts_with("AtomAtomicNum")
                && !redundant(q)
            {
                return true;
            }
        }
    }
    mol.bonds().iter().any(|b| b.bond_type() == BondType::Zero)
}

// ------------------------------------------------------------------ atoms

/// `AtomGetMolFileSymbol`
fn atom_get_mol_file_symbol(
    mol: &ROMol,
    idx: u32,
    pad: bool,
    query_list_atoms: &mut [bool],
) -> String {
    let a = mol.atom(idx);
    let mut res = if a.has_prop(cp::_MolFileRLabel) {
        "R#".to_string()
    } else if a.atomic_num() != 0 {
        a.symbol()
    } else if !a.has_prop(cp::dummyLabel) {
        let q = a.query();
        let label = q.map_or("", |q| q.type_label.as_str());
        let is_a = q.is_some_and(|q| {
            q.type_label == "A" || (q.negate && q.description == "AtomAtomicNum" && q.val == 1)
        });
        let is_q = q.is_some_and(|q| {
            q.type_label == "Q"
                || (q.negate
                    && q.description == "AtomOr"
                    && q.children.len() == 2
                    && q.children[0].description == "AtomAtomicNum"
                    && q.children[0].val == 6
                    && q.children[1].description == "AtomAtomicNum"
                    && q.children[1].val == 1)
        });
        if is_a {
            query_list_atoms[idx as usize] = true;
            "A".into()
        } else if is_q {
            query_list_atoms[idx as usize] = true;
            "Q".into()
        } else if q.is_some() && matches!(label, "X" | "M" | "AH" | "QH" | "XH" | "MH") {
            query_list_atoms[idx as usize] = true;
            label.into()
        } else if has_complex_query(mol, idx) {
            if is_list_query(mol, idx) {
                "L".into()
            } else {
                "*".into()
            }
        } else {
            "R".into()
        }
    } else {
        let symb = str_prop(a, cp::dummyLabel).unwrap_or_default();
        match symb.as_str() {
            "*" | "X" => "R".into(),
            "Xa" => "R1".into(),
            "Xb" => "R2".into(),
            "Xc" => "R3".into(),
            "Xd" => "R4".into(),
            "Xf" => "R5".into(),
            "Xg" => "R6".into(),
            "Xh" => "R7".into(),
            "Xi" => "R8".into(),
            "Xj" => "R9".into(),
            _ => symb,
        }
    };
    if pad {
        while res.len() < 3 {
            res.push(' ');
        }
    }
    res
}

/// `getAtomParityFlag`
fn atom_parity_flag(mol: &ROMol, idx: u32, conf: &Conformer) -> Result<u32> {
    if !conf.is_3d() || !(mol.atom_degree(idx) >= 3 && mol.total_degree(idx)? == 4) {
        return Ok(0);
    }
    let pos = conf.atom_pos(idx)?;
    let mut vs: Vec<(u32, Point3D)> = Vec::new();
    for nb in mol.atom_neighbors(idx) {
        let v = conf.atom_pos(nb)? - pos;
        let key = if mol.atom(nb).atomic_num() == 1 {
            nb + mol.num_atoms()
        } else {
            nb
        };
        vs.push((key, v));
    }
    vs.sort_by_key(|p| p.0);
    let vol = if vs.len() == 4 {
        vs[0].1.cross_product(&vs[1].1).dot_product(&vs[3].1)
    } else {
        -vs[0].1.cross_product(&vs[1].1).dot_product(&vs[2].1)
    };
    Ok(if vol < 0.0 {
        2
    } else if vol > 0.0 {
        1
    } else {
        0
    })
}

/// `hasNonDefaultValence`
fn has_non_default_valence(mol: &ROMol, idx: u32) -> Result<bool> {
    let a = mol.atom(idx);
    if a.num_radical_electrons() != 0 {
        return Ok(true);
    }
    if a.has_query() || a.needs_update_property_cache() {
        return Ok(false);
    }
    if a.atomic_num() == 1 || in_organic_subset(a.atomic_num() as i32) {
        let eff = a.atomic_num() as i32 - a.formal_charge();
        return Ok(a.no_implicit()
            && a.valence(ValenceType::Explicit)? as i32
                != PeriodicTable::get().get_default_valence(eff.max(0) as u32)?);
    }
    Ok(true)
}

struct AtomProps {
    tot_valence: i32,
    atom_map_number: i32,
    parity: u32,
    pos: Point3D,
}

/// `GetMolFileAtomProperties`
fn atom_properties(mol: &ROMol, idx: u32, conf: Option<&Conformer>) -> Result<AtomProps> {
    let a = mol.atom(idx);
    let mut r = AtomProps {
        tot_valence: 0,
        atom_map_number: int_prop(a, cp::molAtomMapNumber).unwrap_or(0),
        parity: 0,
        pos: Point3D::default(),
    };
    if let Some(conf) = conf {
        r.pos = conf.atom_pos(idx)?;
        if conf.is_3d()
            && !matches!(a.chiral_tag(), ChiralType::Unspecified | ChiralType::Other)
            && mol.atom_degree(idx) >= 3
            && mol.total_degree(idx)? == 4
        {
            r.parity = atom_parity_flag(mol, idx, conf)?;
        }
    }
    if has_non_default_valence(mol, idx)? {
        r.tot_valence = if mol.total_degree(idx)? == 0 {
            15
        } else {
            (a.total_valence()? % 15) as i32
        };
    }
    Ok(r)
}

fn out_of_v2000_range(p: &Point3D) -> bool {
    [p.x, p.y, p.z]
        .iter()
        .any(|&c| c >= MAX_V2000_COORD || c <= MIN_V2000_COORD)
}

/// `GetMolFileAtomLine`
fn mol_file_atom_line(
    mol: &ROMol,
    idx: u32,
    conf: Option<&Conformer>,
    qla: &mut [bool],
) -> Result<String> {
    let p = atom_properties(mol, idx, conf)?;
    if out_of_v2000_range(&p.pos) {
        return Err(RDKitError::Value(
            "MolFile coordinates must be in (-100000, 1000000)".into(),
        ));
    }
    let a = mol.atom(idx);
    let rxn_type = int_prop(a, cp::molRxnRole).unwrap_or(0);
    let rxn_num = int_prop(a, cp::molRxnComponent).unwrap_or(0);
    let symbol = atom_get_mol_file_symbol(mol, idx, true, qla);
    let mut s = format!(
        "{:10.4}{:10.4}{:10.4} {:>3}{:2}{:3}{:3}{:3}{:3}{:3}  0{:3}{:3}{:3}{:3}{:3}",
        p.pos.x,
        p.pos.y,
        p.pos.z,
        symbol,
        0,
        0,
        p.parity,
        0,
        0,
        p.tot_valence,
        rxn_type,
        rxn_num,
        p.atom_map_number,
        0,
        0
    );
    s.truncate(127);
    Ok(s)
}

/// `BondGetMolFileSymbol` / `GetV3000BondCode`
fn bond_symbol(mol: &ROMol, bond: u32, v3000: bool) -> i32 {
    let b = mol.bond(bond);
    let mut res = b.query().map_or(0, query_bond_symbol);
    if res == 0 {
        res = match b.bond_type() {
            BondType::Single | BondType::Double if b.is_aromatic() => 4,
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Aromatic => 4,
            BondType::Zero => 1,
            BondType::Dative => 9,
            BondType::Hydrogen if v3000 => 10,
            _ => 0,
        };
    }
    res
}

/// `GetMolFileBondLine`
fn mol_file_bond_line(
    mol: &ROMol,
    bond: u32,
    wedge: &WedgeBonds,
    conf: Option<&Conformer>,
    was_aromatic: bool,
) -> Result<String> {
    let (mut dir, reverse) = get_mol_file_bond_stereo_info(mol, bond, wedge, conf)?;
    if was_aromatic && dir == 3 {
        dir = 0;
    }
    let b = mol.bond(bond);
    let (a1, a2) = if reverse {
        (b.end_atom_idx() + 1, b.begin_atom_idx() + 1)
    } else {
        (b.begin_atom_idx() + 1, b.end_atom_idx() + 1)
    };
    let mut s = format!("{a1:3}{a2:3}{:3} {dir:2}", bond_symbol(mol, bond, false));
    if let Some(q) = b.query() {
        let topol = query_bond_topology(q);
        if topol != 0 {
            s += &format!(" {:2} {topol:2}", 0);
        }
    }
    Ok(s)
}

/// `M  CHG` / `M  RAD` / `M  ISO` lines, 8 entries per line
fn charge_info(mol: &ROMol) -> Result<String> {
    let mut res = String::new();
    let (mut chg, mut rad, mut iso) = (String::new(), String::new(), String::new());
    let (mut nc, mut nr, mut ni) = (0, 0, 0);
    for a in mol.atoms() {
        let i = a.idx() + 1;
        if a.formal_charge() != 0 {
            nc += 1;
            chg += &format!(" {i:3} {:3}", a.formal_charge());
            if nc == 8 {
                res += &format!("M  CHG{nc:3}{chg}\n");
                chg.clear();
                nc = 0;
            }
        }
        let n_rad = a.num_radical_electrons();
        if n_rad != 0 && mol.total_degree(a.idx())? != 0 {
            nr += 1;
            let code = if n_rad % 2 == 1 { 2 } else { 3 };
            rad += &format!(" {i:3} {code:3}");
            if nr == 8 {
                res += &format!("M  RAD{nr:3}{rad}\n");
                rad.clear();
                nr = 0;
            }
        }
        if !is_atom_r_group(mol, a.idx()) && a.isotope() != 0 {
            ni += 1;
            iso += &format!(" {i:3} {:3}", a.isotope());
            if ni == 8 {
                res += &format!("M  ISO{ni:3}{iso}\n");
                iso.clear();
                ni = 0;
            }
        }
    }
    if nc > 0 {
        res += &format!("M  CHG{nc:3}{chg}\n");
    }
    if nr > 0 {
        res += &format!("M  RAD{nr:3}{rad}\n");
    }
    if ni > 0 {
        res += &format!("M  ISO{ni:3}{iso}\n");
    }
    Ok(res)
}

/// `GetMolFileRGroupInfo`
fn rgroup_info(mol: &ROMol) -> String {
    let mut ss = String::new();
    let mut n = 0;
    for a in mol.atoms() {
        if let Some(lbl) = uint_prop(a, cp::_MolFileRLabel) {
            ss += &format!(" {:3} {lbl:3}", a.idx() + 1);
            n += 1;
        }
    }
    if n > 0 {
        format!("M  RGP{n:3}{ss}\n")
    } else {
        String::new()
    }
}

/// `GetMolFileQueryInfo`
fn query_info(mol: &mut ROMol, qla: &[bool]) -> Result<String> {
    let n = mol.num_atoms();
    let list_qs: Vec<bool> = (0..n)
        .map(|i| is_list_query(mol, i) && !qla[i as usize])
        .collect();
    let mut ss = String::new();
    for i in 0..n {
        let mut wrote = false;
        if !list_qs[i as usize] && !qla[i as usize] && has_complex_query(mol, i) {
            let sma = crate::graphmol::smiles_parse::smarts_write::get_atom_smarts(
                mol,
                i,
                &Default::default(),
            )?;
            ss += &format!("V  {:3} {sma}\n", i + 1);
            wrote = true;
        }
        if !wrote && let Some(v) = str_prop(mol.atom(i), cp::molFileValue) {
            ss += &format!("V  {:3} {v}\n", i + 1);
        }
    }
    for i in 0..n {
        if list_qs[i as usize] {
            let q = mol.atom(i).query().expect("list query");
            let mut vals = Vec::new();
            get_atom_list_query_vals(q, &mut vals)?;
            ss += &format!(
                "M  ALS {:3} {:2}{}",
                i + 1,
                vals.len(),
                if q.negate { " T " } else { " F " }
            );
            for v in vals {
                let sym = PeriodicTable::get().get_element_symbol(v as u32)?;
                ss += &format!("{sym:<4}");
            }
            ss.push('\n');
        }
    }
    Ok(ss)
}

/// `GetMolFileAliasInfo`
fn alias_info(mol: &ROMol) -> String {
    let mut ss = String::new();
    for a in mol.atoms() {
        if let Some(lbl) = str_prop(a, cp::molFileAlias)
            && !lbl.is_empty()
        {
            ss += &format!("A  {:3}\n{lbl}\n", a.idx() + 1);
        }
    }
    ss
}

/// `GetMolFilePXAInfo`
fn pxa_info(mol: &ROMol) -> String {
    let mut res = String::new();
    for a in mol.atoms() {
        if let Some(v) = str_prop(a, "_MolFile_PXA") {
            res += &format!("M  PXA {:>3}{v}\n", format!(" {}", a.idx() + 1));
        }
    }
    res
}

// ------------------------------------------------------------------ V3000

/// `GetV3000MolFileAtomLine`
fn v3000_atom_line(
    mol: &ROMol,
    idx: u32,
    conf: Option<&Conformer>,
    qla: &mut [bool],
    precision: u32,
) -> Result<String> {
    let p = atom_properties(mol, idx, conf)?;
    let a = mol.atom(idx);
    let mut ss = format!("M  V30 {}", idx + 1);
    let symbol = atom_get_mol_file_symbol(mol, idx, false, qla);
    if !is_list_query(mol, idx) || qla[idx as usize] {
        ss += &format!(" {symbol}");
    } else {
        let q = a.query().expect("list query");
        let mut vals = Vec::new();
        get_atom_list_query_vals(q, &mut vals)?;
        if q.negate {
            ss += " \"NOT";
        }
        ss += " [";
        for (i, v) in vals.iter().enumerate() {
            if i != 0 {
                ss.push(',');
            }
            ss += PeriodicTable::get().get_element_symbol(*v as u32)?;
        }
        ss.push(']');
        if q.negate {
            ss.push('"');
        }
    }
    let pr = precision as usize;
    ss += &format!(
        " {:.pr$} {:.pr$} {:.pr$} {}",
        p.pos.x, p.pos.y, p.pos.z, p.atom_map_number
    );
    if p.parity != 0 {
        ss += &format!(" CFG={}", p.parity);
    }
    if a.formal_charge() != 0 {
        ss += &format!(" CHG={}", a.formal_charge());
    }
    if a.isotope() != 0 && !is_atom_r_group(mol, idx) {
        let mut mass = round_half_away(a.mass()) as i32;
        if mass == 0 {
            mass = a.isotope() as i32;
        }
        ss += &format!(" MASS={mass}");
    }
    let n_rad = a.num_radical_electrons();
    if n_rad != 0 && mol.total_degree(idx)? != 0 {
        ss += &format!(" RAD={}", if n_rad % 2 == 1 { 2 } else { 3 });
    }
    if p.tot_valence != 0 {
        if p.tot_valence == 15 {
            ss += " VAL=-1";
        } else {
            ss += &format!(" VAL={}", p.tot_valence);
        }
    }
    if symbol == "R#" {
        let lbl = uint_prop(a, cp::_MolFileRLabel).unwrap_or(1);
        ss += &format!(" RGROUPS=(1 {lbl})");
    }
    let nonzero = |k: &str| int_prop(a, k).filter(|v| *v != 0);
    if let Some(v) = nonzero(cp::molAttachOrder) {
        ss += &format!(" ATTCHORD={v}");
    }
    if let Some(v) = nonzero(cp::molAttachPoint) {
        ss += &format!(" ATTCHPT={v}");
    }
    if let Some(v) = nonzero(cp::molAtomSeqId) {
        ss += &format!(" SEQID={v}");
    }
    if let Some(v) = str_prop(a, cp::molAtomSeqName) {
        ss += &format!(" SEQNAME={v}");
    }
    if let Some(v) = nonzero(cp::molRxnExactChange) {
        ss += &format!(" EXACHG={v}");
    }
    if let Some(v) = nonzero(cp::molInversionFlag)
        && (v == 1 || v == 2)
    {
        ss += &format!(" INVRET={v}");
    }
    if let Some(v) = nonzero(cp::molStereoCare) {
        ss += &format!(" STBOX={v}");
    }
    if let Some(v) = nonzero(cp::molSubstCount) {
        ss += &format!(" SUBST={v}");
    }
    if let Some(v) = nonzero(cp::molRingBondCount) {
        ss += &format!(" RBCNT={v}");
    }
    if let Some(v) = str_prop(a, cp::molAtomClass) {
        ss += &format!(" CLASS={v}");
    }
    Ok(ss)
}

/// `std::round`
fn round_half_away(x: f64) -> f64 {
    x.round()
}

/// `BondStereoCodeV2000ToV3000`
fn v2000_to_v3000_stereo(dir: i32) -> i32 {
    match dir {
        1 => 1,
        3 | 4 => 2,
        6 => 3,
        _ => 0,
    }
}

/// `GetV3000MolFileBondLine`
fn v3000_bond_line(
    mol: &ROMol,
    bond: u32,
    wedge: &WedgeBonds,
    conf: Option<&Conformer>,
    was_aromatic: bool,
) -> Result<String> {
    let (mut dir, reverse) = get_mol_file_bond_stereo_info(mol, bond, wedge, conf)?;
    if was_aromatic && dir == 3 {
        dir = 0;
    }
    let b = mol.bond(bond);
    let mut ss = format!("M  V30 {} {}", bond + 1, bond_symbol(mol, bond, true));
    if reverse {
        ss += &format!(" {} {}", b.end_atom_idx() + 1, b.begin_atom_idx() + 1);
    } else {
        ss += &format!(" {} {}", b.begin_atom_idx() + 1, b.end_atom_idx() + 1);
    }
    if dir != 0 {
        ss += &format!(" CFG={}", v2000_to_v3000_stereo(dir));
    }
    if let Some(q) = b.query() {
        let topol = query_bond_topology(q);
        if topol != 0 {
            ss += &format!(" TOPO={topol}");
        }
    }
    if let Some(v) = int_prop(b, cp::molReactStatus).filter(|v| *v != 0) {
        ss += &format!(" RXCTR={v}");
    }
    for (k, tag) in [
        (cp::molStereoCare, "STBOX"),
        (cp::_MolFileBondEndPts, "ENDPTS"),
        (cp::_MolFileBondAttach, "ATTACH"),
    ] {
        if let Some(v) = str_prop(b, k)
            && v != "0"
        {
            ss += &format!(" {tag}={v}");
        }
    }
    Ok(ss)
}

/// `assignStereoGroupIds`
pub fn assign_stereo_group_ids(groups: &mut [StereoGroup]) {
    fn store(ids: &mut Vec<bool>, sg: &mut StereoGroup) {
        let id = sg.write_id() as usize;
        if id == 0 {
            return;
        }
        if id >= ids.len() {
            ids.resize(id + 1, false);
        }
        if ids[id] {
            crate::rd_log_warning!(
                "StereoGroup ID {id} is used by more than one group, and will be reassined"
            );
            sg.set_write_id(0);
        } else {
            ids[id] = true;
        }
    }
    fn assign(ids: &[bool], next: &mut u32, sg: &mut StereoGroup) {
        if sg.write_id() == 0 {
            *next += 1;
            while (*next as usize) < ids.len() && ids[*next as usize] {
                *next += 1;
            }
            sg.set_write_id(*next);
        }
    }
    let (mut and_ids, mut or_ids) = (Vec::new(), Vec::new());
    for sg in groups.iter_mut() {
        match sg.group_type() {
            StereoGroupType::And => store(&mut and_ids, sg),
            StereoGroupType::Or => store(&mut or_ids, sg),
            StereoGroupType::Absolute => {}
        }
    }
    let (mut and_id, mut or_id) = (0, 0);
    for sg in groups.iter_mut() {
        match sg.group_type() {
            StereoGroupType::And => assign(&and_ids, &mut and_id, sg),
            StereoGroupType::Or => assign(&or_ids, &mut or_id, sg),
            StereoGroupType::Absolute => {}
        }
    }
}

/// `Atropisomers::getAllAtomIdsForStereoGroup`
fn all_atom_ids_for_stereo_group(mol: &ROMol, g: &StereoGroup, wedge: &WedgeBonds) -> Vec<u32> {
    let mut ids: Vec<u32> = g.atoms().to_vec();
    for &bond in g.bonds() {
        let b = mol.bond(bond);
        for atom in [b.begin_atom_idx(), b.end_atom_idx()] {
            for ab in mol.atom_bond_indices(atom) {
                if ab == bond {
                    continue;
                }
                let d = mol.bond(ab).bond_dir();
                if (matches!(d, BondDir::BeginWedge | BondDir::BeginDash)
                    || matches!(wedge.get(&ab), Some(WedgeInfo::Atropisomer(_))))
                    && !ids.contains(&atom)
                {
                    ids.push(atom);
                }
            }
        }
    }
    ids
}

/// `appendEnhancedStereoGroups`
fn enhanced_stereo_groups(mol: &ROMol, wedge: &WedgeBonds) -> String {
    let mut res = String::new();
    if mol.stereo_groups().is_empty() {
        return res;
    }
    let mut groups = mol.stereo_groups().to_vec();
    assign_stereo_group_ids(&mut groups);
    res += "M  V30 BEGIN COLLECTION\n";
    let mut tmp = String::new();
    for g in &groups {
        tmp += "M  V30 MDLV30/";
        match g.group_type() {
            StereoGroupType::Absolute => tmp += "STEABS",
            StereoGroupType::Or => tmp += &format!("STEREL{}", g.write_id()),
            StereoGroupType::And => tmp += &format!("STERAC{}", g.write_id()),
        }
        tmp += " ATOMS=(";
        let ids = all_atom_ids_for_stereo_group(mol, g, wedge);
        tmp += &ids.len().to_string();
        for a in ids {
            tmp.push(' ');
            let s = (a + 1).to_string();
            if tmp.len() + s.len() >= 78 {
                res += &tmp;
                res += "-\n";
                tmp = "M  V30 ".into();
            }
            tmp += &s;
        }
        res += &tmp;
        res += ")\n";
        tmp.clear();
    }
    res += &tmp;
    res += "M  V30 END COLLECTION\n";
    res
}

/// `FileParserUtils::getV3000CTAB`
pub fn get_v3000_ctab(
    mol: &mut ROMol,
    was_aromatic: &[bool],
    conf_id: i32,
    precision: u32,
) -> Result<String> {
    let chiral_flag = uint_prop(mol, cp::_MolFileChiralFlag).unwrap_or(0);
    let conf = if conf_id >= 0 || mol.num_conformers() > 0 {
        Some(mol.conformer(conf_id)?.clone())
    } else {
        None
    };
    let mut res = String::from("M  V30 BEGIN CTAB\n");
    res += &format!(
        "M  V30 COUNTS {} {} {} {} {}\n",
        mol.num_atoms(),
        mol.num_bonds(),
        0,
        0,
        chiral_flag
    );
    let mut qla = vec![false; mol.num_atoms() as usize];
    res += "M  V30 BEGIN ATOM\n";
    for i in 0..mol.num_atoms() {
        res += &v3000_atom_line(mol, i, conf.as_ref(), &mut qla, precision)?;
        res.push('\n');
    }
    res += "M  V30 END ATOM\n";
    let wedge = pick_bonds_to_wedge(mol, conf.as_ref())?;
    if mol.num_bonds() > 0 {
        res += "M  V30 BEGIN BOND\n";
        for b in 0..mol.num_bonds() {
            res += &v3000_bond_line(mol, b, &wedge, conf.as_ref(), was_aromatic[b as usize])?;
            res.push('\n');
        }
        res += "M  V30 END BOND\n";
    }
    if let Some(v) = str_prop(mol, cp::molFileLinkNodes) {
        for ln in v.split('|') {
            res += &format!("M  V30 LINKNODE {ln}\n");
        }
    }
    res += &enhanced_stereo_groups(mol, &wedge);
    res += "M  V30 END CTAB\n";
    Ok(res)
}

#[derive(Clone, Copy, PartialEq)]
enum Format {
    V2000,
    V3000,
    Unspecified,
}

/// `outputMolToMolBlock`
fn output_mol_block(
    mol: &mut RWMol,
    conf_id: i32,
    format: Format,
    precision: u32,
    aromatic_bonds: &[bool],
) -> Result<String> {
    let (n_atoms, n_bonds) = (mol.num_atoms(), mol.num_bonds());
    if format == Format::V2000 && (n_atoms > 999 || n_bonds > 999) {
        return Err(RDKitError::Value(
            "V2000 format does not support more than 999 atoms, bonds or SGroups.".into(),
        ));
    }
    let chiral_flag = uint_prop(mol, cp::_MolFileChiralFlag).unwrap_or(0);
    let conf = if conf_id < 0 && mol.num_conformers() == 0 {
        None
    } else {
        Some(mol.conformer(conf_id)?.clone())
    };
    let too_large = conf
        .as_ref()
        .is_some_and(|c| c.positions().iter().any(out_of_v2000_range));
    if format == Format::V2000 && too_large {
        return Err(RDKitError::Value(format!(
            "V2000 format does not support atom positions <= {} or >= {}",
            MIN_V2000_COORD as i32, MAX_V2000_COORD as i32
        )));
    }
    let mut res = String::new();
    if let Some(n) = str_prop(mol, cp::_Name) {
        res += &n;
    }
    res.push('\n');
    if let Some(info) = str_prop(mol, cp::MolFileInfo) {
        res += &info;
    } else {
        res += &format!("  {:>8}{:>10}", "RDKit", "");
        if let Some(c) = &conf {
            res += if c.is_3d() { "3D" } else { cp::TWOD };
        }
    }
    res.push('\n');
    if let Some(c) = str_prop(mol, cp::MolFileComments) {
        res += &c;
    }
    res.push('\n');

    let has_dative = mol
        .bonds()
        .iter()
        .any(|b| b.bond_type() == BondType::Dative);
    let is_v3000 = format == Format::V3000
        || (format == Format::Unspecified
            && (too_large
                || has_dative
                || n_atoms > 999
                || n_bonds > 999
                || !mol.stereo_groups().is_empty()));
    if is_v3000 {
        res += &format!("{:>3}", 0).repeat(10);
        res += "999 V3000\n";
    } else {
        for v in [n_atoms, n_bonds, 0, 0, chiral_flag, 0, 0, 0, 0, 0] {
            res += &format!("{v:3}");
        }
        res += "999 V2000\n";
    }
    if !is_v3000 {
        let mut qla = vec![false; n_atoms as usize];
        for i in 0..n_atoms {
            res += &mol_file_atom_line(mol, i, conf.as_ref(), &mut qla)?;
            res.push('\n');
        }
        let wedge = pick_bonds_to_wedge(mol, conf.as_ref())?;
        for b in 0..n_bonds {
            res += &mol_file_bond_line(mol, b, &wedge, conf.as_ref(), aromatic_bonds[b as usize])?;
            res.push('\n');
        }
        res += &charge_info(mol)?;
        res += &rgroup_info(mol);
        res += &query_info(mol, &qla)?;
        res += &alias_info(mol);
        // GetMolFileZBOInfo: zero-order bonds are rejected in prepare_mol
        res += &pxa_info(mol);
    } else {
        res += &get_v3000_ctab(mol, aromatic_bonds, conf_id, precision)?;
    }
    res += "M  END\n";
    Ok(res)
}

/// `prepareMol`: returns the bonds that were aromatic before kekulization.
fn prepare_mol(mol: &mut RWMol, params: &MolWriterParams) -> Result<Vec<bool>> {
    if mol.needs_update_property_cache() {
        mol.update_property_cache(false)?;
    }
    let mut aromatic = vec![false; mol.num_bonds() as usize];
    if params.kekulize && mol.num_bonds() > 0 {
        for b in mol.bonds() {
            if b.is_aromatic() {
                aromatic[b.idx() as usize] = true;
            }
        }
        kekulize(mol, true, true, 100)?;
    }
    if params.include_stereo && mol.num_conformers() == 0 {
        return Err(not_ported(
            "computing 2D coordinates (RDDepict::compute2DCoords)",
        ));
    }
    if needs_additional_sgroups(mol) {
        return Err(not_ported("writing SubstanceGroups"));
    }
    Ok(aromatic)
}

/// `MolToMolBlock(mol, params, confId)`
pub fn mol_to_mol_block(mol: &ROMol, params: &MolWriterParams, conf_id: i32) -> Result<String> {
    let mut m = mol.clone();
    let arom = prepare_mol(&mut m, params)?;
    let format = if params.force_v3000 {
        Format::V3000
    } else {
        Format::Unspecified
    };
    output_mol_block(&mut m, conf_id, format, params.precision, &arom)
}

/// `MolToV2KMolBlock(mol, params, confId)`
pub fn mol_to_v2k_mol_block(mol: &ROMol, params: &MolWriterParams, conf_id: i32) -> Result<String> {
    let mut m = mol.clone();
    let arom = prepare_mol(&mut m, params)?;
    output_mol_block(&mut m, conf_id, Format::V2000, params.precision, &arom)
}

/// `MolToV3KMolBlock(mol, params, confId)`
pub fn mol_to_v3k_mol_block(mol: &ROMol, params: &MolWriterParams, conf_id: i32) -> Result<String> {
    let p = MolWriterParams {
        force_v3000: true,
        ..params.clone()
    };
    mol_to_mol_block(mol, &p, conf_id)
}

#[cfg(test)]
mod tests {
    #[test]
    fn printf_style_rounding() {
        // ties round to even like printf("%.4f")
        assert_eq!(format!("{:.4}", 0.03125), "0.0312");
        assert_eq!(format!("{:10.4}", -0.0), "   -0.0000");
        assert_eq!(format!("{:.1}", 0.25), "0.2");
    }
}
