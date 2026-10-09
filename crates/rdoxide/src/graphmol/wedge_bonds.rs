//! Bond wedging: choosing the bonds to wedge around chiral centers and
//! setting the wedge directions from a conformer (`Chirality::pickBondsToWedge`,
//! `wedgeMolBonds`, `wedgeBond`), plus the molfile wedging bookkeeping
//! (`reapplyMolBlockWedging`, `clearMolBlockWedgingInfo`,
//! `invertMolBlockWedgingInfo`).
//!
//! Ported from RDKit `Code/GraphMol/WedgeBonds.cpp` (2026.09.1).
//!
//! Copyright (C) 2023-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::f64::consts::PI;

use super::atom::ChiralType;
use super::atropisomers::can_have_direction;
use super::bond::{BondDir, BondStereo, BondType};
use super::conformer::Conformer;
use super::find_rings::find_sssr;
use super::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

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
    crate::graphmol::atropisomers::wedge_bonds_from_atropisomers(mol, conf, &mut wedge)?;
    Ok(wedge)
}

/// `detail::determineBondWedgeState(bond, fromAtomIdx, conf)`
pub fn determine_bond_wedge_state_from(
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
        return Err(RDKitError::test_assert_at(
            "chiralType == Atom::CHI_TETRAHEDRAL_CW || chiralType == Atom::CHI_TETRAHEDRAL_CCW",
            SOURCE_FILE,
            133,
        ));
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

/// `Chirality::BondWedgingParameters`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BondWedgingParameters {
    /// wedge two bonds at degree-four centers when possible: no ring bonds,
    /// no bonds to chiral atoms and no bonds more than 120 degrees apart
    pub wedge_two_bonds_if_possible: bool,
}

fn int_prop(props: &crate::rdgeneral::Dict, key: &str) -> Option<i32> {
    props.get(key).and_then(RDValue::as_int)
}

/// The conformer to use: the given one or the default conformer of `mol`.
fn conformer_or_default(mol: &ROMol, conf: Option<&Conformer>) -> Result<Conformer> {
    match conf {
        Some(c) => Ok(c.clone()),
        None if mol.num_conformers() > 0 => Ok(mol.conformer(-1)?.clone()),
        None => Err(RDKitError::precondition_at(
            "no conformer available",
            "conf || mol.getNumConformers()",
            SOURCE_FILE,
            473,
        )),
    }
}

const SOURCE_FILE: &str = "Code/GraphMol/WedgeBonds.cpp";

/// `addSecondWedgeAroundAtom`: wedges a second single bond around the begin
/// atom of the wedged `ref_bond`, if there is a suitable one.
fn add_second_wedge_around_atom(mol: &mut ROMol, ref_bond: u32, conf: &Conformer) -> Result<()> {
    let atom = mol.bond(ref_bond).begin_atom_idx();
    // only degree four atoms (per IUPAC recommendation)
    if mol.atom_degree(atom) < 4 {
        return Ok(());
    }
    let flat = |idx: u32| -> Result<_> {
        let mut p = conf.atom_pos(idx)?;
        p.z = 0.0;
        Ok(p)
    };
    let aloc = flat(atom)?;
    let ref_vect = aloc.direction_vector(&flat(mol.bond(ref_bond).end_atom_idx())?)?;
    let mut min_angle = 10000.0;
    let mut best_degree = 100;
    let mut bond_to_wedge = None;
    for bid in mol.atom_bond_indices(atom).collect::<Vec<_>>() {
        let bond = mol.bond(bid);
        let other = bond.other_atom_idx(atom)?;
        if bid == ref_bond
            || bond.bond_type() != BondType::Single
            || bond.bond_dir() != BondDir::None
            || mol.atom(other).chiral_tag() != ChiralType::Unspecified
            || mol.ring_info().num_bond_rings(bid)? > 0
        {
            continue;
        }
        let b_vect = aloc.direction_vector(&flat(other)?)?;
        let angle = ref_vect.angle_to(&b_vect);
        if (angle - min_angle) < 5.0 * PI / 180.0 && mol.atom_degree(other) <= best_degree {
            bond_to_wedge = Some(bid);
            min_angle = angle;
            best_degree = mol.atom_degree(other);
        }
    }
    // only wedge when the angle is below 120 degrees (quasi-arbitrary)
    if let Some(bid) = bond_to_wedge
        && min_angle < 2.0 * PI / 3.0
    {
        let dir = if mol.bond(ref_bond).bond_dir() == BondDir::BeginDash {
            BondDir::BeginWedge
        } else {
            BondDir::BeginDash
        };
        let bond = mol.bond_mut(bid);
        bond.set_bond_dir(dir);
        if bond.begin_atom_idx() != atom {
            let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
            bond.set_begin_atom_idx_unchecked(end);
            bond.set_end_atom_idx_unchecked(begin);
        }
    }
    Ok(())
}

/// `Chirality::wedgeMolBonds`: sets wedge/dash directions on the single bonds
/// around chiral centers (and atropisomer bonds) from the coordinates of
/// `conf`, or of the default conformer.
pub fn wedge_mol_bonds(
    mol: &mut ROMol,
    conf: Option<&Conformer>,
    params: &BondWedgingParameters,
) -> Result<()> {
    let conf = conformer_or_default(mol, conf)?;
    if !mol.ring_info().is_sssr_or_better() {
        find_sssr(mol, false, false)?;
    }
    let wedge_bonds = pick_bonds_to_wedge(mol, Some(&conf))?;
    for (&bid, info) in &wedge_bonds {
        match *info {
            WedgeInfo::Atropisomer(dir) => mol.bond_mut(bid).set_bond_dir(dir),
            WedgeInfo::Chiral(from_atom) => {
                let dir = determine_bond_wedge_state_from(mol, bid, from_atom, Some(&conf))?;
                if !matches!(dir, BondDir::BeginWedge | BondDir::BeginDash) {
                    continue;
                }
                let bond = mol.bond_mut(bid);
                bond.set_bond_dir(dir);
                // the wedge may come from a chiral atom at the end of the
                // bond: reverse the bond so the wedge starts at that atom
                if from_atom != bond.begin_atom_idx() {
                    let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
                    bond.set_begin_atom_idx_unchecked(end);
                    bond.set_end_atom_idx_unchecked(begin);
                }
            }
        }
    }
    if params.wedge_two_bonds_if_possible {
        for atom in 0..mol.num_atoms() {
            if !is_tetra(mol.atom(atom).chiral_tag()) {
                continue;
            }
            let wedged: Vec<u32> = mol
                .atom_bond_indices(atom)
                .filter(|&b| {
                    let bond = mol.bond(b);
                    bond.begin_atom_idx() == atom
                        && bond.bond_type() == BondType::Single
                        && matches!(bond.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash)
                })
                .collect();
            if let [only] = wedged[..] {
                add_second_wedge_around_atom(mol, only, &conf)?;
            }
        }
    }
    Ok(())
}

/// `Chirality::wedgeBond`: wedges `bond` as seen from `from_atom`.
pub fn wedge_bond(mol: &mut ROMol, bond: u32, from_atom: u32, conf: &Conformer) -> Result<()> {
    if mol.bond(bond).bond_type() != BondType::Single {
        return Ok(());
    }
    let dir = determine_bond_wedge_state_from(mol, bond, from_atom, Some(conf))?;
    if matches!(dir, BondDir::BeginWedge | BondDir::BeginDash) {
        mol.bond_mut(bond).set_bond_dir(dir);
    }
    Ok(())
}

/// `Chirality::wedgingHasChirality`: whether a wedge on `bond` can mean
/// something, i.e. the begin atom is chiral or part of an atropisomer.
fn wedging_has_chirality(mol: &ROMol, bond: u32) -> bool {
    let atom = mol.bond(bond).begin_atom_idx();
    if mol.atom(atom).chiral_tag() != ChiralType::Unspecified {
        return true;
    }
    mol.atom_bond_indices(atom).any(|other| {
        let b = mol.bond(other);
        other != bond
            && b.bond_type() == BondType::Single
            && matches!(
                b.stereo(),
                BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW
            )
    })
}

/// `Chirality::reapplyMolBlockWedging`: restores the wedging read from the
/// original mol block, replacing whatever directions are set now. With
/// `verify`, wedges are only applied where they mean something.
pub fn reapply_mol_block_wedging(
    mol: &mut ROMol,
    all_bond_types: bool,
    verify: bool,
) -> Result<()> {
    super::chirality::clear_dir_flags(mol, true);
    for bid in 0..mol.num_bonds() {
        let explicit_unknown = int_prop(&mol.bond(bid).props, cp::_UnknownStereo);
        if explicit_unknown.is_some_and(|v| v != 0) {
            mol.bond_mut(bid).set_bond_dir(BondDir::Unknown);
        }
        let bond = mol.bond(bid);
        // a bond that is neither directional nor part of an atropisomer, or
        // whose begin atom is not chiral, should not carry a wedge
        if verify && (!can_have_direction(bond.bond_type()) || !wedging_has_chirality(mol, bid)) {
            continue;
        }
        let may_wedge = all_bond_types || can_have_direction(bond.bond_type());
        let is_double = bond.bond_type() == BondType::Double;
        let bond_dir = int_prop(&bond.props, cp::_MolFileBondStereo);
        let cfg = int_prop(&bond.props, cp::_MolFileBondCfg);
        let stereo_any = bond.stereo() == BondStereo::StereoAny;
        let can_dir = can_have_direction(bond.bond_type());
        let bond = mol.bond_mut(bid);
        if let Some(code) = bond_dir {
            if may_wedge {
                match code {
                    1 => bond.set_bond_dir(BondDir::BeginWedge),
                    6 => bond.set_bond_dir(BondDir::BeginDash),
                    _ => {}
                }
            }
            if is_double {
                if code == 0 && stereo_any {
                    bond.set_bond_dir(BondDir::None);
                    bond.set_stereo(BondStereo::StereoNone)?;
                } else if code == 3 {
                    bond.set_bond_dir(BondDir::EitherDouble);
                    bond.set_stereo(BondStereo::StereoAny)?;
                }
            }
        }
        match cfg.unwrap_or(-1) {
            1 if may_wedge => bond.set_bond_dir(BondDir::BeginWedge),
            2 if can_dir => bond.set_bond_dir(BondDir::Unknown),
            2 if is_double => {
                bond.set_bond_dir(BondDir::EitherDouble);
                bond.set_stereo(BondStereo::StereoAny)?;
            }
            3 if may_wedge => bond.set_bond_dir(BondDir::BeginDash),
            0 | -1 if bond_dir.is_none() && is_double && stereo_any => {
                bond.set_bond_dir(BondDir::None);
                bond.set_stereo(BondStereo::StereoNone)?;
            }
            _ => {}
        }
    }
    Ok(())
}

/// `Chirality::clearMolBlockWedgingInfo`
pub fn clear_mol_block_wedging_info(mol: &mut ROMol) {
    for bid in 0..mol.num_bonds() {
        let bond = mol.bond_mut(bid);
        bond.clear_prop(cp::_MolFileBondStereo);
        bond.clear_prop(cp::_MolFileBondCfg);
    }
}

/// `Chirality::invertMolBlockWedgingInfo`: swaps wedges and dashes in the
/// stored mol block information.
pub fn invert_mol_block_wedging_info(mol: &mut ROMol) {
    let swap = |code: i32, a: i32, b: i32| match code {
        c if c == a => b,
        c if c == b => a,
        c => c,
    };
    for bid in 0..mol.num_bonds() {
        let bond = mol.bond_mut(bid);
        if let Some(code) = int_prop(&bond.props, cp::_MolFileBondStereo) {
            bond.props.set_val(cp::_MolFileBondStereo, swap(code, 1, 6));
        }
        if let Some(code) = int_prop(&bond.props, cp::_MolFileBondCfg) {
            bond.props.set_val(cp::_MolFileBondCfg, swap(code, 1, 3));
        }
    }
}
