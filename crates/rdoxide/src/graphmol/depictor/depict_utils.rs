//! Helpers for 2D depiction (`RDDepict` utilities).
//!
//! Ported from RDKit `Code/GraphMol/Depictor/DepictUtils.{h,cpp}`
//! (2026.09.1).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::f64::consts::PI;

use crate::geometry::{Point2D, Transform2D};
use crate::graphmol::atom::HybridizationType;
use crate::graphmol::bond::BondStereo;
use crate::graphmol::matrices::get_shortest_path;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::fp::fmuladd;
use crate::{RDKitError, Result};

/// `RDDepict::BOND_LEN` (a global upstream; `Compute2DCoords(bondLength=)`
/// changes it temporarily).
static BOND_LEN_BITS: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0x3FF8_0000_0000_0000);

pub fn bond_len() -> f64 {
    f64::from_bits(BOND_LEN_BITS.load(std::sync::atomic::Ordering::Relaxed))
}

pub fn set_bond_len(v: f64) {
    BOND_LEN_BITS.store(v.to_bits(), std::sync::atomic::Ordering::Relaxed);
}
pub const COLLISION_THRES: f64 = 0.70;
pub const BOND_THRES: f64 = 0.50;
pub const ANGLE_OPEN: f64 = 0.1222;
pub const MAX_COLL_ITERS: u32 = 15;
pub const HETEROATOM_COLL_SCALE: f64 = 1.3;
pub const NUM_BONDS_FLIPS: u32 = 3;
pub const ANGLE_EXPANSION_INCREMENT: f64 = 0.0873;
pub const MAX_ANGLE_EXPANSION: f64 = 3. * ANGLE_EXPANSION_INCREMENT;
pub const MAX_ANGLE_EXPANSION_ITERS: u32 = 20;

/// `RDGeom::INT_POINT2D_MAP`
pub type IntPoint2DMap = BTreeMap<u32, Point2D>;

/// `embedRing`: a regular polygon centred at the origin.
pub fn embed_ring(ring: &[i32]) -> IntPoint2DMap {
    let na = ring.len();
    let ang = 2.0 * PI / na as f64;
    let al = bond_len() / (2.0 * (1.0 - ang.cos())).sqrt();
    let mut res = IntPoint2DMap::new();
    for (i, &a) in ring.iter().enumerate() {
        let x = al * (i as f64 * ang).cos();
        let y = al * (i as f64 * ang).sin();
        res.insert(a as u32, Point2D::new(x, y));
    }
    res
}

/// `computeBisectPoint`
pub fn compute_bisect_point(rcr: &Point2D, ang: f64, nb1: &Point2D, nb2: &Point2D) -> Point2D {
    let mut cloc = *nb1;
    cloc += *nb2;
    cloc *= 0.5;
    if ang > PI {
        cloc -= *rcr;
        cloc *= -1.0;
        cloc += *rcr;
    }
    cloc
}

/// `reflectPoint`: mirror `point` in the line through `loc1` and `loc2`.
pub fn reflect_point(point: &Point2D, loc1: &Point2D, loc2: &Point2D) -> Point2D {
    let org = Point2D::new(0.0, 0.0);
    let xaxis = Point2D::new(1.0, 0.0);
    let mut cent = *loc1 + *loc2;
    cent *= 0.5;
    let mut trans = Transform2D::identity();
    trans.set_transform(&org, &xaxis, &cent, loc1);
    let mut itrans = Transform2D::identity();
    itrans.set_transform(&cent, loc1, &org, &xaxis);
    let mut res = *point;
    trans.transform_point(&mut res);
    res.y = -res.y;
    itrans.transform_point(&mut res);
    res
}

/// `getAtomDepictRank`
pub fn get_atom_depict_rank(mol: &ROMol, aid: u32) -> i32 {
    const MAX_AT_NUM: i32 = 1000;
    const MAX_DEG: i32 = 100;
    let a = mol.atom(aid);
    let anum = if a.atomic_num() == 1 {
        MAX_AT_NUM
    } else {
        a.atomic_num() as i32
    };
    MAX_DEG * anum + mol.atom_degree(aid) as i32
}

/// `rankAtomsByRank(mol, commAtms, ascending)`
pub fn rank_atoms_by_rank(mol: &ROMol, atoms: &[i32], ascending: bool) -> Vec<i32> {
    let na = mol.num_atoms();
    let mut rank_aid: Vec<(i32, i32)> = atoms
        .iter()
        .map(|&aid| {
            let at = mol.atom(aid as u32);
            let rank: u32 = match at.props.get(cp::_CIPRank).and_then(RDValue::as_uint) {
                Some(r) => r,
                None => {
                    let mut r = aid as u32;
                    if let Some(cr) = at.props.get(cp::_ChiralAtomRank).and_then(RDValue::as_uint) {
                        r = na.wrapping_sub(cr);
                    }
                    r.wrapping_add(na.wrapping_mul(get_atom_depict_rank(mol, aid as u32) as u32))
                }
            };
            (rank as i32, aid)
        })
        .collect();
    if ascending {
        rank_aid.sort();
    } else {
        rank_aid.sort_by(|a, b| b.cmp(a));
    }
    rank_aid.into_iter().map(|p| p.1).collect()
}

/// `setNbrOrder`
pub fn set_nbr_order(aid: u32, nbrs: &[i32], mol: &ROMol) -> Result<Vec<i32>> {
    let mut reference = -1;
    for nb in mol.atom_neighbors(aid) {
        if !nbrs.contains(&(nb as i32)) {
            reference = nb as i32;
        }
    }
    let mut thold = nbrs.to_vec();
    if reference >= 0 {
        thold.push(reference);
    }
    if thold.len() <= 3 {
        return Err(RDKitError::invariant("setNbrOrder needs more than 3 atoms"));
    }
    let mut thold = rank_atoms_by_rank(mol, &thold, true);
    let ln = thold.len();
    thold.swap(ln - 3, ln - 2);
    let mut res = Vec::with_capacity(thold.len());
    let pos = thold.iter().position(|&a| a == reference);
    match pos {
        Some(p) => {
            res.extend_from_slice(&thold[p + 1..]);
            res.extend_from_slice(&thold[..p]);
        }
        None => res.extend_from_slice(&thold),
    }
    if res.len() != nbrs.len() {
        return Err(RDKitError::postcondition("setNbrOrder size mismatch"));
    }
    Ok(res)
}

/// `pickFirstRingToEmbed`
pub fn pick_first_ring_to_embed(mol: &ROMol, fused_rings: &[Vec<i32>]) -> usize {
    let mut res = 0usize;
    let mut max_size = 0usize;
    let mut minsubs = 100_000_000i32;
    for (cnt, ring) in fused_rings.iter().enumerate() {
        let subs = ring
            .iter()
            .filter(|&&a| mol.atom_degree(a as u32) > 2)
            .count() as i32;
        if subs < minsubs {
            res = cnt;
            minsubs = subs;
            max_size = ring.len();
        } else if subs == minsubs && ring.len() > max_size {
            res = cnt;
            max_size = ring.len();
        }
    }
    res
}

/// `RDKit::Union(rings, res, exclude)`
pub fn union_rings(rings: &[Vec<i32>], exclude: Option<&[usize]>) -> Vec<i32> {
    let mut res = Vec::new();
    for (id, ring) in rings.iter().enumerate() {
        if exclude.is_some_and(|e| e.contains(&id)) {
            continue;
        }
        for &a in ring {
            if !res.contains(&a) {
                res.push(a);
            }
        }
    }
    res
}

/// `findNextRingToEmbed`: returns `(commonAtoms, nextId)`.
pub fn find_next_ring_to_embed(
    done_rings: &[usize],
    fused_rings: &[Vec<i32>],
) -> Result<(Vec<i32>, usize)> {
    if done_rings.is_empty() || fused_rings.len() <= 1 {
        return Err(RDKitError::precondition("findNextRingToEmbed"));
    }
    let not_done: Vec<usize> = (0..fused_rings.len())
        .filter(|i| !done_rings.contains(i))
        .collect();
    let done_atoms = union_rings(fused_rings, Some(&not_done));
    let mut max_common = 0;
    let mut next_id = 0usize;
    let mut res: Vec<i32> = Vec::new();
    for (cur, ring) in fused_rings.iter().enumerate() {
        if done_rings.contains(&cur) {
            continue;
        }
        let common: Vec<i32> = ring
            .iter()
            .copied()
            .filter(|a| done_atoms.contains(a))
            .collect();
        if common.len() == 2 {
            return Ok((common, cur));
        }
        if common.len() > max_common {
            max_common = common.len();
            next_id = cur;
            res = common;
        }
    }
    let n_cmn = res.len();
    let mut cmn_lst = 0;
    for i in 0..n_cmn {
        if res[i] == fused_rings[next_id][i] {
            cmn_lst += 1;
        } else {
            break;
        }
    }
    if cmn_lst > 0 && cmn_lst < res.len() {
        let temp = res.clone();
        res[..n_cmn - cmn_lst].copy_from_slice(&temp[cmn_lst..n_cmn]);
        let n_mov = n_cmn - cmn_lst;
        res[n_mov..n_mov + cmn_lst].copy_from_slice(&temp[..cmn_lst]);
    }
    if res.is_empty() {
        return Err(RDKitError::postcondition("no next ring found"));
    }
    Ok((res, next_id))
}

/// `getAllRotatableBonds`
pub fn get_all_rotatable_bonds(mol: &ROMol) -> Vec<u32> {
    mol.bonds()
        .iter()
        .filter(|b| {
            b.stereo() as u8 <= BondStereo::StereoAny as u8
                && mol.ring_info().num_bond_rings(b.idx()).unwrap_or(0) == 0
        })
        .map(|b| b.idx())
        .collect()
}

/// `getRotatableBonds(mol, aid1, aid2)`: flippable bonds on the shortest
/// path between the atoms.
pub fn get_rotatable_bonds(mol: &ROMol, aid1: u32, aid2: u32) -> Result<Vec<u32>> {
    let mut path = get_shortest_path(mol, aid1, aid2)?;
    let mut res = Vec::new();
    if path.len() >= 4 {
        path.remove(0);
        path.pop();
        let mut pid = path[0];
        for &aid in &path {
            if aid == pid {
                continue;
            }
            let bid = mol
                .bond_idx_between(pid, aid)
                .ok_or_else(|| RDKitError::invariant("no bond on path"))?;
            let b = mol.bond(bid);
            if b.stereo() as u8 <= BondStereo::StereoAny as u8
                && mol.ring_info().num_bond_rings(bid)? == 0
            {
                res.push(bid);
            }
            pid = aid;
        }
    }
    Ok(res)
}

/// `isSpiroCenter`
pub fn is_spiro_center(aid: u32, mol: &ROMol) -> Result<bool> {
    if mol.atom_degree(aid) != 4 || mol.ring_info().num_atom_rings(aid)? != 2 {
        return Ok(false);
    }
    let rings: Vec<&Vec<i32>> = mol
        .ring_info()
        .atom_rings()
        .iter()
        .filter(|r| r.contains(&(aid as i32)))
        .collect();
    if rings.len() != 2 {
        return Ok(false);
    }
    let n = mol.num_atoms() as usize;
    let mut r1 = vec![false; n];
    let mut r2 = vec![false; n];
    for &i in rings[0] {
        r1[i as usize] = true;
    }
    for &i in rings[1] {
        r2[i as usize] = true;
    }
    if (0..n).filter(|&i| r1[i] && r2[i]).count() != 1 {
        return Ok(false);
    }
    let (mut n1, mut n2) = (0, 0);
    for nb in mol.atom_neighbors(aid) {
        let (a, b) = (r1[nb as usize], r2[nb as usize]);
        if a && !b {
            n1 += 1;
        } else if b && !a {
            n2 += 1;
        } else {
            return Ok(false);
        }
    }
    Ok(n1 == 2 && n2 == 2)
}

/// `computeSubAngle`
pub fn compute_sub_angle(degree: u32, htype: HybridizationType) -> f64 {
    match htype {
        HybridizationType::Unspecified | HybridizationType::SP3 => {
            if degree == 4 {
                PI / 2.0
            } else {
                2.0 * PI / 3.0
            }
        }
        HybridizationType::SP2 => 2.0 * PI / 3.0,
        _ => 2. * PI / degree as f64,
    }
}

/// `rotationDir`
pub fn rotation_dir(center: &Point2D, loc1: &Point2D, loc2: &Point2D, rem_angle: f64) -> i32 {
    let pt1 = *loc1 - *center;
    let pt2 = *loc2 - *center;
    let mut cross = fmuladd(pt1.x, pt2.y, -(pt1.y * pt2.x));
    let diff = PI - rem_angle;
    cross *= diff;
    if cross >= 0.0 { -1 } else { 1 }
}

/// `computeNormal`
pub fn compute_normal(center: &Point2D, other: &Point2D) -> Result<Point2D> {
    let mut res = *other - *center;
    res.normalize()?;
    std::mem::swap(&mut res.x, &mut res.y);
    res.x *= -1.0;
    Ok(res)
}

/// `computeAngle`
pub fn compute_angle(center: &Point2D, loc1: &Point2D, loc2: &Point2D) -> Result<f64> {
    let v1 = *loc1 - *center;
    let v2 = *loc2 - *center;
    v1.angle_to(&v2)
}
