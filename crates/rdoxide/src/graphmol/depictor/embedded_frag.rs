//! Fragments with 2D coordinates that grow into a full depiction
//! (`RDDepict::EmbeddedAtom`, `RDDepict::EmbeddedFrag`).
//!
//! Ported from RDKit `Code/GraphMol/Depictor/EmbeddedFrag.{h,cpp}`
//! (2026.09.1). Random flip sampling is not
//! ported yet. Upstream's `std::map::operator[]` default-inserts missing
//! atoms; [`EmbeddedFrag::ea`] reproduces that where upstream relies on it.
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::f64::consts::PI;
use std::sync::Arc;

use super::depict_utils::*;
use super::templates::CoordinateTemplates;
use crate::geometry::{Point2D, Transform2D};
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::matrices::{get_distance_mat, get_shortest_path};
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::fp::fmuladd;
use crate::rdgeneral::random::with_global_generator;
use crate::{RDKitError, Result};

const NEIGH_RADIUS: f64 = 2.5;

/// `RDDepict::EmbeddedAtom`
#[derive(Debug, Clone, PartialEq)]
pub struct EmbeddedAtom {
    pub aid: u32,
    pub angle: f64,
    pub nbr1: i32,
    pub nbr2: i32,
    pub cis_trans_nbr: i32,
    pub ccw: bool,
    pub rot_dir: i32,
    pub loc: Point2D,
    pub normal: Point2D,
    pub neighs: Vec<i32>,
    pub density: f64,
    pub fixed: bool,
}

impl Default for EmbeddedAtom {
    fn default() -> Self {
        EmbeddedAtom {
            aid: 0,
            angle: -1.0,
            nbr1: -1,
            nbr2: -1,
            cis_trans_nbr: -1,
            ccw: true,
            rot_dir: 0,
            loc: Point2D::default(),
            normal: Point2D::default(),
            neighs: Vec::new(),
            density: -1.0,
            fixed: false,
        }
    }
}

impl EmbeddedAtom {
    /// `EmbeddedAtom(aid, pos)`
    pub fn at(aid: u32, pos: Point2D) -> Self {
        EmbeddedAtom {
            aid,
            loc: pos,
            ..Default::default()
        }
    }

    /// `Transform`
    pub fn transform(&mut self, trans: &Transform2D) {
        let mut temp = self.loc + self.normal;
        trans.transform_point(&mut self.loc);
        trans.transform_point(&mut temp);
        self.normal = temp - self.loc;
    }

    /// `Reflect`
    pub fn reflect(&mut self, loc1: &Point2D, loc2: &Point2D) {
        let temp = self.loc + self.normal;
        self.loc = reflect_point(&self.loc, loc1, loc2);
        let temp = reflect_point(&temp, loc1, loc2);
        self.normal = temp - self.loc;
        self.ccw = !self.ccw;
    }
}

/// `RDDepict::EmbeddedFrag`
#[derive(Debug, Clone)]
pub struct EmbeddedFrag<'a> {
    pub(crate) done: bool,
    px: f64,
    nx: f64,
    py: f64,
    ny: f64,
    pub(crate) eatoms: BTreeMap<u32, EmbeddedAtom>,
    attach_pts: VecDeque<i32>,
    mol: &'a ROMol,
    /// position in upstream's `std::list` (fragments are only appended)
    pub(crate) seq: usize,
}

fn missing(aid: impl std::fmt::Display) -> RDKitError {
    RDKitError::Index(format!("map::at: embedded atom {aid} not found"))
}

/// `_recurseAtomOneSide`
fn recurse_atom_one_side(end_aid: u32, beg_aid: u32, mol: &ROMol, flip: &mut Vec<i32>) {
    flip.push(end_aid as i32);
    for nb in mol.atom_neighbors(end_aid) {
        if nb != beg_aid && !flip.contains(&(nb as i32)) {
            recurse_atom_one_side(nb, beg_aid, mol, flip);
        }
    }
}

fn cross_val(v1: &Point2D, v2: &Point2D) -> f64 {
    fmuladd(v1.x, v2.y, -(v2.x * v1.y))
}

const SEG_EPS: f64 = 1e-6;

fn orientation(a: &Point2D, b: &Point2D, c: &Point2D) -> f64 {
    cross_val(&(*b - *a), &(*c - *a))
}

fn point_on_segment(a: &Point2D, b: &Point2D, p: &Point2D) -> bool {
    orientation(a, b, p).abs() <= SEG_EPS
        && p.x >= a.x.min(b.x) - SEG_EPS
        && p.x <= a.x.max(b.x) + SEG_EPS
        && p.y >= a.y.min(b.y) - SEG_EPS
        && p.y <= a.y.max(b.y) + SEG_EPS
}

/// `std::max` / `std::min` of two doubles
fn cmax(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}
fn cmin(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

fn segments_intersect(a: &Point2D, b: &Point2D, c: &Point2D, d: &Point2D) -> bool {
    if cmax(a.x, b.x) + SEG_EPS < cmin(c.x, d.x)
        || cmax(c.x, d.x) + SEG_EPS < cmin(a.x, b.x)
        || cmax(a.y, b.y) + SEG_EPS < cmin(c.y, d.y)
        || cmax(c.y, d.y) + SEG_EPS < cmin(a.y, b.y)
    {
        return false;
    }
    let o1 = orientation(a, b, c);
    let o2 = orientation(a, b, d);
    let o3 = orientation(c, d, a);
    let o4 = orientation(c, d, b);
    let first = (o1 > SEG_EPS && o2 < -SEG_EPS) || (o1 < -SEG_EPS && o2 > SEG_EPS);
    let second = (o3 > SEG_EPS && o4 < -SEG_EPS) || (o3 < -SEG_EPS && o4 > SEG_EPS);
    if first && second {
        return true;
    }
    point_on_segment(a, b, c)
        || point_on_segment(a, b, d)
        || point_on_segment(c, d, a)
        || point_on_segment(c, d, b)
}

/// `_findFarthestPair`
fn find_farthest_pair(
    beg1: u32,
    end1: u32,
    beg2: u32,
    end2: u32,
    na: usize,
    dmat: &[f64],
) -> (u32, u32) {
    let d = |a: u32, b: u32| dmat[a as usize * na + b as usize];
    let cands = [
        (d(beg1, beg2), (beg1, beg2)),
        (d(beg1, end2), (beg1, end2)),
        (d(end1, beg2), (end1, beg2)),
        (d(end1, end2), (end1, end2)),
    ];
    // std::max(a, b, comp) returns a unless comp(a, b)
    let mut max = cands[0];
    for c in &cands[1..] {
        if max.0 < c.0 {
            max = *c;
        }
    }
    max.1
}

/// `_findDeg1Neighbor`
fn find_deg1_neighbor(mol: &ROMol, aid: u32) -> Result<u32> {
    if mol.atom_degree(aid) != 1 {
        return Err(RDKitError::invariant("degree must be 1"));
    }
    mol.atom_neighbors(aid)
        .next()
        .ok_or_else(|| RDKitError::invariant("no neighbor"))
}

/// `_findClosestNeighbor`
fn find_closest_neighbor(mol: &ROMol, dmat: &[f64], aid1: u32, aid2: u32) -> u32 {
    let mut res = 0;
    let mut mdist = 1.0e8;
    let naid = aid1 as usize * mol.num_atoms() as usize;
    for nb in mol.atom_neighbors(aid2) {
        let d = dmat[naid + nb as usize];
        if d < mdist {
            mdist = d;
            res = nb;
        }
    }
    res
}

/// `_recurseDegTwoRingAtoms`
fn recurse_deg_two_ring_atoms(
    aid: u32,
    mol: &ROMol,
    r_path: &mut Vec<i32>,
    nbr_map: &mut BTreeMap<i32, Vec<i32>>,
) -> Result<()> {
    let mut nbrs = Vec::new();
    for b in mol.atom_bonds(aid) {
        if mol.ring_info().num_bond_rings(b.idx())? > 0 {
            nbrs.push(b.other_atom_idx(aid)? as i32);
        }
    }
    if nbrs.len() != 2 {
        return Ok(());
    }
    r_path.push(aid as i32);
    nbr_map.insert(aid as i32, nbrs.clone());
    for nb in nbrs {
        if !r_path.contains(&nb) {
            recurse_deg_two_ring_atoms(nb as u32, mol, r_path, nbr_map)?;
        }
    }
    Ok(())
}

/// `_anyNonRingBonds`
fn any_non_ring_bonds(aid: u32, path: &[u32], mol: &ROMol) -> Result<u32> {
    let mut prev = aid;
    let mut n_open = 0;
    for &pi in path {
        let b = mol
            .bond_idx_between(prev, pi)
            .ok_or_else(|| RDKitError::invariant("no bond found"))?;
        if mol.ring_info().num_bond_rings(b)? == 0 {
            n_open += 1;
        }
        prev = pi;
    }
    Ok(n_open)
}

/// A template and the way a ring system matches it: (template atom, atom) pairs.
type TemplateMatch = (Arc<ROMol>, Vec<(i32, i32)>);

impl<'a> EmbeddedFrag<'a> {
    fn empty(mol: &'a ROMol) -> Self {
        EmbeddedFrag {
            done: false,
            px: 0.0,
            nx: 0.0,
            py: 0.0,
            ny: 0.0,
            eatoms: BTreeMap::new(),
            attach_pts: VecDeque::new(),
            mol,
            seq: 0,
        }
    }

    /// `d_eatoms[aid]` (default-inserting)
    fn ea(&mut self, aid: impl TryInto<u32>) -> &mut EmbeddedAtom {
        let aid = aid.try_into().unwrap_or(u32::MAX);
        self.eatoms.entry(aid).or_default()
    }

    /// `d_eatoms.at(aid)`
    fn at(&self, aid: impl TryInto<u32> + Copy + std::fmt::Display) -> Result<&EmbeddedAtom> {
        let k = aid.try_into().map_err(|_| missing(aid))?;
        self.eatoms.get(&k).ok_or_else(|| missing(aid))
    }

    fn has(&self, aid: impl TryInto<u32>) -> bool {
        aid.try_into().is_ok_and(|k| self.eatoms.contains_key(&k))
    }

    pub fn embedded_atoms(&self) -> &BTreeMap<u32, EmbeddedAtom> {
        &self.eatoms
    }

    pub fn size(&self) -> usize {
        self.eatoms.len()
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn mark_done(&mut self) {
        self.done = true;
    }

    /// `EmbeddedFrag(aid, mol)`: one atom at the origin.
    pub fn single_atom(aid: u32, mol: &'a ROMol) -> Result<Self> {
        let mut f = EmbeddedFrag::empty(mol);
        let eatm = EmbeddedAtom {
            aid,
            loc: Point2D::new(0.0, 0.0),
            normal: Point2D::new(1.0, 0.0),
            ..Default::default()
        };
        f.eatoms.insert(aid, eatm);
        f.update_new_neighs(aid)?;
        Ok(f)
    }

    /// `EmbeddedFrag(mol, fusedRings, useRingTemplates)`
    pub fn fused_rings(
        mol: &'a ROMol,
        fused: &[Vec<i32>],
        use_ring_templates: bool,
    ) -> Result<Self> {
        let mut f = EmbeddedFrag::empty(mol);
        f.embed_fused_rings(fused, use_ring_templates)?;
        Ok(f)
    }

    /// `EmbeddedFrag(mol, coordMap)`
    pub fn from_coord_map(mol: &'a ROMol, coord_map: &IntPoint2DMap) -> Result<Self> {
        let mut f = EmbeddedFrag::empty(mol);
        let na = mol.num_atoms();
        for (&aid, p) in coord_map {
            if aid >= na {
                return Err(RDKitError::invariant("coordMap atom out of range"));
            }
            let mut ea = EmbeddedAtom::at(aid, *p);
            ea.fixed = true;
            f.eatoms.insert(aid, ea);
        }
        f.setup_new_neighs()?;
        f.setup_attachment_points()?;
        Ok(f)
    }

    /// `EmbeddedFrag(const Bond *dblBond)`
    pub fn double_bond(mol: &'a ROMol, bond: u32) -> Result<Self> {
        let b = mol.bond(bond);
        let stype = b.stereo();
        if b.bond_type() != BondType::Double
            || stype as u8 <= BondStereo::StereoAny as u8
            || b.stereo_atoms().len() != 2
        {
            return Err(RDKitError::precondition("bad cis/trans bond"));
        }
        let nbr = b.stereo_atoms().to_vec();
        let (beg, end) = (b.begin_atom_idx(), b.end_atom_idx());
        let mut f = EmbeddedFrag::empty(mol);
        let beatm = EmbeddedAtom {
            aid: beg,
            loc: Point2D::new(0.0, 0.0),
            nbr1: end as i32,
            normal: Point2D::new(0.0, -1.0),
            ccw: false,
            cis_trans_nbr: nbr[0],
            ..Default::default()
        };
        f.eatoms.insert(beg, beatm);
        let cis = matches!(stype, BondStereo::StereoZ | BondStereo::StereoCis);
        let eeatm = EmbeddedAtom {
            aid: end,
            loc: Point2D::new(bond_len(), 0.0),
            nbr1: beg as i32,
            cis_trans_nbr: nbr[1],
            normal: Point2D::new(0.0, if cis { -1.0 } else { 1.0 }),
            ccw: cis,
            ..Default::default()
        };
        f.eatoms.insert(end, eeatm);
        Ok(f)
    }

    /// `computeNbrsAndAng`
    fn compute_nbrs_and_ang(&mut self, aid: u32, done_nbrs: &[i32]) -> Result<()> {
        if done_nbrs.len() < 3 {
            return Err(RDKitError::precondition("computeNbrsAndAng"));
        }
        let mut dirs: Vec<(f64, i32)> = Vec::with_capacity(done_nbrs.len());
        for &nbr in done_nbrs {
            let delta = self.ea(nbr).loc - self.ea(aid).loc;
            dirs.push((delta.y.atan2(delta.x), nbr));
        }
        dirs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut pairs: Vec<(f64, (i32, i32))> = Vec::with_capacity(dirs.len());
        for i in 0..dirs.len() {
            let start = dirs[i];
            let end = dirs[(i + 1) % dirs.len()];
            let mut angle = end.0 - start.0;
            if i + 1 == dirs.len() {
                angle += 2.0 * PI;
            }
            pairs.push((angle, (start.1, end.1)));
        }
        pairs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut winner = *pairs.last().expect("pairs");
        let ri = self.mol.ring_info();
        for pr in pairs.iter().rev() {
            if ri.num_atom_rings(pr.1.0 as u32)? <= 1 && ri.num_atom_rings(pr.1.1 as u32)? <= 1 {
                winner = *pr;
                break;
            }
        }
        let e = self.ea(aid);
        e.rot_dir = 1;
        e.nbr1 = winner.1.1;
        e.nbr2 = winner.1.0;
        e.angle = 2.0 * PI - winner.0;
        Ok(())
    }

    /// `findNumNeigh`
    fn find_num_neigh(&self, pt: &Point2D, radius: f64) -> i32 {
        self.eatoms
            .values()
            .filter(|e| (e.loc - *pt).length() < radius)
            .count() as i32
    }

    /// `updateNewNeighs`
    fn update_new_neighs(&mut self, aid: u32) -> Result<()> {
        let mol = self.mol;
        self.ea(aid).neighs.clear();
        let mut neighs = Vec::new();
        let mut h_indices = Vec::new();
        for nb in mol.atom_neighbors(aid) {
            if !self.eatoms.contains_key(&nb) {
                if mol.atom(nb).atomic_num() != 1 {
                    neighs.push(nb as i32);
                } else {
                    h_indices.push(nb as i32);
                }
            }
        }
        neighs.extend(h_indices);
        let deg = mol.atom_degree(aid);
        if !neighs.is_empty() && (deg < 4 || neighs.len() < 3) {
            neighs = rank_atoms_by_rank(mol, &neighs, true);
        } else if deg >= 4 && neighs.len() >= 3 {
            neighs = set_nbr_order(aid, &neighs, mol)?;
        }
        let has = !neighs.is_empty();
        self.ea(aid).neighs = neighs;
        if has && !self.attach_pts.contains(&(aid as i32)) {
            self.attach_pts.push_back(aid as i32);
        }
        Ok(())
    }

    /// `setupNewNeighs`
    pub fn setup_new_neighs(&mut self) -> Result<()> {
        self.attach_pts.clear();
        let keys: Vec<u32> = self.eatoms.keys().copied().collect();
        for k in keys {
            self.update_new_neighs(k)?;
        }
        let v: Vec<i32> = self.attach_pts.iter().copied().collect();
        self.attach_pts = rank_atoms_by_rank(self.mol, &v, true).into();
        Ok(())
    }

    /// `setupAttachmentPoints`
    fn setup_attachment_points(&mut self) -> Result<()> {
        let pts: Vec<i32> = self.attach_pts.iter().copied().collect();
        for dai in pts {
            let enbrs = self.ea(dai).neighs.clone();
            let done: Vec<i32> = self
                .mol
                .atom_neighbors(dai as u32)
                .map(|n| n as i32)
                .filter(|n| !enbrs.contains(n))
                .collect();
            match done.len() {
                0 => {
                    let e = self.ea(dai);
                    e.normal = Point2D::new(1., 0.);
                    e.angle = -1.;
                }
                1 => {
                    let nbid = done[0];
                    let (l0, l1) = (self.ea(dai).loc, self.ea(nbid).loc);
                    let n = compute_normal(&l0, &l1)?;
                    let e = self.ea(dai);
                    e.nbr1 = nbid;
                    e.normal = n;
                }
                2 => {
                    let (nb1, nb2) = (done[0], done[1]);
                    let (l0, l1, l2) = (self.ea(dai).loc, self.ea(nb1).loc, self.ea(nb2).loc);
                    let a = compute_angle(&l0, &l1, &l2)?;
                    let e = self.ea(dai);
                    e.nbr1 = nb1;
                    e.nbr2 = nb2;
                    e.angle = a;
                }
                _ => self.compute_nbrs_and_ang(dai as u32, &done)?,
            }
        }
        Ok(())
    }

    /// `mirrorTransRingAtoms`
    fn mirror_trans_ring_atoms(
        mol: &ROMol,
        ring: &[i32],
        coords: &mut IntPoint2DMap,
    ) -> Result<()> {
        let n = ring.len();
        for i in 0..n {
            let atom1 = ring[i];
            let atom2 = ring[(i + 1) % n];
            let bid = mol
                .bond_idx_between(atom1 as u32, atom2 as u32)
                .ok_or_else(|| RDKitError::invariant("ring bond missing"))?;
            let b = mol.bond(bid);
            if b.bond_type() != BondType::Double {
                continue;
            }
            let stype = b.stereo();
            if stype as u8 <= BondStereo::StereoAny as u8 || b.stereo_atoms().len() != 2 {
                continue;
            }
            let nb = b.stereo_atoms();
            let left_in = ring.contains(&nb[0]);
            let right_in = ring.contains(&nb[1]);
            let is_trans = if matches!(stype, BondStereo::StereoTrans | BondStereo::StereoE) {
                left_in == right_in
            } else {
                left_in != right_in
            };
            if !is_trans {
                continue;
            }
            let left = ring[(i + n - 1) % n] as u32;
            let right = atom2 as u32;
            let last = *coords.entry(left).or_default();
            let refp = *coords.entry(right).or_default();
            let interest = *coords.entry(atom1 as u32).or_default();
            let d = last - refp;
            let a = fmuladd(d.x, d.x, -(d.y * d.y)) / d.dot_product(&d);
            let bb = 2.0 * d.x * d.y / d.dot_product(&d);
            let x = fmuladd(a, interest.x - refp.x, bb * (interest.y - refp.y)) + refp.x;
            let y = fmuladd(bb, interest.x - refp.x, -(a * (interest.y - refp.y))) + refp.y;
            coords.insert(atom1 as u32, Point2D::new(x, y));
        }
        Ok(())
    }

    /// `checkStereoChemistry`: whether the double bond stereochemistry of
    /// the ring system matches that of the template, as laid out.
    fn check_stereo_chemistry(
        mol: &ROMol,
        template: &ROMol,
        matched: &[(i32, i32)],
    ) -> Result<bool> {
        for bond in mol.bonds() {
            if bond.bond_type() != BondType::Double
                || bond.stereo() == BondStereo::StereoAny
                || bond.stereo() == BondStereo::StereoNone
            {
                continue;
            }
            // the four atoms around the double bond
            let neighbors = bond.stereo_atoms();
            if neighbors.len() != 2 {
                continue;
            }
            let atom1_neighbor1 = neighbors[0];
            let atom2_neighbor1 = neighbors[1];
            let atom1 = bond.begin_atom_idx() as i32;
            let atom2 = bond.end_atom_idx() as i32;
            // now the other two atoms that are not part of the double bond (if any)
            let other_neighbor = |atom: i32, first: i32, other_end: i32| -> i32 {
                if mol.atom_degree(atom as u32) > 2 {
                    for n in mol.atom_neighbors(atom as u32) {
                        if n as i32 != first && n as i32 != other_end {
                            return n as i32;
                        }
                    }
                }
                -1
            };
            let atom1_neighbor2 = other_neighbor(atom1, atom1_neighbor1, atom2);
            let atom2_neighbor2 = other_neighbor(atom2, atom2_neighbor1, atom1);

            // find the template atoms that correspond to the four atoms
            let (mut t_atom1, mut t_atom2) = (-1, -1);
            let (mut t_a1_n1, mut t_a1_n2, mut t_a2_n1, mut t_a2_n2) = (-1, -1, -1, -1);
            for &(template_aidx, rs_aidx) in matched {
                if rs_aidx == atom1 {
                    t_atom1 = template_aidx;
                } else if rs_aidx == atom2 {
                    t_atom2 = template_aidx;
                } else if rs_aidx == atom1_neighbor1 {
                    t_a1_n1 = template_aidx;
                } else if rs_aidx == atom2_neighbor1 {
                    t_a2_n1 = template_aidx;
                } else if rs_aidx == atom1_neighbor2 {
                    t_a1_n2 = template_aidx;
                } else if rs_aidx == atom2_neighbor2 {
                    t_a2_n2 = template_aidx;
                }
            }
            // the atoms controlling the double bond stereochemistry in the
            // molecule may not be the atoms that matched the template: swap to
            // the other atom
            let mut swap_stereo = false;
            if t_a1_n1 == -1 {
                t_a1_n1 = t_a1_n2;
                swap_stereo = !swap_stereo;
            }
            if t_a2_n1 == -1 {
                t_a2_n1 = t_a2_n2;
                swap_stereo = !swap_stereo;
            }
            if t_atom1 == -1 || t_atom2 == -1 || t_a1_n1 == -1 || t_a2_n1 == -1 {
                return Ok(false);
            }
            let conf = template.conformer(-1)?;
            let atom1_loc = conf.atom_pos(t_atom1 as u32)?;
            let atom2_loc = conf.atom_pos(t_atom2 as u32)?;
            let atom1_neighbor_loc = conf.atom_pos(t_a1_n1 as u32)?;
            let atom2_neighbor_loc = conf.atom_pos(t_a2_n1 as u32)?;
            // are the two neighbors on the same side of the bond?
            let v12 = atom1_neighbor_loc - atom1_loc;
            let v42 = atom2_neighbor_loc - atom1_loc;
            let v32 = atom2_loc - atom1_loc;
            let cross1 = v32.x * v12.y - v32.y * v12.x;
            let cross2 = v32.x * v42.y - v32.y * v42.x;
            let mut is_cis = cross1 * cross2 > 0.0;
            if swap_stereo {
                is_cis = !is_cis;
            }
            if is_cis != matches!(bond.stereo(), BondStereo::StereoZ | BondStereo::StereoCis) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// `matchToTemplate`: lays out a ring system as its template, if it has
    /// one. Only an exact match to the ring system is looked for, because
    /// completing rings from a template isn't reliably better than not using
    /// a template at all.
    fn match_to_template(&mut self, ring_system_atoms: &[i32]) -> Result<bool> {
        const DUMMY_ATOMIC_NUM: u32 = 200;
        let candidates = {
            let templates = CoordinateTemplates::ring_system_templates()
                .read()
                .unwrap_or_else(|e| e.into_inner());
            if !templates.has_template_of_size(ring_system_atoms.len() as u32) {
                return Ok(false);
            }
            templates.matching_templates(ring_system_atoms.len() as u32)
        };

        // the induced subgraph of the ring system: everything else is dummies
        let mut rs_mol = self.mol.clone();
        rs_mol.clear_conformers();
        let in_system = |a: u32| ring_system_atoms.contains(&(a as i32));
        for idx in 0..rs_mol.num_atoms() {
            if !in_system(idx) {
                rs_mol.atom_mut(idx).set_atomic_num(DUMMY_ATOMIC_NUM as i32);
            }
        }
        let num_bonds = rs_mol
            .bonds()
            .iter()
            .filter(|b| in_system(b.begin_atom_idx()) && in_system(b.end_atom_idx()))
            .count() as u32;

        // how many atoms have 0 to 4 neighbors that are not dummies
        let degree_counts = |mol: &ROMol| {
            let mut counts = [0i32; 5];
            for atom in mol.atoms() {
                if atom.atomic_num() == DUMMY_ATOMIC_NUM {
                    continue;
                }
                let degree = mol
                    .atom_neighbors(atom.idx())
                    .filter(|&n| mol.atom(n).atomic_num() != DUMMY_ATOMIC_NUM)
                    .count()
                    .min(4);
                counts[degree] += 1;
            }
            counts
        };
        let rs_degrees = degree_counts(&rs_mol);

        // find the template that this molecule matches, if any
        let params = SubstructMatchParameters {
            max_matches: 1,
            ..SubstructMatchParameters::default()
        };
        let mut found: Option<TemplateMatch> = None;
        for template in candidates {
            // to reduce how often a substructure match is needed, check the
            // bond count and the connectivity first
            if template.num_bonds() != num_bonds || degree_counts(&template) != rs_degrees {
                continue;
            }
            let matches = substruct_match(&rs_mol, &template, &params)?;
            if let Some(first) = matches.first()
                && Self::check_stereo_chemistry(&rs_mol, &template, first)?
            {
                found = Some((template, first.clone()));
                break;
            }
        }
        let Some((template, matched)) = found else {
            return Ok(false);
        };

        // copy over the new coordinates
        let conf = template.conformer(-1)?;
        for &(template_aidx, rs_aidx) in &matched {
            let pos = conf.atom_pos(template_aidx as u32)?;
            let mut new_at = EmbeddedAtom::at(rs_aidx as u32, Point2D::new(pos.x, pos.y));
            new_at.fixed = true;
            self.eatoms.entry(rs_aidx as u32).or_insert(new_at);
        }
        self.setup_new_neighs()?;
        self.setup_attachment_points()?;
        Ok(true)
    }

    /// `embedFusedRings`
    fn embed_fused_rings(&mut self, fused: &[Vec<i32>], use_ring_templates: bool) -> Result<()> {
        // Look for a template for the whole system. Failing that simplify the
        // system to a set of core atoms and look for a template for those. If
        // that fails, start from a single ring. Then add rings one by one.
        //
        // Look for a template that matches the entire fused ring system. For
        // single rings, only use templates for macrocycles (size > 8).
        if use_ring_templates
            && (fused.len() > 1 || (fused.len() == 1 && fused[0].len() > 8))
            && self.match_to_template(&union_rings(fused, None))?
        {
            return Ok(());
        }
        let mut coords: Vec<IntPoint2DMap> = Vec::with_capacity(fused.len());
        for ring in fused {
            let mut rc = embed_ring(ring);
            Self::mirror_trans_ring_atoms(self.mol, ring, &mut rc)?;
            coords.push(rc);
        }
        let mut done_rings: Vec<usize> = Vec::new();
        if use_ring_templates {
            let (core_rings, core_ids) = find_core_rings(fused, self.mol);
            if core_rings.len() > 1
                && core_rings.len() < fused.len()
                && self.match_to_template(&union_rings(&core_rings, None))?
            {
                done_rings = core_ids;
            }
        }
        // if not embedded, find a ring as a starting point
        if done_rings.is_empty() {
            let first = pick_first_ring_to_embed(self.mol, fused);
            self.init_from_ring_coords(&fused[first], &coords[first])?;
            done_rings.push(first);
        }
        let funion = union_rings(fused, None);
        while self.eatoms.len() < funion.len() {
            let (common, next) = find_next_ring_to_embed(&done_rings, fused)?;
            let mut emb_ring = EmbeddedFrag::empty(self.mol);
            emb_ring.init_from_ring_coords(&fused[next], &coords[next])?;
            let mut pins = Vec::new();
            if common.len() == 1 {
                let trans = self.compute_one_atom_trans(common[0] as u32, &emb_ring)?;
                emb_ring.transform(&trans);
                pins.push(common[0]);
            } else {
                let (a1, a2) = (common[0], *common.last().expect("common"));
                pins.push(a1);
                pins.push(a2);
                let trans = self.compute_two_atom_trans(a1 as u32, a2 as u32, &coords[next])?;
                emb_ring.transform(&trans);
                self.reflect_if_necessary_density(&mut emb_ring, a1 as u32, a2 as u32)?;
            }
            self.merge_ring(&emb_ring, common.len(), &pins);
            done_rings.push(next);
        }
        Ok(())
    }

    /// `computeOneAtomTrans`
    fn compute_one_atom_trans(
        &mut self,
        comm_aid: u32,
        other: &EmbeddedFrag,
    ) -> Result<Transform2D> {
        let rcr = self.ea(comm_aid).loc;
        let oeatm = other.at(comm_aid)?;
        let ccr = oeatm.loc;
        let (onb1, onb2) = (oeatm.nbr1, oeatm.nbr2);
        if onb1 < 0 || onb2 < 0 {
            return Err(RDKitError::invariant("bad ring neighbors"));
        }
        let mut mid = other.at(onb1)?.loc;
        mid += other.at(onb2)?.loc;
        mid *= 0.5;
        let (nb1, nb2) = (self.ea(comm_aid).nbr1, self.ea(comm_aid).nbr2);
        let nbp1 = self.ea(nb1).loc;
        let nbp2 = self.ea(nb2).loc;
        let ang = self.ea(comm_aid).angle;
        let largest = 2.0 * PI - ang;
        let bpt = compute_bisect_point(&rcr, largest, &nbp1, &nbp2);
        let mut t = Transform2D::identity();
        t.set_transform(&rcr, &bpt, &ccr, &mid);
        Ok(t)
    }

    /// `computeTwoAtomTrans`
    fn compute_two_atom_trans(
        &self,
        aid1: u32,
        aid2: u32,
        nring: &IntPoint2DMap,
    ) -> Result<Transform2D> {
        let loc1 = nring.get(&aid1).ok_or_else(|| missing(aid1))?;
        let loc2 = nring.get(&aid2).ok_or_else(|| missing(aid2))?;
        let ref1 = self.at(aid1)?.loc;
        let ref2 = self.at(aid2)?.loc;
        let mut t = Transform2D::identity();
        t.set_transform(&ref1, &ref2, loc1, loc2);
        Ok(t)
    }

    /// `Reflect`
    fn reflect(&mut self, loc1: &Point2D, loc2: &Point2D) {
        for e in self.eatoms.values_mut() {
            e.reflect(loc1, loc2);
        }
    }

    /// `reflectIfNecessaryCisTrans`
    fn reflect_if_necessary_cis_trans(
        &mut self,
        emb: &mut EmbeddedFrag,
        ct_case: u32,
        aid1: u32,
        aid2: u32,
    ) -> Result<()> {
        let p1 = self.ea(aid1).loc;
        let (mut r_loc, p1norm);
        if ct_case == 1 {
            p1norm = emb.ea(aid1).normal;
            let ring_atm = emb.ea(aid1).cis_trans_nbr;
            if self.has(ring_atm) {
                r_loc = self.ea(ring_atm).loc;
            } else {
                crate::rd_log_warning!(
                    "Warning: stereochemistry around double bond may be incorrect in depiction."
                );
                return Ok(());
            }
        } else {
            p1norm = self.ea(aid1).normal;
            let ring_atm = self.ea(aid1).cis_trans_nbr;
            r_loc = emb.ea(ring_atm).loc;
        }
        r_loc -= p1;
        let dot = r_loc.dot_product(&p1norm);
        let p2 = self.ea(aid2).loc;
        if dot < 0.0 {
            emb.reflect(&p1, &p2);
        }
        Ok(())
    }

    /// `reflectIfNecessaryThirdPt`
    fn reflect_if_necessary_third_pt(
        &mut self,
        emb: &mut EmbeddedFrag,
        aid1: u32,
        aid2: u32,
        aid3: u32,
    ) -> Result<()> {
        let pt1 = self.ea(aid1).loc;
        let pt2 = self.ea(aid2).loc;
        let mut normal = pt2;
        normal -= pt1;
        normal.rotate90();
        let oth3 = emb.at(aid3)?.loc - pt1;
        let pt3 = self.ea(aid3).loc - pt1;
        let dot1 = normal.dot_product(&pt3);
        let dot2 = normal.dot_product(&oth3);
        if dot1 * dot2 < 0.0 {
            emb.reflect(&pt1, &pt2);
        }
        Ok(())
    }

    /// `reflectIfNecessaryDensity`
    fn reflect_if_necessary_density(
        &mut self,
        emb: &mut EmbeddedFrag,
        aid1: u32,
        aid2: u32,
    ) -> Result<()> {
        let pin1 = self.ea(aid1).loc;
        let pin2 = self.ea(aid2).loc;
        let mut dn = 0.0;
        let mut dr = 0.0;
        for (oid, oe) in &emb.eatoms {
            if self.eatoms.contains_key(oid) {
                continue;
            }
            let loc1 = oe.loc;
            let rloc1 = reflect_point(&loc1, &pin1, &pin2);
            for te in self.eatoms.values() {
                let mut t1 = te.loc;
                t1 -= loc1;
                let td = t1.length();
                let mut rt1 = te.loc;
                rt1 -= rloc1;
                let rtd = rt1.length();
                if td > 1.0e-3 {
                    dn += 1.0 / td;
                } else {
                    dn += 1000.0;
                }
                if rtd > 1.0e-3 {
                    dr += 1.0 / rtd;
                } else {
                    dr += 1000.0;
                }
            }
        }
        if dn - dr > 1.0e-4 {
            emb.reflect(&pin1, &pin2);
        }
        Ok(())
    }

    /// `initFromRingCoords`
    fn init_from_ring_coords(&mut self, ring: &[i32], nring: &IntPoint2DMap) -> Result<()> {
        let largest = PI * (1.0 - (2.0 / ring.len() as f64));
        let mut prev = *ring.last().expect("ring");
        for (cnt, &ai) in ring.iter().enumerate() {
            let eatm = EmbeddedAtom {
                loc: *nring.get(&(ai as u32)).ok_or_else(|| missing(ai))?,
                aid: ai as u32,
                angle: largest,
                nbr1: prev,
                ..Default::default()
            };
            if cnt > 0 {
                self.ea(prev).nbr2 = ai;
            }
            self.eatoms.insert(ai as u32, eatm);
            prev = ai;
        }
        self.ea(prev).nbr2 = ring[0];
        Ok(())
    }

    /// `mergeRing`
    fn merge_ring(&mut self, emb_ring: &EmbeddedFrag, n_common: usize, pins: &[i32]) {
        for (&aid, o) in &emb_ring.eatoms {
            if !self.eatoms.contains_key(&aid) {
                self.eatoms.insert(aid, o.clone());
            } else if n_common <= 2 && pins.contains(&(aid as i32)) {
                let e = self.ea(aid);
                e.angle += o.angle;
                if e.nbr1 == o.nbr1 {
                    e.nbr1 = o.nbr2;
                } else if e.nbr1 == o.nbr2 {
                    e.nbr1 = o.nbr1;
                } else if e.nbr2 == o.nbr1 {
                    e.nbr2 = o.nbr2;
                } else if e.nbr2 == o.nbr2 {
                    e.nbr2 = o.nbr1;
                }
            }
        }
    }

    /// `addNonRingAtom`
    fn add_non_ring_atom(&mut self, aid: u32, to_aid: u32) -> Result<()> {
        if self.eatoms.contains_key(&aid) || !self.eatoms.contains_key(&to_aid) {
            return Err(RDKitError::precondition("addNonRingAtom"));
        }
        if self.ea(to_aid).angle > 0.0 {
            self.add_atom_to_atom_with_ang(aid, to_aid)?;
        } else {
            self.add_atom_to_atom_with_no_ang(aid, to_aid)?;
        }
        // neighs.erase(std::remove(...)); aid occurs once
        let neighs = &mut self.ea(to_aid).neighs;
        if let Some(p) = neighs.iter().position(|&x| x == aid as i32) {
            neighs.remove(p);
        }
        self.update_new_neighs(aid)
    }

    /// `addAtomToAtomWithAng`
    fn add_atom_to_atom_with_ang(&mut self, aid: u32, to_aid: u32) -> Result<()> {
        let r = self.ea(to_aid).clone();
        let ref_loc = r.loc;
        if r.angle <= 0.0 {
            return Err(RDKitError::precondition("addAtomToAtomWithAng"));
        }
        let nnbr = r.neighs.len();
        let rem_angle = 2.0 * PI - r.angle;
        let mut curr_angle = rem_angle / (1 + nnbr) as f64;
        self.ea(to_aid).angle += curr_angle;
        let nb1 = self.at(r.nbr1)?.loc;
        let nb2 = self.at(r.nbr2)?.loc;
        let can_choose = r.rot_dir == 0;
        if can_choose {
            self.ea(to_aid).rot_dir = rotation_dir(&ref_loc, &nb1, &nb2, rem_angle);
        }
        curr_angle *= self.ea(to_aid).rot_dir as f64;
        let mut rtrans = Transform2D::identity();
        rtrans.set_transform_rotation(&ref_loc, curr_angle);
        let mut curr_loc = nb2;
        rtrans.transform_point(&mut curr_loc);
        if can_choose && (rem_angle - PI).abs() < 1e-3 {
            let mut curr_loc2 = nb2;
            rtrans.set_transform_rotation(&ref_loc, -curr_angle);
            rtrans.transform_point(&mut curr_loc2);
            if self.find_num_neigh(&curr_loc, 0.5) > self.find_num_neigh(&curr_loc2, 0.5) {
                curr_loc = curr_loc2;
                self.ea(to_aid).rot_dir *= -1;
            }
        }
        self.ea(to_aid).nbr2 = aid as i32;
        let tpt = curr_loc - ref_loc;
        let mut norm = Point2D::new(-tpt.y, tpt.x);
        let tp1 = curr_loc + norm;
        let tp2 = curr_loc - norm;
        let nccw = self.find_num_neigh(&tp1, NEIGH_RADIUS);
        let ncw = self.find_num_neigh(&tp2, NEIGH_RADIUS);
        norm.normalize()?;
        let (normal, ccw) = if nccw < ncw {
            (norm, false)
        } else {
            (-norm, true)
        };
        let eatm = EmbeddedAtom {
            aid,
            loc: curr_loc,
            nbr1: to_aid as i32,
            angle: -1.0,
            normal,
            ccw,
            ..Default::default()
        };
        self.eatoms.insert(aid, eatm);
        Ok(())
    }

    /// `addAtomToAtomWithNoAng`
    fn add_atom_to_atom_with_no_ang(&mut self, aid: u32, to_aid: u32) -> Result<()> {
        let r = self.at(to_aid)?.clone();
        if r.angle > 0.0 {
            return Err(RDKitError::precondition("addAtomToAtomWithNoAng"));
        }
        let ref_loc = r.loc;
        let origin = Point2D::new(0.0, 0.0);
        let mut ref_ccw = r.ccw;
        let mut curr_loc = r.normal;
        if r.cis_trans_nbr >= 0 && r.cis_trans_nbr as u32 != aid {
            ref_ccw = !ref_ccw;
            curr_loc *= -1.0;
        }
        if curr_loc.length_sq() <= 1.0e-8 {
            return Err(RDKitError::invariant("zero normal"));
        }
        let atm = self.mol.atom(to_aid);
        let deg = self.mol.atom_degree(to_aid);
        let mut angle = compute_sub_angle(deg, atm.hybridization());
        let mut flip_norm = false;
        if self.ea(to_aid).nbr1 >= 0 {
            let e = self.ea(to_aid);
            e.angle = angle;
            e.nbr2 = aid as i32;
        } else {
            let mut norm = self.at(to_aid)?.normal;
            let mut rt = Transform2D::identity();
            rt.set_transform_rotation(&origin, angle);
            rt.transform_point(&mut norm);
            let e = self.ea(to_aid);
            e.normal = norm;
            e.nbr1 = aid as i32;
            flip_norm = true;
        }
        angle -= PI / 2.0;
        if !ref_ccw {
            angle *= -1.0;
        }
        let mut trans = Transform2D::identity();
        trans.set_transform_rotation(&origin, angle);
        trans.transform_point(&mut curr_loc);
        curr_loc *= bond_len();
        curr_loc += ref_loc;
        let tpt = ref_loc - curr_loc;
        let mut norm = Point2D::new(-tpt.y, tpt.x);
        if ref_ccw ^ flip_norm {
            norm *= -1.0;
        }
        norm.normalize()?;
        let eatm = EmbeddedAtom {
            loc: curr_loc,
            normal: norm,
            nbr1: to_aid as i32,
            angle: -1.0,
            ccw: (!ref_ccw) ^ flip_norm,
            ..Default::default()
        };
        self.eatoms.insert(aid, eatm);
        Ok(())
    }

    /// `findCommonAtoms`
    fn find_common_atoms(&self, other: &EmbeddedFrag) -> Vec<i32> {
        self.eatoms
            .keys()
            .filter(|k| other.eatoms.contains_key(k))
            .map(|&k| k as i32)
            .collect()
    }

    /// `mergeNoCommon`
    fn merge_no_common(&mut self, emb: &mut EmbeddedFrag, to_aid: u32, nbr_aid: u32) -> Result<()> {
        self.add_non_ring_atom(nbr_aid, to_aid)?;
        emb.add_non_ring_atom(to_aid, nbr_aid)?;
        let mut comm = vec![to_aid as i32, nbr_aid as i32];
        self.merge_with_common(emb, &mut comm)
    }

    /// `mergeWithCommon`
    fn merge_with_common(&mut self, emb: &mut EmbeddedFrag, comm: &mut Vec<i32>) -> Result<()> {
        if comm.is_empty() {
            return Err(RDKitError::precondition("mergeWithCommon"));
        }
        let mut ct_case = 0;
        if comm.len() == 1 {
            let c = comm[0];
            let other_atom;
            if self.ea(c).cis_trans_nbr >= 0 {
                ct_case = 2;
                other_atom = self.ea(c).nbr1;
                emb.add_non_ring_atom(other_atom as u32, c as u32)?;
            } else if emb.ea(c).cis_trans_nbr >= 0 {
                ct_case = 1;
                other_atom = emb.ea(c).nbr1;
                self.add_non_ring_atom(other_atom as u32, c as u32)?;
            } else {
                other_atom = self.ea(c).nbr1;
                if other_atom >= 0 {
                    emb.add_non_ring_atom(other_atom as u32, c as u32)?;
                }
            }
            if other_atom >= 0 {
                comm.push(other_atom);
            }
        }
        let rtrans = if comm.len() == 1 {
            self.compute_one_atom_trans(comm[0] as u32, emb)?
        } else {
            let (c1, c2) = (comm[0], comm[1]);
            let ref1 = self.at(c1)?.loc;
            let ref2 = self.at(c2)?.loc;
            let oth1 = emb.at(c1)?.loc;
            let oth2 = emb.at(c2)?.loc;
            let mut t = Transform2D::identity();
            t.set_transform(&ref1, &ref2, &oth1, &oth2);
            t
        };
        emb.transform(&rtrans);
        if comm.len() >= 2 {
            if ct_case > 0 {
                self.reflect_if_necessary_cis_trans(emb, ct_case, comm[0] as u32, comm[1] as u32)?;
            } else if comm.len() == 2 {
                self.reflect_if_necessary_density(emb, comm[0] as u32, comm[1] as u32)?;
            } else {
                self.reflect_if_necessary_third_pt(
                    emb,
                    comm[0] as u32,
                    comm[1] as u32,
                    comm[2] as u32,
                )?;
            }
        }
        for (&aid, o) in &emb.eatoms {
            if !comm.contains(&(aid as i32)) {
                self.eatoms.insert(aid, o.clone());
                if !o.neighs.is_empty() && !self.attach_pts.contains(&(aid as i32)) {
                    self.attach_pts.push_back(aid as i32);
                }
            } else {
                let e = self.ea(aid);
                if o.cis_trans_nbr >= 0 {
                    e.cis_trans_nbr = o.cis_trans_nbr;
                    e.normal = o.normal;
                    e.ccw = o.ccw;
                }
                if o.angle > 0.0 {
                    e.angle = o.angle;
                    e.nbr1 = o.nbr1;
                    e.nbr2 = o.nbr2;
                }
            }
        }
        for &c in comm.iter() {
            self.update_new_neighs(c as u32)?;
        }
        Ok(())
    }

    fn erase_attach_pt(&mut self, aid: i32) {
        if let Some(p) = self.attach_pts.iter().position(|&a| a == aid) {
            self.attach_pts.remove(p);
        }
    }

    /// `mergeFragsWithComm`
    fn merge_frags_with_comm(&mut self, efrags: &mut Vec<EmbeddedFrag<'a>>) -> Result<()> {
        loop {
            let mut found: Option<(usize, Vec<i32>)> = None;
            for (i, f) in efrags.iter().enumerate() {
                if !f.done {
                    let comm = self.find_common_atoms(f);
                    if !comm.is_empty() {
                        found = Some((i, comm));
                        break;
                    }
                }
            }
            let Some((i, mut comm)) = found else {
                break;
            };
            let mut other = efrags.remove(i);
            self.merge_with_common(&mut other, &mut comm)?;
            for &c in &comm {
                if self.at(c)?.neighs.is_empty() && self.attach_pts.contains(&c) {
                    self.erase_attach_pt(c);
                }
            }
        }
        Ok(())
    }

    /// `expandEfrag`
    pub fn expand_efrag(
        &mut self,
        nratms: &mut Vec<u32>,
        efrags: &mut Vec<EmbeddedFrag<'a>>,
    ) -> Result<()> {
        self.merge_frags_with_comm(efrags)?;
        while let Some(&aid) = self.attach_pts.front() {
            let nbrs = self.ea(aid).neighs.clone();
            if nbrs.is_empty() {
                return Err(RDKitError::invariant("attachment point without neighbors"));
            }
            for nbri in nbrs {
                if let Some(p) = nratms.iter().position(|&a| a as i32 == nbri) {
                    self.add_non_ring_atom(nbri as u32, aid as u32)?;
                    nratms.remove(p);
                } else if let Some(fi) = efrags
                    .iter()
                    .position(|f| !f.done && f.eatoms.contains_key(&(nbri as u32)))
                {
                    let mut other = efrags.remove(fi);
                    self.merge_no_common(&mut other, aid as u32, nbri as u32)?;
                    if self.at(nbri)?.neighs.is_empty() && self.attach_pts.contains(&nbri) {
                        self.erase_attach_pt(nbri);
                    }
                }
            }
            self.attach_pts.pop_front();
            self.ea(aid).neighs.clear();
            self.merge_frags_with_comm(efrags)?;
        }
        Ok(())
    }

    /// `Transform`
    pub fn transform(&mut self, trans: &Transform2D) {
        for e in self.eatoms.values_mut() {
            e.transform(trans);
        }
    }

    /// `Translate`
    pub fn translate(&mut self, shift: &Point2D) {
        for e in self.eatoms.values_mut() {
            e.loc += *shift;
        }
    }

    /// `computeBox`; returns `(px, nx, py, ny)`
    pub fn compute_box(&mut self) -> (f64, f64, f64, f64) {
        self.px = -1.0e8;
        self.nx = 1.0e8;
        self.py = -1.0e8;
        self.ny = 1.0e8;
        for e in self.eatoms.values() {
            self.px = cmax(self.px, e.loc.x);
            self.nx = cmin(self.nx, e.loc.x);
            self.py = cmax(self.py, e.loc.y);
            self.ny = cmin(self.ny, e.loc.y);
        }
        self.nx *= -1.0;
        self.ny *= -1.0;
        (self.px, self.nx, self.py, self.ny)
    }

    /// `canonicalizeOrientation`
    pub fn canonicalize_orientation(&mut self) -> Result<()> {
        if self.eatoms.len() <= 1 {
            return Ok(());
        }
        let mut cent = Point2D::new(0.0, 0.0);
        for e in self.eatoms.values() {
            cent += e.loc;
        }
        cent *= 1.0 / self.eatoms.len() as f64;
        let (mut xx, mut xy, mut yy) = (0.0, 0.0, 0.0);
        for e in self.eatoms.values_mut() {
            e.loc -= cent;
            xx = fmuladd(e.loc.x, e.loc.x, xx);
            xy = fmuladd(e.loc.x, e.loc.y, xy);
            yy = fmuladd(e.loc.y, e.loc.y, yy);
        }
        let d = fmuladd(xx - yy, xx - yy, 4.0 * xy * xy).sqrt();
        let mut eig1 = Point2D::new(2.0 * xy, (yy - xx) + d);
        if eig1.length() <= 1e-4 {
            return Ok(());
        }
        let e_val1 = (xx + yy + d) / 2.0;
        eig1.normalize()?;
        let mut eig2 = Point2D::new(2.0 * xy, (yy - xx) - d);
        let e_val2 = (xx + yy - d) / 2.0;
        if eig2.length() > 1e-4 {
            eig2.normalize()?;
            if e_val2 > e_val1 {
                std::mem::swap(&mut eig1, &mut eig2);
            }
        }
        let mut trans = Transform2D::identity();
        trans.set_val(0, 0, eig1.x);
        trans.set_val(1, 0, -eig1.y);
        trans.set_val(0, 1, eig1.y);
        trans.set_val(1, 1, eig1.x);
        self.transform(&trans);
        Ok(())
    }

    /// `computeDistMat`: the packed lower triangle of the inter-atom
    /// distances, indexed `i * (i - 1) / 2 + j` for `i > j`.
    fn compute_dist_mat(&self) -> Vec<f64> {
        let na = self.mol.num_atoms() as usize;
        let mut dmat = vec![0.0; na * na.saturating_sub(1) / 2];
        for (&ai, ei) in &self.eatoms {
            for (&aj, ej) in self.eatoms.range(..ai) {
                let (i, j) = (ai as usize, aj as usize);
                dmat[i * (i - 1) / 2 + j] = (ej.loc - ei.loc).length();
            }
        }
        dmat
    }

    /// `mimicDistMatAndDensityCostFunc`
    fn mimic_dist_mat_and_density_cost(&self, dmat: Option<&[f64]>, mimic_dmat_wt: f64) -> f64 {
        let na = self.mol.num_atoms() as usize;
        if na < 2 {
            return 0.0;
        }
        let dmat_2d = self.compute_dist_mat();
        let (mut res1, mut res2) = (0.0, 0.0);
        for (i, &d) in dmat_2d.iter().enumerate() {
            let d2 = d * d;
            res1 += if d2 > 1.0e-3 { 1.0 / d2 } else { 1000.0 };
            if let Some(target) = dmat.map(|m| m[i]).filter(|&t| t >= 0.0) {
                let dd = d - target;
                res2 += dd * dd;
            }
        }
        let wt = mimic_dmat_wt.clamp(0.0, 1.0);
        (1.0 - wt) * res1 + wt * res2
    }

    /// `permuteBonds`: swaps the positions of the branches at `aid1` and
    /// `aid2` by reflecting both through the line bisecting the angle they
    /// make at `aid`.
    fn permute_bonds(&mut self, aid: u32, aid1: u32, aid2: u32) -> Result<()> {
        let rl1 = self.at(aid)?.loc;
        let rl2 = (self.at(aid1)?.loc + self.at(aid2)?.loc) * 0.5;
        let mut frag_a = Vec::new();
        recurse_atom_one_side(aid1, aid, self.mol, &mut frag_a);
        let mut frag_b = Vec::new();
        recurse_atom_one_side(aid2, aid, self.mol, &mut frag_b);
        for fi in frag_a.into_iter().chain(frag_b) {
            self.ea(fi).reflect(&rl1, &rl2);
        }
        Ok(())
    }

    /// `randomSampleFlipsAndPermutations`: keeps the best of `n_samples`
    /// random sets of `n_bonds_per_sample` rotatable-bond flips (and, with
    /// `permute_deg4_nodes`, degree 4 permutations), scored by the density
    /// and, when `dmat` is given, the distance to it.
    pub fn random_sample_flips_and_permutations(
        &mut self,
        n_bonds_per_sample: u32,
        n_samples: u32,
        seed: i32,
        dmat: Option<&[f64]>,
        mimic_dmat_wt: f64,
        permute_deg4_nodes: bool,
    ) -> Result<()> {
        let mol = self.mol;
        let rot_bonds: Vec<u32> = get_all_rotatable_bonds(mol)
            .into_iter()
            .filter(|&bid| {
                let b = mol.bond(bid);
                self.has(b.begin_atom_idx()) && self.has(b.end_atom_idx())
            })
            .collect();
        let nb = rot_bonds.len();
        let mut deg4: Vec<(u32, Vec<u32>, Vec<u32>)> = Vec::new();
        if permute_deg4_nodes {
            for atom in 0..mol.num_atoms() {
                if mol.atom_degree(atom) == 4 && mol.ring_info().num_atom_rings(atom)? == 0 {
                    let (aids, bids) = get_nbr_atom_and_bond_ids(atom, mol)?;
                    let all_in = aids
                        .iter()
                        .all(|&a| self.eatoms.get(&a).is_some_and(|e| !e.fixed));
                    if all_in {
                        deg4.push((atom, aids, bids));
                    }
                }
            }
        }
        let nt = nb + deg4.len();
        if nt == 0 {
            return Ok(());
        }
        let per_sample = nt.min(n_bonds_per_sample as usize);
        let snapshot = |f: &Self| -> BTreeMap<u32, Point2D> {
            f.eatoms.iter().map(|(&k, e)| (k, e.loc)).collect()
        };
        let mut best_crd = snapshot(self);
        let mut best_dens = self.mimic_dist_mat_and_density_cost(dmat, mimic_dmat_wt);
        with_global_generator(seed, |rng| -> Result<()> {
            for _ in 0..n_samples {
                for _ in 0..per_sample {
                    let ri = rng.uniform_int(0, nt as i32 - 1) as usize;
                    if ri < nb {
                        self.flip_about_bond(rot_bonds[ri], true)?;
                        continue;
                    }
                    let (ai, aids, bids) = &deg4[ri - nb];
                    let locs: Vec<Point2D> = aids
                        .iter()
                        .map(|&a| self.at(a).map(|e| e.loc))
                        .collect::<Result<_>>()?;
                    let pairs = find_bonds_pairs_to_permute_deg4(&self.at(*ai)?.loc, bids, &locs)?;
                    let fbi = usize::from(rng.next_f64() > 0.5);
                    let aid1 = mol.bond(pairs[fbi].0).other_atom_idx(*ai)?;
                    let aid2 = mol.bond(pairs[fbi].1).other_atom_idx(*ai)?;
                    self.permute_bonds(*ai, aid1, aid2)?;
                }
                let density = self.mimic_dist_mat_and_density_cost(dmat, mimic_dmat_wt);
                if best_dens - density > 1.0e-4 {
                    best_dens = density;
                    best_crd = snapshot(self);
                }
            }
            Ok(())
        })?;
        for (k, e) in self.eatoms.iter_mut() {
            e.loc = best_crd[k];
        }
        Ok(())
    }

    /// `flipAboutBond`
    fn flip_about_bond(&mut self, bond_id: u32, flip_end: bool) -> Result<()> {
        let mol = self.mol;
        if mol.ring_info().num_bond_rings(bond_id)? > 0 {
            return Err(RDKitError::invariant("cannot flip about a ring bond"));
        }
        let b = mol.bond(bond_id);
        let (mut beg, mut end) = (b.begin_atom_idx(), b.end_atom_idx());
        if !flip_end {
            std::mem::swap(&mut beg, &mut end);
        }
        self.at(beg)?;
        self.at(end)?;
        let mut end_side = Vec::new();
        recurse_atom_one_side(end, beg, mol, &mut end_side);
        let n_fixed = self.eatoms.values().filter(|e| e.fixed).count();
        let mut n_end_fixed = 0;
        if n_fixed > 0 {
            for &a in &end_side {
                if self.ea(a).fixed {
                    n_end_fixed += 1;
                }
            }
        }
        if n_end_fixed > 0 {
            return Ok(());
        }
        let nats = self.eatoms.len();
        // size_t arithmetic upstream: the difference wraps around
        let end_side_flip = nats.wrapping_sub(end_side.len()) >= end_side.len();
        let keys: Vec<u32> = self.eatoms.keys().copied().collect();
        for k in keys {
            let in_end = end_side.contains(&(k as i32));
            if end_side_flip ^ !in_end {
                // upstream passes references to the bond atoms' locations:
                // reflecting one of those atoms changes the mirror line seen
                // by its own normal and by every later atom
                let bl = self.eatoms[&beg].loc;
                let el = self.eatoms[&end].loc;
                let e = self.ea(k);
                let temp = e.loc + e.normal;
                e.loc = reflect_point(&e.loc, &bl, &el);
                let l1 = if k == beg { e.loc } else { bl };
                let l2 = if k == end { e.loc } else { el };
                let temp = reflect_point(&temp, &l1, &l2);
                e.normal = temp - e.loc;
                e.ccw = !e.ccw;
            }
        }
        Ok(())
    }

    /// `flipAboutSpiroCenter`
    fn flip_about_spiro_center(&mut self, spiro: u32) -> Result<()> {
        let mol = self.mol;
        let rings: Vec<&Vec<i32>> = mol
            .ring_info()
            .atom_rings()
            .iter()
            .filter(|r| r.contains(&(spiro as i32)))
            .collect();
        if rings.len() != 2 {
            return Err(RDKitError::postcondition("Spiro must have exactly 2 rings"));
        }
        let target: HashSet<i32> = rings[0].iter().copied().collect();
        let ring_nbrs: Vec<u32> = mol
            .atom_neighbors(spiro)
            .filter(|n| target.contains(&(*n as i32)))
            .collect();
        if ring_nbrs.len() != 2 {
            return Err(RDKitError::invariant(
                "Spiro atom should have exactly 2 neighbors in each ring",
            ));
        }
        let mut to_flip = Vec::new();
        recurse_atom_one_side(ring_nbrs[0], spiro, mol, &mut to_flip);
        let spiro_loc = self.at(spiro)?.loc;
        let mid = (self.at(ring_nbrs[0])?.loc + self.at(ring_nbrs[1])?.loc) * 0.5;
        for &a in &to_flip {
            if self.at(a)?.fixed {
                return Ok(());
            }
        }
        for &a in &to_flip {
            self.ea(a).reflect(&spiro_loc, &mid);
        }
        Ok(())
    }

    /// `findCollisions`
    pub fn find_collisions(
        &mut self,
        dmat: &[f64],
        include_bonds: bool,
        include_bonded: bool,
    ) -> Vec<(u32, u32)> {
        let mol = self.mol;
        let mut res = Vec::new();
        for e in self.eatoms.values_mut() {
            e.density = 0.0;
        }
        let col_thres2 = COLLISION_THRES * COLLISION_THRES;
        let keys: Vec<u32> = self.eatoms.keys().copied().collect();
        let locs: Vec<Point2D> = self.eatoms.values().map(|e| e.loc).collect();
        let mut dens = vec![0.0f64; keys.len()];
        for i in 0..keys.len() {
            let f1 = if mol.atom(keys[i]).atomic_num() != 6 {
                HETEROATOM_COLL_SCALE
            } else {
                1.0
            };
            for j in 0..i {
                let f2 = if mol.atom(keys[j]).atomic_num() != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                let mut ptj = locs[j];
                ptj -= locs[i];
                let mut d2 = ptj.length_sq();
                if d2 > 1.0e-3 {
                    dens[i] += 1.0 / d2;
                    dens[j] += 1.0 / d2;
                } else {
                    dens[i] += 1000.0;
                    dens[j] += 1000.0;
                }
                d2 /= f1 * f2;
                if d2 < col_thres2
                    && (include_bonded || mol.bond_idx_between(keys[i], keys[j]).is_none())
                {
                    res.push((keys[i], keys[j]));
                }
            }
        }
        for (e, d) in self.eatoms.values_mut().zip(dens) {
            e.density = d;
        }
        if include_bonds {
            let na = mol.num_atoms() as usize;
            for b1 in mol.bonds() {
                let (beg1, end1) = (b1.begin_atom_idx(), b1.end_atom_idx());
                if !(self.eatoms.contains_key(&beg1) && self.eatoms.contains_key(&end1)) {
                    continue;
                }
                for b2 in mol.bonds() {
                    if b2.idx() <= b1.idx() {
                        continue;
                    }
                    let (beg2, end2) = (b2.begin_atom_idx(), b2.end_atom_idx());
                    if beg1 == beg2 || beg1 == end2 || end1 == beg2 || end1 == end2 {
                        continue;
                    }
                    if self.eatoms.contains_key(&beg2)
                        && self.eatoms.contains_key(&end2)
                        && segments_intersect(
                            &self.eatoms[&beg1].loc,
                            &self.eatoms[&end1].loc,
                            &self.eatoms[&beg2].loc,
                            &self.eatoms[&end2].loc,
                        )
                    {
                        res.push(find_farthest_pair(beg1, end1, beg2, end2, na, dmat));
                    }
                }
            }
        }
        res
    }

    /// `totalDensity`
    fn total_density(&self) -> f64 {
        self.eatoms.values().fold(0.0, |acc, e| e.density + acc)
    }

    /// `openAngles`
    fn open_angles(&mut self, dmat: &[f64], aid1: u32, aid2: u32) -> Result<()> {
        let mol = self.mol;
        let deg1 = mol.atom_degree(aid1);
        let deg2 = mol.atom_degree(aid2);
        let fixed1 = self.at(aid1)?.fixed;
        let fixed2 = self.at(aid2)?.fixed;
        if (deg1 > 1 || fixed1) && (deg2 > 1 || fixed2) {
            return Ok(());
        }
        let (aid_a, aid_b, kind);
        if (deg1 == 1 && !fixed1) && (deg2 == 1 && !fixed2) {
            aid_a = find_deg1_neighbor(mol, aid1)?;
            aid_b = find_deg1_neighbor(mol, aid2)?;
            kind = 1;
        } else if (deg1 == 1 && !fixed1) && (deg2 > 1 || fixed2) {
            aid_a = find_deg1_neighbor(mol, aid1)?;
            aid_b = find_closest_neighbor(mol, dmat, aid_a, aid2);
            kind = 2;
        } else {
            aid_b = find_deg1_neighbor(mol, aid2)?;
            aid_a = find_closest_neighbor(mol, dmat, aid_b, aid1);
            kind = 3;
        }
        let v2 = self.at(aid1)?.loc - self.at(aid_a)?.loc;
        let v1 = self.at(aid_b)?.loc - self.at(aid_a)?.loc;
        let cross = fmuladd(v1.x, v2.y, -(v1.y * v2.x));
        let mut t1 = Transform2D::identity();
        let mut t2 = Transform2D::identity();
        match kind {
            1 => {
                let mut angle = ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let la = self.ea(aid_a).loc;
                t1.set_transform_rotation(&la, angle);
                let lb = self.ea(aid_b).loc;
                t2.set_transform_rotation(&lb, -angle);
                t1.transform_point(&mut self.ea(aid1).loc);
                t2.transform_point(&mut self.ea(aid2).loc);
            }
            2 => {
                let mut angle = 2.0 * ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let la = self.ea(aid_a).loc;
                t1.set_transform_rotation(&la, angle);
                t1.transform_point(&mut self.ea(aid1).loc);
            }
            _ => {
                let mut angle = -2.0 * ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let lb = self.ea(aid_b).loc;
                t2.set_transform_rotation(&lb, angle);
                t2.transform_point(&mut self.ea(aid2).loc);
            }
        }
        Ok(())
    }

    /// `tryResolvingCollisionWithBondFlip`
    fn try_bond_flip(
        &mut self,
        c: (u32, u32),
        ncols: usize,
        prev_density: f64,
        done_bonds: &mut BTreeMap<u32, u32>,
        dmat: &[f64],
    ) -> Result<bool> {
        let rot_bonds = get_rotatable_bonds(self.mol, c.0, c.1)?;
        for ri in rot_bonds {
            let cur = done_bonds.get(&ri).copied();
            if cur.is_none_or(|v| v < NUM_BONDS_FLIPS) {
                *done_bonds.entry(ri).or_insert(0) += 1;
                self.flip_about_bond(ri, true)?;
                let colls = self.find_collisions(dmat, true, true);
                let nd = self.total_density();
                if colls.len() < ncols {
                    done_bonds.insert(ri, NUM_BONDS_FLIPS);
                    return Ok(true);
                } else if colls.len() == ncols && nd < prev_density {
                    return Ok(true);
                }
                self.flip_about_bond(ri, true)?;
                self.flip_about_bond(ri, false)?;
                let colls = self.find_collisions(dmat, true, true);
                let nd = self.total_density();
                if colls.len() < ncols {
                    done_bonds.insert(ri, NUM_BONDS_FLIPS);
                    return Ok(true);
                } else if colls.len() == ncols && nd < prev_density {
                    return Ok(true);
                }
                self.flip_about_bond(ri, false)?;
            }
        }
        Ok(false)
    }

    /// `tryResolvingCollisionWithSpiroFlip`
    fn try_spiro_flip(
        &mut self,
        c: (u32, u32),
        ncols: usize,
        prev_density: f64,
        done_spiros: &mut BTreeMap<u32, u32>,
        spiro_centers: &[bool],
        dmat: &[f64],
    ) -> Result<bool> {
        let path = get_shortest_path(self.mol, c.0, c.1)?;
        let spiros: Vec<u32> = path
            .into_iter()
            .filter(|&a| spiro_centers[a as usize])
            .collect();
        for s in spiros {
            let done = done_spiros.get(&s).copied();
            if done.is_some_and(|v| v >= NUM_BONDS_FLIPS) {
                continue;
            }
            self.flip_about_spiro_center(s)?;
            let colls = self.find_collisions(dmat, true, true);
            let nd = self.total_density();
            if colls.len() < ncols {
                done_spiros.insert(s, NUM_BONDS_FLIPS);
                return Ok(true);
            } else if colls.len() == ncols && nd < prev_density {
                *done_spiros.entry(s).or_insert(0) += 1;
                return Ok(true);
            }
            self.flip_about_spiro_center(s)?;
        }
        Ok(false)
    }

    /// `removeCollisionsBondAndSpiroFlip`
    pub fn remove_collisions_bond_and_spiro_flip(&mut self) -> Result<()> {
        let mol = self.mol;
        let mut spiro = vec![false; mol.num_atoms() as usize];
        for a in 0..mol.num_atoms() {
            spiro[a as usize] = is_spiro_center(a, mol)?;
        }
        let dmat = get_distance_mat(mol);
        let mut colls = self.find_collisions(&dmat, true, true);
        let mut done_bonds = BTreeMap::new();
        let mut done_spiros = BTreeMap::new();
        let mut iter = 0;
        while iter < MAX_COLL_ITERS && !colls.is_empty() {
            let ncols = colls.len();
            let c = colls[0];
            let prev = self.total_density();
            let resolved = self.try_bond_flip(c, ncols, prev, &mut done_bonds, &dmat)?;
            if !resolved {
                self.try_spiro_flip(c, ncols, prev, &mut done_spiros, &spiro, &dmat)?;
            }
            colls = self.find_collisions(&dmat, true, true);
            iter += 1;
        }
        Ok(())
    }

    /// `removeCollisionsOpenAngles`
    pub fn remove_collisions_open_angles(&mut self) -> Result<()> {
        let dmat = get_distance_mat(self.mol);
        for (a, b) in self.find_collisions(&dmat, false, true) {
            self.open_angles(&dmat, a, b)?;
        }
        Ok(())
    }

    /// `removeCollisionsShortenBonds`
    pub fn remove_collisions_shorten_bonds(&mut self) -> Result<()> {
        let mol = self.mol;
        let dmat = get_distance_mat(mol);
        let mut colls = self.find_collisions(&dmat, false, true);
        let mut ncols = colls.len();
        let mut iter = 0;
        while ncols > 0 && iter < MAX_COLL_ITERS {
            let (mut aid1, mut aid2) = colls[0];
            let mut fixed1 = self.at(aid1)?.fixed;
            let mut fixed2 = self.at(aid2)?.fixed;
            if fixed1 && fixed2 {
                colls.remove(0);
                ncols = colls.len();
                iter += 1;
                continue;
            }
            let mut deg1 = mol.atom_degree(aid1);
            let mut deg2 = mol.atom_degree(aid2);
            if fixed1 || (deg2 > deg1 && !fixed2) {
                std::mem::swap(&mut deg1, &mut deg2);
                std::mem::swap(&mut aid1, &mut aid2);
                std::mem::swap(&mut fixed1, &mut fixed2);
            }
            let mut path = get_shortest_path(mol, aid1, aid2)?;
            if path.is_empty() {
                colls.remove(0);
            } else {
                if path[0] != aid1 {
                    return Err(RDKitError::invariant("bad path head"));
                }
                path.remove(0);
                let n_open = any_non_ring_bonds(aid1, &path, mol)?;
                if n_open > 0 {
                    for (aid, deg, fixed) in [(aid1, deg1, fixed1), (aid2, deg2, fixed2)] {
                        if deg == 1 && (aid == aid1 || !fixed) {
                            let mut loc = self.at(aid)?.loc;
                            let aid_a = find_deg1_neighbor(mol, aid)?;
                            loc -= self.ea(aid_a).loc;
                            loc *= 0.9;
                            if loc.length() > 0.75 {
                                loc += self.ea(aid_a).loc;
                                self.ea(aid).loc = loc;
                            }
                        }
                    }
                } else {
                    let mut r_path = Vec::new();
                    let mut nbr_map = BTreeMap::new();
                    recurse_deg_two_ring_atoms(aid1, mol, &mut r_path, &mut nbr_map)?;
                    if r_path.is_empty() {
                        recurse_deg_two_ring_atoms(aid2, mol, &mut r_path, &mut nbr_map)?;
                    }
                    let mut move_map: BTreeMap<i32, Point2D> = BTreeMap::new();
                    for &rpi in &r_path {
                        if self.at(rpi)?.fixed {
                            continue;
                        }
                        let nb = nbr_map[&rpi].clone();
                        let mut mv = self.ea(nb[0]).loc;
                        mv += self.ea(nb[1]).loc;
                        mv *= 0.5;
                        mv -= self.at(rpi)?.loc;
                        mv.normalize()?;
                        mv *= COLLISION_THRES;
                        move_map.insert(rpi, mv);
                    }
                    for &rpi in &r_path {
                        let mv = *move_map.entry(rpi).or_default();
                        self.ea(rpi).loc += mv;
                    }
                }
                colls = self.find_collisions(&dmat, false, true);
            }
            ncols = colls.len();
            iter += 1;
        }
        Ok(())
    }

    /// `getAtomsOnSide`
    fn get_atoms_on_side(&self, center: u32, side_start: u32, exclude: u32) -> Vec<u32> {
        let mut result = Vec::new();
        let mut queue = VecDeque::from([side_start]);
        let mut visited: HashSet<u32> = HashSet::from([center, exclude]);
        while let Some(cur) = queue.pop_front() {
            if !visited.insert(cur) {
                continue;
            }
            if self.eatoms.contains_key(&cur) {
                result.push(cur);
            }
            for b in self.mol.atom_bonds(cur) {
                let nb = if b.begin_atom_idx() == cur {
                    b.end_atom_idx()
                } else {
                    b.begin_atom_idx()
                };
                if !visited.contains(&nb) {
                    queue.push_back(nb);
                }
            }
        }
        result
    }

    /// `openAngleByIncrement`
    fn open_angle_by_increment(
        &mut self,
        prev: u32,
        center: u32,
        next: u32,
        inc: f64,
    ) -> Result<bool> {
        let cl = self.ea(center).loc;
        let v1 = self.ea(prev).loc - cl;
        let v2 = self.ea(next).loc - cl;
        let side1 = self.get_atoms_on_side(center, prev, next);
        let side2 = self.get_atoms_on_side(center, next, prev);
        let s1: HashSet<u32> = side1.iter().copied().collect();
        if side2.iter().any(|a| s1.contains(a)) {
            return Ok(false);
        }
        let contains_fixed = |this: &Self, atoms: &[u32]| -> Result<bool> {
            for &a in atoms {
                if this.at(a)?.fixed {
                    return Ok(true);
                }
            }
            Ok(false)
        };
        let f1 = contains_fixed(self, &side1)?;
        let f2 = contains_fixed(self, &side2)?;
        if f1 && f2 {
            return Ok(false);
        }
        let rotate_side1 = f2 || (!f1 && side1.len() < side2.len());
        let to_move = if rotate_side1 { side1 } else { side2 };
        let vlen = |v: &Point2D| fmuladd(v.x, v.x, v.y * v.y).sqrt();
        let angle_between = |a: &Point2D, b: &Point2D| -> f64 {
            let den = vlen(a) * vlen(b);
            if den < 1e-8 {
                return -1.0;
            }
            (fmuladd(a.x, b.x, a.y * b.y) / den).clamp(-1.0, 1.0).acos()
        };
        let original = angle_between(&v1, &v2);
        if original < 0.0 {
            return Ok(false);
        }
        let mut rotation = 0.0;
        for cand in [inc, -inc] {
            let mut trial = Transform2D::identity();
            let cl = self.ea(center).loc;
            trial.set_transform_rotation(&cl, cand);
            let mut moved = self.ea(if rotate_side1 { prev } else { next }).loc;
            trial.transform_point(&mut moved);
            let mv = moved - self.ea(center).loc;
            let fixed = if rotate_side1 { v2 } else { v1 };
            if angle_between(&mv, &fixed) > original + 1e-8 {
                rotation = cand;
                break;
            }
        }
        if rotation == 0.0 {
            return Ok(false);
        }
        let mut trans = Transform2D::identity();
        let cl = self.ea(center).loc;
        trans.set_transform_rotation(&cl, rotation);
        for a in to_move {
            trans.transform_point(&mut self.ea(a).loc);
        }
        Ok(true)
    }

    fn path_angle_collisions(&mut self, dmat: &[f64]) -> Vec<(u32, u32)> {
        let c = self.find_collisions(dmat, false, false);
        if c.is_empty() {
            self.find_collisions(dmat, true, false)
        } else {
            c
        }
    }

    /// `removeCollisionsPathAngleExpansion`
    pub fn remove_collisions_path_angle_expansion(&mut self) -> Result<()> {
        let mol = self.mol;
        let mut atom_rings: HashMap<i32, HashSet<usize>> = HashMap::new();
        for (ri, ring) in mol.ring_info().atom_rings().iter().enumerate() {
            for &a in ring {
                atom_rings.entry(a).or_default().insert(ri);
            }
        }
        let canon = |c: (u32, u32)| (c.0.min(c.1), c.0.max(c.1));
        let mut skipped: HashSet<(u32, u32)> = HashSet::new();
        let filter = |skipped: &HashSet<(u32, u32)>, cs: &[(u32, u32)]| -> Vec<(u32, u32)> {
            cs.iter()
                .copied()
                .filter(|c| !skipped.contains(&canon(*c)))
                .collect()
        };
        let dmat = get_distance_mat(mol);
        let mut all = self.path_angle_collisions(&dmat);
        let mut colls = filter(&skipped, &all);
        let mut angle_totals: HashMap<(u32, u32), HashMap<u32, f64>> = HashMap::new();
        let mut iter = 0;
        while iter < MAX_ANGLE_EXPANSION_ITERS && !colls.is_empty() {
            let collision = colls[0];
            let path = get_shortest_path(mol, collision.0, collision.1)?;
            if path.len() < 3 {
                skipped.insert(canon(collision));
                colls.remove(0);
                iter += 1;
                continue;
            }
            let mut expandable = Vec::new();
            for i in 1..path.len() - 1 {
                let (p, c, n) = (path[i - 1] as i32, path[i] as i32, path[i + 1] as i32);
                let same_ring = match (atom_rings.get(&p), atom_rings.get(&c), atom_rings.get(&n)) {
                    (Some(pr), Some(cr), Some(nr)) => {
                        pr.iter().any(|r| cr.contains(r) && nr.contains(r))
                    }
                    _ => false,
                };
                if !same_ring {
                    expandable.push(i);
                }
            }
            if expandable.is_empty() {
                skipped.insert(canon(collision));
                colls.remove(0);
                iter += 1;
                continue;
            }
            let prev_count = all.len();
            let saved: HashMap<u32, Point2D> =
                self.eatoms.iter().map(|(&k, e)| (k, e.loc)).collect();
            let mut any = false;
            for &pos in &expandable {
                let center = path[pos];
                let tot = *angle_totals
                    .entry(canon(collision))
                    .or_default()
                    .entry(center)
                    .or_insert(0.0);
                if tot >= MAX_ANGLE_EXPANSION {
                    continue;
                }
                if self.open_angle_by_increment(
                    path[pos - 1],
                    center,
                    path[pos + 1],
                    ANGLE_EXPANSION_INCREMENT,
                )? {
                    *angle_totals
                        .entry(canon(collision))
                        .or_default()
                        .entry(center)
                        .or_insert(0.0) += ANGLE_EXPANSION_INCREMENT;
                    any = true;
                }
            }
            if !any {
                skipped.insert(canon(collision));
                colls.remove(0);
                iter += 1;
                continue;
            }
            let dmat = get_distance_mat(mol);
            let new_colls = self.path_angle_collisions(&dmat);
            if new_colls.len() > prev_count {
                for (k, e) in self.eatoms.iter_mut() {
                    e.loc = saved.get(k).copied().unwrap_or_default();
                }
                angle_totals.remove(&canon(collision));
                skipped.insert(canon(collision));
                all = self.path_angle_collisions(&dmat);
                colls = filter(&skipped, &all);
            } else {
                all = new_colls;
                colls = filter(&skipped, &all);
            }
            iter += 1;
        }
        Ok(())
    }
}
