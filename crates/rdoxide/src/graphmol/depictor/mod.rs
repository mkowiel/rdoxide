//! 2D coordinate generation (`RDDepict::compute2DCoords`).
//!
//! Ported from RDKit `Code/GraphMol/Depictor/RDDepictor.cpp` (2026.09.1).
//! The RDKit depictor is used (not CoordGen); random flip sampling is not
//! ported yet.
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

pub mod depict_utils;
pub mod embedded_frag;
pub mod templates;

use crate::geometry::{Point2D, Point3D};
use crate::graphmol::aromaticity::{make_ring_neighbor_map, pick_fused_rings};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::chirality::nontetrahedral::{
    get_chiral_across_atom, get_ideal_angle_between_ligands, get_trigonal_bipyramidal_axial_bond,
    has_non_tetrahedral_stereo,
};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::{RDKitError, Result};
use depict_utils::{IntPoint2DMap, bond_len, get_atom_depict_rank};
use embedded_frag::EmbeddedFrag;
use std::sync::OnceLock;

/// `RDDepict::DepictException`, which Python sees as a `ValueError`.
pub fn depict_error(msg: impl std::fmt::Display) -> RDKitError {
    RDKitError::Value(format!("Depict error: {msg}"))
}

/// `RDDepict::Compute2DCoordParameters`
#[derive(Debug, Clone)]
pub struct Compute2DCoordParameters<'m> {
    pub coord_map: Option<&'m IntPoint2DMap>,
    pub canon_orient: bool,
    pub clear_confs: bool,
    pub n_flips_per_sample: u32,
    pub n_samples: u32,
    pub sample_seed: i32,
    pub permute_deg4_nodes: bool,
    pub force_rdkit: bool,
    pub use_ring_templates: bool,
    pub use_path_angle_expansion: bool,
}

impl Default for Compute2DCoordParameters<'_> {
    fn default() -> Self {
        Compute2DCoordParameters {
            coord_map: None,
            canon_orient: false,
            clear_confs: true,
            n_flips_per_sample: 0,
            n_samples: 0,
            sample_seed: 0,
            permute_deg4_nodes: false,
            force_rdkit: false,
            use_ring_templates: false,
            use_path_angle_expansion: true,
        }
    }
}

// upstream's rounded constants (not FRAC_1_SQRT_2 / sqrt(3)/2)
#[allow(clippy::approx_constant)]
const ISQRT2: f64 = 0.707107;
const SQRT3_2: f64 = 0.866025;

/// `getRankedAtomNeighbors` (an unstable `std::sort` by depiction rank)
fn ranked_neighbors(mol: &ROMol, atom: u32, ranks: &[i32]) -> Vec<u32> {
    let mut nbrs: Vec<u32> = mol.atom_neighbors(atom).collect();
    libcxx_sort_by(&mut nbrs, |a, b| ranks[*a as usize] < ranks[*b as usize]);
    nbrs
}

fn push_coord_frag<'a>(
    mol: &'a ROMol,
    efrags: &mut Vec<EmbeddedFrag<'a>>,
    map: &IntPoint2DMap,
) -> Result<()> {
    let mut f = EmbeddedFrag::from_coord_map(mol, map)?;
    f.seq = efrags.len();
    efrags.push(f);
    Ok(())
}

/// `embedSquarePlanar`
fn embed_square_planar<'a>(
    mol: &'a ROMol,
    atom: u32,
    efrags: &mut Vec<EmbeddedFrag<'a>>,
    ranks: &[i32],
) -> Result<()> {
    // upstream: function-local statics, fixed at the first call
    static IDEAL: OnceLock<[Point2D; 4]> = OnceLock::new();
    let ideal = IDEAL.get_or_init(|| {
        let b = bond_len();
        [
            Point2D::new(ISQRT2 * b, ISQRT2 * b),
            Point2D::new(ISQRT2 * b, -ISQRT2 * b),
            Point2D::new(-ISQRT2 * b, -ISQRT2 * b),
            Point2D::new(-ISQRT2 * b, ISQRT2 * b),
        ]
    });
    let nbrs = ranked_neighbors(mol, atom, ranks);
    let mut map = IntPoint2DMap::new();
    map.insert(atom, Point2D::new(0., 0.));
    map.insert(nbrs[0], ideal[0]);
    let mut q2_full = false;
    for &nbr in &nbrs[1..] {
        let angle = get_ideal_angle_between_ligands(mol, atom, nbrs[0], nbr);
        if (angle - 180.0).abs() < 0.1 {
            map.insert(nbr, ideal[2]);
        } else if !q2_full {
            map.insert(nbr, ideal[1]);
            q2_full = true;
        } else {
            map.insert(nbr, ideal[3]);
        }
    }
    push_coord_frag(mol, efrags, &map)
}

fn axial_atom(mol: &ROMol, cen: u32, axial: i32) -> Option<u32> {
    get_trigonal_bipyramidal_axial_bond(mol, cen, axial)
        .and_then(|b| mol.bond(b).other_atom_idx(cen).ok())
}

/// `embedTBP`
fn embed_tbp<'a>(
    mol: &'a ROMol,
    atom: u32,
    efrags: &mut Vec<EmbeddedFrag<'a>>,
    ranks: &[i32],
) -> Result<()> {
    static IDEAL: OnceLock<[Point2D; 5]> = OnceLock::new();
    let ideal = IDEAL.get_or_init(|| {
        let b = bond_len();
        [
            Point2D::new(0.0, b),
            Point2D::new(0.0, -b),
            Point2D::new(-SQRT3_2 * b, b / 2.0),
            Point2D::new(-SQRT3_2 * b, -b / 2.0),
            Point2D::new(b, 0.0),
        ]
    });
    let nbrs = ranked_neighbors(mol, atom, ranks);
    let mut map = IntPoint2DMap::new();
    map.insert(atom, Point2D::new(0., 0.));
    let ax1 = axial_atom(mol, atom, 1);
    let ax2 = axial_atom(mol, atom, -1);
    if let Some(a) = ax1 {
        map.insert(a, ideal[0]);
    }
    if let Some(a) = ax2 {
        map.insert(a, ideal[1]);
    }
    let mut which = 2;
    for &nbr in &nbrs {
        if Some(nbr) != ax1 && Some(nbr) != ax2 {
            let p = *ideal
                .get(which)
                .ok_or_else(|| RDKitError::Index("too many equatorial ligands".into()))?;
            map.insert(nbr, p);
            which += 1;
        }
    }
    push_coord_frag(mol, efrags, &map)
}

/// `embedOctahedral`
fn embed_octahedral<'a>(
    mol: &'a ROMol,
    atom: u32,
    efrags: &mut Vec<EmbeddedFrag<'a>>,
    ranks: &[i32],
) -> Result<()> {
    static IDEAL: OnceLock<[Point2D; 6]> = OnceLock::new();
    let ideal = IDEAL.get_or_init(|| {
        let b = bond_len();
        [
            Point2D::new(0.0, b),
            Point2D::new(0.0, -b),
            Point2D::new(SQRT3_2 * b, b / 2.0),
            Point2D::new(SQRT3_2 * b, -b / 2.0),
            Point2D::new(-SQRT3_2 * b, -b / 2.0),
            Point2D::new(-SQRT3_2 * b, b / 2.0),
        ]
    });
    let nbrs = ranked_neighbors(mol, atom, ranks);
    let mut map = IntPoint2DMap::new();
    map.insert(atom, Point2D::new(0., 0.));
    let (mut ax1, mut ax2): (Option<u32>, Option<u32>) = (None, None);
    for i in 0..nbrs.len() {
        let mut all90 = true;
        for j in i + 1..nbrs.len() {
            let ang = get_ideal_angle_between_ligands(mol, atom, nbrs[i], nbrs[j]);
            if (ang - 180.0).abs() < 0.1 {
                ax1 = Some(nbrs[i]);
                ax2 = Some(nbrs[j]);
                all90 = false;
                break;
            } else if (ang - 90.0).abs() > 0.1 {
                all90 = false;
            }
        }
        if all90 {
            ax1 = Some(nbrs[i]);
        }
        if ax1.is_some() {
            break;
        }
    }
    if let Some(a) = ax1 {
        map.insert(a, ideal[0]);
    }
    if let Some(a) = ax2 {
        map.insert(a, ideal[1]);
    }
    let (mut ref1, mut ref2): (Option<u32>, Option<u32>) = (None, None);
    for &nbr in &nbrs {
        if Some(nbr) == ax1 || Some(nbr) == ax2 {
            continue;
        }
        if ref1.is_none() {
            ref1 = Some(nbr);
            map.insert(nbr, ideal[2]);
            ref2 = get_chiral_across_atom(mol, atom, nbr);
            if let Some(r) = ref2 {
                map.insert(r, ideal[4]);
            }
        } else {
            if Some(nbr) == ref2 || Some(nbr) == ref1 {
                continue;
            }
            map.insert(nbr, ideal[3]);
            if let Some(a2) = get_chiral_across_atom(mol, atom, nbr) {
                map.insert(a2, ideal[5]);
            }
            break;
        }
    }
    push_coord_frag(mol, efrags, &map)
}

/// `embedNontetrahedralStereo`
fn embed_nontetrahedral<'a>(
    mol: &'a ROMol,
    efrags: &mut Vec<EmbeddedFrag<'a>>,
    ranks: &[i32],
) -> Result<()> {
    for a in 0..mol.num_atoms() {
        if !has_non_tetrahedral_stereo(mol, a) {
            continue;
        }
        match mol.atom(a).chiral_tag() {
            ChiralType::SquarePlanar => embed_square_planar(mol, a, efrags, ranks)?,
            ChiralType::TrigonalBipyramidal => embed_tbp(mol, a, efrags, ranks)?,
            ChiralType::Octahedral => embed_octahedral(mol, a, efrags, ranks)?,
            _ => {}
        }
    }
    Ok(())
}

/// `embedFusedSystems`
fn embed_fused_systems<'a>(
    mol: &'a ROMol,
    arings: &[Vec<i32>],
    efrags: &mut Vec<EmbeddedFrag<'a>>,
    coord_map: Option<&IntPoint2DMap>,
    use_ring_templates: bool,
) -> Result<()> {
    let neigh_map = make_ring_neighbor_map(arings, 0, 0);
    let cnrs = arings.len();
    let mut fus_done = vec![false; cnrs];
    let mut curr = 0usize;
    while curr < cnrs {
        let mut fused = Vec::new();
        pick_fused_rings(curr as i32, &neigh_map, &mut fused, &mut fus_done);
        let frings: Vec<Vec<i32>> = fused.iter().map(|&r| arings[r as usize].clone()).collect();
        let mut allow = use_ring_templates;
        if use_ring_templates && let Some(cm) = coord_map {
            let mut atoms = std::collections::BTreeSet::new();
            for ring in &frings {
                for &a in ring {
                    if cm.contains_key(&(a as u32)) {
                        atoms.insert(a);
                    }
                }
            }
            allow = atoms.len() < 2;
        }
        let mut f = EmbeddedFrag::fused_rings(mol, &frings, allow)?;
        f.setup_new_neighs()?;
        f.seq = efrags.len();
        efrags.push(f);
        match (0..cnrs).find(|&r| !fus_done[r]) {
            Some(r) => curr = r,
            None => break,
        }
    }
    Ok(())
}

/// `embedCisTransSystems`
fn embed_cis_trans<'a>(mol: &'a ROMol, efrags: &mut Vec<EmbeddedFrag<'a>>) -> Result<()> {
    for b in mol.bonds() {
        if b.bond_type() == BondType::Double
            && b.stereo() as u8 > BondStereo::StereoAny as u8
            && mol.ring_info().num_bond_rings(b.idx())? == 0
        {
            if b.stereo_atoms().len() != 2 {
                crate::rd_log_warning!("WARNING: bond found with stereo spec but no stereo atoms");
                continue;
            }
            let mut f = EmbeddedFrag::double_bond(mol, b.idx())?;
            f.setup_new_neighs()?;
            f.seq = efrags.len();
            efrags.push(f);
        }
    }
    Ok(())
}

/// `_findLargestFrag`: index of the first largest unfinished fragment.
fn find_largest_frag(efrags: &[EmbeddedFrag]) -> Option<usize> {
    let mut best = None;
    let mut msiz = 0;
    for (i, f) in efrags.iter().enumerate() {
        if !f.is_done() && f.size() > msiz {
            msiz = f.size();
            best = Some(i);
        }
    }
    best
}

/// Runs `f` on `efrags[idx]` with the other fragments, keeping upstream's
/// list order (fragments are ordered by creation).
fn with_frag<'a>(
    efrags: &mut Vec<EmbeddedFrag<'a>>,
    idx: usize,
    f: impl FnOnce(&mut EmbeddedFrag<'a>, &mut Vec<EmbeddedFrag<'a>>) -> Result<()>,
) -> Result<()> {
    let mut frag = efrags.remove(idx);
    let r = f(&mut frag, efrags);
    let pos = efrags
        .iter()
        .position(|o| o.seq > frag.seq)
        .unwrap_or(efrags.len());
    efrags.insert(pos, frag);
    r
}

/// `computeInitialCoords`
fn compute_initial_coords<'a>(
    mol: &'a ROMol,
    arings: &[Vec<i32>],
    coord_map: Option<&IntPoint2DMap>,
    use_ring_templates: bool,
) -> Result<Vec<EmbeddedFrag<'a>>> {
    let ranks: Vec<i32> = (0..mol.num_atoms())
        .map(|i| get_atom_depict_rank(mol, i))
        .collect();
    let mut efrags: Vec<EmbeddedFrag<'a>> = Vec::new();
    let mut pre_spec = false;
    if let Some(cm) = coord_map
        && cm.len() > 1
    {
        push_coord_frag(mol, &mut efrags, cm)?;
        pre_spec = true;
    }
    if !arings.is_empty() {
        embed_fused_systems(mol, arings, &mut efrags, coord_map, use_ring_templates)?;
    }
    embed_nontetrahedral(mol, &mut efrags, &ranks)?;
    embed_cis_trans(mol, &mut efrags)?;
    let mut done = vec![false; mol.num_atoms() as usize];
    for f in &efrags {
        for &a in f.embedded_atoms().keys() {
            done[a as usize] = true;
        }
    }
    let mut nratms: Vec<u32> = (0..mol.num_atoms())
        .filter(|&a| !done[a as usize])
        .collect();
    let mut mri = if pre_spec {
        Some(0)
    } else {
        find_largest_frag(&efrags)
    };
    while mri.is_some() || !nratms.is_empty() {
        let idx = match mri {
            Some(i) => i,
            None => {
                let na = mol.num_atoms() as i64;
                let mut mrank = i32::MAX as i64;
                let mut pick = 0;
                for (i, &a) in nratms.iter().enumerate() {
                    let rank = ranks[a as usize] as i64 * na + a as i64;
                    if rank < mrank {
                        mrank = rank;
                        pick = i;
                    }
                }
                let aid = nratms.remove(pick);
                let mut f = EmbeddedFrag::single_atom(aid, mol)?;
                f.seq = efrags.iter().map(|e| e.seq + 1).max().unwrap_or(0);
                efrags.push(f);
                efrags.len() - 1
            }
        };
        with_frag(&mut efrags, idx, |frag, others| {
            frag.mark_done();
            frag.expand_efrag(&mut nratms, others)
        })?;
        mri = find_largest_frag(&efrags);
    }
    Ok(efrags)
}

/// `_shiftCoords`
fn shift_coords(efrags: &mut [EmbeddedFrag]) {
    if efrags.is_empty() {
        return;
    }
    let boxes: Vec<(f64, f64, f64, f64)> = efrags.iter_mut().map(|f| f.compute_box()).collect();
    let (mut xmax, xmin, mut ymax, ymin) = boxes[0];
    for (i, f) in efrags.iter_mut().enumerate().skip(1) {
        let xshift = xmax + xmin <= ymax + ymin;
        let (xp, xn, yp, yn) = boxes[i];
        let mut shift = Point2D::new(0.0, 0.0);
        if xshift {
            shift.x = xmax + xn + 1.0;
            xmax += xp + xn + 1.0;
        } else {
            shift.y = ymax + yn + 1.0;
            ymax += yp + yn + 1.0;
        }
        f.translate(&shift);
    }
}

/// `copyCoordinate`: adds a conformer holding the fragments' coordinates to
/// `mol` and returns its id.
fn copy_coordinate(mol: &mut ROMol, efrags: Vec<EmbeddedFrag>, clear_confs: bool) -> Result<u32> {
    let mut conf = Conformer::new(mol.num_atoms() as usize);
    conf.set_3d(false);
    for f in &efrags {
        for (&aid, e) in f.embedded_atoms() {
            conf.set_atom_pos(aid, Point3D::new(e.loc.x, e.loc.y, 0.0))?;
        }
    }
    drop(efrags);
    if clear_confs {
        mol.clear_conformers();
        conf.set_id(0);
        mol.add_conformer(conf, false)?;
        Ok(0)
    } else {
        mol.add_conformer(conf, true)
    }
}

/// `RDDepict::compute2DCoords(mol, params)`: returns the conformer id.
pub fn compute_2d_coords(mol: &mut ROMol, params: &Compute2DCoordParameters) -> Result<u32> {
    if mol.needs_update_property_cache() {
        mol.update_property_cache(false)?;
    }
    let mut cp = mol.clone();
    let arings = symmetrize_sssr(&mut cp, SymmetrizeSSSRAlgorithm::Default, true, true, false)?;
    assign_stereochemistry(&mut cp, false, false, false)?;
    let mut efrags =
        compute_initial_coords(&cp, &arings, params.coord_map, params.use_ring_templates)?;
    // RDOXIDE_DEPICT_TRACE=1 prints the coordinates after each stage (to
    // compare with a C++ trace of the reference library)
    let trace = std::env::var_os("RDOXIDE_DEPICT_TRACE").is_some();
    let dump = |stage: &str, efrags: &[EmbeddedFrag]| {
        if trace {
            for (fi, f) in efrags.iter().enumerate() {
                for (a, e) in f.embedded_atoms() {
                    eprintln!("{stage} f{fi} {a} {:?} {:?}", e.loc.x, e.loc.y);
                }
            }
        }
    };
    dump("initial", &efrags);
    for f in efrags.iter_mut() {
        // either sample the 2D space by randomly flipping rotatable bonds
        // or flip only the bonds between colliding atoms, not both
        if params.n_samples > 0 && params.n_flips_per_sample > 0 {
            f.random_sample_flips_and_permutations(
                params.n_flips_per_sample,
                params.n_samples,
                params.sample_seed,
                None,
                0.0,
                params.permute_deg4_nodes,
            )?;
        } else {
            f.remove_collisions_bond_and_spiro_flip()?;
        }
    }
    dump("bondflip", &efrags);
    for i in 0..efrags.len() {
        efrags[i].remove_collisions_open_angles()?;
        dump("openangles", &efrags);
        efrags[i].remove_collisions_shorten_bonds()?;
        dump("shorten", &efrags);
        if params.use_path_angle_expansion {
            efrags[i].remove_collisions_path_angle_expansion()?;
        }
        dump("pathangle", &efrags);
    }
    if params.coord_map.is_none_or(|m| m.is_empty()) && params.canon_orient {
        for f in efrags.iter_mut() {
            f.canonicalize_orientation()?;
        }
    }
    shift_coords(&mut efrags);

    let cid = copy_coordinate(mol, efrags, params.clear_confs)?;
    if let Some(cm) = params.coord_map
        && cm.len() == 1
    {
        let (&aid, &target) = cm.iter().next().expect("one entry");
        let c = mol.conformer_mut(cid as i32)?;
        let cp0 = c.atom_pos(aid)?;
        let (dx, dy) = (target.x - cp0.x, target.y - cp0.y);
        for i in 0..c.num_atoms() {
            let mut p = c.atom_pos(i)?;
            p.x += dx;
            p.y += dy;
            c.set_atom_pos(i, p)?;
        }
    }
    Ok(cid)
}

/// `DepictorLocal::copySign`: the sign test has a tolerance, unlike
/// `f64::copysign`.
fn copy_sign(to: f64, from: f64, tol: f64) -> f64 {
    if from < -tol { -to.abs() } else { to.abs() }
}

#[derive(Default)]
struct ThetaBin {
    theta_avg: f64,
    theta_values: Vec<f64>,
}

/// `straightenDepiction`: rotates the conformer so that most bonds lie at
/// multiples of 30 degrees to the X axis.
pub fn straighten_depiction(mol: &mut ROMol, conf_id: i32, minimize_rotation: bool) -> Result<()> {
    use crate::geometry::transform3d::Transform3D;
    use crate::graphmol::mol_transforms::transform_conformer;

    const ALMOST_ZERO: f64 = 1.0e-5;
    const INCR_DEG: f64 = 30.0;
    const HALF_INCR_DEG: f64 = 0.5 * INCR_DEG;
    const QUARTER_INCR_DEG: f64 = 0.25 * INCR_DEG;
    if mol.num_bonds() == 0 {
        return Ok(());
    }
    // bins in insertion order; libc++'s unordered_map iterates the
    // (distinct-bucket) nodes newest-first, which we emulate below
    let mut keys: Vec<i32> = Vec::new();
    let mut bins: Vec<ThetaBin> = Vec::new();
    {
        let pos = mol.conformer(conf_id)?.positions();
        for b in mol.bonds() {
            let bv = pos[b.begin_atom_idx() as usize] - pos[b.end_atom_idx() as usize];
            let x = if bv.x < 0.0 {
                bv.x.min(-ALMOST_ZERO)
            } else {
                bv.x.max(ALMOST_ZERO)
            };
            let theta = (bv.y / x).atan().to_degrees();
            let mut d_theta = (-theta) % INCR_DEG;
            if d_theta.abs() > HALF_INCR_DEG {
                d_theta -= copy_sign(INCR_DEG, d_theta, ALMOST_ZERO);
            }
            let key = (d_theta + copy_sign(0.5, d_theta, ALMOST_ZERO)) as i32;
            let idx = keys.iter().position(|&k| k == key).unwrap_or_else(|| {
                keys.push(key);
                bins.push(ThetaBin::default());
                bins.len() - 1
            });
            bins[idx].theta_avg += d_theta;
            bins[idx].theta_values.push(theta);
        }
    }
    bins.reverse();
    let mut d_theta_smallest = f64::MAX;
    for bin in &mut bins {
        bin.theta_avg /= bin.theta_values.len() as f64;
        if bin.theta_avg.abs() < d_theta_smallest.abs() {
            d_theta_smallest = bin.theta_avg;
        }
    }
    // std::max_element keeps the first of equivalent elements
    let mut best = &bins[0];
    for bin in &bins[1..] {
        if best.theta_values.len() < bin.theta_values.len()
            || (best.theta_values.len() == bin.theta_values.len()
                && best.theta_avg.abs() > bin.theta_avg.abs())
        {
            best = bin;
        }
    }
    let mut d_theta_min = best.theta_avg;
    if !minimize_rotation {
        let mut count60_vs30 = [0u32; 2];
        for &theta in &best.theta_values {
            let abs_theta = (theta + d_theta_min).abs();
            if abs_theta < ALMOST_ZERO {
                continue;
            }
            count60_vs30[((abs_theta + 0.5) / INCR_DEG) as usize % 2] += 1;
        }
        if count60_vs30[0] > count60_vs30[1] {
            d_theta_min -= copy_sign(INCR_DEG, d_theta_min, ALMOST_ZERO);
        }
    } else if d_theta_smallest.abs() < ALMOST_ZERO
        || (d_theta_smallest.abs() < d_theta_min.abs() && d_theta_min.abs() > QUARTER_INCR_DEG)
    {
        d_theta_min = d_theta_smallest;
    }
    if d_theta_min.abs() > ALMOST_ZERO {
        let mut trans = Transform3D::identity();
        trans.set_rotation(d_theta_min.to_radians(), &Point3D::new(0.0, 0.0, 1.0));
        transform_conformer(mol.conformer_mut(conf_id)?, &trans);
    }
    Ok(())
}

/// `normalizeDepiction`: centres the conformer (optionally aligning its
/// principal axes) and scales it to the standard bond length. Returns the
/// scale factor applied, or -1 if none could be determined.
pub fn normalize_depiction(
    mol: &mut ROMol,
    conf_id: i32,
    canonicalize: i32,
    mut scale_factor: f64,
) -> Result<f64> {
    use crate::geometry::transform3d::Transform3D;
    use crate::graphmol::mol_transforms::{
        compute_canonical_transform, compute_centroid, get_bond_length, transform_conformer,
    };

    const SCALE_FACTOR_THRESHOLD: f64 = 1.0e-5;
    if mol.num_bonds() == 0 {
        return Ok(-1.0);
    }
    let conf = mol.conformer(conf_id)?;
    if scale_factor < 0.0 {
        const RDKIT_BOND_LEN: f64 = 1.5;
        let mut most_common = -1;
        let mut max_count = 0u32;
        let mut binned: std::collections::HashMap<i32, u32> = std::collections::HashMap::new();
        for b in mol.bonds() {
            let len =
                (get_bond_length(conf, b.begin_atom_idx(), b.end_atom_idx())? * 10.0 + 0.5) as i32;
            let count = binned.entry(len).or_insert(0);
            *count += 1;
            if *count > max_count {
                max_count = *count;
                most_common = len;
            }
        }
        if most_common > 0 {
            scale_factor = RDKIT_BOND_LEN / (f64::from(most_common) * 0.1);
        }
    }
    let ctd = compute_centroid(mol, conf, true, None)?;
    let mut canon_trans = if canonicalize != 0 {
        let (mut t, _) = compute_canonical_transform(mol, conf, Some(ctd), false, true)?;
        if canonicalize < 0 {
            let mut rotate90 = Transform3D::identity();
            rotate90.set_rotation_cs(0.0, 1.0, &Point3D::new(0.0, 0.0, 1.0));
            t.mul_assign(&rotate90);
        }
        t
    } else {
        let mut t = Transform3D::identity();
        t.set_translation(&-ctd);
        t
    };
    let sane = scale_factor > SCALE_FACTOR_THRESHOLD;
    if sane && (scale_factor - 1.0).abs() > SCALE_FACTOR_THRESHOLD {
        let mut trans = Transform3D::identity();
        trans.set_val(0, 0, scale_factor);
        trans.set_val(1, 1, scale_factor);
        trans.mul_assign(&canon_trans);
        canon_trans = trans;
    }
    transform_conformer(mol.conformer_mut(conf_id)?, &canon_trans);
    Ok(if sane { scale_factor } else { -1.0 })
}

/// `compute2DCoordsMimicDistMat`: 2D coordinates whose inter-atom distances
/// mimic `dmat` (the packed lower triangle, negative entries are ignored),
/// found by random sampling of the rotatable bonds. `weight_dist_mat` in
/// [0, 1] balances mimicking `dmat` against spreading the structure out.
#[allow(clippy::too_many_arguments)]
pub fn compute_2d_coords_mimic_dist_mat(
    mol: &mut ROMol,
    dmat: &[f64],
    canon_orient: bool,
    clear_confs: bool,
    weight_dist_mat: f64,
    n_flips_per_sample: u32,
    n_samples: u32,
    sample_seed: i32,
    permute_deg4_nodes: bool,
) -> Result<u32> {
    let mut cp = mol.clone();
    let arings = symmetrize_sssr(&mut cp, SymmetrizeSSSRAlgorithm::Default, true, true, false)?;
    assign_stereochemistry(&mut cp, false, false, false)?;
    let mut efrags = compute_initial_coords(&cp, &arings, None, false)?;
    for f in efrags.iter_mut() {
        f.random_sample_flips_and_permutations(
            n_flips_per_sample,
            n_samples,
            sample_seed,
            Some(dmat),
            weight_dist_mat,
            permute_deg4_nodes,
        )?;
    }
    if canon_orient {
        for f in efrags.iter_mut() {
            f.canonicalize_orientation()?;
        }
    }
    shift_coords(&mut efrags);
    copy_coordinate(mol, efrags, clear_confs)
}
