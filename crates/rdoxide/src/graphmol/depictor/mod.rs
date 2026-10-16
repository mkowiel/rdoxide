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
        if params.n_samples > 0 && params.n_flips_per_sample > 0 {
            return Err(RDKitError::Runtime(
                "random flip sampling in compute2DCoords is not supported by rdoxide yet".into(),
            ));
        }
        f.remove_collisions_bond_and_spiro_flip()?;
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

    // copyCoordinate
    let mut conf = Conformer::new(mol.num_atoms() as usize);
    conf.set_3d(false);
    for f in &efrags {
        for (&aid, e) in f.embedded_atoms() {
            conf.set_atom_pos(aid, Point3D::new(e.loc.x, e.loc.y, 0.0))?;
        }
    }
    drop(efrags);
    let cid = if params.clear_confs {
        mol.clear_conformers();
        conf.set_id(0);
        mol.add_conformer(conf, false)?;
        0
    } else {
        mol.add_conformer(conf, true)?
    };
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
