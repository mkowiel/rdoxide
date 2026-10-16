//! Encoding the shapes of molecules on grids and comparing them
//! (`MolShapes` of RDKit `Code/GraphMol/ShapeHelpers`).
//!
//! Ported from `ShapeEncoder.cpp` and `ShapeUtils.cpp` (2026.09.1). The
//! functions that take a conformer in upstream need the molecule that owns it
//! (for the atoms), so here they take both.
//!
//! Copyright (C) 2005-2006 Rational Discovery LLC. Licensed under the BSD
//! 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::conformer::Conformer;
use super::mol_transforms::{compute_canonical_transform, compute_centroid};
use super::periodic_table::PeriodicTable;
use super::romol::ROMol;
use crate::Result;
use crate::geometry::grid::{self, UniformGrid3D};
use crate::geometry::{Point3D, Transform3D};

/// `computeConfBox`: the lower-left and upper-right corners of the box around
/// the conformer, enlarged by `padding`.
pub fn compute_conf_box(
    conf: &Conformer,
    trans: Option<&Transform3D>,
    padding: f64,
) -> Result<(Point3D, Point3D)> {
    let (mut lo, mut hi) = (
        Point3D::new(1.0e8, 1.0e8, 1.0e8),
        Point3D::new(-1.0e8, -1.0e8, -1.0e8),
    );
    for i in 0..conf.num_atoms() {
        let mut loc = conf.atom_pos(i)?;
        if let Some(t) = trans {
            loc = t.transform_point(&loc);
        }
        hi.x = hi.x.max(loc.x);
        lo.x = lo.x.min(loc.x);
        hi.y = hi.y.max(loc.y);
        lo.y = lo.y.min(loc.y);
        hi.z = hi.z.max(loc.z);
        lo.z = lo.z.min(loc.z);
    }
    let pad = Point3D::new(padding, padding, padding);
    Ok((lo - pad, hi + pad))
}

/// `computeConfDimsAndOffset`: (dimensions, offset) of the box.
pub fn compute_conf_dims_and_offset(
    conf: &Conformer,
    trans: Option<&Transform3D>,
    padding: f64,
) -> Result<(Point3D, Point3D)> {
    let (offset, hi) = compute_conf_box(conf, trans, padding)?;
    Ok((hi - offset, offset))
}

/// `getConfDimensions`: the dimensions of a box centered on `center` (the
/// centroid by default) that holds the conformer.
pub fn get_conf_dimensions(
    mol: &ROMol,
    conf: &Conformer,
    padding: f64,
    center: Option<&Point3D>,
    ignore_hs: bool,
) -> Result<[f64; 3]> {
    let (mut lb, mut rb) = compute_conf_box(conf, None, padding)?;
    let c = match center {
        Some(c) => *c,
        None => compute_centroid(mol, conf, ignore_hs, None)?,
    };
    rb -= c;
    lb -= c;
    lb *= -1.0;
    Ok([
        2.0 * rb.x.max(lb.x),
        2.0 * rb.y.max(lb.y),
        2.0 * rb.z.max(lb.z),
    ])
}

/// `computeUnionBox`
pub fn compute_union_box(
    left_bottom1: &Point3D,
    right_top1: &Point3D,
    left_bottom2: &Point3D,
    right_top2: &Point3D,
) -> (Point3D, Point3D) {
    (
        Point3D::new(
            left_bottom1.x.min(left_bottom2.x),
            left_bottom1.y.min(left_bottom2.y),
            left_bottom1.z.min(left_bottom2.z),
        ),
        Point3D::new(
            right_top1.x.max(right_top2.x),
            right_top1.y.max(right_top2.y),
            right_top1.z.max(right_top2.z),
        ),
    )
}

/// `EncodeShape`: writes the shape of the conformer (atom spheres of the
/// van der Waals radius times `vdw_scale`, plus `max_layers` shells of
/// `step_size`) onto the grid. The molecule itself is not changed; `trans`
/// is applied to the atom positions.
#[allow(clippy::too_many_arguments)]
pub fn encode_shape(
    mol: &ROMol,
    conf: &Conformer,
    grid: &mut UniformGrid3D,
    trans: Option<&Transform3D>,
    vdw_scale: f64,
    step_size: f64,
    max_layers: i32,
    ignore_hs: bool,
) -> Result<()> {
    let table = PeriodicTable::get();
    for atom in mol.atoms() {
        let anum = atom.atomic_num();
        if anum == 1 && ignore_hs {
            continue;
        }
        let rad = table.get_rvdw(anum)?;
        let mut loc = conf.atom_pos(atom.idx())?;
        if let Some(t) = trans {
            loc = t.transform_point(&loc);
        }
        grid.set_sphere_occupancy(&loc, vdw_scale * rad, step_size, max_layers, true)?;
    }
    Ok(())
}

/// The options shared by the shape comparison functions.
#[derive(Debug, Clone, Copy)]
pub struct ShapeGridOptions {
    pub grid_spacing: f64,
    /// `DiscreteValueType`: 0 = one bit .. 4 = sixteen bits per grid point
    pub bits_per_point: u32,
    pub vdw_scale: f64,
    pub step_size: f64,
    pub max_layers: i32,
    pub ignore_hs: bool,
}

impl Default for ShapeGridOptions {
    fn default() -> Self {
        ShapeGridOptions {
            grid_spacing: 0.5,
            bits_per_point: 1,
            vdw_scale: 0.8,
            step_size: 0.25,
            max_layers: -1,
            ignore_hs: true,
        }
    }
}

/// `createShapeGrids`: both shapes encoded on grids that cover the union of
/// their boxes, in the frame of the first conformer's principal axes.
fn create_shape_grids(
    mol1: &ROMol,
    conf1: &Conformer,
    mol2: &ROMol,
    conf2: &Conformer,
    o: &ShapeGridOptions,
) -> Result<(UniformGrid3D, UniformGrid3D)> {
    // first compute the canonical transform for the first conformer so that
    // we get an efficient grid. The same transform is applied to the second
    // conformer, so the relative orientations/positions are preserved
    let (trans, _) = compute_canonical_transform(mol1, conf1, None, false, true)?;
    let (lb1, rt1) = compute_conf_box(conf1, Some(&trans), 2.5)?;
    let (lb2, rt2) = compute_conf_box(conf2, Some(&trans), 2.5)?;
    let (lb, mut rt) = compute_union_box(&lb1, &rt1, &lb2, &rt2);
    // rt now has the dimensions of the grid
    rt -= lb;
    let mut g1 = UniformGrid3D::new(rt.x, rt.y, rt.z, o.grid_spacing, o.bits_per_point, Some(lb))?;
    let mut g2 = UniformGrid3D::new(rt.x, rt.y, rt.z, o.grid_spacing, o.bits_per_point, Some(lb))?;
    encode_shape(
        mol1,
        conf1,
        &mut g1,
        Some(&trans),
        o.vdw_scale,
        o.step_size,
        o.max_layers,
        o.ignore_hs,
    )?;
    encode_shape(
        mol2,
        conf2,
        &mut g2,
        Some(&trans),
        o.vdw_scale,
        o.step_size,
        o.max_layers,
        o.ignore_hs,
    )?;
    Ok((g1, g2))
}

/// `tverskyIndex` of the shapes of two conformers (which are assumed to be
/// aligned already).
pub fn tversky_index(
    mol1: &ROMol,
    conf1: &Conformer,
    mol2: &ROMol,
    conf2: &Conformer,
    alpha: f64,
    beta: f64,
    o: &ShapeGridOptions,
) -> Result<f64> {
    let (g1, g2) = create_shape_grids(mol1, conf1, mol2, conf2, o)?;
    grid::tversky_index(&g1, &g2, alpha, beta)
}

/// `tanimotoDistance` of the shapes of two conformers.
pub fn tanimoto_distance(
    mol1: &ROMol,
    conf1: &Conformer,
    mol2: &ROMol,
    conf2: &Conformer,
    o: &ShapeGridOptions,
) -> Result<f64> {
    let (g1, g2) = create_shape_grids(mol1, conf1, mol2, conf2, o)?;
    grid::tanimoto_distance(&g1, &g2)
}

/// `protrudeDistance` of the shapes of two conformers; with
/// `allow_reordering` the smaller shape is protruded out of the larger.
pub fn protrude_distance(
    mol1: &ROMol,
    conf1: &Conformer,
    mol2: &ROMol,
    conf2: &Conformer,
    o: &ShapeGridOptions,
    allow_reordering: bool,
) -> Result<f64> {
    let (g1, g2) = create_shape_grids(mol1, conf1, mol2, conf2, o)?;
    if allow_reordering && g2.occupancy_vect().total_val() < g1.occupancy_vect().total_val() {
        grid::protrude_distance(&g2, &g1)
    } else {
        grid::protrude_distance(&g1, &g2)
    }
}

/// The molecule versions of the comparison functions take conformer ids.
pub mod mol {
    use super::*;

    pub fn tversky_index(
        mol1: &ROMol,
        mol2: &ROMol,
        alpha: f64,
        beta: f64,
        conf_id1: i32,
        conf_id2: i32,
        o: &ShapeGridOptions,
    ) -> Result<f64> {
        super::tversky_index(
            mol1,
            mol1.conformer(conf_id1)?,
            mol2,
            mol2.conformer(conf_id2)?,
            alpha,
            beta,
            o,
        )
    }

    pub fn tanimoto_distance(
        mol1: &ROMol,
        mol2: &ROMol,
        conf_id1: i32,
        conf_id2: i32,
        o: &ShapeGridOptions,
    ) -> Result<f64> {
        super::tanimoto_distance(
            mol1,
            mol1.conformer(conf_id1)?,
            mol2,
            mol2.conformer(conf_id2)?,
            o,
        )
    }

    /// (upstream's molecule version ignores the grid spacing it is given and
    /// always uses 0.5)
    pub fn protrude_distance(
        mol1: &ROMol,
        mol2: &ROMol,
        conf_id1: i32,
        conf_id2: i32,
        o: &ShapeGridOptions,
        allow_reordering: bool,
    ) -> Result<f64> {
        let o = ShapeGridOptions {
            grid_spacing: 0.5,
            ..*o
        };
        super::protrude_distance(
            mol1,
            mol1.conformer(conf_id1)?,
            mol2,
            mol2.conformer(conf_id2)?,
            &o,
            allow_reordering,
        )
    }
}
