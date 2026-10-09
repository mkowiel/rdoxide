//! Uniform three-dimensional grids (upstream `Geometry/UniformGrid3D.cpp`,
//! `UniformRealValueGrid3D.cpp` and `GridUtils.cpp`).
//!
//! [`UniformGrid3D`] stores small unsigned integers in a
//! [`DiscreteValueVect`], [`UniformRealValueGrid3D`] stores doubles in a
//! [`RealValueVect`]. Grid points are numbered `z * nx * ny + y * nx + x`.
//!
//! Copyright (C) 2005-2024 Greg Landrum, Novartis Institutes for BioMedical
//! Research and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::fmt::Write as _;

use super::Point3D;
use crate::datastructs::value_vects::{
    DiscreteValueVect, RealValueVect, compute_l1_norm_discrete, cxx_max, cxx_min,
};
use crate::error::{RDKitError, Result};

/// Pickle version of both grid kinds.
const GRID_PICKLE_VERSION: i32 = 0x1;

/// `GridException` ("Center outside the grid boundary").
fn grid_error(msg: &str) -> RDKitError {
    RDKitError::Runtime(msg.to_string())
}

fn index_error(idx: u32) -> RDKitError {
    RDKitError::Index(format!("Index Error: {}", idx as i32))
}

/// The dimensions, spacing and offset shared by both grid kinds.
#[derive(Clone, Debug, PartialEq)]
struct GridGeometry {
    num_x: u32,
    num_y: u32,
    num_z: u32,
    spacing: f64,
    offset: Point3D,
}

impl GridGeometry {
    /// `initGrid`: the number of points per axis is `floor(dim / spacing + 0.5)`.
    fn new(dim_x: f64, dim_y: f64, dim_z: f64, spacing: f64, offset: Point3D) -> Result<Self> {
        if dim_x <= 0.0 {
            return Err(RDKitError::precondition("Invalid x-dimension for grid"));
        }
        if dim_y <= 0.0 {
            return Err(RDKitError::precondition("Invalid y-dimension for grid"));
        }
        if dim_z <= 0.0 {
            return Err(RDKitError::precondition("Invalid z-dimension for grid"));
        }
        if spacing <= 0.0 {
            return Err(RDKitError::precondition("Invalid spacing for grid"));
        }
        Ok(GridGeometry {
            num_x: (dim_x / spacing + 0.5).floor() as u32,
            num_y: (dim_y / spacing + 0.5).floor() as u32,
            num_z: (dim_z / spacing + 0.5).floor() as u32,
            spacing,
            offset,
        })
    }

    /// Dimensions of the grid in grid units (what the constructor was given).
    fn dims(&self) -> (f64, f64, f64) {
        (
            f64::from(self.num_x) * self.spacing,
            f64::from(self.num_y) * self.spacing,
            f64::from(self.num_z) * self.spacing,
        )
    }

    fn size(&self) -> u32 {
        self.num_x.wrapping_mul(self.num_y).wrapping_mul(self.num_z)
    }

    fn grid_index(&self, xi: u32, yi: u32, zi: u32) -> i32 {
        if xi >= self.num_x || yi >= self.num_y || zi >= self.num_z {
            return -1;
        }
        (zi * self.num_x * self.num_y + yi * self.num_x + xi) as i32
    }

    fn grid_indices(&self, idx: u32) -> Result<(u32, u32, u32)> {
        if idx >= self.size() {
            return Err(index_error(idx));
        }
        Ok((
            idx % self.num_x,
            (idx % (self.num_x * self.num_y)) / self.num_x,
            idx / (self.num_x * self.num_y),
        ))
    }

    /// Index of the grid point closest to `point`, -1 outside the grid.
    fn point_index(&self, point: &Point3D) -> i32 {
        let t = (*point - self.offset) / self.spacing;
        let (xi, yi, zi) = (
            (t.x + 0.5).floor() as i32,
            (t.y + 0.5).floor() as i32,
            (t.z + 0.5).floor() as i32,
        );
        if xi < 0 || xi >= self.num_x as i32 {
            return -1;
        }
        if yi < 0 || yi >= self.num_y as i32 {
            return -1;
        }
        if zi < 0 || zi >= self.num_z as i32 {
            return -1;
        }
        (zi as u32 * self.num_x * self.num_y + yi as u32 * self.num_x + xi as u32) as i32
    }

    fn point_loc(&self, id: u32) -> Result<Point3D> {
        if id >= self.size() {
            return Err(index_error(id));
        }
        let loc = Point3D::new(
            f64::from(id % self.num_x) * self.spacing,
            f64::from((id % (self.num_x * self.num_y)) / self.num_x) * self.spacing,
            f64::from(id / (self.num_x * self.num_y)) * self.spacing,
        );
        Ok(loc + self.offset)
    }

    /// `compareParams` with the given tolerances.
    fn compatible_with(&self, other: &Self, spacing_tol: f64, offset_tol: f64) -> bool {
        if self.num_x != other.num_x || self.num_y != other.num_y || self.num_z != other.num_z {
            return false;
        }
        if (self.spacing - other.spacing).abs() > spacing_tol {
            return false;
        }
        (self.offset - other.offset).length_sq() <= offset_tol
    }

    /// The pickle header shared by both grid kinds.
    fn write_header(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&(-GRID_PICKLE_VERSION).to_le_bytes());
        for n in [self.num_x, self.num_y, self.num_z] {
            out.extend_from_slice(&n.to_le_bytes());
        }
        for v in [self.spacing, self.offset.x, self.offset.y, self.offset.z] {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }

    /// Reads the pickle header; returns the geometry and the position after
    /// the length-prefixed storage pickle.
    fn read_header<'a>(pkl: &'a [u8], what: &str) -> Result<(Self, &'a [u8])> {
        let mut pos = 0;
        let mut take = |n: usize| -> [u8; 8] {
            // a short pickle leaves the (zero-initialized) value untouched
            let mut b = [0u8; 8];
            if let Some(s) = pkl.get(pos..pos + n) {
                b[..n].copy_from_slice(s);
            }
            pos += n;
            b
        };
        let vers = i32::from_le_bytes(take(4)[..4].try_into().expect("4 bytes")).wrapping_neg();
        if vers != GRID_PICKLE_VERSION {
            return Err(RDKitError::Value(format!("bad version in {what} pickle")));
        }
        let mut u32s = [0u32; 3];
        for n in &mut u32s {
            *n = u32::from_le_bytes(take(4)[..4].try_into().expect("4 bytes"));
        }
        let mut f64s = [0f64; 4];
        for v in &mut f64s {
            *v = f64::from_le_bytes(take(8));
        }
        let size = u32::from_le_bytes(take(4)[..4].try_into().expect("4 bytes")) as usize;
        let storage = pkl.get(pos..).unwrap_or(&[]);
        let storage = &storage[..size.min(storage.len())];
        Ok((
            GridGeometry {
                num_x: u32s[0],
                num_y: u32s[1],
                num_z: u32s[2],
                spacing: f64s[0],
                offset: Point3D::new(f64s[1], f64s[2], f64s[3]),
            },
            storage,
        ))
    }
}

/// A grid of unsigned integers of 1, 2, 4, 8 or 16 bits per point.
#[derive(Clone, Debug)]
pub struct UniformGrid3D {
    geometry: GridGeometry,
    storage: DiscreteValueVect,
}

const GRID_OFFSET_TOL: f64 = 1.0e-8;
const GRID_SPACING_TOL: f64 = 1.0e-8;

impl UniformGrid3D {
    /// `UniformGrid3D(dimX, dimY, dimZ, spacing, valType, offset)`; the default
    /// offset centers the grid on the origin.
    pub fn new(
        dim_x: f64,
        dim_y: f64,
        dim_z: f64,
        spacing: f64,
        val_type: u32,
        offset: Option<Point3D>,
    ) -> Result<Self> {
        let offset =
            offset.unwrap_or_else(|| Point3D::new(-0.5 * dim_x, -0.5 * dim_y, -0.5 * dim_z));
        let geometry = GridGeometry::new(dim_x, dim_y, dim_z, spacing, offset)?;
        let storage = DiscreteValueVect::new(val_type, geometry.size());
        Ok(UniformGrid3D { geometry, storage })
    }

    /// Reads a grid from its pickle (`toString()`).
    pub fn from_binary(pkl: &[u8]) -> Result<Self> {
        let (geometry, storage) = GridGeometry::read_header(pkl, "UniformGrid3D")?;
        Ok(UniformGrid3D {
            geometry,
            storage: DiscreteValueVect::from_binary(storage)?,
        })
    }

    /// `toString()`
    pub fn to_binary(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.geometry.write_header(&mut out);
        let storage = self.storage.to_binary();
        out.extend_from_slice(&(storage.len() as u32).to_le_bytes());
        out.extend_from_slice(&storage);
        out
    }

    pub fn num_x(&self) -> u32 {
        self.geometry.num_x
    }
    pub fn num_y(&self) -> u32 {
        self.geometry.num_y
    }
    pub fn num_z(&self) -> u32 {
        self.geometry.num_z
    }
    pub fn offset(&self) -> Point3D {
        self.geometry.offset
    }
    pub fn spacing(&self) -> f64 {
        self.geometry.spacing
    }
    /// Number of grid points.
    pub fn size(&self) -> u32 {
        self.geometry.size()
    }
    pub fn occupancy_vect(&self) -> &DiscreteValueVect {
        &self.storage
    }

    pub fn grid_index(&self, xi: u32, yi: u32, zi: u32) -> i32 {
        self.geometry.grid_index(xi, yi, zi)
    }
    pub fn grid_indices(&self, idx: u32) -> Result<(u32, u32, u32)> {
        self.geometry.grid_indices(idx)
    }
    pub fn grid_point_index(&self, point: &Point3D) -> i32 {
        self.geometry.point_index(point)
    }
    pub fn grid_point_loc(&self, id: u32) -> Result<Point3D> {
        self.geometry.point_loc(id)
    }

    /// The value at the grid point closest to `point`, -1 outside the grid.
    pub fn val_at_point(&self, point: &Point3D) -> Result<i64> {
        let id = self.grid_point_index(point);
        if id < 0 {
            return Ok(-1);
        }
        Ok(i64::from(self.storage.get_val(id as u32)?))
    }
    pub fn val(&self, id: u32) -> Result<u32> {
        self.storage.get_val(id)
    }
    pub fn set_val(&mut self, id: u32, val: u32) -> Result<()> {
        self.storage.set_val(id, val)
    }
    /// Sets the value at the grid point closest to `point`; nothing happens
    /// outside the grid.
    pub fn set_val_at_point(&mut self, point: &Point3D, val: u32) -> Result<()> {
        let id = self.grid_point_index(point);
        if id < 0 {
            return Ok(());
        }
        self.storage.set_val(id as u32, val)
    }

    /// `compareParams`
    pub fn compare_params(&self, other: &Self) -> bool {
        self.geometry
            .compatible_with(&other.geometry, GRID_SPACING_TOL, GRID_OFFSET_TOL)
    }

    fn check_compatible(&self, other: &Self) -> Result<()> {
        if self.compare_params(other) {
            Ok(())
        } else {
            Err(RDKitError::precondition("incompatible grids"))
        }
    }

    /// `operator|=`
    pub fn union_assign(&mut self, other: &Self) -> Result<()> {
        self.check_compatible(other)?;
        self.storage = self.storage.or(&other.storage)?;
        Ok(())
    }
    /// `operator&=`
    pub fn intersect_assign(&mut self, other: &Self) -> Result<()> {
        self.check_compatible(other)?;
        self.storage = self.storage.and(&other.storage)?;
        Ok(())
    }
    /// `operator+=`
    pub fn add_assign(&mut self, other: &Self) -> Result<()> {
        self.check_compatible(other)?;
        self.storage.add_assign(&other.storage)
    }
    /// `operator-=`
    pub fn sub_assign(&mut self, other: &Self) -> Result<()> {
        self.check_compatible(other)?;
        self.storage.sub_assign(&other.storage)
    }

    /// `setSphereOccupancy`: encodes a sphere and `max_num_layers` shells of
    /// `step_size` around it with decreasing occupancy values.
    pub fn set_sphere_occupancy(
        &mut self,
        center: &Point3D,
        radius: f64,
        step_size: f64,
        max_num_layers: i32,
        ignore_out_of_bound: bool,
    ) -> Result<()> {
        if self.grid_point_index(center) == -1 {
            return if ignore_out_of_bound {
                Ok(())
            } else {
                Err(grid_error("Center outside the grid boundary"))
            };
        }
        // the center in grid coordinates
        let g_pt = (*center - self.geometry.offset) / self.geometry.spacing;
        let bits_per_val = 1u32 << self.storage.value_type();
        let max_val = (1u32 << bits_per_val) - 1;
        let mut n_layers = max_val;
        let mut val_step = 1u32;
        if max_num_layers > 0 && max_num_layers <= n_layers as i32 {
            n_layers = max_num_layers as u32;
            val_step = (max_val + 1) / n_layers;
        }
        let bg_rad = radius / self.geometry.spacing; // base radius in grid coords
        let g_step_size = step_size / self.geometry.spacing;
        let g_radius = bg_rad + f64::from(n_layers) * g_step_size; // largest radius
        let g_rad2 = g_radius * g_radius;
        let bg_rad2 = bg_rad * bg_rad;
        let (nx, ny, nz) = (
            self.geometry.num_x as i32,
            self.geometry.num_y as i32,
            self.geometry.num_z as i32,
        );
        let (xmax, xmin) = (
            (g_pt.x + g_radius).floor() as i32,
            (g_pt.x - g_radius).ceil() as i32,
        );
        let (ymax, ymin) = (
            (g_pt.y + g_radius).floor() as i32,
            (g_pt.y - g_radius).ceil() as i32,
        );
        let (zmax, zmin) = (
            (g_pt.z + g_radius).floor() as i32,
            (g_pt.z - g_radius).ceil() as i32,
        );

        for k in zmin..=zmax {
            if !(0..nz).contains(&k) {
                continue;
            }
            let dz = f64::from(k) - g_pt.z;
            let dz2 = dz * dz;
            let pt_id1 = k * nx * ny;
            for j in ymin..=ymax {
                if !(0..ny).contains(&j) {
                    continue;
                }
                let dy = f64::from(j) - g_pt.y;
                let dy2z2 = dy * dy + dz2;
                if dy2z2 >= g_rad2 {
                    continue;
                }
                let pt_id2 = pt_id1 + j * nx;
                for i in xmin..=xmax {
                    if !(0..nx).contains(&i) {
                        continue;
                    }
                    let id = (pt_id2 + i) as u32;
                    let oval = self.storage.get_val(id)?;
                    if oval >= max_val {
                        continue; // already at maxVal, nothing to change
                    }
                    let dx = f64::from(i) - g_pt.x;
                    let d2 = dx * dx + dy2z2;
                    if d2 >= g_rad2 {
                        continue;
                    }
                    let val = if d2 < bg_rad2 {
                        max_val
                    } else {
                        let d = d2.sqrt();
                        let value_change = ((d - bg_rad) / g_step_size + 1.0) as u32 * val_step;
                        max_val.saturating_sub(value_change)
                    };
                    if val > oval {
                        self.storage.set_val(id, val)?;
                    }
                }
            }
        }
        Ok(())
    }

    /// `writeGridToStream`: the grid in GRD format.
    pub fn to_grd(&self) -> Result<String> {
        let (dim_x, dim_y, dim_z) = (
            self.num_x() as i32,
            self.num_y() as i32,
            self.num_z() as i32,
        );
        let (len_x, len_y, len_z) = self.geometry.dims();
        let offset = self.offset() / self.spacing();
        let mut out = String::from("Grid file representing a Shape \n\n");
        let _ = writeln!(
            out,
            "{} {} {} 90.0 90.0 90.0",
            cxx_f64(len_x),
            cxx_f64(len_y),
            cxx_f64(len_z)
        );
        let _ = writeln!(out, "{} {} {}", dim_x - 1, dim_y - 1, dim_z - 1);
        let out_x1 = (offset.x + 0.5).floor() as i32;
        let out_x2 = out_x1 + (dim_x - 1);
        // upstream shifts the y and z ranges by one more than x, which usually
        // displays the grid closer to the molecule in PyMol
        let out_y1 = (offset.y + 0.5).floor() as i32;
        let out_y2 = out_y1 + dim_y;
        let out_z1 = (offset.z + 0.5).floor() as i32;
        let out_z2 = out_z1 + dim_z;
        let _ = writeln!(
            out,
            "1 {out_x1} {out_x2} {out_y1} {out_y2} {out_z1} {out_z2}"
        );
        for i in 0..self.size() {
            let _ = writeln!(out, "{}", cxx_f64(f64::from(self.val(i)?)));
        }
        Ok(out)
    }

    /// `writeGridToFile`
    pub fn write_to_file(&self, filename: &str) -> Result<()> {
        let text = self.to_grd()?;
        // upstream does not check whether the file could be opened
        let _ = std::fs::write(filename, text);
        Ok(())
    }
}

/// `operator<<(std::ostream&, double)` with the default precision (6
/// significant digits, `%g` style).
fn cxx_f64(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let exp = v.abs().log10().floor() as i32;
    if !(-4..6).contains(&exp) {
        let s = format!("{v:.5e}");
        let (mant, e) = s.split_once('e').expect("exponent");
        let mant = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        let e: i32 = e.parse().expect("exponent value");
        return format!("{mant}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs());
    }
    let decimals = (5 - exp).max(0) as usize;
    let s = format!("{v:.decimals$}");
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

/// `tverskyIndex`
pub fn tversky_index(g1: &UniformGrid3D, g2: &UniformGrid3D, alpha: f64, beta: f64) -> Result<f64> {
    let (v1, v2, dist) = shapes(g1, g2)?;
    let (totv1, totv2) = (f64::from(v1.total_val()), f64::from(v2.total_val()));
    let inter = 0.5 * (totv1 + totv2 - f64::from(dist));
    Ok(inter / (alpha * (totv1 - inter) + beta * (totv2 - inter) + inter))
}

/// `tanimotoDistance`
pub fn tanimoto_distance(g1: &UniformGrid3D, g2: &UniformGrid3D) -> Result<f64> {
    let (v1, v2, dist) = shapes(g1, g2)?;
    let inter = 0.5 * (f64::from(v1.total_val()) + f64::from(v2.total_val()) - f64::from(dist));
    Ok(f64::from(dist) / (f64::from(dist) + inter))
}

/// `protrudeDistance`
pub fn protrude_distance(g1: &UniformGrid3D, g2: &UniformGrid3D) -> Result<f64> {
    let (v1, v2, total_protrude) = shapes(g1, g2)?;
    let (totv1, totv2) = (v1.total_val(), v2.total_val());
    let intersect_volume = (totv1.wrapping_add(totv2).wrapping_sub(total_protrude)) / 2;
    Ok((f64::from(totv1) - f64::from(intersect_volume)) / f64::from(totv1))
}

/// The occupancy vectors of two compatible grids and their L1 distance.
fn shapes<'a>(
    g1: &'a UniformGrid3D,
    g2: &'a UniformGrid3D,
) -> Result<(&'a DiscreteValueVect, &'a DiscreteValueVect, u32)> {
    if !g1.compare_params(g2) {
        return Err(RDKitError::Value("Grid parameters do not match".into()));
    }
    let (v1, v2) = (g1.occupancy_vect(), g2.occupancy_vect());
    let dist = compute_l1_norm_discrete(v1, v2)?;
    Ok((v1, v2, dist))
}

/// `computeGridIndices`: grid-index offsets of the points within `window_radius`.
fn grid_window_offsets(grid: &UniformGrid3D, window_radius: f64) -> Vec<i32> {
    let rad_in_grid = (window_radius / grid.spacing()).ceil() as i32;
    let (dx, dy) = (grid.num_x() as i32, grid.num_y() as i32);
    let mut res = Vec::new();
    for i in -rad_in_grid..=rad_in_grid {
        for j in -rad_in_grid..=rad_in_grid {
            for k in -rad_in_grid..=rad_in_grid {
                let r2 = f64::from(i * i + j * j + k * k);
                if r2.sqrt() as i32 <= rad_in_grid {
                    res.push((i * dy + j) * dx + k);
                }
            }
        }
    }
    res
}

/// `computeGridCentroid`: the weighted centroid of the grid points around
/// `pt` and the sum of the weights.
pub fn compute_grid_centroid(
    grid: &UniformGrid3D,
    pt: &Point3D,
    window_radius: f64,
) -> Result<(Point3D, f64)> {
    let mut weight_sum = 0.0;
    let v1 = grid.occupancy_vect();
    let mut centroid = Point3D::new(0.0, 0.0, 0.0);
    // an index of -1 wraps around, as in upstream's unsigned arithmetic
    let idx_i = grid.grid_point_index(pt) as u32;
    for offset in grid_window_offsets(grid, window_radius) {
        let idx = idx_i.wrapping_add(offset as u32) as i32;
        if idx >= 0 && (idx as u32) < v1.length() {
            let wt = v1.get_val(idx as u32)?;
            centroid += grid.grid_point_loc(idx as u32)? * f64::from(wt);
            weight_sum += f64::from(wt);
        }
    }
    Ok((centroid / weight_sum, weight_sum))
}

/// `findGridTerminalPoints` (part of the subshape algorithm).
pub fn find_grid_terminal_points(
    grid: &UniformGrid3D,
    window_radius: f64,
    inclusion_fraction: f64,
) -> Result<Vec<Point3D>> {
    let mut res = Vec::new();
    let offsets = grid_window_offsets(grid, window_radius);
    let storage = grid.occupancy_vect();
    let max_grid_val = (1u32 << (1u32 << storage.value_type())) - 1;
    for i in 0..storage.length() {
        if storage.get_val(i)? < max_grid_val {
            continue;
        }
        // the weighted volume of the shape inside the sphere
        let mut vol_in_sphere = 0.0;
        let mut n_pts_here = 0u32;
        for &offset in &offsets {
            let idx = i.wrapping_add(offset as u32) as i32;
            if idx >= 0 && (idx as u32) < storage.length() {
                vol_in_sphere += f64::from(storage.get_val(idx as u32)?);
                n_pts_here += 1;
            }
        }
        // the shape may be cut off by the edge of the grid, so the actual max
        // volume in the sphere may well be less than the theoretical max
        let max_poss_val_in_sphere = f64::from(n_pts_here) * f64::from(max_grid_val);
        if vol_in_sphere / max_poss_val_in_sphere <= inclusion_fraction {
            let pt_i = grid.grid_point_loc(i)?;
            let (centroid, _) = compute_grid_centroid(grid, &pt_i, window_radius)?;
            res.push(centroid);
        }
    }
    Ok(res)
}

const REAL_OFFSET_TOL: f64 = 1.0e-4;
const REAL_SPACING_TOL: f64 = 1.0e-4;
/// `RealValueVect::compareVectors` tolerance
const REAL_VAL_TOL: f64 = 0.01;

/// A grid of doubles.
#[derive(Clone, Debug)]
pub struct UniformRealValueGrid3D {
    geometry: GridGeometry,
    storage: RealValueVect,
}

impl UniformRealValueGrid3D {
    /// `UniformRealValueGrid3D(dimX, dimY, dimZ, spacing, offset)`; the default
    /// offset centers the grid on the origin.
    pub fn new(
        dim_x: f64,
        dim_y: f64,
        dim_z: f64,
        spacing: f64,
        offset: Option<Point3D>,
    ) -> Result<Self> {
        let offset =
            offset.unwrap_or_else(|| Point3D::new(-0.5 * dim_x, -0.5 * dim_y, -0.5 * dim_z));
        let geometry = GridGeometry::new(dim_x, dim_y, dim_z, spacing, offset)?;
        let storage = RealValueVect::new(geometry.size());
        Ok(UniformRealValueGrid3D { geometry, storage })
    }

    /// The default constructor: a unit grid.
    pub fn unit() -> Result<Self> {
        Self::new(1.0, 1.0, 1.0, 1.0, Some(Point3D::default()))
    }

    /// Reads a grid from its pickle (`toString()`).
    pub fn from_binary(pkl: &[u8]) -> Result<Self> {
        let (geometry, storage) = GridGeometry::read_header(pkl, "UniformRealValueGrid3D")?;
        Ok(UniformRealValueGrid3D {
            geometry,
            storage: RealValueVect::from_binary(storage)?,
        })
    }

    /// `toString()`
    pub fn to_binary(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.geometry.write_header(&mut out);
        let storage = self.storage.to_binary();
        out.extend_from_slice(&(storage.len() as u32).to_le_bytes());
        out.extend_from_slice(&storage);
        out
    }

    pub fn num_x(&self) -> u32 {
        self.geometry.num_x
    }
    pub fn num_y(&self) -> u32 {
        self.geometry.num_y
    }
    pub fn num_z(&self) -> u32 {
        self.geometry.num_z
    }
    pub fn offset(&self) -> Point3D {
        self.geometry.offset
    }
    pub fn spacing(&self) -> f64 {
        self.geometry.spacing
    }
    pub fn size(&self) -> u32 {
        self.geometry.size()
    }
    pub fn occupancy_vect(&self) -> &RealValueVect {
        &self.storage
    }

    pub fn grid_index(&self, xi: u32, yi: u32, zi: u32) -> i32 {
        self.geometry.grid_index(xi, yi, zi)
    }
    pub fn grid_indices(&self, idx: u32) -> Result<(u32, u32, u32)> {
        self.geometry.grid_indices(idx)
    }
    pub fn grid_point_index(&self, point: &Point3D) -> i32 {
        self.geometry.point_index(point)
    }
    pub fn grid_point_loc(&self, id: u32) -> Result<Point3D> {
        self.geometry.point_loc(id)
    }

    /// The value at the grid point closest to `point`, -1 outside the grid.
    pub fn val_at_point(&self, point: &Point3D) -> Result<f64> {
        let id = self.grid_point_index(point);
        if id < 0 {
            return Ok(-1.0);
        }
        self.storage.get_val(id as u32)
    }
    pub fn val(&self, id: u32) -> Result<f64> {
        self.storage.get_val(id)
    }
    pub fn set_val(&mut self, id: u32, val: f64) -> Result<()> {
        self.storage.set_val(id, val)
    }
    pub fn set_val_at_point(&mut self, point: &Point3D, val: f64) -> Result<()> {
        let id = self.grid_point_index(point);
        if id < 0 {
            return Ok(());
        }
        self.storage.set_val(id as u32, val)
    }

    /// `compareParams`
    pub fn compare_params(&self, other: &Self) -> bool {
        self.geometry
            .compatible_with(&other.geometry, REAL_SPACING_TOL, REAL_OFFSET_TOL)
    }

    /// `compareVectors`: entries agree to a relative 1%.
    pub fn compare_vectors(&self, other: &Self) -> Result<bool> {
        let (a, b) = (self.storage.data(), other.storage.data());
        if a.len() != b.len() {
            return Err(RDKitError::Value(
                "Comparing vectors of different lengths".into(),
            ));
        }
        Ok(a.iter()
            .zip(b)
            .all(|(&v1, &v2)| ((v1 - v2) / if v1 != 0.0 { v1 } else { 1.0 }).abs() <= REAL_VAL_TOL))
    }

    /// `compareGrids`
    pub fn compare_grids(&self, other: &Self) -> Result<bool> {
        if !self.compare_params(other) {
            return Ok(false);
        }
        self.compare_vectors(other)
    }

    fn apply(&mut self, other: &Self, f: impl Fn(f64, f64) -> f64) -> Result<()> {
        if !self.compare_params(other) {
            return Err(RDKitError::precondition("incompatible grids"));
        }
        self.storage.apply(&other.storage, f)
    }

    /// `operator|=`: element-wise maximum
    pub fn union_assign(&mut self, other: &Self) -> Result<()> {
        self.apply(other, cxx_max)
    }
    /// `operator&=`: element-wise minimum
    pub fn intersect_assign(&mut self, other: &Self) -> Result<()> {
        self.apply(other, cxx_min)
    }
    /// `operator+=`
    pub fn add_assign(&mut self, other: &Self) -> Result<()> {
        self.apply(other, |a, b| a + b)
    }
    /// `operator-=`
    pub fn sub_assign(&mut self, other: &Self) -> Result<()> {
        self.apply(other, |a, b| a - b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dimensions_and_indices() {
        let g = UniformGrid3D::new(20.0, 18.0, 15.0, 0.5, 1, None).unwrap();
        assert_eq!((g.num_x(), g.num_y(), g.num_z()), (40, 36, 30));
        let idx = g.grid_index(3, 2, 1);
        assert_eq!(g.grid_indices(idx as u32).unwrap(), (3, 2, 1));
        assert_eq!(g.grid_index(40, 0, 0), -1);
        assert!(g.grid_indices(g.size()).is_err());
        assert_eq!(g.grid_point_index(&Point3D::new(100.0, 0.0, 0.0)), -1);
    }

    #[test]
    fn pickle_round_trip() {
        let mut g = UniformGrid3D::new(10.0, 9.0, 8.0, 0.5, 1, None).unwrap();
        g.set_sphere_occupancy(&Point3D::new(1.0, 1.0, 0.0), 1.5, 0.25, -1, true)
            .unwrap();
        let g2 = UniformGrid3D::from_binary(&g.to_binary()).unwrap();
        assert_eq!(tanimoto_distance(&g, &g2).unwrap(), 0.0);
        let r = UniformRealValueGrid3D::new(4.0, 4.0, 4.0, 1.0, None).unwrap();
        let r2 = UniformRealValueGrid3D::from_binary(&r.to_binary()).unwrap();
        assert!(r.compare_grids(&r2).unwrap());
    }

    #[test]
    fn mismatched_grids_are_rejected() {
        let a = UniformGrid3D::new(10.0, 10.0, 10.0, 0.5, 1, None).unwrap();
        let b = UniformGrid3D::new(10.0, 10.0, 10.0, 1.0, 1, None).unwrap();
        assert!(tanimoto_distance(&a, &b).is_err());
    }

    #[test]
    fn cxx_double_formatting() {
        assert_eq!(cxx_f64(10.0), "10");
        assert_eq!(cxx_f64(0.5), "0.5");
        assert_eq!(cxx_f64(1.0 / 3.0), "0.333333");
        assert_eq!(cxx_f64(123456.0), "123456");
        assert_eq!(cxx_f64(1234567.0), "1.23457e+06");
        assert_eq!(cxx_f64(0.00001234), "1.234e-05");
    }
}
