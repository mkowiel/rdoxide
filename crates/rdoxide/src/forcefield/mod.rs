//! The force-field framework (upstream `Code/ForceField/ForceField.{h,cpp}`
//! and `Contrib.h`).
//!
//! A [`ForceField`] owns the flat position array (`dimension` values per
//! point) instead of pointing at external points; callers copy coordinates in
//! and out. Upstream caches pairwise distances between calls; the values are
//! the same as computing them on demand, so the cache is not kept.

pub mod constraints;
pub mod mmff;
pub mod uff;

use crate::error::{RDKitError, Result};
use crate::geometry::Point3D;
use crate::numerics::bfgs::{self, Snapshot};

/// One term of the energy (`ForceFieldContrib`).
pub trait ForceFieldContrib: Send + Sync {
    /// the contribution to the energy at `pos`
    fn energy(&self, pos: &[f64]) -> f64;
    /// adds the contribution to the gradient at `pos` into `grad`
    fn add_grad(&self, pos: &[f64], grad: &mut [f64]);
    fn box_clone(&self) -> Box<dyn ForceFieldContrib>;
    /// divides the force constant of a torsion term (`scaleForceConstant`);
    /// other terms ignore it
    fn scale_torsion_force_constant(&mut self, _count: u32) {}
}

impl Clone for Box<dyn ForceFieldContrib> {
    fn clone(&self) -> Self {
        self.box_clone()
    }
}

/// `ForceFields::ForceField`
#[derive(Clone)]
pub struct ForceField {
    dimension: usize,
    positions: Vec<f64>,
    contribs: Vec<Box<dyn ForceFieldContrib>>,
    fixed_points: Vec<usize>,
}

impl std::fmt::Debug for ForceField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ForceField")
            .field("dimension", &self.dimension)
            .field("points", &self.num_points())
            .field("contribs", &self.contribs.len())
            .field("fixed_points", &self.fixed_points)
            .finish()
    }
}

impl ForceField {
    pub fn new(dimension: usize) -> Self {
        ForceField {
            dimension,
            positions: Vec::new(),
            contribs: Vec::new(),
            fixed_points: Vec::new(),
        }
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn num_points(&self) -> usize {
        self.positions.len() / self.dimension
    }

    /// Appends a point; `coords` must have `dimension` values.
    pub fn push_point(&mut self, coords: &[f64]) {
        assert_eq!(coords.len(), self.dimension, "wrong point dimension");
        self.positions.extend_from_slice(coords);
    }

    /// The flat positions (`dimension` values per point).
    pub fn positions(&self) -> &[f64] {
        &self.positions
    }

    pub fn set_positions(&mut self, pos: &[f64]) -> Result<()> {
        if pos.len() != self.positions.len() {
            return Err(RDKitError::invariant("size mismatch"));
        }
        self.positions.copy_from_slice(pos);
        Ok(())
    }

    pub fn push_contrib(&mut self, contrib: Box<dyn ForceFieldContrib>) {
        self.contribs.push(contrib);
    }

    pub fn contribs(&self) -> &[Box<dyn ForceFieldContrib>] {
        &self.contribs
    }

    /// Divides the force constants of the torsion terms from index `first`
    /// on by `count`.
    pub fn scale_torsions_from(&mut self, first: usize, count: u32) {
        for c in &mut self.contribs[first..] {
            c.scale_torsion_force_constant(count);
        }
    }

    /// Appends a point and returns its index.
    pub fn add_point(&mut self, coords: &[f64]) -> usize {
        self.push_point(coords);
        self.num_points() - 1
    }

    /// Replaces the first `pos.len()` coordinates (the points that are not
    /// extra points).
    pub fn set_leading_positions(&mut self, pos: &[f64]) -> Result<()> {
        if pos.len() > self.positions.len() {
            return Err(RDKitError::invariant("size mismatch"));
        }
        self.positions[..pos.len()].copy_from_slice(pos);
        Ok(())
    }

    pub fn fixed_points(&self) -> &[usize] {
        &self.fixed_points
    }

    pub fn fixed_points_mut(&mut self) -> &mut Vec<usize> {
        &mut self.fixed_points
    }

    /// The energy (kcal/mol) at the stored positions; with `contribs` the
    /// energy of every term is returned as well.
    pub fn calc_energy(&self) -> f64 {
        self.energy_at(&self.positions)
    }

    pub fn calc_energy_contribs(&self) -> (f64, Vec<f64>) {
        let parts: Vec<f64> = self
            .contribs
            .iter()
            .map(|c| c.energy(&self.positions))
            .collect();
        let mut total = 0.0;
        for e in &parts {
            total += e;
        }
        (total, parts)
    }

    /// The energy at `pos` (`3 * numPoints` long).
    pub fn energy_at(&self, pos: &[f64]) -> f64 {
        let mut res = 0.0;
        for c in &self.contribs {
            res += c.energy(pos);
        }
        res
    }

    /// The gradient at the stored positions; fixed points get zero.
    pub fn calc_grad(&self) -> Vec<f64> {
        let mut grad = vec![0.0; self.positions.len()];
        self.grad_at(&self.positions, &mut grad);
        grad
    }

    /// Adds the gradient at `pos` into `grad` and zeroes the fixed points.
    pub fn grad_at(&self, pos: &[f64], grad: &mut [f64]) {
        if self.contribs.is_empty() {
            return;
        }
        for c in &self.contribs {
            c.add_grad(pos, grad);
        }
        for &p in &self.fixed_points {
            let idx = self.dimension * p;
            for g in &mut grad[idx..idx + self.dimension] {
                *g = 0.0;
            }
        }
    }

    /// Minimizes the energy, updating the positions. Returns 0 on success and
    /// 1 if the minimization did not converge in `max_its` iterations.
    pub fn minimize(&mut self, max_its: u32, force_tol: f64, energy_tol: f64) -> Result<i32> {
        Ok(self
            .minimize_with_snapshots(0, max_its, force_tol, energy_tol)?
            .0)
    }

    /// Like [`minimize`](Self::minimize), also recording a snapshot of the
    /// positions every `snapshot_freq` iterations (0: none).
    pub fn minimize_with_snapshots(
        &mut self,
        snapshot_freq: u32,
        max_its: u32,
        force_tol: f64,
        _energy_tol: f64,
    ) -> Result<(i32, Vec<Snapshot>)> {
        if self.contribs.is_empty() {
            return Ok((0, Vec::new()));
        }
        let mut points = self.positions.clone();
        let result = {
            let ff = &*self;
            bfgs::minimize(
                &mut points,
                force_tol,
                |pos| ff.energy_at(pos),
                |pos, grad| ff.scaled_grad(pos, grad),
                snapshot_freq,
                max_its,
            )?
        };
        self.positions = points;
        Ok((i32::from(!result.converged), result.snapshots))
    }

    /// The gradient the minimizer works with: scaled down (upstream's
    /// "hack" for efficiency and numeric stability). Returns the scale.
    fn scaled_grad(&self, pos: &[f64], grad: &mut [f64]) -> f64 {
        grad.fill(0.0);
        self.grad_at(pos, grad);
        let mut max_grad = -1e8f64;
        let mut grad_scale = 0.1;
        for g in grad.iter_mut() {
            *g *= grad_scale;
            if g.abs() > max_grad {
                max_grad = g.abs();
            }
        }
        if max_grad > 10.0 {
            while max_grad * grad_scale > 10.0 {
                grad_scale *= 0.5;
            }
            for g in grad.iter_mut() {
                *g *= grad_scale;
            }
        }
        grad_scale
    }

    /// The distance between points `i` and `j` at `pos` (default: stored).
    pub fn distance(&self, i: usize, j: usize, pos: Option<&[f64]>) -> Result<f64> {
        Ok(self.distance2(i, j, pos)?.sqrt())
    }

    pub fn distance2(&self, i: usize, j: usize, pos: Option<&[f64]>) -> Result<f64> {
        let n = self.num_points();
        if i >= n || j >= n {
            return Err(RDKitError::Index("Index out of bounds".into()));
        }
        let pos = pos.unwrap_or(&self.positions);
        let d = self.dimension;
        let mut res = 0.0;
        for k in 0..d {
            let t = pos[d * i + k] - pos[d * j + k];
            res += t * t;
        }
        Ok(res)
    }
}

/// `ForceFieldsHelper::normalizeAngleDeg`: maps an angle to [-180, 180]
pub fn normalize_angle_deg(angle: f64) -> f64 {
    let mut a = angle % 360.0;
    if a < -180.0 {
        a += 360.0;
    } else if a > 180.0 {
        a -= 360.0;
    }
    a
}

/// The point of `pos` with index `idx` (3 dimensions).
pub fn point3(pos: &[f64], idx: usize) -> Point3D {
    Point3D {
        x: pos[3 * idx],
        y: pos[3 * idx + 1],
        z: pos[3 * idx + 2],
    }
}

/// Intermediate vectors of a dihedral computation.
#[derive(Debug, Clone)]
pub struct DihedralParts {
    pub r: [Point3D; 4],
    pub t: [Point3D; 2],
    pub d: [f64; 2],
    pub cos_phi: f64,
    pub dihedral: f64,
}

/// `ForceFieldsHelper::computeDihedral`
pub fn compute_dihedral(p1: &Point3D, p2: &Point3D, p3: &Point3D, p4: &Point3D) -> DihedralParts {
    let r0 = *p1 - *p2;
    let r1 = *p3 - *p2;
    let r2 = -r1;
    let r3 = *p4 - *p3;
    let mut t0 = r0.cross_product(&r1);
    let d0 = t0.length().max(1.0e-5);
    t0 /= d0;
    let mut t1 = r2.cross_product(&r3);
    let d1 = t1.length().max(1.0e-5);
    t1 /= d1;
    let cos_phi = t0.dot_product(&t1).clamp(-1.0, 1.0);
    // a signed dihedral, hence atan2
    let m = t0.cross_product(&r1);
    let m_length = m.length().max(1.0e-5);
    let dihedral = -(m.dot_product(&t1) / m_length).atan2(cos_phi);
    DihedralParts {
        r: [r0, r1, r2, r3],
        t: [t0, t1],
        d: [d0, d1],
        cos_phi,
        dihedral,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn angles_are_normalized() {
        assert_eq!(normalize_angle_deg(190.0), -170.0);
        assert_eq!(normalize_angle_deg(-190.0), 170.0);
        assert_eq!(normalize_angle_deg(540.0), 180.0);
    }

    #[test]
    fn dihedral_is_signed() {
        let a = Point3D {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };
        let b = Point3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let c = Point3D {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        };
        let d = Point3D {
            x: 0.0,
            y: 1.0,
            z: 1.0,
        };
        let phi = compute_dihedral(&a, &b, &c, &d).dihedral;
        assert!((phi.abs() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }
}
