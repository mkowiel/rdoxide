//! Distance, angle, torsion and position constraints (upstream
//! `ForceField/{Distance,Angle,Torsion,Position}Constraint.{h,cpp}`); shared
//! by the UFF and MMFF force fields.

use super::{ForceField, ForceFieldContrib, compute_dihedral, normalize_angle_deg, point3};
use crate::error::{RDKitError, Result};
use crate::geometry::Point3D;

const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

fn check_idx(field: &ForceField, idx: usize, file: &str, line: u32) -> Result<()> {
    if idx >= field.num_points() {
        return Err(RDKitError::range_at(
            "idx",
            idx,
            field.num_points(),
            file,
            line,
        ));
    }
    Ok(())
}

fn range_check(lo: f64, x: f64, hi: f64, name: &str, file: &str, line: u32) -> Result<()> {
    if x < lo || x > hi {
        return Err(RDKitError::Invariant {
            kind: "Range Error",
            message: format!(
                "{name}\n\tViolation occurred on line {line} in file {file}\n\tFailed Expression: {lo} <= {x} <= {hi}"
            ),
        });
    }
    Ok(())
}

// ----------------------------------------------------------------- distance

/// The distance over all the dimensions of the points.
fn distance(pos: &[f64], dim: usize, i: usize, j: usize) -> f64 {
    (0..dim)
        .map(|k| {
            let d = pos[dim * i + k] - pos[dim * j + k];
            d * d
        })
        .sum::<f64>()
        .sqrt()
}

#[derive(Debug, Clone)]
pub struct DistanceConstraintContrib {
    dim: usize,
    end1: usize,
    end2: usize,
    min_len: f64,
    max_len: f64,
    force_constant: f64,
}

impl DistanceConstraintContrib {
    /// With `relative` the bounds are relative to the current distance.
    pub fn new(
        field: &ForceField,
        idx1: usize,
        idx2: usize,
        relative: bool,
        mut min_len: f64,
        mut max_len: f64,
        force_const: f64,
    ) -> Result<Self> {
        let file = "Code/ForceField/DistanceConstraint.cpp";
        check_idx(field, idx1, file, 37)?;
        check_idx(field, idx2, file, 38)?;
        if max_len < min_len {
            return Err(RDKitError::precondition("bad bounds"));
        }
        let dim = field.dimension();
        if relative {
            let dist = distance(field.positions(), dim, idx1, idx2);
            min_len = (dist + min_len).max(0.0);
            max_len = (dist + max_len).max(0.0);
        }
        Ok(DistanceConstraintContrib {
            dim,
            end1: idx1,
            end2: idx2,
            min_len,
            max_len,
            force_constant: force_const,
        })
    }
}

impl ForceFieldContrib for DistanceConstraintContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let dist = distance(pos, self.dim, self.end1, self.end2);
        let mut term = 0.0;
        if dist < self.min_len {
            term = self.min_len - dist;
        } else if dist > self.max_len {
            term = dist - self.max_len;
        }
        0.5 * self.force_constant * term * term
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let dist = distance(pos, self.dim, self.end1, self.end2);
        let mut pre_factor = if dist < self.min_len {
            dist - self.min_len
        } else if dist > self.max_len {
            dist - self.max_len
        } else {
            return;
        };
        pre_factor *= self.force_constant;
        for i in 0..self.dim {
            let d_grad = pre_factor
                * (pos[self.dim * self.end1 + i] - pos[self.dim * self.end2 + i])
                / dist.max(1.0e-8);
            grad[self.dim * self.end1 + i] += d_grad;
            grad[self.dim * self.end2 + i] -= d_grad;
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// -------------------------------------------------------------------- angle

#[derive(Debug, Clone)]
pub struct AngleConstraintContrib {
    dim: usize,
    at1: usize,
    at2: usize,
    at3: usize,
    min_angle_deg: f64,
    max_angle_deg: f64,
    force_constant: f64,
}

/// the angle (degrees) at the middle point
fn angle_deg(
    pos: &[f64],
    dim: usize,
    at1: usize,
    at2: usize,
    at3: usize,
) -> (f64, [Point3D; 2], [f64; 2]) {
    let point = |i: usize| Point3D::new(pos[dim * i], pos[dim * i + 1], pos[dim * i + 2]);
    let p2 = point(at2);
    let r = [point(at1) - p2, point(at3) - p2];
    let len_sq = [r[0].length_sq().max(1.0e-5), r[1].length_sq().max(1.0e-5)];
    let cos_theta = (r[0].dot_product(&r[1]) / (len_sq[0] * len_sq[1]).sqrt()).clamp(-1.0, 1.0);
    (RAD2DEG * cos_theta.acos(), r, len_sq)
}

impl AngleConstraintContrib {
    /// With `relative` the bounds are relative to the current angle.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        relative: bool,
        mut min_angle_deg: f64,
        mut max_angle_deg: f64,
        force_const: f64,
    ) -> Result<Self> {
        let file = "Code/ForceField/AngleConstraint.cpp";
        check_idx(field, idx1, file, 69)?;
        check_idx(field, idx2, file, 70)?;
        check_idx(field, idx3, file, 71)?;
        if min_angle_deg > max_angle_deg {
            return Err(RDKitError::precondition(
                "minAngleDeg must be <= maxAngleDeg",
            ));
        }
        if relative {
            let (angle, ..) = angle_deg(field.positions(), field.dimension(), idx1, idx2, idx3);
            min_angle_deg += angle;
            max_angle_deg += angle;
        }
        range_check(0.0, min_angle_deg, 180.0, "minAngleDeg", file, 87)?;
        range_check(0.0, max_angle_deg, 180.0, "maxAngleDeg", file, 88)?;
        Ok(AngleConstraintContrib {
            dim: field.dimension(),
            at1: idx1,
            at2: idx2,
            at3: idx3,
            min_angle_deg,
            max_angle_deg,
            force_constant: force_const,
        })
    }

    fn angle_term(&self, angle: f64) -> f64 {
        if angle < self.min_angle_deg {
            angle - self.min_angle_deg
        } else if angle > self.max_angle_deg {
            angle - self.max_angle_deg
        } else {
            0.0
        }
    }
}

impl ForceFieldContrib for AngleConstraintContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let (angle, ..) = angle_deg(pos, self.dim, self.at1, self.at2, self.at3);
        let term = self.angle_term(angle);
        self.force_constant * term * term
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let (angle, r, len_sq) = angle_deg(pos, self.dim, self.at1, self.at2, self.at3);
        let term = self.angle_term(angle);
        let de_dtheta = 2.0 * RAD2DEG * self.force_constant * term;
        let rp = r[1].cross_product(&r[0]);
        let prefactor = de_dtheta / rp.length().max(1.0e-5);
        let t = [-prefactor / len_sq[0], prefactor / len_sq[1]];
        let dedp0 = r[0].cross_product(&rp) * t[0];
        let dedp2 = r[1].cross_product(&rp) * t[1];
        let dedp1 = -dedp0 - dedp2;
        for (idx, d) in [(self.at1, dedp0), (self.at2, dedp1), (self.at3, dedp2)] {
            grad[self.dim * idx] += d.x;
            grad[self.dim * idx + 1] += d.y;
            grad[self.dim * idx + 2] += d.z;
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ------------------------------------------------------------------ torsion

#[derive(Debug, Clone)]
pub struct TorsionConstraintContrib {
    at: [usize; 4],
    min_dihedral_deg: f64,
    max_dihedral_deg: f64,
    force_constant: f64,
}

impl TorsionConstraintContrib {
    /// With `relative` the bounds are relative to the current dihedral.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: &ForceField,
        idx: [usize; 4],
        relative: bool,
        mut min_dihedral_deg: f64,
        mut max_dihedral_deg: f64,
        force_const: f64,
    ) -> Result<Self> {
        if min_dihedral_deg > max_dihedral_deg {
            return Err(RDKitError::precondition(
                "minDihedralDeg must be <= maxDihedralDeg",
            ));
        }
        let file = "Code/ForceField/TorsionConstraint.cpp";
        for (k, &i) in idx.iter().enumerate() {
            check_idx(field, i, file, 28 + k as u32)?;
        }
        if relative {
            let p = field.positions();
            let dihedral = compute_dihedral(
                &point3(p, idx[0]),
                &point3(p, idx[1]),
                &point3(p, idx[2]),
                &point3(p, idx[3]),
            )
            .dihedral
                * RAD2DEG;
            min_dihedral_deg += dihedral;
            max_dihedral_deg += dihedral;
        }
        Ok(TorsionConstraintContrib {
            at: idx,
            min_dihedral_deg: normalize_angle_deg(min_dihedral_deg),
            max_dihedral_deg: normalize_angle_deg(max_dihedral_deg),
            force_constant: force_const,
        })
    }

    fn dihedral_term(&self, dihedral: f64) -> f64 {
        let (lo, hi) = (self.min_dihedral_deg, self.max_dihedral_deg);
        let mut target = dihedral;
        if !(dihedral > lo && dihedral < hi)
            && !(dihedral > lo && lo > hi)
            && !(dihedral < hi && lo > hi)
        {
            let min_target = normalize_angle_deg(dihedral - lo);
            let max_target = normalize_angle_deg(dihedral - hi);
            target = if min_target.abs() < max_target.abs() {
                lo
            } else {
                hi
            };
        }
        normalize_angle_deg(dihedral - target)
    }
}

impl ForceFieldContrib for TorsionConstraintContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let [a, b, c, d] = self.at;
        let dihedral = compute_dihedral(
            &point3(pos, a),
            &point3(pos, b),
            &point3(pos, c),
            &point3(pos, d),
        )
        .dihedral
            * RAD2DEG;
        let term = self.dihedral_term(dihedral);
        self.force_constant * term * term
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let [a, b, c, d] = self.at;
        let parts = compute_dihedral(
            &point3(pos, a),
            &point3(pos, b),
            &point3(pos, c),
            &point3(pos, d),
        );
        let r = parts.r;
        let dihedral = parts.dihedral * RAD2DEG;
        let term = self.dihedral_term(dihedral);
        let de_dphi = 2.0 * RAD2DEG * self.force_constant * term;
        let d23 = (point3(pos, b) - point3(pos, c)).length();
        let r31 = point3(pos, c) - point3(pos, a);
        let r42 = point3(pos, d) - point3(pos, b);
        let prefactor = de_dphi / d23;
        let tt = [r[0].cross_product(&r[1]), r[2].cross_product(&r[3])];
        let dedt = [
            tt[0].cross_product(&r[2]) / tt[0].length_sq() * prefactor,
            tt[1].cross_product(&r[1]) / tt[1].length_sq() * prefactor,
        ];
        let dedp = [
            r[2].cross_product(&dedt[0]),
            r31.cross_product(&dedt[0]) - r[3].cross_product(&dedt[1]),
            r[0].cross_product(&dedt[0]) + r42.cross_product(&dedt[1]),
            r[2].cross_product(&dedt[1]),
        ];
        for (idx, g) in self.at.iter().zip(dedp) {
            grad[3 * idx] += g.x;
            grad[3 * idx + 1] += g.y;
            grad[3 * idx + 2] += g.z;
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ----------------------------------------------------------------- position

#[derive(Debug, Clone)]
pub struct PositionConstraintContrib {
    at: usize,
    max_displ: f64,
    pos0: Point3D,
    force_constant: f64,
}

impl PositionConstraintContrib {
    pub fn new(field: &ForceField, idx: usize, max_displ: f64, force_const: f64) -> Result<Self> {
        check_idx(field, idx, "Code/ForceField/PositionConstraint.cpp", 22)?;
        Ok(PositionConstraintContrib {
            at: idx,
            max_displ,
            pos0: point3(field.positions(), idx),
            force_constant: force_const,
        })
    }
}

impl ForceFieldContrib for PositionConstraintContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let dist = (point3(pos, self.at) - self.pos0).length();
        let term = (dist - self.max_displ).max(0.0);
        0.5 * self.force_constant * term * term
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let p = point3(pos, self.at);
        let dist = (p - self.pos0).length();
        if dist <= self.max_displ {
            return;
        }
        let pre_factor = (dist - self.max_displ) * self.force_constant;
        let (dp, d0) = ([p.x, p.y, p.z], [self.pos0.x, self.pos0.y, self.pos0.z]);
        for i in 0..3 {
            grad[3 * self.at + i] += pre_factor * (dp[i] - d0[i]) / dist.max(1.0e-8);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn field() -> ForceField {
        let mut ff = ForceField::new(3);
        for p in [
            [0.1, 0.2, 0.0],
            [1.5, 0.3, 0.2],
            [2.1, 1.6, -0.1],
            [3.4, 1.9, 0.9],
        ] {
            ff.push_point(&p);
        }
        ff
    }

    fn check_gradient(c: &dyn ForceFieldContrib, pos: &[f64]) {
        let mut grad = vec![0.0; pos.len()];
        c.add_grad(pos, &mut grad);
        for i in 0..pos.len() {
            let h = 1e-6;
            let (mut plus, mut minus) = (pos.to_vec(), pos.to_vec());
            plus[i] += h;
            minus[i] -= h;
            let numeric = (c.energy(&plus) - c.energy(&minus)) / (2.0 * h);
            assert!(
                (numeric - grad[i]).abs() < 1e-4 * (1.0 + numeric.abs()),
                "{i}: {} {numeric}",
                grad[i]
            );
        }
    }

    #[test]
    fn gradients_match_finite_differences() {
        let ff = field();
        let pos = ff.positions().to_vec();
        let d = DistanceConstraintContrib::new(&ff, 0, 3, false, 0.5, 1.0, 100.0).unwrap();
        check_gradient(&d, &pos);
        let a = AngleConstraintContrib::new(&ff, 0, 1, 2, false, 150.0, 170.0, 10.0).unwrap();
        check_gradient(&a, &pos);
        let t = TorsionConstraintContrib::new(&ff, [0, 1, 2, 3], false, 10.0, 20.0, 5.0).unwrap();
        check_gradient(&t, &pos);
        let mut moved = pos.clone();
        moved[0] += 0.7;
        let p = PositionConstraintContrib::new(&ff, 0, 0.2, 50.0).unwrap();
        check_gradient(&p, &moved);
    }

    #[test]
    fn relative_bounds_shift_with_the_current_value() {
        let ff = field();
        let d = DistanceConstraintContrib::new(&ff, 0, 1, true, -0.0, 0.0, 10.0).unwrap();
        assert!(d.energy(ff.positions()).abs() < 1e-12);
    }

    #[test]
    fn bad_bounds_are_rejected() {
        let ff = field();
        assert!(DistanceConstraintContrib::new(&ff, 0, 1, false, 2.0, 1.0, 1.0).is_err());
        assert!(AngleConstraintContrib::new(&ff, 0, 1, 2, false, 10.0, 200.0, 1.0).is_err());
        assert!(PositionConstraintContrib::new(&ff, 9, 1.0, 1.0).is_err());
    }
}
