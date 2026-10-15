//! The UFF energy terms (upstream `ForceField/UFF/*`): bond stretch, angle
//! bend, torsion, inversion (out of plane bend) and van der Waals.

mod params;

use std::f64::consts::PI;

pub use params::{AMIDE_BOND_ORDER, AtomicParams, DEG2RAD, G, LAMBDA, ParamCollection, RAD2DEG};

use super::{ForceFieldContrib, compute_dihedral, point3};
use crate::geometry::Point3D;
use params::{clip_to_one, is_double_zero};

const ANGLE_CORRECTION_THRESHOLD: f64 = 0.8660;

fn dist(pos: &[f64], i: usize, j: usize) -> f64 {
    let mut s = 0.0;
    for k in 0..3 {
        let t = pos[3 * i + k] - pos[3 * j + k];
        s += t * t;
    }
    s.sqrt()
}

fn add3(grad: &mut [f64], idx: usize, v: [f64; 3], sign: f64) {
    for k in 0..3 {
        grad[3 * idx + k] += sign * v[k];
    }
}

// ------------------------------------------------------------- bond stretch

/// The rest length of a bond from the Pauling and O'Keefe-Breese corrections.
pub fn calc_bond_rest_length(bond_order: f64, end1: &AtomicParams, end2: &AtomicParams) -> f64 {
    let (ri, rj) = (end1.r1, end2.r1);
    // Pauling correction
    let r_bo = -LAMBDA * (ri + rj) * bond_order.ln();
    // O'Keefe and Breese electronegativity correction
    let (xi, xj) = (end1.gmp_xi, end2.gmp_xi);
    let r_en = ri * rj * (xi.sqrt() - xj.sqrt()) * (xi.sqrt() - xj.sqrt()) / (xi * ri + xj * rj);
    ri + rj + r_bo - r_en
}

pub fn calc_bond_force_constant(rest_length: f64, end1: &AtomicParams, end2: &AtomicParams) -> f64 {
    2.0 * G * end1.z1 * end2.z1 / (rest_length * rest_length * rest_length)
}

#[derive(Debug, Clone)]
pub struct BondStretchContrib {
    end1: usize,
    end2: usize,
    rest_len: f64,
    force_constant: f64,
}

impl BondStretchContrib {
    pub fn new(
        idx1: usize,
        idx2: usize,
        bond_order: f64,
        end1: &AtomicParams,
        end2: &AtomicParams,
    ) -> Self {
        let rest_len = calc_bond_rest_length(bond_order, end1, end2);
        BondStretchContrib {
            end1: idx1,
            end2: idx2,
            rest_len,
            force_constant: calc_bond_force_constant(rest_len, end1, end2),
        }
    }
}

impl ForceFieldContrib for BondStretchContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let dist_term = dist(pos, self.end1, self.end2) - self.rest_len;
        0.5 * self.force_constant * dist_term * dist_term
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let d = dist(pos, self.end1, self.end2);
        let pre_factor = self.force_constant * (d - self.rest_len);
        for i in 0..3 {
            let d_grad = if d > 0.0 {
                pre_factor * (pos[3 * self.end1 + i] - pos[3 * self.end2 + i]) / d
            } else {
                // move a small amount in an arbitrary direction
                self.force_constant * 0.01
            };
            grad[3 * self.end1 + i] += d_grad;
            grad[3 * self.end2 + i] -= d_grad;
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// --------------------------------------------------------------- angle bend

pub fn calc_angle_force_constant(
    theta0: f64,
    bond_order12: f64,
    bond_order23: f64,
    at1: &AtomicParams,
    at2: &AtomicParams,
    at3: &AtomicParams,
) -> f64 {
    let cos_theta0 = theta0.cos();
    let r12 = calc_bond_rest_length(bond_order12, at1, at2);
    let r23 = calc_bond_rest_length(bond_order23, at2, at3);
    let r13 = (r12 * r12 + r23 * r23 - 2. * r12 * r23 * cos_theta0).sqrt();
    let beta = 2. * G / (r12 * r23);
    let pre_factor = beta * at1.z1 * at3.z1 / r13.powi(5);
    let r_term = r12 * r23;
    let inner_bit = 3. * r_term * (1. - cos_theta0 * cos_theta0) - r13 * r13 * cos_theta0;
    pre_factor * r_term * inner_bit
}

#[derive(Debug, Clone)]
pub struct AngleBendContrib {
    at1: usize,
    at2: usize,
    at3: usize,
    order: u32,
    force_constant: f64,
    theta0: f64,
    c0: f64,
    c1: f64,
    c2: f64,
}

impl AngleBendContrib {
    /// `order` is the periodicity of the cosine term (0: the general
    /// expansion); orders of 30 and above are a hack that fixes the angle
    /// for 3- and 4-membered rings with sp2 atoms.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        idx1: usize,
        idx2: usize,
        idx3: usize,
        bond_order12: f64,
        bond_order23: f64,
        at1: &AtomicParams,
        at2: &AtomicParams,
        at3: &AtomicParams,
        order: u32,
    ) -> Self {
        assert!(
            idx1 != idx2 && idx2 != idx3 && idx1 != idx3,
            "degenerate points"
        );
        let mut theta0 = at2.theta0;
        let mut order = order;
        if order >= 30 {
            match order {
                30 => theta0 = 150.0 / 180.0 * PI,
                35 => theta0 = 60.0 / 180.0 * PI,
                40 => theta0 = 135.0 / 180.0 * PI,
                45 => theta0 = 90.0 / 180.0 * PI,
                _ => {}
            }
            order = 0;
        }
        let force_constant =
            calc_angle_force_constant(theta0, bond_order12, bond_order23, at1, at2, at3);
        let (mut c0, mut c1, mut c2) = (0.0, 0.0, 0.0);
        if order == 0 {
            let sin_theta0 = theta0.sin();
            let cos_theta0 = theta0.cos();
            c2 = 1. / (4. * (sin_theta0 * sin_theta0).max(1e-8));
            c1 = -4. * c2 * cos_theta0;
            c0 = c2 * (2. * cos_theta0 * cos_theta0 + 1.);
        }
        AngleBendContrib {
            at1: idx1,
            at2: idx2,
            at3: idx3,
            order,
            force_constant,
            theta0,
            c0,
            c1,
            c2,
        }
    }

    fn energy_term(&self, cos_theta: f64, sin_theta_sq: f64) -> f64 {
        let cos2_theta = cos_theta * cos_theta - sin_theta_sq;
        if self.order == 0 {
            return self.c0 + self.c1 * cos_theta + self.c2 * cos2_theta;
        }
        let res = match self.order {
            1 => -cos_theta,
            2 => cos2_theta,
            3 => cos_theta * (cos_theta * cos_theta - 3. * sin_theta_sq),
            4 => {
                cos_theta.powi(4) - 6. * cos_theta * cos_theta * sin_theta_sq
                    + sin_theta_sq * sin_theta_sq
            }
            _ => panic!("bad order"),
        };
        (1. - res) / f64::from(self.order * self.order)
    }

    fn theta_deriv(&self, cos_theta: f64, sin_theta: f64) -> f64 {
        let sin2_theta = 2. * sin_theta * cos_theta;
        if self.order == 0 {
            return -self.force_constant * (self.c1 * sin_theta + 2. * self.c2 * sin2_theta);
        }
        let d = match self.order {
            1 => -sin_theta,
            2 => sin2_theta,
            3 => sin_theta * (3. - 4. * sin_theta * sin_theta),
            4 => cos_theta * sin_theta * (4. - 8. * sin_theta * sin_theta),
            _ => panic!("bad order"),
        };
        d * (self.force_constant / f64::from(self.order))
    }
}

impl ForceFieldContrib for AngleBendContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let dist1 = dist(pos, self.at1, self.at2);
        let dist2 = dist(pos, self.at2, self.at3);
        let p12 = point3(pos, self.at1) - point3(pos, self.at2);
        let p32 = point3(pos, self.at3) - point3(pos, self.at2);
        let cos_theta = clip_to_one(p12.dot_product(&p32) / (dist1 * dist2));
        // sin^2(theta) gives cos(2*theta)
        let sin_theta_sq = 1. - cos_theta * cos_theta;
        let mut res = self.force_constant * self.energy_term(cos_theta, sin_theta_sq);
        // UFF has no penalty for zero-degree angles, which can overlap 1-3
        // atoms (Github #7901); an empirical penalty borrowed from OpenBabel
        // is added below about 30 degrees
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            let theta = cos_theta.acos();
            res += (-20.0 * (theta - self.theta0 + 0.25)).exp();
        }
        res
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let d = [dist(pos, self.at1, self.at2), dist(pos, self.at2, self.at3)];
        let r = [
            (point3(pos, self.at1) - point3(pos, self.at2)) / d[0],
            (point3(pos, self.at3) - point3(pos, self.at2)) / d[1],
        ];
        let cos_theta = clip_to_one(r[0].dot_product(&r[1]));
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = sin_theta_sq.sqrt().max(1.0e-8);
        // chain rule: dE/dx = dE/dTheta * dTheta/dx
        let mut de_dtheta = self.theta_deriv(cos_theta, sin_theta);
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            let theta = cos_theta.acos();
            de_dtheta += -20.0 * (-20.0 * (theta - self.theta0 + 0.25)).exp();
        }
        let dcos = [
            1.0 / d[0] * (r[1].x - cos_theta * r[0].x),
            1.0 / d[0] * (r[1].y - cos_theta * r[0].y),
            1.0 / d[0] * (r[1].z - cos_theta * r[0].z),
            1.0 / d[1] * (r[0].x - cos_theta * r[1].x),
            1.0 / d[1] * (r[0].y - cos_theta * r[1].y),
            1.0 / d[1] * (r[0].z - cos_theta * r[1].z),
        ];
        for k in 0..3 {
            grad[3 * self.at1 + k] += de_dtheta * dcos[k] / (-sin_theta);
            grad[3 * self.at2 + k] += de_dtheta * (-dcos[k] - dcos[3 + k]) / (-sin_theta);
            grad[3 * self.at3 + k] += de_dtheta * dcos[3 + k] / (-sin_theta);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ------------------------------------------------------------------ torsion

/// The cosine of the torsion p1-p2-p3-p4 (0 when degenerate).
pub fn calculate_cos_torsion(p1: &Point3D, p2: &Point3D, p3: &Point3D, p4: &Point3D) -> f64 {
    let r1 = *p1 - *p2;
    let r2 = *p3 - *p2;
    let r3 = *p2 - *p3;
    let r4 = *p4 - *p3;
    let t1 = r1.cross_product(&r2);
    let t2 = r3.cross_product(&r4);
    let d1 = t1.length();
    let d2 = t2.length();
    if is_double_zero(d1) || is_double_zero(d2) {
        return 0.0;
    }
    clip_to_one(t1.dot_product(&t2) / (d1 * d2))
}

pub fn is_in_group6(num: u32) -> bool {
    matches!(num, 8 | 16 | 34 | 52 | 84)
}

/// Equation 17 of the UFF paper: the barrier of sp2-sp2 and some sp2-sp3 bonds.
pub fn equation17(bond_order23: f64, at2: &AtomicParams, at3: &AtomicParams) -> f64 {
    5. * (at2.u1 * at3.u1).sqrt() * (1. + 4.18 * bond_order23.ln())
}

/// The hybridizations the torsion parameters distinguish.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TorsionHybridization {
    Sp2,
    Sp3,
}

#[derive(Debug, Clone)]
pub struct TorsionAngleContrib {
    at1: usize,
    at2: usize,
    at3: usize,
    at4: usize,
    order: u32,
    force_constant: f64,
    cos_term: f64,
}

impl TorsionAngleContrib {
    /// The torsion atom1 - atom2 - atom3 - atom4 around the bond 2-3.
    /// `end_atom_is_sp2` signals that atom 1 or 4 is sp2 hybridized, which
    /// triggers a special case for sp3-sp2 bonds.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        idx: [usize; 4],
        bond_order23: f64,
        at_num2: u32,
        at_num3: u32,
        hyb2: TorsionHybridization,
        hyb3: TorsionHybridization,
        at2: &AtomicParams,
        at3: &AtomicParams,
        end_atom_is_sp2: bool,
    ) -> Self {
        use TorsionHybridization::{Sp2, Sp3};
        assert!(
            idx[0] != idx[1]
                && idx[0] != idx[2]
                && idx[0] != idx[3]
                && idx[1] != idx[2]
                && idx[1] != idx[3]
                && idx[2] != idx[3],
            "degenerate points"
        );
        let (force_constant, order, cos_term);
        if hyb2 == Sp3 && hyb3 == Sp3 {
            // general case
            let mut fc = (at2.v1 * at3.v1).sqrt();
            let mut ord = 3;
            // single bonds between group 6 elements
            if bond_order23 == 1.0 && is_in_group6(at_num2) && is_in_group6(at_num3) {
                let (mut v2, mut v3) = (6.8f64, 6.8f64);
                if at_num2 == 8 {
                    v2 = 2.0;
                }
                if at_num3 == 8 {
                    v3 = 2.0;
                }
                fc = (v2 * v3).sqrt();
                ord = 2;
            }
            force_constant = fc;
            order = ord;
            cos_term = -1.0;
        } else if hyb2 == Sp2 && hyb3 == Sp2 {
            force_constant = equation17(bond_order23, at2, at3);
            order = 2;
            cos_term = 1.0;
        } else {
            // sp2 - sp3, by default independent of the atom type in UFF
            let (mut fc, mut ord, mut ct) = (1.0, 6, 1.0);
            if bond_order23 == 1.0 {
                // group 6 sp3 with non-group 6 sp2
                if (hyb2 == Sp3 && is_in_group6(at_num2) && !is_in_group6(at_num3))
                    || (hyb3 == Sp3 && is_in_group6(at_num3) && !is_in_group6(at_num2))
                {
                    fc = equation17(bond_order23, at2, at3);
                    ord = 2;
                    ct = -1.0;
                } else if end_atom_is_sp2 {
                    // sp3 - sp2 - sp2 (the sp2 has another sp2 neighbor)
                    fc = 2.0;
                    ord = 3;
                    ct = -1.0;
                }
            }
            force_constant = fc;
            order = ord;
            cos_term = ct;
        }
        TorsionAngleContrib {
            at1: idx[0],
            at2: idx[1],
            at3: idx[2],
            at4: idx[3],
            order,
            force_constant,
            cos_term,
        }
    }

    /// Divides the force constant by `count` (torsions around a bond that is
    /// shared by several of them).
    pub fn scale_force_constant(&mut self, count: u32) {
        self.force_constant /= f64::from(count);
    }

    fn theta_deriv(&self, cos_theta: f64, sin_theta: f64) -> f64 {
        let sin_theta_sq = sin_theta * sin_theta;
        let res = match self.order {
            2 => 2. * sin_theta * cos_theta,
            3 => sin_theta * (3. - 4. * sin_theta_sq),
            6 => cos_theta * sin_theta * (32. * sin_theta_sq * (sin_theta_sq - 1.) + 6.),
            _ => panic!("bad order"),
        };
        res * (-(self.force_constant / 2.0 * self.cos_term) * f64::from(self.order))
    }
}

impl ForceFieldContrib for TorsionAngleContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let cos_phi = calculate_cos_torsion(
            &point3(pos, self.at1),
            &point3(pos, self.at2),
            &point3(pos, self.at3),
            &point3(pos, self.at4),
        );
        let sin_phi_sq = 1. - cos_phi * cos_phi;
        // E(phi) = V/2 * (1 - cos(n*phi_0)*cos(n*phi))
        let cos_n_phi = match self.order {
            2 => 1. - 2. * sin_phi_sq,
            3 => cos_phi * (cos_phi * cos_phi - 3. * sin_phi_sq),
            6 => 1. + sin_phi_sq * (-32. * sin_phi_sq * sin_phi_sq + 48. * sin_phi_sq - 18.),
            _ => panic!("bad order"),
        };
        self.force_constant / 2.0 * (1. - self.cos_term * cos_n_phi)
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let parts = compute_dihedral(
            &point3(pos, self.at1),
            &point3(pos, self.at2),
            &point3(pos, self.at3),
            &point3(pos, self.at4),
        );
        let (r, t, d, cos_phi) = (parts.r, parts.t, parts.d, parts.cos_phi);
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        let sin_phi = if sin_phi_sq > 0.0 {
            sin_phi_sq.sqrt()
        } else {
            0.0
        };
        // dE/dPhi is independent of cartesians
        let de_dphi = self.theta_deriv(cos_phi, sin_phi);
        let sin_term = de_dphi
            * if is_double_zero(sin_phi) {
                1.0 / cos_phi
            } else {
                1.0 / sin_phi
            };
        let dcos = [
            1.0 / d[0] * (t[1].x - cos_phi * t[0].x),
            1.0 / d[0] * (t[1].y - cos_phi * t[0].y),
            1.0 / d[0] * (t[1].z - cos_phi * t[0].z),
            1.0 / d[1] * (t[0].x - cos_phi * t[1].x),
            1.0 / d[1] * (t[0].y - cos_phi * t[1].y),
            1.0 / d[1] * (t[0].z - cos_phi * t[1].z),
        ];
        let g1 = [
            sin_term * (dcos[2] * r[1].y - dcos[1] * r[1].z),
            sin_term * (dcos[0] * r[1].z - dcos[2] * r[1].x),
            sin_term * (dcos[1] * r[1].x - dcos[0] * r[1].y),
        ];
        let g2 = [
            sin_term
                * (dcos[1] * (r[1].z - r[0].z)
                    + dcos[2] * (r[0].y - r[1].y)
                    + dcos[4] * (-r[3].z)
                    + dcos[5] * r[3].y),
            sin_term
                * (dcos[0] * (r[0].z - r[1].z)
                    + dcos[2] * (r[1].x - r[0].x)
                    + dcos[3] * r[3].z
                    + dcos[5] * (-r[3].x)),
            sin_term
                * (dcos[0] * (r[1].y - r[0].y)
                    + dcos[1] * (r[0].x - r[1].x)
                    + dcos[3] * (-r[3].y)
                    + dcos[4] * r[3].x),
        ];
        let g3 = [
            sin_term
                * (dcos[1] * r[0].z
                    + dcos[2] * (-r[0].y)
                    + dcos[4] * (r[3].z - r[2].z)
                    + dcos[5] * (r[2].y - r[3].y)),
            sin_term
                * (dcos[0] * (-r[0].z)
                    + dcos[2] * r[0].x
                    + dcos[3] * (r[2].z - r[3].z)
                    + dcos[5] * (r[3].x - r[2].x)),
            sin_term
                * (dcos[0] * r[0].y
                    + dcos[1] * (-r[0].x)
                    + dcos[3] * (r[3].y - r[2].y)
                    + dcos[4] * (r[2].x - r[3].x)),
        ];
        let g4 = [
            sin_term * (dcos[4] * r[2].z - dcos[5] * r[2].y),
            sin_term * (dcos[5] * r[2].x - dcos[3] * r[2].z),
            sin_term * (dcos[3] * r[2].y - dcos[4] * r[2].x),
        ];
        add3(grad, self.at1, g1, 1.0);
        add3(grad, self.at2, g2, 1.0);
        add3(grad, self.at3, g3, 1.0);
        add3(grad, self.at4, g4, 1.0);
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ---------------------------------------------------------------- inversion

/// The cosine of the angle between the plane (i, j, k) normal and the bond
/// j-l (0 when degenerate).
pub fn calculate_cos_y(i: &Point3D, j: &Point3D, k: &Point3D, l: &Point3D) -> f64 {
    const ZERO_TOL: f64 = 1.0e-16;
    let r_ji = *i - *j;
    let r_jk = *k - *j;
    let r_jl = *l - *j;
    let (l2_ji, l2_jk, l2_jl) = (r_ji.length_sq(), r_jk.length_sq(), r_jl.length_sq());
    if l2_ji < ZERO_TOL || l2_jk < ZERO_TOL || l2_jl < ZERO_TOL {
        return 0.0;
    }
    let mut n = r_ji.cross_product(&r_jk);
    n /= l2_ji.sqrt() * l2_jk.sqrt();
    let l2n = n.length_sq();
    if l2n < ZERO_TOL {
        return 0.0;
    }
    n.dot_product(&r_jl) / (l2_jl.sqrt() * l2n.sqrt())
}

/// The inversion force constant and coefficients (C0, C1, C2).
pub fn calc_inversion_coefficients_and_force_constant(
    at2_atomic_num: u32,
    is_c_bound_to_o: bool,
) -> (f64, f64, f64, f64) {
    let (mut res, c0, c1, c2);
    // the central atom is sp2 carbon, nitrogen or oxygen
    if (6..=8).contains(&at2_atomic_num) {
        c0 = 1.0;
        c1 = -1.0;
        c2 = 0.0;
        res = if is_c_bound_to_o { 50.0 } else { 6.0 };
    } else {
        // group 5 elements are not clearly explained in the UFF paper; this
        // follows MCCCS Towhee's ffuff.F
        let mut w0 = PI / 180.0;
        match at2_atomic_num {
            15 => w0 *= 84.4339,
            33 => w0 *= 86.9735,
            51 => w0 *= 87.7047,
            83 => w0 *= 90.0,
            _ => {}
        }
        c2 = 1.0;
        c1 = -4.0 * w0.cos();
        c0 = -(c1 * w0.cos() + c2 * (2.0 * w0).cos());
        res = 22.0 / (c0 + c1 + c2);
    }
    res /= 3.0;
    (res, c0, c1, c2)
}

#[derive(Debug, Clone)]
pub struct InversionContrib {
    at1: usize,
    at2: usize,
    at3: usize,
    at4: usize,
    force_constant: f64,
    c0: f64,
    c1: f64,
    c2: f64,
}

impl InversionContrib {
    pub fn new(
        idx: [usize; 4],
        at2_atomic_num: u32,
        is_c_bound_to_o: bool,
        oob_force_scaling_factor: f64,
    ) -> Self {
        let (fc, c0, c1, c2) =
            calc_inversion_coefficients_and_force_constant(at2_atomic_num, is_c_bound_to_o);
        InversionContrib {
            at1: idx[0],
            at2: idx[1],
            at3: idx[2],
            at4: idx[3],
            force_constant: oob_force_scaling_factor * fc,
            c0,
            c1,
            c2,
        }
    }
}

impl ForceFieldContrib for InversionContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let cos_y = calculate_cos_y(
            &point3(pos, self.at1),
            &point3(pos, self.at2),
            &point3(pos, self.at3),
            &point3(pos, self.at4),
        );
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y = if sin_y_sq > 0.0 { sin_y_sq.sqrt() } else { 0.0 };
        // cos(2 * W) = 2 * sin(W) * sin(W) - 1
        let cos2_w = 2.0 * sin_y * sin_y - 1.0;
        self.force_constant * (self.c0 + self.c1 * sin_y + self.c2 * cos2_w)
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let r_ji = point3(pos, self.at1) - point3(pos, self.at2);
        let r_jk = point3(pos, self.at3) - point3(pos, self.at2);
        let r_jl = point3(pos, self.at4) - point3(pos, self.at2);
        let (d_ji, d_jk, d_jl) = (r_ji.length(), r_jk.length(), r_jl.length());
        if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
            return;
        }
        let (r_ji, r_jk, r_jl) = (r_ji / d_ji, r_jk / d_jk, r_jl / d_jl);
        let mut n = (-r_ji).cross_product(&r_jk);
        n /= n.length();
        let cos_y = clip_to_one(n.dot_product(&r_jl));
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y = sin_y_sq.sqrt().max(1.0e-8);
        let cos_theta = clip_to_one(r_ji.dot_product(&r_jk));
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = sin_theta_sq.sqrt().max(1.0e-8);
        // sin(2 * W) = 2 * sin(W) * cos(W) = 2 * cos(Y) * sin(Y)
        let de_dw = -self.force_constant * (self.c1 * cos_y + 4.0 * self.c2 * cos_y * sin_y);
        let t1 = r_jl.cross_product(&r_jk);
        let t2 = r_ji.cross_product(&r_jl);
        let t3 = r_jk.cross_product(&r_ji);
        let term1 = sin_y * sin_theta;
        let term2 = cos_y / (sin_y * sin_theta_sq);
        let tg1 = [
            (t1.x / term1 - (r_ji.x - r_jk.x * cos_theta) * term2) / d_ji,
            (t1.y / term1 - (r_ji.y - r_jk.y * cos_theta) * term2) / d_ji,
            (t1.z / term1 - (r_ji.z - r_jk.z * cos_theta) * term2) / d_ji,
        ];
        let tg3 = [
            (t2.x / term1 - (r_jk.x - r_ji.x * cos_theta) * term2) / d_jk,
            (t2.y / term1 - (r_jk.y - r_ji.y * cos_theta) * term2) / d_jk,
            (t2.z / term1 - (r_jk.z - r_ji.z * cos_theta) * term2) / d_jk,
        ];
        let tg4 = [
            (t3.x / term1 - r_jl.x * cos_y / sin_y) / d_jl,
            (t3.y / term1 - r_jl.y * cos_y / sin_y) / d_jl,
            (t3.z / term1 - r_jl.z * cos_y / sin_y) / d_jl,
        ];
        for i in 0..3 {
            grad[3 * self.at1 + i] += de_dw * tg1[i];
            grad[3 * self.at2 + i] += -de_dw * (tg1[i] + tg3[i] + tg4[i]);
            grad[3 * self.at3 + i] += de_dw * tg3[i];
            grad[3 * self.at4 + i] += de_dw * tg4[i];
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ------------------------------------------------------------------- vdW

/// The position of the vdW minimum for a pair of atoms.
pub fn calc_nonbonded_minimum(at1: &AtomicParams, at2: &AtomicParams) -> f64 {
    (at1.x1 * at2.x1).sqrt()
}

/// The depth of the vdW well for a pair of atoms.
pub fn calc_nonbonded_depth(at1: &AtomicParams, at2: &AtomicParams) -> f64 {
    (at1.d1 * at2.d1).sqrt()
}

/// A Lennard-Jones 12-6 term; beyond `thresh` (a multiple of the preferred
/// distance) it contributes neither energy nor gradient.
#[derive(Debug, Clone)]
pub struct VdwContrib {
    at1: usize,
    at2: usize,
    xij: f64,
    well_depth: f64,
    thresh: f64,
}

impl VdwContrib {
    pub fn new(
        idx1: usize,
        idx2: usize,
        at1: &AtomicParams,
        at2: &AtomicParams,
        thresh_multiplier: f64,
    ) -> Self {
        // UFF uses the geometric mean of the vdW parameters
        let xij = calc_nonbonded_minimum(at1, at2);
        VdwContrib {
            at1: idx1,
            at2: idx2,
            xij,
            well_depth: calc_nonbonded_depth(at1, at2),
            thresh: thresh_multiplier * xij,
        }
    }
}

impl ForceFieldContrib for VdwContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let d = dist(pos, self.at1, self.at2);
        if d > self.thresh || d <= 0.0 {
            return 0.0;
        }
        let r = self.xij / d;
        let r6 = r.powi(6);
        let r12 = r6 * r6;
        self.well_depth * (r12 - 2.0 * r6)
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let d = dist(pos, self.at1, self.at2);
        if d > self.thresh {
            return;
        }
        if d <= 0.0 {
            for i in 0..3 {
                // move in an arbitrary direction
                grad[3 * self.at1 + i] += 100.0;
                grad[3 * self.at2 + i] -= 100.0;
            }
            return;
        }
        let r = self.xij / d;
        let r7 = r.powi(7);
        let r13 = r.powi(13);
        let pre_factor = 12. * self.well_depth / self.xij * (r7 - r13);
        for i in 0..3 {
            let d_grad = pre_factor * (pos[3 * self.at1 + i] - pos[3 * self.at2 + i]) / d;
            grad[3 * self.at1 + i] += d_grad;
            grad[3 * self.at2 + i] -= d_grad;
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}
