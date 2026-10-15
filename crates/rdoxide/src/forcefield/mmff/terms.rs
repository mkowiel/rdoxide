//! The MMFF energy terms (upstream `ForceField/MMFF/{BondStretch,AngleBend,
//! StretchBend,OopBend,TorsionAngle,Nonbonded}.cpp`).
//!
//! Like upstream, each contribution holds all the terms of its kind.

use super::params::{Angle, Bond, Oop, Prop, Stbn, Tor, VdW, VdWConstants};
use crate::forcefield::uff::{DEG2RAD, RAD2DEG};
use crate::forcefield::{ForceFieldContrib, compute_dihedral, point3};
use crate::geometry::Point3D;

/// conversion of millidyne-angstrom to kcal/mol
pub const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;

fn clip_to_one(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

fn is_double_zero(x: f64) -> bool {
    x < 1.0e-10 && x > -1.0e-10
}

fn dist(pos: &[f64], i: usize, j: usize) -> f64 {
    (point3(pos, i) - point3(pos, j)).length()
}

fn add3_dim(grad: &mut [f64], dim: usize, idx: usize, v: [f64; 3]) {
    for (k, x) in v.iter().enumerate() {
        grad[dim * idx + k] += x;
    }
}

fn add3(grad: &mut [f64], idx: usize, v: [f64; 3]) {
    for (k, x) in v.iter().enumerate() {
        grad[3 * idx + k] += x;
    }
}

// ------------------------------------------------------------ bond stretch

/// `Utils::calcBondStretchEnergy`
pub fn calc_bond_stretch_energy(r0: f64, kb: f64, distance: f64) -> f64 {
    let dist_term = distance - r0;
    let dist_term2 = dist_term * dist_term;
    let c1 = MDYNE_A_TO_KCAL_MOL;
    let cs = -2.0;
    let c3 = 7.0 / 12.0;
    0.5 * c1 * kb * dist_term2 * (1.0 + cs * dist_term + c3 * cs * cs * dist_term2)
}

#[derive(Debug, Clone)]
struct BondTerm {
    at1: usize,
    at2: usize,
    r0: f64,
    kb: f64,
}

#[derive(Debug, Clone, Default)]
pub struct BondStretchContrib {
    terms: Vec<BondTerm>,
}

impl BondStretchContrib {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(&mut self, idx1: usize, idx2: usize, params: &Bond) {
        self.terms.push(BondTerm {
            at1: idx1,
            at2: idx2,
            r0: params.r0,
            kb: params.kb,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for BondStretchContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|t| calc_bond_stretch_energy(t.r0, t.kb, dist(pos, t.at1, t.at2)))
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let cs = -2.0;
        let c1 = MDYNE_A_TO_KCAL_MOL;
        let c3 = 7.0 / 12.0;
        for t in &self.terms {
            let d = dist(pos, t.at1, t.at2);
            let dist_term = d - t.r0;
            let de_dr = c1
                * t.kb
                * dist_term
                * (1.0 + 1.5 * cs * dist_term + 2.0 * c3 * cs * cs * dist_term * dist_term);
            for i in 0..3 {
                let d_grad = if d > 0.0 {
                    de_dr * (pos[3 * t.at1 + i] - pos[3 * t.at2 + i]) / d
                } else {
                    t.kb * 0.01
                };
                grad[3 * t.at1 + i] += d_grad;
                grad[3 * t.at2 + i] -= d_grad;
            }
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// -------------------------------------------------------------- angle bend

/// `Utils::calcCosTheta`
pub fn calc_cos_theta(p1: &Point3D, p2: &Point3D, p3: &Point3D, dist1: f64, dist2: f64) -> f64 {
    let p12 = *p1 - *p2;
    let p32 = *p3 - *p2;
    clip_to_one(p12.dot_product(&p32) / (dist1 * dist2))
}

/// `Utils::calcAngleBendEnergy`
pub fn calc_angle_bend_energy(theta0: f64, ka: f64, is_linear: bool, cos_theta: f64) -> f64 {
    let angle = RAD2DEG * cos_theta.acos() - theta0;
    let cb = -0.006981317;
    let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    if is_linear {
        MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cos_theta)
    } else {
        0.5 * c2 * ka * angle * angle * (1.0 + cb * angle)
    }
}

#[derive(Debug, Clone)]
struct AngleTerm {
    at: [usize; 3],
    is_linear: bool,
    theta0: f64,
    ka: f64,
}

#[derive(Debug, Clone, Default)]
pub struct AngleBendContrib {
    terms: Vec<AngleTerm>,
}

impl AngleBendContrib {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(&mut self, at: [usize; 3], params: &Angle, central_atom: &Prop) {
        assert!(
            at[0] != at[1] && at[1] != at[2] && at[0] != at[2],
            "degenerate points"
        );
        self.terms.push(AngleTerm {
            at,
            is_linear: central_atom.linh > 0,
            theta0: params.theta0,
            ka: params.ka,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for AngleBendContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|t| {
                let [a, b, c] = t.at;
                let (d1, d2) = (dist(pos, a, b), dist(pos, b, c));
                let cos_theta =
                    calc_cos_theta(&point3(pos, a), &point3(pos, b), &point3(pos, c), d1, d2);
                calc_angle_bend_energy(t.theta0, t.ka, t.is_linear, cos_theta)
            })
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for t in &self.terms {
            let [a, b, c] = t.at;
            let d = [dist(pos, a, b), dist(pos, b, c)];
            let (p1, p2, p3) = (point3(pos, a), point3(pos, b), point3(pos, c));
            let r = [(p1 - p2) / d[0], (p3 - p2) / d[1]];
            let cos_theta = clip_to_one(r[0].dot_product(&r[1]));
            let sin_theta_sq = 1.0 - cos_theta * cos_theta;
            let sin_theta = (if sin_theta_sq > 0.0 {
                sin_theta_sq.sqrt()
            } else {
                0.0
            })
            .max(1.0e-8);
            let angle_term = RAD2DEG * cos_theta.acos() - t.theta0;
            let cb = -0.006981317;
            let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
            let de_dtheta = if t.is_linear {
                -MDYNE_A_TO_KCAL_MOL * t.ka * sin_theta
            } else {
                RAD2DEG * c2 * t.ka * angle_term * (1.0 + 1.5 * cb * angle_term)
            };
            let d_cos = [
                1.0 / d[0] * (r[1].x - cos_theta * r[0].x),
                1.0 / d[0] * (r[1].y - cos_theta * r[0].y),
                1.0 / d[0] * (r[1].z - cos_theta * r[0].z),
                1.0 / d[1] * (r[0].x - cos_theta * r[1].x),
                1.0 / d[1] * (r[0].y - cos_theta * r[1].y),
                1.0 / d[1] * (r[0].z - cos_theta * r[1].z),
            ];
            let f = |x: f64| de_dtheta * x / (-sin_theta);
            add3(grad, a, [f(d_cos[0]), f(d_cos[1]), f(d_cos[2])]);
            add3(
                grad,
                b,
                [
                    f(-d_cos[0] - d_cos[3]),
                    f(-d_cos[1] - d_cos[4]),
                    f(-d_cos[2] - d_cos[5]),
                ],
            );
            add3(grad, c, [f(d_cos[3]), f(d_cos[4]), f(d_cos[5])]);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ----------------------------------------------------------- stretch-bend

#[derive(Debug, Clone)]
struct StretchBendTerm {
    at: [usize; 3],
    rest_len: [f64; 2],
    theta0: f64,
    force_constants: [f64; 2],
}

#[derive(Debug, Clone, Default)]
pub struct StretchBendContrib {
    terms: Vec<StretchBendTerm>,
}

impl StretchBendContrib {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(
        &mut self,
        at: [usize; 3],
        stbn: &Stbn,
        angle: &Angle,
        bond1: &Bond,
        bond2: &Bond,
    ) {
        assert!(
            at[0] != at[1] && at[1] != at[2] && at[0] != at[2],
            "degenerate points"
        );
        self.terms.push(StretchBendTerm {
            at,
            rest_len: [bond1.r0, bond2.r0],
            theta0: angle.theta0,
            force_constants: [stbn.kba_ijk, stbn.kba_kji],
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for StretchBendContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|t| {
                let [a, b, c] = t.at;
                let (d1, d2) = (dist(pos, a, b), dist(pos, b, c));
                let cos_theta =
                    calc_cos_theta(&point3(pos, a), &point3(pos, b), &point3(pos, c), d1, d2);
                let delta_theta = RAD2DEG * cos_theta.acos() - t.theta0;
                let factor = MDYNE_A_TO_KCAL_MOL * DEG2RAD * delta_theta;
                factor * t.force_constants[0] * (d1 - t.rest_len[0])
                    + factor * t.force_constants[1] * (d2 - t.rest_len[1])
            })
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let c5 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
        for t in &self.terms {
            let [a, b, c] = t.at;
            let (dist1, dist2) = (dist(pos, a, b), dist(pos, b, c));
            let (p1, p2, p3) = (point3(pos, a), point3(pos, b), point3(pos, c));
            let p12 = (p1 - p2) / dist1;
            let p32 = (p3 - p2) / dist2;
            let (fc1, fc2) = (t.force_constants[0], t.force_constants[1]);
            let cos_theta = clip_to_one(p12.dot_product(&p32));
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1.0e-8);
            let angle_term = RAD2DEG * cos_theta.acos() - t.theta0;
            let dist_term =
                RAD2DEG * (fc1 * (dist1 - t.rest_len[0]) + fc2 * (dist2 - t.rest_len[1]));
            let d_cos1 = [
                1.0 / dist1 * (p32.x - cos_theta * p12.x),
                1.0 / dist1 * (p32.y - cos_theta * p12.y),
                1.0 / dist1 * (p32.z - cos_theta * p12.z),
            ];
            let d_cos2 = [
                1.0 / dist2 * (p12.x - cos_theta * p32.x),
                1.0 / dist2 * (p12.y - cos_theta * p32.y),
                1.0 / dist2 * (p12.z - cos_theta * p32.z),
            ];
            let u12 = [p12.x, p12.y, p12.z];
            let u32 = [p32.x, p32.y, p32.z];
            let mut g1 = [0.0; 3];
            let mut g2 = [0.0; 3];
            let mut g3 = [0.0; 3];
            for k in 0..3 {
                g1[k] = c5 * (u12[k] * fc1 * angle_term + d_cos1[k] / (-sin_theta) * dist_term);
                g2[k] = c5
                    * ((-u12[k] * fc1 - u32[k] * fc2) * angle_term
                        + (-d_cos1[k] - d_cos2[k]) / (-sin_theta) * dist_term);
                g3[k] = c5 * (u32[k] * fc2 * angle_term + d_cos2[k] / (-sin_theta) * dist_term);
            }
            add3(grad, a, g1);
            add3(grad, b, g2);
            add3(grad, c, g3);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// -------------------------------------------------------- out-of-plane bend

/// `Utils::calcOopChi`: the out-of-plane angle (degrees) of l with respect to
/// the plane i-j-k, j being the center.
pub fn calc_oop_chi(i: &Point3D, j: &Point3D, k: &Point3D, l: &Point3D) -> f64 {
    let mut r_ji = *i - *j;
    let mut r_jk = *k - *j;
    let mut r_jl = *l - *j;
    r_ji /= r_ji.length();
    r_jk /= r_jk.length();
    r_jl /= r_jl.length();
    let mut n = r_ji.cross_product(&r_jk);
    n /= n.length();
    let sin_chi = clip_to_one(n.dot_product(&r_jl));
    RAD2DEG * sin_chi.asin()
}

/// `Utils::calcOopBendEnergy`
pub fn calc_oop_bend_energy(chi: f64, koop: f64) -> f64 {
    let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    0.5 * c2 * koop * chi * chi
}

#[derive(Debug, Clone)]
struct OopTerm {
    at: [usize; 4],
    koop: f64,
}

#[derive(Debug, Clone, Default)]
pub struct OopBendContrib {
    terms: Vec<OopTerm>,
}

impl OopBendContrib {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(&mut self, at: [usize; 4], params: &Oop) {
        for a in 0..4 {
            for b in a + 1..4 {
                assert!(at[a] != at[b], "degenerate points");
            }
        }
        self.terms.push(OopTerm {
            at,
            koop: params.koop,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for OopBendContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|t| {
                let p = t.at.map(|a| point3(pos, a));
                calc_oop_bend_energy(calc_oop_chi(&p[0], &p[1], &p[2], &p[3]), t.koop)
            })
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
        for t in &self.terms {
            let [i, j, k, l] = t.at;
            let (ip, jp, kp, lp) = (
                point3(pos, i),
                point3(pos, j),
                point3(pos, k),
                point3(pos, l),
            );
            let mut r_ji = ip - jp;
            let mut r_jk = kp - jp;
            let mut r_jl = lp - jp;
            let (d_ji, d_jk, d_jl) = (r_ji.length(), r_jk.length(), r_jl.length());
            if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
                continue;
            }
            r_ji /= d_ji;
            r_jk /= d_jk;
            r_jl /= d_jl;
            let mut n = (-r_ji).cross_product(&r_jk);
            n /= n.length();
            let sin_chi = clip_to_one(r_jl.dot_product(&n));
            let cos_chi_sq = 1.0 - sin_chi * sin_chi;
            let cos_chi = (if cos_chi_sq > 0.0 {
                cos_chi_sq.sqrt()
            } else {
                0.0
            })
            .max(1.0e-8);
            let chi = RAD2DEG * sin_chi.asin();
            let cos_theta = clip_to_one(r_ji.dot_product(&r_jk));
            let sin_theta_sq = (1.0 - cos_theta * cos_theta).max(1.0e-8);
            let sin_theta = (if sin_theta_sq > 0.0 {
                sin_theta_sq.sqrt()
            } else {
                0.0
            })
            .max(1.0e-8);

            let de_dchi = RAD2DEG * c2 * t.koop * chi;
            let t1 = r_jl.cross_product(&r_jk);
            let t2 = r_ji.cross_product(&r_jl);
            let t3 = r_jk.cross_product(&r_ji);
            let term1 = cos_chi * sin_theta;
            let term2 = sin_chi / (cos_chi * sin_theta_sq);
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
                (t3.x / term1 - r_jl.x * sin_chi / cos_chi) / d_jl,
                (t3.y / term1 - r_jl.y * sin_chi / cos_chi) / d_jl,
                (t3.z / term1 - r_jl.z * sin_chi / cos_chi) / d_jl,
            ];
            for c in 0..3 {
                grad[3 * i + c] += de_dchi * tg1[c];
                grad[3 * j + c] += -de_dchi * (tg1[c] + tg3[c] + tg4[c]);
                grad[3 * k + c] += de_dchi * tg3[c];
                grad[3 * l + c] += de_dchi * tg4[c];
            }
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ---------------------------------------------------------------- torsion

/// `Utils::calcTorsionCosPhi`
pub fn calc_torsion_cos_phi(i: &Point3D, j: &Point3D, k: &Point3D, l: &Point3D) -> f64 {
    let r1 = *i - *j;
    let r2 = *k - *j;
    let r3 = *j - *k;
    let r4 = *l - *k;
    let t1 = r1.cross_product(&r2);
    let t2 = r3.cross_product(&r4);
    let (t1_len, t2_len) = (t1.length(), t2.length());
    if is_double_zero(t1_len) || is_double_zero(t2_len) {
        return 0.0;
    }
    clip_to_one(t1.dot_product(&t2) / (t1_len * t2_len))
}

/// `Utils::calcTorsionEnergy`
pub fn calc_torsion_energy(v1: f64, v2: f64, v3: f64, cos_phi: f64) -> f64 {
    let cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
    let cos3_phi = cos_phi * (2.0 * cos2_phi - 1.0);
    0.5 * (v1 * (1.0 + cos_phi) + v2 * (1.0 - cos2_phi) + v3 * (1.0 + cos3_phi))
}

/// Adds the gradient of a torsion angle, given the pieces of its dihedral and
/// `sin_term = -dE/dPhi / sin(phi)`, to the four atoms `at`.
pub(crate) fn add_torsion_grad(
    grad: &mut [f64],
    at: [usize; 4],
    r: &[Point3D; 4],
    tv: &[Point3D; 2],
    d: &[f64; 2],
    cos_phi: f64,
    s: f64,
    dim: usize,
) {
    let dc = [
        1.0 / d[0] * (tv[1].x - cos_phi * tv[0].x),
        1.0 / d[0] * (tv[1].y - cos_phi * tv[0].y),
        1.0 / d[0] * (tv[1].z - cos_phi * tv[0].z),
        1.0 / d[1] * (tv[0].x - cos_phi * tv[1].x),
        1.0 / d[1] * (tv[0].y - cos_phi * tv[1].y),
        1.0 / d[1] * (tv[0].z - cos_phi * tv[1].z),
    ];
    add3_dim(
        grad,
        dim,
        at[0],
        [
            s * (dc[2] * r[1].y - dc[1] * r[1].z),
            s * (dc[0] * r[1].z - dc[2] * r[1].x),
            s * (dc[1] * r[1].x - dc[0] * r[1].y),
        ],
    );
    add3_dim(
        grad,
        dim,
        at[1],
        [
            s * (dc[1] * (r[1].z - r[0].z)
                + dc[2] * (r[0].y - r[1].y)
                + dc[4] * (-r[3].z)
                + dc[5] * r[3].y),
            s * (dc[0] * (r[0].z - r[1].z)
                + dc[2] * (r[1].x - r[0].x)
                + dc[3] * r[3].z
                + dc[5] * (-r[3].x)),
            s * (dc[0] * (r[1].y - r[0].y)
                + dc[1] * (r[0].x - r[1].x)
                + dc[3] * (-r[3].y)
                + dc[4] * r[3].x),
        ],
    );
    add3_dim(
        grad,
        dim,
        at[2],
        [
            s * (dc[1] * r[0].z
                + dc[2] * (-r[0].y)
                + dc[4] * (r[3].z - r[2].z)
                + dc[5] * (r[2].y - r[3].y)),
            s * (dc[0] * (-r[0].z)
                + dc[2] * r[0].x
                + dc[3] * (r[2].z - r[3].z)
                + dc[5] * (r[3].x - r[2].x)),
            s * (dc[0] * r[0].y
                + dc[1] * (-r[0].x)
                + dc[3] * (r[3].y - r[2].y)
                + dc[4] * (r[2].x - r[3].x)),
        ],
    );
    add3_dim(
        grad,
        dim,
        at[3],
        [
            s * (dc[4] * r[2].z - dc[5] * r[2].y),
            s * (dc[5] * r[2].x - dc[3] * r[2].z),
            s * (dc[3] * r[2].y - dc[4] * r[2].x),
        ],
    );
}

#[derive(Debug, Clone)]
struct TorsionTerm {
    at: [usize; 4],
    v: [f64; 3],
}

#[derive(Debug, Clone, Default)]
pub struct TorsionAngleContrib {
    terms: Vec<TorsionTerm>,
}

impl TorsionAngleContrib {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(&mut self, at: [usize; 4], params: &Tor) {
        for a in 0..4 {
            for b in a + 1..4 {
                assert!(at[a] != at[b], "degenerate points");
            }
        }
        self.terms.push(TorsionTerm {
            at,
            v: [params.v1, params.v2, params.v3],
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for TorsionAngleContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|t| {
                let p = t.at.map(|a| point3(pos, a));
                calc_torsion_energy(
                    t.v[0],
                    t.v[1],
                    t.v[2],
                    calc_torsion_cos_phi(&p[0], &p[1], &p[2], &p[3]),
                )
            })
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for t in &self.terms {
            let p = t.at.map(|a| point3(pos, a));
            let parts = compute_dihedral(&p[0], &p[1], &p[2], &p[3]);
            let (r, tv, d, cos_phi) = (&parts.r, &parts.t, &parts.d, parts.cos_phi);
            let sin_phi_sq = 1.0 - cos_phi * cos_phi;
            let sin_phi = if sin_phi_sq > 0.0 {
                sin_phi_sq.sqrt()
            } else {
                0.0
            };
            let sin2_phi = 2.0 * sin_phi * cos_phi;
            let sin3_phi = 3.0 * sin_phi - 4.0 * sin_phi * sin_phi_sq;
            // dE/dPhi is independent of the cartesians
            let de_dphi =
                0.5 * (-t.v[0] * sin_phi + 2.0 * t.v[1] * sin2_phi - 3.0 * t.v[2] * sin3_phi);
            // the tolerance is the recommendation of Niketic and Rasmussen
            let sin_term = -de_dphi
                * (if is_double_zero(sin_phi) {
                    1.0 / cos_phi
                } else {
                    1.0 / sin_phi
                });
            add_torsion_grad(grad, t.at, r, tv, d, cos_phi, sin_term, 3);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ------------------------------------------------------------- nonbonded

/// The van der Waals parameters of an atom pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VdWRijstarEps {
    pub r_ij_star_unscaled: f64,
    pub epsilon_unscaled: f64,
    pub r_ij_star: f64,
    pub epsilon: f64,
}

/// `Utils::calcUnscaledVdWMinimum`
pub fn calc_unscaled_vdw_minimum(vdw: &VdWConstants, i: &VdW, j: &VdW) -> f64 {
    let gamma_ij = (i.r_star - j.r_star) / (i.r_star + j.r_star);
    0.5 * (i.r_star + j.r_star)
        * (1.0
            + if i.da == b'D' || j.da == b'D' {
                0.0
            } else {
                vdw.b * (1.0 - (-vdw.beta * gamma_ij * gamma_ij).exp())
            })
}

/// `Utils::calcUnscaledVdWWellDepth`
pub fn calc_unscaled_vdw_well_depth(r_star_ij: f64, i: &VdW, j: &VdW) -> f64 {
    let r_star_ij2 = r_star_ij * r_star_ij;
    let c4 = 181.16;
    c4 * i.g_i * j.g_i * i.alpha_i * j.alpha_i
        / (((i.alpha_i / i.n_i).sqrt() + (j.alpha_i / j.n_i).sqrt())
            * r_star_ij2
            * r_star_ij2
            * r_star_ij2)
}

/// `Utils::scaleVdWParams`: donor-acceptor pairs get a shorter, deeper minimum.
pub fn scale_vdw_params(
    r_star_ij: &mut f64,
    well_depth: &mut f64,
    vdw: &VdWConstants,
    i: &VdW,
    j: &VdW,
) {
    if (i.da == b'D' && j.da == b'A') || (i.da == b'A' && j.da == b'D') {
        *r_star_ij *= vdw.darad;
        *well_depth *= vdw.daeps;
    }
}

impl VdWRijstarEps {
    /// The parameters of the pair, scaled for donor-acceptor pairs.
    pub fn new(vdw: &VdWConstants, i: &VdW, j: &VdW) -> Self {
        let r_ij_star_unscaled = calc_unscaled_vdw_minimum(vdw, i, j);
        let epsilon_unscaled = calc_unscaled_vdw_well_depth(r_ij_star_unscaled, i, j);
        let (mut r_ij_star, mut epsilon) = (r_ij_star_unscaled, epsilon_unscaled);
        scale_vdw_params(&mut r_ij_star, &mut epsilon, vdw, i, j);
        VdWRijstarEps {
            r_ij_star_unscaled,
            epsilon_unscaled,
            r_ij_star,
            epsilon,
        }
    }
}

/// `Utils::calcVdWEnergy`
pub fn calc_vdw_energy(dist: f64, r_star_ij: f64, well_depth: f64) -> f64 {
    let vdw1 = 1.07;
    let vdw1m1 = vdw1 - 1.0;
    let vdw2 = 1.12;
    let vdw2m1 = vdw2 - 1.0;
    let dist2 = dist * dist;
    let dist7 = dist2 * dist2 * dist2 * dist;
    let a_term = vdw1 * r_star_ij / (dist + vdw1m1 * r_star_ij);
    let a_term2 = a_term * a_term;
    let a_term7 = a_term2 * a_term2 * a_term2 * a_term;
    let r_star_ij2 = r_star_ij * r_star_ij;
    let r_star_ij7 = r_star_ij2 * r_star_ij2 * r_star_ij2 * r_star_ij;
    let b_term = vdw2 * r_star_ij7 / (dist7 + vdw2m1 * r_star_ij7) - 2.0;
    well_depth * a_term7 * b_term
}

/// `RDKit::MMFF::CONSTANT` / `DISTANCE`: the dielectric model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DielModel {
    Constant = 1,
    Distance = 2,
}

/// `Utils::calcEleEnergy`
pub fn calc_ele_energy(dist: f64, charge_term: f64, diel_model: DielModel, is_1_4: bool) -> f64 {
    let mut corr_dist = dist + 0.05;
    let diel = 332.0716;
    let sc1_4 = 0.75;
    if diel_model == DielModel::Distance {
        corr_dist *= corr_dist;
    }
    diel * charge_term / corr_dist * if is_1_4 { sc1_4 } else { 1.0 }
}

/// dE/dr of the van der Waals term
fn vdw_de_dr(dist: f64, r_star: f64, well_depth: f64) -> f64 {
    let vdw1 = 1.07;
    let vdw1m1 = vdw1 - 1.0;
    let vdw2 = 1.12;
    let vdw2m1 = vdw2 - 1.0;
    let vdw2t7 = vdw2 * 7.0;
    let q = dist / r_star;
    let q2 = q * q;
    let q6 = q2 * q2 * q2;
    let q7 = q6 * q;
    let q7pvdw2m1 = q7 + vdw2m1;
    let t = vdw1 / (q + vdw1 - 1.0);
    let t2 = t * t;
    let t7 = t2 * t2 * t2 * t;
    well_depth / r_star
        * t7
        * (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1) + ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)))
}

/// dE/dr of the electrostatic term
fn ele_de_dr(dist: f64, charge_term: f64, diel_model: DielModel, is_1_4: bool) -> f64 {
    let mut corr_dist = dist + 0.05;
    corr_dist *= if diel_model == DielModel::Distance {
        corr_dist * corr_dist
    } else {
        corr_dist
    };
    -332.0716 * f64::from(diel_model as u8) * charge_term / corr_dist
        * if is_1_4 { 0.75 } else { 1.0 }
}

/// Applies a radial gradient `de_dr` along the vector between two atoms, or
/// `fallback` along each axis when the atoms coincide.
fn add_radial_grad(
    pos: &[f64],
    grad: &mut [f64],
    at1: usize,
    at2: usize,
    d: f64,
    de_dr: f64,
    fallback: f64,
) {
    for i in 0..3 {
        let d_grad = if d > 0.0 {
            de_dr * (pos[3 * at1 + i] - pos[3 * at2 + i]) / d
        } else {
            fallback
        };
        grad[3 * at1 + i] += d_grad;
        grad[3 * at2 + i] -= d_grad;
    }
}

#[derive(Debug, Clone)]
struct VdWTerm {
    at1: usize,
    at2: usize,
    r_ij_star: f64,
    well_depth: f64,
}

#[derive(Debug, Clone, Default)]
pub struct VdWContrib {
    terms: Vec<VdWTerm>,
}

impl VdWContrib {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(&mut self, idx1: usize, idx2: usize, constants: &VdWRijstarEps) {
        self.terms.push(VdWTerm {
            at1: idx1,
            at2: idx2,
            r_ij_star: constants.r_ij_star,
            well_depth: constants.epsilon,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for VdWContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|t| calc_vdw_energy(dist(pos, t.at1, t.at2), t.r_ij_star, t.well_depth))
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for t in &self.terms {
            let d = dist(pos, t.at1, t.at2);
            let de_dr = vdw_de_dr(d, t.r_ij_star, t.well_depth);
            add_radial_grad(pos, grad, t.at1, t.at2, d, de_dr, t.r_ij_star * 0.01);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

#[derive(Debug, Clone)]
struct EleTerm {
    at1: usize,
    at2: usize,
    charge_term: f64,
    diel_model: DielModel,
    is_1_4: bool,
}

#[derive(Debug, Clone, Default)]
pub struct EleContrib {
    terms: Vec<EleTerm>,
}

impl EleContrib {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_term(
        &mut self,
        idx1: usize,
        idx2: usize,
        charge_term: f64,
        diel_model: DielModel,
        is_1_4: bool,
    ) {
        self.terms.push(EleTerm {
            at1: idx1,
            at2: idx2,
            charge_term,
            diel_model,
            is_1_4,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for EleContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|t| {
                calc_ele_energy(
                    dist(pos, t.at1, t.at2),
                    t.charge_term,
                    t.diel_model,
                    t.is_1_4,
                )
            })
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for t in &self.terms {
            let d = dist(pos, t.at1, t.at2);
            let de_dr = ele_de_dr(d, t.charge_term, t.diel_model, t.is_1_4);
            add_radial_grad(pos, grad, t.at1, t.at2, d, de_dr, 0.02);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

#[derive(Debug, Clone)]
struct NonbondedTerm {
    at1: usize,
    at2: usize,
    /// van der Waals minimum and well depth, if the pair has a vdW term
    vdw: Option<(f64, f64)>,
    /// charge term, dielectric model and 1-4 flag, if the pair is charged
    ele: Option<(f64, DielModel, bool)>,
}

/// The van der Waals and electrostatic terms of atom pairs in one pass.
#[derive(Debug, Clone, Default)]
pub struct NonbondedContrib {
    terms: Vec<NonbondedTerm>,
}

impl NonbondedContrib {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a pair; nothing is added when it has neither term.
    pub fn add_term(
        &mut self,
        idx1: usize,
        idx2: usize,
        vdw: Option<&VdWRijstarEps>,
        ele: Option<(f64, DielModel, bool)>,
    ) {
        if vdw.is_none() && ele.is_none() {
            return;
        }
        self.terms.push(NonbondedTerm {
            at1: idx1,
            at2: idx2,
            vdw: vdw.map(|v| (v.r_ij_star, v.epsilon)),
            ele,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl ForceFieldContrib for NonbondedContrib {
    fn energy(&self, pos: &[f64]) -> f64 {
        let mut sum = 0.0;
        for t in &self.terms {
            let d = dist(pos, t.at1, t.at2);
            if let Some((r_star, depth)) = t.vdw {
                sum += calc_vdw_energy(d, r_star, depth);
            }
            if let Some((charge_term, model, is_1_4)) = t.ele {
                sum += calc_ele_energy(d, charge_term, model, is_1_4);
            }
        }
        sum
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for t in &self.terms {
            let d = dist(pos, t.at1, t.at2);
            if d <= 0.0 {
                // coincident atoms: nudge them apart (upstream returns from
                // the whole loop here and skips the remaining pairs)
                if let Some((r_star, _)) = t.vdw {
                    for i in 0..3 {
                        grad[3 * t.at1 + i] += r_star * 0.01;
                        grad[3 * t.at2 + i] -= r_star * 0.01;
                    }
                }
                if t.ele.is_some() {
                    for i in 0..3 {
                        grad[3 * t.at1 + i] += 0.02;
                        grad[3 * t.at2 + i] -= 0.02;
                    }
                }
                continue;
            }
            let mut de_dr = 0.0;
            if let Some((r_star, depth)) = t.vdw {
                de_dr += vdw_de_dr(d, r_star, depth) / d;
            }
            if let Some((charge_term, model, is_1_4)) = t.ele {
                de_dr += ele_de_dr(d, charge_term, model, is_1_4) / d;
            }
            for i in 0..3 {
                let d_grad = de_dr * (pos[3 * t.at1 + i] - pos[3 * t.at2 + i]);
                grad[3 * t.at1 + i] += d_grad;
                grad[3 * t.at2 + i] -= d_grad;
            }
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::ForceField;

    const POS: [f64; 15] = [
        0.1, 0.2, 0.0, 1.5, 0.3, 0.2, 2.1, 1.6, -0.1, 3.4, 1.9, 0.9, 1.0, -1.1, 0.7,
    ];

    fn check_gradient(contrib: &dyn ForceFieldContrib) {
        let mut grad = vec![0.0; POS.len()];
        contrib.add_grad(&POS, &mut grad);
        for i in 0..POS.len() {
            let h = 1e-6;
            let (mut plus, mut minus) = (POS.to_vec(), POS.to_vec());
            plus[i] += h;
            minus[i] -= h;
            let numeric = (contrib.energy(&plus) - contrib.energy(&minus)) / (2.0 * h);
            assert!(
                (numeric - grad[i]).abs() < 1e-4 * (1.0 + numeric.abs()),
                "component {i}: analytic {} numeric {numeric}",
                grad[i]
            );
        }
    }

    fn prop(linh: u8) -> Prop {
        Prop {
            atno: 6,
            crd: 4,
            val: 4,
            pilp: 0,
            mltb: 0,
            arom: 0,
            linh,
            sbmb: 0,
        }
    }

    #[test]
    fn bond_stretch() {
        let mut c = BondStretchContrib::new();
        c.add_term(
            0,
            1,
            &Bond {
                kb: 4.258,
                r0: 1.508,
            },
        );
        c.add_term(1, 2, &Bond { kb: 5.0, r0: 1.4 });
        check_gradient(&c);
        // zero at the rest length
        assert_eq!(calc_bond_stretch_energy(1.5, 4.0, 1.5), 0.0);
    }

    #[test]
    fn angle_bend() {
        let mut c = AngleBendContrib::new();
        let angle = Angle {
            ka: 0.851,
            theta0: 109.6,
        };
        c.add_term([0, 1, 2], &angle, &prop(0));
        c.add_term([1, 2, 3], &angle, &prop(0));
        check_gradient(&c);
        let mut linear = AngleBendContrib::new();
        linear.add_term(
            [0, 1, 2],
            &Angle {
                ka: 0.2,
                theta0: 180.0,
            },
            &prop(1),
        );
        check_gradient(&linear);
    }

    #[test]
    fn stretch_bend() {
        let mut c = StretchBendContrib::new();
        c.add_term(
            [0, 1, 2],
            &Stbn {
                kba_ijk: 0.2,
                kba_kji: 0.1,
            },
            &Angle {
                ka: 0.8,
                theta0: 105.0,
            },
            &Bond { kb: 4.0, r0: 1.4 },
            &Bond { kb: 4.0, r0: 1.3 },
        );
        check_gradient(&c);
    }

    #[test]
    fn oop_bend() {
        let mut c = OopBendContrib::new();
        c.add_term([0, 1, 2, 3], &Oop { koop: 0.03 });
        c.add_term([4, 1, 2, 3], &Oop { koop: 0.02 });
        check_gradient(&c);
    }

    #[test]
    fn torsion() {
        let mut c = TorsionAngleContrib::new();
        c.add_term(
            [0, 1, 2, 3],
            &Tor {
                v1: 0.2,
                v2: -0.8,
                v3: 1.5,
            },
        );
        c.add_term(
            [1, 2, 3, 4],
            &Tor {
                v1: 0.0,
                v2: 0.5,
                v3: 0.1,
            },
        );
        check_gradient(&c);
    }

    #[test]
    fn nonbonded() {
        let vdw = VdWRijstarEps {
            r_ij_star_unscaled: 3.9,
            epsilon_unscaled: 0.07,
            r_ij_star: 3.9,
            epsilon: 0.07,
        };
        let mut v = VdWContrib::new();
        v.add_term(0, 3, &vdw);
        check_gradient(&v);
        for model in [DielModel::Constant, DielModel::Distance] {
            let mut e = EleContrib::new();
            e.add_term(0, 3, 0.3 * -0.2, model, false);
            e.add_term(1, 4, 0.5 * 0.4, model, true);
            check_gradient(&e);
        }
        // the combined contribution is the sum of the two
        let mut n = NonbondedContrib::new();
        n.add_term(0, 3, Some(&vdw), Some((-0.06, DielModel::Distance, true)));
        n.add_term(1, 4, None, None);
        check_gradient(&n);
        let mut e = EleContrib::new();
        e.add_term(0, 3, -0.06, DielModel::Distance, true);
        assert!((n.energy(&POS) - v.energy(&POS) - e.energy(&POS)).abs() < 1e-12);
    }

    #[test]
    fn vdw_parameters_of_carbons() {
        let p = super::super::params::MmffParams::get(false);
        let c = p.vdw(1).unwrap();
        let r = VdWRijstarEps::new(p.vdw_constants(), c, c);
        assert!((r.r_ij_star - c.r_star).abs() < 1e-12);
        assert!(r.epsilon > 0.0);
    }

    #[test]
    fn minimizing_a_bond_reaches_its_rest_length() {
        let mut ff = ForceField::new(3);
        ff.push_point(&[0.0, 0.0, 0.0]);
        ff.push_point(&[2.0, 0.0, 0.0]);
        let mut c = BondStretchContrib::new();
        c.add_term(
            0,
            1,
            &Bond {
                kb: 4.258,
                r0: 1.508,
            },
        );
        ff.push_contrib(Box::new(c));
        assert_eq!(ff.minimize(200, 1e-6, 1e-6).unwrap(), 0);
        assert!((ff.distance(0, 1, None).unwrap() - 1.508).abs() < 1e-3);
    }
}
