//! The torsion and planarity terms of the ETKDG force field (upstream
//! `ForceFieldHelpers/CrystalFF/{TorsionAngleContribs,GaussianTorsionAngleContribs,
//! PlanarityContribs}`). `TorsionAngleM6` of upstream is not used anywhere and
//! computes the same as one entry of [`TorsionAngleContribs`].

use std::f64::consts::PI;

use super::mmff::terms::{add_torsion_grad, calc_torsion_cos_phi};
use super::{ForceFieldContrib, point3};
use crate::geometry::Point3D;

const TWO_PI: f64 = 2.0 * PI;

/// The number of points in the energy and gradient lookup tables of Gaussian
/// torsion terms (angles from 0 to pi).
pub const LOOKUP_GRID_SIZE: usize = 180;

fn is_double_zero(x: f64) -> bool {
    x < 1.0e-10 && x > -1.0e-10
}

fn all_different(at: &[usize; 4]) -> bool {
    (0..4).all(|a| (a + 1..4).all(|b| at[a] != at[b]))
}

/// Piecewise-linear lookup of a table sampled on `[0, pi]`.
fn interpolate(table: &[f64], phi: f64) -> f64 {
    let norm = phi * (table.len() - 1) as f64 / PI;
    let idx = (norm as usize).min(table.len() - 2);
    let t = norm - idx as f64;
    table[idx] * (1.0 - t) + table[idx + 1] * t
}

/// The energy of a sum of periodic Gaussians at the torsion angle `phi`.
pub fn gaussian_torsion_energy(
    heights: &[f64],
    positions: &[f64],
    widths: &[f64],
    phi: f64,
) -> f64 {
    let mut result = 0.0;
    for i in 0..heights.len() {
        let (a, b, c) = (heights[i], positions[i], widths[i]);
        let cinv = 1.0 / (c * c);
        let d1 = phi - b;
        let d2 = TWO_PI - (phi + b);
        let d3 = TWO_PI + d1;
        let d4 = TWO_PI + d2;
        let t1 = (-d1 * d1 * cinv).exp();
        let t2 = (-d2 * d2 * cinv).exp();
        let t3 = (-d3 * d3 * cinv).exp();
        let t4 = (-d4 * d4 * cinv).exp();
        result += a * (t1 + t2 + t3 + t4);
    }
    -result.ln()
}

/// The derivative of [`gaussian_torsion_energy`] with respect to `phi`.
pub fn gaussian_torsion_de_dphi(
    heights: &[f64],
    positions: &[f64],
    widths: &[f64],
    phi: f64,
) -> f64 {
    let mut inner_derivs = 0.0;
    let mut denominator = 0.0;
    for i in 0..heights.len() {
        let (a, b, c) = (heights[i], positions[i], widths[i]);
        let cinv = 1.0 / (c * c);
        let d1 = phi - b;
        let d2 = TWO_PI - (phi + b);
        let d3 = TWO_PI + d1;
        let d4 = TWO_PI + d2;
        let t1 = (-d1 * d1 * cinv).exp();
        let t2 = (-d2 * d2 * cinv).exp();
        let t3 = (-d3 * d3 * cinv).exp();
        let t4 = (-d4 * d4 * cinv).exp();
        inner_derivs += a * cinv * (-t1 * d1 + t2 * d2 - t3 * d3 + t4 * d4);
        denominator += a * (t1 + t2 + t3 + t4);
    }
    if denominator < f64::MIN_POSITIVE {
        denominator = f64::MIN_POSITIVE;
    }
    -2.0 * inner_derivs / denominator
}

// ------------------------------------------------- Gaussian torsion terms

#[derive(Debug, Clone)]
struct GaussianTerm {
    at: [usize; 4],
    energies: Vec<f64>,
    gradients: Vec<f64>,
    scaling: f64,
}

/// Torsion terms whose energy and gradient are looked up in tables.
#[derive(Debug, Clone)]
pub struct GaussianTorsionAngleContribs {
    dim: usize,
    terms: Vec<GaussianTerm>,
}

impl GaussianTorsionAngleContribs {
    pub fn new(dim: usize) -> Self {
        GaussianTorsionAngleContribs {
            dim,
            terms: Vec::new(),
        }
    }

    /// Adds a term; the tables sample `[0, pi]` and need at least 2 values.
    pub fn add_contrib(
        &mut self,
        at: [usize; 4],
        energies: Vec<f64>,
        gradients: Vec<f64>,
        scaling: f64,
    ) {
        assert!(all_different(&at), "degenerate points");
        assert!(
            energies.len() >= 2,
            "Energy lookup table needs to have at least 2 elements."
        );
        assert!(
            gradients.len() >= 2,
            "Gradient lookup table needs to have at least 2 elements."
        );
        self.terms.push(GaussianTerm {
            at,
            energies,
            gradients,
            scaling,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }
}

impl ForceFieldContrib for GaussianTorsionAngleContribs {
    fn energy(&self, pos: &[f64]) -> f64 {
        let mut accum = 0.0;
        for t in &self.terms {
            let p = t.at.map(|a| point3(pos, a));
            let cos_phi = calc_torsion_cos_phi(&p[0], &p[1], &p[2], &p[3]);
            accum += interpolate(&t.energies, cos_phi.acos()) * t.scaling;
        }
        accum
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        debug_assert_eq!(self.dim, 3);
        for t in &self.terms {
            let p = t.at.map(|a| point3(pos, a));
            let r = [p[0] - p[1], p[2] - p[1], p[1] - p[2], p[3] - p[2]];
            let mut tv = [r[0].cross_product(&r[1]), r[2].cross_product(&r[3])];
            let d = [tv[0].length(), tv[1].length()];
            if is_double_zero(d[0]) || is_double_zero(d[1]) {
                continue;
            }
            tv[0] /= d[0];
            tv[1] /= d[1];
            let cos_phi = tv[0].dot_product(&tv[1]).clamp(-1.0, 1.0);
            let phi = cos_phi.acos();
            let sin_phi = (1.0 - cos_phi * cos_phi).max(0.0).sqrt();
            let safe_sin_phi = sin_phi.max(1e-8);
            let de_dphi = interpolate(&t.gradients, phi) * t.scaling;
            let sin_term = -de_dphi / safe_sin_phi;
            add_torsion_grad(grad, t.at, &r, &tv, &d, cos_phi, sin_term);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ------------------------------------------------ cosine torsion terms

/// `calcTorsionEnergy`: six cosine terms with force constants and signs.
pub fn calc_cosine_torsion_energy(
    force_constants: &[f64; 6],
    signs: &[i32; 6],
    cos_phi: f64,
) -> f64 {
    let cos_phi2 = cos_phi * cos_phi;
    let cos_phi3 = cos_phi * cos_phi2;
    let cos_phi4 = cos_phi * cos_phi3;
    let cos_phi5 = cos_phi * cos_phi4;
    let cos_phi6 = cos_phi * cos_phi5;

    let cos2 = 2.0 * cos_phi2 - 1.0;
    let cos3 = 4.0 * cos_phi3 - 3.0 * cos_phi;
    let cos4 = 8.0 * cos_phi4 - 8.0 * cos_phi2 + 1.0;
    let cos5 = 16.0 * cos_phi5 - 20.0 * cos_phi3 + 5.0 * cos_phi;
    let cos6 = 32.0 * cos_phi6 - 48.0 * cos_phi4 + 18.0 * cos_phi2 - 1.0;
    let s = signs.map(f64::from);
    let k = force_constants;
    k[0] * (1.0 + s[0] * cos_phi)
        + k[1] * (1.0 + s[1] * cos2)
        + k[2] * (1.0 + s[2] * cos3)
        + k[3] * (1.0 + s[3] * cos4)
        + k[4] * (1.0 + s[4] * cos5)
        + k[5] * (1.0 + s[5] * cos6)
}

#[derive(Debug, Clone)]
struct CosineTerm {
    at: [usize; 4],
    force_constants: [f64; 6],
    signs: [i32; 6],
}

/// Experimental torsion terms: sums of up to six cosines.
#[derive(Debug, Clone, Default)]
pub struct TorsionAngleContribs {
    terms: Vec<CosineTerm>,
}

impl TorsionAngleContribs {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_contrib(&mut self, at: [usize; 4], force_constants: [f64; 6], signs: [i32; 6]) {
        assert!(all_different(&at), "degenerate points");
        self.terms.push(CosineTerm {
            at,
            force_constants,
            signs,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }
}

impl ForceFieldContrib for TorsionAngleContribs {
    fn energy(&self, pos: &[f64]) -> f64 {
        let mut accum = 0.0;
        for t in &self.terms {
            let p = t.at.map(|a| point3(pos, a));
            accum += calc_cosine_torsion_energy(
                &t.force_constants,
                &t.signs,
                calc_torsion_cos_phi(&p[0], &p[1], &p[2], &p[3]),
            );
        }
        accum
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for t in &self.terms {
            let p = t.at.map(|a| point3(pos, a));
            let r = [p[0] - p[1], p[2] - p[1], p[1] - p[2], p[3] - p[2]];
            let mut tv = [r[0].cross_product(&r[1]), r[2].cross_product(&r[3])];
            let d = [tv[0].length(), tv[1].length()];
            if is_double_zero(d[0]) || is_double_zero(d[1]) {
                continue;
            }
            tv[0] /= d[0];
            tv[1] /= d[1];
            let cos_phi = tv[0].dot_product(&tv[1]).clamp(-1.0, 1.0);
            let sin_phi_sq = 1.0 - cos_phi * cos_phi;
            let sin_phi = if sin_phi_sq > 0.0 {
                sin_phi_sq.sqrt()
            } else {
                0.0
            };
            let cos_phi2 = cos_phi * cos_phi;
            let cos_phi3 = cos_phi * cos_phi2;
            let cos_phi4 = cos_phi * cos_phi3;
            let cos_phi5 = cos_phi * cos_phi4;
            let k = &t.force_constants;
            let s = t.signs.map(f64::from);
            // dE/dPhi is independent of the cartesians. Like upstream, the
            // sixth term uses the constant and sign of the fifth.
            let de_dphi = -k[0] * s[0] * sin_phi
                - 2.0 * k[1] * s[1] * (2.0 * cos_phi * sin_phi)
                - 3.0 * k[2] * s[2] * (4.0 * cos_phi2 * sin_phi - sin_phi)
                - 4.0 * k[3] * s[3] * (8.0 * cos_phi3 * sin_phi - 4.0 * cos_phi * sin_phi)
                - 5.0
                    * k[4]
                    * s[4]
                    * (16.0 * cos_phi4 * sin_phi - 12.0 * cos_phi2 * sin_phi + sin_phi)
                - 6.0
                    * k[4]
                    * s[4]
                    * (32.0 * cos_phi5 * sin_phi - 32.0 * cos_phi3 * sin_phi + 6.0 * sin_phi);
            // the tolerance is the recommendation of Niketic and Rasmussen
            let sin_term = -de_dphi
                * if is_double_zero(sin_phi) {
                    1.0 / cos_phi
                } else {
                    1.0 / sin_phi
                };
            add_torsion_grad(grad, t.at, &r, &tv, &d, cos_phi, sin_term);
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

// ------------------------------------------------------ planarity terms

fn calc_oop_chi(i: &Point3D, j: &Point3D, k: &Point3D, l: &Point3D) -> f64 {
    let mut r_ji = *i - *j;
    let mut r_jk = *k - *j;
    let mut r_jl = *l - *j;
    r_ji /= r_ji.length();
    r_jk /= r_jk.length();
    r_jl /= r_jl.length();
    let mut n = r_ji.cross_product(&r_jk);
    n /= n.length();
    let sin_chi = n.dot_product(&r_jl).clamp(-1.0, 1.0);
    sin_chi.asin() * 180.0 / PI
}

#[derive(Debug, Clone)]
struct PlanarityTerm {
    at: [usize; 4],
    force_constant: f64,
}

/// Penalizes the out-of-plane angle of the fourth atom of each term.
#[derive(Debug, Clone, Default)]
pub struct PlanarityContribs {
    terms: Vec<PlanarityTerm>,
}

impl PlanarityContribs {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_contrib(&mut self, at: [usize; 4], force_constant: f64) {
        self.terms.push(PlanarityTerm { at, force_constant });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }
}

impl ForceFieldContrib for PlanarityContribs {
    fn energy(&self, pos: &[f64]) -> f64 {
        let mut total = 0.0;
        for t in &self.terms {
            let p = t.at.map(|a| point3(pos, a));
            let chi = calc_oop_chi(&p[0], &p[1], &p[2], &p[3]);
            total += 0.5 * t.force_constant * chi * chi;
        }
        total
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for t in &self.terms {
            let p = t.at.map(|a| point3(pos, a));
            let mut r_ji = p[0] - p[1];
            let mut r_jk = p[2] - p[1];
            let mut r_jl = p[3] - p[1];
            let (d_ji, d_jk, d_jl) = (r_ji.length(), r_jk.length(), r_jl.length());
            if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
                continue;
            }
            r_ji /= d_ji;
            r_jk /= d_jk;
            r_jl /= d_jl;

            let mut n = (-r_ji).cross_product(&r_jk);
            n /= n.length();
            let sin_chi = r_jl.dot_product(&n).clamp(-1.0, 1.0);
            let cos_chi_sq = 1.0 - sin_chi * sin_chi;
            let cos_chi = (if cos_chi_sq > 0.0 {
                cos_chi_sq.sqrt()
            } else {
                0.0
            })
            .max(1.0e-8);
            let chi = sin_chi.asin() * 180.0 / PI;
            let cos_theta = r_ji.dot_product(&r_jk).clamp(-1.0, 1.0);
            let sin_theta_sq = (1.0 - cos_theta * cos_theta).max(1.0e-8);
            let sin_theta = (if sin_theta_sq > 0.0 {
                sin_theta_sq.sqrt()
            } else {
                0.0
            })
            .max(1.0e-8);

            let de_dchi = t.force_constant * chi * 180.0 / PI;
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
            for i in 0..3 {
                grad[3 * t.at[0] + i] += de_dchi * tg1[i];
                grad[3 * t.at[1] + i] += -de_dchi * (tg1[i] + tg3[i] + tg4[i]);
                grad[3 * t.at[2] + i] += de_dchi * tg3[i];
                grad[3 * t.at[3] + i] += de_dchi * tg4[i];
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

    fn numeric_grad(c: &dyn ForceFieldContrib, pos: &[f64]) -> Vec<f64> {
        let h = 1e-6;
        (0..pos.len())
            .map(|i| {
                let mut p = pos.to_vec();
                p[i] += h;
                let up = c.energy(&p);
                p[i] -= 2.0 * h;
                (up - c.energy(&p)) / (2.0 * h)
            })
            .collect()
    }

    fn check(c: &dyn ForceFieldContrib, pos: &[f64], tol: f64) {
        let mut g = vec![0.0; pos.len()];
        c.add_grad(pos, &mut g);
        for (a, n) in g.iter().zip(numeric_grad(c, pos)) {
            assert!((a - n).abs() < tol, "analytic {a} numeric {n}");
        }
    }

    const POS: [f64; 12] = [0.1, 1.2, 0.3, 0.0, 0.0, 0.0, 1.5, 0.1, -0.2, 2.1, 1.1, 0.7];

    #[test]
    fn cosine_torsion_gradient_matches_energy() {
        let mut c = TorsionAngleContribs::new();
        // the fifth and sixth terms share a constant in upstream's gradient,
        // so only the first four terms can be compared with finite differences
        c.add_contrib(
            [0, 1, 2, 3],
            [0.5, 1.2, 0.7, 0.3, 0.0, 0.0],
            [1, -1, 1, -1, 1, 1],
        );
        check(&c, &POS, 1e-5);
    }

    #[test]
    fn gaussian_gradient_matches_tabulated_energy() {
        let (h, p, w) = ([1.0, 0.5], [0.5, 2.5], [0.3, 0.6]);
        let n = 20 * LOOKUP_GRID_SIZE;
        let phis: Vec<f64> = (0..n).map(|i| PI * i as f64 / (n - 1) as f64).collect();
        let energies: Vec<f64> = phis
            .iter()
            .map(|&x| gaussian_torsion_energy(&h, &p, &w, x))
            .collect();
        let grads: Vec<f64> = phis
            .iter()
            .map(|&x| gaussian_torsion_de_dphi(&h, &p, &w, x))
            .collect();
        let mut c = GaussianTorsionAngleContribs::new(3);
        c.add_contrib([0, 1, 2, 3], energies, grads, 1.0);
        check(&c, &POS, 2e-2);
    }

    #[test]
    fn planarity_gradient_matches_energy() {
        let mut c = PlanarityContribs::new();
        c.add_contrib([0, 1, 2, 3], 0.01);
        check(&c, &POS, 1e-5);
    }
}
