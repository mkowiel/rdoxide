//! The terms of the distance-geometry error function (upstream
//! `DistGeom/{DistViolationContribs,ChiralViolationContribs,FourthDimContribs}`).

use super::ChiralSet;
use crate::forcefield::ForceFieldContrib;
use crate::geometry::Point3D;

fn point(pos: &[f64], idx: usize, dim: usize) -> Point3D {
    Point3D::new(pos[idx * dim], pos[idx * dim + 1], pos[idx * dim + 2])
}

#[derive(Debug, Clone)]
struct DistTerm {
    idx1: usize,
    idx2: usize,
    ub2: f64,
    lb2: f64,
    weight: f64,
}

/// Penalizes violations of the distance bounds between pairs of points.
#[derive(Debug, Clone)]
pub struct DistViolationContribs {
    dim: usize,
    terms: Vec<DistTerm>,
}

impl DistViolationContribs {
    pub fn new(dim: usize) -> Self {
        DistViolationContribs {
            dim,
            terms: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, idx1: usize, idx2: usize, ub: f64, lb: f64, weight: f64) {
        self.terms.push(DistTerm {
            idx1,
            idx2,
            ub2: ub * ub,
            lb2: lb * lb,
            weight,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }

    fn distance2(&self, pos: &[f64], a: usize, b: usize) -> f64 {
        (0..self.dim)
            .map(|i| {
                let d = pos[self.dim * a + i] - pos[self.dim * b + i];
                d * d
            })
            .sum()
    }
}

impl ForceFieldContrib for DistViolationContribs {
    fn energy(&self, pos: &[f64]) -> f64 {
        let mut accum = 0.0;
        for c in &self.terms {
            let d2 = self.distance2(pos, c.idx1, c.idx2);
            let mut val = 0.0;
            if d2 > c.ub2 {
                val = d2 / c.ub2 - 1.0;
            } else if d2 < c.lb2 {
                val = (2.0 * c.lb2) / (c.lb2 + d2) - 1.0;
            }
            if val > 0.0 {
                accum += c.weight * val * val;
            }
        }
        accum
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let dim = self.dim;
        for c in &self.terms {
            let d2 = self.distance2(pos, c.idx1, c.idx2);
            let (d, pre_factor);
            if d2 > c.ub2 {
                d = d2.sqrt();
                pre_factor = 4. * ((d * d) / c.ub2 - 1.0) * (d / c.ub2);
            } else if d2 < c.lb2 {
                d = d2.sqrt();
                let l2d2 = d2 + c.lb2;
                pre_factor = 8. * c.lb2 * d * (1. - 2. * c.lb2 / l2d2) / (l2d2 * l2d2);
            } else {
                continue;
            }
            for i in 0..dim {
                let (p1, p2) = (dim * c.idx1 + i, dim * c.idx2 + i);
                let d_grad = if d > 0.0 {
                    c.weight * pre_factor * (pos[p1] - pos[p2]) / d
                } else {
                    // upstream notes that this is likely not right
                    c.weight * pre_factor * (pos[p1] - pos[p2])
                };
                grad[p1] += d_grad;
                grad[p2] -= d_grad;
            }
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

/// The chiral volume of four points, from the first three dimensions.
pub fn calc_chiral_volume(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    pos: &[f64],
    dim: usize,
) -> f64 {
    let p4 = point(pos, idx4, dim);
    let v1 = point(pos, idx1, dim) - p4;
    let v2 = point(pos, idx2, dim) - p4;
    let v3 = point(pos, idx3, dim) - p4;
    v1.dot_product(&v2.cross_product(&v3))
}

#[derive(Debug, Clone)]
struct ChiralTerm {
    idx: [usize; 4],
    vol_upper: f64,
    vol_lower: f64,
    weight: f64,
}

/// Penalizes chiral volumes outside their bounds.
#[derive(Debug, Clone)]
pub struct ChiralViolationContribs {
    dim: usize,
    terms: Vec<ChiralTerm>,
}

impl ChiralViolationContribs {
    pub fn new(dim: usize) -> Self {
        ChiralViolationContribs {
            dim,
            terms: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, cset: &ChiralSet, weight: f64) {
        self.terms.push(ChiralTerm {
            idx: [cset.idx1, cset.idx2, cset.idx3, cset.idx4],
            vol_upper: cset.volume_upper_bound,
            vol_lower: cset.volume_lower_bound,
            weight,
        });
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }
}

impl ForceFieldContrib for ChiralViolationContribs {
    fn energy(&self, pos: &[f64]) -> f64 {
        let mut res = 0.0;
        for c in &self.terms {
            let vol = calc_chiral_volume(c.idx[0], c.idx[1], c.idx[2], c.idx[3], pos, self.dim);
            if vol < c.vol_lower {
                res += c.weight * (vol - c.vol_lower) * (vol - c.vol_lower);
            } else if vol > c.vol_upper {
                res += c.weight * (vol - c.vol_upper) * (vol - c.vol_upper);
            }
        }
        res
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let dim = self.dim;
        for c in &self.terms {
            let [i1, i2, i3, i4] = c.idx;
            let p4 = point(pos, i4, dim);
            let v1 = point(pos, i1, dim) - p4;
            let v2 = point(pos, i2, dim) - p4;
            let v3 = point(pos, i3, dim) - p4;
            let vol = v1.dot_product(&v2.cross_product(&v3));
            let pre_factor = if vol < c.vol_lower {
                c.weight * (vol - c.vol_lower)
            } else if vol > c.vol_upper {
                c.weight * (vol - c.vol_upper)
            } else {
                continue;
            };
            // twelve variables: four points in three dimensions
            grad[dim * i1] += pre_factor * (v2.y * v3.z - v3.y * v2.z);
            grad[dim * i1 + 1] += pre_factor * (v3.x * v2.z - v2.x * v3.z);
            grad[dim * i1 + 2] += pre_factor * (v2.x * v3.y - v3.x * v2.y);

            grad[dim * i2] += pre_factor * (v3.y * v1.z - v3.z * v1.y);
            grad[dim * i2 + 1] += pre_factor * (v3.z * v1.x - v3.x * v1.z);
            grad[dim * i2 + 2] += pre_factor * (v3.x * v1.y - v3.y * v1.x);

            grad[dim * i3] += pre_factor * (v2.z * v1.y - v2.y * v1.z);
            grad[dim * i3 + 1] += pre_factor * (v2.x * v1.z - v2.z * v1.x);
            grad[dim * i3 + 2] += pre_factor * (v2.y * v1.x - v2.x * v1.y);

            let p = |i: usize, k: usize| pos[i * dim + k];
            grad[dim * i4] += pre_factor
                * (p(i1, 2) * (p(i2, 1) - p(i3, 1))
                    + p(i2, 2) * (p(i3, 1) - p(i1, 1))
                    + p(i3, 2) * (p(i1, 1) - p(i2, 1)));
            grad[dim * i4 + 1] += pre_factor
                * (p(i1, 0) * (p(i2, 2) - p(i3, 2))
                    + p(i2, 0) * (p(i3, 2) - p(i1, 2))
                    + p(i3, 0) * (p(i1, 2) - p(i2, 2)));
            grad[dim * i4 + 2] += pre_factor
                * (p(i1, 1) * (p(i2, 0) - p(i3, 0))
                    + p(i2, 1) * (p(i3, 0) - p(i1, 0))
                    + p(i3, 1) * (p(i1, 0) - p(i2, 0)));
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

/// Penalizes the fourth coordinate, to move from 4D to 3D.
#[derive(Debug, Clone, Default)]
pub struct FourthDimContribs {
    terms: Vec<(usize, f64)>,
}

impl FourthDimContribs {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_contrib(&mut self, idx: usize, weight: f64) {
        self.terms.push((idx, weight));
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }
}

impl ForceFieldContrib for FourthDimContribs {
    fn energy(&self, pos: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|&(idx, w)| w * pos[idx * 4 + 3] * pos[idx * 4 + 3])
            .sum()
    }

    fn add_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for &(idx, w) in &self.terms {
            grad[idx * 4 + 3] += w * pos[idx * 4 + 3];
        }
    }

    fn box_clone(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const POS4: [f64; 20] = [
        0.1, 0.2, 0.0, 0.3, 1.5, 0.3, 0.2, -0.2, 2.1, 1.6, -0.1, 0.1, 3.4, 1.9, 0.9, 0.4, 1.0,
        -1.1, 0.7, -0.3,
    ];

    /// `scale` is the ratio between the numeric derivative of the energy and
    /// the analytic gradient: upstream's chiral and fourth-dimension
    /// gradients leave out the factor 2 of the squared terms.
    fn check_gradient(c: &dyn ForceFieldContrib, pos: &[f64], scale: f64) {
        let mut grad = vec![0.0; pos.len()];
        c.add_grad(pos, &mut grad);
        for i in 0..pos.len() {
            let h = 1e-6;
            let (mut plus, mut minus) = (pos.to_vec(), pos.to_vec());
            plus[i] += h;
            minus[i] -= h;
            let numeric = (c.energy(&plus) - c.energy(&minus)) / (2.0 * h);
            assert!(
                (numeric - scale * grad[i]).abs() < 1e-5 * (1.0 + numeric.abs()),
                "component {i}: analytic {} numeric {numeric}",
                grad[i]
            );
        }
    }

    #[test]
    fn distance_violations() {
        let mut c = DistViolationContribs::new(4);
        // violated from above, violated from below, and satisfied
        c.add_contrib(0, 1, 1.0, 0.5, 1.0);
        c.add_contrib(1, 2, 5.0, 4.0, 2.0);
        c.add_contrib(2, 3, 10.0, 0.1, 1.0);
        assert!(c.energy(&POS4) > 0.0);
        check_gradient(&c, &POS4, 1.0);
    }

    #[test]
    fn chiral_violations() {
        let cs = ChiralSet::new([0, 1, 2, 3, 4], 2.0, 3.0, 0).unwrap();
        let mut c = ChiralViolationContribs::new(4);
        c.add_contrib(&cs, 1.5);
        assert!(c.energy(&POS4) > 0.0);
        check_gradient(&c, &POS4, 2.0);
        assert!(ChiralSet::new([0, 1, 2, 3, 4], 3.0, 2.0, 0).is_err());
    }

    #[test]
    fn fourth_dimension() {
        let mut c = FourthDimContribs::new();
        c.add_contrib(0, 0.1);
        c.add_contrib(3, 0.2);
        check_gradient(&c, &POS4, 2.0);
    }
}
