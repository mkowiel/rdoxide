//! Least-squares alignment of two point sets (quaternion method).
//!
//! Ported from RDKit `Code/Numerics/Alignment/AlignPoints.cpp` (2026.09.1).
//! The 4x4 eigenproblem is solved with the Jacobi routine of the
//! `quatfit.c` program by David J. Heisterberg (Ohio Supercomputer Center,
//! 1990), as upstream does.
//!
//! Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::geometry::{Point3D, Transform3D};
use crate::{RDKitError, Result};

const TOLERANCE: f64 = 1.0e-6;

type Mat4 = [[f64; 4]; 4];

fn weighted_sum_of_points(points: &[Point3D], weights: Option<&[f64]>) -> Point3D {
    let mut res = Point3D::default();
    for (i, p) in points.iter().enumerate() {
        res += weights.map_or(*p, |w| *p * w[i]);
    }
    res
}

fn weighted_sum_of_len_sq(points: &[Point3D], weights: Option<&[f64]>) -> f64 {
    points
        .iter()
        .enumerate()
        .map(|(i, p)| p.length_sq() * weights.map_or(1.0, |w| w[i]))
        .sum()
}

/// `covMat[a][b] = sum(w * probe[a] * ref[b])`
fn covariance(refs: &[Point3D], probes: &[Point3D], weights: Option<&[f64]>) -> [[f64; 3]; 3] {
    let mut cov = [[0.0; 3]; 3];
    for (i, (r, p)) in refs.iter().zip(probes).enumerate() {
        let w = weights.map_or(1.0, |w| w[i]);
        let (pv, rv) = ([p.x, p.y, p.z], [r.x, r.y, r.z]);
        for a in 0..3 {
            for b in 0..3 {
                cov[a][b] += w * pv[a] * rv[b];
            }
        }
    }
    cov
}

/// The 4x4 matrix whose lowest eigenvector is the optimal rotation quaternion.
fn covariance_to_quad(
    cov: &[[f64; 3]; 3],
    rpt_sum: &Point3D,
    ppt_sum: &Point3D,
    wts_sum: f64,
) -> Mat4 {
    let r = [rpt_sum.x, rpt_sum.y, rpt_sum.z];
    let p = [ppt_sum.x, ppt_sum.y, ppt_sum.z];
    // PaRb: centered covariance
    let m = |a: usize, b: usize| cov[a][b] - p[a] / wts_sum * r[b];
    let (pxrx, pxry, pxrz) = (m(0, 0), m(0, 1), m(0, 2));
    let (pyrx, pyry, pyrz) = (m(1, 0), m(1, 1), m(1, 2));
    let (pzrx, pzry, pzrz) = (m(2, 0), m(2, 1), m(2, 2));

    let mut quad = [[0.0; 4]; 4];
    quad[0][0] = -2.0 * (pxrx + pyry + pzrz);
    quad[1][1] = -2.0 * (pxrx - pyry - pzrz);
    quad[2][2] = -2.0 * (pyry - pzrz - pxrx);
    quad[3][3] = -2.0 * (pzrz - pxrx - pyry);

    let mut sym = |i: usize, j: usize, v: f64| {
        quad[i][j] = v;
        quad[j][i] = v;
    };
    sym(0, 1, 2.0 * (pyrz - pzry));
    sym(0, 2, 2.0 * (pzrx - pxrz));
    sym(0, 3, 2.0 * (pxry - pyrx));
    sym(1, 2, -2.0 * (pxry + pyrx));
    sym(1, 3, -2.0 * (pzrx + pxrz));
    sym(2, 3, -2.0 * (pyrz + pzry));
    quad
}

/// Jacobi diagonalization of the symmetric 4x4 `quad` (upper triangle used):
/// eigenvalues ascending with the matching eigenvector columns.
fn jacobi(mut quad: Mat4, max_iter: u32) -> ([f64; 4], Mat4) {
    let mut vecs = [[0.0; 4]; 4];
    let mut vals = [0.0; 4];
    for j in 0..4 {
        vecs[j][j] = 1.0;
        vals[j] = quad[j][j];
    }

    for _ in 0..max_iter {
        let diag_norm: f64 = vals.iter().map(|v| v.abs()).sum();
        let mut off_diag_norm = 0.0;
        for j in 0..4 {
            for i in 0..j {
                off_diag_norm += quad[i][j].abs();
            }
        }
        if diag_norm.abs() > 1.0e-16 && off_diag_norm / diag_norm <= TOLERANCE {
            break;
        }
        for j in 1..4 {
            for i in 0..j {
                let b = quad[i][j];
                if b.abs() <= 0.0 {
                    continue;
                }
                let dma = vals[j] - vals[i];
                let t = if dma.abs() + b.abs() <= dma.abs() {
                    b / dma
                } else {
                    let q = 0.5 * dma / b;
                    let t = 1.0 / (q.abs() + (1.0 + q * q).sqrt());
                    if q < 0.0 { -t } else { t }
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                quad[i][j] = 0.0;
                for k in 0..i {
                    let atemp = c * quad[k][i] - s * quad[k][j];
                    quad[k][j] = s * quad[k][i] + c * quad[k][j];
                    quad[k][i] = atemp;
                }
                for k in i + 1..j {
                    let atemp = c * quad[i][k] - s * quad[k][j];
                    quad[k][j] = s * quad[i][k] + c * quad[k][j];
                    quad[i][k] = atemp;
                }
                for k in j + 1..4 {
                    let atemp = c * quad[i][k] - s * quad[j][k];
                    quad[j][k] = s * quad[i][k] + c * quad[j][k];
                    quad[i][k] = atemp;
                }
                for row in &mut vecs {
                    let vtemp = c * row[i] - s * row[j];
                    row[j] = s * row[i] + c * row[j];
                    row[i] = vtemp;
                }
                let dtemp = c * c * vals[i] + s * s * vals[j] - 2.0 * c * s * b;
                vals[j] = s * s * vals[i] + c * c * vals[j] + 2.0 * c * s * b;
                vals[i] = dtemp;
            }
        }
    }

    // selection sort of the eigenvalues, ascending
    for j in 0..3 {
        let mut k = j;
        let mut dtemp = vals[k];
        for i in j + 1..4 {
            if vals[i] < dtemp {
                k = i;
                dtemp = vals[k];
            }
        }
        if k > j {
            vals[k] = vals[j];
            vals[j] = dtemp;
            for row in &mut vecs {
                row.swap(j, k);
            }
        }
    }
    (vals, vecs)
}

/// `AlignPoints`: the transformation taking `probe_points` onto `ref_points`
/// with the least weighted sum of squared residuals, and that sum. With
/// `reflect` the probe is mirrored first.
pub fn align_points(
    ref_points: &[Point3D],
    probe_points: &[Point3D],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<(f64, Transform3D)> {
    let npt = ref_points.len();
    if npt != probe_points.len() {
        return Err(RDKitError::precondition("Mismatch in number of points"));
    }
    if let Some(w) = weights {
        if w.len() != npt {
            return Err(RDKitError::precondition("Mismatch in number of points"));
        }
    }
    let wts_sum = match weights {
        Some(w) => {
            if w.iter().any(|&x| x <= 0.0) {
                return Err(RDKitError::invariant(
                    "Negative weight specified for a point",
                ));
            }
            w.iter().sum()
        }
        None => npt as f64,
    };

    let mut rpt_sum = weighted_sum_of_points(ref_points, weights);
    let ppt_sum = weighted_sum_of_points(probe_points, weights);
    let rpt_sum_len_sq = weighted_sum_of_len_sq(ref_points, weights);
    let ppt_sum_len_sq = weighted_sum_of_len_sq(probe_points, weights);

    let mut cov = covariance(ref_points, probe_points, weights);
    if reflect {
        rpt_sum *= -1.0;
        for row in &mut cov {
            for v in row {
                *v = -*v;
            }
        }
    }

    let quad = covariance_to_quad(&cov, &rpt_sum, &ppt_sum, wts_sum);
    let (vals, vecs) = jacobi(quad, max_iterations);
    let quaternion = [vecs[0][0], vecs[1][0], vecs[2][0], vecs[3][0]];

    let mut trans = Transform3D::identity();
    trans.set_rotation_from_quaternion(&quaternion);
    if reflect {
        trans.reflect();
    }
    let mut ssr = vals[0] - (ppt_sum.length_sq() + rpt_sum.length_sq()) / wts_sum
        + rpt_sum_len_sq
        + ppt_sum_len_sq;
    if ssr < 0.0 && ssr.abs() < TOLERANCE {
        ssr = 0.0;
    }
    if reflect {
        rpt_sum *= -1.0;
    }

    let rotated = trans.transform_point(&ppt_sum);
    let mut shift = rpt_sum - rotated;
    shift /= wts_sum;
    trans.set_translation(&shift);
    Ok((ssr, trans))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square() -> Vec<Point3D> {
        vec![
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(2.0, 0.0, 0.0),
            Point3D::new(2.0, 1.0, 0.5),
            Point3D::new(0.0, 1.0, -0.5),
        ]
    }

    #[test]
    fn rigid_motion_is_undone() {
        let reference = square();
        let mut motion = Transform3D::identity();
        motion.set_rotation(0.7, &Point3D::new(0.0, 0.0, 1.0));
        motion.set_translation(&Point3D::new(1.0, -2.0, 3.0));
        let probe: Vec<_> = reference
            .iter()
            .map(|p| motion.transform_point(p))
            .collect();
        let (ssr, trans) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!(ssr < 1e-12);
        for (p, r) in probe.iter().zip(&reference) {
            assert!((trans.transform_point(p) - *r).length_sq() < 1e-12);
        }
    }

    #[test]
    fn mirror_image_needs_reflection() {
        let reference = square();
        let probe: Vec<_> = reference
            .iter()
            .map(|p| Point3D::new(-p.x, p.y, p.z))
            .collect();
        let (plain, _) = align_points(&reference, &probe, None, false, 50).unwrap();
        let (mirrored, _) = align_points(&reference, &probe, None, true, 50).unwrap();
        assert!(plain > 0.1);
        assert!(mirrored < 1e-12);
    }

    #[test]
    fn bad_input_is_rejected() {
        let pts = square();
        assert!(align_points(&pts, &pts[..3], None, false, 50).is_err());
        assert!(align_points(&pts, &pts, Some(&[1.0, 1.0, -1.0, 1.0]), false, 50).is_err());
    }
}
