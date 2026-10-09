//! Singular value decomposition of real square matrices by two-sided Jacobi
//! rotations, reproducing Eigen's `JacobiSVD<MatrixXd>` (singular values in
//! descending order and the matrices `U` and `V` of singular vectors, with
//! Eigen's choice of signs and of the basis inside degenerate subspaces).
//!
//! Products are fused the way clang contracts them for the reference build
//! (`rdgeneral::fp`); results agree with Eigen to the last digits.
//!
//! Translated from Eigen 5.0.1: `SVD/JacobiSVD.h` (`compute_impl`),
//! `misc/RealSvd2x2.h` (`real_2x2_jacobi_svd`) and `Jacobi/Jacobi.h`
//! (`makeJacobi`, rotation products, `apply_rotation_in_the_plane`).
//!
//! Copyright (C) 2009-2010 Benoit Jacob <jacob.benoit.1@gmail.com>.
//!
//! This Source Code Form is subject to the terms of the Mozilla Public
//! License v. 2.0 (see LICENSES/Eigen-MPL-2.0.txt).

use super::dense::Matrix;
use crate::rdgeneral::fp::fmuladd;

/// `2 * epsilon`: Eigen stops when off-diagonal entries are this small
/// relative to the largest diagonal entry.
const PRECISION: f64 = 2.0 * f64::EPSILON;
/// Doubles per SIMD register of the reference build (NEON).
const PACKET: usize = 2;

/// A plane rotation (`Eigen::JacobiRotation<double>`).
#[derive(Clone, Copy)]
struct Rotation {
    c: f64,
    s: f64,
}

impl Rotation {
    fn transpose(self) -> Rotation {
        Rotation {
            c: self.c,
            s: -self.s,
        }
    }

    fn is_identity(self) -> bool {
        self.c == 1.0 && self.s == 0.0
    }

    /// `operator*`: the concatenation of two rotations.
    fn compose(self, other: Rotation) -> Rotation {
        Rotation {
            c: fmuladd(self.c, other.c, -(self.s * other.s)),
            s: fmuladd(self.c, other.s, self.s * other.c),
        }
    }

    /// `makeJacobi(x, y, z)`: the rotation that diagonalizes the symmetric
    /// 2x2 matrix `[[x, y], [y, z]]`.
    fn make_jacobi(x: f64, y: f64, z: f64) -> Rotation {
        let deno = 2.0 * y.abs();
        if deno < f64::MIN_POSITIVE {
            return Rotation { c: 1.0, s: 0.0 };
        }
        let tau = (x - z) / deno;
        let w = (tau * tau + 1.0).sqrt();
        let t = if tau > 0.0 {
            1.0 / (tau + w)
        } else {
            1.0 / (tau - w)
        };
        let sign_t = if t > 0.0 { 1.0 } else { -1.0 };
        let n = 1.0 / (t * t + 1.0).sqrt();
        Rotation {
            c: n,
            s: -sign_t * (y / y.abs()) * t.abs() * n,
        }
    }
}

/// The scalar update of `apply_rotation_in_the_plane`; clang contracts it
/// into fused multiply-adds.
fn rotate_scalar(x: f64, y: f64, rot: Rotation) -> (f64, f64) {
    (fmuladd(rot.c, x, rot.s * y), fmuladd(-rot.s, x, rot.c * y))
}

/// The update of the vectorized path: products and sums are rounded separately.
fn rotate_packet(x: f64, y: f64, rot: Rotation) -> (f64, f64) {
    (rot.c * x + rot.s * y, rot.c * y - rot.s * x)
}

/// `applyOnTheLeft(p, q, rot)`: rotates rows `p` and `q` (strided, so always
/// the scalar path).
fn apply_on_the_left(m: &mut Matrix, p: usize, q: usize, rot: Rotation) {
    if rot.is_identity() {
        return;
    }
    for col in 0..m.cols() {
        let (x, y) = rotate_scalar(m.get(p, col), m.get(q, col), rot);
        m.set(p, col, x);
        m.set(q, col, y);
    }
}

/// `applyOnTheRight(p, q, rot)`: rotates columns `p` and `q`. Columns of at
/// least two registers are contiguous and go through the vectorized path
/// between the first aligned element of column `q` and the last full register.
fn apply_on_the_right(m: &mut Matrix, p: usize, q: usize, rot: Rotation) {
    let rot = rot.transpose();
    if rot.is_identity() {
        return;
    }
    let n = m.rows();
    let vectorized = n >= 2 * PACKET;
    let start = (q * n) % PACKET;
    let end = start + (n - start.min(n)) / PACKET * PACKET;
    for i in 0..n {
        let (x, y) = (m.get(i, p), m.get(i, q));
        let (x, y) = if vectorized && (start..end).contains(&i) {
            rotate_packet(x, y, rot)
        } else {
            rotate_scalar(x, y, rot)
        };
        m.set(i, p, x);
        m.set(i, q, y);
    }
}

/// `real_2x2_jacobi_svd`: the left and right rotations that diagonalize the
/// 2x2 block `(p, q)` of `a`.
fn real_2x2_jacobi_svd(a: &Matrix, p: usize, q: usize) -> (Rotation, Rotation) {
    let mut m = [[a.get(p, p), a.get(p, q)], [a.get(q, p), a.get(q, q)]];
    let t = m[0][0] + m[1][1];
    let d = m[1][0] - m[0][1];
    let rot1 = if d.abs() < f64::MIN_POSITIVE {
        Rotation { c: 1.0, s: 0.0 }
    } else {
        let u = t / d;
        let tmp = (1.0 + u * u).sqrt();
        Rotation {
            c: u / tmp,
            s: 1.0 / tmp,
        }
    };
    if !rot1.is_identity() {
        for col in 0..2 {
            let (x, y) = rotate_scalar(m[0][col], m[1][col], rot1);
            (m[0][col], m[1][col]) = (x, y);
        }
    }
    let right = Rotation::make_jacobi(m[0][0], m[0][1], m[1][1]);
    (rot1.compose(right.transpose()), right)
}

/// The decomposition `A = U diag(singular_values) V^T`.
#[derive(Debug, Clone, PartialEq)]
pub struct Svd {
    /// descending
    pub singular_values: Vec<f64>,
    /// left singular vectors as columns (when requested)
    pub u: Option<Matrix>,
    /// right singular vectors as columns
    pub v: Matrix,
}

/// `JacobiSVD<MatrixXd>(a, ComputeThinV | ComputeThinU)` of a square matrix;
/// `U` is only accumulated when `compute_u` is set.
pub fn jacobi_svd(a: &Matrix, compute_u: bool) -> Svd {
    assert_eq!(a.rows(), a.cols(), "only square matrices are supported");
    let n = a.rows();
    let mut scale = 0.0_f64;
    for col in 0..n {
        for &x in a.column(col) {
            if x.is_nan() {
                scale = f64::NAN;
            } else if !scale.is_nan() {
                scale = scale.max(x.abs());
            }
        }
    }
    if !scale.is_finite() {
        // Eigen reports invalid input and leaves the result unset
        return Svd {
            singular_values: vec![f64::NAN; n],
            u: compute_u.then(|| Matrix::identity(n)),
            v: Matrix::identity(n),
        };
    }
    if scale == 0.0 {
        scale = 1.0;
    }
    let mut work = Matrix::from_fn(n, n, |r, c| a.get(r, c) / scale);
    let mut u = compute_u.then(|| Matrix::identity(n));
    let mut v = Matrix::identity(n);

    let mut max_diag = (0..n).fold(0.0_f64, |acc, i| acc.max(work.get(i, i).abs()));
    let mut finished = false;
    while !finished {
        finished = true;
        for p in 1..n {
            for q in 0..p {
                let threshold = f64::MIN_POSITIVE.max(PRECISION * max_diag);
                if work.get(p, q).abs() > threshold || work.get(q, p).abs() > threshold {
                    finished = false;
                    let (left, right) = real_2x2_jacobi_svd(&work, p, q);
                    apply_on_the_left(&mut work, p, q, left);
                    if let Some(u) = u.as_mut() {
                        apply_on_the_right(u, p, q, left.transpose());
                    }
                    apply_on_the_right(&mut work, p, q, right);
                    apply_on_the_right(&mut v, p, q, right);
                    max_diag = max_diag.max(work.get(p, p).abs().max(work.get(q, q).abs()));
                }
            }
        }
    }

    // the diagonal is made positive; the sign moves into U
    let mut values = vec![0.0; n];
    for i in 0..n {
        let d = work.get(i, i);
        values[i] = d.abs();
        if let Some(u) = u.as_mut()
            && d < 0.0
        {
            for x in u.column_mut(i) {
                *x = -*x;
            }
        }
        values[i] *= scale;
    }
    // descending order, the first maximum wins ties
    for i in 0..n {
        let mut pos = i;
        for j in i + 1..n {
            if values[j] > values[pos] {
                pos = j;
            }
        }
        if values[pos] == 0.0 {
            break;
        }
        if pos != i {
            values.swap(i, pos);
            swap_columns(&mut v, i, pos);
            if let Some(u) = u.as_mut() {
                swap_columns(u, i, pos);
            }
        }
    }
    Svd {
        singular_values: values,
        u,
        v,
    }
}

fn swap_columns(m: &mut Matrix, a: usize, b: usize) {
    for r in 0..m.rows() {
        let (x, y) = (m.get(r, a), m.get(r, b));
        m.set(r, a, y);
        m.set(r, b, x);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reconstruct(svd: &Svd, i: usize, j: usize) -> f64 {
        let u = svd.u.as_ref().expect("U was requested");
        (0..svd.singular_values.len())
            .map(|k| u.get(i, k) * svd.singular_values[k] * svd.v.get(j, k))
            .sum()
    }

    #[test]
    fn diagonal_matrix_sorts_its_values() {
        let a = Matrix::from_fn(3, 3, |r, c| if r == c { [1.0, 3.0, 2.0][r] } else { 0.0 });
        let svd = jacobi_svd(&a, false);
        assert_eq!(svd.singular_values, vec![3.0, 2.0, 1.0]);
        assert_eq!(svd.v.get(1, 0), 1.0);
        assert_eq!(svd.v.get(2, 1), 1.0);
        assert_eq!(svd.v.get(0, 2), 1.0);
    }

    #[test]
    fn matrices_are_reconstructed() {
        // a general 5 x 5 matrix exercises the vectorized column rotations
        let a = Matrix::from_fn(5, 5, |r, c| {
            ((r * 7 + c * 3) % 11) as f64 - 4.5 + (r == c) as u8 as f64
        });
        let svd = jacobi_svd(&a, true);
        for i in 0..5 {
            for j in 0..5 {
                assert!((reconstruct(&svd, i, j) - a.get(i, j)).abs() < 1e-12);
            }
        }
        assert!(svd.singular_values.windows(2).all(|w| w[0] >= w[1]));
    }

    #[test]
    fn zero_matrix_has_zero_singular_values() {
        let svd = jacobi_svd(&Matrix::zeros(3, 3), false);
        assert_eq!(svd.singular_values, vec![0.0; 3]);
    }

    #[test]
    fn one_by_one_matrix() {
        let svd = jacobi_svd(&Matrix::from_fn(1, 1, |_, _| -2.0), true);
        assert_eq!(svd.singular_values, vec![2.0]);
        assert_eq!(svd.u.unwrap().get(0, 0), -1.0);
    }
}
