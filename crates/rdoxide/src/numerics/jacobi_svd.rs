//! Singular value decomposition of real 3x3 matrices by two-sided Jacobi
//! rotations, reproducing Eigen's `JacobiSVD<MatrixXd>` (singular values in
//! descending order and the matrix `V` of right singular vectors, with Eigen's
//! choice of signs and of the basis inside degenerate subspaces).
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

use crate::rdgeneral::fp::fmuladd;

const N: usize = 3;
/// `2 * epsilon`: Eigen stops when off-diagonal entries are this small
/// relative to the largest diagonal entry.
const PRECISION: f64 = 2.0 * f64::EPSILON;

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

/// Rotates `x` and `y` in place (`apply_rotation_in_the_plane`, scalar path).
fn rotate(x: &mut f64, y: &mut f64, rot: Rotation) {
    let (xi, yi) = (*x, *y);
    *x = fmuladd(rot.c, xi, rot.s * yi);
    *y = fmuladd(-rot.s, xi, rot.c * yi);
}

fn is_identity(rot: Rotation) -> bool {
    rot.c == 1.0 && rot.s == 0.0
}

/// `applyOnTheLeft(p, q, rot)`: rotates rows `p` and `q`.
fn apply_on_the_left(m: &mut [[f64; N]; N], p: usize, q: usize, rot: Rotation) {
    if is_identity(rot) {
        return;
    }
    for col in 0..N {
        let (mut x, mut y) = (m[p][col], m[q][col]);
        rotate(&mut x, &mut y, rot);
        (m[p][col], m[q][col]) = (x, y);
    }
}

/// `applyOnTheRight(p, q, rot)`: rotates columns `p` and `q`.
fn apply_on_the_right(m: &mut [[f64; N]; N], p: usize, q: usize, rot: Rotation) {
    let rot = rot.transpose();
    if is_identity(rot) {
        return;
    }
    for row in m.iter_mut() {
        let (mut x, mut y) = (row[p], row[q]);
        rotate(&mut x, &mut y, rot);
        (row[p], row[q]) = (x, y);
    }
}

/// `real_2x2_jacobi_svd`: the left and right rotations that diagonalize the
/// 2x2 block `(p, q)` of `a`.
fn real_2x2_jacobi_svd(a: &[[f64; N]; N], p: usize, q: usize) -> (Rotation, Rotation) {
    let mut m = [[a[p][p], a[p][q]], [a[q][p], a[q][q]]];
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
    if !is_identity(rot1) {
        for col in 0..2 {
            let (mut x, mut y) = (m[0][col], m[1][col]);
            rotate(&mut x, &mut y, rot1);
            (m[0][col], m[1][col]) = (x, y);
        }
    }
    let right = Rotation::make_jacobi(m[0][0], m[0][1], m[1][1]);
    (rot1.compose(right.transpose()), right)
}

/// The decomposition `A = U diag(singular_values) V^T`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Svd3 {
    /// descending
    pub singular_values: [f64; N],
    /// `v[row][col]`; column `i` belongs to `singular_values[i]`
    pub v: [[f64; N]; N],
}

/// `JacobiSVD<MatrixXd>(m, ComputeThinV)` of a 3x3 matrix given as `m[row][col]`.
pub fn jacobi_svd_3x3(m: &[[f64; N]; N]) -> Svd3 {
    let mut scale = 0.0_f64;
    for row in m {
        for &x in row {
            if x.is_nan() {
                scale = f64::NAN;
            } else if !scale.is_nan() {
                scale = scale.max(x.abs());
            }
        }
    }
    if !scale.is_finite() {
        // Eigen reports invalid input and leaves the result unset
        return Svd3 {
            singular_values: [f64::NAN; N],
            v: [[f64::NAN; N]; N],
        };
    }
    if scale == 0.0 {
        scale = 1.0;
    }
    let mut work = *m;
    for row in &mut work {
        for x in row.iter_mut() {
            *x /= scale;
        }
    }
    let mut v = [[0.0; N]; N];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    let mut max_diag = (0..N).fold(0.0_f64, |acc, i| acc.max(work[i][i].abs()));
    let mut finished = false;
    while !finished {
        finished = true;
        for p in 1..N {
            for q in 0..p {
                let threshold = f64::MIN_POSITIVE.max(PRECISION * max_diag);
                if work[p][q].abs() > threshold || work[q][p].abs() > threshold {
                    finished = false;
                    let (left, right) = real_2x2_jacobi_svd(&work, p, q);
                    apply_on_the_left(&mut work, p, q, left);
                    apply_on_the_right(&mut work, p, q, right);
                    apply_on_the_right(&mut v, p, q, right);
                    max_diag = max_diag.max(work[p][p].abs().max(work[q][q].abs()));
                }
            }
        }
    }

    let mut values = [0.0; N];
    for i in 0..N {
        values[i] = work[i][i].abs() * scale;
    }
    // descending order, the first maximum wins ties
    for i in 0..N {
        let mut pos = i;
        for j in i + 1..N {
            if values[j] > values[pos] {
                pos = j;
            }
        }
        if values[pos] == 0.0 {
            break;
        }
        if pos != i {
            values.swap(i, pos);
            for row in v.iter_mut() {
                row.swap(i, pos);
            }
        }
    }
    Svd3 {
        singular_values: values,
        v,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_matrix_sorts_its_values() {
        let svd = jacobi_svd_3x3(&[[1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 2.0]]);
        assert_eq!(svd.singular_values, [3.0, 2.0, 1.0]);
        assert_eq!(svd.v, [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    }

    #[test]
    fn symmetric_matrix_is_reconstructed() {
        let a = [[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]];
        let svd = jacobi_svd_3x3(&a);
        for i in 0..N {
            for j in 0..N {
                let r: f64 = (0..N)
                    .map(|k| svd.v[i][k] * svd.singular_values[k] * svd.v[j][k])
                    .sum();
                assert!((r - a[i][j]).abs() < 1e-12, "({i},{j}): {r} vs {}", a[i][j]);
            }
        }
        assert!(svd.singular_values[0] >= svd.singular_values[1]);
    }

    #[test]
    fn zero_matrix_has_zero_singular_values() {
        let svd = jacobi_svd_3x3(&[[0.0; N]; N]);
        assert_eq!(svd.singular_values, [0.0; N]);
    }
}
