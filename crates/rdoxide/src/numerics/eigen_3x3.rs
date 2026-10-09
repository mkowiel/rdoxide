//! Eigen decomposition of real symmetric 3x3 matrices that reproduces Eigen's
//! `SelfAdjointEigenSolver<Matrix3d>`, including the signs and order of the
//! eigenvectors (RDKit's principal axes and canonical transforms depend on
//! them).
//!
//! Products are fused the way clang contracts them for the reference build
//! (`rdgeneral::fp`); results agree with Eigen to the last digits, not always
//! bit for bit (degenerate spectra amplify the differences).
//!
//! Translated from Eigen 5.0.1, the version of the reference build:
//! `Eigenvalues/SelfAdjointEigenSolver.h` (`compute`,
//! `computeFromTridiagonal_impl`, `tridiagonal_qr_step`),
//! `Eigenvalues/Tridiagonalization.h` (the 3x3 specialization) and
//! `Jacobi/Jacobi.h` (`makeGivens`, `applyOnTheRight`).
//!
//! Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>,
//! (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>,
//! (C) 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>.
//!
//! This Source Code Form is subject to the terms of the Mozilla Public
//! License v. 2.0 (see LICENSES/Eigen-MPL-2.0.txt).

use crate::rdgeneral::fp::fmuladd;
use crate::{RDKitError, Result};

/// Iteration limit factor of Eigen's `m_maxIterations` (times the size).
const MAX_ITERATIONS: usize = 30;
const N: usize = 3;

/// Eigenvalues in ascending order and the matching eigenvectors.
#[derive(Debug, Clone, PartialEq)]
pub struct SymmetricEigen3 {
    pub values: [f64; N],
    /// `vectors[row][col]`: column `col` is the eigenvector of `values[col]`.
    pub vectors: [[f64; N]; N],
}

/// A Givens rotation (`Eigen::JacobiRotation<double>`).
struct Givens {
    c: f64,
    s: f64,
}

impl Givens {
    /// `makeGivens(p, q)` for reals.
    fn new(p: f64, q: f64) -> Self {
        if q == 0.0 {
            Givens {
                c: if p < 0.0 { -1.0 } else { 1.0 },
                s: 0.0,
            }
        } else if p == 0.0 {
            Givens {
                c: 0.0,
                s: if q < 0.0 { 1.0 } else { -1.0 },
            }
        } else if p.abs() > q.abs() {
            let t = q / p;
            let mut u = (1.0 + t * t).sqrt();
            if p < 0.0 {
                u = -u;
            }
            let c = 1.0 / u;
            Givens { c, s: -t * c }
        } else {
            let t = p / q;
            let mut u = (1.0 + t * t).sqrt();
            if q < 0.0 {
                u = -u;
            }
            let s = -1.0 / u;
            Givens { c: -t * s, s }
        }
    }
}

/// Reduction of the lower triangle of `mat` to tridiagonal form: the diagonal,
/// the sub-diagonal and the orthogonal matrix `Q` with `mat = Q T Q^T`.
fn tridiagonalize(mat: &[[f64; N]; N]) -> ([f64; N], [f64; N - 1], [[f64; N]; N]) {
    let mut diag = [mat[0][0], 0.0, 0.0];
    let v1_norm2 = mat[2][0] * mat[2][0];
    if v1_norm2 <= f64::MIN_POSITIVE {
        diag[1] = mat[1][1];
        diag[2] = mat[2][2];
        let q = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        return (diag, [mat[1][0], mat[2][1]], q);
    }
    let beta = (mat[1][0] * mat[1][0] + v1_norm2).sqrt();
    let inv_beta = 1.0 / beta;
    let m01 = mat[1][0] * inv_beta;
    let m02 = mat[2][0] * inv_beta;
    let q = fmuladd(2.0 * m01, mat[2][1], m02 * (mat[2][2] - mat[1][1]));
    diag[1] = fmuladd(m02, q, mat[1][1]);
    diag[2] = fmuladd(-m02, q, mat[2][2]);
    let sub = [beta, fmuladd(-m01, q, mat[2][1])];
    (
        diag,
        sub,
        [[1.0, 0.0, 0.0], [0.0, m01, m02], [0.0, m02, -m01]],
    )
}

/// One implicit symmetric QR step with Wilkinson shift on the block
/// `start..=end` (Golub and Van Loan, algorithm 8.3.2).
fn tridiagonal_qr_step(
    diag: &mut [f64; N],
    sub: &mut [f64; N - 1],
    start: usize,
    end: usize,
    q: &mut [[f64; N]; N],
) {
    // Wilkinson shift
    let td = (diag[end - 1] - diag[end]) * 0.5;
    let e = sub[end - 1];
    let mut mu = diag[end];
    if td == 0.0 {
        mu -= e.abs();
    } else if e != 0.0 {
        let e2 = e * e;
        let h = td.hypot(e);
        let signed_h = if td > 0.0 { h } else { -h };
        if e2 == 0.0 {
            mu -= e / ((td + signed_h) / e);
        } else {
            mu -= e2 / (td + signed_h);
        }
    }

    let mut x = diag[start] - mu;
    let mut z = sub[start];
    let mut k = start;
    while k < end && z != 0.0 {
        let rot = Givens::new(x, z);
        let (c, s) = (rot.c, rot.s);
        // T = G' T G
        let sdk = fmuladd(s, diag[k], c * sub[k]);
        let dkp1 = fmuladd(s, sub[k], c * diag[k + 1]);
        diag[k] = fmuladd(
            c,
            fmuladd(c, diag[k], -(s * sub[k])),
            -(s * fmuladd(c, sub[k], -(s * diag[k + 1]))),
        );
        diag[k + 1] = fmuladd(s, sdk, c * dkp1);
        sub[k] = fmuladd(c, sdk, -(s * dkp1));
        if k > start {
            sub[k - 1] = fmuladd(c, sub[k - 1], -(s * z));
        }
        // chasing the bulge
        x = sub[k];
        if k < end - 1 {
            z = -s * sub[k + 1];
            sub[k + 1] *= c;
        }
        // Q = Q G
        for row in q.iter_mut() {
            let (xi, yi) = (row[k], row[k + 1]);
            row[k] = fmuladd(c, xi, -(s * yi));
            row[k + 1] = fmuladd(s, xi, c * yi);
        }
        k += 1;
    }
}

/// `SelfAdjointEigenSolver<Matrix3d>(m)`: only the lower triangle of `m` is
/// read.
pub fn symmetric_eigen_3x3(m: &[[f64; N]; N]) -> Result<SymmetricEigen3> {
    // map the coefficients of the lower triangle to [-1, 1]
    let mut mat = [[0.0; N]; N];
    for i in 0..N {
        mat[i][..=i].copy_from_slice(&m[i][..=i]);
    }
    let mut scale = mat.iter().flatten().fold(0.0_f64, |a, x| a.max(x.abs()));
    if scale == 0.0 {
        scale = 1.0;
    }
    for i in 0..N {
        for j in 0..=i {
            mat[i][j] /= scale;
        }
    }
    let (mut diag, mut sub, mut q) = tridiagonalize(&mat);

    // computeFromTridiagonal_impl
    let precision_inv = 1.0 / f64::EPSILON;
    let mut end = N - 1;
    let mut start = 0;
    let mut iter = 0;
    while end > 0 {
        for i in start..end {
            if sub[i].abs() < f64::MIN_POSITIVE {
                sub[i] = 0.0;
            } else {
                let scaled = precision_inv * sub[i];
                if scaled * scaled <= diag[i].abs() + diag[i + 1].abs() {
                    sub[i] = 0.0;
                }
            }
        }
        // the largest unreduced block at the end of the matrix
        while end > 0 && sub[end - 1] == 0.0 {
            end -= 1;
        }
        if end == 0 {
            break;
        }
        iter += 1;
        if iter > MAX_ITERATIONS * N {
            return Err(RDKitError::Runtime(
                "eigenvalue calculation did not converge".into(),
            ));
        }
        start = end - 1;
        while start > 0 && sub[start - 1] != 0.0 {
            start -= 1;
        }
        tridiagonal_qr_step(&mut diag, &mut sub, start, end, &mut q);
    }

    // selection sort into ascending order, swapping the eigenvector columns
    for i in 0..N - 1 {
        let mut k = i;
        for j in i + 1..N {
            if diag[j] < diag[k] {
                k = j;
            }
        }
        if k > i {
            diag.swap(i, k);
            for row in q.iter_mut() {
                row.swap(i, k);
            }
        }
    }
    Ok(SymmetricEigen3 {
        values: diag.map(|d| d * scale),
        vectors: q,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residual(m: &[[f64; N]; N], e: &SymmetricEigen3) -> f64 {
        let mut worst = 0.0_f64;
        for col in 0..N {
            for row in 0..N {
                let mv: f64 = (0..N).map(|k| m[row][k] * e.vectors[k][col]).sum();
                worst = worst.max((mv - e.values[col] * e.vectors[row][col]).abs());
            }
        }
        worst
    }

    #[test]
    fn diagonal_matrix_keeps_identity_vectors() {
        let m = [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]];
        let e = symmetric_eigen_3x3(&m).unwrap();
        assert_eq!(e.values, [1.0, 2.0, 3.0]);
        assert!(residual(&m, &e) < 1e-15);
    }

    #[test]
    fn general_symmetric_matrix() {
        let m = [[4.0, 1.0, -2.0], [1.0, 2.0, 0.5], [-2.0, 0.5, 3.0]];
        let e = symmetric_eigen_3x3(&m).unwrap();
        assert!(e.values[0] <= e.values[1] && e.values[1] <= e.values[2]);
        assert!(residual(&m, &e) < 1e-13);
        // orthonormal columns
        for a in 0..N {
            for b in 0..N {
                let dot: f64 = (0..N).map(|r| e.vectors[r][a] * e.vectors[r][b]).sum();
                assert!((dot - if a == b { 1.0 } else { 0.0 }).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn zero_matrix() {
        let e = symmetric_eigen_3x3(&[[0.0; N]; N]).unwrap();
        assert_eq!(e.values, [0.0; N]);
    }
}
