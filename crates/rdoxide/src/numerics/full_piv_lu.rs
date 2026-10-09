//! LU decomposition with full pivoting of a square matrix that reproduces
//! Eigen's `FullPivLU<MatrixXd>` and its `solve`: the same pivot choice, the
//! same panel structure of the triangular solves (panels of eight rows) and
//! the same order of the floating-point operations, including the products
//! that Eigen's vectorized matrix-vector kernel fuses on aarch64
//! (`rdgeneral::fp::fmuladd`).
//!
//! Translated from Eigen 5.0.1: `LU/FullPivLU.h` (`computeInPlace`,
//! `_solve_impl`, `rank`), `Core/products/TriangularSolverVector.h` and the
//! column-major kernel of `Core/products/GeneralMatrixVector.h`.
//!
//! Copyright (C) 2006-2009 Benoit Jacob <jacob.benoit.1@gmail.com>.
//!
//! This Source Code Form is subject to the terms of the Mozilla Public
//! License v. 2.0 (see LICENSES/Eigen-MPL-2.0.txt).

use crate::rdgeneral::fp::fmuladd;

/// Eigen's `EIGEN_TUNE_TRIANGULAR_PANEL_WIDTH`.
const PANEL_WIDTH: usize = 8;
/// Doubles per SIMD register of the reference build (NEON); the matrix-vector
/// kernel fuses the products of whole register pairs of rows only.
const PACKET: usize = 2;

/// `FullPivLU` of a square matrix.
#[derive(Debug)]
pub struct FullPivLu {
    n: usize,
    /// column-major, `L` (unit diagonal, below) and `U` (on and above)
    lu: Vec<f64>,
    /// row permutation `P`: `(P x)[p[i]] = x[i]`
    p: Vec<usize>,
    /// column permutation `Q`
    q: Vec<usize>,
    nonzero_pivots: usize,
    max_pivot: f64,
}

impl FullPivLu {
    /// Decomposes the `n x n` column-major matrix `matrix`.
    pub fn new(n: usize, matrix: Vec<f64>) -> Self {
        assert_eq!(matrix.len(), n * n, "matrix must be n x n");
        let mut lu = FullPivLu {
            n,
            lu: matrix,
            p: (0..n).collect(),
            q: (0..n).collect(),
            nonzero_pivots: n,
            max_pivot: 0.0,
        };
        lu.compute_in_place();
        lu
    }

    fn at(&self, row: usize, col: usize) -> f64 {
        self.lu[col * self.n + row]
    }

    fn compute_in_place(&mut self) {
        let n = self.n;
        let mut row_swaps = (0..n).collect::<Vec<_>>();
        let mut col_swaps = (0..n).collect::<Vec<_>>();
        for k in 0..n {
            // the biggest coefficient of the bottom-right corner, the first
            // one in column-major order on ties
            let (mut biggest, mut pivot_row, mut pivot_col) = (-1.0, k, k);
            for col in k..n {
                for row in k..n {
                    let score = self.at(row, col).abs();
                    if score > biggest {
                        (biggest, pivot_row, pivot_col) = (score, row, col);
                    }
                }
            }
            if biggest == 0.0 {
                self.nonzero_pivots = k;
                break;
            }
            self.max_pivot = self.max_pivot.max(biggest);
            row_swaps[k] = pivot_row;
            col_swaps[k] = pivot_col;
            if k != pivot_row {
                for col in 0..n {
                    self.lu.swap(col * n + k, col * n + pivot_row);
                }
            }
            if k != pivot_col {
                for row in 0..n {
                    self.lu.swap(k * n + row, pivot_col * n + row);
                }
            }
            // eliminate below the pivot
            let pivot = self.at(k, k);
            for row in k + 1..n {
                self.lu[k * n + row] /= pivot;
            }
            for col in k + 1..n {
                let u = self.at(k, col);
                for row in k + 1..n {
                    // Eigen's outer-product update multiplies and subtracts separately
                    let product = u * self.lu[k * n + row];
                    self.lu[col * n + row] -= product;
                }
            }
        }
        for k in (0..n).rev() {
            self.p.swap(k, row_swaps[k]);
        }
        for k in 0..n {
            self.q.swap(k, col_swaps[k]);
        }
    }

    /// Number of pivots above Eigen's default threshold.
    pub fn rank(&self) -> usize {
        let threshold = self.max_pivot.abs() * (f64::EPSILON * self.n as f64);
        (0..self.nonzero_pivots)
            .filter(|&i| self.at(i, i).abs() > threshold)
            .count()
    }

    /// `res[..rows] -= A[row0.., col0..col0+cols] * x` where `x` is
    /// `v[col0..col0+cols]` and the result is `v[row0..row0+rows]`
    /// (`general_matrix_vector_product` with alpha = -1).
    fn gemv_subtract(&self, v: &mut [f64], row0: usize, rows: usize, col0: usize, cols: usize) {
        let fused_rows = rows - rows % PACKET;
        for i in 0..rows {
            let mut acc = 0.0;
            for j in 0..cols {
                let a = self.at(row0 + i, col0 + j);
                acc = if i < fused_rows {
                    fmuladd(a, v[col0 + j], acc)
                } else {
                    acc + a * v[col0 + j]
                };
            }
            v[row0 + i] -= acc;
        }
    }

    /// Solves `A x = rhs`.
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        assert_eq!(rhs.len(), self.n, "right-hand side must have n entries");
        let rank = self.rank();
        let mut dst = vec![0.0; self.n];
        if rank == 0 {
            return dst;
        }
        let mut c = vec![0.0; self.n];
        for (i, &x) in rhs.iter().enumerate() {
            c[self.p[i]] = x;
        }
        self.solve_unit_lower(&mut c);
        self.solve_upper(&mut c[..rank]);
        for i in 0..rank {
            dst[self.q[i]] = c[i];
        }
        dst
    }

    /// Forward substitution with the unit lower triangle.
    fn solve_unit_lower(&self, c: &mut [f64]) {
        let size = self.n;
        let mut pi = 0;
        while pi < size {
            let width = PANEL_WIDTH.min(size - pi);
            for k in 0..width {
                let i = pi + k;
                if c[i] != 0.0 {
                    let r = width - k - 1;
                    for t in 1..=r {
                        let product = c[i] * self.at(i + t, i);
                        c[i + t] -= product;
                    }
                }
            }
            let end = pi + width;
            if size > end {
                self.gemv_subtract(c, end, size - end, pi, width);
            }
            pi += PANEL_WIDTH;
        }
    }

    /// Back substitution with the leading `c.len() x c.len()` upper triangle.
    fn solve_upper(&self, c: &mut [f64]) {
        let size = c.len();
        let mut pi = size;
        while pi > 0 {
            let width = PANEL_WIDTH.min(pi);
            let start = pi - width;
            for k in 0..width {
                let i = pi - k - 1;
                if c[i] != 0.0 {
                    c[i] /= self.at(i, i);
                    let r = width - k - 1;
                    for t in 1..=r {
                        let product = c[i] * self.at(i - t, i);
                        c[i - t] -= product;
                    }
                }
            }
            if start > 0 {
                self.gemv_subtract(c, 0, start, start, width);
            }
            pi = start;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Column-major copy of a row-major matrix.
    fn col_major(n: usize, rows: &[f64]) -> Vec<f64> {
        (0..n * n).map(|m| rows[(m % n) * n + m / n]).collect()
    }

    #[test]
    fn solves_a_system_that_needs_both_pivot_directions() {
        let a = [0.0, 2.0, 1.0, 4.0, 1.0, 0.0, 1.0, 0.0, 3.0];
        let lu = FullPivLu::new(3, col_major(3, &a));
        assert_eq!(lu.rank(), 3);
        let x = lu.solve(&[3.0, 5.0, 4.0]);
        for (row, b) in a.chunks(3).zip([3.0, 5.0, 4.0]) {
            let r: f64 = row.iter().zip(&x).map(|(a, x)| a * x).sum();
            assert!((r - b).abs() < 1e-12, "{r} != {b}");
        }
    }

    #[test]
    fn large_systems_cross_the_panel_boundary() {
        let n = 20;
        let a: Vec<f64> = (0..n * n)
            .map(|m| {
                let (i, j) = (m / n, m % n);
                if i == j {
                    10.0 + i as f64
                } else {
                    1.0 / (1.0 + (i + 2 * j) as f64)
                }
            })
            .collect();
        let b: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
        let x = FullPivLu::new(n, col_major(n, &a)).solve(&b);
        for i in 0..n {
            let r: f64 = (0..n).map(|j| a[i * n + j] * x[j]).sum();
            assert!((r - b[i]).abs() < 1e-12, "row {i}: {r} != {}", b[i]);
        }
    }

    #[test]
    fn singular_matrix_has_reduced_rank() {
        let lu = FullPivLu::new(2, col_major(2, &[1.0, 2.0, 2.0, 4.0]));
        assert_eq!(lu.rank(), 1);
    }
}
