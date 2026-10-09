//! Dense column-major matrices with products evaluated in the order and with
//! the fused multiply-adds of Eigen's `MatrixXd` (see `rdgeneral::fp`).
//!
//! Eigen evaluates a product coefficient by coefficient when
//! `rhs.rows() + dst.rows() + dst.cols() < 20` and with its blocked kernel
//! otherwise. The small products are computed in registers of two rows at a
//! time, so which coefficients use fused multiply-adds depends on their row and
//! on the alignment of their column; the blocked kernel fuses everything.
//! Results agree with Eigen bit for bit for the products used by RDKit's
//! descriptors, and to rounding error in general.
//!
//! Translated from the behaviour of Eigen 5.0.1 (`Core/products/*.h`,
//! `Core/Redux.h`).
//!
//! This Source Code Form is subject to the terms of the Mozilla Public
//! License v. 2.0 (see LICENSES/Eigen-MPL-2.0.txt).

use crate::rdgeneral::fp::fmuladd;

/// Doubles per SIMD register of the reference build (NEON).
const PACKET: usize = 2;
/// `EIGEN_GEMM_TO_COEFFBASED_THRESHOLD`.
pub const SMALL_PRODUCT: usize = 20;

/// A dense matrix stored column by column.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Matrix {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = Matrix::zeros(n, n);
        for i in 0..n {
            m.set(i, i, 1.0);
        }
        m
    }

    /// The matrix whose `(row, col)` entry is `f(row, col)`.
    pub fn from_fn(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> f64) -> Self {
        let mut m = Matrix::zeros(rows, cols);
        for col in 0..cols {
            for row in 0..rows {
                m.set(row, col, f(row, col));
            }
        }
        m
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[col * self.rows + row]
    }

    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        self.data[col * self.rows + row] = value;
    }

    /// Column `col` as a contiguous slice.
    pub fn column(&self, col: usize) -> &[f64] {
        &self.data[col * self.rows..(col + 1) * self.rows]
    }

    pub(crate) fn column_mut(&mut self, col: usize) -> &mut [f64] {
        &mut self.data[col * self.rows..(col + 1) * self.rows]
    }

    pub fn transpose(&self) -> Matrix {
        Matrix::from_fn(self.cols, self.rows, |r, c| self.get(c, r))
    }

    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.rows.min(self.cols))
            .map(|i| self.get(i, i))
            .collect()
    }

    /// Sum of every row (`rowwise().sum()`): the columns are added in order.
    pub fn row_sums(&self) -> Vec<f64> {
        (0..self.rows)
            .map(|r| (1..self.cols).fold(self.get(r, 0), |acc, c| acc + self.get(r, c)))
            .collect()
    }

    /// Mean of every column (`colwise().mean()`).
    pub fn column_means(&self) -> Vec<f64> {
        (0..self.cols)
            .map(|c| packet_sum(self.column(c), c * self.rows) / self.rows as f64)
            .collect()
    }

    /// `self * rhs` as Eigen evaluates `MatrixXd * MatrixXd`.
    pub fn mul(&self, rhs: &Matrix) -> Matrix {
        assert_eq!(self.cols, rhs.rows, "inner dimensions must agree");
        let (rows, cols, depth) = (self.rows, rhs.cols, self.cols);
        let small = depth > 0 && depth + rows + cols < SMALL_PRODUCT;
        Matrix::from_fn(rows, cols, |i, j| {
            let fused = !small || in_packet(i, rows, j * rows);
            dot(
                (0..depth).map(|k| self.get(i, k)),
                (0..depth).map(|k| rhs.get(k, j)),
                fused,
            )
        })
    }

    /// `self^T * rhs` as Eigen evaluates `A.transpose() * B`: the small
    /// products reduce each pair of contiguous columns like `sum()`.
    pub fn transpose_mul(&self, rhs: &Matrix) -> Matrix {
        assert_eq!(self.rows, rhs.rows, "inner dimensions must agree");
        let (rows, cols, depth) = (self.cols, rhs.cols, self.rows);
        let small = depth > 0 && depth + rows + cols < SMALL_PRODUCT;
        Matrix::from_fn(rows, cols, |i, j| {
            if small {
                let products: Vec<f64> = self
                    .column(i)
                    .iter()
                    .zip(rhs.column(j))
                    .map(|(a, b)| a * b)
                    .collect();
                packet_sum(&products, 0)
            } else {
                dot(
                    self.column(i).iter().copied(),
                    rhs.column(j).iter().copied(),
                    true,
                )
            }
        })
    }
}

/// Dot product of two vectors the way Eigen sums the products: `fused` for
/// the multiply-add chains of its vectorized kernels, otherwise a plain sum of
/// separately rounded products.
pub fn dot(a: impl Iterator<Item = f64>, b: impl Iterator<Item = f64>, fused: bool) -> f64 {
    let mut pairs = a.zip(b);
    let Some((x, y)) = pairs.next() else {
        return 0.0;
    };
    pairs.fold(x * y, |acc, (x, y)| {
        if fused {
            fmuladd(x, y, acc)
        } else {
            acc + x * y
        }
    })
}

/// Whether row `row` of a column of `rows` entries starting `offset` elements
/// after a register boundary is computed by Eigen's vectorized path: the rows
/// between the first aligned element and the last full register.
pub fn in_packet(row: usize, rows: usize, offset: usize) -> bool {
    let start = offset % PACKET;
    let end = start + (rows - start.min(rows)) / PACKET * PACKET;
    row >= start && row < end
}

/// `sum()` of a contiguous vector whose first element is `offset` elements
/// past a register boundary (Eigen's vectorized reduction: two interleaved
/// register accumulators, the lanes added at the end, then the unaligned head
/// and tail one by one).
pub fn packet_sum(v: &[f64], offset: usize) -> f64 {
    let start = (offset % PACKET).min(v.len());
    let size = v.len() - start;
    let aligned = size / PACKET * PACKET;
    if aligned == 0 {
        return v.iter().skip(1).fold(v[0], |acc, x| acc + x);
    }
    let end = start + aligned;
    let end2 = start + size / (2 * PACKET) * (2 * PACKET);
    let load = |i: usize| [v[i], v[i + 1]];
    let mut r0 = load(start);
    if aligned > PACKET {
        let mut r1 = load(start + PACKET);
        let mut i = start + 2 * PACKET;
        while i < end2 {
            let (a, b) = (load(i), load(i + PACKET));
            r0 = [r0[0] + a[0], r0[1] + a[1]];
            r1 = [r1[0] + b[0], r1[1] + b[1]];
            i += 2 * PACKET;
        }
        r0 = [r0[0] + r1[0], r0[1] + r1[1]];
        if end > end2 {
            let a = load(end2);
            r0 = [r0[0] + a[0], r0[1] + a[1]];
        }
    }
    let mut res = r0[0] + r0[1];
    for x in &v[..start] {
        res += x;
    }
    for x in &v[end..] {
        res += x;
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn product_of_small_matrices() {
        let a = Matrix::from_fn(2, 3, |r, c| (r * 3 + c + 1) as f64);
        let b = Matrix::from_fn(3, 2, |r, c| (r * 2 + c + 1) as f64);
        let p = a.mul(&b);
        assert_eq!(
            (p.get(0, 0), p.get(0, 1), p.get(1, 0), p.get(1, 1)),
            (22.0, 28.0, 49.0, 64.0)
        );
        let t = a.transpose().transpose_mul(&b);
        assert_eq!(t, p);
    }

    #[test]
    fn packet_sum_adds_everything() {
        let v: Vec<f64> = (1..=11).map(f64::from).collect();
        for offset in 0..2 {
            assert_eq!(packet_sum(&v, offset), 66.0);
            assert_eq!(packet_sum(&v[..3], offset), 6.0);
        }
        assert_eq!(packet_sum(&[5.0], 1), 5.0);
    }

    #[test]
    fn packet_membership_follows_the_alignment() {
        // a column of five rows starting on a register boundary: rows 0-3
        assert!(
            (0..5)
                .map(|r| in_packet(r, 5, 0))
                .eq([true, true, true, true, false])
        );
        // starting one element later: rows 1-4
        assert!(
            (0..5)
                .map(|r| in_packet(r, 5, 1))
                .eq([false, true, true, true, true])
        );
    }

    #[test]
    fn means_and_row_sums() {
        let m = Matrix::from_fn(3, 2, |r, c| (r + 3 * c) as f64);
        assert_eq!(m.column_means(), vec![1.0, 4.0]);
        assert_eq!(m.row_sums(), vec![3.0, 5.0, 7.0]);
    }
}
