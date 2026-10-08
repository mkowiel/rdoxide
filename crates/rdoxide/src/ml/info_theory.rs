//! Information-theory helpers (upstream `ML/InfoTheory/InfoGainFuncs.h`):
//! entropy, information gain and chi-squared of occurrence counts.
//!
//! Upstream instantiates these templates for the element types of the numpy
//! array it is given, so the sums run in that type (`f32` and the integer
//! types accumulate in their own precision). The same is done here with the
//! [`Count`] trait. Floating-point expressions keep upstream's operation order
//! and fused multiply-adds (see [`crate::rdgeneral::fp`]).

use crate::rdgeneral::fp::fmuladd;

/// An element type of the count arrays.
pub trait Count: Copy + PartialEq {
    const ZERO: Self;
    /// `a + b` (integers wrap, as the C++ arithmetic does in practice)
    fn add(self, other: Self) -> Self;
    /// `(double)x`
    fn to_f64(self) -> f64;
    /// `(int)x`
    fn to_i32(self) -> i32;
}

macro_rules! impl_count {
    ($t:ty, $zero:expr, |$a:ident, $b:ident| $add:expr) => {
        impl Count for $t {
            const ZERO: Self = $zero;
            fn add(self, other: Self) -> Self {
                let ($a, $b) = (self, other);
                $add
            }
            fn to_f64(self) -> f64 {
                self as f64
            }
            fn to_i32(self) -> i32 {
                self as i32
            }
        }
    };
}

impl_count!(f32, 0.0, |a, b| a + b);
impl_count!(f64, 0.0, |a, b| a + b);
impl_count!(i32, 0, |a, b| a.wrapping_add(b));
impl_count!(i64, 0, |a, b| a.wrapping_add(b));

fn sum<T: Count>(values: impl Iterator<Item = T>) -> T {
    values.fold(T::ZERO, T::add)
}

/// `InfoEntropy`: the entropy (in bits) of the distribution given by the counts.
pub fn info_entropy<T: Count>(counts: &[T]) -> f64 {
    let n_instances = sum(counts.iter().copied());
    let mut accum = 0.0;
    if n_instances != T::ZERO {
        for &count in counts {
            let d = count.to_f64() / n_instances.to_f64();
            if d != 0.0 {
                accum = fmuladd(-d, d.ln(), accum);
            }
        }
    }
    accum / std::f64::consts::LN_2
}

fn row_sums<T: Count>(mat: &[T], dim2: usize) -> Vec<T> {
    mat.chunks(dim2.max(1))
        .map(|row| sum(row.iter().copied()))
        .collect()
}

fn col_sums<T: Count>(mat: &[T], dim1: usize, dim2: usize) -> Vec<T> {
    (0..dim2)
        .map(|i| sum((0..dim1).map(|j| mat[j * dim2 + i])))
        .collect()
}

/// `InfoEntropyGain`: the information gain of a variable. `mat` is
/// `dim1 x dim2` (row-major): the number of occurrences of each of the
/// `dim2` results for each of the `dim1` values of the variable.
pub fn info_entropy_gain<T: Count>(mat: &[T], dim1: usize, dim2: usize) -> f64 {
    let variable_res = row_sums(mat, dim2);
    let overall_res = col_sums(mat, dim1, dim2);

    let mut term2 = 0.0;
    for (i, &res) in variable_res.iter().enumerate().take(dim1) {
        term2 = fmuladd(
            res.to_f64(),
            info_entropy(&mat[i * dim2..(i + 1) * dim2]),
            term2,
        );
    }
    let t_sum = overall_res
        .iter()
        .fold(0i32, |acc, &r| acc.wrapping_add(r.to_i32()));
    if t_sum != 0 {
        term2 /= f64::from(t_sum);
        info_entropy(&overall_res) - term2
    } else {
        0.0
    }
}

/// `ChiSquare` of a contingency matrix with a column per class and a row per
/// state of the descriptor: `chi = sum_i (N / R_i) * sum_j (N_ij^2 / C_j) - N`.
pub fn chi_square<T: Count>(mat: &[T], dim1: usize, dim2: usize) -> f64 {
    let row_sums = row_sums(mat, dim2);
    let col_sums = col_sums(mat, dim1, dim2);
    let t_sum = row_sums
        .iter()
        .take(dim1)
        .fold(0i32, |acc, &r| acc.wrapping_add(r.to_i32()));

    let mut chi = 0.0;
    for (i, &row_sum) in row_sums.iter().enumerate().take(dim1) {
        let mut rchi = 0.0;
        for (j, &col_sum) in col_sums.iter().enumerate() {
            let n_ij = mat[i * dim2 + j].to_f64();
            rchi += (n_ij * n_ij) / col_sum.to_f64();
        }
        chi = fmuladd(f64::from(t_sum) / row_sum.to_f64(), rchi, chi);
    }
    chi - f64::from(t_sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entropy_of_uniform_and_degenerate_distributions() {
        assert_eq!(info_entropy(&[1i32, 1]), 1.0);
        assert_eq!(info_entropy(&[4.0f64, 4.0, 4.0, 4.0]), 2.0);
        assert_eq!(info_entropy(&[7i64, 0]), 0.0);
        assert_eq!(info_entropy::<i32>(&[0, 0]), 0.0);
    }

    #[test]
    fn perfectly_informative_variable_gains_all_the_entropy() {
        // two values of the variable, each seeing a single result
        assert_eq!(info_entropy_gain(&[3i32, 0, 0, 3], 2, 2), 1.0);
        // an uninformative variable gains nothing
        assert_eq!(info_entropy_gain(&[2i32, 2, 2, 2], 2, 2), 0.0);
        assert_eq!(info_entropy_gain(&[0i32, 0, 0, 0], 2, 2), 0.0);
    }

    #[test]
    fn chi_square_of_independent_and_dependent_tables() {
        assert!(chi_square(&[2i32, 2, 2, 2], 2, 2).abs() < 1e-12);
        assert!((chi_square(&[3i32, 0, 0, 3], 2, 2) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn float32_counts_accumulate_in_single_precision() {
        let counts = [16_777_216.0f32, 1.0];
        // 2^24 + 1 is not representable as f32, so the total stays 2^24
        let entropy = info_entropy(&counts);
        let in_double = info_entropy(&[16_777_216.0f64, 1.0]);
        assert_ne!(entropy, in_double);
    }
}
