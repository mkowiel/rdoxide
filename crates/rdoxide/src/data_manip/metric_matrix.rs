//! Metric matrices (upstream `DataManip/MetricMatrixCalc`): the Euclidean
//! distance matrix of a descriptor table and Tanimoto similarity or distance
//! matrices of bit vectors.
//!
//! A matrix over `n` items is returned as its lower triangle, row by row
//! without the diagonal (`n * (n - 1) / 2` values; entry `(i, j)` with
//! `j < i` sits at `i * (i - 1) / 2 + j`), the layout the diversity pickers
//! and clustering code take. Rows are independent, so they are computed on
//! up to `num_threads` threads; the result does not depend on the count.

use crate::datastructs::bit_ops::{BitOpsVect, Metric, similarity};
use crate::error::Result;
use crate::rdgeneral::threads::par_map;

/// `MetricMatrixCalc::calcMetricMatrix`: the lower triangle of the matrix of
/// `metric(i, j)` over `n` items, computed on up to `num_threads` threads
/// (zero or negative: all hardware threads minus `|num_threads|`).
pub fn lower_triangle<F>(n: usize, num_threads: i32, metric: F) -> Result<Vec<f64>>
where
    F: Fn(usize, usize) -> Result<f64> + Sync,
{
    let rows: Vec<usize> = (1..n).collect();
    let per_row = par_map(&rows, num_threads, |&i| {
        (0..i).map(|j| metric(i, j)).collect::<Result<Vec<f64>>>()
    });
    let mut matrix = Vec::with_capacity(n * n.saturating_sub(1) / 2);
    for row in per_row {
        matrix.extend(row?);
    }
    Ok(matrix)
}

/// `EuclideanDistanceMetric`
fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .fold(0.0, |sum, (x, y)| sum + (x - y) * (x - y))
        .sqrt()
}

/// The Euclidean distance matrix of `n_rows` descriptor vectors of
/// `n_cols` values each (`table` is row-major).
pub fn euclidean_distance_matrix(
    table: &[f64],
    n_rows: usize,
    n_cols: usize,
    num_threads: i32,
) -> Result<Vec<f64>> {
    if table.len() < n_rows * n_cols {
        return Err(crate::RDKitError::Value(
            "the descriptor table is smaller than its shape".into(),
        ));
    }
    let row = |i: usize| &table[i * n_cols..(i + 1) * n_cols];
    lower_triangle(n_rows, num_threads, |i, j| {
        Ok(euclidean_distance(row(i), row(j)))
    })
}

/// The Tanimoto similarity matrix of the bit vectors (`as_distance`: the
/// distance, `1 - similarity`, instead).
pub fn tanimoto_matrix<T>(items: &[&T], as_distance: bool, num_threads: i32) -> Result<Vec<f64>>
where
    T: BitOpsVect + Sync,
{
    lower_triangle(items.len(), num_threads, |i, j| {
        let sim = similarity(Metric::Tanimoto, items[i], items[j])?;
        Ok(if as_distance { 1.0 - sim } else { sim })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::datastructs::{BitVect, ExplicitBitVect};

    #[test]
    fn euclidean_matrix_of_a_triangle() {
        let table = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0];
        let d = euclidean_distance_matrix(&table, 3, 2, 1).unwrap();
        assert_eq!(d, [1.0, 2f64.sqrt(), 1.0]);
        assert!(euclidean_distance_matrix(&table, 4, 2, 1).is_err());
    }

    #[test]
    fn tanimoto_matrices() {
        let make = |bits: &[u32]| {
            let mut bv = ExplicitBitVect::new(8, false);
            bits.iter().for_each(|&b| {
                bv.set_bit(b).unwrap();
            });
            bv
        };
        let (a, b, c) = (make(&[0, 1]), make(&[1, 2]), make(&[5]));
        let items = [&a, &b, &c];
        let sims = tanimoto_matrix(&items, false, 1).unwrap();
        assert_eq!(sims, [1.0 / 3.0, 0.0, 0.0]);
        let dists = tanimoto_matrix(&items, true, 1).unwrap();
        assert_eq!(dists, [1.0 - 1.0 / 3.0, 1.0, 1.0]);
        let other_size = ExplicitBitVect::new(4, false);
        assert!(tanimoto_matrix(&[&a, &other_size], false, 1).is_err());
    }

    #[test]
    fn threads_do_not_change_the_result() {
        let table: Vec<f64> = (0..4000).map(|i| f64::from(i * 31 % 97)).collect();
        let serial = euclidean_distance_matrix(&table, 1000, 4, 1).unwrap();
        assert_eq!(serial.len(), 1000 * 999 / 2);
        assert_eq!(
            serial,
            euclidean_distance_matrix(&table, 1000, 4, 0).unwrap()
        );
    }
}
