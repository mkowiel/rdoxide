//! Small dense linear-algebra helpers.
//!
//! RDKit gets these from Eigen; the routines here are Householder reduction
//! plus implicit QL iteration (the classic EISPACK `tred1`/`tql1` pair), which
//! agree with Eigen to rounding error but not bit for bit — Eigen's own
//! results already depend on the SIMD width it was built with.

pub mod alignment;
pub mod bfgs;
pub mod conrec;
pub mod dense;
pub mod eigen_3x3;
pub mod full_piv_lu;
pub mod jacobi_svd;
pub mod power_eigen;
pub mod symm_matrix;

use crate::{RDKitError, Result};

const MAX_QL_ITERATIONS: usize = 60;

/// Eigenvalues, in ascending order, of the symmetric `n x n` matrix stored
/// row-major in `a` (only the lower triangle is read).
pub fn symmetric_eigenvalues(n: usize, a: &[f64]) -> Result<Vec<f64>> {
    assert_eq!(a.len(), n * n, "matrix must be n x n");
    if n == 0 {
        return Ok(Vec::new());
    }
    let (mut d, mut e) = tridiagonalize(n, a);
    ql_implicit(&mut d, &mut e)?;
    d.sort_by(f64::total_cmp);
    Ok(d)
}

/// Householder reduction to tridiagonal form: the diagonal `d` and the
/// sub-diagonal `e` (with `e[0]` unused).
fn tridiagonalize(n: usize, lower: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut z = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            z[i * n + j] = lower[i * n + j];
        }
    }
    let at = |z: &[f64], i: usize, j: usize| z[i * n + j];
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    for i in (1..n).rev() {
        let l = i - 1;
        let mut h = 0.0;
        if l > 0 {
            let scale: f64 = (0..=l).map(|k| at(&z, i, k).abs()).sum();
            if scale == 0.0 {
                e[i] = at(&z, i, l);
            } else {
                for k in 0..=l {
                    z[i * n + k] /= scale;
                    h += at(&z, i, k) * at(&z, i, k);
                }
                let f = at(&z, i, l);
                let g = if f >= 0.0 { -h.sqrt() } else { h.sqrt() };
                e[i] = scale * g;
                h -= f * g;
                z[i * n + l] = f - g;
                let mut f_sum = 0.0;
                for j in 0..=l {
                    // the lower triangle is all that is kept, so row j is read up
                    // to the diagonal and column j below it
                    let mut g = 0.0;
                    for k in 0..=j {
                        g += at(&z, j, k) * at(&z, i, k);
                    }
                    for k in j + 1..=l {
                        g += at(&z, k, j) * at(&z, i, k);
                    }
                    e[j] = g / h;
                    f_sum += e[j] * at(&z, i, j);
                }
                let hh = f_sum / (h + h);
                for j in 0..=l {
                    let f = at(&z, i, j);
                    let g = e[j] - hh * f;
                    e[j] = g;
                    for k in 0..=j {
                        z[j * n + k] -= f * e[k] + g * at(&z, i, k);
                    }
                }
            }
        } else {
            e[i] = at(&z, i, l);
        }
        d[i] = h;
    }
    for (i, di) in d.iter_mut().enumerate() {
        *di = z[i * n + i];
    }
    (d, e)
}

/// Implicit QL iteration on the tridiagonal matrix `(d, e)`; the eigenvalues
/// end up in `d` (unsorted).
fn ql_implicit(d: &mut [f64], e: &mut [f64]) -> Result<()> {
    let n = d.len();
    // renumber the sub-diagonal so that e[i] couples d[i] and d[i+1]
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    for l in 0..n {
        let mut iterations = 0;
        loop {
            let mut m = l;
            while m + 1 < n {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iterations += 1;
            if iterations > MAX_QL_ITERATIONS {
                return Err(RDKitError::Runtime(
                    "eigenvalue iteration did not converge".into(),
                ));
            }
            // Wilkinson shift
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + r.copysign(g));
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut underflow = false;
            for i in (l..m).rev() {
                let f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    // recover from underflow
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    underflow = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
            }
            if underflow {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() < 1e-12, "{a:?} != {b:?}");
        }
    }

    #[test]
    fn trivial_sizes() {
        assert!(symmetric_eigenvalues(0, &[]).unwrap().is_empty());
        close(&symmetric_eigenvalues(1, &[3.5]).unwrap(), &[3.5]);
    }

    #[test]
    fn two_by_two() {
        // [[2, 1], [1, 2]] has eigenvalues 1 and 3
        close(
            &symmetric_eigenvalues(2, &[2.0, 1.0, 1.0, 2.0]).unwrap(),
            &[1.0, 3.0],
        );
    }

    #[test]
    fn diagonal_matrix_is_sorted() {
        let a = [3.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.0];
        close(&symmetric_eigenvalues(3, &a).unwrap(), &[-1.0, 2.0, 3.0]);
    }

    #[test]
    fn path_graph_adjacency() {
        // adjacency of the 4-path: 2cos(k pi / 5), k = 1..4
        let a = [
            0.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0, //
            0.0, 1.0, 0.0, 1.0, //
            0.0, 0.0, 1.0, 0.0,
        ];
        let expected: Vec<f64> = (1..=4)
            .rev()
            .map(|k| 2.0 * (f64::from(k) * std::f64::consts::PI / 5.0).cos())
            .collect();
        close(&symmetric_eigenvalues(4, &a).unwrap(), &expected);
    }

    #[test]
    fn trace_and_frobenius_norm_are_preserved() {
        let n = 12;
        let mut a = vec![0.0; n * n];
        let mut x = 0.3_f64;
        for i in 0..n {
            for j in 0..=i {
                x = (x * 7.31 + 0.17).fract();
                a[i * n + j] = x - 0.5;
                a[j * n + i] = x - 0.5;
            }
        }
        let ev = symmetric_eigenvalues(n, &a).unwrap();
        let trace: f64 = (0..n).map(|i| a[i * n + i]).sum();
        let frob: f64 = a.iter().map(|v| v * v).sum();
        assert!((ev.iter().sum::<f64>() - trace).abs() < 1e-12);
        assert!((ev.iter().map(|v| v * v).sum::<f64>() - frob).abs() < 1e-11);
        assert!(ev.windows(2).all(|w| w[0] <= w[1]));
    }
}
