//! Power iteration for the leading eigenpairs of a symmetric matrix (upstream
//! `Numerics/EigenSolvers/PowerEigenSolver.cpp`).

use super::symm_matrix::SymmMatrix;
use crate::rdgeneral::random::MinStdRand;

/// A random unit vector (`Vector::setToRandom`), from a private
/// `minstd_rand` seeded with `seed` (the clock when it is not positive).
fn random_unit_vector(n: usize, seed: i64) -> Vec<f64> {
    let seed = if seed > 0 {
        seed as u32
    } else {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(1, |d| d.subsec_nanos())
            .wrapping_add(1)
    };
    let mut rng = MinStdRand::new(seed);
    let mut v: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in &mut v {
        *x /= norm;
    }
    v
}

/// Computes the `num_eig` largest eigenvalues (by magnitude) of `mat`, which
/// is deflated in place. `eigen_vectors` receives `num_eig` rows of `n`
/// values. Returns whether the last iteration converged.
pub fn power_eigen_solver(
    num_eig: usize,
    mat: &mut SymmMatrix,
    eigen_values: &mut [f64],
    mut eigen_vectors: Option<&mut [f64]>,
    seed: i32,
) -> bool {
    const MAX_ITERATIONS: usize = 1000;
    const TOLERANCE: f64 = 0.001;
    const HUGE_EIGVAL: f64 = 1.0e10;
    const TINY_EIGVAL: f64 = 1.0e-10;

    let n = mat.num_rows();
    assert!(
        eigen_values.len() >= num_eig && num_eig <= n,
        "size mismatch"
    );
    if let Some(vecs) = &eigen_vectors {
        assert!(vecs.len() >= num_eig * n, "size mismatch");
    }
    let mut seed = i64::from(seed);
    let mut z = vec![0.0; n];
    let mut converged = false;
    for ei in 0..num_eig {
        let mut eig_val = -HUGE_EIGVAL;
        seed += ei as i64;
        let mut v = random_unit_vector(n, seed);
        converged = false;
        for _ in 0..MAX_ITERATIONS {
            mat.multiply(&v, &mut z);
            let prev_val = eig_val;
            // the entry with the largest magnitude (the first one on ties)
            let mut best = -1.0;
            let mut id = n;
            for (i, x) in z.iter().enumerate() {
                if x.abs() > best {
                    best = x.abs();
                    id = i;
                }
            }
            eig_val = z[id];
            if eig_val.abs() < TINY_EIGVAL {
                break;
            }
            v.copy_from_slice(&z);
            for x in &mut v {
                *x /= eig_val;
            }
            if (eig_val - prev_val).abs() < TOLERANCE {
                converged = true;
                break;
            }
        }
        if !converged {
            break;
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in &mut v {
            *x /= norm;
        }
        if let Some(vecs) = eigen_vectors.as_deref_mut() {
            vecs[ei * n..(ei + 1) * n].copy_from_slice(&v);
        }
        eigen_values[ei] = eig_val;
        // remove this eigenvector's space from the matrix
        let data = mat.data_mut();
        for i in 0..n {
            let id = i * (i + 1) / 2;
            for j in 0..=i {
                data[id + j] -= eig_val * v[i] * v[j];
            }
        }
    }
    converged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_the_leading_eigenpair() {
        // the all-ones matrix has the eigenvalue 3 for (1, 1, 1) / sqrt(3)
        let mut m = SymmMatrix::filled(3, 1.0);
        let mut vals = [0.0; 1];
        let mut vecs = [0.0; 3];
        assert!(power_eigen_solver(
            1,
            &mut m,
            &mut vals,
            Some(&mut vecs),
            17
        ));
        assert!((vals[0] - 3.0).abs() < 1e-2, "{vals:?}");
        for v in vecs {
            assert!((v - 1.0 / 3f64.sqrt()).abs() < 1e-2);
        }
        // the deflated matrix has no other eigenvalues
        assert!(m.data().iter().all(|x| x.abs() < 1e-2));
    }

    #[test]
    fn symmetric_multiply() {
        let mut m = SymmMatrix::new(3);
        for (i, j, v) in [
            (0, 0, 1.0),
            (1, 0, 2.0),
            (1, 1, 3.0),
            (2, 0, 4.0),
            (2, 1, 5.0),
            (2, 2, 6.0),
        ] {
            m.set_val(i, j, v);
        }
        let mut y = [0.0; 3];
        m.multiply(&[1.0, 1.0, 1.0], &mut y);
        assert_eq!(y, [7.0, 10.0, 15.0]);
    }
}
