//! Embedding helpers (upstream `DistGeom/DistGeomUtils.cpp`): picking a
//! random distance matrix, the initial coordinates from it, and the
//! distance-geometry error function.

use std::collections::BTreeMap;

use super::{
    BoundsMatrix, ChiralSet, ChiralViolationContribs, DistViolationContribs, FourthDimContribs,
};
use crate::error::{RDKitError, Result};
use crate::forcefield::ForceField;
use crate::numerics::power_eigen::power_eigen_solver;
use crate::numerics::symm_matrix::SymmMatrix;
use crate::rdgeneral::random::MinStdRand;

const EIGVAL_TOL: f64 = 0.001;

/// Fills `dist_mat` with random distances between the lower and upper
/// bounds and returns the largest one.
pub fn pick_random_dist_mat(
    mmat: &BoundsMatrix,
    dist_mat: &mut SymmMatrix,
    rng: &mut MinStdRand,
) -> Result<f64> {
    let npt = mmat.num_rows();
    if npt != dist_mat.num_rows() {
        return Err(RDKitError::invariant("Size mismatch"));
    }
    let mut largest = -1.0;
    for i in 1..npt {
        for j in 0..i {
            let ub = mmat.upper_bound(i, j);
            let lb = mmat.lower_bound(i, j);
            if ub < lb {
                return Err(RDKitError::invariant(""));
            }
            let d = lb + rng.next_f64() * (ub - lb);
            dist_mat.set_val(i, j, d);
            if d > largest {
                largest = d;
            }
        }
    }
    Ok(largest)
}

/// Computes coordinates (`dim` values per point, flat) that reproduce the
/// distances as well as possible by metrization. Fails (false) when the
/// matrix has too many (near-)zero eigenvalues or, unless `rand_neg_eig`,
/// negative ones; coordinates for negative eigenvalues are then random.
pub fn compute_initial_coords(
    dist_mat: &SymmMatrix,
    positions: &mut [f64],
    dim: usize,
    rng: &mut MinStdRand,
    rand_neg_eig: bool,
    num_zero_fail: usize,
) -> Result<bool> {
    let n = dist_mat.num_rows();
    if positions.len() != n * dim {
        return Err(RDKitError::precondition("Size mismatch"));
    }
    let mut sq_mat = SymmMatrix::new(n);
    let mut sum_sq_d2 = 0.0;
    for (s, d) in sq_mat.data_mut().iter_mut().zip(dist_mat.data()) {
        *s = d * d;
        sum_sq_d2 += *s;
    }
    sum_sq_d2 /= (n * n) as f64;

    let mut sq_d0i = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            sq_d0i[i] += sq_mat.val(i, j);
        }
        sq_d0i[i] /= n as f64;
        sq_d0i[i] -= sum_sq_d2;
        if sq_d0i[i] < EIGVAL_TOL && n > 3 {
            return Ok(false);
        }
    }

    let mut t = SymmMatrix::new(n);
    for i in 0..n {
        for j in 0..=i {
            t.set_val(i, j, 0.5 * (sq_d0i[i] + sq_d0i[j] - sq_mat.val(i, j)));
        }
    }
    let n_eigs = dim.min(n);
    let mut eig_vals = vec![0.0; dim];
    let mut eig_vecs = vec![0.0; dim * n];
    power_eigen_solver(
        n_eigs,
        &mut t,
        &mut eig_vals,
        Some(&mut eig_vecs),
        (sum_sq_d2 * n as f64) as i32,
    );

    let mut found_neg = false;
    let mut zero_eigs = 0;
    for e in &mut eig_vals {
        if *e > EIGVAL_TOL {
            *e = e.sqrt();
        } else if e.abs() < EIGVAL_TOL {
            *e = 0.0;
            zero_eigs += 1;
        } else {
            found_neg = true;
        }
    }
    if found_neg && !rand_neg_eig {
        return Ok(false);
    }
    if zero_eigs >= num_zero_fail && n > 3 {
        return Ok(false);
    }
    for i in 0..n {
        for j in 0..dim {
            positions[i * dim + j] = if eig_vals[j] >= 0.0 {
                eig_vals[j] * eig_vecs[j * n + i]
            } else {
                1.0 - 2.0 * rng.next_f64()
            };
        }
    }
    Ok(true)
}

/// Random coordinates in a cube of side `box_size` centered on the origin.
pub fn compute_random_coords(
    positions: &mut [f64],
    box_size: f64,
    rng: &mut MinStdRand,
) -> Result<bool> {
    if box_size <= 0.0 {
        return Err(RDKitError::precondition("bad boxSize"));
    }
    for p in positions {
        *p = box_size * (rng.next_f64() - 0.5);
    }
    Ok(true)
}

/// Builds the error function of the embedding: distance violations (for
/// the pairs with a basin no wider than `basin_size_tol` or an entry in
/// `extra_weights`), chiral volume violations and the penalty of the fourth
/// dimension. `fixed_pts` marks points whose mutual distances are skipped.
#[allow(clippy::too_many_arguments)]
pub fn construct_force_field(
    mmat: &BoundsMatrix,
    positions: &[f64],
    dim: usize,
    csets: &[ChiralSet],
    weight_distance: f64,
    weight_chiral: f64,
    weight_fourth_dim: f64,
    extra_weights: Option<&BTreeMap<(usize, usize), f64>>,
    basin_size_tol: f64,
    fixed_pts: Option<&[bool]>,
) -> Result<ForceField> {
    let n = mmat.num_rows();
    if positions.len() != n * dim {
        return Err(RDKitError::precondition("bad positions"));
    }
    let mut field = ForceField::new(dim);
    for p in positions.chunks_exact(dim) {
        field.push_point(p);
    }
    let mut contrib = DistViolationContribs::new(dim);
    for i in 1..n {
        for j in 0..i {
            if fixed_pts.is_some_and(|f| f[i] && f[j]) {
                continue;
            }
            let mut w = weight_distance;
            let l = mmat.lower_bound(i, j);
            let u = mmat.upper_bound(i, j);
            let mut include_it = false;
            if let Some(&extra) = extra_weights.and_then(|e| e.get(&(i, j))) {
                w = extra;
                include_it = true;
            }
            if u - l <= basin_size_tol {
                include_it = true;
            }
            if include_it {
                contrib.add_contrib(i, j, u, l, w);
            }
        }
    }
    if !contrib.is_empty() {
        field.push_contrib(Box::new(contrib));
    }
    if weight_chiral > 1.0e-8 {
        let mut contrib = ChiralViolationContribs::new(dim);
        for cset in csets {
            contrib.add_contrib(cset, weight_chiral);
        }
        if !contrib.is_empty() {
            field.push_contrib(Box::new(contrib));
        }
    }
    if dim == 4 && weight_fourth_dim > 1.0e-8 {
        let mut contrib = FourthDimContribs::new();
        for i in 0..n {
            contrib.add_contrib(i, weight_fourth_dim);
        }
        if !contrib.is_empty() {
            field.push_contrib(Box::new(contrib));
        }
    }
    Ok(field)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// bounds of a regular tetrahedron-ish cluster: all distances in [0.9, 1.1]
    fn tight_bounds(n: usize) -> BoundsMatrix {
        let mut m = BoundsMatrix::new(n);
        for i in 1..n {
            for j in 0..i {
                m.set_upper_bound(i, j, 1.1).unwrap();
                m.set_lower_bound(i, j, 0.9).unwrap();
            }
        }
        m
    }

    #[test]
    fn embedding_satisfies_the_bounds() {
        let n = 5;
        let bounds = tight_bounds(n);
        let mut rng = MinStdRand::new(42);
        let mut dist = SymmMatrix::new(n);
        let largest = pick_random_dist_mat(&bounds, &mut dist, &mut rng).unwrap();
        assert!((0.9..=1.1).contains(&largest));
        let mut pos = vec![0.0; n * 4];
        // a few seeds may give negative eigenvalues; random fill is allowed
        assert!(compute_initial_coords(&dist, &mut pos, 4, &mut rng, true, 4).unwrap());
        let mut ff =
            construct_force_field(&bounds, &pos, 4, &[], 1.0, 0.0, 0.1, None, 5.0, None).unwrap();
        let before = ff.calc_energy();
        assert_eq!(ff.minimize(500, 1e-4, 1e-6).unwrap(), 0);
        assert!(ff.calc_energy() <= before);
        assert!(ff.calc_energy() < 0.5);
    }

    #[test]
    fn random_coordinates_stay_in_the_box() {
        let mut rng = MinStdRand::new(7);
        let mut pos = vec![0.0; 30];
        compute_random_coords(&mut pos, 10.0, &mut rng).unwrap();
        assert!(pos.iter().all(|x| x.abs() <= 5.0));
        assert!(pos.iter().any(|x| *x != 0.0));
        assert!(compute_random_coords(&mut pos, 0.0, &mut rng).is_err());
    }
}
