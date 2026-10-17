//! BFGS minimizer with a backtracking line search (upstream
//! `Numerics/Optimizer/BFGSOpt.h`, after Numerical Recipes in C, sections 9.7
//! and 10.7). The scalar code path is used; upstream's SVE path computes the
//! same sums in another order.

use crate::error::{RDKitError, Result};

/// Default tolerance for function convergence in the minimizer
pub const FUNCTOL: f64 = 1e-4;
/// Default tolerance for x changes in the minimizer
pub const MOVETOL: f64 = 1e-7;
/// Default maximum number of iterations
pub const MAXITS: u32 = 200;
/// Default gradient tolerance in the minimizer
pub const EPS: f64 = 3e-8;
/// Default direction vector tolerance in the minimizer
pub const TOLX: f64 = 4. * EPS;
/// Default maximum step size in the minimizer
pub const MAXSTEP: f64 = 100.0;

const MAX_ITER_LINEAR_SEARCH: u32 = 1000;

/// A position and energy recorded during a minimization.
#[derive(Debug, Clone, PartialEq)]
pub struct Snapshot {
    pub coords: Vec<f64>,
    pub energy: f64,
}

/// Outcome of a line search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LineSearch {
    /// the step was accepted
    Success,
    /// the step size got too small, which probably indicates success
    TinyStep,
    /// the direction is bad (orthogonal to the gradient)
    BadDirection,
}

/// Backtracking line search along `dir` from `old_pt`. `new_pt` receives the
/// final position and the new function value is returned.
fn linear_search(
    old_pt: &[f64],
    old_val: f64,
    grad: &[f64],
    dir: &mut [f64],
    new_pt: &mut [f64],
    new_val: &mut f64,
    func: &mut impl FnMut(&[f64]) -> f64,
    max_step: f64,
) -> LineSearch {
    let dim = old_pt.len();
    let mut sum = 0.0;
    for d in dir.iter() {
        sum += d * d;
    }
    let sum = sum.sqrt();
    // rescale if we're trying to move too far
    if sum > max_step {
        for d in dir.iter_mut() {
            *d *= max_step / sum;
        }
    }
    // the direction must have some component along -grad
    let mut slope = 0.0;
    for i in 0..dim {
        slope += dir[i] * grad[i];
    }
    if slope >= 0.0 {
        return LineSearch::BadDirection;
    }
    let mut test = 0.0f64;
    for i in 0..dim {
        let temp = dir[i].abs() / old_pt[i].abs().max(1.0);
        if temp > test {
            test = temp;
        }
    }
    let lambda_min = MOVETOL / test;
    let mut lambda = 1.0f64;
    let (mut lambda2, mut val2) = (0.0f64, 0.0f64);
    let mut it = 0;
    while it < MAX_ITER_LINEAR_SEARCH {
        if lambda < lambda_min {
            // step size is below the position-scaled threshold: converged;
            // upstream still falls through to restoring the old position
            new_pt.copy_from_slice(old_pt);
            return LineSearch::TinyStep;
        }
        for i in 0..dim {
            new_pt[i] = old_pt[i] + lambda * dir[i];
        }
        *new_val = func(new_pt);
        if *new_val - old_val <= FUNCTOL * lambda * slope {
            // sufficient decrease: accept the step
            return LineSearch::Success;
        }
        // backtrack
        let mut tmp_lambda;
        if it == 0 {
            // quadratic model: only one prior function value is available
            tmp_lambda = -slope / (2.0 * (*new_val - old_val - slope));
        } else {
            let rhs1 = *new_val - old_val - lambda * slope;
            let rhs2 = val2 - old_val - lambda2 * slope;
            let a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
            let b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2))
                / (lambda - lambda2);
            if a == 0.0 {
                tmp_lambda = -slope / (2.0 * b);
            } else {
                let disc = b * b - 3.0 * a * slope;
                if disc < 0.0 {
                    tmp_lambda = 0.5 * lambda;
                } else if b <= 0.0 {
                    tmp_lambda = (-b + disc.sqrt()) / (3.0 * a);
                } else {
                    tmp_lambda = -slope / (b + disc.sqrt());
                }
            }
            if tmp_lambda > 0.5 * lambda {
                tmp_lambda = 0.5 * lambda;
            }
        }
        lambda2 = lambda;
        val2 = *new_val;
        lambda = tmp_lambda.max(0.1 * lambda);
        it += 1;
    }
    // nothing was done
    new_pt.copy_from_slice(old_pt);
    LineSearch::BadDirection
}

/// `out[i] = sum_j mat[i][j] * v[j]`, or `out[i] = 0 - sum_j mat[i][j] * v[j]`
/// when `NEG`, accumulating each row in index order exactly like the scalar
/// loop. Four rows are processed at a time: their sums are independent, so the
/// processor overlaps the dependent additions, and every result stays
/// bit-identical.
fn mat_vec<const NEG: bool>(mat: &[f64], dim: usize, v: &[f64], out: &mut [f64]) {
    #[inline(always)]
    fn step<const NEG: bool>(acc: f64, p: f64) -> f64 {
        if NEG { acc - p } else { acc + p }
    }
    let mut i = 0;
    while i + 4 <= dim {
        let r0 = &mat[i * dim..(i + 1) * dim];
        let r1 = &mat[(i + 1) * dim..(i + 2) * dim];
        let r2 = &mat[(i + 2) * dim..(i + 3) * dim];
        let r3 = &mat[(i + 3) * dim..(i + 4) * dim];
        let (mut a0, mut a1, mut a2, mut a3) = (0.0, 0.0, 0.0, 0.0);
        for j in 0..dim {
            let x = v[j];
            a0 = step::<NEG>(a0, r0[j] * x);
            a1 = step::<NEG>(a1, r1[j] * x);
            a2 = step::<NEG>(a2, r2[j] * x);
            a3 = step::<NEG>(a3, r3[j] * x);
        }
        out[i] = a0;
        out[i + 1] = a1;
        out[i + 2] = a2;
        out[i + 3] = a3;
        i += 4;
    }
    while i < dim {
        let row = &mat[i * dim..(i + 1) * dim];
        let mut a = 0.0;
        for j in 0..dim {
            a = step::<NEG>(a, row[j] * v[j]);
        }
        out[i] = a;
        i += 1;
    }
}

/// Copies the upper triangle of the square matrix onto its lower triangle,
/// tile by tile so that the strided writes stay in cache.
fn mirror_upper_to_lower(mat: &mut [f64], dim: usize) {
    const TILE: usize = 16;
    for bi in (0..dim).step_by(TILE) {
        for bj in (bi..dim).step_by(TILE) {
            for i in bi..(bi + TILE).min(dim) {
                for j in bj.max(i + 1)..(bj + TILE).min(dim) {
                    mat[j * dim + i] = mat[i * dim + j];
                }
            }
        }
    }
}

/// Result of [`minimize`].
#[derive(Debug, Clone, PartialEq)]
pub struct Minimization {
    /// false if too many iterations were required
    pub converged: bool,
    pub num_iters: u32,
    pub func_val: f64,
    pub snapshots: Vec<Snapshot>,
}

/// Minimizes `func` starting from (and updating) `pos`.
///
/// `grad_func` fills the gradient at a position and returns the gradient
/// scale. A snapshot is recorded every `snapshot_freq` iterations (0: none).
pub fn minimize(
    pos: &mut [f64],
    grad_tol: f64,
    mut func: impl FnMut(&[f64]) -> f64,
    mut grad_func: impl FnMut(&[f64], &mut [f64]) -> f64,
    snapshot_freq: u32,
    max_its: u32,
) -> Result<Minimization> {
    if grad_tol <= 0.0 {
        return Err(RDKitError::precondition("bad tolerance"));
    }
    let dim = pos.len();
    let snapshot_freq = snapshot_freq.min(max_its);
    let mut snapshots = Vec::new();
    let mut grad = vec![0.0; dim];
    let mut d_grad = vec![0.0; dim];
    let mut hess_d_grad = vec![0.0; dim];
    let mut xi = vec![0.0; dim];
    let mut inv_hessian = vec![0.0; dim * dim];
    let mut new_pos = vec![0.0; dim];

    let mut fp = func(pos);
    grad_func(pos, &mut grad);

    // the inverse Hessian starts as the identity, the first direction is the
    // steepest descent step
    let mut sum = 0.0;
    for i in 0..dim {
        inv_hessian[i * dim + i] = 1.0;
        xi[i] = -grad[i];
        sum += pos[i] * pos[i];
    }
    let max_step = MAXSTEP * sum.sqrt().max(dim as f64);

    // like upstream, the value of an iteration whose line search takes no
    // step is the one of the previous search
    let mut func_val = fp;
    let mut num_iters = 0;
    for iter in 1..=max_its {
        num_iters = iter;
        let status = linear_search(
            pos,
            fp,
            &grad,
            &mut xi,
            &mut new_pos,
            &mut func_val,
            &mut func,
            max_step,
        );
        if status == LineSearch::BadDirection {
            return Err(RDKitError::invariant("bad direction in linearSearch"));
        }
        // save the function value for the next search
        fp = func_val;
        // set the direction of this line and save the gradient
        let mut test = 0.0f64;
        for i in 0..dim {
            xi[i] = new_pos[i] - pos[i];
            pos[i] = new_pos[i];
            let temp = xi[i].abs() / pos[i].abs().max(1.0);
            if temp > test {
                test = temp;
            }
            d_grad[i] = grad[i];
        }
        if test < TOLX {
            if snapshot_freq != 0 {
                snapshots.push(Snapshot {
                    coords: new_pos.clone(),
                    energy: fp,
                });
            }
            return Ok(Minimization {
                converged: true,
                num_iters,
                func_val,
                snapshots,
            });
        }

        // update the gradient
        let grad_scale = grad_func(pos, &mut grad);
        let mut test = 0.0f64;
        // |funcVal| keeps negative energies from clamping the denominator
        let term = (func_val.abs() * grad_scale).max(1.0);
        for i in 0..dim {
            let temp = grad[i].abs() * pos[i].abs().max(1.0);
            test = test.max(temp);
            d_grad[i] = grad[i] - d_grad[i];
        }
        test /= term;
        if test < grad_tol {
            if snapshot_freq != 0 {
                snapshots.push(Snapshot {
                    coords: new_pos.clone(),
                    energy: fp,
                });
            }
            return Ok(Minimization {
                converged: true,
                num_iters,
                func_val,
                snapshots,
            });
        }

        // BFGS inverse Hessian update
        let (mut fac, mut fae, mut sum_d_grad, mut sum_xi) = (0.0, 0.0, 0.0, 0.0);
        mat_vec::<false>(&inv_hessian, dim, &d_grad, &mut hess_d_grad);
        for i in 0..dim {
            fac += d_grad[i] * xi[i];
            fae += d_grad[i] * hess_d_grad[i];
            sum_d_grad += d_grad[i] * d_grad[i];
            sum_xi += xi[i] * xi[i];
        }
        if fac > (EPS * sum_d_grad * sum_xi).sqrt() {
            fac = 1.0 / fac;
            let fad = 1.0 / fae;
            for i in 0..dim {
                d_grad[i] = fac * xi[i] - fad * hess_d_grad[i];
            }
            // upper triangle (contiguous, so it vectorizes), then mirrored
            for i in 0..dim {
                let pxi = fac * xi[i];
                let hdgi = fad * hess_d_grad[i];
                let dgi = fae * d_grad[i];
                let row = &mut inv_hessian[i * dim + i..(i + 1) * dim];
                for (k, h) in row.iter_mut().enumerate() {
                    let j = i + k;
                    *h = *h + pxi * xi[j] - hdgi * hess_d_grad[j] + dgi * d_grad[j];
                }
            }
            mirror_upper_to_lower(&mut inv_hessian, dim);
        }
        mat_vec::<true>(&inv_hessian, dim, &grad, &mut xi);
        if snapshot_freq != 0 && iter % snapshot_freq == 0 {
            snapshots.push(Snapshot {
                coords: new_pos.clone(),
                energy: fp,
            });
        }
    }
    Ok(Minimization {
        converged: false,
        num_iters,
        func_val,
        snapshots,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimizes_a_quadratic() {
        let mut pos = [3.0, -2.0, 5.0];
        let f = |p: &[f64]| p.iter().map(|x| (x - 1.0) * (x - 1.0)).sum::<f64>();
        let g = |p: &[f64], out: &mut [f64]| {
            for i in 0..p.len() {
                out[i] = 2.0 * (p[i] - 1.0);
            }
            1.0
        };
        let r = minimize(&mut pos, 1e-5, f, g, 0, MAXITS).unwrap();
        assert!(r.converged);
        assert!(pos.iter().all(|x| (x - 1.0).abs() < 1e-6));
    }

    #[test]
    fn rosenbrock() {
        let mut pos = [-1.2, 1.0];
        let f = |p: &[f64]| 100.0 * (p[1] - p[0] * p[0]).powi(2) + (1.0 - p[0]).powi(2);
        let g = |p: &[f64], out: &mut [f64]| {
            out[0] = -400.0 * p[0] * (p[1] - p[0] * p[0]) - 2.0 * (1.0 - p[0]);
            out[1] = 200.0 * (p[1] - p[0] * p[0]);
            1.0
        };
        let r = minimize(&mut pos, 1e-8, f, g, 0, 1000).unwrap();
        assert!(r.converged);
        assert!((pos[0] - 1.0).abs() < 1e-4 && (pos[1] - 1.0).abs() < 1e-4);
    }
}
