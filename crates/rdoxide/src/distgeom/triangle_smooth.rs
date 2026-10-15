//! Triangle smoothing of distance bounds (upstream `DistGeom/TriangleSmooth.cpp`).

use super::BoundsMatrix;

/// Tightens the bounds with the triangle inequality. Returns false when the
/// bounds turn out to be inconsistent; a lower bound that exceeds its upper
/// bound by a relative amount below `tol` is accepted and the upper bound is
/// raised.
pub fn triangle_smooth_bounds(bounds: &mut BoundsMatrix, tol: f64) -> bool {
    let npt = bounds.num_rows();
    for k in 0..npt {
        for i in 0..npt.saturating_sub(1) {
            if i == k {
                continue;
            }
            let (mut ii, mut ik) = (i, k);
            if ii > ik {
                std::mem::swap(&mut ii, &mut ik);
            }
            let uik = bounds.val(ii, ik);
            let lik = bounds.val(ik, ii);
            for j in i + 1..npt {
                if j == k {
                    continue;
                }
                let (mut jj, mut jk) = (j, k);
                if jj > jk {
                    std::mem::swap(&mut jj, &mut jk);
                }
                let ukj = bounds.val(jj, jk);
                let sum_uik_ukj = uik + ukj;
                if bounds.val(i, j) > sum_uik_ukj {
                    bounds.set_val(i, j, sum_uik_ukj);
                }
                let diff_lik_ujk = lik - ukj;
                let diff_ljk_uik = bounds.val(jk, jj) - uik;
                if bounds.val(j, i) < diff_lik_ujk {
                    bounds.set_val(j, i, diff_lik_ujk);
                } else if bounds.val(j, i) < diff_ljk_uik {
                    bounds.set_val(j, i, diff_ljk_uik);
                }
                let l_bound = bounds.val(j, i);
                let u_bound = bounds.val(i, j);
                if tol > 0.
                    && (l_bound - u_bound) / l_bound > 0.
                    && (l_bound - u_bound) / l_bound < tol
                {
                    bounds.set_val(i, j, l_bound);
                } else if l_bound - u_bound > 0. {
                    return false;
                }
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upper_bounds_follow_the_triangle_inequality() {
        let mut m = BoundsMatrix::new(3);
        for (i, j, l, u) in [(0, 1, 1.0, 1.5), (1, 2, 1.0, 1.5), (0, 2, 0.5, 10.0)] {
            m.set_upper_bound(i, j, u).unwrap();
            m.set_lower_bound(i, j, l).unwrap();
        }
        assert!(triangle_smooth_bounds(&mut m, 0.0));
        assert_eq!(m.upper_bound(0, 2), 3.0);
    }

    #[test]
    fn inconsistent_bounds_fail() {
        let mut m = BoundsMatrix::new(3);
        for (i, j, l, u) in [(0, 1, 1.0, 1.5), (1, 2, 1.0, 1.5), (0, 2, 5.0, 10.0)] {
            m.set_upper_bound(i, j, u).unwrap();
            m.set_lower_bound(i, j, l).unwrap();
        }
        assert!(!triangle_smooth_bounds(&mut m, 0.0));
    }
}
