//! Distance helpers for internal-coordinate geometry (upstream
//! `Geometry/Utils.h`).

use super::{Point3D, Transform3D};

/// The 1-3 distance from the two 1-2 distances and the angle between them.
pub fn compute_13_dist(d1: f64, d2: f64, angle: f64) -> f64 {
    let res = d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * angle.cos();
    res.sqrt()
}

/// The 1-4 distance from the three 1-2 distances, the two bond angles and an
/// arbitrary torsion angle.
///
/// Atom 2 is put at the origin with the 2-3 axis along x and atom 1 in the xy
/// plane with a positive y; atom 4 is then rotated about x by the torsion.
/// For 0 and 180 degrees prefer [`compute_14_dist_cis`] and
/// [`compute_14_dist_trans`], which are faster.
pub fn compute_14_dist_3d(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64, tor_ang: f64) -> f64 {
    let p1 = Point3D::new(d1 * ang12.cos(), d1 * ang12.sin(), 0.0);
    // atom 4 for a torsion of 0
    let p4 = Point3D::new(d2 - d3 * ang23.cos(), d3 * ang23.sin(), 0.0);
    let mut trans = Transform3D::identity();
    trans.set_rotation(tor_ang, &Point3D::new(1.0, 0.0, 0.0));
    let mut p4 = trans.transform_point(&p4);
    p4 -= p1;
    p4.length()
}

/// The 1-4 distance for a cis (0 degree) torsion.
pub fn compute_14_dist_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() - d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

/// The 1-4 distance for a trans (180 degree) torsion.
pub fn compute_14_dist_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() + d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn special_cases_match_general_torsion() {
        let (d1, d2, d3, a12, a23) = (1.5, 1.4, 1.3, 2.0, 1.9);
        let cis = compute_14_dist_3d(d1, d2, d3, a12, a23, 0.0);
        let trans = compute_14_dist_3d(d1, d2, d3, a12, a23, PI);
        assert!((cis - compute_14_dist_cis(d1, d2, d3, a12, a23)).abs() < 1e-12);
        assert!((trans - compute_14_dist_trans(d1, d2, d3, a12, a23)).abs() < 1e-12);
        assert!(cis < trans);
    }

    #[test]
    fn right_angle_13_distance() {
        assert!((compute_13_dist(3.0, 4.0, PI / 2.0) - 5.0).abs() < 1e-12);
    }
}
