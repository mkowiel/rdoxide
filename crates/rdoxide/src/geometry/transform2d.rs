//! 2D affine transforms (`RDGeom::Transform2D`, a 3x3 `SquareMatrix`).
//!
//! Ported from RDKit `Code/Geometry/Transform2D.{h,cpp}` and the
//! multiplication in `Code/Numerics/SquareMatrix.h` (2026.09.1); the
//! arithmetic follows upstream's evaluation order.
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::Point2D;
use crate::rdgeneral::fp::fmuladd;

/// Row-major 3x3 matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform2D {
    pub data: [f64; 9],
}

impl Default for Transform2D {
    fn default() -> Self {
        Transform2D::identity()
    }
}

impl Transform2D {
    pub fn identity() -> Self {
        let mut data = [0.0; 9];
        data[0] = 1.0;
        data[4] = 1.0;
        data[8] = 1.0;
        Transform2D { data }
    }

    /// `setToIdentity`
    pub fn set_to_identity(&mut self) {
        *self = Transform2D::identity();
    }

    /// `getVal(i, j)` / `setVal(i, j, v)`
    pub fn get_val(&self, i: usize, j: usize) -> f64 {
        self.data[i * 3 + j]
    }
    pub fn set_val(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * 3 + j] = v;
    }

    /// `TransformPoint`
    pub fn transform_point(&self, pt: &mut Point2D) {
        let d = &self.data;
        let x = fmuladd(d[0], pt.x, d[1] * pt.y) + d[2];
        let y = fmuladd(d[3], pt.x, d[4] * pt.y) + d[5];
        pt.x = x;
        pt.y = y;
    }

    /// `SetTranslation`
    pub fn set_translation(&mut self, pt: &Point2D) {
        self.data[2] = pt.x;
        self.data[5] = pt.y;
        self.data[8] = 1.0;
    }

    /// `SquareMatrix::operator*=`: `self = self * b`
    pub fn mul_assign(&mut self, b: &Transform2D) {
        let a = self.data;
        let mut n = [0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                let mut v = 0.0;
                for k in 0..3 {
                    v = fmuladd(a[i * 3 + k], b.data[k * 3 + j], v);
                }
                n[i * 3 + j] = v;
            }
        }
        self.data = n;
    }

    /// `SetTransform(pt, angle)`: rotation by `angle` about `pt`.
    pub fn set_transform_rotation(&mut self, pt: &Point2D, angle: f64) {
        self.set_to_identity();
        let mut trans1 = Transform2D::identity();
        trans1.set_translation(&-*pt);
        self.data[0] = angle.cos();
        self.data[1] = -angle.sin();
        self.data[3] = angle.sin();
        self.data[4] = angle.cos();
        self.mul_assign(&trans1);
        let mut trans2 = Transform2D::identity();
        trans2.set_translation(pt);
        trans2.mul_assign(self);
        *self = trans2;
    }

    /// `SetTransform(ref1, ref2, pt1, pt2)`: maps `pt1` onto `ref1` and
    /// rotates the line (pt1, pt2) onto (ref1, ref2).
    pub fn set_transform(&mut self, ref1: &Point2D, ref2: &Point2D, pt1: &Point2D, pt2: &Point2D) {
        let rvec = *ref2 - *ref1;
        let pvec = *pt2 - *pt1;
        let dp = rvec.dot_product(&pvec);
        let lp = rvec.length() * pvec.length();
        if lp <= 0.0 {
            self.set_to_identity();
            return;
        }
        let cval = (dp / lp).clamp(-1.0, 1.0);
        let mut ang = cval.acos();
        let cross = fmuladd(pvec.x, rvec.y, -(pvec.y * rvec.x));
        if cross < 0.0 {
            ang *= -1.0;
        }
        self.set_to_identity();
        self.data[0] = ang.cos();
        self.data[1] = -ang.sin();
        self.data[3] = ang.sin();
        self.data[4] = ang.cos();
        let mut npt1 = *pt1;
        self.transform_point(&mut npt1);
        self.data[2] = ref1.x - npt1.x;
        self.data[5] = ref1.y - npt1.y;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotation_about_point() {
        let mut t = Transform2D::identity();
        t.set_transform_rotation(&Point2D::new(1.0, 0.0), std::f64::consts::PI / 2.0);
        let mut p = Point2D::new(2.0, 0.0);
        t.transform_point(&mut p);
        assert!((p.x - 1.0).abs() < 1e-12 && (p.y - 1.0).abs() < 1e-12);
    }

    #[test]
    fn align_points() {
        let mut t = Transform2D::identity();
        let (r1, r2) = (Point2D::new(1.0, 1.0), Point2D::new(1.0, 3.0));
        let (p1, p2) = (Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0));
        t.set_transform(&r1, &r2, &p1, &p2);
        let mut a = p1;
        let mut b = p2;
        t.transform_point(&mut a);
        t.transform_point(&mut b);
        assert!((a - r1).length() < 1e-12);
        assert!((b - Point2D::new(1.0, 3.0)).length() < 1e-12);
    }
}
