//! `RDGeom::Transform3D` (upstream `Geometry/Transform3D.{h,cpp}`): the parts
//! needed so far (rotation about an axis, point transformation).
//!
//! Products are fused like clang's `-ffp-contract=on` does for upstream's
//! expressions (see `rdgeneral::fp`).

use super::point::Point3D;
use crate::rdgeneral::fp::fmuladd;

/// A 4x4 row-major transformation matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct Transform3D {
    pub data: [f64; 16],
}

impl Default for Transform3D {
    fn default() -> Self {
        Self::identity()
    }
}

impl Transform3D {
    pub fn identity() -> Self {
        let mut data = [0.0; 16];
        data[0] = 1.0;
        data[5] = 1.0;
        data[10] = 1.0;
        data[15] = 1.0;
        Transform3D { data }
    }

    /// `getVal(row, col)`
    pub fn get_val(&self, row: usize, col: usize) -> f64 {
        self.data[row * 4 + col]
    }

    /// `setVal(row, col, val)`
    pub fn set_val(&mut self, row: usize, col: usize, val: f64) {
        self.data[row * 4 + col] = val;
    }

    /// `SetTranslation(move)`: the last column of the matrix
    pub fn set_translation(&mut self, translation: &Point3D) {
        self.data[3] = translation.x;
        self.data[7] = translation.y;
        self.data[11] = translation.z;
        self.data[15] = 1.0;
    }

    /// `SetRotation(angle, axis)`: resets to the identity, then a rotation by
    /// `angle` (radians) about `axis` (assumed normalized)
    pub fn set_rotation(&mut self, angle: f64, axis: &Point3D) {
        *self = Self::identity();
        self.set_rotation_cs(angle.cos(), angle.sin(), axis);
    }

    /// `SetRotation(cosT, sinT, axis)`
    pub fn set_rotation_cs(&mut self, cos_t: f64, sin_t: f64, axis: &Point3D) {
        let t = 1.0 - cos_t;
        let (x, y, z) = (axis.x, axis.y, axis.z);
        let d = &mut self.data;
        d[0] = fmuladd(t * x, x, cos_t);
        d[1] = fmuladd(t * x, y, -(sin_t * z));
        d[2] = fmuladd(t * x, z, sin_t * y);

        d[4] = fmuladd(t * x, y, sin_t * z);
        d[5] = fmuladd(t * y, y, cos_t);
        d[6] = fmuladd(t * y, z, -(sin_t * x));

        d[8] = fmuladd(t * x, z, -(sin_t * y));
        d[9] = fmuladd(t * y, z, sin_t * x);
        d[10] = fmuladd(t * z, z, cos_t);
    }

    /// `operator*=(other)`: this becomes the combination of `other` and the
    /// current state, i.e. `new(point) = old(other(point))`.
    pub fn mul_assign(&mut self, other: &Transform3D) {
        let mut res = [0.0; 16];
        for i in 0..4 {
            for j in 0..4 {
                res[i * 4 + j] = (0..4).fold(0.0, |acc, k| {
                    fmuladd(self.data[i * 4 + k], other.data[k * 4 + j], acc)
                });
            }
        }
        self.data = res;
    }

    /// `TransformPoint` / `operator*(Transform3D, Point3D)`
    pub fn transform_point(&self, p: &Point3D) -> Point3D {
        let d = &self.data;
        let row = |o: usize| fmuladd(d[o + 2], p.z, fmuladd(d[o], p.x, d[o + 1] * p.y)) + d[o + 3];
        Point3D {
            x: row(0),
            y: row(4),
            z: row(8),
        }
    }
}
