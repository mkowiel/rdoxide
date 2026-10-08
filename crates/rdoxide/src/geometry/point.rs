//! Ported from RDKit `Code/Geometry/point.{h,cpp}` (Point3D / Point2D).
//!
//! Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::error::{RDKitError, Result};
use crate::rdgeneral::fp::fmuladd;

pub const ZERO_TOLERANCE: f64 = 1.0e-16;

/// A point (or vector) in 3D space.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Point3D { x, y, z }
    }

    /// `x * x + y * y + z * z` (contracted like the reference build)
    pub fn length_sq(&self) -> f64 {
        fmuladd(self.z, self.z, fmuladd(self.x, self.x, self.y * self.y))
    }

    pub fn length(&self) -> f64 {
        self.length_sq().sqrt()
    }

    pub fn normalize(&mut self) -> Result<()> {
        let l = self.length();
        if l < ZERO_TOLERANCE {
            return Err(RDKitError::Value(
                "Cannot normalize a zero length vector".into(),
            ));
        }
        self.x /= l;
        self.y /= l;
        self.z /= l;
        Ok(())
    }

    pub fn dot_product(&self, o: &Point3D) -> f64 {
        fmuladd(self.z, o.z, fmuladd(self.x, o.x, self.y * o.y))
    }

    /// Unsigned angle in [0, pi] between the two vectors.
    pub fn angle_to(&self, o: &Point3D) -> f64 {
        let lsq = self.length_sq() * o.length_sq();
        let dot = self.dot_product(o) / lsq.sqrt();
        if dot <= -1.0 {
            return PI;
        }
        if dot >= 1.0 {
            return 0.0;
        }
        dot.acos()
    }

    /// Signed angle in [0, 2pi), sign from the z-component of the cross product.
    pub fn signed_angle_to(&self, o: &Point3D) -> f64 {
        let res = self.angle_to(o);
        if fmuladd(self.x, o.y, -(self.y * o.x)) < -ZERO_TOLERANCE {
            2.0 * PI - res
        } else {
            res
        }
    }

    pub fn direction_vector(&self, o: &Point3D) -> Result<Point3D> {
        let mut res = *o - *self;
        res.normalize()?;
        Ok(res)
    }

    pub fn cross_product(&self, o: &Point3D) -> Point3D {
        Point3D {
            x: fmuladd(self.y, o.z, -(self.z * o.y)),
            y: fmuladd(-self.x, o.z, self.z * o.x),
            z: fmuladd(self.x, o.y, -(self.y * o.x)),
        }
    }

    pub fn get_perpendicular(&self) -> Result<Point3D> {
        let mut res = Point3D::default();
        if self.x != 0.0 {
            if self.y != 0.0 {
                res.y = -self.x;
                res.x = self.y;
            } else if self.z != 0.0 {
                res.z = -self.x;
                res.x = self.z;
            } else {
                res.y = 1.0;
            }
        } else if self.y != 0.0 {
            if self.z != 0.0 {
                res.z = -self.y;
                res.y = self.z;
            } else {
                res.x = 1.0;
            }
        } else if self.z != 0.0 {
            res.x = 1.0;
        }
        res.normalize()?;
        Ok(res)
    }

    pub fn get(&self, i: usize) -> Option<f64> {
        match i {
            0 => Some(self.x),
            1 => Some(self.y),
            2 => Some(self.z),
            _ => None,
        }
    }
}

impl Add for Point3D {
    type Output = Point3D;
    fn add(self, o: Point3D) -> Point3D {
        Point3D::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl Sub for Point3D {
    type Output = Point3D;
    fn sub(self, o: Point3D) -> Point3D {
        Point3D::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl Neg for Point3D {
    type Output = Point3D;
    fn neg(self) -> Point3D {
        Point3D::new(-self.x, -self.y, -self.z)
    }
}
impl Mul<f64> for Point3D {
    type Output = Point3D;
    fn mul(self, s: f64) -> Point3D {
        Point3D::new(self.x * s, self.y * s, self.z * s)
    }
}
impl Div<f64> for Point3D {
    type Output = Point3D;
    fn div(self, s: f64) -> Point3D {
        Point3D::new(self.x / s, self.y / s, self.z / s)
    }
}
impl AddAssign for Point3D {
    fn add_assign(&mut self, o: Point3D) {
        *self = *self + o;
    }
}
impl SubAssign for Point3D {
    fn sub_assign(&mut self, o: Point3D) {
        *self = *self - o;
    }
}
impl MulAssign<f64> for Point3D {
    fn mul_assign(&mut self, s: f64) {
        *self = *self * s;
    }
}
impl DivAssign<f64> for Point3D {
    fn div_assign(&mut self, s: f64) {
        *self = *self / s;
    }
}

/// A point in 2D space.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

impl Point2D {
    pub const fn new(x: f64, y: f64) -> Self {
        Point2D { x, y }
    }
    pub fn length(&self) -> f64 {
        self.length_sq().sqrt()
    }
    pub fn length_sq(&self) -> f64 {
        fmuladd(self.x, self.x, self.y * self.y)
    }
    pub fn dot_product(&self, o: &Point2D) -> f64 {
        fmuladd(self.x, o.x, self.y * o.y)
    }
    pub fn normalize(&mut self) -> Result<()> {
        let l = self.length();
        if l < ZERO_TOLERANCE {
            return Err(RDKitError::Runtime(
                "Cannot normalize a zero length vector".into(),
            ));
        }
        self.x /= l;
        self.y /= l;
        Ok(())
    }
    pub fn rotate90(&mut self) {
        let t = self.x;
        self.x = -self.y;
        self.y = t;
    }
    /// Unsigned angle in [0, pi] between the two vectors.
    pub fn angle_to(&self, o: &Point2D) -> Result<f64> {
        let mut t1 = *self;
        let mut t2 = *o;
        t1.normalize()?;
        t2.normalize()?;
        let d = t1.dot_product(&t2).clamp(-1.0, 1.0);
        Ok(d.acos())
    }
    /// Angle in [0, 2pi), sign from the cross product.
    pub fn signed_angle_to(&self, o: &Point2D) -> Result<f64> {
        let mut res = self.angle_to(o)?;
        if fmuladd(self.x, o.y, -(self.y * o.x)) < -ZERO_TOLERANCE {
            res = 2.0 * PI - res;
        }
        Ok(res)
    }
    pub fn direction_vector(&self, o: &Point2D) -> Result<Point2D> {
        let mut res = Point2D::new(o.x - self.x, o.y - self.y);
        res.normalize()?;
        Ok(res)
    }
}

impl Add for Point2D {
    type Output = Point2D;
    fn add(self, o: Point2D) -> Point2D {
        Point2D::new(self.x + o.x, self.y + o.y)
    }
}
impl Sub for Point2D {
    type Output = Point2D;
    fn sub(self, o: Point2D) -> Point2D {
        Point2D::new(self.x - o.x, self.y - o.y)
    }
}
impl Neg for Point2D {
    type Output = Point2D;
    fn neg(self) -> Point2D {
        Point2D::new(-self.x, -self.y)
    }
}
impl Mul<f64> for Point2D {
    type Output = Point2D;
    fn mul(self, s: f64) -> Point2D {
        Point2D::new(self.x * s, self.y * s)
    }
}
impl AddAssign for Point2D {
    fn add_assign(&mut self, o: Point2D) {
        self.x += o.x;
        self.y += o.y;
    }
}
impl SubAssign for Point2D {
    fn sub_assign(&mut self, o: Point2D) {
        self.x -= o.x;
        self.y -= o.y;
    }
}
impl MulAssign<f64> for Point2D {
    fn mul_assign(&mut self, s: f64) {
        self.x *= s;
        self.y *= s;
    }
}

/// `RDGeom::computeDihedralAngle` (point.cpp): unsigned dihedral, in [0, PI].
pub fn compute_dihedral_angle(pt1: &Point3D, pt2: &Point3D, pt3: &Point3D, pt4: &Point3D) -> f64 {
    let beg_end = *pt3 - *pt2;
    let beg_nbr = *pt1 - *pt2;
    let crs1 = beg_nbr.cross_product(&beg_end);
    let end_nbr = *pt4 - *pt3;
    let crs2 = end_nbr.cross_product(&beg_end);
    crs1.angle_to(&crs2)
}

/// `RDGeom::computeSignedDihedralAngle` (point.cpp): in [-PI, PI].
pub fn compute_signed_dihedral_angle(
    pt1: &Point3D,
    pt2: &Point3D,
    pt3: &Point3D,
    pt4: &Point3D,
) -> f64 {
    let beg_end = *pt3 - *pt2;
    let beg_nbr = *pt1 - *pt2;
    let crs1 = beg_nbr.cross_product(&beg_end);
    let end_nbr = *pt4 - *pt3;
    let crs2 = end_nbr.cross_product(&beg_end);
    let mut ang = crs1.angle_to(&crs2);
    if crs1.cross_product(&crs2).dot_product(&beg_end) < 0.0 {
        ang *= -1.0;
    }
    ang
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        let a = Point3D::new(1.0, 0.0, 0.0);
        let b = Point3D::new(0.0, 1.0, 0.0);
        assert_eq!(a.cross_product(&b), Point3D::new(0.0, 0.0, 1.0));
        assert!((a.angle_to(&b) - PI / 2.0).abs() < 1e-12);
        assert!((b.signed_angle_to(&a) - 1.5 * PI).abs() < 1e-12);
        assert_eq!((a + b) * 2.0, Point3D::new(2.0, 2.0, 0.0));
        assert!(Point3D::default().normalize().is_err());
        assert_eq!(a.get_perpendicular().unwrap(), Point3D::new(0.0, 1.0, 0.0));
    }
}
