//! `sketcherMinimizerMaths.h`: 2D points and the geometry helpers of the sketcher. All arithmetic follows the C++
//! `float` semantics (f32 storage and operators, `double` where a literal or `M_PI` is involved).

// the if/else chains keep the C++ comparisons (and their NaN behaviour) literally
#![allow(clippy::manual_clamp)]

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// `MACROCYCLE`, `BONDLENGTH`, `FRACTION_OF_BONDLENGTH_FOR_CLASH`, `SKETCHER_EPSILON`
pub const MACROCYCLE: usize = 9;
pub const BONDLENGTH: f32 = 50.0;
pub const FRACTION_OF_BONDLENGTH_FOR_CLASH: f64 = 0.25;
pub const SKETCHER_EPSILON: f32 = 0.0001;
const M_PI: f64 = std::f64::consts::PI;

/// `roundToTwoDecimalDigits`
pub fn round_to_two_decimal_digits(f: f32) -> f32 {
    (((f * 100.0) as f64 + 0.5).floor() * 0.01) as f32
}

/// `roundToPrecision`
pub fn round_to_precision(f: f32, precision: i32) -> f32 {
    let scaled = f64::from(f) * 10f64.powi(precision);
    ((scaled + 0.5).floor() * 0.1f64.powi(precision)) as f32
}

/// `sketcherMinimizerPointF`
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PointF {
    pub x: f32,
    pub y: f32,
}

impl PointF {
    pub const fn new(x: f32, y: f32) -> Self {
        PointF { x, y }
    }

    pub fn square_length(&self) -> f32 {
        self.x * self.x + self.y * self.y
    }

    /// 0 for very short vectors
    pub fn length(&self) -> f32 {
        let dd = self.square_length();
        if dd > SKETCHER_EPSILON {
            dd.sqrt()
        } else {
            0.0
        }
    }

    pub fn normalize(&mut self) {
        let q = self.length();
        if q > SKETCHER_EPSILON {
            self.x /= q;
            self.y /= q;
        }
    }

    /// rotates by the angle with sine `s` and cosine `c`
    pub fn rotate(&mut self, s: f32, c: f32) {
        let (x, y) = (self.x, self.y);
        self.x = x * c + y * s;
        self.y = -x * s + y * c;
    }

    pub fn parallel_component(&self, axis: PointF) -> PointF {
        let dot = self.x * axis.x + self.y * axis.y;
        axis * dot / axis.square_length()
    }

    pub fn round(&mut self, precision: i32) {
        if precision == 2 {
            self.x = round_to_two_decimal_digits(self.x);
            self.y = round_to_two_decimal_digits(self.y);
        } else {
            self.x = round_to_precision(self.x, precision);
            self.y = round_to_precision(self.y, precision);
        }
    }
}

impl Add for PointF {
    type Output = PointF;
    fn add(self, o: PointF) -> PointF {
        PointF::new(self.x + o.x, self.y + o.y)
    }
}
impl Sub for PointF {
    type Output = PointF;
    fn sub(self, o: PointF) -> PointF {
        PointF::new(self.x - o.x, self.y - o.y)
    }
}
impl AddAssign for PointF {
    fn add_assign(&mut self, o: PointF) {
        self.x += o.x;
        self.y += o.y;
    }
}
impl SubAssign for PointF {
    fn sub_assign(&mut self, o: PointF) {
        self.x -= o.x;
        self.y -= o.y;
    }
}
impl Mul<f32> for PointF {
    type Output = PointF;
    fn mul(self, c: f32) -> PointF {
        PointF::new(self.x * c, self.y * c)
    }
}
impl Mul<PointF> for f32 {
    type Output = PointF;
    fn mul(self, p: PointF) -> PointF {
        PointF::new(p.x * self, p.y * self)
    }
}
impl Div<f32> for PointF {
    type Output = PointF;
    fn div(self, c: f32) -> PointF {
        PointF::new(self.x / c, self.y / c)
    }
}
impl Neg for PointF {
    type Output = PointF;
    fn neg(self) -> PointF {
        PointF::new(-self.x, -self.y)
    }
}

pub fn dot_product(a: PointF, b: PointF) -> f32 {
    a.x * b.x + a.y * b.y
}

pub fn cross_product(a: PointF, b: PointF) -> f32 {
    a.x * b.y - a.y * b.x
}

/// `intersectionOfSegments`: the intersection point if the segments cross
pub fn intersection_of_segments(
    s1p1: PointF,
    s1p2: PointF,
    s2p1: PointF,
    s2p2: PointF,
) -> Option<PointF> {
    let p = s1p1;
    let r = s1p2 - s1p1;
    let q = s2p1;
    let s = s2p2 - s2p1;
    let rxs = cross_product(r, s);
    if rxs > -SKETCHER_EPSILON && rxs < SKETCHER_EPSILON {
        return None;
    }
    let qminusp = q - p;
    let t = cross_product(qminusp, s) / rxs;
    if !(0.0..=1.0).contains(&t) {
        return None;
    }
    let u = cross_product(qminusp, r) / rxs;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    Some(p + t * r)
}

/// `signedAngle`: angle between p1p2 and p2p3 in degrees
pub fn signed_angle(p1: PointF, p2: PointF, p3: PointF) -> f32 {
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    let a = (v1.x * v2.y - v1.y * v2.x).atan2(v1.x * v2.x + v1.y * v2.y);
    (f64::from(a * 180.0) / M_PI) as f32
}

/// `unsignedAngle`
pub fn unsigned_angle(p1: PointF, p2: PointF, p3: PointF) -> f32 {
    let (v1x, v1y) = (p1.x - p2.x, p1.y - p2.y);
    let (v2x, v2y) = (p3.x - p2.x, p3.y - p2.y);
    let mut d = (v1x * v1x + v1y * v1y).sqrt() * (v2x * v2x + v2y * v2y).sqrt();
    if d < SKETCHER_EPSILON {
        d = SKETCHER_EPSILON;
    }
    let mut cosine = (v1x * v2x + v1y * v2y) / d;
    if cosine < -1.0 {
        cosine = -1.0;
    } else if cosine > 1.0 {
        cosine = 1.0;
    }
    (f64::from(cosine.acos() * 180.0) / M_PI) as f32
}

pub fn points_coincide(p1: PointF, p2: PointF) -> bool {
    (p1 - p2).square_length() < SKETCHER_EPSILON * SKETCHER_EPSILON
}

/// `sameSide`: p1 and p2 are in the same half plane of the line through the segment
pub fn same_side(p1: PointF, p2: PointF, line_p1: PointF, line_p2: PointF) -> bool {
    let x = line_p2.x - line_p1.x;
    let y = line_p2.y - line_p1.y;
    if x.abs() > y.abs() {
        let m = y / x;
        let d1 = p1.y - line_p1.y - m * (p1.x - line_p1.x);
        let d2 = p2.y - line_p1.y - m * (p2.x - line_p1.x);
        d2 * d1 > 0.0
    } else {
        let m = x / y;
        let d1 = p1.x - line_p1.x - m * (p1.y - line_p1.y);
        let d2 = p2.x - line_p1.x - m * (p2.y - line_p1.y);
        d2 * d1 > 0.0
    }
}

pub fn project_point_on_line(p: PointF, sp1: PointF, sp2: PointF) -> PointF {
    let l1 = p - sp1;
    let l3 = sp2 - sp1;
    let mut segmentl2 = l3.square_length();
    if segmentl2 < SKETCHER_EPSILON {
        segmentl2 = SKETCHER_EPSILON;
    }
    let t = dot_product(l1, l3) / segmentl2;
    sp1 + t * l3
}

/// `squaredDistancePointSegment`: the squared distance and the clamped position `t` on the segment
pub fn squared_distance_point_segment(p: PointF, sp1: PointF, sp2: PointF) -> (f32, f32) {
    let l1 = p - sp1;
    let l2 = sp2 - p;
    let l3 = sp2 - sp1;
    let mut segmentl2 = l3.x * l3.x + l3.y * l3.y;
    if segmentl2 < SKETCHER_EPSILON {
        segmentl2 = SKETCHER_EPSILON;
    }
    let t = (l1.x * l3.x + l1.y * l3.y) / segmentl2;
    let return_t = if t < 0.0 {
        0.0
    } else if t > 1.0 {
        1.0
    } else {
        t
    };
    let mut squaredistance = if t < 0.0 {
        l1.x * l1.x + l1.y * l1.y
    } else if t > 1.0 {
        l2.x * l2.x + l2.y * l2.y
    } else {
        let proj = sp1 + t * l3;
        let l5 = p - proj;
        l5.x * l5.x + l5.y * l5.y
    };
    if squaredistance < SKETCHER_EPSILON {
        squaredistance = SKETCHER_EPSILON;
    }
    (squaredistance, return_t)
}

pub fn squared_distance(p1: PointF, p2: PointF) -> f32 {
    (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)
}

/// `tridiagonalSolve`
pub fn tridiagonal_solve(a: &[f32], b: &[f32], c: &[f32], rhs: &[f32]) -> Vec<f32> {
    let n = rhs.len();
    let mut u = vec![0.0f32; n];
    let mut gam = vec![0.0f32; n];
    let mut bet = b[0];
    u[0] = rhs[0] / bet;
    for j in 1..n {
        gam[j] = c[j - 1] / bet;
        bet = b[j] - a[j] * gam[j];
        u[j] = (rhs[j] - a[j] * u[j - 1]) / bet;
    }
    for j in 1..n {
        u[n - j - 1] -= gam[n - j] * u[n - j];
    }
    u
}

/// `cyclicSolve`
pub fn cyclic_solve(
    a: &[f32],
    b: &[f32],
    c: &[f32],
    alpha: f32,
    beta: f32,
    rhs: &[f32],
) -> Vec<f32> {
    let n = b.len();
    let gamma = -b[0];
    let mut bb = vec![0.0f32; n];
    bb[0] = b[0] - gamma;
    bb[n - 1] = b[n - 1] - alpha * beta / gamma;
    bb[1..n - 1].copy_from_slice(&b[1..n - 1]);
    let mut x = tridiagonal_solve(a, &bb, c, rhs);
    let mut u = vec![0.0f32; n];
    u[0] = gamma;
    u[n - 1] = alpha;
    let z = tridiagonal_solve(a, &bb, c, &u);
    let fact = f64::from((x[0] + beta * x[n - 1] / gamma) / (1.0 + z[0] + beta * z[n - 1] / gamma));
    for i in 0..n {
        x[i] -= (fact * f64::from(z[i])) as f32;
    }
    x
}

/// `pointOnCubicBezier` (de Casteljau)
pub fn point_on_cubic_bezier(p1: PointF, cp1: PointF, cp2: PointF, p2: PointF, t: f32) -> PointF {
    let v1 = (1.0 - t) * p1 + t * cp1;
    let v2 = (1.0 - t) * cp1 + t * cp2;
    let v3 = (1.0 - t) * cp2 + t * p2;
    let v4 = (1.0 - t) * v1 + t * v2;
    let v5 = (1.0 - t) * v2 + t * v3;
    (1.0 - t) * v4 + t * v5
}

/// `ClosedBezierControlPoints`: the two control point lists of a closed curve through the knots (empty for <= 2 knots)
pub fn closed_bezier_control_points(knots: &[PointF]) -> (Vec<PointF>, Vec<PointF>) {
    let n = knots.len();
    if n <= 2 {
        return (Vec::new(), Vec::new());
    }
    let a = vec![1.0f32; n];
    let b = vec![4.0f32; n];
    let c = vec![1.0f32; n];
    let mut rhs = vec![0.0f32; n];
    for i in 0..n {
        let j = if i + 1 > n - 1 { 0 } else { i + 1 };
        rhs[i] = 4.0 * knots[i].x + 2.0 * knots[j].x;
    }
    let x = cyclic_solve(&a, &b, &c, 1.0, 1.0, &rhs);
    for i in 0..n {
        let j = if i + 1 > n - 1 { 0 } else { i + 1 };
        rhs[i] = 4.0 * knots[i].y + 2.0 * knots[j].y;
    }
    let y = cyclic_solve(&a, &b, &c, 1.0, 1.0, &rhs);
    let first = (0..n).map(|i| PointF::new(x[i], y[i])).collect();
    let second = (0..n)
        .map(|i| PointF::new(2.0 * knots[i].x - x[i], 2.0 * knots[i].y - y[i]))
        .collect();
    (first, second)
}

/// `mirrorPoint`: the mirror image of `point` at the line through the segment
pub fn mirror_point(point: PointF, segment_point1: PointF, segment_point2: PointF) -> PointF {
    let segment_v = segment_point2 - segment_point1;
    let v2 = point - segment_point1;
    let parallel = v2.parallel_component(segment_v);
    let normal = v2 - parallel;
    segment_point1 + parallel - normal
}

/// `cannonBallDistance`: how far a ball of radius `ball_r` shot from `origin` travels before it hits a ball of
/// radius `target_r` at `target` (direction normalized)
pub fn cannon_ball_distance(
    origin: [f32; 3],
    direction: [f32; 3],
    target: [f32; 3],
    ball_r: f32,
    target_r: f32,
    cut_off: f32,
) -> f32 {
    let mut tdx = target[0] - origin[0];
    let mut tdy = target[1] - origin[1];
    let mut tdz = target[2] - origin[2];
    let rr = ball_r + target_r;
    let d2 = (tdx * tdx) + (tdy * tdy) + (tdz * tdz);
    if d2 > (cut_off + rr) * (cut_off + rr) {
        return cut_off;
    }
    if d2 < rr * rr {
        return 0.0;
    }
    let d = d2.sqrt();
    if d > SKETCHER_EPSILON {
        tdx /= d;
        tdy /= d;
        tdz /= d;
    }
    let cos = tdx * direction[0] + tdy * direction[1] + tdz * direction[2];
    if cos < 0.0 {
        return cut_off;
    }
    let sin = (1.0 - (cos * cos)).sqrt();
    let f = d * sin;
    if f > rr {
        return cut_off;
    }
    let result = (d2 - (f * f)).sqrt() - ((rr * rr) - (f * f)).sqrt();
    if result > cut_off { cut_off } else { result }
}

pub fn length_3d(x: f32, y: f32, z: f32) -> f32 {
    let m = x * x + y * y + z * z;
    if m > SKETCHER_EPSILON { m.sqrt() } else { m }
}

pub fn dot_product_3d(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> f32 {
    x1 * x2 + y1 * y2 + z1 * z2
}

pub fn cross_product_3d(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> (f32, f32, f32) {
    (y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2)
}

pub fn distance_3d(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> f32 {
    length_3d(x2 - x1, y2 - y1, z2 - z1)
}

/// `angle3D`: the angle at the second point in degrees
#[allow(clippy::too_many_arguments)]
pub fn angle_3d(
    x1: f32,
    y1: f32,
    z1: f32,
    x2: f32,
    y2: f32,
    z2: f32,
    x3: f32,
    y3: f32,
    z3: f32,
) -> f32 {
    let (xa, ya, za) = (x1 - x2, y1 - y2, z1 - z2);
    let (xb, yb, zb) = (x3 - x2, y3 - y2, z3 - z2);
    let l1 = length_3d(xa, ya, za);
    let l2 = length_3d(xb, yb, zb);
    let dp = dot_product_3d(xa, ya, za, xb, yb, zb);
    (f64::from((dp / (l1 * l2)).acos() * 180.0) / M_PI) as f32
}

/// `dihedral3D`
#[allow(clippy::too_many_arguments)]
pub fn dihedral_3d(
    x1: f32,
    y1: f32,
    z1: f32,
    x2: f32,
    y2: f32,
    z2: f32,
    x3: f32,
    y3: f32,
    z3: f32,
    x4: f32,
    y4: f32,
    z4: f32,
) -> f32 {
    let (xa, ya, za) = cross_product_3d(x1 - x2, y1 - y2, z1 - z2, x3 - x2, y3 - y2, z3 - z2);
    let (xb, yb, zb) = cross_product_3d(x2 - x3, y2 - y3, z2 - z3, x4 - x3, y4 - y3, z4 - z3);
    angle_3d(xa, ya, za, 0.0, 0.0, 0.0, xb, yb, zb)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f(s: &str) -> f32 {
        f32::from_bits(u32::from_str_radix(s, 16).unwrap())
    }

    fn p(v: &[f32], i: usize) -> PointF {
        PointF::new(v[2 * i], v[2 * i + 1])
    }

    /// every function must give the bits the C++ header gives (`tests/data/coordgen_maths.txt`, driver
    /// `tools/coordgen_drv/maths.cpp`)
    #[test]
    fn maths_matches_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_maths.txt");
        let mut bad = Vec::new();
        let mut checked = 0;
        for line in text.lines() {
            let (head, tail) = line.split_once(" ;").unwrap();
            let mut it = head.split(' ');
            let name = it.next().unwrap();
            let v: Vec<f32> = it.map(f).collect();
            let mut res: Vec<&str> = tail.split_whitespace().collect();
            // a trailing small decimal is the boolean result of the predicates
            let flag = matches!(name, "intersect" | "coincide" | "sameSide")
                .then(|| res.pop().unwrap().parse::<i32>().unwrap());
            let want: Vec<f32> = res.iter().map(|s| f(s)).collect();
            let mut got: Vec<f32> = Vec::new();
            let mut got_flag = None;
            match name {
                "round2" => got.push(round_to_two_decimal_digits(v[0])),
                "roundp" => got.push(round_to_precision(v[0], v[1] as i32)),
                "normalize" => {
                    let a = p(&v, 0);
                    let mut t = a;
                    t.normalize();
                    got = vec![a.length(), t.x, t.y];
                }
                "rotate" => {
                    let mut t = p(&v, 0);
                    t.rotate(v[2], v[3]);
                    got = vec![t.x, t.y];
                }
                "parallel" => {
                    let r = p(&v, 0).parallel_component(p(&v, 1));
                    got = vec![r.x, r.y];
                }
                "pround" => {
                    let (mut t, mut t3) = (p(&v, 0), p(&v, 0));
                    t.round(2);
                    t3.round(3);
                    got = vec![t.x, t.y, t3.x, t3.y];
                }
                "intersect" => {
                    let r = intersection_of_segments(p(&v, 0), p(&v, 1), p(&v, 2), p(&v, 3));
                    got_flag = Some(i32::from(r.is_some()));
                    let r = r.unwrap_or_default();
                    got = vec![r.x, r.y];
                    // the C++ result stays at its default when there is no crossing
                }
                "signedAngle" => got.push(signed_angle(p(&v, 0), p(&v, 1), p(&v, 2))),
                "unsignedAngle" => got.push(unsigned_angle(p(&v, 0), p(&v, 1), p(&v, 2))),
                "coincide" => got_flag = Some(i32::from(points_coincide(p(&v, 0), p(&v, 1)))),
                "sameSide" => {
                    got_flag = Some(i32::from(same_side(p(&v, 0), p(&v, 1), p(&v, 2), p(&v, 3))))
                }
                "project" => {
                    let r = project_point_on_line(p(&v, 0), p(&v, 1), p(&v, 2));
                    got = vec![r.x, r.y];
                }
                "sqdistseg" => {
                    let (d, t) = squared_distance_point_segment(p(&v, 0), p(&v, 1), p(&v, 2));
                    got = vec![d, t];
                }
                "sqdist" => got.push(squared_distance(p(&v, 0), p(&v, 1))),
                "mirror" => {
                    let r = mirror_point(p(&v, 0), p(&v, 1), p(&v, 2));
                    got = vec![r.x, r.y];
                }
                "bezier" => {
                    let r = point_on_cubic_bezier(p(&v, 0), p(&v, 1), p(&v, 2), p(&v, 3), v[8]);
                    got = vec![r.x, r.y];
                }
                "length3D" => got.push(length_3d(v[0], v[1], v[2])),
                "dot3D" => got.push(dot_product_3d(v[0], v[1], v[2], v[3], v[4], v[5])),
                "cross3D" => {
                    let r = cross_product_3d(v[0], v[1], v[2], v[3], v[4], v[5]);
                    got = vec![r.0, r.1, r.2];
                }
                "dist3D" => got.push(distance_3d(v[0], v[1], v[2], v[3], v[4], v[5])),
                "angle3D" => got.push(angle_3d(
                    v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                )),
                "dihedral3D" => got.push(dihedral_3d(
                    v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11],
                )),
                "cannon" => got.push(cannon_ball_distance(
                    [v[0], v[1], v[2]],
                    [v[3], v[4], v[5]],
                    [v[6], v[7], v[8]],
                    v[9],
                    v[10],
                    4.0,
                )),
                "closedbezier" => {
                    let knots: Vec<PointF> = (0..v.len() / 2).map(|i| p(&v, i)).collect();
                    let (a, b) = closed_bezier_control_points(&knots);
                    got = a.iter().chain(&b).flat_map(|q| [q.x, q.y]).collect();
                }
                "tridiag" | "cyclic" => {
                    let n = v.len() / 4;
                    let col = |k: usize| -> Vec<f32> { (0..n).map(|i| v[4 * i + k]).collect() };
                    let (a, b, c, rhs) = (col(0), col(1), col(2), col(3));
                    got = if name == "tridiag" {
                        tridiagonal_solve(&a, &b, &c, &rhs)
                    } else {
                        cyclic_solve(&a, &b, &c, 1.0, 1.0, &rhs)
                    };
                }
                _ => continue,
            }
            checked += 1;
            // acos/atan2 come from the platform libm, which may differ in the last bit between systems
            let libm = matches!(
                name,
                "signedAngle" | "unsignedAngle" | "angle3D" | "dihedral3D"
            );
            let same = got.len() == want.len()
                && got.iter().zip(&want).all(|(a, b)| {
                    a.to_bits() == b.to_bits()
                        || (a.is_nan() && b.is_nan())
                        || (libm && (a - b).abs() <= 2e-6 * b.abs().max(1.0))
                })
                && got_flag == flag;
            // `intersect` leaves the point at its default without a crossing
            let same = same || (name == "intersect" && flag == Some(0) && got_flag == flag);
            if !same {
                bad.push(format!("{line}\n got {got:?} {got_flag:?}"));
            }
        }
        assert!(checked > 7900, "{checked}");
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
