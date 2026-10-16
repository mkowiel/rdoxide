//! Geometry helpers of the 2D molecule drawing code
//! (`MolDraw2DDetails.cpp` of RDKit 2026.09.1).
//!
//! The comic-mode `handdrawnLine` (driven by `std::rand`) and the
//! substance-group data labels are ported together with the code that uses
//! them.
//!
//! Copyright (C) 2015-2020 Greg Landrum. Licensed under the BSD 3-Clause
//! license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::string_rect::StringRect;
use crate::Result;
use crate::geometry::Point2D;

fn scaled(p: Point2D, f: f64) -> Point2D {
    Point2D::new(p.x * f, p.y * f)
}

/// Format used for coordinates in SVG output (`formatDouble`).
pub fn format_double(val: f64) -> String {
    format!("{val:.1}")
}

/// Points on an elliptical arc (from `sping/pid.py`). The calculation is done
/// in `f32` like upstream so the output matches.
pub fn arc_points(cds1: Point2D, cds2: Point2D, start_ang: f32, extent: f32) -> Vec<Point2D> {
    let x_scale = ((cds2.x - cds1.x) / 2.0).abs() as f32;
    let y_scale = ((cds2.y - cds1.y) / 2.0).abs() as f32;
    let x = (cds1.x.min(cds2.x) + x_scale as f64) as f32;
    let y = (cds1.y.min(cds2.y) + y_scale as f64) as f32;
    let steps = ((extent * 2.0) as i32).max(5);
    let step = (PI * extent as f64 / (180 * steps) as f64) as f32;
    let mut angle = (PI * start_ang as f64 / 180.0) as f32;
    let mut res = Vec::with_capacity(steps as usize + 1);
    for _ in 0..=steps {
        res.push(Point2D::new(
            (x + x_scale * angle.cos()) as f64,
            (y - y_scale * angle.sin()) as f64,
        ));
        angle += step;
    }
    res
}

/// Approximate segment intersection test, good enough for drawing.
fn line_segments_intersect(s1: Point2D, s2: Point2D, s3: Point2D, s4: Point2D) -> bool {
    let mut d1x = s1.x - s2.x;
    let d1y = s1.y - s2.y;
    let mut d2x = s3.x - s4.x;
    let d2y = s3.y - s4.y;
    if d1x.abs() < 1e-4 {
        d1x = 1e-4;
    }
    if d2x.abs() < 1e-4 {
        d2x = 1e-4;
    }
    let m1 = d1y / d1x;
    let m2 = d2y / d2x;
    if m1 == m2 || m1 == -m2 {
        return false;
    }
    let b1 = (s1.x * s2.y - s2.x * s1.y) / d1x;
    let b2 = (s3.x * s4.y - s4.x * s3.y) / d2x;
    let ix = (b2 - b1) / (m1 - m2);
    ((ix < s1.x) ^ (ix < s2.x)) && ((ix < s3.x) ^ (ix < s4.x))
}

/// The four points of a bracket between `p1` and `p2`, opening away from
/// `ref_pt` (or from the first bond segment the bracket crosses).
pub fn get_bracket_points(
    p1: Point2D,
    p2: Point2D,
    ref_pt: Point2D,
    bond_segments: &[(Point2D, Point2D)],
    bracket_frac: f64,
) -> Vec<Point2D> {
    let v = p2 - p1;
    let mut bracket_dir = scaled(Point2D::new(v.y, -v.x), bracket_frac);
    let mut ref_vect = p2 - ref_pt;
    for seg in bond_segments {
        if line_segments_intersect(p1, p2, seg.0, seg.1) {
            ref_vect = p2 - seg.0;
        }
    }
    if bracket_dir.dot_product(&ref_vect) > 0.0 {
        bracket_dir = scaled(bracket_dir, -1.0);
    }
    vec![p1 + bracket_dir, p1, p2, p2 + bracket_dir]
}

/// Intersection of the segments `l1s`-`l1f` and `l2s`-`l2f`, if any.
pub fn do_lines_intersect(
    l1s: Point2D,
    l1f: Point2D,
    l2s: Point2D,
    l2f: Point2D,
) -> Option<Point2D> {
    let s1_x = l1f.x - l1s.x;
    let s1_y = l1f.y - l1s.y;
    let s2_x = l2f.x - l2s.x;
    let s2_y = l2f.y - l2s.y;
    let d = -s2_x * s1_y + s1_x * s2_y;
    if d == 0.0 {
        return None;
    }
    let s = (-s1_y * (l1s.x - l2s.x) + s1_x * (l1s.y - l2s.y)) / d;
    let t = (s2_x * (l1s.y - l2s.y) - s2_y * (l1s.x - l2s.x)) / d;
    if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
        Some(Point2D::new(l1s.x + t * s1_x, l1s.y + t * s1_y))
    } else {
        None
    }
}

pub fn is_point_in_triangle(pt: Point2D, t1: Point2D, t2: Point2D, t3: Point2D) -> bool {
    let d = (t2.y - t3.y) * (t1.x - t3.x) + (t3.x - t2.x) * (t1.y - t3.y);
    let a = ((t2.y - t3.y) * (pt.x - t3.x) + (t3.x - t2.x) * (pt.y - t3.y)) / d;
    let b = ((t3.y - t1.y) * (pt.x - t3.x) + (t1.x - t3.x) * (pt.y - t3.y)) / d;
    let c = 1.0 - a - b;
    (0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b) && (0.0..=1.0).contains(&c)
}

/// The sides of a rectangle as `(tl,tr), (tr,br), (br,bl), (bl,tl)`.
fn rect_sides(rect: &StringRect, padding: f64) -> [(Point2D, Point2D); 4] {
    let [tl, tr, br, bl] = rect.calc_corners(padding);
    [(tl, tr), (tr, br), (br, bl), (bl, tl)]
}

/// Does the line `end1`-`end2` cross a side of the (padded) rectangle?
pub fn does_line_intersect(rect: &StringRect, end1: Point2D, end2: Point2D, padding: f64) -> bool {
    rect_sides(rect, padding)
        .iter()
        .any(|&(a, b)| do_lines_intersect(end2, end1, a, b).is_some())
}

/// Does the triangle overlap the (padded) rectangle in any way?
pub fn does_triangle_intersect(
    rect: &StringRect,
    pt1: Point2D,
    pt2: Point2D,
    pt3: Point2D,
    padding: f64,
) -> bool {
    if [pt1, pt2, pt3]
        .iter()
        .any(|&p| rect.is_point_inside(p, padding))
    {
        return true;
    }
    let corners = rect.calc_corners(padding);
    if corners
        .iter()
        .any(|&c| is_point_in_triangle(c, pt1, pt2, pt3))
    {
        return true;
    }
    let tri = [(pt1, pt2), (pt2, pt3), (pt3, pt1)];
    rect_sides(rect, padding).iter().any(|&(a, b)| {
        tri.iter()
            .any(|&(c, d)| do_lines_intersect(a, b, c, d).is_some())
    })
}

/// True if either end of the line is inside the (padded) ellipse.
pub fn does_line_intersect_ellipse(
    centre: Point2D,
    xradius: f64,
    yradius: f64,
    padding: f64,
    end1: Point2D,
    end2: Point2D,
) -> bool {
    let xr2 = (xradius + padding) * (xradius + padding);
    let yr2 = (yradius + padding) * (yradius + padding);
    let inside = |p: Point2D| {
        (p.x - centre.x) * (p.x - centre.x) / xr2 + (p.y - centre.y) * (p.y - centre.y) / yr2 <= 1.0
    };
    inside(end1) || inside(end2)
}

/// True if either end of the line is inside the arc (angles in degrees; the
/// ends are given relative to the centre by the caller of upstream's
/// `doesLineIntersectArc`, which subtracts `centre` itself).
pub fn does_line_intersect_arc(
    centre: Point2D,
    xradius: f64,
    yradius: f64,
    start_ang: f64,
    stop_ang: f64,
    padding: f64,
    end1: Point2D,
    end2: Point2D,
) -> bool {
    let xr = xradius + padding;
    let yr = yradius + padding;
    let (xr2, yr2) = (xr * xr, yr * yr);
    let point_in_arc = |p: Point2D| -> bool {
        if p.x * p.x / xr2 + p.y * p.y / yr2 > 1.0 {
            return false;
        }
        let mut th = p.x.atan2(p.y) * 180.0 / PI;
        if th < 0.0 {
            th += 360.0;
        }
        if start_ang < stop_ang {
            th >= start_ang && th <= stop_ang
        } else {
            (th >= start_ang && th <= 360.0) || (th >= 0.0 && th <= stop_ang)
        }
    };
    point_in_arc(end1 - centre) || point_in_arc(end2 - centre)
}

/// Segments `(start, control1, control2, end)` of a wavy line.
pub fn get_wavy_line_segments(
    p1: Point2D,
    p2: Point2D,
    n_segments: u32,
    vert_offset: f64,
) -> Result<Vec<(Point2D, Point2D, Point2D, Point2D)>> {
    if n_segments <= 1 {
        return Err(crate::RDKitError::precondition("too few segments"));
    }
    // an even number of segments is assumed
    let n_segments = n_segments + n_segments % 2;
    let delta = p2 - p1;
    let mut perp = Point2D::new(delta.y, -delta.x);
    perp.normalize()?;
    let perp = scaled(perp, vert_offset);
    let delta = Point2D::new(delta.x / n_segments as f64, delta.y / n_segments as f64);
    Ok((0..n_segments)
        .map(|i| {
            let start = p1 + scaled(delta, i as f64);
            let seg = start + delta;
            let sign = if i % 2 == 1 { -1.0 } else { 1.0 };
            (
                start,
                start + scaled(perp, sign),
                seg + scaled(perp, sign),
                seg,
            )
        })
        .collect())
}

/// Shortens the arrow so the mitred head fits and returns the two barbs.
/// `arrow_end` is updated in place.
pub fn calc_arrow_head(
    arrow_end: &mut Point2D,
    arrow_begin: Point2D,
    mut frac: f64,
    line_width: f64,
    angle: f64,
) -> Option<(Point2D, Point2D)> {
    if angle < 1.0e-6 {
        return None;
    }
    let adjuster = 0.5 * line_width / angle.sin();
    let len = (arrow_begin - *arrow_end).length();
    if len > 1.0e-6 {
        let adj_len = len - adjuster;
        arrow_end.x = arrow_begin.x + (arrow_end.x - arrow_begin.x) * adj_len / len;
        arrow_end.y = arrow_begin.y + (arrow_end.y - arrow_begin.y) * adj_len / len;
    }
    let delta = arrow_begin - *arrow_end;
    let (cos_a, sin_a) = (angle.cos(), angle.sin());
    frac /= cos_a;
    let arrow1 = Point2D::new(
        arrow_end.x + frac * (delta.x * cos_a + delta.y * sin_a),
        arrow_end.y + frac * (delta.y * cos_a - delta.x * sin_a),
    );
    let arrow2 = Point2D::new(
        arrow_end.x + frac * (delta.x * cos_a - delta.y * sin_a),
        arrow_end.y + frac * (delta.y * cos_a + delta.x * sin_a),
    );
    Some((arrow1, arrow2))
}

/// Moves `p2` onto the ellipse where the line from `p1` crosses it.
pub fn adjust_line_end_for_ellipse(
    centre: Point2D,
    xradius: f64,
    yradius: f64,
    p1: Point2D,
    p2: &mut Point2D,
) {
    let p1 = p1 - centre;
    let q2 = *p2 - centre;
    let a2 = xradius * xradius;
    let b2 = yradius * yradius;
    let a = (q2.x - p1.x) * (q2.x - p1.x) / a2 + (q2.y - p1.y) * (q2.y - p1.y) / b2;
    let b = 2.0 * p1.x * (q2.x - p1.x) / a2 + 2.0 * p1.y * (q2.y - p1.y) / b2;
    let c = p1.x * p1.x / a2 + p1.y * p1.y / b2 - 1.0;
    let t_to_point = |t: f64| {
        Point2D::new(
            p1.x + (q2.x - p1.x) * t + centre.x,
            p1.y + (q2.y - p1.y) * t + centre.y,
        )
    };
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        *p2 = q2 + centre;
    } else if disc.abs() < 1.0e-6 {
        *p2 = t_to_point(-b / (2.0 * a));
    } else {
        let disc_rt = disc.sqrt();
        let t1 = (-b + disc_rt) / (2.0 * a);
        let t2 = (-b - disc_rt) / (2.0 * a);
        let t1_ok = (0.0..=1.0).contains(&t1);
        let t2_ok = (0.0..=1.0).contains(&t2);
        let t = match (t1_ok, t2_ok) {
            (true, false) => t1,
            (false, true) => t2,
            (true, true) => t1.min(t2),
            // both intersections are outside the segment: leave p2 alone
            (false, false) => {
                *p2 = q2 + centre;
                return;
            }
        };
        *p2 = t_to_point(t);
    }
}
