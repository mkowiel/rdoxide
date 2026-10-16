//! Contouring of data on a grid (`Numerics/Conrec.h` of RDKit 2026.09.1),
//! after Paul Bourke's CONREC.
//!
//! Copyright (C) 2020 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::HashMap;

use crate::geometry::Point2D;
use crate::{RDKitError, Result};

/// A segment of a contour at an iso-value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConrecSegment {
    pub p1: Point2D,
    pub p2: Point2D,
    pub iso_val: f64,
}

const IM: [usize; 4] = [0, 1, 1, 0];
const JM: [usize; 4] = [0, 0, 1, 1];
const CASTAB: [[[u8; 3]; 3]; 3] = [
    [[0, 0, 8], [0, 2, 5], [7, 6, 9]],
    [[0, 3, 4], [1, 3, 1], [4, 3, 0]],
    [[9, 6, 7], [5, 2, 0], [8, 0, 0]],
];

/// The contour segments of the data `d`, indexed `d[i * ny + j]` where `ny` is
/// `jub - jlb + 1`, over the index bounds `ilb..=iub` and `jlb..=jub`. `x`
/// and `y` are the column and row coordinates and `z` the contour levels in
/// increasing order.
pub fn contour(
    d: &[f64],
    ilb: usize,
    iub: usize,
    jlb: usize,
    jub: usize,
    x: &[f64],
    y: &[f64],
    nc: usize,
    z: &[f64],
) -> Result<Vec<ConrecSegment>> {
    if d.is_empty() || x.is_empty() || y.is_empty() || z.is_empty() {
        return Err(RDKitError::precondition("no data"));
    }
    if nc == 0 {
        return Err(RDKitError::precondition("no contours"));
    }
    if iub <= ilb || jub <= jlb {
        return Err(RDKitError::precondition("bad bounds"));
    }
    let ny = jub - jlb + 1;
    let mut res = Vec::new();
    let mut h = [0.0f64; 5];
    let mut sh = [0i32; 5];
    let mut xh = [0.0f64; 5];
    let mut yh = [0.0f64; 5];

    for j in (jlb..jub).rev() {
        for i in ilb..iub {
            let temp1 = d[i * ny + j].min(d[i * ny + j + 1]);
            let temp2 = d[(i + 1) * ny + j].min(d[(i + 1) * ny + j + 1]);
            let dmin = temp1.min(temp2);
            let temp1 = d[i * ny + j].max(d[i * ny + j + 1]);
            let temp2 = d[(i + 1) * ny + j].max(d[(i + 1) * ny + j + 1]);
            let dmax = temp1.max(temp2);
            if dmax < z[0] || dmin > z[nc - 1] {
                continue;
            }
            for &zk in &z[..nc] {
                if zk < dmin || zk > dmax {
                    continue;
                }
                for m in (0..=4).rev() {
                    if m > 0 {
                        h[m] = d[(i + IM[m - 1]) * ny + j + JM[m - 1]] - zk;
                        xh[m] = x[i + IM[m - 1]];
                        yh[m] = y[j + JM[m - 1]];
                    } else {
                        h[0] = 0.25 * (h[1] + h[2] + h[3] + h[4]);
                        xh[0] = 0.50 * (x[i] + x[i + 1]);
                        yh[0] = 0.50 * (y[j] + y[j + 1]);
                    }
                    sh[m] = if h[m] > 0.0 {
                        1
                    } else if h[m] < 0.0 {
                        -1
                    } else {
                        0
                    };
                }
                let xsect =
                    |p1: usize, p2: usize| (h[p2] * xh[p1] - h[p1] * xh[p2]) / (h[p2] - h[p1]);
                let ysect =
                    |p1: usize, p2: usize| (h[p2] * yh[p1] - h[p1] * yh[p2]) / (h[p2] - h[p1]);
                // Scan each triangle in the box: the centre is vertex 0 and
                // the corners 1 to 4.
                for m in 1..=4usize {
                    let m1 = m;
                    let m2 = 0;
                    let m3 = if m != 4 { m + 1 } else { 1 };
                    let case_value =
                        CASTAB[(sh[m1] + 1) as usize][(sh[m2] + 1) as usize][(sh[m3] + 1) as usize];
                    let (x1, y1, x2, y2) = match case_value {
                        0 => continue,
                        // line between vertices 1 and 2
                        1 => (xh[m1], yh[m1], xh[m2], yh[m2]),
                        // between vertices 2 and 3
                        2 => (xh[m2], yh[m2], xh[m3], yh[m3]),
                        // between vertices 3 and 1
                        3 => (xh[m3], yh[m3], xh[m1], yh[m1]),
                        // vertex 1 and side 2-3
                        4 => (xh[m1], yh[m1], xsect(m2, m3), ysect(m2, m3)),
                        // vertex 2 and side 3-1
                        5 => (xh[m2], yh[m2], xsect(m3, m1), ysect(m3, m1)),
                        // vertex 3 and side 1-2
                        6 => (xh[m3], yh[m3], xsect(m1, m2), ysect(m1, m2)),
                        // sides 1-2 and 2-3
                        7 => (xsect(m1, m2), ysect(m1, m2), xsect(m2, m3), ysect(m2, m3)),
                        // sides 2-3 and 3-1
                        8 => (xsect(m2, m3), ysect(m2, m3), xsect(m3, m1), ysect(m3, m1)),
                        // sides 3-1 and 1-2
                        9 => (xsect(m3, m1), ysect(m3, m1), xsect(m1, m2), ysect(m1, m2)),
                        _ => continue,
                    };
                    res.push(ConrecSegment {
                        p1: Point2D::new(x1, y1),
                        p2: Point2D::new(x2, y2),
                        iso_val: zk,
                    });
                }
            }
        }
    }
    Ok(res)
}

type PointKey = (i32, i32, i64);

/// Joins segments that share end points (to the quantization given) into
/// continuous lines, each with its iso-value.
///
/// Upstream walks an `unordered_map` to find the starting points, so the order
/// of the lines depends on the C++ library; here they come in the order the
/// end points were first seen.
pub fn connect_line_segments(
    segments: &[ConrecSegment],
    coord_multiplier: f64,
    iso_val_multiplier: f64,
) -> Vec<(Vec<Point2D>, f64)> {
    let key = |pt: &Point2D, iso: f64| -> PointKey {
        (
            (coord_multiplier * pt.x).round() as i32,
            (coord_multiplier * pt.y).round() as i32,
            (iso_val_multiplier * iso).round() as i64,
        )
    };
    // the end point keys in first-seen order, and the segments at each
    let mut order: Vec<PointKey> = Vec::new();
    let mut end_points: HashMap<PointKey, Vec<usize>> = HashMap::new();
    for (i, seg) in segments.iter().enumerate() {
        for pt in [&seg.p1, &seg.p2] {
            let k = key(pt, seg.iso_val);
            end_points
                .entry(k)
                .or_insert_with(|| {
                    order.push(k);
                    Vec::new()
                })
                .push(i);
        }
    }
    let mut done = vec![false; segments.len()];
    let mut res = Vec::new();
    // a candidate end point has a segment that hasn't been used already
    let candidate = |end_points: &HashMap<PointKey, Vec<usize>>, done: &[bool], k: &PointKey| {
        end_points[k].first().is_some_and(|&s| !done[s])
    };
    let mut start = 0;
    loop {
        while start < order.len() && !candidate(&end_points, &done, &order[start]) {
            start += 1;
        }
        let Some(&first_key) = order.get(start) else {
            break;
        };
        let mut seg_id = end_points[&first_key][0];
        let mut curr_key = first_key;
        let curr_val = segments[seg_id].iso_val;
        let mut contour: Vec<Point2D> = Vec::new();
        loop {
            done[seg_id] = true;
            // move onto the next segment
            let seg = segments[seg_id];
            let k1 = key(&seg.p1, seg.iso_val);
            let k2 = key(&seg.p2, seg.iso_val);
            let mut end_key = k2;
            if k1 == curr_key {
                contour.push(seg.p1);
            } else if k2 == curr_key {
                contour.push(seg.p2);
                end_key = k1;
            }
            // remove this segment from the two lists
            if let Some(l) = end_points.get_mut(&curr_key) {
                if let Some(p) = l.iter().position(|&s| s == seg_id) {
                    l.remove(p);
                }
            }
            let segs = end_points.entry(end_key).or_default();
            if let Some(p) = segs.iter().position(|&s| s == seg_id) {
                segs.remove(p);
            }
            if segs.is_empty() {
                // we're at the end, push on the last point
                if k1 == curr_key {
                    contour.push(seg.p2);
                } else if k2 == curr_key {
                    contour.push(seg.p1);
                }
                break;
            }
            seg_id = segs[0];
            curr_key = end_key;
        }
        res.push((contour, curr_val));
    }
    res
}
