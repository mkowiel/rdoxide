//! The shapes the 2D drawing code builds before rendering them
//! (`DrawShape.cpp` of RDKit 2026.09.1).
//!
//! Each shape keeps its points in drawing coordinates so the drawing can be
//! scaled and moved as a whole; [`DrawShape::draw`] renders it on any
//! [`ShapeCanvas`].
//!
//! Copyright (C) 2021-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::details::{
    calc_arrow_head, do_lines_intersect, does_line_intersect, does_line_intersect_arc,
    does_line_intersect_ellipse, does_triangle_intersect,
};
use super::helpers::{DashPattern, DrawColour};
use super::string_rect::StringRect;
use crate::geometry::Point2D;
use crate::{RDKitError, Result};

pub fn no_dash() -> DashPattern {
    Vec::new()
}
pub fn dots() -> DashPattern {
    vec![2.0, 6.0]
}
pub fn dashes() -> DashPattern {
    vec![6.0, 4.0]
}
pub fn short_dashes() -> DashPattern {
    vec![2.0, 2.0]
}

/// The drawing operations shapes need from a drawer (the part of `MolDraw2D`
/// that `DrawShape` uses). Points passed to the `draw_*` methods are already
/// in drawing coordinates.
pub trait ShapeCanvas {
    fn colour(&self) -> DrawColour;
    fn set_colour(&mut self, col: DrawColour);
    fn line_width(&self) -> f64;
    fn set_line_width(&mut self, width: f64);
    fn fill_polys(&self) -> bool;
    fn set_fill_polys(&mut self, fill: bool);
    fn scale_bond_width(&self) -> bool;
    fn set_scale_bond_width(&mut self, scale: bool);
    fn split_bonds(&self) -> bool;
    fn dash(&self) -> DashPattern;
    fn set_dash(&mut self, dash: DashPattern);
    /// Sets the atoms (or -1) written as active in the output.
    fn set_active_atm_idx(&mut self, atom1: i32, atom2: i32);
    fn set_active_bnd_idx(&mut self, bond: i32);
    #[allow(clippy::too_many_arguments)]
    fn draw_arrow(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        as_polygon: bool,
        frac: f64,
        angle: f64,
        col: DrawColour,
        raw_coords: bool,
    ) -> Result<()>;
    fn draw_ellipse(&mut self, cds1: Point2D, cds2: Point2D, raw_coords: bool) -> Result<()>;
    fn draw_line(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col1: DrawColour,
        col2: DrawColour,
        raw_coords: bool,
    ) -> Result<()>;
    fn draw_polygon(&mut self, cds: &[Point2D], raw_coords: bool) -> Result<()>;
    fn draw_triangle(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        cds3: Point2D,
        raw_coords: bool,
    ) -> Result<()>;
    #[allow(clippy::too_many_arguments)]
    fn draw_wavy_line(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col1: DrawColour,
        col2: DrawColour,
        n_segments: u32,
        vert_offset: f64,
        raw_coords: bool,
    ) -> Result<()>;
    #[allow(clippy::too_many_arguments)]
    fn draw_arc(
        &mut self,
        centre: Point2D,
        xradius: f64,
        yradius: f64,
        angle1: f64,
        angle2: f64,
        raw_coords: bool,
    ) -> Result<()>;
}

/// What distinguishes the kinds of shape.
#[derive(Clone, Debug)]
pub enum ShapeKind {
    /// Two points; `orig_pts` are the ends and arrowhead ends, unadjusted
    /// for the line width, used for extremes and clashes.
    Arrow {
        frac: f64,
        angle: f64,
        orig_pts: [Point2D; 4],
    },
    /// Centre and x/y radii.
    Ellipse,
    SimpleLine {
        dash_pattern: DashPattern,
    },
    PolyLine {
        dash_pattern: DashPattern,
    },
    SolidWedge {
        col2: DrawColour,
        split_bonds: bool,
        other_bond_vecs: Vec<Point2D>,
    },
    DashedWedge {
        col2: DrawColour,
        one_less_dash: bool,
        line_colours: Vec<DrawColour>,
        /// The initial points, kept to rebuild the lines after scaling.
        at1: Point2D,
        end1: Point2D,
        end2: Point2D,
    },
    WavyLine {
        col2: DrawColour,
        offset: f64,
    },
    /// Centre and x/y radii, angles in degrees.
    Arc {
        ang1: f64,
        ang2: f64,
    },
}

/// A shape to be drawn.
#[derive(Clone, Debug)]
pub struct DrawShape {
    pub points: Vec<Point2D>,
    pub line_width: f64,
    pub scale_line_width: bool,
    pub line_colour: DrawColour,
    pub fill: bool,
    pub atom1: i32,
    pub atom2: i32,
    pub bond: i32,
    pub kind: ShapeKind,
}

fn mul(p: Point2D, f: f64) -> Point2D {
    Point2D::new(p.x * f, p.y * f)
}

fn half_sum(a: Point2D, b: Point2D) -> Point2D {
    Point2D::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0)
}

fn check_points(points: &[Point2D], n: usize, msg: &str) -> Result<()> {
    if points.len() == n {
        Ok(())
    } else {
        Err(RDKitError::precondition(msg))
    }
}

impl DrawShape {
    fn base(
        points: Vec<Point2D>,
        line_width: f64,
        scale_line_width: bool,
        line_colour: DrawColour,
        fill: bool,
        atoms_bond: (i32, i32, i32),
        kind: ShapeKind,
    ) -> Self {
        DrawShape {
            points,
            line_width,
            scale_line_width,
            line_colour,
            fill,
            atom1: atoms_bond.0,
            atom2: atoms_bond.1,
            bond: atoms_bond.2,
            kind,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn arrow(
        points: Vec<Point2D>,
        line_width: f64,
        scale_line_width: bool,
        line_colour: DrawColour,
        fill: bool,
        atoms_bond: (i32, i32, i32),
        frac: f64,
        angle: f64,
    ) -> Result<Self> {
        check_points(&points, 2, "arrow bad points size")?;
        let mut ab = points[1];
        let (p1, p2) = calc_arrow_head(&mut ab, points[1], frac, 0.0, angle)
            .unwrap_or((Point2D::default(), Point2D::default()));
        let orig_pts = [points[0], points[1], p1, p2];
        Ok(Self::base(
            points,
            line_width,
            scale_line_width,
            line_colour,
            fill,
            atoms_bond,
            ShapeKind::Arrow {
                frac,
                angle,
                orig_pts,
            },
        ))
    }

    /// Default arrowhead of the upstream constructor: `frac` 0.2, 30 degrees.
    pub fn default_arrow_head() -> (f64, f64) {
        (0.2, PI / 6.0)
    }

    pub fn ellipse(
        points: Vec<Point2D>,
        line_width: f64,
        scale_line_width: bool,
        line_colour: DrawColour,
        fill: bool,
        atom1: i32,
    ) -> Result<Self> {
        check_points(&points, 2, "ellipse wrong points")?;
        Ok(Self::base(
            points,
            line_width,
            scale_line_width,
            line_colour,
            fill,
            (atom1, -1, -1),
            ShapeKind::Ellipse,
        ))
    }

    pub fn simple_line(
        points: Vec<Point2D>,
        line_width: f64,
        scale_line_width: bool,
        line_colour: DrawColour,
        atoms_bond: (i32, i32, i32),
        dash_pattern: DashPattern,
    ) -> Result<Self> {
        check_points(&points, 2, "simple line wrong number of points")?;
        Ok(Self::base(
            points,
            line_width,
            scale_line_width,
            line_colour,
            false,
            atoms_bond,
            ShapeKind::SimpleLine { dash_pattern },
        ))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn poly_line(
        points: Vec<Point2D>,
        line_width: f64,
        scale_line_width: bool,
        line_colour: DrawColour,
        fill: bool,
        atoms_bond: (i32, i32, i32),
        dash_pattern: DashPattern,
    ) -> Result<Self> {
        if points.len() <= 2 {
            return Err(RDKitError::precondition("polyline not enough points"));
        }
        Ok(Self::base(
            points,
            line_width,
            scale_line_width,
            line_colour,
            fill,
            atoms_bond,
            ShapeKind::PolyLine { dash_pattern },
        ))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn solid_wedge(
        points: Vec<Point2D>,
        col1: DrawColour,
        col2: DrawColour,
        split_bonds: bool,
        other_bond_vecs: Vec<Point2D>,
        line_width: f64,
        atoms_bond: (i32, i32, i32),
    ) -> Result<Self> {
        check_points(&points, 3, "solid wedge wrong points")?;
        let mut shape = Self::base(
            points,
            line_width / 2.0,
            false,
            col1,
            false,
            atoms_bond,
            ShapeKind::SolidWedge {
                col2,
                split_bonds,
                other_bond_vecs,
            },
        );
        shape.trim_other_bond_vecs()?;
        shape.order_other_bond_vecs()?;
        shape.build_triangles();
        Ok(shape)
    }

    pub fn dashed_wedge(
        points: Vec<Point2D>,
        col1: DrawColour,
        col2: DrawColour,
        one_less_dash: bool,
        line_width: f64,
        atoms_bond: (i32, i32, i32),
    ) -> Result<Self> {
        check_points(&points, 3, "dashed wedge wrong points")?;
        let (at1, end1, end2) = (points[0], points[1], points[2]);
        let mut shape = Self::base(
            points,
            line_width,
            false,
            col1,
            false,
            atoms_bond,
            ShapeKind::DashedWedge {
                col2,
                one_less_dash,
                line_colours: Vec::new(),
                at1,
                end1,
                end2,
            },
        );
        shape.build_lines();
        Ok(shape)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn wavy_line(
        points: Vec<Point2D>,
        line_width: f64,
        scale_line_width: bool,
        col1: DrawColour,
        col2: DrawColour,
        offset: f64,
        atoms_bond: (i32, i32, i32),
    ) -> Result<Self> {
        check_points(&points, 2, "wavy line wrong points")?;
        Ok(Self::base(
            points,
            line_width,
            scale_line_width,
            col1,
            false,
            atoms_bond,
            ShapeKind::WavyLine { col2, offset },
        ))
    }

    /// Angles are in degrees between 0 and 360, `ang2` > `ang1` unless the
    /// arc crosses 0. Points are the centre and the x/y radii.
    #[allow(clippy::too_many_arguments)]
    pub fn arc(
        points: Vec<Point2D>,
        ang1: f64,
        ang2: f64,
        line_width: f64,
        scale_line_width: bool,
        col1: DrawColour,
        fill: bool,
        atom1: i32,
    ) -> Result<Self> {
        check_points(&points, 2, "arc wrong points")?;
        for ang in [ang1, ang2] {
            if !(0.0..=360.0).contains(&ang) {
                return Err(RDKitError::Value(format!(
                    "Range Error\n\tang\n\tViolation occurred on {ang} not in [0, 360]"
                )));
            }
        }
        Ok(Self::base(
            points,
            line_width,
            scale_line_width,
            col1,
            fill,
            (atom1, -1, -1),
            ShapeKind::Arc { ang1, ang2 },
        ))
    }

    /// Renders the shape, restoring the canvas state afterwards.
    pub fn draw(&self, canvas: &mut dyn ShapeCanvas) -> Result<()> {
        let ocolour = canvas.colour();
        canvas.set_colour(self.line_colour);
        let olw = canvas.line_width();
        canvas.set_line_width(self.line_width);
        let ofill = canvas.fill_polys();
        canvas.set_fill_polys(self.fill);
        let line_scale = canvas.scale_bond_width();
        canvas.set_scale_bond_width(self.scale_line_width);
        canvas.set_active_atm_idx(self.atom1, self.atom2);
        canvas.set_active_bnd_idx(self.bond);

        let res = self.my_draw(canvas);

        canvas.set_active_atm_idx(-1, -1);
        canvas.set_active_bnd_idx(-1);
        canvas.set_colour(ocolour);
        canvas.set_line_width(olw);
        canvas.set_fill_polys(ofill);
        canvas.set_scale_bond_width(line_scale);
        res
    }

    fn my_draw(&self, d: &mut dyn ShapeCanvas) -> Result<()> {
        let p = &self.points;
        match &self.kind {
            ShapeKind::Arrow { frac, angle, .. } => {
                if d.split_bonds() {
                    d.set_active_atm_idx(self.atom2, -1);
                } else {
                    d.set_active_atm_idx(self.atom1, self.atom2);
                }
                d.set_active_bnd_idx(self.bond);
                d.draw_arrow(p[0], p[1], self.fill, *frac, *angle, self.line_colour, true)
            }
            ShapeKind::Ellipse => {
                if self.fill {
                    d.set_line_width(1.0);
                    d.set_scale_bond_width(false);
                }
                d.draw_ellipse(p[0] - p[1], p[0] + p[1], true)
            }
            ShapeKind::SimpleLine { dash_pattern } => {
                let od = d.dash();
                let mut dp = dash_pattern.clone();
                if dp == short_dashes() {
                    // roughly equivalent to the original checks on scale
                    let sq_len = (p[0] - p[1]).length_sq();
                    if sq_len < 55.0 {
                        dp[0] /= 4.0;
                        dp[1] /= 3.0;
                    } else if sq_len < 900.0 {
                        dp[0] /= 2.0;
                        dp[1] /= 1.5;
                    }
                }
                d.set_dash(dp);
                d.set_active_atm_idx(self.atom1, self.atom2);
                d.set_active_bnd_idx(self.bond);
                let res = d.draw_line(p[0], p[1], self.line_colour, self.line_colour, true);
                d.set_active_atm_idx(-1, -1);
                d.set_active_bnd_idx(-1);
                d.set_dash(od);
                res
            }
            ShapeKind::PolyLine { .. } => d.draw_polygon(p, true),
            ShapeKind::SolidWedge {
                col2, split_bonds, ..
            } => {
                let _ = split_bonds;
                d.set_fill_polys(true);
                if d.split_bonds() {
                    d.set_active_atm_idx(self.atom1, -1);
                } else {
                    d.set_active_atm_idx(self.atom1, self.atom2);
                }
                d.set_active_bnd_idx(self.bond);
                if p.len() == 3 || p.len() == 9 {
                    d.draw_triangle(p[0], p[1], p[2], true)?;
                }
                if p.len() == 6 {
                    // two triangles where the wedge goes onto a Y shape
                    d.draw_polygon(&[p[0], p[1], p[2], p[5]], true)?;
                } else if p.len() == 9 {
                    // a conventional two-colour wedge drawn as 3 triangles
                    if d.split_bonds() {
                        d.set_active_atm_idx(self.atom2, -1);
                    }
                    d.set_colour(*col2);
                    d.draw_polygon(&[p[4], p[5], p[6], p[7]], true)?;
                }
                Ok(())
            }
            ShapeKind::DashedWedge { line_colours, .. } => {
                d.set_fill_polys(false);
                d.set_active_atm_idx(self.atom1, self.atom2);
                d.set_active_bnd_idx(self.bond);
                let mut j = 0;
                let mut i = 0;
                while i < p.len() {
                    if d.split_bonds() {
                        if i < p.len() / 2 {
                            d.set_active_atm_idx(self.atom1, -1);
                        } else {
                            d.set_active_atm_idx(self.atom2, -1);
                        }
                    }
                    d.draw_line(p[i], p[i + 1], line_colours[j], line_colours[j], true)?;
                    i += 2;
                    j += 1;
                }
                Ok(())
            }
            ShapeKind::WavyLine { col2, offset } => {
                // negative offset because of inverted y coords
                let nsegs = ((p[0] - p[1]).length() / (offset * 2.0 / 3.0)).round() as i32;
                d.draw_wavy_line(
                    p[0],
                    p[1],
                    self.line_colour,
                    *col2,
                    nsegs as u32,
                    -offset / 2.0,
                    true,
                )
            }
            ShapeKind::Arc { ang1, ang2 } => {
                if self.fill {
                    d.set_line_width(1.0);
                    d.set_scale_bond_width(false);
                }
                let start_ang = if ang1 > ang2 { ang1 - 360.0 } else { *ang1 };
                d.draw_arc(p[0], p[1].x, p[1].y, start_ang, *ang2, true)
            }
        }
    }

    /// Grows `(xmin, xmax, ymin, ymax)` to include the shape.
    pub fn find_extremes(&self, ext: &mut (f64, f64, f64, f64)) {
        let (xmin, xmax, ymin, ymax) = ext;
        let mut grow = |p: &Point2D| {
            *xmin = xmin.min(p.x);
            *xmax = xmax.max(p.x);
            *ymin = ymin.min(p.y);
            *ymax = ymax.max(p.y);
        };
        match &self.kind {
            ShapeKind::Arrow { orig_pts, .. } => orig_pts.iter().for_each(grow),
            ShapeKind::Ellipse | ShapeKind::Arc { .. } => {
                // points[0] is the centre, points[1] the radii
                let (c, r) = (self.points[0], self.points[1]);
                *xmin = xmin.min(c.x - r.x);
                *xmax = xmax.max(c.x + r.x);
                *ymin = ymin.min(c.y - r.y);
                *ymax = ymax.max(c.y + r.y);
            }
            ShapeKind::DashedWedge {
                at1, end1, end2, ..
            } => {
                for p in [at1, end1, end2] {
                    grow(p);
                }
            }
            _ => self.points.iter().for_each(grow),
        }
    }

    pub fn scale(&mut self, f: Point2D) {
        for p in &mut self.points {
            p.x *= f.x;
            p.y *= f.y;
        }
        match &mut self.kind {
            ShapeKind::Arrow { orig_pts, .. } => {
                for p in orig_pts {
                    p.x *= f.x;
                    p.y *= f.y;
                }
            }
            ShapeKind::DashedWedge {
                at1, end1, end2, ..
            } => {
                for p in [at1, end1, end2] {
                    p.x *= f.x;
                    p.y *= f.y;
                }
                self.build_lines();
            }
            ShapeKind::WavyLine { offset, .. } => *offset *= f.x,
            _ => {}
        }
    }

    pub fn translate(&mut self, trans: Point2D) {
        match &mut self.kind {
            // only the centre moves
            ShapeKind::Ellipse | ShapeKind::Arc { .. } => self.points[0] += trans,
            kind => {
                for p in &mut self.points {
                    *p += trans;
                }
                match kind {
                    ShapeKind::Arrow { orig_pts, .. } => {
                        for p in orig_pts {
                            *p += trans;
                        }
                    }
                    ShapeKind::DashedWedge {
                        at1, end1, end2, ..
                    } => {
                        for p in [at1, end1, end2] {
                            *p += trans;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    pub fn does_rect_clash(&self, rect: &StringRect, padding: f64) -> bool {
        let padding = if self.scale_line_width {
            padding * self.line_width
        } else {
            padding
        };
        let p = &self.points;
        match &self.kind {
            ShapeKind::Arrow { orig_pts: o, .. } => {
                does_line_intersect(rect, o[0], o[1], padding)
                    || does_line_intersect(rect, o[1], o[2], padding)
                    || does_line_intersect(rect, o[1], o[3], padding)
                    || does_line_intersect(rect, o[2], o[3], padding)
            }
            ShapeKind::Ellipse => {
                let w = p[1].x.abs();
                let h = p[1].y.abs();
                let centre = p[0];
                let [tl, tr, br, bl] = rect.calc_corners(padding);
                [(tl, tr), (tr, br), (br, bl), (bl, tl)]
                    .iter()
                    .any(|&(a, b)| {
                        does_line_intersect_ellipse(centre, w / 2.0, h / 2.0, padding, a, b)
                    })
            }
            ShapeKind::SimpleLine { .. } => does_line_intersect(rect, p[0], p[1], padding),
            ShapeKind::PolyLine { .. } => {
                (1..p.len()).any(|i| does_line_intersect(rect, p[i - 1], p[i], padding))
                    || does_line_intersect(rect, p[0], p[p.len() - 1], padding)
            }
            ShapeKind::SolidWedge { .. } => {
                does_triangle_intersect(rect, p[0], p[1], p[2], padding)
                    || (p.len() >= 6
                        && (does_triangle_intersect(rect, p[3], p[4], p[5], padding)
                            || (p.len() >= 9
                                && does_triangle_intersect(rect, p[6], p[7], p[8], padding))))
            }
            ShapeKind::DashedWedge {
                at1, end1, end2, ..
            } => does_triangle_intersect(rect, *at1, *end1, *end2, padding),
            ShapeKind::WavyLine { offset, .. } => {
                does_line_intersect(rect, p[0], p[1], padding + offset)
            }
            ShapeKind::Arc { ang1, ang2 } => {
                let [tl, tr, br, bl] = rect.calc_corners(padding);
                [(tl, tr), (tr, br), (br, bl), (bl, tl)]
                    .iter()
                    .any(|&(a, b)| {
                        does_line_intersect_arc(p[0], p[1].x, p[1].y, *ang1, *ang2, padding, a, b)
                    })
            }
        }
    }

    // ---- solid wedge construction ----

    fn other_vecs(&self) -> &[Point2D] {
        match &self.kind {
            ShapeKind::SolidWedge {
                other_bond_vecs, ..
            } => other_bond_vecs,
            _ => &[],
        }
    }

    /// With more than 2 other bonds only the two with the widest angle
    /// between them are kept.
    fn trim_other_bond_vecs(&mut self) -> Result<()> {
        let ShapeKind::SolidWedge {
            other_bond_vecs, ..
        } = &mut self.kind
        else {
            return Ok(());
        };
        if other_bond_vecs.len() < 3 {
            return Ok(());
        }
        let (mut first, mut second) = (0, 1);
        let mut largest = -361.0;
        for i in 0..other_bond_vecs.len() - 1 {
            for j in i + 1..other_bond_vecs.len() {
                let ang = other_bond_vecs[i].angle_to(&other_bond_vecs[j])?;
                if ang > largest {
                    first = i;
                    second = j;
                    largest = ang;
                }
            }
        }
        *other_bond_vecs = vec![other_bond_vecs[first], other_bond_vecs[second]];
        Ok(())
    }

    /// With 2 other bonds, makes the first one be on the `points[1]` side.
    fn order_other_bond_vecs(&mut self) -> Result<()> {
        let (p1, p2) = (self.points[1], self.points[2]);
        let ShapeKind::SolidWedge {
            other_bond_vecs, ..
        } = &mut self.kind
        else {
            return Ok(());
        };
        if other_bond_vecs.len() != 2 {
            return Ok(());
        }
        let mid = half_sum(p1, p2);
        let midp1 = mid.direction_vector(&p1)?;
        let dot1 = midp1.dot_product(&other_bond_vecs[0]);
        let dot2 = midp1.dot_product(&other_bond_vecs[1]);
        if dot1 < dot2 {
            other_bond_vecs.swap(0, 1);
        }
        Ok(())
    }

    fn build_triangles(&mut self) {
        let ShapeKind::SolidWedge {
            col2, split_bonds, ..
        } = &self.kind
        else {
            return;
        };
        if self.line_colour != *col2 || *split_bonds {
            self.build_two_colour_triangles();
        } else {
            self.build_single_colour_triangles();
        }
    }

    /// Where the line from `point` through `end` meets the other bond line
    /// through `mid_end`; `end` itself if it doesn't.
    fn adjust_end(point: Point2D, end: Point2D, mid_end: Point2D, other: Point2D) -> Point2D {
        let side = mul(end - point, 2.0);
        do_lines_intersect(point, point + side, mid_end - other, mid_end + other).unwrap_or(end)
    }

    fn build_single_colour_triangles(&mut self) {
        let (point, end1, end2) = (self.points[0], self.points[1], self.points[2]);
        let mid_end = half_sum(end1, end2);
        let other = self.other_vecs().to_vec();
        self.points.clear();
        match other.len() {
            0 => self.points.extend([point, end1, end2]),
            1 => {
                let a1 = Self::adjust_end(point, end1, mid_end, other[0]);
                let a2 = Self::adjust_end(point, end2, mid_end, other[0]);
                self.points.extend([point, a1, a2]);
            }
            2 => {
                let a1 = Self::adjust_end(point, end1, mid_end, other[0]);
                self.points.extend([point, a1, mid_end]);
                let a2 = Self::adjust_end(point, end2, mid_end, other[1]);
                self.points.extend([point, mid_end, a2]);
            }
            _ => {}
        }
    }

    fn build_two_colour_triangles(&mut self) {
        let (point, end1, end2) = (self.points[0], self.points[1], self.points[2]);
        let mid_end = half_sum(end1, end2);
        let other = self.other_vecs().to_vec();
        self.points.clear();
        let mid1 = point + mul(end1 - point, 0.5);
        let mid2 = point + mul(end2 - point, 0.5);
        self.points.extend([point, mid1, mid2]);
        match other.len() {
            0 => self.points.extend([mid1, end2, end1, mid1, mid2, end2]),
            1 => {
                let a1 = Self::adjust_end(point, end1, mid_end, other[0]);
                let a2 = Self::adjust_end(point, end2, mid_end, other[0]);
                self.points.extend([mid1, a2, a1, mid1, mid2, a2]);
            }
            2 => {
                let a1 = Self::adjust_end(point, end1, mid_end, other[0]);
                let a2 = Self::adjust_end(point, end2, mid_end, other[1]);
                self.points
                    .extend([mid1, a1, mid_end, mid_end, mid2, mid1, mid_end, a2, mid2]);
            }
            _ => {}
        }
    }

    // ---- dashed wedge construction ----

    fn build_lines(&mut self) {
        let line_colour = self.line_colour;
        let line_width = self.line_width;
        let ShapeKind::DashedWedge {
            col2,
            one_less_dash,
            line_colours,
            at1,
            end1,
            end2,
        } = &mut self.kind
        else {
            return;
        };
        let (at1, end1, end2) = (*at1, *end1, *end2);
        let mid_end = mul(end1 + end2, 0.5);
        self.points.clear();
        line_colours.clear();
        let mut e1 = at1.direction_vector(&end1).unwrap_or_default();
        let mut e2 = at1.direction_vector(&end2).unwrap_or_default();
        // ACS1996: a 2.5px gap between lines, so allow for the line width
        let mut dash_sep = 2.5 + line_width;
        let central_len = (at1 - mid_end).length();
        let mut n_dashes = (central_len / dash_sep).round() as u32;
        // at least 3 dashes so the direction of the wedge can be seen
        let needed = if *one_less_dash { 4 } else { 3 };
        n_dashes = n_dashes.max(needed);
        // re-adjust so the last dash is on the end of the wedge
        dash_sep = central_len / f64::from(n_dashes);
        if *one_less_dash {
            // a shorter wedge that is just as wide at the end
            let endlenb2 = (end1 - end2).length() / 2.0;
            let central_line = at1.direction_vector(&mid_end).unwrap_or_default();
            let central_perp = Point2D::new(-central_line.y, central_line.x);
            let base = at1 + mul(central_line, central_len - dash_sep);
            let new_end1 = base + mul(central_perp, endlenb2);
            let new_end2 = base - mul(central_perp, endlenb2);
            e1 = at1.direction_vector(&new_end1).unwrap_or_default();
            e2 = at1.direction_vector(&new_end2).unwrap_or_default();
        }
        // similar triangles: the separation down the sides
        dash_sep *= (end1 - at1).length() / central_len;
        let extra = u32::from(!*one_less_dash);
        for i in 1..n_dashes + extra {
            self.points.push(at1 + mul(e1, f64::from(i) * dash_sep));
            self.points.push(at1 + mul(e2, f64::from(i) * dash_sep));
            line_colours.push(if i > n_dashes / 2 { *col2 } else { line_colour });
        }
    }
}
