//! The SVG writing half of `MolDraw2DSVG` (`MolDraw2DSVG.cpp` of RDKit
//! 2026.09.1): the drawing state and the primitives that turn into SVG
//! elements. All coordinates given to it are drawing coordinates.
//!
//! Copyright (C) 2015-2021 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::details::{calc_arrow_head, format_double, get_wavy_line_segments};
use super::draw_mol_render::MolCanvas;
use super::draw_shape::ShapeCanvas;
use super::helpers::{
    DashPattern, DrawColour, LINE_WIDTH_SCALE_FACTOR, MolDrawOptions, colour_to_svg,
};
use crate::geometry::Point2D;
use crate::rdgeneral::strings::format_g;
use crate::{RDKitError, Result};

/// A `double` as a C++ output stream writes it.
pub(super) fn stream_double(val: f64) -> String {
    format_g(val, 6)
}

/// The state of a drawing and the SVG written so far.
#[derive(Debug)]
pub struct SvgCanvas {
    pub options: MolDrawOptions,
    pub(super) out: String,
    pub(super) colour: DrawColour,
    pub(super) dash: DashPattern,
    pub(super) fill_polys: bool,
    pub(super) scale: f64,
    pub(super) active_class: String,
    pub(super) active_atm_idx1: i32,
    pub(super) active_atm_idx2: i32,
    pub(super) active_bnd_idx: i32,
}

impl SvgCanvas {
    pub fn new(options: MolDrawOptions) -> Self {
        SvgCanvas {
            options,
            out: String::new(),
            colour: DrawColour::new(0.0, 0.0, 0.0),
            dash: Vec::new(),
            fill_polys: true,
            scale: 1.0,
            active_class: String::new(),
            active_atm_idx1: -1,
            active_atm_idx2: -1,
            active_bnd_idx: -1,
        }
    }

    /// `getDrawLineWidth`: the line width in pixels, scaled if asked.
    pub fn draw_line_width(&self) -> f64 {
        let mut width = self.options.bond_line_width;
        if self.options.scale_bond_width {
            width *= self.scale * LINE_WIDTH_SCALE_FACTOR;
            if width < 0.0 {
                width = 0.0;
            }
        }
        width
    }

    pub(super) fn has_active_atm_idx(&self) -> bool {
        self.active_atm_idx1 >= 0
    }

    /// The `class` attribute for the active class, bond and atoms, if any.
    fn output_classes(&mut self) {
        if self.active_class.is_empty() && !self.has_active_atm_idx() {
            return;
        }
        self.out.push_str("class='");
        self.out.push_str(&self.active_class);
        if self.active_bnd_idx >= 0 {
            self.out.push_str(&format!("bond-{}", self.active_bnd_idx));
        }
        if !self.has_active_atm_idx() {
            self.out.push_str("' ");
            return;
        }
        if self.active_bnd_idx >= 0 {
            self.out.push(' ');
        }
        if !self.active_class.is_empty() {
            self.out.push(' ');
        }
        self.out.push_str(&format!("atom-{}", self.active_atm_idx1));
        if self.active_atm_idx2 >= 0 && self.active_atm_idx2 != self.active_atm_idx1 {
            self.out
                .push_str(&format!(" atom-{}", self.active_atm_idx2));
        }
        self.out.push_str("' ");
    }

    fn dash_string(&self) -> String {
        match self.dash.split_last() {
            None => String::new(),
            Some((last, rest)) => {
                let mut res = String::from(";stroke-dasharray:");
                for d in rest {
                    res.push_str(&stream_double(*d));
                    res.push(',');
                }
                res.push_str(&stream_double(*last));
                res
            }
        }
    }

    fn line(&mut self, c1: Point2D, c2: Point2D) -> Result<()> {
        let col = colour_to_svg(&self.colour)?;
        let width = self.draw_line_width();
        let dash = self.dash_string();
        self.out.push_str("<path ");
        self.output_classes();
        self.out.push_str(&format!(
            "d='M {},{} L {},{}' ",
            format_double(c1.x),
            format_double(c1.y),
            format_double(c2.x),
            format_double(c2.y)
        ));
        self.out.push_str(&format!(
            "style='fill:none;fill-rule:evenodd;stroke:{col};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1{dash}'",
            format_double(width)
        ));
        self.out.push_str(" />\n");
        Ok(())
    }

    fn polygon(&mut self, cds: &[Point2D]) -> Result<()> {
        if cds.len() < 3 {
            return Err(RDKitError::precondition("must have at least three points"));
        }
        let col = colour_to_svg(&self.colour)?;
        let width = self.draw_line_width();
        let dash = self.dash_string();
        let alpha = stream_double(self.colour.a);
        self.out.push_str("<path ");
        self.output_classes();
        self.out.push('d');
        self.out.push_str("='M");
        self.out.push_str(&format!(
            " {},{}",
            format_double(cds[0].x),
            format_double(cds[0].y)
        ));
        for c in &cds[1..] {
            self.out
                .push_str(&format!(" L {},{}", format_double(c.x), format_double(c.y)));
        }
        if self.fill_polys {
            // the Z closes the path which we don't want for unfilled polygons
            self.out.push_str(&format!(
                " Z' style='fill:{col};fill-rule:evenodd;fill-opacity:{alpha};"
            ));
        } else {
            self.out.push_str("' style='fill:none;");
        }
        self.out.push_str(&format!(
            "stroke:{col};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:{alpha};{dash}'",
            format_double(width)
        ));
        self.out.push_str(" />\n");
        Ok(())
    }

    fn ellipse(&mut self, c1: Point2D, c2: Point2D) -> Result<()> {
        let mut w = c2.x - c1.x;
        let mut h = c2.y - c1.y;
        let cx = c1.x + w / 2.0;
        let cy = c1.y + h / 2.0;
        w = w.abs();
        h = h.abs();
        let col = colour_to_svg(&self.colour)?;
        let width = self.draw_line_width();
        let dash = self.dash_string();
        self.out.push_str(&format!(
            "<ellipse cx='{}' cy='{}' rx='{}' ry='{}' ",
            format_double(cx),
            format_double(cy),
            format_double(w / 2.0),
            format_double(h / 2.0)
        ));
        self.output_classes();
        self.out.push_str(" style='");
        if self.fill_polys {
            self.out.push_str(&format!("fill:{col};fill-rule:evenodd;"));
        } else {
            self.out.push_str("fill:none;");
        }
        self.out.push_str(&format!(
            "stroke:{col};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1{dash}'",
            format_double(width)
        ));
        self.out.push_str(" />\n");
        Ok(())
    }

    fn wavy_line(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col: DrawColour,
        n_segments: u32,
        vert_offset: f64,
    ) -> Result<()> {
        if n_segments <= 1 {
            return Err(RDKitError::precondition("too few segments"));
        }
        self.colour = col;
        let segments = get_wavy_line_segments(cds1, cds2, n_segments, vert_offset)?;
        let col = colour_to_svg(&self.colour)?;
        let width = self.draw_line_width();
        self.out.push_str("<path ");
        self.output_classes();
        let c1 = segments[0].0;
        self.out.push_str(&format!(
            "d='M{},{}",
            stream_double(c1.x),
            stream_double(c1.y)
        ));
        for (_, cpt1, cpt2, segpt) in segments.iter().take(n_segments as usize) {
            self.out.push_str(&format!(
                " C{},{} {},{} {},{}",
                stream_double(cpt1.x),
                stream_double(cpt1.y),
                stream_double(cpt2.x),
                stream_double(cpt2.y),
                stream_double(segpt.x),
                stream_double(segpt.y)
            ));
        }
        self.out.push_str("' ");
        self.out.push_str(&format!(
            "style='fill:none;stroke:{col};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'",
            format_double(width)
        ));
        self.out.push_str(" />\n");
        Ok(())
    }
}

impl ShapeCanvas for SvgCanvas {
    fn colour(&self) -> DrawColour {
        self.colour
    }
    fn set_colour(&mut self, col: DrawColour) {
        self.colour = col;
    }
    fn line_width(&self) -> f64 {
        self.options.bond_line_width
    }
    fn set_line_width(&mut self, width: f64) {
        self.options.bond_line_width = width;
    }
    fn fill_polys(&self) -> bool {
        self.fill_polys
    }
    fn set_fill_polys(&mut self, fill: bool) {
        self.fill_polys = fill;
    }
    fn scale_bond_width(&self) -> bool {
        self.options.scale_bond_width
    }
    fn set_scale_bond_width(&mut self, scale: bool) {
        self.options.scale_bond_width = scale;
    }
    fn split_bonds(&self) -> bool {
        self.options.split_bonds
    }
    fn dash(&self) -> DashPattern {
        self.dash.clone()
    }
    fn set_dash(&mut self, dash: DashPattern) {
        self.dash = dash;
    }
    fn set_active_atm_idx(&mut self, atom1: i32, atom2: i32) {
        let atom1 = atom1.max(-1);
        let atom2 = atom2.max(-1);
        if atom2 >= 0 && atom1 < 0 {
            self.active_atm_idx1 = atom2;
            self.active_atm_idx2 = atom1;
        } else {
            self.active_atm_idx1 = atom1;
            self.active_atm_idx2 = atom2;
        }
    }
    fn set_active_bnd_idx(&mut self, bond: i32) {
        self.active_bnd_idx = bond.max(-1);
    }

    fn draw_arrow(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        as_polygon: bool,
        frac: f64,
        angle: f64,
        col: DrawColour,
        raw_coords: bool,
    ) -> Result<()> {
        if angle <= 1.0e-6 {
            return Err(RDKitError::precondition(
                "Arrow angle must be positive and greater than 0.0.",
            ));
        }
        let mut ae = cds2;
        let (p1, p2) = calc_arrow_head(&mut ae, cds1, frac, self.draw_line_width(), angle)
            .ok_or_else(|| RDKitError::precondition("bad arrow"))?;
        self.draw_line(cds1, ae, col, col, raw_coords)?;
        let fps = self.fill_polys;
        let dc = self.colour;
        self.fill_polys = as_polygon;
        self.colour = col;
        let res = self.polygon(&[p1, ae, p2]);
        self.fill_polys = fps;
        self.colour = dc;
        res
    }

    fn draw_ellipse(&mut self, cds1: Point2D, cds2: Point2D, _raw_coords: bool) -> Result<()> {
        self.ellipse(cds1, cds2)
    }

    fn draw_line(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col1: DrawColour,
        col2: DrawColour,
        _raw_coords: bool,
    ) -> Result<()> {
        if col1 == col2 {
            self.colour = col1;
            self.line(cds1, cds2)
        } else {
            let mid = (cds1 + cds2) * 0.5;
            self.colour = col1;
            self.line(cds1, mid)?;
            self.colour = col2;
            self.line(mid, cds2)
        }
    }

    fn draw_polygon(&mut self, cds: &[Point2D], _raw_coords: bool) -> Result<()> {
        self.polygon(cds)
    }

    fn draw_triangle(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        cds3: Point2D,
        _raw_coords: bool,
    ) -> Result<()> {
        self.polygon(&[cds1, cds2, cds3])
    }

    fn draw_wavy_line(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col1: DrawColour,
        _col2: DrawColour,
        n_segments: u32,
        vert_offset: f64,
        _raw_coords: bool,
    ) -> Result<()> {
        self.wavy_line(cds1, cds2, col1, n_segments, vert_offset)
    }

    fn draw_arc(
        &mut self,
        centre: Point2D,
        xradius: f64,
        yradius: f64,
        angle1: f64,
        angle2: f64,
        _raw_coords: bool,
    ) -> Result<()> {
        // 5 degree increments should be plenty, as the circles are probably
        // going to be small.
        let num_steps = 1 + ((angle2 - angle1) / 5.0) as i32;
        let ang_incr = (angle2 - angle1) / f64::from(num_steps) * std::f64::consts::PI / 180.0;
        let start = angle1 * std::f64::consts::PI / 180.0;
        let mut pts: Vec<Point2D> = (0..=num_steps)
            .map(|i| {
                let ang = start + f64::from(i) * ang_incr;
                Point2D::new(
                    centre.x + xradius * ang.cos(),
                    centre.y + yradius * ang.sin(),
                )
            })
            .collect();
        if self.fill_polys {
            // otherwise it draws an arc back to the first point rather than
            // filling in the sector.
            pts.push(centre);
        }
        // Very short arcs can be drawn as lines.
        if pts.len() == 2 {
            self.line(pts[0], pts[1])
        } else {
            self.polygon(&pts)
        }
    }
}

impl MolCanvas for SvgCanvas {
    fn scale(&self) -> f64 {
        self.scale
    }
    fn set_scale(&mut self, scale: f64) {
        self.scale = scale;
    }
    fn active_class(&self) -> String {
        self.active_class.clone()
    }
    fn write_text(&mut self, svg: &str) {
        self.out.push_str(svg);
    }
}
