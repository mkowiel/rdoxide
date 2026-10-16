//! The rectangle round a glyph or string used when laying out text
//! (`StringRect.h` of RDKit 2026.09.1).
//!
//! Copyright (C) 2021-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::geometry::Point2D;

/// Dimensions of the rectangle round a string.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StringRect {
    /// Where to draw the char relative to other chars in the string.
    pub trans: Point2D,
    /// Offset for draw coords so the char is centred correctly.
    pub offset: Point2D,
    /// Glyph centre relative to the origin of the char.
    pub g_centre: Point2D,
    /// Shift the whole thing in y by this. For multi-line text.
    pub y_shift: f64,
    /// Width and height of the glyph itself, not of the character cell.
    pub width: f64,
    pub height: f64,
    /// If a char is moved one way, the rectangle must move the other.
    pub rect_corr: f64,
    /// Rough measure of how badly it clashed with other things; lower is
    /// better, 0 is no clash.
    pub clash_score: i32,
}

type Corners = [Point2D; 4];

impl StringRect {
    pub fn new(offset: Point2D, g_centre: Point2D, width: f64, height: f64) -> Self {
        StringRect {
            trans: Point2D::new(0.0, 0.0),
            offset,
            g_centre,
            y_shift: 0.0,
            width,
            height,
            rect_corr: 0.0,
            clash_score: 0,
        }
    }

    /// The corners `[tl, tr, br, bl]` of the glyph box grown by `padding`.
    pub fn calc_corners(&self, padding: f64) -> Corners {
        let wb2 = padding + self.width / 2.0;
        let hb2 = padding + self.height / 2.0;
        let c = self.calc_centre();
        [
            Point2D::new(c.x - wb2, c.y - hb2),
            Point2D::new(c.x + wb2, c.y - hb2),
            Point2D::new(c.x + wb2, c.y + hb2),
            Point2D::new(c.x - wb2, c.y + hb2),
        ]
    }

    pub fn calc_centre(&self) -> Point2D {
        let mut c = self.trans + self.g_centre - self.offset;
        c.y -= self.y_shift;
        c
    }

    /// Corners ordered so that "top" has the larger y, whichever way y points.
    fn oriented_corners(&self, padding: f64) -> Corners {
        let [mut tl, mut tr, mut br, mut bl] = self.calc_corners(padding);
        if tl.y < bl.y {
            std::mem::swap(&mut tl, &mut bl);
            std::mem::swap(&mut tr, &mut br);
        }
        [tl, tr, br, bl]
    }

    pub fn is_point_inside(&self, pt: Point2D, padding: f64) -> bool {
        let [tl, _, br, _] = self.oriented_corners(padding);
        pt.x >= tl.x && pt.x <= br.x && pt.y >= br.y && pt.y <= tl.y
    }

    pub fn does_it_intersect(&self, other: &StringRect, padding: f64) -> bool {
        let t = self.oriented_corners(padding);
        let o = other.oriented_corners(padding);
        // Corner containment in either direction. `c[0]` is the top left and
        // `c[3]` the bottom left, `c[1]` the top right.
        let inside = |p: &Point2D, r: &Corners| {
            p.x >= r[0].x && p.x <= r[1].x && p.y >= r[3].y && p.y <= r[0].y
        };
        o.iter().any(|p| inside(p, &t)) || t.iter().any(|p| inside(p, &o))
    }
}
