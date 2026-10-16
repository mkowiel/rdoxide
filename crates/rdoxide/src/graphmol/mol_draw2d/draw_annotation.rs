//! Annotations and atom symbols as laid out text
//! (`DrawAnnotation.cpp` and `AtomSymbol.cpp` of RDKit 2026.09.1).
//!
//! Upstream objects keep a reference to the text drawer; here the drawer is
//! passed to the methods that need it. The `draw` methods return the SVG for
//! the glyphs, with the active class composed from the enclosing one.
//!
//! Copyright (C) 2021-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::draw_text::{DrawTextFT, StringLayout};
use super::helpers::{DrawColour, OrientType, TextAlignType};
use super::string_rect::StringRect;
use crate::Result;
use crate::geometry::Point2D;

/// `(xmin, xmax, ymin, ymax)`, to be initialised by the caller.
pub type Extremes = (f64, f64, f64, f64);

/// The initial value for [`Extremes`] that any real rectangle improves on.
pub const EMPTY_EXTREMES: Extremes = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);

fn compose_class(outer: &str, inner: &str) -> String {
    if outer.is_empty() {
        inner.to_string()
    } else {
        format!("{outer} {inner}")
    }
}

fn grow_by_rects(rects: &[StringRect], pos: Point2D, padding: f64, ext: &mut Extremes) {
    for r in rects {
        let mut r = *r;
        r.trans += pos;
        let [_, tr, _, bl] = r.calc_corners(padding);
        // the rect may be in a frame where +ve y is down or up: only the
        // ordering of ymin and ymax matters
        ext.0 = ext.0.min(tr.x).min(bl.x);
        ext.2 = ext.2.min(tr.y).min(bl.y);
        ext.1 = ext.1.max(tr.x).max(bl.x);
        ext.3 = ext.3.max(tr.y).max(bl.y);
    }
}

fn clashes(rects: &[StringRect], pos: Point2D, rect: &StringRect, padding: f64) -> bool {
    rects.iter().any(|r| {
        let mut r = *r;
        r.trans += pos;
        r.does_it_intersect(rect, padding)
    })
}

/// A note (atom or bond annotation, legend line, ...) on a single line.
#[derive(Clone, Debug)]
pub struct DrawAnnotation {
    pub text: String,
    pub align: TextAlignType,
    /// For SVG output.
    pub class: String,
    pub font_scale: f64,
    pub pos: Point2D,
    pub colour: DrawColour,
    pub rects: Vec<StringRect>,
}

impl DrawAnnotation {
    pub fn new(
        note: &str,
        align: TextAlignType,
        cls: &str,
        rel_font_scale: f64,
        pos: Point2D,
        colour: DrawColour,
        text: &mut DrawTextFT,
    ) -> Result<Self> {
        let mut a = DrawAnnotation {
            text: note.to_string(),
            align,
            class: cls.to_string(),
            font_scale: rel_font_scale * text.font_scale(),
            pos,
            colour,
            rects: Vec::new(),
        };
        a.extract_rects(text)?;
        Ok(a)
    }

    /// Notes are always one line of plain text, so no draw modes are kept.
    pub fn extract_rects(&mut self, text: &mut DrawTextFT) -> Result<()> {
        let ofs = text.font_scale();
        text.set_font_scale(self.font_scale, true);
        self.font_scale = text.font_scale();
        let layout = text.string_rects(&self.text, OrientType::C, true, self.align);
        text.set_font_scale(ofs, true);
        self.rects = layout?.rects;
        Ok(())
    }

    pub fn find_extremes(&self, ext: &mut Extremes, padding: f64) {
        if self.text.is_empty() {
            return;
        }
        grow_by_rects(&self.rects, self.pos, padding, ext);
    }

    pub fn dimensions(&self) -> (f64, f64) {
        let mut ext = EMPTY_EXTREMES;
        self.find_extremes(&mut ext, 0.0);
        (ext.1 - ext.0, ext.3 - ext.2)
    }

    /// The SVG for the note; `active_class` is the class already active.
    pub fn draw_svg(&self, text: &mut DrawTextFT, active_class: &str) -> Result<String> {
        let class = compose_class(active_class, &self.class);
        text.set_colour(self.colour);
        let ofs = text.font_scale();
        text.set_font_scale(self.font_scale, true);
        let svg = text.draw_string_svg(&self.text, self.pos, self.align, &class);
        text.set_font_scale(ofs, true);
        svg
    }

    pub fn scale(&mut self, scale_factor: Point2D, text: &mut DrawTextFT) -> Result<()> {
        self.pos.x *= scale_factor.x;
        self.pos.y *= scale_factor.y;
        // x scale for the fonts: x and y are very unlikely to differ
        self.font_scale *= scale_factor.x;
        // the widths needn't scale by the same amount as the font scale
        self.rects.clear();
        self.extract_rects(text)
    }

    /// The rects are relative to `pos`, so only that moves.
    pub fn translate(&mut self, trans: Point2D) {
        self.pos += trans;
    }

    pub fn does_rect_clash(&self, rect: &StringRect, padding: f64) -> bool {
        clashes(&self.rects, self.pos, rect, padding)
    }
}

/// The symbol of an atom, with the extra bits such as isotope labels.
#[derive(Clone, Debug)]
pub struct AtomSymbol {
    pub symbol: String,
    pub at_idx: i32,
    pub orient: OrientType,
    pub cds: Point2D,
    pub colour: DrawColour,
    pub layout: StringLayout,
}

impl AtomSymbol {
    pub fn new(
        symbol: &str,
        at_idx: i32,
        orient: OrientType,
        cds: Point2D,
        colour: DrawColour,
        text: &DrawTextFT,
    ) -> Result<Self> {
        let mut s = AtomSymbol {
            symbol: symbol.to_string(),
            at_idx,
            orient,
            cds,
            colour,
            layout: StringLayout::default(),
        };
        if !s.symbol.is_empty() {
            s.layout = text.string_rects(&s.symbol, orient, false, TextAlignType::Middle)?;
            s.adjust_colons()?;
        }
        Ok(s)
    }

    pub fn rects(&self) -> &[StringRect] {
        &self.layout.rects
    }

    pub fn find_extremes(&self, ext: &mut Extremes) {
        grow_by_rects(&self.layout.rects, self.cds, 0.0, ext);
    }

    /// Scales the position; the rects are rebuilt because the font scale may
    /// differ from that of the drawing as a whole.
    pub fn scale(&mut self, scale_factor: Point2D, text: &DrawTextFT) -> Result<()> {
        self.cds.x *= scale_factor.x;
        self.cds.y *= scale_factor.y;
        self.recalculate_rects(text)
    }

    pub fn translate(&mut self, trans: Point2D) {
        self.cds += trans;
    }

    pub fn recalculate_rects(&mut self, text: &DrawTextFT) -> Result<()> {
        self.layout = text.string_rects(&self.symbol, self.orient, false, TextAlignType::Middle)?;
        Ok(())
    }

    /// The SVG for the symbol, in class `atom-<idx>`.
    pub fn draw_svg(&self, text: &mut DrawTextFT, active_class: &str) -> Result<String> {
        let class = compose_class(active_class, &format!("atom-{}", self.at_idx));
        text.set_colour(self.colour);
        text.draw_label_svg(&self.symbol, self.cds, self.orient, &class)
    }

    pub fn does_rect_clash(&self, rect: &StringRect, padding: f64) -> bool {
        clashes(&self.layout.rects, self.cds, rect, padding)
    }

    /// A colon is a lot shorter than other characters, so a vertical bond can
    /// slip between an atom symbol and its map number (`C:8`). This makes the
    /// colon as tall as the smaller of the characters on either side.
    fn adjust_colons(&mut self) -> Result<()> {
        if self.symbol.is_empty() || !self.symbol.contains(':') {
            return Ok(());
        }
        // allow for markup such as <lit>[CH2;X2:4]</lit>
        let mut tmp = self.symbol.clone();
        while let Some(lt) = tmp.find('<') {
            let Some(gt) = tmp.find('>') else {
                return Ok(());
            };
            tmp = format!("{}{}", &tmp[..lt], tmp.get(gt + 1..).unwrap_or(""));
        }
        let Some(colon) = tmp.find(':') else {
            return Ok(());
        };
        let rects = &mut self.layout.rects;
        if colon >= rects.len() {
            return Err(crate::RDKitError::invariant("bad rects_ size"));
        }
        let left = if colon > 0 {
            rects[colon - 1].height
        } else {
            0.0
        };
        let right = rects.get(colon + 1).map_or(0.0, |r| r.height);
        rects[colon].height = left.min(right);
        Ok(())
    }
}
