//! Text layout for the molecule drawing code and its FreeType based text
//! drawer (`DrawText.cpp`, `DrawTextFT.cpp` and `DrawTextFTSVG.cpp` of RDKit
//! 2026.09.1).
//!
//! `DrawTextFT` lays strings out with the metrics of a [`TrueTypeFont`] and
//! reports each glyph as outline commands; [`DrawTextFT::svg_path_for_char`]
//! formats them as the SVG `<path>` elements RDKit writes.
//!
//! Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::fmt::Write;

use super::details::{do_lines_intersect, format_double};
use super::font::{PathCommand, TrueTypeFont};
use super::helpers::{DrawColour, OrientType, TextAlignType, colour_to_svg};
use super::string_rect::StringRect;
use crate::geometry::Point2D;
use crate::{RDKitError, Result, rd_log_warning};

/// How a character is drawn within a label.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TextDrawType {
    Normal,
    Superscript,
    Subscript,
}

/// Sub and superscripts are scaled by this.
const SUBS_SCALE: f64 = 0.66;
const SUPER_SCALE: f64 = 0.66;
/// A good default base font size.
pub const DEFAULT_FONT_SCALE: f64 = 0.6;

/// The result of laying out a string: one rectangle, draw mode and
/// character per glyph.
#[derive(Clone, Debug, Default)]
pub struct StringLayout {
    pub rects: Vec<StringRect>,
    pub draw_modes: Vec<TextDrawType>,
    pub draw_chars: Vec<u8>,
}

/// `setStringDrawMode`: if the text at `i` is `<sub>`, `<sup>` or their
/// closing tags, updates the mode, advances `i` to the last byte of the tag
/// and returns true.
pub fn set_string_draw_mode(text: &[u8], mode: &mut TextDrawType, i: &mut usize) -> bool {
    let rest = &text[*i..];
    if rest.starts_with(b"<sub>") {
        *mode = TextDrawType::Subscript;
        *i += 4;
    } else if rest.starts_with(b"<sup>") {
        *mode = TextDrawType::Superscript;
        *i += 4;
    } else if rest.starts_with(b"</sub>") || rest.starts_with(b"</sup>") {
        *mode = TextDrawType::Normal;
        *i += 5;
    } else {
        return false;
    }
    true
}

/// `atomLabelToPieces`: splits an atom label into the pieces that are laid
/// out separately for N and S orientations.
pub fn atom_label_to_pieces(label: &str, orient: OrientType) -> Vec<String> {
    let mut pieces: Vec<String> = Vec::new();
    if label.is_empty() {
        return pieces;
    }
    // <lit>XX</lit>: the symbol is used without modification
    if let Some(rest) = label.strip_prefix("<lit>") {
        let lit = rest.find("</lit>").map_or(rest, |idx| &rest[..idx]);
        pieces.push(lit.to_string());
        return pieces;
    }
    let bytes = label.as_bytes();
    let mut next_piece = String::new();
    let mut i = 0;
    loop {
        if i == bytes.len() {
            if !next_piece.is_empty() {
                pieces.push(next_piece);
            }
            break;
        }
        if bytes[i..].starts_with(b"<s") || bytes[i] == b':' || bytes[i].is_ascii_uppercase() {
            if !next_piece.is_empty() {
                pieces.push(std::mem::take(&mut next_piece));
            }
        }
        // labels are ASCII in practice; non-ASCII bytes stay in their char
        let ch_len = label[i..].chars().next().map_or(1, char::len_utf8);
        next_piece.push_str(&label[i..i + ch_len]);
        i += ch_len;
    }
    if pieces.len() < 2 {
        return pieces;
    }
    // for orientations S and E any charge flag goes at the end
    if matches!(orient, OrientType::E | OrientType::S) {
        if let Some(idx) = pieces
            .iter()
            .position(|p| p == "<sup>+</sup>" || p == "<sup>-</sup>")
        {
            let flag = std::mem::take(&mut pieces[idx]);
            pieces.push(flag);
        }
    }
    // Group the pieces: each atom symbol is flanked by its <sub>/<sup>.
    let mut final_pieces = Vec::new();
    let mut curr = String::new();
    let mut had_symbol = false;
    for p in pieces.iter().filter(|p| !p.is_empty()) {
        if !p.as_bytes()[0].is_ascii_uppercase() {
            curr.push_str(p);
        } else if had_symbol {
            final_pieces.push(std::mem::replace(&mut curr, p.clone()));
        } else {
            curr.push_str(p);
            had_symbol = true;
        }
    }
    if !curr.is_empty() {
        final_pieces.push(curr);
    }
    final_pieces
}

/// Font metrics and string layout, as in `DrawText` and `DrawTextFT`.
#[derive(Debug)]
pub struct DrawTextFT {
    font: TrueTypeFont,
    font_file: String,
    colour: DrawColour,
    font_scale: f64,
    max_font_size: f64,
    min_font_size: f64,
    base_font_size: f64,
}

impl DrawTextFT {
    /// `font_file` is empty or `BuiltinTelexRegular` for Telex,
    /// `BuiltinRobotoRegular` for Roboto, otherwise a TrueType file name.
    pub fn new(max_font_size: f64, min_font_size: f64, font_file: &str) -> Result<Self> {
        Ok(DrawTextFT {
            font: Self::load_font(font_file)?,
            font_file: font_file.to_string(),
            colour: DrawColour::new(0.0, 0.0, 0.0),
            font_scale: 1.0,
            max_font_size,
            min_font_size,
            base_font_size: DEFAULT_FONT_SCALE,
        })
    }

    fn load_font(font_file: &str) -> Result<TrueTypeFont> {
        match font_file {
            "" | "BuiltinTelexRegular" => Ok(TrueTypeFont::builtin_telex_regular()),
            "BuiltinRobotoRegular" => Ok(TrueTypeFont::builtin_roboto_regular()),
            path => TrueTypeFont::from_file(path),
        }
    }

    pub fn set_font_file(&mut self, font_file: &str) -> Result<()> {
        if font_file != self.font_file {
            self.font = Self::load_font(font_file)?;
            self.font_file = font_file.to_string();
        }
        Ok(())
    }

    pub fn font_file(&self) -> &str {
        &self.font_file
    }

    pub fn colour(&self) -> DrawColour {
        self.colour
    }
    pub fn set_colour(&mut self, col: DrawColour) {
        self.colour = col;
    }

    pub fn font_size(&self) -> f64 {
        self.font_scale * self.base_font_size
    }

    pub fn set_font_size(&mut self, new_size: f64) {
        if new_size < self.min_font_size {
            rd_log_warning!(
                "The new font size {} is below the current minimum ({}).\n",
                new_size,
                self.min_font_size
            );
        } else if new_size > self.max_font_size {
            rd_log_warning!(
                "The new font size {} is above the current maximum ({}).\n",
                new_size,
                self.max_font_size
            );
        }
        self.font_scale = new_size / self.base_font_size;
    }

    pub fn base_font_size(&self) -> f64 {
        self.base_font_size
    }
    pub fn set_base_font_size(&mut self, val: f64) {
        self.base_font_size = val;
    }
    pub fn max_font_size(&self) -> f64 {
        self.max_font_size
    }
    pub fn set_max_font_size(&mut self, val: f64) {
        self.max_font_size = val;
    }
    pub fn min_font_size(&self) -> f64 {
        self.min_font_size
    }
    pub fn set_min_font_size(&mut self, val: f64) {
        self.min_font_size = val;
    }
    pub fn font_scale(&self) -> f64 {
        self.font_scale
    }

    /// Sets the font scale; unless `ignore_limits`, the resulting size is kept
    /// within the min and max font sizes. Returns false if it was clamped.
    pub fn set_font_scale(&mut self, new_scale: f64, ignore_limits: bool) -> bool {
        self.font_scale = new_scale;
        if ignore_limits {
            return true;
        }
        let nfs = self.font_size();
        if self.max_font_size > 0.0 && nfs > self.max_font_size {
            self.font_scale = self.max_font_size / self.base_font_size;
            return false;
        }
        if self.min_font_size > 0.0 && nfs < self.min_font_size {
            self.font_scale = self.min_font_size / self.base_font_size;
            return false;
        }
        true
    }

    /// `fontCoordToDrawCoord`.
    fn font_coord_to_draw_coord(&self, fc: i32) -> f64 {
        self.font_size() * f64::from(fc) * self.font.em_scale()
    }

    fn select_scale_factor(c: u8, draw_type: TextDrawType) -> f64 {
        match draw_type {
            TextDrawType::Normal => 1.0,
            TextDrawType::Subscript => SUBS_SCALE,
            TextDrawType::Superscript if c == b'-' || c == b'+' => SUBS_SCALE,
            TextDrawType::Superscript => SUPER_SCALE,
        }
    }

    /// `FT_Load_Char` for a byte of the string: bytes above 127 are negative
    /// `char`s upstream and give the `.notdef` glyph.
    fn load_char(&self, c: u8) -> Result<super::font::Glyph> {
        self.font
            .load_char(if c < 128 { u32::from(c) } else { u32::MAX })
    }

    /// The rectangles of the characters of an unsplit string (`getStringRects`
    /// of `DrawTextFT`).
    fn raw_string_rects(&self, text: &str) -> Result<StringLayout> {
        let bytes = text.as_bytes();
        let mut layout = StringLayout::default();
        let mut draw_mode = TextDrawType::Normal;
        let mut max_y = 0.0f64;
        let mut mean_width = 0.0f64;
        let mut extras: Vec<f64> = Vec::new();
        let mut i = 0;
        while i < bytes.len() {
            if bytes[i] == b'<' && set_string_draw_mode(bytes, &mut draw_mode, &mut i) {
                i += 1;
                continue;
            }
            let c = bytes[i];
            i += 1;
            layout.draw_chars.push(c);
            let glyph = self.load_char(c)?;
            let bb = glyph.outline.bbox();
            let oscale = Self::select_scale_factor(c, draw_mode);
            let p_x_min = oscale * self.font_coord_to_draw_coord(bb.x_min);
            let p_y_min = oscale * self.font_coord_to_draw_coord(bb.y_min);
            let p_x_max = oscale * self.font_coord_to_draw_coord(bb.x_max);
            let p_y_max = oscale * self.font_coord_to_draw_coord(bb.y_max);
            let p_advance = oscale * self.font_coord_to_draw_coord(glyph.advance);
            let mut width = p_x_max - p_x_min;
            // the mean width defines the spacing, so use full size characters
            mean_width += width / oscale;
            // reduce the horizontal offset, otherwise spacing is uneven
            extras.push(p_advance - p_x_max);
            if bb.x_max == 0 {
                // probably a space: keep it small
                width = p_advance / 3.0;
            }
            let height = p_y_max - p_y_min;
            let offset = Point2D::new(p_x_min + width / 2.0, p_y_max / 2.0);
            let g_centre = Point2D::new(offset.x, p_y_max - height / 2.0);
            layout
                .rects
                .push(StringRect::new(offset, g_centre, width, height));
            layout.draw_modes.push(draw_mode);
            max_y = max_y.max(p_y_max);
        }
        if layout.rects.is_empty() {
            return Ok(layout);
        }
        mean_width /= layout.rects.len() as f64;
        for i in 0..layout.rects.len() {
            extras[i] += mean_width / 10.0;
            layout.rects[i].g_centre.y = max_y - layout.rects[i].g_centre.y;
            layout.rects[i].offset.y = max_y / 2.0;
            if i > 0 {
                let prev = layout.rects[i - 1];
                layout.rects[i].trans.x =
                    prev.trans.x + prev.width / 2.0 + layout.rects[i].width / 2.0 + extras[i - 1];
            }
        }
        Self::adjust_for_super_sub_script(&layout.draw_modes, &mut layout.rects);
        Ok(layout)
    }

    fn adjust_for_super_sub_script(draw_modes: &[TextDrawType], rects: &mut [StringRect]) {
        let mut last_char: Option<usize> = None;
        for i in 0..draw_modes.len() {
            match draw_modes[i] {
                TextDrawType::Superscript => {
                    // superscripts may come before a letter: look ahead
                    if last_char.is_none() {
                        last_char = (i + 1..draw_modes.len())
                            .find(|&j| draw_modes[j] == TextDrawType::Normal);
                    }
                    let h = rects[last_char.unwrap_or(0)].height;
                    rects[i].rect_corr = h;
                    rects[i].trans.y -= h / 2.0;
                }
                TextDrawType::Subscript => {
                    let h = rects[last_char.unwrap_or(0)].height;
                    rects[i].rect_corr = -h;
                    rects[i].trans.y -= rects[i].rect_corr / 2.0;
                }
                TextDrawType::Normal => last_char = Some(i),
            }
        }
        // contiguous sub- and superscripts go at the same x
        for i in 1..rects.len() {
            if draw_modes[i] != draw_modes[i - 1]
                && draw_modes[i] != TextDrawType::Normal
                && draw_modes[i - 1] != TextDrawType::Normal
            {
                let move_back = rects[i].trans.x - rects[i - 1].trans.x;
                for r in &mut rects[i..] {
                    r.trans.x -= move_back;
                }
            }
        }
    }

    /// `alignString`: shifts the rectangles so the alignment point is at 0.
    pub fn align_string(
        &self,
        mut talign: TextAlignType,
        draw_modes: &[TextDrawType],
        rects: &mut [StringRect],
    ) {
        if talign == TextAlignType::Middle
            && draw_modes
                .iter()
                .filter(|&&m| m == TextDrawType::Normal)
                .count()
                == 1
        {
            talign = TextAlignType::Start;
        }
        let mut align_trans = Point2D::new(0.0, 0.0);
        if talign == TextAlignType::Start || talign == TextAlignType::End {
            let mut align_char = 0;
            for (i, m) in draw_modes.iter().enumerate().take(rects.len()) {
                if *m == TextDrawType::Normal {
                    align_char = i;
                    if talign == TextAlignType::Start {
                        break;
                    }
                }
            }
            align_trans.x = rects[0].trans.x - rects[align_char].trans.x;
        } else {
            // centre on the middle of the Normal text; the super- or
            // subscripts should be at the ends
            let (mut x_min, mut x_max) = (f64::MAX, f64::MIN);
            let mut align_char: Option<usize> = None;
            for (i, r) in rects.iter().enumerate() {
                if draw_modes[i] == TextDrawType::Normal {
                    align_char = align_char.or(Some(i));
                    let [_, tr, _, bl] = r.calc_corners(0.0);
                    x_min = x_min.min(bl.x).min(tr.x);
                    x_max = x_max.max(bl.x).max(tr.x);
                }
            }
            let align_char = align_char.unwrap_or(0);
            align_trans.x = -(x_max - x_min - rects[align_char].width) / 2.0;
        }
        for r in rects.iter_mut() {
            r.trans += align_trans;
        }
    }

    /// `DrawText::getStringRects` with orientation: splits the label into
    /// pieces and lays them out for the orientation.
    pub fn string_rects(
        &self,
        text: &str,
        orient: OrientType,
        dont_split: bool,
        text_align: TextAlignType,
    ) -> Result<StringLayout> {
        if text.is_empty() {
            return Err(RDKitError::precondition("empty string"));
        }
        let bits = if dont_split {
            vec![text.to_string()]
        } else {
            atom_label_to_pieces(text, orient)
        };
        let ta = if orient == OrientType::C {
            text_align
        } else {
            TextAlignType::Middle
        };
        match orient {
            OrientType::W => {
                // stick the pieces together backwards and draw as one
                let new_lab: String = bits.iter().rev().map(String::as_str).collect();
                let mut l = self.raw_string_rects(&new_lab)?;
                self.align_string(TextAlignType::End, &l.draw_modes, &mut l.rects);
                Ok(l)
            }
            OrientType::E => {
                let new_lab: String = bits.concat();
                let mut l = self.raw_string_rects(&new_lab)?;
                self.align_string(TextAlignType::Start, &l.draw_modes, &mut l.rects);
                Ok(l)
            }
            _ => {
                let mut layout = StringLayout::default();
                let mut running_y = 0.0;
                for tb in &bits {
                    let mut t = self.raw_string_rects(tb)?;
                    self.align_string(ta, &t.draw_modes, &mut t.rects);
                    let mut max_height = f64::MIN;
                    for r in &mut t.rects {
                        max_height = max_height.max(r.height);
                        r.y_shift = running_y;
                    }
                    layout.rects.extend(t.rects);
                    layout.draw_modes.extend(t.draw_modes);
                    layout.draw_chars.extend(t.draw_chars);
                    if orient == OrientType::N {
                        running_y -= 1.1 * max_height;
                    } else if orient == OrientType::S {
                        running_y += 1.1 * max_height;
                    }
                }
                Ok(layout)
            }
        }
    }

    /// `getStringExtremes`: `(x_min, y_min, x_max, y_max)` of a label.
    pub fn string_extremes(
        &self,
        label: &str,
        orient: OrientType,
        dont_split: bool,
    ) -> Result<(f64, f64, f64, f64)> {
        if label.is_empty() {
            return Ok((0.0, 0.0, 0.0, 0.0));
        }
        let layout = self.string_rects(label, orient, dont_split, TextAlignType::Middle)?;
        let (mut x_min, mut y_min) = (f64::MAX, f64::MAX);
        let (mut x_max, mut y_max) = (f64::MIN, f64::MIN);
        for r in &layout.rects {
            let [_, tr, _, bl] = r.calc_corners(0.0);
            x_min = x_min.min(bl.x).min(tr.x);
            y_min = y_min.min(bl.y).min(tr.y);
            x_max = x_max.max(bl.x).max(tr.x);
            y_max = y_max.max(bl.y).max(tr.y);
        }
        Ok((x_min, y_min, x_max, y_max))
    }

    /// `getStringSize`: width and height of a label.
    pub fn string_size(&self, label: &str) -> Result<(f64, f64)> {
        let (x_min, y_min, x_max, y_max) = self.string_extremes(label, OrientType::E, false)?;
        Ok((x_max - x_min, y_max - y_min))
    }

    /// `adjustLineForString`: moves `end2` to where the line from `end1`
    /// meets the label drawn at `end2`.
    pub fn adjust_line_for_string(
        &self,
        label: &str,
        orient: OrientType,
        end1: Point2D,
        end2: &mut Point2D,
    ) -> Result<()> {
        let lab_pos = *end2;
        let layout = self.string_rects(label, orient, false, TextAlignType::Middle)?;
        let bond_len = (end1 - *end2).length();
        for r in &layout.rects {
            let mut r = *r;
            r.trans += lab_pos;
            let [tl, tr, br, bl] = r.calc_corners(0.025 * bond_len);
            for (a, b) in [(tl, tr), (tr, br), (br, bl), (bl, tl)] {
                if let Some(ip) = do_lines_intersect(*end2, end1, a, b) {
                    *end2 = ip;
                }
            }
        }
        Ok(())
    }

    /// `doesRectIntersect` for a label drawn at `cds`.
    pub fn does_rect_intersect(
        &self,
        label: &str,
        orient: OrientType,
        cds: Point2D,
        rect: &StringRect,
    ) -> Result<bool> {
        if label.is_empty() {
            return Ok(false);
        }
        let layout = self.string_rects(label, orient, false, TextAlignType::Middle)?;
        Ok(layout.rects.iter().any(|r| {
            let mut nr = *r;
            nr.trans += cds;
            nr.does_it_intersect(rect, 0.0)
        }))
    }

    /// `doesLineIntersect` for a label drawn at `cds`.
    pub fn does_line_intersect(
        &self,
        label: &str,
        orient: OrientType,
        cds: Point2D,
        end1: Point2D,
        end2: Point2D,
        padding: f64,
    ) -> Result<bool> {
        let layout = self.string_rects(label, orient, false, TextAlignType::Middle)?;
        Ok(layout.rects.iter().any(|r| {
            let mut nr = *r;
            nr.trans += cds;
            let [tl, tr, br, bl] = nr.calc_corners(padding);
            [(tl, tr), (tr, br), (br, bl), (bl, tl)]
                .iter()
                .any(|&(a, b)| do_lines_intersect(end2, end1, a, b).is_some())
        }))
    }

    /// `doesStringIntersect`: do two labels drawn at the given points overlap?
    pub fn does_string_intersect(
        &self,
        label1: &str,
        orient1: OrientType,
        cds1: Point2D,
        label2: &str,
        orient2: OrientType,
        cds2: Point2D,
    ) -> Result<bool> {
        if label1.is_empty() || label2.is_empty() {
            return Ok(false);
        }
        let l1 = self.string_rects(label1, orient1, false, TextAlignType::Middle)?;
        let l2 = self.string_rects(label2, orient2, false, TextAlignType::Middle)?;
        for r1 in &l1.rects {
            let mut nr1 = *r1;
            nr1.trans += cds1;
            for r2 in &l2.rects {
                let mut nr2 = *r2;
                nr2.trans += cds2;
                if nr1.does_it_intersect(&nr2, 0.0) {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// The position at which each character of a laid out string is drawn
    /// (`drawChars`), with the font scale to draw it at.
    pub fn char_positions(&self, a_cds: Point2D, layout: &StringLayout) -> Vec<(u8, Point2D, f64)> {
        let full_scale = self.font_scale;
        layout
            .rects
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let mut p = Point2D::new(
                    a_cds.x + r.trans.x - r.offset.x,
                    a_cds.y - r.trans.y + r.offset.y,
                );
                p.y -= r.rect_corr + r.y_shift;
                let c = layout.draw_chars[i];
                let scale = full_scale * Self::select_scale_factor(c, layout.draw_modes[i]);
                (c, p, scale)
            })
            .collect()
    }

    /// The outline commands of the character `c` drawn at `cds` with the
    /// current font size (`fontPosToDrawPos` applied to every point), and the
    /// advance.
    pub fn char_outline(&self, c: u8, cds: Point2D) -> Result<(Vec<PathCommand2D>, f64)> {
        let glyph = self.load_char(c)?;
        let to_draw = |x: i32, y: i32| {
            // FreeType has the origin bottom left, we want it top left
            (
                cds.x + self.font_coord_to_draw_coord(x),
                cds.y + self.font_coord_to_draw_coord(-y),
            )
        };
        let cmds = glyph
            .outline
            .decompose()
            .into_iter()
            .map(|cmd| match cmd {
                PathCommand::MoveTo(x, y) => PathCommand2D::MoveTo(to_draw(x, y)),
                PathCommand::LineTo(x, y) => PathCommand2D::LineTo(to_draw(x, y)),
                PathCommand::ConicTo(c, t) => {
                    PathCommand2D::ConicTo(to_draw(c.0, c.1), to_draw(t.0, t.1))
                }
            })
            .collect();
        Ok((cmds, self.font_coord_to_draw_coord(glyph.advance)))
    }

    /// The SVG `<path>` element for a character (`DrawTextFTSVG`).
    pub fn svg_path_for_char(&self, c: u8, cds: Point2D, active_class: &str) -> Result<String> {
        let (cmds, _) = self.char_outline(c, cds)?;
        let mut s = String::from("<path ");
        if active_class.is_empty() {
            s.push_str("d='");
        } else {
            let _ = write!(s, "class='{active_class}' d='");
        }
        for cmd in cmds {
            match cmd {
                PathCommand2D::MoveTo(p) => {
                    let _ = writeln!(s, "M {} {}", format_double(p.0), format_double(p.1));
                }
                PathCommand2D::LineTo(p) => {
                    let _ = writeln!(s, "L {} {}", format_double(p.0), format_double(p.1));
                }
                PathCommand2D::ConicTo(ctl, p) => {
                    let _ = writeln!(
                        s,
                        "Q {} {}, {} {}",
                        format_double(ctl.0),
                        format_double(ctl.1),
                        format_double(p.0),
                        format_double(p.1)
                    );
                }
            }
        }
        let _ = writeln!(s, "' fill='{}'/>", colour_to_svg(&self.colour)?);
        Ok(s)
    }

    /// `drawString` with an alignment: the SVG for all characters of `text`.
    pub fn draw_string_svg(
        &mut self,
        text: &str,
        cds: Point2D,
        talign: TextAlignType,
        active_class: &str,
    ) -> Result<String> {
        let mut layout = self.raw_string_rects(text)?;
        self.align_string(talign, &layout.draw_modes, &mut layout.rects);
        self.layout_to_svg(cds, &layout, active_class)
    }

    /// `drawString` with an orientation.
    pub fn draw_label_svg(
        &mut self,
        label: &str,
        cds: Point2D,
        orient: OrientType,
        active_class: &str,
    ) -> Result<String> {
        let layout = self.string_rects(label, orient, false, TextAlignType::Middle)?;
        self.layout_to_svg(cds, &layout, active_class)
    }

    fn layout_to_svg(
        &mut self,
        cds: Point2D,
        layout: &StringLayout,
        active_class: &str,
    ) -> Result<String> {
        let full_scale = self.font_scale;
        let mut out = String::new();
        for (c, p, scale) in self.char_positions(cds, layout) {
            self.set_font_scale(scale, true);
            let path = self.svg_path_for_char(c, p, active_class);
            self.set_font_scale(full_scale, true);
            out.push_str(&path?);
        }
        Ok(out)
    }
}

/// A point of a drawn glyph outline, in draw coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PathCommand2D {
    MoveTo((f64, f64)),
    LineTo((f64, f64)),
    ConicTo((f64, f64), (f64, f64)),
}
