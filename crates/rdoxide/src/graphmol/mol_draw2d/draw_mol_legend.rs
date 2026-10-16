//! The legend of a [`DrawMol`] (`partitionForLegend`, `extractLegend` and
//! helpers of `DrawMol.cpp`, RDKit 2026.09.1).
//!
//! Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::draw_annotation::{DrawAnnotation, EMPTY_EXTREMES};
use super::draw_mol::DrawMol;
use super::draw_text::DrawTextFT;
use super::helpers::{DrawColour, LegendPosition, TextAlignType};
use crate::Result;
use crate::geometry::Point2D;

fn legend_note(
    bit: &str,
    rel_font_scale: f64,
    pos: Point2D,
    colour: DrawColour,
    text: &mut DrawTextFT,
) -> Result<DrawAnnotation> {
    DrawAnnotation::new(
        bit,
        TextAlignType::Middle,
        "legend",
        rel_font_scale,
        pos,
        colour,
        text,
    )
}

/// `(width, height)` of one line of a legend at the relative font scale.
fn bit_dimensions(
    bit: &str,
    rel_font_scale: f64,
    colour: DrawColour,
    text: &mut DrawTextFT,
) -> Result<(f64, f64)> {
    Ok(legend_note(bit, rel_font_scale, Point2D::new(0.0, 0.0), colour, text)?.dimensions())
}

/// The widest line and total height of the lines.
fn calc_legend_dims(
    bits: &[String],
    rel_font_scale: f64,
    colour: DrawColour,
    text: &mut DrawTextFT,
) -> Result<(f64, f64)> {
    let (mut total_width, mut total_height) = (0.0f64, 0.0);
    for bit in bits {
        let (w, h) = bit_dimensions(bit, rel_font_scale, colour, text)?;
        total_height += h;
        total_width = total_width.max(w);
    }
    Ok((total_width, total_height))
}

fn calc_line_gap(text: &DrawTextFT, vert_text: bool, rel_font_scale: f64) -> f64 {
    let font_px = text.font_size() * rel_font_scale;
    if vert_text {
        // vertical text needs more inter-glyph room
        (font_px * 0.35).max(1.0)
    } else {
        // horizontal lines can be tighter, with a 2 px minimum
        (font_px * 0.15).max(2.0)
    }
}

fn split_legend_bits(legend: &str, vert_text: bool) -> Vec<String> {
    if vert_text {
        legend
            .chars()
            .filter(|&c| c != '\n')
            .map(String::from)
            .collect()
    } else {
        legend
            .split('\n')
            .filter(|s| !s.is_empty())
            .map(String::from)
            .collect()
    }
}

struct LegendPlacement {
    base_x: f64,
    base_y: f64,
    draw_width: f64,
    draw_height: f64,
    rel_font_scale: f64,
    colour: DrawColour,
}

fn place_bottom_legend(
    bits: &[String],
    p: &LegendPlacement,
    text: &mut DrawTextFT,
    legends: &mut Vec<DrawAnnotation>,
) -> Result<()> {
    let loc = Point2D::new(p.base_x + p.draw_width / 2.0, p.base_y + p.draw_height);
    for bit in bits {
        legends.push(legend_note(bit, p.rel_font_scale, loc, p.colour, text)?);
    }
    let mut last_above = 0.0;
    for i in (0..legends.len()).rev() {
        let mut ext = EMPTY_EXTREMES;
        legends[i].find_extremes(&mut ext, 0.0);
        let this_below = legends[i].pos.y - ext.3;
        let this_above = legends[i].pos.y - ext.2;
        if i == legends.len() - 1 {
            legends[i].pos.y += this_below;
        } else {
            legends[i].pos.y = legends[i + 1].pos.y + this_below - last_above;
        }
        last_above = this_above;
    }
    Ok(())
}

fn place_top_legend(
    bits: &[String],
    p: &LegendPlacement,
    legend_height: i32,
    text: &mut DrawTextFT,
    legends: &mut Vec<DrawAnnotation>,
) -> Result<()> {
    let gap = calc_line_gap(text, false, p.rel_font_scale);
    let mut total_h = 0.0;
    for (i, bit) in bits.iter().enumerate() {
        total_h += bit_dimensions(bit, p.rel_font_scale, p.colour, text)?.1;
        if i + 1 < bits.len() {
            total_h += gap;
        }
    }
    let band_h = if legend_height > 0 {
        f64::from(legend_height)
    } else {
        p.draw_height
    };
    let mut y = p.base_y + (band_h - total_h) / 2.0;
    for (i, bit) in bits.iter().enumerate() {
        let h = bit_dimensions(bit, p.rel_font_scale, p.colour, text)?.1;
        y += h / 2.0;
        let loc = Point2D::new(p.base_x + p.draw_width / 2.0, y);
        legends.push(legend_note(bit, p.rel_font_scale, loc, p.colour, text)?);
        y += h / 2.0;
        if i + 1 < bits.len() {
            y += gap;
        }
    }
    Ok(())
}

struct SideLegend {
    vert_text: bool,
    is_left: bool,
    legend_width: i32,
    mol_width: i32,
    margin_padding: f64,
}

fn place_side_legend(
    bits: &[String],
    p: &LegendPlacement,
    side: &SideLegend,
    text: &mut DrawTextFT,
    legends: &mut Vec<DrawAnnotation>,
) -> Result<()> {
    if bits.is_empty() {
        return Ok(());
    }
    let mut strip_centre_x = 0.0;
    if side.vert_text {
        strip_centre_x = if side.is_left {
            p.base_x + f64::from(side.legend_width) / 2.0
        } else {
            p.base_x + p.draw_width - f64::from(side.legend_width) / 2.0
        };
    }
    let strip_centre_y = p.base_y + p.draw_height / 2.0;
    let gap = calc_line_gap(text, side.vert_text, p.rel_font_scale);
    if !side.vert_text {
        let mut max_line_w = 0.0f64;
        for bit in bits {
            max_line_w = max_line_w.max(bit_dimensions(bit, p.rel_font_scale, p.colour, text)?.0);
        }
        let pad = 8.0f64.max(side.margin_padding * p.draw_width);
        if side.is_left {
            let mol_panel_left = p.base_x + f64::from(side.legend_width);
            strip_centre_x = mol_panel_left - pad - max_line_w / 2.0;
            strip_centre_x = strip_centre_x.max(p.base_x + max_line_w / 2.0 + 2.0);
        } else {
            let mol_panel_right = p.base_x + f64::from(side.mol_width);
            strip_centre_x = mol_panel_right + pad + max_line_w / 2.0;
            strip_centre_x = strip_centre_x.min(p.base_x + p.draw_width - max_line_w / 2.0 - 2.0);
        }
    }
    // the stack of lines: uniform line height for vertical text, each line's
    // own height otherwise
    let mut heights = Vec::with_capacity(bits.len());
    for bit in bits {
        heights.push(bit_dimensions(bit, p.rel_font_scale, p.colour, text)?.1);
    }
    if side.vert_text {
        let max_h = heights.iter().copied().fold(0.0f64, f64::max);
        heights.iter_mut().for_each(|h| *h = max_h);
    }
    let n = bits.len();
    let total_h: f64 = heights.iter().sum::<f64>() + (n - 1) as f64 * gap;
    let mut y = strip_centre_y - total_h / 2.0;
    for (i, bit) in bits.iter().enumerate() {
        y += heights[i] / 2.0;
        let loc = Point2D::new(strip_centre_x, y);
        legends.push(legend_note(bit, p.rel_font_scale, loc, p.colour, text)?);
        y += heights[i] / 2.0;
        if i + 1 < n {
            y += gap;
        }
    }
    Ok(())
}

impl DrawMol {
    /// If there is a legend, partitions the height or width for it.
    pub fn partition_for_legend(&mut self) {
        self.legend_width = 0;
        self.mol_width = self.draw_width;
        if self.legend.is_empty() {
            self.mol_height = self.draw_height;
            self.legend_height = 0;
            return;
        }
        match self.opts.legend_position {
            LegendPosition::Bottom | LegendPosition::Top => {
                if !self.flexi_canvas_y {
                    self.legend_height =
                        (self.opts.legend_fraction as f32 * self.draw_height as f32) as i32;
                    self.mol_height = self.draw_height - self.legend_height;
                } else {
                    self.mol_height = self.draw_height;
                }
            }
            LegendPosition::Left | LegendPosition::Right => {
                self.legend_width =
                    (self.opts.legend_fraction as f32 * self.draw_width as f32) as i32;
                if self.legend_width <= 0 && self.draw_width > 0 {
                    self.legend_width = 1;
                }
                self.mol_width = self.draw_width - self.legend_width;
                if self.mol_width < 1 {
                    self.mol_width = self.draw_width;
                    self.legend_width = 0;
                }
                self.mol_height = self.draw_height;
                self.legend_height = 0;
            }
        }
    }

    /// Makes the legend annotations. This must be called after
    /// [`DrawMol::calculate_scale`] because the legend font size, given in
    /// pixels in the options, depends on the final font size, and is scaled
    /// down to fit the canvas if necessary.
    pub fn extract_legend(&mut self, text: &mut DrawTextFT) -> Result<()> {
        if self.legend.is_empty() {
            return Ok(());
        }
        let pos = self.opts.legend_position;
        let side = matches!(pos, LegendPosition::Left | LegendPosition::Right);
        let vert_text = side && self.opts.legend_vertical_text;
        let bits = split_legend_bits(&self.legend, vert_text);
        if bits.is_empty() {
            return Ok(());
        }
        let colour = self.opts.legend_colour;
        let fsize = text.font_size();
        let mut rel = f64::from(self.opts.legend_font_size) / fsize;
        let (mut total_width, mut total_height) = calc_legend_dims(&bits, rel, colour, text)?;

        let mut max_width = f64::from(self.draw_width);
        let mut max_height = f64::from(self.draw_height);
        if side {
            if self.legend_width > 0 {
                max_width = f64::from(self.legend_width);
            }
        } else if self.legend_height > 0 {
            max_height = f64::from(self.legend_height);
        }

        // Horizontal text on the sides is not scaled down by the strip width,
        // so its font matches Top/Bottom.
        let side_horizontal = side && !vert_text;
        if total_width > max_width && max_width > 0.0 && !self.flexi_canvas_x && !side_horizontal {
            rel *= max_width / total_width;
            (total_width, total_height) = calc_legend_dims(&bits, rel, colour, text)?;
        }
        if total_height > max_height && max_height > 0.0 {
            if !self.flexi_canvas_y {
                rel *= max_height / total_height;
                (total_width, total_height) = calc_legend_dims(&bits, rel, colour, text)?;
            } else if matches!(pos, LegendPosition::Bottom | LegendPosition::Top) {
                let extra_padding = (total_height * self.margin_padding).max(2.0);
                self.legend_height = (total_height + extra_padding) as i32;
                self.draw_height += self.legend_height;
                self.height += self.legend_height;
                if pos == LegendPosition::Top {
                    self.mol_height = self.draw_height - self.legend_height;
                }
                (total_width, total_height) = calc_legend_dims(&bits, rel, colour, text)?;
            }
        }
        let _ = (total_width, total_height);

        // never shrink horizontal side text below the default legend scale
        let initial = f64::from(self.opts.legend_font_size) / fsize;
        if side_horizontal && rel < initial {
            rel = initial;
        }

        // Vertical side legends use a uniform line height plus gaps: scale so
        // the stack fits the molecule panel.
        if vert_text && self.draw_height > 0 {
            let margin_y = 4.0f64.max(self.margin_padding * f64::from(self.height));
            let avail = 1.0f64.max(f64::from(self.draw_height) - 2.0 * margin_y);
            for _ in 0..8 {
                let mut max_h = 0.0f64;
                for bit in &bits {
                    max_h = max_h.max(bit_dimensions(bit, rel, colour, text)?.1);
                }
                let gap = calc_line_gap(text, vert_text, rel);
                let n = bits.len();
                let span = n as f64 * max_h + if n > 1 { (n - 1) as f64 * gap } else { 0.0 };
                if span <= avail || span <= 0.0 {
                    break;
                }
                rel *= avail / span;
            }
        }

        let placement = LegendPlacement {
            base_x: self.x_offset + f64::from(self.width) * self.margin_padding,
            base_y: self.margin_padding * f64::from(self.height) + self.y_offset,
            draw_width: f64::from(self.draw_width),
            draw_height: f64::from(self.draw_height),
            rel_font_scale: rel,
            colour,
        };
        let mut legends = std::mem::take(&mut self.legends);
        let res = match pos {
            LegendPosition::Bottom => place_bottom_legend(&bits, &placement, text, &mut legends),
            LegendPosition::Top => {
                place_top_legend(&bits, &placement, self.legend_height, text, &mut legends)
            }
            LegendPosition::Left | LegendPosition::Right => {
                let s = SideLegend {
                    vert_text,
                    is_left: pos == LegendPosition::Left,
                    legend_width: self.legend_width,
                    mol_width: self.mol_width,
                    margin_padding: self.margin_padding,
                };
                place_side_legend(&bits, &placement, &s, text, &mut legends)
            }
        };
        self.legends = legends;
        res
    }
}
