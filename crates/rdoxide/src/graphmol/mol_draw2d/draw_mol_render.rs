//! Building the drawing objects of a `DrawMol` from start to finish, and
//! rendering them on a canvas (`createDrawObjects`, `setScale`,
//! `setTransformation`, `draw` and `drawRadicals` of `DrawMol.cpp`).
//!
//! Copyright (C) 2021-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::draw_mol::DrawMol;
use super::draw_shape::ShapeCanvas;
use super::draw_text::DrawTextFT;
use super::helpers::{DrawColour, OrientType};
use crate::Result;
use crate::geometry::Point2D;

/// What a `DrawMol` needs from a drawer to render itself: the shape drawing
/// operations plus the scale in force and a place for text, which comes
/// ready-made as SVG.
pub trait MolCanvas: ShapeCanvas {
    /// The scale of the molecule being drawn, used for the line widths.
    fn scale(&self) -> f64;
    fn set_scale(&mut self, scale: f64);
    /// The class the drawer has been asked to put on everything.
    fn active_class(&self) -> String;
    /// Adds SVG text, as made by the text drawer, to the output.
    fn write_text(&mut self, svg: &str);
}

impl DrawMol {
    /// Builds everything to be drawn, to the size and scale the options and
    /// canvas call for, ending with all in drawing coordinates.
    pub fn create_draw_objects(&mut self, text: &mut DrawTextFT) -> Result<()> {
        text.set_font_scale(self.font_scale, true);
        self.partition_for_legend();
        self.extract_all(self.scale, text)?;
        self.calculate_scale(text);

        let ignore_font_limits = self.opts.fixed_font_size != -1;
        if !text.set_font_scale(self.font_scale, ignore_font_limits) || ignore_font_limits {
            // The relative font size isn't what we were expecting, so
            // everything must be rebuilt. Furthermore, on a fully flexible
            // canvas where the font scale is greater than the global scale,
            // characters at the edge of the image would not fit (Github6111),
            // so rebuild with an appropriate relative font size.
            if self.flexi_canvas_x && self.flexi_canvas_y && (self.font_scale - self.scale) > 1e-4 {
                self.width = -1;
                self.height = -1;
                let curr_scale = text.font_scale();
                let rel_scale = self.font_scale / self.scale;
                self.reset_everything(text);
                self.font_scale = rel_scale;
                text.set_font_scale(rel_scale, true);
                self.extract_all(self.scale, text)?;
                self.calculate_scale(text);
                text.set_font_scale(curr_scale, true);
            }
            self.set_scale(self.scale, text.font_scale(), ignore_font_limits, text)
        } else {
            self.finish_create_draw_objects(text)
        }
    }

    /// Rebuilds everything at a given scale and font scale.
    pub fn set_scale(
        &mut self,
        new_scale: f64,
        new_font_scale: f64,
        ignore_font_limits: bool,
        text: &mut DrawTextFT,
    ) -> Result<()> {
        self.reset_everything(text);
        self.font_scale = new_font_scale / new_scale;
        text.set_font_scale(self.font_scale, true);

        self.extract_all(new_scale, text)?;
        self.find_extremes();

        text.set_font_scale(new_font_scale, ignore_font_limits);
        self.scale = new_scale;
        self.font_scale = text.font_scale();
        self.finish_create_draw_objects(text)
    }

    /// Rebuilds everything with the scale and extremes of another, so that
    /// molecules in the same picture share one transformation.
    pub fn set_transformation(&mut self, source: &DrawMol, text: &mut DrawTextFT) -> Result<()> {
        self.reset_everything(text);
        let rel_font_scale = source.font_scale / source.scale;
        text.set_font_scale(rel_font_scale, true);
        self.x_min = source.x_min;
        self.x_max = source.x_max;
        self.y_min = source.y_min;
        self.y_max = source.y_max;
        self.x_range = source.x_range;
        self.y_range = source.y_range;

        self.extract_all(self.scale, text)?;
        self.scale = source.scale;
        self.font_scale = source.font_scale;
        text.set_font_scale(self.font_scale, true);
        self.finish_create_draw_objects(text)
    }

    /// Renders everything, in the order of the layers.
    pub fn draw(&self, canvas: &mut dyn MolCanvas, text: &mut DrawTextFT) -> Result<()> {
        if !self.drawing_initialised {
            return Err(crate::RDKitError::precondition(
                "you must call createDrawingObjects before calling draw",
            ));
        }
        if self.at_cds.is_empty() {
            return Ok(());
        }
        let keep_scale = canvas.scale();
        canvas.set_scale(self.scale);
        let keep_font_scale = text.font_scale();
        text.set_font_scale(self.font_scale, true);

        let res = self.draw_layers(canvas, text);

        canvas.set_scale(keep_scale);
        text.set_font_scale(keep_font_scale, true);
        res
    }

    fn draw_layers(&self, canvas: &mut dyn MolCanvas, text: &mut DrawTextFT) -> Result<()> {
        for shape in &self.pre_shapes {
            shape.draw(canvas)?;
        }
        for shape in &self.highlights {
            shape.draw(canvas)?;
        }
        for shape in &self.bonds {
            shape.draw(canvas)?;
        }
        for label in self.atom_labels.iter().flatten() {
            let svg = label.draw_svg(text, &canvas.active_class())?;
            canvas.write_text(&svg);
        }
        if self.include_annotations {
            for annot in &self.annotations {
                let svg = annot.draw_svg(text, &canvas.active_class())?;
                canvas.write_text(&svg);
            }
        }
        if self.opts.include_radicals {
            self.draw_radicals(canvas)?;
        }
        for shape in &self.post_shapes {
            shape.draw(canvas)?;
        }
        for legend in &self.legends {
            let svg = legend.draw_svg(text, &canvas.active_class())?;
            canvas.write_text(&svg);
        }
        Ok(())
    }

    /// The radical spots, next to their atoms.
    fn draw_radicals(&self, canvas: &mut dyn MolCanvas) -> Result<()> {
        // take account of differing font scale and main scale if we've hit
        // max or min font size.
        let spot_rad = 0.2 * self.opts.multiple_bond_offset * self.font_scale;
        canvas.set_colour(DrawColour::new(0.0, 0.0, 0.0));
        let draw_spot = |canvas: &mut dyn MolCanvas, cds: Point2D| -> Result<()> {
            let ofp = canvas.fill_polys();
            canvas.set_fill_polys(true);
            let olw = canvas.line_width();
            canvas.set_line_width(0.0);
            let res = canvas.draw_arc(cds, spot_rad, spot_rad, 0.0, 360.0, true);
            canvas.set_line_width(olw);
            canvas.set_fill_polys(ofp);
            res
        };
        // offsets the spot along x (dir 0) or y
        let shifted = |cds: Point2D, by: f64, dir: bool| {
            if dir {
                Point2D::new(cds.x, cds.y + by)
            } else {
                Point2D::new(cds.x + by, cds.y)
            }
        };
        let draw_spots = |canvas: &mut dyn MolCanvas,
                          cds: Point2D,
                          num_spots: u32,
                          width: f64,
                          dir: bool|
         -> Result<()> {
            match num_spots {
                1 | 3 => {
                    if num_spots == 3 {
                        draw_spot(canvas, shifted(cds, -0.6 * width + spot_rad, dir))?;
                        draw_spot(canvas, shifted(cds, 0.6 * width - spot_rad, dir))?;
                    }
                    draw_spot(canvas, cds)
                }
                2 | 4 => {
                    if num_spots == 4 {
                        draw_spot(canvas, shifted(cds, 6.0 * spot_rad, dir))?;
                        draw_spot(canvas, shifted(cds, -6.0 * spot_rad, dir))?;
                    }
                    draw_spot(canvas, shifted(cds, 2.0 * spot_rad, dir))?;
                    draw_spot(canvas, shifted(cds, -2.0 * spot_rad, dir))
                }
                _ => Ok(()),
            }
        };

        for (rect, orient, at_idx) in &self.radicals {
            let num_rade = self.draw_mol.atom(*at_idx as u32).num_radical_electrons();
            canvas.set_active_atm_idx(*at_idx, -1);
            let res = if matches!(orient, OrientType::N | OrientType::S | OrientType::C) {
                draw_spots(canvas, rect.trans, num_rade, rect.width, false)
            } else {
                draw_spots(canvas, rect.trans, num_rade, rect.height, true)
            };
            canvas.set_active_atm_idx(-1, -1);
            res?;
        }
        Ok(())
    }
}
