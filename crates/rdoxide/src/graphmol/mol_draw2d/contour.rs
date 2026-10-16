//! Contouring and shading of gridded data and sums of Gaussians under a
//! molecule (the contour functions of `MolDraw2DUtils.cpp` of RDKit
//! 2026.09.1).
//!
//! Copyright (C) 2016-2022 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::drawer::MolDraw2DSVG;
use super::helpers::DrawColour;
use crate::geometry::Point2D;
use crate::graphmol::ROMol;
use crate::numerics::conrec::{connect_line_segments, contour};
use crate::{RDKitError, Result};

/// Parameters for drawing contours.
#[derive(Clone, Debug)]
pub struct ContourParams {
    /// assumes the grid is drawn first
    pub set_scale: bool,
    /// use dashed lines for negative contours
    pub dash_negative: bool,
    /// shade the grid
    pub fill_grid: bool,
    /// spacing between elements of the grid
    pub grid_resolution: f64,
    /// line width for drawing contours
    pub contour_width: f64,
    /// extra padding (in molecule coordinates)
    pub extra_grid_padding: f64,
    /// colour for drawing contours
    pub contour_colour: DrawColour,
    /// similarity map colour scheme
    pub colour_map: Vec<DrawColour>,
    /// draws the contours as continuous lines instead of line segments
    pub draw_as_lines: bool,
    /// scaling factor used to convert coordinates to ints when forming the
    /// continuous lines
    pub coord_scale_for_quantization: f64,
    /// scaling factor used to convert isovalues to ints when forming the
    /// continuous lines
    pub isoval_scale_for_quantization: f64,
    /// use a magnitude threshold to determine if a grid box is filled
    pub use_fill_threshold: bool,
    /// magnitude threshold for filling grid boxes
    pub fill_threshold: f64,
    /// if true, the fill threshold is a fraction of the range of the data
    pub fill_threshold_is_fraction: bool,
}

impl Default for ContourParams {
    fn default() -> Self {
        ContourParams {
            set_scale: true,
            dash_negative: true,
            fill_grid: false,
            grid_resolution: 0.15,
            contour_width: 1.0,
            extra_grid_padding: 0.0,
            contour_colour: DrawColour::with_alpha(0.5, 0.5, 0.5, 0.5),
            colour_map: vec![
                DrawColour::with_alpha(0.557, 0.004, 0.322, 0.5),
                DrawColour::with_alpha(1.0, 1.0, 1.0, 0.5),
                DrawColour::with_alpha(0.153, 0.392, 0.098, 0.5),
            ],
            draw_as_lines: true,
            coord_scale_for_quantization: 1000.0,
            isoval_scale_for_quantization: 1e6,
            use_fill_threshold: false,
            fill_threshold: 0.01,
            fill_threshold_is_fraction: true,
        }
    }
}

/// `contourAndDrawGrid`: the data are `grid[i * ycoords.len() + j]`. If
/// `levels` is empty it is filled with `n_contours` levels spanning the data.
pub fn contour_and_draw_grid(
    drawer: &mut MolDraw2DSVG,
    grid: &[f64],
    xcoords: &[f64],
    ycoords: &[f64],
    n_contours: usize,
    levels: &mut Vec<f64>,
    params: &ContourParams,
    mol: Option<&ROMol>,
) -> Result<()> {
    if grid.is_empty() {
        return Err(RDKitError::precondition("no data"));
    }
    if params.colour_map.len() <= 1 {
        return Err(RDKitError::precondition(
            "colourMap must have at least two entries",
        ));
    }
    let (n_x, n_y) = (xcoords.len(), ycoords.len());
    if grid.len() < n_x * n_y || n_x == 0 || n_y == 0 {
        return Err(RDKitError::precondition(
            "grid and coordinates do not match",
        ));
    }
    if params.set_scale {
        let min_p = Point2D::new(xcoords[0], ycoords[0]);
        let max_p = Point2D::new(xcoords[n_x - 1], ycoords[n_y - 1]);
        drawer.set_scale_for_range(drawer.width(), drawer.height(), min_p, max_p, mol)?;
    }

    let mut min_v = f64::MAX;
    let mut max_v = f64::MIN;
    if levels.is_empty() || params.fill_grid {
        for i in 0..n_x {
            for j in 0..n_y {
                min_v = min_v.min(grid[i * n_y + j]);
                max_v = max_v.max(grid[i * n_y + j]);
            }
        }
        if levels.is_empty() {
            levels.extend(
                (0..n_contours)
                    .map(|i| min_v + i as f64 * (max_v - min_v) / (n_contours as f64 - 1.0)),
            );
        }
    }
    if max_v <= min_v {
        return Ok(());
    }

    let olw = drawer.line_width();
    let odash = drawer.dash().to_vec();
    let ocolour = drawer.colour();
    let ofill = drawer.fill_polys();
    if params.fill_grid {
        drawer.set_fill_polys(true);
        drawer.set_line_width(1.0);
        let mut delta = max_v - min_v;
        if params.colour_map.len() > 2 {
            // need to find how fractionally far we are from zero, not the min
            delta = if -min_v > max_v { -min_v } else { max_v };
        }
        for i in 0..n_x - 1 {
            for j in 0..n_y - 1 {
                let grid_v = grid[i * n_y + j];
                let mut thresh_test = grid_v;
                if params.fill_threshold_is_fraction {
                    thresh_test /= delta;
                }
                if params.use_fill_threshold && thresh_test.abs() < params.fill_threshold {
                    continue;
                }
                let mut frac_v = (grid_v - min_v) / delta;
                if params.colour_map.len() > 2 {
                    // need to find how fractionally far we are from zero
                    frac_v = grid_v / delta;
                    if frac_v < 0.0 {
                        frac_v *= -1.0;
                    }
                }
                let negative_or_two = grid_v < 0.0 || params.colour_map.len() == 2;
                let c1 = params.colour_map[1];
                let c2 = if negative_or_two {
                    params.colour_map[0]
                } else {
                    params.colour_map[2]
                };
                let c = c1 + (c2 - c1) * frac_v;
                // don't bother drawing boxes that are the same as the
                // background colour
                if c.feq(&drawer.draw_options().background_colour, 0.01, true) {
                    continue;
                }
                drawer.set_colour(c);
                let p1 = Point2D::new(xcoords[i], ycoords[j]);
                let p2 = Point2D::new(xcoords[i + 1], ycoords[j + 1]);
                drawer.draw_rect(p1, p2, false)?;
            }
        }
    }

    if n_contours > 0 {
        if n_contours > levels.len() {
            drawer_restore(drawer, olw, &odash, ocolour, ofill);
            return Err(RDKitError::Value(
                "nContours larger than the size of the level list".into(),
            ));
        }
        let segs = contour(
            grid,
            0,
            n_x - 1,
            0,
            n_y - 1,
            xcoords,
            ycoords,
            n_contours,
            levels,
        )?;
        let neg_dash = vec![2.0, 2.0];
        drawer.set_colour(params.contour_colour);
        drawer.set_line_width(params.contour_width);
        if !params.draw_as_lines {
            for seg in &segs {
                if params.dash_negative && seg.iso_val < 0.0 {
                    drawer.set_dash(neg_dash.clone());
                } else {
                    drawer.set_dash(Vec::new());
                }
                drawer.draw_line(seg.p1, seg.p2, false)?;
            }
        } else {
            drawer.set_fill_polys(false);
            let lines = connect_line_segments(
                &segs,
                params.coord_scale_for_quantization,
                params.isoval_scale_for_quantization,
            );
            for (pts, val) in lines {
                if params.dash_negative && val < 0.0 {
                    drawer.set_dash(neg_dash.clone());
                } else {
                    drawer.set_dash(Vec::new());
                }
                if pts.len() > 2 {
                    drawer.draw_polygon(&pts, false)?;
                } else if pts.len() == 2 {
                    drawer.draw_line(pts[0], pts[1], false)?;
                }
            }
        }
        drawer_restore(drawer, olw, &odash, ocolour, ofill);
    }
    Ok(())
}

fn drawer_restore(
    drawer: &mut MolDraw2DSVG,
    line_width: f64,
    dash: &[f64],
    colour: DrawColour,
    fill: bool,
) {
    drawer.set_dash(dash.to_vec());
    drawer.set_line_width(line_width);
    drawer.set_colour(colour);
    drawer.set_fill_polys(fill);
}

/// `contourAndDrawGaussians`: contours the sum of Gaussians at `locs` with
/// the given heights (weights) and widths (standard deviations).
pub fn contour_and_draw_gaussians(
    drawer: &mut MolDraw2DSVG,
    locs: &[Point2D],
    weights: &[f64],
    widths: &[f64],
    n_contours: usize,
    levels: &mut Vec<f64>,
    params: &ContourParams,
    mol: Option<&ROMol>,
) -> Result<()> {
    if locs.len() != weights.len() || locs.len() != widths.len() {
        return Err(RDKitError::precondition("size mismatch"));
    }
    // start by setting up the grid
    if params.set_scale {
        let mut min_p = Point2D::new(f64::MAX, f64::MAX);
        let mut max_p = Point2D::new(f64::MIN, f64::MIN);
        for loc in locs {
            min_p.x = min_p.x.min(loc.x);
            min_p.y = min_p.y.min(loc.y);
            max_p.x = max_p.x.max(loc.x);
            max_p.y = max_p.y.max(loc.y);
        }
        let dims = max_p - min_p;
        // Here, the padding is just used to extend the grid beyond the
        // molecule. The actual padding round the image is added later.
        let padding = drawer.draw_options().padding;
        min_p.x -= padding * dims.x;
        min_p.y -= padding * dims.y;
        max_p.x += padding * dims.x;
        max_p.y += padding * dims.y;
        if params.extra_grid_padding > 0.0 {
            min_p.x -= params.extra_grid_padding;
            min_p.y -= params.extra_grid_padding;
            max_p.x += params.extra_grid_padding;
            max_p.y += params.extra_grid_padding;
        }
        drawer.set_scale_for_range(drawer.width(), drawer.height(), min_p, max_p, mol)?;
    }

    let range = drawer.range()?;
    let min_pt = drawer.min_pt()?;
    let nx = (range.x / params.grid_resolution).ceil() as usize + 1;
    let ny = (range.y / params.grid_resolution).ceil() as usize + 1;
    let xcoords: Vec<f64> = (0..nx)
        .map(|i| min_pt.x + i as f64 * params.grid_resolution)
        .collect();
    let ycoords: Vec<f64> = (0..ny)
        .map(|i| min_pt.y + i as f64 * params.grid_resolution)
        .collect();
    let mut grid = vec![0.0f64; nx * ny];
    // populate the grid from the gaussians
    for (ix, &px) in xcoords.iter().enumerate() {
        for (iy, &py) in ycoords.iter().enumerate() {
            let pt = Point2D::new(px, py);
            let mut accum = 0.0;
            for ig in 0..locs.len() {
                let d2 = (pt - locs[ig]).length_sq();
                accum += weights[ig] / widths[ig] * (-0.5 * d2 / (widths[ig] * widths[ig])).exp();
            }
            grid[ix * ny + iy] = accum / (2.0 * PI);
        }
    }

    // and render it; if scaling was needed, it was done already
    let mut params_copy = params.clone();
    params_copy.set_scale = false;
    contour_and_draw_grid(
        drawer,
        &grid,
        &xcoords,
        &ycoords,
        n_contours,
        levels,
        &params_copy,
        None,
    )
}
