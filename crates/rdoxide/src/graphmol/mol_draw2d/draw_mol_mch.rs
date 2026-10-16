//! Multi-colour highlights drawn as circles (or arcs) round atoms and lines
//! beside bonds (`DrawMolMCH.cpp` and `DrawMolMCHCircleAndLine.cpp` of RDKit
//! 2026.09.1).
//!
//! Copyright (C) 2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeSet;

use super::details::{adjust_line_end_for_ellipse, does_line_intersect_arc};
use super::draw_annotation::EMPTY_EXTREMES;
use super::draw_mol::{DrawMol, McHighlights};
use super::draw_mol_utils::{calc_perpendicular, get_highlight_bond_width};
use super::draw_shape::{DrawShape, no_dash};
use crate::geometry::Point2D;
use crate::{RDKitError, Result};

impl DrawMol {
    /// The radii of the highlight of an atom: the option, unless the atom has
    /// a radius of its own.
    pub fn mch_atom_radius(&self, atom_idx: u32) -> (f64, f64) {
        let r = self
            .highlight_radii
            .get(&(atom_idx as i32))
            .copied()
            .unwrap_or(self.opts.highlight_radius);
        (r, r)
    }

    /// `extractMCHighlights` of `DrawMolMCHCircleAndLine`.
    pub fn extract_mc_circle_and_line(&mut self, mch: &McHighlights) -> Result<()> {
        let mut bond_highlights = self.make_mc_bond_highlights(mch)?;
        let mut atom_highlights = self.make_mc_atom_highlights(mch)?;
        Self::fix_highlight_join_problems(&mut atom_highlights, &mut bond_highlights);
        self.highlights.extend(bond_highlights);
        self.highlights.extend(atom_highlights);
        Ok(())
    }

    fn make_mc_bond_highlights(&self, mch: &McHighlights) -> Result<Vec<DrawShape>> {
        let mut res = Vec::new();
        for (&bond_idx, colours) in &mch.bond_map {
            if bond_idx < 0 || bond_idx as u32 >= self.draw_mol.num_bonds() {
                return Err(RDKitError::Value("Bond index out of range".into()));
            }
            let mut line_width = self.opts.bond_line_width;
            if !self.opts.fill_highlights {
                line_width = get_highlight_bond_width(
                    &self.opts,
                    bond_idx,
                    Some(&mch.linewidth_multipliers),
                );
            }
            let bond = self.draw_mol.bond(bond_idx as u32);
            let at1_idx = bond.begin_atom_idx();
            let at2_idx = bond.end_atom_idx();
            let at1_cds = self.at_cds[at1_idx as usize];
            let at2_cds = self.at_cds[at2_idx as usize];
            let perp = calc_perpendicular(at1_cds, at2_cds);
            let rad = 0.7 * self.opts.highlight_radius;
            let ids = (at1_idx as i32, at2_idx as i32, bond_idx);
            let scale_lw = self.opts.scale_bond_width;

            let adjusted_line =
                |res: &mut Vec<DrawShape>, mut p1: Point2D, mut p2: Point2D, col| -> Result<()> {
                    self.adjust_line_end_for_highlight(at1_idx, p2, &mut p1);
                    self.adjust_line_end_for_highlight(at2_idx, p1, &mut p2);
                    res.push(DrawShape::simple_line(
                        vec![p1, p2],
                        line_width,
                        scale_lw,
                        col,
                        ids,
                        no_dash(),
                    )?);
                    Ok(())
                };

            if colours.len() < 2 {
                let col = colours
                    .first()
                    .copied()
                    .unwrap_or(self.opts.highlight_colour);
                if self.opts.fill_highlights {
                    let pts = vec![
                        at1_cds + perp * rad,
                        at2_cds + perp * rad,
                        at2_cds - perp * rad,
                        at1_cds - perp * rad,
                    ];
                    res.push(DrawShape::poly_line(
                        pts,
                        line_width,
                        scale_lw,
                        col,
                        true,
                        ids,
                        no_dash(),
                    )?);
                } else {
                    adjusted_line(&mut res, at1_cds + perp * rad, at2_cds + perp * rad, col)?;
                    adjusted_line(&mut res, at1_cds - perp * rad, at2_cds - perp * rad, col)?;
                }
            } else {
                let col_rad = 2.0 * rad / colours.len() as f64;
                if self.opts.fill_highlights {
                    let mut p1 = at1_cds - perp * rad;
                    let mut p2 = at2_cds - perp * rad;
                    for col in colours {
                        let pts = vec![p1, p1 + perp * col_rad, p2 + perp * col_rad, p2];
                        res.push(DrawShape::poly_line(
                            pts,
                            line_width,
                            scale_lw,
                            *col,
                            true,
                            ids,
                            no_dash(),
                        )?);
                        p1 += perp * col_rad;
                        p2 += perp * col_rad;
                    }
                } else {
                    let mut cols: &[_] = colours;
                    if cols.len() % 2 == 1 {
                        adjusted_line(&mut res, at1_cds, at2_cds, cols[0])?;
                        cols = &cols[1..];
                    }
                    let mut step = 0;
                    for (i, col) in cols.iter().enumerate() {
                        // draw even numbers from the bottom, odd from the top
                        let offset = perp * (rad - f64::from(step) * col_rad);
                        if i % 2 == 0 {
                            adjusted_line(&mut res, at1_cds - offset, at2_cds - offset, *col)?;
                        } else {
                            adjusted_line(&mut res, at1_cds + offset, at2_cds + offset, *col)?;
                            step += 1;
                        }
                    }
                }
            }
        }
        Ok(res)
    }

    fn make_mc_atom_highlights(&self, mch: &McHighlights) -> Result<Vec<DrawShape>> {
        let mut res = Vec::new();
        for (&atom_idx, colours) in &mch.atom_map {
            if atom_idx < 0 || atom_idx as u32 >= self.draw_mol.num_atoms() {
                return Err(RDKitError::Value("Atom index out of range".into()));
            }
            // the width is an int upstream
            let line_width = get_highlight_bond_width(&self.opts, -1, None).trunc();
            let (centre, xradius, yradius) = self.calc_symbol_ellipse(atom_idx as u32)?;
            if colours.len() == 1 {
                res.push(DrawShape::ellipse(
                    vec![centre, Point2D::new(xradius, yradius)],
                    line_width,
                    true,
                    colours[0],
                    self.opts.fill_highlights,
                    atom_idx,
                )?);
            } else {
                let arc_size = 360.0 / colours.len() as f64;
                let mut arc_start = 270.0;
                for col in colours {
                    let mut arc_stop = arc_start + arc_size;
                    if arc_stop >= 360.0 {
                        arc_stop -= 360.0;
                    }
                    res.push(DrawShape::arc(
                        vec![centre, Point2D::new(xradius, yradius)],
                        arc_start,
                        arc_stop,
                        line_width,
                        true,
                        *col,
                        self.opts.fill_highlights,
                        atom_idx,
                    )?);
                    arc_start += arc_size;
                    if arc_start >= 360.0 {
                        arc_start -= 360.0;
                    }
                }
            }
        }
        Ok(res)
    }

    /// Moves `p2` onto the ellipse round the atom's symbol where the line
    /// from `p1` crosses it.
    fn adjust_line_end_for_highlight(&self, at_idx: u32, p1: Point2D, p2: &mut Point2D) {
        let Ok((centre, xradius, yradius)) = self.calc_symbol_ellipse(at_idx) else {
            return;
        };
        if xradius < 1.0e-6 || yradius < 1.0e-6 {
            return;
        }
        adjust_line_end_for_ellipse(centre, xradius, yradius, p1, p2);
    }

    /// The ellipse round an atom's symbol: its centre and radii.
    fn calc_symbol_ellipse(&self, atom_idx: u32) -> Result<(Point2D, f64, f64)> {
        let idx = atom_idx as usize;
        if idx >= self.at_cds.len() || idx >= self.atom_labels.len() {
            return Err(RDKitError::precondition("atom index out of range"));
        }
        let mut centre = self.at_cds[idx];
        let (mut xradius, mut yradius) = self.mch_atom_radius(atom_idx);
        let label = match &self.atom_labels[idx] {
            Some(l) if !self.opts.atom_highlights_are_circles && !l.symbol.is_empty() => l,
            _ => return Ok((centre, xradius, yradius)),
        };
        let mut ext = EMPTY_EXTREMES;
        label.find_extremes(&mut ext);
        let (x_min, x_max, y_min, y_max) = ext;
        let root_2 = 2.0f64.sqrt();
        xradius = xradius.max(root_2 * 0.5 * (x_max - x_min));
        yradius = yradius.max(root_2 * 0.5 * (y_max - y_min));
        centre.x = 0.5 * (x_max + x_min);
        centre.y = 0.5 * (y_max + y_min);
        Ok((centre, xradius, yradius))
    }

    /// Sometimes the bond highlights went outside the atom highlight ellipse
    /// they were supposed to be intersecting, particularly when the atom
    /// highlights were circles rather than ellipses. This detects these cases
    /// and makes the atom highlight radius larger.
    fn fix_highlight_join_problems(
        atom_highlights: &mut [DrawShape],
        bond_highlights: &mut [DrawShape],
    ) {
        let mut fettled: BTreeSet<usize> = BTreeSet::new();
        for (i, atom_hl) in atom_highlights.iter_mut().enumerate() {
            for bond_hl in bond_highlights.iter() {
                if atom_hl.atom1 != bond_hl.atom1 && atom_hl.atom1 != bond_hl.atom2 {
                    continue;
                }
                // A multicoloured bond highlight is a polyline with 4 points,
                // so 2 is opposite 0. Normally a bond highlight is a single
                // line, with only 2 points.
                let p2 = if bond_hl.points.len() == 2 { 1 } else { 2 };
                let intersects = |a: &DrawShape| {
                    does_line_intersect_arc(
                        a.points[0],
                        a.points[1].x,
                        a.points[1].y,
                        0.0,
                        360.0,
                        0.0,
                        bond_hl.points[0],
                        bond_hl.points[p2],
                    )
                };
                let mut ins = intersects(atom_hl);
                if ins {
                    continue;
                }
                fettled.insert(i);
                while !ins {
                    let dist1 = (atom_hl.points[0] - bond_hl.points[0]).length();
                    let mut dist2 = (atom_hl.points[0] - bond_hl.points[p2]).length();
                    let maxrad = atom_hl.points[1].x.max(atom_hl.points[1].y);
                    // If it's a filled highlight, the bond highlight won't
                    // have been truncated at the outer radius of the circle,
                    // because there's no need. Set dist2 to the radius now.
                    if dist2 < 1.0e-4 {
                        dist2 = maxrad;
                    }
                    // This scales the radii of the ellipse by 1.25 times the
                    // ratio of the shorter distance of the bond highlight ends
                    // from the atom and the larger radius of the ellipse, to
                    // make sure the ellipse is wide enough that the bond
                    // highlight intersects it.
                    let mut upscaler = 1.25 * dist1.min(dist2) / maxrad;
                    // We certainly don't want to make the ellipse smaller.
                    if upscaler <= 1.0 {
                        upscaler = 1.1;
                    }
                    atom_hl.points[1] *= upscaler;
                    ins = intersects(atom_hl);
                }
            }
        }
        for fa in fettled {
            // now adjust all the other bond highlights that end on this atom
            let atom_hl = &atom_highlights[fa];
            for bond_hl in bond_highlights.iter_mut() {
                if atom_hl.atom1 != bond_hl.atom1 && atom_hl.atom1 != bond_hl.atom2 {
                    continue;
                }
                let (a0, r) = (atom_hl.points[0], atom_hl.points[1]);
                if (a0 - bond_hl.points[0]).length_sq() < (a0 - bond_hl.points[1]).length_sq() {
                    let other = bond_hl.points[1];
                    adjust_line_end_for_ellipse(a0, r.x, r.y, other, &mut bond_hl.points[0]);
                } else {
                    let other = bond_hl.points[0];
                    adjust_line_end_for_ellipse(a0, r.x, r.y, other, &mut bond_hl.points[1]);
                }
            }
        }
    }
}
