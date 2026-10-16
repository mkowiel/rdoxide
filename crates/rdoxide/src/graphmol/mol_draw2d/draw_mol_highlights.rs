//! Highlights, regions and attachment marks of a [`DrawMol`] (the matching
//! methods of `DrawMol.cpp`, RDKit 2026.09.1).
//!
//! Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::details::do_lines_intersect;
use super::draw_annotation::EMPTY_EXTREMES;
use super::draw_mol::DrawMol;
use super::draw_mol_utils::{
    calc_inner_perpendicular, calc_perpendicular, get_highlight_bond_colour,
    get_highlight_bond_width,
};
use super::draw_shape::{DrawShape, no_dash};
use super::helpers::{DrawColour, LINE_WIDTH_SCALE_FACTOR, MultiColourHighlightStyle};
use crate::Result;
use crate::geometry::Point2D;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::{RDProps, common_properties as cp};

fn mul(p: Point2D, f: f64) -> Point2D {
    Point2D::new(p.x * f, p.y * f)
}

impl DrawMol {
    /// Highlights drawn behind the molecule: fat bond lines and atom ellipses
    /// for continuous highlights, circles otherwise.
    pub fn extract_highlights(&mut self, scale: f64) -> Result<()> {
        // the lasso style makes all the highlights itself
        if let Some(mch) = self.mch.clone() {
            if mch.style == MultiColourHighlightStyle::Lasso {
                return self.extract_mc_lasso(&mch);
            }
        }
        if self.opts.continuous_highlight {
            self.make_continuous_highlights(scale)?;
        } else if self.opts.circle_atoms && !self.highlight_atoms.is_empty() {
            self.make_atom_circle_highlights()?;
        }
        if let Some(mch) = self.mch.clone() {
            self.extract_mc_circle_and_line(&mch)?;
        }
        Ok(())
    }

    /// Grey boxes behind the regions of the options.
    pub fn extract_regions(&mut self) -> Result<()> {
        for region in &self.opts.atom_regions {
            if region.len() <= 1 {
                continue;
            }
            let first = self.at_cds[region[0] as usize];
            let (mut minv, mut maxv) = (first, first);
            for &idx in region {
                let pt = self.at_cds[idx as usize];
                minv.x = minv.x.min(pt.x);
                minv.y = minv.y.min(pt.y);
                maxv.x = maxv.x.max(pt.x);
                maxv.y = maxv.y.max(pt.y);
            }
            let size = mul(maxv - minv, 0.2);
            minv -= Point2D::new(size.x / 2.0, size.y / 2.0);
            maxv += Point2D::new(size.x / 2.0, size.y / 2.0);
            let pts = vec![
                minv,
                Point2D::new(minv.x, maxv.y),
                maxv,
                Point2D::new(maxv.x, minv.y),
            ];
            self.highlights.push(DrawShape::poly_line(
                pts,
                1.0,
                false,
                DrawColour::new(0.8, 0.8, 0.8),
                true,
                (-1, -1, -1),
                no_dash(),
            )?);
        }
        Ok(())
    }

    /// With `dummies_are_attachments`, a wavy line across the bond of each
    /// terminal dummy atom.
    pub fn extract_attachments(&mut self) -> Result<()> {
        if !self.opts.dummies_are_attachments {
            return Ok(());
        }
        for idx in 0..self.draw_mol.num_atoms() {
            let atom = self.draw_mol.atom(idx);
            if atom.has_prop(cp::atomLabel) || self.opts.atom_labels.contains_key(&(idx as i32)) {
                // dummies with an explicit label are skipped
                continue;
            }
            if atom.atomic_num() != 0 || self.draw_mol.atom_degree(idx) != 1 {
                continue;
            }
            let at2 = self
                .draw_mol
                .atom_neighbors(idx)
                .next()
                .expect("a degree 1 atom has a neighbour");
            let at1_cds = self.at_cds[idx as usize];
            let at2_cds = self.at_cds[at2 as usize];
            let perp = calc_perpendicular(at1_cds, at2_cds);
            let p1 = Point2D::new(at1_cds.x - perp.x * 0.5, at1_cds.y - perp.y * 0.5);
            let p2 = Point2D::new(at1_cds.x + perp.x * 0.5, at1_cds.y + perp.y * 0.5);
            let col = DrawColour::new(0.5, 0.5, 0.5);
            let offset = self.opts.multiple_bond_offset * self.mean_bond_length / 2.0;
            self.bonds.push(DrawShape::wavy_line(
                vec![p1, p2],
                self.opts.bond_line_width,
                false,
                col,
                col,
                offset,
                (at2 as i32 + self.active_atm_idx_offset, -1, -1),
            )?);
        }
        Ok(())
    }

    fn make_continuous_highlights(&mut self, scale: f64) -> Result<()> {
        let mut tgt_lw = get_highlight_bond_width(&self.opts, -1, None);
        if tgt_lw < 2.0 {
            tgt_lw = 2.0;
        }
        if !self.opts.continuous_highlight {
            tgt_lw /= 4.0;
        }
        if !self.highlight_bonds.is_empty() {
            self.make_bond_highlight_lines(tgt_lw, scale)?;
        }
        if !self.highlight_atoms.is_empty() {
            self.make_atom_ellipse_highlights(tgt_lw)?;
        }
        Ok(())
    }

    fn make_atom_circle_highlights(&mut self) -> Result<()> {
        for idx in 0..self.draw_mol.num_atoms() {
            let key = idx as i32;
            if !self.highlight_atoms.contains(&key) {
                continue;
            }
            let col = self
                .highlight_atom_map
                .get(&key)
                .copied()
                .unwrap_or(self.opts.highlight_colour);
            let radius = self
                .highlight_radii
                .get(&key)
                .copied()
                .unwrap_or(self.opts.highlight_radius);
            self.highlights.push(DrawShape::ellipse(
                vec![self.at_cds[idx as usize], Point2D::new(radius, radius)],
                2.0,
                false,
                col,
                true,
                key + self.active_atm_idx_offset,
            )?);
        }
        Ok(())
    }

    fn make_atom_ellipse_highlights(&mut self, mut line_width: f64) -> Result<()> {
        if !self.opts.fill_highlights {
            // a narrower circle is needed
            line_width /= 2.0;
        }
        for idx in 0..self.draw_mol.num_atoms() {
            let key = idx as i32;
            if !self.highlight_atoms.contains(&key) {
                continue;
            }
            let col = self
                .highlight_atom_map
                .get(&key)
                .copied()
                .unwrap_or(self.opts.highlight_colour);
            let mut centre = self.at_cds[idx as usize];
            let mut xradius = self
                .highlight_radii
                .get(&key)
                .copied()
                .unwrap_or(self.opts.highlight_radius);
            let mut yradius = xradius;
            if !self.opts.atom_highlights_are_circles {
                if let Some(label) = &self.atom_labels[idx as usize] {
                    let mut ext = EMPTY_EXTREMES;
                    label.find_extremes(&mut ext);
                    let root_2 = 2.0f64.sqrt();
                    xradius = xradius.max(root_2 * 0.5 * (ext.1 - ext.0));
                    yradius = yradius.max(root_2 * 0.5 * (ext.3 - ext.2));
                    centre.x = 0.5 * (ext.1 + ext.0);
                    centre.y = 0.5 * (ext.3 + ext.2);
                }
            }
            self.highlights.push(DrawShape::ellipse(
                vec![centre, Point2D::new(xradius, yradius)],
                line_width,
                true,
                col,
                true,
                key + self.active_atm_idx_offset,
            )?);
        }
        Ok(())
    }

    /// The neighbours of `atom` other than `other_atom` that are bonded to it
    /// by a highlighted bond.
    fn find_high_bond_nbrs(&self, atom: u32, other_atom: u32) -> Vec<u32> {
        self.draw_mol
            .atom_bonds(atom)
            .filter_map(|b| {
                let nbr = b.other_atom_idx(atom).ok()?;
                (nbr != other_atom && self.highlight_bonds.contains(&(b.idx() as i32)))
                    .then_some(nbr)
            })
            .collect()
    }

    /// Highlighted bonds are 4-6 sided polygons under the bond: a rectangle
    /// for an isolated highlight, mitred where it joins another highlighted
    /// bond so there are no gaps.
    fn make_bond_highlight_lines(&mut self, mut line_width: f64, scale: f64) -> Result<()> {
        // this is essentially the inverse of MolDraw2D::getDrawLineWidth,
        // which ignores scaleHighlightBondWidth and only uses scaleBondWidth
        if !self.opts.scale_highlight_bond_width {
            // so that scaling it up again gives the right size
            line_width /= scale;
        } else {
            // same conversion factor as in MolDraw2D::getDrawLineWidth
            line_width *= LINE_WIDTH_SCALE_FACTOR;
        }
        for this_idx in 0..self.draw_mol.num_atoms() {
            let bonds: Vec<(u32, u32)> = self
                .draw_mol
                .atom_bonds(this_idx)
                .filter_map(|b| Some((b.idx(), b.other_atom_idx(this_idx).ok()?)))
                .collect();
            for (bond_idx, nbr_idx) in bonds {
                if !((nbr_idx as usize) < self.at_cds.len() && nbr_idx > this_idx) {
                    continue;
                }
                if !self.highlight_bonds.contains(&(bond_idx as i32)) {
                    continue;
                }
                let col = get_highlight_bond_colour(
                    &self.draw_mol,
                    bond_idx,
                    &self.opts,
                    &self.highlight_bonds,
                    &self.highlight_bond_map,
                    &self.highlight_atom_map,
                );
                let this_high_nbrs = self.find_high_bond_nbrs(this_idx, nbr_idx);
                let nbr_high_nbrs = self.find_high_bond_nbrs(nbr_idx, this_idx);
                let mut points =
                    self.make_highlight_end(this_idx, nbr_idx, line_width, &this_high_nbrs)?;
                points.extend(self.make_highlight_end(
                    nbr_idx,
                    this_idx,
                    line_width,
                    &nbr_high_nbrs,
                )?);
                // The end points are sometimes swapped round, giving a
                // butterfly rather than a rectangle (Github5592). Make a convex
                // hull with a simplified Graham's scan: all points are in the
                // hull, so the usual second step that removes inner points
                // isn't needed. Sort so the lowest y point is first, with
                // lowest x as tie-breaker.
                libcxx_sort_by(&mut points, |p1, p2| {
                    if p1.y < p2.y {
                        true
                    } else if p1.y == p2.y {
                        p1.x < p2.x
                    } else {
                        false
                    }
                });
                // then sort points 1.. anti-clockwise around points[0] by
                // cross products
                let p0 = points[0];
                libcxx_sort_by(&mut points[1..], |p1, p2| {
                    let val = (p1.y - p0.y) * (p2.x - p1.x) - (p1.x - p0.x) * (p2.y - p1.y);
                    if val == 0.0 {
                        (p0 - *p2).length_sq() < (p0 - *p1).length_sq()
                    } else {
                        val < 0.0
                    }
                });
                self.highlights.push(DrawShape::poly_line(
                    points,
                    0.0,
                    false,
                    col,
                    true,
                    (
                        this_idx as i32 + self.active_atm_idx_offset,
                        nbr_idx as i32 + self.active_atm_idx_offset,
                        bond_idx as i32 + self.active_bnd_idx_offset,
                    ),
                    no_dash(),
                )?);
            }
        }
        Ok(())
    }

    /// The 2 or 3 points of the `end1` end of a bond highlight. The final
    /// highlight is formed by calling this for both ends.
    fn make_highlight_end(
        &self,
        end1: u32,
        end2: u32,
        line_width: f64,
        end1_high_nbrs: &[u32],
    ) -> Result<Vec<Point2D>> {
        let half = line_width / 2.0;
        // The intersection of two lines parallel to the lines from e2 to e1 and
        // e3 to e1, half a line width from them: inside the angle they make if
        // pm is 1, outside if -1. If the lines don't intersect, e1.
        let inner_point = |e1: Point2D, e2: Point2D, e3: Point2D, pm: f64| -> Result<Point2D> {
            let perp1 = calc_inner_perpendicular(e2, e1, e3);
            let perp2 = calc_inner_perpendicular(e3, e1, e2);
            let line12 = e2 + mul(perp1, pm * half);
            let mut line11 = e1 + mul(perp1, pm * half);
            line11 = line12 + mul(line12.direction_vector(&line11)?, 2.0 * (e1 - e2).length());
            let line22 = e3 + mul(perp2, pm * half);
            let mut line21 = e1 + mul(perp2, pm * half);
            line21 = line22 + mul(line22.direction_vector(&line21)?, 2.0 * (e1 - e3).length());
            Ok(do_lines_intersect(line12, line11, line22, line21).unwrap_or(e1))
        };
        let end1_cds = self.at_cds[end1 as usize];
        let end2_cds = self.at_cds[end2 as usize];
        let mut points = Vec::new();
        match end1_high_nbrs.len() {
            // no highlighted neighbour bonds: a flat end
            0 => {
                let perp = calc_perpendicular(end1_cds, end2_cds);
                points.push(end1_cds + mul(perp, half));
                points.push(end1_cds - mul(perp, half));
            }
            // one intersection: just a slanted end
            1 => {
                let mut end3_cds = self.at_cds[end1_high_nbrs[0] as usize];
                let b1 = end2_cds.direction_vector(&end1_cds)?;
                let b2 = end2_cds.direction_vector(&end3_cds)?;
                if 1.0 - b1.dot_product(&b2).abs() < 1.0e-4 {
                    // move end3 a little to create an inner and an outer; it
                    // is only used to construct ins1 and ins2
                    let d32 = end3_cds - end2_cds;
                    end3_cds += mul(Point2D::new(d32.y, -d32.x), 0.1);
                }
                points.push(inner_point(end1_cds, end2_cds, end3_cds, 1.0)?);
                points.push(inner_point(end1_cds, end2_cds, end3_cds, -1.0)?);
            }
            // a junction of at least 3 highlights: a triangular end whose point
            // is end1; the other points come from the first and last bond
            // vectors going round from the end1->end2 vector
            _ => {
                let bvec = end1_cds.direction_vector(&end2_cds)?;
                let mut angs: Vec<(usize, f64)> = Vec::new();
                for (i, &n) in end1_high_nbrs.iter().enumerate() {
                    let ovec = end1_cds.direction_vector(&self.at_cds[n as usize])?;
                    angs.push((i, bvec.signed_angle_to(&ovec)?));
                }
                libcxx_sort_by(&mut angs, |a, b| a.1 < b.1);
                let pi = std::f64::consts::PI;
                let (first, last) = (angs[0], angs[angs.len() - 1]);
                // both angles on the same side as end1->end2: the other way
                // round
                if first.1 > pi && last.1 > pi {
                    angs.reverse();
                }
                let (first, last) = (angs[0], angs[angs.len() - 1]);
                let end3_cds = self.at_cds[end1_high_nbrs[first.0] as usize];
                points.push(inner_point(end1_cds, end2_cds, end3_cds, 1.0)?);
                points.push(end1_cds);
                let end4_cds = self.at_cds[end1_high_nbrs[last.0] as usize];
                let pm = if (first.1 > pi && last.1 > pi) || (first.1 < pi && last.1 < pi) {
                    -1.0
                } else {
                    1.0
                };
                points.push(inner_point(end1_cds, end2_cds, end4_cds, pm)?);
            }
        }
        Ok(points)
    }
}
