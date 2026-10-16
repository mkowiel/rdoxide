//! Bonds of a [`DrawMol`] (the bond methods of `DrawMol.cpp`, RDKit
//! 2026.09.1): standard, query, double, triple, wedged, wavy, dative and zero
//! order bonds, the geometry of multiple bonds, and the clean-up of wedge
//! ends and bond joins.
//!
//! Bond shapes are built in molecule-derived coordinates (the atom
//! coordinates of the `DrawMol`, with y inverted); [`DrawMol::finish_create_draw_objects`]
//! moves them to drawing coordinates later.
//!
//! Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::details::do_lines_intersect;
use super::draw_mol::DrawMol;
use super::draw_mol_utils::{
    adjust_bond_end_for_string, are_bonds_parallel, are_bonds_trans, calc_inner_perpendicular,
    calc_perpendicular, get_highlight_bond_width, is_linear_atom, other_neighbor,
};
use super::draw_shape::{DrawShape, ShapeKind, dashes, dots, no_dash, short_dashes};
use super::helpers::DrawColour;
use crate::geometry::Point2D;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::{RDKitError, Result};

fn mul(p: Point2D, f: f64) -> Point2D {
    Point2D::new(p.x * f, p.y * f)
}

fn div(p: Point2D, f: f64) -> Point2D {
    Point2D::new(p.x / f, p.y / f)
}

type Colours = (DrawColour, DrawColour);

/// The four points of a double bond's two lines (start and finish of each).
type Lines = (Point2D, Point2D, Point2D, Point2D);

impl DrawMol {
    /// Extracts all the bonds, then tidies wedge ends and bond joins.
    pub fn extract_bonds(&mut self) -> Result<()> {
        self.calc_mean_bond_length()?;
        let double_bond_offset = self.opts.multiple_bond_offset * self.mean_bond_length;
        for b in 0..self.draw_mol.num_bonds() {
            let mut is_complex = false;
            if let Some(q) = self.draw_mol.bond(b).query() {
                if q.negate || q.description != "BondOrder" {
                    is_complex = true;
                    self.make_query_bond(b, double_bond_offset)?;
                }
            }
            if !is_complex {
                self.make_standard_bond(b, double_bond_offset)?;
            }
        }
        self.adjust_bonds_on_solid_wedge_ends()?;
        self.smooth_bond_joins()?;
        Ok(())
    }

    /// The colours of the two halves of a bond.
    pub fn get_bond_colours(&self, bond: u32) -> Result<Colours> {
        let highlight_bond = self.highlight_bonds.contains(&(bond as i32));
        if !self.bond_colours.is_empty() {
            return Ok(self.bond_colours[bond as usize]);
        }
        let b = self.draw_mol.bond(bond);
        if !highlight_bond || self.opts.continuous_highlight {
            if self.opts.single_colour_bonds {
                return Ok((self.opts.symbol_colour, self.opts.symbol_colour));
            }
            return Ok((
                self.get_colour(b.begin_atom_idx() as i32)?,
                self.get_colour(b.end_atom_idx() as i32)?,
            ));
        }
        let c = self
            .highlight_bond_map
            .get(&(bond as i32))
            .copied()
            .unwrap_or(self.opts.highlight_colour);
        Ok((c, c))
    }

    fn atom_ids(&self, bond: u32) -> (u32, u32) {
        let b = self.draw_mol.bond(bond);
        (b.begin_atom_idx(), b.end_atom_idx())
    }

    /// Adds a line for a bond; two-colour lines (or split bonds) become two
    /// half lines.
    #[allow(clippy::too_many_arguments)]
    pub fn new_bond_line(
        &mut self,
        pt1: Point2D,
        pt2: Point2D,
        col1: DrawColour,
        col2: DrawColour,
        atom1_idx: i32,
        atom2_idx: i32,
        bond_idx: i32,
        dash_pattern: Vec<f64>,
    ) -> Result<()> {
        let mut scale_width = self.opts.scale_bond_width;
        let mut line_width = self.opts.bond_line_width;
        if !self.opts.continuous_highlight && self.highlight_bonds.contains(&bond_idx) {
            scale_width = self.opts.scale_highlight_bond_width;
            line_width = get_highlight_bond_width(&self.opts, bond_idx, None) / 4.0;
        }
        let a_off = self.active_atm_idx_offset;
        let b_off = self.active_bnd_idx_offset;
        let plain = dash_pattern == no_dash();
        if col1 == col2 && !self.opts.split_bonds {
            self.bonds.push(DrawShape::simple_line(
                vec![pt1, pt2],
                line_width,
                scale_width,
                col1,
                (atom1_idx + a_off, atom2_idx + a_off, bond_idx + b_off),
                dash_pattern,
            )?);
            if plain {
                self.single_bond_lines.push(self.bonds.len() as i32 - 1);
            }
        } else {
            let mid = div(pt1 + pt2, 2.0);
            let at2_idx = if self.opts.split_bonds { -1 } else { atom2_idx };
            self.bonds.push(DrawShape::simple_line(
                vec![pt1, mid],
                line_width,
                scale_width,
                col1,
                (atom1_idx + a_off, at2_idx + a_off, bond_idx + b_off),
                dash_pattern.clone(),
            )?);
            if plain {
                self.single_bond_lines.push(self.bonds.len() as i32 - 1);
            }
            let at1_idx = if self.opts.split_bonds {
                atom2_idx
            } else {
                atom1_idx
            };
            self.bonds.push(DrawShape::simple_line(
                vec![mid, pt2],
                line_width,
                scale_width,
                col2,
                (at1_idx + a_off, at2_idx + a_off, bond_idx + b_off),
                dash_pattern,
            )?);
            if plain {
                self.single_bond_lines.push(self.bonds.len() as i32 - 1);
            }
        }
        Ok(())
    }

    /// `adjustBondEndsForLabels`: the ends of a bond, moved clear of the atom
    /// labels at either end.
    pub fn adjust_bond_ends_for_labels(&self, beg: u32, end: u32) -> (Point2D, Point2D) {
        // the scale factor is empirical
        let mut padding = 0.033 * self.mean_bond_length;
        if self.opts.additional_atom_label_padding > 0.0 {
            padding += self.opts.additional_atom_label_padding;
        }
        let mut beg_cds = self.at_cds[beg as usize];
        let mut end_cds = self.at_cds[end as usize];
        if let Some(l) = &self.atom_labels[beg as usize] {
            adjust_bond_end_for_string(end_cds, padding, l.rects(), &mut beg_cds);
        }
        if let Some(l) = &self.atom_labels[end as usize] {
            adjust_bond_end_for_string(beg_cds, padding, l.rects(), &mut end_cds);
        }
        (beg_cds, end_cds)
    }

    fn atoms_overlap(&self, beg: u32, end: u32) -> bool {
        // If the 2 atoms are on top of each other nothing is drawn: this can
        // give NaN points for things like chiral atoms (issue 6569).
        (self.at_cds[beg as usize] - self.at_cds[end as usize]).length_sq() < 0.0001
    }

    pub fn make_standard_bond(&mut self, bond: u32, double_bond_offset: f64) -> Result<()> {
        let (beg, end) = self.atom_ids(bond);
        if self.atoms_overlap(beg, end) {
            return Ok(());
        }
        let cols = self.get_bond_colours(bond)?;
        let b = self.draw_mol.bond(bond);
        let (bt, dir) = (b.bond_type(), b.bond_dir());
        match bt {
            BondType::Double | BondType::Aromatic => {
                self.make_double_bond_lines(bond, double_bond_offset, cols)
            }
            BondType::Single if matches!(dir, BondDir::BeginWedge | BondDir::BeginDash) => {
                self.make_wedged_bond(bond, cols)
            }
            BondType::Single if dir == BondDir::Unknown => {
                self.make_wavy_bond(bond, double_bond_offset, cols)
            }
            BondType::Dative | BondType::DativeL | BondType::DativeR => {
                self.make_dative_bond(bond, double_bond_offset, cols)
            }
            BondType::Zero => self.make_zero_bond(bond, cols, short_dashes()),
            BondType::Hydrogen => self.make_zero_bond(bond, cols, dots()),
            _ => {
                // in all other cases a line between the two atoms is wanted
                let (end1, end2) = self.adjust_bond_ends_for_labels(beg, end);
                self.new_bond_line(
                    end1,
                    end2,
                    cols.0,
                    cols.1,
                    beg as i32,
                    end as i32,
                    bond as i32,
                    no_dash(),
                )?;
                if bt == BondType::Triple {
                    self.make_triple_bond_lines(bond, double_bond_offset, cols)?;
                }
                Ok(())
            }
        }
    }

    pub fn make_query_bond(&mut self, bond: u32, double_bond_offset: f64) -> Result<()> {
        let q = self
            .draw_mol
            .bond(bond)
            .query()
            .ok_or_else(|| RDKitError::precondition("no query"))?
            .clone();
        let (beg, end) = self.atom_ids(bond);
        let at1_cds = self.at_cds[beg as usize];
        let at2_cds = self.at_cds[end as usize];
        if self.atoms_overlap(beg, end) {
            return Ok(());
        }
        let (end1, end2) = self.adjust_bond_ends_for_labels(beg, end);
        let (sat1, sat2) = (self.at_cds[beg as usize], self.at_cds[end as usize]);
        self.at_cds[beg as usize] = end1;
        self.at_cds[end as usize] = end2;

        let mut midp = div(at2_cds + at1_cds, 2.0);
        let tdash = short_dashes();
        let qc = self.opts.query_colour;
        let split = self.opts.split_bonds;
        let (bi, ei, bd) = (beg as i32, end as i32, bond as i32);
        let mut draw_generic_query = false;
        let (mut at1_idx, mut at2_idx) = (bi, ei);
        let desc = q.description.as_str();
        if desc == "SingleOrDoubleBond" || desc == "SingleOrAromaticBond" {
            at1_idx = bi;
            at2_idx = if split { -1 } else { ei };
            self.new_bond_line(at1_cds, midp, qc, qc, at1_idx, at2_idx, bd, no_dash())?;
            let (l1s, l1f, l2s, l2f) = self.calc_double_bond_lines(double_bond_offset, bond)?;
            at1_idx = if split { ei } else { bi };
            at2_idx = if split { -1 } else { ei };
            midp = div(l1s + l1f, 2.0);
            self.new_bond_line(midp, l1f, qc, qc, at1_idx, at2_idx, bd, no_dash())?;
            midp = div(l2s + l2f, 2.0);
            let second = if desc == "SingleOrDoubleBond" {
                no_dash()
            } else {
                tdash.clone()
            };
            self.new_bond_line(midp, l2f, qc, qc, at1_idx, at2_idx, bd, second)?;
        } else if desc == "DoubleOrAromaticBond" {
            at1_idx = bi;
            at2_idx = if split { -1 } else { ei };
            let (l1s, l1f, l2s, l2f) = self.calc_double_bond_lines(double_bond_offset, bond)?;
            midp = div(l1s + l1f, 2.0);
            self.new_bond_line(l1s, midp, qc, qc, at1_idx, at2_idx, bd, no_dash())?;
            self.new_bond_line(midp, l1f, qc, qc, at1_idx, at2_idx, bd, no_dash())?;
            at1_idx = if split { ei } else { bi };
            at2_idx = if split { -1 } else { ei };
            midp = div(l2s + l2f, 2.0);
            self.new_bond_line(l2s, midp, qc, qc, at1_idx, at2_idx, bd, no_dash())?;
            self.new_bond_line(midp, l2f, qc, qc, at1_idx, at2_idx, bd, tdash)?;
        } else if desc == "BondNull" {
            at1_idx = bi;
            at2_idx = ei;
            self.new_bond_line(at1_cds, at2_cds, qc, qc, at1_idx, at2_idx, bd, tdash)?;
        } else if desc == "BondAnd" && q.children.len() == 2 {
            let (mut q1, mut q2) = (&q.children[0], &q.children[1]);
            if q2.description == "BondOrder" {
                std::mem::swap(&mut q1, &mut q2);
            }
            if q1.description == "BondOrder" && q2.description == "BondInRing" {
                let curr_num_bonds = self.bonds.len();
                self.make_standard_bond(bond, double_bond_offset)?;
                for b in &mut self.bonds[curr_num_bonds..] {
                    b.line_colour = qc;
                }
                let mut segment = at2_cds - at1_cds;
                let (a_off, b_off) = (self.active_atm_idx_offset, self.active_bnd_idx_offset);
                if !q2.negate {
                    segment = div(segment, segment.length() * 6.0);
                    let r1 = Point2D::new(
                        0.5 * segment.x - 0.866 * segment.y,
                        0.866 * segment.x + 0.5 * segment.y,
                    );
                    let r2 = Point2D::new(0.5 * r1.x - 0.866 * r1.y, 0.866 * r1.x + 0.5 * r1.y);
                    let pts = vec![
                        midp + segment,
                        midp + r1,
                        midp + r2,
                        midp - segment,
                        midp - r1,
                        midp - r2,
                        midp + segment,
                    ];
                    self.bonds.push(DrawShape::poly_line(
                        pts,
                        1.0,
                        false,
                        qc,
                        false,
                        (bi + a_off, ei + a_off, bd + b_off),
                        no_dash(),
                    )?);
                } else {
                    segment = div(segment, segment.length() * 10.0);
                    let l = segment.length();
                    let p1 = midp + segment;
                    self.bonds.push(DrawShape::ellipse(
                        vec![p1, Point2D::new(l, l)],
                        1.0,
                        false,
                        qc,
                        false,
                        -1,
                    )?);
                    let p1 = midp - segment;
                    self.bonds.push(DrawShape::ellipse(
                        vec![p1, Point2D::new(l, l)],
                        1.0,
                        false,
                        qc,
                        false,
                        -1,
                    )?);
                }
            } else {
                draw_generic_query = true;
            }
        } else {
            draw_generic_query = true;
        }
        if draw_generic_query {
            self.new_bond_line(at1_cds, at2_cds, qc, qc, at1_idx, at2_idx, bd, dots())?;
            let lw = self.opts.bond_line_width;
            let last = self.bonds.last_mut().expect("a line was just added");
            last.line_width = lw;
            last.scale_line_width = false;
        }
        self.at_cds[beg as usize] = sat1;
        self.at_cds[end as usize] = sat2;
        Ok(())
    }

    pub fn make_double_bond_lines(
        &mut self,
        bond: u32,
        double_bond_offset: f64,
        cols: Colours,
    ) -> Result<()> {
        let (at1, at2) = self.atom_ids(bond);
        let (end1, end2) = self.adjust_bond_ends_for_labels(at1, at2);
        let mut skip_outer_line = false;
        if matches!(
            self.draw_mol.bond(bond).bond_dir(),
            BondDir::BeginWedge | BondDir::BeginDash
        ) {
            self.make_wedged_bond(bond, cols)?;
            skip_outer_line = true;
        }
        let (sat1, sat2) = (self.at_cds[at1 as usize], self.at_cds[at2 as usize]);
        self.at_cds[at1 as usize] = end1;
        self.at_cds[at2 as usize] = end2;
        let (l1s, l1f, l2s, l2f) = self.calc_double_bond_lines(double_bond_offset, bond)?;
        let (a1, a2, bi) = (at1 as i32, at2 as i32, bond as i32);
        if !skip_outer_line {
            self.new_bond_line(l1s, l1f, cols.0, cols.1, a1, a2, bi, no_dash())?;
        }
        if self.draw_mol.bond(bond).bond_type() == BondType::Aromatic {
            self.new_bond_line(l2s, l2f, cols.0, cols.1, a1, a2, bi, dashes())?;
        } else {
            // A two colour simple line is exactly half one colour and half
            // the other, but the second line to a terminal atom (an aldehyde,
            // say) might be asymmetrically shorter, so don't change colour at
            // halfway.
            let l1 = (l1s - l1f).length_sq();
            let l2 = (l2s - l2f).length_sq();
            let beg_deg1 = self.draw_mol.atom_degree(at1) == 1;
            let end_deg1 = self.draw_mol.atom_degree(at2) == 1;
            if (beg_deg1 || end_deg1) && cols.0 != cols.1 && (l1 - l2).abs() > 0.01 {
                let midlen = l1.sqrt() / 2.0;
                let not_mid = if beg_deg1 {
                    l2s + mul(l2s.direction_vector(&l2f)?, midlen)
                } else {
                    l2f + mul(l2f.direction_vector(&l2s)?, midlen)
                };
                self.new_bond_line(l2s, not_mid, cols.0, cols.0, a1, a2, bi, no_dash())?;
                self.new_bond_line(not_mid, l2f, cols.1, cols.1, a1, a2, bi, no_dash())?;
            } else {
                self.new_bond_line(l2s, l2f, cols.0, cols.1, a1, a2, bi, no_dash())?;
            }
        }
        self.at_cds[at1 as usize] = sat1;
        self.at_cds[at2 as usize] = sat2;
        Ok(())
    }

    pub fn make_triple_bond_lines(
        &mut self,
        bond: u32,
        double_bond_offset: f64,
        cols: Colours,
    ) -> Result<()> {
        let (at1, at2) = self.atom_ids(bond);
        let (end1, end2) = self.adjust_bond_ends_for_labels(at1, at2);
        let (sat1, sat2) = (self.at_cds[at1 as usize], self.at_cds[at2 as usize]);
        self.at_cds[at1 as usize] = end1;
        self.at_cds[at2 as usize] = end2;
        let (l1s, l1f, l2s, l2f) = self.calc_triple_bond_lines(double_bond_offset, bond);
        let (a1, a2, bi) = (at1 as i32, at2 as i32, bond as i32);
        self.new_bond_line(l1s, l1f, cols.0, cols.1, a1, a2, bi, no_dash())?;
        self.new_bond_line(l2s, l2f, cols.0, cols.1, a1, a2, bi, no_dash())?;
        self.at_cds[at1 as usize] = sat1;
        self.at_cds[at2 as usize] = sat2;
        Ok(())
    }

    pub fn make_wedged_bond(&mut self, bond: u32, cols: Colours) -> Result<()> {
        let (at1, at2) = self.atom_ids(bond);
        let (mut col1, mut col2) = cols;
        if self.opts.single_colour_wedge_bonds || self.opts.single_colour_bonds {
            col1 = self.opts.symbol_colour;
            col2 = self.opts.symbol_colour;
        }
        // If either atom has a label the padding is made a bit bigger so the
        // end of the wedge doesn't run up to the atom symbol. (This relies on
        // the label padding being a function of the mean bond length.)
        let labelled =
            self.atom_labels[at1 as usize].is_some() || self.atom_labels[at2 as usize].is_some();
        if labelled {
            self.mean_bond_length *= 2.0;
        }
        let (end1, end2) = self.adjust_bond_ends_for_labels(at1, at2);
        if labelled {
            self.mean_bond_length /= 2.0;
        }
        let at1_cds = self.at_cds[at1 as usize];
        let at2_cds = self.at_cds[at2 as usize];
        let perp = calc_perpendicular(at1_cds, at2_cds);
        // The 'fatness' of the wedge is a fraction of the mean bond length, so
        // something is always visible.
        let disp = div(
            mul(perp, self.opts.multiple_bond_offset * self.mean_bond_length),
            2.0,
        );
        let pts = vec![end1, end2 + disp, end2 - disp];
        // deliberately not scaling the highlighted bond width
        let line_width = if self.opts.bond_line_width < 1.0 {
            self.opts.bond_line_width
        } else {
            self.opts.bond_line_width / 2.0
        };
        let ids = (
            at1 as i32 + self.active_atm_idx_offset,
            at2 as i32 + self.active_atm_idx_offset,
            bond as i32 + self.active_bnd_idx_offset,
        );
        let shape = if self.draw_mol.bond(bond).bond_dir() == BondDir::BeginWedge {
            let other_bond_vecs = self.find_other_bond_vecs(at2, at1)?;
            DrawShape::solid_wedge(
                pts,
                col1,
                col2,
                self.opts.split_bonds,
                other_bond_vecs,
                line_width,
                ids,
            )?
        } else {
            let one_less_dash = self.draw_mol.atom_degree(at2) > 1;
            DrawShape::dashed_wedge(pts, col1, col2, one_less_dash, line_width, ids)?
        };
        self.bonds.push(shape);
        Ok(())
    }

    pub fn make_wavy_bond(&mut self, bond: u32, offset: f64, cols: Colours) -> Result<()> {
        let (at1, at2) = self.atom_ids(bond);
        let (end1, end2) = self.adjust_bond_ends_for_labels(at1, at2);
        self.bonds.push(DrawShape::wavy_line(
            vec![end1, end2],
            self.opts.bond_line_width,
            false,
            cols.0,
            cols.1,
            offset,
            (
                at1 as i32 + self.active_atm_idx_offset,
                at2 as i32 + self.active_atm_idx_offset,
                bond as i32 + self.active_bnd_idx_offset,
            ),
        )?);
        Ok(())
    }

    pub fn make_dative_bond(&mut self, bond: u32, offset: f64, cols: Colours) -> Result<()> {
        let (at1, at2) = self.atom_ids(bond);
        let (end1, end2) = self.adjust_bond_ends_for_labels(at1, at2);
        let mid = mul(end1 + end2, 0.5);
        let atid2 = if self.opts.split_bonds { at1 } else { at2 };
        self.new_bond_line(
            end1,
            mid,
            cols.0,
            cols.0,
            at1 as i32,
            atid2 as i32,
            bond as i32,
            no_dash(),
        )?;
        // the fraction of the line that is arrowhead gives a consistent
        // number of pixels
        let frac = 2.0 * offset / (end2 - end1).length();
        self.bonds.push(DrawShape::arrow(
            vec![mid, end2],
            self.opts.bond_line_width,
            false,
            cols.1,
            true,
            (
                at1 as i32 + self.active_atm_idx_offset,
                atid2 as i32 + self.active_atm_idx_offset,
                bond as i32 + self.active_bnd_idx_offset,
            ),
            frac,
            std::f64::consts::PI / 12.0,
        )?);
        Ok(())
    }

    pub fn make_zero_bond(&mut self, bond: u32, cols: Colours, dash: Vec<f64>) -> Result<()> {
        let (at1, at2) = self.atom_ids(bond);
        let (end1, end2) = self.adjust_bond_ends_for_labels(at1, at2);
        self.new_bond_line(
            end1,
            end2,
            cols.0,
            cols.1,
            at1 as i32,
            at2 as i32,
            bond as i32,
            dash,
        )
    }

    /// The lines of a double bond: `(l1s, l1f, l2s, l2f)`.
    pub fn calc_double_bond_lines(&mut self, offset: f64, bond: u32) -> Result<Lines> {
        let (at1, at2) = self.atom_ids(bond);
        let linear = is_linear_atom(&self.draw_mol, at1, &self.at_cds)
            || is_linear_atom(&self.draw_mol, at2, &self.at_cds);
        let (deg1, deg2) = (
            self.draw_mol.atom_degree(at1),
            self.draw_mol.atom_degree(at2),
        );
        let (c1, c2) = (self.at_cds[at1 as usize], self.at_cds[at2 as usize]);
        if linear || (deg1 == 1 && deg2 == 1) {
            let perp = mul(calc_perpendicular(c1, c2), offset * 0.5);
            return Ok((c1 + perp, c2 + perp, c1 - perp, c2 - perp));
        }
        if deg1 == 1 || deg2 == 1 {
            return self.double_bond_terminal(at1, at2, offset);
        }
        let (mut l1s, mut l1f) = (c1, c2);
        let (mut l2s, l2f);
        if self.draw_mol.ring_info().num_bond_rings(bond)? > 0 {
            // in a ring the bond is drawn inside the ring
            (l2s, l2f) = self.bond_inside_ring(bond, offset)?;
        } else if self.atom_labels[at1 as usize].is_some()
            && self.atom_labels[at2 as usize].is_some()
        {
            // atom labels at both ends: straddle the atom-atom vector rather
            // than one line on the vector and the other to the side
            (l1s, l1f, l2s, l2f) = self.double_bond_terminal(at1, at2, offset)?;
        } else {
            (l2s, l2f) = self.bond_non_ring(bond, offset)?;
        }
        // Occasionally (Github6170) bad geometry about a bond can leave the
        // lines crossed, because perpendiculars have become confused; usually
        // when a neighbouring bond is roughly linear with the double bond.
        // This uncrosses them if so.
        if !are_bonds_parallel(l1s, l1f, l2f, l2s, 1.0e-4)? {
            std::mem::swap(&mut l1s, &mut l2s);
        }
        let b = self.draw_mol.bond(bond);
        if b.bond_dir() == BondDir::EitherDouble || b.stereo() == BondStereo::StereoAny {
            // crossed bond
            std::mem::swap(&mut l1s, &mut l2s);
        }
        Ok((l1s, l1f, l2s, l2f))
    }

    /// The inner line of a double bond in a ring: `(l2s, l2f)`.
    fn bond_inside_ring(&self, bond: u32, offset: f64) -> Result<(Point2D, Point2D)> {
        let ring_info = self.draw_mol.ring_info();
        let bond_rings = ring_info.bond_rings();
        let bond_in_rings: Vec<usize> = (0..bond_rings.len())
            .filter(|&i| bond_rings[i].contains(&(bond as i32)))
            .collect();
        let b = self.draw_mol.bond(bond);
        let (beg, end) = (b.begin_atom_idx(), b.end_atom_idx());
        // the ring atom connected to bond_atom that isn't the other end of the
        // bond
        let other_ring_atom = |bond_atom: u32, ring_bonds: &[i32]| -> Option<u32> {
            self.draw_mol
                .atom_bonds(bond_atom)
                .filter(|b2| b2.idx() != bond)
                .find(|b2| ring_bonds.contains(&(b2.idx() as i32)))
                .and_then(|b2| b2.other_atom_idx(bond_atom).ok())
        };
        let ring_to_use: &[i32] = if bond_in_rings.len() > 1 {
            // Choose a ring of the same aromaticity as the bond, so an
            // aromatic bond's double line is inside the aromatic ring (as in
            // morphine, with fused aromatic and aliphatic rings).
            let mut chosen = &bond_rings[bond_in_rings[0]];
            for &i in &bond_in_rings {
                chosen = &bond_rings[i];
                if chosen
                    .iter()
                    .all(|&bi| self.draw_mol.bond(bi as u32).is_aromatic() == b.is_aromatic())
                {
                    break;
                }
            }
            chosen
        } else {
            &bond_rings[bond_in_rings[0]]
        };
        let (Some(third), Some(fourth)) = (
            other_ring_atom(beg, ring_to_use),
            other_ring_atom(end, ring_to_use),
        ) else {
            return Err(RDKitError::invariant("ring atom not found"));
        };
        // Bonds in rings can be trans (#5486) while the default code assumes
        // cis. Treat as a non-ring bond if so; it needn't be on the inside of
        // the ring, but that is complicated to fix.
        let c = &self.at_cds;
        if are_bonds_trans(
            c[third as usize],
            c[beg as usize],
            c[end as usize],
            c[fourth as usize],
        ) {
            return self.bond_non_ring(bond, offset);
        }
        let l2s = self.double_bond_end(
            third,
            beg,
            end,
            offset,
            self.atom_labels[beg as usize].is_none(),
        );
        let l2f = self.double_bond_end(
            fourth,
            end,
            beg,
            offset,
            self.atom_labels[end as usize].is_none(),
        );
        Ok((l2s?, l2f?))
    }

    /// The inner line of a double bond not in a ring: `(l2s, l2f)`.
    fn bond_non_ring(&self, bond: u32, offset: f64) -> Result<(Point2D, Point2D)> {
        let (beg, end) = self.atom_ids(bond);
        let m = &self.draw_mol;
        let c = &self.at_cds;
        let beg_trunc = self.atom_labels[beg as usize].is_none();
        let end_trunc = self.atom_labels[end as usize].is_none();
        let nbr = |a: u32, not: u32, n: usize| -> Result<u32> {
            other_neighbor(m, a, not, n).ok_or_else(|| RDKitError::invariant("no neighbour"))
        };
        // a neighbour of at1 that isn't at2 and if possible isn't directly
        // opposite at1 to at2
        let non_colinear_nbor = |at1: u32, at2: u32| -> Result<u32> {
            let mut third = None;
            for i in 1..m.atom_degree(at1) as usize {
                third = other_neighbor(m, at1, at2, i);
                if let Some(t) = third {
                    if !are_bonds_parallel(
                        c[at1 as usize],
                        c[at2 as usize],
                        c[at1 as usize],
                        c[t as usize],
                        1.0e-4,
                    )? {
                        return Ok(t);
                    }
                }
            }
            match third {
                Some(t) => Ok(t),
                // we need something
                None => nbr(at1, at2, 1),
            }
        };
        let (deg_b, deg_e) = (m.atom_degree(beg), m.atom_degree(end));
        let trans = |third: u32, fourth: u32| {
            are_bonds_trans(
                c[third as usize],
                c[beg as usize],
                c[end as usize],
                c[fourth as usize],
            )
        };
        let (l2s, l2f);
        if deg_b == 2 && deg_e == 2 {
            let third = nbr(beg, end, 0)?;
            let fourth = nbr(end, beg, 0)?;
            l2s = self.double_bond_end(third, beg, end, offset, beg_trunc)?;
            if trans(third, fourth) {
                let perp =
                    calc_inner_perpendicular(c[end as usize], c[beg as usize], c[third as usize]);
                l2f = c[end as usize] + mul(perp, offset);
            } else {
                l2f = self.double_bond_end(fourth, end, beg, offset, end_trunc)?;
            }
        } else if deg_b == 2 && deg_e > 2 {
            let third = nbr(beg, end, 0)?;
            let mut fourth = nbr(end, beg, 0)?;
            l2s = self.double_bond_end(third, beg, end, offset, beg_trunc)?;
            if trans(third, fourth) {
                fourth = non_colinear_nbor(end, beg)?;
            }
            l2f = self.double_bond_end(fourth, end, beg, offset, end_trunc)?;
        } else if deg_b > 2 && deg_e == 2 {
            let mut third = nbr(beg, end, 0)?;
            let fourth = nbr(end, beg, 0)?;
            let mut s = self.double_bond_end(third, beg, end, offset, beg_trunc)?;
            if trans(third, fourth) {
                third = non_colinear_nbor(beg, end)?;
                s = self.double_bond_end(third, beg, end, offset, end_trunc)?;
            }
            l2s = s;
            l2f = self.double_bond_end(fourth, end, beg, offset, end_trunc)?;
        } else if deg_b > 2 && deg_e > 2 {
            let third = nbr(beg, end, 0)?;
            l2s = self.double_bond_end(third, beg, end, offset, beg_trunc)?;
            let mut fourth = nbr(end, beg, 0)?;
            if trans(third, fourth) {
                fourth = non_colinear_nbor(end, beg)?;
            }
            l2f = self.double_bond_end(fourth, end, beg, offset, end_trunc)?;
        } else {
            return Err(RDKitError::invariant(
                "double bond between atoms of low degree",
            ));
        }
        Ok((l2s, l2f))
    }

    /// A terminal double bond between `at1` and `at2`, one of which has
    /// degree 1.
    fn double_bond_terminal(
        &mut self,
        mut at1: u32,
        mut at2: u32,
        mut offset: f64,
    ) -> Result<Lines> {
        let mut swapped = false;
        if self.draw_mol.atom_degree(at1) > 1 && self.draw_mol.atom_degree(at2) == 1 {
            std::mem::swap(&mut at1, &mut at2);
            swapped = true;
        }
        let at1_cds = self.at_cds[at1 as usize];
        let at2_cds = self.at_cds[at2 as usize];
        let (mut l1s, mut l1f, mut l2s, mut l2f);
        if self.atom_labels[at2 as usize].is_some() {
            // either side of the bond line if going to a label
            offset /= 2.0;
            let perp = mul(calc_perpendicular(at1_cds, at2_cds), offset);
            (l1s, l1f, l2s, l2f) = (
                at1_cds + perp,
                at2_cds + perp,
                at1_cds - perp,
                at2_cds - perp,
            );
        } else if self.draw_mol.atom_degree(at2) > 2 {
            // Lines either side of the bond line, extended at the at2 end to
            // the intersection of the other bonds from at2; only half the
            // offset is needed.
            offset /= 2.0;
            let perp = mul(calc_perpendicular(at1_cds, at2_cds), offset);
            (l1s, l1f, l2s, l2f) = (
                at1_cds + perp,
                at2_cds + perp,
                at1_cds - perp,
                at2_cds - perp,
            );
            let bl = (l1s - l1f).length().max((l2s - l2f).length());
            l1f = l1s + mul(l1s.direction_vector(&l1f)?, 2.0 * bl);
            l2f = l2s + mul(l2s.direction_vector(&l2f)?, 2.0 * bl);
            for nbr in self.draw_mol.atom_neighbors(at2) {
                let nbr_cds = self.at_cds[nbr as usize];
                if let Some(ip) = do_lines_intersect(l1s, l1f, at2_cds, nbr_cds) {
                    l1f = ip;
                }
                if let Some(ip) = do_lines_intersect(l2s, l2f, at2_cds, nbr_cds) {
                    l2f = ip;
                }
            }
        } else {
            // one line as normal, the 2nd truncates at the internal end only
            l1s = at1_cds;
            l1f = at2_cds;
            let third = other_neighbor(&self.draw_mol, at2, at1, 0)
                .ok_or_else(|| RDKitError::invariant("no third atom"))?;
            let perp = calc_inner_perpendicular(at1_cds, at2_cds, self.at_cds[third as usize]);
            l2s = at1_cds + mul(perp, offset);
            l2f = self.double_bond_end(at1, at2, third, offset, true)?;
            // If at1->at2->at3 is a straight line, l2f may have ended up on
            // the wrong side of the other bond from l2s, as there is no inner
            // side. Do it again with a negative offset if so.
            if l1s
                .direction_vector(&l1f)?
                .dot_product(&l2s.direction_vector(&l2f)?)
                .abs()
                < 0.9999
            {
                l2f = self.double_bond_end(at1, at2, third, -offset, true)?;
            }
            // if at1 has a label it must be centred between the two lines
            // (Github 5511)
            if let Some(l) = &mut self.atom_labels[at1 as usize] {
                l.cds += mul(perp, offset * 0.5);
            }
        }
        if swapped {
            std::mem::swap(&mut l1s, &mut l1f);
            std::mem::swap(&mut l2s, &mut l2f);
        }
        Ok((l1s, l1f, l2s, l2f))
    }

    /// With at1 bonded to at2 and at2 to at3, the position of the at2 end of
    /// a double bond between at2 and at3. If `trunc`, it is along the vector
    /// that bisects the two bonds on the inside, otherwise perpendicular to
    /// the bond from at1 to at2.
    fn double_bond_end(
        &self,
        at1: u32,
        at2: u32,
        at3: u32,
        offset: f64,
        trunc: bool,
    ) -> Result<Point2D> {
        let c = &self.at_cds;
        let v21 = c[at2 as usize].direction_vector(&c[at1 as usize])?;
        let v23 = c[at2 as usize].direction_vector(&c[at3 as usize])?;
        let mut v23perp = Point2D::new(-v23.y, v23.x);
        v23perp.normalize()?;
        let mut bis = v21 + v23;
        if bis.normalize().is_err() {
            // Colinear bonds give a zero bisector which would normalise to
            // NaN (Github #6027). v23perp is safe so long as it is on the
            // right side of the bond, which is checked on return.
            return Ok(c[at2 as usize] - mul(v23perp, offset));
        }
        if v23perp.dot_product(&bis) < 0.0 {
            v23perp = mul(v23perp, -1.0);
        }
        // With an atom label the bond needn't be stepped back because both
        // ends are shortened for the letters; likewise if the two lines don't
        // intersect, it is stepped back enough already (github 6025).
        let ip = if trunc {
            do_lines_intersect(
                c[at2 as usize],
                c[at2 as usize] + bis,
                c[at2 as usize] + mul(v23perp, offset),
                c[at3 as usize] + mul(v23perp, offset),
            )
        } else {
            None
        };
        Ok(ip.unwrap_or(c[at2 as usize] + mul(v23perp, offset)))
    }

    fn calc_triple_bond_lines(&self, offset: f64, bond: u32) -> Lines {
        let (at1, at2) = self.atom_ids(bond);
        let (c1, c2) = (self.at_cds[at1 as usize], self.at_cds[at2 as usize]);
        let perp = calc_perpendicular(c1, c2);
        (
            c1 + mul(perp, offset),
            c2 + mul(perp, offset),
            c1 - mul(perp, offset),
            c2 - mul(perp, offset),
        )
    }

    /// The vectors of the single bonds on `atom` other than the one to
    /// `other_atom`, used to shape the end of a wedge.
    pub fn find_other_bond_vecs(&self, atom: u32, other_atom: u32) -> Result<Vec<Point2D>> {
        let mut res = Vec::new();
        let m = &self.draw_mol;
        if m.atom_degree(atom) == 1 || self.atom_labels[atom as usize].is_some() {
            return Ok(res);
        }
        for i in 1..m.atom_degree(atom) as usize {
            let third = other_neighbor(m, atom, other_atom, i - 1)
                .ok_or_else(|| RDKitError::invariant("no third atom"))?;
            let bond = m
                .get_bond_between_atoms(atom, third)?
                .ok_or_else(|| RDKitError::invariant("no bond"))?;
            // Nothing for a wedge to a triple bond: it gets really messed up,
            // especially if the two bonds aren't exactly co-linear (Github
            // 7620).
            if bond.bond_type() == BondType::Triple {
                continue;
            }
            // A double bond that straddles the atom-atom vector also looks odd
            // or wrong (Github 7739).
            if bond.bond_type() == BondType::Double
                && m.atom_degree(atom) > 2
                && m.atom_degree(third) == 1
            {
                continue;
            }
            res.push(self.at_cds[atom as usize].direction_vector(&self.at_cds[third as usize])?);
        }
        Ok(res)
    }

    /// Moves the central line of a bond that continues from the end of a solid
    /// wedge to the corner of the wedge.
    pub fn adjust_bonds_on_solid_wedge_ends(&mut self) -> Result<()> {
        for b in 0..self.draw_mol.num_bonds() {
            let bond = self.draw_mol.bond(b);
            let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
            if bond.bond_dir() != BondDir::BeginWedge
                || self.draw_mol.atom_degree(end) != 2
                || self.atom_labels[end as usize].is_some()
            {
                continue;
            }
            // the bond at the end atom
            let Some(third) = other_neighbor(&self.draw_mol, end, beg, 0) else {
                continue;
            };
            let Some(bond1) = self.draw_mol.get_bond_between_atoms(end, third)? else {
                continue;
            };
            // moving the central line of a triple bond to the wedge corner is
            // clearly wrong
            if bond1.bond_type() == BondType::Triple {
                continue;
            }
            let bond1_idx = bond1.idx();
            let (b1_beg, b1_end) = (bond1.begin_atom_idx(), bond1.end_atom_idx());
            // co-linear bonds: nothing to do (Github7036)
            let c = &self.at_cds;
            let v1 = c[end as usize].direction_vector(&c[beg as usize])?;
            let v2 = c[b1_end as usize].direction_vector(&c[b1_beg as usize])?;
            if (1.0 - v1.dot_product(&v2)).abs() < 0.001 {
                continue;
            }
            let end_cds = c[end as usize];
            let mut wedge: Option<usize> = None;
            let mut bond_line: Option<usize> = None;
            let mut closest = 1.0;
            for (i, shape) in self.bonds.iter().enumerate() {
                if shape.bond == b as i32 {
                    wedge = Some(i);
                }
                // there may be several lines for the bond: take one with an
                // end as close as possible to the bond end atom
                if shape.bond == bond1_idx as i32 {
                    if !matches!(shape.kind, ShapeKind::SimpleLine { .. }) {
                        continue;
                    }
                    for p in 0..2 {
                        let d = (shape.points[p] - end_cds).length_sq();
                        if d < closest {
                            closest = d;
                            bond_line = Some(i);
                        }
                    }
                }
            }
            let (Some(w), Some(l)) = (wedge, bond_line) else {
                continue;
            };
            // the points at the top of the wedge (this assumes the order the
            // triangles are made in the solid wedge)
            let (mut p1, mut p2) = match self.bonds[w].points.len() {
                3 => (1, 2),
                9 => (4, 5),
                _ => continue,
            };
            // p1 should be the one furthest from the 3rd atom
            let third_cds = self.at_cds[third as usize];
            if (third_cds - self.bonds[w].points[p1]).length_sq()
                < (third_cds - self.bonds[w].points[p2]).length_sq()
            {
                std::mem::swap(&mut p1, &mut p2);
            }
            let target = self.bonds[w].points[p1];
            let wedge_atom2 = self.bonds[w].atom2;
            if self.bonds[l].atom1 == wedge_atom2 {
                self.bonds[l].points[0] = target;
            } else {
                self.bonds[l].points[1] = target;
            }
        }
        Ok(())
    }

    /// Bonds are drawn as individual lines rather than paths, so they don't
    /// join up nicely: put a little polyline round the join where needed. (They
    /// aren't paths because SVG lines carry classes for their atoms and bond,
    /// which people use to identify them.)
    pub fn smooth_bond_joins(&mut self) -> Result<()> {
        for atom in 0..self.draw_mol.num_atoms() {
            // an atom label means no join
            if self.atom_labels[atom as usize].is_some() {
                continue;
            }
            let degree = self.draw_mol.atom_degree(atom);
            let mut do_it = degree == 2;
            if degree == 3 {
                for nbr in self.draw_mol.atom_neighbors(atom) {
                    let bond = self
                        .draw_mol
                        .get_bond_between_atoms(atom, nbr)?
                        .ok_or_else(|| RDKitError::invariant("no bond"))?;
                    if (self.draw_mol.atom_degree(nbr) == 1 && bond.bond_type() == BondType::Double)
                        || matches!(bond.bond_dir(), BondDir::BeginDash | BondDir::BeginWedge)
                    {
                        do_it = true;
                    }
                }
            }
            if !do_it {
                continue;
            }
            let adj_atom_idx = atom as i32 + self.active_atm_idx_offset;
            let end_of = |s: &DrawShape| -> Option<usize> {
                if adj_atom_idx == s.atom1 {
                    Some(0)
                } else if adj_atom_idx == s.atom2 {
                    Some(1)
                } else {
                    None
                }
            };
            let lines: Vec<usize> = self.single_bond_lines.iter().map(|&i| i as usize).collect();
            let mut done = false;
            for (i, &li) in lines.iter().enumerate() {
                let Some(p1) = end_of(&self.bonds[li]) else {
                    continue;
                };
                for (j, &lj) in lines.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    let Some(p2) = end_of(&self.bonds[lj]) else {
                        continue;
                    };
                    let (s1, s2) = (&self.bonds[li], &self.bonds[lj]);
                    if (s1.points[p1] - s2.points[p2]).length_sq() < 1.0e-6 {
                        // A small polyline papers over the crack. Lines of
                        // different colours go less far round the corner so one
                        // colour doesn't extend into the other; the same colour
                        // goes further to hide the join. The numbers are
                        // empirical.
                        let p12 = 1 - p1;
                        let p22 = 1 - p2;
                        let len = if s1.line_colour == s2.line_colour {
                            0.05
                        } else {
                            0.025
                        };
                        let dv1 = mul(s1.points[p1] - s1.points[p12], len);
                        let dv2 = mul(s1.points[p1] - s2.points[p22], len);
                        let pl_pts = vec![s1.points[p1] - dv1, s1.points[p1], s1.points[p1] - dv2];
                        let pl = DrawShape::poly_line(
                            pl_pts,
                            s1.line_width,
                            s1.scale_line_width,
                            s1.line_colour,
                            false,
                            (-1, -1, -1),
                            no_dash(),
                        )?;
                        self.bonds.push(pl);
                        done = true;
                        break;
                    }
                }
                if done {
                    break;
                }
            }
        }
        Ok(())
    }
}
