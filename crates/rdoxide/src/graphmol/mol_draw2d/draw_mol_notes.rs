//! Notes, CIP codes, stereo group labels and radicals of a [`DrawMol`], and
//! the placement of annotations so they avoid the rest of the drawing (the
//! matching methods of `DrawMol.cpp`, RDKit 2026.09.1).
//!
//! Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::details::does_line_intersect;
use super::draw_annotation::{DrawAnnotation, EMPTY_EXTREMES};
use super::draw_mol::DrawMol;
use super::draw_mol_utils::calc_perpendicular;
use super::draw_text::DrawTextFT;
use super::helpers::{OrientType, TextAlignType};
use super::string_rect::StringRect;
use crate::geometry::Point2D;
use crate::graphmol::atropisomers::get_all_atom_ids_for_stereo_group;
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::stereo_group::StereoGroupType;
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result};

fn mul(p: Point2D, f: f64) -> Point2D {
    Point2D::new(p.x * f, p.y * f)
}

impl DrawMol {
    fn string_prop_of_mol(&self, key: &str) -> Option<String> {
        self.draw_mol.get_prop(key).map(|v| v.to_rd_string())
    }

    /// The molecule-level note: the `molNote` property or, if asked for, an
    /// "ABS" for the chiral flag. It goes in a corner of the drawing that
    /// doesn't clash with the molecule.
    pub fn extract_mol_notes(&mut self, text: &mut DrawTextFT) -> Result<()> {
        let mut note = self.string_prop_of_mol(cp::molNote).unwrap_or_default();
        if !self.draw_mol.has_prop(cp::molNote)
            && self.opts.include_chiral_flag_label
            && self
                .draw_mol
                .get_prop(cp::_MolFileChiralFlag)
                .and_then(|v| v.as_int())
                .is_some_and(|f| f != 0)
        {
            note = "ABS".to_string();
        }
        if note.is_empty() {
            return Ok(());
        }
        let colour = self.opts.annotation_colour;
        // molecule annotations use a full-size font
        let tmp = DrawAnnotation::new(
            &note,
            TextAlignType::Start,
            "note",
            1.0,
            Point2D::new(0.0, 0.0),
            colour,
            text,
        )?;
        let (width, height) = tmp.dimensions();
        // Try the corners until there's no clash with the molecule. Even though
        // the alignment is START, the middle of the first char goes at the
        // location, so that is allowed for.
        let half_first = tmp.rects[0].width / 2.0;
        let locs = [
            Point2D::new(f64::from(self.width) - width, height),
            Point2D::new(half_first, height),
            Point2D::new(half_first, f64::from(self.draw_height) - height),
            Point2D::new(
                f64::from(self.width) - width,
                f64::from(self.draw_height) - height,
            ),
        ];
        for loc in locs.iter().take(3) {
            let pos = Point2D::new(loc.x + self.x_offset, loc.y + self.y_offset);
            let annot =
                DrawAnnotation::new(&note, TextAlignType::Start, "note", 1.0, pos, colour, text)?;
            // This goes in the legends because it is already in drawing
            // coordinates, so mustn't be moved by change_to_draw_coords.
            if self.does_note_clash(&annot) == 0 {
                self.legends.push(annot);
                return Ok(());
            }
        }
        // nowhere that doesn't clash: live with it
        self.legends.push(DrawAnnotation::new(
            &note,
            TextAlignType::Start,
            "note",
            1.0,
            locs[0],
            colour,
            text,
        )?);
        Ok(())
    }

    pub fn extract_atom_notes(&mut self, text: &mut DrawTextFT) -> Result<()> {
        for idx in 0..self.draw_mol.num_atoms() {
            let Some(note) = self
                .draw_mol
                .atom(idx)
                .get_prop(cp::atomNote)
                .map(|v| v.to_rd_string())
                .filter(|n| !n.is_empty())
            else {
                continue;
            };
            let mut annot = DrawAnnotation::new(
                &note,
                TextAlignType::Middle,
                "note",
                self.opts.annotation_font_scale,
                Point2D::new(0.0, 0.0),
                self.opts.atom_note_colour,
                text,
            )?;
            self.calc_annotation_position_atom(idx, &mut annot)?;
            self.annotations.push(annot);
        }
        Ok(())
    }

    pub fn extract_cip_codes(
        &mut self,
        show_all_cip_codes: bool,
        text: &mut DrawTextFT,
    ) -> Result<()> {
        let mut masked_atoms = vec![false; self.draw_mol.num_atoms() as usize];
        let mut masked_bonds = vec![false; self.draw_mol.num_bonds() as usize];
        if !show_all_cip_codes {
            // atoms and bonds whose codes are hidden
            for group in self.draw_mol.stereo_groups() {
                if matches!(
                    group.group_type(),
                    StereoGroupType::Or | StereoGroupType::And
                ) {
                    group
                        .atoms()
                        .iter()
                        .for_each(|&a| masked_atoms[a as usize] = true);
                    group
                        .bonds()
                        .iter()
                        .for_each(|&b| masked_bonds[b as usize] = true);
                }
            }
        }
        for idx in 0..self.draw_mol.num_atoms() {
            if masked_atoms[idx as usize] {
                continue;
            }
            let Some(cip) = self
                .draw_mol
                .atom(idx)
                .get_prop(cp::_CIPCode)
                .map(|v| v.to_rd_string())
            else {
                continue;
            };
            let mut annot = DrawAnnotation::new(
                &format!("({cip})"),
                TextAlignType::Middle,
                "CIP_Code",
                self.opts.annotation_font_scale,
                Point2D::new(0.0, 0.0),
                self.opts.atom_note_colour,
                text,
            )?;
            self.calc_annotation_position_atom(idx, &mut annot)?;
            self.annotations.push(annot);
        }
        for idx in 0..self.draw_mol.num_bonds() {
            if masked_bonds[idx as usize] {
                continue;
            }
            let bond = self.draw_mol.bond(idx);
            // E or Z are added if missing to be compatible with the previous
            // implementation; users should call AssignCIPLabels() before
            // drawing to harmonise this with how R, S, M and P are handled
            let cip = match bond.get_prop(cp::_CIPCode).map(|v| v.to_rd_string()) {
                Some(c) => c,
                None => match bond.stereo() {
                    BondStereo::StereoE => "E".to_string(),
                    BondStereo::StereoZ => "Z".to_string(),
                    _ => String::new(),
                },
            };
            if cip.is_empty() {
                continue;
            }
            let mut annot = DrawAnnotation::new(
                &format!("({cip})"),
                TextAlignType::Middle,
                "CIP_Code",
                self.opts.annotation_font_scale,
                Point2D::new(0.0, 0.0),
                self.opts.bond_note_colour,
                text,
            )?;
            self.calc_annotation_position_bond(idx, &mut annot)?;
            self.annotations.push(annot);
        }
        Ok(())
    }

    pub fn extract_stereo_groups(&mut self, text: &mut DrawTextFT) -> Result<()> {
        let (mut or_count, mut and_count) = (0, 0);
        let groups = self.draw_mol.stereo_groups().to_vec();
        for group in &groups {
            let label = match group.group_type() {
                StereoGroupType::Absolute => self.opts.stereo_group_abs_label.clone(),
                StereoGroupType::Or => {
                    or_count += 1;
                    format!("{}{or_count}", self.opts.stereo_group_or_label)
                }
                StereoGroupType::And => {
                    and_count += 1;
                    format!("{}{and_count}", self.opts.stereo_group_and_label)
                }
            };
            for atom in get_all_atom_ids_for_stereo_group(&self.draw_mol, group) {
                let mut annot = DrawAnnotation::new(
                    &label,
                    TextAlignType::Middle,
                    "stereoGroup",
                    self.opts.annotation_font_scale,
                    Point2D::new(0.0, 0.0),
                    self.opts.annotation_colour,
                    text,
                )?;
                self.calc_annotation_position_atom(atom, &mut annot)?;
                self.annotations.push(annot);
            }
        }
        Ok(())
    }

    pub fn extract_bond_notes(&mut self, text: &mut DrawTextFT) -> Result<()> {
        for idx in 0..self.draw_mol.num_bonds() {
            let Some(note) = self
                .draw_mol
                .bond(idx)
                .get_prop(cp::bondNote)
                .map(|v| v.to_rd_string())
                .filter(|n| !n.is_empty())
            else {
                continue;
            };
            let mut annot = DrawAnnotation::new(
                &note,
                TextAlignType::Middle,
                "note",
                self.opts.annotation_font_scale,
                Point2D::new(0.0, 0.0),
                self.opts.bond_note_colour,
                text,
            )?;
            self.calc_annotation_position_bond(idx, &mut annot)?;
            self.annotations.push(annot);
        }
        Ok(())
    }

    pub fn extract_radicals(&mut self) -> Result<()> {
        if !self.opts.include_radicals {
            return Ok(());
        }
        for idx in 0..self.draw_mol.num_atoms() {
            if self.draw_mol.atom(idx).num_radical_electrons() == 0 {
                continue;
            }
            let (rect, orient) = self.calc_radical_rect(idx)?;
            self.radicals.push((rect, orient, idx as i32));
        }
        Ok(())
    }

    /// Scans around the atom for somewhere clear for the annotation, falling
    /// back on the least bad place.
    pub fn calc_annotation_position_atom(
        &self,
        atom: u32,
        annot: &mut DrawAnnotation,
    ) -> Result<()> {
        let start_ang = self.get_note_start_angle(atom)?;
        let at_cds = self.at_cds[atom as usize];
        let rad_step = 0.25;
        let mut least_worst_pos = at_cds;
        let mut least_worst_score = 100;
        for j in 1..4 {
            let note_rad = f64::from(j) * rad_step;
            // experience suggests that with an atom symbol the close in radius
            // won't work
            if j == 1 && self.atom_labels[atom as usize].is_some() {
                continue;
            }
            // scan at 30 degree intervals around the atom
            for i in 0..12 {
                let ang = start_ang + f64::from(i) * 30.0 * PI / 180.0;
                annot.pos.x = at_cds.x + ang.cos() * note_rad;
                annot.pos.y = at_cds.y + ang.sin() * note_rad;
                let score = self.does_note_clash(annot);
                if score == 0 {
                    return Ok(());
                }
                if score < least_worst_score {
                    least_worst_score = score;
                    least_worst_pos = annot.pos;
                }
            }
        }
        annot.pos = least_worst_pos;
        Ok(())
    }

    pub fn calc_annotation_position_bond(
        &self,
        bond: u32,
        annot: &mut DrawAnnotation,
    ) -> Result<()> {
        let b = self.draw_mol.bond(bond);
        let at1_cds = self.at_cds[b.begin_atom_idx() as usize];
        let mut at2_cds = self.at_cds[b.end_atom_idx() as usize];
        // Atoms on top of each other give a NaN perpendicular, with very bad
        // effects (issue 6569): move at2 a little in an arbitrary direction.
        if (at1_cds - at2_cds).length_sq() < 0.0001 {
            at2_cds.x += 0.1;
            at2_cds.y += 0.1;
        }
        let perp = calc_perpendicular(at1_cds, at2_cds);
        let bond_vec = at1_cds.direction_vector(&at2_cds)?;
        let bond_len = (at1_cds - at2_cds).length();
        let offset_step = self.opts.multiple_bond_offset;
        let mut least_worst_pos = mul(at1_cds + at2_cds, 0.5);
        let mut least_worst_score = 100;
        for mo in [0.5, 0.33, 0.66, 0.25, 0.75] {
            let mid = at1_cds + mul(bond_vec, bond_len * mo);
            for j in 1..6 {
                if j == 1 && b.bond_type() as i32 > 1 {
                    continue; // multiple bonds need a bigger offset
                }
                let offset = f64::from(j) * offset_step;
                for sign in [1.0, -1.0] {
                    annot.pos = mid + mul(perp, sign * offset);
                    let score = self.does_note_clash(annot);
                    if score == 0 {
                        return Ok(());
                    }
                    if score < least_worst_score {
                        least_worst_pos = annot.pos;
                        least_worst_score = score;
                    }
                }
            }
        }
        annot.pos = least_worst_pos;
        Ok(())
    }

    /// The angle to start looking for a place for a note at an atom: away
    /// from its bonds.
    // upstream's literal 0.7071 is kept
    #[allow(clippy::approx_constant)]
    pub fn get_note_start_angle(&self, atom: u32) -> Result<f64> {
        let degree = self.draw_mol.atom_degree(atom);
        if degree == 0 {
            return Ok(PI / 2.0);
        }
        let at_cds = self.at_cds[atom as usize];
        let mut bond_vecs: Vec<Point2D> = Vec::new();
        for nbr in self.draw_mol.atom_neighbors(atom) {
            // a neighbour with the same coords gives NaN (issue 6559): use a
            // short arbitrary vector instead
            bond_vecs.push(
                at_cds
                    .direction_vector(&self.at_cds[nbr as usize])
                    .unwrap_or(Point2D::new(0.7071, 0.7071)),
            );
        }
        let mut ret_vec = Point2D::default();
        if bond_vecs.len() == 1 {
            if self.atom_labels[atom as usize].is_none() {
                // perpendicular to the bond: mostly to avoid a zero at the end
                // of a bond to carbon, which looks like a black oxygen in the
                // default font
                ret_vec = Point2D::new(bond_vecs[0].y, -bond_vecs[0].x);
            } else {
                // the opposite end
                ret_vec = mul(bond_vecs[0], -1.0);
            }
        } else if bond_vecs.len() == 2 {
            ret_vec = bond_vecs[0] + bond_vecs[1];
            if ret_vec.length_sq() > 1.0e-6 {
                let atom_ref = self.draw_mol.atom(atom);
                // prefer outside the angle, unless Hs are probably in the way
                if atom_ref.num_implicit_hs().unwrap_or(0) == 0 || atom_ref.atomic_num() == 6 {
                    ret_vec = mul(ret_vec, -1.0);
                }
            } else {
                // a -# or == or some such: take a perpendicular to one of them
                ret_vec = Point2D::new(-bond_vecs[0].y, bond_vecs[0].x);
                ret_vec.normalize()?;
            }
        } else {
            // take 2 that are probably adjacent
            let mut discrim = 4.0 * PI / bond_vecs.len() as f64;
            'outer: for i in 0..bond_vecs.len() - 1 {
                for j in i + 1..bond_vecs.len() {
                    let ang = bond_vecs[i].dot_product(&bond_vecs[j]).acos();
                    if ang < discrim {
                        ret_vec = bond_vecs[i] + bond_vecs[j];
                        // normalize fails on a zero-length sum
                        if ret_vec.normalize().is_err() {
                            continue;
                        }
                        discrim = -1.0;
                        break 'outer;
                    }
                }
            }
            if discrim > 0.0 {
                ret_vec = mul(bond_vecs[0] + bond_vecs[1], -1.0);
            }
        }
        Ok(ret_vec.y.atan2(ret_vec.x))
    }

    /// 0 if the note doesn't clash with anything drawn, else 1-3 for what it
    /// clashes with. Only for finding where to put an annotation, so `annot`
    /// should not yet be in the annotations.
    pub fn does_note_clash(&self, annot: &DrawAnnotation) -> i32 {
        for rect in &annot.rects {
            let mut r = *rect;
            r.trans += annot.pos;
            // if the padding is less than this the letters can fit between the
            // 2 lines of a double bond, giving sub-optimal placements
            let padding = self.scale * 0.04;
            let score = self.does_rect_clash(&r, padding);
            if score != 0 {
                return score;
            }
        }
        0
    }

    /// 0 for no clash, 1 for a bond, 2 for an atom label, 3 for an annotation.
    /// Highlights are no longer checked: that gave bad pictures, and notes look
    /// fine on top of them (issues 5269 and 5195).
    pub fn does_rect_clash(&self, rect: &StringRect, padding: f64) -> i32 {
        // the double bonds themselves, as opposed to the shapes derived from
        // them: atom indices can just fit between the lines of a double bond
        // (Github 5185)
        for b in self.draw_mol.bonds() {
            if b.bond_type() == BondType::Double
                && does_line_intersect(
                    rect,
                    self.at_cds[b.begin_atom_idx() as usize],
                    self.at_cds[b.end_atom_idx() as usize],
                    0.0,
                )
            {
                return 1;
            }
        }
        if self.bonds.iter().any(|b| b.does_rect_clash(rect, padding)) {
            return 1;
        }
        if self
            .atom_labels
            .iter()
            .flatten()
            .any(|l| l.does_rect_clash(rect, padding))
        {
            return 2;
        }
        if self
            .annotations
            .iter()
            .any(|a| a.does_rect_clash(rect, padding))
        {
            return 3;
        }
        0
    }

    /// The rectangle and orientation for the radical spots of an atom.
    pub fn calc_radical_rect(&self, atom: u32) -> Result<(StringRect, OrientType)> {
        let num_rade = self.draw_mol.atom(atom).num_radical_electrons();
        let spot_rad = 0.2 * self.opts.multiple_bond_offset * self.font_scale;
        let mut at_cds = self.at_cds[atom as usize];
        if self.scale != 1.0 {
            at_cds = self.atom_draw_coords(atom as usize)?;
        }
        let orient = self
            .atom_syms
            .get(atom as usize)
            .map(|s| s.1)
            .ok_or_else(|| RDKitError::precondition("bad atom number"))?;
        let rad_size = f64::from(4 * num_rade as i32 - 2) * spot_rad / self.font_scale;
        let (x_min, y_min, x_max, y_max);
        if let Some(label) = &self.atom_labels[atom as usize] {
            let mut ext = EMPTY_EXTREMES;
            label.find_extremes(&mut ext);
            (x_min, x_max, y_min, y_max) = ext;
        } else {
            x_min = at_cds.x - 3.0 * spot_rad;
            x_max = at_cds.x + 3.0 * spot_rad;
            y_min = at_cds.y - 3.0 * spot_rad;
            y_max = at_cds.y + 3.0 * spot_rad;
        }
        let mut rect = StringRect::default();
        let fs = self.font_scale;
        let try_orient = |o: OrientType, rect: &mut StringRect| -> bool {
            match o {
                OrientType::N | OrientType::C => {
                    rect.width = rad_size * fs;
                    rect.height = spot_rad * 3.0;
                    rect.trans.x = at_cds.x;
                    rect.trans.y = y_max + 0.5 * rect.height;
                }
                OrientType::S => {
                    rect.width = rad_size * fs;
                    rect.height = spot_rad * 3.0;
                    rect.trans.x = at_cds.x;
                    rect.trans.y = y_min - 0.5 * rect.height;
                }
                OrientType::E => {
                    rect.trans.x = x_max + 3.0 * spot_rad;
                    rect.trans.y = at_cds.y;
                    rect.width = spot_rad * 1.5;
                    rect.height = rad_size * fs;
                }
                OrientType::W => {
                    rect.trans.x = x_min - 3.0 * spot_rad;
                    rect.trans.y = at_cds.y;
                    rect.width = spot_rad * 1.5;
                    rect.height = rad_size * fs;
                }
            }
            self.does_rect_clash(rect, 0.0) == 0
        };
        if try_orient(orient, &mut rect) {
            return Ok((rect, orient));
        }
        for o in [OrientType::N, OrientType::E, OrientType::S, OrientType::W] {
            if orient != o && try_orient(o, &mut rect) {
                return Ok((rect, o));
            }
        }
        // stick them N irrespective of a clash
        try_orient(OrientType::N, &mut rect);
        Ok((rect, OrientType::N))
    }
}
