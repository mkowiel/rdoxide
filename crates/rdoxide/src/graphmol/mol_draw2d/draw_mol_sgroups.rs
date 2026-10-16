//! S-group data, brackets, link nodes, variable attachment bonds, close
//! contacts and the resolution of clashing atom labels of a [`DrawMol`] (the
//! matching methods of `DrawMol.cpp` and `getSGroupDataLabels` of
//! `MolDraw2DDetails.cpp`, RDKit 2026.09.1).
//!
//! Copyright (C) 2015-2026 Greg Landrum, David Cosgrove and other RDKit
//! contributors. Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::details::get_bracket_points;
use super::draw_annotation::{AtomSymbol, DrawAnnotation};
use super::draw_mol::DrawMol;
use super::draw_shape::{DrawShape, no_dash};
use super::draw_text::DrawTextFT;
use super::helpers::{DrawColour, OrientType, TextAlignType};
use crate::geometry::{Point2D, Transform2D};
use crate::graphmol::mol_enumerator::link_node::get_mol_link_nodes;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};
use crate::{RDKitError, Result, rd_log_warning};

fn mul(p: Point2D, f: f64) -> Point2D {
    Point2D::new(p.x * f, p.y * f)
}

/// The text and position of a DAT S-group label.
#[derive(Clone, Debug, PartialEq)]
pub struct SGroupDataLabel {
    pub text: String,
    /// In molecule coordinates.
    pub pos: Point2D,
    /// True if `pos` came from FIELDDISP, false if it is the position of the
    /// associated atom (the fallback).
    pub positioned: bool,
    /// The associated atom, -1 if none.
    pub atom_idx: i32,
}

/// `stripSpacesAndCast<double>`
fn cast_double(s: &str) -> Result<f64> {
    s.trim()
        .parse::<f64>()
        .map_err(|_| RDKitError::Value(format!("Cannot convert '{s}' to a double")))
}

/// `getSGroupDataLabels`: the text and positions of the DAT S-group labels of
/// a molecule, rotated by `rotate` degrees.
pub fn get_sgroup_data_labels(mol: &ROMol, rotate: f64) -> Result<Vec<SGroupDataLabel>> {
    let mut result = Vec::new();
    if mol.num_conformers() == 0 || mol.substance_groups().is_empty() {
        return Ok(result);
    }
    let rot = rotate * PI / 180.0;
    let mut tform = Transform2D::identity();
    tform.set_transform_rotation(&Point2D::new(0.0, 0.0), rot);
    let conf = mol.conformer(-1)?;
    for sg in mol.substance_groups() {
        if sg.get_prop("TYPE").map(|v| v.to_rd_string()).as_deref() != Some("DAT") {
            continue;
        }
        let mut text = String::new();
        if let Some(RDValue::VecString(dfs)) = sg.get_prop("DATAFIELDS") {
            text = dfs.iter().map(|df| format!("{df}|")).collect();
            text.pop();
        }
        if text.is_empty() {
            continue;
        }
        let atom_idx = sg.atoms().first().map_or(-1, |&a| a as i32);
        let mut located = false;
        let mut pos = Point2D::new(0.0, 0.0);
        if let Some(fd) = sg.get_prop("FIELDDISP").map(|v| v.to_rd_string()) {
            let field = |a: usize, b: usize| fd.get(a..b.min(fd.len())).unwrap_or("");
            let xp = cast_double(field(0, 10))?;
            let yp = cast_double(field(10, 20))?;
            // the y of the molecule coordinates is always inverted
            pos = Point2D::new(xp, -yp);
            if fd.as_bytes().get(25) == Some(&b'R') {
                if atom_idx < 0 {
                    // no atom to anchor the relative position to
                    continue;
                } else if xp.abs() > 1e-3 || yp.abs() > 1e-3 {
                    // opposite sign for y
                    let ap = conf.atom_pos(atom_idx as u32)?;
                    pos.x += ap.x;
                    pos.y -= ap.y;
                    located = true;
                }
            } else {
                // absolute position: allow for a centroid offset set by the
                // drawing pipeline
                if mol.has_prop("_centroidx") {
                    let cx = mol.get_prop("_centroidx").and_then(RDValue::as_double);
                    let cy = mol.get_prop("_centroidy").and_then(RDValue::as_double);
                    pos.x += cx.unwrap_or(0.0);
                    pos.y -= cy.unwrap_or(0.0);
                }
                located = true;
            }
            tform.transform_point(&mut pos);
        }
        if !located {
            if atom_idx >= 0 {
                let p = conf.atom_pos(atom_idx as u32)?;
                pos = Point2D::new(p.x, p.y);
            } else {
                rd_log_warning!(
                    "FIELDDISP info not found for DAT SGroup which isn't associated with an atom. SGroup will not be included.\n"
                );
                continue;
            }
        }
        result.push(SGroupDataLabel {
            text,
            pos,
            positioned: located,
            atom_idx,
        });
    }
    Ok(result)
}

/// The label at the bottom (really the top, after the y inversion) of a
/// bracket.
fn draw_bottom_label(
    label: &str,
    brk: &DrawShape,
    annotation_font_scale: f64,
    text: &mut DrawTextFT,
    horizontal: bool,
) -> Result<DrawAnnotation> {
    // annotations go on the last bracket of an S-group
    let p = &brk.points;
    let (mut top_pt, mut brk_pt) = (p[1], p[0]);
    if (!horizontal && p[2].y > top_pt.y) || (horizontal && p[2].x < top_pt.x) {
        top_pt = p[2];
        brk_pt = p[3];
    }
    let mut da = DrawAnnotation::new(
        label,
        TextAlignType::Middle,
        "connect",
        annotation_font_scale,
        top_pt + (top_pt - brk_pt),
        DrawColour::new(0.0, 0.0, 0.0),
        text,
    )?;
    if brk_pt.x < top_pt.x {
        da.align = TextAlignType::Start;
    }
    Ok(da)
}

/// The values of a V3000 array `(n v1 v2 ...)`.
fn parse_v3000_array(s: &str) -> Result<Vec<u32>> {
    let inner = s.trim().trim_start_matches('(').trim_end_matches(')');
    let mut it = inner.split_whitespace();
    let count: usize = it
        .next()
        .and_then(|c| c.parse().ok())
        .ok_or_else(|| RDKitError::Value("invalid count value in V3000 array".into()))?;
    let values: std::result::Result<Vec<u32>, _> = it.take(count).map(str::parse).collect();
    values.map_err(|_| RDKitError::Value("invalid value in V3000 array".into()))
}

impl DrawMol {
    pub fn extract_sgroup_data(&mut self, text: &mut DrawTextFT) -> Result<()> {
        if !self.include_annotations {
            return Ok(());
        }
        // FIELDNAME could be rendered, but Marvin Sketch, Biovia Draw and
        // ChemDraw don't
        for lbl in get_sgroup_data_labels(&self.draw_mol, self.opts.rotate)? {
            // everybody renders these left justified
            let mut annot = DrawAnnotation::new(
                &lbl.text,
                TextAlignType::Start,
                "note",
                self.opts.annotation_font_scale,
                Point2D::new(0.0, 0.0),
                self.opts.annotation_colour,
                text,
            )?;
            if lbl.positioned {
                annot.pos = lbl.pos;
            } else {
                self.calc_annotation_position_atom(lbl.atom_idx as u32, &mut annot)?;
            }
            self.annotations.push(annot);
        }
        Ok(())
    }

    /// Grey blobs and fat lines for the atoms and bonds of variable
    /// attachment points.
    pub fn extract_variable_bonds(&mut self) -> Result<()> {
        let n_atoms = self.draw_mol.num_atoms();
        for b in 0..self.draw_mol.num_bonds() {
            let bond = self.draw_mol.bond(b);
            let (Some(endpts), Some(_attach)) = (
                bond.get_prop(cp::_MolFileBondEndPts)
                    .map(|v| v.to_rd_string()),
                bond.get_prop(cp::_MolFileBondAttach)
                    .map(|v| v.to_rd_string()),
            ) else {
                continue;
            };
            let beg = bond.begin_atom_idx();
            // FIX upstream: maybe distinguish "ANY" and "ALL" attach values
            let mut involved = vec![false; n_atoms as usize];
            for oat in parse_v3000_array(&endpts)? {
                // decrement the indices, with error checking
                if oat == 0 || oat > n_atoms {
                    return Err(RDKitError::Value("Bad variation point index".into()));
                }
                let oat = oat - 1;
                involved[oat as usize] = true;
                let r = self.opts.variable_atom_radius;
                self.pre_shapes.push(DrawShape::ellipse(
                    vec![self.at_cds[oat as usize], Point2D::new(r, r)],
                    1.0,
                    true,
                    self.opts.variable_attachment_colour,
                    true,
                    oat as i32,
                )?);
            }
            for b2 in 0..self.draw_mol.num_bonds() {
                let bond2 = self.draw_mol.bond(b2);
                let (b2_beg, b2_end) = (bond2.begin_atom_idx(), bond2.end_atom_idx());
                if involved[b2_beg as usize] && involved[b2_end as usize] {
                    self.pre_shapes.push(DrawShape::simple_line(
                        vec![self.at_cds[b2_beg as usize], self.at_cds[b2_end as usize]],
                        f64::from(self.opts.variable_bond_width_multiplier),
                        true,
                        self.opts.variable_attachment_colour,
                        (
                            b2_beg as i32 + self.active_atm_idx_offset,
                            b2_end as i32 + self.active_atm_idx_offset,
                            b2 as i32 + self.active_bnd_idx_offset,
                        ),
                        no_dash(),
                    )?);
                }
            }
            // correct the symbol of the end atom (remove the *)
            if self.draw_mol.atom(beg).atomic_num() == 0 {
                self.atom_syms[beg as usize] = (String::new(), OrientType::C);
                self.atom_labels[beg as usize] = None;
            }
        }
        Ok(())
    }

    /// Brackets of S-groups, with their CONNECT and LABEL/TYPE annotations.
    pub fn extract_brackets(&mut self, text: &mut DrawTextFT) -> Result<()> {
        let sgs = self.draw_mol.substance_groups().to_vec();
        if sgs.is_empty() {
            return Ok(());
        }
        // details of this transformation are in extract_atom_coords
        let rot = self.opts.rotate * PI / 180.0;
        let mut trans = Transform2D::identity();
        trans.set_transform_rotation(&Point2D::new(0.0, 0.0), rot);
        let bracket_pts = |brk: &[crate::geometry::Point3D; 3]| {
            // the atom coords have been inverted in y, so the bracket coords
            // must be, too
            let mut p1 = Point2D::new(brk[0].x, -brk[0].y);
            let mut p2 = Point2D::new(brk[1].x, -brk[1].y);
            trans.transform_point(&mut p1);
            trans.transform_point(&mut p2);
            (p1, p2)
        };
        for sg in &sgs {
            if sg.brackets().is_empty() {
                continue;
            }
            // The reference point that shows which way the bracket points.
            // (Thanks to John Mayfield:
            // http://efficientbits.blogspot.com/2015/11/bringing-molfile-sgroups-to-cdk.html)
            let mut ref_pt = Point2D::new(0.0, 0.0);
            if !sg.atoms().is_empty() {
                // The average position of the S-group atoms is a bit simplistic
                // (Github5768: a long chain stretching outside the bracket area
                // turns the last bracket the wrong way), so only the atoms
                // inside the brackets are used, rather crudely.
                let (mut x_min, mut y_min) = (f64::MAX / 2.0, f64::MAX / 2.0);
                let (mut x_max, mut y_max) = (f64::MIN / 2.0, f64::MIN / 2.0);
                for brk in sg.brackets() {
                    let (p1, p2) = bracket_pts(brk);
                    x_min = x_min.min(p1.x).min(p2.x);
                    y_min = y_min.min(p1.y).min(p2.y);
                    x_max = x_max.max(p1.x).max(p2.x);
                    y_max = y_max.max(p1.y).max(p2.y);
                }
                let mut num_in = 0;
                for &a in sg.atoms() {
                    let c = self.at_cds[a as usize];
                    if c.x >= x_min && c.x <= x_max && c.y >= y_min && c.y <= y_max {
                        ref_pt += c;
                        num_in += 1;
                    }
                }
                if num_in > 0 {
                    ref_pt = mul(ref_pt, 1.0 / f64::from(num_in));
                } else {
                    // all of them, and live with the consequences
                    for &a in sg.atoms() {
                        ref_pt += self.at_cds[a as usize];
                    }
                    ref_pt = mul(ref_pt, 1.0 / sg.atoms().len() as f64);
                }
            }
            let mut sg_bond_segments: Vec<(Point2D, Point2D)> = Vec::new();
            for &bnd_idx in sg.bonds() {
                let bnd = self.draw_mol.bond(bnd_idx);
                let (beg, end) = (bnd.begin_atom_idx(), bnd.end_atom_idx());
                if sg.atoms().contains(&beg) {
                    sg_bond_segments.push((self.at_cds[beg as usize], self.at_cds[end as usize]));
                } else if sg.atoms().contains(&end) {
                    sg_bond_segments.push((self.at_cds[end as usize], self.at_cds[beg as usize]));
                }
            }
            let mut num_brackets = 0usize;
            for brk in sg.brackets() {
                num_brackets += 1;
                let (p1, p2) = bracket_pts(brk);
                let points = get_bracket_points(p1, p2, ref_pt, &sg_bond_segments, 0.1);
                self.post_shapes.push(DrawShape::poly_line(
                    points,
                    self.opts.bond_line_width,
                    false,
                    DrawColour::new(0.0, 0.0, 0.0),
                    false,
                    (-1, -1, -1),
                    no_dash(),
                )?);
            }
            if !self.include_annotations {
                continue;
            }
            // The bottom-most or right-most bracket: first work out whether
            // the bracket is largely horizontal or vertical.
            let last = self.post_shapes.len() - 1;
            let mut longline = self.post_shapes[last].points[1] - self.post_shapes[last].points[2];
            if longline.normalize().is_err() {
                // the bracket had no length so can be ignored
                continue;
            }
            let cos45 = 1.0 / 2.0f64.sqrt();
            let horizontal = longline.x.abs() > cos45;
            let mut label_brk = last;
            for i in 1..num_brackets {
                let idx = self.post_shapes.len() - i - 1;
                let (cand, cur) = (
                    &self.post_shapes[idx].points[2],
                    &self.post_shapes[label_brk].points[2],
                );
                if (horizontal && cand.y > cur.y) || (!horizontal && cand.x > cur.x) {
                    label_brk = idx;
                }
            }
            if let Some(connect) = sg.get_prop("CONNECT").map(|v| v.to_rd_string()) {
                // annotations go on the last bracket of an S-group; CONNECT goes
                // at the top, which is now the bottom due to the y inversion
                let p = &self.post_shapes[label_brk].points;
                let (mut bot_pt, mut brk_pt) = (p[2], p[3]);
                if (!horizontal && p[1].y < bot_pt.y) || (horizontal && p[1].x > bot_pt.x) {
                    bot_pt = p[1];
                    brk_pt = p[0];
                }
                let mut da = DrawAnnotation::new(
                    &connect,
                    TextAlignType::Middle,
                    "connect",
                    self.opts.annotation_font_scale,
                    bot_pt + (bot_pt - brk_pt),
                    DrawColour::new(0.0, 0.0, 0.0),
                    text,
                )?;
                // to the right of the bracket: left justify, otherwise things
                // work as they are
                if brk_pt.x < bot_pt.x {
                    da.align = TextAlignType::Start;
                }
                self.annotations.push(da);
            }
            let label = sg.get_prop("LABEL").map(|v| v.to_rd_string());
            let label = match label {
                Some(l) => Some(l),
                None => match sg.get_prop("TYPE").map(|v| v.to_rd_string()) {
                    // ChemDraw doesn't draw the GEN (generic) label
                    Some(t) if t == "GEN" => None,
                    // the lowercase type if there's no label to go there
                    Some(t) => Some(t.to_lowercase()),
                    None => None,
                },
            };
            if let Some(label) = label {
                let da = draw_bottom_label(
                    &label,
                    &self.post_shapes[label_brk],
                    self.opts.annotation_font_scale,
                    text,
                    horizontal,
                )?;
                self.annotations.push(da);
            }
        }
        Ok(())
    }

    pub fn extract_link_nodes(&mut self, text: &mut DrawTextFT) -> Result<()> {
        if !self.draw_mol.has_prop(cp::molFileLinkNodes) {
            return Ok(());
        }
        let link_nodes = get_mol_link_nodes(&self.draw_mol, false)?;
        for node in &link_nodes {
            let crossing_frac = 0.333;
            let length_frac = 0.333;
            let mut label_pt = Point2D::new(-1000.0, -1000.0);
            let mut label_perp = Point2D::new(0.0, 0.0);
            for &(first, second) in &node.bond_atoms {
                // unlike brackets, we know how these point
                let start_loc = self.at_cds[first as usize];
                let end_loc = self.at_cds[second as usize];
                let vect = end_loc - start_loc;
                let crossing_pt = start_loc + mul(vect, crossing_frac);
                let perp = mul(Point2D::new(vect.y, -vect.x), length_frac);
                let p1 = crossing_pt + mul(perp, 0.5);
                let p2 = crossing_pt - mul(perp, 0.5);
                let points = get_bracket_points(p1, p2, start_loc, &[], 0.1);
                self.post_shapes.push(DrawShape::poly_line(
                    points,
                    self.opts.bond_line_width,
                    false,
                    DrawColour::new(0.0, 0.0, 0.0),
                    false,
                    (-1, -1, -1),
                    no_dash(),
                )?);
                if p1.x > label_pt.x {
                    label_pt = p1;
                    label_perp = crossing_pt - start_loc;
                }
                if p2.x > label_pt.x {
                    label_pt = p2;
                    label_perp = crossing_pt - start_loc;
                }
            }
            if self.include_annotations {
                let label = format!("({}-{})", node.min_rep, node.max_rep);
                let perp = mul(label_perp, 1.0 / (label_perp.length() * 5.0));
                self.annotations.push(DrawAnnotation::new(
                    &label,
                    TextAlignType::Start,
                    "linknode",
                    self.opts.annotation_font_scale,
                    label_pt + perp,
                    DrawColour::new(0.0, 0.0, 0.0),
                    text,
                )?);
            }
        }
        Ok(())
    }

    /// Red boxes where two atoms are drawn too close together. This needs the
    /// final drawing coords, so it is called at the end of
    /// [`DrawMol::finish_create_draw_objects`], and the boxes go in the post
    /// shapes in drawing coordinates.
    pub fn extract_close_contacts(&mut self) -> Result<()> {
        if self.opts.flag_close_contacts_dist < 0 {
            return Ok(());
        }
        let tol =
            f64::from(self.opts.flag_close_contacts_dist * self.opts.flag_close_contacts_dist);
        let n = self.at_cds.len();
        let mut flagged = vec![false; n];
        let (trans, scale, to_centre) = self.draw_transformers();
        let tp =
            |p: Point2D| DrawMol::transform_point(p, Some(trans), Some(scale), Some(to_centre));
        for i in 0..n {
            if flagged[i] {
                continue;
            }
            let ci = tp(self.at_cds[i]);
            for j in i + 1..n {
                if flagged[j] {
                    continue;
                }
                let cj = tp(self.at_cds[j]);
                if (cj - ci).length_sq() <= tol {
                    flagged[i] = true;
                    flagged[j] = true;
                    break;
                }
            }
            if flagged[i] {
                let offset = Point2D::new(0.1 * self.scale, 0.1 * self.scale);
                let p1 = ci - offset;
                let p2 = ci + offset;
                let points = vec![
                    p1,
                    Point2D::new(p1.x, p2.y),
                    p2,
                    Point2D::new(p2.x, p1.y),
                    p1,
                ];
                self.post_shapes.push(DrawShape::poly_line(
                    points,
                    self.opts.bond_line_width,
                    false,
                    DrawColour::new(1.0, 0.0, 0.0),
                    false,
                    (i as i32, -1, -1),
                    no_dash(),
                )?);
            }
        }
        Ok(())
    }

    /// Moves multi-character atom labels that overlap another label to a
    /// different side of their atom.
    pub fn resolve_atom_symbol_clashes(&mut self, text: &DrawTextFT) -> Result<()> {
        let n = self.draw_mol.num_atoms() as usize;
        for i1 in 0..n {
            for i2 in i1 + 1..n {
                let (Some(l1), Some(l2)) = (&self.atom_labels[i1], &self.atom_labels[i2]) else {
                    continue;
                };
                let (n1, n2) = (l1.rects().len(), l2.rects().len());
                if !(n1 > 1 || n2 > 1) || !do_labels_clash(l1, l2) {
                    continue;
                }
                // prefer moving the smaller label that isn't a single character
                let mut idxs: [i32; 2] = [-1, -1];
                if n1 > 1 {
                    idxs[0] = i1 as i32;
                }
                if n2 > 1 {
                    idxs[1] = i2 as i32;
                }
                if n1 > n2 {
                    idxs.swap(0, 1);
                }
                let first_ok = idxs[0] != -1 && self.orient_atom_label(idxs[0] as usize, text)?;
                if !first_ok && idxs[1] != -1 {
                    self.orient_atom_label(idxs[1] as usize, text)?;
                }
            }
        }
        Ok(())
    }

    fn label_clashes_with_labels(&self, idx: usize) -> bool {
        let Some(this) = &self.atom_labels[idx] else {
            return false;
        };
        self.atom_labels
            .iter()
            .enumerate()
            .any(|(i, l)| i != idx && l.as_ref().is_some_and(|l| do_labels_clash(this, l)))
    }

    fn label_clashes_with_bonds(&self, idx: usize) -> bool {
        let Some(sym) = &self.atom_labels[idx] else {
            return false;
        };
        sym.rects().iter().any(|r| {
            let mut r = *r;
            r.trans += sym.cds;
            self.bonds.iter().any(|b| b.does_rect_clash(&r, 0.0))
        })
    }

    /// Tries the other orientations of a label, preferring the opposite of
    /// what it was; leaves it alone if none is clear.
    fn orient_atom_label(&mut self, idx: usize, text: &DrawTextFT) -> Result<bool> {
        const NEW_ORIENTS: [[OrientType; 3]; 4] = [
            [OrientType::S, OrientType::E, OrientType::W],
            [OrientType::N, OrientType::W, OrientType::E],
            [OrientType::E, OrientType::S, OrientType::N],
            [OrientType::W, OrientType::N, OrientType::S],
        ];
        let sym = self.atom_labels[idx].as_ref().expect("a label to orient");
        // A 2 character label whose second character is lower case is left
        // alone. (This is only called if there's more than one character.)
        if sym.rects().len() == 2 && sym.layout.draw_chars[1].is_ascii_lowercase() {
            return Ok(false);
        }
        let orig = sym.orient;
        let pref = match orig {
            OrientType::S => 1,
            OrientType::W => 2,
            OrientType::E => 3,
            _ => 0,
        };
        for orient in NEW_ORIENTS[pref] {
            let sym = self.atom_labels[idx].as_mut().expect("a label to orient");
            sym.orient = orient;
            sym.recalculate_rects(text)?;
            if !self.label_clashes_with_labels(idx) && !self.label_clashes_with_bonds(idx) {
                return Ok(true);
            }
        }
        // leave it as it was and hope a better solution is found when the other
        // label is moved
        self.atom_labels[idx]
            .as_mut()
            .expect("a label to orient")
            .orient = orig;
        Ok(false)
    }
}

/// Do any rectangles of the first label clash with the second label?
fn do_labels_clash(a: &AtomSymbol, b: &AtomSymbol) -> bool {
    a.rects().iter().any(|r| {
        let mut r = *r;
        r.trans += a.cds;
        b.does_rect_clash(&r, 0.0)
    })
}
