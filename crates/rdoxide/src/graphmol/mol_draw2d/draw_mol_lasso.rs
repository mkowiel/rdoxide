//! Multi-colour highlights drawn as lassos round the atoms of each colour
//! (`DrawMolMCHLasso.cpp` of RDKit 2026.09.1).
//!
//! Copyright (C) 2023 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::details::{adjust_line_end_for_ellipse, do_lines_intersect};
use super::draw_mol::{DrawMol, McHighlights};
use super::draw_mol_utils::calc_perpendicular;
use super::draw_shape::{DrawShape, ShapeKind, no_dash};
use super::helpers::DrawColour;
use crate::Result;
use crate::geometry::Point2D;

const LINE_WIDTH: f64 = 3.0;
const SCALE_LINE_WIDTH: bool = true;

/// The two lines either side of a bond at one of its atoms, as indices into
/// the lines of the lasso, and the angles at which they leave the atom's arc.
#[derive(Clone, Copy, Debug, Default)]
struct LinePair {
    line1: usize,
    line2: usize,
    radius: f64,
    atom: i32,
    angle1: f64,
    angle2: f64,
}

fn arc_angles(shape: &DrawShape) -> (f64, f64) {
    match shape.kind {
        ShapeKind::Arc { ang1, ang2 } => (ang1, ang2),
        _ => (0.0, 0.0),
    }
}

/// The angles at which lines to `pt1` and `pt2` leave the arc round `at1`, in
/// degrees, ordered anti-clockwise round the bond to `at2`; the angle of the
/// bond, and whether the points were swapped.
fn calc_subtended_angles(
    pt1: Point2D,
    pt2: Point2D,
    at1: Point2D,
    at2: Point2D,
) -> Result<(f64, f64, f64, bool)> {
    let index = Point2D::new(1.0, 0.0);
    let rad1 = at1.direction_vector(&pt1)?;
    let rad2 = at1.direction_vector(&pt2)?;
    let brad = at1.direction_vector(&at2)?;

    let wrap = |a: f64| if a >= 360.0 { a - 360.0 } else { a };
    let mut ang1 = wrap(360.0 - rad1.signed_angle_to(&index)? * 180.0 / PI);
    let mut ang2 = wrap(360.0 - rad2.signed_angle_to(&index)? * 180.0 / PI);
    let bang = wrap(360.0 - brad.signed_angle_to(&index)? * 180.0 / PI);
    let cross = rad1.x * rad2.y - rad1.y * rad2.x;
    let mut swapped = false;
    if cross > 0.0 {
        std::mem::swap(&mut ang1, &mut ang2);
        swapped = true;
    }
    let min_ang = ang1.min(ang2).min(bang);
    let (rang1, rang2);
    if ang1 - min_ang < bang - min_ang && bang - min_ang < ang2 - min_ang {
        rang1 = ang1;
        rang2 = ang2;
    } else {
        rang1 = ang2;
        rang2 = ang1;
        swapped = !swapped;
    }
    Ok((rang1, rang2, bang, swapped))
}

/// The points where two circles of the same size intersect.
fn intersecting_circle_points_same_size(
    at1: Point2D,
    at2: Point2D,
    rad: f64,
) -> (Point2D, Point2D) {
    // the intersection points are where the perpendicular at the mid-point
    // between the circles intersects with either circle
    let mid = (at1 + at2) * 0.5;
    let perp = calc_perpendicular(at1, at2);
    let mut isect1 = mid + perp * (2.0 * rad);
    adjust_line_end_for_ellipse(at1, rad, rad, mid, &mut isect1);
    let isect2 = mid - perp * (mid - isect1).length();
    (isect1, isect2)
}

/// The points where two circles of different sizes intersect (after
/// ambrnet.com/TrigoCalc/Circles2/circle2intersection).
fn intersecting_circle_points_diff_size(
    at1: Point2D,
    rad1: f64,
    at2: Point2D,
    rad2: f64,
) -> (Point2D, Point2D) {
    let dist = (at1 - at2).length();
    let dist2 = dist * dist;
    // We already know the circles intersect. Triangle area via Heron's formula
    let a1 = dist + rad1 + rad2;
    let a2 = dist + rad1 - rad2;
    let a3 = dist - rad1 + rad2;
    let a4 = rad1 + rad2 - dist;
    let area = (a1 * a2 * a3 * a4).sqrt() / 4.0;
    let c = (rad1 * rad1 - rad2 * rad2) / (2.0 * dist2);
    let mut isect1 = Point2D::new(0.0, 0.0);
    let mut isect2 = Point2D::new(0.0, 0.0);
    let val1 = (at1.x + at2.x) / 2.0 + (at2.x - at1.x) * c;
    let val2 = 2.0 * (at1.y - at2.y) * area / dist2;
    isect1.x = val1 + val2;
    isect2.x = val1 - val2;
    let val1 = (at1.y + at2.y) / 2.0 + (at2.y - at1.y) * c;
    let val2 = 2.0 * (at1.x - at2.x) * area / dist2;
    isect1.y = val1 + val2;
    isect2.y = val1 - val2;
    // check that the 2 combinations of y for x taken are the correct 2 by
    // making sure isect1 is on the circle
    if ((isect1 - at1).length() - rad1 * rad1).abs() > 1.0e-8 {
        std::mem::swap(&mut isect1.y, &mut isect2.y);
    }
    (isect1, isect2)
}

/// The angles of the arc of the circle round `at1` between the points where
/// it meets the circle round `at2`, and the angle of the line between them.
fn calc_intersecting_arc_angles(
    at1: Point2D,
    rad1: f64,
    at2: Point2D,
    rad2: f64,
) -> Result<(f64, f64, f64)> {
    let (isect1, isect2) = if (rad1 - rad2).abs() < 1.0e-8 {
        intersecting_circle_points_same_size(at1, at2, rad1)
    } else {
        intersecting_circle_points_diff_size(at1, rad1, at2, rad2)
    };
    let (ang1, ang2, bang, _) = calc_subtended_angles(isect1, isect2, at1, at2)?;
    Ok((ang1, ang2, bang))
}

impl DrawMol {
    /// The radius of the lasso of the given number round an atom: the areas
    /// double for successive lassos.
    fn lasso_width(&self, at_num: u32, lasso_num: usize) -> f64 {
        const RATS: [f64; 5] = [1.0, 1.414, 2.0, 2.828, 4.0];
        let (xrad, _) = self.mch_atom_radius(at_num);
        if lasso_num > 4 {
            // It's going to look horrible, probably, but it's a lot of lassos.
            xrad * (1 + lasso_num) as f64 * 0.75
        } else {
            xrad * RATS[lasso_num]
        }
    }

    /// `extractMCHighlights` of `DrawMolMCHLasso`.
    pub fn extract_mc_lasso(&mut self, mch: &McHighlights) -> Result<()> {
        let (colours, colour_atoms, colour_lists) = self.extract_atom_colour_lists(mch);
        for list in &colour_lists {
            for (i, &c) in list.iter().enumerate() {
                self.draw_lasso(mch, i, colours[c], &colour_atoms[c])?;
            }
        }
        Ok(())
    }

    /// The atoms of each colour, and the lists of colours whose atoms
    /// overlap, each of which gets nested lassos.
    fn extract_atom_colour_lists(
        &self,
        mch: &McHighlights,
    ) -> (Vec<DrawColour>, Vec<Vec<i32>>, Vec<Vec<usize>>) {
        let n_atoms = self.draw_mol.num_atoms() as usize;
        let mut colours: Vec<DrawColour> = Vec::new();
        let mut colour_atoms: Vec<Vec<i32>> = Vec::new();
        let mut in_colour_atoms: Vec<Vec<bool>> = Vec::new();
        for (&atom, cols) in &mch.atom_map {
            for col in cols {
                match colours.iter().position(|c| c == col) {
                    None => {
                        colours.push(*col);
                        colour_atoms.push(vec![atom]);
                        let mut set = vec![false; n_atoms];
                        if let Some(s) = set.get_mut(atom as usize) {
                            *s = true;
                        }
                        in_colour_atoms.push(set);
                    }
                    Some(ln) => {
                        // each atom must be in the list only once (Github6749)
                        let present = in_colour_atoms[ln]
                            .get(atom as usize)
                            .copied()
                            .unwrap_or(false);
                        if !present {
                            colour_atoms[ln].push(atom);
                            if let Some(s) = in_colour_atoms[ln].get_mut(atom as usize) {
                                *s = true;
                            }
                        }
                    }
                }
            }
        }
        let mut colour_lists: Vec<Vec<usize>> = (0..colour_atoms.len()).map(|i| vec![i]).collect();
        if colour_lists.len() < 2 {
            return (colours, colour_atoms, colour_lists);
        }
        let lists_intersect = |v1: &[usize], v2: &[usize]| -> bool {
            v1.iter().any(|&m1| {
                colour_atoms[m1]
                    .iter()
                    .any(|ca| v2.iter().any(|&m2| colour_atoms[m2].contains(ca)))
            })
        };
        let mut did_something = true;
        while did_something && colour_lists.len() > 1 {
            did_something = false;
            'outer: for i in 0..colour_lists.len() - 1 {
                for j in i + 1..colour_lists.len() {
                    if lists_intersect(&colour_lists[i], &colour_lists[j]) {
                        let moved = std::mem::take(&mut colour_lists[j]);
                        colour_lists[i].extend(moved);
                        did_something = true;
                        break 'outer;
                    }
                }
            }
            colour_lists.retain(|v| !v.is_empty());
        }
        (colours, colour_atoms, colour_lists)
    }

    /// Extracts the arcs and lines for the given atoms in the given colour;
    /// `lasso_num` dictates the radii of the arcs and the displacements of
    /// the lines.
    fn draw_lasso(
        &mut self,
        mch: &McHighlights,
        lasso_num: usize,
        col: DrawColour,
        col_atoms: &[i32],
    ) -> Result<()> {
        if col_atoms.is_empty() {
            return Ok(());
        }
        let n_atoms = self.draw_mol.num_atoms() as usize;
        let mut lines: Vec<DrawShape> = Vec::new();
        let mut atom_lines: Vec<Vec<LinePair>> = vec![Vec::new(); n_atoms];
        self.extract_bond_lines(mch, lasso_num, col, col_atoms, &mut lines, &mut atom_lines)?;
        let mut arcs = self.extract_atom_arcs(&atom_lines, &mut lines)?;
        self.add_no_line_arcs(col_atoms, lasso_num, col, &lines, &mut arcs)?;
        self.highlights.extend(arcs);
        self.highlights.extend(lines);
        Ok(())
    }

    fn extract_bond_lines(
        &self,
        mch: &McHighlights,
        lasso_num: usize,
        col: DrawColour,
        col_atoms: &[i32],
        lines: &mut Vec<DrawShape>,
        atom_lines: &mut [Vec<LinePair>],
    ) -> Result<()> {
        let n_atoms = self.draw_mol.num_atoms() as i32;
        if col_atoms.len() > 1 {
            for i in 0..col_atoms.len() - 1 {
                if col_atoms[i] < 0 || col_atoms[i] >= n_atoms {
                    // there's an error in the colour map
                    continue;
                }
                let lasso_width_i = self.lasso_width(col_atoms[i] as u32, lasso_num);
                let disp_i = lasso_width_i * 0.75;
                for j in i + 1..col_atoms.len() {
                    if col_atoms[j] < 0 || col_atoms[j] >= n_atoms {
                        continue;
                    }
                    let lasso_width_j = self.lasso_width(col_atoms[j] as u32, lasso_num);
                    let disp_j = lasso_width_j * 0.75;
                    let Some(bond) = self
                        .draw_mol
                        .get_bond_between_atoms(col_atoms[i] as u32, col_atoms[j] as u32)?
                    else {
                        continue;
                    };
                    let bond_idx = bond.idx() as i32;
                    if !mch.bond_map.is_empty() && !mch.bond_map.contains_key(&bond_idx) {
                        continue;
                    }
                    let at_cds_i = self.at_cds[col_atoms[i] as usize];
                    let at_cds_j = self.at_cds[col_atoms[j] as usize];
                    let perp = calc_perpendicular(at_cds_i, at_cds_j);
                    let mut pair = LinePair::default();
                    let mut skip = false;
                    for m in [-1.0f64, 1.0] {
                        let mut p1 = at_cds_i + perp * (disp_i * m);
                        let mut p2 = at_cds_j + perp * (disp_j * m);
                        // If the mid point of the line is inside one of the
                        // circles it gets messy because of the intersecting
                        // arcs, so skip it.
                        let mid = (p1 + p2) * 0.5;
                        if (at_cds_i - mid).length_sq() < lasso_width_i * lasso_width_i
                            || (at_cds_i - mid).length_sq() < lasso_width_j * lasso_width_j
                        {
                            skip = true;
                            continue;
                        }
                        adjust_line_end_for_ellipse(
                            at_cds_i,
                            lasso_width_i,
                            lasso_width_i,
                            p2,
                            &mut p1,
                        );
                        adjust_line_end_for_ellipse(
                            at_cds_j,
                            lasso_width_j,
                            lasso_width_j,
                            p1,
                            &mut p2,
                        );
                        lines.push(DrawShape::simple_line(
                            vec![p1, p2],
                            LINE_WIDTH,
                            SCALE_LINE_WIDTH,
                            col,
                            (col_atoms[i], col_atoms[j], bond_idx),
                            no_dash(),
                        )?);
                        if m < 0.0 {
                            pair.line1 = lines.len() - 1;
                        } else {
                            pair.line2 = lines.len() - 1;
                        }
                    }
                    if !skip {
                        pair.radius = lasso_width_i;
                        pair.atom = col_atoms[i];
                        atom_lines[col_atoms[i] as usize].push(pair);
                        pair.radius = lasso_width_j;
                        pair.atom = col_atoms[j];
                        atom_lines[col_atoms[j] as usize].push(pair);
                    }
                }
            }
        }
        self.order_atom_lines(lines, atom_lines)
    }

    /// Puts the start points of each pair of lines at an atom in order so that
    /// going anti-clockwise the bond the lines straddle is between them, and
    /// the arc between the first and second is the bit that isn't drawn; then
    /// orders the pairs so that the bonds move round anti-clockwise.
    fn order_atom_lines(
        &self,
        lines: &[DrawShape],
        atom_lines: &mut [Vec<LinePair>],
    ) -> Result<()> {
        for i in 0..self.draw_mol.num_atoms() as usize {
            if atom_lines[i].is_empty() {
                continue;
            }
            let mut bond_angles: Vec<(f64, usize)> = Vec::new();
            let mut min_bond_angle = 720.0f64;
            for j in 0..atom_lines[i].len() {
                let lp = atom_lines[i][j];
                let l1 = &lines[lp.line1];
                let l2 = &lines[lp.line2];
                // because the same line is used for both atoms, points[0] may
                // not be the nearer point to atom i.
                let oatom = if l1.atom1 == i as i32 {
                    l1.atom2
                } else {
                    l1.atom1
                };
                let pt = if (self.at_cds[i] - l1.points[0]).length_sq()
                    > (self.at_cds[i] - l1.points[1]).length_sq()
                {
                    1
                } else {
                    0
                };
                let (a1, a2, bang, swapped) = calc_subtended_angles(
                    l1.points[pt],
                    l2.points[pt],
                    self.at_cds[i],
                    self.at_cds[oatom as usize],
                )?;
                atom_lines[i][j].angle1 = a1;
                atom_lines[i][j].angle2 = a2;
                bond_angles.push((bang, j));
                min_bond_angle = min_bond_angle.min(bang);
                if swapped {
                    let lp = &mut atom_lines[i][j];
                    std::mem::swap(&mut lp.line1, &mut lp.line2);
                }
            }
            for ba in &mut bond_angles {
                ba.0 -= min_bond_angle;
            }
            bond_angles.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            let reordered: Vec<LinePair> =
                bond_angles.iter().map(|ba| atom_lines[i][ba.1]).collect();
            atom_lines[i] = reordered;
        }
        Ok(())
    }

    /// An arc between line 2 of `lp1` and line 1 of `lp2`. If the two lines
    /// intersect outside the radius of the arc, both are trimmed back to the
    /// point of intersection and there is no arc.
    fn make_lasso_arc(
        lp1: &LinePair,
        lp2: &LinePair,
        at_cds: Point2D,
        lines: &mut [DrawShape],
    ) -> Result<Option<DrawShape>> {
        let adjust_line = |line: &mut DrawShape, pt: Point2D| {
            let d1 = (line.points[0] - pt).length_sq();
            let d2 = (line.points[1] - pt).length_sq();
            if d1 < d2 {
                line.points[0] = pt;
            } else {
                line.points[1] = pt;
            }
        };
        let (l2a, l2b) = (lines[lp1.line2].points[0], lines[lp1.line2].points[1]);
        let (l1a, l1b) = (lines[lp2.line1].points[0], lines[lp2.line1].points[1]);
        if let Some(ip) = do_lines_intersect(l2a, l2b, l1a, l1b) {
            adjust_line(&mut lines[lp1.line2], ip);
            adjust_line(&mut lines[lp2.line1], ip);
            return Ok(None);
        }
        let l = &lines[lp1.line1];
        Ok(Some(DrawShape::arc(
            vec![at_cds, Point2D::new(lp1.radius, lp1.radius)],
            lp1.angle2,
            lp2.angle1,
            l.line_width,
            l.scale_line_width,
            l.line_colour,
            false,
            lp1.atom,
        )?))
    }

    fn extract_atom_arcs(
        &self,
        atom_lines: &[Vec<LinePair>],
        lines: &mut [DrawShape],
    ) -> Result<Vec<DrawShape>> {
        let mut arcs = Vec::new();
        for atom_line in atom_lines {
            if atom_line.is_empty() {
                continue;
            }
            if atom_line.len() == 1 {
                let lp = &atom_line[0];
                let l = &lines[lp.line1];
                arcs.push(DrawShape::arc(
                    vec![
                        self.at_cds[lp.atom as usize],
                        Point2D::new(lp.radius, lp.radius),
                    ],
                    lp.angle2,
                    lp.angle1,
                    l.line_width,
                    l.scale_line_width,
                    l.line_colour,
                    false,
                    lp.atom,
                )?);
            } else {
                for i in 0..atom_line.len() - 1 {
                    let at = self.at_cds[atom_line[i].atom as usize];
                    if let Some(arc) =
                        Self::make_lasso_arc(&atom_line[i], &atom_line[i + 1], at, lines)?
                    {
                        arcs.push(arc);
                    }
                }
                let last = atom_line.last().expect("not empty");
                let first = &atom_line[0];
                let at = self.at_cds[first.atom as usize];
                if let Some(arc) = Self::make_lasso_arc(last, first, at, lines)? {
                    arcs.push(arc);
                }
            }
        }
        Ok(arcs)
    }

    /// Circles for the atoms with no line to another atom of the colour, and
    /// arcs where the circles of atoms without a line between them overlap.
    fn add_no_line_arcs(
        &self,
        col_atoms: &[i32],
        lasso_num: usize,
        col: DrawColour,
        lines: &[DrawShape],
        arcs: &mut Vec<DrawShape>,
    ) -> Result<()> {
        let n = self.draw_mol.num_atoms() as usize;
        let mut in_col_atoms = vec![false; n];
        for &ca in col_atoms {
            if let Some(s) = in_col_atoms.get_mut(ca as usize) {
                *s = true;
            }
        }
        let mut in_lines = vec![vec![false; n]; n];
        for l in lines {
            in_lines[l.atom2 as usize][l.atom1 as usize] = true;
            in_lines[l.atom1 as usize][l.atom2 as usize] = true;
        }
        let mut intersects: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            if in_col_atoms[i] && !in_lines[i].iter().any(|&b| b) {
                intersects[i].push(i);
            }
        }
        for i in 0..n.saturating_sub(1) {
            if !in_col_atoms[i] {
                continue;
            }
            let rad_i = self.lasso_width(i as u32, lasso_num);
            for j in i + 1..n {
                if !in_col_atoms[j] {
                    continue;
                }
                let rad_j = self.lasso_width(j as u32, lasso_num);
                if !in_lines[i][j] && (self.at_cds[i] - self.at_cds[j]).length() < rad_i + rad_j {
                    if intersects[i].is_empty() {
                        intersects[i].push(i);
                    }
                    intersects[i].push(j);
                    if intersects[j].is_empty() {
                        intersects[j].push(j);
                    }
                    intersects[j].push(i);
                }
            }
        }
        let new_arcs = self.make_intersecting_arcs(&intersects, lasso_num, col, arcs)?;
        arcs.extend(new_arcs);
        Ok(())
    }

    fn are_arc_ends_in_other_circle(
        &self,
        centre: Point2D,
        radius: f64,
        ang1: f64,
        ang2: f64,
        at_num: usize,
        other_atoms: &[usize],
        lasso_num: usize,
    ) -> bool {
        let a1 = ang1 * PI / 180.0;
        let end1 = Point2D::new(centre.x + radius * a1.cos(), centre.y + radius * a1.sin());
        let a2 = ang2 * PI / 180.0;
        let end2 = Point2D::new(centre.x + radius * a2.cos(), centre.y + radius * a2.sin());
        for &other in other_atoms {
            if other == at_num {
                continue;
            }
            let other_rad = self.lasso_width(other as u32, lasso_num);
            for end in [end1, end2] {
                if (end - self.at_cds[other]).length_sq() - other_rad * other_rad < -1.0e-4 {
                    return true;
                }
            }
        }
        false
    }

    fn make_intersecting_arcs(
        &self,
        intersects: &[Vec<usize>],
        lasso_num: usize,
        col: DrawColour,
        curr_arcs: &mut Vec<DrawShape>,
    ) -> Result<Vec<DrawShape>> {
        let mut arcs = Vec::new();
        for (i, inter) in intersects.iter().enumerate() {
            if inter.is_empty() {
                continue;
            }
            let rad_i = self.lasso_width(i as u32, lasso_num);
            let pts = vec![self.at_cds[i], Point2D::new(rad_i, rad_i)];
            if inter.len() == 1 && inter[0] == i {
                arcs.push(DrawShape::arc(
                    pts,
                    0.0,
                    360.0,
                    LINE_WIDTH,
                    SCALE_LINE_WIDTH,
                    col,
                    false,
                    i as i32,
                )?);
                continue;
            }
            let mut arc_angs: Vec<f64> = Vec::new();
            let mut bangs: Vec<(f64, usize)> = Vec::new();
            // add the arcs already made for this atom, removing them
            curr_arcs.retain(|arc| {
                if arc.atom1 != i as i32 {
                    return true;
                }
                let (ang1, ang2) = arc_angles(arc);
                arc_angs.push(ang2);
                arc_angs.push(ang1);
                let bang = if ang2 < ang1 {
                    (ang2 + ang1) / 2.0
                } else {
                    -360.0 + (ang2 + ang1 + 360.0) / 2.0
                };
                bangs.push((bang, bangs.len()));
                false
            });
            // the first entry is the atom itself, that's how it's set up
            for &j in &inter[1..] {
                let rad_j = self.lasso_width(j as u32, lasso_num);
                let (ang1, ang2, bang) =
                    calc_intersecting_arc_angles(self.at_cds[i], rad_i, self.at_cds[j], rad_j)?;
                arc_angs.push(ang1);
                arc_angs.push(ang2);
                bangs.push((bang, bangs.len()));
            }
            bangs.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            let sorted: Vec<f64> = bangs
                .iter()
                .flat_map(|b| [arc_angs[2 * b.1], arc_angs[2 * b.1 + 1]])
                .collect();
            // Now make the arcs from i.second to (i + 1).first, and finally
            // from last.first to first.second
            let mut j = 0;
            while j + 2 < sorted.len() {
                // if either end of the arc is inside another of the
                // intersecting rings, skip it
                if !self.are_arc_ends_in_other_circle(
                    self.at_cds[i],
                    rad_i,
                    sorted[j + 1],
                    sorted[j + 2],
                    i,
                    inter,
                    lasso_num,
                ) {
                    arcs.push(DrawShape::arc(
                        pts.clone(),
                        sorted[j + 1],
                        sorted[j + 2],
                        LINE_WIDTH,
                        SCALE_LINE_WIDTH,
                        col,
                        false,
                        i as i32,
                    )?);
                }
                j += 2;
            }
            let (last, first) = (sorted[sorted.len() - 1], sorted[0]);
            if !self.are_arc_ends_in_other_circle(
                self.at_cds[i],
                rad_i,
                last,
                first,
                i,
                inter,
                lasso_num,
            ) {
                arcs.push(DrawShape::arc(
                    pts,
                    last,
                    first,
                    LINE_WIDTH,
                    SCALE_LINE_WIDTH,
                    col,
                    false,
                    i as i32,
                )?);
            }
        }
        Ok(arcs)
    }
}
