//! The object that holds everything drawn for one molecule (`DrawMol.cpp` of
//! RDKit 2026.09.1): atom coordinates, atom labels, bonds, highlights,
//! annotations and the legend, with the scale and transformations that take
//! them from molecule coordinates to drawing coordinates.
//!
//! Upstream keeps references to the owning drawer's options and text drawer;
//! here the options are copied in and the text drawer is passed to the
//! methods that use it. The `extract*` steps are added piece by piece: see
//! [`DrawMol::extract_all`] for the ones present.
//!
//! Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::f64::consts::PI;

use super::draw_annotation::{AtomSymbol, DrawAnnotation};
use super::draw_mol_utils::{
    center_mol_for_drawing, find_radical_extremes, get_atom_symbol_and_orientation,
    get_colour_by_atomic_num, get_highlight_bond_colour,
};
use super::draw_shape::DrawShape;
use super::draw_text::DrawTextFT;
use super::helpers::{DrawColour, LegendPosition, MolDrawOptions, OrientType, draw_element};
use super::string_rect::StringRect;
use super::utils::{mean_bond_length, prepare_mol_for_drawing};
use crate::geometry::{Point2D, Transform2D};
use crate::graphmol::chirality::stereo_groups::simplify_enhanced_stereo;
use crate::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::wedge_bonds::reapply_mol_block_wedging;
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result};

/// The highlighting and options a [`DrawMol`] is created with.
#[derive(Clone, Debug, Default)]
pub struct DrawMolHighlights {
    pub atoms: Vec<i32>,
    pub bonds: Vec<i32>,
    pub atom_map: BTreeMap<i32, DrawColour>,
    pub bond_map: BTreeMap<i32, DrawColour>,
    pub bond_colours: Vec<(DrawColour, DrawColour)>,
    pub radii: BTreeMap<i32, f64>,
}

/// Everything needed to draw one molecule.
#[derive(Debug)]
pub struct DrawMol {
    pub opts: MolDrawOptions,
    /// For reactions padding needs to be 0 whatever the options say, so it is
    /// copied on construction.
    pub margin_padding: f64,
    pub highlight_atoms: Vec<i32>,
    pub highlight_bonds: Vec<i32>,
    pub highlight_atom_map: BTreeMap<i32, DrawColour>,
    pub highlight_bond_map: BTreeMap<i32, DrawColour>,
    pub bond_colours: Vec<(DrawColour, DrawColour)>,
    pub highlight_radii: BTreeMap<i32, f64>,
    pub include_annotations: bool,
    pub is_reaction_mol: bool,
    pub legend: String,

    pub draw_mol: RWMol,
    pub conf_id: i32,
    /// Atom coordinates as in the molecule, but with y inverted and
    /// `opts.rotate` applied.
    pub at_cds: Vec<Point2D>,
    pub bonds: Vec<DrawShape>,
    pub pre_shapes: Vec<DrawShape>,
    pub post_shapes: Vec<DrawShape>,
    pub atomic_nums: Vec<i32>,
    pub atom_syms: Vec<(String, OrientType)>,
    pub atom_labels: Vec<Option<AtomSymbol>>,
    pub highlights: Vec<DrawShape>,
    pub annotations: Vec<DrawAnnotation>,
    pub legends: Vec<DrawAnnotation>,
    pub radicals: Vec<(StringRect, OrientType, i32)>,
    pub single_bond_lines: Vec<i32>,

    /// Total size of the canvas; -1 means it is worked out from the molecule
    /// (the "flexiCanvas").
    pub width: i32,
    pub height: i32,
    /// The drawing area within the canvas: the canvas less the padding.
    pub draw_width: i32,
    pub draw_height: i32,
    /// The font scale is independent of the scale to allow for min and max
    /// font sizes.
    pub scale: f64,
    pub font_scale: f64,
    pub x_min: f64,
    pub y_min: f64,
    pub x_max: f64,
    pub y_max: f64,
    pub x_range: f64,
    pub y_range: f64,
    /// For drawing molecules in grids.
    pub x_offset: f64,
    pub y_offset: f64,
    pub mean_bond_length: f64,
    /// With a legend, `mol_height`/`mol_width` are the parts for the molecule
    /// and `legend_height`/`legend_width` the parts for the legend, in pixels.
    pub mol_height: i32,
    pub legend_height: i32,
    pub mol_width: i32,
    pub legend_width: i32,
    pub drawing_initialised: bool,
    /// Atoms and bonds in an SVG are given a class from the drawer's active
    /// indices, which needn't start from 0.
    pub active_atm_idx_offset: i32,
    pub active_bnd_idx_offset: i32,
    pub flexi_canvas_x: bool,
    pub flexi_canvas_y: bool,
}

impl DrawMol {
    /// The extremes before anything has been found.
    fn initial_limits() -> (f64, f64, f64, f64, f64, f64) {
        (
            f64::MAX / 2.0,
            f64::MAX / 2.0,
            f64::MIN / 2.0,
            f64::MIN / 2.0,
            f64::MAX,
            f64::MAX,
        )
    }

    /// Makes the object for a molecule, scaled to the given pixel size (-1 for
    /// either to have the canvas size set automatically).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mol: &ROMol,
        legend: &str,
        width: i32,
        height: i32,
        opts: &MolDrawOptions,
        text: &mut DrawTextFT,
        highlights: DrawMolHighlights,
        include_annotations: bool,
        conf_id: i32,
        is_reaction_mol: bool,
    ) -> Result<Self> {
        let (x_min, y_min, x_max, y_max, x_range, y_range) = Self::initial_limits();
        let mut dm = DrawMol {
            opts: opts.clone(),
            margin_padding: opts.padding,
            highlight_atoms: highlights.atoms,
            highlight_bonds: highlights.bonds,
            highlight_atom_map: highlights.atom_map,
            highlight_bond_map: highlights.bond_map,
            bond_colours: highlights.bond_colours,
            highlight_radii: highlights.radii,
            include_annotations,
            is_reaction_mol,
            legend: legend.to_string(),
            draw_mol: RWMol::new(),
            conf_id,
            at_cds: Vec::new(),
            bonds: Vec::new(),
            pre_shapes: Vec::new(),
            post_shapes: Vec::new(),
            atomic_nums: Vec::new(),
            atom_syms: Vec::new(),
            atom_labels: Vec::new(),
            highlights: Vec::new(),
            annotations: Vec::new(),
            legends: Vec::new(),
            radicals: Vec::new(),
            single_bond_lines: Vec::new(),
            width,
            height,
            draw_width: width,
            draw_height: height,
            scale: 1.0,
            font_scale: 1.0,
            x_min,
            y_min,
            x_max,
            y_max,
            x_range,
            y_range,
            x_offset: 0.0,
            y_offset: 0.0,
            mean_bond_length: 0.0,
            mol_height: 0,
            legend_height: 0,
            mol_width: 0,
            legend_width: 0,
            drawing_initialised: false,
            active_atm_idx_offset: 0,
            active_bnd_idx_offset: 0,
            flexi_canvas_x: width < 0,
            flexi_canvas_y: height < 0,
        };
        text.set_font_scale(dm.font_scale, true);
        dm.init_draw_molecule(mol)?;
        Ok(dm)
    }

    /// Makes an object with nothing to draw, but with a scale and extremes
    /// for converting coordinates (used by the drawer for its own state).
    #[allow(clippy::too_many_arguments)]
    pub fn without_mol(
        width: i32,
        height: i32,
        opts: &MolDrawOptions,
        text: &mut DrawTextFT,
        xmin: f64,
        xmax: f64,
        ymin: f64,
        ymax: f64,
        scale: f64,
        font_scale: f64,
    ) -> Self {
        let mut dm = DrawMol::new_empty(width, height, opts);
        dm.scale = scale;
        dm.font_scale = font_scale;
        dm.x_min = xmin;
        dm.y_min = ymin;
        dm.x_max = xmax;
        dm.y_max = ymax;
        dm.x_range = xmax - xmin;
        dm.y_range = ymax - ymin;
        dm.mol_height = height;
        text.set_font_scale(dm.font_scale, true);
        // the y coords of everything are reversed, so do that here, too
        let (lo, hi) = (-dm.y_max, -dm.y_min);
        dm.y_min = lo;
        dm.y_max = hi;
        dm
    }

    fn new_empty(width: i32, height: i32, opts: &MolDrawOptions) -> Self {
        let (x_min, y_min, x_max, y_max, x_range, y_range) = Self::initial_limits();
        DrawMol {
            opts: opts.clone(),
            margin_padding: opts.padding,
            highlight_atoms: Vec::new(),
            highlight_bonds: Vec::new(),
            highlight_atom_map: BTreeMap::new(),
            highlight_bond_map: BTreeMap::new(),
            bond_colours: Vec::new(),
            highlight_radii: BTreeMap::new(),
            include_annotations: false,
            is_reaction_mol: false,
            legend: String::new(),
            draw_mol: RWMol::new(),
            conf_id: -1,
            at_cds: Vec::new(),
            bonds: Vec::new(),
            pre_shapes: Vec::new(),
            post_shapes: Vec::new(),
            atomic_nums: Vec::new(),
            atom_syms: Vec::new(),
            atom_labels: Vec::new(),
            highlights: Vec::new(),
            annotations: Vec::new(),
            legends: Vec::new(),
            radicals: Vec::new(),
            single_bond_lines: Vec::new(),
            width,
            height,
            draw_width: width,
            draw_height: height,
            scale: 1.0,
            font_scale: 1.0,
            x_min,
            y_min,
            x_max,
            y_max,
            x_range,
            y_range,
            x_offset: 0.0,
            y_offset: 0.0,
            mean_bond_length: 0.0,
            mol_height: 0,
            legend_height: 0,
            mol_width: 0,
            legend_width: 0,
            drawing_initialised: false,
            active_atm_idx_offset: 0,
            active_bnd_idx_offset: 0,
            flexi_canvas_x: false,
            flexi_canvas_y: false,
        }
    }

    /// Copies the molecule and prepares the copy as the options ask.
    pub fn init_draw_molecule(&mut self, mol: &ROMol) -> Result<()> {
        self.draw_mol = mol.clone();
        if self.opts.centre_molecules_before_drawing && self.draw_mol.num_conformers() > 0 {
            center_mol_for_drawing(&mut self.draw_mol, self.conf_id)?;
        }
        if !self.is_reaction_mol && self.draw_mol.num_conformers() == 0 {
            let params = Compute2DCoordParameters {
                canon_orient: true,
                ..Default::default()
            };
            compute_2d_coords(&mut self.draw_mol, &params)?;
        }
        if self.opts.unspecified_stereo_is_unknown {
            return Err(RDKitError::UnhandledFeature(
                "unspecifiedStereoIsUnknown is not supported yet".into(),
            ));
        }
        if self.opts.use_mol_block_wedging {
            reapply_mol_block_wedging(&mut self.draw_mol, true, false)?;
        }
        if !self.is_reaction_mol && self.opts.prepare_mols_before_drawing {
            prepare_mol_for_drawing(&mut self.draw_mol, true, true, true, false, false)?;
        }
        if self.opts.simplified_stereo_group_label && !mol.has_prop(cp::molNote) {
            simplify_enhanced_stereo(&mut self.draw_mol, true);
        }
        if self.opts.add_atom_indices {
            add_atom_indices(&mut self.draw_mol)?;
        }
        if self.opts.add_bond_indices {
            add_bond_indices(&mut self.draw_mol)?;
        }
        Ok(())
    }

    /// Extracts everything to be drawn. Present so far: atom coordinates, labels,
    /// bonds, regions, highlights and attachments.
    pub fn extract_all(&mut self, scale: f64, text: &mut DrawTextFT) -> Result<()> {
        self.extract_atom_coords()?;
        self.extract_atom_symbols(text)?;
        self.extract_bonds()?;
        self.extract_regions()?;
        self.extract_highlights(scale)?;
        self.extract_attachments()?;
        Ok(())
    }

    /// Atom coordinates with y inverted (so molecule and draw coordinates both
    /// go down the page as y increases) and the rotation of the options
    /// applied. The rotation is anti-clockwise as is conventional, but because
    /// y is inverted it is applied in the positive direction.
    pub fn extract_atom_coords(&mut self) -> Result<()> {
        if self.draw_mol.num_conformers() == 0 {
            return Err(RDKitError::precondition("no coords"));
        }
        let locs = self.draw_mol.conformer(self.conf_id)?.positions().to_vec();
        let rot = self.opts.rotate * PI / 180.0;
        let mut trans = Transform2D::identity();
        trans.set_transform_rotation(&Point2D::new(0.0, 0.0), rot);
        self.at_cds = locs
            .iter()
            .map(|p| {
                let mut pt = Point2D::new(p.x, -p.y);
                if rot != 0.0 {
                    trans.transform_point(&mut pt);
                }
                pt
            })
            .collect();
        Ok(())
    }

    /// Makes the atom labels, with their colours.
    pub fn extract_atom_symbols(&mut self, text: &mut DrawTextFT) -> Result<()> {
        self.atomic_nums.clear();
        for idx in 0..self.draw_mol.num_atoms() {
            let atom = self.draw_mol.atom(idx);
            let complex = crate::graphmol::queries::is_complex_atom_query(atom.query());
            self.atomic_nums
                .push(if complex { 0 } else { atom.atomic_num() as i32 });
            let (symbol, orient) =
                get_atom_symbol_and_orientation(&self.draw_mol, idx, &self.at_cds, &self.opts)?;
            self.atom_syms.push((symbol.clone(), orient));
            if symbol.is_empty() {
                self.atom_labels.push(None);
            } else {
                let col = self.get_colour(idx as i32)?;
                self.atom_labels.push(Some(AtomSymbol::new(
                    &symbol,
                    idx as i32,
                    orient,
                    self.at_cds[idx as usize],
                    col,
                    text,
                )?));
            }
        }
        Ok(())
    }

    /// The colour an atom's symbol is drawn in.
    pub fn get_colour(&self, atom_idx: i32) -> Result<DrawColour> {
        if atom_idx < 0 || atom_idx as usize >= self.atomic_nums.len() {
            return Err(RDKitError::precondition("bad atom_idx"));
        }
        let o = &self.opts;
        let mut retval = get_colour_by_atomic_num(self.atomic_nums[atom_idx as usize], o);
        let mut highlighted = self.highlight_atoms.contains(&atom_idx);
        if !o.circle_atoms && !o.continuous_highlight {
            if highlighted {
                retval = o.highlight_colour;
            }
            // an explicit colour from the highlight map overrides
            if let Some(c) = self.highlight_atom_map.get(&atom_idx) {
                highlighted = true;
                retval = *c;
            }
            // If it's not highlighted itself but all its bonds are, the atom
            // adopts the colour of the bonds (it doesn't look right if only
            // some are).
            if !highlighted {
                let (mut num_bonds, mut num_high_bonds) = (0, 0);
                let mut high_col: Option<DrawColour> = None;
                for nbr in self.draw_mol.atom_bonds(atom_idx as u32) {
                    num_bonds += 1;
                    let nbr_idx = nbr.idx() as i32;
                    if self.highlight_bonds.contains(&nbr_idx)
                        || self.highlight_bond_map.contains_key(&nbr_idx)
                    {
                        let hc = get_highlight_bond_colour(
                            &self.draw_mol,
                            nbr.idx(),
                            o,
                            &self.highlight_bonds,
                            &self.highlight_bond_map,
                            &self.highlight_atom_map,
                        );
                        match high_col {
                            None => high_col = Some(hc),
                            Some(c) if c != hc => {
                                num_high_bonds = 0;
                                break;
                            }
                            Some(_) => {}
                        }
                        num_high_bonds += 1;
                    }
                }
                if num_bonds == num_high_bonds {
                    if let Some(c) = high_col {
                        retval = c;
                    }
                }
            }
        } else if highlighted && !o.standard_colours_for_highlighted_atoms {
            // There's a colour behind the atom, so a symbol is drawn in the
            // colour of carbon.
            retval = o
                .atom_colour_palette
                .get(&6)
                .or_else(|| o.atom_colour_palette.get(&-1))
                .copied()
                .unwrap_or(DrawColour::new(0.0, 0.0, 0.0));
        }
        Ok(retval)
    }

    /// Mean bond length of the molecule, found once.
    pub fn calc_mean_bond_length(&mut self) -> Result<()> {
        if self.mean_bond_length == 0.0 {
            self.mean_bond_length = mean_bond_length(&self.draw_mol, -1)?;
        }
        Ok(())
    }

    /// Finds the extremes of everything included by
    /// `opts.drawing_extents_include`.
    pub fn find_extremes(&mut self) {
        let inc = self.opts.drawing_extents_include;
        let mut ext = (self.x_min, self.x_max, self.y_min, self.y_max);
        if inc & draw_element::PRESHAPES != 0 {
            self.pre_shapes
                .iter()
                .for_each(|s| s.find_extremes(&mut ext));
        }
        if inc & draw_element::BONDS != 0 {
            self.bonds.iter().for_each(|s| s.find_extremes(&mut ext));
        }
        if inc & draw_element::ATOMLABELS != 0 {
            self.atom_labels
                .iter()
                .flatten()
                .for_each(|l| l.find_extremes(&mut ext));
        }
        if inc & draw_element::HIGHLIGHTS != 0 {
            self.highlights
                .iter()
                .for_each(|s| s.find_extremes(&mut ext));
        }
        if self.include_annotations && inc & draw_element::ANNOTATIONS != 0 {
            self.annotations
                .iter()
                .for_each(|a| a.find_extremes(&mut ext, 0.0));
        }
        if inc & draw_element::RADICALS != 0 {
            // radicals are drawn in black, so they can extend the extents even
            // if the atom is not drawn
            find_radical_extremes(&self.radicals, &mut ext);
        }
        if inc & draw_element::POSTSHAPES != 0 {
            self.post_shapes
                .iter()
                .for_each(|s| s.find_extremes(&mut ext));
        }
        (self.x_min, self.x_max, self.y_min, self.y_max) = ext;
        if self.at_cds.is_empty() {
            self.x_min = -1.0;
            self.y_min = -1.0;
            self.x_max = 1.0;
            self.y_max = 1.0;
        }
        // The spans don't include the padding: that is in draw_width and
        // draw_height.
        self.x_range = self.x_max - self.x_min;
        self.y_range = self.y_max - self.y_min;
        if self.x_range < 1e-4 {
            self.x_range = 2.0;
            self.x_min -= 1.0;
            self.x_max += 1.0;
        }
        if self.y_range < 1e-4 {
            self.y_range = 2.0;
            self.y_min -= 1.0;
            self.y_max += 1.0;
        }
    }

    /// Works out the scale (and canvas size for flexible canvases).
    pub fn calculate_scale(&mut self, text: &mut DrawTextFT) {
        self.find_extremes();
        let o = &self.opts;
        // If width < 0 the scale comes off the y range, and likewise for
        // height and the x range. If both are negative, the scaling factor of
        // the options is used.
        let mut new_scale = 1.0;
        if self.width < 0 && self.height < 0 {
            self.width =
                (o.scaling_factor * self.x_range * (1.0 + 2.0 * self.margin_padding)) as i32;
            self.mol_height =
                (o.scaling_factor * self.y_range * (1.0 + 2.0 * self.margin_padding)) as i32;
        } else if self.width < 0 && self.y_range > 1.0e-4 {
            new_scale = f64::from(self.height) / self.y_range;
            // very wide and short molecules (HO-NH2) mustn't get long bonds
            let mbl = self.mean_bond_length * new_scale;
            if mbl > f64::from(self.mol_height / 2) {
                new_scale *= f64::from(self.mol_height / 2) / mbl;
            }
            self.width = (new_scale * self.x_range) as i32;
        } else if self.height < 0 && self.x_range > 1.0e-4 {
            new_scale = f64::from(self.width) / self.x_range;
            let mbl = self.mean_bond_length * new_scale;
            if mbl > f64::from(self.width / 2) {
                new_scale *= f64::from(self.width / 2) / mbl;
            }
            self.mol_height = (new_scale * self.y_range) as i32;
        }
        if self.height < 0 {
            self.height = self.mol_height;
            if self.legend.is_empty() {
                self.legend_height = 0;
            }
        }
        self.draw_width = (f64::from(self.width) * (1.0 - 2.0 * self.margin_padding)) as i32;
        self.draw_height = (f64::from(self.height) * (1.0 - 2.0 * self.margin_padding)) as i32;
        self.partition_for_legend();

        let o = &self.opts;
        if self.x_range > 1e-4 || self.y_range > 1e-4 {
            let mut width_for_scale = if self.mol_width > 0 {
                f64::from(self.mol_width)
            } else {
                f64::from(self.draw_width)
            };
            let side_legend_horiz = !self.legend.is_empty()
                && matches!(
                    o.legend_position,
                    LegendPosition::Left | LegendPosition::Right
                )
                && !o.legend_vertical_text;
            if side_legend_horiz && self.mol_width > 0 {
                // keep a minimum side gap of 8 px so labels don't crowd the
                // horizontal side legend on small canvases
                let g = 8.0f64.max(self.margin_padding * f64::from(self.draw_width));
                // and at least ~24 px of drawable molecule width after the gap
                if f64::from(self.mol_width) > g + 24.0 {
                    width_for_scale = f64::from(self.mol_width) - g;
                }
            }
            new_scale =
                (width_for_scale / self.x_range).min(f64::from(self.mol_height) / self.y_range);
            let mut fix_scale = new_scale;
            // Use the fixed scale unless it's too big, in which case scale the
            // drawing down to fit. fixedScale wins over fixedBondLength.
            if o.fixed_bond_length > 0.0 {
                fix_scale = o.fixed_bond_length;
            }
            if o.fixed_scale > 0.0 {
                fix_scale = width_for_scale * o.fixed_scale;
            }
            if new_scale > fix_scale {
                new_scale = fix_scale;
            }
        }
        let scale_mult = new_scale / self.scale;
        self.scale *= scale_mult;
        if o.fixed_font_size != -1 {
            self.font_scale = f64::from(o.fixed_font_size) / text.base_font_size();
        } else {
            self.font_scale *= scale_mult;
        }
    }

    /// Where the molecule coordinates go: the translation, scale and
    /// centring offset of the transformation to drawing coordinates.
    pub fn draw_transformers(&self) -> (Point2D, Point2D, Point2D) {
        let trans = Point2D::new(-self.x_min, -self.y_min);
        let scale = Point2D::new(self.scale, self.scale);
        let scaled_ranges = Point2D::new(self.scale * self.x_range, self.scale * self.y_range);
        let mol_w = if self.mol_width > 0 {
            f64::from(self.mol_width)
        } else {
            f64::from(self.draw_width)
        };
        let mut x_centre = (mol_w - scaled_ranges.x) / 2.0
            + self.x_offset
            + f64::from(self.width) * self.margin_padding;
        let mut y_centre = (f64::from(self.mol_height) - scaled_ranges.y) / 2.0
            + self.y_offset
            + f64::from(self.height) * self.margin_padding;
        let base_panel_x = self.x_offset + f64::from(self.width) * self.margin_padding;
        let o = &self.opts;
        // the side gap used by horizontal side legends, as in calculate_scale
        let side_gap = || 8.0f64.max(self.margin_padding * f64::from(self.draw_width));
        let horizontal_side_legend =
            !self.legend.is_empty() && !o.legend_vertical_text && self.mol_width > 0;
        match o.legend_position {
            LegendPosition::Bottom => {}
            LegendPosition::Top => y_centre += f64::from(self.legend_height),
            LegendPosition::Left => {
                if horizontal_side_legend {
                    x_centre = base_panel_x + f64::from(self.legend_width) + side_gap();
                } else {
                    x_centre += f64::from(self.legend_width);
                }
            }
            LegendPosition::Right => {
                if horizontal_side_legend {
                    x_centre =
                        base_panel_x + f64::from(self.mol_width) - side_gap() - scaled_ranges.x;
                }
            }
        }
        (trans, scale, Point2D::new(x_centre, y_centre))
    }

    /// Molecule coordinates to drawing coordinates with the given
    /// transformation (y is inverted).
    pub fn draw_coords_with(
        at_cds: Point2D,
        trans: Point2D,
        scale_factor: Point2D,
        to_centre: Point2D,
    ) -> Point2D {
        let mut p = Point2D::new(at_cds.x, -at_cds.y);
        p += trans;
        p.x *= scale_factor.x;
        p.y *= scale_factor.y;
        p += to_centre;
        p
    }

    /// Molecule coordinates to drawing coordinates.
    pub fn draw_coords(&self, at_cds: Point2D) -> Point2D {
        let (trans, scale, to_centre) = self.draw_transformers();
        Self::draw_coords_with(at_cds, trans, scale, to_centre)
    }

    /// Drawing coordinates of an atom.
    pub fn atom_draw_coords(&self, atnum: usize) -> Result<Point2D> {
        let p = self
            .at_cds
            .get(atnum)
            .ok_or_else(|| RDKitError::precondition("bad atom number"))?;
        Ok(self.draw_coords(Point2D::new(p.x, -p.y)))
    }

    /// Screen coordinates back to molecule coordinates.
    pub fn atom_coords_from_screen(&self, screen: Point2D) -> Point2D {
        let (trans, scale, to_centre) = self.draw_transformers();
        let mut p = screen;
        p -= to_centre;
        p.x /= scale.x;
        p.y /= scale.y;
        p -= trans;
        Point2D::new(p.x, -p.y)
    }

    /// Coordinates of an atom as found by [`DrawMol::extract_atom_coords`].
    pub fn atom_coords(&self, atnum: usize) -> Result<Point2D> {
        self.at_cds
            .get(atnum)
            .copied()
            .ok_or_else(|| RDKitError::precondition("bad atom number"))
    }

    /// Applies translation, scale and centring (in that order, each optional)
    /// to a point.
    pub fn transform_point(
        pt: Point2D,
        trans: Option<Point2D>,
        scale: Option<Point2D>,
        to_centre: Option<Point2D>,
    ) -> Point2D {
        let mut p = pt;
        if let Some(t) = trans {
            p += t;
        }
        if let Some(s) = scale {
            p.x *= s.x;
            p.y *= s.y;
        }
        if let Some(c) = to_centre {
            p += c;
        }
        p
    }

    /// Applies the transformations to everything.
    pub fn transform_all(
        &mut self,
        trans: Option<Point2D>,
        scale: Option<Point2D>,
        to_centre: Option<Point2D>,
        text: &mut DrawTextFT,
    ) -> Result<()> {
        let shape_all = |shapes: &mut Vec<DrawShape>| {
            for s in shapes {
                if let Some(t) = trans {
                    s.translate(t);
                }
                if let Some(f) = scale {
                    s.scale(f);
                }
                if let Some(c) = to_centre {
                    s.translate(c);
                }
            }
        };
        shape_all(&mut self.pre_shapes);
        shape_all(&mut self.bonds);
        shape_all(&mut self.highlights);
        for annot in &mut self.annotations {
            if let Some(t) = trans {
                annot.translate(t);
            }
            if let Some(f) = scale {
                annot.scale(f, text)?;
            }
            if let Some(c) = to_centre {
                annot.translate(c);
            }
        }
        for label in self.atom_labels.iter_mut().flatten() {
            if let Some(t) = trans {
                label.translate(t);
            }
            if let Some(f) = scale {
                label.scale(f, text)?;
            }
            if let Some(c) = to_centre {
                label.translate(c);
            }
        }
        // Radicals are StringRects without their own class: they are moved by
        // the scale and scaled by the font scale.
        for rad in &mut self.radicals {
            rad.0.trans = Self::transform_point(rad.0.trans, trans, scale, to_centre);
            rad.0.width *= self.font_scale;
            rad.0.height *= self.font_scale;
        }
        shape_all(&mut self.post_shapes);
        Ok(())
    }

    /// Moves everything to drawing coordinates.
    pub fn change_to_draw_coords(&mut self, text: &mut DrawTextFT) -> Result<()> {
        let (trans, scale, to_centre) = self.draw_transformers();
        self.transform_all(Some(trans), Some(scale), Some(to_centre), text)
    }

    /// Sets the offsets for drawing in a grid, say. The previous offsets are
    /// removed first, so errors accumulate a little if this is done a lot.
    pub fn set_offsets(
        &mut self,
        x_offset: f64,
        y_offset: f64,
        text: &mut DrawTextFT,
    ) -> Result<()> {
        if self.x_offset.abs() > 1e-4 || self.y_offset.abs() > 1e-4 {
            let back = Point2D::new(-self.x_offset, -self.y_offset);
            self.transform_all(Some(back), None, None, text)?;
        }
        self.x_offset = x_offset;
        self.y_offset = y_offset;
        self.transform_all(Some(Point2D::new(x_offset, y_offset)), None, None, text)
    }

    /// Clears everything extracted, back to the state after construction.
    pub fn reset_everything(&mut self, text: &mut DrawTextFT) {
        self.scale = 1.0;
        self.font_scale = 1.0;
        text.set_font_scale(1.0, true);
        let (x_min, y_min, x_max, y_max, x_range, y_range) = Self::initial_limits();
        self.x_min = x_min;
        self.y_min = y_min;
        self.x_max = x_max;
        self.y_max = y_max;
        self.x_range = x_range;
        self.y_range = y_range;
        self.mean_bond_length = 0.0;
        self.at_cds.clear();
        self.bonds.clear();
        self.pre_shapes.clear();
        self.post_shapes.clear();
        self.atomic_nums.clear();
        self.atom_syms.clear();
        self.atom_labels.clear();
        self.highlights.clear();
        self.annotations.clear();
        self.legends.clear();
        self.radicals.clear();
        self.single_bond_lines.clear();
    }

    /// Reduces the canvas to just accommodate the ranges at the current scale
    /// and recentres everything.
    pub fn shrink_to_fit(&mut self, with_padding: bool, text: &mut DrawTextFT) -> Result<()> {
        let padding = if with_padding {
            self.margin_padding
        } else {
            0.0
        };
        let new_width = ((2.0 * padding + 1.0) * self.x_range * self.scale).ceil() as i32;
        let new_height = ((2.0 * padding + 1.0) * self.y_range * self.scale).ceil() as i32;
        let corr = Point2D::new(
            f64::from(new_width - self.width) / 2.0,
            f64::from(new_height - self.height) / 2.0,
        );
        self.transform_all(Some(corr), None, None, text)?;
        self.width = new_width;
        self.draw_width = (f64::from(self.width) * (1.0 - 2.0 * padding)) as i32;
        self.height = new_height;
        if !self.legend.is_empty() {
            self.partition_for_legend();
            self.legends.clear();
            self.extract_legend(text)?;
        } else {
            self.legend_height = 0;
            self.mol_height = self.height;
            self.mol_width = self.draw_width;
            self.draw_height = (f64::from(self.height) * (1.0 - 2.0 * padding)) as i32;
        }
        Ok(())
    }

    /// Everything of the drawing is in drawing coordinates after this.
    pub fn finish_create_draw_objects(&mut self, text: &mut DrawTextFT) -> Result<()> {
        // the legend needs the final scale to get its font the right size
        self.extract_legend(text)?;
        self.change_to_draw_coords(text)?;
        self.drawing_initialised = true;
        Ok(())
    }
}

/// `addAtomIndices`: annotations with atom indices.
pub fn add_atom_indices(mol: &mut ROMol) -> Result<()> {
    // not in the global set of tags: only used here
    if mol.has_prop("_atomIndicesAdded") {
        return Ok(());
    }
    mol.set_prop("_atomIndicesAdded", 1, true)?;
    for idx in 0..mol.num_atoms() {
        let mut lab = idx.to_string();
        if let Some(v) = mol.atom(idx).get_prop(cp::atomNote) {
            lab = format!("{lab},{}", v.to_rd_string());
        }
        mol.atom_mut(idx).set_prop(cp::atomNote, lab, false)?;
    }
    Ok(())
}

/// `addBondIndices`: annotations with bond indices.
pub fn add_bond_indices(mol: &mut ROMol) -> Result<()> {
    if mol.has_prop("_bondIndicesAdded") {
        return Ok(());
    }
    mol.set_prop("_bondIndicesAdded", 1, true)?;
    for idx in 0..mol.num_bonds() {
        let mut lab = idx.to_string();
        if let Some(v) = mol.bond(idx).get_prop(cp::bondNote) {
            lab = format!("{lab},{}", v.to_rd_string());
        }
        mol.bond_mut(idx).set_prop(cp::bondNote, lab, false)?;
    }
    Ok(())
}
