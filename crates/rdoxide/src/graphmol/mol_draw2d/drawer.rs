//! The 2D molecule drawer (`MolDraw2D.cpp` of RDKit 2026.09.1) with the SVG
//! output of `MolDraw2DSVG`.
//!
//! Upstream splits the drawer into an abstract `MolDraw2D` and one subclass
//! per output format. Only SVG is ported, so the two are one type here:
//! [`MolDraw2DSVG`] holds the drawer state and the [`SvgCanvas`] that writes
//! the SVG. Reactions, multi-colour highlights and the comic style are not
//! ported yet.
//!
//! Copyright (C) 2015-2021 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;

use super::draw_mol::{DrawMol, DrawMolHighlights};
use super::draw_mol_utils::get_bond_highlights_for_atoms;
use super::draw_shape::ShapeCanvas;
use super::draw_text::{DrawTextFT, atom_label_to_pieces};
use super::helpers::{DrawColour, MolDrawOptions, OrientType, TextAlignType};
use super::svg_canvas::SvgCanvas;
use crate::geometry::Point2D;
use crate::graphmol::ROMol;
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::{RDKitError, Result};

/// Highlighting of one molecule for [`MolDraw2DSVG::draw_molecules`].
#[derive(Clone, Debug, Default)]
pub struct MoleculeHighlights {
    pub atoms: Option<Vec<i32>>,
    pub bonds: Option<Vec<i32>>,
    pub atom_map: Option<BTreeMap<i32, DrawColour>>,
    pub bond_map: Option<BTreeMap<i32, DrawColour>>,
    pub radii: Option<BTreeMap<i32, f64>>,
}

/// A drawer that makes SVG.
#[derive(Debug)]
pub struct MolDraw2DSVG {
    pub(super) canvas: SvgCanvas,
    pub(super) text: DrawTextFT,
    pub(super) width: i32,
    pub(super) height: i32,
    pub(super) panel_width: i32,
    pub(super) panel_height: i32,
    legend_height: i32,
    force_scale: bool,
    flexi_mode: bool,
    scale: f64,
    font_scale: f64,
    x_offset: i32,
    y_offset: i32,
    active_mol_idx: i32,
    pub(super) draw_mols: Vec<DrawMol>,
    global_draw_trans: Option<DrawMol>,
    pub(super) needs_init: bool,
}

impl MolDraw2DSVG {
    /// A drawer for a canvas of the given size, in pixels. Set a size to -1
    /// to have it follow from the molecule. The panel size is that of a
    /// single molecule when drawing several; a size of 0 or less means the
    /// whole canvas.
    pub fn new(width: i32, height: i32, panel_width: i32, panel_height: i32) -> Result<Self> {
        Self::with_options(
            width,
            height,
            panel_width,
            panel_height,
            MolDrawOptions::default(),
        )
    }

    pub fn with_options(
        width: i32,
        height: i32,
        panel_width: i32,
        panel_height: i32,
        options: MolDrawOptions,
    ) -> Result<Self> {
        let mut text = Self::make_text_drawer(&options);
        if options.base_font_size > 0.0 {
            text.set_base_font_size(options.base_font_size);
        }
        let mut drawer = MolDraw2DSVG {
            canvas: SvgCanvas::new(options),
            text,
            width,
            height,
            panel_width: if panel_width > 0 { panel_width } else { width },
            panel_height: if panel_height > 0 {
                panel_height
            } else {
                height
            },
            legend_height: 0,
            force_scale: false,
            flexi_mode: false,
            scale: 1.0,
            font_scale: 1.0,
            x_offset: 0,
            y_offset: 0,
            active_mol_idx: -1,
            draw_mols: Vec::new(),
            global_draw_trans: None,
            needs_init: true,
        };
        if width > 0 && height > 0 {
            drawer.init_drawing();
            drawer.needs_init = false;
        }
        Ok(drawer)
    }

    /// The text drawer for the font in the options; the built-in font is
    /// used if the file cannot be read.
    fn make_text_drawer(options: &MolDrawOptions) -> DrawTextFT {
        let max = f64::from(options.max_font_size);
        let min = f64::from(options.min_font_size);
        DrawTextFT::new(max, min, &options.font_file)
            .or_else(|_| DrawTextFT::new(max, min, ""))
            .expect("the built-in font is always available")
    }

    pub fn draw_options(&self) -> &MolDrawOptions {
        &self.canvas.options
    }

    pub fn draw_options_mut(&mut self) -> &mut MolDrawOptions {
        &mut self.canvas.options
    }

    pub fn width(&self) -> i32 {
        self.width
    }

    pub fn height(&self) -> i32 {
        self.height
    }

    pub fn panel_width(&self) -> i32 {
        self.panel_width
    }

    pub fn panel_height(&self) -> i32 {
        self.panel_height
    }

    pub fn draw_height(&self) -> i32 {
        self.panel_height - self.legend_height
    }

    pub fn flexi_mode(&self) -> bool {
        self.flexi_mode
    }

    pub fn set_flexi_mode(&mut self, mode: bool) {
        self.flexi_mode = mode;
        if mode {
            self.panel_width = -1;
            self.panel_height = -1;
        }
    }

    pub fn set_offset(&mut self, x: i32, y: i32) {
        self.x_offset = x;
        self.y_offset = y;
    }

    pub fn offset(&self) -> Point2D {
        Point2D::new(f64::from(self.x_offset), f64::from(self.y_offset))
    }

    // ------------------------------------------------------------------
    // state
    // ------------------------------------------------------------------

    pub fn colour(&self) -> DrawColour {
        self.canvas.colour
    }

    pub fn set_colour(&mut self, col: DrawColour) {
        self.canvas.colour = col;
    }

    pub fn dash(&self) -> &[f64] {
        &self.canvas.dash
    }

    pub fn set_dash(&mut self, patt: Vec<f64>) {
        self.canvas.dash = patt;
    }

    pub fn line_width(&self) -> f64 {
        self.canvas.options.bond_line_width
    }

    pub fn set_line_width(&mut self, width: f64) {
        self.canvas.options.bond_line_width = width;
    }

    pub fn fill_polys(&self) -> bool {
        self.canvas.fill_polys
    }

    pub fn set_fill_polys(&mut self, val: bool) {
        self.canvas.fill_polys = val;
    }

    pub fn font_size(&self) -> f64 {
        self.text.font_size()
    }

    pub fn set_font_size(&mut self, new_size: f64) {
        self.text.set_font_size(new_size);
    }

    pub fn active_class(&self) -> &str {
        &self.canvas.active_class
    }

    pub fn set_active_class(&mut self, class: &str) {
        self.canvas.active_class = class.to_string();
    }

    pub fn has_active_atm_idx(&self) -> bool {
        self.canvas.active_atm_idx1 >= 0
    }

    pub fn active_atm_idx1(&self) -> i32 {
        self.canvas.active_atm_idx1
    }

    pub fn active_atm_idx2(&self) -> i32 {
        self.canvas.active_atm_idx2
    }

    pub fn set_active_atm_idx(&mut self, atom1: i32, atom2: i32) {
        ShapeCanvas::set_active_atm_idx(&mut self.canvas, atom1, atom2);
    }

    pub fn has_active_bnd_idx(&self) -> bool {
        self.canvas.active_bnd_idx >= 0
    }

    pub fn active_bnd_idx(&self) -> i32 {
        self.canvas.active_bnd_idx
    }

    pub fn set_active_bnd_idx(&mut self, bond: i32) {
        self.canvas.active_bnd_idx = bond.max(-1);
    }

    pub fn set_active_mol_idx(&mut self, new_idx: i32) -> Result<()> {
        if new_idx < -1 || new_idx >= self.draw_mols.len() as i32 {
            return Err(RDKitError::precondition("bad new activeMolIdx_"));
        }
        self.active_mol_idx = new_idx;
        if let Some(dm) = self.active_draw_mol() {
            self.canvas.scale = dm.scale;
        }
        Ok(())
    }

    fn active_draw_mol(&self) -> Option<&DrawMol> {
        usize::try_from(self.active_mol_idx)
            .ok()
            .and_then(|i| self.draw_mols.get(i))
    }

    /// The scale of the molecule being drawn, as set by
    /// [`MolDraw2DSVG::set_scale`] if that has been called.
    pub fn scale(&self) -> f64 {
        self.canvas.scale
    }

    /// `getDrawLineWidth`
    pub fn draw_line_width(&self) -> f64 {
        self.canvas.draw_line_width()
    }

    // ------------------------------------------------------------------
    // coordinates
    // ------------------------------------------------------------------

    fn transforming_mol(&self) -> Result<&DrawMol> {
        self.global_draw_trans
            .as_ref()
            .or_else(|| self.active_draw_mol())
            .ok_or_else(|| RDKitError::precondition("no scaling info"))
    }

    /// Molecule coordinates to drawing coordinates, preferring the transform
    /// set with [`MolDraw2DSVG::set_scale_for_range`] if there is one.
    pub fn get_draw_coords(&self, mol_cds: Point2D) -> Result<Point2D> {
        Ok(self.transforming_mol()?.draw_coords(mol_cds))
    }

    /// Drawing coordinates of an atom of the active molecule.
    pub fn get_draw_coords_atom(&self, at_num: usize) -> Result<Point2D> {
        self.active_draw_mol()
            .ok_or_else(|| RDKitError::precondition("bad active mol index"))?
            .atom_draw_coords(at_num)
    }

    /// Screen coordinates back to molecule coordinates.
    pub fn get_atom_coords(&self, screen_cds: Point2D) -> Result<Point2D> {
        Ok(self.transforming_mol()?.atom_coords_from_screen(screen_cds))
    }

    /// Coordinates of an atom of the active molecule in the drawing's frame.
    pub fn get_atom_coords_atom(&self, at_num: usize) -> Result<Point2D> {
        self.active_draw_mol()
            .ok_or_else(|| RDKitError::precondition("bad active mol index"))?
            .atom_coords(at_num)
    }

    /// The atom coordinates of the active molecule.
    pub fn atom_coords(&self) -> Result<&[Point2D]> {
        Ok(&self
            .active_draw_mol()
            .ok_or_else(|| RDKitError::precondition("no index"))?
            .at_cds)
    }

    /// The atom symbols of the active molecule, with the orientation of each.
    pub fn atom_syms(&self) -> Result<&[(String, OrientType)]> {
        Ok(&self
            .active_draw_mol()
            .ok_or_else(|| RDKitError::precondition("no index"))?
            .atom_syms)
    }

    pub fn min_pt(&self) -> Result<Point2D> {
        let dm = self.transforming_mol()?;
        // the ys are inverted in the DrawMol
        Ok(Point2D::new(dm.x_min, -dm.y_max))
    }

    pub fn range(&self) -> Result<Point2D> {
        let dm = self.transforming_mol()?;
        Ok(Point2D::new(dm.x_range, dm.y_range))
    }

    /// `setScale(newScale)`: only affects the comic style upstream.
    pub fn set_scale_value(&mut self, new_scale: f64) {
        self.scale = new_scale;
    }

    /// `setScale(width, height, minv, maxv, mol)`: fixes the transformation
    /// for everything drawn from now on, from the extremes of the range and
    /// of a molecule if one is given.
    pub fn set_scale_for_range(
        &mut self,
        width: i32,
        height: i32,
        minv: Point2D,
        maxv: Point2D,
        mol: Option<&ROMol>,
    ) -> Result<()> {
        if width <= 0 {
            return Err(RDKitError::precondition("bad width"));
        }
        if height <= 0 {
            return Err(RDKitError::precondition("bad height"));
        }
        let (mut x_min, mut y_min, mut x_max, mut y_max);
        let mut set_font_scale = false;
        if let Some(mol) = mol {
            self.setup_text_drawer();
            let mut dm = DrawMol::new(
                mol,
                "",
                self.panel_width,
                self.panel_height,
                &self.canvas.options,
                &mut self.text,
                DrawMolHighlights::default(),
                true,
                -1,
                false,
            )?;
            dm.create_draw_objects(&mut self.text)?;
            // in the DrawMol, the ys are all inverted.
            x_min = minv.x.min(dm.x_min);
            y_min = minv.y.min(-dm.y_max);
            x_max = maxv.x.max(dm.x_max);
            y_max = maxv.y.max(-dm.y_min);
            self.font_scale = dm.font_scale;
            set_font_scale = true;
        } else {
            x_min = minv.x;
            y_min = minv.y;
            x_max = maxv.x;
            y_max = maxv.y;
        }
        let mut x_range = x_max - x_min;
        let mut y_range = y_max - y_min;
        if x_range < 1.0e-4 {
            x_range = 1.0;
            x_min = -0.5;
        }
        if y_range < 1.0e-4 {
            y_range = 1.0;
            y_min = -0.5;
        }

        // put a buffer round the drawing and calculate a final scale
        let padding = self.canvas.options.padding;
        x_min -= padding * x_range;
        x_range *= 1.0 + 2.0 * padding;
        x_max = x_min + x_range;
        y_min -= padding * y_range;
        y_range *= 1.0 + 2.0 * padding;
        y_max = y_min + y_range;

        let draw_width = f64::from(width) * (1.0 - 2.0 * padding);
        let draw_height = f64::from(height) * (1.0 - 2.0 * padding);

        self.scale = (draw_width / x_range).min(draw_height / y_range);
        // Absent any other information, we'll have to go with the font scale
        // the same as the scale.
        if !set_font_scale {
            // get the text drawer to decide on a suitable font scale.
            let ofs = self.text.font_scale();
            self.text.set_font_scale(self.scale, false);
            self.font_scale = self.text.font_scale();
            self.text.set_font_scale(ofs, true);
        } else {
            self.font_scale *= self.scale / self.font_scale;
        }
        self.force_scale = true;

        let mut global = DrawMol::without_mol(
            self.panel_width,
            self.panel_height,
            &self.canvas.options,
            &mut self.text,
            x_min,
            x_max,
            y_min,
            y_max,
            self.scale,
            self.font_scale,
        );
        // the padding is now a no-go area around the image, so it is applied
        // as an offset before each drawing move. It needs to be taken off
        // here to compensate.
        global.x_offset -= f64::from(width) * padding;
        global.y_offset -= f64::from(height) * padding;
        self.global_draw_trans = Some(global);
        Ok(())
    }

    // ------------------------------------------------------------------
    // text sizes
    // ------------------------------------------------------------------

    /// The size of a string in molecule units.
    pub fn get_string_size(&self, label: &str) -> Result<(f64, f64)> {
        let (w, h) = self.text.string_size(label)?;
        let scale = self.active_scale()?;
        Ok((w / scale, h / scale))
    }

    pub(super) fn active_scale_for_tags(&self) -> Result<f64> {
        self.active_scale()
    }

    fn active_scale(&self) -> Result<f64> {
        self.active_draw_mol()
            .map(|dm| dm.scale)
            .ok_or_else(|| RDKitError::precondition("bad active mol index"))
    }

    /// The size of an atom label, which may be in several pieces.
    pub fn get_label_size(&self, label: &str, orient: OrientType) -> Result<(f64, f64)> {
        if matches!(orient, OrientType::N | OrientType::S) {
            let mut width = 0.0f64;
            let mut height = 0.0;
            for bit in atom_label_to_pieces(label, orient) {
                let (w, h) = self.get_string_size(&bit)?;
                width = width.max(w);
                height += h;
            }
            Ok((width, height))
        } else {
            self.get_string_size(label)
        }
    }

    /// The extremes of a label in molecule coordinates, as
    /// `(x_min, y_min, x_max, y_max)`.
    pub fn get_string_extremes(
        &self,
        label: &str,
        orient: OrientType,
        cds: Point2D,
    ) -> Result<(f64, f64, f64, f64)> {
        let (mut x_min, mut y_min, mut x_max, mut y_max) =
            self.text.string_extremes(label, orient, false)?;
        let draw_cds = self.get_draw_coords(cds)?;
        x_min += draw_cds.x;
        x_max += draw_cds.x;
        y_min += draw_cds.y;
        y_max += draw_cds.y;
        let new_mins = self.get_atom_coords(Point2D::new(x_min, y_min))?;
        let new_maxs = self.get_atom_coords(Point2D::new(x_max, y_max))?;
        // draw coords to atom coords reverses y
        let (y_min, y_max) = if new_mins.y > new_maxs.y {
            (new_maxs.y, new_mins.y)
        } else {
            (new_mins.y, new_maxs.y)
        };
        Ok((new_mins.x, y_min, new_maxs.x, y_max))
    }

    // ------------------------------------------------------------------
    // drawing molecules
    // ------------------------------------------------------------------

    fn setup_text_drawer(&mut self) {
        let o = &self.canvas.options;
        self.text.set_max_font_size(f64::from(o.max_font_size));
        self.text.set_min_font_size(f64::from(o.min_font_size));
        if o.base_font_size > 0.0 {
            self.text.set_base_font_size(o.base_font_size);
        }
        if self.text.set_font_file(&o.font_file).is_err() {
            // fall back on the built-in font
            let _ = self.text.set_font_file("");
        }
    }

    fn check_supported(&self) -> Result<()> {
        if self.canvas.options.comic_mode {
            return Err(RDKitError::UnhandledFeature(
                "comic mode is not supported yet".into(),
            ));
        }
        Ok(())
    }

    /// Draws a molecule with a legend and highlights. The molecule is copied,
    /// as it needs ring information.
    pub fn draw_molecule_with(
        &mut self,
        mol: &ROMol,
        legend: &str,
        highlights: DrawMolHighlights,
        conf_id: i32,
    ) -> Result<()> {
        self.check_supported()?;
        self.setup_text_drawer();
        // we need ring info for drawing, so copy the molecule in order to
        // add it
        let mut lmol = mol.clone();
        if !mol.ring_info().is_symm_sssr() {
            symmetrize_sssr(
                &mut lmol,
                SymmetrizeSSSRAlgorithm::Default,
                true,
                false,
                false,
            )?;
        }
        let mut dm = DrawMol::new(
            &lmol,
            legend,
            self.panel_width,
            self.panel_height,
            &self.canvas.options,
            &mut self.text,
            highlights,
            true,
            conf_id,
            false,
        )?;
        dm.set_offsets(
            f64::from(self.x_offset),
            f64::from(self.y_offset),
            &mut self.text,
        )?;
        dm.create_draw_objects(&mut self.text)?;
        self.fix_variable_dimensions(&dm);
        self.draw_mols.push(dm);
        self.active_mol_idx += 1;
        self.start_drawing();
        self.draw_the_molecule(self.draw_mols.len() - 1)
    }

    /// Draws a molecule.
    pub fn draw_molecule(&mut self, mol: &ROMol, legend: &str, conf_id: i32) -> Result<()> {
        self.draw_molecule_with(mol, legend, DrawMolHighlights::default(), conf_id)
    }

    /// Draws a molecule with the given atoms highlighted, and the bonds between
    /// them.
    pub fn draw_molecule_atoms(
        &mut self,
        mol: &ROMol,
        legend: &str,
        highlight_atoms: &[i32],
        highlight_atom_map: Option<&BTreeMap<i32, DrawColour>>,
        highlight_radii: Option<&BTreeMap<i32, f64>>,
        conf_id: i32,
    ) -> Result<()> {
        let hl = DrawMolHighlights {
            atoms: highlight_atoms.to_vec(),
            bonds: get_bond_highlights_for_atoms(mol, highlight_atoms),
            atom_map: highlight_atom_map.cloned().unwrap_or_default(),
            radii: highlight_radii.cloned().unwrap_or_default(),
            ..Default::default()
        };
        self.draw_molecule_with(mol, legend, hl, conf_id)
    }

    /// Draws several molecules in a grid of panels. The canvas must have a
    /// fixed size; the panel size says how many fit in a row and column.
    pub fn draw_molecules(
        &mut self,
        mols: &[Option<&ROMol>],
        legends: Option<&[String]>,
        highlights: Option<&[MoleculeHighlights]>,
        conf_ids: Option<&[i32]>,
    ) -> Result<()> {
        let n = mols.len();
        for (what, len) in [
            ("legends", legends.map(<[_]>::len)),
            ("highlights", highlights.map(<[_]>::len)),
            ("confIds", conf_ids.map(<[_]>::len)),
        ] {
            if len.is_some_and(|l| l != n) {
                return Err(RDKitError::precondition(format!("bad size of {what}")));
            }
        }
        if self.panel_width == 0 {
            return Err(RDKitError::precondition("panel width cannot be zero"));
        }
        if self.panel_height == 0 {
            return Err(RDKitError::precondition("panel height cannot be zero"));
        }
        if self.width <= 0 || self.height <= 0 {
            return Err(RDKitError::precondition(
                "drawMolecules() needs a fixed canvas size",
            ));
        }
        if mols.is_empty() {
            return Ok(());
        }
        self.check_supported()?;
        self.setup_text_drawer();
        let mut min_scale_mol = 0usize;
        let mut min_font_scale_mol = 0usize;
        let n_cols = self.width / self.panel_width;
        let n_rows = self.height / self.panel_height;

        for (i, mol) in mols.iter().enumerate() {
            let Some(mol) = mol else { continue };
            let legend = legends.map_or("", |l| l[i].as_str());
            let mut hl = DrawMolHighlights::default();
            if let Some(h) = highlights.map(|h| &h[i]) {
                hl.atoms = h.atoms.clone().unwrap_or_default();
                hl.atom_map = h.atom_map.clone().unwrap_or_default();
                hl.bond_map = h.bond_map.clone().unwrap_or_default();
                hl.radii = h.radii.clone().unwrap_or_default();
                if let Some(bonds) = &h.bonds {
                    hl.bonds = bonds.clone();
                } else if self.canvas.options.continuous_highlight {
                    if let Some(atoms) = &h.atoms {
                        hl.bonds = get_bond_highlights_for_atoms(mol, atoms);
                    }
                }
            }
            let conf_id = conf_ids.map_or(-1, |c| c[i]);
            let mut dm = DrawMol::new(
                mol,
                legend,
                self.panel_width,
                self.panel_height,
                &self.canvas.options,
                &mut self.text,
                hl,
                true,
                conf_id,
                false,
            )?;
            let row = if n_rows > 1 { i as i32 / n_cols } else { 0 };
            let col = if n_cols > 1 { i as i32 % n_cols } else { 0 };
            dm.set_offsets(
                f64::from(col * self.panel_width),
                f64::from(row * self.panel_height),
                &mut self.text,
            )?;
            dm.create_draw_objects(&mut self.text)?;
            let idx = self.draw_mols.len();
            if idx > 0 {
                if dm.scale < self.draw_mols[min_scale_mol].scale {
                    min_scale_mol = idx;
                }
                if dm.font_scale < self.draw_mols[min_font_scale_mol].font_scale {
                    min_font_scale_mol = idx;
                }
            }
            self.draw_mols.push(dm);
        }

        if self.canvas.options.draw_mols_same_scale {
            let scale = self.draw_mols[min_scale_mol].scale;
            let font_scale = self.draw_mols[min_font_scale_mol].font_scale;
            for dm in &mut self.draw_mols {
                dm.set_scale(scale, font_scale, true, &mut self.text)?;
            }
        }

        if self.draw_mols.is_empty() {
            self.active_mol_idx = -1;
            return Ok(());
        }
        self.active_mol_idx = 0;
        self.start_drawing();
        for i in 0..self.draw_mols.len() {
            self.active_mol_idx = i as i32;
            self.draw_the_molecule(i)?;
        }
        Ok(())
    }

    /// The size a molecule would be drawn at, with the panel size as the
    /// target.
    pub fn get_mol_size(
        &mut self,
        mol: &ROMol,
        legend: &str,
        highlights: DrawMolHighlights,
        conf_id: i32,
    ) -> Result<(i32, i32)> {
        self.setup_text_drawer();
        let mut dm = DrawMol::new(
            mol,
            legend,
            self.panel_width,
            self.panel_height,
            &self.canvas.options,
            &mut self.text,
            highlights,
            true,
            conf_id,
            false,
        )?;
        dm.create_draw_objects(&mut self.text)?;
        Ok((dm.width, dm.height))
    }

    /// Takes the canvas size from a molecule where the drawer was asked to.
    fn fix_variable_dimensions(&mut self, dm: &DrawMol) {
        if self.panel_width == -1 {
            self.width = dm.width;
            if !self.flexi_mode {
                self.panel_width = self.width;
            }
        }
        if self.panel_height == -1 {
            self.height = dm.height;
            if !self.flexi_mode {
                self.panel_height = self.height;
            }
        }
    }

    fn start_drawing(&mut self) {
        if self.needs_init {
            self.init_drawing();
            self.needs_init = false;
        }
        if self.active_mol_idx <= 0 && self.canvas.options.clear_background {
            self.clear_drawing();
        }
    }

    fn draw_the_molecule(&mut self, idx: usize) -> Result<()> {
        let dm = &mut self.draw_mols[idx];
        if let Some(global) = &self.global_draw_trans {
            dm.set_transformation(global, &mut self.text)?;
        } else if self.force_scale {
            dm.set_scale(self.scale, self.font_scale, false, &mut self.text)?;
        }
        self.canvas.scale = dm.scale;
        dm.draw(&mut self.canvas, &mut self.text)
    }

    // ------------------------------------------------------------------
    // primitives; coordinates are molecule coordinates unless `raw_coords`
    // ------------------------------------------------------------------

    fn to_draw(&self, cds: Point2D, raw_coords: bool) -> Result<Point2D> {
        if raw_coords {
            Ok(cds)
        } else {
            self.get_draw_coords(cds)
        }
    }

    pub fn draw_line(&mut self, cds1: Point2D, cds2: Point2D, raw_coords: bool) -> Result<()> {
        let c = self.canvas.colour;
        self.draw_line_coloured(cds1, cds2, c, c, raw_coords)
    }

    pub fn draw_line_coloured(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col1: DrawColour,
        col2: DrawColour,
        raw_coords: bool,
    ) -> Result<()> {
        let (c1, c2) = (
            self.to_draw(cds1, raw_coords)?,
            self.to_draw(cds2, raw_coords)?,
        );
        self.canvas.draw_line(c1, c2, col1, col2, true)
    }

    pub fn draw_polygon(&mut self, cds: &[Point2D], raw_coords: bool) -> Result<()> {
        let pts = cds
            .iter()
            .map(|&c| self.to_draw(c, raw_coords))
            .collect::<Result<Vec<_>>>()?;
        self.canvas.draw_polygon(&pts, true)
    }

    pub fn draw_triangle(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        cds3: Point2D,
        raw_coords: bool,
    ) -> Result<()> {
        self.draw_polygon(&[cds1, cds2, cds3], raw_coords)
    }

    pub fn draw_ellipse(&mut self, cds1: Point2D, cds2: Point2D, raw_coords: bool) -> Result<()> {
        let (c1, c2) = (
            self.to_draw(cds1, raw_coords)?,
            self.to_draw(cds2, raw_coords)?,
        );
        self.canvas.draw_ellipse(c1, c2, true)
    }

    pub fn draw_arc(
        &mut self,
        centre: Point2D,
        xradius: f64,
        yradius: f64,
        ang1: f64,
        ang2: f64,
        raw_coords: bool,
    ) -> Result<()> {
        let c = self.to_draw(centre, raw_coords)?;
        self.canvas.draw_arc(c, xradius, yradius, ang1, ang2, true)
    }

    pub fn draw_rect(&mut self, cds1: Point2D, cds2: Point2D, raw_coords: bool) -> Result<()> {
        let mut pts = vec![
            cds1,
            Point2D::new(cds1.x, cds2.y),
            cds2,
            Point2D::new(cds2.x, cds1.y),
        ];
        // if fillPolys() is false, it doesn't close the polygon because of
        // its use for drawing filled or open ellipse segments.
        if !self.canvas.fill_polys {
            pts.push(cds1);
        }
        self.draw_polygon(&pts, raw_coords)
    }

    pub fn draw_wavy_line(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col1: DrawColour,
        col2: DrawColour,
        n_segments: u32,
        vert_offset: f64,
        raw_coords: bool,
    ) -> Result<()> {
        let (c1, c2) = (
            self.to_draw(cds1, raw_coords)?,
            self.to_draw(cds2, raw_coords)?,
        );
        self.canvas
            .draw_wavy_line(c1, c2, col1, col2, n_segments, vert_offset, true)
    }

    /// A wavy line across the end of `cds1`-`cds2`, perpendicular to it.
    pub fn draw_attachment_line(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        col: DrawColour,
        len: f64,
        n_segments: u32,
        raw_coords: bool,
    ) -> Result<()> {
        let perp = super::draw_mol_utils::calc_perpendicular(cds1, cds2);
        let p1 = Point2D::new(cds2.x - perp.x * len / 2.0, cds2.y - perp.y * len / 2.0);
        let p2 = Point2D::new(cds2.x + perp.x * len / 2.0, cds2.y + perp.y * len / 2.0);
        self.draw_wavy_line(p1, p2, col, col, n_segments, 0.05, raw_coords)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn draw_arrow(
        &mut self,
        cds1: Point2D,
        cds2: Point2D,
        as_polygon: bool,
        frac: f64,
        angle: f64,
        col: DrawColour,
        raw_coords: bool,
    ) -> Result<()> {
        let (c1, c2) = (
            self.to_draw(cds1, raw_coords)?,
            self.to_draw(cds2, raw_coords)?,
        );
        self.canvas
            .draw_arrow(c1, c2, as_polygon, frac, angle, col, true)
    }

    pub fn draw_plus(
        &mut self,
        cds: Point2D,
        plus_width: i32,
        col: DrawColour,
        raw_coords: bool,
    ) -> Result<()> {
        let w = f64::from(plus_width);
        self.draw_line_coloured(
            Point2D::new(cds.x, cds.y - w),
            Point2D::new(cds.x, cds.y + w),
            col,
            col,
            raw_coords,
        )?;
        self.draw_line_coloured(
            Point2D::new(cds.x - w, cds.y),
            Point2D::new(cds.x + w, cds.y),
            col,
            col,
            raw_coords,
        )
    }

    /// Draws a string centred on the point.
    pub fn draw_string(&mut self, text: &str, cds: Point2D, raw_coords: bool) -> Result<()> {
        let draw_cds = self.to_draw(cds, raw_coords)?;
        self.text.set_colour(self.canvas.colour);
        let svg =
            self.text
                .draw_label_svg(text, draw_cds, OrientType::N, &self.canvas.active_class)?;
        self.canvas.out.push_str(&svg);
        Ok(())
    }

    /// Draws a string aligned to the point.
    pub fn draw_string_aligned(
        &mut self,
        text: &str,
        cds: Point2D,
        align: TextAlignType,
        raw_coords: bool,
    ) -> Result<()> {
        let draw_cds = self.to_draw(cds, raw_coords)?;
        self.text.set_colour(self.canvas.colour);
        let svg = self
            .text
            .draw_string_svg(text, draw_cds, align, &self.canvas.active_class)?;
        self.canvas.out.push_str(&svg);
        Ok(())
    }
}
