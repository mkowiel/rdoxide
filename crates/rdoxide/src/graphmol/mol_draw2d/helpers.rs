//! Colours, palettes, draw elements and drawing options of the 2D molecule
//! drawing code (`MolDraw2DHelpers.h` of RDKit 2026.09.1).
//!
//! Copyright (C) 2015-2026 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::ops::{Add, Div, Mul, Sub};

use crate::{RDKitError, Result};

/// How text is oriented relative to the atom it labels.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrientType {
    C,
    N,
    E,
    S,
    W,
}

/// Alignment of a string relative to the point it is drawn at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TextAlignType {
    Middle,
    Start,
    End,
}

/// RGBA colour with components in `[0, 1]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrawColour {
    pub r: f64,
    pub g: f64,
    pub b: f64,
    pub a: f64,
}

impl Default for DrawColour {
    fn default() -> Self {
        DrawColour::new(0.0, 0.0, 0.0)
    }
}

impl DrawColour {
    pub const fn new(r: f64, g: f64, b: f64) -> Self {
        DrawColour { r, g, b, a: 1.0 }
    }

    pub const fn with_alpha(r: f64, g: f64, b: f64, a: f64) -> Self {
        DrawColour { r, g, b, a }
    }

    /// Approximate equality with tolerance `tol`, optionally ignoring alpha.
    pub fn feq(&self, other: &DrawColour, tol: f64, ignore_alpha: bool) -> bool {
        (self.r - other.r).abs() <= tol
            && (self.g - other.g).abs() <= tol
            && (self.b - other.b).abs() <= tol
            && (ignore_alpha || (self.a - other.a).abs() <= tol)
    }

    /// Division by a scalar; upstream requires it to be non-zero.
    pub fn div_checked(&self, v: f64) -> Result<DrawColour> {
        if v == 0.0 {
            return Err(RDKitError::precondition("divide by zero"));
        }
        Ok(*self / v)
    }
}

impl Add for DrawColour {
    type Output = DrawColour;
    fn add(self, o: DrawColour) -> DrawColour {
        DrawColour::with_alpha(self.r + o.r, self.g + o.g, self.b + o.b, self.a + o.a)
    }
}

impl Sub for DrawColour {
    type Output = DrawColour;
    fn sub(self, o: DrawColour) -> DrawColour {
        DrawColour::with_alpha(self.r - o.r, self.g - o.g, self.b - o.b, self.a - o.a)
    }
}

impl Mul<f64> for DrawColour {
    type Output = DrawColour;
    fn mul(self, v: f64) -> DrawColour {
        DrawColour::with_alpha(self.r * v, self.g * v, self.b * v, self.a * v)
    }
}

impl Div<f64> for DrawColour {
    type Output = DrawColour;
    fn div(self, v: f64) -> DrawColour {
        DrawColour::with_alpha(self.r / v, self.g / v, self.b / v, self.a / v)
    }
}

/// Atomic number (or the special keys -1, 0, 201) to colour.
pub type ColourPalette = BTreeMap<i32, DrawColour>;
pub type DashPattern = Vec<f64>;

/// Converts the line width into something SVG and Cairo use.
pub const LINE_WIDTH_SCALE_FACTOR: f64 = 0.02;

fn palette_from(entries: &[(i32, DrawColour)]) -> ColourPalette {
    entries.iter().copied().collect()
}

const BLACK: DrawColour = DrawColour::new(0.0, 0.0, 0.0);

/// RDKit's default palette. Key 201 is for hydrogens drawn without symbols.
pub fn default_palette() -> ColourPalette {
    palette_from(&[
        (-1, BLACK),
        (0, DrawColour::new(0.1, 0.1, 0.1)),
        (1, BLACK),
        (6, BLACK),
        (7, DrawColour::new(0.0, 0.0, 1.0)),
        (8, DrawColour::new(1.0, 0.0, 0.0)),
        (9, DrawColour::new(0.2, 0.8, 0.8)),
        (15, DrawColour::new(1.0, 0.5, 0.0)),
        (16, DrawColour::new(0.8, 0.8, 0.0)),
        (17, DrawColour::new(0.0, 0.802, 0.0)),
        (35, DrawColour::new(0.5, 0.3, 0.1)),
        (53, DrawColour::new(0.63, 0.12, 0.94)),
        (201, DrawColour::new(0.68, 0.85, 0.90)),
    ])
}

/// The palette of the Avalon renderer.
pub fn avalon_palette() -> ColourPalette {
    palette_from(&[
        (-1, BLACK),
        (0, DrawColour::new(0.1, 0.1, 0.1)),
        (1, BLACK),
        (6, BLACK),
        (7, DrawColour::new(0.0, 0.0, 1.0)),
        (8, DrawColour::new(1.0, 0.0, 0.0)),
        (9, DrawColour::new(0.0, 0.498, 0.0)),
        (15, DrawColour::new(0.498, 0.0, 0.498)),
        (16, DrawColour::new(0.498, 0.247, 0.0)),
        (17, DrawColour::new(0.0, 0.498, 0.0)),
        (35, DrawColour::new(0.0, 0.498, 0.0)),
        (53, DrawColour::new(0.247, 0.0, 0.498)),
        (201, DrawColour::new(0.68, 0.85, 0.90)),
    ])
}

/// Part of the CDK colour palette (CDK2DAtomColors).
pub fn cdk_palette() -> ColourPalette {
    palette_from(&[
        (-1, BLACK),
        (0, DrawColour::new(0.1, 0.1, 0.1)),
        (1, BLACK),
        (6, BLACK),
        (7, DrawColour::new(0.188, 0.314, 0.972)),
        (8, DrawColour::new(1.0, 0.051, 0.051)),
        (9, DrawColour::new(0.565, 0.878, 0.314)),
        (15, DrawColour::new(1.0, 0.5, 0.0)),
        (16, DrawColour::new(0.776, 0.776, 0.173)),
        (17, DrawColour::new(0.122, 0.498, 0.122)),
        (35, DrawColour::new(0.651, 0.161, 0.161)),
        (53, DrawColour::new(0.580, 0.0, 0.580)),
        (5, DrawColour::new(1.000, 0.710, 0.710)),
        (201, DrawColour::new(0.68, 0.85, 0.90)),
    ])
}

pub fn dark_mode_palette() -> ColourPalette {
    palette_from(&[
        (-1, DrawColour::new(0.8, 0.8, 0.8)),
        (0, DrawColour::new(0.9, 0.9, 0.9)),
        (1, DrawColour::new(0.9, 0.9, 0.9)),
        (6, DrawColour::new(0.9, 0.9, 0.9)),
        (7, DrawColour::new(0.33, 0.41, 0.92)),
        (8, DrawColour::new(1.0, 0.2, 0.2)),
        (9, DrawColour::new(0.2, 0.8, 0.8)),
        (15, DrawColour::new(1.0, 0.5, 0.0)),
        (16, DrawColour::new(0.8, 0.8, 0.0)),
        (17, DrawColour::new(0.0, 0.802, 0.0)),
        (35, DrawColour::new(0.71, 0.4, 0.07)),
        (53, DrawColour::new(0.89, 0.004, 1.0)),
        (201, DrawColour::new(0.68, 0.85, 0.90)),
    ])
}

pub fn bw_palette() -> ColourPalette {
    palette_from(&[(-1, BLACK)])
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MultiColourHighlightStyle {
    CircleAndLine,
    Lasso,
}

/// Bit flags selecting the parts of a drawing (`DrawElement`).
pub mod draw_element {
    pub const NONE: u32 = 0;
    pub const PRESHAPES: u32 = 1 << 0;
    pub const BONDS: u32 = 1 << 1;
    pub const ATOMLABELS: u32 = 1 << 2;
    pub const HIGHLIGHTS: u32 = 1 << 3;
    pub const ANNOTATIONS: u32 = 1 << 4;
    pub const RADICALS: u32 = 1 << 5;
    pub const POSTSHAPES: u32 = 1 << 6;
    pub const ALL: u32 = 0x7fff_ffff;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LegendPosition {
    Bottom,
    Top,
    Left,
    Right,
}

/// Options controlling how molecules are drawn (`MolDrawOptions`).
#[derive(Clone, Debug)]
pub struct MolDrawOptions {
    pub atom_label_deuterium_tritium: bool,
    pub dummies_are_attachments: bool,
    pub circle_atoms: bool,
    pub split_bonds: bool,
    pub highlight_colour: DrawColour,
    pub continuous_highlight: bool,
    pub fill_highlights: bool,
    pub highlight_radius: f64,
    pub flag_close_contacts_dist: i32,
    pub include_atom_tags: bool,
    pub clear_background: bool,
    pub background_colour: DrawColour,
    pub query_colour: DrawColour,
    pub legend_font_size: i32,
    pub legend_fraction: f64,
    pub max_font_size: i32,
    pub min_font_size: i32,
    pub fixed_font_size: i32,
    pub annotation_font_scale: f64,
    pub font_file: String,
    pub legend_colour: DrawColour,
    pub legend_position: LegendPosition,
    pub legend_vertical_text: bool,
    pub multiple_bond_offset: f64,
    pub padding: f64,
    pub component_padding: f64,
    pub additional_atom_label_padding: f64,
    pub atom_labels: BTreeMap<i32, String>,
    pub no_atom_labels: bool,
    pub atom_regions: Vec<Vec<i32>>,
    pub symbol_colour: DrawColour,
    pub annotation_colour: DrawColour,
    pub atom_note_colour: DrawColour,
    pub bond_note_colour: DrawColour,
    pub bond_line_width: f64,
    pub scale_bond_width: bool,
    pub scale_highlight_bond_width: bool,
    pub highlight_bond_width_multiplier: i32,
    pub prepare_mols_before_drawing: bool,
    pub highlight_colour_palette: Vec<DrawColour>,
    pub atom_colour_palette: ColourPalette,
    pub fixed_scale: f64,
    pub fixed_bond_length: f64,
    pub rotate: f64,
    pub add_atom_indices: bool,
    pub add_bond_indices: bool,
    pub isotope_labels: bool,
    pub dummy_isotope_labels: bool,
    pub add_stereo_annotation: bool,
    pub show_all_cip_codes: bool,
    pub atom_highlights_are_circles: bool,
    pub multi_colour_highlight_style: MultiColourHighlightStyle,
    pub centre_molecules_before_drawing: bool,
    pub drawing_extents_include: u32,
    pub explicit_methyl: bool,
    pub include_radicals: bool,
    pub include_metadata: bool,
    pub comic_mode: bool,
    pub variable_bond_width_multiplier: i32,
    pub variable_atom_radius: f64,
    pub variable_attachment_colour: DrawColour,
    pub include_chiral_flag_label: bool,
    pub simplified_stereo_group_label: bool,
    pub stereo_group_and_label: String,
    pub stereo_group_or_label: String,
    pub stereo_group_abs_label: String,
    pub add_stereo_group_annotation: bool,
    pub unspecified_stereo_is_unknown: bool,
    pub single_colour_wedge_bonds: bool,
    pub single_colour_bonds: bool,
    pub use_mol_block_wedging: bool,
    pub scaling_factor: f64,
    pub base_font_size: f64,
    pub draw_mols_same_scale: bool,
    pub use_complex_query_atom_symbols: bool,
    pub brackets_around_atom_lists: bool,
    pub standard_colours_for_highlighted_atoms: bool,
}

impl Default for MolDrawOptions {
    fn default() -> Self {
        let c = DrawColour::new;
        MolDrawOptions {
            atom_label_deuterium_tritium: false,
            dummies_are_attachments: false,
            circle_atoms: true,
            split_bonds: false,
            highlight_colour: c(1.0, 0.5, 0.5),
            continuous_highlight: true,
            fill_highlights: true,
            highlight_radius: 0.3,
            flag_close_contacts_dist: 3,
            include_atom_tags: false,
            clear_background: true,
            background_colour: c(1.0, 1.0, 1.0),
            query_colour: c(0.5, 0.5, 0.5),
            legend_font_size: 16,
            legend_fraction: 0.1,
            max_font_size: 40,
            min_font_size: 6,
            fixed_font_size: -1,
            annotation_font_scale: 0.5,
            font_file: String::new(),
            legend_colour: c(0.0, 0.0, 0.0),
            legend_position: LegendPosition::Bottom,
            legend_vertical_text: true,
            multiple_bond_offset: 0.15,
            padding: 0.05,
            component_padding: 0.0,
            additional_atom_label_padding: 0.0,
            atom_labels: BTreeMap::new(),
            no_atom_labels: false,
            atom_regions: Vec::new(),
            symbol_colour: c(0.0, 0.0, 0.0),
            annotation_colour: c(0.0, 0.0, 0.0),
            atom_note_colour: c(0.0, 0.0, 0.0),
            bond_note_colour: c(0.5, 0.5, 1.0),
            bond_line_width: 2.0,
            scale_bond_width: false,
            scale_highlight_bond_width: true,
            highlight_bond_width_multiplier: 8,
            prepare_mols_before_drawing: true,
            highlight_colour_palette: vec![
                c(1.0, 1.0, 0.67),   // popcorn yellow
                c(1.0, 0.8, 0.6),    // sand
                c(1.0, 0.71, 0.76),  // light pink
                c(0.8, 1.0, 0.8),    // offwhitegreen
                c(0.87, 0.63, 0.87), // plum
                c(0.76, 0.94, 0.96), // pastel blue
                c(0.67, 0.67, 1.0),  // periwinkle
                c(0.64, 0.76, 0.34), // avocado
                c(0.56, 0.93, 0.56), // light green
                c(0.20, 0.63, 0.79), // peacock
            ],
            atom_colour_palette: default_palette(),
            fixed_scale: -1.0,
            fixed_bond_length: -1.0,
            rotate: 0.0,
            add_atom_indices: false,
            add_bond_indices: false,
            isotope_labels: true,
            dummy_isotope_labels: true,
            add_stereo_annotation: false,
            show_all_cip_codes: false,
            atom_highlights_are_circles: false,
            multi_colour_highlight_style: MultiColourHighlightStyle::CircleAndLine,
            centre_molecules_before_drawing: false,
            drawing_extents_include: draw_element::ALL,
            explicit_methyl: false,
            include_radicals: true,
            include_metadata: true,
            comic_mode: false,
            variable_bond_width_multiplier: 16,
            variable_atom_radius: 0.4,
            variable_attachment_colour: c(0.8, 0.8, 0.8),
            include_chiral_flag_label: false,
            simplified_stereo_group_label: false,
            stereo_group_and_label: "and".into(),
            stereo_group_or_label: "or".into(),
            stereo_group_abs_label: "abs".into(),
            add_stereo_group_annotation: true,
            unspecified_stereo_is_unknown: false,
            single_colour_wedge_bonds: false,
            single_colour_bonds: false,
            use_mol_block_wedging: false,
            scaling_factor: 20.0,
            base_font_size: -1.0,
            draw_mols_same_scale: true,
            use_complex_query_atom_symbols: true,
            brackets_around_atom_lists: true,
            standard_colours_for_highlighted_atoms: false,
        }
    }
}

impl MolDrawOptions {
    pub fn use_default_atom_palette(&mut self) {
        self.atom_colour_palette = default_palette();
    }
    pub fn use_avalon_atom_palette(&mut self) {
        self.atom_colour_palette = avalon_palette();
    }
    pub fn use_cdk_atom_palette(&mut self) {
        self.atom_colour_palette = cdk_palette();
    }
    pub fn use_bw_atom_palette(&mut self) {
        self.atom_colour_palette = bw_palette();
    }
    /// `setDarkMode`.
    pub fn set_dark_mode(&mut self) {
        let light = DrawColour::new(0.9, 0.9, 0.9);
        self.atom_colour_palette = dark_mode_palette();
        self.background_colour = DrawColour::new(0.0, 0.0, 0.0);
        self.atom_note_colour = light;
        self.annotation_colour = light;
        self.legend_colour = light;
        self.symbol_colour = light;
        self.variable_attachment_colour = DrawColour::new(0.3, 0.3, 0.3);
    }

    /// `setMonochromeMode`.
    pub fn set_monochrome_mode(&mut self, fg: DrawColour, bg: DrawColour) {
        self.atom_colour_palette = palette_from(&[(-1, fg)]);
        self.background_colour = bg;
        self.atom_note_colour = fg;
        self.bond_note_colour = fg;
        self.annotation_colour = fg;
        self.legend_colour = fg;
        self.symbol_colour = fg;
        self.variable_attachment_colour = fg;
    }
}

/// `DrawColourToSVG`: `#RRGGBB`, or `#RRGGBBAA` when the colour is not opaque.
pub fn colour_to_svg(col: &DrawColour) -> Result<String> {
    let mut res = String::from("#");
    let has_alpha = 1.0 - col.a > 1e-3;
    let comps: &[f64] = if has_alpha {
        &[col.r, col.g, col.b, col.a]
    } else {
        &[col.r, col.g, col.b]
    };
    for &c in comps {
        // the C++ conversion goes through `unsigned`, so negatives wrap
        let v = (255.0 * c) as i32 as u32;
        if v > 255 {
            return Err(RDKitError::Value(
                "elements of the color should be between 0 and 1".into(),
            ));
        }
        res.push_str(&format!("{v:02X}"));
    }
    Ok(res)
}
