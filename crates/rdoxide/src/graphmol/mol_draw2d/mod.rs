//! 2D molecule drawing (`Code/GraphMol/MolDraw2D` of RDKit).
//!
//! This is the drawing foundation: options and colours, the geometry helpers,
//! the TrueType reader that stands in for FreeType, and text layout.

pub mod contour;
pub mod details;
pub mod draw_annotation;
pub mod draw_mol;
pub mod draw_mol_bonds;
pub mod draw_mol_highlights;
pub mod draw_mol_lasso;
pub mod draw_mol_legend;
pub mod draw_mol_mch;
pub mod draw_mol_notes;
pub mod draw_mol_render;
pub mod draw_mol_sgroups;
pub mod draw_mol_utils;
pub mod draw_shape;
pub mod draw_text;
pub mod drawer;
pub mod drawer_reaction;
pub mod font;
pub mod helpers;
pub mod mol_draw2d_svg;
pub mod options_json;
pub mod string_rect;
pub mod svg_canvas;
pub mod utils;
