//! 2D molecule drawing (`Code/GraphMol/MolDraw2D` of RDKit).
//!
//! This is the drawing foundation: options and colours, the geometry helpers,
//! the TrueType reader that stands in for FreeType, and text layout.

pub mod details;
pub mod draw_annotation;
pub mod draw_mol;
pub mod draw_mol_bonds;
pub mod draw_mol_highlights;
pub mod draw_mol_legend;
pub mod draw_mol_notes;
pub mod draw_mol_sgroups;
pub mod draw_mol_utils;
pub mod draw_shape;
pub mod draw_text;
pub mod font;
pub mod helpers;
pub mod string_rect;
pub mod utils;
