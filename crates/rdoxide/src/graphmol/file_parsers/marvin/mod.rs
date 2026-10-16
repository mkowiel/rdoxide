//! Marvin documents (`.mrv`, CML), after RDKit's `MarvinParse`.

pub mod defs;
pub mod reader;
mod transform;
pub mod xml;

pub use reader::{MrvParserParams, mol_from_mrv_block, mol_from_mrv_file};
