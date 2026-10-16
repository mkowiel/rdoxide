//! Marvin documents (`.mrv`, CML), after RDKit's `MarvinParse`.

pub mod defs;
mod emit;
pub mod reader;
mod transform;
pub mod writer;
pub mod xml;

pub use reader::{MrvParserParams, mol_from_mrv_block, mol_from_mrv_file};
pub use writer::{MrvWriterParams, mol_to_mrv_block, mol_to_mrv_file};
