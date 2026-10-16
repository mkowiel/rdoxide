//! Marvin documents (`.mrv`, CML), after RDKit's `MarvinParse`.

pub mod defs;
mod emit;
pub mod reaction;
pub mod reader;
mod transform;
pub mod writer;
pub mod xml;

pub use reaction::{
    mrv_block_is_reaction, mrv_file_is_reaction, reaction_from_mrv_block, reaction_from_mrv_file,
    reaction_to_mrv_block, reaction_to_mrv_file,
};
pub use reader::{MrvParserParams, mol_from_mrv_block, mol_from_mrv_file};
pub use writer::{MrvWriterParams, mol_to_mrv_block, mol_to_mrv_file};
