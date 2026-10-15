//! Port of RDKit `Code/GraphMol/FMCS`: the maximum common substructure of a
//! set of molecules.

mod cache;
mod final_match;
mod graph;
mod mcs;
mod params;
mod seed;
mod target;

pub use final_match::{final_chirality_check, final_match_check};
pub use graph::{Graph, NOT_SET};
pub use mcs::{MaximumCommonSubgraph, find_mcs, find_mcs_json, find_mcs_with_options};
pub use params::*;
pub use seed::{MatchV, Seed};
