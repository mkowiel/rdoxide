//! Port of RDKit `Code/GraphMol/RGroupDecomposition`: decomposes molecules
//! into a core and its R groups.
//!
//! Not ported: enumeration of the input cores (needs
//! `MolEnumerator::enumerate`); requesting it is an error.

mod core;
mod data;
mod decomp;
mod decomp_data;
mod fingerprint_score;
mod ga;
mod params;
mod score;
mod utils;

pub use data::{RGroupData, RGroupMatch};
pub use decomp::{
    RGroupColumn, RGroupColumns, RGroupDecomposition, RGroupRow, RGroupRows,
    rgroup_decompose_columns, rgroup_decompose_rows,
};
pub use decomp_data::RGroupDecompositionProcessResult;
pub use params::{
    RGroupDecompositionParameters, rgroup_core_alignment, rgroup_labelling, rgroup_labels,
    rgroup_matching, rgroup_score,
};
pub use score::check_for_timeout;
pub use utils::{Labelling, relabel_mapped_dummies};
