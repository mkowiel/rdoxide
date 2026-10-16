//! RASCAL maximum common edge substructure (MCES) and the clustering built
//! on it.
//!
//! Ported from RDKit `Code/GraphMol/RascalMCES` (2026.09.1), which implements
//! Raymond, Gardiner & Willett, "RASCAL: Calculation of Graph Similarity using
//! Maximum Common Edge Subgraphs", Comput. J. 45, 631 (2002).
//!
//! Copyright (c) 2023, David Cosgrove and other RDKit contributors. Licensed
//! under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt). The
//! linear assignment solver is translated from SciPy (see
//! LICENSES/SciPy-BSD-3-Clause.txt).

mod cluster;
mod lap;
mod mces;
mod options;
mod partition_set;
mod result;

pub use cluster::{rascal_butina_cluster, rascal_cluster};
pub use mces::rascal_mces;
pub use options::{RascalClusterOptions, RascalOptions};
pub use result::RascalResult;

/// Pieces of the implementation that the tests and the clustering reach into
/// (upstream's `details` namespace).
pub mod details {
    pub use super::cluster::{ClusNode, build_proximity_graph};
    pub use super::mces::{get_bond_labels, tier1_sim, tier2_sim};
    pub use super::result::{clean_smarts, johnson_similarity, result_compare};
}
