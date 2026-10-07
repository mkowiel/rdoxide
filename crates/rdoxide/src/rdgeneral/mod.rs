//! Port of RDKit `Code/RDGeneral`: properties, errors, logging, versions.

pub mod dict;
pub mod log;
pub mod random;
pub mod strings;

pub use dict::{Dict, RDProps, RDValue};

/// Upstream RDKit version this crate is a port of (`rdBase.rdkitVersion`).
pub const RDKIT_VERSION: &str = "2026.09.1";

/// Version of rdoxide itself.
pub const RDOXIDE_VERSION: &str = env!("CARGO_PKG_VERSION");
