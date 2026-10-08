//! # rdoxide — RDKit, rusted
//!
//! A pure-Rust port of the [RDKit](https://www.rdkit.org) cheminformatics
//! toolkit (upstream version tracked in `UPSTREAM.toml`). The module tree
//! mirrors upstream's `Code/` directory so changes in new RDKit releases can be
//! mapped onto this crate file by file (see `porting/manifest.toml`).
//!
//! Portions of this crate are translations of RDKit C++ sources, Copyright (C)
//! Greg Landrum, Rational Discovery LLC, Novartis and other RDKit
//! contributors, used under the BSD 3-Clause license (see `LICENSES/`).

pub mod datastructs;
pub mod error;
pub mod external;
pub mod geometry;
pub mod graphmol;
pub mod ml;
pub mod rdgeneral;

pub use error::{RDKitError, Result, SanitizeError};
