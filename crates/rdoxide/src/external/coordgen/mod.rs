//! Port of coordgenlibs 3.0.2 (Schrödinger, Inc., BSD 3-clause; see `LICENSES/CoordGen-BSD-3-Clause-Schrodinger.txt`
//! and `docs/COORDGEN_PORT.md`), the 2D coordinate generator RDKit uses for `rdCoordGen`.

pub mod fragment;
pub mod fragmenter;
pub mod maths;
pub mod model;
pub mod templates;
#[cfg(test)]
mod testutil;
