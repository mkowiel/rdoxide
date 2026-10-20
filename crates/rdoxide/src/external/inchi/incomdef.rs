//! Constants of `INCHI_BASE/src/incomdef.h` (added as the ported code needs them).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

/// number of hydrogen isotopes: protium, deuterium, tritium
pub const NUM_H_ISOTOPES: usize = 3;

pub const RADICAL_SINGLET: i32 = 1;
pub const RADICAL_DOUBLET: i32 = 2;
pub const RADICAL_TRIPLET: i32 = 3;

/// element: lowest valence
pub const METAL: i32 = 1;
/// element: lowest and next to it valence
pub const METAL2: i32 = 3;
/// metal bitmap
pub const IS_METAL: i32 = 3;
