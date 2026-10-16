//! Gaussian shape overlay: molecules are represented as sums of atom-centred
//! Gaussians (plus optional pharmacophore "colour" features) and overlaid by
//! maximising their Tversky similarity.
//!
//! Ported from RDKit `Code/GraphMol/GaussianShape` (2026.09.1), itself based
//! on the PubChem shape-alignment code and Grant, Gallardo & Pickup, J. Comput.
//! Chem. 17, 1653 (1996).
//!
//! Copyright (c) 2025-2026, Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-BSD-3-Clause.txt).

mod alignment;
mod options;
mod overlay;
mod shape_input;

pub use alignment::{SingleConformerAlignment, calc_color_vol_and_grads, calc_shape_vol_and_grads};
pub use options::{OptimMode, ShapeOverlayOptions, StartMode};
pub use overlay::{
    AllConformerScores, align_molecule, align_molecule_to_shape, align_shape, score_molecule,
    score_molecule_all_conformers, score_molecule_to_shape, score_shape,
};
pub use shape_input::{
    CARBON_RAD, CustomFeature, DUMMY_RAD, KAPPA, ShapeInput, ShapeInputOptions, ShapePair,
    apply_transform_to_shape, compute_feature_pos, find_features, max_score, translate_shape,
};
