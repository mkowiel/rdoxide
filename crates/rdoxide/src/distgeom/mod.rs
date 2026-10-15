//! Distance geometry (upstream `Code/DistGeom`).

mod bounds_matrix;
mod chiral_set;
mod contribs;
mod etkdg_ff;
mod triangle_smooth;
mod utils;
mod zmatrix;

pub use bounds_matrix::BoundsMatrix;
pub use chiral_set::{ChiralSet, ChiralSetStructureFlags};
pub use contribs::{
    ChiralViolationContribs, DistViolationContribs, FourthDimContribs, calc_chiral_volume,
};
pub use etkdg_ff::{
    ChargeProducts, add_torsion_terms, construct_3d_force_field, construct_3d_improper_force_field,
    construct_all_in_one_force_field, construct_plain_3d_force_field, improper_term_energies,
};
pub use triangle_smooth::triangle_smooth_bounds;
pub use utils::{
    compute_initial_coords, compute_random_coords, construct_force_field, pick_random_dist_mat,
};
pub use zmatrix::{
    InternalCoordinateReference, TorsionCandidates, TorsionDependence, TorsionRange, TorsionValues,
    ZMatrix, ZMatrixRow, compute_z_matrix_coords, quantize,
};
