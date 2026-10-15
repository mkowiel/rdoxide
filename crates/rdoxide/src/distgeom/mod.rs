//! Distance geometry (upstream `Code/DistGeom`).

mod bounds_matrix;
mod chiral_set;
mod contribs;
mod triangle_smooth;

pub use bounds_matrix::BoundsMatrix;
pub use chiral_set::{ChiralSet, ChiralSetStructureFlags};
pub use contribs::{
    ChiralViolationContribs, DistViolationContribs, FourthDimContribs, calc_chiral_volume,
};
pub use triangle_smooth::triangle_smooth_bounds;
