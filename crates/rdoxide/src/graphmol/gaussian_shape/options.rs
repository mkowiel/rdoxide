//! `ShapeOverlayOptions.h`: the options controlling an overlay.

use std::fmt;

/// How the optimisation is started.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StartMode {
    /// No rotation, just normalization if requested
    Rotate0,
    /// `Rotate0` plus rotate by 180 degrees about each of x, y, z
    Rotate180,
    /// `Rotate180` plus, as the PubChem code does, rotate by a small amount
    /// (~25 degrees) about each axis and use the highest scoring of those as
    /// the start point for that rotation angle
    Rotate180Wiggle,
    /// `Rotate180` plus rotate by 45 degrees about pairs of each of x, y, z
    Rotate45,
    /// No rotation, translate probe to each end of ref
    Rotate0Fragment,
    /// Translate probe to each end of ref and then `Rotate180`
    Rotate180Fragment,
    /// Translate probe to each end of ref and then `Rotate45`
    Rotate45Fragment,
    /// Uses the eigenvalues of the principal vectors to decide whether to do
    /// `Rotate180Wiggle` or `Rotate45`
    ALaPubchem,
}

impl fmt::Display for StartMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            StartMode::Rotate0 => "ROTATE_0",
            StartMode::Rotate0Fragment => "ROTATE_0_FRAGMENT",
            StartMode::Rotate180 => "ROTATE_180",
            StartMode::Rotate180Fragment => "ROTATE_180_FRAGMENT",
            StartMode::Rotate180Wiggle => "ROTATE_180_WIGGLE",
            StartMode::Rotate45 => "ROTATE_45",
            StartMode::Rotate45Fragment => "ROTATE_45_FRAGMENT",
            StartMode::ALaPubchem => "A_LA_PUBCHEM",
        })
    }
}

/// What drives the optimisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptimMode {
    /// Drive the optimisation by shape overlap only.
    ShapeOnly,
    /// Drive the optimisation by shape, but include colour in the score to
    /// determine the best solution. Colour is never used in the optimisation
    /// stage.
    ShapePlusColorScore,
    /// Drive the optimisation by overlap of shape and colour features.
    ShapePlusColor,
}

/// `ShapeOverlayOptions`
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeOverlayOptions {
    /// Different modes for starting the optimisation. The default is as used
    /// by the PubChem code.
    pub start_mode: StartMode,
    pub optim_mode: OptimMode,
    /// When doing a Tversky similarity, the alpha value. If alpha and beta are
    /// both the default 1.0, it's a Tanimoto similarity. Tversky is
    /// `O / (A * (R - O) + B * (F - O) + O)` where `O` is the overlap volume,
    /// `R` is the reference's volume and `F` is the fit's volume.
    pub sim_alpha: f64,
    /// When doing a Tversky similarity, the beta value.
    pub sim_beta: f64,
    /// If using colours, the relative weight of the colour score.
    pub opt_param: f64,
    /// Maximum number of steps for the optimiser to take.
    pub n_steps: u32,
    /// Whether to put the shapes into their canonical conformations (centred
    /// at the origin, aligned along the principal axes) before starting.
    pub normalize: bool,
    /// Whether to use a distance cutoff for the volume calculation.
    pub use_dist_cutoff: bool,
    /// If 2 atoms are more than this distance apart, they are not included in
    /// the volume calculation.
    pub dist_cutoff: f64,
    /// Optimisation stops when the shape score changes by less than this.
    pub shape_convergence_criterion: f64,
}

impl Default for ShapeOverlayOptions {
    fn default() -> Self {
        ShapeOverlayOptions {
            start_mode: StartMode::ALaPubchem,
            optim_mode: OptimMode::ShapePlusColorScore,
            sim_alpha: 1.0,
            sim_beta: 1.0,
            opt_param: 0.5,
            n_steps: 100,
            normalize: true,
            use_dist_cutoff: true,
            dist_cutoff: 4.5,
            shape_convergence_criterion: 0.001,
        }
    }
}
