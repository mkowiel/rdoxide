//! Internal coordinates keyed by the unified ids of the involved bonds
//! (upstream `ZMatrixBuilder.h`).

use std::collections::HashMap;

use crate::distgeom::TorsionCandidates;

/// Bond lengths by bond index, and bond angles and torsion candidates by the
/// unified id of the bond pair or triple.
#[derive(Debug, Clone, Default)]
pub struct InternalCoordinates {
    pub torsion_range: HashMap<usize, TorsionCandidates>,
    pub angles: HashMap<usize, f64>,
    pub lengths: Vec<f64>,
}

impl InternalCoordinates {
    pub fn new(num_bonds: usize) -> Self {
        InternalCoordinates {
            torsion_range: HashMap::new(),
            angles: HashMap::new(),
            lengths: vec![0.0; num_bonds],
        }
    }
}
