//! Parameters controlling an embedding (upstream `Embedder.h`).

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::distgeom::BoundsMatrix;
use crate::geometry::Point3D;

/// How the initial coordinates are generated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InitialEmbeddingMode {
    #[default]
    DgEmbedding = 0,
    InternalCoordinateEmbedding = 1,
    RandomCoordinateEmbedding = 2,
}

/// Force field used to determine the ideal 1-2 and 1-3 distances.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EmbedFF {
    #[default]
    Uff = 0,
    Mmff = 1,
}

/// `EmbedParameters`
#[derive(Debug, Clone)]
pub struct EmbedParameters {
    pub max_iterations: u32,
    pub num_threads: i32,
    pub random_seed: i32,
    pub clear_confs: bool,
    pub use_random_coords: bool,
    pub box_size_mult: f64,
    pub rand_neg_eig: bool,
    pub num_zero_fail: u32,
    pub coord_map: Option<BTreeMap<i32, Point3D>>,
    pub optimizer_force_tol: f64,
    pub ignore_smoothing_failures: bool,
    pub enforce_chirality: bool,
    pub use_exp_torsion_angle_prefs: bool,
    pub use_basic_knowledge: bool,
    pub verbose: bool,
    pub basin_thresh: f64,
    pub prune_rms_thresh: f64,
    pub only_heavy_atoms_for_rms: bool,
    pub et_version: u32,
    pub bounds_mat: Option<Arc<BoundsMatrix>>,
    pub embed_fragments_separately: bool,
    pub use_small_ring_torsions: bool,
    pub use_macrocycle_torsions: bool,
    pub use_macrocycle_14_config: bool,
    pub timeout: u32,
    pub use_legacy_implementation: bool,
    pub initial_embedding_mode: InitialEmbeddingMode,
    pub cpci: Option<Arc<BTreeMap<(u32, u32), f64>>>,
    pub callback: Option<fn(u32)>,
    pub force_trans_amides: bool,
    pub use_symmetry_for_pruning: bool,
    pub bounds_mat_force_scaling: f64,
    pub only_initial_embedding: bool,
    pub track_failures: bool,
    pub failures: Vec<u32>,
    pub enable_sequential_random_seeds: bool,
    pub symmetrize_conjugated_terminal_groups_for_pruning: bool,
    pub embed_force_field: EmbedFF,
}

impl Default for EmbedParameters {
    fn default() -> Self {
        EmbedParameters {
            max_iterations: 0,
            num_threads: 1,
            random_seed: -1,
            clear_confs: true,
            use_random_coords: false,
            box_size_mult: 2.0,
            rand_neg_eig: true,
            num_zero_fail: 1,
            coord_map: None,
            optimizer_force_tol: 1e-3,
            ignore_smoothing_failures: false,
            enforce_chirality: true,
            use_exp_torsion_angle_prefs: false,
            use_basic_knowledge: false,
            verbose: false,
            basin_thresh: 5.0,
            prune_rms_thresh: -1.0,
            only_heavy_atoms_for_rms: true,
            et_version: 2,
            bounds_mat: None,
            embed_fragments_separately: true,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: false,
            use_macrocycle_14_config: false,
            timeout: 0,
            use_legacy_implementation: true,
            initial_embedding_mode: InitialEmbeddingMode::DgEmbedding,
            cpci: None,
            callback: None,
            force_trans_amides: true,
            use_symmetry_for_pruning: true,
            bounds_mat_force_scaling: 1.0,
            only_initial_embedding: false,
            track_failures: false,
            failures: Vec::new(),
            enable_sequential_random_seeds: false,
            symmetrize_conjugated_terminal_groups_for_pruning: true,
            embed_force_field: EmbedFF::Uff,
        }
    }
}
