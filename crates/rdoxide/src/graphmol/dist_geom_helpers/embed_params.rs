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
    /// reads and serializes as 1 until set (upstream's default), but while
    /// [`num_threads_is_default`](Self::num_threads_is_default) holds the
    /// embedding uses the default thread count (all cores, or
    /// `RDOXIDE_NUM_THREADS`)
    pub num_threads: i32,
    pub num_threads_is_default: bool,
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

impl EmbedParameters {
    /// The thread target handed to `num_threads_to_use`.
    pub fn thread_target(&self) -> i32 {
        if self.num_threads_is_default {
            0
        } else {
            self.num_threads
        }
    }

    /// Sets an explicit thread count.
    pub fn set_num_threads(&mut self, n: i32) {
        self.num_threads = n;
        self.num_threads_is_default = false;
    }
}

impl Default for EmbedParameters {
    fn default() -> Self {
        EmbedParameters {
            max_iterations: 0,
            num_threads: 1,
            num_threads_is_default: true,
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

impl EmbedParameters {
    /// plain distance geometry
    pub fn dg() -> Self {
        Self::default()
    }

    /// distance geometry with basic knowledge terms
    pub fn kdg() -> Self {
        EmbedParameters {
            use_basic_knowledge: true,
            et_version: 1,
            ..Self::default()
        }
    }

    /// distance geometry with experimental torsion preferences (version 1)
    pub fn etdg() -> Self {
        EmbedParameters {
            use_exp_torsion_angle_prefs: true,
            et_version: 1,
            ..Self::default()
        }
    }

    pub fn etdg_v2() -> Self {
        EmbedParameters {
            use_exp_torsion_angle_prefs: true,
            et_version: 2,
            ..Self::default()
        }
    }

    pub fn etkdg() -> Self {
        EmbedParameters {
            use_exp_torsion_angle_prefs: true,
            use_basic_knowledge: true,
            et_version: 1,
            ..Self::default()
        }
    }

    pub fn etkdg_v2() -> Self {
        EmbedParameters {
            et_version: 2,
            ..Self::etkdg()
        }
    }

    pub fn etkdg_v3() -> Self {
        EmbedParameters {
            use_macrocycle_torsions: true,
            use_macrocycle_14_config: true,
            ..Self::etkdg_v2()
        }
    }

    pub fn sr_etkdg_v3() -> Self {
        EmbedParameters {
            use_small_ring_torsions: true,
            ..Self::etkdg_v2()
        }
    }

    pub fn etkdg_v4() -> Self {
        EmbedParameters {
            et_version: 4,
            use_small_ring_torsions: true,
            use_macrocycle_torsions: true,
            use_macrocycle_14_config: true,
            use_legacy_implementation: false,
            ..Self::etkdg()
        }
    }
}

/// The reasons an embedding attempt can fail (`EmbedFailureCauses`); the
/// failure counters of [`EmbedParameters`] are indexed by them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum EmbedFailureCause {
    InitialCoords = 0,
    FirstMinimization = 1,
    CheckTetrahedralCenters = 2,
    CheckChiralCenters = 3,
    MinimizeFourthDimension = 4,
    EtkMinimization = 5,
    FinalChiralBounds = 6,
    FinalCenterInVolume = 7,
    LinearDoubleBond = 8,
    BadDoubleBondStereo = 9,
    CheckChiralCenters2 = 10,
    ExceededTimeout = 11,
    Minimization = 12,
    KtermViolation = 13,
    Clash = 14,
}

/// the number of failure causes (`END_OF_ENUM`)
pub const NUM_EMBED_FAILURE_CAUSES: usize = 15;
