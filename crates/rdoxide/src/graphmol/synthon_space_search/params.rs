//! `SynthonSpaceSearchHelpers.h`: the search and shape-build parameters.

use std::sync::Arc;

use crate::graphmol::enumerate_stereoisomers::StereoEnumerationOptions;
use crate::graphmol::gaussian_shape::{ShapeInput, ShapeOverlayOptions};
use crate::graphmol::romol::RWMol;

/// The atoms that stand for the connectors in the text format.
pub const CONNECTOR_SYMBOLS: [&str; 4] = ["[U]", "[Np]", "[Pu]", "[Am]"];
pub const MAX_CONNECTOR_NUM: usize = 4;

/// User supplied conformer generator: takes a SMILES and the maximum number
/// of conformers and returns a multi-conformer molecule (`None` on failure).
pub type UserConfGenerator = Arc<dyn Fn(&str, u32) -> Option<RWMol> + Send + Sync>;

/// How the threads are used in a search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThreadMode {
    ThreadReactions,
    ThreadFragments,
}

/// How the fragment sets of a query are made unique.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FragSetUniquifyMode {
    BySmiles,
    ByAtoms,
}

/// `SynthonSpaceSearchParams`
#[derive(Clone)]
pub struct SynthonSpaceSearchParams {
    pub max_hits: i64,
    pub max_num_frag_sets: u64,
    pub to_try_chunk_size: u64,
    pub hit_start: i64,
    pub random_sample: bool,
    pub random_seed: i32,
    pub build_hits: bool,
    pub similarity_cutoff: f64,
    pub frag_similarity_adjuster: f64,
    pub approx_similarity_adjuster: f64,
    pub min_hit_heavy_atoms: u32,
    pub max_hit_heavy_atoms: i32,
    pub min_hit_mol_wt: f64,
    pub max_hit_mol_wt: f64,
    pub min_hit_chiral_atoms: u32,
    pub max_hit_chiral_atoms: i32,
    pub num_conformers: u32,
    pub conf_rms_threshold: f64,
    pub shape_overlay_options: ShapeOverlayOptions,
    pub best_hit: bool,
    pub enumerate_unspecified_stereo: bool,
    pub stereo_enum_opts: StereoEnumerationOptions,
    pub time_out: u64,
    pub num_threads: i32,
    pub use_progress_bar: u32,
    pub user_conformer_generator: Option<UserConfGenerator>,
    pub excluded_volume: Option<Arc<ShapeInput>>,
    pub max_excluded_volume: f64,
    pub max_mean_excluded_volume: f64,
    pub possible_hits_file: String,
    pub max_possible_hits_to_write: u64,
    pub write_possible_hits_and_stop: bool,
}

impl Default for SynthonSpaceSearchParams {
    fn default() -> Self {
        SynthonSpaceSearchParams {
            max_hits: 1000,
            max_num_frag_sets: 100000,
            to_try_chunk_size: 2500000,
            hit_start: 0,
            random_sample: false,
            random_seed: -1,
            build_hits: true,
            similarity_cutoff: 0.5,
            frag_similarity_adjuster: 0.1,
            approx_similarity_adjuster: 0.1,
            min_hit_heavy_atoms: 0,
            max_hit_heavy_atoms: -1,
            min_hit_mol_wt: 0.0,
            max_hit_mol_wt: -1.0,
            min_hit_chiral_atoms: 0,
            max_hit_chiral_atoms: -1,
            num_conformers: 100,
            conf_rms_threshold: 0.5,
            shape_overlay_options: ShapeOverlayOptions::default(),
            best_hit: false,
            enumerate_unspecified_stereo: false,
            stereo_enum_opts: StereoEnumerationOptions {
                try_embedding: true,
                only_unassigned: true,
                only_stereo_groups: false,
                unique: true,
                max_isomers: 0,
                random_seed: 0xdac,
            },
            time_out: 600,
            num_threads: 0,
            use_progress_bar: 0,
            user_conformer_generator: None,
            excluded_volume: None,
            max_excluded_volume: -1.0,
            max_mean_excluded_volume: -1.0,
            possible_hits_file: String::new(),
            max_possible_hits_to_write: 10_000_000,
            write_possible_hits_and_stop: false,
        }
    }
}

/// `ShapeBuildParams`
#[derive(Clone)]
pub struct ShapeBuildParams {
    pub num_confs: u32,
    pub rms_threshold: f64,
    pub shape_sim_threshold: f64,
    pub num_threads: i32,
    pub random_seed: i32,
    pub stereo_enum_opts: StereoEnumerationOptions,
    pub use_progress_bar: u32,
    pub max_synthon_atoms: u32,
    pub max_embed_attempts: u32,
    pub time_out: u32,
    pub interim_file: String,
    pub interim_writes: u64,
    pub user_conformer_generator: Option<UserConfGenerator>,
}

impl Default for ShapeBuildParams {
    fn default() -> Self {
        ShapeBuildParams {
            num_confs: 100,
            rms_threshold: 0.5,
            shape_sim_threshold: 0.95,
            num_threads: 0,
            random_seed: 0xdac,
            stereo_enum_opts: StereoEnumerationOptions {
                try_embedding: true,
                only_unassigned: true,
                only_stereo_groups: false,
                unique: true,
                max_isomers: 0,
                random_seed: 0xdac,
            },
            use_progress_bar: 0,
            max_synthon_atoms: 0,
            max_embed_attempts: 10,
            time_out: 60,
            interim_file: String::new(),
            interim_writes: 1000,
            user_conformer_generator: None,
        }
    }
}
