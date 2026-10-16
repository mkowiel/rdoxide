//! SynthonSpaceSearch (upstream `GraphMol/SynthonSpaceSearch`): searching
//! spaces of molecules defined by synthons and the reactions joining them.
#![allow(missing_debug_implementations)]

pub mod details;
pub mod hit_set;
pub mod params;
pub mod search_results;
pub mod searcher;
pub mod space;
mod stream;
pub mod substructure;
pub mod synthon;
pub mod synthon_set;

pub use params::{
    CONNECTOR_SYMBOLS, FragSetUniquifyMode, MAX_CONNECTOR_NUM, ShapeBuildParams,
    SynthonSpaceSearchParams, ThreadMode, UserConfGenerator,
};
pub use search_results::SearchResults;
pub use space::{SynthonSpace, convert_text_to_db_file, formatted_integer_string};
pub use synthon::Synthon;
pub use synthon_set::SynthonSet;
