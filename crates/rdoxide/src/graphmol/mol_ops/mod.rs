//! Port of RDKit `Code/GraphMol/MolOps.{h,cpp}` and friends (sanitization,
//! H handling, ring finding entry points).

pub use crate::graphmol::chirality::{clear_dir_flags, clear_single_bond_dir_flags};

use crate::graphmol::romol::RWMol;
use crate::{RDKitError, Result};

/// The sanitize / removeHs / stereo steps `MolFromSmiles` runs after
/// parsing. Filled in by the stage-3 perception port.
pub fn smiles_post_parse_sanitize(
    _mol: &mut RWMol,
    _sanitize: bool,
    _remove_hs: bool,
) -> Result<()> {
    Err(RDKitError::Value("sanitization is not ported yet".into()))
}
