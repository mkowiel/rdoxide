//! Hooks into perception code that is ported in stage 3 (ring finding,
//! sanitization). Until then these report the upstream precondition failure.

use rdoxide::graphmol::ROMol;
use rdoxide::{RDKitError, Result};

/// `MolOps::findSSSR` as called lazily by `Atom.IsInRing()` & co.
pub fn find_sssr_if_needed(mol: &mut ROMol) -> Result<()> {
    if mol.ring_info().is_sssr_or_better() {
        return Ok(());
    }
    Err(RDKitError::precondition("RingInfo not initialized"))
}
