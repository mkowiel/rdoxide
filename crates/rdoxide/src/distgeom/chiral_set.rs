//! A quartet of points with bounds on its chiral volume (upstream
//! `DistGeom/ChiralSet.h`).

use crate::error::{RDKitError, Result};

/// Special cases that may make embedding more difficult.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u64)]
pub enum ChiralSetStructureFlags {
    /// a chiral center involved in fusing two or more small rings
    InFusedSmallRings = 1 << 0,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChiralSet {
    /// the center
    pub idx0: usize,
    pub idx1: usize,
    pub idx2: usize,
    pub idx3: usize,
    pub idx4: usize,
    pub volume_lower_bound: f64,
    pub volume_upper_bound: f64,
    pub structure_flags: u64,
}

impl ChiralSet {
    pub fn new(
        pids: [usize; 5],
        lower_vol_bound: f64,
        upper_vol_bound: f64,
        structure_flags: u64,
    ) -> Result<Self> {
        if lower_vol_bound > upper_vol_bound {
            return Err(RDKitError::invariant("Inconsistent bounds\n"));
        }
        Ok(ChiralSet {
            idx0: pids[0],
            idx1: pids[1],
            idx2: pids[2],
            idx3: pids[3],
            idx4: pids[4],
            volume_lower_bound: lower_vol_bound,
            volume_upper_bound: upper_vol_bound,
            structure_flags,
        })
    }
}
