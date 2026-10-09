//! 3D autocorrelation descriptors (upstream `GraphMol/Descriptors/AUTOCORR3D.cpp`),
//! after the Dragon 6 `TDB` descriptors.
//!
//! For each topological distance (lag) 1 to 10 and each atom weighting the
//! value is `w^T (B_lag .* D) w / (n (n - 1))`, where `B_lag` flags the atom
//! pairs separated by `lag` bonds and `D` holds the 3D distances. The default
//! result has 8 weightings x 10 lags = 80 values; a custom atom property gives
//! 10. Every value is rounded to three decimals.
//!
//! Copyright (c) 2016, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use super::mol_data_3d as weights;
use super::radial::distance_mat_3d;
use crate::graphmol::matrices::get_distance_mat_cached;
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

pub const AUTOCORR3D_VERSION: &str = "1.0.0";

/// Number of lags (topological distances) reported.
const NUM_LAGS: usize = 10;

fn round3(x: f64) -> f64 {
    (1000.0 * x).round() / 1000.0
}

/// `w^T (B_lag .* D) w` summed over the atom pairs at topological distance
/// `lag`; NaN results are reported as 0.
fn lag_sum(topo: &[f64], dist3d: &[f64], w: &[f64], lag: usize) -> f64 {
    let n = w.len();
    let lag = lag as f64;
    let mut total = 0.0;
    for i in 0..n {
        let mut row = 0.0;
        for j in 0..n {
            if topo[i * n + j] == lag {
                row += dist3d[i * n + j] * w[j];
            }
        }
        total += w[i] * row;
    }
    if total.is_nan() { 0.0 } else { total }
}

/// The rounded, normalized value of each lag for one weighting.
fn profile(topo: &[f64], dist3d: &[f64], w: &[f64]) -> impl Iterator<Item = f64> {
    let n = w.len() as i32;
    // upstream divides by the integer pair count
    let pairs = f64::from(n.wrapping_mul(n - 1));
    (1..=NUM_LAGS).map(move |lag| round3(lag_sum(topo, dist3d, w, lag) / pairs))
}

/// `AUTOCORR3D(mol, res, confId, customAtomPropName)`.
pub fn autocorr3d(mol: &mut ROMol, conf_id: i32, custom_atom_prop_name: &str) -> Result<Vec<f64>> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    let topo = get_distance_mat_cached(mol, false, false, false, "")?;
    let dist3d = distance_mat_3d(mol, conf_id)?;
    let columns = if custom_atom_prop_name.is_empty() {
        vec![
            weights::unit_weights(mol),
            weights::relative_mw(mol),
            weights::relative_vdw(mol),
            weights::relative_eneg(mol),
            weights::relative_pol(mol),
            weights::relative_ion_pol(mol),
            weights::intrinsic_state(mol)?,
            weights::relative_rcov(mol),
        ]
    } else {
        vec![weights::custom_atom_prop(mol, custom_atom_prop_name)?]
    };
    Ok(columns
        .iter()
        .flat_map(|w| profile(&topo, &dist3d, w))
        .collect())
}
