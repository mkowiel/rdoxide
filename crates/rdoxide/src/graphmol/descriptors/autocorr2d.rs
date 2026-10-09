//! Broto-Moreau 2D autocorrelation descriptors (upstream
//! `GraphMol/Descriptors/AUTOCORR2D.cpp`), after the Dragon TDB names.
//!
//! For topological distances 1 to 8 and each atom weighting, four families are
//! computed: plain (`ATS`, log1p of the weight products), centered (`ATSC`-like,
//! absolute deviations), Moran (`MATS`) and Geary (`GATS`). The default result
//! has 6 weightings x 4 families x 8 distances = 192 values; a custom atom
//! property gives 4 x 8 = 32. Every value is rounded to three decimals.
//!
//! Copyright (c) 2016, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use super::mol_data_3d as weights;
use crate::Result;
use crate::graphmol::matrices::get_distance_mat_cached;
use crate::graphmol::romol::ROMol;

pub const AUTOCORR2D_VERSION: &str = "1.0.0";

/// Largest topological distance considered.
const MAX_LAG: usize = 8;
/// Number of built-in atom weightings.
const NUM_WEIGHTS: usize = 6;

/// Per-lag sums of one weighting.
#[derive(Default, Clone, Copy)]
struct LagSums {
    /// sum of w_i * w_j
    plain: f64,
    /// sum of |w_i - mean| * |w_j - mean|
    centered: f64,
    /// sum of (w_i - mean) * (w_j - mean)
    moran: f64,
    /// sum of the squared differences (w_i - w_j)^2
    geary: f64,
}

/// The four rounded families of one lag, or zeros without atom pairs.
struct LagValues {
    plain: f64,
    centered: f64,
    moran: f64,
    geary: f64,
}

fn round3(x: f64) -> f64 {
    (1000.0 * x).round() / 1000.0
}

/// Mean and sum of squared deviations of `w`.
fn mean_and_spread(w: &[f64]) -> (f64, f64) {
    let n = w.len() as f64;
    let mean = w.iter().map(|x| x / n).sum::<f64>();
    let spread = w.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>();
    (mean, spread)
}

/// Sums over the atom pairs at topological distance `lag` and their count.
/// `geary_second` selects the weight of the second atom in the Geary term:
/// upstream's custom-property variant subtracts a zero-filled array there
/// instead of the weights.
fn lag_sums(
    dist: &[f64],
    w: &[f64],
    mean: f64,
    lag: usize,
    geary_second: impl Fn(usize) -> f64,
) -> (LagSums, u32) {
    let n = w.len();
    let mut sums = LagSums::default();
    let mut pairs = 0;
    for i in 0..n {
        for j in i + 1..n {
            // the matrix is read as dist[j * n + i], as upstream does
            if dist[j * n + i] == (lag + 1) as f64 {
                sums.moran += (w[i] - mean) * (w[j] - mean);
                sums.geary += (w[i] - w[j]) * (w[i] - geary_second(j));
                sums.plain += w[i] * w[j];
                sums.centered += (w[i] - mean).abs() * (w[j] - mean).abs();
                pairs += 1;
            }
        }
    }
    (sums, pairs)
}

/// Normalizes the sums of one lag; `valid` is the condition under which the
/// lag has a value, otherwise all four families are zero.
fn lag_values(sums: LagSums, spread: f64, pairs: u32, n_atoms: usize, valid: bool) -> LagValues {
    if !valid {
        return LagValues {
            plain: 0.0,
            centered: 0.0,
            moran: 0.0,
            geary: 0.0,
        };
    }
    let n = n_atoms as f64;
    LagValues {
        plain: round3(sums.plain.ln_1p()),
        centered: round3(sums.centered),
        moran: round3(sums.moran / spread / f64::from(pairs) * n),
        geary: round3(sums.geary / spread / f64::from(pairs) * (n - 1.0) / 2.0),
    }
}

/// The 192 values for the six built-in weightings: mass, van der Waals
/// volume, electronegativity, polarizability, ionization potential and
/// intrinsic state, each with the families plain, centered, Moran, Geary.
fn builtin_weights(mol: &ROMol, dist: &[f64]) -> Result<Vec<f64>> {
    let n_atoms = mol.num_atoms() as usize;
    let per_weight = [
        weights::relative_mw(mol),
        weights::relative_vdw(mol),
        weights::relative_eneg(mol),
        weights::relative_pol(mol),
        weights::relative_ion_pol(mol),
        weights::intrinsic_state(mol)?,
    ];
    let mut res = vec![0.0; 4 * NUM_WEIGHTS * MAX_LAG];
    for (t, w) in per_weight.iter().enumerate() {
        let (mean, spread) = mean_and_spread(w);
        for k in 0..MAX_LAG {
            let (sums, pairs) = lag_sums(dist, w, mean, k, |j| w[j]);
            let v = lag_values(sums, spread, pairs, n_atoms, pairs > 0 && spread > 0.0);
            let at = t * MAX_LAG + k;
            res[at] = v.plain;
            res[at + NUM_WEIGHTS * MAX_LAG] = v.centered;
            res[at + 2 * NUM_WEIGHTS * MAX_LAG] = v.moran;
            res[at + 3 * NUM_WEIGHTS * MAX_LAG] = v.geary;
        }
    }
    Ok(res)
}

/// The 32 values for one custom atom property.
fn custom_weights(mol: &ROMol, dist: &[f64], prop: &str) -> Result<Vec<f64>> {
    let n_atoms = mol.num_atoms() as usize;
    let w = weights::custom_atom_prop(mol, prop)?;
    let (mean, spread) = mean_and_spread(&w);
    let mut res = vec![0.0; 4 * MAX_LAG];
    for k in 0..MAX_LAG {
        let (sums, pairs) = lag_sums(dist, &w, mean, k, |_| 0.0);
        let v = lag_values(sums, spread, pairs, n_atoms, pairs > 0);
        res[k] = v.plain;
        res[k + MAX_LAG] = v.centered;
        res[k + 2 * MAX_LAG] = v.moran;
        res[k + 3 * MAX_LAG] = v.geary;
    }
    Ok(res)
}

/// `AUTOCORR2D(mol, result, customAtomPropName)`: 192 values, or 32 for a
/// non-empty custom atom property name.
pub fn autocorr2d(mol: &mut ROMol, custom_atom_prop_name: &str) -> Result<Vec<f64>> {
    let dist = get_distance_mat_cached(mol, false, false, false, "")?;
    if custom_atom_prop_name.is_empty() {
        builtin_weights(mol, &dist)
    } else {
        custom_weights(mol, &dist, custom_atom_prop_name)
    }
}
