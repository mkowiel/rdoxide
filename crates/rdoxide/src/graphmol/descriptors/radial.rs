//! Shared kernel of the radial-distribution style descriptors (RDF and
//! MORSE): for each sample point the kernel value of every atom pair is summed
//! with the product of the pair's weights.

use crate::Result;
use crate::graphmol::matrices::get_3d_distance_mat_cached;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::fp::fmuladd;

/// `get3DDistanceMat(mol, confId, false, true)`: upstream leaves the property
/// prefix unset, which selects `_` and so caches the matrix on the molecule
/// as `_3DDistanceMatrix_Conf<id>`.
pub(super) fn distance_mat_3d(mol: &mut ROMol, conf_id: i32) -> Result<Vec<f64>> {
    get_3d_distance_mat_cached(mol, conf_id, false, true, "_")
}

fn round3(x: f64) -> f64 {
    (1000.0 * x).round() / 1000.0
}

/// For every entry of `columns` (`None` sums the plain kernel values, `Some`
/// weights each pair with the product of its two atom weights) and every
/// sample point `0..points`: the sum over the atom pairs `j < k` of
/// `kernel(point, dist[j * n + k])`, rounded to three decimals. The result is
/// laid out column after column.
pub(super) fn pair_profiles(
    dist: &[f64],
    n: usize,
    columns: &[Option<&[f64]>],
    points: usize,
    kernel: impl Fn(usize, f64) -> f64,
) -> Vec<f64> {
    let mut out = vec![0.0; columns.len() * points];
    for point in 0..points {
        let mut sums = vec![0.0; columns.len()];
        for j in 0..n.saturating_sub(1) {
            for k in j + 1..n {
                let p = kernel(point, dist[j * n + k]);
                for (sum, column) in sums.iter_mut().zip(columns) {
                    // the weighted terms are contracted into fused multiply-adds
                    *sum = match column {
                        None => *sum + p,
                        Some(w) => fmuladd(w[j] * w[k], p, *sum),
                    };
                }
            }
        }
        for (c, sum) in sums.into_iter().enumerate() {
            out[c * points + point] = round3(sum);
        }
    }
    out
}
