//! GEometry, Topology and Atom-Weights AssemblY (GETAWAY) descriptors
//! (upstream `GraphMol/Descriptors/GETAWAY.cpp`), after Consonni, Todeschini
//! and Pavan, J. Chem. Inf. Comput. Sci. 2002, 42, 682.
//!
//! The influence matrix `H = X (X^T X)^+ X^T` of the centered coordinates and
//! the influence/distance matrix `R_ij = sqrt(h_ii h_jj) / r_ij` are combined
//! with the topological distances and seven atom weightings (unit, mass, van
//! der Waals volume, electronegativity, polarizability, ionization potential,
//! intrinsic state): 273 values; a custom atom property gives 45. The
//! first four values are weighting independent (`ITH`, `ISH`, `HIC`, `HGM`);
//! then every weighting contributes `H0-H8`, `HT`, `HATS0-HATS8`, `HATST`,
//! followed by `RCON`, `RARS`, `REIG` and the `R1-R8`, `RT`, `R1+-R8+`, `RT+`
//! block of every weighting. Every value is rounded to three decimals.
//!
//! Copyright (c) 2016-2017, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use std::collections::VecDeque;

use super::mol_data_3d as weights;
use super::pbf::pbf;
use crate::graphmol::matrices::{
    get_3d_distance_mat_cached, get_adjacency_matrix_cached, get_distance_mat_cached,
};
use crate::graphmol::romol::ROMol;
use crate::numerics::dense::{Matrix, packet_sum};
use crate::numerics::jacobi_svd::jacobi_svd;
use crate::{RDKitError, Result};

pub const GETAWAY_VERSION: &str = "1.0.0";

/// Topological distances reported individually (`H0-H8`, `HATS0-HATS8`); the
/// longer ones only enter the totals.
const NUM_LAGS: usize = 9;
/// Singular values at or below this are treated as zero in the pseudo-inverse.
const PINV_TOLERANCE: f64 = 1.0e-3;
/// Planar molecules (plane of best fit below this) have two dimensions.
const PLANAR_PBF: f64 = 1.0e-5;

fn roundn(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

// ---------------------------------------------------------------------------
// the leverage matrix H and the influence/distance matrix R
// ---------------------------------------------------------------------------

/// Pseudo-inverse of a matrix from its SVD, dropping tiny singular values.
fn pseudo_inverse(a: &Matrix) -> Matrix {
    let svd = jacobi_svd(a, true);
    let inverse = Matrix::from_fn(a.cols(), a.cols(), |r, c| {
        let s = svd.singular_values[r];
        if r == c && s > PINV_TOLERANCE {
            1.0 / s
        } else {
            0.0
        }
    });
    let u = svd.u.as_ref().expect("U was requested");
    svd.v.mul(&inverse).mul(&u.transpose())
}

/// `H = X (X^T X)^+ X^T` for the centered coordinates `x`.
fn leverage_matrix(x: &Matrix) -> Matrix {
    let weighted = x.transpose_mul(x);
    x.mul(&pseudo_inverse(&weighted)).mul(&x.transpose())
}

/// `R_ij = sqrt(h_ii h_jj) / d_ij` off the diagonal, 0 on it.
fn influence_matrix(h: &Matrix, dist3d: &[f64]) -> Matrix {
    let n = h.rows();
    let mut r = Matrix::zeros(n, n);
    for i in 0..n.saturating_sub(1) {
        for j in i + 1..n {
            let v = (h.get(i, i) * h.get(j, j)).sqrt() / dist3d[i * n + j];
            r.set(i, j, v);
            r.set(j, i, v);
        }
    }
    r
}

/// Coordinates minus their column means, as an `n x 3` matrix.
fn centered_coordinates(points: &[crate::geometry::Point3D]) -> Matrix {
    let raw = Matrix::from_fn(points.len(), 3, |i, axis| {
        let p = points[i];
        [p.x, p.y, p.z][axis]
    });
    let mean = raw.column_means();
    Matrix::from_fn(points.len(), 3, |i, axis| raw.get(i, axis) - mean[axis])
}

// ---------------------------------------------------------------------------
// the information indices (leverage clustering)
// ---------------------------------------------------------------------------

/// `round_to_n_digits_`: `x` rounded to `n` significant digits.
fn round_to_digits(x: f64, n: u32) -> f64 {
    let digits = n.max(1) as usize;
    format!("{:.*e}", digits - 1, x).parse().unwrap_or(x)
}

/// A `double` printed by a default C++ output stream (`%g` with six
/// significant digits).
fn ostream_double(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.into();
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let sci = format!("{x:.5e}");
    let (mantissa, exponent) = sci.split_once('e').expect("scientific notation");
    let exponent: i32 = exponent.parse().expect("exponent digits");
    let strip = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if !(-4..6).contains(&exponent) {
        format!(
            "{}e{}{:02}",
            strip(mantissa),
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        )
    } else {
        strip(&format!("{:.*}", (5 - exponent) as usize, x))
    }
}

/// C's `atoi`.
fn atoi(s: &str) -> i64 {
    let s = s.trim_start();
    let (negative, digits) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };
    let value = digits
        .chars()
        .take_while(char::is_ascii_digit)
        .fold(0_i64, |acc, d| {
            acc.saturating_mul(10)
                .saturating_add(i64::from(d as u8 - b'0'))
        });
    if negative { -value } else { value }
}

fn count_leading_zeros(s: &str) -> usize {
    s.chars().take_while(|&c| c == '0').count()
}

/// `IsClose2`: whether the printed fractional digits of `a` and `b` differ by
/// at most one unit of the last digit.
fn is_close(a: f64, b: f64) -> bool {
    let fraction = |x: f64| {
        let s = ostream_double(x);
        match s.find('.') {
            Some(i) => s[i + 1..].to_string(),
            None => s,
        }
    };
    let (mut ta, mut tb) = (fraction(a), fraction(b));
    if count_leading_zeros(&ta) == count_leading_zeros(&tb) {
        if ta.len() > tb.len() {
            tb.push('0');
        }
        if ta.len() < tb.len() {
            ta.push('0');
        }
    }
    ta.len() == tb.len() && (atoi(&ta) - atoi(&tb)).abs() < 2
}

/// `clusterArray2`: sizes of the groups of nearly equal values of `data`,
/// including the way upstream drops elements after a group that reaches the
/// end of the sorted values.
fn cluster_sizes(mut data: Vec<f64>) -> Vec<f64> {
    data.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let mut queue: VecDeque<f64> = data.into();
    let mut sizes = Vec::new();
    let mut count = 0.0;
    while let Some(dk) = queue.pop_front() {
        count += 1.0;
        if queue.is_empty() {
            sizes.push(count);
        }
        for i in 0..queue.len() {
            if is_close(dk, queue[i]) {
                count += 1.0;
            } else {
                sizes.push(count);
                queue.drain(..i);
                count = 0.0;
                break;
            }
            if i == queue.len() - 1 {
                sizes.push(count);
                queue.pop_front();
                count = 0.0;
            }
        }
    }
    sizes
}

/// `ITH` and `ISH`, the information content and standardized information
/// content of the leverages of the heavy atoms.
fn information_indices(h: &Matrix, heavy: &[bool], precision: u32) -> (f64, f64) {
    let leverages: Vec<f64> = (0..h.rows())
        .filter(|&i| heavy[i])
        .map(|i| round_to_digits(h.get(i, i), precision))
        .collect();
    let sizes = cluster_sizes(leverages.clone());
    let num_heavy = leverages.len() as f64;
    let ith0 = num_heavy * num_heavy.ln() / 2.0_f64.ln();
    let ith = sizes
        .iter()
        .fold(ith0, |acc, c| acc - c * c.ln() / 2.0_f64.ln());
    (ith, ith / ith0)
}

// ---------------------------------------------------------------------------
// the weighted sums over atom pairs
// ---------------------------------------------------------------------------

/// `getHATS`
fn hats(w1: f64, w2: f64, h1: f64, h2: f64) -> f64 {
    w1 * h1 * w2 * h2
}

/// `getH`
fn weighted(w1: f64, w2: f64, v: f64) -> f64 {
    w1 * v * w2
}

/// Atom pairs `j < k` by topological distance: the pairs of each lag below
/// `NUM_LAGS` (index 0 is unused) and, in order of increasing distance, the
/// pairs at longer distances.
struct LagPairs {
    short: Vec<Vec<(usize, usize)>>,
    long: Vec<(usize, usize)>,
}

fn lag_pairs(topological: &[f64], n: usize) -> LagPairs {
    let mut by_distance: Vec<(f64, usize, usize)> = Vec::new();
    for j in 0..n.saturating_sub(1) {
        for k in j + 1..n {
            by_distance.push((topological[j * n + k], j, k));
        }
    }
    // stable: the pairs of one distance stay in row-major order
    by_distance.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut pairs = LagPairs {
        short: vec![Vec::new(); NUM_LAGS],
        long: Vec::new(),
    };
    for (d, j, k) in by_distance {
        if d >= 1.0 && d < NUM_LAGS as f64 && d.fract() == 0.0 {
            pairs.short[d as usize].push((j, k));
        } else if d >= NUM_LAGS as f64 && d.fract() == 0.0 {
            pairs.long.push((j, k));
        }
    }
    pairs
}

/// The sums of one atom weighting.
#[derive(Default)]
struct WeightedSums {
    /// `H0-H8`
    hk: [f64; NUM_LAGS],
    /// `HATS0-HATS8`
    hats_k: [f64; NUM_LAGS],
    /// `HT`
    ht: f64,
    /// `HATST`
    hats_t: f64,
    /// `R1-R8`
    rk: [f64; NUM_LAGS - 1],
    /// `R1+-R8+`
    rp: [f64; NUM_LAGS - 1],
    /// `RT`
    rt: f64,
    /// `RT+`
    rtm: f64,
}

fn weighted_sums(w: &[f64], h: &Matrix, r: &Matrix, pairs: &LagPairs) -> WeightedSums {
    let n = h.rows();
    let mut s = WeightedSums::default();

    // lag 0: the atoms themselves
    for j in 0..n {
        let hjj = h.get(j, j);
        if hjj > 0.0 {
            s.hats_k[0] += hats(w[j], w[j], hjj, hjj);
            s.hk[0] += weighted(w[j], w[j], hjj);
        }
    }
    // lags 1-8: sums per lag
    for lag in 1..NUM_LAGS {
        let (mut hats_sum, mut h_sum, mut r_sum, mut r_max) = (0.0, 0.0, 0.0, 0.0);
        for &(j, k) in &pairs.short[lag] {
            hats_sum += hats(w[j], w[k], h.get(j, j), h.get(k, k));
            if h.get(j, k) > 0.0 {
                h_sum += weighted(w[j], w[k], h.get(j, k));
            }
            let tmp = weighted(w[j], w[k], r.get(j, k));
            r_sum += tmp;
            if tmp > r_max {
                r_max = tmp;
            }
        }
        s.hats_k[lag] = hats_sum;
        s.hk[lag] = h_sum;
        s.rk[lag - 1] = r_sum;
        s.rp[lag - 1] = r_max;
    }
    // longer distances only enter the totals
    let (mut hats_long, mut h_long, mut r_long) = (0.0, 0.0, 0.0);
    for &(j, k) in &pairs.long {
        hats_long += 2.0 * hats(w[j], w[k], h.get(j, j), h.get(k, k));
        if h.get(j, k) > 0.0 {
            h_long += 2.0 * weighted(w[j], w[k], h.get(j, k));
        }
        r_long += 2.0 * weighted(w[j], w[k], r.get(j, k));
    }

    s.hats_t = s.hats_k[0] + hats_long;
    s.ht = s.hk[0] + h_long;
    for lag in 1..NUM_LAGS {
        s.hats_t += 2.0 * s.hats_k[lag];
        s.ht += 2.0 * s.hk[lag];
    }
    let rt = s.rk.iter().fold(0.0, |acc, v| acc + 2.0 * v);
    s.rt = rt + r_long;
    s.rtm =
        s.rp.iter()
            .fold(0.0, |acc, &v| if v > acc { v } else { acc });
    s
}

// ---------------------------------------------------------------------------
// the descriptor
// ---------------------------------------------------------------------------

/// `GETAWAY(mol, res, confId, precision, customAtomPropName)`: 273 values, or
/// 45 for a custom atom property. `precision` is the number of significant
/// digits of the leverages when they are clustered.
pub fn getaway(
    mol: &mut ROMol,
    conf_id: i32,
    precision: u32,
    custom_atom_prop_name: &str,
) -> Result<Vec<f64>> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    let n = mol.num_atoms() as usize;
    let points = mol.conformer(conf_id)?.positions().to_vec();
    let heavy: Vec<bool> = mol.atoms().iter().map(|a| a.atomic_num() > 1).collect();
    let dist3d = get_3d_distance_mat_cached(mol, conf_id, false, false, "_")?;
    let adjacency = get_adjacency_matrix_cached(mol, false, 0, false, "")?;

    let x = centered_coordinates(&points);
    let h = leverage_matrix(&x);
    let r = influence_matrix(&h, &dist3d);

    let (ith, ish) = information_indices(&h, &heavy, precision);
    let dimension = if pbf(mol, -1)? < PLANAR_PBF { 2.0 } else { 3.0 };
    let hic = (0..n).fold(0.0, |acc, i| {
        let v = h.get(i, i) / dimension;
        acc - v * v.ln() / 2.0_f64.ln()
    });
    let product = (0..n).fold(1.0, |acc, i| acc * h.get(i, i));
    let hgm = 100.0 * product.powf(1.0 / n as f64);

    let row_sums = r.row_sums();
    let rars = packet_sum(&row_sums, 0) / n as f64;
    let reig = jacobi_svd(&r, false).singular_values[0];
    let mut rcon = 0.0;
    for i in 0..n.saturating_sub(1) {
        for j in i + 1..n {
            if adjacency[i * n + j] > 0.0 {
                rcon += (row_sums[i] * row_sums[j]).sqrt();
            }
        }
    }

    let topological = get_distance_mat_cached(mol, false, false, false, "")?;
    let pairs = lag_pairs(&topological, n);
    let columns = if custom_atom_prop_name.is_empty() {
        vec![
            weights::unit_weights(mol),
            weights::relative_mw(mol),
            weights::relative_vdw(mol),
            weights::relative_eneg(mol),
            weights::relative_pol(mol),
            weights::relative_ion_pol(mol),
            weights::intrinsic_state(mol)?,
        ]
    } else {
        vec![weights::custom_atom_prop(mol, custom_atom_prop_name)?]
    };
    let sums: Vec<WeightedSums> = columns
        .iter()
        .map(|w| weighted_sums(w, &h, &r, &pairs))
        .collect();

    let mut res = vec![ith, ish, hic, hgm];
    for s in &sums {
        res.extend(s.hk);
        res.push(s.ht);
        res.extend(s.hats_k);
        res.push(s.hats_t);
    }
    res.extend([rcon, rars, reig]);
    for s in &sums {
        res.extend(s.rk);
        res.push(s.rt);
        res.extend(s.rp);
        res.push(s.rtm);
    }
    Ok(res.into_iter().map(roundn).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stream_formatting_matches_c_plus_plus() {
        assert_eq!(ostream_double(0.123456789), "0.123457");
        assert_eq!(ostream_double(1.5), "1.5");
        assert_eq!(ostream_double(12.0), "12");
        assert_eq!(ostream_double(0.00001234), "1.234e-05");
        assert_eq!(ostream_double(1234567.0), "1.23457e+06");
        assert_eq!(ostream_double(-0.25), "-0.25");
    }

    #[test]
    fn atoi_stops_at_the_first_non_digit() {
        assert_eq!(atoi("123"), 123);
        assert_eq!(atoi("5e-05"), 5);
        assert_eq!(atoi("e5"), 0);
        assert_eq!(atoi("-7x"), -7);
    }

    #[test]
    fn leverages_round_to_significant_digits() {
        assert_eq!(round_to_digits(0.123456, 2), 0.12);
        assert_eq!(round_to_digits(0.0123456, 3), 0.0123);
        assert_eq!(round_to_digits(123456.0, 2), 120000.0);
    }

    #[test]
    fn clusters_of_nearly_equal_values() {
        // 0.12 and 0.13 are neighbours, 0.15 is not
        assert_eq!(cluster_sizes(vec![0.15, 0.12, 0.13]), vec![1.0, 2.0]);
        assert_eq!(cluster_sizes(vec![0.5]), vec![1.0]);
        assert!(cluster_sizes(vec![]).is_empty());
    }
}
