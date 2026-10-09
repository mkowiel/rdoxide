//! Weighted Holistic Invariant Molecular (WHIM) descriptors (upstream
//! `GraphMol/Descriptors/WHIM.cpp`), after Todeschini and Gramatica and the
//! Dragon 6 selection.
//!
//! For an atom weighting the weighted covariance matrix of the centered
//! coordinates is diagonalized (by SVD); its eigenvalues give the directional
//! WHIMs `L1-L3`, their combinations the global shape WHIMs (`T`, `A`, `V`,
//! `P1-P3`, `K`), the kurtosis along the axes `E1-E3` and `D`, and the
//! symmetry of the atom positions along the axes `G1-G3`, `G`. Seven weightings
//! (unit, mass, van der Waals volume, electronegativity, polarizability,
//! ionization potential, intrinsic state) give 114 reported values; a custom
//! atom property gives 17. Every value is rounded to three decimals.
//!
//! Copyright (c) 2016-2017, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use super::mol_data_3d as weights;
use crate::geometry::Point3D;
use crate::graphmol::romol::ROMol;
use crate::numerics::jacobi_svd::jacobi_svd_3x3;
use crate::rdgeneral::fp::fmuladd;
use crate::{RDKitError, Result};

pub const WHIM_VERSION: &str = "1.0.0";

/// Number of intermediate values per weighting, in the order
/// L1 L2 L3 T A V P1 P2 P3 K E1 E2 E3 D G1 G2 G3 G.
const NUM_INTERMEDIATE: usize = 18;
/// Indices of the Dragon selection L1 L2 L3 P1 P2 G1 G2 G3 E1 E2 E3.
const DRAGON_SELECTION: [usize; 11] = [0, 1, 2, 6, 7, 14, 15, 16, 10, 11, 12];
/// Selection of the single custom-property weighting:
/// L1 L2 L3 P1 P2 G1 G2 G3 E1 E2 E3 T A G K D V.
const CUSTOM_SELECTION: [usize; 17] = [0, 1, 2, 6, 7, 14, 15, 16, 10, 11, 12, 3, 4, 17, 9, 13, 5];

fn roundn(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

/// Eigen evaluates products of matrices with `rows + cols + depth` below this
/// bound coefficient by coefficient (`EIGEN_GEMM_TO_COEFFBASED_THRESHOLD`).
const SMALL_PRODUCT: usize = 20;

/// Doubles per SIMD register of the reference build (NEON).
const PACKET: usize = 2;

/// Dot product of two vectors the way Eigen sums the products: `fused` for
/// the multiply-add chains of its vectorized kernels, otherwise a plain sum of
/// separately rounded products.
fn dot(a: impl Iterator<Item = f64>, b: impl Iterator<Item = f64>, fused: bool) -> f64 {
    let mut pairs = a.zip(b);
    let Some((x, y)) = pairs.next() else {
        return 0.0;
    };
    pairs.fold(x * y, |acc, (x, y)| {
        if fused {
            fmuladd(x, y, acc)
        } else {
            acc + x * y
        }
    })
}

/// Whether row `row` of a column of `rows` entries starting `offset` elements
/// after a register boundary is computed by Eigen's vectorized path: the rows
/// between the first aligned element and the last full register.
fn in_packet(row: usize, rows: usize, offset: usize) -> bool {
    let start = offset % PACKET;
    let end = start + (rows - start.min(rows)) / PACKET * PACKET;
    row >= start && row < end
}

/// `sum()` of a contiguous vector whose first element is `offset` elements
/// past a register boundary (Eigen's vectorized reduction: two interleaved
/// register accumulators, the lanes added at the end, then the unaligned head
/// and tail one by one).
fn packet_sum(v: &[f64], offset: usize) -> f64 {
    let start = (offset % PACKET).min(v.len());
    let size = v.len() - start;
    let aligned = size / PACKET * PACKET;
    if aligned == 0 {
        return v.iter().skip(1).fold(v[0], |acc, x| acc + x);
    }
    let end = start + aligned;
    let end2 = start + size / (2 * PACKET) * (2 * PACKET);
    let load = |i: usize| [v[i], v[i + 1]];
    let mut r0 = load(start);
    if aligned > PACKET {
        let mut r1 = load(start + PACKET);
        let mut i = start + 2 * PACKET;
        while i < end2 {
            let (a, b) = (load(i), load(i + PACKET));
            r0 = [r0[0] + a[0], r0[1] + a[1]];
            r1 = [r1[0] + b[0], r1[1] + b[1]];
            i += 2 * PACKET;
        }
        r0 = [r0[0] + r1[0], r0[1] + r1[1]];
        if end > end2 {
            let a = load(end2);
            r0 = [r0[0] + a[0], r0[1] + a[1]];
        }
    }
    let mut res = r0[0] + r0[1];
    for x in &v[..start] {
        res += x;
    }
    for x in &v[end..] {
        res += x;
    }
    res
}

/// Centered coordinates, as `x[atom][axis]`.
fn centered(points: &[Point3D]) -> Vec<[f64; 3]> {
    let n = points.len();
    let mean: Vec<f64> = (0..3)
        .map(|axis| {
            let column: Vec<f64> = points.iter().map(|p| [p.x, p.y, p.z][axis]).collect();
            // column `axis` of the n x 3 coordinate matrix starts axis * n elements in
            packet_sum(&column, axis * n) / n as f64
        })
        .collect();
    points
        .iter()
        .map(|p| [p.x - mean[0], p.y - mean[1], p.z - mean[2]])
        .collect()
}

/// `x^4` the way Eigen's array `pow(4)` computes it.
fn pow4(x: f64) -> f64 {
    let x2 = x * x;
    x2 * x2
}

/// The 18 WHIM values of one atom weighting (`getWhimD`).
fn whim_values(weight: &[f64], points: &[Point3D], th: f64) -> [f64; NUM_INTERMEDIATE] {
    let n = points.len();
    let mut total_weight = weight.iter().fold(0.0, |acc, w| acc + w);
    // the sum is close to zero for charge weightings
    if total_weight.abs() < 1e-4 {
        total_weight = 1.0;
    }
    let x = centered(points);

    // weighted covariance matrix of the coordinates
    let cov_small = 3 + n + 3 < SMALL_PRODUCT;
    let mut cov = [[0.0; 3]; 3];
    for (i, row) in cov.iter_mut().enumerate() {
        for (l, c) in row.iter_mut().enumerate() {
            let weighted = (0..n).map(|k| x[k][i] * weight[k]);
            // the 3 x 3 result is stored column-major: column l starts 3 * l elements in
            let fused = !cov_small || in_packet(i, 3, 3 * l);
            *c = dot(weighted, (0..n).map(|k| x[k][l]), fused) / total_weight;
        }
    }
    let svd = jacobi_svd_3x3(&cov);
    let sv = svd.singular_values;
    let scores_small = n + 3 + 3 < SMALL_PRODUCT;
    let mut scores = vec![[0.0; 3]; n];
    for (atom, xi) in x.iter().enumerate() {
        for c in 0..3 {
            // the n x 3 result is stored column-major: column c starts c * n elements in
            let fused = !scores_small || in_packet(atom, n, c * n);
            scores[atom][c] = dot(xi.iter().copied(), (0..3).map(|k| svd.v[k][c]), fused);
        }
    }

    let mut w = [0.0; NUM_INTERMEDIATE];
    w[0] = sv[0];
    w[1] = sv[1];
    w[2] = sv[2];
    w[3] = sv[0] + sv[1] + sv[2]; // T
    w[4] = fmuladd(sv[1], sv[2], fmuladd(sv[0], sv[1], sv[0] * sv[2])); // A
    w[5] = fmuladd(sv[0] * sv[1], sv[2], w[3] + w[4]); // V
    let sum = sv[0] + sv[1] + sv[2];
    w[6] = sv[0] / sum; // P1
    w[7] = sv[1] / sum; // P2
    w[8] = sv[2] / sum; // P3
    let res: f64 = (0..3).fold(0.0, |acc, i| acc + (w[i] / w[3] - 1.0 / 3.0).abs());
    w[9] = 3.0 / 4.0 * res; // K

    // inverse of the kurtosis along the principal axes
    for axis in 0..3 {
        let fourth_powers: Vec<f64> = scores.iter().map(|s| pow4(s[axis])).collect();
        let fourth = packet_sum(&fourth_powers, 0);
        w[10 + axis] = if fourth > 0.0 {
            n as f64 * (w[axis] * w[axis]) / fourth
        } else {
            0.0
        };
    }
    w[13] = (w[10] + w[11] + w[12]) / 3.0; // D

    // symmetry of the positions along the principal axes
    for s in &mut scores {
        for v in s.iter_mut() {
            *v = roundn(*v);
        }
    }
    let n_atoms = n as f64;
    let mut gamma = [0.0; 3];
    for (axis, g) in gamma.iter_mut().enumerate() {
        let mut symmetric = vec![0.0; 2 * n];
        let (mut ns, mut na) = (0.0_f64, 0.0_f64);
        for j in 0..n {
            let mut matched = false;
            for k in 0..n {
                if j == k {
                    continue;
                }
                // opposite positions that are not on the axis
                if (scores[j][axis] + scores[k][axis]).abs() <= th {
                    ns += 1.0;
                    matched = true;
                    symmetric[j] = 1.0;
                    symmetric[j + n] = 2.0;
                    symmetric[k] = 1.0;
                    symmetric[k + n] = 2.0;
                    break;
                }
            }
            if !matched {
                na += 1.0;
                symmetric[j] = 0.0;
                symmetric[j + n] = scores[j][axis].abs();
            }
        }
        // atoms on the axis count as symmetric
        for aj in 0..n {
            if symmetric[aj + n] < th && symmetric[aj] < 1.0 {
                ns += 1.0;
                na -= 1.0;
            }
        }
        let mut gamma_inv = 1.0;
        if ns == 0.0 {
            gamma_inv = 1.0 - (na / n_atoms) * (1.0 / n_atoms).ln() / 2.0_f64.ln();
        }
        if ns > 0.0 {
            gamma_inv = 1.0
                - ((ns / n_atoms) * (ns / n_atoms).ln() / 2.0_f64.ln()
                    + (na / n_atoms) * (1.0 / n_atoms).ln() / 2.0_f64.ln());
        }
        *g = 1.0 / gamma_inv;
    }
    w[14] = gamma[0];
    w[15] = gamma[1];
    w[16] = gamma[2];
    w[17] = (gamma[0] * gamma[1] * gamma[2]).powf(1.0 / 3.0);
    w
}

/// `WHIM(mol, res, confId, th, customAtomPropName)`: 114 values, or 17 for a
/// custom atom property.
pub fn whim(mol: &ROMol, conf_id: i32, th: f64, custom_atom_prop_name: &str) -> Result<Vec<f64>> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    let points = mol.conformer(conf_id)?.positions().to_vec();
    if !custom_atom_prop_name.is_empty() {
        let weight = weights::custom_atom_prop(mol, custom_atom_prop_name)?;
        let w = whim_values(&weight, &points, th);
        return Ok(CUSTOM_SELECTION.iter().map(|&i| roundn(w[i])).collect());
    }

    let weightings = [
        weights::unit_weights(mol),
        weights::relative_mw(mol),
        weights::relative_vdw(mol),
        weights::relative_eneg(mol),
        weights::relative_pol(mol),
        weights::relative_ion_pol(mol),
        weights::intrinsic_state(mol)?,
    ];
    let w: Vec<[f64; NUM_INTERMEDIATE]> = weightings
        .iter()
        .map(|weight| whim_values(weight, &points, th))
        .collect();

    // Dragon order: 7 x (L1 L2 L3 P1 P2 G1 G2 G3 E1 E2 E3), then T, A, G
    // (unit and mass weighting only), K, D and V of every weighting
    let mut res = vec![0.0; 114];
    for (k, wk) in w.iter().enumerate() {
        for (i, &sel) in DRAGON_SELECTION.iter().enumerate() {
            res[i + 11 * k] = roundn(wk[sel]);
        }
        res[k + 11 * 7] = roundn(wk[3]);
        res[k + 12 * 7] = roundn(wk[4]);
        res[k + 13 * 7 + 2] = roundn(wk[9]);
        res[k + 14 * 7 + 2] = roundn(wk[13]);
        res[k + 15 * 7 + 2] = roundn(wk[5]);
    }
    for (i, wi) in w.iter().take(2).enumerate() {
        res[i + 13 * 7] = roundn(wi[17]);
    }
    Ok(res)
}
