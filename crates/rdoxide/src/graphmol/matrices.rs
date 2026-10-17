//! Topological and 3D distance matrices, adjacency matrix and shortest paths
//! (`MolOps::getDistanceMat`, `getAdjacencyMatrix`, `get3DDistanceMat`,
//! `getShortestPath`).
//!
//! Ported from RDKit `Code/GraphMol/Matrices.cpp` (2026.09.1). With unit
//! bond weights, upstream's Floyd-Warshall gives the breadth-first-search
//! distances computed here (unreachable pairs stay at 1e8).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::VecDeque;

use crate::graphmol::romol::ROMol;
use crate::rdgeneral::dict::{RDProps, RDValue};
use crate::rdgeneral::threads::{par_row_blocks_mut, par_row_blocks_mut2};
use crate::{RDKitError, Result};

const LOCAL_INF: f64 = 1e8;

/// Rows of a distance matrix a worker thread must have to be worth starting.
const PAR_MIN_ROWS: usize = 64;

/// Rows per worker thread in a Floyd-Warshall pass: the threads are started
/// again for every pass, so they need enough rows to pay for that.
const FW_MIN_ROWS: usize = 128;

/// `getDistanceMat(mol)` (no bond orders, no atom weights): row-major
/// `n x n` matrix of bond-count distances.
pub fn get_distance_mat(mol: &ROMol) -> Vec<f64> {
    let n = mol.num_atoms() as usize;
    let mut d = vec![LOCAL_INF; n * n];
    // every row is a breadth-first search from its own atom
    par_row_blocks_mut(&mut d, n, PAR_MIN_ROWS, |first, block| {
        for (r, row) in block.chunks_mut(n).enumerate() {
            let i = first + r;
            row[i] = 0.0;
            let mut q = VecDeque::from([i as u32]);
            while let Some(a) = q.pop_front() {
                let da = row[a as usize];
                for nb in mol.atom_neighbors(a) {
                    if row[nb as usize] == LOCAL_INF {
                        row[nb as usize] = da + 1.0;
                        q.push_back(nb);
                    }
                }
            }
        }
    });
    d
}

/// Floyd-Warshall all-pairs shortest paths on the row-major `n x n` matrix
/// `d` (overwritten with the distances), keeping upstream's tie handling
/// and operation order so that summed bond weights come out bit-identical.
/// Returns the predecessor matrix (`-1` on the diagonal and for unreachable
/// pairs).
fn floyd_warshall(n: usize, d: &mut [f64]) -> Vec<i32> {
    // unit bond weights fit 16 bits, which cuts the memory traffic of the n
    // passes over the matrices to a third and doubles the SIMD width
    if n < usize::from(NO_PATH) && d.iter().all(|&w| w == 0.0 || w == 1.0 || w == LOCAL_INF) {
        return floyd_warshall_unit(n, d);
    }
    floyd_warshall_f64(n, d)
}

fn floyd_warshall_f64(n: usize, d: &mut [f64]) -> Vec<i32> {
    let mut cur_p: Vec<i32> = (0..n * n)
        .map(|k| {
            let (i, j) = (k / n, k % n);
            if i == j || d[k] == LOCAL_INF {
                -1
            } else {
                i as i32
            }
        })
        .collect();
    let mut cur_d = d.to_vec();
    let mut next_d = cur_d.clone();
    let mut next_p = cur_p.clone();
    for k in 0..n {
        // each pass reads the previous one only, so its rows are independent
        let (prev_d, prev_p) = (&cur_d, &cur_p);
        let (k_row_d, k_row_p) = (&prev_d[k * n..(k + 1) * n], &prev_p[k * n..(k + 1) * n]);
        par_row_blocks_mut2(&mut next_d, &mut next_p, n, FW_MIN_ROWS, |first, nd, np| {
            for (r, (row_d, row_p)) in nd.chunks_mut(n).zip(np.chunks_mut(n)).enumerate() {
                let i = first + r;
                let (old_d, old_p) = (&prev_d[i * n..(i + 1) * n], &prev_p[i * n..(i + 1) * n]);
                let dik = old_d[k];
                // branch-free so that the loop vectorizes; ties keep the old path
                for ((((out_d, out_p), &v1), (&p1, &p2)), &dkj) in row_d
                    .iter_mut()
                    .zip(row_p.iter_mut())
                    .zip(old_d)
                    .zip(old_p.iter().zip(k_row_p))
                    .zip(k_row_d)
                {
                    let v2 = dik + dkj;
                    let keep = v1 <= v2;
                    *out_d = if keep { v1 } else { v2 };
                    *out_p = if keep { p1 } else { p2 };
                }
            }
        });
        // the pass just written becomes the one the next pass reads
        std::mem::swap(&mut cur_d, &mut next_d);
        std::mem::swap(&mut cur_p, &mut next_p);
    }
    d.copy_from_slice(&cur_d);
    cur_p
}

/// "Unreachable" distance and "no predecessor" in [`floyd_warshall_unit`].
const NO_PATH: u16 = u16::MAX;

/// [`floyd_warshall`] for weights 0, 1 and `LOCAL_INF`: integer distances and
/// predecessors in `u16`, updated in place (row `k` does not change in pass
/// `k`, so a copy of it is all a pass has to read besides the row it writes).
fn floyd_warshall_unit(n: usize, d: &mut [f64]) -> Vec<i32> {
    let mut dist: Vec<u16> = d
        .iter()
        .map(|&w| if w == LOCAL_INF { NO_PATH } else { w as u16 })
        .collect();
    let mut pred: Vec<u16> = dist
        .iter()
        .enumerate()
        .map(|(k, &w)| {
            let (i, j) = (k / n, k % n);
            if i == j || w == NO_PATH {
                NO_PATH
            } else {
                i as u16
            }
        })
        .collect();
    let (mut k_dist, mut k_pred) = (vec![0u16; n], vec![0u16; n]);
    for k in 0..n {
        k_dist.copy_from_slice(&dist[k * n..(k + 1) * n]);
        k_pred.copy_from_slice(&pred[k * n..(k + 1) * n]);
        let (k_dist, k_pred) = (&k_dist[..], &k_pred[..]);
        par_row_blocks_mut2(&mut dist, &mut pred, n, FW_MIN_ROWS, |_, bd, bp| {
            for (row_d, row_p) in bd.chunks_mut(n).zip(bp.chunks_mut(n)) {
                let dik = row_d[k];
                if dik == NO_PATH {
                    continue;
                }
                // branch-free so that the loop vectorizes; ties keep the old path
                for (((d, p), &dkj), &pkj) in row_d
                    .iter_mut()
                    .zip(row_p.iter_mut())
                    .zip(k_dist)
                    .zip(k_pred)
                {
                    let via_k = dik.saturating_add(dkj);
                    let better = via_k < *d;
                    *d = if better { via_k } else { *d };
                    *p = if better { pkj } else { *p };
                }
            }
        });
    }
    for (out, &v) in d.iter_mut().zip(&dist) {
        *out = if v == NO_PATH {
            LOCAL_INF
        } else {
            f64::from(v)
        };
    }
    pred.iter()
        .map(|&p| if p == NO_PATH { -1 } else { i32::from(p) })
        .collect()
}

/// The bond-weight matrix `getDistanceMat` starts from: `1 / bond order`
/// (2/3 when aromatic) with `use_bo`, otherwise 1 per bond.
fn bond_weight_matrix(mol: &ROMol, use_bo: bool) -> Result<Vec<f64>> {
    let n = mol.num_atoms() as usize;
    let mut d = vec![LOCAL_INF; n * n];
    for i in 0..n {
        d[i * n + i] = 0.0;
    }
    for bond in mol.bonds() {
        let contrib = if !use_bo {
            1.0
        } else if bond.is_aromatic() {
            2.0 / 3.0
        } else {
            1.0 / bond.bond_type_as_double()?
        };
        let (i, j) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
        d[i * n + j] = contrib;
        d[j * n + i] = contrib;
    }
    Ok(d)
}

/// Sets the diagonal to `6 / atomic number` (`useAtomWts`).
fn apply_atom_weights(mol: &ROMol, d: &mut [f64]) {
    let n = mol.num_atoms() as usize;
    for i in 0..n {
        d[i * n + i] = 6.0 / f64::from(mol.atoms()[i].atomic_num());
    }
}

/// `getDistanceMat(mol, useBO, useAtomWts)`: with `use_bo` the bonds weigh
/// `1 / bond order` (2/3 when aromatic), with `use_atom_wts` the diagonal holds
/// `6 / atomic number`. The result is not cached on the molecule.
pub fn get_distance_mat_weighted(
    mol: &ROMol,
    use_bo: bool,
    use_atom_wts: bool,
) -> Result<Vec<f64>> {
    let n = mol.num_atoms() as usize;
    let mut d = if use_bo {
        let mut d = bond_weight_matrix(mol, true)?;
        floyd_warshall(n, &mut d);
        d
    } else {
        get_distance_mat(mol)
    };
    if use_atom_wts {
        apply_atom_weights(mol, &mut d);
    }
    Ok(d)
}

/// `getAdjacencyMatrix(mol, useBO, emptyVal)`. Upstream fills the matrix with
/// `memset(res, emptyVal, ...)`, so every byte of an empty entry is `emptyVal`.
pub fn get_adjacency_matrix(mol: &ROMol, use_bo: bool, empty_val: i32) -> Vec<f64> {
    let n = mol.num_atoms() as usize;
    let empty = f64::from_ne_bytes([empty_val as u8; 8]);
    let mut res = vec![empty; n * n];
    for bond in mol.bonds() {
        let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let (b, e) = (beg as usize, end as usize);
        if use_bo {
            res[b * n + e] = bond.valence_contrib(mol, beg);
            res[e * n + b] = bond.valence_contrib(mol, end);
        } else {
            res[b * n + e] = 1.0;
            res[e * n + b] = 1.0;
        }
    }
    res
}

/// Whether `name` is the key of a matrix cached by the `*_cached` functions.
/// Upstream stores these as `boost::shared_array` values, which the Python
/// property dictionaries cannot represent, so they are left out there.
pub fn is_matrix_cache_prop(name: &str) -> bool {
    name.contains("DistanceMatrix")
        || name.ends_with("AdjacencyMatrix")
        || name.ends_with("AdjacencyMatrixBO")
}

/// The cached matrix under `key` if it is there, usable and `force` is off.
fn cache_hit(mol: &ROMol, key: &str, force: bool) -> Option<Vec<f64>> {
    let n = mol.num_atoms() as usize;
    match mol.get_prop(key) {
        // upstream hands out a stale matrix of another size (reading past its
        // end); recomputing is the safe equivalent
        Some(RDValue::VecDouble(cached)) if !force && cached.len() == n * n => Some(cached.clone()),
        _ => None,
    }
}

/// `getDistanceMat(mol, useBO, useAtomWts, force, propNamePrefix)`. Caches
/// the distances and the predecessor matrix (`<key>_Paths`) as computed
/// properties.
pub fn get_distance_mat_cached(
    mol: &mut ROMol,
    use_bo: bool,
    use_atom_wts: bool,
    force: bool,
    prefix: &str,
) -> Result<Vec<f64>> {
    // the BO and non-BO matrices must not share a cache entry
    let key = format!(
        "{prefix}DistanceMatrix{}{}",
        if use_bo { "BO" } else { "" },
        if use_atom_wts { "AtomWts" } else { "" }
    );
    if let Some(d) = cache_hit(mol, &key, force) {
        return Ok(d);
    }
    let n = mol.num_atoms() as usize;
    let mut d = bond_weight_matrix(mol, use_bo)?;
    let paths = floyd_warshall(n, &mut d);
    if use_atom_wts {
        apply_atom_weights(mol, &mut d);
    }
    mol.set_prop(&key, d.clone(), true)?;
    mol.set_prop(&format!("{key}_Paths"), paths, true)?;
    Ok(d)
}

/// `getAdjacencyMatrix(mol, useBO, emptyVal, force, propNamePrefix)`.
pub fn get_adjacency_matrix_cached(
    mol: &mut ROMol,
    use_bo: bool,
    empty_val: i32,
    force: bool,
    prefix: &str,
) -> Result<Vec<f64>> {
    let key = format!("{prefix}AdjacencyMatrix{}", if use_bo { "BO" } else { "" });
    if let Some(a) = cache_hit(mol, &key, force) {
        return Ok(a);
    }
    let a = get_adjacency_matrix(mol, use_bo, empty_val);
    mol.set_prop(&key, a.clone(), true)?;
    Ok(a)
}

/// `get3DDistanceMat(mol, confId, useAtomWts, force, propNamePrefix)`; an
/// empty prefix disables the cache, as upstream.
pub fn get_3d_distance_mat_cached(
    mol: &mut ROMol,
    conf_id: i32,
    use_atom_wts: bool,
    force: bool,
    prefix: &str,
) -> Result<Vec<f64>> {
    let conf = mol.conformer(conf_id)?.id();
    if prefix.is_empty() {
        return get_3d_distance_mat_weighted(mol, conf_id, use_atom_wts);
    }
    let key = format!(
        "{prefix}3DDistanceMatrix_Conf{conf}{}",
        if use_atom_wts { "_AtomWeights" } else { "" }
    );
    if let Some(d) = cache_hit(mol, &key, force) {
        return Ok(d);
    }
    let d = get_3d_distance_mat_weighted(mol, conf_id, use_atom_wts)?;
    mol.set_prop(&key, d.clone(), true)?;
    Ok(d)
}

/// `MolOps::getShortestPath(mol, aid1, aid2)`: atoms on a shortest path
/// (both ends included), empty when the atoms are not connected.
pub fn get_shortest_path(mol: &ROMol, aid1: u32, aid2: u32) -> Result<Vec<u32>> {
    let nats = mol.num_atoms();
    if aid1 >= nats || aid2 >= nats {
        return Err(RDKitError::range("aid", aid1.max(aid2), nats));
    }
    if aid1 == aid2 {
        return Err(RDKitError::invariant("aid1 == aid2"));
    }
    let mut pred = vec![-1i64; nats as usize];
    pred[aid1 as usize] = -2;
    pred[aid2 as usize] = -3;
    let mut q = VecDeque::from([aid1]);
    let mut done = false;
    while !done {
        let Some(&cur) = q.front() else {
            break;
        };
        for nb in mol.atom_neighbors(cur) {
            match pred[nb as usize] {
                -1 => {
                    pred[nb as usize] = cur as i64;
                    q.push_back(nb);
                }
                -3 => {
                    pred[nb as usize] = cur as i64;
                    done = true;
                    break;
                }
                _ => {}
            }
        }
        q.pop_front();
    }
    let mut res = VecDeque::new();
    if done {
        let mut prev = aid2 as i64;
        res.push_back(aid2);
        loop {
            prev = pred[prev as usize];
            if prev != aid1 as i64 {
                res.push_front(prev as u32);
            } else {
                break;
            }
        }
        res.push_front(aid1);
    }
    Ok(res.into())
}

/// `get3DDistanceMat(mol, confId)` (no atom weights): row-major `n x n`
/// matrix of Euclidean distances in conformer `conf_id`.
pub fn get_3d_distance_mat(mol: &ROMol, conf_id: i32) -> Result<Vec<f64>> {
    get_3d_distance_mat_weighted(mol, conf_id, false)
}

/// `get3DDistanceMat(mol, confId, useAtomWts)`: with `use_atom_wts` the
/// diagonal holds `6 / atomic number`.
pub fn get_3d_distance_mat_weighted(
    mol: &ROMol,
    conf_id: i32,
    use_atom_wts: bool,
) -> Result<Vec<f64>> {
    let conf = mol.conformer(conf_id)?;
    let pos = conf.positions();
    let n = mol.num_atoms() as usize;
    let mut d = vec![0.0; n * n];
    if use_atom_wts {
        apply_atom_weights(mol, &mut d);
    }
    for i in 0..n {
        for j in i + 1..n {
            let dist = (pos[i] - pos[j]).length();
            d[i * n + j] = dist;
            d[j * n + i] = dist;
        }
    }
    Ok(d)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A graph with many equal-length shortest paths (rings, chords and a
    /// disconnected piece), where tie handling decides the predecessors.
    fn tie_prone_weights(n: usize, seed: u64) -> Vec<f64> {
        let mut d = vec![LOCAL_INF; n * n];
        let mut state = seed;
        let mut edge = |d: &mut Vec<f64>, a: usize, b: usize| {
            if a != b {
                d[a * n + b] = 1.0;
                d[b * n + a] = 1.0;
            }
        };
        for i in 0..n {
            d[i * n + i] = 0.0;
            if i + 1 < n - 3 {
                edge(&mut d, i, i + 1);
            }
        }
        for _ in 0..n / 2 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let a = (state >> 33) as usize % (n - 3);
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let b = (state >> 33) as usize % (n - 3);
            edge(&mut d, a, b);
        }
        d
    }

    #[test]
    fn unit_weight_floyd_warshall_matches_the_float_version() {
        for (n, seed) in [(1, 1), (5, 2), (12, 3), (40, 4), (150, 5), (300, 6)] {
            let w = tie_prone_weights(n.max(4), seed);
            let n = n.max(4);
            let (mut a, mut b) = (w.clone(), w);
            let pa = floyd_warshall_unit(n, &mut a);
            let pb = floyd_warshall_f64(n, &mut b);
            assert_eq!(a, b, "distances, n={n}");
            assert_eq!(pa, pb, "predecessors, n={n}");
        }
    }
}
