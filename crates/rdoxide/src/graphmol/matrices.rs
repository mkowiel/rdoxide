//! Topological distance matrix and shortest paths (`MolOps::getDistanceMat`,
//! `MolOps::getShortestPath`).
//!
//! Ported from RDKit `Code/GraphMol/Matrices.cpp` (2026.09.1). With unit
//! bond weights, upstream's Floyd-Warshall gives the breadth-first-search
//! distances computed here (unreachable pairs stay at 1e8).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::VecDeque;

use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

const LOCAL_INF: f64 = 1e8;

/// `getDistanceMat(mol)` (no bond orders, no atom weights): row-major
/// `n x n` matrix of bond-count distances.
pub fn get_distance_mat(mol: &ROMol) -> Vec<f64> {
    let n = mol.num_atoms() as usize;
    let mut d = vec![LOCAL_INF; n * n];
    for i in 0..n {
        d[i * n + i] = 0.0;
        let mut q = VecDeque::from([i as u32]);
        while let Some(a) = q.pop_front() {
            let da = d[i * n + a as usize];
            for nb in mol.atom_neighbors(a) {
                let k = i * n + nb as usize;
                if d[k] == LOCAL_INF {
                    d[k] = da + 1.0;
                    q.push_back(nb);
                }
            }
        }
    }
    d
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
    let conf = mol.conformer(conf_id)?;
    let pos = conf.positions();
    let n = mol.num_atoms() as usize;
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        for j in i + 1..n {
            let dist = (pos[i] - pos[j]).length();
            d[i * n + j] = dist;
            d[j * n + i] = dist;
        }
    }
    Ok(d)
}
