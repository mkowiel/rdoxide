//! Coulomb matrix (upstream `GraphMol/Descriptors/CoulombMat.cpp`): the
//! diagonal holds `0.5 Z^2.4`, the off-diagonal `Z_i Z_j / r_ij`.
//!
//! Copyright (c) 2018, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

pub const COULOMB_MAT_VERSION: &str = "1.0.0";

/// `CoulombMat(mol, res, confId)`: the symmetric `n x n` matrix, row by row.
pub fn coulomb_mat(mol: &ROMol, conf_id: i32) -> Result<Vec<Vec<f64>>> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    let n = mol.num_atoms() as usize;
    let pos = mol.conformer(conf_id)?.positions();
    let z: Vec<f64> = mol
        .atoms()
        .iter()
        .map(|a| f64::from(a.atomic_num()))
        .collect();
    let mut res = vec![vec![0.0; n]; n];
    for i in 0..n {
        res[i][i] = 0.5 * z[i].powf(2.4);
        for j in 0..i {
            let v = z[i] * z[j] / (pos[i] - pos[j]).length();
            res[i][j] = v;
            res[j][i] = v;
        }
    }
    Ok(res)
}
