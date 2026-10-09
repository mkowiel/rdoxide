//! Electronegativity equalization method (EEM) atomic partial charges
//! (upstream `GraphMol/Descriptors/EEM.cpp`), with the parameters of the NEEMP
//! article (J. Cheminf. 2016, 8, 17), `CCD_gen_DE_RMSD_B3LYP_6311G_NPA`,
//! trained for H, C, N, O, F, P, S, Cl and Br (not iodine).
//!
//! The charges solve the linear system of the method with a full-pivoting LU
//! decomposition: `B_i q_i + sum_j kappa / r_ij q_j - chi = -A_i` for every
//! atom and `sum q_i = total charge`. `A` and `B` depend on the element and on
//! the highest bond order of the atom in the Kekulé form of the molecule.
//!
//! Elements beyond bromine have no parameters (upstream reads past its
//! tables); they get zero parameters here.
//!
//! Copyright (c) 2016, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use crate::graphmol::kekulize::kekulize;
use crate::graphmol::matrices::get_3d_distance_mat_cached;
use crate::graphmol::romol::ROMol;
use crate::numerics::full_piv_lu::FullPivLu;
use crate::{RDKitError, Result};

pub const EEM_VERSION: &str = "1.0.0";

const KAPPA: f64 = 0.5125;

/// Parameter tables indexed by atomic number, one per highest bond order 1-3.
const A: [&[f64]; 3] = [&A1, &A2, &A3];
const B: [&[f64]; 3] = [&B1, &B2, &B3];

#[rustfmt::skip]
const A1: [f64; 36] = [
    0.0, 2.5473, 0.0, 0.0, 0.0, 0.0, 2.7221, 2.9750, 3.1503, 2.9976, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 2.6511, 2.7026, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6263,
];
#[rustfmt::skip]
const A2: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.7667, 2.8895, 3.0486, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 2.2933, 2.6471, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];
#[rustfmt::skip]
const A3: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6944, 3.0240, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];
#[rustfmt::skip]
const B1: [f64; 36] = [
    0.0, 1.1641, 0.0, 0.0, 0.0, 0.0, 0.6403, 0.9083, 1.0577, 0.9983, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.4897, 1.1537, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1105,
];
#[rustfmt::skip]
const B2: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6513, 0.6647, 0.8410, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.5759, 0.4512, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];
#[rustfmt::skip]
const B3: [f64; 36] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6776, 1.4240, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];

/// `getAtomtype`: the highest bond order of the atom (aromatic counts as 2).
fn atom_type(mol: &ROMol, idx: u32) -> Result<usize> {
    let mut t = 1;
    for bond in mol.atom_bonds(idx) {
        let mut order = bond.bond_type_as_double()?;
        if order == 1.5 {
            order = 2.0;
        }
        t = t.max(order as usize);
    }
    Ok(t)
}

/// The parameter of `table` for an atom of highest bond order `t`; 0 outside
/// the tabulated bond orders and elements.
fn parameter(table: [&[f64]; 3], t: usize, atomic_num: u32) -> f64 {
    t.checked_sub(1)
        .and_then(|i| table.get(i))
        .and_then(|row| row.get(atomic_num as usize))
        .copied()
        .unwrap_or(0.0)
}

/// `EEM(mol, res, confId)`: one charge per atom.
pub fn eem(mol: &ROMol, conf_id: i32) -> Result<Vec<f64>> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    let n = mol.num_atoms() as usize;
    // the Kekulé form is required for the bond orders; work on a copy
    let mut wmol = mol.clone();
    kekulize(&mut wmol, true, false, 100)?;
    let dist = get_3d_distance_mat_cached(&mut wmol, conf_id, false, true, "_")?;

    let types = (0..n as u32)
        .map(|i| atom_type(&wmol, i))
        .collect::<Result<Vec<_>>>()?;
    let nums: Vec<u32> = wmol.atoms().iter().map(|a| a.atomic_num()).collect();

    // The system is laid out like upstream's buffer: entry (i, j) lives at
    // `i * (n + 1) + j` and the buffer is read column-major.
    let dim = n + 1;
    let mut a = vec![0.0; dim * dim];
    for i in 0..n {
        a[i * dim + i] = parameter(B, types[i], nums[i]);
        for j in i + 1..n {
            let v = KAPPA / dist[i * n + j];
            a[i * dim + j] = v;
            a[j * dim + i] = v;
        }
    }
    for i in 0..n {
        a[n * dim + i] = 1.0;
        a[i * dim + n] = -1.0;
    }
    let mut b: Vec<f64> = (0..n).map(|j| -parameter(A, types[j], nums[j])).collect();
    b.push(f64::from(
        wmol.atoms().iter().map(|a| a.formal_charge()).sum::<i32>(),
    ));

    let mut solution = FullPivLu::new(dim, a).solve(&b);
    solution.truncate(n);
    Ok(solution)
}
