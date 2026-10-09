//! Per-atom weights shared by the autocorrelation-style descriptors (upstream
//! `GraphMol/Descriptors/MolData3Ddescriptors.cpp`).
//!
//! Dummy atoms (atomic number 0) are looked up at index -1 upstream, which
//! reads the last entry of the table declared before; that is reproduced here
//! (the mass table reads 0). Atomic numbers above 110 get a weight of 0.
//!
//! Copyright (c) 2016-2017, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use super::data_3d::{self, NUM_ELEMENTS};
use crate::RDKitError;
use crate::Result;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::dict::{RDProps, RDValue};

/// Weight of every atom looked up by atomic number; `dummy` is the weight of
/// atomic number 0.
fn weights(mol: &ROMol, table: &[f64; NUM_ELEMENTS], dummy: f64) -> Vec<f64> {
    mol.atoms()
        .iter()
        .map(|atom| match atom.atomic_num() as usize {
            0 => dummy,
            z @ 1..=NUM_ELEMENTS => table[z - 1],
            _ => 0.0,
        })
        .collect()
}

fn last(table: &[f64; NUM_ELEMENTS]) -> f64 {
    table[NUM_ELEMENTS - 1]
}

/// `GetRelativeMW`
pub fn relative_mw(mol: &ROMol) -> Vec<f64> {
    weights(mol, &data_3d::MW, 0.0)
}

/// `GetRelativePol`
pub fn relative_pol(mol: &ROMol) -> Vec<f64> {
    weights(mol, &data_3d::POL, last(&data_3d::NEG))
}

/// `GetRelativeVdW`
pub fn relative_vdw(mol: &ROMol) -> Vec<f64> {
    weights(mol, &data_3d::VDW, last(&data_3d::MW))
}

/// `GetRelativeENeg`
pub fn relative_eneg(mol: &ROMol) -> Vec<f64> {
    weights(mol, &data_3d::NEG, last(&data_3d::VDW))
}

/// `GetRelativeIonPol`
pub fn relative_ion_pol(mol: &ROMol) -> Vec<f64> {
    weights(mol, &data_3d::ION_POL, last(&data_3d::POL))
}

/// `GetUn`: a weight of one for every atom.
pub fn unit_weights(mol: &ROMol) -> Vec<f64> {
    vec![1.0; mol.num_atoms() as usize]
}

/// `GetRelativeRcov`: covalent radius relative to carbon.
pub fn relative_rcov(mol: &ROMol) -> Vec<f64> {
    let carbon = data_3d::RCOV[5];
    weights(mol, &data_3d::RCOV, last(&data_3d::ION_POL))
        .into_iter()
        .map(|r| r / carbon)
        .collect()
}

/// `GetCustomAtomProp`: the double property of every atom, 1 where it is
/// missing.
pub fn custom_atom_prop(mol: &ROMol, name: &str) -> Result<Vec<f64>> {
    mol.atoms()
        .iter()
        .map(|atom| match atom.get_prop(name) {
            None => Ok(1.0),
            // `getPropIfPresent<double>` casts doubles and floats, and parses strings
            Some(RDValue::String(s)) => s
                .trim_end()
                .parse::<f64>()
                .map_err(|_| RDKitError::Runtime("bad any cast".into())),
            Some(v) => v
                .as_double()
                .ok_or_else(|| RDKitError::Runtime("bad any cast".into())),
        })
        .collect()
}

/// `GetPrincipalQuantumNumber`
fn principal_quantum_number(atomic_num: u32) -> u32 {
    match atomic_num {
        0..=2 => 1,
        3..=10 => 2,
        11..=18 => 3,
        19..=36 => 4,
        37..=54 => 5,
        55..=86 => 6,
        _ => 7,
    }
}

/// `GetIState`: the intrinsic state of every atom (1 for hydrogens and
/// isolated atoms), rounded to three decimals.
pub fn intrinsic_state(mol: &ROMol) -> Result<Vec<f64>> {
    let pt = PeriodicTable::get();
    let mut states = vec![1.0; mol.num_atoms() as usize];
    for (i, state) in states.iter_mut().enumerate() {
        let idx = i as u32;
        let atomic_num = mol.atom(idx).atomic_num();
        let degree = mol.atom_degree(idx) as i32;
        if degree > 0 && atomic_num > 1 {
            let h = mol.total_num_hs(idx, true)? as i32;
            let dv = pt.element(atomic_num).n_outer_elecs - h;
            let n = principal_quantum_number(atomic_num) as i32;
            let d = f64::from(degree - h);
            if d > 0.0 {
                *state =
                    (1000.0 * (4.0 / f64::from(n * n) * f64::from(dv) + 1.0) / d).round() / 1000.0;
            }
        }
    }
    Ok(states)
}
