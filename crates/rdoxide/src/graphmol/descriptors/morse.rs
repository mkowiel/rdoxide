//! Molecule Representation of Structures based on Electron diffraction
//! (MoRSE) descriptors (upstream `GraphMol/Descriptors/MORSE.cpp`): for the 32
//! scattering parameters `s = 0..31` the sums `sum w_j w_k sin(s r_jk) /
//! (s r_jk)` for seven atom weightings (unweighted, mass, van der Waals
//! volume, electronegativity, polarizability, ionization potential and
//! intrinsic state), 224 values; a custom atom property gives 32.
//!
//! Atomic property values after Todeschini and Consonni, Molecular Descriptors
//! for Chemoinformatics, 2009.
//!
//! Copyright (c) 2016, Guillaume GODIN. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-Descriptors3D-Godin-BSD-3-Clause.txt).

use super::mol_data_3d as weights;
use super::radial::{distance_mat_3d, pair_profiles};
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

pub const MORSE_VERSION: &str = "1.0.0";

const NUM_PARAMETERS: usize = 32;

/// `MORSE(mol, res, confId, customAtomPropName)`.
pub fn morse(mol: &mut ROMol, conf_id: i32, custom_atom_prop_name: &str) -> Result<Vec<f64>> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    let dist = distance_mat_3d(mol, conf_id)?;
    let n = mol.num_atoms() as usize;
    let kernel = |i: usize, d: f64| {
        if i == 0 {
            1.0
        } else {
            let s = i as f64;
            (s * d).sin() / (s * d)
        }
    };
    if custom_atom_prop_name.is_empty() {
        let mass = weights::relative_mw(mol);
        let vdw = weights::relative_vdw(mol);
        let eneg = weights::relative_eneg(mol);
        let pol = weights::relative_pol(mol);
        let ion = weights::relative_ion_pol(mol);
        let state = weights::intrinsic_state(mol)?;
        let columns = [
            None,
            Some(&mass[..]),
            Some(&vdw[..]),
            Some(&eneg[..]),
            Some(&pol[..]),
            Some(&ion[..]),
            Some(&state[..]),
        ];
        Ok(pair_profiles(&dist, n, &columns, NUM_PARAMETERS, kernel))
    } else {
        let custom = weights::custom_atom_prop(mol, custom_atom_prop_name)?;
        Ok(pair_profiles(
            &dist,
            n,
            &[Some(&custom)],
            NUM_PARAMETERS,
            kernel,
        ))
    }
}
