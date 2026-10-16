//! Wavy bonds for double bonds of unknown stereo (`addWavyBondsForStereoAny`
//! of `MolFileStereochem.cpp`, RDKit 2026.09.1).
//!
//! Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;

use crate::Result;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::romol::ROMol;

/// `StereoBondThresholds`: penalties for choosing a single bond as the wavy
/// one.
#[derive(Debug)]
pub struct StereoBondThresholds;

impl StereoBondThresholds {
    /// neighbouring double bond without stereo info
    pub const DBL_BOND_NO_STEREO: u32 = 1000;
    /// neighbouring double bond with stereo specified
    pub const DBL_BOND_SPECIFIED_STEREO: u32 = 10000;
    /// atom with specified chirality
    pub const CHIRAL_ATOM: u32 = 100000;
    /// single bond with the direction already set
    pub const DIRECTION_SET: u32 = 1000000;
}

/// The bonds sharing an atom with `bond`: those of the begin atom, then of
/// the end atom.
fn bond_neighbors(mol: &ROMol, bond: u32) -> Vec<u32> {
    let b = mol.bond(bond);
    [b.begin_atom_idx(), b.end_atom_idx()]
        .into_iter()
        .flat_map(|a| mol.atom_bonds(a).map(|nb| nb.idx()))
        .filter(|&i| i != bond)
        .collect()
}

/// The atom of `bond1` that `bond2` doesn't touch.
fn nonshared_atom(mol: &ROMol, bond1: u32, bond2: u32) -> Result<u32> {
    let (b1, b2) = (mol.bond(bond1), mol.bond(bond2));
    if b1.begin_atom_idx() == b2.begin_atom_idx() || b1.begin_atom_idx() == b2.end_atom_idx() {
        Ok(b1.end_atom_idx())
    } else if b1.end_atom_idx() == b2.begin_atom_idx() || b1.end_atom_idx() == b2.end_atom_idx() {
        Ok(b1.begin_atom_idx())
    } else {
        Err(crate::RDKitError::invariant("bonds don't share an atom"))
    }
}

/// Sets wavy bonds around double bonds with `STEREOANY` stereo.
///
/// The penalties are all a tenth of the corresponding threshold, so the
/// penalties of one single bond can be added up and still compared with a
/// single set of thresholds.
///
/// * `clear_double_bond_flags`: also remove the flags for unknown double bond
///   stereo.
/// * `add_when_impossible`: if non-zero a neighbouring single bond is made
///   wavy even if it connects to a chiral centre or double bond with
///   specified stereo (as long as its penalty is not above this value).
pub fn add_wavy_bonds_for_stereo_any(
    mol: &mut ROMol,
    clear_double_bond_flags: bool,
    add_when_impossible: u32,
) -> Result<()> {
    let n_bonds = mol.num_bonds() as usize;
    let mut single_bond_scores = vec![0i64; n_bonds];
    let mut single_bond_neighbors: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    let mut double_bonds_to_set = vec![false; n_bonds];
    for dbl in 0..n_bonds as u32 {
        if mol.bond(dbl).bond_type() != BondType::Double {
            continue;
        }
        let dbl_stereo = mol.bond(dbl).stereo();
        if dbl_stereo == BondStereo::StereoAny {
            double_bonds_to_set[dbl as usize] = true;
        }
        for single in bond_neighbors(mol, dbl) {
            if mol.bond(single).bond_type() != BondType::Single {
                continue;
            }
            let mut score = single_bond_scores[single as usize] + 1;
            if mol.bond(single).bond_dir() != BondDir::None {
                score += i64::from(StereoBondThresholds::DIRECTION_SET / 10);
            }
            if dbl_stereo == BondStereo::StereoNone {
                score += i64::from(StereoBondThresholds::DBL_BOND_NO_STEREO / 10);
            } else if dbl_stereo > BondStereo::StereoAny {
                score += i64::from(StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO / 10);
            }
            let other = nonshared_atom(mol, single, dbl)?;
            // favour atoms with fewer neighbours
            score += 10 * i64::from(mol.atom_degree(other));
            let tag = mol.atom(other).chiral_tag();
            if tag != ChiralType::Unspecified && tag != ChiralType::Other {
                score += i64::from(StereoBondThresholds::CHIRAL_ATOM / 10);
            }
            single_bond_scores[single as usize] = score;
            if dbl_stereo == BondStereo::StereoAny {
                single_bond_neighbors.entry(single).or_default().push(dbl);
            }
        }
    }
    // (-number of double bonds served, score, bond index)
    let mut sorted: Vec<(i64, i64, usize)> = (0..n_bonds)
        .filter(|&i| single_bond_scores[i] != 0)
        .map(|i| {
            let n = single_bond_neighbors.get(&(i as u32)).map_or(0, Vec::len);
            (-(n as i64), single_bond_scores[i], i)
        })
        .collect();
    sorted.sort();
    for &(_, score, bond) in &sorted {
        let dbls = single_bond_neighbors
            .get(&(bond as u32))
            .cloned()
            .unwrap_or_default();
        for dbl in dbls {
            if !double_bonds_to_set[dbl as usize] {
                continue;
            }
            if add_when_impossible != 0 {
                if score > i64::from(add_when_impossible) {
                    continue;
                }
            } else if score > i64::from(StereoBondThresholds::DBL_BOND_NO_STEREO) {
                crate::rd_log_warning!(
                    "Setting wavy bond flag on bond {bond} which may make other stereo info ambiguous\n"
                );
            }
            mol.bond_mut(bond as u32).set_bond_dir(BondDir::Unknown);
            if clear_double_bond_flags {
                let d = mol.bond_mut(dbl);
                if d.bond_dir() == BondDir::EitherDouble {
                    d.set_bond_dir(BondDir::None);
                }
                d.set_stereo(BondStereo::StereoNone)?;
            }
            double_bonds_to_set[dbl as usize] = false;
        }
    }
    if add_when_impossible != 0 && double_bonds_to_set.iter().any(|&b| b) {
        let ids: String = double_bonds_to_set
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| format!(" {i}"))
            .collect();
        crate::rd_log_warning!(" unable to set wavy bonds for double bonds:{ids}\n");
    }
    Ok(())
}
