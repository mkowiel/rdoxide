//! The new (non-legacy) stereo perception (upstream `Chirality.cpp`:
//! `stereoPerception`, `updateDoubleBondStereo`, `assignBondCisTrans` and
//! the post-sanitization atropisomer detection).

use super::find_stereo::{
    NOATOM, StereoDescriptor, StereoInfo, StereoSpecified, StereoType, find_potential_stereo,
};
use crate::graphmol::atropisomers::{
    cleanup_atropisomer_stereo_groups, detect_atropisomer_chirality,
};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::mol_ops::cleanup_stereo_groups;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result, rd_log_warning};

/// Detects atropisomers and discards the ones that cannot actually rotate:
/// the ring check of `cleanupAtropisomers`, which has already had its turn
/// before sanitization finished, is redone here.
fn detect_atropisomers_post_sanitization(mol: &mut ROMol, clean_it: bool) -> Result<()> {
    let conf = (mol.num_conformers() > 0).then_some(0);
    detect_atropisomer_chirality(mol, conf, clean_it)?;
    if !clean_it || !mol.ring_info().is_sssr_or_better() {
        return Ok(());
    }
    let mut removed_any = false;
    for idx in 0..mol.num_bonds() {
        // bonds in macrocycles (rings with 9 or more members) are left
        // alone, since they can link actual atropisomeric portions
        let stereo = mol.bond(idx).stereo();
        if matches!(
            stereo,
            BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW
        ) && mol.ring_info().num_bond_rings(idx)? > 0
            && mol.ring_info().min_bond_ring_size(idx)? < 9
        {
            mol.bond_mut(idx).set_stereo(BondStereo::StereoNone)?;
            removed_any = true;
        }
    }
    if removed_any {
        cleanup_atropisomer_stereo_groups(mol)?;
    }
    Ok(())
}

fn is_wedge_dir(dir: BondDir) -> bool {
    matches!(dir, BondDir::EndDownRight | BondDir::EndUpRight)
}

fn flipped(dir: BondDir) -> BondDir {
    if dir == BondDir::EndDownRight {
        BondDir::EndUpRight
    } else {
        BondDir::EndDownRight
    }
}

/// `Chirality::assignBondCisTrans`: cis/trans from the single-bond
/// directions, for a double bond the stereo search found unspecified.
fn assign_bond_cis_trans(mol: &mut ROMol, sinfo: &StereoInfo) -> Result<()> {
    let c = &sinfo.controlling_atoms;
    if sinfo.stereo_type != StereoType::BondDouble
        || sinfo.specified != StereoSpecified::Unspecified
        || c.len() != 4
        || (c[0] == NOATOM && c[1] == NOATOM)
        || (c[2] == NOATOM && c[3] == NOATOM)
    {
        return Ok(());
    }
    let dbl = mol.bond(sinfo.centered_on);
    let (dbl_begin, dbl_end) = (dbl.begin_atom_idx(), dbl.end_atom_idx());

    let mut beg_first_neighbor = true;
    let mut beg_bond = mol
        .bond_idx_between(dbl_begin, c[0])
        .ok_or_else(|| RDKitError::invariant("no initial bond found"))?;
    let mut beg_dir = mol.bond(beg_bond).bond_dir();
    if !is_wedge_dir(beg_dir) {
        beg_first_neighbor = false;
        if c[1] != NOATOM {
            beg_bond = mol
                .bond_idx_between(dbl_begin, c[1])
                .ok_or_else(|| RDKitError::invariant("no initial bond found"))?;
            beg_dir = mol.bond(beg_bond).bond_dir();
        }
    }
    if !is_wedge_dir(beg_dir) {
        return Ok(());
    }
    if mol.bond(beg_bond).begin_atom_idx() != dbl_begin {
        beg_dir = flipped(beg_dir);
    }

    let mut end_first_neighbor = true;
    let mut end_bond = mol
        .bond_idx_between(dbl_end, c[2])
        .ok_or_else(|| RDKitError::invariant("no final bond found"))?;
    let mut end_dir = mol.bond(end_bond).bond_dir();
    if !is_wedge_dir(end_dir) {
        end_first_neighbor = false;
        if c[3] != NOATOM {
            end_bond = mol
                .bond_idx_between(dbl_end, c[3])
                .ok_or_else(|| RDKitError::invariant("no final bond found"))?;
            end_dir = mol.bond(end_bond).bond_dir();
        }
    }
    if !is_wedge_dir(end_dir) {
        return Ok(());
    }
    if mol.bond(end_bond).begin_atom_idx() != dbl_end {
        end_dir = flipped(end_dir);
    }

    // both single bonds are viewed as starting at the double bond: when
    // they point the same way the double bond is cis
    let mut same_dir = beg_dir == end_dir;
    // unless exactly one side used its second neighbor
    if beg_first_neighbor ^ end_first_neighbor {
        same_dir = !same_dir;
    }
    mol.set_bond_stereo_atoms(sinfo.centered_on, c[0], c[2])?;
    mol.bond_mut(sinfo.centered_on).set_stereo(if same_dir {
        BondStereo::StereoCis
    } else {
        BondStereo::StereoTrans
    })?;
    Ok(())
}

/// `updateDoubleBondStereo`: stores the double bond stereo that was found.
fn update_double_bond_stereo(mol: &mut ROMol, sinfo: &[StereoInfo], clean_it: bool) -> Result<()> {
    let mut touched = vec![false; mol.num_bonds() as usize];
    for si in sinfo {
        if si.stereo_type != StereoType::BondDouble {
            continue;
        }
        let idx = si.centered_on;
        touched[idx as usize] = true;
        mol.bond_mut(idx).set_stereo(BondStereo::StereoNone)?;
        match si.specified {
            StereoSpecified::Specified => {
                if si.controlling_atoms.len() != 4 {
                    return Err(RDKitError::invariant("expected four controlling atoms"));
                }
                mol.set_bond_stereo_atoms(idx, si.controlling_atoms[0], si.controlling_atoms[2])?;
                match si.descriptor {
                    StereoDescriptor::BondCis => {
                        mol.bond_mut(idx).set_stereo(BondStereo::StereoCis)?;
                    }
                    StereoDescriptor::BondTrans => {
                        mol.bond_mut(idx).set_stereo(BondStereo::StereoTrans)?;
                    }
                    _ => rd_log_warning!("unrecognized bond stereo type"),
                }
            }
            StereoSpecified::Unknown => {
                // imines without explicit Hs have an unknown double bond
                // with a single neighbor on the N: clear those
                let c = &si.controlling_atoms;
                if c.len() == 4 && c[0] != NOATOM && c[2] != NOATOM {
                    mol.set_bond_stereo_atoms(idx, c[0], c[2])?;
                    mol.bond_mut(idx).set_stereo(BondStereo::StereoAny)?;
                } else {
                    mol.bond_mut(idx).set_stereo(BondStereo::StereoNone)?;
                }
                mol.bond_mut(idx).set_bond_dir(BondDir::None);
            }
            StereoSpecified::Unspecified => assign_bond_cis_trans(mol, si)?,
        }
    }
    if clean_it {
        for idx in 0..mol.num_bonds() {
            if touched[idx as usize] || mol.bond(idx).bond_type() != BondType::Double {
                continue;
            }
            // not seen above, so it cannot have stereo
            let bond = mol.bond_mut(idx);
            bond.set_stereo(BondStereo::StereoNone)?;
            bond.set_bond_dir(BondDir::None);
            bond.stereo_atoms_mut().clear();
        }
    }
    Ok(())
}

/// `Chirality::stereoPerception`: the new stereo perception, based on
/// `findPotentialStereo`.
pub fn stereo_perception(
    mol: &mut ROMol,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
) -> Result<()> {
    if clean_it {
        for idx in 0..mol.num_atoms() {
            let atom = mol.atom_mut(idx);
            atom.clear_prop(cp::_CIPCode);
            atom.clear_prop(cp::_ChiralityPossible);
        }
        for idx in 0..mol.num_bonds() {
            let bond = mol.bond_mut(idx);
            bond.clear_prop(cp::_CIPCode);
            if bond.bond_dir() == BondDir::EitherDouble {
                bond.set_stereo(BondStereo::StereoAny)?;
                bond.stereo_atoms_mut().clear();
                bond.set_bond_dir(BondDir::None);
            }
        }
    }
    // the double bonds need cis/trans markers
    super::bond_dirs::set_bond_stereo_from_directions(mol)?;

    let sinfo = find_potential_stereo(mol, clean_it, flag_possible_stereo_centers)?;
    if flag_possible_stereo_centers {
        for si in &sinfo {
            if matches!(
                si.stereo_type,
                StereoType::AtomTetrahedral
                    | StereoType::AtomSquarePlanar
                    | StereoType::AtomTrigonalBipyramidal
                    | StereoType::AtomOctahedral
            ) {
                mol.atom_mut(si.centered_on)
                    .dict_mut()
                    .set_val(cp::_ChiralityPossible, 1i32);
            }
        }
    }
    update_double_bond_stereo(mol, &sinfo, clean_it)?;
    detect_atropisomers_post_sanitization(mol, clean_it)?;
    if clean_it {
        cleanup_atropisomer_stereo_groups(mol)?;
        cleanup_stereo_groups(mol)?;
    }
    Ok(())
}
