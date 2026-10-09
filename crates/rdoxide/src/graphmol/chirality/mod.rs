//! Port of RDKit `Code/GraphMol/Chirality.{h,cpp}` (in progress) and
//! `NontetrahedralStereo.cpp`.

pub mod bond_dirs;
pub mod from_coords;
pub mod legacy;
pub mod nontetrahedral;
pub mod stereo_groups;

use std::sync::atomic::{AtomicU8, Ordering};

/// Double bond stereo is ignored/removed for rings smaller than this.
pub const MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO: u32 = 8;

const NON_TETRAHEDRAL_STEREO_ENV_VAR: &str = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";
const USE_LEGACY_STEREO_ENV_VAR: &str = "RDK_USE_LEGACY_STEREO_PERCEPTION";

// 0: not set (consult the environment), 1: false, 2: true
static ALLOW_NONTETRAHEDRAL: AtomicU8 = AtomicU8::new(0);
static USE_LEGACY: AtomicU8 = AtomicU8::new(0);

/// `getValFromEnvironment`: any value other than "0" means true.
fn val_from_environment(var: &str, default: bool) -> bool {
    match std::env::var_os(var) {
        Some(v) => v != "0",
        None => default,
    }
}

fn get(flag: &AtomicU8, var: &str, default: bool) -> bool {
    match flag.load(Ordering::Relaxed) {
        1 => false,
        2 => true,
        _ => val_from_environment(var, default),
    }
}

/// `setAllowNontetrahedralChirality`. Upstream stores the value in the
/// process environment; here it is kept in-process (setting environment
/// variables is not thread safe).
pub fn set_allow_nontetrahedral_chirality(val: bool) {
    ALLOW_NONTETRAHEDRAL.store(if val { 2 } else { 1 }, Ordering::Relaxed);
}

/// `getAllowNontetrahedralChirality` (default: true)
pub fn get_allow_nontetrahedral_chirality() -> bool {
    get(&ALLOW_NONTETRAHEDRAL, NON_TETRAHEDRAL_STEREO_ENV_VAR, true)
}

/// `setUseLegacyStereoPerception`
pub fn set_use_legacy_stereo_perception(val: bool) {
    USE_LEGACY.store(if val { 2 } else { 1 }, Ordering::Relaxed);
}

/// `getUseLegacyStereoPerception` (default: true)
pub fn get_use_legacy_stereo_perception() -> bool {
    get(&USE_LEGACY, USE_LEGACY_STEREO_ENV_VAR, true)
}

use crate::graphmol::bond::{BondDir, BondType};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// `MolOps::clearSingleBondDirFlags`
pub fn clear_single_bond_dir_flags(mol: &mut ROMol, only_wedge_flags: bool) {
    for i in 0..mol.num_bonds() {
        let b = mol.bond_mut(i);
        if b.bond_type() == BondType::Single {
            if b.bond_dir() == BondDir::Unknown {
                b.props.set_val(cp::_UnknownStereo, 1i32);
            }
            if !only_wedge_flags
                || !matches!(b.bond_dir(), BondDir::EndDownRight | BondDir::EndUpRight)
            {
                b.set_bond_dir(BondDir::None);
            }
        }
    }
}

/// `MolOps::clearDirFlags`
pub fn clear_dir_flags(mol: &mut ROMol, only_wedge_type_bond_dirs: bool) {
    for i in 0..mol.num_bonds() {
        let b = mol.bond_mut(i);
        if matches!(b.bond_dir(), BondDir::Unknown | BondDir::EitherDouble) {
            b.props.set_val(cp::_UnknownStereo, 1i32);
        }
        if !only_wedge_type_bond_dirs
            || !matches!(b.bond_dir(), BondDir::EndDownRight | BondDir::EndUpRight)
        {
            b.set_bond_dir(BondDir::None);
        }
    }
}

/// `Chirality::removeNonExplicit3DChirality`: clears the tags that were
/// derived from a 3D structure rather than marked in the mol block.
pub fn remove_non_explicit_3d_chirality(mol: &mut ROMol) {
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom_mut(idx);
        if atom.props.has_val(cp::_NonExplicit3DChirality) {
            atom.clear_prop(cp::_NonExplicit3DChirality);
            atom.set_chiral_tag(crate::graphmol::atom::ChiralType::Unspecified);
        }
    }
}
