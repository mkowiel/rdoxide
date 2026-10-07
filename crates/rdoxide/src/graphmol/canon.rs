//! Canonical depth-first traversal used by the SMILES writer
//! (`Canon::canonicalizeFragment` and helpers) and enhanced-stereo
//! canonicalization.
//!
//! Ported from RDKit `Code/GraphMol/Canon.{h,cpp}` (2026.09.1).
//!
//! Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Upstream's molecular stack holds `Atom*`/`Bond*`; here it holds indices.

use std::collections::VecDeque;

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::legacy::{
    assign_stereochemistry, is_atom_potential_tetrahedral_center,
};
use crate::graphmol::chirality::nontetrahedral::{
    get_chiral_permutation, get_max_nbors, has_non_tetrahedral_stereo,
};
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::new_canon::{RankOptions, rank_mol_atoms};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::ops::chiral_atom_needs_tag_inversion;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::random::with_global_generator;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// used in the canonical traversal code
pub const MAX_NATOMS: i64 = 5000;
/// used in the canonical traversal code
pub const MAX_CYCLES: usize = 1024;
/// used in the canonical traversal code
pub const MAX_BONDTYPE: i64 = 32;

/// `Canon::AtomColors`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomColor {
    /// not visited
    White,
    /// visited, but not finished
    Grey,
    /// visited and finished
    Black,
}

/// `Canon::MolStackElem`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolStackElem {
    Atom(u32),
    /// a bond and the index of the atom traversed before it
    Bond(u32, i32),
    /// a ring closure (ring index)
    Ring(i32),
    BranchOpen(i32),
    BranchClose(i32),
}

pub type MolStack = Vec<MolStackElem>;

/// `gboost::hash_range` over the bytes of a string (32-bit result; `char` is
/// hashed as a signed value, as on the platforms RDKit is built for).
pub fn gboost_hash_str(s: &str) -> u32 {
    let mut seed: u32 = 0;
    for &b in s.as_bytes() {
        let h = b as i8 as i32 as u32;
        seed ^= h
            .wrapping_add(0x9e37_79b9)
            .wrapping_add(seed << 6)
            .wrapping_add(seed >> 2);
    }
    seed
}

fn is_canonicalizable_stereo_double_bond(mol: &ROMol, bond: u32) -> bool {
    let b = mol.bond(bond);
    let s = b.stereo();
    b.bond_type() == BondType::Double
        && s >= BondStereo::StereoZ
        && s <= BondStereo::StereoTrans
        && b.stereo_atoms().len() == 2
}

fn flip_stereo_bond_dir(d: BondDir) -> BondDir {
    match d {
        BondDir::EndUpRight => BondDir::EndDownRight,
        BondDir::EndDownRight => BondDir::EndUpRight,
        other => other,
    }
}

/// `isDative(bond)`
pub fn is_dative(t: BondType) -> bool {
    matches!(
        t,
        BondType::Dative | BondType::DativeL | BondType::DativeR | BondType::DativeOne
    )
}

/// `canSetDoubleBondStereo`
pub fn can_set_double_bond_stereo(t: BondType) -> bool {
    matches!(t, BondType::Single | BondType::Aromatic) || is_dative(t)
}

/// `canHaveDirection`
pub fn can_have_direction(t: BondType) -> bool {
    matches!(t, BondType::Single | BondType::Aromatic)
}

fn log_inconsistent_bond_dirs_warning(idx: u32) {
    crate::rd_log_warning!(
        "Conflicting single bond directions around double bond at index {idx}.\n  Stereochemistry may be incorrect."
    );
}

fn set_direction_from_neighboring_bond(
    mol: &mut ROMol,
    source: u32,
    source_flipped: bool,
    target: u32,
    target_flipped: bool,
) {
    let mut dir = mol.bond(source).bond_dir();
    if source_flipped == target_flipped {
        dir = flip_stereo_bond_dir(dir);
    }
    mol.bond_mut(target).set_bond_dir(dir);
}

#[allow(clippy::too_many_arguments)]
fn get_reference_direction(
    mol: &ROMol,
    dbl: u32,
    ref_atom: u32,
    target_atom: u32,
    ref_controlling_bond: u32,
    ref_is_flipped: bool,
    target_bond: u32,
    target_is_flipped: bool,
) -> Result<BondDir> {
    let d = mol.bond(dbl);
    let ref_dir = mol.bond(ref_controlling_bond).bond_dir();
    let mut dir = match d.stereo() {
        BondStereo::StereoE | BondStereo::StereoTrans => ref_dir,
        BondStereo::StereoZ | BondStereo::StereoCis => flip_stereo_bond_dir(ref_dir),
        _ => BondDir::None,
    };
    if dir == BondDir::None {
        return Err(RDKitError::invariant("stereo not set"));
    }
    let stereo_atoms = d.stereo_atoms();
    let ref_other = mol.bond(ref_controlling_bond).other_atom_idx(ref_atom)? as i32;
    if mol.atom_degree(ref_atom) == 3 && !stereo_atoms.contains(&ref_other) {
        dir = flip_stereo_bond_dir(dir);
    }
    let tgt_other = mol.bond(target_bond).other_atom_idx(target_atom)? as i32;
    if mol.atom_degree(target_atom) == 3 && !stereo_atoms.contains(&tgt_other) {
        dir = flip_stereo_bond_dir(dir);
    }
    if ref_is_flipped != target_is_flipped {
        dir = flip_stereo_bond_dir(dir);
    }
    Ok(dir)
}

/// `fixConflictAcrossDoubleBond`.
///
/// Upstream iterates over `{make_pair(firstBond, ..), make_pair(secondBond, ..)}`,
/// which *copies* the bonds, so its `&bond == &firstBond` test is never true
/// and `otherBond` is always `firstBond`. That behaviour is reproduced here.
#[allow(clippy::too_many_arguments)]
fn fix_conflict_across_double_bond(
    mol: &ROMol,
    dbl: u32,
    atom: u32,
    first: u32,
    first_flipped: bool,
    second: u32,
    second_flipped: bool,
    ref_atom: u32,
    ref_bond: u32,
    ref_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> Result<bool> {
    for (bond, is_flipped) in [(first, first_flipped), (second, second_flipped)] {
        let other_bond = first;
        let other_idx = mol.bond(other_bond).other_atom_idx(atom)?;
        if atom_dir_counts[other_idx as usize] != 2 {
            continue;
        }
        let expected = get_reference_direction(
            mol,
            dbl,
            ref_atom,
            atom,
            ref_bond,
            ref_flipped,
            bond,
            is_flipped,
        )?;
        if expected == mol.bond(bond).bond_dir() {
            bond_dir_counts[other_bond as usize] = 0;
            atom_dir_counts[atom as usize] -= 1;
            atom_dir_counts[other_idx as usize] -= 1;
            return Ok(true);
        }
    }
    Ok(false)
}

/// One end of a stereo double bond: (atom, first bond, flipped, second bond, flipped)
#[derive(Clone, Copy)]
struct DblEnd {
    atom: u32,
    first: u32,
    first_flipped: bool,
    second: Option<u32>,
    second_flipped: bool,
}

/// `handleDirConflictsAcrossDoubleBond`
fn handle_dir_conflicts_across_double_bond(
    mol: &ROMol,
    dbl: u32,
    e1: DblEnd,
    atom1_consistent: bool,
    e2: DblEnd,
    atom2_consistent: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> Result<bool> {
    if atom1_consistent && atom2_consistent {
        let expected = get_reference_direction(
            mol,
            dbl,
            e1.atom,
            e2.atom,
            e1.first,
            e1.first_flipped,
            e2.first,
            e2.first_flipped,
        )?;
        Ok(expected == mol.bond(e2.first).bond_dir())
    } else if !atom2_consistent && atom1_consistent {
        let s2 = e2
            .second
            .ok_or_else(|| RDKitError::precondition("bad second bond"))?;
        fix_conflict_across_double_bond(
            mol,
            dbl,
            e2.atom,
            e2.first,
            e2.first_flipped,
            s2,
            e2.second_flipped,
            e1.atom,
            e1.first,
            e1.first_flipped,
            bond_dir_counts,
            atom_dir_counts,
        )
    } else if !atom1_consistent && atom2_consistent {
        let s1 = e1
            .second
            .ok_or_else(|| RDKitError::precondition("bad second bond"))?;
        fix_conflict_across_double_bond(
            mol,
            dbl,
            e1.atom,
            e1.first,
            e1.first_flipped,
            s1,
            e1.second_flipped,
            e2.atom,
            e2.first,
            e2.first_flipped,
            bond_dir_counts,
            atom_dir_counts,
        )
    } else {
        let s1 = e1
            .second
            .ok_or_else(|| RDKitError::precondition("bad second bond"))?;
        let s2 = e2
            .second
            .ok_or_else(|| RDKitError::precondition("bad second bond"))?;
        let a1_bonds = [(e1.first, e1.first_flipped), (s1, e1.second_flipped)];
        let a2_bonds = [(e2.first, e2.first_flipped), (s2, e2.second_flipped)];
        for &(b1, f1) in &a1_bonds {
            for &(b2, f2) in &a2_bonds {
                let expected = get_reference_direction(mol, dbl, e1.atom, e2.atom, b1, f1, b2, f2)?;
                if expected != mol.bond(b2).bond_dir() {
                    continue;
                }
                let a1_other_bond = if b1 == e1.first { s1 } else { e1.first };
                let a1_other_idx = mol.bond(a1_other_bond).other_atom_idx(e1.atom)?;
                if atom_dir_counts[a1_other_idx as usize] != 2 {
                    continue;
                }
                let a2_other_bond = if b2 == e2.first { s2 } else { e2.first };
                let a2_other_idx = mol.bond(a2_other_bond).other_atom_idx(e2.atom)?;
                if a1_other_idx == a2_other_idx {
                    continue;
                }
                if atom_dir_counts[a2_other_idx as usize] != 2 {
                    continue;
                }
                bond_dir_counts[a1_other_bond as usize] = 0;
                atom_dir_counts[e1.atom as usize] -= 1;
                atom_dir_counts[a1_other_idx as usize] -= 1;
                bond_dir_counts[a2_other_bond as usize] = 0;
                atom_dir_counts[e2.atom as usize] -= 1;
                atom_dir_counts[a2_other_idx as usize] -= 1;
                return Ok(true);
            }
        }
        Ok(false)
    }
}

fn same_side_dirs_are_compatible(
    mol: &ROMol,
    first: u32,
    second: u32,
    first_flipped: bool,
    second_flipped: bool,
) -> bool {
    let should_match = first_flipped != second_flipped;
    let matches = mol.bond(first).bond_dir() == mol.bond(second).bond_dir();
    matches == should_match
}

/// `Canon::canonicalizeDoubleBond`
fn canonicalize_double_bond(
    mol: &mut ROMol,
    dbl: u32,
    bond_visit_orders: &[u32],
    atom_visit_orders: &[u32],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> Result<()> {
    let (b0, e0) = (mol.bond(dbl).begin_atom_idx(), mol.bond(dbl).end_atom_idx());
    if mol.bond(dbl).bond_type() != BondType::Double
        || mol.bond(dbl).stereo() <= BondStereo::StereoAny
        || mol.bond(dbl).stereo_atoms().len() < 2
        || (atom_visit_orders[b0 as usize] == 0 && atom_visit_orders[e0 as usize] == 0)
    {
        return Err(RDKitError::precondition(
            "bad bond in canonicalizeDoubleBond",
        ));
    }
    let (mut atom1, mut atom2) = (b0, e0);
    let deg_ok = |d: u32| d == 2 || d == 3;
    if !deg_ok(mol.atom_degree(atom1)) || !deg_ok(mol.atom_degree(atom2)) {
        return Ok(());
    }
    if atom_visit_orders[b0 as usize] >= atom_visit_orders[e0 as usize] {
        std::mem::swap(&mut atom1, &mut atom2);
    }

    let find_neighbor_bonds = |mol: &ROMol, atom: u32, counts: &[i8]| {
        let mut first: Option<u32> = None;
        let mut second: Option<u32> = None;
        let mut dir_set = false;
        let mut first_visit = mol.num_bonds() + 1;
        for b in mol.atom_bond_indices(atom) {
            if b == dbl || !can_set_double_bond_stereo(mol.bond(b).bond_type()) {
                continue;
            }
            if counts[b as usize] > 0 {
                dir_set = true;
            }
            if first.is_none() || bond_visit_orders[b as usize] < first_visit {
                if first.is_some() {
                    second = first;
                }
                first = Some(b);
                first_visit = bond_visit_orders[b as usize];
            } else {
                second = Some(b);
            }
        }
        (first, second, dir_set)
    };
    let (f1, second_from_atom1, dir1_set) = find_neighbor_bonds(mol, atom1, bond_dir_counts);
    let (f2, second_from_atom2, dir2_set) = find_neighbor_bonds(mol, atom2, bond_dir_counts);
    let (Some(first_from_atom1), Some(first_from_atom2)) = (f1, f2) else {
        return Ok(());
    };

    let has_rc = |mol: &ROMol, b: u32| mol.bond(b).has_prop(cp::_TraversalRingClosureBond);
    let anchor_after = |mol: &ROMol, b: u32, at: u32| -> Result<bool> {
        let anchor = mol.bond(b).other_atom_idx(at)?;
        Ok((atom_visit_orders[at as usize] < atom_visit_orders[anchor as usize]) != has_rc(mol, b))
    };
    let anchor_before = |mol: &ROMol, b: u32, at: u32| -> Result<bool> {
        let anchor = mol.bond(b).other_atom_idx(at)?;
        Ok((atom_visit_orders[anchor as usize] < atom_visit_orders[at as usize]) != has_rc(mol, b))
    };
    let is_first_from_atom1_flipped = anchor_after(mol, first_from_atom1, atom1)?;
    let is_second_from_atom1_flipped = match second_from_atom1 {
        Some(b) => anchor_after(mol, b, atom1)?,
        None => false,
    };
    let is_first_from_atom2_flipped = anchor_before(mol, first_from_atom2, atom2)?;
    let is_second_from_atom2_flipped = match second_from_atom2 {
        Some(b) => anchor_before(mol, b, atom2)?,
        None => false,
    };

    if dir1_set && dir2_set {
        let mut atom1_consistent = true;
        if let Some(s1) = second_from_atom1 {
            if bond_dir_counts[first_from_atom1 as usize] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    s1,
                    is_second_from_atom1_flipped,
                    first_from_atom1,
                    is_first_from_atom1_flipped,
                );
            } else if bond_dir_counts[s1 as usize] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    first_from_atom1,
                    is_first_from_atom1_flipped,
                    s1,
                    is_second_from_atom1_flipped,
                );
            } else {
                atom1_consistent = same_side_dirs_are_compatible(
                    mol,
                    first_from_atom1,
                    s1,
                    is_first_from_atom1_flipped,
                    is_second_from_atom1_flipped,
                );
            }
            bond_dir_counts[s1 as usize] += 1;
            atom_dir_counts[atom1 as usize] += 1;
        }
        bond_dir_counts[first_from_atom1 as usize] += 1;
        atom_dir_counts[atom1 as usize] += 1;

        let mut atom2_consistent = true;
        if let Some(s2) = second_from_atom2 {
            if bond_dir_counts[first_from_atom2 as usize] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    s2,
                    is_second_from_atom2_flipped,
                    first_from_atom2,
                    is_first_from_atom2_flipped,
                );
            } else if bond_dir_counts[s2 as usize] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    first_from_atom2,
                    is_first_from_atom2_flipped,
                    s2,
                    is_second_from_atom2_flipped,
                );
            } else {
                atom2_consistent = same_side_dirs_are_compatible(
                    mol,
                    first_from_atom2,
                    s2,
                    is_first_from_atom2_flipped,
                    is_second_from_atom2_flipped,
                );
            }
            bond_dir_counts[s2 as usize] += 1;
            atom_dir_counts[atom2 as usize] += 1;
        }
        bond_dir_counts[first_from_atom2 as usize] += 1;
        atom_dir_counts[atom2 as usize] += 1;

        let e1 = DblEnd {
            atom: atom1,
            first: first_from_atom1,
            first_flipped: is_first_from_atom1_flipped,
            second: second_from_atom1,
            second_flipped: is_second_from_atom1_flipped,
        };
        let e2 = DblEnd {
            atom: atom2,
            first: first_from_atom2,
            first_flipped: is_first_from_atom2_flipped,
            second: second_from_atom2,
            second_flipped: is_second_from_atom2_flipped,
        };
        if !handle_dir_conflicts_across_double_bond(
            mol,
            dbl,
            e1,
            atom1_consistent,
            e2,
            atom2_consistent,
            bond_dir_counts,
            atom_dir_counts,
        )? {
            log_inconsistent_bond_dirs_warning(dbl);
        }
        return Ok(());
    }

    let mut set_from_bond1 = true;
    let mut atom1_controlling = first_from_atom1;
    let mut atom2_controlling = first_from_atom2;
    if !dir1_set && !dir2_set {
        mol.bond_mut(first_from_atom1)
            .set_bond_dir(BondDir::EndUpRight);
        bond_dir_counts[first_from_atom1 as usize] += 1;
        atom_dir_counts[atom1 as usize] += 1;
    } else if !dir2_set {
        if bond_dir_counts[first_from_atom1 as usize] > 0 {
            bond_dir_counts[first_from_atom1 as usize] += 1;
            atom_dir_counts[atom1 as usize] += 1;
            if let Some(s1) = second_from_atom1 {
                if bond_dir_counts[s1 as usize] != 0 {
                    if !same_side_dirs_are_compatible(
                        mol,
                        first_from_atom1,
                        s1,
                        is_first_from_atom1_flipped,
                        is_second_from_atom1_flipped,
                    ) {
                        log_inconsistent_bond_dirs_warning(dbl);
                    }
                    bond_dir_counts[s1 as usize] += 1;
                    atom_dir_counts[atom1 as usize] += 1;
                }
            }
        } else {
            let s1 = second_from_atom1
                .filter(|&s| bond_dir_counts[s as usize] > 0)
                .ok_or_else(|| RDKitError::invariant("inconsistent state"))?;
            set_direction_from_neighboring_bond(
                mol,
                s1,
                is_second_from_atom1_flipped,
                first_from_atom1,
                is_first_from_atom1_flipped,
            );
            bond_dir_counts[s1 as usize] += 1;
            bond_dir_counts[first_from_atom1 as usize] += 1;
            atom_dir_counts[atom1 as usize] += 2;
            atom1_controlling = s1;
        }
    } else {
        set_from_bond1 = false;
        if bond_dir_counts[first_from_atom2 as usize] > 0 {
            bond_dir_counts[first_from_atom2 as usize] += 1;
            atom_dir_counts[atom2 as usize] += 1;
            if let Some(s2) = second_from_atom2 {
                if bond_dir_counts[s2 as usize] != 0 {
                    if !same_side_dirs_are_compatible(
                        mol,
                        first_from_atom2,
                        s2,
                        is_first_from_atom2_flipped,
                        is_second_from_atom2_flipped,
                    ) {
                        log_inconsistent_bond_dirs_warning(dbl);
                    }
                    bond_dir_counts[s2 as usize] += 1;
                    atom_dir_counts[atom2 as usize] += 1;
                }
            }
        } else {
            let s2 = second_from_atom2
                .filter(|&s| bond_dir_counts[s as usize] > 0)
                .ok_or_else(|| RDKitError::invariant("inconsistent state"))?;
            set_direction_from_neighboring_bond(
                mol,
                s2,
                is_second_from_atom2_flipped,
                first_from_atom2,
                is_first_from_atom2_flipped,
            );
            bond_dir_counts[s2 as usize] += 1;
            bond_dir_counts[first_from_atom2 as usize] += 1;
            atom_dir_counts[atom2 as usize] += 2;
            atom2_controlling = s2;
        }
    }

    if set_from_bond1 {
        let ctrl_flipped = if atom1_controlling == first_from_atom1 {
            is_first_from_atom1_flipped
        } else {
            is_second_from_atom1_flipped
        };
        let d = get_reference_direction(
            mol,
            dbl,
            atom1,
            atom2,
            atom1_controlling,
            ctrl_flipped,
            first_from_atom2,
            is_first_from_atom2_flipped,
        )?;
        mol.bond_mut(first_from_atom2).set_bond_dir(d);
        bond_dir_counts[first_from_atom2 as usize] += 1;
        atom_dir_counts[atom2 as usize] += 1;
    } else {
        let ctrl_flipped = if atom2_controlling == first_from_atom2 {
            is_first_from_atom2_flipped
        } else {
            is_second_from_atom2_flipped
        };
        let d = get_reference_direction(
            mol,
            dbl,
            atom2,
            atom1,
            atom2_controlling,
            ctrl_flipped,
            first_from_atom1,
            is_first_from_atom1_flipped,
        )?;
        mol.bond_mut(first_from_atom1).set_bond_dir(d);
        bond_dir_counts[first_from_atom1 as usize] += 1;
        atom_dir_counts[atom1 as usize] += 1;
    }

    if mol.atom_degree(atom1) == 3 {
        if let Some(s1) = second_from_atom1 {
            if bond_dir_counts[s1 as usize] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    first_from_atom1,
                    is_first_from_atom1_flipped,
                    s1,
                    is_second_from_atom1_flipped,
                );
                bond_dir_counts[s1 as usize] += 1;
                atom_dir_counts[atom1 as usize] += 1;
            }
        }
    }
    if mol.atom_degree(atom2) == 3 {
        if let Some(s2) = second_from_atom2 {
            if bond_dir_counts[s2 as usize] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    first_from_atom2,
                    is_first_from_atom2_flipped,
                    s2,
                    is_second_from_atom2_flipped,
                );
                bond_dir_counts[s2 as usize] += 1;
                atom_dir_counts[atom2 as usize] += 1;
            }
        }
    }
    Ok(())
}

/// `Canon::canonicalizeDoubleBonds`
fn canonicalize_double_bonds(
    mol: &mut ROMol,
    bond_visit_orders: &[u32],
    atom_visit_orders: &[u32],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    mol_stack: &MolStack,
) -> Result<()> {
    let nb = mol.num_bonds() as usize;
    let mut stereo_bond_nbrs: Vec<Vec<u32>> = vec![Vec::new(); nb];
    let mut queued: Vec<u32> = Vec::new();

    let neighboring_stereo_bond = |mol: &ROMol, dbl_atom: u32, nbr_bnd: u32| -> Option<u32> {
        let other = mol.bond(nbr_bnd).other_atom_idx(dbl_atom).ok()?;
        mol.atom_bond_indices(other)
            .find(|&b| b != nbr_bnd && is_canonicalizable_stereo_double_bond(mol, b))
    };

    for elem in mol_stack {
        let MolStackElem::Bond(bond, _) = *elem else {
            continue;
        };
        if matches!(
            mol.bond(bond).bond_dir(),
            BondDir::EndDownRight | BondDir::EndUpRight
        ) {
            mol.bond_mut(bond).set_bond_dir(BondDir::None);
        }
        if !is_canonicalizable_stereo_double_bond(mol, bond)
            || mol
                .bond(bond)
                .stereo_atoms()
                .iter()
                .any(|&a| atom_visit_orders[a as usize] == 0)
        {
            let b = mol.bond_mut(bond);
            b.set_stereo(BondStereo::StereoNone)?;
            b.stereo_atoms_mut().clear();
            continue;
        }
        let mut current = Vec::new();
        for dbl_atom in [
            mol.bond(bond).begin_atom_idx(),
            mol.bond(bond).end_atom_idx(),
        ] {
            for nb in mol.atom_bond_indices(dbl_atom) {
                if !can_have_direction(mol.bond(nb).bond_type()) {
                    continue;
                }
                if let Some(n) = neighboring_stereo_bond(mol, dbl_atom, nb) {
                    current.push(n);
                }
            }
        }
        current.sort_by_key(|&b| bond_visit_orders[b as usize]);
        stereo_bond_nbrs[bond as usize] = current;
        queued.push(bond);
    }
    // priority queue order: most stereo neighbors first, then lowest visit order
    queued.sort_by(|&a, &b| {
        stereo_bond_nbrs[b as usize]
            .len()
            .cmp(&stereo_bond_nbrs[a as usize].len())
            .then(bond_visit_orders[a as usize].cmp(&bond_visit_orders[b as usize]))
    });

    let mut seen = vec![false; nb];
    for bond in queued {
        if seen[bond as usize] {
            continue;
        }
        let mut q = VecDeque::new();
        q.push_back(bond);
        while let Some(cur) = q.pop_front() {
            if seen[cur as usize] || bond_visit_orders[cur as usize] == 0 {
                continue;
            }
            if !is_canonicalizable_stereo_double_bond(mol, cur) {
                let b = mol.bond_mut(cur);
                b.set_stereo(BondStereo::StereoNone)?;
                b.stereo_atoms_mut().clear();
                seen[cur as usize] = true;
                continue;
            }
            canonicalize_double_bond(
                mol,
                cur,
                bond_visit_orders,
                atom_visit_orders,
                bond_dir_counts,
                atom_dir_counts,
            )?;
            seen[cur as usize] = true;
            for &n in &stereo_bond_nbrs[cur as usize] {
                if !seen[n as usize] {
                    q.push_back(n);
                }
            }
        }
    }
    Ok(())
}

struct Traversal<'a> {
    ranks: &'a [u32],
    bonds_in_play: Option<&'a [bool]>,
    bond_symbols: Option<&'a [String]>,
    do_random: bool,
}

impl Traversal<'_> {
    fn symbol_hash(&self, bond: u32) -> Option<i64> {
        self.bond_symbols
            .map(|s| (gboost_hash_str(&s[bond as usize]) % MAX_NATOMS as u32) as i64)
    }

    /// The traversal priority of a neighbor; upstream does this arithmetic on
    /// `unsigned int` and stores the result in an `int`.
    fn rank(&self, mol: &ROMol, bond: u32, other: u32, grey: bool) -> i32 {
        if self.do_random {
            return with_global_generator(0, |g| g.next()) as i32;
        }
        let mut rank = self.ranks[other as usize] as i64;
        let bt = mol.bond(bond).bond_type() as u8 as i64;
        if grey {
            rank -= (MAX_BONDTYPE + 1) * MAX_NATOMS * MAX_NATOMS;
            match self.symbol_hash(bond) {
                None => rank += (MAX_BONDTYPE - bt) * MAX_NATOMS,
                Some(h) => rank += h * MAX_NATOMS,
            }
        } else if mol.ring_info().num_bond_rings(bond).unwrap_or(0) > 0 {
            match self.symbol_hash(bond) {
                None => rank += (MAX_BONDTYPE - bt) * MAX_NATOMS * MAX_NATOMS,
                Some(h) => rank += h * MAX_NATOMS * MAX_NATOMS,
            }
        }
        rank as u32 as i32
    }

    fn in_play(&self, bond: u32) -> bool {
        self.bonds_in_play.is_none_or(|b| b[bond as usize])
    }

    /// `Canon::dfsFindCycles`
    fn dfs_find_cycles(
        &self,
        mol: &ROMol,
        atom_idx: u32,
        in_bond_idx: i32,
        colors: &mut [AtomColor],
        atom_ring_closures: &mut [Vec<u32>],
    ) -> Result<()> {
        colors[atom_idx as usize] = AtomColor::Grey;
        let mut possibles: Vec<(i32, u32, u32)> =
            Vec::with_capacity(mol.atom_degree(atom_idx) as usize);
        for b in mol.atom_bond_indices(atom_idx) {
            if !self.in_play(b) {
                continue;
            }
            if in_bond_idx < 0 || b != in_bond_idx as u32 {
                let other = mol.bond(b).other_atom_idx(atom_idx)?;
                let grey = colors[other as usize] == AtomColor::Grey;
                possibles.push((self.rank(mol, b, other, grey), other, b));
            }
        }
        possibles.sort_by_key(|p| p.0);
        for (_, possible_idx, bond) in possibles {
            match colors[possible_idx as usize] {
                AtomColor::White => self.dfs_find_cycles(
                    mol,
                    possible_idx,
                    bond as i32,
                    colors,
                    atom_ring_closures,
                )?,
                AtomColor::Grey => {
                    atom_ring_closures[possible_idx as usize].push(bond);
                    atom_ring_closures[atom_idx as usize].push(bond);
                }
                AtomColor::Black => {}
            }
        }
        colors[atom_idx as usize] = AtomColor::Black;
        Ok(())
    }

    /// `Canon::dfsBuildStack`
    #[allow(clippy::too_many_arguments)]
    fn dfs_build_stack(
        &self,
        mol: &mut ROMol,
        atom_idx: u32,
        in_bond_idx: i32,
        colors: &mut [AtomColor],
        cycles_available: &mut [bool],
        mol_stack: &mut MolStack,
        atom_ring_closures: &[Vec<u32>],
        atom_traversal_bond_order: &mut [Vec<i32>],
    ) -> Result<()> {
        let mut seen_from_here = vec![false; mol.num_atoms() as usize];
        seen_from_here[atom_idx as usize] = true;
        mol_stack.push(MolStackElem::Atom(atom_idx));
        colors[atom_idx as usize] = AtomColor::Grey;

        let mut trav_list: Vec<i32> = Vec::new();
        if in_bond_idx >= 0 {
            trav_list.push(in_bond_idx);
        }

        if !atom_ring_closures[atom_idx as usize].is_empty() {
            let mut rings_closed = Vec::new();
            for &b_idx in &atom_ring_closures[atom_idx as usize] {
                trav_list.push(b_idx as i32);
                let other = mol.bond(b_idx).other_atom_idx(atom_idx)?;
                seen_from_here[other as usize] = true;
                let ring_idx = mol
                    .bond(b_idx)
                    .get_prop(cp::_TraversalRingClosureBond)
                    .and_then(RDValue::as_uint);
                if let Some(ring_idx) = ring_idx {
                    mol_stack.push(MolStackElem::Bond(b_idx, atom_idx as i32));
                    mol_stack.push(MolStackElem::Ring(ring_idx as i32));
                    rings_closed.push(ring_idx as usize - 1);
                } else {
                    let Some(lowest) = cycles_available.iter().position(|&x| x) else {
                        return Err(RDKitError::Value(
                            "Too many rings open at once. SMILES cannot be generated.".into(),
                        ));
                    };
                    cycles_available[lowest] = false;
                    let ring_idx = lowest as u32 + 1;
                    mol.bond_mut(b_idx)
                        .set_prop(cp::_TraversalRingClosureBond, ring_idx, false)?;
                    mol_stack.push(MolStackElem::Ring(ring_idx as i32));
                }
            }
            for r in rings_closed {
                cycles_available[r] = true;
            }
        }

        let mut possibles: Vec<(i32, u32, u32)> =
            Vec::with_capacity(mol.atom_degree(atom_idx) as usize);
        for b in mol.atom_bond_indices(atom_idx) {
            if !self.in_play(b) {
                continue;
            }
            if in_bond_idx < 0 || b != in_bond_idx as u32 {
                let other = mol.bond(b).other_atom_idx(atom_idx)?;
                if colors[other as usize] != AtomColor::White || seen_from_here[other as usize] {
                    continue;
                }
                possibles.push((self.rank(mol, b, other, false), other, b));
            }
        }
        possibles.sort_by_key(|p| p.0);

        let n = possibles.len();
        for (i, &(_, possible_idx, bond)) in possibles.iter().enumerate() {
            if colors[possible_idx as usize] != AtomColor::White {
                continue;
            }
            mol.atom_mut(possible_idx)
                .clear_prop(cp::_TraversalBondIndexOrder);
            trav_list.push(bond as i32);
            if i + 1 != n {
                mol_stack.push(MolStackElem::BranchOpen(i as i32));
            }
            mol_stack.push(MolStackElem::Bond(bond, atom_idx as i32));
            self.dfs_build_stack(
                mol,
                possible_idx,
                bond as i32,
                colors,
                cycles_available,
                mol_stack,
                atom_ring_closures,
                atom_traversal_bond_order,
            )?;
            if i + 1 != n {
                mol_stack.push(MolStackElem::BranchClose(i as i32));
            }
        }
        atom_traversal_bond_order[atom_idx as usize] = trav_list;
        colors[atom_idx as usize] = AtomColor::Black;
        Ok(())
    }
}

/// `Canon::clearBondDirs`
fn clear_bond_dirs(
    mol: &mut ROMol,
    ref_bond: u32,
    from_atom: u32,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> Result<()> {
    let clear_direction =
        |mol: &mut ROMol, bond: u32, bdc: &mut [i8], adc: &mut [i8]| -> Result<()> {
            bdc[bond as usize] -= 1;
            if bdc[bond as usize] == 0 {
                mol.bond_mut(bond).set_bond_dir(BondDir::None);
                adc[from_atom as usize] -= 1;
                let other = mol.bond(bond).other_atom_idx(from_atom)?;
                if adc[other as usize] != 0 {
                    adc[other as usize] -= 1;
                }
            }
            Ok(())
        };
    let bonds: Vec<u32> = mol.atom_bond_indices(from_atom).collect();
    for o in bonds {
        if o != ref_bond && can_have_direction(mol.bond(o).bond_type()) {
            let ob = mol.bond(o);
            let rb = mol.bond(ref_bond);
            if bond_dir_counts[o as usize] >= bond_dir_counts[ref_bond as usize]
                && atom_dir_counts[ob.begin_atom_idx() as usize] != 1
                && atom_dir_counts[ob.end_atom_idx() as usize] != 1
            {
                clear_direction(mol, o, bond_dir_counts, atom_dir_counts)?;
            } else if atom_dir_counts[rb.begin_atom_idx() as usize] != 1
                && atom_dir_counts[rb.end_atom_idx() as usize] != 1
            {
                clear_direction(mol, ref_bond, bond_dir_counts, atom_dir_counts)?;
            }
            break;
        }
    }
    Ok(())
}

/// `Canon::removeUnwantedBondDirSpecs`
fn remove_unwanted_bond_dir_specs(
    mol: &mut ROMol,
    mol_stack: &MolStack,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    bond_visit_orders: &[u32],
) -> Result<()> {
    for elem in mol_stack {
        let MolStackElem::Bond(bond, _) = *elem else {
            continue;
        };
        let b = mol.bond(bond);
        if b.bond_type() != BondType::Double || b.stereo() > BondStereo::StereoAny {
            continue;
        }
        let (first_atom, second_atom) = (b.begin_atom_idx(), b.end_atom_idx());
        if mol.atom_degree(first_atom) == 1 || mol.atom_degree(second_atom) == 1 {
            continue;
        }
        let mut candidates: Vec<u32> = mol
            .atom_bond_indices(first_atom)
            .filter(|&x| bond_dir_counts[x as usize] != 0)
            .collect();
        if candidates.is_empty() {
            continue;
        }
        if atom_dir_counts[first_atom as usize] != 0 {
            candidates.clear();
        }
        let mut on_second = 0;
        for x in mol.atom_bond_indices(second_atom) {
            if bond_dir_counts[x as usize] != 0 {
                candidates.push(x);
                on_second += 1;
            }
        }
        if on_second == 0 {
            continue;
        }
        if atom_dir_counts[second_atom as usize] != 0 {
            continue;
        }
        candidates.sort_by_key(|&x| bond_visit_orders[x as usize]);
        for cand in candidates {
            let cb = mol.bond(cand);
            let other = if cb.begin_atom_idx() == first_atom || cb.end_atom_idx() == first_atom {
                cb.other_atom_idx(first_atom)?
            } else if cb.begin_atom_idx() == second_atom || cb.end_atom_idx() == second_atom {
                cb.other_atom_idx(second_atom)?
            } else {
                return Err(RDKitError::invariant("inconsistent bond ends"));
            };
            if atom_dir_counts[other as usize] == 2 {
                bond_dir_counts[cand as usize] = 0;
                mol.bond_mut(cand).set_bond_dir(BondDir::None);
                atom_dir_counts[other as usize] -= 1;
                break;
            }
        }
    }
    Ok(())
}

/// `Canon::removeRedundantBondDirSpecs`
fn remove_redundant_bond_dir_specs(
    mol: &mut ROMol,
    mol_stack: &MolStack,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> Result<()> {
    let clear_from_atom =
        |mol: &mut ROMol, t_bond: u32, atom: u32, bdc: &mut [i8], adc: &mut [i8]| -> Result<()> {
            if adc[atom as usize] < 2 {
                return Ok(());
            }
            let found = mol.atom_bond_indices(atom).any(|b| {
                b != t_bond
                    && mol.bond(b).bond_type() == BondType::Double
                    && mol.bond(b).stereo() > BondStereo::StereoAny
            });
            if found {
                clear_bond_dirs(mol, t_bond, atom, bdc, adc)?;
            }
            Ok(())
        };
    for elem in mol_stack {
        let MolStackElem::Bond(t_bond, number) = *elem else {
            continue;
        };
        let begin = number as u32;
        let end = mol.bond(t_bond).other_atom_idx(begin)?;
        if can_have_direction(mol.bond(t_bond).bond_type()) && bond_dir_counts[t_bond as usize] != 0
        {
            clear_from_atom(mol, t_bond, begin, bond_dir_counts, atom_dir_counts)?;
            clear_from_atom(mol, t_bond, end, bond_dir_counts, atom_dir_counts)?;
        } else if mol.bond(t_bond).bond_dir() != BondDir::None {
            mol.bond_mut(t_bond).set_bond_dir(BondDir::None);
        }
    }
    Ok(())
}

/// `insertImplicitNbors`
fn insert_implicit_nbors(bonds: &mut Vec<i32>, tag: ChiralType, first_atom: bool) {
    let ref_max = get_max_nbors(tag) as usize;
    if bonds.len() < ref_max {
        let n = ref_max - bonds.len();
        let at = if first_atom { 0 } else { 1.min(bonds.len()) };
        bonds.splice(at..at, std::iter::repeat_n(-1, n));
    }
}

/// Options of `Canon::canonicalizeFragment`.
#[derive(Debug, Clone, Copy)]
pub struct FragmentOptions<'a> {
    pub atoms_in_play: Option<&'a [bool]>,
    pub bonds_in_play: Option<&'a [bool]>,
    pub bond_symbols: Option<&'a [String]>,
    pub do_isomeric_smiles: bool,
    pub do_random: bool,
    pub do_chiral_inversions: bool,
}

impl Default for FragmentOptions<'_> {
    fn default() -> Self {
        FragmentOptions {
            atoms_in_play: None,
            bonds_in_play: None,
            bond_symbols: None,
            do_isomeric_smiles: false,
            do_random: false,
            do_chiral_inversions: true,
        }
    }
}

/// `Canon::canonicalizeFragment`. If `atoms_in_play` is `None`, the atoms in
/// play are derived as in upstream's overload without `atomsInPlay` (all
/// atoms when there is no `bonds_in_play`, else the atoms of those bonds).
pub fn canonicalize_fragment(
    mol: &mut ROMol,
    atom_idx: u32,
    colors: &mut [AtomColor],
    ranks: &[u32],
    mol_stack: &mut MolStack,
    opts: FragmentOptions<'_>,
) -> Result<()> {
    let n_atoms = mol.num_atoms() as usize;
    let derived_atoms: Vec<bool>;
    let atoms_in_play: &[bool] = match (opts.atoms_in_play, opts.bonds_in_play) {
        (Some(a), _) => a,
        (None, None) => {
            derived_atoms = vec![true; n_atoms];
            &derived_atoms
        }
        (None, Some(bip)) => {
            let mut v = vec![false; n_atoms];
            for b in mol.bonds() {
                if bip[b.idx() as usize] {
                    v[b.begin_atom_idx() as usize] = true;
                    v[b.end_atom_idx() as usize] = true;
                }
            }
            derived_atoms = v;
            &derived_atoms
        }
    };
    if colors.len() < n_atoms || ranks.len() < n_atoms {
        return Err(RDKitError::precondition("vector too small"));
    }

    if !mol.has_prop(cp::_StereochemDone) {
        assign_stereochemistry(mol, false, false, false)?;
    }
    if !mol.ring_info().is_symm_sssr() {
        find_sssr(mol, false, false)?;
    }
    mol.atom_mut(atom_idx)
        .set_prop(cp::_TraversalStartPoint, true, false)?;

    let mut atom_ring_closures: Vec<Vec<u32>> = vec![Vec::new(); n_atoms];
    let mut atom_traversal_bond_order: Vec<Vec<i32>> = vec![Vec::new(); n_atoms];
    let trav = Traversal {
        ranks,
        bonds_in_play: opts.bonds_in_play,
        bond_symbols: opts.bond_symbols,
        do_random: opts.do_random,
    };
    // `canonicalDFSTraversal`
    {
        let mut tcolors = colors.to_vec();
        trav.dfs_find_cycles(mol, atom_idx, -1, &mut tcolors, &mut atom_ring_closures)?;
        let mut cycles_available = vec![true; MAX_CYCLES];
        trav.dfs_build_stack(
            mol,
            atom_idx,
            -1,
            colors,
            &mut cycles_available,
            mol_stack,
            &atom_ring_closures,
            &mut atom_traversal_bond_order,
        )?;
    }
    let Some(&MolStackElem::Atom(first_idx)) = mol_stack.first() else {
        return Err(RDKitError::invariant(
            "Corrupted stack. First element should be an atom.",
        ));
    };

    let mut num_swaps_chiral_atoms = vec![false; n_atoms];
    let mut atom_permutation_indices = vec![0i32; n_atoms];
    if opts.do_isomeric_smiles {
        for ai in 0..mol.num_atoms() {
            if !atoms_in_play[ai as usize] {
                continue;
            }
            if mol.atom(ai).chiral_tag() == ChiralType::Unspecified {
                continue;
            }
            if let Some(bip) = opts.bonds_in_play {
                if mol.atom_bond_indices(ai).any(|b| !bip[b as usize]) {
                    mol.atom_mut(ai)
                        .set_prop(cp::_brokenChirality, true, false)?;
                }
            }
            if mol.atom(ai).has_prop(cp::_brokenChirality) {
                continue;
            }
            let non_tet = has_non_tetrahedral_stereo(mol, ai);
            if !(is_atom_potential_tetrahedral_center(mol, ai)? || non_tet) {
                continue;
            }
            let mut perm: i32 = 0;
            if non_tet {
                if let Some(p) = mol
                    .atom(ai)
                    .get_prop(cp::_chiralPermutation)
                    .and_then(RDValue::as_int)
                {
                    perm = p;
                }
            }
            let first_in_part = ai == first_idx;
            let true_order = &atom_traversal_bond_order[ai as usize];
            let mut n_swaps = 0;
            let mut t_order = true_order.clone();
            if (true_order.len() as u32) < mol.atom_degree(ai) {
                for b in mol.atom_bond_indices(ai) {
                    if !true_order.contains(&(b as i32)) {
                        t_order.push(b as i32);
                    }
                }
            }
            let tag = mol.atom(ai).chiral_tag();
            if perm == 0 {
                n_swaps = mol.perturbation_order(ai, &t_order)?;
            } else {
                insert_implicit_nbors(&mut t_order, tag, first_in_part);
                perm = get_chiral_permutation(mol, ai, &t_order, false) as i32;
            }
            if opts.do_chiral_inversions
                && chiral_atom_needs_tag_inversion(
                    mol,
                    ai,
                    first_in_part,
                    atom_ring_closures[ai as usize].len(),
                )
            {
                n_swaps += 1;
            }
            if n_swaps % 2 != 0 {
                num_swaps_chiral_atoms[ai as usize] = true;
            }
            atom_permutation_indices[ai as usize] = perm;
        }
    }

    let mut atom_visit_orders = vec![0u32; n_atoms];
    let mut bond_visit_orders = vec![0u32; mol.num_bonds() as usize];
    for (pos, elem) in mol_stack.iter().enumerate() {
        let pos = pos as u32 + 1;
        match *elem {
            MolStackElem::Atom(a) => atom_visit_orders[a as usize] = pos,
            MolStackElem::Bond(b, _) => {
                bond_visit_orders[b as usize] = pos;
                if matches!(
                    mol.bond(b).bond_dir(),
                    BondDir::EndDownRight | BondDir::EndUpRight
                ) {
                    mol.bond_mut(b).set_bond_dir(BondDir::None);
                }
            }
            _ => {}
        }
    }

    let mut bond_dir_counts = vec![0i8; mol.num_bonds() as usize];
    let mut atom_dir_counts = vec![0i8; n_atoms];
    canonicalize_double_bonds(
        mol,
        &bond_visit_orders,
        &atom_visit_orders,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        mol_stack,
    )?;

    if opts.do_isomeric_smiles {
        let mut ring_stereo_adjusted = vec![false; n_atoms];
        for elem in mol_stack.iter() {
            let MolStackElem::Atom(ai) = *elem else {
                continue;
            };
            let a = mol.atom(ai);
            if a.chiral_tag() == ChiralType::Unspecified || a.has_prop(cp::_brokenChirality) {
                continue;
            }
            if let Some(RDValue::VecInt(ring_stereo_atoms)) = a.get_prop(cp::_ringStereoAtoms) {
                let ring_stereo_atoms = ring_stereo_atoms.clone();
                if !ring_stereo_adjusted[ai as usize] {
                    mol.atom_mut(ai).set_chiral_tag(ChiralType::TetrahedralCCW);
                    ring_stereo_adjusted[ai as usize] = true;
                }
                for nbr_v in ring_stereo_atoms {
                    let nbr = (nbr_v.unsigned_abs() - 1) as usize;
                    if !ring_stereo_adjusted[nbr]
                        && atom_visit_orders[nbr] > atom_visit_orders[ai as usize]
                    {
                        let tag = mol.atom(ai).chiral_tag();
                        mol.atom_mut(nbr as u32).set_chiral_tag(tag);
                        if nbr_v < 0 {
                            mol.atom_mut(nbr as u32).invert_chirality();
                        }
                        if num_swaps_chiral_atoms[ai as usize] != num_swaps_chiral_atoms[nbr] {
                            mol.atom_mut(nbr as u32).invert_chirality();
                        }
                        ring_stereo_adjusted[nbr] = true;
                    }
                }
            } else if let Some(sgidx) = a
                .get_prop("_stereoGroup")
                .and_then(|v| match v {
                    RDValue::UInt64(x) => Some(*x as usize),
                    other => other.as_uint().map(|x| x as usize),
                })
                .filter(|&i| i < mol.stereo_groups().len())
            {
                let swap_it = a.chiral_tag() == ChiralType::TetrahedralCW;
                if swap_it {
                    mol.atom_mut(ai).invert_chirality();
                }
                if swap_it || num_swaps_chiral_atoms[ai as usize] {
                    let others: Vec<u32> = mol.stereo_groups()[sgidx].atoms().to_vec();
                    for at in others {
                        if at != ai {
                            mol.atom_mut(at).invert_chirality();
                        }
                    }
                }
            } else if matches!(
                a.chiral_tag(),
                ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
            ) {
                if num_swaps_chiral_atoms[ai as usize] {
                    mol.atom_mut(ai).invert_chirality();
                }
            } else if atom_permutation_indices[ai as usize] != 0 {
                let p = atom_permutation_indices[ai as usize];
                mol.atom_mut(ai)
                    .set_prop(cp::_chiralPermutation, p, false)?;
            }
        }
    }
    remove_unwanted_bond_dir_specs(
        mol,
        mol_stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        &bond_visit_orders,
    )?;
    remove_redundant_bond_dir_specs(mol, mol_stack, &mut bond_dir_counts, &mut atom_dir_counts)
}

/// `Canon::canonicalizeEnhancedStereo`
pub fn canonicalize_enhanced_stereo(mol: &mut ROMol, atom_ranks: Option<&[u32]>) -> Result<()> {
    if mol.stereo_groups().is_empty() {
        return Ok(());
    }
    let lranks;
    let ranks: &[u32] = match atom_ranks {
        Some(r) => r,
        None => {
            lranks = rank_mol_atoms(mol, RankOptions::default())?;
            &lranks
        }
    };
    let sgs = mol.stereo_groups().to_vec();
    let mut new_sgs: Vec<StereoGroup> = Vec::new();
    for sg in sgs {
        if sg.group_type() == StereoGroupType::Absolute {
            new_sgs.push(sg);
            continue;
        }
        let mut sg_atoms = sg.atoms().to_vec();
        sg_atoms.sort_by_key(|&a| ranks[a as usize]);
        let mut sg_bonds = sg.bonds().to_vec();
        let bond_key = |b: u32| {
            let bd = mol.bond(b);
            let (x, y) = (
                ranks[bd.begin_atom_idx() as usize],
                ranks[bd.end_atom_idx() as usize],
            );
            (x.max(y), x.min(y))
        };
        sg_bonds.sort_by_key(|&b| bond_key(b));
        let found_ref_state = if let Some(&a0) = sg_atoms.first() {
            mol.atom(a0).chiral_tag()
        } else if let Some(&b0) = sg_bonds.first() {
            if mol.bond(b0).stereo() == BondStereo::StereoAtropCCW {
                ChiralType::TetrahedralCCW
            } else {
                ChiralType::TetrahedralCW
            }
        } else {
            ChiralType::TetrahedralCCW
        };
        if found_ref_state != ChiralType::TetrahedralCCW {
            for &a in &sg_atoms {
                mol.atom_mut(a).invert_chirality();
            }
            for &b in &sg_bonds {
                mol.bond_mut(b).invert_chirality();
            }
        }
        if let Some(&a0) = sg_atoms.first() {
            mol.atom_mut(a0)
                .set_prop("_stereoGroup", new_sgs.len() as u64, true)?;
        }
        new_sgs.push(StereoGroup::new(sg.group_type(), sg_atoms, sg_bonds, 0)?);
    }
    mol.set_stereo_groups(new_sgs)
}
