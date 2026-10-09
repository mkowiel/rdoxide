//! Single-bond directions around double bonds
//! (`MolOps::setDoubleBondNeighborDirections` and helpers) and
//! `Chirality::detail::setStereoForBond`.
//!
//! Ported from RDKit `Code/GraphMol/Chirality.cpp` (2026.09.1).
//!
//! Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::legacy::should_detect_double_bond_stereo;
use crate::geometry::Point3D;
use crate::geometry::point::compute_dihedral_angle;
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `Chirality::detail::setStereoForBond`
pub fn set_stereo_for_bond(
    mol: &mut ROMol,
    bond: u32,
    stereo: BondStereo,
    use_cx_smiles_ordering: bool,
) -> Result<()> {
    let b = mol.bond(bond);
    let (mut beg, mut end) = (b.begin_atom_idx(), b.end_atom_idx());
    if beg > end {
        std::mem::swap(&mut beg, &mut end);
    }
    if mol.atom_degree(beg) > 1 && mol.atom_degree(end) > 1 {
        let mut beg_control = mol.num_atoms();
        for nbr in mol.atom_neighbors(beg) {
            if nbr != end {
                beg_control = beg_control.min(nbr);
            }
        }
        let mut end_control = if use_cx_smiles_ordering {
            mol.num_atoms()
        } else {
            0
        };
        for nbr in mol.atom_neighbors(end) {
            if nbr == beg {
                continue;
            }
            end_control = if use_cx_smiles_ordering {
                end_control.min(nbr)
            } else {
                end_control.max(nbr)
            };
        }
        if beg != mol.bond(bond).begin_atom_idx() {
            std::mem::swap(&mut beg_control, &mut end_control);
        }
        mol.set_bond_stereo_atoms(bond, beg_control, end_control)?;
        mol.bond_mut(bond).set_stereo(stereo)?;
        mol.set_prop("_needsDetectBondStereo", 1i32, false)?;
    }
    Ok(())
}

fn is_dir(d: BondDir) -> bool {
    matches!(d, BondDir::EndDownRight | BondDir::EndUpRight)
}

fn opposite(d: BondDir) -> BondDir {
    if d == BondDir::EndUpRight {
        BondDir::EndDownRight
    } else {
        BondDir::EndUpRight
    }
}

/// anonymous-namespace `setBondDirRelativeToAtom`
fn set_bond_dir_relative_to_atom(
    mol: &mut ROMol,
    bond: u32,
    atom: u32,
    mut dir: BondDir,
    mut reverse: bool,
) -> Result<()> {
    if !is_dir(dir) {
        return Err(RDKitError::precondition("bad dir"));
    }
    if mol.bond(bond).begin_atom_idx() != atom {
        reverse = !reverse;
    }
    if reverse {
        dir = opposite(dir);
    }
    mol.bond_mut(bond).set_bond_dir(dir);
    Ok(())
}

/// anonymous-namespace `isLinearArrangement`
fn is_linear_arrangement(v1: &Point3D, v2: &Point3D) -> bool {
    let lsq = v1.length_sq() * v2.length_sq();
    if lsq < 1.0e-6 {
        return true;
    }
    let dot = v1.dot_product(v2);
    let cos178 = -0.999388;
    dot < cos178 * lsq.sqrt()
}

/// anonymous-namespace `controllingBondFromAtom`; returns (bond, obond).
#[allow(clippy::too_many_arguments)]
fn controlling_bond_from_atom(
    mol: &ROMol,
    needs_dir: &[bool],
    single_bond_counts: &[u32],
    dbl_bond: u32,
    atom: u32,
    squiggle_bond_seen: &mut bool,
    double_bond_seen: &mut bool,
) -> (Option<u32>, Option<u32>) {
    let mut bond: Option<u32> = None;
    let mut obond: Option<u32> = None;
    for t in mol.atom_bond_indices(atom) {
        if t == dbl_bond {
            continue;
        }
        let tb = mol.bond(t);
        let single_like = matches!(tb.bond_type(), BondType::Single | BondType::Aromatic);
        if single_like && (tb.bond_dir() == BondDir::None || is_dir(tb.bond_dir())) {
            match bond {
                None => bond = Some(t),
                Some(cur) => {
                    if needs_dir[t as usize] {
                        if single_bond_counts[t as usize] > single_bond_counts[cur as usize] {
                            obond = Some(cur);
                            bond = Some(t);
                        } else {
                            obond = Some(t);
                        }
                    } else {
                        obond = Some(cur);
                        bond = Some(t);
                    }
                }
            }
        } else if tb.bond_type() == BondType::Double {
            *double_bond_seen = true;
        }
        if single_like
            && (tb.bond_dir() == BondDir::Unknown
                || tb
                    .get_prop(cp::_UnknownStereo)
                    .and_then(RDValue::as_int)
                    .is_some_and(|v| v != 0))
        {
            *squiggle_bond_seen = true;
            break;
        }
    }
    (bond, obond)
}

struct DirState<'a> {
    conf: Option<&'a [Point3D]>,
    needs_dir: Vec<bool>,
    single_bond_counts: Vec<u32>,
    single_bond_nbrs: Vec<Vec<u32>>,
}

/// anonymous-namespace `updateDoubleBondNeighbors`
fn update_double_bond_neighbors(mol: &mut ROMol, dbl: u32, st: &mut DirState) -> Result<()> {
    if mol.bond(dbl).bond_type() != BondType::Double {
        return Err(RDKitError::precondition("not a double bond"));
    }
    if !st.needs_dir[dbl as usize] {
        return Ok(());
    }
    st.needs_dir[dbl as usize] = false;
    let mut followup = Vec::new();
    let mut squiggle = false;
    let mut double_seen = false;
    let atom1 = mol.bond(dbl).begin_atom_idx();
    let atom2 = mol.bond(dbl).end_atom_idx();

    let (b1, ob1) = controlling_bond_from_atom(
        mol,
        &st.needs_dir,
        &st.single_bond_counts,
        dbl,
        atom1,
        &mut squiggle,
        &mut double_seen,
    );
    if squiggle {
        return set_stereo_for_bond(mol, dbl, BondStereo::StereoAny, false);
    }
    let Some(mut bond1) = b1 else {
        return Ok(());
    };
    let mut obond1 = ob1;
    let (b2, ob2) = controlling_bond_from_atom(
        mol,
        &st.needs_dir,
        &st.single_bond_counts,
        dbl,
        atom2,
        &mut squiggle,
        &mut double_seen,
    );
    if squiggle {
        return set_stereo_for_bond(mol, dbl, BondStereo::StereoAny, false);
    }
    let Some(mut bond2) = b2 else {
        return Ok(());
    };
    let mut obond2 = ob2;

    let same_torsion_dir;
    if let Some(pos) = st.conf {
        let begin_p = pos[atom1 as usize];
        let end_p = pos[atom2 as usize];
        let mut bond1_p = pos[mol.bond(bond1).other_atom_idx(atom1)? as usize];
        let mut bond2_p = pos[mol.bond(bond2).other_atom_idx(atom2)? as usize];
        let mut linear = false;
        let mut p1 = bond1_p - begin_p;
        let mut p2 = end_p - begin_p;
        if is_linear_arrangement(&p1, &p2) {
            match obond1 {
                None => linear = true,
                Some(o) => {
                    obond1 = Some(bond1);
                    bond1 = o;
                    bond1_p = pos[mol.bond(bond1).other_atom_idx(atom1)? as usize];
                    p1 = bond1_p - begin_p;
                    if is_linear_arrangement(&p1, &p2) {
                        linear = true;
                    }
                }
            }
        }
        if !linear {
            p1 = bond2_p - end_p;
            p2 = begin_p - end_p;
            if is_linear_arrangement(&p1, &p2) {
                match obond2 {
                    None => linear = true,
                    Some(o) => {
                        obond2 = Some(bond2);
                        bond2 = o;
                        bond2_p = pos[mol.bond(bond2).other_atom_idx(atom2)? as usize];
                        // note: upstream measures from beginP here
                        p1 = bond2_p - begin_p;
                        if is_linear_arrangement(&p1, &p2) {
                            linear = true;
                        }
                    }
                }
            }
        }
        if linear {
            return set_stereo_for_bond(mol, dbl, BondStereo::StereoAny, false);
        }
        let ang = compute_dihedral_angle(&bond1_p, &begin_p, &end_p, &bond2_p);
        same_torsion_dir = ang >= PI / 2.0;
    } else {
        let mut s = match mol.bond(dbl).stereo() {
            BondStereo::StereoCis | BondStereo::StereoZ => false,
            BondStereo::StereoTrans | BondStereo::StereoE => true,
            _ => return Ok(()),
        };
        let sa = mol.bond(dbl).stereo_atoms().to_vec();
        if sa.len() < 2 {
            return Err(RDKitError::range("idx", sa.len(), 2));
        }
        let b1a = mol.bond(bond1).other_atom_idx(atom1)? as i32;
        if b1a != sa[0] && b1a != sa[1] {
            s = !s;
        }
        let b2a = mol.bond(bond2).other_atom_idx(atom2)? as i32;
        if b2a != sa[0] && b2a != sa[1] {
            s = !s;
        }
        same_torsion_dir = s;
    }

    let mut reverse = same_torsion_dir;
    for b in [bond1, bond2] {
        if st.needs_dir[b as usize] {
            for &bidx in &st.single_bond_nbrs[b as usize] {
                if st.needs_dir[bidx as usize] {
                    followup.push(bidx);
                }
            }
        }
    }
    if !st.needs_dir[bond1 as usize] {
        if st.needs_dir[bond2 as usize] {
            if mol.bond(bond1).begin_atom_idx() != atom1 {
                reverse = !reverse;
            }
            let d = mol.bond(bond1).bond_dir();
            set_bond_dir_relative_to_atom(mol, bond2, atom2, d, reverse)?;
        }
    } else if !st.needs_dir[bond2 as usize] {
        if mol.bond(bond2).begin_atom_idx() != atom2 {
            reverse = !reverse;
        }
        let d = mol.bond(bond2).bond_dir();
        set_bond_dir_relative_to_atom(mol, bond1, atom1, d, reverse)?;
    } else {
        set_bond_dir_relative_to_atom(mol, bond1, atom1, BondDir::EndDownRight, false)?;
        set_bond_dir_relative_to_atom(mol, bond2, atom2, BondDir::EndDownRight, reverse)?;
    }
    st.needs_dir[bond1 as usize] = false;
    st.needs_dir[bond2 as usize] = false;
    if let Some(o) = obond1 {
        if st.needs_dir[o as usize] {
            let d = mol.bond(bond1).bond_dir();
            let rev = mol.bond(bond1).begin_atom_idx() == atom1;
            set_bond_dir_relative_to_atom(mol, o, atom1, d, rev)?;
            st.needs_dir[o as usize] = false;
        }
    }
    if let Some(o) = obond2 {
        if st.needs_dir[o as usize] {
            let d = mol.bond(bond2).bond_dir();
            let rev = mol.bond(bond2).begin_atom_idx() == atom2;
            set_bond_dir_relative_to_atom(mol, o, atom2, d, rev)?;
            st.needs_dir[o as usize] = false;
        }
    }
    for b in followup {
        update_double_bond_neighbors(mol, b, st)?;
    }
    Ok(())
}

/// anonymous-namespace `isBondCandidateForStereo`
fn is_bond_candidate_for_stereo(mol: &ROMol, bond: u32) -> bool {
    let b = mol.bond(bond);
    b.bond_type() == BondType::Double
        && b.stereo() != BondStereo::StereoAny
        && b.bond_dir() != BondDir::EitherDouble
        && mol.atom_degree(b.begin_atom_idx()) > 1
        && mol.atom_degree(b.end_atom_idx()) > 1
        && should_detect_double_bond_stereo(mol, bond)
}

/// `MolOps::setDoubleBondNeighborDirections(mol, conf)`; `conf` is the
/// position of the conformer in `mol.conformers()`.
pub fn set_double_bond_neighbor_directions(mol: &mut ROMol, conf: Option<usize>) -> Result<()> {
    let nb = mol.num_bonds() as usize;
    let mut single_bond_counts = vec![0u32; nb];
    let mut bonds_in_play = Vec::new();
    let mut dbl_bond_nbrs: Vec<Vec<u32>> = vec![Vec::new(); nb];
    let mut single_bond_nbrs: Vec<Vec<u32>> = vec![Vec::new(); nb];
    let mut needs_dir = vec![false; nb];

    if !mol.ring_info().is_symm_sssr() {
        symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    }
    for bi in 0..mol.num_bonds() {
        if !is_bond_candidate_for_stereo(mol, bi) {
            continue;
        }
        let mut is_candidate = true;
        let ends = [mol.bond(bi).begin_atom_idx(), mol.bond(bi).end_atom_idx()];
        'outer: for bond_atom in ends {
            for nbi in mol.atom_bond_indices(bond_atom) {
                let nbr = mol.bond(nbi);
                if matches!(nbr.bond_type(), BondType::Single | BondType::Aromatic) {
                    single_bond_counts[nbi as usize] += 1;
                    let nbr_dir = nbr.bond_dir();
                    if nbr.begin_atom_idx() == bond_atom
                        && nbr_dir == BondDir::Unknown
                        && nbr
                            .get_prop(cp::_UnknownStereo)
                            .and_then(RDValue::as_int)
                            .is_some_and(|v| v != 0)
                    {
                        is_candidate = false;
                    } else {
                        needs_dir[bi as usize] = true;
                        if nbr_dir == BondDir::None || is_dir(nbr_dir) {
                            needs_dir[nbi as usize] = true;
                            dbl_bond_nbrs[bi as usize].push(nbi);
                            if !single_bond_nbrs[nbi as usize].contains(&bi) {
                                single_bond_nbrs[nbi as usize].push(bi);
                            }
                        }
                    }
                }
                if !is_candidate {
                    break 'outer;
                }
            }
        }
        if is_candidate {
            bonds_in_play.push(bi);
        }
    }
    if bonds_in_play.is_empty() {
        return Ok(());
    }
    let mut ordered: Vec<(u32, u32)> = bonds_in_play
        .iter()
        .map(|&b| {
            // note: upstream accumulates the neighbor *bond indices*
            let mut count: u32 = dbl_bond_nbrs[b as usize].iter().sum();
            if mol.ring_info().num_bond_rings(b).unwrap_or(0) == 0 {
                count *= 10;
            }
            (count, b)
        })
        .collect();
    ordered.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));

    let positions: Option<Vec<Point3D>> = conf.map(|c| mol.conformers()[c].positions().to_vec());
    let mut st = DirState {
        conf: positions.as_deref(),
        needs_dir,
        single_bond_counts,
        single_bond_nbrs,
    };
    for (_, b) in ordered {
        update_double_bond_neighbors(mol, b, &mut st)?;
    }
    Ok(())
}

/// `Chirality::hasStereoBondDir`
pub fn has_stereo_bond_dir(bond: &Bond) -> bool {
    is_dir(bond.bond_dir())
}

/// `Chirality::getNeighboringDirectedBond`: the first non-double bond of
/// `atom` that carries a stereo direction.
pub fn neighboring_directed_bond(mol: &ROMol, atom: u32) -> Option<u32> {
    mol.atom_bond_indices(atom).find(|&b| {
        let bd = mol.bond(b);
        bd.bond_type() != BondType::Double && is_dir(bd.bond_dir())
    })
}

/// `findHighestCIPNeighbor`: the neighbor of `atom` (other than `skip`) with
/// the highest `_CIPRank`; `None` if a rank is missing or the best is tied.
fn highest_cip_neighbor(mol: &ROMol, atom: u32, skip: u32) -> Option<u32> {
    let mut best_rank = 0;
    let mut best: Option<u32> = None;
    for nbr in mol.atom_neighbors(atom).filter(|&n| n != skip) {
        let cip = mol
            .atom(nbr)
            .get_prop(cp::_CIPRank)
            .and_then(RDValue::as_uint)?;
        if cip > best_rank || best.is_none() {
            best_rank = cip;
            best = Some(nbr);
        } else if cip == best_rank {
            crate::rd_log_warning!(
                "Warning: duplicate CIP ranks found in findHighestCIPNeighbor()"
            );
            best = None;
        }
    }
    best
}

/// `Chirality::findStereoAtoms`: the stereo atoms of a double bond with
/// defined stereo, taken from the bond or derived from E/Z labels and CIP
/// ranks (empty when they cannot be determined).
pub fn find_stereo_atoms(mol: &ROMol, bond_idx: u32) -> Result<Vec<i32>> {
    let bond = mol.bond(bond_idx);
    if bond.bond_type() != BondType::Double {
        return Err(RDKitError::precondition("not double bond"));
    }
    if bond.stereo() <= BondStereo::StereoAny {
        return Err(RDKitError::precondition("no defined stereo"));
    }
    if !bond.stereo_atoms().is_empty() {
        return Ok(bond.stereo_atoms().to_vec());
    }
    if matches!(bond.stereo(), BondStereo::StereoE | BondStereo::StereoZ) {
        let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        return Ok(
            match (
                highest_cip_neighbor(mol, begin, end),
                highest_cip_neighbor(mol, end, begin),
            ) {
                (Some(a), Some(b)) => vec![a as i32, b as i32],
                _ => Vec::new(),
            },
        );
    }
    crate::rd_log_warning!("Unable to assign stereo atoms for bond {}", bond_idx);
    Ok(Vec::new())
}

/// `MolOps::setBondStereoFromDirections`
pub fn set_bond_stereo_from_directions(mol: &mut ROMol) -> Result<()> {
    mol.clear_prop("_needsDetectBondStereo");
    for bi in 0..mol.num_bonds() {
        let b = mol.bond(bi);
        if b.bond_type() != BondType::Double || b.stereo() == BondStereo::StereoAny {
            continue;
        }
        let (beg, end) = (b.begin_atom_idx(), b.end_atom_idx());
        let (Some(at_begin), Some(at_end)) = (
            neighboring_directed_bond(mol, beg),
            neighboring_directed_bond(mol, end),
        ) else {
            continue;
        };
        let begin_side_atom = mol.bond(at_begin).other_atom_idx(beg)?;
        let end_side_atom = mol.bond(at_end).other_atom_idx(end)?;
        let mut begin_dir = mol.bond(at_begin).bond_dir();
        if mol.bond(at_begin).begin_atom_idx() == beg {
            begin_dir = opposite(begin_dir);
        }
        let mut end_dir = mol.bond(at_end).bond_dir();
        if mol.bond(at_end).end_atom_idx() == end {
            end_dir = opposite(end_dir);
        }
        mol.set_bond_stereo_atoms(bi, begin_side_atom, end_side_atom)?;
        mol.bond_mut(bi).set_stereo(if begin_dir == end_dir {
            BondStereo::StereoTrans
        } else {
            BondStereo::StereoCis
        })?;
    }
    Ok(())
}
