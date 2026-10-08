//! Legacy stereochemistry perception (`Chirality::legacyStereoPerception`
//! and its helpers) and `MolOps::assignStereochemistry`.
//!
//! Ported from RDKit `Code/GraphMol/Chirality.cpp` (2026.09.1).
//!
//! Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeSet, VecDeque};

use super::nontetrahedral::has_non_tetrahedral_stereo;
use crate::graphmol::atom::{ChiralType, HybridizationType, ValenceType};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::conjug_hybrid::atom_has_conjugated_bond;
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, fast_find_rings, symmetrize_sssr};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::query_ops::query_is_atom_bridgehead;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

type IntPair = (i32, i32);

fn int_prop(props: &impl RDProps, key: &str) -> Option<i32> {
    props.get_prop(key).and_then(RDValue::as_int)
}

fn vec_int_prop(props: &impl RDProps, key: &str) -> Option<Vec<i32>> {
    match props.get_prop(key) {
        Some(RDValue::VecInt(v)) => Some(v.clone()),
        _ => None,
    }
}

/// anonymous-namespace `shouldDetectDoubleBondStereo`
pub(crate) fn should_detect_double_bond_stereo(mol: &ROMol, bond: u32) -> bool {
    let ri = mol.ring_info();
    ri.num_bond_rings(bond).unwrap_or(0) == 0
        || ri.min_bond_ring_size(bond).unwrap_or(0) >= super::MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO
}

/// `Chirality::detail::bondAffectsAtomChirality`
pub fn bond_affects_atom_chirality(mol: &ROMol, bond: u32, atom: u32) -> bool {
    let b = mol.bond(bond);
    !(matches!(b.bond_type(), BondType::Unspecified | BondType::Zero)
        || (b.bond_type() == BondType::Dative && b.begin_atom_idx() == atom))
}

/// `Chirality::detail::getAtomNonzeroDegree`
pub fn get_atom_nonzero_degree(mol: &ROMol, atom: u32) -> u32 {
    mol.atom_bond_indices(atom)
        .filter(|&b| bond_affects_atom_chirality(mol, b, atom))
        .count() as u32
}

/// `Chirality::detail::has_protium_neighbor`
fn has_protium_neighbor(mol: &ROMol, atom: u32) -> bool {
    mol.atom_neighbors(atom).any(|n| {
        let a = mol.atom(n);
        a.atomic_num() == 1 && a.isotope() == 0
    })
}

// ------------------------------------------------------------------ CIP ranks

/// `buildCIPInvariants`
fn build_cip_invariants(mol: &ROMol) -> Result<Vec<f64>> {
    const N_MASS_BITS: u32 = 10;
    const MAX_MASS: i32 = 1 << N_MASS_BITS;
    let mut res = Vec::with_capacity(mol.num_atoms() as usize);
    for atom in mol.atoms() {
        let num = (atom.atomic_num() % 128) as u64;
        let mut mass: i32 = 0;
        if atom.isotope() != 0 {
            mass = atom.isotope() as i32
                - PeriodicTable::get().get_most_common_isotope(atom.atomic_num())?;
            if mass >= 0 {
                mass += 1;
            }
        }
        mass += MAX_MASS / 2;
        if mass < 0 {
            mass = 0;
        } else {
            mass %= MAX_MASS;
        }
        let mut invariant: u64 = num;
        invariant = (invariant << N_MASS_BITS) | mass as u64;
        let mut mapnum = int_prop(atom, cp::molAtomMapNumber).unwrap_or(-1);
        mapnum = (mapnum + 1) % 1024;
        invariant = (invariant << 10) | (mapnum as i64 as u64);
        res.push(invariant as f64);
    }
    Ok(res)
}

#[derive(Clone, Copy)]
struct SortableCipRef {
    cip: usize,
    atom_idx: usize,
    curr_rank: i32,
}

/// `findSegmentsToResort`
fn find_segments_to_resort(
    sorted: &mut [SortableCipRef],
    entries: &[Vec<i32>],
    res: &mut Vec<(usize, usize)>,
) -> u32 {
    res.clear();
    let mut num_independent = sorted.len() as u32;
    let mut current = 0usize;
    let mut running_rank = 0;
    sorted[0].curr_rank = running_rank;
    let mut in_equal = false;
    for i in 1..sorted.len() {
        if entries[sorted[current].cip] == entries[sorted[i].cip] {
            sorted[i].curr_rank = running_rank;
            num_independent -= 1;
            if !in_equal {
                in_equal = true;
                res.push((i - 1, 0));
            }
        } else {
            running_rank += 1;
            sorted[i].curr_rank = running_rank;
            current = i;
            if in_equal {
                res.last_mut().unwrap().1 = i;
                in_equal = false;
            }
        }
    }
    if in_equal {
        res.last_mut().unwrap().1 = sorted.len() - 1;
    }
    num_independent
}

const K_MAX_BONDS: usize = 16;

struct PrecomputedBondFeatures {
    /// (count, neighbor index), strided by `K_MAX_BONDS` per atom
    counts_and_nbrs: Vec<(u8, i32)>,
    num_neighbors: Vec<u8>,
}

/// `computeBondFeatures`
fn compute_bond_features(mol: &ROMol) -> Result<PrecomputedBondFeatures> {
    let n = mol.num_atoms() as usize;
    let mut f = PrecomputedBondFeatures {
        counts_and_nbrs: vec![(0, 0); n * K_MAX_BONDS],
        num_neighbors: vec![0; n],
    };
    for atom_idx in 0..n {
        for (offset, b) in (atom_idx * K_MAX_BONDS..).zip(mol.atom_bond_indices(atom_idx as u32)) {
            let bond = mol.bond(b);
            let nbr = bond.other_atom_idx(atom_idx as u32)?;
            f.num_neighbors[nbr as usize] = f.num_neighbors[nbr as usize].wrapping_add(1);
            let slot = f
                .counts_and_nbrs
                .get_mut(offset)
                .ok_or_else(|| RDKitError::range("idx", offset, n * K_MAX_BONDS))?;
            slot.1 = nbr as i32;
            let mut is_p_special = false;
            if bond.bond_type() == BondType::Double {
                let na = mol.atom(nbr);
                if na.atomic_num() == 15 {
                    let d = mol.atom_degree(nbr);
                    is_p_special = d == 3 || d == 4;
                }
            }
            if is_p_special {
                slot.0 = slot.0.wrapping_add(1);
            } else {
                slot.0 = slot.0.wrapping_add(bond.bond_type().twice()?);
            }
        }
    }
    Ok(f)
}

fn recompute_ranks(sorted: &[SortableCipRef], ranks: &mut [u32]) {
    for e in sorted.iter().take(ranks.len()) {
        ranks[e.atom_idx] = e.curr_rank as u32;
    }
}

/// `iterateCIPRanks`
fn iterate_cip_ranks(
    mol: &ROMol,
    invars: &[f64],
    ranks: &mut [u32],
    seed_with_invars: bool,
) -> Result<()> {
    let num_atoms = mol.num_atoms() as usize;
    if num_atoms == 0 {
        return Ok(());
    }
    let mut entries: Vec<Vec<i32>> = vec![Vec::with_capacity(16); num_atoms];
    let mut sortable: Vec<SortableCipRef> = (0..num_atoms)
        .map(|i| SortableCipRef {
            cip: i,
            atom_idx: i,
            curr_rank: -1,
        })
        .collect();
    for i in 0..num_atoms {
        entries[i].push(invars[i] as i32);
    }
    sortable.sort_by(|a, b| entries[a.cip].cmp(&entries[b.cip]));
    let mut needs_sorting = Vec::new();
    let mut num_ranks = find_segments_to_resort(&mut sortable, &entries, &mut needs_sorting);
    recompute_ranks(&sortable, ranks);

    for i in 0..num_atoms {
        if seed_with_invars {
            entries[i][0] = invars[i] as i32;
        } else {
            entries[i][0] = mol.atom(i as u32).atomic_num() as i32;
            entries[i].push(ranks[i] as i32);
        }
    }
    let cip_rank_index = if seed_with_invars { 1 } else { 2 };

    let max_its = num_atoms / 2 + 1;
    let mut num_its = 0;
    let mut last_num_ranks: i64 = -1;
    let mut features = compute_bond_features(mol)?;

    while !needs_sorting.is_empty()
        && num_its < max_its
        && (last_num_ranks < 0 || (last_num_ranks as u32) < num_ranks)
    {
        for index in 0..num_atoms {
            let offset = K_MAX_BONDS * index;
            let num_nbrs = features.num_neighbors[index] as usize;
            let end = (offset + num_nbrs + 1).min(features.counts_and_nbrs.len());
            let slice = &mut features.counts_and_nbrs[offset..end];
            if num_nbrs > 1 {
                slice.sort_by(|a, b| ranks[b.1 as usize].cmp(&ranks[a.1 as usize]));
            }
            let entry = &mut entries[index];
            for &(count, idx) in slice.iter() {
                let v = ranks[idx as usize] as i32 + 1;
                entry.extend(std::iter::repeat_n(v, count as usize));
            }
            if !mol.atom(index as u32).has_query() {
                let nhs = mol.total_num_hs(index as u32, false)?;
                entry.extend(std::iter::repeat_n(0, nhs as usize));
            }
        }
        last_num_ranks = num_ranks as i64;
        for &(first, last) in &needs_sorting {
            sortable[first..=last].sort_by(|a, b| entries[a.cip].cmp(&entries[b.cip]));
        }
        num_ranks = find_segments_to_resort(&mut sortable, &entries, &mut needs_sorting);
        recompute_ranks(&sortable, ranks);
        if last_num_ranks as u32 != num_ranks {
            for i in 0..num_atoms {
                entries[i].resize(cip_rank_index + 1, 0);
                entries[i][cip_rank_index] = ranks[i] as i32;
            }
        }
        num_its += 1;
    }
    Ok(())
}

/// `Chirality::assignAtomCIPRanks`
pub fn assign_atom_cip_ranks(mol: &mut ROMol, ranks: &mut Vec<u32>) -> Result<()> {
    if ranks.is_empty() {
        ranks.resize(mol.num_atoms() as usize, 0);
    }
    let invars = build_cip_invariants(mol)?;
    iterate_cip_ranks(mol, &invars, ranks, false)?;
    for i in 0..mol.num_atoms() {
        mol.atom_mut(i)
            .set_prop(cp::_CIPRank, ranks[i as usize], true)?;
    }
    Ok(())
}

/// `Chirality::rerankAtoms`
fn rerank_atoms(mol: &mut ROMol, ranks: &mut [u32]) -> Result<()> {
    let mut factor: u32 = 100;
    while factor < mol.num_atoms() {
        factor *= 10;
    }
    let mut invars = vec![0f64; mol.num_atoms() as usize];
    for i in 0..mol.num_atoms() {
        let mut v = ranks[i as usize] as f64 * factor as f64;
        let atom = mol.atom(i);
        if let Some(RDValue::String(code)) = atom.get_prop(cp::_CIPCode) {
            if code == "S" {
                v += 10.0;
            } else if code == "R" {
                v += 20.0;
            }
        }
        for b in mol.atom_bonds(i) {
            if b.bond_type() == BondType::Double {
                if b.stereo() == BondStereo::StereoE {
                    v += 1.0;
                } else if b.stereo() == BondStereo::StereoZ {
                    v += 2.0;
                }
            }
        }
        invars[i as usize] = v;
    }
    iterate_cip_ranks(mol, &invars, ranks, true)?;
    for i in 0..mol.num_atoms() {
        mol.atom_mut(i)
            .set_prop(cp::_CIPRank, ranks[i as usize], false)?;
    }
    Ok(())
}

// ------------------------------------------------------------------ atoms

/// `Chirality::isAtomPotentialChiralCenter`. Returns (legalCenter, hasDupes)
/// and fills `nbrs` with (rank, bond index) pairs when ranks are available.
pub fn is_atom_potential_chiral_center(
    mol: &ROMol,
    atom: u32,
    ranks: &[u32],
    nbrs: &mut Vec<IntPair>,
) -> Result<(bool, bool)> {
    let mut legal = true;
    let mut has_dupes = false;
    let a = mol.atom(atom);
    let nz_degree = get_atom_nonzero_degree(mol, atom);
    let total_hs = mol.total_num_hs(atom, false)?;
    let tnz_degree = nz_degree + total_hs;
    if tnz_degree > 4 {
        legal = false;
    } else {
        if tnz_degree < 3 {
            legal = false;
        } else if nz_degree < 3 && a.atomic_num() != 15 && a.atomic_num() != 33 {
            legal = false;
        } else if nz_degree == 3 {
            if total_hs == 1 {
                if has_protium_neighbor(mol, atom) {
                    legal = false;
                }
            } else {
                legal = false;
                match a.atomic_num() {
                    7 => {
                        if a.hybridization() == HybridizationType::SP3
                            && !atom_has_conjugated_bond(mol, atom)
                            && (mol
                                .ring_info()
                                .is_atom_in_ring_of_size(atom, 3)
                                .unwrap_or(false)
                                || query_is_atom_bridgehead(mol, atom))
                        {
                            legal = true;
                        }
                    }
                    15 | 33 => legal = true,
                    16 | 34 => {
                        let ev = a.valence(ValenceType::Explicit)?;
                        if ev == 4 || (ev == 3 && a.formal_charge() == 1) {
                            legal = true;
                        }
                    }
                    _ => {}
                }
            }
        }
        if legal && !ranks.is_empty() {
            let mut codes_seen = vec![false; mol.num_atoms() as usize];
            for b in mol.atom_bond_indices(atom) {
                let other = mol.bond(b).other_atom_idx(atom)?;
                nbrs.push((ranks[other as usize] as i32, b as i32));
                if !bond_affects_atom_chirality(mol, b, atom) {
                    continue;
                }
                let r = ranks[other as usize];
                if r >= mol.num_atoms() {
                    return Err(RDKitError::invariant(
                        "CIP rank higher than the number of atoms.",
                    ));
                }
                if codes_seen[r as usize] {
                    has_dupes = true;
                    break;
                }
                codes_seen[r as usize] = true;
            }
        }
    }
    Ok((legal, has_dupes))
}

/// `Chirality::assignAtomChiralCodes`: returns (unassigned atoms remain,
/// any atom was assigned).
fn assign_atom_chiral_codes(
    mol: &mut ROMol,
    ranks: &mut Vec<u32>,
    flag_possible: bool,
) -> Result<(bool, bool)> {
    let mut atom_changed = false;
    let mut unassigned: u32 = 0;
    for i in 0..mol.num_atoms() {
        let mut tag = mol.atom(i).chiral_tag();
        if !(flag_possible || (tag != ChiralType::Unspecified && tag != ChiralType::Other)) {
            continue;
        }
        if mol.atom(i).has_prop(cp::_CIPCode) {
            continue;
        }
        if ranks.is_empty() {
            assign_atom_cip_ranks(mol, ranks)?;
        }
        let mut nbrs = Vec::new();
        let (legal, has_dupes) = is_atom_potential_chiral_center(mol, i, ranks, &mut nbrs)?;
        if legal {
            unassigned += 1;
        }
        if legal && !has_dupes && flag_possible {
            mol.atom_mut(i)
                .set_prop(cp::_ChiralityPossible, 1i32, false)?;
        }
        if legal && !has_dupes && tag != ChiralType::Unspecified && tag != ChiralType::Other {
            atom_changed = true;
            unassigned -= 1;
            nbrs.sort_by_key(|p| p.0);
            let nbr_indices: Vec<i32> = nbrs.iter().map(|p| p.1).collect();
            let mut n_swaps = mol.perturbation_order(i, &nbr_indices)?;
            if nbr_indices.len() == 3 && mol.total_num_hs(i, false)? == 1 {
                n_swaps += 1;
            }
            if n_swaps % 2 != 0 {
                tag = if tag == ChiralType::TetrahedralCCW {
                    ChiralType::TetrahedralCW
                } else {
                    ChiralType::TetrahedralCCW
                };
            }
            let code = if tag == ChiralType::TetrahedralCCW {
                "S"
            } else {
                "R"
            };
            mol.atom_mut(i)
                .set_prop(cp::_CIPCode, code.to_string(), false)?;
        }
    }
    Ok((unassigned > 0, atom_changed))
}

// ------------------------------------------------------------------ bonds

/// `Chirality::findAtomNeighborDirHelper`
fn find_atom_neighbor_dir_helper(
    mol: &ROMol,
    atom: u32,
    ref_bond: u32,
    ranks: &[u32],
    neighbors: &mut Vec<(i32, BondDir)>,
    has_explicit_unknown_stereo: &mut bool,
) -> Result<()> {
    let mut seen_dir = false;
    for b in mol.atom_bond_indices(atom) {
        let bond = mol.bond(b);
        if !*has_explicit_unknown_stereo
            && (bond.bond_dir() == BondDir::Unknown
                || int_prop(bond, cp::_UnknownStereo).is_some_and(|v| v != 0))
        {
            *has_explicit_unknown_stereo = true;
        }
        let mut dir = bond.bond_dir();
        if b != ref_bond {
            if matches!(dir, BondDir::EndDownRight | BondDir::EndUpRight) {
                seen_dir = true;
                if atom != bond.begin_atom_idx() {
                    dir = if dir == BondDir::EndDownRight {
                        BondDir::EndUpRight
                    } else {
                        BondDir::EndDownRight
                    };
                }
            }
            neighbors.push((bond.other_atom_idx(atom)? as i32, dir));
        }
    }
    if !seen_dir {
        neighbors.clear();
    } else if neighbors.len() == 2
        && ranks[neighbors[0].0 as usize] == ranks[neighbors[1].0 as usize]
    {
        neighbors.clear();
    } else if !matches!(neighbors[0].1, BondDir::EndDownRight | BondDir::EndUpRight) {
        if neighbors.len() <= 1 {
            return Err(RDKitError::invariant("too few neighbors"));
        }
        neighbors[0].1 = if neighbors[1].1 == BondDir::EndDownRight {
            BondDir::EndUpRight
        } else {
            BondDir::EndDownRight
        };
    } else if neighbors.len() > 1
        && !matches!(neighbors[1].1, BondDir::EndDownRight | BondDir::EndUpRight)
    {
        neighbors[1].1 = if neighbors[0].1 == BondDir::EndDownRight {
            BondDir::EndUpRight
        } else {
            BondDir::EndDownRight
        };
    }
    Ok(())
}

/// `Chirality::assignBondStereoCodes`: returns (unassigned bonds remain,
/// any bond was assigned).
fn assign_bond_stereo_codes(mol: &mut ROMol, ranks: &mut Vec<u32>) -> Result<(bool, bool)> {
    let mut assigned = false;
    let mut unassigned: u32 = 0;
    let mut bonds_to_clear = vec![false; mol.num_bonds() as usize];
    for bi in 0..mol.num_bonds() {
        if mol.bond(bi).bond_type() != BondType::Double {
            continue;
        }
        if mol.bond(bi).stereo() != BondStereo::StereoNone {
            continue;
        }
        if ranks.is_empty() {
            assign_atom_cip_ranks(mol, ranks)?;
        }
        mol.bond_mut(bi).stereo_atoms_mut().clear();
        if !should_detect_double_bond_stereo(mol, bi) {
            continue;
        }
        let beg = mol.bond(bi).begin_atom_idx();
        let end = mol.bond(bi).end_atom_idx();
        let bd = mol.atom_degree(beg);
        let ed = mol.atom_degree(end);
        if !((bd == 2 || bd == 3) && (ed == 2 || ed == 3)) {
            continue;
        }
        unassigned += 1;
        let mut beg_nbrs = Vec::new();
        let mut end_nbrs = Vec::new();
        let mut has_unknown = int_prop(mol.atom(beg), cp::_UnknownStereo).is_some_and(|v| v != 0)
            || int_prop(mol.atom(end), cp::_UnknownStereo).is_some_and(|v| v != 0);
        find_atom_neighbor_dir_helper(mol, beg, bi, ranks, &mut beg_nbrs, &mut has_unknown)?;
        find_atom_neighbor_dir_helper(mol, end, bi, ranks, &mut end_nbrs, &mut has_unknown)?;
        if beg_nbrs.is_empty() || end_nbrs.is_empty() {
            continue;
        }
        let pick = |n: &[(i32, BondDir)]| {
            if n.len() == 1 || ranks[n[0].0 as usize] > ranks[n[1].0 as usize] {
                n[0]
            } else {
                n[1]
            }
        };
        let (beg_nbr_aid, beg_dir) = pick(&beg_nbrs);
        let (end_nbr_aid, end_dir) = pick(&end_nbrs);
        let conflicting_begin = beg_nbrs.len() == 2 && beg_nbrs[0].1 == beg_nbrs[1].1;
        let conflicting_end = end_nbrs.len() == 2 && end_nbrs[0].1 == end_nbrs[1].1;
        if conflicting_begin || conflicting_end {
            mol.bond_mut(bi).set_stereo(BondStereo::StereoNone)?;
            crate::rd_log_warning!(
                "Conflicting single bond directions around double bond at index {}.",
                bi
            );
            crate::rd_log_warning!(
                "  BondStereo set to STEREONONE and single bond directions set to NONE."
            );
            assigned = true;
            let mut clear = |nbrs: &[(i32, BondDir)], center: u32| {
                for n in &nbrs[..2] {
                    if let Some(b) = mol.bond_idx_between(n.0 as u32, center) {
                        bonds_to_clear[b as usize] = true;
                    }
                }
            };
            if conflicting_begin {
                clear(&beg_nbrs, beg);
            }
            if conflicting_end {
                clear(&end_nbrs, end);
            }
        } else {
            let bond = mol.bond_mut(bi);
            bond.stereo_atoms_mut().push(beg_nbr_aid);
            bond.stereo_atoms_mut().push(end_nbr_aid);
            if has_unknown {
                bond.set_stereo(BondStereo::StereoAny)?;
            } else if beg_dir == end_dir {
                bond.set_stereo(BondStereo::StereoZ)?;
            } else {
                bond.set_stereo(BondStereo::StereoE)?;
            }
            assigned = true;
        }
        unassigned -= 1;
    }
    for (i, &c) in bonds_to_clear.iter().enumerate() {
        if c {
            mol.bond_mut(i as u32).set_bond_dir(BondDir::None);
        }
    }
    Ok((unassigned > 0, assigned))
}

// ------------------------------------------------------------------ ring stereo

/// `Chirality::atomIsCandidateForRingStereochem`
fn atom_is_candidate_for_ring_stereochem(
    mol: &mut ROMol,
    atom: u32,
    atom_ranks: &[u32],
) -> Result<bool> {
    if let Some(v) = mol
        .atom(atom)
        .get_prop(cp::_ringStereochemCand)
        .and_then(RDValue::as_bool)
    {
        return Ok(v);
    }
    let mut res = false;
    let ri = mol.ring_info();
    if ri.is_initialized() && ri.num_atom_rings(atom).unwrap_or(0) > 0 {
        let a = mol.atom(atom);
        if a.atomic_num() == 7
            && mol.total_degree(atom)? == 3
            && !ri.is_atom_in_ring_of_size(atom, 3).unwrap_or(false)
            && !query_is_atom_bridgehead(mol, atom)
        {
            // note: upstream returns without caching the result here
            return Ok(false);
        }
        let mut non_ring_nbrs = Vec::new();
        let mut ring_nbrs = 0usize;
        let mut ring_nbr_ranks = BTreeSet::new();
        for b in mol.atom_bond_indices(atom) {
            let nbr = mol.bond(b).other_atom_idx(atom)?;
            if ri.num_bond_rings(b).unwrap_or(0) == 0 {
                non_ring_nbrs.push(nbr);
            } else {
                ring_nbrs += 1;
                ring_nbr_ranks.insert(atom_ranks[nbr as usize]);
            }
        }
        res = match non_ring_nbrs.len() {
            2 => {
                atom_ranks[non_ring_nbrs[0] as usize] != atom_ranks[non_ring_nbrs[1] as usize]
                    && ring_nbrs != ring_nbr_ranks.len()
            }
            1 => ring_nbrs > ring_nbr_ranks.len(),
            0 => {
                (ring_nbrs == 4 && ring_nbr_ranks.len() == 3)
                    || (ring_nbrs == 3 && ring_nbr_ranks.len() == 2)
            }
            _ => false,
        };
    }
    mol.atom_mut(atom)
        .set_prop(cp::_ringStereochemCand, res, true)?;
    Ok(res)
}

fn entry_idx(e: i32) -> usize {
    if e < 0 {
        (-e - 1) as usize
    } else {
        (e - 1) as usize
    }
}

/// `Chirality::findChiralAtomSpecialCases`
fn find_chiral_atom_special_cases(
    mol: &mut ROMol,
    possible: &mut [bool],
    atom_ranks: &[u32],
) -> Result<()> {
    possible.iter_mut().for_each(|x| *x = false);
    if !mol.ring_info().is_symm_sssr() {
        symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    }
    let n = mol.num_atoms() as usize;
    let mut atoms_seen = vec![false; n];
    let mut atoms_used = vec![false; n];
    let mut bonds_seen = vec![false; mol.num_bonds() as usize];

    for atom in 0..mol.num_atoms() {
        if atoms_seen[atom as usize] {
            continue;
        }
        if mol.atom(atom).chiral_tag() == ChiralType::Unspecified
            || mol.atom(atom).has_prop(cp::_CIPCode)
            || mol.ring_info().num_atom_rings(atom).unwrap_or(0) == 0
            || !atom_is_candidate_for_ring_stereochem(mol, atom, atom_ranks)?
        {
            continue;
        }
        let mut next_atoms: VecDeque<u32> = VecDeque::new();
        let bonds: Vec<u32> = mol.atom_bond_indices(atom).collect();
        for b in bonds {
            if !bonds_seen[b as usize] {
                bonds_seen[b as usize] = true;
                if mol.ring_info().num_bond_rings(b).unwrap_or(0) > 0 {
                    let o = mol.bond(b).other_atom_idx(atom)?;
                    if !atoms_seen[o as usize] {
                        next_atoms.push_back(o);
                        atoms_used[o as usize] = true;
                    }
                }
            }
        }
        let mut ring_stereo_atoms: Vec<i32> = Vec::new();
        if !next_atoms.is_empty() {
            if let Some(v) = vec_int_prop(mol.atom(atom), cp::_ringStereoAtoms) {
                ring_stereo_atoms = v;
            }
        }
        while let Some(ratom) = next_atoms.pop_front() {
            atoms_seen[ratom as usize] = true;
            if mol.atom(ratom).chiral_tag() != ChiralType::Unspecified
                && !mol.atom(ratom).has_prop(cp::_CIPCode)
                && atom_is_candidate_for_ring_stereochem(mol, ratom, atom_ranks)?
            {
                let same = if mol.atom(ratom).chiral_tag() == mol.atom(atom).chiral_tag() {
                    1
                } else {
                    -1
                };
                ring_stereo_atoms.push(same * (ratom as i32 + 1));
                let mut oring =
                    vec_int_prop(mol.atom(ratom), cp::_ringStereoAtoms).unwrap_or_default();
                oring.push(same * (atom as i32 + 1));
                mol.atom_mut(ratom)
                    .set_prop(cp::_ringStereoAtoms, oring, true)?;
                possible[ratom as usize] = true;
                possible[atom as usize] = true;
            }
            let rbonds: Vec<u32> = mol.atom_bond_indices(ratom).collect();
            for b in rbonds {
                if !bonds_seen[b as usize] {
                    bonds_seen[b as usize] = true;
                    if mol.ring_info().num_bond_rings(b).unwrap_or(0) > 0 {
                        let o = mol.bond(b).other_atom_idx(ratom)?;
                        if !atoms_seen[o as usize] && !atoms_used[o as usize] {
                            next_atoms.push_back(o);
                            atoms_used[o as usize] = true;
                        }
                    }
                }
            }
        }
        if !ring_stereo_atoms.is_empty() {
            mol.atom_mut(atom)
                .set_prop(cp::_ringStereoAtoms, ring_stereo_atoms.clone(), true)?;
            for (ri, &entry) in ring_stereo_atoms.iter().enumerate() {
                let ring_atom_idx = entry_idx(entry);
                let mut lring = vec_int_prop(mol.atom(ring_atom_idx as u32), cp::_ringStereoAtoms)
                    .unwrap_or_default();
                if lring.is_empty() {
                    return Err(RDKitError::invariant("no other ring atoms found."));
                }
                for &oentry in &ring_stereo_atoms[ri + 1..] {
                    let oidx = entry_idx(oentry);
                    let different = (entry < 0) ^ (oentry < 0);
                    lring.push(if different {
                        -(oidx as i32 + 1)
                    } else {
                        oidx as i32 + 1
                    });
                    let mut olring = vec_int_prop(mol.atom(oidx as u32), cp::_ringStereoAtoms)
                        .unwrap_or_default();
                    if olring.is_empty() {
                        return Err(RDKitError::invariant("no other ring atoms found."));
                    }
                    olring.push(if different {
                        -(ring_atom_idx as i32 + 1)
                    } else {
                        ring_atom_idx as i32 + 1
                    });
                    mol.atom_mut(oidx as u32)
                        .set_prop(cp::_ringStereoAtoms, olring, false)?;
                }
                mol.atom_mut(ring_atom_idx as u32)
                    .set_prop(cp::_ringStereoAtoms, lring, false)?;
            }
        } else {
            possible[atom as usize] = false;
        }
        atoms_seen[atom as usize] = true;
    }
    Ok(())
}

// ------------------------------------------------------------------ atropisomers

fn is_atrop(s: BondStereo) -> bool {
    matches!(s, BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW)
}

/// `Atropisomers::detectAtropisomerChirality(mol, conf, replaceExistingTags)`
///
/// Only the candidate-selection part is ported so far; it returns early for
/// every molecule without wedged bonds (i.e. everything SMILES can express).
fn detect_atropisomer_chirality(mol: &mut ROMol, replace_existing_tags: bool) -> Result<()> {
    let has_conf = mol.num_conformers() > 0;
    let is_3d = has_conf && mol.conformers()[0].is_3d();
    let mut bonds_to_try = BTreeSet::new();
    if replace_existing_tags && is_3d {
        for b in mol.bonds() {
            if is_atrop(b.stereo()) {
                bonds_to_try.insert(b.idx());
            }
        }
    }
    for b in mol.bonds() {
        if can_have_direction(b.bond_type())
            && matches!(b.bond_dir(), BondDir::BeginDash | BondDir::BeginWedge)
        {
            for nb in mol.atom_bond_indices(b.begin_atom_idx()) {
                if nb == b.idx() {
                    continue;
                }
                if !replace_existing_tags && is_atrop(mol.bond(nb).stereo()) {
                    continue;
                }
                bonds_to_try.insert(nb);
            }
        }
    }
    if replace_existing_tags {
        for &b in &bonds_to_try {
            if is_atrop(mol.bond(b).stereo()) {
                mol.bond_mut(b).set_stereo(BondStereo::StereoNone)?;
            }
        }
    }
    if bonds_to_try.is_empty() {
        return Ok(());
    }
    let mut any = false;
    for &b in &bonds_to_try {
        let (bg, en) = (mol.bond(b).begin_atom_idx(), mol.bond(b).end_atom_idx());
        for a in [bg, en] {
            if mol.atom(a).needs_update_property_cache() {
                mol.update_atom_property_cache(a, false)?;
            }
        }
        let bond = mol.bond(b);
        let (bd, ed) = (mol.total_degree(bg)?, mol.total_degree(en)?);
        if bond.bond_type() == BondType::Single
            && bond.stereo() != BondStereo::StereoAny
            && (2..=3).contains(&bd)
            && (2..=3).contains(&ed)
        {
            any = true;
            break;
        }
    }
    if !any {
        return Ok(());
    }
    Err(RDKitError::invariant(
        "atropisomer perception from wedged bonds is not ported yet",
    ))
}

/// `canHaveDirection` (Atropisomers.cpp)
fn can_have_direction(t: BondType) -> bool {
    matches!(t, BondType::Single | BondType::Aromatic)
}

/// `Atropisomers::cleanupAtropisomerStereoGroups`
fn cleanup_atropisomer_stereo_groups(mol: &mut ROMol) -> Result<()> {
    let mut new = Vec::new();
    for sg in mol.stereo_groups() {
        let mut ok_atoms = Vec::new();
        let mut ok_bonds: Vec<u32> = Vec::new();
        for &a in sg.atoms() {
            let mut found = false;
            for b in mol.atom_bond_indices(a) {
                if is_atrop(mol.bond(b).stereo()) {
                    found = true;
                    if !ok_bonds.contains(&b) {
                        ok_bonds.push(b);
                    }
                }
            }
            if !found {
                ok_atoms.push(a);
            }
        }
        if ok_bonds.is_empty() {
            new.push(sg.clone());
        } else {
            new.push(crate::graphmol::stereo_group::StereoGroup::new(
                sg.group_type(),
                ok_atoms,
                ok_bonds,
                0,
            )?);
        }
    }
    mol.set_stereo_groups(new)
}

/// `detectAtropisomersPostSanitization`
fn detect_atropisomers_post_sanitization(mol: &mut ROMol, clean_it: bool) -> Result<()> {
    detect_atropisomer_chirality(mol, clean_it)?;
    if !clean_it || !mol.ring_info().is_sssr_or_better() {
        return Ok(());
    }
    let mut removed = false;
    for b in 0..mol.num_bonds() {
        if is_atrop(mol.bond(b).stereo())
            && mol.ring_info().num_bond_rings(b).unwrap_or(0) > 0
            && mol.ring_info().min_bond_ring_size(b).unwrap_or(0) < 9
        {
            mol.bond_mut(b).set_stereo(BondStereo::StereoNone)?;
            removed = true;
        }
    }
    if removed {
        cleanup_atropisomer_stereo_groups(mol)?;
    }
    Ok(())
}

// ------------------------------------------------------------------ driver

/// `Chirality::legacyStereoPerception`
pub fn legacy_stereo_perception(
    mol: &mut ROMol,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
) -> Result<()> {
    mol.clear_prop("_needsDetectBondStereo");
    if !mol.ring_info().is_find_fast_or_better() {
        fast_find_rings(mol)?;
    }

    let mut has_stereo_atoms = false;
    let mut has_potential_stereo_atoms = false;
    for i in 0..mol.num_atoms() {
        if clean_it {
            let a = mol.atom_mut(i);
            a.clear_prop(cp::_CIPCode);
            a.clear_prop(cp::_ChiralityPossible);
        }
        let tag = mol.atom(i).chiral_tag();
        if !has_stereo_atoms && tag != ChiralType::Unspecified && tag != ChiralType::Other {
            has_stereo_atoms = true;
        } else if !has_potential_stereo_atoms {
            let mut nbrs = Vec::new();
            has_potential_stereo_atoms = is_atom_potential_chiral_center(mol, i, &[], &mut nbrs)?.0;
        }
    }

    let mut has_stereo_bonds = false;
    let mut has_potential_stereo_bonds = false;
    for bi in 0..mol.num_bonds() {
        if clean_it {
            mol.bond_mut(bi).clear_prop(cp::_CIPCode);
            let bt = mol.bond(bi).bond_type();
            if matches!(bt, BondType::Double | BondType::Aromatic)
                && !should_detect_double_bond_stereo(mol, bi)
            {
                let b = mol.bond_mut(bi);
                if b.bond_dir() == BondDir::EitherDouble {
                    b.set_bond_dir(BondDir::None);
                }
                if b.stereo() != BondStereo::StereoNone {
                    b.set_stereo(BondStereo::StereoNone)?;
                    b.stereo_atoms_mut().clear();
                }
                continue;
            } else if bt == BondType::Double {
                let b = mol.bond_mut(bi);
                if b.bond_dir() == BondDir::EitherDouble {
                    b.set_stereo(BondStereo::StereoAny)?;
                    b.stereo_atoms_mut().clear();
                    b.set_bond_dir(BondDir::None);
                } else if b.stereo() != BondStereo::StereoAny {
                    b.set_stereo(BondStereo::StereoNone)?;
                    b.stereo_atoms_mut().clear();
                }
            }
        }
        if !has_stereo_bonds && mol.bond(bi).bond_type() == BondType::Double {
            let is_dir = |d: BondDir| matches!(d, BondDir::EndDownRight | BondDir::EndUpRight);
            let (bg, en) = (mol.bond(bi).begin_atom_idx(), mol.bond(bi).end_atom_idx());
            let mut is_specified = false;
            if mol.atom_bonds(bg).any(|b| is_dir(b.bond_dir()))
                || mol.atom_bonds(en).any(|b| is_dir(b.bond_dir()))
            {
                has_stereo_bonds = true;
                is_specified = true;
            }
            if !has_potential_stereo_bonds
                && !is_specified
                && should_detect_double_bond_stereo(mol, bi)
            {
                has_potential_stereo_bonds = true;
            }
        }
        if !clean_it && has_stereo_bonds && has_potential_stereo_bonds {
            break;
        }
    }

    let mut atom_ranks: Vec<u32> = Vec::new();
    let mut keep_going = has_stereo_atoms | has_stereo_bonds;
    if !keep_going {
        keep_going = flag_possible_stereo_centers
            && (has_potential_stereo_atoms || has_potential_stereo_bonds);
    }
    while keep_going {
        let changed_atoms;
        let changed_bonds;
        if has_stereo_atoms || has_potential_stereo_atoms {
            (has_stereo_atoms, changed_atoms) =
                assign_atom_chiral_codes(mol, &mut atom_ranks, flag_possible_stereo_centers)?;
        } else {
            changed_atoms = false;
        }
        if has_stereo_bonds || has_potential_stereo_bonds {
            (has_stereo_bonds, changed_bonds) = assign_bond_stereo_codes(mol, &mut atom_ranks)?;
        } else {
            changed_bonds = false;
        }
        keep_going = (has_stereo_atoms || has_stereo_bonds) && (changed_atoms || changed_bonds);
        if keep_going {
            rerank_atoms(mol, &mut atom_ranks)?;
        }
    }

    if clean_it {
        for i in 0..mol.num_atoms() {
            let a = mol.atom_mut(i);
            a.clear_prop(cp::_ringStereochemCand);
            a.clear_prop(cp::_ringStereoAtoms);
        }
        let mut possible = vec![false; mol.num_atoms() as usize];
        find_chiral_atom_special_cases(mol, &mut possible, &atom_ranks)?;
        for i in 0..mol.num_atoms() {
            let a = mol.atom(i);
            if a.chiral_tag() != ChiralType::Unspecified
                && !has_non_tetrahedral_stereo(mol, i)
                && !a.has_prop(cp::_CIPCode)
                && (!possible[i as usize] || !a.has_prop(cp::_ringStereoAtoms))
            {
                mol.atom_mut(i).set_chiral_tag(ChiralType::Unspecified);
                let a = mol.atom(i);
                if a.num_explicit_hs() == 1 && a.formal_charge() == 0 && !a.is_aromatic() {
                    let a = mol.atom_mut(i);
                    a.set_num_explicit_hs(0);
                    a.set_no_implicit(false);
                    mol.calc_explicit_valence(i, false)?;
                    mol.calc_implicit_valence(i, false)?;
                }
            }
        }
    }
    detect_atropisomers_post_sanitization(mol, clean_it)?;
    if clean_it {
        let mut found_atropisomer = false;
        for bi in 0..mol.num_bonds() {
            let (bg, en) = (mol.bond(bi).begin_atom_idx(), mol.bond(bi).end_atom_idx());
            if matches!(
                mol.bond(bi).bond_dir(),
                BondDir::BeginWedge | BondDir::BeginDash
            ) && mol.atom(bg).chiral_tag() == ChiralType::Unspecified
                && mol.atom(en).chiral_tag() == ChiralType::Unspecified
            {
                let has_atrop = mol.atom_bonds(bg).any(|b| is_atrop(b.stereo()));
                if has_atrop {
                    found_atropisomer = true;
                } else {
                    mol.bond_mut(bi).set_bond_dir(BondDir::None);
                }
            }

            let b = mol.bond(bi);
            if b.bond_type() == BondType::Double
                && (b.bond_dir() == BondDir::EitherDouble || b.stereo() == BondStereo::StereoAny)
                && (mol.atom_degree(bg) == 1 || mol.atom_degree(en) == 1)
            {
                let b = mol.bond_mut(bi);
                if b.bond_dir() == BondDir::EitherDouble {
                    b.set_bond_dir(BondDir::None);
                }
                b.set_stereo(BondStereo::StereoNone)?;
            }

            let b = mol.bond(bi);
            if b.bond_type() == BondType::Double
                && matches!(b.stereo(), BondStereo::StereoAny | BondStereo::StereoNone)
            {
                for batom in [bg, en] {
                    let nbr_bonds: Vec<u32> = mol.atom_bond_indices(batom).collect();
                    for nbi in nbr_bonds {
                        if nbi == bi {
                            continue;
                        }
                        let nb = mol.bond(nbi);
                        if matches!(nb.bond_dir(), BondDir::EndDownRight | BondDir::EndUpRight)
                            && matches!(nb.bond_type(), BondType::Single | BondType::Aromatic)
                        {
                            let other = nb.other_atom_idx(batom)?;
                            let ok_to_clear = !mol.atom_bonds(other).any(|j| {
                                j.bond_type() == BondType::Double
                                    && !matches!(
                                        j.stereo(),
                                        BondStereo::StereoAny | BondStereo::StereoNone
                                    )
                            });
                            if ok_to_clear {
                                mol.bond_mut(nbi).set_bond_dir(BondDir::None);
                            }
                        }
                    }
                }
            }
        }
        if found_atropisomer || mol.bonds().iter().any(|b| is_atrop(b.stereo())) {
            cleanup_atropisomer_stereo_groups(mol)?;
        }
        crate::graphmol::mol_ops::cleanup_stereo_groups(mol)?;
    }
    Ok(())
}

/// `MolOps::assignStereochemistry(mol, cleanIt, force, flagPossibleStereoCenters)`
pub fn assign_stereochemistry(
    mol: &mut ROMol,
    clean_it: bool,
    force: bool,
    flag_possible_stereo_centers: bool,
) -> Result<()> {
    if !force && mol.has_prop(cp::_StereochemDone) {
        return Ok(());
    }
    if mol.needs_update_property_cache() {
        mol.update_property_cache(false)?;
    }
    if !super::get_use_legacy_stereo_perception() {
        return Err(RDKitError::invariant(
            "the new stereo perception algorithm is not ported yet; \
             use legacy stereo perception",
        ));
    }
    legacy_stereo_perception(mol, clean_it, flag_possible_stereo_centers)?;
    mol.set_prop(cp::_StereochemDone, 1i32, true)?;
    Ok(())
}

/// `Chirality::detail::isAtomPotentialTetrahedralCenter` (FindStereo.cpp)
pub fn is_atom_potential_tetrahedral_center(mol: &ROMol, atom: u32) -> Result<bool> {
    let a = mol.atom(atom);
    let nz_degree = get_atom_nonzero_degree(mol, atom);
    let total_hs = mol.total_num_hs(atom, false)?;
    if nz_degree + total_hs > 4 {
        return Ok(false);
    }
    if nz_degree == 4 {
        return Ok(true);
    }
    if nz_degree <= 1 {
        return Ok(false);
    }
    let anum = a.atomic_num();
    if nz_degree < 3 && anum != 15 && anum != 33 {
        return Ok(false);
    }
    if anum == 15 || anum == 33 {
        return Ok(true);
    }
    if nz_degree == 3 {
        if total_hs == 1 {
            return Ok(!has_protium_neighbor(mol, atom));
        }
        let ev = a.valence(ValenceType::Explicit)?;
        if (anum == 16 || anum == 34) && (ev == 4 || (ev == 3 && a.formal_charge() == 1)) {
            return Ok(true);
        }
        if anum == 7
            && a.hybridization() == HybridizationType::SP3
            && !atom_has_conjugated_bond(mol, atom)
            && (mol
                .ring_info()
                .is_atom_in_ring_of_size(atom, 3)
                .unwrap_or(false)
                || query_is_atom_bridgehead(mol, atom))
        {
            return Ok(true);
        }
        return Ok(false);
    }
    Ok(false)
}
