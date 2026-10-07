//! Ported from RDKit `Code/GraphMol/Aromaticity.cpp` (aromaticity models and
//! the `RingUtils` fused-ring helpers).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, BTreeSet, HashSet};

use super::atom::ValenceType;
use super::bond::BondType;
use super::find_rings::{SymmetrizeSSSRAlgorithm, convert_to_bonds, symmetrize_sssr};
use super::periodic_table::PeriodicTable;
use super::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::utils::{intersect, next_combination};
use crate::{RDKitError, Result};

const MAX_FUSED_AROMATIC_RING_SIZE: usize = 24;

/// `INT_INT_VECT_MAP` used for ring-neighbor maps.
pub type RingNeighborMap = BTreeMap<i32, Vec<i32>>;

// ------------------------------------------------------------------ RingUtils

/// `RingUtils::pickFusedRings`: depth-first collection of the fused system
/// containing ring `curr`.
pub fn pick_fused_rings(
    curr: i32,
    neigh_map: &RingNeighborMap,
    res: &mut Vec<i32>,
    done: &mut [bool],
) {
    done[curr as usize] = true;
    res.push(curr);
    let mut stack: Vec<(i32, usize)> = vec![(curr, 0)];
    while let Some(&mut (ring, ref mut next)) = stack.last_mut() {
        let neighs = &neigh_map[&ring];
        if *next == neighs.len() {
            stack.pop();
            continue;
        }
        let neigh = neighs[*next];
        *next += 1;
        if done[neigh as usize] {
            continue;
        }
        done[neigh as usize] = true;
        res.push(neigh);
        stack.push((neigh, 0));
    }
}

/// `RingUtils::checkFused`
pub fn check_fused(rids: &[i32], ring_neighs: &RingNeighborMap) -> bool {
    let mut done = vec![false; ring_neighs.len()];
    for &rid in ring_neighs.keys() {
        if !rids.contains(&rid) {
            done[rid as usize] = true;
        }
    }
    let mut fused = Vec::new();
    pick_fused_rings(rids[0], ring_neighs, &mut fused, &mut done);
    fused.len() == rids.len()
}

/// `RingUtils::makeRingNeighborMap`
pub fn make_ring_neighbor_map(
    brings: &[Vec<i32>],
    max_size: usize,
    max_overlap_size: usize,
) -> RingNeighborMap {
    let mut map = RingNeighborMap::new();
    let n = brings.len();
    for i in 0..n {
        map.entry(i as i32).or_default();
        if max_size != 0 && brings[i].len() > max_size {
            continue;
        }
        for j in i + 1..n {
            if max_size != 0 && brings[j].len() > max_size {
                continue;
            }
            let inter = intersect(&brings[i], &brings[j]);
            if !inter.is_empty() && (max_overlap_size == 0 || inter.len() <= max_overlap_size) {
                map.entry(i as i32).or_default().push(j as i32);
                map.entry(j as i32).or_default().push(i as i32);
            }
        }
    }
    map
}

// ------------------------------------------------------------------ electron donors

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)] // OneOrTwo exists upstream but is never assigned
enum ElectronDonorType {
    Vacant,
    One,
    Two,
    OneOrTwo,
    Any,
    No,
}

fn min_max_atom_elecs(d: ElectronDonorType) -> (i32, i32) {
    match d {
        ElectronDonorType::Any | ElectronDonorType::OneOrTwo => (1, 2),
        ElectronDonorType::One => (1, 1),
        ElectronDonorType::Two => (2, 2),
        ElectronDonorType::No | ElectronDonorType::Vacant => (0, 0),
    }
}

fn valence_contrib(mol: &ROMol, bond: u32, atom: u32) -> f64 {
    mol.bond(bond).valence_contrib(mol, atom)
}

fn incident_non_cyclic_multiple_bond(mol: &ROMol, at: u32) -> Option<u32> {
    for b in mol.atom_bond_indices(at) {
        if mol.ring_info().num_bond_rings(b).unwrap_or(0) == 0 && valence_contrib(mol, b, at) >= 2.0
        {
            return Some(mol.bond(b).other_atom_idx(at).expect("in bond"));
        }
    }
    None
}

fn incident_cyclic_multiple_bond(mol: &ROMol, at: u32) -> bool {
    mol.atom_bond_indices(at).any(|b| {
        mol.ring_info().num_bond_rings(b).unwrap_or(0) != 0 && valence_contrib(mol, b, at) >= 2.0
    })
}

fn incident_multiple_bond(mol: &ROMol, at: u32) -> bool {
    let a = mol.atom(at);
    let mut deg = mol.atom_degree(at) as i32 + a.num_explicit_hs() as i32;
    for b in mol.atom_bond_indices(at) {
        if valence_contrib(mol, b, at).round() as i64 == 0 {
            deg -= 1;
        }
    }
    a.valence(ValenceType::Explicit).unwrap_or(0) as i32 != deg
}

/// `MolOps::isBondOrderQuery` (no query bonds yet)
fn is_bond_order_query(_mol: &ROMol, _bond: u32) -> bool {
    false
}

/// `MolOps::countAtomElec`: electrons an atom can donate to a conjugated
/// system (-1 if it cannot take part).
pub fn count_atom_elec(mol: &ROMol, at: u32) -> i32 {
    let pt = PeriodicTable::get();
    let a = mol.atom(at);
    let dv = pt.element(a.atomic_num()).valences[0];
    if dv <= 1 {
        return -1;
    }
    let mut degree = mol.atom_degree(at) as i32 + mol.total_num_hs(at, false).unwrap_or(0) as i32;
    for b in mol.atom_bond_indices(at) {
        if valence_contrib(mol, b, at) == 0.0 && !is_bond_order_query(mol, b) {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }
    let nlp = pt.element(a.atomic_num()).n_outer_elecs - dv;
    let nlp = (nlp - a.formal_charge()).max(0);
    let n_radicals = a.num_radical_electrons() as i32;
    let mut res = (dv - degree) + nlp - n_radicals;
    if res > 1 {
        let n_unsat =
            a.valence(ValenceType::Explicit).unwrap_or(0) as i32 - mol.atom_degree(at) as i32;
        if n_unsat > 1 {
            res = 1;
        }
    }
    res
}

fn get_atom_donor_type_arom(mol: &ROMol, at: u32, exocyclic_steal: bool) -> ElectronDonorType {
    use ElectronDonorType::*;
    let a = mol.atom(at);
    if a.atomic_num() == 0 {
        return if incident_cyclic_multiple_bond(mol, at) {
            One
        } else {
            Any
        };
    }
    let mut nelec = count_atom_elec(mol, at);
    let pt = PeriodicTable::get();
    if nelec < 0 {
        No
    } else if nelec == 0 {
        if incident_non_cyclic_multiple_bond(mol, at).is_some() {
            Vacant
        } else if incident_cyclic_multiple_bond(mol, at) {
            One
        } else {
            No
        }
    } else if nelec == 1 {
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, at) {
            let at2 = mol.atom(who);
            if exocyclic_steal
                && pt
                    .more_electronegative(at2.atomic_num(), a.atomic_num())
                    .unwrap_or(false)
            {
                Vacant
            } else {
                One
            }
        } else if incident_multiple_bond(mol, at) {
            One
        } else if a.formal_charge() == 1 {
            Vacant
        } else {
            No
        }
    } else {
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, at) {
            let at2 = mol.atom(who);
            if exocyclic_steal
                && pt
                    .more_electronegative(at2.atomic_num(), a.atomic_num())
                    .unwrap_or(false)
            {
                nelec -= 1;
            }
        }
        if nelec % 2 == 1 { One } else { Two }
    }
}

#[allow(clippy::fn_params_excessive_bools)]
fn is_atom_cand_for_arom(
    mol: &ROMol,
    at: u32,
    edon: ElectronDonorType,
    allow_third_row: bool,
    allow_triple_bonds: bool,
    allow_higher_exceptions: bool,
    only_c_or_n: bool,
    allow_exocyclic_multiple_bonds: bool,
) -> bool {
    let a = mol.atom(at);
    let an = a.atomic_num();
    if only_c_or_n && an != 6 && an != 7 {
        return false;
    }
    if !allow_third_row && an > 10 {
        return false;
    }
    if an > 18 && (!allow_higher_exceptions || (an != 34 && an != 52)) {
        return false;
    }
    if edon == ElectronDonorType::No {
        return false;
    }
    let pt = PeriodicTable::get();
    let def_val = pt.element(an).valences[0];
    if def_val > 0 {
        let eff = an as i32 - a.formal_charge();
        let eff_dv = if eff >= 0 && eff as u32 <= pt.get_max_atomic_number() {
            pt.element(eff as u32).valences[0]
        } else {
            return false;
        };
        if a.total_valence().unwrap_or(0) as i32 > eff_dv {
            return false;
        }
    }
    if a.num_radical_electrons() != 0 && (an != 6 || a.formal_charge() != 0) {
        return false;
    }
    let n_unsat = a.valence(ValenceType::Explicit).unwrap_or(0) as i32 - mol.atom_degree(at) as i32;
    if n_unsat > 1 {
        let mut n_mult = 0;
        for b in mol.atom_bonds(at) {
            match b.bond_type() {
                BondType::Double => n_mult += 1,
                BondType::Triple => {
                    if !allow_triple_bonds {
                        return false;
                    }
                    n_mult += 1;
                }
                _ => {}
            }
            if n_mult > 1 {
                break;
            }
        }
        if n_mult > 1 {
            return false;
        }
    }
    if !allow_exocyclic_multiple_bonds {
        for b in mol.atom_bonds(at) {
            if matches!(b.bond_type(), BondType::Double | BondType::Triple)
                && mol.ring_info().num_bond_rings(b.idx()).unwrap_or(0) == 0
            {
                return false;
            }
        }
    }
    true
}

// ------------------------------------------------------------------ Hückel

fn apply_huckel(ring: &[i32], edon: &[ElectronDonorType], min_ring_size: usize) -> bool {
    if ring.len() < min_ring_size {
        return false;
    }
    let (mut rlw, mut rup) = (0, 0);
    let mut n_any = 0;
    for &idx in ring {
        let t = edon[idx as usize];
        if t == ElectronDonorType::Any {
            n_any += 1;
            if n_any > 1 {
                return false;
            }
        }
        let (lw, up) = min_max_atom_elecs(t);
        rlw += lw;
        rup += up;
    }
    if rup >= 6 {
        (rlw..=rup).any(|rie| (rie - 2) % 4 == 0)
    } else {
        rup == 2
    }
}

fn mark_atoms_bonds_arom(
    mol: &mut ROMol,
    brings: &[Vec<i32>],
    ring_ids: &[i32],
    done_bonds: &mut BTreeSet<u32>,
) {
    let mut cntr: BTreeMap<i32, i32> = BTreeMap::new();
    for &ri in ring_ids {
        for &bi in &brings[ri as usize] {
            *cntr.entry(bi).or_default() += 1;
        }
    }
    for (&bi, &count) in &cntr {
        if count == 1 {
            let (b, e, bt) = {
                let bond = mol.bond(bi as u32);
                (bond.begin_atom_idx(), bond.end_atom_idx(), bond.bond_type())
            };
            mol.bond_mut(bi as u32).set_is_aromatic(true);
            if matches!(bt, BondType::Single | BondType::Double) {
                mol.bond_mut(bi as u32).set_bond_type(BondType::Aromatic);
                mol.atom_mut(b).set_is_aromatic(true);
                mol.atom_mut(e).set_is_aromatic(true);
            }
            done_bonds.insert(bi as u32);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn apply_huckel_to_fused(
    mol: &mut ROMol,
    srings: &[Vec<i32>],
    brings: &[Vec<i32>],
    fused: &[i32],
    edon: &[ElectronDonorType],
    ring_neighs: &RingNeighborMap,
    narom: &mut i32,
    max_num_fused_rings: usize,
    min_ring_size: usize,
) {
    let mut arom_rings: HashSet<i32> = HashSet::new();
    let nrings = fused.len();
    let mut pos: i32 = -1;
    let mut cur_size = 0usize;
    let mut comb: Vec<i32> = Vec::new();
    let n_ring_bonds = {
        let mut fb = vec![false; mol.num_bonds() as usize];
        for &r in fused {
            for &b in &brings[r as usize] {
                fb[b as usize] = true;
            }
        }
        fb.iter().filter(|&&x| x).count()
    };
    let mut done_bonds: BTreeSet<u32> = BTreeSet::new();
    loop {
        if pos == -1 {
            if cur_size == 2 && nrings > 300 {
                crate::rd_log_warning!(
                    "Aromaticity detection halted on some rings due to ring system size."
                );
                break;
            }
            cur_size += 1;
            if cur_size > nrings.min(max_num_fused_rings) || done_bonds.len() >= n_ring_bonds {
                break;
            }
            comb = (0..cur_size as i32).collect();
            pos = 0;
        } else {
            pos = next_combination(&mut comb, nrings as i32);
        }
        if pos == -1 {
            continue;
        }
        let cur_rs: Vec<i32> = comb.iter().map(|&i| fused[i as usize]).collect();
        if !ring_neighs.is_empty() && !check_fused(&cur_rs, ring_neighs) {
            continue;
        }
        let mut ats = vec![0i32; mol.num_atoms() as usize];
        for &ridx in &cur_rs {
            for &a in &srings[ridx as usize] {
                ats[a as usize] += 1;
            }
        }
        let unon: Vec<i32> = (0..ats.len())
            .filter(|&i| ats[i] == 1 || ats[i] == 2)
            .map(|i| i as i32)
            .collect();
        if apply_huckel(&unon, edon, min_ring_size) {
            mark_atoms_bonds_arom(mol, brings, &cur_rs, &mut done_bonds);
            arom_rings.extend(cur_rs.iter().copied());
        }
    }
    *narom += arom_rings.len() as i32;
}

/// Iterate the fused systems of `c_rings` (`pickFusedRings` loop).
#[allow(clippy::too_many_arguments)]
fn huckel_over_fused_systems(
    mol: &mut ROMol,
    c_rings: &[Vec<i32>],
    brings: &[Vec<i32>],
    edon: &[ElectronDonorType],
    neigh_map: &RingNeighborMap,
    narom: &mut i32,
    max_fused: usize,
    min_ring_size: usize,
) {
    let cnrs = c_rings.len();
    let mut fus_done = vec![false; cnrs];
    let mut curr = 0usize;
    while curr < cnrs {
        let mut fused = Vec::new();
        pick_fused_rings(curr as i32, neigh_map, &mut fused, &mut fus_done);
        apply_huckel_to_fused(
            mol,
            c_rings,
            brings,
            &fused,
            edon,
            neigh_map,
            narom,
            max_fused,
            min_ring_size,
        );
        match fus_done.iter().position(|&d| !d) {
            Some(r) => curr = r,
            None => break,
        }
    }
}

fn is_atom_dummy(mol: &ROMol, at: u32) -> bool {
    let a = mol.atom(at);
    !a.has_query() && a.atomic_num() == 0
}

fn aromaticity_helper(
    mol: &mut RWMol,
    srings: &[Vec<i32>],
    min_ring_size: usize,
    max_ring_size: usize,
    include_fused: bool,
) -> Result<i32> {
    let natoms = mol.num_atoms() as usize;
    let mut acands = vec![false; natoms];
    let mut aseen = vec![false; natoms];
    let mut edon = vec![ElectronDonorType::Vacant; natoms];
    let mut c_rings: Vec<Vec<i32>> = Vec::new();
    for sring in srings {
        let sz = sring.len();
        if (min_ring_size != 0 && sz < min_ring_size) || (max_ring_size != 0 && sz > max_ring_size)
        {
            continue;
        }
        let mut all_aromatic = true;
        let mut all_dummy = true;
        for &first in sring {
            let fi = first as usize;
            if all_dummy && !is_atom_dummy(mol, first as u32) {
                all_dummy = false;
            }
            if aseen[fi] {
                if !acands[fi] {
                    all_aromatic = false;
                }
                continue;
            }
            aseen[fi] = true;
            let mut donor = get_atom_donor_type_arom(mol, first as u32, true);
            if donor == ElectronDonorType::Two && sz >= 9 {
                let a = mol.atom(first as u32);
                if (a.atomic_num() == 8 || a.atomic_num() == 16)
                    && mol.atom_degree(first as u32) == 2
                    && a.formal_charge() == 0
                    && !mol
                        .atom_bonds(first as u32)
                        .any(|b| matches!(b.bond_type(), BondType::Double | BondType::Triple))
                {
                    donor = ElectronDonorType::No;
                }
            }
            edon[fi] = donor;
            acands[fi] =
                is_atom_cand_for_arom(mol, first as u32, donor, true, true, true, false, true);
            if !acands[fi] {
                all_aromatic = false;
            }
        }
        if all_aromatic && !all_dummy {
            c_rings.push(sring.clone());
        }
    }
    let brings: Vec<Vec<i32>> = c_rings
        .iter()
        .map(|r| convert_to_bonds(r, mol))
        .collect::<Result<_>>()?;
    let mut narom = 0;
    if !include_fused {
        let empty = RingNeighborMap::new();
        for ri in 0..c_rings.len() {
            apply_huckel_to_fused(
                mol,
                &c_rings,
                &brings,
                &[ri as i32],
                &edon,
                &empty,
                &mut narom,
                6,
                0,
            );
        }
    } else {
        let neigh_map = make_ring_neighbor_map(&brings, MAX_FUSED_AROMATIC_RING_SIZE, 1);
        huckel_over_fused_systems(mol, &c_rings, &brings, &edon, &neigh_map, &mut narom, 6, 0);
    }
    mol.set_prop(cp::numArom, narom, true)?;
    Ok(narom)
}

fn mdl_aromaticity_helper(mol: &mut RWMol, srings: &[Vec<i32>]) -> Result<i32> {
    let natoms = mol.num_atoms() as usize;
    let mut acands = vec![false; natoms];
    let mut aseen = vec![false; natoms];
    let mut edon = vec![ElectronDonorType::Vacant; natoms];
    let mut c_rings: Vec<Vec<i32>> = Vec::new();
    for sring in srings {
        let mut all_aromatic = true;
        let mut all_dummy = true;
        for &first in sring {
            let fi = first as usize;
            if all_dummy && mol.atom(first as u32).atomic_num() != 0 {
                all_dummy = false;
            }
            if aseen[fi] {
                if !acands[fi] {
                    all_aromatic = false;
                }
                continue;
            }
            aseen[fi] = true;
            edon[fi] = get_atom_donor_type_arom(mol, first as u32, false);
            if edon[fi] != ElectronDonorType::One {
                all_aromatic = false;
                continue;
            }
            acands[fi] = is_atom_cand_for_arom(
                mol,
                first as u32,
                edon[fi],
                false,
                false,
                false,
                true,
                false,
            );
            if !acands[fi] {
                all_aromatic = false;
            }
        }
        if all_aromatic && !all_dummy {
            c_rings.push(sring.clone());
        }
    }
    let brings: Vec<Vec<i32>> = c_rings
        .iter()
        .map(|r| convert_to_bonds(r, mol))
        .collect::<Result<_>>()?;
    let neigh_map = make_ring_neighbor_map(&brings, MAX_FUSED_AROMATIC_RING_SIZE, 1);
    let mut narom = 0;
    huckel_over_fused_systems(mol, &c_rings, &brings, &edon, &neigh_map, &mut narom, 6, 6);
    mol.set_prop(cp::numArom, narom, true)?;
    Ok(narom)
}

/// `MolOps::setMMFFAromaticity`
pub fn set_mmff_aromaticity(mol: &mut RWMol) -> Result<()> {
    let atom_rings: Vec<Vec<i32>> = mol.ring_info().atom_rings().to_vec();
    let n = atom_rings.len();
    let mut arom_bit = vec![false; mol.num_atoms() as usize];
    let mut arom_ring_bit = vec![false; n];
    let mut n_arom_set: i64 = 0;
    let mut old_n_arom_set: i64 = -1;
    let mut all_set = false;
    let bt = |mol: &RWMol, a: i32, b: i32| {
        mol.bond(mol.bond_idx_between(a as u32, b as u32).expect("ring bond"))
            .bond_type()
    };
    while !all_set && n > 0 && n_arom_set > old_n_arom_set {
        for i in 0..n {
            let ring = &atom_rings[i];
            let mut pi_e: u32 = 0;
            let mut move_next = false;
            let mut nos_in_ring = false;
            let mut exo_double = false;
            let mut j = 0;
            while !move_next && j < ring.len() {
                let a = mol.atom(ring[j] as u32);
                let an = a.atomic_num();
                if an == 7 || an == 8 || (an == 16 && mol.atom_degree(ring[j] as u32) == 2) {
                    nos_in_ring = true;
                }
                let next_in_ring = if j == ring.len() - 1 {
                    ring[0]
                } else {
                    ring[j + 1]
                };
                if bt(mol, ring[j], next_in_ring) == BondType::Double {
                    pi_e += 2;
                } else {
                    let four_n = an == 7
                        && a.valence(ValenceType::Explicit).unwrap_or(0)
                            + a.num_implicit_hs().unwrap_or(0)
                            == 4;
                    if an != 6 && !four_n {
                        j += 1;
                        continue;
                    }
                    for nbr in mol.atom_neighbors(ring[j] as u32) {
                        if ring.contains(&(nbr as i32)) {
                            continue;
                        }
                        let nb_bt = bt(mol, ring[j], nbr as i32);
                        if nb_bt == BondType::Single {
                            continue;
                        }
                        if mol.ring_info().num_atom_rings(nbr).unwrap_or(0) != 0
                            && !arom_bit[nbr as usize]
                        {
                            move_next = true;
                            break;
                        }
                        if nb_bt == BondType::Double {
                            if mol.atom(nbr).is_aromatic() {
                                pi_e += 1;
                            } else {
                                exo_double = true;
                            }
                        }
                    }
                }
                j += 1;
            }
            if move_next {
                continue;
            }
            let mut can_be_aromatic = true;
            for &a in ring {
                arom_bit[a as usize] = true;
                let at = mol.atom(a as u32);
                if (at.atomic_num() == 6 || at.atomic_num() == 7)
                    && at.hybridization() != super::atom::HybridizationType::SP2
                {
                    can_be_aromatic = false;
                }
            }
            if !can_be_aromatic {
                continue;
            }
            if nos_in_ring && !exo_double && ring.len() % 2 == 1 {
                pi_e += 2;
            }
            if pi_e > 2 && (pi_e - 2) % 4 == 0 {
                arom_ring_bit[i] = true;
                for &a in ring {
                    mol.atom_mut(a as u32).set_is_aromatic(true);
                }
            }
        }
        old_n_arom_set = n_arom_set;
        n_arom_set = 0;
        all_set = true;
        for ring in &atom_rings {
            for &a in ring {
                if arom_bit[a as usize] {
                    n_arom_set += 1;
                } else {
                    all_set = false;
                }
            }
        }
    }
    for i in 0..n {
        if !arom_ring_bit[i] {
            continue;
        }
        let ring = &atom_rings[i];
        for j in 0..ring.len() {
            let next = if j == ring.len() - 1 {
                ring[0]
            } else {
                ring[j + 1]
            };
            let b = mol
                .bond_idx_between(ring[j] as u32, next as u32)
                .expect("ring bond");
            mol.bond_mut(b).set_bond_type(BondType::Aromatic);
            mol.bond_mut(b).set_is_aromatic(true);
        }
    }
    for i in 0..n {
        if !arom_ring_bit[i] {
            continue;
        }
        for &a in &atom_rings[i] {
            if mol.atom(a as u32).atomic_num() != 6 {
                let iv = mol.calc_implicit_valence(a as u32, false)?;
                mol.calc_explicit_valence(a as u32, false)?;
                if iv != 0 {
                    mol.atom_mut(a as u32).set_num_explicit_hs(iv as u32);
                    mol.calc_implicit_valence(a as u32, false)?;
                }
            }
        }
    }
    Ok(())
}

fn mmff94_aromaticity_helper(mol: &mut RWMol, srings: &[Vec<i32>]) -> Result<i32> {
    if !mol.has_prop(cp::_MMFFSanitized) {
        if mol.atoms().iter().any(|a| a.is_aromatic()) {
            super::kekulize::kekulize(mol, true, true, 100)?;
        }
        mol.set_prop(cp::_MMFFSanitized, 1i32, true)?;
    }
    set_mmff_aromaticity(mol)?;
    let mut narom = 1;
    for sring in srings {
        if sring.iter().all(|&a| mol.atom(a as u32).is_aromatic()) {
            narom += 1;
        }
    }
    Ok(narom)
}

/// `MolOps::AromaticityModel`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    Default = 0x0,
    Rdkit = 0x1,
    Simple = 0x2,
    Mdl = 0x4,
    Mmff94 = 0x8,
    Custom = 0xFFFFFFF,
}

/// `MolOps::setAromaticity(mol, model)`
pub fn set_aromaticity(mol: &mut RWMol, model: AromaticityModel) -> Result<i32> {
    let srings: Vec<Vec<i32>> = if mol.ring_info().is_initialized() {
        mol.ring_info().atom_rings().to_vec()
    } else {
        symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?
    };
    match model {
        AromaticityModel::Default | AromaticityModel::Rdkit => {
            aromaticity_helper(mol, &srings, 0, 0, true)
        }
        AromaticityModel::Simple => aromaticity_helper(mol, &srings, 5, 6, false),
        AromaticityModel::Mdl => mdl_aromaticity_helper(mol, &srings),
        AromaticityModel::Mmff94 => mmff94_aromaticity_helper(mol, &srings),
        AromaticityModel::Custom => Err(RDKitError::precondition(
            "function must be set when aromaticity model is AROMATICITY_CUSTOM",
        )),
    }
}
