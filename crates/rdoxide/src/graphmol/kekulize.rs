//! Ported from RDKit `Code/GraphMol/Kekulize.cpp`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, VecDeque};

use super::aromaticity::{make_ring_neighbor_map, pick_fused_rings};
use super::atom::{ValenceType, is_aromatic_atom, is_early_atom};
use super::bond::{BondDir, BondType};
use super::find_rings::{convert_to_bonds, find_sssr};
use super::periodic_table::PeriodicTable;
use super::queries::has_bond_type_query;
use super::romol::RWMol;
use crate::rdgeneral::utils::union_rings;
use crate::{RDKitError, Result, SanitizeError};

/// Default maximum number of backtracking steps.
pub const DEFAULT_MAX_BACKTRACKS: u32 = 100;

fn back_track(
    mol: &mut RWMol,
    last_opt: i32,
    done: &mut Vec<i32>,
    aqueue: &mut VecDeque<i32>,
    d_bnd_cands: &mut [bool],
    d_bnd_adds: &mut [bool],
) {
    let ei = done
        .iter()
        .position(|&d| d == last_opt)
        .unwrap_or(done.len());
    let tdone: Vec<i32> = done[..ei].to_vec();
    // put every atom done since the last occurrence of lastOpt (inclusive) back
    // at the front of the queue, keeping their order
    let eri = done.iter().rposition(|&d| d == last_opt).unwrap_or(0);
    for &a in done[eri..].iter().rev() {
        aqueue.push_front(a);
    }
    for bi in 0..mol.num_bonds() {
        if d_bnd_adds[bi as usize] {
            let (a1, a2) = (
                mol.bond(bi).begin_atom_idx() as i32,
                mol.bond(bi).end_atom_idx() as i32,
            );
            if !tdone.contains(&a1) && !tdone.contains(&a2) {
                d_bnd_adds[bi as usize] = false;
                mol.bond_mut(bi).set_bond_type(BondType::Single);
                d_bnd_cands[a1 as usize] = true;
                d_bnd_cands[a2 as usize] = true;
            }
        }
    }
    *done = tdone;
}

fn mark_dbond_cands(
    mol: &mut RWMol,
    all_atms: &[i32],
    d_bnd_cands: &mut [bool],
    questions: &mut Vec<i32>,
    done: &mut Vec<i32>,
) {
    let has_arom_or_dummy = all_atms
        .iter()
        .any(|&a| mol.atom(a as u32).atomic_num() == 0 || is_aromatic_atom(mol, a as u32));
    if !has_arom_or_dummy {
        return;
    }
    let ri = mol.ring_info();
    let atom_rings = ri.atom_rings().to_vec();
    let mut is_ring_not_cand = vec![false; atom_rings.len()];
    for (r, aring) in atom_rings.iter().enumerate() {
        is_ring_not_cand[r] = true;
        for &ai in aring {
            if is_aromatic_atom(mol, ai as u32)
                && mol.ring_info().num_atom_rings(ai as u32).unwrap_or(0) == 1
            {
                is_ring_not_cand[r] = false;
                break;
            }
        }
    }
    let pt = PeriodicTable::get();
    let mut make_single: Vec<u32> = Vec::new();
    let mut in_all = vec![false; mol.num_atoms() as usize];
    for &aa in all_atms {
        in_all[aa as usize] = true;
        let at_idx = aa as u32;
        let at = mol.atom(at_idx);
        if at.atomic_num() != 0 && !is_aromatic_atom(mol, at_idx) {
            done.push(aa);
            continue;
        }
        let mut sbo: i32 = 0;
        let mut n_to_ignore: i32 = 0;
        let mut non_ar_non_dummy_nbr: u32 = 0;
        for bi in mol.atom_bond_indices(at_idx) {
            let bond = mol.bond(bi);
            let other = bond.other_atom_idx(at_idx).expect("in bond");
            let oa = mol.atom(other);
            if oa.atomic_num() != 0 && !oa.is_aromatic() && in_all[other as usize] {
                non_ar_non_dummy_nbr += 1;
            }
            if bond.is_aromatic()
                && matches!(
                    bond.bond_type(),
                    BondType::Single | BondType::Double | BondType::Aromatic
                )
            {
                sbo += 1;
                make_single.push(bi);
            } else {
                let contrib = bond.valence_contrib(mol, at_idx).round() as i32;
                sbo += contrib;
                if contrib == 0 {
                    n_to_ignore += 1;
                }
            }
        }
        let num_atom_rings = mol.ring_info().num_atom_rings(at_idx).unwrap_or(0);
        let num_non_cand = mol
            .ring_info()
            .atom_members(at_idx)
            .unwrap_or(&[])
            .iter()
            .filter(|&&r| is_ring_not_cand[r as usize])
            .count() as u32;
        if at.atomic_num() == 0
            && non_ar_non_dummy_nbr < num_atom_rings
            && num_non_cand < num_atom_rings
        {
            d_bnd_cands[aa as usize] = true;
            questions.push(aa);
        } else {
            sbo += mol.total_num_hs(at_idx, false).unwrap_or(0) as i32;
            let mut dv = pt.element(at.atomic_num()).valences[0];
            let mut chrg = at.formal_charge();
            if is_early_atom(at.atomic_num() as i32) {
                chrg = -chrg; // GitHub #65
            }
            if at.atomic_num() == 6 && chrg > 0 {
                chrg = -chrg;
            }
            dv += chrg;
            let tbo = at.total_valence().unwrap_or(0) as i32;
            let n_radicals = at.num_radical_electrons() as i32;
            let total_degree = mol.atom_degree(at_idx) as i32
                + at.valence(ValenceType::Implicit).unwrap_or(0) as i32
                - n_to_ignore;
            let val_list = &pt.element(at.atomic_num()).valences;
            let mut vi = 1;
            while tbo > dv && vi < val_list.len() && val_list[vi] > 0 {
                dv = val_list[vi] + chrg;
                vi += 1;
            }
            if tbo == 5
                && sbo == 4
                && dv == 3
                && total_degree == 3
                && n_radicals == 0
                && chrg == 0
                && mol.total_num_hs(at_idx, false).unwrap_or(0) == 0
                && matches!(at.atomic_num(), 7 | 15 | 33)
            {
                dv = 5;
            }
            if total_degree + n_radicals >= dv {
                continue;
            }
            if dv == sbo + 1 + n_radicals {
                d_bnd_cands[aa as usize] = true;
            } else if n_radicals == 0 && at.no_implicit() && dv == sbo + 2 {
                d_bnd_cands[aa as usize] = true;
            }
        }
    }
    for bi in make_single {
        mol.bond_mut(bi).set_bond_type(BondType::Single);
    }
}

#[allow(clippy::too_many_arguments)]
fn kekulize_worker(
    mol: &mut RWMol,
    all_atms: &[i32],
    mut d_bnd_cands: Vec<bool>,
    mut d_bnd_adds: Vec<bool>,
    mut done: Vec<i32>,
    atom_ranks: &[u32],
    max_back_tracks: u32,
) -> bool {
    let mut astack: VecDeque<i32> = VecDeque::new();
    let mut options: BTreeMap<i32, VecDeque<i32>> = BTreeMap::new();
    let mut last_opt: i32 = -1;
    let mut local_bonds_added = vec![false; mol.num_bonds() as usize];
    let mut in_all = vec![false; mol.num_atoms() as usize];
    for &a in all_atms {
        in_all[a as usize] = true;
    }
    let less_by_rank = |a: i32, b: i32| -> std::cmp::Ordering {
        let (ra, rb) = (atom_ranks[a as usize], atom_ranks[b as usize]);
        ra.cmp(&rb).then(a.cmp(&b))
    };
    let mut wedge_end = vec![false; mol.num_atoms() as usize];
    for b in mol.bonds() {
        if matches!(b.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash)
            && in_all[b.end_atom_idx() as usize]
        {
            wedge_end[b.end_atom_idx() as usize] = true;
        }
    }
    let mut sorted: Vec<i32> = all_atms.to_vec();
    sorted.sort_by(|&a, &b| {
        let (wa, wb) = (wedge_end[a as usize], wedge_end[b as usize]);
        if wa != wb {
            return wb.cmp(&wa);
        }
        less_by_rank(a, b)
    });
    let mut btmoves: VecDeque<i32> = VecDeque::new();
    let mut num_bt = 0;
    while done.len() < sorted.len() || !astack.is_empty() {
        let curr = if let Some(c) = astack.pop_front() {
            c
        } else {
            match sorted.iter().find(|a| !done.contains(a)) {
                Some(&c) => c,
                None => return false,
            }
        };
        done.push(curr);
        let c_cand = d_bnd_cands[curr as usize];
        let mut opts: VecDeque<i32>;
        if let Some(o) = options.get(&curr) {
            opts = o.clone();
        } else {
            opts = VecDeque::new();
            let mut lstack: Vec<i32> = Vec::new();
            let mut opts_v = Vec::new();
            let mut wedged_opts_v = Vec::new();
            let mut nbrs: Vec<i32> = mol
                .atom_neighbors(curr as u32)
                .map(|n| n as i32)
                .filter(|&n| in_all[n as usize] && !done.contains(&n))
                .collect();
            nbrs.sort_by(|&a, &b| less_by_rank(a, b));
            for nbr in nbrs {
                let bi = mol
                    .bond_idx_between(curr as u32, nbr as u32)
                    .expect("bonded");
                let nbond = mol.bond(bi);
                if !astack.contains(&nbr) {
                    lstack.push(nbr);
                }
                if c_cand
                    && d_bnd_cands[nbr as usize]
                    && (nbond.is_aromatic()
                        || mol.atom(curr as u32).atomic_num() == 0
                        || mol.atom(nbr as u32).atomic_num() == 0)
                {
                    if matches!(nbond.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash) {
                        wedged_opts_v.push(nbr);
                    } else {
                        opts_v.push(nbr);
                    }
                }
            }
            opts.extend(opts_v);
            opts.extend(wedged_opts_v);
            astack.extend(lstack);
        }
        if c_cand {
            if let Some(ncnd) = opts.pop_front() {
                let bi = mol
                    .bond_idx_between(curr as u32, ncnd as u32)
                    .expect("bonded");
                let bnd = mol.bond_mut(bi);
                bnd.set_bond_type(BondType::Double);
                if bnd.bond_dir() != BondDir::None {
                    bnd.set_bond_dir(BondDir::None);
                }
                d_bnd_cands[curr as usize] = false;
                d_bnd_cands[ncnd as usize] = false;
                d_bnd_adds[bi as usize] = true;
                local_bonds_added[bi as usize] = true;
                if options.contains_key(&curr) {
                    if opts.is_empty() {
                        options.remove(&curr);
                        btmoves.pop_back();
                        last_opt = btmoves.back().copied().unwrap_or(-1);
                    } else {
                        options.insert(curr, opts);
                    }
                } else if !opts.is_empty() {
                    last_opt = curr;
                    btmoves.push_back(last_opt);
                    options.insert(curr, opts);
                }
            } else if mol.atom(curr as u32).atomic_num() != 0 {
                if last_opt >= 0 && num_bt < max_back_tracks {
                    back_track(
                        mol,
                        last_opt,
                        &mut done,
                        &mut astack,
                        &mut d_bnd_cands,
                        &mut d_bnd_adds,
                    );
                    num_bt += 1;
                } else {
                    for (bi, &added) in local_bonds_added.iter().enumerate() {
                        if added {
                            mol.bond_mut(bi as u32).set_bond_type(BondType::Single);
                        }
                    }
                    return false;
                }
            }
        }
    }
    true
}

/// Enumerates subsets of the "question" (dummy) atoms to switch off.
struct QuestionEnumerator {
    questions: Vec<i32>,
    state: Vec<bool>,
    done: bool,
}

impl QuestionEnumerator {
    fn new(questions: Vec<i32>) -> Self {
        let mut state = vec![false; questions.len()];
        if let Some(s) = state.first_mut() {
            *s = true;
        }
        QuestionEnumerator {
            questions,
            state,
            done: false,
        }
    }

    fn next(&mut self) -> Vec<i32> {
        if self.done {
            return Vec::new();
        }
        let res: Vec<i32> = self
            .questions
            .iter()
            .zip(&self.state)
            .filter(|&(_, &s)| s)
            .map(|(&q, _)| q)
            .collect();
        let mut pos = 0;
        while pos < self.state.len() && self.state[pos] {
            self.state[pos] = false;
            pos += 1;
        }
        if pos == self.state.len() {
            self.done = true;
        } else {
            self.state[pos] = true;
        }
        res
    }
}

fn permute_dummies_and_kekulize(
    mol: &mut RWMol,
    all_atms: &[i32],
    d_bnd_cands: &[bool],
    questions: &[i32],
    atom_ranks: &[u32],
    max_back_tracks: u32,
) -> bool {
    let mut in_play = vec![false; mol.num_atoms() as usize];
    for &a in all_atms {
        in_play[a as usize] = true;
    }
    let mut kekulized = false;
    let mut q = QuestionEnumerator::new(questions.to_vec());
    while !kekulized && !questions.is_empty() {
        let d_bnd_adds = vec![false; mol.num_bonds() as usize];
        for bi in 0..mol.num_bonds() {
            let b = mol.bond(bi);
            if b.is_aromatic()
                && b.bond_type() != BondType::Single
                && in_play[b.begin_atom_idx() as usize]
                && in_play[b.end_atom_idx() as usize]
            {
                mol.bond_mut(bi).set_bond_type(BondType::Single);
            }
        }
        let switch_off = q.next();
        if switch_off.is_empty() {
            break;
        }
        let mut t_cands = d_bnd_cands.to_vec();
        for &s in &switch_off {
            t_cands[s as usize] = false;
        }
        kekulized = kekulize_worker(
            mol,
            all_atms,
            t_cands,
            d_bnd_adds,
            Vec::new(),
            atom_ranks,
            max_back_tracks,
        );
    }
    kekulized
}

fn kekulize_fused(
    mol: &mut RWMol,
    arings: &[Vec<i32>],
    atom_ranks: &[u32],
    max_back_tracks: u32,
) -> Result<()> {
    let all_atms = union_rings(arings);
    let mut done = Vec::new();
    let mut questions = Vec::new();
    let mut d_bnd_cands = vec![false; mol.num_atoms() as usize];
    let d_bnd_adds = vec![false; mol.num_bonds() as usize];
    mark_dbond_cands(mol, &all_atms, &mut d_bnd_cands, &mut questions, &mut done);
    let mut kekulized = kekulize_worker(
        mol,
        &all_atms,
        d_bnd_cands.clone(),
        d_bnd_adds,
        done,
        atom_ranks,
        max_back_tracks,
    );
    if !kekulized && !questions.is_empty() {
        kekulized = permute_dummies_and_kekulize(
            mol,
            &all_atms,
            &d_bnd_cands,
            &questions,
            atom_ranks,
            max_back_tracks,
        );
    }
    if !kekulized {
        let problem: Vec<u32> = (0..mol.num_atoms())
            .filter(|&i| d_bnd_cands[i as usize])
            .collect();
        let mut msg = String::from("Can't kekulize mol.  Unkekulized atoms:");
        for p in &problem {
            msg.push_str(&format!(" {p}"));
        }
        crate::rd_log_error!("{msg}");
        return Err(RDKitError::Sanitize(SanitizeError::Kekulize {
            atom_indices: problem,
            message: msg,
        }));
    }
    Ok(())
}

fn atom_kekulize_error(atom_idx: u32, message: String) -> RDKitError {
    crate::rd_log_error!("{message}");
    RDKitError::Sanitize(SanitizeError::AtomKekulize { atom_idx, message })
}

/// `MolOps::details::KekulizeFragment`
pub fn kekulize_fragment(
    mol: &mut RWMol,
    atoms_to_use: &[bool],
    bonds_to_use: &[bool],
    mark_atoms_bonds: bool,
    canonical: bool,
    max_back_tracks: u32,
) -> Result<()> {
    if !atoms_to_use.iter().any(|&a| a) {
        return Ok(());
    }
    let mut bonds_to_use = bonds_to_use.to_vec();
    // there's no point doing kekulization if there are no aromatic bonds
    // without queries
    let mut found_aromatic = false;
    for bond in mol.bonds() {
        let idx = bond.idx() as usize;
        if !bonds_to_use[idx] {
            continue;
        }
        if bond.query().is_some_and(has_bond_type_query) {
            // we don't kekulize bonds with bond type queries
            bonds_to_use[idx] = false;
        } else if bond.is_aromatic() {
            found_aromatic = true;
        }
    }
    let n = mol.num_atoms();
    let mut valences = vec![0i32; n as usize];
    let mut dummy_ats = vec![false; n as usize];
    for i in 0..n {
        if !atoms_to_use[i as usize] {
            continue;
        }
        mol.calc_implicit_valence(i, false)?;
        valences[i as usize] = mol.atom(i).total_valence()? as i32;
        if is_aromatic_atom(mol, i) {
            found_aromatic = true;
        }
        if mol.atom(i).atomic_num() == 0 {
            dummy_ats[i as usize] = true;
        }
    }
    if !found_aromatic {
        return Ok(());
    }
    let atom_ranks: Vec<u32> = if canonical {
        crate::graphmol::new_canon::rank_fragment_atoms(
            mol,
            atoms_to_use,
            &bonds_to_use,
            &crate::graphmol::new_canon::FragmentRankOptions::default(),
        )?
    } else {
        (0..n).collect()
    };
    if bonds_to_use.iter().any(|&b| b) {
        let mut wedged = vec![false; n as usize];
        for b in mol.bonds() {
            if bonds_to_use[b.idx() as usize]
                && matches!(b.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash)
            {
                wedged[b.begin_atom_idx() as usize] = true;
            }
        }
        let allrings: Vec<Vec<i32>> = if !mol.ring_info().is_initialized() {
            find_sssr(mol, false, false)?
        } else {
            mol.ring_info().atom_rings().to_vec()
        };
        let mut tmp: VecDeque<Vec<i32>> = VecDeque::new();
        for ring in &allrings {
            let mut ok = false;
            let mut all_in = true;
            for &ai in ring {
                if !atoms_to_use[ai as usize] {
                    all_in = false;
                    break;
                }
                if !dummy_ats[ai as usize] {
                    ok = true;
                }
            }
            if !(all_in && ok) {
                continue;
            }
            let start = ring.iter().position(|&a| wedged[a as usize]);
            let sp = start.unwrap_or(0);
            let nring: Vec<i32> = (0..ring.len())
                .map(|ri| ring[(ri + sp) % ring.len()])
                .collect();
            if start.is_none() {
                tmp.push_back(nring);
            } else {
                tmp.push_front(nring);
            }
        }
        let arings_all: Vec<Vec<i32>> = tmp.into_iter().collect();
        let allbrings: Vec<Vec<i32>> = arings_all
            .iter()
            .map(|r| convert_to_bonds(r, mol))
            .collect::<Result<_>>()?;
        let mut brings = Vec::new();
        let mut arings = Vec::new();
        for (i, br) in allbrings.iter().enumerate() {
            if br.iter().all(|&bi| bonds_to_use[bi as usize]) {
                brings.push(br.clone());
                arings.push(arings_all[i].clone());
            }
        }
        let neigh_map = make_ring_neighbor_map(&brings, 0, 0);
        let cnrs = arings.len();
        let mut fus_done = vec![false; cnrs];
        let mut curr = 0;
        while curr < cnrs {
            let mut fused = Vec::new();
            pick_fused_rings(curr as i32, &neigh_map, &mut fused, &mut fus_done);
            let frings: Vec<Vec<i32>> = fused.iter().map(|&r| arings[r as usize].clone()).collect();
            kekulize_fused(mol, &frings, &atom_ranks, max_back_tracks)?;
            match fus_done.iter().position(|&d| !d) {
                Some(r) => curr = r,
                None => break,
            }
        }
    }
    if mark_atoms_bonds {
        if !mol.ring_info().is_initialized() {
            find_sssr(mol, false, false)?;
        }
        for bi in 0..mol.num_bonds() {
            if bonds_to_use[bi as usize] {
                mol.bond_mut(bi).set_is_aromatic(false);
            }
        }
        let all_atoms = atoms_to_use.iter().all(|&a| a);
        let all_bonds = bonds_to_use.iter().all(|&b| b);
        for i in 0..n {
            if !(atoms_to_use[i as usize] && mol.atom(i).is_aromatic()) {
                continue;
            }
            if all_atoms && all_bonds && mol.ring_info().num_atom_rings(i)? == 0 {
                return Err(atom_kekulize_error(
                    i,
                    format!("non-ring atom {i} marked aromatic"),
                ));
            }
            mol.atom_mut(i).set_is_aromatic(false);
            let a = mol.atom(i);
            if (a.atomic_num() == 7 || a.atomic_num() == 15)
                && a.formal_charge() == 0
                && a.num_explicit_hs() == 1
            {
                let a = mol.atom_mut(i);
                a.set_no_implicit(false);
                a.set_num_explicit_hs(0);
                mol.update_atom_property_cache(i, false)?;
            }
        }
    }
    for i in 0..n {
        if !atoms_to_use[i as usize] {
            continue;
        }
        let val = mol.atom(i).total_valence()? as i32;
        if val != valences[i as usize] {
            return Err(atom_kekulize_error(
                i,
                format!(
                    "Kekulization somehow screwed up valence on {i}: {val}!={}\n",
                    valences[i as usize]
                ),
            ));
        }
    }
    Ok(())
}

/// `MolOps::Kekulize`
pub fn kekulize(
    mol: &mut RWMol,
    mark_atoms_bonds: bool,
    canonical: bool,
    max_back_tracks: u32,
) -> Result<()> {
    let atoms = vec![true; mol.num_atoms() as usize];
    let bonds = vec![true; mol.num_bonds() as usize];
    kekulize_fragment(
        mol,
        &atoms,
        &bonds,
        mark_atoms_bonds,
        canonical,
        max_back_tracks,
    )
}

/// `MolOps::KekulizeIfPossible`: on failure restores aromatic flags and
/// returns false.
pub fn kekulize_if_possible(
    mol: &mut RWMol,
    mark_atoms_bonds: bool,
    canonical: bool,
    max_back_tracks: u32,
) -> Result<bool> {
    let arom_bonds: Vec<bool> = mol.bonds().iter().map(|b| b.is_aromatic()).collect();
    let arom_atoms: Vec<bool> = (0..mol.num_atoms())
        .map(|i| is_aromatic_atom(mol, i))
        .collect();
    match kekulize(mol, mark_atoms_bonds, canonical, max_back_tracks) {
        Ok(()) => Ok(true),
        Err(RDKitError::Sanitize(_)) => {
            for (i, &a) in arom_bonds.iter().enumerate() {
                if a {
                    let b = mol.bond_mut(i as u32);
                    b.set_is_aromatic(true);
                    b.set_bond_type(BondType::Aromatic);
                }
            }
            for (i, &a) in arom_atoms.iter().enumerate() {
                if a {
                    mol.atom_mut(i as u32).set_is_aromatic(true);
                }
            }
            Ok(false)
        }
        Err(e) => Err(e),
    }
}
