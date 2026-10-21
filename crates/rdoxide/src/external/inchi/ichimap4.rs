//! Mapping of stereo bonds and stereocenters that finds the smallest stereo
//! part of the canonical connection table: `map_stereo_bonds4` and
//! `map_stereo_atoms4` of `ichimap4.c`.
//!
//! The C code walks `AT_RANK **` rank stacks with pointer arithmetic
//! (`pRankStack1 + nStackPtr[i]`); here a stack is a `Vec<Vec<u16>>` and a
//! deeper level works on the tail split off at that offset (and appended back
//! afterwards). The `CANON_STAT` members used are in [`StereoMapCs`]. The
//! global user-quit callbacks of the C code are not ported, and `nSymmStereo`
//! is always present.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSE/InChI-MIT.txt).

use super::ichierr::{
    CT_MAPCOUNT_ERR, CT_OVERFLOW, CT_STEREO_CANON_ERR, CT_STEREOBOND_ERROR, CT_STEREOCOUNT_ERR,
    CT_TIMEOUT_ERR,
};
use super::ichimap::{
    BEST_PARITY, CanonCounters, CurTree, EqNeigh, NeighList, RemoveStereoCs, STEREO_AT_MARK,
    WORSE_PARITY, all_sb_same, all_sc_same, b_unique_at_nbr_from_mapping_rank, break_all_ties,
    clear_previous_mappings, compare_lin_ct_stereo_atom_to_values,
    compare_lin_ct_stereo_double_to_values, map_an_atom2, might_change_other_atom_parity,
    n_get_mcr, n_join2_mcrs, next_sb_at_canon_ranks2, next_sc_at_canon_rank2,
    next_stereo_parity2_test, parity_of_mapped_atom2, parity_of_mapped_half_bond,
    remove_calculated_non_stereo, set_use_atom_for_stereo,
};
use super::ichister::{
    AB_PARITY_CALC, AB_PARITY_UNDF, AB_PARITY_UNKN, AtStereoCarb, AtStereoDble, KNOWN_PARITIES_EQL,
    SpAtom, atom_parity_well_def, parity_known, parity_val,
};
use super::inpdef::MAX_NUM_STEREO_BONDS;

const SB_DEPTH: usize = 6;
const AB_PARITY_NONE: i32 = 0;
/// `CMODE_REDNDNT_STEREO`
pub const CMODE_REDNDNT_STEREO: i32 = 0x80;

/// The members of `CANON_STAT` used by `map_stereo_bonds4` and `map_stereo_atoms4`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StereoMapCs {
    /// stereo part of the linear CT (lengths are the vector lengths), `nSymmRank`
    pub rm: RemoveStereoCs,
    /// `lNumNeighListIter` and `lNumBreakTies`
    pub counters: CanonCounters,
    pub b_first_ct: i32,
    pub b_keep_symm_rank: i32,
    pub b_stereo_is_better: i32,
    pub l_num_rejected_ct: i64,
    pub l_num_tot_ct: i64,
    pub l_num_equal_ct: i64,
    pub l_num_decreased_ct: i64,
    pub b_rank_used_for_stereo: Vec<i8>,
    pub b_atom_used_for_stereo: Vec<i8>,
    pub n_prev_atom_number: Vec<u16>,
    pub n_mode: i32,
    /// result of `bInchiTimeIsOver( ic, ulTimeOutTime )`
    pub time_is_over: bool,
}

/// `map_an_atom2` on the rank stacks starting at `off`
#[allow(clippy::too_many_arguments)]
fn map_at(
    num_at_tg: usize,
    num_max: usize,
    from: usize,
    to: usize,
    n_temp_rank: &mut [u16],
    n_mapped_in: i32,
    n_mapped_out: &mut i32,
    cs: &mut StereoMapCs,
    neigh_list: &mut NeighList,
    st1: &mut Vec<Vec<u16>>,
    st2: &mut Vec<Vec<u16>>,
    off: usize,
    b_add_stack: &mut i32,
) -> i32 {
    let mut t1 = st1.split_off(off);
    let mut t2 = st2.split_off(off);
    let r = map_an_atom2(
        num_at_tg,
        num_max,
        from,
        to,
        n_temp_rank,
        n_mapped_in,
        n_mapped_out,
        &mut cs.counters,
        neigh_list,
        &mut t1,
        &mut t2,
        b_add_stack,
    );
    st1.append(&mut t1);
    st2.append(&mut t2);
    r
}

fn returned_error(v: i32) -> bool {
    (super::ichierr::CT_ERR_MIN..=super::ichierr::CT_ERR_MAX).contains(&v)
}

/// `map_stereo_atoms4`: maps the stereocenters (after all stereo bonds have
/// been mapped); returns 1 when a full mapping was completed, 2 added if a
/// better stereo CT was found, 4 for a restart, or an error
#[allow(clippy::too_many_arguments, unused_assignments)]
pub fn map_stereo_atoms4(
    at: &mut [SpAtom],
    num_atoms: usize,
    num_at_tg: usize,
    num_max: usize,
    n_canon_rank_from: &[u16],
    n_atom_number_canon_from: &[u16],
    n_canon_rank_to: &mut [u16],
    n_symm_rank: &[u16],
    st1: &mut Vec<Vec<u16>>,
    st2: &mut Vec<Vec<u16>>,
    n_temp_rank: &mut [u16],
    n_num_mapped_ranks_input: i32,
    n_symm_stereo: &mut [u16],
    neigh_list: &mut NeighList,
    cs: &mut StereoMapCs,
    cur_tree: &mut CurTree,
    n_num_mapped_atoms: usize,
    v_ab_parity_unknown: i32,
) -> i32 {
    let mut n_tot_success = 0i32;
    let mut prev_atom = AtStereoCarb::default();
    let tpos1 = cur_tree.get_pos();

    if n_num_mapped_atoms < cs.rm.linear_ct_stereo_carb.len() {
        prev_atom = cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms];
        let mut at_rank_canon1: u16 = if n_num_mapped_atoms != 0 {
            cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms - 1].at_num
        } else {
            0
        };
        let mut canon_rank1_min = 0u16;
        let mut b_first_time = true;
        let mut prev_atom2 = AtStereoCarb::default();
        let mut stereo_center_parity = 0i32;
        let mut sb_parity_calc = 0i32;
        // (these arrays live across the iterations, as in the C code)
        let mut n_stack_ptr = [0usize; 5];
        let mut n_mapped_ranks = [0i32; 5];
        let mut en = [EqNeigh::default(); 5];
        let mut jj = [0usize; 5];
        let mut cr = [0u16; 5];
        let mut at_to = [0u16; 5];
        let (mut n_num_unkn, mut n_num_undf, mut n_num_worse, mut n_num_best, mut n_num_calc) =
            (0i32, 0i32, 0i32, 0i32, 0i32);
        let mut pass = 0i32;
        let mut at_rank1 = 0usize;
        let mut at_from1 = 0usize;
        let mut i_max = 0i32;

        #[derive(Clone, Copy, PartialEq)]
        enum Pc {
            Bypass,
            NextCanonRank,
            RepeatAll,
        }
        let mut pc = Pc::Bypass;
        'sm: loop {
            match pc {
                Pc::NextCanonRank => {
                    if cs.b_stereo_is_better == 0
                        && at_rank_canon1 >= cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms].at_num
                    {
                        if n_tot_success == 0 {
                            cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms] = prev_atom;
                        }
                        cur_tree.set_pos(tpos1);
                        return n_tot_success;
                    }
                    pc = Pc::Bypass;
                    continue 'sm;
                }
                Pc::Bypass => {
                    cur_tree.set_pos(tpos1);
                    if !next_sc_at_canon_rank2(
                        &mut at_rank_canon1,
                        &mut canon_rank1_min,
                        &mut b_first_time,
                        &cs.b_atom_used_for_stereo,
                        st1,
                        st2,
                        n_atom_number_canon_from,
                        num_atoms as i32,
                    ) || (cs.b_stereo_is_better == 0
                        && at_rank_canon1 > cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms].at_num)
                    {
                        if n_tot_success == 0 {
                            cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms] = prev_atom;
                        }
                        return n_tot_success;
                    }
                    n_num_unkn = 0;
                    n_num_undf = 0;
                    n_num_best = 0;
                    n_num_worse = 0;
                    n_num_calc = 0;
                    let mut n_num_choices = 0i32;
                    pass = 0;
                    at_from1 = n_atom_number_canon_from[at_rank_canon1 as usize - 1] as usize;
                    at_rank1 = st1[0][at_from1] as usize;
                    i_max = at_rank1 as i32 - 1;
                    if at_rank1 as u16 != st2[0][st2[1][at_rank1 - 1] as usize] {
                        return CT_STEREOCOUNT_ERR;
                    }
                    for j1 in 0..=i_max {
                        let at_to1 = st2[1][(i_max - j1) as usize] as usize;
                        if at_rank1 as u16 != st2[0][at_to1] {
                            break;
                        }
                        if at[at_to1].stereo_bond_neighbor[0] == 0
                            && cs.b_atom_used_for_stereo[at_to1] == STEREO_AT_MARK
                        {
                            let mut no_choice = 0;
                            stereo_center_parity =
                                parity_val(i32::from(at[at_to1].stereo_atom_parity));
                            match stereo_center_parity {
                                AB_PARITY_UNDF => n_num_undf += 1,
                                AB_PARITY_UNKN => n_num_unkn += 1,
                                BEST_PARITY => n_num_best += 1,
                                WORSE_PARITY => n_num_worse += 1,
                                AB_PARITY_CALC => n_num_calc += 1,
                                AB_PARITY_NONE => no_choice += 1,
                                _ => {}
                            }
                            n_num_choices += i32::from(no_choice == 0);
                        }
                    }
                    if n_num_choices
                        != n_num_calc + n_num_undf + n_num_unkn + n_num_best + n_num_worse
                    {
                        return CT_STEREOCOUNT_ERR;
                    }
                    if n_num_choices == 0 {
                        pc = Pc::NextCanonRank;
                        continue 'sm;
                    }
                    sb_parity_calc = if n_num_calc > 0 { BEST_PARITY } else { 0 };
                    pc = Pc::RepeatAll;
                    continue 'sm;
                }
                Pc::RepeatAll => {
                    let first = pass == 0;
                    pass += 1;
                    if first {
                        stereo_center_parity = if sb_parity_calc != 0 {
                            BEST_PARITY
                        } else if n_num_best != 0 {
                            BEST_PARITY
                        } else if n_num_worse != 0 {
                            WORSE_PARITY
                        } else if n_num_unkn != 0 {
                            AB_PARITY_UNKN
                        } else if n_num_undf != 0 {
                            AB_PARITY_UNDF
                        } else {
                            AB_PARITY_NONE
                        };
                    } else {
                        let j1 = next_stereo_parity2_test(
                            &mut stereo_center_parity,
                            &mut sb_parity_calc,
                            n_num_best,
                            n_num_worse,
                            n_num_unkn,
                            n_num_undf,
                            n_num_calc,
                            v_ab_parity_unknown,
                        );
                        match j1 {
                            0 => {}
                            1 => {
                                pc = Pc::NextCanonRank;
                                continue 'sm;
                            }
                            _ => return j1,
                        }
                    }
                    if stereo_center_parity == AB_PARITY_NONE {
                        return CT_STEREOCOUNT_ERR;
                    }
                    if cs.b_stereo_is_better == 0 {
                        let cur: &AtStereoCarb = if n_tot_success != 0 {
                            &cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms]
                        } else {
                            &prev_atom
                        };
                        let c = compare_lin_ct_stereo_atom_to_values(
                            cur,
                            at_rank_canon1,
                            stereo_center_parity as u8,
                        );
                        if c < 0 {
                            if n_tot_success == 0 {
                                cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms] = prev_atom;
                            }
                            cur_tree.set_pos(tpos1);
                            return n_tot_success;
                        }
                    }
                    let mut b_stereo_is_better_was_set_here = false;
                    cur_tree.set_pos(tpos1);
                    for j1 in 0..=i_max {
                        let at_to1 = st2[1][(i_max - j1) as usize] as usize;
                        if at_rank1 as u16 != st2[0][at_to1] {
                            break;
                        }
                        let mut n_num_at_to1_success = 0i32;
                        if at[at_to1].stereo_atom_parity == 0
                            || at[at_to1].stereo_bond_neighbor[0] != 0
                            || cs.b_atom_used_for_stereo[at_to1] != STEREO_AT_MARK
                        {
                            continue;
                        }
                        let sap = i32::from(at[at_to1].stereo_atom_parity);
                        if parity_known(sap) {
                            if stereo_center_parity == sb_parity_calc {
                                continue;
                            }
                            if stereo_center_parity != parity_val(sap) {
                                continue;
                            }
                        } else if parity_val(sap) == AB_PARITY_CALC {
                            if stereo_center_parity != sb_parity_calc {
                                continue;
                            }
                        } else {
                            return CT_STEREOCOUNT_ERR;
                        }
                        let mut b_all_parities_identical =
                            (sap & KNOWN_PARITIES_EQL) != 0 && parity_known(sap);
                        if !b_all_parities_identical
                            && n_num_calc == 0
                            && (i32::from(n_num_undf == 0)
                                + i32::from(n_num_unkn == 0)
                                + i32::from(n_num_best == 0)
                                + i32::from(n_num_worse == 0))
                                == 3
                        {
                            let r =
                                all_sc_same(at_rank_canon1, st1, st2, n_atom_number_canon_from, at);
                            if r < 0 {
                                return CT_STEREOCOUNT_ERR;
                            }
                            b_all_parities_identical = r != 0;
                        }
                        if tpos1 < cur_tree.get_pos()
                            && cur_tree.is_last_rank(at_rank_canon1) == 1
                            && cur_tree.is_last_atom_equ(at_to1 as i32, n_symm_stereo) == 1
                        {
                            continue;
                        }
                        let mut istk = 0usize;
                        n_mapped_ranks[0] = n_num_mapped_ranks_input;
                        let mut b_add_stack = 0i32;
                        if !b_all_parities_identical {
                            let mut out = 0;
                            let ret1 = map_at(
                                num_at_tg,
                                num_max,
                                at_from1,
                                at_to1,
                                n_temp_rank,
                                n_mapped_ranks[istk],
                                &mut out,
                                cs,
                                neigh_list,
                                st1,
                                st2,
                                n_stack_ptr[istk],
                                &mut b_add_stack,
                            );
                            n_mapped_ranks[istk + 1] = out;
                            if returned_error(ret1) {
                                return ret1;
                            }
                            n_stack_ptr[istk + 1] = n_stack_ptr[istk] + b_add_stack as usize;
                            istk += 1;
                        } else {
                            if let Some(t) = st1.get_mut(2..) {
                                clear_previous_mappings(t);
                            }
                        }
                        let mut parity1;
                        if stereo_center_parity == sb_parity_calc {
                            parity1 = parity_of_mapped_atom2(
                                at_from1,
                                at_to1,
                                at,
                                Some(&mut en[istk]),
                                n_canon_rank_from,
                                &st1[n_stack_ptr[istk]],
                                &st2[n_stack_ptr[istk]],
                            );
                            if parity1 == 0 {
                                return CT_STEREOCOUNT_ERR;
                            }
                            if en[istk].num_to == 0 && parity1 != sb_parity_calc {
                                continue;
                            }
                        } else {
                            parity1 = stereo_center_parity;
                            en[istk].num_to = 0;
                        }
                        if (stereo_center_parity == sb_parity_calc && en[istk].num_to == 0)
                            || stereo_center_parity != sb_parity_calc
                        {
                            let c = compare_lin_ct_stereo_atom_to_values(
                                &cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms],
                                at_rank_canon1,
                                parity1 as u8,
                            );
                            if c < 0 && cs.b_stereo_is_better == 0 {
                                cs.l_num_rejected_ct += 1;
                                continue;
                            }
                            if b_add_stack != 0 {
                                if tpos1 == cur_tree.get_pos()
                                    || cur_tree.is_last_rank(at_rank_canon1) == 0
                                {
                                    cur_tree.add_rank(at_rank_canon1);
                                }
                                cur_tree.add_atom(at_to1 as i32);
                            }
                            if c > 0 && cs.b_stereo_is_better == 0 {
                                cs.b_stereo_is_better = 1;
                                b_stereo_is_better_was_set_here = true;
                                prev_atom2 = cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms];
                            }
                            cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms].parity = parity1 as u8;
                            cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms].at_num = at_rank_canon1;
                            cs.b_rank_used_for_stereo[at_from1] = 3;
                            if !b_all_parities_identical {
                                cs.b_atom_used_for_stereo[at_to1] -= STEREO_AT_MARK;
                            }
                            let ret = {
                                let off = n_stack_ptr[istk];
                                let mut t1 = st1.split_off(off);
                                let mut t2 = st2.split_off(off);
                                let r = map_stereo_atoms4(
                                    at,
                                    num_atoms,
                                    num_at_tg,
                                    num_max,
                                    n_canon_rank_from,
                                    n_atom_number_canon_from,
                                    n_canon_rank_to,
                                    n_symm_rank,
                                    &mut t1,
                                    &mut t2,
                                    n_temp_rank,
                                    n_mapped_ranks[istk],
                                    n_symm_stereo,
                                    neigh_list,
                                    cs,
                                    cur_tree,
                                    n_num_mapped_atoms + 1,
                                    v_ab_parity_unknown,
                                );
                                st1.append(&mut t1);
                                st2.append(&mut t2);
                                r
                            };
                            cs.b_rank_used_for_stereo[at_from1] = 0;
                            if !b_all_parities_identical {
                                cs.b_atom_used_for_stereo[at_to1] += STEREO_AT_MARK;
                            }
                            if ret == 4 {
                                return ret;
                            }
                            if returned_error(ret) {
                                return ret;
                            }
                            if ret > 0 {
                                n_tot_success |= 1;
                                n_num_at_to1_success += 1;
                                if b_stereo_is_better_was_set_here || (ret & 2) != 0 {
                                    cur_tree.keep_last_atoms_only(tpos1, 1);
                                    n_tot_success |= 2;
                                }
                            } else {
                                if b_stereo_is_better_was_set_here {
                                    cs.b_stereo_is_better = 0;
                                    cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms] = prev_atom2;
                                }
                                if tpos1 < cur_tree.get_pos()
                                    && cur_tree.is_last_rank(at_rank_canon1) == 1
                                {
                                    cur_tree.remove_if_last_atom(at_to1 as i32);
                                    cur_tree.remove_last_rank_if_no_atoms();
                                }
                            }
                            b_stereo_is_better_was_set_here = false;
                            if b_all_parities_identical {
                                break;
                            }
                            continue;
                        }
                        // need to map the neighbors
                        if stereo_center_parity != sb_parity_calc {
                            return CT_STEREOCOUNT_ERR;
                        }
                        if parity1 == 0 {
                            return CT_STEREOCOUNT_ERR;
                        }
                        if b_add_stack != 0 {
                            if tpos1 == cur_tree.get_pos()
                                || cur_tree.is_last_rank(at_rank_canon1) == 0
                            {
                                cur_tree.add_rank(at_rank_canon1);
                            }
                            cur_tree.add_atom(at_to1 as i32);
                        }
                        let mut lvl = 0usize;
                        let mut b_last_lvl_failed = false;
                        #[derive(Clone, Copy)]
                        enum Lv {
                            NextLvl,
                            NextJ,
                            Backup,
                        }
                        let mut lv = Lv::NextLvl;
                        'lv: loop {
                            match lv {
                                Lv::NextLvl => {
                                    if en[istk + lvl].num_to != 0 {
                                        jj[lvl] = 0;
                                        lv = Lv::NextJ;
                                        continue 'lv;
                                    }
                                    cr[lvl] = 0;
                                    break 'lv;
                                }
                                Lv::NextJ => {
                                    let e = en[istk + lvl];
                                    cr[lvl] = n_canon_rank_from[e.from_at as usize];
                                    at_to[lvl] = e.to_at[jj[lvl]];
                                    if jj[lvl] != 0
                                        && tpos1 < cur_tree.get_pos()
                                        && cur_tree.is_last_rank(cr[lvl]) == 1
                                        && cur_tree
                                            .is_last_atom_equ(i32::from(at_to[lvl]), n_symm_stereo)
                                            == 1
                                    {
                                        lvl += 1;
                                        b_last_lvl_failed = false;
                                        lv = Lv::Backup;
                                        continue 'lv;
                                    }
                                    let mut out = 0;
                                    let ret2 = map_at(
                                        num_at_tg,
                                        num_max,
                                        e.from_at as usize,
                                        e.to_at[jj[lvl]] as usize,
                                        n_temp_rank,
                                        n_mapped_ranks[istk + lvl],
                                        &mut out,
                                        cs,
                                        neigh_list,
                                        st1,
                                        st2,
                                        n_stack_ptr[istk + lvl],
                                        &mut b_add_stack,
                                    );
                                    n_mapped_ranks[istk + lvl + 1] = out;
                                    if returned_error(ret2) {
                                        return ret2;
                                    }
                                    if b_add_stack != 0 {
                                        if tpos1 == cur_tree.get_pos()
                                            || cur_tree.is_last_rank(cr[lvl]) == 0
                                        {
                                            cur_tree.add_rank(cr[lvl]);
                                        }
                                        cur_tree.add_atom(i32::from(at_to[lvl]));
                                    }
                                    n_stack_ptr[istk + lvl + 1] =
                                        n_stack_ptr[istk + lvl] + b_add_stack as usize;
                                    lvl += 1;
                                    parity1 = parity_of_mapped_atom2(
                                        at_from1,
                                        at_to1,
                                        at,
                                        Some(&mut en[istk + lvl]),
                                        n_canon_rank_from,
                                        &st1[n_stack_ptr[istk + lvl]],
                                        &st2[n_stack_ptr[istk + lvl]],
                                    );
                                    if parity1 == 0 {
                                        return CT_STEREOCOUNT_ERR;
                                    }
                                    if parity1 < 0 {
                                        lv = Lv::NextLvl;
                                        continue 'lv;
                                    }
                                    let c = compare_lin_ct_stereo_atom_to_values(
                                        &cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms],
                                        at_rank_canon1,
                                        parity1 as u8,
                                    );
                                    if sb_parity_calc != parity1
                                        || (c < 0 && cs.b_stereo_is_better == 0)
                                    {
                                        cs.l_num_rejected_ct += 1;
                                        b_last_lvl_failed = true;
                                    } else {
                                        if c > 0 && cs.b_stereo_is_better == 0 {
                                            cs.b_stereo_is_better = 1;
                                            b_stereo_is_better_was_set_here = true;
                                            prev_atom2 =
                                                cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms];
                                        }
                                        cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms].parity =
                                            parity1 as u8;
                                        cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms].at_num =
                                            at_rank_canon1;
                                        cs.b_rank_used_for_stereo[at_from1] = 3;
                                        cs.b_atom_used_for_stereo[at_to1] -= STEREO_AT_MARK;
                                        let ret = {
                                            let off = n_stack_ptr[istk + lvl];
                                            let mut t1 = st1.split_off(off);
                                            let mut t2 = st2.split_off(off);
                                            let r = map_stereo_atoms4(
                                                at,
                                                num_atoms,
                                                num_at_tg,
                                                num_max,
                                                n_canon_rank_from,
                                                n_atom_number_canon_from,
                                                n_canon_rank_to,
                                                n_symm_rank,
                                                &mut t1,
                                                &mut t2,
                                                n_temp_rank,
                                                n_mapped_ranks[istk + lvl],
                                                n_symm_stereo,
                                                neigh_list,
                                                cs,
                                                cur_tree,
                                                n_num_mapped_atoms + 1,
                                                v_ab_parity_unknown,
                                            );
                                            st1.append(&mut t1);
                                            st2.append(&mut t2);
                                            r
                                        };
                                        cs.b_rank_used_for_stereo[at_from1] = 0;
                                        cs.b_atom_used_for_stereo[at_to1] += STEREO_AT_MARK;
                                        if ret == 4 {
                                            return ret;
                                        }
                                        if returned_error(ret) {
                                            return ret;
                                        }
                                        if ret > 0 {
                                            n_tot_success |= 1;
                                            n_num_at_to1_success += 1;
                                            if b_stereo_is_better_was_set_here || (ret & 2) != 0 {
                                                cur_tree.keep_last_atoms_only(tpos1, 1);
                                                n_tot_success |= 2;
                                            }
                                        } else {
                                            if b_stereo_is_better_was_set_here {
                                                cs.b_stereo_is_better = 0;
                                                cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms] =
                                                    prev_atom2;
                                            }
                                            b_last_lvl_failed = true;
                                        }
                                        b_stereo_is_better_was_set_here = false;
                                        // (`!nSymmStereo` is never true here)
                                        let _ = might_change_other_atom_parity;
                                    }
                                    lv = Lv::Backup;
                                    continue 'lv;
                                }
                                Lv::Backup => loop {
                                    if lvl == 0 {
                                        break 'lv;
                                    }
                                    lvl -= 1;
                                    jj[lvl] += 1;
                                    if (jj[lvl] as i32) < en[istk + lvl].num_to {
                                        if b_last_lvl_failed {
                                            if tpos1 < cur_tree.get_pos()
                                                && cur_tree.is_last_rank(cr[lvl]) == 1
                                            {
                                                cur_tree.remove_if_last_atom(i32::from(at_to[lvl]));
                                                cur_tree.remove_last_rank_if_no_atoms();
                                            }
                                            b_last_lvl_failed = false;
                                        }
                                        lv = Lv::NextJ;
                                        continue 'lv;
                                    }
                                    if tpos1 < cur_tree.get_pos()
                                        && cur_tree.is_last_rank(cr[lvl]) == 1
                                    {
                                        cur_tree.remove_last_rank();
                                    }
                                },
                            }
                        }
                        // done:
                        if n_num_at_to1_success == 0
                            && tpos1 < cur_tree.get_pos()
                            && cur_tree.is_last_rank(at_rank_canon1) == 1
                        {
                            cur_tree.remove_if_last_atom(at_to1 as i32);
                            cur_tree.remove_last_rank_if_no_atoms();
                        }
                    }
                    if tpos1 < cur_tree.get_pos() && cur_tree.is_last_rank(at_rank_canon1) == 1 {
                        cur_tree.remove_last_rank();
                    } else if tpos1 != cur_tree.get_pos() {
                        return CT_STEREOCOUNT_ERR;
                    }
                    if n_tot_success == 0 || stereo_center_parity == sb_parity_calc {
                        pc = Pc::RepeatAll;
                        continue 'sm;
                    }
                    break 'sm;
                }
            }
        }
    } else {
        // all stereogenic atoms and bonds have been mapped
        if cs.b_stereo_is_better != 0 || cs.b_first_ct != 0 {
            let ret = break_all_ties(
                num_at_tg,
                num_max,
                st1,
                neigh_list,
                n_temp_rank,
                &mut cs.counters,
            );
            if returned_error(ret) {
                return ret;
            }
            let ret = break_all_ties(
                num_at_tg,
                num_max,
                st2,
                neigh_list,
                n_temp_rank,
                &mut cs.counters,
            );
            if returned_error(ret) {
                return ret;
            }
            // the stacks advance by two: the arrays with all ties broken
            cs.n_prev_atom_number[..num_at_tg].fill(0);
            for i1 in 0..num_at_tg {
                let n1 = st1[3][i1] as usize;
                let rc = n_canon_rank_from[n1];
                let n2 = st2[3][i1] as usize;
                n_canon_rank_to[n2] = rc;
                cs.n_prev_atom_number[rc as usize - 1] = n2 as u16;
                n_symm_stereo[i1] = i1 as u16;
                if st1[2][n1] != st2[2][n2] || n_symm_rank[n1] != n_symm_rank[n2] {
                    return CT_STEREO_CANON_ERR;
                }
            }
            cs.l_num_tot_ct += 1;
            cs.l_num_equal_ct = 1;
            cs.l_num_decreased_ct += 1;
            cs.b_stereo_is_better = 0;
            n_tot_success = 1;
            cs.b_first_ct = 0;
            if cs.n_mode & CMODE_REDNDNT_STEREO == 0 {
                let prev = cs.n_prev_atom_number.clone();
                cs.rm.l_num_neigh_list_iter = cs.counters.l_num_neigh_list_iter;
                let (a, b) = (&mut st1[2..], &mut st2[2..]);
                let mut i1 = remove_calculated_non_stereo(
                    at,
                    num_atoms,
                    num_at_tg,
                    a,
                    b,
                    n_temp_rank,
                    neigh_list,
                    n_symm_rank,
                    n_canon_rank_to,
                    &prev,
                    &mut cs.rm,
                    v_ab_parity_unknown,
                );
                cs.counters.l_num_neigh_list_iter = cs.rm.l_num_neigh_list_iter;
                if returned_error(i1) {
                    return i1;
                }
                if i1 < 0 {
                    i1 = -(1 + i1);
                }
                if i1 > 0 {
                    return 4;
                }
            }
        } else {
            // the current stereo name is the same as the previous one
            let mut num_changes = 0;
            let r_max = num_at_tg as u16;
            for r in 1..=r_max {
                if let Some(n1) = b_unique_at_nbr_from_mapping_rank(st1, r) {
                    if let Some(n2) = b_unique_at_nbr_from_mapping_rank(st2, r) {
                        let cr = n_canon_rank_from[n1 as usize] as usize - 1;
                        num_changes += n_join2_mcrs(n_symm_stereo, cs.n_prev_atom_number[cr], n2);
                    } else {
                        return CT_MAPCOUNT_ERR;
                    }
                }
            }
            if num_changes != 0 {
                let mut r = r_max - 1;
                while r != 0 {
                    n_get_mcr(n_symm_stereo, r);
                    r -= 1;
                }
            }
            cs.l_num_equal_ct += 1;
            cs.l_num_tot_ct += 1;
            n_tot_success = 1;
        }
        if cs.time_is_over {
            return CT_TIMEOUT_ERR;
        }
    }
    if n_tot_success == 0 && n_num_mapped_atoms < cs.rm.linear_ct_stereo_carb.len() {
        cs.rm.linear_ct_stereo_carb[n_num_mapped_atoms] = prev_atom;
        cur_tree.set_pos(tpos1);
    }
    n_tot_success
}

/// `map_stereo_bonds4`: maps the stereo bonds and cumulenes (and then the
/// stereocenters) to find the smallest stereo part of the canonical CT
#[allow(clippy::too_many_arguments, unused_assignments)]
pub fn map_stereo_bonds4(
    at: &mut [SpAtom],
    num_atoms: usize,
    num_at_tg: usize,
    num_max: usize,
    mut b_allene: bool,
    n_canon_rank_from: &[u16],
    n_atom_number_canon_from: &[u16],
    n_canon_rank_to: &mut [u16],
    n_symm_rank: &[u16],
    st1: &mut Vec<Vec<u16>>,
    st2: &mut Vec<Vec<u16>>,
    n_temp_rank: &mut [u16],
    n_num_mapped_ranks_input: i32,
    n_symm_stereo: &mut [u16],
    neigh_list: &mut NeighList,
    cs: &mut StereoMapCs,
    cur_tree: &mut CurTree,
    n_num_mapped_bonds: usize,
    v_ab_parity_unknown: i32,
) -> i32 {
    let mut n_tot_success = 0i32;
    let mut prev_bond = AtStereoDble::default();
    let mut tpos1 = cur_tree.get_pos();

    'total_restart: loop {
        if n_num_mapped_bonds == 0 {
            cs.b_rank_used_for_stereo[..num_atoms].fill(0);
            set_use_atom_for_stereo(&mut cs.b_atom_used_for_stereo, at, num_atoms);
            if cs.b_first_ct != 0 && cs.b_keep_symm_rank == 0 {
                for (i, v) in n_symm_stereo[..num_at_tg].iter_mut().enumerate() {
                    *v = i as u16;
                }
            }
            if cs.b_first_ct != 0 {
                tpos1 = 0;
                cur_tree.set_pos(0);
            }
        }

        if n_num_mapped_bonds < cs.rm.linear_ct_stereo_dble.len() {
            let mut b_first_canon_rank = true;
            let mut canon_min1 = 0u16;
            let mut canon_min2 = 0u16;
            let mut prev_bond2 = AtStereoDble::default();
            prev_bond = cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds];
            let (mut at_rank_canon1, mut at_rank_canon2) = if n_num_mapped_bonds != 0 {
                let d = cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds - 1];
                (d.at_num1, d.at_num2)
            } else {
                (0, 0)
            };
            let (mut n_num_choices, mut n_num_unkn, mut n_num_undf, mut n_num_best) =
                (0i32, 0, 0, 0);
            let (mut n_num_worse, mut n_num_calc) = (0i32, 0i32);
            let mut sb_parity_calc = 0i32;
            // (these arrays live across the iterations, as in the C code)
            let mut n_stack_ptr = [0usize; SB_DEPTH];
            let mut n_mapped_ranks = [0i32; SB_DEPTH];
            let mut stereo_bond_parity = 0i32;
            let mut prev_stereo_bond_parity = 0i32;
            let mut pass = 0i32;
            let (mut at_from1, mut at_from2) = (0usize, 0usize);
            let (mut at_rank1, mut at_rank2) = (0usize, 0usize);
            let mut i_max = 0i32;

            #[derive(Clone, Copy, PartialEq)]
            enum Pc {
                Bypass,
                NextCanonRanks,
                RepeatAll,
            }
            let mut pc = Pc::Bypass;
            'sm: loop {
                match pc {
                    Pc::NextCanonRanks => {
                        let cur = cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds];
                        if (cs.b_stereo_is_better == 0 && at_rank_canon1 > cur.at_num1)
                            || (at_rank_canon1 == cur.at_num1 && at_rank_canon2 >= cur.at_num2)
                        {
                            if n_tot_success == 0 {
                                cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds] = prev_bond;
                            }
                            cur_tree.set_pos(tpos1);
                            return n_tot_success;
                        }
                        pc = Pc::Bypass;
                        continue 'sm;
                    }
                    Pc::Bypass => {
                        cur_tree.set_pos(tpos1);
                        let mut ret1 = next_sb_at_canon_ranks2(
                            &mut at_rank_canon1,
                            &mut at_rank_canon2,
                            &mut canon_min1,
                            &mut canon_min2,
                            &mut b_first_canon_rank,
                            &cs.b_atom_used_for_stereo,
                            st1,
                            st2,
                            n_canon_rank_from,
                            n_atom_number_canon_from,
                            at,
                            num_atoms as i32,
                            b_allene,
                        );
                        if !ret1 && !b_allene && b_first_canon_rank {
                            let (mut a1, mut a2, mut m1, mut m2) = (0u16, 0u16, 0u16, 0u16);
                            ret1 = next_sb_at_canon_ranks2(
                                &mut a1,
                                &mut a2,
                                &mut m1,
                                &mut m2,
                                &mut b_first_canon_rank,
                                &cs.b_atom_used_for_stereo,
                                st1,
                                st2,
                                n_canon_rank_from,
                                n_atom_number_canon_from,
                                at,
                                num_atoms as i32,
                                true,
                            );
                            if ret1 {
                                at_rank_canon1 = a1;
                                at_rank_canon2 = a2;
                                canon_min1 = m1;
                                canon_min2 = m2;
                                b_allene = true;
                            }
                        }
                        let cur = cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds];
                        if !ret1
                            || (cs.b_stereo_is_better == 0
                                && (at_rank_canon1 > cur.at_num1
                                    || (at_rank_canon1 == cur.at_num1
                                        && at_rank_canon2 > cur.at_num2)))
                        {
                            if n_tot_success == 0 {
                                cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds] = prev_bond;
                            }
                            return n_tot_success;
                        }
                        n_num_choices = 0;
                        n_num_unkn = 0;
                        n_num_undf = 0;
                        n_num_best = 0;
                        n_num_worse = 0;
                        n_num_calc = 0;
                        pass = 0;
                        prev_stereo_bond_parity = 0;
                        at_from1 = n_atom_number_canon_from[at_rank_canon1 as usize - 1] as usize;
                        at_rank1 = st1[0][at_from1] as usize;
                        at_from2 = n_atom_number_canon_from[at_rank_canon2 as usize - 1] as usize;
                        at_rank2 = st1[0][at_from2] as usize;
                        i_max = at_rank1 as i32 - 1;
                        if at_rank1 as u16 != st2[0][st2[1][at_rank1 - 1] as usize]
                            || at_rank2 as u16 != st2[0][st2[1][at_rank2 - 1] as usize]
                        {
                            return CT_STEREOCOUNT_ERR;
                        }
                        for j1 in 0..=i_max {
                            let at_to1 = st2[1][(i_max - j1) as usize] as usize;
                            if at_rank1 as u16 != st2[0][at_to1] {
                                break;
                            }
                            if at[at_to1].stereo_bond_neighbor[0] == 0 {
                                continue;
                            }
                            for j2 in 0..MAX_NUM_STEREO_BONDS {
                                let s = at[at_to1].stereo_bond_neighbor[j2];
                                if s == 0 {
                                    break;
                                }
                                let at_to2 = s as usize - 1;
                                if st1[0][at_from2] != st2[0][at_to2] {
                                    continue;
                                }
                                stereo_bond_parity =
                                    parity_val(i32::from(at[at_to1].stereo_bond_parity[j2]));
                                let mut i = 0;
                                match stereo_bond_parity {
                                    AB_PARITY_UNDF => n_num_undf += 1,
                                    AB_PARITY_UNKN => n_num_unkn += 1,
                                    BEST_PARITY => n_num_best += 1,
                                    WORSE_PARITY => n_num_worse += 1,
                                    AB_PARITY_CALC => n_num_calc += 1,
                                    AB_PARITY_NONE => i += 1,
                                    _ => {}
                                }
                                n_num_choices += i32::from(i == 0);
                            }
                        }
                        if n_num_choices
                            != n_num_calc + n_num_undf + n_num_unkn + n_num_best + n_num_worse
                        {
                            return CT_STEREOCOUNT_ERR;
                        }
                        if n_num_choices == 0 {
                            pc = Pc::NextCanonRanks;
                            continue 'sm;
                        }
                        sb_parity_calc = if n_num_calc > 0 { BEST_PARITY } else { 0 };
                        pc = Pc::RepeatAll;
                        continue 'sm;
                    }
                    Pc::RepeatAll => {
                        if n_num_mapped_bonds == 0 {
                            cs.b_stereo_is_better = 0;
                        }
                        let first = pass == 0;
                        pass += 1;
                        if first {
                            stereo_bond_parity = if sb_parity_calc != 0 {
                                BEST_PARITY
                            } else if n_num_best != 0 {
                                BEST_PARITY
                            } else if n_num_worse != 0 {
                                WORSE_PARITY
                            } else if n_num_unkn != 0 {
                                AB_PARITY_UNKN
                            } else if n_num_undf != 0 {
                                AB_PARITY_UNDF
                            } else {
                                AB_PARITY_NONE
                            };
                        } else {
                            prev_stereo_bond_parity = stereo_bond_parity;
                            let i = next_stereo_parity2_test(
                                &mut stereo_bond_parity,
                                &mut sb_parity_calc,
                                n_num_best,
                                n_num_worse,
                                n_num_unkn,
                                n_num_undf,
                                n_num_calc,
                                v_ab_parity_unknown,
                            );
                            match i {
                                0 => {}
                                1 => {
                                    pc = Pc::NextCanonRanks;
                                    continue 'sm;
                                }
                                _ => return i,
                            }
                        }
                        if stereo_bond_parity == AB_PARITY_NONE {
                            return CT_STEREOCOUNT_ERR;
                        }
                        if cs.b_stereo_is_better == 0 {
                            let cur: &AtStereoDble = if n_tot_success != 0 {
                                &cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds]
                            } else {
                                &prev_bond
                            };
                            let c = compare_lin_ct_stereo_double_to_values(
                                cur,
                                at_rank_canon1,
                                at_rank_canon2,
                                stereo_bond_parity as u8,
                            );
                            if c < 0 {
                                if n_tot_success == 0 {
                                    cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds] = prev_bond;
                                }
                                cur_tree.set_pos(tpos1);
                                return n_tot_success;
                            }
                        }
                        let mut b_all_parities_identical = false;
                        let mut last_mapped_to1: i32 = -1;
                        let mut b_stereo_is_better_was_set_here = false;
                        if n_num_mapped_bonds == 0 && prev_stereo_bond_parity != stereo_bond_parity
                        {
                            cs.b_stereo_is_better = 0;
                        }
                        if prev_stereo_bond_parity != stereo_bond_parity {
                            cur_tree.set_pos(tpos1);
                        }
                        for j1 in 0..=i_max {
                            let at_to1 = st2[1][(i_max - j1) as usize] as usize;
                            if at_rank1 as u16 != st2[0][at_to1] {
                                break;
                            }
                            let mut n_num_at_to1_success = 0i32;
                            if at[at_to1].stereo_bond_neighbor[0] == 0 {
                                continue;
                            }
                            if tpos1 < cur_tree.get_pos()
                                && cur_tree.is_last_rank(at_rank_canon1) == 1
                                && cur_tree.is_last_atom_equ(at_to1 as i32, n_symm_stereo) == 1
                            {
                                continue;
                            }
                            let mut b_all_parities_identical2 = false;
                            for j2 in 0..MAX_NUM_STEREO_BONDS {
                                let s = at[at_to1].stereo_bond_neighbor[j2];
                                if s == 0 {
                                    break;
                                }
                                let at_to2 = s as usize - 1;
                                let mut en1 = [EqNeigh::default(); 2];
                                let mut en2 = [EqNeigh::default(); 2];
                                if st1[0][at_from2] != st2[0][at_to2] {
                                    continue;
                                }
                                let sbp = i32::from(at[at_to1].stereo_bond_parity[j2]);
                                if parity_known(sbp) {
                                    if stereo_bond_parity == sb_parity_calc {
                                        continue;
                                    }
                                    if stereo_bond_parity != parity_val(sbp) {
                                        continue;
                                    }
                                } else if parity_val(sbp) == AB_PARITY_CALC {
                                    if stereo_bond_parity != sb_parity_calc {
                                        continue;
                                    }
                                } else {
                                    return CT_STEREOCOUNT_ERR;
                                }
                                let mut istk = 0usize;
                                n_mapped_ranks[0] = n_num_mapped_ranks_input;
                                let mut b_add_stack = 0i32;
                                b_all_parities_identical =
                                    (sbp & KNOWN_PARITIES_EQL) != 0 && parity_known(sbp);
                                if !b_all_parities_identical
                                    && n_num_calc == 0
                                    && (i32::from(n_num_undf == 0)
                                        + i32::from(n_num_unkn == 0)
                                        + i32::from(n_num_best == 0)
                                        + i32::from(n_num_worse == 0))
                                        == 3
                                {
                                    let r = all_sb_same(
                                        at_rank_canon1,
                                        at_rank_canon2,
                                        st1,
                                        st2,
                                        n_atom_number_canon_from,
                                        at,
                                    );
                                    if r < 0 {
                                        return CT_STEREOCOUNT_ERR;
                                    }
                                    b_all_parities_identical = r != 0;
                                }
                                if !b_all_parities_identical {
                                    if last_mapped_to1 != at_to1 as i32 {
                                        let mut out = 0;
                                        let ret1 = map_at(
                                            num_at_tg,
                                            num_max,
                                            at_from1,
                                            at_to1,
                                            n_temp_rank,
                                            n_mapped_ranks[istk],
                                            &mut out,
                                            cs,
                                            neigh_list,
                                            st1,
                                            st2,
                                            n_stack_ptr[istk],
                                            &mut b_add_stack,
                                        );
                                        n_mapped_ranks[istk + 1] = out;
                                        if returned_error(ret1) {
                                            return ret1;
                                        }
                                        n_stack_ptr[istk + 1] =
                                            n_stack_ptr[istk] + b_add_stack as usize;
                                        last_mapped_to1 = at_to1 as i32;
                                        if b_add_stack != 0 {
                                            if tpos1 == cur_tree.get_pos()
                                                || cur_tree.is_last_rank(at_rank_canon1) == 0
                                            {
                                                cur_tree.add_rank(at_rank_canon1);
                                            }
                                            cur_tree.add_atom(at_to1 as i32);
                                        }
                                    }
                                    istk += 1;
                                    if st1[n_stack_ptr[istk]][at_from2]
                                        != st2[n_stack_ptr[istk]][at_to2]
                                    {
                                        if sb_parity_calc != stereo_bond_parity {
                                            n_num_choices -= 1;
                                            n_num_undf -=
                                                i32::from(stereo_bond_parity == AB_PARITY_UNDF);
                                            n_num_unkn -=
                                                i32::from(stereo_bond_parity == AB_PARITY_UNKN);
                                            n_num_best -=
                                                i32::from(stereo_bond_parity == BEST_PARITY);
                                            n_num_worse -=
                                                i32::from(stereo_bond_parity == WORSE_PARITY);
                                        } else if sb_parity_calc == BEST_PARITY {
                                            n_num_choices -= 1;
                                            n_num_calc -= 1;
                                        }
                                        continue;
                                    }
                                    if n_stack_ptr[istk] > n_stack_ptr[istk - 1] {
                                        let off = n_stack_ptr[istk];
                                        let r = all_sb_same(
                                            at_rank_canon1,
                                            at_rank_canon2,
                                            &st1[off..],
                                            &st2[off..],
                                            n_atom_number_canon_from,
                                            at,
                                        );
                                        if r < 0 {
                                            return CT_STEREOBOND_ERROR;
                                        }
                                        b_all_parities_identical2 = r != 0;
                                    } else {
                                        b_all_parities_identical2 = false;
                                    }
                                    if b_all_parities_identical2 {
                                        let off = n_stack_ptr[istk] + 2;
                                        if let Some(t) = st1.get_mut(off..) {
                                            clear_previous_mappings(t);
                                        }
                                    } else {
                                        if tpos1 < cur_tree.get_pos()
                                            && cur_tree.is_last_rank(at_rank_canon2) == 1
                                            && cur_tree
                                                .is_last_atom_equ(at_to2 as i32, n_symm_stereo)
                                                == 1
                                        {
                                            continue;
                                        }
                                        let mut out = 0;
                                        let ret2 = map_at(
                                            num_at_tg,
                                            num_max,
                                            at_from2,
                                            at_to2,
                                            n_temp_rank,
                                            n_mapped_ranks[istk],
                                            &mut out,
                                            cs,
                                            neigh_list,
                                            st1,
                                            st2,
                                            n_stack_ptr[istk],
                                            &mut b_add_stack,
                                        );
                                        n_mapped_ranks[istk + 1] = out;
                                        if returned_error(ret2) {
                                            return ret2;
                                        }
                                        n_stack_ptr[istk + 1] =
                                            n_stack_ptr[istk] + b_add_stack as usize;
                                        istk += 1;
                                        if b_add_stack != 0 {
                                            if tpos1 == cur_tree.get_pos()
                                                || cur_tree.is_last_rank(at_rank_canon2) == 0
                                            {
                                                cur_tree.add_rank(at_rank_canon2);
                                            }
                                            cur_tree.add_atom(at_to2 as i32);
                                        }
                                    }
                                } else if let Some(t) = st1.get_mut(2..) {
                                    clear_previous_mappings(t);
                                }

                                // known bond parity: do not map the stereo bond neighbors
                                if stereo_bond_parity != sb_parity_calc {
                                    let bond_parity = stereo_bond_parity;
                                    let c = compare_lin_ct_stereo_double_to_values(
                                        &cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds],
                                        at_rank_canon1,
                                        at_rank_canon2,
                                        bond_parity as u8,
                                    );
                                    if c < 0 && cs.b_stereo_is_better == 0 {
                                        cs.l_num_rejected_ct += 1;
                                        if tpos1 < cur_tree.get_pos()
                                            && cur_tree.is_last_rank(at_rank_canon2) == 1
                                        {
                                            cur_tree.remove_if_last_atom(at_to2 as i32);
                                            cur_tree.remove_last_rank_if_no_atoms();
                                        }
                                        continue;
                                    }
                                    if c > 0 && cs.b_stereo_is_better == 0 {
                                        cs.b_stereo_is_better = 1;
                                        b_stereo_is_better_was_set_here = true;
                                        prev_bond2 =
                                            cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds];
                                    }
                                    {
                                        let d =
                                            &mut cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds];
                                        d.at_num1 = at_rank_canon1;
                                        d.at_num2 = at_rank_canon2;
                                        d.parity = bond_parity as u8;
                                    }
                                    cs.b_rank_used_for_stereo[at_from1] += 1;
                                    cs.b_rank_used_for_stereo[at_from2] += 1;
                                    if !b_all_parities_identical {
                                        cs.b_atom_used_for_stereo[at_to1] -= 1;
                                        cs.b_atom_used_for_stereo[at_to2] -= 1;
                                    }
                                    let ret2 = {
                                        let off = n_stack_ptr[istk];
                                        let mut t1 = st1.split_off(off);
                                        let mut t2 = st2.split_off(off);
                                        let r = map_stereo_bonds4(
                                            at,
                                            num_atoms,
                                            num_at_tg,
                                            num_max,
                                            b_allene,
                                            n_canon_rank_from,
                                            n_atom_number_canon_from,
                                            n_canon_rank_to,
                                            n_symm_rank,
                                            &mut t1,
                                            &mut t2,
                                            n_temp_rank,
                                            n_mapped_ranks[istk],
                                            n_symm_stereo,
                                            neigh_list,
                                            cs,
                                            cur_tree,
                                            n_num_mapped_bonds + 1,
                                            v_ab_parity_unknown,
                                        );
                                        st1.append(&mut t1);
                                        st2.append(&mut t2);
                                        r
                                    };
                                    if !b_all_parities_identical {
                                        cs.b_atom_used_for_stereo[at_to1] += 1;
                                        cs.b_atom_used_for_stereo[at_to2] += 1;
                                    }
                                    cs.b_rank_used_for_stereo[at_from1] -= 1;
                                    cs.b_rank_used_for_stereo[at_from2] -= 1;
                                    if ret2 == 4 {
                                        if n_num_mapped_bonds != 0 {
                                            return ret2;
                                        }
                                        cs.b_first_ct = 1;
                                        continue 'total_restart;
                                    }
                                    if returned_error(ret2) {
                                        return ret2;
                                    }
                                    if ret2 > 0 {
                                        n_tot_success |= 1;
                                        n_num_at_to1_success += 1;
                                        if b_stereo_is_better_was_set_here || (ret2 & 2) != 0 {
                                            cur_tree.keep_last_atoms_only(tpos1, 1);
                                            n_tot_success |= 2;
                                        }
                                    } else {
                                        if b_stereo_is_better_was_set_here {
                                            cs.b_stereo_is_better = 0;
                                            cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds] =
                                                prev_bond2;
                                        }
                                        if tpos1 < cur_tree.get_pos()
                                            && cur_tree.is_last_rank(at_rank_canon2) == 1
                                        {
                                            cur_tree.remove_if_last_atom(at_to2 as i32);
                                            cur_tree.remove_last_rank_if_no_atoms();
                                        }
                                    }
                                    b_stereo_is_better_was_set_here = false;
                                    if b_all_parities_identical || b_all_parities_identical2 {
                                        break;
                                    }
                                    continue;
                                }

                                // unknown bond parity: may need to map the stereo bond neighbors
                                let mut parity1 = parity_of_mapped_half_bond(
                                    at_from1,
                                    at_to1,
                                    at_from2,
                                    at_to2,
                                    at,
                                    Some(&mut en1[0]),
                                    n_canon_rank_from,
                                    &st1[n_stack_ptr[istk]],
                                    &st2[n_stack_ptr[istk]],
                                );
                                if parity1 == 0 {
                                    return CT_STEREOCOUNT_ERR;
                                }
                                let num1 = if parity1 > 0 { 1 } else { 2 };
                                let mut at_rank_canon_n1 = 0u16;
                                for i in 0..num1 {
                                    let mut at_no_n1_num_success = 0;
                                    let mut at_to_n1 = 0usize;
                                    let mut istk2 = istk;
                                    if num1 == 2 {
                                        at_rank_canon_n1 =
                                            n_canon_rank_from[en1[0].from_at as usize];
                                        let at_from_n1 = en1[0].from_at as usize;
                                        at_to_n1 = en1[0].to_at[i] as usize;
                                        if tpos1 < cur_tree.get_pos()
                                            && cur_tree.is_last_rank(at_rank_canon_n1) == 1
                                            && cur_tree
                                                .is_last_atom_equ(at_to_n1 as i32, n_symm_stereo)
                                                == 1
                                        {
                                            continue;
                                        }
                                        let mut out = 0;
                                        let ret1 = map_at(
                                            num_at_tg,
                                            num_max,
                                            at_from_n1,
                                            at_to_n1,
                                            n_temp_rank,
                                            n_mapped_ranks[istk2],
                                            &mut out,
                                            cs,
                                            neigh_list,
                                            st1,
                                            st2,
                                            n_stack_ptr[istk2],
                                            &mut b_add_stack,
                                        );
                                        n_mapped_ranks[istk2 + 1] = out;
                                        if returned_error(ret1) {
                                            return ret1;
                                        }
                                        n_stack_ptr[istk2 + 1] =
                                            n_stack_ptr[istk2] + b_add_stack as usize;
                                        istk2 += 1;
                                        if istk2 >= SB_DEPTH {
                                            return CT_OVERFLOW;
                                        }
                                        if b_add_stack != 0 {
                                            if tpos1 == cur_tree.get_pos()
                                                || cur_tree.is_last_rank(at_rank_canon_n1) == 0
                                            {
                                                cur_tree.add_rank(at_rank_canon_n1);
                                            }
                                            cur_tree.add_atom(at_to_n1 as i32);
                                        }
                                        parity1 = parity_of_mapped_half_bond(
                                            at_from1,
                                            at_to1,
                                            at_from2,
                                            at_to2,
                                            at,
                                            Some(&mut en1[1]),
                                            n_canon_rank_from,
                                            &st1[n_stack_ptr[istk2]],
                                            &st2[n_stack_ptr[istk2]],
                                        );
                                        if parity1 <= 0 {
                                            return CT_STEREOCOUNT_ERR;
                                        }
                                    } else {
                                        n_mapped_ranks[istk2 + 1] = n_mapped_ranks[istk2];
                                        n_stack_ptr[istk2 + 1] = n_stack_ptr[istk2];
                                        istk2 += 1;
                                    }

                                    let mut parity2 = parity_of_mapped_half_bond(
                                        at_from2,
                                        at_to2,
                                        at_from1,
                                        at_to1,
                                        at,
                                        Some(&mut en2[0]),
                                        n_canon_rank_from,
                                        &st1[n_stack_ptr[istk2]],
                                        &st2[n_stack_ptr[istk2]],
                                    );
                                    if parity2 == 0 {
                                        return CT_STEREOCOUNT_ERR;
                                    }
                                    let num2 = if parity2 > 0 { 1 } else { 2 };
                                    let mut at_rank_canon_n2 = 0u16;
                                    for j in 0..num2 {
                                        let mut at_to_n2 = 0usize;
                                        let mut istk3 = istk2;
                                        if num2 == 2 {
                                            at_rank_canon_n2 =
                                                n_canon_rank_from[en2[0].from_at as usize];
                                            let at_from_n2 = en2[0].from_at as usize;
                                            at_to_n2 = en2[0].to_at[j] as usize;
                                            if tpos1 < cur_tree.get_pos()
                                                && cur_tree.is_last_rank(at_rank_canon_n2) == 1
                                                && cur_tree.is_last_atom_equ(
                                                    at_to_n2 as i32,
                                                    n_symm_stereo,
                                                ) == 1
                                            {
                                                continue;
                                            }
                                            let mut out = 0;
                                            let ret1 = map_at(
                                                num_at_tg,
                                                num_max,
                                                at_from_n2,
                                                at_to_n2,
                                                n_temp_rank,
                                                n_mapped_ranks[istk3],
                                                &mut out,
                                                cs,
                                                neigh_list,
                                                st1,
                                                st2,
                                                n_stack_ptr[istk3],
                                                &mut b_add_stack,
                                            );
                                            n_mapped_ranks[istk3 + 1] = out;
                                            if returned_error(ret1) {
                                                return ret1;
                                            }
                                            n_stack_ptr[istk3 + 1] =
                                                n_stack_ptr[istk3] + b_add_stack as usize;
                                            istk3 += 1;
                                            if b_add_stack != 0 {
                                                if tpos1 == cur_tree.get_pos()
                                                    || cur_tree.is_last_rank(at_rank_canon_n2) == 0
                                                {
                                                    cur_tree.add_rank(at_rank_canon_n2);
                                                }
                                                cur_tree.add_atom(at_to_n2 as i32);
                                            }
                                            parity2 = parity_of_mapped_half_bond(
                                                at_from2,
                                                at_to2,
                                                at_from1,
                                                at_to1,
                                                at,
                                                Some(&mut en2[1]),
                                                n_canon_rank_from,
                                                &st1[n_stack_ptr[istk3]],
                                                &st2[n_stack_ptr[istk3]],
                                            );
                                            if parity2 <= 0 {
                                                return CT_STEREOCOUNT_ERR;
                                            }
                                        } else {
                                            n_mapped_ranks[istk3 + 1] = n_mapped_ranks[istk3];
                                            n_stack_ptr[istk3 + 1] = n_stack_ptr[istk3];
                                            istk3 += 1;
                                        }
                                        if parity1 <= 0 || parity2 <= 0 {
                                            return CT_STEREOCOUNT_ERR;
                                        }
                                        let mut bond_parity = if atom_parity_well_def(parity1)
                                            && atom_parity_well_def(parity2)
                                        {
                                            2 - (parity1 + parity2) % 2
                                        } else {
                                            parity1.max(parity2)
                                        };
                                        if atom_parity_well_def(bond_parity)
                                            && at[at_to1].stereo_bond_z_prod[j2] < 0
                                        {
                                            bond_parity = 2 - (bond_parity + 1) % 2;
                                        }
                                        let c = compare_lin_ct_stereo_double_to_values(
                                            &cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds],
                                            at_rank_canon1,
                                            at_rank_canon2,
                                            bond_parity as u8,
                                        );
                                        if sb_parity_calc != bond_parity
                                            || (c < 0 && cs.b_stereo_is_better == 0)
                                        {
                                            cs.l_num_rejected_ct += 1;
                                            if tpos1 < cur_tree.get_pos()
                                                && cur_tree.is_last_rank(at_rank_canon_n2) == 1
                                            {
                                                cur_tree.remove_if_last_atom(at_to_n2 as i32);
                                                cur_tree.remove_last_rank_if_no_atoms();
                                            }
                                            continue;
                                        }
                                        if c > 0 && cs.b_stereo_is_better == 0 {
                                            cs.b_stereo_is_better = 1;
                                            b_stereo_is_better_was_set_here = true;
                                            prev_bond2 =
                                                cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds];
                                        }
                                        {
                                            let d = &mut cs.rm.linear_ct_stereo_dble
                                                [n_num_mapped_bonds];
                                            d.at_num1 = at_rank_canon1;
                                            d.at_num2 = at_rank_canon2;
                                            d.parity = bond_parity as u8;
                                        }
                                        cs.b_rank_used_for_stereo[at_from1] += 1;
                                        cs.b_rank_used_for_stereo[at_from2] += 1;
                                        cs.b_atom_used_for_stereo[at_to1] -= 1;
                                        cs.b_atom_used_for_stereo[at_to2] -= 1;
                                        let ret2 = {
                                            let off = n_stack_ptr[istk3];
                                            let mut t1 = st1.split_off(off);
                                            let mut t2 = st2.split_off(off);
                                            let r = map_stereo_bonds4(
                                                at,
                                                num_atoms,
                                                num_at_tg,
                                                num_max,
                                                b_allene,
                                                n_canon_rank_from,
                                                n_atom_number_canon_from,
                                                n_canon_rank_to,
                                                n_symm_rank,
                                                &mut t1,
                                                &mut t2,
                                                n_temp_rank,
                                                n_mapped_ranks[istk3],
                                                n_symm_stereo,
                                                neigh_list,
                                                cs,
                                                cur_tree,
                                                n_num_mapped_bonds + 1,
                                                v_ab_parity_unknown,
                                            );
                                            st1.append(&mut t1);
                                            st2.append(&mut t2);
                                            r
                                        };
                                        cs.b_rank_used_for_stereo[at_from1] -= 1;
                                        cs.b_rank_used_for_stereo[at_from2] -= 1;
                                        cs.b_atom_used_for_stereo[at_to1] += 1;
                                        cs.b_atom_used_for_stereo[at_to2] += 1;
                                        if ret2 == 4 {
                                            if n_num_mapped_bonds != 0 {
                                                return ret2;
                                            }
                                            cs.b_first_ct = 1;
                                            continue 'total_restart;
                                        }
                                        if returned_error(ret2) {
                                            return ret2;
                                        }
                                        if ret2 > 0 {
                                            n_tot_success |= 1;
                                            n_num_at_to1_success += 1;
                                            if b_stereo_is_better_was_set_here || (ret2 & 2) != 0 {
                                                cur_tree.keep_last_atoms_only(tpos1, 1);
                                                n_tot_success |= 2;
                                            }
                                            at_no_n1_num_success += 1;
                                        } else {
                                            if b_stereo_is_better_was_set_here {
                                                cs.b_stereo_is_better = 0;
                                                cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds] =
                                                    prev_bond2;
                                            }
                                            if tpos1 < cur_tree.get_pos()
                                                && cur_tree.is_last_rank(at_rank_canon_n2) == 1
                                            {
                                                cur_tree.remove_if_last_atom(at_to_n2 as i32);
                                                cur_tree.remove_last_rank_if_no_atoms();
                                            }
                                        }
                                        b_stereo_is_better_was_set_here = false;
                                    }
                                    if tpos1 < cur_tree.get_pos()
                                        && cur_tree.is_last_rank(at_rank_canon_n2) == 1
                                    {
                                        cur_tree.remove_last_rank();
                                    }
                                    if at_no_n1_num_success == 0
                                        && tpos1 < cur_tree.get_pos()
                                        && cur_tree.is_last_rank(at_rank_canon_n1) == 1
                                    {
                                        cur_tree.remove_if_last_atom(at_to_n1 as i32);
                                    }
                                }
                                if tpos1 < cur_tree.get_pos()
                                    && cur_tree.is_last_rank(at_rank_canon_n1) == 1
                                {
                                    cur_tree.remove_last_rank();
                                }
                            }
                            if tpos1 < cur_tree.get_pos()
                                && cur_tree.is_last_rank(at_rank_canon2) == 1
                            {
                                cur_tree.remove_last_rank();
                            }
                            if n_num_at_to1_success == 0
                                && tpos1 < cur_tree.get_pos()
                                && cur_tree.is_last_rank(at_rank_canon1) == 1
                            {
                                cur_tree.remove_if_last_atom(at_to1 as i32);
                                cur_tree.remove_last_rank_if_no_atoms();
                            }
                            if b_all_parities_identical {
                                break;
                            }
                        }
                        if tpos1 < cur_tree.get_pos() && cur_tree.is_last_rank(at_rank_canon1) == 1
                        {
                            cur_tree.remove_last_rank();
                        } else if tpos1 != cur_tree.get_pos() {
                            return CT_STEREOCOUNT_ERR;
                        }
                        if n_tot_success == 0 || stereo_bond_parity == sb_parity_calc {
                            pc = Pc::RepeatAll;
                            continue 'sm;
                        }
                        break 'sm;
                    }
                }
            }
        } else {
            if n_num_mapped_bonds == 0 {
                cs.b_stereo_is_better = 0;
            }
            if n_num_mapped_bonds < cs.rm.linear_ct_stereo_dble.len() {
                prev_bond = cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds];
            }
            let ret = map_stereo_atoms4(
                at,
                num_atoms,
                num_at_tg,
                num_max,
                n_canon_rank_from,
                n_atom_number_canon_from,
                n_canon_rank_to,
                n_symm_rank,
                st1,
                st2,
                n_temp_rank,
                n_num_mapped_ranks_input,
                n_symm_stereo,
                neigh_list,
                cs,
                cur_tree,
                0,
                v_ab_parity_unknown,
            );
            if ret == 4 {
                if n_num_mapped_bonds != 0 {
                    return ret;
                }
                cs.b_first_ct = 1;
                continue 'total_restart;
            }
            if returned_error(ret) {
                return ret;
            }
            if ret > 0 {
                n_tot_success |= 1;
                if ret & 2 != 0 {
                    cur_tree.keep_last_atoms_only(tpos1, 1);
                    n_tot_success |= 2;
                }
            }
        }
        if n_tot_success == 0
            && !cs.rm.linear_ct_stereo_dble.is_empty()
            && n_num_mapped_bonds < cs.rm.linear_ct_stereo_dble.len()
        {
            cs.rm.linear_ct_stereo_dble[n_num_mapped_bonds] = prev_bond;
        }
        return n_tot_success;
    }
}
