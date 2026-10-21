//! Rank refinement and neighbor lists of the InChI canonicalizer:
//! comparators and list builders of `ichisort.c` (`CompNeighListRanks`,
//! `CreateNeighList`, `BreakAllTies`, ...) and the equitable-partition code of
//! `ichimap2.c` (`SortedEquInfoToRanks`, `SetNewRanksFromNeighLists*`,
//! `DifferentiateRanks*`).
//!
//! A `NEIGH_LIST` is a `Vec<u16>` whose first element is the number of
//! neighbors; the `CANON_GLOBALS` sort state of the C code is passed
//! explicitly to the comparators.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSE/InChI-MIT.txt).

use super::ichisort::{
    compare_neigh_list_lex, compare_neigh_list_lex_up_to_max_rank, inchi_qsort, insertions_sort,
    insertions_sort_at_numbers, insertions_sort_neigh_list_at_numbers,
    insertions_sort_neigh_list_at_numbers3, insertions_sort_neigh_list_by_symm_and_canon_rank,
};
use super::ichister::SpAtom;
use super::ichitaut::{T_NUM_ISOTOPIC, T_NUM_NO_ISOTOPIC, TGroupInfo};

/// `NEIGH_LIST` array: `list[i][0]` is the number of neighbors of vertex `i`
pub type NeighList = Vec<Vec<u16>>;

// ichicant.h
pub const AT_INV_HILL_ORDER: usize = 0;
pub const AT_INV_NUM_CONNECTIONS: usize = 1;
pub const AT_INV_NUM_H: usize = 2;
pub const AT_INV_NUM_TG_ENDPOINTS: usize = 3;
pub const AT_INV_TG_NUMBERS: usize = 4;
pub const AT_INV_NUM_H_FIX: usize = AT_INV_TG_NUMBERS + T_NUM_NO_ISOTOPIC;
pub const AT_INV_BREAK1: usize = AT_INV_NUM_H_FIX + 1;
pub const AT_INV_TAUT_ISO: usize = AT_INV_BREAK1;
pub const AT_INV_LENGTH: usize = AT_INV_TAUT_ISO + T_NUM_ISOTOPIC;

/// `ATOM_INVARIANT2`
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AtomInvariant2 {
    pub val: [u16; AT_INV_LENGTH],
    pub iso_sort_key: i64,
    pub iso_aux_key: i8,
}

/// `CompAtomInvariants2Only`
pub fn comp_atom_invariants2_only(a1: &AtomInvariant2, a2: &AtomInvariant2) -> i32 {
    let mut i = 0;
    while i < AT_INV_BREAK1 {
        if a1.val[i] != a2.val[i] {
            return i32::from(a1.val[i]) - i32::from(a2.val[i]);
        }
        i += 1;
    }
    if a1.iso_sort_key != a2.iso_sort_key {
        return if a1.iso_sort_key > a2.iso_sort_key {
            1
        } else {
            -1
        };
    }
    while i < AT_INV_LENGTH {
        // (sic) the C code skips the *different* values and returns the
        // difference of the first equal pair, which is 0
        if a1.val[i] != a2.val[i] {
            i += 1;
            continue;
        }
        return i32::from(a1.val[i]) - i32::from(a2.val[i]);
    }
    if a1.iso_aux_key != a2.iso_aux_key {
        return if a1.iso_aux_key > a2.iso_aux_key {
            1
        } else {
            -1
        };
    }
    0
}

/// `CompAtomInvariants2` on the atoms `a` and `b` of `inv`
pub fn comp_atom_invariants2(inv: &[AtomInvariant2], a: u16, b: u16) -> i32 {
    let ret = comp_atom_invariants2_only(&inv[a as usize], &inv[b as usize]);
    if ret == 0 {
        i32::from(a) - i32::from(b)
    } else {
        ret
    }
}

/// `CompRank`
pub fn comp_rank(rank: &[u16], a: u16, b: u16) -> i32 {
    i32::from(rank[a as usize]) - i32::from(rank[b as usize])
}

/// `CompRanksOrd`
pub fn comp_ranks_ord(rank: &[u16], a: u16, b: u16) -> i32 {
    let ret = comp_rank(rank, a, b);
    if ret == 0 {
        i32::from(a) - i32::from(b)
    } else {
        ret
    }
}

/// `compare_NeighLists`
pub fn compare_neigh_lists(nl: &NeighList, a: u16, b: u16, rank: &[u16]) -> i32 {
    compare_neigh_list_lex(&nl[a as usize], &nl[b as usize], rank)
}

/// `CompNeighListRanks`
pub fn comp_neigh_list_ranks(nl: &NeighList, rank: &[u16], a: u16, b: u16) -> i32 {
    let ret = comp_rank(rank, a, b);
    if ret == 0 {
        compare_neigh_lists(nl, a, b, rank)
    } else {
        ret
    }
}

/// `CompNeighLists`
pub fn comp_neigh_lists(nl: &NeighList, rank: &[u16], a: u16, b: u16) -> i32 {
    compare_neigh_lists(nl, a, b, rank)
}

/// `CompNeighListsUpToMaxRank`
pub fn comp_neigh_lists_up_to_max_rank(
    nl: &NeighList,
    rank: &[u16],
    max_rank: u16,
    a: u16,
    b: u16,
) -> i32 {
    compare_neigh_list_lex_up_to_max_rank(&nl[a as usize], &nl[b as usize], rank, max_rank)
}

/// `CompNeighListRanksOrd`
pub fn comp_neigh_list_ranks_ord(nl: &NeighList, rank: &[u16], a: u16, b: u16) -> i32 {
    let ret = comp_neigh_list_ranks(nl, rank, a, b);
    if ret == 0 {
        i32::from(a) - i32::from(b)
    } else {
        ret
    }
}

/// `CompRanksInvOrd`
pub fn comp_ranks_inv_ord(a: u16, b: u16) -> i32 {
    i32::from(b) - i32::from(a)
}

/// `CompNeighborsRanksCountEql` (`CT_NEIGH_INCREASE`); counts the equal pairs in `count`
pub fn comp_neighbors_ranks_count_eql(rank: &[u16], count: &mut i32, a: u16, b: u16) -> i32 {
    let ret = i32::from(rank[a as usize]) - i32::from(rank[b as usize]);
    *count += i32::from(ret == 0);
    ret
}

/// `CreateNeighListFromLinearCT`: the neighbor lists of a linear connection
/// table (canonical numbers); `None` if the table is inconsistent
pub fn create_neigh_list_from_linear_ct(
    linear_ct: &[u16],
    n_len_ct: usize,
    num_atoms: usize,
) -> Option<NeighList> {
    if linear_ct[0] as usize > num_atoms {
        return None;
    }
    let mut valence = vec![0i32; num_atoms + 1];
    let mut n_vertex = linear_ct[0];
    for &n_neigh in &linear_ct[1..n_len_ct] {
        if n_neigh < n_vertex {
            valence[n_neigh as usize] += 1;
            valence[n_vertex as usize] += 1;
        } else {
            n_vertex = n_neigh;
            if n_vertex as usize > num_atoms {
                return None;
            }
        }
    }
    if n_vertex as usize != num_atoms {
        return None;
    }
    let mut pp: NeighList = (1..=num_atoms)
        .map(|i| vec![0u16; valence[i] as usize + 1])
        .collect();
    let mut n_vertex = linear_ct[0] as i32 - 1;
    for &ct in &linear_ct[1..n_len_ct] {
        let n_neigh = ct as i32 - 1;
        if n_neigh < n_vertex {
            pp[n_vertex as usize][0] += 1;
            let j = pp[n_vertex as usize][0] as usize;
            pp[n_vertex as usize][j] = n_neigh as u16;
            pp[n_neigh as usize][0] += 1;
            let j = pp[n_neigh as usize][0] as usize;
            pp[n_neigh as usize][j] = n_vertex as u16;
        } else {
            n_vertex = n_neigh;
            if n_vertex as usize >= num_atoms {
                return None;
            }
        }
    }
    Some(pp)
}

/// `CreateNeighList`: neighbor lists of atoms and t-groups (the t-group
/// vertices follow the atoms); with `b_double_bond_square` neighbors on double
/// bonds appear twice
pub fn create_neigh_list(
    num_atoms: usize,
    num_at_tg: usize,
    at: &[SpAtom],
    b_double_bond_square: bool,
    t_group_info: Option<&TGroupInfo>,
) -> NeighList {
    const BOND_DOUBLE: u8 = 2;
    let mut pp: NeighList = vec![Vec::new(); num_at_tg];
    let mut num_t_groups = 0;
    let mut t_groups = &[][..];
    let mut endpoints = &[][..];
    if num_at_tg > num_atoms {
        if let Some(tgi) = t_group_info {
            num_t_groups = tgi.num_t_groups;
            t_groups = tgi.t_group.as_deref().unwrap_or(&[]);
            endpoints = &tgi.n_endpoint_atom_number;
        }
    }
    for i in 0..num_atoms {
        let mut list = vec![0u16];
        for j in 0..at[i].valence as usize {
            list.push(at[i].neighbor[j]);
            if b_double_bond_square && at[i].bond_type[j] == BOND_DOUBLE {
                list.push(at[i].neighbor[j]);
            }
        }
        if num_t_groups != 0 && at[i].endpoint != 0 {
            list.push((num_atoms + at[i].endpoint as usize - 1) as u16);
        }
        list[0] = (list.len() - 1) as u16;
        pp[i] = list;
    }
    for i in 0..num_t_groups {
        let val = t_groups[i].n_num_endpoints as usize;
        let first = t_groups[i].n_first_endpoint_at_no_pos as usize;
        let mut list = vec![val as u16];
        list.extend_from_slice(&endpoints[first..first + val]);
        pp[num_atoms + i] = list;
    }
    pp
}

/// `SortedEquInfoToRanks`: converts sorted equivalence information to ranks;
/// returns the number of different ranks and whether a rank has changed
pub fn sorted_equ_info_to_ranks(
    n_symm_rank: &[u16],
    n_rank: &mut [u16],
    n_atom_number: &[u16],
    num_atoms: usize,
) -> (i32, bool) {
    let mut num_diff_ranks = 1;
    let mut num_changes = 0;
    if num_atoms > 0 {
        let i_init = num_atoms - 1;
        let j = n_atom_number[i_init] as usize;
        let mut r_old = n_symm_rank[j];
        let mut r_new = num_atoms as u16;
        n_rank[j] = num_atoms as u16;
        for i in (1..=i_init).rev() {
            let j = n_atom_number[i - 1] as usize;
            if n_symm_rank[j] != r_old {
                num_diff_ranks += 1;
                r_new = i as u16;
                num_changes += i32::from(r_old != r_new + 1);
                r_old = n_symm_rank[j];
            }
            n_rank[j] = r_new;
        }
    }
    (num_diff_ranks, num_changes != 0)
}

/// `SortedEquInfoToRanks` with `nSymmRank` and `nRank` being the same array
pub fn sorted_equ_info_to_ranks_in_place(
    n_rank: &mut [u16],
    n_atom_number: &[u16],
    num_atoms: usize,
) -> (i32, bool) {
    let copy = n_rank.to_vec();
    sorted_equ_info_to_ranks(&copy, n_rank, n_atom_number, num_atoms)
}

/// `SortedRanksToEquInfo`
pub fn sorted_ranks_to_equ_info(
    n_symm_rank: &mut [u16],
    n_rank: &[u16],
    n_atom_number: &[u16],
    num_atoms: usize,
) -> i32 {
    let mut j = n_atom_number[0] as usize;
    let mut r_old = n_rank[j];
    let mut r_new = 1u16;
    n_symm_rank[j] = 1;
    let mut num_diff_ranks = 1;
    for i in 1..num_atoms {
        j = n_atom_number[i] as usize;
        if n_rank[j] != r_old {
            num_diff_ranks += 1;
            r_new = (i + 1) as u16;
            r_old = n_rank[j];
        }
        n_symm_rank[j] = r_new;
    }
    num_diff_ranks
}

/// `SetNewRanksFromNeighLists3`: negative result if new ranks were found
pub fn set_new_ranks_from_neigh_lists3(
    num_atoms: usize,
    neigh_list: &NeighList,
    n_rank: &[u16],
    n_new_rank: &mut [u16],
    n_atom_number: &mut [u16],
) -> i32 {
    let (mut num_diff_ranks, mut num_new_ranks) = (0i32, 0i32);
    n_new_rank[..num_atoms].fill(0);
    let mut i = 0usize;
    let mut r1 = 1u16;
    while i < num_atoms {
        let mut j = n_atom_number[i] as usize;
        let mut r2 = n_rank[j];
        if r1 == r2 {
            n_new_rank[j] = r2;
            num_diff_ranks += 1;
            i += 1;
            r1 += 1;
            continue;
        }
        r1 = r2;
        insertions_sort_at_numbers(&mut n_atom_number[i..r2 as usize], |a, b| {
            comp_neigh_lists(neigh_list, n_rank, *a, *b)
        });
        j = r2 as usize - 1;
        let k = n_atom_number[j] as usize;
        n_new_rank[k] = r2;
        num_diff_ranks += 1;
        while j > i {
            if compare_neigh_list_lex(
                &neigh_list[n_atom_number[j - 1] as usize],
                &neigh_list[n_atom_number[j] as usize],
                n_rank,
            ) != 0
            {
                r2 = j as u16;
                num_diff_ranks += 1;
                num_new_ranks += 1;
            }
            j -= 1;
            n_new_rank[n_atom_number[j] as usize] = r2;
        }
        i = r1 as usize;
        r1 += 1;
    }
    if num_new_ranks != 0 {
        -num_diff_ranks
    } else {
        num_diff_ranks
    }
}

/// `SetNewRanksFromNeighLists4`: like version 3, ignoring neighbors with rank > `n_max_at_rank`
pub fn set_new_ranks_from_neigh_lists4(
    num_atoms: usize,
    neigh_list: &NeighList,
    n_rank: &[u16],
    n_new_rank: &mut [u16],
    n_atom_number: &mut [u16],
    n_max_at_rank: u16,
) -> i32 {
    let (mut num_diff_ranks, mut num_new_ranks) = (0i32, 0i32);
    n_new_rank[..num_atoms].fill(0);
    let mut i = 0usize;
    let mut r1 = 1u16;
    while i < num_atoms {
        let mut j = n_atom_number[i] as usize;
        let mut r2 = n_rank[j];
        if r1 == r2 {
            n_new_rank[j] = r2;
            num_diff_ranks += 1;
            i += 1;
            r1 += 1;
            continue;
        }
        r1 = r2;
        insertions_sort_at_numbers(&mut n_atom_number[i..r2 as usize], |a, b| {
            comp_neigh_lists_up_to_max_rank(neigh_list, n_rank, n_max_at_rank, *a, *b)
        });
        j = r2 as usize - 1;
        n_new_rank[n_atom_number[j] as usize] = r2;
        num_diff_ranks += 1;
        while j > i {
            if compare_neigh_list_lex_up_to_max_rank(
                &neigh_list[n_atom_number[j - 1] as usize],
                &neigh_list[n_atom_number[j] as usize],
                n_rank,
                n_max_at_rank,
            ) != 0
            {
                r2 = j as u16;
                num_diff_ranks += 1;
                num_new_ranks += 1;
            }
            j -= 1;
            n_new_rank[n_atom_number[j] as usize] = r2;
        }
        i = r1 as usize;
        r1 += 1;
    }
    if num_new_ranks != 0 {
        -num_diff_ranks
    } else {
        num_diff_ranks
    }
}

/// Which comparator `SetNewRanksFromNeighLists` sorts with
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RankSortCmp {
    /// `CompNeighListRanks`
    NeighListRanks,
    /// `CompNeighListRanksOrd`
    NeighListRanksOrd,
}

/// `SetNewRanksFromNeighLists`: `b_use_alt_sort & 1` selects the stable
/// insertion sort instead of `inchi_qsort`; returns the number of different ranks
pub fn set_new_ranks_from_neigh_lists(
    num_atoms: usize,
    neigh_list: &NeighList,
    n_rank: &[u16],
    n_new_rank: &mut [u16],
    n_atom_number: &mut [u16],
    b_use_alt_sort: i32,
    comp: RankSortCmp,
) -> i32 {
    let cmp = |a: &u16, b: &u16| match comp {
        RankSortCmp::NeighListRanks => comp_neigh_list_ranks(neigh_list, n_rank, *a, *b),
        RankSortCmp::NeighListRanksOrd => comp_neigh_list_ranks_ord(neigh_list, n_rank, *a, *b),
    };
    if b_use_alt_sort & 1 != 0 {
        insertions_sort(&mut n_atom_number[..num_atoms], cmp);
    } else {
        inchi_qsort(&mut n_atom_number[..num_atoms], cmp);
    }
    let mut num_diff_ranks = 1;
    if num_atoms > 0 {
        let mut current_rank = num_atoms as u16;
        n_new_rank[n_atom_number[num_atoms - 1] as usize] = current_rank;
        for i in (1..num_atoms).rev() {
            if comp_neigh_list_ranks(neigh_list, n_rank, n_atom_number[i - 1], n_atom_number[i])
                != 0
            {
                num_diff_ranks += 1;
                current_rank = i as u16;
            }
            n_new_rank[n_atom_number[i - 1] as usize] = current_rank;
        }
    }
    num_diff_ranks
}

/// `SortNeighListsBySymmAndCanonRank`
pub fn sort_neigh_lists_by_symm_and_canon_rank(
    num_atoms: usize,
    neigh_list: &mut NeighList,
    n_symm_rank: &[u16],
    n_canon_rank: &[u16],
) {
    for list in neigh_list.iter_mut().take(num_atoms) {
        insertions_sort_neigh_list_by_symm_and_canon_rank(list, n_symm_rank, n_canon_rank);
    }
}

/// `SortNeighLists2`
pub fn sort_neigh_lists2(
    num_atoms: usize,
    n_rank: &[u16],
    neigh_list: &mut NeighList,
    n_atom_number: &[u16],
) -> i32 {
    // (FIX_STEREOCOUNT_ERR is defined: every list with more than 1 neighbor is sorted)
    for &i in n_atom_number.iter().take(num_atoms) {
        if neigh_list[i as usize][0] > 1 {
            insertions_sort_neigh_list_at_numbers(&mut neigh_list[i as usize], n_rank);
        }
    }
    0
}

/// `SortNeighLists3`
pub fn sort_neigh_lists3(
    num_atoms: usize,
    n_rank: &[u16],
    neigh_list: &mut NeighList,
    n_atom_number: &[u16],
) -> i32 {
    let mut n_prev_rank = 0u16;
    for k in 0..num_atoms {
        let i = n_atom_number[k] as usize;
        if (n_rank[i] as usize != k + 1 || n_rank[i] == n_prev_rank) && neigh_list[i][0] > 1 {
            insertions_sort_neigh_list_at_numbers3(&mut neigh_list[i], n_rank);
        }
        n_prev_rank = n_rank[i];
    }
    0
}

/// `DifferentiateRanks2`; both rank arrays hold the result on return
#[allow(clippy::too_many_arguments)]
pub fn differentiate_ranks2(
    num_atoms: usize,
    neigh_list: &mut NeighList,
    _n_num_curr_ranks: i32,
    pn_curr_rank: &mut [u16],
    pn_prev_rank: &mut [u16],
    n_atom_number: &mut [u16],
    l_num_iter: &mut i64,
    b_use_alt_sort: i32,
) -> i32 {
    let mut n_num_curr_ranks;
    {
        let cr: &[u16] = pn_curr_rank;
        if b_use_alt_sort & 1 != 0 {
            insertions_sort(&mut n_atom_number[..num_atoms], |a, b| {
                comp_rank(cr, *a, *b)
            });
        } else {
            inchi_qsort(&mut n_atom_number[..num_atoms], |a, b| {
                comp_ranks_ord(cr, *a, *b)
            });
        }
    }
    let (mut cur, mut prev) = (pn_curr_rank, pn_prev_rank);
    loop {
        *l_num_iter += 1;
        std::mem::swap(&mut cur, &mut prev);
        sort_neigh_lists2(num_atoms, prev, neigh_list, n_atom_number);
        n_num_curr_ranks = set_new_ranks_from_neigh_lists(
            num_atoms,
            neigh_list,
            prev,
            cur,
            n_atom_number,
            1,
            RankSortCmp::NeighListRanksOrd,
        );
        if prev[..num_atoms] == cur[..num_atoms] {
            break;
        }
    }
    n_num_curr_ranks
}

/// `DifferentiateRanks3`
pub fn differentiate_ranks3(
    num_atoms: usize,
    neigh_list: &mut NeighList,
    _n_num_curr_ranks: i32,
    pn_curr_rank: &mut [u16],
    pn_prev_rank: &mut [u16],
    n_atom_number: &mut [u16],
    l_num_iter: &mut i64,
) -> i32 {
    let mut n_num_curr_ranks;
    let (mut cur, mut prev) = (pn_curr_rank, pn_prev_rank);
    loop {
        *l_num_iter += 1;
        std::mem::swap(&mut cur, &mut prev);
        sort_neigh_lists3(num_atoms, prev, neigh_list, n_atom_number);
        n_num_curr_ranks =
            set_new_ranks_from_neigh_lists3(num_atoms, neigh_list, prev, cur, n_atom_number);
        if n_num_curr_ranks >= 0 {
            break;
        }
    }
    n_num_curr_ranks
}

/// `DifferentiateRanks4`: ignores neighbors with rank > `n_max_at_rank`
#[allow(clippy::too_many_arguments)]
pub fn differentiate_ranks4(
    num_atoms: usize,
    neigh_list: &mut NeighList,
    _n_num_curr_ranks: i32,
    pn_curr_rank: &mut [u16],
    pn_prev_rank: &mut [u16],
    n_atom_number: &mut [u16],
    n_max_at_rank: u16,
    l_num_iter: &mut i64,
) -> i32 {
    let mut n_num_curr_ranks;
    let (mut cur, mut prev) = (pn_curr_rank, pn_prev_rank);
    loop {
        *l_num_iter += 1;
        std::mem::swap(&mut cur, &mut prev);
        sort_neigh_lists3(num_atoms, prev, neigh_list, n_atom_number);
        n_num_curr_ranks = set_new_ranks_from_neigh_lists4(
            num_atoms,
            neigh_list,
            prev,
            cur,
            n_atom_number,
            n_max_at_rank,
        );
        if n_num_curr_ranks >= 0 {
            break;
        }
    }
    n_num_curr_ranks
}

/// `DifferentiateRanksBasic`
pub fn differentiate_ranks_basic(
    num_atoms: usize,
    neigh_list: &mut NeighList,
    mut n_num_curr_ranks: i32,
    pn_curr_rank: &mut [u16],
    pn_prev_rank: &mut [u16],
    n_atom_number: &mut [u16],
    l_num_iter: &mut i64,
    b_use_alt_sort: i32,
) -> i32 {
    {
        let cr: &[u16] = pn_curr_rank;
        if b_use_alt_sort & 1 != 0 {
            insertions_sort(&mut n_atom_number[..num_atoms], |a, b| {
                comp_rank(cr, *a, *b)
            });
        } else {
            inchi_qsort(&mut n_atom_number[..num_atoms], |a, b| {
                comp_rank(cr, *a, *b)
            });
        }
    }
    let (mut cur, mut prev) = (pn_curr_rank, pn_prev_rank);
    loop {
        *l_num_iter += 1;
        let n_num_prev_ranks = n_num_curr_ranks;
        std::mem::swap(&mut cur, &mut prev);
        sort_neigh_lists2(num_atoms, prev, neigh_list, n_atom_number);
        n_num_curr_ranks = set_new_ranks_from_neigh_lists(
            num_atoms,
            neigh_list,
            prev,
            cur,
            n_atom_number,
            b_use_alt_sort,
            RankSortCmp::NeighListRanks,
        );
        if !(n_num_prev_ranks != n_num_curr_ranks || prev[..num_atoms] != cur[..num_atoms]) {
            break;
        }
    }
    n_num_curr_ranks
}

/// The counters of `CANON_STAT` that `BreakAllTies` updates
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct CanonCounters {
    pub l_num_neigh_list_iter: i64,
    pub l_num_break_ties: i64,
}

/// `BreakAllTies`: `rank_stack[0..2]` are the previous ranks and atom numbers;
/// `rank_stack[2..4]` (allocated here when empty) receive the new ones.
/// Returns the number of broken ties.
pub fn break_all_ties(
    num_atoms: usize,
    num_max: usize,
    rank_stack: &mut Vec<Vec<u16>>,
    neigh_list: &mut NeighList,
    n_temp_rank: &mut [u16],
    cs: &mut CanonCounters,
) -> i32 {
    while rank_stack.len() < 4 {
        rank_stack.push(Vec::new());
    }
    if rank_stack[2].is_empty() {
        rank_stack[2] = vec![0; num_max];
    }
    if rank_stack[3].is_empty() {
        rank_stack[3] = vec![0; num_max];
    }
    let prev_rank = rank_stack[0][..num_atoms].to_vec();
    let prev_atom_number = rank_stack[1][..num_atoms].to_vec();
    let mut n_new_atom_number = std::mem::take(&mut rank_stack[3]);
    let mut n_new_rank = std::mem::take(&mut rank_stack[2]);
    n_new_atom_number[..num_atoms].copy_from_slice(&prev_atom_number);
    n_new_rank[..num_atoms].copy_from_slice(&prev_rank);
    let mut n_ret = 0;
    for i in 1..num_atoms {
        if n_new_rank[n_new_atom_number[i - 1] as usize]
            == n_new_rank[n_new_atom_number[i] as usize]
        {
            n_new_rank[n_new_atom_number[i - 1] as usize] = i as u16;
            differentiate_ranks2(
                num_atoms,
                neigh_list,
                1,
                &mut n_new_rank,
                n_temp_rank,
                &mut n_new_atom_number,
                &mut cs.l_num_neigh_list_iter,
                1,
            );
            cs.l_num_break_ties += 1;
            n_ret += 1;
        }
    }
    rank_stack[2] = n_new_rank;
    rank_stack[3] = n_new_atom_number;
    n_ret
}

// ---------------------------------------------------------------------------
// ichimap1.c: mapping of stereo atoms and bonds, stereo CT comparison, CUR_TREE
// ---------------------------------------------------------------------------

use super::ichierr::CT_STEREOCOUNT_ERR;
use super::ichister::{
    AB_PARITY_UNDF, AB_PARITY_UNKN, AtStereoCarb, AtStereoDble, bond_chain_len, parity_known,
    parity_val,
};
use super::inpdef::MAX_NUM_STEREO_BONDS;

pub const STEREO_AT_MARK: i8 = 8;
/// `BEST_PARITY` / `WORSE_PARITY` (`CT_SMALLEST`)
pub const BEST_PARITY: i32 = 1;
pub const WORSE_PARITY: i32 = 2;

/// `All_SC_Same`: do all possibly stereogenic atoms that `canon_rank1` can be
/// mapped on have the same KNOWN parity? Returns how many such atoms there are, or 0.
pub fn all_sc_same(
    canon_rank1: u16,
    stack1: &[Vec<u16>],
    stack2: &[Vec<u16>],
    n_atom_number_canon_from: &[u16],
    at: &[SpAtom],
) -> i32 {
    let n1 = n_atom_number_canon_from[canon_rank1 as usize - 1] as usize;
    let r1 = stack1[0][n1];
    let i_max1 = r1 as usize;
    let mut b_found = 0;
    let mut stereo_atom_parity = -1;
    let mut i1 = 1;
    while i1 <= i_max1 {
        let s1 = stack2[1][i_max1 - i1] as usize;
        if r1 != stack2[0][s1] {
            break;
        }
        if at[s1].stereo_bond_neighbor[0] != 0 {
            b_found = 0;
            break;
        } else if i1 == 1 {
            stereo_atom_parity = parity_val(i32::from(at[s1].stereo_atom_parity));
            if !(1..=4).contains(&stereo_atom_parity) {
                b_found = 0;
                break;
            }
        } else if stereo_atom_parity != parity_val(i32::from(at[s1].stereo_atom_parity)) {
            b_found = 0;
            break;
        }
        b_found += 1;
        i1 += 1;
    }
    b_found
}

/// `Next_SC_At_CanonRank2`: the next canonical rank that can be mapped on a
/// not yet mapped stereocenter
pub fn next_sc_at_canon_rank2(
    canon_rank1: &mut u16,
    canon_rank1_min: &mut u16,
    b_first_time: &mut bool,
    b_atom_used_for_stereo: &[i8],
    stack1: &[Vec<u16>],
    stack2: &[Vec<u16>],
    n_atom_number_canon_from: &[u16],
    num_atoms: i32,
) -> bool {
    let mut canon_rank1_inp = *canon_rank1;
    if canon_rank1_inp < *canon_rank1_min {
        canon_rank1_inp = *canon_rank1_min;
    } else if canon_rank1_inp < 1 {
        canon_rank1_inp = 1;
    } else {
        canon_rank1_inp += 1;
    }
    let mut cr1 = canon_rank1_inp;
    let mut b_found = false;
    while i32::from(cr1) <= num_atoms {
        let n1 = n_atom_number_canon_from[cr1 as usize - 1] as usize;
        let r1 = stack1[0][n1];
        let i_max1 = r1 as usize;
        let mut i1 = 1;
        while i1 <= i_max1 {
            let s1 = stack2[1][i_max1 - i1] as usize;
            if r1 != stack2[0][s1] {
                break;
            }
            if b_atom_used_for_stereo[s1] == STEREO_AT_MARK {
                b_found = true;
                break;
            }
            i1 += 1;
        }
        if b_found {
            if *b_first_time {
                *canon_rank1_min = cr1;
                *b_first_time = false;
            }
            break;
        }
        cr1 += 1;
    }
    if b_found {
        *canon_rank1 = cr1;
        return true;
    }
    false
}

/// `CompareLinCtStereoDble`
pub fn compare_lin_ct_stereo_dble(
    d1: Option<&[AtStereoDble]>,
    len1: i32,
    d2: Option<&[AtStereoDble]>,
    len2: i32,
) -> i32 {
    if let (Some(a), Some(b)) = (d1, d2) {
        let num = len1.min(len2) as usize;
        let mut ret = 0;
        for i in 0..num {
            ret = i32::from(a[i].at_num1) - i32::from(b[i].at_num1);
            if ret != 0 {
                break;
            }
            ret = i32::from(a[i].at_num2) - i32::from(b[i].at_num2);
            if ret != 0 {
                break;
            }
            ret = i32::from(a[i].parity) - i32::from(b[i].parity);
            if ret != 0 {
                break;
            }
        }
        if ret == 0 {
            ret = len1 - len2;
        }
        ret
    } else if d1.is_some() && len1 > 0 {
        1
    } else if d2.is_some() && len2 > 0 {
        -1
    } else {
        0
    }
}

/// `CompareLinCtStereoCarb`
pub fn compare_lin_ct_stereo_carb(
    c1: Option<&[AtStereoCarb]>,
    len1: i32,
    c2: Option<&[AtStereoCarb]>,
    len2: i32,
) -> i32 {
    if let (Some(a), Some(b)) = (c1, c2) {
        let num = len1.min(len2) as usize;
        let mut ret = 0;
        for i in 0..num {
            ret = i32::from(a[i].at_num) - i32::from(b[i].at_num);
            if ret != 0 {
                break;
            }
            ret = i32::from(a[i].parity) - i32::from(b[i].parity);
            if ret != 0 {
                break;
            }
        }
        if ret == 0 {
            ret = len1 - len2;
        }
        ret
    } else if c1.is_some() && len1 > 0 {
        1
    } else if c2.is_some() && len2 > 0 {
        -1
    } else {
        0
    }
}

/// `CompareLinCtStereo`
#[allow(clippy::too_many_arguments)]
pub fn compare_lin_ct_stereo(
    d1: Option<&[AtStereoDble]>,
    dlen1: i32,
    c1: Option<&[AtStereoCarb]>,
    clen1: i32,
    d2: Option<&[AtStereoDble]>,
    dlen2: i32,
    c2: Option<&[AtStereoCarb]>,
    clen2: i32,
) -> i32 {
    let ret = compare_lin_ct_stereo_dble(d1, dlen1, d2, dlen2);
    if ret == 0 {
        compare_lin_ct_stereo_carb(c1, clen1, c2, clen2)
    } else {
        ret
    }
}

/// `CompareLinCtStereoAtomToValues` (`CT_SMALLEST`: greater than is `>`)
pub fn compare_lin_ct_stereo_atom_to_values(
    c: &AtStereoCarb,
    at_rank_canon1: u16,
    parity: u8,
) -> i32 {
    if c.at_num > at_rank_canon1 {
        return 1;
    }
    if c.at_num != at_rank_canon1 {
        return -1;
    }
    if c.parity > parity {
        return 1;
    }
    if c.parity != parity {
        return -1;
    }
    0
}

/// `CompareLinCtStereoDoubleToValues`
pub fn compare_lin_ct_stereo_double_to_values(
    d: &AtStereoDble,
    at_rank_canon1: u16,
    at_rank_canon2: u16,
    bond_parity: u8,
) -> i32 {
    if d.at_num1 > at_rank_canon1 {
        return 1;
    }
    if d.at_num1 != at_rank_canon1 {
        return -1;
    }
    if d.at_num2 > at_rank_canon2 {
        return 1;
    }
    if d.at_num2 != at_rank_canon2 {
        return -1;
    }
    if d.parity > bond_parity {
        return 1;
    }
    if d.parity != bond_parity {
        return -1;
    }
    0
}

/// `bUniqueAtNbrFromMappingRank`: the atom with mapping rank `n_at_rank` if it is unique
pub fn b_unique_at_nbr_from_mapping_rank(rank_stack: &[Vec<u16>], n_at_rank: u16) -> Option<u16> {
    let r = n_at_rank as usize - 1;
    let i = rank_stack[1][r];
    if n_at_rank == rank_stack[0][i as usize]
        && (r == 0 || n_at_rank != rank_stack[0][rank_stack[1][r - 1] as usize])
    {
        Some(i)
    } else {
        None
    }
}

/// `nGetMcr`: minimal class representative, compressing the path
pub fn n_get_mcr(n_eq_array: &mut [u16], n: u16) -> u16 {
    let mut n1 = n_eq_array[n as usize];
    if n == n1 {
        return n;
    }
    loop {
        let n2 = n_eq_array[n1 as usize];
        if n1 == n2 {
            break;
        }
        n1 = n2;
    }
    let mcr = n1;
    n1 = n;
    loop {
        let n2 = n_eq_array[n1 as usize];
        if mcr == n2 {
            break;
        }
        n_eq_array[n1 as usize] = mcr;
        n1 = n2;
    }
    mcr
}

/// `nJoin2Mcrs`: 1 if the classes were joined
pub fn n_join2_mcrs(n_eq_array: &mut [u16], n1: u16, n2: u16) -> i32 {
    let n1 = n_get_mcr(n_eq_array, n1);
    let n2 = n_get_mcr(n_eq_array, n2);
    if n1 < n2 {
        n_eq_array[n2 as usize] = n1;
        return 1;
    }
    if n2 < n1 {
        n_eq_array[n1 as usize] = n2;
        return 1;
    }
    0
}

fn find_sb_pos(at: &[SpAtom], a: usize, partner: usize) -> Option<usize> {
    for k in 0..MAX_NUM_STEREO_BONDS {
        let m = at[a].stereo_bond_neighbor[k] as usize;
        if m == 0 {
            return None;
        }
        if m - 1 == partner {
            return Some(k);
        }
    }
    None
}

/// `All_SB_Same`: the number of stereo bonds `canon_rank1`-`canon_rank2` can
/// be mapped on if all have the same known parity; 0 if not; -1 on an error
pub fn all_sb_same(
    canon_rank1: u16,
    canon_rank2: u16,
    stack1: &[Vec<u16>],
    stack2: &[Vec<u16>],
    n_atom_number_canon_from: &[u16],
    at: &[SpAtom],
) -> i32 {
    let n1 = n_atom_number_canon_from[canon_rank1 as usize - 1] as usize;
    let n2 = n_atom_number_canon_from[canon_rank2 as usize - 1] as usize;
    let r1 = stack1[0][n1];
    let r2 = stack1[0][n2];
    let i_max1 = r1 as usize;
    let (mut s1, mut s2, mut k1) = (0usize, 0usize, 0usize);
    let mut b_not_found = true;
    let mut i1 = 1;
    while i1 <= i_max1 {
        s1 = stack2[1][i_max1 - i1] as usize;
        if r1 != stack2[0][s1] {
            break;
        }
        k1 = 0;
        s2 = 0;
        b_not_found = true;
        while k1 < MAX_NUM_STEREO_BONDS {
            let t = at[s1].stereo_bond_neighbor[k1] as usize;
            if t == 0 {
                s2 = 0;
                break;
            }
            s2 = t - 1;
            b_not_found = r2 != stack2[0][s2];
            if !b_not_found {
                break;
            }
            k1 += 1;
        }
        if !b_not_found {
            break;
        }
        i1 += 1;
    }
    if b_not_found {
        return -1;
    }
    let Some(k2) = find_sb_pos(at, s2, s1) else {
        return -1;
    };
    let stereo_bond_parity = i32::from(at[s1].stereo_bond_parity[k1]);
    if !parity_known(stereo_bond_parity) {
        return 0;
    }
    let cumulene_len = bond_chain_len(stereo_bond_parity);
    let r_neigh1 = stack2[0][at[s1].neighbor[at[s1].stereo_bond_ord[k1] as usize] as usize];
    let r_neigh2 = stack2[0][at[s2].neighbor[at[s2].stereo_bond_ord[k2] as usize] as usize];
    let mut num_equal = 0;
    let mut sc = 1;
    while sc <= i_max1 {
        let i1 = stack2[1][i_max1 - sc] as usize;
        if r1 != stack2[0][i1] {
            break;
        }
        for k in 0..at[i1].valence as usize {
            let n1 = at[i1].neighbor[k] as usize;
            if r_neigh1 != stack2[0][n1] {
                continue;
            }
            let i2 = if cumulene_len != 0 {
                let (mut prev, mut len, mut next) = (i1, 0, n1);
                while len < cumulene_len {
                    if at[next].valence == 2 && at[next].num_h == 0 {
                        let j = usize::from(at[next].neighbor[0] as usize == prev);
                        prev = next;
                        next = at[next].neighbor[j] as usize;
                    } else {
                        break;
                    }
                    len += 1;
                }
                if len != cumulene_len || r2 != stack2[0][next] || r_neigh2 != stack2[0][prev] {
                    continue;
                }
                next
            } else {
                n1
            };
            let Some(kk1) = find_sb_pos(at, i1, i2) else {
                return 0;
            };
            let Some(kk2) = find_sb_pos(at, i2, i1) else {
                return 0;
            };
            if at[i1].stereo_bond_parity[kk1] != at[i2].stereo_bond_parity[kk2] {
                return -1;
            }
            if stereo_bond_parity != i32::from(at[i1].stereo_bond_parity[kk1]) {
                return 0;
            }
            num_equal += 1;
        }
        sc += 1;
    }
    num_equal
}

/// `Next_SB_At_CanonRanks2`: the next pair of canonical ranks that can be
/// mapped on a not yet mapped stereo bond
#[allow(clippy::too_many_arguments)]
pub fn next_sb_at_canon_ranks2(
    canon_rank1: &mut u16,
    canon_rank2: &mut u16,
    canon_rank1_min: &mut u16,
    canon_rank2_min: &mut u16,
    b_first_time: &mut bool,
    b_atom_used_for_stereo: &[i8],
    stack1: &[Vec<u16>],
    stack2: &[Vec<u16>],
    n_canon_rank_from: &[u16],
    n_atom_number_canon_from: &[u16],
    at: &[SpAtom],
    num_atoms: i32,
    b_allene: bool,
) -> bool {
    let mut canon_rank1_inp = *canon_rank1;
    let mut canon_rank2_inp = *canon_rank2;
    if canon_rank1_inp < *canon_rank1_min
        || (canon_rank1_inp == *canon_rank1_min && canon_rank2_inp < *canon_rank2_min)
    {
        canon_rank1_inp = *canon_rank1_min;
        canon_rank2_inp = *canon_rank2_min;
    } else if canon_rank1_inp < 2 {
        canon_rank1_inp = 2;
        canon_rank2_inp = 0;
    }
    let mut cr1 = canon_rank1_inp;
    let mut cr2 = num_atoms as u16;
    while i32::from(cr1) <= num_atoms {
        cr2 = cr1;
        let n1 = n_atom_number_canon_from[cr1 as usize - 1] as usize;
        let r1 = stack1[0][n1];
        let i_max1 = r1 as usize;
        let mut i1 = 1;
        while i1 <= i_max1 {
            let s1 = stack2[1][i_max1 - i1] as usize;
            if r1 != stack2[0][s1] {
                break;
            }
            if b_atom_used_for_stereo[s1] != 0 && b_atom_used_for_stereo[s1] < STEREO_AT_MARK {
                let mut k = 0;
                while k < MAX_NUM_STEREO_BONDS && at[s1].stereo_bond_neighbor[k] != 0 {
                    let s2 = at[s1].stereo_bond_neighbor[k] as usize - 1;
                    if b_atom_used_for_stereo[s2] != 0 {
                        let cumulene_len = bond_chain_len(i32::from(at[s1].stereo_bond_parity[k]));
                        if (cumulene_len % 2 != 0 && !b_allene)
                            || (cumulene_len % 2 == 0 && b_allene)
                        {
                            k += 1;
                            continue;
                        }
                        let r2 = stack2[0][s2];
                        let i_max2 = r2 as usize;
                        let mut i2 = 1;
                        while i2 <= i_max2 {
                            let n2 = stack1[1][i_max2 - i2] as usize;
                            if r2 != stack1[0][n2] {
                                break;
                            }
                            let mut m;
                            if cumulene_len != 0 {
                                m = 0;
                                while m < at[n1].valence as usize {
                                    let (mut prev, mut len, mut next) =
                                        (n1, 0, at[n1].neighbor[m] as usize);
                                    while len < cumulene_len {
                                        if at[next].valence == 2 && at[next].num_h == 0 {
                                            let j =
                                                usize::from(at[next].neighbor[0] as usize == prev);
                                            prev = next;
                                            next = at[next].neighbor[j] as usize;
                                        } else {
                                            break;
                                        }
                                        len += 1;
                                    }
                                    if len == cumulene_len && n2 == next {
                                        break;
                                    }
                                    m += 1;
                                }
                            } else {
                                m = 0;
                                while m < at[n1].valence as usize
                                    && n2 != at[n1].neighbor[m] as usize
                                {
                                    m += 1;
                                }
                            }
                            if m < at[n1].valence as usize
                                && n_canon_rank_from[n2] < cr2
                                && n_canon_rank_from[n2] > canon_rank2_inp
                            {
                                cr2 = n_canon_rank_from[n2];
                            }
                            i2 += 1;
                        }
                    }
                    k += 1;
                }
            }
            i1 += 1;
        }
        if cr2 >= cr1 {
            cr1 += 1;
            canon_rank2_inp = 0;
        } else {
            if *b_first_time {
                *canon_rank1_min = cr1;
                *canon_rank2_min = cr2;
                *b_first_time = false;
            }
            break;
        }
    }
    if cr1 > cr2 && i32::from(cr1) <= num_atoms {
        *canon_rank1 = cr1;
        *canon_rank2 = cr2;
        return true;
    }
    false
}

/// `NextStereoParity2Test`: steps through the (parity, calculated parity) pairs;
/// 1 means "go to the next pair of canonical ranks"
#[allow(clippy::too_many_arguments)]
pub fn next_stereo_parity2_test(
    stereo_bond_parity: &mut i32,
    sb_parity_calc: &mut i32,
    n_num_best: i32,
    n_num_worse: i32,
    n_num_unkn: i32,
    n_num_undf: i32,
    n_num_calc: i32,
    v_ab_parity_unknown: i32,
) -> i32 {
    loop {
        match *stereo_bond_parity {
            BEST_PARITY => match *sb_parity_calc {
                0 => {
                    *stereo_bond_parity = WORSE_PARITY;
                    if n_num_worse == 0 {
                        continue;
                    }
                }
                BEST_PARITY => {
                    *sb_parity_calc = WORSE_PARITY;
                    if n_num_best == 0 {
                        continue;
                    }
                }
                WORSE_PARITY => {
                    *stereo_bond_parity = WORSE_PARITY;
                    if n_num_calc == 0 {
                        continue;
                    }
                }
                _ => {}
            },
            WORSE_PARITY => match *sb_parity_calc {
                0 => {
                    *stereo_bond_parity = v_ab_parity_unknown;
                    if n_num_unkn == 0 {
                        continue;
                    }
                }
                BEST_PARITY => return CT_STEREOCOUNT_ERR,
                WORSE_PARITY => {
                    *sb_parity_calc = 0;
                    if n_num_worse == 0 {
                        continue;
                    }
                }
                _ => {}
            },
            AB_PARITY_UNKN => {
                if *sb_parity_calc != 0 {
                    return CT_STEREOCOUNT_ERR;
                }
                *stereo_bond_parity = AB_PARITY_UNDF;
                if n_num_undf == 0 {
                    return 1;
                }
            }
            AB_PARITY_UNDF => {
                if *sb_parity_calc != 0 {
                    return CT_STEREOCOUNT_ERR;
                }
                return 1;
            }
            _ => {}
        }
        return 0;
    }
}

/// `SetUseAtomForStereo`: 0 for atoms without parity, `STEREO_AT_MARK` for
/// stereocenters, otherwise the number of adjacent stereo bonds
pub fn set_use_atom_for_stereo(b_atom_used_for_stereo: &mut [i8], at: &[SpAtom], num_atoms: usize) {
    b_atom_used_for_stereo[..num_atoms].fill(0);
    for i in 0..num_atoms {
        if at[i].parity != 0 {
            let mut k = 0;
            while k < MAX_NUM_STEREO_BONDS && at[i].stereo_bond_neighbor[k] != 0 {
                k += 1;
            }
            b_atom_used_for_stereo[i] = if k != 0 { k as i8 } else { STEREO_AT_MARK };
        }
    }
}

/// `CUR_TREE`: a stack of (rank, atoms..., count) segments
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CurTree {
    pub tree: Vec<u16>,
    pub cur_len: i32,
    pub max_len: i32,
    pub incr_len: i32,
}

impl CurTree {
    /// `CurTreeAlloc`
    pub fn alloc(&mut self, num_atoms: i32) -> i32 {
        if !self.tree.is_empty() && self.max_len > 0 && self.max_len % num_atoms == 0 {
            self.cur_len = 0;
            self.incr_len = num_atoms;
            let n = self.max_len as usize;
            self.tree[..n].fill(0);
            return 0;
        }
        *self = CurTree::default();
        if num_atoms > 0 {
            self.tree = vec![0; num_atoms as usize];
            self.incr_len = num_atoms;
            self.max_len = num_atoms;
            return 0;
        }
        -1
    }

    /// `CurTreeReAlloc`
    pub fn realloc(&mut self) -> i32 {
        if !self.tree.is_empty() && self.max_len > 0 && self.incr_len > 0 {
            self.tree.resize((self.max_len + self.incr_len) as usize, 0);
            self.max_len += self.incr_len;
            return 0;
        }
        -1
    }

    /// `CurTreeFree`
    pub fn free(&mut self) {
        *self = CurTree::default();
    }

    /// `CurTreeAddRank`
    pub fn add_rank(&mut self, rank: u16) -> i32 {
        if self.cur_len + 2 > self.max_len && self.realloc() != 0 {
            return -1;
        }
        self.tree[self.cur_len as usize] = rank;
        self.tree[self.cur_len as usize + 1] = 1;
        self.cur_len += 2;
        0
    }

    /// `CurTreeIsLastRank`
    pub fn is_last_rank(&self, rank: u16) -> i32 {
        if self.cur_len > 0 {
            let mut rank_pos = self.cur_len - 1;
            rank_pos -= i32::from(self.tree[rank_pos as usize]);
            if rank_pos >= 0 {
                return i32::from(rank == self.tree[rank_pos as usize]);
            }
        }
        0
    }

    /// `CurTreeRemoveLastRankIfNoAtoms`
    pub fn remove_last_rank_if_no_atoms(&mut self) -> i32 {
        if !self.tree.is_empty() && self.cur_len >= 2 {
            if self.tree[self.cur_len as usize - 1] == 1 {
                return self.remove_last_rank();
            }
            return 1;
        }
        -1
    }

    /// `CurTreeAddAtom`
    pub fn add_atom(&mut self, at_no: i32) -> i32 {
        if self.cur_len + 1 > self.max_len && self.realloc() != 0 {
            return -1;
        }
        if self.cur_len > 0 {
            self.cur_len -= 1;
            let new_len = self.tree[self.cur_len as usize] + 1;
            self.tree[self.cur_len as usize] = at_no as u16;
            self.cur_len += 1;
            self.tree[self.cur_len as usize] = new_len;
            self.cur_len += 1;
            return 0;
        }
        -1
    }

    /// `CurTreeKeepLastAtomsOnly`
    pub fn keep_last_atoms_only(&mut self, tpos: i32, mut shift: i32) {
        if self.tree.is_empty() {
            return;
        }
        let cur_length_pos = self.cur_len - shift;
        if cur_length_pos > tpos {
            let seg = i32::from(self.tree[cur_length_pos as usize]);
            if seg > 2 {
                self.cur_len -= seg - 2;
                let dst = (cur_length_pos - seg + 1) as usize;
                let src = (cur_length_pos - 1) as usize;
                let count = (shift + 1) as usize;
                self.tree.copy_within(src..src + count, dst);
                self.tree[(self.cur_len - shift) as usize] = 2;
                shift += 3;
            } else {
                shift += seg + 1;
            }
            self.keep_last_atoms_only(tpos, shift);
        }
    }

    /// `CurTreeRemoveIfLastAtom`
    pub fn remove_if_last_atom(&mut self, at_no: i32) -> i32 {
        if !self.tree.is_empty() && self.cur_len > 2 {
            let len = self.tree[self.cur_len as usize - 1];
            if len >= 2 && i32::from(self.tree[self.cur_len as usize - 2]) == at_no {
                self.cur_len -= 1;
                self.tree[self.cur_len as usize - 1] = len - 1;
                return 0;
            }
            return 1;
        }
        -1
    }

    /// `CurTreeGetPos`
    pub fn get_pos(&self) -> i32 {
        self.cur_len
    }

    /// `CurTreeSetPos`
    pub fn set_pos(&mut self, len: i32) -> i32 {
        self.cur_len = len;
        0
    }

    /// `CurTreeRemoveLastRank`
    pub fn remove_last_rank(&mut self) -> i32 {
        if self.cur_len > 0 {
            self.cur_len -= i32::from(self.tree[self.cur_len as usize - 1]) + 1;
            if self.cur_len >= 0 {
                return 0;
            }
        }
        -1
    }

    /// `CurTreeIsLastAtomEqu`
    pub fn is_last_atom_equ(&self, at_no: i32, n_symm_stereo: &[u16]) -> i32 {
        if !self.tree.is_empty() && self.cur_len > 1 {
            let n_eq = n_symm_stereo[at_no as usize];
            let end = self.cur_len - 1;
            let mut len = i32::from(self.tree[end as usize]) - 1;
            while len > 0 {
                if n_symm_stereo[self.tree[(end - len) as usize] as usize] == n_eq {
                    return 1;
                }
                len -= 1;
            }
            return 0;
        }
        -1
    }
}
