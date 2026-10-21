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
