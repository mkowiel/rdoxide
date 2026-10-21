//! The canonical numbering search (the McKay-style partition refinement
//! algorithm) of `ichican2.c`: partitions, cells, node sets and connection
//! table containers.
//!
//! Memory management of the C code (`...Create`/`...Free`) is replaced by
//! vectors; the `static` rank mark/mask bits are constants (`AT_NUMB` is 16
//! bits wide, so the mark bit is `0x8000`). `CANON_DATA` is represented by the
//! few values the helpers need.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi_bns::NodeSet;
use super::ichierr::{CT_CANON_ERR, CT_OUT_OF_RAM};
use super::ichimap::{NeighList, differentiate_ranks3, differentiate_ranks4};

pub const INCHI_CANON_INFINITY: u16 = 0x7FFF;
/// `rank_mark_bit` / `rank_mask_bit` (set by `SetBitCreate`)
pub const RANK_MARK_BIT: u16 = 0x8000;
pub const RANK_MASK_BIT: u16 = 0x7FFF;
/// bits in a `bitWord`
const NUM_BIT: usize = 16;

pub type Node = u16;

/// `Cell`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Cell {
    pub first: i32,
    pub next: i32,
    pub prev: i32,
}

/// `Partition`
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Partition {
    pub rank: Vec<u16>,
    pub at_number: Vec<u16>,
}

/// `UnorderedPartition`
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct UnorderedPartition {
    pub equ2: Vec<u16>,
}

/// `Transposition`
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Transposition {
    pub n_at_numb: Vec<u16>,
}

/// `ConTable`: a connection table (a prefix of the one being searched for)
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ConTable {
    pub ctbl: Vec<u16>,
    pub len_ct: i32,
    pub n_len_ct_at_only: i32,
    pub maxlen_ct: i32,
    pub max_pos: i32,
    pub max_vert: i32,
    pub len_pos: i32,
    pub next_at_rank: Vec<u16>,
    pub next_ctbl_pos: Vec<u16>,
    pub num_h: Vec<i16>,
    pub len_num_h: i32,
    pub maxlen_num_h: i32,
    pub num_h_fixed: Vec<i16>,
    pub iso_sort_key: Vec<i64>,
    pub len_iso_sort_key: i32,
    pub maxlen_iso_sort_key: i32,
    pub iso_exchg_atnos: Vec<i8>,
    pub len_iso_exchg_atnos: i32,
    pub maxlen_iso_exchg_atnos: i32,
}

/// `kLeast`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct KLeast {
    pub k: i32,
    pub i: i32,
}

/// `CTableCreate`; the `Option` sizes are the `pCD->maxlen...` values of the
/// arrays that exist in `CANON_DATA` (`NumH`, `NumHfixed`, `iso_sort_key`,
/// `iso_exchg_atnos`)
pub fn c_table_create(
    n: usize,
    n_max_len_linear_ct: i32,
    n_len_ct_at_only: i32,
    maxlen_num_h: Option<i32>,
    maxlen_num_h_fixed: Option<i32>,
    maxlen_iso_sort_key: Option<i32>,
    maxlen_iso_exchg_atnos: Option<i32>,
) -> ConTable {
    let maxlen_ct = n_max_len_linear_ct + 1;
    let maxlen_num_h = maxlen_num_h.map_or(0, |m| m + 1);
    let maxlen_num_h_fixed = maxlen_num_h_fixed.map_or(0, |m| m + 1);
    let maxlen_iso = maxlen_iso_sort_key.map_or(0, |m| m + 1);
    let maxlen_iso_exchg = maxlen_iso_exchg_atnos.map_or(0, |m| m + 1);
    let _ = maxlen_iso_exchg;
    ConTable {
        ctbl: vec![0; maxlen_ct as usize],
        next_ctbl_pos: vec![0; n + 1],
        next_at_rank: vec![0; n + 1],
        num_h: vec![0; maxlen_num_h as usize],
        num_h_fixed: vec![0; maxlen_num_h_fixed as usize],
        iso_sort_key: vec![0; maxlen_iso as usize],
        iso_exchg_atnos: vec![0; maxlen_iso_exchg as usize],
        len_ct: 0,
        n_len_ct_at_only,
        maxlen_ct,
        len_num_h: 0,
        maxlen_num_h,
        len_iso_sort_key: 0,
        maxlen_iso_sort_key: maxlen_iso,
        len_iso_exchg_atnos: 0,
        // (sic) the C code copies maxlenIso here
        maxlen_iso_exchg_atnos: maxlen_iso,
        max_pos: n as i32 + 1,
        max_vert: n as i32,
        len_pos: 0,
    }
}

/// `UnorderedPartitionMakeDiscrete`
pub fn unordered_partition_make_discrete(p: &mut UnorderedPartition, n: usize) {
    for (i, e) in p.equ2[..n].iter_mut().enumerate() {
        *e = i as u16;
    }
}

/// `PartitionIsDiscrete`
pub fn partition_is_discrete(p: &Partition, n: usize) -> bool {
    let mut r = 1u16;
    for i in 0..n {
        if r != (RANK_MASK_BIT & p.rank[p.at_number[i] as usize]) {
            return false;
        }
        r = r.wrapping_add(1);
    }
    true
}

/// `PartitionGetFirstCell`: finds the first non-trivial cell after cell `k-1`
/// (stored in `base_w[k-1]`); returns its size or 0
pub fn partition_get_first_cell(p: &Partition, base_w: &mut [Cell], k: usize, n: usize) -> i32 {
    let mut i = if k > 1 {
        base_w[k - 2].first as usize + 1
    } else {
        0
    };
    if i < n {
        let mut r = (i + 1) as u16;
        while i < n && r == (RANK_MASK_BIT & p.rank[p.at_number[i] as usize]) {
            i += 1;
            r = r.wrapping_add(1);
        }
    }
    let w = &mut base_w[k - 1];
    if i < n {
        w.first = i as i32;
        let r = RANK_MASK_BIT & p.rank[p.at_number[i] as usize];
        i += 1;
        while i < n && r == (RANK_MASK_BIT & p.rank[p.at_number[i] as usize]) {
            i += 1;
        }
        w.next = i as i32;
        return w.next - w.first;
    }
    w.first = INCHI_CANON_INFINITY as i32;
    w.next = 0;
    0
}

/// `CellMakeEmpty`
pub fn cell_make_empty(base_w: &mut [Cell], k: usize) {
    let c = &mut base_w[k - 1];
    c.first = INCHI_CANON_INFINITY as i32;
    c.next = 0;
    c.prev = -1;
}

/// `NodeSetFromVertices`: the set `l` (1-based) becomes `v[..num_v]` (1-based vertices)
pub fn node_set_from_vertices(cur_nodes: &mut NodeSet, l: usize, v: &[Node], num_v: usize) {
    let bits = &mut cur_nodes.bitword[l - 1];
    bits.iter_mut().for_each(|b| *b = 0);
    for &x in &v[..num_v] {
        let j = x as usize - 1;
        bits[j / NUM_BIT] |= 1 << (j % NUM_BIT);
    }
}

/// `AllNodesAreInSet`: is every node of set `lcur_nodes` of `cur_nodes` in set `lset` of `set`?
pub fn all_nodes_are_in_set(
    cur_nodes: &NodeSet,
    lcur_nodes: usize,
    set: &NodeSet,
    lset: usize,
) -> bool {
    let a = &cur_nodes.bitword[lcur_nodes - 1];
    let b = &set.bitword[lset - 1];
    (0..cur_nodes.len_set).all(|i| a[i] & !b[i] == 0)
}

/// `PartitionGetMcrAndFixSet`: the minimal class representatives and the fixed
/// vertices (single vertex cells) of a partition, as sets number `l` (1-based)
pub fn partition_get_mcr_and_fix_set(
    p: &Partition,
    mcr: &mut NodeSet,
    fix: &mut NodeSet,
    n: usize,
    l: usize,
) {
    mcr.bitword[l - 1].iter_mut().for_each(|b| *b = 0);
    fix.bitword[l - 1].iter_mut().for_each(|b| *b = 0);
    let mut i = 0usize;
    let mut r = 1u16;
    while i < n {
        let mut j1 = p.at_number[i] as usize;
        let r1 = RANK_MASK_BIT & p.rank[j1];
        if r == r1 {
            fix.bitword[l - 1][j1 / NUM_BIT] |= 1 << (j1 % NUM_BIT);
            mcr.bitword[l - 1][j1 / NUM_BIT] |= 1 << (j1 % NUM_BIT);
        } else {
            r = r1;
            while i + 1 < n {
                let j2 = p.at_number[i + 1] as usize;
                if r != (RANK_MASK_BIT & p.rank[j2]) {
                    break;
                }
                if j1 > j2 {
                    j1 = j2;
                }
                i += 1;
            }
            mcr.bitword[l - 1][j1 / NUM_BIT] |= 1 << (j1 % NUM_BIT);
        }
        i += 1;
        r = r.wrapping_add(1);
    }
}

/// `PartitionGetTransposition`
pub fn partition_get_transposition(
    p_from: &Partition,
    p_to: &Partition,
    n: usize,
    gamma: &mut Transposition,
) {
    for i in 0..n {
        gamma.n_at_numb[p_from.at_number[i] as usize] = p_to.at_number[i];
    }
}

/// `nGetMcr2`: the minimal class representative, partially compressing the classes
pub fn n_get_mcr2(n_equ_array: &mut [u16], n: u16) -> u16 {
    let mut n1 = n_equ_array[n as usize];
    if n == n1 {
        return n;
    }
    loop {
        let n2 = n_equ_array[n1 as usize];
        if n1 == n2 {
            break;
        }
        n1 = n2;
    }
    let mcr = n1;
    let mut n1 = n;
    loop {
        let n2 = n_equ_array[n1 as usize];
        if mcr == n2 {
            break;
        }
        n_equ_array[n1 as usize] = mcr;
        n1 = n2;
    }
    mcr
}

/// `nJoin2Mcrs2`: joins the classes of `n1` and `n2`; 1 if a change was made
pub fn n_join2_mcrs2(n_equ_array: &mut [u16], n1: u16, n2: u16) -> i32 {
    let n1 = n_get_mcr2(n_equ_array, n1);
    let n2 = n_get_mcr2(n_equ_array, n2);
    if n1 < n2 {
        n_equ_array[n2 as usize] = n1;
        return 1;
    }
    if n2 < n1 {
        n_equ_array[n1 as usize] = n2;
        return 1;
    }
    0
}

/// `GetUnorderedPartitionMcrNode`
pub fn get_unordered_partition_mcr_node(p1: &mut UnorderedPartition, v: Node) -> Node {
    1 + n_get_mcr2(&mut p1.equ2, v.wrapping_sub(1))
}

/// `UnorderedPartitionJoin`: changes `p2` to `p2 v p1`; returns the number of changes
pub fn unordered_partition_join(
    p1: &UnorderedPartition,
    p2: &mut UnorderedPartition,
    n: usize,
) -> i32 {
    let mut changes = 0;
    for i in 0..n {
        let j = p1.equ2[i] as usize;
        if j == i || p2.equ2[i] == p2.equ2[j] {
            continue;
        }
        changes += n_join2_mcrs2(&mut p2.equ2, i as u16, j as u16);
    }
    changes
}

/// `PartitionSatisfiesLemma_2_25`
pub fn partition_satisfies_lemma_2_25(p: &Partition, n: usize) -> bool {
    let mut partition_size = 0usize;
    let mut num_non_trivial_cells = 0usize;
    let mut num = 0;
    let mut r = 1u16;
    for i in 0..n {
        if (RANK_MASK_BIT & p.rank[p.at_number[i] as usize]) == r {
            partition_size += 1;
            if num != 0 {
                num_non_trivial_cells += 1;
                num = 0;
            }
        } else {
            num += 1;
        }
        r = r.wrapping_add(1);
    }
    n <= partition_size + 4
        || n == partition_size + num_non_trivial_cells
        || n == partition_size + num_non_trivial_cells + 1
}

/// `PartitionCopy`
pub fn partition_copy(to: &mut Partition, from: &Partition, n: usize) {
    to.at_number[..n].copy_from_slice(&from.at_number[..n]);
    for i in 0..n {
        to.rank[i] = from.rank[i] & RANK_MASK_BIT;
    }
}

/// `PartitionColorVertex`: makes the equitable partition `p[1]` out of `p[0]`
/// by reducing the rank of vertex `v` (1-based); `p[2]` is a temporary
/// rank array. Returns the number of different ranks or an error.
#[allow(clippy::too_many_arguments)]
pub fn partition_color_vertex(
    g: &mut NeighList,
    p: &mut [Partition],
    v: Node,
    n: usize,
    n_tg: usize,
    n_max: usize,
    b_digraph: bool,
    n_num_prev_ranks: i32,
) -> i32 {
    for q in p.iter_mut().take(3).skip(1) {
        if q.at_number.is_empty() {
            q.at_number = vec![0; n_max];
        }
        if q.rank.is_empty() {
            q.rank = vec![0; n_max];
        }
        if q.at_number.is_empty() || q.rank.is_empty() {
            return CT_OUT_OF_RAM;
        }
    }
    {
        let (a, b) = p.split_at_mut(1);
        partition_copy(&mut b[0], &a[0], n_tg);
    }
    let sv = v.wrapping_sub(1);
    if sv as usize >= n_tg {
        return CT_CANON_ERR;
    }
    let (_, rest) = p.split_at_mut(1);
    let (p1, p2) = rest.split_at_mut(1);
    let p1 = &mut p1[0];
    let p2 = &mut p2[0];
    let rv = p1.rank[sv as usize];
    let mut s = (n_max + 1) as u16;
    let mut j = rv as i32 - 1;
    while j >= 0 {
        s = p1.at_number[j as usize];
        if rv != p1.rank[s as usize] || s == sv {
            break;
        }
        j -= 1;
    }
    if s != sv {
        return CT_CANON_ERR;
    }
    let mut r = rv.wrapping_sub(1);
    let mut i = j;
    j -= 1;
    while j >= 0 {
        s = p1.at_number[j as usize];
        r = p1.rank[s as usize];
        if rv != r {
            break;
        }
        p1.at_number[i as usize] = s;
        i = j;
        j -= 1;
    }
    let r = if i > 0 { r.wrapping_add(1) } else { 1 };
    p1.at_number[i as usize] = sv;
    p1.rank[sv as usize] = r;
    let mut l_num_iter = 0i64;
    if b_digraph {
        differentiate_ranks4(
            n_tg,
            g,
            n_num_prev_ranks + 1,
            &mut p1.rank,
            &mut p2.rank,
            &mut p1.at_number,
            n as u16,
            &mut l_num_iter,
        )
    } else {
        differentiate_ranks3(
            n_tg,
            g,
            n_num_prev_ranks + 1,
            &mut p1.rank,
            &mut p2.rank,
            &mut p1.at_number,
            &mut l_num_iter,
        )
    }
}

/// `CellGetMinNode`: the minimal unmarked node (1-based) of the cell `w` that
/// is greater than `v` (`v` = previous node, 0 on the first call), or
/// `INCHI_CANON_INFINITY`. With `aux_rank` the nodes are ordered by it first.
pub fn cell_get_min_node(p: &Partition, w: &Cell, v: Node, aux_rank: Option<&[u16]>) -> Node {
    if w.first > w.next {
        return INCHI_CANON_INFINITY;
    }
    let mut u_min_at_numb = INCHI_CANON_INFINITY;
    if let Some(aux) = aux_rank {
        let mut i = w.first;
        while i < w.next {
            let cur = p.at_number[i as usize];
            if p.rank[cur as usize] & RANK_MARK_BIT == 0 {
                break;
            }
            i += 1;
        }
        if i == w.next {
            return INCHI_CANON_INFINITY;
        }
        let mut u_min_aux_rank = INCHI_CANON_INFINITY as i32;
        let mut n_min_at_numb = INCHI_CANON_INFINITY as i32;
        let (n_inp_at_numb, u_inp_aux_rank) = if v != 0 {
            (v as i32 - 1, aux[v as usize - 1] as i32)
        } else {
            (-1, 0)
        };
        while i < w.next {
            let n_cur = p.at_number[i as usize] as i32;
            if p.rank[n_cur as usize] & RANK_MARK_BIT == 0 {
                let u_cur_aux_rank = aux[n_cur as usize] as i32;
                if (u_cur_aux_rank == u_inp_aux_rank && n_cur > n_inp_at_numb)
                    || u_cur_aux_rank > u_inp_aux_rank
                {
                    if u_cur_aux_rank == u_min_aux_rank && n_cur < n_min_at_numb {
                        n_min_at_numb = n_cur;
                    } else if u_cur_aux_rank < u_min_aux_rank {
                        u_min_aux_rank = u_cur_aux_rank;
                        n_min_at_numb = n_cur;
                    }
                }
            }
            i += 1;
        }
        u_min_at_numb = if n_min_at_numb == INCHI_CANON_INFINITY as i32 {
            INCHI_CANON_INFINITY
        } else {
            n_min_at_numb as u16
        };
    } else {
        for i in w.first..w.next {
            let cur = p.at_number[i as usize];
            if cur >= v && p.rank[cur as usize] & RANK_MARK_BIT == 0 && cur < u_min_at_numb {
                u_min_at_numb = cur;
            }
        }
    }
    if u_min_at_numb != INCHI_CANON_INFINITY {
        u_min_at_numb += 1;
    }
    u_min_at_numb
}

/// `CellGetNumberOfNodes`: the number of unmarked nodes in the cell
pub fn cell_get_number_of_nodes(p: &Partition, w: &Cell) -> i32 {
    (w.first..w.next)
        .filter(|&i| RANK_MARK_BIT & p.rank[p.at_number[i as usize] as usize] == 0)
        .count() as i32
}

/// `CellIntersectWithSet`: marks the nodes of the cell that are not in the
/// set `l` (1-based) of `mcr`; returns the number of newly marked nodes
pub fn cell_intersect_with_set(p: &mut Partition, w: &Cell, mcr: &NodeSet, l: usize) -> i32 {
    let mcr_bits = &mcr.bitword[l - 1];
    if w.first >= w.next {
        return 0;
    }
    let mut k = 0;
    for i in w.first..w.next {
        let j = p.at_number[i as usize] as usize;
        if mcr_bits[j / NUM_BIT] & (1 << (j % NUM_BIT)) == 0 {
            k += i32::from(p.rank[j] & RANK_MARK_BIT == 0);
            p.rank[j] |= RANK_MARK_BIT;
        }
    }
    k
}

/// `CtPartClear`: clears the part of the connection table from position `k` on
pub fn ct_part_clear(ct: &mut ConTable, k: usize) {
    let start = if k > 1 {
        ct.next_ctbl_pos[k - 1] as i32
    } else {
        0
    };
    let len = ct.len_ct - start;
    if len > 0 {
        ct.ctbl[start as usize..ct.len_ct as usize].fill(0);
    }
    ct.len_ct = start;
    ct.len_pos = k as i32;
}

/// `insertions_sort_NeighList_AT_NUMBERS2`: sorts the neighbors of the list
/// `base` (`base[0]` = count) by ascending rank (the mark bit is ignored)
pub fn insertions_sort_neigh_list_at_numbers2(base: &mut [u16], n_rank: &[u16], max_rj: u16) {
    let num = base[0] as usize;
    let nb = &mut base[1..];
    for k in 1..num {
        // i = k - 1, j = k
        let mut j = k;
        let rj = RANK_MASK_BIT & n_rank[nb[j] as usize];
        if rj < max_rj {
            while j > 0 && rj < (RANK_MASK_BIT & n_rank[nb[j - 1] as usize]) {
                nb.swap(j - 1, j);
                j -= 1;
            }
        }
    }
}
