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

use super::ichi_bns::{NodeSet, node_set_create};
use super::ichierr::{CT_CANON_ERR, CT_OUT_OF_RAM, CT_TIMEOUT_ERR};
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

// ---------------------------------------------------------------------------
// connection table parts: filling, comparing and copying
// ---------------------------------------------------------------------------

use super::ichitaut::T_NUM_NO_ISOTOPIC;

pub const EMPTY_CT: u16 = 0;
pub const MAX_LAYERS: usize = 100;

/// The parts of `CANON_DATA` that the connection table helpers use. An
/// array is present when the C pointer is not NULL.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CanonData {
    /// `nCanonFlags`
    pub n_canon_flags: i32,
    /// `LinearCT`
    pub linear_ct: Vec<u16>,
    pub n_max_len_linear_ct: i32,
    pub n_len_linear_ct: i32,
    pub n_len_ct_at_only: i32,
    pub num_h: Option<Vec<i16>>,
    pub maxlen_num_h: i32,
    pub num_h_fixed: Option<Vec<i16>>,
    pub maxlen_num_h_fixed: i32,
    pub iso_sort_key: Option<Vec<i64>>,
    pub maxlen_iso_sort_key: i32,
    pub iso_exchg_atnos: Option<Vec<i8>>,
    pub maxlen_iso_exchg_atnos: i32,
    pub n_aux_rank: Option<Vec<u16>>,
}

/// `CtPartFill`: appends the part of the connection table for the cells of
/// the partition `p` that became singletons since level `k - 1`
/// (`k` = new level + 1). The neighbor lists of `g` are sorted in place.
#[allow(clippy::too_many_arguments)]
pub fn ct_part_fill(
    g: &mut NeighList,
    cd: &CanonData,
    p: &Partition,
    ct: &mut ConTable,
    k: usize,
    n: usize,
    n_tg: usize,
    n_max: usize,
) {
    let k = k.wrapping_sub(1);
    let (mut start_ctbl, start_at_ord): (usize, i32) = if k != 0 {
        (
            ct.next_ctbl_pos[k - 1] as usize,
            i32::from(ct.next_at_rank[k - 1]) - 1,
        )
    } else {
        (0, 0)
    };
    if start_at_ord < 0 {
        return;
    }
    let start_at_ord = start_at_ord as usize;
    if start_at_ord >= n_max {
        return;
    }
    let an_sao = p.at_number[start_at_ord] as usize;
    let mut r = RANK_MASK_BIT & p.rank[an_sao];
    let mut i = start_at_ord;
    while i < n_tg {
        let m = p.at_number[i] as usize;
        if r != (RANK_MASK_BIT & p.rank[m]) {
            break;
        }
        ct.ctbl[start_ctbl] = r;
        start_ctbl += 1;
        insertions_sort_neigh_list_at_numbers2(&mut g[m], &p.rank, r);
        let nn = g[m][0] as usize;
        let mut j = 1;
        while j <= nn {
            let rj = RANK_MASK_BIT & p.rank[g[m][j] as usize];
            if rj >= r {
                break;
            }
            ct.ctbl[start_ctbl] = rj;
            start_ctbl += 1;
            j += 1;
        }
        i += 1;
        r = r.wrapping_add(1);
    }
    // base hydrogen atoms
    match (&cd.num_h, ct.num_h.is_empty()) {
        (Some(src), false) => {
            let mut nn = n.min(i);
            let mut j = start_at_ord;
            while j < nn {
                ct.num_h[j] = src[p.at_number[j] as usize];
                j += 1;
            }
            while j < i {
                let data_pos = n + T_NUM_NO_ISOTOPIC * (p.at_number[j] as usize - n);
                for d in 0..T_NUM_NO_ISOTOPIC {
                    ct.num_h[nn] = src[data_pos + d];
                    nn += 1;
                }
                j += 1;
            }
            ct.len_num_h = nn as i32;
        }
        _ => ct.len_num_h = 0,
    }
    // fixed hydrogen atoms
    if let (Some(src), false) = (&cd.num_h_fixed, ct.num_h_fixed.is_empty()) {
        let nn = n.min(i);
        for j in start_at_ord..nn {
            ct.num_h_fixed[j] = src[p.at_number[j] as usize];
        }
    }
    // isotopic keys
    match (&cd.iso_sort_key, ct.iso_sort_key.is_empty()) {
        (Some(src), false) => {
            for j in start_at_ord..i {
                ct.iso_sort_key[j] = src[p.at_number[j] as usize];
            }
            ct.len_iso_sort_key = i as i32;
        }
        _ => ct.len_iso_sort_key = 0,
    }
    match (&cd.iso_exchg_atnos, ct.iso_exchg_atnos.is_empty()) {
        (Some(src), false) => {
            for j in start_at_ord..i {
                ct.iso_exchg_atnos[j] = src[p.at_number[j] as usize];
            }
            ct.len_iso_exchg_atnos = i as i32;
        }
        _ => ct.len_iso_exchg_atnos = 0,
    }
    ct.len_ct = start_ctbl as i32;
    ct.next_ctbl_pos[k] = start_ctbl as u16;
    ct.next_at_rank[k] = r;
    ct.len_pos = k as i32 + 1;
}

/// `CtPartINCHI_CANON_INFINITY`: terminates the connection table part of level `k`
pub fn ct_part_inchi_canon_infinity(ct: &mut ConTable, cmp: Option<&mut [i8]>, k: usize) {
    let k = k - 1;
    let start_ctbl = if k != 0 {
        if let Some(c) = cmp {
            c[..k].fill(0);
        }
        ct.next_ctbl_pos[k - 1] as usize
    } else {
        0
    };
    if start_ctbl == 0 || ct.ctbl[start_ctbl - 1] != EMPTY_CT {
        ct.ctbl[start_ctbl] = EMPTY_CT;
    }
}

/// `CtPartCompare`: compares the parts of two connection tables at level
/// `k` (= level + 1); returns `+-(layer + 1)` for the first difference or 0
pub fn ct_part_compare(
    ct1: &ConTable,
    ct2: &ConTable,
    mut cmp: Option<&mut [i8]>,
    k_least_for_layer: Option<&mut [KLeast]>,
    k: usize,
    b_only_common: bool,
    b_split_taut_compare: bool,
) -> i32 {
    let k = k - 1;
    let mut i: i32 = -1;
    if let Some(c) = cmp.as_deref_mut() {
        let mut ii = 0usize;
        while ii <= k && c[ii] == 0 {
            ii += 1;
        }
        i = ii as i32;
        if ii < k {
            c[k] = c[ii];
            return i32::from(c[ii]);
        }
    }
    let (mut start_ct1, start_ct2, mut start_at1, mut start_at2);
    if k != 0 {
        start_ct1 = i32::from(ct1.next_ctbl_pos[k - 1]);
        start_ct2 = i32::from(ct2.next_ctbl_pos[k - 1]);
        start_at1 = i32::from(ct1.next_at_rank[k - 1]) - 1;
        start_at2 = i32::from(ct2.next_at_rank[k - 1]) - 1;
    } else {
        start_ct1 = 0;
        start_ct2 = 0;
        start_at1 = 0;
        start_at2 = 0;
    }
    let mut end_ct1 = i32::from(ct1.next_ctbl_pos[k]);
    let mut end_ct2 = i32::from(ct2.next_ctbl_pos[k]);
    let mut end_at1 = i32::from(ct1.next_at_rank[k]) - 1;
    let mut end_at2 = i32::from(ct2.next_at_rank[k]) - 1;
    let max_vert = ct1.max_vert.min(ct2.max_vert);
    let mut n_layer = 0usize;
    let mut diff: i32;
    let mut mid_num_h = 0i32;
    'done: {
        // lengths
        diff = -(start_ct1 - start_ct2);
        if diff != 0 {
            if b_only_common
                && start_ct1 >= ct1.n_len_ct_at_only
                && start_ct2 >= ct2.n_len_ct_at_only
                && ct1.ctbl[start_ct1 as usize] == EMPTY_CT
                && ct2.ctbl[start_ct2 as usize] == EMPTY_CT
            {
                return 0;
            }
            if b_only_common {
                start_ct1 = start_ct1.min(start_ct2);
                start_at1 = start_at1.min(start_at2);
                start_at2 = start_at1;
                if ct1.len_ct == ct2.len_ct {
                    end_ct1 = end_ct1.max(end_ct2);
                    end_ct2 = end_ct1;
                    end_at1 = end_at1.max(end_at2);
                    end_at2 = end_at1;
                }
            } else if start_ct1 > start_ct2
                && ct1.max_vert > ct2.max_vert
                && start_at2 == ct2.max_vert
            {
                return 0;
            } else {
                break 'done;
            }
        }
        let mut len_num_h = ct1.len_num_h;
        diff = -(end_ct1 - end_ct2);
        if diff != 0 {
            if b_only_common {
                end_ct1 = end_ct1.min(end_ct2);
                end_at1 = end_at1.min(end_at2);
                end_at2 = end_at1;
                len_num_h = ct1.len_num_h.min(ct2.len_num_h);
            } else if end_ct1 > end_ct2 && ct1.maxlen_ct > ct2.maxlen_ct {
                if end_at2 == ct2.max_vert + 1 {
                    end_at2 -= 1;
                    end_at1 = end_at2;
                    len_num_h = end_at1;
                    end_ct2 -= 1;
                    end_ct1 = end_ct2;
                } else if end_at2 == ct2.max_vert {
                    end_at1 = end_at2;
                    len_num_h = end_at1;
                    end_ct1 = end_ct2;
                } else {
                    break 'done;
                }
            } else {
                break 'done;
            }
        }
        let (mid_ct, mid_at);
        if b_split_taut_compare {
            let mut m = ct1.n_len_ct_at_only.min(ct2.n_len_ct_at_only);
            if m > end_ct1 {
                m = end_ct1;
            }
            mid_ct = m;
            mid_at = max_vert.min(end_at1);
        } else {
            mid_ct = end_ct1;
            mid_at = end_at1;
        }
        // layer 0: the connection table without tautomeric groups
        i = start_ct1;
        while i < mid_ct && ct1.ctbl[i as usize] == ct2.ctbl[i as usize] {
            i += 1;
        }
        if i < mid_ct {
            diff = i32::from(ct1.ctbl[i as usize]) - i32::from(ct2.ctbl[i as usize]);
            break 'done;
        }
        // layer 1: H atoms without tautomeric H
        n_layer += 1;
        diff = -(start_at1 - start_at2);
        if diff != 0 {
            break 'done;
        }
        diff = -(end_at1 - end_at2);
        if diff != 0 {
            break 'done;
        }
        if !ct1.num_h.is_empty() && !ct2.num_h.is_empty() {
            if end_at1 < max_vert {
                len_num_h = end_at1;
                mid_num_h = end_at1;
            } else if b_split_taut_compare {
                mid_num_h = max_vert;
            } else {
                mid_num_h = len_num_h;
            }
            i = start_at1;
            while i < mid_num_h && ct1.num_h[i as usize] == ct2.num_h[i as usize] {
                i += 1;
            }
            if i < mid_num_h {
                diff = i32::from(ct1.num_h[i as usize]) - i32::from(ct2.num_h[i as usize]);
                break 'done;
            }
        }
        // layer 2: tautomeric part of the CT and tautomeric H
        n_layer += 1;
        i = mid_ct;
        while i < end_ct1 && ct1.ctbl[i as usize] == ct2.ctbl[i as usize] {
            i += 1;
        }
        if i < end_ct1 {
            diff = i32::from(ct1.ctbl[i as usize]) - i32::from(ct2.ctbl[i as usize]);
            break 'done;
        }
        if !ct1.num_h.is_empty() && !ct2.num_h.is_empty() {
            i = mid_num_h;
            while i < len_num_h && ct1.num_h[i as usize] == ct2.num_h[i as usize] {
                i += 1;
            }
            if i < len_num_h {
                diff = i32::from(ct1.num_h[i as usize]) - i32::from(ct2.num_h[i as usize]);
                i += end_ct1 - mid_ct;
                break 'done;
            }
        }
        // layer 3: fixed H atoms
        n_layer += 1;
        if !ct1.num_h_fixed.is_empty() && !ct2.num_h_fixed.is_empty() {
            i = start_at1;
            while i < mid_at && ct1.num_h_fixed[i as usize] == ct2.num_h_fixed[i as usize] {
                i += 1;
            }
            if i < mid_at {
                diff =
                    i32::from(ct1.num_h_fixed[i as usize]) - i32::from(ct2.num_h_fixed[i as usize]);
                break 'done;
            }
        }
        // layer 4: isotopic atoms, H and tautomeric groups
        n_layer += 1;
        if !ct1.iso_sort_key.is_empty() && !ct2.iso_sort_key.is_empty() {
            i = start_at1;
            while i < end_at1 && ct1.iso_sort_key[i as usize] == ct2.iso_sort_key[i as usize] {
                i += 1;
            }
            if i < end_at1 {
                diff = if ct1.iso_sort_key[i as usize] > ct2.iso_sort_key[i as usize] {
                    1
                } else {
                    -1
                };
                break 'done;
            }
        }
        // (sic) the C code tests `Ct2->len_iso_exchg_atnos` here
        if !ct1.iso_exchg_atnos.is_empty() && ct2.len_iso_exchg_atnos != 0 {
            i = start_at1;
            while i < end_at1 && ct1.iso_exchg_atnos[i as usize] == ct2.iso_exchg_atnos[i as usize]
            {
                i += 1;
            }
            if i < end_at1 {
                diff = if ct1.iso_exchg_atnos[i as usize] > ct2.iso_exchg_atnos[i as usize] {
                    1
                } else {
                    -1
                };
                break 'done;
            }
        }
    }
    // INCHI_CANON_MIN
    diff = -diff;
    if diff != 0 {
        diff = if diff > 0 {
            n_layer as i32 + 1
        } else {
            -(n_layer as i32 + 1)
        };
        if let Some(kl) = k_least_for_layer {
            if kl[n_layer].k == 0 {
                kl[n_layer].k = if diff > 0 {
                    k as i32 + 1
                } else {
                    -(k as i32 + 1)
                };
                kl[n_layer].i = i;
            }
            if n_layer != 0 {
                diff = 0;
            }
        }
    }
    if let Some(c) = cmp {
        c[k] = if diff > 0 {
            1
        } else if diff < 0 {
            -1
        } else {
            0
        };
    }
    diff
}

/// `CtFullCompare`: compares two complete connection tables
pub fn ct_full_compare(
    ct1: &ConTable,
    ct2: &ConTable,
    b_only_common: bool,
    b_split_taut_compare: bool,
) -> i32 {
    let k1 = (ct1.len_pos - 1) as usize;
    let k2 = (ct2.len_pos - 1) as usize;
    let start_ct1 = 0i32;
    let start_at1 = 0i32;
    let mut end_ct1 = i32::from(ct1.next_ctbl_pos[k1]);
    let mut end_ct2 = i32::from(ct2.next_ctbl_pos[k2]);
    let mut end_at1 = i32::from(ct1.next_at_rank[k1]) - 1;
    let end_at2 = i32::from(ct2.next_at_rank[k2]) - 1;
    let max_vert = ct1.max_vert.min(ct2.max_vert);
    let (len_num_h2, len_iso1, len_iso2): (i32, i32, i32);
    let mut len_num_h1: i32;
    let mut n_layer = 0;
    let mut diff: i32;
    let mut mid_num_h = 0i32;
    if b_only_common {
        end_ct1 = end_ct1.min(end_ct2);
        end_ct1 = end_ct1.min(ct1.len_ct);
        end_ct2 = end_ct1;
        end_at1 = end_at1.min(end_at2);
        if ct1.ctbl[end_ct1 as usize] == EMPTY_CT || ct2.ctbl[end_ct1 as usize] == EMPTY_CT {
            end_ct1 -= 1;
            end_ct2 = end_ct1;
        }
        len_num_h1 = ct1.len_num_h.min(ct2.len_num_h);
        len_num_h2 = len_num_h1;
        // (sic) the C code takes the minimum of `Ct1->len_iso_sort_key` twice
        len_iso1 = ct1.len_iso_sort_key.min(ct1.len_iso_sort_key);
        len_iso2 = len_iso1;
    } else {
        if ct1.ctbl[(end_ct1 - 1) as usize] == EMPTY_CT {
            end_ct1 -= 1;
        }
        if ct2.ctbl[(end_ct2 - 1) as usize] == EMPTY_CT {
            end_ct2 -= 1;
        }
        len_num_h1 = ct1.len_num_h;
        len_num_h2 = ct2.len_num_h;
        len_iso1 = ct1.len_iso_sort_key;
        len_iso2 = ct2.len_iso_sort_key;
    }
    'done: {
        diff = -(end_ct1 - end_ct2);
        if diff != 0 {
            break 'done;
        }
        let (mid_ct, _mid_at);
        if b_split_taut_compare {
            let mut m = ct1.n_len_ct_at_only.min(ct2.n_len_ct_at_only);
            if m > end_ct1 {
                m = end_ct1;
            }
            mid_ct = m;
            _mid_at = max_vert.min(end_at1);
        } else {
            mid_ct = end_ct1;
            _mid_at = end_at1;
        }
        let mut i = start_ct1;
        while i < mid_ct && ct1.ctbl[i as usize] == ct2.ctbl[i as usize] {
            i += 1;
        }
        if i < mid_ct {
            diff = i32::from(ct1.ctbl[i as usize]) - i32::from(ct2.ctbl[i as usize]);
            break 'done;
        }
        n_layer += 1;
        if !ct1.num_h.is_empty() && !ct2.num_h.is_empty() {
            diff = -(len_num_h1 - len_num_h2);
            if diff != 0 {
                break 'done;
            }
            if end_at1 < max_vert {
                len_num_h1 = end_at1;
                mid_num_h = end_at1;
            } else if b_split_taut_compare {
                mid_num_h = max_vert;
            } else {
                mid_num_h = len_num_h1;
            }
            i = start_at1;
            while i < mid_num_h && ct1.num_h[i as usize] == ct2.num_h[i as usize] {
                i += 1;
            }
            if i < mid_num_h {
                diff = i32::from(ct1.num_h[i as usize]) - i32::from(ct2.num_h[i as usize]);
                break 'done;
            }
        }
        n_layer += 1;
        i = mid_ct;
        while i < end_ct1 && ct1.ctbl[i as usize] == ct2.ctbl[i as usize] {
            i += 1;
        }
        if i < end_ct1 {
            diff = i32::from(ct1.ctbl[i as usize]) - i32::from(ct2.ctbl[i as usize]);
            break 'done;
        }
        if !ct1.num_h.is_empty() && !ct2.num_h.is_empty() {
            i = mid_num_h;
            while i < len_num_h1 && ct1.num_h[i as usize] == ct2.num_h[i as usize] {
                i += 1;
            }
            if i < len_num_h1 {
                diff = i32::from(ct1.num_h[i as usize]) - i32::from(ct2.num_h[i as usize]);
                break 'done;
            }
        }
        n_layer += 1;
        if !ct1.num_h_fixed.is_empty() && !ct2.num_h_fixed.is_empty() {
            i = start_at1;
            while i < end_at1 && ct1.num_h_fixed[i as usize] == ct2.num_h_fixed[i as usize] {
                i += 1;
            }
            if i < end_at1 {
                diff =
                    i32::from(ct1.num_h_fixed[i as usize]) - i32::from(ct2.num_h_fixed[i as usize]);
                break 'done;
            }
        }
        n_layer += 1;
        if !ct1.iso_sort_key.is_empty() && !ct2.iso_sort_key.is_empty() {
            diff = -(len_iso1 - len_iso2);
            if diff != 0 {
                break 'done;
            }
            i = start_at1;
            while i < end_at1 && ct1.iso_sort_key[i as usize] == ct2.iso_sort_key[i as usize] {
                i += 1;
            }
            if i < end_at1 {
                diff = if ct1.iso_sort_key[i as usize] > ct2.iso_sort_key[i as usize] {
                    1
                } else {
                    -1
                };
                break 'done;
            }
        }
    }
    diff = -diff;
    if diff != 0 {
        diff = if diff > 0 {
            n_layer + 1
        } else {
            -(n_layer + 1)
        };
    }
    diff
}

/// `CtFullCompareLayers`
pub fn ct_full_compare_layers(k_least_for_layer: &[KLeast]) -> i32 {
    for (i_layer, kl) in k_least_for_layer.iter().enumerate().take(MAX_LAYERS) {
        if kl.k != 0 {
            return if kl.k > 0 {
                i_layer as i32 + 1
            } else {
                -(i_layer as i32 + 1)
            };
        }
    }
    0
}

/// `CtCompareLayersGetFirstDiff`: returns 1 if a real layer differs, 0 if none
/// (with the extra layer in `l_rho`), -1 without input
pub fn ct_compare_layers_get_first_diff(
    k_least_rho: Option<&[KLeast]>,
    n_one_additional_layer: i32,
    l_rho: &mut i32,
    i_rho: &mut i32,
    k_rho: &mut i32,
) -> i32 {
    let Some(kl) = k_least_rho else {
        return -1;
    };
    for (i_layer, e) in kl.iter().enumerate().take(MAX_LAYERS) {
        if e.k != 0 {
            *l_rho = i_layer as i32;
            *i_rho = e.i;
            *k_rho = e.k;
            return 1;
        }
    }
    if n_one_additional_layer != 0 {
        *l_rho = n_one_additional_layer;
        *i_rho = -1;
        *k_rho = 0;
        0
    } else {
        *l_rho = INCHI_CANON_INFINITY as i32;
        *i_rho = -1;
        *k_rho = 0;
        0
    }
}

/// `CtPartCompareLayers`
pub fn ct_part_compare_layers(
    k_least_rho: Option<&[KLeast]>,
    l_rho_fix_prev: i32,
    n_one_additional_layer: i32,
) -> i32 {
    let (mut l, mut i, mut k) = (0, 0, 0);
    if 0 < ct_compare_layers_get_first_diff(
        k_least_rho,
        n_one_additional_layer,
        &mut l,
        &mut i,
        &mut k,
    ) && l <= l_rho_fix_prev
    {
        return if k > 0 { l + 1 } else { -(l + 1) };
    }
    0
}

/// `UpdateCompareLayers`
pub fn update_compare_layers(k_least_for_layer: Option<&mut [KLeast]>, hzz: i32) {
    if let Some(kl) = k_least_for_layer {
        for e in kl.iter_mut().take(MAX_LAYERS) {
            if e.k.abs() >= hzz {
                e.k = 0;
                e.i = 0;
            }
        }
    }
}

/// `CtPartCopy`: copies level `k` (= level + 1) of `ct2` into `ct1`
pub fn ct_part_copy(ct1: &mut ConTable, ct2: &ConTable, k: usize) {
    let k = k - 1;
    let (start_ct1, start_ct2, start_at1, start_at2) = if k != 0 {
        (
            ct1.next_ctbl_pos[k - 1] as usize,
            ct2.next_ctbl_pos[k - 1] as usize,
            ct1.next_at_rank[k - 1] as i32 - 1,
            ct2.next_at_rank[k - 1] as i32 - 1,
        )
    } else {
        (0, 0, 0, 0)
    };
    let end_ct2 = ct2.next_ctbl_pos[k] as usize;
    let end_at2 = i32::from(ct2.next_at_rank[k]) - 1;
    let len2 = end_ct2 - start_ct2;
    for i in 0..len2 {
        ct1.ctbl[start_ct1 + i] = ct2.ctbl[start_ct2 + i];
    }
    let mut len2_h = 0i32;
    if !ct1.num_h.is_empty() && !ct2.num_h.is_empty() {
        len2_h = end_at2 - start_at2;
        if end_at2 > ct2.max_vert {
            len2_h = ct2.len_num_h - start_at2;
        }
        for i in 0..len2_h {
            ct1.num_h[(start_at1 + i) as usize] = ct2.num_h[(start_at2 + i) as usize];
        }
    }
    if !ct1.num_h_fixed.is_empty() && !ct2.num_h_fixed.is_empty() {
        for i in 0..(end_at2 - start_at2) {
            ct1.num_h_fixed[(start_at1 + i) as usize] = ct2.num_h_fixed[(start_at2 + i) as usize];
        }
    }
    let mut len2_iso = 0i32;
    if !ct1.iso_sort_key.is_empty() && !ct2.iso_sort_key.is_empty() {
        len2_iso = end_at2 - start_at2;
        for i in 0..len2_iso {
            ct1.iso_sort_key[(start_at1 + i) as usize] = ct2.iso_sort_key[(start_at2 + i) as usize];
        }
    }
    let mut len2_exchg = 0i32;
    if !ct1.iso_exchg_atnos.is_empty() && !ct2.iso_exchg_atnos.is_empty() {
        len2_exchg = end_at2 - start_at2;
        for i in 0..len2_exchg {
            ct1.iso_exchg_atnos[(start_at1 + i) as usize] =
                ct2.iso_exchg_atnos[(start_at2 + i) as usize];
        }
    }
    ct1.len_ct = (start_ct1 + len2) as i32;
    ct1.next_ctbl_pos[k] = (start_ct1 + len2) as u16;
    ct1.next_at_rank[k] = ct2.next_at_rank[k];
    if len2_h != 0 {
        ct1.len_num_h = start_at1 + len2_h;
    }
    if len2_iso != 0 {
        ct1.len_iso_sort_key = start_at1 + len2_iso;
    }
    if len2_exchg != 0 {
        ct1.len_iso_exchg_atnos = start_at1 + len2_exchg;
    }
    ct1.len_pos = k as i32 + 1;
}

/// `CtFullCopy`
pub fn ct_full_copy(ct1: &mut ConTable, ct2: &ConTable) {
    for k in 0..ct2.len_pos as usize {
        ct_part_copy(ct1, ct2, k + 1);
    }
}

/// `TranspositionGetMcrAndFixSetAndUnorderedPartition`: the orbit minimal
/// representatives, fixed points and orbit partition of the permutation `gamma`
pub fn transposition_get_mcr_and_fix_set_and_unordered_partition(
    gamma: &mut Transposition,
    mcr_set: &mut NodeSet,
    fix_set: &mut NodeSet,
    n: usize,
    l: usize,
    p: &mut UnorderedPartition,
) {
    mcr_set.bitword[l - 1].iter_mut().for_each(|b| *b = 0);
    fix_set.bitword[l - 1].iter_mut().for_each(|b| *b = 0);
    for e in p.equ2[..n].iter_mut() {
        *e = INCHI_CANON_INFINITY;
    }
    for i in 0..n {
        let next = gamma.n_at_numb[i];
        let mut j = next as usize;
        if j == i {
            fix_set.bitword[l - 1][i / NUM_BIT] |= 1 << (i % NUM_BIT);
            mcr_set.bitword[l - 1][i / NUM_BIT] |= 1 << (i % NUM_BIT);
            p.equ2[i] = next;
        } else if RANK_MARK_BIT & next == 0 {
            gamma.n_at_numb[i] |= RANK_MARK_BIT;
            let mut mcr = j.min(i);
            loop {
                let nx = gamma.n_at_numb[j];
                if RANK_MARK_BIT & nx != 0 {
                    break;
                }
                gamma.n_at_numb[j] |= RANK_MARK_BIT;
                j = nx as usize;
                if mcr > j {
                    mcr = j;
                }
            }
            mcr_set.bitword[l - 1][mcr / NUM_BIT] |= 1 << (mcr % NUM_BIT);
            p.equ2[mcr] = mcr as u16;
            let mut k = mcr;
            loop {
                j = (RANK_MASK_BIT & gamma.n_at_numb[k]) as usize;
                if mcr == j {
                    break;
                }
                p.equ2[j] = mcr as u16;
                k = j;
            }
        }
    }
    for e in gamma.n_at_numb[..n].iter_mut() {
        *e &= RANK_MASK_BIT;
    }
}

/// `GetOneAdditionalLayer`: the number of the only layer that `pCD` has and
/// `pzb_rho_fix` does not, or 0
pub fn get_one_additional_layer(cd: Option<&CanonData>, pzb_rho_fix: Option<&ConTable>) -> i32 {
    let (Some(cd), Some(z)) = (cd, pzb_rho_fix) else {
        return 0;
    };
    let (mut last, mut num_last, mut layer) = (-1, 0, 0);
    layer += 1;
    if cd.num_h.is_some() && z.num_h.is_empty() {
        last = layer;
        num_last += 1;
    }
    layer += 1;
    if cd.n_len_ct_at_only < cd.n_len_linear_ct && z.n_len_ct_at_only == z.len_ct {
        last = layer;
        num_last += 1;
    }
    layer += 1;
    if cd.num_h_fixed.is_some() && z.num_h_fixed.is_empty() {
        last = layer;
        num_last += 1;
    }
    layer += 1;
    if cd.iso_sort_key.is_some() && z.iso_sort_key.is_empty() {
        last = layer;
        num_last += 1;
    }
    if num_last == 1 { last } else { 0 }
}

// ---------------------------------------------------------------------------
// initial ranks and atom invariants, clean-up of the canonical CT arrays
// ---------------------------------------------------------------------------

use super::ichimap::{
    AT_INV_HILL_ORDER, AT_INV_NUM_CONNECTIONS, AT_INV_NUM_H, AT_INV_NUM_H_FIX,
    AT_INV_NUM_TG_ENDPOINTS, AT_INV_TAUT_ISO, AT_INV_TG_NUMBERS, AtomInvariant2,
    comp_atom_invariants2, comp_atom_invariants2_only,
};
use super::ichisort::{comp_chem_elem_lex, inchi_qsort};
use super::ichister::SpAtom;
use super::ichitaut::{T_NUM_ISOTOPIC, TGroupInfo};
use super::inpdef::MAXVAL;

pub const EMPTY_H_NUMBER: i16 = INCHI_CANON_INFINITY as i16 - 1;
pub const BASE_H_NUMBER: i16 = (INCHI_CANON_INFINITY as i16 - 1) / 2;
pub const EMPTY_ISO_SORT_KEY: i64 = i64::MAX;
const ELEM_NAME_LEN: usize = 2;
const NUM_CHEM_ELEMENTS: usize = 127;

/// `SetInitialRanks2`: sorts the atoms by their invariants and assigns the
/// initial ranks; returns the number of different ranks
pub fn set_initial_ranks2(
    num_atoms: usize,
    atom_invariant2: &[AtomInvariant2],
    n_new_rank: &mut [u16],
    n_atom_number: &mut [u16],
) -> i32 {
    for (i, a) in n_atom_number[..num_atoms].iter_mut().enumerate() {
        *a = i as u16;
    }
    inchi_qsort(&mut n_atom_number[..num_atoms], |a, b| {
        comp_atom_invariants2(atom_invariant2, *a, *b)
    });
    let mut n_num_diff_ranks = 1;
    if num_atoms > 0 {
        let mut n_current_rank = num_atoms as u16;
        n_new_rank[n_atom_number[num_atoms - 1] as usize] = n_current_rank;
        for i in (1..num_atoms).rev() {
            if comp_atom_invariants2_only(
                &atom_invariant2[n_atom_number[i - 1] as usize],
                &atom_invariant2[n_atom_number[i] as usize],
            ) != 0
            {
                n_num_diff_ranks += 1;
                n_current_rank = i as u16;
            }
            n_new_rank[n_atom_number[i - 1] as usize] = n_current_rank;
        }
    }
    n_num_diff_ranks
}

/// the part of a NUL terminated byte string before the NUL
fn c_str(b: &[u8]) -> &[u8] {
    let n = b.iter().position(|&c| c == 0).unwrap_or(b.len());
    &b[..n]
}

/// `strstr`: the position of `needle` in the C string `hay`
fn c_strstr(hay: &[u8], needle: &[u8]) -> Option<usize> {
    let hay = c_str(hay);
    let needle = c_str(needle);
    if needle.is_empty() {
        return Some(0);
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

/// `FillOutAtomInvariant2`
#[allow(clippy::too_many_arguments)]
pub fn fill_out_atom_invariant2(
    at: &[SpAtom],
    num_atoms: usize,
    num_at_tg: usize,
    atom_invariant: &mut [AtomInvariant2],
    b_ignore_isotopic: bool,
    b_hydrogens_in_ranks: bool,
    b_hydrogens_fixed_in_ranks: bool,
    b_digraph: bool,
    b_taut_groups_only: bool,
    t_group_info: Option<&TGroupInfo>,
) {
    let mut chem_elements = [0u8; ELEM_NAME_LEN * NUM_CHEM_ELEMENTS + ELEM_NAME_LEN];
    let mut num_chem_elements = 0usize;
    let mut num_hydrogen_atoms = 0usize;
    let mut num_carbon_atoms = 0usize;
    let mut t_group: &[super::ichitaut::TGroup] = &[];
    let mut have_t_group = false;
    let mut num_t_groups = 0usize;
    let mut num_tautomer_iso = 0usize;
    if num_at_tg > num_atoms
        && let Some(tgi) = t_group_info
    {
        if let Some(tg) = tgi.t_group.as_deref() {
            t_group = tg;
            have_t_group = true;
        }
        num_t_groups = tgi.num_t_groups;
        num_tautomer_iso = if tgi.b_ignore_isotopic != 0 {
            0
        } else {
            T_NUM_ISOTOPIC
        };
    }
    if !b_taut_groups_only {
        for a in &at[..num_atoms] {
            let el = c_str(&a.elname);
            if el == b"C" {
                num_carbon_atoms += 1;
            } else if el == b"H" || el == b"D" || el == b"T" {
                num_hydrogen_atoms += 1;
            } else {
                let cur = [
                    a.elname[0],
                    if a.elname[1] != 0 { a.elname[1] } else { b' ' },
                ];
                if c_strstr(&chem_elements, &cur).is_none() {
                    // strcat
                    let len = c_str(&chem_elements).len();
                    chem_elements[len] = cur[0];
                    chem_elements[len + 1] = cur[1];
                    num_chem_elements += 1;
                }
            }
        }
        if num_chem_elements > 1 {
            let mut recs: Vec<[u8; 2]> = (0..num_chem_elements)
                .map(|i| [chem_elements[2 * i], chem_elements[2 * i + 1]])
                .collect();
            recs.sort_by(|a, b| comp_chem_elem_lex(a, b).cmp(&0));
            for (i, r) in recs.iter().enumerate() {
                chem_elements[2 * i] = r[0];
                chem_elements[2 * i + 1] = r[1];
            }
        }
        if num_carbon_atoms != 0 {
            if num_chem_elements != 0 {
                chem_elements.copy_within(0..ELEM_NAME_LEN * num_chem_elements, ELEM_NAME_LEN);
            }
            chem_elements[0] = b'C';
            chem_elements[1] = b' ';
            num_chem_elements += 1;
        }
        if num_hydrogen_atoms != 0 {
            chem_elements[ELEM_NAME_LEN * num_chem_elements] = b'H';
            chem_elements[ELEM_NAME_LEN * num_chem_elements + 1] = b' ';
            num_chem_elements += 1;
        }
        for i in 0..num_atoms {
            atom_invariant[i] = AtomInvariant2::default();
            let cur = [
                at[i].elname[0],
                if at[i].elname[1] != 0 {
                    at[i].elname[1]
                } else {
                    b' '
                },
            ];
            let j = match c_strstr(&chem_elements, &cur) {
                Some(p) => p / ELEM_NAME_LEN + 1,
                None => num_chem_elements,
            };
            atom_invariant[i].val[AT_INV_HILL_ORDER] = j as u16;
            atom_invariant[i].val[AT_INV_NUM_CONNECTIONS] = at[i].valence as u16;
            if b_hydrogens_in_ranks {
                atom_invariant[i].val[AT_INV_NUM_H] = if have_t_group && at[i].endpoint > 0 {
                    0
                } else {
                    at[i].num_h as u16
                };
            }
            if b_hydrogens_fixed_in_ranks {
                atom_invariant[i].val[AT_INV_NUM_H_FIX] = if have_t_group && at[i].endpoint > 0 {
                    at[i].num_h as u16
                } else {
                    0
                };
            }
            let i_t_group = i32::from(at[i].endpoint) - 1;
            if !b_digraph && have_t_group && i_t_group >= 0 && (i_t_group as usize) < num_t_groups {
                let g = &t_group[i_t_group as usize];
                atom_invariant[i].val[AT_INV_NUM_TG_ENDPOINTS] = g.n_num_endpoints;
                for j in 0..T_NUM_NO_ISOTOPIC {
                    atom_invariant[i].val[AT_INV_TG_NUMBERS + j] = g.num[j];
                }
                for j in 0..num_tautomer_iso {
                    atom_invariant[i].val[AT_INV_TAUT_ISO + j] = g.num[j + T_NUM_NO_ISOTOPIC];
                }
            }
            atom_invariant[i].iso_sort_key = if b_ignore_isotopic {
                0
            } else {
                at[i].iso_sort_key
            };
        }
    } else {
        for a in atom_invariant[..num_at_tg].iter_mut() {
            *a = AtomInvariant2::default();
        }
    }
    // tautomeric groups
    for i in num_atoms..num_at_tg {
        let k = i - num_atoms;
        atom_invariant[i] = AtomInvariant2::default();
        if !have_t_group {
            continue;
        }
        atom_invariant[i].val[AT_INV_HILL_ORDER] = if b_taut_groups_only {
            num_at_tg as u16
        } else {
            num_chem_elements as u16 + 1
        };
        atom_invariant[i].val[AT_INV_NUM_CONNECTIONS] = MAXVAL as u16 + 1;
        if k < num_t_groups {
            let g = &t_group[k];
            atom_invariant[i].val[AT_INV_NUM_TG_ENDPOINTS] = g.n_num_endpoints;
            for j in 0..T_NUM_NO_ISOTOPIC {
                atom_invariant[i].val[AT_INV_TAUT_ISO + j] = g.num[j];
            }
            for j in 0..num_tautomer_iso {
                atom_invariant[i].val[AT_INV_TAUT_ISO + j] = g.num[j + T_NUM_NO_ISOTOPIC];
            }
        }
    }
}

/// `CleanNumH`
pub fn clean_num_h(num_h: &mut [i16], len: usize) {
    for h in num_h[..len].iter_mut() {
        if *h == EMPTY_H_NUMBER {
            *h = 0;
        } else {
            *h = h.wrapping_sub(BASE_H_NUMBER);
        }
    }
}

/// `CleanCt`: replaces the terminating `EMPTY_CT` at `ct[len]` by 0
pub fn clean_ct(ct: &mut [u16], len: usize) -> i32 {
    if ct[len] == EMPTY_CT {
        ct[len] = 0;
        1
    } else {
        0
    }
}

/// `CleanIsoSortKeys`
pub fn clean_iso_sort_keys(isk: &mut [i64], len: usize) {
    for k in isk[..len].iter_mut() {
        if *k == EMPTY_ISO_SORT_KEY {
            *k = 0;
        }
    }
}

// ---------------------------------------------------------------------------
// CanonGraph
// ---------------------------------------------------------------------------

pub const MAX_SET_SIZE: usize = 32766;
pub const NORMALLY_ALLOWED_MAX_SET_SIZE: usize = 2048;

/// `CANON_COUNTS`
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct CanonCounts {
    pub l_num_break_ties: i64,
    pub l_num_decreased_ct: i64,
    pub l_num_rejected_ct: i64,
    pub l_num_equal_ct: i64,
    pub l_num_tot_ct: i64,
    pub d_group_size: f64,
    pub l_num_generators: i64,
    pub l_num_stored_isomorphisms: i64,
}

/// `QZFIX_OK`
fn qzfix_ok(x: i32) -> bool {
    x <= 0
}

/// the labels of the `goto` structure of `CanonGraph`
#[derive(Clone, Copy, PartialEq, Eq)]
enum Lb {
    L2,
    L6,
    L7,
    L8,
    L9,
    L10,
    L11,
    L12,
    L13,
    L14,
    L15,
    L16,
    L17,
}

/// `CanonGraph`: the canonical numbering search (McKay's algorithm).
/// `pzb_rho_inp` is `*pp_zb_rho_inp`; `pzb_rho_out` is `pp_zb_rho_out`.
/// `time_is_over` replaces `bInchiTimeIsOver`. Returns 0 or the last number
/// of ranks (as the C code) or a negative error code.
#[allow(clippy::too_many_arguments, unused_assignments, unused_variables)]
pub fn canon_graph(
    time_is_over: &mut dyn FnMut() -> bool,
    n: usize,
    n_tg: usize,
    n_max: usize,
    b_digraph: bool,
    g: &mut NeighList,
    pi: &mut [Partition],
    n_symm_rank: &mut [u16],
    n_canon_rank: &mut [u16],
    n_atom_number_canon: &mut [u16],
    cd: &mut CanonData,
    cc: &mut CanonCounts,
    pzb_rho_inp: Option<&ConTable>,
    pzb_rho_out: Option<&mut Option<ConTable>>,
    large_molecules: bool,
) -> i32 {
    let n_num_layers = i32::from(cd.num_h.is_some())
        + i32::from(cd.num_h_fixed.is_some())
        + i32::from(cd.iso_sort_key.is_some());
    let dig = b_digraph || n_num_layers != 0;
    let b_split_taut_compare = dig;
    let lab = true;
    let l_curr_max_set_size = if large_molecules {
        MAX_SET_SIZE
    } else {
        NORMALLY_ALLOWED_MAX_SET_SIZE
    } as i32;
    let pzb_rho_fix = pzb_rho_inp;

    let n_one_additional_layer = get_one_additional_layer(Some(cd), pzb_rho_fix);
    if let Some(z) = pzb_rho_fix
        && z.n_len_ct_at_only != cd.n_len_ct_at_only
    {
        return -2;
    }
    let create = |cd: &CanonData| {
        c_table_create(
            n,
            cd.n_max_len_linear_ct,
            cd.n_len_ct_at_only,
            cd.num_h.as_ref().map(|_| cd.maxlen_num_h),
            cd.num_h_fixed.as_ref().map(|_| cd.maxlen_num_h_fixed),
            cd.iso_sort_key.as_ref().map(|_| cd.maxlen_iso_sort_key),
            cd.iso_exchg_atnos
                .as_ref()
                .map(|_| cd.maxlen_iso_exchg_atnos),
        )
    };
    let mut theta = UnorderedPartition {
        equ2: vec![0; n_tg],
    };
    let mut theta_from_gamma = UnorderedPartition {
        equ2: vec![0; n_tg],
    };
    let mut w = vec![Cell::default(); n_tg];
    let mut v: Vec<Node> = vec![0; n_tg];
    let mut e: Vec<i8> = vec![0; n_tg];
    let mut lambda = create(cd);
    let mut zf_zeta = create(cd);
    let mut pzb_rho = create(cd);
    let mut omega = node_set_create(n_tg, l_curr_max_set_size as usize);
    let mut phi = node_set_create(n_tg, l_curr_max_set_size as usize);
    let mut cur_nodes = node_set_create(n_tg, 1);
    let mut zeta = Partition {
        rank: vec![0; n_tg],
        at_number: vec![0; n_tg],
    };
    let mut rho = zeta.clone();
    let mut gamma = Transposition {
        n_at_numb: vec![0; n_tg],
    };

    let mut ret = 0i32;
    let mut k: i32 = 1;
    let mut size = 1.0f64;
    let mut hz_rho = 0i32;
    let mut index = 0i32;
    let mut l = 0i32;
    unordered_partition_make_discrete(&mut theta, n_tg);
    let mut t_lemma = 2i32;
    cc.l_num_break_ties = 0;
    cc.l_num_decreased_ct = 0;
    cc.l_num_rejected_ct = 0;
    cc.l_num_equal_ct = 1;
    cc.l_num_tot_ct = 0;

    let mut hzb_rho_fix = 1i32;
    let mut k_least_rho = vec![KLeast::default(); MAX_LAYERS];
    let mut k_least_rho_fix = vec![KLeast::default(); MAX_LAYERS];
    let mut qzb_rho_fix = 0i32;
    let mut qzb_rho = 0i32;
    let mut b_rho_is_discrete = true;
    let mut pzb_rho_fix_reached = false;
    let (mut l_rho_fix_prev, mut i_rho_fix_prev, mut k_rho_fix_prev) = (0i32, -1i32, 0i32);
    let mut n_num_found_generators = 0i64;
    let mut b_zeta_eq_rho = lab;
    let mut b_zeta_isomorph = false;
    let (mut h_zeta, mut t_eq_zeta, mut hz_zeta, mut h_rho, mut r) = (0i32, 0i32, 0i32, 0i32, 0i32);
    let (mut tvc, mut tvh): (Node, Node) = (0, 0);

    macro_rules! ku {
        () => {
            (k - 1) as usize
        };
    }
    macro_rules! aux {
        () => {
            cd.n_aux_rank.as_deref()
        };
    }

    let mut skip_search = false;
    if partition_is_discrete(&pi[ku!()], n_tg) {
        partition_copy(&mut rho, &pi[ku!()], n_tg);
        ct_part_fill(g, cd, &pi[ku!()], &mut pzb_rho, 1, n, n_tg, n_max);
        ct_part_inchi_canon_infinity(&mut pzb_rho, None, 2);
        cc.l_num_tot_ct += 1;
        skip_search = true;
    }
    if !skip_search {
        if !dig && partition_satisfies_lemma_2_25(&pi[0], n) {
            t_lemma = 1;
        }
        ct_part_clear(&mut lambda, 1);
        // L2: reach the first leaf and save it in zeta and rho
        while k != 0 {
            partition_get_first_cell(&pi[ku!()], &mut w, k as usize, n);
            v[ku!()] = cell_get_min_node(&pi[ku!()], &w[ku!()], 0, aux!());
            e[ku!()] = 0;
            if dig || !partition_satisfies_lemma_2_25(&pi[ku!()], n) {
                t_lemma = k + 1;
            }
            ret =
                partition_color_vertex(g, &mut pi[ku!()..], v[ku!()], n, n_tg, n_max, b_digraph, 0);
            if ret < 0 {
                return ret;
            }
            cc.l_num_break_ties += 1;
            k += 1;
            ct_part_fill(
                g,
                cd,
                &pi[ku!()],
                &mut lambda,
                (k - 1) as usize,
                n,
                n_tg,
                n_max,
            );
            if let Some(z) = pzb_rho_fix
                && qzfix_ok(qzb_rho_fix)
            {
                qzb_rho_fix = ct_part_compare(
                    &lambda,
                    z,
                    None,
                    Some(&mut k_least_rho_fix),
                    (k - 1) as usize,
                    true,
                    b_split_taut_compare,
                );
                if qzfix_ok(qzb_rho_fix) {
                    hzb_rho_fix = k;
                }
            }
            if lab && qzfix_ok(qzb_rho_fix) {
                ct_part_copy(&mut pzb_rho, &lambda, (k - 1) as usize);
            }
            ct_part_copy(&mut zf_zeta, &lambda, (k - 1) as usize);
            if partition_is_discrete(&pi[ku!()], n) {
                break;
            }
        }
        cc.l_num_tot_ct += 1;
        h_zeta = k;
        t_eq_zeta = k;
        hz_zeta = k;
        ct_part_inchi_canon_infinity(&mut zf_zeta, None, k as usize);
        partition_copy(&mut zeta, &pi[ku!()], n_tg);
        if let Some(z) = pzb_rho_fix {
            if qzb_rho_fix == 0 {
                qzb_rho_fix = ct_full_compare(&lambda, z, true, b_split_taut_compare);
                if qzb_rho_fix > 0 {
                    hzb_rho_fix = 1;
                }
            }
            if hzb_rho_fix > 1 {
                partition_copy(&mut rho, &pi[(hzb_rho_fix - 1) as usize], n_tg);
            }
            hz_rho = hzb_rho_fix;
            h_rho = hzb_rho_fix;
            b_rho_is_discrete = hzb_rho_fix == k;
            if b_rho_is_discrete {
                ct_part_inchi_canon_infinity(&mut pzb_rho, None, k as usize);
                pzb_rho_fix_reached = qzb_rho_fix == 0;
                ct_compare_layers_get_first_diff(
                    Some(&k_least_rho_fix),
                    n_one_additional_layer,
                    &mut l_rho_fix_prev,
                    &mut i_rho_fix_prev,
                    &mut k_rho_fix_prev,
                );
            }
        } else {
            partition_copy(&mut rho, &pi[ku!()], n_tg);
            hz_rho = k;
            h_rho = k;
            ct_part_inchi_canon_infinity(&mut pzb_rho, None, k as usize);
        }
        qzb_rho = 0;
        r = k;
        v[ku!()] = INCHI_CANON_INFINITY;
        cell_make_empty(&mut w, k as usize);
        k -= 1;

        let mut lb = Lb::L13;
        'sm: loop {
            match lb {
                Lb::L2 => {
                    ret = partition_color_vertex(
                        g,
                        &mut pi[ku!()..],
                        v[ku!()],
                        n,
                        n_tg,
                        n_max,
                        b_digraph,
                        0,
                    );
                    if ret < 0 {
                        return ret;
                    }
                    cc.l_num_break_ties += 1;
                    k += 1;
                    ct_part_fill(
                        g,
                        cd,
                        &pi[ku!()],
                        &mut lambda,
                        (k - 1) as usize,
                        n,
                        n_tg,
                        n_max,
                    );
                    e[ku!()] = 0;
                    v[ku!()] = INCHI_CANON_INFINITY;
                    cell_make_empty(&mut w, k as usize);
                    if hz_zeta == k - 1
                        && ct_part_compare(
                            &lambda,
                            &zf_zeta,
                            None,
                            None,
                            (k - 1) as usize,
                            false,
                            b_split_taut_compare,
                        ) == 0
                    {
                        hz_zeta = k;
                    }
                    if let Some(z) = pzb_rho_fix
                        && qzb_rho_fix == 0
                    {
                        qzb_rho_fix = ct_part_compare(
                            &lambda,
                            z,
                            None,
                            Some(&mut k_least_rho_fix),
                            (k - 1) as usize,
                            true,
                            b_split_taut_compare,
                        );
                        if qzb_rho_fix == 0 && b_rho_is_discrete {
                            qzb_rho_fix = ct_part_compare_layers(
                                Some(&k_least_rho_fix),
                                l_rho_fix_prev,
                                n_one_additional_layer,
                            );
                            if qzb_rho_fix != 0 {
                                let d = qzb_rho_fix.abs() - 1;
                                // (sic) C compares with the result of `&&`
                                if d < l_rho_fix_prev
                                    || d == i32::from(
                                        l_rho_fix_prev != 0
                                            && k_least_rho_fix[d as usize].i < i_rho_fix_prev,
                                    )
                                {
                                    qzb_rho_fix = d + 1;
                                }
                            }
                        }
                        if !qzfix_ok(qzb_rho_fix) {
                            cc.l_num_rejected_ct += 1;
                        }
                    }
                    if pzb_rho_fix.is_some() && qzfix_ok(qzb_rho_fix) {
                        hzb_rho_fix = k;
                    }
                    if lab && qzfix_ok(qzb_rho_fix) {
                        if hz_rho == k - 1 && qzb_rho == 0 && b_rho_is_discrete {
                            qzb_rho = ct_part_compare(
                                &lambda,
                                &pzb_rho,
                                None,
                                Some(&mut k_least_rho),
                                (k - 1) as usize,
                                false,
                                b_split_taut_compare,
                            );
                            if qzb_rho == 0 && b_rho_is_discrete {
                                qzb_rho =
                                    ct_part_compare_layers(Some(&k_least_rho), l_rho_fix_prev, 0);
                                if qzb_rho != 0 {
                                    let d = qzb_rho.abs() - 1;
                                    if d < l_rho_fix_prev
                                        || d == i32::from(
                                            l_rho_fix_prev != 0
                                                && k_least_rho[d as usize].i < i_rho_fix_prev,
                                        )
                                    {
                                        qzb_rho = -(d + 1);
                                    }
                                }
                            }
                            if qzb_rho == 0 {
                                hz_rho = k;
                            } else if qzb_rho < 0 {
                                cc.l_num_rejected_ct += 1;
                            }
                        }
                        if (qzb_rho > 0 || (qzb_rho == 0 && !b_rho_is_discrete))
                            && n_num_layers == 0
                        {
                            ct_part_copy(&mut pzb_rho, &lambda, (k - 1) as usize);
                        }
                    }
                    // L3
                    if hz_zeta == k || hz_rho == k || (lab && qzb_rho >= 0 && qzfix_ok(qzb_rho_fix))
                    {
                        // L4
                        if partition_is_discrete(&pi[ku!()], n) {
                            cc.l_num_tot_ct += 1;
                            lb = Lb::L7;
                            continue 'sm;
                        }
                        partition_get_first_cell(&pi[ku!()], &mut w, k as usize, n);
                        v[ku!()] = cell_get_min_node(&pi[ku!()], &w[ku!()], 0, aux!());
                        if dig || !partition_satisfies_lemma_2_25(&pi[ku!()], n) {
                            t_lemma = k + 1;
                        }
                        e[ku!()] = 0;
                        lb = Lb::L2;
                        continue 'sm;
                    }
                    lb = Lb::L6;
                }
                Lb::L6 => {
                    let k2 = k;
                    k = (t_lemma - 1).min((t_eq_zeta - 1).max(hz_rho));
                    if k2 == t_lemma {
                        lb = Lb::L13;
                        continue 'sm;
                    }
                    if dig {
                        lb = Lb::L13;
                        continue 'sm;
                    }
                    l = (l + 1).min(l_curr_max_set_size);
                    partition_get_mcr_and_fix_set(
                        &pi[(t_lemma - 1) as usize],
                        &mut omega,
                        &mut phi,
                        n_tg,
                        l as usize,
                    );
                    lb = Lb::L12;
                }
                Lb::L7 => {
                    if h_zeta == 0 {
                        return CT_CANON_ERR;
                    }
                    if k != hz_zeta {
                        lb = Lb::L8;
                        continue 'sm;
                    }
                    if ct_full_compare(&lambda, &zf_zeta, false, b_split_taut_compare) == 0 {
                        partition_get_transposition(&zeta, &pi[ku!()], n_tg, &mut gamma);
                        b_zeta_isomorph = true;
                        lb = Lb::L10;
                        continue 'sm;
                    } else if n_num_layers == 0 {
                        return -2;
                    }
                    lb = Lb::L8;
                }
                Lb::L8 => {
                    if !lab || (qzb_rho < 0 && (pzb_rho_fix.is_none() || qzb_rho_fix > 0)) {
                        lb = Lb::L6;
                        continue 'sm;
                    }
                    if let Some(z) = pzb_rho_fix
                        && qzb_rho_fix == 0
                    {
                        qzb_rho_fix = ct_full_compare(&lambda, z, true, b_split_taut_compare);
                        if !pzb_rho_fix_reached {
                            pzb_rho_fix_reached = qzb_rho_fix == 0;
                        }
                        if 0 < qzb_rho_fix {
                            qzb_rho_fix = 0;
                            lb = Lb::L6;
                            continue 'sm;
                        }
                        qzb_rho_fix = 0;
                    }
                    if qzb_rho < 0 {
                        lb = Lb::L6;
                        continue 'sm;
                    }
                    if qzb_rho > 0 || !b_rho_is_discrete {
                        lb = Lb::L9;
                        continue 'sm;
                    }
                    if k < r {
                        lb = Lb::L9;
                        continue 'sm;
                    }
                    qzb_rho = ct_full_compare(&lambda, &pzb_rho, false, b_split_taut_compare);
                    if 0 < qzb_rho {
                        lb = Lb::L9;
                        continue 'sm;
                    }
                    if 0 > qzb_rho {
                        qzb_rho = 0;
                        lb = Lb::L6;
                        continue 'sm;
                    }
                    if r != k {
                        r = k;
                    }
                    partition_get_transposition(&pi[ku!()], &rho, n_tg, &mut gamma);
                    b_zeta_isomorph = false;
                    cc.l_num_equal_ct += 1;
                    lb = Lb::L10;
                }
                Lb::L9 => {
                    partition_copy(&mut rho, &pi[ku!()], n_tg);
                    if n_num_layers != 0 {
                        ct_full_copy(&mut pzb_rho, &lambda);
                    }
                    b_zeta_eq_rho = false;
                    qzb_rho = 0;
                    ct_compare_layers_get_first_diff(
                        Some(&k_least_rho_fix),
                        n_one_additional_layer,
                        &mut l_rho_fix_prev,
                        &mut i_rho_fix_prev,
                        &mut k_rho_fix_prev,
                    );
                    k_least_rho.iter_mut().for_each(|x| *x = KLeast::default());
                    h_rho = k;
                    hz_rho = k;
                    ct_part_inchi_canon_infinity(&mut pzb_rho, None, k as usize);
                    cc.l_num_decreased_ct += 1;
                    cc.l_num_equal_ct = 1;
                    b_rho_is_discrete = true;
                    lb = Lb::L6;
                }
                Lb::L10 => {
                    cc.l_num_equal_ct +=
                        i64::from(b_zeta_eq_rho || !(b_zeta_isomorph || qzb_rho != 0));
                    l = (l + 1).min(l_curr_max_set_size);
                    transposition_get_mcr_and_fix_set_and_unordered_partition(
                        &mut gamma,
                        &mut omega,
                        &mut phi,
                        n_tg,
                        l as usize,
                        &mut theta_from_gamma,
                    );
                    if unordered_partition_join(&theta_from_gamma, &mut theta, n_tg) == 0 {
                        lb = Lb::L11;
                        continue 'sm;
                    }
                    n_num_found_generators += 1;
                    if tvc == get_unordered_partition_mcr_node(&mut theta, tvc) {
                        lb = Lb::L11;
                        continue 'sm;
                    }
                    k = h_zeta;
                    lb = Lb::L13;
                }
                Lb::L11 => {
                    k = if lab { h_rho } else { h_zeta };
                    lb = Lb::L12;
                }
                Lb::L12 => {
                    if e[ku!()] == 1 && v[ku!()] != INCHI_CANON_INFINITY {
                        let wk = w[ku!()];
                        cell_intersect_with_set(&mut pi[ku!()], &wk, &omega, l as usize);
                    }
                    lb = Lb::L13;
                }
                Lb::L13 => {
                    if time_is_over() {
                        return CT_TIMEOUT_ERR;
                    }
                    if k == 0 {
                        break 'sm;
                    }
                    if lab && k < h_rho {
                        h_rho = k;
                    }
                    if k > h_zeta {
                        if v[ku!()] == INCHI_CANON_INFINITY {
                            k -= 1;
                            continue 'sm;
                        }
                        lb = Lb::L17;
                        continue 'sm;
                    }
                    if k == h_zeta {
                        lb = Lb::L14;
                        continue 'sm;
                    }
                    h_zeta = k;
                    tvh = cell_get_min_node(&pi[ku!()], &w[ku!()], 0, aux!());
                    tvc = tvh;
                    lb = Lb::L14;
                }
                Lb::L14 => {
                    if get_unordered_partition_mcr_node(&mut theta, v[ku!()])
                        == get_unordered_partition_mcr_node(&mut theta, tvh)
                    {
                        index += 1;
                    }
                    v[ku!()] = cell_get_min_node(&pi[ku!()], &w[ku!()], v[ku!()], aux!());
                    if v[ku!()] == INCHI_CANON_INFINITY {
                        lb = Lb::L16;
                        continue 'sm;
                    }
                    if v[ku!()] != get_unordered_partition_mcr_node(&mut theta, v[ku!()]) {
                        continue 'sm;
                    }
                    lb = Lb::L15;
                }
                Lb::L15 => {
                    t_lemma = t_lemma.min(k + 1);
                    hz_zeta = hz_zeta.min(k);
                    if lab {
                        if hz_rho >= k {
                            qzb_rho = 0;
                        }
                        if hz_rho > k {
                            hz_rho = k;
                        }
                        update_compare_layers(Some(&mut k_least_rho), hz_rho);
                    }
                    if pzb_rho_fix.is_some() {
                        if hzb_rho_fix >= k {
                            qzb_rho_fix = 0;
                        }
                        if hzb_rho_fix > k {
                            hzb_rho_fix = k;
                        }
                        update_compare_layers(Some(&mut k_least_rho_fix), hzb_rho_fix);
                    }
                    lb = Lb::L2;
                }
                Lb::L16 => {
                    if t_eq_zeta == k + 1
                        && index == cell_get_number_of_nodes(&pi[ku!()], &w[ku!()])
                    {
                        t_eq_zeta = k;
                    }
                    size *= f64::from(index);
                    index = 0;
                    k -= 1;
                    lb = Lb::L13;
                }
                Lb::L17 => {
                    if e[ku!()] == 0 && v[ku!()] != INCHI_CANON_INFINITY {
                        node_set_from_vertices(&mut cur_nodes, 1, &v, (k - 1) as usize);
                        for i in 1..=l {
                            if all_nodes_are_in_set(&cur_nodes, 1, &phi, i as usize) {
                                let wk = w[ku!()];
                                cell_intersect_with_set(&mut pi[ku!()], &wk, &omega, i as usize);
                            }
                        }
                    }
                    e[ku!()] = 1;
                    v[ku!()] = cell_get_min_node(&pi[ku!()], &w[ku!()], v[ku!()], aux!());
                    if v[ku!()] != INCHI_CANON_INFINITY {
                        lb = Lb::L15;
                        continue 'sm;
                    }
                    k -= 1;
                    lb = Lb::L13;
                }
            }
        }
    }

    // exit_function
    if !b_rho_is_discrete {
        return CT_CANON_ERR;
    }
    if let Some(z) = pzb_rho_fix {
        qzb_rho_fix = ct_full_compare(z, &pzb_rho, true, b_split_taut_compare);
        if qzb_rho_fix != 0 {
            return CT_CANON_ERR;
        }
    }
    n_symm_rank[..n_tg].iter_mut().for_each(|x| *x = 0);
    for i in 0..n_tg {
        let k = rho.at_number[i] as usize;
        let k2 = get_unordered_partition_mcr_node(&mut theta, (k + 1) as u16) as usize - 1;
        if n_symm_rank[k2] == 0 || n_symm_rank[k2] > rho.rank[k] {
            n_symm_rank[k2] = rho.rank[k];
        }
    }
    for i in 0..n_tg {
        let k = rho.at_number[i] as usize;
        let k2 = get_unordered_partition_mcr_node(&mut theta, (k + 1) as u16) as usize - 1;
        n_symm_rank[k] = n_symm_rank[k2];
    }
    n_canon_rank[..n_tg].copy_from_slice(&rho.rank[..n_tg]);
    n_atom_number_canon[..n_tg].copy_from_slice(&rho.at_number[..n_tg]);
    cd.n_len_linear_ct = pzb_rho.len_ct - 1;
    let len = cd.n_len_linear_ct.max(0) as usize;
    if !cd.linear_ct.is_empty() {
        cd.linear_ct[..len].copy_from_slice(&pzb_rho.ctbl[..len]);
    }
    cc.l_num_tot_ct = cc.l_num_decreased_ct + cc.l_num_rejected_ct + cc.l_num_equal_ct;
    cc.d_group_size = size;
    cc.l_num_generators = n_num_found_generators;
    cc.l_num_stored_isomorphisms = i64::from(l);
    if let Some(out) = pzb_rho_out
        && out.is_none()
    {
        *out = Some(pzb_rho);
    }
    ret
}

// ---------------------------------------------------------------------------
// helpers of the base canonical numbering (`GetBaseCanonRanking`)
// ---------------------------------------------------------------------------

use super::ichimap::{comp_ranks_ord, sorted_equ_info_to_ranks, sorted_ranks_to_equ_info};

/// `iOther`: second canonicalization of the no H structure
pub const CANON_FLAG_NO_H_RECANON: i32 = 0x0001;
/// `iOther`: NoTautH equivalence partition differs from NoH
pub const CANON_FLAG_NO_TAUT_H_DIFF: i32 = 0x0002;
/// `iOther`: equivalence partition in isotopic only non-taut differs from non-isotopic
pub const CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF: i32 = 0x0004;
/// `iBase`: isotopic equivalence partition in isotopic taut differs from non-isotopic taut
pub const CANON_FLAG_ISO_TAUT_DIFF: i32 = 0x0008;
/// `iOther`: isotopic equivalence partition in fixed H non-taut differs from non-isotopic fixed H
pub const CANON_FLAG_ISO_FIXED_H_DIFF: i32 = 0x0010;

/// `AT_ISO_SORT_KEY_MULT`: up to 32 identical hydrogen isotopes
pub const AT_ISO_SORT_KEY_MULT: i64 = 32;

/// `ATOM_SIZES`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AtomSizes {
    pub n_max_num_stereo_atoms: i32,
    pub n_max_num_stereo_bonds: i32,
    pub num_isotopic_atoms: i32,
    pub n_len_ct: i32,
    pub n_len_bonds: i32,
    pub n_len_isotopic: i32,
    pub n_len_ct_at_only: i32,
    pub n_len_linear_ct_stereo_dble: i32,
    pub n_len_linear_ct_stereo_carb: i32,
    pub b_may_have_stereo: i32,
    pub b_ignore_isotopic: i32,
    pub n_len_linear_ct_tautomer: i32,
    pub n_len_linear_ct_isotopic_tautomer: i32,
    pub b_has_isotopic_taut_groups: i32,
    pub n_len_isotopic_endpoints: i32,
}

/// `make_iso_sort_key`
pub fn make_iso_sort_key(iso_atw_diff: i32, num_1h: i32, num_2h: i32, num_3h: i32) -> i64 {
    let mut key = 0i64;
    let mut mult = 1i64;
    key += mult * i64::from(num_1h);
    mult *= AT_ISO_SORT_KEY_MULT;
    key += mult * i64::from(num_2h);
    mult *= AT_ISO_SORT_KEY_MULT;
    key += mult * i64::from(num_3h);
    mult *= AT_ISO_SORT_KEY_MULT;
    key += mult * i64::from(iso_atw_diff);
    key
}

/// `FixCanonEquivalenceInfo`: makes the current ranks agree with the
/// equivalence information; returns the number of different ranks and
/// `bChanged` (bit 0: symmetry ranks changed, bit 1: current ranks changed).
/// `n_temp_rank` is always the same as `n_curr_rank` on exit.
pub fn fix_canon_equivalence_info(
    num_at_tg: usize,
    n_symm_rank: &mut [u16],
    n_curr_rank: &mut [u16],
    n_temp_rank: &mut [u16],
    n_atom_number: &mut [u16],
) -> (i32, i32) {
    {
        let sr: &[u16] = n_symm_rank;
        inchi_qsort(&mut n_atom_number[..num_at_tg], |a, b| {
            comp_ranks_ord(sr, *a, *b)
        });
    }
    let (n_num_diff_ranks, b_change_symm_rank) =
        sorted_equ_info_to_ranks(n_symm_rank, n_temp_rank, n_atom_number, num_at_tg);
    let b_change_curr_rank = n_curr_rank[..num_at_tg] != n_temp_rank[..num_at_tg];
    if b_change_curr_rank {
        n_curr_rank[..num_at_tg].copy_from_slice(&n_temp_rank[..num_at_tg]);
    }
    if b_change_symm_rank {
        sorted_ranks_to_equ_info(n_symm_rank, n_temp_rank, n_atom_number, num_at_tg);
    }
    (
        n_num_diff_ranks,
        i32::from(b_change_symm_rank) | (2 * i32::from(b_change_curr_rank)),
    )
}

// ---------------------------------------------------------------------------
// GetBaseCanonRanking
// ---------------------------------------------------------------------------

use super::ichierr::{CT_ERR_MAX, CT_ERR_MIN, CT_ISOCOUNT_ERR, CT_UNKNOWN_ERR};
use super::ichimap::{create_neigh_list, differentiate_ranks2};
use super::ichitaut::{FLAG_NORM_CONSIDER_TAUT, T_NUM_NO_ISOTOPIC as T_NUM_NO_ISO};
use super::inpdef::AT_FLAG_ISO_H_POINT;
use super::runichi3::{TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE};

pub const TAUT_NON: usize = 0;
pub const TAUT_YES: usize = 1;

/// `FTCN`: the canonical numbering of a fixed-H or tautomeric structure
#[derive(Debug, Clone, Default)]
pub struct Ftcn {
    pub num_at_tg: usize,
    pub num_atoms: usize,
    pub n_canon_flags: i32,
    pub neigh_list: Option<NeighList>,
    pub linear_ct: Vec<u16>,
    pub n_len_linear_ct_at_only: i32,
    pub n_len_linear_ct: i32,
    pub n_max_len_linear_ct: i32,
    pub partition_ct: Partition,
    pub n_symm_rank_ct: Vec<u16>,
    pub n_num_h_orig: Vec<i16>,
    pub n_num_h: Vec<i16>,
    pub n_len_num_h: i32,
    pub n_num_h_orig_fix_h: Vec<i16>,
    pub n_num_h_fix_h: Vec<i16>,
    pub n_len_num_h_fix_h: i32,
    pub partition_ct_iso: Partition,
    pub n_symm_rank_ct_iso: Vec<u16>,
    pub iso_sort_keys: Vec<i64>,
    pub iso_sort_keys_orig: Vec<i64>,
    pub len_iso_sort_keys: i32,
    pub iso_exchg_atnos: Vec<i8>,
    pub iso_exchg_atnos_orig: Vec<i8>,
}

/// `BCN`: the base canonical numbering
#[derive(Debug, Clone, Default)]
pub struct Bcn {
    /// `pRankStack` as pairs (ranks, atom numbers)
    pub p_rank_stack: Vec<Partition>,
    pub n_max_len_rank_stack: i32,
    pub num_max: usize,
    pub num_at_tg: usize,
    pub num_atoms: usize,
    pub ftcn: [Ftcn; 2],
}

fn returned_error(v: i32) -> bool {
    (CT_ERR_MIN..=CT_ERR_MAX).contains(&v)
}

fn clean_iso_keys_if_any(keys: &mut [i64], len: usize) {
    if !keys.is_empty() {
        clean_iso_sort_keys(keys, len);
    }
}

fn clean_num_h_if_any(v: &mut [i16], len: usize) {
    if !v.is_empty() {
        clean_num_h(v, len);
    }
}

fn num_h_plus_base(num_h: i8) -> i16 {
    if num_h != 0 {
        i16::from(num_h) + BASE_H_NUMBER
    } else {
        EMPTY_H_NUMBER
    }
}

/// `GetBaseCanonRanking`: the base canonical numbering of a (possibly
/// tautomeric) structure with its fixed-H and isotopic variants. `at` holds
/// the tautomeric and non-tautomeric atoms. `time_is_over` replaces
/// `bInchiTimeIsOver`. Returns 0 or an error code; on success `bcn` has the
/// results.
#[allow(clippy::too_many_arguments, unused_assignments, unused_variables)]
pub fn get_base_canon_ranking(
    time_is_over: &mut dyn FnMut() -> bool,
    num_atoms: usize,
    mut num_at_tg: usize,
    at: [Option<&[SpAtom]>; 2],
    mut t_group_info: Option<&mut TGroupInfo>,
    s: &[AtomSizes; 2],
    bcn: &mut Bcn,
    b_fix_iso_fixed_h: bool,
    large_molecules: bool,
) -> i32 {
    let num_max = num_at_tg;
    let mut pi: Vec<Partition> = Vec::new();
    let mut pcd: [CanonData; 2] = [CanonData::default(), CanonData::default()];
    let mut cc = CanonCounts::default();
    let mut b_use_iso_aux_base = [false; 2];
    let mut b_taut_ignore_isotopic = 0i32;

    // select the base structure
    let (i_base, i_other, b_req_non_taut, b_req_taut);
    let at_yes = at[TAUT_YES].is_some() && s[TAUT_YES].n_len_ct != 0;
    let at_non = at[TAUT_NON].is_some() && s[TAUT_NON].n_len_ct != 0;
    let taut_cond = at_yes
        && t_group_info.as_deref().is_some_and(|t| {
            s[TAUT_YES].n_len_linear_ct_tautomer > 0
                && ((t.t_group.is_some() && t.num_t_groups > 0)
                    || t.tni.b_normalization_flags & FLAG_NORM_CONSIDER_TAUT != 0
                    || (t.n_num_isotopic_endpoints > 1
                        && t.b_taut_flags_done
                            & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                            != 0))
        });
    if taut_cond {
        i_base = TAUT_YES;
        b_req_taut = true;
        let t = t_group_info.as_deref().unwrap();
        b_use_iso_aux_base[i_base] = s[i_base].n_len_isotopic_endpoints > 1
            && t.b_taut_flags_done
                & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                != 0;
        if at_non {
            i_other = TAUT_NON;
            b_req_non_taut = true;
        } else {
            i_other = i_base;
            b_req_non_taut = false;
        }
    } else if at_non {
        i_base = TAUT_NON;
        b_req_taut = false;
        i_other = i_base;
        b_req_non_taut = true;
        num_at_tg = num_atoms;
    } else if at_yes {
        i_base = TAUT_YES;
        b_req_taut = false;
        b_use_iso_aux_base[i_base] = s[i_base].n_len_isotopic_endpoints > 1;
        i_other = i_base;
        b_req_non_taut = true;
        num_at_tg = num_atoms;
    } else {
        return CT_UNKNOWN_ERR;
    }
    if b_req_taut {
        let t = t_group_info.as_deref_mut().unwrap();
        b_taut_ignore_isotopic = t.b_ignore_isotopic;
        t.b_ignore_isotopic = 1;
    }
    let at_base = at[i_base].unwrap();
    let at_other = at[i_other];

    let ret: i32 = 'run: {
        macro_rules! fail {
            ($r:expr) => {{
                let r = $r;
                break 'run if returned_error(r) { r } else { CT_CANON_ERR };
            }};
        }
        let mut p_atom_invariant = vec![AtomInvariant2::default(); num_max];
        let mut n_symm_rank_no_h = vec![0u16; num_max];
        let mut n_canon_rank_no_h = vec![0u16; num_max];
        let mut n_atom_number_canon_no_h = vec![0u16; num_max];
        let mut n_rank = vec![0u16; num_max];
        let mut n_atom_number = vec![0u16; num_max];
        let mut n_temp_rank = vec![0u16; num_max];
        let mut n_rank_aux = vec![0u16; num_max];
        let mut n_atom_number_aux = vec![0u16; num_max];
        let mut p_atom_invariant_aux = vec![AtomInvariant2::default(); num_max];

        let mut neigh_yes: Option<NeighList> = None;
        let mut neigh_non_opt: Option<NeighList>;
        if b_req_taut {
            neigh_yes = Some(create_neigh_list(
                num_atoms,
                num_at_tg,
                at_base,
                false,
                t_group_info.as_deref(),
            ));
            neigh_non_opt = Some(create_neigh_list(
                num_atoms, num_atoms, at_base, false, None,
            ));
        } else {
            neigh_non_opt = Some(create_neigh_list(
                num_atoms, num_atoms, at_base, false, None,
            ));
        }
        bcn.n_max_len_rank_stack = 0;
        bcn.num_max = num_max;
        bcn.num_at_tg = num_at_tg;
        bcn.num_atoms = num_atoms;

        let mut l_count = 0i64;
        fill_out_atom_invariant2(
            at_base,
            num_atoms,
            num_atoms,
            &mut p_atom_invariant,
            true,
            false,
            false,
            false,
            false,
            None,
        );
        let mut n_num_curr_ranks = set_initial_ranks2(
            num_atoms,
            &p_atom_invariant,
            &mut n_rank,
            &mut n_atom_number,
        );
        n_num_curr_ranks = differentiate_ranks2(
            num_atoms,
            neigh_non_opt.as_mut().unwrap(),
            n_num_curr_ranks,
            &mut n_rank,
            &mut n_temp_rank,
            &mut n_atom_number,
            &mut l_count,
            0,
        );
        let n_max_len_rank_stack = 2 * (num_at_tg as i32 - n_num_curr_ranks) + 8;
        pi = vec![Partition::default(); (n_max_len_rank_stack / 2).max(1) as usize];
        bcn.n_max_len_rank_stack = n_max_len_rank_stack;
        pi[0] = Partition {
            rank: std::mem::take(&mut n_rank),
            at_number: std::mem::take(&mut n_atom_number),
        };

        let io = i_other;
        let ib = i_base;
        // the canonical search with the current `pcd[idx]`
        macro_rules! canon {
            ($which:expr, $ntg:expr, $dig:expr, $sym:expr, $can:expr, $atn:expr, $idx:expr, $inp:expr, $out:expr) => {{
                let g: &mut NeighList = if $which == TAUT_NON {
                    neigh_non_opt.as_mut().unwrap()
                } else {
                    neigh_yes.as_mut().unwrap()
                };
                canon_graph(
                    &mut *time_is_over,
                    num_atoms,
                    $ntg,
                    num_max,
                    $dig,
                    g,
                    &mut pi,
                    &mut $sym,
                    &mut $can,
                    &mut $atn,
                    &mut pcd[$idx],
                    &mut cc,
                    $inp,
                    Some(&mut $out),
                    large_molecules,
                )
            }};
        }
        macro_rules! new_ct {
            ($idx:expr) => {{
                let cd = &pcd[$idx];
                c_table_create(
                    num_atoms,
                    cd.n_max_len_linear_ct,
                    cd.n_len_ct_at_only,
                    cd.num_h.as_ref().map(|_| cd.maxlen_num_h),
                    cd.num_h_fixed.as_ref().map(|_| cd.maxlen_num_h_fixed),
                    cd.iso_sort_key.as_ref().map(|_| cd.maxlen_iso_sort_key),
                    cd.iso_exchg_atnos
                        .as_ref()
                        .map(|_| cd.maxlen_iso_exchg_atnos),
                )
            }};
        }
        // the stage set-up shared by the atoms-only stages of `iOther`
        macro_rules! other_stage {
            () => {{
                let cd = &mut pcd[io];
                cd.linear_ct = Vec::new();
                cd.n_max_len_linear_ct = s[io].n_len_ct_at_only + 1;
                cd.n_len_linear_ct = s[io].n_len_ct_at_only;
                cd.n_len_ct_at_only = s[io].n_len_ct_at_only;
                cd.n_aux_rank = None;
            }};
        }

        // ---- NoH ----
        pcd[io].num_h = None;
        pcd[io].num_h_fixed = None;
        pcd[io].iso_sort_key = None;
        pcd[io].iso_exchg_atnos = None;
        pcd[io].maxlen_iso_sort_key = 0;
        pcd[io].maxlen_iso_exchg_atnos = 0;
        other_stage!();
        pcd[io].n_aux_rank = None;
        let mut ct_no_h: Option<ConTable> = None;
        let mut ct_no_taut_h: Option<ConTable> = None;
        let mut sym_no_h = std::mem::take(&mut n_symm_rank_no_h);
        let mut can_no_h = std::mem::take(&mut n_canon_rank_no_h);
        let mut atn_no_h = std::mem::take(&mut n_atom_number_canon_no_h);
        let ret = canon!(
            TAUT_NON, num_atoms, false, sym_no_h, can_no_h, atn_no_h, io, None, ct_no_h
        );
        if ret < 0 {
            fail!(ret);
        }
        let mut b_changed = 0;
        {
            let p0 = &mut pi[0];
            let (nn, bc) = fix_canon_equivalence_info(
                num_atoms,
                &mut sym_no_h,
                &mut p0.rank,
                &mut n_temp_rank,
                &mut p0.at_number,
            );
            n_num_curr_ranks = nn;
            b_changed = bc;
        }
        if b_changed & 3 != 0 {
            ct_no_h = None;
            pcd[io].n_canon_flags |= CANON_FLAG_NO_H_RECANON;
            let ret = canon!(
                TAUT_NON, num_atoms, false, sym_no_h, can_no_h, atn_no_h, io, None, ct_no_h
            );
            if ret < 0 {
                fail!(ret);
            }
        }

        // ---- NoTautH ----
        let maxlen_num_h_no_taut_h = num_atoms as i32 + 1;
        let mut sym_no_taut_h = vec![0u16; num_max];
        let mut can_no_taut_h = vec![0u16; num_max];
        let mut atn_no_taut_h = vec![0u16; num_max];
        let mut num_h_no_taut_h = vec![0i16; maxlen_num_h_no_taut_h as usize];
        for i in 0..num_atoms {
            num_h_no_taut_h[i] = if at_base[i].endpoint == 0 && at_base[i].num_h != 0 {
                i16::from(at_base[i].num_h) + BASE_H_NUMBER
            } else {
                EMPTY_H_NUMBER
            };
        }
        other_stage!();
        pcd[io].num_h = Some(num_h_no_taut_h.clone());
        pcd[io].num_h_fixed = None;
        pcd[io].iso_sort_key = None;
        pcd[io].iso_exchg_atnos = None;
        pcd[io].maxlen_num_h = maxlen_num_h_no_taut_h;
        pcd[io].maxlen_iso_sort_key = 0;
        pcd[io].maxlen_iso_exchg_atnos = 0;

        n_temp_rank.iter_mut().for_each(|x| *x = 0);
        for i in 0..num_atoms {
            let c = sym_no_h[i] as usize - 1;
            if (n_temp_rank[c] as usize) < i {
                n_temp_rank[c] = i as u16;
            }
        }
        let mut i = 0;
        while i < num_atoms {
            if num_h_no_taut_h[i] != num_h_no_taut_h[n_temp_rank[sym_no_h[i] as usize - 1] as usize]
            {
                pcd[io].n_canon_flags |= CANON_FLAG_NO_TAUT_H_DIFF;
                break;
            }
            i += 1;
        }
        if i < num_atoms {
            for x in p_atom_invariant_aux.iter_mut() {
                *x = AtomInvariant2::default();
            }
            for i in 0..num_atoms {
                p_atom_invariant_aux[i].val[0] = sym_no_h[i];
                p_atom_invariant_aux[i].val[1] = num_h_no_taut_h[i] as u16;
            }
            n_num_curr_ranks = set_initial_ranks2(
                num_atoms,
                &p_atom_invariant_aux,
                &mut n_rank_aux,
                &mut n_atom_number_aux,
            );
            n_num_curr_ranks = differentiate_ranks2(
                num_atoms,
                neigh_non_opt.as_mut().unwrap(),
                n_num_curr_ranks,
                &mut n_rank_aux,
                &mut n_temp_rank,
                &mut n_atom_number_aux,
                &mut l_count,
                0,
            );
            pcd[io].n_aux_rank = Some(n_rank_aux.clone());
            let ret = canon!(
                TAUT_NON,
                num_atoms,
                true,
                sym_no_taut_h,
                can_no_taut_h,
                atn_no_taut_h,
                io,
                ct_no_h.as_ref(),
                ct_no_taut_h
            );
            if ret < 0 {
                fail!(ret);
            }
        } else {
            let mut tmp = new_ct!(io);
            ct_full_copy(&mut tmp, ct_no_h.as_ref().unwrap());
            for i in 0..num_atoms {
                tmp.num_h[can_no_h[i] as usize - 1] = num_h_no_taut_h[i];
            }
            tmp.len_num_h = num_atoms as i32;
            ct_no_taut_h = Some(tmp);
            sym_no_taut_h[..num_atoms].copy_from_slice(&sym_no_h[..num_atoms]);
            can_no_taut_h[..num_atoms].copy_from_slice(&can_no_h[..num_atoms]);
            atn_no_taut_h[..num_atoms].copy_from_slice(&atn_no_h[..num_atoms]);
        }

        // ---- isotopic non-tautomeric (only when tautomerism was not found) ----
        let mut ct_no_taut_h_iso: Option<ConTable> = None;
        let mut sym_no_taut_h_iso = Vec::new();
        let mut can_no_taut_h_iso = Vec::new();
        let mut atn_no_taut_h_iso = Vec::new();
        let mut iso_sort_key_no_taut_h: Vec<i64> = Vec::new();
        let mut len_iso_sort_key_no_taut_h = 0i32;
        let mut maxlen_iso_sort_key_no_taut_h = 0i32;
        let mut num_iso_no_taut_h = 0i32;
        if s[io].num_isotopic_atoms != 0
            && s[io].b_ignore_isotopic == 0
            && !b_req_taut
            && b_req_non_taut
        {
            maxlen_iso_sort_key_no_taut_h = num_atoms as i32 + 1;
            sym_no_taut_h_iso = vec![0u16; num_max];
            can_no_taut_h_iso = vec![0u16; num_max];
            atn_no_taut_h_iso = vec![0u16; num_max];
            iso_sort_key_no_taut_h = vec![0i64; maxlen_iso_sort_key_no_taut_h as usize];
            num_iso_no_taut_h = 0;
            for i in 0..num_atoms {
                let a = &at_base[i];
                let key = if a.endpoint != 0 {
                    make_iso_sort_key(i32::from(a.iso_atw_diff), 0, 0, 0)
                } else {
                    make_iso_sort_key(
                        i32::from(a.iso_atw_diff),
                        i32::from(a.num_iso_h[0]),
                        i32::from(a.num_iso_h[1]),
                        i32::from(a.num_iso_h[2]),
                    )
                };
                if key != 0 {
                    iso_sort_key_no_taut_h[i] = key;
                    num_iso_no_taut_h += 1;
                } else {
                    iso_sort_key_no_taut_h[i] = EMPTY_ISO_SORT_KEY;
                }
            }
            other_stage!();
            pcd[io].num_h = Some(num_h_no_taut_h.clone());
            pcd[io].num_h_fixed = None;
            pcd[io].iso_sort_key = Some(iso_sort_key_no_taut_h.clone());
            pcd[io].iso_exchg_atnos = None;
            pcd[io].maxlen_num_h = maxlen_num_h_no_taut_h;
            len_iso_sort_key_no_taut_h = num_atoms as i32;
            pcd[io].maxlen_iso_sort_key = maxlen_iso_sort_key_no_taut_h;
            pcd[io].maxlen_iso_exchg_atnos = 0;
            let mut i = 0;
            if num_iso_no_taut_h != 0 {
                n_temp_rank.iter_mut().for_each(|x| *x = 0);
                for i in 0..num_atoms {
                    let c = sym_no_taut_h[i] as usize - 1;
                    if (n_temp_rank[c] as usize) < i {
                        n_temp_rank[c] = i as u16;
                    }
                }
                while i < num_atoms {
                    if iso_sort_key_no_taut_h[i]
                        != iso_sort_key_no_taut_h
                            [n_temp_rank[sym_no_taut_h[i] as usize - 1] as usize]
                    {
                        pcd[io].n_canon_flags |= CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF;
                        break;
                    }
                    i += 1;
                }
            } else {
                i = num_atoms;
            }
            if i < num_atoms {
                for x in p_atom_invariant_aux.iter_mut() {
                    *x = AtomInvariant2::default();
                }
                for i in 0..num_atoms {
                    p_atom_invariant_aux[i].val[0] = sym_no_taut_h[i];
                    p_atom_invariant_aux[i].iso_sort_key = iso_sort_key_no_taut_h[i];
                }
                n_num_curr_ranks = set_initial_ranks2(
                    num_atoms,
                    &p_atom_invariant_aux,
                    &mut n_rank_aux,
                    &mut n_atom_number_aux,
                );
                n_num_curr_ranks = differentiate_ranks2(
                    num_atoms,
                    neigh_non_opt.as_mut().unwrap(),
                    n_num_curr_ranks,
                    &mut n_rank_aux,
                    &mut n_temp_rank,
                    &mut n_atom_number_aux,
                    &mut l_count,
                    0,
                );
                pcd[io].n_aux_rank = Some(n_rank_aux.clone());
                let ret = canon!(
                    TAUT_NON,
                    num_atoms,
                    true,
                    sym_no_taut_h_iso,
                    can_no_taut_h_iso,
                    atn_no_taut_h_iso,
                    io,
                    ct_no_taut_h.as_ref(),
                    ct_no_taut_h_iso
                );
                if ret < 0 {
                    fail!(ret);
                }
            } else {
                let mut tmp = new_ct!(io);
                ct_full_copy(&mut tmp, ct_no_taut_h.as_ref().unwrap());
                for i in 0..num_atoms {
                    tmp.iso_sort_key[can_no_taut_h[i] as usize - 1] = iso_sort_key_no_taut_h[i];
                }
                tmp.len_iso_sort_key = num_atoms as i32;
                ct_no_taut_h_iso = Some(tmp);
                sym_no_taut_h_iso[..num_atoms].copy_from_slice(&sym_no_taut_h[..num_atoms]);
                can_no_taut_h_iso[..num_atoms].copy_from_slice(&can_no_taut_h[..num_atoms]);
                atn_no_taut_h_iso[..num_atoms].copy_from_slice(&atn_no_taut_h[..num_atoms]);
            }
        }

        // ---- tautomeric base ----
        let mut ct_base: Option<ConTable> = None;
        let mut ct_base_iso: Option<ConTable> = None;
        let mut sym_base = Vec::new();
        let mut can_base = Vec::new();
        let mut atn_base = Vec::new();
        let mut sym_base_iso = Vec::new();
        let mut can_base_iso = Vec::new();
        let mut atn_base_iso = Vec::new();
        let mut num_h: Vec<i16> = Vec::new();
        let mut maxlen_num_h = 0i32;
        let mut iso_sort_key_base: Vec<i64> = Vec::new();
        let mut len_iso_sort_key_base = 0i32;
        let mut iso_exchg_atnos: Vec<i8> = Vec::new();
        let mut num_iso_base = 0i32;
        if b_req_taut {
            maxlen_num_h = (num_atoms + T_NUM_NO_ISO * (num_at_tg - num_atoms) + 1) as i32;
            sym_base = vec![0u16; num_max];
            can_base = vec![0u16; num_max];
            atn_base = vec![0u16; num_max];
            num_h = vec![0i16; maxlen_num_h as usize];
            let tgi = t_group_info.as_deref().unwrap();
            for i in 0..num_atoms {
                num_h[i] = if at_base[i].endpoint == 0 && at_base[i].num_h != 0 {
                    i16::from(at_base[i].num_h) + BASE_H_NUMBER
                } else {
                    EMPTY_H_NUMBER
                };
            }
            let mut k = num_atoms;
            for i in num_atoms..num_at_tg {
                let m = i - num_atoms;
                for j in 0..T_NUM_NO_ISO {
                    let v = tgi.t_group.as_ref().unwrap()[m].num[j];
                    num_h[k] = if v != 0 {
                        (i32::from(v) + i32::from(BASE_H_NUMBER)) as i16
                    } else {
                        EMPTY_H_NUMBER
                    };
                    k += 1;
                }
            }
            let len_num_h = k;
            let cd = &mut pcd[ib];
            cd.linear_ct = Vec::new();
            cd.num_h = Some(num_h.clone());
            cd.num_h_fixed = None;
            cd.iso_sort_key = None;
            cd.iso_exchg_atnos = None;
            cd.n_max_len_linear_ct = s[ib].n_len_ct + 1;
            cd.maxlen_num_h = maxlen_num_h;
            cd.n_len_linear_ct = s[ib].n_len_ct;
            cd.n_len_ct_at_only = s[ib].n_len_ct_at_only;
            cd.maxlen_iso_sort_key = 0;
            cd.maxlen_iso_exchg_atnos = 0;
            cd.n_aux_rank = None;

            fill_out_atom_invariant2(
                at_base,
                num_atoms,
                num_at_tg,
                &mut p_atom_invariant,
                true,
                false,
                false,
                true,
                true,
                t_group_info.as_deref(),
            );
            for i in 0..num_atoms {
                p_atom_invariant[i].val[0] = pi[0].rank[i];
            }
            {
                let p0 = &mut pi[0];
                n_num_curr_ranks = set_initial_ranks2(
                    num_at_tg,
                    &p_atom_invariant,
                    &mut p0.rank,
                    &mut p0.at_number,
                );
                n_num_curr_ranks = differentiate_ranks4(
                    num_at_tg,
                    neigh_yes.as_mut().unwrap(),
                    n_num_curr_ranks,
                    &mut p0.rank,
                    &mut n_temp_rank,
                    &mut p0.at_number,
                    num_atoms as u16,
                    &mut l_count,
                );
            }
            for x in p_atom_invariant_aux.iter_mut() {
                *x = AtomInvariant2::default();
            }
            for i in 0..num_atoms {
                p_atom_invariant_aux[i].val[0] = sym_no_taut_h[i];
                p_atom_invariant_aux[i].val[1] = num_h[i] as u16;
            }
            for j in num_atoms..num_at_tg {
                p_atom_invariant_aux[j].val[0] = pi[0].rank[j];
            }
            n_num_curr_ranks = set_initial_ranks2(
                num_at_tg,
                &p_atom_invariant_aux,
                &mut n_rank_aux,
                &mut n_atom_number_aux,
            );
            n_num_curr_ranks = differentiate_ranks4(
                num_at_tg,
                neigh_yes.as_mut().unwrap(),
                n_num_curr_ranks,
                &mut n_rank_aux,
                &mut n_temp_rank,
                &mut n_atom_number_aux,
                num_atoms as u16,
                &mut l_count,
            );
            pcd[ib].n_aux_rank = Some(n_rank_aux.clone());
            let ret = canon!(
                TAUT_YES,
                num_at_tg,
                true,
                sym_base,
                can_base,
                atn_base,
                ib,
                ct_no_taut_h.as_ref(),
                ct_base
            );
            if ret < 0 {
                fail!(ret);
            }

            // tautomeric isotopic structure
            let tgi_now = t_group_info.as_deref().unwrap();
            if (s[ib].num_isotopic_atoms != 0 && s[ib].b_ignore_isotopic == 0)
                || (s[ib].b_has_isotopic_taut_groups != 0 && b_taut_ignore_isotopic == 0)
                || (b_use_iso_aux_base[ib] && b_taut_ignore_isotopic == 0)
            {
                t_group_info.as_deref_mut().unwrap().b_ignore_isotopic = b_taut_ignore_isotopic;
                sym_base_iso = vec![0u16; num_max];
                can_base_iso = vec![0u16; num_max];
                atn_base_iso = vec![0u16; num_max];
                let mut maxlen_iso_exchg_atnos = 0i32;
                let mut len_iso_exchg_atnos = 0i32;
                if b_use_iso_aux_base[ib] {
                    maxlen_iso_exchg_atnos = num_max as i32 + 1;
                    iso_exchg_atnos = vec![0i8; maxlen_iso_exchg_atnos as usize];
                }
                let mut maxlen_iso_sort_key_base = num_max as i32 + 1;
                iso_sort_key_base = vec![0i64; maxlen_iso_sort_key_base as usize];
                num_iso_no_taut_h = 0;
                let mut num_iso_no_aux_base = 0i32;
                if !iso_exchg_atnos.is_empty() {
                    len_iso_exchg_atnos = num_at_tg as i32;
                }
                let tgi = t_group_info.as_deref().unwrap();
                for i in 0..num_atoms {
                    let a = &at_base[i];
                    let key;
                    if a.endpoint != 0
                        || (!iso_exchg_atnos.is_empty() && a.c_flags & AT_FLAG_ISO_H_POINT != 0)
                    {
                        key = make_iso_sort_key(i32::from(a.iso_atw_diff), 0, 0, 0);
                        if !iso_exchg_atnos.is_empty() {
                            num_iso_no_aux_base += i32::from(a.endpoint == 0);
                        }
                    } else {
                        key = make_iso_sort_key(
                            i32::from(a.iso_atw_diff),
                            i32::from(a.num_iso_h[0]),
                            i32::from(a.num_iso_h[1]),
                            i32::from(a.num_iso_h[2]),
                        );
                        if !iso_exchg_atnos.is_empty() {
                            iso_exchg_atnos[i] = 1;
                        }
                    }
                    if key != 0 {
                        num_iso_no_taut_h += 1;
                        iso_sort_key_base[i] = key;
                    } else {
                        iso_sort_key_base[i] = EMPTY_ISO_SORT_KEY;
                    }
                }
                if !iso_exchg_atnos.is_empty() {
                    if num_iso_no_aux_base != i32::from(tgi.n_isotopic_endpoint_atom_number[0]) {
                        fail!(CT_ISOCOUNT_ERR);
                    }
                    for i in 1..=num_iso_no_aux_base as usize {
                        let j = tgi.n_isotopic_endpoint_atom_number[i] as usize;
                        if at_base[j].endpoint != 0 || at_base[j].c_flags & AT_FLAG_ISO_H_POINT == 0
                        {
                            fail!(CT_ISOCOUNT_ERR);
                        }
                    }
                }
                num_iso_base = 0;
                if !iso_exchg_atnos.is_empty() {
                    for i in num_atoms..num_at_tg {
                        iso_sort_key_base[i] = EMPTY_ISO_SORT_KEY;
                    }
                } else {
                    for i in num_atoms..num_at_tg {
                        let m = i - num_atoms;
                        let key = tgi.t_group.as_ref().unwrap()[m].i_weight;
                        if key != 0 {
                            num_iso_base += 1;
                            iso_sort_key_base[i] = key;
                        } else {
                            iso_sort_key_base[i] = EMPTY_ISO_SORT_KEY;
                        }
                    }
                }
                if num_iso_no_aux_base == 0 && !iso_exchg_atnos.is_empty() {
                    iso_exchg_atnos = Vec::new();
                    len_iso_exchg_atnos = 0;
                    maxlen_iso_exchg_atnos = 0;
                }
                if num_iso_no_taut_h == 0 && num_iso_base == 0 && !iso_sort_key_base.is_empty() {
                    iso_sort_key_base = Vec::new();
                    maxlen_iso_sort_key_base = 0;
                } else {
                    len_iso_sort_key_base = num_at_tg as i32;
                }
                if iso_exchg_atnos.is_empty() && iso_sort_key_base.is_empty() {
                    sym_base_iso = Vec::new();
                    can_base_iso = Vec::new();
                    atn_base_iso = Vec::new();
                } else {
                    let cd = &mut pcd[ib];
                    cd.linear_ct = Vec::new();
                    cd.num_h = Some(num_h.clone());
                    cd.num_h_fixed = None;
                    cd.iso_sort_key =
                        (!iso_sort_key_base.is_empty()).then(|| iso_sort_key_base.clone());
                    cd.iso_exchg_atnos =
                        (!iso_exchg_atnos.is_empty()).then(|| iso_exchg_atnos.clone());
                    cd.n_max_len_linear_ct = s[ib].n_len_ct + 1;
                    cd.maxlen_num_h = maxlen_num_h;
                    cd.n_len_linear_ct = s[ib].n_len_ct;
                    cd.n_len_ct_at_only = s[ib].n_len_ct_at_only;
                    cd.maxlen_iso_sort_key = maxlen_iso_sort_key_base;
                    cd.maxlen_iso_exchg_atnos = maxlen_iso_exchg_atnos;
                    cd.n_aux_rank = None;
                    let mut i;
                    if num_iso_no_taut_h != 0 || num_iso_base != 0 || num_iso_no_aux_base != 0 {
                        n_temp_rank.iter_mut().for_each(|x| *x = 0);
                        for i in 0..num_at_tg {
                            let c = sym_base[i] as usize - 1;
                            if (n_temp_rank[c] as usize) < i {
                                n_temp_rank[c] = i as u16;
                            }
                        }
                        i = 0;
                        while i < num_at_tg {
                            let rep = n_temp_rank[sym_base[i] as usize - 1] as usize;
                            if (!iso_sort_key_base.is_empty()
                                && iso_sort_key_base[i] != iso_sort_key_base[rep])
                                || (!iso_exchg_atnos.is_empty()
                                    && iso_exchg_atnos[i] != iso_exchg_atnos[rep])
                            {
                                pcd[ib].n_canon_flags |= CANON_FLAG_ISO_TAUT_DIFF;
                                break;
                            }
                            i += 1;
                        }
                    } else {
                        i = num_at_tg;
                    }
                    if i < num_at_tg {
                        for x in p_atom_invariant_aux.iter_mut() {
                            *x = AtomInvariant2::default();
                        }
                        for i in 0..num_at_tg {
                            p_atom_invariant_aux[i].val[0] = sym_base[i];
                            p_atom_invariant_aux[i].iso_sort_key = if iso_sort_key_base.is_empty() {
                                0
                            } else {
                                iso_sort_key_base[i]
                            };
                            p_atom_invariant_aux[i].iso_aux_key = if iso_exchg_atnos.is_empty() {
                                0
                            } else {
                                iso_exchg_atnos[i]
                            };
                        }
                        n_num_curr_ranks = set_initial_ranks2(
                            num_at_tg,
                            &p_atom_invariant_aux,
                            &mut n_rank_aux,
                            &mut n_atom_number_aux,
                        );
                        n_num_curr_ranks = differentiate_ranks2(
                            num_at_tg,
                            neigh_yes.as_mut().unwrap(),
                            n_num_curr_ranks,
                            &mut n_rank_aux,
                            &mut n_temp_rank,
                            &mut n_atom_number_aux,
                            &mut l_count,
                            0,
                        );
                        pcd[ib].n_aux_rank = Some(n_rank_aux.clone());
                        let ret = canon!(
                            TAUT_YES,
                            num_at_tg,
                            true,
                            sym_base_iso,
                            can_base_iso,
                            atn_base_iso,
                            ib,
                            ct_base.as_ref(),
                            ct_base_iso
                        );
                        if ret < 0 {
                            fail!(ret);
                        }
                    } else {
                        let mut tmp = new_ct!(ib);
                        ct_full_copy(&mut tmp, ct_base.as_ref().unwrap());
                        if !iso_sort_key_base.is_empty() {
                            for i in 0..num_at_tg {
                                tmp.iso_sort_key[can_base[i] as usize - 1] = iso_sort_key_base[i];
                            }
                            tmp.len_iso_sort_key = num_at_tg as i32;
                        } else {
                            tmp.len_iso_sort_key = 0;
                        }
                        if !iso_exchg_atnos.is_empty() {
                            for i in 0..num_atoms {
                                tmp.iso_exchg_atnos[can_base[i] as usize - 1] = iso_exchg_atnos[i];
                            }
                            tmp.len_iso_exchg_atnos = num_at_tg as i32;
                        } else {
                            tmp.len_iso_exchg_atnos = 0;
                        }
                        ct_base_iso = Some(tmp);
                        sym_base_iso[..num_at_tg].copy_from_slice(&sym_base[..num_at_tg]);
                        can_base_iso[..num_at_tg].copy_from_slice(&can_base[..num_at_tg]);
                        atn_base_iso[..num_at_tg].copy_from_slice(&atn_base[..num_at_tg]);
                    }
                    t_group_info.as_deref_mut().unwrap().b_ignore_isotopic = 1;
                }
            }
        }

        // ---- "fixed H" ----
        let mut ct_fix_h: Option<ConTable> = None;
        let mut ct_fix_h_iso: Option<ConTable> = None;
        let mut sym_fix_h = Vec::new();
        let mut can_fix_h = Vec::new();
        let mut atn_fix_h = Vec::new();
        let mut sym_fix_h_iso = Vec::new();
        let mut can_fix_h_iso = Vec::new();
        let mut atn_fix_h_iso = Vec::new();
        let mut num_h_fixed: Vec<i16> = Vec::new();
        if b_req_taut && b_req_non_taut {
            let maxlen_num_h_fixed = num_atoms as i32 + 1;
            sym_fix_h = vec![0u16; num_max];
            can_fix_h = vec![0u16; num_max];
            atn_fix_h = vec![0u16; num_max];
            num_h_fixed = vec![0i16; maxlen_num_h_fixed as usize];
            let at_other = at_other.unwrap();
            for i in 0..num_atoms {
                num_h_fixed[i] = if at_base[i].endpoint != 0 {
                    num_h_plus_base(at_other[i].num_h)
                } else if at_other[i].num_h != at_base[i].num_h {
                    i16::from(at_other[i].num_h) - i16::from(at_base[i].num_h) + BASE_H_NUMBER
                } else {
                    EMPTY_H_NUMBER
                };
            }
            other_stage!();
            pcd[io].num_h = Some(num_h_no_taut_h.clone());
            pcd[io].num_h_fixed = Some(num_h_fixed.clone());
            pcd[io].iso_sort_key = None;
            pcd[io].iso_exchg_atnos = None;
            pcd[io].maxlen_num_h = maxlen_num_h_no_taut_h;
            pcd[io].maxlen_num_h_fixed = maxlen_num_h_fixed;
            pcd[io].maxlen_iso_sort_key = 0;
            pcd[io].maxlen_iso_exchg_atnos = 0;
            for x in p_atom_invariant_aux.iter_mut() {
                *x = AtomInvariant2::default();
            }
            for i in 0..num_atoms {
                p_atom_invariant_aux[i].val[0] = sym_base[i];
                p_atom_invariant_aux[i].val[1] = num_h_fixed[i] as u16;
            }
            n_num_curr_ranks = set_initial_ranks2(
                num_atoms,
                &p_atom_invariant_aux,
                &mut n_rank_aux,
                &mut n_atom_number_aux,
            );
            n_num_curr_ranks = differentiate_ranks2(
                num_atoms,
                neigh_non_opt.as_mut().unwrap(),
                n_num_curr_ranks,
                &mut n_rank_aux,
                &mut n_temp_rank,
                &mut n_atom_number_aux,
                &mut l_count,
                0,
            );
            pcd[io].n_aux_rank = Some(n_rank_aux.clone());
            let ret = canon!(
                TAUT_NON,
                num_atoms,
                false,
                sym_fix_h,
                can_fix_h,
                atn_fix_h,
                io,
                ct_no_taut_h.as_ref(),
                ct_fix_h
            );
            if ret < 0 {
                fail!(ret);
            }

            // "fixed H" isotopic
            let mut iflag = (s[ib].num_isotopic_atoms != 0 && s[ib].b_ignore_isotopic == 0)
                || (s[ib].b_has_isotopic_taut_groups != 0 && b_taut_ignore_isotopic == 0);
            if b_fix_iso_fixed_h {
                iflag = iflag || (s[io].num_isotopic_atoms != 0 && s[io].b_ignore_isotopic == 0);
            }
            if iflag {
                maxlen_iso_sort_key_no_taut_h = num_atoms as i32 + 1;
                sym_fix_h_iso = vec![0u16; num_max];
                can_fix_h_iso = vec![0u16; num_max];
                atn_fix_h_iso = vec![0u16; num_max];
                iso_sort_key_no_taut_h = vec![0i64; maxlen_iso_sort_key_no_taut_h as usize];
                for i in 0..num_atoms {
                    let a = if b_fix_iso_fixed_h {
                        &at_other[i]
                    } else {
                        &at_base[i]
                    };
                    let key = make_iso_sort_key(
                        i32::from(a.iso_atw_diff),
                        i32::from(a.num_iso_h[0]),
                        i32::from(a.num_iso_h[1]),
                        i32::from(a.num_iso_h[2]),
                    );
                    if key != 0 {
                        iso_sort_key_no_taut_h[i] = key;
                        num_iso_no_taut_h += 1;
                    } else {
                        iso_sort_key_no_taut_h[i] = EMPTY_ISO_SORT_KEY;
                    }
                }
                other_stage!();
                pcd[io].num_h = Some(num_h_no_taut_h.clone());
                pcd[io].num_h_fixed = Some(num_h_fixed.clone());
                pcd[io].iso_sort_key = Some(iso_sort_key_no_taut_h.clone());
                pcd[io].iso_exchg_atnos = None;
                pcd[io].maxlen_num_h = maxlen_num_h_no_taut_h;
                pcd[io].maxlen_num_h_fixed = maxlen_num_h_fixed;
                len_iso_sort_key_no_taut_h = num_atoms as i32;
                pcd[io].maxlen_iso_sort_key = maxlen_iso_sort_key_no_taut_h;
                pcd[io].maxlen_iso_exchg_atnos = 0;
                let mut i;
                if num_iso_no_taut_h != 0 {
                    n_temp_rank.iter_mut().for_each(|x| *x = 0);
                    for i in 0..num_atoms {
                        let c = sym_fix_h[i] as usize - 1;
                        if (n_temp_rank[c] as usize) < i {
                            n_temp_rank[c] = i as u16;
                        }
                    }
                    i = 0;
                    while i < num_atoms {
                        if iso_sort_key_no_taut_h[i]
                            != iso_sort_key_no_taut_h
                                [n_temp_rank[sym_fix_h[i] as usize - 1] as usize]
                        {
                            break;
                        }
                        i += 1;
                    }
                } else {
                    i = num_atoms;
                }
                if i < num_atoms {
                    pcd[io].n_canon_flags |= CANON_FLAG_ISO_FIXED_H_DIFF;
                    for x in p_atom_invariant_aux.iter_mut() {
                        *x = AtomInvariant2::default();
                    }
                    for i in 0..num_atoms {
                        p_atom_invariant_aux[i].val[0] = sym_fix_h[i];
                        p_atom_invariant_aux[i].iso_sort_key = iso_sort_key_no_taut_h[i];
                    }
                    n_num_curr_ranks = set_initial_ranks2(
                        num_atoms,
                        &p_atom_invariant_aux,
                        &mut n_rank_aux,
                        &mut n_atom_number_aux,
                    );
                    n_num_curr_ranks = differentiate_ranks2(
                        num_atoms,
                        neigh_non_opt.as_mut().unwrap(),
                        n_num_curr_ranks,
                        &mut n_rank_aux,
                        &mut n_temp_rank,
                        &mut n_atom_number_aux,
                        &mut l_count,
                        0,
                    );
                    pcd[io].n_aux_rank = Some(n_rank_aux.clone());
                    let ret = canon!(
                        TAUT_NON,
                        num_atoms,
                        true,
                        sym_fix_h_iso,
                        can_fix_h_iso,
                        atn_fix_h_iso,
                        io,
                        ct_fix_h.as_ref(),
                        ct_fix_h_iso
                    );
                    if ret < 0 {
                        fail!(ret);
                    }
                } else {
                    let mut tmp = new_ct!(io);
                    ct_full_copy(&mut tmp, ct_fix_h.as_ref().unwrap());
                    for i in 0..num_atoms {
                        tmp.iso_sort_key[can_fix_h[i] as usize - 1] = iso_sort_key_no_taut_h[i];
                    }
                    tmp.len_iso_sort_key = num_atoms as i32;
                    ct_fix_h_iso = Some(tmp);
                    sym_fix_h_iso[..num_atoms].copy_from_slice(&sym_fix_h[..num_atoms]);
                    can_fix_h_iso[..num_atoms].copy_from_slice(&can_fix_h[..num_atoms]);
                    atn_fix_h_iso[..num_atoms].copy_from_slice(&atn_fix_h[..num_atoms]);
                }
            }
        }

        // ---- consistency check ----
        let ctb = |c: &ConTable, len: i32| c.ctbl[..len as usize].to_vec();
        let nh = |c: &ConTable, len: i32| c.num_h[..len as usize].to_vec();
        let c_no_h = ct_no_h.as_ref().unwrap();
        let c_no_taut_h = ct_no_taut_h.as_ref().unwrap();
        let mut cons = false;
        cons |= c_no_h.len_ct != c_no_taut_h.len_ct
            || ctb(c_no_h, c_no_h.len_ct) != ctb(c_no_taut_h, c_no_h.len_ct);
        if b_req_taut {
            let c_base = ct_base.as_ref().unwrap();
            if let Some(c_fix) = ct_fix_h.as_ref() {
                cons |= c_no_taut_h.len_ct != c_fix.len_ct
                    || ctb(c_no_taut_h, c_no_taut_h.len_ct) != ctb(c_fix, c_no_taut_h.len_ct);
                cons |= c_no_taut_h.len_num_h != c_fix.len_num_h
                    || nh(c_no_taut_h, c_no_taut_h.len_num_h) != nh(c_fix, c_no_taut_h.len_num_h);
            }
            cons |= c_no_taut_h.len_ct > c_base.len_ct
                || ctb(c_no_taut_h, c_no_taut_h.len_ct) != ctb(c_base, c_no_taut_h.len_ct);
            cons |= c_no_taut_h.len_num_h > c_base.len_num_h
                || nh(c_no_taut_h, c_no_taut_h.len_num_h) != nh(c_base, c_no_taut_h.len_num_h);
        }
        if let Some(c_iso) = ct_no_taut_h_iso.as_ref() {
            cons |= c_no_h.len_ct != c_iso.len_ct
                || ctb(c_no_h, c_no_h.len_ct) != ctb(c_iso, c_no_h.len_ct);
            cons |= c_no_taut_h.len_num_h != c_iso.len_num_h
                || nh(c_no_taut_h, c_no_taut_h.len_num_h) != nh(c_iso, c_no_taut_h.len_num_h);
        } else if let Some(c_biso) = ct_base_iso.as_ref() {
            let c_base = ct_base.as_ref().unwrap();
            cons |= c_biso.len_ct != c_base.len_ct
                || ctb(c_biso, c_biso.len_ct) != ctb(c_base, c_biso.len_ct);
            cons |= c_biso.len_num_h != c_base.len_num_h
                || nh(c_biso, c_biso.len_num_h) != nh(c_base, c_biso.len_num_h);
            if let Some(c_fiso) = ct_fix_h_iso.as_ref() {
                cons |= c_fiso.len_ct > c_biso.len_ct
                    || ctb(c_fiso, c_fiso.len_ct) != ctb(c_biso, c_fiso.len_ct);
                cons |= c_fiso.len_num_h > c_biso.len_num_h
                    || nh(c_fiso, c_fiso.len_num_h) != nh(c_biso, c_fiso.len_num_h);
            }
        }
        if cons {
            fail!(1);
        }
        if b_req_taut {
            t_group_info.unwrap().b_ignore_isotopic = b_taut_ignore_isotopic;
        }

        // ---- output ----
        bcn.num_max = num_max;
        bcn.num_at_tg = num_at_tg;
        bcn.num_atoms = num_atoms;
        bcn.ftcn[TAUT_NON].neigh_list = neigh_non_opt.take();
        bcn.ftcn[TAUT_YES].neigh_list = neigh_yes.take();
        if b_req_taut {
            let f = &mut bcn.ftcn[TAUT_YES];
            let mut c_base = ct_base.take().unwrap();
            f.num_at_tg = num_at_tg;
            f.num_atoms = num_atoms;
            f.linear_ct = std::mem::take(&mut c_base.ctbl);
            f.n_len_linear_ct_at_only = s[ib].n_len_ct_at_only;
            f.n_max_len_linear_ct = s[ib].n_len_ct + 1;
            f.n_len_linear_ct = s[ib].n_len_ct;
            f.partition_ct = Partition {
                rank: std::mem::take(&mut can_base),
                at_number: std::mem::take(&mut atn_base),
            };
            f.n_symm_rank_ct = std::mem::take(&mut sym_base);
            f.n_num_h_orig = std::mem::take(&mut num_h);
            f.n_num_h = std::mem::take(&mut c_base.num_h);
            f.n_len_num_h = maxlen_num_h.min(c_base.maxlen_num_h);
            f.n_num_h_orig_fix_h = Vec::new();
            f.n_num_h_fix_h = Vec::new();
            f.n_len_num_h_fix_h = 0;
            f.n_canon_flags |= pcd[ib].n_canon_flags;
            let l = f.n_len_num_h as usize;
            clean_num_h_if_any(&mut f.n_num_h_orig, l);
            clean_num_h_if_any(&mut f.n_num_h, l);
            clean_ct(&mut f.linear_ct, f.n_len_linear_ct as usize);
            if let Some(mut c_biso) = ct_base_iso.take() {
                f.partition_ct_iso = Partition {
                    rank: std::mem::take(&mut can_base_iso),
                    at_number: std::mem::take(&mut atn_base_iso),
                };
                f.n_symm_rank_ct_iso = std::mem::take(&mut sym_base_iso);
                f.iso_sort_keys = std::mem::take(&mut c_biso.iso_sort_key);
                f.iso_sort_keys_orig = std::mem::take(&mut iso_sort_key_base);
                f.len_iso_sort_keys = len_iso_sort_key_base;
                f.iso_exchg_atnos = std::mem::take(&mut c_biso.iso_exchg_atnos);
                f.iso_exchg_atnos_orig = std::mem::take(&mut iso_exchg_atnos);
                let l = f.len_iso_sort_keys as usize;
                clean_iso_keys_if_any(&mut f.iso_sort_keys, l);
                clean_iso_keys_if_any(&mut f.iso_sort_keys_orig, l);
            }
        }
        if b_req_non_taut {
            if !b_req_taut {
                sym_fix_h = std::mem::take(&mut sym_no_taut_h);
                can_fix_h = std::mem::take(&mut can_no_taut_h);
                atn_fix_h = std::mem::take(&mut atn_no_taut_h);
                ct_fix_h = ct_no_taut_h.take();
                sym_fix_h_iso = std::mem::take(&mut sym_no_taut_h_iso);
                can_fix_h_iso = std::mem::take(&mut can_no_taut_h_iso);
                atn_fix_h_iso = std::mem::take(&mut atn_no_taut_h_iso);
                ct_fix_h_iso = ct_no_taut_h_iso.take();
                if io == TAUT_YES
                    && bcn.ftcn[TAUT_NON].neigh_list.is_some()
                    && bcn.ftcn[TAUT_YES].neigh_list.is_none()
                {
                    bcn.ftcn[TAUT_YES].neigh_list = bcn.ftcn[TAUT_NON].neigh_list.take();
                }
            }
            let f = &mut bcn.ftcn[io];
            let mut c_fix = ct_fix_h.take().unwrap();
            f.num_at_tg = num_atoms;
            f.num_atoms = num_atoms;
            f.linear_ct = std::mem::take(&mut c_fix.ctbl);
            f.n_len_linear_ct_at_only = s[io].n_len_ct_at_only;
            f.n_max_len_linear_ct = s[io].n_len_ct_at_only + 1;
            f.n_len_linear_ct = s[io].n_len_ct_at_only;
            f.partition_ct = Partition {
                rank: std::mem::take(&mut can_fix_h),
                at_number: std::mem::take(&mut atn_fix_h),
            };
            f.n_symm_rank_ct = std::mem::take(&mut sym_fix_h);
            f.n_num_h_orig = std::mem::take(&mut num_h_no_taut_h);
            f.n_num_h = std::mem::take(&mut c_fix.num_h);
            f.n_len_num_h = maxlen_num_h_no_taut_h.min(c_fix.maxlen_num_h);
            f.n_num_h_orig_fix_h = std::mem::take(&mut num_h_fixed);
            f.n_num_h_fix_h = std::mem::take(&mut c_fix.num_h_fixed);
            f.n_len_num_h_fix_h = num_atoms as i32;
            f.n_canon_flags |= pcd[io].n_canon_flags;
            let l = f.n_len_num_h as usize;
            clean_num_h_if_any(&mut f.n_num_h_orig, l);
            clean_num_h_if_any(&mut f.n_num_h_orig_fix_h, l);
            clean_num_h_if_any(&mut f.n_num_h, l);
            clean_num_h_if_any(&mut f.n_num_h_fix_h, l);
            clean_ct(&mut f.linear_ct, f.n_len_linear_ct as usize);
            if let Some(mut c_fiso) = ct_fix_h_iso.take() {
                f.partition_ct_iso = Partition {
                    rank: std::mem::take(&mut can_fix_h_iso),
                    at_number: std::mem::take(&mut atn_fix_h_iso),
                };
                f.n_symm_rank_ct_iso = std::mem::take(&mut sym_fix_h_iso);
                f.iso_sort_keys = std::mem::take(&mut c_fiso.iso_sort_key);
                f.iso_sort_keys_orig = std::mem::take(&mut iso_sort_key_no_taut_h);
                f.len_iso_sort_keys = len_iso_sort_key_no_taut_h;
                let l = f.len_iso_sort_keys as usize;
                clean_iso_keys_if_any(&mut f.iso_sort_keys, l);
                clean_iso_keys_if_any(&mut f.iso_sort_keys_orig, l);
            }
        }
        0
    };
    // the first two arrays of the rank stack belong to this function
    if !pi.is_empty() {
        pi[0] = Partition::default();
    }
    bcn.p_rank_stack = pi;
    ret
}
