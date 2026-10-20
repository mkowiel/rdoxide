//! Sorting primitives of the InChI library (`INCHI_BASE/src/ichisort.c`).
//!
//! InChI does not use the C library's `qsort`: the canonical numbering
//! depends on the exact order in which an unstable sort leaves equal
//! elements, so `inchi_qsort` and the insertion sorts are translated
//! literally (and tested against the C code, see `tests/inchi_utils.rs`).
//!
//! Pointers into the arrays become indices; the comparison functions return
//! an `i32` like the C ones. `NEIGH_LIST` (an `AT_NUMB *` whose first element
//! is the number of neighbors) is a `&[u16]`/`&mut [u16]` with that layout.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

/// Number of entries of the explicit stack of `inchi_qsort`
/// (`8 * sizeof(void *) - 2` on a 64-bit platform).
const STKSIZ: usize = 62;

/// `inchi_qsort`: quicksort (median of three, three-way partition) of
/// `base`; `comp(a, b)` is negative, zero or positive like the C comparison
/// callback. The sort is not stable.
pub fn inchi_qsort<T, F>(base: &mut [T], mut comp: F)
where
    F: FnMut(&T, &T) -> i32,
{
    let num = base.len();
    if num < 2 {
        return;
    }
    // indices stand for the pointers `lo`, `hi`, `mid`, `loguy` and `higuy`;
    // `loguy` can reach `hi + 1` and `higuy` can reach `lo - 1`, so they are signed
    let mut lostk = [0isize; STKSIZ];
    let mut histk = [0isize; STKSIZ];
    let mut stkptr: isize = 0;
    let mut lo: isize = 0;
    let mut hi: isize = num as isize - 1;

    macro_rules! cmp {
        ($a:expr, $b:expr) => {
            comp(&base[$a as usize], &base[$b as usize])
        };
    }
    macro_rules! swap {
        ($a:expr, $b:expr) => {
            if $a != $b {
                base.swap($a as usize, $b as usize);
            }
        };
    }

    loop {
        // `recurse:`
        let size = hi - lo + 1;

        // median of the first, middle and last elements
        let mut mid = lo + size / 2;
        if cmp!(lo, mid) > 0 {
            swap!(lo, mid);
        }
        if cmp!(lo, hi) > 0 {
            swap!(lo, hi);
        }
        if cmp!(mid, hi) > 0 {
            swap!(mid, hi);
        }

        let mut loguy = lo;
        let mut higuy = hi;
        loop {
            if mid > loguy {
                loop {
                    loguy += 1;
                    if !(loguy < mid && cmp!(loguy, mid) <= 0) {
                        break;
                    }
                }
            }
            if mid <= loguy {
                loop {
                    loguy += 1;
                    if !(loguy <= hi && cmp!(loguy, mid) <= 0) {
                        break;
                    }
                }
            }
            loop {
                higuy -= 1;
                if !(higuy > mid && cmp!(higuy, mid) > 0) {
                    break;
                }
            }
            if higuy < loguy {
                break;
            }
            swap!(loguy, higuy);
            // follow the partition element if it was moved
            if mid == higuy {
                mid = loguy;
            }
        }

        // adjacent elements equal to the partition element
        higuy += 1;
        if mid < higuy {
            loop {
                higuy -= 1;
                if !(higuy > mid && cmp!(higuy, mid) == 0) {
                    break;
                }
            }
        }
        if mid >= higuy {
            loop {
                higuy -= 1;
                if !(higuy > lo && cmp!(higuy, mid) == 0) {
                    break;
                }
            }
        }

        // sort the smaller sub-array first, stack the bigger one
        let mut again = false;
        if higuy - lo >= hi - loguy {
            if lo < higuy {
                lostk[stkptr as usize] = lo;
                histk[stkptr as usize] = higuy;
                stkptr += 1;
            }
            if loguy < hi {
                lo = loguy;
                again = true;
            }
        } else {
            if loguy < hi {
                lostk[stkptr as usize] = loguy;
                histk[stkptr as usize] = hi;
                stkptr += 1;
            }
            if lo < higuy {
                hi = higuy;
                again = true;
            }
        }
        if again {
            continue;
        }
        stkptr -= 1;
        if stkptr >= 0 {
            lo = lostk[stkptr as usize];
            hi = histk[stkptr as usize];
        } else {
            return;
        }
    }
}

/// `insertions_sort`: insertion sort by swaps; returns the number of swaps.
pub fn insertions_sort<T, F>(base: &mut [T], mut compare: F) -> i32
where
    F: FnMut(&T, &T) -> i32,
{
    let mut num_trans = 0;
    for k in 1..base.len() {
        // `i` and `j` are neighbors, `j` moves down while it is smaller
        let mut j = k;
        while j > 0 && compare(&base[j - 1], &base[j]) > 0 {
            base.swap(j - 1, j);
            num_trans += 1;
            j -= 1;
        }
    }
    num_trans
}

/// `insertions_sort_AT_NUMBERS`: insertion sort by shifting, with `compare`
/// receiving the element and the element being inserted.
pub fn insertions_sort_at_numbers<F>(base: &mut [u16], mut compare: F) -> i32
where
    F: FnMut(&u16, &u16) -> i32,
{
    let mut num_trans = 0;
    for k in 1..base.len() {
        let tmp = base[k];
        let mut j = k;
        while j > 0 && compare(&base[j - 1], &tmp) > 0 {
            base[j] = base[j - 1];
            num_trans += 1;
            j -= 1;
        }
        base[j] = tmp;
    }
    num_trans
}

/// `insertions_sort_NeighList_AT_NUMBERS`: sorts the neighbors of a
/// `NEIGH_LIST` (`base[0]` is their number) by ascending rank.
pub fn insertions_sort_neigh_list_at_numbers(base: &mut [u16], n_rank: &[u16]) {
    let num = base[0] as usize;
    let nb = &mut base[1..];
    for k in 1..num {
        let rj = n_rank[nb[k] as usize];
        let mut j = k;
        while j > 0 && n_rank[nb[j - 1] as usize] > rj {
            nb.swap(j - 1, j);
            j -= 1;
        }
    }
}

/// `insertions_sort_AT_RANK`: ascending insertion sort; returns the number of shifts.
pub fn insertions_sort_at_rank(base: &mut [u16]) -> i32 {
    let mut num_trans = 0;
    for k in 1..base.len() {
        let tmp = base[k];
        let mut j = k;
        while j > 0 && base[j - 1] > tmp {
            base[j] = base[j - 1];
            num_trans += 1;
            j -= 1;
        }
        base[j] = tmp;
    }
    num_trans
}

/// `insertions_sort_NeighList_AT_NUMBERS3`: like
/// [`insertions_sort_neigh_list_at_numbers`], shifting instead of swapping,
/// and returning the number of shifts.
pub fn insertions_sort_neigh_list_at_numbers3(base: &mut [u16], n_rank: &[u16]) -> i32 {
    let num = base[0] as usize;
    let nb = &mut base[1..];
    let mut n = 0;
    for k in 1..num {
        let tmp = nb[k];
        let rj = n_rank[tmp as usize];
        let mut j = k;
        while j > 0 && n_rank[nb[j - 1] as usize] > rj {
            nb[j] = nb[j - 1];
            n += 1;
            j -= 1;
        }
        nb[j] = tmp;
    }
    n
}

/// `insertions_sort_NeighListBySymmAndCanonRank`: sorts the neighbors by
/// descending symmetry rank, then ascending canonical rank.
pub fn insertions_sort_neigh_list_by_symm_and_canon_rank(
    base: &mut [u16],
    n_symm_rank: &[u16],
    n_canon_rank: &[u16],
) {
    let num = base[0] as usize;
    let nb = &mut base[1..];
    for k in 1..num {
        let mut j = k;
        while j > 0 {
            let (a, b) = (nb[j - 1] as usize, nb[j] as usize);
            let diff = i32::from(n_symm_rank[a]) - i32::from(n_symm_rank[b]);
            if diff < 0 || (diff == 0 && n_canon_rank[a] < n_canon_rank[b]) {
                nb.swap(j - 1, j);
                j -= 1;
            } else {
                break;
            }
        }
    }
}

/// `comp_AT_RANK`
pub fn comp_at_rank(a1: &u16, a2: &u16) -> i32 {
    i32::from(*a1) - i32::from(*a2)
}

/// `CompChemElemLex`: lexicographic comparison of two element symbols
/// (`memcmp` of 2 bytes).
pub fn comp_chem_elem_lex(a1: &[u8; 2], a2: &[u8; 2]) -> i32 {
    match a1.cmp(a2) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// `CompareNeighListLex`: lexicographic comparison of two `NEIGH_LIST`s by
/// the ranks of their neighbors.
pub fn compare_neigh_list_lex(pp1: &[u16], pp2: &[u16], n_rank: &[u16]) -> i32 {
    let (len1, len2) = (i32::from(pp1[0]), i32::from(pp2[0]));
    let mut len = len1.min(len2);
    let (mut i1, mut i2) = (1usize, 1usize);
    let mut diff = 0;
    while len > 0 && diff == 0 {
        len -= 1;
        diff = i32::from(n_rank[pp1[i1] as usize]) - i32::from(n_rank[pp2[i2] as usize]);
        i1 += 1;
        i2 += 1;
    }
    if diff != 0 { diff } else { len1 - len2 }
}

/// `CompareNeighListLexUpToMaxRank`: like [`compare_neigh_list_lex`], ignoring
/// the trailing neighbors whose rank exceeds `n_max_at_neigh_rank`.
pub fn compare_neigh_list_lex_up_to_max_rank(
    pp1: &[u16],
    pp2: &[u16],
    n_rank: &[u16],
    n_max_at_neigh_rank: u16,
) -> i32 {
    let (mut len1, mut len2) = (i32::from(pp1[0]), i32::from(pp2[0]));
    let (a, b) = (&pp1[1..], &pp2[1..]);
    while 0 < len1 && n_rank[a[(len1 - 1) as usize] as usize] > n_max_at_neigh_rank {
        len1 -= 1;
    }
    while 0 < len2 && n_rank[b[(len2 - 1) as usize] as usize] > n_max_at_neigh_rank {
        len2 -= 1;
    }
    let mut len = len1.min(len2);
    let (mut i1, mut i2) = (0usize, 0usize);
    let mut diff = 0;
    while len > 0 {
        len -= 1;
        diff = i32::from(n_rank[a[i1] as usize]) - i32::from(n_rank[b[i2] as usize]);
        i1 += 1;
        i2 += 1;
        if diff != 0 {
            break;
        }
    }
    if diff != 0 { diff } else { len1 - len2 }
}

/// `iisort`: ascending insertion sort of a list of ints.
pub fn iisort(list: &mut [i32]) {
    for k in 1..list.len() {
        let tmp = list[k];
        let mut j = k;
        while j > 0 && list[j - 1] > tmp {
            list[j] = list[j - 1];
            j -= 1;
        }
        list[j] = tmp;
    }
}
