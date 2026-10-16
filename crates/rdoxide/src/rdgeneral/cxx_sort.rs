//! `std::sort` as implemented by LLVM's libc++ (17 and later), for a
//! comparator that is not `std::less` (the classic, non-branchless path).
//!
//! `std::sort` is not stable, so where RDKit sorts with ties the result
//! depends on the C++ standard library. The reference builds (conda-forge,
//! macOS) use libc++; this reproduces its exact element order.
//!
//! Ported from libc++ `__algorithm/sort.h`, `partial_sort.h`, `make_heap.h`,
//! `sort_heap.h`, `pop_heap.h`, `sift_down.h` and `push_heap.h`.
//! Part of the LLVM Project, under the Apache License v2.0 with LLVM
//! Exceptions (see LICENSES/LLVM-Apache-2.0-with-LLVM-exception.txt).

/// `__sort3`
fn sort3<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], x: usize, y: usize, z: usize, c: &mut F) {
    if !c(&v[y], &v[x]) {
        if !c(&v[z], &v[y]) {
            return;
        }
        v.swap(y, z);
        if c(&v[y], &v[x]) {
            v.swap(x, y);
        }
        return;
    }
    if c(&v[z], &v[y]) {
        v.swap(x, z);
        return;
    }
    v.swap(x, y);
    if c(&v[z], &v[y]) {
        v.swap(y, z);
    }
}

/// `__sort4`
fn sort4<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    x1: usize,
    x2: usize,
    x3: usize,
    x4: usize,
    c: &mut F,
) {
    sort3(v, x1, x2, x3, c);
    if c(&v[x4], &v[x3]) {
        v.swap(x3, x4);
        if c(&v[x3], &v[x2]) {
            v.swap(x2, x3);
            if c(&v[x2], &v[x1]) {
                v.swap(x1, x2);
            }
        }
    }
}

/// `__sort5`
fn sort5<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], x: [usize; 5], c: &mut F) {
    sort4(v, x[0], x[1], x[2], x[3], c);
    if c(&v[x[4]], &v[x[3]]) {
        v.swap(x[3], x[4]);
        if c(&v[x[3]], &v[x[2]]) {
            v.swap(x[2], x[3]);
            if c(&v[x[2]], &v[x[1]]) {
                v.swap(x[1], x[2]);
                if c(&v[x[1]], &v[x[0]]) {
                    v.swap(x[0], x[1]);
                }
            }
        }
    }
}

/// `__insertion_sort` on `[first, last)`; `guarded == false` is
/// `__insertion_sort_unguarded` (an element at `first - 1` bounds the scan).
fn insertion_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    c: &mut F,
    guarded: bool,
) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        let mut j = i - 1;
        if c(&v[i], &v[j]) {
            let t = v[i].clone();
            let mut k = j;
            j = i;
            loop {
                v[j] = v[k].clone();
                j = k;
                if guarded && j == first {
                    break;
                }
                k -= 1;
                if !c(&t, &v[k]) {
                    break;
                }
            }
            v[j] = t;
        }
    }
}

/// `__insertion_sort_incomplete`: gives up after 8 moved elements.
fn insertion_sort_incomplete<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    c: &mut F,
) -> bool {
    match last - first {
        0 | 1 => return true,
        2 => {
            if c(&v[last - 1], &v[first]) {
                v.swap(first, last - 1);
            }
            return true;
        }
        3 => {
            sort3(v, first, first + 1, last - 1, c);
            return true;
        }
        4 => {
            sort4(v, first, first + 1, first + 2, last - 1, c);
            return true;
        }
        5 => {
            sort5(v, [first, first + 1, first + 2, first + 3, last - 1], c);
            return true;
        }
        _ => {}
    }
    let mut j = first + 2;
    sort3(v, first, first + 1, j, c);
    let mut count = 0;
    let mut i = j + 1;
    while i != last {
        if c(&v[i], &v[j]) {
            let t = v[i].clone();
            let mut k = j;
            j = i;
            loop {
                v[j] = v[k].clone();
                j = k;
                if j == first {
                    break;
                }
                k -= 1;
                if !c(&t, &v[k]) {
                    break;
                }
            }
            v[j] = t;
            count += 1;
            if count == 8 {
                return i + 1 == last;
            }
        }
        j = i;
        i += 1;
    }
    true
}

/// `__partition_with_equals_on_right`: returns (pivot position, already
/// partitioned).
fn partition_with_equals_on_right<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    c: &mut F,
) -> (usize, bool) {
    let mut first = begin;
    let mut last = end;
    let pivot = v[first].clone();
    loop {
        first += 1;
        if !c(&v[first], &pivot) {
            break;
        }
    }
    if begin == first - 1 {
        while first < last {
            last -= 1;
            if c(&v[last], &pivot) {
                break;
            }
        }
    } else {
        loop {
            last -= 1;
            if c(&v[last], &pivot) {
                break;
            }
        }
    }
    let already = first >= last;
    while first < last {
        v.swap(first, last);
        loop {
            first += 1;
            if !c(&v[first], &pivot) {
                break;
            }
        }
        loop {
            last -= 1;
            if c(&v[last], &pivot) {
                break;
            }
        }
    }
    let pivot_pos = first - 1;
    if begin != pivot_pos {
        v[begin] = v[pivot_pos].clone();
    }
    v[pivot_pos] = pivot;
    (pivot_pos, already)
}

/// `__partition_with_equals_on_left`: returns the new `first`.
fn partition_with_equals_on_left<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    c: &mut F,
) -> usize {
    let mut first = begin;
    let mut last = end;
    let pivot = v[first].clone();
    if c(&pivot, &v[last - 1]) {
        loop {
            first += 1;
            if c(&pivot, &v[first]) {
                break;
            }
        }
    } else {
        loop {
            first += 1;
            if !(first < last && !c(&pivot, &v[first])) {
                break;
            }
        }
    }
    if first < last {
        loop {
            last -= 1;
            if !c(&pivot, &v[last]) {
                break;
            }
        }
    }
    while first < last {
        v.swap(first, last);
        loop {
            first += 1;
            if c(&pivot, &v[first]) {
                break;
            }
        }
        loop {
            last -= 1;
            if !c(&pivot, &v[last]) {
                break;
            }
        }
    }
    let pivot_pos = first - 1;
    if begin != pivot_pos {
        v[begin] = v[pivot_pos].clone();
    }
    v[pivot_pos] = pivot;
    first
}

/// `__sift_down` (heap rooted at `first`, element `start`)
fn sift_down<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    c: &mut F,
    len: usize,
    start: usize,
) {
    let mut start = start;
    let mut child = start - first;
    if len < 2 || (len - 2) / 2 < child {
        return;
    }
    child = 2 * child + 1;
    let mut child_i = first + child;
    if child + 1 < len && c(&v[child_i], &v[child_i + 1]) {
        child_i += 1;
        child += 1;
    }
    if c(&v[child_i], &v[start]) {
        return;
    }
    let top = v[start].clone();
    loop {
        v[start] = v[child_i].clone();
        start = child_i;
        if (len - 2) / 2 < child {
            break;
        }
        child = 2 * child + 1;
        child_i = first + child;
        if child + 1 < len && c(&v[child_i], &v[child_i + 1]) {
            child_i += 1;
            child += 1;
        }
        if c(&v[child_i], &top) {
            break;
        }
    }
    v[start] = top;
}

/// `__floyd_sift_down`: returns the hole.
fn floyd_sift_down<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    c: &mut F,
    len: usize,
) -> usize {
    let mut hole = first;
    let mut child_i = first;
    let mut child = 0usize;
    loop {
        child_i += child + 1;
        child = 2 * child + 1;
        if child + 1 < len && c(&v[child_i], &v[child_i + 1]) {
            child_i += 1;
            child += 1;
        }
        v[hole] = v[child_i].clone();
        hole = child_i;
        if child > (len - 2) / 2 {
            return hole;
        }
    }
}

/// `__sift_up` of the element before `last`
fn sift_up<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    c: &mut F,
    len: usize,
) {
    if len > 1 {
        let mut len = (len - 2) / 2;
        let mut ptr = first + len;
        let mut last = last - 1;
        if c(&v[ptr], &v[last]) {
            let t = v[last].clone();
            loop {
                v[last] = v[ptr].clone();
                last = ptr;
                if len == 0 {
                    break;
                }
                len = (len - 1) / 2;
                ptr = first + len;
                if !c(&v[ptr], &t) {
                    break;
                }
            }
            v[last] = t;
        }
    }
}

/// `std::priority_queue` as implemented by libc++: the order in which
/// elements that compare equal come out depends on the heap algorithms.
/// `less` is the `Compare` of the queue, so the top is the element that no
/// other element is "greater" than by it.
pub struct LibcxxPriorityQueue<T, F> {
    heap: Vec<T>,
    less: F,
}

impl<T: std::fmt::Debug, F> std::fmt::Debug for LibcxxPriorityQueue<T, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LibcxxPriorityQueue")
            .field("heap", &self.heap)
            .finish_non_exhaustive()
    }
}

impl<T: Clone, F: FnMut(&T, &T) -> bool> LibcxxPriorityQueue<T, F> {
    pub fn new(less: F) -> Self {
        LibcxxPriorityQueue {
            heap: Vec::new(),
            less,
        }
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }

    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    pub fn top(&self) -> Option<&T> {
        self.heap.first()
    }

    /// `push`: `push_back` and `push_heap`
    pub fn push(&mut self, value: T) {
        self.heap.push(value);
        let len = self.heap.len();
        sift_up(&mut self.heap, 0, len, &mut self.less, len);
    }

    /// `pop`: `pop_heap` and `pop_back`
    pub fn pop(&mut self) {
        let len = self.heap.len();
        if len > 1 {
            let top = self.heap[0].clone();
            let hole = floyd_sift_down(&mut self.heap, 0, &mut self.less, len);
            let last = len - 1;
            if hole == last {
                self.heap[hole] = top;
            } else {
                self.heap[hole] = self.heap[last].clone();
                self.heap[last] = top;
                sift_up(&mut self.heap, 0, hole + 1, &mut self.less, hole + 1);
            }
        }
        self.heap.pop();
    }
}

/// `__partial_sort(first, last, last)`: heap sort.
fn heap_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    c: &mut F,
) {
    if first == last {
        return;
    }
    let n = last - first;
    if n > 1 {
        let mut start = (n - 2) / 2;
        loop {
            sift_down(v, first, c, n, first + start);
            if start == 0 {
                break;
            }
            start -= 1;
        }
    }
    // __sort_heap
    let mut last = last;
    let mut len = n;
    while len > 1 {
        let top = v[first].clone();
        let hole = floyd_sift_down(v, first, c, len);
        let l = last - 1;
        if hole == l {
            v[hole] = top;
        } else {
            v[hole] = v[l].clone();
            v[l] = top;
            sift_up(v, first, hole + 1, c, hole + 1 - first);
        }
        last -= 1;
        len -= 1;
    }
}

/// `__introsort`
fn introsort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    c: &mut F,
    depth: usize,
    leftmost: bool,
) {
    const LIMIT: usize = 24;
    const NINTHER_THRESHOLD: usize = 128;
    let (mut first, mut last, mut depth, mut leftmost) = (first, last, depth, leftmost);
    loop {
        let len = last - first;
        match len {
            0 | 1 => return,
            2 => {
                if c(&v[last - 1], &v[first]) {
                    v.swap(first, last - 1);
                }
                return;
            }
            3 => {
                sort3(v, first, first + 1, last - 1, c);
                return;
            }
            4 => {
                sort4(v, first, first + 1, first + 2, last - 1, c);
                return;
            }
            5 => {
                sort5(v, [first, first + 1, first + 2, first + 3, last - 1], c);
                return;
            }
            _ => {}
        }
        if len < LIMIT {
            insertion_sort(v, first, last, c, leftmost);
            return;
        }
        if depth == 0 {
            heap_sort(v, first, last, c);
            return;
        }
        depth -= 1;
        let half = len / 2;
        if len > NINTHER_THRESHOLD {
            sort3(v, first, first + half, last - 1, c);
            sort3(v, first + 1, first + (half - 1), last - 2, c);
            sort3(v, first + 2, first + (half + 1), last - 3, c);
            sort3(v, first + (half - 1), first + half, first + (half + 1), c);
            v.swap(first, first + half);
        } else {
            sort3(v, first + half, first, last - 1, c);
        }
        if !leftmost && !c(&v[first - 1], &v[first]) {
            first = partition_with_equals_on_left(v, first, last, c);
            continue;
        }
        let (i, already) = partition_with_equals_on_right(v, first, last, c);
        if already {
            let fs = insertion_sort_incomplete(v, first, i, c);
            if insertion_sort_incomplete(v, i + 1, last, c) {
                if fs {
                    return;
                }
                last = i;
                continue;
            } else if fs {
                first = i + 1;
                continue;
            }
        }
        introsort(v, first, i, c, depth, leftmost);
        leftmost = false;
        first = i + 1;
    }
}

/// `std::sort(v.begin(), v.end(), less)` with libc++'s element order.
pub fn libcxx_sort_by<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut less: F) {
    let n = v.len();
    let depth = if n == 0 {
        0
    } else {
        2 * (usize::BITS - 1 - n.leading_zeros()) as usize
    };
    introsort(v, 0, n, &mut less, depth, true);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(n: usize, seed: u64) {
        let mut x = seed;
        let mut v: Vec<(u32, usize)> = (0..n)
            .map(|i| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (((x >> 33) % 7) as u32, i)
            })
            .collect();
        libcxx_sort_by(&mut v, |a, b| a.0 < b.0);
        assert!(v.windows(2).all(|w| w[0].0 <= w[1].0), "n={n}");
        let mut idx: Vec<usize> = v.iter().map(|p| p.1).collect();
        idx.sort_unstable();
        assert_eq!(idx, (0..n).collect::<Vec<_>>());
    }

    #[test]
    fn sorts() {
        for n in [0, 1, 2, 3, 4, 5, 6, 23, 24, 25, 100, 129, 500, 3000] {
            for seed in 0..5 {
                check(n, seed);
            }
        }
    }

    #[test]
    fn small_inputs_are_stable() {
        // insertion sort below 24 elements keeps ties in input order
        let mut v: Vec<(u32, usize)> = (0..20).map(|i| ((i % 3) as u32, i)).collect();
        libcxx_sort_by(&mut v, |a, b| a.0 < b.0);
        let mut w = v.clone();
        w.sort_by_key(|p| p.0);
        assert_eq!(v, w);
    }
}
