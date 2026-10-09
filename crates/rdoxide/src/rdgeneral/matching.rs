//! Maximum cardinality matching in general graphs (Edmonds' blossom algorithm).
//!
//! Replaces `boost::edmonds_maximum_cardinality_matching`, which RDKit uses
//! in DetermineBonds. A maximum matching is not unique: this implementation
//! returns *a* maximum matching, which may differ from the one boost finds.

const NONE: usize = usize::MAX;

struct Blossoms<'a> {
    adj: &'a [Vec<usize>],
    mate: Vec<usize>,
    parent: Vec<usize>,
    base: Vec<usize>,
    in_queue: Vec<bool>,
    in_blossom: Vec<bool>,
}

impl Blossoms<'_> {
    fn lowest_common_ancestor(&self, mut a: usize, mut b: usize) -> usize {
        let mut seen = vec![false; self.adj.len()];
        loop {
            a = self.base[a];
            seen[a] = true;
            if self.mate[a] == NONE {
                break;
            }
            a = self.parent[self.mate[a]];
        }
        loop {
            b = self.base[b];
            if seen[b] {
                return b;
            }
            b = self.parent[self.mate[b]];
        }
    }

    fn mark_path(&mut self, mut v: usize, base: usize, mut child: usize) {
        while self.base[v] != base {
            self.in_blossom[self.base[v]] = true;
            self.in_blossom[self.base[self.mate[v]]] = true;
            self.parent[v] = child;
            child = self.mate[v];
            v = self.parent[self.mate[v]];
        }
    }

    /// An augmenting path from the free vertex `root`, as its last vertex.
    fn find_path(&mut self, root: usize) -> Option<usize> {
        let n = self.adj.len();
        self.in_queue.fill(false);
        self.parent.fill(NONE);
        for (i, b) in self.base.iter_mut().enumerate() {
            *b = i;
        }
        self.in_queue[root] = true;
        let mut queue = std::collections::VecDeque::from([root]);
        while let Some(v) = queue.pop_front() {
            for k in 0..self.adj[v].len() {
                let to = self.adj[v][k];
                if self.base[v] == self.base[to] || self.mate[v] == to {
                    continue;
                }
                if to == root || (self.mate[to] != NONE && self.parent[self.mate[to]] != NONE) {
                    // an odd cycle: contract the blossom
                    let current_base = self.lowest_common_ancestor(v, to);
                    self.in_blossom.fill(false);
                    self.mark_path(v, current_base, to);
                    self.mark_path(to, current_base, v);
                    for i in 0..n {
                        if self.in_blossom[self.base[i]] {
                            self.base[i] = current_base;
                            if !self.in_queue[i] {
                                self.in_queue[i] = true;
                                queue.push_back(i);
                            }
                        }
                    }
                } else if self.parent[to] == NONE {
                    self.parent[to] = v;
                    if self.mate[to] == NONE {
                        return Some(to);
                    }
                    let next = self.mate[to];
                    self.in_queue[next] = true;
                    queue.push_back(next);
                }
            }
        }
        None
    }
}

/// A maximum matching of the graph with `n` vertices and the given edges:
/// `mate[v]` is the partner of `v`, if any. Self loops are ignored.
pub fn max_cardinality_matching(n: usize, edges: &[(usize, usize)]) -> Vec<Option<usize>> {
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in edges {
        if a != b {
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    let mut state = Blossoms {
        adj: &adj,
        mate: vec![NONE; n],
        parent: vec![NONE; n],
        base: (0..n).collect(),
        in_queue: vec![false; n],
        in_blossom: vec![false; n],
    };
    // greedy initial matching
    for v in 0..n {
        if state.mate[v] == NONE
            && let Some(&u) = adj[v].iter().find(|&&u| state.mate[u] == NONE)
        {
            state.mate[v] = u;
            state.mate[u] = v;
        }
    }
    for root in 0..n {
        if state.mate[root] != NONE {
            continue;
        }
        let Some(mut v) = state.find_path(root) else {
            continue;
        };
        while v != NONE {
            let pv = state.parent[v];
            let next = state.mate[pv];
            state.mate[v] = pv;
            state.mate[pv] = v;
            v = next;
        }
    }
    state
        .mate
        .into_iter()
        .map(|m| (m != NONE).then_some(m))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn size(mate: &[Option<usize>]) -> usize {
        mate.iter().flatten().count() / 2
    }

    #[test]
    fn path_needs_augmentation() {
        // 0-1-2-3: a greedy start (0-1, 2-3) is already perfect; 1-0 / 1-2 / 2-3 forces augmenting
        let mate = max_cardinality_matching(4, &[(1, 2), (0, 1), (2, 3)]);
        assert_eq!(size(&mate), 2);
    }

    #[test]
    fn odd_cycle_with_tail_uses_a_blossom() {
        // triangle 0-1-2 plus pendant vertices 3 (on 0) and 4 (on 1)
        let mate = max_cardinality_matching(5, &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 4)]);
        assert_eq!(size(&mate), 2);
        assert!(mate[3].is_some() || mate[4].is_some());
    }

    #[test]
    fn five_cycle() {
        let mate = max_cardinality_matching(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
        assert_eq!(size(&mate), 2);
    }

    #[test]
    fn matching_is_consistent() {
        let edges: Vec<(usize, usize)> = (0..9).map(|i| (i, (i * 7 + 3) % 10)).collect();
        let mate = max_cardinality_matching(10, &edges);
        for (v, m) in mate.iter().enumerate() {
            if let Some(m) = m {
                assert_eq!(mate[*m], Some(v));
            }
        }
    }
}
