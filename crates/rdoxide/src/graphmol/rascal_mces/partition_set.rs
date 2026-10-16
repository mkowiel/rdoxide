//! `PartitionSet`: the vertices of the modular product graph grouped into
//! partitions of vertices that cannot be in a clique together, which gives the
//! upper bound used to prune the clique search.

use std::sync::Arc;

/// A symmetric boolean matrix stored as bits.
#[derive(Debug, Clone)]
pub struct BitMatrix {
    rows: Vec<Vec<u64>>,
}

impl BitMatrix {
    pub fn new(n: usize) -> Self {
        BitMatrix {
            rows: vec![vec![0; n.div_ceil(64)]; n],
        }
    }

    pub fn get(&self, i: usize, j: usize) -> bool {
        self.rows[i][j / 64] >> (j % 64) & 1 == 1
    }

    pub fn set(&mut self, i: usize, j: usize) {
        self.rows[i][j / 64] |= 1 << (j % 64);
    }
}

/// What all the partition sets of a search share.
#[derive(Debug)]
pub struct SharedGraph {
    pub mod_prod: BitMatrix,
    pub vtx_pairs: Vec<(i32, i32)>,
    pub vtx1_labels: Vec<u32>,
    pub vtx2_labels: Vec<u32>,
}

#[derive(Debug, Clone)]
pub struct PartitionSet {
    graph: Arc<SharedGraph>,
    parts: Vec<Vec<u32>>,
    /// counts of the number of times each vertex appears in the partitions
    vtx1_counts: Vec<i32>,
    vtx2_counts: Vec<i32>,
    /// counts of the number of times the labels appear in the partitions
    vtx1_type_counts: Vec<i32>,
    vtx2_type_counts: Vec<i32>,
}

impl PartitionSet {
    /// Makes a partition set from the modular product and the labels of the
    /// vertices of the 2 line graphs. Each pair has a row/column in the
    /// product. The partitions are sorted into descending order of size.
    pub fn new(graph: Arc<SharedGraph>, lower_bound: u32) -> Self {
        let mut ps = PartitionSet {
            parts: Vec::new(),
            vtx1_counts: vec![0; graph.vtx1_labels.len()],
            vtx2_counts: vec![0; graph.vtx2_labels.len()],
            vtx1_type_counts: Vec::new(),
            vtx2_type_counts: Vec::new(),
            graph: Arc::clone(&graph),
        };
        let mut first_vtx = -1;
        // Clearly, a vertex in one of the line graphs can only match one
        // vertex in the other. Thus, the initial partitions can be set up so
        // that all vertices in a partition have the same vertex in the first
        // line graph.
        for (i, vp) in graph.vtx_pairs.iter().enumerate() {
            if vp.0 != first_vtx {
                ps.parts.push(vec![i as u32]);
                first_vtx = vp.0;
            } else {
                ps.parts
                    .last_mut()
                    .expect("a partition was started")
                    .push(i as u32);
            }
            ps.vtx1_counts[vp.0 as usize] += 1;
            ps.vtx2_counts[vp.1 as usize] += 1;
        }
        if ps.parts.is_empty() {
            return ps;
        }
        // Now sort the partitions by size. This means that the vertices at the
        // top of the partition set, above the lower bound (or Pex as Raymond
        // calls it in the paper), are the ones that match the least number of
        // vertices in the other line graph. This has a dramatic effect on the
        // speed compared with other things tried.
        ps.sort_partitions();
        // Now reassign vertices from above Pex to below it if possible. This
        // also improves the speed of finding a large clique early. A vertex
        // is moved to a partition where it isn't connected to a vertex in the
        // modular product graph that is in the partition.
        let lb = lower_bound as usize;
        let mut i = ps.parts.len() - 1;
        while i > lb {
            let mut reassigned = false;
            for vi in 0..ps.parts[i].len() {
                let iv = ps.parts[i][vi];
                for k in 0..=lb {
                    let conn = ps.parts[k]
                        .iter()
                        .any(|&kv| graph.mod_prod.get(iv as usize, kv as usize));
                    if !conn {
                        ps.parts[k].push(iv);
                        ps.parts[i][vi] = u32::MAX;
                        reassigned = true;
                        break;
                    }
                }
            }
            if reassigned {
                ps.parts[i].retain(|&v| v != u32::MAX);
            }
            i -= 1;
        }
        ps.parts.retain(|p| !p.is_empty());
        // Sort again, to make sure the large partitions are dealt with as late
        // as possible.
        ps.sort_partitions();

        // Get the info together for the upper bound calculation.
        ps.calc_vtx_type_counts();
        ps
    }

    /// The vertex pairs (bonds of the 2 molecules) of the product.
    pub fn graph_pairs(&self) -> Vec<(i32, i32)> {
        self.graph.vtx_pairs.clone()
    }

    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }

    pub fn num_parts(&self) -> usize {
        self.parts.len()
    }

    /// The upper bound on the clique that can be extracted from the current
    /// partitions.
    pub fn upper_bound(&self) -> i32 {
        self.vtx1_type_counts
            .iter()
            .zip(&self.vtx2_type_counts)
            .map(|(a, b)| (*a).min(*b))
            .sum()
    }

    /// Removes the last element of the last partition and returns it.
    pub fn pop_last_vertex(&mut self) -> u32 {
        let part = self.parts.last_mut().expect("PartitionSet set is empty.");
        let ret = part.pop().expect("partitions are never empty");
        if part.is_empty() {
            self.parts.pop();
        }
        self.decrement_vertex_counts(ret as usize);
        ret
    }

    /// Removes from the partitions any vertex not connected to the given one.
    pub fn prune_vertices(&mut self, vtx_num: u32) {
        for p in 0..self.parts.len() {
            let mut i = 0;
            while i < self.parts[p].len() {
                let v = self.parts[p][i];
                if !self.graph.mod_prod.get(v as usize, vtx_num as usize) {
                    self.decrement_vertex_counts(v as usize);
                    let part = &mut self.parts[p];
                    part[i] = *part.last().expect("not empty");
                    part.pop();
                } else {
                    i += 1;
                }
            }
        }
        self.parts.retain(|p| !p.is_empty());
        self.sort_partitions();
    }

    /// Sorts by size, descending. Upstream breaks ties on the first value
    /// (the vertices are distinct) so the order is the same whatever the sort.
    fn sort_partitions(&mut self) {
        self.parts.sort_by(|v1, v2| {
            if v1.len() == v2.len() && !v1.is_empty() {
                v1[0].cmp(&v2[0])
            } else {
                v2.len().cmp(&v1.len())
            }
        });
    }

    fn calc_vtx_type_counts(&mut self) {
        let g = &self.graph;
        let max_label = g
            .vtx1_labels
            .iter()
            .chain(&g.vtx2_labels)
            .copied()
            .max()
            .unwrap_or(0) as usize;
        let count = |counts: &[i32], labels: &[u32]| {
            let mut type_counts = vec![0; max_label + 1];
            for (c, &l) in counts.iter().zip(labels) {
                if *c != 0 {
                    type_counts[l as usize] += 1;
                }
            }
            type_counts
        };
        self.vtx1_type_counts = count(&self.vtx1_counts, &g.vtx1_labels);
        self.vtx2_type_counts = count(&self.vtx2_counts, &g.vtx2_labels);
    }

    fn decrement_vertex_counts(&mut self, vtx: usize) {
        let (v1, v2) = self.graph.vtx_pairs[vtx];
        self.vtx1_counts[v1 as usize] -= 1;
        if self.vtx1_counts[v1 as usize] == 0 {
            self.vtx1_type_counts[self.graph.vtx1_labels[v1 as usize] as usize] -= 1;
        }
        self.vtx2_counts[v2 as usize] -= 1;
        if self.vtx2_counts[v2 as usize] == 0 {
            self.vtx2_type_counts[self.graph.vtx2_labels[v2 as usize] as usize] -= 1;
        }
    }
}
