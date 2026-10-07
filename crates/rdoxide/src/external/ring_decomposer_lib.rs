//! Rust port of RingDecomposerLib v1.1.3_rdkit (the version RDKit 2026.09.1
//! builds against, with `qsort` replaced by a stable sort).
//!
//! Copyright (c) 2016 University of Hamburg, ZBH - Center for Bioinformatics,
//! Niek Andresen, Florian Flachsenberg, Matthias Rarey.
//! Licensed under the BSD 3-Clause license
//! (see LICENSES/RingDecomposerLib-BSD-3-Clause.txt).
//!
//! Computes Unique Ring Families (URFs), relevant cycles and an SSSR of an
//! undirected graph. The structure follows the C sources function by function
//! (`RDLgraph.c`, `RDLtarjan.c`, `RDLapsp.c`, `RDLcycleFams.c`,
//! `RDLrelation.c`, `RDLinfo.c`, `RDLhandler.c`, `RingDecomposerLib.c`) so the
//! results — including the order of SSSR rings and relevant cycles — are
//! identical.

const NONE: u32 = u32::MAX;

// ------------------------------------------------------------------ RDLgraph.c

/// `RDL_graph`. Adjacency entries are `(neighbor, edge id)`; for the
/// shortest-path DAGs (`dPaths`) only the neighbor is used.
#[derive(Debug, Clone, Default)]
pub struct Graph {
    v: usize,
    e: usize,
    adj: Vec<Vec<(u32, u32)>>,
    /// edge id -> (smaller node, larger node); empty for non-owning graphs
    edges: Vec<[u32; 2]>,
}

/// Result of [`Graph::add_u_edge`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AddEdge {
    Added(u32),
    Duplicate,
    Invalid,
}

impl Graph {
    /// `RDL_initNewGraph`
    pub fn new(v: usize) -> Self {
        Graph {
            v,
            e: 0,
            adj: vec![Vec::new(); v],
            edges: Vec::new(),
        }
    }

    pub fn num_nodes(&self) -> usize {
        self.v
    }

    pub fn num_edges(&self) -> usize {
        self.e
    }

    pub fn edge(&self, i: usize) -> [u32; 2] {
        self.edges[i]
    }

    fn degree(&self, n: usize) -> usize {
        self.adj[n].len()
    }

    /// `RDL_addEdge`: directed edge, ignored if present.
    fn add_edge(&mut self, from: u32, to: u32) {
        if from == to {
            return;
        }
        if self.adj[from as usize].iter().any(|&(n, _)| n == to) {
            return;
        }
        self.e += 1;
        self.adj[from as usize].push((to, NONE));
    }

    /// `RDL_addUEdge`
    pub fn add_u_edge(&mut self, from: u32, to: u32) -> AddEdge {
        if from as usize >= self.v || to as usize >= self.v || from == to {
            return AddEdge::Invalid;
        }
        if self.adj[from as usize].iter().any(|&(n, _)| n == to) {
            return AddEdge::Duplicate;
        }
        self.add_edge(from, to);
        self.add_edge(to, from);
        self.e -= 1;
        let id = self.e as u32 - 1;
        let (a, b) = if from > to { (to, from) } else { (from, to) };
        self.edges.push([a, b]);
        self.adj[from as usize].last_mut().expect("just added").1 = id;
        self.adj[to as usize].last_mut().expect("just added").1 = id;
        AddEdge::Added(id)
    }

    /// `RDL_edgeId`
    fn edge_id(&self, from: u32, to: u32) -> u32 {
        let (a, b) = if from > to { (to, from) } else { (from, to) };
        self.adj[a as usize]
            .iter()
            .find(|&&(n, _)| n == b)
            .map_or(NONE, |&(_, e)| e)
    }
}

// ------------------------------------------------------------------ RDLtarjan.c

#[derive(Debug)]
struct BccGraph {
    bcc_graphs: Vec<Graph>,
    /// per node: list of (bcc, local node index); used by the ring-system API
    #[allow(dead_code)]
    node_to_bcc: Vec<Vec<(u32, u32)>>,
    edge_from_bcc: Vec<Vec<u32>>,
    node_from_bcc: Vec<Vec<u32>>,
}

#[derive(Debug, Clone, Copy)]
struct DfsElem {
    u: u32,
    parent: u32,
    j: usize,
    expected_time: u32,
}

fn tarjan_visit(
    g: &Graph,
    start: u32,
    d: &mut [u32],
    low: &mut [u32],
    time: &mut u32,
    curr_bcc: &mut u32,
    bcc: &mut [u32],
    edge_stack: &mut Vec<u32>,
) {
    let mut dfs: Vec<DfsElem> = Vec::new();
    *time += 1;
    d[start as usize] = *time;
    low[start as usize] = d[start as usize];
    dfs.push(DfsElem {
        u: start,
        parent: NONE,
        j: 0,
        expected_time: NONE,
    });
    while let Some(top) = dfs.last().copied() {
        let u = top.u as usize;
        if top.j < g.degree(u) {
            let v = g.adj[u][top.j].0;
            let uv_edge = g.edge_id(top.u, v);
            if d[v as usize] == 0 {
                edge_stack.push(uv_edge);
                *time += 1;
                d[v as usize] = *time;
                low[v as usize] = d[v as usize];
                dfs.last_mut().expect("top").expected_time = *time;
                dfs.push(DfsElem {
                    u: v,
                    parent: top.u,
                    j: 0,
                    expected_time: NONE,
                });
            } else {
                dfs.last_mut().expect("top").j += 1;
                if d[v as usize] == top.expected_time {
                    low[u] = low[u].min(low[v as usize]);
                    if low[v as usize] >= d[u] {
                        loop {
                            let curr_edge = edge_stack.pop().expect("edge stack");
                            bcc[curr_edge as usize] = *curr_bcc;
                            if curr_edge == uv_edge {
                                break;
                            }
                        }
                        *curr_bcc += 1;
                    }
                } else if d[v as usize] < d[u] && v != top.parent {
                    edge_stack.push(uv_edge);
                    low[u] = low[u].min(d[v as usize]);
                }
            }
        } else {
            dfs.pop();
        }
    }
}

/// `RDL_tarjanBCC`
fn tarjan_bcc(g: &Graph) -> BccGraph {
    let mut bcc = vec![0u32; g.e];
    let mut d = vec![0u32; g.v];
    let mut low = vec![0u32; g.v];
    let mut edge_stack = Vec::new();
    let mut time = 0;
    let mut curr_bcc = 1;
    for u in 0..g.v {
        if d[u] == 0 {
            tarjan_visit(
                g,
                u as u32,
                &mut d,
                &mut low,
                &mut time,
                &mut curr_bcc,
                &mut bcc,
                &mut edge_stack,
            );
        }
    }
    let nof_bcc = (curr_bcc - 1) as usize;
    let mut bcc_size = vec![0u32; nof_bcc];
    for &b in &bcc {
        bcc_size[b as usize - 1] += 1;
    }
    let mut non_trivial = vec![NONE; nof_bcc];
    let mut n_nontrivial = 0u32;
    for i in 0..nof_bcc {
        if bcc_size[i] > 1 {
            non_trivial[i] = n_nontrivial;
            n_nontrivial += 1;
        }
    }
    let nb = n_nontrivial as usize;
    let mut edge_to_bcc = vec![(NONE, NONE); g.e];
    let mut edge_from_bcc: Vec<Vec<u32>> = vec![Vec::new(); nb];
    let mut node_to_bcc: Vec<Vec<(u32, u32)>> = vec![Vec::new(); g.v];
    let mut node_from_bcc: Vec<Vec<u32>> = vec![Vec::new(); nb];
    for i in 0..g.e {
        let cb = bcc[i] as usize;
        if bcc_size[cb - 1] <= 1 {
            continue;
        }
        let u = non_trivial[cb - 1];
        edge_to_bcc[i] = (u, edge_from_bcc[u as usize].len() as u32);
        edge_from_bcc[u as usize].push(i as u32);
        for k in 0..2 {
            let node = g.edges[i][k] as usize;
            if !node_to_bcc[node].iter().any(|&(b, _)| b == u) {
                let local = node_from_bcc[u as usize].len() as u32;
                node_to_bcc[node].push((u, local));
                node_from_bcc[u as usize].push(node as u32);
            }
        }
    }
    let mut bcc_graphs: Vec<Graph> = node_from_bcc.iter().map(|n| Graph::new(n.len())).collect();
    for i in 0..g.e {
        let u = edge_to_bcc[i].0;
        if u == NONE {
            continue;
        }
        let mut ii = [NONE; 2];
        for k in 0..2 {
            let node = g.edges[i][k] as usize;
            for &(b, local) in &node_to_bcc[node] {
                if b == u {
                    ii[k] = local;
                }
            }
        }
        bcc_graphs[u as usize].add_u_edge(ii[0], ii[1]);
    }
    BccGraph {
        bcc_graphs,
        node_to_bcc,
        edge_from_bcc,
        node_from_bcc,
    }
}

// ------------------------------------------------------------------ RDLapsp.c

#[derive(Debug)]
struct SPathInfo {
    pred: Vec<Vec<u32>>,
    dist: Vec<Vec<u32>>,
    reachable: Vec<Vec<bool>>,
    d_paths: Vec<Graph>,
}

/// `RDL_AllPairsShortestPaths`
fn all_pairs_shortest_paths(g: &Graph) -> SPathInfo {
    let n = g.v;
    let mut spi = SPathInfo {
        pred: vec![vec![NONE; n]; n],
        dist: vec![vec![NONE; n]; n],
        reachable: vec![vec![false; n]; n],
        d_paths: (0..n).map(|_| Graph::new(n)).collect(),
    };
    let mut color = vec![false; n]; // true = black
    let mut queue = vec![0u32; n];
    for run in 1..3 {
        for i in 0..n {
            let (mut head, mut next, mut size) = (0usize, 0usize, 0usize);
            for c in color.iter_mut() {
                *c = false;
            }
            spi.dist[i][i] = 0;
            spi.pred[i][i] = i as u32;
            color[i] = true;
            queue[next] = i as u32;
            next += 1;
            size += 1;
            while size > 0 {
                let w = queue[head] as usize;
                head += 1;
                size -= 1;
                for &(j, _) in &g.adj[w] {
                    let j = j as usize;
                    let allowed = run == 2
                        || (g.degree(j) == g.degree(i) && j < i)
                        || g.degree(j) < g.degree(i);
                    if allowed && !color[j] {
                        let dw1 = spi.dist[i][w].saturating_add(1);
                        if run == 2 && dw1 < spi.dist[i][j] {
                            spi.reachable[i][j] = false;
                        }
                        if run == 1 || dw1 < spi.dist[i][j] {
                            spi.pred[i][j] = w as u32;
                        }
                        if spi.dist[i][j] > dw1 {
                            spi.dist[i][j] = dw1;
                        }
                        color[j] = true;
                        queue[next] = j as u32;
                        next += 1;
                        size += 1;
                        if run == 1 {
                            spi.reachable[i][j] = true;
                        }
                    }
                }
            }
        }
    }
    spi
}

// ------------------------------------------------------------------ RDLcycleFams.c

/// A cycle family (`RDL_cfam`).
#[derive(Debug, Clone)]
struct CFam {
    weight: u32,
    r: u32,
    p: u32,
    q: u32,
    x: u32,
    prototype: Vec<bool>,
    mark: bool,
}

fn paths_share_only_start(r: u32, y: u32, z: u32, g: &Graph, spi: &SPathInfo) -> bool {
    let mut in_ry = vec![false; g.v];
    let mut in_rz = vec![false; g.v];
    in_ry[y as usize] = true;
    in_rz[z as usize] = true;
    let mut pnt = y;
    loop {
        pnt = spi.pred[r as usize][pnt as usize];
        in_ry[pnt as usize] = true;
        if pnt == r {
            break;
        }
    }
    pnt = z;
    loop {
        pnt = spi.pred[r as usize][pnt as usize];
        in_rz[pnt as usize] = true;
        if pnt == r {
            break;
        }
    }
    let mut count = 0;
    for i in 0..g.v {
        if count >= 2 {
            break;
        }
        if in_ry[i] && in_rz[i] {
            count += 1;
        }
    }
    count == 1 && in_ry[r as usize] && in_rz[r as usize]
}

fn find_prototype(r: u32, y: u32, z: u32, x: u32, g: &Graph, spi: &SPathInfo) -> Vec<bool> {
    let mut proto = vec![false; g.e];
    for start in [y, z] {
        let mut v1 = start;
        loop {
            let v2 = v1;
            v1 = spi.pred[r as usize][v1 as usize];
            proto[g.edge_id(v1, v2) as usize] = true;
            if v1 == r {
                break;
            }
        }
    }
    if x == NONE {
        proto[g.edge_id(y, z) as usize] = true;
    } else {
        proto[g.edge_id(y, x) as usize] = true;
        proto[g.edge_id(z, x) as usize] = true;
    }
    proto
}

/// Vismara's algorithm (`RDL_vismara`).
fn vismara(fams: &mut Vec<CFam>, g: &Graph, spi: &mut SPathInfo) {
    let mut even_cand: Vec<u32> = Vec::with_capacity(g.v);
    for rv in 0..g.v as u32 {
        for yv in 0..g.v as u32 {
            if !spi.reachable[rv as usize][yv as usize] {
                continue;
            }
            even_cand.clear();
            for i in 0..g.degree(yv as usize) {
                let zv = g.adj[yv as usize][i].0;
                if spi.reachable[rv as usize][zv as usize] {
                    let dz = spi.dist[rv as usize][zv as usize];
                    let dy = spi.dist[rv as usize][yv as usize];
                    if dz.wrapping_add(1) == dy {
                        even_cand.push(zv);
                    } else if dz != dy.wrapping_add(1)
                        && (g.degree(zv as usize) < g.degree(yv as usize)
                            || (g.degree(zv as usize) == g.degree(yv as usize) && zv < yv))
                        && paths_share_only_start(rv, yv, zv, g, spi)
                    {
                        fams.push(CFam {
                            r: rv,
                            p: yv,
                            q: zv,
                            x: NONE,
                            mark: false,
                            prototype: find_prototype(rv, yv, zv, NONE, g, spi),
                            weight: spi.dist[rv as usize][yv as usize]
                                + spi.dist[rv as usize][zv as usize]
                                + 1,
                        });
                    }
                    if spi.dist[rv as usize][zv as usize] == 1 {
                        spi.d_paths[rv as usize].add_edge(zv, rv);
                    }
                }
            }
            for i in 0..even_cand.len() {
                let pv = even_cand[i];
                for &qv in &even_cand[i + 1..] {
                    if paths_share_only_start(rv, pv, qv, g, spi) {
                        fams.push(CFam {
                            r: rv,
                            p: pv,
                            q: qv,
                            x: yv,
                            mark: false,
                            prototype: find_prototype(rv, pv, qv, yv, g, spi),
                            weight: spi.dist[rv as usize][pv as usize]
                                + spi.dist[rv as usize][qv as usize]
                                + 2,
                        });
                    }
                }
            }
            for &pv in &even_cand {
                spi.d_paths[rv as usize].add_edge(yv, pv);
            }
        }
    }
}

/// `RDL_findCycleFams`: families sorted (stably) by weight.
fn find_cycle_fams(g: &Graph, spi: &mut SPathInfo) -> Vec<CFam> {
    let mut fams = Vec::new();
    vismara(&mut fams, g, spi);
    fams.sort_by_key(|f| f.weight);
    fams
}

// ------------------------------------------------------------------ RDLinfo.c / RDLrelation.c

#[derive(Debug)]
struct UrfInfo {
    nof_protos: Vec<usize>,
    urf_rel: Vec<Vec<Vec<bool>>>,
    /// URF -> indices into the family list
    urfs: Vec<Vec<usize>>,
}

impl UrfInfo {
    fn idx_weights(&self, weight: usize, j: usize) -> usize {
        self.nof_protos[..weight].iter().sum::<usize>() + j
    }
}

fn init_urf_info(fams: &[CFam]) -> UrfInfo {
    let mut nof_protos: Vec<usize> = Vec::new();
    let mut curr = fams[0].weight;
    nof_protos.push(0);
    for f in fams {
        if f.weight != curr {
            nof_protos.push(0);
            curr = f.weight;
        }
        *nof_protos.last_mut().expect("non-empty") += 1;
    }
    let urf_rel = nof_protos
        .iter()
        .map(|&n| vec![vec![false; n]; n])
        .collect();
    UrfInfo {
        nof_protos,
        urf_rel,
        urfs: Vec::new(),
    }
}

/// Column swap on a set of bit rows (`RDL_swap_columns`).
fn swap_columns(rows: &mut [Vec<bool>], c1: usize, c2: usize) {
    for r in rows.iter_mut() {
        r.swap(c1, c2);
    }
}

fn xor_inplace(dst: &mut [bool], src: &[bool]) {
    for (d, s) in dst.iter_mut().zip(src) {
        *d ^= *s;
    }
}

/// `RDL_checkDependencies`: marks relevant families and potential
/// URF-relations via Gaussian elimination over GF(2).
fn check_dependencies(fams: &mut [CFam], g: &Graph, u: &mut UrfInfo) {
    if fams.len() < 3 {
        for i in 0..u.nof_protos.len() {
            for j in 0..u.nof_protos[i] {
                u.urf_rel[i][j][j] = true;
            }
        }
        for f in fams.iter_mut() {
            f.mark = true;
        }
        return;
    }
    let mut basis: Vec<Vec<bool>> = Vec::new();
    let mut relevant: Vec<Vec<bool>> = Vec::new();
    let mut relevant_map: Vec<usize> = Vec::new();
    let mut protos: Vec<Vec<bool>> = fams.iter().map(|f| f.prototype.clone()).collect();
    for i in 0..u.nof_protos.len() {
        let basis_smaller_end = basis.len();
        let relevant_smaller_end = relevant.len();
        for j in 0..u.nof_protos[i] {
            let index = u.idx_weights(i, j);
            let mut with_smaller = protos[index].clone();
            for k in 0..basis_smaller_end {
                if with_smaller[k] {
                    xor_inplace(&mut with_smaller, &basis[k]);
                }
            }
            if !with_smaller.iter().any(|&b| b) {
                continue;
            }
            relevant.push(with_smaller.clone());
            relevant_map.push(j);
            fams[index].mark = true;
            u.urf_rel[i][j][j] = true;
            let mut with_equal = with_smaller.clone();
            for k in basis_smaller_end..basis.len() {
                if with_equal[k] {
                    xor_inplace(&mut with_equal, &basis[k]);
                }
            }
            if !with_equal.iter().any(|&b| b) {
                for k in relevant_smaller_end..relevant.len() - 1 {
                    let mut t = with_smaller.clone();
                    xor_inplace(&mut t, &relevant[k]);
                    if !t.iter().any(|&b| b) {
                        let m = relevant_map[k];
                        u.urf_rel[i][j][m] = true;
                        u.urf_rel[i][m][j] = true;
                    }
                }
            } else {
                let old = basis.len();
                let has_pivot = with_equal[old];
                basis.push(with_equal);
                if !has_pivot {
                    let last = basis.len() - 1;
                    if let Some(k) = (old + 1..g.e).find(|&k| basis[last][k]) {
                        swap_columns(&mut basis, old, k);
                        swap_columns(&mut relevant, old, k);
                        swap_columns(&mut protos, old, k);
                    }
                }
            }
        }
    }
}

/// `RDL_giveEdges`: mark all edges on shortest paths from `b` back to `a`.
fn give_edges(a: u32, b: u32, out: &mut [bool], g: &Graph, spi: &SPathInfo, visited: &mut [bool]) {
    if a == b {
        return;
    }
    visited[b as usize] = true;
    let nbrs: Vec<u32> = spi.d_paths[a as usize].adj[b as usize]
        .iter()
        .map(|&(n, _)| n)
        .collect();
    for vertex in nbrs {
        out[g.edge_id(b, vertex) as usize] = true;
        if !visited[vertex as usize] {
            give_edges(a, vertex, out, g, spi, visited);
        }
    }
}

/// `RDL_giveVertices`
fn give_vertices(a: u32, b: u32, out: &mut [bool], spi: &SPathInfo, visited: &mut [bool]) {
    visited[b as usize] = true;
    if a == b {
        out[a as usize] = true;
        return;
    }
    out[b as usize] = true;
    let nbrs: Vec<u32> = spi.d_paths[a as usize].adj[b as usize]
        .iter()
        .map(|&(n, _)| n)
        .collect();
    for vertex in nbrs {
        if !visited[vertex as usize] {
            give_vertices(a, vertex, out, spi, visited);
        }
    }
}

/// `RDL_findEdges`
fn find_edges(edges: &mut [bool], f: &CFam, g: &Graph, spi: &SPathInfo) {
    let mut visited = vec![false; g.e.max(g.v)];
    give_edges(f.r, f.p, edges, g, spi, &mut visited);
    visited.iter_mut().for_each(|v| *v = false);
    give_edges(f.r, f.q, edges, g, spi, &mut visited);
    if f.x < NONE {
        edges[g.edge_id(f.x, f.p) as usize] = true;
        edges[g.edge_id(f.x, f.q) as usize] = true;
    } else {
        edges[g.edge_id(f.p, f.q) as usize] = true;
    }
}

fn make_edge_list(g: &Graph, f: &CFam, spi: &SPathInfo) -> Vec<u32> {
    let mut edges = vec![false; g.e];
    find_edges(&mut edges, f, g, spi);
    edges
        .iter()
        .enumerate()
        .filter(|&(_, &b)| b)
        .map(|(i, _)| i as u32)
        .collect()
}

/// `RDL_checkEdges`: potential URF relations need a shared edge.
fn check_edges(fams: &[CFam], g: &Graph, u: &mut UrfInfo, spi: &SPathInfo) {
    for i in 0..u.nof_protos.len() {
        let n = u.nof_protos[i];
        let mut lists: Vec<Option<Vec<u32>>> = vec![None; n];
        for j in 0..n {
            for k in j + 1..n {
                if !u.urf_rel[i][j][k] {
                    continue;
                }
                if lists[j].is_none() {
                    lists[j] = Some(make_edge_list(g, &fams[u.idx_weights(i, j)], spi));
                }
                if lists[k].is_none() {
                    lists[k] = Some(make_edge_list(g, &fams[u.idx_weights(i, k)], spi));
                }
                let (l1, l2) = (
                    lists[j].as_ref().expect("set"),
                    lists[k].as_ref().expect("set"),
                );
                let (mut a, mut b) = (0, 0);
                let mut shared = false;
                while a < l1.len() && b < l2.len() && !shared {
                    if l1[a] < l2[b] {
                        a += 1;
                    } else if l1[a] > l2[b] {
                        b += 1;
                    } else {
                        shared = true;
                    }
                }
                u.urf_rel[i][j][k] = shared;
                u.urf_rel[i][k][j] = shared;
            }
        }
    }
}

/// `RDL_findTransitiveClosure`
fn find_transitive_closure(u: &mut UrfInfo) {
    for i in 0..u.nof_protos.len() {
        let n = u.nof_protos[i];
        let mut visited = vec![false; n];
        let mut ccs: Vec<Vec<usize>> = Vec::new();
        let mut stack: Vec<usize> = Vec::new();
        for j in 0..n {
            if visited[j] {
                continue;
            }
            let mut cc = Vec::new();
            stack.push(j);
            visited[j] = true;
            while let Some(k) = stack.pop() {
                cc.push(k);
                for l in 0..n {
                    if !visited[l] && u.urf_rel[i][k][l] {
                        visited[l] = true;
                        stack.push(l);
                    }
                }
            }
            ccs.push(cc);
        }
        for cc in ccs {
            for a in 0..cc.len() {
                for b in a + 1..cc.len() {
                    u.urf_rel[i][cc[a]][cc[b]] = true;
                    u.urf_rel[i][cc[b]][cc[a]] = true;
                }
            }
        }
    }
}

/// `RDL_fillURFs`
fn fill_urfs(u: &mut UrfInfo) {
    let mut already: Vec<Vec<bool>> = u.nof_protos.iter().map(|&n| vec![false; n]).collect();
    let mut urfs = Vec::new();
    for i in 0..u.nof_protos.len() {
        for j in 0..u.nof_protos[i] {
            if u.urf_rel[i][j][j] && !already[i][j] {
                let mut urf = Vec::new();
                for k in j..u.nof_protos[i] {
                    if u.urf_rel[i][j][k] {
                        urf.push(u.idx_weights(i, k));
                        already[i][k] = true;
                    }
                }
                urfs.push(urf);
            }
        }
    }
    u.urfs = urfs;
}

/// `RDL_checkURFRelation`
fn check_urf_relation(fams: &mut [CFam], g: &Graph, spi: &SPathInfo) -> UrfInfo {
    let mut u = init_urf_info(fams);
    check_dependencies(fams, g, &mut u);
    check_edges(fams, g, &mut u, spi);
    find_transitive_closure(&mut u);
    fill_urfs(&mut u);
    u
}

// ------------------------------------------------------------------ RingDecomposerLib.c

#[derive(Debug)]
struct BccData {
    spi: SPathInfo,
    fams: Vec<CFam>,
    urf_info: Option<UrfInfo>,
}

/// `RDL_data`: the result of [`calculate`].
#[derive(Debug)]
pub struct Data {
    graph: Graph,
    bcc: BccGraph,
    per_bcc: Vec<BccData>,
    urf_to_bcc: Vec<(usize, usize)>,
}

/// An edge-defined cycle (`RDL_cycle`): node pairs in the original graph.
pub type Cycle = Vec<[u32; 2]>;

/// `RDL_calculate`. Returns `None` for a graph without nodes.
pub fn calculate(graph: Graph) -> Option<Data> {
    if graph.v == 0 {
        return None;
    }
    let bcc = tarjan_bcc(&graph);
    let mut per_bcc = Vec::with_capacity(bcc.bcc_graphs.len());
    let mut urf_to_bcc = Vec::new();
    for (i, g) in bcc.bcc_graphs.iter().enumerate() {
        let mut spi = all_pairs_shortest_paths(g);
        let mut fams = find_cycle_fams(g, &mut spi);
        let urf_info = if fams.is_empty() {
            None
        } else {
            Some(check_urf_relation(&mut fams, g, &spi))
        };
        if let Some(u) = &urf_info {
            for j in 0..u.urfs.len() {
                urf_to_bcc.push((i, j));
            }
        }
        per_bcc.push(BccData {
            spi,
            fams,
            urf_info,
        });
    }
    Some(Data {
        graph,
        bcc,
        per_bcc,
        urf_to_bcc,
    })
}

impl Data {
    fn urf_fams(&self, bcc: usize, urf: usize) -> &[usize] {
        &self.per_bcc[bcc].urf_info.as_ref().expect("URF info").urfs[urf]
    }

    /// `RDL_getNofURF`
    pub fn num_urfs(&self) -> usize {
        self.urf_to_bcc.len()
    }

    /// `RDL_getNodesForURF`: original node indices (in BCC-local order).
    pub fn nodes_for_urf(&self, index: usize) -> Vec<u32> {
        let (b, u) = self.urf_to_bcc[index];
        let g = &self.bcc.bcc_graphs[b];
        let spi = &self.per_bcc[b].spi;
        let mut atoms = vec![false; g.v];
        for &fi in self.urf_fams(b, u) {
            let f = &self.per_bcc[b].fams[fi];
            let mut visited = vec![false; g.v];
            give_vertices(f.r, f.q, &mut atoms, spi, &mut visited);
            let mut visited = vec![false; g.v];
            give_vertices(f.r, f.p, &mut atoms, spi, &mut visited);
            if f.x < NONE {
                atoms[f.x as usize] = true;
            }
        }
        atoms
            .iter()
            .enumerate()
            .filter(|&(_, &a)| a)
            .map(|(i, _)| self.bcc.node_from_bcc[b][i])
            .collect()
    }

    /// `RDL_getEdgesForURF`: edges as node pairs of the original graph.
    pub fn edges_for_urf(&self, index: usize) -> Vec<[u32; 2]> {
        let (b, u) = self.urf_to_bcc[index];
        let g = &self.bcc.bcc_graphs[b];
        let spi = &self.per_bcc[b].spi;
        let mut edges = vec![false; g.e];
        for &fi in self.urf_fams(b, u) {
            let f = &self.per_bcc[b].fams[fi];
            let mut visited = vec![false; g.v];
            give_edges(f.r, f.q, &mut edges, g, spi, &mut visited);
            let mut visited = vec![false; g.v];
            give_edges(f.r, f.p, &mut edges, g, spi, &mut visited);
            if f.x < NONE {
                edges[g.edge_id(f.q, f.x) as usize] = true;
                edges[g.edge_id(f.p, f.x) as usize] = true;
            } else {
                edges[g.edge_id(f.q, f.p) as usize] = true;
            }
        }
        edges
            .iter()
            .enumerate()
            .filter(|&(_, &e)| e)
            .map(|(i, _)| self.graph.edges[self.bcc.edge_from_bcc[b][i] as usize])
            .collect()
    }

    /// `RDL_getSSSR`
    pub fn sssr(&self) -> Vec<Cycle> {
        let mut result = Vec::new();
        for (b, g) in self.bcc.bcc_graphs.iter().enumerate() {
            let Some(u) = &self.per_bcc[b].urf_info else {
                continue;
            };
            let size = g.e + 1 - g.v;
            let mut protos: Vec<Vec<bool>> = u
                .urfs
                .iter()
                .map(|urf| self.per_bcc[b].fams[urf[0]].prototype.clone())
                .collect();
            let mut basis: Vec<Vec<bool>> = Vec::new();
            let mut added = 0;
            for internal in 0..u.urfs.len() {
                let mut cur = protos[internal].clone();
                for k in 0..basis.len() {
                    if cur[k] {
                        xor_inplace(&mut cur, &basis[k]);
                    }
                }
                if !cur.iter().any(|&x| x) {
                    continue;
                }
                let old = basis.len();
                let has_pivot = cur[old];
                basis.push(cur);
                if !has_pivot {
                    let last = basis.len() - 1;
                    if let Some(k) = (old + 1..g.e).find(|&k| basis[last][k]) {
                        swap_columns(&mut basis, old, k);
                        swap_columns(&mut protos, old, k);
                    }
                }
                let proto = &self.per_bcc[b].fams[u.urfs[internal][0]].prototype;
                let cycle: Cycle = proto
                    .iter()
                    .enumerate()
                    .filter(|&(_, &x)| x)
                    .map(|(j, _)| self.graph.edges[self.bcc.edge_from_bcc[b][j] as usize])
                    .collect();
                result.push(cycle);
                added += 1;
                if added == size {
                    break;
                }
            }
        }
        result
    }

    /// `RDL_getNofRC`: number of relevant cycles (as a float, like upstream).
    pub fn num_relevant_cycles(&self) -> f64 {
        let mut total = 0.0;
        for &(b, u) in &self.urf_to_bcc {
            for &fi in self.urf_fams(b, u) {
                let f = &self.per_bcc[b].fams[fi];
                let g = &self.bcc.bcc_graphs[b];
                let spi = &self.per_bcc[b].spi;
                total += count_paths(f.r, f.q, g.v, spi) * count_paths(f.r, f.p, g.v, spi);
            }
        }
        total
    }

    /// All relevant cycles, in `RDL_getRCyclesIterator` order.
    pub fn relevant_cycles(&self) -> Vec<Cycle> {
        let mut out = Vec::new();
        for (b, g) in self.bcc.bcc_graphs.iter().enumerate() {
            let Some(u) = &self.per_bcc[b].urf_info else {
                continue;
            };
            let spi = &self.per_bcc[b].spi;
            for urf in &u.urfs {
                for &fi in urf {
                    let f = &self.per_bcc[b].fams[fi];
                    let paths_p = list_paths(f.r, f.p, g, spi);
                    let paths_q = list_paths(f.r, f.q, g, spi);
                    // the C iterator advances the q-side fastest
                    for pp in &paths_p {
                        for pq in &paths_q {
                            let mut cyc: Vec<bool> =
                                pp.iter().zip(pq).map(|(a, b)| *a || *b).collect();
                            if f.x < NONE {
                                cyc[g.edge_id(f.p, f.x) as usize] = true;
                                cyc[g.edge_id(f.q, f.x) as usize] = true;
                            } else {
                                cyc[g.edge_id(f.p, f.q) as usize] = true;
                            }
                            out.push(
                                cyc.iter()
                                    .enumerate()
                                    .filter(|&(_, &x)| x)
                                    .map(|(j, _)| {
                                        self.graph.edges[self.bcc.edge_from_bcc[b][j] as usize]
                                    })
                                    .collect(),
                            );
                        }
                    }
                }
            }
        }
        out
    }
}

/// `RDL_countPaths`: number of shortest paths from `b` back to `a`.
fn count_paths(a: u32, b: u32, v: usize, spi: &SPathInfo) -> f64 {
    let mut indegree = vec![NONE; v];
    fn find_indegree(a: u32, b: u32, indegree: &mut [u32], spi: &SPathInfo) {
        indegree[b as usize] = 0;
        if a == b {
            return;
        }
        for &(vertex, _) in &spi.d_paths[a as usize].adj[b as usize] {
            if indegree[vertex as usize] == NONE {
                find_indegree(a, vertex, indegree, spi);
            }
            indegree[vertex as usize] += 1;
        }
    }
    find_indegree(a, b, &mut indegree, spi);
    let mut nof_paths = vec![0.0f64; v];
    let mut empty: Vec<usize> = Vec::new();
    for i in 0..v {
        if indegree[i] == 0 {
            empty.push(i);
            nof_paths[i] = 1.0;
        }
    }
    while let Some(curr) = empty.pop() {
        for &(vertex, _) in &spi.d_paths[a as usize].adj[curr] {
            let vtx = vertex as usize;
            indegree[vtx] -= 1;
            nof_paths[vtx] += nof_paths[curr];
            if indegree[vtx] == 0 {
                empty.push(vtx);
            }
        }
    }
    nof_paths[a as usize]
}

/// All shortest paths from `p` back to `r` as edge bitsets, in the order the
/// C `RDL_pathIterator` (mode 'b') produces them.
fn list_paths(r: u32, p: u32, g: &Graph, spi: &SPathInfo) -> Vec<Vec<bool>> {
    struct Elem {
        p: u32,
        ai_index: usize,
        split_copy: Option<Vec<bool>>,
    }
    let mut paths = Vec::new();
    let mut path = vec![false; g.e];
    let mut stack = vec![Elem {
        p,
        ai_index: 0,
        split_copy: None,
    }];
    loop {
        // one call of RDL_pathIteratorNext
        if stack.is_empty() {
            break;
        }
        while let Some(top) = stack.last_mut() {
            if top.p == r {
                stack.pop();
                break;
            }
            let adj = &spi.d_paths[r as usize].adj[top.p as usize];
            if adj.len() > 1 {
                let vertex = adj[top.ai_index].0;
                if top.ai_index == 0 {
                    top.split_copy = Some(path.clone());
                } else {
                    path.clone_from(top.split_copy.as_ref().expect("split copy"));
                }
                path[g.edge_id(top.p, vertex) as usize] = true;
                top.ai_index += 1;
                if top.ai_index >= adj.len() {
                    stack.pop();
                }
                stack.push(Elem {
                    p: vertex,
                    ai_index: 0,
                    split_copy: None,
                });
            } else {
                let vertex = adj[0].0;
                path[g.edge_id(top.p, vertex) as usize] = true;
                top.p = vertex;
                top.ai_index = 0;
            }
        }
        paths.push(path.clone());
    }
    paths
}

/// `RingUtils::rdlCycleToAtomRing` is in `graphmol::find_rings`.
#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, edges: &[(u32, u32)]) -> Graph {
        let mut g = Graph::new(n);
        for &(a, b) in edges {
            g.add_u_edge(a, b);
        }
        g
    }

    #[test]
    fn benzene() {
        let g = graph(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]);
        let d = calculate(g).unwrap();
        assert_eq!(d.num_urfs(), 1);
        assert_eq!(d.sssr().len(), 1);
        assert_eq!(d.sssr()[0].len(), 6);
        assert_eq!(d.num_relevant_cycles(), 1.0);
        assert_eq!(d.nodes_for_urf(0).len(), 6);
    }

    #[test]
    fn naphthalene_and_chain() {
        // two fused 6-rings plus a pendant chain (no BCC for the chain)
        let g = graph(
            12,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 0),
                (4, 6),
                (6, 7),
                (7, 8),
                (8, 9),
                (9, 5),
                (9, 10),
                (10, 11),
            ],
        );
        let d = calculate(g).unwrap();
        assert_eq!(d.num_urfs(), 2);
        let s = d.sssr();
        assert_eq!(s.len(), 2);
        assert!(s.iter().all(|c| c.len() == 6));
    }

    #[test]
    fn cubane_relevant_cycles() {
        // cube graph: 6 four-membered rings, SSSR has 5
        let e = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let d = calculate(graph(8, &e)).unwrap();
        assert_eq!(d.sssr().len(), 5);
        assert_eq!(d.relevant_cycles().len(), 6);
        assert_eq!(d.num_relevant_cycles(), 6.0);
    }

    #[test]
    fn no_rings() {
        let d = calculate(graph(3, &[(0, 1), (1, 2)])).unwrap();
        assert_eq!(d.num_urfs(), 0);
        assert!(d.sssr().is_empty());
        assert!(calculate(Graph::new(0)).is_none());
    }
}
