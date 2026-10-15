//! The topology of a seed or target (upstream `FMCS/Graph.h`): an undirected
//! `boost::adjacency_list<vecS, vecS, undirectedS, AtomIdx, BondIdx>` whose
//! vertices carry a molecule atom index and whose edges a bond index.

use crate::graphmol::substruct::Vf2Graph;

/// Sentinel for "no index" (upstream `NotSet`).
pub const NOT_SET: u32 = u32::MAX;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Graph {
    /// atom index of each vertex
    atoms: Vec<u32>,
    /// per vertex: (neighbor vertex, edge) in insertion order
    adj: Vec<Vec<(u32, u32)>>,
    /// per edge: (source vertex, target vertex, bond index)
    edges: Vec<(u32, u32, u32)>,
}

impl Graph {
    pub fn new() -> Self {
        Graph::default()
    }

    pub fn add_atom(&mut self, atom: u32) {
        self.atoms.push(atom);
        self.adj.push(Vec::new());
    }

    pub fn add_bond(&mut self, bond: u32, begin: u32, end: u32) {
        let edge = self.edges.len() as u32;
        self.edges.push((begin, end, bond));
        self.adj[begin as usize].push((end, edge));
        self.adj[end as usize].push((begin, edge));
    }

    pub fn num_vertices(&self) -> usize {
        self.atoms.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// The molecule atom index of a vertex.
    pub fn atom(&self, vertex: u32) -> u32 {
        self.atoms[vertex as usize]
    }

    /// The molecule bond index of an edge.
    pub fn bond(&self, edge: u32) -> u32 {
        self.edges[edge as usize].2
    }

    /// `(source, target)` vertices of an edge.
    pub fn edge_vertices(&self, edge: u32) -> (u32, u32) {
        let (s, t, _) = self.edges[edge as usize];
        (s, t)
    }

    /// Edges in insertion order as `(source, target, bond)`.
    pub fn edges(&self) -> &[(u32, u32, u32)] {
        &self.edges
    }

    pub fn out_degree(&self, vertex: u32) -> u32 {
        self.adj[vertex as usize].len() as u32
    }
}

impl Vf2Graph for Graph {
    fn num_nodes(&self) -> u32 {
        self.atoms.len() as u32
    }

    fn adjacency(&self, node: u32) -> &[(u32, u32)] {
        &self.adj[node as usize]
    }
}
