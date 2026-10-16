//! The hierarchical digraph of an atom (upstream `CIPLabeler/{Digraph,Node,Edge}.{h,cpp}`).
//!
//! Nodes and edges live in arenas owned by the [`Digraph`] and refer to each
//! other by index. Nodes are expanded lazily the first time their edges are
//! asked for, hence the accessors that take `&mut self`.

use super::cip_mol::CipMol;
use super::descriptor::Descriptor;
use super::mancude::Rational;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::{RDKitError, Result};

pub type NodeId = usize;
pub type EdgeId = usize;

/// `Atom::NOATOM`, the index reported for implicit hydrogens
pub const NOATOM: u32 = u32::MAX;

const MAX_NODE_COUNT: usize = 100_000;

pub const EXPANDED: u8 = 0x1;
pub const RING_DUPLICATE: u8 = 0x2;
pub const BOND_DUPLICATE: u8 = 0x4;
pub const DUPLICATE: u8 = RING_DUPLICATE | BOND_DUPLICATE;
pub const IMPL_HYDROGEN: u8 = 0x8;
pub const DUPLICATE_OR_H: u8 = RING_DUPLICATE | BOND_DUPLICATE | IMPL_HYDROGEN;

#[derive(Debug)]
pub struct Node {
    /// `None` for implicit hydrogens
    pub atom: Option<u32>,
    pub dist: i32,
    pub atomic_num: Rational,
    pub atomic_mass: f64,
    pub aux: Descriptor,
    flags: u8,
    edges: Vec<EdgeId>,
    visit: Vec<u32>,
}

impl Node {
    fn new(
        mol: &CipMol,
        visit: Vec<u32>,
        atom: Option<u32>,
        frac: Rational,
        dist: i32,
        flags: u8,
    ) -> Self {
        let mut node = Node {
            atom,
            dist,
            atomic_num: frac,
            atomic_mass: 0.0,
            aux: Descriptor::None,
            flags,
            edges: Vec::new(),
            visit,
        };
        if node.flags & DUPLICATE == 0 {
            let table = PeriodicTable::get();
            let atomic_number = node.atomic_number(mol);
            let isotope = node.mass_num(mol);
            if isotope == 0 {
                node.atomic_mass = table.get_atomic_weight(atomic_number).unwrap_or(0.0);
            } else {
                node.atomic_mass = table
                    .get_mass_for_isotope(atomic_number, isotope)
                    .unwrap_or(0.0);
                if atomic_number != 0 && node.atomic_mass == 0.0 {
                    node.atomic_mass = f64::from(isotope);
                }
            }
        }
        if node.visit.is_empty() || node.flags & DUPLICATE != 0 {
            node.flags |= EXPANDED;
        }
        node
    }

    /// `getAtomicNum`: implicit hydrogens are hydrogens
    pub fn atomic_number(&self, mol: &CipMol) -> u32 {
        self.atom.map_or(1, |a| mol.atom(a).atomic_num())
    }

    /// `getMassNum`
    pub fn mass_num(&self, mol: &CipMol) -> u32 {
        match self.atom {
            Some(a) if !self.is_duplicate() => mol.atom(a).isotope(),
            _ => 0,
        }
    }

    /// `getAtomIdx`
    pub fn atom_idx(&self) -> u32 {
        if self.is_set(IMPL_HYDROGEN) {
            NOATOM
        } else {
            self.atom.unwrap_or(NOATOM)
        }
    }

    pub fn is_set(&self, mask: u8) -> bool {
        mask & self.flags != 0
    }

    pub fn is_duplicate(&self) -> bool {
        self.flags & DUPLICATE != 0
    }

    pub fn is_duplicate_or_h(&self) -> bool {
        self.flags & DUPLICATE_OR_H != 0
    }

    pub fn is_expanded(&self) -> bool {
        self.flags & EXPANDED != 0
    }

    pub fn is_terminal(&self) -> bool {
        self.visit.is_empty() || (self.is_expanded() && self.edges.len() == 1)
    }

    fn is_visited(&self, idx: u32) -> bool {
        self.visit[idx as usize] != 0
    }
}

#[derive(Debug)]
pub struct Edge {
    pub beg: NodeId,
    pub end: NodeId,
    pub bond: Option<u32>,
    pub aux: Descriptor,
}

impl Edge {
    pub fn is_beg(&self, node: NodeId) -> bool {
        self.beg == node
    }

    pub fn is_end(&self, node: NodeId) -> bool {
        self.end == node
    }

    /// `getOther`
    pub fn other(&self, node: NodeId) -> Result<NodeId> {
        if self.is_beg(node) {
            Ok(self.end)
        } else if self.is_end(node) {
            Ok(self.beg)
        } else {
            Err(RDKitError::Runtime("Not an end-point of this edge!".into()))
        }
    }

    fn flip(&mut self) {
        std::mem::swap(&mut self.beg, &mut self.end);
    }
}

pub struct Digraph<'a> {
    mol: &'a CipMol,
    nodes: Vec<Node>,
    edges: Vec<Edge>,
    origin: NodeId,
    root: NodeId,
    atropisomer_mode: bool,
    rule6_ref: Option<u32>,
}

impl<'a> Digraph<'a> {
    pub fn new(mol: &'a CipMol, atom: u32, atropisomer_mode: bool) -> Result<Self> {
        let mut visit = vec![0u32; mol.num_atoms() as usize];
        visit[atom as usize] = 1;
        let atomic_num = Rational::new(mol.atom(atom).atomic_num() as i32, 1);
        let mut g = Digraph {
            mol,
            nodes: Vec::new(),
            edges: Vec::new(),
            origin: 0,
            root: 0,
            atropisomer_mode,
            rule6_ref: None,
        };
        let root = g.add_node(visit, Some(atom), atomic_num, 1, 0)?;
        g.origin = root;
        g.root = root;
        Ok(g)
    }

    pub fn mol(&self) -> &'a CipMol {
        self.mol
    }

    pub fn original_root(&self) -> NodeId {
        self.origin
    }

    pub fn current_root(&self) -> NodeId {
        self.root
    }

    pub fn rule6_ref(&self) -> Option<u32> {
        self.rule6_ref
    }

    pub fn set_rule6_ref(&mut self, atom: Option<u32>) {
        self.rule6_ref = atom;
    }

    pub fn node(&self, id: NodeId) -> &Node {
        &self.nodes[id]
    }

    pub fn edge(&self, id: EdgeId) -> &Edge {
        &self.edges[id]
    }

    pub fn set_node_aux(&mut self, id: NodeId, aux: Descriptor) {
        self.nodes[id].aux = aux;
    }

    /// The end node of an edge.
    pub fn end_of(&self, edge: EdgeId) -> &Node {
        &self.nodes[self.edges[edge].end]
    }

    /// The start node of an edge.
    pub fn beg_of(&self, edge: EdgeId) -> &Node {
        &self.nodes[self.edges[edge].beg]
    }

    fn add_node(
        &mut self,
        visit: Vec<u32>,
        atom: Option<u32>,
        frac: Rational,
        dist: i32,
        flags: u8,
    ) -> Result<NodeId> {
        if self.nodes.len() >= MAX_NODE_COUNT {
            return Err(RDKitError::Runtime(format!(
                "Digraph generation failed: more than {MAX_NODE_COUNT} nodes found."
            )));
        }
        self.nodes
            .push(Node::new(self.mol, visit, atom, frac, dist, flags));
        Ok(self.nodes.len() - 1)
    }

    /// `seenAtom`
    pub fn seen_atom(&self, atom: u32) -> bool {
        self.nodes.iter().any(|n| n.atom == Some(atom))
    }

    fn add_edge(&mut self, beg: NodeId, bond: Option<u32>, end: NodeId) {
        let id = self.edges.len();
        self.edges.push(Edge {
            beg,
            end,
            bond,
            aux: Descriptor::None,
        });
        self.nodes[beg].edges.push(id);
        self.nodes[end].edges.push(id);
    }

    /// `Node::getEdges`: the edges of a node, expanding it on first use.
    pub fn edges(&mut self, node: NodeId) -> Result<Vec<EdgeId>> {
        if !self.nodes[node].is_expanded() {
            self.nodes[node].flags |= EXPANDED;
            self.expand(node)?;
        }
        Ok(self.nodes[node].edges.clone())
    }

    /// `Node::getEdges(end)`: the non-duplicate edges that touch the atom.
    pub fn edges_to_atom(&mut self, node: NodeId, end: u32) -> Result<Vec<EdgeId>> {
        let mut res = Vec::new();
        for edge in self.edges(node)? {
            if self.end_of(edge).is_duplicate() {
                continue;
            }
            if Some(end) == self.beg_of(edge).atom || Some(end) == self.end_of(edge).atom {
                res.push(edge);
            }
        }
        Ok(res)
    }

    /// `Node::getNonTerminalOutEdges`
    pub fn non_terminal_out_edges(&mut self, node: NodeId) -> Result<Vec<EdgeId>> {
        let mut res = Vec::new();
        for edge in self.edges(node)? {
            if self.edges[edge].is_beg(node) && !self.end_of(edge).is_terminal() {
                res.push(edge);
            }
        }
        Ok(res)
    }

    /// `getNodes`: the nodes of an atom, in breadth first order.
    pub fn nodes_of_atom(&mut self, atom: u32) -> Result<Vec<NodeId>> {
        let mut result = Vec::new();
        let mut queue = vec![self.root];
        let mut i = 0;
        while i < queue.len() {
            let node = queue[i];
            i += 1;
            if self.nodes[node].atom == Some(atom) {
                result.push(node);
            }
            for e in self.edges(node)? {
                if self.edges[e].is_beg(node) {
                    queue.push(self.edges[e].end);
                }
            }
        }
        Ok(result)
    }

    /// `changeRoot`
    pub fn change_root(&mut self, new_root: NodeId) -> Result<()> {
        let mut to_flip = Vec::new();
        let mut queue = vec![new_root];
        let mut i = 0;
        while i < queue.len() {
            let node = queue[i];
            i += 1;
            for e in self.edges(node)? {
                if self.edges[e].is_end(node) {
                    to_flip.push(e);
                    queue.push(self.edges[e].beg);
                }
            }
        }
        for e in to_flip {
            self.edges[e].flip();
        }
        self.root = new_root;
        Ok(())
    }

    fn new_child(&mut self, parent: NodeId, atom: u32) -> Result<NodeId> {
        let mut visit = self.nodes[parent].visit.clone();
        let dist = self.nodes[parent].dist + 1;
        visit[atom as usize] = dist as u32;
        let atomic_num = Rational::new(self.mol.atom(atom).atomic_num() as i32, 1);
        self.add_node(visit, Some(atom), atomic_num, dist, 0)
    }

    fn new_terminal_child(
        &mut self,
        parent: NodeId,
        idx: Option<u32>,
        atom: Option<u32>,
        flags: u8,
    ) -> Result<NodeId> {
        let new_dist = if flags & DUPLICATE != 0 {
            self.nodes[parent].visit[idx.unwrap() as usize] as i32
        } else {
            self.nodes[parent].dist + 1
        };
        if flags & BOND_DUPLICATE != 0 {
            let frac = self
                .mol
                .fractional_atomic_num(self.nodes[parent].atom.unwrap())?;
            if frac.is_averaged() {
                return self.add_node(Vec::new(), atom, frac.value(), new_dist, flags);
            }
        }
        let atomic_num = atom.map_or(1, |a| self.mol.atom(a).atomic_num() as i32);
        self.add_node(
            Vec::new(),
            atom,
            Rational::new(atomic_num, 1),
            new_dist,
            flags,
        )
    }

    fn add_bond_duplicate(&mut self, beg: NodeId, bond: u32, nbr: u32) -> Result<()> {
        let end = self.new_terminal_child(beg, Some(nbr), Some(nbr), BOND_DUPLICATE)?;
        self.add_edge(beg, Some(bond), end);
        Ok(())
    }

    /// Adds the duplicate nodes standing for the higher order of a bond (one
    /// for an averaged anion, `virtual_nodes` otherwise).
    fn add_bond_order_duplicates(
        &mut self,
        beg: NodeId,
        atom: u32,
        bond: u32,
        nbr: u32,
        virtual_nodes: i32,
    ) -> Result<()> {
        if self.mol.atom(atom).formal_charge() < 0
            && self.mol.fractional_atomic_num(atom)?.is_averaged()
        {
            self.add_bond_duplicate(beg, bond, nbr)?;
        } else {
            for _ in 0..virtual_nodes {
                self.add_bond_duplicate(beg, bond, nbr)?;
            }
        }
        Ok(())
    }

    /// `expand`
    fn expand(&mut self, beg: NodeId) -> Result<()> {
        let mol = self.mol;
        let atom = self.nodes[beg].atom.unwrap();
        let prev = match self.nodes[beg].edges.first() {
            Some(&e) if !self.edges[e].is_beg(beg) => self.edges[e].bond,
            _ => None,
        };
        for bond in mol.bonds_of(atom).collect::<Vec<_>>() {
            let nbr = mol.other_atom(bond, atom)?;
            let virtual_nodes = mol.bond_order(bond)? - 1;
            if !self.nodes[beg].is_visited(nbr) {
                let end = self.new_child(beg, nbr)?;
                self.add_edge(beg, Some(bond), end);
                if self.origin != beg || self.atropisomer_mode {
                    self.add_bond_order_duplicates(beg, atom, bond, nbr, virtual_nodes)?;
                }
            } else if Some(bond) == prev {
                // bond order expansion (backwards)
                if self.nodes[self.origin].atom != Some(nbr) || self.atropisomer_mode {
                    for _ in 0..virtual_nodes {
                        self.add_bond_duplicate(beg, bond, nbr)?;
                    }
                }
            } else {
                // ring closure
                let end = self.new_terminal_child(beg, Some(nbr), Some(nbr), RING_DUPLICATE)?;
                self.add_edge(beg, Some(bond), end);
                self.add_bond_order_duplicates(beg, atom, bond, nbr, virtual_nodes)?;
            }
        }
        for _ in 0..mol.total_num_hs(atom)? {
            let end = self.new_terminal_child(beg, None, None, IMPL_HYDROGEN)?;
            self.add_edge(beg, None, end);
        }
        Ok(())
    }
}
