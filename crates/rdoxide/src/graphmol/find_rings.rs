//! Ported from RDKit `Code/GraphMol/FindRings.cpp` (+ `RingUtils` from
//! `RingInfo.cpp` and `MolOps::getMolFrags`).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::sync::Arc;
use std::sync::atomic::{AtomicU8, Ordering};

use super::bond::BondType;
use super::ring_info::{FindRingType, normalize_ring};
use super::romol::ROMol;
use crate::external::ring_decomposer_lib::{self as rdl, Cycle};
use crate::{RDKitError, Result};

const MAX_BFSQ_SIZE: usize = 200_000;

fn is_dative(bt: BondType) -> bool {
    matches!(
        bt,
        BondType::Dative | BondType::DativeL | BondType::DativeR | BondType::DativeOne
    )
}

fn excluded_bond(bt: BondType, include_dative: bool, include_hydrogen: bool) -> bool {
    bt == BondType::Zero
        || (!include_dative && is_dative(bt))
        || (!include_hydrogen && bt == BondType::Hydrogen)
}

/// `ringComparer`: by size, then lexicographically.
fn ring_cmp(a: &Vec<i32>, b: &Vec<i32>) -> std::cmp::Ordering {
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

// ------------------------------------------------------------------ ring invariants

/// `boost::dynamic_bitset` used as ring invariant; ordered like boost (as an
/// unsigned integer, highest bit most significant).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct RingInvar(Vec<u64>);

impl RingInvar {
    fn new(ring: &[i32], n: usize) -> Self {
        let mut v = vec![0u64; n.div_ceil(64).max(1)];
        for &i in ring {
            v[i as usize / 64] |= 1 << (i as usize % 64);
        }
        RingInvar(v)
    }
}

impl PartialOrd for RingInvar {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RingInvar {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.iter().rev().cmp(other.0.iter().rev())
    }
}

// ------------------------------------------------------------------ helpers

/// `RingUtils::convertToBonds`
pub fn convert_to_bonds(ring: &[i32], mol: &ROMol) -> Result<Vec<i32>> {
    let n = ring.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (a, b) = (ring[i] as u32, ring[(i + 1) % n] as u32);
        let bidx = mol
            .bond_idx_between(a, b)
            .ok_or_else(|| RDKitError::Value("expected bond not found".into()))?;
        out.push(bidx as i32);
    }
    Ok(out)
}

/// `RingUtils::rdlCycleToAtomRing`
pub fn rdl_cycle_to_atom_ring(cycle: &Cycle) -> Vec<i32> {
    let weight = cycle.len();
    let mut ring: Vec<i32> = Vec::with_capacity(weight);
    let mut unseen = vec![true; weight];
    ring.push(cycle[0][0] as i32);
    ring.push(cycle[0][1] as i32);
    unseen[0] = false;
    while ring.len() < weight {
        for edge_idx in 0..weight {
            if !unseen[edge_idx] {
                continue;
            }
            let edge = cycle[edge_idx];
            for j in 0..2 {
                if *ring.last().expect("non-empty") as u32 == edge[j] {
                    ring.push(edge[1 - j] as i32);
                    unseen[edge_idx] = false;
                    break;
                }
            }
            if ring.len() == weight {
                break;
            }
        }
    }
    normalize_ring(&mut ring);
    ring
}

/// `MolOps::getMolFrags(mol, mapping)`: connected-component id per atom,
/// numbered in order of the lowest atom index (boost::connected_components).
pub fn get_mol_frags_mapping(mol: &ROMol) -> (u32, Vec<i32>) {
    let n = mol.num_atoms() as usize;
    let mut mapping = vec![-1i32; n];
    let mut ncomp = 0;
    for start in 0..n {
        if mapping[start] != -1 {
            continue;
        }
        let mut stack = vec![start as u32];
        mapping[start] = ncomp;
        while let Some(a) = stack.pop() {
            for nbr in mol.atom_neighbors(a) {
                if mapping[nbr as usize] == -1 {
                    mapping[nbr as usize] = ncomp;
                    stack.push(nbr);
                }
            }
        }
        ncomp += 1;
    }
    (ncomp as u32, mapping)
}

/// `MolOps::getMolFrags(mol, frags)`
pub fn get_mol_frags(mol: &ROMol) -> Vec<Vec<i32>> {
    let (n, mapping) = get_mol_frags_mapping(mol);
    let mut frags = vec![Vec::new(); n as usize];
    for (i, &m) in mapping.iter().enumerate() {
        frags[m as usize].push(i as i32);
    }
    frags
}

fn store_rings_info(mol: &mut ROMol, rings: &[Vec<i32>]) -> Result<()> {
    for ring in rings {
        let bonds = convert_to_bonds(ring, mol)?;
        mol.ring_info_mut().add_ring(ring.clone(), bonds)?;
    }
    Ok(())
}

// ------------------------------------------------------------------ fast (DFS) ring finding

fn dfs(
    mol: &ROMol,
    atom: u32,
    colors: &mut [u8],
    order: &mut Vec<u32>,
    res: &mut Vec<Vec<i32>>,
    from: Option<u32>,
) {
    colors[atom as usize] = 1;
    order.push(atom);
    let nbrs: Vec<u32> = mol.atom_neighbors(atom).collect();
    for nbr in nbrs {
        match colors[nbr as usize] {
            0 => {
                if mol.atom_degree(nbr) < 2 {
                    colors[nbr as usize] = 2;
                } else {
                    dfs(mol, nbr, colors, order, res, Some(atom));
                }
            }
            1 if from.is_some_and(|f| f != nbr) => {
                let pos = order.iter().rposition(|&a| a == atom).expect("on stack");
                let mut cycle = Vec::new();
                for &a in order[..=pos].iter().rev() {
                    if a == nbr {
                        break;
                    }
                    cycle.push(a as i32);
                }
                cycle.push(nbr as i32);
                res.push(cycle);
            }
            _ => {}
        }
    }
    colors[atom as usize] = 2;
    order.pop();
}

/// `FindRings::fastFindRingsHeuristic`
fn fast_find_rings_heuristic(mol: &ROMol) -> Vec<Vec<i32>> {
    let n = mol.num_atoms() as usize;
    let mut colors = vec![0u8; n];
    let mut res = Vec::new();
    for i in 0..n as u32 {
        if colors[i as usize] != 0 {
            continue;
        }
        if mol.atom_degree(i) < 2 {
            colors[i as usize] = 2;
            continue;
        }
        let mut order = Vec::new();
        dfs(mol, i, &mut colors, &mut order, &mut res, None);
    }
    res
}

/// `MolOps::fastFindRings`
pub fn fast_find_rings(mol: &mut ROMol) -> Result<()> {
    if mol.ring_info().is_initialized() {
        mol.ring_info_mut().reset(true);
    }
    mol.ring_info_mut().initialize(FindRingType::Fast);
    let res = fast_find_rings_heuristic(mol);
    store_rings_info(mol, &res)
}

// ------------------------------------------------------------------ ring families / SSSR

/// `MolOps::findRingFamilies`
pub fn find_ring_families(
    mol: &mut ROMol,
    include_dative: bool,
    include_hydrogen: bool,
) -> Result<()> {
    if mol.ring_info().is_initialized() {
        if mol.ring_info().are_ring_families_initialized() {
            return Ok(());
        }
    } else {
        mol.ring_info_mut().initialize(FindRingType::OtherOrUnknown);
    }
    let num_atoms = (mol.num_atoms() as usize).max(1);
    let mut graph = rdl::Graph::new(num_atoms);
    for b in mol.bonds() {
        if excluded_bond(b.bond_type(), include_dative, include_hydrogen) {
            continue;
        }
        graph.add_u_edge(b.begin_atom_idx(), b.end_atom_idx());
    }
    let Some(data) = rdl::calculate(graph) else {
        mol.ring_info_mut().set_urf_data(None);
        return Err(RDKitError::Value("Cannot get URFs".into()));
    };
    let mut families = Vec::with_capacity(data.num_urfs());
    for i in 0..data.num_urfs() {
        let nodes: Vec<i32> = data
            .nodes_for_urf(i)
            .into_iter()
            .map(|n| n as i32)
            .collect();
        let edges: Result<Vec<i32>> = data
            .edges_for_urf(i)
            .into_iter()
            .map(|[a, b]| {
                mol.bond_idx_between(a, b)
                    .map(|x| x as i32)
                    .ok_or_else(|| RDKitError::Value("Cannot get URF edges".into()))
            })
            .collect();
        families.push((nodes, edges?));
    }
    let ri = mol.ring_info_mut();
    ri.set_urf_data(Some(Arc::new(data)));
    for (n, e) in families {
        ri.add_ring_family(n, e)?;
    }
    Ok(())
}

/// `MolOps::findSSSR`: returns the rings (sorted by size, then indices).
pub fn find_sssr(
    mol: &mut ROMol,
    include_dative: bool,
    include_hydrogen: bool,
) -> Result<Vec<Vec<i32>>> {
    if mol.ring_info().is_initialized() {
        mol.ring_info_mut().reset(true);
    }
    mol.ring_info_mut().initialize(FindRingType::Sssr);
    let (na, nb) = (mol.num_atoms() as usize, mol.num_bonds() as usize);
    mol.ring_info_mut().preallocate(na, nb);
    find_ring_families(mol, include_dative, include_hydrogen)?;
    let data = mol.ring_info().urf_data().expect("ring families computed");
    let mut res: Vec<Vec<i32>> = data.sssr().iter().map(rdl_cycle_to_atom_ring).collect();
    res.sort_by(ring_cmp);
    store_rings_info(mol, &res)?;
    Ok(res)
}

// ------------------------------------------------------------------ legacy Figueras

struct Figueras<'a> {
    mol: &'a ROMol,
    n_atoms: usize,
}

impl Figueras<'_> {
    fn trim_bonds(
        &self,
        cand: u32,
        changed: &mut VecDeque<u32>,
        degrees: &mut [i32],
        active: &mut [bool],
    ) {
        for b in self.mol.atom_bonds(cand) {
            let bi = b.idx() as usize;
            if !active[bi] {
                continue;
            }
            let o = b.other_atom_idx(cand).expect("in bond");
            if degrees[o as usize] <= 2 {
                changed.push_back(o);
            }
            active[bi] = false;
            degrees[o as usize] -= 1;
            degrees[cand as usize] -= 1;
        }
    }

    fn smallest_rings_bfs(
        &self,
        root: u32,
        rings: &mut Vec<Vec<i32>>,
        active: &[bool],
        forbidden: Option<&[i32]>,
    ) -> Result<usize> {
        const WHITE: u8 = 0;
        const GRAY: u8 = 1;
        const BLACK: u8 = 2;
        let n = self.n_atoms;
        let mut done = vec![WHITE; n];
        if let Some(f) = forbidden {
            for &i in f {
                done[i as usize] = BLACK;
            }
        }
        let mut parents = vec![-1i32; n];
        let mut depths = vec![0u32; n];
        let mut bfsq: VecDeque<i32> = VecDeque::from([root as i32]);
        let mut cur_size = u32::MAX;
        while let Some(curr) = bfsq.front().copied() {
            if bfsq.len() >= MAX_BFSQ_SIZE {
                let msg = "Maximum BFS search size exceeded.\nThis is likely due to a highly symmetric fused ring system.";
                crate::rd_log_error!("{msg}");
                return Err(RDKitError::Value(msg.into()));
            }
            bfsq.pop_front();
            done[curr as usize] = BLACK;
            let depth = depths[curr as usize] + 1;
            if depth > cur_size {
                break;
            }
            for b in self.mol.atom_bonds(curr as u32) {
                if !active[b.idx() as usize] {
                    continue;
                }
                let nbr = b.other_atom_idx(curr as u32).expect("in bond") as i32;
                if done[nbr as usize] == BLACK || parents[curr as usize] == nbr {
                    continue;
                }
                if done[nbr as usize] == WHITE {
                    parents[nbr as usize] = curr;
                    done[nbr as usize] = GRAY;
                    depths[nbr as usize] = depth;
                    bfsq.push_back(nbr);
                } else {
                    let mut ring = vec![nbr];
                    let mut parent = parents[nbr as usize];
                    while parent != -1 && parent != root as i32 {
                        ring.push(parent);
                        parent = parents[parent as usize];
                    }
                    ring.insert(0, curr);
                    parent = parents[curr as usize];
                    while parent != -1 {
                        if ring.contains(&parent) {
                            ring.clear();
                            break;
                        }
                        ring.insert(0, parent);
                        parent = parents[parent as usize];
                    }
                    if ring.len() > 1 {
                        if ring.len() as u32 <= cur_size {
                            cur_size = ring.len() as u32;
                            rings.push(ring);
                        } else {
                            return Ok(rings.len());
                        }
                    }
                }
            }
        }
        Ok(rings.len())
    }

    fn mark_useless_d2s(&self, root: u32, forb: &mut [bool], degrees: &[i32], active: &[bool]) {
        for b in self.mol.atom_bonds(root) {
            if !active[b.idx() as usize] {
                continue;
            }
            let o = b.other_atom_idx(root).expect("in bond");
            if !forb[o as usize] && degrees[o as usize] == 2 {
                forb[o as usize] = true;
                self.mark_useless_d2s(o, forb, degrees, active);
            }
        }
    }

    fn pick_d2_nodes(&self, frag: &[i32], degrees: &[i32], active: &[bool]) -> Vec<i32> {
        let mut d2 = Vec::new();
        let mut forb = vec![false; self.n_atoms];
        while let Some(&root) = frag
            .iter()
            .find(|&&a| degrees[a as usize] == 2 && !forb[a as usize])
        {
            d2.push(root);
            forb[root as usize] = true;
            self.mark_useless_d2s(root as u32, &mut forb, degrees, active);
        }
        d2
    }

    fn add_ring_bonds(&self, ring: &[i32], ring_bonds: &mut [bool], ring_atoms: &mut [bool]) {
        for i in 0..ring.len() - 1 {
            let b = self
                .mol
                .bond_idx_between(ring[i] as u32, ring[i + 1] as u32)
                .expect("ring bond");
            ring_bonds[b as usize] = true;
            ring_atoms[ring[i] as usize] = true;
        }
        let b = self
            .mol
            .bond_idx_between(ring[0] as u32, ring[ring.len() - 1] as u32)
            .expect("ring bond");
        ring_bonds[b as usize] = true;
        ring_atoms[ring[ring.len() - 1] as usize] = true;
    }

    #[allow(clippy::too_many_arguments)]
    fn find_rings_d2_nodes(
        &self,
        res: &mut Vec<Vec<i32>>,
        invars: &mut BTreeSet<RingInvar>,
        d2nodes: &[i32],
        degrees: &mut [i32],
        active: &mut [bool],
        ring_bonds: &mut [bool],
        ring_atoms: &mut [bool],
    ) -> Result<()> {
        let mut dup_d2: BTreeMap<RingInvar, Vec<i32>> = BTreeMap::new();
        let mut dup_map: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
        for &cand in d2nodes {
            let mut srings = Vec::new();
            self.smallest_rings_bfs(cand as u32, &mut srings, active, None)?;
            for nring in &srings {
                let invr = RingInvar::new(nring, self.n_atoms);
                let dups = dup_d2.entry(invr.clone()).or_default();
                if !invars.contains(&invr) {
                    res.push(nring.clone());
                    invars.insert(invr);
                    self.add_ring_bonds(nring, ring_bonds, ring_atoms);
                } else {
                    for &other in dups.iter() {
                        dup_map.entry(cand).or_default().push(other);
                        dup_map.entry(other).or_default().push(cand);
                    }
                }
                dups.push(cand);
            }
            if srings.is_empty() {
                let mut changed = VecDeque::from([cand as u32]);
                while let Some(c) = changed.pop_front() {
                    self.trim_bonds(c, &mut changed, degrees, active);
                }
            }
        }
        // findSSSRforDupCands
        for dup_cands in dup_d2.values() {
            if dup_cands.len() <= 1 {
                continue;
            }
            let mut nrings: Vec<Vec<i32>> = Vec::new();
            let mut min_siz = i32::MAX as usize;
            for &dc in dup_cands {
                let mut deg_copy = degrees.to_vec();
                let mut active_copy = active.to_vec();
                let mut changed = VecDeque::new();
                let dni = dup_map
                    .get(&dc)
                    .ok_or_else(|| RDKitError::invariant("duplicate could not be found"))?;
                for &d in dni {
                    self.trim_bonds(d as u32, &mut changed, &mut deg_copy, &mut active_copy);
                }
                let mut srings = Vec::new();
                self.smallest_rings_bfs(dc as u32, &mut srings, &active_copy, None)?;
                for s in srings {
                    min_siz = min_siz.min(s.len());
                    nrings.push(s);
                }
            }
            for nring in nrings {
                if nring.len() == min_siz {
                    let invr = RingInvar::new(&nring, self.n_atoms);
                    if !invars.contains(&invr) {
                        res.push(nring);
                        invars.insert(invr);
                    }
                }
            }
        }
        Ok(())
    }

    fn push_new(
        &self,
        rings: Vec<Vec<i32>>,
        res: &mut Vec<Vec<i32>>,
        invars: &mut BTreeSet<RingInvar>,
    ) {
        for nring in rings {
            let invr = RingInvar::new(&nring, self.n_atoms);
            if invars.insert(invr) {
                res.push(nring);
            }
        }
    }

    fn find_rings_d3_node(
        &self,
        res: &mut Vec<Vec<i32>>,
        invars: &mut BTreeSet<RingInvar>,
        cand: i32,
        active: &[bool],
    ) -> Result<()> {
        let mut srings = Vec::new();
        let nsmall = self.smallest_rings_bfs(cand as u32, &mut srings, active, None)?;
        self.push_new(srings.clone(), res, invars);
        if nsmall >= 3 {
            return Ok(());
        }
        let nbrs: Vec<i32> = self
            .mol
            .atom_bonds(cand as u32)
            .filter(|b| active[b.idx() as usize])
            .take(3)
            .map(|b| b.other_atom_idx(cand as u32).expect("in bond") as i32)
            .collect();
        if nbrs.len() < 3 {
            return Err(RDKitError::invariant("neighbor not found"));
        }
        let (n1, n2, n3) = (nbrs[0], nbrs[1], nbrs[2]);
        if nsmall == 2 {
            let in_both = |n: i32| srings[0].contains(&n) && srings[1].contains(&n);
            let f = if in_both(n1) {
                n1
            } else if in_both(n2) {
                n2
            } else if in_both(n3) {
                n3
            } else {
                return Err(RDKitError::invariant("third ring not found"));
            };
            let mut trings = Vec::new();
            self.smallest_rings_bfs(cand as u32, &mut trings, active, Some(&[f]))?;
            self.push_new(trings, res, invars);
        }
        if nsmall == 1 {
            let (f1, f2) = if !srings[0].contains(&n1) {
                (n2, n3)
            } else if !srings[0].contains(&n2) {
                (n1, n3)
            } else if !srings[0].contains(&n3) {
                (n1, n2)
            } else {
                return Err(RDKitError::invariant("rings not found"));
            };
            let mut trings = Vec::new();
            self.smallest_rings_bfs(cand as u32, &mut trings, active, Some(&[f2]))?;
            self.push_new(trings, res, invars);
            let mut trings = Vec::new();
            self.smallest_rings_bfs(cand as u32, &mut trings, active, Some(&[f1]))?;
            self.push_new(trings, res, invars);
        }
        Ok(())
    }

    fn atom_search_bfs(
        &self,
        start: u32,
        end: u32,
        ring_atoms: &[bool],
        invars: &BTreeSet<RingInvar>,
    ) -> Result<Option<Vec<i32>>> {
        let mut bfsq: VecDeque<Vec<i32>> = VecDeque::from([vec![start as i32]]);
        while let Some(tv) = bfsq.pop_front() {
            if bfsq.len() + 1 >= MAX_BFSQ_SIZE {
                let msg = "Maximum BFS search size exceeded.\nThis is likely due to a highly symmetric fused ring system.";
                crate::rd_log_error!("{msg}");
                return Err(RDKitError::Value(msg.into()));
            }
            let curr = *tv.last().expect("non-empty") as u32;
            for nbr in self.mol.atom_neighbors(curr) {
                if nbr == end {
                    if curr != start {
                        let mut nv = tv.clone();
                        nv.push(nbr as i32);
                        if !invars.contains(&RingInvar::new(&nv, self.n_atoms)) {
                            return Ok(Some(nv));
                        }
                    }
                } else if ring_atoms[nbr as usize] && !tv.contains(&(nbr as i32)) {
                    let mut nv = tv.clone();
                    nv.push(nbr as i32);
                    bfsq.push_back(nv);
                }
            }
        }
        Ok(None)
    }

    fn run(&self, include_dative: bool, include_hydrogen: bool) -> Result<Vec<Vec<i32>>> {
        let mol = self.mol;
        let mut res: Vec<Vec<i32>> = Vec::new();
        let nb_total = mol.num_bonds() as usize;
        let mut active = vec![true; nb_total];
        for b in mol.bonds() {
            if excluded_bond(b.bond_type(), include_dative, include_hydrogen) {
                active[b.idx() as usize] = false;
            }
        }
        let nats = self.n_atoms;
        let mut degrees = vec![0i32; nats];
        let mut degrees_with_zero = vec![0i32; nats];
        for i in 0..nats {
            let deg = mol.atom_degree(i as u32) as i32;
            degrees[i] = deg;
            degrees_with_zero[i] = deg;
            for b in mol.atom_bonds(i as u32) {
                if excluded_bond(b.bond_type(), include_dative, include_hydrogen) {
                    degrees[i] -= 1;
                }
            }
        }
        let mut invars: BTreeSet<RingInvar> = BTreeSet::new();
        let mut ring_atoms = vec![false; nats];
        let mut ring_bonds = vec![false; nb_total];
        for frag in get_mol_frags(mol) {
            if frag.len() < 3 {
                continue;
            }
            let mut changed: VecDeque<u32> = VecDeque::new();
            let mut bnd_with_zero = 0i32;
            let mut nbnds = 0u32;
            for &a in &frag {
                bnd_with_zero += degrees_with_zero[a as usize];
                let deg = degrees[a as usize];
                nbnds += deg as u32;
                if deg < 2 {
                    changed.push_back(a as u32);
                }
            }
            bnd_with_zero /= 2;
            if bnd_with_zero - frag.len() as i32 + 1 < 1 {
                continue;
            }
            let nbnds = nbnds / 2;
            let mut done_ats = vec![false; nats];
            let mut n_done = 0usize;
            let mut frag_res: Vec<Vec<i32>> = Vec::new();
            while n_done + 3 <= frag.len() {
                while let Some(cand) = changed.pop_front() {
                    if !done_ats[cand as usize] {
                        done_ats[cand as usize] = true;
                        n_done += 1;
                        self.trim_bonds(cand, &mut changed, &mut degrees, &mut active);
                    }
                }
                let d2 = self.pick_d2_nodes(&frag, &degrees, &active);
                if !d2.is_empty() {
                    self.find_rings_d2_nodes(
                        &mut frag_res,
                        &mut invars,
                        &d2,
                        &mut degrees,
                        &mut active,
                        &mut ring_bonds,
                        &mut ring_atoms,
                    )?;
                    for &d in &d2 {
                        done_ats[d as usize] = true;
                        n_done += 1;
                        self.trim_bonds(d as u32, &mut changed, &mut degrees, &mut active);
                    }
                } else if n_done + 3 <= frag.len() {
                    let Some(&cand) = frag.iter().find(|&&a| degrees[a as usize] == 3) else {
                        break;
                    };
                    self.find_rings_d3_node(&mut frag_res, &mut invars, cand, &active)?;
                    done_ats[cand as usize] = true;
                    n_done += 1;
                    self.trim_bonds(cand as u32, &mut changed, &mut degrees, &mut active);
                }
            }
            let nexpt = nbnds as i64 - frag.len() as i64 + 1;
            if (frag_res.len() as i64) < nexpt {
                // NB: upstream scans bond indices below this fragment's bond count
                let first_possible =
                    |ring_bonds: &[bool], ring_atoms: &[bool], dead: &[bool]| -> Option<u32> {
                        (0..nbnds).find(|&i| {
                            let b = mol.bond(i);
                            !ring_bonds[i as usize]
                                && !dead[i as usize]
                                && ring_atoms[b.begin_atom_idx() as usize]
                                && ring_atoms[b.end_atom_idx() as usize]
                        })
                    };
                let mut dead = vec![false; nb_total];
                let mut possible = first_possible(&ring_bonds, &ring_atoms, &dead);
                while let Some(bi) = possible {
                    let b = mol.bond(bi);
                    let found = match self.atom_search_bfs(
                        b.begin_atom_idx(),
                        b.end_atom_idx(),
                        &ring_atoms,
                        &invars,
                    )? {
                        Some(nring) => {
                            let invr = RingInvar::new(&nring, nats);
                            if invars.insert(invr) {
                                self.add_ring_bonds(&nring, &mut ring_bonds, &mut ring_atoms);
                                frag_res.push(nring);
                            }
                            true
                        }
                        None => false,
                    };
                    if !found {
                        dead[bi as usize] = true;
                    }
                    possible = first_possible(&ring_bonds, &ring_atoms, &dead);
                }
                if (frag_res.len() as i64) < nexpt {
                    crate::rd_log_warning!(
                        "WARNING: could not find number of expected rings. Switching to an approximate ring finding algorithm."
                    );
                    return Ok(fast_find_rings_heuristic(mol));
                }
            }
            res.extend(frag_res);
        }
        for r in &mut res {
            normalize_ring(r);
        }
        res.sort_by(ring_cmp);
        Ok(res)
    }
}

// ------------------------------------------------------------------ symmetrizeSSSR

/// `MolOps::SymmetrizeSSSRAlgorithm`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymmetrizeSSSRAlgorithm {
    Default,
    Legacy,
    Rdl,
}

static USE_LEGACY_RING_FINDING: AtomicU8 = AtomicU8::new(0);

/// `MolOps::setUseLegacyRingFinding`
pub fn set_use_legacy_ring_finding(val: bool) {
    USE_LEGACY_RING_FINDING.store(if val { 2 } else { 1 }, Ordering::Relaxed);
}

/// `MolOps::getUseLegacyRingFinding` (default false; env var
/// `RDK_USE_LEGACY_RING_FINDING`)
pub fn get_use_legacy_ring_finding() -> bool {
    match USE_LEGACY_RING_FINDING.load(Ordering::Relaxed) {
        1 => false,
        2 => true,
        _ => std::env::var_os("RDK_USE_LEGACY_RING_FINDING").is_some_and(|v| v != "0"),
    }
}

fn legacy_symmetrize_sssr(
    mol: &mut ROMol,
    include_dative: bool,
    include_hydrogen: bool,
) -> Result<Vec<Vec<i32>>> {
    let mut res = find_sssr(mol, include_dative, include_hydrogen)?;
    let fig = Figueras {
        mol,
        n_atoms: mol.num_atoms() as usize,
    }
    .run(include_dative, include_hydrogen)?;
    // std::set_difference on sorted ranges with ringComparer
    let mut extras = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < fig.len() {
        if j >= res.len() || ring_cmp(&fig[i], &res[j]).is_lt() {
            extras.push(fig[i].clone());
            i += 1;
        } else if ring_cmp(&res[j], &fig[i]).is_lt() {
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    if extras.is_empty() {
        return Ok(res);
    }
    let bondsssrs: Vec<Vec<i32>> = mol.ring_info().bond_rings().to_vec();
    let mut bond_counts = vec![0i32; mol.num_bonds() as usize];
    for r in &bondsssrs {
        for &b in r {
            bond_counts[b as usize] += 1;
        }
    }
    for extra in extras {
        let extra_bonds = convert_to_bonds(&extra, mol)?;
        for ring in &bondsssrs {
            if ring.len() != extra_bonds.len() {
                continue;
            }
            let mut share = false;
            let mut replaces_all_unique = true;
            for &bid in ring {
                let count = bond_counts[bid as usize];
                if count == 1 || !share {
                    if extra_bonds.contains(&bid) {
                        share = true;
                    } else if count == 1 {
                        replaces_all_unique = false;
                        break;
                    }
                }
            }
            if share && replaces_all_unique {
                res.push(extra.clone());
                mol.ring_info_mut()
                    .add_ring(extra.clone(), extra_bonds.clone())?;
                break;
            }
        }
    }
    Ok(res)
}

/// `MolOps::symmetrizeSSSR`
pub fn symmetrize_sssr(
    mol: &mut ROMol,
    algorithm: SymmetrizeSSSRAlgorithm,
    recalc_sssr: bool,
    include_dative: bool,
    include_hydrogen: bool,
) -> Result<Vec<Vec<i32>>> {
    if mol.ring_info().is_initialized() {
        mol.ring_info_mut().reset(recalc_sssr);
    }
    mol.ring_info_mut().initialize(FindRingType::SymmSssr);
    let algorithm = match algorithm {
        SymmetrizeSSSRAlgorithm::Default => {
            if get_use_legacy_ring_finding() {
                SymmetrizeSSSRAlgorithm::Legacy
            } else {
                SymmetrizeSSSRAlgorithm::Rdl
            }
        }
        a => a,
    };
    if algorithm != SymmetrizeSSSRAlgorithm::Legacy {
        if !mol.ring_info().are_ring_families_initialized() {
            let (na, nb) = (mol.num_atoms() as usize, mol.num_bonds() as usize);
            mol.ring_info_mut().preallocate(na, nb);
            find_ring_families(mol, include_dative, include_hydrogen)?;
        }
        let res = mol.ring_info().atom_relevant_cycles()?;
        for ring in &res {
            let bonds = convert_to_bonds(ring, mol)?;
            mol.ring_info_mut().add_ring(ring.clone(), bonds)?;
        }
        Ok(res)
    } else {
        legacy_symmetrize_sssr(mol, include_dative, include_hydrogen)
    }
}

/// Atom rings of the relevant cycles stored in RDL data (`RingInfo::atomRelevantCycles`).
pub fn atom_relevant_cycles(data: &rdl::Data) -> Vec<Vec<i32>> {
    data.relevant_cycles()
        .iter()
        .map(rdl_cycle_to_atom_ring)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    fn mol(s: &str) -> ROMol {
        let p = SmilesParserParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        mol_from_smiles(s, &p).unwrap().unwrap()
    }

    #[test]
    fn sssr_basics() {
        let mut m = mol("c1ccc2ccccc2c1");
        let r = find_sssr(&mut m, false, false).unwrap();
        assert_eq!(r.len(), 2);
        assert!(r.iter().all(|x| x.len() == 6));
        assert_eq!(m.ring_info().num_rings().unwrap(), 2);
        let mut m = mol("CCO");
        assert!(find_sssr(&mut m, false, false).unwrap().is_empty());
    }

    #[test]
    fn symmetrized_cubane() {
        let mut m = mol("C12C3C4C1C5C2C3C45");
        let r =
            symmetrize_sssr(&mut m, SymmetrizeSSSRAlgorithm::Default, true, false, false).unwrap();
        assert_eq!(r.len(), 6);
        let mut m2 = mol("C12C3C4C1C5C2C3C45");
        let r2 =
            symmetrize_sssr(&mut m2, SymmetrizeSSSRAlgorithm::Legacy, true, false, false).unwrap();
        assert_eq!(r2.len(), 6);
        let mut m3 = mol("C12C3C4C1C5C2C3C45");
        assert_eq!(find_sssr(&mut m3, false, false).unwrap().len(), 5);
    }

    #[test]
    fn fast_rings_and_frags() {
        let mut m = mol("C1CC1.C1CCC1");
        fast_find_rings(&mut m).unwrap();
        assert_eq!(m.ring_info().num_rings().unwrap(), 2);
        assert_eq!(get_mol_frags(&m), vec![vec![0, 1, 2], vec![3, 4, 5, 6]]);
    }

    #[test]
    fn invar_order_matches_boost() {
        // boost::dynamic_bitset compares as an unsigned integer
        let a = RingInvar::new(&[0], 70);
        let b = RingInvar::new(&[65], 70);
        assert!(a < b);
    }
}
