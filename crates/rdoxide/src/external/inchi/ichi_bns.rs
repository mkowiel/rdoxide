//! The balanced network structure (BNS) of InChI (`INCHI_BASE/src/ichi_bns.c`,
//! `ichi_bns.h`): data structures, allocation, edge/vertex access, flow
//! augmentation and residual capacities. The search itself
//! (`BalancedNetworkSearch`) follows in the next stage.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! Build facts: `BNS_RAD_SEARCH` is 1, `FIX_NUM_TG`, `FIX_TACN_POSSIBLE_BUG`,
//! `FIX_KEEP_H_ON_NH_ANION` and `FIX_AVOID_ADP` are 0,
//! `RESET_EDGE_FORBIDDEN_MASK` is 1 in `ichi_bns.c` (so `IS_FORBIDDEN(x)` is
//! `x`), `ADD_CAPACITY_RADICAL` is 1. The `iedge` pointers of the vertices
//! are offsets into [`BnStruct::iedge`]; the `alt_path` pointer is an index
//! into [`BnStruct::altp`].

use super::ichierr::{BNS_ALTPATH_OVFL, BNS_WRONG_PARMS, is_bns_error};
use super::ichitaut::is_centerpoint_elem;
use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, RADICAL_DOUBLET};
use super::inpdef::{BOND_TYPE_MASK, InpAtom};
use super::util::get_endpoint_valence;

pub type Vertex = i32;
pub type EdgeIndex = i32;
/// `Edge`: `[vertex1, iedge or -(1+vertex1) if vertex2 is s or t]`
pub type Edge = [i32; 2];

pub const BN_MAX_ALTP: usize = 16;
pub const NUM_KINDS_OF_GROUPS: i32 = 2;
pub const BNS_ADD_ATOMS: i32 = 2;
pub const BNS_ADD_EDGES: i32 = 1;

pub const VERTEX_S: Vertex = 0;
pub const VERTEX_T: Vertex = 1;
pub const NO_VERTEX: Vertex = -2;
pub const BLOSSOM_BASE: Vertex = -1;

pub const MAX_BOND_EDGE_CAP: i32 = 2;
pub const AROM_BOND_EDGE_CAP: i32 = 1;

pub const EDGE_FLOW_ST_MASK: i32 = 0x3fff;
pub const EDGE_FLOW_ST_PATH: i32 = 0x4000;
pub const EDGE_FLOW_MASK: i32 = 0x3fff;
pub const EDGE_FLOW_PATH: i32 = 0x4000;

pub const BNS_EF_CHNG_FLOW: i32 = 1;
pub const BNS_EF_RSTR_FLOW: i32 = 2;
pub const BNS_EF_CHNG_BONDS: i32 = 4;
pub const BNS_EF_ALTR_BONDS: i32 = 8;
pub const BNS_EF_UPD_RAD_ORI: i32 = 16;
pub const BNS_EF_SET_NOSTEREO: i32 = 32;
pub const BNS_EF_UPD_H_CHARGE: i32 = 64;
pub const BNS_EF_RAD_SRCH: i32 = 128;

pub const BNS_EDGE_FORBIDDEN_MASK: i8 = 1;
pub const BNS_EDGE_FORBIDDEN_TEMP: i8 = 2;
pub const BNS_EDGE_FORBIDDEN_TEST: i8 = 4;

pub const BNS_VERT_TYPE_ATOM: u16 = 0x0001;
pub const BNS_VERT_TYPE_ENDPOINT: u16 = 0x0002;
pub const BNS_VERT_TYPE_TGROUP: u16 = 0x0004;
pub const BNS_VERT_TYPE_C_POINT: u16 = 0x0008;
pub const BNS_VERT_TYPE_C_GROUP: u16 = 0x0010;
pub const BNS_VERT_TYPE_SUPER_TGROUP: u16 = 0x0020;
pub const BNS_VERT_TYPE_TEMP: u16 = 0x0040;
pub const BNS_VERT_TYPE_AUX: u16 = 0x0080;
pub const BNS_VERT_TYPE_C_NEGATIVE: u16 = 0x0100;
pub const BNS_VERT_TYPE_ACID: u16 = 0x0200;
pub const BNS_VERT_TYPE_CARBON_GR: u16 = 0x0400;
pub const BNS_VERT_TYPE_METAL_GR: u16 = 0x0800;

// `ALT_CONST`
pub const I_ALTP_MAX_LEN: usize = 0;
pub const I_ALTP_FLOW: usize = 1;
pub const I_ALTP_PATH_LEN: usize = 2;
pub const I_ALTP_START_ATOM: usize = 3;
pub const I_ALTP_END_ATOM: usize = 4;
pub const I_ALTP_NEIGHBOR: usize = 5;
pub const I_ALTP_HDR_LEN: usize = I_ALTP_NEIGHBOR;

pub const TREE_NOT_IN_M: i8 = 0;
pub const TREE_IN_2: i8 = 1;
pub const TREE_IN_2BLOSS: i8 = 2;

/// `BNS_EDGE`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BnsEdge {
    pub neighbor1: u16,
    pub neighbor12: u16,
    pub neigh_ord: [u16; 2],
    pub cap: i32,
    pub cap0: i32,
    pub flow: i32,
    pub flow0: i32,
    pub pass: i8,
    pub forbidden: i8,
}

/// `BNS_ST_EDGE`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BnsStEdge {
    pub cap: i32,
    pub cap0: i32,
    pub flow: i32,
    pub flow0: i32,
    pub pass: i8,
}

/// `BNS_VERTEX`; `iedge` is the offset of the vertex's edge indexes in
/// [`BnStruct::iedge`]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BnsVertex {
    pub st_edge: BnsStEdge,
    pub type_: u16,
    pub num_adj_edges: u16,
    pub max_adj_edges: u16,
    pub iedge: usize,
}

/// `BNS_ALT_PATH`: a union of `flow[2]`, `number` and `ineigh[2]` (the
/// `AT_NUMB` pair overlays `flow[0]`, little endian)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BnsAltPath {
    pub flow: [i32; 2],
}

impl BnsAltPath {
    pub fn number(&self) -> i32 {
        self.flow[0]
    }
    pub fn set_number(&mut self, v: i32) {
        self.flow[0] = v;
    }
    pub fn ineigh(&self, i: usize) -> u16 {
        ((self.flow[0] as u32) >> (16 * i)) as u16
    }
    pub fn set_ineigh(&mut self, i: usize, v: u16) {
        let mut u = self.flow[0] as u32;
        u &= !(0xffffu32 << (16 * i));
        u |= u32::from(v) << (16 * i);
        self.flow[0] = u as i32;
    }
}

/// `BN_STRUCT`
#[derive(Debug, Clone, Default)]
pub struct BnStruct {
    pub num_atoms: i32,
    pub num_added_atoms: i32,
    pub n_max_add_atoms: i32,
    pub num_c_groups: i32,
    pub num_t_groups: i32,
    pub num_vertices: i32,
    pub num_bonds: i32,
    pub num_edges: i32,
    pub num_iedges: i32,
    pub num_added_edges: i32,
    pub n_max_add_edges: i32,
    pub max_vertices: i32,
    pub max_edges: i32,
    pub max_iedges: i32,
    pub tot_st_cap: i32,
    pub tot_st_flow: i32,
    pub len_alt_path: i32,
    pub b_not_a_simple_path: i32,
    pub b_change_flow: i32,
    pub vert: Vec<BnsVertex>,
    pub edge: Vec<BnsEdge>,
    pub iedge: Vec<i32>,
    /// index into `altp` of the current alt. path (`alt_path`)
    pub alt_path: Option<usize>,
    pub altp: Vec<Vec<BnsAltPath>>,
    pub max_altp: i32,
    pub num_altp: i32,
    pub type_tacn: u16,
    pub type_t: u16,
    pub type_cn: u16,
    pub edge_forbidden_mask: i8,
}

impl BnStruct {
    fn ap(&mut self) -> &mut Vec<BnsAltPath> {
        let k = self.alt_path.expect("alt_path is set");
        &mut self.altp[k]
    }
}

/// `BN_DATA`
#[derive(Debug, Clone, Default)]
pub struct BnData {
    pub base_ptr: Vec<Vertex>,
    pub switch_edge: Vec<Edge>,
    pub tree: Vec<i8>,
    pub scan_q: Vec<Vertex>,
    pub q_size: i32,
    pub pu: Vec<Vertex>,
    pub pv: Vec<Vertex>,
    pub max_num_vertices: i32,
    pub max_len_pu_pv: i32,
    pub rad_endpoints: Vec<Vertex>,
    pub n_num_rad_endpoints: i32,
    pub rad_edges: Vec<EdgeIndex>,
    pub n_num_rad_edges: i32,
    pub n_num_radicals: i32,
    pub b_rad_srch_mode: i32,
}

fn max_at_flow(a: &InpAtom) -> i32 {
    let rad = if is_centerpoint_elem(a.el_number) || get_endpoint_valence(a.el_number) != 0 {
        i32::from(i32::from(a.radical) == RADICAL_DOUBLET)
    } else {
        0
    };
    (i32::from(a.chem_bonds_valence) - i32::from(a.valence)) + rad
}

/// `AllocateAndInitBnStruct`: builds the network of the atoms (each bond is an
/// edge; the group vertices are added later). Alternating and unknown bonds
/// of `at` are made single; their number is returned as the second value.
pub fn allocate_and_init_bn_struct(
    at: &mut [InpAtom],
    num_atoms: usize,
    n_max_add_atoms: i32,
    mut n_max_add_edges: i32,
    max_altp: i32,
) -> (BnStruct, i32) {
    let mut num_changed_bonds = 0;
    n_max_add_edges += 1; // BNS_RAD_SEARCH
    let max_tg = num_atoms as i32; // FIX_NUM_TG == 0
    let mut num_bonds: i32 = 0;
    let mut num_rad = 0;
    for a in &at[..num_atoms] {
        num_bonds += i32::from(a.valence);
        num_rad += i32::from(i32::from(a.radical) == RADICAL_DOUBLET);
    }
    let num_edges = num_bonds / 2;
    num_bonds = num_edges;
    let mut max_vertices = num_atoms as i32 + n_max_add_atoms + max_tg + 1;
    let mut max_edges = num_edges + (n_max_add_edges + NUM_KINDS_OF_GROUPS) * max_vertices + max_tg;
    if num_rad != 0 {
        max_vertices *= 2;
        max_edges *= 2;
    }
    let max_iedges = 2 * max_edges;
    let len_alt_path = max_vertices + I_ALTP_HDR_LEN as i32 + 1;

    let mut p = BnStruct {
        edge: vec![BnsEdge::default(); max_edges as usize],
        vert: vec![BnsVertex::default(); max_vertices as usize],
        iedge: vec![0; max_iedges as usize],
        ..Default::default()
    };
    let mut num_altp = 0;
    while num_altp < max_altp && (num_altp as usize) < BN_MAX_ALTP {
        let mut ap = vec![BnsAltPath::default(); len_alt_path as usize];
        ap[I_ALTP_MAX_LEN].set_number(len_alt_path);
        p.len_alt_path = len_alt_path;
        ap[I_ALTP_FLOW].flow[0] = 0;
        ap[I_ALTP_START_ATOM].set_number(NO_VERTEX);
        ap[I_ALTP_END_ATOM].set_number(NO_VERTEX);
        ap[I_ALTP_PATH_LEN].set_number(0);
        p.altp.push(ap);
        num_altp += 1;
    }
    p.alt_path = None;
    p.num_altp = 0;
    p.max_altp = num_altp;

    p.vert[0].iedge = 0;
    for i in 0..num_atoms {
        let k = i32::from(at[i].valence) + (n_max_add_edges + NUM_KINDS_OF_GROUPS);
        p.vert[i].max_adj_edges = k as u16;
        let next = p.vert[i].iedge + k as usize;
        p.vert[i + 1].iedge = next;
    }
    p.num_atoms = num_atoms as i32;
    p.num_added_atoms = 0;
    p.num_t_groups = 0;
    p.num_c_groups = 0;
    p.n_max_add_atoms = n_max_add_atoms;
    p.n_max_add_edges = n_max_add_edges;
    p.num_vertices = num_atoms as i32;
    p.max_vertices = max_vertices;
    p.num_bonds = num_bonds;
    p.max_edges = max_edges;
    p.max_iedges = max_iedges;

    let (mut tot_st_cap, mut tot_st_flow) = (0, 0);
    let mut n_edges: usize = 0;
    for i in 0..num_atoms {
        let mut st_flow = 0;
        let valence = at[i].valence as usize;
        for j in 0..valence {
            let neigh = at[i].neighbor[j] as usize;
            let mut k = 0usize;
            while k < at[neigh].valence as usize {
                if at[neigh].neighbor[k] as usize == i {
                    break;
                }
                k += 1;
            }
            let mut bond_type = i32::from(at[i].bond_type[j]) & BOND_TYPE_MASK;
            let bond_mark = i32::from(at[i].bond_type[j]) & !BOND_TYPE_MASK;
            if bond_type != BOND_TYPE_SINGLE
                && bond_type != BOND_TYPE_DOUBLE
                && bond_type != BOND_TYPE_TRIPLE
            {
                bond_type = 1;
                at[i].bond_type[j] = (bond_mark | bond_type) as u8;
                num_changed_bonds += 1;
            }
            let edge_flow;
            if neigh > i {
                let f1 = max_at_flow(&at[i]);
                let f2 = max_at_flow(&at[neigh]);
                let mut ef = bond_type - 1;
                let edge_cap;
                if ef > MAX_BOND_EDGE_CAP {
                    ef = 0;
                    edge_cap = AROM_BOND_EDGE_CAP;
                } else {
                    edge_cap = f1.min(f2).min(MAX_BOND_EDGE_CAP);
                }
                edge_flow = ef;
                let e = &mut p.edge[n_edges];
                e.neighbor1 = i as u16;
                e.neighbor12 = (i ^ neigh) as u16;
                e.flow = edge_flow;
                e.flow0 = edge_flow;
                e.cap = edge_cap;
                e.cap0 = edge_cap;
                e.neigh_ord = [j as u16, k as u16];
                e.pass = 0;
                e.forbidden = 0;
                let o1 = p.vert[i].iedge + j;
                let o2 = p.vert[neigh].iedge + k;
                p.iedge[o1] = n_edges as i32;
                p.iedge[o2] = n_edges as i32;
                n_edges += 1;
            } else {
                let iedge = p.iedge[p.vert[neigh].iedge + k] as usize;
                edge_flow = p.edge[iedge].flow;
            }
            st_flow += edge_flow;
        }
        let v = &mut p.vert[i];
        v.num_adj_edges = valence as u16;
        let cap = max_at_flow(&at[i]);
        v.st_edge.cap = cap;
        v.st_edge.cap0 = cap;
        v.st_edge.flow = st_flow;
        v.st_edge.flow0 = st_flow;
        v.type_ = BNS_VERT_TYPE_ATOM;
        tot_st_cap += cap;
        tot_st_flow += st_flow;
    }
    p.num_edges = n_edges as i32;
    p.num_added_edges = 0;
    p.tot_st_cap = tot_st_cap;
    p.tot_st_flow = tot_st_flow;
    (p, num_changed_bonds / 2)
}

/// `ReInitBnStructAltPaths`
pub fn re_init_bn_struct_alt_paths(p: &mut BnStruct) -> i32 {
    let mut i = 0;
    while i < p.max_altp && (i as usize) < BN_MAX_ALTP {
        let ap = &mut p.altp[i as usize];
        ap[I_ALTP_FLOW].flow[0] = 0;
        ap[I_ALTP_PATH_LEN].set_number(0);
        ap[I_ALTP_START_ATOM].set_number(NO_VERTEX);
        ap[I_ALTP_END_ATOM].set_number(NO_VERTEX);
        i += 1;
    }
    p.alt_path = None;
    p.num_altp = 0;
    i
}

/// `ReInitBnStruct`: strips the group vertices and restores caps and flows
pub fn re_init_bn_struct(
    p: &mut BnStruct,
    at: &mut [InpAtom],
    num_at: usize,
    b_remove_groups_from_atoms: bool,
) -> i32 {
    let mut ret = 0;
    if p.vert.is_empty() || p.edge.is_empty() {
        ret += 2;
    } else {
        let i = p.edge[..p.num_edges as usize]
            .iter()
            .filter(|e| e.pass != 0)
            .count();
        ret += i as i32 * 100;
        for vfict in p.num_atoms..p.num_vertices {
            let nadj = p.vert[vfict as usize].num_adj_edges as usize;
            for kfict in 0..nadj {
                let iedgefict = p.iedge[p.vert[vfict as usize].iedge + kfict] as usize;
                let endpoint = (i32::from(p.edge[iedgefict].neighbor12)) ^ vfict;
                if b_remove_groups_from_atoms && (endpoint as usize) < num_at {
                    at[endpoint as usize].c_point = 0;
                    at[endpoint as usize].endpoint = 0;
                }
                let ne = p.vert[endpoint as usize].num_adj_edges as usize;
                for k in 0..ne {
                    let iedge = p.iedge[p.vert[endpoint as usize].iedge + k] as usize;
                    let centerpoint = (i32::from(p.edge[iedge].neighbor12) ^ endpoint) as usize;
                    p.edge[iedge].cap = p.edge[iedge].cap0;
                    p.edge[iedge].flow = p.edge[iedge].flow0;
                    p.edge[iedge].pass = 0;
                    p.vert[centerpoint].st_edge.cap = p.vert[centerpoint].st_edge.cap0;
                    p.vert[centerpoint].st_edge.flow = p.vert[centerpoint].st_edge.flow0;
                }
                let v = &mut p.vert[endpoint as usize];
                v.st_edge.cap = v.st_edge.cap0;
                v.st_edge.flow = v.st_edge.flow0;
                v.type_ &= BNS_VERT_TYPE_ATOM;
            }
        }
        if p.num_edges > p.num_bonds {
            for i in 0..p.num_atoms as usize {
                p.vert[i].num_adj_edges = (i32::from(p.vert[i].max_adj_edges)
                    - p.n_max_add_edges
                    - NUM_KINDS_OF_GROUPS) as u16;
            }
        }
    }
    if p.edge.is_empty() {
        ret += 4;
    }
    if p.iedge.is_empty() {
        ret += 8;
    }
    re_init_bn_struct_alt_paths(p);
    p.num_vertices = p.num_atoms;
    p.num_edges = p.num_bonds;
    p.num_added_atoms = 0;
    p.num_t_groups = 0;
    p.num_c_groups = 0;
    p.num_added_edges = 0;
    ret
}

/// `AllocateAndInitBnData`
pub fn allocate_and_init_bn_data(max_num_vertices: i32) -> BnData {
    let max_num_vertices = 2 * max_num_vertices + 2;
    let mut max_len_pu_pv = max_num_vertices / 2 + 1;
    max_len_pu_pv += max_len_pu_pv % 2;
    let n = max_num_vertices as usize;
    let m = max_len_pu_pv as usize;
    BnData {
        base_ptr: vec![NO_VERTEX; n],
        switch_edge: vec![[NO_VERTEX, 0]; n],
        tree: vec![TREE_NOT_IN_M; n],
        scan_q: vec![0; n],
        q_size: -1,
        pu: vec![0; m],
        pv: vec![0; m],
        max_num_vertices,
        max_len_pu_pv,
        rad_endpoints: vec![0; m],
        rad_edges: vec![0; m],
        n_num_rad_endpoints: 0,
        ..Default::default()
    }
}

/// `prim(v)`: the complementary vertex
pub fn prim(v: Vertex) -> Vertex {
    v ^ 1
}

/// `ReInitBnData`
pub fn re_init_bn_data(bd: &mut BnData) -> i32 {
    let mut ret = 0;
    if bd.scan_q.is_empty() {
        ret += 2;
    }
    if bd.base_ptr.is_empty() {
        ret += 4;
    }
    if bd.switch_edge.is_empty() {
        ret += 8;
    }
    if bd.tree.is_empty() {
        ret += 16;
    }
    if ret == 0 {
        for i in 0..=bd.q_size {
            let u = bd.scan_q[i as usize] as usize;
            let v = prim(u as Vertex) as usize;
            bd.base_ptr[u] = NO_VERTEX;
            bd.base_ptr[v] = NO_VERTEX;
            bd.switch_edge[u][0] = NO_VERTEX;
            bd.switch_edge[v][0] = NO_VERTEX;
            bd.tree[u] = TREE_NOT_IN_M;
            bd.tree[v] = TREE_NOT_IN_M;
        }
    }
    bd.q_size = -1;
    if bd.pu.is_empty() {
        ret += 32;
    }
    if bd.pv.is_empty() {
        ret += 64;
    }
    ret
}

/// `GetVertexDegree`
pub fn get_vertex_degree(p: &BnStruct, v: Vertex) -> i32 {
    let i = v / 2 - 1;
    if i >= 0 {
        if p.vert[i as usize].st_edge.cap > 0 {
            i32::from(p.vert[i as usize].num_adj_edges) + 1
        } else {
            0
        }
    } else {
        p.num_vertices
    }
}

/// `Get2ndEdgeVertex`
pub fn get_2nd_edge_vertex(p: &BnStruct, uv: Edge) -> Vertex {
    if uv[1] >= 0 {
        return ((uv[0] - 2) ^ (2 * i32::from(p.edge[uv[1] as usize].neighbor12) + 1)) + 2;
    }
    if uv[0] <= 1 { -(1 + uv[1]) } else { uv[0] % 2 }
}

/// `GetVertexNeighbor`: the `neigh`-th neighbor of `v` and the edge index
/// (`NO_VERTEX` if the edge has no capacity or is forbidden)
pub fn get_vertex_neighbor(p: &BnStruct, v: Vertex, neigh: i32, iedge: &mut EdgeIndex) -> Vertex {
    let i = v - 2;
    if i >= 0 {
        if neigh != 0 {
            let neigh = neigh - 1;
            *iedge = p.iedge[p.vert[(i / 2) as usize].iedge + neigh as usize];
            let e = &p.edge[*iedge as usize];
            if e.cap & EDGE_FLOW_MASK == 0 || e.forbidden != 0 {
                return NO_VERTEX;
            }
            (i ^ (2 * i32::from(e.neighbor12) + 1)) + 2
        } else {
            let neighbor = v & 1;
            *iedge = -(neighbor + 1);
            neighbor
        }
    } else {
        if p.vert[neigh as usize].st_edge.cap & EDGE_FLOW_ST_MASK == 0 {
            return NO_VERTEX;
        }
        let neighbor = 2 * neigh + 2 + (v & 1);
        *iedge = -(neighbor + 1);
        neighbor
    }
}

/// What `GetEdgePointer` points to
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeRef {
    None,
    /// an edge between atoms (index into `edge`)
    Edge(usize),
    /// the edge between an atom and s or t (the atom's vertex index)
    St(usize),
}

/// `GetEdgePointer`: returns `bBackward` (or an error) and sets the pointer
/// and `s_or_t` (0 = ordinary edge; 1, 2: u = s, t; 3, 4: v = s, t)
pub fn get_edge_pointer(
    _p: &BnStruct,
    u: Vertex,
    v: Vertex,
    iuv: EdgeIndex,
    uv: &mut EdgeRef,
    s_or_t: &mut i8,
) -> i32 {
    let i = u / 2 - 1;
    let j = v / 2 - 1;
    let mut b_backward = BNS_WRONG_PARMS;
    *uv = EdgeRef::None;
    if i >= 0 {
        if j >= 0 {
            if (u + v) % 2 != 0 {
                *uv = EdgeRef::Edge(iuv as usize);
                b_backward = u & 1;
                *s_or_t = 0;
            }
        } else if v >= 0 && (u + v) % 2 == 0 {
            *uv = EdgeRef::St(i as usize);
            b_backward = i32::from(v & 1 == 0);
            *s_or_t = (v + 3) as i8;
        }
    } else if j >= 0 && u >= 0 && (u + v) % 2 == 0 {
        *uv = EdgeRef::St(j as usize);
        b_backward = u & 1;
        *s_or_t = (u + 1) as i8;
    }
    b_backward
}

/// `AugmentEdge`: changes the flow of the edge u-v by `delta` (in the
/// direction `b_reverse`) and extends the alternating path; returns the new
/// flow or a negative error code
pub fn augment_edge(
    p: &mut BnStruct,
    u: Vertex,
    v: Vertex,
    iuv: EdgeIndex,
    delta: i32,
    b_reverse: bool,
    b_change_flow: i32,
) -> i32 {
    let mut pedge = EdgeRef::None;
    let mut s_or_t: i8 = 0;
    let b_backward = get_edge_pointer(p, u, v, iuv, &mut pedge, &mut s_or_t);
    if is_bns_error(b_backward) {
        return b_backward;
    }
    let delta = if b_backward != 0 { -delta } else { delta };
    let mut ret = 0;
    let f = if s_or_t != 0 {
        let EdgeRef::St(vi) = pedge else {
            unreachable!()
        };
        let st = p.vert[vi].st_edge;
        let flow = st.flow;
        let mut ff = (flow & EDGE_FLOW_ST_MASK) + delta;
        if delta == 0 {
            p.vert[vi].st_edge.flow = st.flow & !EDGE_FLOW_ST_PATH;
        } else {
            let cap = st.cap;
            if ff < 0 || ff > cap {
                ret = BNS_WRONG_PARMS;
            } else {
                if b_change_flow & BNS_EF_CHNG_FLOW == 0 {
                    ff -= delta;
                } else if delta != 0 {
                    p.vert[vi].st_edge.pass += 1;
                }
                let nf = (flow & !(EDGE_FLOW_ST_PATH | EDGE_FLOW_ST_MASK)) + ff;
                p.vert[vi].st_edge.flow = nf;
                if b_reverse {
                    match s_or_t {
                        1 => p.ap()[I_ALTP_END_ATOM].set_number(v / 2 - 1),
                        2 | 3 => ret = BNS_WRONG_PARMS,
                        4 => {
                            p.ap()[I_ALTP_START_ATOM].set_number(u / 2 - 1);
                            p.ap()[I_ALTP_FLOW].flow[0] = delta;
                        }
                        _ => ret = BNS_WRONG_PARMS,
                    }
                } else {
                    match s_or_t {
                        1 => {
                            p.ap()[I_ALTP_START_ATOM].set_number(v / 2 - 1);
                            p.ap()[I_ALTP_FLOW].flow[0] = delta;
                        }
                        2 | 3 => ret = BNS_WRONG_PARMS,
                        4 => p.ap()[I_ALTP_END_ATOM].set_number(u / 2 - 1),
                        _ => ret = BNS_WRONG_PARMS,
                    }
                }
            }
        }
        ff
    } else {
        let EdgeRef::Edge(ei) = pedge else {
            unreachable!()
        };
        let mut ff = (p.edge[ei].flow & EDGE_FLOW_MASK) + delta;
        if delta == 0 {
            p.edge[ei].flow &= !EDGE_FLOW_PATH;
        } else if ff < 0 || ff > p.edge[ei].cap {
            ret = BNS_WRONG_PARMS;
        } else {
            let iu = (u / 2 - 1) as u16;
            let iv = (v / 2 - 1) as u16;
            if b_change_flow & BNS_EF_CHNG_FLOW == 0 {
                ff -= delta;
            } else if delta != 0 {
                p.edge[ei].pass += 1;
            }
            p.edge[ei].flow = (p.edge[ei].flow & !(EDGE_FLOW_PATH | EDGE_FLOW_MASK)) | ff;
            let (path_len, max_len) = {
                let ap = p.ap();
                (ap[I_ALTP_PATH_LEN].number(), ap[I_ALTP_MAX_LEN].number())
            };
            // ALTP_MAY_ADD
            if I_ALTP_NEIGHBOR as i32 + path_len < max_len {
                let e = p.edge[ei];
                let indx = if b_reverse {
                    usize::from(e.neighbor1 == iv)
                } else {
                    usize::from(e.neighbor1 == iu)
                };
                let slot = I_ALTP_NEIGHBOR + path_len as usize;
                let ap = p.ap();
                ap[slot].set_ineigh(0, e.neigh_ord[1 - indx]);
                ap[slot].set_ineigh(1, e.neigh_ord[indx]);
                let n = ap[I_ALTP_PATH_LEN].number() + 1;
                ap[I_ALTP_PATH_LEN].set_number(n);
            } else {
                p.ap()[I_ALTP_FLOW].flow[1] = 1;
                ret = BNS_ALTPATH_OVFL;
            }
        }
        ff
    };
    if ret != 0 { ret } else { f }
}

/// `rescap_mark`: residual capacity of the edge; marks it as belonging to
/// the augmenting path
pub fn rescap_mark(p: &mut BnStruct, u: Vertex, v: Vertex, iuv: EdgeIndex) -> i32 {
    let mut pedge = EdgeRef::None;
    let mut s_or_t: i8 = 0;
    let b_backward = get_edge_pointer(p, u, v, iuv, &mut pedge, &mut s_or_t);
    if is_bns_error(b_backward) {
        return b_backward;
    }
    if s_or_t != 0 {
        let EdgeRef::St(vi) = pedge else {
            unreachable!()
        };
        let st = p.vert[vi].st_edge;
        let flow = st.flow;
        let mut f = flow & EDGE_FLOW_ST_MASK;
        if b_backward == 0 {
            f = st.cap - f;
        }
        if flow & EDGE_FLOW_ST_PATH != 0 {
            p.b_not_a_simple_path += 1;
            f /= 2;
        } else {
            p.vert[vi].st_edge.flow |= EDGE_FLOW_ST_PATH;
        }
        f
    } else {
        let EdgeRef::Edge(ei) = pedge else {
            unreachable!()
        };
        let flow = p.edge[ei].flow;
        let mut f = flow & EDGE_FLOW_MASK;
        if b_backward == 0 {
            f = p.edge[ei].cap - f;
        }
        if flow & EDGE_FLOW_PATH != 0 {
            f /= 2;
            p.b_not_a_simple_path += 1;
        } else {
            p.edge[ei].flow |= EDGE_FLOW_PATH;
        }
        f
    }
}

/// `rescap`: residual capacity of the edge u-v
pub fn rescap(p: &BnStruct, u: Vertex, v: Vertex, iuv: EdgeIndex) -> i32 {
    let mut pedge = EdgeRef::None;
    let mut s_or_t: i8 = 0;
    let b_backward = get_edge_pointer(p, u, v, iuv, &mut pedge, &mut s_or_t);
    if is_bns_error(b_backward) {
        return b_backward;
    }
    match pedge {
        EdgeRef::St(vi) => {
            let st = p.vert[vi].st_edge;
            let mut f = st.flow & EDGE_FLOW_ST_MASK;
            if b_backward == 0 {
                f = st.cap - f;
            }
            f
        }
        EdgeRef::Edge(ei) => {
            let e = p.edge[ei];
            let mut f = e.flow & EDGE_FLOW_MASK;
            if b_backward == 0 {
                f = e.cap - f;
            }
            f
        }
        EdgeRef::None => unreachable!(),
    }
}

/// `GetPrevVertex`: the previous vertex on the path found by the search
/// (`SwitchEdge` is the search table of [`BnData`])
pub fn get_prev_vertex(
    p: &BnStruct,
    y: Vertex,
    switch_edge: &[Edge],
    iuv: &mut EdgeIndex,
) -> Vertex {
    let mut w = switch_edge[y as usize][0];
    let mut z = get_2nd_edge_vertex(p, switch_edge[y as usize]);
    let mut iwy = switch_edge[y as usize][1];
    if z == y {
        *iuv = iwy;
        return w;
    }
    let x2 = prim(y);
    let mut y2 = prim(z);
    while y2 != NO_VERTEX {
        w = switch_edge[y2 as usize][0];
        z = get_2nd_edge_vertex(p, switch_edge[y2 as usize]);
        iwy = switch_edge[y2 as usize][1];
        if w == x2 {
            *iuv = iwy;
            return if (y + z) % 2 != 0 { z } else { prim(z) };
        }
        if w == y2 {
            return NO_VERTEX;
        }
        y2 = w;
    }
    y2
}

/// `bIgnoreVertexNonTACN_atom` (`FIX_TACN_POSSIBLE_BUG == 0`)
pub fn ignore_vertex_non_tacn_atom(p: &BnStruct, u: Vertex, v: Vertex) -> bool {
    if p.type_tacn == 0 || u <= 1 || v <= 1 || p.vert[(v / 2 - 1) as usize].type_ & p.type_tacn != 0
    {
        return false;
    }
    if p.type_t == 0 || p.type_cn == 0 {
        return false;
    }
    let ut = p.vert[(u / 2 - 1) as usize].type_;
    let u_is_taut = ut & p.type_t == p.type_t;
    if u_is_taut || ut & p.type_cn == p.type_cn {
        let degree = get_vertex_degree(p, v);
        let (mut num_allowed, mut num_found_groups) = (0, 0);
        for i in 0..degree {
            let mut ivw = 0;
            let w = get_vertex_neighbor(p, v, i, &mut ivw);
            if w == NO_VERTEX || w <= 1 {
                continue;
            }
            if w != u && rescap(p, v, w, ivw) > 0 {
                num_allowed += 1;
                let wt = p.vert[(w / 2 - 1) as usize].type_;
                let hit = if u_is_taut {
                    wt & p.type_cn == p.type_cn
                } else {
                    wt & p.type_t == p.type_t
                };
                if hit {
                    num_found_groups += 1;
                }
            }
        }
        if num_found_groups != 0 && num_allowed == 1 {
            return true;
        }
    }
    false
}

/// `bIgnoreVertexNonTACN_group` (`FIX_TACN_POSSIBLE_BUG == 0`)
pub fn ignore_vertex_non_tacn_group(
    p: &BnStruct,
    v: Vertex,
    w: Vertex,
    switch_edge: &[Edge],
) -> bool {
    if v <= 1 || w <= 1 {
        return false;
    }
    if p.type_tacn == 0 || p.vert[(v / 2 - 1) as usize].type_ & p.type_tacn != 0 {
        return false;
    }
    if p.type_t == 0 || p.type_cn == 0 {
        return false;
    }
    let mut iuv = 0;
    let u = get_prev_vertex(p, v, switch_edge, &mut iuv);
    if u == NO_VERTEX || iuv < 0 {
        return false;
    }
    let e = &p.edge[iuv as usize];
    if (i32::from(e.neighbor1) != u / 2 - 1 && i32::from(e.neighbor1) != v / 2 - 1)
        || ((i32::from(e.neighbor12) ^ (u / 2 - 1)) != v / 2 - 1)
    {
        return false;
    }
    let ut = p.vert[(u / 2 - 1) as usize].type_;
    let wt = p.vert[(w / 2 - 1) as usize].type_;
    let u_is_taut = ut & p.type_t == p.type_t;
    if u_is_taut || ut & p.type_cn == p.type_cn {
        let w_is_taut = wt & p.type_t == p.type_t;
        if (w_is_taut || wt & p.type_cn == p.type_cn)
            && i32::from(u_is_taut) + i32::from(w_is_taut) == 1
        {
            return true;
        }
    }
    false
}
