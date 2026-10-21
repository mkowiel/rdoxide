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

use super::ichierr::{
    BNS_ALTPATH_OVFL, BNS_BOND_ERR, BNS_PROGRAM_ERR, BNS_WRONG_PARMS, is_bns_error,
};
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
    /// `*pbTautFlags`
    pub pb_taut_flags: u32,
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

// ---- the balanced network search (Kocay & Stone) ----

use super::ichierr::BNS_VERT_EDGE_OVFL;

pub const TREE_IN_1: i8 = 3;
pub const RAD_SRCH_NORM: i32 = 0;
pub const RAD_SRCH_FROM_FICT: i32 = 1;

fn tree_mark(tree: &mut [i8], x: Vertex, mark: i8) {
    if tree[x as usize] < mark {
        tree[x as usize] = mark;
    }
}

fn tree_is_s_reachable(tree: &[i8], x: Vertex) -> bool {
    tree[x as usize] >= TREE_IN_2BLOSS
}

/// `FindBase`: the base of the blossom of `u` (with path compression)
pub fn find_base(u: Vertex, base_ptr: &mut [Vertex]) -> Vertex {
    let bp = base_ptr[u as usize];
    if bp == NO_VERTEX {
        NO_VERTEX
    } else if bp == BLOSSOM_BASE {
        u
    } else {
        let b = find_base(bp, base_ptr);
        base_ptr[u as usize] = b;
        b
    }
}

/// `FindPathToVertex_s`: the path of blossom bases from `x` to s; returns the
/// index of the last element
pub fn find_path_to_vertex_s(
    mut x: Vertex,
    switch_edge: &[Edge],
    base_ptr: &mut [Vertex],
    path: &mut [Vertex],
    max_path_len: i32,
) -> i32 {
    let mut i = 0i32;
    path[0] = x;
    while x != VERTEX_S {
        x = find_base(switch_edge[x as usize][0], base_ptr);
        i += 1;
        if i < max_path_len {
            path[i as usize] = x;
        } else {
            return BNS_WRONG_PARMS;
        }
    }
    i
}

/// `MakeBlossom`: merges the blossoms of `u` and `v`; returns the base `w` of
/// the new blossom or an error code
#[allow(clippy::too_many_arguments)]
pub fn make_blossom(
    p: &BnStruct,
    scan_q: &mut [Vertex],
    p_q_size: &mut i32,
    pu: &mut [Vertex],
    pv: &mut [Vertex],
    max_len_pu_pv: i32,
    switch_edge: &mut [Edge],
    base_ptr: &mut [Vertex],
    u: Vertex,
    v: Vertex,
    iuv: EdgeIndex,
    b_u: Vertex,
    b_v: Vertex,
    tree: &mut [i8],
) -> Vertex {
    let len_pu = find_path_to_vertex_s(b_u, switch_edge, base_ptr, pu, max_len_pu_pv);
    if is_bns_error(len_pu) {
        return len_pu;
    }
    let len_pv = find_path_to_vertex_s(b_v, switch_edge, base_ptr, pv, max_len_pu_pv);
    if is_bns_error(len_pv) {
        return len_pv;
    }
    let mut i = len_pu;
    let mut j = len_pv;
    while i >= 0 && j >= 0 && pu[i as usize] == pv[j as usize] {
        i -= 1;
        j -= 1;
    }
    i += 1;
    let mut w = pu[i as usize];
    let mut z = switch_edge[w as usize][0];
    let mut izw = switch_edge[w as usize][1];
    while w != VERTEX_S && rescap(p, z, w, izw) >= 2 {
        i += 1;
        w = pu[i as usize];
        z = switch_edge[w as usize][0];
        izw = switch_edge[w as usize][1];
    }
    let mut absorb = |z: Vertex, other: Vertex, switch_edge: &mut [Edge]| {
        base_ptr[z as usize] = w;
        base_ptr[prim(z) as usize] = w;
        if !tree_is_s_reachable(tree, prim(z)) {
            switch_edge[prim(z) as usize][0] = other;
            switch_edge[prim(z) as usize][1] = iuv;
            *p_q_size += 1;
            scan_q[*p_q_size as usize] = prim(z);
            tree_mark(tree, prim(z), TREE_IN_2BLOSS);
        }
    };
    let mut k = i - 1;
    while k >= 0 {
        absorb(pu[k as usize], prim(v), switch_edge);
        k -= 1;
    }
    let mut k = j;
    while k >= 0 {
        absorb(pv[k as usize], u, switch_edge);
        k -= 1;
    }
    if !tree_is_s_reachable(tree, prim(w)) {
        switch_edge[prim(w) as usize][0] = u;
        switch_edge[prim(w) as usize][1] = iuv;
        *p_q_size += 1;
        scan_q[*p_q_size as usize] = prim(w);
        tree_mark(tree, prim(w), TREE_IN_2BLOSS);
    }
    w
}

/// `PullFlow`: augments the flow by `delta` along the path between `x` and `y`
pub fn pull_flow(
    p: &mut BnStruct,
    switch_edge: &[Edge],
    x: Vertex,
    y: Vertex,
    delta: i32,
    b_reverse: bool,
    b_change_flow: i32,
) -> i32 {
    let w = switch_edge[y as usize][0];
    let z = get_2nd_edge_vertex(p, switch_edge[y as usize]);
    let iwz = switch_edge[y as usize][1];
    let mut ret = 0;
    if b_reverse {
        if z != y {
            ret = pull_flow(
                p,
                switch_edge,
                prim(y),
                prim(z),
                delta,
                !b_reverse,
                b_change_flow,
            );
        }
        if !is_bns_error(ret) {
            ret = augment_edge(p, w, z, iwz, delta, b_reverse, b_change_flow);
        }
        if w != x && !is_bns_error(ret) {
            ret = pull_flow(p, switch_edge, x, w, delta, b_reverse, b_change_flow);
        }
    } else {
        if w != x && !is_bns_error(ret) {
            ret = pull_flow(p, switch_edge, x, w, delta, b_reverse, b_change_flow);
        }
        if !is_bns_error(ret) {
            ret = augment_edge(p, w, z, iwz, delta, b_reverse, b_change_flow);
        }
        if z != y && !is_bns_error(ret) {
            ret = pull_flow(
                p,
                switch_edge,
                prim(y),
                prim(z),
                delta,
                !b_reverse,
                b_change_flow,
            );
        }
    }
    ret
}

/// `FindPathCap`: the minimum residual capacity along the path between x and y
pub fn find_path_cap(
    p: &mut BnStruct,
    switch_edge: &[Edge],
    x: Vertex,
    y: Vertex,
    mut delta: i32,
) -> i32 {
    let w = switch_edge[y as usize][0];
    let z = get_2nd_edge_vertex(p, switch_edge[y as usize]);
    let iwz = switch_edge[y as usize][1];
    let cap = rescap_mark(p, w, z, iwz);
    if is_bns_error(cap) {
        return cap;
    }
    if cap < delta {
        delta = cap;
    }
    if w != x {
        let delta2 = find_path_cap(p, switch_edge, x, w, delta);
        delta = delta2.min(delta);
    }
    if z != y {
        let delta2 = find_path_cap(p, switch_edge, prim(y), prim(z), delta);
        delta = delta2.min(delta);
    }
    delta
}

/// `bRadChangesAtomType`
pub fn rad_changes_atom_type(
    p: &BnStruct,
    bd: &BnData,
    v: Vertex,
    mut v_1: Vertex,
    mut v_2: Vertex,
) -> bool {
    let mut iuv = 0;
    if v_1 == NO_VERTEX {
        v_1 = get_prev_vertex(p, v, &bd.switch_edge, &mut iuv);
    }
    let v_o = v_1 / 2 - 1;
    if v_o < 0 || v_o >= p.num_atoms {
        return false;
    }
    let second = p.iedge[p.vert[v_o as usize].iedge + 1] as usize;
    let nb = (i32::from(p.edge[second].neighbor12) ^ v_o) as usize;
    if p.vert[nb].type_ & BNS_VERT_TYPE_ATOM != 0 {
        return false;
    }
    if v_2 == NO_VERTEX {
        v_2 = get_prev_vertex(p, v_1, &bd.switch_edge, &mut iuv);
    }
    let v_chg_or_h = v_2 / 2 - 1;
    if v_chg_or_h < p.num_atoms {
        return false;
    }
    p.vert[v_chg_or_h as usize].type_ & (BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP) != 0
}

/// `RegisterRadEndpoint`: registers a (radical, endpoint) pair found by the search
pub fn register_rad_endpoint(p: &BnStruct, bd: &mut BnData, mut u: Vertex) -> i32 {
    let mut iuv = 0;
    match bd.b_rad_srch_mode {
        RAD_SRCH_NORM => {
            let mut num_found = 0;
            while u > VERTEX_T && (u % 2 != 0 || u / 2 > p.num_atoms) {
                u = get_prev_vertex(p, u, &bd.switch_edge, &mut iuv);
            }
            let w = u / 2 - 1;
            if VERTEX_T < u
                && w < p.num_atoms
                && p.vert[w as usize].st_edge.cap
                    == (p.vert[w as usize].st_edge.flow & EDGE_FLOW_ST_MASK)
            {
                let mut v = u;
                let u_last = u;
                while v > VERTEX_T {
                    u = v;
                    v = get_prev_vertex(p, u, &bd.switch_edge, &mut iuv);
                }
                if u % 2 == 0 && VERTEX_T < u {
                    u = u / 2 - 1;
                    if u < p.num_atoms
                        && p.vert[u as usize].st_edge.cap
                            > (p.vert[u as usize].st_edge.flow & EDGE_FLOW_ST_MASK)
                    {
                        let mut v = u_last;
                        while v > VERTEX_T {
                            let v2 = v / 2 - 1;
                            if v % 2 == 0
                                && v2 < p.num_atoms
                                && p.vert[v2 as usize].st_edge.cap
                                    == (p.vert[v2 as usize].st_edge.flow & EDGE_FLOW_ST_MASK)
                            {
                                if p.type_tacn != 0
                                    && rad_changes_atom_type(p, bd, v, NO_VERTEX, NO_VERTEX)
                                {
                                    v = get_prev_vertex(p, v, &bd.switch_edge, &mut iuv);
                                    continue;
                                }
                                let mut i = 0;
                                while i < bd.n_num_rad_endpoints {
                                    if u == bd.rad_endpoints[i as usize]
                                        && v2 == bd.rad_endpoints[i as usize + 1]
                                    {
                                        break;
                                    }
                                    i += 2;
                                }
                                if i >= bd.n_num_rad_endpoints {
                                    if bd.n_num_rad_endpoints + 2 <= bd.max_num_vertices {
                                        bd.rad_endpoints[bd.n_num_rad_endpoints as usize] = u;
                                        bd.n_num_rad_endpoints += 1;
                                        bd.rad_endpoints[bd.n_num_rad_endpoints as usize] = v2;
                                        bd.n_num_rad_endpoints += 1;
                                        num_found += 1;
                                    } else {
                                        return BNS_VERT_EDGE_OVFL;
                                    }
                                }
                            }
                            v = get_prev_vertex(p, v, &bd.switch_edge, &mut iuv);
                        }
                        if num_found != 0 {
                            return 1;
                        }
                    }
                }
            }
        }
        RAD_SRCH_FROM_FICT => {
            let mut v = u;
            let mut w = NO_VERTEX;
            u = NO_VERTEX;
            while v > VERTEX_T {
                u = v;
                if v % 2 == 0 && v / 2 <= p.num_atoms {
                    let st = p.vert[(v / 2 - 1) as usize].st_edge;
                    if st.cap - st.flow < 2 {
                        w = v;
                    }
                }
                v = get_prev_vertex(p, u, &bd.switch_edge, &mut iuv);
            }
            v = u / 2 - 1;
            if w == NO_VERTEX
                || u == NO_VERTEX
                || w % 2 != 0
                || u == w
                || v < p.num_atoms
                || p.vert[v as usize].st_edge.cap == p.vert[v as usize].st_edge.flow
            {
                return 0;
            }
            w = w / 2 - 1;
            if w >= p.num_atoms {
                return 0;
            }
            u = v;
            let mut i = 0;
            while i < bd.n_num_rad_endpoints {
                if u == bd.rad_endpoints[i as usize] && w == bd.rad_endpoints[i as usize + 1] {
                    return 0;
                }
                i += 2;
            }
            if bd.n_num_rad_endpoints + 2 <= bd.max_num_vertices {
                bd.rad_endpoints[bd.n_num_rad_endpoints as usize] = u;
                bd.n_num_rad_endpoints += 1;
                bd.rad_endpoints[bd.n_num_rad_endpoints as usize] = w;
                bd.n_num_rad_endpoints += 1;
                return 1;
            }
            return BNS_VERT_EDGE_OVFL;
        }
        _ => {}
    }
    0
}

/// `BalancedNetworkSearch`: searches for an augmenting path and, if one is
/// found, augments the flow along it. Returns the augmentation (>0), 0 or a
/// negative error code.
pub fn balanced_network_search(p: &mut BnStruct, bd: &mut BnData, b_change_flow: i32) -> i32 {
    let mut q_size;
    let mut k: i32 = 0;
    let mut ret = 0;
    let b_rad_search = (BNS_EF_RAD_SRCH & b_change_flow) != 0 && !bd.rad_endpoints.is_empty();
    let mut b_rad_srch_mode = RAD_SRCH_NORM;
    let mut b_rad_search_prelim = false;
    if b_rad_search {
        bd.n_num_rad_endpoints = 0;
        b_rad_srch_mode = bd.b_rad_srch_mode;
        b_rad_search_prelim = p.type_tacn != 0 && b_rad_srch_mode == RAD_SRCH_NORM;
    }

    q_size = 0;
    bd.scan_q[0] = VERTEX_S;
    bd.base_ptr[VERTEX_T as usize] = VERTEX_S;
    bd.base_ptr[VERTEX_S as usize] = BLOSSOM_BASE;
    bd.tree[VERTEX_S as usize] = TREE_IN_1;

    loop {
        let u = bd.scan_q[k as usize];
        let mut b_u = find_base(u, &mut bd.base_ptr);
        let degree = get_vertex_degree(p, u);
        let mut n = 0;
        for i in 0..degree {
            let mut iuv: EdgeIndex = 0;
            let v = get_vertex_neighbor(p, u, i, &mut iuv);
            if v == NO_VERTEX {
                continue;
            }
            if k == 0 && b_rad_srch_mode == RAD_SRCH_FROM_FICT && v / 2 <= p.num_atoms {
                continue;
            }
            if b_rad_search_prelim && v / 2 > p.num_atoms {
                continue;
            }
            let sw_u = bd.switch_edge[u as usize];
            if (sw_u[0] != v || get_2nd_edge_vertex(p, sw_u) != u) && {
                ret = rescap(p, u, v, iuv);
                ret > 0
            } {
                if p.type_tacn != 0 {
                    if ignore_vertex_non_tacn_atom(p, u, v) {
                        continue;
                    }
                    if ignore_vertex_non_tacn_group(p, u, v, &bd.switch_edge) {
                        continue;
                    }
                }
                let b_v = find_base(v, &mut bd.base_ptr);
                if b_v == NO_VERTEX {
                    q_size += 1;
                    bd.scan_q[q_size as usize] = v;
                    tree_mark(&mut bd.tree, v, TREE_IN_1);
                    tree_mark(&mut bd.tree, prim(v), TREE_IN_2);
                    bd.switch_edge[v as usize][0] = u;
                    bd.switch_edge[v as usize][1] = iuv;
                    bd.base_ptr[prim(v) as usize] = v;
                    bd.base_ptr[v as usize] = BLOSSOM_BASE;
                    n += 1;
                } else {
                    let sw_pu = bd.switch_edge[prim(u) as usize];
                    if tree_is_s_reachable(&bd.tree, prim(v))
                        && (sw_pu[0] != prim(v) || get_2nd_edge_vertex(p, sw_pu) != prim(u))
                        && b_u != b_v
                        && !(p.type_tacn != 0
                            && ignore_vertex_non_tacn_group(p, prim(v), u, &bd.switch_edge))
                    {
                        n += 1;
                        let w = make_blossom(
                            p,
                            &mut bd.scan_q,
                            &mut q_size,
                            &mut bd.pu,
                            &mut bd.pv,
                            bd.max_len_pu_pv,
                            &mut bd.switch_edge,
                            &mut bd.base_ptr,
                            u,
                            v,
                            iuv,
                            b_u,
                            b_v,
                            &mut bd.tree,
                        );
                        if is_bns_error(w) {
                            bd.q_size = q_size;
                            return w;
                        }
                        b_u = w;
                        if prim(w) == VERTEX_T {
                            let sw = bd.switch_edge.clone();
                            let delta = find_path_cap(p, &sw, VERTEX_S, VERTEX_T, 10000);
                            if is_bns_error(delta) {
                                bd.q_size = q_size;
                                return delta;
                            }
                            if delta != 0 {
                                p.b_change_flow |= b_change_flow & BNS_EF_CHNG_FLOW;
                            }
                            let r =
                                pull_flow(p, &sw, VERTEX_S, VERTEX_T, delta, false, b_change_flow);
                            bd.q_size = q_size;
                            return if is_bns_error(r) { r } else { delta };
                        }
                    }
                }
            } else if is_bns_error(ret) {
                bd.q_size = q_size;
                return ret;
            }
        }
        if b_rad_search && n == 0 {
            let r = register_rad_endpoint(p, bd, u);
            if is_bns_error(r) {
                bd.q_size = q_size;
                return r;
            }
        }
        k += 1;
        if k > q_size {
            break;
        }
    }
    bd.q_size = q_size;
    0
}

/// `RunBalancedNetworkSearch`: runs the search until no augmenting path is
/// found (the timeout check of the C code is not ported: no time limit is set)
pub fn run_balanced_network_search(p: &mut BnStruct, bd: &mut BnData, b_change_flow: i32) -> i32 {
    let mut delta = 0;
    let mut n_sum_delta = 0;
    let mut pass = 0;
    while pass < p.max_altp {
        p.alt_path = Some(pass as usize);
        p.b_change_flow = 0;
        delta = balanced_network_search(p, bd, b_change_flow);
        re_init_bn_data(bd);
        if delta > 0 {
            p.num_altp += 1;
            n_sum_delta += delta.abs();
        } else {
            break;
        }
        pass += 1;
    }
    if is_bns_error(delta) {
        return delta;
    }
    n_sum_delta
}

/// `bIsBnsEndpoint`: is atom `v` an endpoint connected to a t-group by an edge
/// that is not forbidden?
pub fn is_bns_endpoint(p: &BnStruct, v: i32) -> bool {
    if v >= 0 && v < p.num_atoms && p.vert[v as usize].type_ & BNS_VERT_TYPE_ENDPOINT != 0 {
        let vert = &p.vert[v as usize];
        for i in (0..vert.num_adj_edges as usize).rev() {
            let e = &p.edge[p.iedge[vert.iedge + i] as usize];
            let vt = (i32::from(e.neighbor12) ^ v) as usize;
            if p.vert[vt].type_ & BNS_VERT_TYPE_TGROUP != 0 {
                return e.forbidden == 0;
            }
        }
    }
    false
}

pub const BNS_VERT_TYPE_ANY_GROUP: u16 =
    BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_SUPER_TGROUP;

/// `bAddNewVertex`: adds a temporary vertex joined to `v2` by a new edge;
/// returns the new vertex or an error code
pub fn add_new_vertex(
    p: &mut BnStruct,
    v2: Vertex,
    cap: i32,
    flow: i32,
    max_adj_edges: i32,
    n_dots: &mut i32,
) -> i32 {
    let vlast = (p.num_vertices - 1) as usize;
    let vnew = p.num_vertices;
    let iedge = p.num_edges;
    if iedge >= p.max_edges || vnew >= p.max_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    if p.vert[vlast].iedge as i32 + i32::from(p.vert[vlast].max_adj_edges) + max_adj_edges
        >= p.max_iedges
    {
        return BNS_VERT_EDGE_OVFL;
    }
    let old = p.vert[v2 as usize];
    if old.num_adj_edges >= old.max_adj_edges || max_adj_edges <= 0 {
        return BNS_VERT_EDGE_OVFL;
    }
    let new_iedge = p.vert[vlast].iedge + usize::from(p.vert[vlast].max_adj_edges);
    let mut e = BnsEdge {
        cap,
        cap0: cap,
        flow,
        flow0: flow,
        pass: 0,
        neighbor1: v2 as u16,
        neighbor12: (v2 ^ vnew) as u16,
        forbidden: 0,
        ..Default::default()
    };
    e.neigh_ord[usize::from(v2 > vnew)] = old.num_adj_edges;
    e.neigh_ord[usize::from(v2 < vnew)] = 0;
    p.edge[iedge as usize] = e;
    p.vert[vnew as usize] = BnsVertex {
        max_adj_edges: max_adj_edges as u16,
        num_adj_edges: 1,
        st_edge: BnsStEdge {
            cap,
            cap0: cap,
            flow,
            flow0: flow,
            pass: 0,
        },
        iedge: new_iedge,
        type_: BNS_VERT_TYPE_TEMP,
    };
    p.iedge[new_iedge] = iedge;
    *n_dots += cap - flow;
    let o = p.vert[v2 as usize].iedge + usize::from(old.num_adj_edges);
    p.iedge[o] = iedge;
    let pv = &mut p.vert[v2 as usize];
    pv.num_adj_edges += 1;
    *n_dots -= pv.st_edge.cap - pv.st_edge.flow;
    pv.st_edge.flow += flow;
    if pv.st_edge.cap < pv.st_edge.flow {
        pv.st_edge.cap = pv.st_edge.flow;
    }
    *n_dots += pv.st_edge.cap - pv.st_edge.flow;
    p.num_edges += 1;
    p.num_vertices += 1;
    vnew
}

/// `AddNewEdge`: connects the vertices `ip1` and `ip2`; returns the edge index
/// or an error code
pub fn add_new_edge(p: &mut BnStruct, ip1: i32, ip2: i32, cap: i32, flow: i32) -> i32 {
    let ie = p.num_edges;
    if ip1 >= p.max_vertices || ip1 < 0 || ip2 >= p.max_vertices || ip2 < 0 || ie >= p.max_edges {
        return BNS_VERT_EDGE_OVFL;
    }
    let (a, b) = (p.vert[ip1 as usize], p.vert[ip2 as usize]);
    if a.iedge as i32 + i32::from(a.max_adj_edges) > p.max_iedges
        || b.iedge as i32 + i32::from(b.max_adj_edges) > p.max_iedges
        || a.num_adj_edges >= a.max_adj_edges
        || b.num_adj_edges >= b.max_adj_edges
    {
        return BNS_VERT_EDGE_OVFL;
    }
    let mut e = BnsEdge {
        neighbor1: ip1.min(ip2) as u16,
        neighbor12: (ip1 ^ ip2) as u16,
        cap,
        cap0: cap,
        flow,
        flow0: flow,
        ..Default::default()
    };
    p.iedge[a.iedge + usize::from(a.num_adj_edges)] = ie;
    p.iedge[b.iedge + usize::from(b.num_adj_edges)] = ie;
    e.neigh_ord[usize::from(ip1 > ip2)] = a.num_adj_edges;
    e.neigh_ord[usize::from(ip1 < ip2)] = b.num_adj_edges;
    p.edge[ie as usize] = e;
    for i in [ip1, ip2] {
        let v = &mut p.vert[i as usize];
        v.num_adj_edges += 1;
        v.st_edge.flow += flow;
        if v.st_edge.cap < v.st_edge.flow {
            v.st_edge.cap = v.st_edge.flow;
        }
    }
    p.num_edges += 1;
    ie
}

/// `GetEdgeToGroupVertex`
pub fn get_edge_to_group_vertex(p: &BnStruct, v1: Vertex, type_: u16) -> i32 {
    if v1 < p.num_atoms {
        let pv = &p.vert[v1 as usize];
        for i in (0..pv.num_adj_edges as usize).rev() {
            let ie = p.iedge[pv.iedge + i];
            let e = &p.edge[ie as usize];
            let v2 = i32::from(e.neighbor12) ^ v1;
            if p.vert[v2 as usize].type_ == type_ {
                return if e.forbidden != 0 { NO_VERTEX } else { ie };
            }
        }
        return NO_VERTEX;
    }
    if v1 < p.num_vertices {
        return NO_VERTEX;
    }
    BNS_VERT_EDGE_OVFL
}

/// `GetGroupVertex`
pub fn get_group_vertex(p: &BnStruct, v1: Vertex, type_: u16) -> Vertex {
    if v1 < p.num_atoms {
        let pv = &p.vert[v1 as usize];
        let type2 = if type_ == BNS_VERT_TYPE_ENDPOINT {
            BNS_VERT_TYPE_TGROUP
        } else if type_ == BNS_VERT_TYPE_C_POINT {
            BNS_VERT_TYPE_C_GROUP
        } else {
            0
        };
        if pv.type_ & type_ == type_ {
            for i in (0..pv.num_adj_edges as usize).rev() {
                let e = &p.edge[p.iedge[pv.iedge + i] as usize];
                let v2 = i32::from(e.neighbor12) ^ v1;
                if p.vert[v2 as usize].type_ == type2 {
                    return if e.forbidden != 0 { NO_VERTEX } else { v2 };
                }
            }
        }
        return super::ichierr::BNS_BOND_ERR;
    }
    if v1 < p.num_vertices {
        return NO_VERTEX;
    }
    BNS_VERT_EDGE_OVFL
}

/// `bAddStCapToAVertex`: increments the st-cap of `v1`, saving the old caps in
/// `old_caps`; returns the number of saved values
pub fn add_st_cap_to_a_vertex(
    p: &mut BnStruct,
    v1: Vertex,
    v2: Vertex,
    old_caps: &mut [i32],
    n_dots: &mut i32,
    b_adjacent_donors: bool,
) -> i32 {
    let mut n = 0usize;
    old_caps[n] = p.vert[v1 as usize].st_edge.cap;
    n += 1;
    p.vert[v1 as usize].st_edge.cap += 1;
    *n_dots += 1;
    let v1v = p.vert[v1 as usize];
    if v1v.type_ & BNS_VERT_TYPE_ANY_GROUP == 0 {
        for i in 0..v1v.num_adj_edges as usize {
            let ie = p.iedge[v1v.iedge + i] as usize;
            old_caps[n] = p.edge[ie].cap;
            n += 1;
            let v = i32::from(p.edge[ie].neighbor12) ^ v1;
            if v == v2 && !b_adjacent_donors {
                continue;
            }
            let pv = p.vert[v as usize];
            if pv.type_ & BNS_VERT_TYPE_ANY_GROUP != 0 {
                continue;
            }
            let new_cap = pv.st_edge.cap.min(v1v.st_edge.cap).min(MAX_BOND_EDGE_CAP);
            p.edge[ie].cap = new_cap;
        }
    }
    n as i32
}

// ---- temporary changes of the network made to check an alternating path ----

pub const ALT_PATH_MODE_TAUTOM: i32 = 1;
pub const ALT_PATH_MODE_CHARGE: i32 = 2;
pub const ALT_PATH_MODE_4_SALT: i32 = 3;
pub const ALT_PATH_MODE_4_SALT2: i32 = 4;
pub const ALT_PATH_MODE_REM2H_CHG: i32 = 5;
pub const ALT_PATH_MODE_ADD2H_CHG: i32 = 6;
pub const ALT_PATH_MODE_REM2H_TST: i32 = 7;
pub const ALT_PATH_MODE_ADD2H_TST: i32 = 8;
pub const ALT_PATH_MODE_REM_PROTON: i32 = 9;

pub const BNS_CHK_ALTP_NO_ALTPATH: i32 = 0;
pub const BNS_CHK_ALTP_SAME_TGROUP: i32 = 1;
pub const BNS_CHK_ALTP_SAME_VERTEX: i32 = 2;
pub const BNS_CHK_ALTP_SET_SUCCESS: i32 = 4;

pub const TG_FLAG_TEST_TAUT2_SALTS: u32 = 0x10;

/// `ALT_PATH_CHANGES`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AltPathChanges {
    pub n_old_caps_vert: [[i32; super::inpdef::MAXVAL + 1]; 2],
    pub v_old_vert: [Vertex; 2],
    pub b_set_old_caps_vert: [i32; 2],
    pub v_new_vertex: [Vertex; 2],
    pub b_set_new: [i32; 2],
}

/// `BNS_FLOW_CHANGES`; a list of them ends with `iedge == NO_VERTEX`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BnsFlowChanges {
    pub iedge: i32,
    pub flow: i32,
    pub cap: i32,
    pub v1: Vertex,
    pub cap_st1: i32,
    pub flow_st1: i32,
    pub v2: Vertex,
    pub cap_st2: i32,
    pub flow_st2: i32,
}

impl Default for BnsFlowChanges {
    fn default() -> Self {
        Self {
            iedge: NO_VERTEX,
            flow: 0,
            cap: 0,
            v1: NO_VERTEX,
            cap_st1: 0,
            flow_st1: 0,
            v2: NO_VERTEX,
            cap_st2: 0,
            flow_st2: 0,
        }
    }
}

fn st_dec(x: i32, d: i32) -> i32 {
    ((x & EDGE_FLOW_ST_MASK) - d) | (x & !EDGE_FLOW_ST_MASK)
}

/// `bRestoreFlowAfterCheckOneBond`
pub fn restore_flow_after_check_one_bond(p: &mut BnStruct, fcd: &[BnsFlowChanges]) -> i32 {
    let mut ifcd = 0usize;
    while fcd[ifcd].iedge != NO_VERTEX {
        ifcd += 1;
    }
    while ifcd > 0 {
        ifcd -= 1;
        let f = fcd[ifcd];
        let e = &mut p.edge[f.iedge as usize];
        e.flow = f.flow;
        e.cap = f.cap;
        e.pass = 0;
        if f.v1 != NO_VERTEX {
            let v = &mut p.vert[f.v1 as usize];
            v.st_edge.flow = f.flow_st1;
            v.st_edge.cap = f.cap_st1;
            v.st_edge.pass = 0;
        }
        if f.v2 != NO_VERTEX {
            let v = &mut p.vert[f.v2 as usize];
            v.st_edge.flow = f.flow_st2;
            v.st_edge.cap = f.cap_st2;
            v.st_edge.pass = 0;
        }
    }
    0
}

/// `bSetFlowToCheckOneBond`: forces the flow of an edge to `flow`, recording
/// the changes in `fcd`; returns the number of dots created or an error
pub fn set_flow_to_check_one_bond(
    p: &mut BnStruct,
    iedge: i32,
    flow: i32,
    fcd: &mut [BnsFlowChanges],
) -> i32 {
    let f12 = p.edge[iedge as usize].flow & EDGE_FLOW_MASK;
    let mut ifcd = 0usize;
    let mut n_dots = 0;
    fcd[ifcd].iedge = NO_VERTEX;
    if f12 < flow {
        let v1 = i32::from(p.edge[iedge as usize].neighbor1);
        let v2 = v1 ^ i32::from(p.edge[iedge as usize].neighbor12);
        if (p.vert[v1 as usize].st_edge.cap & EDGE_FLOW_ST_MASK) < flow
            || (p.vert[v2 as usize].st_edge.cap & EDGE_FLOW_ST_MASK) < flow
        {
            return super::ichierr::BNS_CANT_SET_BOND;
        }
        if (p.vert[v1 as usize].st_edge.flow & EDGE_FLOW_ST_MASK) < f12
            || (p.vert[v2 as usize].st_edge.flow & EDGE_FLOW_ST_MASK) < f12
        {
            return super::ichierr::BNS_CAP_FLOW_ERR;
        }
        {
            let e = p.edge[iedge as usize];
            fcd[ifcd] = BnsFlowChanges {
                iedge,
                flow: e.flow,
                cap: e.cap,
                v1,
                flow_st1: p.vert[v1 as usize].st_edge.flow,
                cap_st1: p.vert[v1 as usize].st_edge.cap,
                v2,
                flow_st2: p.vert[v2 as usize].st_edge.flow,
                cap_st2: p.vert[v2 as usize].st_edge.cap,
            };
            ifcd += 1;
            fcd[ifcd].iedge = NO_VERTEX;
            p.edge[iedge as usize].pass |= 64;
        }
        let mut delta1 = flow - f12;
        let mut delta2 = delta1;
        if f12 > 0 {
            for v in [v1, v2] {
                let st = &mut p.vert[v as usize].st_edge;
                st.cap = st_dec(st.cap, f12);
                st.flow = st_dec(st.flow, f12);
            }
            let e = &mut p.edge[iedge as usize];
            e.flow &= !EDGE_FLOW_MASK;
        }
        p.edge[iedge as usize].cap &= !EDGE_FLOW_MASK;
        for (v, delta) in [(v1, &mut delta1), (v2, &mut delta2)] {
            let st = &mut p.vert[v as usize].st_edge;
            let mut rescap = (st.cap & EDGE_FLOW_ST_MASK) - (st.flow & EDGE_FLOW_ST_MASK);
            while rescap != 0 && *delta != 0 {
                rescap -= 1;
                *delta -= 1;
                st.cap = st_dec(st.cap, 1);
                n_dots -= 1;
            }
        }
        for (v, delta) in [(v1, &mut delta1), (v2, &mut delta2)] {
            let nadj = p.vert[v as usize].num_adj_edges as usize;
            let mut i = 0;
            while *delta != 0 && i < nadj {
                let iedge_i = p.iedge[p.vert[v as usize].iedge + i];
                i += 1;
                if iedge_i == iedge {
                    continue;
                }
                if p.edge[iedge_i as usize].forbidden != 0 {
                    continue;
                }
                let mut f = p.edge[iedge_i as usize].flow & EDGE_FLOW_MASK;
                if f != 0 {
                    let v_i = i32::from(p.edge[iedge_i as usize].neighbor12) ^ v;
                    fcd[ifcd] = BnsFlowChanges {
                        iedge: iedge_i,
                        flow: p.edge[iedge_i as usize].flow,
                        cap: p.edge[iedge_i as usize].cap,
                        v1: v_i,
                        flow_st1: p.vert[v_i as usize].st_edge.flow,
                        cap_st1: p.vert[v_i as usize].st_edge.cap,
                        v2: NO_VERTEX,
                        flow_st2: 0,
                        cap_st2: 0,
                    };
                    ifcd += 1;
                    fcd[ifcd].iedge = NO_VERTEX;
                    p.edge[iedge_i as usize].pass |= 64;
                    while f != 0 && *delta != 0 {
                        f -= 1;
                        *delta -= 1;
                        let ef = p.edge[iedge_i as usize].flow;
                        p.edge[iedge_i as usize].flow =
                            ((ef & EDGE_FLOW_MASK) - 1) | (ef & !EDGE_FLOW_MASK);
                        let sv = &mut p.vert[v_i as usize].st_edge;
                        sv.flow = st_dec(sv.flow, 1);
                        let s0 = &mut p.vert[v as usize].st_edge;
                        s0.cap = st_dec(s0.cap, 1);
                        s0.flow = st_dec(s0.flow, 1);
                        n_dots += 1;
                    }
                }
            }
        }
        if delta1 != 0 || delta2 != 0 {
            return super::ichierr::BNS_CANT_SET_BOND;
        }
    }
    if f12 >= flow {
        let v1 = i32::from(p.edge[iedge as usize].neighbor1);
        let v2 = v1 ^ i32::from(p.edge[iedge as usize].neighbor12);
        if (p.vert[v1 as usize].st_edge.flow & EDGE_FLOW_ST_MASK) < f12
            || (p.vert[v2 as usize].st_edge.flow & EDGE_FLOW_ST_MASK) < f12
            || (p.vert[v1 as usize].st_edge.cap & EDGE_FLOW_ST_MASK) < flow
            || (p.vert[v2 as usize].st_edge.cap & EDGE_FLOW_ST_MASK) < flow
        {
            return super::ichierr::BNS_CAP_FLOW_ERR;
        }
        let e = p.edge[iedge as usize];
        fcd[ifcd] = BnsFlowChanges {
            iedge,
            flow: e.flow,
            cap: e.cap,
            v1,
            flow_st1: p.vert[v1 as usize].st_edge.flow,
            cap_st1: p.vert[v1 as usize].st_edge.cap,
            v2,
            flow_st2: p.vert[v2 as usize].st_edge.flow,
            cap_st2: p.vert[v2 as usize].st_edge.cap,
        };
        ifcd += 1;
        fcd[ifcd].iedge = NO_VERTEX;
        p.edge[iedge as usize].pass |= 64;
        let delta = f12 - flow;
        for v in [v1, v2] {
            let st = &mut p.vert[v as usize].st_edge;
            st.flow = st_dec(st.flow, f12);
            st.cap = st_dec(st.cap, flow);
        }
        let e = &mut p.edge[iedge as usize];
        e.flow &= !EDGE_FLOW_MASK;
        e.cap &= !EDGE_FLOW_MASK;
        n_dots = 2 * delta;
    }
    n_dots
}

/// `bRestoreBnsAfterCheckAltPath`
pub fn restore_bns_after_check_alt_path(
    p: &mut BnStruct,
    apc: &AltPathChanges,
    b_change_flow: i32,
) -> i32 {
    fn remove_new(p: &mut BnStruct, apc: &AltPathChanges, with_flow: bool) {
        for i in (0..2).rev() {
            if apc.b_set_new[i] != 0 {
                let v_new = apc.v_new_vertex[i];
                let nv = p.vert[v_new as usize];
                for j in 0..nv.num_adj_edges as usize {
                    let ie = p.iedge[nv.iedge + j] as usize;
                    let e = p.edge[ie];
                    let v_old = (i32::from(e.neighbor12) ^ v_new) as usize;
                    if with_flow {
                        p.vert[v_old].st_edge.flow -= e.flow;
                        p.vert[v_old].st_edge.cap -= e.flow;
                    }
                    p.vert[v_old].num_adj_edges -= 1;
                    let at = p.vert[v_old].iedge + usize::from(p.vert[v_old].num_adj_edges);
                    p.iedge[at] = 0;
                    p.edge[ie] = BnsEdge::default();
                    p.num_edges -= 1;
                }
                p.vert[v_new as usize] = BnsVertex::default();
                p.num_vertices -= 1;
            }
        }
    }
    fn restore_caps(p: &mut BnStruct, apc: &AltPathChanges, check_flow: bool) {
        for i in (0..2).rev() {
            let mut n = apc.b_set_old_caps_vert[i];
            if n != 0 {
                let ov = p.vert[apc.v_old_vert[i] as usize];
                if !check_flow || ov.st_edge.flow <= apc.n_old_caps_vert[i][0] {
                    p.vert[apc.v_old_vert[i] as usize].st_edge.cap = apc.n_old_caps_vert[i][0];
                    n -= 1;
                    let mut j = 0;
                    while j < n && j < i32::from(ov.num_adj_edges) {
                        let ie = p.iedge[ov.iedge + j as usize] as usize;
                        p.edge[ie].cap = apc.n_old_caps_vert[i][j as usize + 1];
                        j += 1;
                    }
                }
            }
        }
    }
    if b_change_flow & BNS_EF_UPD_H_CHARGE != 0 {
        remove_new(p, apc, true);
        restore_caps(p, apc, true);
    } else {
        restore_caps(p, apc, false);
        remove_new(p, apc, false);
    }
    0
}

/// `bSetBnsToCheckAltPath`: temporarily modifies the network so that an
/// augmenting path between two vertices means that H (or a charge) can move
/// from `v_single` to `v_double`
#[allow(clippy::too_many_arguments)]
pub fn set_bns_to_check_alt_path(
    p: &mut BnStruct,
    v_double: i32,
    v_single: i32,
    type_: u16,
    path_type: i32,
    apc: &mut AltPathChanges,
    fcd: &mut [BnsFlowChanges],
    n_dots: &mut i32,
) -> i32 {
    use super::ichierr::BNS_BOND_ERR as BOND_ERR;
    let hmode = matches!(
        path_type,
        ALT_PATH_MODE_REM2H_CHG
            | ALT_PATH_MODE_ADD2H_CHG
            | ALT_PATH_MODE_REM2H_TST
            | ALT_PATH_MODE_ADD2H_TST
    );
    if p.vert[v_double as usize].st_edge.flow == 0 && !hmode {
        return BNS_CHK_ALTP_NO_ALTPATH;
    }
    let (v1, v2) = (v_single, v_double);
    let pv1 = p.vert[v1 as usize];
    let pv2 = p.vert[v2 as usize];
    let neighbor_of = |p: &BnStruct, pv: &BnsVertex, i: usize, v: Vertex| -> (Vertex, i32) {
        let ie = p.iedge[pv.iedge + i];
        (i32::from(p.edge[ie as usize].neighbor12) ^ v, ie)
    };

    if p.pb_taut_flags & TG_FLAG_TEST_TAUT2_SALTS != 0
        && path_type == ALT_PATH_MODE_4_SALT2
        && type_ & BNS_VERT_TYPE_ENDPOINT != 0
    {
        *apc = AltPathChanges::default();
        fcd[0].iedge = NO_VERTEX;
        *n_dots = 0;
        if v1 == v2 {
            return BNS_CHK_ALTP_SAME_VERTEX;
        }
        let mut n = 0;
        for i in 0..pv1.num_adj_edges as usize {
            let (v, _) = neighbor_of(p, &pv1, i, v1);
            if v == v2 {
                continue;
            }
            n += i32::from(p.vert[v as usize].st_edge.cap > 0);
        }
        if n == 0 {
            return BNS_CHK_ALTP_NO_ALTPATH;
        }
        let (mut v1_act, mut v2_act) = (v1, v2);
        let (mut v1t, mut v2t) = (NO_VERTEX, NO_VERTEX);
        if pv1.type_ & type_ == type_ {
            v1t = get_group_vertex(p, v1, type_);
            if is_bns_error(v1t) {
                return v1t;
            }
            if v1t != NO_VERTEX {
                v1_act = v1t;
            }
        }
        if pv2.type_ & type_ == type_ {
            v2t = get_group_vertex(p, v2, type_);
            if is_bns_error(v2t) {
                return v2t;
            }
            if v2t != NO_VERTEX {
                v2_act = v2t;
            }
        }
        if v1t != NO_VERTEX && v1t == v2t {
            return BNS_CHK_ALTP_SAME_TGROUP;
        }
        let (mut set1, mut set2) = (false, false);
        let mut iapc = 0;
        if v1t != NO_VERTEX {
            let v_new = add_new_vertex(p, v1t, 1, 0, 1, n_dots);
            if is_bns_error(v_new) {
                return v_new;
            }
            apc.v_new_vertex[iapc] = v_new;
            apc.b_set_new[iapc] = 1;
            set1 = true;
            iapc += 1;
        }
        if v2t == NO_VERTEX {
            let v_new = add_new_vertex(p, v2, 1, 0, 1, n_dots);
            if is_bns_error(v_new) {
                return v_new;
            }
            apc.v_new_vertex[iapc] = v_new;
            apc.b_set_new[iapc] = 1;
            set2 = true;
        }
        let mut iapc = 0;
        if !set1 {
            if v1t != NO_VERTEX {
                return BOND_ERR;
            }
            let n = add_st_cap_to_a_vertex(
                p,
                v1,
                v2_act,
                &mut apc.n_old_caps_vert[iapc],
                n_dots,
                false,
            );
            apc.b_set_old_caps_vert[iapc] = n;
            apc.v_old_vert[iapc] = v1;
            iapc += 1;
        }
        if !set2 {
            if v2t == NO_VERTEX {
                return BOND_ERR;
            }
            let n = add_st_cap_to_a_vertex(
                p,
                v2t,
                v1_act,
                &mut apc.n_old_caps_vert[iapc],
                n_dots,
                false,
            );
            apc.b_set_old_caps_vert[iapc] = n;
            apc.v_old_vert[iapc] = v2t;
        }
        if *n_dots < 0 || *n_dots % 2 != 0 {
            return super::ichierr::BNS_SET_ALTP_ERR;
        }
        return BNS_CHK_ALTP_SET_SUCCESS;
    }

    if hmode {
        let b_donors = path_type == ALT_PATH_MODE_REM2H_CHG || path_type == ALT_PATH_MODE_REM2H_TST;
        *apc = AltPathChanges::default();
        fcd[0].iedge = NO_VERTEX;
        *n_dots = 0;
        let mut cap = 1;
        let mut adjacent = [0i32; 2];
        for (k, (pv, v, other)) in [(&pv1, v1, v2), (&pv2, v2, v1)].into_iter().enumerate() {
            let mut n = 0;
            let mut adj = 0;
            for i in 0..pv.num_adj_edges as usize {
                let (w, ie) = neighbor_of(p, pv, i, v);
                let e = &p.edge[ie as usize];
                n += if b_donors {
                    i32::from(p.vert[w as usize].st_edge.cap > 0)
                } else {
                    i32::from(e.flow > 0)
                };
                if b_donors && w == other && e.flow < MAX_BOND_EDGE_CAP {
                    adj += 1;
                }
            }
            adjacent[k] = adj;
            if n == 0 && adj == 0 {
                return BNS_CHK_ALTP_NO_ALTPATH;
            }
        }
        let b_adjacent_donors = adjacent[1];
        let (mut v1_act, mut v2_act) = (v1, v2);
        let (mut v1t, mut v2t) = (NO_VERTEX, NO_VERTEX);
        if pv1.type_ & type_ == type_ {
            v1t = get_group_vertex(p, v1, type_);
            if v1t == BOND_ERR {
                v1t = NO_VERTEX;
            } else if is_bns_error(v1t) {
                return v1t;
            } else if v1t != NO_VERTEX {
                v1_act = v1t;
            }
        }
        if pv2.type_ & type_ == type_ {
            v2t = get_group_vertex(p, v2, type_);
            if v2t == BOND_ERR {
                v2t = NO_VERTEX;
            } else if is_bns_error(v2t) {
                return v2t;
            } else if v2t != NO_VERTEX {
                v2_act = v2t;
            }
        }
        if v1t != NO_VERTEX && v1t == v2t {
            cap = 2;
        }
        let (mut set1, mut set2) = (false, false);
        let mut iapc = 0;
        if b_donors == (v1t != NO_VERTEX) {
            let v_new = add_new_vertex(p, v1_act, cap, 0, 1, n_dots);
            if is_bns_error(v_new) {
                return v_new;
            }
            apc.v_new_vertex[iapc] = v_new;
            apc.b_set_new[iapc] = 1;
            set1 = true;
            iapc += 1;
        }
        if b_donors == (v2t != NO_VERTEX) && cap == 1 {
            let v_new = add_new_vertex(p, v2_act, cap, 0, 1, n_dots);
            if is_bns_error(v_new) {
                return v_new;
            }
            apc.v_new_vertex[iapc] = v_new;
            apc.b_set_new[iapc] = 1;
            set2 = true;
        } else if b_donors == (v2t != NO_VERTEX) {
            set2 = true;
        }
        let mut iapc = 0;
        let adj_donors = b_adjacent_donors != 0;
        if !set1 {
            if b_donors == (v1t != NO_VERTEX) {
                return BOND_ERR;
            }
            let n = add_st_cap_to_a_vertex(
                p,
                v1_act,
                v2_act,
                &mut apc.n_old_caps_vert[iapc],
                n_dots,
                adj_donors,
            );
            apc.b_set_old_caps_vert[iapc] = n;
            apc.v_old_vert[iapc] = v1_act;
            iapc += 1;
        }
        if !set2 {
            if b_donors == (v2t != NO_VERTEX) {
                return BOND_ERR;
            }
            let n = add_st_cap_to_a_vertex(
                p,
                v2_act,
                v1_act,
                &mut apc.n_old_caps_vert[iapc],
                n_dots,
                adj_donors,
            );
            apc.b_set_old_caps_vert[iapc] = n;
            apc.v_old_vert[iapc] = v2_act;
        }
        if *n_dots < 0 || *n_dots % 2 != 0 {
            return super::ichierr::BNS_SET_ALTP_ERR;
        }
        return BNS_CHK_ALTP_SET_SUCCESS;
    }

    if path_type == ALT_PATH_MODE_REM_PROTON
        && v1 >= 0
        && v2 >= 0
        && pv1.type_ & BNS_VERT_TYPE_ANY_GROUP != 0
        && pv2.type_ & BNS_VERT_TYPE_ANY_GROUP != 0
    {
        if pv1.type_ & BNS_VERT_TYPE_C_GROUP != 0
            && pv1.st_edge.flow == 2 * i32::from(pv1.num_adj_edges)
        {
            return BNS_CHK_ALTP_NO_ALTPATH;
        }
        *apc = AltPathChanges::default();
        fcd[0].iedge = NO_VERTEX;
        *n_dots = 0;
        let v_new = add_new_vertex(p, v2, 1, 0, 1, n_dots);
        if is_bns_error(v_new) {
            return v_new;
        }
        apc.v_new_vertex[0] = v_new;
        apc.b_set_new[0] = 1;
        let n = add_st_cap_to_a_vertex(p, v1, v2, &mut apc.n_old_caps_vert[0], n_dots, false);
        apc.b_set_old_caps_vert[0] = n;
        apc.v_old_vert[0] = v1;
        return BNS_CHK_ALTP_SET_SUCCESS;
    }

    // NEUTRALIZE_ENDPOINTS == 1, FIX_H_CHECKING_TAUT == 1
    *n_dots = 0;
    *apc = AltPathChanges::default();
    fcd[0].iedge = NO_VERTEX;
    if type_ & BNS_VERT_TYPE_ENDPOINT != 0 {
        let mut ifcd = 0usize;
        for (v, type2, flow) in [
            (v1, BNS_VERT_TYPE_C_GROUP, 1),
            (v2, BNS_VERT_TYPE_C_GROUP, 1),
            (v1, BNS_VERT_TYPE_TGROUP, 1),
            (v2, BNS_VERT_TYPE_TGROUP, 0),
        ] {
            let iedge = get_edge_to_group_vertex(p, v, type2);
            if iedge != NO_VERTEX {
                let ret2 = set_flow_to_check_one_bond(p, iedge, flow, &mut fcd[ifcd..]);
                if is_bns_error(ret2) {
                    return ret2;
                }
                *n_dots += ret2;
                while fcd[ifcd].iedge != NO_VERTEX {
                    ifcd += 1;
                }
            }
        }
        let v_new = add_new_vertex(p, v2, 1, 0, 1, n_dots);
        if is_bns_error(v_new) {
            return v_new;
        }
        apc.v_new_vertex[0] = v_new;
        apc.b_set_new[0] = 1;
        let n = add_st_cap_to_a_vertex(p, v1, v2, &mut apc.n_old_caps_vert[0], n_dots, false);
        apc.b_set_old_caps_vert[0] = n;
        apc.v_old_vert[0] = v1;
    } else {
        let v_new = add_new_vertex(p, v2, 1, 0, 1, n_dots);
        if is_bns_error(v_new) {
            return v_new;
        }
        apc.v_new_vertex[0] = v_new;
        apc.b_set_new[0] = 1;
        let n = add_st_cap_to_a_vertex(p, v1, v2, &mut apc.n_old_caps_vert[0], n_dots, false);
        apc.b_set_old_caps_vert[0] = n;
        apc.v_old_vert[0] = v1;
    }
    if *n_dots < 0 || *n_dots % 2 != 0 {
        return super::ichierr::BNS_SET_ALTP_ERR;
    }
    BNS_CHK_ALTP_SET_SUCCESS
}

// ---- transferring the flow of the network to the bonds of the structure ----

pub const BOND_ALT_123: i32 = 5;
pub const BOND_ALT_13: i32 = 6;
pub const BOND_ALT_23: i32 = 7;
pub const BOND_MARK_ALT12: i32 = 0x10;
pub const BOND_MARK_ALT123: i32 = 0x20;
pub const BOND_MARK_ALT13: i32 = 0x30;
pub const BOND_MARK_ALT23: i32 = 0x40;
pub const BOND_MARK_ALT12NS: i32 = 0x50;
pub const BOND_MARK_MASK: i32 = 0x70;
pub const BNS_EF_SAVE_ALL: i32 = BNS_EF_CHNG_FLOW | BNS_EF_CHNG_BONDS | BNS_EF_UPD_RAD_ORI;
pub const BNS_EF_ALTR_NS: i32 = BNS_EF_ALTR_BONDS | BNS_EF_SET_NOSTEREO;
pub const BNS_EF_CHNG_RSTR: i32 = BNS_EF_CHNG_FLOW | BNS_EF_RSTR_FLOW;

/// `RestoreEdgeFlow`
pub fn restore_edge_flow(edge: &mut BnsEdge, delta: i32, b_change_flow: i32) -> i32 {
    match b_change_flow & BNS_EF_CHNG_RSTR {
        BNS_EF_CHNG_FLOW => edge.flow0 = edge.flow,
        BNS_EF_CHNG_RSTR => edge.flow -= delta,
        _ => {}
    }
    0
}

/// `SetAtomBondType`: sets the types of the two halves of a bond from the
/// flow of its edge; returns the number of changed bonds or an error
pub fn set_atom_bond_type(
    edge: &BnsEdge,
    bond_type12: &mut u8,
    bond_type21: &mut u8,
    delta: i32,
    b_change_flow: i32,
) -> i32 {
    use super::ichitaut::{BOND_ALT12NS, BOND_ALTERN, BOND_TAUTOM};
    let mut ret = 0;
    if edge.pass == 0 {
        return 0;
    }
    let (mut flow1, mut flow2);
    match b_change_flow & BNS_EF_CHNG_RSTR {
        0 | BNS_EF_RSTR_FLOW => {
            flow1 = edge.flow0;
            flow2 = flow1 + delta;
        }
        _ => {
            flow2 = edge.flow;
            flow1 = edge.flow0;
        }
    }
    if b_change_flow & BNS_EF_CHNG_BONDS != 0 && b_change_flow & BNS_EF_ALTR_NS != BNS_EF_ALTR_NS {
        let new_bond_type = (flow2 + BOND_TYPE_SINGLE) as u8;
        if *bond_type12 != new_bond_type {
            *bond_type12 = new_bond_type;
            *bond_type21 = new_bond_type;
            ret += 1;
        }
    } else if b_change_flow & BNS_EF_ALTR_BONDS != 0 {
        if flow1 == flow2 {
            return ret;
        }
        if flow1 > flow2 {
            std::mem::swap(&mut flow1, &mut flow2);
        }
        let mut bond_type = i32::from(*bond_type12) & BOND_TYPE_MASK;
        let bond_mark;
        match bond_type {
            t if t == BOND_TYPE_SINGLE || t == BOND_TYPE_DOUBLE || t == BOND_TYPE_TRIPLE => {
                if flow1 == 0 && flow2 == 1 {
                    if b_change_flow & BNS_EF_SET_NOSTEREO != 0 {
                        bond_mark = BOND_MARK_ALT12NS;
                        bond_type = BOND_ALT12NS;
                    } else {
                        bond_mark = BOND_MARK_ALT12;
                        bond_type = BOND_ALTERN;
                    }
                } else if flow1 == 0 && flow2 == 2 {
                    bond_mark = BOND_MARK_ALT13;
                    bond_type = BOND_ALT_13;
                } else if flow1 == 1 && flow2 == 2 {
                    bond_mark = BOND_MARK_ALT23;
                    bond_type = BOND_ALT_23;
                } else {
                    return BNS_BOND_ERR;
                }
            }
            BOND_TAUTOM => {
                if flow1 == 0 && flow2 == 1 {
                    bond_mark = BOND_MARK_ALT12NS;
                } else if flow1 == 1 && flow2 == 2 {
                    bond_mark = BOND_MARK_ALT23;
                } else {
                    return BNS_BOND_ERR;
                }
            }
            _ => {
                let mut new_bond_type = bond_type;
                let mut mark = i32::from(*bond_type12) & BOND_MARK_MASK;
                match mark {
                    BOND_MARK_ALT12 | BOND_MARK_ALT12NS => {
                        if mark == BOND_MARK_ALT12
                            && b_change_flow & BNS_EF_SET_NOSTEREO != 0
                            && flow1 == 0
                            && flow2 == 1
                        {
                            mark = BOND_MARK_ALT12NS;
                            new_bond_type = BOND_ALT12NS;
                        } else if flow1 == 2 || flow2 == 2 {
                            mark = BOND_MARK_ALT123;
                            new_bond_type = BOND_ALT_123;
                        }
                    }
                    BOND_MARK_ALT13 => {
                        if flow1 == 1 || flow2 == 1 {
                            mark = BOND_MARK_ALT123;
                            new_bond_type = BOND_ALT_123;
                        }
                    }
                    BOND_MARK_ALT23 => {
                        if flow1 == 0 || flow2 == 0 {
                            mark = BOND_MARK_ALT123;
                            new_bond_type = BOND_ALT_123;
                        }
                    }
                    BOND_MARK_ALT123 => {}
                    0 => {
                        if flow1 == 0 && flow2 == 1 {
                            mark = BOND_MARK_ALT12;
                        } else if flow1 == 0 && flow2 == 2 {
                            mark = BOND_MARK_ALT13;
                        } else if flow1 == 1 && flow2 == 2 {
                            mark = BOND_MARK_ALT23;
                        } else {
                            return BNS_BOND_ERR;
                        }
                    }
                    _ => return BNS_BOND_ERR,
                }
                match bond_type {
                    BOND_TAUTOM => {}
                    BOND_ALTERN | BOND_ALT12NS | BOND_ALT_123 | BOND_ALT_13 | BOND_ALT_23 => {
                        bond_type = new_bond_type;
                    }
                    _ => return BNS_BOND_ERR,
                }
                bond_mark = mark;
            }
        }
        let nb = (bond_type | bond_mark) as u8;
        if nb != *bond_type12 {
            *bond_type12 = nb;
            *bond_type21 = nb;
            ret += 1;
        }
    }
    ret
}

/// `RestoreBnStructFlow`
pub fn restore_bn_struct_flow(p: &mut BnStruct, b_change_flow: i32) -> i32 {
    let mut err = 0;
    for pass in (0..p.num_altp as usize).rev() {
        p.alt_path = Some(pass);
        let (mut v1, n, mut delta, v_last) = {
            let ap = &p.altp[pass];
            (
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_PATH_LEN].number(),
                ap[I_ALTP_FLOW].flow[0],
                ap[I_ALTP_END_ATOM].number(),
            )
        };
        let mut v2 = NO_VERTEX;
        if b_change_flow & BNS_EF_CHNG_RSTR == BNS_EF_CHNG_RSTR {
            p.vert[v1 as usize].st_edge.flow -= delta;
        } else if b_change_flow & BNS_EF_SAVE_ALL == BNS_EF_SAVE_ALL {
            p.vert[v1 as usize].st_edge.flow0 = p.vert[v1 as usize].st_edge.flow;
        }
        for i in 0..n as usize {
            let ineigh1 = p.altp[pass][I_ALTP_NEIGHBOR + i].ineigh(0) as usize;
            let ie = p.iedge[p.vert[v1 as usize].iedge + ineigh1] as usize;
            v2 = i32::from(p.edge[ie].neighbor12) ^ v1;
            restore_edge_flow(&mut p.edge[ie], delta, b_change_flow);
            p.edge[ie].pass = 0;
            delta = -delta;
            v1 = v2;
        }
        if v2 != v_last {
            err = BNS_PROGRAM_ERR;
        } else if b_change_flow & BNS_EF_CHNG_RSTR == BNS_EF_CHNG_RSTR {
            p.vert[v2 as usize].st_edge.flow += delta;
        } else if b_change_flow & BNS_EF_SAVE_ALL == BNS_EF_SAVE_ALL {
            p.vert[v2 as usize].st_edge.flow0 = p.vert[v2 as usize].st_edge.flow;
        }
    }
    err
}

/// `SetAtomRadAndChemValFromVertexCapFlow`
pub fn set_atom_rad_and_chem_val_from_vertex_cap_flow(
    p: &BnStruct,
    atom: &mut [InpAtom],
    v1: usize,
) -> i32 {
    let vert = &p.vert[v1];
    let at = &mut atom[v1];
    let mut n_changes = 0;
    if vert.st_edge.pass == 0 {
        return 0;
    }
    let c = i32::from(at.chem_bonds_valence) - i32::from(at.valence);
    if c >= 0 && c != vert.st_edge.flow {
        at.chem_bonds_valence = (i32::from(at.valence) + vert.st_edge.flow) as i8;
        n_changes += 1;
    }
    let radical = match vert.st_edge.cap - vert.st_edge.flow {
        0 => 0,
        1 => super::incomdef::RADICAL_DOUBLET,
        2 => super::incomdef::RADICAL_TRIPLET,
        _ => return BNS_BOND_ERR,
    } as i8;
    if radical != at.radical {
        at.radical = radical;
        n_changes += 1;
    }
    n_changes
}

/// `SetBondsFromBnStructFlow`
pub fn set_bonds_from_bn_struct_flow(
    p: &mut BnStruct,
    at: &mut [InpAtom],
    num_atoms: i32,
    b_change_flow0: i32,
) -> i32 {
    let b_change_flow = b_change_flow0 & !BNS_EF_SET_NOSTEREO;
    let mut err = 0;
    let mut ret = 0;
    for pass in (0..p.num_altp as usize).rev() {
        p.alt_path = Some(pass);
        let (mut v1, n, mut delta, v_last) = {
            let ap = &p.altp[pass];
            (
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_PATH_LEN].number(),
                ap[I_ALTP_FLOW].flow[0],
                ap[I_ALTP_END_ATOM].number(),
            )
        };
        let mut b_change_flow_add;
        if b_change_flow0 & BNS_EF_SET_NOSTEREO != 0
            && (p.vert[v1 as usize].st_edge.cap0 > p.vert[v1 as usize].st_edge.flow0
                || p.vert[v_last as usize].st_edge.cap0 > p.vert[v_last as usize].st_edge.flow0)
        {
            b_change_flow_add = BNS_EF_SET_NOSTEREO;
            ret |= 2;
        } else {
            b_change_flow_add = 0;
        }
        if b_change_flow & BNS_EF_CHNG_RSTR != BNS_EF_CHNG_RSTR
            && b_change_flow & BNS_EF_SAVE_ALL == BNS_EF_SAVE_ALL
            && v1 < num_atoms
        {
            let r = set_atom_rad_and_chem_val_from_vertex_cap_flow(p, at, v1 as usize);
            if r < 0 {
                err = BNS_PROGRAM_ERR;
            } else {
                ret |= i32::from(r > 0);
            }
        }
        p.vert[v1 as usize].st_edge.pass = 0;
        let mut v0 = NO_VERTEX;
        let mut v2 = NO_VERTEX;
        for i in 0..n as usize {
            let nb = p.altp[pass][I_ALTP_NEIGHBOR + i];
            let (ineigh1, ineigh2) = (nb.ineigh(0) as usize, nb.ineigh(1) as usize);
            let ie = p.iedge[p.vert[v1 as usize].iedge + ineigh1] as usize;
            v2 = i32::from(p.edge[ie].neighbor12) ^ v1;
            if b_change_flow & BNS_EF_CHNG_BONDS != 0 && v1 < num_atoms {
                if v0 >= num_atoms && v2 < num_atoms {
                    at[v1 as usize].chem_bonds_valence += delta as i8;
                } else if v0 < num_atoms && v2 >= num_atoms && v0 != NO_VERTEX {
                    at[v1 as usize].chem_bonds_valence -= delta as i8;
                }
            }
            if p.edge[ie].pass != 0 {
                if v1 < num_atoms
                    && ineigh1 < at[v1 as usize].valence as usize
                    && v2 < num_atoms
                    && ineigh2 < at[v2 as usize].valence as usize
                {
                    if b_change_flow0 & BNS_EF_ALTR_NS == BNS_EF_ALTR_NS
                        && b_change_flow0 & BNS_EF_SAVE_ALL == BNS_EF_SAVE_ALL
                    {
                        b_change_flow_add =
                            if at[v1 as usize].n_ring_system != at[v2 as usize].n_ring_system {
                                BNS_EF_ALTR_NS
                            } else {
                                0
                            };
                    }
                    let mut t12 = at[v1 as usize].bond_type[ineigh1];
                    let mut t21 = at[v2 as usize].bond_type[ineigh2];
                    let r = set_atom_bond_type(
                        &p.edge[ie],
                        &mut t12,
                        &mut t21,
                        delta,
                        b_change_flow | b_change_flow_add,
                    );
                    at[v1 as usize].bond_type[ineigh1] = t12;
                    at[v2 as usize].bond_type[ineigh2] = t21;
                    if r < 0 {
                        err = BNS_PROGRAM_ERR;
                    } else {
                        ret |= i32::from(r > 0);
                    }
                }
                p.edge[ie].pass = 0;
            }
            delta = -delta;
            v0 = v1;
            v1 = v2;
        }
        if v2 != v_last {
            err = BNS_PROGRAM_ERR;
        } else if b_change_flow & BNS_EF_CHNG_RSTR != BNS_EF_CHNG_RSTR
            && b_change_flow & BNS_EF_SAVE_ALL == BNS_EF_SAVE_ALL
            && v2 < num_atoms
        {
            let r = set_atom_rad_and_chem_val_from_vertex_cap_flow(p, at, v2 as usize);
            if r < 0 {
                err = BNS_PROGRAM_ERR;
            } else {
                ret |= i32::from(r > 0);
            }
        }
        p.vert[v2 as usize].st_edge.pass = 0;
    }
    if err != 0 { err } else { ret }
}

/// `bSetBondsAfterCheckOneBond`
pub fn set_bonds_after_check_one_bond(
    p: &mut BnStruct,
    fcd: &[BnsFlowChanges],
    n_test_flow: i32,
    at: &mut [InpAtom],
    num_atoms: i32,
    b_change_flow0: i32,
) -> i32 {
    let b_change_flow0 = b_change_flow0 & !BNS_EF_CHNG_RSTR;
    if b_change_flow0 == 0 {
        return 0;
    }
    let mut b_change_flow = b_change_flow0 & !BNS_EF_SET_NOSTEREO;
    let mut n_changes = 0;
    let mut b_error = 0;
    let mut n = 0usize;
    while fcd[n].iedge != NO_VERTEX {
        n += 1;
    }
    if b_change_flow0 & BNS_EF_SET_NOSTEREO != 0 {
        for (ifcd, f) in fcd.iter().enumerate().take(n) {
            let e = p.edge[f.iedge as usize];
            if e.pass == 0 {
                continue;
            }
            let new_flow = if ifcd == 0 && n_test_flow >= 0 {
                n_test_flow
            } else {
                e.flow
            };
            let v1 = i32::from(e.neighbor1);
            let v2 = i32::from(e.neighbor12) ^ v1;
            if v1 < num_atoms && v2 < num_atoms && new_flow != e.flow0 {
                let (a, b) = (&p.vert[v1 as usize].st_edge, &p.vert[v2 as usize].st_edge);
                if (a.cap0 == a.flow0) != (a.cap == a.flow)
                    || (b.cap0 == b.flow0) != (b.cap == b.flow)
                {
                    b_change_flow |= BNS_EF_SET_NOSTEREO;
                    n_changes |= BNS_EF_SET_NOSTEREO;
                }
            }
        }
    }
    for ifcd in (0..n).rev() {
        let ie = fcd[ifcd].iedge as usize;
        if p.edge[ie].pass == 0 {
            continue;
        }
        let new_flow = if ifcd == 0 && n_test_flow >= 0 {
            n_test_flow
        } else {
            p.edge[ie].flow
        };
        let e = p.edge[ie];
        let v1 = i32::from(e.neighbor1);
        let v2 = i32::from(e.neighbor12) ^ v1;
        if v1 < num_atoms && v2 < num_atoms && b_change_flow != 0 && new_flow != e.flow0 {
            let (i1, i2) = (e.neigh_ord[0] as usize, e.neigh_ord[1] as usize);
            let mut t12 = at[v1 as usize].bond_type[i1];
            let mut t21 = at[v2 as usize].bond_type[i2];
            let r = set_atom_bond_type(&e, &mut t12, &mut t21, new_flow - e.flow0, b_change_flow);
            at[v1 as usize].bond_type[i1] = t12;
            at[v2 as usize].bond_type[i2] = t21;
            if !is_bns_error(r) {
                n_changes |= i32::from(r > 0);
            } else {
                b_error = r;
            }
        }
        p.edge[ie].pass = 0;
    }
    if b_error != 0 { b_error } else { n_changes }
}

// ---- radical endpoints: radicals are replaced by temporary alternating cycles ----

pub const MAX_NUM_RAD: usize = 256;
const NUM_BIT: usize = 16; // `bitWord` is `U_SHORT`

/// `NodeSet`: a list of bit sets of vertices
#[derive(Debug, Clone, Default)]
pub struct NodeSet {
    pub bitword: Vec<Vec<u16>>,
    pub len_set: usize,
}

/// `NodeSetCreate`
pub fn node_set_create(n: usize, l: usize) -> NodeSet {
    let len = n.div_ceil(NUM_BIT);
    NodeSet {
        bitword: vec![vec![0; len]; l],
        len_set: len,
    }
}

/// `NodeSetFromRadEndpoints`
pub fn node_set_from_rad_endpoints(
    s: &mut NodeSet,
    k: usize,
    rad_endpoints: &[Vertex],
    num_v: usize,
) {
    let bits = &mut s.bitword[k];
    bits.iter_mut().for_each(|b| *b = 0);
    let mut i = 1;
    while i < num_v {
        let j = rad_endpoints[i] as usize;
        bits[j / NUM_BIT] |= 1 << (j % NUM_BIT);
        i += 2;
    }
}

/// `RemoveFromNodeSet`
pub fn remove_from_node_set(s: &mut NodeSet, k: usize, v: &[Vertex]) {
    if !s.bitword.is_empty() {
        for &j in v {
            let j = j as usize;
            s.bitword[k][j / NUM_BIT] &= !(1 << (j % NUM_BIT));
        }
    }
}

/// `DoNodeSetsIntersect`
pub fn do_node_sets_intersect(s: &NodeSet, k1: usize, k2: usize) -> bool {
    !s.bitword.is_empty() && (0..s.len_set).any(|i| s.bitword[k1][i] & s.bitword[k2][i] != 0)
}

/// `IsNodeSetEmpty`
pub fn is_node_set_empty(s: &NodeSet, k: usize) -> bool {
    s.bitword.is_empty() || s.bitword[k].iter().all(|&b| b == 0)
}

/// `AddNodeSet2ToNodeSet1`
pub fn add_node_set2_to_node_set1(s: &mut NodeSet, k1: usize, k2: usize) {
    if !s.bitword.is_empty() {
        for i in 0..s.len_set {
            let b = s.bitword[k2][i];
            s.bitword[k1][i] |= b;
        }
    }
}

/// `AddNodesToRadEndpoints`; returns the new length or -1 on overflow
pub fn add_nodes_to_rad_endpoints(
    s: &NodeSet,
    k: usize,
    rad_endpoints: &mut [Vertex],
    v_rad: Vertex,
    n_start: i32,
    n_len: i32,
) -> i32 {
    let mut n = n_start;
    if !s.bitword.is_empty() {
        let mut v: Vertex = 0;
        for i in 0..s.len_set {
            let bits = s.bitword[k][i];
            if bits != 0 {
                for j in 0..NUM_BIT {
                    if bits & (1 << j) != 0 {
                        if n >= n_len {
                            return -1;
                        }
                        rad_endpoints[n as usize] = v_rad;
                        n += 1;
                        rad_endpoints[n as usize] = v;
                        n += 1;
                    }
                    v += 1;
                }
            } else {
                v += NUM_BIT as Vertex;
            }
        }
    }
    n
}

/// `RemoveRadEndpoints`
pub fn remove_rad_endpoints(
    p: &mut BnStruct,
    bd: &mut BnData,
    mut at: Option<&mut [InpAtom]>,
) -> i32 {
    for i in (0..bd.n_num_rad_edges).rev() {
        let ie = bd.rad_edges[i as usize];
        if ie < 0 || ie >= p.num_edges {
            return BNS_PROGRAM_ERR;
        }
        let e = p.edge[ie as usize];
        let v1 = i32::from(e.neighbor1);
        let v2 = i32::from(e.neighbor12) ^ v1;
        if ie + 1 != p.num_edges || v1 < 0 || v1 >= p.num_vertices || v2 < 0 || v2 >= p.num_vertices
        {
            return BNS_PROGRAM_ERR;
        }
        let (p1, p2) = (p.vert[v1 as usize], p.vert[v2 as usize]);
        if p.iedge[p2.iedge + usize::from(p2.num_adj_edges) - 1] != ie
            || p.iedge[p1.iedge + usize::from(p1.num_adj_edges) - 1] != ie
        {
            return BNS_PROGRAM_ERR;
        }
        for v in [v2, v1] {
            let pv = &mut p.vert[v as usize];
            pv.num_adj_edges -= 1;
            let at_ = pv.iedge + usize::from(pv.num_adj_edges);
            pv.st_edge.flow -= e.flow;
            p.iedge[at_] = 0;
        }
        if p.vert[v2 as usize].num_adj_edges == 0 && v2 >= p.num_atoms {
            if v2 + 1 != p.num_vertices {
                return BNS_PROGRAM_ERR;
            }
            p.vert[v2 as usize] = BnsVertex::default();
            p.num_vertices -= 1;
        }
        if p.vert[v1 as usize].num_adj_edges == 0 && v1 >= p.num_atoms {
            if v1 + 1 != p.num_vertices {
                return BNS_PROGRAM_ERR;
            }
            p.vert[v1 as usize] = BnsVertex::default();
            p.num_vertices -= 1;
        }
        if let Some(at) = at.as_deref_mut() {
            if v1 < p.num_atoms {
                let st = p.vert[v1 as usize].st_edge;
                let delta = st.cap - st.flow;
                let mut rad = i32::from(at[v1 as usize].radical);
                match delta {
                    0 if rad == RADICAL_DOUBLET => rad = 0,
                    1 if rad != RADICAL_DOUBLET => rad = RADICAL_DOUBLET,
                    _ => {}
                }
                at[v1 as usize].radical = rad as i8;
            }
        }
        p.edge[ie as usize] = BnsEdge::default();
        p.num_edges -= 1;
    }
    bd.n_num_rad_edges = 0;
    bd.n_num_radicals = 0;
    bd.b_rad_srch_mode = RAD_SRCH_NORM;
    0
}

/// `RestoreRadicalsOnly`
pub fn restore_radicals_only(p: &BnStruct, bd: &BnData, mut at: Option<&mut [InpAtom]>) -> i32 {
    for i in (0..bd.n_num_rad_edges).rev() {
        let ie = bd.rad_edges[i as usize];
        if ie < 0 || ie >= p.num_edges {
            return BNS_PROGRAM_ERR;
        }
        let e = &p.edge[ie as usize];
        let v1 = i32::from(e.neighbor1);
        let v2 = i32::from(e.neighbor12) ^ v1;
        if v1 < 0 || v1 >= p.num_atoms || v2 < p.num_atoms || v2 >= p.num_vertices {
            return BNS_PROGRAM_ERR;
        }
        let (p1, p2) = (&p.vert[v1 as usize], &p.vert[v2 as usize]);
        if p.iedge[p2.iedge + usize::from(e.neigh_ord[1])] != ie
            || p.iedge[p1.iedge + usize::from(e.neigh_ord[0])] != ie
        {
            return BNS_PROGRAM_ERR;
        }
        if let Some(at) = at.as_deref_mut() {
            let delta = p1.st_edge.cap - p1.st_edge.flow + e.flow;
            let mut rad = i32::from(at[v1 as usize].radical);
            match delta {
                0 if rad == RADICAL_DOUBLET => rad = 0,
                1 if rad != RADICAL_DOUBLET => rad = RADICAL_DOUBLET,
                _ => {}
            }
            at[v1 as usize].radical = rad as i8;
        }
    }
    0
}

fn sort_rad_pairs(r: &mut [Vertex], num: usize) {
    let mut pairs: Vec<[Vertex; 2]> = (0..num / 2).map(|i| [r[2 * i], r[2 * i + 1]]).collect();
    pairs.sort();
    for (i, pr) in pairs.iter().enumerate() {
        r[2 * i] = pr[0];
        r[2 * i + 1] = pr[1];
    }
}

/// `SetRadEndpoints`
pub fn set_rad_endpoints(p: &mut BnStruct, bd: &mut BnData, rad_srch_mode: i32) -> i32 {
    if p.tot_st_cap <= p.tot_st_flow {
        return 0;
    }
    bd.n_num_rad_endpoints = 0;
    bd.n_num_rad_edges = 0;
    bd.b_rad_srch_mode = rad_srch_mode;
    p.alt_path = Some(0);
    p.b_change_flow = 0;
    let mut ret = balanced_network_search(p, bd, BNS_EF_RAD_SRCH);
    re_init_bn_data(bd);
    re_init_bn_struct_alt_paths(p);
    if ret == 0 && bd.n_num_rad_endpoints >= 2 {
        let n = bd.n_num_rad_endpoints as usize;
        sort_rad_pairs(&mut bd.rad_endpoints, n);
        let mut n_num_radicals = 0;
        let mut n_dots = 0;
        let mut i = 0usize;
        let mut j;
        while i < n {
            let w_rad = bd.rad_endpoints[i];
            let pr = p.vert[w_rad as usize];
            let mut delta = pr.st_edge.cap - (pr.st_edge.flow & EDGE_FLOW_ST_MASK);
            if delta <= 0 {
                delta = 1;
            }
            let mut n_edges = 0;
            j = i;
            while j < n && w_rad == bd.rad_endpoints[j] {
                n_edges += 1;
                j += 2;
            }
            let v_rad = add_new_vertex(p, w_rad, delta, delta, n_edges + 1, &mut n_dots);
            if is_bns_error(v_rad) {
                remove_rad_endpoints(p, bd, None);
                return v_rad;
            }
            let pr = p.vert[v_rad as usize];
            bd.rad_edges[bd.n_num_rad_edges as usize] =
                p.iedge[pr.iedge + usize::from(pr.num_adj_edges) - 1];
            bd.n_num_rad_edges += 1;
            let mut k = i;
            while k < j {
                bd.rad_endpoints[k] = v_rad;
                k += 2;
            }
            n_num_radicals += 1;
            i = j;
        }
        i = 0;
        while i < n {
            let v_rad = bd.rad_endpoints[i];
            j = i;
            while j < n && v_rad == bd.rad_endpoints[j] {
                let v_endp = bd.rad_endpoints[j + 1];
                ret = add_new_edge(p, v_rad, v_endp, 1, 0);
                if is_bns_error(ret) {
                    remove_rad_endpoints(p, bd, None);
                    return ret;
                }
                bd.rad_edges[bd.n_num_rad_edges as usize] = ret;
                bd.n_num_rad_edges += 1;
                j += 2;
            }
            i = j;
        }
        bd.n_num_radicals = n_num_radicals;
        return n_num_radicals;
    }
    0
}

/// `SetRadEndpoints2`
pub fn set_rad_endpoints2(p: &mut BnStruct, bd: &mut BnData, rad_srch_mode: i32) -> i32 {
    use super::ichierr::{BNS_CAP_FLOW_ERR, BNS_RADICAL_ERR};
    let mut delta = 1;
    let mut ret = 0;
    let mut v_rad_list: Vec<Vertex> = Vec::new();
    let mut v_rad_equl: Vec<Vertex> = Vec::new();
    let mut vert_set = NodeSet::default();
    if p.tot_st_cap <= p.tot_st_flow {
        return 0;
    }
    for i in 0..p.num_atoms as usize {
        if p.vert[i].st_edge.cap - delta == (p.vert[i].st_edge.flow & EDGE_FLOW_ST_MASK)
            && v_rad_list.len() < MAX_NUM_RAD
        {
            p.vert[i].st_edge.cap -= delta;
            p.tot_st_cap -= delta;
            v_rad_equl.push(v_rad_list.len() as Vertex);
            v_rad_list.push(i as Vertex);
        }
    }
    let n_num_rad = v_rad_list.len();
    if p.tot_st_cap - p.tot_st_flow > n_num_rad as i32 {
        return BNS_CAP_FLOW_ERR;
    }
    macro_rules! fail {
        ($r:expr) => {{
            remove_rad_endpoints(p, bd, None);
            return $r;
        }};
    }
    for j in 0..n_num_rad {
        let i = v_rad_list[j] as usize;
        bd.n_num_rad_endpoints = 0;
        bd.n_num_rad_edges = 0;
        bd.b_rad_srch_mode = rad_srch_mode;
        p.alt_path = Some(0);
        p.b_change_flow = 0;
        p.vert[i].st_edge.cap += delta;
        p.tot_st_cap += delta;
        ret = balanced_network_search(p, bd, BNS_EF_RAD_SRCH);
        re_init_bn_data(bd);
        re_init_bn_struct_alt_paths(p);
        bd.b_rad_srch_mode = RAD_SRCH_NORM;
        p.vert[i].st_edge.cap -= delta;
        p.tot_st_cap -= delta;
        if is_bns_error(ret) {
            fail!(ret);
        } else if ret != 0 {
            fail!(BNS_RADICAL_ERR);
        }
        if bd.n_num_rad_endpoints >= 2 {
            let n = bd.n_num_rad_endpoints as usize;
            sort_rad_pairs(&mut bd.rad_endpoints, n);
            if bd.rad_endpoints[0] != i as Vertex || bd.rad_endpoints[n - 2] != i as Vertex {
                fail!(BNS_RADICAL_ERR);
            }
            if n_num_rad > 1 {
                if vert_set.bitword.is_empty() {
                    vert_set = node_set_create(p.num_atoms as usize, n_num_rad);
                }
                node_set_from_rad_endpoints(&mut vert_set, j, &bd.rad_endpoints, n);
                remove_from_node_set(&mut vert_set, j, &v_rad_list);
            }
        }
    }
    for &i in &v_rad_list {
        p.vert[i as usize].st_edge.cap += delta;
        p.tot_st_cap += delta;
    }
    if n_num_rad > 1 {
        for i in 0..n_num_rad {
            if v_rad_equl[i] != i as Vertex {
                continue;
            }
            loop {
                let mut n = 0;
                for j in i + 1..n_num_rad {
                    if v_rad_equl[j] != j as Vertex {
                        continue;
                    }
                    if do_node_sets_intersect(&vert_set, i, j) {
                        add_node_set2_to_node_set1(&mut vert_set, i, j);
                        v_rad_equl[j] = i as Vertex;
                        n += 1;
                    }
                }
                if n == 0 {
                    break;
                }
            }
        }
        let mut n: i32 = 0;
        for i in 0..n_num_rad {
            if i as Vertex == v_rad_equl[i] {
                if !is_node_set_empty(&vert_set, i) {
                    for j in i + 1..n_num_rad {
                        if i as Vertex == v_rad_equl[j] {
                            bd.rad_endpoints[n as usize] = v_rad_list[i];
                            n += 1;
                            bd.rad_endpoints[n as usize] = -v_rad_list[j] - 2;
                            n += 1;
                        }
                    }
                    n = add_nodes_to_rad_endpoints(
                        &vert_set,
                        i,
                        &mut bd.rad_endpoints,
                        v_rad_list[i],
                        n,
                        bd.max_len_pu_pv,
                    );
                    if n < 0 {
                        fail!(BNS_RADICAL_ERR);
                    }
                } else {
                    bd.rad_endpoints[n as usize] = v_rad_list[i];
                    n += 1;
                    bd.rad_endpoints[n as usize] = -1;
                    n += 1;
                }
            }
        }
        bd.n_num_rad_endpoints = n;
    } else if n_num_rad == 1 && bd.n_num_rad_endpoints == 0 {
        bd.rad_endpoints[0] = v_rad_list[0];
        bd.rad_endpoints[1] = -1;
        bd.n_num_rad_endpoints = 2;
    }
    if ret == 0 && bd.n_num_rad_endpoints >= 2 {
        let n = bd.n_num_rad_endpoints as usize;
        let mut n_num_radicals = 0;
        let mut n_dots = 0;
        let mut i = 0usize;
        let mut j;
        while i < n {
            let w_rad = bd.rad_endpoints[i];
            let pr = p.vert[w_rad as usize];
            delta = pr.st_edge.cap - (pr.st_edge.flow & EDGE_FLOW_ST_MASK);
            if delta <= 0 {
                delta = 1;
            }
            let mut n_edges = 0;
            j = i;
            while j < n && w_rad == bd.rad_endpoints[j] {
                n_edges += i32::from(bd.rad_endpoints[j + 1] != -1);
                j += 2;
            }
            let v_rad = add_new_vertex(p, w_rad, delta, delta, n_edges + 1, &mut n_dots);
            if is_bns_error(v_rad) {
                fail!(v_rad);
            }
            let pr = p.vert[v_rad as usize];
            bd.rad_edges[bd.n_num_rad_edges as usize] =
                p.iedge[pr.iedge + usize::from(pr.num_adj_edges) - 1];
            bd.n_num_rad_edges += 1;
            let mut k = i;
            while k < j {
                bd.rad_endpoints[k] = v_rad;
                k += 2;
            }
            n_num_radicals += 1;
            i = j;
        }
        i = 0;
        while i < n {
            let v_rad = bd.rad_endpoints[i];
            j = i;
            while j < n && v_rad == bd.rad_endpoints[j] {
                let mut v_endp = bd.rad_endpoints[j + 1];
                j += 2;
                if v_endp == -1 {
                    continue;
                }
                let edge_flow = if v_endp < 0 {
                    v_endp = -v_endp - 2;
                    1
                } else {
                    0
                };
                ret = add_new_edge(p, v_rad, v_endp, 1, edge_flow);
                if is_bns_error(ret) {
                    fail!(ret);
                }
                bd.rad_edges[bd.n_num_rad_edges as usize] = ret;
                bd.n_num_rad_edges += 1;
            }
            i = j;
        }
        bd.n_num_radicals = n_num_radicals;
        return n_num_radicals;
    }
    0
}

// ---- atom charge/H types (`GetAtomChargeType`) ----

use super::inpdef::{
    ATT_ACIDIC_CO, ATT_ACIDIC_S, ATT_ATOM_N, ATT_ATOM_P, ATT_HAL_ACID, ATT_HAL_ANION, ATT_N_O,
    ATT_NO, ATT_NONE, ATT_NP_MINUS_V23, ATT_O_PLUS, ATT_OH_MINUS, ATT_OO, ATT_OTHER_NEG_O,
    ATT_OTHER_ZO, ATT_PROTON, ATT_ZOO,
};

/// `AT_TYPE_TOTALS`
pub mod attot {
    pub const NUM_NP_PLUS: usize = 0;
    pub const NUM_NP_PROTON: usize = 1;
    pub const NUM_NP_H: usize = 2;
    pub const NUM_N_MINUS: usize = 3;
    pub const NUM_NP: usize = 4;
    pub const NUM_ON: usize = 5;
    pub const NUM_COH: usize = 6;
    pub const NUM_CSH: usize = 7;
    pub const NUM_ZOH: usize = 8;
    pub const NUM_OOH: usize = 9;
    pub const NUM_ZOOH: usize = 10;
    pub const NUM_NOH: usize = 11;
    pub const NUM_N_OH: usize = 12;
    pub const NUM_CO: usize = 13;
    pub const NUM_ZO: usize = 14;
    pub const NUM_NO: usize = 15;
    pub const NUM_N_O: usize = 16;
    pub const NUM_CO_MINUS: usize = 17;
    pub const NUM_CS_MINUS: usize = 18;
    pub const NUM_ZO_MINUS: usize = 19;
    pub const NUM_OO_MINUS: usize = 20;
    pub const NUM_ZOO_MINUS: usize = 21;
    pub const NUM_NO_MINUS: usize = 22;
    pub const NUM_N_O_MINUS: usize = 23;
    pub const NUM_O_MINUS: usize = 24;
    pub const NUM_OH_PLUS: usize = 25;
    pub const NUM_O_PLUS: usize = 26;
    pub const NUM_PROTON: usize = 27;
    pub const NUM_HAL_ANION: usize = 28;
    pub const NUM_HAL_ACID: usize = 29;
    pub const NUM_ERRORS: usize = 30;
    pub const TOT_CHARGE: usize = 31;
    pub const NUM_CHARGES: usize = 32;
    pub const ARRAY_LEN: usize = 33;
}

/// `ATBIT_*`
pub mod atbit {
    use super::attot::*;
    pub const NP_PLUS: i32 = 1 << NUM_NP_PLUS;
    pub const NP_PROTON: i32 = 1 << NUM_NP_PROTON;
    pub const NP_H: i32 = 1 << NUM_NP_H;
    pub const N_MINUS: i32 = 1 << NUM_N_MINUS;
    pub const NP: i32 = 1 << NUM_NP;
    pub const ON: i32 = 1 << NUM_ON;
    pub const COH: i32 = 1 << NUM_COH;
    pub const CSH: i32 = 1 << NUM_CSH;
    pub const ZOH: i32 = 1 << NUM_ZOH;
    pub const OOH: i32 = 1 << NUM_OOH;
    pub const ZOOH: i32 = 1 << NUM_ZOOH;
    pub const NOH: i32 = 1 << NUM_NOH;
    pub const N_OH: i32 = 1 << NUM_N_OH;
    pub const CO: i32 = 1 << NUM_CO;
    pub const ZO: i32 = 1 << NUM_ZO;
    pub const NO: i32 = 1 << NUM_NO;
    pub const N_O: i32 = 1 << NUM_N_O;
    pub const CO_MINUS: i32 = 1 << NUM_CO_MINUS;
    pub const CS_MINUS: i32 = 1 << NUM_CS_MINUS;
    pub const ZO_MINUS: i32 = 1 << NUM_ZO_MINUS;
    pub const OO_MINUS: i32 = 1 << NUM_OO_MINUS;
    pub const ZOO_MINUS: i32 = 1 << NUM_ZOO_MINUS;
    pub const NO_MINUS: i32 = 1 << NUM_NO_MINUS;
    pub const N_O_MINUS: i32 = 1 << NUM_N_O_MINUS;
    pub const O_MINUS: i32 = 1 << NUM_O_MINUS;
    pub const OH_PLUS: i32 = 1 << NUM_OH_PLUS;
    pub const O_PLUS: i32 = 1 << NUM_O_PLUS;
    pub const PROTON: i32 = 1 << NUM_PROTON;
    pub const HAL_ANION: i32 = 1 << NUM_HAL_ANION;
    pub const HAL_ACID: i32 = 1 << NUM_HAL_ACID;
    pub const ERRORS: i32 = 1 << NUM_ERRORS;
}

/// `is_Z_atom`
pub fn is_z_atom(el: u8) -> bool {
    use super::util::*;
    matches!(
        el,
        EL_NUMBER_C
            | EL_NUMBER_N
            | EL_NUMBER_P
            | EL_NUMBER_AS
            | EL_NUMBER_SB
            | EL_NUMBER_S
            | EL_NUMBER_SE
            | EL_NUMBER_TE
            | EL_NUMBER_CL
            | EL_NUMBER_BR
            | EL_NUMBER_I
    )
}

fn is_chalcogen(el: u8) -> bool {
    use super::util::*;
    el == EL_NUMBER_O || el == EL_NUMBER_S || el == EL_NUMBER_SE || el == EL_NUMBER_TE
}

/// `IsZOX`: the number of terminal `=O` (O, S, Se, Te) on the Z neighbor `ord`
/// of `at_x`
pub fn is_zox(atom: &[InpAtom], at_x: usize, ord: usize) -> i32 {
    let z = &atom[atom[at_x].neighbor[ord] as usize];
    let mut num_o = 0;
    for i in 0..z.valence as usize {
        let neigh = z.neighbor[i] as usize;
        if neigh == at_x {
            continue;
        }
        let n = &atom[neigh];
        if n.valence == 1
            && n.chem_bonds_valence == 2
            && n.charge == 0
            && n.radical == 0
            && is_chalcogen(n.el_number)
        {
            num_o += 1;
        }
    }
    num_o
}

/// `GetAtomChargeType`: the type of an atom for proton/charge moving and its
/// bit mask; the totals are updated when given. `b_subtract`: 0 or -2 add,
/// 1 or 2 subtract (-2 and 2 ignore adjacent ions). The mask is `None` when
/// the C code leaves `*pMask` unset.
pub fn get_atom_charge_type(
    atom: &[InpAtom],
    at_no: usize,
    mut totals: Option<&mut [i32]>,
    b_subtract: i32,
) -> (i32, Option<i32>) {
    use super::util::*;
    let at = &atom[at_no];
    let delta = if b_subtract > 0 { -1 } else { 1 };
    let b_no_adj_ion = b_subtract == 0 || b_subtract == 1;
    let mut b_unsat_n_has_terminal_o = 0;
    let mut type_ = i32::from(ATT_NONE);
    let mut mask = 0i32;
    let mut mask_out: Option<i32> = None;
    let el = at.el_number;
    let charge = i32::from(at.charge);
    'exit: {
        'count: {
            if at.radical != 0 && i32::from(at.radical) != super::incomdef::RADICAL_SINGLET {
                break 'exit;
            }
            if is_el_a_metal(i32::from(el)) {
                break 'exit;
            }
            if !(-1..=1).contains(&charge) {
                break 'exit;
            }
            if at.valence == 0
                && charge == 1
                && at.num_h == 0
                && at.radical == 0
                && el == EL_NUMBER_H
            {
                type_ = i32::from(ATT_PROTON);
                mask = atbit::PROTON;
                break 'count;
            }
            if at.valence == 0
                && charge == -1
                && at.num_h == 0
                && at.radical == 0
                && (el == EL_NUMBER_F
                    || el == EL_NUMBER_CL
                    || el == EL_NUMBER_BR
                    || el == EL_NUMBER_I)
            {
                type_ = i32::from(ATT_HAL_ANION);
                mask = atbit::HAL_ANION;
                break 'count;
            }
            if (at.valence == 0
                && charge == 0
                && at.num_h == 1
                && at.radical == 0
                && (el == EL_NUMBER_F
                    || el == EL_NUMBER_CL
                    || el == EL_NUMBER_BR
                    || el == EL_NUMBER_I))
                || (at.valence == 0
                    && charge == 0
                    && at.num_h == 2
                    && at.radical == 0
                    && is_chalcogen(el))
            {
                type_ = i32::from(ATT_HAL_ACID);
                mask = atbit::HAL_ACID;
                break 'count;
            }
            if detect_unusual_el_valence(
                i32::from(el),
                charge,
                i32::from(at.radical),
                i32::from(at.chem_bonds_valence),
                i32::from(at.num_h),
                i32::from(at.valence),
            ) != 0
            {
                break 'exit;
            }
            let (mut num_z, mut num_m, mut num_o) = (0, 0, 0);
            for i in 0..at.valence as usize {
                let neigh = at.neighbor[i] as usize;
                let nb = &atom[neigh];
                if !(-1..=1).contains(&i32::from(nb.charge)) {
                    break 'exit;
                }
                if nb.charge != 0 && at.charge != 0 {
                    if b_no_adj_ion {
                        // FIX_RENUM_BUG_FOR_CASE_OF_ACIDIC_OH_AT_P_PLUS
                        break 'count;
                    }
                    type_ = i32::from(ATT_NONE);
                    mask = 0;
                    break 'count;
                }
                if detect_unusual_el_valence(
                    i32::from(nb.el_number),
                    i32::from(nb.charge),
                    i32::from(nb.radical),
                    i32::from(nb.chem_bonds_valence),
                    i32::from(nb.num_h),
                    i32::from(nb.valence),
                ) != 0
                {
                    break 'exit;
                }
                if is_z_atom(nb.el_number) {
                    num_z += 1;
                }
                if is_el_a_metal(i32::from(nb.el_number)) {
                    num_m += 1;
                }
                num_o += i32::from(nb.el_number == EL_NUMBER_O);
                if el == EL_NUMBER_N
                    && at.valence == 2
                    && at.charge == 0
                    && nb.valence == 1
                    && nb.chem_bonds_valence == 2
                    && is_chalcogen(nb.el_number)
                {
                    b_unsat_n_has_terminal_o += 1;
                }
            }
            if is_chalcogen(el) {
                if charge == 1 {
                    type_ = i32::from(ATT_O_PLUS);
                    mask |= if at.num_h != 0 {
                        atbit::OH_PLUS
                    } else {
                        atbit::O_PLUS
                    };
                } else if at.valence > 1 {
                    break 'exit;
                } else if at.valence != 0 && !(num_z != 0 || num_o != 0) {
                    if num_m == i32::from(at.valence) {
                        break 'exit;
                    }
                    break 'count;
                } else if at.valence != 0 {
                    let neigh = at.neighbor[0] as usize;
                    let nb = &atom[neigh];
                    if nb.charge == 0
                        && nb.el_number == EL_NUMBER_C
                        && nb.chem_bonds_valence > nb.valence
                    {
                        type_ = i32::from(ATT_ACIDIC_CO);
                        if at.num_h == 1 {
                            mask |= atbit::COH;
                        } else if charge == -1 {
                            mask |= atbit::CO_MINUS;
                        } else if at.num_h == 0 && charge == 0 {
                            mask |= atbit::CO;
                        } else {
                            mask |= atbit::ERRORS;
                        }
                    } else if nb.charge == 0
                        && is_chalcogen(nb.el_number)
                        && nb.chem_bonds_valence == nb.valence
                    {
                        type_ = i32::from(ATT_OO);
                        if at.num_h == 1 {
                            mask |= atbit::OOH;
                        } else if charge == -1 {
                            mask |= atbit::OO_MINUS;
                        } else {
                            mask |= atbit::ERRORS;
                        }
                    } else if nb.charge == 0
                        && nb.chem_bonds_valence == nb.valence
                        && nb.el_number == EL_NUMBER_C
                        && el != EL_NUMBER_O
                    {
                        type_ = i32::from(ATT_ACIDIC_S);
                        if at.num_h == 1 {
                            mask |= atbit::CSH;
                        } else if charge == -1 {
                            mask |= atbit::CS_MINUS;
                        } else {
                            mask |= atbit::ERRORS;
                        }
                    } else if nb.el_number == EL_NUMBER_N
                        && nb.valence == 2
                        && (nb.num_h == 0 || (nb.num_h == 1 && nb.charge == 1))
                    {
                        type_ = i32::from(ATT_NO);
                        if at.num_h == 1 {
                            mask |= atbit::NOH;
                        } else if charge == -1 {
                            mask |= atbit::NO_MINUS;
                        } else if nb.charge == 1 || nb.charge == 0 {
                            mask |= atbit::NO;
                        } else {
                            mask |= atbit::ERRORS;
                        }
                    } else if nb.el_number == EL_NUMBER_N {
                        type_ = i32::from(ATT_N_O);
                        if at.num_h == 1 {
                            mask |= atbit::N_OH;
                        } else if charge == -1 {
                            mask |= atbit::N_O_MINUS;
                        } else if nb.charge == 1 {
                            mask |= atbit::N_O;
                        } else {
                            mask |= atbit::ERRORS;
                        }
                    } else if nb.el_number != EL_NUMBER_C
                        && nb.el_number != EL_NUMBER_O
                        && !is_el_a_metal(i32::from(nb.el_number))
                        && nb.chem_bonds_valence > nb.valence
                    {
                        if at.chem_bonds_valence == 1 && is_zox(atom, at_no, 0) != 0 {
                            type_ = i32::from(ATT_ZOO);
                            if at.num_h == 1 {
                                mask |= atbit::ZOOH;
                            } else if charge == -1 {
                                mask |= atbit::ZOO_MINUS;
                            } else {
                                mask |= atbit::ERRORS;
                            }
                        } else {
                            type_ = i32::from(ATT_OTHER_ZO);
                            if at.num_h == 1 {
                                mask |= atbit::ZOH;
                            } else if charge == -1 {
                                mask |= atbit::ZO_MINUS;
                            } else if at.num_h == 0 {
                                mask |= atbit::ZO;
                            } else {
                                mask |= atbit::ERRORS;
                            }
                        }
                    } else if charge == -1 && !is_el_a_metal(i32::from(nb.el_number)) {
                        type_ = i32::from(ATT_OTHER_NEG_O);
                        mask |= atbit::O_MINUS;
                    }
                } else if charge == -1 && at.num_h == 1 {
                    type_ = i32::from(ATT_OH_MINUS);
                    mask |= atbit::O_MINUS;
                }
            } else if (el == EL_NUMBER_N || el == EL_NUMBER_P)
                && (0..=3).contains(&at.valence)
                && i32::from(at.chem_bonds_valence) + i32::from(at.num_h) == 3 + charge
            {
                if at.valence != 0 && num_z == 0 {
                    if num_m == i32::from(at.valence) {
                        break 'exit;
                    }
                    break 'count;
                }
                type_ = i32::from(if el == EL_NUMBER_N {
                    ATT_ATOM_N
                } else {
                    ATT_ATOM_P
                });
                match charge {
                    -1 => {
                        if el == EL_NUMBER_N {
                            mask |= atbit::N_MINUS;
                            if at.num_h != 0 {
                                mask |= atbit::NP_H;
                            } else if at.valence == 1
                                && at.chem_bonds_valence >= 2
                                && i32::from(at.bond_type[0]) & 0xf0 != 0
                            {
                                type_ |= i32::from(ATT_NP_MINUS_V23);
                            }
                        }
                    }
                    0 => {
                        if at.num_h != 0 {
                            mask |= atbit::NP_H;
                        } else if b_unsat_n_has_terminal_o == 1 {
                            mask |= atbit::ON;
                        } else {
                            mask |= atbit::NP;
                        }
                    }
                    1 => {
                        if at.num_h != 0 {
                            mask |= atbit::NP_PROTON;
                        } else if at.chem_bonds_valence > at.valence {
                            mask |= atbit::NP_PLUS;
                        } else {
                            type_ = 0;
                        }
                    }
                    _ => mask |= 1 << attot::NUM_ERRORS,
                }
            }
        }
        // count_mask_bits
        if let Some(t) = totals.as_mut() {
            if mask != 0 && mask & atbit::ERRORS == 0 {
                for (i, slot) in t.iter_mut().enumerate().take(attot::ARRAY_LEN) {
                    if 1i32.checked_shl(i as u32).unwrap_or(0) & mask != 0 {
                        *slot += delta;
                    }
                }
            }
            if charge != 0 {
                t[attot::TOT_CHARGE] += delta * charge;
                t[attot::NUM_CHARGES] += delta;
            }
        }
        mask_out = Some(mask);
    }
    if mask & atbit::ERRORS != 0 {
        type_ = 0;
        if let Some(t) = totals.as_mut() {
            t[attot::NUM_ERRORS] += 1;
        }
    }
    (type_, mask_out)
}
