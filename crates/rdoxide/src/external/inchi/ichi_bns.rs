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

/// The bit sets of a `NodeSet`, `rows` of `len` words each in one allocation (the C code allocates every row
/// separately; a canonicalization allocates thousands of them, almost all never touched)
#[derive(Debug, Clone, Default)]
pub struct BitRows {
    words: Vec<u16>,
    len: usize,
    rows: usize,
}

impl BitRows {
    pub fn is_empty(&self) -> bool {
        self.rows == 0
    }
}

impl std::ops::Index<usize> for BitRows {
    type Output = [u16];
    fn index(&self, i: usize) -> &[u16] {
        assert!(i < self.rows, "node set row out of range");
        &self.words[i * self.len..(i + 1) * self.len]
    }
}

impl std::ops::IndexMut<usize> for BitRows {
    fn index_mut(&mut self, i: usize) -> &mut [u16] {
        assert!(i < self.rows, "node set row out of range");
        &mut self.words[i * self.len..(i + 1) * self.len]
    }
}

/// `NodeSet`: a list of bit sets of vertices
#[derive(Debug, Clone, Default)]
pub struct NodeSet {
    pub bitword: BitRows,
    pub len_set: usize,
}

/// `NodeSetCreate`
pub fn node_set_create(n: usize, l: usize) -> NodeSet {
    let len = n.div_ceil(NUM_BIT);
    NodeSet {
        bitword: BitRows {
            words: vec![0; len * l],
            len,
            rows: l,
        },
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

// ---- applying the flow changes to charges and H counts ----

use super::ichister::{AB_PARITY_UNDF, atom_parity_well_def, get_opposite_sb_atom};
use super::ichitaut::TGroupInfo;
use super::incomdef::NUM_H_ISOTOPES;

pub const MAX_ALT_AATG_ARRAY_LEN: usize = 127;
pub const AATG_MARK_IN_PATH: i8 = 1;
const MAX_NUM_STEREO_BONDS_: usize = 3;
const MIN_NUM_STEREO_BOND_NEIGH: i32 = 2;

/// `BN_AATG` (the members used by the flow analysis)
#[derive(Debug, Clone, Default)]
pub struct BnAatg {
    pub n_num_found: i32,
    pub n_marked_atom: Vec<i8>,
    pub n_at_type_totals: Option<Vec<i32>>,
    pub t_group_info: Option<TGroupInfo>,
}

/// `MarkAtomsAtTautGroups`
pub fn mark_atoms_at_taut_groups(
    p: &mut BnStruct,
    num_atoms: i32,
    aatg: &mut BnAatg,
    mut n_end1: i32,
    mut n_end2: i32,
) -> i32 {
    let mut err = 0;
    let mut c_delta = [0i8; MAX_ALT_AATG_ARRAY_LEN];
    let mut n_vertex = [0u16; MAX_ALT_AATG_ARRAY_LEN];
    let mut n_len_delta = 0usize;
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
        let v_first = v1;
        let mut v2 = NO_VERTEX;
        aatg.n_num_found = 0;
        if n_end1 != v_first && n_end1 != v_last {
            n_end1 = -1;
        }
        if n_end2 != v_first && n_end2 != v_last {
            n_end2 = -1;
        }
        for i in 0..n as usize {
            let ineigh1 = p.altp[pass][I_ALTP_NEIGHBOR + i].ineigh(0) as usize;
            let ie = p.iedge[p.vert[v1 as usize].iedge + ineigh1] as usize;
            v2 = i32::from(p.edge[ie].neighbor12) ^ v1;
            let t1 = if v1 >= 0 {
                p.vert[v1 as usize].type_
            } else {
                0
            };
            let tg1 = t1 & (BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_TEMP) != 0;
            let tg2 = p.vert[v2 as usize].type_ & (BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_TEMP) != 0;
            let mut push = |d: i32, v: i32| {
                if n_len_delta < MAX_ALT_AATG_ARRAY_LEN {
                    c_delta[n_len_delta] = d as i8;
                    n_vertex[n_len_delta] = v as u16;
                    n_len_delta += 1;
                }
            };
            if v1 >= num_atoms
                && tg1
                && (0..num_atoms).contains(&v2)
                && p.vert[v2 as usize].type_ & BNS_VERT_TYPE_ATOM != 0
            {
                push(delta, v2);
            } else if v2 >= num_atoms
                && tg2
                && (0..num_atoms).contains(&v1)
                && p.vert[v1 as usize].type_ & BNS_VERT_TYPE_ATOM != 0
            {
                push(delta, v1);
            } else if (0 <= v1 && v1 == n_end1) || (v1 == n_end2 && (0..num_atoms).contains(&v2)) {
                push(-delta, v1);
            } else if (0 <= v2 && v2 == n_end1) || (v2 == n_end2 && (0..num_atoms).contains(&v1)) {
                push(-delta, v2);
            }
            delta = -delta;
            v1 = v2;
        }
        if v2 != v_last {
            err = BNS_PROGRAM_ERR;
        } else {
            let mut last_i: i32 = -1;
            let mut n_num_found = 0;
            for i in 1..n_len_delta {
                let j = i - 1;
                if (c_delta[j] > 0 && c_delta[i] > 0) || (c_delta[j] < 0 && c_delta[i] < 0) {
                    if j as i32 == last_i {
                        return 0;
                    }
                    if aatg.n_marked_atom[n_vertex[j] as usize] & AATG_MARK_IN_PATH == 0 {
                        n_num_found += 1;
                    }
                    if aatg.n_marked_atom[n_vertex[i] as usize] & AATG_MARK_IN_PATH == 0 {
                        n_num_found += 1;
                    }
                    last_i = i as i32;
                }
            }
            if n_num_found != 4 {
                return 0;
            }
            for i in 1..n_len_delta {
                let j = i - 1;
                if (c_delta[j] > 0 && c_delta[i] > 0) || (c_delta[j] < 0 && c_delta[i] < 0) {
                    for v in [n_vertex[i - 1], n_vertex[i]] {
                        if aatg.n_marked_atom[v as usize] & AATG_MARK_IN_PATH == 0 {
                            aatg.n_marked_atom[v as usize] |= AATG_MARK_IN_PATH;
                            aatg.n_num_found += 1;
                        }
                    }
                }
            }
        }
    }
    if err != 0 { err } else { aatg.n_num_found }
}

/// `EliminatePlusMinusChargeAmbiguity`: replaces (+)--N==(-) with (+)==N--(-)
pub fn eliminate_plus_minus_charge_ambiguity(p: &mut BnStruct, num_atoms: i32) -> i32 {
    let mut err = 0;
    let mut ret = 0;
    for pass in (0..p.num_altp as usize).rev() {
        p.alt_path = Some(pass);
        let (mut v1, n, v_last) = {
            let ap = &p.altp[pass];
            (
                ap[I_ALTP_START_ATOM].number(),
                ap[I_ALTP_PATH_LEN].number(),
                ap[I_ALTP_END_ATOM].number(),
            )
        };
        let mut v0 = NO_VERTEX;
        let mut v2 = NO_VERTEX;
        for i in 0..n as usize {
            let ineigh1 = p.altp[pass][I_ALTP_NEIGHBOR + i].ineigh(0) as usize;
            let ie = p.iedge[p.vert[v1 as usize].iedge + ineigh1] as usize;
            v2 = i32::from(p.edge[ie].neighbor12) ^ v1;
            if v1 < num_atoms
                && ((v0 >= num_atoms && p.vert[v0 as usize].type_ & BNS_VERT_TYPE_C_GROUP != 0)
                    || (v2 >= num_atoms && p.vert[v2 as usize].type_ & BNS_VERT_TYPE_C_GROUP != 0))
            {
                let (mut cg_pos, mut cg_neg) = (0i32, 0i32);
                let (mut neigh_pos, mut neigh_neg) = (-1i32, -1i32);
                let mut n_found = 0;
                let mut k = i32::from(p.vert[v1 as usize].num_adj_edges) - 1;
                while k >= 0 && (neigh_pos < 0 || neigh_neg < 0) {
                    let ne = p.iedge[p.vert[v1 as usize].iedge + k as usize] as usize;
                    let v = i32::from(p.edge[ne].neighbor12) ^ v1;
                    let vt = p.vert[v as usize].type_;
                    if vt & BNS_VERT_TYPE_C_GROUP != 0 {
                        if vt & BNS_VERT_TYPE_C_NEGATIVE != 0 {
                            cg_neg = v;
                            neigh_neg = k;
                        } else {
                            cg_pos = v;
                            neigh_pos = k;
                        }
                        n_found += 1;
                    }
                    k -= 1;
                }
                if n_found == 2 && neigh_pos >= 0 && neigh_neg >= 0 {
                    let base = p.vert[v1 as usize].iedge;
                    let ep = p.iedge[base + neigh_pos as usize] as usize;
                    let en = p.iedge[base + neigh_neg as usize] as usize;
                    if p.edge[ep].flow < p.edge[en].flow {
                        let dflow = p.edge[en].flow - p.edge[ep].flow;
                        p.edge[ep].flow += dflow;
                        p.vert[cg_pos as usize].st_edge.cap += dflow;
                        p.vert[cg_pos as usize].st_edge.flow += dflow;
                        p.edge[en].flow -= dflow;
                        p.vert[cg_neg as usize].st_edge.cap -= dflow;
                        p.vert[cg_neg as usize].st_edge.flow -= dflow;
                        ret += 1;
                    }
                }
            }
            v0 = v1;
            v1 = v2;
        }
        if v2 != v_last {
            err = BNS_PROGRAM_ERR;
        }
    }
    if err != 0 { err } else { ret }
}

/// `AddChangedAtHChargeBNS`
pub fn add_changed_at_h_charge_bns(
    at: &mut [InpAtom],
    num_atoms: i32,
    totals: Option<&mut [i32]>,
    mark: &mut [i8],
) -> i32 {
    let mut totals = totals;
    let mut num = 0;
    for i in 0..num_atoms as usize {
        if mark[i] != 0 {
            mark[i] = 0;
            let (t, _) = get_atom_charge_type(at, i, totals.as_deref_mut(), -2);
            at[i].at_type = t as u16;
            num += 1;
        }
    }
    num
}

/// `AddOrRemoveExplOrImplH`: adds or removes implicit/explicit H of an atom;
/// returns the part of `n_delta` that was done
pub fn add_or_remove_expl_or_impl_h(
    n_delta: i32,
    at: &mut [InpAtom],
    num_atoms: i32,
    at_no: usize,
    tgi: &mut TGroupInfo,
) -> i32 {
    if n_delta == 0 {
        return 0;
    }
    if n_delta > 0 {
        at[at_no].num_h += n_delta as i8;
        tgi.tni.n_num_removed_protons -= 1;
        return n_delta;
    }
    let mut n_num2remove = -n_delta;
    let mut n_removed_expl = tgi.tni.n_num_removed_explicit_h;
    let tot_num_iso_h: i32 = at[at_no].num_iso_h.iter().map(|&x| i32::from(x)).sum();
    let mut num_h = i32::from(at[at_no].num_h);
    let mut num_iso_h = at[at_no].num_iso_h;
    let nat = num_atoms as usize;
    let mut i = 0usize;
    while (i as i32) < n_removed_expl {
        if at[nat + i].neighbor[0] as usize == at_no {
            let orig_no = at[nat + i].orig_at_number;
            n_removed_expl -= 1;
            if n_removed_expl > i as i32 {
                at[nat + i..=nat + n_removed_expl as usize].rotate_left(1);
            }
            if at[at_no].sb_parity[0] != 0 {
                let mut m = 0;
                while m < MAX_NUM_STEREO_BONDS_ && at[at_no].sb_parity[m] != 0 {
                    if at[at_no].sn_orig_at_num[m] == orig_no {
                        if i32::from(at[at_no].valence) >= MIN_NUM_STEREO_BOND_NEIGH {
                            let k = usize::from(at[at_no].sb_ord[m] == 0);
                            at[at_no].sn_ord[m] = k as i8;
                            let nb = at[at_no].neighbor[k] as usize;
                            at[at_no].sn_orig_at_num[m] = at[nb].orig_at_number;
                            if atom_parity_well_def(i32::from(at[at_no].sb_parity[m])) {
                                at[at_no].sb_parity[m] = 3 - at[at_no].sb_parity[m];
                            }
                        } else {
                            at[at_no].sn_ord[m] = -99;
                            at[at_no].sn_orig_at_num[m] = 0;
                            if atom_parity_well_def(i32::from(at[at_no].sb_parity[m])) {
                                let (len, pnxt, _pi2c, pord) =
                                    get_opposite_sb_atom(at, at_no, at[at_no].sb_ord[m] as usize);
                                if len > 0 {
                                    at[at_no].sb_parity[m] = AB_PARITY_UNDF as i8;
                                    at[pnxt].sb_parity[pord] = AB_PARITY_UNDF as i8;
                                }
                            }
                        }
                    }
                    m += 1;
                }
            }
        } else {
            i += 1;
        }
    }
    let mut iso: i32 = -1;
    while iso < NUM_H_ISOTOPES as i32 && n_num2remove > 0 {
        if iso < 0 {
            while tot_num_iso_h < num_h && n_num2remove > 0 {
                num_h -= 1;
                tgi.tni.n_num_removed_protons += 1;
                n_num2remove -= 1;
            }
        } else {
            let u = iso as usize;
            while num_iso_h[u] != 0 && num_h != 0 && n_num2remove > 0 {
                num_h -= 1;
                num_iso_h[u] -= 1;
                tgi.tni.n_num_removed_protons_isotopic[u] += 1;
                tgi.tni.n_num_removed_protons += 1;
                n_num2remove -= 1;
            }
        }
        iso += 1;
    }
    if n_delta + n_num2remove < 0 {
        at[at_no].num_h = num_h as i8;
        at[at_no].num_iso_h = num_iso_h;
        tgi.tni.n_num_removed_explicit_h = n_removed_expl;
    }
    n_delta + n_num2remove
}

/// `SubtractOrChangeAtHChargeBNS`
pub fn subtract_or_change_at_h_charge_bns(
    p: &mut BnStruct,
    at: &mut [InpAtom],
    num_atoms: i32,
    mut totals: Option<&mut [i32]>,
    mut mark: Option<&mut [i8]>,
    mut tgi: Option<&mut TGroupInfo>,
    b_subtract: bool,
) -> i32 {
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
        let mut v0 = NO_VERTEX;
        let mut v2 = NO_VERTEX;
        for i in 0..n as usize {
            let ineigh1 = p.altp[pass][I_ALTP_NEIGHBOR + i].ineigh(0) as usize;
            let ie = p.iedge[p.vert[v1 as usize].iedge + ineigh1] as usize;
            v2 = i32::from(p.edge[ie].neighbor12) ^ v1;
            if v1 < num_atoms && (v0 >= num_atoms || v2 >= num_atoms) {
                let (mut n_delta_h, mut n_delta_charge) = (0, 0);
                if v0 >= num_atoms {
                    let t = p.vert[v0 as usize].type_;
                    if t & BNS_VERT_TYPE_TGROUP != 0 {
                        n_delta_h -= delta;
                    } else if t & BNS_VERT_TYPE_C_GROUP != 0 {
                        n_delta_charge += delta;
                    }
                }
                if v2 >= num_atoms {
                    let t = p.vert[v2 as usize].type_;
                    if t & BNS_VERT_TYPE_TGROUP != 0 {
                        n_delta_h += delta;
                    } else if t & BNS_VERT_TYPE_C_GROUP != 0 {
                        n_delta_charge -= delta;
                    }
                }
                if n_delta_h != 0 || n_delta_charge != 0 {
                    if b_subtract {
                        let m = mark.as_deref_mut().expect("mark");
                        if m[v1 as usize] == 0 {
                            get_atom_charge_type(at, v1 as usize, totals.as_deref_mut(), 2);
                            ret += 1;
                            m[v1 as usize] += 1;
                        }
                    } else {
                        at[v1 as usize].charge += n_delta_charge as i8;
                        if n_delta_h != 0 {
                            add_or_remove_expl_or_impl_h(
                                n_delta_h,
                                at,
                                num_atoms,
                                v1 as usize,
                                tgi.as_deref_mut().expect("t_group_info"),
                            );
                        }
                        ret += 1;
                    }
                }
            }
            delta = -delta;
            v0 = v1;
            v1 = v2;
        }
        if v2 != v_last {
            err = BNS_PROGRAM_ERR;
        }
    }
    if err != 0 { err } else { ret }
}

// ---- does an alternating path exist between two atoms? ----

use super::ichitaut::{
    EndpointInfo, n_get_endpoint_info, n_get_endpoint_info_ket, n_get_endpoint_info_pt_06_00,
    n_get_endpoint_info_pt_13_00, n_get_endpoint_info_pt_16_00, n_get_endpoint_info_pt_18_00,
    n_get_endpoint_info_pt_22_00, n_get_endpoint_info_pt_39_00,
};

pub const ALT_PATH_MODE_TAUTOM_KET: i32 = 10;
pub const ALT_PATH_MODE_TAUTOM_PT_22_00: i32 = 11;
pub const ALT_PATH_MODE_TAUTOM_PT_16_00: i32 = 12;
pub const ALT_PATH_MODE_TAUTOM_PT_06_00: i32 = 13;
pub const ALT_PATH_MODE_TAUTOM_PT_39_00: i32 = 14;
pub const ALT_PATH_MODE_TAUTOM_PT_13_00: i32 = 15;
pub const ALT_PATH_MODE_TAUTOM_PT_18_00: i32 = 16;
pub const BNS_MAX_NUM_FLOW_CHANGES: usize = 1 + 2 * MAX_BOND_EDGE_CAP as usize;

type EndpointInfoFn = fn(&[InpAtom], usize, &mut EndpointInfo) -> i32;

/// `bExistsAltPath`: is there an alternating path allowing to move H (or a
/// charge) from `v_single` to `v_double`? Returns 0, or a positive code
/// (`bSuccess + 4 * nDelta`), or a negative error
#[allow(clippy::too_many_arguments)]
pub fn exists_alt_path(
    p: &mut BnStruct,
    bd: &mut BnData,
    mut aatg: Option<&mut BnAatg>,
    at: &mut [InpAtom],
    num_atoms: i32,
    v_double: i32,
    v_single: i32,
    path_type: i32,
) -> i32 {
    let mut apc = AltPathChanges::default();
    let mut fcd = vec![BnsFlowChanges::default(); 4 * BNS_MAX_NUM_FLOW_CHANGES + 1];
    let mut eif = EndpointInfo::default();
    let mut eif2 = EndpointInfo::default();
    let mut b_change_flow = 0;
    let (vs, vd) = (v_single as usize, v_double as usize);
    let type_: u16;
    let tautomer_check = |at: &[InpAtom],
                          f: EndpointInfoFn,
                          eif: &mut EndpointInfo,
                          eif2: &mut EndpointInfo|
     -> bool {
        if at[vs].endpoint == 0 && (f(at, vs, eif) == 0 || eif.c_donor == 0) {
            return false;
        }
        if at[vd].endpoint == 0 && (f(at, vd, eif2) == 0 || eif2.c_acceptor == 0) {
            return false;
        }
        true
    };
    match path_type {
        ALT_PATH_MODE_TAUTOM
        | ALT_PATH_MODE_TAUTOM_PT_22_00
        | ALT_PATH_MODE_TAUTOM_PT_16_00
        | ALT_PATH_MODE_TAUTOM_PT_06_00
        | ALT_PATH_MODE_TAUTOM_PT_39_00
        | ALT_PATH_MODE_TAUTOM_PT_13_00
        | ALT_PATH_MODE_TAUTOM_PT_18_00
        | ALT_PATH_MODE_TAUTOM_KET => {
            type_ = BNS_VERT_TYPE_ENDPOINT;
            b_change_flow = BNS_EF_CHNG_RSTR;
            let f: EndpointInfoFn = match path_type {
                ALT_PATH_MODE_TAUTOM_PT_22_00 => n_get_endpoint_info_pt_22_00,
                ALT_PATH_MODE_TAUTOM_PT_16_00 => n_get_endpoint_info_pt_16_00,
                ALT_PATH_MODE_TAUTOM_PT_06_00 => n_get_endpoint_info_pt_06_00,
                ALT_PATH_MODE_TAUTOM_PT_39_00 => n_get_endpoint_info_pt_39_00,
                ALT_PATH_MODE_TAUTOM_PT_13_00 => n_get_endpoint_info_pt_13_00,
                ALT_PATH_MODE_TAUTOM_PT_18_00 => n_get_endpoint_info_pt_18_00,
                ALT_PATH_MODE_TAUTOM_KET => n_get_endpoint_info_ket,
                _ => n_get_endpoint_info,
            };
            if !tautomer_check(at, f, &mut eif, &mut eif2) {
                return 0;
            }
        }
        ALT_PATH_MODE_CHARGE => {
            type_ = BNS_VERT_TYPE_C_POINT;
            b_change_flow = BNS_EF_CHNG_RSTR | BNS_EF_ALTR_BONDS;
        }
        ALT_PATH_MODE_4_SALT | ALT_PATH_MODE_4_SALT2 => {
            type_ = BNS_VERT_TYPE_ENDPOINT;
            b_change_flow = BNS_EF_CHNG_RSTR | BNS_EF_ALTR_BONDS;
            if !is_bns_endpoint(p, v_single)
                && (n_get_endpoint_info(at, vs, &mut eif) == 0 || eif.c_donor == 0)
            {
                return 0;
            }
            if !is_bns_endpoint(p, v_double)
                && (n_get_endpoint_info(at, vd, &mut eif) == 0 || eif.c_acceptor == 0)
            {
                return 0;
            }
        }
        ALT_PATH_MODE_REM2H_CHG | ALT_PATH_MODE_REM2H_TST => {
            if path_type == ALT_PATH_MODE_REM2H_CHG {
                b_change_flow |= BNS_EF_ALTR_BONDS;
            }
            b_change_flow |= BNS_EF_CHNG_RSTR;
            type_ = BNS_VERT_TYPE_ENDPOINT;
            if !is_bns_endpoint(p, v_single)
                && (n_get_endpoint_info(at, vs, &mut eif) == 0 || eif.c_donor == 0)
            {
                return 0;
            }
            if !is_bns_endpoint(p, v_double)
                && (n_get_endpoint_info(at, vd, &mut eif) == 0 || eif.c_donor == 0)
            {
                return 0;
            }
        }
        ALT_PATH_MODE_ADD2H_CHG | ALT_PATH_MODE_ADD2H_TST => {
            if path_type == ALT_PATH_MODE_ADD2H_CHG {
                b_change_flow |= BNS_EF_ALTR_BONDS;
            }
            b_change_flow |= BNS_EF_CHNG_RSTR;
            type_ = BNS_VERT_TYPE_ENDPOINT;
            if !is_bns_endpoint(p, v_single)
                && (n_get_endpoint_info(at, vs, &mut eif) == 0 || eif.c_acceptor == 0)
            {
                return 0;
            }
            if !is_bns_endpoint(p, v_double)
                && (n_get_endpoint_info(at, vd, &mut eif) == 0 || eif.c_acceptor == 0)
            {
                return 0;
            }
        }
        ALT_PATH_MODE_REM_PROTON => {
            type_ = 0;
            b_change_flow = (BNS_EF_SAVE_ALL | BNS_EF_UPD_H_CHARGE) | BNS_EF_ALTR_NS;
        }
        _ => {
            type_ = 0;
            b_change_flow = BNS_EF_CHNG_RSTR;
        }
    }

    let mut b_error = 0;
    let mut b_success = 0;
    let mut n_delta = 0;
    let mut ret = set_rad_endpoints2(p, bd, RAD_SRCH_NORM);
    if is_bns_error(ret) {
        return ret;
    }
    let mut n_dots = 0;
    ret = set_bns_to_check_alt_path(
        p,
        v_double,
        v_single,
        type_,
        path_type,
        &mut apc,
        &mut fcd,
        &mut n_dots,
    );
    let mut b_adjust_radicals = false;
    let mut ret_val;
    let mut skip_search = false;
    match ret {
        BNS_CHK_ALTP_NO_ALTPATH => {
            return remove_rad_endpoints(p, bd, None);
        }
        BNS_CHK_ALTP_SAME_TGROUP => {
            b_success = 1;
            skip_search = true;
        }
        BNS_CHK_ALTP_SAME_VERTEX => {
            let r = remove_rad_endpoints(p, bd, None);
            return if r != 0 { r } else { 1 };
        }
        BNS_CHK_ALTP_SET_SUCCESS => {}
        super::ichierr::BNS_CANT_SET_BOND => skip_search = true,
        _ => {
            let _ = remove_rad_endpoints(p, bd, None);
            if is_bns_error(ret) {
                return ret;
            }
            return BNS_PROGRAM_ERR;
        }
    }
    if !skip_search {
        b_adjust_radicals =
            b_change_flow & BNS_EF_UPD_RAD_ORI != 0 && b_change_flow & BNS_EF_RSTR_FLOW == 0;
        ret = run_balanced_network_search(p, bd, b_change_flow);
        if is_bns_error(ret) {
            b_error = ret;
        } else if ret > 0 {
            if 2 * ret >= n_dots {
                n_delta = 2 * ret - n_dots;
                let mut b_do_mark_changed_bonds = 1;
                if let Some(ag) = aatg.as_deref_mut().filter(|a| !a.n_marked_atom.is_empty()) {
                    if ag.n_at_type_totals.is_some() && b_change_flow & BNS_EF_UPD_H_CHARGE != 0 {
                        for m in ag.n_marked_atom.iter_mut().take(num_atoms as usize) {
                            *m = 0;
                        }
                        subtract_or_change_at_h_charge_bns(
                            p,
                            at,
                            num_atoms,
                            ag.n_at_type_totals.as_deref_mut(),
                            Some(&mut ag.n_marked_atom),
                            None,
                            true,
                        );
                        subtract_or_change_at_h_charge_bns(
                            p,
                            at,
                            num_atoms,
                            None,
                            None,
                            ag.t_group_info.as_mut(),
                            false,
                        );
                    } else if ag.n_at_type_totals.is_none() {
                        b_do_mark_changed_bonds =
                            mark_atoms_at_taut_groups(p, num_atoms, ag, v_single, v_double);
                        if b_do_mark_changed_bonds < 0 {
                            b_error = b_do_mark_changed_bonds;
                            b_do_mark_changed_bonds = 0;
                        }
                    }
                }
                if b_do_mark_changed_bonds != 0 {
                    ret_val =
                        set_bonds_after_check_one_bond(p, &fcd, -1, at, num_atoms, b_change_flow);
                    if is_bns_error(ret_val) {
                        b_error = ret_val;
                    }
                    ret = set_bonds_from_bn_struct_flow(p, at, num_atoms, b_change_flow);
                    if is_bns_error(ret) {
                        b_error = ret;
                    } else if ret & 1 == 0 && ret_val & 1 == 0 {
                        b_success = 1;
                    } else if ((ret & 1 != 0 || ret_val & 1 != 0)
                        && b_change_flow & BNS_EF_ALTR_BONDS != 0)
                        || b_change_flow & BNS_EF_UPD_H_CHARGE != 0
                    {
                        b_success = 3;
                    } else {
                        b_error = BNS_BOND_ERR;
                    }
                    if b_error == 0 {
                        if let Some(ag) = aatg.filter(|a| !a.n_marked_atom.is_empty()) {
                            if b_change_flow & BNS_EF_UPD_H_CHARGE != 0 {
                                if b_adjust_radicals {
                                    ret_val = restore_radicals_only(p, bd, Some(at));
                                    if is_bns_error(ret_val) {
                                        b_error = ret_val;
                                    }
                                }
                                add_changed_at_h_charge_bns(
                                    at,
                                    num_atoms,
                                    ag.n_at_type_totals.as_deref_mut(),
                                    &mut ag.n_marked_atom,
                                );
                                if b_change_flow & BNS_EF_CHNG_FLOW != 0 {
                                    eliminate_plus_minus_charge_ambiguity(p, num_atoms);
                                }
                            }
                        }
                    }
                }
            }
            ret = restore_bn_struct_flow(p, b_change_flow & BNS_EF_CHNG_RSTR);
            if is_bns_error(ret) {
                b_error = ret;
            }
        }
    }
    // reinit_BNS
    restore_bns_after_check_alt_path(p, &apc, b_change_flow & BNS_EF_UPD_H_CHARGE);
    restore_flow_after_check_one_bond(p, &fcd);
    ret_val = remove_rad_endpoints(p, bd, if b_adjust_radicals { Some(at) } else { None });
    re_init_bn_struct_alt_paths(p);
    if b_error != 0 {
        b_error
    } else if ret_val != 0 {
        ret_val
    } else {
        b_success + 4 * n_delta
    }
}

/// `bExistsAnyAltPath`: tries both directions
pub fn exists_any_alt_path(
    p: &mut BnStruct,
    bd: &mut BnData,
    at: &mut [InpAtom],
    num_atoms: i32,
    v2: i32,
    v1: i32,
    path_type: i32,
) -> i32 {
    let r1 = exists_alt_path(p, bd, None, at, num_atoms, v2, v1, path_type);
    if r1 > 0 {
        return r1;
    }
    let r2 = exists_alt_path(p, bd, None, at, num_atoms, v1, v2, path_type);
    if r2 > 0 {
        return r2;
    }
    if is_bns_error(r1) {
        return r1;
    }
    if is_bns_error(r2) {
        return r2;
    }
    0
}

// ---- adding t-groups and c-groups to the network ----

use super::ichiparm::{
    TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_MOVE_POS_CHARGES, TG_FLAG_PT_06_00, TG_FLAG_PT_13_00,
    TG_FLAG_PT_16_00, TG_FLAG_PT_18_00, TG_FLAG_PT_22_00, TG_FLAG_PT_39_00,
};
use super::ichitaut::CGroupInfo;

pub const BNS_ADD_SUPER_TGROUP: i32 = 1;

/// `AddTGroups2BnStruct`: adds a fictitious vertex for each t-group and joins
/// the endpoints to it
pub fn add_t_groups_to_bn_struct(
    p: &mut BnStruct,
    at: &[InpAtom],
    num_atoms: i32,
    tgi: Option<&mut TGroupInfo>,
) -> i32 {
    use super::ichitaut::{BOND_ALT12NS, BOND_ALTERN, BOND_TAUTOM};
    let mut ret = 0;
    let Some(tgi) = tgi else { return 0 };
    if tgi.num_t_groups == 0 || tgi.t_group.is_none() {
        return 0;
    }
    let num_tg = tgi.num_t_groups;
    let mut num_edges = p.num_edges;
    let num_vertices = p.num_vertices;
    if num_vertices + num_tg as i32 >= p.max_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    let flags = tgi.b_taut_flags;
    let tg = tgi.t_group.as_mut().expect("t_group");
    let n_max = tg[..num_tg]
        .iter()
        .map(|t| i32::from(t.n_group_number))
        .fold(0, i32::max);
    for v in &mut p.vert[num_vertices as usize..(num_vertices + n_max) as usize] {
        *v = BnsVertex::default();
    }
    if n_max != i32::from(tg[num_tg - 1].n_group_number) {
        tg[..num_tg].sort_by_key(|t| t.n_group_number);
    }
    let mut prev = (num_vertices - 1) as usize;
    for t in tg.iter().take(num_tg) {
        let cur = (num_vertices + i32::from(t.n_group_number) - 1) as usize;
        let iedge = p.vert[prev].iedge + usize::from(p.vert[prev].max_adj_edges);
        let v = &mut p.vert[cur];
        v.iedge = iedge;
        v.max_adj_edges =
            (i32::from(t.n_num_endpoints) + BNS_ADD_EDGES + BNS_ADD_SUPER_TGROUP) as u16;
        v.num_adj_edges = 0;
        v.st_edge.flow = 0;
        v.st_edge.flow0 = 0;
        v.st_edge.cap = 0;
        v.st_edge.cap0 = 0;
        v.type_ = BNS_VERT_TYPE_TGROUP;
        prev = cur;
    }
    for endpoint in 0..num_atoms {
        let ep = endpoint as usize;
        if at[ep].endpoint == 0 {
            continue;
        }
        let fictpoint = i32::from(at[ep].endpoint) + num_vertices - 1;
        let fp = fictpoint as usize;
        if fictpoint >= p.max_vertices
            || num_edges >= p.max_edges
            || p.vert[fp].num_adj_edges >= p.vert[fp].max_adj_edges
            || p.vert[ep].num_adj_edges >= p.vert[ep].max_adj_edges
        {
            ret = BNS_VERT_EDGE_OVFL;
            break;
        }
        let mut eif = EndpointInfo::default();
        let pt = |f: u64, g: EndpointInfoFn, eif: &mut EndpointInfo| {
            flags & f != 0 && g(at, ep, eif) != 0
        };
        if n_get_endpoint_info(at, ep, &mut eif) == 0
            && !pt(TG_FLAG_PT_22_00, n_get_endpoint_info_pt_22_00, &mut eif)
            && !pt(TG_FLAG_PT_16_00, n_get_endpoint_info_pt_16_00, &mut eif)
            && !pt(TG_FLAG_PT_06_00, n_get_endpoint_info_pt_06_00, &mut eif)
            && !pt(TG_FLAG_PT_39_00, n_get_endpoint_info_pt_39_00, &mut eif)
            && !pt(TG_FLAG_PT_13_00, n_get_endpoint_info_pt_13_00, &mut eif)
            && !pt(TG_FLAG_PT_18_00, n_get_endpoint_info_pt_18_00, &mut eif)
            && !pt(TG_FLAG_KETO_ENOL_TAUT, n_get_endpoint_info_ket, &mut eif)
        {
            ret = BNS_BOND_ERR;
            break;
        }
        p.vert[ep].type_ |= BNS_VERT_TYPE_ENDPOINT;
        for k in 0..p.vert[ep].num_adj_edges as usize {
            let iedge = p.iedge[p.vert[ep].iedge + k] as usize;
            if p.edge[iedge].cap == 0 {
                let centerpoint = i32::from(p.edge[iedge].neighbor12) ^ endpoint;
                if centerpoint < p.num_atoms && p.vert[centerpoint as usize].st_edge.cap >= 1 {
                    let bond_type = i32::from(at[ep].bond_type[k]) & BOND_TYPE_MASK;
                    if bond_type == BOND_TAUTOM
                        || bond_type == BOND_ALTERN
                        || bond_type == BOND_ALT12NS
                        || bond_type == BOND_TYPE_SINGLE
                    {
                        p.edge[iedge].cap = 1;
                    }
                }
            }
        }
        let ne = num_edges as usize;
        p.edge[ne].cap = 1;
        p.edge[ne].flow = 0;
        p.edge[ne].pass = 0;
        p.edge[ne].forbidden &= p.edge_forbidden_mask;
        if eif.c_donor != 0 {
            p.edge[ne].flow += 1;
            p.vert[fp].st_edge.flow += 1;
            p.vert[fp].st_edge.cap += 1;
            p.vert[ep].st_edge.flow += 1;
            p.vert[ep].st_edge.cap += 1;
        }
        p.edge[ne].neighbor1 = endpoint as u16;
        p.edge[ne].neighbor12 = (endpoint ^ fictpoint) as u16;
        let o1 = p.vert[ep].num_adj_edges;
        let o2 = p.vert[fp].num_adj_edges;
        let (i1, i2) = (p.vert[ep].iedge, p.vert[fp].iedge);
        p.iedge[i1 + usize::from(o1)] = num_edges;
        p.iedge[i2 + usize::from(o2)] = num_edges;
        num_edges += 1;
        p.edge[ne].neigh_ord[0] = o1;
        p.edge[ne].neigh_ord[1] = o2;
        p.vert[ep].num_adj_edges += 1;
        p.vert[fp].num_adj_edges += 1;
        p.edge[ne].cap0 = p.edge[ne].cap;
        p.edge[ne].flow0 = p.edge[ne].flow;
    }
    p.num_edges = num_edges;
    p.num_vertices += n_max;
    p.num_t_groups = num_tg as i32;
    ret
}

/// `AddCGroups2BnStruct`: adds a fictitious vertex for each charge group and
/// joins the c-points to it
pub fn add_c_groups_to_bn_struct(
    p: &mut BnStruct,
    at: &[InpAtom],
    num_atoms: i32,
    cgi: Option<&mut CGroupInfo>,
) -> i32 {
    let mut ret = 0;
    let Some(cgi) = cgi else { return 0 };
    if cgi.num_c_groups == 0 || cgi.c_group.is_none() {
        return 0;
    }
    let num_cg = cgi.num_c_groups;
    let mut num_edges = p.num_edges;
    let num_vertices = p.num_vertices;
    if num_vertices + num_cg as i32 >= p.max_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    let cg = cgi.c_group.as_mut().expect("c_group");
    let n_max = cg[..num_cg]
        .iter()
        .map(|c| i32::from(c.n_group_number))
        .fold(0, i32::max);
    for v in &mut p.vert[num_vertices as usize..(num_vertices + n_max) as usize] {
        *v = BnsVertex::default();
    }
    if n_max != i32::from(cg[num_cg - 1].n_group_number) {
        cg[..num_cg].sort_by_key(|c| c.n_group_number);
    }
    let mut prev = (num_vertices - 1) as usize;
    for c in cg.iter().take(num_cg) {
        let cur = (num_vertices + i32::from(c.n_group_number) - 1) as usize;
        let iedge = p.vert[prev].iedge + usize::from(p.vert[prev].max_adj_edges);
        let v = &mut p.vert[cur];
        v.iedge = iedge;
        v.max_adj_edges = (i32::from(c.num_c_points) + BNS_ADD_EDGES) as u16;
        v.num_adj_edges = 0;
        v.st_edge.flow = 0;
        v.st_edge.flow0 = 0;
        v.st_edge.cap = 0;
        v.st_edge.cap0 = 0;
        v.type_ = BNS_VERT_TYPE_C_GROUP;
        prev = cur;
    }
    for c_point in 0..num_atoms {
        let cp = c_point as usize;
        if at[cp].c_point == 0 {
            continue;
        }
        let fictpoint = i32::from(at[cp].c_point) + num_vertices - 1;
        let fp = fictpoint as usize;
        if fictpoint >= p.max_vertices
            || num_edges >= p.max_edges
            || p.vert[fp].num_adj_edges >= p.vert[fp].max_adj_edges
            || p.vert[cp].num_adj_edges >= p.vert[cp].max_adj_edges
        {
            ret = BNS_VERT_EDGE_OVFL;
            break;
        }
        p.vert[cp].type_ |= BNS_VERT_TYPE_C_POINT;
        let ne = num_edges as usize;
        p.edge[ne].cap = 1;
        p.edge[ne].flow = 0;
        p.edge[ne].pass = 0;
        p.edge[ne].forbidden &= p.edge_forbidden_mask;
        if at[cp].charge != 1 {
            p.edge[ne].flow += 1;
            p.vert[fp].st_edge.flow += 1;
            p.vert[fp].st_edge.cap += 1;
            p.vert[cp].st_edge.flow += 1;
            p.vert[cp].st_edge.cap += 1;
        }
        for k in 0..p.vert[cp].num_adj_edges as usize {
            let iedge = p.iedge[p.vert[cp].iedge + k] as usize;
            let mut new_cap = p.vert[cp].st_edge.cap;
            let centerpoint = i32::from(p.edge[iedge].neighbor12) ^ c_point;
            if p.edge[iedge].cap == 0
                && centerpoint < p.num_atoms
                && p.vert[centerpoint as usize].st_edge.cap >= 1
            {
                new_cap = new_cap
                    .min(p.vert[centerpoint as usize].st_edge.cap)
                    .min(MAX_BOND_EDGE_CAP);
                p.edge[iedge].cap = new_cap;
            }
        }
        p.edge[ne].neighbor1 = c_point as u16;
        p.edge[ne].neighbor12 = (c_point ^ fictpoint) as u16;
        let o1 = p.vert[cp].num_adj_edges;
        let o2 = p.vert[fp].num_adj_edges;
        let (i1, i2) = (p.vert[cp].iedge, p.vert[fp].iedge);
        p.iedge[i1 + usize::from(o1)] = num_edges;
        p.iedge[i2 + usize::from(o2)] = num_edges;
        num_edges += 1;
        p.edge[ne].neigh_ord[0] = o1;
        p.edge[ne].neigh_ord[1] = o2;
        p.vert[cp].num_adj_edges += 1;
        p.vert[fp].num_adj_edges += 1;
        p.edge[ne].cap0 = p.edge[ne].cap;
        p.edge[ne].flow0 = p.edge[ne].flow;
    }
    p.num_edges = num_edges;
    p.num_vertices += n_max;
    p.num_c_groups = num_cg as i32;
    ret
}

/// `ReInitBnStructAddGroups`: strips the groups and adds the given ones
pub fn re_init_bn_struct_add_groups(
    p: &mut BnStruct,
    at: &mut [InpAtom],
    num_atoms: usize,
    tgi: Option<&mut TGroupInfo>,
    cgi: Option<&mut CGroupInfo>,
) -> i32 {
    let ret = re_init_bn_struct(p, at, num_atoms, false);
    if ret != 0 {
        return super::ichierr::BNS_REINIT_ERR;
    }
    if u64::from(p.pb_taut_flags) & TG_FLAG_MOVE_POS_CHARGES != 0 {
        let r = add_c_groups_to_bn_struct(p, at, num_atoms as i32, cgi);
        if is_bns_error(r) {
            return r;
        }
    }
    add_t_groups_to_bn_struct(p, at, num_atoms as i32, tgi)
}

// ---- creating and removing single groups ----

pub const KNOWN_ACIDIC_TYPE: i32 =
    (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO) as i32;

/// `CreateCGroupInBnStruct`: adds one c-group of the atoms of type `n_type` and
/// mask `n_mask`; returns the new vertex number, 0 if there are no atoms
/// or an error
pub fn create_c_group_in_bn_struct(
    at: &[InpAtom],
    num_atoms: i32,
    p: &mut BnStruct,
    n_type: i32,
    n_mask: i32,
    n_charge: i32,
) -> i32 {
    let mut num_edges = p.num_edges;
    let num_vertices = p.num_vertices;
    if num_vertices + 1 >= p.max_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    let selected = |i: usize| -> Option<i32> {
        let (t, m) = get_atom_charge_type(at, i, None, 0);
        (n_type & t != 0 && m.unwrap_or(0) & n_mask != 0).then_some(t)
    };
    let num_c_points = (0..num_atoms as usize)
        .filter(|&i| selected(i).is_some())
        .count() as i32;
    if num_c_points == 0 {
        return 0;
    }
    let prev = (num_vertices - 1) as usize;
    let cur = num_vertices as usize;
    let iedge = p.vert[prev].iedge + usize::from(p.vert[prev].max_adj_edges);
    p.vert[cur] = BnsVertex {
        iedge,
        max_adj_edges: (num_c_points + BNS_ADD_EDGES) as u16,
        type_: BNS_VERT_TYPE_C_GROUP
            | if n_charge < 0 {
                BNS_VERT_TYPE_C_NEGATIVE
            } else {
                0
            },
        ..Default::default()
    };
    let fictpoint = 1 + num_vertices - 1;
    let fp = fictpoint as usize;
    for c_point in 0..num_atoms {
        let cp = c_point as usize;
        let Some(type_) = selected(cp) else { continue };
        if fictpoint >= p.max_vertices
            || num_edges >= p.max_edges
            || p.vert[fp].num_adj_edges >= p.vert[fp].max_adj_edges
            || p.vert[cp].num_adj_edges >= p.vert[cp].max_adj_edges
        {
            break;
        }
        p.vert[cp].type_ |= BNS_VERT_TYPE_C_POINT;
        if KNOWN_ACIDIC_TYPE & type_ != 0 && n_charge < 0 {
            p.vert[cp].type_ |= p.type_tacn;
        }
        let ne = num_edges as usize;
        p.edge[ne].cap = 1;
        p.edge[ne].flow = 0;
        p.edge[ne].pass = 0;
        p.edge[ne].forbidden &= p.edge_forbidden_mask;
        let charge = i32::from(at[cp].charge);
        if (n_charge == 1 && charge != 1) || (n_charge == -1 && charge == -1) {
            p.edge[ne].flow += 1;
            p.vert[fp].st_edge.flow += 1;
            p.vert[fp].st_edge.cap += 1;
            p.vert[cp].st_edge.flow += 1;
            p.vert[cp].st_edge.cap += 1;
        }
        for k in 0..p.vert[cp].num_adj_edges as usize {
            let iedge = p.iedge[p.vert[cp].iedge + k] as usize;
            let mut new_cap = p.vert[cp].st_edge.cap;
            let centerpoint = i32::from(p.edge[iedge].neighbor12) ^ c_point;
            if p.edge[iedge].cap == 0
                && centerpoint < p.num_atoms
                && p.vert[centerpoint as usize].st_edge.cap >= 1
            {
                new_cap = new_cap
                    .min(p.vert[centerpoint as usize].st_edge.cap)
                    .min(MAX_BOND_EDGE_CAP);
                p.edge[iedge].cap = new_cap;
            }
        }
        p.edge[ne].neighbor1 = c_point as u16;
        p.edge[ne].neighbor12 = (c_point ^ fictpoint) as u16;
        let o1 = p.vert[cp].num_adj_edges;
        let o2 = p.vert[fp].num_adj_edges;
        let (i1, i2) = (p.vert[cp].iedge, p.vert[fp].iedge);
        p.iedge[i1 + usize::from(o1)] = num_edges;
        p.iedge[i2 + usize::from(o2)] = num_edges;
        num_edges += 1;
        p.edge[ne].neigh_ord[0] = o1;
        p.edge[ne].neigh_ord[1] = o2;
        p.vert[cp].num_adj_edges += 1;
        p.vert[fp].num_adj_edges += 1;
        p.edge[ne].cap0 = p.edge[ne].cap;
        p.edge[ne].flow0 = p.edge[ne].flow;
    }
    p.num_edges = num_edges;
    p.num_vertices += 1;
    p.num_c_groups += 1;
    num_vertices
}

/// `CreateTGroupInBnStruct`: adds one t-group of the atoms of type `n_type` and
/// mask `n_mask`; returns the new vertex number, 0 if there are no atoms
/// or an error
pub fn create_t_group_in_bn_struct(
    at: &[InpAtom],
    num_atoms: i32,
    p: &mut BnStruct,
    n_type: i32,
    n_mask: i32,
) -> i32 {
    let mut num_edges = p.num_edges;
    let num_vertices = p.num_vertices;
    if num_vertices + 1 >= p.max_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    let selected = |i: usize| -> bool {
        let (t, m) = get_atom_charge_type(at, i, None, 0);
        n_type & t != 0 && m.unwrap_or(0) & n_mask != 0
    };
    let num_endpoints = (0..num_atoms as usize).filter(|&i| selected(i)).count() as i32;
    if num_endpoints == 0 {
        return 0;
    }
    p.vert[num_vertices as usize] = BnsVertex::default();
    let prev = (num_vertices - 1) as usize;
    let cur = num_vertices as usize;
    let iedge = p.vert[prev].iedge + usize::from(p.vert[prev].max_adj_edges);
    p.vert[cur].iedge = iedge;
    p.vert[cur].max_adj_edges = (num_endpoints + BNS_ADD_EDGES + BNS_ADD_SUPER_TGROUP) as u16;
    p.vert[cur].type_ |= BNS_VERT_TYPE_TGROUP;
    let fictpoint = 1 + num_vertices - 1;
    let fp = fictpoint as usize;
    for endpoint in 0..num_atoms {
        let ep = endpoint as usize;
        if !selected(ep) {
            continue;
        }
        if fictpoint >= p.max_vertices
            || num_edges >= p.max_edges
            || p.vert[fp].num_adj_edges >= p.vert[fp].max_adj_edges
            || p.vert[ep].num_adj_edges >= p.vert[ep].max_adj_edges
        {
            break;
        }
        let a = &at[ep];
        let neutral_valence =
            i32::from(a.chem_bonds_valence) + i32::from(a.num_h) - i32::from(a.charge);
        if neutral_valence != 2 && neutral_valence != 3 {
            break;
        }
        let num_bonds = i32::from(a.valence);
        let mut edge_cap = neutral_valence - num_bonds;
        if neutral_valence == 3 && num_bonds > 1 {
            edge_cap += 1;
        }
        let edge_flow = i32::from(a.num_h).min(edge_cap);
        p.vert[ep].type_ |= BNS_VERT_TYPE_ENDPOINT;
        let ne = num_edges as usize;
        p.edge[ne].cap = edge_cap;
        p.edge[ne].flow = edge_flow;
        p.edge[ne].pass = 0;
        p.edge[ne].forbidden &= p.edge_forbidden_mask;
        p.vert[fp].st_edge.flow += edge_flow;
        p.vert[fp].st_edge.cap += edge_flow;
        p.vert[ep].st_edge.flow += edge_flow;
        p.vert[ep].st_edge.cap += edge_flow;
        for k in 0..p.vert[ep].num_adj_edges as usize {
            let iedge = p.iedge[p.vert[ep].iedge + k] as usize;
            let mut new_cap = p.vert[ep].st_edge.cap;
            if p.edge[iedge].cap == 0 {
                let centerpoint = i32::from(p.edge[iedge].neighbor12) ^ endpoint;
                if centerpoint < p.num_atoms && p.vert[centerpoint as usize].st_edge.cap >= 1 {
                    new_cap = new_cap
                        .min(p.vert[centerpoint as usize].st_edge.cap)
                        .min(MAX_BOND_EDGE_CAP);
                    p.edge[iedge].cap = new_cap;
                }
            }
        }
        p.edge[ne].neighbor1 = endpoint as u16;
        p.edge[ne].neighbor12 = (endpoint ^ fictpoint) as u16;
        let o1 = p.vert[ep].num_adj_edges;
        let o2 = p.vert[fp].num_adj_edges;
        let (i1, i2) = (p.vert[ep].iedge, p.vert[fp].iedge);
        p.iedge[i1 + usize::from(o1)] = num_edges;
        p.iedge[i2 + usize::from(o2)] = num_edges;
        num_edges += 1;
        p.edge[ne].neigh_ord[0] = o1;
        p.edge[ne].neigh_ord[1] = o2;
        p.vert[ep].num_adj_edges += 1;
        p.vert[fp].num_adj_edges += 1;
        p.edge[ne].cap0 = p.edge[ne].cap;
        p.edge[ne].flow0 = p.edge[ne].flow;
    }
    p.num_edges = num_edges;
    p.num_vertices += 1;
    p.num_t_groups += 1;
    num_vertices
}

/// `RemoveLastGroupFromBnStruct`: removes the last t-group or c-group vertex
pub fn remove_last_group_from_bn_struct(
    at: &mut [InpAtom],
    num_atoms: i32,
    tg: i32,
    p: &mut BnStruct,
) -> i32 {
    let mut num_edges = p.num_edges;
    let num_vertices = p.num_vertices;
    if p.num_added_atoms + p.num_c_groups + p.num_t_groups + num_atoms >= p.max_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    if tg + 1 != num_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    let fp = tg as usize;
    let ftype = p.vert[fp].type_;
    let is_t_group = ftype & BNS_VERT_TYPE_TGROUP != 0;
    let mut is_c_group = 0;
    if ftype & BNS_VERT_TYPE_C_GROUP != 0 {
        is_c_group = if ftype & BNS_VERT_TYPE_C_NEGATIVE != 0 {
            2
        } else {
            1
        };
    }
    let mut k = i32::from(p.vert[fp].num_adj_edges) - 1;
    while k >= 0 {
        let iedge = p.iedge[p.vert[fp].iedge + k as usize];
        if iedge + 1 != num_edges {
            return BNS_VERT_EDGE_OVFL;
        }
        let e = p.edge[iedge as usize];
        let endpoint = i32::from(e.neighbor12) ^ tg;
        let ve = &mut p.vert[endpoint as usize];
        ve.st_edge.cap -= e.flow;
        ve.st_edge.cap0 = ve.st_edge.cap;
        ve.st_edge.flow -= e.flow;
        ve.st_edge.flow0 = ve.st_edge.flow;
        if p.type_tacn != 0 && ve.type_ & p.type_tacn == p.type_tacn {
            ve.type_ ^= p.type_tacn;
        }
        if is_t_group {
            ve.type_ ^= ftype & BNS_VERT_TYPE_ENDPOINT;
        }
        if is_c_group != 0 {
            ve.type_ ^= ftype & BNS_VERT_TYPE_C_POINT;
        }
        if i32::from(e.neigh_ord[0]) + 1 != i32::from(ve.num_adj_edges) {
            return BNS_VERT_EDGE_OVFL;
        }
        ve.num_adj_edges -= 1;
        p.edge[iedge as usize] = BnsEdge::default();
        num_edges -= 1;
        if is_t_group && endpoint < num_atoms {
            at[0].endpoint = 0;
        }
        if is_c_group == 1 && endpoint < num_atoms {
            at[0].c_point = 0;
        }
        k -= 1;
    }
    p.vert[fp] = BnsVertex::default();
    p.num_edges = num_edges;
    p.num_vertices = num_vertices - 1;
    if is_t_group {
        p.num_t_groups -= 1;
    }
    if is_c_group != 0 {
        p.num_c_groups -= 1;
    }
    0
}

/// `SetInitCapFlowToCurrent`
pub fn set_init_cap_flow_to_current(p: &mut BnStruct) -> i32 {
    for i in 0..p.num_vertices as usize {
        let st = &mut p.vert[i].st_edge;
        st.flow0 = st.flow;
        st.cap0 = st.cap;
        let v = p.vert[i];
        for j in 0..v.num_adj_edges as usize {
            let e = &mut p.edge[p.iedge[v.iedge + j] as usize];
            e.cap0 = e.cap;
            e.flow0 = e.flow;
        }
    }
    0
}

// ---- simple proton addition and removal ----

use super::ichitaut::{SALT_ACCEPTOR, SALT_DONOR_H, SALT_DONOR_NEG};

const PR_SIMPLE_MSK: i32 = atbit::NP_PROTON | atbit::OH_PLUS;
const PR_SIMPLE_TYP: i32 = (ATT_ATOM_N | ATT_ATOM_P | ATT_O_PLUS) as i32;
const ATBIT_MSK_NP: i32 =
    atbit::NP_PLUS | atbit::NP_PROTON | atbit::NP_H | atbit::N_MINUS | atbit::NP;
const ATBIT_MSK_OS: i32 = atbit::COH
    | atbit::CSH
    | atbit::ZOH
    | atbit::OOH
    | atbit::ZOOH
    | atbit::NOH
    | atbit::N_OH
    | atbit::CO
    | atbit::ZO
    | atbit::NO
    | atbit::N_O
    | atbit::CO_MINUS
    | atbit::CS_MINUS
    | atbit::ZO_MINUS
    | atbit::OO_MINUS
    | atbit::ZOO_MINUS
    | atbit::NO_MINUS
    | atbit::N_O_MINUS;
const ATTYP_OS: i32 =
    (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO | ATT_OTHER_ZO) as i32;
const ATTYP_N: i32 = ATT_ATOM_N as i32;
const PR_HARD_MSK_NEG: i32 = ATBIT_MSK_NP | ATBIT_MSK_OS;
const PR_HARD_TYP_NEG: i32 = ATTYP_N | ATTYP_OS;
const PR_HARD_MSK_H: i32 = ATBIT_MSK_NP | ATBIT_MSK_OS;
const PR_HARD_TYP_H: i32 = ATTYP_N | ATTYP_OS;

/// `ArTypMask` (type, mask) pairs of the acidic H groups (`AR_SIMPLE_*`;
/// the pairs 2 and 3 are 0 because `AR_ANY_OH` is 0)
const AR_TYP_MASK: [(i32, i32); 1] = [(
    (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO | ATT_HAL_ACID) as i32,
    atbit::COH | atbit::CSH | atbit::OOH | atbit::ZOOH | atbit::NOH | atbit::HAL_ACID,
)];

/// `AaTypMask` (type, mask) pairs of the acidic anions (`AA_SIMPLE_*`; the
/// pairs 2 and 3 are 0 because `AA_ANY_O_Minus` is 0)
const AA_TYP_MASK: [(i32, i32); 2] = [
    (
        (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO | ATT_OH_MINUS | ATT_HAL_ANION)
            as i32,
        atbit::CO_MINUS
            | atbit::CS_MINUS
            | atbit::OO_MINUS
            | atbit::ZOO_MINUS
            | atbit::NO_MINUS
            | atbit::O_MINUS
            | atbit::HAL_ANION,
    ),
    (ATT_NP_MINUS_V23 as i32, atbit::N_MINUS),
];

/// `SimpleRemoveHplusNPO`: removes H(+) from the N, P and O atoms that carry
/// it; returns the number of removed protons
pub fn simple_remove_hplus_npo(
    at: &mut [InpAtom],
    num_atoms: i32,
    mut totals: Option<&mut [i32]>,
    tgi: &mut TGroupInfo,
) -> i32 {
    let mut num_removed = 0;
    for i in 0..num_atoms as usize {
        let (type_, mask) = get_atom_charge_type(at, i, None, 0);
        if PR_SIMPLE_TYP & type_ != 0 && PR_SIMPLE_MSK & mask.unwrap_or(0) != 0 {
            get_atom_charge_type(at, i, totals.as_deref_mut(), 1);
            at[i].charge = 0;
            add_or_remove_expl_or_impl_h(-1, at, num_atoms, i, tgi);
            num_removed += 1;
            get_atom_charge_type(at, i, totals.as_deref_mut(), 0);
        }
    }
    num_removed
}

/// `bIsAtomTypeHard` (`OPPOSITE_CHARGE_IN_CGROUP` is 1: the charge is not
/// checked)
pub fn b_is_atom_type_hard(at: &[InpAtom], endpoint: usize, n_type: i32, n_mask: i32) -> bool {
    let (type_, mask) = get_atom_charge_type(at, endpoint, None, 0);
    n_type & type_ != 0 && mask.unwrap_or(0) & n_mask != 0
}

/// `bIsHDonorAccAtomType`
pub fn b_is_h_donor_acc_atom_type(at: &[InpAtom], endpoint: usize, c_sub_type: &mut i32) -> i32 {
    if b_is_atom_type_hard(at, endpoint, PR_HARD_TYP_H, PR_HARD_MSK_H) {
        let a = &at[endpoint];
        let neutral_valence =
            i32::from(a.chem_bonds_valence) + i32::from(a.num_h) - i32::from(a.charge);
        if neutral_valence != 2 && neutral_valence != 3 {
            return -1;
        }
        let mut edge_flow = i32::from(a.num_h);
        let edge_cap = neutral_valence - i32::from(a.valence);
        edge_flow = edge_flow.min(edge_cap);
        if edge_cap != 0 {
            if edge_cap > edge_flow {
                *c_sub_type |= SALT_ACCEPTOR;
            }
            if edge_flow != 0 {
                *c_sub_type |= SALT_DONOR_H;
            }
            return 4;
        }
    }
    -1
}

/// `bIsNegAtomType`
pub fn b_is_neg_atom_type(at: &[InpAtom], endpoint: usize, c_sub_type: &mut i32) -> i32 {
    let mut sub_type = 0;
    if b_is_atom_type_hard(at, endpoint, PR_HARD_TYP_NEG, PR_HARD_MSK_NEG) {
        let a = &at[endpoint];
        let neutral_valence =
            i32::from(a.chem_bonds_valence) + i32::from(a.num_h) - i32::from(a.charge);
        if neutral_valence != 2 && neutral_valence != 3 {
            return -1;
        }
        let mut edge_flow = i32::from(a.charge == -1);
        let edge_cap = neutral_valence - i32::from(a.valence) - i32::from(a.num_h);
        edge_flow = edge_flow.min(edge_cap);
        if edge_cap != 0 {
            if edge_cap > edge_flow {
                sub_type |= SALT_ACCEPTOR;
            }
            if edge_flow != 0 {
                sub_type |= SALT_DONOR_NEG;
            }
            if sub_type != 0 {
                *c_sub_type |= sub_type;
                return 4;
            }
        }
    }
    -1
}

/// `bIsHardRemHCandidate`
pub fn b_is_hard_rem_h_candidate(at: &[InpAtom], i: usize, c_sub_type: &mut i32) -> i32 {
    let mut sub_type = 0;
    let ret1 = b_is_h_donor_acc_atom_type(at, i, &mut sub_type);
    let ret2 = b_is_neg_atom_type(at, i, &mut sub_type);
    let ret = ret1.max(ret2);
    if ret > 0 && sub_type != 0 {
        *c_sub_type |= sub_type;
        return ret;
    }
    -1
}

/// `SimpleRemoveAcidicProtons` / `SimpleAddAcidicProtons`: the first
/// `num` atoms of the earliest matching class lose (`remove`) a proton or
/// get one; the C code tests the mask with `&&`, i.e. for being non-zero
fn simple_move_acidic_protons(
    at: &mut [InpAtom],
    num_atoms: i32,
    aatg: &mut BnAatg,
    num: i32,
    classes: &[(i32, i32)],
    remove: bool,
) -> i32 {
    let selected = |a: &InpAtom| {
        if remove {
            a.charge == 0 && a.num_h != 0
        } else {
            a.charge == -1
        }
    };
    let mut count = vec![0i32; classes.len()];
    let mut max_j = classes.len() - 1;
    for i in 0..num_atoms as usize {
        if selected(&at[i]) {
            let (type_, mask) = get_atom_charge_type(at, i, None, 0);
            if type_ != 0 {
                for (j, c) in classes.iter().enumerate().take(max_j + 1) {
                    if type_ & c.0 != 0 && mask.unwrap_or(0) != 0 {
                        count[j] += 1;
                        break;
                    }
                }
            }
        }
    }
    let mut num_tot = 0;
    for j in 0..=max_j {
        num_tot += count[j];
        if num_tot >= num {
            max_j = j;
            break;
        }
    }
    if num_tot == 0 {
        return 0;
    }
    let mut num_done = 0;
    let mut i = 0usize;
    while i < num_atoms as usize && num_done < num {
        if selected(&at[i]) {
            let (type_, mask) = get_atom_charge_type(at, i, None, 0);
            if type_ != 0 {
                for j in 0..=max_j {
                    if count[j] != 0 && type_ & classes[j].0 != 0 && mask.unwrap_or(0) != 0 {
                        get_atom_charge_type(at, i, aatg.n_at_type_totals.as_deref_mut(), 1);
                        count[j] -= 1;
                        at[i].charge += if remove { -1 } else { 1 };
                        let tgi = aatg.t_group_info.as_mut().expect("t_group_info");
                        add_or_remove_expl_or_impl_h(
                            if remove { -1 } else { 1 },
                            at,
                            num_atoms,
                            i,
                            tgi,
                        );
                        num_done += 1;
                        get_atom_charge_type(at, i, aatg.n_at_type_totals.as_deref_mut(), 0);
                        break;
                    }
                }
            }
        }
        i += 1;
    }
    num_done
}

/// `SimpleRemoveAcidicProtons`: removes up to `num2remove` protons from acids
pub fn simple_remove_acidic_protons(
    at: &mut [InpAtom],
    num_atoms: i32,
    aatg: &mut BnAatg,
    num2remove: i32,
) -> i32 {
    simple_move_acidic_protons(at, num_atoms, aatg, num2remove, &AR_TYP_MASK, true)
}

/// `SimpleAddAcidicProtons`: adds up to `num2add` protons to acidic anions
pub fn simple_add_acidic_protons(
    at: &mut [InpAtom],
    num_atoms: i32,
    aatg: &mut BnAatg,
    num2add: i32,
) -> i32 {
    simple_move_acidic_protons(at, num_atoms, aatg, num2add, &AA_TYP_MASK, false)
}

fn has_class(at: &[InpAtom], i: usize, classes: &[(i32, i32)]) -> i32 {
    let (type_, mask) = get_atom_charge_type(at, i, None, 0);
    if type_ != 0 {
        let mask = mask.unwrap_or(0);
        for c in classes {
            if type_ & c.0 != 0 && mask & c.1 != 0 {
                return 1;
            }
        }
    }
    0
}

/// `bHasAcidicHydrogen`
pub fn b_has_acidic_hydrogen(at: &[InpAtom], i: usize) -> i32 {
    if at[i].charge == 0 && at[i].num_h != 0 {
        has_class(at, i, &AR_TYP_MASK)
    } else {
        0
    }
}

/// `bHasOtherExchangableH`
pub fn b_has_other_exchangable_h(at: &[InpAtom], i: usize) -> i32 {
    if at[i].num_h != 0 {
        let (type_, mask) = get_atom_charge_type(at, i, None, 0);
        if type_ != 0 && type_ & i32::from(ATT_ATOM_N) != 0 && mask.unwrap_or(0) & atbit::NP_H != 0
        {
            return 1;
        }
    }
    0
}

/// `bHasAcidicMinus`
pub fn b_has_acidic_minus(at: &[InpAtom], i: usize) -> i32 {
    if at[i].charge == -1 {
        has_class(at, i, &AA_TYP_MASK)
    } else {
        0
    }
}

// ---- hard proton addition and removal ----

use super::ichierr::BNS_RADICAL_ERR;
use super::ichiparm::{
    TG_FLAG_HARD_ADD_REM_PROTONS, TG_FLAG_TEST_TAUT__SALTS,
    TG_FLAG_TEST_TAUT2_SALTS as TG_FLAG_TEST_TAUT2_SALTS64,
};

pub const FLAG_PROTON_NPO_SIMPLE_REMOVED: u64 = 0x0001;
pub const FLAG_PROTON_NP_HARD_REMOVED: u64 = 0x0002;
pub const FLAG_PROTON_AC_SIMPLE_ADDED: u64 = 0x0004;
pub const FLAG_PROTON_AC_HARD_REMOVED: u64 = 0x0010;
pub const FLAG_PROTON_AC_HARD_ADDED: u64 = 0x0020;
pub const FLAG_PROTON_AC_SIMPLE_REMOVED: u64 = 0x0008;
pub const FLAG_PROTON_CHARGE_CANCEL: u64 = 0x0040;

const HARD_TYP_POS: i32 = ATTYP_N;
const HARD_MSK_POS: i32 = ATBIT_MSK_NP;
const AR_HARD_MSK_HA: i32 = atbit::CO | atbit::NO;
const AR_HARD_TYP_HA: i32 = (ATT_ACIDIC_CO | ATT_NO) as i32;
const AR_HARD_MSK_HN: i32 = (ATBIT_MSK_NP | ATBIT_MSK_OS) & !AR_HARD_MSK_HA;
const AR_HARD_TYP_HN: i32 = ATTYP_N | ATTYP_OS;
const AA_HARD_MSK_NEG: i32 = (ATBIT_MSK_NP | ATBIT_MSK_OS) & !(atbit::CO | atbit::NO);

fn totals(aatg: &BnAatg) -> &[i32] {
    aatg.n_at_type_totals.as_deref().expect("nAtTypeTotals")
}

/// charges (positive, negative) from the totals
fn charge_counts(aatg: &BnAatg) -> (i32, i32) {
    let t = totals(aatg);
    (
        (t[attot::NUM_CHARGES] + t[attot::TOT_CHARGE]) / 2,
        (t[attot::NUM_CHARGES] - t[attot::TOT_CHARGE]) / 2,
    )
}

fn set_group_types(p: &mut BnStruct) {
    p.type_cn = BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE;
    p.type_t = BNS_VERT_TYPE_TGROUP;
    p.type_tacn = BNS_VERT_TYPE_ACID;
}

/// removes the groups (the first of `groups` that were created) in order;
/// returns the first error
fn remove_groups(at: &mut [InpAtom], num_atoms: i32, p: &mut BnStruct, groups: &[i32]) -> i32 {
    let mut ret = 0;
    for &g in groups {
        if g >= num_atoms {
            let ret2 = remove_last_group_from_bn_struct(at, num_atoms, g, p);
            if ret == 0 && ret2 != 0 {
                ret = ret2;
            }
        }
    }
    ret
}

/// the repeated alternating-path search of the hard proton functions: runs
/// `bExistsAltPath` until no path is found (or `limit` paths were found);
/// returns the number of found paths, `Err` on a BNS error. `neutralize`
/// counts neutralized ion pairs the way the C code does for the two kinds of
/// search.
#[allow(clippy::too_many_arguments)]
fn hard_alt_path_loop(
    p: &mut BnStruct,
    bd: &mut BnData,
    aatg: &mut BnAatg,
    at: &mut [InpAtom],
    num_atoms: i32,
    v_double: i32,
    v_single: i32,
    limit: Option<i32>,
    count_proton: bool,
    check_removed_protons: bool,
    n_neutralized: &mut i32,
) -> Result<i32, i32> {
    let mut found = 0;
    loop {
        let prev_charges = totals(aatg)[attot::NUM_CHARGES];
        let prev_removed = aatg
            .t_group_info
            .as_ref()
            .map_or(0, |t| t.tni.n_num_removed_protons);
        let ret = exists_alt_path(
            p,
            bd,
            Some(&mut *aatg),
            at,
            num_atoms,
            v_double,
            v_single,
            ALT_PATH_MODE_REM_PROTON,
        );
        if is_bns_error(ret) {
            return Err(ret);
        }
        if check_removed_protons {
            let curr = aatg
                .t_group_info
                .as_ref()
                .map_or(0, |t| t.tni.n_num_removed_protons);
            let want = if count_proton {
                prev_removed + (ret & 1)
            } else {
                prev_removed
            };
            if curr != want {
                return Err(BNS_RADICAL_ERR);
            }
        }
        if ret & 1 != 0 {
            let cur = totals(aatg)[attot::NUM_CHARGES];
            if count_proton {
                found += 1;
                if prev_charges > cur + 1 {
                    *n_neutralized += (prev_charges - (cur - 1)) / 2;
                }
            } else if prev_charges > cur {
                *n_neutralized += (prev_charges - cur) / 2;
            }
        }
        let go_on = ret & 1 != 0 && limit.is_none_or(|l| found < l);
        if !go_on {
            return Ok(found);
        }
    }
}

/// the neutralizing pass of the hard proton functions: ion pairs like
/// >N(+)=-O(-) => >N-=O are removed along (-) -> (+) paths
#[allow(clippy::too_many_arguments)]
fn hard_neutralize_loop(
    p: &mut BnStruct,
    bd: &mut BnData,
    aatg: &mut BnAatg,
    at: &mut [InpAtom],
    num_atoms: i32,
    cg_minus: i32,
    cg_plus: i32,
    check_removed_protons: bool,
    n_neutralized: &mut i32,
) -> Result<(), i32> {
    hard_alt_path_loop(
        p,
        bd,
        aatg,
        at,
        num_atoms,
        cg_minus,
        cg_plus,
        None,
        false,
        check_removed_protons,
        n_neutralized,
    )
    .map(|_| ())
}

fn finish_hard(
    at: &mut [InpAtom],
    num_atoms: i32,
    p: &mut BnStruct,
    groups: &[i32],
    aatg: &BnAatg,
    (pos, neg): (i32, i32),
    expected_diff: i32,
    n_canceled: Option<&mut i32>,
    n_neutralized: i32,
) -> Option<i32> {
    let ret = remove_groups(at, num_atoms, p, groups);
    p.type_cn = 0;
    p.type_t = 0;
    p.type_tacn = 0;
    if ret != 0 {
        return Some(ret);
    }
    let (pos2, neg2) = charge_counts(aatg);
    if (pos - neg) - (pos2 - neg2) != expected_diff {
        return Some(BNS_PROGRAM_ERR);
    }
    if let Some(c) = n_canceled {
        // FIX_CANCEL_CHARGE_COUNT_BUG is 0: the count is replaced
        *c = 2 * n_neutralized;
    }
    None
}

/// `HardRemoveAcidicProtons`: moves up to `num2remove` protons to acidic
/// groups along alternating paths; returns the number moved
#[allow(clippy::too_many_arguments)]
pub fn hard_remove_acidic_protons(
    at: &mut [InpAtom],
    num_atoms: i32,
    aatg: &mut BnAatg,
    num2remove: i32,
    n_canceled: Option<&mut i32>,
    p: &mut BnStruct,
    bd: &mut BnData,
) -> i32 {
    let counts = charge_counts(aatg);
    set_group_types(p);
    let cg_plus = create_c_group_in_bn_struct(at, num_atoms, p, ATTYP_N, ATBIT_MSK_NP, 1);
    let cg_minus =
        create_c_group_in_bn_struct(at, num_atoms, p, PR_HARD_TYP_NEG, PR_HARD_MSK_NEG, -1);
    set_group_types(p);
    let tg_other = create_t_group_in_bn_struct(at, num_atoms, p, AR_HARD_TYP_HN, AR_HARD_MSK_HN);
    let tg_acid = create_t_group_in_bn_struct(at, num_atoms, p, AR_HARD_TYP_HA, AR_HARD_MSK_HA);
    let mut n_neutralized = 0;
    let mut moved = 0;
    if tg_other >= num_atoms && tg_acid >= num_atoms {
        match hard_alt_path_loop(
            p,
            bd,
            aatg,
            at,
            num_atoms,
            tg_other,
            tg_acid,
            Some(num2remove),
            true,
            false,
            &mut n_neutralized,
        ) {
            Ok(n) => moved = n,
            Err(e) => return e,
        }
        if moved != 0
            && cg_minus >= num_atoms
            && cg_plus >= num_atoms
            && totals(aatg)[attot::NUM_CHARGES] > totals(aatg)[attot::TOT_CHARGE].abs()
        {
            if let Err(e) = hard_neutralize_loop(
                p,
                bd,
                aatg,
                at,
                num_atoms,
                cg_minus,
                cg_plus,
                false,
                &mut n_neutralized,
            ) {
                return e;
            }
        }
    }
    match finish_hard(
        at,
        num_atoms,
        p,
        &[tg_acid, tg_other, cg_minus, cg_plus],
        aatg,
        counts,
        0,
        n_canceled,
        n_neutralized,
    ) {
        Some(e) => e,
        None => moved,
    }
}

/// `HardAddAcidicProtons`: moves up to `num2add` protons to acidic anions
#[allow(clippy::too_many_arguments)]
pub fn hard_add_acidic_protons(
    at: &mut [InpAtom],
    num_atoms: i32,
    aatg: &mut BnAatg,
    num2add: i32,
    n_canceled: Option<&mut i32>,
    p: &mut BnStruct,
    bd: &mut BnData,
) -> i32 {
    let counts = charge_counts(aatg);
    set_group_types(p);
    let cg_plus = create_c_group_in_bn_struct(at, num_atoms, p, ATTYP_N, ATBIT_MSK_NP, 1);
    let cg_minus_co = create_c_group_in_bn_struct(
        at,
        num_atoms,
        p,
        (ATT_ACIDIC_CO | ATT_NO) as i32,
        atbit::CO | atbit::NO,
        -1,
    );
    let cg_minus_other =
        create_c_group_in_bn_struct(at, num_atoms, p, ATTYP_N | ATTYP_OS, AA_HARD_MSK_NEG, -1);
    set_group_types(p);
    let tg_h = create_t_group_in_bn_struct(
        at,
        num_atoms,
        p,
        ATTYP_N | ATTYP_OS,
        ATBIT_MSK_NP | ATBIT_MSK_OS,
    );
    let mut n_neutralized = 0;
    let mut moved = 0;
    if cg_minus_other >= num_atoms && cg_minus_co >= num_atoms {
        match hard_alt_path_loop(
            p,
            bd,
            aatg,
            at,
            num_atoms,
            cg_minus_other,
            cg_minus_co,
            Some(num2add),
            true,
            false,
            &mut n_neutralized,
        ) {
            Ok(n) => moved = n,
            Err(e) => return e,
        }
        if moved != 0
            && cg_minus_other >= num_atoms
            && cg_plus >= num_atoms
            && totals(aatg)[attot::NUM_CHARGES] > totals(aatg)[attot::TOT_CHARGE].abs()
        {
            if let Err(e) = hard_neutralize_loop(
                p,
                bd,
                aatg,
                at,
                num_atoms,
                cg_minus_other,
                cg_plus,
                false,
                &mut n_neutralized,
            ) {
                return e;
            }
        }
    }
    match finish_hard(
        at,
        num_atoms,
        p,
        &[tg_h, cg_minus_other, cg_minus_co, cg_plus],
        aatg,
        counts,
        0,
        n_canceled,
        n_neutralized,
    ) {
        Some(e) => e,
        None => moved,
    }
}

/// `HardRemoveHplusNP`: removes H(+) from cationic N/P along alternating
/// paths to a t-group; returns the number of removed protons
#[allow(clippy::too_many_arguments)]
pub fn hard_remove_hplus_np(
    at: &mut [InpAtom],
    num_atoms: i32,
    b_cancel_charges_always: bool,
    n_canceled: Option<&mut i32>,
    aatg: &mut BnAatg,
    p: &mut BnStruct,
    bd: &mut BnData,
) -> i32 {
    let counts = charge_counts(aatg);
    set_group_types(p);
    let cg_plus = create_c_group_in_bn_struct(at, num_atoms, p, HARD_TYP_POS, HARD_MSK_POS, 1);
    let cg_minus =
        create_c_group_in_bn_struct(at, num_atoms, p, PR_HARD_TYP_NEG, PR_HARD_MSK_NEG, -1);
    let tg_h = create_t_group_in_bn_struct(at, num_atoms, p, PR_HARD_TYP_H, PR_HARD_MSK_H);
    let mut n_neutralized = 0;
    let mut removed = 0;
    if tg_h >= num_atoms && cg_plus >= num_atoms {
        // FIX_N_MINUS_NORN_BUG is 0: no neutralizing pass before the removal
        match hard_alt_path_loop(
            p,
            bd,
            aatg,
            at,
            num_atoms,
            tg_h,
            cg_plus,
            None,
            true,
            true,
            &mut n_neutralized,
        ) {
            Ok(n) => removed = n,
            Err(e) => return e,
        }
        if (removed != 0 || b_cancel_charges_always)
            && cg_minus >= num_atoms
            && cg_plus >= num_atoms
            && totals(aatg)[attot::NUM_CHARGES] > totals(aatg)[attot::TOT_CHARGE].abs()
        {
            if let Err(e) = hard_neutralize_loop(
                p,
                bd,
                aatg,
                at,
                num_atoms,
                cg_minus,
                cg_plus,
                true,
                &mut n_neutralized,
            ) {
                return e;
            }
        }
    }
    match finish_hard(
        at,
        num_atoms,
        p,
        &[tg_h, cg_minus, cg_plus],
        aatg,
        counts,
        removed,
        n_canceled,
        n_neutralized,
    ) {
        Some(e) => e,
        None => removed,
    }
}

/// `mark_at_type`: sets `at_type` of all atoms and the totals; returns the
/// number of charges
pub fn mark_at_type(atom: &mut [InpAtom], num_atoms: i32, totals: Option<&mut [i32]>) -> i32 {
    let mut totals = totals;
    if let Some(t) = totals.as_deref_mut() {
        t[..attot::ARRAY_LEN].fill(0);
    }
    for i in 0..num_atoms as usize {
        let (type_, _) = get_atom_charge_type(atom, i, totals.as_deref_mut(), 0);
        atom[i].at_type = type_ as _;
    }
    totals.map_or(0, |t| t[attot::NUM_CHARGES])
}

/// `RemoveNPProtonsAndAcidCharges`: the removal and addition of protons that
/// neutralizes the structure (TECHMAN-5.1, 5.2); returns 0 or an error code
pub fn remove_np_protons_and_acid_charges(
    at: &mut [InpAtom],
    num_atoms: i32,
    aatg: &mut BnAatg,
    p: &mut BnStruct,
    bd: &mut BnData,
) -> i32 {
    let mut n_canceled = 0;
    let taut_flags = aatg
        .t_group_info
        .as_ref()
        .expect("t_group_info")
        .b_taut_flags;
    let allow_hard = taut_flags & TG_FLAG_TEST_TAUT__SALTS != 0
        && taut_flags & TG_FLAG_TEST_TAUT2_SALTS64 != 0
        && taut_flags & TG_FLAG_MOVE_POS_CHARGES != 0
        && taut_flags & TG_FLAG_HARD_ADD_REM_PROTONS != 0;
    if aatg.n_marked_atom.is_empty() {
        aatg.n_marked_atom = vec![0; num_atoms as usize];
        aatg.n_num_found = 0;
    }
    let mut ret = 0;
    let mut b_error = 0;
    'body: {
        let num = totals(aatg)[attot::NUM_NP_PROTON] + totals(aatg)[attot::NUM_OH_PLUS];
        if num != 0 {
            let mut tgi = aatg.t_group_info.take().expect("t_group_info");
            ret = simple_remove_hplus_npo(
                at,
                num_atoms,
                aatg.n_at_type_totals.as_deref_mut(),
                &mut tgi,
            );
            aatg.t_group_info = Some(tgi);
            if ret != num {
                b_error = BNS_PROGRAM_ERR;
                break 'body;
            }
            if ret > 0 {
                aatg.t_group_info
                    .as_mut()
                    .unwrap()
                    .tni
                    .b_normalization_flags |= FLAG_PROTON_NPO_SIMPLE_REMOVED;
            }
        }
        if totals(aatg)[attot::NUM_NP_PLUS] != 0 && allow_hard {
            ret = hard_remove_hplus_np(at, num_atoms, true, Some(&mut n_canceled), aatg, p, bd);
            if is_bns_error(ret) {
                b_error = ret;
                break 'body;
            }
            if ret > 0 {
                aatg.t_group_info
                    .as_mut()
                    .unwrap()
                    .tni
                    .b_normalization_flags |= FLAG_PROTON_NP_HARD_REMOVED;
            }
        }
        let tot_charge = totals(aatg)[attot::TOT_CHARGE];
        if tot_charge > 0 {
            ret = simple_remove_acidic_protons(at, num_atoms, aatg, tot_charge);
            if is_bns_error(ret) {
                b_error = ret;
                break 'body;
            }
            if ret > 0 {
                aatg.t_group_info
                    .as_mut()
                    .unwrap()
                    .tni
                    .b_normalization_flags |= FLAG_PROTON_AC_SIMPLE_REMOVED;
            }
            if totals(aatg)[attot::TOT_CHARGE] > 0 && allow_hard {
                let n = totals(aatg)[attot::TOT_CHARGE];
                ret = hard_remove_acidic_protons(
                    at,
                    num_atoms,
                    aatg,
                    n,
                    Some(&mut n_canceled),
                    p,
                    bd,
                );
                if is_bns_error(ret) {
                    b_error = ret;
                    break 'body;
                }
                if ret > 0 {
                    let ret2 = simple_remove_acidic_protons(at, num_atoms, aatg, ret);
                    if ret2 != ret {
                        b_error = BNS_PROGRAM_ERR;
                        break 'body;
                    }
                    aatg.t_group_info
                        .as_mut()
                        .unwrap()
                        .tni
                        .b_normalization_flags |= FLAG_PROTON_AC_HARD_REMOVED;
                }
            }
        } else if tot_charge < 0 {
            ret = simple_add_acidic_protons(at, num_atoms, aatg, -tot_charge);
            if is_bns_error(ret) {
                b_error = ret;
                break 'body;
            }
            if ret > 0 {
                aatg.t_group_info
                    .as_mut()
                    .unwrap()
                    .tni
                    .b_normalization_flags |= FLAG_PROTON_AC_SIMPLE_ADDED;
            }
            if totals(aatg)[attot::TOT_CHARGE] < 0 && allow_hard {
                let n = -totals(aatg)[attot::TOT_CHARGE];
                ret = hard_add_acidic_protons(at, num_atoms, aatg, n, Some(&mut n_canceled), p, bd);
                if is_bns_error(ret) {
                    b_error = ret;
                    break 'body;
                }
                if ret > 0 {
                    let ret2 = simple_add_acidic_protons(at, num_atoms, aatg, ret);
                    if ret2 != ret {
                        b_error = BNS_PROGRAM_ERR;
                        break 'body;
                    }
                    aatg.t_group_info
                        .as_mut()
                        .unwrap()
                        .tni
                        .b_normalization_flags |= FLAG_PROTON_AC_HARD_ADDED;
                }
            }
        }
        if n_canceled != 0 {
            aatg.t_group_info
                .as_mut()
                .unwrap()
                .tni
                .b_normalization_flags |= FLAG_PROTON_CHARGE_CANCEL;
        }
    }
    if b_error != 0 {
        ret = if is_bns_error(b_error) {
            b_error
        } else {
            BNS_PROGRAM_ERR
        };
    }
    ret
}

// ---- testing the bonds of a structure for alternation ----

use super::ichierr::{BNS_ALTBOND_ERR, BNS_CANT_SET_BOND};

/// `nMaxFlow2Check`
pub fn n_max_flow2_check(p: &BnStruct, iedge: i32) -> i32 {
    (p.edge[iedge as usize].cap & EDGE_FLOW_MASK).min(MAX_BOND_EDGE_CAP)
}

/// `nCurFlow2Check`
pub fn n_cur_flow2_check(p: &BnStruct, iedge: i32) -> i32 {
    p.edge[iedge as usize].flow & EDGE_FLOW_MASK
}

/// `nMinFlow2Check`: the flow of the edge cannot be smaller than this
pub fn n_min_flow2_check(p: &BnStruct, iedge: i32) -> i32 {
    let e = &p.edge[iedge as usize];
    let v1 = e.neighbor1 as usize;
    let v2 = v1 ^ e.neighbor12 as usize;
    let f12 = e.flow & EDGE_FLOW_MASK;
    if f12 <= 0 {
        return 0;
    }
    let rescap = |v: usize| -> i32 {
        let vert = &p.vert[v];
        let mut r = 0;
        for i in 0..vert.num_adj_edges as usize {
            let ie = p.iedge[vert.iedge + i];
            if ie == iedge {
                continue;
            }
            let ed = &p.edge[ie as usize];
            r += (ed.cap & EDGE_FLOW_MASK) - (ed.flow & EDGE_FLOW_MASK);
        }
        r
    };
    let rescap12 = rescap(v1).min(rescap(v2)).min(f12);
    f12 - rescap12
}

/// `bNeedToTestTheFlow`: false if the bond is already known to have the
/// bond order `n_test_flow + 1`
pub fn b_need_to_test_the_flow(
    bond_type: i32,
    n_test_flow: i32,
    b_test_for_non_stereo_bond: bool,
) -> bool {
    use super::ichitaut::{BOND_ALT12NS, BOND_ALTERN};
    let n_bond_type = BOND_TYPE_MASK & bond_type;
    let n_attrib = BOND_MARK_MASK & bond_type;
    if n_attrib != 0 || n_bond_type == BOND_ALTERN || n_bond_type == BOND_ALT12NS {
        let known: &[i32] = match (b_test_for_non_stereo_bond, n_test_flow) {
            (true, 0) => &[BOND_MARK_ALT12NS, BOND_MARK_ALT123, BOND_MARK_ALT13],
            (true, 1) => &[BOND_MARK_ALT12NS, BOND_MARK_ALT123, BOND_MARK_ALT23],
            (false, 0) => &[
                BOND_MARK_ALT12,
                BOND_MARK_ALT12NS,
                BOND_MARK_ALT123,
                BOND_MARK_ALT13,
            ],
            (false, 1) => &[
                BOND_MARK_ALT12,
                BOND_MARK_ALT12NS,
                BOND_MARK_ALT123,
                BOND_MARK_ALT23,
            ],
            (_, 2) => &[BOND_MARK_ALT13, BOND_MARK_ALT123, BOND_MARK_ALT23],
            _ => &[],
        };
        if known.contains(&n_attrib) {
            return false;
        }
    }
    true
}

/// `nBondsValenceInpAt`: (sum of the bond orders, number of alternating
/// bonds, number of wrong bonds)
pub fn n_bonds_valence_inp_at(at: &InpAtom) -> (i32, i32, i32) {
    use super::ichitaut::BOND_ALTERN;
    let (mut valence, mut alt, mut wrong) = (0, 0, 0);
    for j in 0..at.valence as usize {
        let bt = i32::from(at.bond_type[j]) & BOND_TYPE_MASK;
        if bt == 0 || bt == BOND_TYPE_SINGLE || bt == BOND_TYPE_DOUBLE || bt == BOND_TYPE_TRIPLE {
            valence += bt;
        } else if bt == BOND_ALTERN {
            alt += 1;
        } else {
            wrong += 1;
        }
    }
    match alt {
        0 => {}
        1 => {
            valence += 1;
            wrong += 1;
        }
        _ => valence += alt + 1,
    }
    (valence, alt, wrong)
}

/// `BnsAdjustFlowBondsRad`: matches the bonds to the valences (aromatic
/// bonds and radicals are augmented to the lowest multiplicity); returns the
/// total flow change or an error
pub fn bns_adjust_flow_bonds_rad(
    p: &mut BnStruct,
    bd: &mut BnData,
    at: &mut [InpAtom],
    num_atoms: i32,
) -> i32 {
    let mut b_error = 0;
    let mut orig_delta = 0;
    // atoms that have an extra unsatisfied valence due to aromatic bonds
    let mut excess_before: Option<Vec<i32>> = None;
    for i in 0..num_atoms as usize {
        let (valen, alt, mut ignore) = n_bonds_valence_inp_at(&at[i]);
        let val_minus = i32::from(at[i].chem_bonds_valence) - valen;
        ignore += i32::from(alt > 3);
        if ignore == 0 && val_minus > 0 {
            let v = excess_before.get_or_insert_with(|| vec![0; num_atoms as usize]);
            v[i] = val_minus + i32::from(i32::from(at[i].radical) == RADICAL_DOUBLET);
        }
    }
    loop {
        let mut num_removed = 0;
        let mut ret = run_balanced_network_search(p, bd, BNS_EF_CHNG_FLOW);
        if is_bns_error(ret) {
            b_error = ret;
        } else {
            orig_delta += ret;
            num_removed = p.num_altp;
            if ret > 0 {
                ret = set_bonds_from_bn_struct_flow(p, at, num_atoms, BNS_EF_SAVE_ALL);
                if is_bns_error(ret) {
                    b_error = ret;
                }
                ret = restore_bn_struct_flow(p, BNS_EF_SAVE_ALL);
                if is_bns_error(ret) {
                    b_error = ret;
                }
            }
            re_init_bn_struct_alt_paths(p);
        }
        if !(num_removed != 0 && num_removed == p.max_altp && b_error == 0) {
            break;
        }
    }
    if b_error == 0 {
        if let Some(before) = &excess_before {
            for i in 0..num_atoms as usize {
                if before[i] == 0 {
                    continue;
                }
                let (valen, _, ignore) = n_bonds_valence_inp_at(&at[i]);
                let val_minus = i32::from(at[i].chem_bonds_valence) - valen;
                let excess = val_minus + i32::from(i32::from(at[i].radical) == RADICAL_DOUBLET);
                if ignore != 0 || before[i] - excess != 1 {
                    b_error = BNS_ALTBOND_ERR;
                    break;
                }
            }
        }
    }
    if b_error != 0 { b_error } else { orig_delta }
}

/// `BnsTestAndMarkAltBonds`: tests every bond (of the given type, or all
/// if `n_bond_type_to_test` is 0) for the other bond orders it can take and
/// marks the alternating bonds; returns the number of changes or an error
pub fn bns_test_and_mark_alt_bonds(
    p: &mut BnStruct,
    bd: &mut BnData,
    at: &mut [InpAtom],
    num_atoms: i32,
    fcd: &mut [BnsFlowChanges],
    b_change_flow: i32,
    n_bond_type_to_test: i32,
) -> i32 {
    let mut b_error = 0;
    let mut n_changes = 0;
    let b_test_non_stereo = p.tot_st_cap > p.tot_st_flow;
    let mut iat = 0usize;
    while iat < num_atoms as usize && b_error == 0 {
        let mut ineigh = 0usize;
        while ineigh < at[iat].valence as usize && b_error == 0 {
            let neigh = at[iat].neighbor[ineigh] as usize;
            if neigh < iat {
                ineigh += 1;
                continue;
            }
            let iedge = p.iedge[p.vert[iat].iedge + ineigh];
            if p.edge[iedge as usize].forbidden != 0 {
                ineigh += 1;
                continue;
            }
            if n_bond_type_to_test != 0
                && (i32::from(at[iat].bond_type[ineigh]) & BOND_TYPE_MASK) != n_bond_type_to_test
            {
                ineigh += 1;
                continue;
            }
            let mut n_min = n_min_flow2_check(p, iedge);
            let n_max = n_max_flow2_check(p, iedge);
            let n_cur = n_cur_flow2_check(p, iedge);
            if n_min == n_max {
                if n_max != 0 && b_test_non_stereo {
                    n_min = 0.max(n_max - (p.tot_st_cap - p.tot_st_flow));
                } else {
                    ineigh += 1;
                    continue;
                }
            }
            let mut n_test = n_min;
            while n_test <= n_max && b_error == 0 {
                let this_test = n_test;
                n_test += 1;
                if this_test == n_cur {
                    continue;
                }
                if !b_need_to_test_the_flow(
                    i32::from(at[iat].bond_type[ineigh]),
                    this_test,
                    b_test_non_stereo,
                ) {
                    continue;
                }
                let n_dots = set_flow_to_check_one_bond(p, iedge, this_test, fcd);
                if is_bns_error(n_dots) {
                    if n_dots == BNS_CANT_SET_BOND {
                        let ret = restore_flow_after_check_one_bond(p, fcd);
                        if !is_bns_error(ret) {
                            continue;
                        }
                    }
                    b_error = n_dots;
                } else {
                    if n_dots > 0 {
                        let mut ret = run_balanced_network_search(p, bd, b_change_flow);
                        if is_bns_error(ret) {
                            b_error = ret;
                        } else if ret > 0 {
                            if 2 * ret == n_dots {
                                ret = set_bonds_after_check_one_bond(
                                    p,
                                    fcd,
                                    this_test,
                                    at,
                                    num_atoms,
                                    b_change_flow,
                                );
                                if is_bns_error(ret) {
                                    b_error = ret;
                                } else {
                                    n_changes += ret & 1;
                                    ret = set_bonds_from_bn_struct_flow(
                                        p,
                                        at,
                                        num_atoms,
                                        b_change_flow,
                                    );
                                    if is_bns_error(ret) || ret < 0 {
                                        b_error = ret;
                                    } else {
                                        n_changes += ret & 1;
                                    }
                                }
                            }
                            ret = restore_bn_struct_flow(p, b_change_flow & BNS_EF_CHNG_RSTR);
                            if is_bns_error(ret) {
                                b_error = ret;
                            }
                        }
                        re_init_bn_struct_alt_paths(p);
                    } else if n_dots == 0 {
                        let ret = set_bonds_after_check_one_bond(
                            p,
                            fcd,
                            this_test,
                            at,
                            num_atoms,
                            b_change_flow,
                        );
                        if is_bns_error(ret) {
                            b_error = ret;
                        } else {
                            n_changes += ret & 1;
                        }
                    }
                }
                let ret = restore_flow_after_check_one_bond(p, fcd);
                if is_bns_error(ret) {
                    b_error = ret;
                }
            }
            ineigh += 1;
        }
        iat += 1;
    }
    if b_error != 0 { b_error } else { n_changes }
}

// ---- fixing bonds around removed ion pairs ----

use super::util::{
    EL_NUMBER_C, EL_NUMBER_N, EL_NUMBER_O, get_el_valence, ion_el_group, is_el_a_metal,
    n_no_metal_bonds_valence, n_no_metal_neigh_index, n_no_metal_num_bonds,
    n_no_metal_other_neigh_index, n_no_metal_other_neigh_index2, num_of_h,
};

fn forbid_bond(p: &mut BnStruct, atom: usize, bond: usize, mask: i8) {
    let ie = p.iedge[p.vert[atom].iedge + bond] as usize;
    p.edge[ie].forbidden |= mask;
}

/// `fix_special_bonds`: forbids the edges of the bonds that were set by
/// `remove_ion_pairs` around N(V), S(IV), S(VI) and S(VIII); returns the
/// number of fixed bonds
pub fn fix_special_bonds(
    p: &mut BnStruct,
    at: &[InpAtom],
    num_atoms: i32,
    forbidden_mask: i32,
) -> i32 {
    let mask = forbidden_mask as i8;
    let mut num_changes = 0;
    p.edge_forbidden_mask |= mask;
    let bt = |i: usize, k: usize| i32::from(at[i].bond_type[k]) & BOND_TYPE_MASK;
    let nm_bonds = |i: usize| n_no_metal_num_bonds(at, i);
    let nm_val = |i: usize| n_no_metal_bonds_valence(at, i);
    let group = |i: usize| ion_el_group(at[i].el_number);
    for i in 0..num_atoms as usize {
        let a = &at[i];
        if a.charge != 0 || a.radical != 0 {
            continue;
        }
        let ev = get_el_valence(i32::from(a.el_number), 0, 0);
        let excess = i32::from(a.chem_bonds_valence) + a.numh() - ev;
        let nm_excess = nm_val(i) + a.numh() - ev;
        let no_h = num_of_h(at, i) == 0;
        if excess >= 2 && no_h && nm_excess == 2 && group(i) == EL_NUMBER_N {
            // N(V), no H
            let nb = nm_bonds(i);
            if nb == 2 {
                let i1 = n_no_metal_neigh_index(at, i);
                if i1 < 0 {
                    continue;
                }
                let n1 = i32::from(a.neighbor[i1 as usize]);
                let i2 = n_no_metal_other_neigh_index(at, i, n1);
                if i2 < 0 {
                    continue;
                }
                let n2 = i32::from(a.neighbor[i2 as usize]);
                forbid_bond(p, i, i1 as usize, mask);
                forbid_bond(p, i, i2 as usize, mask);
                num_changes += 1;
                // forbid the single bond beyond the neighboring =N- as in #N=N-
                let n3 = if bt(i, i1 as usize) == BOND_TYPE_DOUBLE {
                    n1
                } else {
                    n2
                } as usize;
                if at[n3].numh() == 0
                    && nm_bonds(n3) == 2
                    && nm_val(n3) == 3
                    && group(n3) == EL_NUMBER_N
                {
                    let k = n_no_metal_other_neigh_index(at, n3, i as i32);
                    if k >= 0 {
                        forbid_bond(p, n3, k as usize, mask);
                        num_changes += 1;
                    }
                }
            } else if nb == 3 {
                let j0 = n_no_metal_neigh_index(at, i);
                if j0 < 0 {
                    continue;
                }
                let m0 = i32::from(a.neighbor[j0 as usize]);
                let j1 = n_no_metal_other_neigh_index(at, i, m0);
                if j1 < 0 {
                    continue;
                }
                let m1 = i32::from(a.neighbor[j1 as usize]);
                let j2 = n_no_metal_other_neigh_index2(at, i, m0, m1);
                if j2 < 0 {
                    continue;
                }
                let j = [j0 as usize, j1 as usize, j2 as usize];
                let m = [m0 as usize, m1 as usize, a.neighbor[j2 as usize] as usize];
                let (mut num_o, mut num_n) = (0, 0);
                let (mut k_o, mut k_n) = (0usize, 0usize);
                for k in 0..3 {
                    let n1 = m[k];
                    if group(n1) == EL_NUMBER_N {
                        k_n = k;
                        num_n += 1;
                    } else if group(n1) == EL_NUMBER_O && nm_bonds(n1) == 1 {
                        k_o = k;
                        num_o += 1;
                    }
                }
                let mut num_other = 0;
                if num_o == 1
                    && at[m[k_o]].charge == 0
                    && at[m[k_o]].radical == 0
                    && bt(i, j[k_o]) == BOND_TYPE_DOUBLE
                {
                    forbid_bond(p, i, j[k_o], mask);
                    num_changes += 1;
                    num_other += 1;
                }
                if num_other == 0
                    && num_o <= 1
                    && num_n == 1
                    && at[m[k_n]].charge == 0
                    && at[m[k_n]].radical == 0
                    && bt(i, j[k_n]) == BOND_TYPE_DOUBLE
                {
                    forbid_bond(p, i, j[k_n], mask);
                    num_changes += 1;
                }
            } else if nb == 4 {
                // -N=N- : fix the double bond and the single bond beyond
                let mut num_n = 0;
                let (mut n2, mut i2) = (0usize, 0usize);
                for i1 in 0..a.valence as usize {
                    let n1 = a.neighbor[i1] as usize;
                    if bt(i, i1) == BOND_TYPE_DOUBLE
                        && !is_el_a_metal(i32::from(at[n1].el_number))
                        && group(n1) == EL_NUMBER_N
                    {
                        num_n += 1;
                        n2 = n1;
                        i2 = i1;
                    }
                }
                if num_n == 1 && at[n2].numh() == 0 && nm_bonds(n2) == 2 && nm_val(n2) == 3 {
                    let i3 = n_no_metal_other_neigh_index(at, n2, i as i32);
                    if i3 >= 0 && bt(n2, i3 as usize) == BOND_TYPE_SINGLE {
                        forbid_bond(p, n2, i3 as usize, mask);
                        num_changes += 1;
                        forbid_bond(p, i, i2, mask);
                        num_changes += 1;
                    }
                }
            }
        } else if excess >= 2
            && no_h
            && nm_excess == 2
            && a.el_number != EL_NUMBER_O
            && group(i) == EL_NUMBER_O
            && nm_bonds(i) == 3
        {
            // S(IV), no H, one double bond, total 3 bonds
            let (mut num_n, mut num_oh, mut num_om, mut num_o, mut num_other) = (0, 0, 0, 0, 0);
            let (mut n2, mut i2, mut n3, mut i4) = (0usize, 0usize, 0usize, 0usize);
            for i1 in 0..a.valence as usize {
                let n1 = a.neighbor[i1] as usize;
                if is_el_a_metal(i32::from(at[n1].el_number)) {
                    continue;
                }
                let b = bt(i, i1);
                if b == BOND_TYPE_DOUBLE {
                    num_n += 1;
                    n2 = n1;
                    i2 = i1;
                    if group(n1) == EL_NUMBER_O {
                        num_o += 1;
                    }
                } else if b == BOND_TYPE_SINGLE
                    && nm_bonds(n1) == 1
                    && nm_val(n1) == 1
                    && group(n1) == EL_NUMBER_O
                {
                    if at[n1].charge == 0 {
                        num_oh += 1;
                        n3 = n1;
                    } else {
                        num_om += 1;
                    }
                } else {
                    num_other += 1;
                    i4 = i1;
                }
            }
            let _ = n3;
            if num_n == 1 && num_o == 1 && num_oh + num_om == 1 {
                if num_other == 1 {
                    forbid_bond(p, i, i4, mask);
                    num_changes += 1;
                }
            } else if num_n == 1 && num_oh == 0 && num_om == 0 {
                let mut found = 0;
                let mut do_not_fix = 0;
                let mut n3 = 0i32;
                if group(n2) == EL_NUMBER_N {
                    let el = at[n2].el_number;
                    for i1 in 0..a.valence as usize {
                        let n1 = a.neighbor[i1] as usize;
                        if bt(i, i1) == BOND_TYPE_SINGLE
                            && (at[n1].numh() != 0 || at[n1].charge == -1)
                            && el == at[n1].el_number
                        {
                            n3 = n1 as i32;
                            found += 1;
                        }
                    }
                }
                for i1 in 0..a.valence as usize {
                    if i1 != i2 {
                        let n1 = a.neighbor[i1] as usize;
                        if at[n2].el_number == at[n1].el_number
                            && at[n2].n_ring_system == at[n1].n_ring_system
                        {
                            do_not_fix += 1;
                        }
                    }
                }
                if do_not_fix != 0 {
                } else if found != 0 {
                    if found == 1 {
                        let i4 = n_no_metal_other_neigh_index2(at, i, n2 as i32, n3);
                        if i4 >= 0 {
                            forbid_bond(p, i, i4 as usize, mask);
                            num_changes += 1;
                        }
                    }
                } else {
                    forbid_bond(p, i, i2, mask);
                    num_changes += 1;
                }
            }
        } else if excess >= 4
            && no_h
            && nm_excess == 4
            && a.el_number != EL_NUMBER_O
            && group(i) == EL_NUMBER_O
            && nm_bonds(i) == 4
        {
            num_changes += fix_hypervalent_sulfur(p, at, i, 2, mask);
        } else if excess >= 6
            && no_h
            && nm_excess == 6
            && a.el_number != EL_NUMBER_O
            && group(i) == EL_NUMBER_O
            && nm_bonds(i) == 5
        {
            num_changes += fix_hypervalent_sulfur(p, at, i, 3, mask);
        }
    }
    num_changes
}

/// the S(VI) (`num_double` = 2) and S(VIII) (3) cases of `fix_special_bonds`:
/// the single bonds to non-terminal neighbors are fixed
fn fix_hypervalent_sulfur(
    p: &mut BnStruct,
    at: &[InpAtom],
    i: usize,
    num_double: i32,
    mask: i8,
) -> i32 {
    let a = &at[i];
    let (mut num_n, mut num_oh, mut num_om, mut num_o, mut num_x, mut num_other) =
        (0, 0, 0, 0, 0, 0);
    let mut idx: Vec<usize> = Vec::new();
    for i1 in 0..a.valence as usize {
        let n1 = a.neighbor[i1] as usize;
        if is_el_a_metal(i32::from(at[n1].el_number)) {
            continue;
        }
        let b = i32::from(a.bond_type[i1]) & BOND_TYPE_MASK;
        let n = &at[n1];
        let nh = num_of_h(at, n1);
        let nb = n_no_metal_num_bonds(at, n1);
        let is_o_group = ion_el_group(n.el_number) == EL_NUMBER_O;
        if b == BOND_TYPE_DOUBLE {
            num_n += 1;
            // S_VI_O_PLUS_METAL_FIX_BOND is 1
            if (n.charge == 0 || (n.charge == 1 && n.valence == 2))
                && n.radical == 0
                && nh == 0
                && is_o_group
                && nb == 1
            {
                num_o += 1;
            }
        } else if b == BOND_TYPE_SINGLE
            && nb == 1
            && is_o_group
            && nh <= 1
            && i32::from(n.charge == 0 && nh == 1) + i32::from(n.charge == -1 && nh == 0) == 1
        {
            num_oh += 1;
        } else if b == BOND_TYPE_SINGLE && nb > 1 {
            idx.push(i1);
            num_om += 1;
        } else if b == BOND_TYPE_SINGLE && nb == 1 {
            num_x += 1;
        } else {
            num_other += 1;
        }
    }
    let mut changes = 0;
    if num_n == num_o && num_o == num_double && num_oh + num_om + num_x == 2 && num_other == 0 {
        for &i1 in &idx {
            forbid_bond(p, i, i1, mask);
            changes += 1;
        }
    }
    changes
}

/// `SetForbiddenEdges` (without the disabled explicit-bond fixing): forbids
/// the single bonds of acetyl-like C and nitro N and the bonds of
/// `fix_special_bonds`; returns the number of forbidden bonds
pub fn set_forbidden_edges(
    p: &mut BnStruct,
    at: &[InpAtom],
    num_atoms: i32,
    forbidden_mask: i32,
) -> i32 {
    let mask = forbidden_mask as i8;
    p.edge_forbidden_mask |= mask;
    let mut num_found = 0;
    for i in 0..num_atoms as usize {
        let a = &at[i];
        let acetyl = a.el_number == EL_NUMBER_C && a.valence == 3 && a.chem_bonds_valence == 4;
        let nitro = a.el_number == EL_NUMBER_N
            && a.valence == 3
            && (a.chem_bonds_valence == 4 || a.chem_bonds_valence == 5);
        if !acetyl && !nitro {
            continue;
        }
        let mut num_o = 0;
        let mut bond_to_o = 0;
        let mut pos: i32 = -1;
        let mut val = -1;
        for j in 0..a.valence as usize {
            let nb = &at[a.neighbor[j] as usize];
            let b = i32::from(a.bond_type[j]) & BOND_TYPE_MASK;
            if nb.el_number == EL_NUMBER_O && nb.valence == 1 {
                num_o += 1;
                bond_to_o += b;
            } else {
                pos = j as i32;
                val = b;
            }
        }
        let ok = if acetyl {
            num_o == 2 && bond_to_o == 3 && val == 1
        } else {
            num_o == 2 && (bond_to_o == 3 || bond_to_o == 4) && val == 1
        };
        if ok {
            forbid_bond(p, i, pos as usize, mask);
            num_found += 1;
        }
    }
    // REMOVE_ION_PAIRS_FIX_BONDS is 1; RESET_EDGE_FORBIDDEN_MASK is 1
    num_found + fix_special_bonds(p, at, num_atoms, forbidden_mask)
}

// ---- blocks and non-stereo bonds of alternating-bond graphs ----

/// `BT_*` bond types of `nBondTypeInpAltBns` (stored in `BnsEdge::pass`)
pub const BT_ALTERN_BOND: i8 = 1;
pub const BT_OTHER_ALTERN_BOND: i8 = 2;
pub const BT_ALTERN_NS_BOND: i8 = 4;
pub const BT_TAUTOM_BOND: i8 = 8;
pub const BT_IGNORE_BOND: i8 = 0;
pub const BT_ALT_BOND_MASK: i8 = BT_ALTERN_BOND | BT_OTHER_ALTERN_BOND;

/// `MarkRingSystemsAltBns` (`BNS_MARK_ONLY_BLOCKS` is 1): finds the blocks
/// of the graph of alternating bonds; `edge.flow` gets the block number and
/// `edge.cap` the number of atoms in the block. Returns the number of
/// blocks.
pub fn mark_ring_systems_alt_bns(p: &mut BnStruct) -> i32 {
    let num_atoms = p.num_atoms as usize;
    let num_bonds = p.num_bonds as usize;
    let mut stack_atom = vec![0usize; num_atoms];
    let mut ring_stack = vec![0usize; num_atoms];
    let mut bond_stack = vec![0usize; num_bonds];
    let mut dfs_number = vec![0u16; num_atoms];
    let mut low_number = vec![0u16; num_atoms];
    let mut neigh_numb = vec![0usize; num_atoms];
    let mut n_blocks = 0;
    let valence = |p: &BnStruct, v: usize| p.vert[v].num_adj_edges as usize;
    let iedge = |p: &BnStruct, v: usize, k: usize| p.iedge[p.vert[v].iedge + k] as usize;
    for start in 0..num_atoms {
        if dfs_number[start] != 0 {
            continue;
        }
        let found =
            (0..valence(p, start)).any(|i| p.edge[iedge(p, start, i)].pass & BT_ALTERN_BOND != 0);
        if !found {
            continue;
        }
        let mut n_dfs: u16 = 0;
        let mut top_stack: isize = -1;
        let mut top_ring: isize = -1;
        let mut top_bond: isize = -1;
        neigh_numb.fill(0);
        n_dfs += 1;
        low_number[start] = n_dfs;
        dfs_number[start] = n_dfs;
        top_stack += 1;
        stack_atom[top_stack as usize] = start;
        top_ring += 1;
        ring_stack[top_ring as usize] = start;
        loop {
            // advance
            let mut i;
            loop {
                i = stack_atom[top_stack as usize];
                let j = neigh_numb[i];
                if valence(p, i) <= j {
                    break;
                }
                neigh_numb[i] += 1;
                let w = iedge(p, i, j);
                if p.edge[w].pass & BT_ALT_BOND_MASK == 0 {
                    continue;
                }
                let u = (p.edge[w].neighbor12 as usize) ^ i;
                if dfs_number[u] == 0 && num_bonds != 0 {
                    top_stack += 1;
                    stack_atom[top_stack as usize] = u;
                    top_ring += 1;
                    ring_stack[top_ring as usize] = u;
                    top_bond += 1;
                    bond_stack[top_bond as usize] = w;
                    n_dfs += 1;
                    low_number[u] = n_dfs;
                    dfs_number[u] = n_dfs;
                } else if (top_stack == 0 || u != stack_atom[top_stack as usize - 1])
                    && dfs_number[u] < dfs_number[i]
                    && num_bonds != 0
                {
                    // back edge, first visit
                    top_bond += 1;
                    bond_stack[top_bond as usize] = w;
                    if low_number[i] > dfs_number[u] {
                        low_number[i] = dfs_number[u];
                    }
                }
            }
            neigh_numb[i] = 0;
            // back up
            if i != start {
                let u = stack_atom[top_stack as usize - 1];
                if low_number[i] >= dfs_number[u] {
                    n_blocks += 1;
                    let mut n_at_in_block: i32 = 1;
                    while top_ring >= 0 {
                        let j = ring_stack[top_ring as usize];
                        top_ring -= 1;
                        n_at_in_block += 1;
                        if i == j {
                            break;
                        }
                    }
                    while top_bond >= 0 {
                        let w = bond_stack[top_bond as usize];
                        top_bond -= 1;
                        let e = &mut p.edge[w];
                        e.flow = n_blocks;
                        e.cap = n_at_in_block;
                        let n1 = e.neighbor1 as usize;
                        let n12 = e.neighbor12 as usize;
                        if (i == n1 && u == (i ^ n12)) || (u == n1 && i == (u ^ n12)) {
                            break;
                        }
                    }
                } else if low_number[u] > low_number[i] {
                    low_number[u] = low_number[i];
                }
            }
            top_stack -= 1;
            if top_stack < 0 {
                break;
            }
        }
    }
    n_blocks
}

/// `ReInitBnStructForAltBns`: strips the groups and classifies the bonds by
/// the `BT_*` types; returns the number of alternating bonds to test or an
/// error
pub fn re_init_bn_struct_for_alt_bns(
    p: &mut BnStruct,
    at: &mut [InpAtom],
    num_atoms: i32,
    b_unkn_alt_as_no_stereo: bool,
) -> i32 {
    use super::ichitaut::{BOND_ALT12NS, BOND_ALTERN, BOND_TAUTOM};
    let mut num_to_test = 0;
    if b_unkn_alt_as_no_stereo {
        for j in 0..p.num_edges as usize {
            p.edge[j].pass = 0;
        }
    }
    let ret = re_init_bn_struct(p, at, num_atoms as usize, false);
    if ret != 0
        || p.num_atoms != num_atoms
        || p.num_vertices != num_atoms
        || p.num_bonds != p.num_edges
    {
        return super::ichierr::BNS_REINIT_ERR;
    }
    for v in 0..num_atoms as usize {
        let nadj = p.vert[v].num_adj_edges as usize;
        for j in 0..nadj {
            let ie = p.iedge[p.vert[v].iedge + j] as usize;
            if p.edge[ie].neighbor1 as usize != v {
                continue;
            }
            let mut bond_type = i32::from(at[v].bond_type[j]) & BOND_TYPE_MASK;
            let v2 = (p.edge[ie].neighbor12 as usize) ^ v;
            if at[v].endpoint != 0 || at[v2].endpoint != 0 {
                bond_type = 0;
            }
            // FIX_EITHER_DB_AS_NONSTEREO is 0
            let bt = if bond_type == BOND_ALTERN {
                num_to_test += 1;
                BT_ALTERN_BOND
            } else if bond_type == BOND_ALT_123
                || bond_type == BOND_ALT_13
                || bond_type == BOND_ALT_23
            {
                BT_OTHER_ALTERN_BOND
            } else if bond_type == BOND_TAUTOM {
                BT_TAUTOM_BOND
            } else if bond_type == BOND_ALT12NS {
                BT_ALTERN_NS_BOND
            } else {
                BT_IGNORE_BOND
            };
            let e = &mut p.edge[ie];
            e.pass = bt;
            e.cap = 0;
            e.flow = 0;
            // nBondNonStereoAltBns, nBlockNumberAltBns and nNumAtInBlockAltBns
            // share cap and flow; RESET_EDGE_FORBIDDEN_MASK is 1
            e.forbidden &= p.edge_forbidden_mask;
        }
        let st = &mut p.vert[v].st_edge;
        st.cap0 = 0;
        st.cap = 0;
        st.flow0 = 0;
        st.flow = 0;
    }
    num_to_test
}

/// `MarkNonStereoAltBns`: marks the alternating bonds outside rings as
/// non-stereo; returns the number of marked bonds
pub fn mark_non_stereo_alt_bns(
    p: &BnStruct,
    at: &mut [InpAtom],
    num_atoms: i32,
    b_unkn_alt_as_no_stereo: bool,
) -> i32 {
    use super::ichitaut::{BOND_ALT12NS, BOND_ALTERN};
    if p.num_atoms != num_atoms || p.num_vertices != num_atoms || p.num_bonds != p.num_edges {
        return super::ichierr::BNS_REINIT_ERR;
    }
    let mut ret = 0;
    for ibond in 0..p.num_bonds as usize {
        let e = &p.edge[ibond];
        if e.pass != BT_ALTERN_BOND && e.pass != BT_IGNORE_BOND {
            continue;
        }
        let iat1 = e.neighbor1 as usize;
        let iat2 = (e.neighbor12 as usize) ^ iat1;
        let ib1 = e.neigh_ord[0] as usize;
        let ib2 = e.neigh_ord[1] as usize;
        let is_alt = (i32::from(at[iat1].bond_type[ib1]) & BOND_TYPE_MASK) == BOND_ALTERN;
        if (e.pass == BT_ALTERN_BOND && e.cap <= 3) || (e.pass == BT_IGNORE_BOND && is_alt) {
            if b_unkn_alt_as_no_stereo {
                if is_alt {
                    at[iat1].bond_stereo[ib1] = super::ichister::STEREO_DBLE_EITHER;
                    at[iat2].bond_stereo[ib2] = super::ichister::STEREO_DBLE_EITHER;
                    ret += 1;
                }
            } else {
                at[iat1].bond_type[ib1] = BOND_ALT12NS as u8;
                at[iat2].bond_type[ib2] = BOND_ALT12NS as u8;
                ret += 1;
            }
        }
    }
    ret
}

// ---- adding and removing protons when reconstructing a structure ----

/// `RI_ERR_PROGR`
pub const RI_ERR_PROGR: i32 = -3;

/// `bHasChargedNeighbor`
pub fn b_has_charged_neighbor(at: &[InpAtom], iat: usize) -> bool {
    (0..at[iat].valence as usize).any(|i| at[at[iat].neighbor[i] as usize].charge != 0)
}

/// the test "the atom would get the proton back" shared by the two
/// directions of `AddRemoveProtonsRestr`
fn type_matches_classes(type_: i32, mask: Option<i32>, classes: &[(i32, i32)]) -> bool {
    type_ != 0
        && classes
            .iter()
            .any(|c| type_ & c.0 != 0 && mask.unwrap_or(0) != 0)
}

/// `AddRemoveProtonsRestr`: adds (`*num_protons_to_add` > 0) or removes
/// (< 0) protons so that the normalization of the reconstructed structure
/// gives the original one; returns the number of changed atoms or an error
#[allow(clippy::too_many_arguments)]
pub fn add_remove_protons_restr(
    at: &mut [InpAtom],
    num_atoms: i32,
    num_protons_to_add: &mut i32,
    n_num_prot_added_by_restr: i32,
    b_normalization_flags: u64,
    num_tg: i32,
    _n_charge_revrs: i32,
    _n_charge_in_chi: i32,
) -> i32 {
    use super::util::get_endpoint_valence;
    let na = num_atoms as usize;
    let mut totals = vec![0i32; attot::ARRAY_LEN];
    let mut num_prot = *num_protons_to_add;
    let mut n_success = 0;
    mark_at_type(at, num_atoms, Some(&mut totals));
    let acid_ok = |t: &[i32]| t[attot::TOT_CHARGE] - n_num_prot_added_by_restr;

    if num_prot < 0 && acid_ok(&totals) <= 0 {
        // remove protons from atoms that the normalization would protonate back
        let mut i = 0;
        while i < na && num_prot != 0 {
            let a = &at[i];
            if a.sb_parity[0] != 0
                || a.p_parity != 0
                || a.charge != 0
                || a.num_h == 0
                || a.radical != 0
                || b_has_charged_neighbor(at, i)
            {
                i += 1;
                continue;
            }
            at[i].charge -= 1;
            at[i].num_h -= 1;
            let (type_, mask) = get_atom_charge_type(at, i, None, 0);
            at[i].charge += 1;
            at[i].num_h += 1;
            if type_matches_classes(type_, mask, &AA_TYP_MASK) {
                get_atom_charge_type(at, i, Some(&mut totals), 1);
                at[i].charge -= 1;
                at[i].num_h -= 1;
                get_atom_charge_type(at, i, Some(&mut totals), 0);
                num_prot += 1;
                n_success += 1;
            }
            i += 1;
        }
    }

    if num_prot < 0 && num_tg != 0 && acid_ok(&totals) <= 0 {
        // O=C-NH => (-)O-C=N, O and N are tautomeric endpoints
        for i in 0..na {
            let a = &at[i];
            if a.endpoint == 0
                || a.sb_parity[0] != 0
                || a.p_parity != 0
                || a.radical != 0
                || a.charge != 0
                || b_has_charged_neighbor(at, i)
            {
                continue;
            }
            if a.valence != 1
                || i32::from(a.bond_type[0]) != BOND_TYPE_DOUBLE
                || a.num_h != 0
                || get_endpoint_valence(a.el_number) != 2
            {
                continue;
            }
            let centp = a.neighbor[0] as usize;
            if at[centp].sb_parity[0] != 0
                || at[centp].p_parity != 0
                || !is_centerpoint_elem(at[centp].el_number)
            {
                continue;
            }
            for k in 0..at[centp].valence as usize {
                if i32::from(at[centp].bond_type[k]) != BOND_TYPE_SINGLE {
                    continue;
                }
                let endp2 = at[centp].neighbor[k] as usize;
                let e2 = &at[endp2];
                if e2.endpoint != at[i].endpoint
                    || e2.num_h == 0
                    || e2.charge != 0
                    || e2.sb_parity[0] != 0
                    || e2.p_parity != 0
                    || e2.valence != e2.chem_bonds_valence
                    || 3 != i32::from(e2.chem_bonds_valence) + i32::from(e2.num_h)
                    || get_endpoint_valence(e2.el_number) != 3
                {
                    continue;
                }
                let i0 = (0..at[centp].valence as usize)
                    .find(|&q| at[centp].neighbor[q] as usize == i)
                    .unwrap_or(at[centp].valence as usize);
                let k0 = (0..at[endp2].valence as usize)
                    .find(|&q| at[endp2].neighbor[q] as usize == centp)
                    .unwrap_or(at[endp2].valence as usize);
                if i0 == at[centp].valence as usize || k0 == at[endp2].valence as usize {
                    return RI_ERR_PROGR;
                }
                get_atom_charge_type(at, i, Some(&mut totals), 1);
                get_atom_charge_type(at, endp2, Some(&mut totals), 1);
                at[i].bond_type[0] = at[i].bond_type[0].wrapping_sub(1);
                at[centp].bond_type[i0] = at[centp].bond_type[i0].wrapping_sub(1);
                at[i].chem_bonds_valence -= 1;
                at[i].charge -= 1;
                at[endp2].bond_type[k0] = at[endp2].bond_type[k0].wrapping_add(1);
                at[centp].bond_type[k] = at[centp].bond_type[k].wrapping_add(1);
                at[endp2].chem_bonds_valence += 1;
                at[endp2].num_h -= 1;
                num_prot += 1;
                n_success += 1;
                get_atom_charge_type(at, i, Some(&mut totals), 0);
                get_atom_charge_type(at, endp2, Some(&mut totals), 0);
            }
        }
    }

    if num_prot > 0 {
        // 1. the test of SimpleRemoveAcidicProtons: H(+) to =C-O(-), etc.
        let mut i = 0;
        while i < na && num_prot != 0 && acid_ok(&totals) >= 0 {
            let a = &at[i];
            if a.sb_parity[0] != 0
                || a.p_parity != 0
                || a.num_h != 0
                || a.charge != -1
                || a.radical != 0
                || b_has_charged_neighbor(at, i)
            {
                i += 1;
                continue;
            }
            at[i].charge += 1;
            at[i].num_h += 1;
            let (type_, mask) = get_atom_charge_type(at, i, None, 0);
            at[i].charge -= 1;
            at[i].num_h -= 1;
            if type_matches_classes(type_, mask, &AR_TYP_MASK) {
                get_atom_charge_type(at, i, Some(&mut totals), 1);
                at[i].charge += 1;
                at[i].num_h += 1;
                get_atom_charge_type(at, i, Some(&mut totals), 0);
                num_prot -= 1;
                n_success += 1;
            }
            i += 1;
        }
        // 2. the test of SimpleRemoveHplusNPO
        let mut i = 0;
        while i < na && num_prot != 0 {
            let a = &at[i];
            if a.sb_parity[0] != 0
                || a.p_parity != 0
                || a.charge != 0
                || a.radical != 0
                || b_has_charged_neighbor(at, i)
            {
                i += 1;
                continue;
            }
            at[i].num_h += 1;
            at[i].charge += 1;
            let (type_, mask) = get_atom_charge_type(at, i, None, 0);
            let success = PR_SIMPLE_TYP & type_ != 0 && PR_SIMPLE_MSK & mask.unwrap_or(0) != 0;
            at[i].num_h -= 1;
            at[i].charge -= 1;
            if success {
                get_atom_charge_type(at, i, Some(&mut totals), 1);
                at[i].num_h += 1;
                at[i].charge += 1;
                get_atom_charge_type(at, i, Some(&mut totals), 0);
                num_prot -= 1;
                n_success += 1;
            }
            i += 1;
        }
    }

    if num_prot < 0
        && b_normalization_flags & FLAG_PROTON_AC_HARD_ADDED != 0
        && num_tg == 1
        && acid_ok(&totals) <= 0
    {
        // remove protons from tautomeric N (a specific ADP must be present)
        let (mut n_acceptors, mut n_donors) = (0, 0);
        let aa_co_typ = (ATT_ACIDIC_CO | ATT_NO) as i32;
        let aa_co_msk = atbit::CO | atbit::NO;
        let nh_donor = |at: &[InpAtom], i: usize, t: i32, m: Option<i32>| {
            t == i32::from(ATT_ATOM_N)
                && m == Some(atbit::NP_H)
                && at[i].charge == 0
                && at[i].valence == at[i].chem_bonds_valence
        };
        for i in 0..na {
            let a = &at[i];
            if a.endpoint == 0
                || a.radical != 0
                || a.sb_parity[0] != 0
                || a.p_parity != 0
                || b_has_charged_neighbor(at, i)
            {
                continue;
            }
            let (t, m) = get_atom_charge_type(at, i, None, 0);
            if t & aa_co_typ != 0 && m.unwrap_or(0) & aa_co_msk != 0 {
                n_acceptors += 1;
            } else if nh_donor(at, i, t, m) {
                n_donors += 1;
            }
        }
        let num_max = n_acceptors.min(n_donors);
        let mut num_success = 0;
        let mut i = 0;
        while i < na && num_success < num_max && num_prot < 0 {
            let a = &at[i];
            if a.endpoint == 0
                || a.radical != 0
                || a.sb_parity[0] != 0
                || a.p_parity != 0
                || b_has_charged_neighbor(at, i)
            {
                i += 1;
                continue;
            }
            let (t, m) = get_atom_charge_type(at, i, None, 0);
            if nh_donor(at, i, t, m) {
                get_atom_charge_type(at, i, Some(&mut totals), 1);
                at[i].num_h -= 1;
                at[i].charge -= 1;
                get_atom_charge_type(at, i, Some(&mut totals), 0);
                num_prot += 1;
                num_success += 1;
                n_success += 1;
            }
            i += 1;
        }
    }
    *num_protons_to_add = num_prot;
    n_success
}

/// `AddRemoveIsoProtonsRestr`: distributes isotopic protons over the atoms;
/// returns the number of placed isotopes or an error
pub fn add_remove_iso_protons_restr(
    at: &mut [InpAtom],
    num_atoms: i32,
    num_protons_to_add: &mut [i16; 3],
    num_tg: i32,
) -> i32 {
    use super::strutil::hetero_atom_may_have_xchg_iso_h;
    let na = num_atoms as usize;
    let mut n_success = 0;
    let mode2 = if num_tg != 0 { 1 } else { 0 };
    for mode in 0..=mode2 {
        for iso in (0..=2usize).rev() {
            if num_protons_to_add[iso] == 0 {
                continue;
            }
            if num_protons_to_add[iso] < 0 {
                return RI_ERR_PROGR;
            }
            let mut i = 0;
            while i < na && num_protons_to_add[iso] > 0 {
                let j;
                if mode != 0 {
                    if at[i].endpoint != 0 {
                        j = i;
                    } else {
                        i += 1;
                        continue;
                    }
                } else if at[i].endpoint == 0 && hetero_atom_may_have_xchg_iso_h(at, i) == 1 {
                    j = i;
                } else if at[i].el_number == super::util::EL_NUMBER_H
                    && at[i].charge == 1
                    && at[i].valence == 0
                    && at[i].radical == 0
                    && at[i].iso_atw_diff == 0
                {
                    at[i].iso_atw_diff = 1 + iso as i8;
                    num_protons_to_add[iso] -= 1;
                    n_success += 1;
                    i += 1;
                    continue;
                } else {
                    i += 1;
                    continue;
                }
                let mut num_h = i32::from(at[j].num_h);
                while num_h > 0 && num_protons_to_add[iso] > 0 {
                    at[j].num_iso_h[iso] += 1;
                    at[j].num_h -= 1;
                    num_protons_to_add[iso] -= 1;
                    num_h -= 1;
                    n_success += 1;
                }
                let (mut num_expl_h, mut num_expl_iso_h) = (0usize, 0usize);
                let mut k = 0;
                while k < at[j].valence as usize && na <= at[j].neighbor[k] as usize {
                    let n = at[j].neighbor[k] as usize;
                    num_expl_h += usize::from(at[n].iso_atw_diff == 0);
                    num_expl_iso_h += usize::from(at[n].iso_atw_diff != 0);
                    k += 1;
                }
                while num_expl_h > 0 && num_protons_to_add[iso] > 0 {
                    // the index is that of the C code
                    let n = at[j].neighbor[num_expl_h] as usize;
                    if at[n].iso_atw_diff != 0 {
                        return RI_ERR_PROGR;
                    }
                    at[n].iso_atw_diff = 1 + iso as i8;
                    num_expl_h -= 1;
                    num_expl_iso_h += 1;
                    num_protons_to_add[iso] -= 1;
                    n_success += 1;
                }
                let _ = num_expl_iso_h;
                i += 1;
            }
        }
    }
    n_success
}
