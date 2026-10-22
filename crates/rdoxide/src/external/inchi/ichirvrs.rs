//! The data types of the InChI -> structure restoration
//! (`INCHI_BASE/src/ichirvrs.h`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi_bns::{
    BNS_VERT_TYPE_AUX, BNS_VERT_TYPE_C_GROUP, BNS_VERT_TYPE_C_NEGATIVE, BNS_VERT_TYPE_CARBON_GR,
    BNS_VERT_TYPE_METAL_GR, BNS_VERT_TYPE_SUPER_TGROUP, BNS_VERT_TYPE_TEMP,
};

pub const NOT_READ_INT: i32 = 9998;
pub const VALUE_OCTET: i32 = 8;
pub const TREAT_ATOM_AS_METAL: i32 = 99;

pub const BNS_VT_C_POS: u16 = BNS_VERT_TYPE_C_GROUP;
pub const BNS_VT_C_NEG: u16 = BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE;
pub const BNS_VT_C_POS_C: u16 = BNS_VT_C_POS | BNS_VERT_TYPE_CARBON_GR;
pub const BNS_VT_C_NEG_C: u16 = BNS_VT_C_NEG | BNS_VERT_TYPE_CARBON_GR;
pub const BNS_VT_C_POS_M: u16 = BNS_VT_C_POS | BNS_VERT_TYPE_METAL_GR;
pub const BNS_VT_C_NEG_M: u16 = BNS_VT_C_NEG | BNS_VERT_TYPE_METAL_GR;
pub const BNS_VT_M_GROUP: u16 = BNS_VERT_TYPE_METAL_GR;
pub const BNS_VT_C_POS_ALL: u16 = BNS_VERT_TYPE_SUPER_TGROUP | BNS_VERT_TYPE_C_GROUP;
pub const BNS_VT_C_NEG_ALL: u16 = BNS_VT_C_POS_ALL | BNS_VERT_TYPE_C_NEGATIVE;
pub const BNS_VT_CHRG_STRUCT: u16 = BNS_VERT_TYPE_AUX | BNS_VERT_TYPE_TEMP;
pub const BNS_VT_YVCONNECTOR: u16 = BNS_VERT_TYPE_AUX;

/// `SRM`: the rules how to handle bonds to metal atoms
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Srm {
    pub b_metal_add_flower: i32,
    pub n_metal_min_bond_order: i32,
    pub n_metal_init_edge_flow: i32,
    pub n_metal_init_bond_order: i32,
    pub n_metal2_endpoint_min_bond_order: i32,
    pub n_metal2_endpoint_init_bond_order: i32,
    pub n_metal2_endpoint_init_edge_flow: i32,
    pub n_metal_flower_param_d: i32,
    pub n_metal_max_charge_d: i32,
    pub b_stereo_removes_metal_flag: i32,
    pub b_fix_stereo_bonds: i32,
}

/// `SetUpSrm`: the structure restore parameters
pub fn set_up_srm() -> Srm {
    // FIX_STEREO_BOND_ORDER = 0, METAL_FREE_CHARGE_VAL = 1
    let b_fix_stereo_bonds = 0;
    let mut s = Srm {
        b_fix_stereo_bonds,
        n_metal2_endpoint_min_bond_order: 1,
        n_metal2_endpoint_init_edge_flow: 0,
        b_metal_add_flower: 1,
        n_metal_min_bond_order: 0,
        n_metal_init_edge_flow: 1,
        n_metal_init_bond_order: 1,
        b_stereo_removes_metal_flag: b_fix_stereo_bonds,
        n_metal_flower_param_d: 16,
        n_metal_max_charge_d: 16,
        ..Default::default()
    };
    s.n_metal2_endpoint_init_bond_order =
        s.n_metal2_endpoint_min_bond_order + s.n_metal2_endpoint_init_edge_flow;
    s
}

/// `VAL_AT`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ValAt {
    pub c_do_not_add_h: i8,
    pub c_metal: i8,
    pub c_num_bonds_to_metal: i8,
    pub c_init_flow_to_metal: i8,
    pub c_init_valence_to_metal: i8,
    pub c_init_orig_valence_to_metal: i8,
    pub c_max_flow_to_metal: i8,
    pub c_init_free_valences: i8,
    pub c_init_charge: i8,
    pub c_num_valence_electrons: i8,
    pub c_periodic_row_number: i8,
    pub c_min_ring_size: i8,
    pub c_periodic_number: u8,
    pub c_cn_list_index: i8,
    pub n_c_minus_group_edge: i32,
    pub n_c_plus_group_edge: i32,
    pub n_metal_group_edge: i32,
    pub n_taut_group_edge: i32,
}

/// `CHARGE_VAL`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ChargeVal {
    pub n_valence: i32,
    pub n_charge: i32,
    pub n_valence_ordering_number: i32,
}

/// `VCF`: the vertex of a charge structure
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Vcf {
    pub type_: u16,
    pub cap: i8,
    pub flow: i8,
    pub valence: i8,
}

/// `ECF`: an edge of a charge structure
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Ecf {
    pub neigh: i16,
    pub cap: i8,
    pub b_forbidden_edge: i8,
    pub flow: i8,
}

/// `C_NODE`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CNode {
    pub v: Vcf,
    pub e: [Ecf; 3],
}

/// `CN_LIST`
#[derive(Debug, Clone, Copy)]
pub struct CnList {
    pub p_cn: &'static [CNode],
    pub bits: i32,
    pub n_initial_charge: i32,
    pub len: i32,
}

pub const CN_BITS_N: i32 = 1;
pub const CN_BITS_P: i32 = 2;
pub const CN_BITS_M: i32 = 4;
pub const CN_BITS_SHIFT: i32 = 3;
pub const CN_BITS_ME: i32 = -1;
pub const CN_LIST_INDEX_ME: i32 = 17;

/// `MAKE_CN_BITS`
pub const fn make_cn_bits(a: i32, b: i32, c: i32, d: i32) -> i32 {
    ((((d << CN_BITS_SHIFT) | c) << CN_BITS_SHIFT) | b) << CN_BITS_SHIFT | a
}

pub const INI_NUM_TCGROUPS: usize = 16;
pub const INC_NUM_TCGROUPS: usize = 16;

/// `TGRF_MINUS_FIRST`
pub const TGRF_MINUS_FIRST: i16 = 1;

/// `TC_GROUP`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TcGroup {
    pub type_: i32,
    pub ord_num: i32,
    pub num_edges: i32,
    pub st_cap: i32,
    pub st_flow: i32,
    pub edges_cap: i32,
    pub edges_flow: i32,
    pub n_vertex_number: i32,
    pub n_forward_edge: i32,
    pub n_backward_edge: i32,
    pub tg_num_h: i16,
    pub tg_num_minus: i16,
    pub tg_set_minus: i16,
    pub tg_restore_flags: i16,
}

/// `TCGR_TYPE`
pub const TCG_NONE: i32 = -1;
pub const TCG_PLUS0: i32 = 0;
pub const TCG_PLUS1: i32 = 1;
pub const TCG_MINUS0: i32 = 2;
pub const TCG_MINUS1: i32 = 3;
pub const TCG_PLUS_C0: i32 = 4;
pub const TCG_PLUS_C1: i32 = 5;
pub const TCG_MINUS_C0: i32 = 6;
pub const TCG_MINUS_C1: i32 = 7;
pub const TCG_PLUS_M0: i32 = 8;
pub const TCG_PLUS_M1: i32 = 9;
pub const TCG_MINUS_M0: i32 = 10;
pub const TCG_MINUS_M1: i32 = 11;
pub const TCG_ME_FLOWER0: i32 = 12;
pub const TCG_ME_FLOWER1: i32 = 13;
pub const TCG_ME_FLOWER2: i32 = 14;
pub const TCG_ME_FLOWER3: i32 = 15;
pub const TCG_PLUS: i32 = 16;
pub const TCG_MINUS: i32 = 17;
pub const NUM_TCGROUP_TYPES: usize = 18;

/// `ALL_TC_GROUPS`
#[derive(Debug, Clone, Default)]
pub struct AllTcGroups {
    pub p_tcg: Vec<TcGroup>,
    pub num_tc_groups: i32,
    pub max_tc_groups: i32,
    pub n_group: [i32; NUM_TCGROUP_TYPES],
    pub n_vertices: i32,
    pub n_edges: i32,
    pub n_add_iedges: i32,
    pub num_atoms: i32,
    pub num_bonds: i32,
    pub num_tgroups: i32,
    pub num_tgroup_edges: i32,
    pub tgroup_charge: i32,
    pub charge_on_atoms: i32,
    pub added_charge: i32,
    pub total_charge: i32,
    pub total_electrons: i32,
    pub total_electrons_metals: i32,
    pub num_metal_atoms: i32,
    pub num_metal_bonds: i32,
    pub n_edge4charge: i32,
    pub n_edge_plus: i32,
    pub n_edge_minus: i32,
    pub i_component: i32,
    pub i_at_no_offset: i32,
}

/// `CC_CAND`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CcCand {
    pub iat: i16,
    pub num_bonds: i8,
    pub chem_valence: i8,
    pub c_metal: i8,
    pub c_num_bonds_to_metal: i8,
    pub c_num_valence_electrons: i8,
    pub c_periodic_row_number: i8,
    pub c_num_charge_states: i8,
    pub el_number: u8,
}
