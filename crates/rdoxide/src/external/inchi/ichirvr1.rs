//! The InChI -> structure restoration, part 1 (`INCHI_BASE/src/ichirvr1.c`):
//! the charge-structure tables, the t-groups of a read InChI, the charge and
//! valence states of an atom.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::InChI;
use super::ichi_bns::{
    BN_MAX_ALTP, BNS_ADD_EDGES, BNS_ADD_SUPER_TGROUP, BNS_EDGE_FORBIDDEN_MASK, BNS_EF_CHNG_FLOW,
    BNS_VERT_TYPE_C_NEGATIVE, BNS_VERT_TYPE_ENDPOINT, BNS_VERT_TYPE_SUPER_TGROUP, BnData, BnStruct,
    BnsAltPath, BnsEdge, BnsVertex, EDGE_FLOW_MASK, EDGE_FLOW_ST_MASK, I_ALTP_END_ATOM,
    I_ALTP_FLOW, I_ALTP_HDR_LEN, I_ALTP_MAX_LEN, I_ALTP_NEIGHBOR, I_ALTP_PATH_LEN,
    I_ALTP_START_ATOM, MAX_BOND_EDGE_CAP, NO_VERTEX, NUM_KINDS_OF_GROUPS, balanced_network_search,
    re_init_bn_data, re_init_bn_struct_alt_paths, run_balanced_network_search,
};
use super::ichi_bns::{
    BNS_VERT_TYPE_ATOM, BNS_VERT_TYPE_AUX, BNS_VERT_TYPE_C_GROUP, BNS_VERT_TYPE_METAL_GR,
    BNS_VERT_TYPE_TEMP, BNS_VERT_TYPE_TGROUP,
};
use super::ichicano::{T_GROUP_HDR_LEN, TGSO_SYMM_IORDER};
use super::ichierr::{
    BNS_BOND_ERR, BNS_CPOINT_ERR, BNS_PROGRAM_ERR, BNS_VERT_EDGE_OVFL, is_bns_error,
};
use super::ichiread::{NO_VALUE_INT, RI_ERR_ALLOC, RI_ERR_PROGR, RI_ERR_SYNTAX};
use super::ichirvrs::*;
use super::ichisort::insertions_sort;
use super::ichitaut::{TGroup, TGroupInfo};
use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, NUM_H_ISOTOPES};
use super::inpdef::{BOND_TYPE_MASK, InpAtom, MAX_NUM_STEREO_BONDS, MAXVAL};
use super::util::{
    EL_NUMBER_N, EL_NUMBER_O, EL_NUMBER_P, EL_NUMBER_S, EL_NUMBER_SE, EL_NUMBER_TE, get_el_valence,
    get_endpoint_valence, if_skip_add_h, is_el_a_metal, n_no_metal_bonds_valence,
    n_no_metal_num_bonds,
};

/// `TGSO_TOTAL_LEN`
pub const TGSO_TOTAL_LEN: usize = 4;

// ChargeStruct fictitious structures `cnList[*]`
pub static CNME: [CNode; 5] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 0,
            flow: 0,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 16,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 16,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 16,
            flow: 16,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 16,
                b_forbidden_edge: 0,
                flow: 16,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS_M,
            cap: 16,
            flow: 16,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 2,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG_M,
            cap: 16,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_M_GROUP,
            cap: 0,
            flow: 0,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 1,
                cap: 3,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNPNPN: [CNode; 5] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 3,
            flow: 3,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 3,
                b_forbidden_edge: 0,
                flow: 3,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 3,
            flow: 3,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 2,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 2,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNNPNP: [CNode; 5] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 3,
            flow: 3,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 3,
                b_forbidden_edge: 0,
                flow: 3,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 3,
            flow: 3,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNNPN: [CNode; 5] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 2,
            flow: 2,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 2,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNPNP: [CNode; 5] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 2,
            flow: 2,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 2,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 2,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNMNP: [CNode; 3] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 2,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNPNM: [CNode; 4] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 2,
            flow: 2,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 2,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 2,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNEN: [CNode; 6] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 6,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS_C,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG_C,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNNMN: [CNode; 5] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 2,
            flow: 2,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 2,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNNE: [CNode; 5] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNNEN: [CNode; 6] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 2,
            flow: 2,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 2,
                b_forbidden_edge: 0,
                flow: 2,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 2,
            flow: 2,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 5,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 4,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 3,
        },
        e: [
            Ecf {
                neigh: 6,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNNP: [CNode; 2] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNPN: [CNode; 3] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 2,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNNM: [CNode; 3] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 1,
            flow: 1,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 1,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_CHRG_STRUCT,
            cap: 1,
            flow: 1,
            valence: 2,
        },
        e: [
            Ecf {
                neigh: 3,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNMN: [CNode; 2] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_NEG,
            cap: 1,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNP_: [CNode; 2] = [
    CNode {
        v: Vcf {
            type_: BNS_VERT_TYPE_ATOM,
            cap: 0,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 2,
                cap: 1,
                b_forbidden_edge: 1,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
    CNode {
        v: Vcf {
            type_: BNS_VT_C_POS,
            cap: 2,
            flow: 0,
            valence: 1,
        },
        e: [
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
            Ecf {
                neigh: 0,
                cap: 0,
                b_forbidden_edge: 0,
                flow: 0,
            },
        ],
    },
];
pub static CNM_: [CNode; 1] = [CNode {
    v: Vcf {
        type_: BNS_VERT_TYPE_ATOM,
        cap: 0,
        flow: 0,
        valence: 0,
    },
    e: [
        Ecf {
            neigh: 0,
            cap: 0,
            b_forbidden_edge: 0,
            flow: 0,
        },
        Ecf {
            neigh: 0,
            cap: 0,
            b_forbidden_edge: 0,
            flow: 0,
        },
        Ecf {
            neigh: 0,
            cap: 0,
            b_forbidden_edge: 0,
            flow: 0,
        },
    ],
}];
pub static CNN_: [CNode; 1] = [CNode {
    v: Vcf {
        type_: BNS_VERT_TYPE_ATOM,
        cap: 0,
        flow: 0,
        valence: 0,
    },
    e: [
        Ecf {
            neigh: 0,
            cap: 0,
            b_forbidden_edge: 0,
            flow: 0,
        },
        Ecf {
            neigh: 0,
            cap: 0,
            b_forbidden_edge: 0,
            flow: 0,
        },
        Ecf {
            neigh: 0,
            cap: 0,
            b_forbidden_edge: 0,
            flow: 0,
        },
    ],
}];

pub const CN_BITS_PNPN: i32 = make_cn_bits(CN_BITS_P, CN_BITS_N, CN_BITS_P, CN_BITS_N);
pub const CN_BITS_NPNP: i32 = make_cn_bits(CN_BITS_N, CN_BITS_P, CN_BITS_N, CN_BITS_P);
pub const CN_BITS_NPN: i32 = make_cn_bits(CN_BITS_N, CN_BITS_P, CN_BITS_N, 0);
pub const CN_BITS_PNP: i32 = make_cn_bits(CN_BITS_P, CN_BITS_N, CN_BITS_P, 0);
pub const CN_BITS_MNP: i32 = make_cn_bits(CN_BITS_M, CN_BITS_N, CN_BITS_P, 0);
pub const CN_BITS_PNM: i32 = make_cn_bits(CN_BITS_P, CN_BITS_N, CN_BITS_M, 0);
pub const CN_BITS_EN: i32 = make_cn_bits(CN_BITS_P | CN_BITS_M, CN_BITS_N, 0, 0);
pub const CN_BITS_NMN: i32 = make_cn_bits(CN_BITS_N, CN_BITS_M, CN_BITS_N, 0);
pub const CN_BITS_NE: i32 = make_cn_bits(CN_BITS_N, CN_BITS_P | CN_BITS_M, 0, 0);
pub const CN_BITS_NEN: i32 = make_cn_bits(CN_BITS_N, CN_BITS_M | CN_BITS_N, CN_BITS_N, 0);
pub const CN_BITS_NP: i32 = make_cn_bits(CN_BITS_N, CN_BITS_P, 0, 0);
pub const CN_BITS_PN: i32 = make_cn_bits(CN_BITS_P, CN_BITS_N, 0, 0);
pub const CN_BITS_NM: i32 = make_cn_bits(CN_BITS_N, CN_BITS_M, 0, 0);
pub const CN_BITS_MN: i32 = make_cn_bits(CN_BITS_M, CN_BITS_N, 0, 0);
pub const CN_BITS_P_: i32 = make_cn_bits(CN_BITS_P, 0, 0, 0);
pub const CN_BITS_M_: i32 = make_cn_bits(CN_BITS_M, 0, 0, 0);
pub const CN_BITS_N_: i32 = make_cn_bits(CN_BITS_N, 0, 0, 0);

macro_rules! cnl {
    ($t:expr, $b:expr, $c:expr) => {
        CnList {
            p_cn: &$t,
            bits: $b,
            n_initial_charge: $c,
            len: $t.len() as i32,
        }
    };
}

/// `cnList[]`
pub static CN_LIST: [CnList; 18] = [
    cnl!(CNPNPN, CN_BITS_PNPN, 0), /*  0 */
    cnl!(CNNPNP, CN_BITS_NPNP, 0), /*  1 */
    cnl!(CNNPN, CN_BITS_NPN, 0),   /*  2 */
    cnl!(CNPNP, CN_BITS_PNP, 0),   /*  3 */
    cnl!(CNMNP, CN_BITS_MNP, 0),   /*  4 */
    cnl!(CNPNM, CN_BITS_PNM, 0),   /*  5 */
    cnl!(CNEN, CN_BITS_EN, 0),     /*  6 */
    cnl!(CNNMN, CN_BITS_NMN, 0),   /*  7 */
    cnl!(CNNE, CN_BITS_NE, 0),     /*  8 */
    cnl!(CNNEN, CN_BITS_NEN, 0),   /*  9 */
    cnl!(CNNP, CN_BITS_NP, 0),     /* 10 */
    cnl!(CNPN, CN_BITS_PN, 0),     /* 11 */
    cnl!(CNNM, CN_BITS_NM, 0),     /* 12 */
    cnl!(CNMN, CN_BITS_MN, 0),     /* 13 */
    cnl!(CNP_, CN_BITS_P_, 0),     /* 14 */
    cnl!(CNM_, CN_BITS_M_, -1),    /* 15 */
    cnl!(CNN_, CN_BITS_N_, 0),     /* 16 */
    cnl!(CNME, CN_BITS_ME, 0),     /* 17 */
];

/// `cnListNumEl`
pub const CN_LIST_NUM_EL: usize = 18;

/// `clear_t_group_info`: zeroes the t-group info keeping the allocations
pub fn clear_t_group_info(ti: &mut TGroupInfo) {
    let mut t_group = ti.t_group.take();
    let mut max_num_t_groups = ti.max_num_t_groups;
    let mut t_group_number = std::mem::take(&mut ti.t_group_number);
    let mut num_t_groups = ti.num_t_groups;
    let mut endpoints = std::mem::take(&mut ti.n_endpoint_atom_number);
    let mut num_endpoints = ti.n_num_endpoints;
    let mut iso_endpoints = std::mem::take(&mut ti.n_isotopic_endpoint_atom_number);
    let mut num_iso_endpoints = ti.n_num_isotopic_endpoints;
    *ti = TGroupInfo::default();
    match t_group.as_mut() {
        Some(v) => {
            for x in v.iter_mut().take(max_num_t_groups) {
                *x = TGroup::default();
            }
        }
        None => max_num_t_groups = 0,
    }
    if t_group_number.is_empty() {
        num_t_groups = 0;
    } else {
        for x in t_group_number.iter_mut().take(num_t_groups) {
            *x = 0;
        }
    }
    if endpoints.is_empty() {
        num_endpoints = 0;
    } else {
        for x in endpoints.iter_mut().take(num_endpoints.max(0) as usize) {
            *x = 0;
        }
    }
    if iso_endpoints.is_empty() {
        num_iso_endpoints = 0;
    } else {
        for x in iso_endpoints
            .iter_mut()
            .take(num_iso_endpoints.max(0) as usize)
        {
            *x = 0;
        }
    }
    ti.t_group = t_group;
    ti.max_num_t_groups = max_num_t_groups;
    ti.t_group_number = t_group_number;
    ti.num_t_groups = num_t_groups;
    ti.n_endpoint_atom_number = endpoints;
    ti.n_num_endpoints = num_endpoints;
    ti.n_isotopic_endpoint_atom_number = iso_endpoints;
    ti.n_num_isotopic_endpoints = num_iso_endpoints;
}

/// `GetTgroupInfoFromInChI`: fills the t-group info of the structure from
/// the tautomeric groups of the InChI and marks the endpoints
pub fn get_tgroup_info_from_inchi(
    ti: &mut TGroupInfo,
    mut at: Option<&mut [InpAtom]>,
    mut endpoint: Option<&mut [u16]>,
    inchi: Option<&InChI>,
) -> i32 {
    clear_t_group_info(ti);
    let Some(p) = inchi else {
        return 0;
    };
    if !(p.len_tautomer > 1 && !p.n_tautomer.is_empty() && p.n_tautomer[0] > 0) {
        return 0;
    }
    let nt = |j: usize| p.n_tautomer.get(j).copied().unwrap_or(0);
    let num_atoms = p.n_number_of_atoms;
    let num_t_groups = i32::from(nt(0));
    let mut len_tg = p.len_tautomer - T_GROUP_HDR_LEN as i32 * num_t_groups - 1;
    if ti.max_num_t_groups as i32 != num_atoms / 2 + 1 || ti.t_group.is_none() {
        ti.max_num_t_groups = (num_atoms / 2 + 1).max(0) as usize;
        ti.t_group = Some(vec![TGroup::default(); ti.max_num_t_groups + 1]);
    }
    if ti.num_t_groups as i32 != num_t_groups || ti.t_group_number.is_empty() {
        ti.num_t_groups = num_t_groups as usize;
        ti.t_group_number = vec![0; (ti.num_t_groups + 1) * TGSO_TOTAL_LEN];
    }
    if len_tg != ti.n_num_endpoints || ti.n_endpoint_atom_number.is_empty() {
        if len_tg < -1 {
            return RI_ERR_ALLOC;
        }
        ti.n_num_endpoints = len_tg;
        ti.n_endpoint_atom_number = vec![0; (len_tg + 1).max(0) as usize];
    }
    let n_t = ti.num_t_groups;
    let mut j = 1usize;
    let mut i = 0usize;
    let mut itg = 0usize;
    while itg < nt(0) as usize && itg <= ti.max_num_t_groups {
        len_tg = i32::from(nt(j));
        {
            let tg = &mut ti.t_group.as_mut().unwrap()[itg];
            tg.num[0] = nt(j + 1).wrapping_add(nt(j + 2));
            tg.num[1] = nt(j + 2);
            tg.n_group_number = (itg + 1) as u16;
        }
        ti.t_group_number[itg] = itg as u16;
        ti.t_group_number[TGSO_SYMM_IORDER * n_t + itg] = itg as u16;
        j += T_GROUP_HDR_LEN;
        len_tg -= T_GROUP_HDR_LEN as i32 - 1;
        {
            let tg = &mut ti.t_group.as_mut().unwrap()[itg];
            tg.n_num_endpoints = len_tg as u16;
            tg.n_first_endpoint_at_no_pos = i as u16;
        }
        while len_tg > 0 {
            let v = (i32::from(nt(j)) - 1) as u16;
            if i < ti.n_endpoint_atom_number.len() {
                ti.n_endpoint_atom_number[i] = v;
            }
            let k = i32::from(v);
            if k < 0 || k > num_atoms {
                return RI_ERR_PROGR;
            }
            if let Some(a) = at.as_deref_mut() {
                if let Some(x) = a.get_mut(k as usize) {
                    x.endpoint = (itg + 1) as u16;
                }
            }
            if let Some(e) = endpoint.as_deref_mut() {
                if let Some(x) = e.get_mut(k as usize) {
                    *x = (itg + 1) as u16;
                }
            }
            len_tg -= 1;
            j += 1;
            i += 1;
        }
        itg += 1;
    }
    if i as i32 != ti.n_num_endpoints {
        return RI_ERR_PROGR;
    }
    0
}

/// `FillOutpStructEndpointFromInChI`
pub fn fill_outp_struct_endpoint_from_inchi(
    inchi: &InChI,
    p_endpoint: &mut Option<Vec<u16>>,
) -> i32 {
    let num_at = inchi.n_number_of_atoms.max(0) as usize;
    let mut endpoint = p_endpoint.take().unwrap_or_else(|| vec![0; num_at]);
    if endpoint.len() < num_at {
        endpoint.resize(num_at, 0);
    }
    for x in endpoint.iter_mut().take(num_at) {
        *x = 0;
    }
    if !(inchi.len_tautomer <= 1 || inchi.n_tautomer.is_empty()) {
        let nt = |j: usize| inchi.n_tautomer.get(j).copied().unwrap_or(0);
        let mut j = 1usize;
        for itg in 0..nt(0) as usize {
            let mut len_tg = i32::from(nt(j));
            j += T_GROUP_HDR_LEN;
            len_tg -= T_GROUP_HDR_LEN as i32 - 1;
            while 0 < len_tg {
                len_tg -= 1;
                let k = (i32::from(nt(j)) - 1) as usize;
                if let Some(x) = endpoint.get_mut(k) {
                    *x = (itg + 1) as u16;
                }
                j += 1;
            }
        }
    }
    *p_endpoint = Some(endpoint);
    0
}

/// `cmp_charge_val`
pub fn cmp_charge_val(p1: &ChargeVal, p2: &ChargeVal) -> i32 {
    let mut diff = p1.n_valence - p2.n_valence; // smaller valence first
    if diff != 0 {
        return diff;
    }
    diff = p1.n_charge.abs() - p2.n_charge.abs(); // smaller abs charge first
    if diff != 0 {
        return diff;
    }
    diff = p2.n_charge - p1.n_charge; // (+) first, (-) second
    if diff != 0 {
        return diff;
    }
    p1.n_valence_ordering_number - p2.n_valence_ordering_number
}

/// `bMayBeACationInMobileHLayer`
pub fn b_may_be_a_cation_in_mobile_h_layer(
    at: &[InpAtom],
    p_va: &[ValAt],
    iat: usize,
    b_mobile_h: i32,
) -> bool {
    if b_mobile_h == 0 || at[iat].num_h == 0 {
        return true;
    }
    let c_val = match at[iat].el_number {
        EL_NUMBER_N | EL_NUMBER_P => 4,
        EL_NUMBER_O | EL_NUMBER_S | EL_NUMBER_SE | EL_NUMBER_TE => 3,
        _ => return true,
    };
    if i32::from(at[iat].valence) + i32::from(at[iat].num_h) <= c_val {
        for j in 0..at[iat].valence.max(0) as usize {
            let neigh = at[iat].neighbor[j] as usize;
            if at[neigh].valence == 4
                && at[neigh].chem_bonds_valence == 4
                && at[neigh].num_h == 0
                && p_va[neigh].c_num_valence_electrons == 3
                && p_va[neigh].c_periodic_row_number == 1
            {
                return true;
            }
        }
        return false;
    }
    true
}

/// `clean_charge_val`: keeps only the suitable charge/valence states;
/// returns their number
#[allow(clippy::too_many_arguments)]
pub fn clean_charge_val(
    p_charge_val: &mut [ChargeVal],
    len: usize,
    atom: &[InpAtom],
    p_va: &[ValAt],
    iat: usize,
    b_is_metal: bool,
    b_mobile_h: i32,
    endpoint: Option<&[u16]>,
) -> usize {
    let at = &atom[iat];
    let n_periodic_num = i32::from(at.el_number);
    let num_bonds = i32::from(at.valence);
    let min_valence = i32::from(at.valence) + i32::from(at.num_h);
    let b_tautomeric = at.endpoint != 0;
    let b_fixed_h_tautomeric = b_mobile_h == 0
        && endpoint.is_some_and(|e| e.get(iat).copied().unwrap_or(0) != 0)
        && p_va[iat].c_num_valence_electrons == 6
        && num_bonds == 1
        && at.num_h == 0
        && !b_is_metal;
    let b_do_not_add_h = if_skip_add_h(n_periodic_num);
    let mut n_period = 0;
    let n_num_valence_el = get_sp_element_type(n_periodic_num, &mut n_period) - 1;
    if len == 0 {
        return len;
    }
    insertions_sort(&mut p_charge_val[..len], cmp_charge_val);
    if b_is_metal && b_do_not_add_h {
        return len.min(1);
    }
    let mut n_num_eq_abs_charges = 0;
    let mut j = 0usize;
    let mut i = 0usize;
    while i < len
        && (j as i32)
            < if n_num_eq_abs_charges != 0 {
                3 + n_num_eq_abs_charges
            } else {
                4
            }
    {
        let cur = p_charge_val[i];
        let skip = 'blk: {
            if cur.n_charge.abs() > 1 {
                break 'blk true;
            }
            if BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * (min_valence - 1) < cur.n_valence {
                break 'blk true;
            }
            if (b_tautomeric || (j != 0 && b_fixed_h_tautomeric)) && cur.n_charge < 0 {
                break 'blk true;
            }
            if (b_tautomeric || b_fixed_h_tautomeric) && cur.n_charge > 0 {
                break 'blk true;
            }
            if j != 0
                && !b_may_be_a_cation_in_mobile_h_layer(atom, p_va, iat, b_mobile_h)
                && cur.n_charge > 0
            {
                if i + 1 < len
                    && cur.n_valence == p_charge_val[i + 1].n_valence
                    && cur.n_charge == -p_charge_val[i + 1].n_charge
                {
                    i += 1;
                }
                break 'blk true;
            }
            if j != 0
                && !b_tautomeric
                && cur.n_valence == p_charge_val[j - 1].n_valence
                && cur.n_charge == -p_charge_val[j - 1].n_charge
            {
                if n_num_valence_el == VALUE_OCTET / 2
                    && cur.n_charge != 0
                    && n_num_eq_abs_charges == 0
                {
                    p_charge_val[j] = cur;
                    j += 1;
                    n_num_eq_abs_charges += 1;
                }
                break 'blk true;
            }
            if n_num_valence_el == 5
                && n_period == 1
                && at.num_h != 0
                && j != 0
                && b_mobile_h == 0
                && cur.n_valence == 5
                && cur.n_charge == 0
            {
                break 'blk true;
            }
            if j != 0 && cur.n_valence > p_charge_val[j - 1].n_valence + 1 {
                // do not accept gaps in allowed valences
                return finish_clean(j, n_num_eq_abs_charges, num_bonds);
            }
            false
        };
        if !skip {
            p_charge_val[j] = cur;
            j += 1;
        }
        i += 1;
    }
    finish_clean(j, n_num_eq_abs_charges, num_bonds)
}

fn finish_clean(len: usize, n_num_eq_abs_charges: i32, num_bonds: i32) -> usize {
    if n_num_eq_abs_charges == 0 && num_bonds < 3 && len == 4 {
        len - 1 // prohibit =S# where # is a triple bond
    } else {
        len
    }
}

/// `get_sp_element_type`: the s or p-element type (and the row)
pub fn get_sp_element_type(n_periodic_number: i32, n_row: &mut i32) -> i32 {
    let (row, type_);
    if n_periodic_number == 1 {
        type_ = 1;
        row = 0;
    } else if n_periodic_number == 2 {
        type_ = 0;
        row = 0;
    } else if n_periodic_number <= 10 {
        type_ = n_periodic_number - 1;
        row = 1;
    } else if n_periodic_number <= 18 {
        type_ = n_periodic_number - 9;
        row = 2;
    } else if n_periodic_number <= 20 {
        type_ = n_periodic_number - 17;
        row = 3;
    } else if n_periodic_number <= 30 {
        type_ = 0;
        row = 3;
    } else if n_periodic_number <= 36 {
        type_ = n_periodic_number - 27;
        row = 3;
    } else if n_periodic_number <= 38 {
        type_ = n_periodic_number - 35;
        row = 4;
    } else if n_periodic_number <= 48 {
        type_ = 0;
        row = 4;
    } else if n_periodic_number <= 54 {
        type_ = n_periodic_number - 45;
        row = 4;
    } else if n_periodic_number <= 56 {
        type_ = n_periodic_number - 53;
        row = 5;
    } else if n_periodic_number <= 80 {
        type_ = 0;
        row = 5;
    } else if n_periodic_number <= 86 {
        type_ = n_periodic_number - 77;
        row = 5;
    } else if n_periodic_number <= 88 {
        type_ = n_periodic_number - 85;
        row = 6;
    } else {
        type_ = 0;
        row = 6;
    }
    *n_row = row;
    if type_ == 9 { 0 } else { type_ }
}

/// `GetAtomRestoreInfo`: finds the charge structure of the atom `iat`.
/// Returns -1/`RI_ERR_*` on an error, 0 when the atom is to be left
/// unchanged, `TREAT_ATOM_AS_METAL`, or 1 on success.
pub fn get_atom_restore_info(
    atom: &mut [InpAtom],
    iat: usize,
    p_varray: &mut [ValAt],
    p_srm: &Srm,
    b_mobile_h: i32,
    endpoint: Option<&[u16]>,
) -> i32 {
    const MIN_ATOM_CHARGE: i32 = -2;
    const MAX_ATOM_CHARGE: i32 = 2;
    const NUM_ATOM_CHARGES: usize = (MAX_ATOM_CHARGE - MIN_ATOM_CHARGE + 1) as usize;
    const MAX_NUM_VALENCES: usize = 5;
    let n_periodic_num = i32::from(atom[iat].el_number);
    let mut charge_val = [ChargeVal::default(); NUM_ATOM_CHARGES * MAX_NUM_VALENCES];

    p_varray[iat].c_do_not_add_h = i8::from(if_skip_add_h(n_periodic_num));

    let (
        mut metal_bonds_chem_valence,
        mut not_metal_bonds_chem_valence,
        mut alt_bonds_delta_valence,
    ) = (0i32, 0i32, 0i32);
    let valence = atom[iat].valence.max(0) as usize;
    let j: i32;
    if p_varray[iat].c_metal != 0 {
        j = i32::from(atom[iat].valence);
        for i in 0..valence {
            let bond_type = i32::from(atom[iat].bond_type[i]) & BOND_TYPE_MASK;
            if bond_type <= BOND_TYPE_TRIPLE {
                metal_bonds_chem_valence += BOND_TYPE_SINGLE.max(bond_type);
            } else {
                metal_bonds_chem_valence += BOND_TYPE_SINGLE;
            }
        }
    } else {
        let mut jj = 0;
        for i in 0..valence {
            let bond_type = i32::from(atom[iat].bond_type[i]) & BOND_TYPE_MASK;
            let neigh = atom[iat].neighbor[i] as usize;
            if p_varray[neigh].c_metal != 0 {
                jj += 1;
                if bond_type <= BOND_TYPE_TRIPLE {
                    metal_bonds_chem_valence += BOND_TYPE_SINGLE.max(bond_type);
                } else {
                    metal_bonds_chem_valence += BOND_TYPE_SINGLE;
                }
            } else if bond_type <= BOND_TYPE_TRIPLE {
                not_metal_bonds_chem_valence += BOND_TYPE_SINGLE.max(bond_type);
            } else {
                not_metal_bonds_chem_valence += BOND_TYPE_SINGLE;
            }
        }
        j = jj;
    }
    let bonds_chem_valence = metal_bonds_chem_valence + not_metal_bonds_chem_valence;
    if i32::from(atom[iat].chem_bonds_valence) > bonds_chem_valence {
        if i32::from(atom[iat].chem_bonds_valence) - bonds_chem_valence > 1 {
            atom[iat].chem_bonds_valence = (bonds_chem_valence + 1) as i8; // should not happen
        }
        alt_bonds_delta_valence = i32::from(atom[iat].chem_bonds_valence) - bonds_chem_valence;
    }
    p_varray[iat].c_num_bonds_to_metal = j as i8;
    if n_periodic_num == 1 {
        return 0; // ignore bridging H
    }
    let num_h = i32::from(atom[iat].num_h);
    let num_bonds = i32::from(atom[iat].valence);
    if num_bonds == 0 && num_h == 0 {
        return 0; // isolated atom
    }
    let min_chem_valence = num_bonds + num_h;
    let cur_chem_valence = bonds_chem_valence + alt_bonds_delta_valence + num_h;
    let pva = &mut p_varray[iat];
    let mut known_chem_valence = get_el_valence(n_periodic_num, 0, 0);
    let cur_chem_valence_fixed;
    if p_srm.b_metal_add_flower != 0 {
        cur_chem_valence_fixed = cur_chem_valence
            - i32::from(pva.c_num_bonds_to_metal) * (1 - p_srm.n_metal_init_bond_order);
        pva.c_init_orig_valence_to_metal = metal_bonds_chem_valence as i8;
        pva.c_init_valence_to_metal = (metal_bonds_chem_valence
            - i32::from(pva.c_num_bonds_to_metal) * (1 - p_srm.n_metal_init_bond_order))
            as i8;
        pva.c_init_flow_to_metal = (i32::from(pva.c_init_valence_to_metal)
            - i32::from(pva.c_num_bonds_to_metal) * p_srm.n_metal_min_bond_order)
            as i8;
        if pva.c_metal != 0 {
            pva.c_init_free_valences =
                (i32::from(pva.c_init_free_valences) + alt_bonds_delta_valence) as i8;
        }
        if p_srm.n_metal_init_edge_flow
            < p_srm.n_metal_init_bond_order - p_srm.n_metal_min_bond_order
        {
            if pva.c_init_flow_to_metal <= pva.c_num_bonds_to_metal {
                if pva.c_metal != 0 {
                    pva.c_init_free_valences = (i32::from(pva.c_init_free_valences)
                        + i32::from(pva.c_init_flow_to_metal))
                        as i8;
                }
                pva.c_init_flow_to_metal = 0;
            } else {
                if pva.c_metal != 0 {
                    pva.c_init_free_valences = (i32::from(pva.c_init_free_valences)
                        + i32::from(pva.c_num_bonds_to_metal) * (1 - p_srm.n_metal_init_edge_flow))
                        as i8;
                }
                pva.c_init_flow_to_metal = (i32::from(pva.c_init_flow_to_metal)
                    - i32::from(pva.c_num_bonds_to_metal) * (1 - p_srm.n_metal_init_edge_flow))
                    as i8;
            }
        }
    } else {
        cur_chem_valence_fixed = cur_chem_valence;
        pva.c_init_flow_to_metal =
            (metal_bonds_chem_valence - i32::from(pva.c_num_bonds_to_metal)) as i8;
        pva.c_init_valence_to_metal = metal_bonds_chem_valence as i8;
        pva.c_init_orig_valence_to_metal = metal_bonds_chem_valence as i8;
    }
    if pva.c_metal != 0 && p_srm.b_metal_add_flower != 0 {
        pva.c_cn_list_index = (CN_LIST_INDEX_ME + 1) as i8;
        return 0; // metal
    }
    if known_chem_valence == 0 {
        // a noble gas like He, Ne, ...
        pva.c_init_free_valences =
            (i32::from(atom[iat].chem_bonds_valence) - i32::from(atom[iat].valence)) as i8;
        return TREAT_ATOM_AS_METAL;
    }
    let _ = known_chem_valence;
    let mut n_found_neutral_valence_excess = NO_VALUE_INT;
    let mut n_last_found_valence_state = NO_VALUE_INT;
    let mut n_num_states = 0usize;
    for cur_charge in MIN_ATOM_CHARGE..=MAX_ATOM_CHARGE {
        for i in 0..MAX_NUM_VALENCES {
            known_chem_valence = get_el_valence(n_periodic_num, cur_charge, i);
            if cur_chem_valence_fixed > known_chem_valence || known_chem_valence == 0 {
                continue;
            }
            if BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * (num_bonds - 1) + num_h < known_chem_valence {
                continue;
            }
            charge_val[n_num_states] = ChargeVal {
                n_valence: known_chem_valence,
                n_charge: cur_charge,
                n_valence_ordering_number: i as i32,
            };
            if cur_charge == 0 && n_found_neutral_valence_excess == NO_VALUE_INT {
                n_found_neutral_valence_excess = known_chem_valence - min_chem_valence;
            }
            if min_chem_valence == known_chem_valence {
                if n_last_found_valence_state == NO_VALUE_INT {
                    n_last_found_valence_state = n_num_states as i32;
                } else if charge_val[n_last_found_valence_state as usize]
                    .n_charge
                    .abs()
                    >= cur_charge.abs()
                {
                    n_last_found_valence_state = n_num_states as i32;
                }
            }
            n_num_states += 1;
        }
    }
    let _ = (n_found_neutral_valence_excess, n_last_found_valence_state);

    // select only appropriate charge & valence so that a suitable ChargeStruct exists
    let c_metal = p_varray[iat].c_metal;
    let mut n_num_selected_states = clean_charge_val(
        &mut charge_val,
        n_num_states,
        atom,
        p_varray,
        iat,
        c_metal != 0,
        b_mobile_h,
        endpoint,
    ) as i32;
    if n_num_selected_states == 0 {
        return TREAT_ATOM_AS_METAL;
    }

    // find an appropriate ChargeStruct index for the ChargeVal found
    let mut cn_bits_array = [0i32; 5];
    let mut len_cn_bits_array = 0usize;
    let mut i = 0i32;
    while i < n_num_selected_states && len_cn_bits_array < 4 {
        match charge_val[i as usize].n_charge {
            -1 => cn_bits_array[len_cn_bits_array] |= CN_BITS_M,
            0 => cn_bits_array[len_cn_bits_array] |= CN_BITS_N,
            1 => cn_bits_array[len_cn_bits_array] |= CN_BITS_P,
            _ => return RI_ERR_PROGR,
        }
        if i + 1 < n_num_selected_states
            && charge_val[i as usize].n_valence == charge_val[i as usize + 1].n_valence
            && charge_val[i as usize].n_charge != 0
            && charge_val[i as usize].n_charge == -charge_val[i as usize + 1].n_charge
        {
            // add the opposite charge to the same element of cn_bits_array[]
        } else {
            len_cn_bits_array += 1;
        }
        i += 1;
    }
    if len_cn_bits_array == 0 || len_cn_bits_array > 4 {
        return RI_ERR_PROGR;
    }
    // accommodate the added 4-state ChargeStruct: +/- cannot be in case of 4 states
    if len_cn_bits_array as i32 + 1 == n_num_selected_states && n_num_selected_states == 4 {
        len_cn_bits_array -= 1;
        n_num_selected_states -= 1;
        cn_bits_array[len_cn_bits_array] = 0;
    }
    // terminal hydrogenless -C as in isocyano or CO: there is no just cnE_[] ChargeStruct
    if len_cn_bits_array == 1
        && cn_bits_array[0] == (CN_BITS_P | CN_BITS_M)
        && charge_val[0].n_valence + 1
            > BOND_TYPE_TRIPLE + BOND_TYPE_DOUBLE * (num_bonds - 1) + num_h
    {
        cn_bits_array[len_cn_bits_array] = CN_BITS_N;
        len_cn_bits_array += 1;
        let prev = charge_val[n_num_selected_states as usize - 1].n_valence;
        charge_val[n_num_selected_states as usize] = ChargeVal {
            n_valence: prev,
            n_charge: 0,
            n_valence_ordering_number: 0,
        };
    }
    let found = loop {
        let cn_bits = make_cn_bits(
            cn_bits_array[0],
            cn_bits_array[1],
            cn_bits_array[2],
            cn_bits_array[3],
        );
        if let Some(j) = CN_LIST.iter().position(|c| c.bits == cn_bits) {
            break j;
        }
        // the ChargeStructure was not found
        if 1 < len_cn_bits_array && len_cn_bits_array as i32 + 1 == n_num_selected_states {
            len_cn_bits_array -= 1;
            cn_bits_array[len_cn_bits_array] = 0;
        } else if n_num_selected_states == 4 {
            len_cn_bits_array -= 1;
            cn_bits_array[len_cn_bits_array] = 0;
            n_num_selected_states -= 1;
        } else {
            return RI_ERR_PROGR;
        }
    };
    let pva = &mut p_varray[iat];
    pva.c_cn_list_index = (found + 1) as i8;
    pva.c_init_charge = CN_LIST[found].n_initial_charge as i8;
    // the "free valence"
    let jfree = charge_val[0].n_valence - cur_chem_valence_fixed;
    if jfree < 0 {
        return RI_ERR_PROGR;
    }
    pva.c_init_free_valences = jfree as i8;
    1
}

// ---------------------------------------------------------------------------
// Small helpers of the restoration: t/c groups, bond flows, edge lists

/// `RegisterTCGroup`: adds the capacities to the group of the type and
/// ordinal number, creating it when absent; returns `index + 1` for a new
/// group, 0 for an existing one
pub fn register_tc_group(
    tcg: &mut AllTcGroups,
    group_type: i32,
    group_ord_num: i32,
    vertex_cap: i32,
    vertex_flow: i32,
    edge_cap: i32,
    edge_flow: i32,
    num_edges: i32,
) -> i32 {
    let mut ret = 0;
    let n = tcg.num_tc_groups as usize;
    let mut i = tcg.p_tcg[..n.min(tcg.p_tcg.len())]
        .iter()
        .position(|g| g.type_ == group_type && g.ord_num == group_ord_num)
        .unwrap_or(n);
    if i == n {
        if tcg.num_tc_groups == tcg.max_tc_groups {
            // ReallocTCGroups
            tcg.p_tcg.resize(
                tcg.max_tc_groups as usize + INC_NUM_TCGROUPS,
                TcGroup::default(),
            );
            tcg.max_tc_groups += INC_NUM_TCGROUPS as i32;
        }
        if tcg.p_tcg.len() <= i {
            tcg.p_tcg.resize(i + 1, TcGroup::default());
        }
        ret = i as i32 + 1;
        tcg.num_tc_groups += 1;
        tcg.p_tcg[i].type_ = group_type;
        tcg.p_tcg[i].ord_num = group_ord_num;
    }
    i = i.min(tcg.p_tcg.len() - 1);
    let g = &mut tcg.p_tcg[i];
    g.num_edges += num_edges;
    g.st_cap += vertex_cap;
    g.st_flow += vertex_flow;
    g.edges_cap += edge_cap;
    g.edges_flow += edge_flow;
    ret
}

/// `nTautEndpointEdgeCap`: the capacity of the t-group edge of an endpoint
pub fn n_taut_endpoint_edge_cap(at: &[InpAtom], p_va: &[ValAt], i: usize) -> i32 {
    let mut cap = i32::from(p_va[i].c_init_free_valences);
    if p_va[i].c_cn_list_index > 0 {
        let cn = &CN_LIST[p_va[i].c_cn_list_index as usize - 1].p_cn;
        cap += i32::from(cn[0].v.cap) - i32::from(cn[0].v.flow);
    }
    let mut excess = 0;
    let a = &at[i];
    let mut j = 0;
    while j < MAX_NUM_STEREO_BONDS && a.sb_parity[j] != 0 {
        let k = a.sb_ord[j] as usize;
        if i32::from(a.bond_type[k]) < BOND_TYPE_TRIPLE {
            excess += i32::from(a.bond_type[k]) - BOND_TYPE_SINGLE;
        }
        j += 1;
    }
    let bonds_valence = (i32::from(a.chem_bonds_valence) - i32::from(a.valence)) - excess;
    if bonds_valence >= 0 {
        cap + bonds_valence
    } else {
        RI_ERR_PROGR
    }
}

/// `BondFlowMaxcapMinorder`: returns the initial flow; `(maxcap, minorder,
/// needs_flower)` go to the tuple
pub fn bond_flow_maxcap_minorder(
    atom: &[InpAtom],
    p_va: &[ValAt],
    srm: &Srm,
    iat: usize,
    ineigh: usize,
) -> (i32, i32, i32, bool) {
    let at = &atom[iat];
    let neigh = at.neighbor[ineigh] as usize;
    let mut bond_type = i32::from(at.bond_type[ineigh]) & BOND_TYPE_MASK;
    let n_metal = i32::from(p_va[iat].c_metal != 0) + i32::from(p_va[neigh].c_metal != 0);
    let n_endpoint = i32::from(at.endpoint != 0) + i32::from(atom[neigh].endpoint != 0);
    let n_stereo = i32::from(at.p_parity != 0 || at.sb_parity[0] != 0)
        + i32::from(atom[neigh].p_parity != 0 || atom[neigh].sb_parity[0] != 0);
    if bond_type > BOND_TYPE_TRIPLE {
        bond_type = BOND_TYPE_SINGLE;
    }
    let (min_order, mut flow, mut needs_flower);
    needs_flower = false;
    if (n_stereo != 0 && srm.b_fix_stereo_bonds != 0) || n_metal == 0 || srm.b_metal_add_flower == 0
    {
        min_order = BOND_TYPE_SINGLE;
        flow = bond_type - min_order;
    } else if n_metal != 0 && n_endpoint == 0 {
        min_order = srm.n_metal_min_bond_order;
        let init = srm.n_metal_init_bond_order + bond_type - BOND_TYPE_SINGLE;
        flow = init - min_order;
        if srm.n_metal_init_edge_flow == 0
            && srm.n_metal_init_bond_order > srm.n_metal_min_bond_order
            && flow > 0
        {
            flow -= 1;
        }
        needs_flower = p_va[iat].c_metal != 0;
    } else {
        min_order = srm.n_metal2_endpoint_min_bond_order;
        let init = srm.n_metal2_endpoint_init_bond_order + bond_type - BOND_TYPE_SINGLE;
        flow = init - min_order;
        if srm.n_metal2_endpoint_init_edge_flow == 0
            && srm.n_metal2_endpoint_init_bond_order > srm.n_metal2_endpoint_min_bond_order
            && flow > 0
        {
            flow -= 1;
        }
        needs_flower = if (p_va[iat].c_metal != 0
            && at.endpoint == 0
            && p_va[neigh].c_metal == 0
            && atom[neigh].endpoint != 0)
            || (p_va[neigh].c_metal != 0
                && atom[neigh].endpoint == 0
                && p_va[iat].c_metal == 0
                && at.endpoint != 0)
        {
            p_va[iat].c_metal != 0
        } else {
            p_va[iat].c_metal != 0 && at.endpoint == 0
        };
    }
    (flow, BOND_TYPE_TRIPLE - min_order, min_order, needs_flower)
}

/// `AtomStcapStflow`: `(num_flowers, stcap, stflow, mgroup_edge_cap,
/// mgroup_edge_flow)`
pub fn atom_stcap_stflow(
    atom: &[InpAtom],
    p_va: &[ValAt],
    srm: &Srm,
    iat: usize,
) -> (i32, i32, i32, i32, i32) {
    let mut stflow = 0;
    let mut stcap = i32::from(atom[iat].chem_bonds_valence);
    let (mut mcap, mut mflow, mut needs) = (0, 0, 0);
    if srm.b_metal_add_flower != 0 {
        stcap -= i32::from(p_va[iat].c_init_orig_valence_to_metal)
            - i32::from(p_va[iat].c_init_valence_to_metal);
    }
    for ineigh in 0..atom[iat].valence as usize {
        let (flow, _maxcap, minorder, flower) =
            bond_flow_maxcap_minorder(atom, p_va, srm, iat, ineigh);
        stflow += flow;
        stcap -= minorder;
        if flower {
            needs += 1;
            mflow += flow;
            mcap += BOND_TYPE_TRIPLE - minorder + srm.n_metal_max_charge_d;
        }
    }
    (
        needs,
        if needs != 0 { stflow } else { stcap },
        stflow,
        mcap,
        mcap - mflow,
    )
}

/// `comp_cc_cand` (with its unfinished "more valence electrons" branch)
pub fn comp_cc_cand(p1: &CcCand, p2: &CcCand) -> i32 {
    let d = |a: i32, b: i32| a - b;
    let mut ret = d(i32::from(p2.c_metal), i32::from(p1.c_metal));
    if ret != 0 {
        return ret;
    }
    ret = d(
        i32::from(p2.c_num_bonds_to_metal),
        i32::from(p1.c_num_bonds_to_metal),
    );
    if ret != 0 {
        return ret;
    }
    ret = d(
        i32::from(p2.c_periodic_row_number),
        i32::from(p1.c_periodic_row_number),
    );
    if ret != 0 {
        return ret;
    }
    ret = d(i32::from(p2.num_bonds), i32::from(p1.num_bonds));
    if ret != 0 {
        return ret;
    }
    ret = d(i32::from(p1.chem_valence), i32::from(p2.chem_valence));
    if ret != 0 {
        return ret;
    }
    if p1.c_num_valence_electrons == 0 && p2.c_num_valence_electrons != 0 {
        return -1;
    }
    if p2.c_num_valence_electrons == 0 && p1.c_num_valence_electrons != 0 {
        return -1;
    }
    if p2.c_num_valence_electrons != p1.c_num_valence_electrons {
        return ret; // sic: ret is 0 here
    }
    d(i32::from(p2.iat), i32::from(p1.iat))
}

pub const EL_TYPE_O: i32 = 0x0001;
pub const EL_TYPE_S: i32 = 0x0002;
pub const EL_TYPE_N: i32 = 0x0004;
pub const EL_TYPE_P: i32 = 0x0008;
pub const EL_TYPE_C: i32 = 0x0010;
pub const EL_TYPE_X: i32 = 0x0020;
pub const EL_TYPE_MASK: i32 = 0x003f;
pub const EL_TYPE_OST: i32 = 0x0100;
pub const EL_TYPE_PT: i32 = 0x0200;

/// `get_pVA_atom_type`
pub fn get_pva_atom_type(p_va: &[ValAt], at: &[InpAtom], iat: usize, bond_type: i32) -> i32 {
    let mut ty = 0;
    let v = &p_va[iat];
    match v.c_num_valence_electrons {
        4 => {
            if v.c_periodic_row_number == 1 {
                ty |= EL_TYPE_C;
            }
        }
        6 => {
            if v.c_periodic_row_number == 1 {
                ty |= EL_TYPE_O;
            } else if v.c_periodic_row_number < 5 {
                ty |= EL_TYPE_S;
            }
            if bond_type == BOND_TYPE_SINGLE
                && ty & (EL_TYPE_O | EL_TYPE_S) != 0
                && n_no_metal_bonds_valence(at, iat) == 1
                && n_no_metal_num_bonds(at, iat) == 1
            {
                ty |= EL_TYPE_OST;
            }
        }
        5 => {
            ty |= if v.c_periodic_row_number == 1 {
                EL_TYPE_N
            } else {
                EL_TYPE_P
            };
        }
        _ => {
            if !is_el_a_metal(i32::from(v.c_periodic_number)) {
                ty |= EL_TYPE_X;
            }
        }
    }
    let val = get_endpoint_valence(at[iat].el_number);
    let a = &at[iat];
    if val != 0
        && val > i32::from(a.valence)
        && a.radical == 0
        && (-1..=0).contains(&a.charge)
        && val == i32::from(a.chem_bonds_valence) - i32::from(a.charge) + i32::from(a.num_h)
    {
        ty |= EL_TYPE_PT;
    }
    ty
}

/// `EDGE_LIST`
#[derive(Debug, Clone, Default)]
pub struct EdgeList {
    pub num_alloc: i32,
    pub num_edges: i32,
    pub pn_edges: Vec<i32>,
}

pub const EDGE_LIST_CLEAR: i32 = -1;
pub const EDGE_LIST_FREE: i32 = -2;

/// `AllocEdgeList`
pub fn alloc_edge_list(e: &mut EdgeList, len: i32) -> i32 {
    match len {
        EDGE_LIST_FREE | EDGE_LIST_CLEAR => *e = EdgeList::default(),
        _ => {
            if len > 0 && len != e.num_alloc {
                let keep = e.num_edges.min(len);
                let mut v = vec![0; len as usize];
                if keep > 0 {
                    v[..keep as usize].copy_from_slice(&e.pn_edges[..keep as usize]);
                    e.num_edges = keep;
                } else {
                    e.num_edges = 0;
                }
                e.pn_edges = v;
                e.num_alloc = len;
            }
        }
    }
    0
}

/// `AddToEdgeList`
pub fn add_to_edge_list(e: &mut EdgeList, iedge: i32, add_len: i32) -> i32 {
    if e.num_alloc == e.num_edges {
        if add_len <= 0 {
            return RI_ERR_PROGR;
        }
        let ret = alloc_edge_list(e, e.num_alloc + add_len);
        if ret != 0 {
            return ret;
        }
    }
    e.pn_edges[e.num_edges as usize] = iedge;
    e.num_edges += 1;
    0
}

/// `RemoveFromEdgeListByIndex`
pub fn remove_from_edge_list_by_index(e: &mut EdgeList, index: i32) -> i32 {
    let len = e.num_edges - index - 1;
    if len >= 0 {
        let i = index as usize;
        e.pn_edges.copy_within(i + 1..i + 1 + len as usize, i);
        e.num_edges -= 1;
        e.pn_edges[e.num_edges as usize] = 0;
        0
    } else {
        -1
    }
}

/// `FindInEdgeList`
pub fn find_in_edge_list(e: &EdgeList, iedge: i32) -> i32 {
    (0..e.num_edges)
        .rev()
        .find(|&i| e.pn_edges[i as usize] == iedge)
        .unwrap_or(-1)
}

/// `RemoveFromEdgeListByValue`
pub fn remove_from_edge_list_by_value(e: &mut EdgeList, iedge: i32) -> i32 {
    let mut n = 0;
    for i in (0..e.num_edges).rev() {
        if e.pn_edges[i as usize] == iedge {
            let ret = remove_from_edge_list_by_index(e, i);
            if ret != 0 {
                return ret;
            }
            n += 1;
        }
    }
    n
}

// ---------------------------------------------------------------------------
// Sizes of the BNS and the super c-groups

const MAX_CN_VAL: usize = 3;

fn is_bns_vt_c_gr(x: u16) -> bool {
    x & BNS_VT_C_POS_ALL == BNS_VERT_TYPE_C_GROUP
}
fn is_bns_vt_m_gr(x: u16) -> bool {
    x == BNS_VERT_TYPE_METAL_GR
}
fn is_bns_vt_chrg_struct(x: u16) -> bool {
    x & BNS_VERT_TYPE_AUX != 0 && x & BNS_VERT_TYPE_TEMP != 0
}

/// `nCountBnsSizes`: fills the totals and registers the t- and c-groups;
/// returns 0, a negative error or (as the C code does) a positive group
/// number left by the last registration
pub fn n_count_bns_sizes(
    at: &[InpAtom],
    num_at: usize,
    ti: &TGroupInfo,
    p_va: &[ValAt],
    srm: &Srm,
    tcg: &mut AllTcGroups,
) -> i32 {
    let (mut n_bonds, mut n_other_edges, mut n_vertices) = (0i32, 0i32, 0i32);
    let (mut n_tgroup_edges, mut n_tgroup_edges_from_tg, mut n_tot_neg_charge) = (0i32, 0i32, 0i32);
    let mut ret = 0;
    for i in 0..num_at {
        tcg.num_metal_atoms += i32::from(p_va[i].c_metal != 0);
        tcg.num_metal_bonds += i32::from(p_va[i].c_num_bonds_to_metal);
        tcg.total_electrons += i32::from(at[i].el_number);
        tcg.total_electrons_metals += if p_va[i].c_metal != 0 {
            i32::from(at[i].el_number)
        } else {
            0
        };
    }
    tcg.total_electrons -= tcg.total_charge;
    tcg.num_metal_bonds /= 2;

    if let Some(tgs) = ti.t_group.as_ref() {
        for g in tgs.iter().take(ti.num_t_groups) {
            ret = register_tc_group(
                tcg,
                i32::from(BNS_VERT_TYPE_TGROUP),
                i32::from(g.n_group_number),
                i32::from(g.num[0]),
                0,
                0,
                0,
                i32::from(g.n_num_endpoints),
            );
            if ret < 0 {
                return ret;
            }
            n_other_edges += i32::from(g.n_num_endpoints);
            n_tgroup_edges_from_tg += i32::from(g.n_num_endpoints);
            n_tot_neg_charge += i32::from(g.num[1]);
            if ret > 0 {
                let j = ret as usize - 1;
                tcg.p_tcg[j].tg_num_h = (i32::from(g.num[0]) - i32::from(g.num[1])) as i16;
                tcg.p_tcg[j].tg_num_minus = g.num[1] as i16;
            }
        }
    }

    let mut b_metal_atoms = false;
    loop {
        for i in 0..num_at {
            n_bonds += i32::from(at[i].valence);
            if p_va[i].c_cn_list_index != 0 {
                let cl = &CN_LIST[p_va[i].c_cn_list_index as usize - 1];
                let (len, bits, p_cn) = (cl.len as usize, cl.bits, cl.p_cn);
                if (bits != CN_BITS_ME) != !b_metal_atoms {
                    continue;
                }
                for j in 0..len {
                    let ty = p_cn[j].v.type_;
                    let mut k = 0;
                    while k < MAX_CN_VAL && p_cn[j].e[k].neigh != 0 {
                        n_other_edges += 1;
                        let n = p_cn[j].e[k].neigh as usize - 1;
                        let neigh_type = p_cn[n].v.type_;
                        let e = &p_cn[j].e[k];
                        let md = srm.n_metal_max_charge_d;
                        let cap = if !b_metal_atoms {
                            i32::from(e.cap)
                        } else if e.cap != 0 {
                            md
                        } else {
                            0
                        };
                        let flow = if !b_metal_atoms {
                            i32::from(e.flow)
                        } else if e.flow != 0 {
                            md
                        } else {
                            0
                        };
                        for t in [neigh_type, ty] {
                            if is_bns_vt_c_gr(t) {
                                ret = register_tc_group(tcg, i32::from(t), 0, 0, 0, cap, flow, 1);
                                if ret < 0 {
                                    return ret;
                                }
                                if ret > 0 {
                                    ret = register_tc_group(tcg, i32::from(t), 0, 0, 0, 0, 0, 1);
                                    if ret < 0 {
                                        return ret;
                                    }
                                    n_other_edges += 1;
                                }
                            }
                        }
                        k += 1;
                    }
                    if ty & BNS_VERT_TYPE_ATOM != 0 {
                        continue;
                    }
                    if is_bns_vt_chrg_struct(ty) {
                        n_vertices += 1;
                        continue;
                    }
                    if srm.b_metal_add_flower != 0 && is_bns_vt_m_gr(ty) {
                        let (needs, _, _, edge_cap, edge_flow) =
                            atom_stcap_stflow(at, p_va, srm, i);
                        if needs == 0 {
                            return RI_ERR_PROGR;
                        }
                        let t = i32::from(ty);
                        ret = register_tc_group(tcg, t, 0, 0, 0, edge_cap, edge_flow, 1);
                        if ret < 0 {
                            return ret;
                        }
                        n_other_edges += 1;
                        if ret > 0 {
                            ret = register_tc_group(tcg, t, 0, 0, 0, 0, 0, 2);
                            if ret < 0 {
                                return ret;
                            }
                            for (ord, ne) in [(1, 3), (2, 3), (3, 2)] {
                                ret = register_tc_group(tcg, t, ord, 0, 0, 0, 0, ne);
                                if ret <= 0 {
                                    return if ret == 0 { RI_ERR_PROGR } else { ret };
                                }
                            }
                            n_vertices += 4;
                            n_other_edges += 5;
                        }
                        continue;
                    }
                    n_vertices += 1;
                    if !is_bns_vt_c_gr(ty) {
                        return RI_ERR_PROGR;
                    }
                    let md = srm.n_metal_max_charge_d;
                    let v = &p_cn[j].v;
                    let cap = if !b_metal_atoms {
                        i32::from(v.cap)
                    } else if v.cap != 0 {
                        md
                    } else {
                        0
                    };
                    let flow = if !b_metal_atoms {
                        i32::from(v.flow)
                    } else if v.flow != 0 {
                        md
                    } else {
                        0
                    };
                    ret = register_tc_group(tcg, i32::from(ty), 0, cap, flow, 0, 0, 0);
                    if ret < 0 {
                        return ret;
                    }
                }
            }
            if at[i].endpoint != 0 {
                let edge_cap = n_taut_endpoint_edge_cap(at, p_va, i);
                n_tgroup_edges += 1;
                if edge_cap < 0 {
                    return edge_cap;
                }
                ret = register_tc_group(
                    tcg,
                    i32::from(BNS_VERT_TYPE_TGROUP),
                    i32::from(at[i].endpoint),
                    0,
                    0,
                    edge_cap,
                    0,
                    0,
                );
                if ret < 0 {
                    return ret;
                }
            }
        }
        if !b_metal_atoms && tcg.num_metal_atoms != 0 {
            b_metal_atoms = true;
            n_bonds = 0;
            continue;
        }
        break;
    }

    n_bonds /= 2;
    tcg.num_atoms = num_at as i32;
    tcg.num_bonds = n_bonds;
    tcg.num_tgroups = ti.num_t_groups as i32;
    tcg.num_tgroup_edges = n_tgroup_edges;
    tcg.tgroup_charge = -n_tot_neg_charge;
    if ret >= 0 && n_tgroup_edges_from_tg != n_tgroup_edges {
        ret = BNS_PROGRAM_ERR;
    }
    n_vertices += num_at as i32;
    n_vertices += ti.num_t_groups as i32;
    n_bonds += n_other_edges;
    tcg.n_vertices = n_vertices;
    tcg.n_edges = n_bonds;
    ret
}

/// `nAddSuperCGroups`
pub fn n_add_super_c_groups(tcg: &mut AllTcGroups) -> i32 {
    for i in 0..tcg.num_tc_groups as usize {
        let ty = tcg.p_tcg[i].type_ as u16;
        if ty & BNS_VERT_TYPE_TGROUP != 0 {
            continue;
        }
        if is_bns_vt_c_gr(ty) || is_bns_vt_m_gr(ty) {
            let k = match ty {
                BNS_VT_C_POS => TCG_PLUS0,
                BNS_VT_C_NEG => TCG_MINUS0,
                BNS_VT_C_POS_C => TCG_PLUS_C0,
                BNS_VT_C_NEG_C => TCG_MINUS_C0,
                BNS_VT_C_POS_M => TCG_PLUS_M0,
                BNS_VT_C_NEG_M => TCG_MINUS_M0,
                BNS_VT_M_GROUP => match tcg.p_tcg[i].ord_num {
                    0 => TCG_ME_FLOWER0,
                    1 => TCG_ME_FLOWER1,
                    2 => TCG_ME_FLOWER2,
                    3 => TCG_ME_FLOWER3,
                    _ => return RI_ERR_PROGR,
                },
                _ => return RI_ERR_PROGR,
            } as usize;
            if tcg.n_group[k] >= 0 || (tcg.p_tcg[i].ord_num != 0 && !is_bns_vt_m_gr(ty)) {
                return RI_ERR_PROGR;
            }
            tcg.n_group[k] = i as i32;
        }
    }
    // (+) supergroup
    let n = [TCG_PLUS0, TCG_PLUS_C0, TCG_PLUS_M0]
        .iter()
        .filter(|&&g| tcg.n_group[g as usize] >= 0)
        .count() as i32;
    if n != 0 {
        let ret = register_tc_group(tcg, i32::from(BNS_VT_C_POS_ALL), 0, 0, 0, 0, 0, 1 + n);
        if ret <= 0 {
            return if ret == 0 { RI_ERR_PROGR } else { ret };
        }
        tcg.n_group[TCG_PLUS as usize] = ret - 1;
        tcg.n_vertices += 2;
        tcg.n_edges += 1 + n;
    }
    // (-) supergroup
    let n = [TCG_MINUS0, TCG_MINUS_C0, TCG_MINUS_M0]
        .iter()
        .filter(|&&g| tcg.n_group[g as usize] >= 0)
        .count() as i32;
    if n != 0 {
        let ret = register_tc_group(tcg, i32::from(BNS_VT_C_NEG_ALL), 0, 0, 0, 0, 0, 1 + n);
        if ret < 0 {
            return ret;
        }
        tcg.n_group[TCG_MINUS as usize] = ret - 1;
        tcg.n_vertices += 2;
        tcg.n_edges += 1 + n;
    }
    let n = i32::from(tcg.n_group[TCG_MINUS as usize] >= 0)
        + i32::from(tcg.n_group[TCG_PLUS as usize] >= 0);
    if n != 0 {
        tcg.n_vertices += 1;
        tcg.n_edges += n;
    }
    0
}

// ---------------------------------------------------------------------------
// Building the BNS: t-groups, metal flower, super c-groups

const MAX_TGROUP_EDGE_CAP: i32 = 2;

/// `ConnectTwoVertices` on the vertices `ip1`, `ip2` and the edge `ie`
pub fn connect_two_vertices(
    p: &mut BnStruct,
    ip1: usize,
    ip2: usize,
    ie: usize,
    clear_edge: bool,
) -> i32 {
    let (v1, v2) = (p.vert[ip1], p.vert[ip2]);
    if ip1 as i32 >= p.max_vertices
        || ip2 as i32 >= p.max_vertices
        || ie as i32 >= p.max_edges
        || v1.iedge + v1.max_adj_edges as usize > p.max_iedges as usize
        || v2.iedge + v2.max_adj_edges as usize > p.max_iedges as usize
        || v1.num_adj_edges >= v1.max_adj_edges
        || v2.num_adj_edges >= v2.max_adj_edges
    {
        return BNS_VERT_EDGE_OVFL;
    }
    if clear_edge {
        p.edge[ie] = BnsEdge::default();
    } else if p.edge[ie].neighbor1 != 0 || p.edge[ie].neighbor12 != 0 {
        return BNS_PROGRAM_ERR;
    }
    p.edge[ie].neighbor1 = ip1.min(ip2) as u16;
    p.edge[ie].neighbor12 = (ip1 ^ ip2) as u16;
    let (n1, n2) = (v1.num_adj_edges, v2.num_adj_edges);
    p.iedge[v1.iedge + n1 as usize] = ie as i32;
    p.iedge[v2.iedge + n2 as usize] = ie as i32;
    p.edge[ie].neigh_ord[usize::from(ip1 > ip2)] = n1;
    p.vert[ip1].num_adj_edges += 1;
    p.edge[ie].neigh_ord[usize::from(ip1 < ip2)] = n2;
    p.vert[ip2].num_adj_edges += 1;
    0
}

/// `AddTGroups2TCGBnStruct`
pub fn add_tgroups_to_tcg_bn_struct(
    bns: &mut BnStruct,
    at: &[InpAtom],
    num_atoms: usize,
    p_va: &mut [ValAt],
    tcg: &mut AllTcGroups,
    max_add_edges: i32,
) -> i32 {
    let mut ret = 0;
    if tcg.num_tgroups == 0 {
        return ret;
    }
    let num_tg = tcg.num_tgroups as usize;
    let mut num_edges = bns.num_edges;
    let num_vertices = bns.num_vertices;
    if num_vertices + num_tg as i32 >= bns.max_vertices
        || num_edges + tcg.num_tgroup_edges >= bns.max_edges
    {
        return BNS_VERT_EDGE_OVFL;
    }
    let mut max_tg_number = 0;
    let mut i = 0usize;
    while i < tcg.num_tc_groups as usize {
        if tcg.p_tcg[i].type_ & i32::from(BNS_VERT_TYPE_TGROUP) != 0 {
            let k = tcg.p_tcg[i].ord_num;
            if k <= 0 || k > tcg.num_tc_groups || k != max_tg_number + 1 {
                return BNS_CPOINT_ERR;
            }
            max_tg_number = k;
        } else {
            break;
        }
        i += 1;
    }
    if i != num_tg {
        return BNS_CPOINT_ERR;
    }
    let nv = num_vertices as usize;
    for v in &mut bns.vert[nv..nv + max_tg_number as usize] {
        *v = BnsVertex::default();
    }
    let mut prev = nv - 1;
    let mut tot_st_cap = 0;
    for i in 0..num_tg {
        let cur = nv + tcg.p_tcg[i].ord_num as usize - 1;
        let (piedge, pmax) = (bns.vert[prev].iedge, bns.vert[prev].max_adj_edges);
        let g = tcg.p_tcg[i];
        let v = &mut bns.vert[cur];
        v.iedge = piedge + pmax as usize;
        v.max_adj_edges = (g.num_edges + max_add_edges + BNS_ADD_SUPER_TGROUP) as u16;
        v.num_adj_edges = 0;
        v.st_edge.flow = 0;
        v.st_edge.flow0 = 0;
        v.st_edge.cap = g.st_cap;
        v.st_edge.cap0 = g.st_cap;
        tot_st_cap += g.st_cap;
        v.type_ = g.type_ as u16;
        tcg.p_tcg[i].n_vertex_number = cur as i32;
        prev = cur;
    }
    for endpoint in 0..num_atoms {
        if at[endpoint].endpoint == 0 {
            continue;
        }
        let fictpoint = at[endpoint].endpoint as usize + nv - 1;
        if fictpoint as i32 >= bns.max_vertices
            || num_edges >= bns.max_edges
            || bns.vert[fictpoint].num_adj_edges >= bns.vert[fictpoint].max_adj_edges
            || bns.vert[endpoint].num_adj_edges >= bns.vert[endpoint].max_adj_edges
        {
            ret = BNS_VERT_EDGE_OVFL;
            break;
        }
        bns.vert[endpoint].type_ |= BNS_VERT_TYPE_ENDPOINT;
        let ne = num_edges as usize;
        let st = bns.vert[endpoint].st_edge;
        let e = &mut bns.edge[ne];
        e.cap = (st.cap - st.flow).clamp(0, MAX_TGROUP_EDGE_CAP);
        e.flow = 0;
        e.pass = 0;
        // RESET_EDGE_FORBIDDEN_MASK is 0 in ichirvrs.h: the mask is not applied
        ret = connect_two_vertices(bns, endpoint, fictpoint, ne, false);
        if is_bns_error(ret) {
            break;
        }
        num_edges += 1;
        let e = &mut bns.edge[ne];
        e.cap0 = e.cap;
        e.flow0 = e.flow;
        p_va[endpoint].n_taut_group_edge = num_edges;
    }
    bns.num_edges = num_edges;
    bns.num_vertices += max_tg_number;
    bns.num_t_groups = num_tg as i32;
    bns.tot_st_cap += tot_st_cap;
    ret
}

/// `AddRadicalToMetal`
pub fn add_radical_to_metal(
    tot_st_cap: &mut i32,
    srm: &Srm,
    bns: &mut BnStruct,
    tcg: &AllTcGroups,
) -> i32 {
    let ig = [
        tcg.n_group[TCG_ME_FLOWER0 as usize],
        tcg.n_group[TCG_ME_FLOWER1 as usize],
        tcg.n_group[TCG_ME_FLOWER2 as usize],
        tcg.n_group[TCG_ME_FLOWER3 as usize],
    ];
    let n = ig.iter().filter(|&&g| g >= 0).count();
    if tcg.num_metal_atoms != 0 && srm.b_metal_add_flower != 0 && *tot_st_cap % 2 != 0 && n == 4 {
        let v3 = tcg.p_tcg[ig[3] as usize].n_vertex_number as usize;
        bns.vert[v3].st_edge.cap += 1;
        bns.vert[v3].st_edge.cap0 += 1;
        *tot_st_cap += 1;
        return 1;
    }
    0
}

/// `SetEdgeCapFlow`
pub fn set_edge_cap_flow(e: &mut BnsEdge, cap: i32, flow: i32) {
    e.cap = cap;
    e.cap0 = cap;
    e.flow = flow;
    e.flow0 = flow;
}

/// `AddStCapFlow`
pub fn add_st_cap_flow(
    v: &mut BnsVertex,
    tot_st_flow: &mut i32,
    tot_st_cap: &mut i32,
    cap: i32,
    flow: i32,
) {
    v.st_edge.flow += flow;
    *tot_st_flow += flow;
    v.st_edge.cap += cap;
    *tot_st_cap += cap;
    v.st_edge.flow0 = v.st_edge.flow;
    v.st_edge.cap0 = v.st_edge.cap;
}

/// `SetStCapFlow`
pub fn set_st_cap_flow(
    v: &mut BnsVertex,
    tot_st_flow: &mut i32,
    tot_st_cap: &mut i32,
    cap: i32,
    flow: i32,
) {
    *tot_st_flow += flow - v.st_edge.flow;
    v.st_edge.flow = flow;
    *tot_st_cap += cap - v.st_edge.cap;
    v.st_edge.cap = cap;
    v.st_edge.flow0 = v.st_edge.flow;
    v.st_edge.cap0 = v.st_edge.cap;
}

/// `AddEdgeFlow`: edge `e01`, source vertex `src`, destination `dst`
pub fn add_edge_flow(
    edge_cap: i32,
    edge_flow: i32,
    bns: &mut BnStruct,
    e01: usize,
    src: usize,
    dst: usize,
    tot_st_cap: &mut i32,
    tot_st_flow: &mut i32,
) -> i32 {
    let e = bns.edge[e01];
    let (d, s) = (bns.vert[dst].st_edge, bns.vert[src].st_edge);
    if e.cap < 0 || edge_cap < 0 || e.cap + edge_cap >= EDGE_FLOW_MASK {
        return BNS_PROGRAM_ERR;
    }
    if d.cap < 0
        || d.cap + edge_cap >= EDGE_FLOW_ST_MASK
        || d.flow < 0
        || d.flow + edge_flow >= EDGE_FLOW_ST_MASK
        || s.cap < 0
        || s.flow < 0
        || s.flow + edge_flow >= EDGE_FLOW_ST_MASK
    {
        return BNS_PROGRAM_ERR;
    }
    let e = &mut bns.edge[e01];
    e.cap += edge_cap;
    e.flow += edge_flow;
    e.cap0 = e.cap;
    e.flow0 = e.flow;
    let v = &mut bns.vert[dst].st_edge;
    v.cap += edge_cap;
    v.cap0 = v.cap;
    *tot_st_cap += edge_cap;
    v.flow += edge_flow;
    v.flow0 = v.flow;
    *tot_st_flow += edge_flow;
    let v = &mut bns.vert[src].st_edge;
    v.flow += edge_flow;
    v.flow0 = v.flow;
    *tot_st_flow += edge_flow;
    0
}

/// `ConnectMetalFlower`
pub fn connect_metal_flower(
    _pcur_num_vertices: &mut i32,
    pcur_num_edges: &mut i32,
    tot_st_cap: &mut i32,
    tot_st_flow: &mut i32,
    srm: &Srm,
    bns: &mut BnStruct,
    tcg: &AllTcGroups,
) -> i32 {
    let ig = [
        tcg.n_group[TCG_ME_FLOWER0 as usize],
        tcg.n_group[TCG_ME_FLOWER1 as usize],
        tcg.n_group[TCG_ME_FLOWER2 as usize],
        tcg.n_group[TCG_ME_FLOWER3 as usize],
    ];
    let n = ig.iter().filter(|&&g| g >= 0).count();
    if n == 0 {
        return 0;
    }
    if n != 4 {
        return RI_ERR_PROGR;
    }
    let vg: Vec<usize> = ig
        .iter()
        .map(|&g| tcg.p_tcg[g as usize].n_vertex_number as usize)
        .collect();
    let (g0, g1, g2, g3) = (vg[0], vg[1], vg[2], vg[3]);
    let mut cur_num_edges = *pcur_num_edges;
    let (mut c, mut f) = (0, 0);
    for i in 0..bns.vert[g0].num_adj_edges as usize {
        let e = &bns.edge[bns.iedge[bns.vert[g0].iedge + i] as usize];
        c += e.cap;
        f += e.flow;
    }
    let t0 = &tcg.p_tcg[ig[0] as usize];
    if !is_bns_vt_m_gr(t0.type_ as u16)
        && (t0.edges_cap != bns.vert[g0].st_edge.cap || t0.edges_flow != bns.vert[g0].st_edge.flow)
    {
        return RI_ERR_PROGR;
    }
    if t0.edges_cap != c || t0.edges_flow != f {
        return RI_ERR_PROGR;
    }
    let ia = cur_num_edges as usize;
    let ib = ia + 1;
    let id = ia + 2;
    let ix = ia + 3;
    let iy = ia + 4;
    cur_num_edges += 5;
    for (a, b, e) in [
        (g0, g1, ib),
        (g0, g2, ia),
        (g1, g2, id),
        (g1, g3, iy),
        (g2, g3, ix),
    ] {
        let ret = connect_two_vertices(bns, a, b, e, true);
        if is_bns_error(ret) {
            return ret;
        }
    }
    let dc = c % 2;
    c /= 2;
    let df = f % 2;
    f /= 2;
    let d = srm.n_metal_flower_param_d;
    let c0 = 2 * c + 2 * d + dc;
    let c1 = c + 2 * d + dc - df;
    let c2 = c + 2 * d;
    let c3 = 0;
    let (ca, fa) = (c + 2 * d, c + d - f);
    let (cb, fb) = (c + 2 * d + dc, c + d + dc - (f + df));
    let (cd, fd) = (c + 2 * d, f + d);
    let (cx, fx, cy, fy) = (d, 0, d, 0);
    for v in [c0, c1, c2, c3] {
        if v >= EDGE_FLOW_ST_MASK {
            return BNS_PROGRAM_ERR;
        }
    }
    for (g, cc) in [(g0, c0), (g1, c1), (g2, c2), (g3, c3)] {
        set_st_cap_flow(&mut bns.vert[g], tot_st_flow, tot_st_cap, cc, cc);
    }
    set_edge_cap_flow(&mut bns.edge[ia], ca, fa);
    set_edge_cap_flow(&mut bns.edge[ib], cb, fb);
    set_edge_cap_flow(&mut bns.edge[id], cd, fd);
    set_edge_cap_flow(&mut bns.edge[ix], cx, fx);
    set_edge_cap_flow(&mut bns.edge[iy], cy, fy);
    *pcur_num_edges = cur_num_edges;
    0
}

/// `ConnectSuperCGroup`; returns the number of connected groups or an error
#[allow(clippy::too_many_arguments)]
pub fn connect_super_c_group(
    n_super_c_group: i32,
    n_add_groups: &[i32],
    pcur_num_vertices: &mut i32,
    pcur_num_edges: &mut i32,
    tot_st_cap: &mut i32,
    tot_st_flow: &mut i32,
    bns: &mut BnStruct,
    tcg: &mut AllTcGroups,
) -> i32 {
    let (i1, fst) = if n_super_c_group >= 0 {
        let i1 = tcg.n_group[n_super_c_group as usize];
        if i1 < 0 {
            return 0;
        }
        (i1, 0usize)
    } else {
        (-1, 1usize)
    };
    let num_groups = n_add_groups
        .iter()
        .filter(|&&g| {
            let x = tcg.n_group[g as usize];
            x >= 0 && x != i1
        })
        .count();
    if num_groups < 1 {
        return 0;
    }
    let mut cur_num_vertices = *pcur_num_vertices;
    let mut cur_num_edges = *pcur_num_edges;
    let mut ix = vec![0i32; num_groups + 1];
    let mut jx = vec![0usize; num_groups + 1];
    let mut e0x = vec![0usize; num_groups + 1];

    let j0 = cur_num_vertices as usize;
    let (piedge, pmax) = (bns.vert[j0 - 1].iedge, bns.vert[j0 - 1].max_adj_edges);
    bns.vert[j0].iedge = piedge + pmax as usize;
    bns.vert[j0].max_adj_edges = (num_groups as i32 + 1 + BNS_ADD_EDGES) as u16;
    bns.vert[j0].num_adj_edges = 0;
    bns.vert[j0].type_ = BNS_VT_YVCONNECTOR;
    cur_num_vertices += 1;

    if fst == 0 {
        ix[0] = i1;
        jx[0] = tcg.p_tcg[i1 as usize].n_vertex_number as usize;
    }
    let mut j = 1;
    for &g in n_add_groups {
        let x = tcg.n_group[g as usize];
        if x >= 0 && x != i1 {
            ix[j] = x;
            jx[j] = tcg.p_tcg[x as usize].n_vertex_number as usize;
            j += 1;
        }
    }
    for i in fst..=num_groups {
        e0x[i] = cur_num_edges as usize;
        let ret = connect_two_vertices(bns, j0, jx[i], cur_num_edges as usize, true);
        if is_bns_error(ret) {
            return ret;
        }
        if i != 0 {
            tcg.p_tcg[ix[i] as usize].n_forward_edge = cur_num_edges;
        } else {
            tcg.p_tcg[ix[i] as usize].n_backward_edge = cur_num_edges;
        }
        cur_num_edges += 1;
    }
    for i in 1..=num_groups {
        let g = tcg.p_tcg[ix[i] as usize];
        let n_delta = g.st_cap - g.edges_cap;
        let edge_cap = g.edges_cap + n_delta;
        let edge_flow = g.edges_cap - g.edges_flow;
        let ret = add_edge_flow(
            edge_cap,
            edge_flow,
            bns,
            e0x[i],
            jx[i],
            j0,
            tot_st_cap,
            tot_st_flow,
        );
        if is_bns_error(ret) {
            return ret;
        }
    }
    if fst == 0 {
        let st = bns.vert[j0].st_edge;
        let (edge_cap, edge_flow) = (st.cap, st.cap - st.flow);
        let ret = add_edge_flow(
            edge_cap,
            edge_flow,
            bns,
            e0x[0],
            j0,
            jx[0],
            tot_st_cap,
            tot_st_flow,
        );
        if is_bns_error(ret) {
            return ret;
        }
        let g = &mut tcg.p_tcg[ix[0] as usize];
        g.edges_cap += edge_cap;
        g.edges_flow += edge_flow;
        g.st_cap += edge_cap;
        g.st_flow += edge_flow;
    } else {
        let st = &mut bns.vert[j0].st_edge;
        *tot_st_cap += st.flow - st.cap;
        st.cap += st.flow - st.cap;
        st.cap0 = st.cap;
    }
    *pcur_num_vertices = cur_num_vertices;
    *pcur_num_edges = cur_num_edges;
    num_groups as i32
}

/// `AddCGroups2TCGBnStruct`: adds the charge-group vertices, the
/// ChargeStruct of every c-point, the metal flower and the (+)/(-) connections
pub fn add_cgroups_to_tcg_bn_struct(
    bns: &mut BnStruct,
    at: &[InpAtom],
    num_atoms: usize,
    srm: &Srm,
    p_va: &mut [ValAt],
    tcg: &mut AllTcGroups,
    max_add_edges: i32,
) -> i32 {
    let num_cg = tcg.num_tc_groups - tcg.num_tgroups;
    let fst_cg_vertex = bns.num_vertices;
    let fst_cg_group = tcg.num_tgroups as usize;
    let num_vertices = bns.num_vertices;
    let num_edges = bns.num_edges;
    if num_cg <= 0 {
        return 0;
    }
    let md = srm.n_metal_max_charge_d;
    if num_vertices >= bns.max_vertices {
        return BNS_VERT_EDGE_OVFL;
    }
    for v in &mut bns.vert[num_vertices as usize..bns.max_vertices as usize] {
        *v = BnsVertex::default();
    }
    let mut tot_st_cap = bns.tot_st_cap;
    let mut tot_st_flow = bns.tot_st_flow;
    let mut prev = fst_cg_vertex as usize - 1;
    for i in 0..num_cg as usize {
        let cur = prev + 1;
        let (piedge, pmax) = (bns.vert[prev].iedge, bns.vert[prev].max_adj_edges);
        let g = tcg.p_tcg[i + fst_cg_group];
        let v = &mut bns.vert[cur];
        v.iedge = piedge + pmax as usize;
        v.max_adj_edges = (g.num_edges + max_add_edges) as u16;
        v.num_adj_edges = 0;
        v.st_edge.flow += g.st_flow;
        tot_st_flow += g.st_flow;
        v.st_edge.cap += g.st_cap;
        tot_st_cap += g.st_cap;
        v.st_edge.flow0 = v.st_edge.flow;
        v.st_edge.cap0 = v.st_edge.cap;
        v.type_ = g.type_ as u16;
        tcg.p_tcg[i + fst_cg_group].n_vertex_number = cur as i32;
        prev = cur;
    }
    let mut cur_num_vertices = prev as i32 + 1;
    let mut cur_num_edges = num_edges;

    for c_point in 0..num_atoms {
        let k = p_va[c_point].c_cn_list_index as usize;
        if k == 0 {
            continue;
        }
        let cl = &CN_LIST[k - 1];
        let (p_cn, cn_len, b_metal) = (cl.p_cn, cl.len as usize, cl.bits == CN_BITS_ME);
        let base = prev;
        let sc = |x: i8| -> i32 {
            if !b_metal {
                i32::from(x)
            } else if x != 0 {
                md
            } else {
                0
            }
        };
        for i1 in 0..cn_len {
            if !is_bns_vt_chrg_struct(p_cn[i1].v.type_) {
                continue;
            }
            let cur = base + i1;
            let (piedge, pmax) = (bns.vert[prev].iedge, bns.vert[prev].max_adj_edges);
            let v = &mut bns.vert[cur];
            v.iedge = piedge + pmax as usize;
            v.max_adj_edges = p_cn[i1].v.valence as u16;
            v.num_adj_edges = 0;
            add_st_cap_flow(
                v,
                &mut tot_st_flow,
                &mut tot_st_cap,
                sc(p_cn[i1].v.cap),
                sc(p_cn[i1].v.flow),
            );
            v.type_ = p_cn[i1].v.type_;
            prev = cur;
            cur_num_vertices = cur as i32 + 1;
            if (bns.vert[cur].iedge + bns.vert[cur].max_adj_edges as usize) as i32 >= bns.max_iedges
            {
                return BNS_VERT_EDGE_OVFL;
            }
            if cur_num_vertices >= bns.max_vertices {
                return BNS_VERT_EDGE_OVFL;
            }
        }
        for i1 in 0..cn_len {
            // find the vertex of the ChargeStruct node `i`
            let find = |bns: &mut BnStruct,
                        tcg: &AllTcGroups,
                        i: usize,
                        tot_st_flow: &mut i32,
                        tot_st_cap: &mut i32|
             -> Result<Option<usize>, i32> {
                let t = p_cn[i].v.type_;
                if t & BNS_VERT_TYPE_ATOM != 0 {
                    add_st_cap_flow(
                        &mut bns.vert[c_point],
                        tot_st_flow,
                        tot_st_cap,
                        sc(p_cn[i].v.cap),
                        sc(p_cn[i].v.flow),
                    );
                    Ok(Some(c_point))
                } else if is_bns_vt_c_gr(t) {
                    let mut found = None;
                    for j in 0..num_cg as usize {
                        if t == bns.vert[fst_cg_vertex as usize + j].type_ {
                            found = Some(j);
                            break;
                        }
                    }
                    if let Some(j) = found {
                        let g = &tcg.p_tcg[j + fst_cg_group];
                        if g.type_ != i32::from(t) || g.ord_num != 0 {
                            return Err(RI_ERR_PROGR);
                        }
                        Ok(Some(fst_cg_vertex as usize + j))
                    } else {
                        Ok(None)
                    }
                } else if is_bns_vt_m_gr(t) {
                    let k = tcg.n_group[TCG_ME_FLOWER0 as usize];
                    if k < 0
                        || tcg.p_tcg[k as usize].type_ != i32::from(t)
                        || tcg.p_tcg[k as usize].ord_num != 0
                        || srm.b_metal_add_flower == 0
                    {
                        return Err(RI_ERR_PROGR);
                    }
                    Ok(Some(tcg.p_tcg[k as usize].n_vertex_number as usize))
                } else if is_bns_vt_chrg_struct(t) {
                    Ok(Some(base + i))
                } else {
                    Ok(None)
                }
            };
            let pv1 = match find(bns, tcg, i1, &mut tot_st_flow, &mut tot_st_cap) {
                Ok(Some(v)) => v,
                Ok(None) => return BNS_BOND_ERR,
                Err(e) => return e,
            };
            let mut k = 0;
            while k < MAX_CN_VAL && p_cn[i1].e[k].neigh != 0 {
                let i2 = p_cn[i1].e[k].neigh as usize - 1;
                let pv2 = match find(bns, tcg, i2, &mut tot_st_flow, &mut tot_st_cap) {
                    Ok(Some(v)) => v,
                    Ok(None) => return BNS_BOND_ERR,
                    Err(e) => return e,
                };
                if pv1 == pv2 {
                    return BNS_BOND_ERR;
                }
                let ne = cur_num_edges as usize;
                let (t1, t2) = (p_cn[i1].v.type_, p_cn[i2].v.type_);
                if (is_bns_vt_m_gr(t1) && t2 & BNS_VERT_TYPE_ATOM != 0)
                    || (is_bns_vt_m_gr(t2) && t1 & BNS_VERT_TYPE_ATOM != 0)
                {
                    let (needs, _, _, ecap, eflow) = atom_stcap_stflow(at, p_va, srm, c_point);
                    bns.edge[ne].cap = ecap;
                    bns.edge[ne].flow = eflow;
                    if needs == 0 {
                        return RI_ERR_PROGR;
                    }
                    p_va[c_point].n_metal_group_edge = cur_num_edges + 1;
                    let ef = bns.edge[ne].flow;
                    let free = i32::from(p_va[c_point].c_init_free_valences);
                    let st = &mut bns.vert[c_point].st_edge;
                    st.flow += ef;
                    st.cap += ef + free;
                    st.flow0 = st.flow;
                    st.cap0 = st.cap;
                    tot_st_flow += ef;
                    tot_st_cap += ef + free;
                } else {
                    bns.edge[ne].cap = sc(p_cn[i1].e[k].cap);
                    bns.edge[ne].flow = sc(p_cn[i1].e[k].flow);
                }
                bns.edge[ne].forbidden = if p_cn[i1].e[k].b_forbidden_edge != 0 {
                    BNS_EDGE_FORBIDDEN_MASK
                } else {
                    0
                };
                bns.edge[ne].pass = 0;
                if bns.vert[pv1].num_adj_edges >= bns.vert[pv1].max_adj_edges
                    || bns.vert[pv2].num_adj_edges >= bns.vert[pv2].max_adj_edges
                    || cur_num_edges >= bns.max_edges
                {
                    return BNS_VERT_EDGE_OVFL;
                }
                let ret = connect_two_vertices(bns, pv1, pv2, ne, false);
                if is_bns_error(ret) {
                    return ret;
                }
                bns.edge[ne].cap0 = bns.edge[ne].cap;
                bns.edge[ne].flow0 = bns.edge[ne].flow;
                let (y1, y2) = (bns.vert[pv1].type_, bns.vert[pv2].type_);
                let ty = if is_bns_vt_c_gr(y1) {
                    y1
                } else if is_bns_vt_c_gr(y2) {
                    y2
                } else {
                    0
                };
                if ty != 0 {
                    if ty & BNS_VERT_TYPE_C_NEGATIVE != 0 {
                        p_va[c_point].n_c_minus_group_edge = cur_num_edges + 1;
                    } else {
                        p_va[c_point].n_c_plus_group_edge = cur_num_edges + 1;
                    }
                }
                cur_num_edges += 1;
                k += 1;
            }
        }
    }

    // pass 2
    for c_point in 0..num_atoms {
        let st_cap = bns.vert[c_point].st_edge.cap;
        for k in 0..bns.vert[c_point].num_adj_edges as usize {
            let ie = bns.iedge[bns.vert[c_point].iedge + k] as usize;
            let c_neigh = (bns.edge[ie].neighbor12 as usize) ^ c_point;
            if c_neigh > c_point || bns.vert[c_neigh].type_ & BNS_VERT_TYPE_ATOM == 0 {
                continue;
            }
            let st_cap2 = bns.vert[c_neigh].st_edge.cap;
            let mut max_edge_flow = st_cap.min(st_cap2);
            if (srm.b_metal_add_flower != 0
                && srm.n_metal_min_bond_order == 0
                && p_va[c_point].c_metal != 0
                && p_va[c_point].c_num_bonds_to_metal != 0)
                || (p_va[c_neigh].c_metal != 0 && p_va[c_neigh].c_num_bonds_to_metal != 0)
            {
                max_edge_flow = max_edge_flow.min(MAX_BOND_EDGE_CAP + 1);
            } else {
                max_edge_flow = max_edge_flow.min(MAX_BOND_EDGE_CAP);
            }
            if i32::from(at[c_point].bond_type[k]) == BOND_TYPE_SINGLE {
                bns.edge[ie].cap = max_edge_flow;
                bns.edge[ie].cap0 = max_edge_flow;
            }
        }
    }

    let mut ret = connect_metal_flower(
        &mut cur_num_vertices,
        &mut cur_num_edges,
        &mut tot_st_cap,
        &mut tot_st_flow,
        srm,
        bns,
        tcg,
    );
    if ret < 0 {
        return ret;
    }
    connect_super_c_group(
        TCG_PLUS,
        &[TCG_PLUS0, TCG_PLUS_C0, TCG_PLUS_M0],
        &mut cur_num_vertices,
        &mut cur_num_edges,
        &mut tot_st_cap,
        &mut tot_st_flow,
        bns,
        tcg,
    );
    connect_super_c_group(
        TCG_MINUS,
        &[TCG_MINUS0, TCG_MINUS_C0, TCG_MINUS_M0],
        &mut cur_num_vertices,
        &mut cur_num_edges,
        &mut tot_st_cap,
        &mut tot_st_flow,
        bns,
        tcg,
    );
    let ret3 = connect_super_c_group(
        -1,
        &[TCG_PLUS, TCG_MINUS],
        &mut cur_num_vertices,
        &mut cur_num_edges,
        &mut tot_st_cap,
        &mut tot_st_flow,
        bns,
        tcg,
    );
    let cg_charge = tcg.total_charge - tcg.tgroup_charge - tcg.charge_on_atoms;
    ret = 1;
    if ret3 > 0 {
        let vpm = cur_num_vertices as usize - 1;
        let n = i32::from(tcg.n_group[TCG_PLUS as usize] >= 0);
        let k = i32::from(tcg.n_group[TCG_MINUS as usize] >= 0);
        let (mut e_plus, mut e_minus): (Option<usize>, Option<usize>) = (None, None);
        let adj = bns.vert[vpm].num_adj_edges;
        let i0 = bns.vert[vpm].iedge;
        if adj == 2 && k + n == 2 {
            e_plus = Some(bns.iedge[i0] as usize);
            e_minus = Some(bns.iedge[i0 + 1] as usize);
        } else if adj == 1 && k + n == 1 {
            if tcg.n_group[TCG_PLUS as usize] >= 0 {
                e_plus = Some(bns.iedge[i0] as usize);
            } else if tcg.n_group[TCG_MINUS as usize] >= 0 {
                e_minus = Some(bns.iedge[i0] as usize);
            }
        } else if k + n != 0 {
            return BNS_BOND_ERR;
        }
        let v_plus = e_plus.map(|e| bns.edge[e].neighbor12 as usize ^ vpm);
        let v_minus = e_minus.map(|e| bns.edge[e].neighbor12 as usize ^ vpm);
        let p_vert = v_plus.or(v_minus);
        let p_edge = e_plus.or(e_minus);
        if let Some(e) = e_minus {
            tcg.n_edge_minus = e as i32;
        }
        if let Some(e) = e_plus {
            tcg.n_edge_plus = e as i32;
        }
        if let Some(e) = p_edge {
            tcg.n_edge4charge = e as i32;
        }
        if let (Some(v), Some(e)) = (p_vert, p_edge) {
            if cg_charge > 0 {
                let st = &mut bns.vert[vpm].st_edge;
                st.cap += cg_charge;
                tot_st_cap += cg_charge;
                st.cap0 = st.cap;
            }
            if cg_charge < 0 {
                let st = &mut bns.vert[v].st_edge;
                st.cap -= cg_charge;
                tot_st_cap -= cg_charge;
                st.cap0 = st.cap;
                let ed = &mut bns.edge[e];
                if ed.cap - ed.flow + cg_charge < 0 {
                    ed.cap = ed.flow - cg_charge;
                }
            }
            tcg.added_charge = cg_charge;
        }
        if cg_charge == 0 || (p_vert.is_some() && p_edge.is_some()) {
            ret = 2;
        }
    }
    add_radical_to_metal(&mut tot_st_cap, srm, bns, tcg);
    bns.num_edges = cur_num_edges;
    bns.num_vertices = cur_num_vertices;
    bns.num_c_groups = num_cg;
    bns.tot_st_cap = tot_st_cap;
    bns.tot_st_flow = tot_st_flow;
    ret
}

/// `nNumEdgesToCnVertex`
pub fn n_num_edges_to_cn_vertex(p_cn: &[CNode], len: usize, v: usize) -> i32 {
    let v1 = v as i32 + 1;
    let mut n = 0;
    for (i, node) in p_cn.iter().enumerate().take(len) {
        let mut j = 0;
        while j < MAX_CN_VAL && node.e[j].neigh != 0 {
            n += i32::from(i == v || i32::from(node.e[j].neigh) == v1);
            j += 1;
        }
    }
    n
}

/// `AllocateAndInitTCGBnStruct`: the BNS with the atom vertices and the bond
/// edges, sized for the t- and c-groups; `None` when the sizes do not fit.
/// Returns also the number of bonds made single.
pub fn allocate_and_init_tcg_bn_struct(
    at: &mut [InpAtom],
    num_atoms: usize,
    srm: &Srm,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    n_max_add_atoms: i32,
    n_max_add_edges: i32,
    max_altp: i32,
) -> Option<(BnStruct, i32)> {
    let max_vertices = tcg.n_vertices + n_max_add_atoms;
    let num_bonds = tcg.num_bonds;
    let num_edges = tcg.n_edges;
    let max_edges = num_edges + (n_max_add_edges + NUM_KINDS_OF_GROUPS) * max_vertices;
    let max_iedges = 2 * max_edges + tcg.n_add_iedges;
    let mut len_alt_path = max_vertices + I_ALTP_HDR_LEN as i32 + 1;
    len_alt_path += (max_vertices / 2).max(16);
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

    let mut iedge_pos = 0usize;
    let mut tot_st_cap = 0;
    let mut tot_st_flow = 0;
    for i in 0..num_atoms {
        let mut k = i32::from(at[i].valence) + i32::from(at[i].endpoint != 0) + n_max_add_edges;
        if p_va[i].c_cn_list_index > 0 {
            let cl = &CN_LIST[p_va[i].c_cn_list_index as usize - 1];
            k += n_num_edges_to_cn_vertex(cl.p_cn, cl.len as usize, 0);
        }
        p.vert[i].max_adj_edges = k as u16;
        p.vert[i].iedge = iedge_pos;
        iedge_pos += k as usize;
        let (needs, c1, _f1, _, _) = atom_stcap_stflow(at, p_va, srm, i);
        let st_cap = c1
            + if needs != 0 {
                0
            } else {
                i32::from(p_va[i].c_init_free_valences)
            };
        p.vert[i].st_edge.cap = st_cap;
        p.vert[i].st_edge.cap0 = st_cap;
        tot_st_cap += st_cap;
    }
    let num_iedges = iedge_pos as i32;
    if max_iedges - num_iedges < (n_max_add_edges + NUM_KINDS_OF_GROUPS) * max_vertices {
        return None;
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

    let mut num_changed_bonds = 0;
    let mut n_edges = 0usize;
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
                bond_type = BOND_TYPE_SINGLE;
                at[i].bond_type[j] = (bond_mark | bond_type) as u8;
                num_changed_bonds += 1;
            }
            let edge_flow;
            if neigh > i {
                let (n1, c1, _, _, _) = atom_stcap_stflow(at, p_va, srm, i);
                let _c1 = c1
                    + if n1 != 0 {
                        0
                    } else {
                        i32::from(p_va[i].c_init_free_valences)
                    };
                let (n2, c2, _, _, _) = atom_stcap_stflow(at, p_va, srm, neigh);
                let _c2 = c2
                    + if n2 != 0 {
                        0
                    } else {
                        i32::from(p_va[neigh].c_init_free_valences)
                    };
                let (ef, edge_cap, _min_order, _) = bond_flow_maxcap_minorder(at, p_va, srm, i, j);
                edge_flow = ef;
                let e = &mut p.edge[n_edges];
                e.neighbor1 = i as u16;
                e.neighbor12 = (i ^ neigh) as u16;
                e.flow = ef;
                e.flow0 = ef;
                e.cap = edge_cap;
                e.cap0 = edge_cap;
                e.neigh_ord = [j as u16, k as u16];
                e.pass = 0;
                e.forbidden = 0;
                if bond_type == BOND_TYPE_DOUBLE {
                    let mut m = 0;
                    while m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m] != 0 {
                        if at[i].sb_ord[m] as usize == j {
                            e.forbidden |= BNS_EDGE_FORBIDDEN_MASK;
                            break;
                        }
                        m += 1;
                    }
                }
                let vi = p.vert[i].iedge + j;
                p.iedge[vi] = n_edges as i32;
                let vn = p.vert[neigh].iedge + k;
                p.iedge[vn] = n_edges as i32;
                n_edges += 1;
            } else {
                let iedge2 = p.iedge[p.vert[neigh].iedge + k] as usize;
                edge_flow = p.edge[iedge2].flow;
            }
            st_flow += edge_flow;
        }
        let v = &mut p.vert[i];
        v.num_adj_edges = valence as u16;
        v.st_edge.flow = st_flow;
        v.st_edge.flow0 = st_flow;
        v.type_ = BNS_VERT_TYPE_ATOM;
        tot_st_flow += st_flow;
    }
    p.num_edges = n_edges as i32;
    p.num_iedges = num_iedges;
    p.num_added_edges = 0;
    p.tot_st_cap = tot_st_cap;
    p.tot_st_flow = tot_st_flow;
    Some((p, num_changed_bonds / 2))
}

// ---------------------------------------------------------------------------
// Charge changes of the BNS paths, BNS runs, forbidden edge masks

/// `IncrZeroBondsAndClearEndpts`
pub fn incr_zero_bonds_and_clear_endpts(at: &mut [InpAtom], num_at: usize, i_component: i32) {
    for a in &mut at[..num_at] {
        a.endpoint = 0;
        a.component = i_component as _;
        for j in 0..a.valence as usize {
            if a.bond_type[j] == 0 {
                a.bond_type[j] = BOND_TYPE_SINGLE as u8;
                a.chem_bonds_valence += BOND_TYPE_SINGLE as i8;
            }
        }
    }
}

/// `IncrZeroBonds`
pub fn incr_zero_bonds(at: &mut [InpAtom], num_at: usize, i_component: i32) {
    for a in &mut at[..num_at] {
        a.component = i_component as _;
        for j in 0..a.valence as usize {
            if a.bond_type[j] == 0 {
                a.bond_type[j] = BOND_TYPE_SINGLE as u8;
                a.chem_bonds_valence += BOND_TYPE_SINGLE as i8;
            }
        }
    }
}

/// `ClearEndpts`
pub fn clear_endpts(at: &mut [InpAtom], num_at: usize) {
    for a in &mut at[..num_at] {
        a.endpoint = 0;
    }
}

fn any_vert_type(x: u16) -> bool {
    x & (BNS_VERT_TYPE_ATOM | BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP) != 0
        && x & BNS_VERT_TYPE_SUPER_TGROUP == 0
}

const VF_USED_IN: i32 = 1;
const VF_USED_OUT: i32 = 2;
const VF_USED_ALL: i32 = VF_USED_IN | VF_USED_OUT;

/// `VF`: a vertex on an alternating path with its in/out edges
#[derive(Debug, Clone, Copy)]
pub struct Vf {
    pub type_: u16,
    pub v: i32,
    pub e_in: i32,
    pub e_out: i32,
    pub delta_in: i32,
    pub delta_out: i32,
    pub b_used: i32,
}

impl Default for Vf {
    fn default() -> Self {
        Vf {
            type_: 0,
            v: NO_VERTEX,
            e_in: NO_VERTEX,
            e_out: NO_VERTEX,
            delta_in: 0,
            delta_out: 0,
            b_used: 0,
        }
    }
}

/// `GetDeltaChargeFromVF`
pub fn get_delta_charge_from_vf(bns: &BnStruct, p_va: &[ValAt], vf: &mut Vf) -> i32 {
    let ie_in1 = if vf.b_used & VF_USED_IN == 0 && vf.e_in >= 0 && vf.delta_in != 0 {
        vf.e_in + 1
    } else {
        NO_VERTEX
    };
    let ie_out1 = if vf.b_used & VF_USED_OUT == 0 && vf.e_out >= 0 && vf.delta_out != 0 {
        vf.e_out + 1
    } else {
        NO_VERTEX
    };
    if vf.type_ & BNS_VERT_TYPE_C_GROUP == 0
        || vf.type_ & BNS_VERT_TYPE_SUPER_TGROUP != 0
        || (ie_in1 == NO_VERTEX && ie_out1 == NO_VERTEX)
    {
        return 0;
    }
    let mut v = NO_VERTEX;
    for i in 0..bns.num_atoms as usize {
        let hit = if vf.type_ & BNS_VERT_TYPE_C_NEGATIVE != 0 {
            p_va[i].n_c_minus_group_edge
        } else {
            p_va[i].n_c_plus_group_edge
        };
        if hit == ie_in1 || hit == ie_out1 {
            v = i as i32;
            break;
        }
    }
    if v == NO_VERTEX {
        return 0;
    }
    let va = &p_va[v as usize];
    let mut n_init_charge = i32::from(va.c_init_charge);
    let (mut plus_flow, mut minus_flow) = (0, 0);
    let e_cplus = va.n_c_plus_group_edge - 1;
    if e_cplus >= 0 {
        let e = &bns.edge[e_cplus as usize];
        plus_flow = e.cap - e.flow;
    }
    let e_cminus = va.n_c_minus_group_edge - 1;
    if e_cminus >= 0 {
        minus_flow = -bns.edge[e_cminus as usize].flow;
    }
    n_init_charge += plus_flow + minus_flow;
    let mut n_delta_charge = 0;
    if vf.b_used & VF_USED_OUT == 0 && (vf.e_out == e_cplus || vf.e_out == e_cminus) {
        n_delta_charge -= vf.delta_out;
        vf.b_used |= VF_USED_OUT;
    }
    if vf.b_used & VF_USED_IN == 0 && (vf.e_in == e_cplus || vf.e_in == e_cminus) {
        n_delta_charge -= vf.delta_in;
        vf.b_used |= VF_USED_IN;
    }
    if n_init_charge == 0 && n_delta_charge != 0 {
        1
    } else if n_init_charge != 0 && n_init_charge + n_delta_charge == 0 {
        -1
    } else {
        0
    }
}

/// `EvaluateChargeChanges`: `(ret, delta_h, delta_charge, num_visited_atoms)`
pub fn evaluate_charge_changes(bns: &mut BnStruct, p_va: &[ValAt]) -> (i32, i32, i32, i32) {
    let (mut pn_delta_h, mut pn_delta_charge, mut pn_visited) = (0, 0, 0);
    let mut err = 0;
    let mut pass = bns.num_altp - 1;
    while pass >= 0 {
        bns.alt_path = Some(pass as usize);
        let ap = &bns.altp[pass as usize];
        let mut v1 = ap[I_ALTP_START_ATOM].number();
        let n = ap[I_ALTP_PATH_LEN].number();
        let mut delta = ap[I_ALTP_FLOW].flow[0];
        let v_last = ap[I_ALTP_END_ATOM].number();
        let mut v2 = NO_VERTEX;
        let mut vf = [Vf::default(); 4];
        let mut i_prev = 0;
        if any_vert_type(bns.vert[v1 as usize].type_) {
            vf[2].type_ = bns.vert[v1 as usize].type_;
            vf[2].v = v1;
            i_prev = 2;
        }
        let mut n_num_delta_h = 0;
        let mut n_num_delta_charge = 0;
        // (C counts the first atom and resets the counter right away)
        let mut n_num_visited = 0;
        for i in 0..n {
            let ineigh1 = bns.altp[pass as usize][I_ALTP_NEIGHBOR + i as usize].ineigh(0) as usize;
            let ie = bns.iedge[bns.vert[v1 as usize].iedge + ineigh1];
            let nb12 = i32::from(bns.edge[ie as usize].neighbor12);
            if i_prev != 0 {
                vf[2].e_out = ie;
                vf[2].delta_out = delta;
            }
            v2 = nb12 ^ v1;
            let t2 = bns.vert[v2 as usize].type_;
            if t2 & BNS_VERT_TYPE_ATOM != 0 {
                n_num_visited += 1;
            }
            if (any_vert_type(t2) || i == n - 1)
                && vf[0].type_ & BNS_VERT_TYPE_C_GROUP != 0
                && vf[0].b_used != VF_USED_ALL
            {
                n_num_delta_charge += get_delta_charge_from_vf(bns, p_va, &mut vf[0]);
            }
            if any_vert_type(t2) {
                vf[0] = vf[1];
                vf[1] = vf[2];
                vf[2] = vf[3];
                vf[2].v = v2;
                vf[2].type_ = t2;
                vf[2].e_in = ie;
                vf[2].delta_in = delta;
                i_prev = 2;
            } else if i == n - 1 {
                vf[0] = vf[1];
                vf[1] = vf[2];
                vf[2] = vf[3];
                i_prev = 1;
            } else {
                i_prev = 0;
            }
            if i_prev != 0 && vf[1].type_ & BNS_VERT_TYPE_ATOM != 0 {
                let v = vf[1].v as usize;
                let mut n_init_charge = i32::from(p_va[v].c_init_charge);
                let (mut plus_flow, mut minus_flow) = (0, 0);
                let e_cplus = p_va[v].n_c_plus_group_edge - 1;
                if e_cplus >= 0 {
                    let e = &bns.edge[e_cplus as usize];
                    plus_flow = e.cap - e.flow;
                }
                let e_cminus = p_va[v].n_c_minus_group_edge - 1;
                if e_cminus >= 0 {
                    minus_flow = -bns.edge[e_cminus as usize].flow;
                }
                n_init_charge += plus_flow + minus_flow;
                let (mut n_delta_h, mut n_delta_charge) = (0, 0);
                if vf[0].type_ & BNS_VERT_TYPE_TGROUP != 0 {
                    n_delta_h -= delta;
                } else if vf[0].type_ & BNS_VERT_TYPE_C_GROUP != 0
                    && vf[0].b_used & VF_USED_OUT == 0
                    && (vf[0].e_out == e_cplus || vf[0].e_out == e_cminus)
                {
                    n_delta_charge -= vf[0].delta_out;
                    vf[0].b_used |= VF_USED_OUT;
                }
                if vf[2].type_ & BNS_VERT_TYPE_TGROUP != 0 {
                    n_delta_h += delta;
                } else if vf[2].type_ & BNS_VERT_TYPE_C_GROUP != 0
                    && vf[2].b_used & VF_USED_IN == 0
                    && (vf[2].e_in == e_cplus || vf[2].e_in == e_cminus)
                {
                    n_delta_charge -= vf[2].delta_in;
                    vf[2].b_used |= VF_USED_IN;
                }
                if n_init_charge == 0 && n_delta_charge != 0 {
                    n_num_delta_charge += 1;
                } else if n_init_charge != 0 && n_init_charge + n_delta_charge == 0 {
                    n_num_delta_charge -= 1;
                }
                n_num_delta_h += n_delta_h.abs();
                vf[1].b_used = VF_USED_ALL;
            }
            delta = -delta;
            v1 = v2;
        }
        for j in 0..3 {
            n_num_delta_charge += get_delta_charge_from_vf(bns, p_va, &mut vf[j]);
        }
        pn_delta_h += n_num_delta_h;
        pn_delta_charge += n_num_delta_charge;
        pn_visited += n_num_visited;
        if v2 != v_last {
            err = BNS_PROGRAM_ERR;
        }
        pass -= 1;
    }
    (err, pn_delta_h, pn_delta_charge, pn_visited)
}

/// result of [`run_bns_test_once`]
#[derive(Debug, Clone, Copy, Default)]
pub struct BnsTestResult {
    pub ret: i32,
    pub v_first: i32,
    pub v_last: i32,
    pub path_len: i32,
    pub delta_h: i32,
    pub delta_charge: i32,
    pub num_visited_atoms: i32,
}

/// `RunBnsTestOnce`: the `delta_h`, `delta_charge` and `num_visited_atoms`
/// of `res` keep their values when no path is found
pub fn run_bns_test_once(
    bns: &mut BnStruct,
    bd: &mut BnData,
    p_va: &[ValAt],
    res: &mut BnsTestResult,
) -> i32 {
    re_init_bn_struct_alt_paths(bns);
    bns.alt_path = Some(0);
    bns.num_altp = 0;
    bns.b_change_flow = 0;
    let delta = balanced_network_search(bns, bd, 0);
    if delta > 0 {
        bns.alt_path = Some(0);
        let ap = &bns.altp[0];
        res.v_first = ap[I_ALTP_START_ATOM].number();
        res.path_len = ap[I_ALTP_PATH_LEN].number();
        res.v_last = ap[I_ALTP_END_ATOM].number();
        bns.num_altp += 1;
        let (_, dh, dc, nv) = evaluate_charge_changes(bns, p_va);
        res.delta_h = dh;
        res.delta_charge = dc;
        res.num_visited_atoms = nv;
    } else {
        res.v_first = NO_VERTEX;
        res.path_len = 0;
        res.v_last = NO_VERTEX;
    }
    re_init_bn_struct_alt_paths(bns);
    let ret = re_init_bn_data(bd);
    res.ret = delta;
    if delta >= 0 && ret > 0 { -ret } else { delta }
}

/// `RunBnsRestoreOnce`
pub fn run_bns_restore_once(bns: &mut BnStruct, bd: &mut BnData) -> i32 {
    let mut total = 0;
    re_init_bn_struct_alt_paths(bns);
    loop {
        let delta = run_balanced_network_search(bns, bd, BNS_EF_CHNG_FLOW);
        if is_bns_error(delta) {
            return delta;
        }
        total += delta;
        re_init_bn_struct_alt_paths(bns);
        let ret = re_init_bn_data(bd);
        if ret > 0 {
            return -ret;
        }
        if !(delta > 0 && ret == 0) {
            break;
        }
    }
    bns.tot_st_flow += 2 * total;
    total
}

/// `RemoveForbiddenEdgeMask`
pub fn remove_forbidden_edge_mask(bns: &mut BnStruct, edges: &EdgeList, mask: i32) {
    for &ie in &edges.pn_edges[..edges.num_edges as usize] {
        bns.edge[ie as usize].forbidden &= !mask as i8;
    }
}

/// `SetForbiddenEdgeMask`
pub fn set_forbidden_edge_mask(bns: &mut BnStruct, edges: &EdgeList, mask: i32) {
    for &ie in &edges.pn_edges[..edges.num_edges as usize] {
        bns.edge[ie as usize].forbidden |= mask as i8;
    }
}

/// `RemoveForbiddenBondFlowBits`
pub fn remove_forbidden_bond_flow_bits(bns: &mut BnStruct, mask: i32) {
    let n = bns.num_bonds as usize;
    for e in &mut bns.edge[..n] {
        e.forbidden &= !mask as i8;
    }
}

/// `GetChargeFlowerUpperEdge`: the index of the upper edge of the charge
/// flower of the charge edge, or `NO_VERTEX`
pub fn get_charge_flower_upper_edge(bns: &BnStruct, charge_edge: i32) -> i32 {
    if charge_edge < 0 {
        return NO_VERTEX;
    }
    let pe = &bns.edge[charge_edge as usize];
    let mut vc = i32::from(pe.neighbor1);
    if !is_bns_vt_c_gr(bns.vert[vc as usize].type_) {
        vc ^= i32::from(pe.neighbor12);
    }
    let v0 = vc ^ i32::from(pe.neighbor12);
    let pv0 = &bns.vert[v0 as usize];
    if pv0.type_ & BNS_VERT_TYPE_ATOM != 0 {
        return NO_VERTEX;
    }
    let mut pe1 = [0usize; 3];
    let mut v1 = [0i32; 3];
    let (mut i, mut j) = (0usize, 0usize);
    while i < pv0.num_adj_edges as usize && j < 3 {
        pe1[j] = bns.iedge[pv0.iedge + i] as usize;
        v1[j] = i32::from(bns.edge[pe1[j]].neighbor12) ^ v0;
        let t = bns.vert[v1[j] as usize].type_;
        if vc != v1[j] && t & BNS_VERT_TYPE_ATOM == 0 && !is_bns_vt_c_gr(t) {
            j += 1;
        }
        i += 1;
    }
    if j != 2 || i != pv0.num_adj_edges as usize {
        return NO_VERTEX;
    }
    let n = |k: usize| bns.vert[v1[k] as usize].num_adj_edges;
    let (i0, i1) = if n(1) == 2 && n(0) == 3 {
        (1, 0)
    } else if n(0) == 2 && n(1) == 3 {
        (0, 1)
    } else {
        return NO_VERTEX;
    };
    let pv_t = &bns.vert[v1[i1] as usize];
    let mut k = 0;
    for i in 0..pv_t.num_adj_edges as usize {
        let pe_t = &bns.edge[bns.iedge[pv_t.iedge + i] as usize];
        let v = i32::from(pe_t.neighbor12) ^ v1[i1];
        if v == v0 {
            k += 1;
        }
        if v == v1[i0] {
            k += 2;
        }
        if bns.vert[v as usize].type_ & BNS_VERT_TYPE_ATOM != 0 {
            k += 4;
        }
    }
    if k != 7 {
        return NO_VERTEX;
    }
    pe1[i0] as i32
}

// ---------------------------------------------------------------------------
// Explicit terminal H

/// `ConnectDisconnectedH`: connects the explicit H atoms stored after the
/// `num_atoms` atoms; returns the total number of atoms or an error
pub fn connect_disconnected_h(at: &mut [InpAtom], num_atoms: usize, num_deleted_h: usize) -> i32 {
    let tot_atoms = num_atoms + num_deleted_h;
    let mut i = num_atoms;
    while i < tot_atoms {
        let k = at[i].neighbor[0] as usize;
        let mut j = i;
        while j < tot_atoms && at[j].neighbor[0] as usize == k {
            j += 1;
        }
        let num_h = j - i;
        if num_h as i32 > i32::from(at[k].num_h) {
            return RI_ERR_PROGR;
        }
        if num_h + at[k].valence as usize > MAXVAL {
            return RI_ERR_SYNTAX;
        }
        let n = at[k].valence as usize;
        at[k].neighbor.copy_within(0..n, num_h);
        at[k].bond_stereo.copy_within(0..n, num_h);
        at[k].bond_type.copy_within(0..n, num_h);
        for q in 0..num_h {
            at[k].neighbor[q] = (i + q) as u16;
            at[k].bond_stereo[q] = 0;
            at[k].bond_type[q] = BOND_TYPE_SINGLE as u8;
        }
        let mut m = 0;
        while m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m] != 0 {
            at[k].sb_ord[m] = at[k].sb_ord[m].wrapping_add(num_h as i8);
            if at[k].sn_ord[m] < 0 {
                let mut n = i;
                while n < j {
                    if at[n].orig_at_number == at[k].sn_orig_at_num[m] {
                        at[k].sn_ord[m] = (n - i) as i8;
                        break;
                    }
                    n += 1;
                }
                if n == j {
                    return RI_ERR_PROGR;
                }
            } else {
                at[k].sn_ord[m] = at[k].sn_ord[m].wrapping_add(num_h as i8);
            }
            m += 1;
        }
        at[k].valence = at[k].valence.wrapping_add(num_h as i8);
        at[k].chem_bonds_valence = at[k].chem_bonds_valence.wrapping_add(num_h as i8);
        at[k].num_h = at[k].num_h.wrapping_sub(num_h as i8);
        for a in &mut at[i..j] {
            a.chem_bonds_valence = BOND_TYPE_SINGLE as i8;
        }
        let mut m = j as isize - 1;
        while i as isize <= m && at[m as usize].iso_atw_diff > 0 {
            let d = at[m as usize].iso_atw_diff as usize;
            if d > NUM_H_ISOTOPES {
                return RI_ERR_PROGR;
            }
            let old = at[k].num_iso_h[d - 1];
            at[k].num_iso_h[d - 1] = old.wrapping_sub(1);
            if old <= 0 {
                return RI_ERR_PROGR;
            }
            m -= 1;
        }
        i = j;
    }
    for a in &mut at[..num_atoms] {
        for m in 0..NUM_H_ISOTOPES {
            a.num_h = a.num_h.wrapping_sub(a.num_iso_h[m]);
        }
        if a.num_h < 0 {
            return RI_ERR_PROGR;
        }
    }
    tot_atoms as i32
}

/// `DisconnectedConnectedH`: the reverse of [`connect_disconnected_h`]
pub fn disconnected_connected_h(at: &mut [InpAtom], num_atoms: usize, num_deleted_h: usize) -> i32 {
    let tot_atoms = num_atoms + num_deleted_h;
    for a in &mut at[..num_atoms] {
        for m in 0..NUM_H_ISOTOPES {
            a.num_h = a.num_h.wrapping_add(a.num_iso_h[m]);
        }
    }
    let mut i = num_atoms;
    while i < tot_atoms {
        let k = at[i].neighbor[0] as usize;
        let mut j = i;
        while j < tot_atoms && at[j].neighbor[0] as usize == k {
            at[j].chem_bonds_valence = 0;
            j += 1;
        }
        let num_h = j - i;
        let mut n = 0usize;
        while n < at[k].valence as usize && at[k].neighbor[n] as usize >= num_atoms {
            n += 1;
        }
        if n != num_h {
            return RI_ERR_PROGR;
        }
        at[k].valence = at[k].valence.wrapping_sub(num_h as i8);
        let n = at[k].valence as usize;
        at[k].chem_bonds_valence = at[k].chem_bonds_valence.wrapping_sub(num_h as i8);
        if n != 0 {
            at[k].neighbor.copy_within(num_h..num_h + n, 0);
            at[k].bond_stereo.copy_within(num_h..num_h + n, 0);
            at[k].bond_type.copy_within(num_h..num_h + n, 0);
        }
        for q in n..n + num_h {
            at[k].neighbor[q] = 0;
            at[k].bond_stereo[q] = 0;
            at[k].bond_type[q] = 0;
        }
        let mut m = 0;
        while m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m] != 0 {
            at[k].sb_ord[m] = at[k].sb_ord[m].wrapping_sub(num_h as i8);
            if 0 <= at[k].sn_ord[m] && (at[k].sn_ord[m] as usize) < num_h {
                at[k].sn_ord[m] = -1;
            }
            m += 1;
        }
        let mut m = j as isize - 1;
        while i as isize <= m && at[m as usize].iso_atw_diff > 0 {
            let d = at[m as usize].iso_atw_diff as usize;
            if d > NUM_H_ISOTOPES {
                return RI_ERR_PROGR;
            }
            at[k].num_iso_h[d - 1] = at[k].num_iso_h[d - 1].wrapping_add(1);
            m -= 1;
        }
        at[k].num_h = at[k].num_h.wrapping_add(num_h as i8);
        i = j;
    }
    tot_atoms as i32
}
