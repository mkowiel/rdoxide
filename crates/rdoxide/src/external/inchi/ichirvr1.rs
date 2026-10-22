//! The InChI -> structure restoration, part 1 (`INCHI_BASE/src/ichirvr1.c`):
//! the charge-structure tables, the t-groups of a read InChI, the charge and
//! valence states of an atom.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::InChI;
use super::ichi_bns::BNS_VERT_TYPE_ATOM;
use super::ichicano::{T_GROUP_HDR_LEN, TGSO_SYMM_IORDER};
use super::ichiread::{NO_VALUE_INT, RI_ERR_ALLOC, RI_ERR_PROGR};
use super::ichirvrs::*;
use super::ichisort::insertions_sort;
use super::ichitaut::{TGroup, TGroupInfo};
use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE};
use super::inpdef::{BOND_TYPE_MASK, InpAtom};
use super::util::{
    EL_NUMBER_N, EL_NUMBER_O, EL_NUMBER_P, EL_NUMBER_S, EL_NUMBER_SE, EL_NUMBER_TE, get_el_valence,
    if_skip_add_h,
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
