//! The InChI -> structure restoration, part 4 (`INCHI_BASE/src/ichirvr4.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi_bns::BnStruct;
use super::ichiread::RI_ERR_PROGR;
use super::ichirvr1::{EdgeList, add_to_edge_list, alloc_edge_list};
use super::ichirvrs::{AllTcGroups, TCG_MINUS_C0, TCG_PLUS_C0};

/// `ForbidCarbonChargeEdges`: forbids the edges that charge carbon atoms;
/// returns their number or an error
pub fn forbid_carbon_charge_edges(
    bns: &mut BnStruct,
    tcg: &AllTcGroups,
    carbon_charge_edges: &mut EdgeList,
    forbidden_edge_mask: i32,
) -> i32 {
    const MAX_NUM_CARBON_CHARGE_EDGES: i32 = 2;
    let ret = alloc_edge_list(carbon_charge_edges, MAX_NUM_CARBON_CHARGE_EDGES);
    if ret != 0 {
        return ret;
    }
    carbon_charge_edges.num_edges = 0;
    for n_type in [TCG_PLUS_C0, TCG_MINUS_C0] {
        let k = tcg.n_group[n_type as usize];
        if k >= 0 {
            let k = tcg.p_tcg[k as usize].n_forward_edge;
            if k > 0 {
                let e = &mut bns.edge[k as usize];
                if i32::from(e.forbidden) & forbidden_edge_mask == 0 {
                    e.forbidden |= forbidden_edge_mask as i8;
                    let ret = add_to_edge_list(carbon_charge_edges, k, 0);
                    if ret != 0 {
                        return ret;
                    }
                }
            } else {
                return RI_ERR_PROGR;
            }
        }
    }
    carbon_charge_edges.num_edges
}
