//! The InChI -> structure restoration, part 5 (`INCHI_BASE/src/ichirvr5.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi_bns::{BnStruct, NO_VERTEX};
use super::ichirvrs::{AllTcGroups, TCG_MINUS, TCG_PLUS};

/// `GetPlusMinusVertex`: the (+/-) vertex connected to the super (+) or (-)
/// charge group, or `NO_VERTEX`. (As in C, the second vertex is computed
/// from the edge index instead of the super vertex number.)
pub fn get_plus_minus_vertex(
    bns: &BnStruct,
    tcg: &AllTcGroups,
    check_forbidden_plus: bool,
    check_forbidden_minus: bool,
) -> i32 {
    let mut v_plus_minus1 = NO_VERTEX;
    let mut v_plus_minus2 = NO_VERTEX;
    let k = tcg.n_group[TCG_PLUS as usize];
    if k >= 0 {
        let g = &tcg.p_tcg[k as usize];
        let e_plus_super = g.n_forward_edge;
        let v_plus_super = g.n_vertex_number;
        if e_plus_super > 0
            && v_plus_super >= bns.num_atoms
            && !(bns.edge[e_plus_super as usize].forbidden != 0 && check_forbidden_plus)
        {
            v_plus_minus1 = i32::from(bns.edge[e_plus_super as usize].neighbor12) ^ v_plus_super;
        }
    }
    let k = tcg.n_group[TCG_MINUS as usize];
    if k >= 0 {
        let g = &tcg.p_tcg[k as usize];
        let e_minus_super = g.n_forward_edge;
        if e_minus_super > 0
            && g.n_vertex_number >= bns.num_atoms
            && !(bns.edge[e_minus_super as usize].forbidden != 0 && check_forbidden_minus)
        {
            v_plus_minus2 = i32::from(bns.edge[e_minus_super as usize].neighbor12) ^ e_minus_super;
        }
    }
    if (check_forbidden_plus && v_plus_minus1 == NO_VERTEX)
        || (check_forbidden_minus && v_plus_minus2 == NO_VERTEX)
    {
        return NO_VERTEX;
    }
    if v_plus_minus1 != NO_VERTEX {
        v_plus_minus1
    } else {
        v_plus_minus2
    }
}
