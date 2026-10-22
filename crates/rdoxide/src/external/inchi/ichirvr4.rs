//! The InChI -> structure restoration, part 4 (`INCHI_BASE/src/ichirvr4.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi_bns::{BnData, BnStruct, NO_VERTEX};
use super::ichimake::Icr;
use super::ichiread::RI_ERR_PROGR;
use super::ichirvr1::{
    BnsTestResult, CN_BITS_MNP, CN_LIST, EdgeList, add_to_edge_list, alloc_edge_list,
    get_charge_flower_upper_edge, remove_forbidden_edge_mask, run_bns_restore_once,
    run_bns_test_once,
};
use super::ichirvr2::StrFromInChI;
use super::ichirvrs::{AllTcGroups, TCG_MINUS_C0, TCG_PLUS_C0, ValAt};
use super::ichitaut::{BOND_ALT12NS, is_centerpoint_elem};
use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE};
use super::inpdef::{BOND_TYPE_MASK, InpAtom};

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

/// `ForbidNintrogenPlus2BondsInSmallRings`: forbids the edges that allow
/// =N(+)= or #N(+)- in a small ring
#[allow(clippy::too_many_arguments)]
pub fn forbid_nitrogen_plus2_bonds_in_small_rings(
    bns: &mut BnStruct,
    at: &[InpAtom],
    num_at: usize,
    p_va: &[ValAt],
    min_ring_size: i32,
    n_plus2_bonds_edges: &mut EdgeList,
    forbidden_edge_mask: i32,
) -> i32 {
    let mask = forbidden_edge_mask as i8;
    for i in 0..num_at {
        let va = &p_va[i];
        if at[i].valence == 2
            && at[i].num_h == 0
            && at[i].endpoint == 0
            && va.c_num_valence_electrons == 5
            && va.c_periodic_row_number == 1
            && va.c_max_flow_to_metal == 0
            && va.n_c_plus_group_edge > 0
            && va.c_cn_list_index > 0
            && CN_LIST[va.c_cn_list_index as usize - 1].bits == CN_BITS_MNP
            && va.c_min_ring_size != 0
            && i32::from(va.c_min_ring_size) <= min_ring_size
        {
            let j = va.n_c_plus_group_edge - 1;
            if bns.edge[j as usize].forbidden & mask == 0 {
                bns.edge[j as usize].forbidden |= mask;
                let r = add_to_edge_list(n_plus2_bonds_edges, j, 128);
                if r != 0 {
                    return r;
                }
            }
        }
    }
    0
}

/// fixes all the charge edges; returns the error of the list
fn fix_all_charges(
    bns: &mut BnStruct,
    p_va: &[ValAt],
    num_at: usize,
    list: &mut EdgeList,
    mask: i8,
) -> i32 {
    let ret = alloc_edge_list(list, 2 * num_at as i32);
    if ret != 0 {
        return ret;
    }
    for va in &p_va[..num_at] {
        for j in [va.n_c_minus_group_edge - 1, va.n_c_plus_group_edge - 1] {
            if j >= 0 {
                let r = add_to_edge_list(list, j, 0);
                if r != 0 {
                    return r;
                }
                bns.edge[j as usize].forbidden |= mask;
            }
        }
    }
    0
}

/// the common test and undo of a one-charge-edge flow removal used by the
/// Fix* functions: `Ok(ret)` of the restore or the rollback
#[allow(clippy::too_many_arguments)]
fn try_charge_move(
    bns: &mut BnStruct,
    bd: &mut BnData,
    p_va: &[ValAt],
    pe: usize,
    max_delta_charge: i32,
    num_run_bns: &mut i32,
    total_delta: &mut i32,
) -> i32 {
    let mut res = BnsTestResult::default();
    let v1 = bns.edge[pe].neighbor1 as usize;
    let v2 = (bns.edge[pe].neighbor12 as usize) ^ v1;
    bns.edge[pe].flow -= 1;
    bns.vert[v1].st_edge.flow -= 1;
    bns.vert[v2].st_edge.flow -= 1;
    bns.tot_st_flow -= 2;
    let mut ret = run_bns_test_once(bns, bd, p_va, &mut res);
    if ret < 0 {
        return ret;
    }
    let (a, b) = (v1 as i32, v2 as i32);
    if ret == 1
        && ((res.v_last == a && res.v_first == b) || (res.v_last == b && res.v_first == a))
        && res.delta_charge <= max_delta_charge
    {
        ret = run_bns_restore_once(bns, bd);
        *num_run_bns += 1;
        if ret < 0 {
            return ret;
        }
        if ret != 0 {
            *total_delta += ret;
        } else {
            ret = RI_ERR_PROGR;
        }
        // marks the success for the callers
        return if ret == RI_ERR_PROGR {
            ret
        } else {
            1_000_000 + ret
        };
    }
    bns.edge[pe].flow += 1;
    bns.vert[v1].st_edge.flow += 1;
    bns.vert[v2].st_edge.flow += 1;
    bns.tot_st_flow += 2;
    0
}

/// `FixLessHydrogenInFormula`: the formula of the restored structure has
/// less H than the input InChI
#[allow(clippy::too_many_arguments)]
pub fn fix_less_hydrogen_in_formula(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &StrFromInChI,
    at: &[InpAtom],
    at2: &mut [InpAtom],
    atf: &[InpAtom],
    p_va: &[ValAt],
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let mask = forbidden_edge_mask as i8;
    let mut fixed = EdgeList::default();
    let mut i_bplus = NO_VERTEX;
    let mut i_nv = NO_VERTEX;
    let mut i_nh = NO_VERTEX;
    let mut res = BnsTestResult::default();
    let mut ret = fix_all_charges(bns, p_va, num_at, &mut fixed, mask);
    'body: {
        if ret != 0 {
            break 'body;
        }
        for i in 0..num_at {
            if at2[i].endpoint == 0
                && p_va[i].c_metal == 0
                && p_va[i].c_num_valence_electrons == 5
                && p_va[i].c_periodic_row_number == 1
                && i32::from(at2[i].num_h) == i32::from(atf[i].num_h) + 1
            {
                i_nh = i as i32;
                break;
            }
        }
        if !(0 <= i_nh && (i_nh as usize) < num_at) {
            break 'body;
        }
        let nh = i_nh as usize;
        for j in 0..at2[nh].valence as usize {
            let neigh = at2[nh].neighbor[j] as usize;
            if at2[neigh].valence == 4
                && at2[neigh].charge == -1
                && at2[neigh].chem_bonds_valence == 4
                && at2[neigh].radical == 0
                && at[neigh].num_h == 0
            {
                i_bplus = neigh as i32;
            }
        }
        let mut delta: i32 = -1;
        let mut nxt = i_nh;
        let mut nxt_is_nplus = false;
        loop {
            let i = nxt as usize;
            nxt = NO_VERTEX;
            delta = -delta;
            for j in 0..at2[i].valence as usize {
                let bt2 = i32::from(at2[i].bond_type[j]) & BOND_TYPE_MASK;
                let btf = i32::from(atf[i].bond_type[j]) & BOND_TYPE_MASK;
                nxt_is_nplus = false;
                let nb = at2[i].neighbor[j] as usize;
                if (btf - bt2 == delta || btf == BOND_ALT12NS)
                    && BOND_TYPE_SINGLE <= bt2 + delta
                    && bt2 + delta <= BOND_TYPE_TRIPLE
                    && at2[nb].c_flags == 0
                {
                    let prv = i;
                    nxt = nb as i32;
                    nxt_is_nplus = at2[nb].charge == 1
                        && atf[nb].charge == 0
                        && p_va[nb].c_num_valence_electrons == 5
                        && p_va[nb].c_periodic_row_number == 1;
                    at2[i].c_flags |= 1;
                    if delta == -1
                        && at2[prv].valence == 4
                        && at2[prv].chem_bonds_valence == 5
                        && at2[prv].charge == 0
                        && at2[prv].radical == 0
                        && p_va[prv].c_num_valence_electrons == 5
                        && p_va[prv].n_c_plus_group_edge > 0
                    {
                        i_nv = prv as i32;
                    }
                    if at2[nb].charge != atf[nb].charge {
                        if (at2[nb].charge == 1 || atf[nb].charge == 1)
                            && p_va[nb].n_c_plus_group_edge > 0
                        {
                            bns.edge[(p_va[nb].n_c_plus_group_edge - 1) as usize].forbidden &=
                                !mask;
                        }
                        if (at2[nb].charge == -1 || atf[nb].charge == -1)
                            && p_va[nb].n_c_minus_group_edge > 0
                        {
                            bns.edge[(p_va[nb].n_c_minus_group_edge - 1) as usize].forbidden &=
                                !mask;
                        }
                    }
                    break;
                }
            }
            if !(nxt >= 0 && !(nxt_is_nplus && delta == -1)) {
                break;
            }
        }
        for a in &mut at2[..num_at] {
            a.c_flags = 0;
        }
        if !(nxt >= 0 && nxt_is_nplus && delta == -1) {
            break 'body;
        }
        if i_bplus == 0 && i_nv == 0 {
            break 'body;
        }
        let mut jsel: Option<usize> = None;
        for cand in [i_nv, i_nh] {
            if cand >= 0 {
                let t = p_va[cand as usize].n_c_plus_group_edge - 1;
                if t > 0 && bns.edge[t as usize].flow > 0 {
                    jsel = Some(t as usize);
                    break;
                }
            }
        }
        let Some(pe) = jsel else { break 'body };
        let v1 = bns.edge[pe].neighbor1 as usize;
        let v2 = (bns.edge[pe].neighbor12 as usize) ^ v1;
        let nxt_u = nxt as usize;
        let fl = get_charge_flower_upper_edge(bns, p_va[nxt_u].n_c_plus_group_edge - 1);
        if fl >= 0 && bns.edge[fl as usize].flow == 1 {
            bns.edge[fl as usize].forbidden |= mask;
            ret = add_to_edge_list(&mut fixed, fl, 0);
            if ret != 0 {
                break 'body;
            }
        }
        bns.edge[pe].forbidden |= mask;
        bns.edge[pe].flow -= 1;
        bns.vert[v1].st_edge.flow -= 1;
        bns.vert[v2].st_edge.flow -= 1;
        bns.tot_st_flow -= 2;
        ret = run_bns_test_once(bns, bd, p_va, &mut res);
        if ret < 0 {
            break 'body;
        }
        let (a, b) = (v1 as i32, v2 as i32);
        if ret == 1
            && ((res.v_last == a && res.v_first == b) || (res.v_last == b && res.v_first == a))
            && res.delta_charge <= 0
        {
            ret = run_bns_restore_once(bns, bd);
            *num_run_bns += 1;
            if ret < 0 {
                break 'body;
            } else if ret == 1 {
                *total_delta += ret;
            } else {
                ret = RI_ERR_PROGR;
                break 'body;
            }
        } else {
            ret = 0;
            bns.edge[pe].flow += 1;
            bns.vert[v1].st_edge.flow += 1;
            bns.vert[v2].st_edge.flow += 1;
            bns.tot_st_flow += 2;
        }
    }
    remove_forbidden_edge_mask(bns, &fixed, forbidden_edge_mask);
    ret
}

/// `FixMoreHydrogenInFormula`: X=Y-O(-) => X(-)-Y=O
#[allow(clippy::too_many_arguments)]
pub fn fix_more_hydrogen_in_formula(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &StrFromInChI,
    _at: &[InpAtom],
    at2: &[InpAtom],
    atf: &[InpAtom],
    p_va: &[ValAt],
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let mask = forbidden_edge_mask as i8;
    let mut fixed = EdgeList::default();
    let mut ret = fix_all_charges(bns, p_va, num_at, &mut fixed, mask);
    'body: {
        if ret != 0 {
            break 'body;
        }
        for i in 0..num_at {
            let neigh = at2[i].neighbor[0] as usize;
            let ep = if ps.b_mobile_h != 0 {
                at2[i].endpoint != 0
            } else {
                ps.endpoint.as_ref().is_some_and(|e| e[i] != 0)
            };
            if ep || p_va[i].c_metal != 0 || i32::from(at2[i].num_h) + 1 != i32::from(atf[i].num_h)
            {
                continue;
            }
            let k = p_va[i].n_c_minus_group_edge - 1;
            if !(k >= 0
                && bns.edge[k as usize].flow == 1
                && at2[i].charge == -1
                && atf[i].charge == 0
                && at2[i].valence == 1
                && at2[i].chem_bonds_valence == 1
                && p_va[i].c_num_valence_electrons == 6
                && at2[neigh].chem_bonds_valence > at2[neigh].valence)
            {
                continue;
            }
            let i_nh = i;
            let mut k2 = 0;
            let mut j = 0;
            while j < at2[neigh].valence as usize {
                let neigh2 = at2[neigh].neighbor[j] as usize;
                if neigh2 != i_nh && at2[neigh2].endpoint == 0 {
                    let ie = bns.iedge[bns.vert[neigh].iedge + j] as usize;
                    if bns.edge[ie].forbidden == 0
                        && (4..=5).contains(&p_va[neigh2].c_num_valence_electrons)
                    {
                        k2 = p_va[neigh2].n_c_minus_group_edge - 1;
                        if k2 >= 0 && bns.edge[k2 as usize].flow == 0 {
                            break;
                        }
                    }
                }
                j += 1;
            }
            if j < at2[neigh].valence as usize {
                bns.edge[k2 as usize].forbidden &= !mask;
                let r = try_charge_move(bns, bd, p_va, k as usize, 1, num_run_bns, total_delta);
                if r < 0 {
                    ret = r;
                    break 'body;
                }
                if r >= 1_000_000 || r == RI_ERR_PROGR {
                    ret = if r == RI_ERR_PROGR { r } else { r - 1_000_000 };
                    break 'body;
                }
                ret = 0;
            }
        }
    }
    remove_forbidden_edge_mask(bns, &fixed, forbidden_edge_mask);
    ret
}

/// `FixRemoveExtraTautEndpoints`: false tautomerism of -NH(=O)OH
#[allow(clippy::too_many_arguments)]
pub fn fix_remove_extra_taut_endpoints(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &StrFromInChI,
    at2: &[InpAtom],
    atn: &[InpAtom],
    p_va: &[ValAt],
    icr: &Icr,
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let mask = forbidden_edge_mask as i8;
    let mut fixed = EdgeList::default();
    let mut ret = fix_all_charges(bns, p_va, num_at, &mut fixed, mask);
    'body: {
        if ret != 0 {
            break 'body;
        }
        for &e in &icr.endp_in1_only {
            let endpoint1 = e as usize - 1;
            if at2[endpoint1].valence == at2[endpoint1].chem_bonds_valence
                || p_va[endpoint1].n_c_minus_group_edge <= 0
            {
                continue;
            }
            for j in 0..at2[endpoint1].valence as usize {
                if BOND_TYPE_DOUBLE != (BOND_TYPE_MASK & i32::from(at2[endpoint1].bond_type[j])) {
                    continue;
                }
                let cp = at2[endpoint1].neighbor[j] as usize;
                if at2[cp].charge != 0
                    || p_va[cp].n_c_plus_group_edge <= 0
                    || !is_centerpoint_elem(at2[cp].el_number)
                {
                    continue;
                }
                let mut k = 0;
                while k < at2[cp].valence as usize {
                    if BOND_TYPE_SINGLE == (BOND_TYPE_MASK & i32::from(at2[cp].bond_type[k])) {
                        let endpoint2 = at2[cp].neighbor[k] as usize;
                        if at2[endpoint2].endpoint == 0 && atn[endpoint2].endpoint != 0 {
                            break;
                        }
                    }
                    k += 1;
                }
                if k == at2[cp].valence as usize {
                    continue;
                }
                let pe = (p_va[cp].n_c_plus_group_edge - 1) as usize;
                if bns.edge[pe].flow == 0 {
                    continue;
                }
                let pe2 = (p_va[endpoint1].n_c_minus_group_edge - 1) as usize;
                if bns.edge[pe2].flow != 0 {
                    continue;
                }
                bns.edge[pe2].forbidden &= !mask;
                let r = try_charge_move(bns, bd, p_va, pe, 1, num_run_bns, total_delta);
                if r < 0 {
                    ret = r;
                    break 'body;
                }
                if r >= 1_000_000 || r == RI_ERR_PROGR {
                    ret = if r == RI_ERR_PROGR { r } else { r - 1_000_000 };
                    break 'body;
                }
                ret = 0;
                bns.edge[pe2].forbidden |= mask;
            }
        }
    }
    remove_forbidden_edge_mask(bns, &fixed, forbidden_edge_mask);
    ret
}
