//! The InChI -> structure restoration, part 4 (`INCHI_BASE/src/ichirvr4.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::InChI;
use super::ichi_bns::{
    BNS_EDGE_FORBIDDEN_MASK, BNS_EDGE_FORBIDDEN_TEST, BnData, BnStruct, NO_VERTEX,
    set_forbidden_edges,
};
use super::ichimake::Icr;
use super::ichiread::RI_ERR_PROGR;
use super::ichirvr1::{
    BnsTestResult, CN_BITS_MNP, CN_LIST, EdgeList, add_to_edge_list, alloc_edge_list,
    get_charge_flower_upper_edge, get_tgroup_info_from_inchi, is_bns_vt_c_gr,
    remove_forbidden_bond_flow_bits, remove_forbidden_edge_mask, run_bns_restore_once,
    run_bns_test_once, set_forbidden_edge_mask,
};
use super::ichirvr2::{StrFromInChI, refresh_at2};
use super::ichirvr5::get_plus_minus_vertex;
use super::ichirvrs::{
    AllTcGroups, CN_BITS_M, CN_BITS_N, CN_BITS_P, CN_BITS_SHIFT, TCG_MINUS_C0, TCG_PLUS_C0, ValAt,
    make_cn_bits,
};
use super::ichitaut::{BOND_ALT12NS, is_centerpoint_elem};
use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE};
use super::inpdef::{BOND_TYPE_MASK, InpAtom};
use super::strutil::is_metal_salt;
use super::util::get_endpoint_valence;

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

/// `MoveChargeToRemoveCenerpoints`: finds and eliminates false Mobile-H
/// groups: Cl(=O)3(-O(-)) => Cl(-)(=O)4
#[allow(clippy::too_many_arguments)]
pub fn move_charge_to_remove_centerpoints(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &mut [ValAt],
    tcg: &AllTcGroups,
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let test = i32::from(BNS_EDGE_FORBIDDEN_TEST);
    let b_possibly_ignore =
        ps.charge >= 0 && (tcg.num_tgroups == 0 || (ps.i_mobile_h == 0 && ps.ti.num_t_groups != 0));
    let mut charge_list = EdgeList::default();
    let mut res = BnsTestResult::default();
    let mut num_success = 0;
    let r2 = refresh_at2(ps, at2, bns, p_va, tcg);
    if r2 < 0 {
        return r2;
    }
    let r2 = super::ichinorm::mark_ring_systems_inp(&mut at2[..num_at], 0);
    if r2 < 0 {
        return r2;
    }
    let ret_forbid_edges = set_forbidden_edges(bns, at2, num_at as i32, test);
    if ret_forbid_edges < 0 {
        return ret_forbid_edges;
    }
    let endp = |ps: &StrFromInChI, v: usize| ps.endpoint.as_ref().is_some_and(|e| e[v] != 0);
    for i in 0..num_at {
        let a = &at2[i];
        let va = p_va[i];
        if !(va.c_num_valence_electrons != 4
            && !(va.n_taut_group_edge != 0 || (ps.i_mobile_h == 0 && endp(ps, i)))
            && a.num_h == 0
            && a.charge == 0
            && a.valence >= 2
            && a.valence < a.chem_bonds_valence
            && is_centerpoint_elem(a.el_number))
        {
            continue;
        }
        let is_centerpoint_n =
            va.c_num_valence_electrons == 5 && (va.c_periodic_row_number == 1 || va.c_metal != 0);
        let mut mobile_charge_neigh: Vec<usize> = Vec::new();
        let mut double_bond_acceptors = 0;
        let mut double_bond_not_o_neigh: Vec<usize> = Vec::new();
        let (mut num_donors, mut num_acceptors, mut num_donors_o, mut num_acceptors_o) =
            (0, 0, 0, 0);
        let (mut num_known_endpoints, mut num_wrong_neigh) = (0, 0);
        let mut j = 0;
        while j < a.valence as usize {
            let neigh = a.neighbor[j] as usize;
            let nb = &at2[neigh];
            if nb.endpoint != 0 || (ps.i_mobile_h == 0 && endp(ps, neigh)) || nb.charge > 0 {
                num_known_endpoints += 1;
                j += 1;
                continue;
            }
            if i32::from(bns.edge[bns.iedge[bns.vert[i].iedge + j] as usize].forbidden) & test != 0
            {
                j += 1;
                continue;
            }
            let bond_type = i32::from(a.bond_type[j]) & BOND_TYPE_MASK;
            if bond_type > BOND_TYPE_DOUBLE {
                num_wrong_neigh += 1;
                j += 1;
                continue;
            }
            if nb.num_h != 0 && bond_type == BOND_TYPE_SINGLE {
                break;
            }
            if i32::from(nb.chem_bonds_valence) - i32::from(nb.charge)
                != get_endpoint_valence(nb.el_number)
            {
                if bond_type == BOND_TYPE_DOUBLE && p_va[neigh].c_num_valence_electrons != 6 {
                    double_bond_not_o_neigh.push(j);
                }
                j += 1;
                continue;
            }
            if nb.charge == -1
                && bond_type == BOND_TYPE_SINGLE
                && (p_va[neigh].n_c_minus_group_edge < 1
                    || bns.edge[(p_va[neigh].n_c_minus_group_edge - 1) as usize].flow != 1)
            {
                break;
            }
            let o_like =
                p_va[neigh].c_num_valence_electrons == 6 && p_va[neigh].c_periodic_row_number <= 4;
            if bond_type == BOND_TYPE_SINGLE {
                if nb.charge != -1 || p_va[neigh].n_c_minus_group_edge <= 0 {
                    num_wrong_neigh += 1;
                    j += 1;
                    continue;
                }
                num_donors += 1;
                num_donors_o += i32::from(o_like);
                mobile_charge_neigh.push(j);
            } else if bond_type == BOND_TYPE_DOUBLE {
                if nb.charge != 0 {
                    num_wrong_neigh += 1;
                    j += 1;
                    continue;
                }
                double_bond_acceptors += 1;
                num_acceptors += 1;
                num_acceptors_o += i32::from(o_like);
            }
            j += 1;
        }
        let _ = double_bond_acceptors;
        if j != a.valence as usize || num_donors == 0 || num_acceptors == 0 {
            continue;
        }
        if is_centerpoint_n && num_donors == num_donors_o && num_acceptors == num_acceptors_o {
            continue;
        }
        if ps.i_mobile_h == 0 && num_donors as usize == double_bond_not_o_neigh.len() {
            if charge_list.num_edges == 0 {
                for k in 0..num_at {
                    let ak = &at2[k];
                    let e0 = bns.iedge[bns.vert[k].iedge] as usize;
                    let em = p_va[k].n_c_minus_group_edge - 1;
                    let ep = p_va[k].n_c_plus_group_edge - 1;
                    let skip = (ak.valence == 1
                        && bns.edge[e0].flow != 0
                        && bns.edge[e0].forbidden == 0
                        && !(em >= 0 && bns.edge[em as usize].flow != 0)
                        && !(ep >= 0 && bns.edge[ep as usize].flow == 0)
                        && p_va[k].c_num_valence_electrons == 6
                        && p_va[k].c_metal == 0
                        && endp(ps, k))
                        || ps.fixed_h.as_ref().is_some_and(|f| f[k] != 0);
                    if !skip
                        && em >= 0
                        && bns.edge[em as usize].flow == 0
                        && bns.edge[em as usize].forbidden == 0
                    {
                        let r = add_to_edge_list(&mut charge_list, em, 64);
                        if r != 0 {
                            return r;
                        }
                    }
                    if ep >= 0 && bns.edge[ep as usize].forbidden == 0 {
                        let r = add_to_edge_list(&mut charge_list, ep, 64);
                        if r != 0 {
                            return r;
                        }
                    }
                }
            }
            for &jj in &double_bond_not_o_neigh {
                let e = bns.iedge[bns.vert[i].iedge + jj];
                if bns.edge[e as usize].forbidden == 0 {
                    let r = add_to_edge_list(&mut charge_list, e, 64);
                    if r != 0 {
                        return r;
                    }
                }
            }
            set_forbidden_edge_mask(bns, &charge_list, forbidden_edge_mask);
            let mut num_moved = 0usize;
            let mut k = 0;
            while k < double_bond_not_o_neigh.len() && num_moved < mobile_charge_neigh.len() {
                let pe = bns.iedge[bns.vert[i].iedge + double_bond_not_o_neigh[k]] as usize;
                k += 1;
                if bns.edge[pe].flow != 1 {
                    continue;
                }
                let v1 = bns.edge[pe].neighbor1 as usize;
                let v2 = (bns.edge[pe].neighbor12 as usize) ^ v1;
                bns.vert[v1].st_edge.flow -= 1;
                bns.vert[v2].st_edge.flow -= 1;
                bns.edge[pe].flow -= 1;
                bns.tot_st_flow -= 2;
                let ret = run_bns_test_once(bns, bd, p_va, &mut res);
                if ret < 0 {
                    return ret;
                }
                let (x, y) = (v1 as i32, v2 as i32);
                if ret == 1
                    && ((res.v_last == x && res.v_first == y)
                        || (res.v_last == y && res.v_first == x))
                    && res.delta_charge == 0
                {
                    let ret = run_bns_restore_once(bns, bd);
                    *num_run_bns += 1;
                    if ret < 0 {
                        return ret;
                    }
                    if ret == 1 {
                        *total_delta += ret;
                        num_moved += 1;
                    } else {
                        return RI_ERR_PROGR;
                    }
                } else {
                    bns.vert[v1].st_edge.flow += 1;
                    bns.vert[v2].st_edge.flow += 1;
                    bns.edge[pe].flow += 1;
                    bns.tot_st_flow += 2;
                }
            }
            remove_forbidden_edge_mask(bns, &charge_list, forbidden_edge_mask);
        } else if !b_possibly_ignore
            || (num_known_endpoints == 0
                && num_wrong_neigh == 0
                && num_acceptors_o + num_donors_o >= 3)
        {
            bns.vert[i].st_edge.cap += num_donors;
            bns.tot_st_cap += num_donors;
            p_va[i].c_init_charge = p_va[i].c_init_charge.wrapping_sub(num_donors as i8);
            for &jm in &mobile_charge_neigh {
                let neigh = at2[i].neighbor[jm] as usize;
                let pm = (p_va[neigh].n_c_minus_group_edge - 1) as usize;
                let v1m = bns.edge[pm].neighbor1 as usize;
                let v2m = (bns.edge[pm].neighbor12 as usize) ^ v1m;
                let delta = bns.edge[pm].flow;
                bns.vert[v1m].st_edge.flow -= delta;
                bns.vert[v2m].st_edge.flow -= delta;
                if is_bns_vt_c_gr(bns.vert[v1m].type_) {
                    bns.vert[v1m].st_edge.cap -= delta;
                } else if is_bns_vt_c_gr(bns.vert[v2m].type_) {
                    bns.vert[v2m].st_edge.cap -= delta;
                } else {
                    return RI_ERR_PROGR;
                }
                bns.tot_st_cap -= delta;
                bns.tot_st_flow -= 2 * delta;
                bns.edge[pm].flow -= delta;
            }
            let ret = run_bns_restore_once(bns, bd);
            *num_run_bns += 1;
            if ret < 0 {
                return ret;
            }
            if ret == num_donors {
                *total_delta += ret;
                num_success += 1;
            } else {
                return RI_ERR_PROGR;
            }
        }
    }
    if ret_forbid_edges != 0 {
        remove_forbidden_bond_flow_bits(bns, test);
    }
    num_success
}

/// `MakeSingleBondsMetal2ChargedHeteroat`: bonds of a metal to a charged
/// heteroatom that are not single get flow removed and the BNS is rerun
#[allow(clippy::too_many_arguments)]
pub fn make_single_bonds_metal2_charged_heteroat(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let len_at = (ps.num_atoms + ps.num_deleted_h) as usize;
    let mask = forbidden_edge_mask as i8;
    let r2 = refresh_at2(ps, at2, bns, p_va, tcg);
    if r2 < 0 {
        return r2;
    }
    let mut fixed: Vec<usize> = Vec::new();
    for i in 0..num_at {
        if p_va[i].c_metal == 0 {
            continue;
        }
        for j in 0..at2[i].valence as usize {
            let neigh = at2[i].neighbor[j] as usize;
            if p_va[neigh].c_num_valence_electrons == 4 && p_va[neigh].c_periodic_row_number == 1 {
                continue;
            }
            if i32::from(at2[i].bond_type[j]) > BOND_TYPE_SINGLE
                && at2[neigh].charge != 0
                && p_va[neigh].c_metal == 0
                && p_va[neigh].c_cn_list_index > 0
            {
                let cn_bits = if at2[neigh].charge > 0 {
                    make_cn_bits(CN_BITS_N, CN_BITS_P, 0, 0)
                } else {
                    make_cn_bits(CN_BITS_N, CN_BITS_M, 0, 0)
                };
                let mut at_bits = CN_LIST[p_va[neigh].c_cn_list_index as usize - 1].bits;
                let mut k = 0;
                while k < 3 {
                    if at_bits & cn_bits == cn_bits {
                        break;
                    }
                    k += 1;
                    at_bits >>= CN_BITS_SHIFT;
                }
                if k == 3 {
                    continue;
                }
                fixed.push(bns.iedge[bns.vert[i].iedge + j] as usize);
            }
        }
    }
    at2.clear();
    at2.extend_from_slice(&ps.at[..len_at]);
    if fixed.is_empty() {
        return 0;
    }
    for &e in &fixed {
        let v1 = bns.edge[e].neighbor1 as usize;
        let v2 = (bns.edge[e].neighbor12 as usize) ^ v1;
        bns.edge[e].flow -= 1;
        bns.edge[e].forbidden |= mask;
        bns.vert[v1].st_edge.flow -= 1;
        bns.vert[v2].st_edge.flow -= 1;
        bns.tot_st_flow -= 2;
        *total_delta -= 2;
    }
    let mut ret = run_bns_restore_once(bns, bd);
    *num_run_bns += 1;
    if ret < 0 {
        return ret;
    }
    *total_delta += ret;
    for &e in &fixed {
        bns.edge[e].forbidden &= !mask;
    }
    if ret < 2 * fixed.len() as i32 {
        ret = run_bns_restore_once(bns, bd);
        *num_run_bns += 1;
        if ret < 0 {
            return ret;
        }
        *total_delta += ret;
    }
    ret
}

/// `SaltBondsToCoordBonds`: in the reconnected structure changes 'salt
/// bonds' to 'coordination bonds', for example M-O-C= -> M(+)-O(-)-C=
#[allow(clippy::too_many_arguments)]
pub fn salt_bonds_to_coord_bonds(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &mut [ValAt],
    tcg: &AllTcGroups,
    num_run_bns: &mut i32,
    _total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let num_edges = bns.num_bonds + 2 * bns.num_atoms;
    let mask = forbidden_edge_mask as i8;
    let mut all_charge = EdgeList::default();
    let mut res = BnsTestResult::default();
    let mut ret = 0;
    if ps.i_inchi_rec == 0 || ps.srm.b_metal_add_flower == 0 || ps.srm.n_metal_min_bond_order != 0 {
        return 0;
    }
    let r2 = refresh_at2(ps, at2, bns, p_va, tcg);
    if r2 < 0 {
        return r2;
    }
    macro_rules! add {
        ($e:expr) => {{
            ret = add_to_edge_list(&mut all_charge, $e, num_edges);
            if ret != 0 {
                return ret;
            }
        }};
    }
    for i in 0..num_at {
        if !is_metal_salt(at2, i) {
            continue;
        }
        if all_charge.num_edges == 0 {
            for j in 0..num_at {
                for k in 0..at2[j].valence as usize {
                    let n = at2[j].neighbor[k] as usize;
                    let ie = bns.iedge[bns.vert[j].iedge + k];
                    if n < j && bns.edge[ie as usize].forbidden == 0 {
                        add!(ie);
                    }
                }
                for ie in [
                    p_va[j].n_c_minus_group_edge - 1,
                    p_va[j].n_c_plus_group_edge - 1,
                ] {
                    if ie >= 0 && bns.edge[ie as usize].forbidden == 0 {
                        add!(ie);
                    }
                }
            }
            for j in 0..tcg.num_tgroups as usize {
                let v1 = tcg.p_tcg[j].n_vertex_number as usize;
                for k in 0..bns.vert[v1].num_adj_edges as usize {
                    add!(bns.iedge[bns.vert[v1].iedge + k]);
                }
            }
        }
        for k in 0..at2[i].valence as usize {
            let n = at2[i].neighbor[k] as usize;
            let pe = bns.iedge[bns.vert[i].iedge + k] as usize;
            if at2[n].charge != 0 || i32::from(at2[i].bond_type[k]) != BOND_TYPE_SINGLE {
                continue;
            }
            let ie_plus_metal = p_va[i].n_c_plus_group_edge - 1;
            let ie_minus_o = p_va[n].n_c_minus_group_edge - 1;
            if bns.edge[pe].flow != 1 || bns.edge[pe].forbidden != 0 || ie_plus_metal < 0 {
                continue;
            }
            let pe_plus_metal = ie_plus_metal as usize;
            if bns.edge[pe_plus_metal].flow <= 0 {
                continue;
            }
            if ie_minus_o >= 0 {
                let pe_minus_o = ie_minus_o as usize;
                if bns.edge[pe_minus_o].flow != 0
                    || bns.edge[pe_plus_metal].forbidden != 0
                    || bns.edge[pe_minus_o].forbidden != 0
                {
                    continue;
                }
                let v1 = bns.edge[pe].neighbor1 as usize;
                let v2 = (bns.edge[pe].neighbor12 as usize) ^ v1;
                bns.edge[pe].flow -= 1;
                bns.vert[v1].st_edge.flow -= 1;
                bns.vert[v2].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                set_forbidden_edge_mask(bns, &all_charge, forbidden_edge_mask);
                bns.edge[pe_plus_metal].forbidden &= !mask;
                bns.edge[pe_minus_o].forbidden &= !mask;
                ret = run_bns_test_once(bns, bd, p_va, &mut res);
                let (a, b) = (v1 as i32, v2 as i32);
                if ret == 1
                    && ((res.v_last == a && res.v_first == b)
                        || (res.v_last == b && res.v_first == a))
                {
                    ret = run_bns_restore_once(bns, bd);
                    if ret > 0 {
                        *num_run_bns += 1;
                    }
                } else {
                    bns.edge[pe].flow += 1;
                    bns.vert[v1].st_edge.flow += 1;
                    bns.vert[v2].st_edge.flow += 1;
                    bns.tot_st_flow += 2;
                }
                remove_forbidden_edge_mask(bns, &all_charge, forbidden_edge_mask);
            } else {
                let vpm = get_plus_minus_vertex(bns, tcg, true, true);
                if vpm == NO_VERTEX {
                    continue;
                }
                let vpm = vpm as usize;
                bns.vert[vpm].st_edge.cap += 1;
                bns.vert[n].st_edge.cap += 1;
                bns.edge[pe].cap += 1;
                bns.tot_st_cap += 2;
                let (v1, v2) = (vpm as i32, n as i32);
                for i1 in 0..bns.vert[vpm].num_adj_edges as usize {
                    let j1 = (bns.edge[bns.iedge[bns.vert[vpm].iedge + i1] as usize].neighbor12
                        as usize)
                        ^ vpm;
                    for k1 in 0..bns.vert[j1].num_adj_edges as usize {
                        let e = bns.iedge[bns.vert[j1].iedge + k1] as usize;
                        bns.edge[e].cap += 1;
                    }
                }
                set_forbidden_edge_mask(bns, &all_charge, forbidden_edge_mask);
                bns.edge[pe_plus_metal].forbidden &= !mask;
                bns.edge[pe].forbidden &= !mask;
                ret = run_bns_test_once(bns, bd, p_va, &mut res);
                let mut cur_success = 0;
                if ret == 1
                    && ((res.v_last == v1 && res.v_first == v2)
                        || (res.v_last == v2 && res.v_first == v1))
                {
                    ret = run_bns_restore_once(bns, bd);
                    if ret > 0 {
                        *num_run_bns += 1;
                        cur_success += 1;
                    }
                }
                if cur_success != 0 {
                    if bns.edge[pe].flow != 2 {
                        return RI_ERR_PROGR;
                    }
                    bns.edge[pe].flow -= 2;
                    bns.vert[n].st_edge.cap -= 2;
                    bns.vert[n].st_edge.flow -= 2;
                    bns.vert[i].st_edge.flow -= 2;
                    bns.vert[i].st_edge.cap -= 2;
                    bns.tot_st_cap -= 3;
                    bns.tot_st_flow -= 4;
                    bns.edge[pe].cap -= 2;
                    p_va[n].c_init_charge -= 1;
                } else {
                    bns.vert[vpm].st_edge.cap -= 1;
                    bns.vert[n].st_edge.cap -= 1;
                    bns.tot_st_cap -= 2;
                    for i1 in 0..bns.vert[vpm].num_adj_edges as usize {
                        let j1 = (bns.edge[bns.iedge[bns.vert[vpm].iedge + i1] as usize].neighbor12
                            as usize)
                            ^ vpm;
                        for k1 in 0..bns.vert[j1].num_adj_edges as usize {
                            let e = bns.iedge[bns.vert[j1].iedge + k1] as usize;
                            bns.edge[e].cap -= 1;
                        }
                    }
                }
                remove_forbidden_edge_mask(bns, &all_charge, forbidden_edge_mask);
            }
        }
    }
    ret
}

/// `FillOutExtraFixedHDataInChI`: allocates and fills the endpoint and
/// fixed-H data of the Fixed-H structure from the input InChI (`p_inchi[0]`
/// is the Fixed-H and `p_inchi[1]` the Mobile-H InChI of the component)
pub fn fill_out_extra_fixed_h_data_inchi(
    ps: &mut StrFromInChI,
    p_inchi: [Option<&InChI>; 2],
) -> i32 {
    let n = ps.num_atoms as usize;
    ps.endpoint = Some(vec![0; n]);
    let fixed = ps.fixed_h.get_or_insert_with(|| vec![0; n]);
    fixed.resize(n, 0);
    get_tgroup_info_from_inchi(&mut ps.ti, None, ps.endpoint.as_deref_mut(), p_inchi[1]);
    match p_inchi[0].and_then(|c| c.n_num_h_fixed.as_ref()) {
        Some(src) => {
            for (d, s) in fixed.iter_mut().zip(src.iter()) {
                *d = *s;
            }
        }
        None => fixed.fill(0),
    }
    0
}

/// `CheckAndRefixStereobonds`: finds A=X< where all bonds to X except A=X
/// are marked as stereogenic; temporarily allows the stereobonds to change
/// and makes the A=X bonds single
#[allow(clippy::too_many_arguments)]
pub fn check_and_refix_stereobonds(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let stereo = i32::from(BNS_EDGE_FORBIDDEN_MASK);
    let num_at = ps.num_atoms as usize;
    let mask = forbidden_edge_mask as i8;
    let r2 = refresh_at2(ps, at2, bns, p_va, tcg);
    if r2 < 0 {
        return 0;
    }
    // the atom has an A=X< double bond and exactly two stereo bonds
    let wrong = |a: &InpAtom| -> Option<(usize, usize)> {
        let (j1, j2) = (a.sb_ord[0] as usize, a.sb_ord[1] as usize);
        (a.valence == 3
            && i32::from(a.chem_bonds_valence) - i32::from(a.valence) == 1
            && a.sb_parity[0] != 0
            && a.sb_parity[1] != 0
            && a.sb_parity[2] == 0
            && j1 < a.bond_type.len()
            && j2 < a.bond_type.len()
            && (i32::from(a.bond_type[j1]) & BOND_TYPE_MASK) == BOND_TYPE_SINGLE
            && (i32::from(a.bond_type[j2]) & BOND_TYPE_MASK) == BOND_TYPE_SINGLE
            && j1 != j2)
            .then_some((j1, j2))
    };
    let num_wrong = at2[..num_at].iter().filter(|a| wrong(a).is_some()).count() as i32;
    if num_wrong == 0 {
        return 0;
    }
    let num_fixed = bns.edge[..bns.num_bonds as usize]
        .iter()
        .filter(|e| i32::from(e.forbidden) & stereo != 0)
        .count() as i32;
    let mut carbon = EdgeList::default();
    let mut fixed = EdgeList::default();
    let mut wrong_edges = EdgeList::default();
    let mut ret;
    'body: {
        ret = forbid_carbon_charge_edges(bns, tcg, &mut carbon, forbidden_edge_mask);
        if ret < 0 {
            break 'body;
        }
        ret = alloc_edge_list(&mut fixed, num_fixed);
        if ret == 0 {
            ret = alloc_edge_list(&mut wrong_edges, num_wrong);
        }
        if ret != 0 {
            break 'body;
        }
        let mut i = 0;
        while i < num_at && wrong_edges.num_edges < num_wrong {
            if let Some((j1, j2)) = wrong(&at2[i]) {
                let k = match j1 + j2 {
                    1 => 2,
                    2 => 1,
                    3 => 0,
                    _ => {
                        ret = RI_ERR_PROGR;
                        break 'body;
                    }
                };
                let ne = bns.iedge[bns.vert[i].iedge + k];
                let pe = ne as usize;
                let v1 = bns.edge[pe].neighbor1 as usize;
                let v2 = (bns.edge[pe].neighbor12 as usize) ^ v1;
                if bns.edge[pe].flow == 0 {
                    ret = RI_ERR_PROGR;
                    break 'body;
                }
                bns.edge[pe].flow -= 1;
                bns.edge[pe].forbidden |= mask;
                bns.vert[v1].st_edge.flow -= 1;
                bns.vert[v2].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                ret = add_to_edge_list(&mut wrong_edges, ne, 0);
                if ret != 0 {
                    break 'body;
                }
            }
            i += 1;
        }
        let mut i = 0;
        while i < bns.num_bonds as usize && fixed.num_edges < num_fixed {
            if i32::from(bns.edge[i].forbidden) & stereo != 0 {
                bns.edge[i].forbidden &= !(stereo as i8);
                fixed.pn_edges[fixed.num_edges as usize] = i as i32;
                fixed.num_edges += 1;
            }
            i += 1;
        }
        let ret_bns = run_bns_restore_once(bns, bd);
        *num_run_bns += 1;
        if ret_bns < 0 {
            break 'body;
        }
        if ret_bns > 0 {
            *total_delta += ret_bns;
        }
        remove_forbidden_edge_mask(bns, &wrong_edges, forbidden_edge_mask);
        remove_forbidden_edge_mask(bns, &carbon, forbidden_edge_mask);
        set_forbidden_edge_mask(bns, &fixed, stereo);
        let ret2 = run_bns_restore_once(bns, bd);
        *num_run_bns += 1;
        if ret2 < 0 {
            break 'body;
        }
        if ret2 > 0 {
            *total_delta += ret_bns;
        }
        ret = ret_bns;
    }
    ret
}
