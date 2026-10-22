//! The InChI -> structure restoration, part 5 (`INCHI_BASE/src/ichirvr5.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::InChI;
use super::ichi_bns::{BnData, BnStruct, NO_VERTEX};
use super::ichiparm::InputParms;
use super::ichiring::is_bond_in_nmax_memb_ring;
use super::ichirvr1::{
    BfsQ, BnsTestResult, EdgeList, add_to_edge_list, find_in_edge_list,
    get_charge_flower_upper_edge, make_one_inchi_out_of_str_from_inchi2,
    remove_forbidden_edge_mask, run_bns_restore_once, run_bns_test_once, set_forbidden_edge_mask,
};
use super::ichirvr2::StrFromInChI;
use super::ichirvr4::{
    Cmp2MhInchi, MAX_DIFF_FIXH, MH_DB_N_MINUS, MH_NH, MH_NH_MINUS, MH_NH2, MH_O_MINUS, MH_OH,
    fill_out_cmp2mh_inchi, fill_out_extra_fixed_h_data_restr,
};
use super::ichirvrs::{AllTcGroups, TCG_MINUS, TCG_PLUS, TGRF_MINUS_FIRST, ValAt};
use super::incomdef::BOND_TYPE_SINGLE;
use super::inpdef::InpAtom;

const INC_ADD_EDGE: i32 = 64;
const RI_ERR_PROGR: i32 = super::ichiread::RI_ERR_PROGR;

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

/// `bIsUnsatCarbonInASmallRing`
pub fn b_is_unsat_carbon_in_a_small_ring(
    at2: &[InpAtom],
    p_va: &[ValAt],
    iat: usize,
    pbfsq: &mut BfsQ,
    min_ring_size: i32,
) -> i32 {
    let a = &at2[iat];
    if min_ring_size < 5 {
        // =C= in a small ring
        if a.valence == 2 && p_va[iat].c_min_ring_size <= 5 && a.chem_bonds_valence == 4 {
            return 1;
        }
    } else {
        if a.valence == 2
            && p_va[iat].c_min_ring_size != 0
            && i32::from(p_va[iat].c_min_ring_size) <= min_ring_size
            && a.chem_bonds_valence == 3
        {
            return 1;
        }
        let mut cur = min_ring_size + 1;
        let mut min = min_ring_size + 1;
        if (a.valence == 2 || a.valence == 3) && a.chem_bonds_valence == a.valence + 1 {
            for j in 0..a.valence as usize {
                let q = pbfsq.q.as_mut().expect("BFS queue");
                cur = is_bond_in_nmax_memb_ring(
                    at2,
                    iat,
                    j,
                    q,
                    &mut pbfsq.n_atom_level,
                    &mut pbfsq.c_source,
                    min as u16,
                );
                if 0 < cur && cur < min {
                    min = cur;
                }
            }
            return if cur >= 0 {
                i32::from(min <= min_ring_size)
            } else {
                cur
            };
        }
    }
    0
}

/// the state of `FixMobileHRestoredStructure`
#[allow(missing_debug_implementations)]
struct FixMhCtx<'a> {
    ip: &'a InputParms,
    bns: &'a mut BnStruct,
    bd: &'a mut BnData,
    ps: &'a mut StrFromInChI,
    at2: &'a mut Vec<InpAtom>,
    at3: &'a mut Vec<InpAtom>,
    p_va: &'a [ValAt],
    tcg: &'a mut AllTcGroups,
    p_inchi: [Option<&'a InChI>; 2],
    mask: i32,
    mask_inv: i8,
    c2i: Cmp2MhInchi,
    all_charge: EdgeList,
    curr: EdgeList,
    curr2: EdgeList,
    taut: EdgeList,
    ret: i32,
    res: BnsTestResult,
}

impl FixMhCtx<'_> {
    fn minus_edge(&self, iat: usize) -> i32 {
        self.p_va[iat].n_c_minus_group_edge - 1
    }

    fn plus_edge(&self, iat: usize) -> i32 {
        self.p_va[iat].n_c_plus_group_edge - 1
    }

    fn free_edge(&self, e: i32) -> bool {
        e >= 0 && self.bns.edge[e as usize].forbidden == 0
    }

    fn vert_edge(&self, v: usize, j: usize) -> usize {
        self.bns.iedge[self.bns.vert[v].iedge + j] as usize
    }

    fn path_between(&self, a: usize, b: usize) -> bool {
        let (a, b) = (a as i32, b as i32);
        (self.res.v_last == a && self.res.v_first == b)
            || (self.res.v_last == b && self.res.v_first == a)
    }

    fn fill_cmp(&mut self) -> bool {
        self.ret = fill_out_cmp2mh_inchi(
            self.ps,
            self.tcg,
            self.at2,
            self.p_va,
            self.p_inchi,
            &mut self.c2i,
        );
        self.ret != 0
    }

    fn no_fixed_h(&self) -> bool {
        self.p_inchi[0].is_none_or(|c| c.n_num_h_fixed.is_none())
            && self.ps.p_one_inchi[0]
                .as_ref()
                .is_none_or(|c| c.n_num_h_fixed.is_none())
    }

    /// recalculates the InChI of the structure and compares it with the
    /// original one; true = exit the function
    fn recalc(&mut self, check_fixed_h: bool) -> bool {
        self.ret = make_one_inchi_out_of_str_from_inchi2(
            self.ip, self.bns, self.ps, self.at2, self.at3, self.p_va, self.tcg,
        );
        if self.ret < 0 {
            return true;
        }
        self.ret = fill_out_extra_fixed_h_data_restr(self.ps);
        if self.ret != 0 {
            return true;
        }
        if check_fixed_h && self.no_fixed_h() {
            return true;
        }
        if self.fill_cmp() {
            return true;
        }
        !self.c2i.b_has_difference
    }

    /// `(Rem)(i) = tg->tg_num_H--, tg_num_Minus++` bookkeeping of the cases
    /// 03 and 04 when an NH is replaced with N(-)
    fn nh_to_n_minus(&mut self, itg: usize, set_minus: Option<usize>) {
        let g = &mut self.tcg.p_tcg[itg];
        g.tg_num_h -= 1;
        g.tg_num_minus += 1;
        g.tg_restore_flags |= TGRF_MINUS_FIRST;
        if let Some(v) = set_minus {
            g.tg_set_minus = (v + 1) as i16;
        }
        self.ps.ti.t_group.as_mut().expect("t-groups")[itg].num[1] += 1;
        self.tcg.total_charge -= 1;
        self.tcg.tgroup_charge -= 1;
        self.ps.n_num_removed_protons_by_revrs += 1;
    }

    /// the reverse of `nh_to_n_minus`
    fn n_minus_to_nh(&mut self, itg: usize, reset_set_minus: bool) {
        let g = &mut self.tcg.p_tcg[itg];
        g.tg_num_h += 1;
        g.tg_num_minus -= 1;
        g.tg_restore_flags &= !TGRF_MINUS_FIRST;
        if reset_set_minus {
            g.tg_set_minus = 0;
        }
        self.ps.ti.t_group.as_mut().expect("t-groups")[itg].num[1] -= 1;
        self.tcg.total_charge += 1;
        self.tcg.tgroup_charge += 1;
        self.ps.n_num_removed_protons_by_revrs -= 1;
    }

    /// the beginning of the C function; true = exit
    fn prologue(&mut self, forbidden_stereo_edge_mask: i32) -> bool {
        if self.no_fixed_h() {
            return true;
        }
        for i in 0..self.tcg.num_tgroups as usize {
            let v1 = self.tcg.p_tcg[i].n_vertex_number as usize;
            for j in 0..self.bns.vert[v1].num_adj_edges as usize {
                let e = self.bns.iedge[self.bns.vert[v1].iedge + j];
                self.ret = add_to_edge_list(&mut self.taut, e, INC_ADD_EDGE);
                if self.ret != 0 {
                    return true;
                }
            }
        }
        let mut n_flower = EdgeList::default();
        let mut other_n_flower = EdgeList::default();
        let mut fixed_large_ring = EdgeList::default();
        for i in 0..self.ps.num_atoms as usize {
            let e = self.minus_edge(i);
            if self.free_edge(e) {
                self.ret = add_to_edge_list(&mut self.all_charge, e, INC_ADD_EDGE);
                if self.ret != 0 {
                    return true;
                }
            }
            let e = self.plus_edge(i);
            if self.free_edge(e) {
                self.ret = add_to_edge_list(&mut self.all_charge, e, INC_ADD_EDGE);
                if self.ret != 0 {
                    return true;
                }
                let va = &self.p_va[i];
                if va.c_num_valence_electrons == 5 && va.c_metal == 0 {
                    let j = get_charge_flower_upper_edge(self.bns, e);
                    if j != NO_VERTEX {
                        let ju = j as usize;
                        if self.bns.edge[ju].forbidden == 0 && self.bns.edge[ju].flow != 0 {
                            self.ret = add_to_edge_list(&mut self.all_charge, j, INC_ADD_EDGE);
                            if self.ret == 0 {
                                self.ret = add_to_edge_list(&mut n_flower, j, INC_ADD_EDGE);
                            }
                        } else {
                            self.ret = add_to_edge_list(&mut other_n_flower, j, INC_ADD_EDGE);
                        }
                        if self.ret != 0 {
                            return true;
                        }
                    }
                }
            }
        }
        if forbidden_stereo_edge_mask != 0 {
            for i in 0..self.ps.num_atoms as usize {
                for j in 0..self.at2[i].valence as usize {
                    let k = self.vert_edge(i, j);
                    if i32::from(self.bns.edge[k].forbidden) == forbidden_stereo_edge_mask {
                        let q = self.ps.pbfsq.q.as_mut().expect("BFS queue");
                        let min_ring = is_bond_in_nmax_memb_ring(
                            self.at2,
                            i,
                            j,
                            q,
                            &mut self.ps.pbfsq.n_atom_level,
                            &mut self.ps.pbfsq.c_source,
                            99,
                        );
                        if 0 < min_ring {
                            self.ret =
                                add_to_edge_list(&mut fixed_large_ring, k as i32, INC_ADD_EDGE);
                            if self.ret != 0 {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        self.ret = fill_out_extra_fixed_h_data_restr(self.ps);
        if self.ret != 0 {
            return true;
        }
        self.fill_cmp()
    }

    /// the common ADP condition of the cases 01-05
    fn adp(&self) -> bool {
        let c = &self.c2i;
        c.n_num_tg_inchi == 1 && (c.n_num_endp_revrs < c.n_num_endp_inchi || c.n_num_tg_revrs > 1)
    }

    fn case_01(&mut self) -> bool {
        let c = &self.c2i;
        if !(self.adp()
            && c.revrs[MH_DB_N_MINUS] + c.revrs[MH_NH_MINUS] == 0
            && c.inchi[MH_O_MINUS] != 0
            && self
                .tcg
                .p_tcg
                .first()
                .is_some_and(|g| g.tg_restore_flags & TGRF_MINUS_FIRST == 0))
        {
            return false;
        }
        self.tcg.p_tcg[0].tg_restore_flags |= TGRF_MINUS_FIRST;
        self.recalc(true)
    }

    fn case_02(&mut self) -> bool {
        let c = &self.c2i;
        if !(self.adp()
            && c.revrs[MH_DB_N_MINUS] + c.revrs[MH_NH_MINUS] == 0
            && c.inchi[MH_O_MINUS] == 0)
        {
            return false;
        }
        const MAX: usize = MAX_DIFF_FIXH;
        let mut sb: Vec<usize> = Vec::new();
        let mut db: Vec<usize> = Vec::new();
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        for iat in 0..self.ps.num_atoms as usize {
            let a = &self.at2[iat];
            let va = &self.p_va[iat];
            if !(va.c_num_valence_electrons == 5
                && va.c_periodic_row_number == 1
                && a.endpoint == 0
                && a.charge == 0
                && a.radical == 0)
            {
                continue;
            }
            let em = self.minus_edge(iat);
            let ep = self.plus_edge(iat);
            let minus_ok = self.free_edge(em)
                && self.bns.edge[em as usize].cap != 0
                && self.bns.edge[em as usize].flow == 0;
            if db.len() < MAX
                && a.num_h == 0
                && a.valence == 2
                && a.chem_bonds_valence == 3
                && a.sb_parity[0] == 0
                && minus_ok
            {
                db.push(iat);
            } else if db.len() < MAX
                && a.num_h == 1
                && a.valence == 1
                && a.chem_bonds_valence == 2
                && a.sb_parity[0] == 0
                && minus_ok
            {
                db.push(iat);
            } else if sb.len() < MAX
                && a.num_h == 0
                && a.valence == 3
                && a.chem_bonds_valence == 3
                && self.free_edge(ep)
                && self.bns.edge[ep as usize].cap != 0
                && self.bns.edge[ep as usize].flow != 0
            {
                sb.push(iat);
                self.ret = add_to_edge_list(&mut self.curr, ep, INC_ADD_EDGE);
                if self.ret != 0 {
                    return true;
                }
            }
        }
        if !db.is_empty() && !sb.is_empty() {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut i = 0;
            while i < db.len() && cur_success == 0 {
                let iat = db[i];
                i += 1;
                let k1 = self.vert_edge(iat, 0);
                let e = if self.bns.edge[k1].flow != 0 {
                    k1 as i32
                } else {
                    let k2 = self.vert_edge(iat, 1);
                    if self.bns.edge[k2].flow != 0 {
                        k2 as i32
                    } else {
                        NO_VERTEX
                    }
                };
                if e == NO_VERTEX {
                    continue;
                }
                let pe_minus = self.minus_edge(iat) as usize;
                let pe = e as usize;
                if self.bns.edge[pe].flow == 0 {
                    continue;
                }
                let v1 = self.bns.edge[pe].neighbor1 as usize;
                let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
                self.bns.edge[pe].forbidden |= self.mask as i8;
                self.bns.edge[pe_minus].forbidden &= self.mask_inv;
                self.bns.edge[pe].flow -= 1;
                self.bns.vert[v1].st_edge.flow -= 1;
                self.bns.vert[v2].st_edge.flow -= 1;
                self.bns.tot_st_flow -= 2;
                self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                if self.ret == 1 && self.path_between(v1, v2) && self.res.delta_charge == 2 {
                    self.ret = run_bns_restore_once(self.bns, self.bd);
                    if self.ret > 0 {
                        cur_success += 1;
                        // eliminate the (-) charge and add H
                        let w1 = self.bns.edge[pe_minus].neighbor1 as usize;
                        let w2 = (self.bns.edge[pe_minus].neighbor12 as usize) ^ w1;
                        self.bns.edge[pe_minus].cap -= 1;
                        self.bns.edge[pe_minus].flow -= 1;
                        self.bns.vert[w1].st_edge.cap -= 1;
                        self.bns.vert[w1].st_edge.flow -= 1;
                        self.bns.vert[w2].st_edge.cap -= 1;
                        self.bns.vert[w2].st_edge.flow -= 1;
                        self.bns.tot_st_flow -= 2;
                        self.bns.tot_st_cap -= 2;
                        self.ps.at[iat].num_h += 1;
                        self.tcg.total_charge += 1;
                        self.ps.n_num_removed_protons_by_revrs -= 1;
                    }
                } else {
                    self.bns.edge[pe].forbidden &= self.mask_inv;
                    self.bns.edge[pe_minus].forbidden |= self.mask as i8;
                    self.bns.edge[pe].flow += 1;
                    self.bns.vert[v1].st_edge.flow += 1;
                    self.bns.vert[v2].st_edge.flow += 1;
                    self.bns.tot_st_flow += 2;
                }
            }
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            self.curr.num_edges = 0;
            if cur_success != 0 {
                return self.recalc(false);
            }
        }
        false
    }

    /// cases 03 and 04
    fn case_03(&mut self) -> bool {
        let c = &self.c2i;
        if !(self.adp()
            && c.n_num_tg_m_inchi == 0
            && c.n_num_tg_n_inchi != 0
            && c.n_num_tg_o_inchi != 0)
        {
            return false;
        }
        let (nh, nh2, oh) = (c.inchi[MH_NH], c.inchi[MH_NH2], c.inchi[MH_OH]);
        let has_o = c.n_num_tg_o_inchi != 0;
        self.curr.num_edges = 0;
        self.curr2.num_edges = 0;
        let mut cur_success = 0;
        let mut b_action = 0;
        let num_tg = self.tcg.num_tgroups as usize;
        // 1st attempt: -NH-=O => -N(-)-=O or -N=-OH => -N(-)-=O
        let mut itg = 0;
        while itg < num_tg && cur_success == 0 {
            let v_tg = self.tcg.p_tcg[itg].n_vertex_number as usize;
            let mut i = 0;
            while i < self.bns.vert[v_tg].num_adj_edges as usize && cur_success == 0 {
                let e_tg = self.vert_edge(v_tg, i);
                i += 1;
                let v_endp = self.bns.edge[e_tg].neighbor1 as usize;
                let mut found: Option<(usize, usize, usize)> = None; // (eTg2, eCent1, eCent2)
                if self.p_va[v_endp].c_num_valence_electrons == 6 && self.bns.edge[e_tg].cap != 0 {
                    let mut j1 = 0;
                    while j1 < self.at2[v_endp].valence as usize && found.is_none() {
                        let e_cent1 = self.vert_edge(v_endp, j1);
                        j1 += 1;
                        let v_cent = (self.bns.edge[e_cent1].neighbor12 as usize) ^ v_endp;
                        let pc1 = self.bns.edge[e_cent1];
                        let ptg = self.bns.edge[e_tg];
                        if self.at2[v_cent].endpoint != 0
                            || pc1.cap == 0
                            || pc1.flow + i32::from(ptg.cap == ptg.flow) != 1
                        {
                            continue;
                        }
                        for j2 in 0..self.at2[v_cent].valence as usize {
                            let e_cent2 = self.vert_edge(v_cent, j2);
                            let v_endp2 = (self.bns.edge[e_cent2].neighbor12 as usize) ^ v_cent;
                            let pc2 = self.bns.edge[e_cent2];
                            if pc2.cap == 0
                                || pc2.flow + pc1.flow != 1
                                || self.at2[v_endp2].endpoint as usize != itg + 1
                                || self.p_va[v_endp2].c_num_valence_electrons != 5
                            {
                                continue;
                            }
                            let j = self.p_va[v_endp2].n_taut_group_edge - 1;
                            if j < 0 {
                                continue;
                            }
                            let pt2 = self.bns.edge[j as usize];
                            if pt2.forbidden != 0 || pc2.flow + i32::from(pt2.cap == pt2.flow) != 1
                            {
                                continue;
                            }
                            found = Some((j as usize, e_cent1, e_cent2));
                            break;
                        }
                    }
                }
                if let Some((e_tg2, e_cent1, e_cent2)) = found {
                    if self.bns.edge[e_tg].cap - self.bns.edge[e_tg].flow == 0
                        && self.bns.edge[e_tg2].cap - self.bns.edge[e_tg2].flow == 1
                        && self.bns.edge[e_cent1].flow == 0
                        && self.bns.edge[e_cent2].flow == 1
                    {
                        self.bns.edge[e_tg].flow -= 1;
                        self.bns.edge[e_tg2].flow += 1;
                        self.bns.edge[e_cent2].flow -= 1;
                        self.bns.edge[e_cent1].flow += 1;
                        b_action |= 1;
                    }
                    if self.bns.edge[e_tg].cap - self.bns.edge[e_tg].flow == 1
                        && self.bns.edge[e_tg2].cap - self.bns.edge[e_tg2].flow == 0
                        && self.bns.edge[e_cent1].flow == 1
                        && self.bns.edge[e_cent2].flow == 0
                    {
                        // vEndp2 is the neighbor1 of the second t-group edge
                        let v_endp2 = self.bns.edge[e_tg2].neighbor1 as usize;
                        self.nh_to_n_minus(itg, Some(v_endp2));
                        b_action |= 2;
                        cur_success += 1;
                    }
                }
            }
            itg += 1;
        }
        if nh + nh2 == 0 && oh != 0 && cur_success == 0 {
            // transfer an attachment to N
            for itg in 0..num_tg {
                let v_tg = self.tcg.p_tcg[itg].n_vertex_number as usize;
                for i in 0..self.bns.vert[v_tg].num_adj_edges as usize {
                    let e_tg = self.vert_edge(v_tg, i);
                    let v_endp = self.bns.edge[e_tg].neighbor1 as usize;
                    let ptg = self.bns.edge[e_tg];
                    let at = &self.at2[v_endp];
                    let va = &self.p_va[v_endp];
                    if va.c_num_valence_electrons == 6
                        && at.valence == at.chem_bonds_valence
                        && ptg.flow != 0
                        && ptg.flow == ptg.cap
                    {
                        self.ret = add_to_edge_list(&mut self.curr, e_tg as i32, INC_ADD_EDGE);
                        if self.ret != 0 {
                            return true;
                        }
                    } else if va.c_num_valence_electrons == 5
                        && va.c_periodic_row_number == 1
                        && i32::from(at.valence) + 1 == i32::from(at.chem_bonds_valence)
                        && ptg.cap != 0
                        && ptg.flow + 1 == ptg.cap
                    {
                        // endpoint -N= or =NH found, check for -N=-OH
                        let mut e: i32 = -1;
                        let mut j1 = 0;
                        while j1 < self.at2[v_endp].valence as usize && e < 0 {
                            let e_cent1 = self.vert_edge(v_endp, j1);
                            j1 += 1;
                            if self.bns.edge[e_cent1].flow != 1 {
                                continue;
                            }
                            let v_cent = (self.bns.edge[e_cent1].neighbor12 as usize) ^ v_endp;
                            if self.at2[v_cent].endpoint != 0 {
                                continue;
                            }
                            for j2 in 0..self.at2[v_cent].valence as usize {
                                let e_cent2 = self.vert_edge(v_cent, j2);
                                let v_endp2 = (self.bns.edge[e_cent2].neighbor12 as usize) ^ v_cent;
                                e = self.p_va[v_endp2].n_taut_group_edge - 1;
                                if self.bns.edge[e_cent2].flow != 0
                                    || self.at2[v_endp2].endpoint as usize != itg + 1
                                    || self.p_va[v_endp2].c_num_valence_electrons != 6
                                    || e <= 0
                                    || self.bns.edge[e as usize].forbidden != 0
                                    || self.bns.edge[e as usize].flow == 0
                                {
                                    e = -1;
                                    continue;
                                }
                                self.ret = add_to_edge_list(&mut self.curr2, e, INC_ADD_EDGE);
                                if self.ret != 0 {
                                    return true;
                                }
                                break;
                            }
                        }
                        if e < 0 {
                            self.ret = add_to_edge_list(&mut self.curr, e_tg as i32, INC_ADD_EDGE);
                            if self.ret != 0 {
                                return true;
                            }
                        }
                    }
                }
            }
            // rearrange the flows
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            set_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            set_forbidden_edge_mask(self.bns, &self.curr2, self.mask);
            let list: Option<Vec<i32>> = if self.curr2.num_edges != 0 {
                Some(self.curr2.pn_edges[..self.curr2.num_edges as usize].to_vec())
            } else if self.curr.num_edges != 0 {
                Some(self.curr.pn_edges[..self.curr.num_edges as usize].to_vec())
            } else {
                None
            };
            if let Some(list) = list {
                let mut i = 0;
                while i < list.len() && cur_success == 0 {
                    let pe = list[i] as usize;
                    i += 1;
                    if self.bns.edge[pe].flow == 0 {
                        continue;
                    }
                    let v1 = self.bns.edge[pe].neighbor1 as usize;
                    let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
                    let Some(itg) =
                        (0..num_tg).find(|&g| v2 as i32 == self.tcg.p_tcg[g].n_vertex_number)
                    else {
                        continue;
                    };
                    self.bns.edge[pe].flow -= 1;
                    self.bns.vert[v1].st_edge.flow -= 1;
                    self.bns.vert[v2].st_edge.flow -= 1;
                    self.bns.tot_st_flow -= 2;
                    self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                    if self.ret == 1 && self.path_between(v1, v2) && self.res.delta_charge == 2 {
                        self.ret = run_bns_restore_once(self.bns, self.bd);
                        if self.ret > 0 {
                            cur_success += 1;
                            self.nh_to_n_minus(itg, None);
                            b_action |= 4;
                        }
                    } else {
                        self.bns.edge[pe].flow += 1;
                        self.bns.vert[v1].st_edge.flow += 1;
                        self.bns.vert[v2].st_edge.flow += 1;
                        self.bns.tot_st_flow += 2;
                    }
                }
            }
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
        } else if nh + nh2 != 0 && has_o && cur_success == 0 {
            // change an attachment to N from H to (-)
            let mut itg = 0;
            while itg < num_tg && cur_success == 0 {
                let v_tg = self.tcg.p_tcg[itg].n_vertex_number as usize;
                let mut i = 0;
                while i < self.bns.vert[v_tg].num_adj_edges as usize && cur_success == 0 {
                    let e_tg = self.vert_edge(v_tg, i);
                    i += 1;
                    let v_endp2 = self.bns.edge[e_tg].neighbor1 as usize;
                    let ptg = self.bns.edge[e_tg];
                    let at = &self.at2[v_endp2];
                    let va = &self.p_va[v_endp2];
                    if va.c_num_valence_electrons == 5
                        && va.c_periodic_row_number == 1
                        && at.valence == at.chem_bonds_valence
                        && ptg.flow != 0
                        && ptg.flow == ptg.cap
                    {
                        cur_success += 1;
                        self.nh_to_n_minus(itg, Some(v_endp2));
                        b_action |= 8;
                    }
                }
                itg += 1;
            }
        }
        if cur_success != 0 {
            if self.recalc(false) {
                return true;
            }
            let one_bad = self
                .ps
                .one_ti
                .t_group
                .as_ref()
                .is_some_and(|g| self.ps.one_ti.num_t_groups == 1 && g[0].num[1] != 0);
            if one_bad {
                // this method did not work: no alt path from N(-) to =O
                if b_action & (8 | 2) != 0 {
                    self.n_minus_to_nh(0, true);
                } else if b_action & 4 != 0 {
                    self.n_minus_to_nh(0, false);
                } else {
                    self.ret = RI_ERR_PROGR;
                    return true;
                }
                // (cur_success is not used afterwards)
                if self.recalc_after_rollback() {
                    return true;
                }
            }
            return !self.c2i.b_has_difference;
        }
        false
    }

    /// the second recalculation of cases 03/04 (no `bHasDifference` test)
    fn recalc_after_rollback(&mut self) -> bool {
        self.ret = make_one_inchi_out_of_str_from_inchi2(
            self.ip, self.bns, self.ps, self.at2, self.at3, self.p_va, self.tcg,
        );
        if self.ret < 0 {
            return true;
        }
        self.ret = fill_out_extra_fixed_h_data_restr(self.ps);
        if self.ret != 0 {
            return true;
        }
        self.fill_cmp()
    }

    fn case_05(&mut self) -> bool {
        let c = &self.c2i;
        if !(self.adp()
            && c.n_num_tg_m_inchi == 0
            && (c.n_num_tg_n_inchi != 0 || c.n_num_tg_o_inchi != 0))
        {
            return false;
        }
        let v_plus_minus = get_plus_minus_vertex(self.bns, self.tcg, true, true);
        if v_plus_minus == NO_VERTEX {
            return false;
        }
        let vpm = v_plus_minus as usize;
        self.curr.num_edges = 0;
        self.curr2.num_edges = 0;
        let mut cur_success = 0;
        'body: {
            for i in 0..self.ps.num_atoms as usize {
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let a = &self.at2[iat];
                let va = &self.p_va[i];
                let e = self.plus_edge(iat);
                if a.endpoint == 0
                    && a.charge == 0
                    && a.radical == 0
                    && a.num_h == 0
                    && va.c_num_valence_electrons == 5
                    && va.c_periodic_row_number == 1
                    && e >= 0
                    && self.bns.edge[e as usize].flow != 0
                    && self.bns.edge[e as usize].forbidden == 0
                {
                    self.ret = add_to_edge_list(&mut self.curr, e, INC_ADD_EDGE);
                    if self.ret != 0 {
                        return true;
                    }
                }
            }
            if self.curr.num_edges == 0 {
                break 'body;
            }
            let mut itg = 0;
            while itg < self.tcg.num_tgroups as usize && cur_success == 0 {
                self.curr2.num_edges = 0;
                let v_tg = self.tcg.p_tcg[itg].n_vertex_number as usize;
                for i in 0..self.bns.vert[v_tg].num_adj_edges as usize {
                    let e_tg = self.vert_edge(v_tg, i);
                    let v_endp = self.bns.edge[e_tg].neighbor1 as usize;
                    let ptg = self.bns.edge[e_tg];
                    if self.at2[v_endp].charge != 0
                        || self.at2[v_endp].radical != 0
                        || ptg.cap - ptg.flow != 1
                    {
                        continue;
                    }
                    self.ret = add_to_edge_list(&mut self.curr2, e_tg as i32, INC_ADD_EDGE);
                    if self.ret != 0 {
                        return true;
                    }
                }
                if self.curr2.num_edges == 0 {
                    break 'body;
                }
                set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                set_forbidden_edge_mask(self.bns, &self.taut, self.mask);
                remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
                remove_forbidden_edge_mask(self.bns, &self.curr2, self.mask);
                self.bns.vert[vpm].st_edge.cap += 1;
                self.bns.vert[v_tg].st_edge.cap += 1;
                self.bns.tot_st_cap += 2;
                // increase capacities of edges to Y
                let cap_change = |s: &mut Self, d: i32| {
                    for i in 0..s.bns.vert[vpm].num_adj_edges as usize {
                        let j = (s.bns.edge[s.vert_edge(vpm, i)].neighbor12 as usize) ^ vpm;
                        for k in 0..s.bns.vert[j].num_adj_edges as usize {
                            let ek = s.vert_edge(j, k);
                            s.bns.edge[ek].cap += d;
                        }
                    }
                };
                cap_change(self, 1);
                self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                if self.ret == 1 && self.path_between(vpm, v_tg) && self.res.delta_charge == 1 {
                    self.ret = run_bns_restore_once(self.bns, self.bd);
                    if self.ret > 0 {
                        cur_success += 1;
                        self.tcg.total_charge += 1;
                        self.tcg.p_tcg[itg].edges_cap += 1;
                        self.tcg.p_tcg[itg].tg_num_h += 1;
                        self.ps.n_num_removed_protons_by_revrs -= 1;
                    }
                } else {
                    self.bns.vert[vpm].st_edge.cap -= 1;
                    self.bns.vert[v_tg].st_edge.cap -= 1;
                    self.bns.tot_st_cap -= 2;
                    cap_change(self, -1);
                }
                remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                remove_forbidden_edge_mask(self.bns, &self.taut, self.mask);
                itg += 1;
            }
            if cur_success != 0 {
                return self.recalc(false);
            }
        }
        false
    }

    fn case_06(&mut self) -> bool {
        while self.c2i.n_num_diff_mob_h != 0
            && self.c2i.n_charge_mob_h_revrs > self.c2i.n_charge_mob_h_inchi
        {
            let mut cur_success = 0;
            let mut i = 0;
            while cur_success == 0 && i < self.c2i.c2at.len() {
                let a = self.c2i.c2at[i];
                i += 1;
                let i_o = a.atom_number as usize;
                let e = self.minus_edge(i_o);
                if !(a.n_mob_h_revrs == a.n_mob_h_inchi + 1
                    && a.n_num_h_revrs == a.n_mob_h_inchi
                    && a.endpt_inchi == 0
                    && a.endpt_revrs == 0
                    && self.at2[i_o].charge == -1
                    && e >= 0
                    && self.bns.edge[e as usize].flow != 0)
                {
                    continue;
                }
                let pe = e as usize;
                // first look for =C= in a small ring
                let mut mode = 4;
                while cur_success == 0 && mode <= 8 {
                    if mode == 8 {
                        mode = 99;
                    }
                    let mut iat = 0;
                    while cur_success == 0 && iat < self.ps.num_atoms as usize {
                        let ic = iat;
                        iat += 1;
                        let a2 = &self.at2[ic];
                        let e2 = self.minus_edge(ic);
                        if !(a2.charge == 0
                            && a2.radical == 0
                            && self.p_va[ic].c_num_valence_electrons == 4
                            && e2 >= 0
                            && self.bns.edge[e2 as usize].flow == 0)
                        {
                            continue;
                        }
                        if b_is_unsat_carbon_in_a_small_ring(
                            self.at2,
                            self.p_va,
                            ic,
                            &mut self.ps.pbfsq,
                            mode,
                        ) <= 0
                        {
                            continue;
                        }
                        let pe2 = e2 as usize;
                        set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                        self.bns.edge[pe2].forbidden &= self.mask_inv;
                        if self.bns.edge[pe].flow == 0 {
                            continue;
                        }
                        let v1 = self.bns.edge[pe].neighbor1 as usize;
                        let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
                        self.bns.edge[pe].flow -= 1;
                        self.bns.vert[v1].st_edge.flow -= 1;
                        self.bns.vert[v2].st_edge.flow -= 1;
                        self.bns.tot_st_flow -= 2;
                        self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                        if self.ret == 1 && self.path_between(v1, v2) && self.res.delta_charge == 1
                        {
                            self.ret = run_bns_restore_once(self.bns, self.bd);
                            if self.ret > 0 {
                                cur_success += 1;
                            }
                        } else {
                            self.bns.edge[pe].forbidden |= self.mask as i8;
                            self.bns.edge[pe].flow += 1;
                            self.bns.vert[v1].st_edge.flow += 1;
                            self.bns.vert[v2].st_edge.flow += 1;
                            self.bns.tot_st_flow += 2;
                        }
                        set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                    }
                    mode += 1;
                }
            }
            if cur_success != 0 {
                if self.recalc_after_rollback() {
                    return true;
                }
                if !self.c2i.b_has_difference {
                    return true;
                }
            } else {
                break;
            }
        }
        false
    }

    fn case_07(&mut self) -> bool {
        const WRONG_N: usize = 0;
        const WRONG_O: usize = 1;
        const WRONG_ALL: usize = 2;
        const O_FIXED: usize = 3;
        if !(!self.c2i.c2at.is_empty()
            && self.c2i.n_charge_mob_h_revrs > self.c2i.n_charge_mob_h_inchi)
        {
            return false;
        }
        let mut ch: [EdgeList; 4] = Default::default();
        let mut db_o: Vec<usize> = Vec::new();
        let mut no2: Vec<usize> = Vec::new();
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        let atm: Option<Vec<u16>> = self.ps.p_one_norm_data[0]
            .as_ref()
            .filter(|d| !d.at.is_empty())
            .map(|d| d.at.iter().map(|a| a.endpoint).collect());
        'body: {
            for i in 0..self.ps.num_atoms as usize {
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let a = &self.at2[iat];
                let e = self.minus_edge(iat);
                if !(db_o.len() < MAX_DIFF_FIXH
                    && self.p_va[iat].c_num_valence_electrons == 6
                    && a.endpoint == 0
                    && self.free_edge(e)
                    && a.num_h == 0
                    && !atm.as_ref().is_some_and(|m| m[iat] != 0)
                    && a.charge == 0
                    && a.valence == 1
                    && a.chem_bonds_valence == 2)
                {
                    continue;
                }
                let i_n = a.neighbor[0] as usize;
                let an = &self.at2[i_n];
                if !(self.p_va[i_n].c_num_valence_electrons == 5
                    && an.valence == 3
                    && (an.charge == 0 || an.charge == 1)
                    && i32::from(an.chem_bonds_valence) == 5 - i32::from(an.charge))
                {
                    continue;
                }
                let (mut num_o, mut num_others) = (0, 0);
                for k in 0..an.valence as usize {
                    let neigh = an.neighbor[k] as usize;
                    if neigh == iat {
                        continue;
                    }
                    let ag = &self.at2[neigh];
                    if self.p_va[neigh].c_num_valence_electrons == 6
                        && ag.endpoint == 0
                        && !atm.as_ref().is_some_and(|m| m[neigh] != 0)
                        && ag.valence == 1
                        && ag.num_h == 0
                        && ag.radical == 0
                        && (ag.charge == 0 || ag.charge == -1)
                        && i32::from(ag.chem_bonds_valence) - i32::from(ag.charge) == 2
                    {
                        num_o += 1;
                    } else if i32::from(an.bond_type[k]) == BOND_TYPE_SINGLE
                        && ag.valence > 1
                        && ag.valence < ag.chem_bonds_valence
                    {
                        num_others += 1;
                    }
                }
                if num_o != 1 || num_others != 1 {
                    continue;
                }
                if !no2.contains(&i_n) {
                    no2.push(i_n);
                    db_o.push(iat);
                }
                self.ret = add_to_edge_list(&mut ch[O_FIXED], e, INC_ADD_EDGE);
                if self.ret != 0 {
                    break 'body;
                }
            }
            if !db_o.is_empty() {
                for ci in 0..self.c2i.c2at.len() {
                    let a = self.c2i.c2at[ci];
                    let iat = a.atom_number as usize;
                    let e = self.minus_edge(iat);
                    if a.endpt_revrs != 0
                        && a.endpt_inchi == 0
                        && a.n_at_charge_revrs == -1
                        && self.free_edge(e)
                        && self.bns.edge[e as usize].flow != 0
                        && find_in_edge_list(&ch[O_FIXED], e) < 0
                    {
                        if a.n_val_electr == 6 {
                            self.ret = add_to_edge_list(&mut ch[WRONG_O], e, INC_ADD_EDGE);
                            if self.ret != 0 {
                                break 'body;
                            }
                        } else if a.n_val_electr == 5 {
                            self.ret = add_to_edge_list(&mut ch[WRONG_N], e, INC_ADD_EDGE);
                            if self.ret != 0 {
                                break 'body;
                            }
                        }
                        self.ret = add_to_edge_list(&mut ch[WRONG_ALL], e, INC_ADD_EDGE);
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                }
                for &iat in &db_o {
                    let pe_db_o_minus = self.minus_edge(iat) as usize;
                    let pe = self.vert_edge(iat, 0);
                    if self.bns.edge[pe].flow == 0 {
                        continue;
                    }
                    let v1 = self.bns.edge[pe].neighbor1 as usize;
                    let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
                    self.bns.edge[pe].forbidden |= self.mask as i8;
                    self.bns.edge[pe].flow -= 1;
                    self.bns.vert[v1].st_edge.flow -= 1;
                    self.bns.vert[v2].st_edge.flow -= 1;
                    self.bns.tot_st_flow -= 2;
                    let mut one_success = 0;
                    let mut k = 0;
                    while one_success == 0 && k <= WRONG_ALL {
                        let kk = k;
                        k += 1;
                        if ch[kk].num_edges == 0 {
                            continue;
                        }
                        set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                        remove_forbidden_edge_mask(self.bns, &ch[kk], self.mask);
                        self.bns.edge[pe_db_o_minus].forbidden &= self.mask_inv;
                        self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                        if self.ret == 1 && self.path_between(v1, v2) && self.res.delta_charge == 0
                        {
                            self.ret = run_bns_restore_once(self.bns, self.bd);
                            if self.ret > 0 {
                                one_success += 1;
                            }
                        }
                    }
                    cur_success += one_success;
                    remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                    self.bns.edge[pe].forbidden &= self.mask_inv;
                    if one_success == 0 {
                        self.bns.edge[pe].flow += 1;
                        self.bns.vert[v1].st_edge.flow += 1;
                        self.bns.vert[v2].st_edge.flow += 1;
                        self.bns.tot_st_flow += 2;
                    }
                }
            }
        }
        self.curr.num_edges = 0;
        if cur_success != 0 {
            return self.recalc_after_rollback() || !self.c2i.b_has_difference;
        }
        false
    }
}

/// `FixMobileHRestoredStructure`: processes the extra or missing Mobile-H
/// attachments. Returns the last BNS or comparison result (negative on an
/// error), as the C function.
#[allow(clippy::too_many_arguments)]
pub fn fix_mobile_h_restored_structure(
    ip: &InputParms,
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    at3: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &mut AllTcGroups,
    p_inchi: [Option<&InChI>; 2],
    forbidden_edge_mask: i32,
    forbidden_stereo_edge_mask: i32,
) -> i32 {
    let mut c = FixMhCtx {
        ip,
        bns,
        bd,
        ps,
        at2,
        at3,
        p_va,
        tcg,
        p_inchi,
        mask: forbidden_edge_mask,
        mask_inv: !(forbidden_edge_mask as i8),
        c2i: Cmp2MhInchi::default(),
        all_charge: EdgeList::default(),
        curr: EdgeList::default(),
        curr2: EdgeList::default(),
        taut: EdgeList::default(),
        ret: 0,
        res: BnsTestResult::default(),
    };
    if !c.prologue(forbidden_stereo_edge_mask) {
        let _ =
            c.case_01() || c.case_02() || c.case_03() || c.case_05() || c.case_06() || c.case_07();
    }
    c.ret
}
