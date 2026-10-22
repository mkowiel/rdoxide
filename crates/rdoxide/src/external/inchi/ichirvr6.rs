//! The InChI -> structure restoration, part 6 (`INCHI_BASE/src/ichirvr6.c`):
//! `FixRestoredStructureStereo`.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::{InChI, InChIStereo};
use super::ichi_bns::{BnData, BnStruct, NO_VERTEX};
use super::ichimake::{
    IDIF_CHARGE, IDIF_DIFF_TG_ENDP, IDIF_EXTRA_TG_ENDP, IDIF_ISO_AT, IDIF_LESS_FH, IDIF_LESS_H,
    IDIF_MISS_TG_ENDP, IDIF_MORE_FH, IDIF_MORE_H, IDIF_MULTIPLE_TG, IDIF_NO_TAUT, IDIF_NUM_ISO_AT,
    IDIF_NUM_TG, IDIF_POSITION_H, IDIF_PROBLEM, IDIF_REM_ISO_H, IDIF_REM_PROT, IDIF_SB_EXTRA,
    IDIF_SB_EXTRA_UNDF, IDIF_SB_MISS, IDIF_SB_MISS_UNDF, IDIF_SB_PARITY, IDIF_SC_EXTRA,
    IDIF_SC_EXTRA_UNDF, IDIF_SC_INV, IDIF_SC_MISS, IDIF_SC_MISS_UNDF, IDIF_SC_PARITY,
    IDIF_SINGLE_TG, IDIF_TG, IDIF_WRONG_TAUT, Icr, compare_icr, compare_reversed_inchi2,
};
use super::ichiparm::InputParms;
use super::ichiread::{RI_ERR_SYNTAX, TAUT_NON, TAUT_YES};
use super::ichirvr1::{
    BnsTestResult, EdgeList, add_to_edge_list, get_charge_flower_upper_edge,
    make_one_inchi_out_of_str_from_inchi2, remove_forbidden_edge_mask, run_bns_restore_once,
    run_bns_test_once, set_forbidden_edge_mask,
};
use super::ichirvr2::StrFromInChI;
use super::ichirvr4::fill_out_extra_fixed_h_data_restr;
use super::ichirvr5::get_plus_minus_vertex;
use super::ichirvrs::{AllTcGroups, ValAt};
use super::inpdef::InpAtom;

const INC_ADD_EDGE: i32 = 64;
const RI_ERR_PROGR: i32 = super::ichiread::RI_ERR_PROGR;

const IDIFF_SB: i32 =
    IDIF_SB_PARITY | IDIF_SB_EXTRA_UNDF | IDIF_SB_EXTRA | IDIF_SB_MISS_UNDF | IDIF_SB_MISS;
const IDIFF_CONSTIT: i32 = IDIF_POSITION_H
    | IDIF_MORE_FH
    | IDIF_LESS_FH
    | IDIF_MORE_H
    | IDIF_LESS_H
    | IDIF_NO_TAUT
    | IDIF_WRONG_TAUT
    | IDIF_SINGLE_TG
    | IDIF_MULTIPLE_TG
    | IDIF_NUM_TG
    | IDIF_EXTRA_TG_ENDP
    | IDIF_MISS_TG_ENDP
    | IDIF_TG
    | IDIF_NUM_ISO_AT
    | IDIF_ISO_AT
    | IDIF_CHARGE
    | IDIF_REM_PROT
    | IDIF_REM_ISO_H
    | IDIF_DIFF_TG_ENDP;
const IDIFF_STEREO: i32 = IDIF_SC_INV
    | IDIF_SC_PARITY
    | IDIF_SC_EXTRA_UNDF
    | IDIF_SC_EXTRA
    | IDIF_SC_MISS_UNDF
    | IDIF_SC_MISS
    | IDIF_SB_PARITY
    | IDIF_SB_EXTRA_UNDF
    | IDIF_SB_EXTRA
    | IDIF_SB_MISS_UNDF
    | IDIF_SB_MISS;

/// the stereo layer that is being compared: isotopic when it has
/// stereo, else the basic one
fn pick_stereo(c: &InChI) -> Option<&InChIStereo> {
    match &c.stereo_isotopic {
        Some(s) if s.n_number_of_stereo_bonds + s.n_number_of_stereo_centers != 0 => Some(s),
        _ => c.stereo.as_ref(),
    }
}

/// the Mobile-H layer in case of Fixed-H
fn pick_stereo2(c: Option<&InChI>, b_mobile_h: i32) -> Option<&InChIStereo> {
    let c = c?;
    if b_mobile_h == TAUT_YES as i32 || c.n_number_of_atoms == 0 || c.b_deleted != 0 {
        None
    } else {
        pick_stereo(c)
    }
}

/// atoms (0-based) of the stereo bond `j`
fn sb(st: &Option<InChIStereo>, j: usize) -> Option<(usize, usize)> {
    let s = st.as_ref()?;
    let a = i32::from(*s.n_bond_atom1.get(j)?) - 1;
    let b = i32::from(*s.n_bond_atom2.get(j)?) - 1;
    (a >= 0 && b >= 0).then_some((a as usize, b as usize))
}

struct Ctx<'a> {
    ip: &'a InputParms,
    bns: &'a mut BnStruct,
    bd: &'a mut BnData,
    ps: &'a mut StrFromInChI,
    at2: &'a mut Vec<InpAtom>,
    at3: &'a mut Vec<InpAtom>,
    p_va: &'a [ValAt],
    tcg: &'a mut AllTcGroups,
    p_inchi: [Option<&'a InChI>; 2],
    icr: &'a mut Icr,
    icr2: &'a mut Icr,
    num_run_bns: &'a mut i32,
    mask: i32,
    stereo_mask: i32,
    cmp1: i32,
    cmp2: i32,
    i_orig: usize,
    i_revr: usize,
    st_inchi: Option<&'a InChIStereo>,
    st2_inchi: Option<&'a InChIStereo>,
    st_revrs: Option<InChIStereo>,
    st2_revrs: Option<InChIStereo>,
    all_charge: EdgeList,
    curr: EdgeList,
    all_rad: EdgeList,
    fixed_stereo: EdgeList,
    taut: [EdgeList; 2],
    ret: i32,
    res: BnsTestResult,
}

impl Ctx<'_> {
    fn refresh_revrs(&mut self) {
        self.st_revrs = self.ps.p_one_inchi[0]
            .as_ref()
            .and_then(pick_stereo)
            .cloned();
        self.st2_revrs = pick_stereo2(self.ps.p_one_inchi[1].as_ref(), self.ps.b_mobile_h).cloned();
    }

    fn vert_edge(&self, v: usize, j: usize) -> usize {
        self.bns.iedge[self.bns.vert[v].iedge + j] as usize
    }

    /// the edge between `v1` and `v2`
    fn find_edge(&self, v1: usize, v2: usize) -> Option<usize> {
        (0..self.at2[v1].valence as usize)
            .map(|k| self.vert_edge(v1, k))
            .find(|&e| v2 == (usize::from(self.bns.edge[e].neighbor12) ^ v1))
    }

    fn minus_edge(&self, j: usize) -> i32 {
        self.p_va[j].n_c_minus_group_edge - 1
    }

    fn plus_edge(&self, j: usize) -> i32 {
        self.p_va[j].n_c_plus_group_edge - 1
    }

    fn free_edge(&self, k: i32) -> bool {
        k >= 0 && self.bns.edge[k as usize].forbidden == 0
    }

    fn is_endpoint(&self, j: usize) -> bool {
        self.at2[j].endpoint != 0 || self.ps.endpoint.as_ref().is_some_and(|e| e[j] != 0)
    }

    /// appends to a list (0 = all charges, 1/2 = `TautMinusEdges[0/1]`);
    /// true on an error
    fn add(&mut self, which: usize, e: i32) -> bool {
        let list = match which {
            0 => &mut self.all_charge,
            1 => &mut self.taut[0],
            _ => &mut self.taut[1],
        };
        self.ret = add_to_edge_list(list, e, INC_ADD_EDGE);
        self.ret != 0
    }

    /// the upper edge of the N charge flower when it has the flow 1 and is
    /// not forbidden (C tests `edge[j].forbidden`, with `j` the atom)
    fn add_flower(&mut self, j: usize, k: i32) -> bool {
        let va = self.p_va[j];
        if va.c_num_valence_electrons == 5 && va.c_metal == 0 {
            let k = get_charge_flower_upper_edge(self.bns, k);
            if k != NO_VERTEX
                && self.bns.edge[j].forbidden == 0
                && self.bns.edge[k as usize].flow != 0
                && self.add(0, k)
            {
                return true;
            }
        }
        false
    }

    /// -NHn terminal N atoms
    fn is_nh(&self, v: usize, num_h: i8) -> bool {
        let a = &self.at2[v];
        let va = &self.p_va[v];
        a.valence == 1
            && a.num_h == num_h
            && a.endpoint == 0
            && va.c_num_valence_electrons == 5
            && va.c_periodic_row_number == 1
    }

    fn compare_main(&mut self, into_icr2: bool) -> i32 {
        let icr = if into_icr2 {
            &mut *self.icr2
        } else {
            &mut *self.icr
        };
        compare_reversed_inchi2(
            self.ps.p_one_inchi[0].as_ref(),
            self.p_inchi[0],
            self.ps.p_one_inchi_aux[0].as_ref(),
            None,
            icr,
        ) as i32
    }

    /// compares the restored InChI with the original one (both layers)
    fn compare_both(&mut self, case06: bool) -> bool {
        self.cmp1 = self.compare_main(case06);
        if self.cmp1 & IDIF_PROBLEM != 0 {
            self.ret = RI_ERR_PROGR;
            return true;
        }
        if !case06 {
            self.cmp2 = 0;
        }
        *self.icr2 = Icr::default();
        if self.i_revr != 0 || self.i_orig != 0 {
            self.cmp2 = compare_reversed_inchi2(
                self.ps.p_one_inchi[self.i_revr].as_ref(),
                self.p_inchi[self.i_orig],
                self.ps.p_one_inchi_aux[self.i_revr].as_ref(),
                None,
                &mut *self.icr2,
            ) as i32;
            // C checks `cmpInChI` here
            if self.cmp1 & IDIF_PROBLEM != 0 {
                self.ret = RI_ERR_PROGR;
                return true;
            }
        }
        false
    }

    /// recalculates the InChI of the structure and compares it with the
    /// original one; true = exit the function (`self.ret` is set).
    /// `case06` reproduces the C code of case 06: the first comparison goes
    /// to `icr2` and `cmpInChI2` is not reset.
    fn recompute(&mut self, case06: bool) -> bool {
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
        if self.compare_both(case06) {
            return true;
        }
        self.refresh_revrs();
        false
    }

    /// decrements the order of the stereobond `pe`, fixes it and checks
    /// whether the BNS can move the charge; restores on success.
    /// Returns true when `cur_success` is incremented.
    fn try_move(&mut self, pe: usize) -> bool {
        let mask = self.mask as i8;
        let v1 = usize::from(self.bns.edge[pe].neighbor1);
        let v2 = usize::from(self.bns.edge[pe].neighbor12) ^ v1;
        self.bns.edge[pe].forbidden |= mask;
        self.bns.edge[pe].flow -= 1;
        self.bns.vert[v1].st_edge.flow -= 1;
        self.bns.vert[v2].st_edge.flow -= 1;
        self.bns.tot_st_flow -= 2;
        self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
        let (a, b) = (v1 as i32, v2 as i32);
        if self.ret == 1
            && ((self.res.v_last == a && self.res.v_first == b)
                || (self.res.v_last == b && self.res.v_first == a))
            && self.res.delta_charge == 0
        {
            self.ret = run_bns_restore_once(self.bns, self.bd);
            if self.ret > 0 {
                *self.num_run_bns += 1;
                return true;
            }
        } else {
            self.bns.edge[pe].forbidden &= !mask;
            self.bns.edge[pe].flow += 1;
            self.bns.vert[v1].st_edge.flow += 1;
            self.bns.vert[v2].st_edge.flow += 1;
            self.bns.tot_st_flow += 2;
        }
        false
    }

    /// case 01: extra stereogenic bond, radical present; true = exit
    fn case_01(&mut self) -> bool {
        let max_success = self.bns.tot_st_cap - self.bns.tot_st_flow;
        let mut cur_success = 0;
        remove_forbidden_edge_mask(self.bns, &self.fixed_stereo, self.stereo_mask);
        let mut i = 0;
        while i < self.icr.sb_in2_only.len() && cur_success < max_success {
            let j12 = self.icr.sb_in2_only[i] as usize;
            i += 1;
            let st = self.st_inchi.cloned();
            let Some((v1, v2)) = sb(&st, j12) else {
                continue;
            };
            let Some(pe) = self.find_edge(v1, v2) else {
                self.ret = RI_ERR_SYNTAX;
                return true;
            };
            // check v1
            self.bns.vert[v1].st_edge.cap -= 1;
            self.bns.vert[v1].st_edge.flow -= 1;
            self.bns.vert[v2].st_edge.flow -= 1;
            self.bns.edge[pe].flow -= 1; // new radical on v2
            self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
            self.bns.vert[v1].st_edge.cap += 1;
            self.bns.vert[v1].st_edge.flow += 1;
            self.bns.vert[v2].st_edge.flow += 1;
            self.bns.edge[pe].flow += 1;
            let (a, b) = (v1 as i32, v2 as i32);
            let mut v_rad = NO_VERTEX;
            if self.ret == 1
                && self.res.delta_charge == 0
                && (b == self.res.v_first || b == self.res.v_last)
            {
                v_rad = if b == self.res.v_first {
                    self.res.v_last
                } else {
                    self.res.v_first
                };
            } else {
                self.bns.vert[v2].st_edge.cap -= 1;
                self.bns.vert[v2].st_edge.flow -= 1;
                self.bns.vert[v1].st_edge.flow -= 1;
                self.bns.edge[pe].flow -= 1; // new radical on v1
                self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                self.bns.vert[v2].st_edge.cap += 1;
                self.bns.vert[v2].st_edge.flow += 1;
                self.bns.vert[v1].st_edge.flow += 1;
                self.bns.edge[pe].flow += 1;
                if self.ret == 1
                    && self.res.delta_charge == 0
                    && (a == self.res.v_first || a == self.res.v_last)
                {
                    v_rad = if a == self.res.v_first {
                        self.res.v_last
                    } else {
                        self.res.v_first
                    };
                }
            }
            if v_rad == NO_VERTEX {
                continue; // radical did not affect this bond
            }
            let v_rad = v_rad as usize;
            // detect =N-O*
            let mut is_no = false;
            let mut a_n = 0usize;
            if self.p_va[v_rad].c_num_valence_electrons == 6 && self.at2[v_rad].valence == 1 {
                let pe_rad = self.vert_edge(v_rad, 0);
                if self.bns.edge[pe_rad].flow == 0 {
                    a_n = usize::from(self.bns.edge[pe_rad].neighbor12) ^ v_rad;
                    is_no =
                        self.p_va[a_n].c_num_valence_electrons == 5 && self.at2[a_n].valence == 2;
                }
            }
            if is_no {
                let pe_no = self.vert_edge(v_rad, 0);
                // N-O*  => N=O
                self.bns.edge[pe_no].flow += 1;
                self.bns.vert[v_rad].st_edge.flow += 1;
                self.bns.vert[a_n].st_edge.cap += 1;
                self.bns.vert[a_n].st_edge.flow += 1;
                self.bns.tot_st_cap += 1;
                self.bns.tot_st_flow += 2;
                cur_success += 1;
            } else {
                // all other radicals that affect stereo
                let delta = self.bns.vert[v_rad].st_edge.cap - self.bns.vert[v_rad].st_edge.flow;
                self.bns.vert[v_rad].st_edge.cap -= delta;
                self.bns.tot_st_cap -= delta;
            }
        }
        set_forbidden_edge_mask(self.bns, &self.fixed_stereo, self.stereo_mask);
        cur_success != 0 && self.recompute(false)
    }

    /// case 02: missing stereogenic bond in Mobile-H only
    fn case_02(&mut self) -> bool {
        let icr_prev = self.icr.clone();
        let icr2_prev = self.icr2.clone();
        self.all_rad.num_edges = 0;
        for i in 0..self.ps.num_atoms {
            let st = self.bns.vert[i as usize].st_edge;
            if st.cap - st.flow == 1 {
                self.ret = add_to_edge_list(&mut self.all_rad, i, INC_ADD_EDGE);
                if self.ret != 0 {
                    return true;
                }
            }
        }
        for i in 0..self.all_rad.num_edges as usize {
            let j = self.all_rad.pn_edges[i] as usize;
            self.bns.vert[j].st_edge.cap -= 1;
            self.bns.tot_st_cap -= 1;
        }
        // re-create InChI and see whether it looks better
        if self.recompute(false) {
            return true;
        }
        let ret_c = compare_icr(self.icr, &icr_prev, IDIFF_CONSTIT as u64).0;
        let ret_s = compare_icr(self.icr, &icr_prev, IDIFF_STEREO as u64).0;
        let ret2_c = compare_icr(self.icr2, &icr2_prev, IDIFF_CONSTIT as u64).0;
        let ret2_s = compare_icr(self.icr2, &icr2_prev, IDIFF_STEREO as u64).0;
        if !(ret_c <= 0 && ret_s <= 0 && ret2_c <= 0 && ret2_s < 0) {
            // reject
            for i in 0..self.all_rad.num_edges as usize {
                let j = self.all_rad.pn_edges[i] as usize;
                self.bns.vert[j].st_edge.cap += 1;
                self.bns.tot_st_cap += 1;
            }
            // re-create InChI: return to the previous state
            if self.recompute(false) {
                return true;
            }
        }
        false
    }

    /// the loop over the tautomeric charges of cases 03 and 04
    fn fix_taut_charges(&mut self, pe: usize) -> bool {
        let mask = self.mask;
        set_forbidden_edge_mask(self.bns, &self.all_charge, mask);
        let mut done = false;
        let mut k = 1;
        while k < 4 && !done {
            let mut num_endpoints = self.taut[0].num_edges + self.taut[1].num_edges;
            if k == 2 {
                set_forbidden_edge_mask(self.bns, &self.taut[0], mask);
                num_endpoints -= self.taut[0].num_edges;
            }
            if k == 1 {
                set_forbidden_edge_mask(self.bns, &self.taut[1], mask);
                num_endpoints -= self.taut[1].num_edges;
            }
            if num_endpoints >= 2 && self.try_move(pe) {
                done = true;
            }
            if k == 2 {
                remove_forbidden_edge_mask(self.bns, &self.taut[0], mask);
            }
            if k == 1 {
                remove_forbidden_edge_mask(self.bns, &self.taut[1], mask);
            }
            k += 1;
        }
        remove_forbidden_edge_mask(self.bns, &self.all_charge, mask);
        done
    }

    /// case 03: extra stereogenic bond in Fixed-H only
    fn case_03(&mut self) -> bool {
        let mut cur_success = 0;
        self.taut[0].num_edges = 0;
        self.taut[1].num_edges = 0;
        self.all_charge.num_edges = 0;
        for i in 0..self.icr.sb_undef_in1_only.len() {
            let j12 = self.icr.sb_undef_in1_only[i] as usize;
            let Some((v1, v2)) = sb(&self.st_revrs, j12) else {
                continue;
            };
            if let Some(s2) = &self.st2_revrs {
                // reject if it is extra in Mobile-H also
                let n = self.icr2.sb_undef_in1_only.len();
                if n != 0 {
                    let mut j = 0;
                    while j < n {
                        let k = self.icr2.sb_undef_in1_only[j] as usize;
                        if v1 == usize::from(s2.n_bond_atom1[k])
                            && v2 == usize::from(s2.n_bond_atom2[k])
                        {
                            break;
                        }
                        j += 1;
                    }
                    if j < self.icr.sb_in1_only.len() {
                        continue; // extra stereobond in Mobile H also
                    }
                }
            }
            // reject if it is a stereobond in Mobile-H also
            if let Some(s2) = self.st2_inchi {
                let n = s2.n_number_of_stereo_bonds as usize;
                if n != 0 {
                    let j = (0..n).find(|&j| {
                        v1 == usize::from(s2.n_bond_atom1[j])
                            && v2 == usize::from(s2.n_bond_atom1[j])
                    });
                    if j.is_some() {
                        continue;
                    }
                }
            }
            let Some(pe) = self.find_edge(v1, v2) else {
                self.ret = RI_ERR_SYNTAX;
                return true;
            };
            // fix all charges except negative charges on tautomeric endpoints
            if self.all_charge.num_edges == 0
                && self.taut[0].num_edges == 0
                && self.taut[1].num_edges == 0
            {
                for j in 0..self.ps.num_atoms as usize {
                    let k = self.minus_edge(j);
                    if self.free_edge(k) {
                        let w = if self.ps.endpoint.as_ref().is_some_and(|e| e[j] == 0) {
                            0
                        } else if self.p_va[j].c_num_valence_electrons == 6 {
                            1
                        } else {
                            2
                        };
                        if self.add(w, k) {
                            return true;
                        }
                    }
                    let k = self.plus_edge(j);
                    if self.free_edge(k) && (self.add(0, k) || self.add_flower(j, k)) {
                        return true;
                    }
                }
            }
            if self.bns.edge[pe].flow == 0 {
                continue;
            }
            if self.fix_taut_charges(pe) {
                cur_success += 1;
            }
        }
        cur_success != 0 && self.recompute(false)
    }

    /// case 04: extra stereogenic bond whose parity is unknown
    fn case_04(&mut self) -> bool {
        let mut cur_success = 0;
        self.taut[0].num_edges = 0;
        self.taut[1].num_edges = 0;
        self.all_charge.num_edges = 0;
        for i in 0..self.icr.sb_undef_in1_only.len() {
            let j12 = self.icr.sb_undef_in1_only[i] as usize;
            let Some((v1, v2)) = sb(&self.st_revrs, j12) else {
                continue;
            };
            let Some(pe) = self.find_edge(v1, v2) else {
                self.ret = RI_ERR_SYNTAX;
                return true;
            };
            if let Some(s2) = &self.st2_revrs {
                // reject if it is not extra in Mobile-H also
                let n = self.icr2.sb_undef_in1_only.len();
                if n != 0 {
                    let mut j = 0;
                    while j < n {
                        let k = self.icr2.sb_undef_in1_only[j] as usize;
                        if v1 == usize::from(s2.n_bond_atom1[k])
                            && v2 == usize::from(s2.n_bond_atom2[k])
                        {
                            break;
                        }
                        j += 1;
                    }
                    if j == self.icr.sb_in1_only.len() {
                        continue; // extra stereobond only in Fixed-H
                    }
                }
            }
            if self.all_charge.num_edges == 0
                && self.taut[0].num_edges == 0
                && self.taut[1].num_edges == 0
            {
                for j in 0..self.ps.num_atoms as usize {
                    let k = self.minus_edge(j);
                    if self.free_edge(k) && self.add(0, k) {
                        return true;
                    }
                    let k = self.plus_edge(j);
                    if self.free_edge(k) {
                        let may_be_unfixed = self.at2[j].num_h == 0
                            && !self.ps.endpoint.as_ref().is_some_and(|e| e[j] != 0);
                        let va = self.p_va[j];
                        let w = if (may_be_unfixed && va.c_num_valence_electrons == 6)
                            || (va.c_num_valence_electrons == 5 && va.c_periodic_row_number > 1)
                        {
                            1 // O & P
                        } else if may_be_unfixed
                            && va.c_num_valence_electrons == 5
                            && va.c_periodic_row_number == 1
                        {
                            2 // N
                        } else {
                            0
                        };
                        if self.add(w, k) || self.add_flower(j, k) {
                            return true;
                        }
                    }
                }
            }
            if self.bns.edge[pe].flow == 0 {
                continue;
            }
            if self.fix_taut_charges(pe) {
                cur_success += 1;
            }
        }
        cur_success != 0 && self.recompute(false)
    }

    /// the charge edges of cases 05 and 06
    fn build_56(&mut self) -> bool {
        for j in 0..self.ps.num_atoms as usize {
            let k = self.minus_edge(j);
            if self.free_edge(k) && self.add(0, k) {
                return true;
            }
            let k = self.plus_edge(j);
            if self.free_edge(k) {
                let va = self.p_va[j];
                let a = &self.at2[j];
                // do not fix -N< or =N(+)<
                let skip = va.c_num_valence_electrons == 5
                    && va.c_periodic_row_number == 1
                    && a.num_h == 0
                    && a.valence == 3
                    && !self.is_endpoint(j);
                if (!skip && self.add(1, k)) || self.add_flower(j, k) {
                    return true;
                }
            }
        }
        false
    }

    /// case 05: extra stereogenic bond on =NH2(+), Mobile-H
    fn case_05(&mut self) -> bool {
        let mask = self.mask;
        let mut cur_success = 0;
        self.all_charge.num_edges = 0;
        for i in 0..self.icr.sb_undef_in1_only.len() {
            let j12 = self.icr.sb_undef_in1_only[i] as usize;
            let Some((v1, v2)) = sb(&self.st_revrs, j12) else {
                continue;
            };
            let (i1, i2) = (self.is_nh(v1, 1), self.is_nh(v2, 1));
            if i1 == i2 {
                continue;
            }
            let Some(pe) = self.find_edge(v1, v2) else {
                self.ret = RI_ERR_SYNTAX;
                return true;
            };
            if self.bns.edge[pe].flow != 1 {
                continue; // already charged
            }
            let a_n = if i1 { v1 } else { v2 };
            let a_c = if i1 { v2 } else { v1 };
            let v_pm = get_plus_minus_vertex(self.bns, self.tcg, true, false);
            if v_pm == NO_VERTEX {
                break; // cannot do anything
            }
            let v_pm = v_pm as usize;
            // increase edges to -Y-(+/-)-Y- capacities
            for i1 in 0..self.bns.vert[v_pm].num_adj_edges as usize {
                let e1 = self.vert_edge(v_pm, i1);
                let i2 = usize::from(self.bns.edge[e1].neighbor12) ^ v_pm;
                for k in 0..self.bns.vert[i2].num_adj_edges as usize {
                    let e2 = self.vert_edge(i2, k);
                    self.bns.edge[e2].cap += 1;
                }
            }
            // fix all charges except (+) on -N<
            if self.all_charge.num_edges == 0 && self.build_56() {
                return true;
            }
            // make bond to =NH single, add radical to aC
            self.bns.edge[pe].flow -= 1;
            self.bns.vert[a_n].st_edge.flow -= 1;
            self.bns.vert[a_n].st_edge.cap -= 1;
            self.bns.vert[a_c].st_edge.flow -= 1;
            self.bns.vert[v_pm].st_edge.cap += 1;
            self.bns.tot_st_flow -= 2;
            // fix C-NH bond
            if self.add(0, pe as i32) {
                return true;
            }
            set_forbidden_edge_mask(self.bns, &self.all_charge, mask);
            self.bns.edge[pe].cap += 1;
            self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
            let (pm, c) = (v_pm as i32, a_c as i32);
            if self.ret == 1
                && ((self.res.v_last == pm && self.res.v_first == c)
                    || (self.res.v_last == c && self.res.v_first == pm))
                && self.res.delta_charge == 1
            {
                self.ret = run_bns_restore_once(self.bns, self.bd);
                if self.ret > 0 {
                    *self.num_run_bns += 1;
                    self.ps.at[a_n].num_h += 1;
                    self.tcg.total_charge += 1;
                    cur_success += 1;
                }
            } else {
                self.bns.edge[pe].flow += 1;
                self.bns.vert[a_n].st_edge.flow += 1;
                self.bns.vert[a_n].st_edge.cap += 1;
                self.bns.vert[a_c].st_edge.flow += 1;
                self.bns.vert[v_pm].st_edge.cap -= 1;
                self.bns.tot_st_flow += 2;
                remove_forbidden_edge_mask(self.bns, &self.all_charge, mask);
                self.all_charge.num_edges -= 1;
                self.curr.num_edges = 0;
                continue; // should not happen
            }
            remove_forbidden_edge_mask(self.bns, &self.all_charge, mask);
            self.all_charge.num_edges -= 1;
            self.curr.num_edges = 0;
        }
        cur_success != 0 && self.recompute(false)
    }

    /// case 06: extra stereogenic bond on =NH2(+), Fixed-H
    fn case_06(&mut self) -> bool {
        let mask = self.mask;
        let mut cur_success = 0;
        self.all_charge.num_edges = 0;
        for i in 0..self.icr2.sb_undef_in1_only.len() {
            let j12 = self.icr2.sb_undef_in1_only[i] as usize;
            let Some((v1, v2)) = sb(&self.st2_revrs, j12) else {
                continue;
            };
            let (i1, i2) = (self.is_nh(v1, 2), self.is_nh(v2, 2));
            if i1 == i2 {
                continue;
            }
            let Some(pe) = self.find_edge(v1, v2) else {
                self.ret = RI_ERR_SYNTAX;
                return true;
            };
            if self.bns.edge[pe].flow != 1 {
                continue; // already charged
            }
            let a_n = if i1 { v1 } else { v2 };
            let e_plus = self.plus_edge(a_n);
            // must be (+) charged
            if e_plus < 0
                || self.bns.edge[e_plus as usize].flow != 0
                || self.bns.edge[e_plus as usize].forbidden != 0
            {
                continue;
            }
            // move (+) from =NH2(+) to some other -N<
            if self.all_charge.num_edges == 0 && self.build_56() {
                return true;
            }
            // the pePlus edge is already fixed; unfix it
            if self.bns.edge[pe].flow == 0 {
                continue;
            }
            self.bns.edge[e_plus as usize].forbidden &= !(mask as i8);
            if self.try_move(pe) {
                cur_success += 1;
            }
            self.bns.edge[pe].forbidden &= !(mask as i8);
            remove_forbidden_edge_mask(self.bns, &self.all_charge, mask);
        }
        cur_success != 0 && self.recompute(true)
    }

    fn run(&mut self) {
        // the first comparison
        if self.compare_both(false) {
            return;
        }
        if self.cmp1 & IDIFF_SB == 0 && self.cmp2 & IDIFF_SB == 0 {
            return;
        }
        // need to temporarily remove fixing of stereogenic bonds
        let stereo_mask = self.stereo_mask as i8;
        for i in 0..self.ps.num_atoms as usize {
            for j in 0..self.at2[i].valence as usize {
                let e = self.vert_edge(i, j);
                let pe = &self.bns.edge[e];
                if j == usize::from(pe.neighbor1) && pe.forbidden & stereo_mask != 0 {
                    self.ret = add_to_edge_list(&mut self.fixed_stereo, e as i32, INC_ADD_EDGE);
                    if self.ret != 0 {
                        return;
                    }
                }
            }
        }
        let free_cap = |c: &Self| c.bns.tot_st_cap - c.bns.tot_st_flow;
        if self.cmp1 & IDIF_SB_MISS != 0
            && (self.cmp2 == 0 || self.cmp2 & IDIF_SB_MISS != 0)
            && free_cap(self) > 0
            && self.case_01()
        {
            return;
        }
        if self.cmp1 & IDIF_SB_MISS == 0
            && self.cmp2 & IDIF_SB_MISS != 0
            && !self.icr2.sb_in2_only.is_empty()
            && free_cap(self) > 0
            && self.case_02()
        {
            return;
        }
        if self.ps.b_mobile_h == TAUT_NON as i32
            && self.cmp1 & IDIF_SB_EXTRA_UNDF != 0
            && self.ps.endpoint.is_some()
            && self.case_03()
        {
            return;
        }
        if self.cmp1 & IDIF_SB_EXTRA_UNDF != 0 && self.case_04() {
            return;
        }
        if self.ps.b_mobile_h == TAUT_YES as i32
            && self.cmp1 & IDIF_SB_EXTRA_UNDF != 0
            && self.ps.ti.num_t_groups == 0
            && self.case_05()
        {
            return;
        }
        if self.ps.b_mobile_h == TAUT_NON as i32
            && self.st2_revrs.is_some()
            && self.cmp2 & IDIF_SB_EXTRA_UNDF != 0
            && self.ps.ti.num_t_groups == 0
        {
            self.case_06();
        }
    }
}

/// `FixRestoredStructureStereo`: processes the extra or missing stereobonds
/// of the restored structure `at2` (the Fixed-H structure when the structure
/// is Fixed-H); the results of the comparison are left in `icr` and `icr2`.
/// Returns the last BNS or comparison result (negative on an error).
#[allow(clippy::too_many_arguments)]
pub fn fix_restored_structure_stereo(
    ip: &InputParms,
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    at3: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &mut AllTcGroups,
    p_inchi: [Option<&InChI>; 2],
    icr: &mut Icr,
    icr2: &mut Icr,
    num_run_bns: &mut i32,
    forbidden_edge_mask: i32,
    forbidden_stereo_edge_mask: i32,
) -> i32 {
    let Some(inchi0) = p_inchi[0] else {
        return RI_ERR_PROGR;
    };
    let st_inchi = pick_stereo(inchi0);
    let st2_inchi = pick_stereo2(p_inchi[1], ps.b_mobile_h);
    let exists = |c: &Option<InChI>| {
        c.as_ref()
            .is_some_and(|c| c.n_number_of_atoms != 0 && c.b_deleted == 0)
    };
    let (i_orig, i_revr) = if ps.b_mobile_h == TAUT_NON as i32 {
        (
            usize::from(p_inchi[1].is_some_and(|c| c.n_number_of_atoms != 0 && c.b_deleted == 0)),
            usize::from(exists(&ps.p_one_inchi[1])),
        )
    } else {
        (0, 0)
    };
    let mut c = Ctx {
        ip,
        bns,
        bd,
        ps,
        at2,
        at3,
        p_va,
        tcg,
        p_inchi,
        icr,
        icr2,
        num_run_bns,
        mask: forbidden_edge_mask,
        stereo_mask: forbidden_stereo_edge_mask,
        cmp1: 0,
        cmp2: 0,
        i_orig,
        i_revr,
        st_inchi,
        st2_inchi,
        st_revrs: None,
        st2_revrs: None,
        all_charge: EdgeList::default(),
        curr: EdgeList::default(),
        all_rad: EdgeList::default(),
        fixed_stereo: EdgeList::default(),
        taut: [EdgeList::default(), EdgeList::default()],
        ret: 0,
        res: BnsTestResult::default(),
    };
    c.refresh_revrs();
    c.run();
    set_forbidden_edge_mask(c.bns, &c.fixed_stereo, forbidden_stereo_edge_mask);
    c.ret
}
