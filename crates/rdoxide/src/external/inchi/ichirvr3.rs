//! The InChI -> structure restoration, part 3 (`INCHI_BASE/src/ichirvr3.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::InChI;
use super::ichi_bns::{BnData, BnStruct, NO_VERTEX};
use super::ichiparm::InputParms;
use super::ichirvr1::{
    BnsTestResult, EdgeList, add_to_edge_list, alloc_edge_list, find_in_edge_list,
    get_charge_flower_upper_edge, make_one_inchi_out_of_str_from_inchi2,
    remove_forbidden_edge_mask, remove_from_edge_list_by_value, run_bns_restore_once,
    run_bns_test_once, set_forbidden_edge_mask,
};
use super::ichirvr2::StrFromInChI;
use super::ichirvr4::{
    Cmp2FhInchi, MAX_DIFF_FIXH, fill_out_cmp2fh_inchi, fill_out_extra_fixed_h_data_restr,
};
use super::ichirvrs::{AllTcGroups, ValAt};
use super::ichitaut::TGroupInfo;
use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE};
use super::inpdef::InpAtom;

const INC_ADD_EDGE: i32 = 64;

/// the edge lists of the context
#[derive(Debug, Clone, Copy)]
pub enum L {
    AllCharge,
    Curr,
    SFlower,
    NFlower,
    OtherNFlower,
    FixedLargeRingStereo,
    AllBond,
}
/// `EL_NUMBER_N`
const EL_NUMBER_N: u8 = 7;

/// `bHas_N_V`: the number of neutral N(V) atoms with three neighbours
pub fn b_has_n_v(at2: &[InpAtom], num_atoms: usize) -> i32 {
    at2[..num_atoms]
        .iter()
        .filter(|a| {
            a.el_number == EL_NUMBER_N
                && a.charge == 0
                && a.num_h == 0
                && a.radical == 0
                && a.chem_bonds_valence == 5
                && a.valence == 3
        })
        .count() as i32
}

/// types of `TgDiffHChgFH`
pub const F_NUM_R_POS_CHG_H: usize = 0;
pub const F_NUM_R_POS_CHG_U: usize = 1;
pub const F_NUM_R_NEG_CHG_O: usize = 2;
pub const F_NUM_R_NEG_CHG_N: usize = 3;
pub const F_NUM_R_NEUTRL_H: usize = 4;
pub const F_NUM_N_POS_CHG_H: usize = 5;
pub const F_NUM_N_POS_CHG_U: usize = 6;
pub const F_NUM_N_NEG_CHG_O: usize = 7;
pub const F_NUM_N_NEG_CHG_N: usize = 8;
pub const F_NUM_N_NEUTRL_H: usize = 9;
pub const F_NUM_ALL_CHG_T: usize = 10;

/// `TgDiffHChgFH`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TgDiffHChgFh {
    pub itg: i16,
    pub n_num_h_inchi: i16,
    pub n_num_h_revrs: i16,
    pub n_num_h_norml: i16,
    pub n_num_m_inchi: i16,
    pub n_num_m_revrs: i16,
    pub n_num_m_norml: i16,
    pub n_num_p_revrs: i16,
    pub n_num_p_norml: i16,
    pub n: [i16; F_NUM_ALL_CHG_T],
    pub i: [i16; F_NUM_ALL_CHG_T],
}

/// `FillTgDiffHChgFH`: describes the t-groups whose numbers of H or (-)
/// differ between the original InChI and the normalized restored structure;
/// the atoms of the groups are appended to `atom_ind_list`. Returns the
/// number of the t-groups found in `tdhc` or an error.
/// (Unlike C, the scan of the sorted index list stops at its end instead of
/// reading the memory behind it.)
pub fn fill_tg_diff_h_chg_fh(
    tdhc: &mut Vec<TgDiffHChgFh>,
    max_tdhc: usize,
    at2: &[InpAtom],
    atf: &[InpAtom],
    canon2atno_revrs: &[u16],
    p_va: &[ValAt],
    ti: &TGroupInfo,
    atom_ind_list: &mut EdgeList,
) -> i32 {
    let mut ind_list = EdgeList::default();
    atom_ind_list.num_edges = 0;
    tdhc.clear();
    let groups = ti.t_group.as_deref().unwrap_or(&[]);
    let mut ret;
    for itg in 0..ti.num_t_groups {
        let g = &groups[itg];
        let mut cur = TgDiffHChgFh {
            itg: itg as i16,
            n_num_h_inchi: (i32::from(g.num[0]) - i32::from(g.num[1])) as i16,
            n_num_m_inchi: g.num[1] as i16,
            ..Default::default()
        };
        let first = g.n_first_endpoint_at_no_pos as usize;
        let cur_len = ind_list.num_edges;
        for j in 0..g.n_num_endpoints as usize {
            let i = ti.n_endpoint_atom_number[first + j] as usize;
            let iat = canon2atno_revrs[i] as usize;
            let (a2, af) = (&at2[iat], &atf[iat]);
            let type_r: i32 = if a2.charge == 1 {
                cur.n_num_p_revrs += 1;
                (if a2.num_h != 0 {
                    F_NUM_R_POS_CHG_H
                } else {
                    F_NUM_R_POS_CHG_U
                }) as i32
            } else if a2.charge == -1 {
                cur.n_num_m_revrs += 1;
                match p_va[iat].c_num_valence_electrons {
                    6 => F_NUM_R_NEG_CHG_O as i32,
                    5 => F_NUM_R_NEG_CHG_N as i32,
                    _ => -1,
                }
            } else if a2.num_h != 0 && a2.valence == a2.chem_bonds_valence {
                F_NUM_R_NEUTRL_H as i32
            } else {
                -1
            };
            cur.n_num_h_revrs += i16::from(a2.num_h);
            let type_n: i32 = if af.charge == 1 {
                cur.n_num_p_norml += 1;
                (if af.num_h != 0 {
                    F_NUM_N_POS_CHG_H
                } else {
                    F_NUM_N_POS_CHG_U
                }) as i32
            } else if af.charge == -1 {
                cur.n_num_m_norml += 1;
                match p_va[iat].c_num_valence_electrons {
                    6 => F_NUM_N_NEG_CHG_O as i32,
                    5 => F_NUM_N_NEG_CHG_N as i32,
                    _ => -1,
                }
            } else if af.num_h != 0 && af.valence == af.chem_bonds_valence {
                F_NUM_N_NEUTRL_H as i32
            } else {
                -1
            };
            cur.n_num_h_norml += i16::from(af.num_h);
            if a2.charge < 0 || 0 < p_va[iat].n_c_plus_group_edge {
                for ty in [type_r, type_n] {
                    if ty >= 0 {
                        for v in [ty, itg as i32, iat as i32] {
                            ret = add_to_edge_list(&mut ind_list, v, INC_ADD_EDGE);
                            if ret != 0 {
                                return ret;
                            }
                        }
                    }
                }
            }
        }
        if cur.n_num_h_norml == cur.n_num_h_inchi && cur.n_num_m_norml == cur.n_num_m_inchi {
            ind_list.num_edges = cur_len;
            continue;
        }
        if tdhc.len() < max_tdhc {
            tdhc.push(cur);
        } else {
            ind_list.num_edges = cur_len;
            break;
        }
    }
    // fill out the atom index list
    if !tdhc.is_empty() {
        let list_len = ind_list.num_edges as usize;
        for ty in 0..F_NUM_ALL_CHG_T {
            let mut j = 0usize;
            for t in tdhc.iter_mut() {
                let mut num = 0;
                let itg = i32::from(t.itg);
                t.i[ty] = -999;
                while j + 1 < list_len && ind_list.pn_edges[j + 1] == itg {
                    if ind_list.pn_edges[j] == ty as i32 {
                        if num == 0 {
                            t.i[ty] = atom_ind_list.num_edges as i16;
                        }
                        num += 1;
                        ret =
                            add_to_edge_list(atom_ind_list, ind_list.pn_edges[j + 2], INC_ADD_EDGE);
                        if ret != 0 {
                            return ret;
                        }
                    }
                    j += 3;
                }
                t.n[ty] = num;
            }
        }
    }
    let _ = alloc_edge_list;
    tdhc.len() as i32
}

/// the state of `FixFixedHRestoredStructure` shared by its cases
#[allow(missing_debug_implementations)]
pub struct FixFhCtx<'a> {
    pub ip: &'a InputParms,
    pub bns: &'a mut BnStruct,
    pub bd: &'a mut BnData,
    pub ps: &'a mut StrFromInChI,
    pub at2: &'a mut Vec<InpAtom>,
    pub at3: &'a mut Vec<InpAtom>,
    pub p_va: &'a [ValAt],
    pub tcg: &'a AllTcGroups,
    pub p_inchi: [Option<&'a InChI>; 2],
    pub num_run_bns: i32,
    pub mask: i32,
    pub mask_inv: i8,
    pub c2i: Cmp2FhInchi,
    pub all_charge: EdgeList,
    pub curr: EdgeList,
    pub s_flower: EdgeList,
    pub n_flower: EdgeList,
    pub other_n_flower: EdgeList,
    pub fixed_large_ring_stereo: EdgeList,
    pub all_bond: EdgeList,
    pub tot_succes: i32,
    pub b_allowed_n_flower_edges: bool,
    pub ret: i32,
    pub stereo_mask: i32,
    res: BnsTestResult,
}

/// the outcome of `try_flow`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tried {
    /// the BNS path was found and the structure restored with a positive result
    Success,
    /// the path matched but the restoration did not return a positive value;
    /// nothing is rolled back (as in C)
    Unrestored,
    /// the attempt was rolled back
    RolledBack,
}

impl<'a> FixFhCtx<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ip: &'a InputParms,
        bns: &'a mut BnStruct,
        bd: &'a mut BnData,
        ps: &'a mut StrFromInChI,
        at2: &'a mut Vec<InpAtom>,
        at3: &'a mut Vec<InpAtom>,
        p_va: &'a [ValAt],
        tcg: &'a AllTcGroups,
        p_inchi: [Option<&'a InChI>; 2],
        mask: i32,
    ) -> Self {
        FixFhCtx {
            ip,
            bns,
            bd,
            ps,
            at2,
            at3,
            p_va,
            tcg,
            p_inchi,
            num_run_bns: 0,
            mask,
            mask_inv: !(mask as i8),
            c2i: Cmp2FhInchi::default(),
            all_charge: EdgeList::default(),
            curr: EdgeList::default(),
            s_flower: EdgeList::default(),
            n_flower: EdgeList::default(),
            other_n_flower: EdgeList::default(),
            fixed_large_ring_stereo: EdgeList::default(),
            all_bond: EdgeList::default(),
            tot_succes: 0,
            b_allowed_n_flower_edges: false,
            ret: 0,
            stereo_mask: 0,
            res: BnsTestResult::default(),
        }
    }

    /// no fixed-H found in the original or the restored InChI
    fn no_fixed_h(&self) -> bool {
        self.p_inchi[0].is_none_or(|c| c.n_num_h_fixed.is_none())
            && self.ps.p_one_inchi[0]
                .as_ref()
                .is_none_or(|c| c.n_num_h_fixed.is_none())
    }

    /// the (-) charge edge index of the atom
    pub fn minus_edge(&self, iat: usize) -> i32 {
        self.p_va[iat].n_c_minus_group_edge - 1
    }

    /// the (+) charge edge index of the atom
    pub fn plus_edge(&self, iat: usize) -> i32 {
        self.p_va[iat].n_c_plus_group_edge - 1
    }

    /// edge `e` exists and is not forbidden
    pub fn free_edge(&self, e: i32) -> bool {
        e >= 0 && self.bns.edge[e as usize].forbidden == 0
    }

    /// appends to `list`; returns true on an error (`self.ret` is set)
    pub fn add_edge(&mut self, which: L, e: i32) -> bool {
        let list = match which {
            L::AllCharge => &mut self.all_charge,
            L::Curr => &mut self.curr,
            L::SFlower => &mut self.s_flower,
            L::NFlower => &mut self.n_flower,
            L::OtherNFlower => &mut self.other_n_flower,
            L::FixedLargeRingStereo => &mut self.fixed_large_ring_stereo,
            L::AllBond => &mut self.all_bond,
        };
        let r = add_to_edge_list(list, e, INC_ADD_EDGE);
        self.ret = r;
        r != 0
    }

    /// removes the flow 1 from the edge `pe`, forbids it and runs the BNS
    /// test; when the path connects the ends of the edge and `ok(delta
    /// charge)` holds, the structure is restored. Otherwise the change is
    /// rolled back.
    pub fn try_flow(&mut self, pe: usize, ok: impl Fn(i32) -> bool) -> Tried {
        self.try_flow_opt(pe, true, ok)
    }

    /// `try_flow`, optionally without forbidding the edge (case 06c, 06d)
    pub fn try_flow_opt(&mut self, pe: usize, forbid: bool, ok: impl Fn(i32) -> bool) -> Tried {
        let mask = self.mask as i8;
        let v1 = self.bns.edge[pe].neighbor1 as usize;
        let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
        if forbid {
            self.bns.edge[pe].forbidden |= mask;
        }
        self.bns.edge[pe].flow -= 1;
        self.bns.vert[v1].st_edge.flow -= 1;
        self.bns.vert[v2].st_edge.flow -= 1;
        self.bns.tot_st_flow -= 2;
        self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
        let (a, b) = (v1 as i32, v2 as i32);
        if self.ret == 1
            && ((self.res.v_last == a && self.res.v_first == b)
                || (self.res.v_last == b && self.res.v_first == a))
            && ok(self.res.delta_charge)
        {
            self.ret = run_bns_restore_once(self.bns, self.bd);
            if self.ret > 0 {
                self.num_run_bns += 1;
                return Tried::Success;
            }
            Tried::Unrestored
        } else {
            if forbid {
                self.bns.edge[pe].forbidden &= self.mask_inv;
            }
            self.bns.edge[pe].flow += 1;
            self.bns.vert[v1].st_edge.flow += 1;
            self.bns.vert[v2].st_edge.flow += 1;
            self.bns.tot_st_flow += 2;
            Tried::RolledBack
        }
    }

    /// the common end of every case: when something has been fixed,
    /// recalculates the InChI of the structure and compares it with the
    /// original one. Returns true when `FixFixedHRestoredStructure` has to
    /// exit (`self.ret` holds its value).
    pub fn recalc(&mut self, cur_success: i32) -> bool {
        if cur_success == 0 {
            return false;
        }
        self.tot_succes += cur_success;
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
        if self.no_fixed_h() {
            return true;
        }
        self.ret = fill_out_cmp2fh_inchi(self.ps, self.at2, self.p_va, self.p_inchi, &mut self.c2i);
        if self.ret != 0 {
            return true;
        }
        !self.c2i.b_has_difference
    }

    /// the beginning of `FixFixedHRestoredStructure` up to the first case;
    /// returns true when the function has to exit
    pub fn prologue(&mut self, forbidden_stereo_edge_mask: i32) -> bool {
        use super::ichiring::is_bond_in_nmax_memb_ring;
        self.stereo_mask = forbidden_stereo_edge_mask;
        if self.no_fixed_h() {
            return true;
        }
        let num_atoms = self.ps.num_atoms as usize;
        for i in 0..num_atoms {
            let e = self.minus_edge(i);
            if self.free_edge(e) && self.add_edge(L::AllCharge, e) {
                return true;
            }
            let e = self.plus_edge(i);
            if self.free_edge(e) {
                if self.add_edge(L::AllCharge, e) {
                    return true;
                }
                let va = self.p_va[i];
                let mut handled = false;
                if va.c_num_valence_electrons == 5 && va.c_metal == 0 {
                    let j = get_charge_flower_upper_edge(self.bns, e);
                    if j != NO_VERTEX {
                        handled = true;
                        if self.bns.edge[j as usize].forbidden != 0 {
                            continue;
                        }
                        if self.bns.edge[j as usize].flow != 0 {
                            if self.add_edge(L::AllCharge, j) || self.add_edge(L::NFlower, j) {
                                return true;
                            }
                        } else if self.add_edge(L::OtherNFlower, j) {
                            return true;
                        }
                    }
                }
                if !handled && va.c_num_valence_electrons == 6 && va.c_metal == 0 {
                    let j = get_charge_flower_upper_edge(self.bns, e);
                    if j != NO_VERTEX {
                        if self.bns.edge[j as usize].forbidden != 0 {
                            continue;
                        }
                        if self.bns.edge[j as usize].flow != 0 && self.add_edge(L::SFlower, j) {
                            return true;
                        }
                    }
                }
            }
            for j in 0..self.at2[i].valence as usize {
                let k = self.at2[i].neighbor[j] as usize;
                let e = self.bns.iedge[self.bns.vert[i].iedge + j];
                if k < i && self.bns.edge[e as usize].forbidden == 0 && self.add_edge(L::AllBond, e)
                {
                    return true;
                }
            }
        }
        if forbidden_stereo_edge_mask != 0 {
            for i in 0..num_atoms {
                for j in 0..self.at2[i].valence as usize {
                    let k = self.bns.iedge[self.bns.vert[i].iedge + j];
                    if i32::from(self.bns.edge[k as usize].forbidden) == forbidden_stereo_edge_mask
                    {
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
                        if 0 < min_ring && self.add_edge(L::FixedLargeRingStereo, k) {
                            return true;
                        }
                    }
                }
            }
        }
        self.ret = fill_out_extra_fixed_h_data_restr(self.ps);
        if self.ret != 0 {
            return true;
        }
        self.ret = fill_out_cmp2fh_inchi(self.ps, self.at2, self.p_va, self.p_inchi, &mut self.c2i);
        if self.ret != 0 {
            return true;
        }
        let c = &self.c2i;
        !c.b_has_difference
            || (c.c2at.is_empty()
                && c.n_num_tg_revrs == c.n_num_tg_inchi
                && c.n_num_endp_revrs == c.n_num_rem_h_inchi
                && c.n_num_endp_revrs == c.n_num_endp_inchi
                && c.n_num_tg_diff_minus == 0
                && c.n_num_tg_diff_h == 0)
    }

    /// the end of `FixFixedHRestoredStructure`
    pub fn result(&self) -> i32 {
        if self.ret < 0 {
            self.ret
        } else {
            i32::from(self.c2i.b_has_difference && self.tot_succes != 0)
        }
    }
}

/// the generic "move a charge" attempt used by cases 01-03a: the edges
/// `cand` (atoms listed in `atoms`) lose their flow one after the other until
/// `num_try` attempts succeed; `ok` tests the charge change of the BNS path
pub fn try_atoms(
    c: &mut FixFhCtx,
    atoms: &[usize],
    use_plus: bool,
    num_try: usize,
    ok: impl Fn(i32) -> bool + Copy,
    cur_success: &mut i32,
    num_zero_ret: &mut i32,
) {
    let mut i = 0;
    while i < atoms.len() && (*cur_success as usize) < num_try {
        let iat = atoms[i];
        i += 1;
        let pe = if use_plus {
            c.plus_edge(iat)
        } else {
            c.minus_edge(iat)
        } as usize;
        if c.bns.edge[pe].flow == 0 {
            continue;
        }
        match c.try_flow(pe, ok) {
            Tried::Success => *cur_success += 1,
            Tried::RolledBack => *num_zero_ret += i32::from(c.ret == 0),
            Tried::Unrestored => {}
        }
    }
}

impl<'a> FixFhCtx<'a> {
    /// case 01: restored O=AB-O(-), original (-)O-AB=O; move (-) from
    /// B-O(-) to O=A. Returns true to exit.
    pub fn case_01(&mut self) -> bool {
        if self.c2i.c2at.len() < 2 {
            return false;
        }
        let mut dbo: Vec<usize> = Vec::new();
        let mut sbo: Vec<usize> = Vec::new();
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.minus_edge(iat);
            if a.n_val_electr == 6 && self.free_edge(e) {
                let at = &self.at2[iat];
                if sbo.len() < MAX_DIFF_FIXH
                    && a.n_fix_h_inchi == 0
                    && a.n_mob_h_inchi == 0
                    && a.n_fix_h_revrs == -1
                    && a.n_mob_h_revrs == 1
                    && a.n_at_charge_revrs == -1
                    && at.num_h == 0
                    && at.valence == 1
                    && at.chem_bonds_valence == 1
                {
                    sbo.push(iat);
                    if self.add_edge(L::Curr, e) {
                        return true;
                    }
                } else if dbo.len() < MAX_DIFF_FIXH
                    && a.n_fix_h_inchi == -1
                    && a.n_mob_h_inchi == 1
                    && a.n_fix_h_revrs == 0
                    && a.n_mob_h_revrs == 0
                    && a.n_at_charge_revrs == 0
                    && at.num_h == 0
                    && at.valence == 1
                    && at.chem_bonds_valence == 2
                {
                    dbo.push(iat);
                    if self.add_edge(L::Curr, e) {
                        return true;
                    }
                }
            }
        }
        let num_try = sbo.len().min(dbo.len());
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut nz = 0;
            try_atoms(
                self,
                &sbo,
                false,
                num_try,
                |d| d == 1,
                &mut cur_success,
                &mut nz,
            );
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            self.curr.num_edges = 0;
        }
        self.recalc(cur_success)
    }

    /// case 02: restored -O(+)=AB-NH2, original -O-AB=NH2(+); move (+) from
    /// O(+) to NH2
    pub fn case_02(&mut self) -> bool {
        if self.c2i.c2at.is_empty() {
            return false;
        }
        let mut db_o_plus: Vec<usize> = Vec::new();
        let mut sb_nh: Vec<usize> = Vec::new();
        let n_mob_h_inchi: Option<Vec<i8>> = self.p_inchi[1]
            .filter(|c| !c.n_num_h.is_empty())
            .map(|c| c.n_num_h.clone())
            .or_else(|| {
                self.p_inchi[0]
                    .filter(|c| !c.n_num_h.is_empty())
                    .map(|c| c.n_num_h.clone())
            });
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.plus_edge(iat);
            let at = &self.at2[iat];
            if sb_nh.len() < MAX_DIFF_FIXH
                && ((a.n_val_electr == 5 && a.n_period_num == 1) || a.n_val_electr == 6)
                && self.free_edge(e)
                && a.n_fix_h_inchi > 0
                && a.n_fix_h_revrs == 0
                && a.n_at_charge_revrs == 0
                && at.num_h != 0
                && at.valence == at.chem_bonds_valence
            {
                sb_nh.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let e = self.plus_edge(iat);
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            if db_o_plus.len() < MAX_DIFF_FIXH
                && at.charge == 1
                && at.num_h == 0
                && at.valence < at.chem_bonds_valence
                && va.c_metal == 0
                && (va.c_num_valence_electrons == 6
                    || va.c_num_valence_electrons == 7
                    || (va.c_num_valence_electrons == 5 && va.c_periodic_row_number > 1))
                && self.ps.endpoint.as_ref().is_none_or(|ep| ep[i] == 0)
                && !self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0)
                && !n_mob_h_inchi.as_ref().is_some_and(|m| m[i] != 0)
                && self.free_edge(e)
            {
                db_o_plus.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = db_o_plus.len().min(sb_nh.len());
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut num_zero_ret = 0;
            loop {
                // as in C, the counters are kept over the repeats
                try_atoms(
                    self,
                    &sb_nh,
                    true,
                    num_try,
                    |d| d == -1,
                    &mut cur_success,
                    &mut num_zero_ret,
                );
                if num_zero_ret as usize == num_try
                    && !self.b_allowed_n_flower_edges
                    && self.n_flower.num_edges != 0
                {
                    remove_forbidden_edge_mask(self.bns, &self.n_flower, self.mask);
                    self.b_allowed_n_flower_edges = true;
                    continue;
                }
                break;
            }
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// the normalized Mobile-H atoms of the restored structure
    fn at_mobile_h_revrs(&self) -> Option<&[InpAtom]> {
        self.ps.p_one_norm_data[1]
            .as_ref()
            .filter(|d| !d.at.is_empty())
            .map(|d| d.at.as_slice())
    }

    /// case 03: restored -N(-)-AB=O, original -N=AB-O(-); move (-) from
    /// N(-) to =O
    pub fn case_03(&mut self) -> bool {
        let c = &self.c2i;
        if !(!c.c2at.is_empty()
            && c.n_num_tg_revrs == 1
            && (c.n_num_endp_revrs > c.n_num_endp_inchi || c.n_num_tg_inchi > 1))
        {
            return false;
        }
        let mut sb_n_minus: Vec<usize> = Vec::new();
        let mut db_o: Vec<usize> = Vec::new();
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            if db_o.len() < MAX_DIFF_FIXH
                && a.n_val_electr == 6
                && a.endpt_inchi == 0
                && self.free_edge(e)
                && a.n_fix_h_inchi == -1
                && a.n_mob_h_inchi == 1
                && a.endpt_revrs != 0
                && a.n_fix_h_revrs == 0
                && a.n_mob_h_revrs == 0
                && a.n_at_charge_revrs == 0
                && at.num_h == 0
                && at.valence == 1
                && at.chem_bonds_valence == 2
            {
                db_o.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let endp_mh = self
                .at_mobile_h_revrs()
                .is_some_and(|a| a[iat].endpoint != 0);
            if sb_n_minus.len() < MAX_DIFF_FIXH
                && at.charge == -1
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 5
                && va.c_periodic_row_number == 1
                && endp_mh
                && self.free_edge(e)
            {
                sb_n_minus.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = sb_n_minus.len().min(db_o.len());
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut nz = 0;
            try_atoms(
                self,
                &sb_n_minus,
                false,
                num_try,
                |d| d == 1,
                &mut cur_success,
                &mut nz,
            );
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 03a: as case 03 for tautomeric atoms of the original InChI
    pub fn case_03a(&mut self) -> bool {
        let c = &self.c2i;
        if !(c.n_num_tg_revrs == 1
            && (c.n_num_endp_revrs > c.n_num_endp_inchi || c.n_num_tg_inchi > 1))
        {
            return false;
        }
        let mut sb_n_minus: Vec<usize> = Vec::new();
        let mut db_o: Vec<usize> = Vec::new();
        let pn_mob_h_inchi: Option<Vec<i8>> = self.p_inchi[1]
            .filter(|c| !c.n_num_h.is_empty())
            .map(|c| c.n_num_h.clone())
            .or_else(|| {
                self.p_inchi[0]
                    .filter(|c| !c.n_num_h.is_empty())
                    .map(|c| c.n_num_h.clone())
            });
        self.curr.num_edges = 0;
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let endp_mh = self
                .at_mobile_h_revrs()
                .is_some_and(|a| a[iat].endpoint != 0);
            if sb_n_minus.len() < MAX_DIFF_FIXH
                && at.charge == -1
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 5
                && va.c_periodic_row_number == 1
                && endp_mh
                && self.free_edge(e)
            {
                sb_n_minus.push(iat);
            } else if db_o.len() < MAX_DIFF_FIXH
                && at.charge == 0
                && i32::from(at.valence) + 1 == i32::from(at.chem_bonds_valence)
                && va.c_metal == 0
                && va.c_num_valence_electrons == 6
                && endp_mh
                && (self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                    || (pn_mob_h_inchi.as_ref().is_some_and(|m| m[i] == 1)
                        && self.ps.fixed_h.as_ref().is_some_and(|f| f[i] == -1)))
                && self.free_edge(e)
            {
                db_o.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = sb_n_minus.len().min(db_o.len());
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut nz = 0;
            try_atoms(
                self,
                &sb_n_minus,
                false,
                num_try,
                |d| d == 1,
                &mut cur_success,
                &mut nz,
            );
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// `nNum_H` of the Mobile-H (or else the Fixed-H) original InChI
    fn mob_h_inchi(&self) -> Option<Vec<i8>> {
        self.p_inchi[1]
            .filter(|c| !c.n_num_h.is_empty())
            .map(|c| c.n_num_h.clone())
            .or_else(|| {
                self.p_inchi[0]
                    .filter(|c| !c.n_num_h.is_empty())
                    .map(|c| c.n_num_h.clone())
            })
    }

    /// case 04: restored OH(+)=AB-O-, original HO-AB=O(+)-; move (+) from
    /// O(+) to the neutral atom
    pub fn case_04(&mut self) -> bool {
        let c = &self.c2i;
        if !(!c.c2at.is_empty()
            && c.n_num_tg_inchi == 1
            && (c.n_num_endp_revrs < c.n_num_endp_inchi || c.n_num_tg_revrs > 1))
        {
            return false;
        }
        let mob = self.mob_h_inchi();
        let mut num_db_charged = 0;
        let mut sb_neutr: Vec<usize> = Vec::new();
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let e = self.plus_edge(iat);
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let fixed = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0);
            let endpoint = self.ps.endpoint.as_ref().is_some_and(|f| f[i] != 0);
            if num_db_charged < MAX_DIFF_FIXH
                && at.charge == 1
                && at.num_h != 0
                && at.valence < at.chem_bonds_valence
                && va.c_metal == 0
                && (va.c_num_valence_electrons == 6
                    || (va.c_num_valence_electrons == 5 && va.c_periodic_row_number == 1))
                && endpoint
                && fixed
                && self.free_edge(e)
            {
                num_db_charged += 1;
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            } else if sb_neutr.len() < MAX_DIFF_FIXH
                && at.charge == 0
                && at.num_h == 0
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && (va.c_num_valence_electrons == 6 || va.c_num_valence_electrons == 5)
                && !fixed
                && !mob.as_ref().is_some_and(|m| m[i] != 0)
                && self.free_edge(e)
            {
                sb_neutr.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = sb_neutr.len().min(num_db_charged);
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut nz = 0;
            try_atoms(
                self,
                &sb_neutr,
                true,
                num_try,
                |d| d == -1,
                &mut cur_success,
                &mut nz,
            );
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 05: restored O=AB-NH, original (-)O-AB=NH(+); separate charges
    pub fn case_05(&mut self) -> bool {
        if self.c2i.c2at.len() <= 1 {
            return false;
        }
        let mut db_o: Vec<usize> = Vec::new();
        let mut sb_nh: Vec<usize> = Vec::new();
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let ep = self.plus_edge(iat);
            let em = self.minus_edge(iat);
            let at = &self.at2[iat];
            if sb_nh.len() < MAX_DIFF_FIXH
                && ((a.n_val_electr == 5 && a.n_period_num == 1) || a.n_val_electr == 6)
                && a.endpt_inchi == 0
                && self.free_edge(ep)
                && a.n_fix_h_inchi == 1
                && a.n_fix_h_revrs == 0
                && a.n_mob_h_revrs != 0
                && a.n_at_charge_revrs == 0
                && at.num_h != 0
                && a.endpt_revrs == 0
                && at.valence == at.chem_bonds_valence
            {
                sb_nh.push(iat);
                if self.add_edge(L::Curr, ep) {
                    return true;
                }
            } else if db_o.len() < MAX_DIFF_FIXH
                && a.n_val_electr == 6
                && a.endpt_inchi == 0
                && self.free_edge(em)
                && a.n_fix_h_inchi == -1
                && a.n_mob_h_inchi == 1
                && a.n_fix_h_revrs == 0
                && a.n_mob_h_revrs == 0
                && a.n_at_charge_revrs == 0
                && at.num_h == 0
                && a.endpt_revrs == 0
                && i32::from(at.valence) + 1 == i32::from(at.chem_bonds_valence)
            {
                db_o.push(iat);
                if self.add_edge(L::Curr, em) {
                    return true;
                }
            }
        }
        let num_try = db_o.len().min(sb_nh.len());
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut nz = 0;
            try_atoms(
                self,
                &sb_nh,
                true,
                num_try,
                |d| d == 1,
                &mut cur_success,
                &mut nz,
            );
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// case 06c: restored -NH- or -NH(+), original -NH-; move the (+)
    /// charges to the -NH- unless it already has one blocked by an adjacent
    /// (+)
    pub fn case_06c(&mut self) -> bool {
        use super::ichirvr1::find_in_edge_list;
        if !(self.ps.fixed_h.is_some()
            && self.ps.endpoint.is_some()
            && self.c2i.n_charge_fix_h_inchi > 0
            && self.c2i.n_charge_fix_h_inchi > self.c2i.n_charge_mob_h_inchi)
        {
            return false;
        }
        let mut cur_charge = EdgeList::default();
        let mut cur_success = 0;
        self.curr.num_edges = 0;
        'work: {
            for ci in 0..self.c2i.c2at.len() {
                let a = self.c2i.c2at[ci];
                let iat = a.atom_number as usize;
                let e = self.plus_edge(iat);
                if (a.n_val_electr == 6 || (a.n_val_electr == 5 && a.n_period_num == 1))
                    && a.endpt_inchi == 0
                    && self.free_edge(e)
                {
                    let at = &self.at2[iat];
                    if a.n_fix_h_inchi == 1
                        && a.n_mob_h_inchi == 0
                        && a.n_fix_h_revrs == 0
                        && a.n_mob_h_revrs == 1
                        && at.num_h != 0
                        && at.valence == at.chem_bonds_valence
                        && self.add_edge(L::Curr, e)
                    {
                        return true;
                    }
                }
            }
            for i in 0..self.ps.num_atoms as usize {
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let at = &self.at2[iat];
                if self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                    || at.charge != 1
                    || at.radical != 0
                    || self.p_va[iat].c_metal != 0
                {
                    continue;
                }
                let e = self.plus_edge(iat);
                if 0 <= e
                    && self.bns.edge[e as usize].forbidden == 0
                    && self.bns.edge[e as usize].flow == 0
                    && self.p_va[iat].c_num_valence_electrons >= 5
                {
                    for j in 0..self.at2[iat].valence as usize {
                        let k = self.at2[iat].neighbor[j] as usize;
                        let e2 = self.plus_edge(k);
                        if self.at2[k].charge == 1
                            && self.p_va[k].c_metal == 0
                            && 0 <= e2
                            && self.bns.edge[e2 as usize].forbidden == 0
                            && self.bns.edge[e2 as usize].flow == 0
                        {
                            for ee in [e, e2] {
                                if find_in_edge_list(&self.curr, ee) < 0
                                    && find_in_edge_list(&cur_charge, ee) < 0
                                {
                                    self.ret = add_to_edge_list(&mut cur_charge, ee, INC_ADD_EDGE);
                                    if self.ret != 0 {
                                        break 'work;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            let num_try = self.curr.num_edges.min(cur_charge.num_edges) as usize;
            if num_try != 0 {
                set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                remove_forbidden_edge_mask(self.bns, &cur_charge, self.mask);
                let mut i = 0;
                while i < self.curr.num_edges as usize && (cur_success as usize) < num_try {
                    let pe = self.curr.pn_edges[i] as usize;
                    i += 1;
                    if self.bns.edge[pe].flow == 0 {
                        continue;
                    }
                    if self.try_flow_opt(pe, false, |d| d == -1) == Tried::Success {
                        cur_success += 1;
                    }
                }
                remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            }
        }
        self.curr.num_edges = 0;
        if self.ret < 0 {
            return true;
        }
        self.recalc(cur_success)
    }

    /// case 06d: restored XH(+)=-AB-NH, original XH-=AB=NH(+); move (+) from
    /// X to NH
    pub fn case_06d(&mut self) -> bool {
        if self.c2i.c2at.len() < 2 {
            return false;
        }
        let mut cur_charge = EdgeList::default();
        let mut cur_success = 0;
        self.curr.num_edges = 0;
        'work: {
            for ci in 0..self.c2i.c2at.len() {
                let a = self.c2i.c2at[ci];
                let iat = a.atom_number as usize;
                let e = self.plus_edge(iat);
                let cond_a = a.n_mob_h_revrs + 1 == a.n_fix_h_revrs
                    && a.n_fix_h_revrs > 0
                    && a.endpt_revrs == 0
                    && a.n_at_charge_revrs == 1
                    && ((a.n_fix_h_inchi == 0 && a.n_mob_h_inchi == a.n_fix_h_revrs)
                        || (a.n_fix_h_inchi == a.n_fix_h_revrs && a.endpt_inchi != 0));
                if cond_a
                    && 0 <= e
                    && self.bns.edge[e as usize].forbidden == 0
                    && self.bns.edge[e as usize].flow == 0
                {
                    self.ret = add_to_edge_list(&mut cur_charge, e, INC_ADD_EDGE);
                    if self.ret != 0 {
                        break 'work;
                    }
                } else {
                    let cond_b = a.n_mob_h_inchi + 1 == a.n_fix_h_inchi
                        && a.n_fix_h_inchi > 0
                        && a.endpt_inchi == 0
                        && a.n_at_charge_revrs == 0
                        && ((a.n_fix_h_revrs == 0 && a.n_mob_h_revrs == a.n_fix_h_inchi)
                            || (a.n_fix_h_revrs == a.n_fix_h_inchi && a.endpt_revrs != 0));
                    if cond_b
                        && 0 <= e
                        && self.bns.edge[e as usize].forbidden == 0
                        && self.bns.edge[e as usize].flow != 0
                    {
                        if self.add_edge(L::Curr, e) {
                            break 'work;
                        }
                    }
                }
            }
            let num_try = self.curr.num_edges.min(cur_charge.num_edges) as usize;
            if num_try != 0 {
                let s_may_be_forbidden = self.s_flower.num_edges > 0;
                let mut pass = i32::from(s_may_be_forbidden);
                while pass >= 0 {
                    if pass != 0 {
                        set_forbidden_edge_mask(self.bns, &self.s_flower, self.mask);
                    }
                    set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                    remove_forbidden_edge_mask(self.bns, &cur_charge, self.mask);
                    let mut i = 0;
                    while i < self.curr.num_edges as usize && (cur_success as usize) < num_try {
                        let pe = self.curr.pn_edges[i] as usize;
                        i += 1;
                        if self.bns.edge[pe].flow == 0 {
                            continue;
                        }
                        if self.try_flow_opt(pe, false, |d| d == -1) == Tried::Success {
                            cur_success += 1;
                        }
                    }
                    remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                    remove_forbidden_edge_mask(self.bns, &self.s_flower, self.mask);
                    pass -= 1;
                }
            }
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// the preprocessed Mobile-H atoms (`at_fixed_bonds`) of the restored
    /// structure
    fn atf_mobile_h_revrs(&self) -> Option<&[InpAtom]> {
        self.ps.p_one_norm_data[1]
            .as_ref()
            .and_then(|d| d.at_fixed_bonds.as_deref())
    }

    /// case 06: restored NHn(+)=AB-NHm, original NHn-AB=NHm(+); move (+)
    /// from NHn(+) to NHn
    pub fn case_06(&mut self) -> bool {
        if self.c2i.c2at.len() < 2 {
            return false;
        }
        let mut sb_nhm_neutr: Vec<usize> = Vec::new();
        let mut num_db_nhn_plus = 0;
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.plus_edge(iat);
            if (a.n_val_electr == 6 || (a.n_val_electr == 5 && a.n_period_num == 1))
                && a.endpt_inchi == 0
                && self.free_edge(e)
            {
                let at = &self.at2[iat];
                if sb_nhm_neutr.len() < MAX_DIFF_FIXH
                    && a.n_fix_h_inchi == 1
                    && a.n_fix_h_revrs == 0
                    && a.n_at_charge_revrs == 0
                    && at.num_h != 0
                    && at.valence == at.chem_bonds_valence
                {
                    sb_nhm_neutr.push(iat);
                    if self.add_edge(L::Curr, e) {
                        return true;
                    }
                } else if num_db_nhn_plus < MAX_DIFF_FIXH
                    && a.n_fix_h_inchi == 0
                    && a.n_fix_h_revrs == 1
                    && a.n_at_charge_revrs == 1
                    && at.num_h != 0
                    && at.valence < at.chem_bonds_valence
                {
                    num_db_nhn_plus += 1;
                    if self.add_edge(L::Curr, e) {
                        return true;
                    }
                }
            }
        }
        let num_try = sb_nhm_neutr.len().min(num_db_nhn_plus);
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut nz = 0;
            try_atoms(
                self,
                &sb_nhm_neutr,
                true,
                num_try,
                |d| d == -1,
                &mut cur_success,
                &mut nz,
            );
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// shared body of cases 06a and 06b: moves (+) from the atoms of
    /// `cur_charge` to the -NH-/-OH atoms of `self.curr`
    fn plus_move(&mut self, cur_charge: &EdgeList, two_pass: bool) -> i32 {
        let num_try = self.curr.num_edges.min(cur_charge.num_edges) as usize;
        let mut cur_success = 0;
        if num_try != 0 {
            let mut pass = i32::from(two_pass && self.s_flower.num_edges > 0);
            while pass >= 0 {
                if pass != 0 {
                    set_forbidden_edge_mask(self.bns, &self.s_flower, self.mask);
                }
                set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                remove_forbidden_edge_mask(self.bns, cur_charge, self.mask);
                let mut i = 0;
                while i < self.curr.num_edges as usize && (cur_success as usize) < num_try {
                    let pe = self.curr.pn_edges[i] as usize;
                    i += 1;
                    if self.bns.edge[pe].flow == 0 {
                        continue;
                    }
                    if self.try_flow_opt(pe, false, |d| d == -1) == Tried::Success {
                        cur_success += 1;
                    }
                }
                remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                if two_pass {
                    remove_forbidden_edge_mask(self.bns, &self.s_flower, self.mask);
                }
                pass -= 1;
            }
        }
        cur_success
    }

    /// the loop of cases 06a/06b collecting the -NH- (`curr`) and the
    /// charged (`cur_charge`) candidates; true = error exit of the case
    fn collect_06ab(&mut self, cur_charge: &mut EdgeList, is_06b: bool) -> bool {
        let mob = self.mob_h_inchi();
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            if self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0) {
                continue;
            }
            let e = self.plus_edge(iat);
            let at = &self.at2[iat];
            let first = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0)
                && mob.as_ref().is_some_and(|m| m[i] == 0)
                && at.charge == 0
                && at.radical == 0
                && 0 <= e
                && self.bns.edge[e as usize].forbidden == 0
                && self.bns.edge[e as usize].flow != 0;
            if first {
                self.ret = add_to_edge_list(&mut self.curr, e, INC_ADD_EDGE);
                if self.ret != 0 {
                    return true;
                }
            }
            // as in C the second test is made also when the first one has
            // passed (it cannot be true then)
            let va = &self.p_va[iat];
            let atf_charge = self.atf_mobile_h_revrs().map(|a| a[iat].charge);
            let kind_ok = if is_06b {
                va.c_num_valence_electrons == 6 || va.c_periodic_row_number == 7
            } else {
                va.c_num_valence_electrons == 5 && va.c_periodic_row_number == 1
            };
            let want_charge = if is_06b { 1 } else { 0 };
            let at = &self.at2[iat];
            if at.charge == 1
                && at.num_h == 0
                && kind_ok
                && atf_charge == Some(want_charge)
                && 0 <= e
                && self.bns.edge[e as usize].forbidden == 0
                && self.bns.edge[e as usize].flow == 0
            {
                self.ret = add_to_edge_list(cur_charge, e, INC_ADD_EDGE);
                if self.ret != 0 {
                    return true;
                }
            }
        }
        false
    }

    /// case 06a: restored N'(+)=-AB-NH, original N'-=AB=NH(+); move (+) from
    /// N' to NH
    pub fn case_06a(&mut self) -> bool {
        let c = &self.c2i;
        if !(((c.n_num_tg_inchi > c.n_num_tg_revrs && c.n_num_tg_revrs == 1)
            || c.n_num_endp_inchi < c.n_num_endp_revrs)
            && self.ps.n_num_removed_protons_mobh_inchi == self.ps.one_ti.tni.n_num_removed_protons
            && self.ps.fixed_h.is_some()
            && self.ps.endpoint.is_some()
            && self.atf_mobile_h_revrs().is_some())
        {
            return false;
        }
        let mut cur_charge = EdgeList::default();
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        if !self.collect_06ab(&mut cur_charge, false) {
            cur_success = self.plus_move(&cur_charge, false);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 06b: restored X(+)=-AB-NH, original X-=AB=NH(+); move (+) from
    /// X to NH
    pub fn case_06b(&mut self) -> bool {
        let c = &self.c2i;
        let removed = (
            self.ps.n_num_removed_protons_mobh_inchi,
            self.ps.one_ti.tni.n_num_removed_protons,
        );
        if !(((c.n_num_tg_inchi > c.n_num_tg_revrs && c.n_num_tg_revrs == 1)
            || c.n_num_endp_inchi < c.n_num_endp_revrs)
            && removed.0 >= removed.1
            && self.ps.fixed_h.is_some()
            && self.ps.endpoint.is_some()
            && self.atf_mobile_h_revrs().is_some())
        {
            return false;
        }
        let mut cur_charge = EdgeList::default();
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        if !self.collect_06ab(&mut cur_charge, true) {
            cur_success = self.plus_move(&cur_charge, true);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// the common tail of cases 07, 07a and 08: tries to move (-) from the
    /// single-bonded atoms `sb` (the C `num_try` is `num_try`)
    fn minus_move(&mut self, sb: &[usize], num_try: usize, stereo: bool) -> i32 {
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            if stereo && self.stereo_mask != 0 {
                remove_forbidden_edge_mask(
                    self.bns,
                    &self.fixed_large_ring_stereo,
                    self.stereo_mask,
                );
            }
            let mut nz = 0;
            try_atoms(
                self,
                sb,
                false,
                num_try,
                |d| d == 1,
                &mut cur_success,
                &mut nz,
            );
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            if stereo && self.stereo_mask != 0 {
                set_forbidden_edge_mask(self.bns, &self.fixed_large_ring_stereo, self.stereo_mask);
            }
        }
        self.curr.num_edges = 0;
        cur_success
    }

    /// the DB=O collection loop of cases 07 and 07a; true = error exit
    fn collect_07_db(&mut self, db_o: &mut Vec<usize>) -> bool {
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            if db_o.len() < MAX_DIFF_FIXH
                && a.n_val_electr == 6
                && a.endpt_inchi == 0
                && self.free_edge(e)
                && a.n_fix_h_inchi == -1
                && a.n_mob_h_inchi == 1
                && a.n_fix_h_revrs == 0
                && a.n_mob_h_revrs == 0
                && a.n_at_charge_revrs == 0
                && at.num_h == 0
                && at.valence < at.chem_bonds_valence
            {
                db_o.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        false
    }

    /// case 07: restored O(-)-AB=O, original O=AB-O(-); move (-) from
    /// O(-)-AB to AB=O
    pub fn case_07(&mut self) -> bool {
        if self.c2i.c2at.is_empty() {
            return false;
        }
        let mut db_o: Vec<usize> = Vec::new();
        let mut sb_o: Vec<usize> = Vec::new();
        let mob = self.mob_h_inchi();
        if self.collect_07_db(&mut db_o) {
            return true;
        }
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let endp_mh = self
                .at_mobile_h_revrs()
                .is_some_and(|a| a[iat].endpoint != 0);
            if sb_o.len() < MAX_DIFF_FIXH
                && at.charge == -1
                && at.num_h == 0
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 6
                && endp_mh
                && !self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0)
                && !mob.as_ref().is_some_and(|m| m[i] != 0)
                && self.free_edge(e)
            {
                sb_o.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = sb_o.len().min(db_o.len());
        let cur_success = self.minus_move(&sb_o, num_try, false);
        self.recalc(cur_success)
    }

    /// case 07a: restored O(-)-N(V)B=O, original O=N(V)B-O(-)
    pub fn case_07a(&mut self) -> bool {
        if self.c2i.c2at.is_empty() {
            return false;
        }
        let mut db_o: Vec<usize> = Vec::new();
        let mut sb_o: Vec<usize> = Vec::new();
        let mob = self.mob_h_inchi();
        if self.collect_07_db(&mut db_o) {
            return true;
        }
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let mut ok = sb_o.len() < MAX_DIFF_FIXH
                && at.charge == -1
                && at.num_h == 0
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 6
                && !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                && !self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0)
                && !mob.as_ref().is_some_and(|m| m[i] != 0)
                && at.valence == 1;
            if ok {
                let n = at.neighbor[0] as usize;
                let an = &self.at2[n];
                ok = an.chem_bonds_valence == 5
                    && an.charge == 0
                    && self.p_va[n].c_num_valence_electrons == 5
                    && self.free_edge(e);
            }
            if ok {
                sb_o.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = sb_o.len().min(db_o.len());
        let cur_success = self.minus_move(&sb_o, num_try, false);
        self.recalc(cur_success)
    }

    /// case 08: restored O(-)-AB=N-, original O=AB-N(-)- (ADP in InChI)
    pub fn case_08(&mut self) -> bool {
        let c = &self.c2i;
        if !(c.n_num_tg_inchi == 1
            && (c.n_num_endp_revrs < c.n_num_endp_inchi || c.n_num_tg_revrs > 1))
        {
            return false;
        }
        let mut db_n: Vec<usize> = Vec::new();
        let mut sb_o: Vec<usize> = Vec::new();
        let mob = self.mob_h_inchi();
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let endp = self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0);
            let fixed = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0);
            let mobh = mob.as_ref().is_some_and(|m| m[i] != 0);
            if sb_o.len() < MAX_DIFF_FIXH
                && at.charge == -1
                && at.num_h == 0
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 6
                && endp
                && !mobh
                && self.free_edge(e)
            {
                sb_o.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            } else if db_n.len() < MAX_DIFF_FIXH
                && at.charge == 0
                && at.num_h == 0
                && at.sb_parity[0] == 0
                && at.valence < at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 5
                && va.c_periodic_row_number == 1
                && endp
                && !fixed
                && !mobh
                && self.free_edge(e)
            {
                db_n.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = db_n.len().min(sb_o.len());
        let cur_success = self.minus_move(&sb_o, num_try, true);
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// the vertex-edge index of the `j`-th edge of vertex `v`
    fn vert_edge(&self, v: usize, j: usize) -> usize {
        self.bns.iedge[self.bns.vert[v].iedge + j] as usize
    }

    /// case 09: restored NH2(+)=C-NH2, original NH2-C(+)-NH2; move (+) from
    /// NH2(+) to C
    pub fn case_09(&mut self) -> bool {
        if self.c2i.c2at.len() < 2 {
            return false;
        }
        let mut cur_success = 0;
        for i in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[i];
            let iat = a.atom_number as usize;
            let e = self.plus_edge(iat);
            let at = &self.at2[iat];
            if !((a.n_val_electr == 5 && a.n_period_num == 1)
                && a.endpt_inchi == 0
                && a.n_fix_h_inchi == 0
                && a.n_mob_h_inchi != 0
                && a.endpt_revrs != 0
                && a.n_fix_h_revrs != 0
                && a.n_mob_h_revrs == 0
                && a.n_at_charge_revrs == 1
                && i32::from(at.valence) + 1 == i32::from(at.chem_bonds_valence)
                && self.free_edge(e))
            {
                continue;
            }
            let val = at.valence as usize;
            let Some(j) = (0..val).find(|&j| i32::from(at.bond_type[j]) == BOND_TYPE_DOUBLE) else {
                continue;
            };
            let e_nc = self.vert_edge(iat, j);
            let iat_c = at.neighbor[j] as usize;
            let ac = &self.at2[iat_c];
            let vc = &self.p_va[iat_c];
            let e_cplus_c = self.plus_edge(iat_c);
            if vc.c_num_valence_electrons != 4
                || vc.c_metal != 0
                || ac.charge != 0
                || ac.valence != 3
                || i32::from(ac.valence) + 1 != i32::from(ac.chem_bonds_valence)
                || e_cplus_c < 0
                || self.bns.edge[e_cplus_c as usize].forbidden != 0
            {
                continue;
            }
            for j in 0..ac.valence as usize {
                let iat_nh2 = self.at2[iat_c].neighbor[j] as usize;
                let vn = &self.p_va[iat_nh2];
                if iat_nh2 == iat
                    || vn.c_num_valence_electrons != 5
                    || vn.c_periodic_row_number != 1
                    || self.at2[iat_nh2].num_h == 0
                    || self.at2[iat_nh2].charge != 0
                {
                    continue;
                }
                let Some(i_nh2) = self
                    .c2i
                    .c2at
                    .iter()
                    .position(|c| c.atom_number as usize == iat_nh2)
                else {
                    continue;
                };
                let b = self.c2i.c2at[i_nh2];
                let an = &self.at2[iat_nh2];
                if (b.n_val_electr == 5 && b.n_period_num == 1)
                    && b.endpt_inchi == 0
                    && b.n_fix_h_inchi == 0
                    && b.n_mob_h_inchi != 0
                    && b.endpt_revrs != 0
                    && b.n_fix_h_revrs != 0
                    && b.n_mob_h_revrs == 0
                    && b.n_at_charge_revrs == 0
                    && an.valence == an.chem_bonds_valence
                {
                    if self.add_edge(L::Curr, e) || self.add_edge(L::Curr, e_cplus_c) {
                        return true;
                    }
                    set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                    remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
                    if self.bns.edge[e_nc].flow == 0 {
                        continue;
                    }
                    if self.try_flow(e_nc, |d| d == 0) == Tried::Success {
                        cur_success += 1;
                    }
                    self.bns.edge[e_nc].forbidden &= self.mask_inv;
                    remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                    self.curr.num_edges = 0;
                    break;
                }
            }
        }
        self.recalc(cur_success)
    }

    /// case 10: restored NH2-X(+)-NH-, original NH2(+)=X-NH-; move (+) from
    /// X(+) to NH2 or NH
    pub fn case_10(&mut self) -> bool {
        if self.c2i.c2at.len() < 2 {
            return false;
        }
        let mut cur_success = 0;
        for i in 0..self.c2i.c2at.len() {
            if self.c2i.c2at[i].n_value != 0 {
                continue;
            }
            let a = self.c2i.c2at[i];
            let iat = a.atom_number as usize;
            let e = self.plus_edge(iat);
            let at = &self.at2[iat];
            if !((a.n_val_electr == 6 || (a.n_val_electr == 5 && a.n_period_num == 1))
                && a.endpt_inchi != 0
                && a.n_fix_h_inchi != 0
                && a.n_mob_h_inchi == 0
                && a.endpt_revrs == 0
                && a.n_fix_h_revrs == 0
                && a.n_mob_h_revrs != 0
                && a.n_at_charge_revrs == 0
                && at.valence == at.chem_bonds_valence
                && self.free_edge(e))
            {
                continue;
            }
            let mut b_continue = true;
            let mut j1 = 0;
            while j1 < self.at2[iat].valence as usize && b_continue {
                let jj1 = j1;
                j1 += 1;
                let at = &self.at2[iat];
                let iat_c = at.neighbor[jj1] as usize;
                let ac = &self.at2[iat_c];
                let nve = self.p_va[iat_c].c_num_valence_electrons;
                let e_cplus_c = self.plus_edge(iat_c);
                if !(i32::from(at.bond_type[jj1]) == BOND_TYPE_SINGLE
                    && ac.charge == 1
                    && (4..=6).contains(&nve)
                    && ac.valence == ac.chem_bonds_valence
                    && self.free_edge(e_cplus_c))
                {
                    continue;
                }
                let mut j2 = 0;
                while j2 < self.at2[iat_c].valence as usize && b_continue {
                    let jj2 = j2;
                    j2 += 1;
                    let ac = &self.at2[iat_c];
                    let iat_nh2 = ac.neighbor[jj2] as usize;
                    let an = &self.at2[iat_nh2];
                    let nve2 = self.p_va[iat_nh2].c_num_valence_electrons;
                    let e_cplus_nh2 = self.plus_edge(iat_nh2);
                    if !(i32::from(ac.bond_type[jj2]) == BOND_TYPE_SINGLE
                        && iat != iat_nh2
                        && an.charge == 0
                        && an.num_h != 0
                        && (nve2 == 5 || nve2 == 6)
                        && an.valence == an.chem_bonds_valence
                        && self.free_edge(e_cplus_nh2))
                    {
                        continue;
                    }
                    for i_nh2 in 0..self.c2i.c2at.len() {
                        let b = self.c2i.c2at[i_nh2];
                        if iat_nh2 != b.atom_number as usize || b.n_value != 0 {
                            continue;
                        }
                        if b.endpt_inchi != 0
                            && b.n_fix_h_inchi != 0
                            && b.n_mob_h_inchi == 0
                            && b.endpt_revrs == 0
                            && b.n_fix_h_revrs == 0
                            && b.n_mob_h_revrs != 0
                            && b.n_at_charge_revrs == 0
                        {
                            let pe_iat = self.vert_edge(iat, jj1);
                            let pe_i_nh2 = self.vert_edge(iat_c, jj2);
                            let pe = if self.bns.edge[pe_iat].forbidden == 0
                                && self.bns.edge[e as usize].flow != 0
                            {
                                e as usize
                            } else if self.bns.edge[pe_i_nh2].forbidden == 0
                                && self.bns.edge[e_cplus_nh2 as usize].flow != 0
                            {
                                e_cplus_nh2 as usize
                            } else {
                                continue;
                            };
                            if self.add_edge(L::Curr, e_cplus_c) {
                                return true;
                            }
                            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
                            if self.try_flow_opt(pe, false, |d| d == -1) == Tried::Success {
                                cur_success += 1;
                                b_continue = false;
                                self.c2i.c2at[i].n_value = 1;
                                self.c2i.c2at[i_nh2].n_value = 1;
                            }
                            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                            self.curr.num_edges = 0;
                            break;
                        }
                    }
                }
            }
        }
        self.recalc(cur_success)
    }

    /// case 11: restored NH(+)=AB-N<, original NH-AB=N(+)<; move (+) from
    /// NH(+) to -N<
    pub fn case_11(&mut self) -> bool {
        let c = &self.c2i;
        if !(c.n_num_tg_inchi == 1
            && (c.n_num_endp_revrs < c.n_num_endp_inchi || c.n_num_tg_revrs > 1))
        {
            return false;
        }
        let mut sb_neutr: Vec<usize> = Vec::new();
        let mut db_charged: Vec<usize> = Vec::new();
        let mob = self.mob_h_inchi();
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let fixed = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0);
            let mobh = mob.as_ref().is_some_and(|m| m[i] != 0);
            let e = self.plus_edge(iat);
            if db_charged.len() < MAX_DIFF_FIXH
                && at.charge == 1
                && at.num_h != 0
                && at.valence < at.chem_bonds_valence
                && va.c_metal == 0
                && ((va.c_num_valence_electrons == 5 && va.c_periodic_row_number == 1)
                    || va.c_num_valence_electrons == 6)
                && fixed
                && self.free_edge(e)
            {
                db_charged.push(iat);
            } else if sb_neutr.len() < MAX_DIFF_FIXH
                && at.charge == 0
                && at.num_h == 0
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 5
                && va.c_periodic_row_number == 1
                && !fixed
                && !mobh
                && self.free_edge(e)
            {
                sb_neutr.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = sb_neutr.len().min(db_charged.len());
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut i = 0;
            while i < db_charged.len() && (cur_success as usize) < num_try {
                let iat = db_charged[i];
                i += 1;
                let pe = self.plus_edge(iat) as usize;
                if self.bns.edge[pe].flow != 0 {
                    continue;
                }
                let v1 = self.bns.edge[pe].neighbor1 as usize;
                let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
                let mut found1 = None;
                for j in (0..self.bns.vert[v1].num_adj_edges as usize).rev() {
                    let pe1n = self.vert_edge(v1, j);
                    if self.bns.edge[pe1n].flow != 0 && self.bns.edge[pe1n].forbidden == 0 {
                        found1 = Some((pe1n, (self.bns.edge[pe1n].neighbor12 as usize) ^ v1));
                        break;
                    }
                }
                let Some((pe1n, v1n)) = found1 else { continue };
                let mut found2 = None;
                for j in (0..(self.bns.vert[v2].num_adj_edges as i32 - 1)).rev() {
                    let pe2n = self.vert_edge(v2, j as usize);
                    if self.bns.edge[pe2n].flow != 0 && self.bns.edge[pe2n].forbidden == 0 {
                        found2 = Some((pe2n, (self.bns.edge[pe2n].neighbor12 as usize) ^ v2));
                        break;
                    }
                }
                let Some((pe2n, v2n)) = found2 else { continue };

                self.bns.edge[pe].flow += 1;
                self.bns.edge[pe1n].flow -= 1;
                self.bns.edge[pe2n].flow -= 1;
                self.bns.vert[v1n].st_edge.flow -= 1;
                self.bns.vert[v2n].st_edge.flow -= 1;
                self.bns.tot_st_flow -= 2;
                self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                let (a, b) = (v1n as i32, v2n as i32);
                if self.ret == 1
                    && ((self.res.v_last == a && self.res.v_first == b)
                        || (self.res.v_last == b && self.res.v_first == a))
                    && (self.res.delta_charge == 0 || self.res.delta_charge == 1)
                {
                    self.ret = run_bns_restore_once(self.bns, self.bd);
                    if self.ret > 0 {
                        self.num_run_bns += 1;
                        cur_success += 1;
                    }
                } else {
                    self.bns.edge[pe].flow -= 1;
                    self.bns.edge[pe1n].flow += 1;
                    self.bns.edge[pe2n].flow += 1;
                    self.bns.vert[v1n].st_edge.flow += 1;
                    self.bns.vert[v2n].st_edge.flow += 1;
                    self.bns.tot_st_flow += 2;
                }
            }
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// the `O=X` double bond attempt shared by cases 12 and 13: the edge to
    /// the first neighbor of `iat` loses its flow, the BNS has to add two
    /// charges
    fn try_o_db(&mut self, db_o: &[usize], num_try: usize) -> i32 {
        let mut cur_success = 0;
        set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
        let mut i = 0;
        while i < db_o.len() && (cur_success as usize) < num_try {
            let iat = db_o[i];
            i += 1;
            let pe_cminus = self.minus_edge(iat) as usize;
            self.bns.edge[pe_cminus].forbidden &= self.mask_inv;
            let pe = self.vert_edge(iat, 0);
            if self.bns.edge[pe].flow == 0 {
                continue;
            }
            if self.try_flow(pe, |d| d == 2) == Tried::Success {
                cur_success += 1;
            }
            self.bns.edge[pe].forbidden &= self.mask_inv;
        }
        remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        cur_success
    }

    /// case 12: restored O=AB-N<, original (-)O-AB=N(+)<
    pub fn case_12(&mut self) -> bool {
        let c = &self.c2i;
        if !(!c.c2at.is_empty()
            && c.n_num_tg_inchi == 1
            && c.n_num_rem_h_inchi >= -1
            && (c.n_num_endp_inchi > c.n_num_endp_revrs || c.n_num_tg_revrs > 1))
        {
            return false;
        }
        let mut only_n_v = true;
        loop {
            let mut sb_n_neutr = 0usize;
            let mut db_o: Vec<usize> = Vec::new();
            let mut num_n_v = 0;
            let mob = self.mob_h_inchi();
            for ci in 0..self.c2i.c2at.len() {
                let a = self.c2i.c2at[ci];
                let iat = a.atom_number as usize;
                let e = self.minus_edge(iat);
                let at = &self.at2[iat];
                if db_o.len() < MAX_DIFF_FIXH
                    && (a.n_val_electr == 6 || (a.n_val_electr == 5 && a.n_period_num == 1))
                    && a.endpt_inchi != 0
                    && self.free_edge(e)
                    && a.n_fix_h_inchi == 0
                    && a.n_mob_h_inchi == 0
                    && a.endpt_revrs == 0
                    && a.n_fix_h_revrs == 0
                    && a.n_mob_h_revrs == 0
                    && a.n_at_charge_revrs == 0
                    && at.num_h == 0
                    && (if a.n_val_electr == 6 {
                        at.valence == 1 && at.chem_bonds_valence == 2
                    } else {
                        a.n_val_electr == 5 && at.valence == 2 && at.chem_bonds_valence == 3
                    })
                {
                    db_o.push(iat);
                }
            }
            for i in 0..self.ps.num_atoms as usize {
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let at = &self.at2[iat];
                let va = &self.p_va[iat];
                let mut b_n_v = false;
                let valence_ok = at.valence == at.chem_bonds_valence || {
                    b_n_v = i32::from(at.valence) + 2 == i32::from(at.chem_bonds_valence);
                    b_n_v
                };
                let e = self.plus_edge(iat);
                if sb_n_neutr < MAX_DIFF_FIXH
                    && at.charge == 0
                    && at.num_h == 0
                    && valence_ok
                    && va.c_metal == 0
                    && va.c_num_valence_electrons == 5
                    && va.c_periodic_row_number == 1
                    && !self
                        .at_mobile_h_revrs()
                        .is_some_and(|a| a[iat].endpoint != 0)
                    && !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                    && !self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0)
                    && !mob.as_ref().is_some_and(|m| m[i] != 0)
                    && self.free_edge(e)
                {
                    let j = if only_n_v && b_n_v {
                        get_charge_flower_upper_edge(self.bns, e)
                    } else {
                        NO_VERTEX
                    };
                    if only_n_v
                        && b_n_v
                        && j != NO_VERTEX
                        && self.bns.edge[j as usize].forbidden == 0
                        && self.bns.edge[j as usize].flow == 0
                    {
                        if num_n_v == 0 {
                            self.curr.num_edges = 0;
                            sb_n_neutr = 0;
                        }
                        sb_n_neutr += 1;
                        num_n_v += 1;
                        if self.add_edge(L::Curr, e) || self.add_edge(L::Curr, j) {
                            return true;
                        }
                    } else if num_n_v == 0 {
                        sb_n_neutr += 1;
                        if self.add_edge(L::Curr, e) {
                            return true;
                        }
                        if b_n_v {
                            let j = get_charge_flower_upper_edge(self.bns, e);
                            if j != NO_VERTEX
                                && self.bns.edge[j as usize].forbidden == 0
                                && self.bns.edge[j as usize].flow == 0
                                && self.add_edge(L::Curr, j)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
            let num_try = sb_n_neutr.min(db_o.len());
            let mut cur_success = 0;
            if num_try != 0 {
                cur_success = self.try_o_db(&db_o, num_try);
            }
            self.curr.num_edges = 0;
            if cur_success != 0 {
                if self.recalc(cur_success) {
                    return true;
                }
                break;
            } else if only_n_v {
                only_n_v = false;
            } else {
                break;
            }
        }
        false
    }

    /// case 13: restored O=AB=N=, original (-)O-AB-N(+)=
    pub fn case_13(&mut self) -> bool {
        if self.c2i.n_num_tg_diff_minus == 0 {
            return false;
        }
        let mob = self.mob_h_inchi();
        let fixed_revrs = self.ps.p_one_inchi[0]
            .as_ref()
            .and_then(|c| c.n_num_h_fixed.clone());
        let mob_revrs: Option<Vec<i8>> = match (&self.ps.p_one_inchi[1], &self.ps.p_one_inchi[0]) {
            (Some(c), _) if !c.n_num_h.is_empty() => Some(c.n_num_h.clone()),
            (_, Some(c)) if !c.n_num_h.is_empty() => Some(c.n_num_h.clone()),
            _ => None,
        };
        let n_tg = self.ps.ti.num_t_groups.min(self.ps.one_ti.num_t_groups);
        for itg in 0..n_tg {
            let (g1, g2) = {
                let a = &self.ps.ti.t_group.as_ref().unwrap()[itg];
                let b = &self.ps.one_ti.t_group.as_ref().unwrap()[itg];
                (*a, *b)
            };
            if !(g1.n_num_endpoints == g2.n_num_endpoints
                && i32::from(g1.num[0]) - i32::from(g1.num[1])
                    == i32::from(g2.num[0]) - i32::from(g2.num[1])
                && g1.num[1] > g2.num[1])
            {
                continue;
            }
            let mut sb_n_neutr = 0usize;
            let mut db_o: Vec<usize> = Vec::new();
            for i in 0..self.ps.num_atoms as usize {
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let e_i = self.minus_edge(i);
                let fixed = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0);
                let mobh = mob.as_ref().is_some_and(|m| m[i] != 0);
                let at = &self.at2[iat];
                if db_o.len() < MAX_DIFF_FIXH
                    && self.p_va[i].c_num_valence_electrons == 6
                    && self
                        .ps
                        .endpoint
                        .as_ref()
                        .is_some_and(|ep| usize::from(ep[i]) == itg + 1)
                    && self.free_edge(e_i)
                    && !fixed
                    && !mobh
                    && !fixed_revrs.as_ref().is_some_and(|f| f[iat] != 0)
                    && !mob_revrs.as_ref().is_some_and(|m| m[iat] != 0)
                    && at.charge == 0
                    && at.num_h == 0
                    && at.valence == 1
                    && at.chem_bonds_valence == 2
                {
                    db_o.push(iat);
                } else {
                    let va = &self.p_va[iat];
                    let e = self.plus_edge(iat);
                    if sb_n_neutr < MAX_DIFF_FIXH
                        && at.charge == 0
                        && at.num_h == 0
                        && (at.valence == 4 && at.chem_bonds_valence == 5)
                        && va.c_metal == 0
                        && va.c_num_valence_electrons == 5
                        && va.c_periodic_row_number >= 1
                        && !self
                            .at_mobile_h_revrs()
                            .is_some_and(|a| a[iat].endpoint != 0)
                        && !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                        && !fixed
                        && !mobh
                        && self.free_edge(e)
                    {
                        sb_n_neutr += 1;
                        if self.add_edge(L::Curr, e) {
                            return true;
                        }
                    }
                }
            }
            let num_try = sb_n_neutr.min(db_o.len());
            let mut cur_success = 0;
            if num_try != 0 {
                cur_success = self.try_o_db(&db_o, num_try);
            }
            self.curr.num_edges = 0;
            if cur_success != 0 {
                if self.recalc(cur_success) {
                    return true;
                }
                break;
            }
        }
        false
    }
}

impl<'a> FixFhCtx<'a> {
    /// the flow edges next to the vertex `v` for case 14: the last edge
    /// with flow that is not forbidden and differs from `skip`
    fn flow_neighbor(&self, v: usize, skip: usize) -> Option<(usize, usize)> {
        for j in (0..self.bns.vert[v].num_adj_edges as usize).rev() {
            let pen = self.vert_edge(v, j);
            if self.bns.edge[pen].flow != 0 && self.bns.edge[pen].forbidden == 0 && pen != skip {
                return Some((pen, (self.bns.edge[pen].neighbor12 as usize) ^ v));
            }
        }
        None
    }

    /// does the last BNS test connect the vertices `a` and `b`
    fn path_between(&self, a: usize, b: usize) -> bool {
        let (a, b) = (a as i32, b as i32);
        (self.res.v_last == a && self.res.v_first == b)
            || (self.res.v_last == b && self.res.v_first == a)
    }

    /// case 14: restored -N=AB=N=CD-XH, original (-)N-AB-N(+)=CD-XH
    pub fn case_14(&mut self) -> bool {
        let c = &self.c2i;
        if !(((c.n_num_tg_inchi <= 1 && c.n_num_rem_h_inchi > c.n_num_rem_h_revrs)
            || !c.c2at.is_empty())
            && b_has_n_v(self.at2, self.ps.num_atoms as usize) != 0)
        {
            return false;
        }
        let mob = self.mob_h_inchi();
        let mut max_success =
            i32::from(self.c2i.n_num_rem_h_inchi) - i32::from(self.c2i.n_num_rem_h_revrs);
        let mut x_list = EdgeList::default();
        let mut n_iii_list = EdgeList::default();
        let mut n_v_array: Vec<i32> = Vec::new();
        let mut cur_success = 0;
        self.ret = 0;
        'body: {
            for i in 0..self.ps.num_atoms as usize {
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let at = &self.at2[iat];
                let va = &self.p_va[iat];
                let endp_i = self.at_mobile_h_revrs().is_some_and(|a| a[i].endpoint != 0);
                let mobh = mob.as_ref().is_some_and(|m| m[i] != 0);
                let fixed_set = self.ps.fixed_h.is_some();
                let fixed = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0);
                let e = self.plus_edge(iat);
                if n_v_array.len() < MAX_DIFF_FIXH
                    && at.chem_bonds_valence == 5
                    && at.valence == 3
                    && at.charge == 0
                    && at.radical == 0
                    && va.c_num_valence_electrons == 5
                    && va.c_periodic_row_number == 1
                    && !endp_i
                    && at.num_h == 0
                    && self.free_edge(e)
                    && self.bns.edge[e as usize].flow != 0
                    && {
                        let j = get_charge_flower_upper_edge(self.bns, e);
                        j != NO_VERTEX
                            && self.bns.edge[j as usize].forbidden == 0
                            && self.bns.edge[j as usize].flow == 0
                    }
                    && !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                    && !mobh
                    && fixed_set
                    && !fixed
                {
                    n_v_array.push(iat as i32);
                } else {
                    let em = self.minus_edge(iat);
                    if at.chem_bonds_valence == 3
                        && at.valence == 2
                        && at.charge == 0
                        && at.radical == 0
                        && va.c_num_valence_electrons == 5
                        && va.c_periodic_row_number == 1
                        && !endp_i
                        && at.num_h == 0
                        && self.free_edge(em)
                        && self.bns.edge[em as usize].flow == 0
                        && !mobh
                        && fixed_set
                        && !fixed
                    {
                        self.ret = add_to_edge_list(&mut n_iii_list, iat as i32, 32);
                        if self.ret != 0 {
                            break 'body;
                        }
                    } else if at.chem_bonds_valence == at.valence
                        && at.charge == 0
                        && at.radical == 0
                        && ((va.c_num_valence_electrons == 5 && va.c_periodic_row_number == 1)
                            || va.c_num_valence_electrons == 6)
                        && at.num_h != 0
                        && self.free_edge(e)
                        && self.bns.edge[e as usize].flow != 0
                        && !mobh
                        && fixed_set
                        && fixed
                    {
                        self.ret = add_to_edge_list(&mut x_list, iat as i32, 32);
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                }
            }
            let num_n_v = n_v_array.len() as i32;
            if max_success == 0 {
                max_success = num_n_v.min(n_iii_list.num_edges).min(x_list.num_edges);
            }
            if num_n_v != 0 && n_iii_list.num_edges != 0 && x_list.num_edges != 0 {
                let mut i1 = 0usize;
                while i1 < n_v_array.len() && cur_success < max_success {
                    let i1c = i1;
                    i1 += 1;
                    let iat_n_v = n_v_array[i1c];
                    if iat_n_v == NO_VERTEX {
                        continue;
                    }
                    let e1_plus = self.plus_edge(iat_n_v as usize);
                    if e1_plus <= 0 {
                        continue;
                    }
                    let e1_flower = get_charge_flower_upper_edge(self.bns, e1_plus);
                    if e1_flower == NO_VERTEX
                        || self.bns.edge[e1_plus as usize].flow != 1
                        || self.bns.edge[e1_flower as usize].flow != 0
                    {
                        continue;
                    }
                    let (e1_plus, e1_flower) = (e1_plus as usize, e1_flower as usize);
                    let mut i2 = n_iii_list.num_edges - 1;
                    while i2 >= 0 && cur_success < max_success {
                        let i2c = i2 as usize;
                        i2 -= 1;
                        let iat_n_iii = n_iii_list.pn_edges[i2c];
                        if iat_n_iii == NO_VERTEX {
                            continue;
                        }
                        let e2_minus = self.minus_edge(iat_n_iii as usize);
                        let e2_plus = self.plus_edge(iat_n_iii as usize);
                        if e2_minus <= 0
                            || e2_plus <= 0
                            || self.bns.edge[e2_minus as usize].flow != 0
                            || self.bns.edge[e2_plus as usize].flow != 1
                        {
                            n_iii_list.pn_edges[i2c] = NO_VERTEX;
                            continue;
                        }
                        let e2_minus = e2_minus as usize;
                        let mut i3 = x_list.num_edges - 1;
                        while i3 >= 0 && cur_success < max_success {
                            let i3c = i3 as usize;
                            i3 -= 1;
                            self.ret = 0;
                            let iat_x = x_list.pn_edges[i3c];
                            if iat_x == NO_VERTEX {
                                continue;
                            }
                            let e3_plus = self.plus_edge(iat_x as usize);
                            if e3_plus <= 0 || self.bns.edge[e3_plus as usize].flow != 1 {
                                x_list.pn_edges[i3c] = NO_VERTEX;
                                continue;
                            }
                            let e3_plus = e3_plus as usize;
                            let pe = e1_flower;
                            let v1 = self.bns.edge[pe].neighbor1 as usize;
                            let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
                            let Some((pe1n, v1n)) = self.flow_neighbor(v1, e1_plus) else {
                                continue;
                            };
                            let Some((pe2n, v2n)) = self.flow_neighbor(v2, e1_plus) else {
                                continue;
                            };
                            // test #1
                            self.bns.edge[pe].flow += 1;
                            self.bns.edge[pe1n].flow -= 1;
                            self.bns.edge[pe2n].flow -= 1;
                            self.bns.vert[v1n].st_edge.flow -= 1;
                            self.bns.vert[v2n].st_edge.flow -= 1;
                            self.bns.tot_st_flow -= 2;
                            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                            set_forbidden_edge_mask(self.bns, &self.other_n_flower, self.mask);
                            self.bns.edge[e2_minus].forbidden &= self.mask_inv;
                            self.bns.edge[e3_plus].forbidden &= self.mask_inv;
                            self.ret =
                                run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                            if self.ret < 0 {
                                break 'body;
                            }
                            if !(self.ret == 1
                                && self.path_between(v1n, v2n)
                                && self.res.delta_charge == 2)
                            {
                                self.ret = 0;
                            }
                            self.bns.edge[e2_minus].forbidden |= self.mask as i8;
                            self.bns.edge[e3_plus].forbidden |= self.mask as i8;
                            self.bns.edge[pe].flow -= 1;
                            self.bns.edge[pe1n].flow += 1;
                            self.bns.edge[pe2n].flow += 1;
                            self.bns.vert[v1n].st_edge.flow += 1;
                            self.bns.vert[v2n].st_edge.flow += 1;
                            self.bns.tot_st_flow += 2;
                            if self.ret == 1 {
                                // test #2
                                self.bns.edge[pe].flow += 1;
                                self.bns.edge[pe1n].flow -= 1;
                                self.bns.edge[pe2n].flow -= 1;
                                self.bns.vert[v1n].st_edge.flow -= 1;
                                self.bns.vert[v2n].st_edge.flow -= 1;
                                self.bns.tot_st_flow -= 2;
                                self.bns.edge[e2_minus].forbidden &= self.mask_inv;
                                self.bns.edge[e1_plus].forbidden &= self.mask_inv;
                                self.ret =
                                    run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                                if self.ret == 1
                                    && self.path_between(v1n, v2n)
                                    && self.res.delta_charge == 2
                                {
                                    self.ret = run_bns_restore_once(self.bns, self.bd);
                                    if self.ret > 0 {
                                        self.num_run_bns += 1;
                                        cur_success += 1;
                                    }
                                }
                                if self.ret <= 0 {
                                    self.bns.edge[pe].flow -= 1;
                                    self.bns.edge[pe1n].flow += 1;
                                    self.bns.edge[pe2n].flow += 1;
                                    self.bns.vert[v1n].st_edge.flow += 1;
                                    self.bns.vert[v2n].st_edge.flow += 1;
                                    self.bns.tot_st_flow += 2;
                                }
                            }
                            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                            remove_forbidden_edge_mask(self.bns, &self.other_n_flower, self.mask);
                            if self.ret > 0 {
                                n_v_array[i1c] = NO_VERTEX;
                                n_iii_list.pn_edges[i2c] = NO_VERTEX;
                                x_list.pn_edges[i3c] = NO_VERTEX;
                            }
                            if self.ret < 0 {
                                break 'body;
                            }
                            if self.ret > 0 {
                                break;
                            }
                        }
                        if self.ret > 0 {
                            break;
                        }
                    }
                }
            }
        }
        remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        remove_forbidden_edge_mask(self.bns, &self.other_n_flower, self.mask);
        self.curr.num_edges = 0;
        if self.ret < 0 {
            return true;
        }
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// the charge-moving attempt of cases 15 and 17: the (+) of the atoms
    /// `db_charged` moves to another atom through two neighbor flows
    fn plus_swap(&mut self, db_charged: &[usize], num_try: i32) -> i32 {
        let mut cur_success = 0;
        set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
        let mut i = 0;
        while i < db_charged.len() && cur_success < num_try {
            let iat = db_charged[i];
            i += 1;
            let pe = self.plus_edge(iat) as usize;
            if self.bns.edge[pe].flow != 0 {
                continue;
            }
            let v1 = self.bns.edge[pe].neighbor1 as usize;
            let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
            let Some((pe1n, v1n)) = self.flow_neighbor(v1, usize::MAX) else {
                continue;
            };
            let Some((pe2n, v2n)) = self.flow_neighbor(v2, usize::MAX) else {
                continue;
            };
            self.bns.edge[pe].flow += 1;
            self.bns.edge[pe1n].flow -= 1;
            self.bns.edge[pe2n].flow -= 1;
            self.bns.vert[v1n].st_edge.flow -= 1;
            self.bns.vert[v2n].st_edge.flow -= 1;
            self.bns.tot_st_flow -= 2;
            self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
            if self.ret == 1
                && self.path_between(v1n, v2n)
                && (self.res.delta_charge == 0 || self.res.delta_charge == 1)
            {
                self.ret = run_bns_restore_once(self.bns, self.bd);
                if self.ret > 0 {
                    self.num_run_bns += 1;
                    cur_success += 1;
                }
            } else {
                self.bns.edge[pe].flow -= 1;
                self.bns.edge[pe1n].flow += 1;
                self.bns.edge[pe2n].flow += 1;
                self.bns.vert[v1n].st_edge.flow += 1;
                self.bns.vert[v2n].st_edge.flow += 1;
                self.bns.tot_st_flow += 2;
            }
        }
        remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        cur_success
    }

    /// case 15: restored -(+)O=AB-N<, original -O-AB=N(+)<
    pub fn case_15(&mut self) -> bool {
        let c = &self.c2i;
        if !(c.n_num_tg_m_revrs > c.n_num_tg_m_inchi
            || c.n_num_rem_h_revrs < c.n_num_rem_h_inchi
            || c.n_num_endp_revrs < c.n_num_endp_inchi
            || (c.n_num_tg_inchi <= 1 && c.n_num_tg_revrs > c.n_num_tg_inchi))
        {
            return false;
        }
        let mob = self.mob_h_inchi();
        let mut sb_neutr: Vec<usize> = Vec::new();
        let mut db_charged: Vec<usize> = Vec::new();
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let fixed = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0);
            let mobh = mob.as_ref().is_some_and(|m| m[i] != 0);
            let e = self.plus_edge(iat);
            if db_charged.len() < MAX_DIFF_FIXH
                && at.charge == 1
                && at.num_h == 0
                && at.valence < at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 6
                && !fixed
                && !mobh
                && self.free_edge(e)
            {
                db_charged.push(iat);
            } else if sb_neutr.len() < MAX_DIFF_FIXH
                && at.charge == 0
                && at.num_h == 0
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && va.c_num_valence_electrons == 5
                && va.c_periodic_row_number == 1
                && !fixed
                && !mobh
                && self.free_edge(e)
            {
                sb_neutr.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let num_try = sb_neutr.len().min(db_charged.len());
        let mut cur_success = 0;
        if num_try != 0 {
            cur_success = self.plus_swap(&db_charged, num_try as i32);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 16: restored O=X-NH(-), original O(-)-X=NH
    pub fn case_16(&mut self) -> bool {
        if self.c2i.n_num_tg_diff_minus == 0 {
            return false;
        }
        let mob = self.mob_h_inchi();
        let n_tg = self.ps.ti.num_t_groups.min(self.ps.one_ti.num_t_groups);
        for itg in 0..n_tg {
            let (g1, g2) = {
                let a = &self.ps.ti.t_group.as_ref().unwrap()[itg];
                let b = &self.ps.one_ti.t_group.as_ref().unwrap()[itg];
                (*a, *b)
            };
            if g1.n_num_endpoints != g2.n_num_endpoints || g1.num[1] >= g2.num[1] {
                continue;
            }
            self.curr.num_edges = 0;
            let mut sb_n_minus: Vec<usize> = Vec::new();
            let mut db_o_neutr = 0usize;
            let k = g2.n_first_endpoint_at_no_pos as usize;
            for j in 0..g2.n_num_endpoints as usize {
                let i = self.ps.one_ti.n_endpoint_atom_number[k + j] as usize;
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let at = &self.at2[iat];
                let va = &self.p_va[iat];
                let e = self.minus_edge(iat);
                let endp = self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0);
                let fixed = self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0);
                let mobh = mob.as_ref().is_some_and(|m| m[i] != 0);
                if db_o_neutr < MAX_DIFF_FIXH
                    && at.charge == 0
                    && at.num_h == 0
                    && at.valence < at.chem_bonds_valence
                    && va.c_metal == 0
                    && va.c_num_valence_electrons == 6
                    && endp
                    && !mobh
                    && self.free_edge(e)
                {
                    db_o_neutr += 1;
                    if self.add_edge(L::Curr, e) {
                        return true;
                    }
                } else if sb_n_minus.len() < MAX_DIFF_FIXH
                    && at.charge == -1
                    && at.num_h != 0
                    && at.valence == at.chem_bonds_valence
                    && va.c_metal == 0
                    && va.c_num_valence_electrons == 5
                    && va.c_periodic_row_number == 1
                    && endp
                    && fixed
                    && !mobh
                    && self.free_edge(e)
                {
                    sb_n_minus.push(iat);
                }
            }
            let num_try = sb_n_minus.len().min(db_o_neutr);
            let mut cur_success = 0;
            if num_try != 0 {
                set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
                let mut i = 0;
                while i < sb_n_minus.len() && (cur_success as usize) < num_try {
                    let iat = sb_n_minus[i];
                    i += 1;
                    let pe = self.minus_edge(iat) as usize;
                    if self.bns.edge[pe].flow == 0 {
                        continue;
                    }
                    match self.try_flow_opt(pe, false, |d| d == 1) {
                        Tried::Success => cur_success += 1,
                        Tried::RolledBack => self.bns.edge[pe].forbidden &= self.mask_inv,
                        Tried::Unrestored => {}
                    }
                }
                remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            }
            self.curr.num_edges = 0;
            if self.recalc(cur_success) {
                return true;
            }
        }
        false
    }

    /// case 17: restored OH(+)=AB-O-, original HO-AB=O(+)-
    pub fn case_17(&mut self) -> bool {
        if self.c2i.n_num_rem_h_inchi >= self.c2i.n_num_rem_h_revrs {
            return false;
        }
        let mob = self.mob_h_inchi();
        let mut sb_neutr: Vec<usize> = Vec::new();
        let mut db_charged: Vec<usize> = Vec::new();
        for i in 0..self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let e = self.plus_edge(iat);
            if db_charged.len() < MAX_DIFF_FIXH
                && at.charge == 1
                && at.num_h != 0
                && at.valence < at.chem_bonds_valence
                && va.c_metal == 0
                && (va.c_num_valence_electrons == 6
                    || (va.c_num_valence_electrons == 5 && va.c_periodic_row_number == 1))
                && self.free_edge(e)
            {
                db_charged.push(iat);
            } else if sb_neutr.len() < MAX_DIFF_FIXH
                && at.charge == 0
                && at.num_h == 0
                && at.valence == at.chem_bonds_valence
                && va.c_metal == 0
                && (va.c_num_valence_electrons == 6
                    || va.c_num_valence_electrons == 7
                    || (va.c_num_valence_electrons == 5 && va.c_periodic_row_number > 1))
                && !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                && !self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0)
                && !mob.as_ref().is_some_and(|m| m[i] != 0)
                && self.free_edge(e)
            {
                sb_neutr.push(iat);
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let mut num_try = sb_neutr.len().min(db_charged.len()) as i32;
        let mut cur_success = 0;
        if num_try != 0 {
            num_try = num_try
                .min(i32::from(self.c2i.n_num_rem_h_revrs) - i32::from(self.c2i.n_num_rem_h_inchi));
            cur_success = self.plus_swap(&db_charged, num_try);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 18: restored -N=AB-X, original -(-)N-AB-X(+)
    pub fn case_18(&mut self) -> bool {
        let c = &self.c2i;
        if !(c.n_num_tg_inchi != 0
            && self.ps.endpoint.is_some()
            && c.n_num_tg_m_inchi > c.n_num_tg_m_revrs
            && c.n_num_endp_inchi > c.n_num_endp_revrs)
        {
            return false;
        }
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        self.ret = 0;
        let n_tg = self.ps.ti.num_t_groups.min(self.ps.one_ti.num_t_groups);
        'body: {
            for itg in 0..n_tg {
                let (g1, g2) = {
                    let a = &self.ps.ti.t_group.as_ref().unwrap()[itg];
                    let b = &self.ps.one_ti.t_group.as_ref().unwrap()[itg];
                    (*a, *b)
                };
                if g1.n_num_endpoints <= g2.n_num_endpoints || g1.num[1] <= g2.num[1] {
                    continue;
                }
                self.curr.num_edges = 0;
                cur_success = 0;
                let k = g1.n_first_endpoint_at_no_pos as usize;
                for j in 0..g1.n_num_endpoints as usize {
                    let i = self.ps.ti.n_endpoint_atom_number[k + j] as usize;
                    let iat = self.ps.n_canon2atno[0][i] as usize;
                    let at = &self.at2[iat];
                    let va = &self.p_va[i];
                    let e1 = self.plus_edge(iat);
                    let e = self.minus_edge(iat);
                    let Some(atm) = self.at_mobile_h_revrs() else {
                        continue;
                    };
                    if self.ps.endpoint.as_ref().is_none_or(|ep| ep[i] == 0)
                        || atm[iat].endpoint != 0
                        || va.c_num_valence_electrons != 5
                        || va.c_periodic_row_number != 1
                        || at.valence != 2
                        || at.num_h != 0
                        || at.radical != 0
                        || (e1 >= 0 && self.bns.edge[e1 as usize].flow == 0)
                        || e < 0
                        || self.bns.edge[e as usize].forbidden != 0
                        || self.bns.edge[e as usize].flow != 0
                    {
                        continue;
                    }
                    if self.add_edge(L::Curr, e) {
                        return true;
                    }
                }
            }
            let max_success = self.curr.num_edges;
            if max_success == 0 {
                break 'body;
            }
            let mut i = 0;
            while i < self.ps.num_atoms as usize && cur_success < max_success {
                let ic = i;
                i += 1;
                let iat = self.ps.n_canon2atno[0][ic] as usize;
                let at = &self.at2[iat];
                let nve = self.p_va[ic].c_num_valence_electrons;
                let e1 = self.minus_edge(iat);
                let e = self.plus_edge(iat);
                if self.ps.endpoint.as_ref().is_some_and(|ep| ep[ic] != 0)
                    || nve == 0
                    || nve == 4
                    || at.num_h != 0
                    || at.radical != 0
                    || (e1 >= 0 && self.bns.edge[e1 as usize].flow == 0)
                    || e < 0
                    || self.bns.edge[e as usize].forbidden != 0
                    || self.bns.edge[e as usize].flow != 1
                {
                    continue;
                }
                set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                set_forbidden_edge_mask(self.bns, &self.other_n_flower, self.mask);
                remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
                let pe = e as usize;
                if self.bns.edge[pe].flow == 0 {
                    continue;
                }
                if self.try_flow_opt(pe, false, |d| d == 1) == Tried::Success {
                    cur_success += 1;
                }
            }
        }
        remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        remove_forbidden_edge_mask(self.bns, &self.other_n_flower, self.mask);
        self.curr.num_edges = 0;
        if self.ret < 0 {
            return true;
        }
        self.recalc(cur_success)
    }

    /// case 19: restored M-OH, original M(-)=OH(+)
    pub fn case_19(&mut self) -> bool {
        if self.c2i.c2at.is_empty() {
            return false;
        }
        let mut cur_success = 0;
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.plus_edge(iat);
            let at = &self.at2[iat];
            if !(a.n_val_electr == 6
                && self.free_edge(e)
                && self.bns.edge[e as usize].flow != 0
                && a.n_fix_h_inchi == 1
                && a.n_mob_h_inchi == 0
                && a.n_fix_h_revrs == 0
                && a.n_mob_h_revrs == 1
                && a.n_at_charge_revrs == 0
                && at.num_h != 0
                && at.valence == 1
                && at.valence == at.chem_bonds_valence)
            {
                continue;
            }
            let metal = at.neighbor[0] as usize;
            let (e_mp, e_mm) = (self.plus_edge(metal), self.minus_edge(metal));
            if self.p_va[metal].c_metal != 0
                && self.free_edge(e_mp)
                && self.free_edge(e_mm)
                && self.bns.edge[self.vert_edge(iat, 0)].forbidden == 0
                && self.add_edge(L::Curr, iat as i32)
            {
                return true;
            }
        }
        if self.curr.num_edges != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            set_forbidden_edge_mask(self.bns, &self.n_flower, self.mask);
            set_forbidden_edge_mask(self.bns, &self.all_bond, self.mask);
            for i in 0..self.curr.num_edges as usize {
                let iat = self.curr.pn_edges[i] as usize;
                let metal = self.at2[iat].neighbor[0] as usize;
                let pe_oh_plus = self.plus_edge(iat) as usize;
                let pe_m_plus = self.plus_edge(metal) as usize;
                let pe_m_minus = self.minus_edge(metal) as usize;
                let pe_om = self.vert_edge(iat, 0);
                self.bns.edge[pe_m_plus].forbidden &= self.mask_inv;
                self.bns.edge[pe_m_minus].forbidden &= self.mask_inv;
                self.bns.edge[pe_om].forbidden &= self.mask_inv;
                let charge_on_metal = (self.bns.edge[pe_m_plus].cap
                    - self.bns.edge[pe_m_plus].flow)
                    - self.bns.edge[pe_m_minus].flow;
                let expected = match charge_on_metal {
                    1 => -1,
                    0 => 1,
                    _ => 0,
                };
                if self.bns.edge[pe_oh_plus].flow == 0 {
                    continue;
                }
                if self.try_flow_opt(pe_oh_plus, false, |d| d == expected) == Tried::Success {
                    cur_success += 1;
                }
                self.bns.edge[pe_m_plus].forbidden |= self.mask as i8;
                self.bns.edge[pe_m_minus].forbidden |= self.mask as i8;
                self.bns.edge[pe_om].forbidden |= self.mask as i8;
            }
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.n_flower, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.all_bond, self.mask);
            self.curr.num_edges = 0;
            return self.recalc(cur_success);
        }
        false
    }
}

impl<'a> FixFhCtx<'a> {
    /// the final (-)-removal attempt of cases 20, 21 and 21a over the single
    /// edges `pes`
    fn minus_remove(&mut self, pes: &[usize], num_try: usize, forbid: bool) -> i32 {
        let mut cur_success = 0;
        set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
        let mut i = 0;
        while i < pes.len() && (cur_success as usize) < num_try {
            let pe = pes[i];
            i += 1;
            if self.bns.edge[pe].flow == 0 {
                continue;
            }
            if self.try_flow_opt(pe, forbid, |d| d == 1) == Tried::Success {
                cur_success += 1;
            }
        }
        remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        cur_success
    }

    /// case 20: restored O(-)-AB=N-, original O=AB-N(-)-
    pub fn case_20(&mut self) -> bool {
        let c = &self.c2i;
        if !(c.c2at.len() > 1 && c.n_num_tg_revrs != 0 && c.n_num_tg_inchi != 0) {
            return false;
        }
        let mut sb_o_minus: Vec<usize> = Vec::new();
        let mut num_db_n = 0usize;
        self.curr.num_edges = 0;
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            if num_db_n < MAX_DIFF_FIXH
                && a.endpt_inchi != 0
                && self.free_edge(e)
                && self.bns.edge[e as usize].flow == 0
                && a.n_fix_h_inchi == 0
                && a.n_mob_h_inchi == 0
                && (self.c2i.n_num_tg_inchi == 1 || a.n_val_electr == 6)
                && a.endpt_revrs == 0
                && a.n_fix_h_revrs == 0
                && a.n_at_charge_revrs == 0
                && at.num_h == 0
                && i32::from(at.valence) + 1 == i32::from(at.chem_bonds_valence)
            {
                num_db_n += 1;
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            } else if sb_o_minus.len() < MAX_DIFF_FIXH
                && a.endpt_inchi == 0
                && self.free_edge(e)
                && self.bns.edge[e as usize].flow == 1
                && a.n_fix_h_inchi == 0
                && a.n_mob_h_inchi == 0
                && a.n_val_electr == 6
                && a.endpt_revrs != 0
                && a.n_fix_h_revrs == 0
                && a.n_mob_h_revrs == 0
                && a.n_at_charge_revrs == -1
                && at.num_h == 0
                && at.valence == 1
                && at.chem_bonds_valence == 1
            {
                sb_o_minus.push(iat);
            }
        }
        if num_db_n == 0 {
            for i in 0..self.ps.num_atoms as usize {
                if !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0) {
                    continue;
                }
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let at = &self.at2[iat];
                let e = self.minus_edge(iat);
                if num_db_n < MAX_DIFF_FIXH
                    && at.charge == 0
                    && at.num_h == 0
                    && i32::from(at.valence) + 1 == i32::from(at.chem_bonds_valence)
                    && self.p_va[iat].c_metal == 0
                    && !self.ps.fixed_h.as_ref().is_some_and(|f| f[i] != 0)
                    && self.free_edge(e)
                    && self.bns.edge[e as usize].flow == 0
                {
                    num_db_n += 1;
                    if self.add_edge(L::Curr, e) {
                        return true;
                    }
                }
            }
        }
        let num_try = sb_o_minus.len().min(num_db_n);
        let mut cur_success = 0;
        if num_try != 0 {
            let pes: Vec<usize> = sb_o_minus
                .iter()
                .map(|&a| self.minus_edge(a) as usize)
                .collect();
            cur_success = self.minus_remove(&pes, num_try, false);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 21: restored R=S(=O)(O(-))O(-), original R-S(=O)(=O)O(-)
    pub fn case_21(&mut self) -> bool {
        let c = &self.c2i;
        if !(!c.c2at.is_empty()
            && c.n_num_tg_revrs != 0
            && c.n_num_tg_h_inchi != 0
            && self.ps.endpoint.is_some())
        {
            return false;
        }
        self.curr.num_edges = 0;
        let mut db_o = 0usize;
        let mut sb_o_minus: Vec<usize> = Vec::new();
        let mut central: Vec<i32> = Vec::new();
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            if db_o < MAX_DIFF_FIXH
                && a.n_val_electr == 6
                && (a.endpt_inchi != 0 || a.n_mob_h_inchi != 0)
                && self.free_edge(e)
                && a.n_fix_h_inchi == 0
                && !(a.endpt_revrs != 0 || a.n_mob_h_revrs != 0)
                && a.n_fix_h_revrs == 0
                && a.n_at_charge_revrs == 0
                && at.num_h == 0
                && at.valence == 1
                && at.chem_bonds_valence == 2
            {
                db_o += 1;
                if self.add_edge(L::Curr, e) {
                    return true;
                }
            }
        }
        let mut i = 0usize;
        while db_o != 0 && i < self.ps.num_atoms as usize {
            let ic = i;
            i += 1;
            let endp = |s: &Self, k: usize| s.ps.endpoint.as_ref().map_or(0, |ep| ep[k]);
            if endp(self, ic) == 0 {
                continue;
            }
            let iat = self.ps.n_canon2atno[0][ic] as usize;
            let at = &self.at2[iat];
            let va = &self.p_va[iat];
            let e = self.minus_edge(iat);
            if !(sb_o_minus.len() < MAX_DIFF_FIXH
                && at.charge == -1
                && at.num_h == 0
                && at.valence == 1
                && at.chem_bonds_valence != 0
                && va.c_metal == 0
                && va.c_num_valence_electrons == 6
                && self
                    .at_mobile_h_revrs()
                    .is_some_and(|a| a[iat].endpoint != 0)
                && !self.ps.fixed_h.as_ref().is_some_and(|f| f[ic] != 0)
                && self.free_edge(e)
                && self.bns.edge[e as usize].flow != 0)
            {
                continue;
            }
            // as in C the neighbor of the atom with the canonical number
            // `ic` (not of `iat`) is taken
            let i_s = self.at2[ic].neighbor[0] as usize;
            if central.contains(&(i_s as i32)) {
                continue;
            }
            let (mut taut_sb, mut taut_db, mut other_db) = (0, 0, 0);
            for j in 0..self.ps.at[i_s].valence as usize {
                let bond_type = i32::from(self.at2[i_s].bond_type[j]);
                let k = self.at2[i_s].neighbor[j] as usize;
                if k == ic {
                    continue;
                }
                if endp(self, k) == endp(self, ic) {
                    taut_sb += i32::from(bond_type == BOND_TYPE_SINGLE);
                    taut_db += i32::from(bond_type == BOND_TYPE_DOUBLE);
                } else if bond_type == BOND_TYPE_DOUBLE {
                    other_db += 1;
                }
            }
            if taut_sb == 0 {
                continue;
            }
            if !(other_db != 0 && taut_db != 0) {
                continue;
            }
            sb_o_minus.push(iat);
            central.push(i_s as i32);
        }
        let num_try = sb_o_minus.len().min(db_o);
        let mut cur_success = 0;
        if num_try != 0 {
            let pes: Vec<usize> = sb_o_minus
                .iter()
                .map(|&a| self.minus_edge(a) as usize)
                .collect();
            cur_success = self.minus_remove(&pes, num_try, true);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 21a: restored R=S-R'=X, original R-S-R'-X(-)
    pub fn case_21a(&mut self) -> bool {
        let c = &self.c2i;
        if !(!c.c2at.is_empty() && c.n_num_tg_revrs != 0 && c.n_num_endp_inchi < c.n_num_endp_revrs)
        {
            return false;
        }
        self.curr.num_edges = 0;
        let mut other_so = EdgeList::default();
        let mut central_s = EdgeList::default();
        let mut so_minus = EdgeList::default();
        let mut cur_success = 0;
        let atm: Option<Vec<u16>> = self
            .at_mobile_h_revrs()
            .map(|a| a.iter().map(|x| x.endpoint).collect());
        'body: {
            let Some(atm) = atm else { break 'body };
            for ci in 0..self.c2i.c2at.len() {
                let a = self.c2i.c2at[ci];
                let iat = a.atom_number as usize;
                let e = self.minus_edge(iat);
                let at = &self.at2[iat];
                if !(a.endpt_inchi == 0
                    && self.free_edge(e)
                    && a.n_fix_h_inchi == 0
                    && (a.endpt_revrs != 0 || a.n_mob_h_revrs != 0)
                    && a.n_fix_h_revrs == 0
                    && at.num_h == 0)
                {
                    continue;
                }
                let e = e as usize;
                let i_s = at.neighbor[0] as usize;
                if self.p_va[iat].c_num_valence_electrons == 6
                    && at.charge == -1
                    && self.bns.edge[e].flow != 0
                    && at.valence == 1
                    && at.chem_bonds_valence == 1
                    && self.p_va[i_s].c_num_valence_electrons == 6
                    && self.p_va[i_s].c_periodic_row_number > 1
                    && self.at2[i_s].valence >= 4
                {
                    if find_in_edge_list(&central_s, i_s as i32) >= 0 {
                        continue;
                    }
                    let (mut taut_sb, mut taut_db, mut other_db, mut num_endp_o) = (0, 0, 0, 0);
                    let val = self.ps.at[i_s].valence as usize;
                    for j in 0..val {
                        let bond_type = i32::from(self.at2[i_s].bond_type[j]);
                        let k = self.at2[i_s].neighbor[j] as usize;
                        if k == iat {
                            continue;
                        }
                        if a.endpt_revrs == atm[k] && self.at2[k].endpoint == 0 {
                            taut_sb += i32::from(bond_type == BOND_TYPE_SINGLE);
                            taut_db += i32::from(bond_type == BOND_TYPE_DOUBLE);
                            num_endp_o += i32::from(
                                self.p_va[k].c_num_valence_electrons == 6
                                    && self.at2[k].valence == 1,
                            );
                        } else if bond_type == BOND_TYPE_DOUBLE {
                            other_db += 1;
                        }
                    }
                    if num_endp_o == 0 {
                        continue;
                    }
                    if taut_sb + taut_db + other_db <= num_endp_o {
                        continue;
                    }
                    for j in 0..val {
                        let bond_type = i32::from(self.at2[i_s].bond_type[j]);
                        let k = self.at2[i_s].neighbor[j] as usize;
                        let ek = self.minus_edge(k);
                        if a.endpt_revrs == atm[k]
                            && self.at2[k].endpoint == 0
                            && self.p_va[k].c_num_valence_electrons == 6
                            && self.at2[k].valence == 1
                            && self.free_edge(ek)
                        {
                            let ek_us = ek as usize;
                            if bond_type == BOND_TYPE_DOUBLE
                                && self.at2[k].charge == 0
                                && self.bns.edge[ek_us].flow == 0
                            {
                                self.ret = add_to_edge_list(&mut other_so, ek, INC_ADD_EDGE);
                                if self.ret != 0 {
                                    break 'body;
                                }
                            } else if bond_type == BOND_TYPE_SINGLE
                                && self.at2[k].charge == -1
                                && self.bns.edge[ek_us].flow != 0
                            {
                                self.ret = add_to_edge_list(&mut so_minus, ek, INC_ADD_EDGE);
                                if self.ret != 0 {
                                    break 'body;
                                }
                            }
                        }
                    }
                    self.ret = add_to_edge_list(&mut central_s, i_s as i32, INC_ADD_EDGE);
                    if self.ret != 0 {
                        break 'body;
                    }
                } else if at.charge == 0
                    && self.bns.edge[e].flow == 0
                    && i32::from(at.valence) + 1 == i32::from(at.chem_bonds_valence)
                    && self.add_edge(L::Curr, e as i32)
                {
                    return true;
                }
            }
            for i in 0..other_so.num_edges as usize {
                remove_from_edge_list_by_value(&mut self.curr, other_so.pn_edges[i]);
            }
            let num_try = (so_minus.num_edges.min(self.curr.num_edges)) as usize;
            if num_try != 0 {
                let pes: Vec<usize> = (0..so_minus.num_edges as usize)
                    .map(|i| so_minus.pn_edges[i] as usize)
                    .collect();
                cur_success = self.minus_remove(&pes, num_try, false);
            }
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// case 22: restored N(-)=N(+)=C...=O, original N#N-N=...-O(-)
    pub fn case_22(&mut self) -> bool {
        if self.c2i.c2at.is_empty() {
            return false;
        }
        self.curr.num_edges = 0;
        let mut db_o: Vec<usize> = Vec::new();
        for ci in 0..self.c2i.c2at.len() {
            let a = self.c2i.c2at[ci];
            let iat = a.atom_number as usize;
            let e = self.minus_edge(iat);
            let at = &self.at2[iat];
            if db_o.len() < MAX_DIFF_FIXH
                && a.n_val_electr == 6
                && (a.endpt_inchi != 0 || a.n_mob_h_inchi != 0)
                && self.free_edge(e)
                && a.n_fix_h_inchi == 0
                && !(a.endpt_revrs != 0 || a.n_mob_h_revrs != 0)
                && a.n_fix_h_revrs == 0
                && a.n_at_charge_revrs == 0
                && at.num_h == 0
                && at.valence == 1
                && at.chem_bonds_valence == 2
            {
                db_o.push(iat);
            }
        }
        let mut i = 0;
        while !db_o.is_empty() && i < self.ps.num_atoms as usize {
            let iat = self.ps.n_canon2atno[0][i] as usize;
            i += 1;
            let at = &self.at2[iat];
            let e = self.minus_edge(iat);
            if !(at.charge == -1
                && at.num_h == 0
                && at.valence == 1
                && at.chem_bonds_valence == 2
                && self.p_va[iat].c_metal == 0
                && self.p_va[iat].c_num_valence_electrons == 5
                && self.free_edge(e)
                && self.bns.edge[e as usize].flow != 0
                && !self
                    .at_mobile_h_revrs()
                    .is_some_and(|a| a[iat].endpoint != 0))
            {
                continue;
            }
            let i_n2 = at.neighbor[0] as usize;
            let an = &self.at2[i_n2];
            if self.p_va[i_n2].c_num_valence_electrons == 5
                && i32::from(at.bond_type[0]) == BOND_TYPE_DOUBLE
                && an.charge == 1
                && an.valence == 2
                && an.chem_bonds_valence == 4
                && {
                    let i_c = an.neighbor[usize::from(an.neighbor[0] as usize == i_n2)] as usize;
                    self.p_va[i_c].c_num_valence_electrons == 4
                }
                && self.add_edge(L::Curr, e)
            {
                return true;
            }
        }
        let num_try = (self.curr.num_edges as usize).min(db_o.len());
        let mut cur_success = 0;
        if num_try != 0 {
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
            let mut i = 0;
            while i < db_o.len() && (cur_success as usize) < num_try {
                let iat = db_o[i];
                i += 1;
                let pe_db_o_minus = self.minus_edge(iat) as usize;
                let pe = self.vert_edge(iat, 0);
                if self.bns.edge[pe].flow == 0 {
                    continue;
                }
                self.bns.edge[pe_db_o_minus].forbidden &= self.mask_inv;
                if self.try_flow(pe, |d| d == 0) == Tried::Success {
                    cur_success += 1;
                }
                self.bns.edge[pe].forbidden &= self.mask_inv;
                self.bns.edge[pe_db_o_minus].forbidden |= self.mask as i8;
            }
            remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// the O/other neighbor count of the N(V) atom `i_n` in case 23 (`oh`:
    /// searching for -OH, otherwise for =O)
    fn count_no2_neighbors(&self, i_n: usize, iat: usize, oh: bool) -> (i32, i32) {
        let (mut num_o, mut num_others) = (0, 0);
        let atm = self.at_mobile_h_revrs();
        for k in 0..self.at2[i_n].valence as usize {
            let neigh = self.at2[i_n].neighbor[k] as usize;
            if neigh == iat {
                continue;
            }
            let an = &self.at2[neigh];
            let endp = self.ps.endpoint.as_ref().is_some_and(|ep| ep[neigh] != 0);
            let o_like = self.p_va[neigh].c_num_valence_electrons == 6
                && endp
                && an.valence == 1
                && an.radical == 0
                && if oh {
                    an.num_h == 1 && an.charge == 0
                } else {
                    !atm.is_some_and(|a| a[neigh].endpoint != 0)
                        && an.num_h == 0
                        && (an.charge == 0 || an.charge == -1)
                        && i32::from(an.chem_bonds_valence) - i32::from(an.charge) == 2
                };
            if o_like {
                num_o += 1;
            } else if oh {
                if i32::from(self.at2[i_n].bond_type[k]) == BOND_TYPE_DOUBLE
                    && an.valence >= 2
                    && an.valence < an.chem_bonds_valence
                {
                    num_others += 1;
                }
            } else if i32::from(self.at2[i_n].bond_type[k]) == BOND_TYPE_SINGLE
                && an.valence > 1
                && an.valence < an.chem_bonds_valence
            {
                num_others += 1;
            }
        }
        (num_o, num_others)
    }

    /// case 23: -NO2 that is to be tautomeric but is not
    pub fn case_23(&mut self) -> bool {
        const NOOH: usize = 0;
        const WRONG_TAUT: usize = 1;
        const TAUT: usize = 2;
        const O_FIXED: usize = 3;
        let c = &self.c2i;
        if !(!c.c2at.is_empty() && c.n_num_tg_inchi == 1) {
            return false;
        }
        let fixed_revrs = self.ps.p_one_inchi[0]
            .as_ref()
            .and_then(|c| c.n_num_h_fixed.clone());
        let mob_revrs: Option<Vec<i8>> = match (&self.ps.p_one_inchi[1], &self.ps.p_one_inchi[0]) {
            (Some(c), _) if !c.n_num_h.is_empty() => Some(c.n_num_h.clone()),
            (_, Some(c)) if !c.n_num_h.is_empty() => Some(c.n_num_h.clone()),
            _ => None,
        };
        let mut changeable: [EdgeList; 4] = Default::default();
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        let mut db_o: Vec<usize> = Vec::new();
        let mut no2: Vec<usize> = Vec::new();
        'body: {
            for ci in 0..self.c2i.c2at.len() {
                let a = self.c2i.c2at[ci];
                let iat = a.atom_number as usize;
                let e = self.minus_edge(iat);
                let at = &self.at2[iat];
                if !(db_o.len() < MAX_DIFF_FIXH
                    && a.n_val_electr == 6
                    && a.endpt_inchi != 0
                    && self.free_edge(e)
                    && a.n_fix_h_inchi == 0
                    && a.endpt_revrs == 0
                    && a.n_fix_h_revrs == 0
                    && a.n_at_charge_revrs == 0
                    && at.num_h == 0
                    && at.valence == 1
                    && at.chem_bonds_valence == 2)
                {
                    continue;
                }
                let i_n = at.neighbor[0] as usize;
                let an = &self.at2[i_n];
                if !(self.p_va[i_n].c_num_valence_electrons == 5
                    && an.valence == 3
                    && (an.charge == 0 || an.charge == 1)
                    && i32::from(an.chem_bonds_valence) == 5 - i32::from(an.charge))
                {
                    continue;
                }
                let (num_o, num_others) = self.count_no2_neighbors(i_n, iat, false);
                if num_o != 1 || num_others != 1 {
                    continue;
                }
                if !no2.contains(&i_n) {
                    no2.push(i_n);
                    db_o.push(iat);
                }
                self.ret = add_to_edge_list(&mut changeable[O_FIXED], e, INC_ADD_EDGE);
                if self.ret != 0 {
                    break 'body;
                }
            }
            if !db_o.is_empty() {
                let n_atoms = self.ps.num_atoms as usize;
                // 1. =N(=O)-OH
                for i in 0..n_atoms {
                    let iat = self.ps.n_canon2atno[0][i] as usize;
                    let at = &self.at2[iat];
                    if !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                        || self.p_va[i].c_num_valence_electrons != 6
                        || at.valence != 1
                        || at.charge != 0
                        || self.minus_edge(iat) < 0
                        || i32::from(at.num_h) + i32::from(at.chem_bonds_valence) != 2
                    {
                        continue;
                    }
                    let i_n = at.neighbor[0] as usize;
                    let e = self.plus_edge(i_n);
                    let an = &self.at2[i_n];
                    if self.p_va[i_n].c_num_valence_electrons != 5
                        || e < 0
                        || self.bns.edge[e as usize].forbidden != 0
                        || self.bns.edge[e as usize].flow == 0
                        || an.charge != 0
                        || an.valence != 3
                        || an.chem_bonds_valence != 5
                    {
                        continue;
                    }
                    let (num_o, num_others) = self.count_no2_neighbors(i_n, iat, true);
                    if num_o != 1 || num_others != 1 {
                        continue;
                    }
                    self.ret = add_to_edge_list(&mut changeable[NOOH], e, INC_ADD_EDGE);
                    if self.ret == 0 {
                        self.ret = add_to_edge_list(&mut changeable[O_FIXED], e, INC_ADD_EDGE);
                    }
                    if self.ret != 0 {
                        break 'body;
                    }
                    let j = get_charge_flower_upper_edge(self.bns, e);
                    if j != NO_VERTEX {
                        self.ret = add_to_edge_list(&mut changeable[NOOH], j, INC_ADD_EDGE);
                        if self.ret == 0 {
                            self.ret = add_to_edge_list(&mut changeable[O_FIXED], e, INC_ADD_EDGE);
                        }
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                }
                // 2. (-) atoms that are tautomeric but should not be
                for i in 0..n_atoms {
                    let iat = self.ps.n_canon2atno[0][i] as usize;
                    let at = &self.at2[iat];
                    let endp = self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0);
                    let e = self.minus_edge(iat);
                    let wrong = at.charge == -1
                        && !endp
                        && self
                            .at_mobile_h_revrs()
                            .is_some_and(|a| a[i].endpoint != 0 || at.num_h < a[i].num_h);
                    let wrong2 = !wrong
                        && at.charge == -1
                        && endp
                        && !self.at_mobile_h_revrs().is_some_and(|a| a[i].endpoint != 0)
                        && fixed_revrs.as_ref().is_some_and(|f| f[i] == -1)
                        && mob_revrs.as_ref().is_some_and(|m| m[i] == 1)
                        && self.ps.fixed_h.as_ref().is_some_and(|f| f[i] == 0);
                    if (wrong || wrong2)
                        && e >= 0
                        && find_in_edge_list(&changeable[O_FIXED], e) < 0
                        && self.bns.edge[e as usize].forbidden == 0
                        && self.bns.edge[e as usize].flow != 0
                    {
                        self.ret = add_to_edge_list(&mut changeable[WRONG_TAUT], e, INC_ADD_EDGE);
                        if self.ret == 0 {
                            self.ret = add_to_edge_list(&mut changeable[O_FIXED], e, INC_ADD_EDGE);
                        }
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                }
                // 3. (-) atoms that are tautomeric
                for i in 0..n_atoms {
                    let iat = self.ps.n_canon2atno[0][i] as usize;
                    let e = self.minus_edge(iat);
                    if self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0)
                        && self.at_mobile_h_revrs().is_some_and(|a| a[i].endpoint != 0)
                        && self.at2[iat].charge == -1
                        && e >= 0
                        && self.bns.edge[e as usize].forbidden == 0
                        && self.bns.edge[e as usize].flow != 0
                        && find_in_edge_list(&changeable[O_FIXED], e) < 0
                    {
                        self.ret = add_to_edge_list(&mut changeable[TAUT], e, INC_ADD_EDGE);
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                }
                // finally, try to move charges from O=N
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
                    while one_success == 0 && k <= TAUT {
                        let kk = k;
                        k += 1;
                        if changeable[kk].num_edges == 0 {
                            continue;
                        }
                        let expected = if kk == NOOH { 2 } else { 0 };
                        set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                        remove_forbidden_edge_mask(self.bns, &changeable[kk], self.mask);
                        self.bns.edge[pe_db_o_minus].forbidden &= self.mask_inv;
                        self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                        if self.ret == 1
                            && self.path_between(v1, v2)
                            && self.res.delta_charge == expected
                        {
                            self.ret = run_bns_restore_once(self.bns, self.bd);
                            if self.ret > 0 {
                                self.num_run_bns += 1;
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
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// the finishing loop of cases 24 and 25: the flow of the edge `pe` is
    /// removed and the BNS test is repeated for the sets
    /// `changeable[0..=last]` until it succeeds with the expected charge
    /// change; `pef` is a flower edge to be unlocked
    fn try_changeable(
        &mut self,
        changeable: &[EdgeList],
        last: usize,
        pe: usize,
        pef: Option<usize>,
        expected: i32,
    ) -> i32 {
        let v1 = self.bns.edge[pe].neighbor1 as usize;
        let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
        self.bns.edge[pe].flow -= 1;
        self.bns.vert[v1].st_edge.flow -= 1;
        self.bns.vert[v2].st_edge.flow -= 1;
        self.bns.tot_st_flow -= 2;
        let mut one_success = 0;
        let mut k = 0;
        while one_success == 0 && k <= last {
            let kk = k;
            k += 1;
            if changeable[kk].num_edges == 0 {
                continue;
            }
            set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
            remove_forbidden_edge_mask(self.bns, &changeable[kk], self.mask);
            if let Some(pef) = pef {
                self.bns.edge[pef].forbidden &= self.mask_inv;
            }
            self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
            if self.ret == 1 && self.path_between(v1, v2) && self.res.delta_charge == expected {
                self.ret = run_bns_restore_once(self.bns, self.bd);
                if self.ret > 0 {
                    self.num_run_bns += 1;
                    one_success += 1;
                }
            }
        }
        remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
        if one_success == 0 {
            self.bns.edge[pe].flow += 1;
            self.bns.vert[v1].st_edge.flow += 1;
            self.bns.vert[v2].st_edge.flow += 1;
            self.bns.tot_st_flow += 2;
        }
        one_success
    }

    /// case 24: InChI normalization -N(-)-N(+)(IV) => -N=N(V) prevents
    /// tautomerism on -N(-)- in case of ADP
    pub fn case_24(&mut self) -> bool {
        const MISSED_TAUT: usize = 0;
        const OTHER_TAUT_O: usize = 1;
        const OTHER_TAUT_N: usize = 2;
        const NN: usize = 3;
        const AVOID: usize = 4;
        let c = &self.c2i;
        if !(!c.c2at.is_empty() && c.n_num_tg_inchi == 1) {
            return false;
        }
        let mut ch: [EdgeList; 5] = Default::default();
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        let atf: Option<Vec<InpAtom>> = self.atf_mobile_h_revrs().map(<[InpAtom]>::to_vec);
        'body: {
            for ci in 0..self.c2i.c2at.len() {
                let a = self.c2i.c2at[ci];
                let iat = a.atom_number as usize;
                let e1 = self.minus_edge(iat);
                let at = &self.at2[iat];
                // -N=N(V)
                let mut done = false;
                if a.n_val_electr == 5
                    && a.endpt_inchi != 0
                    && self.free_edge(e1)
                    && a.n_fix_h_inchi == 0
                    && a.n_mob_h_inchi == 0
                    && a.endpt_revrs == 0
                    && a.n_fix_h_revrs == 0
                    && a.n_at_charge_revrs == 0
                    && at.num_h == 0
                    && at.valence == 2
                    && at.chem_bonds_valence == 3
                {
                    let i_n = at.neighbor
                        [usize::from(i32::from(at.bond_type[0]) != BOND_TYPE_DOUBLE)]
                        as usize;
                    let an = &self.at2[i_n];
                    let e = self.plus_edge(i_n);
                    if self.p_va[i_n].c_num_valence_electrons == 5
                        && an.chem_bonds_valence == 5
                        && an.charge == 0
                        && an.num_h == 0
                        && an.radical == 0
                        && self.free_edge(e)
                        && self.bns.edge[e as usize].flow != 0
                        && find_in_edge_list(&ch[AVOID], e) < 0
                    {
                        done = true;
                        let ef = get_charge_flower_upper_edge(self.bns, e);
                        for (set, val) in [
                            (NN, e),
                            (NN, ef),
                            (NN, 1),
                            (AVOID, e1),
                            (AVOID, e),
                            (AVOID, ef),
                        ] {
                            self.ret = add_to_edge_list(&mut ch[set], val, INC_ADD_EDGE);
                            if self.ret != 0 {
                                break 'body;
                            }
                        }
                        // mark -N= so that (-) will not be moved to it; as in
                        // C the index is not decremented
                        let em = self.p_va[iat].n_c_minus_group_edge;
                        if em >= 0
                            && self
                                .bns
                                .edge
                                .get(em as usize)
                                .is_some_and(|x| x.forbidden == 0)
                            && find_in_edge_list(&ch[AVOID], em) < 0
                        {
                            self.ret = add_to_edge_list(&mut ch[AVOID], em, INC_ADD_EDGE);
                            if self.ret != 0 {
                                break 'body;
                            }
                        }
                    }
                }
                if done {
                    continue;
                }
                // -N(-)N(IV)(+)
                let Some(atf) = atf.as_ref() else { continue };
                let e = self.minus_edge(iat);
                let at = &self.at2[iat];
                if a.n_val_electr == 5
                    && a.endpt_inchi != 0
                    && self.free_edge(e)
                    && a.n_fix_h_inchi == 0
                    && a.n_mob_h_inchi == 0
                    && a.endpt_revrs == 0
                    && a.n_fix_h_revrs == 0
                    && a.n_at_charge_revrs == -1
                    && at.num_h == 0
                    && at.valence == 2
                    && at.chem_bonds_valence == 2
                    && atf[iat].valence == 2
                    && atf[iat].chem_bonds_valence == 3
                {
                    let i_n = atf[iat].neighbor
                        [usize::from(i32::from(atf[iat].bond_type[0]) != BOND_TYPE_DOUBLE)]
                        as usize;
                    let an = &self.at2[i_n];
                    if self.p_va[i_n].c_num_valence_electrons == 5
                        && an.charge == 1
                        && an.chem_bonds_valence == 4
                        && atf[i_n].charge == 0
                        && atf[i_n].chem_bonds_valence == 5
                        && an.num_h == 0
                        && an.radical == 0
                        && self.bns.edge[e as usize].flow != 0
                        && find_in_edge_list(&ch[AVOID], e) < 0
                    {
                        for (set, val) in [(NN, e), (NN, NO_VERTEX), (NN, 1), (AVOID, e)] {
                            self.ret = add_to_edge_list(&mut ch[set], val, INC_ADD_EDGE);
                            if self.ret != 0 {
                                break 'body;
                            }
                        }
                    }
                }
            }
            if ch[NN].num_edges == 0 {
                break 'body;
            }
            // collect all relevant tautomeric atoms
            for i in 0..self.ps.num_atoms as usize {
                if !self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0) {
                    continue;
                }
                let iat = self.ps.n_canon2atno[0][i] as usize;
                let at = &self.at2[iat];
                if at.charge != 0 || at.radical != 0 || at.valence == at.chem_bonds_valence {
                    continue;
                }
                let e = self.minus_edge(iat);
                if e < 0
                    || self.bns.edge[e as usize].forbidden != 0
                    || self.bns.edge[e as usize].flow != 0
                {
                    continue;
                }
                if find_in_edge_list(&ch[AVOID], e) >= 0 {
                    continue;
                }
                let set = if !self
                    .at_mobile_h_revrs()
                    .is_some_and(|a| a[iat].endpoint != 0)
                {
                    Some(MISSED_TAUT)
                } else if self.p_va[iat].c_num_valence_electrons == 6 {
                    Some(OTHER_TAUT_O)
                } else if self.p_va[iat].c_num_valence_electrons == 5 {
                    Some(OTHER_TAUT_N)
                } else {
                    None
                };
                if let Some(set) = set {
                    // the AVOID test is repeated as in C
                    if find_in_edge_list(&ch[AVOID], e) < 0 {
                        self.ret = add_to_edge_list(&mut ch[set], e, INC_ADD_EDGE);
                        if self.ret == 0 {
                            self.ret = add_to_edge_list(&mut ch[AVOID], e, INC_ADD_EDGE);
                        }
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                }
            }
            // finally, try to move charges from -N(-)-N(+) or to N(V)
            let mut i = 0;
            while i < ch[NN].num_edges as usize {
                let pe = ch[NN].pn_edges[i];
                let pef = ch[NN].pn_edges[i + 1];
                let expected = ch[NN].pn_edges[i + 2];
                i += 3;
                if self.bns.edge[pe as usize].flow == 0 {
                    continue;
                }
                let pef = if pef != NO_VERTEX {
                    Some(pef as usize)
                } else {
                    None
                };
                let chs = [ch[0].clone(), ch[1].clone(), ch[2].clone()];
                cur_success += self.try_changeable(&chs, OTHER_TAUT_N, pe as usize, pef, expected);
            }
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }

    /// case 25: the restored InChI does not have 2 or more added protons;
    /// possibly tautomeric endpoints are missing
    pub fn case_25(&mut self) -> bool {
        const MISSED_TAUT_1: usize = 0;
        const MISSED_TAUT_ALL: usize = 1;
        const OTHER_TAUT_1: usize = 2;
        const OTHER_TAUT_ALL: usize = 3;
        const NO_IN_NO2M2: usize = 4;
        const AVOID: usize = 5;
        let c = &self.c2i;
        if !(!c.c2at.is_empty()
            && c.n_num_tg_inchi == 1
            && c.n_num_rem_h_revrs > c.n_num_rem_h_inchi
            && c.n_num_rem_h_inchi < 0
            && (c.n_num_endp_revrs < c.n_num_endp_inchi || c.n_num_tg_revrs > c.n_num_tg_inchi))
        {
            return false;
        }
        let tg_inchi = self.c2i.n_num_tg_inchi;
        let mut ch: [EdgeList; 6] = Default::default();
        self.curr.num_edges = 0;
        let mut cur_success = 0;
        'body: {
            for i in 0..self.ps.num_atoms as usize {
                let iat = self.ps.n_canon2atno[0][i] as usize;
                if self.ps.endpoint.as_ref().is_some_and(|ep| ep[i] != 0) {
                    let e = self.minus_edge(iat);
                    if e < 0
                        || self.bns.edge[e as usize].forbidden != 0
                        || find_in_edge_list(&ch[AVOID], e) >= 0
                    {
                        continue;
                    }
                    let nve = self.p_va[iat].c_num_valence_electrons;
                    let b_first = (nve == 5 && tg_inchi == 1) || (nve == 6 && tg_inchi != 1);
                    let adds: &mut [(usize, bool)] = &mut [
                        (MISSED_TAUT_1, false),
                        (MISSED_TAUT_ALL, false),
                        (OTHER_TAUT_1, false),
                        (OTHER_TAUT_ALL, true),
                        (AVOID, true),
                    ];
                    let missed = !self.at_mobile_h_revrs().is_some_and(|a| a[i].endpoint != 0);
                    for &mut (set, always) in adds.iter_mut() {
                        let want = match set {
                            MISSED_TAUT_1 => missed && b_first,
                            MISSED_TAUT_ALL => missed,
                            OTHER_TAUT_1 => b_first,
                            _ => always,
                        };
                        if want {
                            self.ret = add_to_edge_list(&mut ch[set], e, INC_ADD_EDGE);
                            if self.ret != 0 {
                                break 'body;
                            }
                        }
                    }
                } else {
                    let at = &self.at2[iat];
                    if !(at.valence == 1
                        && at.charge == -1
                        && self.p_va[iat].c_num_valence_electrons == 6)
                    {
                        continue;
                    }
                    let i_n = at.neighbor[0] as usize;
                    let an = &self.at2[i_n];
                    if self.p_va[i_n].c_num_valence_electrons != 5 {
                        continue;
                    }
                    let canon_n = self.ps.n_atno2canon[0][i_n] as usize;
                    let e = self.plus_edge(i_n);
                    let e1 = self.minus_edge(iat);
                    if !(!self.ps.endpoint.as_ref().is_some_and(|ep| ep[canon_n] != 0)
                        && an.valence == 3
                        && an.chem_bonds_valence == 3
                        && an.charge == 0
                        && an.radical == 0
                        && self.free_edge(e)
                        && self.bns.edge[e as usize].flow != 0
                        && self.free_edge(e1)
                        && self.bns.edge[e1 as usize].flow != 0
                        && find_in_edge_list(&ch[AVOID], e) < 0
                        && find_in_edge_list(&ch[AVOID], e1) < 0)
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
                            && !self.ps.endpoint.as_ref().is_some_and(|ep| ep[neigh] != 0)
                            && ag.valence == 1
                            && ag.num_h == 0
                            && ag.radical == 0
                            && ag.charge == -1
                            && ag.chem_bonds_valence == 1
                        {
                            num_o += 1;
                        } else if i32::from(an.bond_type[k]) == BOND_TYPE_SINGLE
                            && ag.valence > 1
                            && ag.valence < ag.chem_bonds_valence
                        {
                            num_others += 1;
                        }
                    }
                    if num_o != 1 && num_others != 1 {
                        continue;
                    }
                    for (set, val) in [(NO_IN_NO2M2, e), (NO_IN_NO2M2, e1), (AVOID, e), (AVOID, e1)]
                    {
                        self.ret = add_to_edge_list(&mut ch[set], val, INC_ADD_EDGE);
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                }
            }
            if ch[NO_IN_NO2M2].num_edges == 0 || ch[OTHER_TAUT_ALL].num_edges == 0 {
                break 'body;
            }
            let mut i = 0;
            while i < ch[NO_IN_NO2M2].num_edges as usize {
                let pe = ch[NO_IN_NO2M2].pn_edges[i + 1] as usize;
                let pef = ch[NO_IN_NO2M2].pn_edges[i] as usize;
                i += 2;
                if self.bns.edge[pe].flow == 0 {
                    continue;
                }
                let chs = [ch[0].clone(), ch[1].clone(), ch[2].clone(), ch[3].clone()];
                cur_success += self.try_changeable(&chs, OTHER_TAUT_ALL, pe, Some(pef), 3);
            }
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

impl<'a> FixFhCtx<'a> {
    /// case 06e: restored XHn(+)=-AB-YHm, original XHn-=AB=YHm(+); move (+)
    /// between the Mobile-H groups
    pub fn case_06e(&mut self) -> bool {
        let norm = self.ps.p_one_norm_data[1].as_ref();
        if !(self.c2i.n_num_tg_inchi > 1
            && (self.ps.n_num_removed_protons_mobh_inchi > 0
                || self.ps.ti.tni.n_num_removed_protons > 0)
            && self.ps.fixed_h.is_some()
            && self.ps.endpoint.is_some()
            && norm.is_some_and(|d| d.at_fixed_bonds.is_some()))
        {
            return false;
        }
        let atf: Option<Vec<InpAtom>> = norm.and_then(|d| {
            d.at_fixed_bonds.clone().or_else(|| {
                if d.at.is_empty() {
                    None
                } else {
                    Some(d.at.clone())
                }
            })
        });
        let mut cur_success = 0;
        let mut cur_charge = EdgeList::default();
        let mut endp_list = EdgeList::default();
        self.curr.num_edges = 0;
        'body: {
            let Some(atf) = atf else { break 'body };
            let mut tdhc: Vec<TgDiffHChgFh> = Vec::new();
            let n_wrong = fill_tg_diff_h_chg_fh(
                &mut tdhc,
                MAX_DIFF_FIXH,
                self.at2,
                &atf,
                &self.ps.n_canon2atno[0],
                self.p_va,
                &self.ps.ti,
                &mut endp_list,
            );
            if n_wrong < 1 {
                break 'body;
            }
            let n_wrong = n_wrong as usize;
            let (mut num_2_remove_plus, mut num_2_add_plus) = (0i32, 0i32);
            for t in tdhc.iter().take(n_wrong) {
                if t.n_num_h_inchi > t.n_num_h_norml
                    && t.n_num_p_revrs > t.n_num_p_norml
                    && t.n[F_NUM_R_POS_CHG_H] != 0
                {
                    let needed = i32::from(t.n_num_h_inchi - t.n_num_h_norml)
                        .min(i32::from(t.n[F_NUM_R_POS_CHG_H]));
                    num_2_remove_plus += needed;
                    let off = t.i[F_NUM_R_POS_CHG_H] as usize;
                    for jj in 0..t.n[F_NUM_R_POS_CHG_H] as usize {
                        let iat = endp_list.pn_edges[off + jj] as usize;
                        let e = self.plus_edge(iat);
                        self.ret = add_to_edge_list(&mut cur_charge, e, INC_ADD_EDGE);
                        if self.ret != 0 {
                            break 'body;
                        }
                    }
                } else if t.n_num_h_inchi < t.n_num_h_norml && t.n[F_NUM_R_NEUTRL_H] != 0 {
                    let needed = i32::from(t.n_num_h_norml - t.n_num_h_inchi)
                        .min(i32::from(t.n[F_NUM_R_NEUTRL_H]));
                    num_2_add_plus += needed;
                    let off = t.i[F_NUM_R_NEUTRL_H] as usize;
                    for jj in 0..t.n[F_NUM_R_NEUTRL_H] as usize {
                        let iat = endp_list.pn_edges[off + jj] as usize;
                        let e = self.plus_edge(iat);
                        if self.add_edge(L::Curr, e) {
                            break 'body;
                        }
                    }
                }
            }
            let mut num_2_move = num_2_remove_plus.min(num_2_add_plus);
            if self.curr.num_edges > 0 && cur_charge.num_edges > 0 {
                let mut i = 0;
                while 0 < num_2_move && i < n_wrong {
                    let t = tdhc[i];
                    i += 1;
                    if !(t.n_num_h_inchi > t.n_num_h_norml
                        && t.n_num_p_revrs > t.n_num_p_norml
                        && t.n[F_NUM_R_POS_CHG_H] != 0)
                    {
                        continue;
                    }
                    let mut num_2_remove = i32::from(t.n_num_h_inchi - t.n_num_h_norml);
                    if num_2_remove < i32::from(t.n[F_NUM_R_POS_CHG_H]) {
                        num_2_remove = i32::from(t.n[F_NUM_R_POS_CHG_H]);
                    }
                    let off = t.i[F_NUM_R_POS_CHG_H] as usize;
                    set_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                    remove_forbidden_edge_mask(self.bns, &self.curr, self.mask);
                    let mut jj = 0usize;
                    while 0 < num_2_move && 0 < num_2_remove && jj < t.n[F_NUM_R_POS_CHG_H] as usize
                    {
                        let iat = endp_list.pn_edges[off + jj] as usize;
                        jj += 1;
                        let e = self.plus_edge(iat);
                        let pe = e as usize;
                        if self.bns.edge[pe].flow != 0 {
                            continue;
                        }
                        let v1 = self.bns.edge[pe].neighbor1 as usize;
                        let v2 = (self.bns.edge[pe].neighbor12 as usize) ^ v1;
                        let Some((pe1n, v1n)) = self.flow_neighbor(v1, usize::MAX) else {
                            continue;
                        };
                        // as in C the last edge of v2 is not considered
                        let mut found2 = None;
                        for j in (0..(self.bns.vert[v2].num_adj_edges as i32 - 1)).rev() {
                            let pen = self.vert_edge(v2, j as usize);
                            if self.bns.edge[pen].flow != 0 && self.bns.edge[pen].forbidden == 0 {
                                found2 = Some((pen, (self.bns.edge[pen].neighbor12 as usize) ^ v2));
                                break;
                            }
                        }
                        let Some((pe2n, v2n)) = found2 else { continue };
                        self.bns.edge[pe].flow += 1;
                        self.bns.edge[pe1n].flow -= 1;
                        self.bns.edge[pe2n].flow -= 1;
                        self.bns.vert[v1n].st_edge.flow -= 1;
                        self.bns.vert[v2n].st_edge.flow -= 1;
                        self.bns.tot_st_flow -= 2;
                        self.ret = run_bns_test_once(self.bns, self.bd, self.p_va, &mut self.res);
                        if self.ret == 1
                            && self.path_between(v1n, v2n)
                            && (self.res.delta_charge == 0 || self.res.delta_charge == 1)
                        {
                            self.ret = run_bns_restore_once(self.bns, self.bd);
                            if self.ret > 0 {
                                self.num_run_bns += 1;
                                num_2_remove -= 1;
                                num_2_move -= 1;
                                cur_success += 1;
                            }
                        } else {
                            self.bns.edge[pe].flow -= 1;
                            self.bns.edge[pe1n].flow += 1;
                            self.bns.edge[pe2n].flow += 1;
                            self.bns.vert[v1n].st_edge.flow += 1;
                            self.bns.vert[v2n].st_edge.flow += 1;
                            self.bns.tot_st_flow += 2;
                        }
                        if self.add_edge(L::Curr, e) {
                            break 'body;
                        }
                    }
                    remove_forbidden_edge_mask(self.bns, &self.all_charge, self.mask);
                }
            }
        }
        self.curr.num_edges = 0;
        self.recalc(cur_success)
    }
}

/// `FixFixedHRestoredStructure`: processes the extra or missing Fixed-H on
/// non-tautomeric atoms; `at2` should be the most recently restored atoms.
/// Returns a negative value on an error, otherwise 1 when something has been
/// fixed and the original and the restored InChI still differ.
#[allow(clippy::too_many_arguments)]
pub fn fix_fixed_h_restored_structure(
    ip: &InputParms,
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    at3: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    p_inchi: [Option<&InChI>; 2],
    forbidden_edge_mask: i32,
    forbidden_stereo_edge_mask: i32,
) -> i32 {
    let mut c = FixFhCtx::new(
        ip,
        bns,
        bd,
        ps,
        at2,
        at3,
        p_va,
        tcg,
        p_inchi,
        forbidden_edge_mask,
    );
    let _ = c.run(forbidden_stereo_edge_mask);
    c.result()
}

impl<'a> FixFhCtx<'a> {
    /// the prologue and all the cases in the order of the C function; true
    /// when the function exited early
    fn run(&mut self, forbidden_stereo_edge_mask: i32) -> bool {
        if self.prologue(forbidden_stereo_edge_mask) {
            return true;
        }
        let cases: [fn(&mut Self) -> bool; 33] = [
            Self::case_01,
            Self::case_02,
            Self::case_03,
            Self::case_03a,
            Self::case_04,
            Self::case_05,
            Self::case_06c,
            Self::case_06d,
            Self::case_06,
            Self::case_06a,
            Self::case_06b,
            Self::case_06e,
            Self::case_07,
            Self::case_07a,
            Self::case_08,
            Self::case_09,
            Self::case_10,
            Self::case_11,
            Self::case_12,
            Self::case_13,
            Self::case_14,
            Self::case_15,
            Self::case_16,
            Self::case_17,
            Self::case_18,
            Self::case_19,
            Self::case_20,
            Self::case_21,
            Self::case_21a,
            Self::case_22,
            Self::case_23,
            Self::case_24,
            Self::case_25,
        ];
        for case in cases {
            if case(self) {
                return true;
            }
        }
        false
    }
}
