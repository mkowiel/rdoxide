//! Port of `ichicano.c` (canonicalization driver): the canonical-length,
//! isotopic / tautomeric linear connection table builders and the
//! `CANON_STAT` structure.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichican2::{AtomSizes, Bcn};
use super::ichierr::{CT_LEN_MISMATCH, CT_OVERFLOW};
use super::ichimap::{NeighList, comp_rank};
use super::ichisort::{inchi_qsort, insertions_sort};
use super::ichister::{AtStereoCarb, AtStereoDble, SpAtom};
use super::ichitaut::{T_NUM_ISOTOPIC, T_NUM_NO_ISOTOPIC, TGroupInfo};
use super::inpdef::{AT_FLAG_ISO_H_POINT, MAX_NUM_STEREO_BONDS};

/// `T_GROUP_HDR_LEN`
pub const T_GROUP_HDR_LEN: usize = 1 + T_NUM_NO_ISOTOPIC;
/// `TGSO_SYMM_RANK`: symmetry ranks (no isotopes)
pub const TGSO_SYMM_RANK: usize = 1;
/// `TGSO_SYMM_IORDER`: isotopic symmetry rank sorting order
pub const TGSO_SYMM_IORDER: usize = 2;
/// `TGSO_SYMM_IRANK`: isotopic symmetry ranks
pub const TGSO_SYMM_IRANK: usize = 3;

/// `AT_ISOTOPIC`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AtIsotopic {
    pub at_num: u16,
    pub num_1h: i16,
    pub num_d: i16,
    pub num_t: i16,
    pub iso_atw_diff: i16,
}

/// `AT_ISO_TGROUP`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AtIsoTgroup {
    pub tgroup_num: u16,
    pub num: [u16; T_NUM_ISOTOPIC],
}

/// `CANON_STAT`
#[derive(Debug, Clone, Default)]
pub struct CanonStat {
    pub l_num_break_ties: i64,
    pub l_num_neigh_list_iter: i64,
    pub l_num_tot_ct: i64,
    pub l_num_decreased_ct: i64,
    pub l_num_rejected_ct: i64,
    pub l_num_equal_ct: i64,
    pub l_total_time: i64,
    pub b_first_ct: i32,
    pub b_keep_symm_rank: i32,
    pub b_stereo_is_better: i32,
    pub n_canon_flags: i32,
    pub linear_ct: Vec<u16>,
    pub linear_ct_isotopic: Vec<AtIsotopic>,
    pub linear_ct_isotopic_tautomer: Vec<AtIsoTgroup>,
    pub linear_ct_stereo_dble: Vec<AtStereoDble>,
    pub linear_ct_stereo_carb: Vec<AtStereoCarb>,
    pub linear_ct_stereo_dble_inv: Vec<AtStereoDble>,
    pub linear_ct_stereo_carb_inv: Vec<AtStereoCarb>,
    pub linear_ct_isotopic_stereo_dble: Vec<AtStereoDble>,
    pub linear_ct_isotopic_stereo_carb: Vec<AtStereoCarb>,
    pub linear_ct_isotopic_stereo_dble_inv: Vec<AtStereoDble>,
    pub linear_ct_isotopic_stereo_carb_inv: Vec<AtStereoCarb>,
    pub linear_ct_tautomer: Vec<u16>,
    pub linear_ct2: Vec<u16>,
    pub n_len_linear_ct_stereo_dble: i32,
    pub n_len_linear_ct_stereo_dble_inv: i32,
    pub n_max_len_linear_ct_stereo_dble: i32,
    pub b_cmp_stereo: i32,
    pub n_len_linear_ct_stereo_carb: i32,
    pub n_len_linear_ct_stereo_carb_inv: i32,
    pub n_max_len_linear_ct_stereo_carb: i32,
    pub n_len_linear_ct_isotopic: i32,
    pub n_max_len_linear_ct_isotopic: i32,
    pub n_len_linear_ct_isotopic_tautomer: i32,
    pub n_max_len_linear_ct_isotopic_tautomer: i32,
    pub n_len_linear_ct: i32,
    pub n_len_linear_ct2: i32,
    pub n_len_linear_ct_at_only: i32,
    pub n_len_linear_ct_at_only2: i32,
    pub n_max_len_linear_ct: i32,
    pub n_len_linear_ct_tautomer: i32,
    pub n_max_len_linear_ct_tautomer: i32,
    pub b_cmp_isotopic_stereo: i32,
    pub n_len_linear_ct_isotopic_stereo_dble: i32,
    pub n_len_linear_ct_isotopic_stereo_dble_inv: i32,
    pub n_max_len_linear_ct_isotopic_stereo_dble: i32,
    pub n_len_linear_ct_isotopic_stereo_carb: i32,
    pub n_len_linear_ct_isotopic_stereo_carb_inv: i32,
    pub n_max_len_linear_ct_isotopic_stereo_carb: i32,
    pub b_rank_used_for_stereo: Vec<i8>,
    pub b_atom_used_for_stereo: Vec<i8>,
    pub n_prev_atom_number: Vec<u16>,
    pub n_canon_ord: Vec<u16>,
    pub n_symm_rank: Vec<u16>,
    pub n_canon_ord_taut: Vec<u16>,
    pub n_symm_rank_taut: Vec<u16>,
    pub n_canon_ord_stereo: Vec<u16>,
    pub n_canon_ord_stereo_inv: Vec<u16>,
    pub n_canon_ord_stereo_taut: Vec<u16>,
    pub n_symm_rank_isotopic: Vec<u16>,
    pub n_canon_ord_isotopic: Vec<u16>,
    pub n_symm_rank_isotopic_taut: Vec<u16>,
    pub n_canon_ord_isotopic_taut: Vec<u16>,
    pub n_canon_ord_isotopic_stereo: Vec<u16>,
    pub n_canon_ord_isotopic_stereo_inv: Vec<u16>,
    pub n_canon_ord_isotopic_stereo_taut: Vec<u16>,
    pub n_len_canon_ord: i32,
    pub n_len_canon_ord_taut: i32,
    pub n_len_canon_ord_isotopic: i32,
    pub n_len_canon_ord_isotopic_taut: i32,
    pub n_len_canon_ord_stereo: i32,
    pub n_len_canon_ord_stereo_taut: i32,
    pub n_len_canon_ord_isotopic_stereo: i32,
    pub n_len_canon_ord_isotopic_stereo_taut: i32,
    pub b_has_isotopic_in_tautomer_groups: i32,
    pub t_group_info: Option<TGroupInfo>,
    pub b_ignore_isotopic: i32,
    pub b_double_bond_square: i32,
    pub n_mode: i64,
    pub neigh_list: Option<NeighList>,
    pub p_bcn: Option<Box<Bcn>>,
    pub n_num_h: Vec<i8>,
    pub n_num_h_fixed: Vec<i8>,
    pub n_exchg_iso_h: Vec<i8>,
    /// which of the four buffers the C pointers `LinearCTStereoDble`,
    /// `LinearCTStereoCarb`, `LinearCTStereoDbleInv` and
    /// `LinearCTStereoCarbInv` currently designate
    pub stereo_slots: StereoSlots,
}

/// One of the four stereo CT buffers of `CANON_STAT`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Slot {
    /// `LinearCTStereo...`
    Base,
    /// `LinearCTStereo...Inv`
    Inv,
    /// `LinearCTIsotopicStereo...`
    Iso,
    /// `LinearCTIsotopicStereo...Inv`
    IsoInv,
}

/// The current targets of the stereo CT pointers of `CANON_STAT`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StereoSlots {
    pub dble: Slot,
    pub carb: Slot,
    pub dble_inv: Slot,
    pub carb_inv: Slot,
}

impl Default for StereoSlots {
    fn default() -> Self {
        StereoSlots {
            dble: Slot::Base,
            carb: Slot::Base,
            dble_inv: Slot::Inv,
            carb_inv: Slot::Inv,
        }
    }
}

impl CanonStat {
    /// the stereo double bond buffer `slot`, read only
    pub fn dble_ref(&self, slot: Slot) -> &Vec<AtStereoDble> {
        match slot {
            Slot::Base => &self.linear_ct_stereo_dble,
            Slot::Inv => &self.linear_ct_stereo_dble_inv,
            Slot::Iso => &self.linear_ct_isotopic_stereo_dble,
            Slot::IsoInv => &self.linear_ct_isotopic_stereo_dble_inv,
        }
    }
    /// the stereo center buffer `slot`, read only
    pub fn carb_ref(&self, slot: Slot) -> &Vec<AtStereoCarb> {
        match slot {
            Slot::Base => &self.linear_ct_stereo_carb,
            Slot::Inv => &self.linear_ct_stereo_carb_inv,
            Slot::Iso => &self.linear_ct_isotopic_stereo_carb,
            Slot::IsoInv => &self.linear_ct_isotopic_stereo_carb_inv,
        }
    }
    /// the stereo double bond buffer `slot`
    pub fn dble_buf(&mut self, slot: Slot) -> &mut Vec<AtStereoDble> {
        match slot {
            Slot::Base => &mut self.linear_ct_stereo_dble,
            Slot::Inv => &mut self.linear_ct_stereo_dble_inv,
            Slot::Iso => &mut self.linear_ct_isotopic_stereo_dble,
            Slot::IsoInv => &mut self.linear_ct_isotopic_stereo_dble_inv,
        }
    }
    /// the stereo center buffer `slot`
    pub fn carb_buf(&mut self, slot: Slot) -> &mut Vec<AtStereoCarb> {
        match slot {
            Slot::Base => &mut self.linear_ct_stereo_carb,
            Slot::Inv => &mut self.linear_ct_stereo_carb_inv,
            Slot::Iso => &mut self.linear_ct_isotopic_stereo_carb,
            Slot::IsoInv => &mut self.linear_ct_isotopic_stereo_carb_inv,
        }
    }
}

/// `CHECK_OVERFLOW`
fn check_overflow(len: i32, max_len: i32) -> bool {
    len >= max_len
}

/// `GetCanonLengths`: the lengths of the canonical representation in
/// `AT_NUMB` units; enlarges the lengths in `s`
pub fn get_canon_lengths(
    num_at: usize,
    at: &[SpAtom],
    s: &mut AtomSizes,
    t_group_info: Option<&TGroupInfo>,
) -> i32 {
    let mut n_num_bonds = 0i32;
    let mut n_num_t_bonds = 0i32;
    let mut n_num_dbl_bonds_stereo = 0i32;
    let mut n_num_asym_carb_stereo = 0i32;
    let mut n_num_isotopic = 0i32;
    let t_group = if s.n_len_linear_ct_tautomer != 0 {
        t_group_info.and_then(|t| t.t_group.as_deref())
    } else {
        None
    };
    for a in &at[..num_at] {
        n_num_bonds += i32::from(a.valence);
        if a.iso_sort_key != 0 {
            n_num_isotopic += 1;
        }
        if a.parity > 0 {
            let mut j = 0;
            let mut n_stereo_bonds_to_atom = 0;
            while j < MAX_NUM_STEREO_BONDS && a.stereo_bond_neighbor[j] != 0 {
                let k = a.stereo_bond_neighbor[j] as usize;
                n_stereo_bonds_to_atom += i32::from(at[k - 1].parity > 0);
                j += 1;
            }
            n_num_dbl_bonds_stereo += n_stereo_bonds_to_atom;
            n_num_asym_carb_stereo += i32::from(j == 0);
        }
    }
    n_num_dbl_bonds_stereo /= 2;
    n_num_bonds /= 2;
    s.n_len_bonds = s.n_len_bonds.max(n_num_bonds);
    let mut n_num_ct = n_num_bonds;
    n_num_ct += num_at as i32;
    s.n_len_ct_at_only = s.n_len_ct_at_only.max(n_num_ct);
    if let (Some(tg), Some(tgi)) = (t_group, t_group_info) {
        for g in &tg[..tgi.num_t_groups] {
            n_num_t_bonds += i32::from(g.n_num_endpoints);
        }
        n_num_ct += n_num_t_bonds;
        n_num_ct += tgi.num_t_groups as i32;
    }
    n_num_ct = n_num_ct.max(1);
    s.n_len_ct = s.n_len_ct.max(n_num_ct);
    s.n_len_isotopic = s.n_len_isotopic.max(n_num_isotopic);
    s.n_len_linear_ct_stereo_dble = s.n_len_linear_ct_stereo_dble.max(n_num_dbl_bonds_stereo);
    s.n_len_linear_ct_stereo_carb = s.n_len_linear_ct_stereo_carb.max(n_num_asym_carb_stereo);
    if let Some(tgi) = t_group_info {
        s.n_len_isotopic_endpoints = s.n_len_isotopic_endpoints.max(tgi.n_num_isotopic_endpoints);
    }
    0
}

/// `FillIsotopicAtLinearCT`: returns the length of the isotopic CT, 0 when
/// there is no room, or an error code
pub fn fill_isotopic_at_linear_ct(
    num_atoms: usize,
    at: &[SpAtom],
    n_atom_number: &[u16],
    linear_ct_isotopic: &mut [AtIsotopic],
    n_max_len_linear_ct_isotopic: i32,
    pn_len_linear_ct_isotopic: &mut i32,
) -> i32 {
    let mut len = 0i32;
    if n_max_len_linear_ct_isotopic > 0 && !linear_ct_isotopic.is_empty() {
        for x in linear_ct_isotopic[..n_max_len_linear_ct_isotopic as usize].iter_mut() {
            *x = AtIsotopic::default();
        }
    } else {
        return 0;
    }
    for rank in 1..=num_atoms {
        let i = n_atom_number[rank - 1] as usize;
        let a = &at[i];
        let k = a.endpoint == 0
            && a.c_flags & AT_FLAG_ISO_H_POINT == 0
            && (a.num_iso_h[0] != 0 || a.num_iso_h[1] != 0 || a.num_iso_h[2] != 0);
        if a.iso_atw_diff != 0 || k {
            if check_overflow(len, n_max_len_linear_ct_isotopic) {
                return CT_OVERFLOW;
            }
            let e = &mut linear_ct_isotopic[len as usize];
            e.at_num = rank as u16;
            e.iso_atw_diff = i16::from(a.iso_atw_diff);
            e.num_1h = if k { i16::from(a.num_iso_h[0]) } else { 0 };
            e.num_d = if k { i16::from(a.num_iso_h[1]) } else { 0 };
            e.num_t = if k { i16::from(a.num_iso_h[2]) } else { 0 };
            len += 1;
        }
    }
    if *pn_len_linear_ct_isotopic != 0 {
        if *pn_len_linear_ct_isotopic != len {
            return CT_LEN_MISMATCH;
        }
    } else {
        *pn_len_linear_ct_isotopic = len;
    }
    len
}

/// `FillTautLinearCT2`: fills the tautomeric (and isotopic tautomeric)
/// linear CT; sorts the t-group endpoints by canonical rank and stores the
/// t-group orders in `tGroupNumber`. Returns the length or an error code
/// (negative length on a program error).
#[allow(clippy::too_many_arguments)]
pub fn fill_taut_linear_ct2(
    num_atoms: usize,
    num_at_tg: usize,
    b_iso_taut: bool,
    n_rank: &[u16],
    n_atom_number: &[u16],
    n_symm_rank: &[u16],
    n_atom_number_iso: &[u16],
    n_symm_rank_iso: &[u16],
    linear_ct_tautomer: &mut [u16],
    n_max_len_linear_ct_tautomer: i32,
    pn_len_linear_ct_tautomer: &mut i32,
    linear_ct_isotopic_tautomer: &mut [AtIsoTgroup],
    n_max_len_linear_ct_isotopic_tautomer: i32,
    pn_len_linear_ct_isotopic_tautomer: &mut i32,
    t_group_info: Option<&mut TGroupInfo>,
) -> i32 {
    let max_num_num = T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC;
    let max_num_iso = T_NUM_ISOTOPIC + T_NUM_NO_ISOTOPIC;
    let mut len = 0i32;
    let mut len_iso = 0i32;
    let mut max_len = 0i32;
    let Some(tgi) = t_group_info else {
        return len;
    };
    if !(num_at_tg > num_atoms && tgi.num_t_groups != 0) {
        return len;
    }
    let num_t_groups = tgi.num_t_groups;
    let n_offset = num_atoms as u16;
    for (j, i) in (num_atoms..num_at_tg).enumerate() {
        tgi.t_group_number[j] = n_atom_number[i].wrapping_sub(n_offset);
        tgi.t_group_number[TGSO_SYMM_RANK * num_t_groups + j] =
            n_symm_rank[i].wrapping_sub(n_offset);
        if b_iso_taut {
            tgi.t_group_number[TGSO_SYMM_IORDER * num_t_groups + j] =
                n_atom_number_iso[i].wrapping_sub(n_offset);
            tgi.t_group_number[TGSO_SYMM_IRANK * num_t_groups + j] =
                n_symm_rank_iso[i].wrapping_sub(n_offset);
        }
    }
    // sort the endpoints within each tautomeric group by the canonical ranks
    for i in 0..num_t_groups {
        let (first, n) = {
            let g = &tgi.t_group.as_ref().unwrap()[i];
            (
                g.n_first_endpoint_at_no_pos as usize,
                g.n_num_endpoints as usize,
            )
        };
        inchi_qsort(&mut tgi.n_endpoint_atom_number[first..first + n], |a, b| {
            comp_rank(n_rank, *a, *b)
        });
    }
    if n_max_len_linear_ct_tautomer != 0 {
        max_len = (T_GROUP_HDR_LEN * num_t_groups + tgi.n_num_endpoints as usize + 1) as i32;
        if max_len > n_max_len_linear_ct_tautomer {
            return CT_OVERFLOW;
        }
    }
    let t_group = tgi.t_group.as_ref().unwrap();
    for i in 0..num_t_groups {
        let g = tgi.t_group_number[i] as usize;
        let tg = &t_group[g];
        if check_overflow(
            len + T_GROUP_HDR_LEN as i32 + i32::from(tg.n_num_endpoints),
            max_len,
        ) {
            return CT_OVERFLOW;
        }
        linear_ct_tautomer[len as usize] = tg.n_num_endpoints;
        len += 1;
        for j in 0..T_NUM_NO_ISOTOPIC {
            linear_ct_tautomer[len as usize] = tg.num[j];
            len += 1;
        }
        let offset = tg.n_first_endpoint_at_no_pos as usize;
        for j in 0..tg.n_num_endpoints as usize {
            linear_ct_tautomer[len as usize] =
                n_rank[tgi.n_endpoint_atom_number[offset + j] as usize];
            len += 1;
        }
    }
    if n_max_len_linear_ct_tautomer != 0 {
        linear_ct_tautomer[len as usize] = 0;
        len += 1;
        if len != max_len {
            len = -len;
        } else if *pn_len_linear_ct_tautomer != 0 && *pn_len_linear_ct_tautomer != len {
            return CT_LEN_MISMATCH;
        } else {
            *pn_len_linear_ct_tautomer = len;
        }
    } else {
        *pn_len_linear_ct_tautomer = 0;
    }
    // isotopic tautomeric mobile groups
    if n_max_len_linear_ct_isotopic_tautomer != 0 && tgi.n_num_isotopic_endpoints == 0 {
        for i in 0..num_t_groups {
            let g = tgi.t_group_number[TGSO_SYMM_IORDER * num_t_groups + i] as usize;
            let tg = &t_group[g];
            if tg.i_weight == 0 {
                continue;
            }
            if check_overflow(len_iso, n_max_len_linear_ct_isotopic_tautomer) {
                return CT_OVERFLOW;
            }
            let mut j = T_NUM_NO_ISOTOPIC;
            while j < max_num_num && j < max_num_iso {
                linear_ct_isotopic_tautomer[len_iso as usize].num[j - T_NUM_NO_ISOTOPIC] =
                    tg.num[j];
                j += 1;
            }
            linear_ct_isotopic_tautomer[len_iso as usize].tgroup_num = (i + 1) as u16;
            len_iso += 1;
        }
    }
    if n_max_len_linear_ct_isotopic_tautomer != 0 {
        if *pn_len_linear_ct_isotopic_tautomer != 0
            && *pn_len_linear_ct_isotopic_tautomer != len_iso
        {
            return CT_LEN_MISMATCH;
        }
        *pn_len_linear_ct_isotopic_tautomer = len_iso;
    } else {
        *pn_len_linear_ct_isotopic_tautomer = 0;
    }
    len
}

/// `UpdateFullLinearCT`: updates the linear connection table `pCS->LinearCT`
/// with the final ranks; returns 1 when the new CT is not minimal, 0 when it
/// is identical, -1 when it is smaller (or the first one), or an error code
pub fn update_full_linear_ct(
    num_atoms: usize,
    num_at_tg: usize,
    at: &[SpAtom],
    n_rank: &[u16],
    n_atom_number: &[u16],
    cs: &mut CanonStat,
    b_first_time: bool,
) -> i32 {
    let mut n_ct_len = 0i32;
    let mut b_compare = !b_first_time;
    macro_rules! compare_with_ct {
        ($value:expr) => {{
            let value: u16 = $value;
            if b_compare {
                if value > cs.linear_ct[n_ct_len as usize] {
                    return 1;
                }
                b_compare = value == cs.linear_ct[n_ct_len as usize];
            }
            cs.linear_ct[n_ct_len as usize] = value;
            n_ct_len += 1;
        }};
    }
    let have_tg = num_at_tg > num_atoms;
    let mut num_neigh = 0usize;
    for rank in 1..=num_atoms {
        let i = n_atom_number[rank - 1] as usize;
        if check_overflow(n_ct_len, cs.n_max_len_linear_ct) {
            return CT_OVERFLOW;
        }
        compare_with_ct!(rank as u16);
        num_neigh = at[i].valence as usize;
        let mut n_neighbor_number: Vec<u16> = (0..num_neigh as u16).collect();
        let nb = &at[i].neighbor;
        insertions_sort(&mut n_neighbor_number, |a, b| {
            i32::from(n_rank[nb[*a as usize] as usize])
                - i32::from(n_rank[nb[*b as usize] as usize])
        });
        for k in 0..num_neigh {
            let r_neigh = n_rank[nb[n_neighbor_number[k] as usize] as usize];
            if (r_neigh as usize) < rank {
                if check_overflow(n_ct_len, cs.n_max_len_linear_ct) {
                    return CT_OVERFLOW;
                }
                compare_with_ct!(r_neigh);
            }
        }
    }
    let n_ct_len_at_only = n_ct_len;
    for rank in num_atoms + 1..=num_at_tg {
        let j = n_atom_number[rank - 1] as usize;
        let i = j - num_atoms;
        if check_overflow(n_ct_len, cs.n_max_len_linear_ct) {
            return CT_OVERFLOW;
        }
        compare_with_ct!(rank as u16);
        let tgi = cs.t_group_info.as_mut().unwrap();
        let (first, n_ep) = {
            let g = &tgi.t_group.as_ref().unwrap()[i];
            (
                g.n_first_endpoint_at_no_pos as usize,
                g.n_num_endpoints as usize,
            )
        };
        if have_tg {
            num_neigh = n_ep;
        }
        insertions_sort(
            &mut tgi.n_endpoint_atom_number[first..first + num_neigh],
            |a, b| comp_rank(n_rank, *a, *b),
        );
        for k in 0..num_neigh {
            let r_neigh = n_rank[tgi.n_endpoint_atom_number[first + k] as usize];
            if (r_neigh as usize) < rank {
                if check_overflow(n_ct_len, cs.n_max_len_linear_ct) {
                    return CT_OVERFLOW;
                }
                compare_with_ct!(r_neigh);
            }
        }
    }
    if cs.n_len_linear_ct != 0 {
        if cs.n_len_linear_ct != n_ct_len {
            return CT_LEN_MISMATCH;
        }
    } else {
        cs.n_len_linear_ct = n_ct_len;
    }
    if cs.n_len_linear_ct != 0 {
        if cs.n_len_linear_ct_at_only != n_ct_len_at_only {
            return CT_LEN_MISMATCH;
        }
    } else {
        cs.n_len_linear_ct_at_only = n_ct_len_at_only;
    }
    i32::from(b_compare) - 1
}

// ---------------------------------------------------------------------------
// AllocateCS / Canon_INChI3
// ---------------------------------------------------------------------------

use super::ichican2::fix_canon_equivalence_info;
use super::ichierr::{
    CT_CALC_STEREO_ERR, CT_CANON_ERR, CT_ERR_MAX, CT_ERR_MIN, CT_OUT_OF_RAM, CT_STEREOCOUNT_ERR,
    CT_TAUCOUNT_ERR,
};
use super::ichimap::{CurTree, RemoveStereoCs, compare_lin_ct_stereo};
use super::ichimap4::{StereoMapCs, map_stereo_bonds4};
use super::ichister::{
    CanonStatStereo, fill_out_stereo_parities, invert_stereo,
    switch_atom_stereo_and_isotopic_stereo,
};
use super::ichitaut::sort_tautomer_groups_and_endpoints;
use super::runichi3::{TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE};

pub const CMODE_CT: i64 = 0x01;
pub const CMODE_ISO: i64 = 0x02;
pub const CMODE_STEREO: i64 = 0x08;
pub const CMODE_ISO_STEREO: i64 = 0x10;
pub const CMODE_TAUT: i64 = 0x20;
pub const CMODE_NOEQ_STEREO: i64 = 0x40;
pub const CANON_MODE_TAUT: i64 = CMODE_CT | CMODE_TAUT;
pub const CANON_MODE_MASK: i64 = 0xFF;
/// `REQ_MODE_DIFF_UU_STEREO`
pub const REQ_MODE_DIFF_UU_STEREO: i64 = 0x4000;

fn returned_error(v: i32) -> bool {
    (CT_ERR_MIN..=CT_ERR_MAX).contains(&v)
}

/// `AllocateCS`: allocates the arrays of `pCS` that the mode requires
/// (`has_bcn`: `pBCN` is not NULL)
#[allow(clippy::too_many_arguments)]
pub fn allocate_cs(
    cs: &mut CanonStat,
    num_at: usize,
    num_at_tg: usize,
    n_len_ct: i32,
    n_len_ct_at_only: i32,
    n_len_linear_ct_stereo_dble: i32,
    n_len_linear_ct_isotopic_stereo_dble: i32,
    n_len_linear_ct_stereo_carb: i32,
    n_len_linear_ct_isotopic_stereo_carb: i32,
    n_len_linear_ct_tautomer: i32,
    n_len_linear_ct_isotopic_tautomer: i32,
    n_len_isotopic: i32,
    n_mode: i64,
    has_bcn: bool,
) -> i32 {
    let num_t_groups = num_at_tg as i64 - num_at as i64;
    cs.n_mode = n_mode;
    if n_mode & CMODE_CT != 0 && n_len_ct > 0 {
        cs.linear_ct = vec![0; n_len_ct as usize];
        cs.n_max_len_linear_ct = n_len_ct;
        cs.n_len_linear_ct = n_len_ct;
        cs.n_len_linear_ct_at_only = n_len_ct_at_only;
        cs.n_canon_ord = vec![0; num_at_tg];
        cs.n_symm_rank = vec![0; num_at_tg];
        if has_bcn {
            cs.n_num_h = vec![0; num_at];
            cs.n_num_h_fixed = vec![0; num_at];
            cs.n_exchg_iso_h = vec![0; num_at];
        }
    }
    if n_mode & CMODE_ISO != 0 && n_len_isotopic > 0 {
        cs.linear_ct_isotopic = vec![AtIsotopic::default(); n_len_isotopic as usize];
        cs.n_max_len_linear_ct_isotopic = n_len_isotopic;
        cs.n_len_linear_ct_isotopic = n_len_isotopic;
    }
    if n_mode & CMODE_ISO != 0 && n_mode & CANON_MODE_TAUT == CANON_MODE_TAUT {
        if n_len_linear_ct_isotopic_tautomer > 0 {
            cs.linear_ct_isotopic_tautomer =
                vec![AtIsoTgroup::default(); n_len_linear_ct_isotopic_tautomer as usize];
            cs.n_max_len_linear_ct_isotopic_tautomer = n_len_linear_ct_isotopic_tautomer;
            cs.n_len_linear_ct_isotopic_tautomer = n_len_linear_ct_isotopic_tautomer;
        }
        if num_t_groups > 0 {
            cs.n_canon_ord_isotopic_taut = vec![0; num_t_groups as usize];
            cs.n_symm_rank_isotopic_taut = vec![0; num_t_groups as usize];
        }
    }
    if n_mode & CMODE_ISO != 0 {
        cs.n_symm_rank_isotopic = vec![0; num_at_tg];
        cs.n_canon_ord_isotopic = vec![0; num_at_tg];
    }
    if n_mode & CMODE_STEREO != 0 && n_len_linear_ct_stereo_dble > 0 {
        cs.linear_ct_stereo_dble =
            vec![AtStereoDble::default(); n_len_linear_ct_stereo_dble as usize];
        cs.linear_ct_stereo_dble_inv =
            vec![AtStereoDble::default(); n_len_linear_ct_stereo_dble as usize];
        cs.n_len_linear_ct_stereo_dble_inv = n_len_linear_ct_stereo_dble;
        cs.n_max_len_linear_ct_stereo_dble = n_len_linear_ct_stereo_dble;
        cs.n_len_linear_ct_stereo_dble = n_len_linear_ct_stereo_dble;
    }
    if n_mode & CMODE_STEREO != 0 && n_len_linear_ct_stereo_carb > 0 {
        cs.linear_ct_stereo_carb =
            vec![AtStereoCarb::default(); n_len_linear_ct_stereo_carb as usize];
        cs.linear_ct_stereo_carb_inv =
            vec![AtStereoCarb::default(); n_len_linear_ct_stereo_carb as usize];
        cs.n_len_linear_ct_stereo_carb_inv = n_len_linear_ct_stereo_carb;
        cs.n_max_len_linear_ct_stereo_carb = n_len_linear_ct_stereo_carb;
        cs.n_len_linear_ct_stereo_carb = n_len_linear_ct_stereo_carb;
    }
    if n_mode & CMODE_STEREO != 0
        && (n_len_linear_ct_stereo_dble > 0 || n_len_linear_ct_stereo_carb > 0)
    {
        cs.n_canon_ord_stereo = vec![0; num_at_tg];
        cs.n_canon_ord_stereo_inv = vec![0; num_at_tg];
        if n_mode & CMODE_TAUT != 0 && n_len_linear_ct_tautomer > 0 && num_t_groups > 0 {
            cs.n_canon_ord_stereo_taut = vec![0; num_t_groups as usize];
        }
    }
    if n_mode & CMODE_ISO_STEREO != 0 && n_len_linear_ct_isotopic_stereo_dble > 0 {
        cs.linear_ct_isotopic_stereo_dble =
            vec![AtStereoDble::default(); n_len_linear_ct_isotopic_stereo_dble as usize];
        cs.linear_ct_isotopic_stereo_dble_inv =
            vec![AtStereoDble::default(); n_len_linear_ct_isotopic_stereo_dble as usize];
        cs.n_len_linear_ct_isotopic_stereo_dble_inv = n_len_linear_ct_isotopic_stereo_dble;
        cs.n_max_len_linear_ct_isotopic_stereo_dble = n_len_linear_ct_isotopic_stereo_dble;
        cs.n_len_linear_ct_isotopic_stereo_dble = n_len_linear_ct_isotopic_stereo_dble;
    }
    if n_mode & CMODE_ISO_STEREO != 0 && n_len_linear_ct_isotopic_stereo_carb > 0 {
        cs.linear_ct_isotopic_stereo_carb =
            vec![AtStereoCarb::default(); n_len_linear_ct_isotopic_stereo_carb as usize];
        cs.linear_ct_isotopic_stereo_carb_inv =
            vec![AtStereoCarb::default(); n_len_linear_ct_isotopic_stereo_carb as usize];
        cs.n_len_linear_ct_isotopic_stereo_carb_inv = n_len_linear_ct_isotopic_stereo_carb;
        cs.n_max_len_linear_ct_isotopic_stereo_carb = n_len_linear_ct_isotopic_stereo_carb;
        cs.n_len_linear_ct_isotopic_stereo_carb = n_len_linear_ct_isotopic_stereo_carb;
    }
    if n_mode & CMODE_ISO_STEREO != 0
        && (n_len_linear_ct_isotopic_stereo_dble > 0 || n_len_linear_ct_isotopic_stereo_carb > 0)
    {
        cs.n_canon_ord_isotopic_stereo = vec![0; num_at_tg];
        cs.n_canon_ord_isotopic_stereo_inv = vec![0; num_at_tg];
        if n_mode & CMODE_TAUT != 0 && n_len_linear_ct_tautomer > 0 && num_t_groups > 0 {
            cs.n_canon_ord_isotopic_stereo_taut = vec![0; num_t_groups as usize];
        }
    }
    if (n_mode & CMODE_STEREO != 0
        && (n_len_linear_ct_stereo_dble > 0 || n_len_linear_ct_stereo_carb > 0))
        || (n_mode & CMODE_ISO_STEREO != 0
            && (n_len_linear_ct_isotopic_stereo_dble > 0
                || n_len_linear_ct_isotopic_stereo_carb > 0))
    {
        cs.b_rank_used_for_stereo = vec![0; num_at];
        cs.b_atom_used_for_stereo = vec![0; num_at];
    }
    if n_mode & CMODE_CT != 0 && n_mode & CMODE_TAUT != 0 && n_len_linear_ct_tautomer > 0 {
        cs.linear_ct_tautomer = vec![0; n_len_linear_ct_tautomer as usize];
        cs.n_max_len_linear_ct_tautomer = n_len_linear_ct_tautomer;
        cs.n_len_linear_ct_tautomer = n_len_linear_ct_tautomer;
        if num_t_groups > 0 {
            cs.n_canon_ord_taut = vec![0; num_t_groups as usize];
            cs.n_symm_rank_taut = vec![0; num_t_groups as usize];
        }
    }
    if n_mode & CMODE_CT != 0 {
        cs.linear_ct2 = vec![0; n_len_ct.max(0) as usize];
    }
    cs.n_prev_atom_number = vec![0; num_at_tg];
    cs.n_len_canon_ord = 0;
    cs.n_len_canon_ord_isotopic = 0;
    cs.n_len_canon_ord_isotopic_taut = 0;
    cs.n_len_canon_ord_stereo = 0;
    cs.n_len_canon_ord_stereo_taut = 0;
    cs.n_len_canon_ord_isotopic_stereo = 0;
    cs.n_len_canon_ord_isotopic_stereo_taut = 0;
    cs.n_len_canon_ord_taut = 0;
    0
}

/// `DeAllocateCS`
pub fn de_allocate_cs(cs: &mut CanonStat) -> i32 {
    cs.linear_ct = Vec::new();
    cs.n_canon_ord = Vec::new();
    cs.n_symm_rank = Vec::new();
    cs.n_num_h = Vec::new();
    cs.n_num_h_fixed = Vec::new();
    cs.n_exchg_iso_h = Vec::new();
    cs.linear_ct_isotopic = Vec::new();
    cs.n_symm_rank_isotopic = Vec::new();
    cs.n_canon_ord_isotopic = Vec::new();
    cs.linear_ct_isotopic_tautomer = Vec::new();
    cs.n_canon_ord_isotopic_taut = Vec::new();
    cs.n_symm_rank_isotopic_taut = Vec::new();
    cs.linear_ct_stereo_dble = Vec::new();
    cs.linear_ct_stereo_carb = Vec::new();
    cs.linear_ct_stereo_dble_inv = Vec::new();
    cs.linear_ct_stereo_carb_inv = Vec::new();
    cs.n_canon_ord_stereo = Vec::new();
    cs.n_canon_ord_stereo_inv = Vec::new();
    cs.n_canon_ord_stereo_taut = Vec::new();
    cs.linear_ct_isotopic_stereo_dble = Vec::new();
    cs.linear_ct_isotopic_stereo_carb = Vec::new();
    cs.linear_ct_isotopic_stereo_dble_inv = Vec::new();
    cs.linear_ct_isotopic_stereo_carb_inv = Vec::new();
    cs.b_rank_used_for_stereo = Vec::new();
    cs.b_atom_used_for_stereo = Vec::new();
    cs.n_canon_ord_isotopic_stereo = Vec::new();
    cs.n_canon_ord_isotopic_stereo_inv = Vec::new();
    cs.n_canon_ord_isotopic_stereo_taut = Vec::new();
    cs.linear_ct_tautomer = Vec::new();
    cs.n_canon_ord_taut = Vec::new();
    cs.n_symm_rank_taut = Vec::new();
    cs.linear_ct2 = Vec::new();
    cs.n_prev_atom_number = Vec::new();
    cs.neigh_list = None;
    cs.n_max_len_linear_ct_stereo_dble = 0;
    cs.n_len_linear_ct_stereo_dble = 0;
    cs.n_max_len_linear_ct_stereo_carb = 0;
    cs.n_len_linear_ct_stereo_carb = 0;
    cs.n_max_len_linear_ct_isotopic_stereo_dble = 0;
    cs.n_len_linear_ct_isotopic_stereo_dble = 0;
    cs.n_max_len_linear_ct_isotopic_stereo_carb = 0;
    cs.n_len_linear_ct_isotopic_stereo_carb = 0;
    cs.n_max_len_linear_ct_tautomer = 0;
    cs.n_len_linear_ct_tautomer = 0;
    cs.n_max_len_linear_ct_isotopic = 0;
    cs.n_len_linear_ct_isotopic = 0;
    cs.n_max_len_linear_ct_isotopic_tautomer = 0;
    cs.n_len_linear_ct_isotopic_tautomer = 0;
    cs.n_len_canon_ord = 0;
    cs.n_len_canon_ord_isotopic = 0;
    cs.n_len_canon_ord_isotopic_taut = 0;
    cs.n_len_canon_ord_stereo = 0;
    cs.n_len_canon_ord_stereo_taut = 0;
    cs.n_len_canon_ord_isotopic_stereo = 0;
    cs.n_len_canon_ord_isotopic_stereo_taut = 0;
    cs.n_len_canon_ord_taut = 0;
    0
}

/// runs `f` on the stereo CT buffers that the `CANON_STAT` pointers
/// `LinearCTStereoDble` / `LinearCTStereoCarb` designate
fn with_stereo_view<R>(cs: &mut CanonStat, f: impl FnOnce(&mut CanonStatStereo) -> R) -> R {
    let (sd, sc) = (cs.stereo_slots.dble, cs.stereo_slots.carb);
    let carb = std::mem::take(cs.carb_buf(sc));
    let dble = std::mem::take(cs.dble_buf(sd));
    let mut v = CanonStatStereo {
        carb: (!carb.is_empty()).then_some(carb),
        n_len_carb: cs.n_len_linear_ct_stereo_carb,
        n_max_len_carb: cs.n_max_len_linear_ct_stereo_carb,
        dble: (!dble.is_empty()).then_some(dble),
        n_len_dble: cs.n_len_linear_ct_stereo_dble,
        n_max_len_dble: cs.n_max_len_linear_ct_stereo_dble,
    };
    let r = f(&mut v);
    *cs.carb_buf(sc) = v.carb.unwrap_or_default();
    *cs.dble_buf(sd) = v.dble.unwrap_or_default();
    cs.n_len_linear_ct_stereo_carb = v.n_len_carb;
    cs.n_len_linear_ct_stereo_dble = v.n_len_dble;
    r
}

/// the snapshot `*pCS2 = *pCS` of the stereo lengths
#[derive(Debug, Clone, Copy, Default)]
struct Cs2 {
    n_len_dble: i32,
    n_len_carb: i32,
    n_len_dble_inv: i32,
    n_len_carb_inv: i32,
    n_len_iso_dble: i32,
    n_len_iso_carb: i32,
    n_len_iso_dble_inv: i32,
    n_len_iso_carb_inv: i32,
    n_max_len_dble: i32,
    n_max_len_carb: i32,
    n_max_len_iso_dble: i32,
    n_max_len_iso_carb: i32,
}

impl Cs2 {
    fn snapshot(cs: &CanonStat) -> Cs2 {
        Cs2 {
            n_len_dble: cs.n_len_linear_ct_stereo_dble,
            n_len_carb: cs.n_len_linear_ct_stereo_carb,
            n_len_dble_inv: cs.n_len_linear_ct_stereo_dble_inv,
            n_len_carb_inv: cs.n_len_linear_ct_stereo_carb_inv,
            n_len_iso_dble: cs.n_len_linear_ct_isotopic_stereo_dble,
            n_len_iso_carb: cs.n_len_linear_ct_isotopic_stereo_carb,
            n_len_iso_dble_inv: cs.n_len_linear_ct_isotopic_stereo_dble_inv,
            n_len_iso_carb_inv: cs.n_len_linear_ct_isotopic_stereo_carb_inv,
            n_max_len_dble: cs.n_max_len_linear_ct_stereo_dble,
            n_max_len_carb: cs.n_max_len_linear_ct_stereo_carb,
            n_max_len_iso_dble: cs.n_max_len_linear_ct_isotopic_stereo_dble,
            n_max_len_iso_carb: cs.n_max_len_linear_ct_isotopic_stereo_carb,
        }
    }
}

/// `SetCtToIsotopicStereo`
fn set_ct_to_isotopic_stereo(cs: &mut CanonStat, cs2: &Cs2) {
    cs.stereo_slots.dble = Slot::Iso;
    cs.stereo_slots.carb = Slot::Iso;
    cs.stereo_slots.dble_inv = Slot::IsoInv;
    cs.stereo_slots.carb_inv = Slot::IsoInv;
    cs.n_max_len_linear_ct_stereo_dble = cs2.n_max_len_iso_dble;
    cs.n_max_len_linear_ct_stereo_carb = cs2.n_max_len_iso_carb;
    cs.n_len_linear_ct_stereo_dble = cs2.n_len_iso_dble;
    cs.n_len_linear_ct_stereo_carb = cs2.n_len_iso_carb;
}

/// `SetCtToNonIsotopicStereo`
fn set_ct_to_non_isotopic_stereo(cs: &mut CanonStat, cs2: &Cs2) {
    cs.stereo_slots.dble = Slot::Base;
    cs.stereo_slots.carb = Slot::Base;
    cs.stereo_slots.dble_inv = Slot::Inv;
    cs.stereo_slots.carb_inv = Slot::Inv;
    cs.n_max_len_linear_ct_stereo_dble = cs2.n_max_len_dble;
    cs.n_max_len_linear_ct_stereo_carb = cs2.n_max_len_carb;
    cs.n_len_linear_ct_stereo_dble = cs2.n_len_dble;
    cs.n_len_linear_ct_stereo_carb = cs2.n_len_carb;
    cs.n_len_linear_ct_isotopic_stereo_dble = cs2.n_len_iso_dble;
    cs.n_len_linear_ct_isotopic_stereo_carb = cs2.n_len_iso_carb;
}

/// `map_stereo_bonds4` on the `CANON_STAT` pointers' current buffers
#[allow(clippy::too_many_arguments)]
fn map_stereo_in_cs(
    time_is_over: &mut dyn FnMut() -> bool,
    at: &mut [SpAtom],
    num_atoms: usize,
    num_at_tg: usize,
    num_max: usize,
    n_canon_rank_from: &[u16],
    n_atom_number_canon_from: &[u16],
    n_canon_rank_to: &mut [u16],
    n_symm_rank: &[u16],
    st1: &mut Vec<Vec<u16>>,
    st2: &mut Vec<Vec<u16>>,
    n_temp_rank: &mut [u16],
    n_num_mapped_ranks: i32,
    n_symm_stereo: &mut [u16],
    neigh_list: &mut NeighList,
    cs: &mut CanonStat,
    cur_tree: &mut CurTree,
    v_ab_parity_unknown: i32,
) -> i32 {
    let (sd, sc) = (cs.stereo_slots.dble, cs.stereo_slots.carb);
    let nd = cs.n_len_linear_ct_stereo_dble.max(0) as usize;
    let nc = cs.n_len_linear_ct_stereo_carb.max(0) as usize;
    let mut dble_buf = std::mem::take(cs.dble_buf(sd));
    let mut carb_buf = std::mem::take(cs.carb_buf(sc));
    let mut m = StereoMapCs {
        rm: RemoveStereoCs {
            linear_ct_stereo_carb: carb_buf[..nc.min(carb_buf.len())].to_vec(),
            linear_ct_stereo_dble: dble_buf[..nd.min(dble_buf.len())].to_vec(),
            n_symm_rank: std::mem::take(&mut cs.n_symm_rank),
            l_num_neigh_list_iter: 0,
        },
        counters: Default::default(),
        b_first_ct: cs.b_first_ct,
        b_keep_symm_rank: cs.b_keep_symm_rank,
        b_stereo_is_better: cs.b_stereo_is_better,
        l_num_rejected_ct: cs.l_num_rejected_ct,
        l_num_tot_ct: cs.l_num_tot_ct,
        l_num_equal_ct: cs.l_num_equal_ct,
        l_num_decreased_ct: cs.l_num_decreased_ct,
        b_rank_used_for_stereo: std::mem::take(&mut cs.b_rank_used_for_stereo),
        b_atom_used_for_stereo: std::mem::take(&mut cs.b_atom_used_for_stereo),
        n_prev_atom_number: std::mem::take(&mut cs.n_prev_atom_number),
        n_mode: cs.n_mode as i32,
        time_is_over: time_is_over(),
    };
    m.counters.l_num_break_ties = cs.l_num_break_ties;
    m.counters.l_num_neigh_list_iter = cs.l_num_neigh_list_iter;
    let ret = map_stereo_bonds4(
        at,
        num_atoms,
        num_at_tg,
        num_max,
        false,
        n_canon_rank_from,
        n_atom_number_canon_from,
        n_canon_rank_to,
        n_symm_rank,
        st1,
        st2,
        n_temp_rank,
        n_num_mapped_ranks,
        n_symm_stereo,
        neigh_list,
        &mut m,
        cur_tree,
        0,
        v_ab_parity_unknown,
    );
    let k = m.rm.linear_ct_stereo_carb.len();
    carb_buf[..k].copy_from_slice(&m.rm.linear_ct_stereo_carb);
    let k = m.rm.linear_ct_stereo_dble.len();
    dble_buf[..k].copy_from_slice(&m.rm.linear_ct_stereo_dble);
    *cs.carb_buf(sc) = carb_buf;
    *cs.dble_buf(sd) = dble_buf;
    cs.n_symm_rank = m.rm.n_symm_rank;
    cs.b_first_ct = m.b_first_ct;
    cs.b_keep_symm_rank = m.b_keep_symm_rank;
    cs.b_stereo_is_better = m.b_stereo_is_better;
    cs.l_num_rejected_ct = m.l_num_rejected_ct;
    cs.l_num_tot_ct = m.l_num_tot_ct;
    cs.l_num_equal_ct = m.l_num_equal_ct;
    cs.l_num_decreased_ct = m.l_num_decreased_ct;
    cs.l_num_break_ties = m.counters.l_num_break_ties;
    cs.l_num_neigh_list_iter = m.counters.l_num_neigh_list_iter;
    cs.b_rank_used_for_stereo = m.b_rank_used_for_stereo;
    cs.b_atom_used_for_stereo = m.b_atom_used_for_stereo;
    cs.n_prev_atom_number = m.n_prev_atom_number;
    ret
}

/// `for (n = 2; n < nRankStackLen && pRankStack[n]; n++) pRankStack[n][0] = 0`
fn mark_ranks_stale(st: &mut [Vec<u16>], len: usize) {
    let mut n = 2;
    while n < len && n < st.len() && !st[n].is_empty() {
        st[n][0] = 0;
        n += 1;
    }
}

fn reset_stat_counters(cs: &mut CanonStat) {
    cs.l_num_break_ties = 0;
    cs.l_num_neigh_list_iter = 0;
    cs.l_num_tot_ct = 0;
    cs.l_num_decreased_ct = 0;
    cs.l_num_rejected_ct = 0;
    cs.l_num_equal_ct = 0;
    cs.b_keep_symm_rank = 0;
    cs.b_first_ct = 1;
}

/// the canonical order of `rank` into `ord`: `ord[rank[i]-1] = i`; the
/// number of failures (ranks that are 0 or larger than `num_at_tg`)
fn canon_ord_from_rank(rank: &[u16], num_at_tg: usize, ord: &mut [u16]) -> i32 {
    let mut failed = 0;
    for i in 0..num_at_tg {
        if rank[i] != 0 && rank[i] as usize <= num_at_tg {
            ord[rank[i] as usize - 1] = i as u16;
        } else {
            failed += 1;
        }
    }
    failed
}

/// `Canon_INChI3`: the final canonical numbering of a component: saves the
/// base canonical numbering (`pCS->pBCN`) in `pCS` and optimizes the isotopic
/// and stereo parts. `b_taut_ftcn` selects the `pBCN->ftcn[]` entry. The
/// user-quit callbacks are not ported; `CMODE_NOEQ_STEREO` is not supported
/// (the equivalence tree is always used). Returns `num_atoms` or an error.
#[allow(clippy::too_many_arguments, unused_assignments)]
pub fn canon_inchi3(
    time_is_over: &mut dyn FnMut() -> bool,
    num_atoms: usize,
    num_at_tg: usize,
    at: &mut [SpAtom],
    cs: &mut CanonStat,
    n_mode: i64,
    b_taut_ftcn: usize,
) -> i32 {
    use super::ichister::{AB_PARITY_UNDF, AB_PARITY_UNKN};
    let mut bcn = cs.p_bcn.take().expect("pBCN");
    let n_rank_stack_len = bcn.n_max_len_rank_stack.max(0) as usize;
    let num_max = bcn.num_max;
    // flatten the rank stack: [0] and [1] are not set, the arrays follow
    let mut st1: Vec<Vec<u16>> = vec![Vec::new(), Vec::new()];
    {
        let flat: Vec<Vec<u16>> = bcn
            .p_rank_stack
            .iter()
            .flat_map(|p| [p.rank.clone(), p.at_number.clone()])
            .collect();
        for v in flat.into_iter().skip(2) {
            if v.is_empty() {
                break;
            }
            st1.push(v);
        }
    }
    let mut neigh_list: Option<NeighList> = bcn.ftcn[b_taut_ftcn].neigh_list.take();
    let mut st2: Option<Vec<Vec<u16>>> = None;
    let mut n_ret: i32 = 0;
    let mut b_switched_atom_to_isotopic = false;
    let mut cs2 = Cs2::default();
    {
        let ftcn = &bcn.ftcn[b_taut_ftcn];
        let b_taut = num_at_tg > num_atoms
            && cs
                .t_group_info
                .as_ref()
                .is_some_and(|t| t.num_t_groups != 0 && t.t_group.is_some());
        let b_iso_xchg_h = cs.t_group_info.as_ref().is_some_and(|t| {
            t.n_num_isotopic_endpoints > 1
                && !t.n_isotopic_endpoint_atom_number.is_empty()
                && t.n_isotopic_endpoint_atom_number[0] != 0
                && t.b_taut_flags_done
                    & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                    != 0
        });
        let b_has_isotopic_canon_data = !ftcn.partition_ct_iso.at_number.is_empty()
            && !ftcn.partition_ct_iso.rank.is_empty()
            && !ftcn.n_symm_rank_ct_iso.is_empty();
        let b_canon_isotopic = b_has_isotopic_canon_data
            && n_mode & CMODE_ISO != 0
            && (!cs.linear_ct_isotopic.is_empty()
                || !cs.linear_ct_isotopic_tautomer.is_empty()
                || b_iso_xchg_h);
        let mut b_canon_stereo = n_mode & CMODE_STEREO != 0
            && (!cs.dble_ref(cs.stereo_slots.dble).is_empty()
                || !cs.carb_ref(cs.stereo_slots.carb).is_empty());
        let mut b_canon_iso_stereo = b_has_isotopic_canon_data
            && n_mode & CMODE_ISO_STEREO != 0
            && (!cs.linear_ct_isotopic_stereo_dble.is_empty()
                || !cs.linear_ct_isotopic_stereo_carb.is_empty())
            && b_canon_isotopic;
        let mut b_iso_taut = b_taut && b_canon_isotopic;

        let mut n_num_curr_ranks;
        let mut n_temp_rank: Vec<u16> = Vec::new();
        let mut n_symm_rank: Vec<u16> = Vec::new();
        let mut n_atom_number: Vec<u16> = Vec::new();
        let mut n_rank: Vec<u16> = Vec::new();
        let mut n_canon_rank_stereo: Vec<u16> = Vec::new();
        let mut n_canon_rank_stereo_inv: Vec<u16> = Vec::new();
        let mut n_symm_stereo: Vec<u16> = Vec::new();
        let mut n_canon_rank_isotopic_stereo: Vec<u16> = Vec::new();
        let mut n_canon_rank_isotopic_stereo_inv: Vec<u16> = Vec::new();
        let mut cur_tree = CurTree::default();

        let v_ab_parity_unknown = if n_mode & REQ_MODE_DIFF_UU_STEREO != 0 {
            AB_PARITY_UNKN
        } else {
            AB_PARITY_UNDF
        };

        'body: {
            macro_rules! exit {
                ($r:expr) => {{
                    n_ret = $r;
                    break 'body;
                }};
            }
            let mut b_ignore_isotopic_input_groups = 1;
            if b_taut {
                let t = cs.t_group_info.as_mut().unwrap();
                b_ignore_isotopic_input_groups = t.b_ignore_isotopic;
                t.b_ignore_isotopic = 1;
            }
            let b_ignore_isotopic_input_atoms = cs.b_ignore_isotopic;
            cs.b_ignore_isotopic = 1;

            // 1. non-isotopic
            cs.n_canon_flags = ftcn.n_canon_flags;
            let l = ftcn.n_len_linear_ct as usize;
            cs.linear_ct[..l].copy_from_slice(&ftcn.linear_ct[..l]);
            if !cs.n_num_h.is_empty() && !ftcn.n_num_h.is_empty() {
                for i in 0..num_atoms {
                    cs.n_num_h[i] = (ftcn.n_num_h[i] & 0xFF) as u8 as i8;
                }
            }
            if !cs.n_num_h_fixed.is_empty() && !ftcn.n_num_h_fix_h.is_empty() {
                for i in 0..num_atoms {
                    cs.n_num_h_fixed[i] = (ftcn.n_num_h_fix_h[i] & 0xFF) as u8 as i8;
                }
            }
            cs.n_len_linear_ct = ftcn.n_len_linear_ct;
            cs.n_len_linear_ct_at_only = ftcn.n_len_linear_ct_at_only;
            if !cs.n_symm_rank.is_empty() {
                cs.n_symm_rank[..num_at_tg].copy_from_slice(&ftcn.n_symm_rank_ct[..num_at_tg]);
            }
            if !cs.n_canon_ord.is_empty() {
                cs.n_canon_ord[..num_at_tg]
                    .copy_from_slice(&ftcn.partition_ct.at_number[..num_at_tg]);
                cs.n_len_canon_ord = num_atoms as i32;
            }
            if !ftcn.iso_exchg_atnos.is_empty() && !cs.n_exchg_iso_h.is_empty() {
                for i in 0..num_atoms {
                    cs.n_exchg_iso_h[i] = i8::from(ftcn.iso_exchg_atnos[i] == 0);
                }
            }
            // 2. isotopic
            if b_canon_isotopic {
                if !cs.n_symm_rank_isotopic.is_empty() {
                    cs.n_symm_rank_isotopic[..num_at_tg]
                        .copy_from_slice(&ftcn.n_symm_rank_ct_iso[..num_at_tg]);
                }
                if !cs.n_canon_ord_isotopic.is_empty() {
                    cs.n_canon_ord_isotopic[..num_at_tg]
                        .copy_from_slice(&ftcn.partition_ct_iso.at_number[..num_at_tg]);
                    cs.n_len_canon_ord_isotopic = num_at_tg as i32;
                }
                let max_iso = cs.n_max_len_linear_ct_isotopic;
                n_ret = fill_isotopic_at_linear_ct(
                    num_atoms,
                    at,
                    &ftcn.partition_ct_iso.at_number,
                    &mut cs.linear_ct_isotopic,
                    max_iso,
                    &mut cs.n_len_linear_ct_isotopic,
                );
                if returned_error(n_ret) {
                    break 'body;
                }
                if n_ret < 0 {
                    exit!(CT_TAUCOUNT_ERR);
                }
            } else {
                cs.n_max_len_linear_ct_isotopic = 0;
                cs.n_max_len_linear_ct_isotopic_tautomer = 0;
            }
            if b_taut {
                b_iso_taut = b_iso_taut
                    && !ftcn.partition_ct_iso.rank.is_empty()
                    && !ftcn.partition_ct_iso.at_number.is_empty()
                    && !ftcn.n_symm_rank_ct_iso.is_empty();
                let (max_t, max_it) = (
                    cs.n_max_len_linear_ct_tautomer,
                    cs.n_max_len_linear_ct_isotopic_tautomer,
                );
                n_ret = fill_taut_linear_ct2(
                    num_atoms,
                    num_at_tg,
                    b_iso_taut,
                    &ftcn.partition_ct.rank,
                    &ftcn.partition_ct.at_number,
                    &ftcn.n_symm_rank_ct,
                    &ftcn.partition_ct_iso.at_number,
                    &ftcn.n_symm_rank_ct_iso,
                    &mut cs.linear_ct_tautomer,
                    max_t,
                    &mut cs.n_len_linear_ct_tautomer,
                    &mut cs.linear_ct_isotopic_tautomer,
                    max_it,
                    &mut cs.n_len_linear_ct_isotopic_tautomer,
                    cs.t_group_info.as_mut(),
                );
                if returned_error(n_ret) {
                    break 'body;
                }
                if n_ret <= 0 {
                    exit!(CT_TAUCOUNT_ERR);
                }
                let tgi = cs.t_group_info.as_ref().unwrap();
                let ntg = tgi.num_t_groups;
                let tgn = &tgi.t_group_number;
                if !cs.n_symm_rank_taut.is_empty() {
                    cs.n_symm_rank_taut[..ntg]
                        .copy_from_slice(&tgn[TGSO_SYMM_RANK * ntg..TGSO_SYMM_RANK * ntg + ntg]);
                }
                if !cs.n_canon_ord_taut.is_empty() {
                    cs.n_canon_ord_taut[..ntg].copy_from_slice(&tgn[..ntg]);
                    cs.n_len_canon_ord_taut = ntg as i32;
                }
                if b_canon_isotopic {
                    if !cs.n_symm_rank_isotopic_taut.is_empty() {
                        cs.n_symm_rank_isotopic_taut[..ntg].copy_from_slice(
                            &tgn[TGSO_SYMM_IRANK * ntg..TGSO_SYMM_IRANK * ntg + ntg],
                        );
                    }
                    cs.n_canon_ord_isotopic_taut[..ntg].copy_from_slice(
                        &tgn[TGSO_SYMM_IORDER * ntg..TGSO_SYMM_IORDER * ntg + ntg],
                    );
                    cs.n_len_canon_ord_isotopic_taut = ntg as i32;
                }
            }
            if !cs.linear_ct2.is_empty() {
                let l = cs.n_len_linear_ct as usize;
                let src = cs.linear_ct[..l].to_vec();
                cs.linear_ct2[..l].copy_from_slice(&src);
                cs.n_len_linear_ct2 = cs.n_len_linear_ct;
                cs.n_len_linear_ct_at_only2 = cs.n_len_linear_ct_at_only;
            }
            if num_atoms <= 1 {
                b_canon_stereo = false;
                if num_atoms < 1 || at[0].parity2 == 0 {
                    b_canon_iso_stereo = false;
                }
            }
            if !b_canon_stereo && !(b_canon_isotopic && b_canon_iso_stereo) {
                break 'body;
            }

            // memory allocation
            n_atom_number = vec![0; num_max];
            n_rank = vec![0; num_max];
            n_temp_rank = vec![0; num_max];
            n_symm_rank = vec![0; num_max];
            if neigh_list.is_none() || cs.linear_ct.is_empty() {
                exit!(CT_OUT_OF_RAM);
            }
            cs2 = Cs2::snapshot(cs);
            cs.b_cmp_stereo = 0;
            cs.b_cmp_isotopic_stereo = 0;
            let mut stack2: Vec<Vec<u16>> = Vec::new();
            if b_canon_stereo || b_canon_iso_stereo {
                let mut nn = 2;
                while nn < n_rank_stack_len && nn < st1.len() && !st1[nn].is_empty() {
                    st1[nn][0] = 0;
                    nn += 1;
                }
                stack2 = vec![Vec::new(), Vec::new()];
                let half = (nn - 2) / 2;
                if half > 0 {
                    let from = 2 + half;
                    let moved: Vec<Vec<u16>> = st1.drain(from..nn.min(st1.len())).collect();
                    stack2.extend(moved);
                    st1.truncate(from);
                }
            }
            st2 = Some(stack2);

            if b_canon_stereo {
                n_atom_number[..num_at_tg]
                    .copy_from_slice(&ftcn.partition_ct.at_number[..num_at_tg]);
                n_rank[..num_at_tg].copy_from_slice(&ftcn.partition_ct.rank[..num_at_tg]);
                n_symm_rank[..num_at_tg].copy_from_slice(&ftcn.n_symm_rank_ct[..num_at_tg]);
                n_num_curr_ranks = fix_canon_equivalence_info(
                    num_at_tg,
                    &mut n_symm_rank,
                    &mut n_rank,
                    &mut n_temp_rank,
                    &mut n_atom_number,
                )
                .0;
                cs.n_prev_atom_number[..num_at_tg]
                    .copy_from_slice(&ftcn.partition_ct.at_number[..num_at_tg]);
                n_ret = with_stereo_view(cs, |v| {
                    fill_out_stereo_parities(
                        at,
                        num_atoms,
                        &ftcn.partition_ct.rank,
                        &ftcn.partition_ct.at_number,
                        &n_rank,
                        &n_atom_number,
                        v,
                        false,
                    )
                });
                if returned_error(n_ret) {
                    break 'body;
                }
                if n_ret < 0 {
                    exit!(CT_STEREOCOUNT_ERR);
                }
                // VI. optimize the non-isotopic stereo descriptors
                n_canon_rank_stereo = vec![0; num_max];
                n_symm_stereo = vec![0; num_max + 1];
                if cur_tree.alloc(num_at_tg as i32) < 0 {
                    exit!(CT_OUT_OF_RAM);
                }
                st1[0] = n_rank.clone();
                st1[1] = n_atom_number.clone();
                let s2 = st2.as_mut().unwrap();
                s2[0] = n_rank.clone();
                s2[1] = n_atom_number.clone();
                mark_ranks_stale(&mut st1, n_rank_stack_len);
                mark_ranks_stale(s2, n_rank_stack_len);
                reset_stat_counters(cs);
                n_ret = map_stereo_in_cs(
                    time_is_over,
                    at,
                    num_atoms,
                    num_at_tg,
                    num_max,
                    &ftcn.partition_ct.rank,
                    &ftcn.partition_ct.at_number,
                    &mut n_canon_rank_stereo,
                    &n_symm_rank,
                    &mut st1,
                    s2,
                    &mut n_temp_rank,
                    n_num_curr_ranks,
                    &mut n_symm_stereo,
                    neigh_list.as_mut().unwrap(),
                    cs,
                    &mut cur_tree,
                    v_ab_parity_unknown,
                );
                if returned_error(n_ret) {
                    break 'body;
                }
                let mut b_failed = 0;
                if n_ret == 0 {
                    cs.n_len_linear_ct_stereo_carb = -cs.n_len_linear_ct_stereo_carb.abs();
                    cs2.n_len_carb = cs.n_len_linear_ct_stereo_carb;
                    cs.n_len_linear_ct_stereo_dble = -cs.n_len_linear_ct_stereo_dble.abs();
                    cs2.n_len_dble = cs.n_len_linear_ct_stereo_dble;
                    exit!(CT_STEREOCOUNT_ERR);
                }
                cs2.n_len_dble = cs.n_len_linear_ct_stereo_dble;
                cs2.n_len_carb = cs.n_len_linear_ct_stereo_carb;
                if !cs.n_canon_ord_stereo.is_empty() {
                    b_failed = canon_ord_from_rank(
                        &n_canon_rank_stereo,
                        num_at_tg,
                        &mut cs.n_canon_ord_stereo,
                    );
                    cs.n_len_canon_ord_stereo = if b_failed != 0 {
                        -(num_atoms as i32)
                    } else {
                        num_atoms as i32
                    };
                }
                if b_taut && !cs.n_canon_ord_stereo_taut.is_empty() {
                    let tgi = cs.t_group_info.as_mut().unwrap();
                    n_ret = sort_tautomer_groups_and_endpoints(
                        tgi,
                        num_atoms,
                        num_at_tg,
                        &n_canon_rank_stereo,
                    );
                    if n_ret > 0 {
                        let ntg = tgi.num_t_groups;
                        cs.n_canon_ord_stereo_taut[..ntg]
                            .copy_from_slice(&tgi.t_group_number[..ntg]);
                        cs.n_len_canon_ord_stereo_taut = if b_failed != 0 {
                            -(ntg as i32)
                        } else {
                            ntg as i32
                        };
                    } else if returned_error(n_ret) {
                        break 'body;
                    }
                }

                // VI-B. optimize the inverted stereo descriptors
                n_canon_rank_stereo_inv = vec![0; num_max];
                cs.stereo_slots.dble = Slot::Inv;
                cs.stereo_slots.carb = Slot::Inv;
                cs.n_len_linear_ct_stereo_dble = cs2.n_len_dble;
                cs.n_len_linear_ct_stereo_dble_inv = cs2.n_len_dble;
                cs2.n_len_dble_inv = cs2.n_len_dble;
                cs.n_len_linear_ct_stereo_carb = cs2.n_len_carb;
                cs.n_len_linear_ct_stereo_carb_inv = cs2.n_len_carb;
                cs2.n_len_carb_inv = cs2.n_len_carb;
                if cs.n_len_linear_ct_stereo_dble > 0 || cs.n_len_linear_ct_stereo_carb > 0 {
                    let nd = cs.n_len_linear_ct_stereo_dble.max(0) as usize;
                    let nc = cs.n_len_linear_ct_stereo_carb.max(0) as usize;
                    let src_d = cs.linear_ct_stereo_dble[..nd].to_vec();
                    let src_c = cs.linear_ct_stereo_carb[..nc].to_vec();
                    cs.dble_buf(Slot::Inv)[..nd].copy_from_slice(&src_d);
                    cs.carb_buf(Slot::Inv)[..nc].copy_from_slice(&src_c);
                }
                n_canon_rank_stereo_inv[..num_max].copy_from_slice(&n_canon_rank_stereo[..num_max]);
                if !cs.n_canon_ord_stereo_inv.is_empty() && !cs.n_canon_ord_stereo.is_empty() {
                    let src = cs.n_canon_ord_stereo[..num_at_tg].to_vec();
                    cs.n_canon_ord_stereo_inv[..num_at_tg].copy_from_slice(&src);
                }
                n_ret = with_stereo_view(cs, |v| {
                    invert_stereo(
                        at,
                        num_at_tg,
                        &n_canon_rank_stereo,
                        &mut n_temp_rank,
                        v,
                        true,
                    )
                });
                if returned_error(n_ret) {
                    break 'body;
                } else if n_ret > 0 {
                    mark_ranks_stale(&mut st1, n_rank_stack_len);
                    mark_ranks_stale(s2, n_rank_stack_len);
                    reset_stat_counters(cs);
                    cur_tree.set_pos(0);
                    n_ret = map_stereo_in_cs(
                        time_is_over,
                        at,
                        num_atoms,
                        num_at_tg,
                        num_max,
                        &ftcn.partition_ct.rank,
                        &ftcn.partition_ct.at_number,
                        &mut n_canon_rank_stereo_inv,
                        &n_symm_rank,
                        &mut st1,
                        s2,
                        &mut n_temp_rank,
                        n_num_curr_ranks,
                        &mut n_symm_stereo,
                        neigh_list.as_mut().unwrap(),
                        cs,
                        &mut cur_tree,
                        v_ab_parity_unknown,
                    );
                    if returned_error(n_ret) {
                        break 'body;
                    }
                    let mut b_failed = 0;
                    if n_ret == 0 {
                        cs.n_len_linear_ct_stereo_carb = -cs.n_len_linear_ct_stereo_carb.abs();
                        cs2.n_len_carb = cs.n_len_linear_ct_stereo_carb;
                        cs.n_len_linear_ct_stereo_dble = -cs.n_len_linear_ct_stereo_dble.abs();
                        cs2.n_len_dble = cs.n_len_linear_ct_stereo_dble;
                        exit!(CT_STEREOCOUNT_ERR);
                    }
                    cs.n_len_linear_ct_stereo_dble_inv = cs.n_len_linear_ct_stereo_dble;
                    cs2.n_len_dble_inv = cs.n_len_linear_ct_stereo_dble;
                    cs.n_len_linear_ct_stereo_carb_inv = cs.n_len_linear_ct_stereo_carb;
                    cs2.n_len_carb_inv = cs.n_len_linear_ct_stereo_carb;
                    cs.stereo_slots.dble = Slot::Base;
                    cs.stereo_slots.carb = Slot::Base;
                    cs.n_len_linear_ct_stereo_dble = cs2.n_len_dble;
                    cs.n_len_linear_ct_stereo_carb = cs2.n_len_carb;
                    if cs.n_len_linear_ct_stereo_dble_inv != cs.n_len_linear_ct_stereo_dble
                        || cs.n_len_linear_ct_stereo_carb_inv != cs.n_len_linear_ct_stereo_carb
                    {
                        exit!(CT_CALC_STEREO_ERR);
                    }
                    n_ret = with_stereo_view(cs, |v| {
                        invert_stereo(
                            at,
                            num_at_tg,
                            &n_canon_rank_stereo,
                            &mut n_temp_rank,
                            v,
                            false,
                        )
                    });
                    if returned_error(n_ret) {
                        break 'body;
                    }
                    n_ret = 0;
                    if !cs.n_canon_ord_stereo_inv.is_empty() {
                        b_failed = canon_ord_from_rank(
                            &n_canon_rank_stereo_inv,
                            num_at_tg,
                            &mut cs.n_canon_ord_stereo_inv,
                        );
                        cs.n_len_canon_ord_stereo = if b_failed != 0 {
                            -(num_atoms as i32)
                        } else {
                            num_atoms as i32
                        };
                    }
                    let (di, ci, d, c) = (
                        cs.stereo_slots.dble_inv,
                        cs.stereo_slots.carb_inv,
                        cs.stereo_slots.dble,
                        cs.stereo_slots.carb,
                    );
                    let r = compare_lin_ct_stereo(
                        (!cs.dble_ref(di).is_empty()).then(|| &cs.dble_ref(di)[..]),
                        cs.n_len_linear_ct_stereo_dble_inv,
                        (!cs.carb_ref(ci).is_empty()).then(|| &cs.carb_ref(ci)[..]),
                        cs.n_len_linear_ct_stereo_carb_inv,
                        (!cs.dble_ref(d).is_empty()).then(|| &cs.dble_ref(d)[..]),
                        cs.n_len_linear_ct_stereo_dble,
                        (!cs.carb_ref(c).is_empty()).then(|| &cs.carb_ref(c)[..]),
                        cs.n_len_linear_ct_stereo_carb,
                    );
                    cs.b_cmp_stereo = 2 + r;
                } else if n_ret == 0 {
                    cs.stereo_slots.dble = Slot::Base;
                    cs.stereo_slots.carb = Slot::Base;
                    cs.n_len_linear_ct_stereo_dble = cs2.n_len_dble;
                    cs.n_len_linear_ct_stereo_carb = cs2.n_len_carb;
                }
            }

            if b_taut {
                cs.t_group_info.as_mut().unwrap().b_ignore_isotopic =
                    b_ignore_isotopic_input_groups;
            }
            cs.b_ignore_isotopic = b_ignore_isotopic_input_atoms;

            if b_canon_iso_stereo && b_canon_isotopic {
                // VII. optimize the isotopic stereo descriptors
                n_atom_number[..num_at_tg]
                    .copy_from_slice(&ftcn.partition_ct_iso.at_number[..num_at_tg]);
                n_rank[..num_at_tg].copy_from_slice(&ftcn.partition_ct_iso.rank[..num_at_tg]);
                n_symm_rank[..num_at_tg].copy_from_slice(&ftcn.n_symm_rank_ct_iso[..num_at_tg]);
                n_num_curr_ranks = fix_canon_equivalence_info(
                    num_at_tg,
                    &mut n_symm_rank,
                    &mut n_rank,
                    &mut n_temp_rank,
                    &mut n_atom_number,
                )
                .0;
                cs.n_prev_atom_number[..num_at_tg]
                    .copy_from_slice(&ftcn.partition_ct_iso.at_number[..num_at_tg]);
                n_canon_rank_isotopic_stereo = vec![0; num_max];
                if n_symm_stereo.is_empty() {
                    n_symm_stereo = vec![0; num_max + 1];
                }
                if cur_tree.alloc(num_at_tg as i32) != 0 {
                    exit!(CT_OUT_OF_RAM);
                }
                st1[0] = n_rank.clone();
                st1[1] = n_atom_number.clone();
                let s2 = st2.as_mut().unwrap();
                s2[0] = n_rank.clone();
                s2[1] = n_atom_number.clone();

                switch_atom_stereo_and_isotopic_stereo(
                    at,
                    num_atoms,
                    &mut b_switched_atom_to_isotopic,
                );
                set_ct_to_isotopic_stereo(cs, &cs2);
                n_ret = with_stereo_view(cs, |v| {
                    fill_out_stereo_parities(
                        at,
                        num_atoms,
                        &ftcn.partition_ct_iso.rank,
                        &ftcn.partition_ct_iso.at_number,
                        &n_rank,
                        &n_atom_number,
                        v,
                        true,
                    )
                });
                if returned_error(n_ret) {
                    break 'body;
                }
                'iso: {
                    if n_ret == 0 {
                        cs.n_len_linear_ct_isotopic_stereo_dble = 0;
                        cs2.n_len_iso_dble = 0;
                        cs.n_len_linear_ct_isotopic_stereo_carb = 0;
                        cs2.n_len_iso_carb = 0;
                        cs.n_len_canon_ord_isotopic_stereo = 0;
                        cs.n_len_canon_ord_isotopic_stereo_taut = 0;
                        cs.n_len_linear_ct_isotopic_stereo_dble_inv = 0;
                        cs2.n_len_iso_dble_inv = 0;
                        cs.n_len_linear_ct_isotopic_stereo_carb_inv = 0;
                        cs2.n_len_iso_carb_inv = 0;
                        break 'iso;
                    }
                    // VII-A. non-inverted isotopic stereo
                    mark_ranks_stale(&mut st1, n_rank_stack_len);
                    mark_ranks_stale(s2, n_rank_stack_len);
                    reset_stat_counters(cs);
                    n_ret = map_stereo_in_cs(
                        time_is_over,
                        at,
                        num_atoms,
                        num_at_tg,
                        num_max,
                        &ftcn.partition_ct_iso.rank,
                        &ftcn.partition_ct_iso.at_number,
                        &mut n_canon_rank_isotopic_stereo,
                        &n_symm_rank,
                        &mut st1,
                        s2,
                        &mut n_temp_rank,
                        n_num_curr_ranks,
                        &mut n_symm_stereo,
                        neigh_list.as_mut().unwrap(),
                        cs,
                        &mut cur_tree,
                        v_ab_parity_unknown,
                    );
                    if returned_error(n_ret) {
                        break 'body;
                    }
                    let mut b_failed = 0;
                    if n_ret == 0 {
                        cs.n_len_linear_ct_isotopic_stereo_dble =
                            -cs.n_len_linear_ct_stereo_dble.abs();
                        cs2.n_len_iso_dble = cs.n_len_linear_ct_isotopic_stereo_dble;
                        cs.n_len_linear_ct_isotopic_stereo_carb =
                            -cs.n_len_linear_ct_stereo_carb.abs();
                        cs2.n_len_iso_carb = cs.n_len_linear_ct_isotopic_stereo_carb;
                        exit!(CT_STEREOCOUNT_ERR);
                    }
                    cs.n_len_linear_ct_isotopic_stereo_dble = cs.n_len_linear_ct_stereo_dble;
                    cs2.n_len_iso_dble = cs.n_len_linear_ct_stereo_dble;
                    cs.n_len_linear_ct_isotopic_stereo_carb = cs.n_len_linear_ct_stereo_carb;
                    cs2.n_len_iso_carb = cs.n_len_linear_ct_stereo_carb;
                    if !cs.n_canon_ord_isotopic_stereo.is_empty() {
                        b_failed = canon_ord_from_rank(
                            &n_canon_rank_isotopic_stereo,
                            num_at_tg,
                            &mut cs.n_canon_ord_isotopic_stereo,
                        );
                        cs.n_len_canon_ord_isotopic_stereo = if b_failed != 0 {
                            -(num_atoms as i32)
                        } else {
                            num_atoms as i32
                        };
                    }
                    if !cs.n_canon_ord_isotopic_stereo_taut.is_empty() {
                        if let Some(tgi) = cs.t_group_info.as_mut().filter(|_| b_taut) {
                            n_ret = sort_tautomer_groups_and_endpoints(
                                tgi,
                                num_atoms,
                                num_at_tg,
                                &n_canon_rank_isotopic_stereo,
                            );
                            if n_ret > 0 {
                                let ntg = tgi.num_t_groups;
                                cs.n_canon_ord_isotopic_stereo_taut[..ntg]
                                    .copy_from_slice(&tgi.t_group_number[..ntg]);
                                cs.n_len_canon_ord_isotopic_stereo_taut = if b_failed != 0 {
                                    -(ntg as i32)
                                } else {
                                    ntg as i32
                                };
                            } else if returned_error(n_ret) {
                                break 'body;
                            }
                        }
                    }

                    // VII-B. inverted isotopic stereo
                    n_canon_rank_isotopic_stereo_inv = vec![0; num_max];
                    cs.stereo_slots.dble = Slot::IsoInv;
                    cs.stereo_slots.carb = Slot::IsoInv;
                    cs.n_len_linear_ct_stereo_dble = cs2.n_len_iso_dble;
                    cs.n_len_linear_ct_stereo_dble_inv = cs2.n_len_iso_dble;
                    cs2.n_len_iso_dble_inv = cs2.n_len_iso_dble;
                    cs.n_len_linear_ct_stereo_carb = cs2.n_len_iso_carb;
                    cs.n_len_linear_ct_stereo_carb_inv = cs2.n_len_iso_carb;
                    cs2.n_len_iso_carb_inv = cs2.n_len_iso_carb;
                    if cs.n_len_linear_ct_stereo_dble > 0 || cs.n_len_linear_ct_stereo_carb > 0 {
                        let nd = cs.n_len_linear_ct_stereo_dble.max(0) as usize;
                        let nc = cs.n_len_linear_ct_stereo_carb.max(0) as usize;
                        let src_d = cs.linear_ct_isotopic_stereo_dble[..nd].to_vec();
                        let src_c = cs.linear_ct_isotopic_stereo_carb[..nc].to_vec();
                        cs.dble_buf(Slot::IsoInv)[..nd].copy_from_slice(&src_d);
                        cs.carb_buf(Slot::IsoInv)[..nc].copy_from_slice(&src_c);
                    }
                    n_canon_rank_isotopic_stereo_inv[..num_max]
                        .copy_from_slice(&n_canon_rank_isotopic_stereo[..num_max]);
                    if !cs.n_canon_ord_isotopic_stereo_inv.is_empty()
                        && !cs.n_canon_ord_isotopic_stereo.is_empty()
                    {
                        let src = cs.n_canon_ord_isotopic_stereo[..num_at_tg].to_vec();
                        cs.n_canon_ord_isotopic_stereo_inv[..num_at_tg].copy_from_slice(&src);
                    }
                    n_ret = with_stereo_view(cs, |v| {
                        invert_stereo(
                            at,
                            num_at_tg,
                            &n_canon_rank_isotopic_stereo,
                            &mut n_temp_rank,
                            v,
                            true,
                        )
                    });
                    if returned_error(n_ret) {
                        break 'body;
                    } else if n_ret > 0 {
                        mark_ranks_stale(&mut st1, n_rank_stack_len);
                        mark_ranks_stale(s2, n_rank_stack_len);
                        reset_stat_counters(cs);
                        n_ret = map_stereo_in_cs(
                            time_is_over,
                            at,
                            num_atoms,
                            num_at_tg,
                            num_max,
                            &ftcn.partition_ct_iso.rank,
                            &ftcn.partition_ct_iso.at_number,
                            &mut n_canon_rank_isotopic_stereo_inv,
                            &n_symm_rank,
                            &mut st1,
                            s2,
                            &mut n_temp_rank,
                            n_num_curr_ranks,
                            &mut n_symm_stereo,
                            neigh_list.as_mut().unwrap(),
                            cs,
                            &mut cur_tree,
                            v_ab_parity_unknown,
                        );
                        if returned_error(n_ret) {
                            break 'body;
                        }
                        let mut b_failed = 0;
                        if n_ret == 0 {
                            cs.n_len_linear_ct_isotopic_stereo_dble =
                                -cs.n_len_linear_ct_stereo_dble.abs();
                            cs2.n_len_iso_dble = cs.n_len_linear_ct_isotopic_stereo_dble;
                            cs.n_len_linear_ct_isotopic_stereo_carb =
                                -cs.n_len_linear_ct_stereo_carb.abs();
                            cs2.n_len_iso_carb = cs.n_len_linear_ct_isotopic_stereo_carb;
                            exit!(CT_STEREOCOUNT_ERR);
                        }
                        cs.n_len_linear_ct_isotopic_stereo_dble_inv =
                            cs.n_len_linear_ct_stereo_dble;
                        cs2.n_len_iso_dble_inv = cs.n_len_linear_ct_stereo_dble;
                        cs.n_len_linear_ct_isotopic_stereo_carb_inv =
                            cs.n_len_linear_ct_stereo_carb;
                        cs2.n_len_iso_carb_inv = cs.n_len_linear_ct_stereo_carb;
                        cs.stereo_slots.dble = Slot::Iso;
                        cs.stereo_slots.carb = Slot::Iso;
                        cs.n_len_linear_ct_stereo_dble = cs2.n_len_iso_dble;
                        cs.n_len_linear_ct_stereo_carb = cs2.n_len_iso_carb;
                        if cs.n_len_linear_ct_isotopic_stereo_dble_inv
                            != cs.n_len_linear_ct_isotopic_stereo_dble
                            || cs.n_len_linear_ct_isotopic_stereo_carb_inv
                                != cs.n_len_linear_ct_isotopic_stereo_carb
                        {
                            exit!(CT_CALC_STEREO_ERR);
                        }
                        n_ret = with_stereo_view(cs, |v| {
                            invert_stereo(
                                at,
                                num_at_tg,
                                &n_canon_rank_isotopic_stereo,
                                &mut n_temp_rank,
                                v,
                                false,
                            )
                        });
                        if returned_error(n_ret) {
                            break 'body;
                        }
                        n_ret = 0;
                        if !cs.n_canon_ord_isotopic_stereo_inv.is_empty() {
                            b_failed = canon_ord_from_rank(
                                &n_canon_rank_isotopic_stereo_inv,
                                num_at_tg,
                                &mut cs.n_canon_ord_isotopic_stereo_inv,
                            );
                            cs.n_len_canon_ord_isotopic_stereo = if b_failed != 0 {
                                -(num_atoms as i32)
                            } else {
                                num_atoms as i32
                            };
                        }
                        let r = compare_lin_ct_stereo(
                            (!cs.linear_ct_isotopic_stereo_dble_inv.is_empty())
                                .then(|| &cs.linear_ct_isotopic_stereo_dble_inv[..]),
                            cs.n_len_linear_ct_isotopic_stereo_dble_inv,
                            (!cs.linear_ct_isotopic_stereo_carb_inv.is_empty())
                                .then(|| &cs.linear_ct_isotopic_stereo_carb_inv[..]),
                            cs.n_len_linear_ct_isotopic_stereo_carb_inv,
                            (!cs.linear_ct_isotopic_stereo_dble.is_empty())
                                .then(|| &cs.linear_ct_isotopic_stereo_dble[..]),
                            cs.n_len_linear_ct_isotopic_stereo_dble,
                            (!cs.linear_ct_isotopic_stereo_carb.is_empty())
                                .then(|| &cs.linear_ct_isotopic_stereo_carb[..]),
                            cs.n_len_linear_ct_isotopic_stereo_carb,
                        );
                        cs.b_cmp_isotopic_stereo = 2 + r;
                    } else if n_ret == 0 {
                        cs.stereo_slots.dble = Slot::Iso;
                        cs.stereo_slots.carb = Slot::Iso;
                        cs.n_len_linear_ct_stereo_dble = cs2.n_len_iso_dble;
                        cs.n_len_linear_ct_stereo_carb = cs2.n_len_iso_carb;
                    }
                }
            }
        }
    }
    // exit_function
    if b_switched_atom_to_isotopic {
        switch_atom_stereo_and_isotopic_stereo(at, num_atoms, &mut b_switched_atom_to_isotopic);
        set_ct_to_non_isotopic_stereo(cs, &cs2);
    }
    if !cs.linear_ct2.is_empty() {
        std::mem::swap(&mut cs.linear_ct, &mut cs.linear_ct2);
        std::mem::swap(&mut cs.n_len_linear_ct, &mut cs.n_len_linear_ct2);
        std::mem::swap(
            &mut cs.n_len_linear_ct_at_only,
            &mut cs.n_len_linear_ct_at_only2,
        );
    }
    // give the arrays of the second stack back to the first one
    st1[0] = Vec::new();
    st1[1] = Vec::new();
    if let Some(s2) = st2.as_mut() {
        let mut i = 2;
        while i < n_rank_stack_len && i < st1.len() && !st1[i].is_empty() {
            i += 1;
        }
        st1.truncate(i);
        let mut n = 2;
        while n < n_rank_stack_len && n < s2.len() && !s2[n].is_empty() {
            if i + 1 < n_rank_stack_len {
                st1.push(std::mem::take(&mut s2[n]));
                i += 1;
            }
            n += 1;
        }
    }
    if st1.len() % 2 == 1 {
        st1.push(Vec::new());
    }
    bcn.p_rank_stack = st1
        .chunks(2)
        .map(|c| super::ichican2::Partition {
            rank: c[0].clone(),
            at_number: c[1].clone(),
        })
        .collect();
    bcn.ftcn[b_taut_ftcn].neigh_list = neigh_list.take();
    cs.neigh_list = None;
    cs.p_bcn = Some(bcn);
    if n_ret >= -1 { num_atoms as i32 } else { n_ret }
}

/// `Canon_INChI`
pub fn canon_inchi(
    time_is_over: &mut dyn FnMut() -> bool,
    num_atoms: usize,
    num_at_tg: usize,
    at: &mut [SpAtom],
    cs: &mut CanonStat,
    n_mode: i64,
    b_taut_ftcn: usize,
) -> i32 {
    if cs.p_bcn.is_some() && cs.neigh_list.is_none() {
        return canon_inchi3(
            time_is_over,
            num_atoms,
            num_at_tg,
            at,
            cs,
            n_mode,
            b_taut_ftcn,
        );
    }
    CT_CANON_ERR
}
