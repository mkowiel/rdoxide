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
