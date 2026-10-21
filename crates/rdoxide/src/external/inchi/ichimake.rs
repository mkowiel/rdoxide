//! Construction of the `INChI` output structures (`INCHI_BASE/src/ichimake.c`):
//! the input-atom to `sp_ATOM` copy and the Hill formula comparison.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::{InChI, InChIAux, InChIStereo};
use super::ichister::SpAtom;
use super::inpdef::InpAtom;
use super::util::get_periodic_table_number;
use std::cmp::Ordering;

fn nul_str(b: &[u8]) -> &str {
    let n = b.iter().position(|&c| c == 0).unwrap_or(b.len());
    std::str::from_utf8(&b[..n]).unwrap_or("")
}

/// `inp2spATOM`: copies the structure fields of the input atoms into `sp_ATOM`
/// records (everything else is zeroed).
pub fn inp2sp_atom(inp_at: &[InpAtom], at: &mut [SpAtom]) {
    for (i, ia) in inp_at.iter().enumerate() {
        let mut a = SpAtom {
            elname: ia.elname,
            ..SpAtom::default()
        };
        a.el_number = get_periodic_table_number(nul_str(&a.elname)) as u8;
        a.valence = ia.valence;
        let val = ia.valence.max(0) as usize;
        for j in 0..val {
            a.neighbor[j] = ia.neighbor[j];
            a.bond_type[j] = ia.bond_type[j];
        }
        a.chem_bonds_valence = ia.chem_bonds_valence;
        a.orig_at_number = ia.orig_at_number;
        a.orig_compt_at_numb = ia.orig_compt_at_numb;
        a.endpoint = ia.endpoint;
        a.iso_atw_diff = ia.iso_atw_diff;
        a.num_h = ia.num_h;
        a.c_flags = ia.c_flags;
        a.num_iso_h = ia.num_iso_h;
        a.charge = ia.charge;
        a.radical = ia.radical;
        a.n_block_system = ia.n_block_system;
        a.b_cut_vertex = ia.b_cut_vertex;
        a.n_ring_system = ia.n_ring_system;
        a.n_num_at_in_ring_system = ia.n_num_at_in_ring_system;
        at[i] = a;
    }
}

/// `GetElementAndCount`: reads the next element and its count from a Hill
/// formula; returns `(ret, element, count)` with `ret` = 1 (read), 0 (end of
/// the formula: "Zzz", 99999) or -1 (not a formula). `f` is advanced.
pub fn get_element_and_count(f: &mut &[u8]) -> (i32, String, i32) {
    let p = *f;
    if let Some(&c0) = p.first().filter(|&&c| c != 0) {
        if c0.is_ascii_uppercase() {
            let mut el = vec![c0];
            let mut k = 1;
            if p.get(1).is_some_and(|c| c.is_ascii_lowercase()) {
                el.push(p[1]);
                k = 2;
            }
            if el == b"C" {
                el[0] = b'A'; // less than any element: carbon compounds first
            }
            let mut q = k;
            let mut count = 1;
            if p.get(q).is_some_and(|c| c.is_ascii_digit()) {
                let mut v: i64 = 0;
                while let Some(d) = p.get(q).filter(|c| c.is_ascii_digit()) {
                    v = (v * 10 + i64::from(d - b'0')).min(i64::from(i32::MAX));
                    q += 1;
                }
                count = v as i32;
            }
            *f = &p[q..];
            return (1, String::from_utf8(el).unwrap_or_default(), count);
        }
        return (-1, String::new(), 0);
    }
    (0, "Zzz".to_string(), 99999)
}

/// `CompareHillFormulas`
pub fn compare_hill_formulas(f1: &str, f2: &str) -> i32 {
    let (mut f1, mut f2) = (f1.as_bytes(), f2.as_bytes());
    loop {
        let (ret1, el1, count1) = get_element_and_count(&mut f1);
        let (ret2, el2, count2) = get_element_and_count(&mut f2);
        if ret1 >= 0 && ret2 >= 0 {
            match el1.cmp(&el2) {
                Ordering::Less => return -1,
                Ordering::Greater => return 1,
                Ordering::Equal => {}
            }
            if count2 != count1 {
                return count2 - count1;
            }
        } else {
            return 0;
        }
        if !(ret1 > 0 && ret2 > 0) {
            return 0;
        }
    }
}

/// `CompareHillFormulasNoH`: like `compare_hill_formulas` ignoring the `H`
/// element; the numbers of H are added to `num_h1`/`num_h2`.
pub fn compare_hill_formulas_no_h(f1: &str, f2: &str, num_h1: &mut i32, num_h2: &mut i32) -> i32 {
    let (mut f1, mut f2) = (f1.as_bytes(), f2.as_bytes());
    loop {
        let (mut ret1, mut el1, mut count1) = get_element_and_count(&mut f1);
        if ret1 > 0 && el1 == "H" {
            *num_h1 += count1;
            (ret1, el1, count1) = get_element_and_count(&mut f1);
        }
        let (mut ret2, mut el2, mut count2) = get_element_and_count(&mut f2);
        if ret2 > 0 && el2 == "H" {
            *num_h2 += count2;
            (ret2, el2, count2) = get_element_and_count(&mut f2);
        }
        if ret1 >= 0 && ret2 >= 0 {
            match el1.cmp(&el2) {
                Ordering::Less => return -1,
                Ordering::Greater => return 1,
                Ordering::Equal => {}
            }
            if count2 != count1 {
                return count2 - count1;
            }
        } else {
            return 0;
        }
        if !(ret1 > 0 && ret2 > 0) {
            return 0;
        }
    }
}

/// `CT_MODE_PREDECESSORS`
pub const CT_MODE_PREDECESSORS: i32 = 8;

/// `mystrrev`: reverses a string in place (up to the first NUL)
pub fn mystrrev(p: &mut [u8]) {
    let n = p.iter().position(|&c| c == 0).unwrap_or(p.len());
    p[..n].reverse();
}

/// `GetDfsOrder4CT`: the DFS order of the connection table output as triplets
/// (canonical number, 16 + number of H, delimiter), zero terminated; `None`
/// for an inconsistent table or when the output buffer would overflow.
pub fn get_dfs_order_4ct(
    linear_ct: &[u16],
    n_len_ct: usize,
    n_num_h: Option<&[i8]>,
    num_atoms: usize,
    n_ct_mode: i32,
) -> Option<Vec<u16>> {
    use super::ichimap::create_neigh_list_from_linear_ct;
    use super::ichisort::insertions_sort;
    use super::inpdef::MAX_ATOMS;
    const DELETED: u16 = MAX_ATOMS as u16 + 1;
    let ct_pred = n_ct_mode & CT_MODE_PREDECESSORS != 0;
    let mut nl = create_neigh_list_from_linear_ct(linear_ct, n_len_ct, num_atoms)?;
    let mut stack = vec![0u16; num_atoms];
    let mut num_desc = vec![0u16; num_atoms];
    let mut dfs_number = vec![0u16; num_atoms];
    let mut neigh_numb = vec![0u8; num_atoms];
    let mut start = 0usize;
    if !ct_pred {
        for i in 1..num_atoms {
            if nl[i][0] < nl[start][0] {
                start = i;
            }
        }
    }
    let mut u = start;
    let mut n_dfs: u16 = 0;
    let mut top: i32 = -1;
    n_dfs += 1;
    dfs_number[u] = n_dfs;
    num_desc[u] = u16::from(!ct_pred);
    top += 1;
    stack[top as usize] = u as u16;
    let mut num_rings = 0usize;
    loop {
        let mut i;
        loop {
            i = stack[top as usize] as usize;
            let j = usize::from(neigh_numb[i]) + 1;
            if (nl[i][0] as usize) < j {
                break;
            }
            neigh_numb[i] += 1;
            u = nl[i][j] as usize;
            if dfs_number[u] == 0 {
                top += 1;
                stack[top as usize] = u as u16;
                n_dfs += 1;
                dfs_number[u] = n_dfs;
                if ct_pred {
                    num_desc[u] = (i + 1) as u16;
                } else {
                    num_desc[u] = num_desc[u].wrapping_add(1);
                }
            } else if top != 0
                && u != stack[top as usize - 1] as usize
                && dfs_number[u] < dfs_number[i]
            {
                if !ct_pred {
                    num_desc[i] = num_desc[i].wrapping_add(1);
                }
                num_rings += 1;
            } else {
                nl[i][j] = DELETED;
            }
        }
        neigh_numb[i] = 0;
        if !ct_pred && top != 0 {
            let p = stack[top as usize - 1] as usize;
            num_desc[p] = num_desc[p].wrapping_add(num_desc[i]);
        }
        top -= 1;
        if top < 0 {
            break;
        }
    }
    for i in 0..num_atoms {
        if nl[i][0] > 1 {
            let cur = i;
            let cnt = nl[i][0] as usize;
            let (dn, nd) = (&dfs_number, &num_desc);
            insertions_sort(&mut nl[i][1..=cnt], |a, b| {
                let (n1, n2) = (i32::from(*a), i32::from(*b));
                if n1 > MAX_ATOMS as i32 {
                    return i32::from(n2 <= MAX_ATOMS as i32);
                }
                if n2 > MAX_ATOMS as i32 {
                    return -1;
                }
                let cur_dfs = dn[cur];
                let d1 = if cur_dfs > dn[n1 as usize] {
                    0
                } else {
                    i32::from(nd[n1 as usize])
                };
                let d2 = if cur_dfs > dn[n2 as usize] {
                    0
                } else {
                    i32::from(nd[n2 as usize])
                };
                if d1 != d2 {
                    return d1 - d2;
                }
                n1 - n2
            });
        }
        let cnt = nl[i][0] as usize;
        let mut k = 0;
        while k < cnt && nl[i][k + 1] <= MAX_ATOMS as u16 {
            k += 1;
        }
        nl[i][0] = k as u16;
    }
    let tot_len = 3 * (num_atoms + num_rings + 1);
    let mut out = vec![0u16; tot_len];
    let h16 = |u: usize| n_num_h.map_or(0, |h| (16 + i32::from(h[u])) as u16);
    if ct_pred {
        let mut k: i64 = -3;
        for u in 0..num_atoms {
            k += 3;
            if k as usize + 6 > tot_len {
                return None;
            }
            let ku = k as usize;
            out[ku] = if num_desc[u] != 0 {
                num_desc[u]
            } else {
                DELETED
            };
            out[ku + 1] = h16(u);
            out[ku + 2] = if k != 0 { u16::from(b',') } else { 0 };
            let mut j = 1;
            while j <= nl[u][0] as usize && dfs_number[u] > dfs_number[nl[u][j] as usize] {
                let i = nl[u][j] as usize;
                k += 3;
                if k as usize + 6 > tot_len {
                    return None;
                }
                let ku = k as usize;
                out[ku] = (i + 1) as u16;
                out[ku + 1] = 0;
                out[ku + 2] = u16::from(b'-');
                j += 1;
            }
        }
    } else {
        let mut u = start;
        top = 0;
        neigh_numb.iter_mut().for_each(|x| *x = 0);
        stack[0] = u as u16;
        let mut k = 0usize;
        out[k] = (u + 1) as u16;
        out[k + 1] = h16(u);
        out[k + 2] = 0;
        loop {
            let mut i;
            loop {
                i = stack[top as usize] as usize;
                let j = usize::from(neigh_numb[i]) + 1;
                if !(i < num_atoms && nl[i][0] as usize >= j) {
                    break;
                }
                k += 3;
                if k + 6 > tot_len {
                    return None;
                }
                neigh_numb[i] += 1;
                u = nl[i][j] as usize;
                out[k] = (u + 1) as u16;
                if dfs_number[u] > dfs_number[i] {
                    top += 1;
                    stack[top as usize] = u as u16;
                    out[k + 1] = h16(u);
                } else {
                    out[k + 1] = 0;
                }
                let n = nl[i][0] as usize;
                let delim = if n > 1 {
                    if j == 1 {
                        b'('
                    } else if j == n {
                        b')'
                    } else {
                        b','
                    }
                } else {
                    b'-'
                };
                out[k + 2] = u16::from(delim);
            }
            neigh_numb[i] = 0;
            top -= 1;
            if top < 0 {
                break;
            }
        }
    }
    Some(out)
}

/// `CheckCanonNumberingCorrectness`: recomputes the connection table with the
/// final (non-isotopic, then isotopic) numbering; `CT_CANON_ERR` if a warning
/// or error bit was set, else 0.
pub fn check_canon_numbering_correctness(
    num_atoms: usize,
    num_at_tg: usize,
    at: &[SpAtom],
    cs: &mut super::ichicano::CanonStat,
) -> i32 {
    use super::ichicano::update_full_linear_ct;
    use super::ichierr::CT_CANON_ERR;
    const WARN_FAILED_STEREO: i32 = 0x0001;
    const WARN_FAILED_ISOTOPIC: i32 = 0x0002;
    const WARN_FAILED_ISOTOPIC_STEREO: i32 = 0x0004;
    const ERR_NO_CANON_RESULTS: i32 = 0x0008;
    let mut n_error_code = 0;
    let mut rank = vec![0u16; num_at_tg + 1];
    let ord = if cs.n_len_canon_ord_stereo > 0 {
        Some(cs.n_canon_ord_stereo.clone())
    } else if cs.n_len_canon_ord > 0 {
        Some(cs.n_canon_ord.clone())
    } else {
        None
    };
    match ord {
        Some(o) => {
            for i in 0..num_at_tg {
                rank[o[i] as usize] = (i + 1) as u16;
            }
            if update_full_linear_ct(num_atoms, num_at_tg, at, &rank, &o, cs, false) != 0 {
                n_error_code |= WARN_FAILED_STEREO;
            }
        }
        None => n_error_code |= ERR_NO_CANON_RESULTS,
    }
    if n_error_code & ERR_NO_CANON_RESULTS == 0 {
        let ord = if cs.n_len_canon_ord_isotopic_stereo > 0 {
            Some(cs.n_canon_ord_isotopic_stereo.clone())
        } else if cs.n_len_canon_ord_isotopic > 0 {
            Some(cs.n_canon_ord_isotopic.clone())
        } else {
            None
        };
        if let Some(o) = ord {
            for i in 0..num_at_tg {
                rank[o[i] as usize] = (i + 1) as u16;
            }
            if update_full_linear_ct(num_atoms, num_at_tg, at, &rank, &o, cs, false) != 0 {
                n_error_code |= if cs.n_len_canon_ord_isotopic_stereo != 0 {
                    WARN_FAILED_ISOTOPIC_STEREO
                } else {
                    WARN_FAILED_ISOTOPIC
                };
            }
        }
    }
    if n_error_code != 0 { CT_CANON_ERR } else { 0 }
}

/// `INP_ATOM_DATA`: the normalized structure of one kind (tautomeric or not)
#[derive(Debug, Clone, Default)]
pub struct InpAtomData {
    pub at: Vec<InpAtom>,
    /// `at_fixed_bonds`: tautomeric case, added or removed H
    pub at_fixed_bonds: Option<Vec<InpAtom>>,
    pub num_at: i32,
    pub num_removed_h: i32,
    pub num_bonds: i32,
    pub num_isotopic: i32,
    pub b_exists: i32,
    pub b_deleted: i32,
    pub b_has_isotopic_layer: i32,
    pub b_tautomeric: i32,
    pub b_taut_preprocessed: i32,
    pub n_num_removed_protons: i32,
    pub n_num_removed_protons_isotopic: [i16; 3],
    pub num_iso_h: [i16; 3],
    pub b_taut_flags: u64,
    pub b_taut_flags_done: u64,
    pub b_normalization_flags: u64,
}

/// `Create_INChI`: normalizes the structure, finds the tautomeric groups and
/// stereo, canonicalizes and fills `p_inchi[..]` / `p_aux[..]` (both slots
/// `TAUT_NON` and `TAUT_YES` of the result arrays must be allocated by the
/// caller when requested). Returns the number of atoms or a negative error code;
/// the t-group info is returned in `ti_out` when given. The time-out checks,
/// the polymer end-bond list and the debug extraction flags of the C code
/// are not ported.
#[allow(clippy::too_many_arguments, unused_assignments)]
pub fn create_inchi(
    p_inchi: &mut [Option<super::ichi::InChI>; 2],
    p_aux: &mut [Option<super::ichi::InChIAux>; 2],
    inp_at: &[InpAtom],
    out_norm_data: &mut [InpAtomData; 2],
    num_inp_at: usize,
    mut n_user_mode: u64,
    pb_taut_flags: &mut u64,
    pb_taut_flags_done: &mut u64,
    ti_out: Option<&mut super::ichitaut::TGroupInfo>,
    str_err: &mut String,
    ip: &super::ichiparm::InputParms,
) -> i32 {
    use super::ichican2::{
        AtomSizes, Bcn, TAUT_NON, TAUT_YES, get_base_canon_ranking, set_atom_iso_sort_keys,
    };
    use super::ichicano::{
        CMODE_CT, CMODE_ISO, CMODE_ISO_STEREO, CMODE_NOEQ_STEREO, CMODE_STEREO, CMODE_TAUT,
        CanonStat, allocate_cs, canon_inchi, get_canon_lengths,
    };
    use super::ichierr::{CT_ERR_MAX, CT_ERR_MIN, CT_TAUCOUNT_ERR};
    use super::ichimak2::fill_out_inchi;
    use super::ichiparm::*;
    use super::ichister::{
        AB_PARITY_UNDF, AB_PARITY_UNKN, PES_BIT_ARSINE_STEREO, PES_BIT_FIX_SP3_BUG,
        PES_BIT_PHOSPHINE_STEREO, PES_BIT_POINT_EDGE_STEREO, set_stereo_parity,
    };
    use super::ichitaut::{
        FLAG_NORM_CONSIDER_TAUT, TGroupInfo, count_tautomer_groups, mark_alt_bonds_and_taut_groups,
        set_tautomer_iso_sort_keys,
    };
    use super::inpdef::ATT_PROTON;
    use super::runichi3::{TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE};
    use super::strutil::{add_dt_to_num_h, remove_terminal_hdt};
    const CMODE_REDNDNT: i64 = 0x80;
    const CMODE_NO_ALT: i64 = 0x100;
    const CMODE_REL: i64 = 0x200;
    const CMODE_RAC: i64 = 0x400;
    const CMODE_SC_UU: i64 = 0x800;
    const CMODE_SB_UU: i64 = 0x1000;
    const CANON_MODE_ISO: i64 = CMODE_CT | CMODE_ISO | 0x04;
    const CANON_MODE_STEREO: i64 = CMODE_CT | CMODE_STEREO;
    const CANON_MODE_ISO_STEREO: i64 = CMODE_CT | CMODE_ISO | 0x04 | CMODE_ISO_STEREO;
    const TG_FLAG_ALL_TAUTOMERIC: u64 = TG_FLAG_TEST_TAUT__ATOMS
        | TG_FLAG_TEST_TAUT__SALTS
        | TG_FLAG_TEST_TAUT2_SALTS
        | TG_FLAG_MERGE_TAUT_SALTS;
    let returned_error = |r: i32| (CT_ERR_MIN..=CT_ERR_MAX).contains(&r);
    let alt_taut = |x: usize| if x > TAUT_YES { TAUT_YES } else { 1 - x };
    let _ = str_err;

    let mut s = [AtomSizes::default(), AtomSizes::default()];
    let mut at: [Option<Vec<super::ichister::SpAtom>>; 2] = [None, None];
    let mut tgi = TGroupInfo::default();
    let mut bcn = Bcn::default();
    let mut ret = 0i32;
    let mut num_at_tg = 0usize;
    let mut num_atoms = 0usize;
    let mut num_taut_at = 0usize;
    let mut num_removed_h_taut = 0usize;
    let mut b_has_isotopic_atoms;
    let b_pointed_edge_stereo = (if *pb_taut_flags & TG_FLAG_POINTED_EDGE_STEREO != 0 {
        PES_BIT_POINT_EDGE_STEREO
    } else {
        0
    }) | (if *pb_taut_flags & TG_FLAG_PHOSPHINE_STEREO != 0 {
        PES_BIT_PHOSPHINE_STEREO
    } else {
        0
    }) | (if *pb_taut_flags & TG_FLAG_ARSINE_STEREO != 0 {
        PES_BIT_ARSINE_STEREO
    } else {
        0
    }) | (if *pb_taut_flags & TG_FLAG_FIX_SP3_BUG != 0 {
        PES_BIT_FIX_SP3_BUG
    } else {
        0
    });
    let b_taut_flags = *pb_taut_flags & !TG_FLAG_ALL_TAUTOMERIC;
    let mut b_taut_flags_done = *pb_taut_flags_done;
    let b_fix_iso_fixed_h = *pb_taut_flags & TG_FLAG_FIX_ISO_FIXEDH_BUG != 0;
    let b_fix_term_h_chrg = *pb_taut_flags & TG_FLAG_FIX_TERM_H_CHRG_BUG != 0;
    let large_molecules = ip.b_large_molecules != 0;
    let v_ab_parity_unknown = if n_user_mode & REQ_MODE_DIFF_UU_STEREO != 0 {
        AB_PARITY_UNKN
    } else {
        AB_PARITY_UNDF
    };

    'body: {
        for i in 0..2 {
            if !out_norm_data[i].at.is_empty() {
                at[i] = Some(vec![super::ichister::SpAtom::default(); num_inp_at]);
            }
        }
        if (out_norm_data[TAUT_NON].at.is_empty() && out_norm_data[TAUT_YES].at.is_empty())
            || inp_at.is_empty()
        {
            ret = -1;
            break 'body;
        }
        let out_idx = if !out_norm_data[TAUT_YES].at.is_empty() {
            TAUT_YES
        } else {
            TAUT_NON
        };
        out_norm_data[out_idx].at[..num_inp_at].clone_from_slice(&inp_at[..num_inp_at]);
        tgi.b_ignore_isotopic = 0;
        tgi.b_taut_flags = *pb_taut_flags;
        tgi.b_taut_flags_done = *pb_taut_flags_done;
        tgi.t_group = None;

        let num_removed_h;
        if TG_FLAG_H_ALREADY_REMOVED & b_taut_flags != 0 {
            let nd = if !out_norm_data[TAUT_YES].at.is_empty() {
                &out_norm_data[TAUT_YES]
            } else {
                &out_norm_data[TAUT_NON]
            };
            num_atoms = (nd.num_at - nd.num_removed_h) as usize;
            num_at_tg = num_atoms;
            num_removed_h = nd.num_removed_h as usize;
            tgi.tni.n_num_removed_explicit_h = num_removed_h as i32;
        } else {
            let out_at = &mut out_norm_data[out_idx].at;
            num_atoms = remove_terminal_hdt(&mut out_at[..num_inp_at], b_fix_term_h_chrg) as usize;
            num_at_tg = num_atoms;
            num_removed_h = num_inp_at - num_atoms;
            tgi.tni.n_num_removed_explicit_h = num_removed_h as i32;
            add_dt_to_num_h(&mut out_at[..num_atoms]);
        }
        super::ichinorm::mark_ring_systems_inp(&mut out_norm_data[out_idx].at[..num_atoms], 0);
        let out_copy = out_norm_data[out_idx].at[..num_inp_at].to_vec();
        if out_idx != TAUT_YES && !out_norm_data[TAUT_YES].at.is_empty() {
            out_norm_data[TAUT_YES].at[..num_inp_at].clone_from_slice(&out_copy);
        }
        if !out_norm_data[TAUT_YES].at.is_empty() {
            if let Some(f) = out_norm_data[TAUT_YES].at_fixed_bonds.as_mut() {
                f[..num_inp_at].clone_from_slice(&out_copy);
            }
        }
        if out_idx != TAUT_NON && !out_norm_data[TAUT_NON].at.is_empty() {
            out_norm_data[TAUT_NON].at[..num_inp_at].clone_from_slice(&out_copy);
        }

        // tautomeric structure
        if !out_norm_data[TAUT_YES].at.is_empty() && at[TAUT_YES].is_some() {
            {
                let nd = &mut out_norm_data[TAUT_YES];
                let mut dummy_done = 0u64;
                let r = mark_alt_bonds_and_taut_groups(
                    &mut nd.at,
                    nd.at_fixed_bonds.as_deref_mut(),
                    num_atoms,
                    Some(&mut tgi),
                    0,
                    &mut dummy_done,
                    &[],
                );
                ret = r;
            }
            if ret < 0 {
                break 'body;
            }
            num_taut_at = ret as usize;
            num_removed_h_taut = tgi.tni.n_num_removed_explicit_h as usize;
            {
                let nd = &mut out_norm_data[TAUT_YES];
                nd.num_at = (num_atoms + num_removed_h_taut) as i32;
                nd.num_removed_h = num_removed_h_taut as i32;
                nd.n_num_removed_protons += tgi.tni.n_num_removed_protons;
                for i in 0..3 {
                    nd.n_num_removed_protons_isotopic[i] +=
                        tgi.tni.n_num_removed_protons_isotopic[i] as i16;
                    nd.num_iso_h[i] += tgi.num_iso_h[i];
                }
                if num_taut_at == 1
                    && nd.at[0].at_type == ATT_PROTON
                    && tgi.tni.n_num_removed_protons == 1
                {
                    nd.b_deleted = 1;
                    nd.at_fixed_bonds = None;
                } else if tgi.tni.b_normalization_flags & FLAG_NORM_CONSIDER_TAUT != 0
                    && nd.at_fixed_bonds.is_some()
                {
                    nd.b_taut_preprocessed = 1;
                }
                *pb_taut_flags = tgi.b_taut_flags;
                nd.b_taut_flags = tgi.b_taut_flags;
                *pb_taut_flags_done = tgi.b_taut_flags_done;
                nd.b_taut_flags_done = tgi.b_taut_flags_done;
                nd.b_normalization_flags = tgi.tni.b_normalization_flags;
            }
            let sp = at[TAUT_YES].as_mut().unwrap();
            inp2sp_atom(&out_norm_data[TAUT_YES].at[..num_inp_at], sp);
            ret = set_stereo_parity(
                &mut out_norm_data[TAUT_YES].at,
                sp,
                num_taut_at,
                num_removed_h_taut,
                Some(&mut s[TAUT_YES].n_max_num_stereo_atoms),
                Some(&mut s[TAUT_YES].n_max_num_stereo_bonds),
                n_user_mode,
                b_pointed_edge_stereo,
                v_ab_parity_unknown,
                ip.b_loose_tsa_check != 0,
                ip.b_stereo_at_zz,
            );
            if returned_error(ret) {
                break 'body;
            }
            s[TAUT_YES].b_may_have_stereo = i32::from(
                s[TAUT_YES].n_max_num_stereo_atoms != 0 || s[TAUT_YES].n_max_num_stereo_bonds != 0,
            );
            s[TAUT_YES].num_isotopic_atoms = set_atom_iso_sort_keys(
                &mut sp[..num_taut_at],
                Some(&tgi),
                Some(&mut s[TAUT_YES].b_has_isotopic_taut_groups),
            );
            let mut eps: Vec<u16> = sp[..num_taut_at].iter().map(|a| a.endpoint).collect();
            let lt = count_tautomer_groups(&mut eps, num_taut_at, &mut tgi);
            for (a, e) in sp.iter_mut().zip(eps) {
                a.endpoint = e;
            }
            s[TAUT_YES].n_len_linear_ct_tautomer = lt;
            if returned_error(lt) {
                ret = lt;
                break 'body;
            } else if lt > 0 {
                num_at_tg = num_taut_at + tgi.num_t_groups;
                let l = set_tautomer_iso_sort_keys(&mut tgi);
                s[TAUT_YES].n_len_linear_ct_isotopic_tautomer = l.max(0);
                out_norm_data[TAUT_YES].b_tautomeric = lt;
            }
            get_canon_lengths(num_taut_at, sp, &mut s[TAUT_YES], Some(&tgi));
        }

        if !out_norm_data[TAUT_NON].at.is_empty()
            && !out_norm_data[TAUT_YES].at.is_empty()
            && at[TAUT_NON].is_some()
            && s[TAUT_YES].n_len_linear_ct_tautomer == 0
        {
            at[TAUT_NON] = None;
        } else if out_norm_data[TAUT_NON].at.is_empty()
            && !out_norm_data[TAUT_YES].at.is_empty()
            && at[TAUT_NON].is_none()
            && at[TAUT_YES].is_some()
            && s[TAUT_YES].n_len_linear_ct_tautomer == 0
        {
            out_norm_data[TAUT_YES].b_tautomeric = 0;
        } else if !out_norm_data[TAUT_NON].at.is_empty() && at[TAUT_NON].is_some() {
            {
                let nd = &mut out_norm_data[TAUT_NON];
                ret = mark_alt_bonds_and_taut_groups(
                    &mut nd.at,
                    None,
                    num_atoms,
                    None,
                    b_taut_flags,
                    &mut b_taut_flags_done,
                    &[],
                );
            }
            if ret < 0 {
                break 'body;
            }
            {
                let nd = &mut out_norm_data[TAUT_NON];
                nd.num_at = (num_atoms + num_removed_h) as i32;
                nd.num_removed_h = num_removed_h as i32;
                nd.b_taut_flags = *pb_taut_flags;
                nd.b_taut_flags_done = *pb_taut_flags_done;
                nd.b_normalization_flags = 0;
            }
            let sp = at[TAUT_NON].as_mut().unwrap();
            inp2sp_atom(&out_norm_data[TAUT_NON].at[..num_inp_at], sp);
            ret = set_stereo_parity(
                &mut out_norm_data[TAUT_NON].at,
                sp,
                num_atoms,
                num_removed_h,
                Some(&mut s[TAUT_NON].n_max_num_stereo_atoms),
                Some(&mut s[TAUT_NON].n_max_num_stereo_bonds),
                n_user_mode,
                b_pointed_edge_stereo,
                v_ab_parity_unknown,
                ip.b_loose_tsa_check != 0,
                ip.b_stereo_at_zz,
            );
            if returned_error(ret) {
                break 'body;
            }
            s[TAUT_NON].b_may_have_stereo = i32::from(
                s[TAUT_NON].n_max_num_stereo_atoms != 0 || s[TAUT_NON].n_max_num_stereo_bonds != 0,
            );
            s[TAUT_NON].num_isotopic_atoms =
                set_atom_iso_sort_keys(&mut sp[..num_atoms], None, None);
            get_canon_lengths(num_atoms, sp, &mut s[TAUT_NON], None);
            out_norm_data[TAUT_NON].b_tautomeric = 0;
        }

        // common
        let b_may_have_stereo =
            s[TAUT_YES].b_may_have_stereo != 0 || s[TAUT_NON].b_may_have_stereo != 0;
        b_has_isotopic_atoms = s[TAUT_NON].num_isotopic_atoms > 0
            || s[TAUT_NON].b_has_isotopic_taut_groups > 0
            || s[TAUT_YES].num_isotopic_atoms > 0
            || s[TAUT_YES].b_has_isotopic_taut_groups > 0;
        if b_fix_iso_fixed_h {
            b_has_isotopic_atoms = b_has_isotopic_atoms
                || (s[TAUT_YES].n_len_linear_ct_tautomer > 0
                    && tgi
                        .tni
                        .n_num_removed_protons_isotopic
                        .iter()
                        .any(|&x| x != 0));
        }
        b_has_isotopic_atoms = b_has_isotopic_atoms
            || (s[TAUT_YES].n_len_isotopic_endpoints > 1
                && tgi.b_taut_flags_done
                    & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                    != 0);
        if n_user_mode & REQ_MODE_DEFAULT == 0 {
            n_user_mode |= REQ_MODE_DEFAULT;
        }
        if n_user_mode & REQ_MODE_ISO != 0 && !b_has_isotopic_atoms {
            n_user_mode ^= REQ_MODE_ISO;
            n_user_mode |= REQ_MODE_NON_ISO;
        }
        if n_user_mode & REQ_MODE_STEREO != 0 && n_user_mode & REQ_MODE_ISO != 0 {
            n_user_mode |= REQ_MODE_ISO_STEREO;
        }
        if n_user_mode & REQ_MODE_STEREO != 0 && n_user_mode & REQ_MODE_NON_ISO == 0 {
            n_user_mode ^= REQ_MODE_STEREO;
        }
        if !b_may_have_stereo {
            if n_user_mode & REQ_MODE_STEREO != 0 {
                n_user_mode ^= REQ_MODE_STEREO;
            }
            if n_user_mode & REQ_MODE_ISO_STEREO != 0 {
                n_user_mode ^= REQ_MODE_ISO_STEREO;
            }
        }
        if n_user_mode & REQ_MODE_BASIC != 0
            && (out_norm_data[TAUT_NON].at.is_empty()
                || p_inchi[TAUT_NON].is_none()
                || p_aux[TAUT_NON].is_none()
                || at[TAUT_NON].is_none())
        {
            n_user_mode ^= REQ_MODE_BASIC;
        }
        if n_user_mode & REQ_MODE_TAUT != 0
            && (out_norm_data[TAUT_YES].at.is_empty()
                || p_inchi[TAUT_YES].is_none()
                || p_aux[TAUT_YES].is_none()
                || at[TAUT_YES].is_none())
        {
            n_user_mode ^= REQ_MODE_TAUT;
        }
        let (mut n1, mut n2) = match n_user_mode & (REQ_MODE_BASIC | REQ_MODE_TAUT) {
            REQ_MODE_BASIC => (TAUT_NON, TAUT_NON),
            REQ_MODE_TAUT => (TAUT_YES, TAUT_YES),
            x if x == REQ_MODE_BASIC | REQ_MODE_TAUT => (TAUT_NON, TAUT_YES),
            _ => {
                ret = -3;
                break 'body;
            }
        };

        if n_user_mode & REQ_MODE_NON_ISO != 0 && n_user_mode & REQ_MODE_ISO == 0 {
            {
                tgi.b_ignore_isotopic = 1;
                if !tgi.n_isotopic_endpoint_atom_number.is_empty() {
                    tgi.n_isotopic_endpoint_atom_number[0] =
                        tgi.n_isotopic_endpoint_atom_number[0].min(1);
                }
                tgi.num_iso_h = [0; 3];
                tgi.tni.n_num_removed_protons_isotopic = [0; 3];
                tgi.b_taut_flags_done &=
                    !(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
            }
            for x in s.iter_mut() {
                x.b_has_isotopic_taut_groups = 0;
                x.b_ignore_isotopic = 1;
                x.n_len_isotopic = 0;
                x.n_len_isotopic_endpoints = 0;
                x.n_len_linear_ct_isotopic_tautomer = 0;
                x.num_isotopic_atoms = 0;
            }
            b_has_isotopic_atoms = false;
        }

        let at_refs: [Option<&[super::ichister::SpAtom]>; 2] =
            [at[TAUT_NON].as_deref(), at[TAUT_YES].as_deref()];
        ret = get_base_canon_ranking(
            &mut || false,
            num_atoms,
            num_at_tg,
            at_refs,
            Some(&mut tgi),
            &s,
            &mut bcn,
            b_fix_iso_fixed_h,
            large_molecules,
        );
        if ret < 0 {
            break 'body;
        }
        if bcn.ftcn[n1].partition_ct.rank.is_empty() {
            n1 = alt_taut(n1);
        }
        if bcn.ftcn[n2].partition_ct.rank.is_empty() {
            n2 = alt_taut(n2);
        }
        if n1 > n2 {
            ret = CT_TAUCOUNT_ERR;
            break 'body;
        }

        let mut i = n2 as i32;
        while i >= n1 as i32 && !returned_error(ret) {
            let iu = i as usize;
            let mut cs = CanonStat::default();
            let mut n_mode: i64 = if s[iu].n_len_linear_ct_tautomer == 0 {
                CMODE_CT
            } else {
                CMODE_CT | CMODE_TAUT
            };
            if b_has_isotopic_atoms && n_user_mode & REQ_MODE_ISO != 0 {
                n_mode |= CANON_MODE_ISO;
            }
            if s[iu].b_may_have_stereo != 0 && n_user_mode & REQ_MODE_STEREO != 0 {
                n_mode |= CANON_MODE_STEREO;
            }
            if b_has_isotopic_atoms
                && s[iu].b_may_have_stereo != 0
                && n_user_mode & REQ_MODE_ISO_STEREO != 0
            {
                n_mode |= CANON_MODE_ISO_STEREO;
            }
            if n_user_mode & REQ_MODE_NOEQ_STEREO != 0 {
                n_mode |= CMODE_NOEQ_STEREO;
            }
            if n_user_mode & REQ_MODE_REDNDNT_STEREO != 0 {
                n_mode |= CMODE_REDNDNT;
            }
            if n_user_mode & REQ_MODE_NO_ALT_SBONDS != 0 {
                n_mode |= CMODE_NO_ALT;
            }
            if n_mode & CANON_MODE_STEREO == CANON_MODE_STEREO
                || n_mode & CANON_MODE_ISO_STEREO == CANON_MODE_ISO_STEREO
            {
                if n_user_mode & REQ_MODE_RELATIVE_STEREO != 0 {
                    n_mode |= CMODE_REL;
                }
                if n_user_mode & REQ_MODE_RACEMIC_STEREO != 0 {
                    n_mode |= CMODE_RAC;
                }
                if n_user_mode & REQ_MODE_SC_IGN_ALL_UU != 0 {
                    n_mode |= CMODE_SC_UU;
                }
                if n_user_mode & REQ_MODE_SB_IGN_ALL_UU != 0 {
                    n_mode |= CMODE_SB_UU;
                }
            }
            let (na, nat) = if iu == TAUT_NON {
                (num_atoms, num_atoms)
            } else {
                (num_atoms, num_at_tg)
            };
            ret = if iu == TAUT_NON {
                allocate_cs(
                    &mut cs,
                    na,
                    nat,
                    s[iu].n_len_ct,
                    s[iu].n_len_ct_at_only,
                    s[iu].n_len_linear_ct_stereo_dble,
                    s[iu].n_max_num_stereo_bonds,
                    s[iu].n_len_linear_ct_stereo_carb,
                    s[iu].n_max_num_stereo_atoms,
                    0,
                    0,
                    s[iu].n_len_isotopic,
                    n_mode,
                    true,
                )
            } else {
                allocate_cs(
                    &mut cs,
                    na,
                    nat,
                    s[iu].n_len_ct,
                    s[iu].n_len_ct_at_only,
                    s[iu].n_len_linear_ct_stereo_dble,
                    s[iu].n_max_num_stereo_bonds,
                    s[iu].n_len_linear_ct_stereo_carb,
                    s[iu].n_max_num_stereo_atoms,
                    s[iu].n_len_linear_ct_tautomer,
                    s[iu].n_len_linear_ct_isotopic_tautomer,
                    s[iu].n_len_isotopic,
                    n_mode,
                    true,
                )
            };
            if ret != 0 {
                break 'body;
            }
            if n_user_mode & REQ_MODE_DIFF_UU_STEREO != 0 {
                n_mode |= REQ_MODE_DIFF_UU_STEREO as i64;
            }
            cs.l_num_decreased_ct = -1;
            cs.b_double_bond_square = 0;
            cs.b_ignore_isotopic = i32::from(
                !((s[TAUT_NON].num_isotopic_atoms != 0
                    || s[TAUT_YES].num_isotopic_atoms != 0
                    || s[TAUT_YES].b_has_isotopic_taut_groups != 0)
                    || n_user_mode & REQ_MODE_NON_ISO != 0
                    || n_user_mode & REQ_MODE_ISO == 0),
            );
            if n_user_mode & REQ_MODE_NON_ISO != 0 && n_user_mode & REQ_MODE_ISO == 0 {
                cs.b_ignore_isotopic = 1;
            }
            let have_tgi_cs = iu == TAUT_YES;
            if have_tgi_cs {
                tgi.b_ignore_isotopic = i32::from(
                    !(s[TAUT_YES].b_has_isotopic_taut_groups != 0
                        || n_user_mode & REQ_MODE_NON_ISO != 0
                        || n_user_mode & REQ_MODE_ISO == 0),
                );
                if n_user_mode & REQ_MODE_NON_ISO != 0 && n_user_mode & REQ_MODE_ISO == 0 {
                    tgi.b_ignore_isotopic = 1;
                }
                cs.t_group_info = Some(std::mem::take(&mut tgi));
            }
            cs.neigh_list = None;
            cs.p_bcn = Some(Box::new(std::mem::take(&mut bcn)));
            let sp = at[iu].as_mut().unwrap();
            ret = canon_inchi(&mut || false, num_atoms, nat, sp, &mut cs, n_mode, iu);
            bcn = *cs.p_bcn.take().unwrap();
            if have_tgi_cs {
                tgi = cs.t_group_info.take().unwrap();
            }
            let p = p_inchi[iu].as_mut().unwrap();
            let a = p_aux[iu].as_mut().unwrap();
            if ret <= 0 {
                p.n_error_code = ret;
                a.n_error_code = ret;
            } else {
                p.n_error_code = 0;
                a.n_error_code = 0;
                p.b_deleted = out_norm_data[iu].b_deleted;
                a.b_deleted = out_norm_data[iu].b_deleted;
                a.n_canon_flags = cs.n_canon_flags;
                a.b_taut_flags = out_norm_data[iu].b_taut_flags;
                a.b_taut_flags_done = out_norm_data[iu].b_taut_flags_done;
                a.b_normalization_flags = out_norm_data[iu].b_normalization_flags;
                let removed = if iu != 0 {
                    num_removed_h_taut
                } else {
                    num_removed_h
                };
                ret = fill_out_inchi(
                    p,
                    a,
                    num_atoms,
                    nat,
                    removed,
                    sp,
                    &mut out_norm_data[iu].at,
                    &mut cs,
                    have_tgi_cs.then_some(&mut tgi),
                    iu != 0,
                    n_user_mode,
                    str_err,
                    ip.b_no_warnings != 0,
                );
                if returned_error(ret) {
                    p.n_error_code = ret;
                    a.n_error_code = ret;
                } else {
                    if s[TAUT_YES].n_len_linear_ct_tautomer != 0 {
                        a.b_is_tautomeric = s[TAUT_YES].n_len_linear_ct_tautomer;
                    }
                    cs.t_group_info = have_tgi_cs.then(|| std::mem::take(&mut tgi));
                    let ret2 = check_canon_numbering_correctness(num_atoms, nat, sp, &mut cs);
                    if let Some(t) = cs.t_group_info.take() {
                        tgi = t;
                    }
                    if ret2 != 0 {
                        p.n_error_code = ret2;
                        a.n_error_code = ret2;
                        ret = ret2;
                    }
                }
            }
            i -= 1;
        }
        if ret == 0 {
            ret = num_atoms as i32;
        }
    }
    if let Some(t) = ti_out {
        *t = tgi;
    }
    ret
}

// ---------------------------------------------------------------------------
// Comparison of the InChI of a reversed structure with the input InChI

pub const IDIF_PROBLEM: i32 = 0x0000_0001;
pub const IDIF_NUM_AT: i32 = 0x0000_0001;
pub const IDIF_ATOMS: i32 = 0x0000_0001;
pub const IDIF_NUM_EL: i32 = 0x0000_0001;
pub const IDIF_CON_LEN: i32 = 0x0000_0001;
pub const IDIF_CON_TBL: i32 = 0x0000_0001;
pub const IDIF_POSITION_H: i32 = 0x0000_0002;
pub const IDIF_MORE_FH: i32 = 0x0000_0004;
pub const IDIF_LESS_FH: i32 = 0x0000_0008;
pub const IDIF_MORE_H: i32 = 0x0000_0010;
pub const IDIF_LESS_H: i32 = 0x0000_0020;
pub const IDIF_NO_TAUT: i32 = 0x0000_0040;
pub const IDIF_WRONG_TAUT: i32 = 0x0000_0080;
pub const IDIF_SINGLE_TG: i32 = 0x0000_0100;
pub const IDIF_MULTIPLE_TG: i32 = 0x0000_0200;
pub const IDIF_NUM_TG: i32 = 0x0000_0400;
pub const IDIF_EXTRA_TG_ENDP: i32 = 0x0000_0800;
pub const IDIF_MISS_TG_ENDP: i32 = 0x0000_1000;
pub const IDIF_DIFF_TG_ENDP: i32 = 0x0000_2000;
pub const IDIF_TG: i32 = 0x0000_4000;
pub const IDIF_NUM_ISO_AT: i32 = 0x0000_8000;
pub const IDIF_ISO_AT: i32 = 0x0001_0000;
pub const IDIF_CHARGE: i32 = 0x0002_0000;
pub const IDIF_REM_PROT: i32 = 0x0004_0000;
pub const IDIF_REM_ISO_H: i32 = 0x0008_0000;
pub const ICR_MAX_ENDP_IN1_ONLY: usize = 32;
pub const ICR_MAX_ENDP_IN2_ONLY: usize = 32;
pub const ICR_MAX_DIFF_FIXED_H: usize = 32;
pub const IDIF_SC_INV: i32 = 0x0010_0000;
pub const IDIF_SC_PARITY: i32 = 0x0020_0000;
pub const IDIF_SC_EXTRA_UNDF: i32 = 0x0040_0000;
pub const IDIF_SC_EXTRA: i32 = 0x0080_0000;
pub const IDIF_SC_MISS_UNDF: i32 = 0x0100_0000;
pub const IDIF_SC_MISS: i32 = 0x0200_0000;
pub const IDIF_SB_PARITY: i32 = 0x0400_0000;
pub const IDIF_SB_EXTRA_UNDF: i32 = 0x0800_0000;
pub const IDIF_SB_EXTRA: i32 = 0x1000_0000;
pub const IDIF_SB_MISS_UNDF: i32 = 0x2000_0000;
pub const IDIF_SB_MISS: i32 = 0x4000_0000;

pub const ICR_MAX_SB_IN1_ONLY: usize = 32;
pub const ICR_MAX_SB_IN2_ONLY: usize = 32;
pub const ICR_MAX_SC_IN1_ONLY: usize = 32;
pub const ICR_MAX_SC_IN2_ONLY: usize = 32;
pub const ICR_MAX_SB_UNDF: usize = 32;
pub const ICR_MAX_SC_UNDF: usize = 32;

/// The stereo part of `ICR` (`tagInChICompareResults`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Icr {
    pub flags: u64,
    pub tot_num_h1: i32,
    pub tot_num_h2: i32,
    pub num_taut_h1: i32,
    pub num_taut_h2: i32,
    pub num_taut_m1: i32,
    pub num_taut_m2: i32,
    pub endp_in1_only: Vec<u16>,
    pub endp_in2_only: Vec<u16>,
    /// (atom index, difference in the number of H)
    pub diff_pos_h: Vec<(u16, i8)>,
    pub fixed_h1_more: Vec<(u16, i8)>,
    pub fixed_h2_more: Vec<(u16, i8)>,
    pub sc_in1_only: Vec<u16>,
    pub sc_in2_only: Vec<u16>,
    pub sb_in1_only: Vec<u16>,
    pub sb_in2_only: Vec<u16>,
    pub sb_undef_in1_only: Vec<u16>,
    pub sb_undef_in2_only: Vec<u16>,
    pub sc_undef_in1_only: Vec<u16>,
    pub sc_undef_in2_only: Vec<u16>,
}

/// `CompareReversedStereoINChI`
pub fn compare_reversed_stereo_inchi(s1: Option<&InChIStereo>, s2: Option<&InChIStereo>) -> i32 {
    let (s1, s2) = match (s1, s2) {
        (None, None) => return 0,
        (Some(s), None) | (None, Some(s)) => {
            return if s.n_number_of_stereo_centers != 0 || s.n_number_of_stereo_bonds != 0 {
                20
            } else {
                0
            };
        }
        (Some(a), Some(b)) => (a, b),
    };
    if s1.n_number_of_stereo_centers != s2.n_number_of_stereo_centers {
        return 21;
    }
    let nc = s1.n_number_of_stereo_centers.max(0) as usize;
    if nc > 0 {
        if s1.n_number[..nc] != s2.n_number[..nc] {
            return 22;
        }
        if s1.t_parity[..nc] != s2.t_parity[..nc] {
            return 23;
        }
        if s1.n_comp_inv2_abs != s2.n_comp_inv2_abs
            && s1.n_comp_inv2_abs != 0
            && s2.n_comp_inv2_abs != 0
        {
            return 24;
        }
    }
    if s1.n_number_of_stereo_bonds != s2.n_number_of_stereo_bonds {
        return 25;
    }
    let nb = s1.n_number_of_stereo_bonds.max(0) as usize;
    if nb > 0 {
        if s1.n_bond_atom1[..nb] != s2.n_bond_atom1[..nb] {
            return 26;
        }
        if s1.n_bond_atom2[..nb] != s2.n_bond_atom2[..nb] {
            return 27;
        }
        if s1.b_parity[..nb] != s2.b_parity[..nb] {
            return 28;
        }
    }
    0
}

fn push_max(v: &mut Vec<u16>, max: usize, x: usize) {
    if v.len() < max {
        v.push(x as u16);
    }
}

/// `CompareReversedStereoINChI2`
pub fn compare_reversed_stereo_inchi2(
    s1: Option<&InChIStereo>,
    s2: Option<&InChIStereo>,
    icr: &mut Icr,
) -> i32 {
    use super::ichister::AB_PARITY_UNDF;
    let undf = AB_PARITY_UNDF as i8;
    let mut ret = 0;
    let add_sb = icr.sb_undef_in1_only.len() + icr.sb_in1_only.len() + icr.sb_in2_only.len() == 0;
    let add_sc = icr.sc_undef_in1_only.len() + icr.sc_in1_only.len() + icr.sc_in2_only.len() == 0;
    let nsc1 = s1.map_or(0, |s| s.n_number_of_stereo_centers.max(0) as usize);
    let nsc2 = s2.map_or(0, |s| s.n_number_of_stereo_centers.max(0) as usize);
    let nsb1 = s1.map_or(0, |s| s.n_number_of_stereo_bonds.max(0) as usize);
    let nsb2 = s2.map_or(0, |s| s.n_number_of_stereo_bonds.max(0) as usize);

    if let (Some(s1), Some(s2)) = (s1, s2) {
        if (nsc1 != 0 || nsc2 != 0)
            && (nsc1 != nsc2
                || s1.n_number[..nsc1] != s2.n_number[..nsc1]
                || s1.t_parity[..nsc1] != s2.t_parity[..nsc1])
        {
            let (mut num_dif, mut extra_undf, mut miss_undf, mut in1_only, mut in2_only) =
                (0, 0, 0, 0, 0);
            let (mut j1, mut j2) = (0, 0);
            while j1 < nsc1 && j2 < nsc2 {
                if s1.n_number[j1] == s2.n_number[j2] {
                    if s1.t_parity[j1] != s2.t_parity[j2] {
                        num_dif += 1;
                    }
                    j1 += 1;
                    j2 += 1;
                } else if s1.n_number[j1] < s2.n_number[j2] {
                    in1_only += 1;
                    if s1.t_parity[j1] == undf {
                        extra_undf += 1;
                    }
                    if add_sc {
                        push_max(&mut icr.sc_in1_only, ICR_MAX_SC_IN1_ONLY, j1);
                        if s1.t_parity[j1] == undf {
                            push_max(&mut icr.sc_undef_in1_only, ICR_MAX_SC_UNDF, j1);
                        }
                    }
                    j1 += 1;
                } else {
                    in2_only += 1;
                    if s2.t_parity[j2] == undf {
                        miss_undf += 1;
                    }
                    if add_sc {
                        push_max(&mut icr.sc_in2_only, ICR_MAX_SC_IN2_ONLY, j2);
                        if s2.t_parity[j2] == undf {
                            // as in the C code: records `j1`
                            push_max(&mut icr.sc_undef_in2_only, ICR_MAX_SC_UNDF, j1);
                        }
                    }
                    j2 += 1;
                }
            }
            while j1 < nsc1 {
                if s1.t_parity[j1] == undf {
                    extra_undf += 1;
                }
                in1_only += 1;
                if add_sc {
                    push_max(&mut icr.sc_in1_only, ICR_MAX_SC_IN1_ONLY, j1);
                    if s1.t_parity[j1] == undf {
                        push_max(&mut icr.sc_undef_in1_only, ICR_MAX_SC_UNDF, j1);
                    }
                }
                j1 += 1;
            }
            while j2 < nsc2 {
                if s2.t_parity[j2] == undf {
                    miss_undf += 1;
                }
                in2_only += 1;
                if add_sc {
                    push_max(&mut icr.sc_in2_only, ICR_MAX_SC_IN2_ONLY, j2);
                }
                j2 += 1;
            }
            if num_dif != 0 {
                ret |= IDIF_SC_PARITY;
            }
            if in1_only != 0 {
                if extra_undf != 0 {
                    ret |= IDIF_SC_EXTRA_UNDF;
                }
                if in1_only != extra_undf {
                    ret |= IDIF_SC_EXTRA;
                }
            }
            if in2_only != 0 {
                if miss_undf != 0 {
                    ret |= IDIF_SC_MISS_UNDF;
                }
                if in2_only != miss_undf {
                    ret |= IDIF_SC_MISS;
                }
            }
        }
        if s1.n_comp_inv2_abs != s2.n_comp_inv2_abs
            && s1.n_comp_inv2_abs != 0
            && s2.n_comp_inv2_abs != 0
        {
            ret |= IDIF_SC_INV;
        }
    }

    let differs = nsb1 != nsb2
        || match (s1, s2) {
            (Some(s1), Some(s2)) => {
                s1.n_bond_atom1[..nsb1] != s2.n_bond_atom1[..nsb1]
                    || s1.n_bond_atom2[..nsb1] != s2.n_bond_atom2[..nsb1]
                    || s1.b_parity[..nsb1] != s2.b_parity[..nsb1]
            }
            _ => false,
        };
    if (nsb1 != 0 || nsb2 != 0) && differs {
        let (mut num_dif, mut extra_undf, mut miss_undf, mut in1_only, mut in2_only) =
            (0, 0, 0, 0, 0);
        let (mut j1, mut j2) = (0, 0);
        let (a1, a2, p1, p2): (&[u16], &[u16], &[i8], &[i8]) = (
            s1.map_or(&[], |s| &s.n_bond_atom1),
            s2.map_or(&[], |s| &s.n_bond_atom1),
            s1.map_or(&[], |s| &s.b_parity),
            s2.map_or(&[], |s| &s.b_parity),
        );
        let (c1, c2): (&[u16], &[u16]) = (
            s1.map_or(&[], |s| &s.n_bond_atom2),
            s2.map_or(&[], |s| &s.n_bond_atom2),
        );
        while j1 < nsb1 && j2 < nsb2 {
            if a1[j1] == a2[j2] && c1[j1] == c2[j2] {
                if p1[j1] != p2[j2] {
                    num_dif += 1;
                }
                j1 += 1;
                j2 += 1;
            } else if a1[j1] < a2[j2] || (a1[j1] == a2[j2] && c1[j1] < c2[j2]) {
                in1_only += 1;
                if p1[j1] == undf {
                    extra_undf += 1;
                }
                if add_sb {
                    push_max(&mut icr.sb_in1_only, ICR_MAX_SB_IN1_ONLY, j1);
                    if p1[j1] == undf {
                        push_max(&mut icr.sb_undef_in1_only, ICR_MAX_SB_UNDF, j1);
                    }
                }
                j1 += 1;
            } else {
                in2_only += 1;
                if p2[j2] == undf {
                    miss_undf += 1;
                }
                if add_sb {
                    push_max(&mut icr.sb_in2_only, ICR_MAX_SB_IN2_ONLY, j2);
                    if p2[j2] == undf {
                        // as in the C code: records `j1`
                        push_max(&mut icr.sb_undef_in2_only, ICR_MAX_SB_UNDF, j1);
                    }
                }
                j2 += 1;
            }
        }
        while j1 < nsb1 {
            in1_only += 1;
            if p1[j1] == undf {
                extra_undf += 1;
            }
            if add_sb {
                push_max(&mut icr.sb_in1_only, ICR_MAX_SB_IN1_ONLY, j1);
                if p1[j1] == undf {
                    push_max(&mut icr.sb_undef_in1_only, ICR_MAX_SB_UNDF, j1);
                }
            }
            j1 += 1;
        }
        while j2 < nsb2 {
            in2_only += 1;
            if p2[j2] == undf {
                miss_undf += 1;
            }
            if add_sb {
                push_max(&mut icr.sb_in2_only, ICR_MAX_SB_IN2_ONLY, j2);
                if p2[j2] == undf {
                    push_max(&mut icr.sb_undef_in2_only, ICR_MAX_SB_UNDF, j1);
                }
            }
            j2 += 1;
        }
        if num_dif != 0 {
            ret |= IDIF_SB_PARITY;
        }
        if in1_only != 0 {
            if extra_undf != 0 {
                ret |= IDIF_SB_EXTRA_UNDF;
            }
            if in1_only != extra_undf {
                ret |= IDIF_SB_EXTRA;
            }
        }
        if in2_only != 0 {
            if miss_undf != 0 {
                ret |= IDIF_SB_MISS_UNDF;
            }
            if in2_only != miss_undf {
                ret |= IDIF_SB_MISS;
            }
        }
    }
    ret
}

/// `CompareReversedINChI`
pub fn compare_reversed_inchi(
    i1: Option<&InChI>,
    i2: Option<&InChI>,
    a1: Option<&InChIAux>,
    a2: Option<&InChIAux>,
) -> i32 {
    let (i1, i2) = match (i1, i2) {
        (None, None) => return 0,
        (Some(_), None) | (None, Some(_)) => return 1,
        (Some(a), Some(b)) => (a, b),
    };
    if i1.n_error_code == i2.n_error_code {
        if i1.n_error_code != 0 {
            return 0;
        }
    } else {
        return 2;
    }
    if i1.b_deleted != i2.b_deleted {
        return 1;
    }
    if i1.n_number_of_atoms != i2.n_number_of_atoms {
        return 3;
    }
    let n = i1.n_number_of_atoms.max(0) as usize;
    if n > 0 {
        if i1.n_atom[..n] != i2.n_atom[..n] {
            return 4;
        }
        if i1.sz_hill_formula != i2.sz_hill_formula {
            return 7;
        }
        if i1.n_num_h[..n] != i2.n_num_h[..n] {
            return if i1.len_conn_table > 1 || i2.len_conn_table > 1 {
                5
            } else {
                6
            };
        }
        if i1.n_num_h_fixed.is_some() || i2.n_num_h_fixed.is_some() {
            let has =
                |f: &Option<Vec<i8>>| f.as_ref().is_some_and(|v| v[..n].iter().any(|&x| x != 0));
            let (h1, h2) = (has(&i1.n_num_h_fixed), has(&i2.n_num_h_fixed));
            let (mut j1, mut j2) = (0, 0);
            if h1 && !h2 {
                return 18;
            } else if !h1 && h2 {
                return 19;
            } else if h1 && h2 {
                let (f1, f2) = (
                    i1.n_num_h_fixed.as_ref().unwrap(),
                    i2.n_num_h_fixed.as_ref().unwrap(),
                );
                if f1[..n] != f2[..n] {
                    for i in 0..n {
                        if f1[i] > f2[i] {
                            j1 += 1;
                        } else if f1[i] < f2[i] {
                            j2 += 1;
                        }
                    }
                }
            }
            let ret = if j1 != 0 && j2 != 0 {
                20
            } else if j1 != 0 {
                18
            } else if j2 != 0 {
                19
            } else {
                0
            };
            if ret != 0 {
                return ret;
            }
        }
    }
    if i1.len_conn_table != i2.len_conn_table {
        return 8;
    }
    let lc = i1.len_conn_table.max(0) as usize;
    if lc > 0 && i1.n_conn_table[..lc] != i2.n_conn_table[..lc] {
        return 9;
    }
    if i1.len_tautomer != i2.len_tautomer && (i1.len_tautomer > 1 || i2.len_tautomer > 1) {
        return 10;
    }
    let lt = i1.len_tautomer.max(0) as usize;
    if i1.len_tautomer > 1 && i2.len_tautomer > 1 && i1.n_tautomer[..lt] != i2.n_tautomer[..lt] {
        return 11;
    }
    if i1.n_number_of_isotopic_atoms != i2.n_number_of_isotopic_atoms {
        return 12;
    }
    let ni = i1.n_number_of_isotopic_atoms.max(0) as usize;
    if ni > 0 && i1.isotopic_atom[..ni] != i2.isotopic_atom[..ni] {
        return 13;
    }
    if i1.n_total_charge != i2.n_total_charge {
        return 14;
    }
    if let (Some(a1), Some(a2)) = (a1, a2) {
        if a1.n_num_removed_protons != a2.n_num_removed_protons {
            return 16;
        }
        if a1.n_num_removed_isotopic_h != a2.n_num_removed_isotopic_h {
            return 17;
        }
    }
    let ret = compare_reversed_stereo_inchi(i1.stereo.as_ref(), i2.stereo.as_ref());
    if ret != 0 {
        return ret + 20;
    }
    let dup = i2.stereo_isotopic.is_none()
        && i2.stereo.is_some()
        && i1
            .stereo_isotopic
            .as_ref()
            .is_some_and(|s| s.n_number_of_stereo_bonds + s.n_number_of_stereo_centers > 0)
        && compare_reversed_stereo_inchi(i1.stereo_isotopic.as_ref(), i2.stereo.as_ref()) == 0;
    if !dup {
        let ret =
            compare_reversed_stereo_inchi(i1.stereo_isotopic.as_ref(), i2.stereo_isotopic.as_ref());
        if ret != 0 {
            return ret + 40;
        }
    }
    0
}

/// `CompareIcr`: returns (result, bits only in `icr1`, bits only in `icr2`)
pub fn compare_icr(icr1: &Icr, icr2: &Icr, mask: u64) -> (i32, u64, u64) {
    let (mut flg1, mut flg2, mut cur_bit) = (icr1.flags, icr2.flags, 1u64);
    let (mut in1, mut in2) = (0u64, 0u64);
    let (mut extra1, mut extra2) = (0, 0);
    let mut i = 0;
    while flg1 != 0 || flg2 != 0 {
        if mask & cur_bit != 0 {
            let (b1, b2) = (flg1 & 1, flg2 & 1);
            if b1 != 0 && b2 == 0 {
                in1 |= 1 << i;
                extra1 += 1;
            } else if b1 == 0 && b2 != 0 {
                in2 |= 1 << i;
                extra2 += 1;
            }
        }
        i += 1;
        flg1 >>= 1;
        flg2 >>= 1;
        cur_bit <<= 1;
    }
    let ret = if extra1 != 0 && extra2 == 0 {
        1
    } else if extra1 == 0 && extra2 != 0 {
        -1
    } else if in1 == 0 && in2 == 0 {
        0
    } else {
        2
    };
    (ret, in1, in2)
}

/// `CompareReversedINChI2`; the error argument (allocation failure) is
/// always 0.
pub fn compare_reversed_inchi2(
    i1: Option<&InChI>,
    i2: Option<&InChI>,
    a1: Option<&InChIAux>,
    a2: Option<&InChIAux>,
    icr: &mut Icr,
) -> u64 {
    *icr = Icr::default();
    let (i1, i2) = match (i1, i2) {
        (None, None) => return 0,
        (Some(a), Some(b)) => (a, b),
        _ => {
            icr.flags = IDIF_PROBLEM as u64;
            return icr.flags;
        }
    };
    let ret = compare_reversed_inchi2_body(i1, i2, a1, a2, icr);
    icr.flags = ret as u64;
    ret as u64
}

fn compare_reversed_inchi2_body(
    i1: &InChI,
    i2: &InChI,
    a1: Option<&InChIAux>,
    a2: Option<&InChIAux>,
    icr: &mut Icr,
) -> i32 {
    let mut ret = 0;
    if i1.n_error_code != i2.n_error_code || i1.n_error_code != 0 {
        return ret | IDIF_PROBLEM;
    }
    if i1.n_number_of_atoms != i2.n_number_of_atoms {
        return ret | IDIF_NUM_AT;
    }
    let n = i1.n_number_of_atoms.max(0) as usize;
    if n > 0 {
        if i1.n_atom[..n] != i2.n_atom[..n] {
            return ret | IDIF_ATOMS;
        }
        if i1.n_num_h[..n] != i2.n_num_h[..n] {
            ret |= IDIF_POSITION_H;
            for j in 0..n {
                if i1.n_num_h[j] != i2.n_num_h[j] && icr.diff_pos_h.len() < ICR_MAX_DIFF_FIXED_H {
                    icr.diff_pos_h
                        .push((j as u16, i1.n_num_h[j].wrapping_sub(i2.n_num_h[j])));
                }
            }
        }
        if i1.n_num_h_fixed.is_some() || i2.n_num_h_fixed.is_some() {
            let has =
                |f: &Option<Vec<i8>>| f.as_ref().is_some_and(|v| v[..n].iter().any(|&x| x != 0));
            let (h1, h2) = (has(&i1.n_num_h_fixed), has(&i2.n_num_h_fixed));
            if h1 && !h2 {
                let f1 = i1.n_num_h_fixed.as_ref().unwrap();
                for (i, &x) in f1[..n].iter().enumerate() {
                    if x != 0 && icr.fixed_h1_more.len() < ICR_MAX_DIFF_FIXED_H {
                        icr.fixed_h1_more.push((i as u16, x));
                    }
                }
                ret |= IDIF_MORE_FH;
            } else if !h1 && h2 {
                let f2 = i2.n_num_h_fixed.as_ref().unwrap();
                for (i, &x) in f2[..n].iter().enumerate() {
                    if x != 0 && icr.fixed_h2_more.len() < ICR_MAX_DIFF_FIXED_H {
                        icr.fixed_h2_more.push((i as u16, x));
                    }
                }
                ret |= IDIF_LESS_FH;
            } else if h1 && h2 {
                let (f1, f2) = (
                    i1.n_num_h_fixed.as_ref().unwrap(),
                    i2.n_num_h_fixed.as_ref().unwrap(),
                );
                if f1[..n] != f2[..n] {
                    for i in 0..n {
                        if f1[i] > f2[i] {
                            if icr.fixed_h1_more.len() < ICR_MAX_DIFF_FIXED_H {
                                icr.fixed_h1_more
                                    .push((i as u16, f1[i].wrapping_sub(f2[i])));
                            }
                        } else if f1[i] < f2[i] && icr.fixed_h2_more.len() < ICR_MAX_DIFF_FIXED_H {
                            icr.fixed_h2_more
                                .push((i as u16, f2[i].wrapping_sub(f1[i])));
                        }
                    }
                    if !icr.fixed_h1_more.is_empty() {
                        ret |= IDIF_MORE_FH;
                    }
                    if !icr.fixed_h2_more.is_empty() {
                        ret |= IDIF_LESS_FH;
                    }
                }
            }
        }
    }

    let (mut num_h1, mut num_h2) = (0, 0);
    let ret2 = compare_hill_formulas_no_h(
        i1.sz_hill_formula.as_deref().unwrap_or(""),
        i2.sz_hill_formula.as_deref().unwrap_or(""),
        &mut num_h1,
        &mut num_h2,
    );
    icr.tot_num_h1 = num_h1;
    icr.tot_num_h2 = num_h2;
    if ret2 != 0 {
        return ret | IDIF_NUM_EL;
    }
    if num_h1 > num_h2 {
        ret |= IDIF_MORE_H;
    }
    if num_h1 < num_h2 {
        ret |= IDIF_LESS_H;
    }
    if i1.len_conn_table != i2.len_conn_table {
        return ret | IDIF_CON_LEN;
    }
    let lc = i1.len_conn_table.max(0) as usize;
    if lc > 0 && i1.n_conn_table[..lc] != i2.n_conn_table[..lc] {
        return ret | IDIF_CON_TBL;
    }

    let n1 = if i1.len_tautomer != 0 {
        i1.n_tautomer[0] as i32
    } else {
        0
    };
    let n2 = if i2.len_tautomer != 0 {
        i2.n_tautomer[0] as i32
    } else {
        0
    };
    if n1 == 0 && n2 != 0 {
        ret |= IDIF_NO_TAUT;
    } else if n1 != 0 && n2 == 0 {
        ret |= IDIF_WRONG_TAUT;
    } else if n1 == 1 && n2 > 1 {
        ret |= IDIF_SINGLE_TG;
    } else if n1 > 1 && n2 == 1 {
        ret |= IDIF_MULTIPLE_TG;
    } else if n1 != n2 {
        ret |= IDIF_NUM_TG;
    }

    if n1 != 0 || n2 != 0 {
        let collect = |t: &[u16], len: i32| {
            let g = |k: usize| i32::from(t.get(k).copied().unwrap_or(0));
            let (mut nh, mut nm, mut pe) = (0, 0, Vec::new());
            let mut m = 1usize;
            while (m as i32) < len {
                let l = g(m) as usize;
                m += 1;
                nh += g(m);
                nm += g(m + 1);
                for j in 2..l {
                    pe.push(t.get(m + j).copied().unwrap_or(0));
                }
                m += l;
            }
            pe.sort_unstable();
            (nh, nm, pe)
        };
        let (h1, m1, pe1) = collect(&i1.n_tautomer, i1.len_tautomer);
        let (h2, m2, pe2) = collect(&i2.n_tautomer, i2.len_tautomer);
        icr.num_taut_h1 = h1;
        icr.num_taut_h2 = h2;
        icr.num_taut_m1 = m1;
        icr.num_taut_m2 = m2;
        let (mut j1, mut j2) = (0, 0);
        let (mut num_eq, mut only1, mut only2) = (0, 0, 0);
        while j1 < pe1.len() && j2 < pe2.len() {
            if pe1[j1] == pe2[j2] {
                j1 += 1;
                j2 += 1;
                num_eq += 1;
            } else if pe1[j1] < pe2[j2] {
                if icr.endp_in1_only.len() < ICR_MAX_ENDP_IN1_ONLY {
                    icr.endp_in1_only.push(pe1[j1]);
                }
                j1 += 1;
                only1 += 1;
            } else {
                if icr.endp_in2_only.len() < ICR_MAX_ENDP_IN2_ONLY {
                    icr.endp_in2_only.push(pe2[j2]);
                }
                j2 += 1;
                only2 += 1;
            }
        }
        for &e in &pe1[j1..] {
            if icr.endp_in1_only.len() < ICR_MAX_ENDP_IN1_ONLY {
                icr.endp_in1_only.push(e);
            }
            only1 += 1;
        }
        for &e in &pe2[j2..] {
            if icr.endp_in2_only.len() < ICR_MAX_ENDP_IN2_ONLY {
                icr.endp_in2_only.push(e);
            }
            only2 += 1;
        }
        if only1 != 0 {
            ret |= IDIF_EXTRA_TG_ENDP;
        }
        if only2 != 0 {
            ret |= IDIF_MISS_TG_ENDP;
        }
        if !(only1 == 0 && only2 == 0 && num_eq != 0) {
            ret |= IDIF_DIFF_TG_ENDP;
        }
    }

    let (lt1, lt2) = (i1.len_tautomer, i2.len_tautomer);
    if lt1 > 1
        && lt2 > 1
        && (lt1 != lt2 || i1.n_tautomer[..lt1 as usize] != i2.n_tautomer[..lt1 as usize])
    {
        ret |= IDIF_TG;
    }
    if i1.n_number_of_isotopic_atoms != i2.n_number_of_isotopic_atoms {
        ret |= IDIF_NUM_ISO_AT;
    } else {
        let ni = i1.n_number_of_isotopic_atoms.max(0) as usize;
        if ni > 0 && i1.isotopic_atom[..ni] != i2.isotopic_atom[..ni] {
            ret |= IDIF_ISO_AT;
        }
    }
    if i1.n_total_charge != i2.n_total_charge {
        ret |= IDIF_CHARGE;
    }
    if let Some(a1) = a1 {
        if a1.n_num_removed_protons != 0
            && a2.is_none_or(|a2| a2.n_num_removed_protons != a1.n_num_removed_protons)
        {
            ret |= IDIF_REM_PROT;
        }
        if a2.is_none_or(|a2| a2.n_num_removed_isotopic_h != a1.n_num_removed_isotopic_h) {
            ret |= IDIF_REM_ISO_H;
        }
    }
    let pick = |i: &InChI| match i.stereo_isotopic.as_ref() {
        Some(s) if s.n_number_of_stereo_bonds + s.n_number_of_stereo_centers != 0 => {
            Some(s.clone())
        }
        _ => i.stereo.clone(),
    };
    let (s1, s2) = (pick(i1), pick(i2));
    ret | compare_reversed_stereo_inchi2(s1.as_ref(), s2.as_ref(), icr)
}
