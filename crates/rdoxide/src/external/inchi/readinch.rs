//! `Extract0DParities` of `INCHI_BASE/src/readinch.c`: the 0D stereo
//! descriptors of the user's input become atom and bond parities of the
//! `inp_ATOM`s.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichierr::add_error_message;
use super::ichister::*;
use super::incomdef::BOND_TYPE_DOUBLE;
use super::inpdef::{
    InpAtom, MAX_ATOMS, MAX_NUM_STEREO_BONDS, SB_PARITY_FLAG, SB_PARITY_MASK, SB_PARITY_SHFT,
};

// inchi_api.h
pub const NO_ATOM: i16 = -1;
pub const INCHI_STEREO_TYPE_NONE: i8 = 0;
pub const INCHI_STEREO_TYPE_DOUBLE_BOND: i8 = 1;
pub const INCHI_STEREO_TYPE_TETRAHEDRAL: i8 = 2;
pub const INCHI_STEREO_TYPE_ALLENE: i8 = 3;
pub const INCHI_PARITY_NONE: i32 = 0;
pub const INCHI_PARITY_ODD: i32 = 1;
pub const INCHI_PARITY_EVEN: i32 = 2;
pub const INCHI_PARITY_UNKNOWN: i32 = 3;
pub const INCHI_PARITY_UNDEFINED: i32 = 4;

/// `inchi_Stereo0D`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InchiStereo0D {
    pub neighbor: [i16; 4],
    /// central tetrahedral atom or a central atom of allene; otherwise `NO_ATOM`
    pub central_atom: i16,
    pub type_: i8,
    /// parity of the connected structure | (parity of the disconnected << 3)
    pub parity: i8,
}

fn pos_in(list: &[u16], next: u16, len: usize) -> Option<usize> {
    list[..len.min(list.len())].iter().position(|&a| a == next)
}

/// `Extract0DParities`. `v_ab_parity_unknown` is the internal value for
/// "unknown" parity (the same as undefined by default, or its own value with
/// `SLUUD`). The messages go to `p_str_err`.
pub fn extract_0d_parities(
    at: &mut [InpAtom],
    stereo_0d: &[InchiStereo0D],
    p_str_err: &mut String,
    err: &mut i32,
    v_ab_parity_unknown: i32,
) -> i32 {
    let n_num_atoms = at.len() as i32;
    if stereo_0d.is_empty() {
        return 0;
    }
    for (i0d, s) in stereo_0d.iter().enumerate() {
        let parity = i32::from(s.parity) & SB_PARITY_MASK;
        let parity_nm = (i32::from(s.parity) & i32::from(SB_PARITY_FLAG)) >> SB_PARITY_SHFT;

        let wrong = |p_str_err: &mut String| {
            add_error_message(p_str_err, "Wrong 0D stereo descriptor(s):");
            add_error_message(p_str_err, &format!("#{}", i0d + 1));
        };
        if parity == INCHI_PARITY_NONE
            || (parity != INCHI_PARITY_ODD
                && parity != INCHI_PARITY_EVEN
                && parity != INCHI_PARITY_UNKNOWN
                && parity != INCHI_PARITY_UNDEFINED)
        {
            wrong(p_str_err);
            continue; // warning
        }

        let type_ = s.type_;
        let a2 = i32::from(s.central_atom);
        let mut j: i32 = -1;
        let (mut sb_ord_from_i1, mut sb_ord_from_i2, mut sn_ord_from_i1, mut sn_ord_from_i2) =
            (-1i32, -1i32, -1i32, -1i32);
        let (mut i1n, mut i2n, mut i1, mut i2) =
            (MAX_ATOMS + 1, MAX_ATOMS + 1, MAX_ATOMS + 1, MAX_ATOMS + 1);

        if type_ == INCHI_STEREO_TYPE_TETRAHEDRAL
            || (type_ == INCHI_STEREO_TYPE_ALLENE && 0 <= a2 && a2 < n_num_atoms)
            || (type_ == INCHI_STEREO_TYPE_DOUBLE_BOND && a2 == i32::from(NO_ATOM))
        {
            // test the quadruplet
            j = 0;
            let mut k_prev = -1i32;
            while j < 4 {
                let k = i32::from(s.neighbor[j as usize]);
                if k < 0 || k >= n_num_atoms || k_prev == k {
                    break;
                }
                // tetrahedral atom connectivity test (a central atom outside
                // the structure is undefined behaviour in C: treated as a
                // failure of the test)
                if type_ == INCHI_STEREO_TYPE_TETRAHEDRAL && k != a2 {
                    if a2 < 0 || a2 >= n_num_atoms {
                        break;
                    }
                    let c = &at[a2 as usize];
                    if pos_in(&c.neighbor, k as u16, c.valence as usize).is_none() {
                        break;
                    }
                }
                k_prev = k;
                j += 1;
            }
        }

        // Find in the adjacency lists the double bond neighbor that leads to the opposite atom
        if j == 4 && (type_ == INCHI_STEREO_TYPE_ALLENE || type_ == INCHI_STEREO_TYPE_DOUBLE_BOND) {
            i1n = s.neighbor[0] as usize;
            i1 = s.neighbor[1] as usize;
            i2 = s.neighbor[2] as usize;
            i2n = s.neighbor[3] as usize;

            let q1 = pos_in(&at[i1].neighbor, i1n as u16, at[i1].valence as usize);
            let q2 = pos_in(&at[i2].neighbor, i2n as u16, at[i2].valence as usize);
            if q1.is_none() || q2.is_none() {
                j = -2;
            } else {
                let (q1, q2) = (q1.unwrap_or(0), q2.unwrap_or(0));
                let p1 = pos_in(&at[i1].neighbor, i2 as u16, at[i1].valence as usize);
                if p1.is_none() {
                    // at[i1] and at[i2] are not connected: allene or cumulene
                    let (mut cur, mut next) = (i1, i1);
                    let mut prev = i1;
                    let mut len = 0;
                    let mut half_len = 0;
                    let mut num_dbond = 0;
                    let mut next_ord = 0usize;
                    while len < 20 {
                        prev = cur;
                        cur = next;
                        num_dbond = 0;
                        for i in 0..at[cur].valence as usize {
                            // follow double bond path && avoid going back
                            if i32::from(at[cur].bond_type[i]) == BOND_TYPE_DOUBLE
                                && prev != at[cur].neighbor[i] as usize
                            {
                                next = at[cur].neighbor[i] as usize;
                                next_ord = i;
                                num_dbond += 1;
                            }
                        }
                        if num_dbond == 1 && next != i1 {
                            len += 1;
                            if len == 1 {
                                sb_ord_from_i1 = next_ord as i32;
                            }
                            if type_ == INCHI_STEREO_TYPE_ALLENE && next as i32 == a2 {
                                half_len = len;
                            }
                        } else {
                            break;
                        }
                    }
                    let p2 = pos_in(&at[i2].neighbor, prev as u16, at[i2].valence as usize);
                    if cur == i2
                        && prev != cur
                        && num_dbond == 0
                        && len > 1
                        && p2.is_some()
                        && (type_ != INCHI_STEREO_TYPE_ALLENE || len == 2 * half_len)
                    {
                        sb_ord_from_i2 = p2.unwrap_or(0) as i32;
                        sn_ord_from_i1 = q1 as i32;
                        sn_ord_from_i2 = q2 as i32;
                    } else {
                        j = -5;
                    }
                } else if type_ == INCHI_STEREO_TYPE_ALLENE {
                    // allene must have been already processed, otherwise error
                    j = -3;
                } else {
                    // double bond only; the bond type is not checked because at
                    // the end of the normalization it may be alternating
                    let p2 = pos_in(&at[i2].neighbor, i1 as u16, at[i2].valence as usize);
                    if type_ == INCHI_STEREO_TYPE_DOUBLE_BOND && p2.is_some() {
                        sb_ord_from_i1 = p1.unwrap_or(0) as i32;
                        sb_ord_from_i2 = p2.unwrap_or(0) as i32;
                        sn_ord_from_i1 = q1 as i32;
                        sn_ord_from_i2 = q2 as i32;
                    } else {
                        j = -4;
                    }
                }
            }
        }

        if j != 4 {
            wrong(p_str_err);
            continue; // error
        }

        match type_ {
            INCHI_STEREO_TYPE_NONE => continue,
            INCHI_STEREO_TYPE_DOUBLE_BOND | INCHI_STEREO_TYPE_ALLENE => {
                let mut j1 = 0;
                while j1 < MAX_NUM_STEREO_BONDS && at[i1].sb_parity[j1] != 0 {
                    j1 += 1;
                }
                let mut j2 = 0;
                while j2 < MAX_NUM_STEREO_BONDS && at[i2].sb_parity[j2] != 0 {
                    j2 += 1;
                }
                if j1 < MAX_NUM_STEREO_BONDS
                    && j2 < MAX_NUM_STEREO_BONDS
                    && sb_ord_from_i1 >= 0
                    && sb_ord_from_i2 >= 0
                    && sn_ord_from_i1 >= 0
                    && sn_ord_from_i2 >= 0
                {
                    let (p_1, p_2): (i32, i32) = match parity {
                        INCHI_PARITY_ODD => (AB_PARITY_ODD, AB_PARITY_EVEN),
                        INCHI_PARITY_EVEN => (AB_PARITY_ODD, AB_PARITY_ODD),
                        INCHI_PARITY_UNDEFINED => (AB_PARITY_UNDF, AB_PARITY_UNDF),
                        INCHI_PARITY_UNKNOWN => (v_ab_parity_unknown, v_ab_parity_unknown),
                        _ => (AB_PARITY_NONE, AB_PARITY_NONE),
                    };
                    at[i1].sb_parity[j1] = p_1 as i8;
                    at[i2].sb_parity[j2] = p_2 as i8;
                    match parity_nm {
                        INCHI_PARITY_ODD => {
                            at[i1].sb_parity[j1] |= (AB_PARITY_ODD << SB_PARITY_SHFT) as i8;
                            at[i2].sb_parity[j2] |= (AB_PARITY_EVEN << SB_PARITY_SHFT) as i8;
                        }
                        INCHI_PARITY_EVEN => {
                            at[i1].sb_parity[j1] |= (AB_PARITY_ODD << SB_PARITY_SHFT) as i8;
                            at[i2].sb_parity[j2] |= (AB_PARITY_ODD << SB_PARITY_SHFT) as i8;
                        }
                        INCHI_PARITY_UNDEFINED => {
                            at[i1].sb_parity[j1] |= (AB_PARITY_UNDF << SB_PARITY_SHFT) as i8;
                            at[i2].sb_parity[j2] |= (AB_PARITY_UNDF << SB_PARITY_SHFT) as i8;
                        }
                        INCHI_PARITY_UNKNOWN => {
                            at[i1].sb_parity[j1] |= (v_ab_parity_unknown << SB_PARITY_SHFT) as i8;
                            at[i2].sb_parity[j2] |= (v_ab_parity_unknown << SB_PARITY_SHFT) as i8;
                        }
                        _ => {}
                    }
                    at[i1].sb_ord[j1] = sb_ord_from_i1 as i8;
                    at[i1].sn_ord[j1] = sn_ord_from_i1 as i8;
                    at[i1].sn_orig_at_num[j1] = at[i1n].orig_at_number;

                    at[i2].sb_ord[j2] = sb_ord_from_i2 as i8;
                    at[i2].sn_ord[j2] = sn_ord_from_i2 as i8;
                    at[i2].sn_orig_at_num[j2] = at[i2n].orig_at_number;
                }
            }
            INCHI_STEREO_TYPE_TETRAHEDRAL => {
                let a2u = a2 as usize;
                match parity {
                    INCHI_PARITY_ODD => at[a2u].p_parity = AB_PARITY_ODD as i8,
                    INCHI_PARITY_EVEN => at[a2u].p_parity = AB_PARITY_EVEN as i8,
                    INCHI_PARITY_UNDEFINED => at[a2u].p_parity = AB_PARITY_UNDF as i8,
                    INCHI_PARITY_UNKNOWN => at[a2u].p_parity = v_ab_parity_unknown as i8,
                    _ => continue,
                }
                for j in 0..4 {
                    let k = s.neighbor[j] as usize;
                    at[a2u].p_orig_at_num[j] = at[k].orig_at_number;
                }
            }
            _ => {}
        }
    }
    // take care of Unknown stereobonds: copy their Unknown stereo descriptors
    // to at->bond_stereo
    fix_unkn_0d_stereo_bonds(at);

    // TARGET_API_LIB
    let k = reconcile_all_cml_bond_parities(at, false);
    if k != 0 {
        add_error_message(p_str_err, "0D Parities Reconciliation failed:");
        add_error_message(p_str_err, &k.to_string());
    }
    let _ = err;
    0
}
