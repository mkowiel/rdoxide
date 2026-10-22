//! The InChI -> structure restoration, part 7 (`INCHI_BASE/src/ichirvr7.c`):
//! the helpers for the component lists and the comparison messages.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichimake::compare_reversed_inchi;
use super::ichiread::{INCHI_BAS, INCHI_REC, InpInChI, TAUT_NON, TAUT_YES};

/// `INCHIDIFF` bits of the comparison of two InChI
pub mod diff {
    pub const PROBLEM: u64 = 0x0000_0001;
    pub const POSITION_H: u64 = 0x0000_0002;
    pub const MORE_FH: u64 = 0x0000_0004;
    pub const LESS_FH: u64 = 0x0000_0004;
    pub const MORE_H: u64 = 0x0000_0008;
    pub const LESS_H: u64 = 0x0000_0008;
    pub const NO_TAUT: u64 = 0x0000_0010;
    pub const WRONG_TAUT: u64 = 0x0000_0020;
    pub const SINGLE_TG: u64 = 0x0000_0040;
    pub const MULTIPLE_TG: u64 = 0x0000_0080;
    pub const EXTRA_TG_ENDP: u64 = 0x0000_0100;
    pub const MISS_TG_ENDP: u64 = 0x0000_0100;
    pub const DIFF_TG_ENDP: u64 = 0x0000_0100;
    pub const NUM_TG: u64 = 0x0000_0200;
    pub const TG: u64 = 0x0000_0200;
    pub const NUM_ISO_AT: u64 = 0x0000_0400;
    pub const ISO_AT: u64 = 0x0000_0400;
    pub const REM_ISO_H: u64 = 0x0000_0800;
    pub const MOB_ISO_H: u64 = 0x0000_1000;
    pub const CHARGE: u64 = 0x0000_2000;
    pub const REM_PROT: u64 = 0x0000_4000;
    pub const MOBH_PROTONS: u64 = 0x0000_8000;
    pub const SC_INV: u64 = 0x0001_0000;
    pub const SC_PARITY: u64 = 0x0002_0000;
    pub const SC_EXTRA_UNDF: u64 = 0x0004_0000;
    pub const SC_EXTRA: u64 = 0x0008_0000;
    pub const SC_MISS_UNDF: u64 = 0x0010_0000;
    pub const SC_MISS: u64 = 0x0020_0000;
    pub const SB_PARITY: u64 = 0x0040_0000;
    pub const SB_EXTRA_UNDF: u64 = 0x0080_0000;
    pub const SB_EXTRA: u64 = 0x0100_0000;
    pub const SB_MISS_UNDF: u64 = 0x0200_0000;
    pub const SB_MISS: u64 = 0x0400_0000;
    pub const COMP_HLAYER: u64 = 0x0800_0000;
    pub const COMP_NUMBER: u64 = 0x1000_0000;
    pub const STR2INCHI_ERR: u64 = 0x2000_0000;
}

/// `CMP_INCHI_MSG_GROUP_ID`
#[allow(clippy::enum_variant_names)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Grp {
    Err = 1,
    H,
    MobGrp,
    IsoAt,
    Charge,
    Protons,
    IsoH,
    Sc,
    Sb,
    HLayer,
    Comp,
    ConvErr,
}

/// `CompareInchiMsgsGroup`
const GROUPS: [(Grp, &str); 12] = [
    (Grp::Err, " Error:"),
    (Grp::H, " Hydrogens:"),
    (Grp::MobGrp, " Mobile-H groups:"),
    (Grp::IsoAt, " Isotopic:"),
    (Grp::Charge, " Charge(s):"),
    (Grp::Protons, " Proton balance:"),
    (Grp::IsoH, " Exchangeable isotopic H:"),
    (Grp::Sc, " Stereo centers/allenes:"),
    (Grp::Sb, " Stereobonds/cumulenes:"),
    (Grp::HLayer, " Fixed-H layer:"),
    (Grp::Comp, " Number of components:"),
    (Grp::ConvErr, " Conversion encountered:"),
];

/// `CompareInchiMsgs`
const MSGS: [(u64, Grp, &str); 35] = [
    (diff::PROBLEM, Grp::Err, " Wrong result"),
    (diff::POSITION_H, Grp::H, " Locations or number"),
    (diff::MORE_FH, Grp::H, " Fixed-H"),
    (diff::LESS_FH, Grp::H, " Fixed-H"),
    (diff::MORE_H, Grp::H, " Number"),
    (diff::LESS_H, Grp::H, " Number"),
    (diff::NO_TAUT, Grp::MobGrp, " Missing"),
    (diff::WRONG_TAUT, Grp::MobGrp, " Falsely present"),
    (diff::SINGLE_TG, Grp::MobGrp, " One instead of multiple"),
    (diff::MULTIPLE_TG, Grp::MobGrp, " Multiple instead of one"),
    (diff::EXTRA_TG_ENDP, Grp::MobGrp, " Attachment points"),
    (diff::MISS_TG_ENDP, Grp::MobGrp, " Attachment points"),
    (diff::DIFF_TG_ENDP, Grp::MobGrp, " Attachment points"),
    (diff::NUM_TG, Grp::MobGrp, " Number"),
    (diff::TG, Grp::MobGrp, " Do not match"),
    (diff::NUM_ISO_AT, Grp::IsoAt, " Atoms do not match"),
    (diff::ISO_AT, Grp::IsoAt, " Atoms do not match"),
    (
        diff::REM_ISO_H,
        Grp::IsoH,
        " Does not match for a component",
    ),
    (diff::MOB_ISO_H, Grp::IsoH, " Do not match"),
    (diff::CHARGE, Grp::Charge, " Do not match"),
    (
        diff::REM_PROT,
        Grp::Protons,
        " Does not match for a component",
    ),
    (diff::MOBH_PROTONS, Grp::Protons, " Does not match"),
    (diff::SC_INV, Grp::Sc, " Falsely inverted"),
    (diff::SC_PARITY, Grp::Sc, " Wrong parity"),
    (diff::SC_EXTRA_UNDF, Grp::Sc, " Extra undefined"),
    (diff::SC_EXTRA, Grp::Sc, " Extra known"),
    (diff::SC_MISS_UNDF, Grp::Sc, " Missing undefined"),
    (diff::SC_MISS, Grp::Sc, " Missing known"),
    (diff::SB_PARITY, Grp::Sb, " Wrong parity"),
    (diff::SB_EXTRA_UNDF, Grp::Sb, " Extra undefined"),
    (diff::SB_EXTRA, Grp::Sb, " Missing known"),
    (diff::SB_MISS_UNDF, Grp::Sb, " Missing undefined"),
    (diff::SB_MISS, Grp::Sb, " Missing known"),
    (diff::COMP_HLAYER, Grp::HLayer, " Missing or extra"),
    (diff::COMP_NUMBER, Grp::Comp, " Does not match"),
];

/// the last message of the table (`INCHIDIFF_STR2INCHI_ERR`) is kept apart
/// to keep the array length a constant
const MSG_LAST: (u64, Grp, &str) = (diff::STR2INCHI_ERR, Grp::ConvErr, " Error");

/// `AddOneMsg`: appends `add` (after `delim` when the message is not empty)
/// to `msg` of the used length `used_len`; the message is cut with an
/// ellipsis when it does not fit into `tot_len`. Returns the new length.
pub fn add_one_msg(
    msg: &mut String,
    used_len: usize,
    tot_len: i32,
    add: &str,
    delim: Option<&str>,
) -> usize {
    const ELLIP: &str = "...";
    msg.truncate(used_len);
    let mut used = used_len as i32;
    let len = add.len() as i32;
    let len_delim = if used != 0 {
        delim.map_or(0, |d| d.len() as i32)
    } else {
        0
    };
    if len + len_delim + used < tot_len {
        if len_delim != 0 {
            msg.push_str(delim.unwrap_or(""));
            used += len_delim;
        }
        msg.push_str(add);
        used += len;
    } else {
        let len_to_copy = tot_len - used - len_delim - (ELLIP.len() as i32 + 1);
        if len_to_copy > 10 {
            if len_delim != 0 {
                msg.push_str(delim.unwrap_or(""));
                used += len_delim;
            }
            let n = (len_to_copy as usize).min(add.len());
            msg.push_str(&add[..n]);
            used += n as i32;
            msg.push_str(ELLIP);
            used += 3;
        }
    }
    used as usize
}

/// `FillOutCompareMessage`: appends the description of the difference bits
/// `bits[TAUT_NON]`, `bits[TAUT_YES]` to `msg` (total length `n_len_msg`).
/// Returns -1 when there are differences.
pub fn fill_out_compare_message(msg: &mut String, n_len_msg: i32, bits: [u64; 2]) -> i32 {
    const HDR: &str = " Problems/mismatches:";
    let mut len = msg.len();
    if bits[TAUT_YES] == 0 && bits[TAUT_NON] == 0 {
        return 0;
    }
    if !msg.contains(HDR) {
        len = add_one_msg(msg, len, n_len_msg, HDR, None);
    }
    for b_mobile_h in [TAUT_YES, TAUT_NON] {
        if bits[b_mobile_h] != 0 {
            let one = if b_mobile_h == TAUT_YES {
                " Mobile-H("
            } else {
                " Fixed-H("
            };
            len = add_one_msg(msg, len, n_len_msg, one, None);
        }
        let mut prev_grp: i32 = -1;
        for sh in 0..64 {
            let bit = 1u64 << sh;
            if bit & bits[b_mobile_h] == 0 {
                continue;
            }
            // search for the message
            for &(m_bit, m_grp, text) in MSGS.iter().chain(std::iter::once(&MSG_LAST)) {
                if bit & m_bit == 0 {
                    continue;
                }
                if let Some(n) = GROUPS.iter().position(|g| g.0 == m_grp) {
                    let cur = n as i32;
                    if cur != prev_grp {
                        if prev_grp >= 0 {
                            len = add_one_msg(msg, len, n_len_msg, ";", None);
                        }
                        len = add_one_msg(msg, len, n_len_msg, GROUPS[n].1, None);
                    }
                    len = add_one_msg(
                        msg,
                        len,
                        n_len_msg,
                        text,
                        if cur == prev_grp { Some(",") } else { None },
                    );
                    prev_grp = cur;
                    break;
                }
            }
        }
        if bits[b_mobile_h] != 0 {
            len = add_one_msg(msg, len, n_len_msg, ")", None);
        }
    }
    -1
}

/// `RemoveFixHInChIIdentical2MobH`: eliminates the Fixed-H InChI that are
/// exactly the same as the corresponding Mobile-H ones
pub fn remove_fix_h_inchi_identical2_mob_h(inp: &mut InpInChI) {
    for rec in 0..2 {
        let cur = inp.n_num_components[rec][TAUT_YES].min(inp.n_num_components[rec][TAUT_NON]);
        for k in 0..cur.max(0) as usize {
            let same = compare_reversed_inchi(
                Some(&inp.p_inp_inchi[rec][TAUT_YES][k]),
                Some(&inp.p_inp_inchi[rec][TAUT_NON][k]),
                None,
                None,
            ) == 0;
            if same {
                inp.p_inp_inchi[rec][TAUT_NON][k] = Default::default();
            }
        }
    }
}

/// `MarkDisconectedIdenticalToReconnected`: marks the Disconnected InChI
/// components that are exactly the same as the Reconnected ones. The
/// Disconnected one gets the negative number of the reconnected component,
/// the Reconnected one the positive number of the disconnected one. Returns
/// the number of the marked components.
pub fn mark_disconnected_identical_to_reconnected(inp: &mut InpInChI) -> i32 {
    let mut num_marked = 0;
    let n = &inp.n_num_components;
    let k1max = n[INCHI_BAS][TAUT_YES].max(n[INCHI_BAS][TAUT_NON]);
    for k1 in 0..k1max.max(0) as usize {
        let n = inp.n_num_components;
        let k2max = n[INCHI_REC][TAUT_YES].max(n[INCHI_REC][TAUT_NON]);
        for k2 in 0..k2max.max(0) as usize {
            let (k1i, k2i) = (k1 as i32, k2 as i32);
            let p = &inp.p_inp_inchi;
            let eq_m = k1i < n[INCHI_BAS][TAUT_YES]
                && k2i < n[INCHI_REC][TAUT_YES]
                && p[INCHI_REC][TAUT_YES][k2].n_link == 0
                && p[INCHI_BAS][TAUT_YES][k1].b_deleted == 0
                && p[INCHI_BAS][TAUT_YES][k1].n_number_of_atoms != 0
                && p[INCHI_BAS][TAUT_YES][k1].n_number_of_atoms
                    == p[INCHI_REC][TAUT_YES][k2].n_number_of_atoms
                && p[INCHI_REC][TAUT_YES][k2].b_deleted == 0
                && compare_reversed_inchi(
                    Some(&p[INCHI_REC][TAUT_YES][k2]),
                    Some(&p[INCHI_BAS][TAUT_YES][k1]),
                    None,
                    None,
                ) == 0;
            let is_f1 = k1i < n[INCHI_BAS][TAUT_NON]
                && p[INCHI_BAS][TAUT_NON][k1].b_deleted == 0
                && p[INCHI_BAS][TAUT_NON][k1].n_number_of_atoms > 0;
            let is_f2 = k2i < n[INCHI_REC][TAUT_NON]
                && p[INCHI_REC][TAUT_NON][k2].b_deleted == 0
                && p[INCHI_REC][TAUT_NON][k2].n_number_of_atoms > 0;
            let eq_f = is_f1
                && is_f2
                && p[INCHI_REC][TAUT_NON][k2].n_link == 0
                && p[INCHI_BAS][TAUT_NON][k1].n_number_of_atoms
                    == p[INCHI_REC][TAUT_NON][k2].n_number_of_atoms
                && compare_reversed_inchi(
                    Some(&p[INCHI_REC][TAUT_NON][k2]),
                    Some(&p[INCHI_BAS][TAUT_NON][k1]),
                    None,
                    None,
                ) == 0;
            if eq_m && ((!is_f1 && !is_f2) || eq_f) {
                let p = &mut inp.p_inp_inchi;
                p[INCHI_BAS][TAUT_YES][k1].n_link = -(k2i + 1);
                p[INCHI_REC][TAUT_YES][k2].n_link = k1i + 1;
                if eq_f {
                    p[INCHI_BAS][TAUT_NON][k1].n_link = -(k2i + 1);
                    p[INCHI_REC][TAUT_NON][k2].n_link = k1i + 1;
                }
                num_marked += 1;
                break;
            }
        }
    }
    num_marked
}
