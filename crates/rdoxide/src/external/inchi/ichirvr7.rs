//! The InChI -> structure restoration, part 7 (`INCHI_BASE/src/ichirvr7.c`):
//! the helpers for the component lists and the comparison messages.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::{InChI, InChIAux, InChIStereo};
use super::ichimake::{compare_hill_formulas_no_h, compare_reversed_inchi};
use super::ichiread::{INCHI_BAS, INCHI_REC, InpInChI, TAUT_NON, TAUT_YES};
use super::ichister::AB_PARITY_UNDF;

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

/// the merge-compare of two sorted index lists of `CompareReversedStereoINChI3`:
/// returns (different parity, only in the first, undefined only in the first,
/// only in the second, undefined only in the second)
fn merge_stereo<K: PartialOrd + Copy>(
    k1: &[K],
    p1: &[i8],
    k2: &[K],
    p2: &[i8],
) -> (i32, i32, i32, i32, i32) {
    let undf = AB_PARITY_UNDF as i8;
    let (mut dif, mut in1, mut undf1, mut in2, mut undf2) = (0, 0, 0, 0, 0);
    let (mut j1, mut j2) = (0, 0);
    while j1 < k1.len() && j2 < k2.len() {
        if k1[j1] == k2[j2] {
            dif += i32::from(p1[j1] != p2[j2]);
            j1 += 1;
            j2 += 1;
        } else if k1[j1] < k2[j2] {
            in1 += 1;
            undf1 += i32::from(p1[j1] == undf);
            j1 += 1;
        } else {
            in2 += 1;
            undf2 += i32::from(p2[j2] == undf);
            j2 += 1;
        }
    }
    for &p in &p1[j1..k1.len()] {
        in1 += 1;
        undf1 += i32::from(p == undf);
    }
    for &p in &p2[j2..k2.len()] {
        in2 += 1;
        undf2 += i32::from(p == undf);
    }
    (dif, in1, undf1, in2, undf2)
}

/// `CompareReversedStereoINChI3` (the lists of the differing centres and
/// bonds collected in the `ICR` are not kept)
pub fn compare_reversed_stereo_inchi3(s1: Option<&InChIStereo>, s2: Option<&InChIStereo>) -> u64 {
    let mut ret = 0;
    let nsc =
        |s: Option<&InChIStereo>| s.map_or(0, |s| s.n_number_of_stereo_centers.max(0) as usize);
    let nsb = |s: Option<&InChIStereo>| s.map_or(0, |s| s.n_number_of_stereo_bonds.max(0) as usize);
    let (c1, c2, b1, b2) = (nsc(s1), nsc(s2), nsb(s1), nsb(s2));
    if c1 != 0 || c2 != 0 {
        let differs = c1 != c2 || {
            let (a, b) = (s1.unwrap(), s2.unwrap());
            a.n_number[..c1] != b.n_number[..c1] || a.t_parity[..c1] != b.t_parity[..c1]
        };
        if differs {
            let (a, b) = (
                s1.map(|s| (&s.n_number[..c1], &s.t_parity[..c1])),
                s2.map(|s| (&s.n_number[..c2], &s.t_parity[..c2])),
            );
            let (k1, p1) = a.unwrap_or((&[], &[]));
            let (k2, p2) = b.unwrap_or((&[], &[]));
            let (dif, in1, undf1, in2, undf2) = merge_stereo(k1, p1, k2, p2);
            if dif != 0 {
                ret |= diff::SC_PARITY;
            }
            if in1 != 0 {
                if undf1 != 0 {
                    ret |= diff::SC_EXTRA_UNDF;
                }
                if in1 != undf1 {
                    ret |= diff::SC_EXTRA;
                }
            }
            if in2 != 0 {
                if undf2 != 0 {
                    ret |= diff::SC_MISS_UNDF;
                }
                if in2 != undf2 {
                    ret |= diff::SC_MISS;
                }
            }
        }
    }
    if let (Some(a), Some(b)) = (s1, s2) {
        if b.n_comp_inv2_abs != 2
            && a.n_comp_inv2_abs != b.n_comp_inv2_abs
            && a.n_comp_inv2_abs != 0
            && b.n_comp_inv2_abs != 0
        {
            ret |= diff::SC_INV;
        }
    }
    if b1 != 0 || b2 != 0 {
        let differs = b1 != b2 || {
            let (a, b) = (s1.unwrap(), s2.unwrap());
            a.n_bond_atom1[..b1] != b.n_bond_atom1[..b1]
                || a.n_bond_atom2[..b1] != b.n_bond_atom2[..b1]
                || a.b_parity[..b1] != b.b_parity[..b1]
        };
        if differs {
            let keys = |s: Option<&InChIStereo>, n: usize| -> Vec<(u16, u16)> {
                s.map_or(Vec::new(), |s| {
                    (0..n)
                        .map(|j| (s.n_bond_atom1[j], s.n_bond_atom2[j]))
                        .collect()
                })
            };
            let (k1, k2) = (keys(s1, b1), keys(s2, b2));
            let p1: &[i8] = s1.map_or(&[], |s| &s.b_parity[..b1]);
            let p2: &[i8] = s2.map_or(&[], |s| &s.b_parity[..b2]);
            // lexicographic order of the pairs, as in C
            let (dif, in1, undf1, in2, undf2) = merge_stereo(&k1, p1, &k2, p2);
            if dif != 0 {
                ret |= diff::SB_PARITY;
            }
            if in1 != 0 {
                if undf1 != 0 {
                    ret |= diff::SB_EXTRA_UNDF;
                }
                if in1 != undf1 {
                    ret |= diff::SB_EXTRA;
                }
            }
            if in2 != 0 {
                if undf2 != 0 {
                    ret |= diff::SB_MISS_UNDF;
                }
                if in2 != undf2 {
                    ret |= diff::SB_MISS;
                }
            }
        }
    }
    ret
}

/// `CompareReversedINChI3`: the bits of `diff` that describe how the InChI
/// `i1` (from the reversed structure) differs from the input InChI `i2`
pub fn compare_reversed_inchi3(
    i1: Option<&InChI>,
    i2: Option<&InChI>,
    a1: Option<&InChIAux>,
    a2: Option<&InChIAux>,
) -> u64 {
    let (i1, i2) = match (i1, i2) {
        (None, None) => return 0,
        (Some(a), Some(b)) => (a, b),
        _ => return diff::PROBLEM,
    };
    if i1.n_error_code != i2.n_error_code || i1.n_error_code != 0 {
        return diff::PROBLEM;
    }
    let mut ret = 0;
    if i1.n_number_of_atoms != i2.n_number_of_atoms {
        return diff::PROBLEM; // INCHIDIFF_NUM_AT
    }
    let n = i1.n_number_of_atoms.max(0) as usize;
    if n > 0 {
        if i1.n_atom[..n] != i2.n_atom[..n] {
            return diff::PROBLEM; // INCHIDIFF_ATOMS
        }
        if i1.n_num_h[..n] != i2.n_num_h[..n] {
            ret |= diff::POSITION_H;
        }
        if i1.n_num_h_fixed.is_some() || i2.n_num_h_fixed.is_some() {
            let has =
                |f: &Option<Vec<i8>>| f.as_ref().is_some_and(|v| v[..n].iter().any(|&x| x != 0));
            let (h1, h2) = (has(&i1.n_num_h_fixed), has(&i2.n_num_h_fixed));
            if h1 && !h2 {
                ret |= diff::MORE_FH;
            } else if !h1 && h2 {
                ret |= diff::LESS_FH;
            } else if h1 && h2 {
                let (f1, f2) = (
                    i1.n_num_h_fixed.as_ref().unwrap(),
                    i2.n_num_h_fixed.as_ref().unwrap(),
                );
                if f1[..n] != f2[..n] {
                    let more1 = (0..n).any(|i| f1[i] > f2[i]);
                    let more2 = (0..n).any(|i| f1[i] < f2[i]);
                    ret |= (if more1 { diff::MORE_FH } else { 0 })
                        | (if more2 { diff::LESS_FH } else { 0 });
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
    if ret2 != 0 {
        return ret | diff::PROBLEM; // INCHIDIFF_NUM_EL
    }
    if num_h1 > num_h2 {
        ret |= diff::MORE_H;
    }
    if num_h1 < num_h2 {
        ret |= diff::LESS_H;
    }
    if i1.len_conn_table != i2.len_conn_table {
        return ret | diff::PROBLEM; // INCHIDIFF_CON_LEN
    }
    let lc = i1.len_conn_table.max(0) as usize;
    if lc > 0 && i1.n_conn_table[..lc] != i2.n_conn_table[..lc] {
        return ret | diff::PROBLEM; // INCHIDIFF_CON_TBL
    }
    let n1 = if i1.len_tautomer != 0 {
        i32::from(i1.n_tautomer[0])
    } else {
        0
    };
    let n2 = if i2.len_tautomer != 0 {
        i32::from(i2.n_tautomer[0])
    } else {
        0
    };
    if n1 == 0 && n2 != 0 {
        ret |= diff::NO_TAUT;
    } else if n1 != 0 && n2 == 0 {
        ret |= diff::WRONG_TAUT;
    } else if n1 == 1 && n2 > 1 {
        ret |= diff::SINGLE_TG;
    } else if n1 > 1 && n2 == 1 {
        ret |= diff::MULTIPLE_TG;
    } else if n1 != n2 {
        ret |= diff::NUM_TG;
    }
    if n1 != 0 || n2 != 0 {
        let collect = |t: &[u16], len: i32| {
            let g = |k: usize| t.get(k).copied().unwrap_or(0);
            let mut pe = Vec::new();
            let mut m = 1usize;
            while (m as i32) < len {
                let l = g(m) as usize;
                m += 1;
                for j in 2..l {
                    pe.push(g(m + j));
                }
                m += l;
            }
            pe.sort_unstable();
            pe
        };
        let (pe1, pe2) = (
            collect(&i1.n_tautomer, i1.len_tautomer),
            collect(&i2.n_tautomer, i2.len_tautomer),
        );
        let (mut j1, mut j2) = (0, 0);
        let (mut num_eq, mut only1, mut only2) = (0, 0, 0);
        while j1 < pe1.len() && j2 < pe2.len() {
            if pe1[j1] == pe2[j2] {
                j1 += 1;
                j2 += 1;
                num_eq += 1;
            } else if pe1[j1] < pe2.get(j1).copied().unwrap_or(0) {
                // as in C the index j1 is used for the second list
                j1 += 1;
                only1 += 1;
            } else {
                j2 += 1;
                only2 += 1;
            }
        }
        only1 += pe1.len() - j1;
        only2 += pe2.len() - j2;
        if only1 != 0 {
            ret |= diff::EXTRA_TG_ENDP;
        }
        if only2 != 0 {
            ret |= diff::MISS_TG_ENDP;
        }
        if !(only1 == 0 && only2 == 0 && num_eq != 0) {
            ret |= diff::DIFF_TG_ENDP;
        }
    }
    let (lt1, lt2) = (i1.len_tautomer, i2.len_tautomer);
    if lt1 > 1
        && lt2 > 1
        && (lt1 != lt2 || i1.n_tautomer[..lt1 as usize] != i2.n_tautomer[..lt1 as usize])
    {
        ret |= diff::TG;
    }
    if i1.n_number_of_isotopic_atoms != i2.n_number_of_isotopic_atoms {
        ret |= diff::NUM_ISO_AT;
    } else {
        let ni = i1.n_number_of_isotopic_atoms.max(0) as usize;
        if ni > 0 && i1.isotopic_atom[..ni] != i2.isotopic_atom[..ni] {
            ret |= diff::ISO_AT;
        }
    }
    if i1.n_total_charge != i2.n_total_charge {
        ret |= diff::CHARGE;
    }
    if let Some(a1) = a1 {
        if a1.n_num_removed_protons != 0
            && a2.is_none_or(|a2| a2.n_num_removed_protons != a1.n_num_removed_protons)
        {
            ret |= diff::REM_PROT;
        }
        if a2.is_none_or(|a2| a2.n_num_removed_isotopic_h != a1.n_num_removed_isotopic_h) {
            ret |= diff::REM_ISO_H;
        }
    }
    let pick = |i: &InChI| match i.stereo_isotopic.as_ref() {
        Some(s) if s.n_number_of_stereo_bonds + s.n_number_of_stereo_centers != 0 => {
            Some(s.clone())
        }
        _ => i.stereo.clone(),
    };
    ret | compare_reversed_stereo_inchi3(pick(i1).as_ref(), pick(i2).as_ref())
}

/// `CompareTwoPairsOfInChI`: ORs the difference bits of the Mobile-H and
/// Fixed-H InChI of two structures into `flags`
pub fn compare_two_pairs_of_inchi(
    p1: [Option<&InChI>; 2],
    p2: [Option<&InChI>; 2],
    flags: &mut [u64; 2],
) {
    for i_mobile_h in 0..2 {
        if p1[i_mobile_h].is_some() != p2[i_mobile_h].is_some() {
            if i_mobile_h == TAUT_NON && p1[TAUT_YES].is_some() && p2[TAUT_YES].is_some() {
                flags[i_mobile_h] |= diff::COMP_HLAYER;
            } else {
                flags[i_mobile_h] |= diff::COMP_NUMBER;
            }
            continue;
        }
        if p1[i_mobile_h].is_some() && p2[i_mobile_h].is_some() {
            flags[i_mobile_h] |=
                compare_reversed_inchi3(p1[i_mobile_h], p2[i_mobile_h], None, None);
        }
    }
}

/// `I2A_FLAG_FIXEDH`, `I2A_FLAG_RECMET`
pub const I2A_FLAG_FIXEDH: i32 = 0x0001;
pub const I2A_FLAG_RECMET: i32 = 0x0002;
/// `EL_NUMBER_ZY`, `EL_NUMBER_ZZ`
const EL_NUMBER_ZY: u8 = 119;
const EL_NUMBER_ZZ: u8 = 120;

/// `InChI2Atom`: restores the component `i_component` of the input InChI
/// (Fixed-H or Mobile-H, disconnected or reconnected as the flag says) as
/// the structure in `ps`. `num_inp_actual` and `polymer` are the
/// `num_inp` and `polymer` members of `OneInput`. Returns 0 when the
/// component does not exist.
#[allow(clippy::too_many_arguments)]
pub fn inchi2_atom(
    ip: &super::ichiparm::InputParms,
    sd: &mut super::runichi3::StructData,
    num_inp: i64,
    ps: &mut super::ichirvr2::StrFromInChI,
    i_component: usize,
    i_at_no_offset: i32,
    i2a_flag: i32,
    b_has_some_fixed_h: bool,
    one_input: &mut InpInChI,
    num_inp_actual: i64,
    polymer: Option<&super::mol2atom::OadPolymer>,
) -> i32 {
    let mut i_inchi = if i2a_flag & I2A_FLAG_RECMET != 0 {
        INCHI_REC
    } else {
        INCHI_BAS
    };
    let mut b_mobile_h = if i2a_flag & I2A_FLAG_FIXEDH != 0 {
        TAUT_NON
    } else {
        TAUT_YES
    };
    let ic = i_component as i32;
    // disconnected or reconnected
    if i_inchi == INCHI_REC && one_input.n_num_components[i_inchi][TAUT_YES] == 0 {
        i_inchi = INCHI_BAS;
    }
    if ic >= one_input.n_num_components[i_inchi][TAUT_YES] {
        return 0; // component does not exist
    }
    // mobile or fixed H
    ps.b_fixed_h_exists = 0;
    if b_mobile_h == TAUT_NON && one_input.n_num_components[i_inchi][b_mobile_h] == 0 {
        // only one InChI exists (no mobile H)
        b_mobile_h = TAUT_YES;
    }
    if ic >= one_input.n_num_components[i_inchi][b_mobile_h] {
        return 0;
    }
    ps.b_mobile_h = b_mobile_h as i32;
    ps.i_inchi = i_inchi as i8;
    // deleted component only in case Mobile-H and compound contains only protons
    if one_input.p_inp_inchi[i_inchi][b_mobile_h][i_component].b_deleted != 0 {
        return 0;
    }
    if b_mobile_h == TAUT_NON
        && let Some(v) = &one_input.n_num_protons[i_inchi][TAUT_YES].p_num_protons
    {
        ps.n_num_removed_protons_mobh_inchi = i32::from(v[i_component][0]);
    }
    let fixed_exists = one_input.p_inp_inchi[i_inchi][TAUT_NON]
        .get(i_component)
        .is_some_and(|c| c.n_number_of_atoms > 0 && c.b_deleted == 0);
    if b_mobile_h == TAUT_NON || (b_mobile_h == TAUT_YES && fixed_exists) {
        ps.b_fixed_h_exists = 1;
    }
    let mobile_exists = b_mobile_h == TAUT_NON
        && ic < one_input.n_num_components[i_inchi][TAUT_YES]
        && one_input.p_inp_inchi[i_inchi][TAUT_YES]
            .get(i_component)
            .is_some_and(|c| c.n_number_of_atoms > 0 && c.b_deleted == 0);
    ps.num_inp_actual = num_inp_actual;
    // non-polymer Zz to Zy
    if let Some(p) = polymer {
        ps.n_pzz = 0;
        ps.n_zy = 0;
        let inchi = &mut one_input.p_inp_inchi[i_inchi][b_mobile_h][i_component];
        for a in 0..inchi.n_number_of_atoms as usize {
            let aglob = i_at_no_offset + a as i32 + 1;
            if inchi.n_atom[a] == EL_NUMBER_ZZ {
                let mut new_num = EL_NUMBER_ZY;
                if p.units
                    .iter()
                    .take(p.n.max(0) as usize)
                    .any(|u| aglob == u.cap1 || aglob == u.cap2)
                {
                    new_num = EL_NUMBER_ZZ;
                }
                inchi.n_atom[a] = new_num;
                if new_num == EL_NUMBER_ZY {
                    ps.n_zy += 1;
                } else {
                    ps.n_pzz += 1;
                }
            }
        }
    }
    let ret = {
        let p0 = &one_input.p_inp_inchi[i_inchi][b_mobile_h][i_component];
        let p1 = mobile_exists.then(|| &one_input.p_inp_inchi[i_inchi][TAUT_YES][i_component]);
        super::ichirvr4::one_inchi2_atom(
            ip,
            sd,
            num_inp,
            ps,
            ic,
            0,
            b_has_some_fixed_h,
            [Some(p0), p1],
        )
    };
    // the Zy atoms are returned to Zz
    if !ps.at.is_empty() && !ps.at2.is_empty() {
        let inchi = &mut one_input.p_inp_inchi[i_inchi][b_mobile_h][i_component];
        for a in 0..inchi.n_number_of_atoms as usize {
            if inchi.n_atom[a] == EL_NUMBER_ZY {
                inchi.n_atom[a] = EL_NUMBER_ZZ;
                for at in [&mut ps.at[a], &mut ps.at2[a]] {
                    at.el_number = EL_NUMBER_ZZ;
                    at.elname = [0; super::inpdef::ATOM_EL_LEN];
                    at.elname[..2].copy_from_slice(b"Zz");
                }
                ps.n_zy -= 1;
                ps.n_pzz += 1;
            }
        }
    }
    ret
}

/// `MergeStructureComponents`: merges the restored components `at2` of
/// the highest available detail level into one atom array (all explicit H
/// are moved to the end). `p_struct[iInchiRec][bMobileH]` are the
/// restored components. Returns the result of the C function and the merged
/// atoms (`OneInput->atom`, `->num_atoms` is their count).
pub fn merge_structure_components(
    p_struct: &[[Vec<super::ichirvr2::StrFromInChI>; 2]; 2],
    one_input: &InpInChI,
) -> (i32, Vec<super::inpdef::InpAtom>) {
    use super::inpdef::{InpAtom, MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS};

    let n_comp = |r: usize, m: usize| one_input.n_num_components[r][m];
    // select the highest detail level
    let (i_inchi_rec, i_mobile_h, num_components) = if n_comp(INCHI_REC, TAUT_NON) != 0 {
        (INCHI_REC, TAUT_NON, n_comp(INCHI_REC, TAUT_NON))
    } else if n_comp(INCHI_REC, TAUT_YES) != 0 {
        (INCHI_REC, TAUT_YES, n_comp(INCHI_REC, TAUT_YES))
    } else if n_comp(INCHI_BAS, TAUT_NON) != 0 {
        (INCHI_BAS, TAUT_NON, n_comp(INCHI_BAS, TAUT_NON))
    } else if n_comp(INCHI_BAS, TAUT_YES) != 0 {
        (INCHI_BAS, TAUT_YES, n_comp(INCHI_BAS, TAUT_YES))
    } else {
        return (0, Vec::new());
    };
    let num_components = num_components as usize;
    let i_altern_h = if i_mobile_h == TAUT_NON && n_comp(i_inchi_rec, TAUT_YES) != 0 {
        Some(TAUT_YES)
    } else {
        None
    };
    // the component structure chosen for the component k
    let pick = |k: usize| {
        let has = |m: usize| p_struct[i_inchi_rec][m].get(k).filter(|c| c.num_atoms != 0);
        has(i_mobile_h).or_else(|| i_altern_h.and_then(has))
    };
    // count number of atoms and removed H
    let mut atom_offs = vec![0usize; num_components + 1];
    let mut del_h_offs = vec![0usize; num_components + 1];
    for k in 0..num_components {
        let (cur_na, cur_nh) = match pick(k) {
            Some(c) if !c.at2.is_empty() && c.b_deleted == 0 => {
                (c.num_atoms as usize, c.num_deleted_h as usize)
            }
            _ => (0, 0),
        };
        atom_offs[k + 1] = atom_offs[k] + cur_na;
        del_h_offs[k + 1] = del_h_offs[k] + cur_nh;
    }
    let tot_just_atoms = atom_offs[num_components];
    // shift all H to the end
    for o in del_h_offs.iter_mut() {
        *o += tot_just_atoms;
    }
    let tot_atoms = del_h_offs[num_components];
    if tot_atoms == 0 {
        return (0, Vec::new());
    }
    // copy
    let mut at = vec![InpAtom::default(); tot_atoms];
    for k in 0..num_components {
        let len = atom_offs[k + 1] - atom_offs[k];
        if len != 0 {
            let c = pick(k).expect("component with atoms");
            at[atom_offs[k]..atom_offs[k] + len].clone_from_slice(&c.at2[..len]);
            let len2 = del_h_offs[k + 1] - del_h_offs[k];
            if len2 != 0 {
                at[del_h_offs[k]..del_h_offs[k] + len2].clone_from_slice(&c.at2[len..len + len2]);
            }
        }
    }
    // update atom numbers
    let mut icomp = 0u16;
    for k in 0..num_components {
        let (cur_at, nxt_at) = (atom_offs[k], atom_offs[k + 1]);
        let (cur_h, nxt_h) = (del_h_offs[k], del_h_offs[k + 1]);
        let len = nxt_at - cur_at;
        let shift_h = cur_h as i64 - len as i64;
        if len == 0 {
            continue;
        }
        icomp += 1;
        let add = |v: u16, base: usize, cond: bool| -> u16 {
            if cond {
                v.wrapping_add(base as u16)
            } else {
                v.wrapping_add(shift_h as u16)
            }
        };
        for a in &mut at[cur_at..nxt_at] {
            a.endpoint = 0;
            a.b_ambiguous_stereo = 0;
            a.at_type = 0;
            a.b_cut_vertex = 0;
            a.b_used_0d_parity = 0;
            a.c_flags = 0;
            a.n_block_system = 0;
            a.n_num_at_in_ring_system = 0;
            a.n_ring_system = 0;
            for j in 0..a.valence as usize {
                let is_atom = (a.neighbor[j] as usize) < len;
                a.neighbor[j] = add(a.neighbor[j], cur_at, is_atom);
            }
            a.orig_at_number = a.orig_at_number.wrapping_add(cur_at as u16);
            a.component = icomp;
            if a.p_parity != 0 {
                for j in 0..MAX_NUM_STEREO_ATOM_NEIGH {
                    // orig_at_num = atom_index + 1, hence <=
                    let is_atom = (a.p_orig_at_num[j] as usize) <= len;
                    a.p_orig_at_num[j] = add(a.p_orig_at_num[j], cur_at, is_atom);
                }
            }
            let mut j = 0;
            while j < MAX_NUM_STEREO_BONDS && a.sb_parity[j] != 0 {
                let is_atom = (a.sn_orig_at_num[j] as usize) <= len;
                a.sn_orig_at_num[j] = add(a.sn_orig_at_num[j], cur_at, is_atom);
                j += 1;
            }
        }
        // fixed-H
        for a in &mut at[cur_h..nxt_h] {
            a.neighbor[0] = a.neighbor[0].wrapping_add(cur_at as u16);
            a.orig_at_number = a.orig_at_number.wrapping_add(shift_h as u16);
        }
    }
    (0, at)
}

/// `AllInchiToStructure`: converts all the components of the input InChI
/// (disconnected/reconnected, Mobile-H/Fixed-H) into structures
/// `p_struct[iInchiRec][bMobileH][k]`. Conversion errors of the single
/// components are stored in `n_error` and counted. Returns the number of
/// the errors or a negative fatal error.
#[allow(clippy::too_many_arguments)]
pub fn all_inchi_to_structure(
    ip_inp: &super::ichiparm::InputParms,
    num_inp: i64,
    srm: &super::ichirvrs::Srm,
    b_has_some_fixed_h: bool,
    p_struct: &mut [[Vec<super::ichirvr2::StrFromInChI>; 2]; 2],
    one_input: &mut InpInChI,
    num_inp_actual: i64,
    polymer: Option<&super::mol2atom::OadPolymer>,
) -> i32 {
    use super::ichiparm::REQ_MODE_BASIC;
    use super::ichirvr2::StrFromInChI;
    use super::runichi3::StructData;

    let mut ip = ip_inp.clone();
    let mut sd = StructData::default();
    let mut num_err = 0;
    for i_inchi_rec in 0..2 {
        for i_mobile_h in 0..2 {
            let cur_num_comp = one_input.n_num_components[i_inchi_rec][i_mobile_h];
            if cur_num_comp == 0 {
                continue;
            }
            p_struct[i_inchi_rec][i_mobile_h] =
                vec![StrFromInChI::default(); cur_num_comp as usize];
            // set conversion mode
            let cur_flag = (if i_mobile_h != 0 { 0 } else { I2A_FLAG_FIXEDH })
                | (if i_inchi_rec != 0 { I2A_FLAG_RECMET } else { 0 });
            if i_mobile_h != 0 {
                ip.n_mode &= !REQ_MODE_BASIC;
            } else {
                ip.n_mode |= REQ_MODE_BASIC;
            }
            // InChI --> structure conversion for all components except duplicated
            let mut at_no_offset = 0;
            for k in 0..cur_num_comp as usize {
                let c = &one_input.p_inp_inchi[i_inchi_rec][i_mobile_h][k];
                if (i_mobile_h == 0 && c.n_number_of_atoms == 0) || c.b_deleted != 0 || c.n_link < 0
                {
                    // no structure out of an unavailable Fixed-H InChI or out
                    // of the one present in the Reconnected layer
                    p_struct[i_inchi_rec][i_mobile_h][k].n_link = c.n_link;
                    p_struct[i_inchi_rec][i_mobile_h][k].b_deleted = c.b_deleted;
                    continue;
                }
                if b_has_some_fixed_h
                    && i_mobile_h != 0
                    && k < one_input.n_num_components[i_inchi_rec][TAUT_NON] as usize
                    && one_input.p_inp_inchi[i_inchi_rec][TAUT_NON][k].n_number_of_atoms != 0
                {
                    continue; // Fixed-H exists: do not process Mobile-H
                }
                let ps = &mut p_struct[i_inchi_rec][i_mobile_h][k];
                ps.srm = *srm;
                ps.i_inchi_rec = i_inchi_rec as i8;
                ps.i_mobile_h = i_mobile_h as i32;
                let ret = inchi2_atom(
                    &ip,
                    &mut sd,
                    num_inp,
                    ps,
                    k,
                    at_no_offset,
                    cur_flag,
                    b_has_some_fixed_h,
                    one_input,
                    num_inp_actual,
                    polymer,
                );
                ps.n_link = one_input.p_inp_inchi[i_inchi_rec][i_mobile_h][k].n_link;
                if ret < 0 {
                    if ret == super::ichierr::CT_USER_QUIT_ERR {
                        return ret;
                    }
                    ps.n_error = ret;
                    num_err += 1; // errors are ignored for now, as in C
                }
                at_no_offset += one_input.p_inp_inchi[i_inchi_rec][i_mobile_h][k].n_number_of_atoms;
            }
        }
    }
    num_err
}

/// `AddProtonAndIsoHBalanceToMobHStruct`: when no Fixed-H layer exists, adds
/// the proton and isotopic exchangeable proton balance of the input InChI to
/// the restored Mobile-H components (the InChI of the changed components is
/// recalculated). Returns a negative error or the result of the last call.
pub fn add_proton_and_iso_h_balance_to_mob_h_struct(
    ip: &super::ichiparm::InputParms,
    sd: &mut super::runichi3::StructData,
    num_inp: i64,
    b_has_some_fixed_h: bool,
    p_struct: &mut [[Vec<super::ichirvr2::StrFromInChI>; 2]; 2],
    one_input: &InpInChI,
) -> i32 {
    use super::ichiread::{NO_VALUE_INT, b_inp_inchi_component_exists};
    use super::ichirvr4::{add_rem_iso_protons_in_restr_struct, add_rem_protons_in_restr_struct};

    if b_has_some_fixed_h {
        return 0;
    }
    // protons removed by the normalization from the original structure
    let mut to_be_removed = [0i32; 2];
    let mut to_be_removed_iso = [[0i16; 3]; 2];
    for i in 0..2 {
        let p = &one_input.n_num_protons[i][TAUT_YES];
        to_be_removed[i] = i32::from(p.n_num_removed_protons);
        to_be_removed_iso[i] = p.n_num_removed_isotopic_h;
    }
    // the protons added/removed by the normalization of the restored
    // structure might have been added by the structure restore
    for i_inchi_rec in 0..2 {
        for k in 0..one_input.n_num_components[i_inchi_rec][TAUT_YES].max(0) as usize {
            if !b_inp_inchi_component_exists(one_input, i_inchi_rec, TAUT_YES, k as i32) {
                continue;
            }
            let mut removed_by_norm_from_revrs = 0;
            let mut removed_by_revrs = 0;
            let n_link = p_struct[i_inchi_rec][TAUT_YES][k].n_link;
            if i_inchi_rec == INCHI_REC || n_link >= 0 {
                let ps = &mut p_struct[i_inchi_rec][TAUT_YES][k];
                let revr = ps.rev_inchi.p_inchi[i_inchi_rec].first();
                let charge_revrs = revr.map_or(NO_VALUE_INT, |r| {
                    r[TAUT_YES]
                        .as_ref()
                        .map_or(NO_VALUE_INT, |c| c.n_total_charge)
                });
                let charge_inchi = one_input.p_inp_inchi[i_inchi_rec][TAUT_YES][k].n_total_charge;
                if let Some(a) = ps.rev_inchi.p_inchi_aux[i_inchi_rec].first() {
                    removed_by_norm_from_revrs = a[TAUT_YES]
                        .as_ref()
                        .map_or(0, |a| i32::from(a.n_num_removed_protons));
                }
                removed_by_revrs = ps.n_num_removed_protons_by_revrs;
                ps.n_charge_revrs = charge_revrs;
                ps.n_charge_inchi = charge_inchi;
            } else {
                let k1 = -(1 + n_link);
                if k1 >= 0 {
                    let k1 = k1 as usize;
                    // as in C: the INCHI_BAS slot of the reconnected structure
                    let (charge_revrs, aux) = {
                        let r = &p_struct[INCHI_REC][TAUT_YES][k1];
                        (
                            r.rev_inchi.p_inchi[INCHI_BAS]
                                .first()
                                .map_or(NO_VALUE_INT, |x| {
                                    x[TAUT_YES]
                                        .as_ref()
                                        .map_or(NO_VALUE_INT, |c| c.n_total_charge)
                                }),
                            r.rev_inchi.p_inchi_aux[INCHI_BAS].first().map(|a| {
                                a[TAUT_YES]
                                    .as_ref()
                                    .map_or(0, |a| i32::from(a.n_num_removed_protons))
                            }),
                        )
                    };
                    let charge_inchi =
                        one_input.p_inp_inchi[INCHI_REC][TAUT_YES][k1].n_total_charge;
                    if let Some(a) = aux {
                        removed_by_norm_from_revrs = a;
                    }
                    // this component is the same as in the reconnected layer
                    removed_by_revrs =
                        p_struct[INCHI_REC][TAUT_YES][k1].n_num_removed_protons_by_revrs;
                    let t = &mut p_struct[i_inchi_rec][TAUT_YES][k1];
                    t.n_charge_revrs = charge_revrs;
                    t.n_charge_inchi = charge_inchi;
                }
            }
            let _ = removed_by_revrs;
            // protons to add = balance in the InChI - protons known to be
            // removed by the normalization from the reconstructed structure
            to_be_removed[i_inchi_rec] -= removed_by_norm_from_revrs;
            p_struct[i_inchi_rec][TAUT_YES][k].n_removed_protons_by_norm_from_revrs =
                removed_by_norm_from_revrs;
        }
    }
    // the fixed-H layer is missing: add the proton balance to the components
    let mut delta_from_disconnected = 0;
    let mut iso_delta_from_disconnected = [0i16; 3];
    let mut ret = 0;
    for i_inchi_rec in [INCHI_REC, INCHI_BAS] {
        let has_rec_mob_h =
            i_inchi_rec == INCHI_BAS && one_input.n_num_components[INCHI_REC][TAUT_YES] != 0;
        if i_inchi_rec == INCHI_BAS {
            // second pass: the common structures have been changed
            to_be_removed[INCHI_BAS] += delta_from_disconnected;
        }
        let (bas, rec) = p_struct.split_at_mut(1);
        let (ps, ps_r) = if i_inchi_rec == INCHI_BAS {
            (&mut bas[0][TAUT_YES], Some(&rec[0][TAUT_YES]))
        } else {
            (&mut rec[0][TAUT_YES], None)
        };
        let n = one_input.n_num_components[i_inchi_rec][TAUT_YES];
        let mut delta = 0;
        ret = add_rem_protons_in_restr_struct(
            ip,
            sd,
            num_inp,
            b_has_some_fixed_h,
            &mut ps[..n.max(0) as usize],
            if has_rec_mob_h {
                ps_r.map(|r| &r[..])
            } else {
                None
            },
            if has_rec_mob_h {
                one_input.n_num_components[INCHI_REC][TAUT_YES].max(0) as usize
            } else {
                0
            },
            &mut to_be_removed[i_inchi_rec],
            (i_inchi_rec == INCHI_REC).then_some(&mut delta),
        );
        if i_inchi_rec == INCHI_REC {
            delta_from_disconnected = delta;
        }
        if ret < 0 {
            return ret;
        }
    }
    // the same for the isotopic exchangeable protons
    for i_inchi_rec in [INCHI_REC, INCHI_BAS] {
        let has_rec_mob_h =
            i_inchi_rec == INCHI_BAS && one_input.n_num_components[INCHI_REC][TAUT_YES] != 0;
        if i_inchi_rec == INCHI_BAS {
            for k in 0..3 {
                to_be_removed_iso[INCHI_BAS][k] += iso_delta_from_disconnected[k];
            }
        }
        let (bas, rec) = p_struct.split_at_mut(1);
        let (ps, ps_r) = if i_inchi_rec == INCHI_BAS {
            (&mut bas[0][TAUT_YES], Some(&rec[0][TAUT_YES]))
        } else {
            (&mut rec[0][TAUT_YES], None)
        };
        let n = one_input.n_num_components[i_inchi_rec][TAUT_YES];
        let mut delta = [0i16; 3];
        ret = add_rem_iso_protons_in_restr_struct(
            ip,
            sd,
            num_inp,
            b_has_some_fixed_h,
            &mut ps[..n.max(0) as usize],
            if has_rec_mob_h {
                ps_r.map(|r| &r[..])
            } else {
                None
            },
            if has_rec_mob_h {
                one_input.n_num_components[INCHI_REC][TAUT_YES].max(0) as usize
            } else {
                0
            },
            &mut to_be_removed_iso[i_inchi_rec],
            (i_inchi_rec == INCHI_REC).then_some(&mut delta),
        );
        if i_inchi_rec == INCHI_REC {
            iso_delta_from_disconnected = delta;
        }
        if ret < 0 {
            return ret;
        }
    }
    ret
}

/// `COMPONENT_REM_PROTONS`: removed protons and exchangeable isotopic H
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ComponentRemProtons {
    pub n_num_removed_protons: i32,
    pub n_num_removed_isotopic_h: [i32; 3],
}

/// `CompareOneOrigInchiToRevInChI`: compares the original InChI of one
/// component (`p_inchi[0]`; `p_inchi[1]` the Mobile-H one for the Fixed-H
/// request) with the InChI created from the restored structure; the
/// difference bits are added to `flags[b_mobile_h]` (or `flags[TAUT_YES]`),
/// the protons removed from the restored Mobile-H layers are returned in
/// `cur_removed`.
pub fn compare_one_orig_inchi_to_rev_inchi(
    ps: &super::ichirvr2::StrFromInChI,
    p_inchi: [Option<&InChI>; 2],
    b_mobile_h: usize,
    cur_removed: &mut ComponentRemProtons,
    flags: &mut [u64; 2],
) -> i32 {
    use super::runichi3::{IS_OKAY, IS_WARNING};

    let ret = ps.rev_inchi.n_ret_val;
    if ret != IS_OKAY && ret != IS_WARNING {
        flags[b_mobile_h] |= diff::STR2INCHI_ERR;
        return 0; // C returns `err`, which is 0 here
    }
    // i0, b0: requested Rec/Disconnected and Mobile/Fixed-H;
    // i, b: what the InChI created out of the restored structure has
    let mut i = ps.i_inchi as usize;
    let i0 = i;
    let mut b = ps.i_mobile_h as usize;
    let b0 = b;
    let rev = &ps.rev_inchi;
    if i == INCHI_REC && rev.num_components[i] == 0 {
        i = INCHI_BAS;
    }
    let slot = |i: usize, k: usize, b: usize| rev.p_inchi[i].get(k).and_then(|c| c[b].as_ref());
    if b == TAUT_NON && !slot(i, 0, b).is_some_and(|c| c.n_number_of_atoms != 0) {
        b = TAUT_YES;
    }
    if ps.b_deleted != 0 && p_inchi[0].is_none_or(|c| c.b_deleted != 0) {
        return 0;
    }
    if !rev.p_inchi[i].is_empty()
        && ((rev.num_components[i] > 1 && !slot(i, 1, b).is_some_and(|c| c.b_deleted != 0))
            || rev.num_components[i] < 1)
    {
        flags[b_mobile_h] |= diff::COMP_NUMBER;
    }
    if b != b0 || b != b_mobile_h || b0 != b_mobile_h || i > i0 {
        // do not print messages about TAUT_YES instead of TAUT_NON
        flags[b_mobile_h] |= diff::COMP_HLAYER;
    }
    if rev.num_components[i] != 0 && !rev.p_inchi[i].is_empty() {
        // '0' in [i][0][b] is the first component
        let mut cmp = 0;
        if b == TAUT_YES
            && slot(i, 0, b).is_some_and(|c| c.b_deleted != 0)
            && p_inchi[0].is_none_or(|c| c.b_deleted != 0)
        {
            // the 1st component is made out of proton(s) and the input
            // component is missing or also a proton
        } else {
            cmp = compare_reversed_inchi3(slot(i, 0, b), p_inchi[0], None, None);
            flags[b_mobile_h] |= cmp;
        }
        let _ = cmp;
        if b == b0 && b == TAUT_NON {
            if slot(i, 0, TAUT_YES).is_some_and(|c| c.b_deleted == 0)
                || p_inchi[1].is_some_and(|c| c.b_deleted == 0)
            {
                // in addition to fixed-H also compare the mobile-H InChI
                flags[TAUT_YES] |=
                    compare_reversed_inchi3(slot(i, 0, TAUT_YES), p_inchi[1], None, None);
            }
            // compare removed H
            let aux0 = rev.p_inchi_aux[i]
                .first()
                .and_then(|c| c[TAUT_YES].as_ref())
                .map_or(0, |a| i32::from(a.n_num_removed_protons));
            if ps.n_num_removed_protons_mobh_inchi != aux0 {
                flags[TAUT_YES] |= diff::MOBH_PROTONS;
            }
        }
        *cur_removed = ComponentRemProtons::default();
        for k in 0..rev.num_components[i].max(0) as usize {
            let deleted = slot(i, k, TAUT_YES).is_some_and(|c| c.b_deleted != 0);
            if k == 0 || deleted {
                // protons of the 1st component; the others only if they are deleted protons
                if let Some(a) = rev.p_inchi_aux[i].get(k).and_then(|c| c[TAUT_YES].as_ref()) {
                    cur_removed.n_num_removed_protons += i32::from(a.n_num_removed_protons);
                    for j in 0..3 {
                        cur_removed.n_num_removed_isotopic_h[j] +=
                            i32::from(a.n_num_removed_isotopic_h[j]);
                    }
                }
            }
        }
    }
    0
}

/// `CompareAllOrigInchiToRevInChI`: compares all components of the input
/// InChI (the Fixed-H layer when requested and present) with the InChI of
/// the restored structures; the bits go to `one_input.compare_inchi_flags[0]`
pub fn compare_all_orig_inchi_to_rev_inchi(
    p_struct: &[[Vec<super::ichirvr2::StrFromInChI>; 2]; 2],
    one_input: &mut InpInChI,
    b_req_non_taut: bool,
) -> i32 {
    one_input.compare_inchi_flags[0] = [0; 2];
    let mut total = ComponentRemProtons::default();
    // do we have reconnected InChI?
    let mut i_inchi_rec = INCHI_REC;
    if one_input.n_num_components[i_inchi_rec][TAUT_YES] == 0
        && one_input.n_num_components[i_inchi_rec][TAUT_NON] == 0
    {
        i_inchi_rec = INCHI_BAS;
    }
    // do we have Mobile or Fixed-H?
    let mut i_mobile_h = TAUT_NON;
    if one_input.n_num_components[i_inchi_rec][TAUT_NON] == 0 || !b_req_non_taut {
        i_mobile_h = TAUT_YES;
    }
    let num_components = one_input.n_num_components[i_inchi_rec][i_mobile_h].max(0) as usize;
    let mut ret = 0;
    for i_component in 0..num_components {
        let mut b_mobile_h = i_mobile_h;
        let mut p_inchi: [Option<&InChI>; 2] = [None, None];
        let cur = &one_input.p_inp_inchi[i_inchi_rec][b_mobile_h][i_component];
        if cur.n_number_of_atoms != 0 && cur.b_deleted == 0 {
            // the requested InChI layer exists
            p_inchi[0] = Some(cur);
            if b_mobile_h == TAUT_NON {
                p_inchi[1] = one_input.p_inp_inchi[i_inchi_rec][TAUT_YES].get(i_component);
            }
        } else if b_mobile_h == TAUT_NON {
            let m = &one_input.p_inp_inchi[i_inchi_rec][TAUT_YES][i_component];
            if m.n_number_of_atoms != 0 && m.b_deleted == 0 {
                // only Mobile-H is available
                b_mobile_h = TAUT_YES;
                p_inchi[0] = Some(m);
            }
        }
        let mut flags = [0u64; 2];
        let mut cur_removed = ComponentRemProtons::default();
        ret = compare_one_orig_inchi_to_rev_inchi(
            &p_struct[i_inchi_rec][b_mobile_h][i_component],
            p_inchi,
            b_mobile_h,
            &mut cur_removed,
            &mut flags,
        );
        if ret < 0 {
            return ret;
        }
        total.n_num_removed_protons += cur_removed.n_num_removed_protons;
        for i in 0..3 {
            total.n_num_removed_isotopic_h[i] += cur_removed.n_num_removed_isotopic_h[i];
        }
        for i in 0..2 {
            one_input.compare_inchi_flags[0][i] |= flags[i];
        }
    }
    if i_mobile_h == TAUT_YES {
        let np = &one_input.n_num_protons[i_inchi_rec][i_mobile_h];
        if np.p_num_protons.is_some() {
            // in Mobile-H case proton balances are split between components
            ret = super::ichiread::RI_ERR_PROGR;
        } else {
            // removed protons in the original and the restored InChI
            if total.n_num_removed_protons != i32::from(np.n_num_removed_protons) {
                one_input.compare_inchi_flags[0][TAUT_YES] |= diff::MOBH_PROTONS;
            }
            let iso = one_input.n_num_protons[i_inchi_rec][TAUT_YES].n_num_removed_isotopic_h;
            for i in 0..3 {
                if total.n_num_removed_isotopic_h[i] != i32::from(iso[i]) {
                    one_input.compare_inchi_flags[0][TAUT_YES] |= diff::MOB_ISO_H;
                }
            }
        }
    }
    ret
}

/// `bRevInchiComponentDeleted`
fn rev_component_deleted(
    ps: &super::ichirvr2::StrFromInChI,
    i_inchi: usize,
    b_mobile_h: usize,
    k: usize,
) -> bool {
    ps.num_atoms != 0
        && (k as i32) < ps.rev_inchi.num_components[i_inchi]
        && ps.rev_inchi.p_inchi[i_inchi]
            .get(k)
            .and_then(|c| c[b_mobile_h].as_ref())
            .is_some_and(|c| c.n_number_of_atoms > 0 && c.b_deleted != 0)
}

/// a fragment of the `INCHI_SORT` arrays of
/// `CompareAllDisconnectedOrigInchiToRevInChI`
#[derive(Default, Clone)]
struct FragSort {
    sort: super::ichi::InChISort,
    n1: i32,
    n2: i32,
}

/// `CompareAllDisconnectedOrigInchiToRevInChI`: compares the fragments of the
/// original InChI that were disconnected from the reconnected (metal)
/// structure with the fragments of the restored reconnected structures; the
/// bits go to `one_input.compare_inchi_flags[1]`. Returns 1 when the numbers
/// of the fragments differ, or a negative error. (The C `qsort` is replaced
/// by a stable sort.)
pub fn compare_all_disconnected_orig_inchi_to_rev_inchi(
    p_struct: &[[Vec<super::ichirvr2::StrFromInChI>; 2]; 2],
    one_input: &mut InpInChI,
    b_has_some_fixed_h: bool,
) -> i32 {
    use super::ichi::comp_inchi_taut2;
    use super::ichiread::{b_inp_inchi_component_deleted, b_inp_inchi_component_exists};

    one_input.compare_inchi_flags[1] = [0; 2];
    // components that are not subject to disconnection
    if one_input.n_num_components[INCHI_REC][TAUT_YES] == 0
        && one_input.n_num_components[INCHI_REC][TAUT_NON] == 0
    {
        return 0;
    }
    let ncomp = |r: usize, m: usize| one_input.n_num_components[r][m];
    let mut rem_d = ComponentRemProtons::default();
    let mut rem_r = ComponentRemProtons::default();
    let mut flags = [0u64; 2];
    let num_components_d = ncomp(INCHI_BAS, TAUT_YES)
        .max(ncomp(INCHI_BAS, TAUT_NON))
        .max(0) as usize;
    let num_components_r = ncomp(INCHI_REC, TAUT_YES)
        .max(ncomp(INCHI_REC, TAUT_NON))
        .max(0) as usize;
    let first_mobile = |r: usize| -> usize {
        if b_has_some_fixed_h {
            usize::from(ncomp(r, TAUT_NON) == 0)
        } else {
            TAUT_YES
        }
    };
    let exists = |o: &InpInChI, r: usize, m: usize, k: usize| {
        b_inp_inchi_component_exists(o, r, m, k as i32)
    };
    // the mobile/fixed-H layer to use for the component k
    let pick_layer = |o: &InpInChI, r: usize, m: usize, k: usize| -> Option<usize> {
        if exists(o, r, m, k) {
            Some(m)
        } else if exists(o, r, TAUT_YES, k) {
            Some(TAUT_YES)
        } else {
            None
        }
    };
    let mut num_non_taut1 = 0;
    let mut num_non_taut2 = 0;

    // ---- original InChI ----
    let mut num_fragments_d = 0;
    let mut have_sep_protons_d = 0;
    let mut have_sep_protons_r = 0;
    {
        let i_inchi = INCHI_BAS;
        let i_mobile_h = first_mobile(i_inchi);
        let np = &one_input.n_num_protons[i_inchi][TAUT_YES];
        rem_d.n_num_removed_protons = i32::from(np.n_num_removed_protons);
        for m in 0..3 {
            rem_d.n_num_removed_isotopic_h[m] = i32::from(np.n_num_removed_isotopic_h[m]);
        }
        for k in 0..num_components_d {
            let Some(b_mobile_h) = pick_layer(one_input, i_inchi, i_mobile_h, k) else {
                continue; // component is missing
            };
            let n_link = one_input.p_inp_inchi[i_inchi][b_mobile_h][k].n_link;
            if n_link < 0 {
                // linked to the identical component in the Reconnected layer
                if let Some(v) = &one_input.n_num_protons[INCHI_REC][TAUT_YES].p_num_protons {
                    have_sep_protons_d += 1;
                    let l = (-(1 + n_link)) as usize;
                    rem_d.n_num_removed_protons += i32::from(v[l][0]);
                    for m in 0..3 {
                        rem_d.n_num_removed_isotopic_h[m] += i32::from(v[l][1 + m]);
                    }
                }
                continue;
            }
            // component in the reconnected layer that was disconnected
            num_non_taut2 += i32::from(b_mobile_h == TAUT_NON);
            if let Some(v) = &one_input.n_num_protons[i_inchi][TAUT_YES].p_num_protons {
                have_sep_protons_d += 1;
                rem_d.n_num_removed_protons += i32::from(v[k][0]);
                for m in 0..3 {
                    rem_d.n_num_removed_isotopic_h[m] += i32::from(v[k][1 + m]);
                }
            }
            num_fragments_d += 1;
        }
    }

    // ---- InChI from the reconstructed reconnected structure ----
    let mut num_fragments_r = 0;
    let i_inchi = INCHI_REC;
    let i_mobile_h = first_mobile(i_inchi);
    let add_aux = |rem: &mut ComponentRemProtons, a: Option<&super::ichi::InChIAux>| {
        if let Some(a) = a {
            rem.n_num_removed_protons += i32::from(a.n_num_removed_protons);
            for m in 0..3 {
                rem.n_num_removed_isotopic_h[m] += i32::from(a.n_num_removed_isotopic_h[m]);
            }
        }
    };
    for k in 0..num_components_r {
        let Some(b_mobile_h) = pick_layer(one_input, i_inchi, i_mobile_h, k) else {
            continue; // deleted proton in the Mobile-H layer
        };
        let ps1 = &p_struct[i_inchi][b_mobile_h][k];
        if one_input.p_inp_inchi[i_inchi][b_mobile_h][k].n_link > 0 {
            // NOT DISCONNECTED: the same component is in the disconnected layer
            if ps1.rev_inchi_component_exists(INCHI_BAS, TAUT_YES, 0) {
                add_aux(
                    &mut rem_r,
                    ps1.rev_inchi.p_inchi_aux[INCHI_BAS]
                        .first()
                        .and_then(|c| c[TAUT_YES].as_ref()),
                );
            }
            have_sep_protons_r += i32::from(ps1.rev_inchi_component_exists(INCHI_BAS, TAUT_NON, 0));
            continue;
        }
        // WAS DISCONNECTED; check its fragments
        let num_fragments = ps1.rev_inchi.num_components[INCHI_BAS].max(0) as usize;
        let if_mobile_h = if b_has_some_fixed_h {
            TAUT_NON
        } else {
            TAUT_YES
        };
        for ifk in 0..num_fragments {
            let bf = if ps1.rev_inchi_component_exists(INCHI_BAS, if_mobile_h, ifk) {
                if_mobile_h
            } else if ps1.rev_inchi_component_exists(INCHI_BAS, TAUT_YES, ifk) {
                TAUT_YES
            } else {
                continue; // fragment does not exist
            };
            num_non_taut1 += i32::from(bf == TAUT_NON);
            have_sep_protons_r += i32::from(bf == TAUT_NON);
            add_aux(
                &mut rem_r,
                ps1.rev_inchi.p_inchi_aux[INCHI_BAS]
                    .get(ifk)
                    .and_then(|c| c[TAUT_YES].as_ref()),
            );
            num_fragments_r += 1;
        }
    }
    // the last reconstructed component may contain separate protons added
    // by the reconstruction
    if num_components_r > 0 {
        let k = num_components_r - 1;
        if let Some(ps1) = p_struct[i_inchi][i_mobile_h].get(k)
            && i_mobile_h == TAUT_YES
            && !b_has_some_fixed_h
            && b_inp_inchi_component_deleted(one_input, i_inchi, i_mobile_h, k as i32)
        {
            for ifk in 0..ps1.rev_inchi.num_components[INCHI_BAS].max(0) as usize {
                if !rev_component_deleted(ps1, INCHI_BAS, TAUT_YES, ifk) {
                    continue;
                }
                add_aux(
                    &mut rem_r,
                    ps1.rev_inchi.p_inchi_aux[INCHI_BAS]
                        .get(ifk)
                        .and_then(|c| c[TAUT_YES].as_ref()),
                );
            }
        }
    }
    let num_fragments_dr = num_fragments_d.max(num_fragments_r);
    if num_fragments_dr == 0 {
        return 0; // no component was disconnected
    }
    if num_fragments_d != num_fragments_r {
        for i in 0..2 {
            if one_input.n_num_components[INCHI_BAS][i] != 0 {
                one_input.compare_inchi_flags[1][i] |= diff::PROBLEM;
            }
        }
        return 1; // severe error
    }
    let ndr = num_fragments_dr as usize;
    let mut sort1 = vec![FragSort::default(); ndr];
    let mut sort2 = vec![FragSort::default(); ndr];

    // the original InChI of the fragments
    {
        let i_inchi = INCHI_BAS;
        let i_mobile_h = first_mobile(i_inchi);
        let mut n = 0;
        for k in 0..num_components_d {
            let Some(b_mobile_h) = pick_layer(one_input, i_inchi, i_mobile_h, k) else {
                continue;
            };
            if one_input.p_inp_inchi[i_inchi][b_mobile_h][k].n_link < 0 {
                continue; // same as reconnected
            }
            let f = &mut sort2[n];
            f.sort.p_inchi[b_mobile_h] =
                Some(one_input.p_inp_inchi[i_inchi][b_mobile_h][k].clone());
            if b_mobile_h == TAUT_NON
                && (exists(one_input, i_inchi, TAUT_YES, k)
                    || b_inp_inchi_component_deleted(one_input, i_inchi, TAUT_YES, k as i32))
            {
                f.sort.p_inchi[TAUT_YES] =
                    Some(one_input.p_inp_inchi[i_inchi][TAUT_YES][k].clone());
            }
            // the last sort key is the number of removed protons
            f.sort.ord_number = one_input.n_num_protons[i_inchi][TAUT_YES]
                .p_num_protons
                .as_ref()
                .map_or(0, |v| v[k][0]);
            f.n1 = k as i32;
            f.n2 = -1;
            n += 1;
        }
    }
    // the fragments from the reconstructed structure
    {
        let i_inchi = INCHI_REC;
        let i_mobile_h = first_mobile(i_inchi);
        let mut n = 0;
        for k in 0..num_components_r {
            let Some(b_mobile_h) = pick_layer(one_input, i_inchi, i_mobile_h, k) else {
                continue;
            };
            if one_input.p_inp_inchi[i_inchi][b_mobile_h][k].n_link > 0 {
                continue; // same as disconnected
            }
            let ps1 = &p_struct[i_inchi][b_mobile_h][k];
            let if_mobile_h = if b_has_some_fixed_h {
                TAUT_NON
            } else {
                TAUT_YES
            };
            for i in 0..ps1.rev_inchi.num_components[INCHI_BAS].max(0) as usize {
                let bf = if ps1.rev_inchi_component_exists(INCHI_BAS, if_mobile_h, i) {
                    if_mobile_h
                } else if ps1.rev_inchi_component_exists(INCHI_BAS, TAUT_YES, i) {
                    TAUT_YES
                } else {
                    continue;
                };
                let slot = |m: usize| ps1.rev_inchi.p_inchi[INCHI_BAS][i][m].clone();
                let f = &mut sort1[n];
                f.sort.p_inchi[bf] = slot(bf);
                if bf == TAUT_NON {
                    f.sort.p_inchi[TAUT_YES] = slot(TAUT_YES);
                    // remove Fixed-H InChI identical to Mobile-H, as for InpInChI
                    if compare_reversed_inchi(
                        f.sort.p_inchi[bf].as_ref(),
                        f.sort.p_inchi[TAUT_YES].as_ref(),
                        None,
                        None,
                    ) == 0
                    {
                        f.sort.p_inchi[bf] = None;
                    } else {
                        f.sort.ord_number = ps1.rev_inchi.p_inchi_aux[INCHI_BAS][i][TAUT_YES]
                            .as_ref()
                            .map_or(0, |a| a.n_num_removed_protons);
                    }
                }
                f.n1 = k as i32;
                f.n2 = i as i32;
                n += 1;
            }
        }
    }
    // sort the fragments before comparing
    let by = |a: &FragSort, b: &FragSort| comp_inchi_taut2(&a.sort, &b.sort).cmp(&0);
    sort1.sort_by(by);
    sort2.sort_by(by);
    // compare fragments -- components present in the disconnected layer only
    for i_component in 0..ndr {
        let a = &sort1[i_component].sort;
        let b = &sort2[i_component].sort;
        compare_two_pairs_of_inchi(
            [a.p_inchi[0].as_ref(), a.p_inchi[1].as_ref()],
            [b.p_inchi[0].as_ref(), b.p_inchi[1].as_ref()],
            &mut flags,
        );
    }
    if b_has_some_fixed_h {
        if have_sep_protons_d != 0 || have_sep_protons_r != 0 {
            // for each component, compare the number of removed protons
            'frag: for i_component in 0..ndr {
                let mut iso1 = [0i16; 3];
                let mut iso2 = [0i16; 3];
                if sort1[i_component].sort.ord_number != sort2[i_component].sort.ord_number {
                    flags[TAUT_YES] |= diff::MOBH_PROTONS;
                }
                // also compare the removed isotopic H
                let k = sort2[i_component].n1 as usize;
                if let Some(v) = &one_input.n_num_protons[INCHI_BAS][TAUT_YES].p_num_protons {
                    iso2.copy_from_slice(&v[k][1..4]);
                }
                // removed protons info of the fragments of the restored structure
                let k = sort1[i_component].n1 as usize;
                let i = sort1[i_component].n2 as usize;
                let i_inchi = INCHI_REC;
                let i_mobile_h = first_mobile(i_inchi);
                let Some(b_mobile_h) = pick_layer(one_input, i_inchi, i_mobile_h, k) else {
                    if iso1 != iso2 {
                        flags[TAUT_YES] |= diff::REM_ISO_H;
                    }
                    continue 'frag;
                };
                if one_input.p_inp_inchi[i_inchi][b_mobile_h][k].n_link != 0 {
                    continue;
                }
                let ps1 = &p_struct[i_inchi][b_mobile_h][k];
                let num_fragments = ps1.rev_inchi.num_components[INCHI_BAS].max(0) as usize;
                let if_mobile_h = if b_has_some_fixed_h {
                    TAUT_NON
                } else {
                    TAUT_YES
                };
                if i < num_fragments
                    && (ps1.rev_inchi_component_exists(INCHI_BAS, if_mobile_h, i)
                        || ps1.rev_inchi_component_exists(INCHI_BAS, TAUT_YES, i))
                    && let Some(a) = ps1.rev_inchi.p_inchi_aux[INCHI_BAS][i][TAUT_YES].as_ref()
                {
                    iso1 = a.n_num_removed_isotopic_h;
                }
                if iso1 != iso2 {
                    flags[TAUT_YES] |= diff::REM_ISO_H;
                }
            }
        }
    } else {
        // compare the totals of removed protons and isotopic H
        let np = &one_input.n_num_protons[INCHI_BAS][TAUT_YES];
        if i32::from(np.n_num_removed_protons) != rem_r.n_num_removed_protons {
            flags[TAUT_YES] |= diff::MOBH_PROTONS;
        }
        if (0..3)
            .any(|m| i32::from(np.n_num_removed_isotopic_h[m]) != rem_r.n_num_removed_isotopic_h[m])
        {
            flags[TAUT_YES] |= diff::REM_ISO_H;
        }
    }
    if (num_non_taut1 == 0) != (num_non_taut2 == 0) {
        if num_non_taut1 != 0 {
            // the restored structure has Fixed-H while the original has not
            flags[TAUT_YES] |= diff::WRONG_TAUT;
        } else {
            // the original InChI has Fixed-H while the restored one has not
            flags[TAUT_YES] |= diff::NO_TAUT;
        }
    }
    for i in 0..2 {
        one_input.compare_inchi_flags[1][i] |= flags[i];
    }
    // (the totals of the removed protons are compared after the flags have
    // been stored, so the bits set here are not reported, as in C)
    let _ = rem_d;
    0
}
