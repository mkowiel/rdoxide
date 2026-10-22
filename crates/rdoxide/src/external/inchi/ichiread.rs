//! The InChI string reader of `INCHI_BASE/src/ichiread.c` (`InChILine2Data`
//! and its segment parsers). Done so far: the error names and states, the
//! `/q` and `/p` layers, delimited number lists, the proton placeholder
//! `INChI`, and the segment reader of the input stream.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! The C parsers walk `const char *` pointers over zero-terminated strings;
//! here a string is a byte slice that reads as zero past its end and a
//! position is an index. Only string input streams are supported.

use super::ichi::InChI;
use super::ichi_io::{INCHI_IOS_TYPE_STRING, InchiIostream};
use super::util::strtol10;

pub const RI_ERR_ALLOC: i32 = -1;
pub const RI_ERR_SYNTAX: i32 = -2;
pub const RI_ERR_PROGR: i32 = -3;
pub const RI_ERR_EOL: i32 = -4;
pub const RI_ERR_EOF: i32 = 0;
pub const RI_ERR_MISMATCH: i32 = -9;
/// `NO_VALUE_INT`
pub const NO_VALUE_INT: i32 = 9999;
pub const MAX_ATOMS: i32 = 32766;

pub const TAUT_NON: usize = 0;
pub const TAUT_YES: usize = 1;

// `INCHI_STATE`
pub const IST_MOBILE_H_FORMULA: i32 = 0;
pub const IST_MOBILE_H_CONNECTIONS: i32 = 1;
pub const IST_MOBILE_H: i32 = 2;
pub const IST_MOBILE_H_CHARGE: i32 = 3;
pub const IST_MOBILE_H_PROTONS: i32 = 4;
pub const IST_MOBILE_H_SP2: i32 = 5;
pub const IST_MOBILE_H_SP3: i32 = 6;
pub const IST_MOBILE_H_SP3_M: i32 = 7;
pub const IST_MOBILE_H_SP3_S: i32 = 8;
pub const IST_MOBILE_H_ISO_LAYER_FORK: i32 = 9;
pub const IST_MOBILE_H_ISO_ATOMS: i32 = 10;
pub const IST_MOBILE_H_ISO_EXCH_H: i32 = 11;
pub const IST_MOBILE_H_ISO_SP2: i32 = 12;
pub const IST_MOBILE_H_ISO_SP3: i32 = 13;
pub const IST_MOBILE_H_ISO_SP3_M: i32 = 14;
pub const IST_MOBILE_H_ISO_SP3_S: i32 = 15;
pub const IST_FIXED_H_LAYER_FORK: i32 = 16;
pub const IST_FIXED_H_FORMULA: i32 = 17;
pub const IST_FIXED_H: i32 = 18;
pub const IST_FIXED_H_CHARGE: i32 = 19;
pub const IST_FIXED_H_SP2: i32 = 20;
pub const IST_FIXED_H_SP3: i32 = 21;
pub const IST_FIXED_H_SP3_M: i32 = 22;
pub const IST_FIXED_H_SP3_S: i32 = 23;
pub const IST_FIXED_H_PERMUTATION: i32 = 24;
pub const IST_FIXED_H_ISO_LAYER_FORK: i32 = 25;
pub const IST_FIXED_H_ISO_ATOMS: i32 = 26;
pub const IST_FIXED_H_ISO_LAYER: i32 = 27;
pub const IST_FIXED_H_ISO_SP2: i32 = 28;
pub const IST_FIXED_H_ISO_SP3: i32 = 29;
pub const IST_FIXED_H_ISO_SP3_M: i32 = 30;
pub const IST_FIXED_H_ISO_SP3_S: i32 = 31;
pub const IST_FIXED_H_ISO_PERMUTATION: i32 = 32;
pub const IST_RECONNECTED_LAYER_FORK: i32 = 33;
pub const IST_RECONNECTED_FORMULA: i32 = 34;
pub const IST_MATERIAL_BALANCE_ERROR: i32 = 35;
pub const IST_MOBILE_H_POLYMER: i32 = 36;
pub const IST_END: i32 = -1;
pub const IST_HAPPENED_IN_RECMET: i32 = 100;

/// `irErrMsg`
const IR_ERR_MSG: [(i32, &str); 38] = [
    (IST_MOBILE_H_FORMULA, "MOBILE_H_FORMULA"),
    (IST_MOBILE_H_CONNECTIONS, "MOBILE_H_CONNECTIONS"),
    (IST_MOBILE_H, "MOBILE_H"),
    (IST_MOBILE_H_CHARGE, "MOBILE_H_CHARGE"),
    (IST_MOBILE_H_PROTONS, "MOBILE_H_PROTONS"),
    (IST_MOBILE_H_SP2, "MOBILE_H_SP2"),
    (IST_MOBILE_H_SP3, "MOBILE_H_SP3"),
    (IST_MOBILE_H_SP3_M, "MOBILE_H_SP3_/m"),
    (IST_MOBILE_H_SP3_S, "MOBILE_H_SP3_/s"),
    (IST_MOBILE_H_ISO_LAYER_FORK, "MOBILE_H_ISO_LAYER_FORK"),
    (IST_MOBILE_H_ISO_ATOMS, "MOBILE_H_ISO_ATOMS"),
    (IST_MOBILE_H_ISO_EXCH_H, "MOBILE_H_ISO_EXCH_H"),
    (IST_MOBILE_H_ISO_SP2, "MOBILE_H_ISO_SP2"),
    (IST_MOBILE_H_ISO_SP3, "MOBILE_H_ISO_SP3"),
    (IST_MOBILE_H_ISO_SP3_M, "MOBILE_H_ISO_SP3_/m"),
    (IST_MOBILE_H_ISO_SP3_S, "MOBILE_H_ISO_SP3_/s"),
    (IST_FIXED_H_LAYER_FORK, "FIXED_H_LAYER_FORK"),
    (IST_FIXED_H_FORMULA, "FIXED_H_FORMULA"),
    (IST_FIXED_H, "FIXED_H"),
    (IST_FIXED_H_CHARGE, "FIXED_H_CHARGE"),
    (IST_FIXED_H_SP2, "FIXED_H_SP2"),
    (IST_FIXED_H_SP3, "FIXED_H_SP3"),
    (IST_FIXED_H_SP3_M, "FIXED_H_SP3_/m"),
    (IST_FIXED_H_SP3_S, "FIXED_H_SP3_/s"),
    (IST_FIXED_H_PERMUTATION, "FIXED_H_PERMUTATION"),
    (IST_FIXED_H_ISO_LAYER_FORK, "FIXED_H_ISO_LAYER_FORK"),
    (IST_FIXED_H_ISO_ATOMS, "FIXED_H_ISO_ATOMS"),
    (IST_FIXED_H_ISO_LAYER, "FIXED_H_ISO_LAYER"),
    (IST_FIXED_H_ISO_SP2, "FIXED_H_ISO_SP2"),
    (IST_FIXED_H_ISO_SP3, "FIXED_H_ISO_SP3"),
    (IST_FIXED_H_ISO_SP3_M, "FIXED_H_ISO_SP3_m"),
    (IST_FIXED_H_ISO_SP3_S, "FIXED_H_ISO_SP3_s"),
    (IST_FIXED_H_ISO_PERMUTATION, "FIXED_H_ISO_PERMUTATION"),
    (IST_RECONNECTED_LAYER_FORK, "RECONNECTED_LAYER_FORK"),
    (IST_RECONNECTED_FORMULA, "RECONNECTED_FORMULA"),
    (IST_MATERIAL_BALANCE_ERROR, "MATERIAL_BALANCE"),
    (IST_MOBILE_H_POLYMER, "POLYMER_LAYER"),
    (IST_END, "Unknown Error"),
];

/// `getInchiStateReadErr`: the name of a reading state
pub fn get_inchi_state_read_err(stat: i32) -> String {
    let mut stat = stat;
    let mut rec_met = false;
    if stat >= IST_HAPPENED_IN_RECMET {
        rec_met = true;
        stat -= IST_HAPPENED_IN_RECMET;
    }
    let mut i = 0;
    while 0 <= IR_ERR_MSG[i].0 && stat != IR_ERR_MSG[i].0 {
        i += 1;
    }
    format!(
        "{}{}",
        IR_ERR_MSG[i].1,
        if rec_met { ", Reconnected layer" } else { "" }
    )
}

/// `getInchiErrName`
pub fn get_inchi_err_name(n_err: i32) -> &'static str {
    match n_err {
        RI_ERR_ALLOC => "Allocation failed",
        RI_ERR_PROGR => "Program error",
        RI_ERR_SYNTAX => "Syntax error",
        RI_ERR_EOL => "End of line",
        _ => "Unknown error",
    }
}

/// a zero-terminated C string as a slice: reads past the end give 0
fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// `strchr(&s[from..], c)` for `c != 0`: the position or `None`
fn strchr(s: &[u8], from: usize, c: u8) -> Option<usize> {
    let mut i = from;
    while at(s, i) != 0 {
        if s[i] == c {
            return Some(i);
        }
        i += 1;
    }
    None
}

/// `strlen`
fn strlen(s: &[u8], from: usize) -> usize {
    let mut i = from;
    while at(s, i) != 0 {
        i += 1;
    }
    i - from
}

/// `inchi_strtol(&s[from..], &end, 10)`: the value as C `long` and the end
/// position
fn strtol(s: &[u8], from: usize) -> (i64, usize) {
    let tail = if from <= s.len() { &s[from..] } else { &[][..] };
    // a zero byte ends the string
    let n = tail.iter().position(|&c| c == 0).unwrap_or(tail.len());
    let (v, used) = strtol10(&tail[..n]);
    (v, from + used)
}

/// `ParseSegmentReadDelimitedNumbers`: reads the numbers of `s[from + 1..]`
/// delimited with `c_delim` into `numlist` until `c_stop` or `p_end`. Returns
/// the end position (`None` for the C NULL) and `ret` (1 on success). (The C
/// code loops forever when the end is reached without `c_stop`; here the
/// position is returned.)
pub fn parse_segment_read_delimited_numbers(
    s: &[u8],
    from: usize,
    p_end: usize,
    numlist: &mut Vec<i32>,
    c_delim: u8,
    c_stop: u8,
) -> (Option<usize>, i32) {
    let p_start = from + 1;
    let mut p = p_start;
    if at(s, p_start) != 0 {
        p = p_start;
        while p < p_end && at(s, p) != c_stop {
            let (v, q) = strtol(s, p);
            p = q;
            let num = i32::from(v as u16);
            if !(0..=MAX_ATOMS).contains(&num) {
                return (Some(p), RI_ERR_SYNTAX);
            }
            numlist.push(num);
            if at(s, p) == c_delim {
                p += 1;
            } else if at(s, p) == c_stop {
                return (Some(p), 1);
            } else {
                return (None, -1);
            }
        }
        if at(s, p) == c_stop {
            return (Some(p), 1);
        }
    }
    (Some(p), 1)
}

/// the part of `REM_PROTONS` that the parsers set
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RemProtons {
    pub n_num_removed_protons: i16,
    pub n_num_removed_isotopic_h: [i16; 3],
    pub p_num_protons: Option<Vec<[i16; 4]>>,
}

fn is_digit(c: u8) -> bool {
    c.is_ascii_digit()
}

/// `ParseSegmentProtons`: the `/p` layer (`s` starts at the `p`)
pub fn parse_segment_protons(s: &[u8], b_mobile_h: usize, n_num_protons: &mut [RemProtons]) -> i32 {
    if at(s, 0) != b'p' {
        return 0;
    }
    let p_start = 1;
    let p_end = strchr(s, p_start, b';').unwrap_or(p_start + strlen(s, p_start));
    let (val, q);
    if at(s, p_start) == b'+' && is_digit(at(s, p_start + 1)) {
        let (v, e) = strtol(s, p_start + 1);
        val = v as i32;
        q = e;
    } else if at(s, p_start) == b'-' && is_digit(at(s, p_start + 1)) {
        let (v, e) = strtol(s, p_start + 1);
        val = -(v as i32);
        q = e;
    } else {
        return RI_ERR_SYNTAX;
    }
    if val == 0 {
        return RI_ERR_SYNTAX;
    }
    n_num_protons[b_mobile_h].n_num_removed_protons = val as i16;
    if at(s, p_end) != 0 || q != p_end {
        return RI_ERR_SYNTAX;
    }
    1
}

/// `ParseSegmentCharge`: the `/q` layer (`s` starts at the `q`);
/// `p_inp_inchi[TAUT_NON/TAUT_YES]` hold the components
pub fn parse_segment_charge(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    num_components: [i32; 2],
) -> i32 {
    if at(s, 0) != b'q' {
        return 0;
    }
    const MULT_TYPE: &[u8] = b"mnMNe";
    let mut p_start = 1;
    let mut i_component: i32 = 0;
    let n_num_components = num_components[b_mobile_h];
    if at(s, p_start) == 0 && b_mobile_h == TAUT_NON {
        for i in 0..n_num_components as usize {
            p_inp_inchi[b_mobile_h][i].n_total_charge = NO_VALUE_INT;
        }
        return n_num_components + 1;
    }
    loop {
        let p_end = strchr(s, p_start, b';').unwrap_or(p_start + strlen(s, p_start));
        let (mut val, mut q, if_cnd);
        if (strtol(s, p_start).0 as i32) > 0 {
            let (v, e) = strtol(s, p_start);
            val = v as i32;
            q = e;
            if_cnd = is_digit(at(s, p_start));
        } else {
            val = 1;
            q = p_start;
            if_cnd = true;
        }
        let mpy_component: i32;
        let qc = at(s, q);
        // `strchr(mult_type, *q)` also finds the terminator for *q == 0
        let in_mult = qc == 0 || MULT_TYPE.contains(&qc);
        let mut abbreviation = false;
        let mut p;
        if if_cnd && in_mult && q + 1 == p_end {
            match b_mobile_h {
                TAUT_YES => return RI_ERR_SYNTAX,
                TAUT_NON => {
                    if qc != b'm'
                        || i_component + val > n_num_components
                        || i_component + val > num_components[TAUT_YES]
                    {
                        return RI_ERR_SYNTAX;
                    }
                    for i in 0..val {
                        let ic = (i_component + i) as usize;
                        let c = p_inp_inchi[TAUT_YES][ic].n_total_charge;
                        p_inp_inchi[TAUT_NON][ic].n_total_charge =
                            if c != 0 { c } else { NO_VALUE_INT };
                    }
                    mpy_component = val;
                    abbreviation = true;
                    p = p_start;
                }
                _ => return RI_ERR_SYNTAX,
            }
        } else if let Some(star) = strchr(s, p_start, b'*').filter(|&x| x < p_end) {
            let (v, e) = strtol(s, p_start);
            let m = v as i32;
            if !(0..=MAX_ATOMS).contains(&m) {
                return RI_ERR_SYNTAX;
            }
            if star != e {
                return RI_ERR_SYNTAX;
            }
            mpy_component = m;
            p = star + 1;
        } else {
            mpy_component = 1;
            p = p_start;
        }
        if !abbreviation {
            if mpy_component + i_component > n_num_components || mpy_component <= 0 {
                return RI_ERR_SYNTAX;
            }
            p_start = p;
            if p_start < p_end {
                if at(s, p_start) == b'+' && is_digit(at(s, p_start + 1)) {
                    let (v, e) = strtol(s, p_start + 1);
                    val = v as i32;
                    p_start = e;
                } else if at(s, p_start) == b'-' && is_digit(at(s, p_start + 1)) {
                    let (v, e) = strtol(s, p_start + 1);
                    val = -(v as i32);
                    p_start = e;
                } else {
                    return RI_ERR_SYNTAX;
                }
                if !(-256..=256).contains(&val) {
                    return RI_ERR_SYNTAX;
                }
                if val == 0 {
                    if p_start != p_end {
                        return RI_ERR_SYNTAX;
                    }
                    if b_mobile_h == TAUT_NON {
                        val = NO_VALUE_INT;
                    }
                }
            } else {
                val = NO_VALUE_INT;
            }
            for i in 0..mpy_component {
                p_inp_inchi[b_mobile_h][(i_component + i) as usize].n_total_charge = val;
            }
        }
        let _ = (&mut q, &mut p);
        i_component += mpy_component;
        if at(s, p_end) != 0 {
            p_start = p_end + 1;
            continue;
        }
        break;
    }
    if n_num_components != i_component {
        return RI_ERR_SYNTAX;
    }
    i_component + 1
}

/// `nFillOutProtonMobileH`: makes `inchi` the placeholder `INChI` of the
/// removed protons (a single "H" atom)
pub fn n_fill_out_proton_mobile_h(inchi: &mut InChI) -> i32 {
    inchi.b_deleted = 1;
    inchi.sz_hill_formula = Some("H".to_string());
    inchi.n_number_of_atoms = 1;
    inchi.n_atom = vec![1, 0];
    inchi.n_total_charge = 1;
    inchi.n_conn_table = vec![1, 0];
    inchi.len_conn_table = 1;
    inchi.n_tautomer = vec![0, 0];
    inchi.n_num_h = vec![0, 0];
    inchi.len_tautomer = 1;
    0
}

/// `nProtonCopyIsotopicInfo`
pub fn n_proton_copy_isotopic_info(to: &mut InChI, from: &InChI) -> i32 {
    if from.n_number_of_isotopic_atoms != 0 {
        if to.n_number_of_isotopic_atoms != 0
            && from.n_number_of_isotopic_atoms > to.n_number_of_isotopic_atoms
        {
            to.isotopic_atom.clear();
            to.n_number_of_isotopic_atoms = 0;
        }
        to.isotopic_atom = from.isotopic_atom[..from.n_number_of_isotopic_atoms as usize].to_vec();
        to.n_number_of_isotopic_atoms = from.n_number_of_isotopic_atoms;
    } else {
        to.isotopic_atom.clear();
        to.n_number_of_isotopic_atoms = 0;
    }
    0
}

/// `insertions_sort_AT_NUMB`: sorts ascending, returns the number of moves
pub fn insertions_sort_at_numb(base: &mut [u16], num: usize) -> i32 {
    let mut num_trans = 0;
    for k in 1..num {
        let tmp = base[k];
        let mut j = k;
        while j > 0 && base[j - 1] > tmp {
            base[j] = base[j - 1];
            num_trans += 1;
            j -= 1;
        }
        base[j] = tmp;
    }
    num_trans
}

/// `SEGM_LINE`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SegmLine {
    pub str: Vec<u8>,
    pub len: usize,
    pub c: i32,
}

/// `getInChIChar` for a string stream
pub fn get_inchi_char(inp: &mut InchiIostream) -> i32 {
    if inp.type_ == INCHI_IOS_TYPE_STRING {
        if inp.s.n_ptr < inp.s.n_used_length {
            let c = inp
                .s
                .p_str
                .as_ref()
                .and_then(|p| p.get(inp.s.n_ptr as usize).copied())
                .unwrap_or(0);
            inp.s.n_ptr += 1;
            return i32::from(c as i8);
        }
        return RI_ERR_EOF;
    }
    RI_ERR_EOF
}

fn put(line: &mut SegmLine, i: usize, c: u8) {
    if line.str.len() <= i {
        line.str.resize(i + 1, 0);
    }
    line.str[i] = c;
}

/// `AddInChIChar`
pub fn add_inchi_char(inp: &mut InchiIostream, line: &mut SegmLine, token: &[u8]) -> i32 {
    let mut c = get_inchi_char(inp);
    if c < 0 {
        put(line, line.len, 0);
        c = RI_ERR_SYNTAX;
        return c;
    }
    if c != 0 && token.contains(&(c as u8)) {
        put(line, line.len, 0);
        return -(c + 2);
    } else if c == 0 && line.len == 0 {
        put(line, line.len, 0);
    } else {
        put(line, line.len, c as u8);
        line.len += 1;
    }
    c
}

/// `nGetInChISegment`
pub fn n_get_inchi_segment(inp: &mut InchiIostream, line: &mut SegmLine, token: &[u8]) -> i32 {
    line.len = 0;
    let mut c;
    loop {
        c = add_inchi_char(inp, line, token);
        if c <= 0 {
            break;
        }
    }
    if c < -2 {
        c = -(c + 2);
    }
    line.c = c;
    c
}
