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

use super::ichi::{
    INCHI_FLAG_RAC_STEREO, INCHI_FLAG_REL_STEREO, InChI, InChIIsotopicAtom, InChIStereo,
};
use super::ichi_io::{INCHI_IOS_TYPE_STRING, InchiIostream};
use super::inpdef::MAXVAL;
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

use super::util::{ERR_ELEM, get_periodic_table_number};

/// `isupper`/`islower`/`isdigit` of a position (reads as 0 past the end)
fn upper(s: &[u8], i: usize) -> bool {
    at(s, i).is_ascii_uppercase()
}
fn lower(s: &[u8], i: usize) -> bool {
    at(s, i).is_ascii_lowercase()
}

/// the elements and counts of a Hill formula string, as parsed by
/// `ParseSegmentFormula` pass 2: `(element symbol, count)`; `Err` is the
/// return code
fn formula_items(f: &[u8]) -> Result<Vec<(String, i32)>, i32> {
    let mut items = Vec::new();
    let mut e = 0;
    while at(f, e) != 0 {
        if !upper(f, e) {
            return Err(RI_ERR_SYNTAX);
        }
        let mut el = String::new();
        el.push(at(f, e) as char);
        e += 1;
        if lower(f, e) {
            el.push(at(f, e) as char);
            e += 1;
        }
        let mpy_atom;
        if is_digit(at(f, e)) {
            let (v, q) = strtol(f, e);
            mpy_atom = v as i32;
            if !(0..=MAX_ATOMS).contains(&mpy_atom) {
                return Err(RI_ERR_SYNTAX);
            }
            e = q;
        } else {
            mpy_atom = 1;
        }
        if mpy_atom == 0 {
            return Err(RI_ERR_SYNTAX);
        }
        items.push((el, mpy_atom));
    }
    Ok(items)
}

/// `ParseSegmentFormula`: the formula layer; `*na_total` is the reading
/// state on entry and the total number of atoms except H on return
pub fn parse_segment_formula(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    pn_num_components: &mut [i32; 2],
    na_total: &mut i32,
) -> i32 {
    let alt = 1 - b_mobile_h;
    let state = *na_total;
    *na_total = 0;
    let mut n_num_components: i32 = 0;

    // Pass 1. Count components
    let mut p_start = 0usize;
    let mut p_end;
    loop {
        p_end = strchr(s, p_start, b'.').unwrap_or(p_start + strlen(s, p_start));
        let mut p = p_start;
        let mpy_component;
        if is_digit(at(s, p)) {
            let (v, q) = strtol(s, p);
            mpy_component = v as i32;
            p = q;
        } else {
            mpy_component = 1;
        }
        if mpy_component == 0 {
            break;
        }
        if !upper(s, p) {
            break;
        }
        if p_end == p {
            break;
        }
        n_num_components = n_num_components.wrapping_add(mpy_component);
        if at(s, p_end) != 0 {
            p_start = p_end + 1;
        } else {
            break;
        }
    }
    pn_num_components[b_mobile_h] = n_num_components;
    if n_num_components > MAX_ATOMS {
        pn_num_components[b_mobile_h] = 0;
        return RI_ERR_SYNTAX;
    }
    if n_num_components == 0 {
        let low_case = at(s, p_start) != 0 && lower(s, p_start);
        if at(s, p_start) == 0 || low_case {
            if low_case
                && state == IST_MOBILE_H_FORMULA
                && at(s, p_start) != b'q'
                && at(s, p_start) != b'p'
            {
                return RI_ERR_SYNTAX;
            }
            let n_alt = pn_num_components[alt];
            if b_mobile_h == TAUT_NON && n_alt > 0 {
                let mut v = vec![InChI::default(); n_alt as usize];
                pn_num_components[b_mobile_h] = n_alt;
                for (i, c) in v.iter_mut().enumerate() {
                    let a = &p_inp_inchi[alt][i];
                    c.n_number_of_atoms = a.n_number_of_atoms;
                    let len = (a.n_number_of_atoms + 1) as usize;
                    let mut na: Vec<u8> = (0..len)
                        .map(|k| a.n_atom.get(k).copied().unwrap_or(0))
                        .collect();
                    na.push(0);
                    c.n_atom = na;
                    let f = a.sz_hill_formula.clone().unwrap_or_default();
                    c.sz_hill_formula = Some(f);
                }
                p_inp_inchi[b_mobile_h] = v;
            } else if b_mobile_h == TAUT_YES {
                pn_num_components[b_mobile_h] = 1;
                let mut c = InChI::default();
                let ret = n_fill_out_proton_mobile_h(&mut c);
                p_inp_inchi[b_mobile_h] = vec![c];
                if ret < 0 {
                    return ret;
                }
            }
            return 0;
        }
        return RI_ERR_SYNTAX;
    }
    if at(s, p_end) != 0 {
        return RI_ERR_SYNTAX;
    }

    // allocate InChI
    p_inp_inchi[b_mobile_h] = vec![InChI::default(); n_num_components as usize];

    // Pass 2. Count elements, save formulas and elements
    let mut p_start = 0usize;
    let mut i_component: i32 = 0;
    let mut n_num_atoms: i32 = -999;
    loop {
        let p_end = strchr(s, p_start, b'.').unwrap_or(p_start + strlen(s, p_start));
        let mut p = p_start;
        let mpy_component;
        if is_digit(at(s, p)) {
            let (v, q) = strtol(s, p);
            mpy_component = v as i32;
            if !(0..=MAX_ATOMS).contains(&mpy_component) {
                return RI_ERR_SYNTAX;
            }
            p = q;
        } else {
            mpy_component = 1;
        }
        if i_component + mpy_component > MAX_ATOMS {
            return RI_ERR_SYNTAX;
        }
        let formula: Vec<u8> = s[p.min(s.len())..p_end.min(s.len())].to_vec();
        for i in 0..mpy_component {
            if i_component + i >= n_num_components {
                return RI_ERR_SYNTAX;
            }
            let ic = i_component as usize;
            let ici = (i_component + i) as usize;
            p_inp_inchi[b_mobile_h][ici].sz_hill_formula =
                Some(String::from_utf8_lossy(&formula).into_owned());
            if i == 0 {
                // Pass 2.1 Parse formula and count atoms except H
                n_num_atoms = 0;
                let mut n_num_h: i32 = 0;
                let items = match formula_items(&formula) {
                    Ok(v) => v,
                    Err(r) => return r,
                };
                for (el, m) in &items {
                    if el == "H" {
                        n_num_h = n_num_h.wrapping_add(*m);
                        continue;
                    }
                    n_num_atoms = n_num_atoms.wrapping_add(*m);
                }
                if n_num_atoms > MAX_ATOMS {
                    return RI_ERR_SYNTAX;
                }
                *na_total = na_total.wrapping_add(mpy_component.wrapping_mul(n_num_atoms));
                let n_num_atoms_and_h = if n_num_atoms != 0 {
                    n_num_atoms
                } else {
                    i32::from(n_num_h > 0)
                };
                p_inp_inchi[b_mobile_h][ici].n_number_of_atoms = n_num_atoms_and_h;
                // Pass 2.2 Store elements; this assumes no bridging H
                p_inp_inchi[b_mobile_h][ici].n_atom =
                    vec![0; n_num_atoms_and_h.max(0) as usize + 1];
                let mut atoms: Vec<u8> = Vec::new();
                if n_num_atoms > 0 {
                    for (el, m) in &items {
                        if el == "H" {
                            continue;
                        }
                        let el_number = get_periodic_table_number(el);
                        if el_number == ERR_ELEM {
                            return RI_ERR_SYNTAX;
                        }
                        for _ in 0..*m {
                            if atoms.len() as i32 >= n_num_atoms {
                                return RI_ERR_PROGR;
                            }
                            p_inp_inchi[b_mobile_h][ici].n_atom[atoms.len()] = el_number as u8;
                            atoms.push(el_number as u8);
                        }
                    }
                } else if n_num_h > 0 {
                    atoms.push(1);
                    n_num_atoms = 1;
                }
                let n_atoms = atoms.len() as i32;
                atoms.push(0);
                p_inp_inchi[b_mobile_h][ici].n_atom = atoms;
                if n_num_atoms != n_atoms {
                    return RI_ERR_PROGR;
                }
            } else {
                // copy the duplicated formula and atoms
                let first = p_inp_inchi[b_mobile_h][ic].clone();
                let c = &mut p_inp_inchi[b_mobile_h][ici];
                c.sz_hill_formula = first.sz_hill_formula;
                c.n_number_of_atoms = n_num_atoms;
                c.n_atom = (0..=n_num_atoms as usize)
                    .map(|k| first.n_atom.get(k).copied().unwrap_or(0))
                    .collect();
            }
        }
        i_component += mpy_component;
        if at(s, p_end) != 0 {
            if at(s, p_end) != b'.' {
                return RI_ERR_SYNTAX;
            }
            p_start = p_end + 1;
        } else {
            break;
        }
    }
    if i_component != n_num_components {
        return RI_ERR_PROGR;
    }
    if b_mobile_h == TAUT_NON {
        // the exact number of atoms including bridging H is known from TAUT_YES
        let n = (n_num_components.min(pn_num_components[alt])).max(0) as usize;
        for i in 0..n {
            let len = p_inp_inchi[alt][i].n_number_of_atoms;
            if p_inp_inchi[b_mobile_h][i].n_number_of_atoms < len {
                let mut na: Vec<u8> = (0..len as usize)
                    .map(|k| p_inp_inchi[alt][i].n_atom.get(k).copied().unwrap_or(0))
                    .collect();
                na.push(0);
                p_inp_inchi[b_mobile_h][i].n_atom = na;
                p_inp_inchi[b_mobile_h][i].n_number_of_atoms = len;
            }
        }
    }
    n_num_components + 1
}

const ALPHA_BASE: i64 = 27;

/// `abctol`: `strtol` for the compressed InChI numbers (`Abc@` digits; the
/// sign is skipped, not applied)
fn abctol(s: &[u8], from: usize) -> (i64, usize) {
    let treshold = (i64::MAX - 1) / ALPHA_BASE;
    let mut val: i64 = 0;
    let mut p = from;
    if at(s, p) == b'-' {
        p += 1;
    }
    if at(s, p) == b'@' {
        return (0, p + 1);
    }
    if !upper(s, p) {
        return (0, from);
    }
    val += i64::from(at(s, p).to_ascii_lowercase() - b'a') + 1;
    p += 1;
    while at(s, p) != 0 {
        if lower(s, p) {
            val = val.wrapping_mul(ALPHA_BASE);
            val = val.wrapping_add(i64::from(at(s, p) - b'a') + 1);
        } else if at(s, p) == b'@' {
            if val > treshold {
                return (0, from);
            }
            val *= ALPHA_BASE;
        } else {
            break;
        }
        p += 1;
    }
    (val, p)
}

/// `inchi_strtol`: base 10 or the compressed base (27)
fn inchi_strtol(s: &[u8], from: usize, base: i64) -> (i64, usize) {
    if base == ALPHA_BASE {
        abctol(s, from)
    } else {
        strtol(s, from)
    }
}

/// `ONE_LINKED_BOND`
#[derive(Debug, Clone, Copy, Default)]
struct OneLinkedBond {
    neigh: u16,
    prev: u16,
}

/// `LINKED_BONDS`
#[derive(Debug, Default)]
struct LinkedBonds {
    p_bond: Option<Vec<OneLinkedBond>>,
    len: usize,
    len_alloc: usize,
}

const LINKED_BOND_ADD: usize = 128;

/// `AddLinkedBond`
fn add_linked_bond(at1: u16, at2: u16, num_at: u16, lb: &mut LinkedBonds) -> i32 {
    let req_len = (2 * usize::from(num_at) + 2).max(lb.len + 2);
    if i32::from(at2) > i32::from(num_at)
        || i32::from(at1) > i32::from(num_at)
        || i32::from(num_at) > MAX_ATOMS
    {
        return RI_ERR_SYNTAX;
    }
    if lb.len_alloc <= req_len {
        let new_len = req_len - req_len % LINKED_BOND_ADD + 2 * LINKED_BOND_ADD;
        let mut p = vec![OneLinkedBond::default(); new_len];
        if let (Some(old), true) = (&lb.p_bond, lb.len > 0) {
            p[..lb.len].copy_from_slice(&old[..lb.len]);
        }
        lb.p_bond = Some(p);
        lb.len_alloc = new_len;
    }
    let pb = lb.p_bond.as_mut().expect("allocated above");
    if lb.len == 0 {
        lb.len = usize::from(num_at) + 1;
        for b in &mut pb[..=usize::from(num_at)] {
            *b = OneLinkedBond::default();
        }
    }
    let (a1, a2) = (usize::from(at1), usize::from(at2));
    let prev = pb[a1].prev;
    if prev == 0 {
        pb[a1].neigh = at2;
        pb[a1].prev = at1;
    } else {
        pb[lb.len] = OneLinkedBond { neigh: at2, prev };
        pb[a1].prev = lb.len as u16;
        lb.len += 1;
    }
    let prev = pb[a2].prev;
    if prev == 0 {
        pb[a2].neigh = at1;
        pb[a2].prev = at2;
    } else {
        pb[lb.len] = OneLinkedBond { neigh: at1, prev };
        pb[a2].prev = lb.len as u16;
        lb.len += 1;
    }
    0
}

/// `ParseSegmentConnections`: the `/c` layer (`s` starts at the `c`, or at
/// whatever follows when the layer is missing). `p_inp_inchi` is empty for
/// the C NULL; `pb_abc` is -1 until the format (compressed or not) is known.
pub fn parse_segment_connections(
    s: &[u8],
    p_inp_inchi: &mut Vec<InChI>,
    pn_num_components: &mut i32,
    pb_abc: &mut i32,
    nb_total: &mut i32,
) -> i32 {
    const LAST_AT_LEN: usize = 256;
    *nb_total = 0;
    let mut i_component: i32 = 0;
    let mut lb = LinkedBonds::default();
    let mut b_prev_version: i32 = -1;

    if at(s, 0) != b'c' {
        if p_inp_inchi.is_empty() && *pn_num_components == 0 {
            let c = InChI {
                sz_hill_formula: Some(String::new()),
                n_conn_table: vec![0],
                len_conn_table: 0,
                ..Default::default()
            };
            *p_inp_inchi = vec![c];
            *pn_num_components = 1;
        } else {
            let n = *pn_num_components;
            for i in 0..n as usize {
                if let Some(c) = p_inp_inchi.get_mut(i) {
                    c.n_conn_table = vec![1, 0];
                    c.len_conn_table = 1;
                }
            }
        }
        return 0;
    }

    // Pass 1. Re-Count atoms, count bonds
    let mut p_start = 1usize;
    let n_num_components = *pn_num_components;
    if n_num_components > MAX_ATOMS {
        return RI_ERR_SYNTAX;
    }
    let ret: i32 = 'main: loop {
        let p_end = strchr(s, p_start, b';').unwrap_or(p_start + strlen(s, p_start));
        let mpy_component: i32;
        let mut p: usize;
        if let Some(star) = strchr(s, p_start, b'*').filter(|&x| x < p_end) {
            let (v, q) = strtol(s, p_start);
            mpy_component = v as i32;
            if !(0..=MAX_ATOMS).contains(&mpy_component) {
                break 'main RI_ERR_SYNTAX;
            }
            if star != q || !is_digit(at(s, p_start)) {
                break 'main RI_ERR_SYNTAX;
            }
            p = star + 1;
        } else {
            mpy_component = 1;
            p = p_start;
        }
        if i_component + mpy_component > MAX_ATOMS {
            break 'main RI_ERR_SYNTAX;
        }
        p_start = p;
        // Pass 1.1 parse a component
        let mut num_open: i32 = 0;
        let mut state: u8 = 0;
        let mut max_atom: u16 = 0;
        let mut n_num_bonds: i32 = 0;
        if p < p_end && *pb_abc == -1 {
            *pb_abc = i32::from(upper(s, p));
        }
        let base: i64 = if *pb_abc != 0 { ALPHA_BASE } else { 10 };
        if *pb_abc == 1 {
            let mut n_num_atoms: i32 = 1;
            while p < p_end {
                if at(s, p) == b'-' {
                    if b_prev_version == -1 {
                        b_prev_version = 1;
                    } else if b_prev_version != 1 {
                        break 'main RI_ERR_SYNTAX;
                    }
                    n_num_atoms = n_num_atoms.wrapping_sub(1);
                    p += 1;
                }
                if is_digit(at(s, p)) {
                    if b_prev_version == -1 {
                        b_prev_version = 0;
                    } else if b_prev_version != 0 {
                        break 'main RI_ERR_SYNTAX;
                    }
                    let (v, q) = strtol(s, p);
                    p = q;
                    n_num_atoms = (i64::from(n_num_atoms).wrapping_sub(v)) as i32;
                }
                let (v, q) = inchi_strtol(s, p, base);
                let cur_atom = v as u16;
                if at(s, p) != b'-' && cur_atom != 0 {
                    n_num_atoms = n_num_atoms.wrapping_add(1);
                    n_num_bonds += 1;
                    p = q;
                    if max_atom < cur_atom {
                        max_atom = cur_atom;
                    }
                } else {
                    break 'main RI_ERR_SYNTAX;
                }
            }
            if i32::from(max_atom) < n_num_atoms && n_num_bonds != 0 {
                max_atom = n_num_atoms as u16;
            }
        } else {
            while p < p_end {
                let c = at(s, p);
                p += 1;
                match c {
                    b'(' | b')' | b',' | b'-' => {
                        if state != b'N' {
                            break 'main RI_ERR_SYNTAX;
                        }
                        state = c;
                        num_open += i32::from(c == b'(') - i32::from(c == b')');
                        if num_open < 0 {
                            break 'main RI_ERR_SYNTAX;
                        }
                    }
                    _ => {
                        let (v, q) = strtol(s, p - 1);
                        let cur_atom = v as u16;
                        if is_digit(c) && cur_atom != 0 {
                            p = q;
                            match state {
                                b'(' | b')' | b',' | b'-' | 0 => {
                                    if state != 0 {
                                        n_num_bonds += 1;
                                    }
                                    if max_atom < cur_atom {
                                        max_atom = cur_atom;
                                    }
                                    state = b'N';
                                }
                                _ => break 'main RI_ERR_SYNTAX,
                            }
                        } else {
                            break 'main RI_ERR_SYNTAX;
                        }
                    }
                }
            }
            if num_open != 0 {
                break 'main RI_ERR_SYNTAX;
            }
        }
        *nb_total += n_num_bonds;

        // Save the results and allocate memory
        let n_num_atoms = i32::from(max_atom);
        let mut len_conn_table = n_num_atoms + n_num_bonds;
        for i in 0..mpy_component {
            let ici = i_component + i;
            if ici > n_num_components - 1 || ici < 0 {
                break 'main RI_ERR_SYNTAX;
            }
            let ici = ici as usize;
            let Some(c) = p_inp_inchi.get_mut(ici) else {
                break 'main RI_ERR_SYNTAX;
            };
            if c.n_number_of_atoms <= 0 || c.n_number_of_atoms > MAX_ATOMS {
                break 'main RI_ERR_SYNTAX;
            }
            let mut j = c.n_number_of_atoms;
            if j < n_num_atoms {
                let mut tmp: Vec<u8> = (0..j as usize)
                    .map(|k| c.n_atom.get(k).copied().unwrap_or(0))
                    .collect();
                while j < n_num_atoms {
                    tmp.push(1);
                    j += 1;
                }
                tmp.push(0);
                c.n_atom = tmp;
                c.n_number_of_atoms = n_num_atoms;
            } else if j > n_num_atoms && (len_conn_table != 0 || j != 1) {
                break 'main RI_ERR_SYNTAX;
            }
            if n_num_atoms == 0 && n_num_bonds == 0 && len_conn_table == 0 {
                len_conn_table = 1;
            }
            c.n_conn_table = vec![0; len_conn_table as usize + 1];
            c.len_conn_table = len_conn_table;
        }

        // Pass 1.2 parse a component and extract the bonds
        let mut num_open: usize = 0;
        let mut last_atom = [0u16; LAST_AT_LEN];
        let mut state: u8 = 0;
        p = p_start;
        lb.len = 0;
        let n_at = n_num_atoms as u16;
        if *pb_abc == 1 {
            last_atom[num_open] = 2;
            while p < p_end {
                if i32::from(last_atom[num_open]) > i32::from(max_atom) {
                    break 'main RI_ERR_SYNTAX;
                }
                if upper(s, p) {
                    let (v, q) = inchi_strtol(s, p, base);
                    let cur_atom = v as u16;
                    let r = add_linked_bond(last_atom[num_open], cur_atom, n_at, &mut lb);
                    if r != 0 {
                        break 'main r;
                    }
                    p = q;
                    if b_prev_version == 1 {
                        while p < p_end && at(s, p) == b'-' {
                            p += 1;
                            let (v, q) = inchi_strtol(s, p, base);
                            let cur_atom = v as u16;
                            if cur_atom != 0 {
                                let r =
                                    add_linked_bond(last_atom[num_open], cur_atom, n_at, &mut lb);
                                if r != 0 {
                                    break 'main r;
                                }
                                p = q;
                            } else {
                                break 'main RI_ERR_SYNTAX;
                            }
                        }
                    } else if b_prev_version == 0 && is_digit(at(s, p)) {
                        let (v, q) = strtol(s, p);
                        let mut num_neigh = v as i32;
                        p = q;
                        loop {
                            let cont = num_neigh != 0;
                            num_neigh = num_neigh.wrapping_sub(1);
                            if !(cont && p < p_end) {
                                break;
                            }
                            let (v, q) = inchi_strtol(s, p, base);
                            let cur_atom = v as u16;
                            if cur_atom != 0 {
                                let r =
                                    add_linked_bond(last_atom[num_open], cur_atom, n_at, &mut lb);
                                if r != 0 {
                                    break 'main r;
                                }
                                p = q;
                            } else {
                                break 'main RI_ERR_SYNTAX;
                            }
                        }
                    }
                    last_atom[num_open] = last_atom[num_open].wrapping_add(1);
                } else {
                    break 'main RI_ERR_SYNTAX;
                }
            }
        } else {
            while p < p_end {
                let c = at(s, p);
                p += 1;
                match c {
                    b'(' | b')' | b',' | b'-' => {
                        if state == b'N' {
                            state = c;
                        } else {
                            break 'main RI_ERR_SYNTAX;
                        }
                    }
                    _ => {
                        let (v, q) = strtol(s, p - 1);
                        let cur_atom = v as u16;
                        if !(is_digit(c) && cur_atom != 0) {
                            break 'main RI_ERR_SYNTAX;
                        }
                        p = q;
                        match state {
                            0 => {
                                last_atom[num_open] = cur_atom;
                                state = b'N';
                            }
                            b'(' => {
                                let r =
                                    add_linked_bond(last_atom[num_open], cur_atom, n_at, &mut lb);
                                if r != 0 {
                                    break 'main r;
                                }
                                num_open += 1;
                                if num_open >= LAST_AT_LEN {
                                    break 'main RI_ERR_PROGR;
                                }
                                last_atom[num_open] = cur_atom;
                                state = b'N';
                            }
                            b')' => {
                                if num_open == 0 {
                                    break 'main RI_ERR_SYNTAX;
                                }
                                num_open -= 1;
                                let r =
                                    add_linked_bond(last_atom[num_open], cur_atom, n_at, &mut lb);
                                if r != 0 {
                                    break 'main r;
                                }
                                last_atom[num_open] = cur_atom;
                                state = b'N';
                            }
                            b',' => {
                                if num_open == 0 {
                                    break 'main RI_ERR_SYNTAX;
                                }
                                let r = add_linked_bond(
                                    last_atom[num_open - 1],
                                    cur_atom,
                                    n_at,
                                    &mut lb,
                                );
                                if r != 0 {
                                    break 'main r;
                                }
                                last_atom[num_open] = cur_atom;
                                state = b'N';
                            }
                            b'-' => {
                                let r =
                                    add_linked_bond(last_atom[num_open], cur_atom, n_at, &mut lb);
                                if r != 0 {
                                    break 'main r;
                                }
                                last_atom[num_open] = cur_atom;
                                state = b'N';
                            }
                            _ => break 'main RI_ERR_SYNTAX,
                        }
                    }
                }
            }
        }

        // Store the bonds in connection table
        let ic = i_component as usize;
        if len_conn_table > 1 {
            let mut m: usize = 0;
            for i in 0..n_num_atoms as usize {
                let Some(pb) = lb.p_bond.as_ref() else {
                    break 'main RI_ERR_SYNTAX;
                };
                let mut neighbor = [0u16; MAXVAL];
                let mut k = 0usize;
                let mut j = usize::from(pb[i + 1].prev);
                if j != 0 {
                    while k < MAXVAL {
                        neighbor[k] = pb[j].neigh;
                        k += 1;
                        if j == i + 1 {
                            break;
                        }
                        j = usize::from(pb[j].prev);
                    }
                }
                if j != i + 1 {
                    break 'main RI_ERR_SYNTAX;
                }
                insertions_sort_at_numb(&mut neighbor, k);
                let ct_len = p_inp_inchi[ic].len_conn_table as usize;
                if m == ct_len {
                    break 'main RI_ERR_SYNTAX;
                }
                for j in 1..k {
                    if neighbor[j] == neighbor[j - 1] {
                        break 'main RI_ERR_SYNTAX;
                    }
                }
                p_inp_inchi[ic].n_conn_table[m] = (i + 1) as u16;
                m += 1;
                let mut j = 0;
                while j < k && usize::from(neighbor[j]) <= i {
                    if m == ct_len {
                        break 'main RI_ERR_SYNTAX;
                    }
                    p_inp_inchi[ic].n_conn_table[m] = neighbor[j];
                    m += 1;
                    j += 1;
                }
            }
            if m as i32 != len_conn_table {
                break 'main RI_ERR_PROGR;
            }
        } else {
            p_inp_inchi[ic].n_conn_table[0] = 1;
        }

        // Duplicate if needed
        let mut i = 1;
        while i < mpy_component {
            let (first, rest) = p_inp_inchi.split_at_mut(ic + 1);
            let dst = &mut rest[(i - 1) as usize];
            if dst.n_conn_table.is_empty() || dst.len_conn_table != len_conn_table {
                break 'main RI_ERR_PROGR;
            }
            let n = len_conn_table as usize;
            dst.n_conn_table[..n].copy_from_slice(&first[ic].n_conn_table[..n]);
            i += 1;
        }
        i_component += i.max(1);
        if at(s, p_end) != 0 {
            p_start = p_end + 1;
        } else {
            break 'main i_component;
        }
        if i_component == n_num_components {
            break 'main i_component;
        }
    };
    ret
}

/// `GetInChIFormulaNumH`: the number of H in the Hill formula including the
/// bridging H; `Err` is the return code
pub fn get_inchi_formula_num_h(inchi: &InChI) -> Result<i32, i32> {
    let mut n = 0i32;
    if let Some(f) = &inchi.sz_hill_formula {
        let f = f.as_bytes();
        let mut from = 0;
        while let Some(h) = strchr(f, from, b'H') {
            let p = h + 1;
            let mut next = p;
            if !lower(f, p) {
                if is_digit(at(f, p)) {
                    let (v, q) = strtol(f, p);
                    let delta = v as i32;
                    if !(0..=MAX_ATOMS).contains(&delta) {
                        return Err(RI_ERR_SYNTAX);
                    }
                    n = n.wrapping_add(delta);
                    next = q;
                } else {
                    n += 1;
                }
            }
            from = next;
        }
    }
    Ok(n)
}

/// `GetInChINumH`: the number of H of the connection table (immobile and
/// bridging H and the H of the tautomeric groups); `Err` is the return code
pub fn get_inchi_num_h(inchi: &InChI) -> Result<i32, i32> {
    let mut n = 0i32;
    for i in 0..inchi.n_number_of_atoms.max(0) as usize {
        n += i32::from(inchi.n_atom.get(i).copied().unwrap_or(0) == 1);
        n += i32::from(inchi.n_num_h.get(i).copied().unwrap_or(0));
    }
    if inchi.len_tautomer > 3 && !inchi.n_tautomer.is_empty() {
        let len = inchi.len_tautomer;
        let g = |j: usize| i32::from(inchi.n_tautomer.get(j).copied().unwrap_or(0));
        let mut j: i32 = 0;
        let n_groups = g(j as usize);
        j += 1;
        let mut i_group = 0;
        while j < len && i_group < n_groups {
            let group_len = g(j as usize) + 1;
            n += g(j as usize + 1);
            i_group += 1;
            j += group_len;
        }
        if i_group != n_groups || j != len {
            return Err(RI_ERR_PROGR);
        }
    }
    if inchi.n_num_h_fixed.is_some() && (inchi.len_tautomer != 0 || !inchi.n_tautomer.is_empty()) {
        return Err(RI_ERR_PROGR);
    }
    Ok(n)
}

const INCHI_T_NUM_MOVABLE: i32 = 2;

/// `strspn(s + from, set)` end position
fn strspn_end(s: &[u8], from: usize, set: &[u8]) -> usize {
    let mut i = from;
    while at(s, i) != 0 && set.contains(&s[i]) {
        i += 1;
    }
    i
}

/// the H counts of the component `ic`: `nNum_H` or `nNum_H_fixed`
fn num_h_mut(c: &mut InChI, mobile: usize) -> Option<&mut Vec<i8>> {
    if mobile == TAUT_YES {
        if c.n_num_h.is_empty() {
            None
        } else {
            Some(&mut c.n_num_h)
        }
    } else {
        c.n_num_h_fixed.as_mut()
    }
}

fn set_num_h(c: &mut InChI, mobile: usize, from: i32, to: i32, num_h: i32, sum: &mut i32) -> bool {
    if let Some(v) = num_h_mut(c, mobile) {
        for i in from..=to {
            let Ok(k) = usize::try_from(i - 1) else {
                return false;
            };
            // the proton placeholder has two entries only; C writes past them
            if k >= v.len() {
                v.resize(k + 1, 0);
            }
            v[k] = num_h as i8;
            *sum = sum.wrapping_add(num_h);
        }
    }
    true
}

/// `ParseSegmentMobileH`: the `/h` layer (`s` starts at the `h`, or at
/// whatever follows when the layer is missing)
pub fn parse_segment_mobile_h(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    pn_num_components: [i32; 2],
    pb_abc: &mut i32,
) -> i32 {
    let alt = 1 - b_mobile_h;
    let n_num_components = pn_num_components[b_mobile_h];

    // number of immobile H is always allocated; immobile H are present in M layer only
    for i in 0..n_num_components.max(0) as usize {
        let Some(c) = p_inp_inchi[b_mobile_h].get(i) else {
            return RI_ERR_SYNTAX;
        };
        let mut len = c.n_number_of_atoms;
        let has_alt = b_mobile_h == TAUT_NON && (i as i32) < pn_num_components[alt];
        if has_alt {
            let alt_n = p_inp_inchi[alt][i].n_number_of_atoms;
            if len < alt_n {
                len = alt_n;
                p_inp_inchi[b_mobile_h][i].n_num_h = Vec::new();
            }
        }
        let len = (len + 1).max(0) as usize;
        if p_inp_inchi[b_mobile_h][i].n_num_h.is_empty() {
            p_inp_inchi[b_mobile_h][i].n_num_h = vec![0; len];
        }
        if has_alt {
            // copy immobile H from the Mobile-H layer to the Fixed-H layer
            let a = &mut p_inp_inchi[alt][i].n_num_h;
            a.resize(len, 0);
            let copy: Vec<i8> = (0..len)
                .map(|k| if k + 1 < len { a[k] } else { 0 })
                .collect();
            p_inp_inchi[b_mobile_h][i].n_num_h = copy;
        }
    }

    if at(s, 0) != b'h' {
        return 0;
    }

    // Read Hydrogen info in 1 pass
    let mut p_start = 1usize;
    let mut i_component: i32 = 0;
    loop {
        let p_end = strchr(s, p_start, b';').unwrap_or(p_start + strlen(s, p_start));
        let mpy_component: i32;
        let mut p: usize;
        if let Some(star) = strchr(s, p_start, b'*').filter(|&x| x < p_end) {
            let (v, q) = strtol(s, p_start);
            mpy_component = v as i32;
            if !(0..=MAX_ATOMS).contains(&mpy_component) {
                return RI_ERR_SYNTAX;
            }
            if star != q || !is_digit(at(s, p_start)) {
                return RI_ERR_SYNTAX;
            }
            p = star + 1;
        } else {
            mpy_component = 1;
            p = p_start;
        }
        p_start = p;
        // Pass 1.1 parse a component
        let ic = i_component as usize;
        if ic >= p_inp_inchi[b_mobile_h].len() {
            return RI_ERR_SYNTAX;
        }
        let mut num_ct_atoms = p_inp_inchi[b_mobile_h][ic].n_number_of_atoms;
        if b_mobile_h == TAUT_NON && i_component < pn_num_components[alt] {
            num_ct_atoms = p_inp_inchi[alt][ic].n_number_of_atoms;
        }
        if p < p_end && *pb_abc == -1 {
            *pb_abc = i32::from(at(s, p) == b',' || upper(s, p));
        }
        let base: i64 = if *pb_abc == 1 { ALPHA_BASE } else { 10 };

        // immobile H: locate the first tautomer group character
        let t0 = if *pb_abc == 1 {
            strchr(s, p, b',')
        } else {
            strchr(s, p, b'(')
        };
        if t0.is_some() && b_mobile_h == TAUT_NON {
            return RI_ERR_SYNTAX;
        }
        let (t, p_taut) = match t0 {
            Some(x) if x <= p_end => (Some(x), x),
            _ => (None, p_end),
        };
        let mut i = 0;
        while i < mpy_component && i_component + i < n_num_components {
            let k = (i_component + i) as usize;
            if b_mobile_h == TAUT_NON {
                if p_inp_inchi[b_mobile_h][k].n_num_h_fixed.is_some() {
                    return RI_ERR_PROGR;
                }
                let own_n = p_inp_inchi[b_mobile_h][k].n_number_of_atoms;
                let in_alt = i_component + i < pn_num_components[alt];
                let alt_n = if in_alt {
                    p_inp_inchi[alt][k].n_number_of_atoms
                } else {
                    0
                };
                let len = if in_alt {
                    own_n.max(alt_n) + 1
                } else {
                    own_n + 1
                };
                p_inp_inchi[b_mobile_h][k].n_num_h_fixed = Some(vec![0; len.max(0) as usize]);
                if in_alt {
                    let len2 = own_n.min(alt_n).max(0) as usize;
                    let own_atoms = p_inp_inchi[b_mobile_h][k].n_atom.clone();
                    if !own_atoms.is_empty() && len2 != 0 {
                        let a: Vec<u8> = (0..len2)
                            .map(|j| own_atoms.get(j).copied().unwrap_or(0))
                            .collect();
                        let b: Vec<u8> = (0..len2)
                            .map(|j| p_inp_inchi[alt][k].n_atom.get(j).copied().unwrap_or(0))
                            .collect();
                        if a != b {
                            return RI_ERR_SYNTAX;
                        }
                    }
                    if own_n < alt_n {
                        let alt_atoms = p_inp_inchi[alt][k].n_atom.clone();
                        let mut na = vec![0u8; len.max(0) as usize];
                        if len > 1 {
                            for j in 0..(len - 1) as usize {
                                na[j] = alt_atoms.get(j).copied().unwrap_or(0);
                            }
                        }
                        p_inp_inchi[b_mobile_h][k].n_atom = na;
                        p_inp_inchi[b_mobile_h][k].n_number_of_atoms = alt_n;
                    }
                }
            }
            i += 1;
        }

        let mut tg_alloc_len: i32 = 0;
        let mut num_h_component: i32 = 0;
        let mut num_taut_h_component: i32 = 0;
        let mob = b_mobile_h;
        if *pb_abc == 1 {
            // read numbers of H: XnYn... or XYn...
            p = p_start;
            while p < p_taut {
                if at(s, p) == 0 || !upper(s, p) {
                    return RI_ERR_SYNTAX;
                }
                let (v, q) = inchi_strtol(s, p, base);
                let cur_atom = v as u16;
                let mut nxt_atom = cur_atom;
                if cur_atom != 0 {
                    p = q;
                    if upper(s, p) {
                        let (v, q) = inchi_strtol(s, p, base);
                        nxt_atom = v as u16;
                        if i32::from(nxt_atom) > MAX_ATOMS {
                            return RI_ERR_SYNTAX;
                        }
                        p = q;
                    }
                }
                if cur_atom > nxt_atom || i32::from(nxt_atom) > num_ct_atoms || p > p_taut {
                    return RI_ERR_SYNTAX;
                }
                let (v, q) = strtol(s, p);
                let num_h = v as i32;
                if num_h == 0 || q > p_taut {
                    return RI_ERR_SYNTAX;
                }
                p = q;
                if !set_num_h(
                    &mut p_inp_inchi[mob][ic],
                    mob,
                    i32::from(cur_atom),
                    i32::from(nxt_atom),
                    num_h,
                    &mut num_h_component,
                ) {
                    return RI_ERR_SYNTAX;
                }
            }
            if p != p_taut {
                return RI_ERR_SYNTAX;
            }
        } else {
            // read numbers of H: 1-2,3H2,4,5H3
            p = p_start;
            let set: &[u8] = if b_mobile_h == TAUT_NON {
                b"0123456789hDHT-,;()*"
            } else {
                b"0123456789DHT-,;()*"
            };
            if at(s, strspn_end(s, 1, set)) != 0 {
                return RI_ERR_SYNTAX;
            }
            while p < p_taut {
                if at(s, p) == 0 || !is_digit(at(s, p)) {
                    return RI_ERR_SYNTAX;
                }
                let mut h = p;
                while at(s, h) != 0 && at(s, h) != b'H' && at(s, h) != b'h' {
                    h += 1;
                }
                if at(s, h) == 0 || h >= p_taut {
                    return RI_ERR_SYNTAX;
                }
                let mut num_h: i32 = match at(s, h) {
                    b'H' => 1,
                    b'h' => -1,
                    _ => 0,
                };
                if num_h == 0 {
                    return RI_ERR_SYNTAX;
                }
                let mut p1;
                if is_digit(at(s, h + 1)) {
                    let (v, e) = strtol(s, h + 1);
                    num_h = num_h.wrapping_mul(v as i32);
                    p1 = e;
                } else {
                    p1 = h + 1;
                }
                if at(s, p1) == b',' {
                    p1 += 1;
                }
                while p < h {
                    if at(s, p) == 0 || !is_digit(at(s, p)) {
                        return RI_ERR_SYNTAX;
                    }
                    let (v, mut q) = strtol(s, p);
                    let cur_atom = v as i32;
                    let mut nxt_atom = cur_atom;
                    if at(s, q) == b'-' {
                        let (v2, q2) = strtol(s, q + 1);
                        nxt_atom = v2 as i32;
                        q = q2;
                        if !(0..=MAX_ATOMS).contains(&nxt_atom) {
                            return RI_ERR_SYNTAX;
                        }
                    }
                    if cur_atom == 0
                        || cur_atom > num_ct_atoms
                        || nxt_atom < cur_atom
                        || nxt_atom > num_ct_atoms
                    {
                        return RI_ERR_SYNTAX;
                    }
                    if !set_num_h(
                        &mut p_inp_inchi[mob][ic],
                        mob,
                        cur_atom,
                        nxt_atom,
                        num_h,
                        &mut num_h_component,
                    ) {
                        return RI_ERR_SYNTAX;
                    }
                    p = q;
                    if at(s, p) == b',' {
                        p += 1;
                    }
                }
                if p == h {
                    p = p1;
                } else if p == p_taut {
                    break;
                } else {
                    return RI_ERR_SYNTAX;
                }
            }
        }

        // tautomeric groups
        let mut len_tautomer: i32 = 0;
        if let Some(t) = t {
            p = t;
            let mut tg_pos_tautomer: i32;
            let mut i_tgroup: i32 = 0;
            let mut num_atoms: i32 = 0;
            let alloc_taut = |p_inp: &mut [Vec<InChI>; 2], tg_alloc_len: i32| -> bool {
                for i in 0..mpy_component {
                    let k = (i_component + i) as usize;
                    let Some(c) = p_inp[mob].get_mut(k) else {
                        return false;
                    };
                    c.n_tautomer = vec![0; tg_alloc_len as usize + 1];
                    c.len_tautomer = 0;
                }
                true
            };
            if *pb_abc == 1 {
                if at(s, p) == b',' {
                    p += 1;
                } else {
                    return RI_ERR_SYNTAX;
                }
                while p < p_end {
                    if !is_digit(at(s, p)) {
                        return RI_ERR_SYNTAX;
                    }
                    let (v, q) = strtol(s, p);
                    let num_h = v as i32;
                    if num_h == 0 || q > p_end {
                        return RI_ERR_SYNTAX;
                    }
                    p = q;
                    let mut num_minus = 0;
                    if at(s, p) == b'-' {
                        p += 1;
                        if is_digit(at(s, p)) {
                            let (v, q) = strtol(s, p);
                            num_minus = v as i32;
                            p = q;
                        } else {
                            num_minus = 1;
                        }
                    }
                    if p >= p_end {
                        return RI_ERR_SYNTAX;
                    }
                    if tg_alloc_len == 0 {
                        tg_alloc_len = ((3 + INCHI_T_NUM_MOVABLE)
                            * p_inp_inchi[mob][ic].n_number_of_atoms)
                            / 2
                            + 1;
                        if !alloc_taut(p_inp_inchi, tg_alloc_len) {
                            return RI_ERR_SYNTAX;
                        }
                        tg_pos_tautomer = 1;
                    } else {
                        tg_pos_tautomer = len_tautomer;
                    }
                    if tg_pos_tautomer + 3 >= tg_alloc_len {
                        return RI_ERR_PROGR;
                    }
                    let c = &mut p_inp_inchi[mob][ic];
                    c.n_tautomer[tg_pos_tautomer as usize + 1] = num_h as u16;
                    c.n_tautomer[tg_pos_tautomer as usize + 2] = num_minus as u16;
                    len_tautomer = tg_pos_tautomer + 3;
                    num_taut_h_component = num_taut_h_component.wrapping_add(num_h);
                    while p < p_end && upper(s, p) {
                        let (v, q) = inchi_strtol(s, p, base);
                        let val = v as i32;
                        if len_tautomer >= tg_alloc_len || val > num_ct_atoms {
                            return RI_ERR_PROGR;
                        }
                        num_atoms += 1;
                        let c = &mut p_inp_inchi[mob][ic];
                        c.n_tautomer[len_tautomer as usize] = val as u16;
                        len_tautomer += 1;
                        p = q;
                    }
                    if num_atoms == 0 || (p < p_end && !is_digit(at(s, p))) {
                        return RI_ERR_PROGR;
                    }
                    i_tgroup += 1;
                    let c = &mut p_inp_inchi[mob][ic];
                    c.n_tautomer[tg_pos_tautomer as usize] =
                        (len_tautomer - tg_pos_tautomer - 1) as u16;
                    c.len_tautomer = len_tautomer;
                }
                if i_tgroup == 0 || p != p_end {
                    return RI_ERR_SYNTAX;
                }
                p_inp_inchi[mob][ic].n_tautomer[0] = i_tgroup as u16;
            } else {
                let mut state: u8 = b')';
                let (mut num_h, mut num_minus) = (-999i32, -999i32);
                tg_pos_tautomer = -999;
                while p < p_end {
                    match at(s, p) {
                        b'(' => {
                            if state != b')' {
                                return RI_ERR_SYNTAX;
                            }
                            state = b'(';
                            p += 1;
                            num_h = 0;
                            num_minus = 0;
                        }
                        b')' => {
                            if state != b'A' {
                                return RI_ERR_SYNTAX;
                            }
                            if tg_alloc_len == 0 {
                                return RI_ERR_SYNTAX;
                            }
                            i_tgroup += 1;
                            state = b')';
                            p += 1;
                            let c = &mut p_inp_inchi[mob][ic];
                            c.n_tautomer[tg_pos_tautomer as usize] =
                                (len_tautomer - tg_pos_tautomer - 1) as u16;
                            c.len_tautomer = len_tautomer;
                        }
                        b'H' => {
                            if state != b'(' {
                                return RI_ERR_SYNTAX;
                            }
                            state = b'H';
                            p += 1;
                            num_h = 1;
                        }
                        b'-' => {
                            if state != b'N' && state != b'H' {
                                return RI_ERR_SYNTAX;
                            }
                            state = b'-';
                            p += 1;
                            num_minus = 1;
                        }
                        b',' => match state {
                            b'N' | b'H' | b'-' | b'M' => {
                                if num_h <= 0 && num_minus <= 0 {
                                    return RI_ERR_SYNTAX;
                                }
                                if tg_alloc_len == 0 {
                                    tg_alloc_len = ((3 + INCHI_T_NUM_MOVABLE)
                                        * p_inp_inchi[mob][ic].n_number_of_atoms)
                                        / 2
                                        + 1;
                                    if !alloc_taut(p_inp_inchi, tg_alloc_len) {
                                        return RI_ERR_SYNTAX;
                                    }
                                    tg_pos_tautomer = 1;
                                } else {
                                    tg_pos_tautomer = len_tautomer;
                                }
                                if tg_pos_tautomer + 3 >= tg_alloc_len {
                                    return RI_ERR_PROGR;
                                }
                                let c = &mut p_inp_inchi[mob][ic];
                                c.n_tautomer[tg_pos_tautomer as usize + 1] = num_h as u16;
                                c.n_tautomer[tg_pos_tautomer as usize + 2] = num_minus as u16;
                                len_tautomer = tg_pos_tautomer + 3;
                                num_taut_h_component = num_taut_h_component.wrapping_add(num_h);
                                state = b',';
                                p += 1;
                            }
                            b'A' => {
                                state = b',';
                                p += 1;
                            }
                            _ => return RI_ERR_SYNTAX,
                        },
                        c if is_digit(c) => {
                            let (v, q) = strtol(s, p);
                            let val = v as i32;
                            if val <= 0 {
                                return RI_ERR_SYNTAX;
                            }
                            p = q;
                            match state {
                                b'H' => {
                                    num_h = val;
                                    state = b'N';
                                }
                                b'-' => {
                                    num_minus = val;
                                    state = b'M';
                                }
                                b',' => {
                                    if len_tautomer >= tg_alloc_len || val > num_ct_atoms {
                                        return RI_ERR_PROGR;
                                    }
                                    num_atoms += 1;
                                    let c = &mut p_inp_inchi[mob][ic];
                                    c.n_tautomer[len_tautomer as usize] = val as u16;
                                    len_tautomer += 1;
                                    state = b'A';
                                }
                                _ => return RI_ERR_SYNTAX,
                            }
                        }
                        _ => return RI_ERR_SYNTAX,
                    }
                }
                let _ = num_atoms;
                if i_tgroup == 0 || state != b')' {
                    return RI_ERR_SYNTAX;
                }
                p_inp_inchi[mob][ic].n_tautomer[0] = i_tgroup as u16;
            }
        }

        // check num_H in components
        let check = |c: &InChI| -> Result<(), i32> {
            let f = get_inchi_formula_num_h(c)?;
            let h = get_inchi_num_h(c)?;
            if f != h + if b_mobile_h == TAUT_NON {
                num_h_component
            } else {
                0
            } {
                return Err(RI_ERR_SYNTAX);
            }
            Ok(())
        };
        if let Err(r) = check(&p_inp_inchi[mob][ic]) {
            return r;
        }

        // duplicate according to multiplication
        let mut i = 1;
        while i < mpy_component {
            let k = i_component + i;
            if k > n_num_components - 1 || k < 0 {
                return RI_ERR_SYNTAX;
            }
            let k = k as usize;
            let dn = p_inp_inchi[mob][k].n_number_of_atoms;
            if dn <= 0 || dn > MAX_ATOMS {
                return RI_ERR_SYNTAX;
            }
            let src_h: Option<Vec<i8>> =
                num_h_mut(&mut p_inp_inchi[mob][ic], mob).map(|v| v.clone());
            if let Some(src) = src_h {
                if let Some(dst) = num_h_mut(&mut p_inp_inchi[mob][k], mob) {
                    for j in 0..dn as usize {
                        if let (Some(d), Some(sv)) = (dst.get_mut(j), src.get(j)) {
                            *d = *sv;
                        }
                    }
                }
            }
            let (src_t, src_len) = {
                let c = &p_inp_inchi[mob][ic];
                (c.n_tautomer.clone(), c.len_tautomer)
            };
            if !p_inp_inchi[mob][k].n_tautomer.is_empty() && !src_t.is_empty() && src_len != 0 {
                let d = &mut p_inp_inchi[mob][k];
                for j in 0..src_len as usize {
                    if let (Some(x), Some(y)) = (d.n_tautomer.get_mut(j), src_t.get(j)) {
                        *x = *y;
                    }
                }
                d.len_tautomer = src_len;
            }
            if let Err(r) = check(&p_inp_inchi[mob][k]) {
                return r;
            }
            i += 1;
        }

        // prepare for the next component
        i_component += i;
        if at(s, p_end) != 0 {
            if i_component >= n_num_components {
                return RI_ERR_SYNTAX;
            }
            p_start = p_end + 1;
        } else {
            break;
        }
    }
    if n_num_components != i_component {
        return RI_ERR_SYNTAX;
    }
    i_component + 1
}

// ---------------------------------------------------------------------------
// The stereo layers: `/b`, `/t`, `/m`, `/s` (and the isotopic ones)

/// `COPY_SEG_TYPE`
pub const CPY_SP2: i32 = 0;
pub const CPY_SP3: i32 = 1;
pub const CPY_SP3_M: i32 = 2;
pub const CPY_SP3_S: i32 = 3;
pub const CPY_ISO_AT: i32 = 4;

const MULT_TYPE: &[u8] = b"mnMNe";
const PARITY_TYPE: &[u8] = b"-+u?";
const AB_MIN_KNOWN_PARITY: i32 = 1;
const AB_MAX_KNOWN_PARITY: i32 = 4;

/// `strchr(set, c)` over a C string constant: the position (the terminating
/// zero is found as well when `c == 0`)
fn strchr_set(set: &[u8], c: u8) -> Option<usize> {
    if c == 0 {
        return Some(set.len());
    }
    set.iter().position(|&x| x == c)
}

fn stereo_slot(c: &mut InChI, iso: bool) -> &mut Option<InChIStereo> {
    if iso {
        &mut c.stereo_isotopic
    } else {
        &mut c.stereo
    }
}

/// `CopySegment`: copies the sp2/sp3/`/m`/`/s`/isotopic atoms segment of
/// `from` (`None` for the C NULL) into `to`. Returns the number of the
/// copied items plus one, 0 when there is nothing to copy, or a `RI_ERR_*`.
pub fn copy_segment(
    to: &mut InChI,
    from: Option<&InChI>,
    segment_type: i32,
    b_isotopic_to: i32,
    b_isotopic_from: i32,
) -> i32 {
    let stereo_from: Option<&InChIStereo> = match (from, b_isotopic_from) {
        (Some(f), 1) => f.stereo_isotopic.as_ref(),
        (Some(f), 0) => f.stereo.as_ref(),
        _ => None,
    };
    let iso_to = b_isotopic_to != 0;
    if matches!(segment_type, CPY_SP2 | CPY_SP3 | CPY_SP3_M | CPY_SP3_S) {
        if stereo_from.is_none() && b_isotopic_from >= 0 {
            return 0;
        }
        match segment_type {
            CPY_SP2 => {
                if b_isotopic_from >= 0 {
                    let f = from.unwrap();
                    if f.n_number_of_atoms > MAX_ATOMS || f.n_number_of_atoms < 0 {
                        return RI_ERR_SYNTAX;
                    }
                }
                let sf = stereo_from;
                if b_isotopic_from < 0
                    || sf.is_some_and(|s| {
                        !s.b_parity.is_empty()
                            && !s.n_bond_atom1.is_empty()
                            && !s.n_bond_atom2.is_empty()
                    })
                {
                    let len = if b_isotopic_from < 0 {
                        0
                    } else {
                        sf.unwrap().n_number_of_stereo_bonds
                    };
                    let src = sf.cloned();
                    let slot = stereo_slot(to, iso_to);
                    let st = slot.get_or_insert_with(Default::default);
                    if st.n_number_of_stereo_bonds > 0
                        || !st.b_parity.is_empty()
                        || !st.n_bond_atom1.is_empty()
                        || !st.n_bond_atom2.is_empty()
                    {
                        return RI_ERR_SYNTAX;
                    }
                    let n = (len.max(0) + 1) as usize;
                    st.b_parity = vec![0; n];
                    st.n_bond_atom1 = vec![0; n];
                    st.n_bond_atom2 = vec![0; n];
                    if b_isotopic_from >= 0 && len != 0 {
                        let s = src.unwrap();
                        for k in 0..n {
                            st.b_parity[k] = s.b_parity.get(k).copied().unwrap_or(0);
                            st.n_bond_atom1[k] = s.n_bond_atom1.get(k).copied().unwrap_or(0);
                            st.n_bond_atom2[k] = s.n_bond_atom2.get(k).copied().unwrap_or(0);
                        }
                    }
                    st.n_number_of_stereo_bonds = len;
                    return len + 1;
                }
                0
            }
            CPY_SP3 => {
                let sf = stereo_from;
                if b_isotopic_from < 0
                    || sf.is_some_and(|s| !s.t_parity.is_empty() && !s.n_number.is_empty())
                {
                    let len = if b_isotopic_from < 0 {
                        0
                    } else {
                        sf.unwrap().n_number_of_stereo_centers
                    };
                    let src = sf.cloned();
                    let slot = stereo_slot(to, iso_to);
                    let st = slot.get_or_insert_with(Default::default);
                    if st.n_number_of_stereo_centers > 0
                        || !st.t_parity.is_empty()
                        || !st.n_number.is_empty()
                    {
                        return RI_ERR_SYNTAX;
                    }
                    let n = (len.max(0) + 1) as usize;
                    st.t_parity = vec![0; n];
                    st.n_number = vec![0; n];
                    if b_isotopic_from >= 0 && len != 0 {
                        let s = src.unwrap();
                        for k in 0..n {
                            st.t_parity[k] = s.t_parity.get(k).copied().unwrap_or(0);
                            st.n_number[k] = s.n_number.get(k).copied().unwrap_or(0);
                        }
                    }
                    st.n_number_of_stereo_centers = len;
                    return len + 1;
                }
                0
            }
            CPY_SP3_M => {
                let inv = stereo_from.map(|s| s.n_comp_inv2_abs);
                let st = stereo_slot(to, iso_to).get_or_insert_with(Default::default);
                if st.n_comp_inv2_abs != 0 && NO_VALUE_INT != st.n_comp_inv2_abs {
                    return RI_ERR_SYNTAX;
                }
                st.n_comp_inv2_abs = if b_isotopic_from < 0 { 0 } else { inv.unwrap() };
                1
            }
            _ => {
                // CPY_SP3_S: uses `bIsotopicFrom` to select the target (as C does)
                let triv = stereo_from.map(|s| s.b_trivial_inv);
                let st = stereo_slot(to, b_isotopic_from != 0).get_or_insert_with(Default::default);
                if st.b_trivial_inv != 0 {
                    return RI_ERR_SYNTAX;
                }
                st.b_trivial_inv = if b_isotopic_from < 0 {
                    0
                } else {
                    triv.unwrap()
                };
                1
            }
        }
    } else if segment_type == CPY_ISO_AT {
        let Some(f) = from else {
            return RI_ERR_SYNTAX;
        };
        if b_isotopic_from < 0 || !f.isotopic_atom.is_empty() {
            let len = if b_isotopic_from < 0 {
                0
            } else {
                f.n_number_of_isotopic_atoms
            };
            if to.n_number_of_isotopic_atoms < 1 {
                to.isotopic_atom = Vec::new();
            }
            let n = (len.max(0) + 1) as usize;
            if to.isotopic_atom.is_empty() {
                to.isotopic_atom = vec![InChIIsotopicAtom::default(); n];
            }
            if to.n_number_of_isotopic_atoms != 0 {
                return RI_ERR_SYNTAX;
            }
            if b_isotopic_from >= 0 && len != 0 {
                for k in 0..n.min(to.isotopic_atom.len()) {
                    to.isotopic_atom[k] = f.isotopic_atom.get(k).cloned().unwrap_or_default();
                }
            }
            to.n_number_of_isotopic_atoms = len;
            return len + 1;
        }
        0
    } else {
        RI_ERR_PROGR
    }
}

/// `inchi_strtol` + the cast to `AT_NUMB`
fn num16(s: &[u8], p: &mut usize, base: i64) -> u16 {
    let (v, q) = inchi_strtol(s, *p, base);
    *p = q;
    v as u16
}

/// `ParseSegmentSp3s`: the `/s` layer (`s` starts at the `s`)
pub fn parse_segment_sp3s(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    sv: &mut [[i32; 2]; 2],
    ppn_num_components: [i32; 2],
    state: i32,
) -> i32 {
    let b_iso = state == IST_MOBILE_H_ISO_SP3_S || state == IST_FIXED_H_ISO_SP3_S;
    if !b_iso && state != IST_MOBILE_H_SP3_S && state != IST_FIXED_H_SP3_S {
        return RI_ERR_PROGR;
    }
    if at(s, 0) != b's' {
        return 0;
    }
    let p_start = 1usize;
    let n_num = ppn_num_components[b_mobile_h];
    if strchr(s, p_start, b'/').is_some() {
        return RI_ERR_SYNTAX;
    }
    let p_end = p_start + strlen(s, p_start);
    if p_end == p_start {
        sv[b_mobile_h][b_iso as usize] = NO_VALUE_INT;
        for i in 0..n_num.max(0) as usize {
            let Some(c) = p_inp_inchi[b_mobile_h].get_mut(i) else {
                return RI_ERR_SYNTAX;
            };
            let st = stereo_slot(c, b_iso).get_or_insert_with(Default::default);
            st.n_comp_inv2_abs = 0;
            if st.n_number_of_stereo_centers != 0 {
                return RI_ERR_SYNTAX;
            }
            if st.t_parity.is_empty() {
                st.t_parity = vec![0; 1];
            }
            if st.n_number.is_empty() {
                st.n_number = vec![0; 1];
            }
        }
        n_num + 1
    } else {
        let (val, q) = strtol(s, p_start);
        if q == p_end && (1..=3).contains(&(val as i32)) {
            sv[b_mobile_h][b_iso as usize] = val as i32;
            n_num + 1
        } else {
            RI_ERR_SYNTAX
        }
    }
}

/// `bIsSp3LayerNotEmpty`
pub fn is_sp3_layer_not_empty(
    p_inp_inchi: &[Vec<InChI>; 2],
    b_mobile_h: usize,
    b_iso: bool,
    n_num_components: i32,
) -> i32 {
    let mut n = 0;
    for i in 0..n_num_components.max(0) as usize {
        let Some(c) = p_inp_inchi[b_mobile_h].get(i) else {
            break;
        };
        if c.b_deleted != 0 || c.n_number_of_atoms == 0 {
            continue;
        }
        let st = if b_iso {
            c.stereo_isotopic.as_ref()
        } else {
            c.stereo.as_ref()
        };
        if st.is_some_and(|s| {
            s.n_number_of_stereo_centers > 0 && !s.n_number.is_empty() && !s.t_parity.is_empty()
        }) {
            n += 1;
        }
    }
    n
}

/// `ParseSegmentSp3m`: the `/m` layer
pub fn parse_segment_sp3m(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    ppn_num_components: [i32; 2],
    state: i32,
) -> i32 {
    let b_iso = state == IST_MOBILE_H_ISO_SP3_M || state == IST_FIXED_H_ISO_SP3_M;
    if !b_iso && state != IST_MOBILE_H_SP3_M && state != IST_FIXED_H_SP3_M {
        return RI_ERR_PROGR;
    }
    let mut n_num = ppn_num_components[b_mobile_h];
    if at(s, 0) != b'm' {
        // /m is missing: inherit from a preceding stereo layer
        let (mut from_m, mut from_iso): (i32, i32) = (-1, -1);
        if b_mobile_h == TAUT_YES && !b_iso {
            return 0;
        }
        if b_mobile_h == TAUT_NON && !b_iso {
            from_m = TAUT_YES as i32;
            from_iso = 0;
        } else if b_mobile_h == TAUT_YES && b_iso {
            from_m = TAUT_YES as i32;
            from_iso = 0;
        } else if b_mobile_h == TAUT_NON && b_iso {
            from_m = TAUT_NON as i32;
            from_iso = 0;
            if is_sp3_layer_not_empty(p_inp_inchi, TAUT_NON, false, ppn_num_components[TAUT_NON])
                == 0
                && is_sp3_layer_not_empty(
                    p_inp_inchi,
                    TAUT_YES,
                    false,
                    ppn_num_components[TAUT_YES],
                ) == 0
            {
                from_m = TAUT_YES as i32;
                from_iso = 1;
            }
        }
        if from_m < 0 || from_iso < 0 {
            return RI_ERR_PROGR;
        }
        if is_sp3_layer_not_empty(
            p_inp_inchi,
            from_m as usize,
            from_iso != 0,
            ppn_num_components[from_m as usize],
        ) == 0
        {
            if (from_m == 0 && from_iso != 0) || (from_m != 0 && from_iso == 0) {
                from_m = TAUT_YES as i32;
                from_iso = 0;
                if is_sp3_layer_not_empty(
                    p_inp_inchi,
                    from_m as usize,
                    false,
                    ppn_num_components[from_m as usize],
                ) == 0
                {
                    return 0;
                }
            } else {
                return 0;
            }
        }
        n_num = ppn_num_components[b_mobile_h].min(ppn_num_components[from_m as usize]);
        for i in 0..n_num.max(0) as usize {
            let (Some(f), Some(t)) = (
                p_inp_inchi[from_m as usize].get(i),
                p_inp_inchi[b_mobile_h].get(i),
            ) else {
                break;
            };
            if f.n_number_of_atoms > 0
                && f.b_deleted == 0
                && t.n_number_of_atoms > 0
                && t.b_deleted == 0
            {
                let sf = if from_iso != 0 {
                    f.stereo_isotopic.as_ref()
                } else {
                    f.stereo.as_ref()
                };
                let inv = sf.map(|x| x.n_comp_inv2_abs);
                let t = &mut p_inp_inchi[b_mobile_h][i];
                if let (Some(inv), Some(st)) = (inv, stereo_slot(t, b_iso).as_mut()) {
                    st.n_comp_inv2_abs = inv;
                }
            }
        }
        return 0;
    }

    let p_start = 1usize;
    if strchr(s, p_start, b'/').is_some() {
        return RI_ERR_SYNTAX;
    }
    let p_end = p_start + strlen(s, p_start);
    let mut p = p_start;
    if p_end == p_start {
        for i in 0..n_num.max(0) as usize {
            let Some(c) = p_inp_inchi[b_mobile_h].get_mut(i) else {
                return RI_ERR_SYNTAX;
            };
            let st = stereo_slot(c, b_iso).get_or_insert_with(Default::default);
            st.n_comp_inv2_abs = NO_VALUE_INT;
            if st.t_parity.is_empty() {
                st.t_parity = vec![0; 1];
            }
            if st.n_number.is_empty() {
                st.n_number = vec![0; 1];
            }
        }
        n_num + 1
    } else {
        let mut i_comp = 0i32;
        while p < p_end && i_comp < n_num {
            let Some(c) = p_inp_inchi[b_mobile_h].get_mut(i_comp as usize) else {
                return RI_ERR_SYNTAX;
            };
            let slot = stereo_slot(c, b_iso);
            if at(s, p) != b'.' && slot.is_none() {
                *slot = Some(InChIStereo::default());
            }
            match at(s, p) {
                b'1' => slot.as_mut().unwrap().n_comp_inv2_abs = -1,
                b'0' => slot.as_mut().unwrap().n_comp_inv2_abs = 1,
                b'.' => {
                    if let Some(st) = slot.as_mut() {
                        st.n_comp_inv2_abs = 0;
                    }
                }
                _ => return RI_ERR_SYNTAX,
            }
            i_comp += 1;
            p += 1;
        }
        if p != p_end || i_comp != n_num {
            return RI_ERR_SYNTAX;
        }
        n_num + 1
    }
}

/// the parsing of the number before a `mnMNe` abbreviation or a `*`
/// multiplier shared by `/b` and `/t`: `(val, q, if_cnd)`
fn leading_count(s: &[u8], p_start: usize) -> (i32, usize, bool) {
    if strtol(s, p_start).0 as i32 > 0 {
        let (v, q) = strtol(s, p_start);
        (v as i32, q, is_digit(at(s, p_start)))
    } else {
        (1, p_start, true)
    }
}

/// `INChI_Stereo` slot `i` of the sp3 arrays (they hold `len + 1` entries)
fn put_sp3(st: &mut InChIStereo, i: usize, parity: i32, n: u16) {
    if i < st.t_parity.len() {
        st.t_parity[i] = parity as i8;
    }
    if i < st.n_number.len() {
        st.n_number[i] = n;
    }
}

/// `SegmentSp3StoreStereoCenters`
fn segment_sp3_store_stereo_centers(
    pb_abc: i32,
    s: &[u8],
    p_start: usize,
    p_end: usize,
    n_atoms: i32,
    st: &mut InChIStereo,
) -> i32 {
    let mut p = p_start;
    let mut i_atom = 0usize;
    if pb_abc == 1 {
        while p < p_end {
            // (the C `base` stays 10 here, so a compressed list is rejected)
            let n1 = num16(s, &mut p, 10);
            if n1 == 0 {
                return RI_ERR_SYNTAX;
            }
            let (v, q) = strtol(s, p);
            p = q;
            let parity = v as i32;
            if !(AB_MIN_KNOWN_PARITY..=AB_MAX_KNOWN_PARITY).contains(&parity) {
                return RI_ERR_SYNTAX;
            }
            if i32::from(n1) > n_atoms {
                return RI_ERR_SYNTAX;
            }
            put_sp3(st, i_atom, parity, n1);
            if i_atom > 0 && st.n_number[i_atom - 1] >= n1 {
                return RI_ERR_SYNTAX;
            }
            i_atom += 1;
        }
    } else {
        while p < p_end {
            let (v, q) = strtol(s, p);
            let n1 = v as u16;
            if i32::from(n1) > MAX_ATOMS {
                return RI_ERR_SYNTAX;
            }
            let Some(r) = strchr_set(PARITY_TYPE, at(s, q)) else {
                return RI_ERR_SYNTAX;
            };
            p = q + 1;
            put_sp3(st, i_atom, r as i32 + 1, n1);
            if i_atom > 0 && st.n_number[i_atom - 1] >= n1 {
                return RI_ERR_SYNTAX;
            }
            i_atom += 1;
            p += usize::from(at(s, p) == b',');
        }
    }
    st.n_number_of_stereo_centers = i_atom as i32;
    st.n_comp_inv2_abs = NO_VALUE_INT;
    if p != p_end {
        return RI_ERR_SYNTAX;
    }
    0
}

/// the `(from, bIsoTo, bIsoFrom)` of a `/t` or `/b` abbreviation, or the
/// syntax error; `from` is 0 for the layer itself, 1 for the other tautomer
/// layer
fn stereo_abbreviation(
    b_mobile_h: usize,
    state: i32,
    q: u8,
    sp2: bool,
) -> Result<(usize, i32, i32), i32> {
    let (st_iso, st_fixed, st_fixed_iso) = if sp2 {
        (IST_MOBILE_H_ISO_SP2, IST_FIXED_H_SP2, IST_FIXED_H_ISO_SP2)
    } else {
        (IST_MOBILE_H_ISO_SP3, IST_FIXED_H_SP3, IST_FIXED_H_ISO_SP3)
    };
    let alt = 1 - b_mobile_h;
    if b_mobile_h == TAUT_YES {
        if state == st_iso {
            return match q {
                b'm' => Ok((b_mobile_h, 1, 0)),
                b'e' => Ok((b_mobile_h, 1, -1)),
                _ => Err(RI_ERR_SYNTAX),
            };
        }
        Err(RI_ERR_SYNTAX)
    } else if state == st_fixed {
        if q == b'm' {
            return Ok((alt, 0, 0));
        }
        // the sp2 layer has no other abbreviation here; the sp3 one has `e`
        if !sp2 && q == b'e' {
            return Ok((b_mobile_h, 1, -1));
        }
        Err(RI_ERR_SYNTAX)
    } else if state == st_fixed_iso {
        match q {
            b'm' => Ok((alt, 1, 0)),
            b'M' => Ok((alt, 1, 1)),
            b'n' => Ok((b_mobile_h, 1, 0)),
            b'e' => Ok((b_mobile_h, 1, -1)),
            _ => Err(RI_ERR_SYNTAX),
        }
    } else {
        Err(RI_ERR_SYNTAX)
    }
}

/// `ParseSegmentSp3`: the `/t` layer
pub fn parse_segment_sp3(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    ppn_num_components: [i32; 2],
    state: i32,
    pb_abc: &mut i32,
) -> i32 {
    let b_iso = state == IST_MOBILE_H_ISO_SP3 || state == IST_FIXED_H_ISO_SP3;
    if !b_iso && state != IST_MOBILE_H_SP3 && state != IST_FIXED_H_SP3 {
        return RI_ERR_PROGR;
    }
    if at(s, 0) != b't' {
        return 0;
    }
    let mut p_start = 1usize;
    let mut i_comp = 0i32;
    let n_num = ppn_num_components[b_mobile_h];

    if at(s, p_start) == 0 {
        // SegmentSp3CreateEmpty
        for i in 0..n_num.max(0) as usize {
            let Some(c) = p_inp_inchi[b_mobile_h].get_mut(i) else {
                return RI_ERR_SYNTAX;
            };
            let st = stereo_slot(c, b_iso).get_or_insert_with(Default::default);
            if st.b_parity.is_empty() {
                st.b_parity = vec![0; 1];
            }
            if st.n_bond_atom1.is_empty() {
                st.n_bond_atom1 = vec![0; 1];
            }
            if st.n_bond_atom2.is_empty() {
                st.n_bond_atom2 = vec![0; 1];
            }
            st.n_comp_inv2_abs = NO_VALUE_INT;
        }
        return n_num + 1;
    }

    loop {
        let p_end = strchr(s, p_start, b';').unwrap_or(p_start + strlen(s, p_start));
        let (val, q, if_cnd) = leading_count(s, p_start);
        let mpy_component: i32;
        let mut p: usize;
        if if_cnd && strchr_set(MULT_TYPE, at(s, q)).is_some() && q + 1 == p_end {
            // SegmentSp3ProcessAbbreviation
            if i_comp.wrapping_add(val) > n_num {
                return RI_ERR_SYNTAX;
            }
            let (from, iso_to, iso_from) =
                match stereo_abbreviation(b_mobile_h, state, at(s, q), false) {
                    Ok(x) => x,
                    Err(e) => return e,
                };
            for i in 0..val {
                let idx = (i_comp + i) as usize;
                let src = p_inp_inchi[from].get(idx).cloned();
                let Some(t) = p_inp_inchi[b_mobile_h].get_mut(idx) else {
                    return RI_ERR_SYNTAX;
                };
                let mut r = copy_segment(t, src.as_ref(), CPY_SP3, iso_to, iso_from);
                if r == 0 {
                    r = RI_ERR_SYNTAX;
                }
                if r < 0 {
                    if r == RI_ERR_SYNTAX {
                        return RI_ERR_SYNTAX;
                    }
                    // other errors are ignored by the caller
                    break;
                }
                if iso_from >= 0 {
                    if let Some(st) = stereo_slot(t, iso_to != 0).as_mut() {
                        st.n_comp_inv2_abs = NO_VALUE_INT;
                    }
                }
            }
            mpy_component = val;
        } else {
            if let Some(star) = strchr(s, p_start, b'*').filter(|&x| x < p_end) {
                let (v, q) = strtol(s, p_start);
                mpy_component = v as i32;
                if star != q {
                    return RI_ERR_SYNTAX;
                }
                p = star + 1;
            } else {
                mpy_component = 1;
                p = p_start;
            }
            if i_comp.wrapping_add(mpy_component) > n_num {
                return RI_ERR_SYNTAX;
            }
            p_start = p;
            if p < p_end && *pb_abc == -1 {
                *pb_abc = i32::from(upper(s, p));
            }
            let base = if *pb_abc == 1 { ALPHA_BASE } else { 10 };
            let n_atoms = match p_inp_inchi[b_mobile_h].get(i_comp.max(0) as usize) {
                Some(c) if i_comp >= 0 => c.n_number_of_atoms,
                _ => 0,
            };
            // pass 1: the number of the stereo centres
            let mut i_atom = 0usize;
            p = p_start;
            if *pb_abc == 1 {
                while p < p_end {
                    let n1 = num16(s, &mut p, base);
                    let ok = n1 != 0 && {
                        let (v, q) = strtol(s, p);
                        p = q;
                        (AB_MIN_KNOWN_PARITY..=AB_MAX_KNOWN_PARITY).contains(&(v as i32))
                    };
                    if !ok || i32::from(n1) > n_atoms {
                        return RI_ERR_SYNTAX;
                    }
                    i_atom += 1;
                }
            } else {
                while p < p_end {
                    let (v, q) = strtol(s, p);
                    let n1 = v as u16;
                    p = q + 1;
                    if n1 == 0
                        || i32::from(n1) > n_atoms
                        || strchr_set(PARITY_TYPE, at(s, q)).is_none()
                    {
                        return RI_ERR_SYNTAX;
                    }
                    i_atom += 1;
                    p += usize::from(at(s, p) == b',');
                }
            }
            if p != p_end {
                return RI_ERR_SYNTAX;
            }
            let len = i_atom;
            if i_comp > n_num - 1 || i_comp < 0 {
                return RI_ERR_SYNTAX;
            }
            if n_atoms <= 0 || n_atoms > MAX_ATOMS {
                return RI_ERR_SYNTAX;
            }
            let Some(c) = p_inp_inchi[b_mobile_h].get_mut(i_comp as usize) else {
                return RI_ERR_SYNTAX;
            };
            let st = stereo_slot(c, b_iso).get_or_insert_with(Default::default);
            if !st.t_parity.is_empty()
                || st.n_number_of_stereo_centers != 0
                || !st.n_number.is_empty()
            {
                return RI_ERR_SYNTAX;
            }
            st.t_parity = vec![0; len + 1];
            st.n_number = vec![0; len + 1];
            if segment_sp3_store_stereo_centers(*pb_abc, s, p_start, p_end, n_atoms, st)
                == RI_ERR_SYNTAX
            {
                return RI_ERR_SYNTAX;
            }
            // SegmentSp3CopyMultiplierCovered
            for i in 1..mpy_component {
                let src = p_inp_inchi[b_mobile_h].get(i_comp as usize).cloned();
                let Some(t) = p_inp_inchi[b_mobile_h].get_mut((i_comp + i) as usize) else {
                    return RI_ERR_SYNTAX;
                };
                for ty in [CPY_SP3, CPY_SP3_M] {
                    let mut r =
                        copy_segment(t, src.as_ref(), ty, i32::from(b_iso), i32::from(b_iso));
                    if r == 0 {
                        r = RI_ERR_SYNTAX;
                    }
                    if r < 0 {
                        if r == RI_ERR_SYNTAX {
                            return RI_ERR_SYNTAX;
                        }
                        break;
                    }
                }
            }
        }
        i_comp = i_comp.wrapping_add(mpy_component);
        if at(s, p_end) != 0 {
            p_start = p_end + 1;
        } else {
            break;
        }
    }
    if n_num != i_comp {
        return RI_ERR_SYNTAX;
    }
    i_comp + 1
}

/// `ParseSegmentSp2`: the `/b` layer
pub fn parse_segment_sp2(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    ppn_num_components: [i32; 2],
    state: i32,
    pb_abc: &mut i32,
) -> i32 {
    let b_iso = state == IST_MOBILE_H_ISO_SP2 || state == IST_FIXED_H_ISO_SP2;
    if !b_iso && state != IST_MOBILE_H_SP2 && state != IST_FIXED_H_SP2 {
        return RI_ERR_PROGR;
    }
    if at(s, 0) != b'b' {
        return 0;
    }
    let mut p_start = 1usize;
    let mut i_comp = 0i32;
    let n_num = ppn_num_components[b_mobile_h];

    if at(s, p_start) == 0 {
        // an empty sp2 segment means no sp2
        for i in 0..n_num.max(0) as usize {
            let Some(c) = p_inp_inchi[b_mobile_h].get_mut(i) else {
                return RI_ERR_SYNTAX;
            };
            if let Some(st) = stereo_slot(c, b_iso).as_ref() {
                if !st.b_parity.is_empty()
                    || st.n_number_of_stereo_bonds != 0
                    || !st.n_bond_atom1.is_empty()
                    || !st.n_bond_atom2.is_empty()
                {
                    return RI_ERR_SYNTAX;
                }
            }
            let r = copy_segment(c, None, CPY_SP2, i32::from(b_iso), -1);
            if r < 0 {
                return r;
            }
        }
        return n_num + 1;
    }

    loop {
        let p_end = strchr(s, p_start, b';').unwrap_or(p_start + strlen(s, p_start));
        let (val, q, if_cnd) = leading_count(s, p_start);
        let mpy_component: i32;
        if if_cnd && strchr_set(MULT_TYPE, at(s, q)).is_some() && q + 1 == p_end {
            if i_comp.wrapping_add(val) > n_num {
                return RI_ERR_SYNTAX;
            }
            let (from, iso_to, iso_from) =
                match stereo_abbreviation(b_mobile_h, state, at(s, q), true) {
                    Ok(x) => x,
                    Err(e) => return e,
                };
            for i in 0..val {
                let idx = i_comp + i;
                if idx > n_num - 1 || idx < 0 {
                    return RI_ERR_SYNTAX;
                }
                let idx = idx as usize;
                let src = p_inp_inchi[from].get(idx).cloned();
                match &src {
                    None => return RI_ERR_SYNTAX,
                    Some(f) if f.n_number_of_atoms <= 0 || f.n_number_of_atoms > MAX_ATOMS => {
                        return RI_ERR_SYNTAX;
                    }
                    _ => {}
                }
                let Some(t) = p_inp_inchi[b_mobile_h].get_mut(idx) else {
                    return RI_ERR_SYNTAX;
                };
                let mut r = copy_segment(t, src.as_ref(), CPY_SP2, iso_to, iso_from);
                if r == 0 {
                    r = RI_ERR_SYNTAX;
                }
                if r < 0 {
                    return r;
                }
            }
            mpy_component = val;
        } else {
            let mut p: usize;
            if let Some(star) = strchr(s, p_start, b'*').filter(|&x| x < p_end) {
                let (v, q) = strtol(s, p_start);
                mpy_component = v as i32;
                if !(0..=MAX_ATOMS).contains(&mpy_component) {
                    return RI_ERR_SYNTAX;
                }
                if star != q {
                    return RI_ERR_SYNTAX;
                }
                p = star + 1;
            } else {
                mpy_component = 1;
                p = p_start;
            }
            if i_comp.wrapping_add(mpy_component) > n_num {
                return RI_ERR_SYNTAX;
            }
            p_start = p;
            if p < p_end && *pb_abc == -1 {
                *pb_abc = i32::from(upper(s, p));
            }
            let base = if *pb_abc == 1 { ALPHA_BASE } else { 10 };
            let n_atoms = match p_inp_inchi[b_mobile_h].get(i_comp.max(0) as usize) {
                Some(c) if i_comp >= 0 => c.n_number_of_atoms,
                _ => 0,
            };
            // pass 1: the number of the stereo bonds
            let mut i_bond = 0usize;
            p = p_start;
            if *pb_abc == 1 {
                while p < p_end {
                    let a1 = num16(s, &mut p, base);
                    let a2 = if a1 != 0 { num16(s, &mut p, base) } else { 0 };
                    let ok = a1 != 0 && a2 != 0 && {
                        let (v, q) = strtol(s, p);
                        p = q;
                        (AB_MIN_KNOWN_PARITY..=AB_MAX_KNOWN_PARITY).contains(&(v as i32))
                    };
                    if !ok || a1 <= a2 || i32::from(a1) > n_atoms {
                        return RI_ERR_SYNTAX;
                    }
                    i_bond += 1;
                }
            } else {
                while p < p_end {
                    let (v1, q) = strtol(s, p);
                    if i32::from(v1 as u16) > MAX_ATOMS {
                        return RI_ERR_SYNTAX;
                    }
                    if at(s, q) != b'-' {
                        return RI_ERR_SYNTAX;
                    }
                    let (v2, q) = strtol(s, q + 1);
                    if i32::from(v2 as u16) > MAX_ATOMS {
                        return RI_ERR_SYNTAX;
                    }
                    let (a1, a2) = (v1 as u16, v2 as u16);
                    if a1 == 0
                        || a2 == 0
                        || a1 <= a2
                        || i32::from(a1) > n_atoms
                        || strchr_set(PARITY_TYPE, at(s, q)).is_none()
                    {
                        return RI_ERR_SYNTAX;
                    }
                    p = q + 1;
                    i_bond += 1;
                    p += usize::from(at(s, p) == b',');
                }
            }
            if p != p_end {
                return RI_ERR_SYNTAX;
            }
            let len = i_bond;
            if i_comp > n_num - 1 || i_comp < 0 {
                return RI_ERR_SYNTAX;
            }
            if n_atoms <= 0 || n_atoms > MAX_ATOMS {
                return RI_ERR_SYNTAX;
            }
            let Some(c) = p_inp_inchi[b_mobile_h].get_mut(i_comp as usize) else {
                return RI_ERR_SYNTAX;
            };
            let st = stereo_slot(c, b_iso).get_or_insert_with(Default::default);
            if !st.b_parity.is_empty()
                || st.n_number_of_stereo_bonds != 0
                || !st.n_bond_atom1.is_empty()
                || !st.n_bond_atom2.is_empty()
            {
                return RI_ERR_SYNTAX;
            }
            st.b_parity = vec![0; len + 1];
            st.n_bond_atom1 = vec![0; len + 1];
            st.n_bond_atom2 = vec![0; len + 1];
            // pass 2: store the stereo bonds
            let mut i_bond = 0usize;
            p = p_start;
            if *pb_abc == 1 {
                while p < p_end {
                    let a1 = num16(s, &mut p, base);
                    let a2 = if a1 != 0 { num16(s, &mut p, base) } else { 0 };
                    let ok = a1 != 0 && a2 != 0 && {
                        let (v, q) = strtol(s, p);
                        p = q;
                        st.b_parity[i_bond] = v as i8;
                        (AB_MIN_KNOWN_PARITY..=AB_MAX_KNOWN_PARITY).contains(&(v as i32))
                    };
                    if !ok {
                        return RI_ERR_SYNTAX;
                    }
                    st.n_bond_atom1[i_bond] = a1;
                    st.n_bond_atom2[i_bond] = a2;
                    if i_bond > 0
                        && !(st.n_bond_atom1[i_bond - 1] < a1
                            || (st.n_bond_atom1[i_bond - 1] == a1
                                && st.n_bond_atom2[i_bond - 1] < a2))
                    {
                        return RI_ERR_SYNTAX;
                    }
                    i_bond += 1;
                }
            } else {
                while p < p_end {
                    let (v1, q) = strtol(s, p);
                    if i32::from(v1 as u16) > MAX_ATOMS || at(s, q) != b'-' {
                        return RI_ERR_SYNTAX;
                    }
                    let (v2, q) = strtol(s, q + 1);
                    if i32::from(v2 as u16) > MAX_ATOMS {
                        return RI_ERR_SYNTAX;
                    }
                    let Some(r) = strchr_set(PARITY_TYPE, at(s, q)) else {
                        return RI_ERR_SYNTAX;
                    };
                    p = q + 1;
                    let (a1, a2) = (v1 as u16, v2 as u16);
                    st.b_parity[i_bond] = (r + 1) as i8;
                    st.n_bond_atom1[i_bond] = a1;
                    st.n_bond_atom2[i_bond] = a2;
                    if i_bond > 0
                        && !(st.n_bond_atom1[i_bond - 1] < a1
                            || (st.n_bond_atom1[i_bond - 1] == a1
                                && st.n_bond_atom2[i_bond - 1] < a2))
                    {
                        return RI_ERR_SYNTAX;
                    }
                    i_bond += 1;
                    p += usize::from(at(s, p) == b',');
                }
            }
            st.n_number_of_stereo_bonds = i_bond as i32;
            if p != p_end {
                return RI_ERR_SYNTAX;
            }
            // multiplier
            for i in 1..mpy_component {
                let src = p_inp_inchi[b_mobile_h].get(i_comp as usize).cloned();
                let Some(t) = p_inp_inchi[b_mobile_h].get_mut((i_comp + i) as usize) else {
                    return RI_ERR_SYNTAX;
                };
                let r = copy_segment(t, src.as_ref(), CPY_SP2, i32::from(b_iso), i32::from(b_iso));
                if r < 0 {
                    return r;
                }
            }
        }
        i_comp = i_comp.wrapping_add(mpy_component);
        if at(s, p_end) != 0 {
            p_start = p_end + 1;
        } else {
            break;
        }
    }
    if n_num != i_comp {
        return RI_ERR_SYNTAX;
    }
    i_comp + 1
}

/// `ParseSegmentIsoExchgH`: the `/i/h` layer (isotopic exchangeable H)
pub fn parse_segment_iso_exchg_h(
    s: &[u8],
    b_mobile_h: usize,
    n_num_protons: &mut [RemProtons],
    state: i32,
    pb_abc: &mut i32,
) -> i32 {
    const ABC_H: &[u8] = b"hdt";
    if at(s, 0) != b'h' {
        return 0;
    }
    let p_start = 1usize;
    if !(b_mobile_h == TAUT_YES && state == IST_MOBILE_H_ISO_EXCH_H) {
        return RI_ERR_PROGR;
    }
    if strchr(s, p_start, b';').is_some() {
        return RI_ERR_SYNTAX;
    }
    let p_end = p_start + strlen(s, p_start);
    let mut p = p_start;
    if p < p_end && *pb_abc == -1 {
        *pb_abc = i32::from(is_digit(at(s, p)));
    }
    if *pb_abc == 1 {
        let mut i_prev = ABC_H.len() as i32 + 1;
        while p < p_end {
            let (v, q) = strtol(s, p);
            let num = v as i32;
            if !(0..=MAX_ATOMS).contains(&num) {
                return RI_ERR_SYNTAX;
            }
            if 0 >= num || p == q || q >= p_end {
                return RI_ERR_SYNTAX;
            }
            match ABC_H.iter().position(|&c| c == at(s, q)).map(|i| i as i32) {
                Some(i) if i < i_prev && i < 3 => {
                    n_num_protons[b_mobile_h].n_num_removed_isotopic_h[i as usize] = num as i16;
                    p = q + 1;
                    i_prev = i;
                }
                _ => return RI_ERR_SYNTAX,
            }
        }
    } else {
        for (c, idx) in [(b'T', 2usize), (b'D', 1), (b'H', 0)] {
            if at(s, p) == c {
                n_num_protons[b_mobile_h].n_num_removed_isotopic_h[idx] = 1;
                p += 1;
                if is_digit(at(s, p)) {
                    let (v, q) = strtol(s, p);
                    n_num_protons[b_mobile_h].n_num_removed_isotopic_h[idx] = v as i16;
                    p = q;
                }
            }
        }
    }
    if p != p_end {
        return RI_ERR_SYNTAX;
    }
    1
}

/// `ParseSegmentPerm`: the `/o` layer (the transposition of the fixed-H
/// components)
pub fn parse_segment_perm(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    ppn_num_components: [i32; 2],
    state: i32,
    pb_abc: &mut i32,
) -> i32 {
    if at(s, 0) != b'o' {
        return 0;
    }
    if p_inp_inchi[b_mobile_h].is_empty() {
        return RI_ERR_ALLOC;
    }
    let mut p_start = 1usize;
    let n_num = ppn_num_components[b_mobile_h];
    if !(b_mobile_h == TAUT_NON
        && (state == IST_FIXED_H_PERMUTATION || state == IST_FIXED_H_ISO_PERMUTATION))
    {
        return RI_ERR_PROGR;
    }
    if strchr(s, p_start, b';').is_some() {
        return RI_ERR_SYNTAX;
    }
    let p_end = p_start + strlen(s, p_start);
    let v = &mut p_inp_inchi[b_mobile_h];
    while p_start < p_end {
        if at(s, p_start) != b'(' {
            return RI_ERR_SYNTAX;
        }
        let p_perm_start = p_start + 1;
        let mut tmp = InChI::default();
        let p_perm_end = match strchr(s, p_perm_start, b')') {
            Some(e) if e != p_perm_start => e,
            _ => return RI_ERR_SYNTAX,
        };
        if p_perm_start < p_perm_end && *pb_abc == -1 {
            *pb_abc = i32::from(upper(s, p_perm_start));
        }
        let base = if *pb_abc == 1 { ALPHA_BASE } else { 10 };
        let mut p = p_perm_start;
        let (mut i_comp2, mut num_trans) = (0i32, 0);
        while p < p_perm_end {
            let (q, i_comp1);
            if *pb_abc == 1 {
                let (val, e) = inchi_strtol(s, p, base);
                q = e;
                i_comp1 = val as i32;
                if 0 >= i_comp1 || i_comp1 > n_num {
                    return RI_ERR_SYNTAX;
                }
            } else {
                if !is_digit(at(s, p)) {
                    return RI_ERR_SYNTAX;
                }
                let (val, e) = strtol(s, p);
                q = e;
                i_comp1 = val as i32;
                if i_comp1 < 1 || i_comp1 > n_num {
                    return RI_ERR_SYNTAX;
                }
            }
            if i_comp2 != 0 {
                v[(i_comp2 - 1) as usize] = v[(i_comp1 - 1) as usize].clone();
                num_trans += 1;
            } else {
                tmp = v[(i_comp1 - 1) as usize].clone();
            }
            i_comp2 = i_comp1;
            p = if *pb_abc == 1 {
                q
            } else {
                q + usize::from(at(s, q) == b',')
            };
        }
        v[(i_comp2 - 1) as usize] = tmp;
        if num_trans == 0 || p != p_perm_end {
            return RI_ERR_SYNTAX;
        }
        p_start = p + 1;
    }
    1
}

/// `ParseSegmentIsoAtoms`: the `/i` layer (isotopic atoms)
pub fn parse_segment_iso_atoms(
    s: &[u8],
    b_mobile_h: usize,
    p_inp_inchi: &mut [Vec<InChI>; 2],
    ppn_num_components: [i32; 2],
    state: i32,
    pb_abc: &mut i32,
) -> i32 {
    const PARITY: &[u8] = b"-+TDH";
    if at(s, 0) != b'i' {
        return 0;
    }
    let mut p_start = 1usize;
    let mut i_comp = 0i32;
    let n_num = ppn_num_components[b_mobile_h];
    if !((b_mobile_h == TAUT_YES && state == IST_MOBILE_H_ISO_ATOMS)
        || (b_mobile_h == TAUT_NON && state == IST_FIXED_H_ISO_ATOMS))
    {
        return RI_ERR_PROGR;
    }
    if at(s, p_start) == 0 {
        return n_num + 1;
    }
    loop {
        let p_end = strchr(s, p_start, b';').unwrap_or(p_start + strlen(s, p_start));
        let mut p: usize;
        let mpy_component: i32;
        if let Some(star) = strchr(s, p_start, b'*').filter(|&x| x < p_end) {
            let (v, q) = strtol(s, p_start);
            mpy_component = v as i32;
            if star != q {
                return RI_ERR_SYNTAX;
            }
            if i_comp.wrapping_add(mpy_component) > n_num {
                return RI_ERR_SYNTAX;
            }
            p = star + 1;
        } else {
            let (val, q, if_cnd) = leading_count(s, p_start);
            if if_cnd && strchr_set(MULT_TYPE, at(s, q)).is_some() && q + 1 == p_end {
                if i_comp.wrapping_add(val) > n_num {
                    return RI_ERR_SYNTAX;
                }
                let (from, iso_from) = match (b_mobile_h, at(s, q)) {
                    (TAUT_YES, _) => return RI_ERR_SYNTAX,
                    (_, b'm') => (1 - b_mobile_h, 0),
                    (_, b'e') => (b_mobile_h, -1),
                    _ => return RI_ERR_SYNTAX,
                };
                for i in 0..val {
                    let idx = i_comp.wrapping_add(i);
                    if idx > n_num || idx < 0 {
                        return RI_ERR_SYNTAX;
                    }
                    let src = p_inp_inchi[from].get(idx as usize).cloned();
                    let Some(t) = p_inp_inchi[b_mobile_h].get_mut(idx as usize) else {
                        return RI_ERR_SYNTAX;
                    };
                    let mut r = copy_segment(t, src.as_ref(), CPY_ISO_AT, 0, iso_from);
                    if r == 0 {
                        r = RI_ERR_SYNTAX;
                    }
                    if r < 0 {
                        return r;
                    }
                }
                i_comp = i_comp.wrapping_add(val);
                if at(s, p_end) != 0 {
                    p_start = p_end + 1;
                    continue;
                }
                break;
            }
            mpy_component = 1;
            p = p_start;
        }
        if i_comp > n_num - 1 || i_comp < 0 {
            return RI_ERR_SYNTAX;
        }
        let n_atoms = p_inp_inchi[b_mobile_h]
            .get(i_comp as usize)
            .map_or(0, |c| c.n_number_of_atoms);
        if n_atoms <= 0 || n_atoms > MAX_ATOMS {
            return RI_ERR_SYNTAX;
        }
        p_start = p;
        if !p_inp_inchi[b_mobile_h][i_comp as usize]
            .isotopic_atom
            .is_empty()
        {
            return RI_ERR_PROGR;
        }
        if p < p_end && *pb_abc == -1 {
            *pb_abc = i32::from(upper(s, p));
        }
        let base = if *pb_abc == 1 { ALPHA_BASE } else { 10 };
        // the 1st pass finds the number of atoms, the 2nd one stores them
        let mut len = 0usize;
        let mut allocated = false;
        loop {
            let mut atoms: Vec<InChIIsotopicAtom> = Vec::new();
            let mut i_atom = 0usize;
            p = p_start;
            while p < p_end {
                let mut iso = InChIIsotopicAtom::default();
                let n1;
                if *pb_abc == 1 {
                    n1 = num16(s, &mut p, base);
                    if i32::from(n1) > MAX_ATOMS || n1 == 0 || i32::from(n1) > n_atoms {
                        return RI_ERR_SYNTAX;
                    }
                    iso.n_atom_number = n1;
                    let (v, q) = strtol(s, p);
                    iso.n_iso_difference = v as i16;
                    if p == q {
                        return RI_ERR_SYNTAX;
                    }
                    p = q;
                    parse_iso_counts(s, &mut p, &mut iso, b"tdh");
                    if p > p_end || iso_is_zero(&iso) {
                        return RI_ERR_SYNTAX;
                    }
                } else {
                    let (v, q) = strtol(s, p);
                    n1 = v as u16;
                    if i32::from(n1) > MAX_ATOMS {
                        return RI_ERR_SYNTAX;
                    }
                    p = q;
                    if n1 == 0 || i32::from(n1) > n_atoms || strchr_set(PARITY, at(s, p)).is_none()
                    {
                        return RI_ERR_SYNTAX;
                    }
                    iso.n_atom_number = n1;
                    if at(s, p) == b'+' && is_digit(at(s, p + 1)) {
                        let (v, q) = strtol(s, p + 1);
                        iso.n_iso_difference = v as i16;
                        if iso.n_iso_difference >= 0 {
                            iso.n_iso_difference = iso.n_iso_difference.wrapping_add(1);
                        }
                        p = q;
                    } else if at(s, p) == b'-' && is_digit(at(s, p + 1)) {
                        let (v, q) = strtol(s, p + 1);
                        iso.n_iso_difference = (-i32::from(v as i16)) as i16;
                        if iso.n_iso_difference == 0 {
                            iso.n_iso_difference = 1;
                        }
                        p = q;
                    }
                    parse_iso_counts(s, &mut p, &mut iso, b"TDH");
                    if iso_is_zero(&iso) {
                        return RI_ERR_SYNTAX;
                    }
                    if p < p_end {
                        if at(s, p) == b',' {
                            p += 1;
                        } else {
                            return RI_ERR_SYNTAX;
                        }
                    }
                }
                atoms.push(iso);
                i_atom += 1;
            }
            if p != p_end {
                return RI_ERR_SYNTAX;
            }
            if !allocated {
                len = i_atom;
                allocated = true;
                // calloc(len + 1); the 2nd pass fills it
                let c = &mut p_inp_inchi[b_mobile_h][i_comp as usize];
                c.isotopic_atom = vec![InChIIsotopicAtom::default(); len + 1];
                continue;
            }
            if len != i_atom {
                return RI_ERR_PROGR;
            }
            let c = &mut p_inp_inchi[b_mobile_h][i_comp as usize];
            for (k, a) in atoms.into_iter().enumerate() {
                c.isotopic_atom[k] = a;
            }
            c.n_number_of_isotopic_atoms = len as i32;
            break;
        }
        // multiplier
        for i in 1..mpy_component {
            let src = p_inp_inchi[b_mobile_h].get(i_comp as usize).cloned();
            let Some(t) = p_inp_inchi[b_mobile_h].get_mut((i_comp + i) as usize) else {
                return RI_ERR_SYNTAX;
            };
            let mut r = copy_segment(t, src.as_ref(), CPY_ISO_AT, 0, 0);
            if r == 0 {
                r = RI_ERR_SYNTAX;
            }
            if r < 0 {
                return r;
            }
        }
        i_comp = i_comp.wrapping_add(mpy_component);
        if at(s, p_end) != 0 {
            p_start = p_end + 1;
        } else {
            break;
        }
    }
    if n_num != i_comp {
        return RI_ERR_SYNTAX;
    }
    i_comp + 1
}

fn iso_is_zero(a: &InChIIsotopicAtom) -> bool {
    a.n_iso_difference == 0 && a.n_num_t == 0 && a.n_num_d == 0 && a.n_num_h == 0
}

/// the optional `T[n]`, `D[n]`, `H[n]` (upper case, or lower case in the
/// compressed form) of an isotopic atom
fn parse_iso_counts(s: &[u8], p: &mut usize, iso: &mut InChIIsotopicAtom, letters: &[u8; 3]) {
    for (k, &c) in letters.iter().enumerate() {
        if at(s, *p) == c {
            let field = match k {
                0 => &mut iso.n_num_t,
                1 => &mut iso.n_num_d,
                _ => &mut iso.n_num_h,
            };
            *field = 1;
            *p += 1;
            if is_digit(at(s, *p)) {
                let (v, q) = strtol(s, *p);
                *field = v as i16;
                *p = q;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The segment state machine: `ReadInChILine`

pub const INCHI_BAS: usize = 0;
pub const INCHI_REC: usize = 1;
const SEG_END: i32 = b'/' as i32;
/// `INCHI_TOKEN`
pub const INCHI_TOKEN: &[u8] = b"/\n\r\t\\";
pub const SAVE_OPT_KET: u8 = 0x10;
pub const SAVE_OPT_15T: u8 = 0x20;

/// `INCHI_INP_EOL`
fn inchi_inp_eol(c: i32) -> bool {
    c == 10 || c == 13 || c == 9
}

/// `strstr(s, pat) != NULL` over a zero-terminated string
fn contains_c(s: &[u8], pat: &[u8]) -> bool {
    let n = strlen(s, 0);
    s[..n].windows(pat.len()).any(|w| w == pat)
}

/// `ParseSegmentPolymer`: only the layer test; the polymer units are not
/// ported, so a `/z` layer is a syntax error
fn parse_segment_polymer_stub(s: &[u8]) -> i32 {
    if at(s, 0) != b'z' {
        return 0;
    }
    RI_ERR_SYNTAX
}

/// the state of the `InChI` structures read by `ReadInChILine`
#[derive(Debug, Clone, Default)]
pub struct InpInChI {
    pub p_inp_inchi: [[Vec<InChI>; 2]; 2],
    pub n_num_components: [[i32; 2]; 2],
    pub n_num_protons: [[RemProtons; 2]; 2],
    pub s: [[[i32; 2]; 2]; 2],
    /// `CompareInchiFlags[INCHI_NUM][TAUT_NUM]`: `INCHIDIFF` bits of the
    /// comparison with the restored structures
    pub compare_inchi_flags: [[u64; 2]; 2],
}

/// `ReadInChILine`: reads the InChI string of the stream into `inp_inchi`;
/// returns `RI_ERR_EOL`/`RI_ERR_EOF` at the end of the line/data, or the
/// error with the failed state in `p_state`. (A `/z` layer is a syntax
/// error.) `p_str` receives a `Structure ...` line.
#[allow(clippy::too_many_arguments)]
pub fn read_inchi_line(
    inp: &mut InchiIostream,
    line: &mut SegmLine,
    p_str: &mut Option<Vec<u8>>,
    p_state: &mut i32,
    d: &mut InpInChI,
    b_std_format: &mut i32,
    input_has_save_opt: &mut i32,
    input_save_opt_bits: &mut u8,
) -> i32 {
    let mut ret: i32;
    let mut b_mobile_h = TAUT_YES;
    let mut b_reconn = INCHI_BAS;
    let mut state: i32;
    let mut prev_state: i32 = -1;
    let mut b_abc = -1i32;
    let mut na_total: i32;
    let mut nb_total = 0i32;
    let (mut let1, mut let2) = (0u8, 0u8);
    const A2P: &[u8] = b"ABCDEFGHIJKLMNOP";
    *p_state = 0;
    let mut c: i32;

    loop {
        // next_line
        c = if line.len != 0 {
            line.c
        } else {
            n_get_inchi_segment(inp, line, INCHI_TOKEN)
        };
        if c == RI_ERR_EOF && line.len == 0 && at(&line.str, 0) == 0 {
            ret = c;
            if ret >= 0 || c == RI_ERR_EOF || c == RI_ERR_EOL {
                line.len = 0;
            }
            return ret;
        }
        if line.len == 0 || (c != SEG_END && c != RI_ERR_EOF) || !contains_c(&line.str, b"InChI=1")
        {
            if line.str.starts_with(b"Structure") {
                let n = strlen(&line.str, 0);
                *p_str = Some(line.str[..n].to_vec());
                *line = SegmLine::default();
                while c != 0 && !inchi_inp_eol(c) {
                    c = get_inchi_char(inp);
                }
                continue;
            }
            while c != RI_ERR_EOF && !inchi_inp_eol(c) {
                c = get_inchi_char(inp);
            }
            line.len = 0;
            continue;
        }
        break;
    }

    // a standard InChI?
    *b_std_format = i32::from(line.len == 8 && at(&line.str, 7) == b'S');

    state = IST_MOBILE_H_FORMULA;
    ret = 1;
    'main: loop {
        if ret < 0 {
            *p_state = prev_state;
            break;
        }
        prev_state = state
            + if b_reconn != 0 {
                IST_HAPPENED_IN_RECMET
            } else {
                0
            };
        if 0 < ret {
            if c != RI_ERR_EOF && c != SEG_END {
                // abnormal reading result unless we got backslash-SaveOpt
                if c == i32::from(b'\\') {
                    *input_has_save_opt = 1;
                }
                let mut k = 0usize;
                while c != RI_ERR_EOF && !inchi_inp_eol(c) {
                    c = get_inchi_char(inp);
                    k += 1;
                    if k == 1 {
                        let1 = c as u8;
                    } else if k == 2 {
                        let2 = c as u8;
                    }
                }
                if k != 3 {
                    *input_has_save_opt = 0;
                } else if (b'A'..=b'D').contains(&let2) {
                    *input_has_save_opt = 0;
                    *input_save_opt_bits = 0;
                    for (kk, &a) in A2P.iter().enumerate() {
                        if a == let1 {
                            *input_save_opt_bits = kk as u8;
                            *input_has_save_opt = 1;
                            break;
                        }
                    }
                    if *input_has_save_opt != 0 {
                        if let2 == b'B' || let2 == b'D' {
                            *input_save_opt_bits |= SAVE_OPT_15T;
                        }
                        if let2 == b'C' || let2 == b'D' {
                            *input_save_opt_bits |= SAVE_OPT_KET;
                        }
                    }
                }
                ret = if c == RI_ERR_EOF {
                    RI_ERR_EOF
                } else {
                    RI_ERR_EOL
                };
                line.len = 0;
                line.c = ret;
                break 'main;
            }
            if c == RI_ERR_EOF {
                ret = RI_ERR_EOF;
                break 'main;
            }
            if c == SEG_END {
                c = n_get_inchi_segment(inp, line, INCHI_TOKEN);
            }
            if c < 0 {
                // exit_error
                break 'main;
            }
            if line.len == 0 {
                ret = RI_ERR_EOL;
                break 'main;
            }
        }
        let rc = b_reconn;
        let mob = b_mobile_h;
        let str_ = line.str.clone();
        let s = &str_[..];
        match state {
            IST_MOBILE_H_FORMULA => {
                b_mobile_h = TAUT_YES;
                na_total = state;
                ret = parse_segment_formula(
                    s,
                    TAUT_YES,
                    &mut d.p_inp_inchi[rc],
                    &mut d.n_num_components[rc],
                    &mut na_total,
                );
                state = IST_MOBILE_H_CONNECTIONS;
            }
            IST_MOBILE_H_CONNECTIONS => {
                ret = parse_segment_connections(
                    s,
                    &mut d.p_inp_inchi[rc][mob],
                    &mut d.n_num_components[rc][mob],
                    &mut b_abc,
                    &mut nb_total,
                );
                state = IST_MOBILE_H;
            }
            IST_MOBILE_H | IST_FIXED_H => {
                ret = parse_segment_mobile_h(
                    s,
                    mob,
                    &mut d.p_inp_inchi[rc],
                    d.n_num_components[rc],
                    &mut b_abc,
                );
                state = if state == IST_MOBILE_H {
                    IST_MOBILE_H_CHARGE
                } else {
                    IST_FIXED_H_CHARGE
                };
            }
            IST_MOBILE_H_CHARGE | IST_FIXED_H_CHARGE => {
                ret = parse_segment_charge(s, mob, &mut d.p_inp_inchi[rc], d.n_num_components[rc]);
                state = if state == IST_MOBILE_H_CHARGE {
                    IST_MOBILE_H_PROTONS
                } else {
                    IST_FIXED_H_SP2
                };
            }
            IST_MOBILE_H_PROTONS => {
                ret = parse_segment_protons(s, mob, &mut d.n_num_protons[rc]);
                state = IST_MOBILE_H_POLYMER;
            }
            IST_MOBILE_H_POLYMER => {
                ret = parse_segment_polymer_stub(s);
                state = IST_MOBILE_H_SP2;
            }
            IST_MOBILE_H_SP2 | IST_MOBILE_H_ISO_SP2 | IST_FIXED_H_SP2 | IST_FIXED_H_ISO_SP2 => {
                ret = parse_segment_sp2(
                    s,
                    mob,
                    &mut d.p_inp_inchi[rc],
                    d.n_num_components[rc],
                    state,
                    &mut b_abc,
                );
                state = match state {
                    IST_MOBILE_H_SP2 => IST_MOBILE_H_SP3,
                    IST_MOBILE_H_ISO_SP2 => IST_MOBILE_H_ISO_SP3,
                    IST_FIXED_H_SP2 => IST_FIXED_H_SP3,
                    _ => IST_FIXED_H_ISO_SP3,
                };
            }
            IST_MOBILE_H_SP3 | IST_MOBILE_H_ISO_SP3 | IST_FIXED_H_SP3 | IST_FIXED_H_ISO_SP3 => {
                ret = parse_segment_sp3(
                    s,
                    mob,
                    &mut d.p_inp_inchi[rc],
                    d.n_num_components[rc],
                    state,
                    &mut b_abc,
                );
                state = match state {
                    IST_MOBILE_H_SP3 => IST_MOBILE_H_SP3_M,
                    IST_MOBILE_H_ISO_SP3 => IST_MOBILE_H_ISO_SP3_M,
                    IST_FIXED_H_SP3 => IST_FIXED_H_SP3_M,
                    _ => IST_FIXED_H_ISO_SP3_M,
                };
            }
            IST_MOBILE_H_SP3_M
            | IST_MOBILE_H_ISO_SP3_M
            | IST_FIXED_H_SP3_M
            | IST_FIXED_H_ISO_SP3_M => {
                ret = parse_segment_sp3m(
                    s,
                    mob,
                    &mut d.p_inp_inchi[rc],
                    d.n_num_components[rc],
                    state,
                );
                state = match state {
                    IST_MOBILE_H_SP3_M => IST_MOBILE_H_SP3_S,
                    IST_MOBILE_H_ISO_SP3_M => IST_MOBILE_H_ISO_SP3_S,
                    IST_FIXED_H_SP3_M => IST_FIXED_H_SP3_S,
                    _ => IST_FIXED_H_ISO_SP3_S,
                };
            }
            IST_MOBILE_H_SP3_S
            | IST_MOBILE_H_ISO_SP3_S
            | IST_FIXED_H_SP3_S
            | IST_FIXED_H_ISO_SP3_S => {
                ret = parse_segment_sp3s(
                    s,
                    mob,
                    &mut d.p_inp_inchi[rc],
                    &mut d.s[rc],
                    d.n_num_components[rc],
                    state,
                );
                state = match state {
                    IST_MOBILE_H_SP3_S => IST_MOBILE_H_ISO_LAYER_FORK,
                    IST_MOBILE_H_ISO_SP3_S => IST_FIXED_H_LAYER_FORK,
                    IST_FIXED_H_SP3_S => IST_FIXED_H_PERMUTATION,
                    _ => IST_FIXED_H_ISO_PERMUTATION,
                };
            }
            IST_MOBILE_H_ISO_LAYER_FORK
            | IST_FIXED_H_LAYER_FORK
            | IST_FIXED_H_ISO_LAYER_FORK
            | IST_RECONNECTED_LAYER_FORK => {
                // find the layer type after M, MI, F or FI
                ret = 0;
                let first = at(s, 0);
                let mut next = state;
                let mut keep_first = false; // `state == ISO_ATOMS`: the layer letter stays
                match (state, first) {
                    (IST_MOBILE_H_ISO_LAYER_FORK, b'i') => {
                        next = IST_MOBILE_H_ISO_ATOMS;
                        keep_first = true;
                    }
                    (IST_MOBILE_H_ISO_LAYER_FORK, b'f') | (IST_FIXED_H_LAYER_FORK, b'f') => {
                        next = IST_FIXED_H_FORMULA;
                    }
                    (IST_MOBILE_H_ISO_LAYER_FORK, b'r')
                    | (IST_FIXED_H_LAYER_FORK, b'r')
                    | (IST_FIXED_H_ISO_LAYER_FORK, b'r')
                    | (IST_RECONNECTED_LAYER_FORK, b'r') => {
                        next = IST_RECONNECTED_FORMULA;
                    }
                    (IST_FIXED_H_ISO_LAYER_FORK, b'i') => {
                        next = IST_FIXED_H_ISO_ATOMS;
                        keep_first = true;
                    }
                    _ => ret = RI_ERR_SYNTAX,
                }
                state = next;
                if inchi_inp_eol(c) && ret == 0 && at(s, 1) == 0 {
                    prev_state = state
                        + if b_reconn != 0 {
                            IST_HAPPENED_IN_RECMET
                        } else {
                            0
                        };
                    ret = RI_ERR_SYNTAX;
                } else if ret == 0 && !keep_first {
                    let len = strlen(s, 0);
                    if len > 1 {
                        // memmove(str, str + 1, len)
                        line.str.remove(0);
                    } else {
                        ret = 1;
                    }
                }
            }
            IST_MOBILE_H_ISO_ATOMS | IST_FIXED_H_ISO_ATOMS => {
                ret = parse_segment_iso_atoms(
                    s,
                    mob,
                    &mut d.p_inp_inchi[rc],
                    d.n_num_components[rc],
                    state,
                    &mut b_abc,
                );
                state = if state == IST_MOBILE_H_ISO_ATOMS {
                    IST_MOBILE_H_ISO_EXCH_H
                } else {
                    IST_FIXED_H_ISO_SP2
                };
            }
            IST_MOBILE_H_ISO_EXCH_H => {
                ret =
                    parse_segment_iso_exchg_h(s, mob, &mut d.n_num_protons[rc], state, &mut b_abc);
                state = IST_MOBILE_H_ISO_SP2;
            }
            IST_FIXED_H_FORMULA => {
                b_mobile_h = TAUT_NON;
                na_total = state;
                ret = parse_segment_formula(
                    s,
                    TAUT_NON,
                    &mut d.p_inp_inchi[rc],
                    &mut d.n_num_components[rc],
                    &mut na_total,
                );
                state = IST_FIXED_H;
            }
            IST_FIXED_H_PERMUTATION | IST_FIXED_H_ISO_PERMUTATION => {
                ret = parse_segment_perm(
                    s,
                    mob,
                    &mut d.p_inp_inchi[rc],
                    d.n_num_components[rc],
                    state,
                    &mut b_abc,
                );
                state = if state == IST_FIXED_H_PERMUTATION {
                    IST_FIXED_H_ISO_LAYER_FORK
                } else {
                    IST_RECONNECTED_LAYER_FORK
                };
            }
            IST_RECONNECTED_FORMULA => {
                if b_reconn == INCHI_REC {
                    // the reconnected layer may appear only once
                    ret = RI_ERR_SYNTAX;
                } else {
                    b_reconn = INCHI_REC;
                    b_mobile_h = TAUT_YES;
                    state = IST_MOBILE_H_FORMULA;
                }
            }
            _ => {}
        }
        if c < 0 {
            break;
        }
    }
    if ret >= 0 || c == RI_ERR_EOF || c == RI_ERR_EOL {
        line.len = 0;
    }
    ret
}

// ---------------------------------------------------------------------------
// `InChILine2Data`: the post-processing of the read layers (omitted layers)

const NUM_H_ISOTOPES: usize = 3;
const REL_RAC: u64 = INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO;
/// `REQ_MODE_RELATIVE_STEREO`, `REQ_MODE_RACEMIC_STEREO`
pub const REQ_MODE_RELATIVE_STEREO: u64 = 0x0200;
pub const REQ_MODE_RACEMIC_STEREO: u64 = 0x0400;

/// `GetInChIIsoH`: the numbers of the isotopic H (H, D, T) of a component
pub fn get_inchi_iso_h(c: &InChI) -> Result<[i32; 3], i32> {
    let mut n = [0i32; 3];
    for a in c
        .isotopic_atom
        .iter()
        .take(c.n_number_of_isotopic_atoms.max(0) as usize)
    {
        if a.n_iso_difference > 0 && i32::from(a.n_iso_difference) <= NUM_H_ISOTOPES as i32 {
            if c.n_atom.is_empty()
                || a.n_atom_number == 0
                || i32::from(a.n_atom_number) > c.n_number_of_atoms
            {
                return Err(RI_ERR_PROGR);
            }
            if c.n_atom
                .get(a.n_atom_number as usize - 1)
                .copied()
                .unwrap_or(0)
                == 1
            {
                n[a.n_iso_difference as usize - 1] += 1;
            }
        }
        n[0] += i32::from(a.n_num_h);
        n[1] += i32::from(a.n_num_d);
        n[2] += i32::from(a.n_num_t);
    }
    Ok(n)
}

fn self_copy(c: &mut InChI, ty: i32, b_to: i32, b_from: i32) -> i32 {
    let f = c.clone();
    copy_segment(c, Some(&f), ty, b_to, b_from)
}

fn nsc(s: &Option<InChIStereo>) -> i32 {
    s.as_ref().map_or(0, |s| s.n_number_of_stereo_centers)
}
fn nsb(s: &Option<InChIStereo>) -> i32 {
    s.as_ref().map_or(0, |s| s.n_number_of_stereo_bonds)
}
fn inv(s: &Option<InChIStereo>) -> i32 {
    s.as_ref().map_or(0, |s| s.n_comp_inv2_abs)
}
fn has_t(s: &Option<InChIStereo>) -> bool {
    s.as_ref().is_some_and(|s| !s.t_parity.is_empty())
}
fn has_b(s: &Option<InChIStereo>) -> bool {
    s.as_ref().is_some_and(|s| !s.b_parity.is_empty())
}
fn set_inv(s: &mut Option<InChIStereo>, v: i32) {
    if let Some(s) = s.as_mut() {
        s.n_comp_inv2_abs = v;
    }
}
/// `!x || NO_VALUE_INT == x`
fn is_unset(v: i32) -> bool {
    v == 0 || v == NO_VALUE_INT
}

/// the part of `InChILine2Data` that restores the omitted stereo segments
/// of one component (`pInChI` of the Fixed-H or Mobile-H layer `j`);
/// `alt` is the component `k` of the other layer when it exists,
/// `alt_atoms` tells that it has atoms (`pInChI_Alt` is not NULL then)
#[allow(clippy::too_many_arguments)]
fn restore_omitted_stereo(
    j: usize,
    cur: &mut InChI,
    mut alt: Option<&mut InChI>,
    alt_atoms: bool,
    s: &[[i32; 2]; 2],
    n_flags: u64,
    n_flags_alt: u64,
    n_comps_j: i32,
    mh2fh_alt_inv: &mut i32,
    fh2ifh_alt_inv: &mut i32,
) -> i32 {
    let ja = 1 - j;
    let rr = |v: i32| (v & 6) != 0;

    if j == TAUT_YES {
        // (4) copy stereo: Mobile H --> isotopic Mobile H
        let b_iso = cur.n_number_of_isotopic_atoms != 0
            || (cur.stereo_isotopic.is_some()
                && nsc(&cur.stereo_isotopic) + nsb(&cur.stereo_isotopic) != 0)
            || (alt_atoms
                && alt
                    .as_ref()
                    .is_some_and(|a| a.n_number_of_isotopic_atoms != 0));
        if b_iso {
            if nsc(&cur.stereo) != 0
                && (cur.stereo_isotopic.is_none() || !has_t(&cur.stereo_isotopic))
            {
                let r = self_copy(cur, CPY_SP3, 1, 0);
                if r < 0 {
                    return r;
                }
                if is_unset(inv(&cur.stereo_isotopic)) {
                    let r = self_copy(cur, CPY_SP3_M, 1, 0);
                    if r < 0 {
                        return r;
                    }
                }
                if (n_flags & REL_RAC) != 0 {
                    if inv(&cur.stereo) == NO_VALUE_INT {
                        set_inv(&mut cur.stereo, if s[j][0] > 0 { 2 } else { 0 });
                    }
                    if inv(&cur.stereo_isotopic) == NO_VALUE_INT {
                        set_inv(&mut cur.stereo_isotopic, if s[j][1] > 0 { 2 } else { 0 });
                    }
                }
            } else if nsc(&cur.stereo) != 0
                && nsc(&cur.stereo_isotopic) != 0
                && inv(&cur.stereo) != 0
            {
                if (n_flags & REL_RAC) != 0
                    && inv(&cur.stereo) == NO_VALUE_INT
                    && inv(&cur.stereo_isotopic) == NO_VALUE_INT
                {
                    set_inv(&mut cur.stereo, if s[j][0] > 0 { 2 } else { 0 });
                    set_inv(&mut cur.stereo_isotopic, if s[j][1] > 0 { 2 } else { 0 });
                } else if is_unset(inv(&cur.stereo_isotopic)) {
                    let v = inv(&cur.stereo);
                    set_inv(&mut cur.stereo_isotopic, v);
                }
            }
        }
        if b_iso
            && nsb(&cur.stereo) != 0
            && (cur.stereo_isotopic.is_none() || !has_b(&cur.stereo_isotopic))
        {
            let r = self_copy(cur, CPY_SP2, 1, 0);
            if r < 0 {
                return r;
            }
        }
    }

    // (0) set nCompInv2Abs to Fixed-H
    if j == TAUT_NON && nsc(&cur.stereo) != 0 && inv(&cur.stereo) == NO_VALUE_INT {
        let alt_has_st = alt_atoms && alt.as_ref().is_some_and(|a| nsc(&a.stereo) != 0);
        if s[j][0] == 0 && s[ja][0] > 0 && alt_has_st {
            let a = alt.as_deref_mut().unwrap();
            if inv(&a.stereo) == 0
                && inv(&cur.stereo) == NO_VALUE_INT
                && a.stereo_isotopic.is_some()
                && inv(&a.stereo_isotopic) != 0
                && inv(&a.stereo_isotopic) != NO_VALUE_INT
                && s[j][0] == 0
                && s[j][1] == 0
                && s[ja][0] == 1
                && s[ja][1] == 1
            {
                let r = copy_segment(cur, Some(a), CPY_SP3_M, 1, 1);
                if r < 0 {
                    return r;
                }
            } else if inv(&a.stereo) != NO_VALUE_INT && s[ja][0] == 1 {
                let v = inv(&a.stereo);
                set_inv(&mut cur.stereo, v);
            } else if inv(&a.stereo) == NO_VALUE_INT
                && inv(&cur.stereo) == NO_VALUE_INT
                && s[j][0] == 0
                && s[j][1] == 0
                && rr(s[ja][0])
                && rr(s[ja][1])
            {
                let b_iso = cur.n_number_of_isotopic_atoms != 0
                    || (cur.stereo_isotopic.is_some()
                        && nsc(&cur.stereo_isotopic) + nsb(&cur.stereo_isotopic) != 0)
                    || a.n_number_of_isotopic_atoms != 0;
                if b_iso {
                    if a.stereo_isotopic.is_none() {
                        let r = self_copy(a, CPY_SP3_M, 1, -1);
                        if r < 0 {
                            return r;
                        }
                    }
                    set_inv(&mut a.stereo_isotopic, 2);
                    set_inv(&mut a.stereo, 0);
                    set_inv(&mut cur.stereo, NO_VALUE_INT + 1);
                } else {
                    set_inv(&mut cur.stereo, 2);
                    set_inv(&mut a.stereo, 2);
                }
            } else {
                set_inv(&mut cur.stereo, 2);
                set_inv(&mut a.stereo, 2);
            }
        } else if s[j][0] == 0
            && s[ja][0] == 0
            && (n_flags & REL_RAC) != 0
            && !(cur.stereo_isotopic.is_some() && nsc(&cur.stereo_isotopic) != 0)
            && !(alt_atoms
                && alt
                    .as_ref()
                    .is_some_and(|a| a.stereo_isotopic.is_some() && nsc(&a.stereo_isotopic) != 0))
        {
            set_inv(&mut cur.stereo, NO_VALUE_INT + 1);
        } else {
            set_inv(&mut cur.stereo, if s[j][0] > 0 { 2 } else { 0 });
        }
    }

    // (1) copy stereo: non-isotopic Fixed H --> isotopic Fixed H
    if j == TAUT_NON {
        let b_iso = cur.n_number_of_isotopic_atoms != 0
            || (cur.stereo_isotopic.is_some()
                && nsc(&cur.stereo_isotopic) + nsb(&cur.stereo_isotopic) != 0)
            || (alt_atoms
                && alt
                    .as_ref()
                    .is_some_and(|a| a.n_number_of_isotopic_atoms != 0));
        if b_iso {
            if nsc(&cur.stereo) != 0
                && (cur.stereo_isotopic.is_none() || !has_t(&cur.stereo_isotopic))
            {
                let r = self_copy(cur, CPY_SP3, 1, 0);
                if r < 0 {
                    return r;
                }
                if inv(&cur.stereo) == NO_VALUE_INT + 1 {
                    set_inv(&mut cur.stereo, 0);
                    set_inv(&mut cur.stereo_isotopic, 2);
                } else if inv(&cur.stereo_isotopic) == 0 {
                    let r = self_copy(cur, CPY_SP3_M, 1, 0);
                    if r < 0 {
                        return r;
                    }
                }
                if (n_flags & REL_RAC) != 0 {
                    if inv(&cur.stereo) == NO_VALUE_INT {
                        set_inv(&mut cur.stereo, if s[j][0] > 0 { 2 } else { 0 });
                    }
                    if inv(&cur.stereo_isotopic) == NO_VALUE_INT {
                        set_inv(&mut cur.stereo_isotopic, if s[j][1] > 0 { 2 } else { 0 });
                    }
                }
            } else if nsc(&cur.stereo) != 0
                && nsc(&cur.stereo_isotopic) != 0
                && inv(&cur.stereo) != 0
            {
                if (n_flags & REL_RAC) != 0
                    && inv(&cur.stereo) == NO_VALUE_INT
                    && inv(&cur.stereo_isotopic) == NO_VALUE_INT
                {
                    set_inv(&mut cur.stereo, if s[j][0] > 0 { 2 } else { 0 });
                    set_inv(&mut cur.stereo_isotopic, if s[j][1] > 0 { 2 } else { 0 });
                } else if is_unset(inv(&cur.stereo_isotopic)) {
                    let v = inv(&cur.stereo);
                    set_inv(&mut cur.stereo_isotopic, v);
                }
            }
        }
        if b_iso
            && nsb(&cur.stereo) != 0
            && (cur.stereo_isotopic.is_none() || !has_b(&cur.stereo_isotopic))
        {
            let r = self_copy(cur, CPY_SP2, 1, 0);
            if r < 0 {
                return r;
            }
        }
    }

    // copy stereo: Mobile H --> Fixed H
    if j == TAUT_NON {
        if let Some(a) = alt {
            let b_iso = cur.n_number_of_isotopic_atoms != 0
                || (cur.stereo_isotopic.is_some()
                    && nsc(&cur.stereo_isotopic) + nsb(&cur.stereo_isotopic) != 0)
                || a.n_number_of_isotopic_atoms != 0
                || (a.stereo_isotopic.is_some()
                    && nsc(&a.stereo_isotopic) + nsb(&a.stereo_isotopic) != 0);
            let b_no_t = !has_t(&cur.stereo);
            let b_no_m = cur.stereo.is_none() || inv(&cur.stereo) == NO_VALUE_INT;

            // (2) non-isotopic Mobile H => non-isotopic Fixed H
            if nsc(&a.stereo) != 0 && !has_t(&cur.stereo) {
                let mut b_have_t_stereo = if cur.stereo.is_some() {
                    nsc(&cur.stereo)
                } else {
                    1
                };
                if b_have_t_stereo < 1 && nsb(&cur.stereo) > 0 {
                    b_have_t_stereo = 1;
                }
                let r = copy_segment(cur, Some(a), CPY_SP3, 0, 0);
                if r < 0 {
                    return r;
                }
                if is_unset(inv(&cur.stereo)) {
                    let r = copy_segment(cur, Some(a), CPY_SP3_M, 0, 0);
                    if r < 0 {
                        return r;
                    }
                }
                if s[j][0] < 1 && b_have_t_stereo < 1 && nsc(&a.stereo) > 0 && s[ja][0] < 1 {
                    let mut all_uu = true;
                    if let Some(st) = a.stereo.as_ref() {
                        for kc in 0..st.n_number_of_stereo_centers.max(0) as usize {
                            let p = i32::from(st.t_parity.get(kc).copied().unwrap_or(0));
                            if p != 3 && p != 4 {
                                all_uu = false;
                                break;
                            }
                        }
                    }
                    if !all_uu {
                        if let Some(st) = cur.stereo.as_mut() {
                            st.n_number_of_stereo_centers = 0;
                        }
                    }
                }
                if inv(&cur.stereo) == NO_VALUE_INT && (n_flags_alt & REL_RAC) != 0 {
                    if s[ja][0] > 0 && s[j][0] > 0 {
                        *mh2fh_alt_inv += 1;
                        set_inv(&mut a.stereo, if *mh2fh_alt_inv % 2 != 0 { 2 } else { 0 });
                        set_inv(&mut cur.stereo, 2);
                    } else if n_comps_j > 1
                        && b_no_t
                        && b_no_m
                        && n_flags == 0
                        && !(a.stereo_isotopic.is_some() && has_t(&a.stereo_isotopic))
                        && !(cur.stereo_isotopic.is_some() && has_t(&cur.stereo_isotopic))
                        && s[j][0] == 0
                        && s[j][1] == 0
                        && rr(s[ja][0])
                        && rr(s[ja][1])
                    {
                        if b_iso {
                            if cur.stereo_isotopic.is_none() {
                                let r = copy_segment(cur, Some(a), CPY_SP3_M, 1, -1);
                                if r < 0 {
                                    return r;
                                }
                            }
                            if a.stereo_isotopic.is_none() {
                                let r = self_copy(a, CPY_SP3_M, 1, -1);
                                if r < 0 {
                                    return r;
                                }
                            }
                            set_inv(&mut cur.stereo_isotopic, 2);
                            set_inv(&mut cur.stereo, 0);
                            set_inv(&mut a.stereo_isotopic, 2);
                            set_inv(&mut a.stereo, 0);
                        } else if inv(&cur.stereo) == NO_VALUE_INT && inv(&a.stereo) == NO_VALUE_INT
                        {
                            set_inv(&mut cur.stereo, 2);
                            set_inv(&mut a.stereo, 2);
                        }
                    } else if (s[ja][0] > 0 || s[j][0] > 0) && s[j][0] >= 0 {
                        set_inv(&mut cur.stereo, 2);
                    } else if b_iso
                        && b_no_t
                        && b_no_m
                        && n_flags == 0
                        && !(a.stereo_isotopic.is_some() && has_t(&a.stereo_isotopic))
                        && !(cur.stereo_isotopic.is_some() && has_t(&cur.stereo_isotopic))
                        && s[ja][0] == 0
                        && s[j][0] == 0
                        && s[j][1] == 0
                        && rr(s[ja][1])
                    {
                        if cur.stereo_isotopic.is_none() {
                            let r = copy_segment(cur, Some(a), CPY_SP3_M, 1, -1);
                            if r < 0 {
                                return r;
                            }
                        }
                        set_inv(&mut cur.stereo_isotopic, 2);
                        set_inv(&mut cur.stereo, 0);
                    } else {
                        set_inv(&mut cur.stereo, 0);
                    }
                    if cur.n_flags & REL_RAC == 0 {
                        cur.n_flags |= (n_flags_alt | n_flags) & REL_RAC;
                    }
                }
            } else if nsc(&a.stereo) != 0
                && nsc(&cur.stereo) != 0
                && inv(&a.stereo) != 0
                && is_unset(inv(&cur.stereo))
                && (!(n_flags_alt != 0 && n_flags == 0) || NO_VALUE_INT == inv(&cur.stereo))
            {
                let v = inv(&a.stereo);
                set_inv(&mut cur.stereo, v);
            }

            // use the same rule to copy the stereobonds
            if nsb(&a.stereo) != 0 && !has_b(&cur.stereo) {
                let r = copy_segment(cur, Some(a), CPY_SP2, 0, 0);
                if r < 0 {
                    return r;
                }
            }
            // (3) isotopic Mobile H -> isotopic Fixed H
            if b_iso {
                if !has_t(&cur.stereo)
                    && !has_t(&a.stereo)
                    && nsc(&a.stereo_isotopic) != 0
                    && (cur.stereo_isotopic.is_none() || !has_t(&cur.stereo_isotopic))
                {
                    let r = copy_segment(cur, Some(a), CPY_SP3, 1, 1);
                    if r < 0 {
                        return r;
                    }
                    if is_unset(inv(&cur.stereo_isotopic)) {
                        let r = copy_segment(cur, Some(a), CPY_SP3_M, 1, 1);
                        if r < 0 {
                            return r;
                        }
                    }
                    if inv(&cur.stereo_isotopic) == NO_VALUE_INT && (n_flags_alt & REL_RAC) != 0 {
                        set_inv(&mut cur.stereo_isotopic, if s[ja][1] > 0 { 2 } else { 0 });
                        if cur.n_flags & REL_RAC == 0 {
                            cur.n_flags |= n_flags_alt & REL_RAC;
                        }
                    }
                } else if !has_t(&cur.stereo)
                    && !has_t(&a.stereo)
                    && nsc(&a.stereo_isotopic) != 0
                    && nsc(&cur.stereo_isotopic) != 0
                    && inv(&a.stereo_isotopic) != 0
                    && is_unset(inv(&cur.stereo_isotopic))
                {
                    let r = copy_segment(cur, Some(a), CPY_SP3_M, 1, 1);
                    if r < 0 {
                        return r;
                    }
                }
                // use the same rule to copy the stereobonds
                if !has_b(&cur.stereo)
                    && !has_b(&a.stereo)
                    && nsb(&a.stereo_isotopic) != 0
                    && (cur.stereo_isotopic.is_none() || !has_b(&cur.stereo_isotopic))
                {
                    let r = copy_segment(cur, Some(a), CPY_SP2, 1, 1);
                    if r < 0 {
                        return r;
                    }
                }

                // (4) Copy Fixed-H -> isotopic Fixed-H
                if nsc(&cur.stereo) != 0
                    && !has_t(&a.stereo_isotopic)
                    && !has_t(&cur.stereo_isotopic)
                {
                    if a.stereo_isotopic.is_some()
                        && inv(&a.stereo_isotopic) != 0
                        && b_no_t
                        && NO_VALUE_INT != inv(&a.stereo_isotopic)
                        && !(cur.stereo_isotopic.is_some()
                            && NO_VALUE_INT != inv(&cur.stereo_isotopic))
                    {
                        let r = copy_segment(cur, Some(a), CPY_SP3_M, 1, 1);
                        if r < 0 {
                            return r;
                        }
                    }
                    let r = self_copy(cur, CPY_SP3, 1, 0);
                    if r < 0 {
                        return r;
                    }
                    if is_unset(inv(&cur.stereo_isotopic)) {
                        let r = self_copy(cur, CPY_SP3_M, 1, 0);
                        if r < 0 {
                            return r;
                        }
                    }
                    if inv(&cur.stereo_isotopic) == NO_VALUE_INT && (n_flags & REL_RAC) != 0 {
                        if s[j][0] > 0 && s[j][1] > 0 {
                            set_inv(&mut cur.stereo_isotopic, 2);
                            *fh2ifh_alt_inv += 1;
                            set_inv(
                                &mut cur.stereo,
                                if *fh2ifh_alt_inv % 2 != 0 { 2 } else { 0 },
                            );
                        } else if (s[j][0] > 0 || s[j][1] > 0) && s[j][1] >= 0 {
                            set_inv(&mut cur.stereo_isotopic, 2);
                        } else {
                            set_inv(&mut cur.stereo_isotopic, 0);
                        }
                        if cur.n_flags & REL_RAC == 0 {
                            cur.n_flags |= n_flags & REL_RAC;
                        }
                    }
                } else if has_t(&cur.stereo)
                    && !has_t(&a.stereo_isotopic)
                    && nsc(&a.stereo_isotopic) != 0
                    && nsc(&cur.stereo_isotopic) != 0
                    && inv(&cur.stereo) != 0
                    && is_unset(inv(&cur.stereo_isotopic))
                {
                    let r = self_copy(cur, CPY_SP3_M, 1, 0);
                    if r < 0 {
                        return r;
                    }
                }
            }
            if b_iso
                && nsb(&cur.stereo) == 0
                && nsb(&a.stereo) == 0
                && nsb(&a.stereo_isotopic) != 0
                && (cur.stereo_isotopic.is_none() || !has_b(&cur.stereo_isotopic))
            {
                let r = copy_segment(cur, Some(a), CPY_SP2, 1, 1);
                if r < 0 {
                    return r;
                }
            }
        }
    }
    0
}

/// `ReadInChICoord` (AuxInfo): the AuxInfo line is not ported; the lines
/// that are not `AuxInfo=` are handled as in C
fn read_inchi_coord(inp: &mut InchiIostream, line: &mut SegmLine, p_state: &mut i32) -> i32 {
    *p_state = 0;
    let c = if line.len != 0 {
        line.c
    } else {
        n_get_inchi_segment(inp, line, INCHI_TOKEN)
    };
    if c == RI_ERR_EOF && line.len == 0 && at(&line.str, 0) == 0 {
        line.len = 0;
        return c;
    }
    if line.len == 0 || (c != SEG_END && c != RI_ERR_EOF && !inchi_inp_eol(c)) {
        *p_state = -1;
        line.len = 0;
        return RI_ERR_PROGR;
    }
    if !line.str.starts_with(b"AuxInfo=") {
        *p_state = -1;
        return c;
    }
    // AuxInfo is not supported yet
    *p_state = -1;
    line.len = 0;
    RI_ERR_SYNTAX
}

/// the stereo type of a layer pair: the maximum of the two `/s` values
fn stereo_type_flags(s: &[[i32; 2]; 2], j: usize) -> (u64, u64) {
    let ja = 1 - j;
    let (mut t, mut ta) = (0, 0);
    for b in 0..2 {
        if t == 0 || t < s[j][b] {
            t = s[j][b];
        }
        if ta == 0 || ta < s[ja][b] {
            ta = s[ja][b];
        }
    }
    let f = |t: i32| match t {
        2 => INCHI_FLAG_REL_STEREO,
        3 => INCHI_FLAG_RAC_STEREO,
        _ => 0,
    };
    (f(t), f(ta))
}

/// the part of `InChILine2Data` after `ReadInChILine`: adds the omitted
/// layers and checks the balances; returns `ret2` (0 or an error)
fn line2data_postprocess(p_state: &mut i32, d: &mut InpInChI) -> i32 {
    *p_state = IST_MATERIAL_BALANCE_ERROR;
    let mut num_elem: [[Option<Vec<i32>>; 2]; 2] = Default::default();
    let mut tot_charge = [[0i32; 2]; 2];
    let el_len = super::util::N_EL_DATA_LEN as usize;
    for i_inchi in 0..2usize {
        for j in 0..2usize {
            let ja = 1 - j;
            let (mut mh2fh_alt_inv, mut fh2ifh_alt_inv) = (0i32, 0i32);
            tot_charge[i_inchi][j] = 0;
            let s = d.s[i_inchi];
            let (n_flags, n_flags_alt) = stereo_type_flags(&s, j);
            let ncj = d.n_num_components[i_inchi][j];
            if ncj != 0 {
                num_elem[i_inchi][j] = Some(vec![0; el_len + 1]);
            }
            for k in 0..ncj.max(0) as usize {
                if d.p_inp_inchi[i_inchi][j].is_empty() {
                    continue;
                }
                let n_alt = d.n_num_components[i_inchi][ja];
                let alt_exists = (k as i32) < n_alt && !d.p_inp_inchi[i_inchi][ja].is_empty();
                let alt_atoms = alt_exists && d.p_inp_inchi[i_inchi][ja][k].n_number_of_atoms != 0;
                let ncomp_yes = d.n_num_components[i_inchi][TAUT_YES];
                let ncomp_non = d.n_num_components[i_inchi][TAUT_NON];
                // the data of the Mobile-H component read before it is changed
                let mob_ct = if j == TAUT_NON
                    && (k as i32) < ncomp_yes
                    && !d.p_inp_inchi[i_inchi][TAUT_YES].is_empty()
                {
                    let m = &d.p_inp_inchi[i_inchi][TAUT_YES][k];
                    Some((m.n_conn_table.clone(), m.len_conn_table, m.n_total_charge))
                } else {
                    None
                };
                let yes_exists = !d.p_inp_inchi[i_inchi][TAUT_YES].is_empty();
                let pp = &mut d.p_inp_inchi[i_inchi];
                {
                    let cur = &mut pp[j][k];
                    if n_flags != 0 {
                        cur.n_flags |= n_flags;
                    } else if j == TAUT_NON && n_flags_alt != 0 {
                        cur.n_flags |= n_flags_alt;
                    }
                    // add empty immobile H (nNum_H) if it is missing
                    if cur.n_num_h.is_empty() {
                        cur.n_num_h = vec![0; cur.n_number_of_atoms.max(0) as usize + 1];
                    }
                    // add a single atom nConnTable if it is missing
                    if cur.n_conn_table.is_empty() && yes_exists {
                        match mob_ct.as_ref() {
                            Some((ct, len, _)) if j == TAUT_NON && !ct.is_empty() && *len > 0 => {
                                let mut v = vec![0u16; *len as usize + 1];
                                for (x, y) in v.iter_mut().zip(ct.iter()).take(*len as usize) {
                                    *x = *y;
                                }
                                cur.n_conn_table = v;
                                cur.len_conn_table = *len;
                            }
                            _ => {
                                if j == TAUT_YES && cur.n_number_of_atoms > 1 {
                                    *p_state = IST_MOBILE_H_CONNECTIONS
                                        + if i_inchi == INCHI_REC {
                                            IST_HAPPENED_IN_RECMET
                                        } else {
                                            0
                                        };
                                    return RI_ERR_SYNTAX;
                                }
                                cur.n_conn_table =
                                    vec![0; cur.n_number_of_atoms.max(0) as usize + 1];
                                cur.len_conn_table = 1;
                                cur.n_conn_table[0] = 1;
                            }
                        }
                    } else if !cur.n_conn_table.is_empty()
                        && cur.len_conn_table == 0
                        && cur.n_number_of_atoms == 1
                    {
                        cur.n_conn_table[0] = 1;
                        cur.len_conn_table = 1;
                    }
                    // copy charge: Mobile H --> Fixed H
                    if j == TAUT_NON && cur.n_total_charge == 0 {
                        if let Some((_, _, ch)) = mob_ct.as_ref() {
                            if *ch != 0 && *ch != NO_VALUE_INT {
                                cur.n_total_charge = *ch;
                            }
                        }
                    }
                    // Fixed H: add nNum_H_fixed to nNum_H
                    if j == TAUT_NON && !cur.n_num_h.is_empty() {
                        if let Some(fx) = cur.n_num_h_fixed.clone() {
                            for m in 0..cur.n_number_of_atoms.max(0) as usize {
                                let add = fx.get(m).copied().unwrap_or(0);
                                if let Some(x) = cur.n_num_h.get_mut(m) {
                                    *x = x.wrapping_add(add);
                                }
                            }
                        }
                    }
                }
                // copy isotopic atoms: Mobile H --> Fixed H
                if j == TAUT_YES
                    && pp[j][k].n_number_of_isotopic_atoms != 0
                    && (k as i32) < ncomp_non
                {
                    let from = pp[TAUT_YES][k].clone();
                    let alt = &mut pp[TAUT_NON][k];
                    if alt.n_number_of_isotopic_atoms == 0 {
                        let r = copy_segment(alt, Some(&from), CPY_ISO_AT, 0, 0);
                        if r < 0 {
                            return r;
                        }
                    }
                }
                // restore the omitted stereo segments
                {
                    let (a0, a1) = pp.split_at_mut(1);
                    let (cur_v, alt_v) = if j == 0 {
                        (&mut a0[0], &mut a1[0])
                    } else {
                        (&mut a1[0], &mut a0[0])
                    };
                    let cur = &mut cur_v[k];
                    let alt = if alt_exists { alt_v.get_mut(k) } else { None };
                    let r = restore_omitted_stereo(
                        j,
                        cur,
                        alt,
                        alt_atoms,
                        &s,
                        n_flags,
                        n_flags_alt,
                        ncj,
                        &mut mh2fh_alt_inv,
                        &mut fh2ifh_alt_inv,
                    );
                    if r < 0 {
                        return r;
                    }
                }
            }
        }

        // replace NO_VALUE_INT with zeroes in all Mobile & Fixed H components
        for j in 0..2usize {
            for k in 0..d.n_num_components[i_inchi][j].max(0) as usize {
                if d.p_inp_inchi[i_inchi][j].is_empty() {
                    continue;
                }
                let s = d.s[i_inchi];
                let c = &mut d.p_inp_inchi[i_inchi][j][k];
                if c.n_total_charge == NO_VALUE_INT {
                    c.n_total_charge = 0;
                }
                if c.stereo.is_some()
                    && c.stereo_isotopic.is_some()
                    && inv(&c.stereo_isotopic) == NO_VALUE_INT
                    && nsc(&c.stereo) != 0
                    && inv(&c.stereo) != NO_VALUE_INT
                {
                    let v = inv(&c.stereo);
                    set_inv(&mut c.stereo_isotopic, v);
                }
                // the special nCompInv2Abs = 2 forces /s2 or /s3 in the InChI output
                if c.stereo.is_some() && inv(&c.stereo) == NO_VALUE_INT {
                    if (c.n_flags & REL_RAC) != 0 && nsc(&c.stereo) != 0 {
                        set_inv(&mut c.stereo, if s[j][0] > 0 { 2 } else { 0 });
                    } else {
                        set_inv(&mut c.stereo, 0);
                    }
                }
                if c.stereo_isotopic.is_some() && inv(&c.stereo_isotopic) == NO_VALUE_INT {
                    if (c.n_flags & REL_RAC) != 0 && nsc(&c.stereo_isotopic) != 0 {
                        set_inv(&mut c.stereo_isotopic, if s[j][1] > 0 { 2 } else { 0 });
                    } else {
                        set_inv(&mut c.stereo_isotopic, 0);
                    }
                }
                if (c.stereo.is_some() && inv(&c.stereo) == NO_VALUE_INT)
                    || (c.stereo_isotopic.is_some() && inv(&c.stereo_isotopic) == NO_VALUE_INT)
                {
                    return RI_ERR_PROGR;
                }
                if c.b_deleted == 0 && c.n_number_of_atoms != 0 {
                    tot_charge[i_inchi][j] += c.n_total_charge;
                    for m in 0..c.n_number_of_atoms.max(0) as usize {
                        let a = i32::from(c.n_atom.get(m).copied().unwrap_or(0));
                        if a < i32::from(super::util::EL_NUMBER_H) || a > super::util::N_EL_DATA_LEN
                        {
                            return RI_ERR_PROGR;
                        }
                    }
                    let nh = match get_inchi_num_h(c) {
                        Ok(v) => v,
                        Err(e) => return e,
                    };
                    if let Some(ne) = num_elem[i_inchi][j].as_mut() {
                        for m in 0..c.n_number_of_atoms.max(0) as usize {
                            let a = usize::from(c.n_atom.get(m).copied().unwrap_or(0));
                            if a > super::util::EL_NUMBER_H as usize {
                                ne[a] += 1;
                            }
                        }
                        ne[super::util::EL_NUMBER_H as usize] += nh;
                    }
                }
            }
        }
        for j in 0..2usize {
            for k in 0..d.n_num_components[i_inchi][j].max(0) as usize {
                if d.p_inp_inchi[i_inchi][j].is_empty() {
                    continue;
                }
                let c = &mut d.p_inp_inchi[i_inchi][j][k];
                if c.stereo.is_some() && nsc(&c.stereo) == 0 {
                    set_inv(&mut c.stereo, 0);
                }
                if c.stereo_isotopic.is_some() && nsc(&c.stereo_isotopic) == 0 {
                    set_inv(&mut c.stereo_isotopic, 0);
                }
            }
        }

        // (2008-04-10) missed repeating /s in FI after F for multi-component case
        if d.n_num_components[i_inchi][TAUT_NON] > 1
            && d.s[i_inchi][TAUT_YES][0] == 0
            && d.s[i_inchi][TAUT_YES][1] == 0
            && d.s[i_inchi][TAUT_NON][0] > 1
            && d.s[i_inchi][TAUT_NON][1] > 1
            && d.s[i_inchi][TAUT_NON][0] == d.s[i_inchi][TAUT_NON][1]
        {
            for k in 0..d.n_num_components[i_inchi][TAUT_NON].max(0) as usize {
                let c = &mut d.p_inp_inchi[i_inchi][TAUT_NON][k];
                if inv(&c.stereo) != 0 {
                    set_inv(&mut c.stereo, 0);
                    break;
                }
            }
        }
        if let Some(ne) = num_elem[i_inchi][TAUT_YES].as_mut() {
            let np = i32::from(d.n_num_protons[i_inchi][TAUT_YES].n_num_removed_protons);
            tot_charge[i_inchi][TAUT_YES] += np;
            ne[super::util::EL_NUMBER_H as usize] += np;
        }

        // count H and isotopic H in Mobile and Fixed H representations
        if d.n_num_components[i_inchi][TAUT_NON] != 0 {
            let mut n_tot_removed_protons = 0i32;
            let mut n_tot_removed_iso = [0i32; 3];
            let len2 =
                d.n_num_components[i_inchi][TAUT_YES].max(d.n_num_components[i_inchi][TAUT_NON]);
            for k in 0..len2.max(0) as usize {
                let mut exists = [false; 2];
                for j in 0..2 {
                    exists[j] = (k as i32) < d.n_num_components[i_inchi][j]
                        && d.p_inp_inchi[i_inchi][j]
                            .get(k)
                            .is_some_and(|c| c.n_number_of_atoms != 0 && c.b_deleted == 0);
                }
                if !exists[TAUT_NON] {
                    return RI_ERR_SYNTAX;
                }
                let mut formula_h = [0i32; 2];
                let mut num_h = [0i32; 2];
                let mut charge = [0i32; 2];
                let mut iso_h = [[0i32; 3]; 2];
                for j in 0..2 {
                    if !exists[j] {
                        continue;
                    }
                    let c = &d.p_inp_inchi[i_inchi][j][k];
                    formula_h[j] = match get_inchi_formula_num_h(c) {
                        Ok(v) => v,
                        Err(e) => return e,
                    };
                    num_h[j] = match get_inchi_num_h(c) {
                        Ok(v) => v,
                        Err(e) => return e,
                    };
                    iso_h[j] = match get_inchi_iso_h(c) {
                        Ok(v) => v,
                        Err(e) => return e,
                    };
                    charge[j] = c.n_total_charge;
                }
                for j in 0..2 {
                    if exists[j] && formula_h[j] != num_h[j] {
                        return RI_ERR_SYNTAX;
                    }
                }
                let mut found_removed_iso = 0;
                let removed_charge = charge[TAUT_NON] - charge[TAUT_YES];
                let removed_h = num_h[TAUT_NON] - num_h[TAUT_YES];
                let mut removed_iso = [0i32; 3];
                for m in 0..3 {
                    removed_iso[m] = iso_h[TAUT_NON][m] - iso_h[TAUT_YES][m];
                    found_removed_iso += i32::from(removed_iso[m] != 0);
                }
                if removed_charge != removed_h {
                    return RI_ERR_SYNTAX;
                }
                if removed_charge != 0 || found_removed_iso != 0 {
                    let np = &mut d.n_num_protons[i_inchi][TAUT_YES];
                    if np.p_num_protons.is_none() {
                        np.p_num_protons = Some(vec![[0; 4]; len2.max(0) as usize]);
                    }
                    let e = &mut np.p_num_protons.as_mut().unwrap()[k];
                    e[0] = removed_h as i16;
                    n_tot_removed_protons += removed_h;
                    for m in 0..3 {
                        e[1 + m] = removed_iso[m] as i16;
                        n_tot_removed_iso[m] += removed_iso[m];
                    }
                    // make sure the Mobile-H InChI has nTautomer
                    if exists[TAUT_YES] {
                        let c = &mut d.p_inp_inchi[i_inchi][TAUT_YES][k];
                        if c.len_tautomer == 0 {
                            c.len_tautomer = 1;
                        }
                        if c.n_tautomer.is_empty() {
                            c.n_tautomer = vec![0; c.len_tautomer.max(0) as usize];
                        }
                    }
                }
            }
            let np = &d.n_num_protons[i_inchi][TAUT_YES];
            if np.p_num_protons.is_some() {
                if n_tot_removed_protons != i32::from(np.n_num_removed_protons) {
                    return RI_ERR_SYNTAX;
                }
                for m in 0..3 {
                    if n_tot_removed_iso[m] != i32::from(np.n_num_removed_isotopic_h[m]) {
                        return RI_ERR_SYNTAX;
                    }
                }
            }
        }

        // make the Mobile H and Fixed H InChI arrays have the same length
        let (nm, nf) = (
            d.n_num_components[i_inchi][TAUT_YES],
            d.n_num_components[i_inchi][TAUT_NON],
        );
        let (j, len1, len2) = if nm < nf {
            (TAUT_YES, nm, nf)
        } else if nm > nf {
            (TAUT_NON, nf, nm)
        } else {
            (0, 0, 0)
        };
        if len2 != 0 && len1 != 0 {
            let pp = &mut d.p_inp_inchi[i_inchi];
            pp[j].truncate(len1 as usize);
            pp[j].resize(len2 as usize, InChI::default());
            d.n_num_components[i_inchi][j] = len2;
            for l in len1 as usize..len2 as usize {
                if j == TAUT_YES {
                    let r = n_fill_out_proton_mobile_h(&mut pp[TAUT_YES][l]);
                    if r < 0 {
                        return r;
                    }
                    let from = pp[TAUT_NON][l].clone();
                    let r = n_proton_copy_isotopic_info(&mut pp[TAUT_YES][l], &from);
                    if r < 0 {
                        return r;
                    }
                } else {
                    pp[l][0].b_deleted = 1;
                }
            }
        }
    }

    // check balances
    for i_inchi in 0..2usize {
        for i in i_inchi..2usize {
            for j in 0..2usize {
                for k in j..2usize {
                    if (i_inchi != i || j != k)
                        && num_elem[i_inchi][j].is_some()
                        && num_elem[i][k].is_some()
                    {
                        if tot_charge[i_inchi][j] != tot_charge[i][k] {
                            return RI_ERR_SYNTAX;
                        }
                        let (a, b) = (
                            num_elem[i_inchi][j].as_ref().unwrap(),
                            num_elem[i][k].as_ref().unwrap(),
                        );
                        if a != b {
                            return RI_ERR_SYNTAX;
                        }
                    }
                }
            }
        }
    }
    0
}

/// `InChILine2Data`: reads an InChI line and adds the omitted layers.
/// Returns the result of `ReadInChILine`; `n_err` is the error of the
/// conversion (0 or a negative `RI_ERR_*`).
#[allow(clippy::too_many_arguments)]
pub fn inchi_line2data(
    inp: &mut InchiIostream,
    line: &mut SegmLine,
    p_str: &mut Option<Vec<u8>>,
    p_state: &mut i32,
    n_err: &mut i32,
    d: &mut InpInChI,
    b_read_coord: bool,
    n_mode: u64,
    b_std_format: &mut i32,
    input_has_save_opt: &mut i32,
    input_save_opt_bits: &mut u8,
) -> i32 {
    let target_stereo_type = if n_mode & REQ_MODE_RELATIVE_STEREO != 0 {
        2
    } else if n_mode & REQ_MODE_RACEMIC_STEREO != 0 {
        3
    } else {
        1
    };
    let ret = read_inchi_line(
        inp,
        line,
        p_str,
        p_state,
        d,
        b_std_format,
        input_has_save_opt,
        input_save_opt_bits,
    );
    if target_stereo_type > 1 {
        for i1 in 0..2 {
            for i2 in 0..2 {
                for i3 in 0..2 {
                    if d.s[i1][i2][i3] != 0 {
                        d.s[i1][i2][i3] = target_stereo_type;
                    }
                }
            }
        }
    }
    *n_err = 0;
    let ncomp = d.n_num_components[INCHI_BAS][TAUT_YES] + d.n_num_components[INCHI_BAS][TAUT_NON];
    if ret == RI_ERR_EOL && ncomp != 0 && b_read_coord {
        let mut st = 0;
        let _ = read_inchi_coord(inp, line, &mut st);
    }
    let ret2 = if (ret == RI_ERR_EOL || ret == RI_ERR_EOF) && ncomp != 0 {
        line2data_postprocess(p_state, d)
    } else {
        ret
    };
    *n_err = if ret2 < 0 && ret2 != RI_ERR_EOL {
        ret2
    } else {
        0
    };
    ret
}

// ---------------------------------------------------------------------------
// Preparing the read InChI for the structure restoration

/// `MODE_PIXH`
pub const MODE_PIXH_UNDEFINED: i32 = 0;
pub const MODE_PIXH_ADD_TO_FIRST: i32 = 1;
pub const MODE_PIXH_ADD_TO_EACH: i32 = 2;
pub const MODE_PIXH_ADD_A_PIXH_COMPONENT: i32 = 3;
pub const MODE_PIXH_KEEP_TOTALS: i32 = 4;

/// `GetNumNeighborsFromInchi`: the number of the neighbours (+1000 for a
/// tautomeric atom) or a negative `RI_ERR_*`
pub fn get_num_neighbors_from_inchi(inchi: &InChI, n_at_number: u16) -> i32 {
    let at_number = i32::from(n_at_number) - 1;
    let mut num_neigh = 0;
    let mut b_taut_atom = 0;
    let num_atoms = inchi.n_number_of_atoms;
    let ct = |i: usize| i32::from(inchi.n_conn_table.get(i).copied().unwrap_or(0));
    let mut n_vertex = ct(0) - 1;
    let mut i = 1;
    while i < inchi.len_conn_table {
        let n_neigh = ct(i as usize) - 1;
        if n_neigh < n_vertex {
            num_neigh += i32::from(at_number == n_vertex || at_number == n_neigh);
        } else {
            n_vertex = n_neigh;
            if n_vertex >= num_atoms {
                return RI_ERR_PROGR;
            }
        }
        i += 1;
    }
    if inchi.len_tautomer > 1 && !inchi.n_tautomer.is_empty() && inchi.n_tautomer[0] > 0 {
        let t = |j: i32| i32::from(inchi.n_tautomer.get(j as usize).copied().unwrap_or(0));
        let tot_len_tg = inchi.len_tautomer - T_GROUP_HDR_LEN_I * t(0) - 1;
        let mut j = 1;
        let mut i = 0;
        for _ in 0..t(0) {
            let mut len_tg = t(j);
            j += T_GROUP_HDR_LEN_I;
            len_tg -= T_GROUP_HDR_LEN_I - 1;
            while 0 < len_tg {
                len_tg -= 1;
                let taut_at_number = t(j) - 1;
                b_taut_atom += i32::from(taut_at_number == at_number);
                j += 1;
                i += 1;
            }
        }
        if i != tot_len_tg {
            return RI_ERR_PROGR;
        }
    }
    let num_h = if inchi.n_num_h.is_empty() {
        0
    } else {
        i32::from(
            inchi
                .n_num_h
                .get(at_number.max(0) as usize)
                .copied()
                .unwrap_or(0),
        )
    };
    let mut tot = num_neigh + num_h;
    if b_taut_atom != 0 {
        tot += 1000;
    }
    tot
}

const T_GROUP_HDR_LEN_I: i32 = super::ichicano::T_GROUP_HDR_LEN as i32;

/// `CountStereoTypes`: returns 0 (no InChI), 1 (no stereo), 2 (has stereo)
/// or a negative error; the counters are added to
pub fn count_stereo_types(
    inchi: &InChI,
    num_known_sb: &mut i32,
    num_known_sc: &mut i32,
    num_unk_und_sb: &mut i32,
    num_unk_und_sc: &mut i32,
    num_sc_piii: &mut i32,
    num_sc_asiii: &mut i32,
) -> i32 {
    if inchi.n_number_of_atoms == 0 || inchi.b_deleted != 0 {
        return 0;
    }
    let tot = |s: &Option<InChIStereo>| nsb(s) + nsc(s);
    let stereo = if inchi.stereo_isotopic.is_some() && tot(&inchi.stereo_isotopic) != 0 {
        inchi.stereo_isotopic.as_ref()
    } else if inchi.stereo.is_some() && tot(&inchi.stereo) != 0 {
        inchi.stereo.as_ref()
    } else {
        None
    };
    let Some(st) = stereo else {
        return 1;
    };
    for i in 0..st.n_number_of_stereo_bonds.max(0) as usize {
        if super::ichister::atom_parity_well_def(i32::from(
            st.b_parity.get(i).copied().unwrap_or(0),
        )) {
            *num_known_sb += 1;
        } else {
            *num_unk_und_sb += 1;
        }
    }
    for i in 0..st.n_number_of_stereo_centers.max(0) as usize {
        let at_number = st.n_number.get(i).copied().unwrap_or(0);
        if at_number == 0 || i32::from(at_number) > inchi.n_number_of_atoms {
            return RI_ERR_PROGR;
        }
        if super::ichister::atom_parity_well_def(i32::from(
            st.t_parity.get(i).copied().unwrap_or(0),
        )) {
            *num_known_sc += 1;
        } else {
            *num_unk_und_sc += 1;
        }
        let el = inchi
            .n_atom
            .get(at_number as usize - 1)
            .copied()
            .unwrap_or(0);
        if el != super::util::EL_NUMBER_P && el != super::util::EL_NUMBER_AS {
            continue;
        }
        let ret = get_num_neighbors_from_inchi(inchi, at_number);
        if ret < 0 {
            return ret;
        }
        if ret == 3 {
            *num_sc_piii += i32::from(el == super::util::EL_NUMBER_P);
            *num_sc_asiii += i32::from(el == super::util::EL_NUMBER_AS);
        }
    }
    2
}

/// `bInpInchiComponentExists`
pub fn b_inp_inchi_component_exists(
    d: &InpInChI,
    i_inchi: usize,
    b_mobile_h: usize,
    k: i32,
) -> bool {
    if i_inchi > 1 || b_mobile_h > 1 || k < 0 {
        return false;
    }
    k < d.n_num_components[i_inchi][b_mobile_h]
        && d.p_inp_inchi[i_inchi][b_mobile_h]
            .get(k as usize)
            .is_some_and(|c| c.n_number_of_atoms > 0 && c.b_deleted == 0)
}

/// `bInpInchiComponentDeleted`
pub fn b_inp_inchi_component_deleted(
    d: &InpInChI,
    i_inchi: usize,
    b_mobile_h: usize,
    k: i32,
) -> bool {
    if i_inchi > 1 || b_mobile_h > 1 || k < 0 {
        return false;
    }
    k < d.n_num_components[i_inchi][b_mobile_h]
        && d.p_inp_inchi[i_inchi][b_mobile_h]
            .get(k as usize)
            .is_some_and(|c| c.n_number_of_atoms > 0 && c.b_deleted != 0)
}

/// `bInChIHasReconnectedMetal`
pub fn b_inchi_has_reconnected_metal(c: &InChI) -> bool {
    if c.b_deleted == 0 && c.n_number_of_atoms != 0 && !c.n_atom.is_empty() {
        for i in 0..c.n_number_of_atoms.max(0) as usize {
            if super::util::is_el_a_metal(i32::from(c.n_atom.get(i).copied().unwrap_or(0)))
                && (c.n_number_of_atoms > 1 || c.n_num_h.first().copied().unwrap_or(0) != 0)
            {
                return true;
            }
        }
    }
    false
}

/// the result of `DetectInpInchiCreationOptions`
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct InpInchiOptions {
    pub has_reconnected: bool,
    pub has_metal: bool,
    pub has_fixed_h: bool,
    pub mode_flags_stereo: u64,
    pub taut_flags_stereo: u64,
}

/// `DetectInpInchiCreationOptions`; `Err` is the return code
pub fn detect_inp_inchi_creation_options(d: &InpInChI) -> Result<InpInchiOptions, i32> {
    use super::ichiparm::{
        REQ_MODE_ISO_STEREO, REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU, REQ_MODE_STEREO,
        TG_FLAG_ARSINE_STEREO, TG_FLAG_PHOSPHINE_STEREO,
    };
    let mut o = InpInchiOptions::default();
    let mut mode = 0u64;
    let mut taut = 0u64;
    let mut has_stereo = 0;
    let (mut nksb, mut nksc, mut nusb, mut nusc, mut npiii, mut nasiii) = (0, 0, 0, 0, 0, 0);
    for i_inchi in 0..2usize {
        for i_mob in 0..2usize {
            let mut b_iso = 1i32;
            while mode == 0 && b_iso >= 0 {
                match d.s[i_inchi][i_mob][b_iso as usize] {
                    1 => mode |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO,
                    2 => mode |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RELATIVE_STEREO,
                    3 => mode |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RACEMIC_STEREO,
                    _ => {}
                }
                b_iso -= 1;
            }
            let max_components = if d.p_inp_inchi[i_inchi][i_mob].is_empty() {
                0
            } else {
                d.n_num_components[i_inchi][i_mob]
            };
            for k in 0..max_components.max(0) as usize {
                let Some(c) = d.p_inp_inchi[i_inchi][i_mob].get(k) else {
                    break;
                };
                let ret = count_stereo_types(
                    c,
                    &mut nksb,
                    &mut nksc,
                    &mut nusb,
                    &mut nusc,
                    &mut npiii,
                    &mut nasiii,
                );
                if ret < 0 {
                    return Err(ret);
                }
                has_stereo += i32::from(ret == 2);
                if ret > 0 {
                    o.has_reconnected |= i_inchi == INCHI_REC;
                    o.has_fixed_h |= i_mob == TAUT_NON;
                }
                o.has_metal |= b_inchi_has_reconnected_metal(c);
            }
        }
    }
    if (mode & REQ_MODE_RELATIVE_STEREO) != 0 && (mode & REQ_MODE_RACEMIC_STEREO) != 0 {
        return Err(RI_ERR_SYNTAX);
    }
    if has_stereo != 0 && mode == 0 {
        mode = REQ_MODE_STEREO | REQ_MODE_ISO_STEREO;
    }
    if !(nksb == 0 && nusb != 0) {
        mode |= REQ_MODE_SB_IGN_ALL_UU;
    }
    if !(nksc == 0 && nusc != 0) {
        mode |= REQ_MODE_SC_IGN_ALL_UU;
    }
    if npiii != 0 {
        taut |= TG_FLAG_PHOSPHINE_STEREO;
    }
    if nasiii != 0 {
        taut |= TG_FLAG_ARSINE_STEREO;
    }
    o.mode_flags_stereo = mode;
    o.taut_flags_stereo = taut;
    Ok(o)
}

/// `RemoveFixHInChIIdentical2MobH`: eliminates the Fixed-H InChI that are
/// exactly the same as the corresponding Mobile-H ones
pub fn remove_fix_h_inchi_identical2_mob_h(d: &mut InpInChI) {
    for i_rec in 0..2usize {
        let n = d.n_num_components[i_rec][TAUT_YES].min(d.n_num_components[i_rec][TAUT_NON]);
        for k in 0..n.max(0) as usize {
            let (Some(m), Some(f)) = (
                d.p_inp_inchi[i_rec][TAUT_YES].get(k),
                d.p_inp_inchi[i_rec][TAUT_NON].get(k),
            ) else {
                break;
            };
            if super::ichimake::compare_reversed_inchi(Some(m), Some(f), None, None) == 0 {
                d.p_inp_inchi[i_rec][TAUT_NON][k] = InChI::default();
            }
        }
    }
}

/// `MarkDisconectedIdenticalToReconnected`: returns the number of the
/// marked components
pub fn mark_disconnected_identical_to_reconnected(d: &mut InpInChI) -> i32 {
    use super::ichimake::compare_reversed_inchi as cmp;
    let mut num_marked = 0;
    let nc = d.n_num_components;
    let k1max = nc[INCHI_BAS][TAUT_YES].max(nc[INCHI_BAS][TAUT_NON]);
    for k1 in 0..k1max.max(0) {
        let k2max = nc[INCHI_REC][TAUT_YES].max(nc[INCHI_REC][TAUT_NON]);
        for k2 in 0..k2max.max(0) {
            let (u1, u2) = (k1 as usize, k2 as usize);
            let pi = &d.p_inp_inchi;
            let eq_m = k1 < nc[INCHI_BAS][TAUT_YES]
                && k2 < nc[INCHI_REC][TAUT_YES]
                && pi[INCHI_REC][TAUT_YES][u2].n_link == 0
                && pi[INCHI_BAS][TAUT_YES][u1].b_deleted == 0
                && pi[INCHI_BAS][TAUT_YES][u1].n_number_of_atoms != 0
                && pi[INCHI_BAS][TAUT_YES][u1].n_number_of_atoms
                    == pi[INCHI_REC][TAUT_YES][u2].n_number_of_atoms
                && pi[INCHI_REC][TAUT_YES][u2].b_deleted == 0
                && cmp(
                    Some(&pi[INCHI_REC][TAUT_YES][u2]),
                    Some(&pi[INCHI_BAS][TAUT_YES][u1]),
                    None,
                    None,
                ) == 0;
            let is_f1 = k1 < nc[INCHI_BAS][TAUT_NON]
                && pi[INCHI_BAS][TAUT_NON][u1].b_deleted == 0
                && pi[INCHI_BAS][TAUT_NON][u1].n_number_of_atoms > 0;
            let is_f2 = k2 < nc[INCHI_REC][TAUT_NON]
                && pi[INCHI_REC][TAUT_NON][u2].b_deleted == 0
                && pi[INCHI_REC][TAUT_NON][u2].n_number_of_atoms > 0;
            let eq_f = is_f1
                && is_f2
                && pi[INCHI_REC][TAUT_NON][u2].n_link == 0
                && pi[INCHI_BAS][TAUT_NON][u1].n_number_of_atoms
                    == pi[INCHI_REC][TAUT_NON][u2].n_number_of_atoms
                && cmp(
                    Some(&pi[INCHI_REC][TAUT_NON][u2]),
                    Some(&pi[INCHI_BAS][TAUT_NON][u1]),
                    None,
                    None,
                ) == 0;
            if eq_m && ((!is_f1 && !is_f2) || eq_f) {
                d.p_inp_inchi[INCHI_BAS][TAUT_YES][u1].n_link = -(k2 + 1);
                d.p_inp_inchi[INCHI_REC][TAUT_YES][u2].n_link = k1 + 1;
                if eq_f {
                    d.p_inp_inchi[INCHI_BAS][TAUT_NON][u1].n_link = -(k2 + 1);
                    d.p_inp_inchi[INCHI_REC][TAUT_NON][u2].n_link = k1 + 1;
                }
                num_marked += 1;
                break;
            }
        }
    }
    num_marked
}

/// `SetProtonsAndXchgIsoH`: decides how to treat the removed protons and
/// the exchangeable isotopic H; `Err`/`Ok` carry `ret2`. Returns the
/// numbers of the components and the modes.
#[allow(clippy::too_many_arguments)]
pub fn set_protons_and_xchg_iso_h(
    b_inchi2_structure: bool,
    b_req_split_output_inchi: bool,
    b_req_protons_for_each_component: bool,
    b_req_non_taut: bool,
    b_req_stereo: bool,
    d: &mut InpInChI,
) -> (i32, [i32; 2], [i32; 2]) {
    let mut num_components = [0i32; 2];
    let mut mode = [MODE_PIXH_UNDEFINED; 2];
    for i_inchi in 0..2usize {
        mode[i_inchi] = MODE_PIXH_UNDEFINED;
        let np = &d.n_num_protons[i_inchi][TAUT_YES];
        let mut avail_total = np.n_num_removed_protons != 0;
        for k1 in 0..NUM_H_ISOTOPES {
            avail_total |= np.n_num_removed_isotopic_h[k1] != 0;
        }
        let mut avail_each = np.p_num_protons.is_some();
        if b_inchi2_structure {
            mode[i_inchi] = if avail_each {
                MODE_PIXH_ADD_TO_EACH
            } else {
                MODE_PIXH_KEEP_TOTALS
            };
        } else if !b_req_split_output_inchi {
            mode[i_inchi] = if avail_each {
                MODE_PIXH_ADD_TO_EACH
            } else {
                MODE_PIXH_ADD_TO_FIRST
            };
        } else if !avail_each {
            mode[i_inchi] = if avail_total {
                MODE_PIXH_ADD_A_PIXH_COMPONENT
            } else {
                MODE_PIXH_ADD_TO_FIRST
            };
        } else if b_req_protons_for_each_component {
            mode[i_inchi] = MODE_PIXH_ADD_TO_EACH;
        } else {
            mode[i_inchi] = if b_req_non_taut {
                MODE_PIXH_ADD_TO_EACH
            } else {
                MODE_PIXH_ADD_A_PIXH_COMPONENT
            };
        }
        // remove unneeded data: protons for each component
        if avail_each && mode[i_inchi] != MODE_PIXH_ADD_TO_EACH {
            d.n_num_protons[i_inchi][TAUT_YES].p_num_protons = None;
            avail_each = false;
        }
        // remove unneeded data: total protons of all components
        if avail_total && mode[i_inchi] == MODE_PIXH_ADD_TO_EACH {
            let np = &mut d.n_num_protons[i_inchi][TAUT_YES];
            np.n_num_removed_protons = 0;
            np.n_num_removed_isotopic_h = [0; 3];
        }
        // remove unneeded data: Fixed-H InChI
        if !b_req_non_taut && d.n_num_components[i_inchi][TAUT_NON] != 0 {
            d.p_inp_inchi[i_inchi][TAUT_NON] = Vec::new();
            d.n_num_components[i_inchi][TAUT_NON] = 0;
        }
        // add one more component containing only /p and /i/h
        if (mode[i_inchi] == MODE_PIXH_ADD_A_PIXH_COMPONENT
            && d.n_num_components[i_inchi][TAUT_YES] != 0)
            || (b_inchi2_structure
                && !avail_each
                && d.n_num_components[i_inchi][TAUT_NON] == 0
                && d.n_num_components[i_inchi][TAUT_YES] != 0)
        {
            let j = TAUT_YES;
            let prev_len = d.n_num_components[i_inchi][j];
            let mut n_len = 0;
            for k in 0..prev_len.max(0) as usize {
                n_len += i32::from(
                    d.p_inp_inchi[i_inchi][j]
                        .get(k)
                        .is_none_or(|c| c.b_deleted == 0),
                );
            }
            let v = &mut d.p_inp_inchi[i_inchi][j];
            if n_len == prev_len {
                v.truncate(n_len.max(0) as usize);
                v.push(InChI::default());
            }
            d.n_num_components[i_inchi][j] = n_len + 1;
            let v = &mut d.p_inp_inchi[i_inchi][j];
            for k in n_len.max(0) as usize..prev_len.max(0) as usize {
                if let Some(c) = v.get_mut(k) {
                    *c = InChI::default();
                }
            }
            let r = n_fill_out_proton_mobile_h(&mut v[n_len as usize]);
            if r < 0 {
                return (r, num_components, mode);
            }
        }
        // remove unneeded Stereo
        if !b_req_stereo {
            for j in 0..2usize {
                for k in 0..d.n_num_components[i_inchi][j].max(0) as usize {
                    if let Some(c) = d.p_inp_inchi[i_inchi][j].get_mut(k) {
                        c.stereo = None;
                        c.stereo_isotopic = None;
                    }
                }
            }
        }
    }
    num_components[INCHI_BAS] =
        d.n_num_components[INCHI_BAS][TAUT_YES].max(d.n_num_components[INCHI_BAS][TAUT_NON]);
    num_components[INCHI_REC] =
        d.n_num_components[INCHI_REC][TAUT_YES].max(d.n_num_components[INCHI_REC][TAUT_NON]);
    (0, num_components, mode)
}
