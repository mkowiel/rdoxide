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
