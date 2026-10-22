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

use super::ichi::{InChI, InChIIsotopicAtom, InChIStereo};
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
        let f = from.unwrap();
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
