//! Port of the SD-file reading helpers of `mol_fmt4.c`
//! (`SDFileSkipExtraData`, `SDFileIdentifyLabel`, `SDFileExtractCASNo`). The
//! writer half of `mol_fmt4.c` is executable-only code and is not ported.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust. Licensed under the MIT license
//! (see LICENSES/InChI-MIT.txt).

use super::ichi_io::{InchiIostream, inchi_fgets_lf};
use super::ichierr::add_error_message;
use super::molfmt::{MOL_FMT_INPLINELEN, MOL_FMT_MAXLINELEN, SD_FMT_END_OF_DATA, treat_err};
use super::util::{
    dotify_non_printable_chars, inchi_memicmp, mystrncpy, normalize_string, remove_trailing_spaces,
};

pub const MAX_SDF_VALUE: usize = 255;

const SDF_START: i32 = 0;
const SDF_DATA_HEADER: i32 = 1;
const SDF_DATA_HEADER_NAME: i32 = 2;
const SDF_DATA_HEADER_COMMENT: i32 = 3;
const SDF_DATA_HEADER_CAS: i32 = 4;
const SDF_DATA_HEADER_USER: i32 = 5;
const SDF_DATA_LINE: i32 = 6;
const SD_FMT_END_OF_DATA_ITEM: i32 = 7;
const SDF_EMPTY_LINE: i32 = 8;
const SD_FMT_END_OF_DATA_BLOCK: i32 = 9;

fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// `inchi_memicmp` over `n` bytes, reading zeros past the end of the slices
fn memicmp(a: &[u8], b: &[u8], n: usize) -> i32 {
    let pad = |s: &[u8]| {
        let mut v = s.to_vec();
        v.resize(n.max(v.len()), 0);
        v
    };
    inchi_memicmp(&pad(a), &pad(b), n)
}

/// `SDFileExtractCASNo`: keeps the digits (and a leading '-') of `line` in
/// place and returns them as a number
pub fn sdfile_extract_cas_no(line: &mut Vec<u8>) -> u64 {
    let mut i = usize::from(at(line, 0) == b'-');
    let mut j = i;
    while at(line, i) != 0 {
        let c = line[i];
        if c.is_ascii_digit() {
            line[j] = c;
            j += 1;
        } else if c != b'-' {
            break;
        }
        i += 1;
    }
    line.truncate(j);
    // strtoul( line, NULL, 10 )
    let mut k = 0;
    let neg = at(line, 0) == b'-';
    if neg || at(line, 0) == b'+' {
        k = 1;
    }
    let mut v: u128 = 0;
    while at(line, k).is_ascii_digit() {
        v = (v * 10 + u128::from(line[k] - b'0')).min(u128::from(u64::MAX));
        k += 1;
    }
    let v = v as u64;
    if neg { v.wrapping_neg() } else { v }
}

/// `SDFileIdentifyLabel`
pub fn sdfile_identify_label(inp_line: &[u8], sdf_label: Option<&[u8]>) -> i32 {
    let s = &inp_line[..inp_line
        .iter()
        .position(|&c| c == 0)
        .unwrap_or(inp_line.len())];
    if let Some(p) = s.iter().position(|&c| c == b'<')
        && let Some(qq) = s[p..].iter().position(|&c| c == b'>')
    {
        let len = qq as i64 - 1;
        if len > 0 && len < MOL_FMT_MAXLINELEN as i64 {
            let line = &s[p + 1..p + 1 + len as usize];
            let len = len as usize;
            let sp = |c: u8| matches!(c, b' ' | 0x09..=0x0d);
            let mut tmp1 = 0usize;
            for (i, &c) in line.iter().enumerate() {
                if sp(c) {
                    tmp1 = i;
                    break;
                }
            }
            let mut tmp2 = 0usize;
            let mut j = len as i64 - 1;
            while j >= tmp1 as i64 {
                if sp(line[j as usize]) {
                    tmp2 = j as usize;
                    break;
                }
                j -= 1;
            }
            let plen = (tmp2 as i64 - tmp1 as i64 + 1).max(0) as usize;
            let pp = &line[tmp1.min(line.len())..];
            if let Some(l) = sdf_label.filter(|l| !l.is_empty())
                && plen == l.len()
                && memicmp(pp, l, plen) == 0
            {
                return SDF_DATA_HEADER_USER;
            }
            if plen == b"NAME".len() && memicmp(pp, b"NAME", plen) == 0 {
                return SDF_DATA_HEADER_NAME;
            }
            if plen == b"COMMENT".len() && memicmp(pp, b"COMMENT", plen) == 0 {
                return SDF_DATA_HEADER_COMMENT;
            }
            if memicmp(pp, b"CAS", 3) == 0 {
                return SDF_DATA_HEADER_CAS;
            }
        }
    }
    SDF_DATA_HEADER
}

fn set_c_string(target: &mut Vec<u8>, line: &[u8], size: usize) {
    let mut buf = vec![0u8; size.max(1)];
    mystrncpy(&mut buf, line, size);
    let n = buf.iter().position(|&c| c == 0).unwrap_or(buf.len());
    buf.truncate(n);
    *target = buf;
}

/// the buffers `SDFileSkipExtraData` may fill (C strings: empty = `\0`)
#[derive(Default, Debug)]
pub struct SdfOut<'a> {
    pub cas_num: Option<&'a mut u64>,
    pub comment: Option<(&'a mut Vec<u8>, usize)>,
    pub name: Option<(&'a mut Vec<u8>, usize)>,
    pub sdf_label: Option<&'a [u8]>,
    pub sdf_value: Option<&'a mut Vec<u8>>,
}

/// `SDFileSkipExtraData`
pub fn sdfile_skip_extra_data(
    inp: &mut InchiIostream,
    out: SdfOut,
    prev_err: i32,
    p_str_err: &mut String,
    no_warnings: bool,
) -> i32 {
    let SdfOut {
        mut cas_num,
        mut comment,
        mut name,
        sdf_label,
        mut sdf_value,
    } = out;
    let mut n_blank_lines = 0;
    let mut n_lines = 0;
    let mut current_state = SDF_START;
    let mut err = 0;
    let mut wait_for_cas = false;
    let mut cas_num_is_user = false;
    let mut wait_for_name = name.as_ref().is_some_and(|(n, l)| *l > 0 && n.is_empty());
    let mut wait_for_comment = comment
        .as_ref()
        .is_some_and(|(n, l)| *l > 0 && n.is_empty());
    let label = sdf_label.filter(|l| !l.is_empty());
    let mut wait_for_user = label.is_some() && sdf_value.is_some();
    if let Some(c) = cas_num.as_deref_mut() {
        wait_for_cas = true;
        *c = 0;
        cas_num_is_user = wait_for_user && memicmp(label.unwrap(), b"CAS", 3) == 0;
    }
    let mut line: Vec<u8> = Vec::new();
    let mut got_p = false;
    let get_line = |inp: &mut InchiIostream| {
        inchi_fgets_lf(MOL_FMT_INPLINELEN, inp).map(|mut l| {
            l.truncate(l.iter().position(|&c| c == 0).unwrap_or(l.len()));
            l
        })
    };
    while err == 0 && current_state != SD_FMT_END_OF_DATA_BLOCK {
        let Some(l) = get_line(inp) else {
            got_p = false;
            break;
        };
        got_p = true;
        line = l;
        if n_lines == 0 && line.starts_with(b"M  END") {
            continue;
        }
        n_lines += 1;
        remove_trailing_spaces(&mut line);
        if at(&line, MOL_FMT_MAXLINELEN) != 0 {
            if !matches!(
                current_state,
                SDF_DATA_HEADER
                    | SDF_DATA_LINE
                    | SDF_DATA_HEADER_NAME
                    | SDF_DATA_HEADER_USER
                    | SDF_DATA_HEADER_COMMENT
            ) {
                line.truncate(MOL_FMT_MAXLINELEN);
                if prev_err == 0 {
                    treat_err(&mut err, 0, p_str_err, "Too long SData line truncated");
                }
            } else {
                line.truncate(MOL_FMT_MAXLINELEN);
            }
        }
        if line.is_empty() {
            n_blank_lines += 1;
        }
        let mut got_data_line = false;
        match current_state {
            SDF_START | SD_FMT_END_OF_DATA_ITEM | SDF_EMPTY_LINE => {
                if line == SD_FMT_END_OF_DATA {
                    current_state = SD_FMT_END_OF_DATA_BLOCK;
                } else if at(&line, 0) == b'>' {
                    current_state =
                        if wait_for_name || wait_for_comment || wait_for_cas || wait_for_user {
                            sdfile_identify_label(&line, sdf_label)
                        } else {
                            SDF_DATA_HEADER
                        };
                } else if line.is_empty() {
                    current_state = SDF_EMPTY_LINE;
                } else if prev_err == 0 {
                    treat_err(&mut err, 3, p_str_err, "Unexpected SData header line:");
                    dotify_non_printable_chars(&mut line);
                    add_error_message(p_str_err, &String::from_utf8_lossy(&line));
                } else {
                    err = 3;
                }
            }
            SDF_DATA_HEADER_NAME => {
                if wait_for_name && normalize_string(&mut line) > 0 {
                    wait_for_name = false;
                    let (n, l) = name.as_mut().unwrap();
                    set_c_string(n, &line, *l);
                }
                got_data_line = true;
            }
            SDF_DATA_HEADER_COMMENT => {
                if wait_for_comment && normalize_string(&mut line) > 0 {
                    wait_for_comment = false;
                    let (n, l) = comment.as_mut().unwrap();
                    set_c_string(n, &line, *l);
                }
                got_data_line = true;
            }
            SDF_DATA_HEADER_USER => {
                if wait_for_user && normalize_string(&mut line) > 0 {
                    wait_for_user = false;
                    set_c_string(sdf_value.as_mut().unwrap(), &line, MAX_SDF_VALUE + 1);
                    if cas_num_is_user && wait_for_cas {
                        let c = sdfile_extract_cas_no(&mut line);
                        **cas_num.as_mut().unwrap() = c;
                        wait_for_cas = c == 0;
                    }
                }
                got_data_line = true;
            }
            SDF_DATA_HEADER_CAS => {
                if wait_for_cas && normalize_string(&mut line) > 0 {
                    let c = sdfile_extract_cas_no(&mut line);
                    **cas_num.as_mut().unwrap() = c;
                    wait_for_cas = c == 0;
                }
                got_data_line = true;
            }
            _ => got_data_line = true,
        }
        if got_data_line {
            current_state = if at(&line, 0) != 0 {
                SDF_DATA_LINE
            } else {
                SD_FMT_END_OF_DATA_ITEM
            };
        }
    }
    if err != 0 && n_blank_lines == n_lines && at(&line, 0) == 0 {
        err = 5;
    }
    if err != 0 && err != 5 && current_state != SD_FMT_END_OF_DATA_BLOCK && got_p {
        let mut reached = false;
        while let Some(l) = get_line(inp) {
            if l.starts_with(SD_FMT_END_OF_DATA) {
                reached = true;
                break;
            }
        }
        if reached {
            err = 9;
            if !no_warnings {
                add_error_message(p_str_err, "Bypassing to next structure");
            }
        }
    }
    err
}
