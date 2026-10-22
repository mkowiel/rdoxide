//! The InChI -> structure entry points of `INCHI_API/libinchi/src/inchi_dll.c`
//! (`CheckINCHI`, `GetStructFromINCHIEx`, `GetStructFromINCHI`,
//! `GetStructFromStdINCHI`) with the single-structure path of
//! `ReadWriteInChI` (`ichiread.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! Polymer InChI (`/z` layer) are not supported: the reader reports them as
//! a syntax error and `DetectAndExposePolymerInternals` is not ported.

use super::ichi_io::{
    INCHI_IOS_TYPE_STRING, InchiIostream, IosFile, inchi_ios_eprint, inchi_ios_init,
};
use super::ichiparm::{
    INCHI_MAX_NUM_ARG, INCHI_OPTION_PREFX, INPUT_INCHI, InputParms, READ_INCHI_TO_STRUCTURE,
    REQ_MODE_BASIC, REQ_MODE_ISO, REQ_MODE_ISO_STEREO, REQ_MODE_STEREO, REQ_MODE_TAUT,
    TG_FLAG_DISCONNECT_COORD, TG_FLAG_RECONNECT_COORD, parse_options_string, print_input_parms,
    read_command_line_parms,
};
use super::ichiread::{
    InpInChI, RI_ERR_EOF, RI_ERR_EOL, RI_ERR_PROGR, RI_ERR_SYNTAX, SegmLine, get_inchi_err_name,
    get_inchi_state_read_err, inchi_line2data, set_protons_and_xchg_iso_h,
};
use super::ichirvr7::convert_inchi2_struct;
use super::inchi_dll::{InchiAtom, inp_atom_0d_to_inchi_atom};
use super::inchi_dll_b::{InchiOutput, copy_corrected_log_tail};
use super::inpdef::InpAtom;
use super::readinch::InchiStereo0D;

/// `INCHI_VALID_*` / `INCHI_INVALID_*` return codes of `CheckINCHI`
pub const INCHI_VALID_STANDARD: i32 = 0;
pub const INCHI_VALID_NON_STANDARD: i32 = 1;
pub const INCHI_VALID_BETA: i32 = 2;
pub const INCHI_INVALID_PREFIX: i32 = 3;
pub const INCHI_INVALID_VERSION: i32 = 4;
pub const INCHI_INVALID_LAYOUT: i32 = 5;
pub const INCHI_FAIL_I2I: i32 = 6;

/// `CheckINCHI` for `strict == false` (the strict mode round-trips the
/// string through `GetINCHIfromINCHI`, which is not ported, and reports
/// `INCHI_FAIL_I2I`).
pub fn check_inchi(sz_inchi: &str, strict: bool) -> i32 {
    const PREFIX: &[u8] = b"InChI=";
    let s = sz_inchi.as_bytes();
    let mut ret = INCHI_VALID_NON_STANDARD;
    if s.len() < PREFIX.len() + 3 || &s[..PREFIX.len()] != PREFIX {
        return INCHI_INVALID_PREFIX;
    }
    if s[PREFIX.len()] != b'1' {
        return INCHI_INVALID_VERSION;
    }
    let mut pos_slash1 = PREFIX.len() + 1;
    if s[pos_slash1] == b'S' {
        ret = INCHI_VALID_STANDARD;
        pos_slash1 += 1;
    } else if s[pos_slash1] == b'B' {
        ret = INCHI_VALID_BETA;
        pos_slash1 += 1;
    }
    if s.get(pos_slash1) != Some(&b'/') {
        return INCHI_INVALID_LAYOUT;
    }
    // do not check the trailing white spaces
    let mut slen = s.len();
    while slen > 0 && s[slen - 1].is_ascii_whitespace() {
        slen -= 1;
    }
    // possible SaveOpt letters
    let mut slen0 = slen;
    if slen >= 3
        && s[slen - 3] == b'\\'
        && s[slen - 2].is_ascii_uppercase()
        && s[slen - 1].is_ascii_uppercase()
    {
        slen0 = slen - 3;
    }
    let mut prev_is_slash = true;
    for (i, &pp) in s.iter().enumerate().take(slen0).skip(pos_slash1 + 1) {
        if prev_is_slash {
            // '0' is never allowed after a slash; only lowercase letters
            // outside the formula layer
            if pp == b'0' {
                return INCHI_INVALID_LAYOUT;
            }
            if i > pos_slash1 + 1 && !pp.is_ascii_lowercase() {
                return INCHI_INVALID_LAYOUT;
            }
        }
        prev_is_slash = pp == b'/';
        if pp.is_ascii_alphanumeric() {
            continue;
        }
        match pp {
            b'(' | b')' | b'*' | b'+' | b',' | b'-' | b'.' | b'/' | b';' | b'?' => {}
            _ => return INCHI_INVALID_LAYOUT,
        }
    }
    if strict {
        return INCHI_FAIL_I2I;
    }
    ret
}

/// `inchi_OutputStructEx` (without the polymer data)
#[derive(Debug, Clone, Default)]
pub struct InchiOutputStructEx {
    pub atom: Vec<InchiAtom>,
    pub stereo_0d: Vec<InchiStereo0D>,
    pub sz_message: Option<String>,
    pub sz_log: Option<String>,
    pub warning_flags: [[u64; 2]; 2],
}

/// `inchi_Ret_*`
pub const INCHI_RET_EOF: i32 = -1;
pub const INCHI_RET_OKAY: i32 = 0;
pub const INCHI_RET_WARNING: i32 = 1;
pub const INCHI_RET_ERROR: i32 = 2;
pub const INCHI_RET_FATAL: i32 = 3;

/// the single-structure path of `ReadWriteInChI` for the conversion InChI ->
/// structure: returns the result and, on success, the restored atoms; the
/// message and the warning flags go to `out`
fn read_write_inchi(
    ip_inp: &InputParms,
    log: &mut InchiIostream,
    sz_inchi: &str,
    out: &mut InchiOutputStructEx,
) -> (i32, Vec<InpAtom>) {
    const MAX_MSG_LEN: usize = 512;
    if ip_inp.b_read_inchi_options & READ_INCHI_TO_STRUCTURE == 0 {
        // InChI -> InChI is not an API mode
        inchi_ios_eprint(
            Some(log),
            "\nWrong command line options: expected Inch2Struct or Inchi2Inchi\n",
        );
        return (0, Vec::new());
    }
    let mut ip = ip_inp.clone();
    let mut one_input = InpInChI::default();
    let mut line = SegmLine::default();
    let mut str_hdr: Option<Vec<u8>> = None;
    let (mut state, mut err) = (0, 0);
    let (mut std_fmt, mut has_save_opt, mut save_opt_bits) = (0, 0, 0u8);
    ip.b_taut_flags |= TG_FLAG_DISCONNECT_COORD | TG_FLAG_RECONNECT_COORD;
    ip.n_mode |=
        REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_ISO;
    // bRELEASE_VERSION == 1: the Fixed-H layer of the input is always used
    let b_req_non_taut = true;

    let mut inp = InchiIostream::default();
    inchi_ios_init(&mut inp, INCHI_IOS_TYPE_STRING, IosFile::None);
    let bytes = sz_inchi.as_bytes().to_vec();
    inp.s.n_used_length = bytes.len() as i32;
    inp.s.p_str = Some(bytes);
    inp.s.n_ptr = 0;

    let read_status = inchi_line2data(
        &mut inp,
        &mut line,
        &mut str_hdr,
        &mut state,
        &mut err,
        &mut one_input,
        true,
        ip_inp.n_mode,
        &mut std_fmt,
        &mut has_save_opt,
        &mut save_opt_bits,
    );
    let end_of_data = read_status == RI_ERR_EOL || read_status == RI_ERR_EOF;
    let read_ok = end_of_data
        && err == 0
        && one_input.n_num_components[0][1] + one_input.n_num_components[0][0] != 0;
    if !read_ok {
        // TreatErrorsInReadInChIString
        if !(read_status == RI_ERR_EOF && err == 0 && state == 0) {
            let hdr = str_hdr
                .as_ref()
                .map(|h| {
                    String::from_utf8_lossy(h)
                        .trim_end_matches('\0')
                        .to_string()
                })
                .filter(|h| !h.is_empty())
                .unwrap_or_else(|| "Structure: 1".to_string());
            inchi_ios_eprint(
                Some(log),
                &format!(
                    "\n{} {} ({}) in {} ({})\n",
                    hdr,
                    get_inchi_err_name(err),
                    err,
                    get_inchi_state_read_err(state),
                    state
                ),
            );
        }
        let ret = if read_status == RI_ERR_SYNTAX || read_status == RI_ERR_PROGR {
            read_status
        } else {
            0
        };
        return (ret, Vec::new());
    }
    // the InChI has been read: how to treat the removed protons and the
    // exchangeable isotopic H
    let (ret, _num_components, _modes) = set_protons_and_xchg_iso_h(
        true,
        true,
        b_req_non_taut,
        b_req_non_taut,
        true,
        &mut one_input,
    );
    if ret < 0 {
        return (ret, Vec::new());
    }
    let conv = convert_inchi2_struct(&ip, &mut one_input, "", b_req_non_taut, 1, None);
    out.sz_message = Some(conv.message.chars().take(MAX_MSG_LEN - 1).collect());
    if conv.ret < 0 {
        // errors (and the mismatches treated as errors) give no structure
        return (conv.ret, Vec::new());
    }
    out.warning_flags = conv.warning_flags;
    (conv.ret, conv.atoms)
}

/// `GetStructFromINCHIEx`: `sz_options` are the user's options (the
/// conversion option `Inchi2Struct` is added). Returns `inchi_Ret_*` and the
/// output structure.
pub fn get_struct_from_inchi_ex(
    sz_inchi: Option<&str>,
    sz_options: Option<&str>,
) -> (i32, InchiOutputStructEx) {
    let mut out = InchiOutputStructEx::default();
    let mut ip = InputParms::default();
    let mut log = InchiIostream::default();
    inchi_ios_init(&mut log, INCHI_IOS_TYPE_STRING, IosFile::None);
    let mut n_ret: i32 = 0;
    'body: {
        let Some(sz_inchi) = sz_inchi else {
            n_ret = INCHI_RET_ERROR;
            break 'body;
        };
        let opts = format!(
            "{} {}InChI2Struct",
            sz_options.unwrap_or(""),
            INCHI_OPTION_PREFX as char
        );
        let argv_bytes = parse_options_string(opts.as_bytes(), INCHI_MAX_NUM_ARG);
        let argv_str: Vec<String> = argv_bytes
            .iter()
            .map(|a| String::from_utf8_lossy(a).into_owned())
            .collect();
        let argv: Vec<&str> = argv_str.iter().map(String::as_str).collect();
        let r1 = read_command_line_parms(&argv, &mut ip, true, &mut log);
        // INChI DLL specific
        ip.b_no_struct_labels = 1;
        if r1 < 0 {
            break 'body;
        }
        if ip.n_input_type != 0 && ip.n_input_type != INPUT_INCHI {
            inchi_ios_eprint(Some(&mut log), "Input type set to INPUT_INCHI\n");
            ip.n_input_type = INPUT_INCHI;
        }
        let chk = check_inchi(sz_inchi, false);
        if !(chk == INCHI_VALID_STANDARD
            || chk == INCHI_VALID_NON_STANDARD
            || chk == INCHI_VALID_BETA)
        {
            n_ret = INCHI_RET_ERROR;
            break 'body;
        }
        print_input_parms(&mut log, &ip);
        let (r, at) = read_write_inchi(&ip, &mut log, sz_inchi, &mut out);
        n_ret = r;
        if n_ret >= 0 && !at.is_empty() {
            // success
            let (r, atoms, stereo) = inp_atom_0d_to_inchi_atom(&at);
            n_ret = r;
            if n_ret >= 0 {
                out.atom = atoms;
                out.stereo_0d = stereo;
            } else {
                inchi_ios_eprint(Some(&mut log), "Final structure conversion failed\n");
            }
        }
    }
    let mut tmp = InchiOutput::default();
    copy_corrected_log_tail(&mut tmp, &mut log);
    out.sz_log = tmp.sz_log;
    let code = match n_ret {
        -3 | -2 => INCHI_RET_ERROR,
        -1 => INCHI_RET_FATAL,
        _ => {
            if out.atom.is_empty() {
                INCHI_RET_EOF
            } else if out.warning_flags.iter().flatten().any(|&f| f != 0) {
                INCHI_RET_WARNING
            } else {
                INCHI_RET_OKAY
            }
        }
    };
    (code, out)
}

/// `GetStructFromINCHI`
pub fn get_struct_from_inchi(
    sz_inchi: Option<&str>,
    sz_options: Option<&str>,
) -> (i32, InchiOutputStructEx) {
    get_struct_from_inchi_ex(sz_inchi, sz_options)
}

/// `GetStructFromStdINCHI`: only a standard InChI string is accepted
pub fn get_struct_from_std_inchi(
    sz_inchi: Option<&str>,
    sz_options: Option<&str>,
) -> (i32, InchiOutputStructEx) {
    match sz_inchi {
        Some(s) if s.len() >= 6 + 3 && s.as_bytes().get(7) == Some(&b'S') => {
            get_struct_from_inchi(sz_inchi, sz_options)
        }
        _ => (INCHI_RET_ERROR, InchiOutputStructEx::default()),
    }
}
