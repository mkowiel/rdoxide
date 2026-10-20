//! Error codes and messages of the InChI library (`INCHI_BASE/src/ichierr.{h,c}`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

pub const CT_ERR_FIRST: i32 = -30000;
pub const CT_OVERFLOW: i32 = CT_ERR_FIRST;
pub const CT_LEN_MISMATCH: i32 = CT_ERR_FIRST - 1;
pub const CT_OUT_OF_RAM: i32 = CT_ERR_FIRST - 2;
pub const CT_RANKING_ERR: i32 = CT_ERR_FIRST - 3;
pub const CT_ISOCOUNT_ERR: i32 = CT_ERR_FIRST - 4;
pub const CT_TAUCOUNT_ERR: i32 = CT_ERR_FIRST - 5;
pub const CT_ISOTAUCOUNT_ERR: i32 = CT_ERR_FIRST - 6;
pub const CT_MAPCOUNT_ERR: i32 = CT_ERR_FIRST - 7;
pub const CT_TIMEOUT_ERR: i32 = CT_ERR_FIRST - 8;
pub const CT_ISO_H_ERR: i32 = CT_ERR_FIRST - 9;
pub const CT_STEREOCOUNT_ERR: i32 = CT_ERR_FIRST - 10;
pub const CT_ATOMCOUNT_ERR: i32 = CT_ERR_FIRST - 11;
pub const CT_STEREOBOND_ERROR: i32 = CT_ERR_FIRST - 12;
pub const CT_USER_QUIT_ERR: i32 = CT_ERR_FIRST - 13;
pub const CT_REMOVE_STEREO_ERR: i32 = CT_ERR_FIRST - 14;
pub const CT_CALC_STEREO_ERR: i32 = CT_ERR_FIRST - 15;
pub const CT_CANON_ERR: i32 = CT_ERR_FIRST - 16;
pub const CT_STEREO_CANON_ERR: i32 = CT_ERR_FIRST - 17;
pub const CT_WRONG_FORMULA: i32 = CT_ERR_FIRST - 18;
pub const CT_UNKNOWN_ERR: i32 = CT_ERR_FIRST - 19;
pub const CT_ERR_MIN: i32 = CT_UNKNOWN_ERR;
pub const CT_ERR_MAX: i32 = CT_ERR_FIRST;

pub const BNS_ERR: i32 = -9999;
pub const BNS_WRONG_PARMS: i32 = BNS_ERR;
pub const BNS_OUT_OF_RAM: i32 = BNS_ERR + 1;
pub const BNS_PROGRAM_ERR: i32 = BNS_ERR + 2;
pub const BNS_ALTPATH_OVFL: i32 = BNS_ERR + 3;
pub const BNS_BOND_ERR: i32 = BNS_ERR + 4;
pub const BNS_VERT_NUM_ERR: i32 = BNS_ERR + 5;
pub const BNS_VERT_EDGE_OVFL: i32 = BNS_ERR + 6;
pub const BNS_SET_ALTP_ERR: i32 = BNS_ERR + 7;
pub const BNS_CPOINT_ERR: i32 = BNS_ERR + 8;
pub const BNS_CANT_SET_BOND: i32 = BNS_ERR + 9;
pub const BNS_CAP_FLOW_ERR: i32 = BNS_ERR + 10;
pub const BNS_RADICAL_ERR: i32 = BNS_ERR + 11;
pub const BNS_REINIT_ERR: i32 = BNS_ERR + 12;
pub const BNS_ALTBOND_ERR: i32 = BNS_ERR + 13;
/// v. 1.05
pub const BNS_TIMEOUT: i32 = BNS_ERR + 14;
pub const BNS_MAX_ERR_VALUE: i32 = BNS_ERR + 19;

/// `IS_BNS_ERROR`
pub fn is_bns_error(x: i32) -> bool {
    (BNS_ERR..=BNS_MAX_ERR_VALUE).contains(&x)
}

pub const INCHI_INP_ERROR_ERR: i32 = 40;
pub const INCHI_INP_ERROR_RET: i32 = -1;
pub const INCHI_INP_FATAL_ERR: i32 = 1;
pub const INCHI_INP_FATAL_RET: i32 = 0;
pub const INCHI_INP_EOF_ERR: i32 = 11;
pub const INCHI_INP_EOF_RET: i32 = 0;

pub const LOG_MASK_WARN: i32 = 1;
pub const LOG_MASK_ERR: i32 = 2;
pub const LOG_MASK_FATAL: i32 = 4;
pub const LOG_MASK_ALL: i32 = LOG_MASK_WARN | LOG_MASK_ERR | LOG_MASK_FATAL;
pub const LOG_MASK_NO_WARN: i32 = LOG_MASK_ERR | LOG_MASK_FATAL;

pub const USER_ACTION_QUIT: i32 = 1;

/// Size of the buffer of the accumulated error messages.
pub const STR_ERR_LEN: usize = 256;

/// `ErrMsg`
pub fn err_msg(error_code: i32) -> String {
    let p = match error_code {
        0 => "",
        CT_OVERFLOW => "ARRAY OVERFLOW",
        CT_LEN_MISMATCH => "LENGTH_MISMATCH",
        CT_OUT_OF_RAM => "Out of RAM",
        CT_RANKING_ERR => "RANKING_ERR",
        CT_ISOCOUNT_ERR => "ISOCOUNT_ERR",
        CT_TAUCOUNT_ERR => "TAUCOUNT_ERR",
        CT_ISOTAUCOUNT_ERR => "ISOTAUCOUNT_ERR",
        CT_MAPCOUNT_ERR => "MAPCOUNT_ERR",
        CT_TIMEOUT_ERR => "Time limit exceeded",
        CT_ISO_H_ERR => "ISO_H_ERR",
        CT_STEREOCOUNT_ERR => "STEREOCOUNT_ERR",
        CT_ATOMCOUNT_ERR => "ATOMCOUNT_ERR",
        CT_STEREOBOND_ERROR => "STEREOBOND_ERR",
        CT_USER_QUIT_ERR => "User requested termination",
        CT_REMOVE_STEREO_ERR => "REMOVE_STEREO_ERR",
        CT_CALC_STEREO_ERR => "CALC_STEREO_ERR",
        CT_STEREO_CANON_ERR => "STEREO_CANON_ERR",
        CT_CANON_ERR => "CANON_ERR",
        CT_WRONG_FORMULA => "Wrong or missing chemical formula",
        CT_UNKNOWN_ERR => "UNKNOWN_ERR",
        BNS_RADICAL_ERR => "Cannot process free radical center",
        BNS_ALTBOND_ERR => "Cannot process aromatic bonds",
        BNS_TIMEOUT => "Structure normalization timeout",
        c if c > CT_UNKNOWN_ERR => return format!("No description({c})"),
        c => return format!("UNKNOWN_ERR({})", CT_UNKNOWN_ERR - c),
    };
    p.to_string()
}

/// `already_have_this_message`: whether `new_message` is already one of the
/// `"; "`- or `": "`-separated messages of `prev_messages`.
fn already_have_this_message(prev_messages: &str, new_message: &str) -> bool {
    let (prev, new) = (prev_messages.as_bytes(), new_message.as_bytes());
    let Some(p) = prev_messages.find(new_message) else {
        return false;
    };
    // the C code reads `*(p - 2)` unguarded only after `*(p - 1) == ' '`
    let at_start =
        p == 0 || (prev[p - 1] == b' ' && p >= 2 && (prev[p - 2] == b';' || prev[p - 2] == b':'));
    if !at_start {
        return false;
    }
    let end = p + new.len();
    let at = |i: usize| prev.get(i).copied().unwrap_or(0);
    end == prev.len()
        || (at(end) == b';' && at(end + 1) == b' ')
        || (new.last() == Some(&b':') && at(end) == b' ')
}

/// `AddErrorMessage`: appends `new_message` to the messages in `all_messages`
/// (a buffer of [`STR_ERR_LEN`] bytes) unless it is already there; returns
/// whether the message is in the buffer afterwards.
pub fn add_error_message(all_messages: &mut String, new_message: &str) -> bool {
    if new_message.is_empty() {
        return false;
    }
    if already_have_this_message(all_messages, new_message) {
        return true;
    }
    let len_all = all_messages.len();
    let len = new_message.len();
    if len_all + len + 2 * usize::from(len_all > 0) < STR_ERR_LEN {
        // enough room
        if len_all > 0 {
            if !all_messages.ends_with(':') {
                all_messages.push(';');
            }
            all_messages.push(' ');
        }
        all_messages.push_str(new_message);
        return true;
    }
    // not enough room: add the marker if it is not there yet
    if all_messages.contains("...") {
        return false;
    }
    if len_all + 3 < STR_ERR_LEN {
        all_messages.push_str("...");
    }
    false
}
