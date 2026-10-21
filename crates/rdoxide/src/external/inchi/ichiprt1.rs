//! Output layer helpers of `INCHI_BASE/src/ichiprt1.c` (`sCompDelim`,
//! `EquString`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

/// component delimiter
pub const S_COMP_DELIM: &str = ";";

pub const II_STEREO: i32 = 0x0001;
pub const II_STEREO_INV: i32 = 0x0002;
pub const II_NUMB: i32 = 0x0004;
pub const II_EQU: i32 = 0x0008;
pub const IIT_ISO: i32 = 0x0010;
pub const IIT_NONTAUT: i32 = 0x0020;
pub const II_EQ2_NONTAUT: i32 = 0x0040;
pub const II_EQ2_ISO: i32 = 0x0080;
pub const II_EQ2_INV: i32 = 0x0100;
pub const II_EMPTY: i32 = 0x0200;

/// `EquString`: the repetition mark of a layer that equals another one
/// (`FIX_EMPTY_LAYER_BUG` is 0)
pub fn equ_string(equ_val: i32) -> &'static str {
    let b_from = equ_val & (II_STEREO | II_STEREO_INV | II_NUMB | II_EQU);
    let b_type = equ_val & (IIT_ISO | IIT_NONTAUT);
    let b_eq2 = equ_val & (II_EQ2_NONTAUT | II_EQ2_ISO | II_EQ2_INV);
    let both = IIT_NONTAUT | IIT_ISO;
    // the plain table shared by several branches: `0` => "m", ISO => "M",
    // NONTAUT => "n", both => "N"
    let m_n = |t: i32, e: i32| -> &'static str {
        match (t, e) {
            (IIT_ISO, 0) | (IIT_NONTAUT, 0) => "m",
            (x, 0) if x == both => "m",
            (x, II_EQ2_ISO) if x == both => "M",
            (x, II_EQ2_NONTAUT) if x == both => "n",
            _ => "??",
        }
    };
    match b_from {
        II_STEREO => m_n(b_type, b_eq2),
        II_STEREO_INV => {
            if b_eq2 & II_EQ2_INV != 0 {
                match (b_type, b_eq2 & !II_EQ2_INV) {
                    (0 | IIT_ISO | IIT_NONTAUT, 0) => "im",
                    (0 | IIT_ISO | IIT_NONTAUT, II_EQ2_ISO) => "iM",
                    (0 | IIT_ISO | IIT_NONTAUT, II_EQ2_NONTAUT) => "in",
                    (0 | IIT_ISO | IIT_NONTAUT, e) if e == (II_EQ2_NONTAUT | II_EQ2_ISO) => "iN",
                    (t, 0) if t == both => "im",
                    (t, II_EQ2_ISO) if t == both => "iM",
                    (t, II_EQ2_NONTAUT) if t == both => "in",
                    (t, e) if t == both && e == (II_EQ2_NONTAUT | II_EQ2_ISO) => "iN",
                    _ => "??",
                }
            } else {
                m_n(b_type, b_eq2)
            }
        }
        x if x == (II_NUMB | II_STEREO_INV) => match b_type {
            0 => {
                if b_eq2 == 0 {
                    "m"
                } else {
                    "??"
                }
            }
            IIT_ISO => match b_eq2 {
                0 => "m",
                II_EQ2_INV => "im",
                II_EQ2_ISO => "M",
                _ => "??",
            },
            IIT_NONTAUT => match b_eq2 {
                0 => "m",
                II_EQ2_NONTAUT => "n",
                II_EQ2_INV => "im",
                _ => "??",
            },
            t if t == both => match b_eq2 {
                0 => "m",
                II_EQ2_ISO => "M",
                e if e == (II_EQ2_ISO | II_EQ2_INV) => "iM",
                II_EQ2_NONTAUT => "n",
                e if e == (II_EQ2_NONTAUT | II_EQ2_ISO) => "N",
                II_EQ2_INV => "im",
                e if e == (II_EQ2_NONTAUT | II_EQ2_INV) => "in",
                _ => "??",
            },
            _ => "??",
        },
        II_NUMB | II_EQU => {
            if b_type == 0 {
                "??"
            } else {
                m_n(b_type, b_eq2)
            }
        }
        _ => "??",
    }
}
