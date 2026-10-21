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

// ---------------------------------------------------------------------------
// Tags, line ends and the layers of the InChI string
// (`INCHI_BASE/src/ichiprt1.c`)

use super::ichi::{
    DIFL_F, DIFS_B_SBONDS, DIFS_H_H_ATOMS, DIFS_I_IATOMS, DIFS_M_SP3INV, DIFS_O_TRANSP,
    DIFS_P_PROTONS, DIFS_Q_CHARGE, DIFS_S_STYPE, DIFS_T_SATOMS, DifSegs, INCHI_SEGM_FILL,
    InChISort, inchi_segment_action,
};
use super::ichi_io::{
    InchiIosString, InchiIostream, inchi_ios_print_nodisplay, inchi_strbuf_printf,
    inchi_strbuf_reset, inchi_strbuf_update,
};
use super::ichican2::TAUT_YES;
use super::ichiprt2::{make_delim, make_iso_h_string};
use super::ichiprt3::{
    OUT_NN, OUT_TN, bin_aux_taut_trans, str_aux_taut_trans, str_charge2, str_connections,
    str_h_atoms, str_hill_formula, str_iso_atoms, str_iso_sp2, str_iso_sp3, str_iso_stereo_abs_inv,
    str_sp2, str_sp3, str_stereo_abs_inv,
};

pub const INCHI_BAS: i32 = 0;
pub const INCHI_REC: i32 = 1;
pub const INCHI_OUT_PLAIN_TEXT: i32 = 0x0040;
pub const INCHI_OUT_PLAIN_TEXT_COMMENTS: i32 = 0x0080;
pub const INCHI_OUT_TABBED_OUTPUT: i32 = 0x0400;
pub const INCHI_OUT_STDINCHI: i32 = 0x4000;
const INCHI_NAME: &str = "InChI";
const X_CURR_VER: &str = "1";

/// `INCHI_TAG`: plain label, plain comment, XML label, always output
#[derive(Debug)]
pub struct InchiTag(
    pub &'static str,
    pub &'static str,
    pub &'static str,
    pub i32,
);

// bits of `IDENT_LBL_BIT`
pub const IL_FIXH: i32 = 1 << 0;
pub const IL_ISOT: i32 = 1 << 1;
pub const IL_STER: i32 = 1 << 2;
pub const IL_VERS: i32 = 1 << 3;
pub const IL_FML_: i32 = 1 << 4;
pub const IL_CONN: i32 = 1 << 5;
pub const IL_ALLH: i32 = 1 << 6;
pub const IL_CHRG: i32 = 1 << 7;
pub const IL_PROT: i32 = 1 << 8;
pub const IL_DBND: i32 = 1 << 9;
pub const IL_SP3S: i32 = 1 << 10;
pub const IL_INVS: i32 = 1 << 11;
pub const IL_TYPS: i32 = 1 << 12;
pub const IL_ATMS: i32 = 1 << 13;
pub const IL_XCGA: i32 = 1 << 14;
pub const IL_FMLF: i32 = 1 << 15;
pub const IL_HFIX: i32 = 1 << 16;
pub const IL_TRNS: i32 = 1 << 17;
pub const IL_REC_: i32 = 1 << 18;

/// `IdentLbl`
pub const IDENT_LBL: [InchiTag; 19] = [
    InchiTag("/", "fixed_H", "fixed-H", 0),
    InchiTag("/", "isotopic", "isotopic", 0),
    InchiTag("/", "stereo", "stereo", 0),
    InchiTag("", "version", "version", 1),
    InchiTag("/", "formula", "formula", 1),
    InchiTag("/c", "connections", "connections", 1),
    InchiTag("/h", "H_atoms", "H", 1),
    InchiTag("/q", "charge", "charge", 1),
    InchiTag("/p", "protons", "protons", 0),
    InchiTag("/b", "dbond", "dbond", 0),
    InchiTag("/t", "sp3", "sp3", 0),
    InchiTag("/m", "sp3:inverted", "abs.inverted", 0),
    InchiTag("/s", "type (1=abs, 2=rel, 3=rac)", "type", 0),
    InchiTag("/i", "atoms", "atoms", 1),
    InchiTag("/h", "exchangeable_H", "H-isotopic", 1),
    InchiTag("/f", "formula", "formula", 1),
    InchiTag("/h", "H_fixed", "H-fixed", 1),
    InchiTag("/o", "transposition", "transposition", 0),
    InchiTag(
        "/r",
        "reconnected bond(s) to metal(s) formula",
        "formula",
        0,
    ),
];

// bits of `AUX_LBL_BIT`
pub const AL_FIXH: i32 = 1 << 0;
pub const AL_ISOT: i32 = 1 << 1;
pub const AL_STER: i32 = 1 << 2;
pub const AL_REVR: i32 = 1 << 3;
pub const AL_VERS: i32 = 1 << 4;
pub const AL_NORM: i32 = 1 << 5;
pub const AL_ANBR: i32 = 1 << 6;
pub const AL_AEQU: i32 = 1 << 7;
pub const AL_GEQU: i32 = 1 << 8;
pub const AL_SP3I: i32 = 1 << 9;
pub const AL_SP3N: i32 = 1 << 10;
pub const AL_CRV_: i32 = 1 << 11;
pub const AL_ATMR: i32 = 1 << 12;
pub const AL_BNDR: i32 = 1 << 13;
pub const AL_XYZR: i32 = 1 << 14;
pub const AL_FIXN: i32 = 1 << 15;
pub const AL_ISON: i32 = 1 << 16;
pub const AL_REC_: i32 = 1 << 17;

/// `AuxLbl`
pub const AUX_LBL: [InchiTag; 18] = [
    InchiTag("/", "fixed_H", "fixed-H", 0),
    InchiTag("/", "isotopic", "isotopic", 0),
    InchiTag("/", "abs_stereo_inverted", "stereo.abs.inverted", 0),
    InchiTag("/", "reversibility", "reversibility", 0),
    InchiTag("", "version", "version", 1),
    InchiTag("/", "normalization_type", "norm-type", 1),
    InchiTag("/N:", "original_atom_numbers", "atom.orig-nbr", 1),
    InchiTag("/E:", "atom_equivalence", "atom.equivalence", 0),
    InchiTag("/gE:", "group_equivalence", "group.equivalence", 0),
    InchiTag("/it:", "sp3", "sp3", 0),
    InchiTag("/iN:", "original_atom_numbers", "atom.orig-nbr", 0),
    InchiTag("/CRV:", "charge_radical_valence", "charges-rad-val", 0),
    InchiTag("/rA:", "atoms", "atoms", 0),
    InchiTag("/rB:", "bonds", "bonds", 0),
    InchiTag("/rC:", "xyz", "xyz", 0),
    InchiTag("/F:", "original_atom_numbers", "atom.orig-nbr", 1),
    InchiTag("/I:", "original_atom_numbers", "atom.orig-nbr", 1),
    InchiTag("/R:", "reconnected bond(s) to metal(s) part", "", 1),
];

/// `MAX_TAG_NUM`
const MAX_TAG_NUM: usize = 19;

/// `szGetTag`
pub fn sz_get_tag(tags: &[InchiTag], n_tag: i32, b_tag: i32, always: &mut i32) -> String {
    let bit_tag = |i: usize| i < tags.len() && b_tag & (1 << i) != 0;
    if 0 < n_tag && n_tag < 3 {
        let j = (0..MAX_TAG_NUM).rev().find(|&i| bit_tag(i));
        if let Some(j) = j {
            if n_tag != 2 {
                *always = tags[j].3;
            }
            return (if n_tag == 1 { tags[j].2 } else { tags[j].0 }).to_string();
        }
    } else if n_tag == 3 {
        let mut s = String::from("{");
        let (mut j, mut num) = (None, 0);
        for i in (0..MAX_TAG_NUM).filter(|&i| bit_tag(i)) {
            j = Some(i);
            if num > 0 {
                s.push(':');
            }
            num += 1;
            s.push_str(tags[i].1);
        }
        return match j {
            Some(j) => {
                s.push('}');
                *always = tags[j].3;
                format!("{}{}", tags[j].0, s)
            }
            None => "???".to_string(),
        };
    }
    "???".to_string()
}

/// the content of the buffer up to the terminating zero
pub fn strbuf_str(buf: &InchiIosString) -> String {
    match &buf.p_str {
        Some(p) => {
            let len = p.iter().position(|&c| c == 0).unwrap_or(p.len());
            String::from_utf8_lossy(&p[..len]).into_owned()
        }
        None => String::new(),
    }
}

/// `str_LineEnd`
pub fn str_line_end(
    tag: &str,
    overflow: &mut i32,
    buf: &mut InchiIosString,
    ind: i32,
    plain_text_tags: i32,
) -> i32 {
    if *overflow != 0 {
        return 1;
    }
    if ind < 0 {
        let first = buf
            .p_str
            .as_ref()
            .and_then(|p| p.first())
            .copied()
            .unwrap_or(0);
        if first != 0 || ind < -1 {
            let tag_len = if plain_text_tags != 0 { tag.len() } else { 0 };
            if tag_len > 0 {
                inchi_strbuf_update(buf, (tag_len + 4) as i32);
                let used = buf.n_used_length as usize;
                let p = buf.p_str.as_mut().expect("allocated by update");
                p.copy_within(0..used + 1, tag_len);
                p[..tag_len].copy_from_slice(tag.as_bytes());
                buf.n_used_length = p.iter().position(|&c| c == 0).unwrap_or(p.len()) as i32;
            }
        }
    }
    0
}

/// `GetSaveOptLetters`
pub fn get_save_opt_letters(save_opt_bits: u8) -> (char, char) {
    let a2p = b"ABCDEFGHIJKLMNOP";
    (
        a2p[usize::from(save_opt_bits & 0x0f)] as char,
        a2p[usize::from((save_opt_bits & 0x30) >> 4)] as char,
    )
}

/// `set_line_separators` (`TARGET_API_LIB`: the tab separator is `\n`)
pub fn set_line_separators(output_options: i32) -> (&'static str, &'static str) {
    let plain_comments = output_options & INCHI_OUT_PLAIN_TEXT_COMMENTS != 0;
    (if plain_comments { "\n" } else { "" }, "\n")
}

/// `INCHI_OUT_CTL`
#[derive(Debug, Clone, Default)]
pub struct InchiOutCtl {
    pub atom_mode: i32,
    pub taut_mode: i32,
    pub b_overflow: i32,
    pub b_always: i32,
    pub b_output_type: i32,
    pub b_out_type: i32,
    pub b_plain_text_tags: i32,
    pub b_omit_repetitions: i32,
    pub b_use_multipliers: i32,
    pub b_non_taut_non_iso_identifier_not_empty: i32,
    pub b_non_taut_iso_identifier_not_empty: i32,
    pub b_second_non_taut_pass: i32,
    pub b_tautomeric_output_allowed: i32,
    pub b_tautomeric: i32,
    pub b_non_tautomeric: i32,
    pub b_non_taut_is_identical_to_taut: i32,
    pub b_fh_tag: i32,
    pub b_rel_rac: i32,
    pub b_abc_numbers: i32,
    pub b_isotopic: i32,
    pub b_polymers: i32,
    pub i_cur_taut_mode: usize,
    pub num_components: usize,
    pub n_num_removed_protons: i32,
    pub n_tag: i32,
    pub b_tag1: i32,
    pub b_tag2: i32,
    pub b_tag3: i32,
    pub tot_len: i32,
    pub tot_len2: i32,
    pub n_cur_inchi_segment: usize,
    pub n_segm_action: i32,
    pub num_comp: [i32; 2],
    pub num_iso_h: [i32; 3],
    pub b_relative_stereo: [i32; 2],
    pub b_racemic_stereo: [i32; 2],
    pub b_isotopic_relative_stereo: [i32; 2],
    pub b_isotopic_racemic_stereo: [i32; 2],
    /// `*pSortPrintINChIFlags`
    pub sort_print_inchi_flags: i32,
    pub s_dif_segs: DifSegs,
    pub sz_tag1: String,
    pub sz_tag2: String,
    pub sz_tag3: String,
    pub n_pzz: i32,
    pub n_zy: i32,
    pub p_inchi_sort: Vec<InChISort>,
    pub p_inchi_sort2: Vec<InChISort>,
}

/// `OutputINCHI_VersionAndKind`
pub fn output_inchi_version_and_kind(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    output_options: i32,
    is_beta: bool,
    lf: &str,
) -> i32 {
    inchi_ios_print_nodisplay(Some(out_file), &format!("{lf}{INCHI_NAME}={lf}"));
    inchi_strbuf_reset(strbuf);
    inchi_strbuf_printf(strbuf, X_CURR_VER);
    if is_beta {
        inchi_strbuf_printf(strbuf, "B");
    } else if output_options & INCHI_OUT_STDINCHI != 0 {
        inchi_strbuf_printf(strbuf, "S");
    }
    inchi_ios_print_nodisplay(Some(out_file), &format!("{}{lf}", strbuf_str(strbuf)));
    0
}

/// `MergeZzInStrHillFormulaComponent`
fn merge_zz_in_str_hill_formula_component(s: &mut Vec<u8>) {
    let find = |s: &[u8], from: usize| -> Option<usize> {
        s.get(from..)?
            .windows(2)
            .position(|w| w == b"Zz")
            .map(|p| p + from)
    };
    let number = |s: &[u8], at: usize| -> (i32, usize) {
        let mut e = at;
        let mut v = 0i32;
        while e < s.len() && s[e].is_ascii_digit() {
            v = v.wrapping_mul(10).wrapping_add(i32::from(s[e] - b'0'));
            e += 1;
        }
        (v, e)
    };
    let Some(pz) = find(s, 0) else {
        return;
    };
    let pd = pz + 2;
    let (mut n, mut pz2) = (1, pd);
    if s.get(pd).is_some_and(u8::is_ascii_digit) {
        (n, pz2) = number(s, pd);
    }
    if let Some(p2) = find(s, pz2) {
        let mut n2 = 1;
        if s.get(p2 + 2).is_some_and(u8::is_ascii_digit) {
            n2 = number(s, p2 + 2).0;
        }
        n = n.wrapping_add(n2);
        s.truncate(pd);
        s.extend_from_slice(n.to_string().as_bytes());
    }
}

/// `MergeZzInHillFormula`
pub fn merge_zz_in_hill_formula(strbuf: &mut InchiIosString) -> i32 {
    if strbuf.p_str.is_none() || strbuf.n_used_length < 1 {
        return 0;
    }
    let text: Vec<u8> = {
        let p = strbuf.p_str.as_ref().expect("checked");
        let len = p.iter().position(|&c| c == 0).unwrap_or(p.len());
        p[..len].to_vec()
    };
    inchi_strbuf_reset(strbuf);
    let mut p = 0usize;
    loop {
        let pend = text[p..]
            .iter()
            .position(|&c| c == b'.')
            .map_or(text.len(), |x| x + p);
        let mut comp = text[p..pend].to_vec();
        merge_zz_in_str_hill_formula_component(&mut comp);
        let delim = text.get(pend).copied().unwrap_or(0);
        let mut out = String::from_utf8_lossy(&comp).into_owned();
        out.push(delim as char);
        inchi_strbuf_printf(strbuf, &out);
        if delim == 0 {
            break;
        }
        p = pend + 1;
    }
    0
}

/// `OutputINCHI_MainLayerFormula`
pub fn output_inchi_main_layer_formula(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    num_components2: [i32; 2],
    basic_or_reconnected: i32,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    if num_components2[0] != 0 || num_components2[1] != 0 {
        io.b_tag1 = if basic_or_reconnected == INCHI_REC {
            IL_REC_
        } else {
            IL_FML_
        };
        io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = str_hill_formula(
            &io.p_inchi_sort,
            strbuf,
            &mut io.b_overflow,
            io.b_out_type,
            io.num_components,
            io.b_use_multipliers != 0,
        );
        if str_line_end(&io.sz_tag1, &mut io.b_overflow, strbuf, -1, 1) != 0 {
            return 1;
        }
        if io.n_pzz > 0 && io.n_zy > 0 && merge_zz_in_hill_formula(strbuf) != 0 {
            return -1;
        }
        inchi_ios_print_nodisplay(Some(out_file), &format!("{}{lf}", strbuf_str(strbuf)));
    }
    0
}

/// `OutputINCHI_MainLayerConnections`
pub fn output_inchi_main_layer_connections(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    io.b_tag1 = IL_CONN;
    io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
    inchi_strbuf_reset(strbuf);
    io.tot_len = 0;
    io.tot_len2 = str_connections(
        &io.p_inchi_sort,
        strbuf,
        &mut io.b_overflow,
        io.b_out_type,
        io.atom_mode,
        io.num_components,
        io.b_use_multipliers != 0,
    );
    if io.tot_len != io.tot_len2 {
        io.tot_len = io.tot_len2;
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -2,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        inchi_ios_print_nodisplay(Some(out_file), &format!("{}{lf}", strbuf_str(strbuf)));
    }
    0
}

/// `OutputINCHI_MainLayerHydrogens`
pub fn output_inchi_main_layer_hydrogens(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    if inchi_segment_action(io.s_dif_segs[io.n_cur_inchi_segment][DIFS_H_H_ATOMS])
        == INCHI_SEGM_FILL
    {
        io.b_tag1 = IL_ALLH;
        io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        io.tot_len2 = str_h_atoms(
            &io.p_inchi_sort,
            strbuf,
            &mut io.b_overflow,
            io.b_out_type,
            io.atom_mode,
            io.taut_mode,
            io.num_components,
            io.b_use_multipliers != 0,
        );
        if io.tot_len != io.tot_len2 {
            io.tot_len = io.tot_len2;
            if str_line_end(&io.sz_tag1, &mut io.b_overflow, strbuf, -2, 1) != 0 {
                return 1;
            }
            inchi_ios_print_nodisplay(Some(out_file), &format!("{}{lf}", strbuf_str(strbuf)));
        }
    }
    0
}

/// `OutputINCHI_ChargeAndRemovedAddedProtonsLayers`
pub fn output_inchi_charge_and_removed_added_protons_layers(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    io.n_segm_action = inchi_segment_action(io.s_dif_segs[io.n_cur_inchi_segment][DIFS_Q_CHARGE]);
    if io.n_segm_action != 0 {
        io.b_tag1 = IL_CHRG | io.b_fh_tag;
        io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        if io.n_segm_action == INCHI_SEGM_FILL {
            io.tot_len = str_charge2(
                &io.p_inchi_sort,
                &io.p_inchi_sort2,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.num_components,
                io.b_second_non_taut_pass != 0,
                io.b_omit_repetitions != 0,
                io.b_use_multipliers != 0,
            );
            io.b_non_taut_non_iso_identifier_not_empty += io.b_second_non_taut_pass;
        }
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -io.n_segm_action,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        inchi_ios_print_nodisplay(Some(out_file), &format!("{}{lf}", strbuf_str(strbuf)));
    }
    if io.i_cur_taut_mode == TAUT_YES && io.b_second_non_taut_pass == 0 {
        io.n_segm_action =
            inchi_segment_action(io.s_dif_segs[io.n_cur_inchi_segment][DIFS_P_PROTONS]);
        if io.n_segm_action != 0 {
            io.b_tag1 = IL_PROT | io.b_fh_tag;
            io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
            inchi_strbuf_reset(strbuf);
            io.tot_len = 0;
            inchi_strbuf_printf(strbuf, &format!("{:+}", io.n_num_removed_protons));
            if str_line_end(
                &io.sz_tag1,
                &mut io.b_overflow,
                strbuf,
                -io.n_segm_action,
                io.b_plain_text_tags,
            ) != 0
            {
                return 1;
            }
            inchi_ios_print_nodisplay(Some(out_file), &format!("{}{lf}", strbuf_str(strbuf)));
        } else if io.b_plain_text_tags == 1 {
            inchi_ios_print_nodisplay(Some(out_file), "/");
        }
    }
    0
}

pub const FLAG_SORT_PRINT_TRANSPOS_BAS: i32 = 1;
pub const FLAG_SORT_PRINT_TRANSPOS_REC: i32 = 2;
const X_REL: &str = "2";
const X_RAC: &str = "3";
const X_ABS: &str = "1";

fn seg_action(io: &InchiOutCtl, seg: usize) -> i32 {
    inchi_segment_action(io.s_dif_segs[io.n_cur_inchi_segment][seg])
}

fn print_layer(out_file: &mut InchiIostream, strbuf: &InchiIosString, lf: &str) {
    inchi_ios_print_nodisplay(Some(out_file), &format!("{}{lf}", strbuf_str(strbuf)));
}

fn any_stereo_segment(io: &InchiOutCtl) -> bool {
    seg_action(io, DIFS_B_SBONDS) != 0
        || seg_action(io, DIFS_T_SATOMS) != 0
        || seg_action(io, DIFS_M_SP3INV) != 0
        || seg_action(io, DIFS_S_STYPE) != 0
}

/// `OutputINCHI_StereoLayer`
pub fn output_inchi_stereo_layer(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    let slash = |out_file: &mut InchiIostream, io: &InchiOutCtl, s: &str| {
        if io.b_plain_text_tags == 1 {
            inchi_ios_print_nodisplay(Some(out_file), s);
        }
    };
    if !any_stereo_segment(io) {
        slash(out_file, io, "////");
        return 0;
    }
    io.b_tag1 = IL_STER | io.b_fh_tag;
    io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
    let cur = io.i_cur_taut_mode;

    // sp2
    io.n_segm_action = seg_action(io, DIFS_B_SBONDS);
    if io.n_segm_action != 0 {
        io.b_tag2 = io.b_tag1 | IL_DBND;
        io.sz_tag2 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag2, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        if io.n_segm_action == INCHI_SEGM_FILL {
            io.tot_len = str_sp2(
                &io.p_inchi_sort,
                &io.p_inchi_sort2,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.taut_mode,
                io.num_components,
                io.b_second_non_taut_pass != 0,
                io.b_omit_repetitions != 0,
                io.b_use_multipliers != 0,
            );
            io.b_non_taut_non_iso_identifier_not_empty += io.b_second_non_taut_pass;
        }
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -io.n_segm_action,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        print_layer(out_file, strbuf, lf);
    } else {
        slash(out_file, io, "/");
    }

    // sp3
    io.n_segm_action = seg_action(io, DIFS_T_SATOMS);
    if io.n_segm_action != 0 {
        io.b_rel_rac = i32::from(io.b_relative_stereo[cur] != 0 || io.b_racemic_stereo[cur] != 0);
        io.b_tag2 = io.b_tag1 | IL_SP3S;
        io.sz_tag2 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag2, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        if io.n_segm_action == INCHI_SEGM_FILL {
            io.tot_len = str_sp3(
                &io.p_inchi_sort,
                &io.p_inchi_sort2,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.taut_mode,
                io.num_components,
                io.b_second_non_taut_pass != 0,
                io.b_omit_repetitions != 0,
                io.b_use_multipliers != 0,
            );
            io.b_non_taut_non_iso_identifier_not_empty += io.b_second_non_taut_pass;
        }
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -io.n_segm_action,
            io.b_plain_text_tags,
        ) != 0
        {
            return 2;
        }
        print_layer(out_file, strbuf, lf);
    } else {
        slash(out_file, io, "/");
    }

    // inverted absolute stereo
    io.n_segm_action = seg_action(io, DIFS_M_SP3INV);
    if io.n_segm_action != 0 {
        io.b_tag2 = io.b_tag1 | IL_INVS;
        io.sz_tag2 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag2, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        if io.n_segm_action == INCHI_SEGM_FILL {
            io.tot_len = str_stereo_abs_inv(
                &io.p_inchi_sort,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.num_components,
            );
            io.b_non_taut_non_iso_identifier_not_empty += io.b_second_non_taut_pass;
        }
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -io.n_segm_action,
            io.b_plain_text_tags,
        ) != 0
        {
            return 3;
        }
        print_layer(out_file, strbuf, lf);
    } else {
        slash(out_file, io, "/");
    }

    // stereo type
    io.n_segm_action = seg_action(io, DIFS_S_STYPE);
    if io.n_segm_action != 0 {
        let p_stereo = if io.b_relative_stereo[cur] != 0 {
            X_REL
        } else if io.b_racemic_stereo[cur] != 0 {
            X_RAC
        } else {
            X_ABS
        };
        io.b_tag2 = io.b_tag1 | IL_TYPS;
        io.sz_tag2 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag2, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        if io.n_segm_action == INCHI_SEGM_FILL {
            io.tot_len += make_delim(Some(p_stereo), strbuf, &mut io.b_overflow);
            io.b_non_taut_non_iso_identifier_not_empty += io.b_second_non_taut_pass;
        }
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -io.n_segm_action,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        print_layer(out_file, strbuf, lf);
    }
    slash(out_file, io, "/");
    0
}

/// `OutputINCHI_IsotopicLayer`
pub fn output_inchi_isotopic_layer(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    basic_or_reconnected: i32,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    let slash = |out_file: &mut InchiIostream, io: &InchiOutCtl, s: &str| {
        if io.b_plain_text_tags == 1 {
            inchi_ios_print_nodisplay(Some(out_file), s);
        }
    };
    let cur = io.i_cur_taut_mode;
    if seg_action(io, DIFS_I_IATOMS) != 0 {
        io.b_tag1 = IL_ISOT | io.b_fh_tag;
        io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
        io.n_segm_action = seg_action(io, DIFS_I_IATOMS);
        if io.n_segm_action != 0 {
            io.b_tag2 = io.b_tag1 | IL_ATMS;
            io.sz_tag2 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag2, &mut io.b_always);
            inchi_strbuf_reset(strbuf);
            io.tot_len = 0;
            if io.n_segm_action == INCHI_SEGM_FILL {
                io.tot_len2 = str_iso_atoms(
                    &io.p_inchi_sort,
                    &io.p_inchi_sort2,
                    strbuf,
                    &mut io.b_overflow,
                    io.b_out_type,
                    io.taut_mode,
                    io.num_components,
                    io.b_abc_numbers != 0,
                    io.b_second_non_taut_pass != 0,
                    io.b_omit_repetitions != 0,
                    io.b_use_multipliers != 0,
                );
                io.b_non_taut_iso_identifier_not_empty += io.b_second_non_taut_pass;
            } else {
                io.tot_len2 = io.tot_len;
            }
            io.tot_len = io.tot_len2;
            if str_line_end(
                &io.sz_tag2,
                &mut io.b_overflow,
                strbuf,
                -io.n_segm_action,
                io.b_plain_text_tags,
            ) != 0
            {
                return 1;
            }
            print_layer(out_file, strbuf, lf);
        }

        // exchangeable isotopic H (mobile H only)
        io.n_segm_action = seg_action(io, DIFS_H_H_ATOMS);
        if io.n_segm_action != 0 {
            io.b_tag2 = io.b_tag1 | IL_XCGA;
            io.sz_tag2 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag2, &mut io.b_always);
            inchi_strbuf_reset(strbuf);
            io.tot_len = 0;
            io.tot_len +=
                make_iso_h_string(&io.num_iso_h, strbuf, io.taut_mode, &mut io.b_overflow);
            io.b_non_taut_iso_identifier_not_empty += io.b_second_non_taut_pass;
            if str_line_end(
                &io.sz_tag2,
                &mut io.b_overflow,
                strbuf,
                -io.n_segm_action,
                io.b_plain_text_tags,
            ) != 0
            {
                return 2;
            }
            print_layer(out_file, strbuf, lf);
        }

        // isotopic stereo
        if any_stereo_segment(io) {
            io.b_tag2 = io.b_tag1 | IL_STER;
            io.sz_tag2 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag2, &mut io.b_always);

            io.n_segm_action = seg_action(io, DIFS_B_SBONDS);
            if io.n_segm_action != 0 {
                io.b_tag3 = io.b_tag2 | IL_DBND;
                io.sz_tag3 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag3, &mut io.b_always);
                inchi_strbuf_reset(strbuf);
                io.tot_len = 0;
                if io.n_segm_action == INCHI_SEGM_FILL {
                    io.tot_len = str_iso_sp2(
                        &io.p_inchi_sort,
                        &io.p_inchi_sort2,
                        strbuf,
                        &mut io.b_overflow,
                        io.b_out_type,
                        io.taut_mode,
                        io.num_components,
                        io.b_second_non_taut_pass != 0,
                        io.b_omit_repetitions != 0,
                        io.b_use_multipliers != 0,
                    );
                    io.b_non_taut_iso_identifier_not_empty += io.b_second_non_taut_pass;
                }
                if str_line_end(
                    &io.sz_tag3,
                    &mut io.b_overflow,
                    strbuf,
                    -io.n_segm_action,
                    io.b_plain_text_tags,
                ) != 0
                {
                    return 3;
                }
                print_layer(out_file, strbuf, lf);
            } else {
                slash(out_file, io, "/");
            }

            io.n_segm_action = seg_action(io, DIFS_T_SATOMS);
            if io.n_segm_action != 0 {
                io.b_rel_rac = i32::from(
                    io.b_isotopic_relative_stereo[cur] != 0
                        || io.b_isotopic_racemic_stereo[cur] != 0,
                );
                io.b_tag3 = io.b_tag2 | IL_SP3S;
                io.sz_tag3 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag3, &mut io.b_always);
                inchi_strbuf_reset(strbuf);
                io.tot_len = 0;
                if io.n_segm_action == INCHI_SEGM_FILL {
                    io.tot_len = str_iso_sp3(
                        &io.p_inchi_sort,
                        &io.p_inchi_sort2,
                        strbuf,
                        &mut io.b_overflow,
                        io.b_out_type,
                        io.taut_mode,
                        io.num_components,
                        io.b_second_non_taut_pass != 0,
                        io.b_omit_repetitions != 0,
                        io.b_use_multipliers != 0,
                    );
                    io.b_non_taut_iso_identifier_not_empty += io.b_second_non_taut_pass;
                }
                if str_line_end(
                    &io.sz_tag3,
                    &mut io.b_overflow,
                    strbuf,
                    -io.n_segm_action,
                    io.b_plain_text_tags,
                ) != 0
                {
                    return 5;
                }
                print_layer(out_file, strbuf, lf);
            } else {
                slash(out_file, io, "/");
            }

            io.n_segm_action = seg_action(io, DIFS_M_SP3INV);
            if io.n_segm_action != 0 {
                io.b_tag3 = io.b_tag2 | IL_INVS;
                io.sz_tag3 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag3, &mut io.b_always);
                inchi_strbuf_reset(strbuf);
                io.tot_len = 0;
                if io.n_segm_action == INCHI_SEGM_FILL {
                    io.tot_len = str_iso_stereo_abs_inv(
                        &io.p_inchi_sort,
                        strbuf,
                        &mut io.b_overflow,
                        io.b_out_type,
                        io.num_components,
                    );
                    io.b_non_taut_iso_identifier_not_empty += io.b_second_non_taut_pass;
                }
                if str_line_end(
                    &io.sz_tag3,
                    &mut io.b_overflow,
                    strbuf,
                    -io.n_segm_action,
                    io.b_plain_text_tags,
                ) != 0
                {
                    return 5;
                }
                print_layer(out_file, strbuf, lf);
            } else {
                slash(out_file, io, "/");
            }

            io.n_segm_action = seg_action(io, DIFS_S_STYPE);
            if io.n_segm_action != 0 {
                let p_stereo = if io.b_isotopic_relative_stereo[cur] != 0 {
                    X_REL
                } else if io.b_isotopic_racemic_stereo[cur] != 0 {
                    X_RAC
                } else {
                    X_ABS
                };
                io.b_tag3 = io.b_tag2 | IL_TYPS;
                io.sz_tag3 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag3, &mut io.b_always);
                inchi_strbuf_reset(strbuf);
                io.tot_len = 0;
                if io.n_segm_action == INCHI_SEGM_FILL {
                    io.tot_len += make_delim(Some(p_stereo), strbuf, &mut io.b_overflow);
                    io.b_non_taut_iso_identifier_not_empty += io.b_second_non_taut_pass;
                }
                if str_line_end(
                    &io.sz_tag3,
                    &mut io.b_overflow,
                    strbuf,
                    -io.n_segm_action,
                    io.b_plain_text_tags,
                ) != 0
                {
                    return 6;
                }
                print_layer(out_file, strbuf, lf);
            }
            slash(out_file, io, "/");
        } else {
            slash(out_file, io, "////");
        }
    } else {
        slash(out_file, io, "///");
        slash(out_file, io, "//");
    }

    // CANON_FIXH_TRANS: the transposition of the non-tautomeric components
    if io.b_out_type == OUT_NN
        && io.b_output_type == OUT_TN
        && io.b_second_non_taut_pass != 0
        && inchi_segment_action(io.s_dif_segs[DIFL_F][DIFS_O_TRANSP]) == INCHI_SEGM_FILL
    {
        let trans = bin_aux_taut_trans(
            &io.p_inchi_sort,
            &io.p_inchi_sort2,
            io.b_out_type,
            io.num_components,
        );
        if trans.is_some() {
            io.b_tag1 = IL_TRNS | io.b_fh_tag;
            io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
            inchi_strbuf_reset(strbuf);
            io.tot_len = str_aux_taut_trans(
                trans,
                strbuf,
                &mut io.b_overflow,
                io.taut_mode,
                io.num_components,
            );
            io.b_non_taut_iso_identifier_not_empty += io.b_second_non_taut_pass;
            if str_line_end(
                &io.sz_tag1,
                &mut io.b_overflow,
                strbuf,
                -1,
                io.b_plain_text_tags,
            ) != 0
            {
                return 7;
            }
            print_layer(out_file, strbuf, lf);
            io.sort_print_inchi_flags |= if basic_or_reconnected == INCHI_BAS {
                FLAG_SORT_PRINT_TRANSPOS_BAS
            } else {
                FLAG_SORT_PRINT_TRANSPOS_REC
            };
        } else {
            slash(out_file, io, "/");
        }
    }
    0
}
