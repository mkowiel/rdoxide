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
    DIFL_F, DIFS_B_SBONDS, DIFS_F_FORMULA, DIFS_H_H_ATOMS, DIFS_I_IATOMS, DIFS_M_SP3INV,
    DIFS_O_TRANSP, DIFS_P_PROTONS, DIFS_Q_CHARGE, DIFS_S_STYPE, DIFS_T_SATOMS, DifSegs,
    INCHI_SEGM_FILL, InChISort, inchi_segment_action,
};
use super::ichi_io::{
    InchiIosString, InchiIostream, inchi_ios_print, inchi_ios_print_nodisplay, inchi_strbuf_printf,
    inchi_strbuf_reset, inchi_strbuf_update,
};
use super::ichican2::TAUT_YES;
use super::ichiprt2::{make_delim, make_iso_h_string};
use super::ichiprt3::{
    OUT_NN, OUT_TN, bin_aux_taut_trans, str_aux_taut_trans, str_charge2, str_connections,
    str_fixed_h_atoms, str_h_atoms, str_hill_formula, str_hill_formula2, str_iso_atoms,
    str_iso_sp2, str_iso_sp3, str_iso_stereo_abs_inv, str_sp2, str_sp3, str_stereo_abs_inv,
};

pub const INCHI_BAS: i32 = 0;
pub const INCHI_REC: i32 = 1;
pub const INCHI_OUT_PLAIN_TEXT: i32 = 0x0040;
pub const INCHI_OUT_PLAIN_TEXT_COMMENTS: i32 = 0x0080;
pub const INCHI_OUT_WINCHI_WINDOW: i32 = 0x0200;
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
    pub b_atom_equ: [i32; 2],
    pub b_taut_equ: [i32; 2],
    pub b_inv_stereo: [i32; 2],
    pub b_inv_stereo_orig_numb: [i32; 2],
    pub b_isotopic_orig_numb: [i32; 2],
    pub b_isotopic_atom_equ: [i32; 2],
    pub b_isotopic_taut_equ: [i32; 2],
    pub b_inv_isotopic_stereo: [i32; 2],
    pub b_inv_isotopic_stereo_orig_numb: [i32; 2],
    pub b_ign_uu_sp3_iso: [i32; 2],
    pub b_ign_uu_sp2_iso: [i32; 2],
    pub b_charges_rad_val: [i32; 2],
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
    /// `pINChISortTautAndNonTaut[TAUT_NUM]`
    pub p_inchi_sort_taut_and_non_taut: [Vec<InChISort>; 2],
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
pub const FLAG_SORT_PRINT_NO_NFIX_H_BAS: i32 = 4;
pub const FLAG_SORT_PRINT_NO_NFIX_H_REC: i32 = 8;
pub const FLAG_SORT_PRINT_NO_IFIX_H_BAS: i32 = 16;
pub const FLAG_SORT_PRINT_NO_IFIX_H_REC: i32 = 32;
const X_REL: &str = "2";
const X_RAC: &str = "3";
const TAUT_NON: usize = 0;
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

/// `OutputINCHI_FixedHLayerWithSublayers`; `then_goto_repeat` is set when
/// the caller has to run the layers again for the non-tautomeric pass
pub fn output_inchi_fixed_h_layer_with_sublayers(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    basic_or_reconnected: i32,
    io: &mut InchiOutCtl,
    lf: &str,
    then_goto_repeat: &mut i32,
) -> i32 {
    let bas = basic_or_reconnected == INCHI_BAS;
    *then_goto_repeat = 0;
    if io.b_out_type == OUT_TN
        && io.b_second_non_taut_pass == 0
        && io.b_non_taut_is_identical_to_taut != 0
        && io.b_tautomeric != 0
        && io.b_non_tautomeric != 0
    {
        io.sort_print_inchi_flags |= if bas {
            FLAG_SORT_PRINT_NO_NFIX_H_BAS
        } else {
            FLAG_SORT_PRINT_NO_NFIX_H_REC
        };
        io.sort_print_inchi_flags |= if bas {
            FLAG_SORT_PRINT_NO_IFIX_H_BAS
        } else {
            FLAG_SORT_PRINT_NO_IFIX_H_REC
        };
    }

    if io.b_out_type == OUT_TN
        && io.b_non_taut_is_identical_to_taut == 0
        && inchi_segment_action(io.s_dif_segs[DIFL_F][DIFS_F_FORMULA]) != 0
    {
        io.b_out_type = OUT_NN;
        io.i_cur_taut_mode = TAUT_NON;
        io.p_inchi_sort = io.p_inchi_sort_taut_and_non_taut[TAUT_NON].clone();
        io.b_second_non_taut_pass = 1;
        io.n_cur_inchi_segment = DIFL_F;
        io.num_components = io.num_comp[io.i_cur_taut_mode] as usize;
        io.b_fh_tag = IL_FIXH;
        io.b_tag1 = io.b_fh_tag;
        io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
        io.b_tag1 = IL_FMLF | io.b_fh_tag;
        io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        io.n_segm_action = seg_action(io, DIFS_F_FORMULA);
        if io.n_segm_action == INCHI_SEGM_FILL {
            io.tot_len2 = str_hill_formula2(
                &io.p_inchi_sort,
                &io.p_inchi_sort2,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.num_components,
                io.b_use_multipliers != 0,
            );
            if io.n_pzz > 0 && io.n_zy > 0 {
                merge_zz_in_hill_formula(strbuf);
            }
            io.b_non_taut_non_iso_identifier_not_empty += io.b_second_non_taut_pass;
        } else {
            io.tot_len2 = io.tot_len;
        }
        io.tot_len = io.tot_len2;
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
        print_layer(out_file, strbuf, lf);

        io.n_segm_action = seg_action(io, DIFS_H_H_ATOMS);
        if io.n_segm_action == INCHI_SEGM_FILL {
            io.b_tag1 = IL_HFIX | io.b_fh_tag;
            io.sz_tag1 = sz_get_tag(&IDENT_LBL, io.n_tag, io.b_tag1, &mut io.b_always);
            inchi_strbuf_reset(strbuf);
            io.tot_len = 0;
            io.tot_len2 = str_fixed_h_atoms(
                &io.p_inchi_sort,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.atom_mode,
                io.num_components,
                io.b_use_multipliers != 0,
            );
            io.tot_len = io.tot_len2;
            if str_line_end(
                &io.sz_tag1,
                &mut io.b_overflow,
                strbuf,
                -io.n_segm_action,
                io.b_plain_text_tags,
            ) != 0
            {
                return 2;
            }
            print_layer(out_file, strbuf, lf);
            io.b_non_taut_non_iso_identifier_not_empty += io.b_second_non_taut_pass;
        }
        *then_goto_repeat = 1;
        return 0;
    }

    if io.b_out_type == OUT_NN && io.b_output_type == OUT_TN && io.b_second_non_taut_pass != 0 {
        io.b_out_type = OUT_TN;
        io.i_cur_taut_mode = TAUT_YES;
        io.p_inchi_sort = io.p_inchi_sort_taut_and_non_taut[TAUT_YES].clone();
        io.b_second_non_taut_pass = 0;
        io.num_components = io.num_comp[io.i_cur_taut_mode] as usize;
        if io.b_non_taut_non_iso_identifier_not_empty == 0 {
            io.sort_print_inchi_flags |= if bas {
                FLAG_SORT_PRINT_NO_NFIX_H_BAS
            } else {
                FLAG_SORT_PRINT_NO_NFIX_H_REC
            };
        }
        if io.b_non_taut_iso_identifier_not_empty == 0 {
            io.sort_print_inchi_flags |= if bas {
                FLAG_SORT_PRINT_NO_IFIX_H_BAS
            } else {
                FLAG_SORT_PRINT_NO_IFIX_H_REC
            };
        }
        io.b_fh_tag = 0;
    }
    0
}

use super::ichiprt3 as p3;

fn aux_print(out_file: &mut InchiIostream, text: &str) {
    inchi_ios_print(Some(out_file), text);
}

fn aux_tag(io: &mut InchiOutCtl, bits: i32) -> String {
    sz_get_tag(&AUX_LBL, io.n_tag, bits, &mut io.b_always)
}

/// `OutputAUXINFO_HeaderAndNormalization_type`
pub fn output_auxinfo_header_and_normalization_type(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    output_options: i32,
    basic_or_reconnected: i32,
    num_components2: [i32; 2],
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    if basic_or_reconnected == INCHI_BAS {
        inchi_strbuf_printf(strbuf, "AuxInfo=");
        let nl = if output_options & INCHI_OUT_WINCHI_WINDOW != 0 {
            "\n"
        } else {
            ""
        };
        aux_print(out_file, &format!("{nl}{}{lf}", strbuf_str(strbuf)));
        io.b_tag1 = AL_VERS;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        inchi_strbuf_printf(strbuf, X_CURR_VER);
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    } else if basic_or_reconnected == INCHI_REC {
        io.b_tag1 = AL_REC_;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        aux_print(out_file, &format!("{}{lf}", io.sz_tag1));
    }
    if num_components2[0] != 0 || num_components2[1] != 0 {
        io.b_tag1 = AL_NORM;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        inchi_strbuf_reset(strbuf);
        io.tot_len = 0;
        let v = if io.b_tautomeric != 0 && io.b_tautomeric_output_allowed != 0 {
            io.b_tautomeric
        } else {
            0
        };
        inchi_strbuf_printf(strbuf, &v.to_string());
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    }
    0
}

/// `OutputAUXINFO_OriginalNumbersAndEquivalenceClasses`
pub fn output_auxinfo_original_numbers_and_equivalence_classes(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    num_components2: [i32; 2],
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    if num_components2[0] != 0 || num_components2[1] != 0 {
        io.b_tag1 = (if io.b_second_non_taut_pass != 0 {
            AL_FIXN
        } else {
            AL_ANBR
        }) | io.b_fh_tag;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_numb(
            &io.p_inchi_sort,
            strbuf,
            &mut io.b_overflow,
            io.b_out_type,
            io.taut_mode,
            io.num_components,
            io.b_second_non_taut_pass != 0,
            io.b_omit_repetitions != 0,
        );
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    }
    if io.b_atom_equ[io.i_cur_taut_mode] != 0 {
        io.b_tag1 = AL_AEQU | io.b_fh_tag;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_equ(
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
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    } else if io.b_plain_text_tags == 1 {
        aux_print(out_file, "/");
    }
    0
}

/// `OutputAUXINFO_TautomericGroupsEquivalence`
pub fn output_auxinfo_tautomeric_groups_equivalence(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    io: &mut InchiOutCtl,
) -> i32 {
    if io.b_tautomeric_output_allowed != 0
        && io.b_tautomeric != 0
        && io.b_taut_equ[io.i_cur_taut_mode] != 0
        && io.b_second_non_taut_pass == 0
    {
        io.b_tag1 = AL_GEQU | io.b_fh_tag;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_tgroup_equ(
            &io.p_inchi_sort,
            strbuf,
            &mut io.b_overflow,
            io.b_out_type,
            io.taut_mode,
            io.num_components,
            io.b_use_multipliers != 0,
        );
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &strbuf_str(strbuf));
    } else if io.b_tautomeric_output_allowed != 0
        && io.b_tautomeric != 0
        && io.b_plain_text_tags == 1
    {
        aux_print(out_file, "/");
    }
    0
}

/// `OutputAUXINFO_Stereo`
pub fn output_auxinfo_stereo(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    let cur = io.i_cur_taut_mode;
    if io.b_inv_stereo[cur] != 0 {
        io.b_tag1 = AL_STER | io.b_fh_tag;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        io.b_tag2 = io.b_tag1 | AL_SP3I;
        io.sz_tag2 = aux_tag(io, io.b_tag2);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_inv_sp3(
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
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
        if io.b_inv_stereo_orig_numb[cur] != 0 {
            io.b_tag2 = io.b_tag1 | AL_SP3N;
            io.sz_tag2 = aux_tag(io, io.b_tag2);
            inchi_strbuf_reset(strbuf);
            io.tot_len = p3::str_aux_inv_sp3_numb(
                &io.p_inchi_sort,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.taut_mode,
                io.num_components,
                io.b_second_non_taut_pass != 0,
                io.b_omit_repetitions != 0,
            );
            if str_line_end(
                &io.sz_tag2,
                &mut io.b_overflow,
                strbuf,
                -1,
                io.b_plain_text_tags,
            ) != 0
            {
                return 1;
            }
            aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
        } else if io.b_plain_text_tags == 1 {
            aux_print(out_file, "/");
        }
    } else if io.b_plain_text_tags == 1 {
        aux_print(out_file, "//");
    }
    0
}

/// `OutputAUXINFO_IsotopicInfo`
pub fn output_auxinfo_isotopic_info(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    basic_or_reconnected: i32,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    let cur = io.i_cur_taut_mode;
    let no_ifix = io.b_second_non_taut_pass != 0
        && io.sort_print_inchi_flags
            & (if basic_or_reconnected == INCHI_BAS {
                FLAG_SORT_PRINT_NO_IFIX_H_BAS
            } else {
                FLAG_SORT_PRINT_NO_IFIX_H_REC
            })
            != 0;
    if !(io.b_isotopic != 0
        && !no_ifix
        && (io.b_isotopic_orig_numb[cur] != 0
            || io.b_isotopic_atom_equ[cur] != 0
            || (io.b_tautomeric_output_allowed != 0
                && io.b_tautomeric != 0
                && io.b_isotopic_taut_equ[cur] != 0)
            || (io.b_inv_isotopic_stereo[cur] != 0 && io.b_ign_uu_sp3_iso[cur] != 0)
            || io.b_ign_uu_sp2_iso[cur] != 0))
    {
        return 0;
    }
    io.b_tag1 = AL_ISOT | io.b_fh_tag;
    io.sz_tag1 = aux_tag(io, io.b_tag1);
    inchi_strbuf_reset(strbuf);
    io.b_tag2 = io.b_tag1 | AL_ISON;
    io.sz_tag2 = aux_tag(io, io.b_tag2);
    if io.b_isotopic_orig_numb[cur] != 0 {
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_iso_numb(
            &io.p_inchi_sort,
            strbuf,
            &mut io.b_overflow,
            io.b_out_type,
            io.taut_mode,
            io.num_components,
            io.b_second_non_taut_pass != 0,
            io.b_omit_repetitions != 0,
        );
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    } else {
        aux_print(out_file, &format!("{}{lf}", io.sz_tag2));
    }

    if io.b_isotopic_atom_equ[cur] != 0 {
        io.b_tag2 = io.b_tag1 | AL_AEQU;
        io.sz_tag2 = aux_tag(io, io.b_tag2);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_iso_equ(
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
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -2,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    } else if io.b_plain_text_tags == 1 {
        aux_print(out_file, "/");
    }

    if io.b_tautomeric_output_allowed != 0
        && io.b_tautomeric != 0
        && io.b_isotopic_taut_equ[cur] != 0
    {
        io.b_tag2 = io.b_tag1 | AL_GEQU;
        io.sz_tag2 = aux_tag(io, io.b_tag2);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_iso_tgroup_equ(
            &io.p_inchi_sort,
            strbuf,
            &mut io.b_overflow,
            io.b_out_type,
            io.taut_mode,
            io.num_components,
            io.b_omit_repetitions != 0,
            io.b_use_multipliers != 0,
        );
        if str_line_end(
            &io.sz_tag2,
            &mut io.b_overflow,
            strbuf,
            -2,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    } else if io.b_tautomeric_output_allowed != 0
        && io.b_tautomeric != 0
        && io.b_plain_text_tags == 1
    {
        aux_print(out_file, "/");
    }

    if io.b_inv_isotopic_stereo[cur] != 0 {
        io.b_tag2 = io.b_tag1 | AL_STER;
        io.sz_tag2 = aux_tag(io, io.b_tag2);
        io.b_tag3 = io.b_tag2 | AL_SP3I;
        io.sz_tag3 = aux_tag(io, io.b_tag3);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_inv_iso_sp3(
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
        if str_line_end(
            &io.sz_tag3,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &strbuf_str(strbuf));
        if io.b_inv_isotopic_stereo_orig_numb[cur] != 0 {
            io.b_tag3 = io.b_tag2 | AL_SP3N;
            io.sz_tag3 = aux_tag(io, io.b_tag3);
            inchi_strbuf_reset(strbuf);
            io.tot_len = p3::str_aux_inv_iso_sp3_numb(
                &io.p_inchi_sort,
                strbuf,
                &mut io.b_overflow,
                io.b_out_type,
                io.taut_mode,
                io.num_components,
                io.b_second_non_taut_pass != 0,
                io.b_omit_repetitions != 0,
            );
            if str_line_end(
                &io.sz_tag3,
                &mut io.b_overflow,
                strbuf,
                -1,
                io.b_plain_text_tags,
            ) != 0
            {
                return 1;
            }
            aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
        } else if io.b_plain_text_tags == 1 {
            aux_print(out_file, "/");
        }
    } else if io.b_plain_text_tags == 1 {
        aux_print(out_file, "//");
    }
    0
}

/// `OutputAUXINFO_ChargesRadicalsAndUnusualValences`
pub fn output_auxinfo_charges_radicals_and_unusual_valences(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    io: &mut InchiOutCtl,
    lf: &str,
) -> i32 {
    if io.b_second_non_taut_pass == 0 && io.b_charges_rad_val[io.i_cur_taut_mode] != 0 {
        io.b_tag1 = AL_CRV_ | io.b_fh_tag;
        io.sz_tag1 = aux_tag(io, io.b_tag1);
        inchi_strbuf_reset(strbuf);
        io.tot_len = p3::str_aux_charge_rad_val(
            &io.p_inchi_sort,
            strbuf,
            &mut io.b_overflow,
            io.b_out_type,
            io.taut_mode,
            io.num_components,
            io.b_use_multipliers != 0,
        );
        if str_line_end(
            &io.sz_tag1,
            &mut io.b_overflow,
            strbuf,
            -1,
            io.b_plain_text_tags,
        ) != 0
        {
            return 1;
        }
        aux_print(out_file, &format!("{}{lf}", strbuf_str(strbuf)));
    }
    0
}

// ---------------------------------------------------------------------------
// The original structure for the reversibility part of the AuxInfo

use super::ichi_bns::n_bonds_valence_inp_at;
use super::ichisort::{comp_at_rank, insertions_sort};
use super::ichister::{
    AB_PARITY_EVEN, AB_PARITY_ODD, AB_PARITY_UNDF, AB_PARITY_UNKN, atom_parity_well_def,
};
use super::inpdef::{InpAtom, MAX_NUM_STEREO_ATOM_NEIGH};
use super::util::{get_atomic_mass_from_elnum, lrtrim, needed_unusual_el_valence};

pub const LEN_COORD: usize = 10;
pub const NUM_COORD: usize = 3;
/// `sizeof(MOL_COORD)`
pub const MOL_COORD_LEN: usize = LEN_COORD * NUM_COORD + NUM_COORD - 1;
pub const FLAG_INP_AT_CHIRAL: i32 = 1;
pub const FLAG_INP_AT_NONCHIRAL: i32 = 2;

pub type MolCoord = [u8; MOL_COORD_LEN];

/// the longest prefix of `s` that `strtod` reads (decimal numbers only)
fn strtod(s: &[u8]) -> (f64, usize) {
    let n = s.len();
    let mut i = 0;
    if i < n && (s[i] == b'+' || s[i] == b'-') {
        i += 1;
    }
    let ds = i;
    while i < n && s[i].is_ascii_digit() {
        i += 1;
    }
    let mut digits = i - ds;
    if i < n && s[i] == b'.' {
        let f = i + 1;
        let mut k = f;
        while k < n && s[k].is_ascii_digit() {
            k += 1;
        }
        if digits + (k - f) > 0 {
            digits += k - f;
            i = k;
        }
    }
    if digits == 0 {
        return (0.0, 0);
    }
    if i < n && (s[i] == b'e' || s[i] == b'E') {
        let mut k = i + 1;
        if k < n && (s[k] == b'+' || s[k] == b'-') {
            k += 1;
        }
        let es = k;
        while k < n && s[k].is_ascii_digit() {
            k += 1;
        }
        if k > es {
            i = k;
        }
    }
    let v = std::str::from_utf8(&s[..i])
        .ok()
        .and_then(|t| t.parse::<f64>().ok());
    (v.unwrap_or(0.0), i)
}

/// `strtol(.., 10)`: (value, bytes read)
fn strtol10(s: &[u8]) -> (i32, usize) {
    let mut i = 0;
    while i < s.len() && matches!(s[i], b' ' | 0x09..=0x0d) {
        i += 1;
    }
    let mut neg = false;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let ds = i;
    let mut v: i64 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        v = (v * 10 + i64::from(s[i] - b'0')).min(i64::from(i32::MAX));
        i += 1;
    }
    if i == ds {
        return (0, 0);
    }
    (if neg { -v } else { v } as i32, i)
}

/// `CleanOrigCoord`: removes the insignificant characters of the three
/// coordinates and joins them with `delim`; returns the number of zero ones
pub fn clean_orig_coord(sz_coord: &mut MolCoord, delim: u8) -> i32 {
    const MIN_BOND_LENGTH: f64 = 1.0e-6;
    let mut sz_buf = [0u8; MOL_COORD_LEN];
    let (mut len_buf, mut num_zer) = (0usize, 0);
    for k in (0..NUM_COORD * LEN_COORD).step_by(LEN_COORD) {
        // the value as a zero terminated string
        let raw = &sz_coord[k..k + LEN_COORD];
        let end = raw.iter().position(|&c| c == 0).unwrap_or(LEN_COORD);
        let mut v = raw[..end].to_vec();
        lrtrim(&mut v);
        let mut sz_val = [0u8; LEN_COORD + 1];
        sz_val[..v.len()].copy_from_slice(&v);
        let (coord, q) = strtod(&v);
        let mut len: usize;
        if MIN_BOND_LENGTH > coord.abs() {
            sz_val = [0u8; LEN_COORD + 1];
            sz_val[0] = b'0';
            len = 1;
            num_zer += 1;
        } else {
            len = q;
            let mut last = len;
            if let Some(p) = v
                .iter()
                .position(|&c| matches!(c, b'e' | b'E' | b'd' | b'D'))
            {
                last = p;
                let (e, _) = strtol10(&sz_val[last + 1..]);
                if e != 0 {
                    let t = e.to_string();
                    sz_val[last + 1..last + 1 + t.len()].copy_from_slice(t.as_bytes());
                    sz_val[last + 1 + t.len()] = 0;
                    len = last + 1 + t.len();
                } else {
                    len = last;
                }
            }
            let fst = usize::from(sz_val[0] != b'.' && !sz_val[0].is_ascii_digit());
            let cstr_len = sz_val.iter().position(|&c| c == 0).unwrap_or(LEN_COORD);
            let dec_pnt: i64 = sz_val[..cstr_len]
                .iter()
                .position(|&c| c == b'.')
                .map_or(last as i64, |p| p as i64);
            let last_i = last as i64 - 1;
            let mut i = last_i;
            while dec_pnt < i && sz_val[i as usize] == b'0' {
                i -= 1;
            }
            if i == dec_pnt {
                i -= 1;
            }
            if i < last_i {
                let from = (last_i + 1) as usize;
                let cnt = len as i64 - last as i64;
                if cnt > 0 {
                    let cnt = cnt as usize;
                    sz_val.copy_within(from..from + cnt, (i + 1) as usize);
                }
                len = (len as i64 - (last_i - i)) as usize;
            }
            let mut i = fst;
            while i < len && sz_val[i] == b'0' {
                i += 1;
            }
            if i > fst {
                sz_val.copy_within(i..len, fst);
                len -= i - fst;
            }
        }
        if len_buf != 0 && len_buf < MOL_COORD_LEN {
            sz_buf[len_buf] = delim;
            len_buf += 1;
        }
        sz_buf[len_buf..len_buf + len].copy_from_slice(&sz_val[..len]);
        len_buf += len;
    }
    if len_buf < MOL_COORD_LEN {
        sz_buf[len_buf..].fill(0);
    }
    *sz_coord = sz_buf;
    num_zer
}

/// `WriteOrigCoord`: the coordinates of the atoms from `*i` on, each ended
/// with `;`; returns the text (at most `buf_len - 1` characters)
pub fn write_orig_coord(
    num_inp_atoms: usize,
    sz_mol_coord: &[MolCoord],
    i: &mut usize,
    buf_len: usize,
) -> String {
    let mut buf: Vec<u8> = Vec::new();
    let mut j = *i;
    while j < num_inp_atoms {
        let mut cur = sz_mol_coord[j];
        let num_zer = clean_orig_coord(&mut cur, b',');
        let len = if NUM_COORD as i32 == num_zer {
            0
        } else {
            cur.iter().position(|&c| c == 0).unwrap_or(MOL_COORD_LEN)
        };
        if len + buf.len() + 1 < buf_len {
            buf.extend_from_slice(&cur[..len]);
            buf.push(b';');
            j += 1;
        } else {
            break;
        }
    }
    *i = j;
    String::from_utf8_lossy(&buf).into_owned()
}

/// `WriteOrigAtoms`: the atoms from `*i` on (element, valence, charge,
/// radical, isotope, parity, isotopic H); the first call (`*i == 0`) starts
/// with the number of atoms and the chirality flag
pub fn write_orig_atoms(
    num_inp_atoms: usize,
    at: &[InpAtom],
    i: &mut usize,
    buf_len: usize,
    chiral_flag: i32,
) -> String {
    const SZ_ISO_H: &[u8; 3] = b"hdt";
    let mut buf = String::new();
    if *i == 0 {
        buf = format!(
            "{num_inp_atoms}{}",
            if chiral_flag & FLAG_INP_AT_CHIRAL != 0 {
                "c"
            } else if chiral_flag & FLAG_INP_AT_NONCHIRAL != 0 {
                "n"
            } else {
                ""
            }
        );
    }
    let mut j = *i;
    while j < num_inp_atoms {
        let a = &at[j];
        let mut parity = 0;
        let mut num_trans = 0;
        if a.p_parity != 0 {
            let mut is_ok = true;
            let mut b_self = 0;
            let mut k = 0;
            let mut neigh_order = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
            for n in 0..MAX_NUM_STEREO_ATOM_NEIGH {
                let neigh = a.p_orig_at_num[n].wrapping_sub(1);
                let in_list = a.neighbor[..a.valence.max(0) as usize].contains(&neigh);
                if in_list && at[usize::from(neigh)].orig_at_number == a.p_orig_at_num[n] {
                    neigh_order[k] = a.p_orig_at_num[n];
                    k += 1;
                } else if usize::from(neigh) == j
                    && at[usize::from(neigh)].orig_at_number == a.p_orig_at_num[n]
                {
                    num_trans = n as i32;
                    b_self += 1;
                } else {
                    is_ok = false;
                    break;
                }
            }
            if is_ok && b_self <= 1 && b_self + k == MAX_NUM_STEREO_ATOM_NEIGH {
                num_trans += insertions_sort(&mut neigh_order[..k], comp_at_rank);
                let p = i32::from(a.p_parity);
                if atom_parity_well_def(p) {
                    parity = 2 - (num_trans + p) % 2;
                } else if p == AB_PARITY_UNKN || p == AB_PARITY_UNDF {
                    parity = p;
                }
            }
        }
        let el_end = a
            .elname
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(a.elname.len());
        let mut cur: String = String::from_utf8_lossy(&a.elname[..el_end]).into_owned();
        let len0 = cur.len();
        let (bonds_val, _, _) = n_bonds_valence_inp_at(a);
        let mut val = needed_unusual_el_valence(
            i32::from(a.el_number),
            i32::from(a.charge),
            i32::from(a.radical),
            i32::from(a.chem_bonds_valence),
            bonds_val,
            i32::from(a.num_h),
            i32::from(a.valence),
        );
        if val != 0
            || a.charge != 0
            || a.radical != 0
            || a.iso_atw_diff != 0
            || a.num_iso_h_total() != 0
            || parity != 0
        {
            if val != 0 {
                cur.push_str(&(if val > 0 { val } else { 0 }).to_string());
            }
            val = i32::from(a.charge);
            if val != 0 {
                cur.push(if val > 0 { '+' } else { '-' });
                if val.abs() > 1 {
                    cur.push_str(&val.abs().to_string());
                }
            }
            val = i32::from(a.radical);
            if val != 0 {
                cur.push_str(&format!(".{val}"));
            }
            val = i32::from(a.iso_atw_diff);
            if val != 0 {
                let mw = get_atomic_mass_from_elnum(i32::from(a.el_number));
                val = if val == 1 {
                    mw
                } else if val > 0 {
                    mw + val - 1
                } else {
                    mw + val
                };
                cur.push_str(&format!(
                    "{}i{val}",
                    if cur.len() == len0 { "." } else { "" }
                ));
            }
            if parity != 0 {
                let s = match parity {
                    AB_PARITY_ODD => "o",
                    AB_PARITY_EVEN => "e",
                    AB_PARITY_UNKN => "u",
                    AB_PARITY_UNDF => "?",
                    _ => "",
                };
                cur.push_str(&format!("{}{s}", if cur.len() == len0 { "." } else { "" }));
            }
            if a.num_iso_h_total() != 0 {
                for k in 0..3 {
                    let v = i32::from(a.num_iso_h[k]);
                    if v != 0 {
                        cur.push_str(&format!(
                            "{}{}",
                            if cur.len() == len0 { "." } else { "" },
                            SZ_ISO_H[k] as char
                        ));
                        if v > 1 {
                            cur.push_str(&v.to_string());
                        }
                    }
                }
            }
        }
        if cur.len() + buf.len() < buf_len {
            buf.push_str(&cur);
            j += 1;
        } else {
            break;
        }
        *i = j;
    }
    buf
}

/// `WriteOrigBonds`: the bonds of the atoms from `*i` on, in ascending order
/// of the neighbor numbers, each atom ended with `;`; `err` receives the
/// structure error messages (`pStrErrStruct`)
pub fn write_orig_bonds(
    num_inp_atoms: usize,
    at: &[InpAtom],
    i: &mut usize,
    buf_len: usize,
    err: &mut String,
) -> String {
    use super::ichister::{
        STEREO_DBLE_EITHER, STEREO_SNGL_DOWN, STEREO_SNGL_EITHER, STEREO_SNGL_UP,
        get_opposite_sb_atom,
    };
    use super::inpdef::{MAX_NUM_STEREO_BONDS, sb_parity_1, sb_parity_2};
    use super::util::is_el_a_metal;
    let mut buf = String::new();
    let mut j = *i;
    while j < num_inp_atoms {
        let aj = &at[j];
        let val = aj.valence.max(0) as usize;
        let mut neigh_order: Vec<u16> = (0..val as u16).collect();
        insertions_sort(&mut neigh_order, |a, b| {
            i32::from(aj.neighbor[usize::from(*a)]) - i32::from(aj.neighbor[usize::from(*b)])
        });
        let mut cur = String::new();
        for &kk in &neigh_order {
            let k = usize::from(kk);
            let j2 = usize::from(aj.neighbor[k]);
            let (mut bond_parity, mut bond_parity_nm) = (0, 0);
            if j2 >= j {
                continue;
            }
            let bond_stereo = i32::from(aj.bond_stereo[k]);
            let bond_char = match i32::from(aj.bond_type[k]) {
                1 => match bond_stereo {
                    x if x == STEREO_SNGL_UP => 'p',
                    x if x == -STEREO_SNGL_UP => 'P',
                    x if x == STEREO_SNGL_DOWN => 'n',
                    x if x == -STEREO_SNGL_DOWN => 'N',
                    x if x == STEREO_SNGL_EITHER => 'v',
                    x if x == -STEREO_SNGL_EITHER => 'V',
                    _ => 's',
                },
                2 => {
                    if bond_stereo == i32::from(STEREO_DBLE_EITHER)
                        || bond_stereo == -i32::from(STEREO_DBLE_EITHER)
                    {
                        'w'
                    } else {
                        'd'
                    }
                }
                3 => 't',
                4 => 'a',
                _ => 's',
            };
            let a2 = &at[j2];
            let k2 = a2.neighbor[..a2.valence.max(0) as usize]
                .iter()
                .position(|&x| usize::from(x) == j)
                .unwrap_or(0);
            let (mut chain_len, mut chain_len2) = (0, 0);
            let (mut pnxt_atom, mut pinxt2cur, mut pinxt_sb_parity_ord) = (0usize, 0usize, 0usize);
            let (mut pnxt_atom2, mut pinxt2cur2, mut pinxt_sb_parity_ord2) =
                (0usize, 0usize, 0usize);
            if aj.sb_parity[0] != 0 {
                let mut m1 = 0;
                while m1 < MAX_NUM_STEREO_BONDS && aj.sb_parity[m1] != 0 {
                    if k == aj.sb_ord[m1] as usize {
                        (chain_len, pnxt_atom, pinxt2cur, pinxt_sb_parity_ord) =
                            get_opposite_sb_atom(at, j, k);
                        break;
                    }
                    m1 += 1;
                }
            }
            if a2.sb_parity[0] != 0 {
                let mut m2 = 0;
                while m2 < MAX_NUM_STEREO_BONDS && a2.sb_parity[m2] != 0 {
                    if k2 == a2.sb_ord[m2] as usize {
                        (chain_len2, pnxt_atom2, pinxt2cur2, pinxt_sb_parity_ord2) =
                            get_opposite_sb_atom(at, j2, k2);
                        break;
                    }
                    m2 += 1;
                }
            }
            let (mut pcur_atom, mut picur2nxt, mut picur_sb_parity_ord) = (0usize, 0usize, -1i32);
            if (chain_len == 1 && chain_len2 == 1) || (chain_len > 1 && j > pnxt_atom) {
                pcur_atom = j;
                picur2nxt = k;
                let mut m = 0;
                while m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m] != 0 {
                    if at[pcur_atom].sb_ord[m] as usize == k {
                        picur_sb_parity_ord = m as i32;
                        break;
                    }
                    m += 1;
                }
            } else if chain_len2 > 1 && j2 > pnxt_atom2 {
                pcur_atom = j2;
                picur2nxt = k2;
                pnxt_atom = pnxt_atom2;
                pinxt2cur = pinxt2cur2;
                pinxt_sb_parity_ord = pinxt_sb_parity_ord2;
                let mut m = 0;
                while m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m] != 0 {
                    if at[pcur_atom].sb_ord[m] as usize == k2 {
                        picur_sb_parity_ord = m as i32;
                    }
                    m += 1;
                }
                chain_len = chain_len2;
            } else {
                chain_len = 0;
            }
            if chain_len != 0 {
                let (mut p1, mut p1_nm) = (0, 0);
                if picur_sb_parity_ord >= 0 {
                    let sp = at[pcur_atom].sb_parity[picur_sb_parity_ord as usize];
                    p1 = sb_parity_1(sp);
                    p1_nm = sb_parity_2(sp);
                }
                let p2 = sb_parity_1(at[pnxt_atom].sb_parity[pinxt_sb_parity_ord]);
                let p2_nm = sb_parity_2(at[pnxt_atom].sb_parity[pinxt_sb_parity_ord]);
                let well_def = atom_parity_well_def(p1) && atom_parity_well_def(p2);
                let has_metal = atom_parity_well_def(p1_nm) && atom_parity_well_def(p2_nm);
                let (mut switched1, mut switched2) = (0, 0);
                if well_def || has_metal {
                    let lowest = |atom: usize, skip: usize| -> usize {
                        let mut n1 = num_inp_atoms;
                        for kc in 0..at[atom].valence.max(0) as usize {
                            if kc == skip {
                                continue;
                            }
                            let neigh = usize::from(at[atom].neighbor[kc]);
                            if has_metal && is_el_a_metal(i32::from(at[neigh].el_number)) {
                                continue;
                            }
                            n1 = n1.min(neigh);
                        }
                        n1
                    };
                    let neigh1 = lowest(pcur_atom, picur2nxt);
                    if neigh1 < num_inp_atoms {
                        let sn = at[pcur_atom].sn_ord[picur_sb_parity_ord as usize] as usize;
                        switched1 = i32::from(neigh1 != usize::from(at[pcur_atom].neighbor[sn]));
                    } else {
                        super::ichierr::add_error_message(
                            err,
                            "Cannot find 0D stereobond neighbor",
                        );
                    }
                    let neigh2 = lowest(pnxt_atom, pinxt2cur);
                    if neigh2 < num_inp_atoms {
                        let sn = at[pnxt_atom].sn_ord[pinxt_sb_parity_ord] as usize;
                        switched2 = i32::from(neigh2 != usize::from(at[pnxt_atom].neighbor[sn]));
                    } else {
                        super::ichierr::add_error_message(
                            err,
                            "Cannot find 0D stereobond neighbor",
                        );
                    }
                    if neigh1 < num_inp_atoms && neigh2 < num_inp_atoms {
                        if atom_parity_well_def(p1) && atom_parity_well_def(p2) {
                            bond_parity = 2 - (p1 + p2 + switched1 + switched2) % 2;
                        } else {
                            bond_parity = p1.min(p2);
                        }
                        if has_metal {
                            bond_parity_nm = 2 - (p1_nm + p2_nm + switched1 + switched2) % 2;
                        } else if p1_nm != 0 && p2_nm != 0 {
                            bond_parity_nm = p1_nm.min(p2_nm);
                        }
                    }
                } else {
                    if p1 != 0 && p2 != 0 {
                        bond_parity = p1.min(p2);
                    }
                    if p1_nm != 0 && p2_nm != 0 {
                        bond_parity_nm = p1_nm.min(p2_nm);
                    }
                    if bond_parity_nm != 0 && bond_parity == 0 {
                        bond_parity = AB_PARITY_UNDF;
                    }
                }
            }
            let ps = |p: i32| match p {
                AB_PARITY_ODD => "-",
                AB_PARITY_EVEN => "+",
                AB_PARITY_UNKN => "u",
                AB_PARITY_UNDF => "?",
                _ => "",
            };
            cur.push_str(&format!(
                "{bond_char}{}{}{}",
                ps(bond_parity),
                ps(bond_parity_nm),
                j2 + 1
            ));
        }
        if cur.len() + buf.len() + 2 < buf_len {
            buf.push_str(&cur);
            buf.push(';');
            j += 1;
        } else {
            break;
        }
    }
    *i = if num_inp_atoms > 0 { j } else { 0 };
    buf
}
