//! Option parsing of InChI (`INCHI_BASE/src/ichiparm.c`, the part shared by
//! the command line and the library: `set_common_options_by_parg`).
//!
//! The C function sets about forty variables through pointers; they are the
//! fields of [`CommonOptions`] here, and the fields of `INPUT_PARMS` it
//! touches are in [`InputParms`] (more are added with the code that needs
//! them). The library build (`USE_STDINCHI_API` undefined) is translated.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::util::{inchi_memicmp, inchi_stricmp, lrtrim, mystrncpy, strtol10};

// ichi.h
pub const REQ_MODE_BASIC: u64 = 0x000001;
pub const REQ_MODE_TAUT: u64 = 0x000002;
pub const REQ_MODE_ISO: u64 = 0x000004;
pub const REQ_MODE_NON_ISO: u64 = 0x000008;
pub const REQ_MODE_STEREO: u64 = 0x000010;
pub const REQ_MODE_ISO_STEREO: u64 = 0x000020;
pub const REQ_MODE_NOEQ_STEREO: u64 = 0x000040;
pub const REQ_MODE_REDNDNT_STEREO: u64 = 0x000080;
pub const REQ_MODE_NO_ALT_SBONDS: u64 = 0x000100;
pub const REQ_MODE_RELATIVE_STEREO: u64 = 0x000200;
pub const REQ_MODE_RACEMIC_STEREO: u64 = 0x000400;
pub const REQ_MODE_SC_IGN_ALL_UU: u64 = 0x000800;
pub const REQ_MODE_SB_IGN_ALL_UU: u64 = 0x001000;
pub const REQ_MODE_CHIR_FLG_STEREO: u64 = 0x002000;
pub const REQ_MODE_DIFF_UU_STEREO: u64 = 0x004000;
pub const REQ_MODE_MIN_SB_RING_MASK: u64 = 0x0F0000;
pub const REQ_MODE_MIN_SB_RING_SHFT: u32 = 16;
pub const REQ_MODE_DEFAULT: u64 =
    REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO;

pub const INCHI_OUT_NO_AUX_INFO: i32 = 0x0001;
pub const INCHI_OUT_SHORT_AUX_INFO: i32 = 0x0002;
pub const INCHI_OUT_ONLY_AUX_INFO: i32 = 0x0004;
pub const INCHI_OUT_EMBED_REC: i32 = 0x0008;
pub const INCHI_OUT_SDFILE_ONLY: i32 = 0x0010;
pub const INCHI_OUT_XML: i32 = 0x0020;
pub const INCHI_OUT_PLAIN_TEXT: i32 = 0x0040;
pub const INCHI_OUT_PLAIN_TEXT_COMMENTS: i32 = 0x0080;
pub const INCHI_OUT_XML_TEXT_COMMENTS: i32 = 0x0100;
pub const INCHI_OUT_WINCHI_WINDOW: i32 = 0x0200;
pub const INCHI_OUT_TABBED_OUTPUT: i32 = 0x0400;
pub const INCHI_OUT_SDFILE_ATOMS_DT: i32 = 0x0800;
pub const INCHI_OUT_SDFILE_SPLIT: i32 = 0x1000;
pub const INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG: i32 = 0x2000;
pub const INCHI_OUT_STDINCHI: i32 = 0x4000;
pub const INCHI_OUT_SAVEOPT: i32 = 0x8000;

// second set of output options
pub const INCHI_OUT_INCHI_GEN_ERROR: i32 = 0x0001;
pub const INCHI_OUT_MISMATCH_AS_ERROR: i32 = 0x0002;
pub const INCHI_OUT_REQ_LARGE_MOL: i32 = 0x0004;
pub const INCHI_OUT_REQ_POLYMERS: i32 = 0x0008;

// inpdef.h
pub const FLAG_SET_INP_AT_CHIRAL: i32 = 4;
pub const FLAG_SET_INP_AT_NONCHIRAL: i32 = 8;

// mode.h
pub const POLYMERS_NO: i32 = 0;
pub const POLYMERS_MODERN: i32 = 1;
pub const POLYMERS_LEGACY: i32 = 2;
pub const POLYMERS_LEGACY_PLUS: i32 = 3;

// ichidrp.h: `tagFrameShiftScheme`
pub const FSS_STARS_CYCLED: i32 = 0;
pub const FSS_NONE: i32 = 1;
pub const FSS_STARS_CYCLED_SORTED: i32 = 2;
pub const FSS_STARS_OPENED: i32 = 3;
pub const FSS_STARS_ENDS_OPENED: i32 = 4;

// ichidrp.h: `tagInputType`
pub const INPUT_NONE: i32 = 0;
pub const INPUT_MOLFILE: i32 = 1;
pub const INPUT_SDFILE: i32 = 2;
pub const INPUT_INCHI_XML: i32 = 3;
pub const INPUT_INCHI_PLAIN: i32 = 4;

/// The `INPUT_PARMS` fields set by [`set_common_options_by_parg`]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InputParms {
    pub n_input_type: i32,
    pub first_struct_number: i64,
    pub last_struct_number: i64,
    pub b_save_all_good_structs_as_problem: i32,
    pub b_filter_ss: i32,
    pub b_fix_non_uniform_draw: i32,
}

/// The variables `set_common_options_by_parg` sets through pointers
/// (`INCHI_MODE *pbVer1DefaultMode`, `int *pnMode`, ...)
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CommonOptions {
    pub b_ver1_default_mode: u64,
    pub n_mode: u64,
    pub b_inchi_output_options: i32,
    pub b_inchi_output_options2: i32,
    pub b_std_format: i32,
    pub b_hash_key: i32,
    pub b_hash_xtra1: i32,
    pub b_hash_xtra2: i32,
    pub b_fix_sp3_bug: i32,
    pub b_fix_fb2: i32,
    pub b_add_phosphine_stereo: i32,
    pub b_add_arsine_stereo: i32,
    pub b_no_struct_labels: i32,
    pub b_pointed_edge_stereo: i32,
    pub b_do_not_add_h: i32,
    pub b_forced_chiral_flag: i32,
    pub b_reconnect_coord: i32,
    pub b_keto_enol_taut: i32,
    pub b15_taut_non_ring: i32,
    pub b_pt_06_00_taut: i32,
    pub b_pt_13_00_taut: i32,
    pub b_pt_16_00_taut: i32,
    pub b_pt_18_00_taut: i32,
    pub b_pt_22_00_taut: i32,
    pub b_pt_39_00_taut: i32,
    pub b_loose_tsa_check: i32,
    pub b_large_molecules: i32,
    pub b_polymers: i32,
    pub b_fold_polymer_sru: i32,
    pub b_frame_shift_scheme: i32,
    pub b_stereo_at_zz: i32,
    pub b_np_zz: i32,
    pub b_no_warnings: i32,
    pub b_merge_hash: i32,
    pub b_hide_inchi: i32,
}

/// `set_common_options_by_parg`: returns whether `p_arg` was recognized and
/// treated. (As in C, the `PT_xx_00` options set their flag but report the
/// argument as not recognized.)
pub fn set_common_options_by_parg(
    p_arg: &str,
    developer_options: bool,
    ip: &mut InputParms,
    o: &mut CommonOptions,
) -> bool {
    let is = |name: &str| inchi_stricmp(p_arg, name) == 0;
    let prefix = |name: &str| {
        p_arg.len() >= name.len()
            && inchi_memicmp(p_arg.as_bytes(), name.as_bytes(), name.len()) == 0
    };
    let after = |n: usize| strtol10(&p_arg.as_bytes()[n..]).0;
    let mut got = false;

    // Input
    if is("INPAUX") {
        if ip.n_input_type == INPUT_NONE {
            ip.n_input_type = INPUT_INCHI_PLAIN;
        }
        got = true;
    } else if prefix("START:") {
        ip.first_struct_number = after(6);
        got = true;
    } else if prefix("END:") {
        ip.last_struct_number = after(4);
        got = true;
    } else if prefix("RECORD:") {
        let num = after(7);
        ip.first_struct_number = num;
        ip.last_struct_number = num;
        got = true;
    } else if is("NOLABELS") {
        o.b_no_struct_labels = 1;
        got = true;
    } else if is("SAVEOPT") {
        o.b_inchi_output_options |= INCHI_OUT_SAVEOPT;
        got = true;
    }
    // Generation
    else if is("AUXNONE") {
        o.b_inchi_output_options |= INCHI_OUT_NO_AUX_INFO;
        o.b_inchi_output_options &= !INCHI_OUT_SHORT_AUX_INFO;
        got = true;
    } else if is("MISMATCHISERROR") {
        o.b_inchi_output_options2 |= INCHI_OUT_MISMATCH_AS_ERROR;
        got = true;
    } else if is("OUTERRINCHI") {
        o.b_inchi_output_options2 |= INCHI_OUT_INCHI_GEN_ERROR;
        got = true;
    }
    // InChIKey/InChI hash
    else if is("Key") {
        o.b_hash_key = 1;
        got = true;
    } else if is("XHash1") {
        o.b_hash_xtra1 = 1;
        got = true;
    } else if is("XHash2") {
        o.b_hash_xtra2 = 1;
        got = true;
    }
    // Options that do not turn off the standard flag
    else if is("SNON") {
        o.b_ver1_default_mode &= !REQ_MODE_STEREO;
        o.n_mode &=
            !(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
        got = true;
    } else if is("NEWPSOFF") {
        o.b_pointed_edge_stereo = 0;
        got = true;
    } else if is("DONOTADDH") {
        o.b_do_not_add_h = 1;
        got = true;
    } else if is("LooseTSACheck") {
        o.b_loose_tsa_check = 1;
        got = true;
    }
    // Options that turn off the standard flag
    else if is("SREL") {
        o.n_mode &= !(REQ_MODE_RACEMIC_STEREO | REQ_MODE_CHIR_FLG_STEREO);
        o.n_mode |= REQ_MODE_RELATIVE_STEREO | REQ_MODE_STEREO;
        o.b_std_format = 0;
        got = true;
    } else if is("SRAC") {
        o.n_mode &= !(REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
        o.n_mode |= REQ_MODE_RACEMIC_STEREO | REQ_MODE_STEREO;
        o.b_std_format = 0;
        got = true;
    } else if is("SUCF") {
        o.n_mode &= !(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
        o.n_mode |= REQ_MODE_CHIR_FLG_STEREO | REQ_MODE_STEREO;
        o.b_std_format = 0;
        got = true;
    } else if is("ChiralFlagON") {
        o.b_forced_chiral_flag &= !FLAG_SET_INP_AT_NONCHIRAL;
        o.b_forced_chiral_flag |= FLAG_SET_INP_AT_CHIRAL;
        got = true;
    } else if is("ChiralFlagOFF") {
        o.b_forced_chiral_flag &= !FLAG_SET_INP_AT_CHIRAL;
        o.b_forced_chiral_flag |= FLAG_SET_INP_AT_NONCHIRAL;
        got = true;
    } else if is("SUU") {
        o.b_ver1_default_mode &= !(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
        o.b_std_format = 0;
        got = true;
    } else if is("SLUUD") {
        o.b_ver1_default_mode |= REQ_MODE_DIFF_UU_STEREO;
        o.b_std_format = 0;
        got = true;
    } else if is("FIXEDH") {
        o.b_ver1_default_mode |= REQ_MODE_BASIC;
        o.b_std_format = 0;
        got = true;
    } else if is("RECMET") {
        o.b_reconnect_coord = 1;
        o.b_std_format = 0;
        got = true;
    } else if is("KET") {
        o.b_keto_enol_taut = 1;
        o.b_std_format = 0;
        got = true;
    } else if is("15T") {
        o.b15_taut_non_ring = 1;
        o.b_std_format = 0;
        got = true;
    } else if is("PT_22_00") {
        o.b_pt_22_00_taut = 1;
        o.b_std_format = 0;
    } else if is("PT_16_00") {
        o.b_pt_16_00_taut = 1;
        o.b_std_format = 0;
    } else if is("PT_06_00") {
        o.b_pt_06_00_taut = 1;
        o.b_std_format = 0;
    } else if is("PT_39_00") {
        o.b_pt_39_00_taut = 1;
        o.b_std_format = 0;
    } else if is("PT_13_00") {
        o.b_pt_13_00_taut = 1;
        o.b_std_format = 0;
    } else if is("PT_18_00") {
        o.b_pt_18_00_taut = 1;
        o.b_std_format = 0;
    } else if is("LargeMolecules") {
        o.b_large_molecules = 1;
        got = true;
    } else if is("Polymers") {
        o.b_polymers = POLYMERS_MODERN;
        got = true;
    } else if is("Polymers105") {
        o.b_polymers = POLYMERS_LEGACY;
        got = true;
    } else if is("NPZz") {
        o.b_np_zz = 1;
        got = true;
    } else if is("NoWarnings") {
        o.b_no_warnings = 1;
        got = true;
    } else if is("MergeHash") {
        o.b_merge_hash = 1;
        got = true;
    } else if is("NoInChI") || is("HideInChI") {
        o.b_hide_inchi = 1;
        got = true;
    } else if is("FoldCRU") || is("FoldSRU") {
        o.b_fold_polymer_sru = 1;
        got = true;
    } else if prefix("FrameShift:") {
        let mut wrd = vec![0u8; 256];
        mystrncpy(&mut wrd, &p_arg.as_bytes()[11..], 256);
        wrd.truncate(wrd.iter().position(|&c| c == 0).unwrap_or(256));
        let k = lrtrim(&mut wrd);
        if k != 0 {
            let w = String::from_utf8_lossy(&wrd).to_string();
            if inchi_stricmp(&w, "None") == 0 {
                o.b_frame_shift_scheme = FSS_NONE;
            } else if inchi_stricmp(&w, "Cyclize") == 0 {
                o.b_frame_shift_scheme = FSS_STARS_CYCLED;
            } else if inchi_stricmp(&w, "MoveStars") == 0 {
                o.b_frame_shift_scheme = FSS_STARS_OPENED;
            } else if inchi_stricmp(&w, "MoveBrackets") == 0 {
                o.b_frame_shift_scheme = FSS_STARS_ENDS_OPENED;
            }
        } else {
            o.b_frame_shift_scheme = FSS_STARS_CYCLED;
        }
        got = true;
    } else if is("NoFrameShift") {
        o.b_frame_shift_scheme = FSS_NONE;
        got = true;
    } else if is("NoEdits") {
        o.b_fold_polymer_sru = 0;
        o.b_frame_shift_scheme = FSS_NONE;
        got = true;
    } else if is("SATZZ") {
        o.b_stereo_at_zz = 1;
        got = true;
    }

    if !got && developer_options {
        if is("PGO") {
            ip.b_save_all_good_structs_as_problem = 1;
            got = true;
        } else if is("FilterSS") {
            ip.b_filter_ss = 1;
            got = true;
        } else if is("InvFilterSS") {
            ip.b_filter_ss = -1;
            got = true;
        }
        // Options that turn off the standard flag (a separate `if` in C, so
        // it also runs after the ones above)
        if is("FNUDOFF") {
            ip.b_fix_non_uniform_draw = 0;
            o.b_std_format = 0;
            got = true;
        } else if is("FixSp3bugOFF") || is("FBOFF") {
            o.b_fix_sp3_bug = 0;
            o.b_std_format = 0;
            got = true;
        } else if is("FB2OFF") {
            o.b_fix_fb2 = 0;
            o.b_std_format = 0;
            got = true;
        } else if is("SPXYZOFF") {
            o.b_add_phosphine_stereo = 0;
            o.b_std_format = 0;
            got = true;
        } else if is("SASXYZOFF") {
            o.b_add_arsine_stereo = 0;
            o.b_std_format = 0;
            got = true;
        } else if is("Polymers105+") {
            o.b_polymers = POLYMERS_LEGACY_PLUS;
            o.b_std_format = 0;
            got = true;
        }
    }
    got
}
