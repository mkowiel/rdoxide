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

/// `MAX_NUM_PATHS`
pub const MAX_NUM_PATHS: usize = 4;
/// `MAX_SDF_HEADER`
pub const MAX_SDF_HEADER: usize = 64;

// ichidrp.h: `tagInputType` (rest)
pub const INPUT_CMLFILE: i32 = 5;
pub const INPUT_INCHI: i32 = 6;

// ichidrp.h: `tagInChIHashCalc`
pub const INCHIHASH_NONE: i32 = 0;
pub const INCHIHASH_KEY: i32 = 1;
pub const INCHIHASH_KEY_XTRA1: i32 = 2;
pub const INCHIHASH_KEY_XTRA2: i32 = 3;
pub const INCHIHASH_KEY_XTRA1_XTRA2: i32 = 4;

// ichi.h: component comparison
pub const CMP_COMPONENTS: i32 = 0x0001;
pub const CMP_COMPONENTS_NONISO: i32 = 0x0002;
pub const CMP_COMPONENTS_NONTAUT: i32 = 0x0004;

// mode.h: `READ_INCHI_*`
pub const READ_INCHI_OUTPUT_INCHI: i32 = 0x00000001;
pub const READ_INCHI_SPLIT_OUTPUT: i32 = 0x00000002;
pub const READ_INCHI_KEEP_BALANCE_P: i32 = 0x00000004;
pub const READ_INCHI_TO_STRUCTURE: i32 = 0x00000008;

// mode.h: `TG_FLAG_*` (`INPUT_PARMS::bTautFlags`)
pub const TG_FLAG_TEST_TAUT__ATOMS: u64 = 0x00000001;
pub const TG_FLAG_DISCONNECT_SALTS: u64 = 0x00000002;
pub const TG_FLAG_TEST_TAUT__SALTS: u64 = 0x00000004;
pub const TG_FLAG_MOVE_POS_CHARGES: u64 = 0x00000008;
pub const TG_FLAG_TEST_TAUT2_SALTS: u64 = 0x00000010;
pub const TG_FLAG_ALLOW_NO_NEGTV_O: u64 = 0x00000020;
pub const TG_FLAG_MERGE_TAUT_SALTS: u64 = 0x00000040;
pub const TG_FLAG_DISCONNECT_COORD: u64 = 0x00000080;
pub const TG_FLAG_RECONNECT_COORD: u64 = 0x00000100;
pub const TG_FLAG_CHECK_VALENCE_COORD: u64 = 0x00000200;
pub const TG_FLAG_MOVE_HPLUS2NEUTR: u64 = 0x00000400;
pub const TG_FLAG_VARIABLE_PROTONS: u64 = 0x00000800;
pub const TG_FLAG_HARD_ADD_REM_PROTONS: u64 = 0x00001000;
pub const TG_FLAG_POINTED_EDGE_STEREO: u64 = 0x00002000;
pub const TG_FLAG_FIX_ADJ_RADICALS: u64 = 0x00004000;
pub const TG_FLAG_PHOSPHINE_STEREO: u64 = 0x00008000;
pub const TG_FLAG_ARSINE_STEREO: u64 = 0x00010000;
pub const TG_FLAG_H_ALREADY_REMOVED: u64 = 0x00020000;
pub const TG_FLAG_FIX_SP3_BUG: u64 = 0x00040000;
pub const TG_FLAG_KETO_ENOL_TAUT: u64 = 0x00080000;
pub const TG_FLAG_1_5_TAUT: u64 = 0x00100000;
pub const TG_FLAG_FIX_ISO_FIXEDH_BUG: u64 = 0x00200000;
pub const TG_FLAG_FIX_TERM_H_CHRG_BUG: u64 = 0x00400000;
pub const TG_FLAG_PT_22_00: u64 = 0x00800000;
pub const TG_FLAG_PT_16_00: u64 = 0x01000000;
pub const TG_FLAG_PT_06_00: u64 = 0x02000000;
pub const TG_FLAG_PT_39_00: u64 = 0x04000000;
pub const TG_FLAG_PT_13_00: u64 = 0x08000000;
pub const TG_FLAG_PT_18_00: u64 = 0x10000000;

/// `INPUT_PARMS` of the library build (`ichidrp.h`; the drawing parameters
/// of the GUI/stand-alone builds are not part of it). `pSdfLabel` points
/// into `szSdfDataHeader` and `pSdfValue` into a caller's buffer, so they are
/// flags here.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InputParms {
    pub sz_sdf_data_header: String,
    pub p_sdf_label: bool,
    pub p_sdf_value: bool,
    pub l_sdf_id: u64,
    pub l_molfile_number: i64,
    pub path: [Option<String>; MAX_NUM_PATHS],
    pub num_paths: usize,
    pub first_struct_number: i64,
    pub last_struct_number: i64,
    pub n_input_type: i32,
    pub n_mode: u64,
    pub b_abc_numbers: i32,
    pub b_inchi_output_options: i32,
    pub b_inchi_output_options2: i32,
    pub b_ct_predecessors: i32,
    pub b_display_each_component_inchi: i32,
    pub msec_max_time: i64,
    pub msec_left_time: i64,
    pub ul_displ_time: i64,
    pub b_display: i32,
    pub b_display_if_restore_warnings: i32,
    pub b_merge_all_input_structures: i32,
    pub b_save_warning_structs_as_problem: i32,
    pub b_save_all_good_structs_as_problem: i32,
    pub b_get_sdfile_id: i32,
    pub b_get_molfile_number: i32,
    pub b_compare_components: i32,
    pub b_display_composite_results: i32,
    pub b_do_not_add_h: i32,
    pub b_no_struct_labels: i32,
    pub b_fix_non_uniform_draw: i32,
    pub b_calc_inchi_hash: i32,
    pub b_chiral_flag: i32,
    pub b_allow_empty_structure: i32,
    pub b_large_molecules: i32,
    pub b_loose_tsa_check: i32,
    pub b_polymers: i32,
    pub b_fold_polymer_sru: i32,
    pub b_frame_shift_scheme: i32,
    pub b_np_zz: i32,
    pub b_stereo_at_zz: i32,
    pub b_merge_hash: i32,
    pub b_no_warnings: i32,
    pub b_hide_inchi: i32,
    pub b_taut_flags: u64,
    pub b_taut_flags_done: u64,
    pub b_read_inchi_options: i32,
    pub b_underivatize: i32,
    pub b_ring2_chain: i32,
    pub b_ignore_unchanged: i32,
    pub b_filter_ss: i32,
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

/// `INCHI_MAX_NUM_ARG` (`inchi_dll.h`)
pub const INCHI_MAX_NUM_ARG: usize = 32;

/// `parse_options_string` (`inchi_dll.c`): splits the option string into
/// arguments with the quoting rules of the Windows command line (a pair of
/// backslashes is one backslash, a backslash before a quote makes it
/// literal, `""` inside quotes is a quote). The C function edits the string
/// in place; here the arguments are returned. `argv[0]` is the empty program
/// name, as in C, and at most `maxargs - 2` arguments follow it.
pub fn parse_options_string(cmd: &[u8], maxargs: usize) -> Vec<Vec<u8>> {
    // the C string with its terminating zero; reads past the end give 0
    let mut buf: Vec<u8> = cmd.iter().copied().take_while(|&c| c != 0).collect();
    buf.push(0);
    let at = |buf: &Vec<u8>, k: usize| buf.get(k).copied().unwrap_or(0);
    let mut argv: Vec<Vec<u8>> = vec![Vec::new()];
    let mut p = 0usize;
    let mut inside_quotes = false;
    while argv.len() + 1 < maxargs {
        while matches!(at(&buf, p), b' ' | b'\t') {
            p += 1;
        }
        if at(&buf, p) == 0 {
            break;
        }
        let start = p;
        let mut cur = p; // pArgCurChar
        loop {
            let mut copy_char = true;
            let mut backslashes = 0usize;
            while at(&buf, p) == b'\\' {
                p += 1;
                backslashes += 1;
            }
            if at(&buf, p) == b'"' {
                if backslashes % 2 == 0 {
                    if inside_quotes {
                        if at(&buf, p + 1) == b'"' {
                            p += 1;
                        } else {
                            copy_char = false;
                        }
                    } else {
                        copy_char = false;
                    }
                    inside_quotes = !inside_quotes;
                }
                backslashes /= 2;
            }
            for _ in 0..backslashes {
                buf[cur] = b'\\';
                cur += 1;
            }
            if at(&buf, p) == 0 {
                break;
            }
            if !inside_quotes && matches!(at(&buf, p), b' ' | b'\t') {
                p += 1;
                break;
            }
            if copy_char {
                buf[cur] = buf[p];
                cur += 1;
            }
            p += 1;
        }
        let arg = buf[start..cur].to_vec();
        if cur < buf.len() {
            buf[cur] = 0;
        }
        argv.push(arg);
    }
    argv
}

use super::ichi_io::{InchiIostream, inchi_ios_eprint};

/// `INCHI_OPTION_PREFX` of the non-Windows build
pub const INCHI_OPTION_PREFX: u8 = b'-';
const LONG_MAX: i64 = i64::MAX;

/// `VER103_DEFAULT_MODE`
const VER103_DEFAULT_MODE: u64 = REQ_MODE_TAUT
    | REQ_MODE_ISO
    | REQ_MODE_STEREO
    | REQ_MODE_SB_IGN_ALL_UU
    | REQ_MODE_SC_IGN_ALL_UU;

/// C's `strtol(s, &end, 10)` with the overflow indication (`ERANGE`)
fn strtol_range(s: &[u8]) -> (i64, usize, bool) {
    let (v, used) = strtol10(s);
    let digits = s[..used].iter().filter(|c| c.is_ascii_digit()).count();
    let overflow = used > 0 && (v == i64::MAX || digits > 18 && v == 0);
    (v, used, overflow)
}

/// C's `strtod(s, &end)` for decimal numbers and `inf`/`nan`: the value, the
/// number of bytes consumed and the range error flag. (Hexadecimal floating
/// point is not recognized.)
fn strtod_prefix(s: &[u8]) -> (f64, usize, bool) {
    let mut i = 0;
    while s
        .get(i)
        .is_some_and(|c| c.is_ascii_whitespace() || *c == 0x0b)
    {
        i += 1;
    }
    let start = i;
    if matches!(s.get(i), Some(b'+' | b'-')) {
        i += 1;
    }
    let lower: Vec<u8> = s[i..].iter().take(8).map(u8::to_ascii_lowercase).collect();
    let neg = s.get(start) == Some(&b'-');
    if lower.starts_with(b"infinity") {
        return (
            if neg {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            },
            i + 8,
            false,
        );
    }
    if lower.starts_with(b"inf") {
        return (
            if neg {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            },
            i + 3,
            false,
        );
    }
    if lower.starts_with(b"nan") {
        return (f64::NAN, i + 3, false);
    }
    let digits_start = i;
    while s.get(i).is_some_and(u8::is_ascii_digit) {
        i += 1;
    }
    let mut n_digits = i - digits_start;
    if s.get(i) == Some(&b'.') {
        let frac = i + 1;
        let mut j = frac;
        while s.get(j).is_some_and(u8::is_ascii_digit) {
            j += 1;
        }
        n_digits += j - frac;
        if n_digits > 0 {
            i = j;
        }
    }
    if n_digits == 0 {
        return (0.0, 0, false);
    }
    if matches!(s.get(i), Some(b'e' | b'E')) {
        let mut j = i + 1;
        if matches!(s.get(j), Some(b'+' | b'-')) {
            j += 1;
        }
        if s.get(j).is_some_and(u8::is_ascii_digit) {
            while s.get(j).is_some_and(u8::is_ascii_digit) {
                j += 1;
            }
            i = j;
        }
    }
    let text = std::str::from_utf8(&s[start..i]).unwrap_or("0");
    let v: f64 = text.parse().unwrap_or(0.0);
    let nonzero_mantissa = s[digits_start..i]
        .iter()
        .take_while(|c| !matches!(c, b'e' | b'E'))
        .any(|c| (b'1'..=b'9').contains(c));
    let erange = v.is_infinite() || (v == 0.0 && nonzero_mantissa);
    (v, i, erange)
}

/// `ReadCommandLineParms` as built into the library (`TARGET_API_LIB`: no
/// GUI options, no developer options, no display; messages go to `log_file`).
/// `argv[0]` is the unused program name. Returns 0, or -1 when the options
/// are contradictory (a message is written to the log).
pub fn read_command_line_parms(
    argv: &[&str],
    ip: &mut InputParms,
    b_release_version: bool,
    log_file: &mut InchiIostream,
) -> i32 {
    let developer_options = false;
    let ext = [".mol", ".txt", ".log", ".prb", ""];

    let mut n_mode: u64 = 0;
    let _n_release_mode =
        n_mode | (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO);
    let mut b_ver1_default_mode: u64 = VER103_DEFAULT_MODE;
    let mut b_name_suffix;
    let sz_name_suffix = String::new(); // only the GUI option O: sets it
    let b_disconnect_salts = true; // only developer options change it
    let b_tg_flag_variable_protons = true; // only developer options change it
    let b_tg_flag_hard_add_ren_protons = true; // only developer options change it
    let b_disconnect_coord = true; // only developer options change it
    let b_disconnect_coord_chk_val = false; // only developer options change it
    let b_move_positive_charges = true; // only developer options change it
    let b_acid_tautomerism = 2; // only developer options change it
    let b_uncharged_acid_taut = false; // only developer options change it
    let b_merge_salt_t_groups = true; // only developer options change it
    let n_min_db_rin_size = 8u64; // MIN_SB_RING_SIZE (the developer option RSB: changes it)
    let mut b_output_style = INCHI_OUT_PLAIN_TEXT;
    let b_display = 0;
    let mut b_output_molfile_only = false;
    let mut b_output_molfile_dt = false;
    let b_output_molfile_split = false;
    let mut timeout_set_warning = false;
    let mut timeout_set_error = false;
    // C's `errno`, which the timeout code tests but never resets: once a
    // conversion has overflowed it stays ERANGE
    let mut errno_erange = false;
    let b_display_if_restore_warnings = 0;

    // The variables of `set_common_options_by_parg`
    let mut o = CommonOptions {
        n_mode: 0,
        b_ver1_default_mode: 0,
        b_inchi_output_options: INCHI_OUT_EMBED_REC,
        b_inchi_output_options2: 0,
        b_std_format: 1,
        b_fix_sp3_bug: 1,
        b_fix_fb2: 1,
        b_add_phosphine_stereo: 1,
        b_add_arsine_stereo: 1,
        b_polymers: POLYMERS_NO,
        b_frame_shift_scheme: FSS_STARS_CYCLED,
        ..Default::default()
    };

    // Init table of parameters
    *ip = InputParms::default();

    // Default are StdInChI generation options
    b_ver1_default_mode &= !REQ_MODE_BASIC;
    o.b_reconnect_coord = 0;
    o.b_pointed_edge_stereo = 1;
    ip.b_fix_non_uniform_draw = 1;

    b_name_suffix = true; // bVer1Options
    ip.msec_max_time = 0; // unlimited in the library
    if b_release_version {
        ip.b_abc_numbers = 0;
        ip.b_ct_predecessors = 0;
    }
    ip.b_filter_ss = 0;

    // Analyze command line switches
    let mut i = 1;
    while i < argv.len() {
        let arg = argv[i];
        let ab = arg.as_bytes();
        if ab.first() == Some(&INCHI_OPTION_PREFX) && ab.len() > 1 {
            let p_arg = &arg[1..];
            let pb = p_arg.as_bytes();
            let is = |name: &str| inchi_stricmp(p_arg, name) == 0;
            let prefix = |name: &str| {
                pb.len() >= name.len() && inchi_memicmp(pb, name.as_bytes(), name.len()) == 0
            };
            o.n_mode = n_mode;
            o.b_ver1_default_mode = b_ver1_default_mode;
            let got = set_common_options_by_parg(p_arg, developer_options, ip, &mut o);
            n_mode = o.n_mode;
            b_ver1_default_mode = o.b_ver1_default_mode;

            if got {
                // handled
            } else if is("STDIO") {
                b_name_suffix = false;
            } else if prefix("SDF:") {
                // SDfile label
                let mut hdr = vec![0u8; MAX_SDF_HEADER + 1];
                mystrncpy(&mut hdr, &pb[4..], MAX_SDF_HEADER + 1);
                hdr.truncate(hdr.iter().position(|&c| c == 0).unwrap_or(hdr.len()));
                let k = lrtrim(&mut hdr);
                ip.sz_sdf_data_header = String::from_utf8_lossy(&hdr).to_string();
                if k != 0 {
                    ip.p_sdf_label = true;
                    ip.p_sdf_value = true;
                    if ip.n_input_type == INPUT_NONE {
                        ip.n_input_type = INPUT_SDFILE;
                    }
                } else {
                    ip.p_sdf_label = false;
                    ip.p_sdf_value = false;
                    if ip.n_input_type == INPUT_NONE {
                        ip.n_input_type = INPUT_MOLFILE;
                    }
                }
            } else if is("OUTPUTSDF") {
                b_output_molfile_only = true;
            } else if is("SdfAtomsDT") {
                b_output_molfile_dt = true;
            } else if is("D") {
                // display the structures: not in the library
            } else if prefix("F") && {
                let (_, used, erange) = strtol_range(&pb[1..]);
                errno_erange |= erange;
                used > 0
            } {
                // structure display font size: not used
            } else if is("EQU") {
                // compared components are not used in the library (see below)
            } else if is("WarnOnEmptyStructure") {
                ip.b_allow_empty_structure = 1;
            } else if prefix("W") {
                let mut timeout_value: i64 = 0;
                let c1 = pb.get(1).copied().unwrap_or(0);
                if c1 == b'M' || c1 == b'm' {
                    // "WMnumber", milliseconds
                    let (v, used, erange) = strtol_range(&pb[2..]);
                    errno_erange |= erange;
                    timeout_value = v;
                    if timeout_value != 0 && used > 0 && 2 + used == pb.len() {
                        if errno_erange || timeout_value < 0 {
                            timeout_value = 0;
                            timeout_set_warning = true;
                        }
                        timeout_set_error = false;
                    } else {
                        timeout_set_error = true;
                    }
                } else {
                    // "Wnumber", seconds
                    let (t, used, erange) = strtod_prefix(&pb[1..]);
                    errno_erange |= erange;
                    if t != 0.0 && used > 0 {
                        if 1 + used != pb.len() {
                            timeout_set_warning = true;
                        }
                        if errno_erange || t < 0.0 || t * 1000.0 > LONG_MAX as f64 {
                            timeout_value = 0;
                            timeout_set_warning = true;
                        } else {
                            timeout_value = if t.is_nan() { 0 } else { (t * 1000.0) as i64 };
                        }
                        timeout_set_error = false;
                    } else {
                        timeout_set_error = true;
                    }
                }
                if !timeout_set_error {
                    ip.msec_max_time = timeout_value;
                }
            } else if is("InChI2InChI") {
                ip.n_input_type = INPUT_INCHI;
                ip.b_read_inchi_options |= READ_INCHI_OUTPUT_INCHI;
                ip.b_read_inchi_options &= !READ_INCHI_TO_STRUCTURE;
            } else if is("InChI2Struct") {
                ip.b_read_inchi_options |= READ_INCHI_TO_STRUCTURE;
                ip.b_read_inchi_options &= !READ_INCHI_OUTPUT_INCHI;
                ip.n_input_type = INPUT_INCHI;
            } else if is("DoDrvReport") {
                ip.b_underivatize = 3;
                o.b_std_format = 0;
            } else {
                inchi_ios_eprint(
                    Some(log_file),
                    &format!("Unrecognized optionQ3: \"{p_arg}\".\n"),
                );
            }
        } else if ip.num_paths < MAX_NUM_PATHS && !ab.is_empty() {
            ip.path[ip.num_paths] = Some(arg.to_string());
            ip.num_paths += 1;
        }
        i += 1;
    }

    // Print messages and set control variables according to the options
    if timeout_set_warning {
        inchi_ios_eprint(
            Some(log_file),
            "Warning: timeout value may have been modified (truncated?) due to number formatting issues;\n",
        );
    }
    if timeout_set_error {
        inchi_ios_eprint(
            Some(log_file),
            "Warning: specified timeout value was ignored due to invalid number format, using the default;\n",
        );
    }

    // InChIKey option(s)
    if o.b_hash_key != 0 {
        if ip.b_abc_numbers == 1 && ip.b_ct_predecessors == 1 {
            inchi_ios_eprint(
                Some(log_file),
                "Terminating: generation of InChIKey is not available with 'Compress' option\n",
            );
            return -1;
        }
        if ip.n_input_type == INPUT_INCHI {
            inchi_ios_eprint(
                Some(log_file),
                "Terminating: generation of InChIKey is not available in InChI conversion mode\n",
            );
            return -1;
        } else if b_output_molfile_only {
            inchi_ios_eprint(
                Some(log_file),
                "Terminating: generation of InChIKey is not available with 'OutputSDF' option\n",
            );
            return -1;
        }
    }

    // output file names from the paths
    if b_name_suffix {
        let p: Option<String> = ip.path.iter().flatten().find(|s| !s.is_empty()).cloned();
        if let Some(p_prev) = p.clone() {
            for k in 0..MAX_NUM_PATHS {
                if ip.path[k].as_ref().is_none_or(|s| s.is_empty()) {
                    let mut sz = p_prev.clone();
                    sz.push_str(&sz_name_suffix);
                    sz.push_str(ext[k]);
                    ip.num_paths += 1;
                    ip.path[k] = Some(sz);
                } else if ip.path[k]
                    .as_ref()
                    .is_some_and(|s| inchi_stricmp(s, "NUL") == 0)
                {
                    ip.path[k] = None;
                }
            }
        }
    }

    // inchi2inchi and inchi2struct option(s)
    if ip.n_input_type == INPUT_INCHI {
        b_output_molfile_only = false;
        o.b_inchi_output_options |= INCHI_OUT_NO_AUX_INFO;
        o.b_inchi_output_options &= !INCHI_OUT_SHORT_AUX_INFO;
        o.b_inchi_output_options &= !INCHI_OUT_ONLY_AUX_INFO;
        ip.b_display_if_restore_warnings = b_display_if_restore_warnings;
        if o.b_inchi_output_options
            & (INCHI_OUT_SDFILE_ONLY | INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS)
            == 0
        {
            o.b_inchi_output_options |= INCHI_OUT_PLAIN_TEXT;
        }
    }

    n_mode |= b_ver1_default_mode;

    // COMPILE_ANSI_ONLY: nothing is displayed or compared
    ip.b_display_each_component_inchi = 0;
    let b_compare_components = 0;
    let _ = b_display;

    ip.b_merge_all_input_structures = 0;
    ip.b_do_not_add_h = o.b_do_not_add_h;

    // Set default options
    if n_mode == 0 || n_mode == REQ_MODE_STEREO {
        n_mode |=
            REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO;
    } else {
        if n_mode & (REQ_MODE_BASIC | REQ_MODE_TAUT) == 0 {
            n_mode |= REQ_MODE_BASIC | REQ_MODE_TAUT;
        }
        if n_mode & REQ_MODE_STEREO != 0 && n_mode & (REQ_MODE_ISO | REQ_MODE_NON_ISO) == 0 {
            n_mode |= REQ_MODE_ISO | REQ_MODE_NON_ISO;
        }
    }
    // isotopic requested: unconditionally add non-isotopic output
    if n_mode & REQ_MODE_ISO != 0 {
        n_mode |= REQ_MODE_NON_ISO;
    }
    n_mode |= (n_min_db_rin_size << REQ_MODE_MIN_SB_RING_SHFT) & REQ_MODE_MIN_SB_RING_MASK;

    // Input file
    if ip.n_input_type == INPUT_NONE && ip.num_paths > 0 {
        ip.n_input_type = INPUT_MOLFILE;
    }
    ip.n_mode = n_mode;
    let mut b_compare_components = b_compare_components;
    if b_compare_components & CMP_COMPONENTS != 0 && n_mode & REQ_MODE_BASIC != 0 {
        b_compare_components |= CMP_COMPONENTS_NONTAUT;
    }
    ip.b_compare_components = b_compare_components;
    // Output
    ip.b_inchi_output_options = o.b_inchi_output_options
        | if b_output_molfile_only {
            INCHI_OUT_SDFILE_ONLY
        } else {
            0
        };
    if b_output_molfile_only {
        b_output_style &=
            !(INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS | INCHI_OUT_TABBED_OUTPUT);
        if b_output_molfile_dt {
            ip.b_inchi_output_options |= INCHI_OUT_SDFILE_ATOMS_DT;
        }
        if b_output_molfile_split {
            ip.b_inchi_output_options |= INCHI_OUT_SDFILE_SPLIT;
        }
    }
    ip.b_inchi_output_options |= b_output_style;
    ip.b_no_struct_labels = o.b_no_struct_labels;

    // Processing options
    if o.b_forced_chiral_flag != 0 {
        ip.b_chiral_flag = o.b_forced_chiral_flag;
    }

    // Tautomeric/salts options
    ip.b_taut_flags = 0;
    ip.b_taut_flags_done = 0;
    let mut f = TG_FLAG_TEST_TAUT__ATOMS;
    let flag = |c: bool, v: u64| if c { v } else { 0 };
    f |= flag(b_disconnect_salts, TG_FLAG_DISCONNECT_SALTS);
    f |= flag(b_acid_tautomerism != 0, TG_FLAG_TEST_TAUT__SALTS);
    f |= flag(b_move_positive_charges, TG_FLAG_MOVE_POS_CHARGES);
    f |= flag(b_acid_tautomerism > 1, TG_FLAG_TEST_TAUT2_SALTS);
    f |= flag(b_uncharged_acid_taut, TG_FLAG_ALLOW_NO_NEGTV_O);
    f |= flag(b_merge_salt_t_groups, TG_FLAG_MERGE_TAUT_SALTS);
    f |= flag(b_disconnect_coord, TG_FLAG_DISCONNECT_COORD);
    f |= flag(
        b_disconnect_coord && o.b_reconnect_coord != 0,
        TG_FLAG_RECONNECT_COORD,
    );
    f |= flag(b_disconnect_coord_chk_val, TG_FLAG_CHECK_VALENCE_COORD);
    f |= flag(b_tg_flag_variable_protons, TG_FLAG_VARIABLE_PROTONS);
    f |= flag(b_tg_flag_hard_add_ren_protons, TG_FLAG_HARD_ADD_REM_PROTONS);
    f |= flag(o.b_keto_enol_taut != 0, TG_FLAG_KETO_ENOL_TAUT);
    f |= flag(o.b15_taut_non_ring != 0, TG_FLAG_1_5_TAUT);
    f |= flag(o.b_pt_22_00_taut != 0, TG_FLAG_PT_22_00);
    f |= flag(o.b_pt_16_00_taut != 0, TG_FLAG_PT_16_00);
    f |= flag(o.b_pt_06_00_taut != 0, TG_FLAG_PT_06_00);
    f |= flag(o.b_pt_39_00_taut != 0, TG_FLAG_PT_39_00);
    f |= flag(o.b_pt_13_00_taut != 0, TG_FLAG_PT_13_00);
    f |= flag(o.b_pt_18_00_taut != 0, TG_FLAG_PT_18_00);
    f |= flag(o.b_pointed_edge_stereo != 0, TG_FLAG_POINTED_EDGE_STEREO);
    f |= flag(o.b_add_phosphine_stereo != 0, TG_FLAG_PHOSPHINE_STEREO);
    f |= flag(o.b_add_arsine_stereo != 0, TG_FLAG_ARSINE_STEREO);
    f |= flag(o.b_fix_sp3_bug != 0, TG_FLAG_FIX_SP3_BUG);
    ip.b_taut_flags = f;

    // Bug fixes
    if o.b_fix_fb2 != 0 {
        ip.b_taut_flags |= TG_FLAG_FIX_ISO_FIXEDH_BUG;
        ip.b_taut_flags |= TG_FLAG_FIX_TERM_H_CHRG_BUG;
        ip.b_inchi_output_options |= INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG;
    }

    if ip.n_input_type == 0 {
        ip.n_input_type = INPUT_MOLFILE;
    }

    // Check if /SNon requested turn OFF SUU/SLUUD
    if ip.n_mode & REQ_MODE_STEREO == 0 {
        ip.n_mode &= !REQ_MODE_DIFF_UU_STEREO;
        ip.n_mode &= !(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
    }

    // Standard InChI ?
    if o.b_std_format != 0 {
        ip.b_inchi_output_options |= INCHI_OUT_STDINCHI;
    }

    // InChIKey ?
    if o.b_hash_key == 0 {
        ip.b_calc_inchi_hash = INCHIHASH_NONE;
        if o.b_hash_xtra1 != 0 || o.b_hash_xtra2 != 0 {
            inchi_ios_eprint(
                Some(log_file),
                "Hash extension(s) not generated: InChIKey not requested",
            );
        }
    } else {
        ip.b_calc_inchi_hash = INCHIHASH_KEY;
        if o.b_hash_xtra1 != 0 {
            ip.b_calc_inchi_hash = if o.b_hash_xtra2 != 0 {
                INCHIHASH_KEY_XTRA1_XTRA2
            } else {
                INCHIHASH_KEY_XTRA1
            };
        } else if o.b_hash_xtra2 != 0 {
            ip.b_calc_inchi_hash = INCHIHASH_KEY_XTRA2;
        }
    }

    ip.b_large_molecules = o.b_large_molecules;
    ip.b_loose_tsa_check = o.b_loose_tsa_check;
    ip.b_np_zz = o.b_np_zz;
    ip.b_stereo_at_zz = o.b_stereo_at_zz;
    ip.b_no_warnings = o.b_no_warnings;
    ip.b_merge_hash = o.b_merge_hash;
    ip.b_hide_inchi = o.b_hide_inchi;
    ip.b_polymers = o.b_polymers;
    ip.b_fold_polymer_sru = o.b_fold_polymer_sru;
    ip.b_frame_shift_scheme = o.b_frame_shift_scheme;
    ip.b_inchi_output_options2 = o.b_inchi_output_options2;
    0
}

/// The "Ensure standardness" block of `GetINCHIEx` (`inchi_dll.c`): when the
/// caller asks for a standard InChI, every option that would make the result
/// non-standard is switched off.
pub fn enforce_std_format(ip: &mut InputParms) {
    ip.b_inchi_output_options &= !INCHI_OUT_SAVEOPT;
    ip.b_taut_flags &= !TG_FLAG_RECONNECT_COORD;
    ip.n_mode &= !REQ_MODE_BASIC;
    if ip.n_mode & REQ_MODE_RELATIVE_STEREO != 0 {
        ip.n_mode &=
            !(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
    }
    if ip.n_mode & REQ_MODE_RACEMIC_STEREO != 0 {
        ip.n_mode &=
            !(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
    }
    if ip.n_mode & REQ_MODE_CHIR_FLG_STEREO != 0 {
        ip.n_mode &=
            !(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
    }
    ip.n_mode &= !REQ_MODE_DIFF_UU_STEREO;
    if ip.n_mode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) == 0 {
        ip.n_mode |= REQ_MODE_SB_IGN_ALL_UU;
        ip.n_mode |= REQ_MODE_SC_IGN_ALL_UU;
    }
    ip.b_taut_flags &= !TG_FLAG_KETO_ENOL_TAUT;
    ip.b_taut_flags &= !TG_FLAG_1_5_TAUT;
    // And anyway...
    ip.b_inchi_output_options |= INCHI_OUT_STDINCHI;
    ip.b_inchi_output_options &= !INCHI_OUT_SAVEOPT;
}
