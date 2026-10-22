//! Reading and component extraction of `INCHI_BASE/src/runichi2.c`
//! (`GetOneComponent`, the error classification of `TreatErrorsInReadTheStructure`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi_io::{InchiIostream, inchi_ios_eprint};
use super::ichierr::{CT_ATOMCOUNT_ERR, CT_UNKNOWN_ERR, add_error_message};
use super::ichimake::InpAtomData;
use super::ichimake::get_inp_struct_error_type;
use super::ichiparm::{
    FLAG_SET_INP_AT_CHIRAL, FLAG_SET_INP_AT_NONCHIRAL, INCHI_OUT_NO_AUX_INFO,
    INCHI_OUT_SDFILE_ONLY, INCHI_OUT_SHORT_AUX_INFO, INPUT_MOLFILE, InputParms,
    REQ_MODE_CHIR_FLG_STEREO, REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO, REQ_MODE_STEREO,
};
use super::mol2atom::{
    FLAG_INP_AT_CHIRAL, FLAG_INP_AT_NONCHIRAL, FLAG_SET_INP_LARGE_MOLS, MolfileAtomsArgs,
    create_orig_inp_data_from_molfile,
};
use super::runichi::create_inp_atom_data;
use super::runichi3::{IS_EOF, IS_ERROR, IS_FATAL, IS_OKAY, IS_SKIP, IS_WARNING, StructData};
use super::strutil::{OrigAtomData, extract_connected_component};

pub const LOG_MASK_WARN: u32 = 1;
pub const LOG_MASK_ERR: u32 = 2;
pub const LOG_MASK_FATAL: u32 = 4;
pub const LOG_MASK_ALL: u32 = 7;
pub const LOG_MASK_NO_WARN: u32 = 6;

/// `GetOneComponent`: extracts the component `i` (from 0) of `orig` into
/// `inp_cur_data`; returns `sd.n_error_type`. (The structure labels of the
/// log line are empty: the library API reads no SDfile labels.)
pub fn get_one_component(
    sd: &mut StructData,
    log: &mut InchiIostream,
    inp_cur_data: &mut InpAtomData,
    orig: &OrigAtomData,
    i: usize,
    num_inp: i64,
) -> i32 {
    let cur_len = orig
        .n_cur_at_len
        .as_ref()
        .and_then(|v| v.get(i).copied())
        .unwrap_or(0) as usize;
    create_inp_atom_data(inp_cur_data, cur_len, false);
    let n = extract_connected_component(
        &orig.at[..orig.num_inp_atoms],
        (i + 1) as u16,
        &mut inp_cur_data.at,
    );
    inp_cur_data.num_at = n;
    if n <= 0 || cur_len as i32 != n {
        add_error_message(&mut sd.str_err_struct, "Cannot extract Component");
        inchi_ios_eprint(
            Some(log),
            &format!("{} #{} structure #{}.\n", sd.str_err_struct, i + 1, num_inp),
        );
        sd.n_error_code = if n < 0 {
            n
        } else if cur_len as i32 != n {
            CT_ATOMCOUNT_ERR
        } else {
            CT_UNKNOWN_ERR
        };
        sd.n_error_type = IS_ERROR;
    }
    sd.n_error_type
}

/// `TreatErrorsInReadTheStructure` (the
/// problem-file copies are not modelled); `num_inp` is the number of the
/// structure.
pub fn treat_errors_in_read_the_structure(
    sd: &mut StructData,
    ip: &InputParms,
    log_mask: u32,
    log: &mut InchiIostream,
    orig: &OrigAtomData,
    num_inp: i64,
) -> i32 {
    let mut ret = IS_OKAY;
    'body: {
        if 10 < sd.n_struct_read_error && sd.n_struct_read_error < 20 {
            // end of file
            if !sd.str_err_struct.is_empty() {
                inchi_ios_eprint(
                    Some(log),
                    &format!(
                        "{} inp structure #{}: End of file.    \n",
                        sd.str_err_struct, num_inp
                    ),
                );
            }
            inchi_ios_eprint(
                Some(log),
                &format!(
                    "End of file detected after structure #{}.   \n",
                    num_inp - 1
                ),
            );
            ret = IS_EOF;
            break 'body;
        }
        if num_inp < ip.first_struct_number {
            sd.n_error_type = IS_SKIP;
            ret = IS_SKIP;
            break 'body;
        }
        sd.n_error_type = get_inp_struct_error_type(
            ip,
            sd.n_struct_read_error,
            &sd.str_err_struct,
            orig.num_inp_atoms,
        );
        if sd.n_error_type == IS_FATAL && log_mask & LOG_MASK_FATAL != 0 {
            inchi_ios_eprint(
                Some(log),
                &format!(
                    "Fatal Error {} (aborted; {}) inp structure #{}.\n",
                    sd.n_struct_read_error, sd.str_err_struct, num_inp
                ),
            );
        }
        if sd.n_error_type == IS_ERROR && log_mask & LOG_MASK_ERR != 0 {
            inchi_ios_eprint(
                Some(log),
                &format!(
                    "Error {} (no {}; {}) inp structure #{}.\n",
                    sd.n_struct_read_error,
                    if ip.b_inchi_output_options & INCHI_OUT_SDFILE_ONLY != 0 {
                        "Molfile"
                    } else {
                        "InChI"
                    },
                    sd.str_err_struct,
                    num_inp
                ),
            );
        }
        if sd.n_error_type == IS_WARNING && log_mask & LOG_MASK_WARN != 0 {
            inchi_ios_eprint(
                Some(log),
                &format!(
                    "Warning: ({}) inp structure #{}.\n",
                    sd.str_err_struct, num_inp
                ),
            );
        }
    }
    if ret <= IS_OKAY && sd.n_error_type > 0 {
        ret = sd.n_error_type;
    }
    ret
}

/// `ReadTheStructure` for a Molfile input (`INPUT_MOLFILE`; any other input
/// type gives `_IS_FATAL`) reading into `orig`
pub fn read_the_structure(
    sd: &mut StructData,
    ip: &mut InputParms,
    inp_file: &mut InchiIostream,
    orig: &mut OrigAtomData,
) -> i32 {
    let get_orig_coord =
        ip.b_inchi_output_options & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO) == 0;
    let mut flags: u64 = if ip.b_large_molecules != 0 {
        FLAG_SET_INP_LARGE_MOLS
    } else {
        0
    };
    *sd = StructData::default();
    if ip.n_input_type != INPUT_MOLFILE {
        return IS_FATAL;
    }
    let mut od = Some(std::mem::take(orig));
    let mut molfile_number = ip.l_molfile_number;
    create_orig_inp_data_from_molfile(
        inp_file,
        &mut od,
        ip.b_merge_all_input_structures != 0,
        get_orig_coord,
        ip.b_do_not_add_h != 0,
        ip.b_polymers,
        ip.b_np_zz,
        MolfileAtomsArgs {
            sdf_label: None,
            sdf_value: None,
            id: Some(&mut ip.l_sdf_id),
            inp_atom_flags: Some(&mut flags),
        },
        &mut molfile_number,
        &mut sd.n_struct_read_error,
        &mut sd.str_err_struct,
        ip.b_no_warnings != 0,
    );
    ip.l_molfile_number = molfile_number;
    *orig = od.unwrap_or_default();
    if ip.b_get_sdfile_id == 0 || ip.l_sdf_id == 999_999 {
        ip.l_sdf_id = 0;
    }
    if ip.b_get_molfile_number == 0 || ip.l_molfile_number < 0 {
        ip.l_molfile_number = 0;
    }
    // Molfile chiral flag mode
    if ip.b_chiral_flag & FLAG_SET_INP_AT_CHIRAL != 0 {
        flags = FLAG_INP_AT_CHIRAL;
    } else if ip.b_chiral_flag & FLAG_SET_INP_AT_NONCHIRAL != 0 {
        flags = FLAG_INP_AT_NONCHIRAL;
    } else if flags & FLAG_INP_AT_CHIRAL != 0 && flags & FLAG_INP_AT_NONCHIRAL != 0 {
        flags &= !FLAG_INP_AT_NONCHIRAL;
    }
    sd.b_chiral_flag &= !((FLAG_INP_AT_CHIRAL | FLAG_INP_AT_NONCHIRAL) as i32);
    sd.b_chiral_flag |= (flags & (FLAG_INP_AT_CHIRAL | FLAG_INP_AT_NONCHIRAL)) as i32;
    if ip.n_mode & REQ_MODE_CHIR_FLG_STEREO != 0 && ip.n_mode & REQ_MODE_STEREO != 0 {
        if flags & FLAG_INP_AT_CHIRAL != 0 {
            ip.n_mode &= !(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
            sd.b_chiral_flag |= FLAG_INP_AT_CHIRAL as i32;
        } else {
            ip.n_mode &= !REQ_MODE_RACEMIC_STEREO;
            ip.n_mode |= REQ_MODE_RELATIVE_STEREO;
            sd.b_chiral_flag |= FLAG_INP_AT_NONCHIRAL as i32;
        }
    }
    0
}

/// `GetOneStructure` for a Molfile input without the file-position table;
/// `num_inp` is updated as in C
pub fn get_one_structure(
    sd: &mut StructData,
    ip: &mut InputParms,
    log: &mut InchiIostream,
    inp_file: &mut InchiIostream,
    orig: &mut OrigAtomData,
    num_inp: &mut i64,
) -> i32 {
    *orig = OrigAtomData::default();
    let mut n_ret = read_the_structure(sd, ip, inp_file, orig);
    if n_ret == 0 {
        if ip.l_molfile_number != 0 {
            *num_inp = ip.l_molfile_number;
        } else {
            *num_inp += 1;
        }
        n_ret = treat_errors_in_read_the_structure(sd, ip, LOG_MASK_ALL, log, orig, *num_inp);
    }
    if n_ret == IS_EOF {
        *num_inp -= 1;
    }
    n_ret
}
