//! Reading and component extraction of `INCHI_BASE/src/runichi2.c`
//! (`GetOneComponent`, the error classification of `TreatErrorsInReadTheStructure`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichierr::{CT_ATOMCOUNT_ERR, CT_UNKNOWN_ERR, add_error_message};
use super::ichimake::InpAtomData;
use super::ichimake::get_inp_struct_error_type;
use super::ichiparm::InputParms;
use super::runichi::create_inp_atom_data;
use super::runichi3::{IS_ERROR, IS_FATAL, IS_OKAY, IS_SKIP, StructData};
use super::strutil::{OrigAtomData, extract_connected_component};

/// `GetOneComponent`: extracts the component `i` (from 0) of `orig` into
/// `inp_cur_data`; returns `sd.n_error_type`.
pub fn get_one_component(
    sd: &mut StructData,
    inp_cur_data: &mut InpAtomData,
    orig: &OrigAtomData,
    i: usize,
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

/// `TreatErrorsInReadTheStructure` (without the logging and the problem-file
/// copies); `num_inp` is the number of the structure.
pub fn treat_errors_in_read_the_structure(
    sd: &mut StructData,
    ip: &InputParms,
    orig: &OrigAtomData,
    num_inp: i64,
) -> i32 {
    let mut ret = IS_OKAY;
    'body: {
        if 10 < sd.n_struct_read_error && sd.n_struct_read_error < 20 {
            ret = super::runichi3::IS_EOF;
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
        let _ = IS_FATAL;
    }
    if ret <= IS_OKAY && sd.n_error_type > 0 {
        ret = sd.n_error_type;
    }
    ret
}
