//! Per-structure result treatment (`INCHI_BASE/src/runichi4.c`): the
//! processing warnings of one component.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::InChI;
use super::ichierr::add_error_message;
use super::ichimak2::{
    AMBIGUOUS_STEREO_ATOM, AMBIGUOUS_STEREO_ATOM_ISO, AMBIGUOUS_STEREO_BOND,
    AMBIGUOUS_STEREO_BOND_ISO,
};
use super::ichimake::InpAtomData;

/// `GetProcessingWarningsOneInChI`: 1 when there is ambiguous stereo
pub fn get_processing_warnings_one_inchi(
    inchi: &InChI,
    inp_norm_data: &InpAtomData,
    str_err_struct: &mut String,
    b_no_warnings: bool,
) -> i32 {
    let (mut atoms, mut bonds) = (0, 0);
    if !inp_norm_data.at.is_empty() {
        for a in &inp_norm_data.at[..inchi.n_number_of_atoms.max(0) as usize] {
            if a.b_ambiguous_stereo & (AMBIGUOUS_STEREO_ATOM | AMBIGUOUS_STEREO_ATOM_ISO) != 0 {
                atoms += 1;
            }
            if a.b_ambiguous_stereo & (AMBIGUOUS_STEREO_BOND | AMBIGUOUS_STEREO_BOND_ISO) != 0 {
                bonds += 1;
            }
        }
        if atoms != 0 && !b_no_warnings {
            add_error_message(str_err_struct, "Ambiguous stereo:");
            add_error_message(str_err_struct, "center(s)");
        }
        if bonds != 0 && !b_no_warnings {
            add_error_message(str_err_struct, "Ambiguous stereo:");
            add_error_message(str_err_struct, "bond(s)");
        }
    }
    i32::from(atoms != 0 || bonds != 0)
}

/// `GetProcessingWarningsOneComponentInChI`
pub fn get_processing_warnings_one_component_inchi(
    cur_inchi: &[Option<InChI>; 2],
    inp_norm_data: &[InpAtomData; 2],
    str_err_struct: &mut String,
    b_no_warnings: bool,
) -> i32 {
    let mut ret = 0;
    for i in 0..2 {
        if let Some(c) = cur_inchi[i].as_ref().filter(|c| c.n_number_of_atoms > 0) {
            ret |= get_processing_warnings_one_inchi(
                c,
                &inp_norm_data[i],
                str_err_struct,
                b_no_warnings,
            );
        }
    }
    ret
}

/// `TreatErrorsInCreateOneComponentINChI` (without the logging): sets the
/// error message and `nErrorType`; returns it
pub fn treat_errors_in_create_one_component_inchi(sd: &mut super::runichi3::StructData) -> i32 {
    use super::ichierr::{CT_OUT_OF_RAM, CT_USER_QUIT_ERR, err_msg};
    use super::runichi3::{IS_ERROR, IS_FATAL};
    if sd.n_error_code != 0 {
        add_error_message(&mut sd.str_err_struct, &err_msg(sd.n_error_code));
        sd.n_error_type = if sd.n_error_code == CT_OUT_OF_RAM || sd.n_error_code == CT_USER_QUIT_ERR
        {
            IS_FATAL
        } else {
            IS_ERROR
        };
    }
    sd.n_error_type
}
