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

// ---------------------------------------------------------------------------
// The main serialization procedure

use super::ichi::{InChIAux, InChISort, comp_inchi_non_taut2, comp_inchi_taut2};
use super::ichi_io::{InchiIosString, InchiIostream, inchi_ios_print};
use super::ichiprt1::{
    INCHI_OUT_EMBED_REC, INCHI_OUT_PLAIN_TEXT, INCHI_OUT_PLAIN_TEXT_COMMENTS,
    INCHI_OUT_TABBED_OUTPUT, OrigStruct, OutputInput, output_inchi2,
};
use super::ichiprt3::OUT_TN;
use super::runichi3::TG_FLAG_DISCONNECT_COORD_DONE;

const INCHI_BAS: usize = 0;
const TAUT_NON: usize = 0;
const TAUT_YES: usize = 1;
/// `INCHI_OUT_PRINT_OPTIONS`
const INCHI_OUT_PRINT_OPTIONS: i32 =
    INCHI_OUT_EMBED_REC | INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
const CMP_COMPONENTS: i32 = 0x0001;
const IS_FATAL: i32 = 3;

/// the per-component results: `pINChI[iINCHI][component][bTaut]`
pub type PInChI2 = Vec<[Option<InChI>; 2]>;
/// `pINChI_Aux[iINCHI][component][bTaut]`
pub type PInChIAux2 = Vec<[Option<InChIAux>; 2]>;

/// the parts of `INPUT_PARMS` used by `SortAndPrintINChI`
#[derive(Debug, Clone, Default)]
pub struct SortPrintParams {
    pub inchi_output_options: i32,
    pub compare_components: i32,
    pub input: OutputInput,
}

/// `bIsStructChiral`
pub fn is_struct_chiral(pinchi: &[PInChI2; 2], num_components: [i32; 2]) -> bool {
    for j in 0..2 {
        for comp in pinchi[j].iter().take(num_components[j].max(0) as usize) {
            for c in comp.iter().flatten() {
                if c.b_deleted == 0 && c.n_number_of_atoms > 0 {
                    for s in [c.stereo.as_ref(), c.stereo_isotopic.as_ref()]
                        .into_iter()
                        .flatten()
                    {
                        if !s.t_parity.is_empty()
                            && s.n_number_of_stereo_centers > 0
                            && s.n_comp_inv2_abs != 0
                        {
                            return true;
                        }
                    }
                }
            }
        }
    }
    false
}

/// `SortAndPrintINChI`: sorts the components of the mobile-H and fixed-H
/// orders and either prints the InChI strings or stores the components in
/// the mobile-H order into `pinchi`/`pinchi_aux`; returns 0 or `_IS_FATAL`.
/// (The component-equivalence labels of `SaveEquComponentsInfoAndSortOrder`
/// are only produced for the GUI builds, where the structure was saved in
/// the library; they are not ported.)
#[allow(clippy::too_many_arguments)]
pub fn sort_and_print_inchi(
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    mut log_file: Option<&mut InchiIostream>,
    ip: &SortPrintParams,
    orig_struct: Option<&OrigStruct>,
    num_components: [i32; 2],
    taut_flags: &mut [u64; 2],
    taut_flags_done: &mut [u64; 2],
    nc_taut_flags: [[u64; 2]; 2],
    nc_taut_flags_done: [[u64; 2]; 2],
    num_inp: i32,
    pinchi: &mut [PInChI2; 2],
    pinchi_aux: &mut [PInChIAux2; 2],
    flags: &mut i32,
    save_opt_bits: u8,
) -> i32 {
    let disconnected_coord = taut_flags_done[0] & TG_FLAG_DISCONNECT_COORD_DONE != 0;
    let mut ret = 1;
    for i in 0..2 {
        for k in 0..2 {
            taut_flags[i] |= nc_taut_flags[i][k];
            taut_flags_done[i] |= nc_taut_flags_done[i][k];
        }
    }
    let mut sorts: [[Vec<InChISort>; 2]; 2] = Default::default();
    for j in 0..2 {
        if num_components[j] == 0 {
            continue;
        }
        let n = num_components[j] as usize;
        for k1 in 0..2 {
            let mut v: Vec<InChISort> = (0..n)
                .map(|i| InChISort {
                    p_inchi: pinchi[j][i].clone(),
                    p_inchi_aux: pinchi_aux[j][i].clone(),
                    ord_number: i as i16,
                })
                .collect();
            v.sort_by(|a, b| {
                let r = if k1 == TAUT_NON {
                    comp_inchi_non_taut2(a, b)
                } else {
                    comp_inchi_taut2(a, b)
                };
                r.cmp(&0)
            });
            sorts[j][k1] = v;
        }
    }

    if ip.inchi_output_options & INCHI_OUT_PRINT_OPTIONS == 0 {
        for j in 0..2 {
            if num_components[j] == 0 {
                continue;
            }
            for i in 0..num_components[j] as usize {
                let s = &sorts[j][TAUT_YES][i];
                if s.p_inchi[TAUT_NON].is_some() && s.p_inchi[TAUT_YES].is_none() {
                    for k in 0..2 {
                        pinchi[j][i][k] = s.p_inchi[1 - k].clone();
                        pinchi_aux[j][i][k] = s.p_inchi_aux[1 - k].clone();
                    }
                } else {
                    for k in 0..2 {
                        pinchi[j][i][k] = s.p_inchi[k].clone();
                        pinchi_aux[j][i][k] = s.p_inchi_aux[k].clone();
                    }
                }
            }
        }
    } else {
        let options0 = ip.inchi_output_options & !INCHI_OUT_PRINT_OPTIONS;
        let embed_reconnected = ip.inchi_output_options & INCHI_OUT_EMBED_REC;
        for i in 1..=2 {
            let cur = if i == 2 {
                INCHI_OUT_PLAIN_TEXT_COMMENTS
            } else {
                INCHI_OUT_PLAIN_TEXT
            };
            if ip.inchi_output_options & cur == 0 {
                continue;
            }
            let mut options_cur = options0 | cur;
            if i == 1 {
                options_cur |= embed_reconnected;
            } else {
                inchi_ios_print(Some(out_file), "\n==== InChI ANNOTATED CONTENTS ====\n");
                options_cur |= embed_reconnected;
                options_cur &= !INCHI_OUT_TABBED_OUTPUT;
            }
            let ret2 = output_inchi2(
                strbuf,
                &sorts,
                INCHI_BAS as i32,
                orig_struct,
                &ip.input,
                disconnected_coord,
                OUT_TN,
                options_cur,
                num_components,
                out_file,
                log_file.as_deref_mut(),
                num_inp,
                flags,
                save_opt_bits,
            );
            ret &= ret2;
            if ret == 0 {
                break;
            }
        }
    }
    let _ = CMP_COMPONENTS;
    if ret != 0 { 0 } else { IS_FATAL }
}
