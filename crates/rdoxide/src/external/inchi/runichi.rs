//! The per-component driver of `INCHI_BASE/src/runichi.c`
//! (`CreateOneComponentINChI`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

//! The elapsed-time bookkeeping and the interactive display of the C code are
//! not modelled; there is no time limit.

use super::ichi::{InChI, InChIAux};
use super::ichierr::{CT_OUT_OF_RAM, CT_USER_QUIT_ERR};
use super::ichimake::{InpAtomData, create_inchi};
use super::ichiparm::{
    INCHI_OUT_NO_AUX_INFO, INCHI_OUT_SHORT_AUX_INFO, InputParms, REQ_MODE_BASIC, REQ_MODE_ISO,
    REQ_MODE_TAUT,
};
use super::inpdef::InpAtom;
use super::runichi3::{
    IS_ERROR, IS_FATAL, StructData, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE,
};
use super::runichi4::get_processing_warnings_one_component_inchi;
use super::strutil::{alloc_inchi, alloc_inchi_aux, set_connected_component_number};

/// `NORM_CANON_FLAGS` for `INCHI_NUM` = 3 structures and the two tautomer modes
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NormCanonFlags {
    pub b_taut_flags: [[u64; 2]; 3],
    pub b_taut_flags_done: [[u64; 2]; 3],
    pub b_normalization_flags: [[u64; 2]; 3],
    pub n_canon_flags: [[i32; 2]; 3],
}

/// `CreateInpAtomData`
pub fn create_inp_atom_data(d: &mut InpAtomData, num_atoms: usize, create_at_fixed_bonds: bool) {
    *d = InpAtomData::default();
    d.at = vec![InpAtom::default(); num_atoms];
    if create_at_fixed_bonds {
        d.at_fixed_bonds = Some(vec![InpAtom::default(); num_atoms]);
    }
    d.num_at = num_atoms as i32;
}

/// `CreateOneComponentINChI`: canonicalizes the component `inp_cur_data` and
/// stores the results in `p_inchi` / `p_aux`; `i_inchi` selects the
/// `NORM_CANON_FLAGS` row (`INCHI_BAS` = 0, ...).
#[allow(clippy::too_many_arguments)]
pub fn create_one_component_inchi(
    sd: &mut StructData,
    ip: &InputParms,
    inp_cur_data: &mut InpAtomData,
    p_inchi: &mut [Option<InChI>; 2],
    p_aux: &mut [Option<InChIAux>; 2],
    i_inchi: usize,
    component_number: usize,
    inp_norm_data: &mut [InpAtomData; 2],
    nc_flags: &mut NormCanonFlags,
) -> i32 {
    let mut taut_flags = ip.b_taut_flags;
    let mut taut_flags_done = ip.b_taut_flags_done | sd.taut_flags_done[0];
    let _b_orig_coord =
        ip.b_inchi_output_options & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO) == 0;
    let mut cur_inchi: [Option<InChI>; 2] = [None, None];
    let mut cur_aux: [Option<InChIAux>; 2] = [None, None];
    for k in 0..2 {
        let mut alloc_mode = if k == 1 { REQ_MODE_TAUT } else { 0 };
        if taut_flags_done & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE) != 0
        {
            alloc_mode |= ip.n_mode & REQ_MODE_ISO;
        }
        if (k == 0 && ip.n_mode & REQ_MODE_BASIC != 0) || (k == 1 && ip.n_mode & REQ_MODE_TAUT != 0)
        {
            let n = inp_cur_data.num_at.max(0) as usize;
            if let Some((p, nb, ni)) = alloc_inchi(&inp_cur_data.at, n, alloc_mode) {
                inp_cur_data.num_bonds = nb as i32;
                inp_cur_data.num_isotopic = ni as i32;
                cur_inchi[k] = Some(p);
                cur_aux[k] = alloc_inchi_aux(n, ni, alloc_mode);
            }
            create_inp_atom_data(&mut inp_norm_data[k], n, k == 1);
        } else {
            inp_norm_data[k] = InpAtomData::default();
        }
    }

    let num_inp = inp_cur_data.num_at.max(0) as usize;
    let num_at = create_inchi(
        &mut cur_inchi,
        &mut cur_aux,
        &inp_cur_data.at,
        inp_norm_data,
        num_inp,
        ip.n_mode,
        &mut taut_flags,
        &mut taut_flags_done,
        None,
        &mut sd.str_err_struct,
        ip,
    );
    set_connected_component_number(&mut inp_cur_data.at, (component_number + 1) as u16);

    for k in 0..2 {
        if let Some(a) = cur_aux[k].as_ref().filter(|a| a.n_number_of_atoms > 0) {
            nc_flags.b_normalization_flags[i_inchi][k] |= a.b_normalization_flags;
            nc_flags.b_taut_flags[i_inchi][k] |= a.b_taut_flags;
            nc_flags.b_taut_flags_done[i_inchi][k] |= a.b_taut_flags_done;
            nc_flags.n_canon_flags[i_inchi][k] |= a.n_canon_flags;
        }
    }

    if num_at < 0 {
        sd.n_error_code = num_at;
    } else if num_at == 0 {
        sd.n_error_code = -1;
    } else if let Some(c) = cur_inchi[0].as_ref().filter(|c| c.n_error_code != 0) {
        sd.n_error_code = c.n_error_code;
    } else if let Some(c) = cur_inchi[1].as_ref().filter(|c| c.n_error_code != 0) {
        sd.n_error_code = c.n_error_code;
    }
    if sd.n_error_code == 0 {
        get_processing_warnings_one_component_inchi(
            &cur_inchi,
            inp_norm_data,
            &mut sd.str_err_struct,
            ip.b_no_warnings != 0,
        );
    }

    *p_inchi = cur_inchi;
    *p_aux = cur_aux;

    if sd.n_error_code == 0 {
        let in_non_taut = p_inchi[0].as_ref().is_some_and(|c| c.n_number_of_atoms > 0);
        let in_taut = p_inchi[1].as_ref().is_some_and(|c| c.n_number_of_atoms > 0);
        let is_non_taut = (in_non_taut && p_inchi[0].as_ref().unwrap().len_tautomer == 0)
            || (in_taut && p_inchi[1].as_ref().unwrap().len_tautomer == 0);
        let is_taut = in_taut && p_inchi[1].as_ref().unwrap().len_tautomer > 0;
        if is_non_taut || is_taut {
            let j1 = if in_non_taut { 0 } else { 1 };
            let j2 = if in_taut { 1 } else { 0 };
            sd.num_non_taut[i_inchi] += i32::from(is_non_taut);
            sd.num_taut[i_inchi] += i32::from(is_taut);
            for j in j1..=j2 {
                let c = p_inchi[j].as_ref().unwrap();
                let mut iso = c.n_number_of_isotopic_atoms != 0
                    || c.n_number_of_isotopic_t_groups != 0
                    || c.n_possible_locations_of_isotopic_h
                        .first()
                        .is_some_and(|&x| x > 1);
                if j == 1 {
                    if let Some(a) = p_aux[j].as_ref() {
                        iso |= a
                            .n_num_removed_isotopic_h
                            .iter()
                            .map(|&x| i32::from(x))
                            .sum::<i32>()
                            > 0;
                    }
                }
                inp_norm_data[j].b_exists = 1;
                inp_norm_data[j].b_has_isotopic_layer = i32::from(iso);
            }
        }
    }

    if sd.n_error_code == CT_OUT_OF_RAM || sd.n_error_code == CT_USER_QUIT_ERR {
        IS_FATAL
    } else if sd.n_error_code != 0 {
        IS_ERROR
    } else {
        0
    }
}
