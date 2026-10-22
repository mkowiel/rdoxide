//! The per-component driver of `INCHI_BASE/src/runichi.c`
//! (`CreateOneComponentINChI`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

//! The elapsed-time bookkeeping and the interactive display of the C code are
//! not modelled; there is no time limit.

use super::ichi::{InChI, InChIAux};
use super::ichi_io::{InchiIosString, InchiIostream, inchi_ios_eprint};
use super::ichierr::{CT_OUT_OF_RAM, CT_USER_QUIT_ERR};
use super::ichimake::{InpAtomData, create_inchi};
use super::ichiparm::{
    INCHI_OUT_NO_AUX_INFO, INCHI_OUT_SHORT_AUX_INFO, InputParms, REQ_MODE_BASIC, REQ_MODE_ISO,
    REQ_MODE_TAUT,
};
use super::inpdef::InpAtom;
use super::runichi3::{
    INCHI_BAS, IS_ERROR, IS_FATAL, StructData, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE,
    TG_FLAG_FOUND_ISOTOPIC_H_DONE,
};
use super::runichi4::{PInChI2, PInChIAux2, get_processing_warnings_one_component_inchi};
use super::strutil::OrigAtomData;
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
                if let Some(a) = cur_aux[k].as_mut() {
                    a.b_is_isotopic = ni as i32;
                }
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

/// The results of `CreateOneStructureINChI`: for each of the three layers
/// (`INCHI_BAS`, `INCHI_REC`, ...) the `INChI` / `INChI_Aux` pairs of its
/// components.
#[derive(Debug, Clone, Default)]
pub struct StructResults {
    pub p_inchi: [Vec<[Option<InChI>; 2]>; 3],
    pub p_aux: [Vec<[Option<InChIAux>; 2]>; 3],
}

/// `CreateOneStructureINChI`: preprocesses the structure (once), splits it
/// into components and creates the `INChI` of each. The displays, the
/// composite-structure data, the logging and the problem-file copies of the C
/// code are not modelled. `INChI` results shared between layers
/// (`bREUSE_INCHI`) are copies here.
pub fn create_one_structure_inchi(
    sd: &mut StructData,
    ip: &InputParms,
    log: &mut InchiIostream,
    num_inp: i64,
    i_inchi: usize,
    orig: &mut OrigAtomData,
    prep: &mut [OrigAtomData; 2],
    res: &mut StructResults,
    nc_flags: &mut NormCanonFlags,
) -> i32 {
    use super::ichierr::add_error_message;
    use super::ichiparm::CMP_COMPONENTS;
    use super::runichi2::{get_one_component, treat_errors_in_read_the_structure};
    use super::runichi3::preprocess_one_structure;
    use super::runichi4::treat_errors_in_create_one_component_inchi;

    let empty_ok = ip.b_allow_empty_structure != 0 && orig.at.is_empty() && orig.num_inp_atoms == 0;
    if !empty_ok && (orig.at.is_empty() || orig.num_inp_atoms == 0) {
        return 0;
    }
    if i_inchi == 1 && orig.b_disconnect_coord <= 0 {
        return 0;
    }
    if i_inchi != INCHI_BAS && i_inchi != 1 {
        add_error_message(&mut sd.str_err_struct, "Fatal undetermined program error");
        sd.n_struct_read_error = 97;
        sd.n_error_type = IS_FATAL;
        return IS_FATAL;
    }
    let mut n_ret = 0;
    'cycle: {
        if (prep[0].at.is_empty() || prep[0].num_inp_atoms == 0) && orig.num_inp_atoms > 0 {
            preprocess_one_structure(sd, ip, orig, prep);
            let f = sd.taut_flags[INCHI_BAS] | ip.b_taut_flags;
            nc_flags.b_taut_flags[i_inchi] = [f, f];
            let f = sd.taut_flags_done[INCHI_BAS] | ip.b_taut_flags_done;
            nc_flags.b_taut_flags_done[i_inchi] = [f, f];
            if sd.n_error_type == IS_ERROR || sd.n_error_type == IS_FATAL {
                n_ret = treat_errors_in_read_the_structure(
                    sd,
                    ip,
                    super::runichi2::LOG_MASK_ALL,
                    log,
                    &prep[0],
                    num_inp,
                );
                break 'cycle;
            }
        }
        let ncomp = prep[i_inchi].num_components.max(0) as usize;
        if sd.num_components[i_inchi] as usize <= ncomp {
            res.p_inchi[i_inchi].resize(ncomp + 1, [None, None]);
            res.p_aux[i_inchi].resize(ncomp + 1, [None, None]);
            sd.num_components[i_inchi] = ncomp as i32;
        }
        let mut inp_cur_data = InpAtomData::default();
        let mut inp_norm_data = [InpAtomData::default(), InpAtomData::default()];
        for i in 0..ncomp {
            if sd.b_user_quit_component != 0 {
                break;
            }
            if i_inchi == 1 && ip.b_compare_components & CMP_COMPONENTS == 0 {
                let m = i_inchi - 1;
                let mut n = 0;
                let old: Vec<u16> = prep[m].n_old_comp_number.clone().unwrap_or_default();
                for j in 0..prep[m].num_components.max(0) as usize {
                    let has = res.p_inchi[m]
                        .get(j)
                        .is_some_and(|p| p[0].is_some() || p[1].is_some());
                    if i + 1 == old[j] as usize && has {
                        n += 1;
                        if n == 1 {
                            let (pi, pa) = (res.p_inchi[m][j].clone(), res.p_aux[m][j].clone());
                            for k in 0..2 {
                                if let Some(c) = pi[k].as_ref().filter(|c| c.n_number_of_atoms > 0)
                                {
                                    if k == 0 {
                                        sd.num_non_taut[i_inchi] += 1;
                                    } else if c.len_tautomer > 0 {
                                        sd.num_taut[i_inchi] += 1;
                                    } else if !pi[0]
                                        .as_ref()
                                        .is_some_and(|c0| c0.n_number_of_atoms != 0)
                                    {
                                        sd.num_non_taut[i_inchi] += 1;
                                    }
                                }
                            }
                            res.p_inchi[i_inchi][i] = pi;
                            res.p_aux[i_inchi][i] = pa;
                        }
                    }
                }
                if n == 1 {
                    continue;
                }
                if n > 1 {
                    add_error_message(&mut sd.str_err_struct, "Cannot distinguish components");
                    sd.n_struct_read_error = 99;
                    sd.n_error_type = IS_ERROR;
                    return n_ret;
                }
            }
            n_ret = get_one_component(sd, log, &mut inp_cur_data, &prep[i_inchi], i, num_inp);
            if n_ret == IS_ERROR || n_ret == IS_FATAL {
                break 'cycle;
            }
            let (mut pi, mut pa) = ([None, None], [None, None]);
            n_ret = create_one_component_inchi(
                sd,
                ip,
                &mut inp_cur_data,
                &mut pi,
                &mut pa,
                i_inchi,
                i,
                &mut inp_norm_data,
                nc_flags,
            );
            res.p_inchi[i_inchi][i] = pi;
            res.p_aux[i_inchi][i] = pa;
            if n_ret != 0 {
                n_ret = treat_errors_in_create_one_component_inchi(sd, log, i, num_inp);
                break;
            }
        }
    }
    n_ret
}

/// `PrepareSaveOptBits` (`SAVE_OPT_*` of `ichi.h`; the bits above 0xFF are
/// lost in the `unsigned char`, as in C)
pub fn prepare_save_opt_bits(ip: &InputParms) -> u8 {
    use super::ichiparm::{
        INCHI_OUT_SAVEOPT, REQ_MODE_DIFF_UU_STEREO, REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU,
        REQ_MODE_STEREO, TG_FLAG_1_5_TAUT, TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_PT_06_00,
        TG_FLAG_PT_13_00, TG_FLAG_PT_16_00, TG_FLAG_PT_18_00, TG_FLAG_PT_22_00, TG_FLAG_PT_39_00,
        TG_FLAG_RECONNECT_COORD,
    };
    let mut bits: u8 = 0;
    if ip.b_inchi_output_options & INCHI_OUT_SAVEOPT == 0 {
        return bits;
    }
    let tf = |f: u64| ip.b_taut_flags & f != 0;
    if tf(TG_FLAG_RECONNECT_COORD) {
        bits |= 0x0008;
    }
    if ip.n_mode & REQ_MODE_BASIC != 0 {
        bits |= 0x0004;
    }
    if ip.n_mode & REQ_MODE_DIFF_UU_STEREO != 0 {
        bits |= 0x0001;
    }
    if ip.n_mode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) == 0 {
        bits |= 0x0002;
    }
    if tf(TG_FLAG_KETO_ENOL_TAUT) {
        bits |= 0x0010;
    }
    if tf(TG_FLAG_1_5_TAUT) {
        bits |= 0x0020;
    }
    if tf(TG_FLAG_PT_22_00) {
        bits |= 0x0040;
    }
    if tf(TG_FLAG_PT_16_00) {
        bits |= 0x0080;
    }
    // SAVE_OPT_PT_06_00 .. SAVE_OPT_PT_18_00 do not fit into the unsigned char
    let _ = (
        TG_FLAG_PT_06_00,
        TG_FLAG_PT_39_00,
        TG_FLAG_PT_13_00,
        TG_FLAG_PT_18_00,
    );
    if ip.n_mode & REQ_MODE_STEREO == 0 {
        bits &= !0x0002;
        bits &= !0x0001;
    }
    bits
}

/// `ProcessOneStructure` for a structure without polymer data: creates the
/// `INChI` of the (disconnected and reconnected) structure and prints the
/// InChI string and the AuxInfo into `out_file`. The Molfile output
/// (`OrigAtData_SaveMolfile`), the displays and the problem-file copies of
/// the C code are not modelled.
#[allow(clippy::too_many_arguments)]
pub fn process_one_structure(
    sd: &mut StructData,
    ip: &InputParms,
    log: &mut InchiIostream,
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    orig: &mut OrigAtomData,
    prep: &mut [OrigAtomData; 2],
    num_inp: i64,
) -> i32 {
    process_one_structure_arrays(sd, ip, log, out_file, strbuf, orig, prep, num_inp, None)
}

/// the `INChI` and `INChI_Aux` arrays returned by `ProcessOneStructure`
/// (`pINChI[iINChI][component][bTaut]`)
pub type InchiArrays = ([PInChI2; 2], [PInChIAux2; 2]);

/// `ProcessOneStructure` that also hands out the `INChI` arrays (the
/// `pINChI`, `pINChI_Aux` arguments of C)
#[allow(clippy::too_many_arguments)]
pub fn process_one_structure_arrays(
    sd: &mut StructData,
    ip: &InputParms,
    log: &mut InchiIostream,
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    orig: &mut OrigAtomData,
    prep: &mut [OrigAtomData; 2],
    num_inp: i64,
    out_arrays: Option<&mut InchiArrays>,
) -> i32 {
    use super::ichierr::add_error_message;
    use super::ichiparm::{
        REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO, REQ_MODE_STEREO, TG_FLAG_RECONNECT_COORD,
    };
    use super::ichiprt1::{FLAG_INP_AT_CHIRAL, OrigStruct, OutputInput, orig_struct_fill_out};
    use super::runichi3::{INCHI_REC, IS_FATAL, TG_FLAG_DISCONNECT_COORD_DONE};
    use super::runichi4::{
        SortPrintParams, is_struct_chiral, sort_and_print_inchi, treat_create_inchi_warning,
    };

    let mut n_ret = 0;
    sd.b_user_quit_component = 0;
    let mut res = StructResults::default();
    let mut nc = NormCanonFlags::default();

    // OrigAtData_StoreNativeInput
    let mut orig_struct = OrigStruct::default();
    match orig_struct_fill_out(
        orig.num_inp_atoms,
        orig.n_zy,
        orig.sz_coord.as_deref(),
        &orig.at,
        sd.b_chiral_flag,
        &mut sd.str_err_struct,
    ) {
        Ok(o) => {
            orig_struct = o;
            orig.sz_coord = None;
        }
        Err(_) => {
            add_error_message(
                &mut sd.str_err_struct,
                "Cannot interpret reversibility information",
            );
            sd.n_struct_read_error = 99;
            sd.n_error_type = IS_ERROR;
            n_ret = IS_ERROR;
        }
    }
    let ok = |n: i32| n != IS_FATAL && n != IS_ERROR;

    // 2. INChI of the whole disconnected or original structure
    if ok(n_ret) {
        let n1 = create_one_structure_inchi(
            sd, ip, log, num_inp, INCHI_BAS, orig, prep, &mut res, &mut nc,
        );
        n_ret = n_ret.max(n1);
    }
    // 3. INChI of the whole metal-reconnected structure
    if ok(n_ret)
        && sd.taut_flags_done[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE != 0
        && ip.b_taut_flags & TG_FLAG_RECONNECT_COORD != 0
    {
        let n1 = create_one_structure_inchi(
            sd, ip, log, num_inp, INCHI_REC, orig, prep, &mut res, &mut nc,
        );
        n_ret = n_ret.max(n1);
    }
    let mut pinchi: [PInChI2; 2] = [
        std::mem::take(&mut res.p_inchi[0]),
        std::mem::take(&mut res.p_inchi[1]),
    ];
    let mut pinchi_aux: [PInChIAux2; 2] = [
        std::mem::take(&mut res.p_aux[0]),
        std::mem::take(&mut res.p_aux[1]),
    ];
    let num_components = [sd.num_components[0], sd.num_components[1]];
    if ok(n_ret) {
        if sd.b_chiral_flag & FLAG_INP_AT_CHIRAL != 0
            && ip.n_mode & REQ_MODE_STEREO != 0
            && ip.n_mode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO) == 0
            && !is_struct_chiral(&pinchi, num_components)
            && ip.b_no_warnings == 0
        {
            add_error_message(&mut sd.str_err_struct, "Not chiral");
        }
        if sd.b_user_quit_component == 0 {
            let n1 = treat_create_inchi_warning(sd, log, num_inp);
            n_ret = n_ret.max(n1);
        }
    }

    // 4. sort and print
    let save_opt_bits = prepare_save_opt_bits(ip);
    if ok(n_ret) {
        let params = SortPrintParams {
            inchi_output_options: ip.b_inchi_output_options,
            compare_components: ip.b_compare_components,
            input: OutputInput {
                abc_numbers: ip.b_abc_numbers != 0,
                ct_predecessors: ip.b_ct_predecessors != 0,
                polymers: ip.b_polymers,
                large_molecules: ip.b_large_molecules != 0,
                sdf_label: String::new(),
                sdf_value: String::new(),
                orig_present: true,
                n_zy: orig.n_zy,
                n_pzz: 0,
            },
        };
        let mut flags = 0;
        let mut taut_flags = [sd.taut_flags[0], sd.taut_flags[1]];
        let mut taut_flags_done = [sd.taut_flags_done[0], sd.taut_flags_done[1]];
        let nc_tf = [nc.b_taut_flags[0], nc.b_taut_flags[1]];
        let nc_tfd = [nc.b_taut_flags_done[0], nc.b_taut_flags_done[1]];
        sort_and_print_inchi(
            out_file,
            strbuf,
            Some(log),
            &params,
            Some(&orig_struct),
            num_components,
            &mut taut_flags,
            &mut taut_flags_done,
            nc_tf,
            nc_tfd,
            num_inp as i32,
            &mut pinchi,
            &mut pinchi_aux,
            &mut flags,
            save_opt_bits,
        );
        sd.taut_flags[0] = taut_flags[0];
        sd.taut_flags[1] = taut_flags[1];
        sd.taut_flags_done[0] = taut_flags_done[0];
        sd.taut_flags_done[1] = taut_flags_done[1];
    }
    if let Some(o) = out_arrays {
        *o = (pinchi, pinchi_aux);
    }
    n_ret
}

/// `ProcessOneStructureExCore` + `ProcessOneStructureEx` for structures
/// without polymer data (the CRU frame-shift and folding code is not ported)
#[allow(clippy::too_many_arguments)]
pub fn process_one_structure_ex(
    sd: &mut StructData,
    ip: &InputParms,
    log: &mut InchiIostream,
    out_file: &mut InchiIostream,
    strbuf: &mut InchiIosString,
    orig: &mut OrigAtomData,
    prep: &mut [OrigAtomData; 2],
    num_inp: i64,
) -> i32 {
    use super::ichiparm::{INCHI_OUT_INCHI_GEN_ERROR, INCHI_OUT_STDINCHI};
    use super::runichi3::oad_validate_polymer_and_pseudo_element_data;
    // PreprocessPolymerCRUData
    orig.valid_polymer = 0;
    if orig.polymer.is_some() {
        // the CRU frame-shift and folding code is not ported
        super::ichierr::add_error_message(
            &mut sd.str_err_struct,
            "Polymer processing is not supported yet",
        );
        sd.n_error_type = IS_ERROR;
        return IS_ERROR;
    }

    // ProcessOneStructureExCore
    let mut res = oad_validate_polymer_and_pseudo_element_data(
        orig,
        ip.b_polymers,
        ip.b_np_zz,
        &mut sd.str_err_struct,
    );
    if res != 0 {
        sd.n_error_code = res;
        inchi_ios_eprint(
            Some(log),
            &format!(
                "Error {} ({}) structure #{}.\n",
                sd.n_error_code, sd.str_err_struct, num_inp
            ),
        );
        res = IS_ERROR;
        orig.num_inp_atoms = 0;
    } else {
        res = process_one_structure(sd, ip, log, out_file, strbuf, orig, prep, num_inp);
    }

    if ip.b_inchi_output_options2 & INCHI_OUT_INCHI_GEN_ERROR != 0
        && (res == IS_FATAL || res == IS_ERROR)
    {
        let t = if ip.b_inchi_output_options & INCHI_OUT_STDINCHI != 0 {
            "InChI=1S//\n"
        } else {
            "InChI=1//\n"
        };
        inchi_ios_eprint(Some(out_file), t);
    }
    res
}
