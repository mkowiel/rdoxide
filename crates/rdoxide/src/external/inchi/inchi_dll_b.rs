//! The Molfile text entry point of `INCHI_API/libinchi/src/inchi_dll_b.c`
//! (`MakeINCHIFromMolfileText`, `PrepareToMakeINCHI`) with the output
//! assembly of `inchi_dll.c` (`produce_generation_output`,
//! `copy_corrected_log_tail`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! Structures with polymer data (Molfile S-groups of the polymer types) are
//! not supported yet: the CRU frame-shift and folding code of the C library
//! is not ported.

use super::ichi_io::{
    INCHI_IOS_TYPE_STRING, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT, InchiIosString,
    InchiIostream, IosFile, inchi_ios_init, inchi_strbuf_init,
};
use super::ichierr::add_error_message;
use super::ichiparm::{
    INCHI_MAX_NUM_ARG, INCHI_OUT_INCHI_GEN_ERROR, INCHI_OUT_SDFILE_ONLY, INCHI_OUT_STDINCHI,
    INPUT_MOLFILE, InputParms, parse_options_string, read_command_line_parms,
};
use super::runichi::process_one_structure_ex;
use super::runichi2::get_one_structure;
use super::runichi3::{IS_EOF, IS_ERROR, IS_FATAL, IS_SKIP, IS_WARNING, StructData};
use super::strutil::OrigAtomData;

/// `mol2inchi_Ret_*`
pub const MOL2INCHI_RET_OKAY: i32 = 0;
pub const MOL2INCHI_RET_WARNING: i32 = 1;
pub const MOL2INCHI_RET_ERROR: i32 = 2;
/// `MOL2INCHI_BAD_COMMAND_LINE`
pub const MOL2INCHI_BAD_COMMAND_LINE: i32 = 1002;

/// `inchi_Output`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InchiOutput {
    pub sz_inchi: Option<String>,
    pub sz_aux_info: Option<String>,
    pub sz_message: Option<String>,
    pub sz_log: Option<String>,
}

fn c_prefix(b: &[u8]) -> &[u8] {
    &b[..b.iter().position(|&c| c == 0).unwrap_or(b.len())]
}

fn text(b: &[u8]) -> String {
    String::from_utf8_lossy(c_prefix(b)).into_owned()
}

/// `copy_corrected_log_tail`: strips the trailing LFs of the log and cuts it
/// at the first ` structure #`
fn copy_corrected_log_tail(out: &mut InchiOutput, log: &mut InchiIostream) {
    let used = log.s.n_used_length.max(0) as usize;
    let Some(p) = log.s.p_str.as_mut() else {
        return;
    };
    if used == 0 {
        return;
    }
    let mut used = used.min(p.len());
    while used > 0 && p[used - 1] == b'\n' {
        used -= 1;
        p[used] = 0;
    }
    log.s.n_used_length = used as i32;
    let mut s = c_prefix(p).to_vec();
    if let Some(k) = s.windows(12).position(|w| w == b" structure #") {
        s.truncate(k);
    }
    out.sz_log = Some(String::from_utf8_lossy(&s).into_owned());
    log.s.p_str = None;
}

/// `produce_generation_output`
fn produce_generation_output(
    out: &mut InchiOutput,
    sd: &StructData,
    ip: &InputParms,
    log: &mut InchiIostream,
    out_file: &mut InchiIostream,
) {
    if !sd.str_err_struct.is_empty() {
        out.sz_message = Some(sd.str_err_struct.clone());
    }
    if out_file.s.n_used_length > 0 && out_file.s.p_str.is_some() {
        let mut b = out_file.s.p_str.take().unwrap_or_default();
        b.push(0);
        let at = |b: &[u8], k: usize| b.get(k).copied().unwrap_or(0);
        let mut aux: Option<usize> = None;
        if ip.b_inchi_output_options & INCHI_OUT_SDFILE_ONLY == 0 {
            let find = |b: &[u8], from: usize| {
                (from..b.len())
                    .take_while(|&k| b[k] != 0 || k < from)
                    .find(|&k| b[k] == b'\n')
            };
            let mut p = find(&b, 0);
            while let Some(pp) = p {
                if (0..8).all(|k| at(&b, pp + k) == b"\nAuxInfo"[k]) {
                    b[pp] = 0;
                    aux = Some(pp + 1);
                } else if aux.is_some() || at(&b, pp + 1) == 0 {
                    b[pp] = 0;
                    break;
                }
                p = find(&b, pp + 1);
            }
        }
        out.sz_inchi = Some(text(&b));
        out.sz_aux_info = aux.map(|k| text(&b[k..]));
    }
    copy_corrected_log_tail(out, log);
}

/// `MakeINCHIFromMolfileText`: returns the `mol2inchi_Ret_*` code and the
/// output strings; `options` is the option string (`None` for NULL)
pub fn make_inchi_from_molfile_text(moltext: &[u8], options: Option<&[u8]>) -> (i32, InchiOutput) {
    let mut result = InchiOutput::default();
    let mut sd = StructData::default();
    let mut ip = InputParms::default();
    let mut orig = OrigAtomData::default();
    let mut prep = [OrigAtomData::default(), OrigAtomData::default()];
    let mut pout = InchiIostream::default();
    let mut plog = InchiIostream::default();
    let mut inp_file = InchiIostream::default();
    let mut strbuf = InchiIosString::default();
    let mut num_inp: i64 = 0;
    let mut retcode;

    // PrepareToMakeINCHI
    let argv: Vec<Vec<u8>> = match options {
        Some(o) => parse_options_string(o, INCHI_MAX_NUM_ARG),
        None => vec![Vec::new()],
    };
    inchi_ios_init(&mut pout, INCHI_IOS_TYPE_STRING, IosFile::None);
    inchi_ios_init(&mut plog, INCHI_IOS_TYPE_STRING, IosFile::None);
    inchi_ios_init(&mut inp_file, INCHI_IOS_TYPE_STRING, IosFile::None);
    let mut mol = c_prefix(moltext).to_vec();
    mol.push(0);
    inp_file.s.n_used_length = mol.len() as i32;
    inp_file.s.n_allocated_length = mol.len() as i32;
    inp_file.s.p_str = Some(mol);
    inp_file.s.n_ptr = 0;
    let argv_str: Vec<String> = argv
        .iter()
        .map(|a| String::from_utf8_lossy(a).into_owned())
        .collect();
    let argv_ref: Vec<&str> = argv_str.iter().map(String::as_str).collect();
    let mut prepare = 0;
    if read_command_line_parms(&argv_ref, &mut ip, true, &mut plog) < 0 {
        prepare = MOL2INCHI_BAD_COMMAND_LINE;
    } else {
        ip.n_input_type = INPUT_MOLFILE;
        ip.b_no_struct_labels = 1;
        inchi_strbuf_init(
            &mut strbuf,
            INCHI_STRBUF_INITIAL_SIZE,
            INCHI_STRBUF_SIZE_INCREMENT,
        );
    }
    let output_error_inchi = ip.b_inchi_output_options2 & INCHI_OUT_INCHI_GEN_ERROR != 0;

    'main: {
        if prepare != 0 {
            if let Some(p) = plog.s.p_str.as_mut() {
                let mut t = text(p);
                add_error_message(&mut t, "Error while preparing to make InChI");
                let b = t.as_bytes();
                if p.len() <= b.len() {
                    p.resize(b.len() + 1, 0);
                }
                p[..b.len()].copy_from_slice(b);
                p[b.len()] = 0;
            }
            retcode = MOL2INCHI_RET_ERROR;
            break 'main;
        }
        retcode = get_one_structure(
            &mut sd,
            &mut ip,
            &mut plog,
            &mut inp_file,
            &mut orig,
            &mut num_inp,
        );
        if retcode == IS_FATAL || retcode == IS_ERROR {
            retcode = MOL2INCHI_RET_ERROR;
            if !output_error_inchi {
                break 'main;
            }
        } else if retcode == IS_EOF || retcode == IS_SKIP {
            break 'main;
        }
        // always enable polymer extensions
        if orig.polymer.is_some() {
            orig.valid_polymer = 1;
        }
        retcode = process_one_structure_ex(
            &mut sd,
            &ip,
            &mut plog,
            &mut pout,
            &mut strbuf,
            &mut orig,
            &mut prep,
            num_inp,
        );
        if retcode != IS_FATAL && retcode != IS_ERROR {
            produce_generation_output(&mut result, &sd, &ip, &mut plog, &mut pout);
            if retcode == IS_WARNING {
                retcode = MOL2INCHI_RET_WARNING;
            }
        } else {
            retcode = MOL2INCHI_RET_ERROR;
            if output_error_inchi {
                produce_generation_output(&mut result, &sd, &ip, &mut plog, &mut pout);
                if result.sz_inchi.is_none() {
                    result.sz_inchi = Some(
                        if ip.b_inchi_output_options & INCHI_OUT_STDINCHI != 0 {
                            "InChI=1S//"
                        } else {
                            "InChI=1//"
                        }
                        .to_string(),
                    );
                }
            }
        }
    }
    copy_corrected_log_tail(&mut result, &mut plog);
    (retcode.max(0), result)
}
