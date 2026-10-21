//! Preprocessing of one structure (`INCHI_BASE/src/runichi3.c`:
//! `PreprocessOneStructure`, `OrigAtData_Duplicate`,
//! `OrigAtData_bCheckUnusualValences`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! Build facts: `bRELEASE_VERSION` is 1 (no `bExtract` bookkeeping),
//! `FIX_ADJ_RAD` is 0 (`FixAdjacentRadicals` is not called) and
//! `REMOVE_ION_PAIRS_DISC_STRU` is 1. Polymer, V3000 and coordinate-string
//! members of `ORIG_ATOM_DATA` are not modelled.

use super::ichierr::add_error_message;
use super::ichiparm::{
    InputParms, TG_FLAG_CHECK_VALENCE_COORD, TG_FLAG_DISCONNECT_COORD, TG_FLAG_DISCONNECT_SALTS,
    TG_FLAG_FIX_SP3_BUG, TG_FLAG_RECONNECT_COORD,
};
use super::ichister::reconcile_all_cml_bond_parities;
use super::incomdef::{RADICAL_DOUBLET, RADICAL_SINGLET, RADICAL_TRIPLET};
use super::inpdef::{MAX_NUM_STEREO_BONDS, SB_PARITY_FLAG, SB_PARITY_MASK, SB_PARITY_SHFT};
use super::strutil::{
    OrigAtomData, disconnect_metals, disconnect_salts, fix_odd_things,
    mark_disconnected_components, may_disconnect_metals, num_hetero_atom_has_isotopic_h,
    post_fix_odd_things, remove_ion_pairs,
};
use super::util::detect_unusual_el_valence;

pub const IS_OKAY: i32 = 0;
pub const IS_WARNING: i32 = 1;
pub const IS_ERROR: i32 = 2;
pub const IS_FATAL: i32 = 3;
pub const IS_EOF: i32 = -1;
pub const IS_SKIP: i32 = -2;

pub const INCHI_BAS: usize = 0;
pub const INCHI_REC: usize = 1;

pub const TG_FLAG_DISCONNECT_SALTS_DONE: u64 = 0x0000_0004;
pub const TG_FLAG_DISCONNECT_COORD_DONE: u64 = 0x0000_0100;
pub const TG_FLAG_FIX_ODD_THINGS_DONE: u64 = 0x0000_0800;
pub const TG_FLAG_FOUND_ISOTOPIC_H_DONE: u64 = 0x0000_4000;
pub const TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE: u64 = 0x0000_8000;

/// The members of `STRUCT_DATA` used by the preprocessing
#[derive(Debug, Clone, Default)]
pub struct StructData {
    pub n_error_type: i32,
    pub n_struct_read_error: i32,
    pub str_err_struct: String,
    pub taut_flags: [u64; 2],
    pub taut_flags_done: [u64; 2],
    pub n_error_code: i32,
    pub num_taut: [i32; 3],
    pub num_non_taut: [i32; 3],
    pub num_components: [i32; 3],
    pub b_user_quit_component: i32,
    /// `bChiralFlag`
    pub b_chiral_flag: i32,
}

/// `OrigAtData_bCheckUnusualValences`: counts atoms with an unusual valence
/// and reports them in `str_err_struct`
pub fn orig_at_data_check_unusual_valences(
    orig: &OrigAtomData,
    b_add_iso_h: bool,
    str_err_struct: &mut String,
    b_no_warnings: bool,
) -> i32 {
    let mut num_found = 0;
    if orig.num_inp_atoms == 0 {
        return 0;
    }
    for a in &orig.at[..orig.num_inp_atoms] {
        let num_h = if b_add_iso_h {
            a.numh()
        } else {
            i32::from(a.num_h)
        };
        let val = detect_unusual_el_valence(
            i32::from(a.el_number),
            i32::from(a.charge),
            i32::from(a.radical),
            i32::from(a.chem_bonds_valence),
            num_h,
            i32::from(a.valence),
        );
        if val != 0 {
            num_found += 1;
            if !b_no_warnings {
                add_error_message(str_err_struct, "Accepted unusual valence(s):");
            }
            let end = a
                .elname
                .iter()
                .position(|&c| c == 0)
                .unwrap_or(a.elname.len());
            let mut msg = String::from_utf8_lossy(&a.elname[..end]).into_owned();
            if a.charge != 0 {
                msg += &format!("{:+}", a.charge);
            }
            if a.radical != 0 {
                let r = i32::from(a.radical);
                msg += &format!(
                    ",{}",
                    if r == RADICAL_SINGLET {
                        "s"
                    } else if r == RADICAL_DOUBLET {
                        "d"
                    } else if r == RADICAL_TRIPLET {
                        "t"
                    } else {
                        "?"
                    }
                );
            }
            msg += &format!("({val})");
            if !b_no_warnings {
                add_error_message(str_err_struct, &msg);
            }
        }
    }
    num_found
}

/// `OrigAtData_Duplicate` (into an empty `new_orig_atom`): copies the atoms
/// and the component arrays; the equivalence data are not copied
pub fn orig_at_data_duplicate(new_orig_atom: &mut OrigAtomData, orig_atom: &OrigAtomData) -> i32 {
    let n = orig_atom.num_inp_atoms;
    let nc = orig_atom.num_components.max(0) as usize;
    let mut at = orig_atom.at[..n.min(orig_atom.at.len())].to_vec();
    at.resize(n + 1, Default::default());
    let copy_len = |src: &Option<Vec<u16>>| {
        let mut v = vec![0u16; nc + 1];
        if let Some(s) = src {
            let k = nc.min(s.len());
            v[..k].copy_from_slice(&s[..k]);
        }
        v
    };
    let cur = copy_len(&orig_atom.n_cur_at_len);
    let old = copy_len(&orig_atom.n_old_comp_number);
    *new_orig_atom = orig_atom.clone();
    new_orig_atom.at = at;
    new_orig_atom.n_cur_at_len = Some(cur);
    new_orig_atom.n_old_comp_number = Some(old);
    0
}

fn strip_disconnected_parities(od: &mut OrigAtomData) {
    let n = od.num_inp_atoms;
    for a in &mut od.at[..n] {
        let mut m = 0;
        while m < MAX_NUM_STEREO_BONDS && a.sb_parity[m] != 0 {
            a.sb_parity[m] &= SB_PARITY_MASK as i8;
            m += 1;
        }
    }
}

/// `PreprocessOneStructure`: copies `orig` to `prep[0]`, fixes odd things,
/// disconnects salts and metals, marks the components. `prep[1]` receives
/// the reconnected structure if requested. Returns `sd.n_error_type`.
pub fn preprocess_one_structure(
    sd: &mut StructData,
    ip: &InputParms,
    orig: &mut OrigAtomData,
    prep: &mut [OrigAtomData; 2],
) -> i32 {
    let b_taut_flags = 0u64;
    let mut b_taut_flags_done = 0u64;
    let no_warn = ip.b_no_warnings != 0;
    let fatal_ram = |sd: &mut StructData| {
        add_error_message(&mut sd.str_err_struct, "Out of RAM");
        sd.n_struct_read_error = 99;
        sd.n_error_type = IS_FATAL;
    };

    'body: {
        // 1. copy
        if orig_at_data_duplicate(&mut prep[0], orig) < 0 {
            fatal_ram(sd);
            break 'body;
        }
        // 2. fix odd things
        let n = prep[0].num_inp_atoms;
        if fix_odd_things(
            &mut prep[0].at[..n],
            ip.b_taut_flags & TG_FLAG_FIX_SP3_BUG != 0,
            ip.b_fix_non_uniform_draw != 0,
        ) > 0
        {
            if !no_warn {
                add_error_message(&mut sd.str_err_struct, "Charges were rearranged");
            }
            if sd.n_error_type < IS_WARNING {
                sd.n_error_type = IS_WARNING;
            }
            sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
        }

        // needs salt disconnection?
        if ip.b_taut_flags & TG_FLAG_DISCONNECT_SALTS != 0 {
            let n = prep[0].num_inp_atoms;
            let mut nb = prep[0].num_inp_bonds;
            prep[0].b_disconnect_salts =
                i32::from(disconnect_salts(&mut prep[0].at[..n], &mut nb, false) > 0);
        } else {
            prep[0].b_disconnect_salts = 0;
        }
        // needs metal disconnection?
        if ip.b_taut_flags & TG_FLAG_DISCONNECT_COORD != 0 {
            let n = prep[0].num_inp_atoms;
            let (_, bdc) = may_disconnect_metals(&prep[0].at[..n]);
            prep[0].b_disconnect_coord = bdc;
            sd.taut_flags_done[INCHI_BAS] |= b_taut_flags_done;
        } else {
            prep[0].b_disconnect_coord = 0;
        }
        orig.b_disconnect_salts = prep[0].b_disconnect_salts;
        orig.b_disconnect_coord = prep[0].b_disconnect_coord;

        // 3. disconnect salts
        let mut salts_done = false;
        if ip.b_taut_flags & TG_FLAG_DISCONNECT_SALTS != 0 && prep[0].b_disconnect_salts != 0 {
            let n = prep[0].num_inp_atoms;
            let mut nb = prep[0].num_inp_bonds;
            let r = disconnect_salts(&mut prep[0].at[..n], &mut nb, true);
            prep[0].num_inp_bonds = nb;
            if r > 0 {
                salts_done = true;
                if !no_warn {
                    add_error_message(&mut sd.str_err_struct, "Salt was disconnected");
                }
                sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_DISCONNECT_SALTS_DONE;
                if sd.n_error_type < IS_WARNING {
                    sd.n_error_type = IS_WARNING;
                }
                let i = reconcile_all_cml_bond_parities(&mut prep[0].at[..n], false);
                if i != 0 {
                    add_error_message(&mut sd.str_err_struct, "0D Parities Reconciliation failed:");
                    add_error_message(&mut sd.str_err_struct, &i.to_string());
                }
            }
        }
        if !salts_done {
            prep[0].b_disconnect_salts = 0;
        }

        // mark the (disconnected) components
        prep[0].num_components = mark_disconnected_components(&mut prep[0], false);
        if prep[0].num_components < 0 {
            fatal_ram(sd);
            break 'body;
        }

        // isotopic H on heteroatoms
        let n = prep[0].num_inp_atoms;
        let i = num_hetero_atom_has_isotopic_h(&prep[0].at[..n]);
        if i & 1 != 0 {
            sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_H_DONE;
        }
        if i & 2 != 0 {
            sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE;
        }

        // 4a. unusual valences
        orig_at_data_check_unusual_valences(&prep[0], true, &mut sd.str_err_struct, no_warn);

        // 5. metals
        if prep[0].b_disconnect_coord != 0 {
            prep[0].num_components = mark_disconnected_components(&mut prep[0], false);
            if prep[0].num_components < 0 {
                fatal_ram(sd);
                break 'body;
            }
            if ip.b_taut_flags & TG_FLAG_RECONNECT_COORD != 0 {
                let (first, second) = prep.split_at_mut(1);
                if orig_at_data_duplicate(&mut second[0], &first[0]) < 0 {
                    fatal_ram(sd);
                    break 'body;
                }
                sd.taut_flags[INCHI_REC] = sd.taut_flags[INCHI_BAS];
                sd.taut_flags_done[INCHI_REC] = sd.taut_flags_done[INCHI_BAS];
                strip_disconnected_parities(&mut prep[1]);
            }
            let check = ip.b_taut_flags & TG_FLAG_CHECK_VALENCE_COORD != 0;
            let i = disconnect_metals(&mut prep[0], check, Some(&mut b_taut_flags_done));
            if i > 0 {
                if !no_warn {
                    add_error_message(&mut sd.str_err_struct, "Metal was disconnected");
                }
                sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_DISCONNECT_COORD_DONE;
                if sd.n_error_type < IS_WARNING {
                    sd.n_error_type = IS_WARNING;
                }
                prep[0].num_components = mark_disconnected_components(&mut prep[0], true);
                if prep[0].num_components < 0 {
                    fatal_ram(sd);
                    break 'body;
                }
                let n = prep[0].num_inp_atoms;
                for a in &mut prep[0].at[..n] {
                    let mut m = 0;
                    while m < MAX_NUM_STEREO_BONDS && a.sb_parity[m] != 0 {
                        let p = a.sb_parity[m];
                        if p & SB_PARITY_FLAG != 0 {
                            a.sb_parity[m] =
                                ((i32::from(p) >> SB_PARITY_SHFT) & SB_PARITY_MASK) as i8;
                        }
                        m += 1;
                    }
                }
                let i = reconcile_all_cml_bond_parities(&mut prep[0].at[..n], true);
                if i != 0 {
                    add_error_message(&mut sd.str_err_struct, "0D Parities Reconciliation failed:");
                    add_error_message(&mut sd.str_err_struct, &i.to_string());
                }
                // REMOVE_ION_PAIRS_DISC_STRU == 1
                if remove_ion_pairs(&mut prep[0].at[..n]) > 0 {
                    if !no_warn {
                        add_error_message(&mut sd.str_err_struct, "Charges were rearranged");
                    }
                    if sd.n_error_type < IS_WARNING {
                        sd.n_error_type = IS_WARNING;
                    }
                    sd.taut_flags_done[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
                    sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
                }
            } else if i < 0 {
                add_error_message(&mut sd.str_err_struct, "Cannot disconnect metal error");
                sd.n_struct_read_error = i;
                sd.n_error_type = IS_ERROR;
                break 'body;
            }
        } else {
            strip_disconnected_parities(&mut prep[0]);
        }
    }

    // exit_function
    if sd.n_error_type < IS_ERROR {
        let n = prep[0].num_inp_atoms;
        if prep[0].at.len() >= n && post_fix_odd_things(&mut prep[0].at[..n]) > 0 {
            if !no_warn {
                add_error_message(&mut sd.str_err_struct, "Charges were rearranged");
            }
            if sd.n_error_type < IS_WARNING {
                sd.n_error_type = IS_WARNING;
            }
            sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
        }
        if sd.taut_flags_done[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE != 0
            && !prep[1].at.is_empty()
            && prep[1].num_inp_atoms > 0
        {
            let n = prep[1].num_inp_atoms;
            if post_fix_odd_things(&mut prep[1].at[..n]) > 0 {
                if !no_warn {
                    add_error_message(&mut sd.str_err_struct, "Charges were rearranged");
                }
                if sd.n_error_type < IS_WARNING {
                    sd.n_error_type = IS_WARNING;
                }
                sd.taut_flags_done[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
                sd.taut_flags_done[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
            }
        }
    }
    sd.taut_flags[INCHI_BAS] |= b_taut_flags;
    sd.taut_flags_done[INCHI_BAS] |= b_taut_flags_done;
    sd.n_error_type
}

/// `OAD_ValidateAndSortOutPseudoElementAtoms` for structures without polymer
/// data: the pseudoelement atoms `*`, `Zz` and `Zy` become `Zy`; sets
/// `orig.n_zy`. Returns the error code (0 for none).
fn oad_validate_and_sort_out_pseudo_element_atoms(
    orig: &mut OrigAtomData,
    treat_polymers: i32,
    b_np_zz: i32,
    p_str_err: &mut String,
) -> i32 {
    use super::molfmt::treat_err;
    let mut err = 0;
    let mut n_pseudo = 0;
    let pseudos_allowed = b_np_zz == 1 || treat_polymers != super::ichiparm::POLYMERS_NO;
    let nat = orig.num_inp_atoms;
    for k in 0..nat {
        let name_len = orig.at[k]
            .elname
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(orig.at[k].elname.len());
        let name = orig.at[k].elname[..name_len].to_vec();
        let is_zy = name == b"Zy";
        let is_star = name == b"*";
        let is_zz = !is_star && name == b"Zz";
        if is_star || is_zz || is_zy {
            n_pseudo += 1;
            if !pseudos_allowed {
                treat_err(&mut err, 70 + 5, p_str_err, "Invalid element(s):");
                treat_err(&mut err, 70 + 5, p_str_err, &String::from_utf8_lossy(&name));
                continue;
            }
            if orig.at[k].valence > 1 || orig.at[k].chem_bonds_valence > 1 {
                treat_err(
                    &mut err,
                    70 + 7,
                    p_str_err,
                    "Invalid pseudo element(s) bonding",
                );
                continue;
            }
            orig.at[k].elname[0] = b'Z';
            orig.at[k].elname[1] = b'y';
            orig.at[k].elname[2] = 0;
        }
    }
    orig.n_zy = n_pseudo;
    if orig.n_zy != 0 && b_np_zz == 0 {
        treat_err(
            &mut err,
            70 + 4,
            p_str_err,
            "Polymer-unrelated pseudoatoms are not allowed",
        );
    }
    if err != 0 {
        orig.valid_polymer = 0;
    }
    err
}

/// `OAD_ValidatePolymerAndPseudoElementData` for structures without polymer
/// data (`orig.polymer` is `None`; the polymer unit checks are not ported).
/// Returns the error code (0 for none).
pub fn oad_validate_polymer_and_pseudo_element_data(
    orig: &mut OrigAtomData,
    treat_polymers: i32,
    b_np_zz: i32,
    p_str_err: &mut String,
) -> i32 {
    use super::molfmt::treat_err;
    orig.valid_polymer = 0;
    let mut err =
        oad_validate_and_sort_out_pseudo_element_atoms(orig, treat_polymers, b_np_zz, p_str_err);
    if err != 0 {
        return err;
    }
    if orig.n_zy > 0 && b_np_zz == 0 {
        treat_err(
            &mut err,
            9,
            p_str_err,
            "Non-polymer-related Zz/star atoms are not allowed",
        );
    }
    err
}
