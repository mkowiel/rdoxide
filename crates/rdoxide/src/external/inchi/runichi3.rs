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
use super::mol2atom::{OadPolymer, OadPolymerUnit};
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

pub const NO_POLYMER: i32 = -1;
pub const POLYMER_REPRESENTATION_SOURCE_BASED: i32 = 1;
pub const POLYMER_REPRESENTATION_STRUCTURE_BASED: i32 = 2;
pub const POLYMER_REPRESENTATION_MIXED: i32 = 3;
pub const POLYMER_REPRESENTATION_UNRECOGNIZED: i32 = 4;

pub const POLYMER_STY_NON: i32 = 0;
pub const POLYMER_STY_SRU: i32 = 1;
pub const POLYMER_STY_MON: i32 = 2;
pub const POLYMER_STY_COP: i32 = 3;
pub const POLYMER_STY_MOD: i32 = 4;
pub const POLYMER_STY_CRO: i32 = 5;
pub const POLYMER_STY_MER: i32 = 6;

pub const POLYMER_SST_NON: i32 = 0;
pub const POLYMER_SST_ALT: i32 = 1;
pub const POLYMER_SST_RAN: i32 = 2;
pub const POLYMER_SST_BLK: i32 = 3;

pub const POLYMER_CONN_NON: i32 = 0;
pub const POLYMER_CONN_HT: i32 = 1;
pub const POLYMER_CONN_HH: i32 = 2;
pub const POLYMER_CONN_EU: i32 = 3;

pub const CLOSING_SRU_NOT_APPLICABLE: i32 = 0;
pub const CLOSING_SRU_RING: i32 = 1;
pub const CLOSING_SRU_HIGHER_ORDER_BOND: i32 = 2;
pub const CLOSING_SRU_DIRADICAL: i32 = 3;
pub const CLOSING_SRU_RING_OPENED: i32 = 4;
pub const CLOSING_SRU_MOVED_BRACKETS: i32 = 5;

// inpdef.h: ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND, ALLOW_CLOSING_SRU_VIA_DIRADICAL,
// CLOSING_STARRED_SRU_IS_A_MUST and ALLOW_MIXED_SRU_AND_MON are all 1

fn elname_is(orig: &OrigAtomData, idx: usize, name: &str) -> bool {
    let e = &orig.at[idx].elname;
    let len = e.iter().position(|&c| c == 0).unwrap_or(e.len());
    &e[..len] == name.as_bytes()
}

fn set_elname(orig: &mut OrigAtomData, idx: usize, name: &str) {
    let e = &mut orig.at[idx].elname;
    e.iter_mut().for_each(|c| *c = 0);
    e[..name.len()].copy_from_slice(name.as_bytes());
}

/// `OAD_PolymerUnit_FindEndsAndCaps`: the aliased output parameters of the C
/// code are the fields of the unit itself. Returns the error code.
fn oad_polymer_unit_find_ends_and_caps(
    unit: &mut OadPolymerUnit,
    orig: &OrigAtomData,
    p_str_err: &mut String,
) -> i32 {
    use super::molfmt::treat_err;
    let num_atoms = orig.num_inp_atoms as i32;
    unit.end_atom1 = 0;
    unit.end_atom2 = 0;
    unit.cap1 = 0;
    unit.cap2 = 0;
    unit.cap1_is_undef = 0;
    unit.cap2_is_undef = 0;
    let mut err = 0;
    if unit.blist.is_empty() || unit.nb < 1 {
        return 0;
    }
    // left crossing bond
    let (i, j) = (unit.blist[0], unit.blist[1]);
    let i_inside = unit.alist.contains(&i);
    let j_inside = unit.alist.contains(&j);
    if i_inside && j_inside {
        treat_err(
            &mut err,
            9032,
            p_str_err,
            "Polymer CRU cap(s) lie inside CRU",
        );
        return err;
    }
    if i_inside {
        unit.end_atom1 = i;
        unit.cap1 = j;
    } else {
        unit.end_atom1 = j;
        unit.cap1 = i;
    }
    if unit.cap1 >= 1
        && (unit.cap1 as usize) <= orig.at.len()
        && elname_is(orig, unit.cap1 as usize - 1, "Zz")
    {
        unit.cap1_is_undef = 1;
    }
    // right crossing bond (an error here does not stop the C code)
    let (i, j) = (unit.blist[2], unit.blist[3]);
    let i_inside = unit.alist.contains(&i);
    let j_inside = unit.alist.contains(&j);
    if i_inside && j_inside {
        treat_err(
            &mut err,
            9032,
            p_str_err,
            "Polymer CRU cap(s) lie inside CRU",
        );
    }
    if i_inside {
        unit.end_atom2 = i;
        unit.cap2 = j;
    } else {
        unit.end_atom2 = j;
        unit.cap2 = i;
    }
    if unit.cap2 >= 1
        && (unit.cap2 as usize) <= orig.at.len()
        && elname_is(orig, unit.cap2 as usize - 1, "Zz")
    {
        unit.cap2_is_undef = 1;
    }
    let bad = |v: i32| v <= 0 || v > num_atoms;
    if bad(unit.end_atom1) || bad(unit.cap1) {
        treat_err(
            &mut err,
            9090,
            p_str_err,
            "Invalid polymer CRU crossing bond",
        );
        return err;
    }
    if bad(unit.end_atom2) || bad(unit.cap2) {
        treat_err(
            &mut err,
            9091,
            p_str_err,
            "Invalid polymer CRU crossing bond",
        );
        return err;
    }
    if unit.cap1 == unit.cap2 {
        treat_err(&mut err, 9090, p_str_err, "Invalid polymer CRU surrounding");
        return err;
    }
    // `*err = 0` at the end of the C function: a 9032 on the right bond is dropped
    0
}

/// `OAD_PolymerUnit_SetEndsAndCaps`: returns the error code.
fn oad_polymer_unit_set_ends_and_caps(
    unit: &mut OadPolymerUnit,
    orig: &OrigAtomData,
    p_str_err: &mut String,
) -> i32 {
    unit.cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    unit.end_atom1 = -1;
    unit.end_atom2 = -1;
    unit.cap1 = -1;
    unit.cap2 = -1;
    unit.cap1_is_undef = 0;
    unit.cap2_is_undef = 0;
    let err = oad_polymer_unit_find_ends_and_caps(unit, orig, p_str_err);
    if err != 0 {
        return err;
    }
    if unit.cap1_is_undef == 0 && unit.cap2_is_undef == 0 {
        return 0;
    }
    // the rest is applicable only to the *---SRU---* case
    if unit.end_atom1 == unit.end_atom2 {
        unit.cyclizable = CLOSING_SRU_DIRADICAL;
        return 0;
    }
    let a1 = &orig.at[unit.end_atom1 as usize - 1];
    for k in 0..a1.valence as usize {
        if i32::from(a1.neighbor[k]) == unit.end_atom2 - 1 {
            unit.cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
            return 0;
        }
    }
    unit.cyclizable = CLOSING_SRU_RING;
    0
}

/// `OAD_Polymer_GetRepresentation`
pub fn oad_polymer_get_representation(p: Option<&mut OadPolymer>) -> i32 {
    let Some(p) = p else {
        return NO_POLYMER;
    };
    let (mut n_source, mut n_struct) = (0, 0);
    for u in p.units.iter_mut().take(p.n.max(0) as usize) {
        if u.nb == 2 || u.nbkbonds > 0 || (u.cap1 > 0 && u.cap2 > 0) {
            u.representation = POLYMER_REPRESENTATION_STRUCTURE_BASED;
            n_struct += 1;
        } else if u.nb == 0 {
            u.representation = POLYMER_REPRESENTATION_SOURCE_BASED;
            n_source += 1;
        }
    }
    if p.n == n_source {
        POLYMER_REPRESENTATION_SOURCE_BASED
    } else if p.n == n_struct {
        POLYMER_REPRESENTATION_STRUCTURE_BASED
    } else if n_source > 0 && n_struct > 0 && n_source + n_struct == p.n {
        POLYMER_REPRESENTATION_MIXED
    } else {
        POLYMER_REPRESENTATION_UNRECOGNIZED
    }
}

/// `OAD_ValidateAndSortOutPseudoElementAtoms`: the pseudoelement atoms `*`,
/// `Zz` and `Zy` become `Zy`; for a valid polymer the caps of units with two
/// undefined-nature caps are turned back into `Zz`. Sets `orig.n_zy` and the
/// number of polymer `Zz` atoms; returns the error code (0 for none).
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
    for k in 0..orig.num_inp_atoms {
        let is_zy = elname_is(orig, k, "Zy");
        let is_star = elname_is(orig, k, "*");
        let is_zz = !is_star && elname_is(orig, k, "Zz");
        if is_star || is_zz || is_zy {
            n_pseudo += 1;
            if !pseudos_allowed {
                let name = {
                    let e = &orig.at[k].elname;
                    String::from_utf8_lossy(&e[..e.iter().position(|&c| c == 0).unwrap_or(e.len())])
                        .into_owned()
                };
                treat_err(&mut err, 70 + 5, p_str_err, "Invalid element(s):");
                treat_err(&mut err, 70 + 5, p_str_err, &name);
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
            set_elname(orig, k, "Zy");
        }
    }
    orig.n_zy = 0;
    let mut nzz = 0;
    if orig.valid_polymer != 0 {
        let caps: Vec<(i32, i32)> = orig
            .polymer
            .as_ref()
            .map(|pd| {
                pd.units
                    .iter()
                    .take(pd.n.max(0) as usize)
                    .filter(|u| u.cap1_is_undef != 0 && u.cap2_is_undef != 0)
                    .map(|u| (u.cap1, u.cap2))
                    .collect()
            })
            .unwrap_or_default();
        for (c1, c2) in caps {
            set_elname(orig, c1 as usize - 1, "Zz");
            set_elname(orig, c2 as usize - 1, "Zz");
            nzz += 2;
        }
        if let Some(pd) = orig.polymer.as_mut() {
            pd.n_pzz = nzz;
        }
    }
    orig.n_zy = n_pseudo - nzz;
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

/// `OAD_ValidatePolymerAndPseudoElementData`: checks the polymer units and
/// pseudoelement atoms and fixes what can be fixed. Returns the error code
/// (0 for none).
pub fn oad_validate_polymer_and_pseudo_element_data(
    orig: &mut OrigAtomData,
    treat_polymers: i32,
    b_np_zz: i32,
    p_str_err: &mut String,
    b_no_warnings: bool,
) -> i32 {
    use super::molfmt::treat_err;
    let nat = orig.num_inp_atoms as i32;
    let mut err = 0;
    orig.valid_polymer = 0;
    if treat_polymers != 0 && orig.polymer.is_some() {
        orig.valid_polymer = 1;
    }
    let nsgroups = if orig.valid_polymer != 0 {
        orig.polymer.as_ref().map_or(0, |p| p.n.max(0) as usize)
    } else {
        0
    };
    let warn = |m: &mut String, text: &str| {
        if !b_no_warnings {
            add_error_message(m, text);
        }
    };
    'exit: {
        if nsgroups == 1 {
            let u = &orig.polymer.as_ref().expect("polymer").units[0];
            if u.type_ == POLYMER_STY_COP {
                treat_err(
                    &mut err,
                    9001,
                    p_str_err,
                    "Copolymer must contain more than one unit",
                );
                break 'exit;
            }
            if matches!(
                u.subtype,
                POLYMER_SST_RAN | POLYMER_SST_ALT | POLYMER_SST_BLK
            ) {
                treat_err(
                    &mut err,
                    9002,
                    p_str_err,
                    "Single polymer unit may not be RAN/ALT/BLO",
                );
                break 'exit;
            }
        }
        // for each CRU
        for i in 0..nsgroups {
            let (nb, na) = {
                let u = &orig.polymer.as_ref().expect("polymer").units[i];
                (u.nb, u.na)
            };
            if nb != 0 && nb != 2 {
                treat_err(
                    &mut err,
                    9003,
                    p_str_err,
                    "Number of crossing bonds in polymer unit is not 0 or 2",
                );
                break 'exit;
            }
            if na < 1 {
                treat_err(&mut err, 9004, p_str_err, "Empty polymer unit");
                break 'exit;
            }
            if na > nat {
                treat_err(&mut err, 9005, p_str_err, "Too large polymer unit");
                break 'exit;
            }
            let bad_atom = orig.polymer.as_ref().expect("polymer").units[i]
                .alist
                .iter()
                .take(na as usize)
                .any(|&a| a < 1 || a > nat);
            if bad_atom {
                treat_err(
                    &mut err,
                    9006,
                    p_str_err,
                    "Invalid atom number in polymer unit",
                );
                break 'exit;
            }
            let mut unit = std::mem::take(&mut orig.polymer.as_mut().expect("polymer").units[i]);
            let e = oad_polymer_unit_set_ends_and_caps(&mut unit, orig, p_str_err);
            if e != 0 {
                err = e;
            }
            if err == 0 {
                unit.nbkbonds = 0;
                unit.cyclizable = CLOSING_SRU_NOT_APPLICABLE;
                unit.cyclized = 0;
            }
            orig.polymer.as_mut().expect("polymer").units[i] = unit;
            if err != 0 {
                break 'exit;
            }
        }

        err = oad_validate_and_sort_out_pseudo_element_atoms(
            orig,
            treat_polymers,
            b_np_zz,
            p_str_err,
        );
        if err != 0 {
            break 'exit;
        }
        if orig.n_zy > 0 && b_np_zz == 0 {
            treat_err(
                &mut err,
                9,
                p_str_err,
                "Non-polymer-related Zz/star atoms are not allowed",
            );
            break 'exit;
        }
        if orig.polymer.is_none() || orig.valid_polymer == 0 {
            break 'exit;
        }

        let n_pzz = orig.polymer.as_ref().expect("polymer").n_pzz;
        if n_pzz > 0 {
            if orig.polymer.as_ref().expect("polymer").treat == super::ichiparm::POLYMERS_NO {
                treat_err(
                    &mut err,
                    9,
                    p_str_err,
                    "Pseudoelement endgroups are not allowed",
                );
                break 'exit;
            }
            let pzz: Vec<i32> = (0..orig.num_inp_atoms)
                .filter(|&k| elname_is(orig, k, "Zz"))
                .map(|k| k as i32 + 1)
                .collect();
            orig.polymer.as_mut().expect("polymer").pzz = Some(pzz);
        }

        // copolymers must include more than one SRU
        let n = orig.polymer.as_ref().expect("polymer").n.max(0) as usize;
        for i in 0..n {
            let pd = orig.polymer.as_ref().expect("polymer");
            let (utype, unb) = (pd.units[i].type_, pd.units[i].nb);
            if utype != POLYMER_STY_COP && utype != POLYMER_STY_SRU {
                continue;
            }
            if unb > 0 {
                if utype == POLYMER_STY_COP {
                    treat_err(
                        &mut err,
                        9026,
                        p_str_err,
                        "Polymer COP unit contains bracket-crossing bonds, not supported",
                    );
                    break 'exit;
                }
                continue;
            }
            let mut in_units = 0;
            for j in 0..n {
                if pd.units[j].type_ == POLYMER_STY_COP {
                    continue;
                }
                let aj = &pd.units[j].alist[..pd.units[j].na as usize];
                let ai = &pd.units[i].alist[..pd.units[i].na as usize];
                if super::util::is_ilist_inside(aj, ai) {
                    in_units += 1;
                    if in_units == 2 {
                        break;
                    }
                }
            }
            if in_units > 1 {
                if utype != POLYMER_STY_COP {
                    orig.polymer.as_mut().expect("polymer").units[i].type_ = POLYMER_STY_COP;
                    warn(p_str_err, "Convert multiple-subunits unit to copolymer");
                }
            } else if utype == POLYMER_STY_COP {
                treat_err(
                    &mut err,
                    9027,
                    p_str_err,
                    "Polymer COP unit contains a single SRU instead of multiple",
                );
                break 'exit;
            }
        }

        let representation = oad_polymer_get_representation(orig.polymer.as_mut());
        let nsg = nsgroups;
        if representation == POLYMER_REPRESENTATION_SOURCE_BASED {
            let pd = orig.polymer.as_mut().expect("polymer");
            for u in pd.units.iter_mut().take(nsg) {
                if u.type_ == POLYMER_STY_SRU {
                    u.type_ = POLYMER_STY_MON;
                    warn(p_str_err, "Converted src-based polymer unit type to MON");
                }
                if u.type_ == POLYMER_STY_COP && u.subtype == POLYMER_SST_NON {
                    u.subtype = POLYMER_SST_RAN;
                    warn(p_str_err, "Set missing copolymer subtype to RAN");
                }
                if u.conn != POLYMER_CONN_NON {
                    u.conn = POLYMER_CONN_NON;
                    warn(
                        p_str_err,
                        "Ignore connection pattern for src-based polymer unit",
                    );
                }
            }
        } else if representation == POLYMER_REPRESENTATION_STRUCTURE_BASED
            || representation == POLYMER_REPRESENTATION_MIXED
        {
            for i in 0..nsg {
                let (utype, unb) = {
                    let u = &orig.polymer.as_ref().expect("polymer").units[i];
                    (u.type_, u.nb)
                };
                if unb == 0 && utype == POLYMER_STY_SRU {
                    orig.polymer.as_mut().expect("polymer").units[i].type_ = POLYMER_STY_COP;
                    warn(p_str_err, "Set copolymer embedding unit mark to COP");
                }
                let utype = orig.polymer.as_ref().expect("polymer").units[i].type_;
                if utype == POLYMER_STY_COP {
                    let u = &mut orig.polymer.as_mut().expect("polymer").units[i];
                    u.cyclizable = CLOSING_SRU_NOT_APPLICABLE;
                    if u.subtype == POLYMER_SST_NON {
                        u.subtype = POLYMER_SST_RAN;
                        warn(p_str_err, "Set missing copolymer subtype to RAN");
                    }
                    continue;
                }
                if utype == POLYMER_STY_MON {
                    continue;
                }
                // SRU with end groups or stars
                let pd = orig.polymer.as_ref().expect("polymer");
                let u = &pd.units[i];
                let pzz: &[i32] = pd.pzz.as_deref().unwrap_or(&[]);
                for k in 0..u.nb as usize {
                    let (a1, a2) = (u.blist[2 * k], u.blist[2 * k + 1]);
                    let is_h = |a: i32| {
                        elname_is(orig, a as usize - 1, "H")
                            || elname_is(orig, a as usize - 1, "D")
                            || elname_is(orig, a as usize - 1, "T")
                    };
                    if is_h(a1) {
                        treat_err(
                            &mut err,
                            9030,
                            p_str_err,
                            "H as polymer end group is not supported",
                        );
                        break 'exit;
                    }
                    if is_h(a2) {
                        treat_err(
                            &mut err,
                            9031,
                            p_str_err,
                            "H as polymer end group is not supported",
                        );
                        break 'exit;
                    }
                    let na = u.na as usize;
                    let out1 = !u.alist[..na].contains(&a1) || pzz.contains(&a1);
                    let out2 = !u.alist[..na].contains(&a2) || pzz.contains(&a2);
                    if out1 && out2 {
                        treat_err(
                            &mut err,
                            9032,
                            p_str_err,
                            "Caps of polymer unit lie inside it",
                        );
                        break 'exit;
                    }
                }
                let (utype, cap1, cap2, na, conn) = (u.type_, u.cap1, u.cap2, u.na, u.conn);
                if matches!(
                    utype,
                    POLYMER_STY_SRU | POLYMER_STY_MOD | POLYMER_STY_CRO | POLYMER_STY_MER
                ) {
                    let mut conn = conn;
                    let mut cyc = orig.polymer.as_ref().expect("polymer").units[i].cyclizable;
                    if conn == POLYMER_CONN_NON {
                        warn(p_str_err, "Set missing copolymer unit connection to EU");
                        conn = POLYMER_CONN_EU;
                    }
                    if cap1 != 0 && cap2 != 0 {
                        cyc = match na {
                            1 => CLOSING_SRU_DIRADICAL,
                            2 => CLOSING_SRU_HIGHER_ORDER_BOND,
                            _ => CLOSING_SRU_RING,
                        };
                    }
                    if conn != POLYMER_CONN_HT {
                        cyc = CLOSING_SRU_NOT_APPLICABLE;
                    }
                    let nb_total = orig.num_inp_bonds;
                    let u = &mut orig.polymer.as_mut().expect("polymer").units[i];
                    u.conn = conn;
                    u.cyclizable = cyc;
                    if cyc != CLOSING_SRU_NOT_APPLICABLE {
                        u.maxbkbonds = nb_total + 2;
                        u.bkbonds = vec![vec![0; 2]; u.maxbkbonds as usize];
                    }
                }
            }
        } else {
            treat_err(
                &mut err,
                9035,
                p_str_err,
                "Invalid kind of polymer representation",
            );
            break 'exit;
        }
        orig.valid_polymer = 1;
    }
    if err != 0 {
        orig.valid_polymer = 0;
    }
    err
}
