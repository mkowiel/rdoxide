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

// ---------------------------------------------------------------------------
// Polymer helpers: bond edits of the original input, ring systems, atom
// seniority and the working set of the polymer layer

const MAXVAL: usize = super::inpdef::MAXVAL;
const INCHI_BOND_TYPE_DOUBLE: i32 = 2;
const INCHI_BOND_TYPE_TRIPLE: i32 = 3;

/// `OrigAtData_RemoveHalfBond` (0-based atoms)
fn orig_at_data_remove_half_bond(
    at: &mut [super::inpdef::InpAtom],
    this_atom: i32,
    other_atom: i32,
    bond_type: &mut i32,
    bond_stereo: &mut i32,
) -> bool {
    if this_atom < 0 || other_atom < 0 {
        return false;
    }
    let a = &mut at[this_atom as usize];
    for k in 0..a.valence as usize {
        if i32::from(a.neighbor[k]) != other_atom {
            continue;
        }
        *bond_type = i32::from(a.bond_type[k]);
        *bond_stereo = i32::from(a.bond_stereo[k]);
        for kk in k + 1..a.valence as usize {
            a.neighbor[kk - 1] = a.neighbor[kk];
            a.bond_type[kk - 1] = a.bond_type[kk];
            a.bond_stereo[kk - 1] = a.bond_stereo[kk];
        }
        for kk in (a.valence as usize - 1)..MAXVAL {
            a.neighbor[kk] = 0;
            a.bond_type[kk] = 0;
            a.bond_stereo[kk] = 0;
        }
        return true;
    }
    false
}

/// `OrigAtData_RemoveBond`
pub fn orig_at_data_remove_bond(
    at: &mut [super::inpdef::InpAtom],
    this_atom: i32,
    other_atom: i32,
    bond_type: &mut i32,
    bond_stereo: &mut i32,
    num_inp_bonds: &mut i32,
) -> bool {
    if this_atom < 0 || other_atom < 0 {
        return false;
    }
    let mut del = i32::from(orig_at_data_remove_half_bond(
        at,
        this_atom,
        other_atom,
        bond_type,
        bond_stereo,
    ));
    del += i32::from(orig_at_data_remove_half_bond(
        at,
        other_atom,
        this_atom,
        bond_type,
        bond_stereo,
    ));
    if del == 2 {
        *num_inp_bonds -= 1;
        at[this_atom as usize].valence -= 1;
        at[this_atom as usize].chem_bonds_valence -= *bond_type as i8;
        at[other_atom as usize].valence -= 1;
        at[other_atom as usize].chem_bonds_valence -= *bond_type as i8;
        return true;
    }
    false
}

/// `OrigAtData_AddBond`
pub fn orig_at_data_add_bond(
    at: &mut [super::inpdef::InpAtom],
    this_atom: i32,
    other_atom: i32,
    mut bond_type: i32,
    bond_stereo: i32,
    num_bonds: &mut i32,
) -> bool {
    let (ta, oa) = (this_atom as usize, other_atom as usize);
    if at[ta].valence as usize >= MAXVAL || at[oa].valence as usize >= MAXVAL {
        return false;
    }
    if bond_type != INCHI_BOND_TYPE_DOUBLE && bond_type != INCHI_BOND_TYPE_TRIPLE {
        bond_type = 1;
    }
    let a = &mut at[ta];
    let k = a.valence as usize;
    if !a.neighbor[..k].iter().any(|&n| i32::from(n) == other_atom) {
        a.neighbor[k] = other_atom as u16;
        a.bond_type[k] = bond_type as u8;
        a.bond_stereo[k] = bond_stereo as i8;
        a.chem_bonds_valence += bond_type as i8;
        a.valence += 1;
    }
    let a = &mut at[oa];
    let k = a.valence as usize;
    if !a.neighbor[..k].iter().any(|&n| i32::from(n) == this_atom) && k < MAXVAL {
        a.neighbor[k] = this_atom as u16;
        a.bond_type[k] = bond_type as u8;
        a.bond_stereo[k] = bond_stereo as i8;
        a.chem_bonds_valence += bond_type as i8;
        a.valence += 1;
    }
    *num_bonds += 1;
    true
}

/// `OrigAtData_AddSingleStereolessBond`
pub fn orig_at_data_add_single_stereoless_bond(
    at: &mut [super::inpdef::InpAtom],
    this_atom: i32,
    other_atom: i32,
    num_bonds: &mut i32,
) -> bool {
    orig_at_data_add_bond(at, this_atom, other_atom, 1, 0, num_bonds)
}

/// `OAD_Polymer_FindRingSystems`: `num_ring_sys[orig_at_number]` (-1 for none)
/// and `size_ring_sys[orig_at_number]` for the atoms in rings systems of more
/// than two atoms; returns the number of such atoms.
pub fn oad_polymer_find_ring_systems(
    pd: &OadPolymer,
    at: &mut [super::inpdef::InpAtom],
    nat: usize,
    num_inp_bonds: &mut i32,
    num_ring_sys: &mut [i32],
    size_ring_sys: &mut [i32],
    start: usize,
) -> i32 {
    let (mut bt, mut bs) = (0, 0);
    let n = pd.n.max(0) as usize;
    for u in pd.units.iter().take(n) {
        if u.cyclized != 0 {
            orig_at_data_remove_bond(
                at,
                u.end_atom1 - 1,
                u.end_atom2 - 1,
                &mut bt,
                &mut bs,
                num_inp_bonds,
            );
        }
    }
    super::ichinorm::mark_ring_systems_inp(&mut at[..nat], start);
    for v in num_ring_sys.iter_mut().take(nat + 1) {
        *v = -1;
    }
    for a in at.iter().take(nat) {
        if a.n_num_at_in_ring_system > 2 {
            let atnum = a.orig_at_number as usize;
            num_ring_sys[atnum] = i32::from(a.n_ring_system);
            size_ring_sys[atnum] = i32::from(a.n_num_at_in_ring_system);
        }
    }
    super::ichinorm::unmark_ring_systems_inp(&mut at[..nat]);
    let nrings = num_ring_sys.iter().take(nat).filter(|&&v| v > -1).count() as i32;
    for u in pd.units.iter().take(n) {
        if u.cyclized != 0 {
            orig_at_data_add_single_stereoless_bond(
                at,
                u.end_atom1 - 1,
                u.end_atom2 - 1,
                num_inp_bonds,
            );
        }
    }
    nrings
}

/// `OAD_AtProps`
#[derive(Debug, Clone, Copy, Default)]
pub struct OadAtProps {
    pub erank: i32,
    pub ring_erank: i32,
    pub ring_size: i32,
    pub ring_num: i32,
}

const ERANK_RULE2: &[i32] = &[
    0, 1, 198, 197, 196, 202, 2, 216, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182, 181,
    180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164, 163, 162,
    161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146, 145, 144, 143,
    142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128, 127, 126, 125, 124,
    123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 105,
    104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83,
    82, 81,
];
const ERANK_RULE4: &[i32] = &[
    0, 1, 198, 197, 196, 202, 2, 211, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182, 181,
    180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164, 163, 162,
    161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146, 145, 144, 143,
    142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128, 127, 126, 125, 124,
    123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 105,
    104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83,
    82, 81,
];

/// `OAD_Polymer_SetAtProps`: `aprops` has `nat + 2` elements; `cano_nums` maps
/// the original atom numbers to canonical numbers (0-based) when given.
pub fn oad_polymer_set_at_props(
    pd: &OadPolymer,
    at: &mut [super::inpdef::InpAtom],
    nat: usize,
    num_inp_bonds: &mut i32,
    aprops: &mut [OadAtProps],
    cano_nums: Option<&[i32]>,
) {
    for k in 0..nat {
        let atnum = at[k].orig_at_number as usize;
        let index = cano_nums.map_or(k as i32, |c| c[atnum]);
        let el = at[k].el_number as usize;
        if index >= 0 && el < ERANK_RULE4.len() {
            let p = &mut aprops[index as usize];
            p.erank = ERANK_RULE4[el];
            p.ring_erank = 0;
            p.ring_size = 0;
            p.ring_num = -1;
        }
    }
    let mut num_ring_sys = vec![0i32; nat + 1];
    let mut size_ring_sys = vec![0i32; nat + 1];
    // break temporarily the first 'cyclizing' bond of the SRUs
    let n = pd.n.max(0) as usize;
    let (mut bt, mut bs, mut dummy) = (-1, -1, 0);
    let breakable = |u: &OadPolymerUnit| {
        u.na > 2 && u.nbkbonds > 0 && u.cyclized == 0 && u.cyclizable == CLOSING_SRU_RING
    };
    for u in pd.units.iter().take(n).filter(|u| breakable(u)) {
        orig_at_data_remove_bond(
            at,
            u.bkbonds[0][0] - 1,
            u.bkbonds[0][1] - 1,
            &mut bt,
            &mut bs,
            &mut dummy,
        );
    }
    let nrings = oad_polymer_find_ring_systems(
        pd,
        at,
        nat,
        num_inp_bonds,
        &mut num_ring_sys,
        &mut size_ring_sys,
        0,
    );
    for u in pd.units.iter().take(n).filter(|u| breakable(u)) {
        orig_at_data_add_bond(
            at,
            u.bkbonds[0][0] - 1,
            u.bkbonds[0][1] - 1,
            bt,
            bs,
            &mut dummy,
        );
    }
    if nrings > 0 {
        let idx = |k: usize, at: &[super::inpdef::InpAtom]| -> usize {
            let atnum = at[k].orig_at_number as usize;
            cano_nums.map_or(k as i32, |c| c[atnum] + 1) as usize
        };
        let mut max_ring_num = 0;
        for k in 0..nat {
            let atnum = at[k].orig_at_number as usize;
            let index = idx(k, at);
            if num_ring_sys[atnum] >= 0 {
                aprops[index].ring_num = num_ring_sys[atnum];
                max_ring_num = max_ring_num.max(aprops[index].ring_num);
                aprops[index].ring_size = size_ring_sys[atnum];
            }
        }
        for i in 0..=max_ring_num {
            let mut max_erank = 0;
            for k in 0..nat {
                if aprops[idx(k, at)].ring_num == i {
                    max_erank = max_erank.max(
                        ERANK_RULE2
                            .get(at[k].el_number as usize)
                            .copied()
                            .unwrap_or(0),
                    );
                }
            }
            for k in 0..nat {
                let index = idx(k, at);
                if aprops[index].ring_num == i && aprops[index].ring_size > 2 {
                    aprops[index].ring_erank = max_erank;
                }
            }
        }
    }
}

/// `OAD_Polymer_CompareRanksOfTwoAtoms`
pub fn oad_polymer_compare_ranks_of_two_atoms(
    atom1: i32,
    atom2: i32,
    aprops: &[OadAtProps],
) -> i32 {
    const HETEROCYC: i32 = 3;
    const HETEROAT: i32 = 2;
    const CARBOCYC: i32 = 1;
    const CARBOAT: i32 = 0;
    let (a1, a2) = (atom1 - 1, atom2 - 1);
    if a1 < 0 || a2 < 0 {
        return 0;
    }
    let (p1, p2) = (&aprops[a1 as usize], &aprops[a2 as usize]);
    let kind = |p: &OadAtProps| {
        if p.ring_size > 2 {
            if p.ring_erank <= 2 {
                CARBOCYC
            } else {
                HETEROCYC
            }
        } else if p.erank == 2 {
            CARBOAT
        } else {
            HETEROAT
        }
    };
    let (t1, t2) = (kind(p1), kind(p2));
    let lower_first = |x: i32, y: i32| -> i32 {
        // smaller value is more senior (returns 1)
        if x < y {
            1
        } else if x > y {
            -1
        } else {
            0
        }
    };
    if t1 == HETEROCYC && t2 == HETEROCYC {
        let r = lower_first(p1.ring_erank, p2.ring_erank);
        if r != 0 {
            return r;
        }
        return lower_first(p1.ring_size, p2.ring_size);
    } else if t1 == HETEROCYC {
        return -1;
    } else if t2 == HETEROCYC {
        return 1;
    }
    if t1 == HETEROAT && t2 == HETEROAT {
        return lower_first(p1.erank, p2.erank);
    } else if t1 == HETEROAT {
        return -1;
    } else if t2 == HETEROAT {
        return 1;
    }
    if t1 == CARBOCYC && t2 == CARBOCYC {
        return lower_first(p1.ring_size, p2.ring_size);
    } else if t1 == CARBOCYC {
        return -1;
    } else if t2 == CARBOCYC {
        return 1;
    }
    0
}

/// `OAD_Polymer_IsFirstAtomRankLower`
pub fn oad_polymer_is_first_atom_rank_lower(atom1: i32, atom2: i32, aprops: &[OadAtProps]) -> i32 {
    let r = oad_polymer_compare_ranks_of_two_atoms(atom1, atom2, aprops);
    if r != 0 {
        return r;
    }
    if atom1 < atom2 {
        1
    } else if atom1 > atom2 {
        -1
    } else {
        0
    }
}

/// `OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves`: returns 1 when a
/// crossing bond is invalid.
fn oad_polymer_unit_order_bond_atoms_and_bonds_themselves(
    u: &mut OadPolymerUnit,
    star_atoms: &[i32],
) -> i32 {
    for k in 0..u.nb as usize {
        let (a1, a2) = (u.blist[2 * k], u.blist[2 * k + 1]);
        let na = u.na as usize;
        let a1_out = !u.alist[..na].contains(&a1) || star_atoms.contains(&a1);
        let a2_out = !u.alist[..na].contains(&a2) || star_atoms.contains(&a2);
        if a1_out && a2_out {
            return 1;
        }
        if a2_out {
            u.blist[2 * k] = a2;
            u.blist[2 * k + 1] = a1;
        }
    }
    if u.nb == 2 {
        let (b1a1, b1a2, b2a1, b2a2) = (u.blist[0], u.blist[1], u.blist[2], u.blist[3]);
        if b1a1 > b2a1 {
            u.blist[0] = b2a1;
            u.blist[1] = b2a2;
            u.blist[2] = b1a1;
            u.blist[3] = b1a2;
        }
    }
    0
}

/// `OAD_PolymerUnit_CompareAtomListsMod`
fn oad_polymer_unit_compare_atom_lists_mod(u1: &OadPolymerUnit, u2: &OadPolymerUnit) -> i32 {
    if u1.na != u2.na {
        return if u1.na < u2.na { -1 } else { 1 };
    }
    for i in 0..u1.na as usize {
        if u1.alist[i] != u2.alist[i] {
            return if u1.alist[i] < u2.alist[i] { -1 } else { 1 };
        }
    }
    0
}

/// `OAD_Polymer_PrepareWorkingSet`: `p.pzz` and the lists of `units2` (copies
/// of the units of `p`) get canonical numbers + 1; fills `unum` with the unit
/// numbers in the order of their atom lists. `cano_nums` is indexed by the
/// original atom number (and holds -1 for removed atoms). Returns the error.
pub fn oad_polymer_prepare_working_set(
    p: &mut OadPolymer,
    cano_nums: &[i32],
    units2: &mut [OadPolymerUnit],
    unum: &mut [usize],
) -> i32 {
    let n_pzz = p.n_pzz.max(0) as usize;
    if let Some(pzz) = p.pzz.as_mut() {
        for k in 0..n_pzz {
            let c = cano_nums[pzz[k] as usize] + 1;
            if c == 0 {
                return 10;
            }
            pzz[k] = c;
        }
    }
    let n = p.n.max(0) as usize;
    for u in units2.iter_mut().take(n) {
        let mut na_new = 0usize;
        for k in 0..u.na as usize {
            let c = cano_nums[u.alist[k] as usize] + 1;
            if c == 0 {
                continue;
            }
            u.alist[na_new] = c;
            na_new += 1;
        }
        u.na = na_new as i32;
        for k in 0..2 * u.nb as usize {
            let c = cano_nums[u.blist[k] as usize] + 1;
            if c == 0 {
                return 11;
            }
            u.blist[k] = c;
        }
        for field in 0..4 {
            let v = match field {
                0 => u.cap1,
                1 => u.cap2,
                2 => u.end_atom1,
                _ => u.end_atom2,
            };
            let c = if v >= 0 && (v as usize) < cano_nums.len() {
                cano_nums[v as usize] + 1
            } else {
                0
            };
            if c == 0 {
                return 11;
            }
            match field {
                0 => u.cap1 = c,
                1 => u.cap2 = c,
                2 => u.end_atom1 = c,
                _ => u.end_atom2 = c,
            }
        }
        for k in 0..u.nbkbonds as usize {
            let c1 = cano_nums[u.bkbonds[k][0] as usize] + 1;
            if c1 == 0 {
                continue;
            }
            let c2 = cano_nums[u.bkbonds[k][1] as usize] + 1;
            if c2 == 0 {
                continue;
            }
            u.bkbonds[k][0] = c1.min(c2);
            u.bkbonds[k][1] = c1.max(c2);
        }
    }
    let pzz: Vec<i32> = p.pzz.clone().unwrap_or_default();
    for u in units2.iter_mut().take(n) {
        let na = u.na as usize;
        super::ichisort::iisort(&mut u.alist[..na]);
        if oad_polymer_unit_order_bond_atoms_and_bonds_themselves(u, &pzz[..n_pzz.min(pzz.len())])
            != 0
        {
            return 12;
        }
    }
    for (i, v) in unum.iter_mut().enumerate().take(n) {
        *v = i;
    }
    for i in 1..n {
        let tmp = unum[i];
        let mut j = i as isize - 1;
        while j >= 0
            && oad_polymer_unit_compare_atom_lists_mod(&units2[unum[j as usize]], &units2[tmp]) > 0
        {
            unum[j as usize + 1] = unum[j as usize];
            j -= 1;
        }
        unum[(j + 1) as usize] = tmp;
    }
    0
}

/// `OrigAtData_IncreaseBondOrder`: number of the bond ends raised
pub fn orig_at_data_increase_bond_order(
    at: &mut [super::inpdef::InpAtom],
    this_atom: usize,
    other_atom: usize,
) -> i32 {
    use super::inpdef::MAXVAL;
    if at[this_atom].valence as usize >= MAXVAL || at[other_atom].valence as usize >= MAXVAL {
        return 0;
    }
    let mut n_up = 0;
    for (from, to) in [(this_atom, other_atom), (other_atom, this_atom)] {
        let a = &mut at[from];
        if a.chem_bonds_valence as usize > MAXVAL - 1 {
            return 0;
        }
        for i in 0..a.valence as usize {
            if a.neighbor[i] as usize != to {
                continue;
            }
            if a.bond_type[i] > 3 {
                return 0;
            }
            a.bond_type[i] += 1;
            a.chem_bonds_valence += 1;
            n_up += 1;
            break;
        }
    }
    n_up
}

/// `OAD_PolymerUnit_HasMetal`
fn oad_polymer_unit_has_metal(u: &OadPolymerUnit, at: &[super::inpdef::InpAtom]) -> bool {
    u.alist
        .iter()
        .take(u.na as usize)
        .any(|&a| super::util::is_el_a_metal(i32::from(at[a as usize - 1].el_number)))
}

/// `OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms`
fn oad_polymer_unit_unlink_caps_and_connect_end_atoms(
    unit: &mut OadPolymerUnit,
    orig: &mut OrigAtomData,
) {
    if unit.cyclizable == CLOSING_SRU_NOT_APPLICABLE {
        return;
    }
    let (mut bt, mut bs) = (0, 0);
    let (c1, c2, e1, e2) = (
        unit.cap1 - 1,
        unit.cap2 - 1,
        unit.end_atom1 - 1,
        unit.end_atom2 - 1,
    );
    let mut nb = orig.num_inp_bonds;
    if unit.cyclizable == CLOSING_SRU_RING {
        orig_at_data_remove_bond(&mut orig.at, c1, e1, &mut bt, &mut bs, &mut nb);
        orig_at_data_remove_bond(&mut orig.at, c2, e2, &mut bt, &mut bs, &mut nb);
        orig_at_data_add_single_stereoless_bond(&mut orig.at, e1, e2, &mut nb);
    } else if unit.cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND {
        orig_at_data_increase_bond_order(&mut orig.at, e1 as usize, e2 as usize);
        orig_at_data_remove_bond(&mut orig.at, c1, e1, &mut bt, &mut bs, &mut nb);
        orig_at_data_remove_bond(&mut orig.at, c2, e2, &mut bt, &mut bs, &mut nb);
    } else if unit.cyclizable == CLOSING_SRU_DIRADICAL {
        orig.at[e1 as usize].radical = RADICAL_TRIPLET as i8;
        orig_at_data_remove_bond(&mut orig.at, c1, e1, &mut bt, &mut bs, &mut nb);
        orig_at_data_remove_bond(&mut orig.at, c2, e2, &mut bt, &mut bs, &mut nb);
    }
    orig.num_inp_bonds = nb;
    unit.cyclized = 1;
}

/// `OAD_Polymer_CyclizeCloseableUnits`: returns the error code
pub fn oad_polymer_cyclize_closeable_units(
    orig: &mut OrigAtomData,
    p_str_err: &mut String,
    b_no_warnings: bool,
) -> i32 {
    let n = orig.polymer.as_ref().map_or(0, |p| p.units.len());
    for i in 0..n {
        let mut unit = std::mem::take(&mut orig.polymer.as_mut().expect("polymer").units[i]);
        let mut err = 0;
        'unit: {
            if unit.cyclizable == CLOSING_SRU_NOT_APPLICABLE {
                break 'unit;
            }
            err = oad_polymer_unit_set_ends_and_caps(&mut unit, orig, p_str_err);
            if err != 0 || unit.cyclizable == 0 {
                break 'unit;
            }
            if oad_polymer_unit_has_metal(&unit, &orig.at)
                && unit.cyclizable == CLOSING_SRU_RING
                && !b_no_warnings
            {
                add_error_message(
                    p_str_err,
                    "Frame shift in metallated polymer unit may be missed",
                );
            }
            oad_polymer_unit_unlink_caps_and_connect_end_atoms(&mut unit, orig);
        }
        orig.polymer.as_mut().expect("polymer").units[i] = unit;
        if err != 0 {
            return err;
        }
    }
    0
}

/// `subgraf`: a light-weight subgraph of the original atoms (`subgraf_new`)
struct Subgraf {
    nodes: Vec<i32>,
    /// original atom number -> node number (-1: not in the subgraph)
    orig2node: Vec<i32>,
    adj: Vec<Vec<usize>>,
}

fn subgraf_new(at: &[super::inpdef::InpAtom], nat: usize, nodes: &[i32]) -> Subgraf {
    let mut orig2node = vec![-1; nat + 1];
    for (i, &n) in nodes.iter().enumerate() {
        orig2node[n as usize] = i as i32;
    }
    let adj = nodes
        .iter()
        .map(|&n| {
            let a = &at[n as usize - 1];
            (0..a.valence as usize)
                .filter_map(|j| {
                    let nbr = orig2node[a.neighbor[j] as usize + 1];
                    (nbr >= 0).then_some(nbr as usize)
                })
                .collect()
        })
        .collect();
    Subgraf {
        nodes: nodes.to_vec(),
        orig2node,
        adj,
    }
}

/// `subgraf_pathfinder_run` (no forbidden edges, no atom list): all simple
/// paths from the last node of `seen` to `end`; collects their bonds as pairs of
/// original atom numbers
fn subgraf_pathfinder_run(
    sg: &Subgraf,
    seen: &mut Vec<usize>,
    end: usize,
    bonds: &mut Vec<[i32; 2]>,
) {
    let Some(&node0) = seen.last() else {
        return;
    };
    for &node in &sg.adj[node0] {
        if seen.contains(&node) {
            continue;
        }
        if node == end {
            seen.push(node);
            for k in 1..seen.len() {
                let (a1, a2) = (sg.nodes[seen[k - 1]], sg.nodes[seen[k]]);
                if !bonds
                    .iter()
                    .any(|b| (b[0] == a1 && b[1] == a2) || (b[0] == a2 && b[1] == a1))
                {
                    bonds.push([a1, a2]);
                }
            }
            seen.pop();
            break;
        }
    }
    for &node in &sg.adj[node0] {
        if node == end || seen.contains(&node) {
            continue;
        }
        seen.push(node);
        subgraf_pathfinder_run(sg, seen, end, bonds);
        seen.pop();
    }
}

/// `OAD_PolymerUnit_RemoveLinkFromCRUChain`
fn oad_polymer_unit_remove_link_from_cru_chain(
    at1: i32,
    at2: i32,
    nbonds: &mut i32,
    bonds: &mut [Vec<i32>],
) {
    for p in 0..*nbonds as usize {
        if bonds[p][0] == at1 && bonds[p][1] == at2 {
            for q in p + 1..*nbonds as usize {
                bonds[q - 1][0] = bonds[q][0];
                bonds[q - 1][1] = bonds[q][1];
            }
            *nbonds -= 1;
            break;
        }
    }
}

/// `OAD_Polymer_FindBackbones`: for all frame-shiftable units, the bonds
/// potentially involved in the frame shift. `at_data` is the preprocessed
/// structure; `taut_bonds` holds the tautomeric bonds of the composite
/// normalized structure as pairs of original atom numbers.
pub fn oad_polymer_find_backbones(
    at_data: &mut OrigAtomData,
    pd: &mut OadPolymer,
    taut_bonds: &std::collections::HashSet<(i32, i32)>,
) {
    let nat = at_data.num_inp_atoms;
    for i in 0..pd.units.len().min(pd.n.max(0) as usize) {
        if pd.units[i].cyclizable == CLOSING_SRU_NOT_APPLICABLE {
            continue;
        }
        {
            let u = &mut pd.units[i];
            let sg = subgraf_new(&at_data.at, nat, &u.alist[..u.na as usize]);
            let (start, end) = (
                sg.orig2node[u.end_atom1 as usize],
                sg.orig2node[u.end_atom2 as usize],
            );
            u.nbkbonds = 0;
            if start < 0 || end < 0 {
                continue;
            }
            let mut bonds = Vec::new();
            subgraf_pathfinder_run(&sg, &mut vec![start as usize], end as usize, &mut bonds);
            u.nbkbonds = bonds.len() as i32;
            for (k, b) in bonds.iter().enumerate() {
                if k < u.bkbonds.len() {
                    u.bkbonds[k] = b.to_vec();
                } else {
                    u.bkbonds.push(b.to_vec());
                }
            }
        }
        if pd.units[i].nbkbonds < 2 {
            continue;
        }
        // OAD_PolymerUnit_DelistIntraRingBackboneBonds
        let mut num_ring_sys = vec![-1; nat + 1];
        let mut size_ring_sys = vec![0; nat + 1];
        let start = (pd.units[i].end_atom1 - 1) as usize;
        let nrings = {
            let snapshot = pd.clone();
            let mut nb = at_data.num_inp_bonds;
            let n = oad_polymer_find_ring_systems(
                &snapshot,
                &mut at_data.at,
                nat,
                &mut nb,
                &mut num_ring_sys,
                &mut size_ring_sys,
                start,
            );
            at_data.num_inp_bonds = nb;
            n
        };
        let u = &mut pd.units[i];
        if nrings != 0 {
            let mut j = 0;
            loop {
                let (a1, a2) = (u.bkbonds[j][0], u.bkbonds[j][1]);
                let (r1, r2) = (num_ring_sys[a1 as usize], num_ring_sys[a2 as usize]);
                if r1 == r2 && r1 != -1 {
                    oad_polymer_unit_remove_link_from_cru_chain(
                        a1,
                        a2,
                        &mut u.nbkbonds,
                        &mut u.bkbonds,
                    );
                } else {
                    j += 1;
                }
                if j as i32 >= u.nbkbonds {
                    break;
                }
            }
        }
        // OAD_PolymerUnit_DelistHighOrderBackboneBonds (tautomeric bonds only)
        if u.na >= 2 && u.nb >= 2 && u.nbkbonds >= 1 {
            let mut j = 0;
            loop {
                let (a1, a2) = (u.bkbonds[j][0], u.bkbonds[j][1]);
                if taut_bonds.contains(&(a1.min(a2), a1.max(a2))) {
                    oad_polymer_unit_remove_link_from_cru_chain(
                        a1,
                        a2,
                        &mut u.nbkbonds,
                        &mut u.bkbonds,
                    );
                } else {
                    j += 1;
                }
                if j as i32 >= u.nbkbonds {
                    break;
                }
            }
        }
        if u.nbkbonds == 0 {
            // keep a single frame-shift bond between the original caps of the CRU ends
            u.cyclizable = 1;
            u.nbkbonds = 1;
            u.bkbonds[0][0] = u.end_atom1;
            u.bkbonds[0][1] = u.end_atom2;
        }
    }
}

/// `OAD_PolymerUnit_SetReopeningDetails`: returns the number of backbone bonds
pub fn oad_polymer_unit_set_reopening_details(
    u: &mut OadPolymerUnit,
    at: &[super::inpdef::InpAtom],
) -> i32 {
    if u.nbkbonds == 1 {
        u.end_atom1 = u.bkbonds[0][0];
        u.end_atom2 = u.bkbonds[0][1];
        if u.end_atom1 == u.end_atom2 {
            u.cyclizable = CLOSING_SRU_DIRADICAL;
        } else {
            let a = &at[u.end_atom1 as usize - 1];
            for k in 0..a.valence as usize {
                if i32::from(a.neighbor[k]) == u.end_atom2 - 1 {
                    if a.bond_type[k] > 1 {
                        u.cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
                    }
                    break;
                }
            }
        }
    }
    u.nbkbonds
}

/// `OAD_Polymer_CompareBackboneBondsSeniority`
fn oad_polymer_compare_backbone_bonds_seniority(
    b1: &[i32],
    b2: &[i32],
    aprops: &[OadAtProps],
) -> i32 {
    let (mut b1max, mut b1min) = (b1[0], b1[1]);
    let (mut b2max, mut b2min) = (b2[0], b2[1]);
    if oad_polymer_is_first_atom_rank_lower(b1min, b1max, aprops) == -1 {
        std::mem::swap(&mut b1max, &mut b1min);
    }
    if oad_polymer_is_first_atom_rank_lower(b2min, b2max, aprops) == -1 {
        std::mem::swap(&mut b2max, &mut b2min);
    }
    match oad_polymer_compare_ranks_of_two_atoms(b1max, b2max, aprops) {
        1 => return 1,
        -1 => return -1,
        _ => {}
    }
    match oad_polymer_compare_ranks_of_two_atoms(b1min, b2min, aprops) {
        1 => return -1,
        -1 => return 1,
        _ => {}
    }
    if b1max < b2max {
        return 1;
    }
    if b1max > b2max {
        return -1;
    }
    if b1min < b2min {
        return -1;
    }
    if b1min > b2min {
        return 1;
    }
    0
}

/// `OAD_PolymerUnit_SortBackboneBondsAndSetSeniors`: returns the index of the
/// senior backbone bond
pub fn oad_polymer_unit_sort_backbone_bonds_and_set_seniors(
    u: &mut OadPolymerUnit,
    aprops: &[OadAtProps],
) -> usize {
    let mut senior = 0;
    if u.nbkbonds > 1 {
        let n = u.nbkbonds as usize;
        let mut bnum: Vec<usize> = (0..n).collect();
        for i in 1..n {
            let tmp = bnum[i];
            let mut j = i as isize - 1;
            while j >= 0
                && oad_polymer_compare_backbone_bonds_seniority(
                    &u.bkbonds[bnum[j as usize]],
                    &u.bkbonds[tmp],
                    aprops,
                ) > 0
            {
                bnum[(j + 1) as usize] = bnum[j as usize];
                j -= 1;
            }
            bnum[(j + 1) as usize] = tmp;
        }
        senior = bnum[0];
    }
    if oad_polymer_is_first_atom_rank_lower(u.bkbonds[senior][0], u.bkbonds[senior][1], aprops) == 1
    {
        u.bkbonds[senior].swap(0, 1);
    }
    u.end_atom1 = u.bkbonds[senior][0];
    u.end_atom2 = u.bkbonds[senior][1];
    senior
}

/// `inchi_strtol`-like: the integer at `pos` (after white space and a sign) and the
/// position after it; `(0, pos)` when there are no digits
fn strtol_at(b: &[u8], pos: usize) -> (i32, usize) {
    let mut i = pos;
    while i < b.len() && b[i].is_ascii_whitespace() {
        i += 1;
    }
    let neg = i < b.len() && b[i] == b'-';
    if i < b.len() && (b[i] == b'-' || b[i] == b'+') {
        i += 1;
    }
    let st = i;
    let mut v: i64 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        v = (v * 10 + i64::from(b[i] - b'0')).min(i64::from(i32::MAX));
        i += 1;
    }
    if i == st {
        return (0, pos);
    }
    (if neg { -v as i32 } else { v as i32 }, i)
}

/// `GetFrameShiftInfoFrom105PlusInChI`: `(unit, atom1, atom2)` of the senior
/// backbone bond of each eligible unit in the `/z` layer
fn get_frame_shift_info_from_105plus_inchi(
    sinchi: &str,
    max_crossing: usize,
) -> Vec<(usize, i32, i32)> {
    let b = sinchi.as_bytes();
    let find = |from: usize, c: u8| {
        b.get(from..)
            .and_then(|s| s.iter().position(|&x| x == c))
            .map(|k| k + from)
    };
    let mut out = Vec::new();
    let mut p = sinchi.find("/z");
    let mut iunit = 0;
    while let Some(pp) = p {
        let Some(open) = find(pp + 2, b'(') else {
            break;
        };
        let mut pos = open + 1;
        let mut num = [-1, -1];
        let mut j = 0;
        let mut skip = false;
        loop {
            let (k, q) = strtol_at(b, pos);
            if k == 0 || j >= 2 {
                break;
            }
            num[j] = k;
            j += 1;
            let c = b.get(q).copied().unwrap_or(0);
            if (j == 1 && c == b'-') || c == b')' {
                skip = true;
                break;
            }
            pos = q + 1;
        }
        if !skip && j >= 2 {
            out.push((iunit, num[0], num[1]));
            if out.len() >= max_crossing {
                break;
            }
        }
        p = find(pos, b';');
        iunit += 1;
    }
    out
}

/// `extract_orig_nums_from_auxinfo_string`: `orig[canonical number]`
fn extract_orig_nums_from_auxinfo_string(saux: &str, orig: &mut [i32]) -> bool {
    let b = saux.as_bytes();
    let Some(p) = saux.find("/N:") else {
        return false;
    };
    if !b.get(p + 3).is_some_and(u8::is_ascii_digit) {
        return false;
    }
    let mut pos = p + 3;
    let mut cano = 1;
    loop {
        let (k, q) = strtol_at(b, pos);
        if k == 0 {
            break;
        }
        if cano < orig.len() {
            orig[cano] = k;
        }
        cano += 1;
        match b.get(q).copied() {
            Some(c) if c != 0 && c != b'/' => pos = q + 1,
            _ => break,
        }
    }
    true
}

/// `OAD_Polymer_PrepareFrameShiftEdits`: the bond modifications
/// `(old_a1, old_a2, new_a1, new_a2)`; `None` is an error
pub fn oad_polymer_prepare_frame_shift_edits(
    orig: &OrigAtomData,
    sinchi: &str,
    saux: &str,
) -> Option<Vec<[i32; 4]>> {
    let p = orig.polymer.as_ref()?;
    let nu = p.n.max(0) as usize;
    let mut cano2orig = vec![0; orig.num_inp_atoms + 1];
    if !extract_orig_nums_from_auxinfo_string(saux, &mut cano2orig) {
        return None;
    }
    let mut mods = Vec::new();
    for (_, c1, c2) in get_frame_shift_info_from_105plus_inchi(sinchi, nu) {
        let at = |c: i32| cano2orig.get(c as usize).copied().unwrap_or(0);
        let (end1, end2) = (at(c1), at(c2));
        let mut iu = None;
        for (k, u) in p.units.iter().enumerate().take(nu) {
            if u.blist.is_empty() || u.nb < 2 {
                continue;
            }
            let present = u
                .alist
                .iter()
                .take(u.na as usize)
                .filter(|&&a| a == end1 || a == end2)
                .count();
            if present >= 2 {
                iu = Some(k);
            }
        }
        let u = &p.units[iu?];
        let mut probe = u.clone();
        let mut dummy = String::new();
        let err = oad_polymer_unit_find_ends_and_caps(&mut probe, orig, &mut dummy);
        let (cap1, cap2) = (probe.cap1, probe.cap2);
        if err == 0
            && probe.cap1_is_undef != 0
            && probe.cap2_is_undef != 0
            && end1 != 0
            && end2 != 0
            && cap1 != 0
            && cap2 != 0
        {
            let bl = &u.blist;
            let old_end = |cap: i32| -> Option<i32> {
                if cap == bl[0] {
                    Some(bl[1])
                } else if cap == bl[1] {
                    Some(bl[0])
                } else if cap == bl[2] {
                    Some(bl[3])
                } else if cap == bl[3] {
                    Some(bl[2])
                } else {
                    None
                }
            };
            let (Some(old_end1), Some(old_end2)) = (old_end(cap1), old_end(cap2)) else {
                continue;
            };
            if old_end1 == 0 || old_end2 == 0 || old_end1 == old_end2 {
                continue;
            }
            let same = |a1: i32, a2: i32, b1: i32, b2: i32| {
                (a1 == b1 && a2 == b2) || (a1 == b2 && a2 == b1)
            };
            if same(old_end1, cap1, end1, cap1) && same(old_end2, cap2, end2, cap2) {
                continue;
            }
            if !same(old_end1, cap1, end1, cap1) {
                mods.push([old_end1, cap1, end1, cap1]);
            }
            if !same(old_end2, cap2, end2, cap2) {
                mods.push([old_end2, cap2, end2, cap2]);
            }
            mods.push([end1, end2, old_end1, old_end2]);
        }
    }
    Some(mods)
}

/// `OAD_StructureEdits_Apply` for bond modifications; returns the number of
/// edits done, `None` at error
pub fn oad_structure_edits_apply_mod_bonds(
    orig: &mut OrigAtomData,
    mods: &[[i32; 4]],
) -> Option<i32> {
    let mut n = 0;
    for m in mods {
        let (o1, o2, n1, n2) = (m[0], m[1], m[2], m[3]);
        if (o1 == n1 && o2 == n2) || (o2 == n1 && o1 == n2) {
            continue;
        }
        let (mut bt, mut bs) = (0, 0);
        let mut nb = orig.num_inp_bonds;
        if !orig_at_data_remove_bond(&mut orig.at, o1 - 1, o2 - 1, &mut bt, &mut bs, &mut nb) {
            return None;
        }
        if !orig_at_data_add_bond(&mut orig.at, n1 - 1, n2 - 1, bt, bs, &mut nb) {
            return None;
        }
        orig.num_inp_bonds = nb;
        if let Some(p) = orig.polymer.as_mut() {
            for u in p.units.iter_mut().take(p.n.max(0) as usize) {
                if u.blist.is_empty() {
                    continue;
                }
                let same = |a1: i32, a2: i32| (a1 == o1 && a2 == o2) || (a1 == o2 && a2 == o1);
                if same(u.blist[0], u.blist[1]) {
                    u.blist[0] = n1;
                    u.blist[1] = n2;
                } else if same(u.blist[2], u.blist[3]) {
                    u.blist[2] = n1;
                    u.blist[3] = n2;
                }
            }
        }
        n += 1;
    }
    Some(n)
}
