//! Structure preprocessing of InChI (`INCHI_BASE/src/strutil.c`), the part
//! that corrects drawing peculiarities of `inp_ATOM` arrays:
//! non-uniformly drawn oxoanions and amidinium cations, the doublet-radical
//! neighbor search and the choice of one terminal O among several.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, RADICAL_DOUBLET, RADICAL_SINGLET};
use super::inpdef::InpAtom;

const D_CL: u8 = 17;
const D_BR: u8 = 35;
const D_I: u8 = 53;
const D_AT: u8 = 85;
const D_O: u8 = 8;
const D_S: u8 = 16;
const D_SE: u8 = 34;
const D_TE: u8 = 52;
const D_P: u8 = 15;
const D_C: u8 = 6;
const D_N: u8 = 7;

/// `cmp_iso_atw_diff_component_no`: order by isotopic shift, then by component
pub fn cmp_iso_atw_diff_component_no(a1: &InpAtom, a2: &InpAtom) -> i32 {
    let ret = i32::from(a1.iso_atw_diff) - i32::from(a2.iso_atw_diff);
    if ret == 0 {
        i32::from(a1.component) - i32::from(a2.component)
    } else {
        ret
    }
}

/// `the_only_doublet_neigh`: if atom `i1` is a doublet radical with exactly
/// one doublet neighbor, and that neighbor has exactly one doublet neighbor
/// (`i1`), returns the neighbor and the positions of the two bonds.
pub fn the_only_doublet_neigh(at: &[InpAtom], i1: usize) -> Option<(usize, usize, usize)> {
    let a = &at[i1];
    if i32::from(a.radical) != RADICAL_DOUBLET {
        return None;
    }
    let mut num_rad1 = 0;
    let mut ineigh1 = 0;
    for i in 0..a.valence as usize {
        let b = &at[a.neighbor[i] as usize];
        if i32::from(b.radical) == RADICAL_DOUBLET {
            num_rad1 += 1;
            ineigh1 = i;
        }
    }
    if num_rad1 == 1 {
        let neigh1 = a.neighbor[ineigh1] as usize;
        let a = &at[neigh1];
        let mut num_rad2 = 0;
        let mut ineigh2 = 0;
        for i in 0..a.valence as usize {
            let b = &at[a.neighbor[i] as usize];
            if i32::from(b.radical) == RADICAL_DOUBLET {
                num_rad2 += 1;
                ineigh2 = i;
            }
        }
        if num_rad2 == 1 {
            return Some((neigh1, ineigh1, ineigh2));
        }
    }
    None
}

/// `fix_non_uniform_drawn_oxoanions`: moves the negative charge of a halogen
/// (or highest-oxidation-state chalcogen) to a terminal =O and makes the bond
/// single. `num_changes` is incremented per fix.
pub fn fix_non_uniform_drawn_oxoanions(at: &mut [InpAtom], num_changes: &mut i32) -> i32 {
    const HALOGEN: [u8; 4] = [D_CL, D_BR, D_I, D_AT];
    const HALCOGEN: [u8; 3] = [D_S, D_SE, D_TE];
    for i in 0..at.len() {
        // charged exactly (-1) ...
        if at[i].charge != -1 {
            continue;
        }
        let en_center = at[i].el_number;
        // from the eligible element list ...
        if !HALOGEN.contains(&en_center) {
            if HALCOGEN.contains(&en_center) {
                if at[i].chem_bonds_valence < 7 {
                    continue;
                }
            } else {
                continue;
            }
        }
        // non-radical ...
        if at[i].radical != 0 && i32::from(at[i].radical) != RADICAL_SINGLET {
            continue;
        }
        // the center is found, examine the adjacent terminals
        let mut kk = 0usize;
        let mut jj: i32 = -1;
        let mut min_en = 999i32;
        let mut min_iso = 999i32;
        for k in 0..at[i].valence as usize {
            let j = at[i].neighbor[k] as usize;
            if at[j].valence != 1 {
                continue;
            }
            if i32::from(at[i].bond_type[k]) != BOND_TYPE_DOUBLE {
                continue;
            }
            if at[j].charge != 0 {
                continue;
            }
            if at[j].radical != 0 && i32::from(at[j].radical) != RADICAL_SINGLET {
                continue;
            }
            // of eligible elements ...
            let en_term = at[j].el_number;
            let eligible = match en_term {
                D_O => true,
                D_S => matches!(en_center, D_SE | D_AT | D_TE),
                D_SE => matches!(en_center, D_AT | D_TE),
                D_TE => en_center == D_AT,
                _ => false,
            };
            if !eligible {
                continue;
            }
            // from several candidates, select one with less el. number
            // (the isotopic shift is read from the center, as in the C code)
            let shift = |a: &InpAtom, b: &InpAtom| -> i32 {
                if a.iso_atw_diff > 0 {
                    i32::from(b.iso_atw_diff) - 1
                } else {
                    i32::from(b.iso_atw_diff)
                }
            };
            if i32::from(en_term) < min_en {
                min_en = i32::from(en_term);
                kk = k;
                jj = j as i32;
                min_iso = shift(&at[j], &at[i]);
                continue;
            } else if i32::from(en_term) == min_en {
                let iso = shift(&at[j], &at[i]);
                if iso < min_iso {
                    min_iso = iso;
                    kk = k;
                    jj = j as i32;
                    continue;
                }
            }
        }
        // if OK, apply changes
        if jj >= 0 {
            let jj = jj as usize;
            at[i].charge = 0;
            at[jj].charge = -1;
            at[i].bond_type[kk] = BOND_TYPE_SINGLE as u8;
            at[jj].bond_type[0] = BOND_TYPE_SINGLE as u8;
            at[i].bond_stereo[kk] = 0;
            at[jj].bond_stereo[0] = 0;
            at[i].chem_bonds_valence -= 1;
            at[jj].chem_bonds_valence -= 1;
            *num_changes += 1;
        }
    }
    0
}

/// `fix_non_uniform_drawn_amidiniums`: moves the positive charge of an
/// amidinium cation to a nitrogen and makes the C-N bond double.
pub fn fix_non_uniform_drawn_amidiniums(at: &mut [InpAtom], num_changes: &mut i32) -> i32 {
    const CENTER: [u8; 3] = [D_C, D_S, D_P];
    for i in 0..at.len() {
        if at[i].charge != 1 {
            continue;
        }
        if !CENTER.contains(&at[i].el_number) {
            continue;
        }
        // exactly 3 neighbors connected by single bonds
        if at[i].valence != 3 || at[i].chem_bonds_valence != 3 {
            continue;
        }
        // non-radical
        if at[i].radical != 0 && i32::from(at[i].radical) != RADICAL_SINGLET {
            continue;
        }
        // the neighbors are neutral, two of them are aliphatic N, at least one bears H
        let mut mismatch = false;
        let mut nu_h = 0;
        let mut nu_n = 0;
        let mut kk = 0usize;
        let mut jj: i32 = -1;
        for k in 0..at[i].valence as usize {
            let j = at[i].neighbor[k] as usize;
            if at[j].charge != 0 {
                mismatch = true;
                break;
            }
            if at[j].el_number == D_N {
                if at[j].valence > 3 || at[j].chem_bonds_valence > 3 {
                    mismatch = true;
                    break;
                }
                nu_h += at[j].numh();
                nu_n += 1;
                if jj < 0 {
                    jj = j as i32;
                    kk = k;
                }
            }
        }
        if mismatch || nu_n != 2 || nu_h < 1 {
            continue;
        }
        if jj >= 0 {
            let jj = jj as usize;
            at[i].charge = 0;
            at[jj].charge = 1;
            at[i].bond_type[kk] = BOND_TYPE_DOUBLE as u8;
            let mut k1 = 0;
            while k1 < at[jj].valence as usize && i != at[jj].neighbor[k1] as usize {
                k1 += 1;
            }
            at[jj].bond_type[k1] = BOND_TYPE_DOUBLE as u8;
            at[i].chem_bonds_valence += 1;
            at[jj].chem_bonds_valence += 1;
            // nothing is done with wedge stereo bonds (the wedge is retained)
            *num_changes += 1;
        }
    }
    0
}

/// `nFindOneOM`: selects one of the neighbors `ord_om` (positions in the
/// neighbor list of `at_no`): fewest bonds, then lowest element number, and
/// among terminal atoms the one without or with the lightest isotope. The
/// list is reordered like the C array. Returns -1 if nothing is selected.
pub fn n_find_one_om(at: &[InpAtom], at_no: usize, ord_om: &mut [i32], num_om: usize) -> i32 {
    if num_om == 1 {
        return ord_om[0];
    }
    if num_om < 1 {
        return -1;
    }
    // select neighbors with the min. number of bonds
    let mut num_best = 1;
    let mut n_om = at[at_no].neighbor[ord_om[0] as usize] as usize;
    let mut best_value = i32::from(at[n_om].valence);
    for i in 1..num_om {
        n_om = at[at_no].neighbor[ord_om[i] as usize] as usize;
        let cur_value = i32::from(at[n_om].valence);
        let diff = cur_value - best_value;
        if diff < 0 {
            best_value = cur_value;
            ord_om[0] = ord_om[i];
            num_best = 1;
        } else if diff == 0 {
            ord_om[num_best] = ord_om[i];
            num_best += 1;
        }
    }
    let num_om = num_best;
    if num_om == 1 {
        return ord_om[0];
    }

    // select neighbors with the min. periodic numbers
    num_best = 1;
    n_om = at[at_no].neighbor[ord_om[0] as usize] as usize;
    best_value = i32::from(at[n_om].el_number);
    for i in 1..num_om {
        n_om = at[at_no].neighbor[ord_om[i] as usize] as usize;
        let cur_value = i32::from(at[n_om].el_number);
        let diff = cur_value - best_value;
        if diff < 0 {
            best_value = cur_value;
            ord_om[0] = ord_om[i];
            num_best = 1;
        } else if diff == 0 {
            ord_om[num_best] = ord_om[i];
            num_best += 1;
        }
    }
    let num_om = num_best;
    if num_om == 1 {
        return ord_om[0];
    }

    // if the neighbors are not terminal atoms then reject
    if at[n_om].valence > 1 {
        return -1;
    }

    // terminal atoms: the one without isotope or with the lightest isotope
    num_best = 1;
    n_om = at[at_no].neighbor[ord_om[0] as usize] as usize;
    best_value = i32::from(at[n_om].iso_atw_diff);
    for i in 1..num_om {
        n_om = at[at_no].neighbor[ord_om[i] as usize] as usize;
        // (the C code compares the element number with the isotopic shift)
        let cur_value = i32::from(at[n_om].el_number);
        let diff = cur_value - best_value;
        if (cur_value == 0 && best_value != 0) || diff < 0 {
            best_value = cur_value;
            ord_om[0] = ord_om[i];
            num_best = 1;
        } else if diff == 0 {
            ord_om[num_best] = ord_om[i];
            num_best += 1;
        }
    }
    // return any
    ord_om[0]
}

// ---- remove_ion_pairs and fix_odd_things ----

use super::incomdef::{BOND_TYPE_TRIPLE, RADICAL_SINGLET as RAD_SINGLET};
use super::util::{
    EL_NUMBER_AS, EL_NUMBER_C, EL_NUMBER_H, EL_NUMBER_N, EL_NUMBER_O, EL_NUMBER_P, EL_NUMBER_S,
    EL_NUMBER_SB, EL_NUMBER_SE, EL_NUMBER_TE, get_el_valence, has_other_ion_in_sphere_2,
    has_other_ion_neigh, ion_el_group, n_no_metal_bonds_valence, n_no_metal_neigh_index,
    n_no_metal_num_bonds, n_no_metal_other_neigh_index, n_no_metal_other_neigh_index2, num_of_h,
};

const SINGLE: u8 = BOND_TYPE_SINGLE as u8;
const DOUBLE: u8 = BOND_TYPE_DOUBLE as u8;
const TRIPLE: u8 = BOND_TYPE_TRIPLE as u8;

/// Position of `of` in the neighbor list of `at[n]` (`is_in_the_list(...) -
/// neighbor` in C; the bond is always present in valid input)
fn pos_of(at: &[InpAtom], n: usize, of: usize) -> usize {
    at[n].neighbor[..at[n].valence as usize]
        .iter()
        .position(|&a| a as usize == of)
        .expect("symmetric adjacency lists")
}

fn nb(at: &[InpAtom], i: usize, k: usize) -> usize {
    at[i].neighbor[k] as usize
}

fn nnm_bonds(at: &[InpAtom], i: usize) -> i32 {
    n_no_metal_num_bonds(at, i)
}

fn nnm_val(at: &[InpAtom], i: usize) -> i32 {
    n_no_metal_bonds_valence(at, i)
}

fn grp(at: &[InpAtom], i: usize) -> u8 {
    ion_el_group(at[i].el_number)
}

/// The common action of the two-atom conversions: removes the charges and
/// raises the bond order between `i` (bond position `i1`) and `n` (`i2`).
fn neutralize_pair(at: &mut [InpAtom], i: usize, i1: usize, n: usize, i2: usize) {
    at[i].charge -= 1;
    at[n].charge += 1;
    at[i].chem_bonds_valence += 1;
    at[n].chem_bonds_valence += 1;
    at[i].bond_type[i1] += 1;
    at[n].bond_type[i2] += 1;
}

/// `remove_ion_pairs`: neutralizes ion pairs of the types 1..18 described in
/// the comments of the C code (`O=N(+)-O(-)`, `N(+)-C(-)`, `C(+)-O-N(-)`, ...).
/// Returns the number of changes.
pub fn remove_ion_pairs(at: &mut [InpAtom]) -> i32 {
    let num_atoms = at.len();
    let mut num_changes = 0;
    let (mut num_c_ii, mut num_c_plus, mut num_c_minus) = (0, 0, 0);
    let (mut num_n_plus, mut num_n_minus, mut num_o_plus, mut num_o_minus) = (0, 0, 0, 0);
    // FIX_P_IV_Plus_O_Minus: P(IV)(+) atoms, counted among the N group
    let mut num_p_iv_plus = 0;

    // count candidates
    for a in at.iter() {
        let chrg = a.charge;
        if chrg == 1 || chrg == -1 {
            match ion_el_group(a.el_number) {
                EL_NUMBER_C => {
                    if chrg > 0 {
                        num_c_plus += 1;
                    } else {
                        num_c_minus += 1;
                    }
                }
                EL_NUMBER_O => {
                    if chrg > 0 {
                        num_o_plus += 1;
                    } else {
                        num_o_minus += 1;
                    }
                }
                EL_NUMBER_N => {
                    if chrg > 0 {
                        num_n_plus += 1;
                    } else {
                        num_n_minus += 1;
                    }
                    num_p_iv_plus += i32::from(
                        a.el_number != EL_NUMBER_N
                            && chrg == 1
                            && a.valence == 4
                            && a.chem_bonds_valence == 4,
                    );
                }
                _ => {}
            }
        } else if chrg == 0
            && i32::from(a.chem_bonds_valence) + a.numh() == 2
            && get_el_valence(i32::from(a.el_number), 0, 0) == 4
            && ion_el_group(a.el_number) == EL_NUMBER_C
        {
            num_c_ii += 1;
        }
    }
    let mut num_all =
        num_c_ii + num_c_plus + num_c_minus + num_n_plus + num_n_minus + num_o_plus + num_o_minus;
    if num_all == 0 {
        return 0;
    }

    for type_ in 1..=18 {
        // ---- terminal ion pair type 1: O=N(+)-O(-) => O=N=O ----
        if type_ == 1 {
            let mut i = 0;
            while i < num_atoms && num_n_plus > 0 && num_o_minus > 0 {
                if at[i].charge == 1
                    && nnm_bonds(at, i) == 3
                    && nnm_val(at, i) == 4
                    && grp(at, i) == EL_NUMBER_N
                {
                    let mut ord_om = [0i32; 3];
                    let mut num_om = 0usize;
                    let mut num_o = 0;
                    let mut num_o_other = 0;
                    for i1 in 0..at[i].valence as usize {
                        let n = nb(at, i, i1);
                        if nnm_bonds(at, n) == 1
                            && num_of_h(at, n) == 0
                            && grp(at, n) == EL_NUMBER_O
                        {
                            if at[i].bond_type[i1] == SINGLE && at[n].charge == -1 {
                                if num_om < 3 {
                                    ord_om[num_om] = i1 as i32;
                                }
                                num_om += 1;
                            } else if at[n].bond_type[0] == DOUBLE && at[n].charge == 0 {
                                num_o += 1;
                            } else {
                                num_o_other += 1;
                            }
                        }
                    }
                    if num_om > 0 && num_o > 0 && num_o_other == 0 {
                        let i1 = n_find_one_om(at, i, &mut ord_om, num_om);
                        if i1 >= 0 {
                            // remove charges and increase bond order
                            let i1 = i1 as usize;
                            let n = nb(at, i, i1);
                            let i2 = pos_of(at, n, i);
                            neutralize_pair(at, i, i1, n, i2);
                            at[i].radical = 0;
                            at[n].radical = 0;
                            num_changes += 1;
                            num_n_plus -= 1;
                            num_o_minus -= 1;
                            num_all -= 2;
                        }
                    }
                }
                i += 1;
            }
        }

        // ---- terminal pair type 1a: Y3P(+)-O(-) => Y3P=O (FIX_P_IV_Plus_O_Minus) ----
        if type_ == 1 {
            let mut i = 0;
            while i < num_atoms && num_p_iv_plus > 0 && num_o_minus > 0 {
                if at[i].charge == 1
                    && nnm_bonds(at, i) == 4
                    && nnm_val(at, i) == 4
                    && at[i].el_number != EL_NUMBER_N
                    && grp(at, i) == EL_NUMBER_N
                {
                    let mut ord_om = [0i32; 4];
                    let mut num_om = 0usize;
                    for i1 in 0..at[i].valence as usize {
                        let n = nb(at, i, i1);
                        if nnm_bonds(at, n) == 1
                            && num_of_h(at, n) == 0
                            && grp(at, n) == EL_NUMBER_O
                            && at[i].bond_type[i1] == SINGLE
                            && at[n].charge == -1
                        {
                            ord_om[num_om] = i1 as i32;
                            num_om += 1;
                        }
                    }
                    if num_om > 0 {
                        let i1 = n_find_one_om(at, i, &mut ord_om, num_om);
                        if i1 >= 0 {
                            let i1 = i1 as usize;
                            let n = nb(at, i, i1);
                            let i2 = pos_of(at, n, i);
                            neutralize_pair(at, i, i1, n, i2);
                            at[i].radical = 0;
                            at[n].radical = 0;
                            num_changes += 1;
                            num_n_plus -= 1;
                            num_o_minus -= 1;
                            num_p_iv_plus -= 1;
                            num_all -= 2;
                        }
                    }
                }
                i += 1;
            }
        }

        // ---- terminal pair types 2..9 ----
        if (2..=9).contains(&type_) {
            let mut i = 0;
            while i < num_atoms && num_all > 0 {
                let i1c = n_no_metal_neigh_index(at, i);
                if at[i].charge == 0
                    && nnm_bonds(at, i) == 1
                    && nnm_val(at, i) == 2
                    && num_of_h(at, i) == 0
                    && grp(at, i) == EL_NUMBER_O
                    && i1c >= 0
                    && at[i].bond_type[i1c as usize] <= TRIPLE
                {
                    // terminal O=
                    let i1 = i1c as usize;
                    let n = nb(at, i, i1);
                    if type_ == 2 && num_c_ii > 0 {
                        // avoid alternating bonds
                        if at[n].charge == 0
                            && nnm_bonds(at, n) == 2
                            && nnm_val(at, n) == 3
                            && num_of_h(at, n) == 0
                            && grp(at, n) == EL_NUMBER_N
                        {
                            let i2c = n_no_metal_other_neigh_index(at, n, i as i32);
                            if i2c >= 0 && at[n].bond_type[i2c as usize] <= TRIPLE {
                                let i2 = i2c as usize;
                                let n2 = nb(at, n, i2);
                                if at[n2].charge == 0
                                    && at[n2].valence == 2
                                    && at[n2].chem_bonds_valence == 2
                                    && num_of_h(at, n2) == 0
                                    && grp(at, n2) == EL_NUMBER_C
                                {
                                    // O=N-C(II)- => O=N#C-
                                    let i3 = usize::from(at[n2].neighbor[0] as usize != n);
                                    at[n].chem_bonds_valence = 5;
                                    at[n2].chem_bonds_valence = 4;
                                    at[n].bond_type[i2] = TRIPLE;
                                    at[n2].bond_type[i3] = TRIPLE;
                                    at[n2].radical = 0;
                                    num_changes += 1;
                                    num_c_ii -= 1;
                                    num_all -= 1;
                                    i += 1;
                                    continue;
                                }
                            }
                        }
                    }
                    if type_ == 3 && num_o_plus > 0 && num_c_minus > 0 {
                        if at[n].charge == 1
                            && nnm_bonds(at, n) == 2
                            && nnm_val(at, n) == 3
                            && num_of_h(at, n) == 0
                            && grp(at, n) == EL_NUMBER_O
                        {
                            let i2c = n_no_metal_other_neigh_index(at, n, i as i32);
                            if i2c >= 0 && at[n].bond_type[i2c as usize] <= TRIPLE {
                                let i2 = i2c as usize;
                                let n2 = nb(at, n, i2);
                                if at[n2].charge == -1
                                    && nnm_bonds(at, n2) <= 3
                                    && nnm_val(at, n2) + at[n2].numh() == 3
                                    && grp(at, n2) == EL_NUMBER_C
                                {
                                    // O=O(+)-C(-)(III) => O=O=C(IV)
                                    let i3 = usize::from(at[n2].neighbor[0] as usize != n);
                                    at[n].charge -= 1;
                                    at[n2].charge += 1;
                                    at[n].chem_bonds_valence += 1;
                                    at[n2].chem_bonds_valence += 1;
                                    at[n].bond_type[i2] = DOUBLE;
                                    at[n2].bond_type[i3] = DOUBLE;
                                    num_changes += 1;
                                    num_o_plus -= 1;
                                    num_c_minus -= 1;
                                    num_all -= 2;
                                    i += 1;
                                    continue;
                                }
                            }
                        }
                    }
                } else if at[i].charge == -1
                    && num_o_minus + num_n_minus > 0
                    && num_n_plus + num_o_plus + num_c_plus > 0
                    && nnm_bonds(at, i) == 1
                    && nnm_val(at, i) == 1
                    && num_of_h(at, i) == 0
                    && grp(at, i) == EL_NUMBER_O
                    && i1c >= 0
                    && at[i].bond_type[i1c as usize] <= TRIPLE
                {
                    // terminal O(-)-
                    let i1 = i1c as usize;
                    let n = nb(at, i, i1);
                    if type_ == 4
                        && num_o_minus > 0
                        && num_n_plus > 0
                        && at[n].charge == 1
                        && nnm_bonds(at, n) <= 3
                        && nnm_val(at, n) == 4
                        && num_of_h(at, n) == 0
                        && grp(at, n) == EL_NUMBER_N
                    {
                        // O(-)-N(+)(IV) => O=N(V)
                        let i2 = pos_of(at, n, i);
                        at[i].charge += 1;
                        at[n].charge -= 1;
                        at[i].chem_bonds_valence += 1;
                        at[n].chem_bonds_valence += 1;
                        at[i].bond_type[i1] += 1;
                        at[n].bond_type[i2] += 1;
                        num_changes += 1;
                        num_o_minus -= 1;
                        num_n_plus -= 1;
                        num_all -= 2;
                        i += 1;
                        continue;
                    }
                    if type_ == 5
                        && num_o_minus > 0
                        && num_o_plus > 0
                        && at[n].charge == 1
                        && nnm_bonds(at, n) <= 3
                        && nnm_val(at, n) == 3
                        && num_of_h(at, n) == 0
                        && grp(at, n) == EL_NUMBER_O
                    {
                        // O(-)-O(+)(III) => O=O(IV)
                        let i2 = pos_of(at, n, i);
                        at[i].charge += 1;
                        at[n].charge -= 1;
                        at[i].chem_bonds_valence += 1;
                        at[n].chem_bonds_valence += 1;
                        at[i].bond_type[i1] += 1;
                        at[n].bond_type[i2] += 1;
                        num_changes += 1;
                        num_o_minus -= 1;
                        num_o_plus -= 1;
                        num_all -= 2;
                        i += 1;
                        continue;
                    }
                    if type_ == 6
                        && num_o_minus > 0
                        && num_c_plus > 0
                        && at[n].charge == 0
                        && nnm_bonds(at, n) == 2
                        && nnm_val(at, n) == 2
                        && num_of_h(at, n) == 0
                        && grp(at, n) == EL_NUMBER_O
                    {
                        let i2c = n_no_metal_other_neigh_index(at, n, i as i32);
                        if i2c >= 0 && at[n].bond_type[i2c as usize] <= TRIPLE {
                            let i2 = i2c as usize;
                            let n2 = nb(at, n, i2);
                            if at[n2].charge == 1
                                && nnm_bonds(at, n2) <= 3
                                && nnm_val(at, n2) + at[n2].numh() == 3
                                && grp(at, n2) == EL_NUMBER_C
                            {
                                // O(-)-O-C(+)(III) => O=O=C(IV)
                                let i3 = pos_of(at, n2, n);
                                let i4 = pos_of(at, n, i);
                                at[i].charge += 1;
                                at[n2].charge -= 1;
                                at[i].chem_bonds_valence += 1;
                                at[n].chem_bonds_valence += 2;
                                at[n2].chem_bonds_valence += 1;
                                at[i].bond_type[i1] = DOUBLE;
                                at[n].bond_type[i4] = DOUBLE;
                                at[n].bond_type[i2] = DOUBLE;
                                at[n2].bond_type[i3] = DOUBLE;
                                num_changes += 1;
                                num_o_minus -= 1;
                                num_c_plus -= 1;
                                num_all -= 2;
                                i += 1;
                                continue;
                            }
                        }
                    }
                } else if at[i].charge == -1
                    && num_n_minus > 0
                    && num_n_plus + num_o_plus + num_c_plus > 0
                    && nnm_bonds(at, i) == 1
                    && nnm_val(at, i) + at[i].numh() == 2
                    && grp(at, i) == EL_NUMBER_N
                    && i1c >= 0
                    && at[i].bond_type[i1c as usize] <= TRIPLE
                {
                    // terminal N(-)=
                    let i1 = 0usize;
                    let n = nb(at, i, i1);
                    let apply = |at: &mut [InpAtom], n: usize| {
                        let i2 = pos_of(at, n, i);
                        at[i].charge += 1;
                        at[n].charge -= 1;
                        at[i].chem_bonds_valence += 1;
                        at[n].chem_bonds_valence += 1;
                        at[i].bond_type[i1] += 1;
                        at[n].bond_type[i2] += 1;
                    };
                    if type_ == 7
                        && num_n_plus > 0
                        && at[n].charge == 1
                        && nnm_bonds(at, n) <= 3
                        && nnm_val(at, n) == 4
                        && num_of_h(at, n) == 0
                        && grp(at, n) == EL_NUMBER_N
                    {
                        // N(-)=N(+)(IV) => N#N(V)
                        apply(at, n);
                        num_changes += 1;
                        num_n_minus -= 1;
                        num_n_plus -= 1;
                        num_all -= 2;
                        i += 1;
                        continue;
                    }
                    if type_ == 8
                        && num_o_plus > 0
                        && at[n].charge == 1
                        && nnm_bonds(at, n) == 2
                        && nnm_val(at, n) == 3
                        && num_of_h(at, n) == 0
                        && grp(at, n) == EL_NUMBER_O
                    {
                        // N(-)=O(+)(III) => N#O(IV)-
                        apply(at, n);
                        num_changes += 1;
                        num_n_minus -= 1;
                        num_o_plus -= 1;
                        num_all -= 2;
                        i += 1;
                        continue;
                    }
                    if type_ == 9
                        && num_c_plus > 0
                        && at[n].charge == 1
                        && at[n].valence == 2
                        && at[n].chem_bonds_valence == 3
                        && num_of_h(at, n) == 0
                        && grp(at, n) == EL_NUMBER_C
                    {
                        // N(-)=C(+)(III) => N#C(IV)-
                        apply(at, n);
                        num_changes += 1;
                        num_n_minus -= 1;
                        num_c_plus -= 1;
                        num_all -= 2;
                        i += 1;
                        continue;
                    }
                }
                i += 1;
            }
        }

        // ---- non-terminal ion pairs: types 10..15 ----
        if (10..=15).contains(&type_) {
            let mut i = 0;
            while i < num_atoms && num_all > 0 {
                if at[i].charge == 1
                    && num_n_plus + num_o_plus > 0
                    && num_c_minus + num_n_minus > 0
                    && nnm_bonds(at, i) <= 4
                    && nnm_val(at, i) == 4
                    && num_of_h(at, i) == 0
                    && grp(at, i) == EL_NUMBER_N
                {
                    // non-terminal N(+)(IV)
                    if type_ == 10 && num_n_plus > 0 && num_c_minus > 0 {
                        let mut num_neigh = 0;
                        let mut pos_neigh = 0usize;
                        for i1 in 0..at[i].valence as usize {
                            let n = nb(at, i, i1);
                            if at[n].charge == -1
                                && at[n].valence <= 3
                                && i32::from(at[n].chem_bonds_valence) + at[n].numh() == 3
                                && at[i].bond_type[i1] == SINGLE
                                && grp(at, n) == EL_NUMBER_C
                            {
                                num_neigh += 1;
                                pos_neigh = i1;
                            }
                        }
                        if num_neigh == 1 {
                            let i1 = pos_neigh;
                            let n = nb(at, i, i1);
                            if at[i].bond_type[i1] <= TRIPLE
                                && !has_other_ion_neigh(at, i, n)
                                && !has_other_ion_neigh(at, n, i)
                            {
                                let i2 = pos_of(at, n, i);
                                neutralize_pair(at, i, i1, n, i2);
                                num_changes += 1;
                                num_c_minus -= 1;
                                num_n_plus -= 1;
                                num_all -= 2;
                                i += 1;
                                continue;
                            }
                        }
                    }
                    if type_ == 11 && num_n_plus > 0 && num_c_minus > 0 {
                        let mut num_neigh = 0;
                        let mut pos_neigh = 0usize;
                        for i1 in 0..at[i].valence as usize {
                            let n = nb(at, i, i1);
                            if at[n].charge == -1
                                && at[n].valence <= 3
                                && i32::from(at[n].chem_bonds_valence) + at[n].numh() == 3
                                && at[i].bond_type[i1] == DOUBLE
                                && grp(at, n) == EL_NUMBER_C
                            {
                                num_neigh += 1;
                                pos_neigh = i1;
                            }
                        }
                        if num_neigh == 1 {
                            let i1 = pos_neigh;
                            let n = nb(at, i, i1);
                            if !has_other_ion_neigh(at, i, n) && !has_other_ion_neigh(at, n, i) {
                                let i2 = pos_of(at, n, i);
                                neutralize_pair(at, i, i1, n, i2);
                                num_changes += 1;
                                num_c_minus -= 1;
                                num_n_plus -= 1;
                                num_all -= 2;
                                i += 1;
                                continue;
                            }
                        }
                    }
                    if type_ == 12 && num_n_plus > 0 && num_n_minus > 0 {
                        let mut num_neigh = 0;
                        let mut pos_neigh = 0usize;
                        for i1 in 0..at[i].valence as usize {
                            let n = nb(at, i, i1);
                            if at[n].charge == -1
                                && nnm_bonds(at, n) <= 2
                                && nnm_val(at, n) + at[n].numh() == 2
                                && at[i].bond_type[i1] == SINGLE
                                && grp(at, n) == EL_NUMBER_N
                            {
                                num_neigh += 1;
                                pos_neigh = i1;
                            }
                        }
                        if num_neigh == 1 {
                            let i1 = pos_neigh;
                            let n = nb(at, i, i1);
                            if !has_other_ion_neigh(at, i, n) && !has_other_ion_neigh(at, n, i) {
                                let i2 = pos_of(at, n, i);
                                neutralize_pair(at, i, i1, n, i2);
                                num_changes += 1;
                                num_n_minus -= 1;
                                num_n_plus -= 1;
                                num_all -= 2;
                                i += 1;
                                continue;
                            }
                        }
                    }
                } else if at[i].charge == 1
                    && num_o_plus > 0
                    && num_c_minus + num_n_minus > 0
                    && nnm_bonds(at, i) <= 3
                    && nnm_val(at, i) == 3
                    && num_of_h(at, i) == 0
                    && grp(at, i) == EL_NUMBER_O
                {
                    // non-terminal O(+)(III)
                    for (tt, bt, want_c) in
                        [(13, SINGLE, true), (14, DOUBLE, true), (15, SINGLE, false)]
                    {
                        if type_ != tt
                            || (want_c && num_c_minus == 0)
                            || (!want_c && num_n_minus == 0)
                        {
                            continue;
                        }
                        let mut num_neigh = 0;
                        let mut pos_neigh = 0usize;
                        for i1 in 0..at[i].valence as usize {
                            let n = nb(at, i, i1);
                            let ok = if want_c {
                                at[n].charge == -1
                                    && at[n].valence <= 3
                                    && i32::from(at[n].chem_bonds_valence) + at[n].numh() == 3
                                    && at[i].bond_type[i1] == bt
                                    && grp(at, n) == EL_NUMBER_C
                            } else {
                                at[n].charge == -1
                                    && nnm_bonds(at, n) <= 2
                                    && nnm_val(at, n) + at[n].numh() == 2
                                    && at[i].bond_type[i1] == bt
                                    && grp(at, n) == EL_NUMBER_N
                            };
                            if ok {
                                num_neigh += 1;
                                pos_neigh = i1;
                            }
                        }
                        if num_neigh == 1 {
                            let i1 = pos_neigh;
                            let n = nb(at, i, i1);
                            if !has_other_ion_neigh(at, i, n) && !has_other_ion_neigh(at, n, i) {
                                let i2 = pos_of(at, n, i);
                                neutralize_pair(at, i, i1, n, i2);
                                num_changes += 1;
                                if want_c {
                                    num_c_minus -= 1;
                                } else {
                                    num_n_minus -= 1;
                                }
                                num_o_plus -= 1;
                                num_all -= 2;
                                break;
                            }
                        }
                    }
                    // (the C code `continue`s after a change; the loop below
                    // advances in both cases, so the next atom is examined)
                }
                i += 1;
            }
        }

        // ---- non-terminal ion triples: types 16, 17, 18 ----
        if type_ == 16 && num_c_plus > 0 && num_n_minus > 0 {
            for i in 0..num_atoms {
                if !(at[i].charge == 0
                    && nnm_bonds(at, i) == 2
                    && nnm_val(at, i) == 2
                    && num_of_h(at, i) == 0)
                {
                    continue;
                }
                let j0 = n_no_metal_neigh_index(at, i);
                if j0 < 0 {
                    continue;
                }
                let m0 = nb(at, i, j0 as usize);
                if at[m0].charge == 0 {
                    continue;
                }
                let j1 = n_no_metal_other_neigh_index(at, i, m0 as i32);
                if j1 < 0 {
                    continue;
                }
                let m1 = nb(at, i, j1 as usize);
                if at[m0].charge + at[m1].charge != 0
                    || nnm_val(at, m0) + nnm_val(at, m1) > 5
                    || grp(at, i) != EL_NUMBER_O
                {
                    continue;
                }
                let m = [m0, m1];
                let j = [j0 as usize, j1 as usize];
                // found non-terminal A(+)-O-B(-); chem_bond_val of A+B <= 5
                let (mut n_n, mut n_c, mut i_c) = (-1i32, -1i32, 0usize);
                for k in 0..2 {
                    let n = m[k];
                    if at[n].charge == -1
                        && nnm_bonds(at, n) + at[n].numh() == 2
                        && grp(at, n) == EL_NUMBER_N
                    {
                        n_n = n as i32;
                    } else if at[n].charge == 1
                        && i32::from(at[n].chem_bonds_valence) + at[n].numh() == 3
                        && grp(at, n) == EL_NUMBER_C
                    {
                        n_c = n as i32;
                        i_c = k;
                    }
                }
                if n_c < 0
                    || n_n < 0
                    || has_other_ion_in_sphere_2(at, n_c as usize, n_n as usize) != 0
                    || has_other_ion_in_sphere_2(at, n_n as usize, n_c as usize) != 0
                {
                    continue;
                }
                // C(+)(III)-O-N(-)(II) => C(IV)=O=N(III)
                for k in 0..2 {
                    let n = if k == 1 { n_c as usize } else { n_n as usize };
                    let i1 = if k == 1 { j[i_c] } else { j[1 - i_c] };
                    let i2 = pos_of(at, n, i);
                    at[i].bond_type[i1] += 1;
                    at[n].bond_type[i2] += 1;
                    at[i].chem_bonds_valence += 1;
                    at[n].chem_bonds_valence += 1;
                    at[n].charge += if k == 1 { -1 } else { 1 };
                }
                num_changes += 1;
                num_n_minus -= 1;
                num_c_plus -= 1;
                num_all -= 2;
            }
        }

        if type_ == 17 && num_c_plus > 0 && num_c_minus > 0 {
            for i in 0..num_atoms {
                if !(at[i].charge == 0
                    && nnm_bonds(at, i) == 3
                    && nnm_val(at, i) == 3
                    && num_of_h(at, i) == 0)
                {
                    continue;
                }
                let j0 = n_no_metal_neigh_index(at, i);
                if j0 < 0 {
                    continue;
                }
                let m0 = nb(at, i, j0 as usize);
                let j1 = n_no_metal_other_neigh_index(at, i, m0 as i32);
                if j1 < 0 {
                    continue;
                }
                let m1 = nb(at, i, j1 as usize);
                let j2 = n_no_metal_other_neigh_index2(at, i, m0 as i32, m1 as i32);
                if j2 < 0 {
                    continue;
                }
                let m2 = nb(at, i, j2 as usize);
                let m = [m0, m1, m2];
                let j = [j0 as usize, j1 as usize, j2 as usize];
                let c = [at[m0].charge, at[m1].charge, at[m2].charge];
                let nzero = c.iter().filter(|&&x| x == 0).count();
                let three = |k: usize| -> i32 {
                    i32::from(
                        (if c[k] != 0 {
                            i32::from(at[m[k]].chem_bonds_valence) + at[m[k]].numh()
                        } else {
                            0
                        }) == 3,
                    )
                };
                if !(nzero == 1
                    && i32::from(c[0]) + i32::from(c[1]) + i32::from(c[2]) == 0
                    && three(0) + three(1) + three(2) == 2
                    && grp(at, i) == EL_NUMBER_N)
                {
                    continue;
                }
                let (mut n_cp, mut n_cm, mut i_cp, mut i_cm) = (-1i32, -1i32, 0usize, 0usize);
                for k in 0..3 {
                    if c[k] != 0 {
                        let n = m[k];
                        if at[n].charge == -1 && grp(at, n) == EL_NUMBER_C {
                            n_cm = n as i32;
                            i_cm = k;
                        } else if at[n].charge == 1 && grp(at, n) == EL_NUMBER_C {
                            n_cp = n as i32;
                            i_cp = k;
                        }
                    }
                }
                if n_cp < 0
                    || n_cm < 0
                    || has_other_ion_in_sphere_2(at, n_cp as usize, n_cm as usize) != 0
                    || has_other_ion_in_sphere_2(at, n_cm as usize, n_cp as usize) != 0
                {
                    continue;
                }
                // C(+)(III)-N-C(-)(III) => C(IV)=N=C(IV)
                for k in 0..2 {
                    let n = if k == 1 { n_cp as usize } else { n_cm as usize };
                    let i1 = if k == 1 { j[i_cp] } else { j[i_cm] };
                    let i2 = pos_of(at, n, i);
                    at[i].bond_type[i1] += 1;
                    at[n].bond_type[i2] += 1;
                    at[i].chem_bonds_valence += 1;
                    at[n].chem_bonds_valence += 1;
                    at[n].charge += if k == 1 { -1 } else { 1 };
                }
                num_changes += 1;
                num_c_minus -= 1;
                num_c_plus -= 1;
                num_all -= 2;
            }
        }

        if type_ == 18 && ((num_c_plus > 0 && num_c_minus > 0) || num_c_ii > 0) {
            for i in 0..num_atoms {
                if !(at[i].charge == 0
                    && nnm_bonds(at, i) == 2
                    && nnm_val(at, i) == 3
                    && num_of_h(at, i) == 0)
                {
                    continue;
                }
                let j0 = n_no_metal_neigh_index(at, i);
                if j0 < 0 {
                    continue;
                }
                let m0 = nb(at, i, j0 as usize);
                let j1 = n_no_metal_other_neigh_index(at, i, m0 as i32);
                if j1 < 0 {
                    continue;
                }
                let m1 = nb(at, i, j1 as usize);
                let v = [
                    i32::from(at[m0].chem_bonds_valence) + at[m0].numh(),
                    i32::from(at[m1].chem_bonds_valence) + at[m1].numh(),
                ];
                if !(at[m0].charge + at[m1].charge == 0
                    && v[0] + v[1] == 6
                    && (v[0] - v[1]).abs() <= 2
                    && grp(at, i) == EL_NUMBER_N
                    && grp(at, m0) == EL_NUMBER_C
                    && grp(at, m1) == EL_NUMBER_C)
                {
                    continue;
                }
                let m = [m0, m1];
                let j = [j0 as usize, j1 as usize];
                // C(-)(III)-N=C(+)(III) or C(IV)=N-C(II) => C(IV)=N#C(IV)
                let (mut n_cp, mut n_cm, mut i_cp, mut i_cm) = (-1i32, -1i32, 0usize, 0usize);
                for k in 0..2 {
                    let n = m[k];
                    if v[k] == 4 || (v[k] == 3 && at[i].bond_type[j[k]] == SINGLE) {
                        n_cm = n as i32;
                        i_cm = k;
                    } else if v[k] == 2 || (v[k] == 3 && at[i].bond_type[j[k]] == DOUBLE) {
                        n_cp = n as i32;
                        i_cp = k;
                    }
                }
                if n_cp < 0 || n_cm < 0 {
                    continue;
                }
                let (n_cp, n_cm) = (n_cp as usize, n_cm as usize);
                if i32::from(at[n_cp].valence) + at[n_cp].numh() != 2 {
                    continue; // guarantees at[n_cp].valence <= 2
                }
                if v[i_cp] == 2 || at[n_cp].charge == 0 {
                    if at[n_cp].valence == 2 {
                        // neighbor of at[n_cp] opposite to at[i]
                        let k = nb(at, n_cp, usize::from(at[n_cp].neighbor[0] as usize == i));
                        if grp(at, k) == EL_NUMBER_N {
                            continue;
                        }
                    }
                } else if at[n_cp].charge != 0 {
                    if has_other_ion_in_sphere_2(at, n_cp, n_cm) != 0
                        || has_other_ion_in_sphere_2(at, n_cm, n_cp) != 0
                    {
                        continue;
                    }
                } else {
                    continue; // unknown case
                }
                if at[n_cp].charge != 0 {
                    num_c_minus -= 1;
                    num_c_plus -= 1;
                    num_all -= 2;
                } else {
                    num_c_ii -= 1;
                    num_all -= 1;
                }
                for k in 0..2 {
                    let n = if k == 1 { n_cp } else { n_cm };
                    let i3 = if k == 1 { i_cp } else { i_cm };
                    let i1 = j[i3];
                    if v[i3] < 4 {
                        let delta = (4 - v[i3]) as u8;
                        let i2 = pos_of(at, n, i);
                        at[i].bond_type[i1] += delta;
                        at[n].bond_type[i2] += delta;
                        at[i].chem_bonds_valence += delta as i8;
                        at[n].chem_bonds_valence += delta as i8;
                        at[n].charge = 0;
                        at[n].radical = 0;
                    }
                }
                at[i].charge = 0;
                at[i].radical = 0;
                num_changes += 1;
            }
        }
    }
    num_changes
}

/// `fix_odd_things`: corrects structure-drawing oddities — H(-)-X, charge
/// separated `XHm(-)--Y==XHn(+)`, `Y++` with two `X-`, doublet pairs — and
/// then removes ion pairs (`REMOVE_ION_PAIRS_EARLY`). `b_fix_bug` is the
/// `bFixBug` flag, `b_fix_non_uniform_draw` the `bFixNonUniformDraw` flag.
pub fn fix_odd_things(at: &mut [InpAtom], b_fix_bug: bool, b_fix_non_uniform_draw: bool) -> i32 {
    const EN: [u8; 8] = [
        EL_NUMBER_N,
        EL_NUMBER_P,
        EL_NUMBER_AS,
        EL_NUMBER_SB,
        EL_NUMBER_O,
        EL_NUMBER_S,
        EL_NUMBER_SE,
        EL_NUMBER_TE,
    ];
    const FIRST_NEIGHB2: usize = 4;
    const FIRST_CENTER2: usize = 5;
    let rad_ok = |r: i8| r == 0 || i32::from(r) == RAD_SINGLET;
    let num_atoms = at.len();
    let mut num_changes = 0;
    let mut c: usize = 0;

    if b_fix_non_uniform_draw {
        fix_non_uniform_drawn_oxoanions(at, &mut num_changes);
        fix_non_uniform_drawn_amidiniums(at, &mut num_changes);
    }

    // H(-)-X -> H-X(-); H(+)-X -> H-X(+)
    for i1 in 0..num_atoms {
        if at[i1].valence == 1
            && at[i1].charge.abs() == 1
            && rad_ok(at[i1].radical)
            && at[i1].bond_type[0] == SINGLE
            && at[i1].el_number == EL_NUMBER_H
        {
            let i2 = nb(at, i1, 0);
            if at[i2].el_number != EL_NUMBER_H && at[i1].numh() == 0 && at[i2].numh() == 0 {
                at[i2].charge += at[i1].charge;
                at[i1].charge = 0;
            }
        }
    }

    // replace XHm(-)--Y==XHn(+) with XHm==Y--XHn
    for i1 in 0..num_atoms {
        if at[i1].charge != 1
            || (at[i1].radical != 0 && i32::from(at[i1].radical) != RAD_SINGLET)
            || at[i1].chem_bonds_valence == at[i1].valence
            || !EN.contains(&at[i1].el_number)
            || get_el_valence(i32::from(at[i1].el_number), i32::from(at[i1].charge), 0)
                != i32::from(at[i1].chem_bonds_valence) + at[i1].numh()
        {
            continue;
        }
        // a candidate for X in XHn(+)
        if at[i1].valence == 1 && at[i1].bond_type[0] == DOUBLE {
            c = nb(at, i1, 0);
            for k2 in 0..at[c].valence as usize {
                let i2 = nb(at, c, k2);
                if at[i2].valence == 1
                    && at[i2].charge == -1
                    && at[i2].el_number == at[i1].el_number
                    && rad_ok(at[i2].radical)
                    && at[i2].bond_type[0] == SINGLE
                    && get_el_valence(i32::from(at[i2].el_number), i32::from(at[i2].charge), 0)
                        == i32::from(at[i2].chem_bonds_valence) + at[i2].numh()
                {
                    // found both X(-) and X(+); change bonds and remove charges
                    let k1 = pos_of(at, c, i1);
                    at[i1].charge = 0;
                    at[i2].charge = 0;
                    at[i1].bond_type[0] = SINGLE;
                    at[c].bond_type[k1] = SINGLE;
                    at[i1].chem_bonds_valence -= 1;
                    at[i2].bond_type[0] = DOUBLE;
                    at[c].bond_type[k2] = DOUBLE;
                    at[i2].chem_bonds_valence += 1;
                    num_changes += 1;
                    break;
                }
            }
        } else {
            // explicit H case: detect H-neighbors and Y
            let mut num_h_i1 = 0;
            let mut i1_c: i32 = -1;
            for ineigh in 0..at[i1].valence as usize {
                let neigh = nb(at, i1, ineigh);
                if at[neigh].el_number == EL_NUMBER_H {
                    if at[neigh].chem_bonds_valence == 1 && rad_ok(at[neigh].radical) {
                        num_h_i1 += 1;
                    } else {
                        break; // wrong neighbor
                    }
                } else if at[i1].bond_type[ineigh] == DOUBLE {
                    // a candidate for Y; the bond must be double
                    i1_c = ineigh as i32;
                    c = neigh;
                }
            }
            if i1_c < 0 || num_h_i1 + 1 != i32::from(at[i1].valence) {
                continue;
            }
            for k2 in 0..at[c].valence as usize {
                let i2 = nb(at, c, k2);
                if at[i2].charge == -1
                    && at[i2].el_number == at[i1].el_number
                    && rad_ok(at[i2].radical)
                    && get_el_valence(i32::from(at[i2].el_number), i32::from(at[i2].charge), 0)
                        == i32::from(at[i2].chem_bonds_valence) + at[i2].numh()
                {
                    let mut num_h_i2 = 0;
                    let mut i2_c: i32 = -1;
                    for ineigh in 0..at[i2].valence as usize {
                        let neigh = nb(at, i2, ineigh);
                        if at[neigh].el_number == EL_NUMBER_H {
                            if at[neigh].chem_bonds_valence == 1 && rad_ok(at[neigh].radical) {
                                num_h_i2 += 1;
                            } else {
                                break;
                            }
                        } else if c == neigh && at[i2].bond_type[ineigh] == SINGLE {
                            i2_c = ineigh as i32; // position of the Y neighbor
                        } else {
                            break;
                        }
                    }
                    if num_h_i2 + i32::from(i2_c >= 0) != i32::from(at[i2].valence) {
                        continue;
                    }
                    // found both X(-) and X(+); change bonds and remove charges
                    let k1 = pos_of(at, c, i1);
                    at[i1].charge = 0;
                    at[i2].charge = 0;
                    at[i1].bond_type[i1_c as usize] = SINGLE;
                    at[c].bond_type[k1] = SINGLE;
                    at[i1].chem_bonds_valence -= 1;
                    at[i2].bond_type[i2_c as usize] = DOUBLE;
                    at[c].bond_type[k2] = DOUBLE;
                    at[i2].chem_bonds_valence += 1;
                    num_changes += 1;
                    break;
                }
            }
        }
    }

    // X- ... >Y++ ... X-  =>  X= ... >Y ... X
    for i1 in 0..num_atoms {
        if !(at[i1].valence == 1
            && at[i1].charge == -1
            && rad_ok(at[i1].radical)
            && at[i1].numh() == 0
            && at[i1].bond_type[0] == SINGLE
            && EN[FIRST_NEIGHB2..].contains(&at[i1].el_number))
        {
            continue;
        }
        // a candidate for X
        c = nb(at, i1, 0); // candidate for Y
        let mut charge = 2;
        let accept = if (charge == i32::from(at[c].charge)
            && EN[FIRST_CENTER2..].contains(&at[c].el_number))
            || {
                charge = 1;
                charge == i32::from(at[c].charge) && at[c].el_number == EL_NUMBER_P
            } {
            at[c].valence == 4
                && rad_ok(at[c].radical)
                && at[c].valence == at[c].chem_bonds_valence
                && at[c].numh() == 0
        } else {
            false
        };
        if !accept {
            continue;
        }
        for k2 in 0..at[c].valence as usize {
            let i2 = nb(at, c, k2);
            if i2 == i1 {
                continue;
            }
            if at[i2].valence == 1
                && at[i2].charge == -1
                && EN[FIRST_NEIGHB2..].contains(&at[i2].el_number)
                && rad_ok(at[i2].radical)
                && at[i2].numh() == 0
                && at[i2].bond_type[0] == SINGLE
            {
                // found both X(-) and X(-); change bonds and remove charges
                let k1 = pos_of(at, c, i1);
                for i in 0..charge {
                    match i {
                        0 => {
                            at[i1].charge += 1;
                            at[i1].bond_type[0] = DOUBLE;
                            at[c].bond_type[k1] = DOUBLE;
                            at[i1].bond_stereo[0] = 0;
                            at[c].bond_stereo[k1] = 0;
                            at[i1].chem_bonds_valence += 1;
                            at[c].chem_bonds_valence += 1;
                            if b_fix_bug {
                                at[c].charge -= 1;
                            }
                            num_changes += 1;
                        }
                        1 => {
                            at[i2].charge += 1;
                            at[i2].bond_type[0] = DOUBLE;
                            at[c].bond_type[k2] = DOUBLE;
                            at[i2].bond_stereo[0] = 0;
                            at[c].bond_stereo[k2] = 0;
                            at[i2].chem_bonds_valence += 1;
                            at[c].chem_bonds_valence += 1;
                            if b_fix_bug {
                                at[c].charge -= 1;
                            }
                            num_changes += 1;
                        }
                        _ => {}
                    }
                }
                break;
            }
        }
    }

    // A(doublet)-B(doublet) -> A=B;  A(doublet)=B(doublet) -> A#B
    for i1 in 0..num_atoms {
        if i32::from(at[i1].radical) == RADICAL_DOUBLET
            && let Some((i2, k1, k2)) = the_only_doublet_neigh(at, i1)
            && i32::from(at[i1].bond_type[k1]) <= BOND_TYPE_DOUBLE
        {
            at[i1].bond_type[k1] += 1;
            at[i1].chem_bonds_valence += 1;
            at[i2].bond_type[k2] += 1;
            at[i2].chem_bonds_valence += 1;
            at[i1].radical = 0;
            at[i2].radical = 0;
        }
    }

    // REMOVE_ION_PAIRS_EARLY == 1
    num_changes += remove_ion_pairs(at);
    num_changes
}

// ---- bond removal for salt and metal disconnection ----

use super::ichister::{AB_PARITY_UNDF, atom_parity_well_def, get_opposite_sb_atom};
use super::inpdef::{MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS};

/// Removes the parity slot `m` of `at` (the following slots move up)
fn remove_sb_slot(a: &mut InpAtom, m: usize) {
    for q in m..MAX_NUM_STEREO_BONDS - 1 {
        a.sb_parity[q] = a.sb_parity[q + 1];
        a.sb_ord[q] = a.sb_ord[q + 1];
        a.sn_ord[q] = a.sn_ord[q + 1];
        a.sn_orig_at_num[q] = a.sn_orig_at_num[q + 1];
    }
    let l = MAX_NUM_STEREO_BONDS - 1;
    a.sb_parity[l] = 0;
    a.sb_ord[l] = 0;
    a.sn_ord[l] = 0;
    a.sn_orig_at_num[l] = 0;
}

/// `RemoveInpAtBond`: removes the bond `k` of atom `iat` (one direction only)
/// and updates the CML stereo of that atom
pub fn remove_inp_at_bond(atom: &mut [InpAtom], iat: usize, k: usize) -> i32 {
    let val = i32::from(atom[iat].valence) - 1;
    if val < 0 {
        return 0;
    }
    let mut bond = atom[iat].bond_type[k];
    if i32::from(bond) > BOND_TYPE_TRIPLE {
        bond = SINGLE;
    }
    let ki = k as i32;

    if atom[iat].p_parity != 0 {
        let own = atom[iat].orig_at_number;
        if atom[iat].p_orig_at_num[..MAX_NUM_STEREO_ATOM_NEIGH].contains(&own) {
            atom[iat].p_parity = 0; // only 3 bonds are present
        }
        if atom[iat].p_parity != 0 {
            let nk_orig = atom[atom[iat].neighbor[k] as usize].orig_at_number;
            let valence = (atom[iat].valence as usize).min(MAX_NUM_STEREO_ATOM_NEIGH);
            let m = (0..valence).find(|&m| nk_orig == atom[iat].p_orig_at_num[m]);
            match m {
                Some(m) => atom[iat].p_orig_at_num[m] = own,
                None => atom[iat].p_parity = 0,
            }
        }
    }

    if atom[iat].sb_parity[0] != 0 {
        let mut m = 0;
        while m < MAX_NUM_STEREO_BONDS && atom[iat].sb_parity[m] != 0 {
            let a = &atom[iat];
            let (sb_ord, sn_ord) = (i32::from(a.sb_ord[m]), i32::from(a.sn_ord[m]));
            if ki == sb_ord
                || (ki == sn_ord && val < 2 && atom_parity_well_def(i32::from(a.sb_parity[m])))
            {
                // stereogenic bond is removed or its only neighbor is removed
                let (len, pnxt, _, pord) = get_opposite_sb_atom(atom, iat, sb_ord as usize);
                let i = if len != 0 { pord } else { MAX_NUM_STEREO_BONDS };
                if i < MAX_NUM_STEREO_BONDS && atom[pnxt].sb_parity[i] != 0 {
                    remove_sb_slot(&mut atom[iat], m);
                    remove_sb_slot(&mut atom[pnxt], i);
                } else {
                    m += 1; // inconsistent stereobond parity
                }
            } else if ki == sn_ord {
                let alt = (0..atom[iat].valence as i32).find(|&j| j != ki && j != sb_ord);
                match alt {
                    None => {
                        if atom_parity_well_def(i32::from(atom[iat].sb_parity[m])) {
                            let (len, pnxt, _, pord) =
                                get_opposite_sb_atom(atom, iat, sb_ord as usize);
                            if len > 0 {
                                atom[iat].sb_parity[m] = AB_PARITY_UNDF as i8;
                                atom[pnxt].sb_parity[pord] = AB_PARITY_UNDF as i8;
                            }
                        }
                        let a = &mut atom[iat];
                        a.sn_ord[m] = -99;
                        a.sb_ord[m] -= i8::from(sb_ord > ki);
                        a.sn_orig_at_num[m] = 0;
                    }
                    Some(i) => {
                        let ni = atom[iat].neighbor[i as usize] as usize;
                        let orig = atom[ni].orig_at_number;
                        let a = &mut atom[iat];
                        if atom_parity_well_def(i32::from(a.sb_parity[m])) {
                            a.sb_parity[m] = 3 - a.sb_parity[m];
                        }
                        a.sn_ord[m] = (i - i32::from(i > ki)) as i8;
                        a.sb_ord[m] -= i8::from(sb_ord > ki);
                        a.sn_orig_at_num[m] = orig;
                    }
                }
                m += 1;
            } else {
                let a = &mut atom[iat];
                if sb_ord > ki {
                    a.sb_ord[m] -= 1;
                }
                if sn_ord > ki {
                    a.sn_ord[m] -= 1;
                }
                m += 1;
            }
        }
    }

    let a = &mut atom[iat];
    let v = val as usize;
    for q in k..v {
        a.neighbor[q] = a.neighbor[q + 1];
        a.bond_stereo[q] = a.bond_stereo[q + 1];
        a.bond_type[q] = a.bond_type[q + 1];
    }
    a.neighbor[v] = 0;
    a.bond_stereo[v] = 0;
    a.bond_type[v] = 0;
    a.valence = val as i8;
    a.chem_bonds_valence -= bond as i8;
    1
}

/// `DisconnectInpAtBond`: removes the bond `neigh_ord` of `iat` in both
/// directions; resets the old component numbers of both ends
pub fn disconnect_inp_at_bond(
    at: &mut [InpAtom],
    n_old_comp_number: Option<&mut [u16]>,
    iat: usize,
    neigh_ord: usize,
) -> bool {
    let neigh = at[iat].neighbor[neigh_ord] as usize;
    let Some(i) = at[neigh].neighbor[..at[neigh].valence as usize]
        .iter()
        .position(|&x| x as usize == iat)
    else {
        return false;
    };
    let ret = remove_inp_at_bond(at, iat, neigh_ord) + remove_inp_at_bond(at, neigh, i);
    if let Some(old) = n_old_comp_number
        && ret != 0
    {
        for x in [iat, neigh] {
            let component = at[x].component;
            if component != 0 {
                old[component as usize - 1] = 0;
            }
        }
    }
    ret == 2
}

// ---- ammonium and metal salt disconnection ----

use super::incomdef::IS_METAL;
use super::incomdef::NUM_H_ISOTOPES;
use super::util::{EL_NUMBER_BR, EL_NUMBER_CL, EL_NUMBER_F, EL_NUMBER_I, get_el_type};

/// `at[i].num_iso_H[idx]`; as in the C code the index `NUM_H_ISOTOPES`
/// addresses the next structure member, `iso_atw_diff`
fn iso_h_slot(a: &mut InpAtom, idx: usize) -> &mut i8 {
    if idx < NUM_H_ISOTOPES {
        &mut a.num_iso_h[idx]
    } else {
        &mut a.iso_atw_diff
    }
}

fn is_halogen(el: u8) -> bool {
    el == EL_NUMBER_F || el == EL_NUMBER_CL || el == EL_NUMBER_BR || el == EL_NUMBER_I
}

/// `bIsAmmoniumSalt`: `NH4(+)-O(-)-C` or `NH4-X`; returns the position of
/// the O (or halogen) atom `iO`, its index `k` in the neighbor list of `i`
/// and fills `num_explicit_h` with the explicit hydrogens by isotope
pub fn is_ammonium_salt(
    at: &[InpAtom],
    i: usize,
    num_explicit_h: &mut [i8; NUM_H_ISOTOPES + 1],
) -> Option<(usize, usize)> {
    if at[i].el_number != EL_NUMBER_N {
        return None;
    }
    let val = at[i].valence as usize;
    let mut num_h = i32::from(at[i].num_h) + at[i].num_iso_h_total();
    if val as i32 + num_h != 5 {
        return None;
    }
    let mut num_o = 0;
    let (mut io, mut k) = (0usize, 0usize);
    *num_explicit_h = [0; NUM_H_ISOTOPES + 1];
    for j in 0..val {
        let neigh = nb(at, i, j);
        let n = &at[neigh];
        if n.num_h != 0
            || (n.charge != 0
                && (n.el_number != EL_NUMBER_O
                    || i32::from(n.charge) + i32::from(at[i].charge) != 0))
            || (n.radical != 0 && i32::from(n.radical) != RAD_SINGLET)
        {
            return None;
        }
        if n.el_number == EL_NUMBER_H && n.valence == 1 && n.charge == 0 && n.radical == 0 {
            num_h += 1;
            num_explicit_h[n.iso_atw_diff as usize] += 1;
        } else if n.el_number == EL_NUMBER_O && n.valence == 2 && num_o == 0 {
            num_o += 1;
            io = neigh;
            k = j;
            let ic = nb(at, io, usize::from(nb(at, io, 0) == i));
            if at[ic].el_number != EL_NUMBER_C
                || at[ic].charge != 0
                || (at[ic].radical != 0 && i32::from(at[ic].radical) != RAD_SINGLET)
            {
                return None;
            }
        } else if is_halogen(n.el_number)
            && n.valence == 1
            && n.chem_bonds_valence == 1
            && n.charge == 0
            && n.numh() == 0
            && num_o == 0
        {
            num_o += 1;
            io = neigh;
            k = j;
        } else {
            return None;
        }
    }
    (num_o == 1 && num_h == 4).then_some((io, k))
}

/// `DisconnectAmmoniumSalt`: removes the bond N-O and moves one hydrogen
/// from N to O
pub fn disconnect_ammonium_salt(
    at: &mut [InpAtom],
    i_n: usize,
    i_o: usize,
    k: usize,
    num_explicit_h: &[i8; NUM_H_ISOTOPES + 1],
) -> i32 {
    let mut n_move_h_iso_diff: i32 = -1;
    let mut val = at[i_n].valence as i32;
    if at[i_n].charge != 0 && at[i_n].charge + at[i_o].charge == 0 {
        at[i_n].charge = 0;
        at[i_o].charge = 0;
    }
    let neigh_pos = if at[i_o].valence == 2 {
        usize::from(nb(at, i_o, 1) == i_n)
    } else {
        0
    };
    remove_inp_at_bond(at, i_o, neigh_pos);
    remove_inp_at_bond(at, i_n, k);
    val -= 1;

    for iso_diff in 0..=NUM_H_ISOTOPES {
        if iso_diff == 0 {
            if at[i_n].num_h != 0 {
                at[i_n].num_h -= 1;
                at[i_o].num_h += 1;
                break;
            } else if num_explicit_h[0] != 0 {
                n_move_h_iso_diff = 0;
                break;
            }
        } else if *iso_h_slot(&mut at[i_n], iso_diff) != 0 {
            *iso_h_slot(&mut at[i_n], iso_diff) -= 1;
            *iso_h_slot(&mut at[i_o], iso_diff) += 1;
            break;
        } else if num_explicit_h[iso_diff] != 0 {
            n_move_h_iso_diff = iso_diff as i32;
            break;
        }
    }

    if n_move_h_iso_diff >= 0 {
        let mut min_dist2: f64 = -1.0;
        let (mut jh, mut ih) = (0usize, 0usize);
        for j in 0..val as usize {
            let neigh = nb(at, i_n, j);
            if at[neigh].el_number == EL_NUMBER_H
                && i32::from(at[neigh].iso_atw_diff) == n_move_h_iso_diff
            {
                let d2 = (at[neigh].x - at[i_o].x) * (at[neigh].x - at[i_o].x)
                    + (at[neigh].y - at[i_o].y) * (at[neigh].y - at[i_o].y)
                    + (at[neigh].z - at[i_o].z) * (at[neigh].z - at[i_o].z);
                if min_dist2 < 0.0 || min_dist2 > d2 {
                    min_dist2 = d2;
                    ih = neigh;
                    jh = j;
                }
            }
        }
        let pos = at[i_o].valence as usize;
        let bt = at[ih].bond_type[0];
        at[i_o].neighbor[pos] = ih as u16;
        at[i_o].bond_stereo[pos] = 0;
        at[i_o].bond_type[pos] = bt;
        at[i_o].chem_bonds_valence += bt as i8;
        at[i_o].valence += 1;
        at[ih].neighbor[0] = i_o as u16;
        at[ih].bond_stereo[0] = 0;
        remove_inp_at_bond(at, i_n, jh);
    }
    1
}

/// `bIsMetalSalt`: a metal with the normal valence whose neighbors are all
/// `-O-C(=)` (without H on C) or halogens
pub fn is_metal_salt(at: &[InpAtom], i: usize) -> bool {
    let val = i32::from(at[i].valence);
    if val == 0 {
        return false;
    }
    let ty = get_el_type(i32::from(at[i].el_number));
    if ty == 0 || ty & IS_METAL == 0 || at[i].num_h != 0 {
        return false;
    }
    let el = i32::from(at[i].el_number);
    let chg = i32::from(at[i].charge);
    let ok = (chg == 0
        && ((ty & 1 != 0 && val == get_el_valence(el, 0, 0))
            || (ty & 2 != 0 && val == get_el_valence(el, 0, 1))))
        || (chg > 0 && ty & 1 != 0 && val == get_el_valence(el, chg, 0));
    if !ok {
        return false;
    }
    for k in 0..at[i].valence as usize {
        let io = nb(at, i, k);
        let o = &at[io];
        if is_halogen(o.el_number)
            && o.valence == 1
            && o.chem_bonds_valence == 1
            && o.charge == 0
            && !(o.radical != 0 && i32::from(o.radical) != RAD_SINGLET)
            && o.numh() == 0
        {
            continue;
        }
        if o.el_number != EL_NUMBER_O
            || o.numh() != 0
            || o.valence != 2
            || o.charge != 0
            || (o.radical != 0 && i32::from(o.radical) != RAD_SINGLET)
            || o.valence != o.chem_bonds_valence
        {
            return false;
        }
        let ic = nb(at, io, usize::from(nb(at, io, 0) == i));
        let c = &at[ic];
        if c.el_number != EL_NUMBER_C
            || c.num_h != 0
            || c.chem_bonds_valence != 4
            || c.charge != 0
            || (c.radical != 0 && i32::from(c.radical) != RAD_SINGLET)
            || c.valence == c.chem_bonds_valence
        {
            return false;
        }
        if (0..c.valence as usize).any(|j| at[c.neighbor[j] as usize].el_number == EL_NUMBER_H) {
            return false;
        }
    }
    true
}

/// `DisconnectMetalSalt`: cuts all bonds of the metal `i`, the neighbors
/// become anions; returns the number of removed bonds
pub fn disconnect_metal_salt(at: &mut [InpAtom], i: usize) -> i32 {
    let n = at[i].valence as usize;
    for k in 0..n {
        let io = nb(at, i, k);
        if at[io].valence == 2 {
            if at[io].neighbor[0] as usize == i {
                at[io].neighbor[0] = at[io].neighbor[1];
                at[io].bond_stereo[0] = at[io].bond_stereo[1];
                at[io].bond_type[0] = at[io].bond_type[1];
            }
            at[io].neighbor[1] = 0;
            at[io].bond_stereo[1] = 0;
            at[io].bond_type[1] = 0;
        } else {
            at[io].neighbor[0] = 0;
            at[io].bond_stereo[0] = 0;
            at[io].bond_type[0] = 0;
        }
        at[io].charge = -1;
        at[io].valence -= 1;
        at[io].chem_bonds_valence -= 1;
        at[i].neighbor[k] = 0;
        at[i].bond_stereo[k] = 0;
        at[i].bond_type[k] = 0;
        at[i].charge += 1;
    }
    at[i].valence = 0;
    at[i].chem_bonds_valence = 0;
    n as i32
}

/// `DisconnectSalts`: `num_inp_bonds` is `orig_inp_data->num_inp_bonds`
pub fn disconnect_salts(at: &mut [InpAtom], num_inp_bonds: &mut i32, b_disconnect: bool) -> i32 {
    let mut num_changes = 0;
    let mut num_explicit_h = [0i8; NUM_H_ISOTOPES + 1];
    for i in 0..at.len() {
        let val = at[i].valence;
        if val == 0
            || val != at[i].chem_bonds_valence
            || (at[i].radical != 0 && i32::from(at[i].radical) != RAD_SINGLET)
        {
            continue;
        }
        if let Some((io, k)) = is_ammonium_salt(at, i, &mut num_explicit_h) {
            if b_disconnect {
                disconnect_ammonium_salt(at, i, io, k, &num_explicit_h);
                *num_inp_bonds -= 1;
            }
            num_changes += 1;
        } else if is_metal_salt(at, i) {
            if b_disconnect {
                *num_inp_bonds -= disconnect_metal_salt(at, i);
            }
            num_changes += 1;
        }
    }
    num_changes
}

/// `bIsMetalToDisconnect`: 0 = nothing to disconnect, 1 = disconnect
/// (the C code's check for a normal valence, result 2, never runs: its loop
/// condition `i & type` is false for i = 0)
pub fn is_metal_to_disconnect(at: &[InpAtom], i: usize) -> i32 {
    let ty = get_el_type(i32::from(at[i].el_number));
    if ty == 0 || ty & IS_METAL == 0 {
        return 0;
    }
    let at_valence = at[i].numh() + i32::from(at[i].chem_bonds_valence);
    i32::from(at_valence != 0)
}

/// `bMayDisconnectMetals`: counts the metals that `DisconnectMetals` would
/// disconnect; returns that count and the `bDisconnectCoord` value
/// (implicit H count + 1, or 0)
pub fn may_disconnect_metals(at: &[InpAtom]) -> (i32, i32) {
    let mut num_changes = 0;
    let mut num_impl_h = 0;
    let mut num_explicit_h = [0i8; NUM_H_ISOTOPES + 1];
    for i in 0..at.len() {
        let val = at[i].valence;
        if val == 0 && at[i].numh() == 0 {
            continue;
        }
        let rad_or_mult = val == 0
            || val != at[i].chem_bonds_valence
            || (at[i].radical != 0 && i32::from(at[i].radical) != RAD_SINGLET);
        if !rad_or_mult && is_ammonium_salt(at, i, &mut num_explicit_h).is_some() {
        } else if !rad_or_mult && is_metal_salt(at, i) {
        } else if is_metal_to_disconnect(at, i) == 1 {
            num_impl_h += at[i].numh();
            num_changes += 1;
        }
    }
    (
        num_changes,
        if num_changes != 0 { num_impl_h + 1 } else { 0 },
    )
}

// ---- metal disconnection ----

/// `TG_FLAG_MOVE_CHARGE_COORD_DONE`
pub const TG_FLAG_MOVE_CHARGE_COORD_DONE: u64 = 0x0000_0400;

/// The members of `ORIG_ATOM_DATA` that the metal disconnection uses
#[derive(Debug, Clone, Default)]
pub struct OrigAtomData {
    pub at: Vec<InpAtom>,
    pub num_inp_atoms: usize,
    pub num_inp_bonds: i32,
    /// `bDisconnectCoord`: number of implicit H on the metals + 1, or 0
    pub b_disconnect_coord: i32,
    /// `nOldCompNumber`
    pub n_old_comp_number: Option<Vec<u16>>,
    /// `nCurAtLen`: number of atoms per component
    pub n_cur_at_len: Option<Vec<u16>>,
    /// `num_components`
    pub num_components: i32,
    /// `bDisconnectSalts`
    pub b_disconnect_salts: i32,
    /// `szCoord`: the original coordinate strings
    pub sz_coord: Option<Vec<super::ichiprt1::MolCoord>>,
    /// `num_dimensions`
    pub num_dimensions: i32,
    /// `polymer` (v. 1.05 extension)
    pub polymer: Option<super::mol2atom::OadPolymer>,
    /// `v3000` (v. 1.05 extension)
    pub v3000: Option<super::mol2atom::OadV3000>,
    /// `n_zy`: number of non-polymer-related pseudoelement atoms
    pub n_zy: i32,
    /// `valid_polymer`
    pub valid_polymer: i32,
}

/// Heteroatoms that may take a charge after a ligand is disconnected: the
/// first `NUM_HALOGENS` are halogens (`elnumber_Heteroat`)
const ELNUMBER_HETEROATOM: [u8; 13] = [
    EL_NUMBER_F,
    EL_NUMBER_CL,
    EL_NUMBER_BR,
    EL_NUMBER_I,
    EL_NUMBER_AT,
    EL_NUMBER_N,
    EL_NUMBER_P,
    EL_NUMBER_AS,
    EL_NUMBER_O,
    EL_NUMBER_S,
    EL_NUMBER_SE,
    EL_NUMBER_TE,
    EL_NUMBER_B,
];
const NUM_HALOGENS: usize = 5;

use super::inpdef::{ATOM_EL_LEN, MAX_ATOMS, MAXVAL};
use super::util::{EL_NUMBER_AT, EL_NUMBER_B};

/// `DisconnectOneLigand`: cuts the bonds between the ligand
/// `at[iMetal].neighbor[jLigand]` and metal atoms, then tries to give the
/// ligand a charge that fits its valence
pub fn disconnect_one_ligand(
    at: &mut [InpAtom],
    mut n_old_comp_number: Option<&mut [u16]>,
    b_metal: &[u8],
    num_atoms: usize,
    i_metal: usize,
    j_ligand: usize,
    taut_flags_done: Option<&mut u64>,
) -> i32 {
    let i_ligand = nb(at, i_metal, j_ligand);
    let mut metal_neigh_ord: Vec<usize> = Vec::new();
    let mut num_neigh_arom_bonds = [0i32; MAXVAL];
    let mut num_disconnections = 0;
    let (mut num_del_arom_bonds, mut num_tot_arom_bonds) = (0i32, 0i32);

    for i in 0..at[i_ligand].valence as usize {
        let neigh = nb(at, i_ligand, i);
        if neigh < num_atoms && b_metal[neigh] != 0 {
            metal_neigh_ord.push(i);
            if i32::from(at[i_ligand].bond_type[i]) > BOND_TYPE_TRIPLE {
                for j in 0..at[neigh].valence as usize {
                    num_neigh_arom_bonds[i] +=
                        i32::from(i32::from(at[neigh].bond_type[j]) > BOND_TYPE_TRIPLE);
                }
                num_del_arom_bonds += 1;
            }
        }
        num_tot_arom_bonds += i32::from(i32::from(at[i_ligand].bond_type[i]) > BOND_TYPE_TRIPLE);
    }

    if num_del_arom_bonds != 0 {
        for &j in &metal_neigh_ord {
            if num_neigh_arom_bonds[j] != 0 {
                let neigh = nb(at, i_ligand, j);
                let d = num_neigh_arom_bonds[j] / 2 - (num_neigh_arom_bonds[j] - 1) / 2;
                at[neigh].chem_bonds_valence -= d as i8;
            }
        }
        let d = num_tot_arom_bonds / 2 - (num_tot_arom_bonds - num_del_arom_bonds) / 2;
        at[i_ligand].chem_bonds_valence -= d as i8;
    }

    for &ord in metal_neigh_ord.iter().rev() {
        num_disconnections += i32::from(disconnect_inp_at_bond(
            at,
            n_old_comp_number.as_deref_mut(),
            i_ligand,
            ord,
        ));
    }

    // attempt to change the ligand charge to make its valence natural
    let i = num_tot_arom_bonds - num_del_arom_bonds;
    let p = ELNUMBER_HETEROATOM
        .iter()
        .position(|&e| e == at[i_ligand].el_number)
        .or(if at[i_ligand].el_number == 0 {
            Some(ELNUMBER_HETEROATOM.len())
        } else {
            None
        });
    if (i != 0 && i != 2 && i != 3)
        || (at[i_ligand].radical != 0 && i32::from(at[i_ligand].radical) != RAD_SINGLET)
        || p.is_none()
    {
        return num_disconnections;
    }
    let p = p.unwrap();
    let val = i32::from(at[i_ligand].chem_bonds_valence) + at[i_ligand].numh();
    let mut new_charge = MAX_ATOMS as i32; // impossible value
    if val == 0 {
        if p < NUM_HALOGENS {
            new_charge = -1;
        }
    } else {
        for c in -1..=1 {
            if val == get_el_valence(i32::from(at[i_ligand].el_number), c, 0) {
                new_charge = c;
                break;
            }
        }
    }
    if new_charge != MAX_ATOMS as i32
        && (new_charge != i32::from(at[i_ligand].charge)
            || (at[i_ligand].radical != 0 && i32::from(at[i_ligand].radical) != RAD_SINGLET))
        && metal_neigh_ord.len() == 1
    {
        let l = &at[i_ligand];
        if new_charge == 1
            && val == 4
            && l.valence == 2
            && l.chem_bonds_valence == 4
            && l.bond_type[0] == l.bond_type[1]
        {
            // do not add +1 charge to disconnected =N=
        } else {
            if let Some(f) = taut_flags_done
                && new_charge != i32::from(l.charge)
            {
                *f |= TG_FLAG_MOVE_CHARGE_COORD_DONE;
            }
            let old = i32::from(at[i_ligand].charge);
            at[i_metal].charge -= (new_charge - old) as i8;
            at[i_ligand].charge = new_charge as i8;
        }
    }
    num_disconnections
}

/// `dist3D`
pub fn dist3d(a: &InpAtom, b: &InpAtom) -> f64 {
    let (dx, dy, dz) = (a.x - b.x, a.y - b.y, a.z - b.z);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

const MIN_BOND_LENGTH: f64 = 1.0e-6;
const MIN_COS: f64 = 1.0e-6;
const MIN_BOND_LENGTH2: f64 = MIN_BOND_LENGTH * MIN_BOND_LENGTH;
const MAX_BOND_LENGTH: f64 = 1.0e30;

/// `GetMinDistDistribution`: distances from `at[iat]` to the nearest bonds in
/// `num_segm` angular sectors; returns the average bond length
pub fn get_min_dist_distribution(
    at: &[InpAtom],
    num_at: usize,
    iat: usize,
    iat_h: usize,
    b_in_all_components: bool,
    min_dist: &mut [f64],
    num_segm: usize,
) -> f64 {
    let one_pi = std::f64::consts::PI;
    let two_pi = 2.0 * one_pi;
    let f_step = two_pi / num_segm as f64;
    let h_step = f_step / 2.0;
    let ns = num_segm as i32;
    for d in min_dist.iter_mut().take(num_segm) {
        *d = MAX_BOND_LENGTH;
    }
    let mut num_bonds = 0;
    let mut ave_bond_len = 0.0;
    let seg = |f: f64| ((f + h_step) / f_step).floor() as i32;

    for i in 0..num_at {
        if i == iat || i == iat_h || !(b_in_all_components || at[i].component == at[iat].component)
        {
            continue;
        }
        for j in 0..at[i].valence as usize {
            let n = nb(at, i, j);
            if (n > i && n != iat) || n == iat_h {
                continue;
            }
            let mut xi = at[i].x - at[iat].x;
            let mut yi = at[i].y - at[iat].y;
            let mut xn = at[n].x - at[iat].x;
            let mut yn = at[n].y - at[iat].y;
            let cross_prod_in = xi * yn - xn * yi;
            if cross_prod_in < -0.01 * MIN_BOND_LENGTH2 {
                std::mem::swap(&mut xi, &mut xn);
                std::mem::swap(&mut yi, &mut yn);
            }
            let xni = xn - xi;
            let yni = yn - yi;
            let rni = xni * xni + yni * yni;
            let tni;
            let rmin;
            if rni > 0.01 * MIN_BOND_LENGTH2 {
                tni = -(xni * xi + yni * yi) / rni;
                rmin = if tni < 0.0 {
                    (xi * xi + yi * yi).sqrt()
                } else if tni > 1.0 {
                    (xn * xn + yn * yn).sqrt()
                } else {
                    (tni * tni * rni).sqrt()
                };
                ave_bond_len += rni.sqrt();
                num_bonds += 1;
            } else {
                tni = 0.5;
                rmin = (xi * xi + yi * yi).sqrt();
            }
            if rmin >= 0.1 * MIN_BOND_LENGTH {
                let mut b_calc_rt = true;
                let (mut ft, mut rt) = (0.0, 0.0);
                let mut fi = yi.atan2(xi);
                let mut fn_ = if n == iat { fi } else { yn.atan2(xn) };
                if fi > fn_ {
                    fn_ += two_pi;
                }
                if fi < 0.0 {
                    fi += two_pi;
                    fn_ += two_pi;
                }
                let ki = seg(fi);
                let kn = seg(fn_);
                for k in ki..=kn {
                    let kk = (k % ns) as usize;
                    if min_dist[kk] < rmin {
                        continue;
                    }
                    if b_calc_rt {
                        if n == iat {
                            ft = fi;
                            rt = rmin;
                        } else {
                            let xt = xi + xni * tni;
                            let yt = yi + yni * tni;
                            ft = yt.atan2(xt);
                            rt = (xt * xt + yt * yt).sqrt();
                        }
                        b_calc_rt = false;
                    }
                    let fk = f_step * kk as f64;
                    let mut c = (fk - ft).cos().abs();
                    if c < MIN_COS {
                        c = MIN_COS;
                    }
                    let rk = rt / c;
                    if min_dist[kk] > rk {
                        min_dist[kk] = rk;
                    }
                }
            } else {
                let ri = xi * xi + yi * yi;
                let rn = xn * xn + yn * yn;
                if ri > MIN_BOND_LENGTH2 && rn > MIN_BOND_LENGTH2 {
                    let dot_prod_in = xn * xi + yn * yi;
                    if dot_prod_in > 0.01 * MIN_BOND_LENGTH2 {
                        let mut f = yi.atan2(xi);
                        if f < 0.0 {
                            f += two_pi;
                        }
                        let kk = (seg(f) % ns) as usize;
                        if min_dist[kk] > rmin {
                            min_dist[kk] = rmin;
                        }
                    } else if dot_prod_in < -0.01 * MIN_BOND_LENGTH2 {
                        let mut f = yi.atan2(xi);
                        if f < 0.0 {
                            f += two_pi;
                        }
                        let kk = (seg(f) % ns) as usize;
                        if min_dist[kk] > rmin {
                            min_dist[kk] = rmin;
                        }
                        f += one_pi;
                        let kk = (seg(f) % ns) as usize;
                        if min_dist[kk] > rmin {
                            min_dist[kk] = rmin;
                        }
                    }
                } else if ri <= MIN_BOND_LENGTH2 && rn <= MIN_BOND_LENGTH2 {
                    // coincides with at[iat]; ignore
                } else {
                    let mut f = if ri > rn { yi.atan2(xi) } else { yn.atan2(xn) };
                    if f < 0.0 {
                        f += two_pi;
                    }
                    let kk = (seg(f) % ns) as usize;
                    if min_dist[kk] > rmin {
                        min_dist[kk] = rmin;
                    }
                }
            }
        }
    }
    if num_bonds != 0 {
        ave_bond_len / f64::from(num_bonds)
    } else {
        0.0
    }
}

/// `move_explicit_Hcation`: attaches the hydrogen `iat_h` to `iat`, placing it
/// at a free position around `iat` in the 2D drawing
pub fn move_explicit_h_cation(
    at: &mut [InpAtom],
    num_at: usize,
    iat: usize,
    iat_h: usize,
    b_in_all_components: bool,
) -> i32 {
    const NUM_SEGM: usize = 20;
    let one_pi = std::f64::consts::PI;
    let two_pi = 2.0 * one_pi;
    let f_step = two_pi / NUM_SEGM as f64;
    let h_step = f_step / 2.0;
    let mut min_dist = [0f64; NUM_SEGM];
    let (mut xd, mut yd, mut zd) = (0.0f64, 0.0f64, 0.0f64);
    let mut r;
    let mut r0;

    if at[iat].valence != 0 {
        let mut nbn = 0.0;
        r = 0.0;
        for i in 0..at[iat].valence as usize {
            let next = nb(at, iat, i);
            xd += at[next].x;
            yd += at[next].y;
            zd += at[next].z;
            r += dist3d(&at[iat], &at[next]);
            nbn += 1.0;
        }
        xd /= nbn;
        yd /= nbn;
        zd /= nbn;
        r /= nbn;
        r0 = ((xd - at[iat].x) * (xd - at[iat].x) + (yd - at[iat].y) * (yd - at[iat].y)).sqrt();
    } else {
        r = if at[iat_h].valence != 0 {
            dist3d(&at[iat_h], &at[nb(at, iat_h, 0)])
        } else {
            0.0
        };
        r0 = 0.0;
    }

    let ave_bond_len = get_min_dist_distribution(
        at,
        num_at,
        iat,
        iat_h,
        b_in_all_components,
        &mut min_dist,
        NUM_SEGM,
    );
    if r < MIN_BOND_LENGTH && ave_bond_len > MIN_BOND_LENGTH {
        r = ave_bond_len;
    }

    let (mut xr, mut yr, mut zr);
    if r > MIN_BOND_LENGTH {
        if 10.0 * r0 < r {
            xr = -r;
            yr = 0.0;
            zr = 0.0;
        } else {
            xr = r * (at[iat].x - xd) / r0;
            yr = r * (at[iat].y - yd) / r0;
            zr = r * (at[iat].z - zd) / r0;
            if xr * xr + yr * yr < 0.04 * r * r {
                xr = -r;
                yr = 0.0;
            }
        }
        r = (xr * xr + yr * yr).sqrt();
        let mut f = yr.atan2(xr);
        if f < 0.0 {
            f += two_pi;
        }
        let kk = (((f + h_step) / f_step).floor() as i32 % NUM_SEGM as i32) as usize;
        if min_dist[kk] < 1.5 * r {
            let mut dist = 1.5 * r;
            let mut start_max: i32 = -1;
            loop {
                let (mut start, mut len, mut len_max) = (0i32, 0i32, 0i32);
                for k in 0..2 * NUM_SEGM as i32 {
                    let kk = (k % NUM_SEGM as i32) as usize;
                    if min_dist[kk] >= dist {
                        if len == 0 {
                            start = k;
                        }
                        len += 1;
                    } else {
                        if len > len_max {
                            len_max = len;
                            start_max = start;
                        }
                        len = 0;
                    }
                }
                if len_max == 0 {
                    if dist > 0.1 * r {
                        dist *= 0.75;
                        continue;
                    }
                    break; // do it anyway
                }
                f = f_step * (f64::from(start_max) + (f64::from(len_max) - 1.0) / 2.0);
                r0 = dist / 1.5;
                xr = r0 * f.cos();
                yr = r0 * f.sin();
                zr = zr / r * r0;
                break;
            }
        }
    } else {
        xr = 0.0;
        yr = 0.0;
        zr = 0.0;
    }

    let (ok, next) = if at[iat_h].valence != 0 {
        let next = nb(at, iat_h, 0);
        let pos = at[next].neighbor[..at[next].valence as usize]
            .iter()
            .position(|&x| x as usize == iat_h);
        if let Some(i) = pos {
            remove_inp_at_bond(at, next, i);
        }
        (pos.is_some(), next)
    } else {
        at[iat_h].valence = 1;
        at[iat_h].chem_bonds_valence = 1;
        at[iat_h].bond_type[0] = SINGLE;
        (true, iat_h)
    };
    if !ok {
        return 0;
    }
    if at[next].charge > 0 && at[iat].charge < 0 {
        at[next].charge -= 1;
        at[iat].charge += 1;
    }
    let val = at[iat].valence as usize;
    if val < MAXVAL {
        let bt = at[iat_h].bond_type[0];
        at[iat].neighbor[val] = iat_h as u16;
        at[iat].bond_type[val] = bt;
        at[iat].bond_stereo[val] = 0;
        at[iat].chem_bonds_valence += bt as i8;
        at[iat].valence = (val + 1) as i8;
    }
    at[iat_h].component = at[iat].component;
    at[iat_h].neighbor[0] = iat as u16;
    at[iat_h].bond_stereo[0] = 0;
    at[iat_h].x = at[iat].x + xr;
    at[iat_h].y = at[iat].y + yr;
    at[iat_h].z = at[iat].z + zr;
    1
}

/// `DisconnectMetals`: returns the number of disconnected bonds, or the
/// negated error code; on success `orig` gets the new atom array in which
/// the implicit H of the metals are explicit
pub fn disconnect_metals(
    orig: &mut OrigAtomData,
    b_check_metal_valence: bool,
    mut taut_flags_done: Option<&mut u64>,
) -> i32 {
    let _ = b_check_metal_valence;
    let num_atoms = orig.num_inp_atoms;
    let n_num_expl_h = if orig.b_disconnect_coord > 0 {
        (orig.b_disconnect_coord - 1) as usize
    } else {
        0
    };
    let mut at: Vec<InpAtom> = orig.at[..num_atoms].to_vec();
    at.resize(num_atoms + n_num_expl_h, InpAtom::default());
    let mut b_metal = vec![0u8; num_atoms + n_num_expl_h];
    let mut num_at = num_atoms;
    let mut num_disconnected = 0;
    let mut num_impl_h = 0usize;
    let mut err = 0;
    let mut num_explicit_h = [0i8; NUM_H_ISOTOPES + 1];

    'body: {
        for i in 0..num_atoms {
            let val = at[i].valence;
            if val == 0 && at[i].numh() == 0 {
                continue;
            }
            let rad_or_mult = val == 0
                || val != at[i].chem_bonds_valence
                || (at[i].radical != 0 && i32::from(at[i].radical) != RAD_SINGLET);
            if !rad_or_mult && is_ammonium_salt(&at, i, &mut num_explicit_h).is_some() {
            } else if !rad_or_mult && is_metal_salt(&at, i) {
            } else if is_metal_to_disconnect(&at, i) == 1 {
                let k = at[i].numh();
                num_impl_h += k as usize;
                b_metal[i] = (1 + k) as u8;
            }
            // (`bIsMetalToDisconnect` never returns 2, so
            // TG_FLAG_CHECK_VALENCE_COORD_DONE is never set)
        }
        if num_impl_h != n_num_expl_h {
            err = 2;
            break 'body;
        }

        // replace implicit H atoms with explicit H atoms
        let mut i = 0;
        while i < num_atoms && num_impl_h > 0 {
            if b_metal[i] <= 1 {
                i += 1;
                continue;
            }
            for k in 0..=NUM_H_ISOTOPES {
                let n = if k != 0 {
                    i32::from(at[i].num_iso_h[k - 1])
                } else {
                    i32::from(at[i].num_h)
                };
                for _ in 0..n {
                    if num_at >= num_atoms + n_num_expl_h {
                        err = 3;
                        break 'body;
                    }
                    at[num_at].elname = [0; ATOM_EL_LEN];
                    at[num_at].elname[0] = b'H';
                    at[num_at].el_number = EL_NUMBER_H;
                    at[num_at].iso_atw_diff = k as i8;
                    at[num_at].component = at[i].component;
                    move_explicit_h_cation(&mut at, num_at + 1, i, num_at, true);
                    at[num_at].orig_at_number = (num_at + 1) as u16;
                    num_at += 1;
                    num_impl_h -= 1;
                    b_metal[i] -= 1;
                    if k != 0 {
                        at[i].num_iso_h[k - 1] -= 1;
                    } else {
                        at[i].num_h -= 1;
                    }
                }
            }
            if b_metal[i] != 1 {
                err = 4;
                break 'body;
            }
            i += 1;
        }
        if num_at != num_atoms + n_num_expl_h {
            err = 5;
            break 'body;
        }

        // metal - ligand bonds, then metal - metal bonds
        for metal_pass in [false, true] {
            for i in 0..num_atoms {
                if b_metal[i] == 0 {
                    continue;
                }
                let mut j = i32::from(at[i].valence) - 1;
                while j >= 0 {
                    if j < i32::from(at[i].valence)
                        && (b_metal[nb(&at, i, j as usize)] != 0) == metal_pass
                    {
                        num_disconnected += disconnect_one_ligand(
                            &mut at,
                            orig.n_old_comp_number.as_deref_mut(),
                            &b_metal,
                            num_atoms,
                            i,
                            j as usize,
                            taut_flags_done.as_deref_mut(),
                        );
                    }
                    j -= 1;
                }
            }
        }
    }

    if num_disconnected == 0 {
        err = 6;
    }
    if err == 0 {
        orig.at = at;
        orig.num_inp_atoms = num_at;
    }
    if err != 0 { -err } else { num_disconnected }
}

// ---- terminal hydrogens ----

use super::ichisort::insertions_sort_at_rank;

/// `ATW_H`
const ATW_H: i32 = 1;

/// `add_DT_to_num_H`: adds the isotopic H counts to `num_H`
pub fn add_dt_to_num_h(at: &mut [InpAtom]) -> i32 {
    for a in at.iter_mut() {
        for j in 0..NUM_H_ISOTOPES {
            a.num_h += a.num_iso_h[j];
        }
    }
    0
}

/// `remove_terminal_HDT`: moves the terminal explicit H, D, T to the end of
/// the array and makes them implicit hydrogens of their neighbors. Returns
/// the new number of atoms (the first `ret` atoms of `at` are the heavy ones)
pub fn remove_terminal_hdt(at: &mut [InpAtom], b_fix_term_h_chrg: bool) -> i32 {
    const HDT: [u8; 4] = [b'H', b'D', b'T', 0];
    const K_MAX: i32 = 4;
    let num_atoms = at.len();
    let na = num_atoms as i32;
    let mut new_ord = vec![0u16; num_atoms];
    let mut new_at = vec![InpAtom::default(); num_atoms];
    let mut num_hydrogens: i32 = 0;
    let mut num_h = 0;

    for i in 0..num_atoms {
        at[i].component = i as u16;
        let mut k = if at[i].elname[1] != 0 {
            K_MAX
        } else {
            HDT.iter()
                .position(|&c| c == at[i].elname[0])
                .map_or(K_MAX, |p| p as i32)
        };
        if k == ATW_H || k == ATW_H + 1 {
            at[i].elname[0] = b'H';
            k += 1;
            at[i].iso_atw_diff = k as i8;
        }
        num_h += i32::from(
            k != K_MAX && at[i].valence == 1 && at[i].chem_bonds_valence == 1 && at[i].numh() == 0,
        );
    }

    if num_h == 2 && num_atoms == 2 && at[0].numh() == 0 && at[1].numh() == 0 {
        if at[0].iso_atw_diff >= at[1].iso_atw_diff {
            new_ord[0] = 0;
            new_ord[1] = 1;
        } else {
            new_ord[0] = 1;
            new_ord[1] = 0;
        }
        let (o0, o1) = (new_ord[0] as usize, new_ord[1] as usize);
        if at[o1].charge != 0 {
            at[o0].charge += at[o1].charge;
            at[o1].charge = 0;
        }
        new_at[o0] = at[0].clone();
        new_at[o1] = at[1].clone();
        num_hydrogens = 1;
    } else {
        for i in 0..num_atoms {
            let k = if at[i].elname[1] != 0 || at[i].numh() != 0 {
                K_MAX
            } else if at[i].elname[0] == b'H' {
                i32::from(at[i].iso_atw_diff)
            } else {
                K_MAX
            };
            let n = at[i].neighbor[0] as usize;
            if k < K_MAX
                && at[i].valence == 1
                && at[i].chem_bonds_valence == 1
                && (n > i || i32::from(new_ord[n]) < na - num_hydrogens)
            {
                num_hydrogens += 1;
                if at[i].charge != 0 {
                    at[n].charge += at[i].charge;
                    at[i].charge = 0;
                    if b_fix_term_h_chrg && i > n {
                        let t = new_ord[n] as usize;
                        new_at[t].charge = at[n].charge;
                    }
                }
                new_ord[i] = (na - num_hydrogens) as u16;
            } else {
                new_ord[i] = (i as i32 - num_hydrogens) as u16;
            }
            new_at[new_ord[i] as usize] = at[i].clone();
        }
    }

    if num_hydrogens == 0 {
        return na;
    }
    let num_others = (na - num_hydrogens) as usize;
    if num_hydrogens > 1 {
        new_at[num_others..].sort_by(|a, b| cmp_iso_atw_diff_component_no(a, b).cmp(&0));
    }
    for i in num_others..num_atoms {
        new_ord[new_at[i].component as usize] = i as u16;
    }

    for i in 0..num_others {
        let mut new_hydrogen_at_order = [0u16; NUM_H_ISOTOPES + 1];
        let mut new_other_neigh_order = [0u16; MAXVAL];
        let mut num_other_neigh = 0usize;
        let mut num_hydrogen_at = 0usize;
        let mut num_h_flag = 0;
        let mut val = 0usize;

        for k in 0..new_at[i].valence as usize {
            let n = new_ord[new_at[i].neighbor[k] as usize] as usize;
            if num_others <= n {
                let iso = i32::from(new_at[n].iso_atw_diff);
                if 0 < iso && iso < ATW_H + NUM_H_ISOTOPES as i32 {
                    new_at[i].num_iso_h[iso as usize - 1] += 1;
                    num_hydrogen_at += usize::from(new_hydrogen_at_order[iso as usize] == 0);
                    new_hydrogen_at_order[iso as usize] = (k + 1) as u16;
                } else {
                    new_at[i].num_h += 1;
                    num_hydrogen_at += usize::from(num_h_flag == 0);
                    num_h_flag += 1;
                    new_hydrogen_at_order[0] = (k + 1) as u16;
                }
                new_at[i].chem_bonds_valence = (new_at[i].chem_bonds_valence - 1).max(0);
                new_at[n].neighbor[0] = i as u16;
                if new_at[i].sb_parity[0] != 0 {
                    let mut m = 0;
                    while m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m] != 0 {
                        if k as i32 == i32::from(new_at[i].sn_ord[m]) {
                            new_at[i].sn_ord[m] = -(new_at[n].iso_atw_diff + 1);
                        }
                        m += 1;
                    }
                }
            } else {
                if new_at[i].sb_parity[0] != 0 {
                    if num_other_neigh < MAX_NUM_STEREO_BONDS {
                        new_other_neigh_order[num_other_neigh] = (k + 1) as u16;
                    }
                    num_other_neigh += 1;
                    if val != k {
                        let mut m = 0;
                        while m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m] != 0 {
                            if k as i32 == i32::from(new_at[i].sb_ord[m]) {
                                new_at[i].sb_ord[m] = val as i8;
                            } else if k as i32 == i32::from(new_at[i].sn_ord[m]) {
                                new_at[i].sn_ord[m] = val as i8;
                            }
                            m += 1;
                        }
                    }
                }
                new_at[i].neighbor[val] = new_ord[new_at[i].neighbor[k] as usize];
                new_at[i].bond_type[val] = new_at[i].bond_type[k];
                new_at[i].bond_stereo[val] = new_at[i].bond_stereo[k];
                val += 1;
            }
        }
        if new_at[i].valence as usize > val
            && new_at[i].sb_parity[0] != 0
            && num_hydrogen_at == new_at[i].valence as usize - val
            && num_hydrogen_at + num_other_neigh <= MAXVAL
        {
            new_other_neigh_order.copy_within(0..num_other_neigh, num_hydrogen_at);
            let mut j = 1usize;
            for k in 0..=NUM_H_ISOTOPES {
                if new_hydrogen_at_order[k] != 0 {
                    new_other_neigh_order[num_hydrogen_at - j] = new_hydrogen_at_order[k];
                    let mut m = 0;
                    while m < MAX_NUM_STEREO_BONDS && new_at[i].sb_parity[m] != 0 {
                        if i32::from(new_at[i].sn_ord[m]) == -(k as i32 + 1) {
                            new_at[i].sn_ord[m] = -(j as i8);
                        }
                        m += 1;
                    }
                    j += 1;
                }
            }
            insertions_sort_at_rank(
                &mut new_other_neigh_order[..num_hydrogen_at + num_other_neigh],
            );
        }
        new_at[i].valence = val as i8;
    }
    at.clone_from_slice(&new_at);
    num_others as i32
}

// ---- connected components ----

/// `get_iat_number`: the position of an element in the `ION_ATOM_TYPE` list
/// (H, C, N, P, O, S, Se, Te, F, Cl, Br, I), -1 for others
pub fn get_iat_number(el_number: u8) -> i32 {
    use super::util::{EL_NUMBER_SE, EL_NUMBER_TE};
    match el_number {
        EL_NUMBER_H => 0,
        EL_NUMBER_C => 1,
        EL_NUMBER_N => 2,
        super::util::EL_NUMBER_P => 3,
        EL_NUMBER_O => 4,
        super::util::EL_NUMBER_S => 5,
        EL_NUMBER_SE => 6,
        EL_NUMBER_TE => 7,
        EL_NUMBER_F => 8,
        EL_NUMBER_CL => 9,
        EL_NUMBER_BR => 10,
        EL_NUMBER_I => 11,
        _ => -1,
    }
}

/// The target valence of an atom whose hydrogens may be exchanged
/// (`None`: the atom is excluded); `is_h` is set for an isolated proton
fn xchg_valence(a: &InpAtom, is_h: &mut bool) -> Option<i32> {
    let iat = get_iat_number(a.el_number);
    if iat < 0 {
        return None;
    }
    if a.charge.abs() > 1 || (a.radical != 0 && i32::from(a.radical) != RAD_SINGLET) {
        return None;
    }
    let chg = i32::from(a.charge);
    let mut val = -1;
    match iat {
        2 | 3 => val = 3 + chg,
        4..=7 => val = 2 + chg,
        8..=11 => {
            if chg == 0 {
                val = 1;
            }
        }
        0 if a.valence == 0 && chg == 1 => {
            *is_h = true;
            val = 0;
        }
        _ => {}
    }
    (val >= 0).then_some(val)
}

/// `bHeteroAtomMayHaveXchgIsoH`: 0 = no, 1 = heteroatom, 2 = isolated proton
pub fn hetero_atom_may_have_xchg_iso_h(atom: &[InpAtom], iat: usize) -> i32 {
    let at = &atom[iat];
    let mut is_h = false;
    let Some(val) = xchg_valence(at, &mut is_h) else {
        return 0;
    };
    if val != i32::from(at.chem_bonds_valence) + at.numh() {
        return 0;
    }
    if is_h {
        return 2;
    }
    for j in 0..at.valence as usize {
        let at2 = &atom[at.neighbor[j] as usize];
        if (at2.charge != 0 && at.charge != 0)
            || (at2.radical != 0 && i32::from(at2.radical) != RAD_SINGLET)
        {
            return 0;
        }
    }
    1
}

/// `bNumHeterAtomHasIsotopicH`: bit 0 = a neutralizable heteroatom has an
/// isotopic H, bit 1 = the structure has isotopic atoms
pub fn num_hetero_atom_has_isotopic_h(atom: &[InpAtom]) -> i32 {
    let mut is_h = false;
    let mut num_iso_h = 0;
    let mut num_iso_atoms = 0i32;
    for at in atom {
        num_iso_atoms += i32::from(at.iso_atw_diff != 0 || at.num_iso_h_total() != 0);
        let Some(val) = xchg_valence(at, &mut is_h) else {
            continue;
        };
        if val != i32::from(at.chem_bonds_valence) + at.numh() {
            continue;
        }
        let mut b_accept = true;
        let mut cur_num_iso_h;
        if is_h {
            cur_num_iso_h = i32::from(at.iso_atw_diff != 0);
        } else {
            cur_num_iso_h = 0;
            for j in 0..at.valence as usize {
                let at2 = &atom[at.neighbor[j] as usize];
                if (at2.charge != 0 && at.charge != 0)
                    || (at2.radical != 0 && i32::from(at2.radical) != RAD_SINGLET)
                {
                    b_accept = false;
                    break;
                } else if at2.el_number == EL_NUMBER_H && at2.valence == 1 && at2.iso_atw_diff != 0
                {
                    cur_num_iso_h += 1;
                }
            }
            if b_accept {
                num_iso_atoms -= cur_num_iso_h;
                cur_num_iso_h += at.num_iso_h_total();
            }
        }
        num_iso_h += i32::from(b_accept && cur_num_iso_h != 0);
    }
    i32::from(num_iso_h != 0) | if num_iso_atoms != 0 { 2 } else { 0 }
}

/// `MarkDisconnectedComponents`: sets `at[*].component` (1-based, components
/// sorted by decreasing number of atoms, stable); returns the number of
/// components
pub fn mark_disconnected_components(
    orig: &mut OrigAtomData,
    mut b_process_old_comp_numbers: bool,
) -> i32 {
    if b_process_old_comp_numbers && orig.n_old_comp_number.is_none() {
        b_process_old_comp_numbers = false;
    }
    let num_at = orig.num_inp_atoms;
    if num_at == 0 {
        return 0;
    }
    let at = &mut orig.at;
    let mut n_new_comp = vec![0u16; num_at];
    let mut n_prev = vec![0u16; num_at];
    let mut i_neigh = vec![0i8; num_at];
    let mut num_components = 0usize;
    for j in 0..num_at {
        if n_new_comp[j] != 0 {
            continue;
        }
        num_components += 1;
        let fst = j;
        let mut cur = j;
        n_new_comp[fst] = num_components as u16;
        loop {
            if i_neigh[cur] < at[cur].valence {
                let nxt = at[cur].neighbor[i_neigh[cur] as usize] as usize;
                i_neigh[cur] += 1;
                if n_new_comp[nxt] == 0 {
                    n_new_comp[nxt] = num_components as u16;
                    n_prev[nxt] = cur as u16;
                    cur = nxt;
                }
            } else if cur == fst {
                break;
            } else {
                cur = n_prev[cur] as usize;
            }
        }
    }

    let size = num_components.max(orig.num_components.max(0) as usize);
    let mut n_old = vec![0u16; size + 1];
    // [0] atoms, [1] original number, [2] new number
    let mut comp: Vec<[u16; 3]> = (0..num_components).map(|i| [0, i as u16, 0]).collect();
    for &c in &n_new_comp {
        comp[c as usize - 1][0] += 1;
    }
    comp.sort_by(|a, b| b[0].cmp(&a[0]).then(a[1].cmp(&b[1])));
    let mut n_cur_len = vec![0u16; num_components + 1];
    for i in 0..num_components {
        n_cur_len[i] = comp[i][0];
        let o = comp[i][1] as usize;
        comp[o][2] = (i + 1) as u16;
    }
    let no_component = (num_at + 1) as u16;
    for j in 0..num_at {
        let new_comp_no = comp[n_new_comp[j] as usize - 1][2] as usize - 1;
        if b_process_old_comp_numbers {
            let old_comp_no = at[j].component;
            if old_comp_no == 0 {
                n_old[new_comp_no] = no_component;
            } else if n_old[new_comp_no] != old_comp_no {
                if n_old[new_comp_no] == 0 {
                    n_old[new_comp_no] = old_comp_no;
                } else {
                    let another = n_old[new_comp_no];
                    for v in n_old.iter_mut().take(num_components) {
                        if *v == old_comp_no || *v == another {
                            *v = no_component;
                        }
                    }
                }
            }
        }
        at[j].component = (new_comp_no + 1) as u16;
    }
    if b_process_old_comp_numbers {
        let prev = orig.n_old_comp_number.as_ref().unwrap();
        for v in n_old.iter_mut().take(num_components) {
            if *v == no_component {
                *v = 0;
            } else if *v != 0 && prev[*v as usize - 1] == 0 {
                *v = 0;
            }
        }
    } else {
        for (j, v) in n_old.iter_mut().take(num_components).enumerate() {
            *v = (j + 1) as u16;
        }
    }
    orig.n_cur_at_len = Some(n_cur_len);
    orig.n_old_comp_number = Some(n_old);
    orig.num_components = num_components as i32;
    num_components as i32
}

/// `ExtractConnectedComponent`: copies the atoms of a component to
/// `component_at` with renumbered neighbors; returns their number
pub fn extract_connected_component(
    at: &[InpAtom],
    component_number: u16,
    component_at: &mut Vec<InpAtom>,
) -> i32 {
    let mut number = vec![0u16; at.len()];
    component_at.clear();
    for (i, a) in at.iter().enumerate() {
        if a.component == component_number {
            number[i] = component_at.len() as u16;
            component_at.push(a.clone());
        }
    }
    for (i, a) in component_at.iter_mut().enumerate() {
        a.orig_compt_at_numb = (i + 1) as u16;
        for j in 0..a.valence as usize {
            a.neighbor[j] = number[a.neighbor[j] as usize];
        }
    }
    component_at.len() as i32
}

/// `SetConnectedComponentNumber`
pub fn set_connected_component_number(at: &mut [InpAtom], component_number: u16) -> i32 {
    for a in at.iter_mut() {
        a.component = component_number;
    }
    0
}

/// `post_fix_odd_things`: currently does nothing in the C code
pub fn post_fix_odd_things(_at: &mut [InpAtom]) -> i32 {
    0
}

// ---- allocation of the `INChI` / `INChI_Aux` result structures ----

/// `Alloc_INChI_Stereo`
pub fn alloc_inchi_stereo(num_at: usize, num_bonds: usize) -> super::ichi::InChIStereo {
    super::ichi::InChIStereo {
        n_number: vec![0; num_at],
        t_parity: vec![0; num_at],
        n_number_inv: vec![0; num_at],
        t_parity_inv: vec![0; num_at],
        n_bond_atom1: vec![0; num_bonds],
        n_bond_atom2: vec![0; num_bonds],
        b_parity: vec![0; num_bonds],
        ..Default::default()
    }
}

/// `Alloc_INChI`: returns the zeroed structure with the numbers of bonds and
/// isotopic atoms found; `None` for no atoms.
pub fn alloc_inchi(
    at: &[InpAtom],
    num_at: usize,
    n_alloc_mode: u64,
) -> Option<(super::ichi::InChI, usize, usize)> {
    use super::ichi::{InChI, InChIIsotopicAtom, InChIIsotopicTGroup};
    use super::ichiparm::REQ_MODE_ISO;
    const INCHI_T_NUM_MOVABLE: usize = 2;
    if num_at == 0 {
        return None;
    }
    let mut num_bonds = 0usize;
    let mut num_isotopic = 0usize;
    for a in &at[..num_at] {
        num_bonds += a.valence.max(0) as usize;
        let name_end = a
            .elname
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(a.elname.len());
        let name = &a.elname[..name_end];
        num_isotopic += usize::from(
            a.iso_atw_diff != 0
                || name == b"D"
                || name == b"T"
                || a.num_iso_h.iter().any(|&x| x != 0),
        );
    }
    num_bonds /= 2;
    let mut p = InChI {
        n_atom: vec![0; num_at],
        n_conn_table: vec![0; num_at + num_bonds],
        n_tautomer: vec![0; (3 + INCHI_T_NUM_MOVABLE) * num_at / 2 + 1],
        n_num_h: vec![0; num_at],
        n_num_h_fixed: Some(vec![0; num_at]),
        stereo: Some(alloc_inchi_stereo(num_at, num_bonds)),
        ..Default::default()
    };
    if n_alloc_mode & REQ_MODE_ISO != 0 {
        if num_isotopic > 0 {
            p.isotopic_atom = vec![InChIIsotopicAtom::default(); num_isotopic];
            p.isotopic_t_group = vec![InChIIsotopicTGroup::default(); num_isotopic];
        }
        p.n_possible_locations_of_isotopic_h = vec![0; num_at + 1];
        p.stereo_isotopic = Some(alloc_inchi_stereo(num_at, num_bonds));
    }
    Some((p, num_bonds, num_isotopic))
}

/// `Alloc_INChI_Aux`
pub fn alloc_inchi_aux(
    num_at: usize,
    _num_isotopic_atoms: usize,
    n_alloc_mode: u64,
) -> Option<super::ichi::InChIAux> {
    use super::ichi::{InChIAux, OrigInfo};
    use super::ichiparm::REQ_MODE_ISO;
    if num_at == 0 {
        return None;
    }
    let num_at_tg = num_at + num_at / 2;
    let mut a = InChIAux {
        n_orig_at_nos_in_canon_ord: vec![0; num_at_tg],
        n_orig_at_nos_in_canon_ord_inv: vec![0; num_at_tg],
        n_constit_equ_numbers: vec![0; num_at_tg],
        orig_info: vec![OrigInfo::default(); num_at],
        ..Default::default()
    };
    if num_at > 1 {
        a.n_constit_equ_t_group_numbers = vec![0; num_at / 2 + 1];
    }
    if n_alloc_mode & REQ_MODE_ISO != 0 {
        a.n_isotopic_orig_at_nos_in_canon_ord = vec![0; num_at_tg];
        a.n_isotopic_orig_at_nos_in_canon_ord_inv = vec![0; num_at_tg];
        a.n_constit_equ_isotopic_numbers = vec![0; num_at_tg];
        a.n_constit_equ_isotopic_t_group_numbers = vec![0; num_at / 2 + 1];
    }
    Some(a)
}
