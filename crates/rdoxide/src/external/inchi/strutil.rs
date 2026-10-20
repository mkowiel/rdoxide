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
