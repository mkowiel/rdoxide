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
