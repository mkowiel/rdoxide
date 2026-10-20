//! Chemical-element helpers of InChI (`INCHI_BASE/src/util.c`, the parts that
//! do not depend on the atom structures): element lookup, standard valences,
//! the number of implicit hydrogens, and the parsing of charges, radicals and
//! hydrogens in element names.
//!
//! The element table is generated, see [`super::el_data`].
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::el_data::EL_DATA;
use super::incomdef::{
    IS_METAL, NUM_H_ISOTOPES, RADICAL_DOUBLET, RADICAL_SINGLET, RADICAL_TRIPLET,
};

pub const EL_NUMBER_H: u8 = 1;
pub const EL_NUMBER_B: u8 = 5;
pub const EL_NUMBER_C: u8 = 6;
pub const EL_NUMBER_N: u8 = 7;
pub const EL_NUMBER_O: u8 = 8;
pub const EL_NUMBER_F: u8 = 9;
pub const EL_NUMBER_SI: u8 = 14;
pub const EL_NUMBER_P: u8 = 15;
pub const EL_NUMBER_S: u8 = 16;
pub const EL_NUMBER_CL: u8 = 17;
pub const EL_NUMBER_GE: u8 = 32;
pub const EL_NUMBER_AS: u8 = 33;
pub const EL_NUMBER_SE: u8 = 34;
pub const EL_NUMBER_BR: u8 = 35;
pub const EL_NUMBER_SB: u8 = 51;
pub const EL_NUMBER_TE: u8 = 52;
pub const EL_NUMBER_I: u8 = 53;
pub const EL_NUMBER_PO: u8 = 84;
pub const EL_NUMBER_AT: u8 = 85;

pub const MIN_ATOM_CHARGE: i32 = -2;
pub const MAX_ATOM_CHARGE: i32 = 2;
/// index of the neutral state in the valence table
pub const NEUTRAL_STATE: i32 = -MIN_ATOM_CHARGE;
/// max. number of valences + 1 (zero termination)
pub const MAX_NUM_VALENCES: usize = 5;

/// `ERR_ELEM`
pub const ERR_ELEM: i32 = 255;

/// `nElDataLen`: the table without its terminating empty entry
const N_EL_DATA_LEN: i32 = EL_DATA.len() as i32 - 1;

/// The index of the element `n_at_num` in the table (which has D and T after H).
fn table_index(n_at_num: i32) -> i32 {
    let mut n = n_at_num - 1;
    if n > 0 {
        // bypass D, T
        n += 2;
    }
    n
}

/// `get_element_chemical_symbol`: the symbol of the element, or `"??"` and
/// `-1` for an unknown number.
pub fn get_element_chemical_symbol(n_at_num: i32) -> (i32, &'static str) {
    let n = table_index(n_at_num);
    if (0..N_EL_DATA_LEN).contains(&n) {
        (0, EL_DATA[n as usize].name)
    } else {
        (-1, "??")
    }
}

/// `get_element_or_pseudoelement_symbol`: like the above, translating the
/// pseudoelement `Zy` to `Zz`.
pub fn get_element_or_pseudoelement_symbol(n_at_num: i32) -> (i32, &'static str) {
    match get_element_chemical_symbol(n_at_num) {
        (0, "Zy") => (0, "Zz"),
        r => r,
    }
}

/// `el_number_in_internal_ref_table`: the index in the table (with D and T),
/// or [`ERR_ELEM`].
pub fn el_number_in_internal_ref_table(elname: &str) -> i32 {
    EL_DATA
        .iter()
        .take_while(|e| !e.name.is_empty())
        .position(|e| e.name == elname)
        .map_or(ERR_ELEM, |i| i as i32)
}

/// `get_periodic_table_number`: the element number by symbol.
pub fn get_periodic_table_number(elname: &str) -> i32 {
    let b = elname.as_bytes();
    if b.len() == 1 {
        match b[0] {
            b'H' => return i32::from(EL_NUMBER_H),
            b'B' => return i32::from(EL_NUMBER_B),
            b'C' => return i32::from(EL_NUMBER_C),
            b'N' => return i32::from(EL_NUMBER_N),
            b'O' => return i32::from(EL_NUMBER_O),
            b'P' => return i32::from(EL_NUMBER_P),
            b'S' => return i32::from(EL_NUMBER_S),
            b'F' => return i32::from(EL_NUMBER_F),
            b'I' => return i32::from(EL_NUMBER_I),
            _ => {}
        }
    }
    let num = el_number_in_internal_ref_table(elname);
    if num < ERR_ELEM {
        // D and T are in the internal table, but not in Mendeleev's
        return 1.max(num - 1);
    }
    num
}

/// The table entry of the periodic number (`ElData[n > 1 ? n + 1 : 0]`).
fn el_by_periodic_num(n_periodic_num: i32) -> &'static super::el_data::ElData {
    &EL_DATA[if n_periodic_num > 1 {
        (n_periodic_num + 1) as usize
    } else {
        0
    }]
}

/// `if_skip_add_H`: whether no hydrogens are added to the element.
pub fn if_skip_add_h(n_periodic_num: i32) -> bool {
    el_by_periodic_num(n_periodic_num).skip_adding_h
}

/// `get_el_valence`: the reference valence `val_num` of the element at the charge.
pub fn get_el_valence(n_periodic_num: i32, charge: i32, val_num: usize) -> i32 {
    if !(MIN_ATOM_CHARGE..=MAX_ATOM_CHARGE).contains(&charge) || val_num >= MAX_NUM_VALENCES {
        return 0;
    }
    i32::from(
        el_by_periodic_num(n_periodic_num).valence[(NEUTRAL_STATE + charge) as usize][val_num],
    )
}

fn radical_adjustment(radical: i32, singlet_counts: bool) -> i32 {
    if radical == RADICAL_DOUBLET {
        1
    } else if radical == RADICAL_TRIPLET || (singlet_counts && radical == RADICAL_SINGLET) {
        2
    } else {
        0
    }
}

/// `get_unusual_el_valence`: the valence needed to unambiguously reconstruct the bonds.
pub fn get_unusual_el_valence(
    n_periodic_num: i32,
    charge: i32,
    radical: i32,
    bonds_valence: i32,
    num_h: i32,
    num_bonds: i32,
) -> i32 {
    if num_bonds == 0 && num_h == 0 {
        return 0;
    }
    if !(MIN_ATOM_CHARGE..=MAX_ATOM_CHARGE).contains(&charge) {
        // all single bonds
        return if bonds_valence == num_bonds {
            0
        } else {
            bonds_valence
        };
    }
    if get_el_valence(n_periodic_num, charge, 0) == 0 && bonds_valence == num_bonds {
        return 0;
    }
    let chem_valence = bonds_valence + num_h;
    let rad_adj = radical_adjustment(radical, false);
    let mut num_found = 0;
    let mut exact_found = false;
    for i in 0..MAX_NUM_VALENCES {
        let known = get_el_valence(n_periodic_num, charge, i) - rad_adj;
        if 0 < known && num_bonds <= known && known <= chem_valence {
            num_found += 1;
            if known == chem_valence {
                exact_found = true;
                break;
            }
        }
    }
    if exact_found && num_found == 1 {
        0
    } else {
        chem_valence
    }
}

/// `needed_unusual_el_valence`: the valence needed to unambiguously
/// reconstruct the number of hydrogens.
pub fn needed_unusual_el_valence(
    n_periodic_num: i32,
    charge: i32,
    radical: i32,
    bonds_valence: i32,
    actual_bonds_valence: i32,
    num_h: i32,
    num_bonds: i32,
) -> i32 {
    let (ret, symbol) = get_element_chemical_symbol(n_periodic_num);
    let num_h_expected = if num_bonds != 0 && ret != -1 {
        get_num_h(
            symbol,
            0,
            None,
            charge,
            radical,
            actual_bonds_valence,
            0,
            false,
            false,
            false,
        )
    } else {
        num_h
    };
    let chem_valence = bonds_valence + num_h;
    // `BUILD_WITH_ENG_OPTIONS` is not set in the library build
    if !(MIN_ATOM_CHARGE..=MAX_ATOM_CHARGE).contains(&charge)
        || get_el_valence(n_periodic_num, charge, 0) == 0
        || if_skip_add_h(n_periodic_num)
        || bonds_valence != actual_bonds_valence
        || num_h_expected != num_h
    {
        if num_h == 0 && num_h_expected == 0 && bonds_valence == actual_bonds_valence {
            // no H
            return 0;
        }
        // needs to add H-atoms
        return chem_valence;
    }
    let rad_adj = radical_adjustment(radical, false);
    let (mut num_found_known, mut num_found, mut exact_found) = (0, 0, false);
    for i in 0..MAX_NUM_VALENCES {
        let mut known = get_el_valence(n_periodic_num, charge, i);
        if 0 < known {
            known -= rad_adj;
            if bonds_valence <= known {
                // known valence that fits without H
                num_found_known += 1;
                if known <= chem_valence {
                    // large enough to accommodate the (implicit) H
                    num_found += 1;
                }
                if known == chem_valence {
                    exact_found = true;
                    break;
                }
            }
        }
    }
    if exact_found && num_found == 1 && num_found_known == 1 {
        0
    } else if chem_valence != 0 {
        chem_valence
    } else {
        // needs zero
        -1
    }
}

/// `detect_unusual_el_valence`: the valence that does not fit any known valence.
pub fn detect_unusual_el_valence(
    n_periodic_num: i32,
    charge: i32,
    radical: i32,
    bonds_valence: i32,
    num_h: i32,
    num_bonds: i32,
) -> i32 {
    if num_bonds == 0 && num_h == 0 {
        return 0;
    }
    if !(MIN_ATOM_CHARGE..=MAX_ATOM_CHARGE).contains(&charge) {
        return if bonds_valence == num_bonds {
            0
        } else {
            bonds_valence
        };
    }
    if get_el_valence(n_periodic_num, charge, 0) == 0 && bonds_valence == num_bonds {
        return 0;
    }
    let chem_valence = bonds_valence + num_h;
    let rad_adj = radical_adjustment(radical, true);
    for i in 0..MAX_NUM_VALENCES {
        let known = get_el_valence(n_periodic_num, charge, i) - rad_adj;
        if 0 < known && known == chem_valence {
            return 0;
        }
    }
    chem_valence
}

/// `get_el_type`: the type of the element after `n_periodic_num`'s entry in
/// the table (as in the C code, which indexes `ElData[n + 1]`).
pub fn get_el_type(n_periodic_num: i32) -> i32 {
    EL_DATA[(n_periodic_num + 1) as usize].el_type
}

/// `is_el_a_metal`
pub fn is_el_a_metal(n_periodic_num: i32) -> bool {
    EL_DATA[(n_periodic_num + 1) as usize].el_type & IS_METAL != 0
}

/// C's `strtol(s, &end, 10)`: the value and the number of bytes consumed
/// (0 when there is no number).
fn strtol10(s: &[u8]) -> (i64, usize) {
    let mut i = 0;
    while s
        .get(i)
        .is_some_and(|c| c.is_ascii_whitespace() || *c == 0x0b)
    {
        i += 1;
    }
    let neg = match s.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let start = i;
    let mut v: i64 = 0;
    while let Some(c) = s.get(i).filter(|c| c.is_ascii_digit()) {
        v = v.saturating_mul(10).saturating_add(i64::from(c - b'0'));
        i += 1;
    }
    if i == start {
        return (0, 0);
    }
    (if neg { -v } else { v }, i)
}

/// `extract_charges_and_radicals`: removes the charges (`+`, `-`, `+2`, `--`),
/// radical marks (`^`, `:`, trailing `.`) from `elname` and returns
/// `(radical, charge)`.
pub fn extract_charges_and_radicals(elname: &mut Vec<u8>) -> (i32, i32) {
    let mut n_charge: i32 = 0;
    let mut n_rad: i32 = 0;
    let mut n_last_sign: i32 = 1;
    while let Some(q) = elname.iter().position(|c| matches!(c, b'+' | b'-' | b'^')) {
        let charge_len;
        match elname[q] {
            b'+' | b'-' => {
                let mut n_val: i32 = 0;
                let mut k = 0;
                let mut len = 0;
                loop {
                    let sign = match elname.get(q + k) {
                        Some(b'+') => 1,
                        Some(b'-') => -1,
                        _ => break,
                    };
                    n_last_sign = sign;
                    n_val += sign;
                    len += 1;
                    k += 1;
                }
                let (v, used) = strtol10(&elname[q + k..]);
                let n_sign = v as i32;
                if n_sign != 0 {
                    n_val += n_last_sign * (n_sign - 1);
                }
                // `charge_len = r - q`: after the digits, or just after the signs
                charge_len = k + used;
                let _ = len;
                n_charge += n_val;
            }
            _ => {
                // '^': a doublet is 1 here
                n_rad = 1;
                let mut k = 1;
                while elname.get(q + k) == Some(&b'^') {
                    n_rad += 1;
                    k += 1;
                }
                charge_len = k;
            }
        }
        elname.drain(q..q + charge_len);
    }
    // radical
    if elname.last() == Some(&b':') {
        n_rad = RADICAL_SINGLET;
        elname.pop();
    } else {
        while elname.last() == Some(&b'.') {
            n_rad += 1;
            elname.pop();
        }
        n_rad = match n_rad {
            1 => RADICAL_DOUBLET,
            2 => RADICAL_TRIPLET,
            _ => 0,
        };
    }
    (n_rad, n_charge)
}

/// `extract_H_atoms`: removes the hydrogen atoms (`H`, `D`, `T` with an
/// optional count) from `elname`; returns the number of protium atoms and
/// adds the numbers of D and T to `num_iso_h[1]` and `num_iso_h[2]`.
pub fn extract_h_atoms(elname: &mut Vec<u8>, num_iso_h: &mut [i8; NUM_H_ISOTOPES]) -> i32 {
    let at = |v: &Vec<u8>, i: usize| v.get(i).copied().unwrap_or(0);
    let mut i = 0usize;
    let mut num_h: i32 = 0;
    let mut len = elname.len();
    let elname1 = if len > 1 { elname[1] } else { 0 };
    let mut c = at(elname, 0);
    while i < len {
        let k: i32 = match c {
            b'H' => 0,
            b'D' => 1,
            b'T' => 2,
            _ => -1,
        };
        // the character after elname[i]
        let mut q = i + 1;
        c = at(elname, q);
        if k >= 0 && !c.is_ascii_lowercase() {
            // a hydrogen atom
            let val = if c.is_ascii_digit() {
                let (v, used) = strtol10(&elname[q..]);
                q += used;
                v as i32
            } else {
                1
            };
            if k != 0 {
                let slot = &mut num_iso_h[k as usize];
                *slot = slot.wrapping_add(val as i8);
            } else {
                num_h += val;
            }
            // remove the hydrogen atom from the string
            len -= q - i;
            elname.drain(i..q);
        } else {
            i += 1;
        }
        c = at(elname, i);
    }
    let len = elname.len();
    if len == 2 && elname[1] != elname1 {
        // the 2nd character changed after removing {H,D,T}: an error
        elname[1] = b'?';
    }
    num_h
}

/// `get_num_H`: the number of attached (implicit) hydrogens.
#[allow(clippy::too_many_arguments)]
pub fn get_num_h(
    elname: &str,
    inp_num_h: i32,
    inp_num_iso_h: Option<&[i8; NUM_H_ISOTOPES]>,
    charge: i32,
    radical: i32,
    chem_bonds_valence: i32,
    atom_input_valence: i32,
    b_aliased: bool,
    b_do_not_add_h: bool,
    b_has_metal_neighbor: bool,
) -> i32 {
    let el_number = el_number_in_internal_ref_table(elname);
    let mut num_h;
    if b_aliased {
        num_h = inp_num_h;
    } else if atom_input_valence != 0 && (atom_input_valence != 15 || chem_bonds_valence != 0) {
        num_h = 0.max(atom_input_valence - chem_bonds_valence);
    } else if atom_input_valence == 15 && chem_bonds_valence == 0 {
        num_h = 0;
    } else if (MIN_ATOM_CHARGE..=MAX_ATOM_CHARGE).contains(&charge)
        && el_number != ERR_ELEM
        && !EL_DATA[el_number as usize].skip_adding_h
        && !b_do_not_add_h
    {
        let valences = &EL_DATA[el_number as usize].valence[(NEUTRAL_STATE + charge) as usize];
        num_h = 0;
        if radical != 0 && radical != RADICAL_SINGLET {
            let mut val = i32::from(valences[0]);
            if val != 0 {
                val -= if radical == RADICAL_DOUBLET {
                    1
                } else if radical == RADICAL_SINGLET || radical == RADICAL_TRIPLET {
                    2
                } else {
                    val
                };
                // an unknown radical adds no H
                num_h = 0.max(val - chem_bonds_valence);
            }
        } else {
            // the smallest valence that is greater than the sum of the bond valences
            let mut i = 0;
            let mut val;
            loop {
                val = valences.get(i).map_or(0, |v| i32::from(*v));
                if val != 0 && val < chem_bonds_valence {
                    i += 1;
                } else {
                    break;
                }
            }
            let n = el_number_in_internal_ref_table("N");
            let s = el_number_in_internal_ref_table("S");
            let c = el_number_in_internal_ref_table("C");
            if el_number == n && charge == 0 && radical == 0 && val == 5 {
                // do not add H to N(IV)
                val = 3;
            } else if el_number == s
                && charge == 0
                && radical == 0
                && val == 4
                && chem_bonds_valence == 3
            {
                // nor to S(III)
                val = 3;
            } else if b_has_metal_neighbor && el_number != c && val > 0 {
                val -= 1;
            }
            num_h = 0.max(val - chem_bonds_valence);
        }
        let num_iso_h: i32 = inp_num_iso_h.map_or(0, |a| a.iter().map(|v| i32::from(*v)).sum());
        if num_iso_h != 0 {
            if num_h >= num_iso_h {
                num_h -= num_iso_h;
            } else {
                // as requested in the alias
                num_h = inp_num_h;
            }
        }
        if inp_num_h > num_h {
            num_h = inp_num_h;
        }
    } else {
        num_h = inp_num_h;
    }
    num_h
}

/// `get_atomic_mass_from_elnum`
pub fn get_atomic_mass_from_elnum(n_at_num: i32) -> i32 {
    let n = table_index(n_at_num);
    if (0..N_EL_DATA_LEN).contains(&n) {
        EL_DATA[n as usize].at_mass
    } else {
        0
    }
}

/// `get_atomic_mass`
pub fn get_atomic_mass(elname: &str) -> i32 {
    let n = el_number_in_internal_ref_table(elname);
    if n != ERR_ELEM {
        EL_DATA[n as usize].at_mass
    } else {
        0
    }
}
