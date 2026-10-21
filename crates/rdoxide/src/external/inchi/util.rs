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
    BOND_TYPE_ALTERN, BOND_TYPE_TRIPLE, IS_METAL, NUM_H_ISOTOPES, RADICAL_DOUBLET, RADICAL_SINGLET,
    RADICAL_TRIPLET,
};
use super::inpdef::{BOND_TYPE_MASK, InpAtom};

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
pub(crate) fn strtol10(s: &[u8]) -> (i64, usize) {
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

// ---- functions on `inp_ATOM` ----

/// `S_VI_O_PLUS_METAL_FIX_BOND == 1`
const S_VI_O_PLUS_METAL_FIX_BOND: bool = true;

/// `get_endpoint_valence`
pub fn get_endpoint_valence(el_number: u8) -> i32 {
    match el_number {
        EL_NUMBER_O | EL_NUMBER_S | EL_NUMBER_SE | EL_NUMBER_TE => 2,
        EL_NUMBER_N => 3,
        _ => 0,
    }
}

/// `is_in_the_ilist`: the position of `next_atom` in the list
pub fn is_in_the_ilist(path_atom: &[i32], next_atom: i32) -> Option<usize> {
    path_atom.iter().position(|&a| a == next_atom)
}

/// `is_ilist_inside`: whether every element of `ilist` is in `ilist2`
pub fn is_ilist_inside(ilist: &[i32], ilist2: &[i32]) -> bool {
    ilist.iter().all(|&a| is_in_the_ilist(ilist2, a).is_some())
}

/// `nBondsValToMetal`: the sum of the orders of the bonds to metals, or -1
/// if such a bond has no well-defined order.
pub fn n_bonds_val_to_metal(at: &[InpAtom], iat: usize) -> i32 {
    let a = &at[iat];
    let mut val_to_metal = 0;
    for i in 0..a.valence as usize {
        if is_el_a_metal(i32::from(at[a.neighbor[i] as usize].el_number)) {
            let bond_type = i32::from(a.bond_type[i]);
            if bond_type <= BOND_TYPE_TRIPLE {
                val_to_metal += bond_type;
            } else {
                return -1;
            }
        }
    }
    val_to_metal
}

/// `num_of_H`: explicit terminal H plus implicit H
pub fn num_of_h(at: &[InpAtom], iat: usize) -> i32 {
    let a = &at[iat];
    let explicit = (0..a.valence as usize)
        .filter(|&i| {
            let n = &at[a.neighbor[i] as usize];
            n.valence == 1 && n.el_number == EL_NUMBER_H
        })
        .count() as i32;
    explicit + a.numh()
}

/// `ion_el_group` (`FIX_REM_ION_PAIRS_Si_BUG == 0`)
pub fn ion_el_group(el: u8) -> u8 {
    match el {
        EL_NUMBER_C => EL_NUMBER_C,
        EL_NUMBER_N | EL_NUMBER_P | EL_NUMBER_AS | EL_NUMBER_SB => EL_NUMBER_N,
        EL_NUMBER_O | EL_NUMBER_S | EL_NUMBER_SE | EL_NUMBER_TE => EL_NUMBER_O,
        _ => 0,
    }
}

/// `has_other_ion_neigh`
pub fn has_other_ion_neigh(at: &[InpAtom], iat: usize, iat_ion_neigh: usize) -> bool {
    let charge = at[iat_ion_neigh].charge;
    (0..at[iat].valence as usize).any(|i| {
        let neigh = at[iat].neighbor[i] as usize;
        neigh != iat_ion_neigh
            && at[neigh].charge == charge
            && ion_el_group(at[neigh].el_number) != 0
    })
}

/// `has_other_ion_in_sphere_2`: BFS to radius 2; `c_flags` is used as the
/// visited mark and restored.
pub fn has_other_ion_in_sphere_2(at: &mut [InpAtom], iat: usize, iat_ion_neigh: usize) -> i32 {
    let mut q: Vec<usize> = vec![iat];
    at[iat].c_flags = 1;
    let mut iq = 0;
    let mut dist = 1;
    let mut ret = 0;
    while dist <= 2 {
        let lenq2 = q.len();
        while iq < lenq2 {
            let i = q[iq];
            for j in 0..at[i].valence as usize {
                let neigh = at[i].neighbor[j] as usize;
                if at[neigh].c_flags == 0
                    && at[neigh].valence <= 3
                    && ion_el_group(at[neigh].el_number) != 0
                {
                    q.push(neigh);
                    at[neigh].c_flags = 1;
                    if neigh != iat_ion_neigh && at[iat_ion_neigh].charge == at[neigh].charge {
                        ret += 1;
                    }
                }
            }
            iq += 1;
        }
        dist += 1;
    }
    for &i in &q {
        at[i].c_flags = 0;
    }
    ret
}

/// The bonds to metals of an atom whose valence exceeds the standard one (or
/// the `O(+)` case): `(count, sum of orders)`, `None` for the C code's
/// "fall back" return.
fn metal_bonds(at: &[InpAtom], a: &InpAtom) -> Option<(i32, i32)> {
    let mut count = 0;
    let mut sum = 0;
    for i in 0..a.valence as usize {
        if is_el_a_metal(i32::from(at[a.neighbor[i] as usize].el_number)) {
            if i32::from(a.bond_type[i]) & BOND_TYPE_MASK >= BOND_TYPE_ALTERN {
                return None;
            }
            count += 1;
            sum += i32::from(a.bond_type[i]) & BOND_TYPE_MASK;
        }
    }
    Some((count, sum))
}

/// `nNoMetalNumBonds`: the number of bonds not counting those to metals when
/// that gives the standard valence
pub fn n_no_metal_num_bonds(at: &[InpAtom], at_no: usize) -> i32 {
    let a = &at[at_no];
    let num_h = a.numh();
    let std = get_el_valence(i32::from(a.el_number), i32::from(a.charge), 0);
    let cbv = i32::from(a.chem_bonds_valence);
    let valence = i32::from(a.valence);
    if cbv + num_h > std {
        let Some((n, v)) = metal_bonds(at, a) else {
            return valence;
        };
        if cbv + num_h - v == std {
            return valence - n;
        }
    } else if S_VI_O_PLUS_METAL_FIX_BOND
        && a.charge == 1
        && get_endpoint_valence(a.el_number) == 2
        && cbv + num_h == std
    {
        let Some((n, v)) = metal_bonds(at, a) else {
            return valence;
        };
        if v == 1 {
            return valence - n;
        }
    }
    valence
}

/// `nNoMetalBondsValence`
pub fn n_no_metal_bonds_valence(at: &[InpAtom], at_no: usize) -> i32 {
    let a = &at[at_no];
    let num_h = a.numh();
    let std = get_el_valence(i32::from(a.el_number), i32::from(a.charge), 0);
    let cbv = i32::from(a.chem_bonds_valence);
    if cbv + num_h > std {
        let Some((_, v)) = metal_bonds(at, a) else {
            return i32::from(a.valence);
        };
        if cbv + num_h - v == std {
            return cbv - v;
        }
    } else if S_VI_O_PLUS_METAL_FIX_BOND
        && a.charge == 1
        && get_endpoint_valence(a.el_number) == 2
        && cbv + num_h == std
    {
        let Some((_, v)) = metal_bonds(at, a) else {
            return i32::from(a.valence);
        };
        if v == 1 {
            return cbv - v;
        }
    }
    cbv
}

/// `nNoMetalNeighIndex`: the index of the first non-metal neighbor, or -1
pub fn n_no_metal_neigh_index(at: &[InpAtom], at_no: usize) -> i32 {
    let a = &at[at_no];
    (0..a.valence as usize)
        .find(|&i| !is_el_a_metal(i32::from(at[a.neighbor[i] as usize].el_number)))
        .map_or(-1, |i| i as i32)
}

/// `nNoMetalOtherNeighIndex`
pub fn n_no_metal_other_neigh_index(at: &[InpAtom], at_no: usize, cur_neigh: i32) -> i32 {
    let a = &at[at_no];
    (0..a.valence as usize)
        .find(|&i| {
            let n = i32::from(a.neighbor[i]);
            n != cur_neigh && !is_el_a_metal(i32::from(at[n as usize].el_number))
        })
        .map_or(-1, |i| i as i32)
}

/// `nNoMetalOtherNeighIndex2`
pub fn n_no_metal_other_neigh_index2(
    at: &[InpAtom],
    at_no: usize,
    cur_neigh: i32,
    cur_neigh2: i32,
) -> i32 {
    let a = &at[at_no];
    (0..a.valence as usize)
        .find(|&i| {
            let n = i32::from(a.neighbor[i]);
            n != cur_neigh && n != cur_neigh2 && !is_el_a_metal(i32::from(at[n as usize].el_number))
        })
        .map_or(-1, |i| i as i32)
}

// ---- strings and text handling ----

/// C `isspace` in the "C" locale
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | 0x09..=0x0d)
}

/// The part of a buffer before the first NUL, like a C string.
fn c_str(s: &[u8]) -> &[u8] {
    &s[..s.iter().position(|&c| c == 0).unwrap_or(s.len())]
}

/// `MakeRemovedProtonsString`: the text and the number of isotopic H kinds
/// removed (`num_removed_iso_H`, only set when there are any).
pub fn make_removed_protons_string(
    n_num_removed_protons: i32,
    n_num_exchg_isotopic_h: Option<&[i32; NUM_H_ISOTOPES]>,
    n_num_removed_protons_isotopic: Option<&[i32; NUM_H_ISOTOPES]>,
    b_isotopic: bool,
) -> (String, Option<i32>) {
    let mut s = String::new();
    let mut num_iso = None;
    if n_num_removed_protons != 0 {
        s = format!(
            "Proton balance: {} {} H+",
            if n_num_removed_protons >= 0 { '+' } else { '-' },
            n_num_removed_protons.abs()
        );
    }
    if b_isotopic && (n_num_removed_protons_isotopic.is_some() || n_num_exchg_isotopic_h.is_some())
    {
        let mut j = 0;
        for i in 0..NUM_H_ISOTOPES {
            let num = n_num_exchg_isotopic_h.map_or(0, |a| a[i])
                + n_num_removed_protons_isotopic.map_or(0, |a| a[i]);
            if num != 0 {
                s += &format!(
                    "{} {}^{}H",
                    if j != 0 { ", " } else { "  [ removed " },
                    num,
                    i + 1
                );
                j += 1;
            }
        }
        if j != 0 {
            s += " ]";
            num_iso = Some(j);
        }
    }
    (s, num_iso)
}

/// `normalize_string`: removes leading and trailing white space and replaces
/// runs of white space with a single space; returns the new length.
pub fn normalize_string(name: &mut Vec<u8>) -> i32 {
    let mut len = name.len() as i32;
    let mut n: i32 = 0;
    let mut i: i32 = 0;
    while i < len {
        if is_c_space(name[i as usize]) {
            name[i as usize] = b' ';
            n += 1;
        } else {
            if n > 0 {
                name.copy_within((i as usize).., (i - n) as usize);
                i -= n;
                len -= n;
            }
            n = -1;
        }
        i += 1;
    }
    if n == len {
        len = 0;
    } else {
        n += 1;
        if n != 0 && n <= len {
            len -= n;
        }
    }
    name.truncate(len as usize);
    len
}

/// `dotify_non_printable_chars`: replaces non-ASCII and control characters
/// with '.', returns the number of replacements
pub fn dotify_non_printable_chars(line: &mut [u8]) -> i32 {
    let mut num = 0;
    for c in line.iter_mut() {
        if *c == 0 {
            break;
        }
        if *c < b' ' || *c >= 0x7f {
            *c = b'.';
            num += 1;
        }
    }
    num
}

/// `is_matching_any_delim`
pub fn is_matching_any_delim(c: u8, delims: &[u8]) -> bool {
    c_str(delims).contains(&c)
}

/// `read_upto_delim`: `s` is the C string and `*pstring` the reading
/// position in it (`None` = NULL). Copies up to the next delimiter or the end
/// to `field`; the position is moved to the delimiter (or `None` at the end).
/// Returns the number of bytes copied, -1 on error.
pub fn read_upto_delim(
    s: &[u8],
    pstring: &mut Option<usize>,
    field: &mut Vec<u8>,
    maxlen: i32,
    delims: &[u8],
) -> i32 {
    let Some(start) = *pstring else {
        return -1;
    };
    let s = c_str(s);
    let at = |k: usize| s.get(k).copied().unwrap_or(0);
    let mut i = 0;
    while at(start + i) != 0 && is_c_space(at(start + i)) {
        i += 1;
    }
    let p = start + i;
    let mut n = 0;
    while at(p + n) != 0 && !is_matching_any_delim(at(p + n), delims) {
        n += 1;
    }
    if n as i32 + 1 > maxlen {
        return -1;
    }
    *field = s[p..p + n].to_vec();
    *pstring = if at(p + n) == 0 {
        None
    } else {
        Some(start + i + n)
    };
    n as i32
}

/// `remove_trailing_spaces`
pub fn remove_trailing_spaces(p: &mut Vec<u8>) {
    while p.last().is_some_and(|&c| is_c_space(c)) {
        p.pop();
    }
}

/// `remove_one_lf`: removes a final "\n" or "\r\n"
pub fn remove_one_lf(p: &mut Vec<u8>) {
    let len = p.len();
    if len > 0 && p[len - 1] == b'\n' {
        p.truncate(len - 1);
        if len >= 2 && p[len - 2] == b'\r' {
            p.truncate(len - 2);
        }
    }
}

/// `mystrncpy`: copies up to `maxlen` bytes including the terminating zero
/// from `source` to the start of `target` and zero-fills up to `maxlen`.
/// Returns false (C: 0) for empty input.
pub fn mystrncpy(target: &mut [u8], source: &[u8], maxlen: usize) -> bool {
    if maxlen == 0 {
        return false;
    }
    let window = &source[..maxlen.min(source.len())];
    let len = window.iter().position(|&c| c == 0).unwrap_or(maxlen - 1);
    let len = len.min(window.len());
    target[..len].copy_from_slice(&source[..len]);
    target[len..maxlen].fill(0);
    true
}

/// `lrtrim`: removes leading and trailing white space; returns the length
pub fn lrtrim(p: &mut Vec<u8>) -> i32 {
    let is_sp = |c: u8| c.is_ascii() && is_c_space(c);
    let lead = p.iter().take_while(|&&c| is_sp(c)).count();
    p.drain(..lead);
    while p.last().is_some_and(|&c| is_sp(c)) {
        p.pop();
    }
    p.len() as i32
}

/// `extract_inchi_substring`: the InChI embedded in `s` (from "InChI=" up to
/// the first character an InChI cannot contain, at most `slen` characters)
pub fn extract_inchi_substring(s: &str, slen: usize) -> Option<String> {
    let p = &s[s.find("InChI=")?..];
    let b = p.as_bytes();
    let mut i = 0;
    while i < slen {
        let c = b.get(i).copied().unwrap_or(0);
        if !(c.is_ascii_alphanumeric() || b"()*+,-./;=?@".contains(&c)) {
            break;
        }
        i += 1;
    }
    Some(p[..i.min(b.len())].to_string())
}

/// `extract_auxinfo_substring`: from "AuxInfo=" to the first white space
pub fn extract_auxinfo_substring(s: &str, slen: usize) -> Option<String> {
    let p = &s[s.find("AuxInfo=")?..];
    let b = p.as_bytes();
    let mut i = 0;
    while i < slen {
        if b.get(i).is_some_and(|&c| is_c_space(c)) {
            break;
        }
        i += 1;
    }
    Some(p[..i.min(b.len())].to_string())
}

fn my_tolower(c: u8) -> i32 {
    i32::from(c.to_ascii_lowercase())
}

/// `inchi_memicmp`
pub fn inchi_memicmp(p1: &[u8], p2: &[u8], length: usize) -> i32 {
    for i in 0..length {
        if p1[i] != p2[i] && my_tolower(p1[i]) != my_tolower(p2[i]) {
            return my_tolower(p1[i]) - my_tolower(p2[i]);
        }
    }
    0
}

/// `inchi_stricmp`
pub fn inchi_stricmp(s1: &str, s2: &str) -> i32 {
    let (a, b) = (s1.as_bytes(), s2.as_bytes());
    let mut i = 0;
    while i < a.len() {
        let c2 = b.get(i).copied().unwrap_or(0);
        if a[i] != c2 && my_tolower(a[i]) != my_tolower(c2) {
            return my_tolower(a[i]) - my_tolower(c2);
        }
        i += 1;
    }
    if i < b.len() { -1 } else { 0 }
}

/// `inchi__strnset`: sets up to `length` characters before the terminating
/// zero to `val`
pub fn inchi_strnset(s: &mut [u8], val: u8, length: usize) {
    for c in s.iter_mut().take(length) {
        if *c == 0 {
            break;
        }
        *c = val;
    }
}

/// `get_endpoint_valence_KET` (`KETO_ENOL_TAUT == 1`): only O and C can be an
/// endpoint of keto-enol tautomerism
pub fn get_endpoint_valence_ket(el_number: u8) -> i32 {
    match el_number {
        EL_NUMBER_C => 4,
        EL_NUMBER_O => 2,
        _ => 0,
    }
}
