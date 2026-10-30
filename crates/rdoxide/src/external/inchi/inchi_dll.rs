//! The input side of the InChI library (`INCHI_API/libinchi/src/inchi_dll.c`):
//! conversion of the user's `inchi_Atom` array into `inp_ATOM`s —
//! `SetAtomProperties`, `SetBondProperties`, `SetAtomAndBondProperties`,
//! `SetNumImplicitH` and the loop of `ExtractOneStructure` that calls them —
//! and `WriteCoord` (`ichimak2.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichierr::add_error_message;
use super::incomdef::{
    BOND_TYPE_ALTERN, BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, NUM_H_ISOTOPES,
    RADICAL_DOUBLET, RADICAL_TRIPLET,
};
use super::inpdef::{ATOM_EL_LEN, InpAtom, MAXVAL};
use super::util::{
    EL_NUMBER_H, ERR_ELEM, detect_unusual_el_valence, extract_charges_and_radicals,
    extract_h_atoms, get_atomic_mass_from_elnum, get_num_h, get_periodic_table_number,
    is_el_a_metal, mystrncpy, n_bonds_val_to_metal,
};

// inchi_api.h
pub const INCHI_RADICAL_NONE: i8 = 0;
pub const INCHI_RADICAL_SINGLET: i8 = 1;
pub const INCHI_RADICAL_DOUBLET: i8 = 2;
pub const INCHI_RADICAL_TRIPLET: i8 = 3;

pub const INCHI_BOND_TYPE_NONE: i8 = 0;
pub const INCHI_BOND_TYPE_SINGLE: i8 = 1;
pub const INCHI_BOND_TYPE_DOUBLE: i8 = 2;
pub const INCHI_BOND_TYPE_TRIPLE: i8 = 3;
pub const INCHI_BOND_TYPE_ALTERN: i8 = 4;

pub const INCHI_BOND_STEREO_NONE: i8 = 0;
pub const INCHI_BOND_STEREO_SINGLE_1UP: i8 = 1;
pub const INCHI_BOND_STEREO_SINGLE_1EITHER: i8 = 4;
pub const INCHI_BOND_STEREO_SINGLE_1DOWN: i8 = 6;
pub const INCHI_BOND_STEREO_SINGLE_2UP: i8 = -1;
pub const INCHI_BOND_STEREO_SINGLE_2EITHER: i8 = -4;
pub const INCHI_BOND_STEREO_SINGLE_2DOWN: i8 = -6;
pub const INCHI_BOND_STEREO_DOUBLE_EITHER: i8 = 3;

pub const ISOTOPIC_SHIFT_FLAG: i32 = 10000;
pub const ISOTOPIC_SHIFT_MAX: i32 = 100;

// incomdef.h
pub const STEREO_SNGL_UP: i8 = 1;
pub const STEREO_SNGL_EITHER: i8 = 4;
pub const STEREO_SNGL_DOWN: i8 = 6;
pub const STEREO_DBLE_EITHER: i8 = 3;

/// `inchi_Atom` (`AT_NUM` is a signed 16 bit integer)
#[derive(Debug, Clone, PartialEq)]
pub struct InchiAtom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub neighbor: [i16; MAXVAL],
    pub bond_type: [i8; MAXVAL],
    pub bond_stereo: [i8; MAXVAL],
    pub elname: [u8; ATOM_EL_LEN],
    pub num_bonds: i16,
    /// `[0]`: implicit non-isotopic H (-1: add automatically), `[1..]`: 1H, 2H, 3H
    pub num_iso_h: [i8; NUM_H_ISOTOPES + 1],
    pub isotopic_mass: i16,
    pub radical: i8,
    pub charge: i8,
}

impl Default for InchiAtom {
    fn default() -> Self {
        InchiAtom {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            neighbor: [0; MAXVAL],
            bond_type: [0; MAXVAL],
            bond_stereo: [0; MAXVAL],
            elname: [0; ATOM_EL_LEN],
            num_bonds: 0,
            num_iso_h: [0; NUM_H_ISOTOPES + 1],
            isotopic_mass: 0,
            radical: 0,
            charge: 0,
        }
    }
}

impl InchiAtom {
    pub fn elname_str(&self) -> String {
        let n = self
            .elname
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(ATOM_EL_LEN);
        String::from_utf8_lossy(&self.elname[..n]).to_string()
    }
}

/// `MOL_COORD`: three 10-character numbers
pub use super::ichiprt1::MolCoord;

const MIN_BOND_LENGTH: f64 = 1.0e-6;
/// `SINGLET_IS_TRIPLET`
const SINGLET_IS_TRIPLET: bool = true;

/// C's `%.{prec}e`
fn c_exp(x: f64, prec: usize) -> String {
    let s = format!("{x:.prec$e}");
    let (mant, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let e: i32 = exp.parse().unwrap_or(0);
    format!("{mant}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
}

/// `format!("{x:.prec$}")` for moderate values without the exact big-number path of the standard formatter: the
/// scaled value is rounded as an integer unless it lies too close to a rounding boundary to be sure
fn fixed(x: f64, prec: usize) -> String {
    const POW: [f64; 5] = [1.0, 10.0, 100.0, 1000.0, 10000.0];
    let v = x.abs() * POW[prec];
    if v < 1e15 {
        let r = v.round();
        if (v - r).abs() < 0.499 {
            let n = r as u64;
            let p = POW[prec] as u64;
            let sign = if x.is_sign_negative() { "-" } else { "" };
            return if prec == 0 {
                format!("{sign}{n}")
            } else {
                format!("{sign}{}.{:0width$}", n / p, n % p, width = prec)
            };
        }
    }
    format!("{x:.prec$}")
}

/// `WriteCoord`: a coordinate in 10 characters (or a little more for huge
/// values)
pub fn write_coord(x: f64) -> String {
    let s = if x.is_nan() || x.is_infinite() {
        let t = if x.is_nan() {
            "nan"
        } else if x < 0.0 {
            "-inf"
        } else {
            "inf"
        };
        t.to_string()
    } else if x < -9999999.9 {
        c_exp(x, 2)
    } else if x < -999999.99 {
        fixed(x, 2)
    } else if x < -99999.999 {
        fixed(x, 3)
    } else if x < 99999.9999 {
        fixed(x, 4)
    } else if x < 999999.999 {
        fixed(x, 3)
    } else if x < 9999999.99 {
        fixed(x, 2)
    } else if x < 99999999.9 {
        fixed(x, 1)
    } else {
        c_exp(x, 3)
    };
    format!("{s:>10}")
}

/// `SetAtomProperties` (the coordinate text is stored when `sz_coord` is
/// given). `err` is the error bit mask, `p_str_err` the message buffer.
pub fn set_atom_properties(
    at: &mut [InpAtom],
    sz_coord: Option<&mut [MolCoord]>,
    ati: &[InchiAtom],
    a1: usize,
    n_dim: &mut i32,
    p_str_err: &mut String,
    err: &mut i32,
) -> i32 {
    // element, check later
    at[a1].elname = ati[a1].elname;
    at[a1].charge = ati[a1].charge;

    let c_radical: i8 = match ati[a1].radical {
        INCHI_RADICAL_NONE => 0,
        INCHI_RADICAL_SINGLET => {
            if SINGLET_IS_TRIPLET {
                RADICAL_TRIPLET as i8
            } else {
                1
            }
        }
        INCHI_RADICAL_DOUBLET => RADICAL_DOUBLET as i8,
        INCHI_RADICAL_TRIPLET => RADICAL_TRIPLET as i8,
        r => {
            let mut n_rad = i32::from(r);
            while n_rad > RADICAL_TRIPLET {
                n_rad -= 2;
            }
            let radical_type = format!("{}->{}", r, n_rad);
            add_error_message(p_str_err, "Radical center type replaced:");
            add_error_message(p_str_err, &radical_type);
            if n_rad < 0 {
                *err |= 8;
            }
            n_rad as i8
        }
    };
    at[a1].radical = c_radical;

    // coordinates
    at[a1].x = ati[a1].x;
    at[a1].y = ati[a1].y;
    at[a1].z = ati[a1].z;

    if let Some(coords) = sz_coord {
        let c = &mut coords[a1];
        for (k, v) in [ati[a1].x, ati[a1].y, ati[a1].z].into_iter().enumerate() {
            let s = write_coord(v);
            c[10 * k..10 * k + 10].copy_from_slice(&s.as_bytes()[..10]);
        }
    }

    if MIN_BOND_LENGTH < ati[a1].x.abs()
        || MIN_BOND_LENGTH < ati[a1].y.abs()
        || MIN_BOND_LENGTH < ati[a1].z.abs()
    {
        if MIN_BOND_LENGTH < ati[a1].z.abs() {
            *n_dim |= 3;
        } else {
            *n_dim |= 2;
        }
    }

    // orig. at. number
    at[a1].orig_at_number = (a1 + 1) as u16;
    0
}

fn is_in_the_list(list: &[u16], next: u16, len: usize) -> Option<usize> {
    list[..len].iter().position(|&a| a == next)
}

/// `SetBondProperties`: returns 1 if the bond was ignored
pub fn set_bond_properties(
    at: &mut [InpAtom],
    ati: &[InchiAtom],
    a1: usize,
    j: usize,
    n_num_atoms: usize,
    n_num_bonds: &mut i32,
    p_str_err: &mut String,
    err: &mut i32,
) -> i32 {
    // bond type
    let c_bond_type: u8 = match ati[a1].bond_type[j] {
        INCHI_BOND_TYPE_SINGLE => BOND_TYPE_SINGLE as u8,
        INCHI_BOND_TYPE_DOUBLE => BOND_TYPE_DOUBLE as u8,
        INCHI_BOND_TYPE_TRIPLE => BOND_TYPE_TRIPLE as u8,
        INCHI_BOND_TYPE_ALTERN => BOND_TYPE_ALTERN as u8,
        other => {
            add_error_message(p_str_err, "Unrecognized bond type:");
            add_error_message(p_str_err, &other.to_string());
            *err |= 8;
            BOND_TYPE_SINGLE as u8
        }
    };

    // 2D stereo
    let (c_stereo1, c_stereo2): (i8, i8) = match ati[a1].bond_stereo[j] {
        INCHI_BOND_STEREO_NONE => (0, 0),
        INCHI_BOND_STEREO_SINGLE_1UP => (STEREO_SNGL_UP, -STEREO_SNGL_UP),
        INCHI_BOND_STEREO_SINGLE_1EITHER => (STEREO_SNGL_EITHER, -STEREO_SNGL_EITHER),
        INCHI_BOND_STEREO_SINGLE_1DOWN => (STEREO_SNGL_DOWN, -STEREO_SNGL_DOWN),
        INCHI_BOND_STEREO_SINGLE_2UP => (-STEREO_SNGL_UP, STEREO_SNGL_UP),
        INCHI_BOND_STEREO_SINGLE_2EITHER => (-STEREO_SNGL_EITHER, STEREO_SNGL_EITHER),
        INCHI_BOND_STEREO_SINGLE_2DOWN => (-STEREO_SNGL_DOWN, STEREO_SNGL_DOWN),
        3 | -3 => (STEREO_DBLE_EITHER, STEREO_DBLE_EITHER),
        other => {
            add_error_message(p_str_err, "Unrecognized bond stereo:");
            add_error_message(p_str_err, &other.to_string());
            *err |= 8;
            (0, 0)
        }
    };

    // neighbor
    if ati[a1].neighbor[j] < 0 || ati[a1].neighbor[j] as usize >= n_num_atoms {
        *err |= 1;
        add_error_message(p_str_err, "Bond to nonexistent atom");
        return 1;
    }
    let a2 = ati[a1].neighbor[j] as usize;
    if a2 == a1 {
        *err |= 1;
        add_error_message(p_str_err, "Atom has a bond to itself");
        return 1;
    }

    // consistency check; locate the bond in the opposite atom
    let v1 = at[a1].valence as usize;
    let v2 = at[a2].valence as usize;
    let p1 = is_in_the_list(&at[a1].neighbor, a2 as u16, v1);
    let p2 = is_in_the_list(&at[a2].neighbor, a1 as u16, v2);
    let (n1, n2);

    if let (Some(i1), Some(i2)) = (p1, p2) {
        n1 = i1;
        n2 = i2;
        let dup_later = (n1 + 1 < v1
            && is_in_the_list(&at[a1].neighbor[n1 + 1..], a2 as u16, v1 - n1 - 1).is_some())
            || (n2 + 1 < v2
                && is_in_the_list(&at[a2].neighbor[n2 + 1..], a1 as u16, v2 - n2 - 1).is_some());
        if dup_later {
            add_error_message(p_str_err, "Multiple bonds between two atoms");
            *err |= 2;
        } else if n1 < v1
            && n2 < v2
            && c_bond_type == at[a2].bond_type[n2]
            && c_bond_type == at[a1].bond_type[n1]
            && c_stereo1 == at[a1].bond_stereo[n1]
            && c_stereo2 == at[a2].bond_stereo[n2]
        {
            // duplicated bond: not reported
        } else {
            add_error_message(p_str_err, "Multiple bonds between two atoms");
            *err |= 2;
        }
    } else if (p1.is_some() || p2.is_some())
        && (p1.is_some() || v1 < MAXVAL)
        && (p2.is_some() || v2 < MAXVAL)
    {
        n1 = p1.unwrap_or_else(|| {
            let k = at[a1].valence as usize;
            at[a1].valence += 1;
            k
        });
        n2 = p2.unwrap_or_else(|| {
            let k = at[a2].valence as usize;
            at[a2].valence += 1;
            k
        });
        // the bond is present in one atom only: possibly program error
        if (p1.is_some()
            && (c_bond_type != at[a1].bond_type[n1] || at[a1].bond_stereo[n1] != c_stereo1))
            || (p2.is_some()
                && (c_bond_type != at[a2].bond_type[n2] || at[a2].bond_stereo[n2] != c_stereo2))
        {
            add_error_message(p_str_err, "Multiple bonds between two atoms");
            *err |= 2;
        } else {
            add_error_message(p_str_err, "Duplicated bond(s) between two atoms");
        }
    } else if p1.is_none() && p2.is_none() && v1 < MAXVAL && v2 < MAXVAL {
        n1 = at[a1].valence as usize;
        at[a1].valence += 1;
        n2 = at[a2].valence as usize;
        at[a2].valence += 1;
        *n_num_bonds += 1;
    } else {
        *err |= 4;
        let name = if v1 >= MAXVAL {
            at[a1].elname_str()
        } else {
            at[a2].elname_str()
        };
        add_error_message(
            p_str_err,
            &format!("Atom '{name}' has more than {MAXVAL} bonds"),
        );
        return 1;
    }

    // store the connection
    at[a1].bond_type[n1] = c_bond_type;
    at[a2].bond_type[n2] = c_bond_type;
    at[a1].neighbor[n1] = a2 as u16;
    at[a2].neighbor[n2] = a1 as u16;
    at[a1].bond_stereo[n1] = c_stereo1;
    at[a2].bond_stereo[n2] = c_stereo2;
    0
}

/// `SetAtomAndBondProperties`: element number, chemical valence, implicit H
/// flags and isotopic shifts. `at_type` gets temporary flags (1: data in the
/// atom name; 2: request to add H).
pub fn set_atom_and_bond_properties(
    at: &mut [InpAtom],
    ati: &[InchiAtom],
    a1: usize,
    b_do_not_add_h: bool,
    p_str_err: &mut String,
    err: &mut i32,
) -> i32 {
    let valence = at[a1].valence as usize;
    let mut chem_valence: i32 = 0;
    let mut num_alt_bonds = 0;
    for j in 0..valence {
        if i32::from(at[a1].bond_type[j]) <= BOND_TYPE_TRIPLE {
            chem_valence += i32::from(at[a1].bond_type[j]);
        } else {
            num_alt_bonds += 1;
        }
    }
    match num_alt_bonds {
        0 => {}
        2 => chem_valence += 3,
        3 => chem_valence += 4,
        _ => {
            *err |= 8;
            add_error_message(
                p_str_err,
                &format!(
                    "Atom '{}' has {} alternating bonds",
                    at[a1].elname_str(),
                    num_alt_bonds
                ),
            );
        }
    }
    at[a1].chem_bonds_valence = chem_valence as i8;

    // aliased hydrogen atoms
    let mut n1 = get_periodic_table_number(at[a1].elname_str());
    if n1 == ERR_ELEM {
        // elname contains more than 1 element: extract number of H if possible
        let mut name: Vec<u8> = at[a1].elname_str().to_string().into_bytes();
        let (n_radical, n_charge) = extract_charges_and_radicals(&mut name);
        if n_radical != 0 || n_charge != 0 {
            if (n_radical != 0 && at[a1].radical != 0 && n_radical != i32::from(at[a1].radical))
                || (n_charge != 0 && at[a1].charge != 0 && n_charge != i32::from(at[a1].charge))
            {
                add_error_message(p_str_err, "Ignored charge/radical redefinition:");
                add_error_message(p_str_err, &ati[a1].elname_str());
            } else {
                if n_radical != 0 {
                    at[a1].radical = n_radical as i8;
                }
                if n_charge != 0 {
                    at[a1].charge = n_charge as i8;
                }
            }
        }
        let mut num_iso = at[a1].num_iso_h;
        at[a1].num_h = extract_h_atoms(&mut name, &mut num_iso) as i8;
        at[a1].num_iso_h = num_iso;
        at[a1].elname = [0; ATOM_EL_LEN];
        for (d, s) in at[a1]
            .elname
            .iter_mut()
            .zip(name.iter().take(ATOM_EL_LEN - 1))
        {
            *d = *s;
        }
        if at[a1].elname[0] == 0 && at[a1].numh() != 0 {
            // alias contains only H: move the heaviest isotope to the "central atom"
            at[a1].set_elname("H");
            if at[a1].num_iso_h_total() != 0 {
                for j in (0..NUM_H_ISOTOPES).rev() {
                    if at[a1].num_iso_h[j] != 0 {
                        at[a1].num_iso_h[j] -= 1;
                        at[a1].iso_atw_diff = (1 + j) as i8;
                        break;
                    }
                }
            } else {
                at[a1].num_h -= 1;
            }
        }
        n1 = get_periodic_table_number(at[a1].elname_str());
        if n1 == ERR_ELEM {
            n1 = 0;
        }
        if n1 != 0 {
            at[a1].at_type |= 1;
            add_error_message(p_str_err, "Parsed compound atom(s):");
            add_error_message(p_str_err, &ati[a1].elname_str());
        }
    }

    at[a1].el_number = n1 as u8;
    if n1 == 0 {
        *err |= 64;
        add_error_message(p_str_err, "Unknown element(s):");
        add_error_message(p_str_err, at[a1].elname_str());
    } else if n1 == i32::from(EL_NUMBER_H) && at[a1].iso_atw_diff == 0 {
        // replace explicit D or T with isotopic H
        match at[a1].elname[0] {
            b'D' => {
                at[a1].iso_atw_diff = 2;
                let mut e = vec![0u8; ATOM_EL_LEN];
                mystrncpy(&mut e, b"H", ATOM_EL_LEN);
                at[a1].elname.copy_from_slice(&e);
            }
            b'T' => {
                at[a1].iso_atw_diff = 3;
                let mut e = vec![0u8; ATOM_EL_LEN];
                mystrncpy(&mut e, b"H", ATOM_EL_LEN);
                at[a1].elname.copy_from_slice(&e);
            }
            b'H' if ati[a1].isotopic_mass >= 1 => {
                // reproduce Bug04: allowed non-terminal H heavier than T
                let diff16 = shift_i16(ati[a1].isotopic_mass, at[a1].el_number);
                if diff16 >= 1
                    && (at[a1].valence != 1 || i32::from(diff16) <= NUM_H_ISOTOPES as i32)
                {
                    at[a1].iso_atw_diff = diff16 as i8;
                }
            }
            _ => {}
        }
    } else if ati[a1].isotopic_mass != 0 {
        // isotopic shift
        at[a1].iso_atw_diff = shift_i16(ati[a1].isotopic_mass, at[a1].el_number) as i8;
    }

    // add implicit hydrogen atoms flag
    if ati[a1].num_iso_h[0] == -1 {
        if !b_do_not_add_h {
            at[a1].at_type |= 2;
        }
    } else {
        at[a1].num_h = ati[a1].num_iso_h[0];
    }
    for j in 0..NUM_H_ISOTOPES {
        at[a1].num_iso_h[j] = ati[a1].num_iso_h[j + 1];
    }

    if num_alt_bonds != 0 {
        // atom has aromatic bonds AND the chemical valence is not known
        let num_h = at[a1].numh();
        let chem_valence_alt = i32::from(at[a1].chem_bonds_valence) + num_h;
        let unusual_arom = detect_unusual_el_valence(
            i32::from(at[a1].el_number),
            i32::from(at[a1].charge),
            i32::from(at[a1].radical),
            chem_valence_alt,
            num_h,
            i32::from(at[a1].valence),
        );
        let unusual_no_arom = detect_unusual_el_valence(
            i32::from(at[a1].el_number),
            i32::from(at[a1].charge),
            i32::from(at[a1].radical),
            chem_valence_alt - 1,
            num_h,
            i32::from(at[a1].valence),
        );
        if unusual_arom != 0 && unusual_no_arom == 0 && n_bonds_val_to_metal(at, a1) == 0 {
            // typically NH in a 5-member aromatic ring
            at[a1].chem_bonds_valence -= 1;
        }
    }
    0
}

/// `AT_NUM iso_atw_diff` of `SetAtomAndBondProperties` (a signed 16 bit
/// value) after the increment of a non-negative shift
fn shift_i16(isotopic_mass: i16, el_number: u8) -> i16 {
    let m = i32::from(isotopic_mass);
    let mut diff: i16 = if (ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX
        ..=ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX)
        .contains(&m)
    {
        (m - ISOTOPIC_SHIFT_FLAG) as i16
    } else {
        (m - get_atomic_mass_from_elnum(i32::from(el_number))) as i16
    };
    if diff >= 0 {
        diff += 1;
    }
    diff
}

/// `SetNumImplicitH`: sets the number of implicit H of all atoms, metals
/// first, and clears the temporary `at_type` flags.
pub fn set_num_implicit_h(at: &mut [InpAtom]) {
    for b_non_metal in 0..2 {
        for a1 in 0..at.len() {
            if b_non_metal != i32::from(is_el_a_metal(i32::from(at[a1].el_number))) {
                continue;
            }
            let iso = at[a1].num_iso_h;
            at[a1].num_h = get_num_h(
                at[a1].elname_str(),
                i32::from(at[a1].num_h),
                Some(&iso),
                i32::from(at[a1].charge),
                i32::from(at[a1].radical),
                i32::from(at[a1].chem_bonds_valence),
                0,
                at[a1].at_type & 1 != 0,
                at[a1].at_type & 2 == 0,
                false,
            ) as i8;
            at[a1].at_type = 0;
        }
    }
}

/// The atom/bond loops of `ExtractOneStructure`: returns the `inp_ATOM`s,
/// the number of dimensions and bonds. `err` and `p_str_err` are the
/// structure's error bit mask and message; processing stops at the first
/// atom after which `err` is non-zero (the atoms are then incomplete).
pub fn extract_atoms_and_bonds(
    ati: &[InchiAtom],
    b_do_not_add_h: bool,
    p_str_err: &mut String,
    err: &mut i32,
) -> (Vec<InpAtom>, i32, i32, Vec<MolCoord>) {
    let n = ati.len();
    let mut at = vec![InpAtom::default(); n];
    let mut sz_coord = vec![[0u8; super::ichiprt1::MOL_COORD_LEN]; n.max(1)];
    let mut n_dim = 0;
    let mut n_num_bonds = 0;
    for a1 in 0..n {
        set_atom_properties(
            &mut at,
            Some(&mut sz_coord),
            ati,
            a1,
            &mut n_dim,
            p_str_err,
            err,
        );
        if *err != 0 {
            return (at, n_dim, n_num_bonds, sz_coord);
        }
        for j in 0..ati[a1].num_bonds.max(0) as usize {
            set_bond_properties(&mut at, ati, a1, j, n, &mut n_num_bonds, p_str_err, err);
        }
        if *err != 0 {
            return (at, n_dim, n_num_bonds, sz_coord);
        }
    }
    for a1 in 0..n {
        set_atom_and_bond_properties(&mut at, ati, a1, b_do_not_add_h, p_str_err, err);
        if *err != 0 {
            return (at, n_dim, n_num_bonds, sz_coord);
        }
    }
    set_num_implicit_h(&mut at);
    (at, n_dim, n_num_bonds, sz_coord)
}

use super::inpdef::MAX_NUM_STEREO_ATOM_NEIGH;
use super::inpdef::MAX_NUM_STEREO_BONDS;
use super::readinch::{
    INCHI_PARITY_EVEN, INCHI_PARITY_ODD, INCHI_STEREO_TYPE_ALLENE, INCHI_STEREO_TYPE_DOUBLE_BOND,
    INCHI_STEREO_TYPE_TETRAHEDRAL, InchiStereo0D, NO_ATOM,
};

/// `InpAtom0DToInchiAtom`: converts `inp_ATOM`s with 0D parities back to the
/// API structures. Returns the result code (bit 1: more stereo descriptors
/// than were counted) and the atoms and 0D stereo descriptors.
pub fn inp_atom_0d_to_inchi_atom(at: &[InpAtom]) -> (i32, Vec<InchiAtom>, Vec<InchiStereo0D>) {
    let num_inp_atoms = at.len();
    let mut ret = 0;
    let mut num_stereo_centers = 0usize;
    let mut num_stereo_bonds = 0usize;
    // count stereobonds, allenes, cumulenes and stereoatoms
    for a in at {
        if a.p_parity != 0 {
            num_stereo_centers += 1;
        } else {
            let mut m = 0;
            while m < MAX_NUM_STEREO_BONDS && a.sb_parity[m] != 0 {
                m += 1;
            }
            num_stereo_bonds += m;
        }
    }
    num_stereo_bonds /= 2;
    let num_inp_stereo_0d = num_stereo_bonds + num_stereo_centers;

    let mut atom = vec![InchiAtom::default(); num_inp_atoms];
    let mut stereo_0d: Vec<InchiStereo0D> = Vec::new();

    // copy atom properties
    for (i, a) in at.iter().enumerate() {
        atom[i].num_bonds = i16::from(a.valence);
        for m in 0..a.valence as usize {
            atom[i].bond_type[m] = a.bond_type[m] as i8;
            atom[i].neighbor[m] = a.neighbor[m] as i16;
        }
        atom[i].charge = a.charge;
        atom[i].elname = a.elname;
        if a.iso_atw_diff != 0 {
            atom[i].isotopic_mass = (ISOTOPIC_SHIFT_FLAG
                + if a.iso_atw_diff > 0 {
                    i32::from(a.iso_atw_diff) - 1
                } else {
                    i32::from(a.iso_atw_diff)
                }) as i16;
        }
        atom[i].num_iso_h[0] = a.num_h;
        for m in 0..NUM_H_ISOTOPES {
            atom[i].num_iso_h[m + 1] = a.num_iso_h[m];
        }
        atom[i].radical = a.radical;
    }

    // stereo
    'atoms: for i in 0..num_inp_atoms {
        if at[i].p_parity != 0 {
            if stereo_0d.len() < num_inp_stereo_0d {
                let mut neighbor = [0i16; 4];
                for m in 0..MAX_NUM_STEREO_ATOM_NEIGH {
                    neighbor[m] = (i32::from(at[i].p_orig_at_num[m]) - 1) as i16;
                }
                stereo_0d.push(InchiStereo0D {
                    neighbor,
                    central_atom: i as i16,
                    type_: INCHI_STEREO_TYPE_TETRAHEDRAL,
                    parity: at[i].p_parity,
                });
            } else {
                ret |= 1;
                break 'atoms;
            }
        } else {
            let mut m1 = 0;
            while m1 < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m1] != 0 {
                // find the opposite atom at the other end of a double bond,
                // allene, or cumulene
                let mut chain = [0usize; 12];
                let mut len = 0usize;
                let mut cur = i;
                chain[len] = i;
                len += 1;
                let mut nxt_neigh = at[cur].sb_ord[m1] as usize;
                loop {
                    let nxt = at[cur].neighbor[nxt_neigh] as usize;
                    chain[len] = nxt;
                    len += 1;
                    nxt_neigh = usize::from(at[nxt].neighbor[0] as usize == cur);
                    cur = nxt;
                    if !(at[cur].sb_parity[0] == 0 && len < 12 && at[cur].valence == 2) {
                        break;
                    }
                }
                if at[cur].sb_parity[0] != 0 && len <= 4 && i < cur {
                    // double bond, cumulene, or allene has been found
                    let mut m2 = 0;
                    while m2 < MAX_NUM_STEREO_BONDS && at[cur].sb_parity[m2] != 0 {
                        if chain[len - 2] == at[cur].neighbor[at[cur].sb_ord[m2] as usize] as usize
                        {
                            if stereo_0d.len() < num_inp_stereo_0d {
                                let parity1 = i32::from(at[i].sb_parity[m1]);
                                let parity2 = i32::from(at[cur].sb_parity[m2]);
                                let parity = if (parity1 == INCHI_PARITY_ODD
                                    || parity1 == INCHI_PARITY_EVEN)
                                    && (parity2 == INCHI_PARITY_ODD || parity2 == INCHI_PARITY_EVEN)
                                {
                                    // well-defined parity
                                    if parity1 == parity2 {
                                        INCHI_PARITY_EVEN
                                    } else {
                                        INCHI_PARITY_ODD
                                    }
                                } else {
                                    parity1.max(parity2)
                                };
                                stereo_0d.push(InchiStereo0D {
                                    neighbor: [
                                        (i32::from(at[i].sn_orig_at_num[m1]) - 1) as i16,
                                        i as i16,
                                        cur as i16,
                                        (i32::from(at[cur].sn_orig_at_num[m2]) - 1) as i16,
                                    ],
                                    central_atom: if len == 3 { chain[1] as i16 } else { NO_ATOM },
                                    type_: if len == 3 {
                                        INCHI_STEREO_TYPE_ALLENE
                                    } else {
                                        INCHI_STEREO_TYPE_DOUBLE_BOND
                                    },
                                    parity: parity as i8,
                                });
                            } else {
                                ret |= 1;
                            }
                            break;
                        }
                        m2 += 1;
                    }
                }
                m1 += 1;
            }
        }
    }
    (ret, atom, stereo_0d)
}

// ---------------------------------------------------------------------------
// GetINCHI / GetStdINCHI

use super::ichi_io::{
    INCHI_IOS_TYPE_STRING, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT, InchiIosString,
    InchiIostream, IosFile, inchi_ios_eprint, inchi_ios_init, inchi_strbuf_init,
};
use super::ichiparm::{
    FLAG_SET_INP_AT_CHIRAL, FLAG_SET_INP_AT_NONCHIRAL, INCHI_MAX_NUM_ARG, INCHI_OUT_SAVEOPT,
    INCHI_OUT_STDINCHI, InputParms, REQ_MODE_BASIC, REQ_MODE_CHIR_FLG_STEREO,
    REQ_MODE_DIFF_UU_STEREO, REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO,
    REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU, REQ_MODE_STEREO, TG_FLAG_1_5_TAUT,
    TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_RECONNECT_COORD, parse_options_string, print_input_parms,
    read_command_line_parms,
};
use super::ichister::{AB_PARITY_UNDF, AB_PARITY_UNKN};
use super::inchi_dll_b::{InchiOutput, copy_corrected_log_tail, produce_generation_output};
use super::readinch::extract_0d_parities;
use super::runichi::process_one_structure_ex;
use super::runichi2::{LOG_MASK_NO_WARN, treat_errors_in_read_the_structure};
use super::runichi3::{IS_EOF, IS_ERROR, IS_FATAL, IS_OKAY, IS_SKIP, IS_WARNING, StructData};
use super::strutil::OrigAtomData;

/// `inchi_Input`
#[derive(Debug, Clone, Default)]
pub struct InchiInput {
    pub atom: Vec<InchiAtom>,
    pub stereo_0d: Vec<InchiStereo0D>,
    pub options: Option<String>,
}

/// `ExtractOneStructure` without polymer and V3000 data
pub fn extract_one_structure(
    sd: &mut StructData,
    ip: &mut InputParms,
    inp: &InchiInput,
    log: &mut InchiIostream,
    orig: &mut OrigAtomData,
    num_inp: &mut i64,
) -> i32 {
    let v_ab_parity_unknown = if ip.n_mode & REQ_MODE_DIFF_UU_STEREO != 0 {
        AB_PARITY_UNKN
    } else {
        AB_PARITY_UNDF
    };
    *orig = OrigAtomData::default();
    let n = inp.atom.len();
    let max_num_at = if ip.b_large_molecules != 0 {
        32766
    } else {
        1024
    };
    'body: {
        if n == 0 {
            add_error_message(&mut sd.str_err_struct, "Empty structure");
            sd.n_struct_read_error = 98;
            break 'body;
        }
        if n >= max_num_at {
            add_error_message(
                &mut sd.str_err_struct,
                "Too many atoms [did you forget 'LargeMolecules' switch?]",
            );
            sd.n_struct_read_error = 70;
            break 'body;
        }
        let (mut at, n_dim, n_bonds, sz_coord) = extract_atoms_and_bonds(
            &inp.atom,
            ip.b_do_not_add_h != 0,
            &mut sd.str_err_struct,
            &mut sd.n_struct_read_error,
        );
        if sd.n_struct_read_error != 0 {
            break 'body;
        }
        extract_0d_parities(
            &mut at,
            &inp.stereo_0d,
            &mut sd.str_err_struct,
            &mut sd.n_struct_read_error,
            v_ab_parity_unknown,
        );
        if sd.n_struct_read_error != 0 {
            break 'body;
        }
        orig.at = at;
        orig.num_dimensions = n_dim;
        orig.num_inp_atoms = n;
        orig.num_inp_bonds = n_bonds;
        orig.sz_coord = Some(sz_coord);

        // chiral flag
        let (chiral, nonchiral) = (
            super::mol2atom::FLAG_INP_AT_CHIRAL as i32,
            super::mol2atom::FLAG_INP_AT_NONCHIRAL as i32,
        );
        if ip.n_mode & REQ_MODE_CHIR_FLG_STEREO != 0 && ip.n_mode & REQ_MODE_STEREO != 0 {
            if ip.b_chiral_flag & FLAG_SET_INP_AT_CHIRAL != 0 {
                ip.n_mode &= !(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
                sd.b_chiral_flag &= !nonchiral;
                sd.b_chiral_flag |= chiral;
            } else {
                ip.n_mode &= !REQ_MODE_RACEMIC_STEREO;
                ip.n_mode |= REQ_MODE_RELATIVE_STEREO;
                sd.b_chiral_flag &= !chiral;
                sd.b_chiral_flag |= nonchiral;
            }
        } else if ip.b_chiral_flag & FLAG_SET_INP_AT_CHIRAL != 0 {
            sd.b_chiral_flag &= !nonchiral;
            sd.b_chiral_flag |= chiral;
        } else if ip.b_chiral_flag & FLAG_SET_INP_AT_NONCHIRAL != 0 {
            sd.b_chiral_flag &= !chiral;
            sd.b_chiral_flag |= nonchiral;
        }
        *num_inp += 1;
    }
    treat_errors_in_read_the_structure(sd, ip, LOG_MASK_NO_WARN, log, orig, *num_inp)
}

/// `GetINCHI1`: `enforce_std_format` is the `GetStdINCHI` mode. Structures
/// with polymer data and the usage text (given for an empty input without
/// options; it contains the build date) are not supported: the usage case
/// returns `inchi_Ret_EOF` with no log.
fn get_inchi1(inp: &InchiInput, enforce_std_format: bool) -> (i32, InchiOutput) {
    let mut out = InchiOutput::default();
    let mut sd = StructData::default();
    let mut ip = InputParms::default();
    let mut orig = OrigAtomData::default();
    let mut prep = [OrigAtomData::default(), OrigAtomData::default()];
    let mut out_file = InchiIostream::default();
    let mut log = InchiIostream::default();
    let mut strbuf = InchiIosString::default();
    let mut num_inp: i64 = 0;
    inchi_ios_init(&mut out_file, INCHI_IOS_TYPE_STRING, IosFile::None);
    inchi_ios_init(&mut log, INCHI_IOS_TYPE_STRING, IosFile::None);

    let mut n_ret;
    let mut translate_only = false;
    'main: {
        let argv: Vec<Vec<u8>> = match &inp.options {
            Some(o) => parse_options_string(o.as_bytes(), INCHI_MAX_NUM_ARG),
            None => vec![Vec::new()],
        };
        let help = (argv.len() == 1 && inp.atom.is_empty())
            || (argv.len() == 2
                && argv[1].first() == Some(&b'-')
                && (&argv[1][1..] == b"?" || argv[1][1..].eq_ignore_ascii_case(b"help")));
        if help {
            n_ret = IS_EOF;
            translate_only = true;
            break 'main;
        }
        let argv_str: Vec<String> = argv
            .iter()
            .map(|a| String::from_utf8_lossy(a).into_owned())
            .collect();
        let argv_ref: Vec<&str> = argv_str.iter().map(String::as_str).collect();
        let n_ret1 = read_command_line_parms(&argv_ref, &mut ip, true, &mut log);
        ip.b_no_struct_labels = 1;
        if n_ret1 < 0 {
            n_ret = IS_FATAL;
            break 'main;
        }
        if enforce_std_format {
            ip.b_inchi_output_options &= !INCHI_OUT_SAVEOPT;
            ip.b_taut_flags &= !TG_FLAG_RECONNECT_COORD;
            ip.n_mode &= !REQ_MODE_BASIC;
            if ip.n_mode & REQ_MODE_RELATIVE_STEREO != 0 {
                ip.n_mode &= !(REQ_MODE_RACEMIC_STEREO
                    | REQ_MODE_RELATIVE_STEREO
                    | REQ_MODE_CHIR_FLG_STEREO);
            }
            if ip.n_mode & REQ_MODE_RACEMIC_STEREO != 0 {
                ip.n_mode &= !(REQ_MODE_RACEMIC_STEREO
                    | REQ_MODE_RELATIVE_STEREO
                    | REQ_MODE_CHIR_FLG_STEREO);
            }
            if ip.n_mode & REQ_MODE_CHIR_FLG_STEREO != 0 {
                ip.n_mode &= !(REQ_MODE_RACEMIC_STEREO
                    | REQ_MODE_RELATIVE_STEREO
                    | REQ_MODE_CHIR_FLG_STEREO);
            }
            ip.n_mode &= !REQ_MODE_DIFF_UU_STEREO;
            if ip.n_mode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) == 0 {
                ip.n_mode |= REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU;
            }
            ip.b_taut_flags &= !(TG_FLAG_KETO_ENOL_TAUT | TG_FLAG_1_5_TAUT);
            ip.b_inchi_output_options |= INCHI_OUT_STDINCHI;
            ip.b_inchi_output_options &= !INCHI_OUT_SAVEOPT;
        }
        print_input_parms(&mut log, &ip);
        inchi_strbuf_init(
            &mut strbuf,
            INCHI_STRBUF_INITIAL_SIZE,
            INCHI_STRBUF_SIZE_INCREMENT,
        );
        if ip.last_struct_number != 0 && num_inp >= ip.last_struct_number {
            n_ret = IS_EOF;
            break 'main;
        }
        n_ret = extract_one_structure(&mut sd, &mut ip, inp, &mut log, &mut orig, &mut num_inp);
        if n_ret == IS_FATAL || n_ret == IS_EOF || n_ret == IS_ERROR {
            break 'main;
        }
        let n_ret1 = process_one_structure_ex(
            &mut sd,
            &ip,
            &mut log,
            &mut out_file,
            &mut strbuf,
            &mut orig,
            &mut prep,
            num_inp,
        );
        n_ret = n_ret.max(n_ret1);
    }
    if !translate_only {
        produce_generation_output(&mut out, &sd, &ip, &mut log, &mut out_file);
    } else {
        copy_corrected_log_tail(&mut out, &mut log);
    }
    let code = match n_ret {
        IS_SKIP => -2,
        IS_EOF => -1,
        IS_OKAY => 0,
        IS_WARNING => 1,
        IS_ERROR => 2,
        IS_FATAL => 3,
        _ => 4,
    };
    let _ = inchi_ios_eprint;
    (code, out)
}

/// `input_erroneously_contains_pseudoatoms` (which, as in C, looks at the
/// first atom only)
fn pseudoatoms_message(inp: &InchiInput) -> Option<InchiOutput> {
    for _ in 0..inp.atom.len() {
        let e = &inp.atom[0].elname;
        let len = e.iter().position(|&c| c == 0).unwrap_or(e.len());
        if &e[..len] == b"Zz" || &e[..len] == b"*" {
            return Some(InchiOutput {
                sz_message: Some("Unsupported in this mode element '*'".to_string()),
                ..Default::default()
            });
        }
    }
    None
}

/// `GetINCHI`
pub fn get_inchi(inp: &InchiInput) -> (i32, InchiOutput) {
    if let Some(o) = pseudoatoms_message(inp) {
        return (IS_ERROR, o);
    }
    get_inchi1(inp, false)
}

/// `GetStdINCHI`
pub fn get_std_inchi(inp: &InchiInput) -> (i32, InchiOutput) {
    if let Some(o) = pseudoatoms_message(inp) {
        return (IS_ERROR, o);
    }
    get_inchi1(inp, true)
}

#[cfg(test)]
mod fixed_format_tests {
    use super::fixed;

    #[test]
    fn fixed_matches_the_standard_formatter() {
        let mut x = -123.456789f64;
        for i in 0..20000 {
            x = (x * 1.000731 + 0.0371 * f64::from(i % 17)) % 99999.0;
            for p in 0..=4usize {
                assert_eq!(fixed(x, p), format!("{x:.p$}"), "{x} at {p}");
            }
        }
        for x in [
            0.0,
            -0.0,
            0.00005,
            0.00015,
            -0.00005,
            1.00005,
            2.5,
            0.5,
            12345.67895,
            -0.00001,
        ] {
            for p in 0..=4usize {
                assert_eq!(fixed(x, p), format!("{x:.p$}"), "{x} at {p}");
            }
        }
    }
}
