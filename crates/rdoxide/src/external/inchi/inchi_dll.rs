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
pub type MolCoord = [u8; 30];

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
        format!("{x:.2}")
    } else if x < -99999.999 {
        format!("{x:.3}")
    } else if x < 99999.9999 {
        format!("{x:.4}")
    } else if x < 999999.999 {
        format!("{x:.3}")
    } else if x < 9999999.99 {
        format!("{x:.2}")
    } else if x < 99999999.9 {
        format!("{x:.1}")
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
    let mut sz_coord = vec![[0u8; 30]; n.max(1)];
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
