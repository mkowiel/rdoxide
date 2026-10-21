//! Port of `mol2atom.c`: conversion of the data read from a Molfile into the
//! InChI input atom array (`inp_ATOM`), the original-structure bookkeeping,
//! and the v. 1.05 polymer / V3000 extension containers.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust. Licensed under the MIT license
//! (see LICENSES/InChI-MIT.txt).

use super::ichierr::add_error_message;
use super::ichiparm::{FSS_STARS_CYCLED, POLYMERS_MODERN};
use super::ichiprt1::MolCoord;
use super::incomdef::{BOND_TYPE_ALTERN, RADICAL_SINGLET};
use super::inpdef::{InpAtom, MAXVAL};
use super::molfmt::molfile_extract_struc_num;
use super::molfmt::{
    MolFmtData, ReadMolfileExtras, ZERO_ATW_DIFF, molfile_get_xyz_dim_and_norm_factors,
    molfile_has_no_chem_struc, read_molfile, treat_err,
};
use super::strutil::OrigAtomData;
use super::util::{
    EL_NUMBER_H, ERR_ELEM, detect_unusual_el_valence, extract_h_atoms, get_num_h,
    get_periodic_table_number, inchi_stricmp, is_el_a_metal, lrtrim, mystrncpy,
    n_bonds_val_to_metal,
};

pub const MIN_INPUT_BOND_TYPE: i32 = 1;
pub const MAX_INPUT_BOND_TYPE: i32 = 4;
pub const BOND_TYPE_SINGLE: u8 = 1;
const STEREO_SNGL_UP: i8 = 1;
const STEREO_SNGL_EITHER: i8 = 4;
const STEREO_SNGL_DOWN: i8 = 6;
const STEREO_DBLE_EITHER: i8 = 3;
pub const MAX_SDF_VALUE: usize = 255;
pub const MAX_ATOMS: i32 = 32766;
pub const NORMALLY_ALLOWED_INP_MAX_ATOMS: i32 = 1024;
pub const FLAG_INP_AT_CHIRAL: u64 = 1;
pub const FLAG_INP_AT_NONCHIRAL: u64 = 2;
pub const FLAG_SET_INP_LARGE_MOLS: u64 = 16;

/// `OAD_PolymerUnit`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OadPolymerUnit {
    pub id: i32,
    pub type_: i32,
    pub subtype: i32,
    pub conn: i32,
    pub label: i32,
    pub na: i32,
    pub nb: i32,
    pub cyclizable: i32,
    pub cyclized: i32,
    pub xbr1: [f64; 4],
    pub xbr2: [f64; 4],
    pub smt: Vec<u8>,
    pub representation: i32,
    pub cap1: i32,
    pub end_atom1: i32,
    pub end_atom2: i32,
    pub cap2: i32,
    pub cap1_is_undef: i32,
    pub cap2_is_undef: i32,
    pub alist: Vec<i32>,
    pub blist: Vec<i32>,
    pub maxbkbonds: i32,
    pub nbkbonds: i32,
    pub bkbonds: Vec<Vec<i32>>,
}

/// `OAD_Polymer`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OadPolymer {
    pub units: Vec<OadPolymerUnit>,
    pub n: i32,
    pub n_pzz: i32,
    pub pzz: Option<Vec<i32>>,
    pub really_do_frame_shift: i32,
    pub frame_shift_scheme: i32,
    pub treat: i32,
    pub representation: i32,
    pub is_in_reconn: i32,
    pub edit_repeats: i32,
}

/// `OAD_V3000`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OadV3000 {
    pub n_non_star_atoms: i32,
    pub n_star_atoms: i32,
    pub atom_index_orig: Vec<i32>,
    pub atom_index_fin: Vec<i32>,
    pub n_sgroups: i32,
    pub n_3d_constraints: i32,
    pub n_collections: i32,
    pub n_non_haptic_bonds: i32,
    pub n_haptic_bonds: i32,
    pub lists_haptic_bonds: Vec<Vec<i32>>,
    pub n_steabs: i32,
    pub lists_steabs: Vec<Vec<i32>>,
    pub n_sterel: i32,
    pub lists_sterel: Vec<Vec<i32>>,
    pub n_sterac: i32,
    pub lists_sterac: Vec<Vec<i32>>,
}

fn c_len(s: &[u8]) -> usize {
    s.iter().position(|&c| c == 0).unwrap_or(s.len())
}

fn is_in_the_list(list: &[u16], next: u16, len: usize) -> Option<usize> {
    list[..len.min(list.len())].iter().position(|&a| a == next)
}

/// `SetExtOrigAtDataByMolfileExtInput`
pub fn set_ext_orig_at_data_by_molfile_ext_input(
    mfdata: &MolFmtData,
    p_str_err: &mut String,
) -> (i32, Option<OadPolymer>, Option<OadV3000>) {
    let _ = p_str_err;
    let mut err = 0;
    let mut polymer = None;
    let nsgroups = mfdata.ctab.sgroups.len();
    if nsgroups > 0 {
        let mut pol = OadPolymer {
            n: nsgroups as i32,
            is_in_reconn: 0,
            pzz: None,
            edit_repeats: 0,
            really_do_frame_shift: 0,
            frame_shift_scheme: FSS_STARS_CYCLED,
            treat: POLYMERS_MODERN,
            ..Default::default()
        };
        for g in &mfdata.ctab.sgroups {
            let mut u = OadPolymerUnit {
                id: g.id,
                type_: g.type_,
                subtype: g.subtype,
                conn: g.conn,
                label: g.label,
                xbr1: g.xbr1,
                xbr2: g.xbr2,
                smt: g.smt.clone(),
                na: g.alist.len() as i32,
                alist: g.alist.clone(),
                nb: g.blist.len() as i32,
                ..Default::default()
            };
            if u.nb > 0 {
                u.blist = vec![0; 2 * u.nb as usize];
                for (m, &ib) in g.blist.iter().enumerate() {
                    if ib < 1 || ib > mfdata.ctab.n_bonds {
                        treat_err(
                            &mut err,
                            9004,
                            p_str_err,
                            "Polymer unit in Molfile refers to invalid bond",
                        );
                        return (err, None, None);
                    }
                    let b = &mfdata.ctab.bonds.as_ref().unwrap()[(ib - 1) as usize];
                    let (ia1, ia2) = (i32::from(b.atnum1), i32::from(b.atnum2));
                    u.blist[2 * m] = ia1;
                    u.blist[2 * m + 1] = ia2;
                    let atoms = mfdata.ctab.atoms.as_ref().unwrap();
                    if atoms[(ia1 - 1) as usize].symbol == b"H"
                        || atoms[(ia2 - 1) as usize].symbol == b"H"
                    {
                        treat_err(
                            &mut err,
                            9002,
                            p_str_err,
                            "Hydrogen as polymer end group is not supported",
                        );
                        return (err, None, None);
                    }
                }
            }
            pol.units.push(u);
        }
        polymer = Some(pol);
    }
    let mut v3000 = None;
    if let Some(mpv) = mfdata.ctab.v3000.as_ref() {
        let mut pv = OadV3000 {
            n_collections: mpv.n_collections,
            n_haptic_bonds: mpv.n_haptic_bonds,
            n_non_haptic_bonds: mpv.n_non_haptic_bonds,
            n_sgroups: mpv.n_sgroups,
            n_non_star_atoms: mpv.n_non_star_atoms,
            n_star_atoms: mpv.n_star_atoms,
            n_steabs: mpv.n_steabs,
            n_sterac: mpv.n_sterac,
            n_sterel: mpv.n_sterel,
            n_3d_constraints: mpv.n_3d_constraints,
            ..Default::default()
        };
        // the C code copies `n_atoms` bytes (not ints) into the zeroed int arrays
        let copy_bytes = |src: &[i32]| {
            let n = mfdata.ctab.n_atoms.max(0) as usize;
            let mut out = vec![0i32; n];
            let bytes: Vec<u8> = src.iter().flat_map(|v| v.to_le_bytes()).collect();
            let mut raw: Vec<u8> = out.iter().flat_map(|v| v.to_le_bytes()).collect();
            let k = n.min(bytes.len()).min(raw.len());
            raw[..k].copy_from_slice(&bytes[..k]);
            for (o, ch) in out.iter_mut().zip(raw.chunks(4)) {
                *o = i32::from_le_bytes([ch[0], ch[1], ch[2], ch[3]]);
            }
            out
        };
        if !mpv.atom_index_orig.is_empty() {
            pv.atom_index_orig = copy_bytes(&mpv.atom_index_orig);
        }
        if !mpv.atom_index_fin.is_empty() {
            pv.atom_index_fin = copy_bytes(&mpv.atom_index_fin);
        }
        let copy_lists = |lists: &[Vec<i32>], n: i32, first_len_at: usize, add: i32| {
            lists
                .iter()
                .take(n.max(0) as usize)
                .map(|l| {
                    let nn = (l[first_len_at] + add).max(0) as usize;
                    l[..nn].to_vec()
                })
                .collect::<Vec<_>>()
        };
        if mpv.n_haptic_bonds != 0 && !mpv.haptic_bonds.is_empty() {
            pv.lists_haptic_bonds = copy_lists(&mpv.haptic_bonds, mpv.n_haptic_bonds, 2, 3);
        }
        if mpv.n_steabs != 0 && !mpv.steabs.is_empty() {
            pv.lists_steabs = copy_lists(&mpv.steabs, mpv.n_steabs, 1, 2);
        }
        if mpv.n_sterac != 0 && !mpv.sterac.is_empty() {
            pv.lists_sterac = copy_lists(&mpv.sterac, mpv.n_sterac, 1, 2);
        }
        if mpv.n_sterel != 0 && !mpv.sterel.is_empty() {
            pv.lists_sterel = copy_lists(&mpv.sterel, mpv.n_sterel, 1, 2);
        }
        v3000 = Some(pv);
    }
    (err, polymer, v3000)
}

fn elname_str(a: &InpAtom) -> String {
    String::from_utf8_lossy(&a.elname[..c_len(&a.elname)]).into_owned()
}

fn set_elname(a: &mut InpAtom, name: &[u8]) {
    a.elname = [0; super::inpdef::ATOM_EL_LEN];
    let n = name.len().min(a.elname.len() - 1);
    a.elname[..n].copy_from_slice(&name[..n]);
}

/// `MakeInpAtomsFromMolfileData`: the atoms (`None` for a structure without
/// atoms), the number of atoms and bonds, and the error code in `err`
pub fn make_inp_atoms_from_molfile_data(
    mfdata: &MolFmtData,
    at_inp: Option<Vec<InpAtom>>,
    num_atoms: &mut i32,
    num_bonds: &mut i32,
    err: &mut i32,
    p_str_err: &mut String,
    do_not_add_h: bool,
) -> Option<Vec<InpAtom>> {
    *err = 0;
    *num_atoms = mfdata.ctab.n_atoms;
    *num_bonds = 0;
    if molfile_has_no_chem_struc(Some(mfdata)) {
        return None;
    }
    let n = *num_atoms as usize;
    let mut at = at_inp.unwrap_or_default();
    if at.len() < n {
        at.resize(n, InpAtom::default());
    }
    let atoms = mfdata.ctab.atoms.as_ref().unwrap();
    for i in 0..n {
        let src = &atoms[i];
        let a = &mut at[i];
        let mut sbuf = src.symbol.clone();
        sbuf.push(0);
        let mut tmp = vec![0u8; a.elname.len()];
        mystrncpy(&mut tmp, &sbuf, a.elname.len());
        a.elname.copy_from_slice(&tmp);
        a.orig_at_number = (i + 1) as u16;
        a.charge = src.charge;
        a.radical = src.radical;
        let iso = i32::from(src.mass_difference);
        a.iso_atw_diff = (if iso == ZERO_ATW_DIFF {
            1
        } else if iso > 0 {
            iso + 1
        } else {
            iso
        }) as i8;
        // SINGLET_IS_TRIPLET
        if i32::from(a.radical) == RADICAL_SINGLET {
            a.radical = 3;
        }
        let mut n1 = get_periodic_table_number(&elname_str(a));
        if n1 == ERR_ELEM {
            // the element name contains more than one element: extract the H atoms
            let mut name = a.elname[..c_len(&a.elname)].to_vec();
            a.num_h = extract_h_atoms(&mut name, &mut a.num_iso_h) as i8;
            set_elname(a, &name);
            if a.elname[0] == 0 && a.numh() != 0 {
                set_elname(a, b"H");
                if a.num_iso_h_total() != 0 {
                    for j in (0..a.num_iso_h.len()).rev() {
                        if a.num_iso_h[j] != 0 {
                            a.num_iso_h[j] -= 1;
                            a.iso_atw_diff = (1 + j) as i8;
                            break;
                        }
                    }
                } else {
                    a.num_h -= 1;
                }
            }
            n1 = get_periodic_table_number(&elname_str(a));
            if n1 == ERR_ELEM {
                n1 = 0;
            }
        }
        a.el_number = n1 as u8;
        if n1 == 0 {
            *err = -2;
            treat_err(err, -2, p_str_err, "Unknown element(s):");
            treat_err(err, -2, p_str_err, &elname_str(a));
            return Some(at);
        }
        if n1 == i32::from(EL_NUMBER_H) && a.iso_atw_diff == 0 {
            match a.elname[0] {
                b'D' => {
                    a.iso_atw_diff = 2;
                    set_elname(a, b"H");
                }
                b'T' => {
                    a.iso_atw_diff = 3;
                    set_elname(a, b"H");
                }
                _ => {}
            }
        }
    }
    let mut bonds = 0;
    if let Some(blist) = mfdata.ctab.bonds.as_ref() {
        for b in blist.iter().take(mfdata.ctab.n_bonds.max(0) as usize) {
            let bond_stereo = b.bond_stereo;
            let mut bond_type = b.bond_type;
            let a1 = i32::from(b.atnum1) - 1;
            let a2 = i32::from(b.atnum2) - 1;
            if a1 < 0 || a1 >= *num_atoms || a2 < 0 || a2 >= *num_atoms || a1 == a2 {
                *err |= 1;
                treat_err(err, 0, p_str_err, "Bond to nonexistent atom");
                continue;
            }
            let (a1, a2) = (a1 as usize, a2 as usize);
            let p1 = is_in_the_list(&at[a1].neighbor, a2 as u16, at[a1].valence as usize);
            let p2 = is_in_the_list(&at[a2].neighbor, a1 as u16, at[a2].valence as usize);
            let maxval = MAXVAL as i8;
            let (n1, n2);
            if (p1.is_some() || p2.is_some())
                && (p1.is_some() || at[a1].valence < maxval)
                && (p2.is_some() || at[a2].valence < maxval)
            {
                n1 = p1.unwrap_or_else(|| {
                    at[a1].valence += 1;
                    at[a1].valence as usize - 1
                });
                n2 = p2.unwrap_or_else(|| {
                    at[a2].valence += 1;
                    at[a2].valence as usize - 1
                });
                treat_err(err, 0, p_str_err, "Multiple bonds between two atoms");
                *err |= 2;
            } else if p1.is_none()
                && p2.is_none()
                && at[a1].valence < maxval
                && at[a2].valence < maxval
            {
                n1 = at[a1].valence as usize;
                at[a1].valence += 1;
                n2 = at[a2].valence as usize;
                at[a2].valence += 1;
                bonds += 1;
            } else {
                *err |= 4;
                let who = if at[a1].valence >= maxval {
                    elname_str(&at[a1])
                } else {
                    elname_str(&at[a2])
                };
                treat_err(
                    err,
                    0,
                    p_str_err,
                    &format!("Atom '{who}' has more than {MAXVAL} bonds"),
                );
                continue;
            }
            if i32::from(bond_type) < MIN_INPUT_BOND_TYPE
                || i32::from(bond_type) > MAX_INPUT_BOND_TYPE
            {
                let msg = format!("{bond_type}");
                bond_type = 1;
                treat_err(err, 0, p_str_err, "Unrecognized bond type:");
                treat_err(err, 0, p_str_err, &msg);
                *err |= 8;
            }
            if n1 < MAXVAL && n2 < MAXVAL {
                at[a1].bond_type[n1] = bond_type as u8;
                at[a2].bond_type[n2] = bond_type as u8;
                at[a1].neighbor[n1] = a2 as u16;
                at[a2].neighbor[n2] = a1 as u16;
                if bond_stereo == 3 {
                    at[a1].bond_stereo[n1] = STEREO_DBLE_EITHER;
                    at[a2].bond_stereo[n2] = STEREO_DBLE_EITHER;
                } else if matches!(bond_stereo, 1 | 4 | 6) {
                    let c = match bond_stereo {
                        1 => STEREO_SNGL_UP,
                        4 => STEREO_SNGL_EITHER,
                        _ => STEREO_SNGL_DOWN,
                    };
                    at[a1].bond_stereo[n1] = c;
                    at[a2].bond_stereo[n2] = -c;
                } else if bond_stereo != 0 {
                    *err |= 16;
                    treat_err(err, 0, p_str_err, "Unrecognized bond stereo");
                    continue;
                }
            }
        }
    }
    *num_bonds = bonds;
    calculate_valences(
        Some(mfdata),
        &mut at,
        *num_atoms,
        do_not_add_h,
        err,
        p_str_err,
    );
    Some(at)
}

/// `calculate_valences`
pub fn calculate_valences(
    mfdata: Option<&MolFmtData>,
    at: &mut [InpAtom],
    num_atoms: i32,
    do_not_add_h: bool,
    err: &mut i32,
    p_str_err: &mut String,
) {
    for b_non_metal in 0..2 {
        for a1 in 0..num_atoms.max(0) as usize {
            let mut num_bond_type =
                [0i32; (MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1) as usize];
            if (b_non_metal == 1) == is_el_a_metal(i32::from(at[a1].el_number)) {
                continue;
            }
            let valence: i32 = match mfdata.and_then(|m| m.ctab.atoms.as_ref()) {
                Some(atoms) => i32::from(atoms[a1].valence),
                None => i32::from(at[a1].chem_bonds_valence),
            };
            at[a1].chem_bonds_valence = 0;
            let has_metal_neighbor = false;
            for n1 in 0..at[a1].valence as usize {
                let mut bond_type = i32::from(at[a1].bond_type[n1]) - MIN_INPUT_BOND_TYPE;
                if !(0..=MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE).contains(&bond_type) {
                    bond_type = 0;
                    treat_err(
                        err,
                        0,
                        p_str_err,
                        "Unknown bond type in MOLfile assigned as a single bond",
                    );
                }
                num_bond_type[bond_type as usize] += 1;
            }
            let mut n1 = 0;
            // MAX_INPUT_BOND_TYPE (4) > 3: only the single, double and triple bonds are added
            while MIN_INPUT_BOND_TYPE + n1 <= 3 {
                at[a1].chem_bonds_valence = at[a1]
                    .chem_bonds_valence
                    .wrapping_add(((MIN_INPUT_BOND_TYPE + n1) * num_bond_type[n1 as usize]) as i8);
                n1 += 1;
            }
            let mut n2 = 0;
            if MIN_INPUT_BOND_TYPE <= BOND_TYPE_ALTERN && BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE {
                n2 = num_bond_type[(BOND_TYPE_ALTERN - MIN_INPUT_BOND_TYPE) as usize];
                if n2 != 0 {
                    match n2 {
                        2 => at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_add(3),
                        3 => at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_add(4),
                        _ => {
                            for n1 in 0..at[a1].valence as usize {
                                if i32::from(at[a1].bond_type[n1]) == BOND_TYPE_ALTERN {
                                    let a2 = at[a1].neighbor[n1] as usize;
                                    match is_in_the_list(
                                        &at[a2].neighbor,
                                        a1 as u16,
                                        at[a2].valence as usize,
                                    ) {
                                        Some(p) => {
                                            at[a1].bond_type[n1] = BOND_TYPE_SINGLE;
                                            at[a2].bond_type[p] = BOND_TYPE_SINGLE;
                                        }
                                        None => {
                                            *err = -2;
                                            treat_err(
                                                err,
                                                0,
                                                p_str_err,
                                                "Program error interpreting MOLfile",
                                            );
                                            return;
                                        }
                                    }
                                }
                            }
                            at[a1].chem_bonds_valence =
                                at[a1].chem_bonds_valence.wrapping_add(n2 as i8);
                            *err |= 32;
                            treat_err(
                                err,
                                0,
                                p_str_err,
                                "Atom has 1 or more than 3 aromatic bonds",
                            );
                            n2 = 0;
                        }
                    }
                }
            }
            if n2 != 0 && valence == 0 {
                let num_h = at[a1].numh();
                let chem_valence = i32::from(at[a1].chem_bonds_valence);
                let arom = detect_unusual_el_valence(
                    i32::from(at[a1].el_number),
                    i32::from(at[a1].charge),
                    i32::from(at[a1].radical),
                    chem_valence,
                    num_h,
                    i32::from(at[a1].valence),
                );
                let no_arom = detect_unusual_el_valence(
                    i32::from(at[a1].el_number),
                    i32::from(at[a1].charge),
                    i32::from(at[a1].radical),
                    chem_valence - 1,
                    num_h,
                    i32::from(at[a1].valence),
                );
                if arom != 0 && no_arom == 0 && n_bonds_val_to_metal(at, a1) == 0 {
                    at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_sub(1);
                }
            } else if n2 != 0 && valence != 0 {
                let num_h = at[a1].numh();
                let chem_valence = i32::from(at[a1].chem_bonds_valence) + num_h;
                if valence == chem_valence - 1 {
                    at[a1].chem_bonds_valence = at[a1].chem_bonds_valence.wrapping_sub(1);
                }
            }
            if let Some(m) = mfdata {
                let src = &m.ctab.atoms.as_ref().unwrap()[a1];
                let a = &at[a1];
                at[a1].num_h = get_num_h(
                    &elname_str(a),
                    i32::from(a.num_h),
                    Some(&a.num_iso_h),
                    i32::from(a.charge),
                    i32::from(a.radical),
                    i32::from(a.chem_bonds_valence),
                    i32::from(src.valence),
                    src.atom_aliased_flag != 0,
                    do_not_add_h,
                    has_metal_neighbor,
                ) as i8;
            }
        }
    }
}

/// `SetInpAtomsXYZ`: copies the coordinates, returns the number of dimensions
pub fn set_inp_atoms_xyz(
    mfdata: &MolFmtData,
    num_atoms: i32,
    at: &mut [InpAtom],
    err: &mut i32,
    p_str_err: &mut String,
) -> i32 {
    // NORMALIZE_INP_COORD == 0: the coordinates are kept unchanged
    let nx = molfile_get_xyz_dim_and_norm_factors(mfdata, false, err, p_str_err);
    if nx.num_dimensions == 0 {
        return 0;
    }
    let atoms = mfdata.ctab.atoms.as_ref().unwrap();
    for i in 0..num_atoms.max(0) as usize {
        at[i].x = atoms[i].fx;
        at[i].y = atoms[i].fy;
        at[i].z = atoms[i].fz;
    }
    nx.num_dimensions
}

/// the result of [`read_molfile_to_inp_atoms`]
#[derive(Debug, Default)]
pub struct MolfileAtoms {
    pub num_atoms: i32,
    pub at: Option<Vec<InpAtom>>,
    pub sz_coord: Option<Vec<MolCoord>>,
    pub polymer: Option<OadPolymer>,
    pub v3000: Option<OadV3000>,
    pub num_dimensions: i32,
    pub num_bonds: i32,
    /// the number of a "Structure #n" name line, when a header was available
    pub molfile_number: Option<i64>,
    pub err: i32,
}

/// the optional in/out arguments of `ReadMolfileToInpAtoms`
#[derive(Debug, Default)]
pub struct MolfileAtomsArgs<'a> {
    pub sdf_label: Option<&'a [u8]>,
    pub sdf_value: Option<&'a mut Vec<u8>>,
    pub id: Option<&'a mut u64>,
    pub inp_atom_flags: Option<&'a mut u64>,
}

/// `ReadMolfileToInpAtoms` (`want_atoms` = false is the mode that only counts
/// the atoms: `at == NULL`)
#[allow(clippy::too_many_arguments)]
pub fn read_molfile_to_inp_atoms(
    inp: &mut super::ichi_io::InchiIostream,
    do_not_add_h: bool,
    want_atoms: bool,
    get_orig_coord: bool,
    treat_polymers: i32,
    treat_np_zz: i32,
    max_num_at: i32,
    args: MolfileAtomsArgs,
    p_str_err: &mut String,
    no_warnings: bool,
) -> MolfileAtoms {
    let MolfileAtomsArgs {
        sdf_label,
        mut sdf_value,
        id,
        mut inp_atom_flags,
    } = args;
    let mut res = MolfileAtoms::default();
    let saved_value = sdf_value.as_deref().cloned();
    if let Some(v) = sdf_value.as_deref_mut() {
        v.clear();
    }
    let mut only = MolFmtData::default();
    let (mfdata, err) = read_molfile(
        inp,
        if want_atoms { None } else { Some(&mut only) },
        get_orig_coord,
        treat_polymers,
        treat_np_zz,
        ReadMolfileExtras {
            name: None,
            id,
            sdf_label,
            sdf_value: sdf_value.as_deref_mut(),
        },
        p_str_err,
        no_warnings,
    );
    res.err = err;
    // `pHdr`: the header of the data, or of the header-only block when the read failed
    let hdr = match (&mfdata, want_atoms) {
        (Some(m), true) => Some(&m.hdr),
        (None, false) => Some(&only.hdr),
        _ => None,
    };
    if let Some(h) = hdr {
        res.molfile_number = Some(molfile_extract_struc_num(h));
    }
    if let (Some(v), Some(label), Some(h)) = (&mut sdf_value, sdf_label, hdr)
        && v.is_empty()
        && !label.is_empty()
    {
        let label = String::from_utf8_lossy(label).into_owned();
        let take = |s: &[u8]| {
            let mut t = s[..c_len(s)].to_vec();
            t.truncate(MAX_SDF_VALUE);
            lrtrim(&mut t);
            t
        };
        if inchi_stricmp(&label, "MolfileName") == 0 {
            **v = take(&h.molname);
        } else if inchi_stricmp(&label, "MolfileLine2") == 0 {
            **v = take(&h.line2);
        } else if inchi_stricmp(&label, "MolfileComment") == 0 {
            **v = take(&h.comment);
        } else if inchi_stricmp(&label, "MolfileIntRegNo") == 0 && h.internal_regno != 0 {
            **v = h.internal_regno.to_string().into_bytes();
        }
        if v.is_empty()
            && let Some(old) = &saved_value
        {
            **v = old.clone();
        }
    }
    if let (Some(mut m), true, 0) = (mfdata.clone(), want_atoms, res.err) {
        let mut num_atoms = 0;
        if m.ctab.n_atoms <= max_num_at {
            let mut err = 0;
            let mut at = make_inp_atoms_from_molfile_data(
                &m,
                None,
                &mut num_atoms,
                &mut res.num_bonds,
                &mut err,
                p_str_err,
                do_not_add_h,
            );
            if err >= 0 {
                res.num_dimensions = match at.as_mut() {
                    Some(a) => set_inp_atoms_xyz(&m, num_atoms, a, &mut err, p_str_err),
                    None => set_inp_atoms_xyz(&m, num_atoms, &mut [], &mut err, p_str_err),
                };
                if get_orig_coord {
                    res.sz_coord = m.ctab.coords.take();
                }
                if let Some(a) = at.as_ref() {
                    for x in a.iter().take(num_atoms.max(0) as usize) {
                        if i32::from(x.num_iso_h[0])
                            + i32::from(x.num_iso_h[1])
                            + i32::from(x.num_iso_h[2])
                            > MAXVAL as i32
                        {
                            err = 70 + 8;
                            treat_err(&mut err, 0, p_str_err, "Too many hydrogens at heavy atom");
                            num_atoms = -1;
                            break;
                        }
                    }
                }
                if err == 0 {
                    let (e, pol, v3) = set_ext_orig_at_data_by_molfile_ext_input(&m, p_str_err);
                    res.polymer = pol;
                    res.v3000 = v3;
                    if e != 0 {
                        err = 80;
                        num_atoms = -1;
                    }
                }
            } else if err == -2 {
                err = 90;
                num_atoms = -1;
            }
            res.at = at;
            res.err = err;
        } else {
            let mut err = 0;
            treat_err(
                &mut err,
                0,
                p_str_err,
                "Too many atoms [did you forget 'LargeMolecules' switch?]",
            );
            res.err = 70;
            num_atoms = -1;
        }
        if res.err > 0 {
            res.err += 100;
        }
        if num_atoms > 0
            && res.at.is_some()
            && let Some(f) = inp_atom_flags.as_mut()
        {
            **f |= if m.ctab.chiral_flag != 0 {
                FLAG_INP_AT_CHIRAL
            } else {
                FLAG_INP_AT_NONCHIRAL
            };
        }
        res.num_atoms = num_atoms;
    } else if !want_atoms {
        res.num_atoms = only.ctab.n_atoms;
    }
    res
}

/// `CreateOrigInpDataFromMolfile`: reads one structure (the legacy merge of
/// several structures is supported as in the C code) into `orig_at_data`;
/// returns the number of atoms, `err` is the error code
#[allow(clippy::too_many_arguments)]
pub fn create_orig_inp_data_from_molfile(
    inp: &mut super::ichi_io::InchiIostream,
    orig_at_data: &mut Option<OrigAtomData>,
    merge_all_input_structures: bool,
    get_orig_coord: bool,
    do_not_add_h: bool,
    treat_polymers: i32,
    treat_np_zz: i32,
    mut args: MolfileAtomsArgs,
    molfile_number: &mut i64,
    err: &mut i32,
    p_str_err: &mut String,
    no_warnings: bool,
) -> i32 {
    p_str_err.clear();
    let flags_now = args.inp_atom_flags.as_deref().copied().unwrap_or(0);
    let max_num_at = if flags_now & FLAG_SET_INP_LARGE_MOLS == 0 {
        NORMALLY_ALLOWED_INP_MAX_ATOMS
    } else {
        MAX_ATOMS
    };
    let mut n_num_atoms = 0;
    loop {
        let r = read_molfile_to_inp_atoms(
            inp,
            do_not_add_h,
            orig_at_data.is_some(),
            get_orig_coord && orig_at_data.is_some(),
            treat_polymers,
            treat_np_zz,
            max_num_at,
            MolfileAtomsArgs {
                sdf_label: args.sdf_label,
                sdf_value: args.sdf_value.as_deref_mut(),
                id: args.id.as_deref_mut(),
                inp_atom_flags: args.inp_atom_flags.as_deref_mut(),
            },
            p_str_err,
            no_warnings,
        );
        *err = r.err;
        if let Some(n) = r.molfile_number {
            *molfile_number = n;
        }
        let num_new = r.num_atoms;
        if num_new <= 0 && *err == 0 {
            treat_err(err, 0, p_str_err, "Empty structure");
            *err = 98;
        } else if orig_at_data.is_some()
            && num_new == 0
            && 10 < *err
            && *err < 20
            && orig_at_data.as_ref().unwrap().num_inp_atoms > 0
            && merge_all_input_structures
        {
            *err = 0;
            break;
        } else if num_new > 0 && orig_at_data.is_some() {
            let od = orig_at_data.as_mut().unwrap();
            let total = num_new + od.num_inp_atoms as i32;
            n_num_atoms = total;
            if total >= max_num_at {
                treat_err(
                    err,
                    0,
                    p_str_err,
                    "Too many atoms [did you forget 'LargeMolecules' switch?]",
                );
                *err = 70;
                od.num_inp_atoms = usize::MAX; // C: -1
            } else if od.at.is_empty() {
                od.at = r.at.unwrap_or_default();
                od.sz_coord = r.sz_coord;
                od.num_inp_atoms = num_new as usize;
                od.num_inp_bonds = r.num_bonds;
                od.num_dimensions = r.num_dimensions;
                od.polymer = r.polymer;
                od.v3000 = r.v3000;
            } else {
                let mut at_new = r.at.unwrap_or_default();
                let old_n = od.num_inp_atoms;
                if old_n > 0 {
                    for a in at_new.iter_mut().take(num_new as usize) {
                        for j in 0..a.valence as usize {
                            a.neighbor[j] += old_n as u16;
                        }
                        a.orig_at_number += old_n as u16;
                    }
                }
                let mut merged: Vec<InpAtom> = od.at[..old_n.min(od.at.len())].to_vec();
                merged.extend(at_new.into_iter().take(num_new as usize));
                od.at = merged;
                if let (Some(old), Some(new)) = (od.sz_coord.as_ref(), r.sz_coord.as_ref()) {
                    let mut c = old[..old_n.min(old.len())].to_vec();
                    c.extend(new.iter().take(num_new as usize).cloned());
                    od.sz_coord = Some(c);
                }
                od.num_inp_atoms += num_new as usize;
                od.num_inp_bonds += r.num_bonds;
                od.num_dimensions = od.num_dimensions.max(r.num_dimensions);
                od.polymer = r.polymer;
                od.v3000 = r.v3000;
            }
        } else if num_new > 0 {
            n_num_atoms += num_new;
        }
        if !(*err == 0 && merge_all_input_structures) {
            break;
        }
    }
    if *err != 0
        && let Some(od) = orig_at_data.as_mut()
    {
        *od = OrigAtomData::default();
    }
    if *err != 0 && !(10 < *err && *err < 20) && p_str_err.is_empty() {
        treat_err(err, 0, p_str_err, "Unknown error");
    }
    match orig_at_data {
        Some(od) => od.num_inp_atoms as i32,
        None => n_num_atoms,
    }
}

#[allow(dead_code)]
fn _unused() {
    let _ = add_error_message;
}
