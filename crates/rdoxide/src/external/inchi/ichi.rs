//! The `INChI` result structures (`INCHI_BASE/src/ichi.h`) and the component
//! comparison functions of `ichimake.c` used to sort the components.
//!
//! NULL C pointers are `None` / empty vectors; the element counts of the C
//! structures are kept as separate fields because the comparison functions
//! read them.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichican2::{TAUT_NON, TAUT_YES};
use super::ichimake::{compare_hill_formulas, compare_hill_formulas_no_h};

pub const INCHI_FLAG_REL_STEREO: u64 = 0x0002;
pub const INCHI_FLAG_RAC_STEREO: u64 = 0x0004;

/// `INChI_IsotopicAtom`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InChIIsotopicAtom {
    pub n_atom_number: u16,
    pub n_iso_difference: i16,
    pub n_num_h: i16,
    pub n_num_d: i16,
    pub n_num_t: i16,
}

/// `INChI_IsotopicTGroup`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InChIIsotopicTGroup {
    pub n_t_group_number: u16,
    pub n_num_h: u16,
    pub n_num_d: u16,
    pub n_num_t: u16,
}

/// `INChI_Stereo`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InChIStereo {
    pub n_number_of_stereo_centers: i32,
    pub n_number: Vec<u16>,
    pub t_parity: Vec<i8>,
    pub n_number_inv: Vec<u16>,
    pub t_parity_inv: Vec<i8>,
    pub n_comp_inv2_abs: i32,
    pub b_trivial_inv: i32,
    pub n_number_of_stereo_bonds: i32,
    pub n_bond_atom1: Vec<u16>,
    pub n_bond_atom2: Vec<u16>,
    pub b_parity: Vec<i8>,
}

/// `INChI`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InChI {
    pub n_error_code: i32,
    pub n_flags: u64,
    pub n_total_charge: i32,
    pub n_number_of_atoms: i32,
    pub sz_hill_formula: Option<String>,
    pub n_atom: Vec<u8>,
    pub len_conn_table: i32,
    pub n_conn_table: Vec<u16>,
    pub len_tautomer: i32,
    pub n_tautomer: Vec<u16>,
    pub n_num_h: Vec<i8>,
    pub n_num_h_fixed: Option<Vec<i8>>,
    pub n_number_of_isotopic_atoms: i32,
    pub isotopic_atom: Vec<InChIIsotopicAtom>,
    pub n_number_of_isotopic_t_groups: i32,
    pub isotopic_t_group: Vec<InChIIsotopicTGroup>,
    pub stereo: Option<InChIStereo>,
    pub stereo_isotopic: Option<InChIStereo>,
    pub n_possible_locations_of_isotopic_h: Vec<u16>,
    pub b_deleted: i32,
}

/// `INCHI_SORT` (the `INChI_Aux` pointers are not needed here)
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InChISort {
    pub p_inchi: [Option<InChI>; 2],
    pub ord_number: i16,
}

fn nz(x: i32) -> Option<i32> {
    (x != 0).then_some(x)
}

/// `CompareTautNonIsoPartOfINChI`
pub fn compare_taut_non_iso_part_of_inchi(i1: &InChI, i2: &InChI) -> i32 {
    let len1 = if i1.len_tautomer > 0 && i1.n_tautomer[0] != 0 {
        i1.len_tautomer
    } else {
        0
    };
    let len2 = if i2.len_tautomer > 0 && i2.n_tautomer[0] != 0 {
        i2.len_tautomer
    } else {
        0
    };
    if let Some(r) = nz(len2 - len1) {
        return r;
    }
    for i in 0..len1 as usize {
        if let Some(r) = nz(i32::from(i2.n_tautomer[i]) - i32::from(i1.n_tautomer[i])) {
            return r;
        }
    }
    0
}

/// `CompareInchiStereo` (`REL_RAC_STEREO_IGN_1_SC` is 0)
pub fn compare_inchi_stereo(
    stereo1: Option<&InChIStereo>,
    n_flags1: u64,
    stereo2: Option<&InChIStereo>,
    n_flags2: u64,
) -> i32 {
    if let (Some(s1), Some(s2)) = (stereo1, stereo2) {
        let num = s2.n_number_of_stereo_bonds.min(s1.n_number_of_stereo_bonds) as usize;
        for i in 0..num {
            if let Some(r) = nz(i32::from(s2.n_bond_atom1[i]) - i32::from(s1.n_bond_atom1[i])) {
                return r;
            }
            if let Some(r) = nz(i32::from(s2.n_bond_atom2[i]) - i32::from(s1.n_bond_atom2[i])) {
                return r;
            }
            if let Some(r) = nz(i32::from(s2.b_parity[i]) - i32::from(s1.b_parity[i])) {
                return r;
            }
        }
        if let Some(r) = nz(s2.n_number_of_stereo_bonds - s1.n_number_of_stereo_bonds) {
            return r;
        }
        let num = s2
            .n_number_of_stereo_centers
            .min(s1.n_number_of_stereo_centers) as usize;
        for i in 0..num {
            if let Some(r) = nz(i32::from(s2.n_number[i]) - i32::from(s1.n_number[i])) {
                return r;
            }
            if let Some(r) = nz(i32::from(s2.t_parity[i]) - i32::from(s1.t_parity[i])) {
                return r;
            }
        }
        if let Some(r) = nz(s2.n_number_of_stereo_centers - s1.n_number_of_stereo_centers) {
            return r;
        }
        if (n_flags1 | n_flags2) & (INCHI_FLAG_RAC_STEREO | INCHI_FLAG_REL_STEREO) == 0 {
            if let Some(r) =
                nz(i32::from(s2.n_comp_inv2_abs < 0) - i32::from(s1.n_comp_inv2_abs < 0))
            {
                return r;
            }
        }
    } else if stereo2
        .is_some_and(|s| s.n_number_of_stereo_bonds > 0 || s.n_number_of_stereo_centers > 0)
    {
        return 1;
    } else if stereo1
        .is_some_and(|s| s.n_number_of_stereo_bonds > 0 || s.n_number_of_stereo_centers > 0)
    {
        return -1;
    }
    0
}

const SP3_NONE: i32 = 0;
const SP3_ONLY: i32 = 1;
const SP3_ABS: i32 = 2;
const SP3_REL: i32 = 4;
const SP3_RAC: i32 = 8;

/// `GetSp3RelRacAbs`: `SP3_*` type of the sp3 stereo layer
pub fn get_sp3_rel_rac_abs(p_inchi: Option<&InChI>, stereo: Option<&InChIStereo>) -> i32 {
    if let (Some(p), Some(s)) = (p_inchi, stereo) {
        if p.b_deleted == 0 && s.n_number_of_stereo_centers > 0 {
            return if s.n_comp_inv2_abs != 0 {
                if p.n_flags & INCHI_FLAG_REL_STEREO != 0 {
                    SP3_REL
                } else if p.n_flags & INCHI_FLAG_RAC_STEREO != 0 {
                    SP3_RAC
                } else {
                    SP3_ABS
                }
            } else {
                SP3_ONLY
            };
        }
    }
    SP3_NONE
}

fn nonempty(p: &Option<InChI>) -> Option<&InChI> {
    p.as_ref()
}

fn compare_isotopic_part(i1: &InChI, i2: &InChI, with_t_groups: bool) -> i32 {
    if let Some(r) = nz(i2.n_number_of_isotopic_atoms - i1.n_number_of_isotopic_atoms) {
        return r;
    }
    let num = i1.n_number_of_isotopic_atoms as usize;
    for i in 0..num {
        let (a, b) = (&i2.isotopic_atom[i], &i1.isotopic_atom[i]);
        if let Some(r) = nz(i32::from(a.n_atom_number) - i32::from(b.n_atom_number)) {
            return r;
        }
        if let Some(r) = nz(i32::from(a.n_iso_difference) - i32::from(b.n_iso_difference)) {
            return r;
        }
    }
    for i in 0..num {
        let (a, b) = (&i2.isotopic_atom[i], &i1.isotopic_atom[i]);
        if let Some(r) = nz(i32::from(a.n_num_t) - i32::from(b.n_num_t)) {
            return r;
        }
        if let Some(r) = nz(i32::from(a.n_num_d) - i32::from(b.n_num_d)) {
            return r;
        }
        if let Some(r) = nz(i32::from(a.n_num_h) - i32::from(b.n_num_h)) {
            return r;
        }
    }
    if with_t_groups {
        if let Some(r) = nz(i2.n_number_of_isotopic_t_groups - i1.n_number_of_isotopic_t_groups) {
            return r;
        }
        for i in 0..i1.n_number_of_isotopic_t_groups as usize {
            let (a, b) = (&i2.isotopic_t_group[i], &i1.isotopic_t_group[i]);
            if let Some(r) = nz(i32::from(a.n_t_group_number) - i32::from(b.n_t_group_number)) {
                return r;
            }
            if let Some(r) = nz(i32::from(a.n_num_t) - i32::from(b.n_num_t)) {
                return r;
            }
            if let Some(r) = nz(i32::from(a.n_num_d) - i32::from(b.n_num_d)) {
                return r;
            }
            if let Some(r) = nz(i32::from(a.n_num_h) - i32::from(b.n_num_h)) {
                return r;
            }
        }
    }
    0
}

/// The common start of `CompINChITautVsNonTaut` and `CompINChI2` after the
/// Hill formulas: non-tautomeric connection table, H and tautomeric part.
fn compare_basic_layers(i1: &InChI, i2: &InChI, num_h1: i32, num_h2: i32) -> i32 {
    if let Some(r) = nz(i2.n_number_of_atoms - i1.n_number_of_atoms) {
        return r;
    }
    for i in 0..i1.n_number_of_atoms as usize {
        if let Some(r) = nz(i32::from(i2.n_atom[i]) - i32::from(i1.n_atom[i])) {
            return r;
        }
    }
    if let Some(r) = nz(i2.len_conn_table - i1.len_conn_table) {
        return r;
    }
    for i in 0..i2.len_conn_table as usize {
        if let Some(r) = nz(i32::from(i2.n_conn_table[i]) - i32::from(i1.n_conn_table[i])) {
            return r;
        }
    }
    if let Some(r) = nz(num_h2 - num_h1) {
        return r;
    }
    for i in 0..i1.n_number_of_atoms as usize {
        if i2.n_num_h[i] != i1.n_num_h[i] {
            return if i2.n_num_h[i] == 0 {
                1
            } else if i1.n_num_h[i] == 0 {
                -1
            } else {
                i32::from(i2.n_num_h[i]) - i32::from(i1.n_num_h[i])
            };
        }
    }
    compare_taut_non_iso_part_of_inchi(i1, i2)
}

fn compare_charges(i1: &InChI, i2: &InChI) -> i32 {
    if i2.n_total_charge != 0 && i1.n_total_charge != 0 {
        return i1.n_total_charge - i2.n_total_charge;
    }
    i32::from(i1.n_total_charge != 0) - i32::from(i2.n_total_charge != 0)
}

/// `CompINChITautVsNonTaut`
pub fn comp_inchi_taut_vs_non_taut(
    p1: &InChISort,
    p2: &InChISort,
    b_compare_isotopic: bool,
) -> i32 {
    let has = |p: &Option<InChI>| p.as_ref().is_some_and(|i| i.n_number_of_atoms != 0);
    let n1 = if has(&p1.p_inchi[TAUT_YES]) {
        TAUT_YES
    } else {
        TAUT_NON
    };
    let i1 = nonempty(&p1.p_inchi[n1]);
    let i2 = if n1 == TAUT_YES && has(&p2.p_inchi[TAUT_NON]) {
        nonempty(&p2.p_inchi[TAUT_NON])
    } else {
        None
    };
    let (Some(i1), Some(i2)) = (i1, i2) else {
        return 0;
    };
    if i1.b_deleted != 0 {
        return 1;
    }
    if i2.b_deleted != 0 {
        return -1;
    }
    if i1.n_number_of_atoms > 0 && i2.n_number_of_atoms == 0 {
        return 0;
    }
    let (mut num_h1, mut num_h2) = (0, 0);
    let r = compare_hill_formulas_no_h(
        i1.sz_hill_formula.as_deref().unwrap_or(""),
        i2.sz_hill_formula.as_deref().unwrap_or(""),
        &mut num_h1,
        &mut num_h2,
    );
    if r != 0 {
        return r;
    }
    let r = compare_basic_layers(i1, i2, num_h1, num_h2);
    if r != 0 {
        return r;
    }
    if let Some(fixed) = &i2.n_num_h_fixed {
        if fixed[..i2.n_number_of_atoms as usize]
            .iter()
            .any(|&h| h != 0)
        {
            return 1;
        }
    }
    let r = compare_inchi_stereo(
        i1.stereo.as_ref(),
        i1.n_flags,
        i2.stereo.as_ref(),
        i2.n_flags,
    );
    if r != 0 {
        return r;
    }
    if b_compare_isotopic {
        let r = compare_isotopic_part(i1, i2, false);
        if r != 0 {
            return r;
        }
        // the C code returns the value of `a || b` here
        if i2.n_number_of_isotopic_t_groups != 0 || i1.n_number_of_isotopic_t_groups != 0 {
            return 1;
        }
        let r = compare_inchi_stereo(
            i1.stereo_isotopic.as_ref(),
            i1.n_flags,
            i2.stereo_isotopic.as_ref(),
            i2.n_flags,
        );
        if r != 0 {
            return r;
        }
    }
    compare_charges(i1, i2)
}

/// `CompINChI2`: sorting in descending order (-1 if `p1` > `p2`)
pub fn comp_inchi2(p1: &InChISort, p2: &InChISort, b_taut: usize, b_compare_isotopic: bool) -> i32 {
    let has = |p: &Option<InChI>| p.as_ref().is_some_and(|i| i.n_number_of_atoms != 0);
    let n1 = if has(&p1.p_inchi[TAUT_YES]) {
        TAUT_YES
    } else {
        TAUT_NON
    };
    let n2 = if has(&p2.p_inchi[TAUT_YES]) {
        TAUT_YES
    } else {
        TAUT_NON
    };
    let i1 = nonempty(&p1.p_inchi[n1]);
    let i1n = if n1 == TAUT_YES && has(&p1.p_inchi[TAUT_NON]) {
        nonempty(&p1.p_inchi[TAUT_NON])
    } else {
        None
    };
    let i2 = nonempty(&p2.p_inchi[n2]);
    let i2n = if n2 == TAUT_YES && has(&p2.p_inchi[TAUT_NON]) {
        nonempty(&p2.p_inchi[TAUT_NON])
    } else {
        None
    };
    let (a, b) = match (i1, i2) {
        (Some(a), Some(b)) => (a, b),
        (Some(_), None) => return -1,
        (None, Some(_)) => return 1,
        (None, None) => return 0,
    };
    if a.b_deleted != 0 && b.b_deleted == 0 {
        return 1;
    }
    if a.b_deleted == 0 && b.b_deleted != 0 {
        return -1;
    }
    let (mut num_h1, mut num_h2) = (0, 0);
    let r = compare_hill_formulas_no_h(
        a.sz_hill_formula.as_deref().unwrap_or(""),
        b.sz_hill_formula.as_deref().unwrap_or(""),
        &mut num_h1,
        &mut num_h2,
    );
    if r != 0 {
        return r;
    }
    let r = compare_basic_layers(a, b, num_h1, num_h2);
    if r != 0 {
        return r;
    }
    let fixed = |x: Option<&InChI>| x.is_some_and(|i| i.n_num_h_fixed.is_some());
    if b_taut == TAUT_NON && (fixed(i1n) || fixed(i2n)) {
        let f1 = i1n.unwrap_or(a).sz_hill_formula.as_deref();
        let f2 = i2n.unwrap_or(b).sz_hill_formula.as_deref();
        if let (Some(f1), Some(f2)) = (f1, f2) {
            let r = compare_hill_formulas(f1, f2);
            if r != 0 {
                return r;
            }
        }
        if fixed(i1n) && fixed(i2n) {
            let (n1, n2) = (i1n.unwrap(), i2n.unwrap());
            let (h1, h2) = (
                n1.n_num_h_fixed.as_ref().unwrap(),
                n2.n_num_h_fixed.as_ref().unwrap(),
            );
            let num = n1.n_number_of_atoms.min(n2.n_number_of_atoms) as usize;
            for i in 0..num {
                if h2[i] != h1[i] {
                    return if h2[i] == 0 {
                        1
                    } else if h1[i] == 0 {
                        -1
                    } else {
                        i32::from(h2[i]) - i32::from(h1[i])
                    };
                }
            }
            if let Some(r) = nz(n2.n_number_of_atoms - n1.n_number_of_atoms) {
                return r;
            }
        } else if fixed(i1n) {
            let n = i1n.unwrap();
            if n.n_num_h_fixed.as_ref().unwrap()[..n.n_number_of_atoms as usize]
                .iter()
                .any(|&h| h != 0)
            {
                return -1;
            }
        } else {
            let n = i2n.unwrap();
            if n.n_num_h_fixed.as_ref().unwrap()[..n.n_number_of_atoms as usize]
                .iter()
                .any(|&h| h != 0)
            {
                return 1;
            }
        }
    }
    let (mut a, mut b) = (a, b);
    if b_taut == TAUT_NON {
        if let Some(x) = i1n {
            a = x;
        }
        if let Some(x) = i2n {
            b = x;
        }
    }
    let r = compare_inchi_stereo(a.stereo.as_ref(), a.n_flags, b.stereo.as_ref(), b.n_flags);
    if r != 0 {
        return r;
    }
    if b_compare_isotopic {
        let r = compare_isotopic_part(a, b, true);
        if r != 0 {
            return r;
        }
        let r = compare_inchi_stereo(
            a.stereo_isotopic.as_ref(),
            a.n_flags,
            b.stereo_isotopic.as_ref(),
            b.n_flags,
        );
        if r != 0 {
            return r;
        }
    }
    compare_charges(a, b)
}

/// `CompINChINonTaut2` (`CANON_FIXH_TRANS` is 1)
pub fn comp_inchi_non_taut2(p1: &InChISort, p2: &InChISort) -> i32 {
    let mut r = comp_inchi2(p1, p2, TAUT_NON, true);
    if r == 0 {
        r = comp_inchi2(p1, p2, TAUT_YES, true);
    }
    if r == 0 {
        r = i32::from(p1.ord_number) - i32::from(p2.ord_number);
    }
    r
}

/// `CompINChITaut2`
pub fn comp_inchi_taut2(p1: &InChISort, p2: &InChISort) -> i32 {
    let mut r = comp_inchi2(p1, p2, TAUT_YES, true);
    if r == 0 {
        r = comp_inchi2(p1, p2, TAUT_NON, true);
    }
    if r == 0 {
        r = i32::from(p1.ord_number) - i32::from(p2.ord_number);
    }
    r
}
