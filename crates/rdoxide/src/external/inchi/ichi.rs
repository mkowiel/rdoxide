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

pub const EQL_EXISTS: i32 = 1;
pub const EQL_SP3: i32 = 2;
pub const EQL_SP3_INV: i32 = 4;
pub const EQL_SP2: i32 = 8;

/// `Eql_INChI_Stereo` (`ichiprt2.c`; `REL_RAC_STEREO_IGN_1_SC` is 0)
pub fn eql_inchi_stereo(
    s1: Option<&InChIStereo>,
    eql1: i32,
    s2: Option<&InChIStereo>,
    eql2: i32,
) -> bool {
    use super::ichister::atom_parity_well_def;
    let Some(s1) = s1 else {
        return false;
    };
    if eql1 == EQL_SP2 {
        let len = s1.n_number_of_stereo_bonds;
        if len > 0
            && !s1.b_parity.is_empty()
            && !s1.n_bond_atom1.is_empty()
            && !s1.n_bond_atom2.is_empty()
        {
            let Some(s2) = s2 else {
                return eql2 == EQL_EXISTS;
            };
            let len = len as usize;
            return eql2 == EQL_SP2
                && len as i32 == s2.n_number_of_stereo_bonds
                && !s2.b_parity.is_empty()
                && !s2.n_bond_atom1.is_empty()
                && !s2.n_bond_atom2.is_empty()
                && s1.n_bond_atom1[..len] == s2.n_bond_atom1[..len]
                && s1.n_bond_atom2[..len] == s2.n_bond_atom2[..len]
                && s1.b_parity[..len] == s2.b_parity[..len];
        }
        return false;
    }
    let inv1 = eql1 == EQL_SP3_INV;
    if (eql1 == EQL_SP3 || inv1) && s1.n_number_of_stereo_centers > 0 {
        let len = s1.n_number_of_stereo_centers as usize;
        let (t_parity1, n_number1) = if inv1 {
            if s1.n_comp_inv2_abs != 0 {
                (&s1.t_parity_inv, &s1.n_number_inv)
            } else {
                return false;
            }
        } else {
            (&s1.t_parity, &s1.n_number)
        };
        if !t_parity1.is_empty() && !n_number1.is_empty() {
            let Some(s2) = s2 else {
                return eql2 == EQL_EXISTS && (!inv1 || s1.n_comp_inv2_abs != 0);
            };
            let inv2 = eql2 == EQL_SP3_INV;
            if (eql2 == EQL_SP3 || inv2) && len as i32 == s2.n_number_of_stereo_centers {
                let (t_parity2, n_number2) = if inv2 {
                    if s2.n_comp_inv2_abs != 0 && s1.n_comp_inv2_abs != 0 {
                        (&s2.t_parity_inv, &s2.n_number_inv)
                    } else {
                        return false;
                    }
                } else {
                    if inv1 && s2.n_comp_inv2_abs == 0 {
                        return false;
                    }
                    (&s2.t_parity, &s2.n_number)
                };
                if !t_parity2.is_empty() && !n_number2.is_empty() {
                    if inv1 ^ inv2 {
                        let mut num_inv = 0;
                        let mut i = 0;
                        while i < len {
                            if n_number1[i] != n_number2[i] {
                                break;
                            }
                            let (p1, p2) = (i32::from(t_parity1[i]), i32::from(t_parity2[i]));
                            if atom_parity_well_def(p1) && atom_parity_well_def(p2) {
                                if p1 + p2 == 3 {
                                    num_inv += 1;
                                } else {
                                    break;
                                }
                            } else if p1 != p2 {
                                break;
                            }
                            i += 1;
                        }
                        return len == i && num_inv > 0;
                    }
                    return t_parity1[..len] == t_parity2[..len]
                        && n_number1[..len] == n_number2[..len];
                }
            }
        }
    }
    false
}

pub const DIFL_M: usize = 0;
pub const DIFL_MI: usize = 1;
pub const DIFL_F: usize = 2;
pub const DIFL_FI: usize = 3;
pub const DIFL_LENGTH: usize = 4;
pub const DIFS_F_FORMULA: usize = 0;
pub const DIFS_H_H_ATOMS: usize = 2;
pub const DIFS_Q_CHARGE: usize = 3;
pub const DIFS_B_SBONDS: usize = 5;
pub const DIFS_T_SATOMS: usize = 6;
pub const DIFS_M_SP3INV: usize = 7;
pub const DIFS_S_STYPE: usize = 8;
pub const DIFS_I_IATOMS: usize = 9;
pub const DIFS_O_TRANSP: usize = 10;
pub const DIFS_LENGTH: usize = 11;
pub const DIFV_BOTH_EMPTY: u8 = 0;
pub const DIFV_EQL2PRECED: u8 = 1;
pub const DIFV_NEQ2PRECED: u8 = 2;
pub const DIFV_IS_EMPTY: u8 = 4;
pub const DIFV_FI_EQ_MI: u8 = 8;
pub const DIFV_OUTPUT_EMPTY_T: u8 = DIFV_IS_EMPTY;
pub const DIFV_OUTPUT_EMPTY_F: u8 = DIFV_EQL2PRECED | DIFV_NEQ2PRECED | DIFV_FI_EQ_MI;
pub const DIFV_OUTPUT_OMIT_F: u8 = DIFV_NEQ2PRECED | DIFV_IS_EMPTY;
pub const DIFV_OUTPUT_FILL_T: u8 = DIFV_EQL2PRECED | DIFV_NEQ2PRECED | DIFV_FI_EQ_MI;
pub const INCHI_SEGM_OMIT: i32 = 0;
pub const INCHI_SEGM_FILL: i32 = 1;
pub const INCHI_SEGM_EMPTY: i32 = 2;

/// The segment difference table `char sDifSegs[DIFL_LENGTH][DIFS_LENGTH]`
pub type DifSegs = [[u8; DIFS_LENGTH]; DIFL_LENGTH];

const SP3_TYPE: i32 = SP3_ABS | SP3_REL | SP3_RAC;
const SP3_ANY: i32 = SP3_TYPE | SP3_ONLY;

fn live(x: Option<&InChI>) -> Option<&InChI> {
    x.filter(|i| i.b_deleted == 0)
}

/// `CompINChILayers`: marks how the layers of the fixed-H and isotopic
/// variants of a component differ from the preceding ones. The quirks of the
/// C code (e.g. the `/c` flags written into the formula segment) are kept.
#[allow(clippy::unnecessary_unwrap)]
pub fn comp_inchi_layers(
    p1: &InChISort,
    p2: &InChISort,
    s: &mut DifSegs,
    b_fix_transp_charge_bug: i32,
) -> i32 {
    let has = |p: &Option<InChI>| p.as_ref().is_some_and(|i| i.n_number_of_atoms != 0);
    let n1 = if has(&p1.p_inchi[TAUT_YES]) {
        TAUT_YES
    } else {
        TAUT_NON
    };
    let i1 = p1.p_inchi[n1].as_ref();
    let i2 = if n1 == TAUT_YES && has(&p2.p_inchi[TAUT_NON]) {
        p2.p_inchi[TAUT_NON].as_ref()
    } else {
        None
    };
    let nonempty_formula = |i: &InChI| i.sz_hill_formula.as_deref().is_some_and(|f| !f.is_empty());
    let (mut num_h1, mut num_h2) = (0, 0);
    let mut b_rel_rac = [DIFV_BOTH_EMPTY as i32; DIFL_LENGTH];

    // /f
    if let Some(a) = live(i1).filter(|i| nonempty_formula(i)) {
        s[DIFL_M][DIFS_F_FORMULA] |= DIFV_NEQ2PRECED;
        if let Some(b) = live(i2).filter(|i| nonempty_formula(i)) {
            if compare_hill_formulas_no_h(
                a.sz_hill_formula.as_deref().unwrap(),
                b.sz_hill_formula.as_deref().unwrap(),
                &mut num_h1,
                &mut num_h2,
            ) == 0
                && num_h1 == num_h2
            {
                s[DIFL_F][DIFS_F_FORMULA] |= DIFV_EQL2PRECED;
            } else {
                s[DIFL_F][DIFS_F_FORMULA] |= DIFV_NEQ2PRECED;
            }
        } else {
            s[DIFL_F][DIFS_F_FORMULA] |= if i2.is_some() {
                DIFV_IS_EMPTY
            } else {
                DIFV_EQL2PRECED
            };
        }
    } else {
        s[DIFL_M][DIFS_F_FORMULA] |= DIFV_BOTH_EMPTY;
        if live(i2).is_some_and(nonempty_formula) {
            s[DIFL_F][DIFS_F_FORMULA] |= DIFV_NEQ2PRECED;
        } else {
            s[DIFL_F][DIFS_F_FORMULA] |= DIFV_BOTH_EMPTY;
        }
    }
    // /c
    if live(i1).is_some_and(|i| i.len_conn_table > 1) {
        s[DIFL_M][DIFS_F_FORMULA] |= DIFV_NEQ2PRECED;
    } else {
        s[DIFL_M][DIFS_F_FORMULA] |= DIFV_BOTH_EMPTY;
    }
    // /h
    if let Some(a) = live(i1) {
        let mut h = i32::from(a.len_tautomer > 0 && a.n_tautomer.first().is_some_and(|&x| x != 0));
        if h == 0 {
            h = i32::from(
                a.n_num_h[..a.n_number_of_atoms as usize]
                    .iter()
                    .any(|&x| x != 0),
            );
        }
        s[DIFL_M][DIFS_H_H_ATOMS] |= if h != 0 {
            DIFV_NEQ2PRECED
        } else {
            DIFV_BOTH_EMPTY
        };
    } else {
        s[DIFL_M][DIFS_H_H_ATOMS] |= DIFV_BOTH_EMPTY;
    }
    if let Some(b) = live(i2).filter(|i| i.n_num_h_fixed.is_some()) {
        let mut h = 0;
        if let Some(a) = live(i1) {
            let fixed = b.n_num_h_fixed.as_ref().unwrap();
            h = i32::from(
                (0..a.n_number_of_atoms as usize).any(|i| fixed.get(i).copied().unwrap_or(0) != 0),
            );
        }
        s[DIFL_F][DIFS_H_H_ATOMS] |= if h != 0 {
            DIFV_NEQ2PRECED
        } else {
            DIFV_BOTH_EMPTY
        };
    } else {
        s[DIFL_F][DIFS_H_H_ATOMS] |= DIFV_BOTH_EMPTY;
    }

    // /q
    if let Some(a) = live(i1) {
        s[DIFL_M][DIFS_Q_CHARGE] |= if a.n_total_charge != 0 {
            DIFV_NEQ2PRECED
        } else {
            DIFV_BOTH_EMPTY
        };
        let ch = |c1: i32, c2: i32| -> u8 {
            if c1 != 0 {
                if c1 == c2 {
                    DIFV_EQL2PRECED
                } else if c2 != 0 {
                    DIFV_NEQ2PRECED
                } else {
                    DIFV_IS_EMPTY
                }
            } else if c2 != 0 {
                DIFV_NEQ2PRECED
            } else {
                DIFV_BOTH_EMPTY
            }
        };
        if let Some(b) = live(i2) {
            s[DIFL_F][DIFS_Q_CHARGE] |= ch(a.n_total_charge, b.n_total_charge);
        } else if i2.is_none() {
            let same = if a.n_total_charge != 0 {
                DIFV_EQL2PRECED
            } else {
                DIFV_BOTH_EMPTY
            };
            if b_fix_transp_charge_bug == 1 {
                let alt = p2.p_inchi[TAUT_YES]
                    .as_ref()
                    .filter(|x| x.b_deleted == 0 && x.n_number_of_atoms != 0);
                if p1.ord_number != p2.ord_number && n1 == TAUT_YES && alt.is_some() {
                    s[DIFL_F][DIFS_Q_CHARGE] |= ch(a.n_total_charge, alt.unwrap().n_total_charge);
                } else {
                    s[DIFL_F][DIFS_Q_CHARGE] |= same;
                }
            } else {
                s[DIFL_F][DIFS_Q_CHARGE] |= same;
            }
        } else {
            s[DIFL_F][DIFS_Q_CHARGE] |= if a.n_total_charge != 0 {
                DIFV_IS_EMPTY
            } else {
                DIFV_BOTH_EMPTY
            };
        }
    } else {
        s[DIFL_M][DIFS_Q_CHARGE] |= DIFV_BOTH_EMPTY;
        if let Some(b) = live(i2) {
            s[DIFL_F][DIFS_Q_CHARGE] |= if b.n_total_charge != 0 {
                DIFV_NEQ2PRECED
            } else {
                DIFV_BOTH_EMPTY
            };
        }
    }

    // stereo
    let (stereo1, iso_stereo1) = match live(i1) {
        Some(a) => (a.stereo.as_ref(), a.stereo_isotopic.as_ref()),
        None => (None, None),
    };
    let (stereo2, iso_stereo2) = match live(i2) {
        Some(b) => (b.stereo.as_ref(), b.stereo_isotopic.as_ref()),
        None => (None, None),
    };
    let sb = |x: Option<&InChIStereo>| x.is_some_and(|s| s.n_number_of_stereo_bonds != 0);

    // /b
    s[DIFL_M][DIFS_B_SBONDS] |= if sb(stereo1) {
        DIFV_NEQ2PRECED
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_F][DIFS_B_SBONDS] |= if sb(stereo2) {
        if sb(stereo1) {
            if eql_inchi_stereo(stereo1, EQL_SP2, stereo2, EQL_SP2) {
                DIFV_EQL2PRECED
            } else {
                DIFV_NEQ2PRECED
            }
        } else {
            DIFV_NEQ2PRECED
        }
    } else if sb(stereo1) {
        if i2.is_some() {
            DIFV_IS_EMPTY
        } else {
            DIFV_EQL2PRECED
        }
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_MI][DIFS_B_SBONDS] |= if sb(iso_stereo1) {
        if eql_inchi_stereo(iso_stereo1, EQL_SP2, stereo1, EQL_SP2) {
            DIFV_EQL2PRECED
        } else {
            DIFV_NEQ2PRECED
        }
    } else if sb(stereo1) {
        DIFV_EQL2PRECED
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_FI][DIFS_B_SBONDS] |= if sb(iso_stereo2) {
        if eql_inchi_stereo(iso_stereo2, EQL_SP2, stereo2, EQL_SP2) {
            DIFV_EQL2PRECED
        } else if !sb(stereo1)
            && !sb(stereo2)
            && eql_inchi_stereo(iso_stereo2, EQL_SP2, iso_stereo1, EQL_SP2)
        {
            DIFV_FI_EQ_MI
        } else {
            DIFV_NEQ2PRECED
        }
    } else if sb(stereo2) {
        DIFV_EQL2PRECED
    } else if sb(iso_stereo1) && !sb(stereo1) {
        if i2.is_some() {
            DIFV_IS_EMPTY
        } else {
            DIFV_EQL2PRECED
        }
    } else {
        DIFV_BOTH_EMPTY
    };

    // /t, /m, /s for M
    b_rel_rac[DIFL_M] = get_sp3_rel_rac_abs(i1, stereo1);
    b_rel_rac[DIFL_MI] = get_sp3_rel_rac_abs(i1, iso_stereo1);
    b_rel_rac[DIFL_F] = get_sp3_rel_rac_abs(i2, stereo2);
    b_rel_rac[DIFL_FI] = get_sp3_rel_rac_abs(i2, iso_stereo2);
    let br = b_rel_rac;
    let seg = |bit: u8, cond: bool| if cond { bit } else { DIFV_BOTH_EMPTY };
    if br[DIFL_M] != SP3_NONE {
        s[DIFL_M][DIFS_T_SATOMS] |= seg(DIFV_NEQ2PRECED, br[DIFL_M] & SP3_ANY != 0);
        s[DIFL_M][DIFS_M_SP3INV] |= seg(DIFV_NEQ2PRECED, br[DIFL_M] & SP3_ABS != 0);
        s[DIFL_M][DIFS_S_STYPE] |= seg(DIFV_NEQ2PRECED, br[DIFL_M] & SP3_TYPE != 0);
    } else {
        s[DIFL_M][DIFS_T_SATOMS] |= DIFV_BOTH_EMPTY;
        s[DIFL_M][DIFS_M_SP3INV] |= DIFV_BOTH_EMPTY;
        s[DIFL_M][DIFS_S_STYPE] |= DIFV_BOTH_EMPTY;
    }
    // /t
    s[DIFL_F][DIFS_T_SATOMS] |= if SP3_ANY & br[DIFL_F] != 0 {
        if eql_inchi_stereo(stereo2, EQL_SP3, stereo1, EQL_SP3) {
            DIFV_EQL2PRECED
        } else {
            DIFV_NEQ2PRECED
        }
    } else if SP3_ANY & br[DIFL_M] != 0 {
        if i2.is_some() {
            DIFV_IS_EMPTY
        } else {
            DIFV_EQL2PRECED
        }
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_MI][DIFS_T_SATOMS] |= if SP3_ANY & br[DIFL_MI] != 0 {
        if eql_inchi_stereo(iso_stereo1, EQL_SP3, stereo1, EQL_SP3) {
            DIFV_EQL2PRECED
        } else {
            DIFV_NEQ2PRECED
        }
    } else if SP3_ANY & br[DIFL_M] != 0 {
        DIFV_EQL2PRECED
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_FI][DIFS_T_SATOMS] |= if SP3_ANY & br[DIFL_FI] != 0 {
        if eql_inchi_stereo(iso_stereo2, EQL_SP3, stereo2, EQL_SP3) {
            DIFV_EQL2PRECED
        } else if SP3_ANY & br[DIFL_M] == 0
            && SP3_ANY & br[DIFL_F] == 0
            && eql_inchi_stereo(iso_stereo2, EQL_SP3, iso_stereo1, EQL_SP3)
        {
            DIFV_FI_EQ_MI
        } else {
            DIFV_NEQ2PRECED
        }
    } else if SP3_ANY & br[DIFL_F] != 0 {
        DIFV_EQL2PRECED
    } else if SP3_ANY & br[DIFL_MI] != 0 && SP3_ANY & br[DIFL_M] == 0 {
        if i2.is_some() {
            DIFV_IS_EMPTY
        } else {
            DIFV_EQL2PRECED
        }
    } else {
        DIFV_BOTH_EMPTY
    };
    // /m
    s[DIFL_F][DIFS_M_SP3INV] |=
        if br[DIFL_F] & SP3_ABS != 0 && stereo1.is_some() && stereo2.is_some() {
            if br[DIFL_M] & SP3_ABS == 0
                || stereo2.unwrap().n_comp_inv2_abs != stereo1.unwrap().n_comp_inv2_abs
            {
                DIFV_NEQ2PRECED
            } else {
                DIFV_EQL2PRECED
            }
        } else if br[DIFL_M] & SP3_ABS != 0 {
            if i2.is_some() {
                DIFV_IS_EMPTY
            } else {
                DIFV_EQL2PRECED
            }
        } else {
            DIFV_BOTH_EMPTY
        };
    s[DIFL_MI][DIFS_M_SP3INV] |= if SP3_ABS & br[DIFL_MI] != 0 && stereo1.is_some() {
        if SP3_ABS & br[DIFL_M] != 0
            && iso_stereo1.unwrap().n_comp_inv2_abs == stereo1.unwrap().n_comp_inv2_abs
        {
            DIFV_EQL2PRECED
        } else {
            DIFV_NEQ2PRECED
        }
    } else if SP3_ABS & br[DIFL_M] != 0 {
        DIFV_EQL2PRECED
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_FI][DIFS_M_SP3INV] |= if SP3_ABS & br[DIFL_FI] != 0
        && iso_stereo1.is_some()
        && iso_stereo2.is_some()
        && stereo2.is_some()
    {
        if SP3_ABS & br[DIFL_F] != 0
            && iso_stereo2.unwrap().n_comp_inv2_abs == stereo2.unwrap().n_comp_inv2_abs
        {
            DIFV_EQL2PRECED
        } else if SP3_ABS & br[DIFL_M] == 0
            && SP3_ABS & br[DIFL_F] == 0
            && SP3_ABS & br[DIFL_MI] != 0
            && iso_stereo2.unwrap().n_comp_inv2_abs == iso_stereo1.unwrap().n_comp_inv2_abs
        {
            DIFV_FI_EQ_MI
        } else {
            DIFV_NEQ2PRECED
        }
    } else if SP3_ABS & br[DIFL_F] != 0 {
        DIFV_EQL2PRECED
    } else if SP3_ABS & br[DIFL_MI] != 0 && SP3_ABS & br[DIFL_M] == 0 {
        if i2.is_some() {
            DIFV_IS_EMPTY
        } else {
            DIFV_EQL2PRECED
        }
    } else {
        DIFV_BOTH_EMPTY
    };
    // /s
    s[DIFL_F][DIFS_S_STYPE] |= if br[DIFL_F] & SP3_TYPE != 0 {
        if br[DIFL_F] & SP3_TYPE == br[DIFL_M] & SP3_TYPE {
            DIFV_EQL2PRECED
        } else {
            DIFV_NEQ2PRECED
        }
    } else if br[DIFL_M] & SP3_TYPE != 0 {
        if i2.is_some() {
            DIFV_IS_EMPTY
        } else {
            DIFV_EQL2PRECED
        }
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_MI][DIFS_S_STYPE] |= if SP3_TYPE & br[DIFL_MI] != 0 {
        if SP3_TYPE & br[DIFL_MI] == SP3_TYPE & br[DIFL_M] {
            DIFV_EQL2PRECED
        } else {
            DIFV_NEQ2PRECED
        }
    } else if SP3_TYPE & br[DIFL_M] != 0 {
        DIFV_EQL2PRECED
    } else {
        DIFV_BOTH_EMPTY
    };
    s[DIFL_FI][DIFS_S_STYPE] |= if SP3_TYPE & br[DIFL_FI] != 0 {
        if SP3_TYPE & br[DIFL_FI] == SP3_TYPE & br[DIFL_F] {
            DIFV_EQL2PRECED
        } else if SP3_TYPE & br[DIFL_M] == 0
            && SP3_TYPE & br[DIFL_F] == 0
            && SP3_TYPE & br[DIFL_MI] != 0
        {
            DIFV_FI_EQ_MI
        } else {
            DIFV_NEQ2PRECED
        }
    } else if SP3_TYPE & br[DIFL_F] != 0 {
        DIFV_EQL2PRECED
    } else if SP3_TYPE & br[DIFL_MI] != 0 && SP3_TYPE & br[DIFL_M] == 0 {
        if i2.is_some() {
            DIFV_IS_EMPTY
        } else {
            DIFV_EQL2PRECED
        }
    } else {
        DIFV_BOTH_EMPTY
    };
    // /o
    if p1.ord_number != p2.ord_number {
        s[DIFL_F][DIFS_O_TRANSP] |= DIFV_NEQ2PRECED;
    }
    // /i
    let iso = |x: Option<&InChI>| {
        x.is_some_and(|i| i.n_number_of_isotopic_atoms != 0 || i.n_number_of_isotopic_t_groups != 0)
    };
    s[DIFL_MI][DIFS_I_IATOMS] |= if live(i1).is_some_and(|i| iso(Some(i))) {
        DIFV_NEQ2PRECED
    } else {
        DIFV_BOTH_EMPTY
    };
    if let Some(b) = live(i2) {
        if iso(Some(b)) {
            let a = live(i1);
            if a.is_none_or(|a| {
                b.n_number_of_isotopic_atoms != a.n_number_of_isotopic_atoms
                    || b.n_number_of_isotopic_t_groups != a.n_number_of_isotopic_t_groups
            }) {
                s[DIFL_FI][DIFS_I_IATOMS] |= DIFV_NEQ2PRECED;
            } else {
                let a = a.unwrap();
                let mut diff = 0;
                for i in 0..a.n_number_of_isotopic_atoms as usize {
                    let (x, y) = (&b.isotopic_atom[i], &a.isotopic_atom[i]);
                    diff = i32::from(x.n_atom_number) - i32::from(y.n_atom_number);
                    if diff == 0 {
                        diff = i32::from(x.n_iso_difference) - i32::from(y.n_iso_difference);
                    }
                    if diff == 0 {
                        diff = i32::from(x.n_num_t) - i32::from(y.n_num_t);
                    }
                    if diff == 0 {
                        diff = i32::from(x.n_num_d) - i32::from(y.n_num_d);
                    }
                    if diff == 0 {
                        diff = i32::from(x.n_num_h) - i32::from(y.n_num_h);
                    }
                    if diff != 0 {
                        break;
                    }
                }
                if diff == 0 {
                    for i in 0..a.n_number_of_isotopic_t_groups as usize {
                        let (x, y) = (&b.isotopic_t_group[i], &a.isotopic_t_group[i]);
                        diff = i32::from(x.n_t_group_number) - i32::from(y.n_t_group_number);
                        if diff != 0 {
                            break;
                        }
                        diff = i32::from(x.n_num_t) - i32::from(y.n_num_t);
                        if diff != 0 {
                            break;
                        }
                        diff = i32::from(x.n_num_d) - i32::from(y.n_num_d);
                        if diff != 0 {
                            // the C code returns here, skipping the final update
                            return diff;
                        }
                        diff = i32::from(x.n_num_h) - i32::from(y.n_num_h);
                        if diff != 0 {
                            break;
                        }
                    }
                }
                s[DIFL_FI][DIFS_I_IATOMS] |= if diff != 0 {
                    DIFV_NEQ2PRECED
                } else {
                    DIFV_FI_EQ_MI
                };
            }
        } else if live(i1).is_some_and(|i| iso(Some(i))) {
            s[DIFL_FI][DIFS_I_IATOMS] |= DIFV_IS_EMPTY;
        }
    } else if i2.is_none() {
        s[DIFL_FI][DIFS_I_IATOMS] |= if live(i1).is_some_and(|i| iso(Some(i))) {
            DIFV_EQL2PRECED
        } else {
            DIFV_BOTH_EMPTY
        };
    }
    0
}

/// `INChI_SegmentAction`
pub fn inchi_segment_action(c_dif_segs: u8) -> i32 {
    if c_dif_segs & DIFV_OUTPUT_OMIT_F == 0 {
        return INCHI_SEGM_OMIT;
    }
    if c_dif_segs & DIFV_OUTPUT_EMPTY_T != 0 && c_dif_segs & DIFV_OUTPUT_EMPTY_F == 0 {
        return INCHI_SEGM_EMPTY;
    }
    if c_dif_segs & DIFV_OUTPUT_FILL_T != 0 {
        return INCHI_SEGM_FILL;
    }
    INCHI_SEGM_OMIT
}

/// `MarkUnusedAndEmptyLayers`
pub fn mark_unused_and_empty_layers(s: &mut DifSegs) -> i32 {
    const FIRST_FML_SEGM: usize = DIFS_F_FORMULA;
    const FIRST_ISO_SEGM: usize = DIFS_I_IATOMS;
    for (layer, first) in [(DIFL_FI, FIRST_ISO_SEGM), (DIFL_MI, FIRST_ISO_SEGM)] {
        let bits = s[layer].iter().fold(0, |a, &b| a | b);
        if bits & DIFV_OUTPUT_OMIT_F == 0 {
            s[layer] = [DIFV_BOTH_EMPTY; DIFS_LENGTH];
        } else if s[layer][first] == DIFV_BOTH_EMPTY || s[layer][first] & DIFV_OUTPUT_OMIT_F == 0 {
            s[layer][first] = DIFV_IS_EMPTY;
        }
    }
    let bits = s[DIFL_F].iter().fold(0, |a, &b| a | b);
    if bits & DIFV_OUTPUT_OMIT_F == 0 && s[DIFL_FI][FIRST_ISO_SEGM] == DIFV_BOTH_EMPTY {
        s[DIFL_F] = [DIFV_BOTH_EMPTY; DIFS_LENGTH];
    } else if s[DIFL_F][FIRST_FML_SEGM] == DIFV_BOTH_EMPTY
        || s[DIFL_F][FIRST_FML_SEGM] & DIFV_OUTPUT_OMIT_F == 0
    {
        s[DIFL_F][FIRST_FML_SEGM] = DIFV_IS_EMPTY;
    }
    0
}

/// `ORIG_INFO`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct OrigInfo {
    pub c_charge: i8,
    /// 0=none, 1=doublet, 2=triplet, 3=unknown
    pub c_radical: i8,
    pub c_unusual_valence: i8,
}

/// `INChI_Aux` (without the coordinate buffer)
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InChIAux {
    pub n_error_code: i32,
    pub n_number_of_atoms: i32,
    pub n_number_of_t_groups: i32,
    pub b_is_isotopic: i32,
    pub b_is_tautomeric: i32,
    pub n_orig_at_nos_in_canon_ord: Vec<u16>,
    pub n_isotopic_orig_at_nos_in_canon_ord: Vec<u16>,
    pub n_orig_at_nos_in_canon_ord_inv: Vec<u16>,
    pub n_isotopic_orig_at_nos_in_canon_ord_inv: Vec<u16>,
    pub n_constit_equ_numbers: Vec<u16>,
    pub n_constit_equ_t_group_numbers: Vec<u16>,
    pub n_constit_equ_isotopic_numbers: Vec<u16>,
    pub n_constit_equ_isotopic_t_group_numbers: Vec<u16>,
    pub orig_info: Vec<OrigInfo>,
    pub n_num_removed_protons: i16,
    pub n_num_removed_isotopic_h: [i16; 3],
    pub b_deleted: i32,
    pub b_taut_flags: u64,
    pub b_taut_flags_done: u64,
    pub b_normalization_flags: u64,
    pub n_canon_flags: i32,
}
