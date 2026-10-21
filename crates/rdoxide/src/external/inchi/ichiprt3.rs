//! Substrings of the whole-structure InChI string (`INCHI_BASE/src/ichiprt3.c`):
//! the formula, connection, hydrogen and charge layers.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::{InChI, InChISort, compare_taut_non_iso_part_of_inchi};
use super::ichi_io::{InchiIosString, inchi_strbuf_printf};
use super::ichiprt1::{II_EQU, IIT_NONTAUT, S_COMP_DELIM, equ_string};
use super::ichiprt2::{make_delim, make_eq_str, make_h_string, make_mult, make_taut_string};

pub const OUT_N1: i32 = 0;
pub const OUT_T1: i32 = 1;
pub const OUT_NT: i32 = 2;
pub const OUT_TN: i32 = 3;
pub const OUT_NN: i32 = 4;

const TAUT_NON: usize = 0;
const TAUT_YES: usize = 1;

fn has_t(s: &InChISort) -> bool {
    s.p_inchi[TAUT_YES]
        .as_ref()
        .is_some_and(|c| c.n_number_of_atoms != 0)
}
fn has_n(s: &InChISort) -> bool {
    s.p_inchi[TAUT_NON]
        .as_ref()
        .is_some_and(|c| c.n_number_of_atoms != 0)
}
fn has_tt(s: &InChISort) -> bool {
    s.p_inchi[TAUT_YES]
        .as_ref()
        .is_some_and(|c| c.n_number_of_atoms != 0 && c.len_tautomer > 0)
}
fn has_tn(s: &InChISort) -> bool {
    s.p_inchi[TAUT_YES]
        .as_ref()
        .is_some_and(|c| c.n_number_of_atoms != 0 && c.len_tautomer == 0)
}
fn has_nn(s: &InChISort) -> bool {
    s.p_inchi[TAUT_NON]
        .as_ref()
        .is_some_and(|c| c.n_number_of_atoms != 0 && c.len_tautomer == 0)
}

/// `GET_II`: the index of the `INChI` to output, -1 for none
pub fn get_ii(mode: i32, s: &InChISort) -> i32 {
    let pick = |a: bool, ia: i32, b: bool, ib: i32| {
        if a {
            ia
        } else if b {
            ib
        } else {
            -1
        }
    };
    match mode {
        OUT_N1 => pick(has_tn(s), 1, has_nn(s), 0),
        OUT_T1 | OUT_TN => pick(has_t(s), 1, has_n(s), 0),
        OUT_NN => pick(has_nn(s), 0, has_tn(s), 1),
        OUT_NT if has_tt(s) && has_nn(s) => 0,
        _ => -1,
    }
}

/// the `INChI` of the component `i` (`NULL` outside `num_components` or when
/// nothing is to be output)
fn pick(is0: &[InChISort], i: usize, num_components: usize, mode: i32) -> Option<&InChI> {
    if i < num_components {
        let ii = get_ii(mode, &is0[i]);
        if ii >= 0 {
            return is0[i].p_inchi[ii as usize].as_ref();
        }
    }
    None
}

fn first(is0: &[InChISort], mode: i32) -> Option<&InChI> {
    let ii = get_ii(mode, &is0[0]);
    if ii >= 0 {
        is0[0].p_inchi[ii as usize].as_ref()
    } else {
        None
    }
}

fn used(b: &InchiIosString) -> i32 {
    b.n_used_length
}

/// truncates the buffer to `len` and terminates it
fn truncate(b: &mut InchiIosString, len: i32) {
    b.n_used_length = len;
    if let Some(p) = b.p_str.as_mut() {
        p[len as usize] = 0;
    }
}

fn make_hill_formula_string(f: &str, strbuf: &mut InchiIosString, overflow: &mut i32) -> i32 {
    let u0 = used(strbuf);
    if *overflow == 0 {
        inchi_strbuf_printf(strbuf, f);
    }
    used(strbuf) - u0
}

fn nonempty_formula(c: Option<&InChI>) -> Option<&str> {
    c.and_then(|c| c.sz_hill_formula.as_deref())
        .filter(|f| !f.is_empty())
}

/// `str_HillFormula`
pub fn str_hill_formula(
    is0: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    num_components: usize,
    use_multipliers: bool,
) -> i32 {
    let u0 = used(strbuf);
    if is0.is_empty() {
        return u0;
    }
    let mut prev = first(is0, out_type);
    let (mut mult, mut next) = (0, 0);
    for i in 1..=num_components {
        let cur = pick(is0, i, num_components, out_type);
        let eq2prev = use_multipliers
            && matches!((cur, prev), (Some(c), Some(p))
                if c.sz_hill_formula.is_some() && p.sz_hill_formula.is_some()
                    && !c.sz_hill_formula.as_ref().unwrap().is_empty()
                    && c.sz_hill_formula == p.sz_hill_formula);
        if eq2prev {
            mult += 1;
            continue;
        }
        if next != 0 {
            make_delim(Some("."), strbuf, overflow);
        }
        next += 1;
        if let Some(f) = nonempty_formula(prev) {
            make_mult(mult + 1, "", strbuf, 0, overflow);
            make_hill_formula_string(f, strbuf, overflow);
        }
        prev = cur;
        mult = 0;
    }
    used(strbuf) - u0
}

/// `str_HillFormula2`: the second (FixedH) appearance
pub fn str_hill_formula2(
    is0: &[InChISort],
    is20: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    num_components: usize,
    use_multipliers: bool,
) -> i32 {
    let u0 = used(strbuf);
    let mut prev = first(is0, out_type);
    let mut taut_prev = first(is20, OUT_T1);
    let (mut mult, mut next) = (0, 0);
    let same = |a: Option<&InChI>, b: Option<&InChI>| {
        matches!((a, b), (Some(a), Some(b)) if !b.b_deleted_nz()
            && a.sz_hill_formula.is_some() && b.sz_hill_formula.is_some()
            && a.sz_hill_formula == b.sz_hill_formula)
    };
    let mut eq_to_taut = same(prev, taut_prev);
    for i in 1..=num_components {
        let cur = pick(is0, i, num_components, out_type);
        let taut = pick(is20, i, num_components, OUT_T1);
        if eq_to_taut && (cur.is_some() || taut.is_some()) {
            eq_to_taut = same(cur, taut);
        }
        let eq2prev = use_multipliers
            && matches!((cur, prev), (Some(c), Some(p))
                if c.sz_hill_formula.is_some() && p.sz_hill_formula.is_some()
                    && !c.sz_hill_formula.as_ref().unwrap().is_empty()
                    && c.sz_hill_formula == p.sz_hill_formula);
        if eq2prev {
            mult += 1;
            continue;
        }
        if next != 0 {
            make_delim(Some("."), strbuf, overflow);
        }
        next += 1;
        if let Some(f) = nonempty_formula(prev) {
            make_mult(mult + 1, "", strbuf, 0, overflow);
            make_hill_formula_string(f, strbuf, overflow);
        }
        prev = cur;
        taut_prev = taut;
        mult = 0;
    }
    let _ = taut_prev;
    if eq_to_taut {
        truncate(strbuf, u0);
    } else {
        truncate(strbuf, used(strbuf));
    }
    used(strbuf) - u0
}

trait Deleted {
    fn b_deleted_nz(&self) -> bool;
}
impl Deleted for InChI {
    fn b_deleted_nz(&self) -> bool {
        self.b_deleted != 0
    }
}

/// `str_Connections`
pub fn str_connections(
    is0: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    atom_mode: i32,
    num_components: usize,
    use_multipliers: bool,
) -> i32 {
    use super::ichiprt2::make_ct_string_new;
    let u0 = used(strbuf);
    if is0.is_empty() {
        return u0;
    }
    let mut prev = first(is0, out_type);
    let (mut mult, mut next, mut num_empty) = (0, 0, 0usize);
    for i in 1..=num_components {
        let cur = pick(is0, i, num_components, out_type);
        let eq2prev = use_multipliers
            && matches!((cur, prev), (Some(c), Some(p))
                if c.len_conn_table > 1 && p.len_conn_table == c.len_conn_table
                    && p.n_conn_table[..p.len_conn_table as usize] == c.n_conn_table[..p.len_conn_table as usize]);
        if eq2prev {
            mult += 1;
            continue;
        }
        if let Some(p) = prev {
            if next != 0 {
                make_delim(Some(S_COMP_DELIM), strbuf, overflow);
            }
            next += 1;
            if p.len_conn_table > 1 {
                make_mult(mult + 1, "*", strbuf, 0, overflow);
                make_ct_string_new(
                    &p.n_conn_table,
                    p.len_conn_table,
                    false,
                    None,
                    p.n_number_of_atoms.max(0) as usize,
                    strbuf,
                    atom_mode,
                    overflow,
                );
            } else {
                num_empty += 1;
            }
        }
        prev = cur;
        mult = 0;
    }
    if num_empty == num_components && used(strbuf) > u0 {
        truncate(strbuf, u0);
    }
    used(strbuf) - u0
}

/// `str_H_atoms`
pub fn str_h_atoms(
    is0: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    atom_mode: i32,
    taut_mode: i32,
    num_components: usize,
    use_multipliers: bool,
) -> i32 {
    let u0 = used(strbuf);
    let mut num_empty = 0usize;
    let mut prev = first(is0, out_type);
    let (mut mult, mut next) = (0, 0);
    for i in 1..=num_components {
        let cur = pick(is0, i, num_components, out_type);
        let mut eq2prev = use_multipliers
            && matches!((cur, prev), (Some(c), Some(p))
                if (c.n_number_of_atoms > 0 || c.len_tautomer > 1)
                    && p.n_number_of_atoms == c.n_number_of_atoms
                    && (p.n_number_of_atoms == 0
                        || p.n_num_h[..p.n_number_of_atoms as usize] == c.n_num_h[..p.n_number_of_atoms as usize])
                    && compare_taut_non_iso_part_of_inchi(p, c) == 0);
        if eq2prev {
            let p = prev.unwrap();
            if p.len_tautomer <= 1 {
                eq2prev = p.n_num_h[..p.n_number_of_atoms.max(0) as usize]
                    .iter()
                    .any(|&x| x != 0);
            }
        }
        if eq2prev {
            mult += 1;
            continue;
        }
        if let Some(p) = prev {
            if next != 0 {
                make_delim(Some(S_COMP_DELIM), strbuf, overflow);
            }
            next += 1;
            let not_empty = p.len_tautomer > 1
                || p.n_num_h[..p.n_number_of_atoms.max(0) as usize]
                    .iter()
                    .any(|&x| x != 0);
            if not_empty {
                make_mult(mult + 1, "*", strbuf, 0, overflow);
                let len_h = make_h_string(
                    false,
                    &p.n_num_h,
                    p.n_number_of_atoms,
                    strbuf,
                    atom_mode,
                    overflow,
                );
                make_taut_string(
                    &p.n_tautomer,
                    p.len_tautomer,
                    len_h != 0,
                    strbuf,
                    taut_mode,
                    overflow,
                );
            } else {
                num_empty += 1;
            }
        }
        prev = cur;
        mult = 0;
    }
    if num_empty == num_components && used(strbuf) > u0 {
        truncate(strbuf, u0);
    }
    used(strbuf) - u0
}

/// `str_Charge2`
#[allow(clippy::too_many_arguments)]
pub fn str_charge2(
    is0: &[InChISort],
    is20_in: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    num_components: usize,
    second_non_taut_pass: bool,
    omit_repetitions: bool,
    use_multipliers: bool,
) -> i32 {
    let is20: &[InChISort] = if second_non_taut_pass { is20_in } else { &[] };
    let mut taut: Option<&InChI> = None;
    let mut prev: Option<&InChI> = None;
    let mut taut_prev: Option<&InChI> = None;
    let (mut mult, mut next) = (0, 0);
    let mut eq2taut_prev = true;
    let mut prev_equ: Option<&str> = None;
    let mut mult_prev_equ = 0;
    let u0 = used(strbuf);
    for i in 0..=num_components {
        let cur = pick(is0, i, num_components, out_type);
        if second_non_taut_pass {
            taut = pick(is20, i, num_components, OUT_T1);
        }
        let mut eq2taut = false;
        if second_non_taut_pass && omit_repetitions {
            eq2taut = matches!((cur, taut), (Some(c), Some(t)) if t.b_deleted == 0
                && c.n_total_charge != 0 && t.n_total_charge != 0 && c.n_total_charge == t.n_total_charge);
        }
        if eq2taut {
            if let Some(p) = prev.filter(|p| p.n_number_of_atoms != 0) {
                if next != 0 {
                    make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                }
                next += 1;
                if p.n_total_charge != 0 {
                    make_mult(mult + 1, "*", strbuf, 0, overflow);
                    inchi_strbuf_printf(strbuf, &format!("{:+}", p.n_total_charge));
                }
            } else if taut_prev.is_some_and(|t| t.n_number_of_atoms != 0 && t.b_deleted == 0) {
                if next != 0 {
                    make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                }
                next += 1;
            }
            let curr_equ = equ_string(II_EQU | IIT_NONTAUT);
            if mult_prev_equ != 0 && prev_equ.is_some() {
                if prev_equ == Some(curr_equ) {
                    mult_prev_equ += 1;
                } else {
                    if next != 0 {
                        make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                    }
                    next += 1;
                    make_eq_str(prev_equ, mult_prev_equ, strbuf, overflow);
                    prev_equ = Some(curr_equ);
                    mult_prev_equ = 1;
                }
            } else {
                prev_equ = Some(curr_equ);
                mult_prev_equ = 1;
            }
            prev = None;
            taut_prev = None;
            mult = 0;
            eq2taut_prev = true;
        } else if eq2taut_prev {
            if mult_prev_equ != 0 && prev_equ.is_some() {
                if next != 0 {
                    make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                }
                next += 1;
                make_eq_str(prev_equ, mult_prev_equ, strbuf, overflow);
                prev_equ = None;
                mult_prev_equ = 0;
            }
            eq2taut_prev = false;
            prev = cur;
            taut_prev = taut;
            mult = 0;
        } else {
            let eq2prev = use_multipliers
                && matches!((cur, prev), (Some(c), Some(p))
                    if c.n_total_charge != 0 && p.n_total_charge != 0 && c.n_total_charge == p.n_total_charge);
            if eq2prev {
                mult += 1;
                continue;
            }
            if next != 0 {
                make_delim(Some(S_COMP_DELIM), strbuf, overflow);
            }
            next += 1;
            if let Some(p) = prev.filter(|p| p.n_number_of_atoms != 0)
                && p.n_total_charge != 0
            {
                make_mult(mult + 1, "*", strbuf, 0, overflow);
                inchi_strbuf_printf(strbuf, &format!("{:+}", p.n_total_charge));
            }
            prev = cur;
            taut_prev = taut;
            mult = 0;
        }
    }
    used(strbuf) - u0
}

/// `str_FixedH_atoms`
pub fn str_fixed_h_atoms(
    is0: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    atom_mode: i32,
    num_components: usize,
    use_multipliers: bool,
) -> i32 {
    let fixed = |c: &InChI| -> Vec<i8> { c.n_num_h_fixed.clone().unwrap_or_default() };
    let mut prev = first(is0, out_type);
    let (mut mult, mut next, mut num_empty) = (0, 0, 0usize);
    let u0 = used(strbuf);
    for i in 1..=num_components {
        let cur = pick(is0, i, num_components, out_type);
        let mut eq2prev = use_multipliers
            && matches!((cur, prev), (Some(c), Some(p))
                if c.n_number_of_atoms > 0 && p.n_number_of_atoms == c.n_number_of_atoms
                    && fixed(p)[..p.n_number_of_atoms as usize] == fixed(c)[..p.n_number_of_atoms as usize]);
        if eq2prev {
            let p = prev.unwrap();
            eq2prev = fixed(p)[..p.n_number_of_atoms as usize]
                .iter()
                .any(|&x| x != 0);
        }
        if eq2prev {
            mult += 1;
            continue;
        }
        if next != 0 {
            make_delim(Some(S_COMP_DELIM), strbuf, overflow);
        }
        next += 1;
        if let Some(p) = prev {
            let f = fixed(p);
            let not_empty = f[..p.n_number_of_atoms.max(0) as usize]
                .iter()
                .any(|&x| x != 0);
            if not_empty {
                make_mult(mult + 1, "*", strbuf, 0, overflow);
                make_h_string(false, &f, p.n_number_of_atoms, strbuf, atom_mode, overflow);
            } else {
                num_empty += 1;
            }
        }
        prev = cur;
        mult = 0;
    }
    if num_empty == num_components && used(strbuf) > u0 {
        truncate(strbuf, u0);
    }
    used(strbuf) - u0
}

// ---------------------------------------------------------------------------
// Layers that share the "same as the tautomeric layer / same as the previous
// component" skeleton (stereo and isotopic layers)

type EqTaut<'a> = dyn Fn(Option<&InChI>, Option<&InChI>) -> i32 + 'a;
type EqPrev<'a> = dyn Fn(&InChI, &InChI) -> bool + 'a;
type PrintOne<'a> = dyn Fn(&InChI, i32, &mut InchiIosString, &mut i32) + 'a;

struct Layer<'a> {
    /// the equivalence value (`iiSTEREO | iitNONTAUT`, ...) of the current
    /// component to the tautomeric layer, 0 for none
    eq2taut: &'a EqTaut<'a>,
    /// the current component equals the previous one
    eq2prev: &'a EqPrev<'a>,
    /// output of a component (the `MakeMult` and the string)
    print: &'a PrintOne<'a>,
}

#[allow(clippy::too_many_arguments)]
fn run_layer(
    is0: &[InChISort],
    is20_in: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    num_components: usize,
    second_non_taut_pass: bool,
    use_multipliers: bool,
    layer: &Layer,
) -> i32 {
    let is20: &[InChISort] = if second_non_taut_pass { is20_in } else { &[] };
    let mut taut: Option<&InChI> = None;
    let mut prev: Option<&InChI> = None;
    let mut taut_prev: Option<&InChI> = None;
    let (mut mult, mut next) = (0, 0);
    let mut eq2taut_prev = true;
    let mut prev_equ: Option<&str> = None;
    let mut mult_prev_equ = 0;
    let u0 = used(strbuf);
    for i in 0..=num_components {
        let cur = pick(is0, i, num_components, out_type);
        if second_non_taut_pass {
            taut = pick(is20, i, num_components, OUT_T1);
        }
        let eq2taut = (layer.eq2taut)(cur, taut);
        if eq2taut != 0 {
            if let Some(p) = prev.filter(|p| p.n_number_of_atoms != 0) {
                if next != 0 {
                    make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                }
                next += 1;
                (layer.print)(p, mult, strbuf, overflow);
            } else if taut_prev.is_some_and(|t| t.n_number_of_atoms != 0) {
                if next != 0 {
                    make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                }
                next += 1;
            }
            let curr_equ = equ_string(eq2taut);
            if mult_prev_equ != 0 && prev_equ.is_some() {
                if prev_equ == Some(curr_equ) {
                    mult_prev_equ += 1;
                } else {
                    if next != 0 {
                        make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                    }
                    next += 1;
                    make_eq_str(prev_equ, mult_prev_equ, strbuf, overflow);
                    prev_equ = Some(curr_equ);
                    mult_prev_equ = 1;
                }
            } else {
                prev_equ = Some(curr_equ);
                mult_prev_equ = 1;
            }
            prev = None;
            taut_prev = None;
            mult = 0;
            eq2taut_prev = true;
        } else if eq2taut_prev {
            if mult_prev_equ != 0 && prev_equ.is_some() {
                if next != 0 {
                    make_delim(Some(S_COMP_DELIM), strbuf, overflow);
                }
                next += 1;
                make_eq_str(prev_equ, mult_prev_equ, strbuf, overflow);
                prev_equ = None;
                mult_prev_equ = 0;
            }
            eq2taut_prev = false;
            prev = cur;
            taut_prev = taut;
            mult = 0;
        } else {
            let eq2prev = use_multipliers
                && matches!((cur, prev), (Some(c), Some(p)) if (layer.eq2prev)(c, p));
            if eq2prev {
                mult += 1;
                continue;
            }
            if next != 0 {
                make_delim(Some(S_COMP_DELIM), strbuf, overflow);
            }
            next += 1;
            if let Some(p) = prev.filter(|p| p.n_number_of_atoms != 0) {
                (layer.print)(p, mult, strbuf, overflow);
            }
            prev = cur;
            taut_prev = taut;
            mult = 0;
        }
    }
    used(strbuf) - u0
}

use super::ichi::{EQL_EXISTS, EQL_SP2, EQL_SP3, InChIStereo, eql_inchi_stereo};
use super::ichiprt1::{II_EQ2_ISO, II_EQ2_NONTAUT, II_NUMB, II_STEREO, IIT_ISO};
use super::ichiprt2::{
    eql_inchi_isotopic, make_iso_atom_string, make_iso_taut_string, make_stereo_string,
};

fn sp2_of(c: &InChI) -> Option<&InChIStereo> {
    c.stereo.as_ref()
}
fn iso_of(c: &InChI) -> Option<&InChIStereo> {
    c.stereo_isotopic.as_ref()
}

fn print_sp2(
    s: Option<&InChIStereo>,
    mult: i32,
    strbuf: &mut InchiIosString,
    taut_mode: i32,
    ovf: &mut i32,
) {
    if let Some(s) = s.filter(|s| s.n_number_of_stereo_bonds > 0) {
        make_mult(mult + 1, "*", strbuf, 0, ovf);
        make_stereo_string(
            Some(&s.n_bond_atom1),
            Some(&s.n_bond_atom2),
            Some(&s.b_parity),
            s.n_number_of_stereo_bonds,
            strbuf,
            taut_mode,
            ovf,
        );
    }
}

fn print_sp3(
    s: Option<&InChIStereo>,
    mult: i32,
    strbuf: &mut InchiIosString,
    taut_mode: i32,
    ovf: &mut i32,
) {
    if let Some(s) = s.filter(|s| s.n_number_of_stereo_centers > 0) {
        make_mult(mult + 1, "*", strbuf, 0, ovf);
        make_stereo_string(
            Some(&s.n_number),
            None,
            Some(&s.t_parity),
            s.n_number_of_stereo_centers,
            strbuf,
            taut_mode,
            ovf,
        );
    }
}

fn eql(a: Option<&InChIStereo>, b: Option<&InChIStereo>, e: i32) -> bool {
    a.is_some() && b.is_some() && eql_inchi_stereo(a, e, b, e)
}

/// `str_Sp2`: double bond stereo
#[allow(clippy::too_many_arguments)]
pub fn str_sp2(
    is0: &[InChISort],
    is20: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    taut_mode: i32,
    num_components: usize,
    second_non_taut_pass: bool,
    omit_repetitions: bool,
    use_multipliers: bool,
) -> i32 {
    let _ = EQL_EXISTS;
    run_layer(
        is0,
        is20,
        strbuf,
        overflow,
        out_type,
        num_components,
        second_non_taut_pass,
        use_multipliers,
        &Layer {
            eq2taut: &|c, t| {
                if second_non_taut_pass
                    && omit_repetitions
                    && matches!((c, t), (Some(c), Some(t)) if eql(sp2_of(c), sp2_of(t), EQL_SP2))
                {
                    II_STEREO | IIT_NONTAUT
                } else {
                    0
                }
            },
            eq2prev: &|c, p| eql(sp2_of(c), sp2_of(p), EQL_SP2),
            print: &|p, m, sb, ovf| print_sp2(sp2_of(p), m, sb, taut_mode, ovf),
        },
    )
}

/// `str_Sp3`: tetrahedral stereo (`bRelRac` is ignored: `REL_RAC_STEREO_IGN_1_SC` is 0)
#[allow(clippy::too_many_arguments)]
pub fn str_sp3(
    is0: &[InChISort],
    is20: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    taut_mode: i32,
    num_components: usize,
    second_non_taut_pass: bool,
    omit_repetitions: bool,
    use_multipliers: bool,
) -> i32 {
    run_layer(
        is0,
        is20,
        strbuf,
        overflow,
        out_type,
        num_components,
        second_non_taut_pass,
        use_multipliers,
        &Layer {
            eq2taut: &|c, t| {
                if second_non_taut_pass
                    && omit_repetitions
                    && matches!((c, t), (Some(c), Some(t)) if eql(sp2_of(c), sp2_of(t), EQL_SP3))
                {
                    II_STEREO | IIT_NONTAUT
                } else {
                    0
                }
            },
            eq2prev: &|c, p| eql(sp2_of(c), sp2_of(p), EQL_SP3),
            print: &|p, m, sb, ovf| print_sp3(sp2_of(p), m, sb, taut_mode, ovf),
        },
    )
}

/// the isotopic stereo equivalences of `str_IsoSp2/3`
fn iso_eq2taut(c: Option<&InChI>, t: Option<&InChI>, e: i32, second: bool, omit: bool) -> i32 {
    if second && omit {
        if let Some(c) = c
            && eql(iso_of(c), sp2_of(c), e)
        {
            return II_STEREO | IIT_ISO | IIT_NONTAUT | II_EQ2_NONTAUT;
        }
        if let (Some(c), Some(t)) = (c, t) {
            if eql(iso_of(c), sp2_of(t), e) {
                return II_STEREO | IIT_ISO | IIT_NONTAUT;
            }
            if eql(iso_of(c), iso_of(t), e) {
                return II_STEREO | IIT_ISO | IIT_NONTAUT | II_EQ2_ISO;
            }
        }
    } else if !second
        && omit
        && let Some(c) = c
        && eql(iso_of(c), sp2_of(c), e)
    {
        return II_STEREO | IIT_ISO;
    }
    0
}

/// `str_IsoSp2`: isotopic double bond stereo
#[allow(clippy::too_many_arguments)]
pub fn str_iso_sp2(
    is0: &[InChISort],
    is20: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    taut_mode: i32,
    num_components: usize,
    second_non_taut_pass: bool,
    omit_repetitions: bool,
    use_multipliers: bool,
) -> i32 {
    run_layer(
        is0,
        is20,
        strbuf,
        overflow,
        out_type,
        num_components,
        second_non_taut_pass,
        use_multipliers,
        &Layer {
            eq2taut: &|c, t| iso_eq2taut(c, t, EQL_SP2, second_non_taut_pass, omit_repetitions),
            eq2prev: &|c, p| eql(iso_of(c), iso_of(p), EQL_SP2),
            print: &|p, m, sb, ovf| print_sp2(iso_of(p), m, sb, taut_mode, ovf),
        },
    )
}

/// `str_IsoSp3`: isotopic tetrahedral stereo
#[allow(clippy::too_many_arguments)]
pub fn str_iso_sp3(
    is0: &[InChISort],
    is20: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    taut_mode: i32,
    num_components: usize,
    second_non_taut_pass: bool,
    omit_repetitions: bool,
    use_multipliers: bool,
) -> i32 {
    run_layer(
        is0,
        is20,
        strbuf,
        overflow,
        out_type,
        num_components,
        second_non_taut_pass,
        use_multipliers,
        &Layer {
            eq2taut: &|c, t| iso_eq2taut(c, t, EQL_SP3, second_non_taut_pass, omit_repetitions),
            eq2prev: &|c, p| eql(iso_of(c), iso_of(p), EQL_SP3),
            print: &|p, m, sb, ovf| print_sp3(iso_of(p), m, sb, taut_mode, ovf),
        },
    )
}

/// `str_IsoAtoms`: isotopic atoms and t-groups
#[allow(clippy::too_many_arguments)]
pub fn str_iso_atoms(
    is0: &[InChISort],
    is20: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    taut_mode: i32,
    num_components: usize,
    abc_numbers: bool,
    second_non_taut_pass: bool,
    omit_repetitions: bool,
    use_multipliers: bool,
) -> i32 {
    run_layer(
        is0,
        is20,
        strbuf,
        overflow,
        out_type,
        num_components,
        second_non_taut_pass,
        use_multipliers,
        &Layer {
            eq2taut: &|c, t| {
                if second_non_taut_pass && omit_repetitions && eql_inchi_isotopic(c, t) {
                    II_NUMB | IIT_NONTAUT
                } else {
                    0
                }
            },
            eq2prev: &|c, p| eql_inchi_isotopic(Some(c), Some(p)),
            print: &|p, m, sb, ovf| {
                if p.n_number_of_isotopic_atoms > 0 || p.n_number_of_isotopic_t_groups > 0 {
                    make_mult(m + 1, "*", sb, 0, ovf);
                    if p.n_number_of_isotopic_atoms > 0 && *ovf == 0 {
                        make_iso_atom_string(
                            &p.isotopic_atom,
                            p.n_number_of_isotopic_atoms,
                            sb,
                            taut_mode,
                            ovf,
                        );
                    }
                    if p.n_number_of_isotopic_t_groups > 0 && *ovf == 0 {
                        make_delim(Some(if abc_numbers { "," } else { "(" }), sb, ovf);
                        make_iso_taut_string(
                            &p.isotopic_t_group,
                            p.n_number_of_isotopic_t_groups,
                            sb,
                            taut_mode,
                            ovf,
                        );
                        if !abc_numbers {
                            make_delim(Some(")"), sb, ovf);
                        }
                    }
                }
            },
        },
    )
}

fn abs_inv(
    is0: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    num_components: usize,
    sel: fn(&InChI) -> Option<&InChIStereo>,
) -> i32 {
    let u0 = used(strbuf);
    let mut i = 0;
    while *overflow == 0 && i < num_components {
        let cur = {
            let ii = get_ii(out_type, &is0[i]);
            if ii >= 0 {
                is0[i].p_inchi[ii as usize].as_ref()
            } else {
                None
            }
        };
        let j = cur.and_then(sel).map_or(0, |s| s.n_comp_inv2_abs);
        if j != 0 {
            make_delim(Some(if j < 0 { "1" } else { "0" }), strbuf, overflow);
        } else {
            make_delim(Some("."), strbuf, overflow);
        }
        i += 1;
    }
    used(strbuf) - u0
}

/// `str_StereoAbsInv`
pub fn str_stereo_abs_inv(
    is0: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    num_components: usize,
) -> i32 {
    abs_inv(is0, strbuf, overflow, out_type, num_components, sp2_of)
}

/// `str_IsoStereoAbsInv`
pub fn str_iso_stereo_abs_inv(
    is0: &[InChISort],
    strbuf: &mut InchiIosString,
    overflow: &mut i32,
    out_type: i32,
    num_components: usize,
) -> i32 {
    abs_inv(is0, strbuf, overflow, out_type, num_components, iso_of)
}
