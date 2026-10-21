//! Filling of the `INChI` output structures, part 1 (`INCHI_BASE/src/ichimak2.c`):
//! Hill formula, stereo copy helpers and ambiguous stereo marking.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichi::{InChI, InChIStereo};
use super::ichimap::compare_lin_ct_stereo_dble;
use super::ichister::{
    AtStereoCarb, AtStereoDble, SpAtom, atom_parity_well_def, bond_chain_len, parity_known,
};
use super::inpdef::{InpAtom, MAX_NUM_STEREO_BONDS};
use super::util::get_element_or_pseudoelement_symbol;

pub const AMBIGUOUS_STEREO_ATOM: i8 = 2;
pub const AMBIGUOUS_STEREO_BOND: i8 = 4;
pub const AMBIGUOUS_STEREO_ATOM_ISO: i8 = 8;
pub const AMBIGUOUS_STEREO_BOND_ISO: i8 = 16;

/// `GetHillFormulaIndexLength`: the number of digits of a count above 1
pub fn get_hill_formula_index_length(count: i32) -> i32 {
    if count > 1 {
        count.to_string().len() as i32
    } else {
        0
    }
}

/// Result of [`get_hill_formula_counts`]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HillCounts {
    pub num_c: i32,
    pub num_h: i32,
    pub len: i32,
    pub num_non_h_atoms: i32,
}

/// `GetHillFormulaCounts`: dimensions of the Hill formula of one component;
/// `None` for an unknown element.
pub fn get_hill_formula_counts(
    n_atom: &[u8],
    n_num_h: &[i8],
    num_atoms: usize,
    n_tautomer: &[u16],
    len_tautomer: i32,
) -> Option<HillCounts> {
    let mut prev_atom: i32 = -2;
    let (mut b_carbon, mut b_hydrogen) = (false, false);
    let (mut elem_len, mut form_len) = (0i32, 0i32);
    let (mut num_h, mut num_c, mut mult) = (0i32, 0i32, 0i32);
    let mut num_non_h = num_atoms as i32;
    let flush = |mult: i32,
                 b_h: bool,
                 b_c: bool,
                 elem_len: i32,
                 num_h: &mut i32,
                 num_c: &mut i32,
                 form_len: &mut i32| {
        if mult != 0 {
            if b_h {
                *num_h += mult;
            } else if b_c {
                *num_c += mult;
            } else {
                *form_len += elem_len + get_hill_formula_index_length(mult);
            }
        }
    };
    for i in 0..num_atoms {
        if prev_atom != i32::from(n_atom[i]) {
            flush(
                mult,
                b_hydrogen,
                b_carbon,
                elem_len,
                &mut num_h,
                &mut num_c,
                &mut form_len,
            );
            let (r, sym) = get_element_or_pseudoelement_symbol(i32::from(n_atom[i]));
            if r == -1 {
                return None;
            }
            mult = 1;
            elem_len = sym.len() as i32;
            prev_atom = i32::from(n_atom[i]);
            b_carbon = sym == "C";
            b_hydrogen = sym == "H";
            if b_hydrogen {
                num_non_h = i as i32;
            }
        } else {
            mult += 1;
        }
        num_h += i32::from(n_num_h[i]);
    }
    if !n_tautomer.is_empty() && len_tautomer > 0 {
        let mut num_groups = i32::from(n_tautomer[0]);
        let mut i = 1usize;
        while (i as i32) < len_tautomer && num_groups > 0 {
            num_h += i32::from(n_tautomer[i + 1]);
            i += n_tautomer[i] as usize + 1;
            num_groups -= 1;
        }
    }
    flush(
        mult,
        b_hydrogen,
        b_carbon,
        elem_len,
        &mut num_h,
        &mut num_c,
        &mut form_len,
    );
    if num_c != 0 {
        form_len += 1 + get_hill_formula_index_length(num_c);
    }
    if num_h != 0 {
        form_len += 1 + get_hill_formula_index_length(num_h);
    }
    Some(HillCounts {
        num_c,
        num_h,
        len: form_len,
        num_non_h_atoms: num_non_h,
    })
}

/// `AddElementAndCount`: appends an element and its count to `out` if it fits
/// into `n_len` characters (including the terminator); returns the length added.
pub fn add_element_and_count(
    element: &str,
    mult: i32,
    out: &mut String,
    n_len: i32,
    b_overflow: &mut i32,
) -> i32 {
    if mult > 0 && *b_overflow == 0 && !element.is_empty() {
        let len1 = element.len() as i32;
        let mult_s = if mult > 1 {
            mult.to_string()
        } else {
            String::new()
        };
        let len2 = mult_s.len() as i32;
        if len1 + len2 < n_len {
            out.push_str(element);
            out.push_str(&mult_s);
            return len1 + len2;
        }
        *b_overflow += 1;
    }
    0
}

/// `MakeHillFormula`: the Hill formula of one component into `out`; the
/// return value is the length, `n_len + 1` on overflow and -1 on error.
pub fn make_hill_formula(
    n_atom: &[u8],
    num_atoms: usize,
    out: &mut String,
    n_len: i32,
    num_c: i32,
    mut num_h: i32,
    b_overflow: &mut i32,
) -> i32 {
    let mut len = 0;
    let mut mult = 0;
    let mut bovfl = 0;
    let mut prev_atom: i32 = -2;
    let mut element = String::new();
    if num_c != 0 {
        len += add_element_and_count("C", num_c, out, n_len - len, &mut bovfl);
        if num_h != 0 {
            len += add_element_and_count("H", num_h, out, n_len - len, &mut bovfl);
            num_h = 0;
        }
    }
    for &a in &n_atom[..num_atoms] {
        if prev_atom != i32::from(a) {
            if mult != 0 {
                len += add_element_and_count(&element, mult, out, n_len - len, &mut bovfl);
            }
            mult = 1;
            let (r, sym) = get_element_or_pseudoelement_symbol(i32::from(a));
            if r == -1 {
                return -1;
            }
            element = sym.to_string();
            prev_atom = i32::from(a);
            if element == "C" {
                return -1;
            }
            let cmp = "H".cmp(element.as_str());
            if cmp == std::cmp::Ordering::Equal {
                return -1;
            }
            if cmp == std::cmp::Ordering::Less && num_h != 0 {
                len += add_element_and_count("H", num_h, out, n_len - len, &mut bovfl);
                num_h = 0;
            }
        } else {
            mult += 1;
        }
    }
    if mult != 0 {
        len += add_element_and_count(&element, mult, out, n_len - len, &mut bovfl);
    }
    if num_h != 0 {
        len += add_element_and_count("H", num_h, out, n_len - len, &mut bovfl);
    }
    *b_overflow |= i32::from(bovfl != 0);
    if bovfl != 0 { n_len + 1 } else { len }
}

/// `AllocateAndFillHillFormula`
pub fn allocate_and_fill_hill_formula(p: &InChI) -> Option<String> {
    let c = get_hill_formula_counts(
        &p.n_atom,
        &p.n_num_h,
        p.n_number_of_atoms as usize,
        &p.n_tautomer,
        p.len_tautomer,
    )?;
    let mut s = String::new();
    let mut ovfl = 0;
    let ret = make_hill_formula(
        &p.n_atom[c.num_c as usize..],
        (c.num_non_h_atoms - c.num_c).max(0) as usize,
        &mut s,
        c.len + 1,
        c.num_c,
        c.num_h,
        &mut ovfl,
    );
    (ret == c.len && ovfl == 0).then_some(s)
}

/// `Copy2StereoBondOrAllene`: 0 => copied to the stereo bonds, 1 => an allene
/// was copied to the stereocenters. `n_bonds` of `None` means the second call
/// (inverted arrays).
#[allow(clippy::too_many_arguments)]
pub fn copy_2_stereo_bond_or_allene(
    stereo: &mut InChIStereo,
    n_centers: &mut i32,
    n_bonds: Option<&mut i32>,
    dble: &AtStereoDble,
    canon_ord: Option<&[u16]>,
    canon_rank: Option<&[u16]>,
    at: &[SpAtom],
    b_isotopic: bool,
) -> i32 {
    if let (Some(canon_ord), Some(canon_rank)) = (canon_ord, canon_rank) {
        let mut j = canon_ord[dble.at_num1 as usize - 1] as usize;
        let mut next_j = 0usize;
        let a = &at[j];
        let (par0, nb1, ord0) = if b_isotopic {
            (
                a.stereo_bond_parity2[0],
                a.stereo_bond_neighbor2[1],
                a.stereo_bond_ord2[0],
            )
        } else {
            (
                a.stereo_bond_parity[0],
                a.stereo_bond_neighbor[1],
                a.stereo_bond_ord[0],
            )
        };
        let mut cumulene_len = bond_chain_len(i32::from(par0));
        if cumulene_len % 2 != 0 && (MAX_NUM_STEREO_BONDS <= 1 || nb1 == 0) {
            next_j = a.neighbor[ord0 as usize] as usize;
            cumulene_len = (cumulene_len - 1) / 2;
            while cumulene_len != 0 && at[next_j].valence == 2 {
                let next_neigh = usize::from(j == at[next_j].neighbor[0] as usize);
                j = next_j;
                next_j = at[next_j].neighbor[next_neigh] as usize;
                cumulene_len -= 1;
            }
        } else {
            cumulene_len = -1;
        }
        if cumulene_len == 0 {
            let at_num = canon_rank[next_j];
            let parity = dble.parity as i8;
            let n = *n_centers as usize;
            let mut j = 0;
            while j < n && stereo.n_number[j] < at_num {
                j += 1;
            }
            let (nn, tp) = if n_bonds.is_some() {
                (&mut stereo.n_number, &mut stereo.t_parity)
            } else {
                (&mut stereo.n_number_inv, &mut stereo.t_parity_inv)
            };
            if j < n {
                nn.copy_within(j..n, j + 1);
                tp.copy_within(j..n, j + 1);
            }
            nn[j] = at_num;
            tp[j] = parity;
            *n_centers += 1;
            return 1;
        }
    }
    if let Some(nb) = n_bonds {
        let j = *nb as usize;
        stereo.b_parity[j] = dble.parity as i8;
        stereo.n_bond_atom1[j] = dble.at_num1;
        stereo.n_bond_atom2[j] = dble.at_num2;
        *nb += 1;
    }
    0
}

/// `CopyLinearCTStereoToINChIStereo`: 0 on success, 1 for no stereo structure,
/// -4/-5 for inconsistent inverted stereo. The stereo arrays must be
/// allocated by the caller.
#[allow(clippy::too_many_arguments)]
pub fn copy_linear_ct_stereo_to_inchi_stereo(
    stereo: Option<&mut InChIStereo>,
    carb: Option<&[AtStereoCarb]>,
    n_len_carb: usize,
    dble: &[AtStereoDble],
    n_len_dble: usize,
    canon_ord: Option<&[u16]>,
    canon_rank: Option<&[u16]>,
    at: &[SpAtom],
    b_isotopic: bool,
    carb_inv: Option<&[AtStereoCarb]>,
    dble_inv: &[AtStereoDble],
    canon_ord_inv: Option<&[u16]>,
    canon_rank_inv: Option<&[u16]>,
) -> i32 {
    let Some(stereo) = stereo else {
        return 1;
    };
    stereo.n_number_of_stereo_centers = n_len_carb as i32;
    for i in 0..n_len_carb {
        if let Some(c) = carb {
            stereo.n_number[i] = c[i].at_num;
            stereo.t_parity[i] = c[i].parity as i8;
        }
        if let Some(c) = carb_inv {
            stereo.n_number_inv[i] = c[i].at_num;
            stereo.t_parity_inv[i] = c[i].parity as i8;
        }
    }
    let mut len_inv = stereo.n_number_of_stereo_centers;
    let mut len = 0;
    for i in 0..n_len_dble {
        let mut nc = stereo.n_number_of_stereo_centers;
        let b_allene = copy_2_stereo_bond_or_allene(
            stereo,
            &mut nc,
            Some(&mut len),
            &dble[i],
            canon_ord,
            canon_rank,
            at,
            b_isotopic,
        );
        stereo.n_number_of_stereo_centers = nc;
        let b_allene_inv = copy_2_stereo_bond_or_allene(
            stereo,
            &mut len_inv,
            None,
            &dble_inv[i],
            canon_ord_inv,
            canon_rank_inv,
            at,
            b_isotopic,
        );
        if b_allene != b_allene_inv
            || (b_allene == 0
                && compare_lin_ct_stereo_dble(Some(&dble[i..]), 1, Some(&dble_inv[i..]), 1) != 0)
        {
            return -4;
        }
    }
    stereo.n_number_of_stereo_bonds = len;
    if len_inv != stereo.n_number_of_stereo_centers {
        return -5;
    }
    let n = stereo.n_number_of_stereo_centers as usize;
    let mut diff = 0;
    for i in 0..n {
        if stereo.n_number_inv[i] != stereo.n_number[i] {
            diff = if stereo.n_number_inv[i] > stereo.n_number[i] {
                2
            } else {
                -2
            };
            break;
        }
        if stereo.t_parity_inv[i] != stereo.t_parity[i] {
            diff = if stereo.t_parity_inv[i] > stereo.t_parity[i] {
                1
            } else {
                -1
            };
            break;
        }
    }
    stereo.n_comp_inv2_abs = (diff > 0) as i32 - (diff < 0) as i32;
    if diff == -1 || diff == 1 {
        let mut diff = 0;
        for i in 0..n {
            if stereo.n_number_inv[i] != stereo.n_number[i] {
                diff = 2;
                break;
            }
            if stereo.t_parity_inv[i] == stereo.t_parity[i]
                || i32::from(stereo.t_parity_inv[i]) + i32::from(stereo.t_parity[i]) != 3
            {
                diff = 1;
                break;
            }
        }
        stereo.b_trivial_inv = i32::from(diff == 0);
    } else {
        stereo.b_trivial_inv = 0;
    }
    0
}

/// `MarkAmbiguousStereo`: marks the ambiguous stereo atoms/bonds in `at` and
/// `norm_at`; returns their number, -1 without the canonical order.
pub fn mark_ambiguous_stereo(
    at: &mut [SpAtom],
    norm_at: &mut [InpAtom],
    b_isotopic: bool,
    canon_ord: Option<&[u16]>,
    carb: &[AtStereoCarb],
    dble: &[AtStereoDble],
) -> i32 {
    let Some(canon_ord) = canon_ord else {
        return -1;
    };
    let mut num = 0;
    let mark_atom = if b_isotopic {
        AMBIGUOUS_STEREO_ATOM_ISO
    } else {
        AMBIGUOUS_STEREO_ATOM
    };
    for c in carb {
        let p = i32::from(c.parity);
        if !(parity_known(p) && p != 3) {
            continue;
        }
        // the C code indexes before the parity test (a stale zero entry reads out of range)
        let j1 = canon_ord[c.at_num as usize - 1] as usize;
        if at[j1].b_ambiguous_stereo != 0 {
            at[j1].b_ambiguous_stereo |= mark_atom;
            norm_at[j1].b_ambiguous_stereo |= mark_atom;
            num += 1;
        }
    }
    let mark_bond = if b_isotopic {
        AMBIGUOUS_STEREO_BOND_ISO
    } else {
        AMBIGUOUS_STEREO_BOND
    };
    for d in dble {
        if !atom_parity_well_def(i32::from(d.parity)) {
            continue;
        }
        let j1 = canon_ord[d.at_num1 as usize - 1] as usize;
        let j2 = canon_ord[d.at_num2 as usize - 1] as usize;
        if at[j1].b_ambiguous_stereo != 0 || at[j2].b_ambiguous_stereo != 0 {
            let j1_parity = if b_isotopic {
                at[j1].stereo_bond_parity2[0]
            } else {
                at[j1].stereo_bond_parity[0]
            };
            let mut cumulene_len = bond_chain_len(i32::from(j1_parity));
            let nb1 = if b_isotopic {
                at[j1].stereo_bond_neighbor2[1]
            } else {
                at[j1].stereo_bond_neighbor[1]
            };
            if cumulene_len % 2 != 0 && (MAX_NUM_STEREO_BONDS <= 1 || nb1 == 0) {
                let mut j = j1;
                let ord = if b_isotopic {
                    at[j].stereo_bond_ord2[0]
                } else {
                    at[j].stereo_bond_ord[0]
                };
                let mut next_j = at[j].neighbor[ord as usize] as usize;
                cumulene_len = (cumulene_len - 1) / 2;
                while cumulene_len != 0 && at[next_j].valence == 2 {
                    let next_neigh = usize::from(j == at[next_j].neighbor[0] as usize);
                    j = next_j;
                    next_j = at[next_j].neighbor[next_neigh] as usize;
                    cumulene_len -= 1;
                }
                if at[next_j].valence == 2 {
                    at[next_j].b_ambiguous_stereo |= mark_atom;
                    norm_at[next_j].b_ambiguous_stereo |= mark_atom;
                    num += 1;
                    continue;
                }
            }
            if at[j1].b_ambiguous_stereo != 0 {
                at[j1].b_ambiguous_stereo |= mark_bond;
                norm_at[j1].b_ambiguous_stereo |= mark_bond;
                num += 1;
            }
            if at[j2].b_ambiguous_stereo != 0 {
                at[j2].b_ambiguous_stereo |= mark_bond;
                norm_at[j2].b_ambiguous_stereo |= mark_bond;
                num += 1;
            }
        }
    }
    num
}

/// `UnmarkAllUndefinedUnknownStereo`: removes all-undefined/unknown stereo
/// layers when requested; returns the `REQ_MODE_*_IGN_ALL_UU` bits applied.
pub fn unmark_all_undefined_unknown_stereo(
    stereo: Option<&mut InChIStereo>,
    n_user_mode: u64,
) -> u64 {
    use super::ichiparm::{REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU};
    let mut ret = 0;
    let Some(s) = stereo else {
        return ret;
    };
    if s.n_number_of_stereo_centers == 0 && s.n_number_of_stereo_bonds == 0 {
        return ret;
    }
    let n = s.n_number_of_stereo_centers;
    if s.n_comp_inv2_abs == 0 && n > 0 && n_user_mode & REQ_MODE_SC_IGN_ALL_UU != 0 {
        let n = n as usize;
        if s.t_parity[..n]
            .iter()
            .all(|&p| !atom_parity_well_def(i32::from(p)))
        {
            s.n_number_of_stereo_centers = 0;
            for i in 0..n {
                s.t_parity[i] = 0;
                s.n_number[i] = 0;
                s.t_parity_inv[i] = 0;
                s.n_number_inv[i] = 0;
            }
            ret |= REQ_MODE_SC_IGN_ALL_UU;
        }
    }
    let n = s.n_number_of_stereo_bonds;
    if n > 0 && n_user_mode & REQ_MODE_SB_IGN_ALL_UU != 0 {
        let n = n as usize;
        if s.b_parity[..n]
            .iter()
            .all(|&p| !atom_parity_well_def(i32::from(p)))
        {
            s.n_number_of_stereo_bonds = 0;
            for i in 0..n {
                s.b_parity[i] = 0;
                s.n_bond_atom1[i] = 0;
                s.n_bond_atom2[i] = 0;
            }
            ret |= REQ_MODE_SB_IGN_ALL_UU;
        }
    }
    ret
}

const WARN_FAILED_STEREO: i32 = 0x0001;
const WARN_FAILED_ISOTOPIC: i32 = 0x0002;
const WARN_FAILED_ISOTOPIC_STEREO: i32 = 0x0004;
const ERR_NO_CANON_RESULTS: i32 = 0x0008;
const INCHI_FLAG_ACID_TAUT: u64 = 0x0001;
const INCHI_FLAG_SC_IGN_ALL_UU: u64 = 0x0008;
const INCHI_FLAG_SB_IGN_ALL_UU: u64 = 0x0010;
const INCHI_FLAG_SC_IGN_ALL_ISO_UU: u64 = 0x0020;
const INCHI_FLAG_HARD_ADD_REM_PROTON: u64 = 0x0080;
/// `INCHI_T_NUM_MOVABLE` (`DISCONNECT_SALTS` is 1, so `REMOVE_TGROUP_CHARGE` is 0)
const INCHI_T_NUM_MOVABLE: usize = 2;
const TG_FLAG_ALL_SALT_DONE: u64 = 0x0000_0008 | 0x0000_0020 | 0x0000_0080;

/// the grouping of equivalent ranks done three times in `FillOutINChI`:
/// `equ[i]` is replaced by the minimal rank of its group, 0 for singletons
fn group_equivalence(equ: &mut [u16], sort_ord: &mut [u16], n: usize) {
    use super::ichimap::comp_ranks_ord;
    use super::ichisort::inchi_qsort;
    {
        let e: &[u16] = equ;
        inchi_qsort(&mut sort_ord[..n], |a, b| comp_ranks_ord(e, *a, *b));
    }
    let mut i = 0usize;
    let mut n_min_ord = sort_ord[0];
    let mut j = 1usize;
    while j <= n {
        if j == n || equ[sort_ord[i] as usize] != equ[sort_ord[j] as usize] {
            n_min_ord = n_min_ord.wrapping_add(1);
            if j - i > 1 {
                while i < j {
                    equ[sort_ord[i] as usize] = n_min_ord;
                    i += 1;
                }
            } else {
                equ[sort_ord[i] as usize] = 0;
                i += 1;
            }
            n_min_ord = sort_ord[j];
        }
        j += 1;
    }
}

fn radical_code(r: i8) -> i8 {
    use super::incomdef::{RADICAL_DOUBLET, RADICAL_SINGLET, RADICAL_TRIPLET};
    match i32::from(r) {
        x if x == RADICAL_SINGLET => 0,
        x if x == RADICAL_DOUBLET => 1,
        x if x == RADICAL_TRIPLET => 2,
        0 => 0,
        _ => 3,
    }
}

/// `FillOutINChI`: serializes one component canonicalization result into
/// `INChI` and `INChI_Aux` (the output arrays must be allocated by the
/// caller). `cs.nCanonOrdStereo`-like arrays are overwritten by the inverted
/// numbering when it is the smaller one, as in the C code.
#[allow(clippy::too_many_arguments)]
pub fn fill_out_inchi(
    p: &mut InChI,
    aux: &mut super::ichi::InChIAux,
    num_atoms: usize,
    num_at_tg: usize,
    num_removed_h: usize,
    at: &mut [SpAtom],
    norm_at: &mut [InpAtom],
    cs: &mut super::ichicano::CanonStat,
    mut tgi: Option<&mut super::ichitaut::TGroupInfo>,
    b_tautomeric: bool,
    n_user_mode: u64,
    str_err: &mut String,
    b_no_warnings: bool,
) -> i32 {
    use super::ichi::{INCHI_FLAG_RAC_STEREO, INCHI_FLAG_REL_STEREO, OrigInfo};
    use super::ichi_bns::FLAG_PROTON_CHARGE_CANCEL;
    use super::ichierr::{CT_OUT_OF_RAM, CT_WRONG_FORMULA, add_error_message};
    use super::ichiparm::{
        REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO, REQ_MODE_SB_IGN_ALL_UU,
        REQ_MODE_SC_IGN_ALL_UU,
    };
    use super::ichitaut::{
        FLAG_FORCE_SALT_TAUT, FLAG_NORM_CONSIDER_TAUT, sort_tautomer_groups_and_endpoints,
    };
    use super::inpdef::MAX_ATOMS;
    use super::runichi3::{TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE};
    use super::util::get_unusual_el_valence;

    let mut n_error_code = 0i32;
    let mut ret = 0i32;
    macro_rules! warn {
        ($m:expr) => {
            if !b_no_warnings {
                add_error_message(str_err, $m);
            }
        };
    }

    'body: {
        // warnings
        if cs.n_len_linear_ct_stereo_carb < 0
            || cs.n_len_linear_ct_stereo_dble < 0
            || cs.n_len_canon_ord_stereo < 0
            || cs.n_len_canon_ord_stereo_taut < 0
        {
            n_error_code |= WARN_FAILED_STEREO;
        }
        if cs.n_len_linear_ct_isotopic < 0
            || cs.n_len_linear_ct_isotopic_tautomer < 0
            || cs.n_len_canon_ord_isotopic < 0
            || cs.n_len_canon_ord_isotopic_taut < 0
        {
            n_error_code |= WARN_FAILED_ISOTOPIC;
        }
        if cs.n_len_linear_ct_isotopic_stereo_carb < 0
            || cs.n_len_linear_ct_isotopic_stereo_dble < 0
            || cs.n_len_canon_ord_isotopic_stereo < 0
            || cs.n_len_canon_ord_isotopic_stereo_taut < 0
        {
            n_error_code |= WARN_FAILED_ISOTOPIC_STEREO;
        }
        let canon_rank_a = vec![0u16; num_at_tg + 1];
        let mut sort_ord = vec![0u16; num_at_tg + 1];

        p.n_total_charge = at[..num_atoms + num_removed_h]
            .iter()
            .map(|a| i32::from(a.charge))
            .sum();
        p.n_number_of_atoms = num_atoms as i32;
        aux.n_number_of_atoms = num_atoms as i32;

        if b_tautomeric {
            if let Some(t) = tgi.as_deref() {
                aux.n_num_removed_protons = t.tni.n_num_removed_protons as i16;
                for i in 0..3 {
                    aux.n_num_removed_isotopic_h[i] = (i32::from(t.num_iso_h[i])
                        + t.tni.n_num_removed_protons_isotopic[i])
                        as i16;
                }
                if aux.b_normalization_flags & FLAG_FORCE_SALT_TAUT != 0 {
                    p.n_flags |= INCHI_FLAG_HARD_ADD_REM_PROTON;
                }
                if aux.b_normalization_flags
                    & (FLAG_NORM_CONSIDER_TAUT & !FLAG_PROTON_CHARGE_CANCEL)
                    != 0
                {
                    warn!("Proton(s) added/removed");
                }
                if aux.b_normalization_flags & FLAG_PROTON_CHARGE_CANCEL != 0 {
                    warn!("Charges neutralized");
                }
            }
        }

        // ------------------------------------------------ non-isotopic numbering
        let nz = |v: i32| v.max(0) as usize;
        let stereo_ok = |c: &super::ichicano::CanonStat,
                         n_carb: i32,
                         n_dble: i32,
                         n_ord: i32,
                         ord: &Vec<u16>,
                         ord_inv: &Vec<u16>,
                         s: [super::ichicano::Slot; 4]| {
            (n_carb > 0 || n_dble > 0)
                && n_ord > 0
                && ((!c.carb_ref(s[0]).is_empty() && !c.carb_ref(s[1]).is_empty())
                    || (!c.dble_ref(s[2]).is_empty() && !c.dble_ref(s[3]).is_empty()))
                && !ord.is_empty()
                && !ord_inv.is_empty()
        };
        let sl = cs.stereo_slots;
        let canon_ord: Option<Vec<u16>>;
        let mut canon_rank = canon_rank_a;
        if stereo_ok(
            cs,
            cs.n_len_linear_ct_stereo_carb,
            cs.n_len_linear_ct_stereo_dble,
            cs.n_len_canon_ord_stereo,
            &cs.n_canon_ord_stereo,
            &cs.n_canon_ord_stereo_inv,
            [sl.carb, sl.carb_inv, sl.dble, sl.dble_inv],
        ) {
            let mut ord = cs.n_canon_ord_stereo.clone();
            let mut ord_inv = cs.n_canon_ord_stereo_inv.clone();
            let mut rank_inv = std::mem::take(&mut sort_ord);
            for i in 0..num_at_tg {
                let v = (i + 1) as u16;
                rank_inv[ord_inv[i] as usize] = v;
                canon_rank[ord[i] as usize] = v;
            }
            let carb = cs.carb_ref(sl.carb);
            let carb_inv = cs.carb_ref(sl.carb_inv);
            let dble = cs.dble_ref(sl.dble);
            let dble_inv = cs.dble_ref(sl.dble_inv);
            let n_err = copy_linear_ct_stereo_to_inchi_stereo(
                p.stereo.as_mut(),
                (!carb.is_empty()).then(|| &carb[..nz(cs.n_len_linear_ct_stereo_carb)]),
                nz(cs.n_len_linear_ct_stereo_carb),
                &dble[..nz(cs.n_len_linear_ct_stereo_dble).min(dble.len())],
                nz(cs.n_len_linear_ct_stereo_dble),
                Some(&ord),
                Some(&canon_rank),
                at,
                false,
                (!carb_inv.is_empty()).then(|| &carb_inv[..nz(cs.n_len_linear_ct_stereo_carb)]),
                dble_inv,
                Some(&ord_inv),
                Some(&rank_inv),
            );
            n_error_code = n_err;
            if n_err == 1 {
                n_error_code = 0;
                ret = CT_OUT_OF_RAM;
                p.n_error_code = CT_OUT_OF_RAM;
                aux.n_error_code = CT_OUT_OF_RAM;
                break 'body;
            }
            let mut use_inv = false;
            {
                let st = p.stereo.as_mut().unwrap();
                if !st.t_parity_inv.is_empty() && !st.n_number_inv.is_empty() {
                    if n_user_mode & REQ_MODE_RELATIVE_STEREO != 0 {
                        p.n_flags |= INCHI_FLAG_REL_STEREO;
                    }
                    if n_user_mode & REQ_MODE_RACEMIC_STEREO != 0 {
                        p.n_flags |= INCHI_FLAG_RAC_STEREO;
                    }
                    if st.n_comp_inv2_abs == -1 {
                        std::mem::swap(&mut st.n_number, &mut st.n_number_inv);
                        std::mem::swap(&mut st.t_parity, &mut st.t_parity_inv);
                        std::mem::swap(&mut canon_rank, &mut rank_inv);
                        std::mem::swap(&mut ord, &mut ord_inv);
                        use_inv = true;
                    }
                }
            }
            for i in 0..num_atoms {
                aux.n_orig_at_nos_in_canon_ord_inv[i] = at[ord_inv[i] as usize].orig_at_number;
                aux.n_orig_at_nos_in_canon_ord[i] = at[ord[i] as usize].orig_at_number;
            }
            if use_inv {
                std::mem::swap(&mut canon_rank, &mut rank_inv);
                std::mem::swap(&mut ord, &mut ord_inv);
                canon_rank[..num_at_tg].copy_from_slice(&rank_inv[..num_at_tg]);
                ord[..num_at_tg].copy_from_slice(&ord_inv[..num_at_tg]);
                cs.n_canon_ord_stereo[..num_at_tg].copy_from_slice(&ord_inv[..num_at_tg]);
            }
            sort_ord = rank_inv;
            canon_ord = Some(ord);
        } else {
            let ord = if cs.n_len_canon_ord_stereo > 0 {
                Some(cs.n_canon_ord_stereo.clone())
            } else if cs.n_len_canon_ord > 0 {
                Some(cs.n_canon_ord.clone())
            } else {
                None
            };
            if let Some(o) = &ord {
                for i in 0..num_atoms {
                    canon_rank[o[i] as usize] = (i + 1) as u16;
                    aux.n_orig_at_nos_in_canon_ord[i] = at[o[i] as usize].orig_at_number;
                }
                for i in num_atoms..num_at_tg {
                    canon_rank[o[i] as usize] = (i + 1) as u16;
                }
            }
            canon_ord = ord;
        }
        let mut canon_ord = canon_ord;

        if !aux.orig_info.is_empty() {
            if let Some(o) = &canon_ord {
                for i in 0..num_atoms {
                    let a = &norm_at[o[i] as usize];
                    if a.valence != 0 || a.num_h != 0 {
                        aux.orig_info[i] = OrigInfo {
                            c_charge: a.charge,
                            c_radical: radical_code(a.radical),
                            c_unusual_valence: get_unusual_el_valence(
                                i32::from(a.el_number),
                                i32::from(a.charge),
                                i32::from(a.radical),
                                i32::from(a.chem_bonds_valence),
                                i32::from(a.num_h),
                                i32::from(a.valence),
                            ) as i8,
                        };
                    } else {
                        aux.orig_info[i].c_radical = radical_code(a.radical);
                    }
                }
            }
        }

        // constitutional equivalence of the atoms
        match &canon_ord {
            Some(o) if !cs.n_symm_rank.is_empty() && !aux.n_constit_equ_numbers.is_empty() => {
                for i in 0..num_atoms {
                    aux.n_constit_equ_numbers[i] = cs.n_symm_rank[o[i] as usize];
                    sort_ord[i] = i as u16;
                }
                for s in sort_ord.iter_mut().take(num_at_tg).skip(num_atoms) {
                    *s = MAX_ATOMS as u16;
                }
                group_equivalence(&mut aux.n_constit_equ_numbers, &mut sort_ord, num_atoms);
            }
            _ => {
                n_error_code |= ERR_NO_CANON_RESULTS;
                ret = -1;
                break 'body;
            }
        }
        {
            let o = canon_ord.as_ref().unwrap();
            for i in 0..num_atoms {
                p.n_atom[i] = at[o[i] as usize].el_number;
            }
        }
        if cs.n_len_linear_ct_at_only <= 0 || cs.linear_ct.is_empty() || p.n_conn_table.is_empty() {
            n_error_code |= ERR_NO_CANON_RESULTS;
            ret = -2;
            break 'body;
        }
        let l = cs.n_len_linear_ct_at_only as usize;
        p.n_conn_table[..l].copy_from_slice(&cs.linear_ct[..l]);
        p.len_conn_table = l as i32;

        // tautomeric groups
        let mut len = 0usize;
        let n_t = if b_tautomeric {
            tgi.as_deref_mut().map_or(0, |t| {
                sort_tautomer_groups_and_endpoints(t, num_atoms, num_at_tg, &canon_rank)
            })
        } else {
            0
        };
        if b_tautomeric && n_t > 0 {
            let t = tgi.as_deref().unwrap();
            if t.b_taut_flags_done & TG_FLAG_ALL_SALT_DONE != 0 {
                p.n_flags |= INCHI_FLAG_ACID_TAUT;
            }
            p.n_tautomer[len] = n_t as u16;
            len += 1;
            for i in 0..n_t as usize {
                let g = t.t_group_number[i] as usize;
                let tg = &t.t_group.as_ref().unwrap()[g];
                p.n_tautomer[len] = tg.n_num_endpoints + INCHI_T_NUM_MOVABLE as u16;
                len += 1;
                for j in 0..INCHI_T_NUM_MOVABLE {
                    p.n_tautomer[len] = tg.num[j];
                    len += 1;
                }
                let first = tg.n_first_endpoint_at_no_pos as usize;
                for j in first..first + tg.n_num_endpoints as usize {
                    p.n_tautomer[len] = canon_rank[t.n_endpoint_atom_number[j] as usize];
                    len += 1;
                }
            }
            p.len_tautomer = len as i32;
            aux.n_number_of_t_groups = n_t;
        } else {
            p.len_tautomer = 0;
            aux.n_number_of_t_groups = 0;
            if let Some(t) = tgi.as_deref() {
                if t.tni.b_normalization_flags & FLAG_NORM_CONSIDER_TAUT != 0
                    || (t.n_num_isotopic_endpoints > 1
                        && t.b_taut_flags_done
                            & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                            != 0)
                {
                    p.len_tautomer = 1;
                    p.n_tautomer[0] = 0;
                }
            }
        }
        if !cs.n_num_h.is_empty() {
            p.n_num_h[..num_atoms].copy_from_slice(&cs.n_num_h[..num_atoms]);
        }
        if !cs.n_num_h_fixed.is_empty() && p.len_tautomer == 0 {
            let f = p.n_num_h_fixed.as_mut().unwrap();
            f[..num_atoms].copy_from_slice(&cs.n_num_h_fixed[..num_atoms]);
            for i in 0..num_atoms {
                p.n_num_h[i] += cs.n_num_h_fixed[i];
            }
        }

        // t-group numbering and symmetry
        let taut_equ = |cs: &super::ichicano::CanonStat,
                        ord_taut: Option<(Vec<u16>, usize)>,
                        symm: &Vec<u16>,
                        equ: &mut Vec<u16>,
                        sort_ord: &mut Vec<u16>| {
            if let Some((o, n)) = ord_taut {
                if !symm.is_empty() && !equ.is_empty() && n > 0 {
                    for i in 0..n {
                        equ[i] = symm[o[i] as usize];
                        sort_ord[i] = i as u16;
                    }
                    group_equivalence(equ, sort_ord, n);
                }
            }
            let _ = cs;
        };
        if p.len_tautomer != 0 {
            let n = aux.n_number_of_t_groups as usize;
            if n != 0 {
                let ord_taut = if cs.n_len_canon_ord_stereo_taut > 0 {
                    Some((cs.n_canon_ord_stereo_taut.clone(), n))
                } else if cs.n_len_canon_ord_taut > 0 {
                    Some((cs.n_canon_ord_taut.clone(), n))
                } else {
                    None
                };
                let symm = cs.n_symm_rank_taut.clone();
                taut_equ(
                    cs,
                    ord_taut,
                    &symm,
                    &mut aux.n_constit_equ_t_group_numbers,
                    &mut sort_ord,
                );
            }
        }

        // Hill formula
        p.sz_hill_formula = allocate_and_fill_hill_formula(p);
        if p.sz_hill_formula.is_none() {
            n_error_code = 0;
            ret = CT_WRONG_FORMULA;
            p.n_error_code = ret;
            aux.n_error_code = ret;
            break 'body;
        }
        let um = unmark_all_undefined_unknown_stereo(p.stereo.as_mut(), n_user_mode);
        if um != 0 {
            if um & REQ_MODE_SC_IGN_ALL_UU != 0 {
                p.n_flags |= INCHI_FLAG_SC_IGN_ALL_UU;
            }
            if um & REQ_MODE_SB_IGN_ALL_UU != 0 {
                p.n_flags |= INCHI_FLAG_SB_IGN_ALL_UU;
            }
            if um & (REQ_MODE_SC_IGN_ALL_UU | REQ_MODE_SB_IGN_ALL_UU) != 0 {
                warn!("Omitted undefined stereo");
            }
        }
        {
            let carb = cs.carb_ref(sl.carb);
            let dble = cs.dble_ref(sl.dble);
            mark_ambiguous_stereo(
                at,
                norm_at,
                false,
                canon_ord.as_deref(),
                &carb[..nz(cs.n_len_linear_ct_stereo_carb).min(carb.len())],
                &dble[..nz(cs.n_len_linear_ct_stereo_dble).min(dble.len())],
            );
        }

        // ------------------------------------------------ isotopic part
        use super::ichicano::Slot;
        if stereo_ok(
            cs,
            cs.n_len_linear_ct_isotopic_stereo_carb,
            cs.n_len_linear_ct_isotopic_stereo_dble,
            cs.n_len_canon_ord_isotopic_stereo,
            &cs.n_canon_ord_isotopic_stereo,
            &cs.n_canon_ord_isotopic_stereo_inv,
            [Slot::Iso, Slot::IsoInv, Slot::Iso, Slot::IsoInv],
        ) {
            let mut ord = cs.n_canon_ord_isotopic_stereo.clone();
            let mut ord_inv = cs.n_canon_ord_isotopic_stereo_inv.clone();
            let mut rank_inv = std::mem::take(&mut sort_ord);
            for i in 0..num_at_tg {
                let v = (i + 1) as u16;
                rank_inv[ord_inv[i] as usize] = v;
                canon_rank[ord[i] as usize] = v;
            }
            let carb = &cs.linear_ct_isotopic_stereo_carb;
            let carb_inv = &cs.linear_ct_isotopic_stereo_carb_inv;
            let dble = &cs.linear_ct_isotopic_stereo_dble;
            let dble_inv = &cs.linear_ct_isotopic_stereo_dble_inv;
            let n_err = copy_linear_ct_stereo_to_inchi_stereo(
                p.stereo_isotopic.as_mut(),
                (!carb.is_empty()).then(|| &carb[..nz(cs.n_len_linear_ct_isotopic_stereo_carb)]),
                nz(cs.n_len_linear_ct_isotopic_stereo_carb),
                &dble[..nz(cs.n_len_linear_ct_isotopic_stereo_dble).min(dble.len())],
                nz(cs.n_len_linear_ct_isotopic_stereo_dble),
                Some(&ord),
                Some(&canon_rank),
                at,
                true,
                (!carb_inv.is_empty())
                    .then(|| &carb_inv[..nz(cs.n_len_linear_ct_isotopic_stereo_carb)]),
                dble_inv,
                Some(&ord_inv),
                Some(&rank_inv),
            );
            n_error_code = n_err;
            if n_err == 1 {
                n_error_code = 0;
                ret = CT_OUT_OF_RAM;
                p.n_error_code = CT_OUT_OF_RAM;
                aux.n_error_code = CT_OUT_OF_RAM;
                break 'body;
            }
            let mut use_inv = false;
            {
                let st = p.stereo_isotopic.as_mut().unwrap();
                if !st.t_parity_inv.is_empty() && !st.n_number_inv.is_empty() {
                    if n_user_mode & REQ_MODE_RELATIVE_STEREO != 0 {
                        p.n_flags |= INCHI_FLAG_REL_STEREO;
                    }
                    if n_user_mode & REQ_MODE_RACEMIC_STEREO != 0 {
                        p.n_flags |= INCHI_FLAG_RAC_STEREO;
                    }
                    if st.n_comp_inv2_abs == -1 {
                        std::mem::swap(&mut st.n_number, &mut st.n_number_inv);
                        std::mem::swap(&mut st.t_parity, &mut st.t_parity_inv);
                        std::mem::swap(&mut canon_rank, &mut rank_inv);
                        std::mem::swap(&mut ord, &mut ord_inv);
                        use_inv = true;
                    }
                }
            }
            for i in 0..num_atoms {
                aux.n_isotopic_orig_at_nos_in_canon_ord_inv[i] =
                    at[ord_inv[i] as usize].orig_at_number;
                aux.n_isotopic_orig_at_nos_in_canon_ord[i] = at[ord[i] as usize].orig_at_number;
            }
            if use_inv {
                std::mem::swap(&mut canon_rank, &mut rank_inv);
                std::mem::swap(&mut ord, &mut ord_inv);
                canon_rank[..num_at_tg].copy_from_slice(&rank_inv[..num_at_tg]);
                ord[..num_at_tg].copy_from_slice(&ord_inv[..num_at_tg]);
                cs.n_canon_ord_isotopic_stereo[..num_at_tg].copy_from_slice(&ord_inv[..num_at_tg]);
            }
            sort_ord = rank_inv;
            canon_ord = Some(ord);
        } else {
            let ord = if cs.n_len_canon_ord_isotopic_stereo > 0 {
                Some(cs.n_canon_ord_isotopic_stereo.clone())
            } else if cs.n_len_canon_ord_isotopic > 0 {
                Some(cs.n_canon_ord_isotopic.clone())
            } else {
                None
            };
            if let Some(o) = &ord {
                if !aux.n_isotopic_orig_at_nos_in_canon_ord.is_empty() {
                    for i in 0..num_atoms {
                        canon_rank[o[i] as usize] = (i + 1) as u16;
                        aux.n_isotopic_orig_at_nos_in_canon_ord[i] =
                            at[o[i] as usize].orig_at_number;
                    }
                    for i in num_atoms..num_at_tg {
                        canon_rank[o[i] as usize] = (i + 1) as u16;
                    }
                }
            }
            canon_ord = ord;
        }
        match &canon_ord {
            Some(o)
                if !aux.n_constit_equ_isotopic_numbers.is_empty()
                    && !cs.n_symm_rank_isotopic.is_empty() =>
            {
                for i in 0..num_atoms {
                    aux.n_constit_equ_isotopic_numbers[i] = cs.n_symm_rank_isotopic[o[i] as usize];
                    sort_ord[i] = i as u16;
                }
                for (i, s) in sort_ord
                    .iter_mut()
                    .enumerate()
                    .take(num_at_tg)
                    .skip(num_atoms)
                {
                    *s = i as u16;
                }
                group_equivalence(
                    &mut aux.n_constit_equ_isotopic_numbers,
                    &mut sort_ord,
                    num_atoms,
                );
            }
            _ => break 'body,
        }

        let n = cs.n_len_linear_ct_isotopic;
        p.n_number_of_isotopic_atoms = n;
        for i in 0..n.max(0) as usize {
            let e = &cs.linear_ct_isotopic[i];
            p.isotopic_atom[i] = super::ichi::InChIIsotopicAtom {
                n_atom_number: e.at_num,
                n_iso_difference: e.iso_atw_diff,
                n_num_h: e.num_1h,
                n_num_d: e.num_d,
                n_num_t: e.num_t,
            };
        }
        let n = cs.n_len_linear_ct_isotopic_tautomer;
        p.n_number_of_isotopic_t_groups = n;
        for i in 0..n.max(0) as usize {
            let e = &cs.linear_ct_isotopic_tautomer[i];
            p.isotopic_t_group[i] = super::ichi::InChIIsotopicTGroup {
                n_t_group_number: e.tgroup_num,
                n_num_h: e.num[2],
                n_num_d: e.num[1],
                n_num_t: e.num[0],
            };
        }
        if !cs.n_exchg_iso_h.is_empty() && !p.n_possible_locations_of_isotopic_h.is_empty() {
            let mut j = 1;
            for i in 0..num_atoms {
                if cs.n_exchg_iso_h[i] != 0 {
                    p.n_possible_locations_of_isotopic_h[j] = (i + 1) as u16;
                    j += 1;
                }
            }
            p.n_possible_locations_of_isotopic_h[0] = j as u16;
        }
        let um = unmark_all_undefined_unknown_stereo(p.stereo_isotopic.as_mut(), n_user_mode);
        if um != 0 {
            if um & REQ_MODE_SC_IGN_ALL_UU != 0 {
                p.n_flags |= INCHI_FLAG_SC_IGN_ALL_ISO_UU;
            }
            if um & REQ_MODE_SB_IGN_ALL_UU != 0 {
                // as in the C code: the stereo bond flag is the stereocenter one
                p.n_flags |= INCHI_FLAG_SC_IGN_ALL_ISO_UU;
            }
            if um & (REQ_MODE_SC_IGN_ALL_UU | REQ_MODE_SB_IGN_ALL_UU) != 0 {
                warn!("Omitted undefined stereo");
            }
        }
        mark_ambiguous_stereo(
            at,
            norm_at,
            true,
            canon_ord.as_deref(),
            &cs.linear_ct_isotopic_stereo_carb[..nz(cs.n_len_linear_ct_isotopic_stereo_carb)
                .min(cs.linear_ct_isotopic_stereo_carb.len())],
            &cs.linear_ct_isotopic_stereo_dble[..nz(cs.n_len_linear_ct_isotopic_stereo_dble)
                .min(cs.linear_ct_isotopic_stereo_dble.len())],
        );

        // isotopic t-groups numbering and symmetry
        if p.len_tautomer != 0
            && !aux.n_constit_equ_isotopic_t_group_numbers.is_empty()
            && !cs.n_symm_rank_isotopic_taut.is_empty()
            && (cs.n_len_linear_ct_isotopic != 0 || cs.n_len_linear_ct_isotopic_tautomer != 0)
            && tgi.as_deref().is_some_and(|t| t.num_t_groups > 0)
        {
            let ord_taut = if cs.n_len_canon_ord_isotopic_stereo_taut > 0 {
                Some((
                    cs.n_canon_ord_isotopic_stereo_taut.clone(),
                    cs.n_len_canon_ord_isotopic_stereo_taut as usize,
                ))
            } else if cs.n_len_canon_ord_isotopic_taut > 0 {
                Some((
                    cs.n_canon_ord_isotopic_taut.clone(),
                    cs.n_len_canon_ord_isotopic_taut as usize,
                ))
            } else {
                None
            };
            let symm = cs.n_symm_rank_isotopic_taut.clone();
            taut_equ(
                cs,
                ord_taut,
                &symm,
                &mut aux.n_constit_equ_isotopic_t_group_numbers,
                &mut sort_ord,
            );
        }
    }
    p.n_error_code |= n_error_code;
    aux.n_error_code |= n_error_code;
    ret
}
