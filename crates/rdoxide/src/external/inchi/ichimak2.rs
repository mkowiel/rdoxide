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
        let j1 = canon_ord[c.at_num as usize - 1] as usize;
        if parity_known(p) && p != 3 && at[j1].b_ambiguous_stereo != 0 {
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
