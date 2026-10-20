//! Stereo helpers of InChI (`INCHI_BASE/src/ichister.c`) that work on bond
//! parities of `inp_ATOM`s: `FixUnkn0DStereoBonds`,
//! `ReconcileAllCmlBondParities`, `ReconcileCmlIncidentBondParities` and
//! `get_opposite_sb_atom`.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::incomdef::BOND_TYPE_DOUBLE;
use super::inpdef::{
    InpAtom, MAX_NUM_STEREO_BONDS, SB_PARITY_FLAG, SB_PARITY_MASK, SB_PARITY_SHFT,
};
use super::util::is_el_a_metal;

// extr_ct.h
pub const AB_PARITY_NONE: i32 = 0;
pub const AB_PARITY_ODD: i32 = 1;
pub const AB_PARITY_EVEN: i32 = 2;
pub const AB_PARITY_UNKN: i32 = 3;
pub const AB_PARITY_UNDF: i32 = 4;
pub const AB_PARITY_IISO: i32 = 5;
pub const AB_PARITY_CALC: i32 = 6;
pub const AB_MIN_WELL_DEFINED_PARITY: i32 = 1;
pub const AB_MAX_WELL_DEFINED_PARITY: i32 = 2;

// incomdef.h
pub const STEREO_DBLE_EITHER: i8 = 3;

/// `ATOM_PARITY_WELL_DEF`
pub fn atom_parity_well_def(x: i32) -> bool {
    (AB_MIN_WELL_DEFINED_PARITY..=AB_MAX_WELL_DEFINED_PARITY).contains(&x)
}

/// `FixUnkn0DStereoBonds`: marks every bond with an Unknown stereobond
/// parity as an "either" double bond; returns how many there are
pub fn fix_unkn_0d_stereo_bonds(at: &mut [InpAtom]) -> i32 {
    let mut num = 0;
    for a in at.iter_mut() {
        let mut m = 0;
        while m < MAX_NUM_STEREO_BONDS && a.sb_parity[m] != 0 {
            if i32::from(a.sb_parity[m]) == AB_PARITY_UNKN {
                a.bond_stereo[a.sb_ord[m] as usize] = STEREO_DBLE_EITHER;
                num += 1;
            }
            m += 1;
        }
    }
    num
}

/// `get_opposite_sb_atom`: follows a cumulene from `cur_atom` through the
/// neighbor `icur2nxt` to the atom at the other end of the stereobond.
/// Returns the length of the path (0 if not found) and the opposite atom,
/// the index of its bond to the path and the index of its parity slot.
pub fn get_opposite_sb_atom(
    at: &[InpAtom],
    mut cur_atom: usize,
    mut icur2nxt: usize,
) -> (i32, usize, usize, usize) {
    let mut len = 0;
    while len < 20 {
        len += 1;
        let nxt_atom = at[cur_atom].neighbor[icur2nxt] as usize;
        let mut j = 0;
        while j < MAX_NUM_STEREO_BONDS && at[nxt_atom].sb_parity[j] != 0 {
            if cur_atom == at[nxt_atom].neighbor[at[nxt_atom].sb_ord[j] as usize] as usize {
                return (len, nxt_atom, at[nxt_atom].sb_ord[j] as usize, j);
            }
            j += 1;
        }
        if j != 0 {
            return (0, 0, 0, 0);
        }
        if at[nxt_atom].valence == 2
            && 2 * BOND_TYPE_DOUBLE == i32::from(at[nxt_atom].chem_bonds_valence)
        {
            // follow cumulene =X= path
            icur2nxt = usize::from(at[nxt_atom].neighbor[0] as usize == cur_atom);
            cur_atom = nxt_atom;
        } else {
            return (0, 0, 0, 0);
        }
    }
    (0, 0, 0, 0)
}

/// `ReconcileCmlIncidentBondParities`
pub fn reconcile_cml_incident_bond_parities(
    at: &mut [InpAtom],
    cur_atom: usize,
    prev_atom: i32,
    visited: &mut [i8],
    b_disconnected: bool,
) -> i32 {
    // visited = 0 or parity => atom has not been visited
    //           10 + parity => currently is on the stack + its final parity
    //           20 + parity => has been visited; not on the stack anymore
    if at[cur_atom].valence as usize > MAX_NUM_STEREO_BONDS {
        return 0; // ignore
    }
    if at[cur_atom].sb_parity[0] == 0 {
        return 1; // wrong call
    }
    if visited[cur_atom] >= 10 {
        return 2; // program error
    }
    let mut cur_parity = i32::from(visited[cur_atom] % 10);
    visited[cur_atom] += 10;

    let mut ret = 0;
    let mut i = 0;
    while i < MAX_NUM_STEREO_BONDS && at[cur_atom].sb_parity[i] != 0 {
        let icur2nxt = at[cur_atom].sb_ord[i] as usize;
        let (len, nxt_atom, inxt2cur, j) = get_opposite_sb_atom(at, cur_atom, icur2nxt);
        if len == 0 {
            return 4; // could not find the opposite atom: bond parity data error
        }
        if nxt_atom as i32 == prev_atom {
            i += 1;
            continue;
        }
        if visited[nxt_atom] >= 20 {
            i += 1;
            continue; // back edge, second visit: ignore
        }
        if at[nxt_atom].valence as usize > MAX_NUM_STEREO_BONDS {
            i += 1;
            continue;
        }

        let (cur_sb_parity, b_cur_mask) =
            if b_disconnected && at[cur_atom].sb_parity[i] & SB_PARITY_FLAG != 0 {
                (
                    i32::from(at[cur_atom].sb_parity[i]) >> SB_PARITY_SHFT,
                    3 << SB_PARITY_SHFT,
                )
            } else {
                (i32::from(at[cur_atom].sb_parity[i]) & SB_PARITY_MASK, 3)
            };
        let (mut nxt_sb_parity, b_nxt_mask) =
            if b_disconnected && at[nxt_atom].sb_parity[j] & SB_PARITY_FLAG != 0 {
                (
                    i32::from(at[nxt_atom].sb_parity[j]) >> SB_PARITY_SHFT,
                    3 << SB_PARITY_SHFT,
                )
            } else {
                (i32::from(at[nxt_atom].sb_parity[j]) & SB_PARITY_MASK, 3)
            };

        if !atom_parity_well_def(cur_sb_parity) || !atom_parity_well_def(nxt_sb_parity) {
            if cur_sb_parity == nxt_sb_parity {
                i += 1;
                continue; // bypass unknown/undefined
            }
            return 3; // sb parities do not match: bond parity data error
        }

        let icur2neigh = i32::from(at[cur_atom].sn_ord[i]);
        let inxt2neigh = i32::from(at[nxt_atom].sn_ord[j]);
        let icur2nxt_i = icur2nxt as i32;
        let inxt2cur_i = inxt2cur as i32;
        // parity of the permutation of at[cur_atom].neighbor[] that gives the
        // order { next_atom, neigh_atom, ... } (4 is added because icur2neigh
        // may be negative due to isotopic H removal)
        let cur_order_parity =
            (4 + icur2nxt_i + icur2neigh + i32::from(icur2neigh > icur2nxt_i)) % 2;
        let nxt_order_parity =
            (4 + inxt2cur_i + inxt2neigh + i32::from(inxt2neigh > inxt2cur_i)) % 2;
        let mut nxt_parity = i32::from(visited[nxt_atom] % 10);

        if cur_parity == 0 {
            cur_parity = 2 - (cur_order_parity + cur_sb_parity) % 2;
            visited[cur_atom] += cur_parity as i8;
        } else if cur_parity != 2 - (cur_order_parity + cur_sb_parity) % 2 {
            // Reconcile bond parities: invert both end-atom parities of the bond
            at[cur_atom].sb_parity[i] ^= b_cur_mask as i8;
            at[nxt_atom].sb_parity[j] ^= b_nxt_mask as i8;
            nxt_sb_parity ^= 3;
        }

        if nxt_parity == 0 {
            nxt_parity = 2 - (nxt_order_parity + nxt_sb_parity) % 2;
            visited[nxt_atom] += nxt_parity as i8;
        } else if nxt_parity != 2 - (nxt_order_parity + nxt_sb_parity) % 2 {
            return 5; // algorithm does not work for Mebius-like structures
        }

        // move forward
        if visited[nxt_atom] < 10 {
            ret = reconcile_cml_incident_bond_parities(
                at,
                nxt_atom,
                cur_atom as i32,
                visited,
                b_disconnected,
            );
            if ret != 0 {
                break;
            }
        }
        i += 1;
    }
    visited[cur_atom] += 10; // all bonds of the current atom were processed
    ret
}

/// `ReconcileAllCmlBondParities`
pub fn reconcile_all_cml_bond_parities(at: &mut [InpAtom], b_disconnected: bool) -> i32 {
    let mut visited = vec![0i8; at.len()];
    let mut ret = 0;
    for i in 0..at.len() {
        if at[i].sb_parity[0] != 0
            && visited[i] == 0
            && !(b_disconnected && is_el_a_metal(i32::from(at[i].el_number)))
        {
            ret = reconcile_cml_incident_bond_parities(at, i, -1, &mut visited, b_disconnected);
            if ret != 0 {
                break;
            }
        }
    }
    ret
}
