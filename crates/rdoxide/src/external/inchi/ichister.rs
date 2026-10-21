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

// ---- stereo geometry helpers ----

use super::ichisort::insertions_sort;

pub const ZTYPE_DOWN: i32 = -1;
pub const ZTYPE_NONE: i32 = 0;
pub const ZTYPE_UP: i32 = 1;
pub const ZTYPE_3D: i32 = 3;
pub const ZTYPE_EITHER: i32 = 9999;

pub const STEREO_SNGL_UP: i32 = 1;
pub const STEREO_SNGL_EITHER: i32 = 4;
pub const STEREO_SNGL_DOWN: i32 = 6;
pub const AMBIGUOUS_STEREO: i32 = 1;

pub const MAX_NUM_STEREO_ATOM_NEIGH: usize = 4;
pub const MIN_ANGLE: f64 = 0.10;
pub const MIN_SINE: f64 = 0.03;
pub const MIN_SINE_OUTSIDE: f64 = 0.06;
pub const MIN_SINE_SQUARE: f64 = 0.125;
pub const MIN_SINE_EDGE: f64 = 0.167;
pub const MIN_LEN_STRAIGHT: f64 = 1.900;
pub const MIN_BOND_LEN: f64 = 0.000001;
pub const ZERO_LENGTH: f64 = MIN_BOND_LEN;
/// `STEREO_CENTER_BONDS_NORM == 1`
pub const MPY_SINE: f64 = 1.00;
pub const MAX_EDGE_RATIO: f64 = 2.50;

pub const T2D_OKAY: i32 = 1;
pub const T2D_WARN: i32 = 2;
pub const T2D_UNDF: i32 = 4;

/// `get_z_coord`: the z difference between the neighbor and the atom and the
/// kind of the stereo information of the bond (`ZTYPE_*`)
pub fn get_z_coord(
    at: &[InpAtom],
    cur_atom: usize,
    neigh_no: usize,
    b_pointed_edge_stereo: i32,
) -> (f64, i32) {
    let stereo_value = i32::from(at[cur_atom].bond_stereo[neigh_no]);
    let stereo_type = stereo_value.abs();
    let neigh = at[cur_atom].neighbor[neigh_no] as usize;
    let z = at[neigh].z - at[cur_atom].z;
    let mut b_flat = z.abs() < ZERO_LENGTH;
    if b_flat {
        for i in 0..at[cur_atom].valence as usize {
            if (at[cur_atom].z - at[at[cur_atom].neighbor[i] as usize].z).abs() > ZERO_LENGTH {
                b_flat = false;
                break;
            }
        }
    }
    let n_type;
    if b_flat {
        if b_pointed_edge_stereo == 0 || b_pointed_edge_stereo * stereo_value >= 0 {
            let mut t = match stereo_type {
                0 => ZTYPE_NONE,
                STEREO_SNGL_UP => ZTYPE_UP,
                STEREO_SNGL_EITHER => ZTYPE_EITHER,
                STEREO_SNGL_DOWN => ZTYPE_DOWN,
                _ => ZTYPE_NONE,
            };
            if stereo_value < 0 && (t == ZTYPE_DOWN || t == ZTYPE_UP) {
                t = -t;
            }
            n_type = t;
        } else {
            n_type = ZTYPE_NONE;
        }
    } else if stereo_type == STEREO_SNGL_EITHER
        && (b_pointed_edge_stereo == 0 || b_pointed_edge_stereo * stereo_value >= 0)
    {
        n_type = ZTYPE_EITHER;
    } else {
        n_type = ZTYPE_3D;
    }
    (z, n_type)
}

pub fn len3(c: &[f64; 3]) -> f64 {
    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}
pub fn len2(c: &[f64; 2]) -> f64 {
    (c[0] * c[0] + c[1] * c[1]).sqrt()
}
pub fn diff3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn add3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn mult3(a: &[f64; 3], b: f64) -> [f64; 3] {
    [a[0] * b, a[1] * b, a[2] * b]
}
pub fn change_sign3(a: &[f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
pub fn dot_prod3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// `dot_prodchar3`
pub fn dot_prodchar3(a: &[i8; 3], b: &[i8; 3]) -> i32 {
    let prod = (i32::from(a[0]) * i32::from(b[0])
        + i32::from(a[1]) * i32::from(b[1])
        + i32::from(a[2]) * i32::from(b[2]))
        / 100;
    prod.clamp(-100, 100)
}
pub fn cross_prod3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        -(a[0] * b[2] - a[2] * b[0]),
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `triple_prod`: the triple product and, when asked, the sine of the angle
/// between `c` and the plane of `a` and `b`
pub fn triple_prod(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3], sine_value: Option<&mut f64>) -> f64 {
    let ab = cross_prod3(a, b);
    let dot_prod_ab_c = dot_prod3(&ab, c);
    if let Some(sv) = sine_value {
        let abs_c = len3(c);
        let abs_ab = len3(&ab);
        *sv = if abs_c > 1.0e-7 && abs_ab > 1.0e-7 {
            MPY_SINE * dot_prod_ab_c / (abs_c * abs_ab)
        } else {
            0.0
        };
    }
    dot_prod_ab_c
}

/// `Get2DTetrahedralAmbiguity`: is the 2D picture of a stereocenter with 3 or
/// 4 bonds (and `at_coord` relative coordinates) unambiguous? Returns
/// `T2D_*` flags or -1 on an error.
pub fn get_2d_tetrahedral_ambiguity(
    at_coord: &[[f64; 3]],
    b_add_explicit_neighbor: bool,
    b_fix_2d_stereo_border_case: bool,
    v_min_angle: f64,
) -> i32 {
    let one_pi = std::f64::consts::PI;
    let two_pi = 2.0 * one_pi;
    let d_angle_and_pi_max_diff = 2.0 * 1.0f64.atan2(7.0f64.sqrt());
    let n_num_neigh = MAX_NUM_STEREO_ATOM_NEIGH - usize::from(b_add_explicit_neighbor);
    let mut n_bond_type = [0i32; MAX_NUM_STEREO_ATOM_NEIGH];
    let mut n_bond_order = [0usize; MAX_NUM_STEREO_ATOM_NEIGH];
    let mut d_bond_direction = [0f64; MAX_NUM_STEREO_ATOM_NEIGH];
    let (mut num_up, mut num_dn) = (0i32, 0i32);
    for i in 0..n_num_neigh {
        let mut d_angle = at_coord[i][1].atan2(at_coord[i][0]);
        if d_angle < 0.0 {
            d_angle += two_pi;
        }
        d_bond_direction[i] = d_angle;
        n_bond_type[i] = if at_coord[i][2] > 0.0 {
            1
        } else if at_coord[i][2] < 0.0 {
            -1
        } else {
            0
        };
        if n_bond_type[i] > 0 {
            num_up += 1;
        } else if n_bond_type[i] < 0 {
            num_dn += 1;
        }
        n_bond_order[i] = i;
    }
    if num_up < num_dn {
        for t in n_bond_type.iter_mut().take(n_num_neigh) {
            *t = -*t;
        }
        std::mem::swap(&mut num_dn, &mut num_up);
    }
    if num_up == 0 {
        return T2D_UNDF;
    }
    let dir = d_bond_direction;
    insertions_sort(&mut n_bond_order[..n_num_neigh], |a, b| {
        let diff = dir[*a] - dir[*b];
        if diff > 0.0 {
            1
        } else if diff < 0.0 {
            -1
        } else {
            0
        }
    });
    let nn = n_num_neigh;
    let (mut cur_len_up, mut len_up, mut first_up);
    if num_up as usize == nn {
        len_up = nn as i32;
        cur_len_up = nn as i32;
        first_up = 0usize;
    } else {
        cur_len_up = 0;
        len_up = 0;
        first_up = 0;
        let mut b_prev_up = false;
        let mut cur_first_up = 0usize;
        let mut i = 0usize;
        while i < nn || b_prev_up {
            if n_bond_type[n_bond_order[i % nn]] > 0 {
                if b_prev_up {
                    cur_len_up += 1;
                } else {
                    b_prev_up = true;
                    cur_len_up = 1;
                    cur_first_up = i % nn;
                }
            } else if b_prev_up {
                if cur_len_up > len_up {
                    first_up = cur_first_up;
                    len_up = cur_len_up;
                }
                b_prev_up = false;
            }
            i += 1;
        }
    }
    // FIX_2D_STEREO_BORDER_CASE == 0: controlled by the runtime flag
    const ZERO_ANGLE: f64 = 0.000001;
    if b_fix_2d_stereo_border_case {
        if nn as i32 - len_up >= 2 {
            for i in 1..(nn as i32 - len_up) as usize {
                let next2 = (first_up + len_up as usize + i) % nn;
                if n_bond_type[n_bond_order[next2]] > 0 {
                    let next1 = (first_up + len_up as usize) % nn;
                    let d_angle = d_bond_direction[n_bond_order[next1]]
                        - d_bond_direction[n_bond_order[next2]];
                    if d_angle.abs() < ZERO_ANGLE {
                        n_bond_order.swap(next1, next2);
                        len_up += 1;
                        break;
                    }
                }
            }
        }
        if nn as i32 - len_up >= 2 {
            for i in 1..(nn as i32 - len_up) as usize {
                let next2 = (first_up + nn - i - 1) % nn;
                if n_bond_type[n_bond_order[next2]] > 0 {
                    let next1 = (first_up + nn - 1) % nn;
                    let d_angle = d_bond_direction[n_bond_order[next1]]
                        - d_bond_direction[n_bond_order[next2]];
                    if d_angle.abs() < ZERO_ANGLE {
                        n_bond_order.swap(next1, next2);
                        first_up = next1;
                        len_up += 1;
                        break;
                    }
                }
            }
        }
    }
    // turn all the bonds so that the 1st Up bond has zero direction
    let d_alpha0 = d_bond_direction[n_bond_order[first_up]];
    for i in 0..nn {
        if i == n_bond_order[first_up] {
            d_bond_direction[i] = 0.0;
        } else {
            let mut d_angle = d_bond_direction[i] - d_alpha0;
            if d_angle < 0.0 {
                d_angle += two_pi;
            }
            d_bond_direction[i] = d_angle;
        }
    }
    let dd = |k: usize| d_bond_direction[n_bond_order[(first_up + k) % nn]];
    let mut ret: i32;
    if nn == 3 {
        match num_up {
            0 => return T2D_UNDF,
            1 => {
                if num_dn != 0 {
                    ret = T2D_UNDF | T2D_WARN;
                } else {
                    let mut d_angle = dd(2) - dd(1);
                    if d_angle < 0.0 {
                        d_angle += two_pi;
                    }
                    ret = if d_angle - one_pi < -v_min_angle || d_angle - one_pi > v_min_angle {
                        T2D_OKAY
                    } else {
                        T2D_UNDF | T2D_WARN
                    };
                }
            }
            2 => {
                if num_dn != 0 {
                    let mut d_alpha = dd(1) - dd(0);
                    if d_alpha < 0.0 {
                        d_alpha += two_pi;
                    }
                    if d_alpha > one_pi - v_min_angle {
                        ret = T2D_OKAY;
                    } else if d_alpha < two_pi / 3.0 - v_min_angle {
                        ret = T2D_UNDF | T2D_WARN;
                    } else {
                        let mut d_bisector = dd(0);
                        d_bisector += dd(1);
                        d_bisector /= 2.0;
                        d_bisector -= one_pi;
                        if d_bisector < 0.0 {
                            d_bisector += two_pi;
                        }
                        let d_limit = if d_alpha < two_pi / 3.0 + v_min_angle {
                            v_min_angle * 3.0 / 2.0
                        } else {
                            d_alpha * 3.0 / 2.0 - one_pi
                        };
                        let d_angle = dd(2);
                        ret = if d_bisector - d_angle < -d_limit || d_bisector - d_angle > d_limit {
                            T2D_UNDF | T2D_WARN
                        } else {
                            T2D_OKAY
                        };
                    }
                } else {
                    ret = T2D_OKAY;
                }
            }
            3 => ret = T2D_OKAY,
            _ => return -1,
        }
    } else if nn == 4 {
        match num_up {
            0 => return T2D_UNDF,
            1 => {
                if num_dn != 0 {
                    if n_bond_type[n_bond_order[(first_up + 2) % nn]] < 0 {
                        let mut d_angle = dd(3) - dd(1);
                        if d_angle < 0.0 {
                            d_angle += two_pi;
                        }
                        ret = if (d_angle - one_pi).abs() < d_angle_and_pi_max_diff + v_min_angle {
                            T2D_UNDF | T2D_WARN
                        } else {
                            T2D_OKAY
                        };
                    } else {
                        ret = T2D_OKAY;
                    }
                } else {
                    ret = T2D_OKAY;
                    let mut d_angle = dd(3) - dd(1);
                    if d_angle < 0.0 {
                        d_angle += two_pi;
                    }
                    if d_angle < one_pi - v_min_angle {
                        ret |= T2D_WARN;
                    }
                }
            }
            2 => {
                if b_fix_2d_stereo_border_case {
                    if len_up == 1 {
                        ret = T2D_OKAY;
                    } else {
                        let d_angle = (two_pi - (dd(3) - dd(0))).abs();
                        let d_alpha = (dd(2) - dd(1)).abs();
                        ret = if (d_angle < 2.0 * ZERO_ANGLE && d_alpha > v_min_angle)
                            || (d_alpha < 2.0 * ZERO_ANGLE && d_angle > v_min_angle)
                        {
                            T2D_OKAY | T2D_WARN
                        } else {
                            T2D_UNDF | T2D_WARN
                        };
                    }
                } else {
                    // the original InChI v. 1 bug
                    ret = if cur_len_up == 1 {
                        T2D_OKAY
                    } else {
                        T2D_UNDF | T2D_WARN
                    };
                }
            }
            3 => {
                ret = T2D_OKAY;
                let mut d_angle = dd(2) - dd(0);
                if d_angle < 0.0 {
                    d_angle += two_pi;
                }
                if d_angle < one_pi - v_min_angle {
                    ret |= T2D_WARN;
                }
            }
            4 => ret = T2D_UNDF | T2D_WARN,
            _ => return -1,
        }
        if ret == T2D_OKAY {
            for i in 0..nn {
                let mut d_angle = d_bond_direction[n_bond_order[(i + nn - 1) % nn]]
                    - d_bond_direction[n_bond_order[i % nn]];
                if d_angle < 0.0 {
                    d_angle += two_pi;
                }
                if d_angle < one_pi - v_min_angle {
                    ret |= T2D_WARN;
                    break;
                }
            }
        }
    } else {
        return -1;
    }
    let _ = &mut cur_len_up;
    ret
}

/// `triple_prod_and_min_abs_sine2`
pub fn triple_prod_and_min_abs_sine2(
    at_coord: &[[f64; 3]],
    central_at_coord: &[f64; 3],
    b_added_explicit_neighbor: bool,
    min_sine: Option<&mut f64>,
    mut b_ambiguous: Option<&mut i32>,
    v_min_sine: f64,
) -> f64 {
    let Some(min_sine) = min_sine else {
        return triple_prod(&at_coord[0], &at_coord[1], &at_coord[2], None);
    };
    let mut min_sine_value = 9999.0f64;
    let mut sine_value = 0.0;
    let ret = triple_prod(
        &at_coord[0],
        &at_coord[1],
        &at_coord[2],
        Some(&mut sine_value),
    );
    sine_value = MPY_SINE * sine_value.abs();
    let e2 = diff3(&at_coord[1], &at_coord[0]);
    let e1 = diff3(&at_coord[0], &at_coord[2]);
    let e0 = diff3(&at_coord[2], &at_coord[1]);
    let e03 = len3(&at_coord[0]);
    let e13 = len3(&at_coord[1]);
    let e23 = len3(&at_coord[2]);
    let e02 = len3(&e1);
    let e12 = len3(&e0);
    let e01 = len3(&e2);
    let mut max_edge_len = e03;
    let mut min_edge_len = e03;
    if min_edge_len > e13 {
        min_edge_len = e13;
    }
    if min_edge_len > e01 {
        min_edge_len = e01;
    }
    let mut min_edge_len_no_expl = min_edge_len;
    if min_edge_len > e23 {
        min_edge_len = e23;
    }
    if min_edge_len > e02 {
        min_edge_len = e02;
    }
    if min_edge_len > e12 {
        min_edge_len = e12;
    }
    if max_edge_len < e13 {
        max_edge_len = e13;
    }
    if max_edge_len < e01 {
        max_edge_len = e01;
    }
    let mut max_edge_len_no_expl = max_edge_len;
    if max_edge_len < e23 {
        max_edge_len = e23;
    }
    if max_edge_len < e02 {
        max_edge_len = e02;
    }
    if max_edge_len < e12 {
        max_edge_len = e12;
    }
    if !b_added_explicit_neighbor {
        min_edge_len_no_expl = min_edge_len;
        max_edge_len_no_expl = max_edge_len;
    }
    let b_long_edges = if b_added_explicit_neighbor {
        max_edge_len_no_expl < MAX_EDGE_RATIO * min_edge_len_no_expl
    } else {
        max_edge_len < MAX_EDGE_RATIO * min_edge_len
    };
    if sine_value > v_min_sine {
        let prod = ret.abs();
        let s0 = prod / len3(&cross_prod3(&at_coord[1], &at_coord[2]));
        let s1 = prod / len3(&cross_prod3(&at_coord[0], &at_coord[2]));
        let s2 = prod / len3(&cross_prod3(&at_coord[0], &at_coord[1]));
        let s3 = prod / len3(&cross_prod3(&e0, &e1));
        for (s, e) in [
            (s0, e01),
            (s0, e02),
            (s0, e03),
            (s1, e01),
            (s1, e12),
            (s1, e13),
            (s2, e02),
            (s2, e12),
            (s2, e23),
            (s3, e03),
            (s3, e13),
            (s3, e23),
        ] {
            sine_value = s / e;
            if sine_value < min_sine_value {
                min_sine_value = sine_value;
            }
        }
        sine_value = MPY_SINE * min_sine_value;
        *min_sine = sine_value;
        if let Some(amb) = b_ambiguous.as_deref_mut() {
            if sine_value >= v_min_sine {
                let tmp = diff3(central_at_coord, &at_coord[0]);
                let central_prod = [
                    triple_prod(&at_coord[0], &at_coord[1], central_at_coord, None),
                    triple_prod(&at_coord[1], &at_coord[2], central_at_coord, None),
                    triple_prod(&at_coord[2], &at_coord[0], central_at_coord, None),
                    triple_prod(&e2, &e1, &tmp, None),
                ];
                for cp in central_prod {
                    if cp / ret < -MIN_SINE_OUTSIDE {
                        *amb |= AMBIGUOUS_STEREO;
                        break;
                    }
                }
            }
        }
        // STEREO_CENTER_BONDS_NORM == 1
        if b_long_edges
            && !b_added_explicit_neighbor
            && max_edge_len >= MIN_LEN_STRAIGHT
            && sine_value < MIN_SINE_SQUARE
        {
            *min_sine = v_min_sine / 2.0;
            if let Some(amb) = b_ambiguous.as_deref_mut() {
                if *amb == 0 {
                    *amb |= AMBIGUOUS_STEREO;
                }
            }
        }
        if b_long_edges
            && sine_value < MIN_SINE_SQUARE
            && sine_value < MIN_SINE_EDGE * min_edge_len_no_expl
        {
            *min_sine = v_min_sine / 2.0;
            if let Some(amb) = b_ambiguous {
                if *amb == 0 {
                    *amb |= AMBIGUOUS_STEREO;
                }
            }
        }
    } else {
        *min_sine = sine_value;
    }
    ret
}

/// `triple_prod_and_min_abs_sine`
pub fn triple_prod_and_min_abs_sine(at_coord: &[[f64; 3]], min_sine: Option<&mut f64>) -> f64 {
    let Some(min_sine) = min_sine else {
        return triple_prod(&at_coord[0], &at_coord[1], &at_coord[2], None);
    };
    let mut min_sine_value = 9999.0f64;
    let mut sine_value = 0.0;
    triple_prod(
        &at_coord[0],
        &at_coord[1],
        &at_coord[2],
        Some(&mut sine_value),
    );
    min_sine_value = min_sine_value.min(sine_value.abs());
    triple_prod(
        &at_coord[1],
        &at_coord[2],
        &at_coord[0],
        Some(&mut sine_value),
    );
    min_sine_value = min_sine_value.min(sine_value.abs());
    let prod = triple_prod(
        &at_coord[2],
        &at_coord[0],
        &at_coord[1],
        Some(&mut sine_value),
    );
    min_sine_value = min_sine_value.min(sine_value.abs());
    *min_sine = min_sine_value;
    prod
}

/// `are_3_vect_in_one_plane`
pub fn are_3_vect_in_one_plane(at_coord: &[[f64; 3]], min_sine: f64) -> bool {
    let mut actual = 0.0;
    triple_prod_and_min_abs_sine(at_coord, Some(&mut actual));
    actual <= min_sine
}

/// `are_4at_in_one_plane`
pub fn are_4at_in_one_plane(at_coord: &[[f64; 3]], min_sine: f64) -> bool {
    let mut min_actual = 0.0;
    for k in 0..4 {
        let mut coord = [[0.0; 3]; 3];
        let mut j = 0;
        for i in 0..4 {
            if i != k {
                coord[j] = diff3(&at_coord[i], &at_coord[k]);
                j += 1;
            }
        }
        let mut actual = 0.0;
        triple_prod_and_min_abs_sine(&coord, Some(&mut actual));
        if k == 0 || actual < min_actual {
            min_actual = actual;
        }
    }
    min_actual <= min_sine
}

/// `triple_prod_char`: the signed triple product (times 100) of the stereo
/// bond directions of two atoms of a stereo bond chain
pub fn triple_prod_char(
    at: &[InpAtom],
    at_1: usize,
    i_next_at_1: usize,
    z_dir1: &[i8; 3],
    at_2: usize,
    i_next_at_2: usize,
    z_dir2: &[i8; 3],
) -> i32 {
    let a1 = &at[at_1];
    let a2 = &at[a1.neighbor[i_next_at_1] as usize];
    let mut pnt0 = [a2.x - a1.x, a2.y - a1.y, a2.z - a1.z];
    let b2 = &at[at_2];
    let b1 = &at[b2.neighbor[i_next_at_2] as usize];
    let mut pnt1 = [b2.x - b1.x, b2.y - b1.y, b2.z - b1.z];
    pnt1 = add3(&pnt0, &pnt1);
    let mut pnt2 = [0.0; 3];
    for i in 0..3 {
        pnt0[i] = f64::from(z_dir1[i]);
        pnt2[i] = f64::from(z_dir2[i]);
    }
    let mut pnts = [pnt0, pnt1, pnt2];
    for i in 0..3 {
        let mut len = len3(&pnts[i]);
        if len < MIN_BOND_LEN {
            if i == 1 && (at[at_1].b_used_0d_parity != 0 || at[at_2].b_used_0d_parity != 0) {
                pnts[i] = [0.0, 1.0, 0.0];
                len = 1.0;
            } else {
                return 0;
            }
        }
        pnts[i] = mult3(&pnts[i], 1.0 / len);
    }
    let len = 100.0 * triple_prod(&pnts[0], &pnts[1], &pnts[2], None);
    if len >= 0.0 {
        (len + 0.5).floor() as i32
    } else {
        -((0.5 - len).floor() as i32)
    }
}

// ---------------------------------------------------------------------------
// stereo center / stereo bond candidate checks
// ---------------------------------------------------------------------------

pub const PES_BIT_POINT_EDGE_STEREO: i32 = 1;
pub const PES_BIT_PHOSPHINE_STEREO: i32 = 2;
pub const PES_BIT_ARSINE_STEREO: i32 = 4;
pub const PES_BIT_FIX_SP3_BUG: i32 = 8;
pub const MAX_NUM_STEREO_BOND_NEIGH: i32 = 3;
pub const MIN_NUM_STEREO_BOND_NEIGH: i32 = 2;
pub const CMODE_NO_ALT_SBONDS: u64 = 0x000100;

use super::ichiring::is_atom_in_3memb_ring;
use super::ichitaut::{BOND_ALT12NS, BOND_ALTERN, BOND_MARK_ALL, BOND_TAUTOM};
use super::incomdef::{BOND_TYPE_SINGLE, RADICAL_SINGLET};
use super::util::get_endpoint_valence;

const BOND_SINGLE: i32 = BOND_TYPE_SINGLE;
const BOND_DOUBLE: i32 = BOND_TYPE_DOUBLE;

fn radical_ok(radical: i8) -> bool {
    radical == 0 || i32::from(radical) == RADICAL_SINGLET
}

fn elem_table_check(table: &[(&str, i8)], elname: &str, charge: i8, radical: i8) -> bool {
    table
        .iter()
        .find(|(e, c)| *e == elname && *c == charge)
        .is_some_and(|_| radical_ok(radical))
}

/// `bCanAtomHaveAStereoBond`
pub fn b_can_atom_have_a_stereo_bond(elname: &str, charge: i8, radical: i8) -> bool {
    const T: [(&str, i8); 5] = [("C", 0), ("Si", 0), ("Ge", 0), ("N", 0), ("N", 1)];
    elem_table_check(&T, elname, charge, radical)
}

/// `bCanAtomBeMiddleAllene`
pub fn b_can_atom_be_middle_allene(elname: &str, charge: i8, radical: i8) -> bool {
    const T: [(&str, i8); 3] = [("C", 0), ("Si", 0), ("Ge", 0)];
    elem_table_check(&T, elname, charge, radical)
}

/// `bCanAtomBeTerminalAllene`
pub fn b_can_atom_be_terminal_allene(elname: &str, charge: i8, radical: i8) -> bool {
    b_can_atom_be_middle_allene(elname, charge, radical)
}

/// `bInpAtomHasRequirdNeigh`
pub fn b_inp_atom_has_requird_neigh(
    at: &[InpAtom],
    cur_at: usize,
    requird_neigh_type: i32,
    num_dble_bonds: i32,
    b_stereo_at_zz: i32,
) -> bool {
    let a = &at[cur_at];
    if a.endpoint != 0 {
        return false;
    }
    if requird_neigh_type & 1 != 0 && a.num_h != 0 {
        return false;
    }
    let val = a.valence as usize;
    if requird_neigh_type & 2 != 0 {
        for i in 0..val {
            let ni = &at[a.neighbor[i] as usize];
            if ni.valence != 1 || get_endpoint_valence(ni.el_number) == 0 {
                continue;
            }
            for j in i + 1..val {
                let nj = &at[a.neighbor[j] as usize];
                if nj.valence != 1
                    || ni.el_number != nj.el_number
                    || get_endpoint_valence(nj.el_number) == 0
                {
                    continue;
                }
                if ni.num_h != nj.num_h && ni.elname_str() != "C" {
                    return false;
                }
            }
        }
    }
    if b_stereo_at_zz == 0 {
        for i in 0..val {
            let n = at[a.neighbor[i] as usize].elname_str();
            if n == "Zz" || n == "Zy" {
                return false;
            }
        }
    }
    let (mut num_1s, mut num_mult, mut num_other) = (0, 0, 0);
    for i in 0..val {
        match i32::from(a.bond_type[i]) & !BOND_MARK_ALL {
            BOND_SINGLE => num_1s += 1,
            BOND_DOUBLE | BOND_ALTERN | BOND_TAUTOM | BOND_ALT12NS => num_mult += 1,
            _ => num_other += 1,
        }
    }
    if num_other != 0 {
        return false;
    }
    if (num_dble_bonds != 0 && num_dble_bonds > num_mult)
        || (num_dble_bonds == 0 && i32::from(a.valence) != num_1s)
    {
        return false;
    }
    true
}

/// `bCanInpAtomBeAStereoCenter`: the number of bonds and implicit H (3 or 4)
/// of a possible stereocenter, otherwise 0
pub fn b_can_inp_atom_be_a_stereo_center(
    at: &[InpAtom],
    cur_at: usize,
    b_pointed_edge_stereo: i32,
    b_stereo_at_zz: i32,
) -> i32 {
    // (element, charge, bonds+H, chemical valence+H, has 3-membered ring, required neighbors)
    const T: [(&str, i8, i32, i32, bool, i32); 21] = [
        ("C", 0, 4, 4, false, 0),
        ("Si", 0, 4, 4, false, 0),
        ("Ge", 0, 4, 4, false, 0),
        ("Sn", 0, 4, 4, false, 0),
        ("As", 1, 4, 4, false, 3),
        ("B", -1, 4, 4, false, 0),
        ("S", 0, 3, 4, false, 3),
        ("S", 0, 4, 6, false, 3),
        ("S", 1, 3, 3, false, 3),
        ("S", 1, 4, 5, false, 3),
        ("Se", 0, 3, 4, false, 3),
        ("Se", 0, 4, 6, false, 3),
        ("Se", 1, 3, 3, false, 3),
        ("Se", 1, 4, 5, false, 3),
        ("N", 0, 4, 5, false, 3),
        ("N", 1, 4, 4, false, 3),
        ("N", 0, 3, 3, true, 1),
        ("P", 1, 4, 4, false, 3),
        ("P", 0, 4, 5, false, 3),
        ("P", 0, 3, 3, false, 2),
        ("As", 0, 3, 3, false, 2),
    ];
    const PHOSPHINE_STEREO: usize = 19;
    const ARSINE_STEREO: usize = 20;
    let a = &at[cur_at];
    let mut ret = 0;
    let mut found = T.len();
    for (i, &(el, charge, nb, cv, ring3, req)) in T.iter().enumerate() {
        // (ALLOW_NO_CHARGE_ON_STEREO_CENTERS is disabled upstream)
        if a.elname_str() == el
            && a.charge == charge
            && (a.radical == 0 || a.radical == 1)
            && i32::from(a.valence) + i32::from(a.num_h) == nb
            && i32::from(a.chem_bonds_valence) + i32::from(a.num_h) == cv
            && (!ring3 || is_atom_in_3memb_ring(at, cur_at))
            && b_inp_atom_has_requird_neigh(at, cur_at, req, cv - nb, b_stereo_at_zz)
        {
            ret = nb;
            found = i;
            break;
        }
    }
    if found == PHOSPHINE_STEREO && b_pointed_edge_stereo & PES_BIT_PHOSPHINE_STEREO == 0 {
        ret = 0;
    }
    if found == ARSINE_STEREO && b_pointed_edge_stereo & PES_BIT_ARSINE_STEREO == 0 {
        ret = 0;
    }
    ret
}

/// `bIsSuitableHeteroInpAtom`: 1 for `=N-`/`=NH`, 0 for `=O`, -1 otherwise
pub fn b_is_suitable_hetero_inp_atom(a: &InpAtom) -> i32 {
    if a.charge == 0 && radical_ok(a.radical) {
        let val = get_endpoint_valence(a.el_number);
        if val > 0 {
            let num_h = i32::from(a.num_h);
            if val == i32::from(a.chem_bonds_valence) + num_h {
                match val {
                    2 if num_h == 0 && a.valence == 1 => return 0,
                    3 if (a.valence == 1 && num_h == 1) || (a.valence == 2 && num_h == 0) => {
                        return 1;
                    }
                    _ => {}
                }
            }
        }
    }
    -1
}

/// `bIsOxide`; like the C code it clears the mark bits of the bond types of `cur_at`
pub fn b_is_oxide(at: &mut [InpAtom], cur_at: usize) -> bool {
    for i in 0..at[cur_at].valence as usize {
        at[cur_at].bond_type[i] &= !(BOND_MARK_ALL as u8);
        let bond_type = i32::from(at[cur_at].bond_type[i]);
        let an = &at[at[cur_at].neighbor[i] as usize];
        if bond_type == BOND_DOUBLE {
            if an.valence == 1
                && an.charge == 0
                && an.num_h == 0
                && an.radical == 0
                && get_endpoint_valence(an.el_number) == 2
            {
                return true;
            }
        } else if (bond_type == BOND_TAUTOM || bond_type == BOND_ALT12NS)
            && an.valence == 1
            && get_endpoint_valence(an.el_number) == 2
        {
            return true;
        }
    }
    false
}

/// `get_allowed_stereo_bond_type`
pub fn get_allowed_stereo_bond_type(bond_type: i32) -> i32 {
    let bt = bond_type & !BOND_MARK_ALL;
    match bt {
        BOND_ALTERN | BOND_DOUBLE | BOND_SINGLE | BOND_TAUTOM => bt,
        _ => 0,
    }
}

/// `half_stereo_bond_action`
pub fn half_stereo_bond_action(
    n_parity: i32,
    b_unknown: bool,
    b_isotopic: bool,
    v_ab_parity_unknown: i32,
) -> i32 {
    const AB_NEGATIVE: i32 = 0x10;
    const AB_UNKNOWN: i32 = 0x20;
    let mut p = n_parity;
    if p == AB_PARITY_NONE {
        return AB_PARITY_NONE;
    }
    if p == v_ab_parity_unknown {
        p = AB_PARITY_ODD | AB_UNKNOWN;
    }
    if p == -v_ab_parity_unknown {
        p = AB_PARITY_ODD | AB_UNKNOWN | AB_NEGATIVE;
    }
    if p < 0 {
        p = (if p == -AB_PARITY_EVEN {
            AB_PARITY_ODD
        } else {
            -p
        }) | AB_NEGATIVE;
    } else if p == AB_PARITY_EVEN {
        p = AB_PARITY_ODD;
    }
    if b_unknown {
        p |= AB_UNKNOWN;
    }
    const ODD: i32 = AB_PARITY_ODD;
    const UNDF: i32 = AB_PARITY_UNDF;
    const IISO: i32 = AB_PARITY_IISO;
    if b_isotopic {
        match p {
            ODD | 0x11 => AB_PARITY_CALC,
            0x21 | 0x24 | 0x31 | 0x34 => v_ab_parity_unknown,
            IISO | 0x25 => AB_PARITY_NONE,
            UNDF | 0x14 => AB_PARITY_UNDF,
            _ => -1,
        }
    } else {
        match p {
            ODD => AB_PARITY_CALC,
            0x21 | 0x24 => v_ab_parity_unknown,
            UNDF => AB_PARITY_UNDF,
            0x31 | 0x11 | IISO | 0x25 | 0x14 | 0x34 => AB_PARITY_NONE,
            _ => -1,
        }
    }
}

/// `can_be_a_stereo_bond_with_isotopic_H`: the number of possible stereo
/// bonds of `cur_at` (to atoms with smaller numbers)
pub fn can_be_a_stereo_bond_with_isotopic_h(at: &mut [InpAtom], cur_at: usize, n_mode: u64) -> i32 {
    fn type_n(a: &InpAtom) -> i32 {
        if a.num_h == 0 && a.charge == 0 && a.radical == 0 && get_endpoint_valence(a.el_number) == 3
        {
            if a.valence == 2 && a.chem_bonds_valence == 3 {
                return 1;
            } else if a.valence == 3 && a.chem_bonds_valence == 5 {
                return 2;
            }
        }
        0
    }
    let nbh = |a: &InpAtom| i32::from(a.valence) + i32::from(a.num_h);
    if nbh(&at[cur_at]) > MAX_NUM_STEREO_BOND_NEIGH || nbh(&at[cur_at]) < MIN_NUM_STEREO_BOND_NEIGH
    {
        return 0;
    }
    if !b_can_atom_have_a_stereo_bond(
        at[cur_at].elname_str(),
        at[cur_at].charge,
        at[cur_at].radical,
    ) {
        return 0;
    }
    let (mut num_2s, mut num_alt, mut num_wrong_bonds_1) = (0, 0, 0);
    let mut num_2s_hetero = [0i32; 2];
    let type_n_cur = type_n(&at[cur_at]);
    let mut num_stereo_bonds = 0;
    for i in 0..at[cur_at].valence as usize {
        let mut b_found = false;
        let next_at = at[cur_at].neighbor[i] as usize;
        let bond_type = get_allowed_stereo_bond_type(i32::from(at[cur_at].bond_type[i]));
        if bond_type == BOND_ALTERN {
            num_alt += 1;
            if cur_at > next_at && n_mode & CMODE_NO_ALT_SBONDS == 0 {
                b_found = true;
            }
        } else if bond_type == BOND_DOUBLE {
            num_2s += 1;
            let n2sh = b_is_suitable_hetero_inp_atom(&at[next_at]);
            if n2sh >= 0 {
                num_2s_hetero[n2sh as usize] += 1;
            }
            if cur_at > next_at {
                b_found = true;
            }
        } else if bond_type != BOND_SINGLE && bond_type != BOND_TAUTOM {
            num_wrong_bonds_1 += 1;
            if num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 != 0 && at[cur_at].valence <= 2) {
                return 0;
            }
            continue;
        }
        if b_found {
            if nbh(&at[next_at]) > MAX_NUM_STEREO_BOND_NEIGH
                || nbh(&at[next_at]) < MIN_NUM_STEREO_BOND_NEIGH
            {
                continue;
            }
            if !b_can_atom_have_a_stereo_bond(
                at[next_at].elname_str(),
                at[next_at].charge,
                at[next_at].radical,
            ) {
                continue;
            }
            let (mut num_2s_next, mut num_alt_next) = (0, 0);
            let mut num_2s_hetero_next = [0i32; 2];
            let type_n_next = type_n(&at[next_at]);
            let nv = at[next_at].valence as usize;
            let mut j = 0;
            while j < nv {
                let bond_type = get_allowed_stereo_bond_type(i32::from(at[next_at].bond_type[j]));
                if bond_type == BOND_ALTERN {
                    num_alt_next += 1;
                } else if bond_type == BOND_DOUBLE {
                    num_2s_next += 1;
                    let n2sh = b_is_suitable_hetero_inp_atom(&at[at[next_at].neighbor[j] as usize]);
                    if n2sh >= 0 {
                        num_2s_hetero_next[n2sh as usize] += 1;
                    }
                } else if bond_type != BOND_SINGLE && bond_type != BOND_TAUTOM {
                    // (sic) the C code tests the counter of the first atom here
                    if num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 != 0 && at[cur_at].valence <= 2)
                    {
                        break;
                    }
                }
                j += 1;
            }
            if 3 == (type_n_cur | type_n_next)
                && ((type_n_cur == 2 && !b_is_oxide(at, cur_at))
                    || (type_n_next == 2 && !b_is_oxide(at, next_at)))
            {
                b_found = false;
            } else if j < nv
                || i32::from(num_alt_next > 0) + i32::from(num_2s_next > 0) != 1
                || num_2s_next > 2
            {
                b_found = false;
            } else if num_2s_next == 2 && at[next_at].valence != 2 {
                let n = &at[next_at];
                if 1 == (num_2s_hetero_next[0] | num_2s_hetero_next[1])
                    && nbh(n) == 3
                    && i32::from(n.chem_bonds_valence) + i32::from(n.num_h) == 5
                    && get_endpoint_valence(n.el_number) == 3
                    && (type_n_cur == 0 || b_is_oxide(at, next_at))
                {
                } else {
                    b_found = false;
                }
            }
        }
        if b_found {
            num_stereo_bonds += 1;
        }
    }
    if i32::from(num_alt > 0) + i32::from(num_2s > 0) != 1 || num_stereo_bonds == 0 {
        return 0;
    }
    if num_2s > 1 {
        let a = &at[cur_at];
        if num_2s == 2
            && 1 == (num_2s_hetero[0] | num_2s_hetero[1])
            && nbh(a) == 3
            && i32::from(a.chem_bonds_valence) + i32::from(a.num_h) == 5
            && get_endpoint_valence(a.el_number) == 3
        {
        } else {
            return 0;
        }
    }
    num_stereo_bonds
}

// ---------------------------------------------------------------------------
// half stereo bond parity
// ---------------------------------------------------------------------------

pub const MIN_ANGLE_DBOND: f64 = 0.087156;
pub const MAX_SINE: f64 = std::f64::consts::FRAC_1_SQRT_2;
pub const AB_PARITY_0D: i32 = 8;

/// `bAtomHasValence3`
pub fn b_atom_has_valence3(elname: &str, charge: i8, radical: i8) -> bool {
    elname == "N" && charge == 0 && radical_ok(radical)
}

/// `GetHalfStereobond0DParity`
pub fn get_half_stereobond_0d_parity(
    at: &mut [InpAtom],
    cur_at: usize,
    sb_neigh_orig_at_numb: &[u16],
    n_num_explict_attachments: usize,
    mut bond_parity: i32,
    n_flag: i8,
) -> i32 {
    let mut last_parity = 0i32;
    let mut m = 0;
    while m < MAX_NUM_STEREO_BONDS && at[cur_at].sb_parity[m] != 0 {
        let a = &at[cur_at];
        let mut cur_parity = 0i32;
        let mut icur2nxt: i32 = -1;
        let mut icur2neigh: i32 = -1;
        let ord = i32::from(a.sb_ord[m]);
        let valid = ord >= 0 && ord < i32::from(a.valence) && {
            let nxt_at = a.neighbor[ord as usize] as usize;
            at[nxt_at].valence as usize <= MAX_NUM_STEREO_BONDS && at[nxt_at].orig_at_number != 0
        };
        if !valid {
            m += 1;
            continue;
        }
        let nxt_at = a.neighbor[ord as usize] as usize;
        let next_orig = at[nxt_at].orig_at_number;
        for (i, &orig) in sb_neigh_orig_at_numb
            .iter()
            .enumerate()
            .take(n_num_explict_attachments)
        {
            if a.sn_orig_at_num[m] == orig {
                icur2neigh = i as i32;
            } else if next_orig == orig {
                icur2nxt = i as i32;
            }
        }
        if icur2neigh >= 0 && icur2nxt >= 0 {
            let sbp = i32::from(a.sb_parity[m]);
            if atom_parity_well_def(sbp) {
                let cur_order_parity =
                    (icur2nxt + icur2neigh + i32::from(icur2nxt > icur2neigh) - 1) % 2;
                cur_parity = 2 - (cur_order_parity + sbp) % 2;
            } else {
                cur_parity = sbp;
            }
        }
        if last_parity == 0 {
            last_parity = cur_parity;
        } else if last_parity != cur_parity && cur_parity != 0 {
            if atom_parity_well_def(last_parity) {
                if atom_parity_well_def(cur_parity) {
                    last_parity = 0;
                    break;
                }
            } else if atom_parity_well_def(cur_parity) {
                last_parity = cur_parity;
            } else {
                last_parity = cur_parity.min(last_parity);
            }
        }
        m += 1;
    }
    if last_parity != 0 {
        bond_parity = last_parity;
        at[cur_at].b_used_0d_parity |= n_flag;
    }
    bond_parity
}

/// `half_stereo_bond_parity`: the parity of one end of a stereo bond.
/// The removed explicit terminal H atoms are `at[removed_h_start..removed_h_start + num_removed_h]`
/// (the C code passes a pointer into the same array).
#[allow(clippy::too_many_arguments)]
pub fn half_stereo_bond_parity(
    at: &mut [InpAtom],
    cur_at: usize,
    removed_h_start: usize,
    num_removed_h: usize,
    mut z_dir: Option<&mut [i8; 3]>,
    b_pointed_edge_stereo: i32,
    v_ab_parity_unknown: i32,
) -> i32 {
    use super::ichierr::CT_ISO_H_ERR;
    use super::incomdef::NUM_H_ISOTOPES;
    use super::inpdef::{AT_FLAG_ISO_H_POINT, FLAG_SB_0D};
    let one_pi = std::f64::consts::PI;
    let two_pi = 2.0 * one_pi;
    let mut at_coord = [[0.0f64; 3]; 3];
    let mut tmp = [0.0f64; 3];
    let mut pnt = [[0.0f64; 3]; 3];
    let mut bond_parity = AB_PARITY_UNDF;
    let mut num_ex_h = 0usize;
    let mut num_either_single = 0;
    let mut num_iso_h = [0i32; NUM_H_ISOTOPES + 1];
    let mut index_h = [0usize; 5];
    let b_ignore_isotopic_h = at[cur_at].c_flags & AT_FLAG_ISO_H_POINT != 0;
    let mut sb_neigh_orig_at_numb = [0u16; 3];

    if let Some(z) = z_dir.as_deref_mut() {
        if z[0] == 0 && z[1] == 0 && z[2] == 0 {
            z[2] = 100;
        }
    }
    let num_h = i32::from(at[cur_at].num_h);
    if num_h > NUM_H_ISOTOPES as i32 {
        return 0;
    }
    let valence = i32::from(at[cur_at].valence);
    if valence + num_h > MAX_NUM_STEREO_BOND_NEIGH || valence + num_h < MIN_NUM_STEREO_BOND_NEIGH {
        return 0;
    }
    if !b_can_atom_have_a_stereo_bond(
        at[cur_at].elname_str(),
        at[cur_at].charge,
        at[cur_at].radical,
    ) {
        return 0;
    }
    let mut num_n_h = 0;
    if !b_ignore_isotopic_h {
        num_n_h = num_h;
        for j in 0..NUM_H_ISOTOPES {
            let k = i32::from(at[cur_at].num_iso_h[j]);
            if k > 1 {
                return AB_PARITY_IISO;
            }
            num_n_h -= k;
        }
    }
    if num_n_h > 1 {
        return AB_PARITY_IISO;
    }
    if num_n_h < 0 {
        return CT_ISO_H_ERR;
    }

    // the "goto exit_function" paths are the `'body` breaks
    'body: {
        if num_removed_h > 0 {
            for j in 0..num_removed_h {
                let h = &at[removed_h_start + j];
                if h.neighbor[0] as usize == cur_at {
                    let k = if b_ignore_isotopic_h {
                        0
                    } else {
                        i32::from(h.iso_atw_diff)
                    };
                    if (0..=NUM_H_ISOTOPES as i32).contains(&k) {
                        num_iso_h[k as usize] += 1;
                        if num_iso_h[k as usize] > 1 {
                            return CT_ISO_H_ERR;
                        }
                        index_h[num_ex_h] = j;
                        num_ex_h += 1;
                    } else {
                        return CT_ISO_H_ERR;
                    }
                }
            }
            let num_i_h = num_h - num_ex_h as i32;
            if num_i_h > 1 {
                bond_parity = -AB_PARITY_UNDF;
                break 'body;
            }
        }
        if valence + (num_ex_h as i32) < MIN_NUM_STEREO_BOND_NEIGH {
            return if num_h == 1 {
                AB_PARITY_UNDF
            } else {
                -AB_PARITY_UNDF
            };
        }
        let b_valence3 = b_atom_has_valence3(
            at[cur_at].elname_str(),
            at[cur_at].charge,
            at[cur_at].radical,
        );
        if valence + num_ex_h as i32 == MIN_NUM_STEREO_BOND_NEIGH
            && !b_valence3
            && !((num_n_h == 1 && num_iso_h[0] == 0) || (num_h == 1 && num_ex_h == 0))
        {
            break 'body;
        }

        // coordinates of the neighbors
        let mut n_attach = 0usize;
        for k in 0..2 {
            let count = if k == 0 { num_ex_h } else { valence as usize };
            for j in 0..count {
                let (next, z0, n_type);
                if k == 0 {
                    let nx = removed_h_start + index_h[j];
                    next = nx;
                    at_coord[n_attach][0] = at[nx].x - at[cur_at].x;
                    at_coord[n_attach][1] = at[nx].y - at[cur_at].y;
                    sb_neigh_orig_at_numb[n_attach] = at[nx].orig_at_number;
                    let (zz, t) = get_z_coord(
                        at,
                        nx,
                        0,
                        -(b_pointed_edge_stereo & PES_BIT_POINT_EDGE_STEREO),
                    );
                    z0 = -zz;
                    n_type = t;
                } else {
                    next = at[cur_at].neighbor[j] as usize;
                    at_coord[n_attach][0] = at[next].x - at[cur_at].x;
                    at_coord[n_attach][1] = at[next].y - at[cur_at].y;
                    sb_neigh_orig_at_numb[n_attach] = at[next].orig_at_number;
                    let (zz, t) = get_z_coord(
                        at,
                        cur_at,
                        j,
                        b_pointed_edge_stereo & PES_BIT_POINT_EDGE_STEREO,
                    );
                    z0 = zz;
                    n_type = t;
                }
                let _ = next;
                let mut z = z0;
                let mut n_type = n_type;
                match n_type {
                    ZTYPE_EITHER => num_either_single += 1,
                    ZTYPE_UP | ZTYPE_DOWN => {
                        if k == 0 {
                            n_type = -n_type;
                        }
                        z = len2(&[at_coord[n_attach][0], at_coord[n_attach][1]]);
                        if n_type == ZTYPE_DOWN {
                            z = -z;
                        }
                    }
                    _ => {}
                }
                at_coord[n_attach][2] = z;
                n_attach += 1;
            }
        }
        if num_either_single != 0 {
            bond_parity = v_ab_parity_unknown;
            break 'body;
        }
        if n_attach == 2 {
            for j in 0..3 {
                at_coord[2][j] = -(at_coord[0][j] + at_coord[1][j]);
            }
            // (nSbNeighOrigAtNumb[2] = 0: implicit H or lone pair)
            sb_neigh_orig_at_numb[2] = 0;
        }
        for j in 0..3 {
            tmp[j] = len3(&at_coord[j]);
        }
        let mut min_tmp = tmp[0].min(tmp[1].min(tmp[2]));
        let mut max_tmp = tmp[0].max(tmp[1].max(tmp[2]));
        if min_tmp < MIN_BOND_LEN || min_tmp < MIN_SINE * max_tmp {
            if at[cur_at].sb_parity[0] != 0 {
                bond_parity = get_half_stereobond_0d_parity(
                    at,
                    cur_at,
                    &sb_neigh_orig_at_numb,
                    n_attach,
                    bond_parity,
                    FLAG_SB_0D,
                );
            }
            break 'body;
        }
        for j in 0..3 {
            at_coord[j] = mult3(&at_coord[j], 1.0 / tmp[j]);
        }
        for j in 0..3 {
            pnt[j] = diff3(&at_coord[j], &at_coord[(j + 1) % 3]);
            tmp[j] = len3(&pnt[j]);
            if tmp[j] < MIN_SINE {
                break 'body;
            }
            pnt[j] = mult3(&pnt[j], 1.0 / tmp[j]);
        }
        let mut tmp1 = len3(&cross_prod3(&pnt[0], &pnt[1]));
        let mut k = 0usize;
        for j in 1..3 {
            let tmp2 = len3(&cross_prod3(&pnt[j], &pnt[(j + 1) % 3]));
            if tmp2 > tmp1 {
                tmp1 = tmp2;
                k = j;
            }
        }
        let p0 = k;
        let p1 = (k + 1) % 3;
        let p2 = (k + 2) % 3;
        pnt[p2] = cross_prod3(&pnt[p0], &pnt[p1]);
        tmp[p2] = len3(&pnt[p2]);
        if tmp[p2] < MIN_SINE * tmp[p0] * tmp[p1] {
            break 'body;
        }
        pnt[p2] = mult3(
            &pnt[p2],
            (if pnt[p2][2] > 0.0 { 1.0 } else { -1.0 }) / tmp[p2],
        );
        min_tmp = dot_prod3(&at_coord[0], &pnt[p2]);
        pnt[p0] = mult3(&pnt[p2], min_tmp);
        pnt[p0] = diff3(&at_coord[0], &pnt[p0]);
        tmp[p0] = len3(&pnt[p0]);
        pnt[p0] = mult3(&pnt[p0], 1.0 / tmp[p0]);
        pnt[p1] = cross_prod3(&pnt[p2], &pnt[p0]);
        for j in 0..3 {
            let temp = at_coord[j];
            for k in 0..3 {
                at_coord[j][k] = dot_prod3(&temp, &pnt[(k + p0) % 3]);
            }
            tmp[j] = (at_coord[j][0] * at_coord[j][0] + at_coord[j][1] * at_coord[j][1]).sqrt();
            at_coord[j] = mult3(&at_coord[j], 1.0 / tmp[j]);
        }
        let s = (at_coord[1][0] * at_coord[2][1] - at_coord[1][1] * at_coord[2][0]).abs();
        let c = at_coord[1][0] * at_coord[2][0] + at_coord[1][1] * at_coord[2][1];
        if s < MIN_SINE && c > 0.5 {
            break 'body;
        }
        let c = at_coord[0][0];
        let s = at_coord[0][1];
        for j in 0..MAX_NUM_STEREO_BOND_NEIGH as usize {
            let t1 = c * at_coord[j][0] + s * at_coord[j][1];
            let t2 = -s * at_coord[j][0] + c * at_coord[j][1];
            at_coord[j][0] = t1;
            at_coord[j][1] = t2;
        }
        tmp1 = at_coord[1][1].atan2(at_coord[1][0]);
        let mut tmp2 = at_coord[2][1].atan2(at_coord[2][0]);
        if tmp1 < 0.0 {
            tmp1 += two_pi;
        }
        if tmp2 < 0.0 {
            tmp2 += two_pi;
        }
        bond_parity = 2 - i32::from(tmp1 < tmp2);
        if let Some(z) = z_dir {
            for j in 0..3 {
                let v = pnt[p2][j];
                z[j] = (if v >= 0.0 {
                    (0.5 + 100.0 * v).floor()
                } else {
                    -(0.5 - 100.0 * v).floor()
                }) as i8;
            }
        }
        if n_attach > 2 {
            min_tmp = tmp1.min(tmp2);
            max_tmp = tmp1.max(tmp2);
            if min_tmp > one_pi - MIN_SINE
                || max_tmp < one_pi + MIN_SINE
                || max_tmp - min_tmp > one_pi - MIN_SINE
            {
                at[cur_at].b_ambiguous_stereo |= AMBIGUOUS_STEREO as i8;
            } else if at_coord[0][2].abs() > MAX_SINE {
                at[cur_at].b_ambiguous_stereo |= AMBIGUOUS_STEREO as i8;
            }
        } else if n_attach == 2 {
            min_tmp = (tmp1 - one_pi).abs();
            if min_tmp < MIN_SINE {
                bond_parity = AB_PARITY_UNDF;
            } else if min_tmp < MIN_ANGLE_DBOND {
                at[cur_at].b_ambiguous_stereo |= AMBIGUOUS_STEREO as i8;
            }
        }
    }

    if num_h > 1 && bond_parity > 0 && bond_parity & AB_PARITY_0D == 0 {
        bond_parity = -bond_parity;
    }
    bond_parity
}

// ---------------------------------------------------------------------------
// stereo bonds
// ---------------------------------------------------------------------------

pub const MAX_CUMULENE_LEN: i32 = 2;
pub const MULT_STEREOBOND: i32 = 0x08;
pub const MASK_CUMULENE_LEN: i32 = 0x38;
pub const MIN_DOT_PROD: i32 = 50;
pub const MIN_SB_RING_SIZE: u16 = 8;

/// The stereo part of `sp_ATOM` (extr_ct.h) used by the stereo detection
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SpAtom {
    pub neighbor: [u16; super::inpdef::MAXVAL],
    pub init_rank: u16,
    pub orig_at_number: u16,
    pub valence: i8,
    pub b_has_stereo_or_equ_to_stereo: i8,
    pub b_has_stereo_or_equ_to_stereo2: i8,
    pub endpoint: u16,
    pub num_h: i8,
    pub num_iso_h: [i8; 3],
    pub n_ring_system: u16,
    pub stereo_bond_neighbor: [u16; MAX_NUM_STEREO_BONDS],
    pub stereo_bond_neighbor2: [u16; MAX_NUM_STEREO_BONDS],
    pub stereo_bond_ord: [i8; MAX_NUM_STEREO_BONDS],
    pub stereo_bond_ord2: [i8; MAX_NUM_STEREO_BONDS],
    pub stereo_bond_z_prod: [i8; MAX_NUM_STEREO_BONDS],
    pub stereo_bond_z_prod2: [i8; MAX_NUM_STEREO_BONDS],
    pub stereo_bond_parity: [i8; MAX_NUM_STEREO_BONDS],
    pub stereo_bond_parity2: [i8; MAX_NUM_STEREO_BONDS],
    pub parity: i8,
    pub parity2: i8,
    pub stereo_atom_parity: i8,
    pub stereo_atom_parity2: i8,
    pub final_parity: i8,
    pub final_parity2: i8,
    pub b_ambiguous_stereo: i8,
    pub z_dir: [i8; 3],
}

/// `save_a_stereo_bond`: stores a stereo bond in both end atoms; 0 if there is no room
#[allow(clippy::too_many_arguments)]
pub fn save_a_stereo_bond(
    z_prod: i32,
    result_action: i32,
    at1: usize,
    ord1: i32,
    nb1: &mut [u16; MAX_NUM_STEREO_BONDS],
    sord1: &mut [i8; MAX_NUM_STEREO_BONDS],
    zp1: &mut [i8; MAX_NUM_STEREO_BONDS],
    par1: &mut [i8; MAX_NUM_STEREO_BONDS],
    at2: usize,
    ord2: i32,
    nb2: &mut [u16; MAX_NUM_STEREO_BONDS],
    sord2: &mut [i8; MAX_NUM_STEREO_BONDS],
    zp2: &mut [i8; MAX_NUM_STEREO_BONDS],
    par2: &mut [i8; MAX_NUM_STEREO_BONDS],
) -> bool {
    let i1 = nb1
        .iter()
        .position(|&x| x == 0)
        .unwrap_or(MAX_NUM_STEREO_BONDS);
    let i2 = nb2
        .iter()
        .position(|&x| x == 0)
        .unwrap_or(MAX_NUM_STEREO_BONDS);
    if i1 == MAX_NUM_STEREO_BONDS || i2 == MAX_NUM_STEREO_BONDS {
        return false;
    }
    par1[i1] = result_action as i8;
    par2[i2] = result_action as i8;
    nb1[i1] = (at2 + 1) as u16;
    sord1[i1] = ord1 as i8;
    nb2[i2] = (at1 + 1) as u16;
    sord2[i2] = ord2 as i8;
    zp1[i1] = z_prod as i8;
    zp2[i2] = z_prod as i8;
    true
}

/// `FixSb0DParities`: reconciles the 0D parities of both ends of a stereo bond
#[allow(clippy::too_many_arguments)]
pub fn fix_sb_0d_parities(
    at: &[InpAtom],
    chain_length: i32,
    at_1: usize,
    i_next_at_1: i32,
    z_dir1: &mut [i8; 3],
    at_2: usize,
    i_next_at_2: i32,
    z_dir2: &mut [i8; 3],
    pparity1: &mut i32,
    pparity2: &mut i32,
) -> i32 {
    use super::inpdef::FLAG_SB_0D;
    let mut parity1 = AB_PARITY_NONE;
    let mut parity2 = AB_PARITY_NONE;
    let (mut j1, mut j2) = (-1i32, -1i32);
    let parity_sign = if *pparity1 < 0 || *pparity2 < 0 {
        -1
    } else {
        1
    };
    let mut abs_parity1 = pparity1.abs();
    let mut abs_parity2 = pparity2.abs();
    let wd = atom_parity_well_def;
    let mut k = 0;
    while k < MAX_NUM_STEREO_BONDS && at[at_1].sb_parity[k] != 0 {
        if i32::from(at[at_1].sb_ord[k]) == i_next_at_1 {
            parity1 = i32::from(at[at_1].sb_parity[k]);
            j1 = k as i32;
        }
        k += 1;
    }
    k = 0;
    while k < MAX_NUM_STEREO_BONDS && at[at_2].sb_parity[k] != 0 {
        if i32::from(at[at_2].sb_ord[k]) == i_next_at_2 {
            parity2 = i32::from(at[at_2].sb_parity[k]);
            j2 = k as i32;
        }
        k += 1;
    }
    let pick = |p: i32, ap: i32| {
        if wd(p) {
            ap
        } else if wd(ap) {
            p
        } else {
            ap.min(p)
        }
    };
    match i32::from(j1 >= 0) + 2 * i32::from(j2 >= 0) {
        0 => {
            *pparity1 = parity_sign * AB_PARITY_UNDF;
            *pparity2 = parity_sign * AB_PARITY_UNDF;
            return 0;
        }
        1 | 2 => {
            *pparity1 = AB_PARITY_NONE;
            *pparity2 = AB_PARITY_NONE;
            return -1;
        }
        _ => match i32::from(!(wd(abs_parity1) && wd(parity1)))
            + 2 * i32::from(!(wd(abs_parity2) && wd(parity2)))
        {
            0 => {}
            1 => {
                *pparity1 = parity_sign * pick(parity1, abs_parity1);
                *pparity2 = parity_sign * abs_parity2;
                return -1;
            }
            2 => {
                *pparity1 = parity_sign * abs_parity1;
                *pparity2 = parity_sign * pick(parity2, abs_parity2);
                return -1;
            }
            _ => {
                abs_parity1 = pick(parity1, abs_parity1);
                abs_parity2 = pick(parity2, abs_parity2);
                let v = parity_sign * abs_parity1.min(abs_parity2);
                *pparity1 = v;
                *pparity2 = v;
                return -1;
            }
        },
    }
    *pparity1 = parity_sign * abs_parity1;
    *pparity2 = parity_sign * abs_parity2;

    if chain_length % 2 != 0 {
        let wrong1 = at[at_1].b_used_0d_parity & FLAG_SB_0D != 0;
        let wrong2 = at[at_2].b_used_0d_parity & FLAG_SB_0D != 0;
        'set_default: {
            if wrong1 && wrong2 {
                break 'set_default;
            }
            if wrong1 || wrong2 {
                let mut r12 = [
                    at[at_2].x - at[at_1].x,
                    at[at_2].y - at[at_1].y,
                    at[at_2].z - at[at_1].z,
                ];
                let abs_r12 = len3(&r12);
                if abs_r12 < MIN_BOND_LEN {
                    break 'set_default;
                }
                let zi1;
                if wrong1 {
                    zi1 = [
                        f64::from(z_dir2[0]),
                        f64::from(z_dir2[1]),
                        f64::from(z_dir2[2]),
                    ];
                    r12 = mult3(&r12, 1.0 / abs_r12);
                } else {
                    zi1 = [
                        f64::from(z_dir1[0]),
                        f64::from(z_dir1[1]),
                        f64::from(z_dir1[2]),
                    ];
                    r12 = mult3(&r12, -1.0 / abs_r12);
                }
                let mut zi2 = cross_prod3(&r12, &zi1);
                let abs_zi2 = len3(&zi2);
                zi2 = mult3(&zi2, 100.0 / abs_zi2);
                let mut z_dir = [0i8; 3];
                for j in 0..3 {
                    z_dir[j] = (if zi2[j] >= 0.0 {
                        (0.5 + zi2[j]).floor()
                    } else {
                        -(0.5 - zi2[j]).floor()
                    }) as i8;
                }
                if wrong1 {
                    *z_dir1 = z_dir;
                } else {
                    *z_dir2 = z_dir;
                }
            }
            return 0;
        }
        *z_dir1 = [100, 0, 0];
        *z_dir2 = [0, 0, 100];
    }
    0
}

/// `set_stereo_bonds_parity`: finds the stereo bonds (and cumulenes) that
/// start at `at_1` and stores their parities in `out_at`. Returns the number
/// of stored isotopic stereo bonds or a negative error code.
#[allow(clippy::too_many_arguments)]
pub fn set_stereo_bonds_parity(
    out_at: &mut [SpAtom],
    at: &mut [InpAtom],
    at_1: usize,
    removed_h_start: usize,
    num_removed_h: usize,
    n_mode: u64,
    q: &mut super::ichiring::Queue,
    n_atom_level: &mut [u16],
    c_source: &mut [i8],
    min_sb_ring_size: u16,
    b_pointed_edge_stereo: i32,
    v_ab_parity_unknown: i32,
) -> i32 {
    use super::ichierr::{CT_CALC_STEREO_ERR, CT_ERR_MAX, CT_ERR_MIN, CT_STEREOBOND_ERROR};
    use super::inpdef::FLAG_SB_0D;
    let nbh = |a: &InpAtom| i32::from(a.valence) + i32::from(a.num_h);
    let type_n_of = |a: &InpAtom| -> i32 {
        if a.num_h == 0 && a.charge == 0 && a.radical == 0 && get_endpoint_valence(a.el_number) == 3
        {
            if a.valence == 2 && a.chem_bonds_valence == 3 {
                return 1;
            } else if a.valence == 3 && a.chem_bonds_valence == 5 {
                return 2;
            }
        }
        0
    };
    if nbh(&at[at_1]) > MAX_NUM_STEREO_BOND_NEIGH || nbh(&at[at_1]) < MIN_NUM_STEREO_BOND_NEIGH {
        return 0;
    }
    if !b_can_atom_have_a_stereo_bond(at[at_1].elname_str(), at[at_1].charge, at[at_1].radical) {
        return 0;
    }
    if at[at_1].c_point != 0 {
        return 0;
    }
    if at[at_1].num_h == 0
        && at[at_1].valence == 2
        && get_allowed_stereo_bond_type(i32::from(at[at_1].bond_type[0])) == BOND_DOUBLE
        && get_allowed_stereo_bond_type(i32::from(at[at_1].bond_type[1])) == BOND_DOUBLE
    {
        return 0;
    }

    let (mut num_2s_1, mut num_alt_1) = (0, 0);
    let mut chain_length = 0;
    let mut num_chains = 0;
    let mut num_wrong_bonds_1 = 0;
    let mut num_2s_hetero = [0i32; 2];
    let type_n = type_n_of(&at[at_1]);
    let mut all_pos_1 = [0i32; MAX_NUM_STEREO_BONDS];
    let mut all_pos_2 = [0i32; MAX_NUM_STEREO_BONDS];
    let mut all_at_2 = [0usize; MAX_NUM_STEREO_BONDS];
    let mut all_unkn = [0i8; MAX_NUM_STEREO_BONDS];
    let mut num_stereo_bonds = 0usize;
    let mut i_next_at_2: i32 = -1;

    for i_next_at_1 in 0..at[at_1].valence as usize {
        let mut n_unknown = i32::from(at[at_1].bond_stereo[i_next_at_1] == STEREO_DBLE_EITHER);
        let mut bond_type =
            get_allowed_stereo_bond_type(i32::from(at[at_1].bond_type[i_next_at_1]));
        let mut at_2: i32 = -1;
        let mut next_at_2 = at_1;
        if bond_type == BOND_ALTERN || bond_type == BOND_DOUBLE {
            at_2 = i32::from(at[at_1].neighbor[i_next_at_1]);
            next_at_2 = at_1;
        }
        match bond_type {
            BOND_ALTERN => {
                num_alt_1 += 1;
                let a2 = at_2 as usize;
                if at[at_1].n_ring_system != at[a2].n_ring_system {
                    continue;
                }
                if n_mode & CMODE_NO_ALT_SBONDS != 0
                    || !b_can_atom_have_a_stereo_bond(
                        at[a2].elname_str(),
                        at[a2].charge,
                        at[a2].radical,
                    )
                {
                    continue;
                }
            }
            BOND_DOUBLE => {
                num_2s_1 += 1;
                let mut cur_chain_length = 0;
                if b_can_atom_be_terminal_allene(
                    at[at_1].elname_str(),
                    at[at_1].charge,
                    at[at_1].radical,
                ) {
                    let mut b_allene;
                    loop {
                        let a2 = at_2 as usize;
                        b_allene = at[a2].num_h == 0
                            && at[a2].valence == 2
                            && get_allowed_stereo_bond_type(i32::from(at[a2].bond_type[0]))
                                == BOND_DOUBLE
                            && get_allowed_stereo_bond_type(i32::from(at[a2].bond_type[1]))
                                == BOND_DOUBLE;
                        if !(b_allene
                            && b_can_atom_be_middle_allene(
                                at[a2].elname_str(),
                                at[a2].charge,
                                at[a2].radical,
                            ))
                        {
                            break;
                        }
                        let k = usize::from(at[a2].neighbor[0] as usize == next_at_2);
                        next_at_2 = a2;
                        n_unknown += i32::from(at[a2].bond_stereo[k] == STEREO_DBLE_EITHER);
                        at_2 = i32::from(at[a2].neighbor[k]);
                        cur_chain_length += 1;
                    }
                    if cur_chain_length != 0 {
                        num_chains += 1;
                        let a2 = at_2 as usize;
                        if b_allene
                            || !b_can_atom_be_terminal_allene(
                                at[a2].elname_str(),
                                at[a2].charge,
                                at[a2].radical,
                            )
                        {
                            continue;
                        }
                        chain_length = cur_chain_length;
                    }
                }
                let a2 = at_2 as usize;
                if cur_chain_length == 0 {
                    let n2sh = b_is_suitable_hetero_inp_atom(&at[a2]);
                    if n2sh >= 0 {
                        num_2s_hetero[n2sh as usize] += 1;
                    }
                    if !b_can_atom_have_a_stereo_bond(
                        at[a2].elname_str(),
                        at[a2].charge,
                        at[a2].radical,
                    ) {
                        continue;
                    }
                }
            }
            BOND_SINGLE | BOND_TAUTOM => continue,
            _ => {
                num_wrong_bonds_1 += 1;
                continue;
            }
        }
        bond_type = 0;
        let _ = bond_type;

        let mut b_found = at_2 >= 0 && at_1 as i32 > at_2;
        if b_found {
            let a2 = at_2 as usize;
            if nbh(&at[a2]) > MAX_NUM_STEREO_BOND_NEIGH || nbh(&at[a2]) < MIN_NUM_STEREO_BOND_NEIGH
            {
                continue;
            }
            let (mut num_2s_2, mut num_alt_2) = (0, 0);
            let mut num_2s_hetero_next = [0i32; 2];
            let type_n_next = type_n_of(&at[a2]);
            i_next_at_2 = -1;
            let mut num_wrong_bonds_2 = 0;
            for j in 0..at[a2].valence as usize {
                let bt = get_allowed_stereo_bond_type(i32::from(at[a2].bond_type[j]));
                if bt == 0 {
                    num_wrong_bonds_2 += 1;
                    continue;
                }
                if bt == BOND_DOUBLE {
                    num_2s_2 += 1;
                    let n2sh = b_is_suitable_hetero_inp_atom(&at[at[a2].neighbor[j] as usize]);
                    if n2sh >= 0 {
                        num_2s_hetero_next[n2sh as usize] += 1;
                    }
                } else {
                    num_alt_2 += i32::from(bt == BOND_ALTERN);
                }
                if at[a2].neighbor[j] as usize == next_at_2 {
                    i_next_at_2 = j as i32;
                }
            }
            if num_wrong_bonds_2 > 1
                || (num_wrong_bonds_2 != 0 && at[a2].valence <= 2)
                || i32::from(num_alt_2 > 0) + i32::from(num_2s_2 > 0) != 1
                || i_next_at_2 < 0
            {
                b_found = false;
            } else if at[a2].c_point != 0 {
                b_found = false;
            } else if num_2s_2 > 2 {
                b_found = false;
            } else if 3 == (type_n | type_n_next)
                && ((type_n == 2 && !b_is_oxide(at, at_1))
                    || (type_n_next == 2 && !b_is_oxide(at, a2)))
            {
                b_found = false;
            } else if num_2s_2 == 2 {
                let a = &at[a2];
                if chain_length == 0
                    && 1 == (num_2s_hetero_next[0] | num_2s_hetero_next[1])
                    && nbh(a) == 3
                    && i32::from(a.chem_bonds_valence) + i32::from(a.num_h) == 5
                    && get_endpoint_valence(a.el_number) == 3
                    && (type_n == 0 || b_is_oxide(at, a2))
                {
                } else {
                    b_found = false;
                }
            }
            if chain_length != 0 && num_alt_2 != 0 {
                return 0;
            }
        }
        if b_found {
            all_pos_1[num_stereo_bonds] = i_next_at_1 as i32;
            all_pos_2[num_stereo_bonds] = i_next_at_2;
            all_at_2[num_stereo_bonds] = at_2 as usize;
            all_unkn[num_stereo_bonds] = n_unknown as i8;
            num_stereo_bonds += 1;
        }
    }
    if num_chains > 1 {
        return 0;
    }
    if num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 != 0 && at[at_1].valence <= 2) {
        return 0;
    }
    if chain_length != 0
        && (num_stereo_bonds != 1 || num_alt_1 != 0 || chain_length > MAX_CUMULENE_LEN)
    {
        return 0;
    }
    if i32::from(num_alt_1 > 0) + i32::from(num_2s_1 > 0) != 1 || num_stereo_bonds == 0 {
        return 0;
    }
    if num_2s_1 > 1 {
        let a = &at[at_1];
        if num_2s_1 == 2
            && type_n == 2
            && 1 == (num_2s_hetero[0] | num_2s_hetero[1])
            && nbh(a) == 3
            && i32::from(a.chem_bonds_valence) + i32::from(a.num_h) == 5
            && get_endpoint_valence(a.el_number) == 3
        {
        } else {
            return 0;
        }
    }

    // ================== calculate the parities ==================
    let mut num_stored_isotopic_stereo_bonds = 0;
    let mut stop = false;
    let returned_error = |v: i32| (CT_ERR_MIN..=CT_ERR_MAX).contains(&v);
    for k in 0..num_stereo_bonds {
        let chain_len_bits = chain_length * MULT_STEREOBOND;
        let at_2 = all_at_2[k];
        let i_next_at_1 = all_pos_1[k];
        if at[at_1].n_ring_system == at[at_2].n_ring_system {
            let j = super::ichiring::is_bond_in_nmax_memb_ring(
                at,
                at_1,
                i_next_at_1 as usize,
                q,
                n_atom_level,
                c_source,
                min_sb_ring_size,
            );
            if j > 0 {
                continue;
            } else if j < 0 {
                return CT_STEREOBOND_ERROR;
            }
        }
        let i_next_at_2 = all_pos_2[k];
        let n_unknown = i32::from(all_unkn[k]);
        let mut z_dir1 = [0i8; 3];
        let mut z_dir2 = [0i8; 3];
        let mut cur_parity = half_stereo_bond_parity(
            at,
            at_1,
            removed_h_start,
            num_removed_h,
            Some(&mut z_dir1),
            b_pointed_edge_stereo,
            v_ab_parity_unknown,
        );
        let mut next_parity = half_stereo_bond_parity(
            at,
            at_2,
            removed_h_start,
            num_removed_h,
            Some(&mut z_dir2),
            b_pointed_edge_stereo,
            v_ab_parity_unknown,
        );
        if returned_error(cur_parity) || returned_error(next_parity) {
            return CT_CALC_STEREO_ERR;
        }
        if (at[at_1].b_used_0d_parity & FLAG_SB_0D) != 0
            || (at[at_2].b_used_0d_parity & FLAG_SB_0D) != 0
        {
            fix_sb_0d_parities(
                at,
                chain_length,
                at_1,
                i_next_at_1,
                &mut z_dir1,
                at_2,
                i_next_at_2,
                &mut z_dir2,
                &mut cur_parity,
                &mut next_parity,
            );
        }
        if cur_parity == AB_PARITY_NONE || cur_parity.abs() == AB_PARITY_IISO {
            continue;
        }
        if next_parity == AB_PARITY_NONE || next_parity.abs() == AB_PARITY_IISO {
            continue;
        }
        let mut cur_action =
            half_stereo_bond_action(cur_parity, n_unknown != 0, false, v_ab_parity_unknown);
        let mut next_action =
            half_stereo_bond_action(next_parity, n_unknown != 0, false, v_ab_parity_unknown);
        let mut result_action = cur_action.min(next_action);
        if result_action == -1 {
            stop = true;
        }
        let abs_cur_parity = cur_parity.abs();
        let abs_next_parity = next_parity.abs();
        let cur_parity_defined = atom_parity_well_def(abs_cur_parity);
        let next_parity_defined = atom_parity_well_def(abs_next_parity);
        let dot_prod_z;
        if cur_parity_defined && next_parity_defined {
            dot_prod_z = if chain_len_bits != 0
                && ((chain_len_bits & MASK_CUMULENE_LEN) / MULT_STEREOBOND) % 2 != 0
            {
                triple_prod_char(
                    at,
                    at_1,
                    i_next_at_1 as usize,
                    &z_dir1,
                    at_2,
                    i_next_at_2 as usize,
                    &z_dir2,
                )
            } else {
                dot_prodchar3(&z_dir1, &z_dir2)
            };
            if dot_prod_z.abs() < MIN_DOT_PROD {
                result_action = result_action.min(AB_PARITY_UNDF);
            }
        } else {
            dot_prod_z = 0;
        }
        let (lo, hi) = if at_1 < at_2 {
            (at_1, at_2)
        } else {
            (at_2, at_1)
        };
        let (left, right) = out_at.split_at_mut(hi);
        let (o_lo, o_hi) = (&mut left[lo], &mut right[0]);
        // (o1, o2) are the records of at_1 and at_2
        let swap = at_1 > at_2;
        macro_rules! recs {
            () => {
                if swap {
                    (&mut *o_hi, &mut *o_lo)
                } else {
                    (&mut *o_lo, &mut *o_hi)
                }
            };
        }
        if result_action != AB_PARITY_NONE
            && result_action != -1
            && cur_parity > 0
            && next_parity > 0
        {
            let (o1, o2) = recs!();
            if save_a_stereo_bond(
                dot_prod_z,
                result_action | chain_len_bits,
                at_1,
                i_next_at_1,
                &mut o1.stereo_bond_neighbor,
                &mut o1.stereo_bond_ord,
                &mut o1.stereo_bond_z_prod,
                &mut o1.stereo_bond_parity,
                at_2,
                i_next_at_2,
                &mut o2.stereo_bond_neighbor,
                &mut o2.stereo_bond_ord,
                &mut o2.stereo_bond_z_prod,
                &mut o2.stereo_bond_parity,
            ) {
                if o1.parity == 0
                    || (cur_parity_defined && !atom_parity_well_def(i32::from(o1.parity).abs()))
                {
                    o1.parity = cur_parity as i8;
                    o1.z_dir = z_dir1;
                }
                if o2.parity == 0
                    || (next_parity_defined && !atom_parity_well_def(i32::from(o2.parity).abs()))
                {
                    o2.parity = next_parity as i8;
                    o2.z_dir = z_dir2;
                }
                o1.b_ambiguous_stereo |= at[at_1].b_ambiguous_stereo;
                o2.b_ambiguous_stereo |= at[at_2].b_ambiguous_stereo;
            }
        }
        cur_action = half_stereo_bond_action(cur_parity, n_unknown != 0, true, v_ab_parity_unknown);
        next_action =
            half_stereo_bond_action(next_parity, n_unknown != 0, true, v_ab_parity_unknown);
        result_action = cur_action.min(next_action);
        cur_parity = abs_cur_parity;
        next_parity = abs_next_parity;
        if result_action != AB_PARITY_NONE && result_action != -1 {
            if cur_parity > 0 && next_parity > 0 {
                let (o1, o2) = recs!();
                if save_a_stereo_bond(
                    dot_prod_z,
                    result_action | chain_len_bits,
                    at_1,
                    i_next_at_1,
                    &mut o1.stereo_bond_neighbor2,
                    &mut o1.stereo_bond_ord2,
                    &mut o1.stereo_bond_z_prod2,
                    &mut o1.stereo_bond_parity2,
                    at_2,
                    i_next_at_2,
                    &mut o2.stereo_bond_neighbor2,
                    &mut o2.stereo_bond_ord2,
                    &mut o2.stereo_bond_z_prod2,
                    &mut o2.stereo_bond_parity2,
                ) {
                    if o1.parity2 == 0
                        || (cur_parity_defined
                            && !atom_parity_well_def(i32::from(o1.parity2).abs()))
                    {
                        o1.parity2 = cur_parity as i8;
                        if o1.parity == 0 {
                            o1.z_dir = z_dir1;
                        }
                    }
                    if o2.parity2 == 0
                        || (next_parity_defined
                            && !atom_parity_well_def(i32::from(o2.parity2).abs()))
                    {
                        o2.parity2 = next_parity as i8;
                        if o2.parity == 0 {
                            o2.z_dir = z_dir2;
                        }
                    }
                    o1.b_ambiguous_stereo |= at[at_1].b_ambiguous_stereo;
                    o2.b_ambiguous_stereo |= at[at_2].b_ambiguous_stereo;
                    num_stored_isotopic_stereo_bonds += 1;
                }
            }
        } else if result_action == -1 {
            stop = true;
        }
    }
    if stop {
        return CT_CALC_STEREO_ERR;
    }
    num_stored_isotopic_stereo_bonds
}

// ---------------------------------------------------------------------------
// stereo atoms (stereocenters)
// ---------------------------------------------------------------------------

pub const MIN_ANGLE_RELAXED: f64 = 0.001;
pub const MIN_SINE_RELAXED: f64 = 0.001;
pub const ZERO_FLOAT: f64 = 1.0e-12;
pub const AMBIGUOUS_STEREO_ERROR: i32 = 32;
pub const FLAG_SC_0D: i8 = 1;
pub const REQ_MODE_MIN_SB_RING_MASK: u64 = 0x0F0000;
pub const REQ_MODE_MIN_SB_RING_SHFT: u32 = 16;
const ADD_EXPLICIT_HYDROGEN_NEIGH: i32 = 1;
const ADD_EXPLICIT_LONE_PAIR_NEIGH: i32 = 2;

/// `can_be_a_stereo_atom_with_isotopic_H`
pub fn can_be_a_stereo_atom_with_isotopic_h(
    at: &[InpAtom],
    cur_at: usize,
    b_pointed_edge_stereo: i32,
    b_stereo_at_zz: i32,
) -> i32 {
    let n = b_can_inp_atom_be_a_stereo_center(at, cur_at, b_pointed_edge_stereo, b_stereo_at_zz);
    i32::from(
        n != 0
            && i32::from(at[cur_at].valence) + i32::from(at[cur_at].num_h) == n
            && at[cur_at].num_h as usize <= super::incomdef::NUM_H_ISOTOPES,
    )
}

/// `GetStereocenter0DParity`; `sb_neigh_orig_at_numb[..j1]` gets sorted
pub fn get_stereocenter_0d_parity(
    at: &mut [InpAtom],
    cur_at: usize,
    j1: usize,
    sb_neigh_orig_at_numb: &mut [u16],
    n_flag: i8,
) -> i32 {
    let mut parity = AB_PARITY_NONE;
    let n4 = MAX_NUM_STEREO_ATOM_NEIGH;
    if at[cur_at].p_parity != 0 && (j1 == n4 - 1 || j1 == n4) {
        let mut inp = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
        for i in 0..n4 {
            inp[i] = at[cur_at].p_orig_at_num[i];
            if inp[i] == at[cur_at].orig_at_number {
                inp[i] = 0;
            }
        }
        let cmp = |a: &u16, b: &u16| i32::from(*a) - i32::from(*b);
        let num_trans_inp = insertions_sort(&mut inp[..], cmp);
        let num_trans_neigh = insertions_sort(&mut sb_neigh_orig_at_numb[..j1], cmp);
        if inp[n4 - j1..] == sb_neigh_orig_at_numb[..j1] {
            let pp = i32::from(at[cur_at].p_parity);
            parity = if atom_parity_well_def(pp) {
                2 - (num_trans_inp + num_trans_neigh + pp) % 2
            } else {
                pp
            };
            at[cur_at].b_used_0d_parity |= n_flag;
        }
    }
    parity
}

/// `set_stereo_atom_parity`: the parity of a stereocenter in `out_at[cur_at]`;
/// returns `out_at[cur_at].parity2`
#[allow(clippy::too_many_arguments)]
pub fn set_stereo_atom_parity(
    out_at: &mut [SpAtom],
    at: &mut [InpAtom],
    cur_at: usize,
    removed_h_start: usize,
    num_removed_h: usize,
    b_pointed_edge_stereo: i32,
    v_ab_parity_unknown: i32,
    loose_tsa_check: bool,
    b_stereo_at_zz: i32,
) -> i32 {
    use super::incomdef::NUM_H_ISOTOPES;
    use super::inpdef::AT_FLAG_ISO_H_POINT;
    let n4 = MAX_NUM_STEREO_ATOM_NEIGH;
    let mut num_explicit_iso_h = [0i32; NUM_H_ISOTOPES + 1];
    let mut index_h = [0usize; MAX_NUM_STEREO_ATOM_NEIGH];
    let mut num_explicit_h = 0usize;
    let mut at_coord = [[0.0f64; 3]; MAX_NUM_STEREO_ATOM_NEIGH];
    let mut bond_len_xy = [0.0f64; 4];
    let (mut rmax, mut rmin) = (0.0f64, 0.0f64);
    let mut sum_xyz = [0.0f64; 3];
    let mut at_coord_center = [0.0f64; 3];
    let mut parity = AB_PARITY_NONE;
    let mut b_ambiguous = 0i32;
    let mut b_add_explicit_neighbor = 0;
    let mut b2d = 0;
    let mut n2d_tetrahedral_ambiguity = 0;
    let b_ignore_isotopic_h = at[cur_at].c_flags & AT_FLAG_ISO_H_POINT != 0;
    let mut sb_neigh_orig_at_numb = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
    let mut v_min_angle = MIN_ANGLE;
    let mut v_min_sine = MIN_SINE;
    if loose_tsa_check && at[cur_at].n_num_at_in_ring_system >= 3 {
        v_min_angle = MIN_ANGLE_RELAXED;
        v_min_sine = MIN_SINE_RELAXED;
    }
    out_at[cur_at].parity = 0;
    out_at[cur_at].parity2 = 0;
    out_at[cur_at].stereo_atom_parity = 0;
    out_at[cur_at].stereo_atom_parity2 = 0;
    let mut tot_num_iso_h = 0;
    let mut j1 = 0usize;

    'body: {
        let n_must_have_num_neigh =
            b_can_inp_atom_be_a_stereo_center(at, cur_at, b_pointed_edge_stereo, b_stereo_at_zz)
                as usize;
        if n_must_have_num_neigh == 0 || at[cur_at].num_h as usize > NUM_H_ISOTOPES {
            break 'body;
        }
        for j in 0..NUM_H_ISOTOPES {
            if at[cur_at].num_iso_h[j] > 1 {
                break 'body;
            }
            tot_num_iso_h += i32::from(at[cur_at].num_iso_h[j]);
        }
        if b_ignore_isotopic_h {
            tot_num_iso_h = 0;
        }
        if i32::from(at[cur_at].num_h) - tot_num_iso_h > 1 {
            break 'body;
        }
        if num_removed_h > 0 {
            for j in 0..num_removed_h {
                let h = &at[removed_h_start + j];
                if h.neighbor[0] as usize == cur_at {
                    let mut k = i32::from(h.iso_atw_diff);
                    if k < 0 || k > NUM_H_ISOTOPES as i32 || b_ignore_isotopic_h {
                        k = 0;
                    }
                    num_explicit_iso_h[k as usize] += 1;
                    // (C: index_H[4] could overflow with more than 4 removed H on one atom)
                    if num_explicit_h < index_h.len() {
                        index_h[num_explicit_h] = j;
                    }
                    num_explicit_h += 1;
                }
            }
        }
        let mut num_z = 0;
        // coordinates of the neighbors
        'coords: for k in 0..2 {
            let count = if k == 0 {
                num_explicit_h
            } else {
                at[cur_at].valence as usize
            };
            for j in 0..count {
                let (z, mut n_type, nbr_idx);
                if k == 0 {
                    let nx = removed_h_start + index_h[j];
                    nbr_idx = nx;
                    let (zz, t) = get_z_coord(
                        at,
                        nx,
                        0,
                        -(b_pointed_edge_stereo & PES_BIT_POINT_EDGE_STEREO),
                    );
                    z = -zz;
                    n_type = t;
                } else {
                    nbr_idx = at[cur_at].neighbor[j] as usize;
                    let (zz, t) = get_z_coord(
                        at,
                        cur_at,
                        j,
                        b_pointed_edge_stereo & PES_BIT_POINT_EDGE_STEREO,
                    );
                    z = zz;
                    n_type = t;
                }
                match n_type {
                    ZTYPE_EITHER => {
                        parity = v_ab_parity_unknown;
                        break 'body;
                    }
                    ZTYPE_UP | ZTYPE_DOWN => {
                        if k == 0 {
                            n_type = -n_type;
                        }
                        b2d += 1;
                        num_z += 1;
                    }
                    ZTYPE_3D => num_z += 1,
                    _ => {}
                }
                if j1 < n4 {
                    sb_neigh_orig_at_numb[j1] = at[nbr_idx].orig_at_number;
                    at_coord[j1][0] = at[nbr_idx].x - at[cur_at].x;
                    at_coord[j1][1] = at[nbr_idx].y - at[cur_at].y;
                    bond_len_xy[j1] = len2(&[at_coord[j1][0], at_coord[j1][1]]);
                    at_coord[j1][2] = if n_type == ZTYPE_3D {
                        z
                    } else if n_type == ZTYPE_UP {
                        bond_len_xy[j1]
                    } else if n_type == ZTYPE_DOWN {
                        -bond_len_xy[j1]
                    } else {
                        0.0
                    };
                    j1 += 1;
                } else {
                    continue 'coords;
                }
            }
        }
        let b2d_flag = b2d == num_z && num_z != 0;
        b2d = i32::from(b2d_flag);
        let valence = at[cur_at].valence as usize;
        if n4 != valence + num_explicit_h && n4 - 1 != valence + num_explicit_h {
            if n_must_have_num_neigh == valence + at[cur_at].num_h as usize && at[cur_at].num_h > 1
            {
                parity = if parity == v_ab_parity_unknown {
                    -v_ab_parity_unknown
                } else {
                    -AB_PARITY_UNDF
                };
            } else {
                parity = AB_PARITY_NONE;
            }
            break 'body;
        }
        for j in 0..j1 {
            let z = len3(&at_coord[j]);
            if z < MIN_BOND_LEN {
                parity = get_stereocenter_0d_parity(
                    at,
                    cur_at,
                    j1,
                    &mut sb_neigh_orig_at_numb,
                    FLAG_SC_0D,
                );
                if parity == AB_PARITY_NONE {
                    parity = AB_PARITY_UNDF;
                }
                break 'body;
            }
            at_coord[j] = mult3(&at_coord[j], 1.0 / z);
            rmax = if j != 0 { rmax.max(z) } else { z };
            rmin = if j != 0 { rmin.min(z) } else { z };
        }
        if rmin / rmax < MIN_SINE {
            parity =
                get_stereocenter_0d_parity(at, cur_at, j1, &mut sb_neigh_orig_at_numb, FLAG_SC_0D);
            if parity == AB_PARITY_NONE {
                parity = AB_PARITY_UNDF;
            }
            break 'body;
        }
        for j in 0..j1 {
            sum_xyz = add3(&sum_xyz, &at_coord[j]);
        }
        let j = j1;
        if j1 == n4 - 1 {
            if n_must_have_num_neigh == n4 - 1 {
                b_add_explicit_neighbor = ADD_EXPLICIT_LONE_PAIR_NEIGH;
            } else if n_must_have_num_neigh == n4
                && 1 == i32::from(at[cur_at].num_h) - num_explicit_h as i32
                && 1 == i32::from(at[cur_at].num_h) - tot_num_iso_h
            {
                b_add_explicit_neighbor = ADD_EXPLICIT_HYDROGEN_NEIGH;
            }
        }
        if b_add_explicit_neighbor != 0 {
            if parity == v_ab_parity_unknown {
                break 'body;
            }
            if num_z == 0 || are_3_vect_in_one_plane(&at_coord, v_min_sine) {
                if num_z > 0 {
                    b_ambiguous |= AMBIGUOUS_STEREO;
                }
                // (APPLY_IMPLICIT_H_DOWN_RULE is 0, ALWAYS_SET_STEREO_PARITY is 0)
                parity = get_stereocenter_0d_parity(
                    at,
                    cur_at,
                    j1,
                    &mut sb_neigh_orig_at_numb,
                    FLAG_SC_0D,
                );
                if parity == AB_PARITY_NONE {
                    parity = AB_PARITY_UNDF;
                }
                break 'body;
            } else {
                at_coord[j] = change_sign3(&sum_xyz);
                let z = len3(&at_coord[j]);
                // (FIX_STEREO_SCALING_BUG is 0)
                if b_pointed_edge_stereo & PES_BIT_FIX_SP3_BUG != 0 {
                    if z > 1.0 {
                        rmax *= z;
                    } else {
                        rmin *= z;
                    }
                } else {
                    rmax = rmax.max(z);
                    rmin = rmin.min(z);
                }
                if z < MIN_BOND_LEN || rmin / rmax < MIN_SINE {
                    parity = get_stereocenter_0d_parity(
                        at,
                        cur_at,
                        j1,
                        &mut sb_neigh_orig_at_numb,
                        FLAG_SC_0D,
                    );
                    if parity == AB_PARITY_NONE {
                        parity = AB_PARITY_UNDF;
                    }
                    break 'body;
                }
            }
        } else if j1 != n4 {
            if parity == v_ab_parity_unknown {
                parity = -AB_PARITY_UNDF;
            }
            break 'body;
        } else if num_z == 0 || are_4at_in_one_plane(&at_coord, v_min_sine) {
            if num_z > 0 {
                b_ambiguous |= AMBIGUOUS_STEREO;
            }
            if parity != v_ab_parity_unknown {
                parity = get_stereocenter_0d_parity(
                    at,
                    cur_at,
                    j1,
                    &mut sb_neigh_orig_at_numb,
                    FLAG_SC_0D,
                );
                if parity == AB_PARITY_NONE {
                    parity = AB_PARITY_UNDF;
                } else if atom_parity_well_def(parity) {
                    b_ambiguous &= !AMBIGUOUS_STEREO;
                }
            }
            break 'body;
        }
        // 4 neighbors; the 2D tetrahedral ambiguity
        if b2d != 0 {
            n2d_tetrahedral_ambiguity = get_2d_tetrahedral_ambiguity(
                &at_coord,
                b_add_explicit_neighbor != 0,
                b_pointed_edge_stereo & PES_BIT_FIX_SP3_BUG != 0,
                v_min_angle,
            );
            if n2d_tetrahedral_ambiguity > 0 {
                if T2D_WARN & n2d_tetrahedral_ambiguity != 0 {
                    b_ambiguous |= AMBIGUOUS_STEREO;
                }
                if T2D_UNDF & n2d_tetrahedral_ambiguity != 0 {
                    if parity != v_ab_parity_unknown {
                        parity = AB_PARITY_UNDF;
                    }
                    break 'body;
                }
            } else if n2d_tetrahedral_ambiguity < 0 {
                b_ambiguous |= AMBIGUOUS_STEREO_ERROR;
                parity = AB_PARITY_UNDF;
                break 'body;
            }
        }
        // origin to the neighbor #0
        let c0 = at_coord[0];
        for j in 1..n4 {
            at_coord[j] = diff3(&at_coord[j], &c0);
        }
        at_coord_center = diff3(&at_coord_center, &c0);
        let mut min_sine = 0.0;
        let mut amb = b_ambiguous;
        let triple_product = triple_prod_and_min_abs_sine2(
            &at_coord[1..],
            &at_coord_center,
            b_add_explicit_neighbor != 0,
            Some(&mut min_sine),
            Some(&mut amb),
            v_min_sine,
        );
        b_ambiguous = amb;
        if triple_product.abs() > ZERO_FLOAT
            && (min_sine > v_min_sine
                || (min_sine.abs() > ZERO_FLOAT && (n2d_tetrahedral_ambiguity & T2D_OKAY) != 0))
        {
            parity = if triple_product > 0.0 {
                AB_PARITY_EVEN
            } else {
                AB_PARITY_ODD
            };
            if j1 == n4 - 1 {
                parity = 3 - parity;
            }
        } else {
            if num_z > 0 {
                b_ambiguous |= AMBIGUOUS_STEREO;
            }
            parity = AB_PARITY_UNDF;
        }
    }

    // exit_function
    if parity != 0 {
        out_at[cur_at].b_ambiguous_stereo |= b_ambiguous as i8;
        at[cur_at].b_ambiguous_stereo |= b_ambiguous as i8;
    }
    if at[cur_at].num_h > 1 || parity <= 0 {
    } else {
        out_at[cur_at].parity = parity as i8;
    }
    if parity == -AB_PARITY_UNDF || parity == -v_ab_parity_unknown {
        parity = -parity;
    }
    if parity < 0 {
        parity = AB_PARITY_NONE;
    }
    out_at[cur_at].parity2 = parity as i8;
    let p = i32::from(out_at[cur_at].parity) & 0x07;
    out_at[cur_at].stereo_atom_parity = (if atom_parity_well_def(p) {
        AB_PARITY_CALC
    } else {
        p
    }) as i8;
    let p = i32::from(out_at[cur_at].parity2) & 0x07;
    out_at[cur_at].stereo_atom_parity2 = (if atom_parity_well_def(p) {
        AB_PARITY_CALC
    } else {
        p
    }) as i8;
    let _ = num_explicit_iso_h;
    i32::from(out_at[cur_at].parity2)
}

/// `set_stereo_parity`: the stereo descriptors of all atoms in `at_output`.
/// `at[num_at..num_at + num_removed_h]` are the removed explicit H atoms.
/// Returns the number of well-defined stereocenters or a negative error code.
#[allow(clippy::too_many_arguments)]
pub fn set_stereo_parity(
    at: &mut [InpAtom],
    at_output: &mut [SpAtom],
    num_at: usize,
    num_removed_h: usize,
    n_max_num_stereo_atoms: Option<&mut i32>,
    n_max_num_stereo_bonds: Option<&mut i32>,
    n_mode: u64,
    b_pointed_edge_stereo: i32,
    v_ab_parity_unknown: i32,
    b_loose_tsa_check: bool,
    b_stereo_at_zz: i32,
) -> i32 {
    use super::ichierr::{CT_ERR_MAX, CT_ERR_MIN};
    use super::ichiring::Queue;
    let mut num_3d_stereo_atoms = 0;
    let mut num_stereo_bonds = 0;
    let (mut max_stereo_atoms, mut max_stereo_bonds) = (0, 0);
    for s in at_output.iter_mut().take(num_at) {
        s.parity = 0;
        s.parity2 = 0;
        s.stereo_bond_neighbor = Default::default();
        s.stereo_bond_neighbor2 = Default::default();
        s.stereo_bond_ord = Default::default();
        s.stereo_bond_ord2 = Default::default();
        s.stereo_bond_z_prod = Default::default();
        s.stereo_bond_z_prod2 = Default::default();
        s.stereo_bond_parity = Default::default();
        s.stereo_bond_parity2 = Default::default();
    }
    let (mut want_atoms, mut want_bonds) = (n_max_num_stereo_atoms, n_max_num_stereo_bonds);
    if want_atoms.is_some() || want_bonds.is_some() {
        for i in 0..num_at {
            let num =
                can_be_a_stereo_atom_with_isotopic_h(at, i, b_pointed_edge_stereo, b_stereo_at_zz);
            if num != 0 {
                max_stereo_atoms += num;
            } else {
                let num = can_be_a_stereo_bond_with_isotopic_h(at, i, n_mode);
                if num != 0 {
                    max_stereo_bonds += num;
                }
            }
        }
        if let Some(a) = want_atoms.as_deref_mut() {
            *a = max_stereo_atoms;
        }
        if let Some(b) = want_bonds.as_deref_mut() {
            *b = max_stereo_bonds;
        }
    }
    let mut min_sb_ring_size =
        ((n_mode & REQ_MODE_MIN_SB_RING_MASK) >> REQ_MODE_MIN_SB_RING_SHFT) as u16;
    let (mut q, mut n_atom_level, mut c_source);
    if min_sb_ring_size >= 3 {
        q = Queue::create(num_at as i32 + 1).unwrap();
        n_atom_level = vec![0u16; num_at];
        c_source = vec![0i8; num_at];
    } else {
        min_sb_ring_size = 2;
        // not used: the ring search returns at once for sizes below 3
        q = Queue::create(1).unwrap();
        n_atom_level = Vec::new();
        c_source = Vec::new();
    }
    for i in 0..num_at {
        let mut is_stereo = set_stereo_atom_parity(
            at_output,
            at,
            i,
            num_at,
            num_removed_h,
            b_pointed_edge_stereo,
            v_ab_parity_unknown,
            b_loose_tsa_check,
            b_stereo_at_zz,
        );
        if is_stereo != 0 {
            num_3d_stereo_atoms += i32::from(atom_parity_well_def(is_stereo));
        } else {
            is_stereo = set_stereo_bonds_parity(
                at_output,
                at,
                i,
                num_at,
                num_removed_h,
                n_mode,
                &mut q,
                &mut n_atom_level,
                &mut c_source,
                min_sb_ring_size,
                b_pointed_edge_stereo,
                v_ab_parity_unknown,
            );
            if (CT_ERR_MIN..=CT_ERR_MAX).contains(&is_stereo) {
                num_3d_stereo_atoms = is_stereo;
                break;
            }
            num_stereo_bonds += i32::from(is_stereo != 0);
        }
    }
    if let Some(a) = want_atoms {
        if max_stereo_atoms < num_3d_stereo_atoms {
            *a = num_3d_stereo_atoms;
        }
    }
    if let Some(b) = want_bonds {
        if max_stereo_bonds < num_stereo_bonds {
            *b = num_stereo_bonds;
        }
    }
    num_3d_stereo_atoms
}

// ---------------------------------------------------------------------------
// ichicans.c / ichimap2.c: parity bookkeeping on sp_ATOM
// ---------------------------------------------------------------------------

pub const AB_MIN_KNOWN_PARITY: i32 = 1;
pub const AB_MAX_KNOWN_PARITY: i32 = 4;
pub const AB_INV_PARITY_BITS: i32 = AB_PARITY_ODD ^ AB_PARITY_EVEN;
pub const BITS_PARITY: i32 = 0x07;

/// `PARITY_VAL`
pub fn parity_val(x: i32) -> i32 {
    x & BITS_PARITY
}
/// `PARITY_KNOWN`
pub fn parity_known(x: i32) -> bool {
    (AB_MIN_KNOWN_PARITY..=AB_MAX_KNOWN_PARITY).contains(&parity_val(x))
}
/// `PARITY_WELL_DEF`
pub fn parity_well_def(x: i32) -> bool {
    atom_parity_well_def(parity_val(x))
}

/// `find_atoms_with_parity`: depth first search for an atom with a parity
pub fn find_atoms_with_parity(
    at: &[SpAtom],
    visited: &mut [i8],
    from_atom: i32,
    cur_atom: usize,
) -> bool {
    if visited[cur_atom] != 0 {
        return false;
    }
    if at[cur_atom].parity != 0 {
        return true;
    }
    visited[cur_atom] = 1;
    for i in 0..at[cur_atom].valence as usize {
        let next_atom = at[cur_atom].neighbor[i] as usize;
        if next_atom as i32 != from_atom
            && find_atoms_with_parity(at, visited, cur_atom as i32, next_atom)
        {
            return true;
        }
    }
    false
}

/// `SetHalfStereoBondIllDefPariy`
pub fn set_half_stereo_bond_ill_def_pariy(
    at: &mut [SpAtom],
    jn: usize,
    k1: usize,
    new_parity: i32,
) -> bool {
    if k1 < MAX_NUM_STEREO_BOND_NEIGH as usize && at[jn].stereo_bond_neighbor[k1] != 0 {
        let sbp = i32::from(at[jn].stereo_bond_parity[k1]);
        let parity = sbp ^ parity_val(sbp);
        at[jn].stereo_bond_parity[k1] = (parity | parity_val(new_parity)) as i8;
        at[jn].parity = parity_val(new_parity) as i8;
        return true;
    }
    false
}

/// `RemoveHalfStereoBond`
pub fn remove_half_stereo_bond(at: &mut [SpAtom], jn: usize, k1: usize) -> bool {
    let n = MAX_NUM_STEREO_BOND_NEIGH as usize;
    if k1 < n && at[jn].stereo_bond_neighbor[k1] != 0 {
        let a = &mut at[jn];
        for k2 in k1..n - 1 {
            a.stereo_bond_neighbor[k2] = a.stereo_bond_neighbor[k2 + 1];
            a.stereo_bond_ord[k2] = a.stereo_bond_ord[k2 + 1];
            a.stereo_bond_z_prod[k2] = a.stereo_bond_z_prod[k2 + 1];
            a.stereo_bond_parity[k2] = a.stereo_bond_parity[k2 + 1];
        }
        a.stereo_bond_neighbor[n - 1] = 0;
        a.stereo_bond_ord[n - 1] = 0;
        a.stereo_bond_z_prod[n - 1] = 0;
        a.stereo_bond_parity[n - 1] = 0;
        if a.stereo_bond_neighbor[0] == 0 {
            a.parity = 0;
            a.stereo_atom_parity = 0;
            a.final_parity = 0;
        }
        return true;
    }
    false
}

fn find_opposite_half_bond(at: &[SpAtom], jc: usize, jn: usize) -> Option<usize> {
    for k1 in 0..MAX_NUM_STEREO_BOND_NEIGH as usize {
        let kn = at[jn].stereo_bond_neighbor[k1];
        if kn == 0 {
            break;
        }
        if kn as usize - 1 == jc {
            return Some(k1);
        }
    }
    None
}

/// `SetOneStereoBondIllDefParity`
pub fn set_one_stereo_bond_ill_def_parity(
    at: &mut [SpAtom],
    jc: usize,
    k: usize,
    new_parity: i32,
) -> bool {
    let jn = at[jc].stereo_bond_neighbor[k] as usize - 1;
    let mut ret = false;
    if let Some(k1) = find_opposite_half_bond(at, jc, jn) {
        ret = set_half_stereo_bond_ill_def_pariy(at, jn, k1, new_parity);
    }
    if ret {
        ret = set_half_stereo_bond_ill_def_pariy(at, jc, k, new_parity);
    }
    ret
}

/// `RemoveOneStereoBond`
pub fn remove_one_stereo_bond(at: &mut [SpAtom], jc: usize, k: usize) -> bool {
    let jn = at[jc].stereo_bond_neighbor[k] as usize - 1;
    let mut ret = false;
    if let Some(k1) = find_opposite_half_bond(at, jc, jn) {
        ret = remove_half_stereo_bond(at, jn, k1);
    }
    if ret {
        ret = remove_half_stereo_bond(at, jc, k);
    }
    ret
}

/// `RemoveOneStereoCenter`
pub fn remove_one_stereo_center(at: &mut [SpAtom], jc: usize) -> bool {
    if at[jc].parity != 0 {
        at[jc].parity = 0;
        at[jc].stereo_atom_parity = 0;
        at[jc].final_parity = 0;
        return true;
    }
    false
}

/// `GetStereoNeighborPos`: the position of the stereo bond `iat1`-`iat2` in `iat1`, or -1
pub fn get_stereo_neighbor_pos(at: &[SpAtom], iat1: usize, iat2: usize) -> i32 {
    let s_neigh = (iat2 + 1) as u16;
    for k1 in 0..MAX_NUM_STEREO_BONDS {
        let s = at[iat1].stereo_bond_neighbor[k1];
        if s == 0 {
            break;
        }
        if s == s_neigh {
            return k1 as i32;
        }
    }
    -1
}

/// `HalfStereoBondParity` (ichimap2.c): 0 if no parity can be found, 1 odd, 2 even
pub fn half_stereo_bond_parity_ranks(
    at: &[SpAtom],
    at_no1: usize,
    i_sb_neigh: i32,
    n_rank: &[u16],
) -> i32 {
    use super::ichierr::CT_STEREOBOND_ERROR;
    let a = &at[at_no1];
    let at1_parity = i32::from(a.parity);
    if a.valence as i32 > MAX_NUM_STEREO_BOND_NEIGH || at1_parity <= 0 {
        return 0;
    }
    if !parity_well_def(at1_parity) {
        return if parity_known(at1_parity) {
            at1_parity
        } else {
            -at1_parity
        };
    }
    if !(0..MAX_NUM_STEREO_BOND_NEIGH).contains(&i_sb_neigh) {
        return CT_STEREOBOND_ERROR;
    }
    let isb = i_sb_neigh as usize;
    for i in 0..=isb {
        if a.stereo_bond_neighbor[i] == 0 {
            return CT_STEREOBOND_ERROR;
        }
    }
    let at_no2 = a.neighbor[a.stereo_bond_ord[isb] as usize] as usize;
    let mut r = [0u16; MAX_NUM_STEREO_BOND_NEIGH as usize];
    let mut j = 0usize;
    let mut i_neigh: i32 = -1;
    for i in 0..a.valence as usize {
        let k = a.neighbor[i] as usize;
        if k == at_no2 {
            i_neigh = i as i32;
        } else {
            r[j] = n_rank[k];
            j += 1;
        }
    }
    if i_neigh < 0 || i_neigh != i32::from(a.stereo_bond_ord[isb]) {
        return CT_STEREOBOND_ERROR;
    }
    if (j > 0 && r[0] == 0) || (j > 1 && r[1] == 0) {
        return 0;
    }
    if j == 2 && r[0] == r[1] {
        AB_PARITY_CALC
    } else {
        2 - (at1_parity + i_neigh + i32::from(r[1] < r[0])) % 2
    }
}

/// `GetStereoBondParity`: -1 on an error, otherwise the parity of the stereo bond `i`-`n`
pub fn get_stereo_bond_parity(at: &[SpAtom], i: usize, n: usize, n_rank: &[u16]) -> i32 {
    use super::ichierr::CT_STEREOBOND_ERROR;
    if at[i].stereo_bond_neighbor[0] == 0 {
        return -1;
    }
    let mut k1 = None;
    for k in 0..MAX_NUM_STEREO_BONDS {
        let s = at[i].stereo_bond_neighbor[k];
        if s == 0 {
            break;
        }
        if s as usize - 1 == n {
            k1 = Some(k);
            break;
        }
    }
    let Some(k1) = k1 else { return -1 };
    if parity_known(i32::from(at[i].stereo_bond_parity[k1])) {
        return parity_val(i32::from(at[i].stereo_bond_parity[k1]));
    }
    let mut k2 = None;
    for k in 0..MAX_NUM_STEREO_BONDS {
        let s = at[n].stereo_bond_neighbor[k];
        if s == 0 {
            break;
        }
        if s as usize - 1 == i {
            k2 = Some(k);
            break;
        }
    }
    let Some(k2) = k2 else { return -1 };
    if atom_parity_well_def(i32::from(at[i].parity))
        && atom_parity_well_def(i32::from(at[n].parity))
        && MIN_DOT_PROD <= i32::from(at[i].stereo_bond_z_prod[k1]).abs()
    {
        let half_parity1 = half_stereo_bond_parity_ranks(at, i, k1 as i32, n_rank);
        let half_parity2 = half_stereo_bond_parity_ranks(at, n, k2 as i32, n_rank);
        if half_parity1 == 0 || half_parity2 == 0 {
            return 0;
        }
        if atom_parity_well_def(half_parity1) && atom_parity_well_def(half_parity2) {
            2 - (half_parity1 + half_parity2 + i32::from(at[i].stereo_bond_z_prod[k1] < 0)) % 2
        } else {
            CT_STEREOBOND_ERROR
        }
    } else {
        let mut parity = i32::from(at[i].parity).max(i32::from(at[n].parity));
        if parity != AB_PARITY_NONE {
            parity = AB_PARITY_UNDF;
        }
        parity
    }
}

/// `GetPermutationParity`: parity of the permutation sorting the canonical
/// ranks of the neighbors (without `n_avoid_neighbor`); -1 on an error
pub fn get_permutation_parity(at: &SpAtom, n_avoid_neighbor: u16, n_canon_rank: &[u16]) -> i32 {
    if at.valence as usize > MAX_NUM_STEREO_ATOM_NEIGH {
        return -1;
    }
    let mut ranks = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
    let mut k = 0;
    for j in 0..at.valence as usize {
        if at.neighbor[j] != n_avoid_neighbor {
            ranks[k] = n_canon_rank[at.neighbor[j] as usize];
            k += 1;
        }
    }
    if k == 0 {
        return 2;
    }
    let parity = insertions_sort(&mut ranks[..k], super::ichisort::comp_at_rank);
    if ranks[0] != 0 { 2 - parity % 2 } else { 0 }
}

/// `GetStereoCenterParity`
pub fn get_stereo_center_parity(at: &[SpAtom], i: usize, n_rank: &[u16]) -> i32 {
    if at[i].parity == 0 {
        return 0;
    }
    if at[i].stereo_bond_neighbor[0] != 0 {
        return -1;
    }
    if atom_parity_well_def(i32::from(at[i].parity)) {
        let num_neigh = at[i].valence as usize;
        let mut order = [0u16; super::inpdef::MAXVAL];
        for k in 0..num_neigh {
            if n_rank[at[i].neighbor[k] as usize] == 0 {
                return 0;
            }
            order[k] = k as u16;
        }
        // CT_NEIGH_INCREASE: neighbors ranks increase
        let num_trans = insertions_sort(&mut order[..num_neigh], |a, b| {
            i32::from(n_rank[at[i].neighbor[*a as usize] as usize])
                - i32::from(n_rank[at[i].neighbor[*b as usize] as usize])
        });
        2 - (i32::from(at[i].parity) + num_trans) % 2
    } else {
        i32::from(at[i].parity)
    }
}

// ---------------------------------------------------------------------------
// ichicans.c: stereo descriptors of the canonical connection table
// ---------------------------------------------------------------------------

/// `AT_STEREO_CARB`
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AtStereoCarb {
    pub at_num: u16,
    pub parity: u8,
}

/// `AT_STEREO_DBLE`
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AtStereoDble {
    pub at_num1: u16,
    pub at_num2: u16,
    pub parity: u8,
}

/// The stereo part of `CANON_STAT`: the (maybe absent) linear stereo CTs
/// and their lengths and capacities
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CanonStatStereo {
    pub carb: Option<Vec<AtStereoCarb>>,
    pub n_len_carb: i32,
    pub n_max_len_carb: i32,
    pub dble: Option<Vec<AtStereoDble>>,
    pub n_len_dble: i32,
    pub n_max_len_dble: i32,
}

pub fn is_allene_chain(x: i32) -> bool {
    ((x & MASK_CUMULENE_LEN) / MULT_STEREOBOND) % 2 != 0
}
pub fn bond_chain_len(x: i32) -> i32 {
    (x & MASK_CUMULENE_LEN) / MULT_STEREOBOND
}

/// `UnmarkNonStereo`: removes the parities of centers and bonds whose
/// constitutionally equivalent neighbors cannot make them stereogenic;
/// returns the number of removed parities
pub fn unmark_non_stereo(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_rank: &[u16],
    n_atom_number: &[u16],
    b_isotopic: bool,
) -> i32 {
    use super::incomdef::NUM_H_ISOTOPES;
    let mut visited = vec![0i8; num_atoms];
    let mut num_removed_parities = 0i32;
    let mut num_no_parity_atoms = 0i32;
    // (declared once in the C code: the values persist between the loops)
    let mut num_in_same_ring_system = 1i32;
    let mut num_with_eq_neigh_in_same_ring_system = 0i32;
    let mut n_neighbor_number = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
    let mut j = 0usize;
    loop {
        let num_removed_parities0 = num_removed_parities;
        let mut i1 = 0usize;
        let mut n_prev_atom_rank = 0u16;
        let mut i = 0usize;
        while i <= num_atoms {
            if i == num_atoms || {
                j = n_atom_number[i] as usize;
                n_prev_atom_rank != n_rank[j]
            } {
                let i2 = i;
                if (i2 - i1) as i32 > num_no_parity_atoms {
                    let mut jc = n_atom_number[i1] as usize;
                    num_no_parity_atoms = 0;
                    let val = at[jc].valence as usize;
                    let num_implicit_h = if at[jc].endpoint != 0 {
                        0
                    } else {
                        i32::from(at[jc].num_h)
                    };
                    if val as i32 + num_implicit_h > MAX_NUM_STEREO_ATOM_NEIGH as i32 {
                        i += 1;
                        continue;
                    }
                    for k in 0..val {
                        n_neighbor_number[k] = k as u16;
                    }
                    for ic in i1..i2 {
                        jc = n_atom_number[ic] as usize;
                        let nbrs = at[jc].neighbor;
                        insertions_sort(&mut n_neighbor_number[..val], |a, b| {
                            i32::from(n_rank[nbrs[*a as usize] as usize])
                                - i32::from(n_rank[nbrs[*b as usize] as usize])
                        });
                        let mut num_neighbors_with_parity = -1i32;
                        let mut k1 = 0usize;
                        let mut n_prev_neigh_rank = 0u16;
                        let mut k = 0usize;
                        while k <= val {
                            if k == val
                                || n_prev_neigh_rank
                                    != n_rank[nbrs[n_neighbor_number[k] as usize] as usize]
                            {
                                let k2 = k;
                                if k2 - k1 > 1 {
                                    num_in_same_ring_system = 0;
                                    let mut n_ring_system = 0u16;
                                    for kn in k1..k2 {
                                        let cur = at[nbrs[n_neighbor_number[kn] as usize] as usize]
                                            .n_ring_system;
                                        if n_ring_system == 0 {
                                            n_ring_system = cur;
                                        } else {
                                            num_in_same_ring_system +=
                                                i32::from(n_ring_system == cur);
                                        }
                                    }
                                    num_neighbors_with_parity = 0;
                                    for kn in k1..k2 {
                                        visited.iter_mut().for_each(|v| *v = 0);
                                        visited[jc] = 1;
                                        num_neighbors_with_parity +=
                                            i32::from(find_atoms_with_parity(
                                                at,
                                                &mut visited,
                                                jc as i32,
                                                nbrs[n_neighbor_number[kn] as usize] as usize,
                                            ));
                                    }
                                }
                                if num_neighbors_with_parity == 0 && num_in_same_ring_system == 0 {
                                    break;
                                }
                                if k + 1 < val {
                                    k1 = k;
                                    n_prev_neigh_rank =
                                        n_rank[nbrs[n_neighbor_number[k] as usize] as usize];
                                } else {
                                    break;
                                }
                            }
                            k += 1;
                        }
                        if num_implicit_h > 1 {
                            let a = &at[jc];
                            if (b_isotopic
                                && (a.num_iso_h[0] > 1 || a.num_iso_h[1] > 1 || a.num_iso_h[2] > 1))
                                || num_implicit_h > NUM_H_ISOTOPES as i32
                                || !b_isotopic
                            {
                                num_neighbors_with_parity = 0;
                            }
                        }
                        num_no_parity_atoms += i32::from(num_neighbors_with_parity == 0);
                        num_with_eq_neigh_in_same_ring_system +=
                            i32::from(num_in_same_ring_system != 0);
                    }
                    let cnt = (i2 - i1) as i32;
                    if num_no_parity_atoms == cnt && num_with_eq_neigh_in_same_ring_system != cnt {
                        for ic in i1..i2 {
                            jc = n_atom_number[ic] as usize;
                            at[jc].parity = 0;
                            at[jc].stereo_atom_parity = 0;
                            at[jc].final_parity = 0;
                            at[jc].b_has_stereo_or_equ_to_stereo = 0;
                            let mut k = 0;
                            while k < MAX_NUM_STEREO_BOND_NEIGH as usize
                                && at[jc].stereo_bond_neighbor[k] != 0
                            {
                                let jn = at[jc].stereo_bond_neighbor[k] as usize - 1;
                                let mut k1 = 0;
                                while k1 < MAX_NUM_STEREO_BOND_NEIGH as usize
                                    && at[jn].stereo_bond_neighbor[k1] != 0
                                {
                                    if at[jn].stereo_bond_neighbor[k1] as usize - 1 == jc {
                                        remove_half_stereo_bond(at, jn, k1);
                                        break;
                                    }
                                    k1 += 1;
                                }
                                at[jc].stereo_bond_neighbor[k] = 0;
                                at[jc].stereo_bond_ord[k] = 0;
                                at[jc].stereo_bond_z_prod[k] = 0;
                                at[jc].stereo_bond_parity[k] = 0;
                                k += 1;
                            }
                        }
                        num_removed_parities += num_no_parity_atoms;
                    }
                }
                if i < num_atoms {
                    n_prev_atom_rank = n_rank[j];
                    i1 = i;
                }
                num_no_parity_atoms = 0;
            }
            num_no_parity_atoms += i32::from(i < num_atoms && at[j].parity == 0);
            i += 1;
        }
        if num_removed_parities == num_removed_parities0 {
            break;
        }
    }
    num_removed_parities
}

/// `FillSingleStereoDescriptors`: adds the stereo descriptors of atom `i` to the linear CTs
#[allow(clippy::too_many_arguments)]
pub fn fill_single_stereo_descriptors(
    at: &[SpAtom],
    i: usize,
    mut num_trans: i32,
    n_rank: &[u16],
    carb: Option<&mut [AtStereoCarb]>,
    n_stereo_carb_len: &mut i32,
    n_max_stereo_carb_len: i32,
    dble: Option<&mut [AtStereoDble]>,
    n_stereo_dble_len: &mut i32,
    n_max_stereo_dble_len: i32,
    b_allene: i32,
) -> i32 {
    use super::ichierr::{CT_OVERFLOW, CT_STEREOBOND_ERROR};
    if dble.is_none() && carb.is_none() {
        return 0;
    }
    if at[i].parity != 0 || at[i].stereo_bond_neighbor[0] != 0 {
        let rank = n_rank[i];
        let mut num_allene = 0;
        if atom_parity_well_def(i32::from(at[i].parity)) && num_trans < 0 {
            let num_neigh = at[i].valence as usize;
            let mut order = [0u16; super::inpdef::MAXVAL];
            for k in 0..num_neigh {
                order[k] = k as u16;
            }
            num_trans = insertions_sort(&mut order[..num_neigh], |a, b| {
                i32::from(n_rank[at[i].neighbor[*a as usize] as usize])
                    - i32::from(n_rank[at[i].neighbor[*b as usize] as usize])
            });
        }
        if let (Some(dble), true) = (dble, at[i].stereo_bond_neighbor[0] != 0) {
            let mut n_stereo_neigh_number = [0u16; MAX_NUM_STEREO_BONDS];
            let mut n_stereo_neigh = [0u16; MAX_NUM_STEREO_BONDS];
            let mut num_stereo = 0usize;
            while num_stereo < MAX_NUM_STEREO_BONDS && at[i].stereo_bond_neighbor[num_stereo] != 0 {
                n_stereo_neigh_number[num_stereo] = num_stereo as u16;
                n_stereo_neigh[num_stereo] = at[i].stereo_bond_neighbor[num_stereo] - 1;
                num_allene += i32::from(is_allene_chain(i32::from(
                    at[i].stereo_bond_parity[num_stereo],
                )));
                num_stereo += 1;
            }
            if (b_allene > 0 && num_allene == 0) || (b_allene == 0 && num_allene != 0) {
                return 0;
            }
            insertions_sort(&mut n_stereo_neigh_number[..num_stereo], |a, b| {
                i32::from(n_rank[n_stereo_neigh[*a as usize] as usize])
                    - i32::from(n_rank[n_stereo_neigh[*b as usize] as usize])
            });
            for k in 0..num_stereo {
                let stereo_neigh_ord = n_stereo_neigh_number[k] as usize;
                let stereo_neigh = n_stereo_neigh[stereo_neigh_ord] as usize;
                let r_neigh = n_rank[stereo_neigh];
                if r_neigh < rank {
                    let stereo_bond_parity =
                        parity_val(i32::from(at[i].stereo_bond_parity[stereo_neigh_ord]));
                    if stereo_bond_parity == AB_PARITY_NONE {
                        continue;
                    }
                    let mut parity;
                    if (AB_MIN_KNOWN_PARITY..=AB_MAX_KNOWN_PARITY).contains(&stereo_bond_parity) {
                        parity = stereo_bond_parity;
                    } else if atom_parity_well_def(i32::from(at[i].parity))
                        && atom_parity_well_def(i32::from(at[stereo_neigh].parity))
                        && MIN_DOT_PROD
                            <= i32::from(at[i].stereo_bond_z_prod[stereo_neigh_ord]).abs()
                    {
                        let mut stereo_neigh_ord2: i32 = -1;
                        for j in 0..MAX_NUM_STEREO_BONDS {
                            let nn = at[stereo_neigh].stereo_bond_neighbor[j] as usize;
                            if nn == 0 {
                                break;
                            }
                            if i + 1 == nn {
                                stereo_neigh_ord2 = j as i32;
                                break;
                            }
                        }
                        if stereo_neigh_ord2 >= 0 {
                            let hp1 = half_stereo_bond_parity_ranks(
                                at,
                                i,
                                stereo_neigh_ord as i32,
                                n_rank,
                            );
                            let hp2 = half_stereo_bond_parity_ranks(
                                at,
                                stereo_neigh,
                                stereo_neigh_ord2,
                                n_rank,
                            );
                            if atom_parity_well_def(hp1) && atom_parity_well_def(hp2) {
                                parity = 2
                                    - (hp1
                                        + hp2
                                        + i32::from(
                                            at[i].stereo_bond_z_prod[stereo_neigh_ord] < 0,
                                        ))
                                        % 2;
                            } else {
                                return CT_STEREOBOND_ERROR;
                            }
                        } else {
                            return CT_STEREOBOND_ERROR;
                        }
                    } else {
                        parity = i32::from(at[i].parity).max(i32::from(at[stereo_neigh].parity));
                        if parity == AB_PARITY_NONE {
                            continue;
                        }
                        if atom_parity_well_def(parity) {
                            parity = AB_PARITY_UNDF;
                        }
                    }
                    if *n_stereo_dble_len >= n_max_stereo_dble_len {
                        return CT_OVERFLOW;
                    }
                    dble[*n_stereo_dble_len as usize] = AtStereoDble {
                        at_num1: rank,
                        at_num2: r_neigh,
                        parity: parity as u8,
                    };
                    *n_stereo_dble_len += 1;
                }
            }
        }
        if b_allene > 0 {
            return 0;
        }
        if let (Some(carb), true) = (carb, at[i].stereo_bond_neighbor[0] == 0) {
            if *n_stereo_carb_len >= n_max_stereo_carb_len {
                return CT_OVERFLOW;
            }
            let parity = if atom_parity_well_def(i32::from(at[i].parity)) {
                2 - (i32::from(at[i].parity) + num_trans) % 2
            } else {
                i32::from(at[i].parity)
            };
            carb[*n_stereo_carb_len as usize] = AtStereoCarb {
                at_num: rank,
                parity: parity as u8,
            };
            *n_stereo_carb_len += 1;
        }
    }
    0
}

/// `SwitchAtomStereoAndIsotopicStereo`
pub fn switch_atom_stereo_and_isotopic_stereo(
    at: &mut [SpAtom],
    num_atoms: usize,
    b_switched: &mut bool,
) {
    for a in at.iter_mut().take(num_atoms) {
        std::mem::swap(&mut a.parity, &mut a.parity2);
        std::mem::swap(&mut a.final_parity, &mut a.final_parity2);
        std::mem::swap(&mut a.stereo_atom_parity, &mut a.stereo_atom_parity2);
        std::mem::swap(
            &mut a.b_has_stereo_or_equ_to_stereo,
            &mut a.b_has_stereo_or_equ_to_stereo2,
        );
        std::mem::swap(&mut a.stereo_bond_neighbor, &mut a.stereo_bond_neighbor2);
        std::mem::swap(&mut a.stereo_bond_ord, &mut a.stereo_bond_ord2);
        std::mem::swap(&mut a.stereo_bond_z_prod, &mut a.stereo_bond_z_prod2);
        std::mem::swap(&mut a.stereo_bond_parity, &mut a.stereo_bond_parity2);
    }
    *b_switched = !*b_switched;
}

/// `FillAllStereoDescriptors`
pub fn fill_all_stereo_descriptors(
    at: &[SpAtom],
    num_atoms: usize,
    n_canon_rank: &[u16],
    n_atom_number_canon: &[u16],
    cs: &mut CanonStatStereo,
) -> i32 {
    cs.n_len_carb = 0;
    cs.n_len_dble = 0;
    let mut ret = 0;
    for b_allene in 0..2 {
        let mut i = 0;
        while ret == 0 && i < num_atoms {
            ret = fill_single_stereo_descriptors(
                at,
                n_atom_number_canon[i] as usize,
                -1,
                n_canon_rank,
                cs.carb.as_deref_mut(),
                &mut cs.n_len_carb,
                cs.n_max_len_carb,
                cs.dble.as_deref_mut(),
                &mut cs.n_len_dble,
                cs.n_max_len_dble,
                b_allene,
            );
            i += 1;
        }
    }
    ret
}

/// `InvertStereo`: inverts the known parities in `at` and (optionally) in the
/// linear stereo CTs; `n_atom_number_canon` is filled in
pub fn invert_stereo(
    at: &mut [SpAtom],
    num_at_tg: usize,
    n_canon_rank: &[u16],
    n_atom_number_canon: &mut [u16],
    cs: &mut CanonStatStereo,
    b_invert_linear_ct_stereo: bool,
) -> i32 {
    use super::ichierr::CT_STEREOCOUNT_ERR;
    let mut num_changes = 0;
    for i in 0..num_at_tg {
        n_atom_number_canon[n_canon_rank[i] as usize - 1] = i as u16;
    }
    for i in 0..cs.n_len_carb as usize {
        let c = cs.carb.as_ref().unwrap()[i];
        if atom_parity_well_def(i32::from(c.parity)) {
            let j = n_atom_number_canon[c.at_num as usize - 1] as usize;
            if parity_well_def(i32::from(at[j].parity)) {
                at[j].parity ^= AB_INV_PARITY_BITS as i8;
            } else {
                return CT_STEREOCOUNT_ERR;
            }
            if b_invert_linear_ct_stereo {
                // FIX_STEREOCOUNT_ERR: deliberately worse
                cs.carb.as_mut().unwrap()[i].parity = AB_PARITY_EVEN as u8;
            }
            num_changes += 1;
            if parity_well_def(i32::from(at[j].stereo_atom_parity)) {
                at[j].stereo_atom_parity ^= AB_INV_PARITY_BITS as i8;
            }
            if parity_well_def(i32::from(at[j].final_parity)) {
                at[j].final_parity ^= AB_INV_PARITY_BITS as i8;
            }
        }
    }
    for i in 0..cs.n_len_dble as usize {
        let d = cs.dble.as_ref().unwrap()[i];
        if atom_parity_well_def(i32::from(d.parity)) {
            let j1 = n_atom_number_canon[d.at_num1 as usize - 1] as usize;
            let cumulene_len = bond_chain_len(i32::from(at[j1].stereo_bond_parity[0]));
            if cumulene_len % 2 != 0 {
                let j2 = n_atom_number_canon[d.at_num2 as usize - 1] as usize;
                if at[j1].stereo_bond_neighbor[1] != 0 || at[j2].stereo_bond_neighbor[1] != 0 {
                    return CT_STEREOCOUNT_ERR;
                }
                if cumulene_len != bond_chain_len(i32::from(at[j2].stereo_bond_parity[0]))
                    || j1 + 1 != at[j2].stereo_bond_neighbor[0] as usize
                    || j2 + 1 != at[j1].stereo_bond_neighbor[0] as usize
                {
                    return CT_STEREOCOUNT_ERR;
                }
                if parity_well_def(i32::from(at[j1].parity))
                    && parity_well_def(i32::from(at[j2].parity))
                {
                    let j = j1.min(j2);
                    at[j].parity ^= AB_INV_PARITY_BITS as i8;
                } else {
                    return CT_STEREOCOUNT_ERR;
                }
                if b_invert_linear_ct_stereo {
                    cs.dble.as_mut().unwrap()[i].parity ^= AB_INV_PARITY_BITS as u8;
                }
                num_changes += 1;
                if parity_well_def(i32::from(at[j1].stereo_bond_parity[0])) {
                    at[j1].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS as i8;
                }
                if parity_well_def(i32::from(at[j2].stereo_bond_parity[0])) {
                    at[j2].stereo_bond_parity[0] ^= AB_INV_PARITY_BITS as i8;
                }
            }
        }
    }
    num_changes
}

// ---------------------------------------------------------------------------
// ichicans.c: stereo bond parities known in advance
// ---------------------------------------------------------------------------

pub const KNOWN_PARITIES_EQL: i32 = 0x40;
const MAX_ATOMS_I: i32 = 32766;

fn all_but_parity(x: i32) -> i32 {
    x & !BITS_PARITY
}

/// the stereo bond `i1`-`i2` as seen from `i2`: the position of `i1` in the list of `i2`
fn find_stereo_neighbor(at: &[SpAtom], i2: usize, i1: usize) -> Option<usize> {
    for n2 in 0..MAX_NUM_STEREO_BONDS {
        let m = at[i2].stereo_bond_neighbor[n2] as usize;
        if m == 0 {
            return None;
        }
        if m - 1 == i1 {
            return Some(n2);
        }
    }
    None
}

/// `SetKnownStereoBondParities`: calculates the parities of the stereo
/// bonds that do not depend on breaking ties; returns the number of bonds set
pub fn set_known_stereo_bond_parities(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_canon_rank: &[u16],
    n_rank: &[u16],
    n_atom_number: &[u16],
) -> i32 {
    use super::ichierr::CT_STEREOCOUNT_ERR;
    use super::ichisort::comp_at_rank;
    let mut num_set = 0;
    for i1 in 0..num_atoms {
        if at[i1].parity == 0 || at[i1].stereo_bond_neighbor[0] == 0 {
            continue;
        }
        if !parity_well_def(i32::from(at[i1].parity)) {
            continue;
        }
        let n_atom_rank1 = n_rank[i1];
        let i_max1 = i32::from(n_atom_rank1) - 1;
        let num_neigh1 = at[i1].valence as usize;
        let mut n1 = 0;
        while n1 < MAX_NUM_STEREO_BONDS && at[i1].stereo_bond_neighbor[n1] != 0 {
            let i2 = at[i1].stereo_bond_neighbor[n1] as usize - 1;
            let Some(n2) = find_stereo_neighbor(at, i2, i1) else {
                return CT_STEREOCOUNT_ERR;
            };
            if at[i1].stereo_bond_parity[n1] != at[i2].stereo_bond_parity[n2] {
                return CT_STEREOCOUNT_ERR;
            }
            n1 += 1;
            if i1 < i2 {
                continue;
            }
            let sbp1 = i32::from(at[i1].stereo_bond_parity[n1 - 1]);
            if parity_known(sbp1) || parity_val(sbp1) == 0 {
                continue;
            }
            if !parity_well_def(i32::from(at[i1].parity))
                || !parity_well_def(i32::from(at[i2].parity))
            {
                continue;
            }
            if parity_val(sbp1) != AB_PARITY_CALC {
                continue;
            }
            let n1 = n1 - 1;
            let cumulene_len = bond_chain_len(sbp1);
            let n_atom_rank2 = n_rank[i2];
            let n_atom1_neigh_rank =
                n_rank[at[i1].neighbor[at[i1].stereo_bond_ord[n1] as usize] as usize];
            let num_neigh2 = at[i2].valence as usize;

            let mut n_neigh_rank1 = [0u16; MAX_NUM_STEREO_BONDS];
            let mut n_neigh_rank2 = [0u16; MAX_NUM_STEREO_BONDS];
            let mut n_neigh_canon_rank1 = [0u16; MAX_NUM_STEREO_BONDS];
            let mut n_neigh_canon_rank2 = [0u16; MAX_NUM_STEREO_BONDS];

            let k = at[i1].stereo_bond_ord[n1] as usize;
            let mut trans_i1 = 0;
            let mut j = 0;
            for i in 0..num_neigh1 {
                if i != k {
                    n_neigh_rank1[j] = n_rank[at[i1].neighbor[i] as usize];
                    j += 1;
                }
            }
            if j == 2 {
                if n_neigh_rank1[0] == n_neigh_rank1[1] {
                    continue;
                }
                trans_i1 = insertions_sort(&mut n_neigh_rank1[..j], comp_at_rank);
            }
            let k = at[i2].stereo_bond_ord[n2] as usize;
            let mut trans_i2 = 0;
            j = 0;
            for i in 0..num_neigh2 {
                if i != k {
                    n_neigh_rank2[j] = n_rank[at[i2].neighbor[i] as usize];
                    j += 1;
                }
            }
            if j == 2 {
                if n_neigh_rank2[0] == n_neigh_rank2[1] {
                    continue;
                }
                trans_i2 = insertions_sort(&mut n_neigh_rank2[..j], comp_at_rank);
            }

            let mut prev_trans: i32 = -1;
            let mut trans_k: i32 = -4;
            let mut j1: i32 = 0;
            'outer: while j1 <= i_max1 {
                let k1 = n_atom_number[(i_max1 - j1) as usize] as usize;
                if n_atom_rank1 != n_rank[k1] {
                    break;
                }
                if at[k1].valence as usize != num_neigh1 {
                    return CT_STEREOCOUNT_ERR;
                }
                for m1 in 0..num_neigh1 {
                    let mut k2 = at[k1].neighbor[m1] as usize;
                    if n_atom1_neigh_rank != n_rank[k2] {
                        continue;
                    }
                    let mut prev = k1;
                    if cumulene_len != 0 {
                        let mut next = at[k1].neighbor[m1] as usize;
                        let mut len = 0;
                        while len < cumulene_len {
                            if at[next].valence == 2 && at[next].num_h == 0 {
                                let jj = usize::from(at[next].neighbor[0] as usize == prev);
                                prev = next;
                                next = at[next].neighbor[jj] as usize;
                            } else {
                                break;
                            }
                            len += 1;
                        }
                        if len != cumulene_len || n_atom_rank2 != n_rank[next] {
                            continue;
                        }
                        k2 = next;
                    }
                    if at[k2].valence as usize != num_neigh2 {
                        return CT_STEREOCOUNT_ERR;
                    }
                    let mut jc = 0;
                    for n in 0..num_neigh1 {
                        if n != m1 {
                            let i = at[k1].neighbor[n] as usize;
                            for m in 0..num_neigh1 - 1 {
                                if n_rank[i] == n_neigh_rank1[m] {
                                    n_neigh_canon_rank1[m] = n_canon_rank[i];
                                    jc += 1;
                                    break;
                                }
                            }
                        }
                    }
                    if jc != num_neigh1 - 1 {
                        return CT_STEREOCOUNT_ERR;
                    }
                    let trans_k1 = if jc == 2 {
                        insertions_sort(&mut n_neigh_canon_rank1[..jc], comp_at_rank)
                    } else {
                        0
                    };
                    jc = 0;
                    let mut m2: i32 = -1;
                    for n in 0..num_neigh2 {
                        let i = at[k2].neighbor[n] as usize;
                        if i == prev {
                            m2 = n as i32;
                        } else {
                            for m in 0..num_neigh2 - 1 {
                                if n_rank[i] == n_neigh_rank2[m] {
                                    n_neigh_canon_rank2[m] = n_canon_rank[i];
                                    jc += 1;
                                    break;
                                }
                            }
                        }
                    }
                    if jc != num_neigh2 - 1 || m2 < 0 {
                        return CT_STEREOCOUNT_ERR;
                    }
                    let trans_k2 = if jc == 2 {
                        insertions_sort(&mut n_neigh_canon_rank2[..jc], comp_at_rank)
                    } else {
                        0
                    };
                    trans_k = (trans_k1 + trans_k2) % 2;
                    if prev_trans < 0 {
                        prev_trans = trans_k;
                    } else if prev_trans != trans_k {
                        break;
                    }
                }
                if prev_trans >= 0 && prev_trans != trans_k {
                    break 'outer;
                }
                j1 += 1;
            }
            if prev_trans == trans_k {
                let mut parity = (i32::from(at[i1].parity)
                    + i32::from(at[i2].parity)
                    + prev_trans
                    + trans_i1
                    + trans_i2
                    + num_neigh1 as i32
                    + num_neigh2 as i32
                    + i32::from(at[i1].stereo_bond_ord[n1])
                    + i32::from(at[i2].stereo_bond_ord[n2]))
                    % 2;
                let z_prod = i32::from(at[i1].stereo_bond_z_prod[n1]);
                if MIN_DOT_PROD > z_prod.abs() {
                    parity = AB_PARITY_UNDF;
                } else {
                    parity = if z_prod > 0 { 2 - parity } else { 1 + parity };
                }
                at[i1].stereo_bond_parity[n1] =
                    (all_but_parity(i32::from(at[i1].stereo_bond_parity[n1])) | parity) as i8;
                at[i2].stereo_bond_parity[n2] =
                    (all_but_parity(i32::from(at[i2].stereo_bond_parity[n2])) | parity) as i8;
                num_set += 1;
            }
        }
    }
    num_set
}

/// `MarkKnownEqualStereoBondParities`
pub fn mark_known_equal_stereo_bond_parities(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_rank: &[u16],
    n_atom_number: &[u16],
) -> i32 {
    use super::ichierr::CT_STEREOCOUNT_ERR;
    let mut num_set = 0;
    for i1 in 0..num_atoms {
        if at[i1].parity == 0 || at[i1].stereo_bond_neighbor[0] == 0 {
            continue;
        }
        let n_atom_rank1 = n_rank[i1];
        let i_max1 = i32::from(n_atom_rank1) - 1;
        let num_neigh1 = at[i1].valence as usize;
        let mut n1 = 0;
        while n1 < MAX_NUM_STEREO_BONDS && at[i1].stereo_bond_neighbor[n1] != 0 {
            let i2 = at[i1].stereo_bond_neighbor[n1] as usize - 1;
            let cur_n1 = n1;
            n1 += 1;
            let n_atom_rank2 = n_rank[i2];
            if n_atom_rank2 < n_atom_rank1 || (n_atom_rank2 == n_atom_rank1 && i1 < i2) {
                continue;
            }
            let mut b_different_parities: i32 = -1;
            if i32::from(at[i1].stereo_bond_parity[cur_n1]) & KNOWN_PARITIES_EQL != 0 {
                continue;
            }
            let Some(n2) = find_stereo_neighbor(at, i2, i1) else {
                return CT_STEREOCOUNT_ERR;
            };
            if at[i1].stereo_bond_parity[cur_n1] != at[i2].stereo_bond_parity[n2] {
                return CT_STEREOCOUNT_ERR;
            }
            let sbp_i1 = i32::from(at[i1].stereo_bond_parity[cur_n1]);
            let stereo_bond_parity = parity_val(sbp_i1);
            let cumulene_len = bond_chain_len(sbp_i1);
            let n_atom1_neigh_rank =
                n_rank[at[i1].neighbor[at[i1].stereo_bond_ord[cur_n1] as usize] as usize];
            let n_atom2_neigh_rank =
                n_rank[at[i2].neighbor[at[i2].stereo_bond_ord[n2] as usize] as usize];
            let num_neigh2 = at[i2].valence as usize;

            let mut j1: i32 = 0;
            while j1 <= i_max1 {
                let k1 = n_atom_number[(i_max1 - j1) as usize] as usize;
                if n_atom_rank1 != n_rank[k1] {
                    break;
                }
                j1 += 1;
                if at[k1].valence as usize != num_neigh1 {
                    return CT_STEREOCOUNT_ERR;
                }
                if at[k1].b_has_stereo_or_equ_to_stereo == 0 {
                    at[k1].b_has_stereo_or_equ_to_stereo = 1;
                }
                for m1 in 0..num_neigh1 {
                    let mut k2 = at[k1].neighbor[m1] as usize;
                    if n_atom1_neigh_rank != n_rank[k2] {
                        continue;
                    }
                    let mut prev = k1;
                    if cumulene_len != 0 {
                        let mut next = at[k1].neighbor[m1] as usize;
                        let mut len = 0;
                        while len < cumulene_len {
                            if at[next].valence == 2 && at[next].num_h == 0 {
                                let jj = usize::from(at[next].neighbor[0] as usize == prev);
                                prev = next;
                                next = at[next].neighbor[jj] as usize;
                            } else {
                                break;
                            }
                            len += 1;
                        }
                        if len != cumulene_len || n_atom_rank2 != n_rank[next] {
                            continue;
                        }
                        if n_atom2_neigh_rank != n_rank[prev] {
                            return CT_STEREOCOUNT_ERR;
                        }
                        k2 = next;
                    }
                    if at[k2].valence as usize != num_neigh2 {
                        return CT_STEREOCOUNT_ERR;
                    }
                    let mut m2 = -1i32;
                    for n in 0..num_neigh2 {
                        if prev == at[k2].neighbor[n] as usize {
                            m2 = n as i32;
                            break;
                        }
                    }
                    if m2 < 0 {
                        return CT_STEREOCOUNT_ERR;
                    }
                    if at[k2].b_has_stereo_or_equ_to_stereo == 0 {
                        at[k2].b_has_stereo_or_equ_to_stereo = 1;
                    }
                    let Some(s1) = find_stereo_neighbor(at, k1, k2) else {
                        b_different_parities = 1;
                        at[k1].b_has_stereo_or_equ_to_stereo = 2;
                        at[k2].b_has_stereo_or_equ_to_stereo = 2;
                        continue;
                    };
                    let Some(s2) = find_stereo_neighbor(at, k2, k1) else {
                        return CT_STEREOCOUNT_ERR;
                    };
                    if at[k1].stereo_bond_parity[s1] != at[k2].stereo_bond_parity[s2] {
                        b_different_parities = 1;
                        continue;
                    }
                    let stereo_bond_parity2 = parity_val(i32::from(at[k1].stereo_bond_parity[s1]));
                    if stereo_bond_parity2 != stereo_bond_parity {
                        b_different_parities = 1;
                        continue;
                    }
                    if stereo_bond_parity2 == stereo_bond_parity && b_different_parities < 0 {
                        b_different_parities = 0;
                    }
                }
            }
            if b_different_parities == 0 && parity_known(stereo_bond_parity) {
                let mut j1: i32 = 0;
                while j1 <= i_max1 {
                    let k1 = n_atom_number[(i_max1 - j1) as usize] as usize;
                    if n_atom_rank1 != n_rank[k1] {
                        break;
                    }
                    j1 += 1;
                    let mut s1 = 0;
                    while s1 < MAX_NUM_STEREO_BONDS && at[k1].stereo_bond_neighbor[s1] != 0 {
                        let k2 = at[k1].stereo_bond_neighbor[s1] as usize - 1;
                        if n_rank[k2] == n_atom_rank2 {
                            let Some(s2) = find_stereo_neighbor(at, k2, k1) else {
                                return CT_STEREOCOUNT_ERR;
                            };
                            let b1 = i32::from(
                                i32::from(at[k1].stereo_bond_parity[s1]) & KNOWN_PARITIES_EQL == 0,
                            );
                            let b2 = i32::from(
                                i32::from(at[k2].stereo_bond_parity[s2]) & KNOWN_PARITIES_EQL == 0,
                            );
                            if b1 + b2 == 2 {
                                at[k1].stereo_bond_parity[s1] |= KNOWN_PARITIES_EQL as i8;
                                at[k2].stereo_bond_parity[s2] |= KNOWN_PARITIES_EQL as i8;
                                num_set += 1;
                            } else if b1 != 0 || b2 != 0 {
                                return CT_STEREOCOUNT_ERR;
                            }
                        }
                        s1 += 1;
                    }
                }
            }
        }
    }
    num_set
}

/// `GetNextNeighborAndRank`
pub fn get_next_neighbor_and_rank(
    at: &[SpAtom],
    cur: u16,
    prev: u16,
    n: &mut u16,
    cr: &mut u16,
    n_canon_rank: &[u16],
) -> bool {
    let mut cr1 = (MAX_ATOMS_I + 1) as u16;
    let mut j1 = cr1;
    for i in 0..at[cur as usize].valence as usize {
        let j = at[cur as usize].neighbor[i];
        if j != prev {
            let crj = n_canon_rank[j as usize];
            if cr1 > crj && crj > *cr {
                cr1 = crj;
                j1 = j;
            }
        }
    }
    if i32::from(cr1) <= MAX_ATOMS_I {
        *cr = cr1;
        *n = j1;
        return true;
    }
    false
}

/// `GetAndCheckNextNeighbors` (one `n_visited` array serves both paths: the
/// only caller passes the same array twice)
#[allow(clippy::too_many_arguments)]
pub fn get_and_check_next_neighbors(
    at: &[SpAtom],
    cur1: u16,
    prev1: u16,
    cur2: u16,
    prev2: u16,
    n1: &mut u16,
    n2: &mut u16,
    n_visited: &[u16],
    n_rank: &[u16],
    n_canon_rank: &[u16],
) -> bool {
    let mut cr1 = if i32::from(*n1) > MAX_ATOMS_I {
        0
    } else {
        n_canon_rank[*n1 as usize]
    };
    let mut cr2 = if i32::from(*n2) > MAX_ATOMS_I {
        0
    } else {
        n_canon_rank[*n2 as usize]
    };
    if !get_next_neighbor_and_rank(at, cur1, prev1, n1, &mut cr1, n_canon_rank)
        || !get_next_neighbor_and_rank(at, cur2, prev2, n2, &mut cr2, n_canon_rank)
        || n_rank[*n1 as usize] != n_rank[*n2 as usize]
        || n_visited[*n1 as usize] != n_visited[*n2 as usize]
    {
        return false;
    }
    let find = |cur: u16, n: u16| -> (usize, bool) {
        let a = &at[cur as usize];
        let mut i = 0;
        let mut k = false;
        while i < MAX_NUM_STEREO_BONDS && a.stereo_bond_neighbor[i] != 0 {
            k = a.neighbor[a.stereo_bond_ord[i] as usize] == n;
            if k {
                break;
            }
            i += 1;
        }
        (i, k)
    };
    let (i1, k1) = find(cur1, *n1);
    let (i2, k2) = find(cur2, *n2);
    if k1 != k2 {
        return false;
    }
    if k1
        && (at[cur1 as usize].stereo_bond_parity[i1] != at[cur2 as usize].stereo_bond_parity[i2]
            || !parity_well_def(i32::from(at[cur1 as usize].stereo_bond_parity[i1])))
    {
        return false;
    }
    true
}

/// `PathsHaveIdenticalKnownParities`
#[allow(clippy::too_many_arguments)]
pub fn paths_have_identical_known_parities(
    at: &[SpAtom],
    prev1: u16,
    cur1: u16,
    prev2: u16,
    cur2: u16,
    n_visited: &mut [u16],
    n_rank: &[u16],
    n_canon_rank: &[u16],
    mut n_length: u16,
) -> u16 {
    n_length += 1;
    n_visited[cur1 as usize] = n_length;
    n_visited[cur2 as usize] = n_length;
    let (a1, a2) = (&at[cur1 as usize], &at[cur2 as usize]);
    if a1.stereo_atom_parity != a2.stereo_atom_parity
        || (a1.stereo_atom_parity != 0 && !parity_well_def(i32::from(a1.stereo_atom_parity)))
    {
        return 0;
    }
    if a1.valence != a2.valence {
        return 0;
    }
    if a1.valence == 1 {
        return n_length;
    }
    let sentinel = (MAX_ATOMS_I + 1) as u16;
    let (mut n1, mut n2) = (sentinel, sentinel);
    for _ in 1..a1.valence {
        if !get_and_check_next_neighbors(
            at,
            cur1,
            prev1,
            cur2,
            prev2,
            &mut n1,
            &mut n2,
            n_visited,
            n_rank,
            n_canon_rank,
        ) {
            return 0;
        }
        if n_visited[n1 as usize] == 0 {
            n_length = paths_have_identical_known_parities(
                at,
                cur1,
                n1,
                cur2,
                n2,
                n_visited,
                n_rank,
                n_canon_rank,
                n_length,
            );
            if n_length == 0 {
                return 0;
            }
        }
    }
    let (mut n1, mut n2) = (sentinel, sentinel);
    for _ in 1..a1.valence {
        if !get_and_check_next_neighbors(
            at,
            cur1,
            prev1,
            cur2,
            prev2,
            &mut n1,
            &mut n2,
            n_visited,
            n_rank,
            n_canon_rank,
        ) {
            return 0;
        }
    }
    n_length
}

/// `RemoveKnownNonStereoBondParities`: removes stereo bonds that are known
/// not to be stereogenic (also from the linear CT); returns how many
pub fn remove_known_non_stereo_bond_parities(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_canon_rank: &[u16],
    n_rank: &[u16],
    cs: &mut CanonStatStereo,
) -> i32 {
    use super::ichierr::{CT_STEREOBOND_ERROR, CT_STEREOCOUNT_ERR};
    let mut ret = 0;
    let mut n_visited: Vec<u16> = Vec::new();
    for i1 in 0..num_atoms {
        if at[i1].valence != 3 || at[i1].stereo_bond_neighbor[0] == 0 {
            continue;
        }
        let mut n1: isize = 0;
        while (n1 as usize) < MAX_NUM_STEREO_BONDS && at[i1].stereo_bond_neighbor[n1 as usize] != 0
        {
            let s2 = at[i1].stereo_bond_neighbor[n1 as usize];
            let sbp = i32::from(at[i1].stereo_bond_parity[n1 as usize]);
            if !(parity_val(sbp) == AB_PARITY_CALC) && parity_well_def(sbp) {
                n1 += 1;
                continue;
            }
            let opposite_atom = (s2 - 1) as usize;
            let m = at[i1].stereo_bond_ord[n1 as usize] as usize;
            let mut neigh = [0u16; 3];
            let mut n = 0;
            for j in 0..at[i1].valence as usize {
                if j != m {
                    neigh[n] = at[i1].neighbor[j];
                    n += 1;
                }
            }
            if n != 2 || n_rank[neigh[0] as usize] != n_rank[neigh[1] as usize] {
                n1 += 1;
                continue;
            }
            if at[i1].n_ring_system == at[neigh[0] as usize].n_ring_system {
                n1 += 1;
                continue;
            }
            n_visited.clear();
            n_visited.resize(num_atoms, 0);
            n_visited[i1] = 1;
            if paths_have_identical_known_parities(
                at,
                i1 as u16,
                neigh[0],
                i1 as u16,
                neigh[1],
                &mut n_visited,
                n_rank,
                n_canon_rank,
                1,
            ) != 0
            {
                if !remove_one_stereo_bond(at, i1, n1 as usize) {
                    return CT_STEREOBOND_ERROR;
                }
                n1 -= 1;
                let n_atom_rank1 = n_canon_rank[i1].max(n_canon_rank[opposite_atom]);
                let n_atom_rank2 = n_canon_rank[i1].min(n_canon_rank[opposite_atom]);
                let mut mm = cs.n_len_dble - 1;
                let mut nn = 0;
                let mut found = false;
                while nn <= mm {
                    let d = cs.dble.as_ref().unwrap()[nn as usize];
                    if d.at_num1 == n_atom_rank1 && d.at_num2 == n_atom_rank2 {
                        if nn < mm {
                            cs.dble
                                .as_mut()
                                .unwrap()
                                .copy_within(nn as usize + 1..mm as usize + 1, nn as usize);
                        }
                        cs.n_len_dble -= 1;
                        found = true;
                        break;
                    }
                    nn += 1;
                }
                mm = if found { -1 } else { mm };
                if mm >= 0 {
                    return CT_STEREOCOUNT_ERR;
                }
                ret += 1;
            }
            n1 += 1;
        }
    }
    ret
}

// ---------------------------------------------------------------------------
// ichicans.c: stereo center parities known in advance
// ---------------------------------------------------------------------------

/// `SetKnownStereoCenterParities`: returns the number of centers set
pub fn set_known_stereo_center_parities(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_canon_rank: &[u16],
    n_rank: &[u16],
    n_atom_number: &[u16],
) -> i32 {
    use super::ichierr::CT_STEREOCOUNT_ERR;
    use super::ichisort::comp_at_rank;
    let mut num_set = 0;
    for i in 0..num_atoms {
        if at[i].parity == 0 || at[i].stereo_bond_neighbor[0] != 0 {
            continue;
        }
        if i32::from(at[i].stereo_atom_parity) != AB_PARITY_CALC
            || !parity_well_def(i32::from(at[i].parity))
        {
            continue;
        }
        let num_neigh = at[i].valence as usize;
        let mut n_neigh_rank = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
        let mut n_neigh_canon_rank = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
        for j in 0..num_neigh {
            n_neigh_rank[j] = n_rank[at[i].neighbor[j] as usize];
        }
        let n_atom_rank = n_rank[i];
        let trans_i;
        if num_neigh == 1 {
            at[i].stereo_atom_parity = at[i].parity;
            trans_i = 0;
        } else {
            trans_i = insertions_sort(&mut n_neigh_rank[..num_neigh], comp_at_rank);
            let mut j = 1;
            while j < num_neigh {
                if n_neigh_rank[j - 1] == n_neigh_rank[j] {
                    break;
                }
                j += 1;
            }
            if j < num_neigh {
                continue;
            }
        }
        let mut prev_trans: i32 = -1;
        let mut trans_k = 0;
        let i_max = i32::from(n_atom_rank) - 1;
        let mut j1 = 0;
        while j1 <= i_max {
            let k = n_atom_number[(i_max - j1) as usize] as usize;
            if n_atom_rank != n_rank[k] {
                break;
            }
            j1 += 1;
            if at[k].valence as usize != num_neigh {
                return CT_STEREOCOUNT_ERR;
            }
            let mut j = 0;
            for m in 0..num_neigh {
                for n in 0..num_neigh {
                    if n_rank[at[k].neighbor[n] as usize] == n_neigh_rank[m] {
                        n_neigh_canon_rank[m] = n_canon_rank[at[k].neighbor[n] as usize];
                        j += 1;
                        break;
                    }
                }
            }
            if j != num_neigh {
                return CT_STEREOCOUNT_ERR;
            }
            trans_k = insertions_sort(&mut n_neigh_canon_rank[..num_neigh], comp_at_rank) % 2;
            if prev_trans < 0 {
                prev_trans = trans_k;
            } else if trans_k != prev_trans {
                break;
            }
        }
        if trans_k == prev_trans {
            at[i].stereo_atom_parity =
                (2 - (i32::from(at[i].parity) + trans_i + prev_trans) % 2) as i8;
            num_set += 1;
        }
    }
    num_set
}

/// `RemoveKnownNonStereoCenterParities`: removes the parities of centers that
/// are known not to be stereogenic (also from the linear CT)
pub fn remove_known_non_stereo_center_parities(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_canon_rank: &[u16],
    n_rank: &[u16],
    cs: &mut CanonStatStereo,
) -> i32 {
    use super::ichierr::CT_STEREOCOUNT_ERR;
    let mut ret = 0;
    let mut n_visited: Vec<u16> = Vec::new();
    for i in 0..num_atoms {
        if at[i].parity == 0 || at[i].stereo_bond_neighbor[0] != 0 {
            continue;
        }
        let sap = i32::from(at[i].stereo_atom_parity);
        if parity_val(sap) != AB_PARITY_CALC && parity_well_def(sap) {
            continue;
        }
        let num_neigh = at[i].valence as usize;
        let mut n_neigh_rank = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
        let mut n_neigh_ord = [0u16; MAX_NUM_STEREO_ATOM_NEIGH];
        for j in 0..num_neigh {
            n_neigh_rank[j] = n_rank[at[i].neighbor[j] as usize];
            n_neigh_ord[j] = j as u16;
        }
        if num_neigh == 1 {
            continue;
        }
        // CompRanksOrd
        insertions_sort(&mut n_neigh_ord[..num_neigh], |a, b| {
            let r = i32::from(n_neigh_rank[*a as usize]) - i32::from(n_neigh_rank[*b as usize]);
            if r == 0 {
                i32::from(*a) - i32::from(*b)
            } else {
                r
            }
        });
        let mut k = 1usize;
        let mut j = 1usize;
        while k != 0 && j < num_neigh {
            if at[i].n_ring_system
                != at[at[i].neighbor[n_neigh_ord[j] as usize] as usize].n_ring_system
                && n_neigh_rank[n_neigh_ord[j - 1] as usize]
                    == n_neigh_rank[n_neigh_ord[j] as usize]
            {
                k = j;
                loop {
                    n_visited.clear();
                    n_visited.resize(num_atoms, 0);
                    n_visited[i] = 1;
                    if paths_have_identical_known_parities(
                        at,
                        i as u16,
                        at[i].neighbor[n_neigh_ord[j - 1] as usize],
                        i as u16,
                        at[i].neighbor[n_neigh_ord[k] as usize],
                        &mut n_visited,
                        n_rank,
                        n_canon_rank,
                        1,
                    ) != 0
                    {
                        at[i].parity = 0;
                        at[i].stereo_atom_parity = 0;
                        at[i].final_parity = 0;
                        let m = cs.n_len_carb - 1;
                        let mut n = 0;
                        while n <= m {
                            if cs.carb.as_ref().unwrap()[n as usize].at_num == n_canon_rank[i] {
                                if n < m {
                                    cs.carb
                                        .as_mut()
                                        .unwrap()
                                        .copy_within(n as usize + 1..m as usize + 1, n as usize);
                                }
                                cs.n_len_carb -= 1;
                                k = 0;
                                break;
                            }
                            n += 1;
                        }
                        if k != 0 {
                            return CT_STEREOCOUNT_ERR;
                        }
                        ret += 1;
                        break;
                    }
                    k += 1;
                    if !(k < num_neigh
                        && n_neigh_rank[n_neigh_ord[j - 1] as usize]
                            == n_neigh_rank[n_neigh_ord[k] as usize])
                    {
                        break;
                    }
                }
            }
            j += 1;
        }
    }
    ret
}

/// `MarkKnownEqualStereoCenterParities`
pub fn mark_known_equal_stereo_center_parities(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_rank: &[u16],
    n_atom_number: &[u16],
) -> i32 {
    let mut num_centers = 0;
    for i in 0..num_atoms {
        if at[i].parity == 0 || at[i].stereo_bond_neighbor[0] != 0 {
            continue;
        }
        if at[i].b_has_stereo_or_equ_to_stereo != 0 {
            continue;
        }
        if i32::from(at[i].stereo_atom_parity) & KNOWN_PARITIES_EQL != 0 {
            continue;
        }
        let parity = parity_val(i32::from(at[i].stereo_atom_parity));
        if parity == AB_PARITY_NONE {
            continue;
        }
        let n_atom_rank = n_rank[i];
        let mut b_different_parities: i32 = -1;
        let i_max = i32::from(n_atom_rank) - 1;
        let mut j1 = 0;
        while j1 <= i_max {
            let k = n_atom_number[(i_max - j1) as usize] as usize;
            if n_atom_rank != n_rank[k] {
                break;
            }
            j1 += 1;
            let parity_k = parity_val(i32::from(at[k].stereo_atom_parity));
            if parity_k != parity {
                b_different_parities = 1;
            } else if b_different_parities < 0 {
                b_different_parities = 0;
            }
            if parity_k == 0 {
                at[k].b_has_stereo_or_equ_to_stereo = 2;
            } else if at[k].b_has_stereo_or_equ_to_stereo == 0 {
                at[k].b_has_stereo_or_equ_to_stereo = 1;
            }
        }
        if b_different_parities == 0 && parity_known(parity) {
            let mut j1 = 0;
            while j1 <= i_max {
                let k = n_atom_number[(i_max - j1) as usize] as usize;
                if n_atom_rank != n_rank[k] {
                    break;
                }
                j1 += 1;
                at[k].stereo_atom_parity |= KNOWN_PARITIES_EQL as i8;
                num_centers += 1;
            }
        }
    }
    num_centers
}

/// `FillOutStereoParities`: makes the stereo descriptors fit the molecular
/// symmetry and removes parities from obviously non-stereo atoms and bonds
#[allow(clippy::too_many_arguments)]
pub fn fill_out_stereo_parities(
    at: &mut [SpAtom],
    num_atoms: usize,
    n_canon_rank: &[u16],
    n_atom_number_canon: &[u16],
    n_rank: &[u16],
    n_atom_number: &[u16],
    cs: &mut CanonStatStereo,
    b_isotopic: bool,
) -> i32 {
    use super::ichierr::{CT_ERR_MAX, CT_ERR_MIN};
    let mut ret = unmark_non_stereo(at, num_atoms, n_rank, n_atom_number, b_isotopic);
    if ret < 0 {
        return ret;
    }
    ret = fill_all_stereo_descriptors(at, num_atoms, n_canon_rank, n_atom_number_canon, cs);
    if ret == 0 {
        ret = cs.n_len_carb + cs.n_len_dble;
    }
    if ret < 0 {
        return ret;
    }
    let mut ret2 =
        set_known_stereo_center_parities(at, num_atoms, n_canon_rank, n_rank, n_atom_number);
    if ret2 >= 0 {
        ret2 = mark_known_equal_stereo_center_parities(at, num_atoms, n_rank, n_atom_number);
    }
    if ret2 >= 0 {
        ret2 = set_known_stereo_bond_parities(at, num_atoms, n_canon_rank, n_rank, n_atom_number);
        if ret2 >= 0 {
            ret2 = mark_known_equal_stereo_bond_parities(at, num_atoms, n_rank, n_atom_number);
        }
    }
    if ret2 >= 0 {
        loop {
            ret2 = remove_known_non_stereo_center_parities(at, num_atoms, n_canon_rank, n_rank, cs);
            if ret2 >= 0 {
                let ret3 =
                    remove_known_non_stereo_bond_parities(at, num_atoms, n_canon_rank, n_rank, cs);
                ret2 = if ret3 >= 0 { ret2 + ret3 } else { ret3 };
            }
            if ret2 <= 0 {
                break;
            }
        }
    }
    if (CT_ERR_MIN..=CT_ERR_MAX).contains(&ret2) {
        ret = ret2;
    }
    ret
}
