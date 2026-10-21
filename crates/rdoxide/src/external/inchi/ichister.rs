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
