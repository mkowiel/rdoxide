//! Tautomer endpoints and charge groups (`INCHI_BASE/src/ichitaut.c`): the
//! endpoint description, `GetChargeType`, `AddAtom2num`, `AddAtom2DA`,
//! `AddEndPoint` and the `nGetEndpointInfo*` family.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! Build facts: `KETO_ENOL_TAUT`, `TAUT_PT_22_00` .. `TAUT_PT_18_00` and
//! `ADD_MOVEABLE_O_PLUS` are all 1.

use super::incomdef::{NUM_H_ISOTOPES, RADICAL_SINGLET};
use super::inpdef::{ATT_ACIDIC_CO, InpAtom};
use super::util::{
    EL_NUMBER_AS, EL_NUMBER_BR, EL_NUMBER_C, EL_NUMBER_CL, EL_NUMBER_I, EL_NUMBER_N, EL_NUMBER_O,
    EL_NUMBER_P, EL_NUMBER_S, EL_NUMBER_SB, EL_NUMBER_SE, EL_NUMBER_TE, get_endpoint_valence,
    get_endpoint_valence_ket,
};

pub const T_NUM_NO_ISOTOPIC: usize = 2;
pub const T_NUM_ISOTOPIC: usize = NUM_H_ISOTOPES;

// `TG_NUM_DA` counters
pub const TG_NUM_DH: usize = 0;
pub const TG_NUM_DM: usize = 1;
pub const TG_NUM_AH: usize = 2;
pub const TG_NUM_AM: usize = 3;
pub const TG_NUM_DO: usize = 4;
pub const TG_NUM_AO: usize = 5;
pub const TG_NUM_DA: usize = 6;

// charge subtype bits
pub const C_SUBTYPE_CHARGED: i8 = 0;
pub const C_SUBTYPE_P_DONOR: i8 = 1;
pub const C_SUBTYPE_P_ACCEPT: i8 = 2;
pub const C_SUBTYPE_H_ACCEPT: i8 = 4;
pub const C_SUBTYPE_H_DONOR: i8 = 8;
pub const C_SUBTYPE_NEUTRAL: i8 = 16;
pub const C_SUBTYPE_CHARGED_NON_TAUT: i8 = C_SUBTYPE_CHARGED;
pub const C_SUBTYPE_CHARGED_P_DONOR: i8 = C_SUBTYPE_CHARGED | C_SUBTYPE_P_DONOR;
pub const C_SUBTYPE_CHARGED_H_ACCEPT: i8 = C_SUBTYPE_CHARGED | C_SUBTYPE_H_ACCEPT;
pub const C_SUBTYPE_CHARGED_H_ACCEPT_P_DONOR: i8 =
    C_SUBTYPE_CHARGED | C_SUBTYPE_H_ACCEPT | C_SUBTYPE_P_DONOR;
pub const C_SUBTYPE_CHARGED_H_DONOR: i8 = C_SUBTYPE_CHARGED | C_SUBTYPE_H_DONOR | C_SUBTYPE_P_DONOR;
pub const C_SUBTYPE_NEUTRAL_NON_TAUT: i8 = C_SUBTYPE_NEUTRAL;
pub const C_SUBTYPE_NEUTRAL_H_ACCEPT: i8 = C_SUBTYPE_NEUTRAL | C_SUBTYPE_H_ACCEPT;
pub const C_SUBTYPE_NEUTRAL_H_ACCEPT_P_ACCEPT: i8 =
    C_SUBTYPE_NEUTRAL | C_SUBTYPE_H_ACCEPT | C_SUBTYPE_P_ACCEPT;
pub const C_SUBTYPE_NEUTRAL_H_DONOR: i8 = C_SUBTYPE_NEUTRAL | C_SUBTYPE_H_DONOR;

/// `ENDPOINT_INFO`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EndpointInfo {
    pub c_moveable_charge: i8,
    pub c_neutral_bonds_valence: i8,
    pub c_mobile: i8,
    pub c_donor: i8,
    pub c_acceptor: i8,
    /// 1 => carbon, 2 => oxygen
    pub c_keto_enol_code: i8,
}

/// `T_ENDPOINT`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TEndpoint {
    pub num: [u16; T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC],
    pub num_da: [u16; TG_NUM_DA],
    pub n_group_number: u16,
    pub n_equ_number: u16,
    pub n_atom_number: u16,
}

/// `is_centerpoint_elem`
pub fn is_centerpoint_elem(el_number: u8) -> bool {
    matches!(
        el_number,
        EL_NUMBER_C
            | EL_NUMBER_N
            | EL_NUMBER_P
            | EL_NUMBER_S
            | EL_NUMBER_I
            | EL_NUMBER_AS
            | EL_NUMBER_SB
            | EL_NUMBER_SE
            | EL_NUMBER_TE
            | EL_NUMBER_CL
            | EL_NUMBER_BR
    )
}

/// `is_centerpoint_elem_KET`
pub fn is_centerpoint_elem_ket(el_number: u8) -> bool {
    el_number == EL_NUMBER_C
}

/// `is_centerpoint_elem_strict`
pub fn is_centerpoint_elem_strict(el_number: u8) -> bool {
    matches!(
        el_number,
        EL_NUMBER_C | EL_NUMBER_N | EL_NUMBER_P | EL_NUMBER_AS | EL_NUMBER_SB
    )
}

/// `AddAtom2num`: `b_subtract`: 0 = add, 1 = subtract, 2 = fill. Returns the
/// number of mobile groups (H and negative charge)
pub fn add_atom2num(num: &mut [u16], at: &InpAtom, b_subtract: i32) -> i32 {
    let mut n_mobile = i32::from(at.charge == -1);
    if b_subtract == 1 {
        num[1] = num[1].wrapping_sub(n_mobile as u16);
        n_mobile += i32::from(at.num_h);
        num[0] = num[0].wrapping_sub(n_mobile as u16);
        for k in 0..T_NUM_ISOTOPIC {
            let v = at.num_iso_h[NUM_H_ISOTOPES - k - 1] as u16;
            num[T_NUM_NO_ISOTOPIC + k] = num[T_NUM_NO_ISOTOPIC + k].wrapping_sub(v);
        }
    } else {
        if b_subtract == 2 {
            for v in num.iter_mut().take(T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC) {
                *v = 0;
            }
        }
        num[1] = num[1].wrapping_add(n_mobile as u16);
        n_mobile += i32::from(at.num_h);
        num[0] = num[0].wrapping_add(n_mobile as u16);
        for k in 0..T_NUM_ISOTOPIC {
            let v = at.num_iso_h[NUM_H_ISOTOPES - k - 1] as u16;
            num[T_NUM_NO_ISOTOPIC + k] = num[T_NUM_NO_ISOTOPIC + k].wrapping_add(v);
        }
    }
    n_mobile
}

/// `AddAtom2DA`: the donor/acceptor counters of an endpoint
pub fn add_atom2da(num_da: &mut [u16], at: &InpAtom, b_subtract: i32) {
    if at.charge < -1 || (at.charge == 1 && at.c_point == 0) || at.charge > 1 {
        return;
    }
    let n_delta: i32 = if b_subtract == 1 { -1 } else { 1 };
    let n_acidic_o = if at.at_type & ATT_ACIDIC_CO != 0 {
        n_delta
    } else {
        0
    };
    let add = |v: &mut u16, d: i32| *v = (i32::from(*v) + d) as u16;
    if b_subtract == 2 {
        for v in num_da.iter_mut().take(TG_NUM_DA) {
            *v = 0;
        }
    }
    let (val, cbv) = (i32::from(at.valence), i32::from(at.chem_bonds_valence));
    if (at.charge <= 0 && val == cbv) || (at.charge > 0 && val + 1 == cbv) {
        if at.charge < 0 {
            add(&mut num_da[TG_NUM_DM], n_delta);
            add(&mut num_da[TG_NUM_DO], n_acidic_o);
        } else if at.num_h != 0 {
            add(&mut num_da[TG_NUM_DH], n_delta);
            add(&mut num_da[TG_NUM_DO], n_acidic_o);
        }
    } else if (at.charge <= 0 && val + 1 == cbv) || (at.charge > 0 && val + 2 == cbv) {
        if at.charge < 0 {
            add(&mut num_da[TG_NUM_AM], n_delta);
        } else if at.num_h != 0 {
            add(&mut num_da[TG_NUM_AH], n_delta);
        } else {
            add(&mut num_da[TG_NUM_AO], n_acidic_o);
        }
    }
}

/// `AddEndPoint`
pub fn add_end_point(p_end_point: &mut TEndpoint, at: &[InpAtom], iat: usize) -> i32 {
    p_end_point.n_atom_number = iat as u16;
    p_end_point.n_equ_number = 0;
    p_end_point.n_group_number = at[iat].endpoint;
    if at[iat].endpoint != 0 {
        p_end_point.num = [0; T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC];
    } else {
        add_atom2num(&mut p_end_point.num, &at[iat], 2);
        add_atom2da(&mut p_end_point.num_da, &at[iat], 2);
    }
    0
}

/// `CHARGE_TYPE`
struct ChargeType {
    elname: &'static str,
    charge: i8,
    neutral_valence: i8,
    neutral_bonds_valence: i8,
    c_change_valence: i8,
    c_charge_type: i8,
    num_bonds: i8,
}

/// `CType[]` (`ADD_MOVEABLE_O_PLUS == 1`)
const C_TYPE: [ChargeType; 6] = [
    ChargeType {
        elname: "N",
        charge: 1,
        neutral_valence: 3,
        neutral_bonds_valence: 3,
        c_change_valence: 1,
        c_charge_type: 0,
        num_bonds: 0,
    },
    ChargeType {
        elname: "P",
        charge: 1,
        neutral_valence: 3,
        neutral_bonds_valence: 3,
        c_change_valence: 1,
        c_charge_type: 1,
        num_bonds: 0,
    },
    ChargeType {
        elname: "O",
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        c_change_valence: 1,
        c_charge_type: 2,
        num_bonds: 2,
    },
    ChargeType {
        elname: "S",
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        c_change_valence: 1,
        c_charge_type: 3,
        num_bonds: 2,
    },
    ChargeType {
        elname: "Se",
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        c_change_valence: 1,
        c_charge_type: 4,
        num_bonds: 2,
    },
    ChargeType {
        elname: "Te",
        charge: 1,
        neutral_valence: 2,
        neutral_bonds_valence: 2,
        c_change_valence: 1,
        c_charge_type: 5,
        num_bonds: 2,
    },
];

/// `bCanBeACPoint`: whether the atom can be a charge point of the given type
/// (see the comments of the C code); fills `c_charge_subtype`
pub fn can_be_ac_point(
    at: &InpAtom,
    c_charge: i8,
    c_change_valence: i8,
    neutral_bonds_valence: i8,
    neutral_valence: i8,
    n_endpoint_valence: i8,
    c_charge_subtype: &mut i8,
) -> bool {
    let b_neg_charge = i32::from(at.charge == -1);
    let (val, cbv, nh) = (
        i32::from(at.valence),
        i32::from(at.chem_bonds_valence),
        i32::from(at.num_h),
    );
    let (nbv0, nv0) = (i32::from(neutral_bonds_valence), i32::from(neutral_valence));
    if at.charge == c_charge && val == cbv && nh != 0 {
        let n_change_valence = i32::from(at.charge) * i32::from(c_change_valence);
        let n_bonds_valence = cbv + nh;
        if n_bonds_valence == nbv0 + n_change_valence && n_endpoint_valence != 0 {
            *c_charge_subtype = C_SUBTYPE_CHARGED_P_DONOR;
        }
        return false;
    }
    if at.charge == c_charge && val < cbv {
        let n_change_valence = i32::from(at.charge) * i32::from(c_change_valence);
        let n_bonds_valence = cbv + nh;
        let n_num_bonds = val + nh;
        if n_bonds_valence == nbv0 + n_change_valence {
            if n_num_bonds == nv0 {
                if val == nv0 || n_endpoint_valence == 0 {
                    *c_charge_subtype = C_SUBTYPE_CHARGED_NON_TAUT;
                } else {
                    *c_charge_subtype = C_SUBTYPE_CHARGED_H_DONOR;
                }
                return true;
            }
            if n_num_bonds == nv0 - 1 {
                if n_endpoint_valence != 0 {
                    *c_charge_subtype = if nh != 0 {
                        C_SUBTYPE_CHARGED_H_ACCEPT_P_DONOR
                    } else {
                        C_SUBTYPE_CHARGED_H_ACCEPT
                    };
                } else {
                    *c_charge_subtype = C_SUBTYPE_CHARGED_NON_TAUT;
                }
                return true;
            }
        }
    } else if at.charge == 0 || b_neg_charge != 0 {
        let n_bonds_valence = cbv + nh + b_neg_charge;
        let n_num_bonds = val + nh + b_neg_charge;
        if n_bonds_valence == nbv0 {
            if n_num_bonds == nv0 {
                if val == nv0 || n_endpoint_valence == 0 {
                    *c_charge_subtype = C_SUBTYPE_NEUTRAL_NON_TAUT;
                } else if val < nv0 {
                    *c_charge_subtype = if nh != 0 {
                        C_SUBTYPE_NEUTRAL_H_DONOR
                    } else {
                        C_SUBTYPE_NEUTRAL_H_ACCEPT
                    };
                } else {
                    return false;
                }
                return true;
            }
            if n_num_bonds == nv0 - 1 {
                if n_endpoint_valence != 0 {
                    *c_charge_subtype = C_SUBTYPE_NEUTRAL_H_ACCEPT_P_ACCEPT;
                } else {
                    *c_charge_subtype = C_SUBTYPE_NEUTRAL_NON_TAUT;
                }
                return true;
            }
        }
    }
    false
}

/// `GetChargeType`: returns the charge type (>= 0) or -1; `c_charge_subtype`
/// receives the subtype
pub fn get_charge_type(atom: &[InpAtom], iat: usize, c_charge_subtype: &mut i8) -> i32 {
    let at = &atom[iat];
    *c_charge_subtype = 0;
    if at.charge.abs() == 1 {
        for i in 0..at.valence as usize {
            let n = at.neighbor[i] as usize;
            if (i32::from(atom[n].charge) + i32::from(at.charge)).abs()
                < (i32::from(atom[n].charge) - i32::from(at.charge)).abs()
                && atom[n].endpoint == 0
            {
                return -1;
            }
        }
    } else if at.charge != 0 {
        return -1;
    }
    let end = at
        .elname
        .iter()
        .position(|&c| c == 0)
        .unwrap_or(at.elname.len());
    let name = &at.elname[..end];
    for ct in &C_TYPE {
        if name == ct.elname.as_bytes()
            && (ct.num_bonds == 0
                || (ct.num_bonds == at.valence && at.n_num_at_in_ring_system >= 5))
        {
            let n_endpoint_valence = get_endpoint_valence(at.el_number) as i8;
            if can_be_ac_point(
                at,
                ct.charge,
                ct.c_change_valence,
                ct.neutral_bonds_valence,
                ct.neutral_valence,
                n_endpoint_valence,
                c_charge_subtype,
            ) {
                return i32::from(ct.c_charge_type);
            }
        }
    }
    -1
}

/// The common positive charge-point branch of the `nGetEndpointInfo*`
/// functions; `strict` = the branch returns 0 when the subtype is neither an
/// H acceptor nor an H donor (always the case but in `PT_18_00`)
fn charge_point_info(
    atom: &[InpAtom],
    iat: usize,
    n_endpoint_valence: i32,
    eif: &mut EndpointInfo,
    keto: i8,
) -> i32 {
    let a = &atom[iat];
    let mut sub = 0i8;
    if a.c_point != 0
        && get_charge_type(atom, iat, &mut sub) >= 0
        && sub & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR) != 0
    {
        if sub & C_SUBTYPE_H_ACCEPT != 0 {
            eif.c_donor = 0;
            eif.c_acceptor = 1;
        } else if sub & C_SUBTYPE_H_DONOR != 0 {
            eif.c_donor = 1;
            eif.c_acceptor = 0;
        } else {
            return 0;
        }
        eif.c_mobile = a.num_h;
        eif.c_neutral_bonds_valence = (n_endpoint_valence - i32::from(a.num_h)) as i8;
        eif.c_moveable_charge = a.charge;
        eif.c_keto_enol_code = keto;
        return n_endpoint_valence;
    }
    0
}

/// The neutral/anion branch shared by the variants (`donor_by_mobile`: the
/// donor/acceptor is decided by the number of mobile groups, as in `PT_13`
/// and `PT_18`; otherwise by `chem_bonds_valence - valence`)
fn neutral_info(a: &InpAtom, n_endpoint_valence: i32, eif: &mut EndpointInfo, keto: i8) -> i32 {
    if n_endpoint_valence < i32::from(a.chem_bonds_valence) {
        return 0;
    }
    let n_mobile = i32::from(a.num_h) + i32::from(a.charge == -1);
    if n_mobile + i32::from(a.chem_bonds_valence) != n_endpoint_valence {
        return 0;
    }
    match i32::from(a.chem_bonds_valence) - i32::from(a.valence) {
        0 => {
            eif.c_donor = 1;
            eif.c_acceptor = 0;
        }
        1 => {
            eif.c_donor = 0;
            eif.c_acceptor = 1;
        }
        _ => return 0,
    }
    eif.c_mobile = n_mobile as i8;
    eif.c_neutral_bonds_valence = (n_endpoint_valence - n_mobile) as i8;
    eif.c_moveable_charge = 0;
    eif.c_keto_enol_code = keto;
    n_endpoint_valence
}

fn is_radical(a: &InpAtom) -> bool {
    a.radical != 0 && i32::from(a.radical) != RADICAL_SINGLET
}

/// Shared skeleton of `nGetEndpointInfo`, `_PT_22_00`, `_PT_16_00`,
/// `_PT_06_00`, `_PT_39_00` (`valence` checks as in the C functions)
fn endpoint_info_generic(
    atom: &[InpAtom],
    iat: usize,
    eif: &mut EndpointInfo,
    n_endpoint_valence: i32,
    extra_reject: impl Fn(&InpAtom, i32) -> bool,
) -> i32 {
    let a = &atom[iat];
    if is_radical(a) {
        return 0;
    }
    if n_endpoint_valence == 0 {
        return 0;
    }
    if n_endpoint_valence <= i32::from(a.valence) {
        return 0;
    }
    if extra_reject(a, n_endpoint_valence) {
        return 0;
    }
    if a.charge == -1 || a.charge == 0 {
        neutral_info(a, n_endpoint_valence, eif, 0)
    } else {
        charge_point_info(atom, iat, n_endpoint_valence, eif, 0)
    }
}

/// `nGetEndpointInfo`
pub fn n_get_endpoint_info(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let v = get_endpoint_valence(atom[iat].el_number);
    endpoint_info_generic(atom, iat, eif, v, |_, _| false)
}

/// `nGetEndpointInfo_PT_22_00`: only carbon
pub fn n_get_endpoint_info_pt_22_00(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let v = if atom[iat].el_number == EL_NUMBER_C {
        4
    } else {
        0
    };
    endpoint_info_generic(atom, iat, eif, v, |_, _| false)
}

/// `nGetEndpointInfo_PT_16_00`
pub fn n_get_endpoint_info_pt_16_00(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let v = get_endpoint_valence_ket(atom[iat].el_number);
    endpoint_info_generic(atom, iat, eif, v, |a, v| {
        (v == 4 && a.valence < 2) || (v == 2 && a.valence > 1)
    })
}

/// `nGetEndpointInfo_PT_06_00`: carbon and the usual endpoints
pub fn n_get_endpoint_info_pt_06_00(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let el = atom[iat].el_number;
    let v = if el == EL_NUMBER_C {
        4
    } else {
        get_endpoint_valence(el)
    };
    endpoint_info_generic(atom, iat, eif, v, |_, _| false)
}

/// `nGetEndpointInfo_PT_39_00`: C and N
pub fn n_get_endpoint_info_pt_39_00(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let el = atom[iat].el_number;
    let v = if el == EL_NUMBER_C {
        4
    } else if el == EL_NUMBER_N {
        3
    } else {
        0
    };
    endpoint_info_generic(atom, iat, eif, v, |_, _| false)
}

/// `nGetEndpointInfo_PT_13_00`: C, S, O, Se, Te; the donor/acceptor is
/// decided by the number of mobile groups
pub fn n_get_endpoint_info_pt_13_00(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let a = &atom[iat];
    if is_radical(a) {
        return 0;
    }
    let v = match a.el_number {
        EL_NUMBER_C => 4,
        EL_NUMBER_S | EL_NUMBER_O | EL_NUMBER_SE | EL_NUMBER_TE => 2,
        _ => 0,
    };
    if v == 0 || v <= i32::from(a.valence) {
        return 0;
    }
    if a.charge == -1 || a.charge == 0 {
        if v < i32::from(a.chem_bonds_valence) {
            return 0;
        }
        let n_mobile = i32::from(a.num_h) + i32::from(a.charge == -1);
        if n_mobile + i32::from(a.chem_bonds_valence) != v {
            return 0;
        }
        if n_mobile > 0 {
            eif.c_donor = 1;
            eif.c_acceptor = 0;
        } else {
            eif.c_donor = 0;
            eif.c_acceptor = 1;
        }
        eif.c_mobile = n_mobile as i8;
        eif.c_neutral_bonds_valence = (v - n_mobile) as i8;
        eif.c_moveable_charge = 0;
        eif.c_keto_enol_code = 0;
        return v;
    }
    charge_point_info(atom, iat, v, eif, 0)
}

/// `nGetEndpointInfo_PT_18_00`: O and N, no valence checks
pub fn n_get_endpoint_info_pt_18_00(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let a = &atom[iat];
    if is_radical(a) {
        return 0;
    }
    let v = match a.el_number {
        EL_NUMBER_O => 2,
        EL_NUMBER_N => 3,
        _ => 0,
    };
    if v == 0 {
        return 0;
    }
    if a.charge == -1 || a.charge == 0 {
        let n_mobile = i32::from(a.num_h);
        if n_mobile > 0 {
            eif.c_donor = 1;
            eif.c_acceptor = 0;
        } else {
            eif.c_donor = 0;
            eif.c_acceptor = 1;
        }
        eif.c_mobile = n_mobile as i8;
        eif.c_neutral_bonds_valence = (v - n_mobile) as i8;
        eif.c_moveable_charge = 0;
        eif.c_keto_enol_code = 0;
        return v;
    }
    let mut sub = 0i8;
    if a.c_point != 0
        && get_charge_type(atom, iat, &mut sub) >= 0
        && sub & (C_SUBTYPE_H_ACCEPT | C_SUBTYPE_H_DONOR) != 0
    {
        if sub & C_SUBTYPE_H_ACCEPT != 0 {
            eif.c_donor = 0;
            eif.c_acceptor = 1;
        } else if sub & C_SUBTYPE_H_DONOR != 0 {
            eif.c_donor = 1;
            eif.c_acceptor = 0;
        }
    } else if a.num_h != 0 {
        eif.c_donor = 1;
        eif.c_acceptor = 0;
    } else {
        eif.c_donor = 0;
        eif.c_acceptor = 1;
    }
    eif.c_mobile = a.num_h;
    eif.c_neutral_bonds_valence = (v - i32::from(a.num_h)) as i8;
    eif.c_moveable_charge = a.charge;
    eif.c_keto_enol_code = 0;
    v
}

/// `nGetEndpointInfo_KET`: O and C for keto-enol tautomerism
pub fn n_get_endpoint_info_ket(atom: &[InpAtom], iat: usize, eif: &mut EndpointInfo) -> i32 {
    let v = get_endpoint_valence_ket(atom[iat].el_number);
    let code: i8 = match v {
        2 => 1,
        4 => 2,
        _ => 0,
    };
    let a = &atom[iat];
    if is_radical(a) || v == 0 || v <= i32::from(a.valence) {
        return 0;
    }
    if (v == 4 && a.valence < 2) || (v == 2 && a.valence > 1) {
        return 0;
    }
    if a.charge == -1 || a.charge == 0 {
        neutral_info(a, v, eif, code)
    } else {
        charge_point_info(atom, iat, v, eif, code)
    }
}

// ---- charge groups and salt charge types ----

use super::ichierr::BNS_CPOINT_ERR;
use super::incomdef::{BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE};
use super::inpdef::BOND_TYPE_MASK;
use super::util::get_el_valence;

pub const SALT_DONOR_H: i32 = 1;
pub const SALT_DONOR_NEG: i32 = 2;
pub const SALT_ACCEPTOR: i32 = 4;
pub const SALT_P_DONOR: i32 = 8;
pub const SALT_P_ACCEPTOR: i32 = 16;

pub const BOND_ALTERN: i32 = 4;
pub const BOND_TAUTOM: i32 = 8;
pub const BOND_ALT12NS: i32 = 9;

/// `T_GROUP` (the members the charge/salt code uses)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TGroup {
    pub num: [u16; T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC],
    pub num_da: [u16; TG_NUM_DA],
    pub n_group_number: u16,
    pub n_num_endpoints: u16,
    pub n_first_endpoint_at_no_pos: u16,
}

/// `TNI`: the protons and explicit H removed by the normalization
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Tni {
    pub n_num_removed_explicit_h: i32,
    pub n_num_removed_protons: i32,
    pub n_num_removed_protons_isotopic: [i32; NUM_H_ISOTOPES],
    pub b_normalization_flags: u64,
}

/// `T_GROUP_INFO` (the t-group list)
#[derive(Debug, Clone, Default)]
pub struct TGroupInfo {
    pub t_group: Option<Vec<TGroup>>,
    pub num_t_groups: usize,
    pub max_num_t_groups: usize,
    /// `tGroupNumber`: index + 1 of the t-group with the given number
    pub t_group_number: Vec<u16>,
    pub tni: Tni,
    pub b_taut_flags: u64,
}

/// `C_GROUP_INFO` (the c-group list)
#[derive(Debug, Clone, Default)]
pub struct CGroupInfo {
    pub c_group: Option<Vec<CGroup>>,
    pub num_c_groups: usize,
    pub max_num_c_groups: usize,
    pub c_candidate: Vec<CCandidate>,
    pub max_num_candidates: i32,
    /// 0 = uninitialized, -1 = no candidates found
    pub num_candidates: i32,
}

/// `C_CANDIDATE`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CCandidate {
    pub atnumber: u16,
    pub type_: i8,
    pub subtype: i8,
}

/// `C_GROUP`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CGroup {
    pub num: [u16; 2],
    pub num_c_points: u16,
    pub n_group_number: u16,
    pub c_group_type: u8,
}

/// `S_CANDIDATE`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SCandidate {
    pub atnumber: u16,
    pub type_: i8,
    pub subtype: i8,
    pub endpoint: u16,
}

/// `CmpCCandidates`
pub fn cmp_c_candidates(c1: &CCandidate, c2: &CCandidate) -> i32 {
    let ret = i32::from(c1.type_) - i32::from(c2.type_);
    if ret != 0 {
        return ret;
    }
    let ret = i32::from(c1.subtype) - i32::from(c2.subtype);
    if ret != 0 {
        return ret;
    }
    i32::from(c1.atnumber) - i32::from(c2.atnumber)
}

/// `RegisterCPoints`: registers the charge points `point1` and `point2` in
/// the charge groups (new group, add to a group, or merge two groups);
/// `c_group` holds `max_num_c` elements, `pnum_c` the number in use.
pub fn register_c_points(
    c_group: &mut [CGroup],
    pnum_c: &mut usize,
    point1: usize,
    point2: usize,
    ctype: u8,
    at: &mut [InpAtom],
) -> i32 {
    let max_num_c = c_group.len();
    let mut num_c = *pnum_c;
    let (mut point1, mut point2) = (point1, point2);
    let charged = |a: &InpAtom| u16::from(a.charge == 1);

    if at[point1].c_point == at[point2].c_point {
        if at[point1].c_point != 0 {
            return 0;
        }
        if num_c < max_num_c {
            c_group[num_c] = CGroup::default();
            c_group[num_c].num[0] = charged(&at[point1]) + charged(&at[point2]);
            c_group[num_c].num_c_points += 2;
            c_group[num_c].c_group_type = ctype;
            let mut n_group_number = 0u16;
            for g in &c_group[..num_c] {
                if n_group_number < g.n_group_number {
                    n_group_number = g.n_group_number;
                }
            }
            n_group_number += 1;
            c_group[num_c].n_group_number = n_group_number;
            at[point1].c_point = n_group_number;
            at[point2].c_point = n_group_number;
            *pnum_c = num_c + 1;
            if at[point1].num_h != 0 || at[point2].num_h != 0 {
                c_group[num_c].num[1] += 1;
            }
            return 1;
        }
        return BNS_CPOINT_ERR;
    }
    if at[point1].c_point > at[point2].c_point {
        std::mem::swap(&mut point1, &mut point2);
    }
    if at[point1].c_point == 0 {
        let n_group_number = at[point2].c_point;
        for g in c_group[..num_c].iter_mut() {
            if n_group_number == g.n_group_number {
                at[point1].c_point = at[point2].c_point;
                g.num_c_points += 1;
                g.num[0] += u16::from(at[point1].charge == 1);
                return 1;
            }
        }
        return BNS_CPOINT_ERR;
    }
    // merge two c-groups
    let n_new_group_number = at[point1].c_point;
    let n_group_number = at[point2].c_point;
    let (mut i1, mut i2): (Option<usize>, Option<usize>) = (None, None);
    let mut i = 0;
    while i < num_c && (i1.is_none() || i2.is_none()) {
        if n_new_group_number == c_group[i].n_group_number {
            i1 = Some(i);
        } else if n_group_number == c_group[i].n_group_number {
            i2 = Some(i);
        }
        i += 1;
    }
    let (Some(i1), Some(i2)) = (i1, i2) else {
        return BNS_CPOINT_ERR;
    };
    let (n0, np) = (c_group[i2].num[0], c_group[i2].num_c_points);
    c_group[i1].num[0] += n0;
    c_group[i1].num_c_points += np;
    num_c -= 1;
    if num_c > i2 {
        c_group.copy_within(i2 + 1..num_c + 1, i2);
    }
    *pnum_c = num_c;
    for g in c_group[..num_c].iter_mut() {
        if g.n_group_number > n_group_number {
            g.n_group_number -= 1;
        }
    }
    for a in at.iter_mut() {
        if a.c_point > n_group_number {
            a.c_point -= 1;
        } else if a.c_point == n_group_number {
            a.c_point = n_new_group_number;
        }
    }
    1
}

/// `GetSaltChargeType`: the "salt" role of a terminal O, S, Se, Te attached to
/// a trigonal carbon: returns the type (0) or -1; `s_subtype` receives
/// `SALT_*` bits
pub fn get_salt_charge_type(
    at: &[InpAtom],
    at_no: usize,
    t_group_info: Option<&TGroupInfo>,
    s_subtype: &mut i32,
) -> i32 {
    *s_subtype = 0;
    let a = &at[at_no];
    if a.valence != 1
        || (a.radical != 0 && i32::from(a.radical) != RADICAL_SINGLET)
        || a.charge < -1
        || (a.charge > 0 && a.c_point == 0)
    {
        return -1;
    }
    let el = a.el_number;
    if !(el == EL_NUMBER_O || el == EL_NUMBER_S || el == EL_NUMBER_SE || el == EL_NUMBER_TE)
        || i32::from(a.chem_bonds_valence) + i32::from(a.num_h)
            != get_el_valence(i32::from(el), i32::from(a.charge), 0)
    {
        return -1;
    }
    let ic = a.neighbor[0] as usize;
    let c = &at[ic];
    // SALT_WITH_PROTONS == 1
    if c.el_number != EL_NUMBER_C
        || i32::from(c.chem_bonds_valence) + i32::from(c.num_h) != 4
        || c.charge != 0
        || (c.radical != 0 && i32::from(c.radical) != RADICAL_SINGLET)
        || c.valence == c.chem_bonds_valence
    {
        return -1;
    }
    let tg = a.endpoint;
    if tg != 0
        && let Some(info) = t_group_info
        && let Some(groups) = &info.t_group
    {
        for g in groups.iter().take(info.num_t_groups) {
            if tg == g.n_group_number {
                if g.num[0] > g.num[1] {
                    *s_subtype |= SALT_DONOR_H;
                }
                if g.num[1] != 0 {
                    *s_subtype |= SALT_DONOR_NEG;
                }
                *s_subtype |= SALT_ACCEPTOR;
                return 0;
            }
        }
        return -1;
    }
    if a.charge == -1 {
        *s_subtype |= SALT_DONOR_NEG;
    }
    if a.charge <= 0 && a.num_h != 0 {
        *s_subtype |= SALT_DONOR_H;
    }
    if a.charge == 0 && a.chem_bonds_valence == 2 {
        *s_subtype |= SALT_ACCEPTOR;
    }
    if a.charge == 1 && a.c_point != 0 && a.chem_bonds_valence == 2 && a.num_h != 0 {
        *s_subtype |= SALT_DONOR_H;
    }
    0
}

/// `bDoNotMergeNonTautAtom`
pub fn do_not_merge_non_taut_atom(at: &[InpAtom], at_no: usize) -> bool {
    at[at_no].el_number == EL_NUMBER_N
}

/// `GetOtherSaltChargeType`: type 1 for a tautomeric endpoint next to a
/// possible centerpoint, -1 to ignore
pub fn get_other_salt_charge_type(
    at: &[InpAtom],
    at_no: usize,
    t_group_info: Option<&TGroupInfo>,
    s_subtype: &mut i32,
    b_accept_o: bool,
) -> i32 {
    *s_subtype = 0;
    let a = &at[at_no];
    if !b_accept_o
        && (a.el_number == EL_NUMBER_O
            || a.el_number == EL_NUMBER_S
            || a.el_number == EL_NUMBER_SE
            || a.el_number == EL_NUMBER_TE)
    {
        return -1;
    }
    let mut eif = EndpointInfo::default();
    if n_get_endpoint_info(at, at_no, &mut eif) == 0 {
        return -1;
    }
    let mut num_centerpoints = 0;
    for j in 0..a.valence as usize {
        let bond_type = i32::from(a.bond_type[j]) & BOND_TYPE_MASK;
        let cp = a.neighbor[j] as usize;
        let c = &at[cp];
        let alt = bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM;
        if ((eif.c_acceptor != 0 && (bond_type == BOND_TYPE_DOUBLE || alt))
            || (eif.c_donor != 0 && (bond_type == BOND_TYPE_SINGLE || alt)))
            && (c.chem_bonds_valence > c.valence
                || (c.chem_bonds_valence == c.valence && (c.endpoint != 0 || c.c_point != 0)))
            && is_centerpoint_elem(c.el_number)
        {
            num_centerpoints += 1;
            break;
        }
    }
    if num_centerpoints == 0 {
        return -1;
    }
    let tg = a.endpoint;
    if tg != 0
        && let Some(info) = t_group_info
        && let Some(groups) = &info.t_group
    {
        for g in groups.iter().take(info.num_t_groups) {
            if tg == g.n_group_number {
                if g.num[0] > g.num[1] {
                    *s_subtype |= SALT_DONOR_H;
                }
                if g.num[1] != 0 {
                    *s_subtype |= SALT_DONOR_NEG;
                }
                *s_subtype |= SALT_ACCEPTOR;
                return 1;
            }
        }
        return -1;
    }
    if eif.c_acceptor != 0 {
        *s_subtype |= SALT_ACCEPTOR;
    }
    if eif.c_donor != 0 {
        if a.charge == -1 {
            *s_subtype |= SALT_DONOR_NEG;
        }
        if a.num_h != 0 {
            *s_subtype |= SALT_DONOR_H;
        }
    }
    1
}

/// `GetOtherSaltType`: type 2 for a terminal C-SH or C-S(-)
pub fn get_other_salt_type(at: &[InpAtom], at_no: usize, s_subtype: &mut i32) -> i32 {
    let a = &at[at_no];
    if a.valence != 1
        || a.chem_bonds_valence != 1
        || 1 != i32::from(a.num_h == 1) + i32::from(a.charge == -1)
    {
        return -1;
    }
    *s_subtype = 0;
    if !(a.el_number == EL_NUMBER_S || a.el_number == EL_NUMBER_SE || a.el_number == EL_NUMBER_TE) {
        return -1;
    }
    let mut eif = EndpointInfo::default();
    if n_get_endpoint_info(at, at_no, &mut eif) == 0
        || (eif.c_moveable_charge != 0 && a.c_point == 0)
        || eif.c_donor == 0
        || eif.c_acceptor != 0
    {
        return -1;
    }
    let cp = a.neighbor[0] as usize;
    let c = &at[cp];
    if c.el_number != EL_NUMBER_C
        || c.charge != 0
        || (c.radical != 0 && i32::from(c.radical) != RADICAL_SINGLET)
        || c.valence != c.chem_bonds_valence
    {
        return -1;
    }
    if a.num_h == 1 {
        *s_subtype |= SALT_P_DONOR;
    } else if a.charge == -1 {
        *s_subtype |= SALT_P_ACCEPTOR;
    } else {
        return -1;
    }
    2
}

/// `comp_candidates`
pub fn comp_candidates(s1: &SCandidate, s2: &SCandidate) -> i32 {
    if s1.type_ >= 0 && s2.type_ < 0 {
        return -1;
    }
    if s1.type_ < 0 && s2.type_ >= 0 {
        return 1;
    }
    if s1.endpoint != 0 && s2.endpoint == 0 {
        return -1;
    }
    if s1.endpoint == 0 && s2.endpoint != 0 {
        return 1;
    }
    if s1.endpoint != 0 && s2.endpoint != 0 {
        let ret = i32::from(s1.endpoint) - i32::from(s2.endpoint);
        if ret != 0 {
            return ret;
        }
    }
    i32::from(s1.atnumber) - i32::from(s2.atnumber)
}

/// `MarkChargeGroups`: finds the alternating paths between charged and
/// neutral charge-group candidates (N, P) and registers the charge groups;
/// returns the number of changes or an error
pub fn mark_charge_groups(
    at: &mut [InpAtom],
    num_atoms: i32,
    cgi: &mut CGroupInfo,
    p: &mut super::ichi_bns::BnStruct,
    bd: &mut super::ichi_bns::BnData,
) -> i32 {
    use super::ichi_bns::{ALT_PATH_MODE_CHARGE, exists_alt_path};
    use super::ichierr::{BNS_VERT_EDGE_OVFL, is_bns_error};
    let mut n_changes = 0;
    if cgi.c_candidate.is_empty() || cgi.max_num_candidates <= 0 {
        return n_changes;
    }
    let max_candidates = cgi.max_num_candidates as usize;
    let mut n_cand = cgi.num_candidates;
    if n_cand == -1 {
        n_cand = 0;
    }
    if n_cand == 0 {
        for i in 0..num_atoms as usize {
            let mut c_subtype = 0i8;
            let c_type = get_charge_type(at, i, &mut c_subtype);
            if c_type >= 0 {
                if n_cand as usize >= max_candidates {
                    return BNS_VERT_EDGE_OVFL;
                }
                cgi.c_candidate[n_cand as usize] = CCandidate {
                    atnumber: i as u16,
                    type_: c_type as i8,
                    subtype: c_subtype,
                };
                n_cand += 1;
            }
        }
        if n_cand <= 1 {
            cgi.num_candidates = -1;
            return 0;
        }
    }
    let n_cand = n_cand as usize;
    cgi.c_candidate[..n_cand].sort_by(|a, b| cmp_c_candidates(a, b).cmp(&0));
    let cand = cgi.c_candidate.clone();
    let mut num_tested = 0;
    let mut i1 = 0usize;
    'outer: while i1 < n_cand {
        while i1 < n_cand && cand[i1].subtype & C_SUBTYPE_NEUTRAL != 0 {
            i1 += 1;
        }
        if i1 == n_cand {
            break;
        }
        let mut i2 = i1 + 1;
        while i2 < n_cand
            && cand[i2].type_ == cand[i1].type_
            && cand[i2].subtype & C_SUBTYPE_NEUTRAL == 0
        {
            i2 += 1;
        }
        if i2 == n_cand {
            break;
        }
        let mut i3 = i2;
        while i3 < n_cand && cand[i3].type_ == cand[i1].type_ {
            i3 += 1;
        }
        if i3 == i2 {
            if i2 < n_cand {
                i1 = i3;
                continue;
            }
            break;
        }
        for i in i1..i2 {
            let iat1 = cand[i].atnumber as usize;
            for j in i2..i3 {
                num_tested += 1;
                let iat2 = cand[j].atnumber as usize;
                if at[iat1].c_point != 0 && at[iat1].c_point == at[iat2].c_point {
                    continue;
                }
                let ret = exists_alt_path(
                    p,
                    bd,
                    None,
                    at,
                    num_atoms,
                    iat1 as i32,
                    iat2 as i32,
                    ALT_PATH_MODE_CHARGE,
                );
                if is_bns_error(ret) {
                    return ret;
                }
                if ret & 1 != 0 {
                    let n_delta = (ret & !3) >> 2;
                    n_changes += ret & 2;
                    let max_cg = cgi.max_num_c_groups;
                    let cg = cgi.c_group.get_or_insert_with(Vec::new);
                    if cg.len() < max_cg {
                        cg.resize(max_cg, CGroup::default());
                    }
                    let r = register_c_points(
                        &mut cg[..max_cg],
                        &mut cgi.num_c_groups,
                        iat1,
                        iat2,
                        cand[i1].type_ as u8,
                        at,
                    );
                    if is_bns_error(r) {
                        return r;
                    }
                    if n_delta != 0 {
                        break 'outer;
                    }
                }
            }
        }
        i1 = i3;
    }
    if cgi.num_candidates == 0 {
        cgi.num_candidates = if num_tested != 0 { n_cand as i32 } else { -1 };
    }
    n_changes
}

/// `T_BONDPOS`
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TBondPos {
    pub n_atom_number: u16,
    pub neighbor_index: u16,
}

/// `BOND_MARK_ALL` (`extr_ct.h`)
pub const BOND_MARK_ALL: i32 = 0xf0;

/// `SetTautomericBonds`: changes the non-tautomeric bonds (single and double
/// ones) at the given positions to tautomeric ones; returns their number
pub fn set_tautomeric_bonds(at: &mut [InpAtom], bond_pos: &[TBondPos]) -> i32 {
    let mut n = 0;
    for bp in bond_pos {
        let neighbor_index = bp.neighbor_index as usize;
        let center = bp.n_atom_number as usize;
        let bond_mark = i32::from(at[center].bond_type[neighbor_index]);
        let bond_type = bond_mark & !BOND_MARK_ALL;
        // REPLACE_ALT_WITH_TAUT is 1
        if bond_type != BOND_TAUTOM {
            let new_type = ((bond_mark & BOND_MARK_ALL) | BOND_TAUTOM) as u8;
            at[center].bond_type[neighbor_index] = new_type;
            let neighbor = at[center].neighbor[neighbor_index] as usize;
            for ii in 0..at[neighbor].valence as usize {
                if at[neighbor].neighbor[ii] as usize == center {
                    at[neighbor].bond_type[ii] = new_type;
                    break;
                }
            }
            n += 1;
        }
    }
    n
}

/// `RegisterEndPoints`: puts the endpoints into the t-groups (creating,
/// extending and merging them) and re-builds the network groups. Returns
/// the number of changes (> 0), 0 for none or -1 for a program error.
pub fn register_end_points(
    tgi: &mut TGroupInfo,
    end_point: &mut [TEndpoint],
    at: &mut [InpAtom],
    num_atoms: i32,
    cgi: Option<&mut CGroupInfo>,
    p: Option<&mut super::ichi_bns::BnStruct>,
) -> i32 {
    let Some(mut tg) = tgi.t_group.take() else {
        return -1;
    };
    let max_num_t = tgi.max_num_t_groups;
    let mut num_t = tgi.num_t_groups;
    let ret = register_end_points_inner(
        tgi, &mut tg, &mut num_t, max_num_t, end_point, at, num_atoms, cgi, p,
    );
    tgi.t_group = Some(tg);
    ret
}

#[allow(clippy::too_many_arguments)]
fn register_end_points_inner(
    tgi: &mut TGroupInfo,
    t_group: &mut [TGroup],
    pnum_t: &mut usize,
    max_num_t: usize,
    end_point: &mut [TEndpoint],
    at: &mut [InpAtom],
    num_atoms: i32,
    cgi: Option<&mut CGroupInfo>,
    p: Option<&mut super::ichi_bns::BnStruct>,
) -> i32 {
    use super::ichi_bns::{
        add_c_groups_to_bn_struct, add_t_groups_to_bn_struct, re_init_bn_struct,
    };
    use super::ichierr::is_bns_error;
    let n_end = end_point.len();
    if n_end == 0 {
        return 0;
    }
    let mut num_t = *pnum_t;
    let mut difference = 0i32;
    let mut next_group: u16 = 0;
    for g in &t_group[..num_t] {
        if next_group < g.n_group_number {
            next_group = g.n_group_number;
        }
    }
    next_group += 1;
    let mut least_group = next_group;
    let prev_group = end_point[0].n_group_number;
    let prev_eqnum = end_point[0].n_equ_number;
    let (mut j, mut k, mut n_zero_equ) = (0usize, 0usize, 0usize);
    for e in end_point.iter() {
        if e.n_group_number != 0 && e.n_group_number < least_group {
            least_group = e.n_group_number;
        }
        j += usize::from(prev_group == e.n_group_number);
        k += usize::from(prev_eqnum == e.n_equ_number);
        n_zero_equ += usize::from(e.n_equ_number == 0);
    }
    if j == n_end && prev_group != 0 && k == n_end {
        return 0;
    }
    let mut n_new_tgroups = 0usize;
    let mut group_new_number = vec![0u16; n_end + 1];
    let mut group_number = vec![0u16; n_end + 1];
    if n_zero_equ == 0 {
        for e in end_point.iter_mut() {
            let group = e.n_equ_number;
            if group >= next_group {
                let jj = (0..n_new_tgroups)
                    .find(|&q| group == group_new_number[q])
                    .unwrap_or(n_new_tgroups);
                if jj == n_new_tgroups {
                    group_new_number[jj] = group;
                    n_new_tgroups += 1;
                }
                e.n_equ_number = next_group + jj as u16;
            }
        }
    } else if n_zero_equ == n_end {
        if least_group == next_group {
            n_new_tgroups = 1;
        }
        for e in end_point.iter_mut() {
            e.n_equ_number = least_group;
        }
    } else {
        return -1;
    }
    if n_new_tgroups != 0 {
        if num_t + n_new_tgroups > max_num_t {
            return -1;
        }
        for i in 0..n_new_tgroups {
            t_group[num_t + i] = TGroup::default();
            t_group[num_t + i].n_group_number = next_group + i as u16;
        }
    }
    let mut n_groups = 0usize;
    for i in 0..n_end {
        let group = end_point[i].n_group_number;
        if group != 0 {
            if group == end_point[i].n_equ_number {
                continue;
            }
            let jj = (0..n_groups)
                .find(|&q| group == group_number[q])
                .unwrap_or(n_groups);
            if jj < n_groups {
                if end_point[i].n_equ_number != group_new_number[jj] {
                    return -1;
                }
            } else {
                group_number[jj] = group;
                group_new_number[jj] = end_point[i].n_equ_number;
                n_groups += 1;
            }
        } else {
            let group = end_point[i].n_equ_number;
            let mut jj;
            if group >= next_group {
                jj = num_t + (group - next_group) as usize;
            } else {
                jj = 0;
                // old t-group
                if jj >= num_t || group != t_group[jj].n_group_number {
                    jj = match (0..num_t).find(|&q| group == t_group[q].n_group_number) {
                        Some(q) => q,
                        None => return -1,
                    };
                }
            }
            t_group[jj].n_num_endpoints += 1;
            for q in 0..t_group[jj].num.len() {
                t_group[jj].num[q] += end_point[i].num[q];
            }
            for q in 0..t_group[jj].num_da.len() {
                t_group[jj].num_da[q] += end_point[i].num_da[q];
            }
            at[end_point[i].n_atom_number as usize].endpoint = group;
            difference += 1;
        }
    }
    difference += n_groups as i32;
    num_t += n_new_tgroups;
    if difference == 0 {
        return 0;
    }
    if n_groups != 0 {
        next_group = 0;
        for g in &t_group[..num_t] {
            if next_group < g.n_group_number {
                next_group = g.n_group_number;
            }
        }
    }
    for i in 0..n_groups {
        let group1 = group_number[i];
        let group2 = group_new_number[i];
        let (mut i1, mut i2) = (None, None);
        let mut jj = 0;
        while jj < num_t && (i1.is_none() || i2.is_none()) {
            if i1.is_none() && group1 == t_group[jj].n_group_number {
                i1 = Some(jj);
            }
            if i2.is_none() && group2 == t_group[jj].n_group_number {
                i2 = Some(jj);
            }
            jj += 1;
        }
        let (Some(i1), Some(i2)) = (i1, i2) else {
            return -1;
        };
        let src = t_group[i1];
        for q in 0..src.num.len() {
            t_group[i2].num[q] += src.num[q];
        }
        for q in 0..src.num_da.len() {
            t_group[i2].num_da[q] += src.num_da[q];
        }
        t_group[i2].n_num_endpoints += src.n_num_endpoints;
        num_t -= 1;
        if num_t > i1 {
            t_group.copy_within(i1 + 1..num_t + 1, i1);
        }
    }
    if n_groups != 0 {
        let mut new_tg_number = vec![0u16; next_group as usize + 1];
        for (i, g) in t_group[..num_t].iter().enumerate() {
            new_tg_number[g.n_group_number as usize] = i as u16 + 1;
        }
        for jj in 0..n_groups {
            let old = group_number[jj] as usize;
            let new = group_new_number[jj] as usize;
            if new_tg_number[old] == 0 && new_tg_number[new] != 0 {
                new_tg_number[old] = new_tg_number[new];
            } else {
                return -1;
            }
        }
        for g in t_group[..num_t].iter_mut() {
            g.n_group_number = new_tg_number[g.n_group_number as usize];
        }
        for a in at.iter_mut().take(num_atoms as usize) {
            let group = a.endpoint;
            if group != 0 {
                a.endpoint = new_tg_number[group as usize];
                if a.endpoint == 0 || next_group <= new_tg_number[group as usize] {
                    return -1;
                }
            }
        }
    }
    tgi.t_group_number = vec![0u16; 2 * max_num_t];
    for (i, g) in t_group[..num_t].iter().enumerate() {
        if g.n_num_endpoints != 0 && g.n_group_number != 0 {
            tgi.t_group_number[g.n_group_number as usize] = i as u16 + 1;
        }
    }
    if let Some(p) = p {
        // ALWAYS_ADD_TG_ON_THE_FLY is 1
        let mut new_tgi = TGroupInfo {
            num_t_groups: num_t,
            t_group: Some(t_group[..num_t].to_vec()),
            ..Default::default()
        };
        use super::ichiparm::{
            TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_PT_06_00, TG_FLAG_PT_13_00, TG_FLAG_PT_16_00,
            TG_FLAG_PT_18_00, TG_FLAG_PT_22_00, TG_FLAG_PT_39_00,
        };
        new_tgi.b_taut_flags |= tgi.b_taut_flags
            & (TG_FLAG_KETO_ENOL_TAUT
                | TG_FLAG_PT_22_00
                | TG_FLAG_PT_16_00
                | TG_FLAG_PT_06_00
                | TG_FLAG_PT_39_00
                | TG_FLAG_PT_13_00
                | TG_FLAG_PT_18_00);
        let ret_bns = re_init_bn_struct(p, at, num_atoms as usize, false);
        if is_bns_error(ret_bns) {
            return ret_bns;
        }
        if u64::from(p.pb_taut_flags) & super::ichiparm::TG_FLAG_MOVE_POS_CHARGES != 0 {
            let r = add_c_groups_to_bn_struct(p, at, num_atoms, cgi);
            if is_bns_error(r) {
                return r;
            }
        }
        let r = add_t_groups_to_bn_struct(p, at, num_atoms, Some(&mut new_tgi));
        if is_bns_error(r) {
            return r;
        }
    }
    *pnum_t = num_t;
    tgi.num_t_groups = num_t;
    difference
}
