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
    pub b_taut_flags_done: u64,
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

/// `GetNeutralRepsIfNeeded`: if `*pri` and `*prj` are in one charge group
/// that has only one neutral atom, replaces them with neutral
/// representatives of their t-groups. (As in the C code, only the first
/// charge group is examined: the loop breaks after its first pass.)
pub fn get_neutral_reps_if_needed(
    pri: &mut u16,
    prj: &mut u16,
    at: &[InpAtom],
    num_atoms: i32,
    end_point: &[TEndpoint],
    cgi: Option<&CGroupInfo>,
) -> i32 {
    let mut ri = *pri as usize;
    let mut rj = *prj as usize;
    let c_point = at[ri].c_point;
    let Some(cgi) = cgi else { return 0 };
    if c_point == 0
        || at[rj].c_point != c_point
        || !(at[ri].charge == 1 || at[rj].charge == 1)
        || cgi.num_c_groups == 0
    {
        return 0;
    }
    let cg = cgi.c_group.as_ref().expect("c_group");
    if cg[0].n_group_number == c_point
        && i32::from(cg[0].num_c_points) - i32::from(cg[0].num[0]) < 2
    {
        let endpoint = at[rj].endpoint;
        if endpoint != 0 {
            for e in end_point {
                let r = e.n_atom_number as usize;
                if r == *prj as usize {
                    continue;
                }
                if at[r].endpoint != endpoint {
                    continue;
                }
                if at[r].c_point == 0 {
                    rj = r;
                    break;
                }
                if at[r].c_point != c_point && c_point == at[rj].c_point {
                    rj = r;
                }
            }
            if rj == *prj as usize {
                for i in 0..num_atoms as usize {
                    if at[i].endpoint != endpoint || i == *prj as usize {
                        continue;
                    }
                    if at[i].c_point == 0 {
                        rj = i;
                        break;
                    }
                    if at[i].c_point != c_point && c_point == at[rj].c_point {
                        rj = i;
                    }
                }
            }
        }
        let endpoint = at[ri].endpoint;
        if endpoint != 0 {
            for e in end_point {
                let r = e.n_atom_number as usize;
                if r == *pri as usize {
                    continue;
                }
                if at[r].endpoint != endpoint {
                    continue;
                }
                if at[r].c_point == 0 {
                    ri = r;
                    break;
                }
                if at[r].c_point != c_point
                    && c_point == at[ri].c_point
                    && at[r].c_point != at[rj].c_point
                {
                    ri = r;
                }
            }
            if ri == *pri as usize && at[rj].endpoint != 0 {
                for i in 0..num_atoms as usize {
                    if at[i].endpoint != endpoint || i == *pri as usize {
                        continue;
                    }
                    if at[i].c_point == 0 {
                        ri = i;
                        break;
                    }
                    if at[i].c_point != c_point
                        && c_point == at[ri].c_point
                        && at[i].c_point != at[rj].c_point
                    {
                        ri = i;
                    }
                }
            }
        }
    }
    *prj = rj as u16;
    *pri = ri as u16;
    0
}

/// `FindAccessibleEndPoints`: finds which of the endpoints are connected by
/// alternating paths (possibly through t-groups), sets their `n_equ_number`
/// and drops the endpoints that are not connected to any other. Returns the
/// number of kept endpoints or an error.
#[allow(clippy::too_many_arguments)]
pub fn find_accessible_end_points(
    end_point: &mut [TEndpoint],
    num_end_points: &mut usize,
    bond_pos: &mut [TBondPos],
    num_bond_pos: &mut usize,
    p: &mut super::ichi_bns::BnStruct,
    bd: &mut super::ichi_bns::BnData,
    at: &mut [InpAtom],
    num_atoms: i32,
    cgi: Option<&CGroupInfo>,
    taut_mode: i32,
) -> i32 {
    use super::ichi_bns::exists_any_alt_path;
    use super::ichierr::is_bns_error;
    if *num_end_points != *num_bond_pos {
        return 0;
    }
    let n_end = *num_end_points;
    let cap = n_end + 1;
    let mut representative = vec![0u16; cap];
    let mut equ = vec![0u16; cap];
    let mut ep_no = vec![0usize; cap];
    let mut n_diff = 0usize;
    for i in 0..n_end {
        end_point[i].n_equ_number = end_point[i].n_group_number;
        let cur = end_point[i].n_equ_number;
        if cur != 0 && !(0..n_diff).any(|j| equ[j] == cur) {
            representative[n_diff] = end_point[i].n_atom_number;
            equ[n_diff] = end_point[i].n_group_number;
            ep_no[n_diff] = i;
            n_diff += 1;
        }
    }
    let mut n_found_equ = 0;
    for i in 0..n_diff {
        for j in i + 1..n_diff {
            let (mut ri, mut rj) = (representative[i], representative[j]);
            get_neutral_reps_if_needed(&mut ri, &mut rj, at, num_atoms, &end_point[..n_end], cgi);
            let err = exists_any_alt_path(
                p,
                bd,
                at,
                num_atoms,
                i32::from(ri),
                i32::from(rj),
                taut_mode,
            );
            if is_bns_error(err) {
                return err;
            }
            if err == 0 {
                continue;
            }
            let cur = equ[i].min(equ[j]);
            let max = equ[i].max(equ[j]);
            for k in 0..n_diff {
                if equ[k] == max {
                    equ[k] = cur;
                    n_found_equ += 1;
                }
            }
            for e in end_point[..n_end].iter_mut() {
                if e.n_equ_number == max {
                    e.n_equ_number = cur;
                }
            }
        }
    }
    if n_found_equ != 0 {
        for i in 1..n_diff {
            if (0..i).any(|j| equ[j] == equ[i]) {
                equ[i] = 0;
            }
        }
        let mut j = 0;
        for i in 0..n_diff {
            if equ[i] != 0 {
                if i != j {
                    equ[j] = equ[i];
                    representative[j] = representative[i];
                    ep_no[j] = ep_no[i];
                }
                j += 1;
            }
        }
        n_diff = j;
    }
    // endpoints that have not been assigned to t-groups
    let mut j = n_diff;
    let mut equ_all = equ.clone();
    equ_all.resize(cap + n_end, 0);
    representative.resize(cap + n_end, 0);
    ep_no.resize(cap + n_end, 0);
    for i in 0..n_end {
        if end_point[i].n_equ_number != 0 {
            continue;
        }
        equ_all[j] = 0;
        representative[j] = end_point[i].n_atom_number;
        ep_no[j] = i;
        j += 1;
    }
    let equ = &mut equ_all;
    let n_tg = j;
    let mut max_equ = (num_atoms + 1) as u16;
    for i in 0..n_tg {
        for j in i + 1..n_tg {
            if (equ[i] != equ[j] && (i >= n_diff || j >= n_diff)) || (equ[i] == 0 && equ[j] == 0) {
                let (mut ri, mut rj) = (representative[i], representative[j]);
                get_neutral_reps_if_needed(
                    &mut ri,
                    &mut rj,
                    at,
                    num_atoms,
                    &end_point[..n_end],
                    cgi,
                );
                let err = exists_any_alt_path(
                    p,
                    bd,
                    at,
                    num_atoms,
                    i32::from(ri),
                    i32::from(rj),
                    taut_mode,
                );
                if is_bns_error(err) {
                    return err;
                }
                if err <= 0 {
                    continue;
                }
                if equ[i] != 0 && equ[j] != 0 {
                    let cur = equ[i].min(equ[j]);
                    let max = equ[i].max(equ[j]);
                    for k in 0..n_tg {
                        if equ[k] == max {
                            equ[k] = cur;
                        }
                    }
                    for e in end_point[..n_end].iter_mut() {
                        if e.n_equ_number == max {
                            e.n_equ_number = cur;
                        }
                    }
                } else if equ[i] != 0 {
                    equ[j] = equ[i];
                    end_point[ep_no[j]].n_equ_number = equ[i];
                } else if equ[j] != 0 {
                    equ[i] = equ[j];
                    end_point[ep_no[i]].n_equ_number = equ[j];
                } else {
                    equ[i] = max_equ;
                    equ[j] = max_equ;
                    end_point[ep_no[i]].n_equ_number = max_equ;
                    end_point[ep_no[j]].n_equ_number = max_equ;
                    max_equ += 1;
                }
            }
        }
    }
    // IGNORE_SINGLE_ENDPOINTS is 1: endpoints with no equivalent are dropped
    let mut j = 0;
    for i in 0..n_end {
        if end_point[i].n_equ_number != 0 {
            let same = end_point[..n_end]
                .iter()
                .filter(|e| e.n_equ_number == end_point[i].n_equ_number)
                .count();
            if same <= 1 {
                continue;
            }
            if i != j {
                end_point[j] = end_point[i];
                bond_pos[j] = bond_pos[i];
            }
            j += 1;
        }
    }
    *num_end_points = j;
    *num_bond_pos = j;
    j as i32
}

/// `TG_FLAG_ALLOW_NO_NEGTV_O_DONE`
pub const TG_FLAG_ALLOW_NO_NEGTV_O_DONE: u64 = 0x0000_0040;
/// `TG_FLAG_FOUND_SALT_CHARGES_DONE`
pub const TG_FLAG_FOUND_SALT_CHARGES_DONE: u64 = 0x0000_2000;
/// `FLAG_FORCE_SALT_TAUT` (`FIX_N_MINUS_NORN_BUG` is 0)
pub const FLAG_FORCE_SALT_TAUT: u64 = 0x0002 | 0x0010 | 0x0020;
/// `SALT_DONOR`
pub const SALT_DONOR: i32 = SALT_DONOR_NEG | SALT_DONOR_H;
/// `SALT_SELECTED`
pub const SALT_SELECTED: i32 = 32;

/// `S_GROUP_INFO`: the "salt" candidates
#[derive(Debug, Clone, Default)]
pub struct SGroupInfo {
    pub s_candidate: Vec<SCandidate>,
    pub max_num_candidates: i32,
    /// 0 = uninitialized, -1 = no candidates found
    pub num_candidates: i32,
    pub num_other_candidates: i32,
    pub num_p_only_candidates: i32,
}

/// `MarkSaltChargeGroups2`: finds acceptor and donor pairs of "salt"
/// candidates connected by alternating paths and merges the candidates into
/// one t-group. Returns 0 (no changes), 1 (new salt tautomerism), 2 (only
/// charges: no movable H), 5 (not only acidic atoms merged) or an error.
pub fn mark_salt_charge_groups2(
    at: &mut [InpAtom],
    num_atoms: i32,
    sgi: &mut SGroupInfo,
    tgi: &mut TGroupInfo,
    cgi: Option<&mut CGroupInfo>,
    p: &mut super::ichi_bns::BnStruct,
    bd: &mut super::ichi_bns::BnData,
) -> i32 {
    let mut n_other_changes = 0;
    let mut tot = 0;
    if sgi.s_candidate.is_empty() || sgi.max_num_candidates <= 0 {
        return 0;
    }
    mark_salt_charge_groups2_body(
        at,
        num_atoms,
        sgi,
        tgi,
        cgi,
        p,
        bd,
        &mut n_other_changes,
        &mut tot,
    );
    if n_other_changes != 0 && tot == 1 {
        tot = 5;
    }
    tot
}

#[allow(clippy::too_many_arguments)]
fn mark_salt_charge_groups2_body(
    at: &mut [InpAtom],
    num_atoms: i32,
    sgi: &mut SGroupInfo,
    tgi: &mut TGroupInfo,
    cgi: Option<&mut CGroupInfo>,
    p: &mut super::ichi_bns::BnStruct,
    bd: &mut super::ichi_bns::BnData,
    n_other_changes: &mut i32,
    tot: &mut i32,
) {
    use super::ichi_bns::{
        ALT_PATH_MODE_ADD2H_CHG, ALT_PATH_MODE_ADD2H_TST, ALT_PATH_MODE_REM2H_CHG,
        ALT_PATH_MODE_REM2H_TST, b_is_hard_rem_h_candidate, exists_alt_path,
    };
    use super::ichierr::{
        BNS_OUT_OF_RAM, BNS_PROGRAM_ERR, BNS_RADICAL_ERR, BNS_VERT_EDGE_OVFL, is_bns_error,
    };
    const DISABLE_CANDIDATE: i8 = 10;
    const ACCEPTOR_PAIR: u8 = 1;
    const DONOR_PAIR: u8 = 2;
    let max_cand = sgi.max_num_candidates as usize;
    let mut n_cand = sgi.num_candidates;
    let hard_added_removed = tgi.tni.b_normalization_flags & FLAG_FORCE_SALT_TAUT != 0;
    let mut s_subtype_all = 0;
    if n_cand <= -2 || tgi.t_group.is_none() {
        return;
    }
    // all candidates, including those of the other types
    n_cand = 0;
    let mut n_other_cand = 0;
    let mut s_subtype = 0i32;
    for i in 0..num_atoms as usize {
        let mut s_type = get_salt_charge_type(at, i, Some(tgi), &mut s_subtype);
        let mut ok = s_type == 0;
        if !ok {
            s_type = get_other_salt_charge_type(at, i, Some(tgi), &mut s_subtype, true);
            ok = s_type == 1;
        }
        if !ok {
            s_type = get_other_salt_type(at, i, &mut s_subtype);
            ok = s_type == 2;
        }
        if !ok && hard_added_removed {
            s_type = b_is_hard_rem_h_candidate(at, i, &mut s_subtype);
            ok = s_type == 4;
        }
        if ok {
            if n_cand as usize >= max_cand {
                *tot = BNS_VERT_EDGE_OVFL;
                return;
            }
            sgi.s_candidate[n_cand as usize] = SCandidate {
                atnumber: i as u16,
                type_: s_type as i8,
                subtype: s_subtype as i8,
                endpoint: at[i].endpoint,
            };
            n_cand += 1;
            n_other_cand += i32::from(s_type == 1);
            s_subtype_all |= s_subtype;
        }
    }
    let allow_no_neg = tgi.b_taut_flags & super::ichiparm::TG_FLAG_ALLOW_NO_NEGTV_O != 0
        || tgi.b_taut_flags_done & TG_FLAG_FOUND_SALT_CHARGES_DONE != 0
        || tgi.tni.b_normalization_flags & FLAG_FORCE_SALT_TAUT != 0;
    let no_donor = if allow_no_neg {
        s_subtype_all & SALT_DONOR == 0
    } else {
        s_subtype_all & SALT_DONOR_NEG == 0 || n_other_cand == n_cand
    };
    if n_cand <= 1 || s_subtype_all & SALT_ACCEPTOR == 0 || no_donor {
        sgi.num_candidates = 0;
        return;
    }
    if s_subtype_all & SALT_DONOR_NEG == 0 {
        tgi.b_taut_flags_done |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
    }
    let n_cand = n_cand as usize;
    // mark redundant candidates: one candidate per t-group is left in
    let mut n_left = 0usize;
    for i in 0..n_cand {
        // oss-fuzz #45059 guard of the C code: reads the element after the last
        if n_cand < max_cand && sgi.s_candidate[n_cand].type_ == 2 {
            sgi.s_candidate[i].type_ -= DISABLE_CANDIDATE;
            n_left += 1;
            continue;
        }
        if sgi.s_candidate[i].endpoint != 0 {
            for j in (0..i).rev() {
                if sgi.s_candidate[i].endpoint == sgi.s_candidate[j].endpoint {
                    sgi.s_candidate[i].type_ -= DISABLE_CANDIDATE;
                    n_left += 1;
                    break;
                }
            }
        }
    }
    let n_left = n_cand - n_left;
    sgi.num_candidates = 0;
    sgi.s_candidate[..n_cand].sort_by(|a, b| comp_candidates(a, b).cmp(&0));
    if n_left == 0 {
        *tot = BNS_OUT_OF_RAM;
        return;
    }
    let mut c_pair = vec![0u8; n_left * n_left];
    let pair = |c: &Vec<u8>, a: usize, b: usize| c[a + b * n_left];
    let (mut n_donor_pairs, mut n_acceptor_pairs) = (0, 0);
    let (mut fi, mut fj) = (0usize, 0usize);
    let mut found_two = false;
    'search: for i in 0..n_left {
        let (mut cur_donor, mut cur_acceptor) = (0, 0);
        for j in 0..=i {
            if i == j && sgi.s_candidate[i].endpoint == 0 {
                continue;
            }
            let (si, sj) = (sgi.s_candidate[i], sgi.s_candidate[j]);
            if i32::from(si.subtype) & SALT_ACCEPTOR != 0
                && i32::from(sj.subtype) & SALT_ACCEPTOR != 0
            {
                let ret = exists_alt_path(
                    p,
                    bd,
                    None,
                    at,
                    num_atoms,
                    i32::from(si.atnumber),
                    i32::from(sj.atnumber),
                    ALT_PATH_MODE_ADD2H_TST,
                );
                if ret != 0 {
                    if is_bns_error(ret) {
                        *tot = ret;
                        return;
                    }
                    if ret & 1 != 0 {
                        if (ret & !3) >> 2 != 0 {
                            *tot = BNS_RADICAL_ERR;
                            return;
                        }
                        c_pair[i + j * n_left] |= ACCEPTOR_PAIR;
                    }
                }
            }
            if i32::from(si.subtype) & SALT_DONOR != 0 && i32::from(sj.subtype) & SALT_DONOR != 0 {
                let ret = exists_alt_path(
                    p,
                    bd,
                    None,
                    at,
                    num_atoms,
                    i32::from(si.atnumber),
                    i32::from(sj.atnumber),
                    ALT_PATH_MODE_REM2H_TST,
                );
                if ret != 0 {
                    if is_bns_error(ret) {
                        *tot = ret;
                        return;
                    }
                    if ret & 1 != 0 {
                        if (ret & !3) >> 2 != 0 {
                            *tot = BNS_RADICAL_ERR;
                            return;
                        }
                        c_pair[i + j * n_left] |= DONOR_PAIR;
                    }
                }
            }
            if pair(&c_pair, i, j) & ACCEPTOR_PAIR != 0 {
                cur_acceptor += 1;
                if n_donor_pairs != 0 {
                    for i1 in 0..i {
                        for j1 in 0..=i1 {
                            if j1 != j && i1 != j && pair(&c_pair, i1, j1) & DONOR_PAIR != 0 {
                                found_two = true;
                                fi = i;
                                fj = j;
                                break 'search;
                            }
                        }
                    }
                }
            }
            if pair(&c_pair, i, j) & DONOR_PAIR != 0 {
                cur_donor += 1;
                if n_acceptor_pairs != 0 {
                    for i1 in 0..i {
                        for j1 in 0..=i1 {
                            if j1 != j && i1 != j && pair(&c_pair, i1, j1) & ACCEPTOR_PAIR != 0 {
                                found_two = true;
                                fi = i;
                                fj = j;
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
        n_donor_pairs += cur_donor;
        n_acceptor_pairs += cur_acceptor;
    }
    if !found_two {
        return;
    }
    // both a donor and an acceptor pair have been found
    let (i1, j1) = (fi, fj);
    let mut n_marked = 0usize;
    for i in 0..n_left {
        for j in 0..=i {
            let already = i < i1 || (i == i1 && j <= j1);
            for (st_mask, mode, flag) in [
                (SALT_ACCEPTOR, ALT_PATH_MODE_ADD2H_CHG, ACCEPTOR_PAIR),
                (SALT_DONOR, ALT_PATH_MODE_REM2H_CHG, DONOR_PAIR),
            ] {
                let run = (already && pair(&c_pair, i, j) & flag != 0) || !already;
                if !run {
                    continue;
                }
                let (si, sj) = (sgi.s_candidate[i], sgi.s_candidate[j]);
                if i32::from(si.subtype) & st_mask != 0 && i32::from(sj.subtype) & st_mask != 0 {
                    let ret = exists_alt_path(
                        p,
                        bd,
                        None,
                        at,
                        num_atoms,
                        i32::from(si.atnumber),
                        i32::from(sj.atnumber),
                        mode,
                    );
                    if ret != 0 {
                        if is_bns_error(ret) {
                            *tot = ret;
                            return;
                        }
                        if ret & 1 != 0 {
                            if (ret & !3) >> 2 != 0 {
                                *tot = BNS_RADICAL_ERR;
                                return;
                            }
                            c_pair[i + j * n_left] |= flag;
                            for k in [i, j] {
                                if i32::from(sgi.s_candidate[k].subtype) & SALT_SELECTED == 0 {
                                    sgi.s_candidate[k].subtype |= SALT_SELECTED as i8;
                                    n_marked += 1;
                                    if sgi.s_candidate[k].endpoint == 0
                                        && sgi.s_candidate[k].type_ != 0
                                    {
                                        *n_other_changes += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    let mut bt_group_has_neg_only = true;
    let mut i1 = 0usize;
    if n_marked != 0 {
        let mut end_points = vec![TEndpoint::default(); n_marked];
        let mut j = 0usize;
        for i in 0..n_left {
            if i32::from(sgi.s_candidate[i].subtype) & SALT_SELECTED != 0 {
                sgi.s_candidate[i].subtype ^= SALT_SELECTED as i8;
                if j < n_marked {
                    i1 = sgi.s_candidate[i].atnumber as usize;
                    add_end_point(&mut end_points[j], at, i1);
                }
                j += 1;
            }
        }
        if j != n_marked {
            *tot = BNS_PROGRAM_ERR;
            return;
        }
        let mut ret = register_end_points(tgi, &mut end_points, at, num_atoms, cgi, Some(p));
        if ret == -1 {
            ret = BNS_PROGRAM_ERR;
        }
        if ret < 0 {
            *tot = ret;
            return;
        }
        *tot += i32::from(ret > 0);
        for i in n_left..n_cand {
            sgi.s_candidate[i].type_ += DISABLE_CANDIDATE;
            let j1 = sgi.s_candidate[i].atnumber as usize;
            if at[j1].endpoint == at[i1].endpoint
                && sgi.s_candidate[i].endpoint == 0
                && sgi.s_candidate[i].type_ != 0
            {
                *n_other_changes += 1;
            }
        }
        bt_group_has_neg_only = false;
        let tg = tgi.t_group.as_ref().unwrap();
        for g in &tg[..tgi.num_t_groups] {
            if g.n_group_number == at[i1].endpoint && g.num[0] == g.num[1] {
                bt_group_has_neg_only = true;
                break;
            }
        }
    }
    *tot = i32::from(*tot > 0);
    if *tot != 0 && bt_group_has_neg_only {
        *tot = 2;
    }
}

/// `MarkSaltChargeGroups` (`SALT_WITH_PROTONS == 1`): looks for alternating
/// paths between "salt" candidates (H or (-) donors and acceptors), one
/// pair at a time, and merges the endpoints into t-groups. Returns the
/// number of changes or an error.
pub fn mark_salt_charge_groups(
    at: &mut [InpAtom],
    num_atoms: i32,
    sgi: &mut SGroupInfo,
    tgi: &mut TGroupInfo,
    mut cgi: Option<&mut CGroupInfo>,
    p: &mut super::ichi_bns::BnStruct,
    bd: &mut super::ichi_bns::BnData,
) -> i32 {
    use super::ichi_bns::{ALT_PATH_MODE_4_SALT, exists_alt_path};
    use super::ichierr::{BNS_PROGRAM_ERR, BNS_RADICAL_ERR, BNS_VERT_EDGE_OVFL, is_bns_error};
    if sgi.s_candidate.is_empty() || sgi.max_num_candidates <= 0 {
        return 0;
    }
    if sgi.num_candidates <= -1 || tgi.t_group.is_none() {
        return 0;
    }
    let max_cand = sgi.max_num_candidates as usize;
    let mut n_cand = 0usize;
    let mut n_other = 0;
    let mut s_subtype_all = 0;
    let mut s_subtype = 0i32;
    let mut num_tested = 0;
    let mut n_changes = 0;
    let mut tot = 0;
    let mut err = 0;
    'body: {
        for i in 0..num_atoms as usize {
            let mut s_type = get_salt_charge_type(at, i, Some(tgi), &mut s_subtype);
            let mut ok = s_type == 0;
            if !ok {
                s_type = get_other_salt_charge_type(at, i, Some(tgi), &mut s_subtype, true);
                ok = s_type == 1;
            }
            if !ok {
                s_type = get_other_salt_type(at, i, &mut s_subtype);
                ok = s_type == 2;
            }
            if ok {
                if n_cand >= max_cand {
                    err = BNS_VERT_EDGE_OVFL;
                    break 'body;
                }
                sgi.s_candidate[n_cand] = SCandidate {
                    atnumber: i as u16,
                    type_: s_type as i8,
                    subtype: s_subtype as i8,
                    endpoint: at[i].endpoint,
                };
                n_cand += 1;
                n_other += i32::from(s_type == 1);
                s_subtype_all |= s_subtype;
            }
        }
        let cond1 = s_subtype_all & SALT_ACCEPTOR != 0;
        let cond2a = tgi.b_taut_flags & super::ichiparm::TG_FLAG_ALLOW_NO_NEGTV_O != 0;
        let cond2b = tgi.b_taut_flags_done & TG_FLAG_FOUND_SALT_CHARGES_DONE != 0;
        let cond2c = tgi.tni.b_normalization_flags & FLAG_FORCE_SALT_TAUT != 0;
        let cond2 = if cond2a || cond2b || cond2c {
            s_subtype_all & (SALT_DONOR_NEG | SALT_DONOR_H) == 0
        } else {
            s_subtype_all & SALT_DONOR_NEG == 0 || n_other as usize == n_cand
        };
        if n_cand <= 1 || !cond1 || cond2 {
            sgi.num_candidates = -1;
            break 'body;
        }
        if s_subtype_all & SALT_DONOR_NEG == 0 {
            tgi.b_taut_flags_done |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
        }
        loop {
            n_changes = 0;
            for i1 in 0..n_cand {
                let j1 = sgi.s_candidate[i1].atnumber as usize;
                for i2 in i1 + 1..n_cand {
                    let j2 = sgi.s_candidate[i2].atnumber as usize;
                    if at[j1].endpoint != 0 && at[j1].endpoint == at[j2].endpoint {
                        continue;
                    }
                    for j in 0..2 {
                        let (ii1, ii2, jj1, jj2) = if j != 0 {
                            (i2, i1, j2, j1)
                        } else {
                            (i1, i2, j1, j2)
                        };
                        if i32::from(sgi.s_candidate[ii1].subtype) & (SALT_DONOR_NEG | SALT_DONOR_H)
                            != 0
                            && i32::from(sgi.s_candidate[ii2].subtype) & SALT_ACCEPTOR != 0
                        {
                            let ret = exists_alt_path(
                                p,
                                bd,
                                None,
                                at,
                                num_atoms,
                                jj2 as i32,
                                jj1 as i32,
                                ALT_PATH_MODE_4_SALT,
                            );
                            num_tested += 1;
                            if is_bns_error(ret) {
                                err = ret;
                                break 'body;
                            }
                            if ret & 1 != 0 {
                                let n_delta = (ret & !3) >> 2;
                                n_changes += ret & 2;
                                let mut end_points = [TEndpoint::default(); 2];
                                add_end_point(&mut end_points[0], at, jj1);
                                add_end_point(&mut end_points[1], at, jj2);
                                let mut ret = register_end_points(
                                    tgi,
                                    &mut end_points,
                                    at,
                                    num_atoms,
                                    cgi.as_deref_mut(),
                                    Some(p),
                                );
                                if ret == -1 {
                                    ret = BNS_PROGRAM_ERR;
                                }
                                if ret < 0 {
                                    err = ret;
                                    break 'body;
                                }
                                if n_delta != 0 {
                                    err = BNS_RADICAL_ERR;
                                    break 'body;
                                }
                                n_changes += i32::from(ret > 0);
                                break;
                            }
                        }
                    }
                }
            }
            tot += n_changes;
            if !(num_tested != 0 && n_changes != 0) {
                break;
            }
        }
    }
    if err == 0 {
        tot += n_changes;
        if sgi.num_candidates == 0 {
            sgi.num_candidates = if num_tested != 0 { n_cand as i32 } else { -1 };
        }
    } else {
        tot = err;
    }
    tot
}

pub const SALT_DONOR_ALL: i32 = SALT_DONOR_NEG | SALT_DONOR_H | SALT_P_ACCEPTOR | SALT_P_DONOR;
pub const SALT_DONOR_NEG2: i32 = SALT_DONOR_NEG | SALT_P_ACCEPTOR;
pub const SALT_DONOR_H2: i32 = SALT_DONOR_H | SALT_P_DONOR;

/// `MergeSaltTautGroups`: collects the acidic donors/acceptors (tautomeric
/// or not) and merges their t-groups into one. Returns the number of
/// changes or an error.
pub fn merge_salt_taut_groups(
    at: &mut [InpAtom],
    num_atoms: i32,
    sgi: &mut SGroupInfo,
    tgi: &mut TGroupInfo,
    cgi: Option<&mut CGroupInfo>,
    p: &mut super::ichi_bns::BnStruct,
) -> i32 {
    use super::ichi_bns::{b_has_acidic_hydrogen, b_has_acidic_minus};
    use super::ichierr::{BNS_PROGRAM_ERR, BNS_VERT_EDGE_OVFL};
    use super::inpdef::MAX_ATOMS;
    if sgi.s_candidate.is_empty() || tgi.t_group.is_none() || cgi.is_none() {
        return 0;
    }
    let max_cand = sgi.max_num_candidates;
    let mut cur_tg: i32 = MAX_ATOMS as i32;
    let mut s_subtype_all = 0;
    let mut n_cand = 0i32;
    for i in 0..num_atoms as usize {
        let mut s_subtype = 0i32;
        let mut s_type = get_salt_charge_type(at, i, Some(tgi), &mut s_subtype);
        let mut ok = s_type == 0;
        if !ok {
            s_type = get_other_salt_charge_type(at, i, Some(tgi), &mut s_subtype, true);
            ok = s_type == 1;
        }
        if !ok {
            s_type = get_other_salt_type(at, i, &mut s_subtype);
            ok = s_type == 2;
        }
        if !ok && b_has_acidic_hydrogen(at, i) != 0 {
            s_type = 3;
            s_subtype = SALT_P_DONOR;
            ok = true;
        }
        if !ok && b_has_acidic_minus(at, i) != 0 {
            s_type = 3;
            s_subtype = SALT_P_ACCEPTOR;
            ok = true;
        }
        if !ok {
            continue;
        }
        if n_cand >= max_cand {
            return BNS_VERT_EDGE_OVFL;
        }
        if at[i].endpoint == 0 && do_not_merge_non_taut_atom(at, i) {
            continue;
        }
        if s_subtype & SALT_DONOR_ALL == 0
            || (s_subtype & SALT_ACCEPTOR != 0 && at[i].endpoint == 0)
        {
            continue;
        }
        sgi.s_candidate[n_cand as usize] = SCandidate {
            atnumber: i as u16,
            type_: s_type as i8,
            subtype: s_subtype as i8,
            endpoint: at[i].endpoint,
        };
        n_cand += 1;
        s_subtype_all |= s_subtype;
    }
    let forced = tgi.b_taut_flags & super::ichiparm::TG_FLAG_ALLOW_NO_NEGTV_O != 0
        || tgi.b_taut_flags_done & TG_FLAG_FOUND_SALT_CHARGES_DONE != 0
        || tgi.tni.b_normalization_flags & FLAG_FORCE_SALT_TAUT != 0;
    let lacks = s_subtype_all & SALT_DONOR_NEG2 == 0 || s_subtype_all & SALT_DONOR_H2 == 0;
    let bail = if forced {
        n_cand <= 1 || (lacks && tgi.num_t_groups == 0)
    } else {
        n_cand <= 1 || lacks
    };
    if bail {
        sgi.num_candidates = -1;
        return 0;
    }
    if s_subtype_all & SALT_DONOR_NEG2 == 0 {
        tgi.b_taut_flags_done |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
    }
    sgi.num_candidates = n_cand;

    let mut min_endpoints = 0;
    for i in 0..sgi.num_candidates as usize {
        let c = sgi.s_candidate[i];
        let iat = c.atnumber as usize;
        if i32::from(c.subtype) & SALT_ACCEPTOR != 0 && at[iat].endpoint == 0 {
            continue;
        }
        if i32::from(at[iat].endpoint) != cur_tg || at[iat].endpoint == 0 {
            min_endpoints += 1;
        }
        cur_tg = i32::from(at[iat].endpoint);
    }
    if min_endpoints <= 1 {
        return 0;
    }
    let mut end_points = vec![TEndpoint::default(); min_endpoints];
    cur_tg = MAX_ATOMS as i32;
    let mut j = 0;
    for i in 0..sgi.num_candidates as usize {
        let c = sgi.s_candidate[i];
        let iat = c.atnumber as usize;
        if i32::from(c.subtype) == SALT_ACCEPTOR && at[iat].endpoint == 0 {
            continue;
        }
        if i32::from(at[iat].endpoint) != cur_tg || at[iat].endpoint == 0 {
            add_end_point(&mut end_points[j], at, iat);
            j += 1;
        }
        cur_tg = i32::from(at[iat].endpoint);
    }
    let ret = register_end_points(tgi, &mut end_points[..j], at, num_atoms, cgi, Some(p));
    if ret == -1 { BNS_PROGRAM_ERR } else { ret }
}
