//! The InChI -> structure restoration, part 2 (`INCHI_BASE/src/ichirvr2.c`):
//! atom connections and 0D stereo from an InChI component, the bond types of
//! the 0D stereo bonds.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

#![allow(clippy::needless_late_init)]

use super::ichi::{InChI, InChIStereo};
use super::ichi_bns::{
    BNS_VERT_TYPE_AUX, BNS_VERT_TYPE_TGROUP, BnData, BnStruct, NO_VERTEX, RAD_SRCH_FROM_FICT,
    remove_rad_endpoints, set_rad_endpoints,
};
use super::ichiread::{RI_ERR_PROGR, RI_ERR_SYNTAX};
use super::ichiring::is_bond_in_nmax_memb_ring;
use super::ichirvr1::BfsQ;
use super::ichirvr1::CN_LIST;
use super::ichirvr1::{
    BnsTestResult, EL_TYPE_C, EL_TYPE_MASK, EL_TYPE_N, EL_TYPE_O, EL_TYPE_OST, EL_TYPE_P,
    EL_TYPE_PT, EL_TYPE_S, EdgeList, add_to_edge_list, alloc_edge_list, bond_flow_maxcap_minorder,
    comp_cc_cand, get_pva_atom_type, get_tgroup_info_from_inchi, remove_forbidden_edge_mask,
    run_bns_restore_once, run_bns_test_once, set_forbidden_edge_mask,
};
use super::ichirvr4::forbid_carbon_charge_edges;
use super::ichirvrs::{
    AllTcGroups, BNS_VT_C_POS, CN_BITS_SHIFT, CcCand, Srm, TCG_ME_FLOWER0, TCG_PLUS,
    TGRF_MINUS_FIRST, ValAt,
};
use super::ichisort::insertions_sort;
use super::ichister::{
    AB_PARITY_EVEN, AB_PARITY_ODD, MAX_CUMULENE_LEN, MAX_NUM_STEREO_BOND_NEIGH,
    MIN_NUM_STEREO_BOND_NEIGH, atom_parity_well_def, b_can_atom_be_middle_allene,
    b_can_atom_be_terminal_allene,
};
use super::ichitaut::TGroupInfo;
use super::ichitaut::is_centerpoint_elem;
use super::incomdef::RADICAL_SINGLET;
use super::incomdef::{
    BOND_TYPE_ALTERN, BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, NUM_H_ISOTOPES,
};
use super::inpdef::{
    BOND_TYPE_MASK, InpAtom, MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS, MAXVAL,
};
use super::util::{
    EL_NUMBER_H, get_el_valence, get_element_chemical_symbol, get_endpoint_valence, is_el_a_metal,
    n_no_metal_bonds_valence, n_no_metal_num_bonds,
};

/// `BOND_TYPE_STEREO`: `BOND_TYPE_SINGLE | BOND_MARK_STEREO`
pub const BOND_TYPE_STEREO: i32 = BOND_TYPE_SINGLE | 0x10;

/// `inp_ATOM_STEREO`: the 0D stereo that exists only in the Mobile-H layer
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct InpAtomStereo {
    pub p_orig_at_num: [u16; MAX_NUM_STEREO_ATOM_NEIGH],
    pub p_parity: i8,
    pub sb_ord: [i8; MAX_NUM_STEREO_BONDS],
    pub sn_ord: [i8; MAX_NUM_STEREO_BONDS],
    pub sb_parity: [i8; MAX_NUM_STEREO_BONDS],
    pub sn_orig_at_num: [u16; MAX_NUM_STEREO_BONDS],
}

/// the used part of `StrFromINChI`: the structure restored from one InChI
/// component
#[derive(Debug, Clone, Default)]
pub struct StrFromInChI {
    /// `pSrm`: the rules for the bonds to metal atoms
    pub srm: Srm,
    /// length = `num_atoms + num_deleted_h`
    pub at: Vec<InpAtom>,
    pub st: Option<Vec<InpAtomStereo>>,
    pub ti: TGroupInfo,
    pub num_atoms: i32,
    pub num_deleted_h: i32,
    pub charge: i8,
    pub b_isotopic: i8,
    /// `endpoint`: the endpoints of the original Fixed-H InChI
    pub endpoint: Option<Vec<u16>>,
    /// `pbfsq`: work arrays of the ring search
    pub pbfsq: BfsQ,
}

/// `CopyAt2St`
pub fn copy_at_to_st(at: &[InpAtom], st: &mut [InpAtomStereo], num_atoms: usize) {
    for i in 0..num_atoms {
        if at[i].p_parity != 0 {
            st[i].p_orig_at_num = at[i].p_orig_at_num;
            st[i].p_parity = at[i].p_parity;
        }
        if at[i].sb_parity[0] != 0 {
            st[i].sb_ord = at[i].sb_ord;
            st[i].sb_parity = at[i].sb_parity;
            st[i].sn_ord = at[i].sn_ord;
            st[i].sn_orig_at_num = at[i].sn_orig_at_num;
        }
    }
}

/// `CopySt2At`
pub fn copy_st_to_at(at: &mut [InpAtom], st: Option<&[InpAtomStereo]>, num_atoms: usize) {
    let Some(st) = st else { return };
    for i in 0..num_atoms {
        if st[i].p_parity != 0 {
            at[i].p_orig_at_num = st[i].p_orig_at_num;
            at[i].p_parity = st[i].p_parity;
        }
        if st[i].sb_parity[0] != 0 {
            at[i].sb_ord = st[i].sb_ord;
            at[i].sb_parity = st[i].sb_parity;
            at[i].sn_ord = st[i].sn_ord;
            at[i].sn_orig_at_num = st[i].sn_orig_at_num;
        }
    }
}

/// `is_in_the_list`: the position of `next` among the first `len` items
fn is_in_the_list(list: &[u16], next: u16, len: usize) -> Option<usize> {
    list[..len].iter().position(|&x| x == next)
}

/// `bFindCumuleneChain`
pub fn b_find_cumulene_chain(
    at: &[InpAtom],
    i1: u16,
    i2: u16,
    cumulene: &mut [u16],
    max_len: i32,
) -> bool {
    cumulene[0] = i1;
    let a1 = &at[i1 as usize];
    for i in 0..a1.valence as usize {
        let mut len: i32 = 0;
        let mut iat = i1;
        let mut nat = a1.neighbor[i];
        if len + 1 == max_len {
            if nat == i2 {
                len += 1;
                cumulene[len as usize] = nat;
                return true;
            }
            continue;
        }
        while at[nat as usize].valence == 2
            && at[nat as usize].num_h == 0
            && b_can_atom_be_middle_allene(at[nat as usize].elname_str(), 0, 0)
        {
            len += 1;
            cumulene[len as usize] = nat;
            nat = at[nat as usize].neighbor[usize::from(at[nat as usize].neighbor[0] == iat)];
            if len + 1 == max_len {
                if nat == i2 {
                    len += 1;
                    cumulene[len as usize] = nat;
                    return true;
                }
                break;
            }
            iat = cumulene[len as usize];
        }
    }
    false
}

/// `set_bond_type`
pub fn set_bond_type(at: &mut [InpAtom], i1: u16, i2: u16, b_type: i32) -> i32 {
    let p1 = is_in_the_list(
        &at[i1 as usize].neighbor,
        i2,
        at[i1 as usize].valence as usize,
    );
    let p2 = is_in_the_list(
        &at[i2 as usize].neighbor,
        i1,
        at[i2 as usize].valence as usize,
    );
    if let (Some(j1), Some(j2)) = (p1, p2) {
        let prev = i32::from(at[i1 as usize].bond_type[j1]);
        at[i1 as usize].bond_type[j1] = b_type as u8;
        at[i2 as usize].bond_type[j2] = b_type as u8;
        if prev != 0 && prev <= BOND_TYPE_TRIPLE && b_type != 0 && b_type <= BOND_TYPE_TRIPLE {
            let d = (b_type - prev) as i8;
            at[i1 as usize].chem_bonds_valence = at[i1 as usize].chem_bonds_valence.wrapping_add(d);
            at[i2 as usize].chem_bonds_valence = at[i2 as usize].chem_bonds_valence.wrapping_add(d);
        }
        return 0;
    }
    RI_ERR_SYNTAX
}

#[derive(Clone, Copy)]
struct Sb {
    ord: [i8; MAX_NUM_STEREO_BONDS],
    par: [i8; MAX_NUM_STEREO_BONDS],
    sn: [i8; MAX_NUM_STEREO_BONDS],
    orig: [u16; MAX_NUM_STEREO_BONDS],
}

fn load_sb(at: &[InpAtom], st: Option<&[InpAtomStereo]>, i: usize) -> Sb {
    match st {
        Some(s) => Sb {
            ord: s[i].sb_ord,
            par: s[i].sb_parity,
            sn: s[i].sn_ord,
            orig: s[i].sn_orig_at_num,
        },
        None => Sb {
            ord: at[i].sb_ord,
            par: at[i].sb_parity,
            sn: at[i].sn_ord,
            orig: at[i].sn_orig_at_num,
        },
    }
}

fn store_sb(at: &mut [InpAtom], st: Option<&mut [InpAtomStereo]>, i: usize, x: &Sb) {
    match st {
        Some(s) => {
            s[i].sb_ord = x.ord;
            s[i].sb_parity = x.par;
            s[i].sn_ord = x.sn;
            s[i].sn_orig_at_num = x.orig;
        }
        None => {
            at[i].sb_ord = x.ord;
            at[i].sb_parity = x.par;
            at[i].sn_ord = x.sn;
            at[i].sn_orig_at_num = x.orig;
        }
    }
}

/// `set_cumulene_0D_parity`
#[allow(clippy::too_many_arguments)]
pub fn set_cumulene_0d_parity(
    at: &mut [InpAtom],
    mut st: Option<&mut [InpAtomStereo]>,
    num_at: usize,
    idel_h1: usize,
    i1: usize,
    i2: usize,
    idel_h2: usize,
    mut parity: i32,
    len: i32,
) -> i32 {
    let mut cumulene = [0u16; MAX_CUMULENE_LEN as usize + 2];
    if !b_find_cumulene_chain(at, i1 as u16, i2 as u16, &mut cumulene, len) {
        return RI_ERR_SYNTAX;
    }
    let p1 = is_in_the_list(&at[i1].neighbor, cumulene[1], at[i1].valence as usize);
    let p2 = is_in_the_list(
        &at[i2].neighbor,
        cumulene[len as usize - 1],
        at[i2].valence as usize,
    );
    let (sb_ord_m1, sb_ord_m2) = match (p1, p2) {
        (Some(a), Some(b)) => (a as i32, b as i32),
        _ => return RI_ERR_PROGR,
    };
    let num_neigh1 = i32::from(at[i1].valence) + i32::from(at[i1].num_h);
    let num_neigh2 = i32::from(at[i2].valence) + i32::from(at[i2].num_h);
    if !(MIN_NUM_STEREO_BOND_NEIGH..=MAX_NUM_STEREO_BOND_NEIGH).contains(&num_neigh1)
        || !(MIN_NUM_STEREO_BOND_NEIGH..=MAX_NUM_STEREO_BOND_NEIGH).contains(&num_neigh2)
    {
        return RI_ERR_SYNTAX;
    }
    let mut x1 = load_sb(at, st.as_deref(), i1);
    let mut x2 = load_sb(at, st.as_deref(), i2);
    let (mut m1, mut m2) = (0usize, 0usize);
    let (mut k1, mut k2) = (false, false);
    while m1 < MAX_NUM_STEREO_BONDS && x1.par[m1] != 0 {
        k1 = i32::from(x1.ord[m1]) == sb_ord_m1;
        if k1 {
            break;
        }
        m1 += 1;
    }
    while m2 < MAX_NUM_STEREO_BONDS && x2.par[m2] != 0 {
        k2 = i32::from(x2.ord[m2]) == sb_ord_m2;
        if k2 {
            break;
        }
        m2 += 1;
    }
    if m1 == MAX_NUM_STEREO_BONDS || m2 == MAX_NUM_STEREO_BONDS {
        return RI_ERR_SYNTAX;
    }
    if k1 && k2 {
        return 0;
    }
    if k1 || k2 {
        return RI_ERR_SYNTAX;
    }
    x1.ord[m1] = sb_ord_m1 as i8;
    x2.ord[m2] = sb_ord_m2 as i8;
    x1.orig[m1] = at[idel_h1].orig_at_number;
    let mut ret = 0;
    if idel_h1 < num_at {
        match is_in_the_list(&at[i1].neighbor, idel_h1 as u16, at[i1].valence as usize) {
            Some(p) => x1.sn[m1] = p as i8,
            None => ret = RI_ERR_PROGR,
        }
    } else {
        x1.sn[m1] = -1;
    }
    if ret == 0 {
        x2.orig[m2] = at[idel_h2].orig_at_number;
        if idel_h2 < num_at {
            match is_in_the_list(&at[i2].neighbor, idel_h2 as u16, at[i2].valence as usize) {
                Some(p) => x2.sn[m2] = p as i8,
                None => ret = RI_ERR_PROGR,
            }
        } else {
            x2.sn[m2] = -1;
        }
    }
    if ret == 0 {
        let (parity1, parity2);
        if atom_parity_well_def(parity) {
            let num_inv = i32::from(num_neigh1 == MIN_NUM_STEREO_BOND_NEIGH)
                + i32::from(num_neigh2 == MIN_NUM_STEREO_BOND_NEIGH);
            if num_inv % 2 != 0 {
                parity = if parity == AB_PARITY_EVEN {
                    AB_PARITY_ODD
                } else {
                    AB_PARITY_EVEN
                };
            }
            parity1 = AB_PARITY_EVEN;
            parity2 = if parity == AB_PARITY_EVEN {
                AB_PARITY_EVEN
            } else {
                AB_PARITY_ODD
            };
        } else {
            parity1 = parity;
            parity2 = parity;
        }
        x1.par[m1] = parity1 as i8;
        x2.par[m2] = parity2 as i8;
    }
    store_sb(at, st.as_deref_mut(), i1, &x1);
    store_sb(at, st, i2, &x2);
    ret
}

/// `set_atom_0D_parity`
pub fn set_atom_0d_parity(
    at: &mut [InpAtom],
    st: Option<&mut [InpAtomStereo]>,
    num_at: usize,
    num_deleted_h: usize,
    i1: usize,
    parity: i32,
) -> i32 {
    if st.is_some() && at[i1].p_parity != 0 {
        return 0;
    }
    let mut p_orig = match st.as_deref() {
        Some(s) => s[i1].p_orig_at_num,
        None => at[i1].p_orig_at_num,
    };
    let mut m1 = 0usize;
    let tot = i32::from(at[i1].valence) + i32::from(at[i1].num_h);
    let ret = (|| {
        if tot == MAX_NUM_STEREO_ATOM_NEIGH as i32 - 1 {
            p_orig[m1] = at[i1].orig_at_number;
            m1 += 1;
        } else if tot != MAX_NUM_STEREO_ATOM_NEIGH as i32 {
            return RI_ERR_PROGR;
        }
        let m2 = m1 as i32 + (MAX_NUM_STEREO_ATOM_NEIGH as i32 - i32::from(at[i1].valence));
        if at[i1].num_h != 0 {
            let mut j = 0;
            while (m1 as i32) < m2 && j < num_deleted_h {
                if at[j + num_at].neighbor[0] as usize == i1 {
                    p_orig[m1] = at[j + num_at].orig_at_number;
                    m1 += 1;
                }
                j += 1;
            }
        }
        if m1 as i32 + i32::from(at[i1].valence) != MAX_NUM_STEREO_ATOM_NEIGH as i32 {
            return RI_ERR_PROGR;
        }
        for i in 0..at[i1].valence as usize {
            let m2 = at[i1].neighbor[i] as usize;
            p_orig[m1] = at[m2].orig_at_number;
            m1 += 1;
        }
        0
    })();
    // the C code writes through the pointer as it goes; on success both are set
    match st {
        Some(s) => {
            s[i1].p_orig_at_num = p_orig;
            if ret == 0 {
                s[i1].p_parity = parity as i8;
            }
        }
        None => {
            at[i1].p_orig_at_num = p_orig;
            if ret == 0 {
                at[i1].p_parity = parity as i8;
            }
        }
    }
    ret
}

/// `AddExplicitDeletedH`
pub fn add_explicit_deleted_h(
    at: &mut [InpAtom],
    jv: usize,
    num_at: usize,
    i_deleted_h: &mut i32,
    i_h: &mut usize,
    num_deleted_h: i32,
    b_two_stereo: bool,
) -> i32 {
    let mut tot_num_iso_h: i32 = at[jv].num_iso_h.iter().map(|&x| i32::from(x)).sum();
    let mut num_h = i32::from(at[jv].num_h);
    let mut iso_h = 0usize;
    if at[jv].at_type == 0 {
        return RI_ERR_PROGR;
    }
    if at[jv].at_type > 1 {
        for i in 0..*i_deleted_h as usize {
            if at[num_at + i].neighbor[0] as usize == jv {
                *i_h = num_at + i;
                return 0;
            }
        }
        return RI_ERR_PROGR;
    }
    *i_h = *i_deleted_h as usize + num_at;
    let mut num_iso_h = [0i8; NUM_H_ISOTOPES];
    num_iso_h.copy_from_slice(&at[jv].num_iso_h);
    while num_h != 0 && *i_deleted_h < num_deleted_h {
        let ch = num_at + *i_deleted_h as usize;
        let v = at[ch].valence as usize;
        at[ch].neighbor[v] = jv as u16;
        at[ch].bond_type[v] = BOND_TYPE_SINGLE as u8;
        at[ch].valence += 1;
        if num_h > tot_num_iso_h {
            num_h -= 1;
            if num_h != tot_num_iso_h {
                if b_two_stereo {
                    *i_deleted_h += 1;
                    continue;
                }
                return RI_ERR_SYNTAX;
            }
        } else {
            while iso_h < NUM_H_ISOTOPES && num_iso_h[iso_h] == 0 {
                iso_h += 1;
            }
            if iso_h < NUM_H_ISOTOPES {
                at[ch].iso_atw_diff = (iso_h + 1) as i8;
                num_h -= 1;
                tot_num_iso_h -= 1;
                num_iso_h[iso_h] -= 1;
                if num_iso_h[iso_h] != 0 {
                    return RI_ERR_SYNTAX;
                }
            } else {
                return RI_ERR_SYNTAX;
            }
        }
        *i_deleted_h += 1;
    }
    if num_h != 0 {
        return RI_ERR_SYNTAX;
    }
    at[jv].at_type += 1;
    0
}

fn invert_parity(parity: i32) -> i32 {
    if parity == AB_PARITY_EVEN {
        AB_PARITY_ODD
    } else if parity == AB_PARITY_ODD {
        AB_PARITY_EVEN
    } else {
        parity
    }
}

fn stereo_of(inchi: &InChI) -> Option<&InChIStereo> {
    let nz = |s: &&InChIStereo| s.n_number_of_stereo_bonds + s.n_number_of_stereo_centers != 0;
    inchi
        .stereo_isotopic
        .as_ref()
        .filter(nz)
        .or_else(|| inchi.stereo.as_ref().filter(nz))
}

/// is `n` an allene centre (two terminal allene neighbors with 3 bonds)
fn is_allene_center(at: &[InpAtom], n: usize) -> Option<(usize, usize)> {
    let a = &at[n];
    if a.valence == 2 && a.num_h == 0 && b_can_atom_be_middle_allene(a.elname_str(), 0, 0) {
        let jv = a.neighbor[0] as usize;
        if i32::from(at[jv].valence) + i32::from(at[jv].num_h) == 3
            && b_can_atom_be_terminal_allene(at[jv].elname_str(), 0, 0)
        {
            let jn = a.neighbor[1] as usize;
            if i32::from(at[jn].valence) + i32::from(at[jn].num_h) == 3
                && b_can_atom_be_terminal_allene(at[jn].elname_str(), 0, 0)
            {
                return Some((jv, jn));
            }
        }
    }
    None
}

struct DelH {
    i_deleted_h: i32,
    num_deleted_h: i32,
    two: bool,
}

/// sets the 0D parity of a stereo centre (an allene or an sp3 atom)
fn set_center_parity(
    at: &mut [InpAtom],
    st: Option<&mut [InpAtomStereo]>,
    num_atoms: usize,
    n_vertex: usize,
    parity: i32,
    d: &mut DelH,
) -> i32 {
    let mut idel = 0usize;
    if let Some((jv, jn)) = is_allene_center(at, n_vertex) {
        let idel_h1;
        let idel_h2;
        if at[jv].num_h != 0 {
            let r = add_explicit_deleted_h(
                at,
                jv,
                num_atoms,
                &mut d.i_deleted_h,
                &mut idel,
                d.num_deleted_h,
                d.two,
            );
            if r < 0 {
                return r;
            }
            idel_h1 = idel;
        } else {
            idel_h1 =
                at[jv].neighbor[usize::from(at[jv].neighbor[0] as usize == n_vertex)] as usize;
        }
        if at[jn].num_h != 0 {
            let r = add_explicit_deleted_h(
                at,
                jn,
                num_atoms,
                &mut d.i_deleted_h,
                &mut idel,
                d.num_deleted_h,
                d.two,
            );
            if r < 0 {
                return r;
            }
            idel_h2 = idel;
        } else {
            idel_h2 =
                at[jn].neighbor[usize::from(at[jn].neighbor[0] as usize == n_vertex)] as usize;
        }
        set_cumulene_0d_parity(at, st, num_atoms, idel_h1, jv, jn, idel_h2, parity, 2)
    } else {
        if at[n_vertex].num_h != 0 {
            let r = add_explicit_deleted_h(
                at,
                n_vertex,
                num_atoms,
                &mut d.i_deleted_h,
                &mut idel,
                d.num_deleted_h,
                d.two,
            );
            if r < 0 {
                return r;
            }
        }
        set_atom_0d_parity(
            at,
            st,
            num_atoms,
            d.num_deleted_h as usize,
            n_vertex,
            parity,
        )
    }
}

/// sets the 0D parity of a stereo bond or cumulene
fn set_bond_parity(
    at: &mut [InpAtom],
    st: Option<&mut [InpAtomStereo]>,
    num_atoms: usize,
    jv: usize,
    jn: usize,
    parity: i32,
    d: &mut DelH,
) -> i32 {
    let mut cumulene = [0u16; MAX_CUMULENE_LEN as usize + 2];
    let len;
    if is_in_the_list(&at[jv].neighbor, jn as u16, at[jv].valence as usize).is_none() {
        if !b_find_cumulene_chain(
            at,
            jv as u16,
            jn as u16,
            &mut cumulene,
            MAX_CUMULENE_LEN + 1,
        ) {
            return RI_ERR_SYNTAX;
        }
        len = MAX_CUMULENE_LEN + 1;
    } else {
        cumulene[0] = jv as u16;
        cumulene[1] = jn as u16;
        len = 1;
    }
    let mut idel = 0usize;
    let idel_h1;
    let idel_h2;
    if at[jv].num_h != 0 {
        let r = add_explicit_deleted_h(
            at,
            jv,
            num_atoms,
            &mut d.i_deleted_h,
            &mut idel,
            d.num_deleted_h,
            d.two,
        );
        if r < 0 {
            return r;
        }
        idel_h1 = idel;
    } else {
        idel_h1 = at[jv].neighbor[usize::from(at[jv].neighbor[0] == cumulene[1])] as usize;
    }
    if at[jn].num_h != 0 {
        let r = add_explicit_deleted_h(
            at,
            jn,
            num_atoms,
            &mut d.i_deleted_h,
            &mut idel,
            d.num_deleted_h,
            d.two,
        );
        if r < 0 {
            return r;
        }
        idel_h2 = idel;
    } else {
        idel_h2 =
            at[jn].neighbor[usize::from(at[jn].neighbor[0] == cumulene[len as usize - 1])] as usize;
    }
    set_cumulene_0d_parity(at, st, num_atoms, idel_h1, jv, jn, idel_h2, parity, len)
}

/// `RestoreAtomConnectionsSetStereo`: restores the atoms, connections, H,
/// t-groups and 0D stereo of one InChI component into `ps`; returns the number
/// of atoms or a negative error
pub fn restore_atom_connections_set_stereo(
    ps: &mut StrFromInChI,
    i_component: i32,
    i_at_no_offset: i32,
    inchi: &InChI,
    inchi_mobh: Option<&InChI>,
) -> i32 {
    let num_atoms = inchi.n_number_of_atoms;
    if num_atoms <= 0 {
        return 0;
    }
    let n = num_atoms as usize;
    let mut at = vec![InpAtom::default(); n];
    ps.num_atoms = num_atoms;
    ps.charge = inchi.n_total_charge as i8;
    let r = restore_inner(ps, &mut at, i_component, i_at_no_offset, inchi, inchi_mobh);
    ps.at = at;
    r
}

fn restore_inner(
    ps: &mut StrFromInChI,
    at_ref: &mut Vec<InpAtom>,
    i_component: i32,
    i_at_no_offset: i32,
    inchi: &InChI,
    inchi_mobh: Option<&InChI>,
) -> i32 {
    let num_atoms = inchi.n_number_of_atoms as usize;
    let mut at = std::mem::take(at_ref);
    let ret = restore_body(
        ps,
        &mut at,
        num_atoms,
        i_component,
        i_at_no_offset,
        inchi,
        inchi_mobh,
    );
    *at_ref = at;
    ret
}

fn restore_body(
    ps: &mut StrFromInChI,
    at: &mut Vec<InpAtom>,
    num_atoms: usize,
    i_component: i32,
    i_at_no_offset: i32,
    inchi: &InChI,
    inchi_mobh: Option<&InChI>,
) -> i32 {
    for i in 0..num_atoms {
        at[i].el_number = inchi.n_atom[i];
        let (r, sym) = get_element_chemical_symbol(i32::from(inchi.n_atom[i]));
        if r == -1 {
            return RI_ERR_PROGR;
        }
        at[i].set_elname(sym);
        at[i].orig_at_number = (i_at_no_offset + i as i32 + 1) as u16;
        at[i].orig_compt_at_numb = (i + 1) as u16;
        at[i].component = (i_component + 1) as u16;
        at[i].num_h = inchi.n_num_h[i];
    }
    // connections
    let mut n_vertex = inchi.n_conn_table[0] as i32 - 1;
    for i in 1..inchi.len_conn_table as usize {
        let n_neigh = inchi.n_conn_table[i] as i32 - 1;
        if n_neigh < n_vertex {
            let (v, nb) = (n_vertex as usize, n_neigh as usize);
            let jv = at[v].valence as usize;
            at[v].valence += 1;
            at[v].neighbor[jv] = nb as u16;
            at[v].bond_type[jv] = BOND_TYPE_SINGLE as u8;
            at[v].chem_bonds_valence += BOND_TYPE_SINGLE as i8;
            let jn = at[nb].valence as usize;
            at[nb].valence += 1;
            at[nb].neighbor[jn] = v as u16;
            at[nb].bond_type[jn] = BOND_TYPE_SINGLE as u8;
            at[nb].chem_bonds_valence += BOND_TYPE_SINGLE as i8;
        } else {
            n_vertex = n_neigh;
            if n_vertex >= num_atoms as i32 {
                return RI_ERR_PROGR;
            }
        }
    }
    // isotopic atoms
    if !inchi.isotopic_atom.is_empty() && inchi.n_number_of_isotopic_atoms != 0 {
        for ia in inchi
            .isotopic_atom
            .iter()
            .take(inchi.n_number_of_isotopic_atoms as usize)
        {
            let v = ia.n_atom_number as i32 - 1;
            if v < num_atoms as i32 {
                let a = &mut at[v as usize];
                a.iso_atw_diff = ia.n_iso_difference as i8;
                a.num_iso_h[0] = ia.n_num_h as i8;
                a.num_iso_h[1] = ia.n_num_d as i8;
                a.num_iso_h[2] = ia.n_num_t as i8;
            }
        }
        ps.b_isotopic |= 1;
    }
    // tautomeric groups
    let r = get_tgroup_info_from_inchi(&mut ps.ti, Some(&mut at[..]), None, Some(inchi));
    if r != 0 {
        return r;
    }

    let p_stereo = stereo_of(inchi);
    let mut p_stereo2: Option<&InChIStereo> = None;
    if let Some(m) = inchi_mobh.filter(|m| m.n_number_of_atoms != 0) {
        p_stereo2 = stereo_of(m);
    }
    let (mut num_stereo_bonds2, mut num_stereo_centers2) = (0, 0);
    let mut num_deleted_h = 0i32;
    let mut d = DelH {
        i_deleted_h: 0,
        num_deleted_h: 0,
        two: p_stereo2.is_some(),
    };

    if p_stereo.is_some() || p_stereo2.is_some() {
        let len1 = p_stereo.map_or(0, |s| s.n_number_of_stereo_centers as usize);
        let len2 = p_stereo2.map_or(0, |s| s.n_number_of_stereo_centers as usize);
        let (mut i, mut i2) = (0usize, 0usize);
        while i < len1 || i2 < len2 {
            let n_vertex;
            if i < len1 && i2 < len2 {
                let diff = i32::from(p_stereo.unwrap().n_number[i])
                    - i32::from(p_stereo2.unwrap().n_number[i2]);
                if diff <= 0 {
                    n_vertex = p_stereo.unwrap().n_number[i] as usize - 1;
                    i += 1;
                    i2 += usize::from(diff == 0);
                } else {
                    n_vertex = p_stereo2.unwrap().n_number[i2] as usize - 1;
                    num_stereo_centers2 += 1;
                    i2 += 1;
                }
            } else if i < len1 {
                n_vertex = p_stereo.unwrap().n_number[i] as usize - 1;
                i += 1;
            } else {
                n_vertex = p_stereo2.unwrap().n_number[i2] as usize - 1;
                num_stereo_centers2 += 1;
                i2 += 1;
            }
            let targets: Vec<usize> = match is_allene_center(at, n_vertex) {
                Some((jv, jn)) => vec![jv, jn],
                None => vec![n_vertex],
            };
            for t in targets {
                if at[t].at_type == 0 && at[t].num_h != 0 {
                    num_deleted_h += i32::from(at[t].num_h);
                    at[t].at_type += 1;
                }
            }
        }
        let len1 = p_stereo.map_or(0, |s| s.n_number_of_stereo_bonds as usize);
        let len2 = p_stereo2.map_or(0, |s| s.n_number_of_stereo_bonds as usize);
        let (mut i, mut i2) = (0usize, 0usize);
        while i < len1 || i2 < len2 {
            let (n_vertex, n_neigh);
            if i < len1 && i2 < len2 {
                let (s1, s2) = (p_stereo.unwrap(), p_stereo2.unwrap());
                let diff = i32::from(s1.n_bond_atom1[i]) - i32::from(s2.n_bond_atom1[i2]);
                let diff2 = i32::from(s1.n_bond_atom2[i]) - i32::from(s2.n_bond_atom2[i2]);
                if diff < 0 || (diff == 0 && diff2 <= 0) {
                    n_vertex = s1.n_bond_atom1[i] as i32 - 1;
                    n_neigh = s1.n_bond_atom2[i] as i32 - 1;
                    i += 1;
                    i2 += usize::from(diff == 0 && diff2 == 0);
                } else {
                    n_vertex = s2.n_bond_atom1[i2] as i32 - 1;
                    n_neigh = s2.n_bond_atom2[i2] as i32 - 1;
                    num_stereo_bonds2 += 1;
                    i2 += 1;
                }
            } else if i < len1 {
                let s1 = p_stereo.unwrap();
                n_vertex = s1.n_bond_atom1[i] as i32 - 1;
                n_neigh = s1.n_bond_atom2[i] as i32 - 1;
                i += 1;
            } else {
                let s2 = p_stereo2.unwrap();
                n_vertex = s2.n_bond_atom1[i2] as i32 - 1;
                n_neigh = s2.n_bond_atom2[i2] as i32 - 1;
                num_stereo_bonds2 += 1;
                i2 += 1;
            }
            if n_vertex < num_atoms as i32 {
                let v = n_vertex as usize;
                if is_in_the_list(&at[v].neighbor, n_neigh as u16, at[v].valence as usize).is_none()
                {
                    let mut c = [0u16; MAX_CUMULENE_LEN as usize + 2];
                    if !b_find_cumulene_chain(
                        at,
                        v as u16,
                        n_neigh as u16,
                        &mut c,
                        MAX_CUMULENE_LEN + 1,
                    ) {
                        return RI_ERR_SYNTAX;
                    }
                }
                if at[v].at_type == 0 && at[v].num_h != 0 {
                    num_deleted_h += i32::from(at[v].num_h);
                    at[v].at_type += 1;
                }
            }
            let nn = n_neigh as usize;
            if at[nn].at_type == 0 && at[nn].num_h != 0 {
                num_deleted_h += i32::from(at[nn].num_h);
                at[nn].at_type += 1;
            }
        }
        if num_deleted_h != 0 {
            ps.num_deleted_h = num_deleted_h;
            at.resize(num_atoms + num_deleted_h as usize, InpAtom::default());
            for i in num_atoms..num_atoms + num_deleted_h as usize {
                at[i].set_elname("H");
                at[i].el_number = EL_NUMBER_H;
                at[i].orig_at_number = (i_at_no_offset + i as i32 + 1) as u16;
                at[i].orig_compt_at_numb = (i + 1) as u16;
                at[i].component = (i_component + 1) as u16;
            }
            d.num_deleted_h = num_deleted_h;
            for i in 0..num_atoms {
                if at[i].at_type == 1 {
                    let mut idel = 0usize;
                    let r = add_explicit_deleted_h(
                        at,
                        i,
                        num_atoms,
                        &mut d.i_deleted_h,
                        &mut idel,
                        num_deleted_h,
                        d.two,
                    );
                    if r < 0 {
                        return r;
                    }
                }
            }
        }
    }
    d.num_deleted_h = num_deleted_h;

    let mut st: Option<Vec<InpAtomStereo>> = None;
    if let Some(s) = p_stereo {
        let inverted = s.n_comp_inv2_abs == -1;
        for i in 0..s.n_number_of_stereo_centers as usize {
            let n_vertex = s.n_number[i] as usize - 1;
            let mut parity = i32::from(s.t_parity[i]);
            if inverted {
                parity = invert_parity(parity);
            }
            let r = set_center_parity(at, None, num_atoms, n_vertex, parity, &mut d);
            if r < 0 {
                return r;
            }
        }
        for i in 0..s.n_number_of_stereo_bonds as usize {
            let jv = s.n_bond_atom1[i] as usize - 1;
            let jn = s.n_bond_atom2[i] as usize - 1;
            let r = set_bond_parity(
                at,
                None,
                num_atoms,
                jv,
                jn,
                i32::from(s.b_parity[i]),
                &mut d,
            );
            if r < 0 {
                return r;
            }
        }
    }
    if num_stereo_centers2 + num_stereo_bonds2 != 0 {
        let mut v = vec![InpAtomStereo::default(); num_atoms];
        copy_at_to_st(at, &mut v, num_atoms);
        st = Some(v);
    }
    if num_stereo_centers2 != 0 {
        let s2 = p_stereo2.unwrap();
        let inverted = s2.n_comp_inv2_abs == -1;
        for i in 0..s2.n_number_of_stereo_centers as usize {
            let n_vertex = s2.n_number[i] as usize - 1;
            let mut parity = i32::from(s2.t_parity[i]);
            if at[n_vertex].p_parity != 0 {
                continue;
            }
            if inverted {
                parity = invert_parity(parity);
            }
            let r = set_center_parity(at, st.as_deref_mut(), num_atoms, n_vertex, parity, &mut d);
            if r < 0 {
                return r;
            }
        }
    }
    if num_stereo_bonds2 != 0 {
        let s2 = p_stereo2.unwrap();
        for i in 0..s2.n_number_of_stereo_bonds as usize {
            let jv = s2.n_bond_atom1[i] as usize - 1;
            let jn = s2.n_bond_atom2[i] as usize - 1;
            let r = set_bond_parity(
                at,
                st.as_deref_mut(),
                num_atoms,
                jv,
                jn,
                i32::from(s2.b_parity[i]),
                &mut d,
            );
            if r < 0 {
                return r;
            }
        }
    }
    ps.st = st;
    num_atoms as i32
}

/// `SetStereoBondTypeFor0DParity`
pub fn set_stereo_bond_type_for_0d_parity(at: &mut [InpAtom], i1: usize, m1: usize) -> i32 {
    let mut cumulene = [0u16; MAX_CUMULENE_LEN as usize + 3];
    let mut len = 0usize;
    let mut k1 = at[i1].sb_ord[m1] as usize;
    let mut n1 = i1;
    cumulene[len] = n1 as u16;
    len += 1;
    let found = loop {
        let n2 = at[n1].neighbor[k1] as usize;
        cumulene[len] = n2 as u16;
        len += 1;
        let mut m2 = 0;
        let mut hit = false;
        while m2 < MAX_NUM_STEREO_BONDS && at[n2].sb_parity[m2] != 0 {
            if n1 == at[n2].neighbor[at[n2].sb_ord[m2] as usize] as usize {
                hit = true;
                break;
            }
            m2 += 1;
        }
        if hit {
            break true;
        }
        if at[n2].num_h != 0 || at[n2].valence != 2 || at[n2].endpoint != 0 {
            break false;
        }
        k1 = usize::from(at[n2].neighbor[0] as usize == n1);
        n1 = n2;
        if !(at[n1].valence == 2
            && at[n1].num_h == 0
            && (len as i32) < MAX_CUMULENE_LEN + 2
            && b_can_atom_be_middle_allene(at[n1].elname_str(), at[n1].charge, at[n1].radical))
        {
            break false;
        }
    };
    if !found {
        return RI_ERR_SYNTAX;
    }
    let bond_type = if len == 2 {
        BOND_TYPE_STEREO
    } else {
        BOND_TYPE_DOUBLE
    };
    for j in 1..len {
        if set_bond_type(at, cumulene[j - 1], cumulene[j], bond_type) < 0 {
            return RI_ERR_PROGR;
        }
    }
    len as i32
}

/// `SetStereoBondTypesFrom0DStereo`
pub fn set_stereo_bond_types_from_0d_stereo(ps: &mut StrFromInChI, inchi: &InChI) -> i32 {
    let num_atoms = ps.num_atoms as usize;
    let at = &mut ps.at;
    if stereo_of(inchi).is_none() {
        return 0;
    }
    let mut num_stereo_bonds = 0;
    for i in 0..num_atoms {
        let mut j = 0;
        while j < MAX_NUM_STEREO_BONDS && at[i].sb_parity[j] != 0 {
            num_stereo_bonds += 1;
            let r = set_stereo_bond_type_for_0d_parity(at, i, j);
            if r < 0 {
                return r;
            }
            j += 1;
        }
    }
    if num_stereo_bonds == 0 {
        return 0;
    }
    let count = |a: &InpAtom| {
        let mut s = 0;
        let mut al = 0;
        for j in 0..a.valence as usize {
            s += i32::from(i32::from(a.bond_type[j]) == BOND_TYPE_STEREO);
            al += i32::from(i32::from(a.bond_type[j]) == BOND_TYPE_ALTERN);
        }
        (s, al)
    };
    for i in 0..num_atoms {
        let (s, al) = count(&at[i]);
        if s + al > 1 && s != 0 {
            for j in 0..at[i].valence as usize {
                if i32::from(at[i].bond_type[j]) == BOND_TYPE_STEREO {
                    let neigh = at[i].neighbor[j];
                    let r = set_bond_type(at, i as u16, neigh, BOND_TYPE_ALTERN);
                    if r < 0 {
                        return r;
                    }
                }
            }
        }
    }
    for i in 0..num_atoms {
        let (s, al) = count(&at[i]);
        if s == 0 && al != 0 {
            at[i].chem_bonds_valence += 1;
        } else if s == 1 {
            for j in 0..at[i].valence as usize {
                if i32::from(at[i].bond_type[j]) == BOND_TYPE_STEREO {
                    let neigh = at[i].neighbor[j];
                    let r = set_bond_type(at, i as u16, neigh, BOND_TYPE_DOUBLE);
                    if r < 0 {
                        return r;
                    }
                    at[i].chem_bonds_valence += 1;
                    at[neigh as usize].chem_bonds_valence += 1;
                }
            }
        } else if s + al != 0 {
            return RI_ERR_PROGR;
        }
    }
    0
}

/// `CopyBnsToAtom`: copies the bond orders, charges, radicals and the
/// tautomeric H and (-) of the solved BNS to the atoms of `ps`
pub fn copy_bns_to_atom(
    ps: &mut StrFromInChI,
    bns: &BnStruct,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    allow_zero_bond_order: bool,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let srm = ps.srm;
    let at = &mut ps.at;
    for i in 0..num_at {
        let pv = &bns.vert[i];
        let mut chem_bonds_valence = 0;
        for j in 0..at[i].valence as usize {
            let pe = &bns.edge[bns.iedge[pv.iedge + j] as usize];
            let (_, _, min_order, _) = bond_flow_maxcap_minorder(at, p_va, &srm, i, j);
            let mut bond_order = pe.flow + min_order;
            if !allow_zero_bond_order && bond_order == 0 {
                bond_order = 1;
            }
            chem_bonds_valence += bond_order;
            at[i].bond_type[j] = bond_order as u8;
        }
        at[i].chem_bonds_valence = chem_bonds_valence as i8;
        at[i].charge = p_va[i].c_init_charge;
        if p_va[i].n_c_minus_group_edge != 0 {
            let pe = &bns.edge[p_va[i].n_c_minus_group_edge as usize - 1];
            if pe.flow != 0 {
                at[i].charge = at[i].charge.wrapping_sub(pe.flow as i8);
            }
        }
        if p_va[i].n_c_plus_group_edge != 0 {
            let pe = &bns.edge[p_va[i].n_c_plus_group_edge as usize - 1];
            let c = pe.cap - pe.flow;
            if c != 0 {
                at[i].charge = at[i].charge.wrapping_add(c as i8);
            }
        }
        if pv.st_edge.cap > pv.st_edge.flow {
            at[i].radical = (RADICAL_SINGLET + (pv.st_edge.cap - pv.st_edge.flow)) as i8;
        }
    }
    for i in 0..bns.num_t_groups as usize {
        let g = &tcg.p_tcg[i];
        let mut num_h = i32::from(g.tg_num_h);
        let mut num_minus = i32::from(g.tg_num_minus);
        let minus_first = g.tg_restore_flags & TGRF_MINUS_FIRST != 0;
        let mut v_minus = NO_VERTEX;
        let pv = &bns.vert[num_at + i];
        if pv.type_ & BNS_VERT_TYPE_TGROUP == 0 {
            return RI_ERR_PROGR;
        }
        if g.tg_set_minus > 0 && num_minus > 0 {
            v_minus = i32::from(g.tg_set_minus) - 1;
            num_minus -= 1;
        }
        let order: Vec<usize> = if minus_first {
            (0..pv.num_adj_edges as usize).collect()
        } else {
            (0..pv.num_adj_edges as usize).rev().collect()
        };
        for j in order {
            let pe = &bns.edge[bns.iedge[pv.iedge + j] as usize];
            let v1 = pe.neighbor1 as usize;
            let mut num_at_add = pe.flow;
            if v1 as i32 == v_minus {
                if num_at_add != 0 {
                    at[v1].charge = -1;
                    num_at_add -= 1;
                } else {
                    num_minus += 1;
                }
                v_minus = NO_VERTEX;
            }
            if num_at_add > 0 {
                if num_minus != 0
                    && at[v1].charge == 0
                    && at[v1].valence == at[v1].chem_bonds_valence
                {
                    at[v1].charge -= 1;
                    num_at_add -= 1;
                    num_minus -= 1;
                }
                if num_at_add > 0 {
                    at[v1].num_h = at[v1].num_h.wrapping_add(num_at_add as i8);
                    num_h -= num_at_add;
                }
            }
            at[v1].endpoint = (i + 1) as u16;
        }
        if num_h + num_minus != pv.st_edge.cap - pv.st_edge.flow
            && (num_h != 0 || num_minus != 0 || v_minus != NO_VERTEX)
        {
            return RI_ERR_PROGR;
        }
    }
    0
}

/// `CheckBnsConsistency` (a no-op in the library build: `TARGET_API_LIB`)
pub fn check_bns_consistency() -> i32 {
    0
}

/// `MoveRadToAtomsAddCharges`: moves the radicals of the charge structure to
/// atoms and adds the charges that this requires. `at2` receives the atoms
/// made from the BNS (the work array of the C function).
#[allow(clippy::too_many_arguments)]
pub fn move_rad_to_atoms_add_charges(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &mut [ValAt],
    tcg: &AllTcGroups,
    forbidden_mask: i32,
) -> i32 {
    let num_atoms = bns.num_atoms as usize;
    let len_at = (ps.num_atoms + ps.num_deleted_h) as usize;
    let mut num_rad_not_atom = 0;
    for i in num_atoms..bns.num_vertices as usize {
        num_rad_not_atom += bns.vert[i].st_edge.cap - bns.vert[i].st_edge.flow;
    }
    if num_rad_not_atom == 0 {
        return 0;
    }
    let nv = bns.num_vertices as usize;
    let mut pn_rad: Vec<i16> = (0..nv)
        .map(|i| (bns.vert[i].st_edge.cap - bns.vert[i].st_edge.flow) as i16)
        .collect();
    let mut pn_delta = vec![0i16; num_atoms];
    let mut num_moved = 0;
    let mut num_candidates = 0;
    let mut extra_charge = 0;
    let mut atom_rad_removed = 0;
    loop {
        for i in 0..num_atoms {
            pn_delta[i] = (bns.vert[i].st_edge.cap - bns.vert[i].st_edge.flow) as i16;
            bns.vert[i].st_edge.cap -= i32::from(pn_delta[i]);
            atom_rad_removed += i32::from(pn_delta[i] != 0);
        }
        let ret = set_rad_endpoints(bns, bd, RAD_SRCH_FROM_FICT);
        if ret == 0 {
            break;
        }
        if ret < 0 {
            return ret;
        }
        for i in 0..ret as usize {
            let ie = bd.rad_edges[i] as usize;
            let (n1, n12, fl) = (
                bns.edge[ie].neighbor1 as usize,
                bns.edge[ie].neighbor12 as usize,
                bns.edge[ie].flow,
            );
            let v2 = n12 ^ n1;
            bns.vert[n1].st_edge.flow -= fl;
            bns.vert[v2].st_edge.flow -= fl;
            bns.tot_st_flow -= 2 * fl;
            bns.edge[ie].flow = 0;
            bns.edge[ie].forbidden |= forbidden_mask as i8;
            bns.edge_forbidden_mask |= forbidden_mask as i8;
        }
        let ret = run_bns_restore_once(bns, bd);
        if ret < 0 {
            return ret;
        }
        num_moved += ret;
        remove_rad_endpoints(bns, bd, None);
        if ret == 0 {
            break;
        }
        if atom_rad_removed != 0 {
            for i in 0..num_atoms {
                bns.vert[i].st_edge.cap += i32::from(pn_delta[i]);
            }
            atom_rad_removed = 0;
        }
    }
    if atom_rad_removed != 0 {
        for i in 0..num_atoms {
            bns.vert[i].st_edge.cap += i32::from(pn_delta[i]);
        }
    }
    bns.edge_forbidden_mask &= !(forbidden_mask as i8);

    if num_moved != 0 {
        extra_charge = 0;
        for i in num_atoms..nv {
            let pv = &bns.vert[i];
            let delta = pv.st_edge.cap - pv.st_edge.flow;
            if delta != 0 {
                if pv.type_ & BNS_VT_C_POS != 0 {
                    extra_charge -= delta;
                } else if BNS_VERT_TYPE_AUX == pv.type_ {
                    extra_charge += delta;
                } else {
                    return RI_ERR_PROGR;
                }
            }
        }
        if extra_charge == 0 {
            return 0;
        }
        num_candidates = 0;
        for i in 0..nv {
            pn_rad[i] = (bns.vert[i].st_edge.cap - bns.vert[i].st_edge.flow) as i16;
            if pn_rad[i] > 0 && i < num_atoms && p_va[i].n_taut_group_edge == 0 {
                num_candidates += 1;
            }
        }
    }
    if num_candidates > 0 {
        at2.clear();
        at2.extend_from_slice(&ps.at[..len_at]);
        std::mem::swap(&mut ps.at, at2);
        let ret2 = copy_bns_to_atom(ps, bns, p_va, tcg, true);
        std::mem::swap(&mut ps.at, at2);
        if ret2 < 0 {
            return ret2;
        }
        let mut cand: Vec<CcCand> = Vec::new();
        for i in 0..nv {
            if pn_rad[i] > 0 && i < num_atoms && p_va[i].n_taut_group_edge == 0 {
                let mut c = CcCand {
                    iat: i as i16,
                    num_bonds: at2[i].valence,
                    chem_valence: at2[i].chem_bonds_valence,
                    c_metal: p_va[i].c_metal,
                    c_num_bonds_to_metal: p_va[i].c_num_bonds_to_metal,
                    c_num_valence_electrons: p_va[i].c_num_valence_electrons,
                    c_periodic_row_number: p_va[i].c_periodic_row_number,
                    c_num_charge_states: 0,
                    el_number: at2[i].el_number,
                };
                let mut cn_bits = if p_va[i].c_cn_list_index > 0 {
                    CN_LIST[p_va[i].c_cn_list_index as usize - 1].bits
                } else {
                    0
                };
                while cn_bits > 0 {
                    c.c_num_charge_states += 1;
                    cn_bits >>= CN_BITS_SHIFT;
                }
                cand.push(c);
            }
        }
        let j = cand.len();
        if j > 1 {
            insertions_sort(&mut cand, comp_cc_cand);
        }
        let mut added_charge = 0;
        for k in 0..j {
            let rest = extra_charge - added_charge;
            let per_left = (rest.abs() + j as i32 - k as i32 - 1) / (j - k) as i32;
            let mut this_add = if rest > 0 { per_left } else { -per_left };
            p_va[cand[k].iat as usize].c_init_charge += this_add as i8;
            added_charge += this_add;
            if this_add != 0 {
                let mut i = nv as isize - 1;
                while this_add != 0 && num_atoms as isize <= i {
                    let pv = &mut bns.vert[i as usize];
                    let mut delta = pv.st_edge.cap - pv.st_edge.flow;
                    if delta != 0 {
                        if this_add < 0 && pv.type_ & BNS_VT_C_POS != 0 {
                            if delta + this_add > 0 {
                                delta = -this_add;
                            }
                            pv.st_edge.cap -= delta;
                            bns.tot_st_cap -= delta;
                            this_add += delta;
                        } else if this_add > 0 && BNS_VERT_TYPE_AUX == pv.type_ {
                            if delta > this_add {
                                delta = this_add;
                            }
                            pv.st_edge.cap -= delta;
                            bns.tot_st_cap -= delta;
                            this_add -= delta;
                        }
                    }
                    i -= 1;
                }
            }
        }
    }
    0
}

#[derive(Clone, Copy, Default)]
#[allow(dead_code)]
struct MobileGr {
    group_number: u16,
    atom_number: usize,
    atom_type_pva: i32,
    ineigh: usize,
    bond_type: i32,
    forbidden: i32,
    endpoint_valence: i32,
    num_bonds: i32,
    bonds_valence: i32,
    num_bonds_non_metal: i32,
    bonds_valence_non_metal: i32,
}

#[derive(Clone, Copy, Default)]
struct MGroups {
    group_number: u16,
    num: i32,
}

/// the edge to the t-group of an endpoint (`nTautGroupEdge - 1`); the C code
/// reads before the edge array when the atom has none, here that case is skipped
fn taut_edge(p_va: &[ValAt], atom: usize) -> Option<usize> {
    let e = p_va[atom].n_taut_group_edge;
    (e > 0).then(|| e as usize - 1)
}

/// `AdjustTgroupsToForbiddenEdges2`: the cases of the InChI Technical
/// Manual, Table 5; returns the number of changes
pub fn adjust_tgroups_to_forbidden_edges2(
    bns: &mut BnStruct,
    at: &[InpAtom],
    p_va: &[ValAt],
    num_atoms: usize,
    forbidden_mask: i32,
) -> i32 {
    let inv_mask = !forbidden_mask;
    let mut num_changes = 0;
    // the edge `j` of the vertex `v`
    fn ei(b: &BnStruct, v: usize, j: usize) -> usize {
        b.iedge[b.vert[v].iedge + j] as usize
    }
    fn other(b: &BnStruct, e: usize, v: usize) -> usize {
        b.edge[e].neighbor12 as usize ^ v
    }
    'atoms: for i in 0..num_atoms {
        let a = &at[i];
        let val0 = get_el_valence(i32::from(a.el_number), 0, 0);
        let center_type = get_pva_atom_type(p_va, at, i, 0);
        if a.chem_bonds_valence == a.valence
            || a.num_h != 0
            || a.endpoint != 0
            || a.charge != 0
            || a.radical != 0
            || !is_centerpoint_elem(a.el_number)
            || center_type == 0
        {
            continue;
        }
        let delta_val = i32::from(a.chem_bonds_valence) - val0;
        if 2 > delta_val {
            continue;
        }
        let bonds_valence_non_metal = n_no_metal_bonds_valence(at, i);
        let delta_met = bonds_valence_non_metal - val0;
        if 2 > delta_met {
            continue;
        }
        let (
            mut num_donors,
            mut num_acceptors,
            mut num_donor_endpoints,
            mut num_acceptor_endpoints,
        ) = (0, 0, 0, 0);
        let (mut num_mgroups, mut num_diff_t_groups) = (0usize, 0);
        let mut mgr: Vec<MobileGr> = Vec::new();
        let mut mgs = [MGroups::default(); MAXVAL];
        let (mut has_mobile_h, mut num_eql, mut num_dif, mut tg_number) = (0, 0, 0, 0u16);
        let mut ind_forbidden: isize = -1;
        let mut num_forbidden = 0;
        let (mut num_n, mut num_o, mut num_p, mut num_s, mut num_ost) = (0, 0, 0, 0, 0);
        let num_bonds_non_metal = n_no_metal_num_bonds(at, i);

        for j in 0..a.valence as usize {
            let neigh = a.neighbor[j] as usize;
            let val = get_endpoint_valence(at[neigh].el_number);
            let forbidden = i32::from(bns.edge[ei(bns, i, j)].forbidden);
            let bond_type = i32::from(a.bond_type[j]) & BOND_TYPE_MASK;
            let neigh_type = get_pva_atom_type(p_va, at, neigh, bond_type);
            if forbidden == 0 && at[neigh].endpoint == 0 {
                if is_el_a_metal(i32::from(at[neigh].el_number)) {
                    continue;
                }
                match bond_type {
                    BOND_TYPE_SINGLE => {
                        if at[neigh].num_h == 0 && at[neigh].charge != -1 {
                            continue;
                        }
                    }
                    BOND_TYPE_DOUBLE => {
                        if neigh_type == 0 {
                            continue;
                        }
                    }
                    _ => continue,
                }
            }
            mgr.push(MobileGr {
                atom_number: neigh,
                ineigh: j,
                bond_type,
                group_number: at[neigh].endpoint,
                endpoint_valence: val,
                forbidden,
                atom_type_pva: neigh_type,
                num_bonds: i32::from(at[neigh].valence),
                bonds_valence: i32::from(at[neigh].chem_bonds_valence),
                num_bonds_non_metal: n_no_metal_num_bonds(at, neigh),
                bonds_valence_non_metal: n_no_metal_bonds_valence(at, neigh),
            });
            if forbidden & forbidden_mask != 0 {
                num_forbidden += 1;
                ind_forbidden = mgr.len() as isize - 1;
            }
            let nb = &at[neigh];
            num_o += i32::from(neigh_type & EL_TYPE_O != 0 && nb.valence == 1);
            num_n += i32::from(
                neigh_type & EL_TYPE_N != 0 && !(nb.valence == 3 && nb.chem_bonds_valence == 3),
            );
            num_s += i32::from(neigh_type & EL_TYPE_S != 0 && nb.valence == 1);
            num_p += i32::from(
                neigh_type & EL_TYPE_P != 0 && !(nb.valence == 3 && nb.chem_bonds_valence == 3),
            );
            num_ost += i32::from(neigh_type & EL_TYPE_OST != 0);
            num_acceptors +=
                i32::from(bond_type == BOND_TYPE_DOUBLE && neigh_type & EL_TYPE_PT != 0);
            num_donors += i32::from(
                bond_type == BOND_TYPE_SINGLE
                    && neigh_type & EL_TYPE_PT != 0
                    && (nb.num_h != 0 || nb.charge == -1 || nb.endpoint != 0),
            );
            if nb.endpoint != 0 {
                num_acceptor_endpoints += i32::from(bond_type == BOND_TYPE_DOUBLE);
                num_donor_endpoints += i32::from(bond_type == BOND_TYPE_SINGLE);
                if tg_number == 0 {
                    tg_number = nb.endpoint;
                    num_eql = 1;
                } else if tg_number == nb.endpoint {
                    num_eql += 1;
                } else {
                    num_dif += 1;
                }
            } else if bond_type == BOND_TYPE_SINGLE && val != 0 {
                has_mobile_h |= i32::from(nb.num_h != 0);
            }
            if nb.endpoint != 0 || neigh_type & EL_TYPE_PT != 0 {
                let mut k = 0;
                while k < num_mgroups {
                    if mgs[k].group_number == nb.endpoint {
                        mgs[k].num += 1;
                        break;
                    }
                    k += 1;
                }
                if k == num_mgroups {
                    mgs[k].group_number = nb.endpoint;
                    mgs[k].num = 1;
                    num_mgroups += 1;
                    num_diff_t_groups += i32::from(nb.endpoint != 0);
                }
            }
        }
        let num_endpoints = mgr.len() as i32;
        let _ = (num_p, num_donors);
        if num_acceptors == 0
            || num_donors == 0
            || (num_eql == num_endpoints && num_forbidden == 0)
            || (tg_number == 0 && has_mobile_h == 0)
        {
            continue;
        }
        let ne = mgr.len();

        // case 5.1
        'c51: {
            if 3 == num_bonds_non_metal
                && 4 == bonds_valence_non_metal
                && center_type == EL_TYPE_C
                && 2 == num_o
                && 1 == num_n + num_s
                && num_ost != 0
                && 1 == num_forbidden
                && 3 == num_eql
            {
                let (mut jxh, mut jo1, mut jo2, mut n) = (-1isize, -1isize, -1isize, 0);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & (EL_TYPE_N | EL_TYPE_S) != 0
                        && m.forbidden == forbidden_mask
                        && m.bond_type == BOND_TYPE_SINGLE
                        && jxh < 0
                    {
                        jxh = j as isize;
                        n += 1;
                    } else if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_O
                        && m.num_bonds_non_metal == 1
                        && m.forbidden == 0
                    {
                        if m.bond_type == BOND_TYPE_DOUBLE && jo1 < 0 {
                            jo1 = j as isize;
                            n += 1;
                        } else if m.bond_type == BOND_TYPE_SINGLE && jo2 < 0 {
                            jo2 = j as isize;
                            n += 1;
                        }
                    }
                }
                if n != 3 {
                    break 'c51;
                }
                let e = ei(bns, i, mgr[jxh as usize].ineigh);
                let ev1 = ei(bns, i, mgr[jo1 as usize].ineigh);
                let Some(ev2) = taut_edge(p_va, mgr[jxh as usize].atom_number) else {
                    break 'c51;
                };
                if bns.edge[ev1].flow == 0 || bns.edge[ev2].flow == 0 {
                    break 'c51;
                }
                bns.edge[e].flow += 1;
                bns.edge[ev1].flow -= 1;
                bns.edge[ev2].flow -= 1;
                let v = other(bns, ev1, i);
                bns.vert[v].st_edge.flow -= 1;
                let v = other(bns, ev2, bns.edge[ev2].neighbor1 as usize);
                bns.vert[v].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.2
        'c52: {
            if 3 == num_bonds_non_metal
                && 5 == bonds_valence_non_metal
                && center_type == EL_TYPE_N
                && 2 == num_o
                && 1 == num_n + num_s
                && 1 == num_forbidden
                && 3 == num_eql
            {
                let (mut jxh, mut jo1, mut jo2, mut n) = (-1isize, -1isize, -1isize, 0);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & (EL_TYPE_N | EL_TYPE_S) != 0
                        && m.forbidden == forbidden_mask
                        && m.bond_type == BOND_TYPE_SINGLE
                        && jxh < 0
                    {
                        jxh = j as isize;
                        n += 1;
                    } else if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_O
                        && m.bond_type == BOND_TYPE_DOUBLE
                        && m.num_bonds_non_metal == 1
                        && m.forbidden == 0
                    {
                        if jo1 < 0 {
                            jo1 = j as isize;
                            n += 1;
                        } else if jo2 < 0 {
                            jo2 = j as isize;
                            n += 1;
                        }
                    }
                }
                if n != 3 {
                    break 'c52;
                }
                let e = ei(bns, i, mgr[jxh as usize].ineigh);
                let ev1 = ei(bns, i, mgr[jo1 as usize].ineigh);
                let Some(ev2) = taut_edge(p_va, mgr[jxh as usize].atom_number) else {
                    break 'c52;
                };
                if bns.edge[ev1].flow == 0 || bns.edge[ev2].flow == 0 {
                    break 'c52;
                }
                bns.edge[e].flow += 1;
                bns.edge[ev1].flow -= 1;
                bns.edge[ev2].flow -= 1;
                let v = other(bns, ev1, i);
                bns.vert[v].st_edge.flow -= 1;
                let v = other(bns, ev2, bns.edge[ev2].neighbor1 as usize);
                bns.vert[v].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.4
        'c54: {
            if 3 == num_bonds_non_metal
                && 5 == bonds_valence_non_metal
                && center_type & (EL_TYPE_N | EL_TYPE_P) != 0
                && 1 == num_o + num_s
                && 0 < num_n
                && 2 == num_n + num_p
                && 1 == num_forbidden
                && num_o + num_s + num_n == num_eql
            {
                let (mut jo, mut jnh, mut jx, mut n) = (-1isize, -1isize, -1isize, 0);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & (EL_TYPE_O | EL_TYPE_S) != 0
                        && m.forbidden == forbidden_mask
                        && m.bond_type == BOND_TYPE_DOUBLE
                        && m.num_bonds_non_metal == 1
                        && jo < 0
                    {
                        jo = j as isize;
                        n += 1;
                    } else if m.atom_type_pva & (EL_TYPE_N | EL_TYPE_P) != 0 && m.forbidden == 0 {
                        if m.bond_type == BOND_TYPE_SINGLE
                            && m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N
                            && jnh < 0
                        {
                            jnh = j as isize;
                            n += 1;
                        } else if m.bond_type == BOND_TYPE_DOUBLE && jx < 0 {
                            jx = j as isize;
                            n += 1;
                        }
                    }
                }
                if n != 3 {
                    break 'c54;
                }
                let e = ei(bns, i, mgr[jo as usize].ineigh);
                let ev1 = ei(bns, i, mgr[jnh as usize].ineigh);
                if bns.edge[e].flow == 0 {
                    break 'c54;
                }
                bns.edge[e].flow -= 1;
                bns.edge[ev1].flow += 1;
                let v = other(bns, e, i);
                bns.vert[v].st_edge.flow -= 1;
                let v = other(bns, ev1, i);
                bns.vert[v].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.6
        'c56: {
            if 2 == delta_met
                && 4 == num_bonds_non_metal
                && 5 == bonds_valence_non_metal
                && 1 == num_forbidden
                && 1 < num_eql
                && num_dif == 0
                && center_type & (EL_TYPE_N | EL_TYPE_P) != 0
                && 1 <= num_n
                && 2 <= num_n + num_o + num_s
                && 1 == num_acceptor_endpoints
                && 0 < num_donor_endpoints
            {
                let mut jn: isize = -1;
                let mut jfix: Vec<usize> = Vec::new();
                let mut n = 0;
                for j in 0..ne {
                    let m = &mgr[j];
                    let neigh = m.atom_number;
                    if m.bond_type == BOND_TYPE_DOUBLE
                        && m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N
                        && m.num_bonds_non_metal == 2
                        && m.bonds_valence_non_metal == 3
                        && at[neigh].endpoint != 0
                        && at[neigh].num_h == 0
                        && at[neigh].charge == 0
                        && at[neigh].radical == 0
                        && m.forbidden & forbidden_mask != 0
                        && jn < 0
                    {
                        jn = j as isize;
                        n += 1;
                    } else if m.bond_type == BOND_TYPE_SINGLE && at[neigh].endpoint != 0 {
                        if m.num_bonds > 1 {
                            jfix.push(j);
                        }
                        n += 1;
                    }
                }
                if jn < 0 || n < 2 || 1 + jfix.len() as i32 == n {
                    break 'c56;
                }
                let e = ei(bns, i, mgr[jn as usize].ineigh);
                if bns.edge[e].flow == 0 {
                    break 'c56;
                }
                bns.edge[e].flow -= 1;
                bns.vert[i].st_edge.flow -= 1;
                let v = other(bns, e, i);
                bns.vert[v].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                for &jf in &jfix {
                    let ev = ei(bns, i, mgr[jf].ineigh);
                    bns.edge[ev].forbidden |= forbidden_mask as i8;
                }
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.7
        'c57: {
            if 3 == num_bonds_non_metal
                && 4 == bonds_valence_non_metal
                && center_type == EL_TYPE_S
                && 2 == num_o + num_s
                && 1 == num_ost
                && 1 == num_n
                && 1 == num_forbidden
                && 3 == num_eql
                && mgr[ind_forbidden as usize].bond_type == BOND_TYPE_SINGLE
            {
                let (mut jo, mut jnh, mut joh, mut n) = (-1isize, -1isize, -1isize, 0);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N
                        && m.forbidden == forbidden_mask
                        && m.bond_type == BOND_TYPE_SINGLE
                        && m.num_bonds_non_metal <= 2
                        && jnh < 0
                    {
                        jnh = j as isize;
                        n += 1;
                    } else if m.atom_type_pva & (EL_TYPE_O | EL_TYPE_S) != 0
                        && m.forbidden == 0
                        && m.num_bonds_non_metal == 1
                    {
                        if m.bond_type == BOND_TYPE_DOUBLE && jo < 0 {
                            jo = j as isize;
                            n += 1;
                        } else if joh < 0 {
                            joh = j as isize;
                            n += 1;
                        }
                    }
                }
                if n != 3 {
                    break 'c57;
                }
                let e = ei(bns, i, mgr[jnh as usize].ineigh);
                let ev1 = ei(bns, i, mgr[jo as usize].ineigh);
                let Some(ev2) = taut_edge(p_va, mgr[jnh as usize].atom_number) else {
                    break 'c57;
                };
                if bns.edge[ev1].flow == 0 || bns.edge[ev2].flow == 0 {
                    break 'c57;
                }
                bns.edge[e].flow += 1;
                bns.edge[ev1].flow -= 1;
                bns.edge[ev2].flow -= 1;
                let v = other(bns, ev1, i);
                bns.vert[v].st_edge.flow -= 1;
                let v = other(bns, ev2, bns.edge[ev2].neighbor1 as usize);
                bns.vert[v].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.9a
        'c59a: {
            if 3 == num_bonds_non_metal
                && 4 == bonds_valence_non_metal
                && center_type == EL_TYPE_S
                && 1 == num_o + num_s
                && num_ost == 0
                && 1 <= num_n
                && 1 == num_forbidden
                && num_o + num_s + num_n == num_eql
                && mgr[ind_forbidden as usize].bond_type == BOND_TYPE_SINGLE
            {
                let (mut jo, mut jnh) = (-1isize, -1isize);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & (EL_TYPE_O | EL_TYPE_S) != 0
                        && m.forbidden == forbidden_mask
                        && m.bond_type == BOND_TYPE_DOUBLE
                        && m.num_bonds_non_metal == 1
                        && jo < 0
                    {
                        jo = j as isize;
                    } else if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N
                        && m.forbidden == 0
                        && m.bond_type == BOND_TYPE_SINGLE
                        && jnh < 0
                    {
                        jnh = j as isize;
                    }
                }
                if jo < 0 || jnh < 0 {
                    break 'c59a;
                }
                let e = ei(bns, i, mgr[ind_forbidden as usize].ineigh);
                if bns.edge[e].flow == 0 {
                    break 'c59a;
                }
                bns.edge[e].flow -= 1;
                let n1 = bns.edge[e].neighbor1 as usize;
                let n12 = bns.edge[e].neighbor12 as usize;
                bns.vert[n1].st_edge.flow -= 1;
                bns.vert[n1 ^ n12].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.8b -> 9b
        'c58b: {
            if 3 == num_bonds_non_metal
                && 4 == bonds_valence_non_metal
                && center_type == EL_TYPE_S
                && 0 == num_o + num_s
                && 2 == num_n
                && 0 == num_p
                && num_ost == 0
                && 1 == num_forbidden
                && 1 == num_eql
                && 0 == num_dif
                && 1 == num_donor_endpoints
                && 0 == num_acceptor_endpoints
                && 1 == num_donors
                && 1 == num_acceptors
                && mgr[ind_forbidden as usize].bond_type == BOND_TYPE_SINGLE
            {
                let (mut jn, mut jnh, mut jx, mut n) = (-1isize, -1isize, -1isize, 0);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N
                        && m.forbidden & forbidden_mask == 0
                    {
                        if m.bond_type == BOND_TYPE_DOUBLE
                            && at[m.atom_number].endpoint == 0
                            && jn < 0
                        {
                            jn = j as isize;
                            n += 1;
                        } else if m.bond_type == BOND_TYPE_SINGLE
                            && m.num_bonds == 2
                            && m.bonds_valence == 2
                            && at[m.atom_number].endpoint != 0
                            && jnh < 0
                        {
                            jnh = j as isize;
                            n += 1;
                        }
                    } else if !(m.atom_type_pva & (EL_TYPE_N | EL_TYPE_P) != 0
                        || (m.atom_type_pva & (EL_TYPE_O | EL_TYPE_S) != 0 && m.num_bonds > 1))
                        && m.forbidden & forbidden_mask != 0
                        && m.bond_type == BOND_TYPE_SINGLE
                        && jx < 0
                    {
                        jx = j as isize;
                        n += 1;
                    }
                }
                if n != 3 {
                    break 'c58b;
                }
                let Some(e) = taut_edge(p_va, mgr[jnh as usize].atom_number) else {
                    break 'c58b;
                };
                if bns.edge[e].flow == 0 {
                    break 'c58b;
                }
                bns.edge[e].flow -= 1;
                let n1 = bns.edge[e].neighbor1 as usize;
                let n12 = bns.edge[e].neighbor12 as usize;
                bns.vert[n1].st_edge.flow -= 1;
                bns.vert[n1 ^ n12].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                bns.edge[e].forbidden |= forbidden_mask as i8;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.8c -> 9c
        'c58c: {
            if 3 == num_bonds_non_metal
                && 4 == bonds_valence_non_metal
                && center_type == EL_TYPE_S
                && 0 == num_o + num_s
                && 3 == num_n
                && 0 == num_p
                && 1 == num_forbidden
                && 3 == num_eql
                && 0 == num_dif
                && 2 == num_donor_endpoints
                && 1 == num_acceptor_endpoints
                && 2 == num_donors
                && 1 == num_acceptors
                && mgr[ind_forbidden as usize].bond_type == BOND_TYPE_SINGLE
            {
                let (mut jn1, mut jn2, mut jx, mut n) = (-1isize, -1isize, -1isize, 0);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N
                        && m.forbidden & forbidden_mask == 0
                    {
                        if m.bond_type == BOND_TYPE_DOUBLE
                            && at[m.atom_number].endpoint != 0
                            && jn1 < 0
                        {
                            jn1 = j as isize;
                            n += 1;
                        } else if m.bond_type == BOND_TYPE_SINGLE
                            && m.num_bonds == 2
                            && m.bonds_valence == 3
                            && m.forbidden == forbidden_mask
                            && at[m.atom_number].endpoint != 0
                            && jn2 < 0
                        {
                            jn2 = j as isize;
                            n += 1;
                        } else if m.bond_type == BOND_TYPE_SINGLE
                            && m.num_bonds <= 2
                            && m.bonds_valence <= 3
                            && at[m.atom_number].endpoint != 0
                            && jx < 0
                        {
                            jx = j as isize;
                            n += 1;
                        }
                    }
                }
                if n != 3 {
                    break 'c58c;
                }
                let Some(e) = taut_edge(p_va, mgr[jn2 as usize].atom_number) else {
                    break 'c58c;
                };
                if bns.edge[e].flow != 0 {
                    break 'c58c;
                }
                let mut pv1 = bns.edge[e].neighbor1 as usize;
                let mut pv2 = pv1 ^ bns.edge[e].neighbor12 as usize;
                let (mut ie1, mut ie2) = (-1isize, -1isize);
                let (mut ev1, mut ev2) = (0usize, 0usize);
                for j in 0..bns.vert[pv1].num_adj_edges as usize {
                    ev1 = ei(bns, pv1, j);
                    if bns.edge[ev1].flow != 0 && bns.edge[ev1].forbidden == 0 {
                        ie1 = ev1 as isize;
                        pv1 = other(bns, ev1, pv1);
                        break;
                    }
                }
                for j in 0..bns.vert[pv2].num_adj_edges as usize {
                    ev2 = ei(bns, pv2, j);
                    if bns.edge[ev2].flow != 0 && bns.edge[ev2].forbidden == 0 {
                        ie2 = ev2 as isize;
                        pv2 = other(bns, ev2, pv2);
                        break;
                    }
                }
                if ie1 < 0 || ie2 < 0 {
                    break 'c58c;
                }
                bns.edge[e].flow += 1;
                bns.edge[e].forbidden |= forbidden_mask as i8;
                bns.edge[ev1].flow -= 1;
                bns.edge[ev2].flow -= 1;
                bns.vert[pv1].st_edge.flow -= 1;
                bns.vert[pv2].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.9b -> 8b
        'c59b: {
            if 3 == num_bonds_non_metal
                && 4 == bonds_valence_non_metal
                && center_type == EL_TYPE_S
                && 0 == num_o + num_s
                && 2 == num_n
                && 0 == num_p
                && 1 == num_forbidden
                && 2 == num_eql
                && 0 == num_dif
                && 1 == num_donor_endpoints
                && 1 == num_acceptor_endpoints
                && 1 == num_donors
                && 1 == num_acceptors
                && mgr[ind_forbidden as usize].bond_type == BOND_TYPE_DOUBLE
            {
                let (mut jn1, mut jn2) = (-1isize, -1isize);
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N {
                        if m.bond_type == BOND_TYPE_DOUBLE
                            && at[m.atom_number].endpoint != 0
                            && m.forbidden == forbidden_mask
                            && jn1 < 0
                        {
                            jn1 = j as isize;
                        } else if m.bond_type == BOND_TYPE_SINGLE
                            && m.num_bonds == 2
                            && m.bonds_valence == 3
                            && at[m.atom_number].endpoint != 0
                            && m.forbidden & forbidden_mask == 0
                            && jn2 < 0
                        {
                            jn2 = j as isize;
                        }
                    }
                }
                if jn1 < 0 || jn2 < 0 {
                    break 'c59b;
                }
                let ev = ei(bns, i, mgr[jn1 as usize].ineigh);
                let Some(e) = taut_edge(p_va, mgr[jn2 as usize].atom_number) else {
                    break 'c59b;
                };
                if bns.edge[e].flow != 0 {
                    break 'c59b;
                }
                let mut pv1 = bns.edge[e].neighbor1 as usize;
                let mut pv2 = pv1 ^ bns.edge[e].neighbor12 as usize;
                let (mut ie1, mut ie2) = (-1isize, -1isize);
                let (mut ev1, mut ev2) = (0usize, 0usize);
                bns.edge[ev].forbidden &= inv_mask as i8;
                for j in 0..bns.vert[pv1].num_adj_edges as usize {
                    ev1 = ei(bns, pv1, j);
                    if bns.edge[ev1].flow != 0 && bns.edge[ev1].forbidden == 0 {
                        ie1 = ev1 as isize;
                        pv1 = other(bns, ev1, pv1);
                        break;
                    }
                }
                for j in 0..bns.vert[pv2].num_adj_edges as usize {
                    ev2 = ei(bns, pv2, j);
                    if bns.edge[ev2].flow != 0 && bns.edge[ev2].forbidden == 0 {
                        ie2 = ev2 as isize;
                        pv2 = other(bns, ev2, pv2);
                        break;
                    }
                }
                if ie1 < 0 || ie2 < 0 {
                    bns.edge[ev].forbidden |= forbidden_mask as i8;
                    break 'c59b;
                }
                bns.edge[e].flow += 1;
                bns.edge[e].forbidden |= forbidden_mask as i8;
                bns.edge[ev1].flow -= 1;
                bns.edge[ev2].flow -= 1;
                bns.vert[pv1].st_edge.flow -= 1;
                bns.vert[pv2].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.9c -> 8c
        'c59c: {
            if 3 == num_bonds_non_metal
                && 4 == bonds_valence_non_metal
                && center_type == EL_TYPE_S
                && 0 == num_o + num_s
                && 3 == num_n
                && 0 == num_p
                && 0 == num_forbidden
                && 2 == num_diff_t_groups
                && 2 == num_mgroups
                && 3 == num_eql + num_dif
                && 2 == num_donor_endpoints
                && 1 == num_acceptor_endpoints
                && 2 == num_donors
                && 1 == num_acceptors
            {
                let (mut jn1, mut jn2, mut jx, mut n) = (-1isize, -1isize, -1isize, 0);
                // (sic) the C loop tests MGroups[k] with k, not j
                let mut k: u16 = 0;
                for _ in 0..num_mgroups {
                    if 1 == mgs[k as usize].num && mgs[k as usize].group_number != 0 {
                        k = mgs[k as usize].group_number;
                        break;
                    }
                }
                if k == 0 {
                    break 'c59c;
                }
                for j in 0..ne {
                    let m = &mgr[j];
                    if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N {
                        let ep = at[m.atom_number].endpoint;
                        if m.bond_type == BOND_TYPE_DOUBLE && ep != 0 && ep != k && jn1 < 0 {
                            jn1 = j as isize;
                            n += 1;
                        } else if m.bond_type == BOND_TYPE_SINGLE
                            && m.num_bonds == 2
                            && m.bonds_valence == 2
                            && ep == k
                            && jn2 < 0
                        {
                            jn2 = j as isize;
                            n += 1;
                        } else if m.bond_type == BOND_TYPE_SINGLE
                            && m.num_bonds <= 2
                            && m.bonds_valence <= 3
                            && ep != 0
                            && ep != k
                            && jx < 0
                        {
                            jx = j as isize;
                            n += 1;
                        }
                    }
                }
                if n != 3 {
                    break 'c59c;
                }
                let Some(e) = taut_edge(p_va, mgr[jn2 as usize].atom_number) else {
                    break 'c59c;
                };
                if bns.edge[e].flow == 0 {
                    break 'c59c;
                }
                bns.edge[e].flow -= 1;
                let n1 = bns.edge[e].neighbor1 as usize;
                let n12 = bns.edge[e].neighbor12 as usize;
                bns.vert[n1].st_edge.flow -= 1;
                bns.vert[n1 ^ n12].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                bns.edge[e].forbidden |= forbidden_mask as i8;
                num_changes += 1;
                continue 'atoms;
            }
        }
        // case 5.9d
        if 3 == num_bonds_non_metal
            && 4 == bonds_valence_non_metal
            && center_type == EL_TYPE_S
            && 0 == num_o + num_s
            && 3 == num_n
            && 0 == num_p
            && 0 == num_forbidden
            && 3 == num_mgroups
            && 2 == num_diff_t_groups
            && 2 == num_donor_endpoints
            && 0 == num_acceptor_endpoints
            && 2 == num_donors
            && 1 == num_acceptors
        {
            let (mut jn1, mut jn2, mut jx, mut n) = (-1isize, -1isize, -1isize, 0);
            for j in 0..ne {
                let m = &mgr[j];
                if m.atom_type_pva & EL_TYPE_MASK == EL_TYPE_N {
                    if m.bond_type == BOND_TYPE_DOUBLE
                        && at[m.atom_number].endpoint == 0
                        && m.forbidden & forbidden_mask == 0
                        && jn1 < 0
                    {
                        jn1 = j as isize;
                        n += 1;
                    } else if m.bond_type == BOND_TYPE_SINGLE
                        && m.num_bonds == 2
                        && m.bonds_valence <= 3
                        && at[m.atom_number].endpoint != 0
                        && m.forbidden & forbidden_mask == 0
                    {
                        if jn2 < 0 {
                            jn2 = j as isize;
                            n += 1;
                        } else if jx < 0 {
                            jx = j as isize;
                            n += 1;
                        }
                    }
                }
            }
            if n != 3 {
                continue;
            }
            let ev = ei(bns, i, mgr[jn1 as usize].ineigh);
            if bns.edge[ev].flow == 0 {
                continue;
            }
            let Some(e1) = taut_edge(p_va, mgr[jn2 as usize].atom_number) else {
                continue;
            };
            if bns.edge[e1].flow == 0 {
                continue;
            }
            let Some(e2) = taut_edge(p_va, mgr[jx as usize].atom_number) else {
                continue;
            };
            if bns.edge[e2].flow == 0 {
                continue;
            }
            for e in [e1, e2] {
                bns.edge[e].flow -= 1;
                let n1 = bns.edge[e].neighbor1 as usize;
                let n12 = bns.edge[e].neighbor12 as usize;
                bns.vert[n1].st_edge.flow -= 1;
                bns.vert[n1 ^ n12].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                bns.edge[e].forbidden |= forbidden_mask as i8;
                num_changes += 1;
            }
            bns.edge[ev].forbidden |= forbidden_mask as i8;
        }
    }
    num_changes
}

/// `RearrangePlusMinusEdgesFlow`: replaces the ambiguous neutral
/// (+)flow=0, (-)flow=1 by (+)flow=1, (-)flow=0 and fixes the charges
pub fn rearrange_plus_minus_edges_flow(
    bns: &mut BnStruct,
    bd: &mut BnData,
    p_va: &[ValAt],
    forbidden_edge_mask: i32,
) -> i32 {
    let num_atoms = bns.num_atoms as usize;
    let mask = forbidden_edge_mask as i8;
    let mut num_found = 0;
    let mut num_tot = 0;
    for va in &p_va[..num_atoms] {
        let e_minus = va.n_c_minus_group_edge - 1;
        let e_plus = va.n_c_plus_group_edge - 1;
        num_tot += i32::from(e_minus >= 0) + i32::from(e_plus >= 0);
        if e_minus >= 0 && e_plus >= 0 {
            let p = &bns.edge[e_plus as usize];
            let m = &bns.edge[e_minus as usize];
            if m.flow > 0 && p.cap - p.flow > 0 {
                num_found += 1;
            }
        }
    }
    if num_found == 0 {
        return 0;
    }
    let mut fixed = EdgeList::default();
    let ret = alloc_edge_list(&mut fixed, num_tot + bns.num_bonds);
    if ret != 0 {
        return ret;
    }
    for va in &p_va[..num_atoms] {
        let e_minus = va.n_c_minus_group_edge - 1;
        let e_plus = va.n_c_plus_group_edge - 1;
        if e_minus >= 0 && e_plus >= 0 {
            let (k1, k2) = (
                bns.edge[e_minus as usize].flow,
                bns.edge[e_plus as usize].cap - bns.edge[e_plus as usize].flow,
            );
            if k1 > 0 && k2 > 0 {
                let v1 = bns.edge[e_minus as usize].neighbor1 as usize;
                let v2 = bns.edge[e_minus as usize].neighbor12 as usize ^ v1;
                let delta = k1.min(k2);
                bns.edge[e_minus as usize].flow -= delta;
                bns.vert[v1].st_edge.flow -= delta;
                bns.vert[v2].st_edge.flow -= delta;
                bns.tot_st_flow -= 2 * delta;
            }
            bns.edge[e_plus as usize].forbidden |= mask;
            bns.edge[e_minus as usize].forbidden |= mask;
            for e in [e_minus, e_plus] {
                let r = add_to_edge_list(&mut fixed, e, 0);
                if r != 0 {
                    return r;
                }
            }
        } else if e_minus >= 0 {
            bns.edge[e_minus as usize].forbidden |= mask;
            let r = add_to_edge_list(&mut fixed, e_minus, 0);
            if r != 0 {
                return r;
            }
        } else if e_plus >= 0 {
            bns.edge[e_plus as usize].forbidden |= mask;
            let r = add_to_edge_list(&mut fixed, e_plus, 0);
            if r != 0 {
                return r;
            }
        }
    }
    for i in 0..bns.num_bonds as usize {
        bns.edge[i].forbidden |= mask;
        let r = add_to_edge_list(&mut fixed, i as i32, 0);
        if r != 0 {
            return r;
        }
    }
    let ret = run_bns_restore_once(bns, bd);
    remove_forbidden_edge_mask(bns, &fixed, forbidden_edge_mask);
    if ret < 0 {
        return ret;
    }
    ret
}

/// `IncrementZeroOrderBondsToHeteroat`: returns the number of changes;
/// `num_run_bns` and `total_delta` are incremented
#[allow(clippy::too_many_arguments)]
pub fn increment_zero_order_bonds_to_heteroat(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    const FIX_BOND_ADD_ALLOC: i32 = 128;
    let num_at = ps.num_atoms as usize;
    let len_at = (ps.num_atoms + ps.num_deleted_h) as usize;
    let mask = forbidden_edge_mask as i8;
    let mut carbon_edges = EdgeList::default();
    let mut fixed = EdgeList::default();
    let mut num_changes = 0;
    let mut fixed_carbon_charges = false;
    let mut res = BnsTestResult::default();

    let ret = (|| -> i32 {
        let k = tcg.n_group[TCG_ME_FLOWER0 as usize];
        if tcg.num_metal_atoms == 0 || k < 0 || tcg.p_tcg[k as usize].n_vertex_number < 0 {
            return 0;
        }
        let v_me_flower0 = tcg.p_tcg[k as usize].n_vertex_number as usize;
        let refresh = |ps: &mut StrFromInChI, at2: &mut Vec<InpAtom>, bns: &BnStruct| -> i32 {
            at2.clear();
            at2.extend_from_slice(&ps.at[..len_at]);
            std::mem::swap(&mut ps.at, at2);
            let r = copy_bns_to_atom(ps, bns, p_va, tcg, true);
            std::mem::swap(&mut ps.at, at2);
            r
        };
        let r = refresh(ps, at2, bns);
        if r < 0 {
            return r;
        }
        for i in 0..num_at {
            if p_va[i].c_metal == 0 || p_va[i].n_metal_group_edge <= 0 {
                continue;
            }
            let pe_me_flower = (p_va[i].n_metal_group_edge - 1) as usize;
            if v_me_flower0 != (bns.edge[pe_me_flower].neighbor12 as usize ^ i) {
                return RI_ERR_PROGR;
            }
            let mut v_neigh_neigh = 0usize;
            let mut pe_neigh_neigh = 0usize;
            for j in 0..at2[i].valence as usize {
                if bns.edge[pe_me_flower].flow == 0 {
                    break;
                }
                if i32::from(at2[i].bond_type[j]) & BOND_TYPE_MASK != 0 {
                    continue;
                }
                if !fixed_carbon_charges {
                    let r = forbid_carbon_charge_edges(
                        bns,
                        tcg,
                        &mut carbon_edges,
                        forbidden_edge_mask,
                    );
                    if r < 0 {
                        return r;
                    }
                    fixed_carbon_charges = true;
                }
                let pe_zero = bns.iedge[bns.vert[i].iedge + j] as usize;
                if bns.edge[pe_zero].flow != 0 {
                    return RI_ERR_PROGR;
                }
                for k in 0..at2[i].valence as usize {
                    let pe = bns.iedge[bns.vert[i].iedge + k] as usize;
                    if bns.edge[pe].flow == 1 && bns.edge[pe].forbidden & mask == 0 {
                        let r = add_to_edge_list(&mut fixed, pe as i32, FIX_BOND_ADD_ALLOC);
                        if r != 0 {
                            return r;
                        }
                        bns.edge[pe].forbidden |= mask;
                    }
                }
                for k in 0..num_at {
                    if p_va[k].c_metal == 0
                        && p_va[k].c_num_valence_electrons == 5
                        && at2[k].valence == 2
                        && at2[k].num_h == 0
                        && p_va[k].c_min_ring_size <= 6
                        && p_va[k].n_c_plus_group_edge > 0
                    {
                        let pe = (p_va[k].n_c_plus_group_edge - 1) as usize;
                        if bns.edge[pe].flow == 1 && bns.edge[pe].forbidden & mask == 0 {
                            let r = add_to_edge_list(&mut fixed, pe as i32, FIX_BOND_ADD_ALLOC);
                            if r != 0 {
                                return r;
                            }
                            bns.edge[pe].forbidden |= mask;
                        }
                    }
                }
                let v_neigh = at2[i].neighbor[j] as usize;
                let mut found = false;
                let mut k = bns.vert[v_neigh].num_adj_edges as isize - 1;
                while k >= 0 {
                    let pe = bns.iedge[bns.vert[v_neigh].iedge + k as usize] as usize;
                    if bns.edge[pe].flow != 0 {
                        let vnm = bns.edge[pe].neighbor12 as usize ^ v_neigh;
                        if vnm != i && vnm != v_me_flower0 {
                            pe_neigh_neigh = pe;
                            v_neigh_neigh = vnm;
                            found = true;
                            break;
                        }
                    }
                    k -= 1;
                }
                if !found {
                    continue;
                }
                bns.edge[pe_zero].flow += 1;
                bns.edge[pe_zero].forbidden |= mask;
                bns.edge[pe_me_flower].flow -= 1;
                bns.edge[pe_neigh_neigh].flow -= 1;
                bns.vert[v_me_flower0].st_edge.flow -= 1;
                bns.vert[v_neigh_neigh].st_edge.flow -= 1;
                bns.tot_st_flow -= 2;
                let mut success = 0;
                let mut ret = run_bns_test_once(bns, bd, p_va, &mut res);
                let vnn = v_neigh_neigh as i32;
                if ret == 1
                    && ((res.v_last == v_me_flower0 as i32 && res.v_first == vnn)
                        || (res.v_last == vnn && res.v_first == v_me_flower0 as i32))
                    && res.delta_charge.abs() <= 2
                {
                    ret = run_bns_restore_once(bns, bd);
                    if ret > 0 {
                        *num_run_bns += 1;
                        *total_delta += ret;
                        num_changes += 1;
                        success = ret;
                    }
                    let r = add_to_edge_list(&mut fixed, pe_zero as i32, FIX_BOND_ADD_ALLOC);
                    if r != 0 {
                        return r;
                    }
                } else {
                    bns.edge[pe_zero].flow -= 1;
                    bns.edge[pe_zero].forbidden &= !mask;
                    bns.edge[pe_me_flower].flow += 1;
                    bns.edge[pe_neigh_neigh].flow += 1;
                    bns.vert[v_me_flower0].st_edge.flow += 1;
                    bns.vert[v_neigh_neigh].st_edge.flow += 1;
                    bns.tot_st_flow += 2;
                }
                remove_forbidden_edge_mask(bns, &fixed, forbidden_edge_mask);
                fixed.num_edges = 0;
                if success != 0 {
                    let r = refresh(ps, at2, bns);
                    if r < 0 {
                        return r;
                    }
                }
            }
        }
        num_changes
    })();
    remove_forbidden_edge_mask(bns, &carbon_edges, forbidden_edge_mask);
    remove_forbidden_edge_mask(bns, &fixed, forbidden_edge_mask);
    ret
}

/// `MovePlusFromS2DiaminoCarbon`: C=S(+) with two non-tautomeric NH2 on C
/// becomes C(+)-S
#[allow(clippy::too_many_arguments)]
pub fn move_plus_from_s2_diamino_carbon(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    num_run_bns: &mut i32,
    _total_delta: &mut i32,
    forbidden_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let len_at = (ps.num_atoms + ps.num_deleted_h) as usize;
    let mut all_charge_edges = EdgeList::default();
    let mut res = BnsTestResult::default();
    let mut ret;
    at2.clear();
    at2.extend_from_slice(&ps.at[..len_at]);
    std::mem::swap(&mut ps.at, at2);
    let ret2 = copy_bns_to_atom(ps, bns, p_va, tcg, true);
    std::mem::swap(&mut ps.at, at2);
    if ret2 < 0 {
        return ret2;
    }
    ret = 0;
    for i in 0..num_at {
        let pv_s = bns.vert[i];
        if p_va[i].c_metal != 0
            || p_va[i].c_num_valence_electrons != 6
            || at2[i].valence != 2
            || pv_s.st_edge.cap != pv_s.st_edge.flow
        {
            continue;
        }
        let e_plus_s = p_va[i].n_c_plus_group_edge - 1;
        if e_plus_s < 0 || bns.edge[e_plus_s as usize].flow != 0 {
            continue;
        }
        let pe1 = bns.iedge[pv_s.iedge] as usize;
        let pe2 = bns.iedge[pv_s.iedge + 1] as usize;
        if bns.edge[pe1].flow + bns.edge[pe2].flow != 1 {
            continue;
        }
        let pe_sc = if bns.edge[pe1].flow != 0 { pe1 } else { pe2 };
        let vc = bns.edge[pe_sc].neighbor12 as usize ^ i;
        if p_va[vc].c_num_valence_electrons != 4 || at2[vc].valence != 3 {
            continue;
        }
        let e_plus_c = p_va[vc].n_c_plus_group_edge - 1;
        if e_plus_c < 0 || bns.edge[e_plus_c as usize].flow == 0 {
            continue;
        }
        let e_minus_c = p_va[vc].n_c_minus_group_edge - 1;
        if e_minus_c >= 0 && bns.edge[e_minus_c as usize].flow != 0 {
            continue;
        }
        let pv_c = bns.vert[vc];
        let mut pe_cn: Vec<usize> = Vec::new();
        for j in 0..ps.at[vc].valence as usize {
            let e = bns.iedge[pv_c.iedge + j] as usize;
            if e != pe_sc && bns.edge[e].flow == 0 {
                pe_cn.push(e);
            }
        }
        if pe_cn.len() != 2 {
            continue;
        }
        let mut j = 0;
        while j < 2 {
            let vn = bns.edge[pe_cn[j]].neighbor12 as usize ^ vc;
            if p_va[vn].c_num_valence_electrons != 5
                || bns.vert[vn].st_edge.cap != bns.vert[vn].st_edge.flow
                || at2[vn].num_h != 2
                || at2[vn].endpoint != 0
                || ps.endpoint.as_ref().is_some_and(|e| e[vn] != 0)
            {
                break;
            }
            j += 1;
        }
        if j != 2 {
            continue;
        }
        if all_charge_edges.num_edges == 0 {
            for jj in 0..num_at {
                let e = p_va[jj].n_c_plus_group_edge - 1;
                if e >= 0 && bns.edge[e as usize].forbidden == 0 {
                    ret = add_to_edge_list(&mut all_charge_edges, e, 2 * num_at as i32);
                    if ret != 0 {
                        return ret;
                    }
                }
                let e = p_va[jj].n_c_minus_group_edge - 1;
                if e >= 0 && bns.edge[e as usize].forbidden == 0 {
                    ret = add_to_edge_list(&mut all_charge_edges, e, 2 * num_at as i32);
                    if ret != 0 {
                        return ret;
                    }
                }
            }
        }
        set_forbidden_edge_mask(bns, &all_charge_edges, forbidden_edge_mask);
        bns.edge[e_plus_s as usize].forbidden &= !forbidden_edge_mask as i8;
        let pe = e_plus_c as usize;
        if bns.edge[pe].flow == 0 {
            continue;
        }
        let delta = 1;
        let v1 = bns.edge[pe].neighbor1 as usize;
        let v2 = bns.edge[pe].neighbor12 as usize ^ v1;
        bns.edge[pe].flow -= delta;
        bns.vert[v1].st_edge.flow -= delta;
        bns.vert[v2].st_edge.flow -= delta;
        bns.tot_st_flow -= 2 * delta;
        ret = run_bns_test_once(bns, bd, p_va, &mut res);
        if ret == 1
            && ((res.v_last == v1 as i32 && res.v_first == v2 as i32)
                || (res.v_last == v2 as i32 && res.v_first == v1 as i32))
            && res.delta_charge == -1
        {
            ret = run_bns_restore_once(bns, bd);
            if ret > 0 {
                *num_run_bns += 1;
            }
        } else {
            bns.edge[pe].flow += delta;
            bns.vert[v1].st_edge.flow += delta;
            bns.vert[v2].st_edge.flow += delta;
            bns.tot_st_flow += 2 * delta;
        }
        remove_forbidden_edge_mask(bns, &all_charge_edges, forbidden_edge_mask);
    }
    ret
}

/// `EliminateChargeSeparationOnHeteroatoms`: cancels pairs of opposite
/// charges; returns the number of pairs left or an error
#[allow(clippy::too_many_arguments)]
pub fn eliminate_charge_separation_on_heteroatoms(
    bns: &mut BnStruct,
    bd: &mut BnData,
    ps: &mut StrFromInChI,
    at2: &mut Vec<InpAtom>,
    p_va: &[ValAt],
    tcg: &AllTcGroups,
    num_run_bns: &mut i32,
    total_delta: &mut i32,
    forbidden_edge_mask: i32,
    forbidden_stereo_edge_mask: i32,
) -> i32 {
    let num_at = ps.num_atoms as usize;
    let len_at = (ps.num_atoms + ps.num_deleted_h) as usize;
    let mask = forbidden_edge_mask as i8;
    let mut ring_edges = EdgeList::default();
    let mut carbon_edges = EdgeList::default();
    let mut fixed_carbon = false;
    let mut num_min = 0;
    let mut res = BnsTestResult::default();

    let mut ret = (|| -> i32 {
        if forbidden_stereo_edge_mask != 0 {
            for i in 0..num_at {
                for j in 0..at2[i].valence as usize {
                    let k = bns.iedge[bns.vert[i].iedge + j];
                    if i32::from(bns.edge[k as usize].forbidden) == forbidden_stereo_edge_mask {
                        let q = ps.pbfsq.q.as_mut().expect("BFS queue");
                        let min_ring = is_bond_in_nmax_memb_ring(
                            at2,
                            i,
                            j,
                            q,
                            &mut ps.pbfsq.n_atom_level,
                            &mut ps.pbfsq.c_source,
                            99,
                        );
                        if 0 < min_ring {
                            let r = add_to_edge_list(&mut ring_edges, k, 64);
                            if r != 0 {
                                return r;
                            }
                        }
                    }
                }
            }
            if ring_edges.num_edges == 0 {
                return 0;
            }
            remove_forbidden_edge_mask(bns, &ring_edges, forbidden_stereo_edge_mask);
        }
        at2.clear();
        at2.extend_from_slice(&ps.at[..len_at]);
        std::mem::swap(&mut ps.at, at2);
        let ret2 = copy_bns_to_atom(ps, bns, p_va, tcg, true);
        std::mem::swap(&mut ps.at, at2);
        if ret2 < 0 {
            return ret2;
        }
        let (mut num_pos, mut num_neg) = (0, 0);
        for i in 0..num_at {
            if p_va[i].c_metal == 0 && at2[i].radical == 0 {
                num_pos += i32::from(at2[i].charge > 0);
                num_neg += i32::from(at2[i].charge < 0);
            }
        }
        num_min = num_pos.min(num_neg);
        let k = tcg.n_group[TCG_PLUS as usize];
        if num_min == 0 || k < 0 {
            return finish_keep(ps, at2, len_at);
        }
        let e_plus_super = tcg.p_tcg[k as usize].n_forward_edge;
        let v_plus_super = tcg.p_tcg[k as usize].n_vertex_number;
        if e_plus_super <= 0
            || v_plus_super < num_at as i32
            || bns.edge[e_plus_super as usize].forbidden != 0
        {
            return finish_keep(ps, at2, len_at);
        }
        let e = e_plus_super as usize;
        let vps = v_plus_super as usize;
        let vpm = bns.edge[e].neighbor12 as usize ^ vps;
        num_min = num_min.min(bns.edge[e].flow);
        let mut d_plus = bns.vert[vps].st_edge.cap - bns.vert[vps].st_edge.flow;
        let mut d_minus = bns.vert[vpm].st_edge.cap - bns.vert[vpm].st_edge.flow;
        if num_min != 0 && d_plus == 0 && d_minus == 0 {
            if !fixed_carbon {
                let r =
                    forbid_carbon_charge_edges(bns, tcg, &mut carbon_edges, forbidden_edge_mask);
                if r < 0 {
                    return r;
                }
                fixed_carbon = true;
            }
            let delta = 1;
            bns.edge[e].forbidden |= mask;
            bns.edge_forbidden_mask |= mask;
            let mut i = 0;
            while i < num_min {
                bns.edge[e].flow -= delta;
                bns.vert[vps].st_edge.flow -= delta;
                bns.vert[vpm].st_edge.flow -= delta;
                bns.tot_st_flow -= 2 * delta;
                let r = run_bns_test_once(bns, bd, p_va, &mut res);
                if r < 0 {
                    return r;
                }
                if r == 1
                    && ((res.v_last == vps as i32 && res.v_first == vpm as i32)
                        || (res.v_last == vpm as i32 && res.v_first == vps as i32))
                    && res.delta_charge < 0
                {
                    let r = run_bns_restore_once(bns, bd);
                    *num_run_bns += 1;
                    if r < 0 {
                        return r;
                    } else if r == 1 {
                        *total_delta += r;
                    } else {
                        return RI_ERR_PROGR;
                    }
                } else {
                    bns.edge[e].flow += delta;
                    bns.vert[vps].st_edge.flow += delta;
                    bns.vert[vpm].st_edge.flow += delta;
                    bns.tot_st_flow += 2 * delta;
                    break;
                }
                i += delta;
            }
            num_min -= i;
            bns.edge[e].forbidden &= !mask;
        }
        d_plus = bns.vert[vps].st_edge.cap - bns.vert[vps].st_edge.flow;
        d_minus = bns.vert[vpm].st_edge.cap - bns.vert[vpm].st_edge.flow;
        if num_min > 1 && d_plus == 0 && d_minus == 0 {
            let delta = 2;
            bns.edge[e].forbidden |= mask;
            bns.edge_forbidden_mask |= mask;
            let mut i = 0;
            while i < num_min {
                bns.vert[vps].st_edge.cap += delta;
                bns.tot_st_cap += delta;
                let r = run_bns_test_once(bns, bd, p_va, &mut res);
                if r < 0 {
                    return r;
                }
                bns.vert[vps].st_edge.cap -= delta;
                bns.tot_st_cap -= delta;
                if r != 1
                    || res.v_last != vps as i32
                    || res.v_first != vps as i32
                    || res.delta_charge >= 0
                {
                    break;
                }
                bns.vert[vpm].st_edge.cap += delta;
                bns.tot_st_cap += delta;
                let r = run_bns_test_once(bns, bd, p_va, &mut res);
                if r < 0 {
                    return r;
                }
                bns.vert[vpm].st_edge.cap -= delta;
                bns.tot_st_cap -= delta;
                if r != 1
                    || res.v_last != vpm as i32
                    || res.v_first != vpm as i32
                    || res.delta_charge >= 0
                {
                    break;
                }
                bns.edge[e].flow -= delta;
                bns.vert[vps].st_edge.flow -= delta;
                bns.vert[vpm].st_edge.flow -= delta;
                bns.tot_st_flow -= 2 * delta;
                let r = run_bns_restore_once(bns, bd);
                *num_run_bns += 1;
                if r < 0 {
                    return r;
                } else if r == 2 {
                    *total_delta += r;
                } else {
                    return RI_ERR_PROGR;
                }
                i += delta;
            }
            num_min -= i;
            bns.edge[e].forbidden &= !mask;
        }
        finish_keep(ps, at2, len_at)
    })();
    if fixed_carbon {
        remove_forbidden_edge_mask(bns, &carbon_edges, forbidden_edge_mask);
    }
    if forbidden_stereo_edge_mask != 0 && ring_edges.num_edges != 0 {
        set_forbidden_edge_mask(bns, &ring_edges, forbidden_stereo_edge_mask);
    }
    if ret >= 0 {
        ret = num_min;
    }
    ret
}

/// the C code copies `at` over `at2` at the end of
/// `EliminateChargeSeparationOnHeteroatoms`
fn finish_keep(ps: &StrFromInChI, at2: &mut Vec<InpAtom>, len_at: usize) -> i32 {
    at2.clear();
    at2.extend_from_slice(&ps.at[..len_at]);
    0
}
