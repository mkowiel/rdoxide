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
use super::ichirvr1::CN_LIST;
use super::ichirvr1::{
    bond_flow_maxcap_minorder, comp_cc_cand, get_tgroup_info_from_inchi, run_bns_restore_once,
};
use super::ichirvrs::{
    AllTcGroups, BNS_VT_C_POS, CN_BITS_SHIFT, CcCand, Srm, TGRF_MINUS_FIRST, ValAt,
};
use super::ichisort::insertions_sort;
use super::ichister::{
    AB_PARITY_EVEN, AB_PARITY_ODD, MAX_CUMULENE_LEN, MAX_NUM_STEREO_BOND_NEIGH,
    MIN_NUM_STEREO_BOND_NEIGH, atom_parity_well_def, b_can_atom_be_middle_allene,
    b_can_atom_be_terminal_allene,
};
use super::ichitaut::TGroupInfo;
use super::incomdef::RADICAL_SINGLET;
use super::incomdef::{
    BOND_TYPE_ALTERN, BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, NUM_H_ISOTOPES,
};
use super::inpdef::{InpAtom, MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS};
use super::util::{EL_NUMBER_H, get_element_chemical_symbol};

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
