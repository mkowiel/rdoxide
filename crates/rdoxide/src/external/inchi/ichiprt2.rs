//! String construction for the InChI output (`INCHI_BASE/src/ichiprt2.c`):
//! "alphabetic" and decimal numbers, the connection-table, H, tautomer,
//! charge/radical/valence, equivalence, isotopic and stereo substrings, and the
//! comparisons that decide whether a layer is repeated.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! The strings of 2048 bytes (`szValue`) of the C code never overflow here
//! (the C overflow branches are unreachable for such values), so the helpers
//! return `String`s.

use super::ichi::{InChI, InChIAux, InChIIsotopicAtom, InChIIsotopicTGroup, OrigInfo};
use super::ichi_io::{InchiIosString, inchi_strbuf_printf};
use super::ichimake::{CT_MODE_PREDECESSORS, get_dfs_order_4ct};
use super::inpdef::MAX_ATOMS;

pub const ITEM_DELIMETER: &str = ",";
pub const EXTRA_SPACE: &str = "";
pub const COMMA_EXTRA_SPACE: &str = ",";
pub const CT_MODE_NO_ORPHANS: i32 = 1;
pub const CT_MODE_ABC_NUMBERS: i32 = 2;
pub const CT_MODE_ATOM_COUNTS: i32 = 4;
pub const CT_MODE_EQL_H_TOGETHER: i32 = 16;
pub const CT_MODE_ABC_NUM_CLOSURES: i32 = 32;

pub const EQL_EQU: i32 = 0;
pub const EQL_EQU_TG: i32 = 1;
pub const EQL_EQU_ISO: i32 = 2;
pub const EQL_NUM: i32 = 0;
pub const EQL_NUM_INV: i32 = 1;
pub const EQL_NUM_ISO: i32 = 2;

const INCHI_T_NUM_MOVABLE: usize = 2;
const NUM_H_ISOTOPES: usize = 3;

const ALPHA_BASE: i32 = 27;

/// `MakeAbcNumber` into a buffer: base-27 "alphabetic" number with the
/// leading delimiter; returns the length or -1
pub fn make_abc_number(
    buf: &mut [u8],
    mut n_string_len: i32,
    delim: Option<&str>,
    mut value: i32,
) -> i32 {
    if n_string_len < 2 {
        return -1;
    }
    let mut p = 0usize;
    for &c in delim.unwrap_or("").as_bytes() {
        n_string_len -= 1;
        if n_string_len == 0 {
            break;
        }
        buf[p] = c;
        p += 1;
    }
    if n_string_len < 2 {
        return -1;
    }
    if value == 0 {
        buf[p] = b'.';
        buf[p + 1] = 0;
        return 1;
    }
    if value < 0 {
        buf[p] = b'-';
        p += 1;
        n_string_len -= 1;
        value = value.wrapping_neg();
    }
    let mut q = p;
    while value != 0 {
        n_string_len -= 1;
        if n_string_len == 0 {
            break;
        }
        let ch = value % ALPHA_BASE;
        buf[q] = if ch != 0 { b'a' + ch as u8 - 1 } else { b'@' };
        q += 1;
        value /= ALPHA_BASE;
    }
    if n_string_len <= 0 {
        return -1;
    }
    buf[q] = 0;
    buf[p..q].reverse();
    buf[p] = buf[p].to_ascii_uppercase();
    q as i32
}

/// `MakeDecNumber` into a buffer
pub fn make_dec_number(
    buf: &mut [u8],
    mut n_string_len: i32,
    delim: Option<&str>,
    mut value: i32,
) -> i32 {
    if n_string_len < 2 {
        return -1;
    }
    let mut p = 0usize;
    for &c in delim.unwrap_or("").as_bytes() {
        n_string_len -= 1;
        if n_string_len == 0 {
            break;
        }
        buf[p] = c;
        p += 1;
    }
    if n_string_len < 2 {
        return -1;
    }
    if value == 0 {
        buf[p] = b'0';
        buf[p + 1] = 0;
        return (p + 1) as i32;
    }
    if value < 0 {
        buf[p] = b'-';
        p += 1;
        n_string_len -= 1;
        value = value.wrapping_neg();
    }
    let mut q = p;
    while value != 0 {
        n_string_len -= 1;
        if n_string_len == 0 {
            break;
        }
        buf[q] = b'0' + (value % 10) as u8;
        q += 1;
        value /= 10;
    }
    if n_string_len <= 0 {
        return -1;
    }
    buf[q] = 0;
    buf[p..q].reverse();
    q as i32
}

/// the number as a string (the 2048-byte `szValue` never overflows)
pub fn abc(delim: Option<&str>, v: i32) -> String {
    let mut b = [0u8; 2048];
    let n = make_abc_number(&mut b, 2048, delim, v);
    String::from_utf8_lossy(&b[..n.max(0) as usize]).into_owned()
}

/// the number as a string
pub fn dec(delim: Option<&str>, v: i32) -> String {
    let mut b = [0u8; 2048];
    let n = make_dec_number(&mut b, 2048, delim, v);
    String::from_utf8_lossy(&b[..n.max(0) as usize]).into_owned()
}

fn num(abc_numbers: bool, delim: Option<&str>, v: i32) -> String {
    if abc_numbers {
        abc(delim, v)
    } else {
        dec(delim, v)
    }
}

fn print(buf: &mut InchiIosString, s: &str) {
    inchi_strbuf_printf(buf, s);
}

fn used(buf: &InchiIosString) -> i32 {
    buf.n_used_length
}

/// `Eql_INChI_Isotopic`
pub fn eql_inchi_isotopic(i1: Option<&InChI>, i2: Option<&InChI>) -> bool {
    let (Some(i1), Some(i2)) = (i1, i2) else {
        return false;
    };
    let na = i1.n_number_of_isotopic_atoms.max(0) as usize;
    let nt = i1.n_number_of_isotopic_t_groups.max(0) as usize;
    i1.b_deleted == 0
        && i2.b_deleted == 0
        && (na > 0 || nt > 0)
        && i1.n_number_of_isotopic_atoms == i2.n_number_of_isotopic_atoms
        && i1.n_number_of_isotopic_t_groups == i2.n_number_of_isotopic_t_groups
        && (na == 0
            || (!i1.isotopic_atom.is_empty()
                && !i2.isotopic_atom.is_empty()
                && i1.isotopic_atom[..na] == i2.isotopic_atom[..na]))
        && (nt == 0
            || (!i1.isotopic_t_group.is_empty()
                && !i2.isotopic_t_group.is_empty()
                // the C code compares `nt` records of the (larger) atom
                // structure; the t-group records are 8 bytes against 10
                && tgroups_prefix_eq(&i1.isotopic_t_group, &i2.isotopic_t_group, nt)))
}

fn tgroups_prefix_eq(a: &[InChIIsotopicTGroup], b: &[InChIIsotopicTGroup], n: usize) -> bool {
    // memcmp over n * sizeof(INChI_IsotopicAtom) = n * 10 bytes of arrays of
    // 8-byte records: covers ceil(n * 10 / 8) records
    let bytes = n * std::mem::size_of::<InChIIsotopicAtomC>();
    let recs = bytes.div_ceil(8).min(a.len()).min(b.len());
    let full = bytes / 8;
    if a[..full.min(recs)] != b[..full.min(recs)] {
        return false;
    }
    // a partial trailing record compares only its leading bytes
    let rem = bytes % 8;
    if rem == 0 || full >= recs {
        return true;
    }
    let tb = |t: &InChIIsotopicTGroup| {
        let mut v = [0u8; 8];
        v[0..2].copy_from_slice(&t.n_t_group_number.to_ne_bytes());
        v[2..4].copy_from_slice(&t.n_num_h.to_ne_bytes());
        v[4..6].copy_from_slice(&t.n_num_d.to_ne_bytes());
        v[6..8].copy_from_slice(&t.n_num_t.to_ne_bytes());
        v
    };
    tb(&a[full])[..rem] == tb(&b[full])[..rem]
}

/// the C `INChI_IsotopicAtom` (10 bytes: `AT_NUMB` + 4 `NUM_H`)
#[allow(dead_code)]
struct InChIIsotopicAtomC([u8; 10]);

/// `memcmp` of the first `len` elements; the C arrays are longer than their
/// meaningful part (zero-filled), so missing elements read as 0
fn same_prefix(a: &[u16], b: &[u16], len: usize) -> bool {
    (0..len).all(|i| a.get(i).copied().unwrap_or(0) == b.get(i).copied().unwrap_or(0))
}

/// `Eql_INChI_Aux_Equ`
pub fn eql_inchi_aux_equ(
    a1: Option<&InChIAux>,
    eql1: i32,
    a2: Option<&InChIAux>,
    eql2: i32,
) -> bool {
    let (Some(a1), Some(a2)) = (a1, a2) else {
        return false;
    };
    let (t1, t2) = (eql1 & EQL_EQU_TG != 0, eql2 & EQL_EQU_TG != 0);
    let mut len = 0usize;
    let (mut n1, mut n2): (Option<&[u16]>, Option<&[u16]>) = (None, None);
    if t1 && t2 {
        if a1.n_number_of_t_groups > 0
            && a1.n_number_of_t_groups == a2.n_number_of_t_groups
            && a1.b_deleted == 0
            && a2.b_deleted == 0
        {
            len = a1.n_number_of_t_groups as usize;
            if eql1 & EQL_EQU_ISO != 0 {
                if a1.b_is_isotopic != 0 {
                    n1 = Some(&a1.n_constit_equ_isotopic_t_group_numbers);
                }
            } else {
                n1 = Some(&a1.n_constit_equ_t_group_numbers);
            }
            if eql2 & EQL_EQU_ISO != 0 {
                if a2.b_is_isotopic != 0 {
                    n2 = Some(&a2.n_constit_equ_isotopic_t_group_numbers);
                }
            } else {
                n2 = Some(&a2.n_constit_equ_t_group_numbers);
            }
        }
    } else if !t1
        && !t2
        && a1.n_number_of_atoms > 0
        && a1.n_number_of_atoms == a2.n_number_of_atoms
        && a1.b_deleted == 0
        && a2.b_deleted == 0
    {
        len = a1.n_number_of_atoms as usize;
        if eql1 & EQL_EQU_ISO != 0 {
            if a1.b_is_isotopic != 0 {
                n1 = Some(&a1.n_constit_equ_isotopic_numbers);
            }
        } else {
            n1 = Some(&a1.n_constit_equ_numbers);
        }
        if eql2 & EQL_EQU_ISO != 0 {
            if a2.b_is_isotopic != 0 {
                n2 = Some(&a2.n_constit_equ_isotopic_numbers);
            }
        } else {
            n2 = Some(&a2.n_constit_equ_numbers);
        }
    }
    match (n1, n2) {
        (Some(n1), Some(n2)) if !n1.is_empty() && !n2.is_empty() => {
            same_prefix(n1, n2, len) && has_equ_string(n1, len as i32)
        }
        _ => false,
    }
}

/// `Eql_INChI_Aux_Num`
pub fn eql_inchi_aux_num(
    a1: Option<&InChIAux>,
    eql1: i32,
    a2: Option<&InChIAux>,
    eql2: i32,
) -> bool {
    let (Some(a1), Some(a2)) = (a1, a2) else {
        return false;
    };
    let len = a1.n_number_of_atoms;
    if len <= 0 || len != a2.n_number_of_atoms || a1.b_deleted != 0 || a2.b_deleted != 0 {
        return false;
    }
    if (eql1 & EQL_NUM_ISO != 0 && a1.b_is_isotopic == 0)
        || (eql2 & EQL_NUM_ISO != 0 && a2.b_is_isotopic == 0)
    {
        return false;
    }
    fn pick(a: &InChIAux, e: i32) -> Option<&[u16]> {
        match e {
            EQL_NUM => Some(&a.n_orig_at_nos_in_canon_ord),
            EQL_NUM_ISO => Some(&a.n_isotopic_orig_at_nos_in_canon_ord),
            EQL_NUM_INV => Some(&a.n_orig_at_nos_in_canon_ord_inv),
            x if x == (EQL_NUM_INV | EQL_NUM_ISO) => {
                Some(&a.n_isotopic_orig_at_nos_in_canon_ord_inv)
            }
            _ => None,
        }
    }
    let (Some(n1), Some(n2)) = (pick(a1, eql1), pick(a2, eql2)) else {
        return false;
    };
    let len = len as usize;
    !n1.is_empty() && !n2.is_empty() && same_prefix(n1, n2, len)
}

/// `bHasOrigInfo`
pub fn has_orig_info(orig_info: &[OrigInfo], num_atoms: i32) -> bool {
    num_atoms > 0
        && !orig_info.is_empty()
        && orig_info[..(num_atoms as usize).min(orig_info.len())]
            .iter()
            .any(|o| o.c_charge != 0 || o.c_radical != 0 || o.c_unusual_valence != 0)
}

/// `EqlOrigInfo`
pub fn eql_orig_info(a1: Option<&InChIAux>, a2: Option<&InChIAux>) -> bool {
    let (Some(a1), Some(a2)) = (a1, a2) else {
        return false;
    };
    let n = a1.n_number_of_atoms;
    a1.n_number_of_atoms == a2.n_number_of_atoms
        && has_orig_info(&a1.orig_info, n)
        && !a2.orig_info.is_empty()
        && (0..n as usize).all(|i| {
            // the C arrays are zero-padded past their end
            a1.orig_info.get(i).copied().unwrap_or_default()
                == a2.orig_info.get(i).copied().unwrap_or_default()
        })
}

/// `bHasEquString`
pub fn has_equ_string(linear_ct: &[u16], n_len_ct: i32) -> bool {
    if linear_ct.is_empty() {
        return false;
    }
    let n = n_len_ct.max(0) as usize;
    for k in 0..n {
        if k as i32 != i32::from(linear_ct[k]) - 1 {
            continue;
        }
        for i in k..n {
            if k as i32 != i32::from(linear_ct[i]) - 1 {
                continue;
            }
            if k < i {
                return true;
            }
        }
    }
    false
}

/// `MakeMult`
pub fn make_mult(
    mult: i32,
    tailing_delim: &str,
    buf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    if mult == 1 || *overflow != 0 {
        return 0;
    }
    let s = num(ct_mode & CT_MODE_ABC_NUMBERS != 0, None, mult);
    let t = format!("{s}{tailing_delim}");
    inchi_strbuf_printf(buf, &t)
}

/// `MakeDelim`
pub fn make_delim(
    tailing_delim: Option<&str>,
    buf: &mut InchiIosString,
    overflow: &mut i32,
) -> i32 {
    match tailing_delim {
        Some(d) if !d.is_empty() && *overflow == 0 => inchi_strbuf_printf(buf, d),
        _ => 0,
    }
}

/// `MakeEqStr`
pub fn make_eq_str(
    tailing_delim: Option<&str>,
    mult: i32,
    buf: &mut InchiIosString,
    overflow: &mut i32,
) -> i32 {
    let Some(d) = tailing_delim.filter(|d| !d.is_empty()) else {
        return 0;
    };
    if *overflow != 0 {
        return 0;
    }
    let n0 = used(buf);
    if mult != 1 {
        let s = dec(None, mult);
        inchi_strbuf_printf(buf, &s);
    }
    inchi_strbuf_printf(buf, d);
    used(buf) - n0
}

/// `MakeCtStringNew`
pub fn make_ct_string_new(
    linear_ct: &[u16],
    n_len_ct: i32,
    add_delim: bool,
    n_num_h: Option<&[i8]>,
    num_atoms: usize,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let no_num_h = n_num_h.is_none();
    let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    let predecessors = ct_mode & CT_MODE_PREDECESSORS != 0;
    let count_ring_closures =
        abc_numbers && predecessors && ct_mode & CT_MODE_ABC_NUM_CLOSURES != 0;
    if n_len_ct <= 1 {
        return 0;
    }
    let used0 = used(strbuf);
    let Some(dfs) = get_dfs_order_4ct(linear_ct, n_len_ct as usize, n_num_h, num_atoms, ct_mode)
    else {
        *overflow += 1;
        return 0;
    };
    if ovfl == 0 && add_delim {
        print(strbuf, ",");
    }
    if ovfl == 0 {
        let mut num_ring_closures = 0;
        let mut sz = [0u8; 2048];
        let mut i = 0usize;
        while dfs.get(i).copied().unwrap_or(0) != 0 {
            let n_value = if dfs[i] as usize > MAX_ATOMS {
                0
            } else {
                i32::from(dfs[i])
            };
            let num_h = if dfs[i + 1] != 0 {
                i32::from(dfs[i + 1]) - 16
            } else {
                0
            };
            let n_delim = dfs[i + 2] as u8;
            // `szValue` is not cleared between the items: a character put without a
            // terminator (`H`, a delimiter) is followed by the rest of the previous item
            let mut len = 0usize;
            let put_char = |sz: &mut [u8; 2048], len: &mut usize, c: u8| {
                sz[*len] = c;
                *len += 1;
            };
            let put_num = |sz: &mut [u8; 2048], len: &mut usize, digits: &str| {
                let b = digits.as_bytes();
                sz[*len..*len + b.len()].copy_from_slice(b);
                sz[*len + b.len()] = 0;
                *len += b.len();
            };
            if predecessors {
                if count_ring_closures {
                    if n_delim == b'-' && i > 3 && no_num_h {
                        if num_ring_closures == 0 {
                            let mut j = i;
                            while dfs.get(j).copied().unwrap_or(0) != 0
                                && dfs[j + 2] == u16::from(b'-')
                            {
                                num_ring_closures += 1;
                                j += 3;
                            }
                            if num_ring_closures != 0 {
                                put_num(&mut sz, &mut len, &dec(None, num_ring_closures));
                            }
                            num_ring_closures -= 1;
                        } else {
                            num_ring_closures -= 1;
                        }
                    } else {
                        num_ring_closures = 0;
                    }
                } else if n_delim != 0
                    && !(abc_numbers && n_delim == b',')
                    && (n_num_h.is_some() || i > 3)
                {
                    put_char(&mut sz, &mut len, n_delim);
                }
            } else if n_delim != 0 && !(abc_numbers && n_delim == b'-') {
                put_char(&mut sz, &mut len, n_delim);
            }
            if abc_numbers {
                if n_value != 0 || i != 0 {
                    put_num(&mut sz, &mut len, &abc(None, n_value));
                }
                if num_h != 0 {
                    put_num(&mut sz, &mut len, &dec(None, num_h));
                }
            } else {
                if n_value != 0 || i != 0 {
                    put_num(&mut sz, &mut len, &dec(None, n_value));
                }
                if num_h != 0 {
                    put_char(&mut sz, &mut len, b'H');
                    if num_h > 1 {
                        put_num(&mut sz, &mut len, &dec(None, num_h));
                    }
                }
            }
            if len > 0 {
                let end = sz.iter().position(|&c| c == 0).unwrap_or(sz.len());
                print(strbuf, &String::from_utf8_lossy(&sz[..end]));
            }
            i += 3;
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeCtStringOld`
pub fn make_ct_string_old(
    linear_ct: &[u16],
    n_len_ct: i32,
    add_delim: bool,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let used0 = used(strbuf);
    let mut n_max: u16 = 0;
    let mut next = 0;
    if ct_mode & CT_MODE_ABC_NUMBERS == 0 && ovfl == 0 && add_delim {
        print(strbuf, ",");
    }
    if ovfl == 0 {
        let n = n_len_ct.max(0) as usize;
        for i in 0..n {
            let mut less_than_prev = false;
            let take = if ct_mode & CT_MODE_NO_ORPHANS == 0 {
                true
            } else {
                less_than_prev = linear_ct[i] < n_max;
                if less_than_prev {
                    true
                } else if i + 1 < n {
                    n_max = linear_ct[i];
                    linear_ct[i + 1] < n_max
                } else {
                    false
                }
            };
            if take {
                let v = i32::from(linear_ct[i]);
                let s = if ct_mode & CT_MODE_ABC_NUMBERS != 0 {
                    abc(
                        if next == 0 && add_delim {
                            Some(ITEM_DELIMETER)
                        } else {
                            None
                        },
                        v,
                    )
                } else if ct_mode & CT_MODE_NO_ORPHANS != 0 {
                    dec(Some(if less_than_prev { "-" } else { ITEM_DELIMETER }), v)
                } else {
                    dec(if i != 0 { Some(ITEM_DELIMETER) } else { None }, v)
                };
                if !s.is_empty() {
                    print(strbuf, &s);
                    next += 1;
                }
            }
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeCtString`
pub fn make_ct_string(
    linear_ct: &[u16],
    n_len_ct: i32,
    add_delim: bool,
    n_num_h: Option<&[i8]>,
    num_atoms: usize,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    if n_num_h.is_none() || ct_mode & CT_MODE_NO_ORPHANS == 0 {
        make_ct_string_old(linear_ct, n_len_ct, add_delim, strbuf, ct_mode, overflow)
    } else {
        make_ct_string_new(
            linear_ct, n_len_ct, add_delim, n_num_h, num_atoms, strbuf, ct_mode, overflow,
        )
    }
}

/// `MakeHString`
pub fn make_h_string(
    add_delim: bool,
    linear_ct: &[i8],
    n_len_ct: i32,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let used0 = used(strbuf);
    let ovfl = *overflow;
    let abc_nums = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    let mut next = 0;
    if !abc_nums && ovfl == 0 && add_delim {
        print(strbuf, ",");
    }
    if ovfl == 0 && n_len_ct > 0 && !linear_ct.is_empty() {
        let n = n_len_ct as usize;
        if ct_mode & CT_MODE_EQL_H_TOGETHER != 0 {
            let (mut cur_min, mut cur_max) = (-4i32, 16i32);
            let mut n_num_h: Vec<i32>;
            let mut grown = false;
            loop {
                let mut out_of_range = 0;
                let mut tot_no_h = 0;
                n_num_h = vec![0; (cur_max - cur_min + 1) as usize];
                for &x in &linear_ct[..n] {
                    let c = i32::from(x);
                    if c < cur_min {
                        cur_min = c;
                        out_of_range += 1;
                    } else if c > cur_max {
                        cur_max = c;
                        out_of_range += 1;
                    } else if out_of_range == 0 {
                        n_num_h[(c - cur_min) as usize] += 1;
                    }
                    tot_no_h += i32::from(c == 0);
                }
                if tot_no_h as usize == n {
                    return 0;
                }
                if out_of_range != 0 {
                    if grown {
                        *overflow |= 1;
                        return 0;
                    }
                    grown = true;
                    continue;
                }
                break;
            }
            for cur in cur_min..=cur_max {
                let mut num_at = n_num_h[(cur - cur_min) as usize];
                if num_at == 0 || cur == 0 {
                    continue;
                }
                let mut j = 0usize;
                while j < n && num_at != 0 {
                    if cur == i32::from(linear_ct[j]) {
                        j += 1;
                        let i_first = j;
                        num_at -= 1;
                        while j < n && cur == i32::from(linear_ct[j]) && num_at != 0 {
                            num_at -= 1;
                            j += 1;
                        }
                        let mut s = if abc_nums {
                            abc(None, i_first as i32)
                        } else {
                            let s = dec(
                                if next != 0 {
                                    Some(ITEM_DELIMETER)
                                } else {
                                    None
                                },
                                i_first as i32,
                            );
                            next += 1;
                            s
                        };
                        if i_first < j {
                            s += &if abc_nums {
                                abc(None, j as i32)
                            } else {
                                dec(Some("-"), j as i32)
                            };
                        }
                        if num_at == 0 || abc_nums {
                            if abc_nums {
                                s += &dec(None, cur);
                            } else {
                                let ph = if cur > 0 { "H" } else { "h" };
                                let v = cur.abs();
                                if v > 1 {
                                    s += &dec(Some(ph), v);
                                } else {
                                    s += ph;
                                }
                            }
                        }
                        if !s.is_empty() {
                            print(strbuf, &s);
                            next += 1;
                        }
                    } else {
                        j += 1;
                    }
                }
            }
        } else {
            let mut i_first = 0usize;
            for i in 1..=n {
                if i < n && linear_ct[i] == linear_ct[i_first] {
                    continue;
                }
                if linear_ct[i_first] != 0 {
                    let n_val = i32::from(linear_ct[i_first]);
                    let first = i_first + 1;
                    let mut s = if abc_nums {
                        abc(None, first as i32)
                    } else {
                        dec(
                            if next != 0 {
                                Some(ITEM_DELIMETER)
                            } else {
                                None
                            },
                            first as i32,
                        )
                    };
                    if first < i {
                        s += &if abc_nums {
                            abc(None, i as i32)
                        } else {
                            dec(Some("-"), i as i32)
                        };
                    }
                    if abc_nums {
                        s += &dec(None, n_val);
                    } else {
                        let ph = if n_val > 0 { "H" } else { "h" };
                        let v = n_val.abs();
                        if v > 1 {
                            s += &dec(Some(ph), v);
                        } else {
                            s += ph;
                        }
                    }
                    if !s.is_empty() {
                        print(strbuf, &s);
                        next += 1;
                    }
                }
                i_first = i;
            }
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeTautString`
pub fn make_taut_string(
    linear_ct: &[u16],
    n_len_ct: i32,
    add_delim: bool,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    if n_len_ct == 0 || linear_ct.is_empty() || linear_ct[0] == 0 {
        return 0;
    }
    let ovfl = *overflow;
    let used0 = used(strbuf);
    let compressed = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    if !compressed && ovfl == 0 && add_delim {
        print(strbuf, COMMA_EXTRA_SPACE);
    }
    let lin = &linear_ct[1..];
    let n = (n_len_ct - 1).max(0) as usize;
    if ovfl == 0 {
        let (mut group_len, mut group_count) = (0i32, 0i32);
        let mut i = 0usize;
        while i < n {
            let v = i32::from(lin[i]);
            let mut s = String::new();
            if group_len == group_count {
                group_len = v;
                group_count = 0;
                if !compressed {
                    s = if i == 0 { "(".into() } else { ")(".into() };
                }
            } else if compressed && group_count >= INCHI_T_NUM_MOVABLE as i32 {
                s = abc(None, v);
                group_count += 1;
            } else {
                if compressed {
                    let mut p: Option<&str> = None;
                    match group_count {
                        0 => s = dec(if i == 1 { Some(ITEM_DELIMETER) } else { None }, v),
                        1 => p = Some("-"),
                        2 => p = Some("+"),
                        _ => {}
                    }
                    if let Some(p) = p {
                        s = match v {
                            0 => String::new(),
                            1 => p.into(),
                            _ => dec(Some(p), v),
                        };
                    }
                } else if group_count >= INCHI_T_NUM_MOVABLE as i32 {
                    s = dec(Some(ITEM_DELIMETER), v);
                } else if v != 0 {
                    let p = match group_count {
                        0 => Some("H"),
                        1 => Some("-"),
                        2 => Some("+"),
                        _ => None,
                    };
                    if let Some(p) = p {
                        s = if v == 1 { p.into() } else { dec(Some(p), v) };
                    }
                }
                group_count += 1;
            }
            if !s.is_empty() {
                print(strbuf, &s);
            }
            i += 1;
        }
        if !compressed && i != 0 {
            print(strbuf, ")");
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeCRVString`
pub fn make_crv_string(
    orig_info: &[OrigInfo],
    n_len_ct: i32,
    add_delim: bool,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    let used0 = used(strbuf);
    let mut next = 0;
    if ovfl == 0 && add_delim {
        print(strbuf, ", ");
    }
    for k in 0..n_len_ct.max(0) as usize {
        if ovfl != 0 {
            break;
        }
        let o = &orig_info[k];
        let mut s = String::new();
        if o.c_charge != 0 || o.c_radical != 0 || o.c_unusual_valence != 0 {
            let ch = i32::from(o.c_charge);
            let uv = i32::from(o.c_unusual_valence);
            let rad = |r: i8| match r {
                1 => 'd',
                2 => 't',
                _ => 'u',
            };
            if abc_numbers {
                s += &abc(None, k as i32 + 1);
                if ch != 0 {
                    s += &if ch > 0 {
                        dec(Some("+"), ch)
                    } else {
                        dec(None, ch)
                    };
                }
                if o.c_radical != 0 {
                    if ch == 0 {
                        s.push('.');
                    }
                    s.push(rad(o.c_radical));
                }
                if uv != 0 {
                    if ch != 0 && o.c_radical == 0 {
                        s.push('.');
                    }
                    s += &dec(None, uv);
                }
            } else {
                s += &dec(
                    if next != 0 {
                        Some(ITEM_DELIMETER)
                    } else {
                        None
                    },
                    k as i32 + 1,
                );
                if ch != 0 {
                    s += &if ch > 0 {
                        dec(Some("+"), ch)
                    } else {
                        dec(None, ch)
                    };
                }
                if o.c_radical != 0 {
                    s.push(rad(o.c_radical));
                }
                if uv != 0 {
                    if o.c_radical == 0 {
                        s.push('.');
                    }
                    s += &dec(None, uv);
                }
            }
        }
        if !s.is_empty() {
            print(strbuf, &s);
            next += 1;
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeEquString`
pub fn make_equ_string(
    linear_ct: &[u16],
    n_len_ct: i32,
    add_delim: bool,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let used0 = used(strbuf);
    let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    let mut next = 0;
    if ovfl == 0 && add_delim {
        print(strbuf, ", ");
    }
    let n = n_len_ct.max(0) as usize;
    for k in 0..n {
        if ovfl != 0 {
            break;
        }
        if k as i32 != i32::from(linear_ct[k]) - 1 {
            continue;
        }
        for i in k..n {
            if k as i32 != i32::from(linear_ct[i]) - 1 {
                continue;
            }
            let s = if abc_numbers {
                abc(
                    if i == k && next != 0 {
                        Some(ITEM_DELIMETER)
                    } else {
                        None
                    },
                    i as i32 + 1,
                )
            } else {
                dec(
                    Some(if i == k { "(" } else { ITEM_DELIMETER }),
                    i as i32 + 1,
                )
            };
            print(strbuf, &s);
            next += 1;
        }
        print(strbuf, ")");
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeIsoAtomString`
pub fn make_iso_atom_string(
    iso: &[InChIIsotopicAtom],
    n: i32,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let used0 = used(strbuf);
    let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    const LETTER: [char; 4] = ['i', 't', 'd', 'h'];
    const H: [&str; 3] = ["T", "D", "H"];
    const SIGN: [&str; 2] = ["-", "+"];
    if ovfl == 0 {
        for i in 0..n.max(0) as usize {
            let a = &iso[i];
            let mut s = String::new();
            for j in 0..5 {
                let v = match j {
                    0 => i32::from(a.n_atom_number),
                    1 => i32::from(a.n_iso_difference),
                    2 => i32::from(a.n_num_t),
                    3 => i32::from(a.n_num_d),
                    _ => i32::from(a.n_num_h),
                };
                if j == 0 {
                    s += &if abc_numbers {
                        abc(None, v)
                    } else {
                        dec(Some(if i != 0 { ITEM_DELIMETER } else { EXTRA_SPACE }), v)
                    };
                } else if abc_numbers {
                    if j == 1 {
                        s += &dec(None, v);
                    } else if v != 0 {
                        s.push(LETTER[j - 1]);
                        if v != 1 {
                            s += &dec(None, v);
                        }
                    }
                } else if v != 0 {
                    if j == 1 {
                        let subtract = i32::from(v > 0);
                        s += &dec(Some(SIGN[usize::from(v >= 0)]), (v - subtract).abs());
                    } else if v != 1 {
                        s += &dec(Some(H[j - 2]), v);
                    } else {
                        s += H[j - 2];
                    }
                }
            }
            print(strbuf, &s);
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeIsoTautString`
pub fn make_iso_taut_string(
    t: &[InChIIsotopicTGroup],
    n: i32,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let used0 = used(strbuf);
    let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    const LETTER: [char; 3] = ['t', 'd', 'h'];
    const H: [&str; 3] = ["T", "D", "H"];
    if ovfl == 0 {
        for i in 0..n.max(0) as usize {
            let g = &t[i];
            let mut s = String::new();
            for j in 0..4 {
                let v = match j {
                    0 => i32::from(g.n_t_group_number),
                    1 => i32::from(g.n_num_t),
                    2 => i32::from(g.n_num_d),
                    _ => i32::from(g.n_num_h),
                };
                if j == 0 {
                    s += &if abc_numbers {
                        abc(None, v)
                    } else {
                        dec(Some(if i != 0 { ITEM_DELIMETER } else { EXTRA_SPACE }), v)
                    };
                } else if v != 0 {
                    if abc_numbers {
                        s += &dec(None, v);
                        s.push(LETTER[j - 1]);
                    } else if v != 1 {
                        s += &dec(Some(H[j - 1]), v);
                    } else {
                        s += H[j - 1];
                    }
                }
            }
            print(strbuf, &s);
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeIsoHString`
pub fn make_iso_h_string(
    num_iso_h: &[i32; NUM_H_ISOTOPES],
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let used0 = used(strbuf);
    let abc_numbers = ct_mode & CT_MODE_ABC_NUMBERS != 0;
    const LETTER: [char; 3] = ['t', 'd', 'h'];
    const H: [&str; 3] = ["T", "D", "H"];
    if ovfl == 0 {
        let mut s = String::new();
        for j in 1..4 {
            let v = num_iso_h[NUM_H_ISOTOPES - j];
            if v == 0 {
                continue;
            }
            if abc_numbers {
                s += &dec(None, v);
                s.push(LETTER[j - 1]);
            } else if v != 1 {
                s += &dec(Some(H[j - 1]), v);
            } else {
                s += H[j - 1];
            }
        }
        print(strbuf, &s);
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `MakeStereoString`
pub fn make_stereo_string(
    at1: Option<&[u16]>,
    at2: Option<&[u16]>,
    parity: Option<&[i8]>,
    n_len_ct: i32,
    strbuf: &mut InchiIosString,
    ct_mode: i32,
    overflow: &mut i32,
) -> i32 {
    let ovfl = *overflow;
    let used0 = used(strbuf);
    const PARITY_CHAR: &[u8; 5] = b"!-+u?";
    if ovfl == 0 {
        for i in 0..n_len_ct.max(0) as usize {
            let mut s = String::new();
            let mut tot_len = 0usize;
            for j in 0..3 {
                let v = match j {
                    0 if at1.is_some() => i32::from(at1.unwrap()[i]),
                    1 if at2.is_some() => i32::from(at2.unwrap()[i]),
                    2 if parity.is_some() => i32::from(parity.unwrap()[i]),
                    _ => continue,
                };
                let piece = if ct_mode & CT_MODE_ABC_NUMBERS != 0 {
                    if j == 2 { dec(None, v) } else { abc(None, v) }
                } else if j < 2 {
                    dec(
                        if tot_len != 0 {
                            Some("-")
                        } else if i != 0 {
                            Some(ITEM_DELIMETER)
                        } else {
                            None
                        },
                        v,
                    )
                } else {
                    let c = if (0..=4).contains(&v) {
                        PARITY_CHAR[v as usize]
                    } else {
                        PARITY_CHAR[0]
                    };
                    (c as char).to_string()
                };
                tot_len += piece.len();
                s += &piece;
            }
            print(strbuf, &s);
        }
    }
    *overflow |= ovfl;
    used(strbuf) - used0
}

/// `print_sequence_of_nums_compressing_ranges`
pub fn print_sequence_of_nums_compressing_ranges(num: &[i32], strbuf: &mut InchiIosString) {
    let n = num.len();
    let mut range = 0;
    for k in 0..n.saturating_sub(1) {
        if num[k + 1] == num[k] + 1 {
            if range != 0 {
                range += 1;
            } else {
                print(strbuf, &format!("{}-", num[k]));
                range = 1;
            }
        } else {
            range = 0;
            print(strbuf, &format!("{},", num[k]));
        }
    }
    print(strbuf, &format!("{}", num[n - 1]));
}
