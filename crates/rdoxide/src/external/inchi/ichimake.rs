//! Construction of the `INChI` output structures (`INCHI_BASE/src/ichimake.c`):
//! the input-atom to `sp_ATOM` copy and the Hill formula comparison.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::ichister::SpAtom;
use super::inpdef::InpAtom;
use super::util::get_periodic_table_number;
use std::cmp::Ordering;

fn nul_str(b: &[u8]) -> &str {
    let n = b.iter().position(|&c| c == 0).unwrap_or(b.len());
    std::str::from_utf8(&b[..n]).unwrap_or("")
}

/// `inp2spATOM`: copies the structure fields of the input atoms into `sp_ATOM`
/// records (everything else is zeroed).
pub fn inp2sp_atom(inp_at: &[InpAtom], at: &mut [SpAtom]) {
    for (i, ia) in inp_at.iter().enumerate() {
        let mut a = SpAtom {
            elname: ia.elname,
            ..SpAtom::default()
        };
        a.el_number = get_periodic_table_number(nul_str(&a.elname)) as u8;
        a.valence = ia.valence;
        let val = ia.valence.max(0) as usize;
        for j in 0..val {
            a.neighbor[j] = ia.neighbor[j];
            a.bond_type[j] = ia.bond_type[j];
        }
        a.chem_bonds_valence = ia.chem_bonds_valence;
        a.orig_at_number = ia.orig_at_number;
        a.orig_compt_at_numb = ia.orig_compt_at_numb;
        a.endpoint = ia.endpoint;
        a.iso_atw_diff = ia.iso_atw_diff;
        a.num_h = ia.num_h;
        a.c_flags = ia.c_flags;
        a.num_iso_h = ia.num_iso_h;
        a.charge = ia.charge;
        a.radical = ia.radical;
        a.n_block_system = ia.n_block_system;
        a.b_cut_vertex = ia.b_cut_vertex;
        a.n_ring_system = ia.n_ring_system;
        a.n_num_at_in_ring_system = ia.n_num_at_in_ring_system;
        at[i] = a;
    }
}

/// `GetElementAndCount`: reads the next element and its count from a Hill
/// formula; returns `(ret, element, count)` with `ret` = 1 (read), 0 (end of
/// the formula: "Zzz", 99999) or -1 (not a formula). `f` is advanced.
pub fn get_element_and_count(f: &mut &[u8]) -> (i32, String, i32) {
    let p = *f;
    if let Some(&c0) = p.first().filter(|&&c| c != 0) {
        if c0.is_ascii_uppercase() {
            let mut el = vec![c0];
            let mut k = 1;
            if p.get(1).is_some_and(|c| c.is_ascii_lowercase()) {
                el.push(p[1]);
                k = 2;
            }
            if el == b"C" {
                el[0] = b'A'; // less than any element: carbon compounds first
            }
            let mut q = k;
            let mut count = 1;
            if p.get(q).is_some_and(|c| c.is_ascii_digit()) {
                let mut v: i64 = 0;
                while let Some(d) = p.get(q).filter(|c| c.is_ascii_digit()) {
                    v = (v * 10 + i64::from(d - b'0')).min(i64::from(i32::MAX));
                    q += 1;
                }
                count = v as i32;
            }
            *f = &p[q..];
            return (1, String::from_utf8(el).unwrap_or_default(), count);
        }
        return (-1, String::new(), 0);
    }
    (0, "Zzz".to_string(), 99999)
}

/// `CompareHillFormulas`
pub fn compare_hill_formulas(f1: &str, f2: &str) -> i32 {
    let (mut f1, mut f2) = (f1.as_bytes(), f2.as_bytes());
    loop {
        let (ret1, el1, count1) = get_element_and_count(&mut f1);
        let (ret2, el2, count2) = get_element_and_count(&mut f2);
        if ret1 >= 0 && ret2 >= 0 {
            match el1.cmp(&el2) {
                Ordering::Less => return -1,
                Ordering::Greater => return 1,
                Ordering::Equal => {}
            }
            if count2 != count1 {
                return count2 - count1;
            }
        } else {
            return 0;
        }
        if !(ret1 > 0 && ret2 > 0) {
            return 0;
        }
    }
}

/// `CompareHillFormulasNoH`: like `compare_hill_formulas` ignoring the `H`
/// element; the numbers of H are added to `num_h1`/`num_h2`.
pub fn compare_hill_formulas_no_h(f1: &str, f2: &str, num_h1: &mut i32, num_h2: &mut i32) -> i32 {
    let (mut f1, mut f2) = (f1.as_bytes(), f2.as_bytes());
    loop {
        let (mut ret1, mut el1, mut count1) = get_element_and_count(&mut f1);
        if ret1 > 0 && el1 == "H" {
            *num_h1 += count1;
            (ret1, el1, count1) = get_element_and_count(&mut f1);
        }
        let (mut ret2, mut el2, mut count2) = get_element_and_count(&mut f2);
        if ret2 > 0 && el2 == "H" {
            *num_h2 += count2;
            (ret2, el2, count2) = get_element_and_count(&mut f2);
        }
        if ret1 >= 0 && ret2 >= 0 {
            match el1.cmp(&el2) {
                Ordering::Less => return -1,
                Ordering::Greater => return 1,
                Ordering::Equal => {}
            }
            if count2 != count1 {
                return count2 - count1;
            }
        } else {
            return 0;
        }
        if !(ret1 > 0 && ret2 > 0) {
            return 0;
        }
    }
}

/// `CT_MODE_PREDECESSORS`
pub const CT_MODE_PREDECESSORS: i32 = 8;

/// `mystrrev`: reverses a string in place (up to the first NUL)
pub fn mystrrev(p: &mut [u8]) {
    let n = p.iter().position(|&c| c == 0).unwrap_or(p.len());
    p[..n].reverse();
}

/// `GetDfsOrder4CT`: the DFS order of the connection table output as triplets
/// (canonical number, 16 + number of H, delimiter), zero terminated; `None`
/// for an inconsistent table or when the output buffer would overflow.
pub fn get_dfs_order_4ct(
    linear_ct: &[u16],
    n_len_ct: usize,
    n_num_h: Option<&[i8]>,
    num_atoms: usize,
    n_ct_mode: i32,
) -> Option<Vec<u16>> {
    use super::ichimap::create_neigh_list_from_linear_ct;
    use super::ichisort::insertions_sort;
    use super::inpdef::MAX_ATOMS;
    const DELETED: u16 = MAX_ATOMS as u16 + 1;
    let ct_pred = n_ct_mode & CT_MODE_PREDECESSORS != 0;
    let mut nl = create_neigh_list_from_linear_ct(linear_ct, n_len_ct, num_atoms)?;
    let mut stack = vec![0u16; num_atoms];
    let mut num_desc = vec![0u16; num_atoms];
    let mut dfs_number = vec![0u16; num_atoms];
    let mut neigh_numb = vec![0u8; num_atoms];
    let mut start = 0usize;
    if !ct_pred {
        for i in 1..num_atoms {
            if nl[i][0] < nl[start][0] {
                start = i;
            }
        }
    }
    let mut u = start;
    let mut n_dfs: u16 = 0;
    let mut top: i32 = -1;
    n_dfs += 1;
    dfs_number[u] = n_dfs;
    num_desc[u] = u16::from(!ct_pred);
    top += 1;
    stack[top as usize] = u as u16;
    let mut num_rings = 0usize;
    loop {
        let mut i;
        loop {
            i = stack[top as usize] as usize;
            let j = usize::from(neigh_numb[i]) + 1;
            if (nl[i][0] as usize) < j {
                break;
            }
            neigh_numb[i] += 1;
            u = nl[i][j] as usize;
            if dfs_number[u] == 0 {
                top += 1;
                stack[top as usize] = u as u16;
                n_dfs += 1;
                dfs_number[u] = n_dfs;
                if ct_pred {
                    num_desc[u] = (i + 1) as u16;
                } else {
                    num_desc[u] = num_desc[u].wrapping_add(1);
                }
            } else if top != 0
                && u != stack[top as usize - 1] as usize
                && dfs_number[u] < dfs_number[i]
            {
                if !ct_pred {
                    num_desc[i] = num_desc[i].wrapping_add(1);
                }
                num_rings += 1;
            } else {
                nl[i][j] = DELETED;
            }
        }
        neigh_numb[i] = 0;
        if !ct_pred && top != 0 {
            let p = stack[top as usize - 1] as usize;
            num_desc[p] = num_desc[p].wrapping_add(num_desc[i]);
        }
        top -= 1;
        if top < 0 {
            break;
        }
    }
    for i in 0..num_atoms {
        if nl[i][0] > 1 {
            let cur = i;
            let cnt = nl[i][0] as usize;
            let (dn, nd) = (&dfs_number, &num_desc);
            insertions_sort(&mut nl[i][1..=cnt], |a, b| {
                let (n1, n2) = (i32::from(*a), i32::from(*b));
                if n1 > MAX_ATOMS as i32 {
                    return i32::from(n2 <= MAX_ATOMS as i32);
                }
                if n2 > MAX_ATOMS as i32 {
                    return -1;
                }
                let cur_dfs = dn[cur];
                let d1 = if cur_dfs > dn[n1 as usize] {
                    0
                } else {
                    i32::from(nd[n1 as usize])
                };
                let d2 = if cur_dfs > dn[n2 as usize] {
                    0
                } else {
                    i32::from(nd[n2 as usize])
                };
                if d1 != d2 {
                    return d1 - d2;
                }
                n1 - n2
            });
        }
        let cnt = nl[i][0] as usize;
        let mut k = 0;
        while k < cnt && nl[i][k + 1] <= MAX_ATOMS as u16 {
            k += 1;
        }
        nl[i][0] = k as u16;
    }
    let tot_len = 3 * (num_atoms + num_rings + 1);
    let mut out = vec![0u16; tot_len];
    let h16 = |u: usize| n_num_h.map_or(0, |h| (16 + i32::from(h[u])) as u16);
    if ct_pred {
        let mut k: i64 = -3;
        for u in 0..num_atoms {
            k += 3;
            if k as usize + 6 > tot_len {
                return None;
            }
            let ku = k as usize;
            out[ku] = if num_desc[u] != 0 {
                num_desc[u]
            } else {
                DELETED
            };
            out[ku + 1] = h16(u);
            out[ku + 2] = if k != 0 { u16::from(b',') } else { 0 };
            let mut j = 1;
            while j <= nl[u][0] as usize && dfs_number[u] > dfs_number[nl[u][j] as usize] {
                let i = nl[u][j] as usize;
                k += 3;
                if k as usize + 6 > tot_len {
                    return None;
                }
                let ku = k as usize;
                out[ku] = (i + 1) as u16;
                out[ku + 1] = 0;
                out[ku + 2] = u16::from(b'-');
                j += 1;
            }
        }
    } else {
        let mut u = start;
        top = 0;
        neigh_numb.iter_mut().for_each(|x| *x = 0);
        stack[0] = u as u16;
        let mut k = 0usize;
        out[k] = (u + 1) as u16;
        out[k + 1] = h16(u);
        out[k + 2] = 0;
        loop {
            let mut i;
            loop {
                i = stack[top as usize] as usize;
                let j = usize::from(neigh_numb[i]) + 1;
                if !(i < num_atoms && nl[i][0] as usize >= j) {
                    break;
                }
                k += 3;
                if k + 6 > tot_len {
                    return None;
                }
                neigh_numb[i] += 1;
                u = nl[i][j] as usize;
                out[k] = (u + 1) as u16;
                if dfs_number[u] > dfs_number[i] {
                    top += 1;
                    stack[top as usize] = u as u16;
                    out[k + 1] = h16(u);
                } else {
                    out[k + 1] = 0;
                }
                let n = nl[i][0] as usize;
                let delim = if n > 1 {
                    if j == 1 {
                        b'('
                    } else if j == n {
                        b')'
                    } else {
                        b','
                    }
                } else {
                    b'-'
                };
                out[k + 2] = u16::from(delim);
            }
            neigh_numb[i] = 0;
            top -= 1;
            if top < 0 {
                break;
            }
        }
    }
    Some(out)
}
