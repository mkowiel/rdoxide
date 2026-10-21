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
