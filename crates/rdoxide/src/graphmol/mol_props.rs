//! Molecular weights and formula (upstream `GraphMol/MolProps.cpp`,
//! `MolOps::getAvgMolWt`, `getExactMolWt` and `getMolFormula`).

use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::error::Result;
use crate::graphmol::periodic_table::{ELECTRON_MASS, PeriodicTable};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::fp::fmuladd;

/// `getAvgMolWt`: average molecular weight, implicit Hs included unless
/// `only_heavy` is set.
pub fn get_avg_mol_wt(mol: &ROMol, only_heavy: bool) -> Result<f64> {
    let table = PeriodicTable::get();
    let h_weight = table.get_atomic_weight(1)?;
    let mut res = 0.0;
    for (i, atom) in mol.atoms().iter().enumerate() {
        if !only_heavy || atom.atomic_num() != 1 {
            res += atom.mass();
        }
        if !only_heavy {
            let n_hs = f64::from(mol.total_num_hs(i as u32, false)?);
            res = fmuladd(n_hs, h_weight, res);
        }
    }
    Ok(res)
}

/// `getExactMolWt`: monoisotopic molecular weight with the electron mass
/// corrected for the formal charges.
pub fn get_exact_mol_wt(mol: &ROMol, only_heavy: bool) -> Result<f64> {
    let table = PeriodicTable::get();
    let mut res = 0.0;
    let mut n_hs_to_count = 0u32;
    for (i, atom) in mol.atoms().iter().enumerate() {
        let atomic_num = atom.atomic_num();
        if atomic_num != 1 || !only_heavy {
            res += if atom.isotope() == 0 {
                table.get_most_common_isotope_mass(atomic_num)?
            } else {
                atom.mass()
            };
            res = fmuladd(-ELECTRON_MASS, f64::from(atom.formal_charge()), res);
        }
        if !only_heavy {
            n_hs_to_count += mol.total_num_hs(i as u32, false)?;
        }
    }
    if !only_heavy {
        res = fmuladd(
            f64::from(n_hs_to_count),
            table.get_most_common_isotope_mass(1)?,
            res,
        );
    }
    Ok(res)
}

/// A formula term: the isotope (0 for the natural mixture) and the symbol.
type FormulaKey = (u32, String);

/// `HillCompare` as a three-way ordering: C first, then H, then D and T, then
/// everything else by isotope and symbol. The keys it sees are distinct, so
/// the order is total and independent of the sorting algorithm.
fn hill_order(a: &FormulaKey, b: &FormulaKey) -> Ordering {
    // special cases: C, H, D and T go at the beginning
    for special in ["C", "H", "D", "T"] {
        match (a.1 == special, b.1 == special) {
            (true, true) => return a.0.cmp(&b.0),
            (true, false) => return Ordering::Less,
            (false, true) => return Ordering::Greater,
            (false, false) => {}
        }
    }
    a.cmp(b)
}

/// `getMolFormula`: Hill-ordered formula with the charge appended.
pub fn get_mol_formula(
    mol: &ROMol,
    separate_isotopes: bool,
    abbreviate_h_isotopes: bool,
) -> Result<String> {
    let table = PeriodicTable::get();
    let mut counts: BTreeMap<FormulaKey, u32> = BTreeMap::new();
    let mut charge = 0i32;
    let mut n_hs = 0u32;
    for (i, atom) in mol.atoms().iter().enumerate() {
        let atomic_num = atom.atomic_num();
        let mut key: FormulaKey = (0, table.get_element_symbol(atomic_num)?.to_string());
        if separate_isotopes {
            let isotope = atom.isotope();
            if abbreviate_h_isotopes && atomic_num == 1 && (isotope == 2 || isotope == 3) {
                key.1 = if isotope == 2 { "D" } else { "T" }.to_string();
            } else {
                key.0 = isotope;
            }
        }
        *counts.entry(key).or_insert(0) += 1;
        n_hs += mol.total_num_hs(i as u32, false)?;
        charge += atom.formal_charge();
    }
    if n_hs > 0 {
        *counts.entry((0, "H".to_string())).or_insert(0) += n_hs;
    }

    let mut keys: Vec<&FormulaKey> = counts.keys().collect();
    keys.sort_by(|a, b| hill_order(a, b));

    let mut res = String::new();
    for key in keys {
        if key.0 > 0 {
            res.push_str(&format!("[{}{}]", key.0, key.1));
        } else {
            res.push_str(&key.1);
        }
        if counts[key] > 1 {
            res.push_str(&counts[key].to_string());
        }
    }
    match charge.cmp(&0) {
        Ordering::Greater => {
            res.push('+');
            if charge > 1 {
                res.push_str(&charge.to_string());
            }
        }
        Ordering::Less => {
            res.push('-');
            if charge < -1 {
                res.push_str(&(-charge).to_string());
            }
        }
        Ordering::Equal => {}
    }
    Ok(res)
}
