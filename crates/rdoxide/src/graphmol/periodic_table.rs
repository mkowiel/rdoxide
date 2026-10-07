//! Ported from RDKit `Code/GraphMol/PeriodicTable.{h,cpp}` and
//! `Code/GraphMol/atomic_data.{h,cpp}`.
//!
//! Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! The element data lives in `data/*.tsv`, extracted verbatim from upstream's
//! `atomic_data.cpp` by `tools/sync_upstream.py regen`.

use std::collections::{BTreeMap, HashMap};
use std::sync::OnceLock;

use crate::error::{RDKitError, Result};

static ATOMIC_DATA: &str = include_str!("../../data/atomic_data.tsv");
static ISOTOPE_DATA: &str = include_str!("../../data/isotopes.tsv");
static ELEMENT_NAMES: &str = include_str!("../../data/element_names.txt");

/// Mass of an electron in amu.
pub const ELECTRON_MASS: f64 = 0.00054857991;

/// Per-element data (`RDKit::atomicData`).
#[derive(Debug, Clone)]
pub struct AtomicData {
    pub anum: u32,
    pub symbol: String,
    pub name: String,
    pub row: u32,
    pub r_cov: f64,
    pub r_b0: f64,
    pub r_vdw: f64,
    pub mass: f64,
    pub n_outer_elecs: i32,
    pub common_isotope: i32,
    pub common_isotope_mass: f64,
    /// First entry is the default valence; a trailing `-1` means any higher
    /// valence is allowed.
    pub valences: Vec<i32>,
    /// isotope -> (mass, abundance)
    pub isotopes: BTreeMap<u32, (f64, f64)>,
}

/// The periodic table singleton (`RDKit::PeriodicTable`).
#[derive(Debug)]
pub struct PeriodicTable {
    by_anum: Vec<AtomicData>,
    by_name: HashMap<String, u32>,
}

fn data_lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
}

impl PeriodicTable {
    /// `PeriodicTable::getTable()`
    pub fn get() -> &'static PeriodicTable {
        static TABLE: OnceLock<PeriodicTable> = OnceLock::new();
        TABLE.get_or_init(PeriodicTable::build)
    }

    fn build() -> PeriodicTable {
        let names: Vec<&str> = data_lines(ELEMENT_NAMES).collect();
        let mut by_anum: Vec<AtomicData> = Vec::new();
        let mut by_name = HashMap::new();
        for line in data_lines(ATOMIC_DATA) {
            let t: Vec<&str> = line.split_whitespace().collect();
            let num = |i: usize| -> f64 { t[i].parse().expect("numeric atomic data") };
            let anum: u32 = t[0].parse().expect("atomic number");
            let data = AtomicData {
                anum,
                symbol: t[1].to_string(),
                name: names[anum as usize].to_string(),
                row: num(2) as u32,
                r_cov: num(3),
                r_b0: num(4),
                r_vdw: num(5),
                mass: num(6),
                n_outer_elecs: num(7) as i32,
                common_isotope: num(8) as i32,
                common_isotope_mass: num(9),
                valences: t[10..]
                    .iter()
                    .map(|v| v.parse().expect("valence"))
                    .collect(),
                isotopes: BTreeMap::new(),
            };
            by_name.insert(data.symbol.clone(), anum);
            // duplicate rows for backwards-compatible symbols (#2784): keep
            // only the first entry per atomic number
            if anum as usize == by_anum.len() {
                by_anum.push(data);
            }
        }
        for line in data_lines(ISOTOPE_DATA) {
            let t: Vec<&str> = line.split_whitespace().collect();
            if t.len() < 5 {
                continue;
            }
            let anum: usize = t[0].parse().expect("atomic number");
            let iso: u32 = t[2].parse().expect("isotope");
            let mass: f64 = t[3].parse().expect("isotope mass");
            let abundance: f64 = t[4].parse().expect("abundance");
            by_anum[anum].isotopes.insert(iso, (mass, abundance));
        }
        PeriodicTable { by_anum, by_name }
    }

    fn data(&self, anum: u32) -> Result<&AtomicData> {
        self.by_anum
            .get(anum as usize)
            .ok_or_else(|| RDKitError::precondition("Atomic number not found"))
    }

    /// Infallible lookup for internal callers that already validated `anum`.
    pub fn element(&self, anum: u32) -> &AtomicData {
        &self.by_anum[anum as usize]
    }

    fn anum_for_symbol(&self, symbol: &str) -> Result<u32> {
        self.by_name
            .get(symbol)
            .copied()
            .ok_or_else(|| RDKitError::precondition(format!("Element '{symbol}' not found")))
    }

    /// `getAtomicNumber`
    pub fn get_atomic_number(&self, symbol: &str) -> Result<u32> {
        match symbol {
            "C" => Ok(6),
            "N" => Ok(7),
            "O" => Ok(8),
            _ => {
                self.by_name.get(symbol).copied().ok_or_else(|| {
                    RDKitError::postcondition(format!("Element '{symbol}' not found"))
                })
            }
        }
    }

    /// Lookup accepting either an atomic number or a symbol.
    pub fn resolve(&self, elem: Element<'_>) -> Result<&AtomicData> {
        match elem {
            Element::Num(n) => self.data(n),
            Element::Symbol(s) => self.data(self.anum_for_symbol(s)?),
        }
    }

    pub fn get_atomic_weight(&self, anum: u32) -> Result<f64> {
        Ok(self.data(anum)?.mass)
    }
    pub fn get_element_symbol(&self, anum: u32) -> Result<&str> {
        Ok(&self.data(anum)?.symbol)
    }
    pub fn get_element_name(&self, anum: u32) -> Result<&str> {
        Ok(&self.data(anum)?.name)
    }
    pub fn get_rvdw(&self, anum: u32) -> Result<f64> {
        Ok(self.data(anum)?.r_vdw)
    }
    pub fn get_rcovalent(&self, anum: u32) -> Result<f64> {
        Ok(self.data(anum)?.r_cov)
    }
    pub fn get_rb0(&self, anum: u32) -> Result<f64> {
        Ok(self.data(anum)?.r_b0)
    }
    pub fn get_default_valence(&self, anum: u32) -> Result<i32> {
        Ok(self.data(anum)?.valences[0])
    }
    pub fn get_valence_list(&self, anum: u32) -> Result<&[i32]> {
        Ok(&self.data(anum)?.valences)
    }
    pub fn get_n_outer_elecs(&self, anum: u32) -> Result<i32> {
        Ok(self.data(anum)?.n_outer_elecs)
    }
    pub fn get_most_common_isotope(&self, anum: u32) -> Result<i32> {
        Ok(self.data(anum)?.common_isotope)
    }
    pub fn get_most_common_isotope_mass(&self, anum: u32) -> Result<f64> {
        Ok(self.data(anum)?.common_isotope_mass)
    }
    pub fn get_row(&self, anum: u32) -> Result<u32> {
        Ok(self.data(anum)?.row)
    }
    /// Returns 0.0 for unknown isotopes, like upstream.
    pub fn get_mass_for_isotope(&self, anum: u32, isotope: u32) -> Result<f64> {
        Ok(self.data(anum)?.isotopes.get(&isotope).map_or(0.0, |x| x.0))
    }
    /// Returns 0.0 for unknown isotopes, like upstream.
    pub fn get_abundance_for_isotope(&self, anum: u32, isotope: u32) -> Result<f64> {
        Ok(self.data(anum)?.isotopes.get(&isotope).map_or(0.0, |x| x.1))
    }
    pub fn get_max_atomic_number(&self) -> u32 {
        (self.by_anum.len() - 1) as u32
    }

    /// `moreElectroNegative`: crude ordering by number of outer-shell
    /// electrons, then by (lower) atomic number.
    pub fn more_electronegative(&self, anum1: u32, anum2: u32) -> Result<bool> {
        let ne1 = self.data(anum1)?.n_outer_elecs;
        let ne2 = self.data(anum2)?.n_outer_elecs;
        Ok(ne1 > ne2 || (ne1 == ne2 && anum1 < anum2))
    }
}

/// An element given either by atomic number or by symbol.
#[derive(Debug, Clone, Copy)]
pub enum Element<'a> {
    Num(u32),
    Symbol(&'a str),
}

#[cfg(test)]
mod tests {
    //! Ported from Code/GraphMol/catch_periodictable.cpp and test1.cpp.
    use super::*;

    #[test]
    fn basics() {
        let pt = PeriodicTable::get();
        assert_eq!(pt.get_atomic_number("C").unwrap(), 6);
        assert_eq!(pt.get_atomic_number("Cl").unwrap(), 17);
        assert_eq!(pt.get_element_symbol(8).unwrap(), "O");
        assert_eq!(pt.get_element_name(8).unwrap(), "Oxygen");
        assert_eq!(pt.get_element_symbol(0).unwrap(), "*");
        assert!((pt.get_atomic_weight(6).unwrap() - 12.011).abs() < 1e-6);
        assert_eq!(pt.get_default_valence(7).unwrap(), 3);
        assert_eq!(pt.get_valence_list(16).unwrap(), &[2, 4, 6]);
        assert_eq!(pt.get_n_outer_elecs(6).unwrap(), 4);
        assert_eq!(pt.get_max_atomic_number(), 118);
        assert!(pt.get_atomic_number("Xx").is_err());
        assert!(pt.get_atomic_weight(200).is_err());
    }

    #[test]
    fn legacy_symbols() {
        // #2784: duplicated rows keep the first symbol per atomic number
        let pt = PeriodicTable::get();
        assert_eq!(pt.get_atomic_number("Uut").unwrap(), 113);
        assert_eq!(pt.get_atomic_number("Nh").unwrap(), 113);
        assert_eq!(pt.get_element_symbol(113).unwrap(), "Nh");
        assert_eq!(pt.get_element_symbol(115).unwrap(), "Mc");
    }

    #[test]
    fn isotopes() {
        let pt = PeriodicTable::get();
        assert!((pt.get_mass_for_isotope(6, 13).unwrap() - 13.003354838).abs() < 1e-8);
        assert!((pt.get_abundance_for_isotope(1, 2).unwrap() - 0.0115).abs() < 1e-8);
        assert_eq!(pt.get_mass_for_isotope(6, 1000).unwrap(), 0.0);
        assert_eq!(pt.get_most_common_isotope(6).unwrap(), 12);
    }

    #[test]
    fn electronegativity() {
        let pt = PeriodicTable::get();
        assert!(pt.more_electronegative(8, 6).unwrap());
        assert!(!pt.more_electronegative(6, 8).unwrap());
        assert!(pt.more_electronegative(9, 17).unwrap());
    }
}
