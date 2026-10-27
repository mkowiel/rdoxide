//! `utilities.c`, `casutils.c`: atom symbol matching, implicit hydrogens, CAS valences.
//! Copyright 2001-2011 Novartis Pharma AG.

use super::tables::*;

fn in_table(sym: &str, table: &[&str]) -> bool {
    table.contains(&sym)
}

/// `AtomSymbolMatch`: `atsym` is in the comma separated list `pattern`, whose entries may be
/// the pseudo symbols `alk`, `gr2`, `hal`, ... (lower case) or a plain symbol; the whole pattern
/// may also be `A`, `Qs`, `G`, `ONS`, `X`, `Q`, `M` or `Ami`
pub fn atom_symbol_match(atsym: &str, pattern: &str) -> bool {
    for tok in pattern.split(',').filter(|t| !t.is_empty()) {
        if tok.starts_with(|c: char| c.is_ascii_lowercase()) {
            let table: Option<&[&str]> = match tok {
                "alk" => Some(ALKALI_METALS),
                "gr2" => Some(GR2),
                "gr3" => Some(GR3),
                "gr4" => Some(GR4),
                "ons" => Some(ONS_TABLE),
                "on2" => Some(ON2),
                "hal" => Some(HALOGENES),
                "ha2" => Some(HA2),
                "trn" => Some(TRANSITION_METALS),
                "tra" => Some(TRA),
                "trb" => Some(TRB),
                "tm1" => Some(TM1),
                "tm2" => Some(TM2),
                "tm3" => Some(TM3),
                "tm4" => Some(TM4),
                "tm5" => Some(TM5),
                "tm6" => Some(TM6),
                "tm7" => Some(TM7),
                "tm8" => Some(TM8),
                "lan" => Some(LANTHANOIDS),
                "ami" => Some(AMINO_ACIDS),
                _ => None,
            };
            if table.is_some_and(|t| in_table(atsym, t)) {
                return true;
            }
        }
        if atsym == tok {
            return true;
        }
    }
    match pattern {
        "A" => atsym != "H",
        "Qs" => in_table(atsym, NON_METAL_SMALL_SOLUTION),
        "G" => in_table(atsym, HC_TABLE),
        "ONS" => in_table(atsym, ONS_TABLE),
        "X" => in_table(atsym, HALOGENES),
        "Q" => in_table(atsym, NON_METAL_HETERO_ELEMENTS),
        "M" => in_table(atsym, METALS),
        "Ami" => in_table(atsym, AMINO_ACIDS),
        _ => false,
    }
}

/// `CASValence`: 0 for unknown elements
pub fn cas_valence(symbol: &str) -> i32 {
    CAS_VALENCE
        .iter()
        .find(|(s, _)| *s == symbol)
        .map_or(0, |e| e.1)
}

/// `valence_table_entry`: symbol, from, to, step, lone pairs, pair deficit
const VALENCE_TABLE: &[(&str, i32, i32, i32, i32, i32)] = &[
    ("C", 4, 4, 1, 0, 0),
    ("H", -1, 1, 2, 0, 0),
    ("N", 3, 5, 2, 1, 0),
    ("O", 2, 2, 1, 2, 0),
    ("Cl", 1, 7, 2, 0, 0),
    ("P", 3, 5, 2, 1, 0),
    ("S", 2, 6, 2, 2, 0),
    ("F", 1, 1, 1, 1, 0),
    ("H", -1, 1, 2, 0, 0),
    ("Li", -1, 1, 2, 0, 0),
    ("Na", -1, 1, 2, 0, 0),
    ("K", -1, 1, 2, 0, 0),
    ("Rb", -1, 1, 2, 0, 0),
    ("Cs", -1, 1, 2, 0, 0),
    ("Be", -2, 2, 2, 0, 0),
    ("Mg", -2, 2, 2, 0, 0),
    ("B", 3, 3, 1, 0, 1),
    ("Al", -3, 3, 2, 0, 1),
    ("Ga", -3, 3, 2, 0, 1),
    ("In", -3, 3, 2, 0, 1),
    ("Tl", -3, 3, 2, 0, 0),
    ("Si", 4, 4, 1, 0, 0),
    ("As", 3, 5, 2, 0, 0),
    ("Sb", 3, 5, 2, 0, 0),
    ("Bi", 3, 5, 2, 0, 0),
    ("Se", 2, 6, 2, 0, 0),
    ("Te", 2, 6, 2, 0, 0),
    ("La", -3, 3, 2, 0, 1),
    ("Ce", -3, 3, 2, 0, 1),
    ("Pr", -3, 3, 2, 0, 1),
    ("Nd", -3, 3, 2, 0, 1),
    ("Pm", -3, 3, 2, 0, 1),
    ("Sm", -3, 3, 2, 0, 1),
    ("Eu", -3, 3, 2, 0, 1),
    ("Gd", -3, 3, 2, 0, 1),
    ("Tb", -3, 3, 2, 0, 1),
    ("Dy", -3, 3, 2, 0, 1),
    ("Ho", -3, 3, 2, 0, 1),
    ("Er", -3, 3, 2, 0, 1),
    ("Tm", -3, 3, 2, 0, 1),
    ("Yb", -3, 3, 2, 0, 1),
    ("Lu", -3, 3, 2, 0, 1),
    ("Br", 1, 7, 2, 0, 0),
    ("I", 1, 7, 2, 0, 0),
];

/// `ImplicitHydrogens`
pub fn implicit_hydrogens(
    symbol: &str,
    nsingle: i32,
    naromatic: i32,
    ndouble: i32,
    ntriple: i32,
    radical: bool,
    charge: i32,
) -> i32 {
    let mut bond_electrons = nsingle + 2 * ndouble + 3 * ntriple;
    if radical {
        bond_electrons += 1;
    }
    bond_electrons += match naromatic {
        0 => 0,
        1 => 2,
        2 => 3,
        3 => 4,
        n => n + 1,
    };
    for &(name, from, to, step, lone_pairs, pair_deficit) in VALENCE_TABLE {
        if name != symbol {
            continue;
        }
        let mut val = from;
        while val <= to {
            if charge == 0 {
                let h = val - bond_electrons;
                if h >= 0 {
                    return h;
                }
            } else if charge > 0 {
                let h = val - bond_electrons + charge;
                if h >= 0 {
                    return if lone_pairs > 0 { h } else { 0 };
                }
            } else {
                let mut h = val - bond_electrons - charge;
                if h >= 0 {
                    if pair_deficit < h {
                        h = pair_deficit;
                    }
                    return if lone_pairs > 0 { 0 } else { h };
                }
            }
            val += step;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbol_groups() {
        assert!(atom_symbol_match("Cl", "hal"));
        assert!(atom_symbol_match("N", "O,S,P,N,L"));
        assert!(atom_symbol_match("C", "A"));
        assert!(!atom_symbol_match("H", "A"));
        assert!(atom_symbol_match("O", "Q"));
        assert_eq!(cas_valence("C"), 4);
        assert_eq!(cas_valence("??x"), 0);
        assert_eq!(implicit_hydrogens("N", 3, 0, 0, 0, false, 0), 0);
        assert_eq!(implicit_hydrogens("O", 1, 0, 0, 0, false, 0), 1);
    }
}
