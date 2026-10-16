//! Atom symbols and the generic atom type lists ("alk", "Q", ...) of the
//! structure checker.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/AtomSymbolMatch.cpp`
//! (2026.09.1).
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

/// Periodic table of the checker, indexed by atomic number (`*` is 0).
const ATOM_SYMBOLS: [&str; 110] = [
    "*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hn", "Mt",
];

/// `getAtomicNumber`: `None` for a symbol that is not in the table (upstream
/// dereferences the end iterator there).
pub fn get_atomic_number(symbol: &str) -> Option<u32> {
    ATOM_SYMBOLS
        .iter()
        .rposition(|s| *s == symbol)
        .map(|n| n as u32)
}

const HC_TABLE: &[&str] = &["H", "C"];
const NON_METAL_HETERO_ELEMENTS: &[&str] = &[
    "He", "B", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Sb", "Te",
    "I", "Xe", "At",
];
const METALS: &[&str] = &[
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu",
];
const NON_METAL_SMALL_SOLUTION: &[&str] = &[
    "H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I",
];
const ALKALI_METALS: &[&str] = &["Li", "Na", "K", "Rb", "Cs", "Fr"];
const GR2: &[&str] = &["Be", "Mg", "Ca", "Sr", "Ba", "Ra"];
const GR3: &[&str] = &["B", "Al", "Ga", "In", "Tl"];
const GR4: &[&str] = &["C", "Si", "Ge", "Sn", "Pb"];
const ONS_TABLE: &[&str] = &["O", "N", "S"];
const ON2: &[&str] = &["O", "N", "S", "P", "Se", "Te", "Po"];
const HALOGENES: &[&str] = &["F", "Cl", "Br", "I", "At"];
const HA2: &[&str] = &["Cl", "Br", "I", "At"];
const TRANSITION_METALS: &[&str] = &[
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
];
const TRA: &[&str] = &[
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt",
];
const TRB: &[&str] = &["Cu", "Zn", "Ag", "Cd", "Au", "Hg"];
const TM1: &[&str] = &["Cu", "Ag", "Au"];
const TM2: &[&str] = &["Zn", "Cd", "Hg"];
const TM3: &[&str] = &["Sc", "Y", "La"];
const TM4: &[&str] = &["Ti", "Zr", "Hf"];
const TM5: &[&str] = &["V", "Nb", "Ta"];
const TM6: &[&str] = &["Cr", "Mo", "W"];
const TM7: &[&str] = &["Mn", "Tc", "Re"];
const TM8: &[&str] = &["Fe", "Co", "Ni", "Ru", "Rh", "Pd", "Os", "Ir", "Pt"];
const LANTHANOIDS: &[&str] = &[
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
];
const AMINO_ACIDS: &[&str] = &[
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
    "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
];

fn is_in_table(symbol: &str, table: &[&str]) -> bool {
    table.contains(&symbol)
}

/// The lower case pseudosymbols that may appear in a comma separated list.
fn pseudo_table(token: &str) -> Option<&'static [&'static str]> {
    Some(match token {
        "alk" => ALKALI_METALS,
        "gr2" => GR2,
        "gr3" => GR3,
        "gr4" => GR4,
        "ons" => ONS_TABLE,
        "on2" => ON2,
        "hal" => HALOGENES,
        "ha2" => HA2,
        "trn" => TRANSITION_METALS,
        "tra" => TRA,
        "trb" => TRB,
        "tm1" => TM1,
        "tm2" => TM2,
        "tm3" => TM3,
        "tm4" => TM4,
        "tm5" => TM5,
        "tm6" => TM6,
        "tm7" => TM7,
        "tm8" => TM8,
        "lan" => LANTHANOIDS,
        "ami" => AMINO_ACIDS,
        _ => return None,
    })
}

/// `AtomSymbolMatch`: `symbol` is in the comma separated list `pattern`
/// (which can contain the generic atom types like "alk" or "Q").
pub fn atom_symbol_match(symbol: &str, pattern: &str) -> bool {
    for token in pattern.split(',').filter(|t| !t.is_empty()) {
        if token.as_bytes()[0].is_ascii_lowercase()
            && let Some(table) = pseudo_table(token)
            && is_in_table(symbol, table)
        {
            return true;
        }
        if symbol == token {
            return true;
        }
    }
    match pattern {
        "A" => symbol != "H",
        "Qs" => is_in_table(symbol, NON_METAL_SMALL_SOLUTION),
        "G" => is_in_table(symbol, HC_TABLE),
        "ONS" => is_in_table(symbol, ONS_TABLE),
        "X" => is_in_table(symbol, HALOGENES),
        "Q" => is_in_table(symbol, NON_METAL_HETERO_ELEMENTS),
        "M" => is_in_table(symbol, METALS),
        "Ami" => is_in_table(symbol, AMINO_ACIDS),
        _ => false,
    }
}
