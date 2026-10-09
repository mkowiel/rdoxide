//! The Hill-order molecular formula of `MolHash` (upstream `mf.h`).

use std::sync::OnceLock;

use crate::error::Result;
use crate::graphmol::romol::ROMol;

/// Element symbols by atomic number; 0 (a dummy) is written `X`.
const SYMBOLS: [&str; 119] = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// Atomic numbers in alphabetical order of their symbols; with `organic`,
/// carbon and hydrogen come first.
fn hill_order(organic: bool) -> &'static [usize] {
    static ORDERS: OnceLock<[Vec<usize>; 2]> = OnceLock::new();
    let orders = ORDERS.get_or_init(|| {
        let alphabetical = |skip: &[usize]| {
            let mut v: Vec<usize> = (0..SYMBOLS.len()).filter(|i| !skip.contains(i)).collect();
            v.sort_by_key(|&i| SYMBOLS[i]);
            v
        };
        let mut organic = vec![6, 1];
        organic.extend(alphabetical(&[6, 1]));
        [alphabetical(&[]), organic]
    });
    &orders[usize::from(organic)]
}

/// `NMMolecularFormula(mol)`
pub(super) fn molecular_formula(mol: &ROMol) -> Result<String> {
    let mut hist = [0u32; 256];
    let mut charge = 0i32;
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom(idx);
        let elem = atom.atomic_num() as usize;
        hist[if elem < 256 { elem } else { 0 }] += 1;
        hist[1] += mol.total_num_hs(idx, false)?;
        charge += atom.formal_charge();
    }
    let mut result = String::new();
    for &elem in hill_order(hist[6] > 0) {
        if hist[elem] > 0 {
            result += SYMBOLS[elem];
            if hist[elem] > 1 {
                result += &hist[elem].to_string();
            }
        }
    }
    match charge {
        0 => {}
        1 => result.push('+'),
        -1 => result.push('-'),
        c if c > 0 => result += &format!("+{c}"),
        c => result += &format!("-{}", -c),
    }
    Ok(result)
}
