//! Port of the SMILES reading helpers of Avalon Toolkit's `smi2mol.c`: the bracket atom
//! reader (`SmilesToMDLAtom`), the valence preprocessing for denormalization
//! (`ComputeSmilesValenceAndHCount`) and the syntax pre-check (`QuickCheck`).

use super::perceive::{DOUBLE, SINGLE, TRIPLE, setup_neighbourhood};
use super::reaccs::Molecule;
use super::tables::PTABLE;

/// the MDL attributes of a SMILES atom in square brackets
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmilesAtom {
    pub symbol: String,
    pub charge: i32,
    pub radical: i32,
    pub isotope: i32,
    pub hydrogen: i32,
    pub stereo: i32,
    pub mapping: i32,
    pub aromatic: bool,
}

/// `SmilesToMDLAtom`: reads the atom at `smiles[pos..]`, which must be at the `[`. Returns the
/// position of the terminal `]` (or where the scan stopped) and the atom; `symbol` is empty
/// on errors. If the text is not at a `[`, `pos` is returned unchanged with an empty atom.
pub fn smiles_to_mdl_atom(smiles: &[u8], pos: usize) -> (usize, SmilesAtom) {
    let at = |i: usize| smiles.get(i).copied().unwrap_or(0);
    let mut a = SmilesAtom::default();
    let mut p = pos;
    if at(p) == b'[' {
        p += 1;
    } else {
        return (p, a);
    }
    while at(p).is_ascii_digit() {
        a.isotope = a
            .isotope
            .wrapping_mul(10)
            .wrapping_add(i32::from(at(p) - b'0'));
        p += 1;
    }
    let mut symbol: Vec<u8> = Vec::new();
    if at(p).is_ascii_uppercase() || at(p) == b'*' {
        symbol.push(at(p));
    } else if at(p).is_ascii_lowercase() {
        a.aromatic = true;
        symbol.push(at(p) - b'a' + b'A');
    } else {
        return (p, a);
    }
    p += 1;
    while at(p).is_ascii_lowercase() {
        symbol.push(at(p));
        p += 1;
    }
    while at(p) == b'@' {
        a.stereo += 1;
        p += 1;
    }
    if a.stereo == 1 && at(p).is_ascii_digit() {
        a.stereo = 0;
        while at(p).is_ascii_digit() {
            a.stereo = a
                .stereo
                .wrapping_mul(10)
                .wrapping_add(i32::from(at(p) - b'0'));
            p += 1;
        }
    }
    let mut sym = String::from_utf8_lossy(&symbol).into_owned();
    if a.isotope > 0 && symbol[0].is_ascii_uppercase() {
        let Some(e) = PTABLE.iter().find(|e| e.0 == sym) else {
            return (p, a);
        };
        if sym == "H" {
            if a.isotope == 2 {
                sym = "D".into();
            } else if a.isotope == 3 {
                sym = "T".into();
            }
            a.isotope = 0;
        } else {
            a.isotope -= (f64::from(e.2) + 0.5) as i32;
        }
    }
    if at(p) == b'H' {
        p += 1;
        if at(p).is_ascii_digit() {
            a.hydrogen = i32::from(at(p) - b'0');
            p += 1;
        } else {
            a.hydrogen = 1;
        }
    }
    a.hydrogen += 1; // hydrogen is always explicit
    while at(p) == b'+' || at(p) == b'-' {
        if a.charge >= 0 && at(p) == b'+' {
            a.charge += 1;
        } else if a.charge <= 0 && at(p) == b'-' {
            a.charge -= 1;
        } else {
            break;
        }
        p += 1;
    }
    if (b'1'..=b'9').contains(&at(p)) && (a.charge == -1 || a.charge == 1) {
        a.charge *= i32::from(at(p) - b'0');
        p += 1;
    }
    if at(p).is_ascii_digit() && a.charge != 0 {
        if a.charge < 0 {
            a.charge = a.charge * 10 - i32::from(at(p) - b'0');
        } else {
            a.charge = a.charge * 10 + i32::from(at(p) - b'0');
        }
        p += 1;
    }
    if at(p) == b':' {
        p += 1;
        while at(p).is_ascii_digit() {
            a.mapping = a
                .mapping
                .wrapping_mul(10)
                .wrapping_add(i32::from(at(p) - b'0'));
            p += 1;
        }
        a.mapping += 1; // MDL mappings start with 1
    }
    if at(p) != b']' {
        sym.clear();
    }
    if sym == "*" || sym == "R" {
        if a.isotope > 0 && a.isotope <= 99 {
            sym = "R#".into();
            a.charge = a.isotope;
            a.isotope = 0;
        } else if a.charge > 0 && a.charge <= 99 {
            sym = "R#".into();
        } else if a.charge == -1 {
            sym = "Frm".into();
            a.charge = 0;
        } else if a.charge == -2 {
            sym = "To".into();
            a.charge = 0;
        } else {
            sym = "R".into();
        }
    }
    a.symbol = sym;
    (p, a)
}

/// `ComputeSmilesValenceAndHCount`: puts the SMILES valence and hydrogen count into the
/// `dummy1` and `dummy2` fields of the atoms
pub fn compute_smiles_valence_and_hcount(mp: &mut Molecule) -> Result<(), &'static str> {
    let nbp = setup_neighbourhood(mp).ok_or("too many neighbours")?;
    for i in 0..mp.atoms.len() {
        let mdl_valence = PTABLE
            .iter()
            .find(|e| e.0 == mp.atoms[i].atom_symbol)
            .map_or(0, |e| e.1);
        let (mut nsingle, mut ndouble, mut ntriple, mut naromatic) = (0, 0, 0, 0);
        for &b in &nbp[i].bonds {
            match mp.bonds[b].bond_type {
                SINGLE => nsingle += 1,
                DOUBLE => ndouble += 1,
                TRIPLE => ntriple += 1,
                _ => naromatic += 1,
            }
        }
        let ap = &mut mp.atoms[i];
        let mut valence = (nsingle * 2 + ndouble * 4 + ntriple * 6 + naromatic * 3) / 2;
        valence += ap.charge.abs();
        if ap.query_h_count > 0 {
            valence += ap.query_h_count - 1;
        }
        if valence == mdl_valence {
            ap.dummy2 = if ap.query_h_count > 0 {
                ap.query_h_count - 1
            } else {
                0
            };
            ap.dummy1 = mdl_valence;
        } else if valence < mdl_valence {
            if ap.query_h_count > 0 {
                ap.dummy2 = ap.query_h_count - 1;
                ap.dummy1 = if naromatic == 2 && nsingle == 0 && ndouble == 0 && mdl_valence == 4 {
                    2
                } else {
                    valence
                };
            } else {
                ap.dummy2 = mdl_valence - valence;
                ap.dummy1 = mdl_valence;
            }
        } else if naromatic > 0 {
            if ap.query_h_count > 0 {
                ap.dummy2 = ap.query_h_count - 1;
                ap.dummy1 = if valence > mdl_valence {
                    mdl_valence + 2 * ((valence - mdl_valence) / 2)
                } else {
                    valence
                };
            } else if valence - mdl_valence == 3 {
                ap.dummy2 = 0;
                ap.dummy1 = mdl_valence + 2;
            } else if valence - mdl_valence > 1 {
                ap.dummy2 = 0;
                ap.dummy1 = valence;
            } else {
                ap.dummy2 = 0;
                ap.dummy1 = mdl_valence;
            }
        } else if ap.query_h_count > 0 {
            ap.dummy2 = ap.query_h_count - 1;
            ap.dummy1 = valence;
        } else {
            ap.dummy2 = (valence - mdl_valence) % 2;
            ap.dummy1 = valence + ap.dummy2;
        }
    }
    Ok(())
}

const SMILES_DIGIT: i32 = 16;
const SMILES_BOND: i32 = 32;
const SMILES_UC: i32 = 1;
const SMILES_CLOSE: i32 = 4;
const SMILES_AROMATIC: i32 = 64;
const SMILES_ARO_BOND: i32 = 128;
const SMILES_BEFORE_RING: i32 =
    SMILES_UC | SMILES_AROMATIC | SMILES_ARO_BOND | SMILES_BOND | SMILES_CLOSE | SMILES_DIGIT;

/// `smiles_class[]` of `InitSMILESClasses` (bytes above 127 have no class)
fn smiles_class(c: u8) -> i32 {
    let mut r = 0;
    if c.is_ascii_uppercase() {
        r |= SMILES_UC;
    }
    if matches!(c, b')' | b']' | b'}') {
        r |= SMILES_CLOSE;
    }
    if c.is_ascii_digit() {
        r |= SMILES_DIGIT;
    }
    if matches!(c, b'-' | b'=' | b'#' | b'\\' | b'/' | b'.') {
        r |= SMILES_BOND;
    }
    if c == b':' {
        r |= SMILES_ARO_BOND;
    }
    if matches!(c, b'c' | b'n' | b'o' | b's' | b'p' | b'*') {
        r |= SMILES_AROMATIC;
    }
    r
}

/// `QuickCheck`: a quick test if `smiles` can be a legal SMILES (matching brackets and
/// ring closures, legal characters before ring digits)
pub fn quick_check(smiles: &str) -> bool {
    let b: Vec<u8> = smiles.bytes().take_while(|&c| c != b'|').collect();
    if smiles.is_empty() {
        return false;
    }
    for (open, close) in [(b'[', b']'), (b'(', b')'), (b'{', b'}')] {
        let (mut bra, mut ket) = (0, 0);
        for &c in &b {
            if c == open {
                bra += 1;
            } else if c == close {
                ket += 1;
            }
            if bra < ket {
                return false;
            }
        }
        if bra != ket {
            return false;
        }
    }
    if b.first().is_some_and(|c| c.is_ascii_digit() || *c == b'%') {
        return false;
    }
    let at = |i: usize| b.get(i).copied().unwrap_or(0);
    let mut closure_count = [0i32; 100];
    let (mut in_atom, mut in_shortcut) = (0, 0);
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c == b'{' {
            in_shortcut += 1;
        } else if c == b'}' {
            in_shortcut -= 1;
        }
        if c == b'[' {
            in_atom += 1;
        } else if c == b']' {
            in_atom -= 1;
        }
        if in_atom > 0 || in_shortcut > 0 {
            i += 1;
            continue;
        }
        if smiles_class(c) & SMILES_DIGIT != 0 {
            if i == 0 || smiles_class(b[i - 1]) & SMILES_BEFORE_RING == 0 {
                return false;
            }
            closure_count[usize::from(c - b'0')] += 1;
        } else if c == b'%' {
            if i == 0 || smiles_class(b[i - 1]) & SMILES_BEFORE_RING == 0 {
                return false;
            }
            if smiles_class(at(i + 1)) & SMILES_DIGIT == 0
                || smiles_class(at(i + 2)) & SMILES_DIGIT == 0
            {
                return false;
            }
            closure_count[usize::from(10 * (at(i + 1) - b'0') + (at(i + 2) - b'0'))] += 1;
            i += 2;
        }
        i += 1;
    }
    closure_count.iter().all(|n| n % 2 == 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::reaccsio::mol_str_to_mol;

    #[test]
    fn atoms_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_smiles_atoms.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let (pos, a) = smiles_to_mdl_atom(input.as_bytes(), 0);
            let got = if !input.starts_with('[') {
                format!("{pos} [] -")
            } else {
                format!(
                    "{pos} [{}] {} {} {} {} {} {} {}",
                    a.symbol,
                    a.charge,
                    a.radical,
                    a.isotope,
                    a.hydrogen,
                    a.stereo,
                    a.mapping,
                    i32::from(a.aromatic)
                )
            };
            let got = format!("{got}|{}", i32::from(quick_check(input)));
            if got != want {
                bad.push(format!("{input:?}\nwant {want}\ngot  {got}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    #[test]
    fn valence_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_smiles_valence.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mut mp = mol_str_to_mol(input).expect("molecule");
            compute_smiles_valence_and_hcount(&mut mp).unwrap();
            let mut got = String::new();
            for a in &mp.atoms {
                got += &format!("{} {}\n", a.dummy1, a.dummy2);
            }
            if got != want {
                bad.push(format!("{input}\nwant\n{want}\ngot\n{got}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
