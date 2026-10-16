//! Fractional atomic numbers of mancude rings (upstream `CIPLabeler/Mancude.{h,cpp}`).

use super::cip_mol::CipMol;
use crate::Result;
use std::cmp::Ordering;
use std::collections::VecDeque;

/// A rational number in lowest terms with a positive denominator
/// (`boost::rational<int>`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rational {
    num: i32,
    den: i32,
}

impl Rational {
    pub fn new(num: i32, den: i32) -> Self {
        fn gcd(a: i32, b: i32) -> i32 {
            if b == 0 { a.abs() } else { gcd(b, a % b) }
        }
        let g = gcd(num, den);
        if g == 0 {
            return Rational { num, den };
        }
        let sign = if den < 0 { -1 } else { 1 };
        Rational {
            num: sign * num / g,
            den: sign * den / g,
        }
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        (i64::from(self.num) * i64::from(other.den))
            .cmp(&(i64::from(other.num) * i64::from(self.den)))
    }
}

/// `FractionalAtomicNum`
#[derive(Debug, Clone, Copy)]
pub struct FractionalAtomicNum {
    pub(super) numerator: i32,
    pub(super) denominator: i32,
}

impl FractionalAtomicNum {
    pub fn new(numerator: i32, denominator: i32) -> Self {
        FractionalAtomicNum {
            numerator,
            denominator,
        }
    }

    pub fn value(&self) -> Rational {
        Rational::new(self.numerator, self.denominator)
    }

    pub fn is_averaged(&self) -> bool {
        self.denominator > 1
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Type {
    /// =C(X)-
    Cv4D3,
    /// =N-
    Nv3D2,
    /// =[N+]<
    Nv4D3Plus,
    /// -[N-]-
    Nv2D2Minus,
    /// -[C(X)-]-
    Cv3D3Minus,
    /// -[O+]=
    Ov3D2Plus,
    Other,
}

fn seed_types(types: &mut [Type], mol: &CipMol) -> Result<bool> {
    let mut result = false;
    for aidx in 0..mol.num_atoms() {
        let atom = mol.atom(aidx);
        let mut btypes = mol.total_num_hs(aidx)? as i32;
        let mut ring = false;
        for bond in mol.bonds_of(aidx) {
            btypes += match mol.bond_order(bond)? {
                1 => 0x0000_0001,
                2 => 0x0000_0100,
                _ => 0x0100_0000,
            };
            if mol.is_in_ring(bond)? {
                ring = true;
            }
        }
        if !ring {
            continue;
        }
        let q = atom.formal_charge();
        match atom.atomic_num() {
            6 | 14 | 32 => {
                if q == 0 && btypes == 0x0102 {
                    types[aidx as usize] = Type::Cv4D3;
                } else if q == -1 && btypes == 0x0003 {
                    types[aidx as usize] = Type::Cv3D3Minus;
                    result = true;
                }
            }
            7 | 15 | 33 => {
                if q == 0 && btypes == 0x0101 {
                    types[aidx as usize] = Type::Nv3D2;
                    result = true;
                } else if q == -1 && btypes == 0x0002 {
                    types[aidx as usize] = Type::Nv2D2Minus;
                    result = true;
                } else if q == 1 && btypes == 0x0102 {
                    types[aidx as usize] = Type::Nv4D3Plus;
                    result = true;
                }
            }
            8 if q == 1 && btypes == 0x0101 => {
                types[aidx as usize] = Type::Ov3D2Plus;
                result = true;
            }
            _ => {}
        }
    }
    Ok(result)
}

fn relax_types(types: &mut [Type], mol: &CipMol) -> Result<()> {
    let mut queue = VecDeque::new();
    let mut counts = vec![0i32; mol.num_atoms() as usize];
    for aidx in 0..mol.num_atoms() {
        if types[aidx as usize] == Type::Other {
            continue;
        }
        for bond in mol.bonds_of(aidx) {
            if !mol.is_in_ring(bond)? {
                continue;
            }
            let nbr = mol.other_atom(bond, aidx)?;
            if types[nbr as usize] != Type::Other {
                counts[aidx as usize] += 1;
            }
        }
        if counts[aidx as usize] <= 1 {
            queue.push_back(aidx);
        }
    }
    while let Some(aidx) = queue.pop_front() {
        if types[aidx as usize] == Type::Other {
            continue;
        }
        types[aidx as usize] = Type::Other;
        for bond in mol.bonds_of(aidx) {
            if !mol.is_in_ring(bond)? {
                continue;
            }
            let nbr = mol.other_atom(bond, aidx)?;
            if types[nbr as usize] == Type::Other {
                continue;
            }
            counts[nbr as usize] -= 1;
            if counts[nbr as usize] == 1 {
                queue.push_back(nbr);
            }
        }
    }
    Ok(())
}

fn visit_part(
    parts: &mut [i32],
    types: &[Type],
    part: i32,
    mut atom: u32,
    mol: &CipMol,
) -> Result<()> {
    loop {
        let mut next = None;
        for bond in mol.bonds_of(atom) {
            if !mol.is_in_ring(bond)? {
                continue;
            }
            let nbr = mol.other_atom(bond, atom)?;
            if parts[nbr as usize] == 0 && types[nbr as usize] != Type::Other {
                parts[nbr as usize] = part;
                if next.is_some() {
                    visit_part(parts, types, part, nbr, mol)?;
                } else {
                    next = Some(nbr);
                }
            }
        }
        match next {
            Some(n) => atom = n,
            None => return Ok(()),
        }
    }
}

fn visit_parts(parts: &mut [i32], types: &[Type], mol: &CipMol) -> Result<i32> {
    let mut numparts = 0;
    for aidx in 0..mol.num_atoms() {
        if parts[aidx as usize] == 0 && types[aidx as usize] != Type::Other {
            numparts += 1;
            parts[aidx as usize] = numparts;
            visit_part(parts, types, numparts, aidx, mol)?;
        }
    }
    Ok(numparts)
}

/// `calcFracAtomNums`
pub fn calc_frac_atom_nums(mol: &CipMol) -> Result<Vec<FractionalAtomicNum>> {
    let num_atoms = mol.num_atoms();
    let mut fractions: Vec<FractionalAtomicNum> = (0..num_atoms)
        .map(|i| FractionalAtomicNum::new(mol.atom(i).atomic_num() as i32, 1))
        .collect();
    let mut types = vec![Type::Other; num_atoms as usize];
    if !seed_types(&mut types, mol)? {
        return Ok(fractions);
    }
    relax_types(&mut types, mol)?;
    let mut parts = vec![0i32; num_atoms as usize];
    let numparts = visit_parts(&mut parts, &types, mol)?;
    if numparts <= 0 {
        return Ok(fractions);
    }
    let mut resonant_parts = vec![false; numparts as usize + 1];
    for i in 0..num_atoms as usize {
        if parts[i] == 0 {
            continue;
        }
        if types[i] == Type::Cv3D3Minus || types[i] == Type::Nv2D2Minus {
            resonant_parts[parts[i] as usize] = true;
            resonant_parts[0] = true;
        }
        let mut numerator = 0;
        let mut denominator = 0;
        for nbr in mol.neighbors_of(i as u32) {
            if parts[nbr as usize] == parts[i] {
                numerator += mol.atom(nbr).atomic_num() as i32;
                denominator += 1;
            }
        }
        fractions[i] = FractionalAtomicNum::new(numerator, denominator);
    }
    if !resonant_parts[0] {
        return Ok(fractions);
    }
    let mut numerators = vec![0i32; numparts as usize + 1];
    let mut denominators = vec![0i32; numparts as usize + 1];
    for i in 0..num_atoms {
        let part = parts[i as usize] as usize;
        if part == 0 || !resonant_parts[part] {
            continue;
        }
        denominators[part] += 1;
        for bond in mol.bonds_of(i) {
            let nbr = mol.other_atom(bond, i)?;
            let bord = mol.bond_order(bond)?;
            if bord > 1 && parts[nbr as usize] as usize == part {
                numerators[part] += (bord - 1) * mol.atom(nbr).atomic_num() as i32;
            }
        }
    }
    for i in 0..num_atoms as usize {
        let part = parts[i] as usize;
        if part != 0 && resonant_parts[part] && denominators[part] != 0 {
            fractions[i] = FractionalAtomicNum::new(numerators[part], denominators[part]);
        }
    }
    Ok(fractions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    fn fractions(smiles: &str) -> Vec<FractionalAtomicNum> {
        let mol = mol_from_smiles(smiles, &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        let cip = CipMol::new(&mol).unwrap();
        calc_frac_atom_nums(&cip).unwrap()
    }

    #[test]
    fn negative_resonance_component_gets_one_final_fraction() {
        for f in fractions("[CH-]1C=CC=C1") {
            assert_eq!((f.numerator, f.denominator), (24, 5));
            assert_eq!(f.value(), Rational::new(24, 5));
            assert!(f.is_averaged());
        }
    }

    #[test]
    fn unreduced_denominator_is_kept() {
        for f in fractions("[CH-]1C=C1") {
            assert_eq!((f.numerator, f.denominator), (12, 3));
            assert_eq!(f.value(), Rational::new(4, 1));
            assert!(f.is_averaged());
        }
    }

    #[test]
    fn typed_atom_outside_the_ring_two_core_is_relaxed() {
        let f = fractions("[N-]1CCC1")[0];
        assert_eq!((f.numerator, f.denominator), (7, 1));
        assert!(!f.is_averaged());
    }

    #[test]
    fn rationals_compare_by_value() {
        assert_eq!(Rational::new(4, 2), Rational::new(2, 1));
        assert!(Rational::new(1, 3) < Rational::new(1, 2));
        assert!(Rational::new(-1, -2) > Rational::new(1, 3));
    }
}
