//! SMILES tokenizer — port of RDKit `Code/GraphMol/SmilesParse/smiles.ll`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Reproduces flex semantics: longest match wins, ties go to the earlier
//! rule, and the `IN_ATOM_STATE` start condition (inside `[...]`) is
//! inclusive, i.e. the unconditional rules stay active inside brackets.

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondType};

/// A lexical token.
#[derive(Debug, Clone, PartialEq)]
pub enum Tok {
    /// `@TH`, `@AL`, `@SP`, `@TB`, `@OH`
    ChiClass(ChiralType),
    At,
    /// element symbol only valid inside brackets
    Atom(u32),
    /// `B C N O P S F Cl Br I *`
    Organic(u32),
    /// `b c n o p s` (and `si as se te` inside brackets)
    Aromatic(u32),
    H,
    Colon,
    Hash,
    /// bond symbols; `unspecified_order` is set for `/` and `\`
    Bond {
        bond_type: BondType,
        dir: BondDir,
        unspecified_order: bool,
        aromatic: bool,
    },
    /// `~` (any bond; a query bond)
    AnyBond,
    Minus,
    Plus,
    GroupOpen,
    GroupClose,
    AtomOpen,
    AtomClose,
    Separator,
    Percent,
    Zero,
    NonZeroDigit(i32),
    Eos,
    Bad,
}

const BRACKET_ELEMENTS: &[(&str, u32)] = &[
    ("He", 2),
    ("Li", 3),
    ("Be", 4),
    ("Ne", 10),
    ("Na", 11),
    ("Mg", 12),
    ("Al", 13),
    ("Si", 14),
    ("Ar", 18),
    ("K", 19),
    ("Ca", 20),
    ("Sc", 21),
    ("Ti", 22),
    ("V", 23),
    ("Cr", 24),
    ("Mn", 25),
    ("Fe", 26),
    ("Co", 27),
    ("Ni", 28),
    ("Cu", 29),
    ("Zn", 30),
    ("Ga", 31),
    ("Ge", 32),
    ("As", 33),
    ("Se", 34),
    ("Kr", 36),
    ("Rb", 37),
    ("Sr", 38),
    ("Y", 39),
    ("Zr", 40),
    ("Nb", 41),
    ("Mo", 42),
    ("Tc", 43),
    ("Ru", 44),
    ("Rh", 45),
    ("Pd", 46),
    ("Ag", 47),
    ("Cd", 48),
    ("In", 49),
    ("Sn", 50),
    ("Sb", 51),
    ("Te", 52),
    ("Xe", 54),
    ("Cs", 55),
    ("Ba", 56),
    ("La", 57),
    ("Ce", 58),
    ("Pr", 59),
    ("Nd", 60),
    ("Pm", 61),
    ("Sm", 62),
    ("Eu", 63),
    ("Gd", 64),
    ("Tb", 65),
    ("Dy", 66),
    ("Ho", 67),
    ("Er", 68),
    ("Tm", 69),
    ("Yb", 70),
    ("Lu", 71),
    ("Hf", 72),
    ("Ta", 73),
    ("W", 74),
    ("Re", 75),
    ("Os", 76),
    ("Ir", 77),
    ("Pt", 78),
    ("Au", 79),
    ("Hg", 80),
    ("Tl", 81),
    ("Pb", 82),
    ("Bi", 83),
    ("Po", 84),
    ("At", 85),
    ("Rn", 86),
    ("Fr", 87),
    ("Ra", 88),
    ("Ac", 89),
    ("Th", 90),
    ("Pa", 91),
    ("U", 92),
    ("Np", 93),
    ("Pu", 94),
    ("Am", 95),
    ("Cm", 96),
    ("Bk", 97),
    ("Cf", 98),
    ("Es", 99),
    ("Fm", 100),
    ("Md", 101),
    ("No", 102),
    ("Lr", 103),
    ("Rf", 104),
    ("Db", 105),
    ("Sg", 106),
    ("Bh", 107),
    ("Hs", 108),
    ("Mt", 109),
    ("Ds", 110),
    ("Rg", 111),
    ("Cn", 112),
    ("Nh", 113),
    ("Fl", 114),
    ("Mc", 115),
    ("Lv", 116),
    ("Ts", 117),
    ("Og", 118),
    ("Uun", 110),
    ("Uuu", 111),
    ("Uub", 112),
    ("Uut", 113),
    ("Uuq", 114),
    ("Uup", 115),
    ("Uuh", 116),
    ("Uus", 117),
    ("Uuo", 118),
    ("'Rf'", 104),
    ("'Db'", 105),
    ("'Sg'", 106),
    ("'Bh'", 107),
    ("'Hs'", 108),
    ("'Mt'", 109),
    ("'Ds'", 110),
    ("'Rg'", 111),
    ("'Cn'", 112),
    ("'Nh'", 113),
    ("'Fl'", 114),
    ("'Mc'", 115),
    ("'Lv'", 116),
    ("'Ts'", 117),
    ("'Og'", 118),
];

const BRACKET_AROMATIC: &[(&str, u32)] = &[("si", 14), ("as", 33), ("se", 34), ("te", 52)];

const ORGANIC: &[(&str, u32)] = &[
    ("B", 5),
    ("C", 6),
    ("N", 7),
    ("O", 8),
    ("P", 15),
    ("S", 16),
    ("F", 9),
    ("Cl", 17),
    ("Br", 35),
    ("I", 53),
];

const AROMATIC: &[(&str, u32)] = &[("b", 5), ("c", 6), ("n", 7), ("o", 8), ("p", 15), ("s", 16)];

/// Streaming tokenizer over a (trimmed) SMILES string.
#[derive(Debug, Clone)]
pub struct Lexer<'a> {
    input: &'a [u8],
    pos: usize,
    in_atom: bool,
    /// `current_token_position` from the flex scanner: end offset of the
    /// last token read (used in error messages).
    pub token_end: usize,
}

impl<'a> Lexer<'a> {
    pub fn new(input: &'a str) -> Self {
        Lexer {
            input: input.as_bytes(),
            pos: 0,
            in_atom: false,
            token_end: 0,
        }
    }

    fn starts_with(&self, s: &str) -> bool {
        self.input[self.pos..].starts_with(s.as_bytes())
    }

    /// Longest match among (text, value) tables; returns (len, value).
    fn best(&self, table: &[(&str, u32)]) -> Option<(usize, u32)> {
        let mut best: Option<(usize, u32)> = None;
        for &(t, v) in table {
            if self.starts_with(t) && best.is_none_or(|(l, _)| t.len() > l) {
                best = Some((t.len(), v));
            }
        }
        best
    }

    /// `@[' ']*XX`
    fn chi_class(&self) -> Option<(usize, ChiralType)> {
        let rest = &self.input[self.pos..];
        if rest.first() != Some(&b'@') {
            return None;
        }
        let mut i = 1;
        while rest.get(i) == Some(&b' ') {
            i += 1;
        }
        let two = rest.get(i..i + 2)?;
        let t = match two {
            b"TH" => ChiralType::Tetrahedral,
            b"AL" => ChiralType::Allene,
            b"SP" => ChiralType::SquarePlanar,
            b"TB" => ChiralType::TrigonalBipyramidal,
            b"OH" => ChiralType::Octahedral,
            _ => return None,
        };
        Some((i + 2, t))
    }

    pub fn next_token(&mut self) -> Tok {
        let (len, tok) = self.scan();
        self.pos += len;
        self.token_end += len;
        match tok {
            Tok::AtomOpen => self.in_atom = true,
            Tok::AtomClose => self.in_atom = false,
            _ => {}
        }
        tok
    }

    /// Peek at the next token without consuming it.
    pub fn peek(&self) -> Tok {
        self.scan().1
    }

    fn scan(&self) -> (usize, Tok) {
        let rest = &self.input[self.pos..];
        let Some(&c) = rest.first() else {
            return (0, Tok::Eos);
        };
        // candidates: (length, rule order, token); longest wins, then the
        // earliest rule in smiles.ll
        let mut cands: Vec<(usize, usize, Tok)> = Vec::with_capacity(4);
        if let Some((l, t)) = self.chi_class() {
            cands.push((l, 0, Tok::ChiClass(t)));
        }
        if c == b'@' {
            cands.push((1, 1, Tok::At));
        }
        if self.in_atom {
            if let Some((l, v)) = self.best(BRACKET_ELEMENTS) {
                cands.push((l, 2, Tok::Atom(v)));
            }
        }
        if let Some((l, v)) = self.best(ORGANIC) {
            cands.push((l, 3, Tok::Organic(v)));
        }
        if c == b'H' {
            cands.push((1, 4, Tok::H));
        }
        if let Some((l, v)) = self.best(AROMATIC) {
            cands.push((l, 5, Tok::Aromatic(v)));
        }
        if self.in_atom {
            if let Some((l, v)) = self.best(BRACKET_AROMATIC) {
                cands.push((l, 6, Tok::Aromatic(v)));
            }
        }
        if c == b'*' {
            cands.push((1, 7, Tok::Organic(0)));
        }
        if self.in_atom && c == b':' {
            cands.push((1, 8, Tok::Colon));
        }
        if self.in_atom && c == b'#' {
            cands.push((1, 9, Tok::Hash));
        }
        let bond = |t: BondType| Tok::Bond {
            bond_type: t,
            dir: BondDir::None,
            unspecified_order: false,
            aromatic: false,
        };
        match c {
            b'=' => cands.push((1, 20, bond(BondType::Double))),
            b'#' => cands.push((1, 21, bond(BondType::Triple))),
            b':' => cands.push((
                1,
                22,
                Tok::Bond {
                    bond_type: BondType::Aromatic,
                    dir: BondDir::None,
                    unspecified_order: false,
                    aromatic: true,
                },
            )),
            b'$' => cands.push((1, 23, bond(BondType::Quadruple))),
            b'~' => cands.push((1, 26, Tok::AnyBond)),
            b'\\' => {
                let l = if rest.get(1) == Some(&b'\\') { 2 } else { 1 };
                cands.push((
                    l,
                    27,
                    Tok::Bond {
                        bond_type: BondType::Unspecified,
                        dir: BondDir::EndDownRight,
                        unspecified_order: true,
                        aromatic: false,
                    },
                ));
            }
            b'/' => cands.push((
                1,
                28,
                Tok::Bond {
                    bond_type: BondType::Unspecified,
                    dir: BondDir::EndUpRight,
                    unspecified_order: true,
                    aromatic: false,
                },
            )),
            b'-' => {
                if rest.get(1) == Some(&b'>') {
                    cands.push((2, 24, bond(BondType::DativeR)));
                }
                cands.push((1, 29, Tok::Minus));
            }
            b'<' if rest.get(1) == Some(&b'-') => cands.push((2, 25, bond(BondType::DativeL))),
            b'+' => cands.push((1, 30, Tok::Plus)),
            b'(' => cands.push((1, 31, Tok::GroupOpen)),
            b')' => cands.push((1, 32, Tok::GroupClose)),
            b'[' => cands.push((1, 33, Tok::AtomOpen)),
            b']' if self.in_atom => cands.push((1, 34, Tok::AtomClose)),
            b'.' => cands.push((1, 35, Tok::Separator)),
            b'%' => cands.push((1, 36, Tok::Percent)),
            b'0' => cands.push((1, 37, Tok::Zero)),
            b'1'..=b'9' => cands.push((1, 38, Tok::NonZeroDigit(i32::from(c - b'0')))),
            b'\n' => return (1, Tok::Eos),
            _ => {}
        }
        cands
            .into_iter()
            .max_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)))
            .map(|(l, _, t)| (l, t))
            .unwrap_or((1, Tok::Bad))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toks(s: &str) -> Vec<Tok> {
        let mut l = Lexer::new(s);
        let mut v = Vec::new();
        loop {
            let t = l.next_token();
            if t == Tok::Eos {
                break;
            }
            v.push(t);
        }
        v
    }

    #[test]
    fn longest_match() {
        assert_eq!(toks("Cl"), vec![Tok::Organic(17)]);
        assert_eq!(toks("Sc"), vec![Tok::Organic(16), Tok::Aromatic(6)]);
        assert_eq!(
            toks("[Sc]"),
            vec![Tok::AtomOpen, Tok::Atom(21), Tok::AtomClose]
        );
        assert_eq!(
            toks("[Cn]"),
            vec![Tok::AtomOpen, Tok::Atom(112), Tok::AtomClose]
        );
        assert_eq!(
            toks("[Hg]"),
            vec![Tok::AtomOpen, Tok::Atom(80), Tok::AtomClose]
        );
        assert_eq!(
            toks("[se]"),
            vec![Tok::AtomOpen, Tok::Aromatic(34), Tok::AtomClose]
        );
        assert_eq!(toks("[C@ TH2]")[2], Tok::ChiClass(ChiralType::Tetrahedral));
        assert_eq!(
            toks("C->N")[1],
            Tok::Bond {
                bond_type: BondType::DativeR,
                dir: BondDir::None,
                unspecified_order: false,
                aromatic: false
            }
        );
        assert_eq!(toks("[O:1]")[2], Tok::Colon);
        assert_eq!(
            toks("C:C")[1],
            Tok::Bond {
                bond_type: BondType::Aromatic,
                dir: BondDir::None,
                unspecified_order: false,
                aromatic: true
            }
        );
        assert_eq!(toks("C\\\\C").len(), 3);
        assert_eq!(toks("X"), vec![Tok::Bad]);
    }
}
