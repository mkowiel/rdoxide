//! SMARTS parser: tokenizer (`smarts.ll`), grammar (`smarts.yy`) and the
//! `MolFromSmarts` / `AtomFromSmarts` / `BondFromSmarts` drivers
//! (SmilesParse.cpp).
//!
//! Ported from RDKit `Code/GraphMol/SmilesParse/` (2026.09.1).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! The LALR grammar is implemented as recursive descent. The precedence
//! that bison derives for atom/bond expressions is reproduced: implicit and
//! explicit `&` bind tightest (left-associative), then `,`, then `;`; `!` is
//! a prefix operator on a single primitive.

use std::collections::BTreeMap;

use super::ops::{self, SmilesParseError};
use crate::graphmol::atom::{Atom, ChiralType, HybridizationType};
use crate::graphmol::bond::{Bond, BondDir, BondType};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::*;
use crate::graphmol::romol::{CI_RIGHTMOST_ATOM, RWMol};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};

const SMARTS_H_MASK: u64 = 0x1;
const SMARTS_CHARGE_MASK: u64 = 0x2;

// ================================================================== lexer

/// The atom carried by an atom-query token.
#[derive(Debug, Clone, PartialEq)]
pub enum Tok {
    ChiClass(ChiralType),
    At,
    /// `ATOM_TOKEN`: an element symbol inside brackets
    Atom(u32),
    /// `COMPLEX_ATOM_QUERY_TOKEN` (D, d, X, v, R)
    Complex(Box<AtomQuery>),
    /// `RINGBOND_ATOM_QUERY_TOKEN` (x)
    RingBond,
    /// `HETERONEIGHBOR_ATOM_QUERY_TOKEN` (z)
    HeteroNbr,
    /// `ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN` (Z)
    AliphHeteroNbr,
    /// `IMPLICIT_H_ATOM_QUERY_TOKEN` (h)
    ImplicitH,
    /// `MIN_RINGSIZE_ATOM_QUERY_TOKEN` (r)
    MinRingSize,
    /// `RINGSIZE_ATOM_QUERY_TOKEN` (k)
    RingSize,
    H,
    Organic(u32),
    Aromatic(u32),
    /// `SIMPLE_ATOM_QUERY_TOKEN`: `*`, `a`, `A`
    Simple(char),
    Colon,
    Underscore,
    Hash,
    /// `BOND_TOKEN`
    Bond(char),
    Minus,
    Plus,
    BeginRecurse,
    GroupOpen,
    GroupClose,
    EndRecurse,
    RangeOpen,
    RangeClose,
    AtomOpen,
    AtomClose,
    Separator,
    Percent,
    Zero,
    NonZeroDigit(i32),
    Not,
    Semi,
    And,
    Or,
    Hyb(HybridizationType),
    Eos,
    Bad,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Initial,
    InAtom,
    InBranch,
    InRecursion,
}

/// Element symbols recognised inside brackets (the `<IN_ATOM_STATE>` list).
const BRACKET_ELEMENTS: &[&str] = &[
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Co", "Fe", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
    "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",
    "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv",
];

const ORGANIC: &[(&str, u32)] = &[
    ("B", 5),
    ("C", 6),
    ("N", 7),
    ("O", 8),
    ("F", 9),
    ("P", 15),
    ("S", 16),
    ("Cl", 17),
    ("Br", 35),
    ("I", 53),
];
const AROMATIC: &[(&str, u32)] = &[("b", 5), ("c", 6), ("n", 7), ("o", 8), ("p", 15), ("s", 16)];
const BRACKET_AROMATIC: &[(&str, u32)] = &[("si", 14), ("as", 33), ("se", 34), ("te", 52)];

/// Streaming SMARTS tokenizer with flex semantics (longest match, then
/// earliest rule; inclusive start conditions kept on a stack).
#[derive(Debug, Clone)]
pub struct Lexer<'a> {
    input: &'a [u8],
    pos: usize,
    states: Vec<State>,
    pub token_end: usize,
}

impl<'a> Lexer<'a> {
    pub fn new(input: &'a str) -> Self {
        Lexer {
            input: input.as_bytes(),
            pos: 0,
            states: vec![State::Initial],
            token_end: 0,
        }
    }

    fn state(&self) -> State {
        *self.states.last().unwrap_or(&State::Initial)
    }

    fn starts_with(&self, s: &str) -> bool {
        self.input[self.pos..].starts_with(s.as_bytes())
    }

    fn best<'t>(&self, table: impl IntoIterator<Item = (&'t str, u32)>) -> Option<(usize, u32)> {
        let mut best: Option<(usize, u32)> = None;
        for (t, v) in table {
            if self.starts_with(t) && best.is_none_or(|(l, _)| t.len() > l) {
                best = Some((t.len(), v));
            }
        }
        best
    }

    fn chi_class(&self) -> Option<(usize, ChiralType)> {
        let rest = &self.input[self.pos..];
        if rest.first() != Some(&b'@') {
            return None;
        }
        let mut i = 1;
        while rest.get(i) == Some(&b' ') {
            i += 1;
        }
        let t = match rest.get(i..i + 2)? {
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
            Tok::AtomOpen => self.states.push(State::InAtom),
            Tok::AtomClose if self.state() == State::InAtom => {
                self.states.pop();
            }
            Tok::GroupOpen => self.states.push(State::InBranch),
            Tok::BeginRecurse => self.states.push(State::InRecursion),
            Tok::GroupClose | Tok::EndRecurse => {
                self.states.pop();
            }
            _ => {}
        }
        tok
    }

    fn scan(&self) -> (usize, Tok) {
        let rest = &self.input[self.pos..];
        let Some(&c) = rest.first() else {
            return (0, Tok::Eos);
        };
        let in_atom = self.state() == State::InAtom;
        let mut cands: Vec<(usize, usize, Tok)> = Vec::with_capacity(4);
        if let Some((l, t)) = self.chi_class() {
            cands.push((l, 0, Tok::ChiClass(t)));
        }
        if c == b'@' {
            cands.push((1, 1, Tok::At));
        }
        if in_atom {
            let elems = BRACKET_ELEMENTS.iter().map(|&s| (s, 0u32));
            if let Some((l, _)) = self.best(elems) {
                let sym = std::str::from_utf8(&rest[..l]).unwrap_or("");
                if let Ok(n) = PeriodicTable::get().get_atomic_number(sym) {
                    cands.push((l, 2, Tok::Atom(n)));
                }
            }
            let q = match c {
                b'D' => Some((
                    3,
                    Tok::Complex(Box::new(make_atom_explicit_degree_query(1))),
                )),
                b'd' => Some((
                    4,
                    Tok::Complex(Box::new(make_atom_non_hydrogen_degree_query(1))),
                )),
                b'X' => Some((5, Tok::Complex(Box::new(make_atom_total_degree_query(1))))),
                b'x' => Some((6, Tok::RingBond)),
                b'v' => Some((7, Tok::Complex(Box::new(make_atom_total_valence_query(1))))),
                b'z' => Some((8, Tok::HeteroNbr)),
                b'Z' => Some((9, Tok::AliphHeteroNbr)),
                b'h' => Some((10, Tok::ImplicitH)),
                b'R' => Some((11, Tok::Complex(Box::new(make_atom_ring_query(-1))))),
                b'r' => Some((12, Tok::MinRingSize)),
                b'k' => Some((13, Tok::RingSize)),
                _ => None,
            };
            if let Some((order, t)) = q {
                cands.push((1, order, t));
            }
        }
        if c == b'H' {
            cands.push((1, 14, Tok::H));
        }
        if let Some((l, v)) = self.best(ORGANIC.iter().copied()) {
            cands.push((l, 15, Tok::Organic(v)));
        }
        if let Some((l, v)) = self.best(AROMATIC.iter().copied()) {
            cands.push((l, 16, Tok::Aromatic(v)));
        }
        if in_atom {
            if let Some((l, v)) = self.best(BRACKET_AROMATIC.iter().copied()) {
                cands.push((l, 17, Tok::Aromatic(v)));
            }
        }
        match c {
            b'*' | b'a' | b'A' => cands.push((1, 18, Tok::Simple(c as char))),
            b':' => cands.push((1, 19, Tok::Colon)),
            b'_' => cands.push((1, 20, Tok::Underscore)),
            b'#' => cands.push((1, 21, Tok::Hash)),
            b'=' | b'~' | b'$' | b'/' => cands.push((1, 22, Tok::Bond(c as char))),
            b'\\' => {
                let l = if rest.get(1) == Some(&b'\\') { 2 } else { 1 };
                cands.push((l, 22, Tok::Bond('\\')));
            }
            b'-' => {
                if rest.get(1) == Some(&b'>') {
                    cands.push((2, 23, Tok::Bond('>')));
                }
                cands.push((1, 25, Tok::Minus));
            }
            b'<' if rest.get(1) == Some(&b'-') => cands.push((2, 24, Tok::Bond('<'))),
            b'+' => cands.push((1, 26, Tok::Plus)),
            b'(' => cands.push((1, 28, Tok::GroupOpen)),
            b')' => match self.state() {
                State::InBranch => cands.push((1, 29, Tok::GroupClose)),
                State::InRecursion => cands.push((1, 30, Tok::EndRecurse)),
                _ => {}
            },
            b'{' => cands.push((1, 31, Tok::RangeOpen)),
            b'}' => cands.push((1, 32, Tok::RangeClose)),
            b'[' => cands.push((1, 33, Tok::AtomOpen)),
            b']' => cands.push((1, 34, Tok::AtomClose)),
            b'.' => cands.push((1, 35, Tok::Separator)),
            b'%' => cands.push((1, 36, Tok::Percent)),
            b'0' => cands.push((1, 37, Tok::Zero)),
            b'1'..=b'9' => cands.push((1, 38, Tok::NonZeroDigit(i32::from(c - b'0')))),
            b'!' => cands.push((1, 39, Tok::Not)),
            b';' => cands.push((1, 40, Tok::Semi)),
            b'&' => cands.push((1, 41, Tok::And)),
            b',' => cands.push((1, 42, Tok::Or)),
            b'^' => {
                let h = match rest.get(1) {
                    Some(b'0') => Some(HybridizationType::S),
                    Some(b'1') => Some(HybridizationType::SP),
                    Some(b'2') => Some(HybridizationType::SP2),
                    Some(b'3') => Some(HybridizationType::SP3),
                    Some(b'4') => Some(HybridizationType::SP3D),
                    Some(b'5') => Some(HybridizationType::SP3D2),
                    _ => None,
                };
                if let Some(h) = h {
                    cands.push((2, 43, Tok::Hyb(h)));
                }
            }
            b'\n' => return (1, Tok::Eos),
            _ => {}
        }
        if in_atom && c == b'$' && rest.get(1) == Some(&b'(') {
            cands.push((2, 27, Tok::BeginRecurse));
        }
        cands
            .into_iter()
            .max_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)))
            .map(|(l, _, t)| (l, t))
            .unwrap_or((1, Tok::Bad))
    }
}

// ================================================================== parser

#[derive(Debug, Clone)]
pub struct ParseFailure;

type PResult<T> = Result<T, ParseFailure>;

/// One molecule being built (`molList` entry) with its ring bookmarks.
#[derive(Debug, Default)]
struct MolBuilder {
    mol: RWMol,
    ring_atoms: BTreeMap<i32, Vec<u32>>,
    ring_bonds: BTreeMap<i32, Vec<Bond>>,
}

impl MolBuilder {
    fn active_atom(&self) -> u32 {
        self.mol
            .get_atom_with_bookmark(CI_RIGHTMOST_ATOM)
            .unwrap_or_else(|_| self.mol.num_atoms().saturating_sub(1))
    }
}

/// A query atom (`QueryAtom`): an atom whose `query` is set.
fn query_atom_num(num: u32) -> Atom {
    let mut a = Atom::new(num);
    a.query = Some(Box::new(make_atom_num_query(num as i32)));
    a
}

fn query_atom(q: AtomQuery) -> Atom {
    let mut a = Atom::new(0);
    a.query = Some(Box::new(q));
    a
}

/// `QueryAtom::expandQuery(what, how, maintainOrder=true)`
fn expand_atom(atom: &mut Atom, what: AtomQuery, how: CompositeQueryType) -> PResult<()> {
    let cur = atom.query.take().ok_or(ParseFailure)?;
    atom.query = Some(Box::new(
        expand_query(*cur, what, how, true, "Atom").map_err(|_| ParseFailure)?,
    ));
    Ok(())
}

/// `QueryBond::expandQuery(what, how, maintainOrder=true)`
fn expand_bond(bond: &mut Bond, what: BondQuery, how: CompositeQueryType) -> PResult<()> {
    let cur = bond.query.take().ok_or(ParseFailure)?;
    bond.query = Some(Box::new(
        expand_query(*cur, what, how, true, "Bond").map_err(|_| ParseFailure)?,
    ));
    Ok(())
}

/// `QueryBond(bT)`
fn query_bond(bt: BondType) -> Bond {
    let mut b = Bond::new(bt);
    b.query = Some(Box::new(if bt != BondType::Unspecified {
        make_bond_order_equals_query(bt)
    } else {
        make_bond_null_query()
    }));
    b
}

/// `QueryBond::setBondType` (replaces the query)
fn query_bond_set_type(b: &mut Bond, bt: BondType) {
    b.bond_type = bt;
    b.query = Some(Box::new(make_bond_order_equals_query(bt)));
}

/// `SmilesParseOps::ClearAtomChemicalProps`
fn clear_atom_chemical_props(a: &mut Atom) {
    a.set_isotope(0);
    a.set_formal_charge(0);
    a.set_num_explicit_hs(0);
}

/// `SmilesParseOps::getUnspecifiedQueryBond`
fn unspecified_query_bond(a1: &Atom, a2: Option<&Atom>) -> Bond {
    let mut b = if !a1.is_aromatic() || a2.is_some_and(|a| !a.is_aromatic()) {
        Bond::new(BondType::Single)
    } else {
        Bond::new(BondType::Aromatic)
    };
    b.query = Some(Box::new(make_single_or_aromatic_bond_query()));
    b.props.set_val(cp::_unspecifiedOrder, 1i32);
    b
}

fn atom_query_of(a: &Atom) -> PResult<AtomQuery> {
    a.query.as_deref().cloned().ok_or(ParseFailure)
}

fn is_number_start(t: &Tok) -> bool {
    matches!(t, Tok::Zero | Tok::NonZeroDigit(_))
}

/// Parser state (bison's parse-params).
#[derive(Debug)]
pub struct SmartsParser<'a> {
    input: &'a str,
    lex: Lexer<'a>,
    look: Tok,
    look_end: usize,
    num_bonds_parsed: u32,
    pub(crate) branch_points: Vec<(u32, usize)>,
}

impl<'a> SmartsParser<'a> {
    pub fn new(input: &'a str) -> Self {
        let mut lex = Lexer::new(input);
        let look = lex.next_token();
        let look_end = lex.token_end;
        SmartsParser {
            input,
            lex,
            look,
            look_end,
            num_bonds_parsed: 0,
            branch_points: Vec::new(),
        }
    }

    fn bump(&mut self) -> Tok {
        let t = std::mem::replace(&mut self.look, self.lex.next_token());
        self.look_end = self.lex.token_end;
        t
    }

    fn error(&self, msg: &str) -> ParseFailure {
        ops::print_syntax_error_message(self.input, msg, self.look_end, "SMARTS");
        ParseFailure
    }

    fn syntax_error(&self) -> ParseFailure {
        self.error("syntax error")
    }

    fn expect(&mut self, t: Tok) -> PResult<()> {
        if self.look == t {
            self.bump();
            Ok(())
        } else {
            Err(self.syntax_error())
        }
    }

    // ----------------------------------------------------------- numbers

    fn digit(&mut self) -> PResult<i32> {
        match self.look {
            Tok::Zero => {
                self.bump();
                Ok(0)
            }
            Tok::NonZeroDigit(d) => {
                self.bump();
                Ok(d)
            }
            _ => Err(self.syntax_error()),
        }
    }

    /// `number: ZERO_TOKEN | nonzero_number`
    fn number(&mut self) -> PResult<i32> {
        match self.look {
            Tok::Zero => {
                self.bump();
                Ok(0)
            }
            Tok::NonZeroDigit(_) => self.nonzero_number(),
            _ => Err(self.syntax_error()),
        }
    }

    fn nonzero_number(&mut self) -> PResult<i32> {
        let Tok::NonZeroDigit(mut v) = self.look else {
            return Err(self.syntax_error());
        };
        self.bump();
        while let Tok::Zero | Tok::NonZeroDigit(_) = self.look {
            let d = match self.look {
                Tok::NonZeroDigit(d) => d,
                _ => 0,
            };
            if v >= i32::MAX / 10 || v * 10 >= i32::MAX - d {
                return Err(self.error("number too large"));
            }
            self.bump();
            v = v * 10 + d;
        }
        Ok(v)
    }

    fn ring_number(&mut self) -> PResult<i32> {
        match self.look {
            Tok::Zero | Tok::NonZeroDigit(_) => self.digit(),
            Tok::Percent => {
                self.bump();
                match self.look {
                    Tok::NonZeroDigit(d1) => {
                        self.bump();
                        let d2 = self.digit()?;
                        Ok(d1 * 10 + d2)
                    }
                    Tok::GroupOpen => {
                        self.bump();
                        let mut v = 0;
                        let mut n = 0;
                        while n < 5 && matches!(self.look, Tok::Zero | Tok::NonZeroDigit(_)) {
                            v = v * 10 + self.digit()?;
                            n += 1;
                        }
                        if n == 0 {
                            return Err(self.syntax_error());
                        }
                        self.expect(Tok::GroupClose)?;
                        Ok(v)
                    }
                    _ => Err(self.syntax_error()),
                }
            }
            _ => Err(self.syntax_error()),
        }
    }

    fn charge_spec(&mut self) -> PResult<i32> {
        let plus = match self.look {
            Tok::Plus => true,
            Tok::Minus => false,
            _ => return Err(self.syntax_error()),
        };
        self.bump();
        let sign = if plus { 1 } else { -1 };
        match self.look {
            Tok::Plus if plus => {
                self.bump();
                Ok(2)
            }
            Tok::Minus if !plus => {
                self.bump();
                Ok(-2)
            }
            Tok::Zero | Tok::NonZeroDigit(_) => Ok(sign * self.number()?),
            _ => Ok(sign),
        }
    }

    // ----------------------------------------------------------- atoms

    /// `simple_atom`
    fn simple_atom(&mut self) -> PResult<Atom> {
        let t = self.bump();
        Ok(match t {
            Tok::Organic(n) => {
                let mut a = Atom::new(n);
                a.query = Some(Box::new(make_atom_type_query(n as i32, false)));
                a
            }
            Tok::Aromatic(n) => {
                let mut a = Atom::new(n);
                a.set_is_aromatic(true);
                a.query = Some(Box::new(make_atom_type_query(n as i32, true)));
                a
            }
            Tok::Simple('*') => query_atom(make_atom_null_query()),
            Tok::Simple('a') => {
                let mut a = query_atom(make_atom_aromatic_query());
                a.set_is_aromatic(true);
                a
            }
            Tok::Simple(_) => query_atom(make_atom_aliphatic_query()),
            _ => return Err(self.syntax_error()),
        })
    }

    fn starts_simple_atom(t: &Tok) -> bool {
        matches!(t, Tok::Organic(_) | Tok::Aromatic(_) | Tok::Simple(_))
    }

    /// `atomd`
    fn atomd(&mut self, mols: &mut Vec<MolBuilder>) -> PResult<Atom> {
        if Self::starts_simple_atom(&self.look) {
            return self.simple_atom();
        }
        if self.look != Tok::AtomOpen {
            return Err(self.syntax_error());
        }
        self.bump();
        if let Some(a) = self.try_hydrogen_atom()? {
            return Ok(a);
        }
        let mut a = self.atom_expr(mols)?;
        if self.look == Tok::Colon {
            self.bump();
            let n = self.number()?;
            a.props.set_val(cp::molAtomMapNumber, n);
        }
        self.expect(Tok::AtomClose)?;
        Ok(a)
    }

    /// The `hydrogen_atom` productions: `[H]`, `[2H]`, `[H+]`, `[2H-:3]`, ...
    /// bison commits to them (shift over reduce) once it has seen
    /// `[ (number)? H` followed by `]`, `+`, `-` or `:`; afterwards a
    /// malformed tail is a syntax error, not a fallback to `atom_expr`.
    fn try_hydrogen_atom(&mut self) -> PResult<Option<Atom>> {
        let mut probe = self.lex.clone();
        let mut t = self.look.clone();
        let mut saw_number = false;
        let mut leading_zero = false;
        while let Tok::Zero | Tok::NonZeroDigit(_) = t {
            if leading_zero {
                return Ok(None);
            }
            leading_zero = !saw_number && t == Tok::Zero;
            saw_number = true;
            t = probe.next_token();
        }
        if t != Tok::H {
            return Ok(None);
        }
        let after = probe.next_token();
        if !matches!(after, Tok::AtomClose | Tok::Plus | Tok::Minus | Tok::Colon) {
            return Ok(None);
        }
        let iso = if saw_number {
            Some(self.number()?)
        } else {
            None
        };
        self.expect(Tok::H)?;
        let charge = if matches!(self.look, Tok::Plus | Tok::Minus) {
            Some(self.charge_spec()?)
        } else {
            None
        };
        let map = if self.look == Tok::Colon {
            self.bump();
            Some(self.number()?)
        } else {
            None
        };
        self.expect(Tok::AtomClose)?;
        let mut a = query_atom_num(1);
        if let Some(n) = iso {
            a.set_isotope(n as u32);
        }
        if let Some(c) = charge {
            a.set_formal_charge(c);
            a.flags |= SMARTS_CHARGE_MASK;
        }
        if let Some(n) = iso {
            expand_atom(&mut a, make_atom_isotope_query(n), CompositeQueryType::And)?;
        }
        if let Some(c) = charge {
            expand_atom(
                &mut a,
                make_atom_formal_charge_query(c),
                CompositeQueryType::And,
            )?;
        }
        if let Some(m) = map {
            a.props.set_val(cp::molAtomMapNumber, m);
        }
        Ok(Some(a))
    }

    /// tokens that can start a `point_query`
    fn starts_point_query(t: &Tok) -> bool {
        Self::starts_simple_atom(t)
            || matches!(
                t,
                Tok::Not
                    | Tok::BeginRecurse
                    | Tok::Atom(_)
                    | Tok::Hash
                    | Tok::Zero
                    | Tok::NonZeroDigit(_)
                    | Tok::H
                    | Tok::At
                    | Tok::ChiClass(_)
                    | Tok::Hyb(_)
                    | Tok::Plus
                    | Tok::Minus
                    | Tok::Complex(_)
                    | Tok::HeteroNbr
                    | Tok::AliphHeteroNbr
                    | Tok::MinRingSize
                    | Tok::RingSize
                    | Tok::RingBond
                    | Tok::ImplicitH
            )
    }

    /// `atom_expr` (lowest precedence: `;`)
    fn atom_expr(&mut self, mols: &mut Vec<MolBuilder>) -> PResult<Atom> {
        let mut lhs = self.or_expr(mols)?;
        while self.look == Tok::Semi {
            self.bump();
            let rhs = self.or_expr(mols)?;
            expand_atom(&mut lhs, atom_query_of(&rhs)?, CompositeQueryType::And)?;
            if lhs.chiral_tag() == ChiralType::Unspecified {
                lhs.set_chiral_tag(rhs.chiral_tag());
            }
            clear_atom_chemical_props(&mut lhs);
        }
        Ok(lhs)
    }

    fn or_expr(&mut self, mols: &mut Vec<MolBuilder>) -> PResult<Atom> {
        let mut lhs = self.and_expr(mols)?;
        while self.look == Tok::Or {
            self.bump();
            let rhs = self.and_expr(mols)?;
            expand_atom(&mut lhs, atom_query_of(&rhs)?, CompositeQueryType::Or)?;
            if lhs.chiral_tag() == ChiralType::Unspecified {
                lhs.set_chiral_tag(rhs.chiral_tag());
            }
            clear_atom_chemical_props(&mut lhs);
            lhs.set_atomic_num(0);
        }
        Ok(lhs)
    }

    /// `atom_expr point_query` and `atom_expr AND_TOKEN point_query`
    fn and_expr(&mut self, mols: &mut Vec<MolBuilder>) -> PResult<Atom> {
        let mut lhs = self.point_query(mols)?;
        loop {
            if self.look == Tok::And {
                self.bump();
                let rhs = self.point_query(mols)?;
                self.atom_expr_and_point_query(&mut lhs, rhs)?;
            } else if Self::starts_point_query(&self.look) {
                let rhs = self.point_query(mols)?;
                self.atom_expr_and_point_query(&mut lhs, rhs)?;
            } else {
                break;
            }
        }
        Ok(lhs)
    }

    /// `atom_expr_and_point_query`
    fn atom_expr_and_point_query(&mut self, ae: &mut Atom, pq: Atom) -> PResult<()> {
        expand_atom(ae, atom_query_of(&pq)?, CompositeQueryType::And)?;
        if ae.chiral_tag() == ChiralType::Unspecified {
            ae.set_chiral_tag(pq.chiral_tag());
            if let Some(p) = pq
                .props
                .get(cp::_chiralPermutation)
                .and_then(RDValue::as_int)
            {
                ae.props.set_val(cp::_chiralPermutation, p);
            }
        }
        if pq.flags & SMARTS_H_MASK != 0 {
            if ae.flags & SMARTS_H_MASK == 0 {
                ae.set_num_explicit_hs(pq.num_explicit_hs());
                ae.set_no_implicit(true);
                ae.flags |= SMARTS_H_MASK;
            } else if ae.num_explicit_hs() != pq.num_explicit_hs() {
                ae.set_num_explicit_hs(0);
                ae.set_no_implicit(true);
            }
        }
        if pq.flags & SMARTS_CHARGE_MASK != 0 {
            if ae.flags & SMARTS_CHARGE_MASK == 0 {
                ae.set_formal_charge(pq.formal_charge());
                ae.flags |= SMARTS_CHARGE_MASK;
            } else if ae.formal_charge() != pq.formal_charge() {
                ae.set_formal_charge(0);
            }
        }
        Ok(())
    }

    /// `point_query`
    fn point_query(&mut self, mols: &mut Vec<MolBuilder>) -> PResult<Atom> {
        match self.look {
            Tok::Not => {
                self.bump();
                let mut a = self.point_query(mols)?;
                if let Some(q) = a.query.as_deref_mut() {
                    q.negate = !q.negate;
                }
                a.set_atomic_num(0);
                clear_atom_chemical_props(&mut a);
                Ok(a)
            }
            Tok::BeginRecurse => self.recursive_query(mols),
            _ => self.atom_query(),
        }
    }

    /// `recursive_query`
    fn recursive_query(&mut self, mols: &mut Vec<MolBuilder>) -> PResult<Atom> {
        self.bump(); // BEGIN_RECURSE
        let idx = mols.len();
        self.mol(mols)?;
        if self.look != Tok::EndRecurse {
            return Err(self.syntax_error());
        }
        self.bump();
        let mut serial = 0u32;
        if self.look == Tok::Underscore {
            self.bump();
            if !matches!(self.look, Tok::NonZeroDigit(_)) {
                return Err(self.syntax_error());
            }
            serial = self.nonzero_number()? as u32;
        }
        let mut mb = mols.remove(idx);
        ops::close_mol_rings(&mut mb.mol, &mut mb.ring_atoms, &mut mb.ring_bonds, false).map_err(
            |SmilesParseError(e)| {
                crate::rd_log_error!("SMARTS Parse Error: {e} for input: '{}'", self.input);
                ParseFailure
            },
        )?;
        Ok(query_atom(make_recursive_structure_query(mb.mol, serial)))
    }

    /// `atom_query` (plus the forms that need a leading number or token)
    fn atom_query(&mut self) -> PResult<Atom> {
        let t = self.look.clone();
        match t {
            ref t if Self::starts_simple_atom(t) => self.simple_atom(),
            Tok::Atom(n) => {
                self.bump();
                Ok(query_atom_num(n))
            }
            Tok::Hash => {
                self.bump();
                let n = self.number()?;
                Ok(query_atom_num(n as u32))
            }
            Tok::Zero | Tok::NonZeroDigit(_) => {
                let n = self.number()?;
                match self.look.clone() {
                    ref t if Self::starts_simple_atom(t) => {
                        let mut a = self.simple_atom()?;
                        a.set_isotope(n as u32);
                        expand_atom(&mut a, make_atom_isotope_query(n), CompositeQueryType::And)?;
                        Ok(a)
                    }
                    Tok::Atom(z) => {
                        self.bump();
                        let mut a = query_atom_num(z);
                        a.set_isotope(n as u32);
                        expand_atom(&mut a, make_atom_isotope_query(n), CompositeQueryType::And)?;
                        Ok(a)
                    }
                    Tok::Hash => {
                        self.bump();
                        let z = self.number()?;
                        let mut a = query_atom_num(z as u32);
                        a.set_isotope(n as u32);
                        expand_atom(&mut a, make_atom_isotope_query(n), CompositeQueryType::And)?;
                        Ok(a)
                    }
                    Tok::H => {
                        self.bump();
                        let hc = if is_number_start(&self.look) {
                            self.number()?
                        } else {
                            1
                        };
                        let mut a = query_atom(make_atom_isotope_query(n));
                        a.set_isotope(n as u32);
                        expand_atom(&mut a, make_atom_h_count_query(hc), CompositeQueryType::And)?;
                        a.set_num_explicit_hs(hc as u32);
                        a.set_no_implicit(true);
                        a.flags |= SMARTS_H_MASK;
                        Ok(a)
                    }
                    _ => Ok(query_atom(make_atom_isotope_query(n))),
                }
            }
            Tok::H => {
                self.bump();
                let hc = if is_number_start(&self.look) {
                    self.number()?
                } else {
                    1
                };
                let mut a = query_atom(make_atom_h_count_query(hc));
                a.set_num_explicit_hs(hc as u32);
                a.set_no_implicit(true);
                a.flags |= SMARTS_H_MASK;
                Ok(a)
            }
            Tok::At => {
                self.bump();
                let mut a = query_atom(make_atom_null_query());
                if self.look == Tok::At {
                    self.bump();
                    a.set_chiral_tag(ChiralType::TetrahedralCW);
                } else {
                    a.set_chiral_tag(ChiralType::TetrahedralCCW);
                }
                Ok(a)
            }
            Tok::ChiClass(tag) => {
                self.bump();
                let mut perm = 0;
                if is_number_start(&self.look) {
                    perm = self.number()?;
                    if perm == 0 {
                        return Err(self.error("chiral permutation cannot be zero"));
                    }
                }
                let mut a = query_atom(make_atom_null_query());
                a.set_chiral_tag(tag);
                a.props.set_val(cp::_chiralPermutation, perm);
                Ok(a)
            }
            Tok::Hyb(h) => {
                self.bump();
                Ok(query_atom(make_atom_hybridization_query(h as i32)))
            }
            Tok::Plus | Tok::Minus => {
                // a range query on charge (`+{..}`), or a charge spec
                if self.peek_is_range_open() {
                    let base = if t == Tok::Plus {
                        make_atom_formal_charge_query(0)
                    } else {
                        make_atom_negative_formal_charge_query(0)
                    };
                    self.bump();
                    return self.range_suffix(query_atom(base));
                }
                let c = self.charge_spec()?;
                let mut a = query_atom(make_atom_formal_charge_query(c));
                a.set_formal_charge(c);
                a.flags |= SMARTS_CHARGE_MASK;
                Ok(a)
            }
            Tok::Complex(q) => {
                self.bump();
                let mut a = query_atom(*q);
                if is_number_start(&self.look) {
                    let n = self.number()?;
                    if let Some(q) = a.query.as_deref_mut() {
                        q.val = n;
                    }
                    Ok(a)
                } else if self.look == Tok::RangeOpen {
                    self.range_suffix(a)
                } else {
                    Ok(a)
                }
            }
            Tok::HeteroNbr
            | Tok::AliphHeteroNbr
            | Tok::MinRingSize
            | Tok::RingBond
            | Tok::ImplicitH => {
                self.bump();
                if is_number_start(&self.look) {
                    let n = self.number()?;
                    let q = match t {
                        Tok::HeteroNbr => make_atom_num_heteroatom_nbrs_query(n),
                        Tok::AliphHeteroNbr => make_atom_num_aliphatic_heteroatom_nbrs_query(n),
                        Tok::MinRingSize => make_atom_min_ring_size_query(n),
                        Tok::RingBond => make_atom_ring_bond_count_query(n),
                        _ => make_atom_implicit_h_count_query(n),
                    };
                    Ok(query_atom(q))
                } else if self.look == Tok::RangeOpen {
                    let q = match t {
                        Tok::HeteroNbr => make_atom_num_heteroatom_nbrs_query(0),
                        Tok::AliphHeteroNbr => make_atom_num_aliphatic_heteroatom_nbrs_query(0),
                        Tok::MinRingSize => make_atom_min_ring_size_query(5),
                        Tok::RingBond => make_atom_ring_bond_count_query(0),
                        _ => make_atom_implicit_h_count_query(0),
                    };
                    self.range_suffix(query_atom(q))
                } else {
                    let q = match t {
                        Tok::HeteroNbr => make_atom_has_heteroatom_nbrs_query(),
                        Tok::AliphHeteroNbr => make_atom_has_aliphatic_heteroatom_nbrs_query(),
                        Tok::MinRingSize => make_atom_in_ring_query(),
                        Tok::RingBond => make_atom_has_ring_bond_query(),
                        _ => make_atom_has_implicit_h_query(),
                    };
                    Ok(query_atom(q))
                }
            }
            Tok::RingSize => {
                self.bump();
                if is_number_start(&self.look) {
                    let n = self.number()?;
                    Ok(query_atom(make_atom_in_ring_of_size_query(n)))
                } else if self.look == Tok::RangeOpen {
                    let (lo, hi) = self.range_bounds()?;
                    let q = match (lo, hi) {
                        (None, Some(uv)) => {
                            let mut q = Query::new(QueryKind::GreaterEqual);
                            q.val = uv;
                            q.data_func = Some(AtomFn::RingSizeRange {
                                lower: -1,
                                upper: uv,
                                lower_open: false,
                                upper_open: false,
                            });
                            q.description = "greater_AtomRingSize".into();
                            q
                        }
                        (Some(lv), None) => {
                            let mut q = Query::new(QueryKind::LessEqual);
                            q.val = lv;
                            q.data_func = Some(AtomFn::RingSizeRange {
                                lower: lv,
                                upper: -1,
                                lower_open: false,
                                upper_open: false,
                            });
                            q.description = "less_AtomRingSize".into();
                            q
                        }
                        (Some(lv), Some(uv)) => make_atom_range_query(
                            lv,
                            uv,
                            false,
                            false,
                            AtomFn::RingSizeRange {
                                lower: lv,
                                upper: uv,
                                lower_open: false,
                                upper_open: false,
                            },
                            "range_AtomRingSize",
                        ),
                        (None, None) => return Err(self.syntax_error()),
                    };
                    Ok(query_atom(q))
                } else {
                    Ok(query_atom(make_atom_in_ring_query()))
                }
            }
            _ => Err(self.syntax_error()),
        }
    }

    fn peek_is_range_open(&self) -> bool {
        let mut probe = self.lex.clone();
        probe.next_token() == Tok::RangeOpen
    }

    /// `RANGE_OPEN (MINUS number | number MINUS | number MINUS number) RANGE_CLOSE`
    fn range_bounds(&mut self) -> PResult<(Option<i32>, Option<i32>)> {
        self.expect(Tok::RangeOpen)?;
        let res = if self.look == Tok::Minus {
            self.bump();
            (None, Some(self.number()?))
        } else {
            let lo = self.number()?;
            self.expect(Tok::Minus)?;
            if self.look == Tok::RangeClose {
                (Some(lo), None)
            } else {
                (Some(lo), Some(self.number()?))
            }
        };
        self.expect(Tok::RangeClose)?;
        Ok(res)
    }

    /// The `possible_range_query RANGE_OPEN ...` productions
    fn range_suffix(&mut self, mut a: Atom) -> PResult<Atom> {
        let (lo, hi) = self.range_bounds()?;
        let oq = atom_query_of(&a)?;
        let func = oq.data_func.clone().ok_or(ParseFailure)?;
        let q = match (lo, hi) {
            (None, Some(v)) => {
                let mut q = Query::new(QueryKind::GreaterEqual);
                q.val = v;
                q.data_func = Some(func);
                q.description = format!("greater_{}", oq.description);
                q
            }
            (Some(v), None) => {
                let mut q = Query::new(QueryKind::LessEqual);
                q.val = v;
                q.data_func = Some(func);
                q.description = format!("less_{}", oq.description);
                q
            }
            (Some(l), Some(u)) => make_atom_range_query(
                l,
                u,
                false,
                false,
                func,
                &format!("range_{}", oq.description),
            ),
            (None, None) => return Err(self.syntax_error()),
        };
        a.query = Some(Box::new(q));
        Ok(a)
    }

    // ----------------------------------------------------------- bonds

    fn starts_bondd(t: &Tok) -> bool {
        matches!(
            t,
            Tok::Bond(_) | Tok::Minus | Tok::Hash | Tok::Colon | Tok::At | Tok::Not
        )
    }

    /// `bondd`
    fn bondd(&mut self) -> PResult<Bond> {
        match self.bump() {
            Tok::Bond(c) => Ok(match c {
                '=' => query_bond(BondType::Double),
                '~' => {
                    let mut b = Bond::new(BondType::Unspecified);
                    b.query = Some(Box::new(make_bond_null_query()));
                    b
                }
                '$' => query_bond(BondType::Quadruple),
                '\\' | '/' => {
                    let mut b = query_bond(BondType::Single);
                    b.set_bond_dir(if c == '\\' {
                        BondDir::EndDownRight
                    } else {
                        BondDir::EndUpRight
                    });
                    b.query = Some(Box::new(make_single_or_aromatic_bond_query()));
                    b
                }
                '>' => query_bond(BondType::DativeR),
                _ => query_bond(BondType::DativeL),
            }),
            Tok::Minus => {
                let mut b = Bond::new(BondType::Unspecified);
                b.query = Some(Box::new(make_bond_null_query()));
                query_bond_set_type(&mut b, BondType::Single);
                Ok(b)
            }
            Tok::Hash => {
                let mut b = Bond::new(BondType::Unspecified);
                query_bond_set_type(&mut b, BondType::Triple);
                Ok(b)
            }
            Tok::Colon => {
                let mut b = Bond::new(BondType::Unspecified);
                query_bond_set_type(&mut b, BondType::Aromatic);
                Ok(b)
            }
            Tok::At => {
                let mut b = Bond::new(BondType::Unspecified);
                b.query = Some(Box::new(make_bond_is_in_ring_query()));
                Ok(b)
            }
            Tok::Not => {
                let mut b = self.bondd()?;
                if let Some(q) = b.query.as_deref_mut() {
                    q.negate = !q.negate;
                }
                Ok(b)
            }
            _ => Err(self.syntax_error()),
        }
    }

    fn bond_query(&mut self) -> PResult<Bond> {
        let mut b = self.bondd()?;
        while Self::starts_bondd(&self.look) {
            let r = self.bondd()?;
            let q = r.query.as_deref().cloned().ok_or(ParseFailure)?;
            expand_bond(&mut b, q, CompositeQueryType::And)?;
        }
        Ok(b)
    }

    fn bond_and(&mut self) -> PResult<Bond> {
        let mut b = self.bond_query()?;
        while self.look == Tok::And {
            self.bump();
            let r = self.bond_query()?;
            let q = r.query.as_deref().cloned().ok_or(ParseFailure)?;
            expand_bond(&mut b, q, CompositeQueryType::And)?;
        }
        Ok(b)
    }

    fn bond_or(&mut self) -> PResult<Bond> {
        let mut b = self.bond_and()?;
        while self.look == Tok::Or {
            self.bump();
            let r = self.bond_and()?;
            let q = r.query.as_deref().cloned().ok_or(ParseFailure)?;
            expand_bond(&mut b, q, CompositeQueryType::Or)?;
        }
        Ok(b)
    }

    /// `bond_expr`
    fn bond_expr(&mut self) -> PResult<Bond> {
        let mut b = self.bond_or()?;
        while self.look == Tok::Semi {
            self.bump();
            let r = self.bond_or()?;
            let q = r.query.as_deref().cloned().ok_or(ParseFailure)?;
            expand_bond(&mut b, q, CompositeQueryType::And)?;
        }
        Ok(b)
    }

    // ----------------------------------------------------------- molecule

    fn starts_atomd(t: &Tok) -> bool {
        Self::starts_simple_atom(t) || *t == Tok::AtomOpen
    }

    fn starts_ring_number(t: &Tok) -> bool {
        matches!(t, Tok::Zero | Tok::NonZeroDigit(_) | Tok::Percent)
    }

    fn add_bonded_atom(&mut self, mb: &mut MolBuilder, mut bond: Bond, atom: Atom) -> PResult<()> {
        let a1 = mb.active_atom();
        let a2 = mb.mol.add_atom(atom, true);
        match bond.bond_type() {
            BondType::DativeR => {
                bond.begin = a1;
                bond.end = a2;
                query_bond_set_type(&mut bond, BondType::Dative);
            }
            BondType::DativeL => {
                bond.begin = a2;
                bond.end = a1;
                query_bond_set_type(&mut bond, BondType::Dative);
            }
            _ => {
                bond.begin = a1;
                bond.end = a2;
            }
        }
        bond.props
            .set_val("_cxsmilesBondIdx", self.num_bonds_parsed);
        self.num_bonds_parsed += 1;
        mb.mol
            .add_bond_object(bond)
            .map_err(|_| self.syntax_error())?;
        Ok(())
    }

    fn add_unspecified(&mut self, mb: &mut MolBuilder, atom: Atom) -> PResult<()> {
        let a1 = mb.active_atom();
        let a2 = mb.mol.add_atom(atom, true);
        let mut b = unspecified_query_bond(mb.mol.atom(a1), Some(mb.mol.atom(a2)));
        b.begin = a1;
        b.end = a2;
        b.props.set_val("_cxsmilesBondIdx", self.num_bonds_parsed);
        self.num_bonds_parsed += 1;
        mb.mol.add_bond_object(b).map_err(|_| self.syntax_error())?;
        Ok(())
    }

    fn ring_closure(
        &mut self,
        mb: &mut MolBuilder,
        num: i32,
        mut partial: Bond,
        count_on_pair: bool,
    ) {
        let atom = mb.active_atom();
        partial.begin = atom;
        let list = mb.ring_bonds.entry(num).or_default();
        list.push(partial);
        if !count_on_pair || list.len() % 2 == 0 {
            let b = list.last_mut().expect("just pushed");
            b.props.set_val("_cxsmilesBondIdx", self.num_bonds_parsed);
            self.num_bonds_parsed += 1;
        }
        mb.ring_atoms.entry(num).or_default().push(atom);
        ops::check_ring_closure_branch_status(&mut mb.mol, atom);
        let mut closures = match mb.mol.atom(atom).props.get(cp::_RingClosures) {
            Some(RDValue::VecInt(v)) => v.clone(),
            _ => Vec::new(),
        };
        closures.push(-(num + 1));
        mb.mol
            .atom_mut(atom)
            .props
            .set_val(cp::_RingClosures, closures);
    }

    /// `mol` (the molecule is pushed onto `mols`)
    fn mol(&mut self, mols: &mut Vec<MolBuilder>) -> PResult<()> {
        let mut first = self.atomd(mols)?;
        first.props.set_val(cp::_SmilesStart, 1i32);
        let mut mb = MolBuilder::default();
        mb.mol.add_atom(first, true);
        loop {
            let t = self.look.clone();
            match t {
                ref t if Self::starts_atomd(t) => {
                    let a = self.atomd(mols)?;
                    self.add_unspecified(&mut mb, a)?;
                }
                ref t if Self::starts_bondd(t) => {
                    let bond = self.bond_expr()?;
                    if Self::starts_atomd(&self.look) {
                        let a = self.atomd(mols)?;
                        self.add_bonded_atom(&mut mb, bond, a)?;
                    } else if Self::starts_ring_number(&self.look) {
                        let num = self.ring_number()?;
                        self.ring_closure(&mut mb, num, bond, false);
                    } else {
                        return Err(self.syntax_error());
                    }
                }
                Tok::Separator => {
                    self.bump();
                    let mut a = self.atomd(mols)?;
                    a.set_prop(cp::_SmilesStart, 1i32, true).expect("valid key");
                    mb.mol.add_atom(a, true);
                }
                ref t if Self::starts_ring_number(t) => {
                    let num = self.ring_number()?;
                    let atom = mb.active_atom();
                    let partial = unspecified_query_bond(mb.mol.atom(atom), None);
                    self.ring_closure(&mut mb, num, partial, true);
                }
                Tok::GroupOpen => {
                    let open_pos = self.look_end;
                    self.bump();
                    let branch_atom = mb.active_atom();
                    if Self::starts_atomd(&self.look) {
                        let a = self.atomd(mols)?;
                        self.add_unspecified(&mut mb, a)?;
                    } else if Self::starts_bondd(&self.look) {
                        let bond = self.bond_expr()?;
                        let a = self.atomd(mols)?;
                        self.add_bonded_atom(&mut mb, bond, a)?;
                    } else {
                        return Err(self.syntax_error());
                    }
                    self.branch_points.push((branch_atom, open_pos));
                }
                Tok::GroupClose => {
                    let Some((atom, _)) = self.branch_points.pop() else {
                        return Err(self.error("extra close parentheses"));
                    };
                    self.bump();
                    mb.mol.replace_atom_bookmark(atom, CI_RIGHTMOST_ATOM);
                }
                _ => break,
            }
        }
        mols.push(mb);
        Ok(())
    }

    /// `meta_start: START_MOL mol` + `EOS_TOKEN`
    fn parse_mol(&mut self) -> PResult<MolBuilder> {
        let mut mols = Vec::new();
        self.mol(&mut mols)?;
        if self.look != Tok::Eos {
            return Err(self.syntax_error());
        }
        mols.pop().ok_or(ParseFailure)
    }

    /// `meta_start: START_ATOM atomd EOS_TOKEN`
    fn parse_atom(&mut self) -> PResult<Atom> {
        let mut mols = Vec::new();
        let a = self.atomd(&mut mols)?;
        if self.look != Tok::Eos {
            return Err(ParseFailure);
        }
        Ok(a)
    }

    /// `meta_start: START_BOND bond_expr EOS_TOKEN`
    fn parse_bond(&mut self) -> PResult<Bond> {
        let b = self.bond_expr()?;
        if self.look != Tok::Eos {
            return Err(ParseFailure);
        }
        Ok(b)
    }
}

// ================================================================== drivers

/// `labelRecursivePatterns`: give every `$(...)` a serial number so that
/// identical recursive queries can share match results.
pub fn label_recursive_patterns(sma: &str) -> String {
    #[derive(PartialEq)]
    enum S {
        Base,
        Branch,
        Recurse,
    }
    let b = sma.as_bytes();
    let mut state = vec![S::Base];
    let mut start_recurse: Vec<usize> = Vec::new();
    let mut patterns: BTreeMap<Vec<u8>, String> = BTreeMap::new();
    let mut res: Vec<u8> = Vec::with_capacity(b.len() + 8);
    let mut pos = 0;
    while pos < b.len() {
        res.push(b[pos]);
        if b[pos] == b'$' && pos + 1 < b.len() && b[pos + 1] == b'(' {
            state.push(S::Recurse);
            start_recurse.push(pos);
            pos += 1;
            res.push(b[pos]);
        } else if b[pos] == b'(' {
            state.push(S::Branch);
        } else if b[pos] == b')' {
            if state.last().is_none_or(|s| *s == S::Base) {
                return sma.to_string();
            }
            if state.pop() == Some(S::Recurse) {
                let dollar = start_recurse.pop().unwrap_or(0);
                if pos + 1 >= b.len() || b[pos + 1] != b'_' {
                    let recurs = b[dollar..=pos].to_vec();
                    let n = patterns.len();
                    let label = patterns
                        .entry(recurs)
                        .or_insert_with(|| (n + 100).to_string())
                        .clone();
                    res.push(b'_');
                    res.extend_from_slice(label.as_bytes());
                }
            }
        }
        pos += 1;
    }
    String::from_utf8(res).unwrap_or_else(|_| sma.to_string())
}

/// `SmartsParserParams`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmartsParserParams {
    pub debug_parse: bool,
    pub allow_cx_smiles: bool,
    pub strict_cx_smiles: bool,
    pub parse_name: bool,
    pub merge_hs: bool,
    pub skip_cleanup: bool,
    pub replacements: BTreeMap<String, String>,
}

impl Default for SmartsParserParams {
    fn default() -> Self {
        SmartsParserParams {
            debug_parse: false,
            allow_cx_smiles: true,
            strict_cx_smiles: true,
            parse_name: true,
            merge_hs: false,
            skip_cleanup: false,
            replacements: BTreeMap::new(),
        }
    }
}

/// `setup_smarts_string`'s whitespace trimming: returns the trimmed input.
fn trim_ws(s: &str) -> &str {
    let b = s.as_bytes();
    let start = b.iter().position(|&c| c > 32).unwrap_or(b.len());
    let end = b.iter().rposition(|&c| c > 32).map_or(start, |e| e + 1);
    &s[start..end.max(start)]
}

fn to_mol(inp: &str, orig: &str) -> Option<RWMol> {
    if inp.is_empty() {
        return Some(RWMol::new());
    }
    let trimmed = trim_ws(inp);
    let fail = |e: &str| {
        crate::rd_log_error!("SMARTS Parse Error: {e} for input: '{orig}'");
        None
    };
    let mut parser = SmartsParser::new(trimmed);
    let parsed = parser.parse_mol();
    if parsed.is_ok() {
        for &(_, pos) in &parser.branch_points {
            ops::print_syntax_error_message(trimmed, "extra open parentheses", pos, "SMARTS");
        }
    }
    let mut mb = match parsed {
        Ok(m) if parser.branch_points.is_empty() => m,
        _ => return fail(&format!("Failed parsing SMARTS '{inp}'")),
    };
    let post = (|| -> std::result::Result<(), SmilesParseError> {
        ops::close_mol_rings(&mut mb.mol, &mut mb.ring_atoms, &mut mb.ring_bonds, false)?;
        ops::check_chirality_specifications(&mut mb.mol, true)?;
        ops::set_unspecified_bond_types(&mut mb.mol);
        ops::adjust_atom_chirality_flags(&mut mb.mol);
        Ok(())
    })();
    if let Err(SmilesParseError(e)) = post {
        return fail(&e);
    }
    if mb.mol.has_atom_bookmark(CI_RIGHTMOST_ATOM) {
        mb.mol.clear_atom_bookmark(CI_RIGHTMOST_ATOM);
    }
    Some(mb.mol)
}

/// `MolFromSmarts(smarts, params)`; `None` on a parse error (logged).
pub fn mol_from_smarts(smarts: &str, params: &SmartsParserParams) -> crate::Result<Option<RWMol>> {
    let sp = super::SmilesParserParams {
        allow_cx_smiles: params.allow_cx_smiles,
        strict_cx_smiles: params.strict_cx_smiles,
        parse_name: params.parse_name,
        replacements: params.replacements.clone(),
        ..Default::default()
    };
    let (lsmarts, mut name, cx_part) = super::preprocess(smarts, &sp);
    let Some(mut res) = to_mol(&label_recursive_patterns(&lsmarts), &lsmarts) else {
        return Ok(None);
    };
    super::handle_cx_part_and_name(
        &mut res,
        params.allow_cx_smiles,
        params.strict_cx_smiles,
        params.parse_name,
        &cx_part,
        &mut name,
    )?;
    if params.merge_hs {
        crate::graphmol::mol_ops::query_hs::merge_query_hs(&mut res, false, false)?;
    }
    crate::graphmol::chirality::bond_dirs::set_bond_stereo_from_directions(&mut res)?;
    if !params.skip_cleanup {
        ops::cleanup_after_parsing(&mut res);
    }
    if !name.is_empty() {
        res.props.set_val(cp::_Name, name);
    }
    Ok(Some(res))
}

/// `AtomFromSmarts`
pub fn atom_from_smarts(s: &str) -> Option<Atom> {
    if s.is_empty() {
        return None;
    }
    let trimmed = trim_ws(s);
    let mut p = SmartsParser::new(trimmed);
    match p.parse_atom() {
        Ok(a) => Some(a),
        Err(_) => {
            crate::rd_log_error!(
                "SMARTS Parse Error: Failed parsing SMARTS '{s}' for input: '{s}'"
            );
            None
        }
    }
}

/// `BondFromSmarts`
pub fn bond_from_smarts(s: &str) -> Option<Bond> {
    if s.is_empty() {
        return None;
    }
    let trimmed = trim_ws(s);
    let mut p = SmartsParser::new(trimmed);
    match p.parse_bond() {
        Ok(b) => Some(b),
        Err(_) => {
            crate::rd_log_error!(
                "SMARTS Parse Error: Failed parsing SMARTS '{s}' for input: '{s}'"
            );
            None
        }
    }
}
