//! SMILES grammar — port of RDKit `Code/GraphMol/SmilesParse/smiles.yy`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! A recursive-descent parser whose actions mirror the bison rules one to
//! one (atom/bond creation order, bookmarks, ring-closure bookkeeping), so
//! that atom and bond indices match upstream exactly.

use std::collections::BTreeMap;

use super::lexer::{Lexer, Tok};
use super::ops;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::romol::{CI_RIGHTMOST_ATOM, RWMol};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};

/// Parse failure; `message` has already been logged where upstream logs.
#[derive(Debug, Clone)]
pub struct ParseFailure;

/// Mutable parser state (bison's parse-params).
pub struct SmilesParser<'a> {
    pub(crate) input: &'a str,
    lex: Lexer<'a>,
    /// lookahead token and the scanner position after reading it
    look: Tok,
    look_end: usize,
    pub(crate) num_bonds_parsed: u32,
    pub(crate) branch_points: Vec<(u32, usize)>,
    /// ring-closure number -> atoms (upstream: atom bookmarks)
    pub(crate) ring_atoms: BTreeMap<i32, Vec<u32>>,
    /// ring-closure number -> partial bonds (upstream: bond bookmarks)
    pub(crate) ring_bonds: BTreeMap<i32, Vec<Bond>>,
}

impl std::fmt::Debug for SmilesParser<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SmilesParser")
            .field("input", &self.input)
            .finish()
    }
}

type PResult<T> = Result<T, ParseFailure>;

impl<'a> SmilesParser<'a> {
    pub fn new(input: &'a str) -> Self {
        let mut lex = Lexer::new(input);
        let look = lex.next_token();
        let look_end = lex.token_end;
        SmilesParser {
            input,
            lex,
            look,
            look_end,
            num_bonds_parsed: 0,
            branch_points: Vec::new(),
            ring_atoms: BTreeMap::new(),
            ring_bonds: BTreeMap::new(),
        }
    }

    fn bump(&mut self) -> Tok {
        let t = std::mem::replace(&mut self.look, self.lex.next_token());
        self.look_end = self.lex.token_end;
        t
    }

    /// `yyerror(..., msg)`: logs the syntax-error message at the position
    /// of the current lookahead token.
    pub(crate) fn error(&self, msg: &str) -> ParseFailure {
        ops::print_syntax_error_message(self.input, msg, self.look_end, "SMILES");
        ParseFailure
    }

    fn syntax_error(&self) -> ParseFailure {
        self.error("syntax error")
    }

    // ----------------------------------------------------------- numbers

    fn is_digit(t: &Tok) -> bool {
        matches!(t, Tok::Zero | Tok::NonZeroDigit(_))
    }

    fn digit(&mut self) -> PResult<i32> {
        let d = match self.look {
            Tok::Zero => 0,
            Tok::NonZeroDigit(d) => d,
            _ => return Err(self.syntax_error()),
        };
        self.bump();
        Ok(d)
    }

    /// `number: ZERO_TOKEN | nonzero_number`
    fn number(&mut self) -> PResult<i32> {
        match self.look {
            Tok::Zero => {
                self.bump();
                Ok(0)
            }
            Tok::NonZeroDigit(d) => {
                self.bump();
                let mut v: i32 = d;
                while Self::is_digit(&self.look) {
                    let dg = match self.look {
                        Tok::NonZeroDigit(x) => x,
                        _ => 0,
                    };
                    if v >= i32::MAX / 10 || v * 10 >= i32::MAX - dg {
                        return Err(self.error("number too large"));
                    }
                    self.bump();
                    v = v * 10 + dg;
                }
                Ok(v)
            }
            _ => Err(self.syntax_error()),
        }
    }

    /// `ring_number`
    fn ring_number(&mut self) -> PResult<i32> {
        if Self::is_digit(&self.look) {
            return self.digit();
        }
        // PERCENT
        self.bump();
        match self.look {
            Tok::NonZeroDigit(d) => {
                self.bump();
                Ok(d * 10 + self.digit()?)
            }
            Tok::GroupOpen => {
                self.bump();
                let mut v = 0;
                let mut n = 0;
                while Self::is_digit(&self.look) && n < 5 {
                    v = v * 10 + self.digit()?;
                    n += 1;
                }
                if n == 0 || self.look != Tok::GroupClose {
                    return Err(self.syntax_error());
                }
                self.bump();
                Ok(v)
            }
            _ => Err(self.syntax_error()),
        }
    }

    // ----------------------------------------------------------- atoms

    fn simple_atom(tok: &Tok) -> Option<Atom> {
        match *tok {
            Tok::Organic(n) => {
                let mut a = Atom::new(n);
                if n == 0 {
                    a.props.set_val(cp::dummyLabel, "*");
                }
                Some(a)
            }
            Tok::Aromatic(n) => {
                let mut a = Atom::new(n);
                a.set_is_aromatic(true);
                Some(a)
            }
            _ => None,
        }
    }

    fn starts_atom(t: &Tok) -> bool {
        matches!(t, Tok::Organic(_) | Tok::Aromatic(_) | Tok::AtomOpen)
    }

    /// `element`
    fn element(&mut self) -> PResult<Atom> {
        let isotope = if Self::is_digit(&self.look) {
            Some(self.number()?)
        } else {
            None
        };
        let mut atom = match self.look {
            Tok::Organic(_) | Tok::Aromatic(_) => {
                let tok = self.bump();
                Self::simple_atom(&tok).expect("simple atom")
            }
            Tok::Atom(n) => {
                self.bump();
                Atom::new(n)
            }
            Tok::Hash => {
                self.bump();
                Atom::new(self.number()? as u32)
            }
            _ => return Err(self.syntax_error()),
        };
        if let Some(iso) = isotope {
            atom.set_isotope(iso as u32);
        }
        Ok(atom)
    }

    /// `chiral_element`
    fn chiral_element(&mut self) -> PResult<Atom> {
        let mut atom = self.element()?;
        match self.look {
            Tok::At => {
                self.bump();
                if self.look == Tok::At {
                    self.bump();
                    atom.set_chiral_tag(ChiralType::TetrahedralCW);
                } else {
                    atom.set_chiral_tag(ChiralType::TetrahedralCCW);
                }
            }
            Tok::ChiClass(t) => {
                self.bump();
                if Self::is_digit(&self.look) {
                    let n = self.number()?;
                    if n == 0 {
                        return Err(self.error("chiral permutation cannot be zero"));
                    }
                    atom.set_chiral_tag(t);
                    atom.props.set_val(cp::_chiralPermutation, n);
                } else {
                    atom.set_chiral_tag(t);
                    atom.props.set_val(cp::_chiralPermutation, 0i32);
                }
            }
            _ => {}
        }
        Ok(atom)
    }

    /// `h_element`
    fn h_element(&mut self) -> PResult<Atom> {
        // forms starting with an (isotope) number followed by H
        let save_is_h = |p: &Self| p.look == Tok::H;
        if save_is_h(self) {
            self.bump();
            let mut a = Atom::new(1);
            if self.look == Tok::H {
                self.bump();
                if Self::is_digit(&self.look) {
                    let n = self.number()?;
                    a.set_num_explicit_hs(n as u32);
                } else {
                    a.set_num_explicit_hs(1);
                }
            }
            return Ok(a);
        }
        if Self::is_digit(&self.look) {
            // need to know whether the number is followed by H
            let (iso, after) = self.peek_number_then_token();
            if after == Tok::H {
                let iso = iso?;
                self.number()?;
                self.bump(); // H
                let mut a = Atom::new(1);
                a.set_isotope(iso as u32);
                if self.look == Tok::H {
                    self.bump();
                    if Self::is_digit(&self.look) {
                        let n = self.number()?;
                        a.set_num_explicit_hs(n as u32);
                    } else {
                        a.set_num_explicit_hs(1);
                    }
                }
                return Ok(a);
            }
        }
        let mut atom = self.chiral_element()?;
        if self.look == Tok::H {
            self.bump();
            if Self::is_digit(&self.look) {
                let n = self.number()?;
                atom.set_num_explicit_hs(n as u32);
            } else {
                atom.set_num_explicit_hs(1);
            }
        }
        Ok(atom)
    }

    /// Look past a run of digits: (value, following token).
    fn peek_number_then_token(&self) -> (PResult<i32>, Tok) {
        let mut l = self.lex.clone();
        let mut v: i64 = 0;
        let mut t = self.look.clone();
        while let Tok::Zero | Tok::NonZeroDigit(_) = t {
            let d = if let Tok::NonZeroDigit(d) = t { d } else { 0 };
            v = (v * 10 + i64::from(d)).min(i64::from(i32::MAX));
            t = l.next_token();
        }
        (Ok(v as i32), t)
    }

    /// `charge_element`
    fn charge_element(&mut self) -> PResult<Atom> {
        let mut atom = self.h_element()?;
        match self.look {
            Tok::Plus => {
                self.bump();
                if self.look == Tok::Plus {
                    self.bump();
                    atom.set_formal_charge(2);
                } else if Self::is_digit(&self.look) {
                    let n = self.number()?;
                    atom.set_formal_charge(n);
                } else {
                    atom.set_formal_charge(1);
                }
            }
            Tok::Minus => {
                self.bump();
                if self.look == Tok::Minus {
                    self.bump();
                    atom.set_formal_charge(-2);
                } else if Self::is_digit(&self.look) {
                    let n = self.number()?;
                    atom.set_formal_charge(-n);
                } else {
                    atom.set_formal_charge(-1);
                }
            }
            _ => {}
        }
        Ok(atom)
    }

    /// `atomd`
    fn atomd(&mut self) -> PResult<Atom> {
        if let Some(a) = Self::simple_atom(&self.look) {
            self.bump();
            return Ok(a);
        }
        if self.look != Tok::AtomOpen {
            return Err(self.syntax_error());
        }
        self.bump();
        let mut atom = self.charge_element()?;
        if self.look == Tok::Colon {
            self.bump();
            let n = self.number()?;
            atom.set_no_implicit(true);
            atom.props.set_val(cp::molAtomMapNumber, n);
        } else {
            atom.set_no_implicit(true);
        }
        if self.look != Tok::AtomClose {
            return Err(self.syntax_error());
        }
        self.bump();
        Ok(atom)
    }

    // ----------------------------------------------------------- bonds

    fn bond_from_token(&self, t: &Tok) -> Option<Bond> {
        match *t {
            Tok::Bond {
                bond_type,
                dir,
                unspecified_order,
                aromatic,
            } => {
                let mut b = Bond::new(bond_type);
                b.set_bond_dir(dir);
                b.set_is_aromatic(aromatic);
                if unspecified_order {
                    b.props.set_val(cp::_unspecifiedOrder, 1i32);
                }
                Some(b)
            }
            _ => None,
        }
    }

    fn active_atom(mol: &RWMol) -> u32 {
        mol.get_atom_with_bookmark(CI_RIGHTMOST_ATOM)
            .unwrap_or_else(|_| mol.num_atoms().saturating_sub(1))
    }

    /// `mol BOND_TOKEN atomd` (and the branch variant)
    fn add_explicit_bond(&mut self, mol: &mut RWMol, mut bond: Bond, atom: Atom) -> PResult<()> {
        let a1 = Self::active_atom(mol);
        let a2 = mol.add_atom(atom, true);
        match bond.bond_type() {
            BondType::DativeR => {
                bond.begin = a1;
                bond.end = a2;
                bond.set_bond_type(BondType::Dative);
            }
            BondType::DativeL => {
                bond.begin = a2;
                bond.end = a1;
                bond.set_bond_type(BondType::Dative);
            }
            _ => {
                bond.begin = a1;
                bond.end = a2;
            }
        }
        self.num_bonds_parsed += 1;
        mol.add_bond_object(bond).map_err(|_| self.syntax_error())?;
        Ok(())
    }

    fn add_implicit_bond(
        &mut self,
        mol: &mut RWMol,
        atom: Atom,
        explicit_single: bool,
    ) -> PResult<()> {
        let a1 = Self::active_atom(mol);
        let a2 = mol.add_atom(atom, true);
        let bt = if explicit_single {
            BondType::Single
        } else {
            ops::get_unspecified_bond_type(mol.atom(a1), mol.atom(a2))
        };
        mol.add_bond(a1, a2, bt).map_err(|_| self.syntax_error())?;
        self.num_bonds_parsed += 1;
        Ok(())
    }

    /// The ring-closure actions of `mol ring_number`, `mol BOND ring_number`,
    /// `mol MINUS ring_number`.
    fn ring_closure(&mut self, mol: &mut RWMol, num: i32, mut partial: Bond) {
        let atom = Self::active_atom(mol);
        self.ring_atoms.entry(num).or_default().push(atom);
        partial.begin = atom;
        let list = self.ring_bonds.entry(num).or_default();
        list.push(partial);
        if list.len() % 2 == 0 {
            let b = list.last_mut().expect("just pushed");
            b.props.set_val("_cxsmilesBondIdx", self.num_bonds_parsed);
            self.num_bonds_parsed += 1;
        }
        ops::check_ring_closure_branch_status(mol, atom);
        let mut closures = match mol.atom(atom).props.get(cp::_RingClosures) {
            Some(RDValue::VecInt(v)) => v.clone(),
            _ => Vec::new(),
        };
        closures.push(-(num + 1));
        mol.atom_mut(atom)
            .props
            .set_val(cp::_RingClosures, closures);
    }

    fn starts_ring_number(t: &Tok) -> bool {
        matches!(t, Tok::Zero | Tok::NonZeroDigit(_) | Tok::Percent)
    }

    // ----------------------------------------------------------- molecule

    /// `meta_start: START_MOL mol` followed by end of input.
    pub fn parse_mol(&mut self) -> PResult<RWMol> {
        let mut mol = RWMol::new();
        // first atom
        let mut first = self.atomd()?;
        first.props.set_val(cp::_SmilesStart, 1i32);
        mol.add_atom(first, true);
        loop {
            let t = self.look.clone();
            match t {
                Tok::Eos => break,
                ref t if Self::starts_atom(t) => {
                    let a = self.atomd()?;
                    self.add_implicit_bond(&mut mol, a, false)?;
                }
                Tok::Bond { .. } => {
                    self.bump();
                    let bond = self.bond_from_token(&t).expect("bond token");
                    if Self::starts_atom(&self.look) {
                        let a = self.atomd()?;
                        self.add_explicit_bond(&mut mol, bond, a)?;
                    } else if Self::starts_ring_number(&self.look) {
                        let num = self.ring_number()?;
                        let mut partial = Bond::new(bond.bond_type());
                        if bond.has_prop(cp::_unspecifiedOrder) {
                            partial.props.set_val(cp::_unspecifiedOrder, 1i32);
                        }
                        partial.set_bond_dir(bond.bond_dir());
                        self.ring_closure(&mut mol, num, partial);
                    } else {
                        return Err(self.syntax_error());
                    }
                }
                Tok::AnyBond => {
                    // query bonds arrive with the SMARTS/query machinery (stage 4)
                    return Err(self.syntax_error());
                }
                Tok::Minus => {
                    self.bump();
                    if Self::starts_atom(&self.look) {
                        let a = self.atomd()?;
                        self.add_implicit_bond(&mut mol, a, true)?;
                    } else if Self::starts_ring_number(&self.look) {
                        let num = self.ring_number()?;
                        self.ring_closure(&mut mol, num, Bond::new(BondType::Single));
                    } else {
                        return Err(self.syntax_error());
                    }
                }
                Tok::Separator => {
                    self.bump();
                    let mut a = self.atomd()?;
                    a.set_prop(cp::_SmilesStart, 1i32, true).expect("valid key");
                    mol.add_atom(a, true);
                }
                Tok::Zero | Tok::NonZeroDigit(_) | Tok::Percent => {
                    let num = self.ring_number()?;
                    let mut partial = Bond::new(BondType::Unspecified);
                    partial.props.set_val(cp::_unspecifiedOrder, 1i32);
                    self.ring_closure(&mut mol, num, partial);
                }
                Tok::GroupOpen => {
                    // branch_open_token: position at the end of '('
                    let open_pos = self.look_end;
                    self.bump();
                    let branch_atom = Self::active_atom(&mol);
                    let t2 = self.look.clone();
                    if Self::starts_atom(&t2) {
                        let a = self.atomd()?;
                        self.add_implicit_bond(&mut mol, a, false)?;
                    } else if let Tok::Bond { .. } = t2 {
                        self.bump();
                        let bond = self.bond_from_token(&t2).expect("bond token");
                        let a = self.atomd()?;
                        self.add_explicit_bond(&mut mol, bond, a)?;
                    } else if t2 == Tok::Minus {
                        self.bump();
                        let a = self.atomd()?;
                        self.add_implicit_bond(&mut mol, a, true)?;
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
                    mol.replace_atom_bookmark(atom, CI_RIGHTMOST_ATOM);
                }
                _ => return Err(self.syntax_error()),
            }
        }
        Ok(mol)
    }
}
