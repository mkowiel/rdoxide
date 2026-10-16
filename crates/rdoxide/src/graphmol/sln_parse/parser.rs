//! The SLN parser: a recursive descent version of the bison grammar in
//! `sln.yy` (2026.09.1) that performs the same actions in the same order.
//!
//! Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::attribs::{
    AttribCombineOp, AttribListType, AttribType, parse_atom_attribs, parse_bond_attribs,
    parse_mol_attribs,
};
use super::lexer::{Lexer, Tok};
use super::ops::{
    MolList, add_atom_to_mol, add_atom_to_mol_single, add_branch_to_mol, add_frag_to_mol,
    close_ring_bond, start_mol,
};
use super::{SlnError, finalize_query_mol};
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::queries::{
    CompositeQueryType, Query, QueryKind, make_bond_null_query, make_recursive_structure_query,
};
use crate::graphmol::query_ops::{bond_expand_query, query_atom_from_atom, query_bond_from_bond};
use crate::rdgeneral::{RDProps, common_properties as cp};

type PResult<T> = std::result::Result<T, SlnError>;

pub struct Parser {
    lexer: Lexer,
    peeked: Option<Tok>,
    pub list: MolList,
    do_queries: bool,
}

fn is_one_bond(t: &Tok) -> bool {
    matches!(t, Tok::Minus | Tok::Equals | Tok::Hash | Tok::Colon)
}

fn is_bond(t: &Tok) -> bool {
    is_one_bond(t) || matches!(t, Tok::Tilde)
}

fn is_atom_start(t: &Tok) -> bool {
    matches!(t, Tok::Atom(_) | Tok::HAsterix | Tok::HBracket)
}

impl Parser {
    pub fn new(input: &str, do_queries: bool) -> Self {
        Parser {
            lexer: Lexer::new(input, do_queries),
            peeked: None,
            list: Vec::new(),
            do_queries,
        }
    }

    fn peek(&mut self) -> PResult<&Tok> {
        if self.peeked.is_none() {
            self.peeked = Some(self.lexer.next_token()?);
        }
        Ok(self.peeked.as_ref().expect("just filled"))
    }

    fn next(&mut self) -> PResult<Tok> {
        self.peek()?;
        Ok(self.peeked.take().expect("just filled"))
    }

    fn expect(&mut self, ok: fn(&Tok) -> bool) -> PResult<Tok> {
        let t = self.next()?;
        if ok(&t) { Ok(t) } else { Err(SlnError::Syntax) }
    }

    fn h_atom(&self) -> PResult<Atom> {
        let a = Atom::new(1);
        if self.do_queries {
            Ok(query_atom_from_atom(&a)?)
        } else {
            Ok(a)
        }
    }

    fn new_bond(&self, t: BondType) -> Bond {
        let b = Bond::new(t);
        if self.do_queries {
            query_bond_from_bond(&b)
        } else {
            b
        }
    }

    /// `cmpd`: the start symbol.
    pub fn cmpd(&mut self) -> PResult<usize> {
        let idx = self.mol()?;
        loop {
            match self.peek()? {
                Tok::OpenAngle => {
                    self.next()?;
                    let attribs = self.ctab_attrib_list()?;
                    self.expect(|t| matches!(t, Tok::CloseAngle))?;
                    // allowing mol<attrs><attrs> seems to be a NIBR thing, I
                    // don't think it's standard SLN (upstream)
                    let m = self.list[idx].as_mut().expect("null molecule");
                    parse_mol_attribs(&mut m.mol, attribs)?;
                }
                Tok::Eos => return Err(SlnError::Accept),
                _ => return Ok(idx),
            }
        }
    }

    fn mol(&mut self) -> PResult<usize> {
        let mut m = self.primmol()?;
        while matches!(self.peek()?, Tok::Separator) {
            self.next()?;
            let p = self.primmol()?;
            m = add_frag_to_mol(&mut self.list, m, p)?;
        }
        Ok(m)
    }

    fn number(&mut self) -> PResult<i32> {
        let Tok::Digit(mut n) = self.expect(|t| matches!(t, Tok::Digit(_)))? else {
            unreachable!()
        };
        while let Tok::Digit(d) = self.peek()? {
            let d = *d;
            self.next()?;
            n = n.wrapping_mul(10).wrapping_add(d);
        }
        Ok(n)
    }

    /// An optional hydrogen atom after a ring closure or a branch.
    fn optional_h(&mut self, idx: usize) -> PResult<()> {
        if matches!(self.peek()?, Tok::H) {
            self.next()?;
            let a = self.h_atom()?;
            add_atom_to_mol_single(&mut self.list, idx, a, self.do_queries)?;
        }
        Ok(())
    }

    fn primmol(&mut self) -> PResult<usize> {
        let idx = match self.peek()? {
            Tok::H => {
                self.next()?;
                let a = self.h_atom()?;
                start_mol(&mut self.list, a, self.do_queries)?
            }
            t if is_atom_start(t) => {
                let a = self.atom()?;
                start_mol(&mut self.list, a, self.do_queries)?
            }
            _ => return Err(SlnError::Syntax),
        };
        let dq = self.do_queries;
        loop {
            let t = self.peek()?;
            if is_atom_start(t) {
                let a = self.atom()?;
                add_atom_to_mol_single(&mut self.list, idx, a, dq)?;
            } else if is_bond(t) {
                let bond = self.bond()?;
                match self.peek()? {
                    Tok::At => {
                        self.next()?;
                        let n = self.number()?;
                        close_ring_bond(&mut self.list, idx, n, bond)?;
                        self.optional_h(idx)?;
                    }
                    Tok::H => {
                        self.next()?;
                        let a = self.h_atom()?;
                        add_atom_to_mol(&mut self.list, idx, a, bond, dq)?;
                    }
                    t if is_atom_start(t) => {
                        let a = self.atom()?;
                        add_atom_to_mol(&mut self.list, idx, a, bond, dq)?;
                    }
                    _ => return Err(SlnError::Syntax),
                }
            } else if matches!(t, Tok::At) {
                self.next()?;
                let n = self.number()?;
                close_ring_bond(&mut self.list, idx, n, Bond::new(BondType::Single))?;
                self.optional_h(idx)?;
            } else if matches!(t, Tok::OpenParen) {
                self.next()?;
                let t = self.peek()?;
                if matches!(t, Tok::At) {
                    self.next()?;
                    let n = self.number()?;
                    self.expect(|t| matches!(t, Tok::CloseParen))?;
                    close_ring_bond(&mut self.list, idx, n, Bond::new(BondType::Single))?;
                } else if is_bond(t) {
                    let bond = self.bond()?;
                    if matches!(self.peek()?, Tok::At) {
                        self.next()?;
                        let n = self.number()?;
                        self.expect(|t| matches!(t, Tok::CloseParen))?;
                        close_ring_bond(&mut self.list, idx, n, bond)?;
                    } else {
                        let branch = self.primmol()?;
                        self.expect(|t| matches!(t, Tok::CloseParen))?;
                        add_branch_to_mol(&mut self.list, idx, branch, bond)?;
                    }
                } else {
                    let branch = self.primmol()?;
                    self.expect(|t| matches!(t, Tok::CloseParen))?;
                    add_branch_to_mol(&mut self.list, idx, branch, Bond::new(BondType::Single))?;
                }
                self.optional_h(idx)?;
            } else {
                return Ok(idx);
            }
        }
    }

    /// `atom`
    fn atom(&mut self) -> PResult<Atom> {
        match self.next()? {
            Tok::HAsterix => {
                let mut a = self.h_atom()?;
                a.set_prop(cp::_starred, 1i32, true)?;
                Ok(a)
            }
            Tok::HBracket => self.hatom_bracket(),
            Tok::Atom(a) => {
                let mut a = self.primatom(a)?;
                if matches!(self.peek()?, Tok::H) {
                    self.next()?;
                    a.set_num_explicit_hs(1);
                    if let Tok::Digit(d) = self.peek()? {
                        let d = *d;
                        self.next()?;
                        a.set_num_explicit_hs(d as u32);
                    }
                }
                Ok(a)
            }
            _ => Err(SlnError::Syntax),
        }
    }

    /// The part of `hatom` after `H[`
    fn hatom_bracket(&mut self) -> PResult<Atom> {
        let mut a = self.h_atom()?;
        if matches!(self.peek()?, Tok::Digit(_)) {
            let n = self.number()?;
            a.set_prop(cp::_AtomID, n as u32, false)?;
            if matches!(self.peek()?, Tok::Colon) {
                self.next()?;
                let attribs = self.attrib_list()?;
                self.expect(|t| matches!(t, Tok::CloseBracket))?;
                parse_atom_attribs(&mut a, attribs, self.do_queries)?;
            } else {
                self.expect(|t| matches!(t, Tok::CloseBracket))?;
            }
        } else {
            let attribs = self.attrib_list()?;
            self.expect(|t| matches!(t, Tok::CloseBracket))?;
            parse_atom_attribs(&mut a, attribs, self.do_queries)?;
        }
        Ok(a)
    }

    /// `primatom` after the atom token
    fn primatom(&mut self, mut a: Atom) -> PResult<Atom> {
        loop {
            match self.peek()? {
                Tok::Asterix => {
                    self.next()?;
                    a.set_prop(cp::_starred, 1i32, true)?;
                }
                Tok::OpenBracket => {
                    self.next()?;
                    if matches!(self.peek()?, Tok::Digit(_)) {
                        let n = self.number()?;
                        a.set_prop(cp::_AtomID, n as u32, false)?;
                        match self.next()? {
                            Tok::Colon => {
                                let attribs = self.attrib_list()?;
                                self.expect(|t| matches!(t, Tok::CloseBracket))?;
                                parse_atom_attribs(&mut a, attribs, self.do_queries)?;
                            }
                            Tok::CloseBracket => {}
                            _ => return Err(SlnError::Syntax),
                        }
                    } else {
                        let attribs = self.attrib_list()?;
                        self.expect(|t| matches!(t, Tok::CloseBracket))?;
                        parse_atom_attribs(&mut a, attribs, self.do_queries)?;
                    }
                }
                _ => return Ok(a),
            }
        }
    }

    /// `bond`
    fn bond(&mut self) -> PResult<Bond> {
        let t = self.next()?;
        let mut bond = if matches!(t, Tok::Tilde) {
            // tildes can't be mixed with regular bonds in expressions
            let mut b = Bond::new(BondType::Unspecified);
            b.set_query(Some(make_bond_null_query()));
            b
        } else {
            let mut bond = self.new_bond(Self::bond_type(&t)?);
            while is_one_bond(self.peek()?) {
                if !self.do_queries {
                    return Err(SlnError::Abort(
                        "sequential bonds not allowed in non-queries".into(),
                    ));
                }
                let t2 = self.next()?;
                let b2 = self.new_bond(Self::bond_type(&t2)?);
                let q2 = b2.query().expect("query bond").clone();
                bond_expand_query(&mut bond, q2, CompositeQueryType::Or, true)?;
            }
            bond
        };
        if matches!(self.peek()?, Tok::OpenBracket) {
            self.next()?;
            let attribs = self.attrib_list()?;
            self.expect(|t| matches!(t, Tok::CloseBracket))?;
            parse_bond_attribs(&mut bond, attribs, self.do_queries)?;
        }
        Ok(bond)
    }

    fn bond_type(t: &Tok) -> PResult<BondType> {
        Ok(match t {
            Tok::Minus => BondType::Single,
            Tok::Equals => BondType::Double,
            Tok::Hash => BondType::Triple,
            Tok::Colon => BondType::Aromatic,
            _ => return Err(SlnError::Syntax),
        })
    }

    /// `attriblist`
    fn attrib_list(&mut self) -> PResult<AttribListType> {
        let first = self.attrib()?;
        let mut list = vec![(AttribCombineOp::LowPriAnd, first)];
        loop {
            let op = match self.peek()? {
                Tok::And => AttribCombineOp::And,
                Tok::Or => AttribCombineOp::Or,
                Tok::Semi => AttribCombineOp::LowPriAnd,
                _ => return Ok(list),
            };
            self.next()?;
            let a = self.attrib()?;
            list.push((op, a));
        }
    }

    /// `ctabattriblist`
    fn ctab_attrib_list(&mut self) -> PResult<AttribListType> {
        let mut list = vec![(AttribCombineOp::And, self.ctab_attrib()?)];
        while matches!(self.peek()?, Tok::Semi) {
            self.next()?;
            let a = self.ctab_attrib()?;
            list.push((AttribCombineOp::And, a));
        }
        Ok(list)
    }

    fn ctab_attrib(&mut self) -> PResult<AttribType> {
        let Tok::Text(first) = self.expect(|t| matches!(t, Tok::Text(_)))? else {
            unreachable!()
        };
        let op = match self.peek()? {
            Tok::Equals => "=",
            Tok::ColonEquals => ":=",
            Tok::CaretEquals => "^=",
            _ => {
                return Ok(AttribType {
                    first: first.to_lowercase(),
                    ..Default::default()
                });
            }
        };
        self.next()?;
        let Tok::Text(second) = self.expect(|t| matches!(t, Tok::Text(_)))? else {
            unreachable!()
        };
        Ok(AttribType {
            first: first.to_lowercase(),
            op: op.into(),
            second: second.to_lowercase(),
            ..Default::default()
        })
    }

    fn simple_attrib(first: &str, op: &str, second: &str) -> AttribType {
        AttribType {
            first: first.into(),
            op: op.into(),
            second: second.into(),
            ..Default::default()
        }
    }

    /// `attrib`
    fn attrib(&mut self) -> PResult<AttribType> {
        match self.next()? {
            Tok::Text(first) => {
                if let Tok::Compare(_) = self.peek()? {
                    let Tok::Compare(op) = self.next()? else {
                        unreachable!()
                    };
                    let Tok::Text(second) = self.expect(|t| matches!(t, Tok::Text(_)))? else {
                        unreachable!()
                    };
                    Ok(Self::simple_attrib(
                        &first.to_lowercase(),
                        &op,
                        &second.to_lowercase(),
                    ))
                } else {
                    Ok(Self::simple_attrib(&first.to_lowercase(), "", ""))
                }
            }
            Tok::Not => {
                let mut a = self.attrib()?;
                a.negated = true;
                Ok(a)
            }
            Tok::Plus => {
                if let Tok::Digit(d) = self.peek()? {
                    let d = *d;
                    self.next()?;
                    Ok(Self::simple_attrib("charge", "=", &d.to_string()))
                } else {
                    Ok(Self::simple_attrib("charge", "=", "+1"))
                }
            }
            Tok::Minus => {
                if let Tok::Digit(d) = self.peek()? {
                    let d = *d;
                    self.next()?;
                    Ok(Self::simple_attrib(
                        "charge",
                        "=",
                        &d.wrapping_neg().to_string(),
                    ))
                } else {
                    Ok(Self::simple_attrib("charge", "=", "-1"))
                }
            }
            Tok::Asterix => Ok(Self::simple_attrib("spin", "=", "d")),
            t @ (Tok::Recurse | Tok::NegRecurse) => {
                let negate = matches!(t, Tok::NegRecurse);
                let idx = self.cmpd()?;
                let mut orq: Query<_> = Query::new(QueryKind::Or);
                orq.negate = negate;
                orq.children.push(self.recursive_child(idx)?);
                while matches!(self.peek()?, Tok::Comma) {
                    self.next()?;
                    let idx = self.cmpd()?;
                    orq.children.push(self.recursive_child(idx)?);
                }
                let mut a = Self::simple_attrib("is", "=", "");
                a.struct_query = Some(orq);
                Ok(a)
            }
            _ => Err(SlnError::Syntax),
        }
    }

    /// The recursive structure query of the molecule at `idx`, which is
    /// removed from the list.
    fn recursive_child(&mut self, idx: usize) -> PResult<crate::graphmol::queries::AtomQuery> {
        let sz = self.list.len();
        let mut mol = self.list[idx].take().expect("null molecule").mol;
        self.list.truncate(sz - 1);
        finalize_query_mol(&mut mol, true)?;
        Ok(make_recursive_structure_query(mol, 0))
    }
}
