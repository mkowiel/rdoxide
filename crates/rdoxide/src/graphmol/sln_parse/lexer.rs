//! The SLN scanner: a port of the flex rules in `sln.ll` (2026.09.1) with
//! their start conditions and "longest match, first rule wins" semantics.
//!
//! Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::graphmol::atom::Atom;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::{
    CompositeQueryType, expand_query, make_atom_null_query, make_atom_num_query,
};
use crate::graphmol::query_ops::query_atom_from_atom;

/// A token of the SLN grammar.
#[derive(Debug, Clone)]
pub enum Tok {
    /// Token 0: the end of the input (a newline)
    End,
    Eos,
    Text(String),
    /// Any other character
    Char,
    Digit(i32),
    H,
    HBracket,
    HAsterix,
    At,
    Atom(Atom),
    Compare(String),
    OpenParen,
    CloseParen,
    OpenBracket,
    CloseBracket,
    OpenAngle,
    CloseAngle,
    Separator,
    Asterix,
    Plus,
    Minus,
    Colon,
    Equals,
    Tilde,
    Hash,
    Comma,
    Not,
    And,
    Or,
    Semi,
    CaretEquals,
    ColonEquals,
    Recurse,
    NegRecurse,
}

/// The start conditions of the scanner.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Initial,
    SlnParam,
    CtabName,
    CtabVal,
    PropVal,
    Recurse,
}

/// What a matching rule does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Act {
    Text,
    AtomSymbol,
    Any,
    Hev,
    Hal,
    Het,
    HBracket,
    HAsterix,
    H,
    Recurse,
    NegRecurse,
    Minus,
    Plus,
    Hash,
    Tilde,
    Colon,
    OpenParen,
    CloseParen,
    Compare,
    CtabEquals,
    CtabColonEquals,
    CtabCaretEquals,
    Equals,
    /// `;`, `&` or `|` that leave the recursion / property value state
    PopThen(Op),
    /// `;`, `&` in the other states
    Plain(Op),
    /// `;` and `&` in the ctab value state
    CtabValOp(Op),
    OpenBracket,
    RecurseClose,
    PropValClose,
    SlnParamClose,
    OpenAngle,
    PopAngle,
    PropValAngle,
    Separator,
    Comma,
    At,
    Asterix,
    Digit,
    Newline,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Op {
    Semi,
    And,
    Or,
    Not,
}

impl Op {
    fn tok(self) -> Tok {
        match self {
            Op::Semi => Tok::Semi,
            Op::And => Tok::And,
            Op::Or => Tok::Or,
            Op::Not => Tok::Not,
        }
    }
}

const TWO_LETTER: &[&str] = &[
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Ar", "Ca", "Sc", "Ti", "Cr", "Mn", "Co", "Fe",
    "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Kr", "Rb", "Sr", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Cl", "Br",
];
const ONE_LETTER: &[u8] = b"KVWYUBCNOPSFI";

fn is_alpha(c: u8) -> bool {
    c.is_ascii_alphabetic()
}
fn is_alnum(c: u8) -> bool {
    c.is_ascii_alphanumeric()
}

/// The length of the longest prefix made of one or more `first` chars
/// followed by any number of `rest` chars.
fn run(s: &[u8], first: impl Fn(u8) -> bool, rest: impl Fn(u8) -> bool) -> usize {
    let mut i = 0;
    while i < s.len() && first(s[i]) {
        i += 1;
    }
    if i == 0 {
        return 0;
    }
    while i < s.len() && rest(s[i]) {
        i += 1;
    }
    i
}

fn starts(s: &[u8], p: &str) -> bool {
    s.starts_with(p.as_bytes())
}

fn starts_nocase(s: &[u8], p: &str) -> bool {
    s.len() >= p.len() && s[..p.len()].eq_ignore_ascii_case(p.as_bytes())
}

pub struct Lexer {
    input: Vec<u8>,
    pos: usize,
    state: State,
    stack: Vec<State>,
    do_queries: bool,
}

impl Lexer {
    pub fn new(input: &str, do_queries: bool) -> Self {
        // (the scanner reads a C string)
        let bytes = input.as_bytes();
        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        Lexer {
            input: bytes[..end].to_vec(),
            pos: 0,
            state: State::Initial,
            stack: Vec::new(),
            do_queries,
        }
    }

    fn push(&mut self, s: State) {
        self.stack.push(self.state);
        self.state = s;
    }

    fn pop(&mut self) {
        if let Some(s) = self.stack.pop() {
            self.state = s;
        }
    }

    /// The matching rules at the current position in the order of the flex
    /// file: `(length, action)`.
    fn candidates(&self, s: &[u8]) -> Vec<(usize, Act)> {
        use State::*;
        let st = self.state;
        let mut c: Vec<(usize, Act)> = Vec::new();
        let mut add = |len: usize, act: Act| {
            if len > 0 {
                c.push((len, act));
            }
        };
        let first = s[0];
        if st == PropVal {
            let n = s.iter().take_while(|&&b| !b";]>&|!".contains(&b)).count();
            add(n, Act::Text);
        }
        if st == SlnParam {
            add(
                run(s, is_alpha, |b| is_alnum(b) || b"_-,".contains(&b)),
                Act::Text,
            );
        }
        if st == CtabVal {
            let mut i = usize::from(first == b'"');
            let body = s[i..]
                .iter()
                .take_while(|&&b| is_alnum(b) || b"_-, .()".contains(&b))
                .count();
            if body > 0 {
                i += body;
                if i < s.len() && s[i] == b'"' {
                    i += 1;
                }
                add(i, Act::Text);
            }
        }
        if st == CtabName {
            add(
                run(s, is_alpha, |b| is_alnum(b) || b == b'_' || b == b'.'),
                Act::Text,
            );
        }
        if matches!(st, Initial | Recurse) {
            let mut best = 0;
            for sym in TWO_LETTER {
                if starts(s, sym) {
                    best = best.max(2);
                }
            }
            if ONE_LETTER.contains(&first) {
                best = best.max(1);
            }
            add(best, Act::AtomSymbol);
            if starts(s, "Any") {
                add(3, Act::Any);
            }
            if starts(s, "Hev") {
                add(3, Act::Hev);
            }
            if starts(s, "Hal") {
                add(3, Act::Hal);
            }
            if starts(s, "Het") {
                add(3, Act::Het);
            }
            if starts(s, "H[") {
                add(2, Act::HBracket);
            }
            if starts(s, "H*") {
                add(2, Act::HAsterix);
            }
            if first == b'H' {
                add(1, Act::H);
            }
        }
        if st == SlnParam {
            if starts_nocase(s, "is=") {
                add(3, Act::Recurse);
            }
            if starts_nocase(s, "not=") {
                add(4, Act::NegRecurse);
            }
        }
        match first {
            b'-' => add(1, Act::Minus),
            b'+' => add(1, Act::Plus),
            b'#' => add(1, Act::Hash),
            b'~' => add(1, Act::Tilde),
            b':' => add(1, Act::Colon),
            b'(' => add(1, Act::OpenParen),
            b')' => add(1, Act::CloseParen),
            _ => {}
        }
        if st == SlnParam {
            for op in ["!=", ">=", "<=", "=", ">", "<"] {
                if starts(s, op) {
                    add(op.len(), Act::Compare);
                    break;
                }
            }
        }
        if st == CtabName {
            if first == b'=' {
                add(1, Act::CtabEquals);
            }
            if starts(s, ":=") {
                add(2, Act::CtabColonEquals);
            }
            if starts(s, "^=") {
                add(2, Act::CtabCaretEquals);
            }
        }
        if matches!(st, Initial | Recurse) && first == b'=' {
            add(1, Act::Equals);
        }
        let op_of = |b: u8| match b {
            b';' => Some(Op::Semi),
            b'&' => Some(Op::And),
            b'|' => Some(Op::Or),
            _ => None,
        };
        if let Some(op) = op_of(first) {
            match st {
                Recurse if op != Op::Or => add(1, Act::PopThen(op)),
                CtabName if op != Op::Or => add(1, Act::Plain(op)),
                CtabVal if op != Op::Or => add(1, Act::CtabValOp(op)),
                PropVal => add(1, Act::PopThen(op)),
                _ => {}
            }
        }
        if first == b';' {
            add(1, Act::Plain(Op::Semi));
        }
        if first == b'&' {
            add(1, Act::Plain(Op::And));
        }
        if first == b'!' {
            add(1, Act::Plain(Op::Not));
        }
        if first == b'[' {
            add(1, Act::OpenBracket);
        }
        if first == b']' {
            match st {
                Recurse => add(1, Act::RecurseClose),
                PropVal => add(1, Act::PropValClose),
                SlnParam => add(1, Act::SlnParamClose),
                _ => {}
            }
        }
        if first == b'<' {
            add(1, Act::OpenAngle);
        }
        if first == b'>' {
            match st {
                CtabName | CtabVal => add(1, Act::PopAngle),
                PropVal => add(1, Act::PropValAngle),
                _ => {}
            }
        }
        if first == b'.' {
            add(1, Act::Separator);
        }
        if st == Recurse && first == b',' {
            add(1, Act::Comma);
        }
        if first == b'@' {
            add(1, Act::At);
        }
        if first == b'*' {
            add(1, Act::Asterix);
        }
        if first.is_ascii_digit() {
            add(
                s.iter().take_while(|b| b.is_ascii_digit()).count(),
                Act::Digit,
            );
        }
        if first == b'\n' {
            add(1, Act::Newline);
        }
        add(1, Act::Other); // `.`: any character but a newline
        if first == b'\n' {
            c.pop();
        }
        c
    }

    fn atom(&self, symbol: &str) -> Result<Atom> {
        let n = PeriodicTable::get().get_atomic_number(symbol)?;
        let mut a = if self.do_queries {
            query_atom_from_atom(&Atom::new(n))?
        } else {
            Atom::new(n)
        };
        // SLN has no concept of implicit Hs: they're either in the SLN or
        // they don't exist
        a.set_no_implicit(true);
        Ok(a)
    }

    fn query_atom(query: crate::graphmol::queries::AtomQuery) -> Tok {
        let mut a = Atom::new(0);
        a.set_query(Some(query));
        a.set_no_implicit(true);
        Tok::Atom(a)
    }

    /// The next token; `Eos` for ever at the end of the input.
    pub fn next_token(&mut self) -> Result<Tok> {
        if self.pos >= self.input.len() {
            return Ok(Tok::Eos);
        }
        let rest = &self.input[self.pos..];
        let mut best: Option<(usize, Act)> = None;
        for (len, act) in self.candidates(rest) {
            if best.is_none_or(|(l, _)| len > l) {
                best = Some((len, act));
            }
        }
        let (len, act) = best.expect("the catch-all rule matches");
        let text = String::from_utf8_lossy(&rest[..len]).into_owned();

        self.pos += len;
        use Act::*;
        Ok(match act {
            Text => Tok::Text(text),
            AtomSymbol => Tok::Atom(self.atom(&text)?),
            Any => {
                let mut a = Atom::new(0);
                a.set_query(Some(make_atom_null_query()));
                a.set_no_implicit(true);
                Tok::Atom(a)
            }
            Hev => {
                let mut q = make_atom_num_query(1);
                // FIX (upstream): are 2H or 3H heavy atoms or Hs?
                q.negate = true;
                Self::query_atom(q)
            }
            Hal => {
                let mut q = make_atom_num_query(9);
                for n in [17, 35, 53] {
                    q = expand_query(
                        q,
                        make_atom_num_query(n),
                        CompositeQueryType::Or,
                        true,
                        "Atom",
                    )?;
                }
                Self::query_atom(q)
            }
            Het => {
                let mut q = expand_query(
                    make_atom_num_query(6),
                    make_atom_num_query(1),
                    CompositeQueryType::Or,
                    true,
                    "Atom",
                )?;
                q.negate = true;
                Self::query_atom(q)
            }
            HBracket => {
                self.push(State::SlnParam);
                Tok::HBracket
            }
            HAsterix => Tok::HAsterix,
            H => Tok::H,
            Recurse => {
                self.push(State::Recurse);
                Tok::Recurse
            }
            NegRecurse => {
                self.push(State::Recurse);
                Tok::NegRecurse
            }
            Minus => Tok::Minus,
            Plus => Tok::Plus,
            Hash => Tok::Hash,
            Tilde => Tok::Tilde,
            Colon => Tok::Colon,
            OpenParen => Tok::OpenParen,
            CloseParen => Tok::CloseParen,
            Compare => {
                self.push(State::PropVal);
                Tok::Compare(text)
            }
            CtabEquals => {
                self.pop();
                self.push(State::CtabVal);
                Tok::Equals
            }
            CtabColonEquals => {
                self.pop();
                self.push(State::CtabVal);
                Tok::ColonEquals
            }
            CtabCaretEquals => {
                self.pop();
                self.push(State::CtabVal);
                Tok::CaretEquals
            }
            Equals => Tok::Equals,
            PopThen(op) => {
                self.pop();
                op.tok()
            }
            Plain(op) => op.tok(),
            CtabValOp(op) => {
                self.pop();
                self.push(State::CtabName);
                op.tok()
            }
            OpenBracket => {
                self.push(State::SlnParam);
                Tok::OpenBracket
            }
            RecurseClose => {
                // we're closing a recursive definition, which means we should
                // also be closing a parameter block
                self.pop();
                if self.state != State::SlnParam {
                    eprintln!(" after closing a recursion, we were not in the appropriate state.");
                } else {
                    self.pop();
                }
                Tok::CloseBracket
            }
            PropValClose => {
                // pop the property value state and, in a property block, the
                // block state
                self.pop();
                if self.state == State::SlnParam {
                    self.pop();
                }
                Tok::CloseBracket
            }
            SlnParamClose => {
                self.pop();
                Tok::CloseBracket
            }
            OpenAngle => {
                self.push(State::CtabName);
                Tok::OpenAngle
            }
            PopAngle => {
                self.pop();
                Tok::CloseAngle
            }
            PropValAngle => {
                self.pop();
                if self.state == State::CtabVal {
                    self.pop();
                }
                Tok::CloseAngle
            }
            Separator => Tok::Separator,
            Comma => Tok::Comma,
            At => Tok::At,
            Asterix => Tok::Asterix,
            Digit => Tok::Digit(text.parse::<i64>().map_or(i32::MAX, |v| {
                v.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
            })),
            Newline => Tok::End,
            Other => Tok::Char,
        })
    }
}
