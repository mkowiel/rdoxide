//! Sybyl line notation (SLN) parsing.
//!
//! Ported from RDKit `Code/GraphMol/SLNParse` (2026.09.1). The flex/bison
//! sources (`sln.ll`, `sln.yy`) are replaced by a hand-written scanner
//! ([`lexer`]) and a recursive descent parser ([`parser`]) with the same
//! semantic actions.
//!
//! Macro atoms (`{Z:O}`) are replaced as plain text; upstream uses the macro
//! name as a regular expression.
//!
//! Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

pub mod attribs;
mod lexer;
mod ops;
mod parser;

use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::queries::CompositeQueryType;
use crate::graphmol::queries::make_atom_h_count_query;
use crate::graphmol::query_ops::atom_expand_query;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result};

/// What can go wrong while parsing.
#[derive(Debug)]
pub enum SlnError {
    /// An `SLNParseException`
    Sln(String),
    /// Any other error
    Other(RDKitError),
    /// A syntax error
    Syntax,
    /// A semantic error that aborts the parse (after being reported)
    Abort(String),
    /// `YYACCEPT`
    Accept,
}

impl From<RDKitError> for SlnError {
    fn from(e: RDKitError) -> Self {
        SlnError::Other(e)
    }
}

/// `sln_parse`: the molecules built by the parse.
fn sln_parse(inp: &str, do_queries: bool) -> std::result::Result<Option<ops::SlnMol>, SlnError> {
    let mut p = parser::Parser::new(inp, do_queries);
    let res = match p.cmpd() {
        Err(SlnError::Accept) => Ok(()),
        Ok(_) | Err(SlnError::Syntax) => Err(SlnError::Syntax),
        Err(e) => Err(e),
    };
    match res {
        Ok(()) => Ok(p.list.into_iter().next().flatten()),
        Err(SlnError::Syntax) => {
            crate::rd_log_error!("SLN Parse Error: syntax error while parsing: {inp}\n");
            Err(SlnError::Sln(format!("Failed parsing SLN '{inp}'")))
        }
        Err(SlnError::Abort(msg)) => {
            crate::rd_log_error!("SLN Parse Error: {msg} while parsing: {inp}\n");
            Err(SlnError::Sln(format!("Failed parsing SLN '{inp}'")))
        }
        Err(e) => Err(e),
    }
}

/// `SLNParse::finalizeQueryMol`
pub fn finalize_query_mol(mol: &mut ROMol, merge_hs: bool) -> std::result::Result<(), SlnError> {
    // do we need to remove the Hs from the molecule?
    if merge_hs {
        for idx in 0..mol.num_atoms() {
            // set a query for the H count
            let n = mol.atom(idx).num_explicit_hs();
            if n != 0 {
                atom_expand_query(
                    mol.atom_mut(idx),
                    make_atom_h_count_query(n as i32),
                    CompositeQueryType::And,
                    true,
                )?;
            }
        }
    }
    // we don't want to sanitize, but we do need to get some ring info
    symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    let mut root_idx: i32 = -1;
    for idx in 0..mol.num_atoms() {
        attribs::parse_final_atom_attribs(mol, idx, true)?;
        if mol.atom(idx).has_prop(cp::_starred) {
            if root_idx > -1 {
                crate::rd_log_error!(
                    "SLN Error: multiple starred atoms in a recursive query. Extra stars ignored\n"
                );
            } else {
                root_idx = idx as i32;
            }
        }
    }
    if root_idx > -1 {
        mol.set_prop(cp::_queryRootAtom, root_idx, false)?;
    }
    Ok(())
}

/// `replaceSLNMacroAtoms`: removes the macro definitions `{name:value}` and
/// replaces the macro names by their values.
fn replace_sln_macro_atoms(inp: &str) -> String {
    // the first `\{(.+?):(.+?)\}` at or after `from` (`.` doesn't match line
    // terminators)
    fn find_macro(s: &[char], from: usize) -> Option<(usize, usize, usize, usize)> {
        let not_nl = |c: char| !matches!(c, '\n' | '\r' | '\u{2028}' | '\u{2029}');
        for start in from..s.len() {
            if s[start] != '{' {
                continue;
            }
            // (start, colon, close)
            let mut colon = start + 2;
            while colon < s.len() && not_nl(s[colon - 1]) {
                if s[colon] == ':' {
                    let mut close = colon + 2;
                    while close < s.len() && not_nl(s[close - 1]) {
                        if s[close] == '}' {
                            return Some((start, colon, close, close + 1));
                        }
                        close += 1;
                    }
                }
                colon += 1;
            }
        }
        None
    }
    let chars: Vec<char> = inp.chars().collect();
    let mut macros: Vec<(String, String)> = Vec::new();
    let mut res = String::new();
    let mut pos = 0;
    while let Some((start, colon, close, end)) = find_macro(&chars, pos) {
        res.extend(&chars[pos..start]);
        macros.push((
            chars[start + 1..colon].iter().collect(),
            chars[colon + 1..close].iter().collect(),
        ));
        pos = end;
    }
    if macros.is_empty() {
        return inp.to_string();
    }
    res.extend(&chars[pos..]);
    for (name, value) in macros {
        res = res.replace(&name, &value);
    }
    res
}

/// `SLNParse::toMol`
fn to_mol(inp: &str, do_queries: bool, debug_parse: i32) -> Result<Option<RWMol>> {
    let inp = inp.trim_matches([' ', '\t', '\r', '\n']);
    let inp = replace_sln_macro_atoms(inp);
    if debug_parse != 0 {
        eprintln!("****** PARSING SLN: ->{inp}<-");
    }
    let first = match sln_parse(&inp, do_queries) {
        Ok(m) => m,
        Err(SlnError::Other(e)) => return Err(e),
        Err(SlnError::Sln(msg)) => {
            crate::rd_log_error!("{msg}\n");
            return Ok(None);
        }
        Err(_) => unreachable!("reported by sln_parse"),
    };
    let Some(first) = first else {
        return Ok(None);
    };
    // ring closures without an opener
    for key in first.open_bond_bookmarks() {
        if key > 0 && key < first.mol.num_atoms() as i32 {
            crate::rd_log_error!(
                "SLN Parser error: Ring closure {key} does not have a corresponding opener.\n"
            );
            return Ok(None);
        }
    }
    let mut res = first.mol;
    // set up the chirality flags as soon as the molecule is finished since
    // the Hs will be removed later and that will break things
    match attribs::adjust_atom_chiralities(&mut res) {
        Ok(()) => Ok(Some(res)),
        Err(e) => Err(into_error(e)),
    }
}

/// `SLNToMol`
pub fn sln_to_mol(sln: &str, sanitize: bool, debug_parse: i32) -> Result<Option<RWMol>> {
    let Some(mut res) = to_mol(sln, false, debug_parse)? else {
        return Ok(None);
    };
    for idx in 0..res.num_atoms() {
        attribs::parse_final_atom_attribs(&mut res, idx, false).map_err(into_error)?;
    }
    if sanitize {
        // removing the explicit Hs from the graph triggers a sanitization
        remove_hs(&mut res, &RemoveHsParameters::default(), true)?;
    }
    Ok(Some(res))
}

/// `SLNQueryToMol`
pub fn sln_query_to_mol(sln: &str, merge_hs: bool, debug_parse: i32) -> Result<Option<RWMol>> {
    let Some(mut res) = to_mol(sln, true, debug_parse)? else {
        return Ok(None);
    };
    finalize_query_mol(&mut res, merge_hs).map_err(into_error)?;
    Ok(Some(res))
}

fn into_error(e: SlnError) -> RDKitError {
    match e {
        SlnError::Other(e) => e,
        SlnError::Sln(m) | SlnError::Abort(m) => RDKitError::Value(m),
        SlnError::Syntax => RDKitError::Value("syntax error".into()),
        SlnError::Accept => RDKitError::Value("unexpected accept".into()),
    }
}
