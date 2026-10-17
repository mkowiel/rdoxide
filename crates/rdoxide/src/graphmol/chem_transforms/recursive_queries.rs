//! Ported from RDKit `Code/GraphMol/ChemTransforms/ChemTransforms.cpp`
//! (`addRecursiveQueries`): attaches recursive-SMARTS queries, looked up by
//! name in an atom property, to the atoms of a molecule.

use std::collections::BTreeMap;

use crate::graphmol::queries::{
    AtomQuery, CompositeQueryType, Query, QueryKind, make_recursive_structure_query,
};
use crate::graphmol::query_ops::{atom_expand_query, query_atom_from_atom};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::rdgeneral::RDProps;
use crate::{RDKitError, Result};

const DELIMITER: char = ',';

/// The query for the (already lower-cased) property value `pval`, or `None`
/// if a name is missing from `queries`.
fn query_from_names(pval: &str, queries: &BTreeMap<String, ROMol>) -> Option<AtomQuery> {
    let recursive = |name: &str| {
        queries
            .get(name)
            .map(|mol| make_recursive_structure_query(mol.clone(), 0))
    };
    if !pval.contains(DELIMITER) {
        return recursive(pval);
    }
    let mut or = Query::new(QueryKind::Or);
    or.description = "AtomOr".into();
    for name in pval.split(DELIMITER).filter(|t| !t.is_empty()) {
        or.children.push(recursive(name)?);
    }
    Some(or)
}

/// `addRecursiveQueries`: for every atom with property `prop_name`, ANDs a
/// recursive query to the atom's query. The property is a (comma separated
/// list of) name(s) in `queries`; anything else must be a SMARTS pattern.
/// Returns the `(atom index, lower-cased value)` pairs of the atoms found.
pub fn add_recursive_queries(
    mol: &mut ROMol,
    queries: &BTreeMap<String, ROMol>,
    prop_name: &str,
) -> Result<Vec<(u32, String)>> {
    let mut labels = Vec::new();
    for idx in 0..mol.num_atoms() {
        let Some(value) = mol.atom(idx).get_prop(prop_name) else {
            continue;
        };
        let original = value.to_rd_string();
        let pval = original.to_lowercase();
        labels.push((idx, pval.clone()));

        let query = match query_from_names(&pval, queries) {
            Some(q) => q,
            None => {
                // maybe the value is a SMARTS expression
                let smarts = mol_from_smarts(&original, &SmartsParserParams::default())
                    .ok()
                    .flatten()
                    .ok_or_else(|| RDKitError::Key(pval.clone()))?;
                make_recursive_structure_query(smarts, 0)
            }
        };
        if !mol.atom(idx).has_query() {
            let query_atom = query_atom_from_atom(mol.atom(idx))?;
            mol.replace_atom(idx, &query_atom, false)?;
        }
        atom_expand_query(mol.atom_mut(idx), query, CompositeQueryType::And, true)?;
    }
    Ok(labels)
}

/// The Python wrapper's `addRecursiveQuery`: ANDs (or, without
/// `preserve_existing_query`, replaces the query with) a recursive query for a
/// copy of `query` on atom `atom_idx`.
pub fn add_recursive_query(
    mol: &mut ROMol,
    query: &ROMol,
    atom_idx: u32,
    preserve_existing_query: bool,
) -> Result<()> {
    if atom_idx >= mol.num_atoms() {
        return Err(RDKitError::Value(
            "atom index exceeds mol.GetNumAtoms()".into(),
        ));
    }
    let q = make_recursive_structure_query(query.clone(), 0);
    if !mol.atom(atom_idx).has_query() {
        let query_atom = query_atom_from_atom(mol.atom(atom_idx))?;
        mol.replace_atom(atom_idx, &query_atom, false)?;
    }
    if preserve_existing_query {
        atom_expand_query(mol.atom_mut(atom_idx), q, CompositeQueryType::And, true)
    } else {
        mol.atom_mut(atom_idx).set_query(Some(q));
        Ok(())
    }
}

/// `parseQueryDefText`: the query definitions of a text with one
/// `name<delimiter>SMARTS` per line. Lines starting with `comment` are
/// skipped, and so are the ones whose SMARTS cannot be parsed (with a
/// warning). Names are lower-cased with `standardize`.
pub fn parse_query_def_text(
    text: &str,
    standardize: bool,
    delimiter: &str,
    comment: &str,
    name_column: usize,
    smarts_column: usize,
) -> BTreeMap<String, ROMol> {
    let mut defs = BTreeMap::new();
    let trim = |s: &str| s.trim_matches([' ', '\t']).to_string();
    for (i, raw) in text.lines().enumerate() {
        let line = raw.trim_end_matches('\r');
        if line.is_empty() || line.starts_with(comment) {
            continue;
        }
        // boost::char_separator: every character of the delimiter separates
        // and empty tokens are dropped
        let tokens: Vec<&str> = line
            .split(|c| delimiter.contains(c))
            .filter(|t| !t.is_empty())
            .collect();
        let name = tokens.get(name_column).map(|t| trim(t)).unwrap_or_default();
        let smarts = tokens
            .get(smarts_column)
            .map(|t| trim(t))
            .unwrap_or_default();
        if name.is_empty() || smarts.is_empty() {
            continue;
        }
        match mol_from_smarts(&smarts, &SmartsParserParams::default())
            .ok()
            .flatten()
        {
            Some(m) => {
                let key = if standardize {
                    name.to_lowercase()
                } else {
                    name
                };
                defs.insert(key, m);
            }
            None => crate::rd_log_warning!(
                "cannot convert SMARTS {smarts} to molecule at line {}",
                i + 1
            ),
        }
    }
    defs
}

/// `parseQueryDefFile`
pub fn parse_query_def_file(
    filename: &str,
    standardize: bool,
    delimiter: &str,
    comment: &str,
    name_column: usize,
    smarts_column: usize,
) -> Result<BTreeMap<String, ROMol>> {
    let text = std::fs::read(filename)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {filename}")))?;
    Ok(parse_query_def_text(
        &String::from_utf8_lossy(&text),
        standardize,
        delimiter,
        comment,
        name_column,
        smarts_column,
    ))
}
