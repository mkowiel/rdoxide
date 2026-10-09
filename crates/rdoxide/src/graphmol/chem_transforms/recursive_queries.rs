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
