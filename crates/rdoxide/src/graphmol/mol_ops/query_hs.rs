//! Merging explicit hydrogens into query atoms (`MolOps::mergeQueryHs`,
//! `MolOps::hasQueryHs`): `C[H]` becomes `[C;!H0]`.
//!
//! Ported from RDKit `Code/GraphMol/AddHs.cpp` (2026.09.1).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::sync::Arc;

use crate::Result;
use crate::graphmol::atom::Atom;
use crate::graphmol::queries::{
    AtomQuery, CompositeQueryType, QueryKind, QueryMol, make_atom_h_count_query,
    make_atom_num_query,
};
use crate::graphmol::query_ops::atom_expand_query;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rd_log_warning;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::log::LogStateSetter;

/// `HydrogenType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HydrogenType {
    NotH,
    Query,
    Unmergeable,
}

/// `queryHasHs`: `(has an H, inside an OR)` for the direct children of a
/// compound query.
fn query_has_hs(query: &AtomQuery, in_or: bool) -> (bool, bool) {
    for child in &query.children {
        match child.description.as_str() {
            "AtomOr" => return query_has_hs(child, true),
            "AtomAtomicNum" if child.val == 1 && !child.negate => return (true, in_or),
            // 1001 is an aromatic hydrogen (not a thing, really), 1 an aliphatic one
            "AtomType" if (child.val == 1001 || child.val == 1) && !child.negate => {
                return (true, in_or);
            }
            _ => {}
        }
    }
    (false, in_or)
}

fn is_query_h(mol: &ROMol, idx: u32) -> HydrogenType {
    let atom = mol.atom(idx);
    let query = atom.query();
    if atom.atomic_num() == 1 && query.is_none_or(|q| !q.negate && q.description == "AtomAtomicNum")
    {
        // the atom is flagged as being an H and has no (other) query
        return HydrogenType::Query;
    }
    if mol.atom_degree(idx) > 1 {
        return HydrogenType::NotH;
    }
    let Some(query) = query.filter(|q| !q.negate) else {
        // negated queries are never merged
        return HydrogenType::NotH;
    };
    let (has_h, in_or) = match query.description.as_str() {
        "AtomOr" => query_has_hs(query, true),
        "AtomAnd" => query_has_hs(query, false),
        _ => (false, false),
    };
    match (has_h, in_or) {
        (true, true) => {
            rd_log_warning!(
                "WARNING: merging explicit H queries involved in ORs is not supported. This query will not be merged"
            );
            HydrogenType::Unmergeable
        }
        (true, false) => HydrogenType::Query,
        _ => HydrogenType::NotH,
    }
}

/// Applies [`merge_query_hs`] to the molecules of the recursive queries
/// anywhere in `query`.
fn merge_recursive_queries(
    query: &mut AtomQuery,
    merge_unmapped_only: bool,
    merge_isotopes: bool,
) -> Result<()> {
    if let QueryKind::Recursive { mol, .. } = &mut query.kind {
        let mut merged = (*mol.0).clone();
        merge_query_hs(&mut merged, merge_unmapped_only, merge_isotopes)?;
        *mol = QueryMol(Arc::new(merged));
    }
    for child in &mut query.children {
        merge_recursive_queries(child, merge_unmapped_only, merge_isotopes)?;
    }
    Ok(())
}

/// `MolOps::mergeQueryHs(mol, mergeUnmappedOnly, mergeIsotopes)`: removes the
/// explicit hydrogens bound to heavy atoms and adds `!H<n>` queries to those
/// atoms instead. Hydrogens that are not bound to a heavy atom stay, and with
/// `merge_unmapped_only` so do atom-mapped ones.
pub fn merge_query_hs(
    mol: &mut RWMol,
    merge_unmapped_only: bool,
    merge_isotopes: bool,
) -> Result<()> {
    let is_h: Vec<bool> = (0..mol.num_atoms())
        .map(|i| is_query_h(mol, i) == HydrogenType::Query)
        .collect();
    let mut atoms_to_remove = Vec::new();
    for idx in 0..mol.num_atoms() {
        if is_h[idx as usize] {
            continue;
        }
        let mut num_hs_to_remove = 0;
        for nbr in mol.atom_neighbors(idx).collect::<Vec<_>>() {
            if !is_h[nbr as usize] {
                continue;
            }
            let h = mol.atom(nbr);
            let unmapped_ok = !merge_unmapped_only || !h.has_prop(cp::molAtomMapNumber);
            let isotope_ok = merge_isotopes || h.isotope() == 0;
            if unmapped_ok && isotope_ok {
                atoms_to_remove.push(nbr);
                num_hs_to_remove += 1;
            }
        }
        if num_hs_to_remove > 0 {
            // C[H] -> [C;!H0], C([H])[H] -> [C;!H0;!H1]
            if !mol.atom(idx).has_query() {
                let mut replacement = Atom::new(0);
                replacement.set_query(Some(make_atom_num_query(mol.atom(idx).atomic_num() as i32)));
                replacement.props.update(&mol.atom(idx).props, false);
                mol.replace_atom(idx, &replacement, false)?;
            }
            for count in 0..num_hs_to_remove {
                let mut no_hs = make_atom_h_count_query(count);
                no_hs.negate = true;
                atom_expand_query(mol.atom_mut(idx), no_hs, CompositeQueryType::And, true)?;
            }
        }
        // recurse into recursive queries (github issue 544)
        if let Some(mut query) = mol.atom_mut(idx).query.take() {
            let merged = merge_recursive_queries(&mut query, merge_unmapped_only, merge_isotopes);
            mol.atom_mut(idx).query = Some(query);
            merged?;
        }
    }
    mol.begin_batch_edit()?;
    for idx in atoms_to_remove {
        mol.remove_atom(idx)?;
    }
    mol.commit_batch_edit()
}

/// `MolOps::hasQueryHs(mol)`: `(has query Hs, has Hs that mergeQueryHs cannot
/// remove)`.
pub fn has_query_hs(mol: &ROMol) -> Result<(bool, bool)> {
    // ORs and the like are not announced here
    let _blocker = LogStateSetter::new();
    let mut query_hs = false;
    for idx in 0..mol.num_atoms() {
        match is_query_h(mol, idx) {
            HydrogenType::Unmergeable => return Ok((true, true)),
            HydrogenType::Query => query_hs = true,
            HydrogenType::NotH => {}
        }
        if let Some(query) = mol.atom(idx).query() {
            let (found, unmergeable) = recursive_query_hs(query)?;
            if unmergeable {
                return Ok((true, true));
            }
            query_hs |= found;
        }
    }
    Ok((query_hs, false))
}

fn recursive_query_hs(query: &AtomQuery) -> Result<(bool, bool)> {
    let mut found = false;
    if let QueryKind::Recursive { mol, .. } = &query.kind {
        let (has, unmergeable) = has_query_hs(&mol.0)?;
        if unmergeable {
            return Ok((true, true));
        }
        found |= has;
    }
    for child in &query.children {
        let (has, unmergeable) = recursive_query_hs(child)?;
        if unmergeable {
            return Ok((true, true));
        }
        found |= has;
    }
    Ok((found, false))
}
