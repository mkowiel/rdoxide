//! A query that matches all the tautomers of a molecule (upstream
//! `Code/GraphMol/TautomerQuery`).
//!
//! The query is the molecule with the atoms and bonds that the tautomer
//! transforms modify relaxed to generic queries; a match of this template is
//! accepted when one of the tautomers is compatible at the modified atoms and
//! bonds. Not ported: Boost serialization (`toStream`/`serialize`), as in an
//! upstream build without it.

use std::collections::BTreeSet;
use std::sync::{Arc, Mutex};

use super::atom::Atom;
use super::bond::{Bond, BondType};
use super::fingerprints::pattern::pattern_fingerprint_mol;
use super::mol_standardize::catalogs::TautomerCatalogParams;
use super::mol_standardize::tautomer::TautomerEnumerator;
use super::queries::{make_atom_num_query, make_single_or_double_or_aromatic_bond_query};
use super::romol::{ROMol, RWMol};
use super::substruct::{
    FinalMatchCheck, MatchVect, SubstructMatchParameters, atom_compat, bond_compat, substruct_match,
};
use crate::Result;
use crate::datastructs::ExplicitBitVect;
use crate::rdgeneral::RDProps;

/// `TautomerQuery`
#[derive(Debug, Clone)]
pub struct TautomerQuery {
    tautomers: Vec<Arc<ROMol>>,
    template_molecule: Arc<ROMol>,
    modified_atoms: Vec<u32>,
    modified_bonds: Vec<u32>,
}

/// Whether `tautomer` agrees with `mol` at the modified atoms and bonds of
/// the match (`TautomerQuery::matchTautomer`).
fn match_tautomer(
    mol: &ROMol,
    tautomer: &ROMol,
    matched: &[u32],
    modified_atoms: &[u32],
    modified_bonds: &[u32],
    params: &SubstructMatchParameters,
) -> bool {
    for &idx in modified_atoms {
        if !atom_compat(tautomer, idx, mol, matched[idx as usize], params).unwrap_or(false) {
            return false;
        }
    }
    for &idx in modified_bonds {
        let query_bond = tautomer.bond(idx);
        let begin = matched[query_bond.begin_atom_idx() as usize];
        let end = matched[query_bond.end_atom_idx() as usize];
        let Some(target_bond) = mol.bond_idx_between(begin, end) else {
            return false;
        };
        if !bond_compat(tautomer, idx, mol, target_bond, params).unwrap_or(false) {
            return false;
        }
    }
    true
}

impl TautomerQuery {
    pub fn new(
        tautomers: Vec<Arc<ROMol>>,
        template_molecule: ROMol,
        modified_atoms: Vec<u32>,
        modified_bonds: Vec<u32>,
    ) -> Self {
        TautomerQuery {
            tautomers,
            template_molecule: Arc::new(template_molecule),
            modified_atoms,
            modified_bonds,
        }
    }

    /// `fromMol`: the tautomer query of `query`, using the default tautomer
    /// transforms or those of `tautomer_transform_file`.
    pub fn from_mol(query: &ROMol, tautomer_transform_file: &str) -> Result<Self> {
        let enumerator = if tautomer_transform_file.is_empty() {
            TautomerEnumerator::default()
        } else {
            TautomerEnumerator::with_params(TautomerCatalogParams::from_file(
                tautomer_transform_file,
            )?)
        };
        let res = enumerator.enumerate(query)?;
        let modified_atoms = res.modified_atom_indices();
        let modified_bonds = res.modified_bond_indices();

        let mut template: RWMol = query.clone();
        for &idx in &modified_atoms {
            let atom = query.atom(idx);
            let mut query_atom = Atom::new(atom.atomic_num());
            query_atom.set_query(Some(match atom.query() {
                Some(q) => q.clone(),
                None => make_atom_num_query(atom.atomic_num() as i32),
            }));
            template.replace_atom(idx, &query_atom, true)?;
        }
        for &idx in &modified_bonds {
            let mut query_bond = Bond::new(BondType::Unspecified);
            query_bond.set_query(Some(make_single_or_double_or_aromatic_bond_query()));
            template.replace_bond(idx, &query_bond, true, false)?;
        }

        let tautomers = res
            .tautomer_mols()
            .into_iter()
            .map(|t| {
                let mut t = (*t).clone();
                t.update_props(query.dict(), false);
                Arc::new(t)
            })
            .collect();
        Ok(TautomerQuery::new(
            tautomers,
            template,
            modified_atoms,
            modified_bonds,
        ))
    }

    /// `getTemplateMolecule`
    pub fn template_molecule(&self) -> &ROMol {
        &self.template_molecule
    }

    /// `getTautomers`
    pub fn tautomers(&self) -> &[Arc<ROMol>] {
        &self.tautomers
    }

    /// `getModifiedAtoms`
    pub fn modified_atoms(&self) -> &[u32] {
        &self.modified_atoms
    }

    /// `getModifiedBonds`
    pub fn modified_bonds(&self) -> &[u32] {
        &self.modified_bonds
    }

    /// `substructOf`: the matches of the query in `mol` and, for each, the
    /// tautomer that matched.
    pub fn substruct_of(
        &self,
        mol: &ROMol,
        params: &SubstructMatchParameters,
    ) -> Result<(Vec<MatchVect>, Vec<Arc<ROMol>>)> {
        let found: Arc<Mutex<Vec<Arc<ROMol>>>> = Arc::default();
        let mut template_params = params.clone();
        template_params.uniquify = false;
        {
            let tautomers = self.tautomers.clone();
            let modified_atoms = self.modified_atoms.clone();
            let modified_bonds = self.modified_bonds.clone();
            let user_check = params.extra_final_check.clone();
            let inner_params = params.clone();
            let found = found.clone();
            template_params.extra_final_check = Some(FinalMatchCheck::new(move |mol, matched| {
                for tautomer in &tautomers {
                    if match_tautomer(
                        mol,
                        tautomer,
                        matched,
                        &modified_atoms,
                        &modified_bonds,
                        &inner_params,
                    ) {
                        let accepted = user_check.as_ref().is_none_or(|c| c.call(mol, matched));
                        if accepted {
                            found.lock().expect("tautomer lock").push(tautomer.clone());
                        }
                        return accepted;
                    }
                }
                false
            }));
        }
        let mut matches = substruct_match(mol, &self.template_molecule, &template_params)?;
        let mut matching = std::mem::take(&mut *found.lock().expect("tautomer lock"));
        if params.uniquify && matches.len() > 1 {
            remove_tautomer_duplicates(&mut matches, &mut matching, mol.num_atoms());
        }
        Ok((matches, matching))
    }

    /// `isSubstructOf`
    pub fn is_substruct_of(&self, mol: &ROMol, params: &SubstructMatchParameters) -> Result<bool> {
        let mut params = params.clone();
        params.max_matches = 1;
        Ok(!self.substruct_of(mol, &params)?.0.is_empty())
    }

    /// `patternFingerprintTemplate`
    pub fn pattern_fingerprint_template(&self, fp_size: u32) -> Result<ExplicitBitVect> {
        pattern_fingerprint_mol(&self.template_molecule, fp_size, true)
    }

    /// `patternFingerprintTarget`
    pub fn pattern_fingerprint_target(target: &ROMol, fp_size: u32) -> Result<ExplicitBitVect> {
        pattern_fingerprint_mol(target, fp_size, true)
    }
}

/// `removeTautomerDuplicates`: keeps the first match for each set of atoms,
/// with its tautomer.
fn remove_tautomer_duplicates(
    matches: &mut Vec<MatchVect>,
    tautomers: &mut Vec<Arc<ROMol>>,
    n_atoms: u32,
) {
    let mut seen: BTreeSet<Vec<bool>> = BTreeSet::new();
    let mut kept_matches = Vec::with_capacity(matches.len());
    let mut kept_tautomers = Vec::with_capacity(tautomers.len());
    for (m, t) in std::mem::take(matches)
        .into_iter()
        .zip(std::mem::take(tautomers))
    {
        let mut val = vec![false; n_atoms as usize];
        for &(_, mol_idx) in &m {
            val[mol_idx as usize] = true;
        }
        if seen.insert(val) {
            kept_matches.push(m);
            kept_tautomers.push(t);
        }
    }
    *matches = kept_matches;
    *tautomers = kept_tautomers;
}

/// `SubstructMatch(mol, TautomerQuery, params)`
pub fn tautomer_substruct_match(
    mol: &ROMol,
    query: &TautomerQuery,
    params: &SubstructMatchParameters,
) -> Result<Vec<MatchVect>> {
    Ok(query.substruct_of(mol, params)?.0)
}
