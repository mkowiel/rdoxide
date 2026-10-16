//! Tautomer standardization of the structure checker.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/Tautomer.cpp` (2026.09.1).
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::options::StructCheckerOptions;
use crate::Result;
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};

/// Applies the tautomer rules of the options to a molecule.
#[derive(Debug)]
pub struct StructCheckTautomer<'a> {
    options: &'a StructCheckerOptions,
}

impl<'a> StructCheckTautomer<'a> {
    pub fn new(options: &'a StructCheckerOptions) -> Self {
        StructCheckTautomer { options }
    }

    /// `applyTautomer`: if the molecule matches the "from" pattern of rule
    /// `it` the bond orders, charges and radicals of the matched atoms are
    /// set to those of the "to" pattern. True if there was a match.
    pub fn apply_tautomer(&self, mol: &mut ROMol, it: usize) -> Result<bool> {
        let options = self.options;
        if options.from_tautomer.len() <= it || options.to_tautomer.len() <= it {
            if options.verbose {
                crate::rd_log_info!("ERROR: incorrect Tautomer index it={it}\n");
            }
            return Ok(false);
        }
        let from_tautomer = &options.from_tautomer[it];
        let to_tautomer = &options.to_tautomer[it];
        if to_tautomer.num_atoms() != from_tautomer.num_atoms() {
            if options.verbose {
                crate::rd_log_info!(
                    "ERROR: incorrect data toTautomer.getNumAtoms() != fromTautomer.getNumAtoms()\n"
                );
            }
            return Ok(false);
        }
        let params = SubstructMatchParameters {
            max_matches: 1,
            ..Default::default()
        };
        let matches = substruct_match(mol, from_tautomer, &params)?;
        let Some(match_) = matches.first() else {
            return Ok(false);
        };
        if options.verbose {
            crate::rd_log_info!(
                "found match for from_tautomer with {} atoms\n",
                to_tautomer.num_atoms()
            );
        }
        let invalid_idx = 1 + mol.num_atoms();
        let mut atom_idx_map = vec![invalid_idx; mol.num_atoms() as usize];
        for &(tai, mai) in match_ {
            atom_idx_map[mai as usize] = tai as u32;
        }
        // replace the bond order of the completely mapped bonds with that of
        // the bond of the "to" tautomer
        for b in 0..mol.num_bonds() {
            let (begin, end) = {
                let bond = mol.bond(b);
                (bond.begin_atom_idx(), bond.end_atom_idx())
            };
            let ti = atom_idx_map[begin as usize];
            let tj = atom_idx_map[end as usize];
            if ti == invalid_idx || tj == invalid_idx {
                continue;
            }
            if let Some(tb) = to_tautomer.get_bond_between_atoms(ti, tj)?
                && mol.bond(b).bond_type() != tb.bond_type()
            {
                let bt = tb.bond_type();
                mol.bond_mut(b).set_bond_type(bt);
            }
        }
        // apply the charge/radical fixes if any
        for &(tai, mai) in match_ {
            let ta = to_tautomer.atom(tai as u32);
            let (charge, radicals) = (ta.formal_charge(), ta.num_radical_electrons());
            let atom = mol.atom_mut(mai as u32);
            atom.set_formal_charge(charge);
            atom.set_num_radical_electrons(radicals);
        }
        Ok(true)
    }
}
