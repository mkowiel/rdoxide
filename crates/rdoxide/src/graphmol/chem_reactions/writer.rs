//! Ported from RDKit `Code/GraphMol/ChemReactions/ReactionWriter.cpp`:
//! reaction SMARTS/SMILES output and the single-molecule reaction encoding.

use super::{ChemicalReaction, ReactionMoleculeType, templates};
use crate::Result;
use crate::graphmol::mol_ops::frags::get_mol_frags_mapping;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

fn mol_to_string(mol: &ROMol, to_smiles: bool, params: &SmilesWriteParams) -> Result<String> {
    let mut copy = mol.clone();
    let mut res = if to_smiles {
        mol_to_smiles(&mut copy, params)?
    } else {
        mol_to_smarts(&mut copy, params)?
    };
    if get_mol_frags_mapping(mol).0 > 1 {
        res = format!("({res})");
    }
    Ok(res)
}

fn templates_to_string(
    rxn: &ChemicalReaction,
    t: ReactionMoleculeType,
    to_smiles: bool,
    params: &SmilesWriteParams,
) -> Result<String> {
    let mut parts = templates(rxn, t)
        .iter()
        .map(|m| mol_to_string(m, to_smiles, params))
        .collect::<Result<Vec<_>>>()?;
    if params.canonical {
        parts.sort();
    }
    Ok(parts.join("."))
}

fn reaction_to_string(
    rxn: &ChemicalReaction,
    to_smiles: bool,
    params: &SmilesWriteParams,
) -> Result<String> {
    Ok(format!(
        "{}>{}>{}",
        templates_to_string(rxn, ReactionMoleculeType::Reactant, to_smiles, params)?,
        templates_to_string(rxn, ReactionMoleculeType::Agent, to_smiles, params)?,
        templates_to_string(rxn, ReactionMoleculeType::Product, to_smiles, params)?,
    ))
}

/// `ChemicalReactionToRxnSmarts`
pub fn reaction_to_smarts(rxn: &ChemicalReaction, params: &SmilesWriteParams) -> Result<String> {
    reaction_to_string(rxn, false, params)
}

/// `ChemicalReactionToRxnSmiles`
pub fn reaction_to_smiles(rxn: &ChemicalReaction, params: &SmilesWriteParams) -> Result<String> {
    reaction_to_string(rxn, true, params)
}

/// `ChemicalReactionToRxnMol`: one molecule holding every template, with each
/// atom's `molRxnRole` set to 1 (reactant), 2 (product) or 3 (agent). Unlike
/// upstream, the reaction's own templates are left untouched.
pub fn reaction_to_rxn_mol(rxn: &ChemicalReaction) -> Result<ROMol> {
    let mut res = RWMol::default();
    for (t, role) in [
        (ReactionMoleculeType::Reactant, 1i32),
        (ReactionMoleculeType::Product, 2),
        (ReactionMoleculeType::Agent, 3),
    ] {
        for tpl in templates(rxn, t) {
            let mut copy = tpl.clone();
            for i in 0..copy.num_atoms() {
                copy.atom_mut(i).dict_mut().set_val(cp::molRxnRole, role);
            }
            res.insert_mol(&copy)?;
        }
    }
    Ok(res)
}
