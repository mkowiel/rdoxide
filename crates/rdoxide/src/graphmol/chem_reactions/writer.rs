//! Ported from RDKit `Code/GraphMol/ChemReactions/ReactionWriter.cpp`:
//! reaction SMARTS/SMILES output (the rxn-file and single-molecule forms are in
//! [`super::mdl`]).

use super::{ChemicalReaction, ReactionMoleculeType, templates};
use crate::Result;
use crate::graphmol::mol_ops::frags::get_mol_frags_mapping;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};

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
