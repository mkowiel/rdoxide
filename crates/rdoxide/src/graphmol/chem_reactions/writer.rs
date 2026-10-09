//! Ported from RDKit `Code/GraphMol/ChemReactions/ReactionWriter.cpp`:
//! reaction SMARTS/SMILES output (the rxn-file and single-molecule forms are in
//! [`super::mdl`]).

use super::{ChemicalReaction, ReactionMoleculeType, templates};
use crate::Result;
use crate::graphmol::mol_ops::frags::get_mol_frags_mapping;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::cx_write::get_cx_extensions_for_mols;
use crate::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};

/// Writes `mol` (a copy, which keeps the SMILES output order properties the
/// CX extensions need).
fn mol_to_string(
    mol: &ROMol,
    to_smiles: bool,
    params: &SmilesWriteParams,
) -> Result<(String, ROMol)> {
    let mut copy = mol.clone();
    let mut res = if to_smiles {
        mol_to_smiles(&mut copy, params)?
    } else {
        mol_to_smarts(&mut copy, params)?
    };
    if get_mol_frags_mapping(mol).0 > 1 {
        res = format!("({res})");
    }
    Ok((res, copy))
}

/// The written templates of one role, with the (unsorted) written molecules
/// appended to `written`.
fn templates_to_string(
    rxn: &ChemicalReaction,
    t: ReactionMoleculeType,
    to_smiles: bool,
    params: &SmilesWriteParams,
    written: &mut Vec<ROMol>,
) -> Result<String> {
    let mut parts = Vec::new();
    for mol in templates(rxn, t) {
        let (text, copy) = mol_to_string(mol, to_smiles, params)?;
        parts.push(text);
        written.push(copy);
    }
    if params.canonical {
        parts.sort();
    }
    Ok(parts.join("."))
}

fn reaction_to_string(
    rxn: &ChemicalReaction,
    to_smiles: bool,
    params: &SmilesWriteParams,
    cx_flags: Option<u32>,
) -> Result<String> {
    let mut written = Vec::new();
    let mut res = String::new();
    for (i, t) in [
        ReactionMoleculeType::Reactant,
        ReactionMoleculeType::Agent,
        ReactionMoleculeType::Product,
    ]
    .into_iter()
    .enumerate()
    {
        if i > 0 {
            res.push('>');
        }
        res += &templates_to_string(rxn, t, to_smiles, params, &mut written)?;
    }
    if let Some(flags) = cx_flags {
        let mols: Vec<&ROMol> = written.iter().collect();
        let ext = get_cx_extensions_for_mols(&mols, flags)?;
        if !ext.is_empty() {
            res.push(' ');
            res += &ext;
        }
    }
    Ok(res)
}

/// `ChemicalReactionToRxnSmarts`
pub fn reaction_to_smarts(rxn: &ChemicalReaction, params: &SmilesWriteParams) -> Result<String> {
    reaction_to_string(rxn, false, params, None)
}

/// `ChemicalReactionToRxnSmiles`
pub fn reaction_to_smiles(rxn: &ChemicalReaction, params: &SmilesWriteParams) -> Result<String> {
    reaction_to_string(rxn, true, params, None)
}

/// `ChemicalReactionToCXRxnSmarts`: reaction SMARTS followed by the CX
/// extensions selected by `flags` (see [`cx_fields`]).
pub fn reaction_to_cx_smarts(
    rxn: &ChemicalReaction,
    params: &SmilesWriteParams,
    flags: u32,
) -> Result<String> {
    reaction_to_string(rxn, false, params, Some(flags))
}

/// `ChemicalReactionToCXRxnSmiles`
pub fn reaction_to_cx_smiles(
    rxn: &ChemicalReaction,
    params: &SmilesWriteParams,
    flags: u32,
) -> Result<String> {
    reaction_to_string(rxn, true, params, Some(flags))
}
