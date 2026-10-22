//! Reactions in PNG metadata.
//!
//! Ported from RDKit `Code/GraphMol/ChemReactions/PNGParser.cpp` (2026.09.1).
//!
//! Copyright (C) 2020 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::ChemicalReaction;
use super::daylight::{ReactionSmartsParserParams, reaction_from_smarts, reaction_from_smiles};
use super::mdl::{reaction_from_rxn_block, reaction_to_rxn_block};
use super::pickler::{pickle_reaction_default, reaction_from_pickle};
use super::writer::{reaction_to_smarts, reaction_to_smiles};
use crate::graphmol::file_parsers::MolFileParserParams;
use crate::graphmol::file_parsers::png::{
    PngMetadata, add_metadata_to_png, augment_tag_name, png_to_metadata,
};
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::{RDKitError, Result};

/// `PNGData::rxnSmilesTag`
pub const RXN_SMILES_TAG: &str = "ReactionSmiles";
/// `PNGData::rxnSmartsTag`
pub const RXN_SMARTS_TAG: &str = "ReactionSmarts";
/// `PNGData::rxnRxnTag`
pub const RXN_RXN_TAG: &str = "ReactionRxn";
/// `PNGData::rxnPklTag`
pub const RXN_PKL_TAG: &str = "rdkitReactionPKL";

/// What `add_chemical_reaction_to_png` stores.
#[derive(Debug, Clone, Copy)]
pub struct ReactionPngOptions {
    /// the reaction pickle
    pub include_pkl: bool,
    /// reaction SMILES
    pub include_smiles: bool,
    /// reaction SMARTS
    pub include_smarts: bool,
    /// an RXN block
    pub include_rxn: bool,
}

impl Default for ReactionPngOptions {
    fn default() -> Self {
        ReactionPngOptions {
            include_pkl: true,
            include_smiles: true,
            include_smarts: false,
            include_rxn: false,
        }
    }
}

/// `addChemicalReactionToPNGStream`: the PNG with the reaction added as metadata.
pub fn add_chemical_reaction_to_png(
    rxn: &ChemicalReaction,
    data: &[u8],
    opts: &ReactionPngOptions,
) -> Result<Vec<u8>> {
    let mut metadata = PngMetadata::new();
    if opts.include_pkl {
        metadata.push((augment_tag_name(RXN_PKL_TAG), pickle_reaction_default(rxn)?));
    }
    if opts.include_smiles {
        let smi = reaction_to_smiles(rxn, &SmilesWriteParams::default())?;
        metadata.push((augment_tag_name(RXN_SMILES_TAG), smi.into_bytes()));
    }
    if opts.include_smarts {
        let sma = reaction_to_smarts(rxn, &SmilesWriteParams::default())?;
        metadata.push((augment_tag_name(RXN_SMARTS_TAG), sma.into_bytes()));
    }
    if opts.include_rxn {
        let mb = reaction_to_rxn_block(rxn, false, false)?;
        metadata.push((augment_tag_name(RXN_RXN_TAG), mb.into_bytes()));
    }
    add_metadata_to_png(data, &metadata, true)
}

/// `v2::ReactionParser::ReactionFromPNGStream`: the reaction of the first
/// reaction tag of the PNG.
pub fn png_to_chemical_reaction(data: &[u8]) -> Result<ChemicalReaction> {
    for (key, value) in png_to_metadata(data)? {
        let text = || String::from_utf8_lossy(&value).into_owned();
        if key.starts_with(RXN_PKL_TAG) {
            return reaction_from_pickle(&value);
        } else if key.starts_with(RXN_SMILES_TAG) {
            return reaction_from_smiles(&text(), &ReactionSmartsParserParams::new());
        } else if key.starts_with(RXN_SMARTS_TAG) {
            return reaction_from_smarts(&text(), &ReactionSmartsParserParams::new());
        } else if key.starts_with(RXN_RXN_TAG) {
            return reaction_from_rxn_block(&text(), &MolFileParserParams::default());
        }
    }
    Err(RDKitError::FileParse("No suitable metadata found.".into()))
}
