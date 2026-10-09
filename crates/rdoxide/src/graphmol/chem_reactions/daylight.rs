//! Ported from RDKit `Code/GraphMol/ChemReactions/DaylightParser.cpp`:
//! reading reactions from reaction SMARTS and reaction SMILES.

use std::collections::BTreeMap;

use super::{ChemicalReaction, update_products_stereochem};
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::cx::parse_cx_extensions_at;
use crate::graphmol::smiles_parse::ops::cleanup_after_parsing;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::{RDKitError, Result};

/// `v2::ReactionParser::ReactionSmartsParserParams`
#[derive(Debug, Clone, Default)]
pub struct ReactionSmartsParserParams {
    /// sanitize the molecules after building them (reaction SMILES only)
    pub sanitize: bool,
    /// SMILES "macros"
    pub replacements: BTreeMap<String, String>,
    /// recognize and parse CXSMILES
    pub allow_cx_smiles: bool,
    /// fail if the CXSMILES parsing fails
    pub strict_cx_smiles: bool,
}

impl ReactionSmartsParserParams {
    /// The upstream defaults (CXSMILES recognized and strict).
    pub fn new() -> Self {
        ReactionSmartsParserParams {
            allow_cx_smiles: true,
            strict_cx_smiles: true,
            ..Default::default()
        }
    }
}

fn parser_error(msg: impl Into<String>) -> RDKitError {
    RDKitError::ChemicalReactionParser(msg.into())
}

/// `DaylightParserUtils::splitSmartsIntoComponents`: splits on top-level
/// dots, unwrapping fully parenthesized components so that `(A.B)` stays one
/// (intramolecular) template.
pub fn split_smarts_into_components(react_text: &str) -> Vec<String> {
    let bytes = react_text.as_bytes();
    let mut res = Vec::new();
    let mut block_start = 0usize;
    let mut level = 0u32;
    // 0: not a block, 1: inside a leading parenthesis, 2: block closed
    let mut in_block = 0u8;
    let element = |start: usize, end: usize, in_block: u8| -> String {
        if in_block == 2 {
            react_text[start + 1..end - 1].to_string()
        } else {
            react_text[start..end].to_string()
        }
    };
    for (pos, &c) in bytes.iter().enumerate() {
        match c {
            b'(' => {
                if pos == block_start {
                    in_block = 1;
                }
                level += 1;
            }
            b')' => {
                if level == 1 && in_block != 0 {
                    in_block = 2;
                }
                level = level.saturating_sub(1);
            }
            b'.' if level == 0 => {
                res.push(element(block_start, pos, in_block));
                block_start = pos + 1;
                in_block = 0;
            }
            _ => {}
        }
    }
    if block_start < bytes.len() {
        res.push(element(block_start, bytes.len(), in_block));
    }
    res
}

fn construct_mol_from_string(
    txt: &str,
    params: &ReactionSmartsParserParams,
    use_smiles: bool,
) -> Result<Option<ROMol>> {
    if use_smiles {
        mol_from_smiles(
            txt,
            &SmilesParserParams {
                replacements: params.replacements.clone(),
                allow_cx_smiles: false,
                parse_name: false,
                sanitize: params.sanitize,
                remove_hs: false,
                skip_cleanup: true,
                ..Default::default()
            },
        )
    } else {
        mol_from_smarts(
            txt,
            &SmartsParserParams {
                replacements: params.replacements.clone(),
                allow_cx_smiles: false,
                parse_name: false,
                merge_hs: false,
                skip_cleanup: true,
                ..Default::default()
            },
        )
    }
}

fn is_space(c: u8) -> bool {
    c == b' ' || c == b'\t'
}

/// `removeSpacesAround`: drops blanks on both sides of `text[pos]`.
fn remove_spaces_around(text: &mut Vec<u8>, pos: usize) {
    let next = pos + 1;
    while next < text.len() && is_space(text[next]) {
        text.remove(next);
    }
    let mut prev = pos;
    while prev > 0 && is_space(text[prev - 1]) {
        text.remove(prev - 1);
        prev -= 1;
    }
}

/// Positions of the `>` characters that are not part of a `->` arrow.
fn arrow_positions(text: &[u8]) -> Vec<usize> {
    (0..text.len())
        .filter(|&i| text[i] == b'>' && (i == 0 || text[i - 1] != b'-'))
        .collect()
}

fn parse_reaction(
    orig_text: &str,
    params: &ReactionSmartsParserParams,
    use_smiles: bool,
) -> Result<ChemicalReaction> {
    let mut text = orig_text;
    let mut cx_part = "";
    if params.allow_cx_smiles {
        if let Some(sidx) = orig_text.find('|').filter(|&s| s != 0) {
            text = &orig_text[..sidx];
            cx_part = orig_text[sidx..].trim();
        }
    }
    let mut text: Vec<u8> = text.trim().as_bytes().to_vec();
    if arrow_positions(&text).len() < 2 {
        return Err(parser_error(
            "a reaction requires at least two > characters",
        ));
    }
    // remove spaces around ">" and "." symbols
    for symbol in *b">." {
        let positions: Vec<usize> = if symbol == b'>' {
            arrow_positions(&text)
        } else {
            (0..text.len()).filter(|&i| text[i] == b'.').collect()
        };
        for &p in positions.iter().rev() {
            remove_spaces_around(&mut text, p);
        }
    }
    // anything after the first remaining blank is a name
    if let Some(sidx) = text.iter().position(|&c| is_space(c)).filter(|&s| s != 0) {
        text.truncate(sidx);
    }
    let pos = arrow_positions(&text);
    if pos.len() < 2 {
        return Err(parser_error(
            "a reaction requires at least two > characters",
        ));
    }
    if pos.len() > 2 {
        return Err(parser_error("multi-step reactions not supported"));
    }
    let text = String::from_utf8_lossy(&text).into_owned();
    let (pos1, pos2) = (pos[0], pos[1]);
    let react_text = &text[..pos1];
    let agent_text = if pos2 != pos1 + 1 {
        &text[pos1 + 1..pos2]
    } else {
        ""
    };
    let product_text = &text[pos2 + 1..];

    let build = |txt: &str, what: &str| -> Result<ROMol> {
        construct_mol_from_string(txt, params, use_smiles)?
            .ok_or_else(|| parser_error(format!("Problems constructing {what} from SMARTS: {txt}")))
    };

    let mut rxn = ChemicalReaction::new();
    for txt in split_smarts_into_components(react_text) {
        rxn.add_reactant_template(build(&txt, "reactant")?);
    }
    for txt in split_smarts_into_components(product_text) {
        rxn.add_product_template(build(&txt, "product")?);
    }
    update_products_stereochem(&mut rxn)?;

    // a reaction template may have no agent
    if !agent_text.is_empty() {
        let agent_mol = build(agent_text, "agent")?;
        let (frags, _) = get_mol_frags_as_mols(&agent_mol, false, true)?;
        for agent in frags {
            rxn.add_agent_template(agent);
        }
    }

    if params.allow_cx_smiles && !cx_part.is_empty() {
        // the indices in the block run over all templates, in order
        let (mut start_atom, mut start_bond) = (0, 0);
        for mol in rxn.templates_mut() {
            let mut pos = 0;
            parse_cx_extensions_at(mol, cx_part, &mut pos, start_atom, start_bond).map_err(
                |_| {
                    RDKitError::Value(
                        "SmilesParseException: failure parsing CXSMILES extensions".into(),
                    )
                },
            )?;
            start_atom += mol.num_atoms();
            start_bond += mol.num_bonds();
        }
    }

    for mol in rxn.templates_mut() {
        cleanup_after_parsing(mol);
    }
    // "SMARTS"-based reactions have implicit properties
    rxn.set_implicit_properties_flag(true);
    Ok(rxn)
}

/// `v2::ReactionParser::ReactionFromSmarts`
pub fn reaction_from_smarts(
    text: &str,
    params: &ReactionSmartsParserParams,
) -> Result<ChemicalReaction> {
    parse_reaction(text, params, false)
}

/// `v2::ReactionParser::ReactionFromSmiles`
pub fn reaction_from_smiles(
    text: &str,
    params: &ReactionSmartsParserParams,
) -> Result<ChemicalReaction> {
    parse_reaction(text, params, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_components() {
        assert_eq!(
            split_smarts_into_components("[C:1].[N:2]"),
            ["[C:1]", "[N:2]"]
        );
        assert_eq!(
            split_smarts_into_components("([C:1].[N:2]).O"),
            ["[C:1].[N:2]", "O"]
        );
        assert_eq!(split_smarts_into_components(""), Vec::<String>::new());
    }

    #[test]
    fn parses_basic_reaction() {
        let rxn = reaction_from_smarts(
            "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]",
            &ReactionSmartsParserParams::new(),
        )
        .unwrap();
        assert_eq!(rxn.num_reactant_templates(), 2);
        assert_eq!(rxn.num_product_templates(), 1);
        assert!(rxn.implicit_properties_flag());
    }

    #[test]
    fn rejects_malformed() {
        let p = ReactionSmartsParserParams::new();
        assert!(reaction_from_smarts("CC>>", &p).is_ok());
        assert!(reaction_from_smarts("CC>", &p).is_err());
        assert!(reaction_from_smarts("C>C>C>C", &p).is_err());
    }
}
