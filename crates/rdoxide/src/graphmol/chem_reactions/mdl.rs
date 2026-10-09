//! Ported from RDKit `Code/GraphMol/ChemReactions/` (`MDLParser.cpp`,
//! `MoleculeParser.cpp`, and the rxn parts of `ReactionWriter.cpp`): MDL
//! rxn files (V2000 and V3000) and the single-molecule encoding of a
//! reaction (`molRxnRole`).

use super::{ChemicalReaction, ReactionMoleculeType, templates, update_products_stereochem};
use crate::graphmol::file_parsers::mol_file_parser::{
    MolFileParserParams, finish_mol_processing, mol_from_mol_data_stream,
};
use crate::graphmol::file_parsers::mol_file_v3000::{get_v3000_line, parse_v3000_ctab};
use crate::graphmol::file_parsers::mol_file_writer::{
    MolWriterParams, get_v3000_ctab, mol_to_mol_block,
};
use crate::graphmol::file_parsers::sgroup_writing::move_additional_properties_to_sgroups;
use crate::graphmol::file_parsers::utils::LineStream;
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::query_ops::replace_atom_with_query_atom;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result, rd_log_warning};

fn parser_error(msg: impl Into<String>) -> RDKitError {
    RDKitError::ChemicalReactionParser(msg.into())
}

fn premature_eof() -> RDKitError {
    parser_error("premature EOF hit.")
}

fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// `FileParserUtils::stripSpacesAndCast<unsigned int>`
fn strip_and_cast(text: &[u8]) -> Result<u32> {
    let s = lossy(text);
    s.trim()
        .parse()
        .map_err(|_| parser_error(format!("Cannot convert {s} to int")))
}

/// Wraps a failure to read template `i` of the given kind like upstream does.
fn template_error(kind: &str, i: u32, err: RDKitError) -> RDKitError {
    match err {
        RDKitError::FileParse(m) => {
            parser_error(format!("Cannot parse {kind} {i}. The error was:\n\t{m}"))
        }
        other => other,
    }
}

fn file_error(msg: &str) -> RDKitError {
    RDKitError::FileParse(msg.to_string())
}

/// Reads one `$MOL` section (V2000 reaction files).
fn read_v2000_template(
    s: &mut LineStream<'_>,
    line: &mut u32,
    params: &MolFileParserParams,
    kind: &str,
    i: u32,
) -> Result<Option<RWMol>> {
    *line += 1;
    let header = s.get_line();
    if s.eof() {
        return Err(premature_eof());
    }
    if !header.starts_with(b"$MOL") {
        return Err(parser_error("$MOL header not found"));
    }
    mol_from_mol_data_stream(s, line, params).map_err(|e| template_error(kind, i, e))
}

fn parse_v2000_rxn_block(
    s: &mut LineStream<'_>,
    line: &mut u32,
    params: &MolFileParserParams,
) -> Result<ChemicalReaction> {
    // the name and comment fields are skipped
    let mut counts = Vec::new();
    for _ in 0..4 {
        *line += 1;
        counts = s.get_line();
    }
    if s.eof() {
        return Err(premature_eof());
    }
    if counts.len() < 6 {
        return Err(parser_error("rxn counts line is too short"));
    }
    let n_reacts = strip_and_cast(&counts[0..3])?;
    let n_prods = strip_and_cast(&counts[3..6])?;
    let mut n_agents = 0;
    if counts.len() > 6 {
        let field = &counts[6..counts.len().min(9)];
        if !lossy(field).trim().is_empty() {
            n_agents = strip_and_cast(field)?;
        }
    }
    let mut rxn = ChemicalReaction::new();
    for i in 0..n_reacts {
        let mol = read_v2000_template(s, line, params, "reactant", i)?
            .ok_or_else(|| parser_error("Null reactant in reaction file."))?;
        rxn.add_reactant_template(mol);
    }
    for i in 0..n_prods {
        let mol = read_v2000_template(s, line, params, "product", i)?
            .ok_or_else(|| parser_error("Null product in reaction file."))?;
        rxn.add_product_template(mol);
    }
    // agents are neither sanitized nor stripped of Hs
    let agent_params = MolFileParserParams {
        sanitize: false,
        remove_hs: false,
        strict_parsing: params.strict_parsing,
        ..MolFileParserParams::default()
    };
    for i in 0..n_agents {
        let mol = read_v2000_template(s, line, &agent_params, "agent", i)?
            .ok_or_else(|| parser_error("Null agent in reaction file."))?;
        rxn.add_agent_template(mol);
    }
    Ok(rxn)
}

/// Reads a V3000 line and checks that it starts with `expected` (case-insensitively).
fn expect_v3000_line(s: &mut LineStream<'_>, line: &mut u32, expected: &str) -> Result<()> {
    let text = lossy(&get_v3000_line(s, line)?).to_uppercase();
    if text.starts_with(expected) {
        Ok(())
    } else {
        Err(file_error(&format!("{expected} line not found")))
    }
}

fn read_v3000_template(
    s: &mut LineStream<'_>,
    line: &mut u32,
    params: &MolFileParserParams,
    strict: bool,
    kind: &str,
    i: u32,
) -> Result<RWMol> {
    let mut mol = RWMol::new();
    let mut chirality_possible = false;
    parse_v3000_ctab(
        s,
        line,
        &mut mol,
        &mut chirality_possible,
        strict,
        false,
        false,
    )
    .and_then(|_| finish_mol_processing(&mut mol, chirality_possible, params))
    .map_err(|e| template_error(kind, i, e))?;
    Ok(mol)
}

fn parse_v3000_rxn_block(
    s: &mut LineStream<'_>,
    line: &mut u32,
    params: &MolFileParserParams,
) -> Result<ChemicalReaction> {
    // skip the header block
    for _ in 0..3 {
        *line += 1;
        s.get_line();
    }
    let counts = lossy(&get_v3000_line(s, line)?).to_uppercase();
    let tokens: Vec<&str> = counts
        .split([' ', '\t'])
        .filter(|t| !t.is_empty())
        .collect();
    if tokens.len() < 3 || tokens[0] != "COUNTS" {
        return Err(parser_error("bad counts line"));
    }
    let n_reacts = strip_and_cast(tokens[1].as_bytes())?;
    let n_prods = strip_and_cast(tokens[2].as_bytes())?;
    let n_agents = match tokens.get(3) {
        Some(t) => strip_and_cast(t.as_bytes())?,
        None => 0,
    };

    let mut rxn = ChemicalReaction::new();
    expect_v3000_line(s, line, "BEGIN REACTANT")?;
    for i in 0..n_reacts {
        let mol = read_v3000_template(s, line, params, params.strict_parsing, "reactant", i)?;
        rxn.add_reactant_template(mol);
    }
    expect_v3000_line(s, line, "END REACTANT")?;
    expect_v3000_line(s, line, "BEGIN PRODUCT")?;
    for i in 0..n_prods {
        let mol = read_v3000_template(s, line, params, params.strict_parsing, "product", i)?;
        rxn.add_product_template(mol);
    }
    expect_v3000_line(s, line, "END PRODUCT")?;
    if n_agents > 0 {
        expect_v3000_line(s, line, "BEGIN AGENT")?;
        for i in 0..n_agents {
            let mol = read_v3000_template(s, line, params, true, "agent", i)?;
            rxn.add_agent_template(mol);
        }
        expect_v3000_line(s, line, "END AGENT")?;
    }
    Ok(rxn)
}

/// `v2::ReactionParser::ReactionFromRxnDataStream`
pub fn reaction_from_rxn_data_stream(
    s: &mut LineStream<'_>,
    line: &mut u32,
    params: &MolFileParserParams,
) -> Result<ChemicalReaction> {
    *line += 1;
    let header = s.get_line();
    if s.eof() {
        return Err(premature_eof());
    }
    if !header.starts_with(b"$RXN") {
        return Err(parser_error("$RXN header not found"));
    }
    let v3000 = header.len() >= 10 && &header[5..10] == b"V3000";
    let mut rxn = if v3000 {
        parse_v3000_rxn_block(s, line, params)?
    } else {
        parse_v2000_rxn_block(s, line, params)?
    };
    // the matching needs query atoms
    for mol in rxn.reactants_and_products_mut() {
        for idx in 0..mol.num_atoms() {
            replace_atom_with_query_atom(mol, idx)?;
        }
    }
    update_products_stereochem(&mut rxn)?;
    // RXN-based reactions do not have implicit properties
    rxn.set_implicit_properties_flag(false);
    Ok(rxn)
}

/// `v2::ReactionParser::ReactionFromRxnBlock`
pub fn reaction_from_rxn_block(
    block: &str,
    params: &MolFileParserParams,
) -> Result<ChemicalReaction> {
    let mut s = LineStream::new(block.as_bytes());
    reaction_from_rxn_data_stream(&mut s, &mut 0, params)
}

/// `v2::ReactionParser::ReactionFromRxnFile`: `Ok(None)` when the file cannot be read.
pub fn reaction_from_rxn_file(
    file_name: &str,
    params: &MolFileParserParams,
) -> Result<Option<ChemicalReaction>> {
    let Ok(data) = std::fs::read(file_name) else {
        return Ok(None);
    };
    if data.is_empty() {
        return Ok(None);
    }
    let mut s = LineStream::new(&data);
    reaction_from_rxn_data_stream(&mut s, &mut 0, params).map(Some)
}

// ------------------------------------------------------------------ writing

fn v3000_template(mol: &ROMol) -> Result<String> {
    let mut copy = mol.clone();
    if copy.needs_update_property_cache() {
        copy.update_property_cache(false)?;
    }
    move_additional_properties_to_sgroups(&mut copy)?;
    let was_aromatic = vec![false; copy.num_bonds() as usize];
    get_v3000_ctab(&mut copy, &was_aromatic, -1, 6)
}

/// `ChemicalReactionToV3KRxnBlock`
pub fn reaction_to_v3k_rxn_block(rxn: &ChemicalReaction, separate_agents: bool) -> Result<String> {
    let (reactants, agents, products) = (rxn.reactants(), rxn.agents(), rxn.products());
    let mut res = String::from("$RXN V3000\n\n      RDKit\n\n");
    if separate_agents {
        res += &format!(
            "M  V30 COUNTS {} {} {}\n",
            reactants.len(),
            products.len(),
            agents.len()
        );
    } else {
        res += &format!(
            "M  V30 COUNTS {} {}\n",
            reactants.len() + agents.len(),
            products.len()
        );
    }
    let mut section = |name: &str, mols: &[&[ROMol]]| -> Result<()> {
        res += &format!("M  V30 BEGIN {name}\n");
        for mol in mols.iter().flat_map(|group| group.iter()) {
            res += &v3000_template(mol)?;
        }
        res += &format!("M  V30 END {name}\n");
        Ok(())
    };
    if separate_agents {
        section("REACTANT", &[reactants])?;
        section("PRODUCT", &[products])?;
        section("AGENT", &[agents])?;
    } else {
        section("REACTANT", &[reactants, agents])?;
        section("PRODUCT", &[products])?;
    }
    res += "M  END\n";
    Ok(res)
}

/// `ChemicalReactionToRxnBlock`
pub fn reaction_to_rxn_block(
    rxn: &ChemicalReaction,
    separate_agents: bool,
    force_v3000: bool,
) -> Result<String> {
    if force_v3000 {
        return reaction_to_v3k_rxn_block(rxn, separate_agents);
    }
    let (reactants, agents, products) = (rxn.reactants(), rxn.agents(), rxn.products());
    let mut res = String::from("$RXN\n\n      RDKit\n\n");
    if separate_agents {
        res += &format!(
            "{:3}{:3}{:3}\n",
            reactants.len(),
            products.len(),
            agents.len()
        );
    } else {
        res += &format!("{:3}{:3}\n", reactants.len() + agents.len(), products.len());
    }
    let params = MolWriterParams {
        kekulize: false,
        ..MolWriterParams::default()
    };
    let groups: [&[ROMol]; 3] = if separate_agents {
        [reactants, products, agents]
    } else {
        [reactants, agents, products]
    };
    for mol in groups.iter().flat_map(|group| group.iter()) {
        res += "$MOL\n";
        res += &mol_to_mol_block(mol, &params, -1)?;
    }
    Ok(res)
}

// ------------------------------------------------------ single-molecule form

fn rxn_role(atom_roles: impl Iterator<Item = Option<i32>>) -> (i32, bool) {
    let mut first = None;
    let mut uniform = true;
    for role in atom_roles.flatten() {
        match first {
            None => first = Some(role),
            Some(f) if f != role => uniform = false,
            _ => {}
        }
    }
    (first.unwrap_or(-1), uniform)
}

/// `RxnMolToChemicalReaction`: builds a reaction from a molecule whose atoms
/// carry a `molRxnRole` (1 reactant, 2 product, 3 agent); each fragment
/// becomes one template.
pub fn reaction_from_rxn_mol(mol: &ROMol) -> Result<ChemicalReaction> {
    let mut rxn = ChemicalReaction::new();
    let (fragments, _) = get_mol_frags_as_mols(mol, true, true)?;
    for (i, frag) in fragments.into_iter().enumerate() {
        let roles = frag
            .atoms()
            .iter()
            .map(|a| a.get_prop(cp::molRxnRole).and_then(|v| v.as_int()));
        let (role, uniform) = rxn_role(roles);
        if !uniform {
            rd_log_warning!(">> Atoms within one molecule have different RXN roles.\n");
            continue;
        }
        match role {
            1 => rxn.add_reactant_template(frag),
            2 => rxn.add_product_template(frag),
            3 => rxn.add_agent_template(frag),
            _ => {
                rd_log_warning!(
                    ">> Fragment {i} not included in the reaction, atoms do not have a correct RXN role.\n"
                );
                continue;
            }
        };
    }
    Ok(rxn)
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
