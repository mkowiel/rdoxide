//! Reactions from CDXML documents.
//!
//! Ported from RDKit `Code/GraphMol/ChemReactions/CDXMLParser.cpp`: the
//! molecules of the document are grouped by scheme step and become the
//! reactant, agent and product templates of one reaction per step.

use std::collections::{BTreeMap, BTreeSet};

use super::ChemicalReaction;
use super::sanitize::{flags, match_only_at_rgroups_adjust_params, sanitize_rxn};
use super::utils::update_products_stereochem;
use crate::geometry::Point3D;
use crate::graphmol::chem_transforms::combine_mols;
use crate::graphmol::file_parsers::cdxml::{CdxmlParserParams, mols_from_cdxml};
use crate::graphmol::query_ops::replace_atom_with_query_atom;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

fn int_prop(mol: &ROMol, key: &str) -> Option<i64> {
    match mol.get_prop(key)? {
        RDValue::Int(i) => Some(i64::from(*i)),
        RDValue::UnsignedInt(u) => Some(i64::from(*u)),
        _ => None,
    }
}

fn make_query_atoms(mol: &mut RWMol) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        replace_atom_with_query_atom(mol, idx)?;
    }
    Ok(())
}

fn add_template(prop: &str, templates: &mut BTreeMap<i64, ROMol>, mol: RWMol) -> Result<()> {
    let idx =
        int_prop(&mol, prop).ok_or_else(|| RDKitError::Value(format!("Key Error: {prop}")))?;
    match templates.get(&idx) {
        Some(existing) => {
            let combined = combine_mols(existing, &mol, Point3D::new(0.0, 0.0, 0.0))?;
            templates.insert(idx, combined);
        }
        None => {
            templates.insert(idx, mol);
        }
    }
    Ok(())
}

/// `CDXMLToChemicalReactions`: one reaction per scheme step of the document.
pub fn cdxml_to_chemical_reactions(
    cdxml: &str,
    sanitize: bool,
    remove_hs: bool,
) -> Result<Vec<ChemicalReaction>> {
    let params = CdxmlParserParams {
        sanitize,
        remove_hs,
        ..Default::default()
    };
    let mut mols: Vec<Option<RWMol>> = mols_from_cdxml(cdxml, &params)?
        .into_iter()
        .map(Some)
        .collect();
    let mut schemes: BTreeMap<(i64, i64), Vec<usize>> = BTreeMap::new();
    for (i, mol) in mols.iter().enumerate() {
        let mol = mol.as_ref().expect("molecule");
        if let (Some(scheme), Some(step)) =
            (int_prop(mol, "CDX_SCHEME_ID"), int_prop(mol, "CDX_STEP_ID"))
        {
            schemes.entry((scheme, step)).or_default().push(i);
        }
    }
    let mut result = Vec::new();
    if schemes.is_empty() {
        return Ok(result);
    }
    let mut used: BTreeSet<usize> = BTreeSet::new();
    // as upstream, the templates are shared by all schemes of the document
    let mut reactant_templates: BTreeMap<i64, ROMol> = BTreeMap::new();
    let mut product_templates: BTreeMap<i64, ROMol> = BTreeMap::new();
    let mut agent_templates: BTreeMap<i64, ROMol> = BTreeMap::new();
    for indices in schemes.values() {
        let mut res = ChemicalReaction::new();
        for &idx in indices {
            if used.contains(&idx) {
                return Err(RDKitError::invariant(
                    "Fragment used in twice in one or more reactions, this shouldn't happen",
                ));
            }
            let Some(mol) = mols[idx].as_ref() else {
                continue;
            };
            let (prop, target) = if mol.has_prop("CDX_REAGENT_ID") {
                ("CDX_REAGENT_ID", &mut reactant_templates)
            } else if mol.has_prop("CDX_AGENT_ID") {
                ("CDX_AGENT_ID", &mut agent_templates)
            } else if mol.has_prop("CDX_PRODUCT_ID") {
                ("CDX_PRODUCT_ID", &mut product_templates)
            } else {
                continue;
            };
            used.insert(idx);
            let mut mol = mols[idx].take().expect("molecule");
            make_query_atoms(&mut mol)?;
            add_template(prop, target, mol)?;
        }
        for t in reactant_templates.values() {
            res.add_reactant_template(t.clone());
        }
        for t in agent_templates.values() {
            res.add_agent_template(t.clone());
        }
        for t in product_templates.values() {
            res.add_product_template(t.clone());
        }
        update_products_stereochem(&mut res)?;
        // CDXML-based reactions do not have implicit properties
        res.set_implicit_properties_flag(false);
        if !sanitize {
            // the reaction still has to be fixed up for SMARTS style matching
            let _ = sanitize_rxn(
                &mut res,
                flags::SANITIZE_ADJUST_REACTANTS | flags::SANITIZE_ADJUST_PRODUCTS,
                &match_only_at_rgroups_adjust_params(),
            );
        }
        result.push(res);
    }
    Ok(result)
}

/// `CDXMLFileToChemicalReactions`: an unreadable file gives no reactions.
pub fn cdxml_file_to_chemical_reactions(
    file_name: &str,
    sanitize: bool,
    remove_hs: bool,
) -> Result<Vec<ChemicalReaction>> {
    let Ok(bytes) = std::fs::read(file_name) else {
        return Ok(Vec::new());
    };
    if bytes.is_empty() {
        return Ok(Vec::new());
    }
    cdxml_to_chemical_reactions(&String::from_utf8_lossy(&bytes), sanitize, remove_hs)
}
