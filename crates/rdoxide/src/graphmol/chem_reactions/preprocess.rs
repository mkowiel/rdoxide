//! Ported from RDKit `Code/GraphMol/ChemReactions/PreprocessRxn.cpp` and
//! `addRecursiveQueriesToReaction` of `Reaction.cpp`: reactant templates that
//! name functional groups (in an atom property) get the groups' queries.

use std::collections::BTreeMap;

use super::{ChemicalReaction, not_initialized};
use crate::Result;
use crate::graphmol::chem_transforms::recursive_queries::add_recursive_queries;
use crate::graphmol::functional_group_hierarchy::flattened_functional_group_hierarchy;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;

/// `(atom index, group name)` pairs per reactant template.
pub type ReactantLabels = Vec<Vec<(u32, String)>>;

/// The outcome of [`preprocess_reaction`].
#[derive(Debug, Clone, Default)]
pub struct PreprocessResult {
    /// whether the reaction validated (the queries were only added then)
    pub valid: bool,
    pub num_warnings: u32,
    pub num_errors: u32,
    pub reactant_labels: ReactantLabels,
}

/// `addRecursiveQueriesToReaction`: the reaction must be initialized.
pub fn add_recursive_queries_to_reaction(
    rxn: &mut ChemicalReaction,
    queries: &BTreeMap<String, ROMol>,
    prop_name: &str,
) -> Result<ReactantLabels> {
    if !rxn.is_initialized() {
        return Err(not_initialized());
    }
    rxn.reactants
        .iter_mut()
        .map(|template| add_recursive_queries(template, queries, prop_name))
        .collect()
}

/// `preprocessReaction(rxn, ..., queries, propName)`: initializes and
/// validates the reaction (as one with implicit properties), then adds the
/// queries if it is valid.
pub fn preprocess_reaction(
    rxn: &mut ChemicalReaction,
    queries: &BTreeMap<String, ROMol>,
    prop_name: &str,
) -> Result<PreprocessResult> {
    rxn.set_implicit_properties_flag(true);
    rxn.init_reactant_matchers(false);
    let (valid, num_warnings, num_errors) = rxn.validate(false);
    let reactant_labels = if valid {
        add_recursive_queries_to_reaction(rxn, queries, prop_name)?
    } else {
        Vec::new()
    };
    Ok(PreprocessResult {
        valid,
        num_warnings,
        num_errors,
        reactant_labels,
    })
}

/// `preprocessReaction(rxn, propName)` with the built-in functional groups.
pub fn preprocess_reaction_with_default_groups(
    rxn: &mut ChemicalReaction,
    prop_name: &str,
) -> Result<PreprocessResult> {
    preprocess_reaction(rxn, flattened_functional_group_hierarchy(true), prop_name)
}

/// The property holding the group names (`common_properties::molFileValue`).
pub const DEFAULT_PROP_NAME: &str = cp::molFileValue;
