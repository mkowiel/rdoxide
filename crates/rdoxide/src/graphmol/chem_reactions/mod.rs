//! Ported from RDKit `Code/GraphMol/ChemReactions/` (`Reaction.{h,cpp}`,
//! `ReactionUtils.{h,cpp}`): reaction templates, their validation, and the
//! queries on reactions (`isMoleculeReactantOfReaction`, `getReactingAtoms`, …).
//!
//! Reactions are run by [`runner`], and read from reaction SMARTS/SMILES by
//! [`daylight`].
//!
//! Templates are owned by the reaction; where upstream hands out shared
//! pointers, the accessors return references (or copies at the Python
//! boundary).

pub mod daylight;
pub mod fingerprints;
pub mod mdl;
pub mod pickler;
pub mod preprocess;
pub mod runner;
pub mod sanitize;
pub mod utils;
pub mod writer;

use crate::graphmol::mol_props::get_avg_mol_wt;
use crate::graphmol::queries::{AtomQuery, MASS_INTEGER_CONVERSION_FACTOR};
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{Dict, RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

pub use utils::*;

fn not_initialized() -> RDKitError {
    RDKitError::ChemicalReaction("initReactantMatchers() must be called first".into())
}

/// `ChemicalReaction`
#[derive(Debug, Clone)]
pub struct ChemicalReaction {
    props: Dict,
    needs_init: bool,
    implicit_properties: bool,
    reactants: Vec<ROMol>,
    products: Vec<ROMol>,
    agents: Vec<ROMol>,
    substruct_params: SubstructMatchParameters,
}

impl Default for ChemicalReaction {
    fn default() -> Self {
        ChemicalReaction {
            props: Dict::new(),
            needs_init: true,
            implicit_properties: false,
            reactants: Vec::new(),
            products: Vec::new(),
            agents: Vec::new(),
            substruct_params: SubstructMatchParameters::default(),
        }
    }
}

impl RDProps for ChemicalReaction {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

impl ChemicalReaction {
    pub fn new() -> Self {
        Self::default()
    }

    /// `addReactantTemplate`; returns the number of reactant templates.
    pub fn add_reactant_template(&mut self, mol: ROMol) -> usize {
        self.needs_init = true;
        self.reactants.push(mol);
        self.reactants.len()
    }
    /// `addAgentTemplate`; returns the number of agent templates.
    pub fn add_agent_template(&mut self, mol: ROMol) -> usize {
        self.agents.push(mol);
        self.agents.len()
    }
    /// `addProductTemplate`; returns the number of product templates.
    pub fn add_product_template(&mut self, mol: ROMol) -> usize {
        self.products.push(mol);
        self.products.len()
    }

    pub fn reactants(&self) -> &[ROMol] {
        &self.reactants
    }
    pub fn products(&self) -> &[ROMol] {
        &self.products
    }
    pub fn agents(&self) -> &[ROMol] {
        &self.agents
    }
    pub fn reactants_mut(&mut self) -> &mut Vec<ROMol> {
        &mut self.reactants
    }
    pub fn products_mut(&mut self) -> &mut Vec<ROMol> {
        &mut self.products
    }
    pub fn agents_mut(&mut self) -> &mut Vec<ROMol> {
        &mut self.agents
    }
    /// All templates (reactants, agents, then products), mutably.
    pub fn templates_mut(&mut self) -> impl Iterator<Item = &mut ROMol> {
        self.reactants
            .iter_mut()
            .chain(self.agents.iter_mut())
            .chain(self.products.iter_mut())
    }
    /// The reactant and product templates (no agents), mutably.
    pub fn reactants_and_products_mut(&mut self) -> impl Iterator<Item = &mut ROMol> {
        self.reactants.iter_mut().chain(self.products.iter_mut())
    }
    pub fn num_reactant_templates(&self) -> usize {
        self.reactants.len()
    }
    pub fn num_product_templates(&self) -> usize {
        self.products.len()
    }
    pub fn num_agent_templates(&self) -> usize {
        self.agents.len()
    }

    pub fn is_initialized(&self) -> bool {
        !self.needs_init
    }
    /// `getImplicitPropertiesFlag`: whether unspecified atomic properties of
    /// the products are copied from the reactants (reaction-SMARTS
    /// semantics).
    pub fn implicit_properties_flag(&self) -> bool {
        self.implicit_properties
    }
    pub fn set_implicit_properties_flag(&mut self, val: bool) {
        self.implicit_properties = val;
    }
    pub fn substruct_params(&self) -> &SubstructMatchParameters {
        &self.substruct_params
    }
    pub fn substruct_params_mut(&mut self) -> &mut SubstructMatchParameters {
        &mut self.substruct_params
    }

    /// `initReactantMatchers`: validates the reaction; it can be run only
    /// when validation finds no errors.
    pub fn init_reactant_matchers(&mut self, silent: bool) {
        let (ok, _, _) = self.validate(silent);
        if !ok {
            rd_log_error!("initialization failed\n");
            self.needs_init = true;
        } else {
            self.needs_init = false;
        }
    }

    /// `validate`: `(valid, number of warnings, number of errors)`. Besides
    /// checking the atom maps, this records the charge, H count, mass and
    /// isotope queries of the product atoms as the implicit properties the
    /// runner uses.
    pub fn validate(&mut self, silent: bool) -> (bool, u32, u32) {
        let mut res = true;
        let mut warnings = 0;
        let mut errors = 0;
        if self.reactants.is_empty() {
            if !silent {
                rd_log_error!("reaction has no reactants\n");
            }
            errors += 1;
            res = false;
        }
        if self.products.is_empty() {
            if !silent {
                rd_log_error!("reaction has no products\n");
            }
            errors += 1;
            res = false;
        }

        // map number -> (reactant template, atom) of the reacting atoms
        let mut map_numbers_seen: Vec<i32> = Vec::new();
        let mut reacting_atoms: std::collections::BTreeMap<i32, (usize, u32)> = Default::default();
        for (mol_idx, mol) in self.reactants.iter().enumerate() {
            let mut this_mapped = false;
            for (ai, atom) in mol.atoms().iter().enumerate() {
                if let Some(map_num) = atom
                    .get_prop(cp::molAtomMapNumber)
                    .and_then(RDValue::as_int)
                {
                    this_mapped = true;
                    if map_numbers_seen.contains(&map_num) {
                        if !silent {
                            rd_log_error!(
                                "reactant atom-mapping number {map_num} found multiple times.\n"
                            );
                        }
                        errors += 1;
                        res = false;
                    } else {
                        map_numbers_seen.push(map_num);
                        reacting_atoms.insert(map_num, (mol_idx, ai as u32));
                    }
                }
            }
            if !this_mapped {
                if !silent {
                    rd_log_warning!("reactant {mol_idx} has no mapped atoms.\n");
                }
                warnings += 1;
            }
        }

        let mut product_numbers_seen: Vec<i32> = Vec::new();
        for mol_idx in 0..self.products.len() {
            for key in [
                cp::_QueryFormalCharge,
                cp::_QueryHCount,
                cp::_QueryMass,
                cp::_QueryIsotope,
            ] {
                for ai in 0..self.products[mol_idx].num_atoms() {
                    self.products[mol_idx].atom_mut(ai).clear_prop(key);
                }
            }
            let mut this_mapped = false;
            for ai in 0..self.products[mol_idx].num_atoms() {
                let map_num = self.products[mol_idx]
                    .atom(ai)
                    .get_prop(cp::molAtomMapNumber)
                    .and_then(RDValue::as_int);
                if let Some(map_num) = map_num {
                    this_mapped = true;
                    let seen_already = product_numbers_seen.contains(&map_num);
                    if seen_already {
                        if !silent {
                            rd_log_warning!(
                                "product atom-mapping number {map_num} found multiple times.\n"
                            );
                        }
                        warnings += 1;
                        if let Some(&(rm, ra)) = reacting_atoms.get(&map_num) {
                            if self.reactants[rm].atom_degree(ra)
                                != self.products[mol_idx].atom_degree(ai)
                            {
                                self.products[mol_idx]
                                    .atom_mut(ai)
                                    .dict_mut()
                                    .set_val(cp::_ReactionDegreeChanged, 1i32);
                            }
                        }
                    } else {
                        product_numbers_seen.push(map_num);
                    }
                    match map_numbers_seen.iter().position(|&m| m == map_num) {
                        None => {
                            if !seen_already {
                                if !silent {
                                    rd_log_warning!(
                                        "product atom-mapping number {map_num} not found in reactants.\n"
                                    );
                                }
                                warnings += 1;
                            }
                        }
                        Some(pos) => {
                            map_numbers_seen.remove(pos);
                            if let Some(&(rm, ra)) = reacting_atoms.get(&map_num) {
                                if self.reactants[rm].atom_degree(ra)
                                    != self.products[mol_idx].atom_degree(ai)
                                {
                                    self.products[mol_idx]
                                        .atom_mut(ai)
                                        .dict_mut()
                                        .set_val(cp::_ReactionDegreeChanged, 1i32);
                                }
                            }
                        }
                    }
                }
                warnings += self.record_query_properties(mol_idx, ai, silent);
            }
            if !this_mapped {
                if !silent {
                    rd_log_warning!("product {mol_idx} has no mapped atoms.\n");
                }
                warnings += 1;
            }
        }
        if !map_numbers_seen.is_empty() {
            if !silent {
                let numbers: String = map_numbers_seen.iter().map(|n| format!("{n} ")).collect();
                rd_log_warning!(
                    "mapped atoms in the reactants were not mapped in the products.\n  unmapped numbers are: {numbers}\n"
                );
            }
            warnings += 1;
        }
        (res, warnings, errors)
    }

    /// The query-walking part of `validate`: copies charge, H count, mass
    /// and isotope specifications of a product atom's query into properties.
    /// Returns the number of conflicting specifications (warnings).
    fn record_query_properties(&mut self, mol_idx: usize, atom_idx: u32, silent: bool) -> u32 {
        let Some(root) = self.products[mol_idx].atom(atom_idx).query() else {
            return 0;
        };
        // breadth first, as upstream's queue
        let mut queue: std::collections::VecDeque<&AtomQuery> = std::collections::VecDeque::new();
        queue.push_back(root);
        let mut found: Vec<(&'static str, RDValue)> = Vec::new();
        let mut warnings = 0;
        let mut seen: Vec<(&'static str, i64)> = Vec::new();
        while let Some(q) = queue.pop_front() {
            queue.extend(q.children.iter());
            let (key, what, value): (&'static str, &str, i64) = match q.description() {
                "AtomFormalCharge" => (cp::_QueryFormalCharge, "charge", i64::from(q.val)),
                "AtomNegativeFormalCharge" => (cp::_QueryFormalCharge, "charge", -i64::from(q.val)),
                "AtomHCount" => (cp::_QueryHCount, "H count", i64::from(q.val)),
                "AtomMass" => (cp::_QueryMass, "mass", i64::from(q.val)),
                "AtomIsotope" => (cp::_QueryIsotope, "isotope", i64::from(q.val)),
                _ => continue,
            };
            // the first specification wins; a different later one is a warning
            let previous = seen.iter().find(|(k, _)| *k == key).map(|&(_, v)| v);
            let check = if key == cp::_QueryFormalCharge {
                // upstream compares the sign-adjusted value with the raw one
                let neg = if q.description() == "AtomNegativeFormalCharge" {
                    -1
                } else {
                    1
                };
                previous.map(|p| neg * p != i64::from(q.val))
            } else {
                previous.map(|p| p != i64::from(q.val))
            };
            if check == Some(true) {
                if !silent {
                    rd_log_warning!(
                        "atom {atom_idx} in product {mol_idx} has multiple {what} specifications.\n"
                    );
                }
                warnings += 1;
            } else {
                seen.push((key, value));
                let rd: RDValue = if key == cp::_QueryMass {
                    RDValue::Double(f64::from(q.val) / f64::from(MASS_INTEGER_CONVERSION_FACTOR))
                } else {
                    RDValue::Int(value as i32)
                };
                found.push((key, rd));
            }
        }
        for (key, value) in found {
            self.products[mol_idx]
                .atom_mut(atom_idx)
                .dict_mut()
                .set_val(key, value);
        }
        warnings
    }

    /// `removeUnmappedReactantTemplates`: templates with (almost) no mapped
    /// atoms are agents; they move to the agent templates and/or `target`.
    pub fn remove_unmapped_reactant_templates(
        &mut self,
        threshold_unmapped_atoms: f64,
        move_to_agent_templates: bool,
        target: Option<&mut Vec<ROMol>>,
    ) {
        let (kept, removed) = split_agents(
            std::mem::take(&mut self.reactants),
            threshold_unmapped_atoms,
        );
        self.reactants = kept;
        self.take_removed(removed, move_to_agent_templates, target);
    }

    /// `removeUnmappedProductTemplates`
    pub fn remove_unmapped_product_templates(
        &mut self,
        threshold_unmapped_atoms: f64,
        move_to_agent_templates: bool,
        target: Option<&mut Vec<ROMol>>,
    ) {
        let (kept, removed) =
            split_agents(std::mem::take(&mut self.products), threshold_unmapped_atoms);
        self.products = kept;
        self.take_removed(removed, move_to_agent_templates, target);
    }

    fn take_removed(
        &mut self,
        removed: Vec<ROMol>,
        to_agents: bool,
        target: Option<&mut Vec<ROMol>>,
    ) {
        if to_agents {
            self.agents.extend(removed.iter().cloned());
        }
        if let Some(t) = target {
            t.extend(removed);
        }
    }

    /// `removeAgentTemplates`
    pub fn remove_agent_templates(&mut self, target: Option<&mut Vec<ROMol>>) {
        if let Some(t) = target {
            t.extend(self.agents.iter().cloned());
        }
        self.agents.clear();
    }

    /// `isMoleculeReactantOfReaction`: the reactant templates `mol` matches
    /// (only the first one with `stop_at_first_match`).
    pub fn reactant_templates_matching(
        &self,
        mol: &ROMol,
        stop_at_first_match: bool,
    ) -> Result<Vec<usize>> {
        self.templates_matching(&self.reactants, mol, stop_at_first_match)
    }

    /// `isMoleculeProductOfReaction`
    pub fn product_templates_matching(
        &self,
        mol: &ROMol,
        stop_at_first_match: bool,
    ) -> Result<Vec<usize>> {
        self.templates_matching(&self.products, mol, stop_at_first_match)
    }

    fn templates_matching(
        &self,
        templates: &[ROMol],
        mol: &ROMol,
        stop_at_first: bool,
    ) -> Result<Vec<usize>> {
        if !self.is_initialized() {
            return Err(not_initialized());
        }
        let mut which = Vec::new();
        for (i, t) in templates.iter().enumerate() {
            if !substruct_match(mol, t, &self.substruct_params)?.is_empty() {
                which.push(i);
                if stop_at_first {
                    break;
                }
            }
        }
        Ok(which)
    }

    /// `isMoleculeAgentOfReaction`: the index of the first agent template
    /// that `mol` matches (templates of a different size or weight are
    /// skipped quickly).
    pub fn agent_template_matching(&self, mol: &ROMol) -> Result<Option<usize>> {
        if !self.is_initialized() {
            return Err(not_initialized());
        }
        for (i, templ) in self.agents.iter().enumerate() {
            if templ.num_heavy_atoms() != mol.num_heavy_atoms()
                || templ.num_bonds() != mol.num_bonds()
                || get_avg_mol_wt(templ, false)? != get_avg_mol_wt(mol, false)?
            {
                continue;
            }
            if !substruct_match(mol, templ, &self.substruct_params)?.is_empty() {
                return Ok(Some(i));
            }
        }
        Ok(None)
    }
}

/// Splits templates into those that take part in the reaction and agents.
fn split_agents(templates: Vec<ROMol>, threshold: f64) -> (Vec<ROMol>, Vec<ROMol>) {
    templates
        .into_iter()
        .partition(|t| !is_reaction_template_molecule_agent(t, threshold))
}
