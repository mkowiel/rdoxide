//! Coordinate templates for ring systems (`Templates.{h,cpp}` of RDKit
//! 2026.09.1): 2D layouts, given as CXSMARTS, that complex ring systems are
//! matched against so that they are drawn the way chemists draw them.
//!
//! The default templates are the 578 of `TemplateSmarts.h`
//! (`template_smarts.txt`, one per line). The templates in use are a process
//! wide singleton, which can be replaced or extended from a file of CXSMILES.
//!
//! Copyright (C) 2023 Schrödinger, LLC.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock, RwLock};

use super::depict_error;
use crate::Result;
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, get_mol_frags, symmetrize_sssr};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::rdgeneral::threads::par_map;

/// `TEMPLATE_SMARTS`
const DEFAULT_TEMPLATES: &str = include_str!("template_smarts.txt");

type TemplatesBySize = BTreeMap<u32, Vec<Arc<ROMol>>>;

/// The templates, by number of atoms.
#[derive(Debug, Default)]
pub struct CoordinateTemplates {
    templates: TemplatesBySize,
}

/// A template from its CXSMARTS, with the ring info of the depictor.
fn parse_template(smarts: &str) -> Option<ROMol> {
    let mut mol = mol_from_smarts(smarts, &SmartsParserParams::default())
        .ok()
        .flatten()?;
    // initialize the ring info using symmetrizeSSSR to match the depictor's
    // ring counting
    symmetrize_sssr(
        &mut mol,
        SymmetrizeSSSRAlgorithm::Default,
        true,
        false,
        false,
    )
    .ok()?;
    Some(mol)
}

impl CoordinateTemplates {
    /// The singleton, which starts with the default templates.
    pub fn ring_system_templates() -> &'static RwLock<CoordinateTemplates> {
        static TEMPLATES: OnceLock<RwLock<CoordinateTemplates>> = OnceLock::new();
        TEMPLATES.get_or_init(|| {
            let mut t = CoordinateTemplates::default();
            t.load_default_templates();
            RwLock::new(t)
        })
    }

    pub fn has_template_of_size(&self, atom_count: u32) -> bool {
        self.templates.contains_key(&atom_count)
    }

    pub fn matching_templates(&self, atom_count: u32) -> Vec<Arc<ROMol>> {
        self.templates.get(&atom_count).cloned().unwrap_or_default()
    }

    /// Replaces the templates with the default ones.
    pub fn load_default_templates(&mut self) {
        self.templates.clear();
        // the templates are independent, so they parse on all cores; the
        // order within a size (which decides the template used) is kept
        let lines: Vec<&str> = DEFAULT_TEMPLATES.lines().collect();
        for mol in par_map(&lines, 0, |smarts| parse_template(smarts))
            .into_iter()
            .flatten()
        {
            self.templates
                .entry(mol.num_atoms())
                .or_default()
                .push(Arc::new(mol));
        }
    }

    /// `assertValidTemplate`: a template needs 2D coordinates and must be a
    /// single ring system (spiro ring systems are fine).
    pub fn assert_valid_template(mol: &ROMol, smiles: &str) -> Result<()> {
        if mol.num_conformers() == 0 {
            return Err(depict_error(format!(
                "Template missing coordinates: {smiles}"
            )));
        }
        if mol.conformer(-1)?.is_3d() {
            return Err(depict_error(format!(
                "Template has 3D coordinates, 2D coordinates required: {smiles}"
            )));
        }
        // every bond must be in a ring and there must be one connected
        // component in the molecular graph
        if get_mol_frags(mol).len() != 1 {
            return Err(depict_error(format!(
                "Template consists of multiple fragments, single fragment required: {smiles}"
            )));
        }
        if mol.num_atoms() == 1 {
            return Err(depict_error(format!(
                "Template is not a ring system: {smiles}"
            )));
        }
        for bond in 0..mol.num_bonds() {
            if mol.ring_info().num_bond_rings(bond)? == 0 {
                return Err(depict_error(format!(
                    "Template is not a ring system: {smiles}"
                )));
            }
        }
        Ok(())
    }

    /// `loadTemplatesFromPath`: the templates of a file with one CXSMILES per line.
    fn load_templates_from_path(template_path: &str) -> Result<TemplatesBySize> {
        let content = std::fs::read_to_string(template_path)
            .map_err(|_| depict_error(format!("Could not open file {template_path}")))?;
        let mut lines: Vec<&str> = content.split('\n').collect();
        if lines.last() == Some(&"") {
            lines.pop(); // a final newline does not start another line
        }
        let mut templates = TemplatesBySize::new();
        let parsed = par_map(&lines, 0, |line| parse_template(line));
        for (line, mol) in lines.into_iter().zip(parsed) {
            let mol = mol.ok_or_else(|| {
                depict_error(format!(
                    "Could not load templates from {template_path}: Invalid smarts"
                ))
            })?;
            Self::assert_valid_template(&mol, line)?;
            templates
                .entry(mol.num_atoms())
                .or_default()
                .push(Arc::new(mol));
        }
        Ok(templates)
    }

    /// Replaces the templates with the ones in a file. If any of them is
    /// invalid the current templates are kept.
    pub fn set_ring_system_templates(&mut self, template_path: &str) -> Result<()> {
        self.templates = Self::load_templates_from_path(template_path)?;
        Ok(())
    }

    /// Adds the templates in a file, which take precedence over the existing
    /// ones. If any of them is invalid the current templates are kept.
    pub fn add_ring_system_templates(&mut self, template_path: &str) -> Result<()> {
        for (size, added) in Self::load_templates_from_path(template_path)? {
            let existing = self.templates.entry(size).or_default();
            existing.splice(0..0, added);
        }
        Ok(())
    }
}

/// `RDDepict::setRingSystemTemplates`
pub fn set_ring_system_templates(template_path: &str) -> Result<()> {
    CoordinateTemplates::ring_system_templates()
        .write()
        .unwrap_or_else(|e| e.into_inner())
        .set_ring_system_templates(template_path)
}

/// `RDDepict::addRingSystemTemplates`
pub fn add_ring_system_templates(template_path: &str) -> Result<()> {
    CoordinateTemplates::ring_system_templates()
        .write()
        .unwrap_or_else(|e| e.into_inner())
        .add_ring_system_templates(template_path)
}

/// `RDDepict::loadDefaultRingSystemTemplates`
pub fn load_default_ring_system_templates() {
    CoordinateTemplates::ring_system_templates()
        .write()
        .unwrap_or_else(|e| e.into_inner())
        .load_default_templates();
}
