//! Ported from RDKit `Code/GraphMol/ChemReactions/SanitizeRxn.{h,cpp}`:
//! clean-up of reactions read from drawings (R-group labels, atom maps,
//! template aromaticity and explicit hydrogens).

use std::collections::BTreeMap;

use super::ChemicalReaction;
use crate::graphmol::adjust_query::{
    AdjustQueryParameters, adjust_flags::*, adjust_query_properties_in_place,
};
use crate::graphmol::mol_ops::query_hs::merge_query_hs;
use crate::graphmol::mol_ops::{sanitize_flags::SANITIZE_SETAROMATICITY, sanitize_mol};
use crate::graphmol::queries::make_atom_null_query;
use crate::graphmol::query_ops::query_atom_from_atom;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result, rd_log_warning};

/// `RxnOps::SanitizeRxnFlags`
pub mod flags {
    pub const SANITIZE_NONE: u32 = 0x0;
    pub const SANITIZE_RGROUP_NAMES: u32 = 0x1;
    pub const SANITIZE_ATOM_MAPS: u32 = 0x2;
    pub const SANITIZE_ADJUST_REACTANTS: u32 = 0x4;
    pub const SANITIZE_MERGEHS: u32 = 0x8;
    pub const SANITIZE_ADJUST_PRODUCTS: u32 = 0x10;
    pub const SANITIZE_ALL: u32 = 0xFFFF_FFFF;
}
use flags::*;

/// `RxnSanitizeException` (a `RuntimeError` in Python)
fn sanitize_error(msg: impl Into<String>) -> RDKitError {
    RDKitError::Runtime(msg.into())
}

/// `RxnOps::DefaultRxnAdjustParams`: aromatizes Kekule templates if possible.
pub fn default_rxn_adjust_params() -> AdjustQueryParameters {
    AdjustQueryParameters {
        adjust_degree: false,
        adjust_degree_flags: ADJUST_IGNOREALL,
        adjust_ring_count: false,
        adjust_ring_count_flags: ADJUST_IGNOREALL,
        make_dummies_queries: false,
        aromatize_if_possible: true,
        ..AdjustQueryParameters::default()
    }
}

/// `RxnOps::MatchOnlyAtRgroupsAdjustParams`
pub fn match_only_at_rgroups_adjust_params() -> AdjustQueryParameters {
    AdjustQueryParameters {
        adjust_degree: true,
        adjust_degree_flags: ADJUST_IGNOREDUMMIES,
        adjust_ring_count: false,
        adjust_ring_count_flags: ADJUST_IGNORENONE,
        make_dummies_queries: false,
        aromatize_if_possible: true,
        ..AdjustQueryParameters::default()
    }
}

/// `RxnOps::ChemDrawRxnAdjustParams` (deprecated alias of the above).
pub fn chem_draw_rxn_adjust_params() -> AdjustQueryParameters {
    rd_log_warning!(" deprecated -- please use MatchOnlyAtRgroupsAdjustParams instead");
    match_only_at_rgroups_adjust_params()
}

fn int_prop(atom: &crate::graphmol::atom::Atom, key: &str) -> Option<i32> {
    atom.get_prop(key).and_then(|v| v.as_int())
}

fn uint_prop(atom: &crate::graphmol::atom::Atom, key: &str) -> Option<u32> {
    atom.get_prop(key).and_then(|v| v.as_uint())
}

/// The largest value of an atom property over all templates.
fn max_prop<T: Default + PartialOrd + Copy>(
    rxn: &ChemicalReaction,
    get: impl Fn(&crate::graphmol::atom::Atom) -> Option<T>,
) -> T {
    rxn.reactants
        .iter()
        .chain(&rxn.agents)
        .chain(&rxn.products)
        .flat_map(|m| m.atoms().iter().filter_map(&get))
        .fold(T::default(), |max, v| if v > max { v } else { max })
}

/// A snapshot of the labels of a template atom (`AtomInfo`).
#[derive(Debug)]
struct AtomInfo {
    template_idx: usize,
    atom_idx: u32,
    rlabel: u32,
    atom_map: i32,
    isotope: u32,
    dummy_label: String,
    /// set once the atom has been paired up
    used: bool,
}

impl AtomInfo {
    fn new(mol: &ROMol, template_idx: usize, atom_idx: u32) -> Self {
        let atom = mol.atom(atom_idx);
        AtomInfo {
            template_idx,
            atom_idx,
            rlabel: uint_prop(atom, cp::_MolFileRLabel).unwrap_or(0),
            atom_map: int_prop(atom, cp::molAtomMapNumber).unwrap_or(0),
            isotope: atom.isotope(),
            dummy_label: atom
                .get_prop(cp::dummyLabel)
                .map(|v| v.to_rd_string())
                .unwrap_or_default(),
            used: false,
        }
    }

    fn needs_rlabel(&self, mol: &ROMol) -> bool {
        mol.atom(self.atom_idx).atomic_num() == 0 && self.rlabel == 0
    }

    fn best_guess_rlabel(&self) -> u32 {
        if self.rlabel != 0 {
            return self.rlabel;
        }
        if self.isotope != 0 {
            return self.isotope;
        }
        if self.atom_map != 0 {
            return self.atom_map as u32;
        }
        self.dummy_label
            .get(1..)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0)
    }
}

fn set_rlabel(mol: &mut ROMol, atom_idx: u32, rlabel: u32) -> Result<()> {
    let mut qatom = query_atom_from_atom(mol.atom(atom_idx))?;
    qatom.dict_mut().set_val(cp::_MolFileRLabel, rlabel);
    qatom
        .dict_mut()
        .set_val(cp::dummyLabel, format!("R{rlabel}"));
    if rlabel > 0 && rlabel < 999 {
        qatom.set_isotope(rlabel);
    }
    qatom.set_query(Some(make_atom_null_query()));
    mol.replace_atom(atom_idx, &qatom, false)
}

fn set_atom_map(mol: &mut ROMol, atom_idx: u32, map: i32) {
    mol.atom_mut(atom_idx)
        .dict_mut()
        .set_val(cp::molAtomMapNumber, map);
}

fn atoms_needing_rlabels(templates: &[ROMol]) -> Vec<AtomInfo> {
    let mut res = Vec::new();
    for (ti, mol) in templates.iter().enumerate() {
        for idx in 0..mol.num_atoms() {
            let info = AtomInfo::new(mol, ti, idx);
            if info.needs_rlabel(mol) {
                res.push(info);
            }
        }
    }
    res
}

/// `RxnOps::fixRGroups`: dummy atoms with a map but no R-group label become
/// R groups when a matching dummy exists on the other side.
pub fn fix_rgroups(rxn: &mut ChemicalReaction) -> Result<()> {
    let reactant_atoms = atoms_needing_rlabels(&rxn.reactants);
    let mut product_atoms = atoms_needing_rlabels(&rxn.products);
    if reactant_atoms.is_empty() && product_atoms.is_empty() {
        return Ok(());
    }
    if reactant_atoms.len() > product_atoms.len() {
        rd_log_warning!(
            "Mismatched potential rlabels: {} unmapped reactant dummy atom rlabels,{} unmappped product dummy atom rlabels",
            reactant_atoms.len(),
            product_atoms.len()
        );
    }
    let max_rlabel: u32 = max_prop(rxn, |a| uint_prop(a, cp::_MolFileRLabel));
    let max_atom_map: i32 = max_prop(rxn, |a| int_prop(a, cp::molAtomMapNumber));

    for rat in &reactant_atoms {
        let guess = rat.best_guess_rlabel();
        if guess == 0 {
            continue;
        }
        let partner = product_atoms
            .iter_mut()
            .find(|p| !p.used && p.best_guess_rlabel() == guess && p.atom_map == rat.atom_map);
        let Some(pat) = partner else {
            rd_log_warning!(
                "Could not find RLabel mapping for atom: {} in template: {}",
                rat.atom_idx,
                rat.template_idx
            );
            continue;
        };
        set_rlabel(
            &mut rxn.reactants[rat.template_idx],
            rat.atom_idx,
            max_rlabel + guess,
        )?;
        set_rlabel(
            &mut rxn.products[pat.template_idx],
            pat.atom_idx,
            max_rlabel + guess,
        )?;
        if rat.atom_map == 0 {
            let map = max_atom_map + guess as i32;
            set_atom_map(&mut rxn.reactants[rat.template_idx], rat.atom_idx, map);
            set_atom_map(&mut rxn.products[pat.template_idx], pat.atom_idx, map);
        }
        pat.used = true;
    }
    Ok(())
}

/// `RxnOps::fixAtomMaps`: gives unmapped R-group atoms an atom map derived
/// from their label and checks the products agree.
pub fn fix_atom_maps(rxn: &mut ChemicalReaction) -> Result<()> {
    let max_atom_map: i32 = max_prop(rxn, |a| int_prop(a, cp::molAtomMapNumber));
    let mut potential: BTreeMap<u32, i32> = BTreeMap::new();

    for mol in &mut rxn.reactants {
        for idx in 0..mol.num_atoms() {
            let info = AtomInfo::new(mol, 0, idx);
            if info.rlabel != 0 && info.atom_map == 0 {
                if potential.contains_key(&info.rlabel) {
                    return Err(sanitize_error("Duplicated RLabels"));
                }
                let map = info.rlabel as i32 + max_atom_map;
                potential.insert(info.rlabel, map);
                set_atom_map(mol, idx, map);
            }
        }
    }
    if potential.is_empty() {
        return Ok(());
    }
    for (ti, mol) in rxn.products.iter_mut().enumerate() {
        for idx in 0..mol.num_atoms() {
            let info = AtomInfo::new(mol, ti, idx);
            if info.rlabel == 0 {
                continue;
            }
            let expected = potential.get(&info.rlabel).copied().unwrap_or(0);
            if info.atom_map == 0 {
                set_atom_map(mol, idx, expected);
            } else if info.atom_map != expected {
                return Err(sanitize_error(format!(
                    "RLabel is mapped in product but not in reactant for product idx: {} atom: {}",
                    info.template_idx, info.atom_idx
                )));
            }
        }
    }
    Ok(())
}

/// `RxnOps::fixHs`: merges query hydrogens of the reactant templates, keeping
/// the maps of hydrogens that are mapped to heavy product atoms.
pub fn fix_hs(rxn: &mut ChemicalReaction) -> Result<()> {
    let mapped_to_heavy: std::collections::BTreeSet<i32> = rxn
        .products
        .iter()
        .flat_map(|m| m.atoms())
        .filter(|a| a.atomic_num() != 1)
        .filter_map(|a| int_prop(a, cp::molAtomMapNumber))
        .filter(|&map| map != 0)
        .collect();
    for mol in &mut rxn.reactants {
        for idx in 0..mol.num_atoms() {
            let atom = mol.atom(idx);
            if atom.atomic_num() != 1 {
                continue;
            }
            match int_prop(atom, cp::molAtomMapNumber) {
                Some(map) if map != 0 => {
                    if mapped_to_heavy.contains(&map) {
                        rd_log_warning!(
                            "Reaction has explicit hydrogens, reactants will need explicit hydrogens (addHs)"
                        );
                    } else {
                        mol.atom_mut(idx).clear_prop(cp::molAtomMapNumber);
                    }
                }
                _ => {}
            }
        }
    }
    for mol in &mut rxn.reactants {
        merge_query_hs(mol, true, false)?;
    }
    Ok(())
}

/// `RxnOps::adjustTemplates`
pub fn adjust_templates(templates: &mut [ROMol], params: &AdjustQueryParameters) -> Result<()> {
    templates
        .iter_mut()
        .try_for_each(|t| adjust_query_properties_in_place(t, params))
}

/// `RxnOps::sanitizeRxn`: carries out the operations selected by `ops`; on
/// failure the error comes with the flag of the operation that failed.
pub fn sanitize_rxn(
    rxn: &mut ChemicalReaction,
    ops: u32,
    params: &AdjustQueryParameters,
) -> std::result::Result<(), (u32, RDKitError)> {
    if ops & SANITIZE_RGROUP_NAMES != 0 {
        fix_rgroups(rxn).map_err(|e| (SANITIZE_RGROUP_NAMES, e))?;
    }
    if ops & SANITIZE_ATOM_MAPS != 0 {
        fix_atom_maps(rxn).map_err(|e| (SANITIZE_ATOM_MAPS, e))?;
    }
    if ops & SANITIZE_ADJUST_REACTANTS != 0 {
        adjust_templates(&mut rxn.reactants, params).map_err(|e| (SANITIZE_ADJUST_REACTANTS, e))?;
    }
    if ops & SANITIZE_ADJUST_PRODUCTS != 0 {
        adjust_templates(&mut rxn.products, params).map_err(|e| (SANITIZE_ADJUST_PRODUCTS, e))?;
    }
    if ops & SANITIZE_MERGEHS != 0 {
        fix_hs(rxn).map_err(|e| (SANITIZE_MERGEHS, e))?;
    }
    Ok(())
}

/// `RxnOps::fixReactantTemplateAromaticity`
pub fn fix_reactant_template_aromaticity(rxn: &mut ChemicalReaction) -> Result<()> {
    for mol in &mut rxn.reactants {
        sanitize_mol(mol, SANITIZE_SETAROMATICITY).map_err(|(_, e)| e)?;
    }
    Ok(())
}

/// `RxnOps::sanitizeRxnAsMols`: the usual molecule sanitization of every template.
pub fn sanitize_rxn_as_mols(rxn: &mut ChemicalReaction, sanitize_ops: u32) -> Result<()> {
    rxn.reactants
        .iter_mut()
        .chain(&mut rxn.agents)
        .chain(&mut rxn.products)
        .try_for_each(|mol| sanitize_mol(mol, sanitize_ops).map_err(|(_, e)| e))
}
