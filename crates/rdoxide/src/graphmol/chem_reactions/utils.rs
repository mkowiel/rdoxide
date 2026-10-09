//! Ported from RDKit `ChemReactions/ReactionUtils.{h,cpp}` and the helper
//! queries of `Reaction.cpp`: substructure queries on reactions, atom-map
//! housekeeping, stereo bookkeeping and `getReactingAtoms`.

use std::collections::BTreeMap;

use super::ChemicalReaction;
use crate::Result;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::substruct_match;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::utils::count_swaps_to_interconvert;
use crate::rdgeneral::{RDProps, RDValue};

/// The three roles a molecule can have in a reaction (`ReactionMoleculeType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionMoleculeType {
    Reactant,
    Product,
    Agent,
}

pub fn templates(rxn: &ChemicalReaction, t: ReactionMoleculeType) -> &[ROMol] {
    match t {
        ReactionMoleculeType::Reactant => rxn.reactants(),
        ReactionMoleculeType::Product => rxn.products(),
        ReactionMoleculeType::Agent => rxn.agents(),
    }
}

fn has_template_match(
    rxn: &ChemicalReaction,
    query: &ChemicalReaction,
    t: ReactionMoleculeType,
) -> Result<bool> {
    let (mine, theirs) = (templates(rxn, t), templates(query, t));
    if mine.len() < theirs.len() {
        return Ok(false);
    }
    if theirs.is_empty() {
        return Ok(true);
    }
    for mol in mine {
        for q in theirs {
            if !substruct_match(mol, q, rxn.substruct_params())?.is_empty() {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// `hasReactantTemplateSubstructMatch`
pub fn has_reactant_template_substruct_match(
    rxn: &ChemicalReaction,
    query: &ChemicalReaction,
) -> Result<bool> {
    has_template_match(rxn, query, ReactionMoleculeType::Reactant)
}

/// `hasProductTemplateSubstructMatch`
pub fn has_product_template_substruct_match(
    rxn: &ChemicalReaction,
    query: &ChemicalReaction,
) -> Result<bool> {
    has_template_match(rxn, query, ReactionMoleculeType::Product)
}

/// `hasAgentTemplateSubstructMatch`
pub fn has_agent_template_substruct_match(
    rxn: &ChemicalReaction,
    query: &ChemicalReaction,
) -> Result<bool> {
    has_template_match(rxn, query, ReactionMoleculeType::Agent)
}

/// `hasReactionSubstructMatch`
pub fn has_reaction_substruct_match(
    rxn: &ChemicalReaction,
    query: &ChemicalReaction,
    include_agents: bool,
) -> Result<bool> {
    Ok(has_reactant_template_substruct_match(rxn, query)?
        && has_product_template_substruct_match(rxn, query)?
        && (!include_agents || has_agent_template_substruct_match(rxn, query)?))
}

fn num_atoms_with_prop(mol: &ROMol, prop: &str) -> usize {
    mol.atoms().iter().filter(|a| a.has_prop(prop)).count()
}

/// `hasReactionAtomMapping`: do any reactant or product templates carry
/// atom-map numbers?
pub fn has_reaction_atom_mapping(rxn: &ChemicalReaction) -> bool {
    rxn.reactants()
        .iter()
        .chain(rxn.products())
        .any(|m| num_atoms_with_prop(m, cp::molAtomMapNumber) > 0)
}

/// `isReactionTemplateMoleculeAgent`: a template is an agent when the
/// fraction of its heavy atoms that are mapped is below `agent_threshold`.
pub fn is_reaction_template_molecule_agent(mol: &ROMol, agent_threshold: f64) -> bool {
    let mapped = num_atoms_with_prop(mol, cp::molAtomMapNumber) as f64;
    let heavy = f64::from(mol.num_heavy_atoms());
    !(heavy > 0.0 && mapped / heavy >= agent_threshold)
}

/// `removeMappingNumbersFromReactions`
pub fn remove_mapping_numbers_from_reactions(rxn: &mut ChemicalReaction) {
    let clear = |mols: &mut Vec<ROMol>| {
        for mol in mols {
            for i in 0..mol.num_atoms() {
                mol.atom_mut(i).clear_prop(cp::molAtomMapNumber);
            }
        }
    };
    clear(rxn.agents_mut());
    clear(rxn.products_mut());
    clear(rxn.reactants_mut());
}

/// The map numbers of an atom's neighbors in bond order (`-1` for unmapped
/// ones, plus one for an extra neighbor when `other` has the higher degree),
/// and the number of unmapped entries.
fn neighbor_order(mol: &ROMol, atom: u32, other_degree: u32) -> (u32, Vec<i32>) {
    let mut order = Vec::with_capacity(mol.atom_degree(atom) as usize);
    let mut unmapped = 0;
    for nbr in mol.atom_neighbors(atom) {
        match mol.atom(nbr).atom_map_num() {
            n if n > 0 => order.push(n),
            _ => {
                order.push(-1);
                unmapped += 1;
            }
        }
    }
    if mol.atom_degree(atom) < other_degree {
        order.push(-1);
        unmapped += 1;
    }
    (unmapped, order)
}

/// Can every mapped entry of `reference` be found in `order` (using one
/// unmapped slot per missing entry)?
fn check_order_overlap(order: &mut [i32], unmapped: u32, reference: &[i32]) -> bool {
    for &elem in reference.iter().filter(|&&e| e >= 0) {
        if order.contains(&elem) {
            continue;
        }
        if unmapped == 0 {
            return false;
        }
        match order.iter_mut().find(|o| **o == -1) {
            Some(slot) => *slot = elem,
            None => return false,
        }
    }
    true
}

/// `countSwapsBetweenReactantAndProduct`: the number of neighbor swaps that
/// turn the reactant atom's neighbor order into the product atom's, or `-1`
/// when they cannot be compared.
pub fn count_swaps_between_reactant_and_product(
    react_mol: &ROMol,
    react_atom: u32,
    prod_mol: &ROMol,
    prod_atom: u32,
) -> i32 {
    let (rd, pd) = (
        react_mol.atom_degree(react_atom),
        prod_mol.atom_degree(prod_atom),
    );
    if rd >= 3 && pd >= 3 && rd.abs_diff(pd) <= 1 {
        let (r_unmapped, mut r_order) = neighbor_order(react_mol, react_atom, pd);
        if r_unmapped <= 1 {
            let (p_unmapped, mut p_order) = neighbor_order(prod_mol, prod_atom, rd);
            if p_unmapped <= 1
                && check_order_overlap(&mut r_order, r_unmapped, &p_order)
                && check_order_overlap(&mut p_order, p_unmapped, &r_order)
            {
                return count_swaps_to_interconvert(&r_order, &p_order).map_or(-1, |n| n as i32);
            }
        }
    }
    -1
}

fn is_chiral(tag: ChiralType) -> bool {
    tag != ChiralType::Unspecified && tag != ChiralType::Other
}

/// `updateProductsStereochem`: sets each product atom's `molInversionFlag`
/// (1 inversion, 2 retention, 3 stereo destroyed, 4 stereo created) from the
/// chirality of the mapped reactant and product template atoms.
pub fn update_products_stereochem(rxn: &mut ChemicalReaction) -> Result<()> {
    // map number -> (reactant template, atom)
    let mut reactant_mapping: BTreeMap<i32, (usize, u32)> = BTreeMap::new();
    for (ri, mol) in rxn.reactants().iter().enumerate() {
        for (ai, atom) in mol.atoms().iter().enumerate() {
            if let Some(n) = atom
                .get_prop(cp::molAtomMapNumber)
                .and_then(RDValue::as_int)
            {
                reactant_mapping.insert(n, (ri, ai as u32));
            }
        }
    }
    let mut flags: Vec<(usize, u32, i32)> = Vec::new();
    for (pi, prod) in rxn.products().iter().enumerate() {
        for ai in 0..prod.num_atoms() {
            let atom = prod.atom(ai);
            if atom.has_prop(cp::molInversionFlag) {
                continue;
            }
            let Some(map_num) = atom
                .get_prop(cp::molAtomMapNumber)
                .and_then(RDValue::as_int)
            else {
                flags.push((pi, ai, 4));
                continue;
            };
            let Some(&(ri, rai)) = reactant_mapping.get(&map_num) else {
                flags.push((pi, ai, 4));
                continue;
            };
            let react_mol = &rxn.reactants()[ri];
            let react_atom = react_mol.atom(rai);
            if is_chiral(atom.chiral_tag()) {
                if is_chiral(react_atom.chiral_tag()) {
                    let mut flag = if react_atom.chiral_tag() == atom.chiral_tag() {
                        2
                    } else {
                        1
                    };
                    let swaps = count_swaps_between_reactant_and_product(react_mol, rai, prod, ai);
                    if swaps >= 0 && swaps % 2 == 1 {
                        flag = if flag == 1 { 2 } else { 1 };
                    }
                    flags.push((pi, ai, flag));
                } else {
                    flags.push((pi, ai, 4));
                }
            } else if is_chiral(react_atom.chiral_tag()) {
                flags.push((pi, ai, 3));
            }
        }
    }
    for (pi, ai, flag) in flags {
        rxn.products_mut()[pi]
            .atom_mut(ai)
            .dict_mut()
            .set_val(cp::molInversionFlag, flag);
    }
    Ok(())
}

/// Does `a` (a bond-query description plus value) describe the same bond
/// query as `b`? Used by [`is_changed_atom`].
fn bond_queries_unchanged(
    r: &crate::graphmol::bond::Bond,
    p: &crate::graphmol::bond::Bond,
) -> bool {
    match (r.query(), p.query()) {
        (Some(rq), Some(pq)) => {
            if pq.description() == "BondNull" {
                true
            } else if r.bond_type() == BondType::Single
                && p.bond_type() == BondType::Single
                && ((rq.description() == "BondOr" && pq.description() == "BondOr")
                    || (rq.description() == "SingleOrAromaticBond"
                        && pq.description() == "SingleOrAromaticBond"))
            {
                true
            } else {
                r.bond_type() == p.bond_type()
                    && rq.description() == "BondOrder"
                    && pq.description() == "BondOrder"
                    && rq.val == pq.val
            }
        }
        (Some(_), None) => false,
        (None, Some(pq)) => pq.description() == "BondNull",
        (None, None) => r.bond_type() == p.bond_type(),
    }
}

/// `isChangedAtom`
fn is_changed_atom(r_mol: &ROMol, r_idx: u32, p_mol: &ROMol, p_idx: u32) -> bool {
    let (r, p) = (r_mol.atom(r_idx), p_mol.atom(p_idx));
    if r.atomic_num() != p.atomic_num() && p.atomic_num() > 0 {
        return true;
    }
    if r_mol.atom_degree(r_idx) != p_mol.atom_degree(p_idx) {
        return true;
    }
    if p.atomic_num() > 0 && crate::graphmol::queries::is_complex_atom_query(r.query()) {
        return true;
    }
    let mut reactant_bonds: BTreeMap<i32, u32> = BTreeMap::new();
    for nbr in r_mol.atom_neighbors(r_idx) {
        match r_mol
            .atom(nbr)
            .get_prop(cp::molAtomMapNumber)
            .and_then(RDValue::as_int)
        {
            Some(n) => {
                if let Some(b) = r_mol.bond_idx_between(r_idx, nbr) {
                    reactant_bonds.insert(n, b);
                }
            }
            None => return true,
        }
    }
    for nbr in p_mol.atom_neighbors(p_idx) {
        if let Some(n) = p_mol
            .atom(nbr)
            .get_prop(cp::molAtomMapNumber)
            .and_then(RDValue::as_int)
        {
            let Some(&rb) = reactant_bonds.get(&n) else {
                return true;
            };
            let Some(pb) = p_mol.bond_idx_between(p_idx, nbr) else {
                return true;
            };
            if !bond_queries_unchanged(r_mol.bond(rb), p_mol.bond(pb)) {
                return true;
            }
        }
    }
    false
}

/// `getReactingAtoms`: for each reactant template the indices of the atoms
/// that change in the reaction (see upstream's documentation for the rules);
/// with `mapped_atoms_only` unmapped atoms are not automatically included.
pub fn get_reacting_atoms(
    rxn: &ChemicalReaction,
    mapped_atoms_only: bool,
) -> Result<Vec<Vec<i32>>> {
    if !rxn.is_initialized() {
        return Err(super::not_initialized());
    }
    let mut mapped_product_atoms: BTreeMap<i32, (usize, u32)> = BTreeMap::new();
    for (pi, mol) in rxn.products().iter().enumerate() {
        for (ai, atom) in mol.atoms().iter().enumerate() {
            if let Some(n) = atom
                .get_prop(cp::molAtomMapNumber)
                .and_then(RDValue::as_int)
            {
                mapped_product_atoms.insert(n, (pi, ai as u32));
            }
        }
    }
    let mut res = Vec::with_capacity(rxn.num_reactant_templates());
    for mol in rxn.reactants() {
        let mut changed = Vec::new();
        for (ai, atom) in mol.atoms().iter().enumerate() {
            match atom
                .get_prop(cp::molAtomMapNumber)
                .and_then(RDValue::as_int)
            {
                None => {
                    if !mapped_atoms_only {
                        changed.push(ai as i32);
                    }
                }
                Some(n) => {
                    if let Some(&(pi, pai)) = mapped_product_atoms.get(&n) {
                        if is_changed_atom(mol, ai as u32, &rxn.products()[pi], pai) {
                            changed.push(ai as i32);
                        }
                    }
                }
            }
        }
        res.push(changed);
    }
    Ok(res)
}
