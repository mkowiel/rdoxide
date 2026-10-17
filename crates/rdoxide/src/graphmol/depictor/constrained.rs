//! Depictions constrained to, or aligned onto, a reference structure.
//!
//! Ported from RDKit `Code/GraphMol/Depictor/RDDepictor.cpp` (the
//! `generateDepictionMatching*` functions) and the R-group helpers of
//! `DepictUtils.cpp` (BSD-3-Clause).

use super::depict_utils::IntPoint2DMap;
use super::{
    Compute2DCoordParameters, compute_2d_coords, compute_2d_coords_mimic_dist_mat, depict_error,
};
use crate::Result;
use crate::geometry::{Point2D, Transform3D};
use crate::graphmol::adjust_query::{AdjustQueryParameters, adjust_query_properties_in_place};
use crate::graphmol::mol_align::{
    BestAlignmentParams, get_alignment_transform, get_best_alignment_transform,
};
use crate::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use crate::graphmol::mol_transforms::transform_conformer;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::substruct::{
    MatchVect, SubstructMatchParameters, get_most_substituted_core_match,
    is_atom_terminal_r_group_or_query_hydrogen, sort_matches_by_degree_of_core_substitution,
    substruct_match,
};
use crate::graphmol::wedge_bonds::{clear_mol_block_wedging_info, invert_mol_block_wedging_info};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use std::collections::BTreeMap;

const FORMER_NBR_INDICES: &str = "__formerNbrIndices";
const FORMER_IDX: &str = "__formerIdx";

/// `ConstrainedDepictionParams`
#[derive(Debug, Clone)]
pub struct ConstrainedDepictionParams {
    /// If false, a depict error is returned when the molecule does not match
    /// the reference; if true, an unconstrained depiction is generated.
    pub accept_failure: bool,
    pub force_rdkit: bool,
    /// Terminal dummy atoms in the reference are ignored if they match an
    /// implicit hydrogen or are attached to a query atom.
    pub allow_rgroups: bool,
    /// Rigid-body align the existing (or an unconstrained) conformer to the
    /// reference instead of rebuilding the coordinates around it.
    pub align_only: bool,
    /// Update or clear existing molblock wedging information as required.
    pub adjust_mol_block_wedging: bool,
    /// The conformer to align (`align_only`) or to check the wedging against.
    pub existing_conf_id: i32,
    pub use_ring_templates: bool,
}

impl Default for ConstrainedDepictionParams {
    fn default() -> Self {
        ConstrainedDepictionParams {
            accept_failure: false,
            force_rdkit: false,
            allow_rgroups: false,
            align_only: false,
            adjust_mol_block_wedging: true,
            existing_conf_id: -1,
            use_ring_templates: false,
        }
    }
}

fn unconstrained_coords(mol: &mut ROMol, p: &ConstrainedDepictionParams) -> Result<u32> {
    let params = Compute2DCoordParameters {
        canon_orient: false,
        clear_confs: true,
        force_rdkit: p.force_rdkit,
        use_ring_templates: p.use_ring_templates,
        ..Default::default()
    };
    compute_2d_coords(mol, &params)
}

/// `hasTerminalRGroupOrQueryHydrogen`
fn has_terminal_r_group_or_query_hydrogen(mol: &ROMol) -> bool {
    (0..mol.num_atoms()).any(|i| is_atom_terminal_r_group_or_query_hydrogen(mol, i))
}

/// `prepareTemplateForRGroups`: adjusts the query bonds of the template and
/// returns a copy without the terminal R-group atoms (remembering their
/// indices and attachment points), if there are any.
fn prepare_template_for_r_groups(template: &mut RWMol) -> Result<Option<RWMol>> {
    let mut params = AdjustQueryParameters::no_adjustments();
    params.adjust_single_bonds_to_degree_one_neighbors = true;
    params.adjust_single_bonds_between_aromatic_atoms = true;
    adjust_query_properties_in_place(template, &params)?;
    let mut removed: BTreeMap<u32, u32> = BTreeMap::new();
    for bond in template.bonds() {
        let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        if is_atom_terminal_r_group_or_query_hydrogen(template, begin)
            && template.atom(end).has_query()
        {
            removed.insert(begin, end);
        } else if is_atom_terminal_r_group_or_query_hydrogen(template, end)
            && template.atom(begin).has_query()
        {
            removed.insert(end, begin);
        }
    }
    if removed.is_empty() {
        return Ok(None);
    }
    let mut reduced = template.clone();
    for idx in 0..reduced.num_atoms() {
        reduced.atom_mut(idx).set_prop(FORMER_IDX, idx, false)?;
        if let Some(&nbr) = removed.get(&idx) {
            let atom = reduced.atom_mut(nbr);
            let mut former = match atom.get_prop(FORMER_NBR_INDICES) {
                Some(RDValue::VecUnsignedInt(v)) => v.clone(),
                _ => Vec::new(),
            };
            former.push(idx);
            atom.set_prop(FORMER_NBR_INDICES, former, false)?;
        }
    }
    reduced.begin_batch_edit()?;
    for &idx in removed.keys() {
        reduced.remove_atom(idx)?;
    }
    reduced.commit_batch_edit()?;
    Ok(Some(reduced))
}

/// `reducedToFullMatches`: maps matches of the reduced query back to the
/// full one, matching each removed R-group atom to a so far unmatched
/// neighbour.
fn reduced_to_full_matches(reduced: &RWMol, mol_hs: &RWMol, matches: &mut [MatchVect]) {
    for m in matches.iter_mut() {
        let mut matched = vec![false; mol_hs.num_atoms() as usize];
        for &(_, a) in m.iter() {
            matched[a as usize] = true;
        }
        let mut new_match: MatchVect = Vec::new();
        for pair in m.iter_mut() {
            let q_atom = reduced.atom(pair.0 as u32);
            let former_idx = q_atom
                .get_prop(FORMER_IDX)
                .and_then(RDValue::as_uint)
                .expect("reduced query atoms have a former index");
            pair.0 = former_idx as i32;
            let mut former_nbrs = match q_atom.get_prop(FORMER_NBR_INDICES) {
                Some(RDValue::VecUnsignedInt(v)) => v.clone(),
                _ => Vec::new(),
            };
            for nbr in mol_hs.atom_neighbors(pair.1 as u32) {
                let Some(&former) = former_nbrs.last() else {
                    break;
                };
                if !matched[nbr as usize] {
                    former_nbrs.pop();
                    new_match.push((former as i32, nbr as i32));
                }
            }
        }
        m.extend(new_match);
    }
}

/// `invertWedgingIfMolHasFlipped`
fn invert_wedging_if_mol_has_flipped(mol: &mut ROMol, trans: &Transform3D) -> bool {
    const FLIP_THRESHOLD: f64 = -0.99;
    let should_flip = trans.get_val(2, 2) < FLIP_THRESHOLD;
    if should_flip {
        invert_mol_block_wedging_info(mol);
    }
    should_flip
}

/// The anonymous-namespace `removeAllConformersButOne`
fn remove_all_conformers_but_one(mol: &mut ROMol, conf_id: i32) -> Result<()> {
    let ids: Vec<u32> = mol.conformers().iter().map(|c| c.id()).collect();
    for (pos, &id) in ids.iter().enumerate() {
        if (conf_id != -1 && id as i32 == conf_id) || (conf_id == -1 && pos == 0) {
            continue;
        }
        mol.remove_conformer(id);
    }
    if mol.num_conformers() != 1 {
        return Err(crate::RDKitError::invariant("num_conformers == 1"));
    }
    mol.conformer_mut(-1)?.set_id(0);
    Ok(())
}

fn match_params(uniquify: bool, max_matches: u32) -> SubstructMatchParameters {
    SubstructMatchParameters {
        uniquify,
        max_matches,
        ..Default::default()
    }
}

fn first_match(mol: &ROMol, query: &ROMol) -> Result<MatchVect> {
    let mut ms = substruct_match(mol, query, &match_params(true, 1))?;
    Ok(if ms.is_empty() {
        Vec::new()
    } else {
        ms.swap_remove(0)
    })
}

/// `generateDepictionMatching2DStructure(mol, reference, refMatchVect, ...)`:
/// 2D coordinates where the atoms of `ref_match_vect` (reference index,
/// molecule index) are hard-constrained to, or aligned onto, the reference.
pub fn generate_depiction_matching_2d_structure_with_map(
    mol: &mut ROMol,
    reference: &ROMol,
    ref_match_vect: &MatchVect,
    conf_id: i32,
    p: &ConstrainedDepictionParams,
) -> Result<()> {
    if ref_match_vect.len() > reference.num_atoms() as usize {
        return Err(depict_error(
            "When a refMatchVect is provided, it must have size <= number of atoms in the reference",
        ));
    }
    for &(r, m) in ref_match_vect {
        if r > reference.num_atoms() as i32 {
            return Err(depict_error(
                "Reference atom index in refMatchVect out of range",
            ));
        }
        if m > mol.num_atoms() as i32 {
            return Err(depict_error(
                "Molecule atom index in refMatchVect out of range",
            ));
        }
    }
    let has_existing_coords = mol.num_conformers() > 0;
    let mut should_clear_wedging = p.adjust_mol_block_wedging && !has_existing_coords;
    let mut should_invert_wedging = false;
    let mut trans = Transform3D::identity();
    if p.align_only {
        if !has_existing_coords {
            unconstrained_coords(mol, p)?;
        }
        let atom_map: MatchVect = ref_match_vect.iter().map(|&(r, m)| (m, r)).collect();
        let (_, t) = get_alignment_transform(
            mol,
            reference,
            p.existing_conf_id,
            conf_id,
            Some(&atom_map),
            None,
            false,
            50,
        )?;
        trans = t;
        transform_conformer(mol.conformer_mut(p.existing_conf_id)?, &trans);
        remove_all_conformers_but_one(mol, p.existing_conf_id)?;
        if !should_clear_wedging {
            should_invert_wedging = p.adjust_mol_block_wedging;
        }
    } else {
        let ref_conf = reference.conformer(conf_id)?;
        let mut coord_map = IntPoint2DMap::new();
        for &(r, m) in ref_match_vect {
            let pt = ref_conf.atom_pos(r as u32)?;
            coord_map.insert(m as u32, Point2D::new(pt.x, pt.y));
        }
        let params = Compute2DCoordParameters {
            coord_map: Some(&coord_map),
            canon_orient: false,
            clear_confs: !(p.adjust_mol_block_wedging && has_existing_coords),
            force_rdkit: p.force_rdkit,
            use_ring_templates: p.use_ring_templates,
            ..Default::default()
        };
        let new_conf_id = compute_2d_coords(mol, &params)?;
        if p.adjust_mol_block_wedging {
            // the existing wedging information has to be cleared if the
            // original molecule had no coordinates, or if wedges are
            // outside the constrained scaffold
            const RMSD_THRESHOLD: f64 = 1.0e-2;
            const MSD_THRESHOLD: f64 = RMSD_THRESHOLD * RMSD_THRESHOLD;
            if !should_clear_wedging {
                let mut in_scaffold = vec![false; mol.num_atoms() as usize];
                for &(_, m) in ref_match_vect {
                    in_scaffold[m as usize] = true;
                }
                should_clear_wedging = mol.bonds().iter().any(|b| {
                    (b.has_prop(cp::_MolFileBondStereo) || b.has_prop(cp::_MolFileBondCfg))
                        && (!in_scaffold[b.begin_atom_idx() as usize]
                            || !in_scaffold[b.end_atom_idx() as usize])
                });
            }
            if !should_clear_wedging {
                // the scaffold coordinates may have changed (CoordGen)
                let mol_pos = mol.conformer(new_conf_id as i32)?.positions();
                let ref_pos = ref_conf.positions();
                should_clear_wedging = ref_match_vect.iter().any(|&(r, m)| {
                    (mol_pos[m as usize] - ref_pos[r as usize]).length_sq() > MSD_THRESHOLD
                });
            }
            if !should_clear_wedging {
                // the molecule may still have flipped to match the scaffold
                let identity: MatchVect = ref_match_vect.iter().map(|&(_, m)| (m, m)).collect();
                let (rmsd, t) = get_alignment_transform(
                    mol,
                    mol,
                    new_conf_id as i32,
                    p.existing_conf_id,
                    Some(&identity),
                    None,
                    false,
                    50,
                )?;
                trans = t;
                if rmsd > RMSD_THRESHOLD {
                    should_clear_wedging = true;
                } else {
                    should_invert_wedging = true;
                }
            }
        }
        if has_existing_coords {
            remove_all_conformers_but_one(mol, new_conf_id as i32)?;
        }
    }
    if should_clear_wedging {
        clear_mol_block_wedging_info(mol);
    } else if should_invert_wedging {
        invert_wedging_if_mol_has_flipped(mol, &trans);
    }
    Ok(())
}

/// `generateDepictionMatching2DStructure(mol, reference, confId,
/// referencePattern, params)`: finds the substructure match between the
/// molecule and the reference (or `reference_pattern`) itself; returns the
/// (reference index, molecule index) pairs used.
pub fn generate_depiction_matching_2d_structure(
    mol: &mut ROMol,
    reference: &ROMol,
    conf_id: i32,
    reference_pattern: Option<&ROMol>,
    params: &ConstrainedDepictionParams,
) -> Result<MatchVect> {
    let mut p = params.clone();
    let query: &ROMol = reference_pattern.unwrap_or(reference);
    // the allowRGroups logic is not needed without terminal dummy atoms
    p.allow_rgroups = p.allow_rgroups && has_terminal_r_group_or_query_hydrogen(query);
    let mut mol_hs: Option<RWMol> = None;
    let mut query_adj: Option<RWMol> = None;
    let mut reduced_query: Option<RWMol> = None;
    if p.allow_rgroups {
        let mut hs = mol.clone();
        add_hs(&mut hs, &AddHsParameters::default(), None)?;
        mol_hs = Some(hs);
        let mut adj = query.clone();
        reduced_query = prepare_template_for_r_groups(&mut adj)?;
        query_adj = Some(adj);
    }
    let prb_mol: &ROMol = mol_hs.as_ref().unwrap_or(mol);
    let ref_mol: &ROMol = match (&reduced_query, &query_adj) {
        (Some(r), _) => r,
        (None, Some(a)) => a,
        (None, None) => query,
    };
    // maps the atoms of the query to those of the reference
    let mut pattern_to_ref: Vec<i32> = vec![-1; query.num_atoms() as usize];
    if let Some(pattern) = reference_pattern {
        let best: MatchVect = if p.allow_rgroups {
            let mut reference_hs = reference.clone();
            add_hs(&mut reference_hs, &AddHsParameters::default(), None)?;
            let adj = query_adj.as_ref().expect("set with allow_rgroups");
            let mut matches = substruct_match(&reference_hs, ref_mol, &match_params(true, 1000))?;
            if let Some(reduced) = &reduced_query {
                reduced_to_full_matches(reduced, &reference_hs, &mut matches);
            }
            if matches.is_empty() {
                Vec::new()
            } else {
                get_most_substituted_core_match(&reference_hs, adj, &matches).clone()
            }
        } else {
            first_match(reference, pattern)?
        };
        if best.is_empty() {
            return Err(depict_error("Reference pattern does not map to reference."));
        }
        let num_ref_atoms = reference.num_atoms() as i32;
        for &(q, r) in &best {
            // skip the indices of added Hs
            if p.allow_rgroups && r >= num_ref_atoms {
                continue;
            }
            pattern_to_ref[q as usize] = r;
        }
    } else {
        for (i, v) in pattern_to_ref.iter_mut().enumerate() {
            *v = i as i32;
        }
    }
    let mut match_vect: MatchVect = Vec::new();
    if p.align_only {
        let mut matches = substruct_match(prb_mol, ref_mol, &match_params(false, 1000))?;
        if !matches.is_empty() {
            if p.allow_rgroups {
                let adj = query_adj.as_ref().expect("set with allow_rgroups");
                // match the maximum number of R-groups to heavy atoms
                if let Some(reduced) = &reduced_query {
                    reduced_to_full_matches(reduced, prb_mol, &mut matches);
                }
                matches = sort_matches_by_degree_of_core_substitution(prb_mol, adj, &matches);
                let mut max_matched_heavies = -1;
                let mut max_pruned_size = -1;
                let mut pruned_matches: Vec<MatchVect> = Vec::new();
                let num_mol_atoms = mol.num_atoms() as i32;
                for m in &matches {
                    // added hydrogens were not originally part of the molecule
                    let mut n_heavies = 0;
                    let mut pruned: MatchVect = Vec::new();
                    for &(q, a) in m {
                        if is_atom_terminal_r_group_or_query_hydrogen(adj, q as u32) {
                            if a >= num_mol_atoms {
                                continue;
                            }
                            n_heavies += 1;
                        }
                        if pattern_to_ref[q as usize] == -1 {
                            continue;
                        }
                        pruned.push((q, a));
                    }
                    if n_heavies < max_matched_heavies {
                        break;
                    }
                    max_matched_heavies = n_heavies;
                    let size = pruned.len() as i32;
                    if size > max_pruned_size {
                        max_pruned_size = size;
                        pruned_matches.clear();
                    }
                    if size == max_pruned_size {
                        pruned_matches.push(pruned);
                    }
                }
                matches = pruned_matches;
            }
            // the alignment wants (molecule index, reference index)
            for m in matches.iter_mut() {
                for pair in m.iter_mut() {
                    let ref_idx = pattern_to_ref[pair.0 as usize];
                    *pair = (pair.1, ref_idx);
                }
            }
            // a molecule without coordinates needs some before the alignment
            if mol.num_conformers() == 0 {
                unconstrained_coords(mol, &p)?;
                if p.adjust_mol_block_wedging {
                    clear_mol_block_wedging_info(mol);
                    p.adjust_mol_block_wedging = false;
                }
            }
            // cap the effort spent on finding the best alignment
            const MAX_MATCHES: u32 = 1000;
            let align_params = BestAlignmentParams {
                max_matches: MAX_MATCHES,
                map: matches,
                ..Default::default()
            };
            let (_, trans, best) = get_best_alignment_transform(
                mol,
                reference,
                &align_params,
                p.existing_conf_id,
                conf_id,
                false,
                50,
            )?;
            match_vect = best.into_iter().map(|(a, b)| (b, a)).collect();
            transform_conformer(mol.conformer_mut(p.existing_conf_id)?, &trans);
            remove_all_conformers_but_one(mol, p.existing_conf_id)?;
            if p.adjust_mol_block_wedging {
                invert_wedging_if_mol_has_flipped(mol, &trans);
            }
        }
    } else {
        // a full rebuild of the coordinates around the constrained reference
        if p.allow_rgroups {
            let mut matches = substruct_match(prb_mol, ref_mol, &match_params(false, 1000))?;
            if !matches.is_empty() {
                let adj = query_adj.as_ref().expect("set with allow_rgroups");
                if let Some(reduced) = &reduced_query {
                    reduced_to_full_matches(reduced, prb_mol, &mut matches);
                }
                let num_mol_atoms = mol.num_atoms() as i32;
                for &(q, a) in get_most_substituted_core_match(prb_mol, adj, &matches) {
                    if a < num_mol_atoms && pattern_to_ref[q as usize] != -1 {
                        match_vect.push((q, a));
                    }
                }
            }
        } else {
            match_vect = first_match(prb_mol, ref_mol)?;
        }
        if !match_vect.is_empty() {
            for pair in match_vect.iter_mut() {
                pair.0 = pattern_to_ref[pair.0 as usize];
            }
            generate_depiction_matching_2d_structure_with_map(
                mol,
                reference,
                &match_vect,
                conf_id,
                &p,
            )?;
        }
    }
    if match_vect.is_empty() {
        if !p.accept_failure {
            return Err(depict_error("Substructure match with reference not found."));
        }
        // a standard depiction, without the old wedging information
        unconstrained_coords(mol, &p)?;
        if p.adjust_mol_block_wedging {
            clear_mol_block_wedging_info(mol);
        }
    }
    Ok(match_vect)
}

/// `generateDepictionMatching3DStructure`: a depiction whose inter-atom
/// distances mimic those of a 3D reference (restricted to the atoms of
/// `reference_pattern`, if given).
pub fn generate_depiction_matching_3d_structure(
    mol: &mut ROMol,
    reference: &ROMol,
    conf_id: i32,
    reference_pattern: Option<&ROMol>,
    accept_failure: bool,
    _force_rdkit: bool,
) -> Result<()> {
    let default_depiction =
        |mol: &mut ROMol| compute_2d_coords(mol, &Compute2DCoordParameters::default()).map(|_| ());
    let num_ats = mol.num_atoms() as usize;
    if reference_pattern.is_none() && (reference.num_atoms() as usize) < num_ats {
        if accept_failure {
            return default_depiction(mol);
        }
        return Err(depict_error(
            "Reference molecule not compatible with target molecule.",
        ));
    }
    let mut mol_to_ref: Vec<i32> = vec![-1; num_ats];
    match reference_pattern.filter(|r| r.num_atoms() > 0) {
        Some(pattern) => {
            let mol_match = first_match(mol, pattern)?;
            let ref_match = first_match(reference, pattern)?;
            if mol_match.is_empty() || ref_match.is_empty() {
                if accept_failure {
                    return default_depiction(mol);
                }
                return Err(depict_error(
                    "Reference pattern didn't match molecule or reference.",
                ));
            }
            for (m, r) in mol_match.iter().zip(&ref_match) {
                mol_to_ref[m.1 as usize] = r.1;
            }
        }
        None => {
            for (i, v) in mol_to_ref.iter_mut().enumerate() {
                *v = i as i32;
            }
        }
    }
    let conf = reference.conformer(conf_id)?;
    // negative distances are ignored, so a pattern selects a subset
    let mut dmat = vec![-1.0; num_ats * num_ats.saturating_sub(1) / 2];
    for i in 0..num_ats {
        if mol_to_ref[i] == -1 {
            continue;
        }
        let cds_i = conf.atom_pos(i as u32)?;
        for j in i + 1..num_ats {
            if mol_to_ref[j] == -1 {
                continue;
            }
            let cds_j = conf.atom_pos(mol_to_ref[j] as u32)?;
            dmat[j * (j - 1) / 2 + i] = (cds_i - cds_j).length();
        }
    }
    compute_2d_coords_mimic_dist_mat(mol, &dmat, false, true, 0.5, 3, 100, 25, true)?;
    Ok(())
}
