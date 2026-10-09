//! Operations on enhanced stereo: annotating atoms and bonds with stereo
//! labels, simplifying a molecule whose only stereo group covers all of its
//! stereocenters, finding meso centers and turning the molfile chiral flag
//! into stereo groups.
//!
//! Ported from RDKit `Code/GraphMol/Chirality.cpp` (`addStereoAnnotations`,
//! `simplifyEnhancedStereo`, `findMesoCenters`) and
//! `Code/GraphMol/FileParsers/MolFileStereochem.cpp`
//! (`translateChiralFlagToStereoGroups`), 2026.09.1.
//!
//! Copyright (C) 2020-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::HashSet;

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::BondStereo;
use crate::graphmol::new_canon::{RankOptions, rank_mol_atoms};
use crate::graphmol::romol::ROMol;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType, assign_stereo_group_ids};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{Result, rd_log_warning};

/// The labels `add_stereo_annotations` writes. `{id}` stands for the index of
/// the stereo group and `{cip}` for the CIP code; an empty label disables the
/// corresponding annotation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoAnnotationLabels {
    /// atoms in an ABS stereo group
    pub abs: String,
    /// atoms in an OR stereo group
    pub or: String,
    /// atoms in an AND stereo group
    pub and: String,
    /// chiral atoms that are in no stereo group
    pub cip: String,
    /// double bonds with CIP stereo
    pub bond: String,
}

impl Default for StereoAnnotationLabels {
    fn default() -> Self {
        StereoAnnotationLabels {
            abs: "abs ({cip})".into(),
            or: "or{id}".into(),
            and: "and{id}".into(),
            cip: "({cip})".into(),
            bond: "({cip})".into(),
        }
    }
}

fn string_prop(props: &impl RDProps, key: &str) -> Option<String> {
    props.get_prop(key).map(RDValue::to_rd_string)
}

/// `Chirality::addStereoAnnotations`: sets the `atomNote` and `bondNote`
/// properties from the stereo groups and CIP codes of `mol`.
pub fn add_stereo_annotations(mol: &mut ROMol, labels: &StereoAnnotationLabels) {
    let mut groups = mol.stereo_groups().to_vec();
    assign_stereo_group_ids(&mut groups);
    let mut done = vec![false; mol.num_atoms() as usize];
    for group in &groups {
        let id = group.write_id().to_string();
        let label = match group.group_type() {
            StereoGroupType::Absolute => &labels.abs,
            StereoGroupType::Or => &labels.or,
            StereoGroupType::And => &labels.and,
        };
        for &atom in group.atoms() {
            if done[atom as usize] {
                rd_log_warning!(
                    "Warning: atom {atom} is in more than one stereogroup. Only the label from the first group will be used."
                );
                continue;
            }
            done[atom as usize] = true;
            if label.is_empty() {
                continue;
            }
            let mut note = label.replace("{id}", &id);
            if let Some(cip) = string_prop(mol.atom(atom), cp::_CIPCode).filter(|c| !c.is_empty()) {
                note = note.replace("{cip}", &cip);
            }
            mol.atom_mut(atom).props.set_val(cp::atomNote, note);
        }
    }
    if !labels.cip.is_empty() {
        for idx in 0..mol.num_atoms() {
            if done[idx as usize] {
                continue;
            }
            if let Some(cip) = string_prop(mol.atom(idx), cp::_CIPCode) {
                let note = labels.cip.replace("{cip}", &cip);
                mol.atom_mut(idx).props.set_val(cp::atomNote, note);
            }
        }
    }
    if !labels.bond.is_empty() {
        for idx in 0..mol.num_bonds() {
            let bond = mol.bond(idx);
            let cip = string_prop(bond, cp::_CIPCode).unwrap_or_else(|| match bond.stereo() {
                BondStereo::StereoE => "E".into(),
                BondStereo::StereoZ => "Z".into(),
                _ => String::new(),
            });
            if !cip.is_empty() {
                let note = labels.bond.replace("{cip}", &cip);
                mol.bond_mut(idx).props.set_val(cp::bondNote, note);
            }
        }
    }
}

fn is_specified_tetrahedral(tag: ChiralType) -> bool {
    tag > ChiralType::Unspecified && tag < ChiralType::Other
}

/// `Chirality::simplifyEnhancedStereo`: when the only stereo group holds all
/// specified stereocenters and is an AND or OR group, the molecule gets a
/// `molNote` ("AND enantiomer"/"OR enantiomer") and, with
/// `remove_affected_stereo_groups`, loses the group; the CIP codes of the
/// group's atoms are cleared.
pub fn simplify_enhanced_stereo(mol: &mut ROMol, remove_affected_stereo_groups: bool) {
    let [group] = mol.stereo_groups() else {
        return;
    };
    let group_type = group.group_type();
    let group_atoms = group.atoms().to_vec();
    let all_covered = (0..mol.num_atoms())
        .filter(|&a| is_specified_tetrahedral(mol.atom(a).chiral_tag()))
        .all(|a| group_atoms.contains(&a));
    if !all_covered {
        return;
    }
    let label = match group_type {
        StereoGroupType::Or => "OR enantiomer",
        StereoGroupType::And => "AND enantiomer",
        StereoGroupType::Absolute => return,
    };
    if remove_affected_stereo_groups {
        // clearing cannot fail: it never creates a second absolute group
        let _ = mol.set_stereo_groups(Vec::new());
    }
    mol.props.set_val(cp::molNote, label);
    for atom in group_atoms {
        mol.atom_mut(atom).clear_prop(cp::_CIPCode);
    }
}

/// `Chirality::findMesoCenters`: pairs of stereocenters that make `mol` a
/// meso compound. The partners are also recorded in the `_mesoOtherAtom`
/// atom property.
pub fn find_meso_centers(
    mol: &mut ROMol,
    include_isotopes: bool,
    include_atom_maps: bool,
) -> Result<Vec<(u32, u32)>> {
    let n = mol.num_atoms();
    let mut specified = vec![false; n as usize];
    let mut ring_stereo = vec![n; n as usize];
    for idx in 0..n {
        let atom = mol.atom_mut(idx);
        atom.clear_prop(cp::_mesoOtherAtom);
        specified[idx as usize] = atom.chiral_tag() > ChiralType::Unspecified;
        if let Some(other) = atom
            .props
            .get(cp::_ringStereoOtherAtom)
            .and_then(RDValue::as_int)
            .filter(|&o| o >= 0)
        {
            ring_stereo[idx as usize] = other as u32;
        }
    }
    let mut res = Vec::new();
    if !specified.contains(&true) {
        return Ok(res);
    }
    // ranks with chirality tell the centers apart; ranks that only know
    // whether chirality is present say which centers are otherwise alike
    let mut ranks = |include_chirality| {
        rank_mol_atoms(
            mol,
            RankOptions {
                break_ties: false,
                include_chirality,
                include_isotopes,
                include_atom_maps,
                include_chiral_presence: true,
                include_stereo_groups: false,
                use_non_stereo_ranks: false,
                ..RankOptions::default()
            },
        )
    };
    let chiral_ranks = ranks(true)?;
    let presence_ranks = ranks(false)?;
    for i in 0..n {
        if !specified[i as usize] {
            continue;
        }
        for j in (i + 1)..n {
            if !specified[j as usize] {
                continue;
            }
            let (iu, ju) = (i as usize, j as usize);
            if chiral_ranks[iu] != chiral_ranks[ju] && presence_ranks[iu] == presence_ranks[ju] {
                res.push((i, j));
            } else if ring_stereo[iu] == j && ring_stereo[ju] == i {
                // both atoms are part of ring stereo: they can have different
                // ranks and still be meso centers, as in N[C@H]1CC[C@@H](O)CC1
                let neighbor_ranks: HashSet<u32> = mol
                    .atom_neighbors(i)
                    .map(|nbr| presence_ranks[nbr as usize])
                    .collect();
                if (neighbor_ranks.len() as u32) < mol.atom_degree(i) {
                    res.push((i, j));
                }
            }
        }
    }
    for &(i, j) in &res {
        mol.atom_mut(i).props.set_val(cp::_mesoOtherAtom, j);
        mol.atom_mut(j).props.set_val(cp::_mesoOtherAtom, i);
    }
    Ok(res)
}

/// `translateChiralFlagToStereoGroups`: uses the molfile chiral flag to put
/// the stereocenters that are in no group yet (and the atropisomer bonds)
/// into an ABS group (flag 1) or into a group of type `zero_flag_group_type`
/// (flag 0). Nothing happens without the `_MolFileChiralFlag` property.
pub fn translate_chiral_flag_to_stereo_groups(
    mol: &mut ROMol,
    zero_flag_group_type: StereoGroupType,
) -> Result<()> {
    let Some(flag) = mol
        .props
        .get(cp::_MolFileChiralFlag)
        .and_then(RDValue::as_int)
    else {
        return Ok(());
    };
    mol.props.clear_val(cp::_MolFileChiralFlag);
    let group_type = if flag != 0 {
        StereoGroupType::Absolute
    } else {
        zero_flag_group_type
    };
    let mut groups = mol.stereo_groups().to_vec();
    let grouped_atoms: HashSet<u32> = groups
        .iter()
        .flat_map(|g| g.atoms().iter().copied())
        .collect();
    let grouped_bonds: HashSet<u32> = groups
        .iter()
        .flat_map(|g| g.bonds().iter().copied())
        .collect();
    // an existing ABS group takes the new members
    let abs_group = (group_type == StereoGroupType::Absolute)
        .then(|| {
            groups
                .iter()
                .position(|g| g.group_type() == StereoGroupType::Absolute)
        })
        .flatten();
    let atoms: Vec<u32> = (0..mol.num_atoms())
        .filter(|a| {
            !grouped_atoms.contains(a)
                && matches!(
                    mol.atom(*a).chiral_tag(),
                    ChiralType::TetrahedralCCW | ChiralType::TetrahedralCW
                )
        })
        .collect();
    let bonds: Vec<u32> = (0..mol.num_bonds())
        .filter(|b| {
            !grouped_bonds.contains(b)
                && matches!(
                    mol.bond(*b).stereo(),
                    BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW
                )
        })
        .collect();
    if atoms.is_empty() && bonds.is_empty() {
        return Ok(());
    }
    match abs_group {
        None => groups.push(StereoGroup::new(group_type, atoms, bonds, 0)?),
        Some(pos) => {
            let mut all_atoms = groups[pos].atoms().to_vec();
            all_atoms.extend(atoms);
            let mut all_bonds = groups[pos].bonds().to_vec();
            all_bonds.extend(bonds);
            groups[pos] = StereoGroup::new(StereoGroupType::Absolute, all_atoms, all_bonds, 0)?;
        }
    }
    mol.set_stereo_groups(groups)
}
