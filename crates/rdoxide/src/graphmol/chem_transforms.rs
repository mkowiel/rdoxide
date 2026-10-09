//! Substructure-based molecule transforms (upstream
//! `GraphMol/ChemTransforms/ChemTransforms.cpp`): currently
//! `deleteSubstructs` and `combineMols`.
//!
//! Ported from RDKit 2026.09.1.
//!
//! Copyright (C) 2006-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeSet;

use crate::geometry::Point3D;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_ops::frags::get_mol_frags_atoms;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::{Result, rd_log_warning};

/// `details::updateSubMolConfs`: rebuilds the conformers of `res` from those
/// of `mol`, leaving out the removed atoms.
fn update_sub_mol_confs(mol: &ROMol, res: &mut RWMol, removed: &[bool]) -> Result<()> {
    res.clear_conformers();
    for conf in mol.conformers() {
        let mut new_conf = Conformer::new(res.num_atoms() as usize);
        new_conf.set_id(conf.id());
        new_conf.set_3d(conf.is_3d());
        let kept = conf
            .positions()
            .iter()
            .zip(removed)
            .filter(|(_, removed)| !**removed)
            .map(|(p, _)| *p);
        for (i, p) in kept.enumerate() {
            new_conf.set_atom_pos(i as u32, p)?;
        }
        res.add_conformer(new_conf, false)?;
    }
    Ok(())
}

/// `deleteSubstructs(mol, query, onlyFrags, useChirality)`: a copy of `mol`
/// without the atoms matched by `query`. With `only_frags` an atom is removed
/// only when a whole fragment is matched.
pub fn delete_substructs(
    mol: &ROMol,
    query: &ROMol,
    only_frags: bool,
    use_chirality: bool,
) -> Result<RWMol> {
    let mut res = mol.clone();
    let params = SubstructMatchParameters {
        use_chirality,
        ..SubstructMatchParameters::default()
    };
    let matches = substruct_match(&res, query, &params)?;
    if matches.is_empty() {
        return Ok(res);
    }
    let matches: Vec<Vec<i32>> = matches
        .into_iter()
        .map(|m| m.into_iter().map(|(_, atom)| atom).collect())
        .collect();

    let mut to_delete: BTreeSet<i32> = BTreeSet::new();
    if only_frags {
        for mut frag in get_mol_frags_atoms(&res) {
            frag.sort_unstable();
            let whole_frag_matched = matches.iter().find(|m| {
                let mut sorted = (*m).clone();
                sorted.sort_unstable();
                sorted == frag
            });
            if let Some(m) = whole_frag_matched {
                to_delete.extend(m);
            }
        }
    } else {
        for m in &matches {
            to_delete.extend(m);
        }
    }
    if to_delete.is_empty() {
        return Ok(res);
    }

    res.begin_batch_edit()?;
    let mut removed = vec![false; mol.num_atoms() as usize];
    for &idx in &to_delete {
        removed[idx as usize] = true;
        res.remove_atom(idx as u32)?;
    }
    res.commit_batch_edit()?;
    update_sub_mol_confs(mol, &mut res, &removed)?;
    res.clear_computed_props(true);
    // update our properties, but allow unhappiness:
    res.update_property_cache(false)?;
    Ok(res)
}

/// `combineMols(mol1, mol2, offset)`: the atoms of both molecules in one
/// molecule. The coordinates of `mol2` are shifted by `offset` in every
/// conformer that both molecules have.
pub fn combine_mols(mol1: &ROMol, mol2: &ROMol, offset: Point3D) -> Result<RWMol> {
    let mut res = mol1.clone();
    let n_atoms1 = res.num_atoms();
    res.insert_mol(mol2)?;

    if mol1.num_conformers() > 0 && mol2.num_conformers() > 0 {
        if mol1.num_conformers() != mol2.num_conformers() {
            rd_log_warning!("combineMols: molecules have unequal numbers of conformers");
        }
        for conf1 in res.conformers_mut() {
            match mol2.conformer(conf1.id() as i32) {
                Ok(conf2) => {
                    for (i, &p) in conf2.positions().iter().enumerate() {
                        conf1.set_atom_pos(i as u32 + n_atoms1, p + offset)?;
                    }
                }
                Err(_) => {
                    rd_log_warning!("combineMols: conformer id {} not found in mol2", conf1.id())
                }
            }
        }
    }
    res.clear_computed_props(true);
    Ok(res)
}
