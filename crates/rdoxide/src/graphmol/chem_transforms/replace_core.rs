//! `replaceCore`: removes the core matched by a query and labels the
//! sidechains' attachment points with dummy atoms.
//!
//! Ported from RDKit 2026.09.1 `GraphMol/ChemTransforms/ChemTransforms.cpp`.
//!
//! Copyright (C) 2006-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;

use super::{copy_stereo_groups, first_match, update_sub_mol_confs};
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::BondStereo;
use crate::graphmol::chirality::MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO;
use crate::graphmol::mol_ops::add_hs::set_terminal_atom_coords;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::substruct::MatchVect;
use crate::rdgeneral::RDProps;
use crate::{RDKitError, Result};

/// Marks the bonds that were added to connect a sidechain to a new dummy.
const REPLACE_CORE_DUMMY_BOND: &str = "_replaceCoreDummyBond";

/// One `(molecule atom, core atom)` pair of the match, or a molecule atom
/// outside the match.
#[derive(Clone, Copy)]
struct SideChainMapping {
    mol_index: i32,
    core_index: i32,
    use_match: bool,
}

/// The neighbor of `bond_atom` through a bond other than `bond` that is
/// either one of the new dummy bonds or an unaffected bond of the ring `bring`.
fn find_nbr_bond(
    mol: &ROMol,
    bond: u32,
    bond_atom: u32,
    bring: &[i32],
    removed: &[bool],
) -> Option<u32> {
    mol.atom_bonds(bond_atom).find_map(|nbr_bond| {
        let other = nbr_bond.other_atom_idx(bond_atom).ok()?;
        let in_unaffected_ring_bond = !removed.get(other as usize).copied().unwrap_or(false)
            && bring.contains(&(nbr_bond.idx() as i32));
        (nbr_bond.idx() != bond
            && (nbr_bond.has_prop(REPLACE_CORE_DUMMY_BOND) || in_unaffected_ring_bond))
            .then_some(other)
    })
}

/// Sets cis/trans on the small-ring double bonds whose ring lost atoms, so
/// that the stereo of the remaining ring pieces is not lost.
fn set_sub_mol_broken_ring_stereo(mol: &mut ROMol, removed: &[bool]) -> Result<()> {
    if !mol.ring_info().is_initialized() {
        return Err(RDKitError::precondition("bad ringinfo"));
    }
    let is_removed = |a: u32| removed.get(a as usize).copied().unwrap_or(false);
    let bond_rings: Vec<Vec<i32>> = mol.ring_info().bond_rings().to_vec();
    for bring in &bond_rings {
        // is this bond ring affected by the removal?
        let affected = bring.iter().any(|&idx| {
            let bond = mol.bond(idx as u32);
            is_removed(bond.begin_atom_idx()) || is_removed(bond.end_atom_idx())
        });
        if !affected {
            continue;
        }
        for &bidx in bring {
            let bidx = bidx as u32;
            let bond = mol.bond(bidx);
            // a bond we can reasonably set cis/trans for, and where neither
            // atom is being removed?
            let candidate = (bond.is_aromatic()
                || bond.bond_type() == crate::graphmol::bond::BondType::Double)
                && bond.stereo() == BondStereo::StereoNone
                && mol.ring_info().min_bond_ring_size(bidx)? < MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO
                && !is_removed(bond.begin_atom_idx())
                && !is_removed(bond.end_atom_idx());
            if !candidate {
                continue;
            }
            // the two neighboring bonds in the ring or to the new dummy
            // atoms, not affected by the atom removal
            let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
            if let Some(begin_nbr) = find_nbr_bond(mol, bidx, begin, bring, removed)
                && let Some(end_nbr) = find_nbr_bond(mol, bidx, end, bring, removed)
            {
                mol.set_bond_stereo_atoms(bidx, begin_nbr, end_nbr)?;
                mol.bond_mut(bidx).set_stereo(BondStereo::StereoCis)?;
            }
        }
    }
    Ok(())
}

/// Validates `match_v` against both molecules and records the core atom each
/// molecule atom was matched to. Returns the mapping triples and, per
/// molecule atom, the core atom it is *used* for (`replaceDummies`) and any
/// core atom it matched.
struct Mappings {
    matches: Vec<(i32, SideChainMapping)>,
    matching_indices: Vec<i32>,
    all_indices: Vec<i32>,
    mapped: Vec<bool>,
    multiple_mapped: Vec<bool>,
}

fn collect_mappings(
    mol: &ROMol,
    core: &ROMol,
    match_v: &MatchVect,
    replace_dummies: bool,
) -> Result<Mappings> {
    let orig_num_atoms = mol.num_atoms() as usize;
    let mut m = Mappings {
        matches: Vec::with_capacity(orig_num_atoms),
        matching_indices: vec![-1; orig_num_atoms],
        all_indices: vec![-1; orig_num_atoms],
        mapped: vec![false; orig_num_atoms],
        multiple_mapped: vec![false; orig_num_atoms],
    };
    for &(core_idx, mol_idx) in match_v {
        if core_idx < 0 || core_idx >= core.num_atoms() as i32 {
            return Err(RDKitError::Value(
                "Supplied MatchVect indices out of bounds of the core molecule".into(),
            ));
        }
        if mol_idx < 0 || mol_idx >= mol.num_atoms() as i32 {
            return Err(RDKitError::Value(
                "Supplied MatchVect indices out of bounds of the target molecule".into(),
            ));
        }
        let use_match = replace_dummies || core.atom(core_idx as u32).atomic_num() > 0;
        if use_match {
            m.matching_indices[mol_idx as usize] = core_idx;
        }
        m.all_indices[mol_idx as usize] = core_idx;
        m.matches.push((
            mol_idx,
            SideChainMapping {
                mol_index: mol_idx,
                core_index: core_idx,
                use_match,
            },
        ));
        if m.mapped[mol_idx as usize] {
            m.multiple_mapped[mol_idx as usize] = true;
        }
        m.mapped[mol_idx as usize] = true;
    }
    Ok(m)
}

/// The bonds that connect an atom matched by several core atoms to the rest
/// of the match, i.e. the ones owned by the extra (terminal) core atoms.
fn multiple_owned_bonds(
    mol: &ROMol,
    core: &ROMol,
    match_v: &MatchVect,
    m: &Mappings,
) -> Result<Vec<bool>> {
    let mut owned = vec![false; mol.num_bonds() as usize];
    if !m.multiple_mapped.iter().any(|&b| b) {
        return Ok(owned);
    }
    for (_, info) in &m.matches {
        if !m.multiple_mapped[info.mol_index as usize] {
            continue;
        }
        let core_atom = info.core_index as u32;
        if core.atom_degree(core_atom) != 1 {
            return Err(RDKitError::invariant(
                "Multiple core atoms match a mol atom, but one of the core atoms has degree > 1 ",
            ));
        }
        let core_nbr = core.atom_neighbors(core_atom).next().expect("degree is 1") as i32;
        let Some(&(_, mol_nbr)) = match_v.iter().find(|&&(c, _)| c == core_nbr) else {
            return Err(RDKitError::Value(
                "Supplied MatchVect is missing the neighbor of a multiply mapped core atom".into(),
            ));
        };
        if mol_nbr > -1 {
            let bond = mol
                .bond_idx_between(info.mol_index as u32, mol_nbr as u32)
                .ok_or_else(|| RDKitError::invariant("expected bond in molecule not found"))?;
            owned[bond as usize] = true;
        }
    }
    Ok(owned)
}

/// `replaceCore(mol, coreQuery, replaceDummies, labelByIndex,
/// requireDummyMatch, useChirality)`: matches `core_query` and removes the
/// match. `None` when the core is not found.
pub fn replace_core(
    mol: &ROMol,
    core_query: &ROMol,
    replace_dummies: bool,
    label_by_index: bool,
    require_dummy_match: bool,
    use_chirality: bool,
) -> Result<Option<RWMol>> {
    match first_match(mol, core_query, use_chirality)? {
        Some(m) if !m.is_empty() => replace_core_with_match(
            mol,
            core_query,
            &m,
            replace_dummies,
            label_by_index,
            require_dummy_match,
        ),
        _ => Ok(None),
    }
}

/// `replaceCore(mol, core, matchV, ...)`: removes the core given by the
/// `(core atom, molecule atom)` pairs `match_v`, replacing each bond to a
/// sidechain by a bond to a dummy atom labelled (as isotope) by the core
/// atom it was attached to, or by order when `label_by_index` is off.
/// `None` when `require_dummy_match` is set and a sidechain is attached to a
/// non-dummy core atom.
pub fn replace_core_with_match(
    mol: &ROMol,
    core: &ROMol,
    match_v: &MatchVect,
    replace_dummies: bool,
    label_by_index: bool,
    require_dummy_match: bool,
) -> Result<Option<RWMol>> {
    let orig_num_atoms = mol.num_atoms();
    let mut m = collect_mappings(mol, core, match_v, replace_dummies)?;
    let owned_bonds = multiple_owned_bonds(mol, core, match_v, &m)?;

    let mut new_mol = mol.clone();
    // atoms to keep, indexed like `new_mol`
    let mut keep = vec![false; orig_num_atoms as usize];
    // new dummy atom -> the molecule atom it replaces
    let mut dummy_atom_map: BTreeMap<u32, u32> = BTreeMap::new();
    // sidechain atom in the molecule -> in `new_mol` (stereo group copying)
    let mut mol_atom_map: BTreeMap<u32, u32> = BTreeMap::new();

    // go through the matches in query order, not target molecule order
    for i in 0..orig_num_atoms as usize {
        if !m.mapped[i] {
            m.matches.push((
                i as i32,
                SideChainMapping {
                    mol_index: i as i32,
                    core_index: -1,
                    use_match: false,
                },
            ));
        }
    }
    m.matches.sort_by(|p1, p2| {
        if p1.1.core_index == p2.1.core_index {
            p1.0.cmp(&p2.0)
        } else {
            p1.1.core_index.cmp(&p2.1.core_index)
        }
    });

    // (label, dummy atom)
    let mut dummies: Vec<(i32, u32)> = Vec::new();
    // (begin, end) atoms of the bonds to the new dummies
    let mut all_new_bonds: Vec<(u32, u32)> = Vec::new();
    for (_, info) in &m.matches {
        if info.use_match {
            continue;
        }
        let sidechain = info.mol_index as u32;
        // we're keeping the sidechain atoms
        keep[sidechain as usize] = true;
        mol_atom_map.insert(sidechain, sidechain);

        // the neighbors that are in the original molecule; atoms and bonds
        // are added below, so collect them now
        let nbr_list: Vec<u32> = new_mol
            .atom_neighbors(sidechain)
            .take_while(|&n| n < orig_num_atoms)
            .collect();
        let mut which_nbr = 0u32;
        let mut new_bonds = Vec::new();
        for nbr in nbr_list {
            let connecting_bond = new_mol
                .bond_idx_between(sidechain, nbr)
                .ok_or_else(|| RDKitError::invariant("expected bond in molecule not found"))?;
            let mut bond_to_core = m.matching_indices[nbr as usize] > -1;
            let core_bond = if bond_to_core
                && m.all_indices[nbr as usize] > -1
                && info.core_index > -1
            {
                core.bond_idx_between(info.core_index as u32, m.all_indices[nbr as usize] as u32)
            } else {
                None
            };
            if bond_to_core && m.multiple_mapped[sidechain as usize] && info.core_index > -1 {
                // the core has several atoms mapping onto this atom: unless
                // we matched the right core bond, the bond can be used only
                // if nobody else owns it
                if core_bond.is_none() && owned_bonds[connecting_bond as usize] {
                    bond_to_core = false;
                }
            }
            if bond_to_core {
                // we've matched an atom in the core
                if require_dummy_match
                    && core
                        .atom(m.matching_indices[nbr as usize] as u32)
                        .atomic_num()
                        != 0
                {
                    return Ok(None);
                }
                let new_at = new_mol.add_atom(Atom::new(0), false);
                keep.push(true);

                // if we were not in the matching list, still keep the
                // original indices (replaceDummies=False); without a core
                // bond the label belongs to the neighbor
                let mut mapping = info.core_index;
                if core_bond.is_none() {
                    mapping = m.all_indices[nbr as usize];
                }
                // dummies are ordered like the mappings; by index they are in
                // arbitrary order right now, so save and sort later
                dummies.push((
                    if mapping != -1 {
                        mapping
                    } else {
                        m.matching_indices[nbr as usize]
                    },
                    new_at,
                ));
                dummy_atom_map.insert(new_at, nbr);

                let mut bnd = new_mol.bond(connecting_bond).copy_detached();
                // stereo of the connecting bond cannot be preserved
                if bnd.stereo() > BondStereo::StereoAny {
                    bnd.set_stereo(BondStereo::StereoAny)?;
                }
                if bnd.begin_atom_idx() == sidechain {
                    bnd.set_end_atom_idx_unchecked(new_at);
                } else {
                    bnd.set_begin_atom_idx_unchecked(new_at);
                }
                all_new_bonds.push((bnd.begin_atom_idx(), bnd.end_atom_idx()));
                new_bonds.push(bnd);

                // are we breaking a stereo bond definition by removing one of
                // its stereo atoms? If so, use the new atom
                let other_bonds: Vec<u32> = new_mol.atom_bond_indices(sidechain).collect();
                for b in other_bonds {
                    if b == connecting_bond || new_mol.bond(b).stereo() <= BondStereo::StereoAny {
                        continue;
                    }
                    for stereo_atom in new_mol.bond_mut(b).stereo_atoms_mut() {
                        if *stereo_atom == nbr as i32 {
                            *stereo_atom = new_at as i32;
                        }
                    }
                }

                // we may be changing the bond ordering at the atom, e.g.
                // replacing the N in C[C@](Cl)(N)F gives an atom ordering of
                // C[C?](Cl)(F)[X], so the chirality flag may have to be
                // adjusted
                if matches!(
                    new_mol.atom(sidechain).chiral_tag(),
                    ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
                ) {
                    let switch_it = match new_mol.atom_degree(sidechain) {
                        // N[C@](F)(Cl)C -> F[C@@](Cl)(C)X
                        4 => which_nbr % 2 == 0,
                        // things are different in the degree three case
                        // because of the implicit H: N[C@H](F)C -> [C@H](F)(C)X
                        3 => which_nbr == 1,
                        _ => false,
                    };
                    if switch_it {
                        new_mol.atom_mut(sidechain).invert_chirality();
                    }
                    // the bond to the core will now be at the end of the
                    // neighbor list
                    which_nbr = which_nbr.wrapping_sub(1);
                }
            }
            which_nbr = which_nbr.wrapping_add(1);
        }
        // add the bonds now, after we've finished the loop over neighbors
        for mut bond in new_bonds {
            bond.set_prop(REPLACE_CORE_DUMMY_BOND, 1, false)?;
            new_mol.add_bond_object(bond)?;
        }
    }

    if label_by_index {
        // don't sort, just label by the index
        for &(label, atom) in &dummies {
            new_mol.atom_mut(atom).set_isotope(label as u32);
        }
    } else {
        // sort the mapping indices, but label from 1..N. Upstream breaks ties
        // by the dummies' addresses, i.e. in allocation order on common
        // allocators; use the creation order.
        dummies.sort();
        for (n, &(_, atom)) in dummies.iter().enumerate() {
            new_mol.atom_mut(atom).set_isotope(n as u32 + 1);
        }
    }

    let total_atoms = new_mol.num_atoms() as usize;
    let mut removed = vec![false; total_atoms];
    let mut removed_ring_atom = false;
    new_mol.begin_batch_edit()?;
    for idx in 0..total_atoms {
        if !keep[idx] {
            new_mol.remove_atom(idx as u32)?;
            removed[idx] = true;
            if !removed_ring_atom
                && mol.ring_info().is_initialized()
                && mol.ring_info().num_atom_rings(idx as u32)? > 0
            {
                removed_ring_atom = true;
            }
        }
    }
    if removed_ring_atom {
        set_sub_mol_broken_ring_stereo(&mut new_mol, &removed)?;
        for bond in 0..new_mol.num_bonds() {
            new_mol.bond_mut(bond).clear_prop(REPLACE_CORE_DUMMY_BOND);
        }
    }
    new_mol.commit_batch_edit()?;

    // atom indices after the removal
    let mut kept_so_far = 0u32;
    let final_idx: Vec<u32> = removed
        .iter()
        .map(|&r| {
            let idx = kept_so_far;
            kept_so_far += u32::from(!r);
            idx
        })
        .collect();

    update_sub_mol_confs(mol, &mut new_mol, &removed[..orig_num_atoms as usize])?;

    // update any terminal dummy atom coordinates after removing the atoms
    // not in the keep list and rebuilding the conformers
    for &(begin, end) in &all_new_bonds {
        let (begin, end) = (final_idx[begin as usize], final_idx[end as usize]);
        if new_mol.atom_degree(begin) != 1 && new_mol.atom_degree(end) != 1 {
            return Err(RDKitError::invariant("neither atom has degree one"));
        }
        if new_mol.num_conformers() > 0 {
            if new_mol.atom(end).atomic_num() == 0 && new_mol.atom_degree(end) == 1 {
                set_terminal_atom_coords(&mut new_mol, end, begin)?;
            } else {
                set_terminal_atom_coords(&mut new_mol, begin, end)?;
            }
        }
    }

    // make a guess at the position of the dummy atoms showing the
    // attachment point
    for conf in mol.conformers() {
        let new_conf = new_mol.conformer_mut(conf.id() as i32)?;
        for (&dummy, &nbr) in &dummy_atom_map {
            new_conf.set_atom_pos(final_idx[dummy as usize], conf.atom_pos(nbr)?)?;
        }
    }

    // copy over stereo groups
    let final_atom_map = mol_atom_map
        .iter()
        .map(|(&old, &new)| (old, final_idx[new as usize]))
        .collect();
    copy_stereo_groups(&final_atom_map, mol, &mut new_mol)?;

    // clear computed props and do basic updates on the resulting molecule,
    // but allow unhappiness
    new_mol.clear_computed_props(true);
    new_mol.update_property_cache(false)?;
    Ok(Some(new_mol))
}
