//! Canonical enhanced stereo groups (`CanonicalizeStereoGroups.cpp` of RDKit
//! 2026.09.1).
//!
//! The molecule is expanded to every stereoisomer its OR (or AND) groups
//! stand for; each is ranked in a SMILES-like atom order, and the groups that
//! reproduce the set of isomers are rebuilt from how the stereocentres vary
//! together.
//!
//! Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::cmp::Ordering;
use std::collections::BTreeSet;

use super::cx_write::{RestoreBondDirOption, cx_fields, mol_to_cx_smiles};
use super::write::SmilesWriteParams;
use super::{SmilesParserParams, mol_from_smiles};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::BondStereo;
use crate::graphmol::chirality::stereo_groups::find_meso_centers;
use crate::graphmol::mol_ops::structure::renumber_atoms;
use crate::graphmol::new_canon::{RankOptions, rank_mol_atoms};
use crate::graphmol::romol::ROMol;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result};

/// `StereoGroupAbsOptions`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoGroupAbsOptions {
    /// ABS groups are only added when other stereo groups are, so a molecule
    /// never ends up with just an ABS group.
    OnlyIncludeWhenOtherGroupsExist = 0,
    /// Only AND and OR groups are produced.
    NeverInclude,
    /// The ABS group is always added.
    AlwaysInclude,
}

/// Why the canonicalization stopped: `RigorousEnhancedStereoException` (the
/// molecule is left as it was) or a real error.
enum Failure {
    Rigorous,
    Error(RDKitError),
}

impl From<RDKitError> for Failure {
    fn from(e: RDKitError) -> Self {
        Failure::Error(e)
    }
}

/// `buildTree`: the atoms in depth-first order, neighbours by rank.
fn build_tree(
    atom: u32,
    mol: &ROMol,
    chosen_order: &mut Vec<u32>,
    reverse_order: &mut [i32],
    ranks: &[u32],
) {
    chosen_order.push(atom);
    reverse_order[atom as usize] = chosen_order.len() as i32 - 1;
    let mut nbr_ranks: Vec<(u32, u32)> = mol
        .atom_neighbors(atom)
        .map(|n| (ranks[n as usize], n))
        .collect();
    nbr_ranks.sort_unstable();
    for (_, nbr) in nbr_ranks {
        if reverse_order[nbr as usize] == -1 {
            build_tree(nbr, mol, chosen_order, reverse_order, ranks);
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ChiralAtomItem {
    atom_id: u32,
    chiral_type: ChiralType,
}

impl ChiralAtomItem {
    fn key(&self) -> (u32, u8) {
        (self.atom_id, self.chiral_type as u8)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ChiralBondItem {
    stereo: BondStereo,
    bond_id: u32,
    atom_id1: u32,
    atom_id2: u32,
}

impl ChiralBondItem {
    fn key(&self) -> (u32, u32, u8) {
        (self.atom_id1, self.atom_id2, self.stereo as u8)
    }
}

/// The chiral centres and atropisomeric bonds outside the stereo groups of
/// one expanded stereoisomer.
#[derive(Debug, Clone)]
struct RankedValue {
    atoms: Vec<ChiralAtomItem>,
    /// in the order they were added; see [`RankedValue::sorted_bonds`]
    bonds: Vec<ChiralBondItem>,
}

impl RankedValue {
    fn sorted_bonds(&self) -> Vec<ChiralBondItem> {
        let mut b = self.bonds.clone();
        b.sort_by_key(ChiralBondItem::key);
        b
    }

    /// Same chiral atoms and bonds, whatever their configuration.
    fn equivalent_to(&self, other: &RankedValue) -> bool {
        self.atoms.len() == other.atoms.len()
            && self.bonds.len() == other.bonds.len()
            && self
                .atoms
                .iter()
                .zip(&other.atoms)
                .all(|(a, b)| a.atom_id == b.atom_id)
            && self
                .bonds
                .iter()
                .zip(&other.bonds)
                .all(|(a, b)| a.atom_id1 == b.atom_id1 && a.atom_id2 == b.atom_id2)
    }
}

impl PartialEq for RankedValue {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for RankedValue {}

impl PartialOrd for RankedValue {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RankedValue {
    fn cmp(&self, other: &Self) -> Ordering {
        self.atoms
            .len()
            .cmp(&other.atoms.len())
            .then(self.bonds.len().cmp(&other.bonds.len()))
            .then_with(|| {
                self.atoms
                    .iter()
                    .map(ChiralAtomItem::key)
                    .cmp(other.atoms.iter().map(ChiralAtomItem::key))
            })
            .then_with(|| {
                self.bonds
                    .iter()
                    .map(ChiralBondItem::key)
                    .cmp(other.bonds.iter().map(ChiralBondItem::key))
            })
    }
}

/// The chiral type of atom `index` in each ranked value.
fn atom_types(all: &BTreeSet<RankedValue>, index: usize) -> impl Iterator<Item = ChiralType> + '_ {
    all.iter().map(move |r| r.atoms[index].chiral_type)
}

/// The stereo of the `index`-th (sorted) bond in each ranked value.
fn bond_stereos(all: &BTreeSet<RankedValue>, index: usize) -> Vec<BondStereo> {
    all.iter().map(|r| r.sorted_bonds()[index].stereo).collect()
}

fn does_atom_chirality_vary(all: &BTreeSet<RankedValue>, index: usize) -> bool {
    let first = all.iter().next().expect("not empty").atoms[index].chiral_type;
    atom_types(all, index).any(|t| t != first)
}

fn does_bond_stereo_vary(all: &BTreeSet<RankedValue>, index: usize) -> bool {
    let stereos = bond_stereos(all, index);
    stereos.iter().any(|&s| s != stereos[0])
}

fn do_two_atoms_vary_the_same(all: &BTreeSet<RankedValue>, i1: usize, i2: usize) -> bool {
    let first = all.iter().next().expect("not empty");
    let (f1, f2) = (first.atoms[i1].chiral_type, first.atoms[i2].chiral_type);
    all.iter()
        .all(|r| (f1 == r.atoms[i1].chiral_type) == (f2 == r.atoms[i2].chiral_type))
}

fn do_atom_and_bond_vary_the_same(all: &BTreeSet<RankedValue>, atom: usize, bond: usize) -> bool {
    let first = all.iter().next().expect("not empty");
    let f1 = first.atoms[atom].chiral_type;
    let f2 = first.sorted_bonds()[bond].stereo;
    all.iter()
        .all(|r| (f1 == r.atoms[atom].chiral_type) == (f2 == r.sorted_bonds()[bond].stereo))
}

fn do_two_bonds_vary_the_same(all: &BTreeSet<RankedValue>, b1: usize, b2: usize) -> bool {
    let first = all.iter().next().expect("not empty").sorted_bonds();
    let (f1, f2) = (first[b1].stereo, first[b2].stereo);
    all.iter().all(|r| {
        let b = r.sorted_bonds();
        (f1 == b[b1].stereo) == (f2 == b[b2].stereo)
    })
}

fn count_swaps(nbrs: &[i32]) -> usize {
    let mut swaps = 0;
    for i in 0..nbrs.len() {
        for j in i + 1..nbrs.len() {
            if nbrs[i] > nbrs[j] {
                swaps += 1;
            }
        }
    }
    swaps
}

/// `getAtomsToInvert2`: the (new) indices of the chiral atoms whose tag a
/// SMILES written in `new_order` would invert. `renumberAtoms` keeps the bond
/// order around every atom, so it does not invert them itself.
fn get_atoms_to_invert(mol: &ROMol, new_order: &[u32], reversed_order: &[i32]) -> Vec<u32> {
    let mut res = Vec::new();
    for (n_idx, &o_idx) in new_order.iter().enumerate() {
        if mol.atom(o_idx).chiral_tag() == ChiralType::Unspecified {
            continue;
        }
        let nbrs: Vec<i32> = mol
            .atom_neighbors(o_idx)
            .map(|n| reversed_order[n as usize])
            .collect();
        if count_swaps(&nbrs) % 2 == 1 {
            res.push(n_idx as u32);
        }
    }
    res
}

fn is_tetrahedral(t: ChiralType) -> bool {
    matches!(t, ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW)
}

fn is_atropisomeric(s: BondStereo) -> bool {
    matches!(s, BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW)
}

/// `addSingleAbsGroup`: all chiral centres go into one ABS group (none if
/// there are no chiral centres).
fn add_single_abs_group(mol: &mut ROMol) -> Result<()> {
    let atoms: Vec<u32> = (0..mol.num_atoms())
        .filter(|&i| is_tetrahedral(mol.atom(i).chiral_tag()))
        .collect();
    let bonds: Vec<u32> = (0..mol.num_bonds())
        .filter(|&i| is_atropisomeric(mol.bond(i).stereo()))
        .collect();
    let mut groups = Vec::new();
    if !atoms.is_empty() || !bonds.is_empty() {
        groups.push(StereoGroup::new(
            StereoGroupType::Absolute,
            atoms,
            bonds,
            0,
        )?);
    }
    mol.set_stereo_groups(groups)
}

fn invert_tetrahedral(t: ChiralType) -> ChiralType {
    match t {
        ChiralType::TetrahedralCW => ChiralType::TetrahedralCCW,
        ChiralType::TetrahedralCCW => ChiralType::TetrahedralCW,
        other => other,
    }
}

/// `canonicalizeStereoGroups_internal`
fn canonicalize_internal(
    mol: &mut ROMol,
    group_type: StereoGroupType,
    output_abs: StereoGroupAbsOptions,
    max_stereo_groups: u32,
) -> std::result::Result<(), Failure> {
    let mut groups_to_process: Vec<StereoGroup> = Vec::new();
    let mut and_groups_to_keep: Vec<StereoGroup> = Vec::new();
    for grp in mol.stereo_groups() {
        if grp.group_type() == group_type {
            groups_to_process.push(grp.clone());
        } else if group_type == StereoGroupType::Or && grp.group_type() == StereoGroupType::And {
            and_groups_to_keep.push(grp.clone());
        }
    }
    mol.set_stereo_groups(and_groups_to_keep.clone())?;

    if max_stereo_groups > 0 && groups_to_process.len() > max_stereo_groups as usize {
        return Err(Failure::Rigorous);
    }

    let mut all_ranked: BTreeSet<RankedValue> = BTreeSet::new();
    let mut best_new_mol: Option<ROMol> = None;
    let new_mol_count = 1u64 << groups_to_process.len();
    for mol_index in 0..new_mol_count {
        let mut new_mol = mol.clone();
        for (grp_index, grp) in groups_to_process.iter().enumerate() {
            if mol_index & (1 << grp_index) == 0 {
                continue;
            }
            for &a in grp.atoms() {
                let tag = mol.atom(a).chiral_tag();
                if is_tetrahedral(tag) {
                    new_mol.atom_mut(a).set_chiral_tag(invert_tetrahedral(tag));
                }
            }
            for &b in grp.bonds() {
                match mol.bond(b).stereo() {
                    BondStereo::StereoAtropCW => {
                        new_mol.bond_mut(b).set_stereo(BondStereo::StereoAtropCCW)?
                    }
                    BondStereo::StereoAtropCCW => {
                        new_mol.bond_mut(b).set_stereo(BondStereo::StereoAtropCW)?
                    }
                    _ => {}
                }
            }
        }
        if !and_groups_to_keep.is_empty() {
            canonicalize_internal(
                &mut new_mol,
                StereoGroupType::And,
                StereoGroupAbsOptions::NeverInclude,
                max_stereo_groups,
            )?;
        }
        let ranks = rank_mol_atoms(
            &mut new_mol,
            RankOptions {
                break_ties: true,
                include_chirality: true,
                include_isotopes: false,
                include_atom_maps: true,
                include_chiral_presence: true,
                include_stereo_groups: true,
                use_non_stereo_ranks: true,
                include_ring_stereo: true,
            },
        )?;

        // the atom order of a SMILES for the whole molecule (not the
        // fragments, so not the same as the one of an actual SMILES)
        let n = new_mol.num_atoms() as usize;
        let mut chosen_order: Vec<u32> = Vec::with_capacity(n);
        let mut reversed_order = vec![-1i32; n];
        loop {
            let start = (0..n)
                .filter(|&i| reversed_order[i] == -1)
                .min_by_key(|&i| (ranks[i], i));
            let Some(start) = start else { break };
            build_tree(
                start as u32,
                &new_mol,
                &mut chosen_order,
                &mut reversed_order,
                &ranks,
            );
        }
        if chosen_order.len() != n {
            return Err(RDKitError::Value("atomOrdering size mismatch".into()).into());
        }
        let atoms_to_invert = get_atoms_to_invert(&new_mol, &chosen_order, &reversed_order);
        let new_mol = renumber_atoms(&new_mol, &chosen_order)?;

        let mut in_groups_atoms = vec![false; new_mol.num_atoms() as usize];
        let mut in_groups_bonds = vec![false; new_mol.num_bonds() as usize];
        for grp in new_mol.stereo_groups() {
            for &a in grp.atoms() {
                in_groups_atoms[a as usize] = true;
            }
            for &b in grp.bonds() {
                in_groups_bonds[b as usize] = true;
            }
        }
        let mut value = RankedValue {
            atoms: Vec::new(),
            bonds: Vec::new(),
        };
        for idx in 0..new_mol.num_atoms() {
            let tag = new_mol.atom(idx).chiral_tag();
            if is_tetrahedral(tag) && !in_groups_atoms[idx as usize] {
                let chiral_type = if atoms_to_invert.contains(&idx) {
                    invert_tetrahedral(tag)
                } else {
                    tag
                };
                value.atoms.push(ChiralAtomItem {
                    atom_id: idx,
                    chiral_type,
                });
            }
        }
        for idx in 0..new_mol.num_bonds() {
            let b = new_mol.bond(idx);
            if is_atropisomeric(b.stereo()) && !in_groups_bonds[idx as usize] {
                let (mut a1, mut a2) = (b.begin_atom_idx(), b.end_atom_idx());
                if a1 > a2 {
                    std::mem::swap(&mut a1, &mut a2);
                }
                value.bonds.push(ChiralBondItem {
                    stereo: b.stereo(),
                    bond_id: idx,
                    atom_id1: a1,
                    atom_id2: a2,
                });
            }
        }
        if let Some(first) = all_ranked.iter().next()
            && !value.equivalent_to(first)
        {
            return Err(Failure::Rigorous);
        }
        let inserted = all_ranked.insert(value.clone());
        if inserted && all_ranked.iter().next() == Some(&value) {
            best_new_mol = Some(new_mol);
        }
    }

    let mut best = best_new_mol.expect("at least one stereoisomer");
    let first = all_ranked.iter().next().expect("not empty").clone();
    let first_bonds = first.sorted_bonds();
    let mut atoms_done = vec![false; first.atoms.len()];
    let mut bonds_done = vec![false; first_bonds.len()];
    let mut abs_atoms: Vec<u32> = Vec::new();
    let mut abs_bonds: Vec<u32> = Vec::new();
    let mut new_groups: Vec<StereoGroup> = Vec::new();

    if all_ranked.len() == 1 {
        // no variation: everything is absolute
        if output_abs == StereoGroupAbsOptions::NeverInclude {
            *mol = best;
            return Ok(());
        }
        abs_atoms.extend(first.atoms.iter().map(|a| a.atom_id));
        abs_bonds.extend(first_bonds.iter().map(|b| b.bond_id));
    } else {
        let mut group_count = 0;
        for index1 in 0..first.atoms.len() {
            if atoms_done[index1] {
                continue;
            }
            let atom_index1 = first.atoms[index1].atom_id;
            if !does_atom_chirality_vary(&all_ranked, index1) {
                if output_abs != StereoGroupAbsOptions::NeverInclude {
                    abs_atoms.push(atom_index1);
                }
                atoms_done[index1] = true;
                continue;
            }
            let mut atoms_to_add = vec![atom_index1];
            atoms_done[index1] = true;
            // the other atoms and bonds that vary the same way
            for index2 in index1 + 1..first.atoms.len() {
                if !atoms_done[index2] && do_two_atoms_vary_the_same(&all_ranked, index1, index2) {
                    atoms_to_add.push(first.atoms[index2].atom_id);
                    atoms_done[index2] = true;
                }
            }
            let mut bonds_to_add: Vec<u32> = Vec::new();
            for index2 in 0..first_bonds.len() {
                if !bonds_done[index2]
                    && do_atom_and_bond_vary_the_same(&all_ranked, index1, index2)
                {
                    bonds_to_add.push(first_bonds[index2].bond_id);
                    bonds_done[index2] = true;
                }
            }
            atoms_to_add.sort_unstable();
            bonds_to_add.sort_unstable();
            group_count += 1;
            new_groups.push(StereoGroup::new(
                group_type,
                atoms_to_add,
                bonds_to_add,
                group_count,
            )?);
        }
        // the groups that only involve bonds
        for index1 in 0..first_bonds.len() {
            if bonds_done[index1] {
                continue;
            }
            let bond_index1 = first_bonds[index1].bond_id;
            if !does_bond_stereo_vary(&all_ranked, index1) {
                if output_abs != StereoGroupAbsOptions::NeverInclude {
                    abs_bonds.push(bond_index1);
                }
                bonds_done[index1] = true;
                continue;
            }
            let mut bonds_to_add = vec![bond_index1];
            bonds_done[index1] = true;
            for index2 in index1 + 1..first_bonds.len() {
                if !bonds_done[index2] && do_two_bonds_vary_the_same(&all_ranked, index1, index2) {
                    bonds_to_add.push(first_bonds[index2].bond_id);
                    bonds_done[index2] = true;
                }
            }
            bonds_to_add.sort_unstable();
            group_count += 1;
            new_groups.push(StereoGroup::new(
                group_type,
                Vec::new(),
                bonds_to_add,
                group_count,
            )?);
        }
    }

    // keep the groups of the best molecule, if it got them from a call for
    // the other kind of stereo groups
    new_groups.extend(best.stereo_groups().iter().cloned());

    let add_abs = output_abs == StereoGroupAbsOptions::AlwaysInclude
        || (output_abs == StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist
            && !new_groups.is_empty());
    if add_abs && (!abs_atoms.is_empty() || !abs_bonds.is_empty()) {
        abs_atoms.sort_unstable();
        abs_bonds.sort_unstable();
        new_groups.push(StereoGroup::new(
            StereoGroupType::Absolute,
            abs_atoms,
            abs_bonds,
            0,
        )?);
    }
    best.set_stereo_groups(new_groups)?;
    *mol = best;
    Ok(())
}

/// `canonicalizeStereoGroups`: rationalizes the enhanced stereo of `mol` to a
/// canonical form. More than `max_stereo_groups` OR (or AND) groups leave the
/// molecule as it was.
pub fn canonicalize_stereo_groups(
    mol: &ROMol,
    output_abs: StereoGroupAbsOptions,
    max_stereo_groups: u32,
) -> Result<ROMol> {
    let mut mol = mol.clone();
    let groups = mol.stereo_groups();
    if groups.is_empty()
        || (groups.len() == 1 && groups[0].group_type() == StereoGroupType::Absolute)
    {
        if output_abs == StereoGroupAbsOptions::AlwaysInclude {
            add_single_abs_group(&mut mol)?;
        } else {
            mol.set_stereo_groups(Vec::new())?;
        }
        return Ok(mol);
    }

    // a single non-ABS group of at most two atoms that holds all the chiral
    // centres is not processed; it is dropped if the molecule is meso
    if groups.len() == 1 {
        let sg = &groups[0];
        if sg.atoms().len() <= 2 && sg.bonds().is_empty() {
            let simple = (0..mol.num_atoms())
                .all(|i| !is_tetrahedral(mol.atom(i).chiral_tag()) || sg.atoms().contains(&i))
                && (0..mol.num_bonds())
                    .all(|i| !is_atropisomeric(mol.bond(i).stereo()) || sg.bonds().contains(&i));
            if simple {
                let meso = find_meso_centers(&mut mol, true, true)?;
                if meso.len() == 1 {
                    if output_abs == StereoGroupAbsOptions::AlwaysInclude {
                        add_single_abs_group(&mut mol)?;
                    } else {
                        mol.set_stereo_groups(Vec::new())?;
                    }
                }
                return Ok(mol);
            }
        }
    }

    let found_or_group = mol
        .stereo_groups()
        .iter()
        .any(|g| g.group_type() == StereoGroupType::Or);
    if mol.needs_update_property_cache() {
        mol.update_property_cache(true)?;
    }

    // the rankings without stereo do not change with the expansion
    let ranks = rank_mol_atoms(
        &mut mol,
        RankOptions {
            break_ties: false,
            include_chirality: false,
            include_isotopes: false,
            include_atom_maps: true,
            include_chiral_presence: true,
            include_stereo_groups: false,
            use_non_stereo_ranks: false,
            include_ring_stereo: false,
        },
    )?;
    for (idx, rank) in ranks.iter().enumerate() {
        mol.atom_mut(idx as u32)
            .set_prop(cp::_CanonicalRankingNumber, *rank, false)?;
    }

    let saved_groups = mol.stereo_groups().to_vec();
    let group_type = if found_or_group {
        StereoGroupType::Or
    } else {
        StereoGroupType::And
    };
    match canonicalize_internal(&mut mol, group_type, output_abs, max_stereo_groups) {
        Ok(()) => {}
        Err(Failure::Rigorous) => {
            mol.set_stereo_groups(saved_groups)?;
            return Ok(mol);
        }
        Err(Failure::Error(e)) => return Err(e),
    }

    // fix up the molecule: round trip through SMILES
    mol.clear_computed_props(true);
    let params = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    let smiles = mol_to_cx_smiles(&mol, &params, cx_fields::ALL, RestoreBondDirOption::Clear)?;
    let parse = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    mol_from_smiles(&smiles, &parse)?
        .ok_or_else(|| RDKitError::Value("could not parse the canonicalized SMILES".into()))
}
