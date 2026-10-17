//! Port of RDKit `Code/GraphMol/MolOps.{h,cpp}`: sanitization and friends.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

pub use crate::graphmol::chirality::{clear_dir_flags, clear_single_bond_dir_flags};

use crate::graphmol::aromaticity::{AromaticityModel, set_aromaticity};
use crate::graphmol::atom::{ChiralType, HybridizationType, ValenceType};
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::conjug_hybrid::{set_conjugation, set_hybridization};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, find_sssr, symmetrize_sssr};
use crate::graphmol::kekulize::{DEFAULT_MAX_BACKTRACKS, kekulize, kekulize_if_possible};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::StereoGroup;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result};

/// `MolOps::SanitizeFlags`
pub mod sanitize_flags {
    pub const SANITIZE_NONE: u32 = 0x0;
    pub const SANITIZE_CLEANUP: u32 = 0x1;
    pub const SANITIZE_PROPERTIES: u32 = 0x2;
    pub const SANITIZE_SYMMRINGS: u32 = 0x4;
    pub const SANITIZE_KEKULIZE: u32 = 0x8;
    pub const SANITIZE_FINDRADICALS: u32 = 0x10;
    pub const SANITIZE_SETAROMATICITY: u32 = 0x20;
    pub const SANITIZE_SETCONJUGATION: u32 = 0x40;
    pub const SANITIZE_SETHYBRIDIZATION: u32 = 0x80;
    pub const SANITIZE_CLEANUPCHIRALITY: u32 = 0x100;
    pub const SANITIZE_ADJUSTHS: u32 = 0x200;
    pub const SANITIZE_CLEANUP_ORGANOMETALLICS: u32 = 0x400;
    pub const SANITIZE_CLEANUPATROPISOMERS: u32 = 0x800;
    pub const SANITIZE_ALL: u32 = 0xFFF_FFFF;
}
use sanitize_flags::*;

pub mod add_hs;
pub mod frags;
pub mod haptic;
pub mod query_hs;
pub mod structure;
use crate::rdgeneral::RDProps;
pub use add_hs::{RemoveHsParameters, remove_hs};
pub use frags::{get_mol_frags_as_mols, get_mol_frags_atoms, get_mol_frags_mapping};
pub use structure::{get_formal_charge, remove_stereochemistry, renumber_atoms};

/// `QueryOps::isMetal` for plain atoms (Marvin's "M" pseudo-element).
pub fn is_metal(atomic_num: u32) -> bool {
    const NON_METALS: [u32; 23] = [
        0, 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86,
    ];
    !NON_METALS.contains(&atomic_num)
}

// ------------------------------------------------------------------ cleanUp

fn nitrogens_cleanup(mol: &mut RWMol) -> Result<()> {
    let mut to_consider = Vec::new();
    for aid in 0..mol.num_atoms() {
        let a = mol.atom(aid);
        if a.atomic_num() != 7 || a.formal_charge() != 0 {
            continue;
        }
        if mol.calc_explicit_valence(aid, false)? != 5 {
            continue;
        }
        to_consider.push(aid);
        let arom = mol.atom(aid).is_aromatic();
        mol.atom_mut(aid).set_is_aromatic(false);
        let mut update = false;
        let nbrs: Vec<u32> = mol.atom_neighbors(aid).collect();
        for nbr in nbrs {
            let b = mol.bond_idx_between(aid, nbr).expect("bonded");
            if mol.atom(nbr).atomic_num() == 8
                && mol.atom(nbr).formal_charge() == 0
                && mol.bond(b).bond_type() == BondType::Double
            {
                mol.bond_mut(b).set_bond_type(BondType::Single);
                mol.atom_mut(aid).set_formal_charge(1);
                mol.atom_mut(nbr).set_formal_charge(-1);
                update = true;
                break;
            }
        }
        mol.atom_mut(aid).set_is_aromatic(arom);
        if update {
            mol.calc_explicit_valence(aid, false)?;
        }
    }
    for aid in to_consider {
        let arom = mol.atom(aid).is_aromatic();
        mol.atom_mut(aid).set_is_aromatic(false);
        let mut update = false;
        let nbrs: Vec<u32> = mol.atom_neighbors(aid).collect();
        for nbr in nbrs {
            let b = mol.bond_idx_between(aid, nbr).expect("bonded");
            if mol.atom(nbr).atomic_num() == 7
                && mol.atom(nbr).formal_charge() == 0
                && mol.bond(b).bond_type() == BondType::Triple
            {
                mol.bond_mut(b).set_bond_type(BondType::Double);
                mol.atom_mut(aid).set_formal_charge(1);
                mol.atom_mut(nbr).set_formal_charge(-1);
                update = true;
                break;
            }
        }
        mol.atom_mut(aid).set_is_aromatic(arom);
        if update {
            mol.calc_explicit_valence(aid, false)?;
        }
    }
    Ok(())
}

fn phosphorus_cleanup(mol: &mut RWMol, aid: u32) -> Result<()> {
    if mol.atom(aid).formal_charge() != 0 {
        return Ok(());
    }
    if mol.calc_explicit_valence(aid, false)? == 5 && mol.atom_degree(aid) == 3 {
        let mut dbl_to_o = None;
        let mut o_atom = None;
        let mut has_double_to_c_or_n = false;
        for nbr in mol.atom_neighbors(aid).collect::<Vec<_>>() {
            let b = mol.bond_idx_between(aid, nbr).expect("bonded");
            let na = mol.atom(nbr);
            if na.atomic_num() == 8
                && na.formal_charge() == 0
                && mol.bond(b).bond_type() == BondType::Double
            {
                dbl_to_o = Some(b);
                o_atom = Some(nbr);
            } else if (na.atomic_num() == 6 || na.atomic_num() == 7)
                && mol.atom_degree(nbr) >= 2
                && mol.bond(b).bond_type() == BondType::Double
            {
                has_double_to_c_or_n = true;
            }
        }
        if let (true, Some(b), Some(o)) = (has_double_to_c_or_n, dbl_to_o, o_atom) {
            mol.atom_mut(o).set_formal_charge(-1);
            mol.bond_mut(b).set_bond_type(BondType::Single);
            mol.atom_mut(aid).set_formal_charge(1);
        }
    }
    mol.calc_explicit_valence(aid, false)?;
    Ok(())
}

fn halogen_cleanup(mol: &mut RWMol, aid: u32) -> Result<()> {
    let ev = mol.calc_explicit_valence(aid, false)?;
    if mol.atom(aid).formal_charge() == 0 && matches!(ev, 7 | 5 | 3) {
        if mol
            .atom_neighbors(aid)
            .all(|n| mol.atom(n).atomic_num() == 8)
        {
            let mut fc = 0;
            for b in mol.atom_bond_indices(aid).collect::<Vec<_>>() {
                if mol.bond(b).bond_type() == BondType::Double {
                    mol.bond_mut(b).set_bond_type(BondType::Single);
                    let other = mol.bond(b).other_atom_idx(aid)?;
                    fc += 1;
                    mol.atom_mut(other).set_formal_charge(-1);
                    mol.calc_explicit_valence(other, false)?;
                }
            }
            mol.atom_mut(aid).set_formal_charge(fc);
            mol.calc_explicit_valence(aid, false)?;
        }
    }
    Ok(())
}

/// `MolOps::cleanUp`
pub fn clean_up(mol: &mut RWMol) -> Result<()> {
    nitrogens_cleanup(mol)?;
    for aid in 0..mol.num_atoms() {
        match mol.atom(aid).atomic_num() {
            15 => phosphorus_cleanup(mol, aid)?,
            17 | 35 | 53 => halogen_cleanup(mol, aid)?,
            _ => {}
        }
    }
    Ok(())
}

fn is_hypervalent_non_metal(mol: &mut RWMol, aid: u32) -> Result<bool> {
    let an = mol.atom(aid).atomic_num();
    if is_metal(an) {
        return Ok(false);
    }
    mol.update_atom_property_cache(aid, false)?;
    let a = mol.atom(aid);
    let ev = a.valence(ValenceType::Explicit)? as i32;
    let eff = an as i32 - a.formal_charge();
    if eff <= 0 {
        return Ok(false);
    }
    let pt = PeriodicTable::get();
    if eff as u32 > pt.get_max_atomic_number() {
        return Err(RDKitError::precondition("Atomic number not found"));
    }
    let max_v = *pt.element(eff as u32).valences.last().expect("valences");
    Ok(
        max_v > 0
            && (ev > max_v || (ev == max_v && a.is_aromatic() && mol.total_degree(aid)? == 4)),
    )
}

fn no_dative(atomic_num: u32) -> bool {
    matches!(atomic_num, 1 | 2 | 9 | 10)
}

/// `MolOps::cleanUpOrganometallics`
pub fn clean_up_organometallics(mol: &mut RWMol) -> Result<()> {
    let mut needs_fixing = false;
    for aid in 0..mol.num_atoms() {
        if is_hypervalent_non_metal(mol, aid)? && !no_dative(mol.atom(aid).atomic_num()) {
            if mol.atom_bonds(aid).any(|b| {
                b.bond_type() == BondType::Single
                    && is_metal(
                        mol.atom(b.other_atom_idx(aid).expect("in bond"))
                            .atomic_num(),
                    )
            }) {
                needs_fixing = true;
                break;
            }
        }
    }
    if !needs_fixing {
        return Ok(());
    }
    mol.update_property_cache(false)?;
    let ranks = crate::graphmol::new_canon::rank_mol_atoms(
        mol,
        crate::graphmol::new_canon::RankOptions::default(),
    )?;
    let mut atom_ranks: Vec<(usize, u32)> = ranks.iter().copied().enumerate().collect();
    // std::sort by rank (ranks are unique after tie breaking)
    atom_ranks.sort_by_key(|&(_, r)| r);
    for (aid, _) in atom_ranks {
        metal_bond_cleanup(mol, aid as u32, &ranks)?;
    }
    Ok(())
}

fn num_dative_bonds(mol: &RWMol, aid: u32) -> i32 {
    mol.atom_bonds(aid)
        .filter(|b| {
            matches!(
                b.bond_type(),
                BondType::Dative | BondType::DativeOne | BondType::DativeL | BondType::DativeR
            )
        })
        .count() as i32
}

fn metal_bond_cleanup(mol: &mut RWMol, aid: u32, ranks: &[u32]) -> Result<()> {
    if !(is_hypervalent_non_metal(mol, aid)? && !no_dative(mol.atom(aid).atomic_num())) {
        return Ok(());
    }
    let mut metals: Vec<u32> = mol
        .atom_bonds(aid)
        .filter(|b| b.bond_type() == BondType::Single)
        .map(|b| b.other_atom_idx(aid).expect("in bond"))
        .filter(|&o| is_metal(mol.atom(o).atomic_num()))
        .collect();
    if metals.is_empty() {
        return Ok(());
    }
    metals.sort_by(|&a1, &a2| {
        let (n1, n2) = (num_dative_bonds(mol, a1), num_dative_bonds(mol, a2));
        if n1 == n2 {
            ranks[a2 as usize].cmp(&ranks[a1 as usize])
        } else {
            n1.cmp(&n2)
        }
    });
    if let Some(b) = mol.bond_idx_between(aid, metals[0]) {
        mol.bond_mut(b).set_bond_type(BondType::Dative);
        mol.set_bond_begin_atom_idx(b, aid)?;
        mol.set_bond_end_atom_idx(b, metals[0])?;
    }
    Ok(())
}

// ------------------------------------------------------------------ radicals / Hs

/// `MolOps::assignRadicals`
pub fn assign_radicals(mol: &mut RWMol) {
    let pt = PeriodicTable::get();
    for aid in 0..mol.num_atoms() {
        let a = mol.atom(aid);
        if !a.no_implicit() || a.atomic_num() == 0 {
            continue;
        }
        let an = a.atomic_num();
        let valens = &pt.element(an).valences;
        let chg = a.formal_charge();
        let n_outer = pt.element(an).n_outer_elecs;

        let num_radicals = if valens.len() != 1 || valens[0] != -1 {
            let mut accum: f64 = mol
                .atom_bonds(aid)
                .map(|b| b.valence_contrib(mol, aid))
                .sum();
            accum += f64::from(a.num_explicit_hs());
            let total_valence = (accum + 0.1) as i32;
            let base_count = if an == 1 || an == 2 { 2 } else { 8 };
            let mut nr = base_count - n_outer - total_valence + chg;
            if nr < 0 {
                nr = 0;
                if valens.len() > 1 {
                    for &v in valens {
                        if v - total_valence + chg >= 0 {
                            nr = v - total_valence + chg;
                            break;
                        }
                    }
                }
            }
            let nr2 = n_outer - total_valence - chg;
            if nr2 >= 0 {
                nr = nr.min(nr2);
            }
            nr
        } else if mol.atom_degree(aid) > 0 {
            0
        } else {
            let mut n_valence = n_outer - chg;
            if n_valence < 0 {
                n_valence = 0;
                crate::rd_log_warning!(
                    "Unusual charge on atom {} number of radical electrons set to zero",
                    aid
                );
            }
            n_valence % 2
        };
        mol.atom_mut(aid)
            .set_num_radical_electrons(num_radicals as u32);
    }
}

fn adjust_hs_impl(mol: &mut RWMol, atoms_to_adjust: Option<&[bool]>) -> Result<()> {
    for aid in 0..mol.num_atoms() {
        if let Some(sel) = atoms_to_adjust {
            if !sel[aid as usize] {
                continue;
            }
        }
        let orig_implicit = mol.atom(aid).valence(ValenceType::Implicit)? as i32;
        mol.calc_explicit_valence(aid, false)?;
        let orig_explicit = mol.atom(aid).num_explicit_hs() as i32;
        let new_implicit = mol.calc_implicit_valence(aid, false)?;
        if new_implicit < orig_implicit {
            mol.atom_mut(aid)
                .set_num_explicit_hs((orig_explicit + orig_implicit - new_implicit) as u32);
            mol.calc_explicit_valence(aid, false)?;
        }
    }
    Ok(())
}

/// `MolOps::adjustHs`
pub fn adjust_hs(mol: &mut RWMol) -> Result<()> {
    adjust_hs_impl(mol, None)
}

fn include_aromatic_atoms(mol: &RWMol, sel: &mut [bool]) {
    for a in mol.atoms() {
        if a.is_aromatic() {
            sel[a.idx() as usize] = true;
        }
    }
    for b in mol.bonds() {
        if b.is_aromatic() {
            sel[b.begin_atom_idx() as usize] = true;
            sel[b.end_atom_idx() as usize] = true;
        }
    }
}

// ------------------------------------------------------------------ stereo cleanups

fn regroup(
    mol: &RWMol,
    sg: &StereoGroup,
    atoms: Vec<u32>,
    bonds: Vec<u32>,
    read_id: u32,
) -> Result<StereoGroup> {
    let _ = mol;
    StereoGroup::new(sg.group_type(), atoms, bonds, read_id)
}

/// `Chirality::cleanupStereoGroups`
pub fn cleanup_stereo_groups(mol: &mut RWMol) -> Result<()> {
    let mut new = Vec::new();
    for sg in mol.stereo_groups() {
        let mut keep = true;
        let mut ok_atoms = Vec::new();
        let mut ok_bonds = Vec::new();
        for &a in sg.atoms() {
            if mol.atom(a).chiral_tag() == ChiralType::Unspecified {
                keep = false;
            } else {
                ok_atoms.push(a);
            }
        }
        for &b in sg.bonds() {
            if !matches!(
                mol.bond(b).stereo(),
                BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW
            ) {
                keep = false;
            } else {
                ok_bonds.push(b);
            }
        }
        if keep {
            new.push(sg.clone());
        } else if !ok_atoms.is_empty() {
            new.push(regroup(mol, sg, ok_atoms, ok_bonds, sg.read_id())?);
        }
    }
    mol.set_stereo_groups(new)
}

/// `MolOps::cleanupChirality`
pub fn cleanup_chirality(mol: &mut RWMol) -> Result<()> {
    let mut need_cleanup = false;
    for aid in 0..mol.num_atoms() {
        let tag = mol.atom(aid).chiral_tag();
        let hyb = mol.atom(aid).hybridization();
        let perm = mol
            .atom(aid)
            .props
            .get(cp::_chiralPermutation)
            .and_then(RDValue::as_uint)
            .unwrap_or(0);
        let limit = match tag {
            ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW => {
                if hyb != HybridizationType::SP3 {
                    mol.atom_mut(aid).set_chiral_tag(ChiralType::Unspecified);
                    need_cleanup = true;
                }
                continue;
            }
            ChiralType::Tetrahedral => {
                if hyb != HybridizationType::SP3 {
                    mol.atom_mut(aid).set_chiral_tag(ChiralType::Unspecified);
                    need_cleanup = true;
                    continue;
                }
                2
            }
            ChiralType::SquarePlanar | ChiralType::TrigonalBipyramidal | ChiralType::Octahedral => {
                let (lo, hi, lim) = match tag {
                    ChiralType::SquarePlanar => (2, 4, 3),
                    ChiralType::TrigonalBipyramidal => (2, 5, 20),
                    _ => (2, 6, 30),
                };
                let degree = mol.total_degree(aid)?;
                if degree < lo || degree > hi {
                    mol.atom_mut(aid).set_chiral_tag(ChiralType::Unspecified);
                    continue;
                }
                lim
            }
            _ => continue,
        };
        if perm > limit {
            mol.atom_mut(aid)
                .props
                .set_val(cp::_chiralPermutation, 0u32);
        }
    }
    if need_cleanup {
        cleanup_stereo_groups(mol)?;
    }
    Ok(())
}

/// `Atropisomers::cleanupAtropisomerStereoGroups`
fn cleanup_atropisomer_stereo_groups(mol: &mut RWMol) -> Result<()> {
    let mut new = Vec::new();
    for sg in mol.stereo_groups() {
        let mut ok_atoms = Vec::new();
        let mut ok_bonds: Vec<u32> = Vec::new();
        for &a in sg.atoms() {
            let mut found = false;
            for b in mol.atom_bonds(a) {
                if matches!(
                    b.stereo(),
                    BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW
                ) {
                    found = true;
                    if !ok_bonds.contains(&b.idx()) {
                        ok_bonds.push(b.idx());
                    }
                }
            }
            if !found {
                ok_atoms.push(a);
            }
        }
        if ok_bonds.is_empty() {
            new.push(sg.clone());
        } else {
            new.push(regroup(mol, sg, ok_atoms, ok_bonds, 0)?);
        }
    }
    mol.set_stereo_groups(new)
}

/// `MolOps::Hybridizations(mol)`: current hybridizations, computed on a
/// copy when unset.
fn hybridizations(mol: &RWMol) -> Result<Vec<HybridizationType>> {
    if mol.num_atoms() == 0 {
        return Ok(Vec::new());
    }
    if mol.atom(0).hybridization() != HybridizationType::Unspecified {
        return Ok(mol.atoms().iter().map(|a| a.hybridization()).collect());
    }
    let mut copy = ROMol::copy_from(mol, false, -1);
    let _ = sanitize_mol(
        &mut copy,
        SANITIZE_SETCONJUGATION | SANITIZE_SETHYBRIDIZATION,
    );
    Ok(copy.atoms().iter().map(|a| a.hybridization()).collect())
}

/// `MolOps::cleanupAtropisomers`
pub fn cleanup_atropisomers(mol: &mut RWMol) -> Result<()> {
    if !mol.bonds().iter().any(|b| {
        matches!(
            b.stereo(),
            BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW
        )
    }) {
        return Ok(());
    }
    let hybs = hybridizations(mol)?;
    let mut need = false;
    for bi in 0..mol.num_bonds() {
        if !matches!(
            mol.bond(bi).stereo(),
            BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW
        ) {
            continue;
        }
        if !mol.ring_info().is_sssr_or_better() {
            find_sssr(mol, false, false)?;
        }
        let (b, e) = (mol.bond(bi).begin_atom_idx(), mol.bond(bi).end_atom_idx());
        let ri = mol.ring_info();
        if hybs[b as usize] != HybridizationType::SP2
            || hybs[e as usize] != HybridizationType::SP2
            || (ri.num_bond_rings(bi)? > 0 && ri.min_bond_ring_size(bi)? < 8)
        {
            mol.bond_mut(bi).set_stereo(BondStereo::StereoNone)?;
            need = true;
        }
    }
    if need {
        cleanup_atropisomer_stereo_groups(mol)?;
    }
    Ok(())
}

// ------------------------------------------------------------------ sanitizeMol

fn kekulize_for_sanitize(mol: &mut RWMol) -> Result<()> {
    if !kekulize_if_possible(mol, true, false, DEFAULT_MAX_BACKTRACKS)? {
        kekulize(mol, true, true, DEFAULT_MAX_BACKTRACKS)?;
    }
    Ok(())
}

/// `MolOps::sanitizeMol(mol, operationThatFailed, sanitizeOps)`. On error,
/// returns the failing operation flag together with the error.
pub fn sanitize_mol(mol: &mut RWMol, ops: u32) -> std::result::Result<(), (u32, RDKitError)> {
    let step = |flag: u32, r: Result<()>| r.map_err(|e| (flag, e));
    mol.clear_computed_props(true);
    if ops & SANITIZE_CLEANUP != 0 {
        step(SANITIZE_CLEANUP, clean_up(mol))?;
    }
    if ops & SANITIZE_CLEANUP_ORGANOMETALLICS != 0 {
        step(
            SANITIZE_CLEANUP_ORGANOMETALLICS,
            clean_up_organometallics(mol),
        )?;
    }
    step(
        SANITIZE_PROPERTIES,
        mol.update_property_cache(ops & SANITIZE_PROPERTIES != 0),
    )?;
    let track =
        ops & SANITIZE_ADJUSTHS != 0 && ops & (SANITIZE_KEKULIZE | SANITIZE_SETAROMATICITY) != 0;
    let mut to_adjust = vec![false; if track { mol.num_atoms() as usize } else { 0 }];
    if track {
        include_aromatic_atoms(mol, &mut to_adjust);
    }
    if ops & SANITIZE_KEKULIZE != 0 {
        step(SANITIZE_KEKULIZE, kekulize_for_sanitize(mol))?;
    }
    if ops & SANITIZE_SYMMRINGS != 0 {
        step(
            SANITIZE_SYMMRINGS,
            symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, false, false, false).map(|_| ()),
        )?;
    }
    if ops & SANITIZE_FINDRADICALS != 0 {
        assign_radicals(mol);
    }
    if ops & SANITIZE_SETAROMATICITY != 0 {
        step(
            SANITIZE_SETAROMATICITY,
            set_aromaticity(mol, AromaticityModel::Default).map(|_| ()),
        )?;
        if track {
            include_aromatic_atoms(mol, &mut to_adjust);
        }
    }
    if ops & SANITIZE_SETCONJUGATION != 0 {
        set_conjugation(mol);
    }
    if ops & SANITIZE_SETHYBRIDIZATION != 0 {
        set_hybridization(mol);
    }
    if ops & SANITIZE_CLEANUPATROPISOMERS != 0 {
        step(SANITIZE_CLEANUPATROPISOMERS, cleanup_atropisomers(mol))?;
    }
    if ops & SANITIZE_CLEANUPCHIRALITY != 0 {
        step(SANITIZE_CLEANUPCHIRALITY, cleanup_chirality(mol))?;
    }
    if ops & SANITIZE_ADJUSTHS != 0 {
        let sel = if track {
            Some(to_adjust.as_slice())
        } else {
            None
        };
        step(SANITIZE_ADJUSTHS, adjust_hs_impl(mol, sel))?;
    }
    if ops & SANITIZE_PROPERTIES != 0 {
        step(SANITIZE_PROPERTIES, mol.update_property_cache(true))?;
    }
    Ok(())
}

/// The sanitize / removeHs / stereo steps `MolFromSmiles` runs after
/// parsing. removeHs and stereochemistry assignment are ported next.
/// The sanitize/removeHs/stereo part of `SmilesToMol` (`MolFromSmiles` in
/// SmilesParse.cpp), run when `sanitize || removeHs`.
pub fn smiles_post_parse_sanitize(mol: &mut RWMol, sanitize: bool, remove_hs: bool) -> Result<()> {
    if remove_hs {
        let rhp = RemoveHsParameters {
            update_explicit_count: true,
            ..Default::default()
        };
        add_hs::remove_hs(mol, &rhp, sanitize)?;
    } else if sanitize {
        sanitize_mol(mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
    }
    if mol.has_prop("_needsDetectBondStereo") {
        // the first 2D conformer, else the first 3D one
        let confs = mol.conformers();
        let conf = confs
            .iter()
            .position(|c| !c.is_3d())
            .or_else(|| confs.iter().position(|c| c.is_3d()));
        if conf.is_some() {
            clear_single_bond_dir_flags(mol, false);
        }
        crate::graphmol::chirality::bond_dirs::set_double_bond_neighbor_directions(mol, conf)?;
    }
    mol.clear_prop("_needsDetectBondStereo");
    crate::graphmol::chirality::legacy::assign_stereochemistry(mol, true, true, true)
}
