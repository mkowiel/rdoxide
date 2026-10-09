//! Hydrogen removal (`MolOps::removeHs` and helpers); `addHs` follows later.
//!
//! Ported from RDKit `Code/GraphMol/AddHs.cpp` (2026.09.1).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! SubstanceGroups are not ported yet; the SGroup-related checks are no-ops
//! until they are (a molecule without SGroups behaves identically).

use std::collections::BTreeMap;

use super::{sanitize_flags::SANITIZE_ALL, sanitize_mol};
use crate::graphmol::atom::{ChiralType, is_aromatic_atom};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::nontetrahedral::has_non_tetrahedral_stereo;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result};

/// `MolOps::RemoveHsParameters`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemoveHsParameters {
    pub remove_degree_zero: bool,
    pub remove_higher_degrees: bool,
    pub remove_only_h_neighbors: bool,
    pub remove_isotopes: bool,
    pub remove_and_track_isotopes: bool,
    pub remove_dummy_neighbors: bool,
    pub remove_defining_bond_stereo: bool,
    pub remove_with_wedged_bond: bool,
    pub remove_with_query: bool,
    pub remove_mapped: bool,
    pub remove_in_sgroups: bool,
    pub show_warnings: bool,
    pub remove_nonimplicit: bool,
    pub update_explicit_count: bool,
    pub remove_hydrides: bool,
    pub remove_nontetrahedral_neighbors: bool,
}

impl Default for RemoveHsParameters {
    fn default() -> Self {
        Self {
            remove_degree_zero: false,
            remove_higher_degrees: false,
            remove_only_h_neighbors: false,
            remove_isotopes: false,
            remove_and_track_isotopes: false,
            remove_dummy_neighbors: false,
            remove_defining_bond_stereo: false,
            remove_with_wedged_bond: true,
            remove_with_query: false,
            remove_mapped: true,
            remove_in_sgroups: true,
            show_warnings: true,
            remove_nonimplicit: true,
            update_explicit_count: false,
            remove_hydrides: false,
            remove_nontetrahedral_neighbors: false,
        }
    }
}

/// anonymous-namespace `getIsoMap`
fn get_iso_map(mol: &mut RWMol) -> BTreeMap<u32, Vec<u32>> {
    for i in 0..mol.num_atoms() {
        mol.atom_mut(i).clear_prop(cp::_isotopicHs);
    }
    let mut iso_map: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for b in mol.bonds() {
        let ba = mol.atom(b.begin_atom_idx());
        let ea = mol.atom(b.end_atom_idx());
        let (ha, iso) = if ba.atomic_num() == 1 && ba.isotope() != 0 && ea.atomic_num() != 1 {
            (ea.idx(), ba.isotope())
        } else if ea.atomic_num() == 1 && ea.isotope() != 0 && ba.atomic_num() != 1 {
            (ba.idx(), ea.isotope())
        } else {
            continue;
        };
        iso_map.entry(ha).or_default().push(iso);
    }
    iso_map
}

/// anonymous-namespace `may_need_extra_H`
fn may_need_extra_h(mol: &RWMol, atom: u32) -> Result<bool> {
    let mut single = 0;
    let mut aromatic = 0;
    for b in mol.atom_bonds(atom) {
        match b.bond_type() {
            BondType::Single => single += 1,
            BondType::Aromatic => aromatic += 1,
            _ => return Ok(false),
        }
    }
    Ok(single == 1 && aromatic == 2 && mol.atom(atom).total_valence()? == 3)
}

/// anonymous-namespace `adjustStereoAtomsIfRequired`
fn adjust_stereo_atoms_if_required(mol: &mut RWMol, atom: u32, heavy: u32) -> Result<bool> {
    if mol.atom_degree(heavy) == 2 {
        return Ok(false);
    }
    if mol.bond_idx_between(atom, heavy).is_none() {
        return Ok(false);
    }
    let bonds: Vec<u32> = mol.atom_bond_indices(heavy).collect();
    for bi in bonds {
        let bnd = mol.bond(bi);
        if bnd.bond_type() == BondType::Double && bnd.stereo() > BondStereo::StereoAny {
            let Some(pos) = bnd.stereo_atoms().iter().position(|&a| a == atom as i32) else {
                continue;
            };
            let dbl_nbr = bnd.other_atom_idx(heavy)?;
            let nbrs: Vec<u32> = mol.atom_neighbors(heavy).collect();
            for nbr in nbrs {
                if nbr == dbl_nbr || nbr == atom {
                    continue;
                }
                let bnd = mol.bond_mut(bi);
                bnd.stereo_atoms_mut()[pos] = nbr as i32;
                return Ok(match bnd.stereo() {
                    BondStereo::StereoCis => {
                        bnd.set_stereo(BondStereo::StereoTrans)?;
                        true
                    }
                    BondStereo::StereoTrans => {
                        bnd.set_stereo(BondStereo::StereoCis)?;
                        true
                    }
                    _ => false,
                });
            }
        }
    }
    Ok(false)
}

/// anonymous-namespace `molRemoveH`
fn mol_remove_h(mol: &mut RWMol, idx: u32, update_explicit_count: bool) -> Result<()> {
    if mol.atom(idx).atomic_num() != 1 {
        return Err(RDKitError::precondition(
            "idx corresponds to a non-Hydrogen",
        ));
    }
    let bonds: Vec<u32> = mol.atom_bond_indices(idx).collect();
    for bond in bonds {
        let heavy = mol.bond(bond).other_atom_idx(idx)?;
        let heavy_num = mol.atom(heavy).atomic_num();

        let h = mol.atom(heavy);
        if update_explicit_count || h.no_implicit() || h.chiral_tag() != ChiralType::Unspecified {
            let n = h.num_explicit_hs();
            mol.atom_mut(heavy).set_num_explicit_hs(n + 1);
        } else {
            let default_vs = PeriodicTable::get().get_valence_list(heavy_num)?;
            let tv = h.total_valence()? as i32;
            if ((heavy_num == 7 || heavy_num == 15 || may_need_extra_h(mol, heavy)?)
                && is_aromatic_atom(mol, heavy))
                || default_vs.iter().skip(1).any(|&v| v == tv)
            {
                let n = mol.atom(heavy).num_explicit_hs();
                mol.atom_mut(heavy).set_num_explicit_hs(n + 1);
            }
        }

        if mol.atom(heavy).chiral_tag() != ChiralType::Unspecified {
            let mut nbr_indices: Vec<i32> = mol
                .atom_bond_indices(heavy)
                .filter(|&b| b != bond)
                .map(|b| b as i32)
                .collect();
            nbr_indices.push(bond as i32);
            let n_swaps = mol.perturbation_order(heavy, &nbr_indices)?;
            if n_swaps % 2 != 0 {
                mol.atom_mut(heavy).invert_chirality();
            }
        }

        if mol.atom_degree(heavy) == 2 {
            let other = mol.atom_bond_indices(heavy).find(|&b| b != bond);
            if let Some(nb) = other {
                let nbnd = mol.bond_mut(nb);
                if nbnd.stereo() > BondStereo::StereoAny {
                    nbnd.set_stereo(BondStereo::StereoNone)?;
                    nbnd.stereo_atoms_mut().clear();
                }
            }
        }

        let dir = mol.bond(bond).bond_dir();
        if dir == BondDir::Unknown && mol.bond(bond).begin_atom_idx() == heavy {
            mol.atom_mut(heavy)
                .set_prop(cp::_UnknownStereo, 1i32, false)?;
        } else if matches!(dir, BondDir::EndDownRight | BondDir::EndUpRight) {
            let mut found_a_dir = false;
            let mut obond = None;
            for nb in mol.atom_bond_indices(heavy) {
                let nbnd = mol.bond(nb);
                if nb != bond && nbnd.bond_type() == BondType::Single {
                    if nbnd.bond_dir() == BondDir::None {
                        obond = Some(nb);
                    } else {
                        found_a_dir = true;
                    }
                }
            }
            if let (false, Some(ob)) = (found_a_dir, obond) {
                let flip_it = mol.bond(ob).begin_atom_idx() == heavy
                    && mol.bond(bond).begin_atom_idx() == heavy;
                let new_dir = if flip_it {
                    if dir == BondDir::EndDownRight {
                        BondDir::EndUpRight
                    } else {
                        BondDir::EndDownRight
                    }
                } else {
                    dir
                };
                mol.bond_mut(ob).set_bond_dir(new_dir);
            }
            adjust_stereo_atoms_if_required(mol, idx, heavy)?;
        } else {
            adjust_stereo_atoms_if_required(mol, idx, heavy)?;
        }
        // TODO(SubstanceGroups): sg.removeBondWithIdx(bond)
    }
    // TODO(SubstanceGroups): remove the atom from SGroups / attach points
    mol.remove_atom_ext(idx, false)
}

/// anonymous-namespace `shouldRemoveH`
fn should_remove_h(mol: &RWMol, atom: u32, ps: &RemoveHsParameters) -> Result<bool> {
    let a = mol.atom(atom);
    if a.atomic_num() != 1 {
        return Ok(false);
    }
    if !ps.remove_with_query && a.has_query() {
        return Ok(false);
    }
    let degree = mol.atom_degree(atom);
    if !ps.remove_degree_zero && degree == 0 {
        if ps.show_warnings {
            crate::rd_log_warning!("WARNING: not removing hydrogen atom without neighbors");
        }
        return Ok(false);
    }
    if !ps.remove_higher_degrees && degree > 1 {
        return Ok(false);
    }
    if !ps.remove_isotopes && !ps.remove_and_track_isotopes && a.isotope() != 0 {
        return Ok(false);
    }
    if !ps.remove_nonimplicit && !a.has_prop(cp::isImplicit) {
        return Ok(false);
    }
    if !ps.remove_mapped && a.atom_map_num() != 0 {
        return Ok(false);
    }
    // TODO(SubstanceGroups): removeInSGroups checks
    if !ps.remove_hydrides && a.formal_charge() == -1 {
        return Ok(false);
    }
    if degree > 0
        && (!ps.remove_dummy_neighbors
            || !ps.remove_defining_bond_stereo
            || !ps.remove_only_h_neighbors
            || !ps.remove_nontetrahedral_neighbors
            || !ps.remove_with_wedged_bond)
    {
        let mut only_h_neighbors = true;
        for nbr in mol.atom_neighbors(atom) {
            let n = mol.atom(nbr);
            if !ps.remove_dummy_neighbors && n.atomic_num() < 1 {
                if ps.show_warnings {
                    crate::rd_log_warning!(
                        "WARNING: not removing hydrogen atom with dummy atom neighbors"
                    );
                }
                return Ok(false);
            }
            if !ps.remove_nontetrahedral_neighbors && has_non_tetrahedral_stereo(mol, nbr) {
                if ps.show_warnings {
                    crate::rd_log_warning!(
                        "WARNING: not removing hydrogen atom with neighbor that has non-tetrahedral stereochemistry"
                    );
                }
                return Ok(false);
            }
            if !ps.remove_only_h_neighbors && n.atomic_num() != 1 {
                only_h_neighbors = false;
            }
            let bnd = mol.bond(mol.bond_idx_between(atom, nbr).expect("neighbor bond"));
            if !ps.remove_with_wedged_bond
                && matches!(bnd.bond_dir(), BondDir::BeginDash | BondDir::BeginWedge)
            {
                if ps.show_warnings {
                    crate::rd_log_warning!("WARNING: not removing hydrogen atom with wedged bond");
                }
                return Ok(false);
            }
            if !ps.remove_defining_bond_stereo && mol.atom_degree(nbr) == 2 {
                for nb in mol.atom_bonds(nbr) {
                    if nb.bond_type() == BondType::Double
                        && (nb.stereo() > BondStereo::StereoAny || bnd.bond_dir() > BondDir::None)
                    {
                        return Ok(false);
                    }
                }
            }
        }
        if !ps.remove_only_h_neighbors && only_h_neighbors {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The parameters `MolOps::removeAllHs` uses: remove every hydrogen.
pub fn remove_all_hs_parameters() -> RemoveHsParameters {
    RemoveHsParameters {
        remove_degree_zero: true,
        remove_higher_degrees: true,
        remove_only_h_neighbors: true,
        remove_isotopes: true,
        remove_dummy_neighbors: true,
        remove_defining_bond_stereo: true,
        remove_with_wedged_bond: true,
        remove_with_query: true,
        remove_nonimplicit: true,
        remove_in_sgroups: true,
        show_warnings: false,
        remove_hydrides: true,
        remove_nontetrahedral_neighbors: true,
        ..Default::default()
    }
}

/// `MolOps::removeAllHs(RWMol &mol, bool sanitize)`
pub fn remove_all_hs(mol: &mut RWMol, sanitize: bool) -> Result<()> {
    remove_hs(mol, &remove_all_hs_parameters(), sanitize)
}

/// `MolOps::removeHs(RWMol &mol, const RemoveHsParameters &ps, bool sanitize)`
pub fn remove_hs(mol: &mut RWMol, ps: &RemoveHsParameters, sanitize: bool) -> Result<()> {
    if ps.remove_and_track_isotopes {
        let need = mol
            .atoms()
            .iter()
            .any(|a| a.atomic_num() == 1 && a.isotope() == 0);
        if need {
            let mut copy = ps.clone();
            copy.remove_and_track_isotopes = false;
            copy.remove_isotopes = false;
            remove_hs(mol, &copy, false)?;
        }
    }
    for i in 0..mol.num_atoms() {
        mol.update_atom_property_cache(i, false)?;
    }
    if ps.remove_and_track_isotopes {
        for (atom, isos) in get_iso_map(mol) {
            mol.atom_mut(atom).set_prop(cp::_isotopicHs, isos, false)?;
        }
    }
    let n = mol.num_atoms();
    let mut to_remove = vec![false; n as usize];
    for i in 0..n {
        if should_remove_h(mol, i, ps)? {
            to_remove[i as usize] = true;
        }
    }
    // TODO(SubstanceGroups): filter_sgroup_emptying_hydrogens
    for idx in (0..n).rev() {
        if to_remove[idx as usize] {
            mol_remove_h(mol, idx, ps.update_explicit_count)?;
        }
    }
    mol.clear_computed_props(true);
    // note: upstream tests `!atomsToRemove.empty()`, i.e. that the molecule
    // had atoms, not that any were removed.
    if n > 0 && ps.remove_nonimplicit && sanitize {
        sanitize_mol(mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
    }
    if n > 0 {
        for i in 0..mol.num_atoms() {
            let a = mol.atom(i);
            if !a.no_implicit()
                && a.chiral_tag() != ChiralType::Unspecified
                && a.num_explicit_hs() > 1
            {
                mol.atom_mut(i).set_num_explicit_hs(0);
                mol.update_atom_property_cache(i, false)?;
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// addHs

use crate::geometry::point::Point3D;
use crate::geometry::transform3d::Transform3D;
use crate::graphmol::atom::{Atom, HybridizationType};

const SQ_DIST_ZERO_TOL: f64 = 1.0e-4;

/// `MolOps::AddHsParameters`
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AddHsParameters {
    pub explicit_only: bool,
    pub add_coords: bool,
    pub add_residue_info: bool,
    pub skip_queries: bool,
}

/// anonymous-namespace `getAtomNeighborNot`: first neighbor of `atom` that is not `other`
fn atom_neighbor_not(mol: &RWMol, atom: u32, other: u32) -> Result<u32> {
    mol.atom_neighbors(atom)
        .find(|&n| n != other)
        .ok_or_else(|| RDKitError::invariant("no neighbor found"))
}

/// anonymous-namespace `pickBisector`
fn pick_bisector(nbr1: &Point3D, nbr2: &Point3D, nbr3: &Point3D) -> Point3D {
    let mut dir = *nbr2 + *nbr3;
    if dir.length_sq() < SQ_DIST_ZERO_TOL {
        dir = *nbr2;
        std::mem::swap(&mut dir.x, &mut dir.y);
        dir.x *= -1.0;
    }
    if dir.dot_product(nbr1) < 0.0 {
        dir *= -1.0;
    }
    dir
}

/// `MolOps::setTerminalAtomCoords`: places terminal atom `idx` bonded to `other_idx`.
pub fn set_terminal_atom_coords(mol: &mut RWMol, idx: u32, other_idx: u32) -> Result<()> {
    const FILE: &str = "Code/GraphMol/AddHs.cpp";
    if other_idx == idx {
        return Err(RDKitError::precondition_at(
            "degenerate atoms",
            "otherIdx != idx",
            FILE,
            183,
        ));
    }
    if idx >= mol.num_atoms() {
        return Err(RDKitError::range_at(
            "idx",
            idx,
            mol.num_atoms(),
            "Code/GraphMol/ROMol.cpp",
            202,
        ));
    }
    if mol.atom_degree(idx) != 1 {
        return Err(RDKitError::precondition_at(
            "bad atom degree",
            "mol.getAtomDegree(atom) == 1",
            FILE,
            185,
        ));
    }
    if other_idx >= mol.num_atoms() {
        return Err(RDKitError::range_at(
            "idx1",
            other_idx,
            mol.num_atoms(),
            "Code/GraphMol/ROMol.cpp",
            344,
        ));
    }
    if mol.get_bond_between_atoms(other_idx, idx)?.is_none() {
        return Err(RDKitError::precondition_at(
            "no bond between atoms",
            "bond",
            FILE,
            187,
        ));
    }
    let pt = PeriodicTable::get();
    let bond_length = pt.get_rb0(1)? + pt.get_rb0(mol.atom(other_idx).atomic_num())?;
    let other_nbrs: Vec<u32> = mol.atom_neighbors(other_idx).collect();
    let hyb = mol.atom(other_idx).hybridization();
    // upstream keeps these across the conformer loop
    let mut dir = Point3D::default();
    let mut perp = Point3D::default();
    let mut other_pos = Point3D::default();
    let scale = |is3d: bool| if is3d { bond_length } else { 1.0 };
    match other_nbrs.len() {
        1 => {
            for conf in mol.conformers_mut() {
                if conf.is_3d() {
                    dir.z = 1.0;
                } else {
                    dir.x = 1.0;
                }
                other_pos = conf.atom_pos(other_idx)?;
                conf.set_atom_pos(idx, other_pos + dir * scale(conf.is_3d()))?;
            }
        }
        2 => {
            let nbr1 = atom_neighbor_not(mol, other_idx, idx)?;
            let nbr2 = if mol.atom_neighbors(nbr1).count() > 1 {
                let nb = mol
                    .get_bond_between_atoms(other_idx, nbr1)?
                    .expect("neighbors are bonded");
                if nb.is_aromatic() || nb.bond_type() == BondType::Double || nb.is_conjugated() {
                    Some(atom_neighbor_not(mol, nbr1, other_idx)?)
                } else {
                    None
                }
            } else {
                None
            };
            for conf in mol.conformers_mut() {
                other_pos = conf.atom_pos(other_idx)?;
                let nbr1_pos = conf.atom_pos(nbr1)?;
                let mut nbr1_vect = nbr1_pos - other_pos;
                if nbr1_vect.length_sq() < SQ_DIST_ZERO_TOL {
                    conf.set_atom_pos(idx, other_pos)?;
                    continue;
                }
                nbr1_vect.normalize()?;
                nbr1_vect *= -1.0;
                let is3d = conf.is_3d();
                match hyb {
                    HybridizationType::SP3 => {
                        if is3d {
                            perp = nbr1_vect.get_perpendicular()?;
                        } else {
                            perp.z = 1.0;
                        }
                        let mut tf = Transform3D::identity();
                        tf.set_rotation((180.0 - 109.471) * std::f64::consts::PI / 180.0, &perp);
                        dir = tf.transform_point(&nbr1_vect);
                    }
                    HybridizationType::SP2 => {
                        if is3d {
                            perp = nbr1_vect.get_perpendicular()?;
                        } else {
                            perp.z = 1.0;
                        }
                        if let Some(n2) = nbr2 {
                            let nbr2_vect = nbr1_pos.direction_vector(&conf.atom_pos(n2)?)?;
                            let cross = nbr2_vect.cross_product(&nbr1_vect);
                            if cross.length_sq() >= SQ_DIST_ZERO_TOL {
                                perp = cross;
                            }
                        }
                        perp.normalize()?;
                        let mut tf = Transform3D::identity();
                        tf.set_rotation(60.0 * std::f64::consts::PI / 180.0, &perp);
                        dir = tf.transform_point(&nbr1_vect);
                    }
                    _ => dir = nbr1_vect,
                }
                conf.set_atom_pos(idx, other_pos + dir * scale(is3d))?;
            }
        }
        3 => {
            let others: Vec<u32> = other_nbrs.iter().copied().filter(|&n| n != idx).collect();
            let (nbr1, nbr2) = (others[0], others[1]);
            for conf in mol.conformers_mut() {
                other_pos = conf.atom_pos(other_idx)?;
                let mut v1 = other_pos - conf.atom_pos(nbr1)?;
                let mut v2 = other_pos - conf.atom_pos(nbr2)?;
                if v1.length_sq() < SQ_DIST_ZERO_TOL || v2.length_sq() < SQ_DIST_ZERO_TOL {
                    conf.set_atom_pos(idx, other_pos)?;
                    continue;
                }
                v1.normalize()?;
                v2.normalize()?;
                dir = v1 + v2;
                if dir.length_sq() < SQ_DIST_ZERO_TOL {
                    continue;
                }
                dir.normalize()?;
                if conf.is_3d() {
                    if hyb == HybridizationType::SP3 {
                        let nbr_perp = v1.cross_product(&v2);
                        let mut axis = nbr_perp.cross_product(&dir);
                        axis.normalize()?;
                        let mut tf = Transform3D::identity();
                        tf.set_rotation((109.471 / 2.0) * std::f64::consts::PI / 180.0, &axis);
                        dir = tf.transform_point(&dir);
                    }
                    conf.set_atom_pos(idx, other_pos + dir * bond_length)?;
                } else {
                    conf.set_atom_pos(idx, other_pos + dir)?;
                }
            }
        }
        4 => {
            let others: Vec<u32> = other_nbrs.iter().copied().filter(|&n| n != idx).collect();
            let (nbr1, nbr2, nbr3) = (others[0], others[1], others[2]);
            let other_atom = mol.atom(other_idx);
            let has_cip = other_atom.get_prop(cp::_CIPCode).is_some();
            let chiral = other_atom.chiral_tag();
            for conf in mol.conformers_mut() {
                other_pos = conf.atom_pos(other_idx)?;
                let mut v1 = other_pos - conf.atom_pos(nbr1)?;
                let mut v2 = other_pos - conf.atom_pos(nbr2)?;
                let mut v3 = other_pos - conf.atom_pos(nbr3)?;
                if v1.length_sq() < SQ_DIST_ZERO_TOL
                    || v2.length_sq() < SQ_DIST_ZERO_TOL
                    || v3.length_sq() < SQ_DIST_ZERO_TOL
                {
                    conf.set_atom_pos(idx, other_pos)?;
                    continue;
                }
                v1.normalize()?;
                v2.normalize()?;
                v3.normalize()?;
                if conf.is_3d() {
                    if v3.dot_product(&v1.cross_product(&v2)).abs() < 0.1 {
                        dir = v1.cross_product(&v2);
                        if dir.length_sq() < SQ_DIST_ZERO_TOL {
                            dir = v1.cross_product(&v3) * -1.0;
                        }
                        if dir.length_sq() < SQ_DIST_ZERO_TOL {
                            dir = v2.cross_product(&v3);
                        }
                        if dir.length_sq() < SQ_DIST_ZERO_TOL {
                            continue;
                        }
                        if has_cip {
                            let a = dir - v3;
                            let b = v1 - v3;
                            let c = v2 - v3;
                            let vol = a.dot_product(&b.cross_product(&c));
                            if (chiral == ChiralType::TetrahedralCCW && vol < 0.0)
                                || (chiral == ChiralType::TetrahedralCW && vol > 0.0)
                            {
                                dir *= -1.0;
                            }
                        }
                    } else {
                        dir = v1 + v2 + v3;
                    }
                } else {
                    let a12 = v1.angle_to(&v2);
                    let a13 = v1.angle_to(&v3);
                    let a23 = v2.angle_to(&v3);
                    let (s1, s2, s3) = (a12 + a13, a12 + a23, a13 + a23);
                    dir = if s1 <= s2 && s1 <= s3 {
                        pick_bisector(&v1, &v2, &v3)
                    } else if s2 <= s1 && s2 <= s3 {
                        pick_bisector(&v2, &v1, &v3)
                    } else {
                        pick_bisector(&v3, &v1, &v2)
                    };
                }
                dir.normalize()?;
                let s = scale(conf.is_3d());
                conf.set_atom_pos(idx, other_pos + dir * s)?;
            }
        }
        _ => {
            let p = other_pos + dir * bond_length;
            for conf in mol.conformers_mut() {
                conf.set_atom_pos(idx, p)?;
            }
        }
    }
    Ok(())
}

/// anonymous-namespace `isQueryAtom`
fn is_query_atom(mol: &RWMol, idx: u32) -> bool {
    mol.atom(idx).has_query() || mol.atom_bonds(idx).any(|b| b.has_query())
}

/// `MolOps::addHs(RWMol&, params, onlyOnAtoms)`
pub fn add_hs(
    mol: &mut RWMol,
    params: &AddHsParameters,
    only_on_atoms: Option<&[u32]>,
) -> Result<()> {
    mol.clear_computed_props(false);
    let n = mol.num_atoms();
    let mut on_atoms = vec![only_on_atoms.is_none(); n as usize];
    if let Some(oa) = only_on_atoms {
        for &a in oa {
            *on_atoms
                .get_mut(a as usize)
                .ok_or_else(|| RDKitError::Index(format!("Index Error: {a}")))? = true;
        }
    }
    let mut n_expl = vec![0u32; n as usize];
    let mut n_impl = vec![0u32; n as usize];
    for i in 0..n {
        n_expl[i as usize] = mol.atom(i).num_explicit_hs();
        n_impl[i as usize] = mol.atom(i).num_implicit_hs()?;
        if on_atoms[i as usize] && params.skip_queries && is_query_atom(mol, i) {
            on_atoms[i as usize] = false;
        }
    }
    for aidx in 0..n {
        if !on_atoms[aidx as usize] {
            continue;
        }
        let iso_hs: Vec<u32> = match mol.atom(aidx).get_prop(cp::_isotopicHs) {
            Some(v) => {
                let v = match v {
                    crate::rdgeneral::RDValue::VecUnsignedInt(v) => v.clone(),
                    _ => Vec::new(),
                };
                mol.atom_mut(aidx).clear_prop(cp::_isotopicHs);
                v
            }
            None => Vec::new(),
        };
        let mut iso = iso_hs.iter();
        mol.atom_mut(aidx).clear_computed_props();
        let n_h = n_expl[aidx as usize]
            + if params.explicit_only {
                0
            } else {
                n_impl[aidx as usize]
            };
        for i in 0..n_h {
            let new_idx = mol.add_atom(Atom::new(1), false);
            mol.add_bond(aidx, new_idx, BondType::Single)?;
            if i >= n_expl[aidx as usize] {
                mol.atom_mut(new_idx)
                    .set_prop(cp::isImplicit, 1i32, false)?;
            }
            mol.update_atom_property_cache(new_idx, true)?;
            if params.add_coords {
                set_terminal_atom_coords(mol, new_idx, aidx)?;
            }
            if let Some(&iv) = iso.next() {
                mol.atom_mut(new_idx).set_isotope(iv);
            }
        }
        mol.atom_mut(aidx).set_num_explicit_hs(0);
        mol.update_atom_property_cache(aidx, false)?;
        if iso.next().is_some() {
            crate::rd_log_warning!("extra H isotope information found on atom {aidx}");
        }
    }
    Ok(())
}

/// `MolOps::needsHs`: whether any atom has implicit hydrogens.
pub fn needs_hs(mol: &crate::graphmol::romol::ROMol) -> Result<bool> {
    for idx in 0..mol.num_atoms() {
        if mol.total_num_hs(idx, false)? > 0 {
            return Ok(true);
        }
    }
    Ok(false)
}
