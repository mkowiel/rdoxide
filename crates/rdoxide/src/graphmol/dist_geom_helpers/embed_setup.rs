//! Preparing an embedding: bounds matrix, ETKDG details, chiral sets and
//! double bonds (upstream `Embedder.cpp`, `EmbeddingOps`).

use std::collections::BTreeMap;

use super::bounds_builder::{
    TopolBoundsOptions, TopolBoundsOutput, init_bounds_mat, set_topol_bounds,
    set_topol_bounds_with_bonds_and_angles,
};
use super::embed_params::EmbedParameters;
use crate::distgeom::{BoundsMatrix, ChiralSet, ChiralSetStructureFlags, triangle_smooth_bounds};
use crate::error::{RDKitError, Result};
use crate::geometry::Point3D;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::atropisomers::get_atropisomer_atoms_and_bonds;
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::forcefield_helpers::crystal_ff::{
    CrystalFFDetails, EtkdgForceConsts, TorsionOptions, get_experimental_torsions,
    populate_ref_table,
};
use crate::graphmol::romol::ROMol;
use crate::rd_log_warning;

/// The atoms whose positions are given: map from atom index to position.
pub type CoordMap = BTreeMap<i32, Point3D>;

/// A double bond end (neighbor, double-bond atom, other double-bond atom) and
/// a stereo double bond (the four controlling atoms; -1 for cis, 1 for trans).
pub(super) type DoubleBondEnd = (u32, u32, u32);
pub(super) type StereoDoubleBond = ([u32; 4], i32);

/// Collects the double bonds that must not become linear and the stereo double
/// bonds with their expected cis/trans arrangement.
pub(super) fn find_double_bonds(
    mol: &ROMol,
    coord_map: Option<&CoordMap>,
) -> Result<(Vec<DoubleBondEnd>, Vec<StereoDoubleBond>)> {
    let mut ends = Vec::new();
    let mut stereo = Vec::new();
    for bnd in mol.bonds() {
        if bnd.bond_type() != BondType::Double {
            continue;
        }
        for atm in [bnd.begin_atom_idx(), bnd.end_atom_idx()] {
            if mol.atom_degree(atm) < 2 {
                continue;
            }
            let oatm = bnd.other_atom_idx(atm)?;
            for nbr in mol.atom_neighbors(atm) {
                if nbr == oatm {
                    continue;
                }
                let Some(obnd) = mol.bond_idx_between(atm, nbr) else {
                    continue;
                };
                if mol.bond(obnd).bond_type() != BondType::Single && mol.atom_degree(atm) == 2 {
                    continue;
                }
                ends.push((nbr, atm, oatm));
            }
        }
        // stereo
        if bnd.stereo() > BondStereo::StereoAny {
            let atoms = bnd.stereo_atoms();
            if atoms.len() < 2 {
                continue;
            }
            // only if the controlling atoms are not in the coord map
            if coord_map.is_some_and(|m| m.contains_key(&atoms[0]) && m.contains_key(&atoms[1])) {
                continue;
            }
            let sign = if matches!(bnd.stereo(), BondStereo::StereoCis | BondStereo::StereoZ) {
                -1
            } else {
                1
            };
            stereo.push((
                [
                    atoms[0] as u32,
                    bnd.begin_atom_idx(),
                    bnd.end_atom_idx(),
                    atoms[1] as u32,
                ],
                sign,
            ));
        }
    }
    Ok((ends, stereo))
}

/// Sets the bounds between the atoms of the coordinate map to their distances.
pub(super) fn adjust_bounds_mat_from_coord_map(
    mmat: &mut BoundsMatrix,
    coord_map: &CoordMap,
) -> Result<()> {
    let entries: Vec<(i32, Point3D)> = coord_map.iter().map(|(&k, &v)| (k, v)).collect();
    for (i, (i_idx, i_pt)) in entries.iter().enumerate() {
        for (j_idx, j_pt) in &entries[i + 1..] {
            let dist = (*i_pt - *j_pt).length();
            mmat.set_upper_bound(*i_idx as usize, *j_idx as usize, dist)?;
            mmat.set_lower_bound(*i_idx as usize, *j_idx as usize, dist)?;
        }
    }
    Ok(())
}

/// Sets the force constants and collects the torsion terms of ETKDG.
pub(super) fn init_etkdg(
    mol: &ROMol,
    params: &EmbedParameters,
    details: &mut CrystalFFDetails,
) -> Result<()> {
    let n_atoms = mol.num_atoms();
    details.force_consts = match (params.et_version == 4, params.use_legacy_implementation) {
        (true, true) => EtkdgForceConsts::SEQ_GAUSSIAN,
        (true, false) => EtkdgForceConsts::AIO_GAUSSIAN,
        (false, true) => EtkdgForceConsts::SEQ_COSINE,
        (false, false) => EtkdgForceConsts::AIO_COSINE,
    };
    if params.use_exp_torsion_angle_prefs || params.use_basic_knowledge {
        get_experimental_torsions(
            mol,
            details,
            &TorsionOptions {
                use_exp_torsions: params.use_exp_torsion_angle_prefs,
                use_small_ring_torsions: params.use_small_ring_torsions,
                use_macrocycle_torsions: params.use_macrocycle_torsions,
                use_basic_knowledge: params.use_basic_knowledge,
                version: params.et_version,
                verbose: params.verbose,
            },
        )?;
        details.atom_nums = (0..n_atoms)
            .map(|i| mol.atom(i).atomic_num() as i32)
            .collect();
    }
    if params.et_version == 4 {
        populate_ref_table(details)?;
    }
    details.bounds_mat_force_scaling = params.bounds_mat_force_scaling;
    Ok(())
}

/// Builds and smooths the bounds matrix from the topology; false if it cannot
/// be smoothed.
pub(super) fn setup_initial_bounds_matrix(
    mol: &ROMol,
    mmat: &mut BoundsMatrix,
    coord_map: Option<&CoordMap>,
    params: &EmbedParameters,
    details: &mut CrystalFFDetails,
) -> Result<bool> {
    let options = TopolBoundsOptions {
        scale_vdw: false,
        set15_bounds: true,
        set14_bounds: true,
        set13_bounds: true,
        embed_force_field: params.embed_force_field,
    };
    if params.use_exp_torsion_angle_prefs || params.use_basic_knowledge {
        let mut paths = Vec::new();
        let (bonds, angles) = set_topol_bounds_with_bonds_and_angles(
            mol,
            mmat,
            params,
            &options,
            TopolBoundsOutput {
                paths14: Some(&mut paths),
                internal_coords: details.internal_coords.as_deref_mut(),
            },
        )?;
        details.bonds = bonds;
        details.angles = angles;
        details.path14_configs = paths;
    } else {
        set_topol_bounds(
            mol,
            mmat,
            params,
            &options,
            TopolBoundsOutput {
                paths14: None,
                internal_coords: details.internal_coords.as_deref_mut(),
            },
        )?;
    }
    let mut tol = 0.0;
    if let Some(map) = coord_map {
        adjust_bounds_mat_from_coord_map(mmat, map)?;
        tol = 0.05;
    }
    if triangle_smooth_bounds(mmat, tol) {
        return Ok(true);
    }
    // the bounds failed to smooth: recompute them without 1-5 bounds and with
    // van der Waals scaling
    let relaxed = TopolBoundsOptions {
        scale_vdw: true,
        set15_bounds: false,
        ..TopolBoundsOptions::default()
    };
    let recompute = |mmat: &mut BoundsMatrix| -> Result<()> {
        init_bounds_mat(mmat, 0.0, 1000.0)?;
        set_topol_bounds(mol, mmat, params, &relaxed, TopolBoundsOutput::default())?;
        if let Some(map) = coord_map {
            adjust_bounds_mat_from_coord_map(mmat, map)?;
        }
        Ok(())
    };
    recompute(mmat)?;
    if triangle_smooth_bounds(mmat, tol) {
        return Ok(true);
    }
    if params.ignore_smoothing_failures {
        // proceed with the more relaxed bounds matrix
        recompute(mmat)?;
        return Ok(true);
    }
    rd_log_warning!("Could not triangle bounds smooth molecule.");
    Ok(false)
}

/// The chiral sets (centers with a specified chirality, and atropisomeric
/// bonds) and the tetrahedral centers whose geometry is checked.
pub(super) fn find_chiral_sets(
    mol: &ROMol,
    coord_map: Option<&CoordMap>,
) -> Result<(Vec<ChiralSet>, Vec<ChiralSet>)> {
    let mut chiral_centers = Vec::new();
    let mut tetrahedral_centers = Vec::new();
    let rinfo = mol.ring_info();
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom(idx);
        if atom.atomic_num() == 1 {
            continue;
        }
        let chiral_type = atom.chiral_tag();
        let is_tet = matches!(
            chiral_type,
            ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
        );
        if !(is_tet
            || ((atom.atomic_num() == 6 || atom.atomic_num() == 7) && mol.atom_degree(idx) == 4))
        {
            continue;
        }
        // the chiral set is made from the neighbors
        let mut nbrs: Vec<usize> = mol
            .atom_bond_indices(idx)
            .map(|b| mol.bond(b).other_atom_idx(idx).map(|a| a as usize))
            .collect::<Result<_>>()?;
        if nbrs.len() < 3 {
            return Err(RDKitError::invariant("Cannot be a chiral center"));
        }
        let mut vol_lower_bound = 5.0;
        let vol_upper_bound = 100.0;
        if nbrs.len() < 4 {
            // lower volumes with three neighbors (github #5883): the center
            // takes the place of the missing one
            vol_lower_bound = 2.0;
            nbrs.push(idx as usize);
        }
        // tetrahedral centers in multiple small rings
        const SMALL_RING_SIZE: i32 = 5;
        let num_small_rings = rinfo
            .atom_ring_sizes(idx)?
            .iter()
            .filter(|&&sz| sz < SMALL_RING_SIZE)
            .count();
        let structure_flags = if num_small_rings > 1 {
            ChiralSetStructureFlags::InFusedSmallRings as u64
        } else {
            0
        };
        let pids = [idx as usize, nbrs[0], nbrs[1], nbrs[2], nbrs[3]];
        match chiral_type {
            ChiralType::TetrahedralCCW => {
                // positive chiral volume
                chiral_centers.push(ChiralSet::new(
                    pids,
                    vol_lower_bound,
                    vol_upper_bound,
                    structure_flags,
                )?);
            }
            ChiralType::TetrahedralCW => {
                chiral_centers.push(ChiralSet::new(
                    pids,
                    -vol_upper_bound,
                    -vol_lower_bound,
                    structure_flags,
                )?);
            }
            _ => {
                let in_coord_map = coord_map.is_some_and(|m| m.contains_key(&(idx as i32)));
                // these tests are only for ring atoms that are not in the
                // coordinate map; 3-rings are a nightmare
                let skip = in_coord_map
                    || (rinfo.is_initialized()
                        && (rinfo.num_atom_rings(idx)? < 2
                            || rinfo.is_atom_in_ring_of_size(idx, 3)?));
                if !skip {
                    tetrahedral_centers.push(ChiralSet::new(pids, 0.0, 0.0, structure_flags)?);
                }
            }
        }
    }

    // atropisomers
    for bond in mol.bonds() {
        if !matches!(
            bond.stereo(),
            BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW
        ) {
            continue;
        }
        let sides = get_atropisomer_atoms_and_bonds(mol, bond.idx());
        // we only manage cases with two exo-substituents on at least one side
        let Some(sides) = sides.filter(|s| s[0].1.len() == 2 || s[1].1.len() == 2) else {
            rd_log_warning!(
                "Atropisomer bond stereochemistry not used for bond {}, which does not have two exo substituents on at least one side.",
                bond.idx()
            );
            continue;
        };
        let idx0 = sides[0].0;
        let idx1 = sides[1].0;
        let other = |b: u32, a: u32| mol.bond(b).other_atom_idx(a);
        let nbr1 = other(sides[0].1[0], idx0)?;
        let (nbr2, nbr3, nbr4);
        if sides[0].1.len() == 2 {
            nbr2 = other(sides[0].1[1], idx0)?;
            nbr3 = other(sides[1].1[0], idx1)?;
            nbr4 = if sides[1].1.len() == 2 {
                other(sides[1].1[1], idx1)?
            } else {
                idx0
            };
        } else {
            nbr2 = other(sides[1].1[0], idx1)?;
            nbr3 = other(sides[1].1[1], idx1)?;
            nbr4 = idx0;
        }
        // FIX (upstream): these numbers are empirical
        let (mut lower, mut upper) = (1.0, 100.0);
        if bond.stereo() == BondStereo::StereoAtropCCW {
            std::mem::swap(&mut lower, &mut upper);
            lower *= -1.0;
            upper *= -1.0;
        }
        chiral_centers.push(ChiralSet::new(
            [
                idx0 as usize,
                nbr1 as usize,
                nbr2 as usize,
                nbr3 as usize,
                nbr4 as usize,
            ],
            lower,
            upper,
            0,
        )?);
    }
    Ok((chiral_centers, tetrahedral_centers))
}
