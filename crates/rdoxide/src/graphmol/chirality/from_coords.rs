//! Stereochemistry from coordinates: wedged bonds in 2D
//! (`MolOps::assignChiralTypesFromBondDirs`,
//! `Chirality::atomChiralTypeFromBondDirPseudo3D`), 3D geometry
//! (`MolOps::assignChiralTypesFrom3D`) and double bonds
//! (`MolOps::detectBondStereochemistry`).
//!
//! Ported from RDKit `Code/GraphMol/Chirality.cpp` (2026.09.1).
//!
//! Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::f64::consts::PI;

use super::bond_dirs::set_double_bond_neighbor_directions;
use super::get_allow_nontetrahedral_chirality;
use super::legacy::{bond_affects_atom_chirality, get_atom_nonzero_degree};
use crate::geometry::Point3D;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::atropisomers::detect_atropisomer_chirality;
use crate::graphmol::bond::{BondDir, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result, rd_log_warning};

fn conf_position(mol: &ROMol, conf_id: i32) -> Result<usize> {
    let c = mol.conformer(conf_id)?;
    let id = c.id();
    Ok(mol
        .conformers()
        .iter()
        .position(|x| x.id() == id)
        .unwrap_or(0))
}

/// `Chirality::atomChiralTypeFromBondDirPseudo3D`. `Ok(None)` is upstream's
/// `std::nullopt` (an ambiguous center, warned about).
pub fn atom_chiral_type_from_bond_dir_pseudo3d(
    mol: &ROMol,
    bond: u32,
    conf: &Conformer,
) -> Result<Option<ChiralType>> {
    const COORD_ZERO_TOL: f64 = 1e-4;
    const ZERO_TOL: f64 = 1e-3;
    const T_SHAPE_TOL: f64 = 0.00031;
    const PSEUDO_3D_OFFSET: f64 = 0.1;
    const VOLUME_TOLERANCE: f64 = 0.00174;

    let b = mol.bond(bond);
    let bond_dir = b.bond_dir();
    if !matches!(bond_dir, BondDir::BeginWedge | BondDir::BeginDash) {
        return Err(RDKitError::precondition("bad bond direction"));
    }
    let atom = b.begin_atom_idx();
    if mol.atom_degree(atom) > 4 {
        return Ok(Some(ChiralType::Unspecified));
    }
    let bond_atom = b.end_atom_idx();
    let mut res = ChiralType::Unspecified;

    let mut center_loc = conf.atom_pos(atom)?;
    center_loc.z = 0.0;
    let mut ref_pt = conf.atom_pos(bond_atom)?;
    let ref_length = (center_loc - ref_pt).length();
    ref_pt.z = if bond_dir == BondDir::BeginWedge {
        PSEUDO_3D_OFFSET
    } else {
        -PSEUDO_3D_OFFSET
    };
    if ref_length != 0.0 {
        ref_pt.z *= ref_length;
    }

    let mut ref_idx = mol.num_bonds() as usize + 1;
    let mut bond_vects: Vec<Point3D> = Vec::new();
    let mut all_single = true;
    for (nbr_idx, nb) in mol.atom_bond_indices(atom).enumerate() {
        let nbr_bond = mol.bond(nb);
        let o_atom = nbr_bond.other_atom_idx(atom)?;
        let mut tmp = conf.atom_pos(o_atom)?;
        if nb == bond {
            ref_idx = nbr_idx;
            tmp = ref_pt;
        } else {
            if nbr_bond.begin_atom_idx() == atom
                && matches!(
                    nbr_bond.bond_dir(),
                    BondDir::BeginWedge | BondDir::BeginDash
                )
            {
                tmp.z = if nbr_bond.bond_dir() == BondDir::BeginWedge {
                    PSEUDO_3D_OFFSET
                } else {
                    -PSEUDO_3D_OFFSET
                };
                if ref_length != 0.0 {
                    tmp.z *= ref_length;
                }
            } else {
                tmp.z = 0.0;
            }
            if (center_loc - tmp).length_sq() < ZERO_TOL {
                rd_log_warning!(
                    "Warning: ambiguous stereochemistry - zero-length (or near zero-length) bond - at atom {atom} ignored."
                );
                return Ok(None);
            }
        }
        if nbr_bond.bond_type() != BondType::Single {
            all_single = false;
        }
        bond_vects.push(center_loc.direction_vector(&tmp)?);
    }
    if ref_idx >= mol.num_bonds() as usize {
        return Err(RDKitError::invariant(
            "could not find reference bond in neighbors",
        ));
    }
    let n_nbrs = bond_vects.len();
    if !(3..=4).contains(&n_nbrs) {
        return Ok(None);
    }
    for i in 0..n_nbrs {
        for j in 0..i {
            if (bond_vects[i] - bond_vects[j]).length_sq() < ZERO_TOL {
                rd_log_warning!(
                    "Warning: ambiguous stereochemistry - overlapping neighbors  - at atom {atom} ignored"
                );
                return Ok(None);
            }
        }
    }

    let anum = mol.atom(atom).atomic_num();
    if all_single || anum == 15 || anum == 16 {
        let mut order = [0usize, 1, 2, 3];
        let mut prefactor = 1.0;
        if ref_idx != 0 {
            order.swap(0, ref_idx);
            prefactor *= -1.0;
        }
        if n_nbrs > 3
            && bond_vects[order[1]]
                .cross_product(&bond_vects[order[2]])
                .length_sq()
                < 10.0 * ZERO_TOL
            && bond_vects[order[1]]
                .cross_product(&bond_vects[order[0]])
                .length_sq()
                > 10.0 * ZERO_TOL
        {
            bond_vects[order[1]].z = -bond_vects[order[0]].z;
        }

        let needs_swap = |cp01: &Point3D, cp02: &Point3D, dp01: f64, dp02: f64| -> bool {
            if dp01.abs() - 1.0 > -ZERO_TOL {
                return cp02.z < 0.0;
            }
            if dp02.abs() - 1.0 > -ZERO_TOL && cp01.z < 0.0 {
                return true;
            }
            if cp01.z * cp02.z < -ZERO_TOL {
                return cp01.z < cp02.z;
            }
            if dp01 * dp02 < -ZERO_TOL {
                return dp01 < dp02;
            }
            dp01.abs() > dp02.abs()
        };

        if n_nbrs == 3 {
            let cp01 = bond_vects[order[0]].cross_product(&bond_vects[order[1]]);
            let cp02 = bond_vects[order[0]].cross_product(&bond_vects[order[2]]);
            let dp01 = bond_vects[order[0]].dot_product(&bond_vects[order[1]]);
            let dp02 = bond_vects[order[0]].dot_product(&bond_vects[order[2]]);
            if needs_swap(&cp01, &cp02, dp01, dp02) {
                order.swap(1, 2);
                prefactor *= -1.0;
            }
        } else {
            let mut ordered: Vec<(f64, f64, usize)> = (1..4)
                .map(|i| {
                    let cp0i = bond_vects[order[0]].cross_product(&bond_vects[order[i]]);
                    let sgn = if cp0i.z < -ZERO_TOL { -1.0 } else { 1.0 };
                    let dp0i = bond_vects[order[0]].dot_product(&bond_vects[order[i]]);
                    (sgn, sgn * dp0i, order[i])
                })
                .collect();
            // std::sort(rbegin, rend): descending lexicographic order
            ordered.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            let mut n_changed = 0;
            for i in 1..4 {
                let ni = ordered[i - 1].2;
                if order[i] != ni {
                    order[i] = ni;
                    n_changed += 1;
                }
            }
            if n_changed == 2 {
                prefactor *= -1.0;
            }
        }

        for i in 0..n_nbrs {
            for j in (i + 1)..n_nbrs {
                if bond_vects[order[i]].z * bond_vects[order[j]].z < -ZERO_TOL {
                    let cp = bond_vects[order[i]]
                        .cross_product(&bond_vects[order[j]])
                        .length_sq();
                    if cp < 0.01 {
                        if n_nbrs == 4
                            && (bond_vects[order[i]].dot_product(&bond_vects[order[j]]) + 1.0).abs()
                                < ZERO_TOL
                            && (j - i == 1 || (i == 0 && j == 3))
                        {
                            bond_vects[order[j]].z = 0.0;
                            continue;
                        }
                        rd_log_warning!(
                            "Warning: ambiguous stereochemistry - opposing bonds have opposite wedging - at atom {atom} ignored."
                        );
                        return Ok(None);
                    }
                }
            }
        }

        if n_nbrs == 3 {
            let bv = &bond_vects;
            let mut conflict = false;
            if bv[order[1]].z * bv[order[0]].z < -COORD_ZERO_TOL
                && bv[order[2]].z.abs() < COORD_ZERO_TOL
            {
                conflict = bv[order[2]].cross_product(&bv[order[0]]).z
                    * bv[order[2]].cross_product(&bv[order[1]]).z
                    < -1e-4;
            } else if bv[order[2]].z * bv[order[0]].z < -COORD_ZERO_TOL
                && bv[order[1]].z.abs() < COORD_ZERO_TOL
            {
                conflict = bv[order[1]].cross_product(&bv[order[0]]).z
                    * bv[order[1]].cross_product(&bv[order[2]]).z
                    < -COORD_ZERO_TOL;
            }
            if conflict {
                rd_log_warning!(
                    "Warning: conflicting stereochemistry - bond wedging contradiction - at atom {atom} ignored"
                );
                return Ok(None);
            }
        }
        let mut bv1 = bond_vects[order[1]];
        bv1.z = 0.0;
        let mut bv2 = bond_vects[order[2]];
        bv2.z = 0.0;
        let mut crossp1 = bv1.cross_product(&bv2);
        if n_nbrs == 3 {
            if crossp1.length_sq() < T_SHAPE_TOL {
                bv1.z = -bond_vects[order[0]].z;
                bv2.z = -bond_vects[order[0]].z;
                crossp1 = bv1.cross_product(&bv2);
            }
        } else if crossp1.length_sq() < 10.0 * ZERO_TOL
            && bond_vects[order[3]].z.abs() < COORD_ZERO_TOL
        {
            bond_vects[order[3]].z = -bond_vects[order[0]].z;
        }
        let mut vol = crossp1.dot_product(&bond_vects[order[0]]);
        if n_nbrs == 4 {
            let dotp1 = bond_vects[order[1]].dot_product(&bond_vects[order[2]]);
            let mut bv3 = bond_vects[order[3]];
            bv3.z = 0.0;
            let crossp2 = bv1.cross_product(&bv3);
            let dotp2 = bond_vects[order[1]].dot_product(&bond_vects[order[3]]);
            let vol2 = crossp2.dot_product(&bond_vects[order[0]]);
            if vol.abs() < ZERO_TOL {
                if vol2.abs() < ZERO_TOL {
                    rd_log_warning!(
                        "Warning: ambiguous stereochemistry - no chiral volume - at atom {atom} ignored"
                    );
                    return Ok(None);
                }
                vol = vol2;
                prefactor *= -1.0;
            } else if vol * vol2 > 0.0 && vol2.abs() > VOLUME_TOLERANCE && dotp1 < dotp2 {
                vol = vol2;
                prefactor *= -1.0;
            } else if vol.abs() < VOLUME_TOLERANCE && vol2.abs() > VOLUME_TOLERANCE {
                if vol * vol2 < 0.0 {
                    prefactor *= -1.0;
                }
                vol = vol2;
            }
        }
        vol *= prefactor;
        if vol > VOLUME_TOLERANCE {
            res = ChiralType::TetrahedralCCW;
        } else if vol < -VOLUME_TOLERANCE {
            res = ChiralType::TetrahedralCW;
        } else {
            rd_log_warning!(
                "Warning: ambiguous stereochemistry - zero final chiral volume - at atom {atom} ignored"
            );
            return Ok(None);
        }
    }
    Ok(Some(res))
}

/// `MolOps::assignChiralTypesFromBondDirs`
pub fn assign_chiral_types_from_bond_dirs(
    mol: &mut ROMol,
    conf_id: i32,
    replace_existing_tags: bool,
) -> Result<()> {
    if mol.num_conformers() == 0 {
        return detect_atropisomer_chirality(mol, None, replace_existing_tags);
    }
    let cpos = conf_position(mol, conf_id)?;
    let conf = mol.conformers()[cpos].clone();
    let mut atoms_set = vec![false; mol.num_atoms() as usize];
    for bi in 0..mol.num_bonds() {
        let dir = mol.bond(bi).bond_dir();
        let atom = mol.bond(bi).begin_atom_idx();
        if dir == BondDir::Unknown {
            if atoms_set[atom as usize] || replace_existing_tags {
                mol.atom_mut(atom).set_chiral_tag(ChiralType::Unspecified);
                atoms_set[atom as usize] = true;
            }
        } else if matches!(dir, BondDir::BeginWedge | BondDir::BeginDash) {
            if atoms_set[atom as usize]
                || (!replace_existing_tags
                    && mol.atom(atom).chiral_tag() != ChiralType::Unspecified)
            {
                continue;
            }
            if mol.atom(atom).needs_update_property_cache() {
                mol.update_atom_property_cache(atom, false)?;
            }
            let code = atom_chiral_type_from_bond_dir_pseudo3d(mol, bi, &conf)?
                .unwrap_or(ChiralType::Unspecified);
            if code != ChiralType::Unspecified {
                atoms_set[atom as usize] = true;
            }
            mol.atom_mut(atom).set_chiral_tag(code);
            let a = mol.atom(atom);
            if mol.atom_degree(atom) == 3 && a.num_explicit_hs() == 0 && a.num_implicit_hs()? == 1 {
                mol.atom_mut(atom).set_num_explicit_hs(1);
                mol.update_atom_property_cache(atom, true)?;
            }
        }
    }
    detect_atropisomer_chirality(mol, Some(cpos), replace_existing_tags)
}

/// `isWigglyBond`
pub fn is_wiggly_bond(mol: &ROMol, bond: u32, atom: u32) -> bool {
    let b = mol.bond(bond);
    b.begin_atom_idx() == atom
        && b.bond_type() == BondType::Single
        && (b.bond_dir() == BondDir::Unknown
            || b.props
                .get(cp::_UnknownStereo)
                .and_then(RDValue::as_int)
                .is_some_and(|v| v != 0))
}

fn voltest(v: &[Point3D; 6], x: usize, y: usize, z: usize) -> bool {
    v[x].dot_product(&v[y].cross_product(&v[z])) >= 0.0
}

/// `OctahedralPermFrom3D`
fn octahedral_perm_from_3d(pair: &[u8; 6], v: &[Point3D; 6]) -> u32 {
    let t = |x, y, z, a, b| if voltest(v, x, y, z) { a } else { b };
    match pair[0] {
        2 => match pair[2] {
            4 => t(0, 3, 4, 28, 27),
            5 => t(0, 2, 3, 25, 30),
            _ => t(0, 2, 3, 26, 29),
        },
        3 => match pair[1] {
            4 => t(0, 3, 4, 22, 21),
            5 => t(0, 1, 3, 19, 24),
            _ => t(0, 1, 3, 20, 23),
        },
        4 => match pair[1] {
            3 => t(0, 2, 4, 13, 12),
            5 => t(0, 1, 2, 6, 18),
            _ => t(0, 1, 2, 7, 17),
        },
        5 => match pair[1] {
            3 => t(0, 2, 3, 11, 9),
            4 => t(0, 1, 2, 3, 16),
            _ => t(0, 1, 2, 5, 15),
        },
        _ => match pair[1] {
            3 => t(0, 2, 3, 10, 8),
            4 => t(0, 1, 2, 1, 2),
            _ => t(0, 1, 2, 4, 14),
        },
    }
}

/// `assignNontetrahedralChiralTypeFrom3D` (tolerance 0.1)
fn assign_nontetrahedral_chiral_type_from_3d(
    mol: &mut ROMol,
    conf: &Conformer,
    atom: u32,
) -> Result<bool> {
    const TOLERANCE: f64 = 0.1;
    if mol.atom(atom).atomic_num() < 15 {
        return Ok(false);
    }
    for b in mol.atom_bond_indices(atom) {
        if is_wiggly_bond(mol, b, atom) {
            return Ok(false);
        }
    }
    let cen = conf.atom_pos(atom)?;
    let mut v = [Point3D::default(); 6];
    let mut count = 0usize;
    for nbr in mol.atom_neighbors(atom) {
        if count == 6 {
            return Ok(false);
        }
        v[count] = cen.direction_vector(&conf.atom_pos(nbr)?)?;
        count += 1;
    }
    if count < 3 {
        return Ok(false);
    }
    let mut pair = [0u8; 6];
    let mut pairs = 0;
    for i in 0..count {
        for j in (i + 1)..count {
            if v[i].dot_product(&v[j]) < -(1.0 - TOLERANCE) {
                if pair[i] != 0 || pair[j] != 0 {
                    return Ok(false);
                }
                pair[i] = (j + 1) as u8;
                pair[j] = (i + 1) as u8;
                pairs += 1;
            }
        }
    }
    let angle_lt_100 = |a: usize, b: usize| v[a].angle_to(&v[b]) < 100.0 * PI / 180.0;
    let t = |x, y, z, a: u32, b: u32| if voltest(&v, x, y, z) { a } else { b };
    let result: Option<(ChiralType, u32)> = match (pairs, count) {
        (1, 3) => {
            let perm = match pair[0] {
                0 => 3,
                2 => 2,
                _ => 1,
            };
            Some((ChiralType::SquarePlanar, perm))
        }
        (1, 4) => {
            use ChiralType::{Octahedral as Oh, TrigonalBipyramidal as Tbp};
            Some(if pair[0] == 2 {
                if angle_lt_100(2, 3) {
                    (Oh, t(0, 2, 3, 25, 29))
                } else {
                    (Tbp, t(0, 2, 3, 7, 8))
                }
            } else if pair[0] == 3 {
                if angle_lt_100(1, 3) {
                    (Oh, t(0, 1, 3, 19, 23))
                } else {
                    (Tbp, t(0, 1, 3, 5, 6))
                }
            } else if pair[0] == 4 {
                if angle_lt_100(1, 2) {
                    (Oh, t(0, 1, 2, 6, 17))
                } else {
                    (Tbp, t(0, 1, 2, 3, 4))
                }
            } else if pair[1] == 3 {
                if angle_lt_100(0, 3) {
                    (Oh, t(0, 1, 3, 10, 8))
                } else {
                    (Tbp, t(1, 0, 3, 13, 14))
                }
            } else if pair[1] == 4 {
                if angle_lt_100(0, 2) {
                    (Oh, t(0, 1, 3, 1, 2))
                } else {
                    (Tbp, t(1, 0, 2, 10, 12))
                }
            } else if angle_lt_100(0, 1) {
                (Oh, t(0, 1, 3, 4, 14))
            } else {
                (Tbp, t(3, 0, 1, 16, 19))
            })
        }
        (1, 5) => {
            let perm = if pair[0] == 2 {
                t(0, 2, 3, 7, 8)
            } else if pair[0] == 3 {
                t(0, 1, 3, 5, 6)
            } else if pair[0] == 4 {
                t(0, 1, 2, 3, 4)
            } else if pair[0] == 5 {
                t(0, 1, 2, 1, 2)
            } else if pair[1] == 3 {
                t(1, 0, 3, 13, 14)
            } else if pair[1] == 4 {
                t(1, 0, 2, 10, 12)
            } else if pair[1] == 5 {
                t(1, 0, 2, 9, 11)
            } else if pair[2] == 4 {
                t(2, 0, 1, 16, 19)
            } else if pair[2] == 5 {
                t(2, 0, 1, 15, 20)
            } else {
                t(3, 0, 1, 17, 18)
            };
            Some((ChiralType::TrigonalBipyramidal, perm))
        }
        (2, 4) => {
            let perm = match pair[0] {
                2 => 2,
                3 => 1,
                _ => 3,
            };
            Some((ChiralType::SquarePlanar, perm))
        }
        (2, 5) | (3, 6) => Some((ChiralType::Octahedral, octahedral_perm_from_3d(&pair, &v))),
        _ => None,
    };
    let Some((tag, perm)) = result else {
        return Ok(false);
    };
    let a = mol.atom_mut(atom);
    a.set_chiral_tag(tag);
    a.props.set_val(cp::_chiralPermutation, perm);
    Ok(true)
}

/// `MolOps::assignChiralTypesFrom3D`
pub fn assign_chiral_types_from_3d(
    mol: &mut ROMol,
    conf_id: i32,
    replace_existing_tags: bool,
) -> Result<()> {
    const ZERO_VOLUME_TOL: f64 = 0.1;
    if mol.num_conformers() == 0 {
        return Ok(());
    }
    let cpos = conf_position(mol, conf_id)?;
    let conf = mol.conformers()[cpos].clone();
    if !conf.is_3d() {
        return Ok(());
    }
    mol.props.clear_val(cp::_StereochemDone);
    let allow_nontet = get_allow_nontetrahedral_chirality();

    let n = mol.num_atoms() as usize;
    let mut explicit_atoms = vec![false; n];
    for b in mol.bonds() {
        if matches!(b.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash) {
            explicit_atoms[b.begin_atom_idx() as usize] = true;
        }
    }
    for a in mol.atoms() {
        if a.chiral_tag() != ChiralType::Unspecified {
            explicit_atoms[a.idx() as usize] = true;
        }
    }

    for ai in 0..mol.num_atoms() {
        if !replace_existing_tags && mol.atom(ai).chiral_tag() != ChiralType::Unspecified {
            continue;
        }
        mol.atom_mut(ai).set_chiral_tag(ChiralType::Unspecified);
        let nz_degree = get_atom_nonzero_degree(mol, ai);
        let tnz_degree = nz_degree + mol.total_num_hs(ai, false)?;
        if nz_degree < 3 || tnz_degree > 6 {
            continue;
        }
        if allow_nontet && assign_nontetrahedral_chiral_type_from_3d(mol, &conf, ai)? {
            if !explicit_atoms[ai as usize] {
                mol.atom_mut(ai)
                    .props
                    .set_val(cp::_NonExplicit3DChirality, 1i32);
            }
            continue;
        }
        if tnz_degree > 4 {
            continue;
        }
        let anum = mol.atom(ai).atomic_num();
        if anum != 16 && anum != 34 && tnz_degree != 4 {
            continue;
        }
        let p0 = conf.atom_pos(ai)?;
        let mut nbrs: Vec<Point3D> = Vec::with_capacity(4);
        let mut wiggly = false;
        for b in mol.atom_bond_indices(ai) {
            if is_wiggly_bond(mol, b, ai) {
                wiggly = true;
                break;
            }
            if !bond_affects_atom_chirality(mol, b, ai) {
                continue;
            }
            nbrs.push(conf.atom_pos(mol.bond(b).other_atom_idx(ai)?)?);
        }
        if wiggly {
            continue;
        }
        let v1 = nbrs[0] - p0;
        let v2 = nbrs[1] - p0;
        let v3 = nbrs[2] - p0;
        let mut chiral_vol = v1.dot_product(&v2.cross_product(&v3));
        let mut set = false;
        let mut tag = ChiralType::Unspecified;
        if chiral_vol < -ZERO_VOLUME_TOL {
            tag = ChiralType::TetrahedralCW;
            set = true;
        } else if chiral_vol > ZERO_VOLUME_TOL {
            tag = ChiralType::TetrahedralCCW;
            set = true;
        } else if nbrs.len() == 4 {
            let v4 = nbrs[3] - p0;
            chiral_vol = -v1.dot_product(&v2.cross_product(&v4));
            if chiral_vol < -ZERO_VOLUME_TOL {
                tag = ChiralType::TetrahedralCW;
                set = true;
            } else if chiral_vol > ZERO_VOLUME_TOL {
                tag = ChiralType::TetrahedralCCW;
                set = true;
            }
        }
        mol.atom_mut(ai).set_chiral_tag(tag);
        if set && !explicit_atoms[ai as usize] {
            mol.atom_mut(ai)
                .props
                .set_val(cp::_NonExplicit3DChirality, 1i32);
        }
    }
    detect_atropisomer_chirality(mol, Some(cpos), replace_existing_tags)
}

/// `MolOps::detectBondStereochemistry(mol, confId)`
pub fn detect_bond_stereochemistry(mol: &mut ROMol, conf_id: i32) -> Result<()> {
    if mol.num_conformers() == 0 {
        return Ok(());
    }
    let cpos = conf_position(mol, conf_id)?;
    set_double_bond_neighbor_directions(mol, Some(cpos))
}

/// `MolOps::assignStereochemistryFrom3D`
pub fn assign_stereochemistry_from_3d(
    mol: &mut ROMol,
    conf_id: i32,
    replace_existing_tags: bool,
) -> Result<()> {
    if mol.num_conformers() == 0 || !mol.conformer(conf_id)?.is_3d() {
        return Ok(());
    }
    if mol.needs_update_property_cache() {
        mol.update_property_cache(false)?;
    }
    detect_bond_stereochemistry(mol, conf_id)?;
    assign_chiral_types_from_3d(mol, conf_id, replace_existing_tags)?;
    super::legacy::assign_stereochemistry(mol, replace_existing_tags, true, true)
}
