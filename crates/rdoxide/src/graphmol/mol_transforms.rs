//! Geometric operations on conformers (upstream
//! `GraphMol/MolTransforms/MolTransforms.cpp`): centroids, principal axes and
//! moments, canonical orientations and internal coordinates (bond lengths,
//! angles and dihedrals, which can also be set by moving the atoms behind the
//! bond).
//!
//! Upstream's `Conformer` knows its molecule; here the functions that need the
//! topology or cache results take the molecule and a conformer id instead.
//!
//! Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::romol::ROMol;
use crate::numerics::eigen_3x3::symmetric_eigen_3x3;
use crate::rdgeneral::fp::fmuladd;
use crate::rdgeneral::{RDProps, RDValue, cxx_sort::libcxx_sort_by};
use crate::{RDKitError, Result};

/// Squared lengths below this count as identical coordinates.
const MIN_SQ_LENGTH: f64 = 1.0e-16;

/// Principal axes (columns) and moments of a conformer.
#[derive(Debug, Clone, PartialEq)]
pub struct PrincipalAxes {
    /// `axes[row][col]`: column `col` is the axis of `moments[col]`.
    pub axes: [[f64; 3]; 3],
    pub moments: [f64; 3],
}

fn check_weights(weights: Option<&[f64]>, conf: &Conformer) -> Result<()> {
    match weights {
        Some(w) if w.len() < conf.num_atoms() as usize => {
            Err(RDKitError::precondition("bad weights vector"))
        }
        _ => Ok(()),
    }
}

fn is_hydrogen(mol: &ROMol, idx: u32) -> bool {
    mol.atom(idx).atomic_num() == 1
}

/// `transformConformer(conf, trans)`
pub fn transform_conformer(conf: &mut Conformer, trans: &Transform3D) {
    for p in conf.positions_mut() {
        *p = trans.transform_point(p);
    }
}

/// `computeCentroid(conf, ignoreHs, weights)`
pub fn compute_centroid(
    mol: &ROMol,
    conf: &Conformer,
    ignore_hs: bool,
    weights: Option<&[f64]>,
) -> Result<Point3D> {
    check_weights(weights, conf)?;
    let mut res = Point3D::default();
    let mut w_sum = 0.0;
    for i in 0..conf.num_atoms() {
        if ignore_hs && is_hydrogen(mol, i) {
            continue;
        }
        let w = weights.map_or(1.0, |w| w[i as usize]);
        w_sum += w;
        res += conf.atom_pos(i)? * w;
    }
    res /= w_sum;
    Ok(res)
}

/// The six distinct terms of a symmetric 3x3 matrix.
#[derive(Default)]
struct SymmTerms {
    xx: f64,
    xy: f64,
    xz: f64,
    yy: f64,
    yz: f64,
    zz: f64,
}

impl SymmTerms {
    fn to_matrix(&self) -> [[f64; 3]; 3] {
        [
            [self.xx, self.xy, self.xz],
            [self.xy, self.yy, self.yz],
            [self.xz, self.yz, self.zz],
        ]
    }
}

/// Atoms taking part, with their position relative to `center` and weight.
fn weighted_offsets<'a>(
    mol: &'a ROMol,
    conf: &'a Conformer,
    center: Point3D,
    ignore_hs: bool,
    weights: Option<&'a [f64]>,
) -> impl Iterator<Item = Result<(Point3D, f64)>> + 'a {
    (0..mol.num_atoms())
        .filter(move |&i| !(ignore_hs && is_hydrogen(mol, i)))
        .map(move |i| {
            let loc = conf.atom_pos(i)? - center;
            Ok((loc, weights.map_or(1.0, |w| w[i as usize])))
        })
}

/// `computeCovarianceTerms`
fn covariance_terms(
    mol: &ROMol,
    conf: &Conformer,
    center: Point3D,
    normalize: bool,
    ignore_hs: bool,
    weights: Option<&[f64]>,
) -> Result<SymmTerms> {
    check_weights(weights, conf)?;
    let mut t = SymmTerms::default();
    let mut w_sum = 0.0;
    for item in weighted_offsets(mol, conf, center, ignore_hs, weights) {
        let (loc, w) = item?;
        w_sum += w;
        t.xx = fmuladd(w * loc.x, loc.x, t.xx);
        t.xy = fmuladd(w * loc.x, loc.y, t.xy);
        t.xz = fmuladd(w * loc.x, loc.z, t.xz);
        t.yy = fmuladd(w * loc.y, loc.y, t.yy);
        t.yz = fmuladd(w * loc.y, loc.z, t.yz);
        t.zz = fmuladd(w * loc.z, loc.z, t.zz);
    }
    if normalize {
        t.xx /= w_sum;
        t.xy /= w_sum;
        t.xz /= w_sum;
        t.yy /= w_sum;
        t.yz /= w_sum;
        t.zz /= w_sum;
    }
    Ok(t)
}

/// `computeInertiaTerms`
fn inertia_terms(
    mol: &ROMol,
    conf: &Conformer,
    center: Point3D,
    ignore_hs: bool,
    weights: Option<&[f64]>,
) -> Result<SymmTerms> {
    check_weights(weights, conf)?;
    let mut t = SymmTerms::default();
    for item in weighted_offsets(mol, conf, center, ignore_hs, weights) {
        let (loc, w) = item?;
        t.xx = fmuladd(w, fmuladd(loc.y, loc.y, loc.z * loc.z), t.xx);
        t.yy = fmuladd(w, fmuladd(loc.x, loc.x, loc.z * loc.z), t.yy);
        t.zz = fmuladd(w, fmuladd(loc.y, loc.y, loc.x * loc.x), t.zz);
        t.xy = fmuladd(-(w * loc.x), loc.y, t.xy);
        t.xz = fmuladd(-(w * loc.x), loc.z, t.xz);
        t.yz = fmuladd(-(w * loc.z), loc.y, t.yz);
    }
    Ok(t)
}

/// `getEigenValEigenVectHelper`: `None` (after logging) if the solver fails.
fn eigen_decomposition(terms: &SymmTerms) -> Option<PrincipalAxes> {
    match symmetric_eigen_3x3(&terms.to_matrix()) {
        Ok(e) => Some(PrincipalAxes {
            axes: e.vectors,
            moments: e.values,
        }),
        Err(_) => {
            crate::rd_log_error!("eigenvalue calculation did not converge");
            None
        }
    }
}

fn flat_axes(res: &PrincipalAxes) -> Vec<f64> {
    res.axes.iter().flatten().copied().collect()
}

/// A `PrincipalAxes` cached on the molecule, if both properties are present.
fn cached_axes(mol: &ROMol, axes_name: &str, moments_name: &str) -> Option<PrincipalAxes> {
    let (Some(RDValue::VecDouble(a)), Some(RDValue::VecDouble(m))) =
        (mol.get_prop(axes_name), mol.get_prop(moments_name))
    else {
        return None;
    };
    if a.len() != 9 || m.len() != 3 {
        return None;
    }
    let mut res = PrincipalAxes {
        axes: [[0.0; 3]; 3],
        moments: [m[0], m[1], m[2]],
    };
    for (i, v) in a.iter().enumerate() {
        res.axes[i / 3][i % 3] = *v;
    }
    Some(res)
}

fn cache_axes(mol: &mut ROMol, names: (&str, &str), res: &PrincipalAxes) -> Result<()> {
    mol.set_prop(names.0, flat_axes(res), true)?;
    mol.set_prop(names.1, res.moments.to_vec(), true)
}

/// `computePrincipalAxesAndMoments(conf, axes, moments, ignoreHs, force, weights)`:
/// from the inertia tensor. `Ok(None)` is upstream's `false` (the eigenvalue
/// calculation failed). Unweighted results are cached on the molecule.
pub fn compute_principal_axes_and_moments(
    mol: &mut ROMol,
    conf_id: i32,
    ignore_hs: bool,
    force: bool,
    weights: Option<&[f64]>,
) -> Result<Option<PrincipalAxes>> {
    check_weights(weights, mol.conformer(conf_id)?)?;
    let names = if ignore_hs {
        ("_principalAxes_noH", "_principalMoments_noH")
    } else {
        ("_principalAxes", "_principalMoments")
    };
    if weights.is_none() && !force {
        if let Some(cached) = cached_axes(mol, names.0, names.1) {
            return Ok(Some(cached));
        }
    }
    let conf = mol.conformer(conf_id)?;
    let origin = compute_centroid(mol, conf, ignore_hs, weights)?;
    let terms = inertia_terms(mol, conf, origin, ignore_hs, weights)?;
    let Some(res) = eigen_decomposition(&terms) else {
        return Ok(None);
    };
    if weights.is_none() {
        cache_axes(mol, names, &res)?;
    }
    Ok(Some(res))
}

/// `computePrincipalAxesAndMomentsFromGyrationMatrix`: from the normalized
/// covariance matrix (the axes need not be right-handed).
pub fn compute_principal_axes_and_moments_from_gyration_matrix(
    mol: &mut ROMol,
    conf_id: i32,
    ignore_hs: bool,
    force: bool,
    weights: Option<&[f64]>,
) -> Result<Option<PrincipalAxes>> {
    check_weights(weights, mol.conformer(conf_id)?)?;
    let names = if ignore_hs {
        ("_principalAxes_noH_cov", "_principalMoments_noH_cov")
    } else {
        ("_principalAxes_cov", "_principalMoments_cov")
    };
    if weights.is_none() && !force {
        if let Some(cached) = cached_axes(mol, names.0, names.1) {
            return Ok(Some(cached));
        }
    }
    let conf = mol.conformer(conf_id)?;
    let origin = compute_centroid(mol, conf, ignore_hs, weights)?;
    let terms = covariance_terms(mol, conf, origin, true, ignore_hs, weights)?;
    let Some(res) = eigen_decomposition(&terms) else {
        return Ok(None);
    };
    if weights.is_none() {
        cache_axes(mol, names, &res)?;
    }
    Ok(Some(res))
}

/// `computeCanonicalTransform(conf, center, normalizeCovar, ignoreHs)`: the
/// transformation that aligns the principal axes with x, y and z; also the
/// eigenvalues (in the order of the eigenvectors of the covariance matrix).
pub fn compute_canonical_transform(
    mol: &ROMol,
    conf: &Conformer,
    center: Option<Point3D>,
    normalize_covar: bool,
    ignore_hs: bool,
) -> Result<(Transform3D, [f64; 3])> {
    let mut origin = match center {
        Some(c) => c,
        None => compute_centroid(mol, conf, ignore_hs, None)?,
    };
    let mut trans = Transform3D::identity();
    let mut eigen_values = [0.0; 3];
    // a single atom needs only the translation
    if conf.num_atoms() > 1 {
        let terms = covariance_terms(mol, conf, origin, normalize_covar, ignore_hs, None)?;
        if let Some(e) = eigen_decomposition(&terms) {
            let mut sorted: Vec<(usize, f64)> = e.moments.iter().copied().enumerate().collect();
            libcxx_sort_by(&mut sorted, |a, b| a.1 > b.1);
            for (col, &(col_sorted, value)) in sorted.iter().enumerate() {
                eigen_values[col_sorted] = value;
                for row in 0..3 {
                    trans.set_val(col, row, e.axes[row][col_sorted]);
                }
            }
        }
    }
    origin *= -1.0;

    // The eigenvector matrix may be left-handed (its determinant is -1); the
    // sign of the z axis is flipped then. This fix came from Andy Grant, RIP.
    let t = |i, j| trans.get_val(i, j);
    let test = t(0, 2) * (t(1, 0) * t(2, 1) - t(2, 0) * t(1, 1))
        + t(1, 2) * (t(2, 0) * t(0, 1) - t(2, 1) * t(0, 0))
        + t(2, 2) * (t(0, 0) * t(1, 1) - t(1, 0) * t(0, 1));
    if test < 0.0 {
        for i in 0..3 {
            trans.set_val(2, i, -trans.get_val(2, i));
        }
    }
    let moved = trans.transform_point(&origin);
    trans.set_translation(&moved);
    Ok((trans, eigen_values))
}

/// `canonicalizeConformer(conf, center, normalizeCovar, ignoreHs)`
pub fn canonicalize_conformer(
    mol: &mut ROMol,
    conf_id: i32,
    center: Option<Point3D>,
    normalize_covar: bool,
    ignore_hs: bool,
) -> Result<()> {
    let (trans, _) = compute_canonical_transform(
        mol,
        mol.conformer(conf_id)?,
        center,
        normalize_covar,
        ignore_hs,
    )?;
    transform_conformer(mol.conformer_mut(conf_id)?, &trans);
    Ok(())
}

/// `canonicalizeMol(mol, normalizeCovar, ignoreHs)`: every conformer.
pub fn canonicalize_mol(mol: &mut ROMol, normalize_covar: bool, ignore_hs: bool) -> Result<()> {
    let ids: Vec<u32> = mol.conformers().iter().map(Conformer::id).collect();
    for id in ids {
        canonicalize_conformer(mol, id as i32, None, normalize_covar, ignore_hs)?;
    }
    Ok(())
}

const SOURCE_FILE: &str = "Code/GraphMol/MolTransforms/MolTransforms.cpp";
const ATOM_ARG_NAMES: [&str; 4] = ["iAtomId", "jAtomId", "kAtomId", "lAtomId"];

/// `URANGE_CHECK(xAtomId, pos.size())` for each atom in turn; `first_line` is
/// the source line of the first check (the error text names the line).
fn check_range(positions: &[Point3D], atoms: &[u32], first_line: u32) -> Result<()> {
    for (n, &a) in atoms.iter().enumerate() {
        if a as usize >= positions.len() {
            return Err(RDKitError::range_at(
                ATOM_ARG_NAMES[n],
                a,
                positions.len(),
                SOURCE_FILE,
                first_line + n as u32,
            ));
        }
    }
    Ok(())
}

fn identical_coordinates(a: &str, b: &str) -> RDKitError {
    RDKitError::Value(format!("atoms {a} and {b} have identical 3D coordinates"))
}

/// `getBondLength`
pub fn get_bond_length(conf: &Conformer, i: u32, j: u32) -> Result<f64> {
    let pos = conf.positions();
    check_range(pos, &[i, j], 511)?;
    Ok((pos[i as usize] - pos[j as usize]).length())
}

/// `_toBeMovedIdxList`: the atoms reachable from `j` without passing through
/// `i`, including `j`, in increasing order.
fn atoms_to_move(mol: &ROMol, i: u32, j: u32) -> Vec<u32> {
    let mut visited = vec![false; mol.num_atoms() as usize];
    visited[i as usize] = true;
    visited[j as usize] = true;
    let mut stack = vec![j];
    while let Some(&top) = stack.last() {
        match mol.atom_neighbors(top).find(|&n| !visited[n as usize]) {
            Some(n) => {
                visited[n as usize] = true;
                stack.push(n);
            }
            None => {
                stack.pop();
            }
        }
    }
    (0..mol.num_atoms())
        .filter(|&a| visited[a as usize] && a != i)
        .collect()
}

/// `queryIsBondInRing`
fn bond_in_ring(mol: &ROMol, a: u32, b: u32) -> Result<bool> {
    let bond = mol
        .get_bond_between_atoms(a, b)?
        .ok_or_else(|| RDKitError::Value("atoms must be bonded".into()))?;
    Ok(mol.ring_info().num_bond_rings(bond.idx())? != 0)
}

/// `setBondLength`: moves `j` and everything attached to it.
pub fn set_bond_length(mol: &mut ROMol, conf_id: i32, i: u32, j: u32, value: f64) -> Result<()> {
    let pos = mol.conformer(conf_id)?.positions();
    check_range(pos, &[i, j], 520)?;
    if mol.get_bond_between_atoms(i, j)?.is_none() {
        return Err(RDKitError::Value("atoms i and j must be bonded".into()));
    }
    if bond_in_ring(mol, i, j)? {
        return Err(RDKitError::Value(
            "bond (i,j) must not belong to a ring".into(),
        ));
    }
    let mut v = pos[i as usize] - pos[j as usize];
    let orig_value = v.length();
    if orig_value <= 1.0e-8 {
        return Err(identical_coordinates("i", "j"));
    }
    let moved = atoms_to_move(mol, i, j);
    v *= value / orig_value - 1.0;
    let pos = mol.conformer_mut(conf_id)?.positions_mut();
    for a in moved {
        pos[a as usize] -= v;
    }
    Ok(())
}

/// `getAngleRad`
pub fn get_angle_rad(conf: &Conformer, i: u32, j: u32, k: u32) -> Result<f64> {
    let pos = conf.positions();
    check_range(pos, &[i, j, k], 548)?;
    let r_ji = pos[i as usize] - pos[j as usize];
    if r_ji.length_sq() <= MIN_SQ_LENGTH {
        return Err(identical_coordinates("i", "j"));
    }
    let r_jk = pos[k as usize] - pos[j as usize];
    if r_jk.length_sq() <= MIN_SQ_LENGTH {
        return Err(identical_coordinates("j", "k"));
    }
    Ok(r_ji.angle_to(&r_jk))
}

/// Rotates the atoms in `moved` by `angle` about the axis through `origin`.
fn rotate_atoms(pos: &mut [Point3D], moved: &[u32], origin: Point3D, axis: &Point3D, angle: f64) {
    let mut rot = Transform3D::identity();
    rot.set_rotation(angle, axis);
    for &a in moved {
        let shifted = pos[a as usize] - origin;
        pos[a as usize] = rot.transform_point(&shifted) + origin;
    }
}

/// `setAngleRad`: rotates `k` and everything attached to it about the normal
/// of the plane of the three atoms.
pub fn set_angle_rad(
    mol: &mut ROMol,
    conf_id: i32,
    i: u32,
    j: u32,
    k: u32,
    value: f64,
) -> Result<()> {
    let pos = mol.conformer(conf_id)?.positions();
    check_range(pos, &[i, j, k], 567)?;
    if mol.get_bond_between_atoms(j, i)?.is_none() {
        return Err(RDKitError::Value("atoms i and j must be bonded".into()));
    }
    if mol.get_bond_between_atoms(j, k)?.is_none() {
        return Err(RDKitError::Value("atoms j and k must be bonded".into()));
    }
    if bond_in_ring(mol, j, i)? && bond_in_ring(mol, j, k)? {
        return Err(RDKitError::Value(
            "bonds (i,j) and (j,k) must not both belong to a ring".into(),
        ));
    }
    let r_ji = pos[i as usize] - pos[j as usize];
    if r_ji.length_sq() <= MIN_SQ_LENGTH {
        return Err(identical_coordinates("i", "j"));
    }
    let r_jk = pos[k as usize] - pos[j as usize];
    if r_jk.length_sq() <= MIN_SQ_LENGTH {
        return Err(identical_coordinates("j", "k"));
    }
    // only the difference to the current angle has to be applied
    let delta = value - r_ji.angle_to(&r_jk);
    let origin = pos[j as usize];
    // upstream builds the axis end point and subtracts the origin again
    let mut axis = (r_ji.cross_product(&r_jk) + origin) - origin;
    axis.normalize()?;
    let moved = atoms_to_move(mol, j, k);
    rotate_atoms(
        mol.conformer_mut(conf_id)?.positions_mut(),
        &moved,
        origin,
        &axis,
        delta,
    );
    Ok(())
}

/// The vectors and normals shared by `getDihedralRad` and `setDihedralRad`.
struct DihedralFrame {
    r_jk: Point3D,
    n_ijk: Point3D,
    n_jkl: Point3D,
}

impl DihedralFrame {
    fn new(pos: &[Point3D], i: u32, j: u32, k: u32, l: u32) -> Result<Self> {
        let at = |a: u32| pos[a as usize];
        let r_ij = at(j) - at(i);
        if r_ij.length_sq() <= MIN_SQ_LENGTH {
            return Err(identical_coordinates("i", "j"));
        }
        let r_jk = at(k) - at(j);
        if r_jk.length_sq() <= MIN_SQ_LENGTH {
            return Err(identical_coordinates("j", "k"));
        }
        let r_kl = at(l) - at(k);
        if r_kl.length_sq() <= MIN_SQ_LENGTH {
            return Err(identical_coordinates("k", "l"));
        }
        Ok(DihedralFrame {
            r_jk,
            n_ijk: r_ij.cross_product(&r_jk),
            n_jkl: r_jk.cross_product(&r_kl),
        })
    }

    /// The signed dihedral angle (atan2 rather than acos).
    fn angle(&self) -> f64 {
        let m = self.n_ijk.cross_product(&self.r_jk);
        let n_ijk_sq = self.n_ijk.length_sq();
        let n_jkl_sq = self.n_jkl.length_sq();
        -(m.dot_product(&self.n_jkl) / (n_jkl_sq * m.length_sq()).sqrt())
            .atan2(self.n_ijk.dot_product(&self.n_jkl) / (n_ijk_sq * n_jkl_sq).sqrt())
    }
}

/// `getDihedralRad`
pub fn get_dihedral_rad(conf: &Conformer, i: u32, j: u32, k: u32, l: u32) -> Result<f64> {
    let pos = conf.positions();
    check_range(pos, &[i, j, k, l], 621)?;
    Ok(DihedralFrame::new(pos, i, j, k, l)?.angle())
}

/// `setDihedralRad`: rotates `k` and everything attached to it about the
/// (j,k) bond.
pub fn set_dihedral_rad(
    mol: &mut ROMol,
    conf_id: i32,
    i: u32,
    j: u32,
    k: u32,
    l: u32,
    value: f64,
) -> Result<()> {
    let pos = mol.conformer(conf_id)?.positions();
    check_range(pos, &[i, j, k, l], 654)?;
    if mol.get_bond_between_atoms(j, k)?.is_none() {
        return Err(RDKitError::Value("atoms j and k must be bonded".into()));
    }
    if bond_in_ring(mol, j, k)? {
        return Err(RDKitError::Value(
            "bond (j,k) must not belong to a ring".into(),
        ));
    }
    let frame = DihedralFrame::new(pos, i, j, k, l)?;
    let delta = value - frame.angle();
    let origin = pos[j as usize];
    let mut axis = pos[k as usize] - origin;
    axis.normalize()?;
    let moved = atoms_to_move(mol, j, k);
    rotate_atoms(
        mol.conformer_mut(conf_id)?.positions_mut(),
        &moved,
        origin,
        &axis,
        delta,
    );
    Ok(())
}

/// `getAngleDeg`
pub fn get_angle_deg(conf: &Conformer, i: u32, j: u32, k: u32) -> Result<f64> {
    Ok(180.0 / std::f64::consts::PI * get_angle_rad(conf, i, j, k)?)
}

/// `setAngleDeg`
pub fn set_angle_deg(
    mol: &mut ROMol,
    conf_id: i32,
    i: u32,
    j: u32,
    k: u32,
    value: f64,
) -> Result<()> {
    set_angle_rad(mol, conf_id, i, j, k, value / 180.0 * std::f64::consts::PI)
}

/// `getDihedralDeg`
pub fn get_dihedral_deg(conf: &Conformer, i: u32, j: u32, k: u32, l: u32) -> Result<f64> {
    Ok(180.0 / std::f64::consts::PI * get_dihedral_rad(conf, i, j, k, l)?)
}

/// `setDihedralDeg`
pub fn set_dihedral_deg(
    mol: &mut ROMol,
    conf_id: i32,
    i: u32,
    j: u32,
    k: u32,
    l: u32,
    value: f64,
) -> Result<()> {
    set_dihedral_rad(
        mol,
        conf_id,
        i,
        j,
        k,
        l,
        value / 180.0 * std::f64::consts::PI,
    )
}
