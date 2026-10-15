//! The ETKDG force fields (upstream `DistGeom/DistGeomUtils.cpp`): the 3D
//! force field that refines the embedded coordinates with torsion
//! preferences and knowledge terms, and the all-in-one field.

use std::collections::BTreeMap;

use super::{
    BoundsMatrix, ChiralSet, ChiralViolationContribs, DistViolationContribs, FourthDimContribs,
};
use crate::error::{RDKitError, Result};
use crate::forcefield::constraints::{AngleConstraintContrib, DistanceConstraintContrib};
use crate::forcefield::crystal_ff::{
    GaussianTorsionAngleContribs, PlanarityContribs, TorsionAngleContribs,
};
use crate::forcefield::mmff::terms::{DielModel, EleContrib};
use crate::forcefield::uff::InversionContrib;
use crate::forcefield::{ContribGroup, ForceField};
use crate::geometry::Point3D;
use crate::graphmol::forcefield_helpers::crystal_ff::{
    CrystalFFDetails, ImproperAtoms, TorsionParamKind, TorsionParams,
};

const KNOWN_DIST_TOL: f64 = 0.01;
const KNOWN_DIST_FORCE_CONSTANT: f64 = 100.0;

/// Partial charge products of atom pairs.
pub type ChargeProducts = BTreeMap<(u32, u32), f64>;

fn point(positions: &[f64], dim: usize, idx: usize) -> Point3D {
    Point3D::new(
        positions[dim * idx],
        positions[dim * idx + 1],
        positions[dim * idx + 2],
    )
}

fn distance(positions: &[f64], dim: usize, i: usize, j: usize) -> f64 {
    (point(positions, dim, i) - point(positions, dim, j)).length()
}

fn new_field(positions: &[f64], dim: usize) -> ForceField {
    let mut field = ForceField::new(dim);
    for p in positions.chunks_exact(dim) {
        field.push_point(p);
    }
    field
}

/// Basic-knowledge improper torsion (inversion) terms. `is_improper_constrained`
/// gets the central atoms.
fn add_improper_torsion_terms(
    ff: &mut ForceField,
    force_scaling_factor: f64,
    improper_atoms: &[ImproperAtoms],
    is_improper_constrained: &mut [bool],
) {
    let mut inversions = ContribGroup::new();
    for imp in improper_atoms {
        let a = imp.atoms;
        // the three permutations of the outer atoms
        for n in [[0, 1, 2, 3], [0, 1, 3, 2], [2, 1, 3, 0]] {
            inversions.push(InversionContrib::new(
                n.map(|i| a[i] as usize),
                imp.center_atomic_num,
                imp.bound_to_sp2_o,
                force_scaling_factor,
            ));
            is_improper_constrained[a[n[1]] as usize] = true;
        }
    }
    if !inversions.is_empty() {
        ff.push_contrib(Box::new(inversions));
    }
}

/// Experimental torsion terms; `atom_pairs` gets the end atoms of the
/// torsions if `exclude_torsions`.
fn add_experimental_torsion_terms(
    ff: &mut ForceField,
    details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    num_atoms: usize,
    exclude_torsions: bool,
) -> Result<()> {
    let mut done_bonds = vec![false; num_atoms * num_atoms];
    let mut gaussian = GaussianTorsionAngleContribs::new(3);
    let mut cosine = TorsionAngleContribs::new();
    let mut non_empty = false;
    for (t, atoms) in details.exp_torsion_atoms.iter().enumerate() {
        let [i, j, k, l] = atoms.map(|a| a as usize);
        let idx = if i < l {
            i * num_atoms + l
        } else {
            l * num_atoms + i
        };
        let bidx = if j < k {
            j * num_atoms + k
        } else {
            k * num_atoms + j
        };
        if done_bonds[bidx] {
            continue;
        }
        atom_pairs[idx] = exclude_torsions;
        match (&details.torsion_param_kind, &details.exp_torsion_angles[t]) {
            (TorsionParamKind::Gaussian, TorsionParams::Gaussian { scaling, .. }) => {
                let table = details.torsion_idx[t];
                gaussian.add_contrib(
                    [i, j, k, l],
                    details.phi_to_energy[table].clone(),
                    details.phi_to_grad[table].clone(),
                    *scaling,
                );
            }
            (TorsionParamKind::Cosine, TorsionParams::Cosine { signs, v }) => {
                let bad = || RDKitError::invariant("cosine torsion needs six terms");
                cosine.add_contrib(
                    [i, j, k, l],
                    v.as_slice().try_into().map_err(|_| bad())?,
                    signs.as_slice().try_into().map_err(|_| bad())?,
                );
            }
            _ => {
                return Err(RDKitError::invariant(
                    "torsion parameters do not match the parameter kind",
                ));
            }
        }
        non_empty = true;
        done_bonds[bidx] = true;
    }
    if non_empty {
        match details.torsion_param_kind {
            TorsionParamKind::Gaussian => ff.push_contrib(Box::new(gaussian)),
            TorsionParamKind::Cosine => ff.push_contrib(Box::new(cosine)),
        }
    }
    Ok(())
}

/// Bond-length constraints with padding at the current positions.
fn add12_terms(
    ff: &mut ForceField,
    details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    positions: &[f64],
    dim: usize,
    force_constant: f64,
    num_atoms: usize,
) -> Result<()> {
    let mut dists = ContribGroup::new();
    for &(a, b) in &details.bonds {
        let (i, j) = (a as usize, b as usize);
        atom_pairs[i.min(j) * num_atoms + i.max(j)] = true;
        let d = distance(positions, dim, i, j);
        dists.push(DistanceConstraintContrib::new(
            ff,
            i,
            j,
            false,
            d - KNOWN_DIST_TOL,
            d + KNOWN_DIST_TOL,
            force_constant,
        )?);
    }
    if !dists.is_empty() {
        ff.push_contrib(Box::new(dists));
    }
    Ok(())
}

/// 1-3 constraints with padding at the current positions, or linear angle
/// constraints for triple bonds; for angles centered on an improper torsion
/// atom the bounds of the matrix are used.
#[allow(clippy::too_many_arguments)]
fn add13_terms(
    ff: &mut ForceField,
    details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    positions: &[f64],
    dim: usize,
    force_constant: f64,
    is_improper_constrained: &[bool],
    use_basic_knowledge: bool,
    mmat: &BoundsMatrix,
    num_atoms: usize,
) -> Result<()> {
    let mut dists = ContribGroup::new();
    let mut angles = ContribGroup::new();
    for angle in &details.angles {
        let (i, j, k) = (angle[0] as usize, angle[1] as usize, angle[2] as usize);
        atom_pairs[i.min(k) * num_atoms + i.max(k)] = true;
        if use_basic_knowledge && angle[3] != 0 {
            // triple bonds
            angles.push(AngleConstraintContrib::new(
                ff,
                i,
                j,
                k,
                false,
                179.0,
                180.0,
                details.force_consts.k_term_angle,
            )?);
        } else if is_improper_constrained[j] {
            dists.push(DistanceConstraintContrib::new(
                ff,
                i,
                k,
                false,
                mmat.lower_bound(i, k),
                mmat.upper_bound(i, k),
                force_constant,
            )?);
        } else {
            let d = distance(positions, dim, i, k);
            dists.push(DistanceConstraintContrib::new(
                ff,
                i,
                k,
                false,
                d - KNOWN_DIST_TOL,
                d + KNOWN_DIST_TOL,
                force_constant,
            )?);
        }
    }
    if !angles.is_empty() {
        ff.push_contrib(Box::new(angles));
    }
    if !dists.is_empty() {
        ff.push_contrib(Box::new(dists));
    }
    Ok(())
}

/// Bounds-matrix constraints for the pairs that are not constrained
/// otherwise, tight ones between constrained atoms.
fn add_long_range_distance_constraints(
    ff: &mut ForceField,
    details: &CrystalFFDetails,
    atom_pairs: &[bool],
    positions: &[f64],
    dim: usize,
    known_distance_force_constant: f64,
    mmat: &BoundsMatrix,
    num_atoms: usize,
) -> Result<()> {
    let mut dists = ContribGroup::new();
    for i in 1..num_atoms {
        for j in 0..i {
            if atom_pairs[j * num_atoms + i] {
                continue;
            }
            let mut fdist = details.bounds_mat_force_scaling * 10.0;
            let mut l = mmat.lower_bound(i, j);
            let mut u = mmat.upper_bound(i, j);
            if !details.constrained_atoms.is_empty()
                && details.constrained_atoms[i]
                && details.constrained_atoms[j]
            {
                // constrained: use very tight bounds
                let d = distance(positions, dim, i, j);
                l = d - KNOWN_DIST_TOL;
                u = d + KNOWN_DIST_TOL;
                fdist = known_distance_force_constant;
            }
            dists.push(DistanceConstraintContrib::new(
                ff, i, j, false, l, u, fdist,
            )?);
        }
    }
    if !dists.is_empty() {
        ff.push_contrib(Box::new(dists));
    }
    Ok(())
}

fn check_inputs(
    mmat: &BoundsMatrix,
    positions: &[f64],
    dim: usize,
    details: &CrystalFFDetails,
) -> Result<usize> {
    let n = mmat.num_rows();
    if dim == 0 || positions.len() != n * dim {
        return Err(RDKitError::precondition("bad positions"));
    }
    if details.exp_torsion_atoms.len() != details.exp_torsion_angles.len() {
        return Err(RDKitError::precondition(
            "torsion atoms and angles differ in number",
        ));
    }
    Ok(n)
}

fn add_charges(ff: &mut ForceField, cpci: &ChargeProducts, divisor: f64) {
    let mut ele = EleContrib::new();
    for (&(i, j), &charge) in cpci {
        ele.add_term(
            i as usize,
            j as usize,
            charge / divisor,
            DielModel::Constant,
            false,
        );
    }
    ff.push_contrib(Box::new(ele));
}

/// The ETKDG force field on 3D `positions` (flat, 3 values per atom): torsion
/// preferences, improper torsions, bond and angle constraints at the current
/// geometry, and the bounds as long-range constraints. With `cpci`
/// Coulomb terms are added.
pub fn construct_3d_force_field(
    mmat: &BoundsMatrix,
    positions: &[f64],
    details: &CrystalFFDetails,
    cpci: Option<&ChargeProducts>,
) -> Result<ForceField> {
    let dim = 3;
    let n = check_inputs(mmat, positions, dim, details)?;
    let mut field = new_field(positions, dim);
    // the atom pairs that are 1-2, 1-3 or 1-4 restrained
    let mut atom_pairs = vec![false; n * n];
    // no 1-3 constraints for angles centered on an improper torsion atom
    let mut is_improper_constrained = vec![false; n];

    add_experimental_torsion_terms(&mut field, details, &mut atom_pairs, n, true)?;
    add_improper_torsion_terms(
        &mut field,
        details.force_consts.k_term_improper,
        &details.improper_atoms,
        &mut is_improper_constrained,
    );
    add12_terms(
        &mut field,
        details,
        &mut atom_pairs,
        positions,
        dim,
        KNOWN_DIST_FORCE_CONSTANT,
        n,
    )?;
    add13_terms(
        &mut field,
        details,
        &mut atom_pairs,
        positions,
        dim,
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        true,
        mmat,
        n,
    )?;
    add_long_range_distance_constraints(
        &mut field,
        details,
        &atom_pairs,
        positions,
        dim,
        KNOWN_DIST_FORCE_CONSTANT,
        mmat,
        n,
    )?;
    if let Some(cpci) = cpci {
        add_charges(&mut field, cpci, 1.0);
    }
    Ok(field)
}

/// Like [`construct_3d_force_field`] without the improper torsions and the
/// linear-angle constraints.
pub fn construct_plain_3d_force_field(
    mmat: &BoundsMatrix,
    positions: &[f64],
    details: &CrystalFFDetails,
) -> Result<ForceField> {
    let dim = 3;
    let n = check_inputs(mmat, positions, dim, details)?;
    let mut field = new_field(positions, dim);
    let mut atom_pairs = vec![false; n * n];
    let is_improper_constrained = vec![false; n];

    add_experimental_torsion_terms(&mut field, details, &mut atom_pairs, n, true)?;
    add12_terms(
        &mut field,
        details,
        &mut atom_pairs,
        positions,
        dim,
        KNOWN_DIST_FORCE_CONSTANT,
        n,
    )?;
    add13_terms(
        &mut field,
        details,
        &mut atom_pairs,
        positions,
        dim,
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        false,
        mmat,
        n,
    )?;
    add_long_range_distance_constraints(
        &mut field,
        details,
        &atom_pairs,
        positions,
        dim,
        KNOWN_DIST_FORCE_CONSTANT,
        mmat,
        n,
    )?;
    Ok(field)
}

/// A force field of only the improper torsions and linear-angle constraints.
pub fn construct_3d_improper_force_field(
    mmat: &BoundsMatrix,
    positions: &[f64],
    improper_atoms: &[ImproperAtoms],
    angles: &[[i32; 4]],
) -> Result<ForceField> {
    let dim = 3;
    let n = mmat.num_rows();
    if positions.len() != n * dim {
        return Err(RDKitError::precondition("bad positions"));
    }
    let mut field = new_field(positions, dim);
    let oob_force_scaling_factor = 10.0;
    let mut is_improper_constrained = vec![false; n];
    add_improper_torsion_terms(
        &mut field,
        oob_force_scaling_factor,
        improper_atoms,
        &mut is_improper_constrained,
    );
    // sp centers need an angle of 180 degrees
    let mut angle_contribs = ContribGroup::new();
    for angle in angles.iter().filter(|a| a[3] != 0) {
        angle_contribs.push(AngleConstraintContrib::new(
            &field,
            angle[0] as usize,
            angle[1] as usize,
            angle[2] as usize,
            false,
            179.0,
            180.0,
            oob_force_scaling_factor,
        )?);
    }
    if !angle_contribs.is_empty() {
        field.push_contrib(Box::new(angle_contribs));
    }
    Ok(field)
}

fn add_angle_terms(ff: &mut ForceField, angles: &[[i32; 4]], force_const: f64) -> Result<()> {
    let mut group = ContribGroup::new();
    for angle in angles.iter().filter(|a| a[3] != 0) {
        group.push(AngleConstraintContrib::new(
            ff,
            angle[0] as usize,
            angle[1] as usize,
            angle[2] as usize,
            false,
            179.0,
            180.0,
            force_const,
        )?);
    }
    if !group.is_empty() {
        ff.push_contrib(Box::new(group));
    }
    Ok(())
}

/// Distance-violation terms between all pairs (1-2 and 1-3 pairs get harmonic
/// constraints instead).
fn add_distance_terms(
    ff: &mut ForceField,
    force_const: f64,
    mmat: &BoundsMatrix,
    num_atoms: usize,
    bounds_mat_force_scaling: f64,
    extra_weights: Option<&BTreeMap<(usize, usize), f64>>,
    dist_mat: &[f64],
    fixed_pts: Option<&[bool]>,
) -> Result<()> {
    let mut dists = DistViolationContribs::new(ff.dimension());
    let mut harmonic = ContribGroup::new();
    for i in 1..num_atoms {
        for j in 0..i {
            if fixed_pts.is_some_and(|f| f[i] && f[j]) {
                continue;
            }
            let l = mmat.lower_bound(i, j);
            let u = mmat.upper_bound(i, j);
            let dist = dist_mat[i * num_atoms + j];
            let is_1213 = (dist - 1.0).abs() < 1e-4 || (dist - 2.0).abs() < 1e-4;
            let mut w = force_const;
            w *= bounds_mat_force_scaling;
            if let Some(&extra) = extra_weights.and_then(|e| e.get(&(i, j))) {
                w = extra;
            }
            if is_1213 {
                harmonic.push(DistanceConstraintContrib::new(ff, i, j, false, l, u, 10.0)?);
                continue;
            }
            dists.add_contrib(i, j, u, l, w);
        }
    }
    if !dists.is_empty() {
        ff.push_contrib(Box::new(dists));
    }
    if !harmonic.is_empty() {
        ff.push_contrib(Box::new(harmonic));
    }
    Ok(())
}

fn add_chirality_terms(
    ff: &mut ForceField,
    csets: &[ChiralSet],
    num_atoms: usize,
    weight_chiral: f64,
    weight_four_d: f64,
) {
    let mut chiral = ChiralViolationContribs::new(ff.dimension());
    let mut fourth = FourthDimContribs::new();
    for cset in csets {
        chiral.add_contrib(cset, weight_chiral);
    }
    for i in 0..num_atoms {
        fourth.add_contrib(i, weight_four_d);
    }
    if !chiral.is_empty() {
        ff.push_contrib(Box::new(chiral));
    }
    if !fourth.is_empty() {
        ff.push_contrib(Box::new(fourth));
    }
}

/// Planarity terms of the improper torsion centers.
fn add_planarity_terms(ff: &mut ForceField, force_const: f64, improper_atoms: &[ImproperAtoms]) {
    let mut planarity = PlanarityContribs::new();
    for imp in improper_atoms {
        let a = imp.atoms;
        for (i0, i2, i3) in [(0, 2, 3), (0, 3, 2), (2, 3, 0)] {
            planarity.add_contrib(
                [
                    a[i0] as usize,
                    a[1] as usize,
                    a[i2] as usize,
                    a[i3] as usize,
                ],
                force_const,
            );
        }
    }
    if !planarity.is_empty() {
        ff.push_contrib(Box::new(planarity));
    }
}

/// The all-in-one force field: distance violation terms between all pairs,
/// linear-angle constraints, chirality and fourth-dimension terms (in 4D) and
/// optionally Coulomb terms; torsions are added by [`add_torsion_terms`].
#[allow(clippy::too_many_arguments)]
pub fn construct_all_in_one_force_field(
    mmat: &BoundsMatrix,
    positions: &[f64],
    dim: usize,
    details: &CrystalFFDetails,
    csets: &[ChiralSet],
    cpci: Option<&ChargeProducts>,
    extra_weights: Option<&BTreeMap<(usize, usize), f64>>,
    fixed_pts: Option<&[bool]>,
) -> Result<ForceField> {
    let n = check_inputs(mmat, positions, dim, details)?;
    let mut field = new_field(positions, dim);
    add_angle_terms(
        &mut field,
        &details.angles,
        details.force_consts.k_term_angle,
    )?;
    if dim == 4 {
        add_chirality_terms(
            &mut field,
            csets,
            n,
            details.force_consts.chiral,
            details.force_consts.fourth_dim,
        );
    }
    add_distance_terms(
        &mut field,
        details.force_consts.distance,
        mmat,
        n,
        details.bounds_mat_force_scaling,
        extra_weights,
        &details.dist_mat,
        fixed_pts,
    )?;
    if let Some(cpci) = cpci {
        add_charges(&mut field, cpci, 500.0);
    }
    Ok(field)
}

/// Adds planarity terms (`do_k`) and experimental torsion terms (`do_k` or
/// `do_et`) to a force field.
pub fn add_torsion_terms(
    field: &mut ForceField,
    details: &CrystalFFDetails,
    do_k: bool,
    do_et: bool,
) -> Result<()> {
    if do_k {
        add_planarity_terms(
            field,
            details.force_consts.k_term_improper,
            &details.improper_atoms,
        );
    }
    if do_k || do_et {
        let n = field.num_points();
        let mut atom_pairs = vec![false; n * n];
        add_experimental_torsion_terms(field, details, &mut atom_pairs, n, false)?;
    }
    Ok(())
}
