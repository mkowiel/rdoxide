//! Principal moments of inertia and the shape descriptors derived from them
//! (upstream `GraphMol/Descriptors/PMI.cpp`): NPR1/NPR2, PMI1-3, radius of
//! gyration, inertial shape factor, eccentricity, asphericity and spherocity.
//!
//! The moments of inertia come from the inertia tensor, the gyration based
//! quantities from the covariance matrix; both are computed with the atomic
//! masses as weights or unweighted, and cached on the molecule.
//!
//! Copyright (C) 2016 Greg Landrum.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::mol_transforms::{
    PrincipalAxes, compute_principal_axes_and_moments,
    compute_principal_axes_and_moments_from_gyration_matrix,
};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::fp::fmuladd;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

pub const NPR1_VERSION: &str = "1.0.0";
pub const NPR2_VERSION: &str = "1.0.0";
pub const PMI1_VERSION: &str = "1.0.0";
pub const PMI2_VERSION: &str = "1.0.0";
pub const PMI3_VERSION: &str = "1.0.0";
pub const RADIUS_OF_GYRATION_VERSION: &str = "1.0.0";
pub const INERTIAL_SHAPE_FACTOR_VERSION: &str = "1.0.0";
pub const ECCENTRICITY_VERSION: &str = "1.0.0";
pub const ASPHERICITY_VERSION: &str = "1.0.0";
pub const SPHEROCITY_INDEX_VERSION: &str = "1.0.0";

/// The three moments, smallest first.
type Moments = [f64; 3];

type Solver = fn(&mut ROMol, i32, bool, bool, Option<&[f64]>) -> Result<Option<PrincipalAxes>>;

fn require_conformer(mol: &ROMol) -> Result<()> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    Ok(())
}

/// The moments from `solver`, cached on the molecule under
/// `_PMI{1,2,3}{suffix}` (`suffix` tells the variants apart). `None` when the
/// eigenvalue calculation failed.
fn cached_moments(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
    solver: Solver,
    suffix: &str,
) -> Result<Option<Moments>> {
    require_conformer(mol)?;
    let names = ["_PMI1", "_PMI2", "_PMI3"].map(|n| {
        if use_atomic_masses {
            format!("{n}_mass{suffix}")
        } else {
            format!("{n}{suffix}")
        }
    });
    if !force {
        let cached = names.each_ref().map(|n| match mol.get_prop(n) {
            Some(RDValue::Double(d)) => Some(*d),
            _ => None,
        });
        if let [Some(a), Some(b), Some(c)] = cached {
            return Ok(Some([a, b, c]));
        }
    }
    let weights: Option<Vec<f64>> =
        use_atomic_masses.then(|| mol.atoms().iter().map(|a| a.mass()).collect());
    let ignore_hs = false;
    let Some(res) = solver(mol, conf_id, ignore_hs, force, weights.as_deref())? else {
        return Ok(None);
    };
    for (name, value) in names.iter().zip(res.moments) {
        mol.set_prop(name, value, true)?;
    }
    Ok(Some(res.moments))
}

fn inertia_moments(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
) -> Result<Option<Moments>> {
    cached_moments(
        mol,
        conf_id,
        use_atomic_masses,
        force,
        compute_principal_axes_and_moments,
        "",
    )
}

fn gyration_moments(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
) -> Result<Option<Moments>> {
    cached_moments(
        mol,
        conf_id,
        use_atomic_masses,
        force,
        compute_principal_axes_and_moments_from_gyration_matrix,
        "_cov",
    )
}

/// Applies `f` to the inertia moments; a failed eigenvalue calculation gives 0.
fn from_inertia(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
    f: impl FnOnce(Moments) -> f64,
) -> Result<f64> {
    require_conformer(mol)?;
    Ok(inertia_moments(mol, conf_id, use_atomic_masses, force)?.map_or(0.0, f))
}

/// Applies `f` to the gyration moments; a failed eigenvalue calculation gives 0.
fn from_gyration(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
    f: impl FnOnce(Moments) -> f64,
) -> Result<f64> {
    require_conformer(mol)?;
    Ok(gyration_moments(mol, conf_id, use_atomic_masses, force)?.map_or(0.0, f))
}

/// Normalized principal moments ratio 1 (I1/I3), from Sauer and Schwarz,
/// JCIM 43:987-1003 (2003).
pub fn npr1(mol: &mut ROMol, conf_id: i32, use_atomic_masses: bool, force: bool) -> Result<f64> {
    from_inertia(mol, conf_id, use_atomic_masses, force, |[pm1, _, pm3]| {
        if pm3 < 1e-8 { 0.0 } else { pm1 / pm3 }
    })
}

/// Normalized principal moments ratio 2 (I2/I3).
pub fn npr2(mol: &mut ROMol, conf_id: i32, use_atomic_masses: bool, force: bool) -> Result<f64> {
    from_inertia(mol, conf_id, use_atomic_masses, force, |[_, pm2, pm3]| {
        if pm3 < 1e-8 { 0.0 } else { pm2 / pm3 }
    })
}

/// First (smallest) principal moment of inertia.
pub fn pmi1(mol: &mut ROMol, conf_id: i32, use_atomic_masses: bool, force: bool) -> Result<f64> {
    from_inertia(mol, conf_id, use_atomic_masses, force, |m| m[0])
}

/// Second principal moment of inertia.
pub fn pmi2(mol: &mut ROMol, conf_id: i32, use_atomic_masses: bool, force: bool) -> Result<f64> {
    from_inertia(mol, conf_id, use_atomic_masses, force, |m| m[1])
}

/// Third (largest) principal moment of inertia.
pub fn pmi3(mol: &mut ROMol, conf_id: i32, use_atomic_masses: bool, force: bool) -> Result<f64> {
    from_inertia(mol, conf_id, use_atomic_masses, force, |m| m[2])
}

/// Radius of gyration, sqrt(t1 + t2 + t3) over the moments of the gyration
/// matrix (G. A. Arteca, Reviews in Computational Chemistry vol 9, eq. A4).
pub fn radius_of_gyration(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
) -> Result<f64> {
    from_gyration(mol, conf_id, use_atomic_masses, force, |[a, b, c]| {
        (a + b + c).sqrt()
    })
}

/// Inertial shape factor, I2 / (I1 * I3); 0 for planar or coordinate-less
/// molecules.
pub fn inertial_shape_factor(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
) -> Result<f64> {
    from_inertia(mol, conf_id, use_atomic_masses, force, |[pm1, pm2, pm3]| {
        if pm1 < 1e-4 || pm3 < 1e-4 {
            0.0
        } else {
            pm2 / (pm1 * pm3)
        }
    })
}

/// Molecular eccentricity, sqrt(I3^2 - I1^2) / I3; 0 without coordinates or
/// near degeneracy.
pub fn eccentricity(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
) -> Result<f64> {
    from_inertia(mol, conf_id, use_atomic_masses, force, |[pm1, _, pm3]| {
        // the difference of squares is fused like clang does for upstream
        let diff = fmuladd(pm3, pm3, -(pm1 * pm1));
        if pm3 < 1e-4 || diff < 1e-4 {
            0.0
        } else {
            diff.sqrt() / pm3
        }
    })
}

/// Asphericity from the gyration moments (Arteca, eq. A14).
pub fn asphericity(
    mol: &mut ROMol,
    conf_id: i32,
    use_atomic_masses: bool,
    force: bool,
) -> Result<f64> {
    from_gyration(mol, conf_id, use_atomic_masses, force, |[pm1, pm2, pm3]| {
        if pm3 < 1e-4 {
            0.0
        } else {
            let denom = pm1 + pm2 + pm3;
            0.5 * ((pm1 - pm2).powi(2) + (pm1 - pm3).powi(2) + (pm2 - pm3).powi(2))
                / (denom * denom)
        }
    })
}

/// Spherocity index, 3 * t1 / (t1 + t2 + t3) over the unweighted gyration
/// moments.
pub fn spherocity_index(mol: &mut ROMol, conf_id: i32, force: bool) -> Result<f64> {
    from_gyration(mol, conf_id, false, force, |[pm1, pm2, pm3]| {
        if pm3 < 1e-4 {
            0.0
        } else {
            3.0 * pm1 / (pm1 + pm2 + pm3)
        }
    })
}
