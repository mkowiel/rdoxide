//! Plane of best fit (upstream `GraphMol/Descriptors/PBF.cpp`), a measure of
//! the three-dimensionality of a molecule: the average distance of the atoms
//! from their least-squares plane. From Plane of Best Fit: A Novel Method to
//! Characterize the Three-Dimensionality of Molecules, Nicholas C. Firth,
//! Nathan Brown, and Julian Blagg, J. Chem. Inf. Model. 2012, 52 (10),
//! 2516-2525.
//!
//! Created by Nicholas Firth, November 2011; modified by Greg Landrum for
//! inclusion in the RDKit.
//!
//! Copyright (c) 2012, Institue of Cancer Research. All rights reserved.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-PBF-InstituteOfCancerResearch-BSD-3-Clause.txt).

use crate::Result;
use crate::geometry::Point3D;
use crate::graphmol::mol_transforms::compute_principal_axes_and_moments_from_gyration_matrix;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::fp::fmuladd;
use crate::{RDKitError, rd_log_error};

pub const PBF_VERSION: &str = "1.0.0";

/// `PBF(mol, confId)`: 0 for molecules with fewer than four atoms or without
/// 3D coordinates.
pub fn pbf(mol: &mut ROMol, conf_id: i32) -> Result<f64> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("molecule has no conformers"));
    }
    let num_atoms = mol.num_atoms();
    if num_atoms < 4 {
        return Ok(0.0);
    }
    let conf = mol.conformer(conf_id)?;
    if !conf.is_3d() {
        return Ok(0.0);
    }
    let points: Vec<Point3D> = conf.positions().to_vec();

    // the plane is spanned by the two largest principal axes of the covariance
    // matrix; its normal is the axis of the smallest moment
    let Some(axes) =
        compute_principal_axes_and_moments_from_gyration_matrix(mol, conf_id, false, false, None)?
    else {
        rd_log_error!("eigenvalue calculation did not converge");
        return Ok(0.0);
    };
    let normal = Point3D::new(axes.axes[0][0], axes.axes[1][0], axes.axes[2][0]);
    let mut origin = Point3D::default();
    let mut w_sum = 0.0;
    for p in &points {
        w_sum += 1.0;
        origin += *p;
    }
    origin /= w_sum;
    let offset = -normal.dot_product(&origin);

    // products are fused like clang does for upstream's expressions
    let denom = [normal.x, normal.y, normal.z]
        .iter()
        .fold(0.0, |acc, &c| fmuladd(c, c, acc))
        .sqrt();
    let distance = |p: &Point3D| {
        let numer = fmuladd(p.z, normal.z, fmuladd(p.x, normal.x, p.y * normal.y)) + offset;
        numer.abs() / denom
    };
    let mut total = 0.0;
    for p in &points {
        total += distance(p);
    }
    Ok(total / f64::from(num_atoms))
}
