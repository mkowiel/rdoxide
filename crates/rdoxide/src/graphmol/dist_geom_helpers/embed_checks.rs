//! Geometry checks that reject bad embeddings (upstream `Embedder.cpp`).
//!
//! Positions are flat arrays with `dim` values (3, or 4 while a fourth
//! dimension is used) per atom; only the first three are looked at.

use crate::distgeom::{BoundsMatrix, ChiralSet, ChiralSetStructureFlags, calc_chiral_volume};
use crate::geometry::Point3D;

pub(super) const ERROR_TOL: f64 = 0.00001;
/// these tolerances, all to detect and filter out bogus conformations, are a
/// delicate balance between sensitive enough to detect obviously bad
/// conformations but not so sensitive that a bunch of ok conformations get
/// filtered out, which slows down the whole conformation generation process
pub(super) const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;
pub(super) const MAX_MINIMIZED_E_PER_ATOM_AIO: f64 = 0.1;
const MIN_TETRAHEDRAL_CHIRAL_VOL: f64 = 0.50;
pub(super) const TETRAHEDRAL_CENTERINVOLUME_TOL: f64 = 0.30;
/// corresponds to the range [-0.1, 0.1]
pub(super) const AMPLITUDE_FOURTH_COORDINATE_IC: f64 = 0.2;

pub(super) fn point(positions: &[f64], dim: usize, idx: usize) -> Point3D {
    Point3D::new(
        positions[dim * idx],
        positions[dim * idx + 1],
        positions[dim * idx + 2],
    )
}

fn have_opposite_sign(a: f64, b: f64) -> bool {
    a.is_sign_negative() ^ b.is_sign_negative()
}

/// Whether the four neighbors span a volume around the center that is large
/// enough in every direction.
pub(super) fn volume_test(chiral_set: &ChiralSet, positions: &[f64], dim: usize) -> bool {
    let p = |i| point(positions, dim, i);
    let p0 = p(chiral_set.idx0);
    let unit = |idx| {
        let mut v = p0 - p(idx);
        // a zero vector stays zero, like upstream
        let _ = v.normalize();
        v
    };
    let v1 = unit(chiral_set.idx1);
    let v2 = unit(chiral_set.idx2);
    let v3 = unit(chiral_set.idx3);
    let v4 = unit(chiral_set.idx4);

    // be more tolerant of tetrahedral centers in multiple small rings
    let vol_scale =
        if chiral_set.structure_flags & (ChiralSetStructureFlags::InFusedSmallRings as u64) != 0 {
            0.25
        } else {
            1.0
        };
    let min_vol = vol_scale * MIN_TETRAHEDRAL_CHIRAL_VOL;
    v1.cross_product(&v2).dot_product(&v3).abs() >= min_vol
        && v1.cross_product(&v2).dot_product(&v4).abs() >= min_vol
        && v1.cross_product(&v3).dot_product(&v4).abs() >= min_vol
        && v2.cross_product(&v3).dot_product(&v4).abs() >= min_vol
}

fn same_side(
    v1: &Point3D,
    v2: &Point3D,
    v3: &Point3D,
    v4: &Point3D,
    p0: &Point3D,
    tol: f64,
) -> bool {
    let normal = (*v2 - *v1).cross_product(&(*v3 - *v1));
    let d1 = normal.dot_product(&(*v4 - *v1));
    let d2 = normal.dot_product(&(*p0 - *v1));
    if d1.abs() < tol || d2.abs() < tol {
        return false;
    }
    !((d1 < 0.0) ^ (d2 < 0.0))
}

/// Whether the center lies inside the tetrahedron of the other four points.
pub(super) fn center_in_volume(
    chiral_set: &ChiralSet,
    positions: &[f64],
    dim: usize,
    tol: f64,
) -> bool {
    // this happens for three-coordinate centers
    if chiral_set.idx0 == chiral_set.idx4 {
        return true;
    }
    let p = |i| point(positions, dim, i);
    let (p0, p1, p2, p3, p4) = (
        p(chiral_set.idx0),
        p(chiral_set.idx1),
        p(chiral_set.idx2),
        p(chiral_set.idx3),
        p(chiral_set.idx4),
    );
    same_side(&p1, &p2, &p3, &p4, &p0, tol)
        && same_side(&p2, &p3, &p4, &p1, &p0, tol)
        && same_side(&p3, &p4, &p1, &p2, &p0, tol)
        && same_side(&p4, &p1, &p2, &p3, &p0, tol)
}

/// Whether the distances between the atoms are (about) within the bounds.
pub(super) fn bounds_fulfilled(
    atoms: &[usize],
    mmat: &BoundsMatrix,
    positions: &[f64],
    dim: usize,
) -> bool {
    for i in 0..atoms.len().saturating_sub(1) {
        for j in i + 1..atoms.len() {
            let (a1, a2) = (atoms[i], atoms[j]);
            let d2 = (point(positions, dim, a1) - point(positions, dim, a2)).length();
            let lb = mmat.lower_bound(a1, a2);
            let ub = mmat.upper_bound(a1, a2);
            if (d2 < lb && (d2 - lb).abs() > 0.1 * ub) || (d2 > ub && (d2 - ub).abs() > 0.1 * ub) {
                return false;
            }
        }
    }
    true
}

/// Whether no atom pair is closer than its lower bound while under the
/// threshold (given as a squared distance).
pub(super) fn clash_check(
    positions: &[f64],
    dim: usize,
    mmat: &BoundsMatrix,
    threshold2: f64,
) -> bool {
    let n = positions.len() / dim;
    for i in 1..n {
        for j in 0..i {
            let dist2: f64 = (0..3)
                .map(|d| {
                    let diff = positions[dim * j + d] - positions[dim * i + d];
                    diff * diff
                })
                .sum();
            if dist2 >= threshold2 {
                continue;
            }
            let lb = mmat.lower_bound(i, j);
            if dist2 < lb * lb {
                return false;
            }
        }
    }
    true
}

/// The tetrahedral carbons keep a minimum volume around them and stay inside
/// it (part of github #971).
pub(super) fn check_tetrahedral_centers(
    tetrahedral_carbons: &[ChiralSet],
    positions: &[f64],
    dim: usize,
    tol: f64,
) -> bool {
    tetrahedral_carbons
        .iter()
        .all(|set| volume_test(set, positions, dim) && center_in_volume(set, positions, dim, tol))
}

/// The chiral volumes are not badly out of range.
pub(super) fn check_chiral_centers(
    chiral_centers: &[ChiralSet],
    positions: &[f64],
    dim: usize,
) -> bool {
    for set in chiral_centers {
        let vol = calc_chiral_volume(set.idx1, set.idx2, set.idx3, set.idx4, positions, dim);
        let lb = set.volume_lower_bound;
        let ub = set.volume_upper_bound;
        if (lb > 0.0 && vol < lb && (vol / lb < 0.8 || have_opposite_sign(vol, lb)))
            || (ub < 0.0 && vol > ub && (vol / ub < 0.8 || have_opposite_sign(vol, ub)))
        {
            return false;
        }
    }
    true
}

/// No double-bond end is part of a linear arrangement (and, optionally, no
/// sp2 center is flat in the wrong way).
pub(super) fn double_bond_geometry_checks(
    positions: &[f64],
    dim: usize,
    double_bond_ends: &[(u32, u32, u32)],
    improper_atoms: &[[u32; 4]],
    linear_tol: f64,
    do_sp2_centers: bool,
) -> bool {
    for &(i0, i1, i2) in double_bond_ends {
        let p0 = point(positions, dim, i0 as usize);
        let p1 = point(positions, dim, i1 as usize);
        let p2 = point(positions, dim, i2 as usize);
        // the arrangement a0 / a1 = a2 must not be a0 - a1 = a2
        let mut v1 = p1 - p0;
        let _ = v1.normalize();
        let mut v2 = p1 - p2;
        let _ = v2.normalize();
        if v1.dot_product(&v2) + 1.0 < linear_tol {
            return false;
        }
    }
    if do_sp2_centers {
        // the three neighbors of an sp2 center must not be arranged as
        // a0 - a1 = a2 with a3 above
        for atoms in improper_atoms {
            for (a, b, c) in [(0, 1, 2), (0, 1, 3), (2, 1, 3)] {
                let p1 = point(positions, dim, atoms[a] as usize);
                let p2 = point(positions, dim, atoms[b] as usize);
                let p3 = point(positions, dim, atoms[c] as usize);
                let mut p12 = p2 - p1;
                let mut p32 = p2 - p3;
                let _ = p12.normalize();
                let _ = p32.normalize();
                if p12.dot_product(&p32) + 1.0 < linear_tol {
                    return false;
                }
            }
        }
    }
    true
}
