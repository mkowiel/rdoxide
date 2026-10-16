//! `rdoxide.DistanceGeometry.DistGeom` (upstream `DistGeom/Wrap/DistGeom.cpp`):
//! basic distance-geometry operations on bounds matrices.

use std::collections::BTreeMap;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::distgeom::{
    BoundsMatrix, compute_initial_coords, construct_force_field, pick_random_dist_mat,
    triangle_smooth_bounds,
};
use rdoxide::graphmol::dist_geom_helpers::seed_rng;
use rdoxide::numerics::symm_matrix::SymmMatrix;

use crate::errors::IntoPyResult;
use crate::pyconv::{float64_array, float64_values};

const ERROR_TOL: f64 = 1e-5;

/// The bounds matrix held by a square, nonempty array of doubles.
fn bounds_matrix(py: Python<'_>, array: &Bound<'_, PyAny>) -> PyResult<BoundsMatrix> {
    let (values, shape) =
        float64_values(py, array).map_err(|_| PyValueError::new_err("Argument isn't an array"))?;
    if shape.len() != 2 || shape[0] != shape[1] {
        return Err(PyValueError::new_err("The array has to be square"));
    }
    let n = shape[0];
    if n == 0 {
        return Err(PyValueError::new_err(
            "The array has to have a nonzero size",
        ));
    }
    let mut mat = BoundsMatrix::new(n);
    for i in 0..n {
        for j in 0..n {
            mat.set_val(i, j, values[i * n + j]);
        }
    }
    Ok(mat)
}

/// Do triangle smoothing on a bounds matrix
///
///  ARGUMENTS:
///
///     - mat: a square Numeric array of doubles containing the bounds matrix, this matrix
///            *is* modified by the smoothing
///
///  RETURNS:
///
///     a boolean indicating whether or not the smoothing worked.
#[pyfunction]
#[pyo3(name = "DoTriangleSmoothing", signature = (boundsMatrix, tol = 0.0))]
#[allow(non_snake_case)]
fn do_triangle_smoothing(
    py: Python<'_>,
    boundsMatrix: &Bound<'_, PyAny>,
    tol: f64,
) -> PyResult<bool> {
    let mut mat = bounds_matrix(py, boundsMatrix)?;
    let ok = triangle_smooth_bounds(&mut mat, tol);
    let n = mat.num_rows();
    let values: Vec<f64> = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j)))
        .map(|(i, j)| mat.val(i, j))
        .collect();
    // the smoothed values are written back into the caller's array
    let smoothed = float64_array(py, &values)?.call_method1("reshape", (n, n))?;
    boundsMatrix.call_method1("__setitem__", (py.Ellipsis(), smoothed))?;
    Ok(ok)
}

/// Embed a bounds matrix and return the coordinates
///
///  ARGUMENTS:
///
///     - boundsMatrix: a square Numeric array of doubles containing the bounds matrix, this matrix
///            should already be smoothed
///     - maxIters: (optional) the maximum number of random distance matrices to try
///     - randomizeOnFailure: (optional) toggles using random coords if a matrix fails to embed
///     - numZeroFail: (optional) sets the number of zero eigenvalues to be considered a failure
///     - weights: (optional) a sequence of 3 sequences (i,j,weight) indicating elements of
///        the bounds matrix whose weights should be adjusted
///     - randomSeed: (optional) sets the random number seed used for embedding
///
///  RETURNS:
///
///     a Numeric array of doubles with the coordinates
#[pyfunction]
#[pyo3(
    name = "EmbedBoundsMatrix",
    signature = (boundsMatrix, maxIters = 10, randomizeOnFailure = false, numZeroFail = 2,
                 weights = None, randomSeed = -1)
)]
#[allow(non_snake_case)]
fn embed_bounds_matrix<'py>(
    py: Python<'py>,
    boundsMatrix: &Bound<'py, PyAny>,
    maxIters: i32,
    randomizeOnFailure: bool,
    numZeroFail: usize,
    weights: Option<Vec<(usize, usize, f64)>>,
    mut randomSeed: i32,
) -> PyResult<Bound<'py, PyAny>> {
    let bm = bounds_matrix(py, boundsMatrix)?;
    let n = bm.num_rows();
    let mut dist_mat = SymmMatrix::new(n);
    let mut positions = vec![0.0; n * 3];
    let mut got_coords = false;
    let mut iter = 0;
    while iter < maxIters && !got_coords {
        // every step restarts the generator from the seed, like upstream
        let mut rng = seed_rng(randomSeed).py()?;
        pick_random_dist_mat(&bm, &mut dist_mat, &mut rng).py()?;
        let mut rng = seed_rng(randomSeed).py()?;
        got_coords = compute_initial_coords(
            &dist_mat,
            &mut positions,
            3,
            &mut rng,
            randomizeOnFailure,
            numZeroFail,
        )
        .py()?;
        if randomSeed >= 0 {
            randomSeed += iter * 999;
        }
        iter += 1;
    }
    if !got_coords {
        return Err(PyValueError::new_err("could not embed matrix"));
    }
    let weight_map: BTreeMap<(usize, usize), f64> = weights
        .unwrap_or_default()
        .into_iter()
        .map(|(i, j, w)| ((i, j), w))
        .collect();
    let mut field = construct_force_field(
        &bm,
        &positions,
        3,
        &[],
        1.0,
        0.0,
        0.0,
        Some(&weight_map),
        5.0,
        None,
    )
    .py()?;
    if field.calc_energy() > ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            need_more = field.minimize(200, 1e-4, 1e-6).py()?;
        }
    }
    float64_array(py, field.positions())?.call_method1("reshape", (n, 3))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions for basic distance geometry operations",
    )?;
    m.add_function(wrap_pyfunction!(do_triangle_smoothing, m)?)?;
    m.add_function(wrap_pyfunction!(embed_bounds_matrix, m)?)?;
    Ok(())
}
