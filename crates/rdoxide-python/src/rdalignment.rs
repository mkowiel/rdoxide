//! `rdoxide.Numerics.rdAlignment` (upstream `Numerics/Alignment/Wrap/rdAlignment.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::geometry::Point3D;
use rdoxide::numerics::alignment::align_points;

use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;
use crate::pyconv::{float64_array, float64_values};

fn value_error(msg: &str) -> PyErr {
    PyValueError::new_err(msg.to_string())
}

/// Points given as an `N x 3` numpy array, a sequence of 3-sequences or a
/// sequence of `Point3D`.
fn points(py: Python<'_>, obj: &Bound<'_, PyAny>) -> PyResult<Vec<Point3D>> {
    let np = py.import("numpy")?;
    if obj.is_instance(&np.getattr("ndarray")?)? {
        let (values, shape) = float64_values(py, obj)?;
        let &[_, 3] = shape.as_slice() else {
            return Err(value_error("Wrong dimension for the points array"));
        };
        return Ok(values
            .chunks_exact(3)
            .map(|c| Point3D::new(c[0], c[1], c[2]))
            .collect());
    }
    let Ok(n) = obj.len() else {
        return Err(value_error("non-sequence argument provided"));
    };
    if n == 0 {
        return Err(value_error("Empty sequence passed in"));
    }
    (0..n)
        .map(|i| {
            let item = obj.get_item(i)?;
            if let Ok(pt) = item.extract::<PyRef<'_, PyPoint3D>>() {
                return Ok(pt.p);
            }
            if obj.get_item(0)?.extract::<PyRef<'_, PyPoint3D>>().is_ok() {
                return Err(value_error("non-Point3D found in sequence of points"));
            }
            let row: Vec<f64> = item
                .try_iter()
                .and_then(|it| it.map(|v| v?.extract::<f64>()).collect())
                .map_err(|_| value_error("Wrong number of entries in the list of lists"))?;
            match row[..] {
                [x, y, z] => Ok(Point3D::new(x, y, z)),
                _ => Err(value_error("Wrong number of entries in the list of lists")),
            }
        })
        .collect()
}

/// Compute the optimal alignment (minimum RMSD) between two set of points using the quaternion algorithm
///
///  ARGUMENTS:
///
///     - refPoints : reference points specified as a N by 3 Numeric array or
///                   sequence of 3-sequences or sequence of Point3Ds
///     - probePoints : probe points to align to reference points - same format
///                   restrictions as reference points apply here
///     - weights : optional numeric vector or list of weights to associate to each pair of points
///     - reflect : reflect the probe points before attempting alignment
///     - maxIteration : maximum number of iterations for the eigen solver
///
///  RETURNS:
///
///     a 2-tuple:
///       - SSD value for the alignment
///       - the 4x4 transform matrix, as a Numeric array
#[pyfunction]
#[pyo3(name = "GetAlignmentTransform", signature = (refPoints, probePoints, weights = None, reflect = false, maxIterations = 50))]
#[allow(non_snake_case)]
fn get_alignment_transform<'py>(
    py: Python<'py>,
    refPoints: &Bound<'py, PyAny>,
    probePoints: &Bound<'py, PyAny>,
    weights: Option<&Bound<'py, PyAny>>,
    reflect: bool,
    maxIterations: u32,
) -> PyResult<(f64, Bound<'py, PyAny>)> {
    let (reference, probe) = (points(py, refPoints)?, points(py, probePoints)?);
    if reference.len() != probe.len() {
        return Err(value_error("Mis-match in number of points"));
    }
    let weights: Option<Vec<f64>> = match weights {
        None => None,
        Some(w) => {
            let np = py.import("numpy")?;
            let values: Vec<f64> = if w.is_instance(&np.getattr("ndarray")?)? {
                float64_values(py, w)?.0
            } else {
                w.try_iter()?
                    .map(|v| v?.extract::<f64>())
                    .collect::<PyResult<_>>()?
            };
            if values.is_empty() {
                None
            } else if values.len() != reference.len() {
                return Err(value_error(
                    "Number of weights supplied do not match the number of points",
                ));
            } else {
                Some(values)
            }
        }
    };
    let (ssd, trans) = align_points(
        &reference,
        &probe,
        weights.as_deref(),
        reflect,
        maxIterations,
    )
    .py()?;
    let flat: Vec<f64> = (0..4)
        .flat_map(|i| (0..4).map(move |j| (i, j)))
        .map(|(i, j)| trans.get_val(i, j))
        .collect();
    let matrix = float64_array(py, &flat)?.call_method1("reshape", ((4, 4),))?;
    Ok((ssd, matrix))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to align pairs of points in 3D",
    )?;
    m.add_function(wrap_pyfunction!(get_alignment_transform, m)?)?;
    Ok(())
}
