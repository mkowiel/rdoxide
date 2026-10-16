//! `rdoxide.Chem.rdShapeHelpers` (upstream `ShapeHelpers/Wrap/rdShapeHelpers.cpp`):
//! encoding and comparing the shapes of molecules.

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::geometry::Transform3D;
use rdoxide::graphmol::shape_helpers::{self as sh, ShapeGridOptions};

use crate::conformer::PyConformer;
use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;
use crate::geometry_grids::PyUniformGrid3D;
use crate::mol::PyMol;

/// The `trans` argument: a 4x4 double numpy array, anything else means no
/// transformation (upstream only checks `PyArray_Check`).
fn trans_arg(py: Python<'_>, trans: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Transform3D>> {
    let Some(t) = trans else { return Ok(None) };
    let np = py.import("numpy")?;
    if !t.is_instance(&np.getattr("ndarray")?)? {
        return Ok(None);
    }
    let shape: Vec<usize> = t.getattr("shape")?.extract()?;
    if shape.len() != 2 || shape[0] != 4 || shape[1] != 4 {
        return Err(PyValueError::new_err(
            "The transform has to be square matrix, of size 4x4",
        ));
    }
    if t.getattr("dtype")?.getattr("name")?.extract::<String>()? != "float64" {
        return Err(PyValueError::new_err(
            "Only double arrays allowed for transform object ",
        ));
    }
    let flat: Vec<f64> = np
        .call_method1("ascontiguousarray", (t, "float64"))?
        .call_method0("ravel")?
        .extract()?;
    let mut x = Transform3D::identity();
    x.data.copy_from_slice(&flat[..16]);
    Ok(Some(x))
}

fn bits(v: i64) -> PyResult<u32> {
    if !(0..=4).contains(&v) {
        return Err(PyTypeError::new_err("invalid DiscreteValueType"));
    }
    Ok(v as u32)
}

#[allow(clippy::too_many_arguments)]
fn options(
    grid_spacing: f64,
    bits_per_point: i64,
    vdw_scale: f64,
    step_size: f64,
    max_layers: i32,
    ignore_hs: bool,
) -> PyResult<ShapeGridOptions> {
    Ok(ShapeGridOptions {
        grid_spacing,
        bits_per_point: bits(bits_per_point)?,
        vdw_scale,
        step_size,
        max_layers,
        ignore_hs,
    })
}

/// Encode the shape of a molecule (one of its conformer) onto a grid
///
///  ARGUMENTS:
///
///    - mol : the molecule of interest
///    - grid : grid onto which the encoding is written
///    - confId : id of the conformation of interest on mol (defaults to the first one)
///    - trans : any transformation that needs to be used to encode onto the grid (note the molecule remains unchanged)
///    - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
///                 used in the encoding - grid points inside this sphere carry the maximum occupancy
///    - setpSize : thickness of the layers outside the base radius, the occupancy value is decreased
///                 from layer to layer from the maximum value
///    - maxLayers : the maximum number of layers - defaults to the number of bits
///                  used per grid point - e.g. two bits per grid point will allow 3 layers
///    - ignoreHs : when set, the contribution of Hs to the shape will be ignored
#[pyfunction]
#[pyo3(
    name = "EncodeShape",
    signature = (mol, grid, confId = -1, trans = None, vdwScale = 0.8, stepSize = 0.25, maxLayers = -1, ignoreHs = true)
)]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn encode_shape(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    grid: &Bound<'_, PyUniformGrid3D>,
    confId: i32,
    trans: Option<&Bound<'_, PyAny>>,
    vdwScale: f64,
    stepSize: f64,
    maxLayers: i32,
    ignoreHs: bool,
) -> PyResult<()> {
    let trans = trans_arg(py, trans)?;
    let conf = mol.mol.conformer(confId).py()?;
    let mut g = grid.borrow_mut();
    sh::encode_shape(
        &mol.mol,
        conf,
        &mut g.g,
        trans.as_ref(),
        vdwScale,
        stepSize,
        maxLayers,
        ignoreHs,
    )
    .py()
}

/// Compute the shape tversky index between two molecule based on a predefined alignment
///
///  ARGUMENTS:
///    - mol1 : The first molecule of interest
///    - mol2 : The second molecule of interest
///    - alpha : first parameter of the Tversky index
///    - beta : second parameter of the Tversky index
///    - confId1 : Conformer in the first molecule (defaults to first conformer)
///    - confId2 : Conformer in the second molecule (defaults to first conformer)
///    - gridSpacing : resolution of the grid used to encode the molecular shapes
///    - bitsPerPoint : number of bits used to encode the occupancy at each grid point
///                          defaults to two bits per grid point
///    - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
///                used in the encoding - grid points inside this sphere carry the maximum occupancy
///    - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
///                 from layer to layer from the maximum value
///    - maxLayers : the maximum number of layers - defaults to the number of bits
///                  used per grid point - e.g. two bits per grid point will allow 3 layers
///    - ignoreHs : when set, the contribution of Hs to the shape will be ignored
#[pyfunction]
#[pyo3(
    name = "ShapeTverskyIndex",
    signature = (mol1, mol2, alpha, beta, confId1 = -1, confId2 = -1, gridSpacing = 0.5, bitsPerPoint = 1, vdwScale = 0.8, stepSize = 0.25, maxLayers = -1, ignoreHs = true)
)]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn shape_tversky_index(
    py: Python<'_>,
    mol1: PyRef<'_, PyMol>,
    mol2: PyRef<'_, PyMol>,
    alpha: f64,
    beta: f64,
    confId1: i32,
    confId2: i32,
    gridSpacing: f64,
    bitsPerPoint: i64,
    vdwScale: f64,
    stepSize: f64,
    maxLayers: i32,
    ignoreHs: bool,
) -> PyResult<f64> {
    let o = options(
        gridSpacing,
        bitsPerPoint,
        vdwScale,
        stepSize,
        maxLayers,
        ignoreHs,
    )?;
    let (m1, m2) = (&mol1.mol, &mol2.mol);
    py.detach(|| sh::mol::tversky_index(m1, m2, alpha, beta, confId1, confId2, &o))
        .py()
}

/// Compute the shape tanimoto distance between two molecule based on a predefined alignment
///
///  ARGUMENTS:
///    - mol1 : The first molecule of interest
///    - mol2 : The second molecule of interest
///    - confId1 : Conformer in the first molecule (defaults to first conformer)
///    - confId2 : Conformer in the second molecule (defaults to first conformer)
///    - gridSpacing : resolution of the grid used to encode the molecular shapes
///    - bitsPerPoint : number of bits used to encode the occupancy at each grid point
///                          defaults to two bits per grid point
///    - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
///                used in the encoding - grid points inside this sphere carry the maximum occupancy
///    - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
///                 from layer to layer from the maximum value
///    - maxLayers : the maximum number of layers - defaults to the number of bits
///                  used per grid point - e.g. two bits per grid point will allow 3 layers
///    - ignoreHs : when set, the contribution of Hs to the shape will be ignored
#[pyfunction]
#[pyo3(
    name = "ShapeTanimotoDist",
    signature = (mol1, mol2, confId1 = -1, confId2 = -1, gridSpacing = 0.5, bitsPerPoint = 1, vdwScale = 0.8, stepSize = 0.25, maxLayers = -1, ignoreHs = true)
)]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn shape_tanimoto_dist(
    py: Python<'_>,
    mol1: PyRef<'_, PyMol>,
    mol2: PyRef<'_, PyMol>,
    confId1: i32,
    confId2: i32,
    gridSpacing: f64,
    bitsPerPoint: i64,
    vdwScale: f64,
    stepSize: f64,
    maxLayers: i32,
    ignoreHs: bool,
) -> PyResult<f64> {
    let o = options(
        gridSpacing,
        bitsPerPoint,
        vdwScale,
        stepSize,
        maxLayers,
        ignoreHs,
    )?;
    let (m1, m2) = (&mol1.mol, &mol2.mol);
    py.detach(|| sh::mol::tanimoto_distance(m1, m2, confId1, confId2, &o))
        .py()
}

/// Compute the shape protrude distance between two molecule based on a predefined alignment
///
///  ARGUMENTS:
///    - mol1 : The first molecule of interest
///    - mol2 : The second molecule of interest
///    - confId1 : Conformer in the first molecule (defaults to first conformer)
///    - confId2 : Conformer in the second molecule (defaults to first conformer)
///    - gridSpacing : resolution of the grid used to encode the molecular shapes
///    - bitsPerPoint : number of bit used to encode the occupancy at each grid point
///                          defaults to two bits per grid point
///    - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
///                used in the encoding - grid points inside this sphere carry the maximum occupancy
///    - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
///                 from layer to layer from the maximum value
///    - maxLayers : the maximum number of layers - defaults to the number of bits
///                  used per grid point - e.g. two bits per grid point will allow 3 layers
///    - ignoreHs : when set, the contribution of Hs to the shape will be ignored
///    - allowReordering : when set, the order will be automatically updated so that the value calculated
///                        is the protrusion of the smaller shape from the larger one.
#[pyfunction]
#[pyo3(
    name = "ShapeProtrudeDist",
    signature = (mol1, mol2, confId1 = -1, confId2 = -1, gridSpacing = 0.5, bitsPerPoint = 1, vdwScale = 0.8, stepSize = 0.25, maxLayers = -1, ignoreHs = true, allowReordering = true)
)]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn shape_protrude_dist(
    py: Python<'_>,
    mol1: PyRef<'_, PyMol>,
    mol2: PyRef<'_, PyMol>,
    confId1: i32,
    confId2: i32,
    gridSpacing: f64,
    bitsPerPoint: i64,
    vdwScale: f64,
    stepSize: f64,
    maxLayers: i32,
    ignoreHs: bool,
    allowReordering: bool,
) -> PyResult<f64> {
    let o = options(
        gridSpacing,
        bitsPerPoint,
        vdwScale,
        stepSize,
        maxLayers,
        ignoreHs,
    )?;
    let (m1, m2) = (&mol1.mol, &mol2.mol);
    py.detach(|| sh::mol::protrude_distance(m1, m2, confId1, confId2, &o, allowReordering))
        .py()
}

fn point_pair(
    a: rdoxide::geometry::Point3D,
    b: rdoxide::geometry::Point3D,
) -> (PyPoint3D, PyPoint3D) {
    (PyPoint3D { p: a }, PyPoint3D { p: b })
}

/// Compute the size of the box that can fit the conformations, and offset
///    of the box from the origin
#[pyfunction]
#[pyo3(name = "ComputeConfDimsAndOffset", signature = (conf, trans = None, padding = 2.0))]
fn compute_conf_dims_and_offset(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    trans: Option<&Bound<'_, PyAny>>,
    padding: f64,
) -> PyResult<(PyPoint3D, PyPoint3D)> {
    let trans = trans_arg(py, trans)?;
    let (dims, offset) = conf.borrow().with(py, |c| {
        sh::compute_conf_dims_and_offset(c, trans.as_ref(), padding).py()
    })??;
    Ok(point_pair(dims, offset))
}

/// Compute the lower and upper corners of a cuboid that will fit the conformer
#[pyfunction]
#[pyo3(name = "ComputeConfBox", signature = (conf, trans = None, padding = 2.0))]
fn compute_conf_box(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    trans: Option<&Bound<'_, PyAny>>,
    padding: f64,
) -> PyResult<(PyPoint3D, PyPoint3D)> {
    let trans = trans_arg(py, trans)?;
    let (lo, hi) = conf.borrow().with(py, |c| {
        sh::compute_conf_box(c, trans.as_ref(), padding).py()
    })??;
    Ok(point_pair(lo, hi))
}

/// Compute the union of two boxes, so that all the points in both boxes are
///     contained in the new box
#[pyfunction]
#[pyo3(name = "ComputeUnionBox", signature = (box1, box2))]
fn compute_union_box(
    box1: &Bound<'_, PyTuple>,
    box2: &Bound<'_, PyTuple>,
) -> PyResult<(PyPoint3D, PyPoint3D)> {
    if box1.len() != 2 || box2.len() != 2 {
        return Err(PyValueError::new_err(
            "In correct format for one of the box: expecting a tuple of two Point3D",
        ));
    }
    let p = |t: &Bound<'_, PyTuple>, i: usize| -> PyResult<rdoxide::geometry::Point3D> {
        Ok(t.get_item(i)?.extract::<PyRef<'_, PyPoint3D>>()?.p)
    };
    let (lo, hi) = sh::compute_union_box(&p(box1, 0)?, &p(box1, 1)?, &p(box2, 0)?, &p(box2, 1)?);
    Ok(point_pair(lo, hi))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to encode and compare the shapes of molecules",
    )?;
    m.add_function(wrap_pyfunction!(encode_shape, m)?)?;
    m.add_function(wrap_pyfunction!(shape_tversky_index, m)?)?;
    m.add_function(wrap_pyfunction!(shape_tanimoto_dist, m)?)?;
    m.add_function(wrap_pyfunction!(shape_protrude_dist, m)?)?;
    m.add_function(wrap_pyfunction!(compute_conf_dims_and_offset, m)?)?;
    m.add_function(wrap_pyfunction!(compute_conf_box, m)?)?;
    m.add_function(wrap_pyfunction!(compute_union_box, m)?)?;
    Ok(())
}
