//! `rdkit.Geometry.rdGeometry` grids: `UniformGrid3D_` (and its
//! `UniformGrid3D` factory function), `UniformRealValueGrid3D` and the grid
//! utilities (upstream `Geometry/Wrap/UniformGrid3D.cpp`,
//! `UniformRealValueGrid3D.cpp`).

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyTuple};
use rdoxide::datastructs::value_vects::DiscreteValueVect;
use rdoxide::geometry::Point3D;
use rdoxide::geometry::grid::{
    UniformGrid3D, UniformRealValueGrid3D, compute_grid_centroid, find_grid_terminal_points,
    protrude_distance, tanimoto_distance, tversky_index,
};

use crate::datastructs_vv::{PyDiscreteValueVect, PyRealValueVect};
use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;

/// The bytes of a pickle given as `bytes` or `str`.
fn pkl_bytes(o: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(b) = o.cast::<PyBytes>() {
        return Ok(b.as_bytes().to_vec());
    }
    Ok(o.extract::<String>()?.into_bytes())
}

/// Class to represent a uniform three-dimensional
/// cubic grid. Each grid point can store a poisitive integer value. For the sake
/// of efficiency these value can either be binary, fit in 2, 4, 8 or 16 bits
#[pyclass(
    name = "UniformGrid3D_",
    module = "rdkit.Geometry.rdGeometry",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyUniformGrid3D {
    g: UniformGrid3D,
}

#[pymethods]
impl PyUniformGrid3D {
    /// pickle constructor
    #[new]
    #[pyo3(signature = (pkl))]
    fn new(pkl: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(PyUniformGrid3D {
            g: UniformGrid3D::from_binary(&pkl_bytes(pkl)?).py()?,
        })
    }

    /// Get the index to the grid point closest to the specified point
    #[pyo3(name = "GetGridPointIndex", signature = (point))]
    fn get_grid_point_index(&self, point: PyRef<'_, PyPoint3D>) -> i32 {
        self.g.grid_point_index(&point.p)
    }
    /// Get the index to the grid point with the three integer indices provided
    #[pyo3(name = "GetGridIndex", signature = (xi, yi, zi))]
    fn get_grid_index(&self, xi: u32, yi: u32, zi: u32) -> i32 {
        self.g.grid_index(xi, yi, zi)
    }
    /// Returns the integer indices of the grid index provided.
    #[pyo3(name = "GetGridIndices", signature = (idx))]
    fn get_grid_indices(&self, idx: u32) -> PyResult<(u32, u32, u32)> {
        self.g.grid_indices(idx).py()
    }
    /// Get the value at the closest grid point
    #[pyo3(name = "GetValPoint", signature = (pt))]
    fn get_val_point(&self, pt: PyRef<'_, PyPoint3D>) -> PyResult<i64> {
        self.g.val_at_point(&pt.p).py()
    }
    /// Get the value at the specified grid point
    #[pyo3(name = "GetVal", signature = (id))]
    fn get_val(&self, id: u32) -> PyResult<u32> {
        self.g.val(id).py()
    }
    /// Set the value at the specified grid point
    #[pyo3(name = "SetVal", signature = (id, val))]
    fn set_val(&mut self, id: u32, val: u32) -> PyResult<()> {
        self.g.set_val(id, val).py()
    }
    /// Set the value at grid point closest to the specified point
    #[pyo3(name = "SetValPoint", signature = (pt, val))]
    fn set_val_point(&mut self, pt: PyRef<'_, PyPoint3D>, val: u32) -> PyResult<()> {
        self.g.set_val_at_point(&pt.p, val).py()
    }
    /// Get the location of the specified grid point
    #[pyo3(name = "GetGridPointLoc", signature = (pointId))]
    #[allow(non_snake_case)]
    fn get_grid_point_loc(&self, pointId: u32) -> PyResult<PyPoint3D> {
        Ok(PyPoint3D {
            p: self.g.grid_point_loc(pointId).py()?,
        })
    }
    /// Get the size of the grid (number of grid points)
    #[pyo3(name = "GetSize")]
    fn get_size(&self) -> u32 {
        self.g.size()
    }
    /// Get the number of grid points along x-axis
    #[pyo3(name = "GetNumX")]
    fn get_num_x(&self) -> u32 {
        self.g.num_x()
    }
    /// Get the number of grid points along y-axis
    #[pyo3(name = "GetNumY")]
    fn get_num_y(&self) -> u32 {
        self.g.num_y()
    }
    /// Get the number of grid points along z-axis
    #[pyo3(name = "GetNumZ")]
    fn get_num_z(&self) -> u32 {
        self.g.num_z()
    }
    /// Get the location of the center of the grid
    #[pyo3(name = "GetOffset")]
    fn get_offset(&self) -> PyPoint3D {
        PyPoint3D { p: self.g.offset() }
    }
    /// Get the grid spacing
    #[pyo3(name = "GetSpacing")]
    fn get_spacing(&self) -> f64 {
        self.g.spacing()
    }
    /// Get the occupancy vector for the grid (a copy: changes to it do not reach the grid)
    #[pyo3(name = "GetOccupancyVect")]
    fn get_occupancy_vect(&self) -> PyDiscreteValueVect {
        PyDiscreteValueVect {
            v: self.g.occupancy_vect().clone(),
        }
    }
    /// Compare the parameters between two grid object
    #[pyo3(name = "CompareParams", signature = (other))]
    fn compare_params(&self, other: PyRef<'_, PyUniformGrid3D>) -> bool {
        self.g.compare_params(&other.g)
    }
    /// Set the occupancy on the grid for a sphere or specified radius
    ///  and multiple layers around this sphere, with decreasing values of
    /// occupancy
    #[pyo3(name = "SetSphereOccupancy", signature = (center, radius, stepSize, maxLayers = -1, ignoreOutOfBound = true))]
    #[allow(non_snake_case)]
    fn set_sphere_occupancy(
        &mut self,
        center: PyRef<'_, PyPoint3D>,
        radius: f64,
        stepSize: f64,
        maxLayers: i32,
        ignoreOutOfBound: bool,
    ) -> PyResult<()> {
        self.g
            .set_sphere_occupancy(&center.p, radius, stepSize, maxLayers, ignoreOutOfBound)
            .py()
    }

    fn __iadd__(&mut self, o: PyRef<'_, PyUniformGrid3D>) -> PyResult<()> {
        self.g.add_assign(&o.g).py()
    }
    fn __isub__(&mut self, o: PyRef<'_, PyUniformGrid3D>) -> PyResult<()> {
        self.g.sub_assign(&o.g).py()
    }
    fn __iand__(&mut self, o: PyRef<'_, PyUniformGrid3D>) -> PyResult<()> {
        self.g.intersect_assign(&o.g).py()
    }
    fn __ior__(&mut self, o: PyRef<'_, PyUniformGrid3D>) -> PyResult<()> {
        self.g.union_assign(&o.g).py()
    }

    fn __getinitargs__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, [PyBytes::new(py, &self.g.to_binary())])
    }
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let bin = PyBytes::new(py, &slf.borrow().g.to_binary());
        PyTuple::new(
            py,
            [
                slf.get_type().into_any(),
                PyTuple::new(py, [bin])?.into_any(),
            ],
        )
    }
}

/// Faking the constructor
#[pyfunction]
#[pyo3(name = "UniformGrid3D", signature = (dimX, dimY, dimZ, spacing = 0.5, valType = None, offSet = None))]
#[allow(non_snake_case)]
fn make_uniform_grid_3d(
    dimX: f64,
    dimY: f64,
    dimZ: f64,
    spacing: f64,
    valType: Option<&Bound<'_, PyAny>>,
    offSet: Option<PyRef<'_, PyPoint3D>>,
) -> PyResult<PyUniformGrid3D> {
    // DiscreteValueType.TWOBITVALUE by default
    let val_type: i64 = match valType {
        Some(v) => v.extract()?,
        None => 1,
    };
    if !(0..=4).contains(&val_type) {
        return Err(PyTypeError::new_err("invalid DiscreteValueType"));
    }
    let g = UniformGrid3D::new(
        dimX,
        dimY,
        dimZ,
        spacing,
        val_type as u32,
        offSet.map(|p| p.p),
    )
    .py()?;
    Ok(PyUniformGrid3D { g })
}

/// Write the grid to a grid file
#[pyfunction]
#[pyo3(name = "WriteGridToFile", signature = (grid, filename))]
fn write_grid_to_file(grid: PyRef<'_, PyUniformGrid3D>, filename: &str) -> PyResult<()> {
    grid.g.write_to_file(filename).py()
}

/// Compute the tversky index between two grid objects
#[pyfunction]
#[pyo3(name = "TverskyIndex", signature = (grid1, grid2, alpha, beta))]
fn py_tversky_index(
    grid1: PyRef<'_, PyUniformGrid3D>,
    grid2: PyRef<'_, PyUniformGrid3D>,
    alpha: f64,
    beta: f64,
) -> PyResult<f64> {
    tversky_index(&grid1.g, &grid2.g, alpha, beta).py()
}

/// Compute the tanimoto distance between two grid objects
#[pyfunction]
#[pyo3(name = "TanimotoDistance", signature = (grid1, grid2))]
fn py_tanimoto_distance(
    grid1: PyRef<'_, PyUniformGrid3D>,
    grid2: PyRef<'_, PyUniformGrid3D>,
) -> PyResult<f64> {
    tanimoto_distance(&grid1.g, &grid2.g).py()
}

/// Compute the protrude distance between two grid objects
#[pyfunction]
#[pyo3(name = "ProtrudeDistance", signature = (grid1, grid2))]
fn py_protrude_distance(
    grid1: PyRef<'_, PyUniformGrid3D>,
    grid2: PyRef<'_, PyUniformGrid3D>,
) -> PyResult<f64> {
    protrude_distance(&grid1.g, &grid2.g).py()
}

/// Compute the grid point at the center of sphere around a Point3D
#[pyfunction]
#[pyo3(name = "ComputeGridCentroid", signature = (grid, pt, windowRadius))]
#[allow(non_snake_case)]
fn py_compute_grid_centroid(
    grid: PyRef<'_, PyUniformGrid3D>,
    pt: PyRef<'_, PyPoint3D>,
    windowRadius: f64,
) -> PyResult<(f64, PyPoint3D)> {
    let (centroid, weight_sum) = compute_grid_centroid(&grid.g, &pt.p, windowRadius).py()?;
    Ok((weight_sum, PyPoint3D { p: centroid }))
}

/// Find a grid's terminal points (defined in the subshape algorithm).
#[pyfunction]
#[pyo3(name = "FindGridTerminalPoints", signature = (grid, windowRadius, inclusionFraction))]
#[allow(non_snake_case)]
fn py_find_grid_terminal_points<'py>(
    py: Python<'py>,
    grid: PyRef<'py, PyUniformGrid3D>,
    windowRadius: f64,
    inclusionFraction: f64,
) -> PyResult<Bound<'py, PyTuple>> {
    let pts: Vec<Point3D> =
        find_grid_terminal_points(&grid.g, windowRadius, inclusionFraction).py()?;
    PyTuple::new(py, pts.into_iter().map(|p| PyPoint3D { p }))
}

/// Class to represent a uniform three-dimensional
/// cubic grid. Each grid point can store a floating point value.
#[pyclass(
    name = "UniformRealValueGrid3D",
    module = "rdkit.Geometry.rdGeometry",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyUniformRealValueGrid3D {
    g: UniformRealValueGrid3D,
}

#[pymethods]
impl PyUniformRealValueGrid3D {
    /// `UniformRealValueGrid3D(pkl)`, `()`, `(grid)` (copy) or
    /// `(dimX, dimY, dimZ, spacing=0.5, offSet=None)`.
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    fn new(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let kw = |name: &str| -> PyResult<Option<Bound<'_, PyAny>>> {
            kwargs.map_or(Ok(None), |k| k.get_item(name))
        };
        let g = match args.len() {
            0 if kwargs.is_none_or(|k| k.is_empty()) => UniformRealValueGrid3D::unit().py()?,
            1 => {
                let a = args.get_item(0)?;
                if let Ok(other) = a.extract::<PyRef<'_, PyUniformRealValueGrid3D>>() {
                    other.g.clone()
                } else {
                    UniformRealValueGrid3D::from_binary(&pkl_bytes(&a)?).py()?
                }
            }
            _ => {
                let arg = |i: usize, name: &str| -> PyResult<Bound<'_, PyAny>> {
                    if i < args.len() {
                        args.get_item(i)
                    } else {
                        kw(name)?.ok_or_else(|| {
                            PyTypeError::new_err(format!("missing required argument '{name}'"))
                        })
                    }
                };
                let dim_x: f64 = arg(0, "dimX")?.extract()?;
                let dim_y: f64 = arg(1, "dimY")?.extract()?;
                let dim_z: f64 = arg(2, "dimZ")?.extract()?;
                let spacing: f64 = if args.len() > 3 {
                    args.get_item(3)?.extract()?
                } else {
                    kw("spacing")?.map_or(Ok(0.5), |v| v.extract())?
                };
                let offset = if args.len() > 4 {
                    Some(args.get_item(4)?)
                } else {
                    kw("offSet")?
                };
                let offset = match offset {
                    Some(o) if !o.is_none() => Some(o.extract::<PyRef<'_, PyPoint3D>>()?.p),
                    _ => None,
                };
                UniformRealValueGrid3D::new(dim_x, dim_y, dim_z, spacing, offset).py()?
            }
        };
        Ok(PyUniformRealValueGrid3D { g })
    }

    /// Get the index to the grid point closest to the specified point
    #[pyo3(name = "GetGridPointIndex", signature = (point))]
    fn get_grid_point_index(&self, point: PyRef<'_, PyPoint3D>) -> i32 {
        self.g.grid_point_index(&point.p)
    }
    /// Get the index to the grid point with the three integer indices provided
    #[pyo3(name = "GetGridIndex", signature = (xi, yi, zi))]
    fn get_grid_index(&self, xi: u32, yi: u32, zi: u32) -> i32 {
        self.g.grid_index(xi, yi, zi)
    }
    /// Returns the integer indices of the grid index provided.
    #[pyo3(name = "GetGridIndices", signature = (idx))]
    fn get_grid_indices(&self, idx: u32) -> PyResult<(u32, u32, u32)> {
        self.g.grid_indices(idx).py()
    }
    /// Get the value at the closest grid point
    #[pyo3(name = "GetValPoint", signature = (pt))]
    fn get_val_point(&self, pt: PyRef<'_, PyPoint3D>) -> PyResult<f64> {
        self.g.val_at_point(&pt.p).py()
    }
    /// Get the value at the specified grid point
    #[pyo3(name = "GetVal", signature = (id))]
    fn get_val(&self, id: u32) -> PyResult<f64> {
        self.g.val(id).py()
    }
    /// Set the value at the specified grid point
    #[pyo3(name = "SetVal", signature = (id, val))]
    fn set_val(&mut self, id: u32, val: f64) -> PyResult<()> {
        self.g.set_val(id, val).py()
    }
    /// Set the value at grid point closest to the specified point
    #[pyo3(name = "SetValPoint", signature = (pt, val))]
    fn set_val_point(&mut self, pt: PyRef<'_, PyPoint3D>, val: f64) -> PyResult<()> {
        self.g.set_val_at_point(&pt.p, val).py()
    }
    /// Get the location of the specified grid point
    #[pyo3(name = "GetGridPointLoc", signature = (pointId))]
    #[allow(non_snake_case)]
    fn get_grid_point_loc(&self, pointId: u32) -> PyResult<PyPoint3D> {
        Ok(PyPoint3D {
            p: self.g.grid_point_loc(pointId).py()?,
        })
    }
    /// Get the size of the grid (number of grid points)
    #[pyo3(name = "GetSize")]
    fn get_size(&self) -> u32 {
        self.g.size()
    }
    /// Get the number of grid points along x-axis
    #[pyo3(name = "GetNumX")]
    fn get_num_x(&self) -> u32 {
        self.g.num_x()
    }
    /// Get the number of grid points along y-axis
    #[pyo3(name = "GetNumY")]
    fn get_num_y(&self) -> u32 {
        self.g.num_y()
    }
    /// Get the number of grid points along z-axis
    #[pyo3(name = "GetNumZ")]
    fn get_num_z(&self) -> u32 {
        self.g.num_z()
    }
    /// Get the location of the center of the grid
    #[pyo3(name = "GetOffset")]
    fn get_offset(&self) -> PyPoint3D {
        PyPoint3D { p: self.g.offset() }
    }
    /// Get the grid spacing
    #[pyo3(name = "GetSpacing")]
    fn get_spacing(&self) -> f64 {
        self.g.spacing()
    }
    /// Get the occupancy vector for the grid (a copy: changes to it do not reach the grid)
    #[pyo3(name = "GetOccupancyVect")]
    fn get_occupancy_vect(&self) -> PyRealValueVect {
        PyRealValueVect {
            v: self.g.occupancy_vect().clone(),
        }
    }
    /// Compare the vector values between two grid objects.
    #[pyo3(name = "CompareVectors", signature = (other))]
    fn compare_vectors(&self, other: PyRef<'_, PyUniformRealValueGrid3D>) -> PyResult<bool> {
        self.g.compare_vectors(&other.g).py()
    }
    /// Compare the parameters between two grid object.
    #[pyo3(name = "CompareParams", signature = (other))]
    fn compare_params(&self, other: PyRef<'_, PyUniformRealValueGrid3D>) -> bool {
        self.g.compare_params(&other.g)
    }
    /// Compare the parameters and values between two grid objects.
    #[pyo3(name = "CompareGrids", signature = (other))]
    fn compare_grids(&self, other: PyRef<'_, PyUniformRealValueGrid3D>) -> PyResult<bool> {
        self.g.compare_grids(&other.g).py()
    }

    fn __iadd__(&mut self, o: PyRef<'_, PyUniformRealValueGrid3D>) -> PyResult<()> {
        self.g.add_assign(&o.g).py()
    }
    fn __isub__(&mut self, o: PyRef<'_, PyUniformRealValueGrid3D>) -> PyResult<()> {
        self.g.sub_assign(&o.g).py()
    }
    fn __iand__(&mut self, o: PyRef<'_, PyUniformRealValueGrid3D>) -> PyResult<()> {
        self.g.intersect_assign(&o.g).py()
    }
    fn __ior__(&mut self, o: PyRef<'_, PyUniformRealValueGrid3D>) -> PyResult<()> {
        self.g.union_assign(&o.g).py()
    }

    fn __getinitargs__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, [PyBytes::new(py, &self.g.to_binary())])
    }
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let bin = PyBytes::new(py, &slf.borrow().g.to_binary());
        PyTuple::new(
            py,
            [
                slf.get_type().into_any(),
                PyTuple::new(py, [bin])?.into_any(),
            ],
        )
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyUniformGrid3D>()?;
    m.add_class::<PyUniformRealValueGrid3D>()?;
    m.add_function(wrap_pyfunction!(make_uniform_grid_3d, m)?)?;
    m.add_function(wrap_pyfunction!(write_grid_to_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_tversky_index, m)?)?;
    m.add_function(wrap_pyfunction!(py_tanimoto_distance, m)?)?;
    m.add_function(wrap_pyfunction!(py_protrude_distance, m)?)?;
    m.add_function(wrap_pyfunction!(py_compute_grid_centroid, m)?)?;
    m.add_function(wrap_pyfunction!(py_find_grid_terminal_points, m)?)?;
    Ok(())
}

// the discrete grid shares its storage type with `DataStructs.DiscreteValueVect`
const _: fn(&DiscreteValueVect) = |_| {};
