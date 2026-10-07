//! `rdkit.Geometry.rdGeometry` — port of upstream `Code/Geometry/Wrap/Point.cpp`
//! (Point3D, Point2D so far).

use pyo3::exceptions::{PyIndexError, PyRuntimeError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::geometry::{Point2D, Point3D};

/// A class to represent a three-dimensional point
/// The x, y, and z coordinates can be read and written using either attributes
/// (i.e. pt.x = 4) or indexing (i.e. pt[0] = 4).
#[pyclass(
    name = "Point3D",
    module = "rdkit.Geometry.rdGeometry",
    subclass,
    skip_from_py_object
)]
#[derive(Debug, Clone, Copy)]
pub struct PyPoint3D {
    pub p: Point3D,
}

fn index_error(idx: i64) -> PyErr {
    PyIndexError::new_err(idx.to_string())
}

#[pymethods]
impl PyPoint3D {
    #[new]
    #[pyo3(signature = (xv = 0.0, yv = 0.0, zv = 0.0))]
    fn py_new(xv: f64, yv: f64, zv: f64) -> Self {
        PyPoint3D {
            p: Point3D::new(xv, yv, zv),
        }
    }
    #[getter]
    fn get_x(&self) -> f64 {
        self.p.x
    }
    #[setter]
    fn set_x(&mut self, v: f64) {
        self.p.x = v;
    }
    #[getter]
    fn get_y(&self) -> f64 {
        self.p.y
    }
    #[setter]
    fn set_y(&mut self, v: f64) {
        self.p.y = v;
    }
    #[getter]
    fn get_z(&self) -> f64 {
        self.p.z
    }
    #[setter]
    fn set_z(&mut self, v: f64) {
        self.p.z = v;
    }
    fn __getitem__(&self, idx: i64) -> PyResult<f64> {
        match idx {
            0 | -3 => Ok(self.p.x),
            1 | -2 => Ok(self.p.y),
            2 | -1 => Ok(self.p.z),
            _ => Err(index_error(idx)),
        }
    }
    fn __len__(&self) -> usize {
        3
    }
    fn __getinitargs__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, [self.p.x, self.p.y, self.p.z])
    }
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<(Bound<'py, PyAny>, (f64, f64, f64))> {
        let p = slf.borrow().p;
        Ok((slf.get_type().into_any(), (p.x, p.y, p.z)))
    }
    fn __add__(&self, other: PyRef<'_, PyPoint3D>) -> Self {
        PyPoint3D {
            p: self.p + other.p,
        }
    }
    fn __sub__(&self, other: PyRef<'_, PyPoint3D>) -> Self {
        PyPoint3D {
            p: self.p - other.p,
        }
    }
    fn __mul__(&self, scale: f64) -> Self {
        PyPoint3D { p: self.p * scale }
    }
    fn __truediv__(&self, scale: f64) -> Self {
        PyPoint3D { p: self.p / scale }
    }
    fn __iadd__(&mut self, other: PyRef<'_, PyPoint3D>) {
        self.p += other.p;
    }
    fn __isub__(&mut self, other: PyRef<'_, PyPoint3D>) {
        self.p -= other.p;
    }
    fn __imul__(&mut self, scale: f64) {
        self.p *= scale;
    }
    fn __itruediv__(&mut self, scale: f64) {
        self.p /= scale;
    }
    fn __neg__(&self) -> Self {
        PyPoint3D { p: -self.p }
    }
    /// Normalize the vector (using L2 norm)
    #[pyo3(name = "Normalize")]
    fn normalize(&mut self) -> PyResult<()> {
        self.p
            .normalize()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }
    /// Length of the vector
    #[pyo3(name = "Length")]
    fn length(&self) -> f64 {
        self.p.length()
    }
    /// Distance from this point to another point
    #[pyo3(name = "Distance")]
    fn distance(&self, pt2: PyRef<'_, PyPoint3D>) -> f64 {
        (self.p - pt2.p).length()
    }
    /// Square of the length
    #[pyo3(name = "LengthSq")]
    fn length_sq(&self) -> f64 {
        self.p.length_sq()
    }
    /// Dot product with another point
    #[pyo3(name = "DotProduct")]
    fn dot_product(&self, other: PyRef<'_, PyPoint3D>) -> f64 {
        self.p.dot_product(&other.p)
    }
    /// determines the angle between a vector to this point
    /// from the origin and a vector to the other point.
    ///
    /// The angle is unsigned: the results of this call will always
    /// be between 0 and M_PI
    #[pyo3(name = "AngleTo")]
    fn angle_to(&self, other: PyRef<'_, PyPoint3D>) -> f64 {
        self.p.angle_to(&other.p)
    }
    /// determines the signed angle between a vector to this point
    /// from the origin and a vector to the other point.
    ///
    /// The results of this call will be between 0 and M_2_PI
    #[pyo3(name = "SignedAngleTo")]
    fn signed_angle_to(&self, other: PyRef<'_, PyPoint3D>) -> f64 {
        self.p.signed_angle_to(&other.p)
    }
    /// return a normalized direction vector from this point to another
    #[pyo3(name = "DirectionVector")]
    fn direction_vector(&self, other: PyRef<'_, PyPoint3D>) -> PyResult<Self> {
        Ok(PyPoint3D {
            p: self
                .p
                .direction_vector(&other.p)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))?,
        })
    }
    /// Get the cross product between two points
    #[pyo3(name = "CrossProduct")]
    fn cross_product(&self, other: PyRef<'_, PyPoint3D>) -> Self {
        PyPoint3D {
            p: self.p.cross_product(&other.p),
        }
    }
}

/// A class to represent a two-dimensional point
#[pyclass(
    name = "Point2D",
    module = "rdkit.Geometry.rdGeometry",
    subclass,
    skip_from_py_object
)]
#[derive(Debug, Clone, Copy)]
pub struct PyPoint2D {
    pub p: Point2D,
}

#[derive(FromPyObject)]
enum Point2DInit<'py> {
    P3(PyRef<'py, PyPoint3D>),
    Xy(f64, f64),
}

#[pymethods]
impl PyPoint2D {
    #[new]
    #[pyo3(signature = (*args))]
    fn py_new(args: &Bound<'_, PyTuple>) -> PyResult<Self> {
        Ok(match args.len() {
            0 => PyPoint2D {
                p: Point2D::default(),
            },
            1 => match args.get_item(0)?.extract::<Point2DInit<'_>>()? {
                Point2DInit::P3(p) => PyPoint2D {
                    p: Point2D::new(p.p.x, p.p.y),
                },
                Point2DInit::Xy(x, y) => PyPoint2D {
                    p: Point2D::new(x, y),
                },
            },
            _ => {
                let (x, y): (f64, f64) = args.extract()?;
                PyPoint2D {
                    p: Point2D::new(x, y),
                }
            }
        })
    }
    #[getter]
    fn get_x(&self) -> f64 {
        self.p.x
    }
    #[setter]
    fn set_x(&mut self, v: f64) {
        self.p.x = v;
    }
    #[getter]
    fn get_y(&self) -> f64 {
        self.p.y
    }
    #[setter]
    fn set_y(&mut self, v: f64) {
        self.p.y = v;
    }
    fn __getitem__(&self, idx: i64) -> PyResult<f64> {
        match idx {
            0 | -2 => Ok(self.p.x),
            1 | -1 => Ok(self.p.y),
            _ => Err(index_error(idx)),
        }
    }
    fn __len__(&self) -> usize {
        2
    }
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<(Bound<'py, PyAny>, (f64, f64))> {
        let p = slf.borrow().p;
        Ok((slf.get_type().into_any(), (p.x, p.y)))
    }
    fn __add__(&self, o: PyRef<'_, PyPoint2D>) -> Self {
        PyPoint2D {
            p: Point2D::new(self.p.x + o.p.x, self.p.y + o.p.y),
        }
    }
    fn __sub__(&self, o: PyRef<'_, PyPoint2D>) -> Self {
        PyPoint2D {
            p: Point2D::new(self.p.x - o.p.x, self.p.y - o.p.y),
        }
    }
    fn __mul__(&self, s: f64) -> Self {
        PyPoint2D {
            p: Point2D::new(self.p.x * s, self.p.y * s),
        }
    }
    fn __truediv__(&self, s: f64) -> Self {
        PyPoint2D {
            p: Point2D::new(self.p.x / s, self.p.y / s),
        }
    }
    fn __neg__(&self) -> Self {
        PyPoint2D {
            p: Point2D::new(-self.p.x, -self.p.y),
        }
    }
    #[pyo3(name = "Length")]
    fn length(&self) -> f64 {
        self.p.length()
    }
    #[pyo3(name = "LengthSq")]
    fn length_sq(&self) -> f64 {
        self.p.length_sq()
    }
    #[pyo3(name = "DotProduct")]
    fn dot_product(&self, other: PyRef<'_, PyPoint2D>) -> f64 {
        self.p.dot_product(&other.p)
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing geometry objects like points, lines, etc\n",
    )?;
    m.add_class::<PyPoint3D>()?;
    m.add_class::<PyPoint2D>()?;
    Ok(())
}
