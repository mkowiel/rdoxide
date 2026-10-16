//! `DiscreteValueVect` / `RealValueVect` (upstream
//! `DataStructs/Wrap/{DiscreteValueVect,RealValueVect}.cpp`).

use pyo3::IntoPyObjectExt;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple};
use rdoxide::datastructs::value_vects::{
    DiscreteValueVect, RealValueVect, compute_l1_norm_discrete, compute_l1_norm_real, cxx_max,
    cxx_min,
};

use crate::enums::{add_enum_exported, to_py};
use crate::errors::IntoPyResult;

const DVV_DOC: &str = "A container class for storing unsigned integer\nvalues within a particular range.\n\nThe length of the vector and type of its elements (determines the maximum value\nthat can be stored) are both set at construction time.\n\nAs you would expect, _DiscreteValueVects_ support a set of binary operations\nso you can do things like:\n  dvv3 = dvv1 & dvv2  the result contains the smallest value in each entry\n  dvv3 = dvv1 | dvv2  the result contains the largest value in each entry\n  dvv1 += dvv2     values are truncated when necessary\n  dvv3 = dvv1 + dvv2    values are truncated when necessary\n  dvv1 -= dvv3    would-be negative values are set to zero\n  dvv3 = dvv1 - dvv2    would-be negative values are set to zero\n\nElements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])\n\n";
const RVV_DOC: &str = "A container class for storing real\nvalues.\n\nThe length of the vector is set at construction time.\n\nAs you would expect, _RealValueVects_ support a set of binary operations\nso you can do things like:\n  rvv3 = rvv1 & rvv2  the result contains the smallest value in each entry\n  rvv3 = rvv1 | rvv2  the result contains the largest value in each entry\n  rvv1 += rvv2     \n  rvv3 = rvv1 + rvv2    \n  rvv1 -= rvv3    \n  rvv3 = rvv1 - rvv2    \n\nElements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])\n\n";

fn pkl_bytes(o: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(b) = o.cast::<PyBytes>() {
        return Ok(b.as_bytes().to_vec());
    }
    Ok(o.extract::<String>()?.into_bytes())
}

#[pyclass(
    name = "DiscreteValueVect",
    module = "rdoxide.DataStructs.cDataStructs",
    dict,
    sequence,
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PyDiscreteValueVect {
    pub v: DiscreteValueVect,
}

#[pymethods]
impl PyDiscreteValueVect {
    /// `DiscreteValueVect(valType, length)` or `DiscreteValueVect(pkl)`
    #[new]
    #[pyo3(signature = (valType, length = None))]
    #[allow(non_snake_case)]
    fn new(valType: &Bound<'_, PyAny>, length: Option<u32>) -> PyResult<Self> {
        match length {
            Some(n) => {
                let t: i64 = valType.extract()?;
                if !(0..=4).contains(&t) {
                    return Err(PyTypeError::new_err("invalid DiscreteValueType"));
                }
                Ok(PyDiscreteValueVect {
                    v: DiscreteValueVect::new(t as u32, n),
                })
            }
            None => Ok(PyDiscreteValueVect {
                v: DiscreteValueVect::from_binary(&pkl_bytes(valType)?).py()?,
            }),
        }
    }
    #[classattr]
    fn __doc__() -> &'static str {
        DVV_DOC
    }
    /// Get the number of entries in the vector
    fn __len__(&self) -> usize {
        self.v.length() as usize
    }
    /// Set the value at a specified location
    fn __setitem__(&mut self, i: u32, val: u32) -> PyResult<()> {
        self.v.set_val(i, val).py()
    }
    /// Get the value at a specified location
    fn __getitem__(&self, i: u32) -> PyResult<u32> {
        self.v.get_val(i).py()
    }
    fn __and__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        Ok(PyDiscreteValueVect {
            v: self.v.and(&o.v).py()?,
        })
    }
    fn __or__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        Ok(PyDiscreteValueVect {
            v: self.v.or(&o.v).py()?,
        })
    }
    fn __add__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        let mut v = self.v.clone();
        v.add_assign(&o.v).py()?;
        Ok(PyDiscreteValueVect { v })
    }
    fn __sub__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        let mut v = self.v.clone();
        v.sub_assign(&o.v).py()?;
        Ok(PyDiscreteValueVect { v })
    }
    fn __iadd__(&mut self, o: PyRef<'_, Self>) -> PyResult<()> {
        self.v.add_assign(&o.v).py()
    }
    fn __isub__(&mut self, o: PyRef<'_, Self>) -> PyResult<()> {
        self.v.sub_assign(&o.v).py()
    }
    /// Get the type of value stored in the vector
    #[pyo3(name = "GetValueType")]
    fn get_value_type(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "DiscreteValueType", self.v.value_type() as i64)
    }
    /// Get the sum of the values in the vector, basically L1 norm
    #[pyo3(name = "GetTotalVal")]
    fn get_total_val(&self) -> u32 {
        self.v.total_val()
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let bin = PyBytes::new(py, &slf.borrow().v.to_binary());
        let dict = slf.getattr("__dict__")?;
        PyTuple::new(
            py,
            [
                slf.get_type().into_any(),
                PyTuple::new(py, [bin])?.into_any(),
                PyTuple::new(py, [dict])?.into_any(),
            ],
        )
    }
    fn __getinitargs__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, [PyBytes::new(py, &self.v.to_binary())])
    }
    fn __getstate__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(slf.py(), [slf.getattr("__dict__")?])
    }
    fn __setstate__(slf: &Bound<'_, Self>, state: &Bound<'_, PyTuple>) -> PyResult<()> {
        if state.len() != 1 {
            return Err(PyValueError::new_err(format!(
                "expected 1-item tuple in call to __setstate__; got {state}"
            )));
        }
        slf.getattr("__dict__")?
            .call_method1("update", (state.get_item(0)?,))?;
        Ok(())
    }
}

#[pyclass(
    name = "RealValueVect",
    module = "rdoxide.DataStructs.cDataStructs",
    dict,
    sequence,
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PyRealValueVect {
    pub v: RealValueVect,
}

#[pymethods]
impl PyRealValueVect {
    /// `RealValueVect(length)` or `RealValueVect(pkl)`
    #[new]
    fn new(arg: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(n) = arg.extract::<u32>() {
            return Ok(PyRealValueVect {
                v: RealValueVect::new(n),
            });
        }
        Ok(PyRealValueVect {
            v: RealValueVect::from_binary(&pkl_bytes(arg)?).py()?,
        })
    }
    #[classattr]
    fn __doc__() -> &'static str {
        RVV_DOC
    }
    /// Get the number of entries in the vector
    fn __len__(&self) -> usize {
        self.v.length() as usize
    }
    /// Set the value at a specified location
    fn __setitem__(&mut self, i: u32, val: f64) -> PyResult<()> {
        self.v.set_val(i, val).py()
    }
    /// Get the value at a specified location
    fn __getitem__(&self, i: u32) -> PyResult<f64> {
        self.v.get_val(i).py()
    }
    fn __and__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        self.op(&o, cxx_min)
    }
    fn __or__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        self.op(&o, cxx_max)
    }
    fn __add__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        self.op(&o, |a, b| a + b)
    }
    fn __sub__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
        self.op(&o, |a, b| a - b)
    }
    fn __iadd__(&mut self, o: PyRef<'_, Self>) -> PyResult<()> {
        self.v.apply(&o.v, |a, b| a + b).py()
    }
    fn __isub__(&mut self, o: PyRef<'_, Self>) -> PyResult<()> {
        self.v.apply(&o.v, |a, b| a - b).py()
    }
    /// Get the sum of the values in the vector, basically L1 norm
    #[pyo3(name = "GetTotalVal")]
    fn get_total_val(&self) -> f64 {
        self.v.total_val()
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let bin = PyBytes::new(py, &slf.borrow().v.to_binary());
        let dict = slf.getattr("__dict__")?;
        PyTuple::new(
            py,
            [
                slf.get_type().into_any(),
                PyTuple::new(py, [bin])?.into_any(),
                PyTuple::new(py, [dict])?.into_any(),
            ],
        )
    }
    fn __getinitargs__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, [PyBytes::new(py, &self.v.to_binary())])
    }
    fn __getstate__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(slf.py(), [slf.getattr("__dict__")?])
    }
    fn __setstate__(slf: &Bound<'_, Self>, state: &Bound<'_, PyTuple>) -> PyResult<()> {
        if state.len() != 1 {
            return Err(PyValueError::new_err(format!(
                "expected 1-item tuple in call to __setstate__; got {state}"
            )));
        }
        slf.getattr("__dict__")?
            .call_method1("update", (state.get_item(0)?,))?;
        Ok(())
    }
}

impl PyRealValueVect {
    fn op(&self, o: &Self, f: impl Fn(f64, f64) -> f64) -> PyResult<Self> {
        let mut v = self.v.clone();
        v.apply(&o.v, f).py()?;
        Ok(PyRealValueVect { v })
    }
}

/// Compute the distance between two discrete vector values
#[pyfunction]
#[pyo3(name = "ComputeL1Norm", signature = (v1, v2))]
fn compute_l1_norm(
    py: Python<'_>,
    v1: &Bound<'_, PyAny>,
    v2: &Bound<'_, PyAny>,
) -> PyResult<Py<PyAny>> {
    if let (Ok(a), Ok(b)) = (
        v1.extract::<PyRef<'_, PyDiscreteValueVect>>(),
        v2.extract::<PyRef<'_, PyDiscreteValueVect>>(),
    ) {
        return Ok(compute_l1_norm_discrete(&a.v, &b.v)
            .py()?
            .into_pyobject(py)?
            .into_any()
            .unbind());
    }
    let a = v1.extract::<PyRef<'_, PyRealValueVect>>()?;
    let b = v2.extract::<PyRef<'_, PyRealValueVect>>()?;
    Ok(compute_l1_norm_real(&a.v, &b.v)
        .py()?
        .into_pyobject(py)?
        .into_any()
        .unbind())
}

/// Values for `ConvertToNumpyArray`; `None` for other types.
pub fn numpy_values(o: &Bound<'_, PyAny>) -> Option<Vec<Py<PyAny>>> {
    let py = o.py();
    if let Ok(d) = o.extract::<PyRef<'_, PyDiscreteValueVect>>() {
        return (0..d.v.length())
            .map(|i| d.v.get_val(i).unwrap_or(0).into_py_any(py))
            .collect::<PyResult<_>>()
            .ok();
    }
    if let Ok(r) = o.extract::<PyRef<'_, PyRealValueVect>>() {
        return r
            .v
            .data()
            .iter()
            .map(|&x| x.into_py_any(py))
            .collect::<PyResult<_>>()
            .ok();
    }
    None
}

pub fn add_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum_exported(
        m,
        "DiscreteValueType",
        &[
            ("ONEBITVALUE", 0),
            ("TWOBITVALUE", 1),
            ("FOURBITVALUE", 2),
            ("EIGHTBITVALUE", 3),
            ("SIXTEENBITVALUE", 4),
        ],
        None,
    )?;
    m.add_class::<PyDiscreteValueVect>()?;
    m.add_class::<PyRealValueVect>()?;
    m.add_function(wrap_pyfunction!(compute_l1_norm, m)?)?;
    Ok(())
}
