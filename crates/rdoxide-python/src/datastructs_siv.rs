//! `IntSparseIntVect` / `LongSparseIntVect` / `UIntSparseIntVect` /
//! `ULongSparseIntVect` (upstream `DataStructs/Wrap/SparseIntVect.cpp`).

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList, PyTuple};
use rdoxide::datastructs::sparse_int_vect::{
    IndexKind, SparseIntVect, dice_similarity, tanimoto_similarity, tversky_similarity,
};

use crate::errors::IntoPyResult;

const SIV_DOC: &str = "A container class for storing integer\nvalues within a particular range.\n\nThe length of the vector is set at construction time.\n\nAs you would expect, _SparseIntVects_ support a set of binary operations\nso you can do things like:\n  Arithmetic:\n  siv1 += siv2\n  siv3 = siv1 + siv2\n  siv1 -= siv3\n  siv3 = siv1 - siv2\n  \"Fuzzy\" binary operations:\n  siv3 = siv1 & siv2  the result contains the smallest value in each entry\n  siv3 = siv1 | siv2  the result contains the largest value in each entry\n\nElements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])\n\n";

macro_rules! siv_class {
    ($py:ident, $name:literal, $kind:expr, $idx:ty) => {
        #[pyclass(name = $name, module = "rdkit.DataStructs.cDataStructs", dict, skip_from_py_object)]
        #[derive(Debug)]
        pub struct $py {
            pub v: SparseIntVect,
        }

        impl $py {
            pub(crate) fn wrap(v: SparseIntVect) -> Self {
                $py { v }
            }
        }

        #[pymethods]
        impl $py {
            #[new]
            #[pyo3(signature = (arg))]
            fn new(arg: &Bound<'_, PyAny>) -> PyResult<Self> {
                if let Ok(b) = arg.cast::<PyBytes>() {
                    return Ok(Self::wrap(
                        SparseIntVect::from_binary($kind, b.as_bytes()).py()?,
                    ));
                }
                if let Ok(s) = arg.extract::<String>() {
                    return Ok(Self::wrap(
                        SparseIntVect::from_binary($kind, s.as_bytes()).py()?,
                    ));
                }
                let n: $idx = arg.extract()?;
                Ok(Self::wrap(SparseIntVect::new($kind, n as i128)))
            }
            #[classattr]
            fn __doc__() -> &'static str {
                SIV_DOC
            }
            /// Set the value at a specified location
            fn __setitem__(&mut self, idx: $idx, val: i32) -> PyResult<()> {
                self.v.set_val(idx as i128, val).py()
            }
            /// Get the value at a specified location
            fn __getitem__(&self, idx: $idx) -> PyResult<i32> {
                self.v.get_val(idx as i128).py()
            }
            fn __and__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
                let mut r = self.v.clone();
                r.and_assign(&o.v).py()?;
                Ok(Self::wrap(r))
            }
            fn __or__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
                let mut r = self.v.clone();
                r.or_assign(&o.v).py()?;
                Ok(Self::wrap(r))
            }
            fn __iand__(&mut self, o: PyRef<'_, Self>) -> PyResult<()> {
                self.v.and_assign(&o.v).py()
            }
            fn __ior__(&mut self, o: PyRef<'_, Self>) -> PyResult<()> {
                self.v.or_assign(&o.v).py()
            }
            fn __add__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
                let mut r = self.v.clone();
                r.add_assign(&o.v).py()?;
                Ok(Self::wrap(r))
            }
            fn __sub__(&self, o: PyRef<'_, Self>) -> PyResult<Self> {
                let mut r = self.v.clone();
                r.sub_assign(&o.v).py()?;
                Ok(Self::wrap(r))
            }
            fn __iadd__(&mut self, o: &Bound<'_, PyAny>) -> PyResult<()> {
                if let Ok(o) = o.extract::<PyRef<'_, Self>>() {
                    return self.v.add_assign(&o.v).py();
                }
                self.v.scalar_op('+', o.extract()?).py()
            }
            fn __isub__(&mut self, o: &Bound<'_, PyAny>) -> PyResult<()> {
                if let Ok(o) = o.extract::<PyRef<'_, Self>>() {
                    return self.v.sub_assign(&o.v).py();
                }
                self.v.scalar_op('-', o.extract()?).py()
            }
            fn __imul__(&mut self, s: i32) -> PyResult<()> {
                self.v.scalar_op('*', s).py()
            }
            fn __eq__(&self, o: &Bound<'_, PyAny>) -> bool {
                o.extract::<PyRef<'_, Self>>().is_ok_and(|o| o.v == self.v)
            }
            fn __ne__(&self, o: &Bound<'_, PyAny>) -> bool {
                !self.__eq__(o)
            }
            /// Get the sum of the values in the vector, basically L1 norm
            #[pyo3(name = "GetTotalVal", signature = (useAbs = false))]
            #[allow(non_snake_case)]
            fn get_total_val(&self, useAbs: bool) -> i32 {
                self.v.total_val(useAbs)
            }
            /// Returns the length of the vector
            #[pyo3(name = "GetLength")]
            fn get_length(&self) -> $idx {
                self.v.length() as $idx
            }
            /// returns a binary (pickle) representation of the vector
            #[pyo3(name = "ToBinary")]
            fn to_binary<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
                PyBytes::new(py, &self.v.to_binary())
            }
            /// update the vector based on the values in the list or tuple
            #[pyo3(name = "UpdateFromSequence", signature = (seq))]
            fn update_from_sequence(&mut self, seq: Vec<$idx>) -> PyResult<()> {
                let s: Vec<i128> = seq.into_iter().map(|i| i as i128).collect();
                self.v.update_from_sequence(&s).py()
            }
            /// returns a dictionary of the nonzero elements
            #[pyo3(name = "GetNonzeroElements")]
            fn get_nonzero_elements<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                let d = PyDict::new(py);
                for (&k, &v) in self.v.nonzero_elements() {
                    d.set_item(k as $idx, v)?;
                }
                Ok(d)
            }
            /// Return the SparseIntVect as a python list
            #[pyo3(name = "ToList")]
            fn to_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
                let mut l = vec![0i32; self.v.length() as usize];
                for (&k, &v) in self.v.nonzero_elements() {
                    l[k as usize] = v;
                }
                PyList::new(py, l)
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
    };
}

siv_class!(PyIntSIV, "IntSparseIntVect", IndexKind::I32, i32);
siv_class!(PyLongSIV, "LongSparseIntVect", IndexKind::I64, i64);
siv_class!(PyUIntSIV, "UIntSparseIntVect", IndexKind::U32, u32);
siv_class!(PyULongSIV, "ULongSparseIntVect", IndexKind::U64, u64);

/// Extracts a copy of the vector if `o` is any SparseIntVect class.
pub fn get_siv(o: &Bound<'_, PyAny>) -> Option<SparseIntVect> {
    if let Ok(x) = o.extract::<PyRef<'_, PyIntSIV>>() {
        return Some(x.v.clone());
    }
    if let Ok(x) = o.extract::<PyRef<'_, PyLongSIV>>() {
        return Some(x.v.clone());
    }
    if let Ok(x) = o.extract::<PyRef<'_, PyUIntSIV>>() {
        return Some(x.v.clone());
    }
    if let Ok(x) = o.extract::<PyRef<'_, PyULongSIV>>() {
        return Some(x.v.clone());
    }
    None
}

fn pair(
    a: &Bound<'_, PyAny>,
    b: &Bound<'_, PyAny>,
) -> PyResult<Option<(SparseIntVect, SparseIntVect)>> {
    let Some(x) = get_siv(a) else { return Ok(None) };
    match get_siv(b) {
        Some(y) if y.kind() == x.kind() => Ok(Some((x, y))),
        _ => Err(PyTypeError::new_err(format!(
            "Python argument types did not match C++ signature: got {} and {}",
            a.get_type().name()?,
            b.get_type().name()?
        ))),
    }
}

/// Which SIV similarity a generic similarity function maps to.
#[derive(Clone, Copy)]
pub enum SivMetric {
    Dice,
    Tanimoto,
    Tversky(f64, f64),
}

fn sim(
    m: SivMetric,
    a: &SparseIntVect,
    b: &SparseIntVect,
    dist: bool,
    bounds: f64,
) -> rdoxide::Result<f64> {
    match m {
        SivMetric::Dice => dice_similarity(a, b, dist, bounds),
        SivMetric::Tanimoto => tanimoto_similarity(a, b, dist),
        SivMetric::Tversky(x, y) => tversky_similarity(a, b, x, y, dist),
    }
}

/// The SparseIntVect overload of `DiceSimilarity`/`TanimotoSimilarity`/`TverskySimilarity`;
/// `None` if the arguments are not SparseIntVects.
pub fn try_similarity(
    m: SivMetric,
    a: &Bound<'_, PyAny>,
    b: &Bound<'_, PyAny>,
    dist: bool,
    bounds: f64,
) -> PyResult<Option<f64>> {
    match pair(a, b)? {
        Some((x, y)) => Ok(Some(sim(m, &x, &y, dist, bounds).py()?)),
        None => Ok(None),
    }
}

/// The SparseIntVect overload of the Bulk* functions.
pub fn try_bulk(
    m: SivMetric,
    a: &Bound<'_, PyAny>,
    list: &Bound<'_, PyAny>,
    dist: bool,
) -> PyResult<Option<Vec<f64>>> {
    let Some(x) = get_siv(a) else { return Ok(None) };
    let mut res = Vec::new();
    for item in list.try_iter()? {
        let item = item?;
        let Some((_, y)) = pair(a, &item)? else {
            unreachable!()
        };
        res.push(sim(m, &x, &y, dist, 0.0).py()?);
    }
    Ok(Some(res))
}

pub fn add_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyIntSIV>()?;
    m.add_class::<PyLongSIV>()?;
    m.add_class::<PyUIntSIV>()?;
    m.add_class::<PyULongSIV>()?;
    Ok(())
}
