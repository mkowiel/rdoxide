//! `rdBase` container classes: Boost.Python `vector_indexing_suite` /
//! `list_indexing_suite` wrappers of `std::vector<int>` (`_vecti`),
//! `std::vector<unsigned>` (`_vectj`), `std::vector<double>` (`_vectd`) and
//! sequences of int vectors (paths).
//!
//! Like upstream they are mutable sequences without value equality. The
//! class names of the nested types are C++-mangled upstream and so differ
//! between standard libraries; the libc++ names are used on macOS and the
//! libstdc++ ones elsewhere.

use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PySlice};

/// Normalizes a (possibly negative) index like the indexing suites do.
fn norm_index(i: isize, len: usize) -> PyResult<usize> {
    let j = if i < 0 { i + len as isize } else { i };
    if j < 0 || j as usize >= len {
        return Err(PyIndexError::new_err("Index out of range"));
    }
    Ok(j as usize)
}

macro_rules! flat_vect {
    ($py:ident, $name:literal, $t:ty) => {
        #[pyclass(name = $name, module = "rdkit.rdBase", sequence, skip_from_py_object)]
        #[derive(Debug, Clone, Default)]
        pub struct $py {
            pub v: Vec<$t>,
        }

        impl $py {
            pub fn new_py(py: Python<'_>, v: Vec<$t>) -> PyResult<Py<PyAny>> {
                Ok(Py::new(py, $py { v })?.into_any())
            }
        }

        #[pymethods]
        impl $py {
            #[new]
            fn new() -> Self {
                Self::default()
            }
            fn __len__(&self) -> usize {
                self.v.len()
            }
            fn __getitem__(&self, py: Python<'_>, idx: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
                if let Ok(s) = idx.cast::<PySlice>() {
                    let ind = s.indices(self.v.len() as isize)?;
                    let mut out = Vec::new();
                    let mut i = ind.start;
                    while (ind.step > 0 && i < ind.stop) || (ind.step < 0 && i > ind.stop) {
                        out.push(self.v[i as usize].clone());
                        i += ind.step;
                    }
                    return Self::new_py(py, out);
                }
                let i = norm_index(idx.extract::<isize>()?, self.v.len())?;
                Ok(self.v[i].clone().into_pyobject(py)?.into_any().unbind())
            }
            fn __setitem__(&mut self, idx: isize, val: $t) -> PyResult<()> {
                let i = norm_index(idx, self.v.len())?;
                self.v[i] = val;
                Ok(())
            }
            fn __delitem__(&mut self, idx: isize) -> PyResult<()> {
                let i = norm_index(idx, self.v.len())?;
                self.v.remove(i);
                Ok(())
            }
            fn __contains__(&self, val: &Bound<'_, PyAny>) -> bool {
                val.extract::<$t>().is_ok_and(|x| self.v.contains(&x))
            }
            fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<PyAny>> {
                let py = slf.py();
                Ok(PyList::new(py, slf.v.iter().cloned())?
                    .into_any()
                    .try_iter()?
                    .into_any()
                    .unbind())
            }
            fn append(&mut self, val: $t) {
                self.v.push(val);
            }
            fn extend(&mut self, vals: &Bound<'_, PyAny>) -> PyResult<()> {
                for x in vals.try_iter()? {
                    self.v.push(x?.extract()?);
                }
                Ok(())
            }
        }
    };
}

flat_vect!(PyVectI, "_vecti", i32);
flat_vect!(PyVectJ, "_vectj", u32);
flat_vect!(PyVectD, "_vectd", f64);
// `std::vector<std::string>`; the class name is C++-mangled like the nested ones below
#[cfg(target_os = "macos")]
flat_vect!(
    PyVectStr,
    "_vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE",
    String
);
#[cfg(not(target_os = "macos"))]
flat_vect!(
    PyVectStr,
    "_vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE",
    String
);

#[cfg(target_os = "macos")]
const LIST_VECTI: &str = "_listNSt3__16vectorIiNS_9allocatorIiEEEE";
#[cfg(not(target_os = "macos"))]
const LIST_VECTI: &str = "_listSt6vectorIiSaIiEE";

/// A sequence of int vectors (`std::list<std::vector<int>>`, e.g. `PATH_LIST`).
/// Items come back as `_vecti` copies.
#[pyclass(module = "rdkit.rdBase", sequence, skip_from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct PyListVectI {
    pub v: Vec<Vec<i32>>,
}

impl PyListVectI {
    pub fn new_py(py: Python<'_>, v: Vec<Vec<i32>>) -> PyResult<Py<PyAny>> {
        Ok(Py::new(py, PyListVectI { v })?.into_any())
    }
}

#[pymethods]
impl PyListVectI {
    fn __len__(&self) -> usize {
        self.v.len()
    }
    fn __getitem__(&self, py: Python<'_>, idx: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(s) = idx.cast::<PySlice>() {
            let ind = s.indices(self.v.len() as isize)?;
            let mut out = Vec::new();
            let mut i = ind.start;
            while (ind.step > 0 && i < ind.stop) || (ind.step < 0 && i > ind.stop) {
                out.push(self.v[i as usize].clone());
                i += ind.step;
            }
            return Self::new_py(py, out);
        }
        let i = norm_index(idx.extract::<isize>()?, self.v.len())?;
        PyVectI::new_py(py, self.v[i].clone())
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        let items = slf
            .v
            .iter()
            .map(|p| PyVectI::new_py(py, p.clone()))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, items)?
            .into_any()
            .try_iter()?
            .into_any()
            .unbind())
    }
    fn __contains__(&self, val: &Bound<'_, PyAny>) -> bool {
        val.extract::<PyRef<'_, PyVectI>>()
            .is_ok_and(|x| self.v.contains(&x.v))
            || val.extract::<Vec<i32>>().is_ok_and(|x| self.v.contains(&x))
    }
    fn __setitem__(&mut self, _idx: isize, _val: &Bound<'_, PyAny>) -> PyResult<()> {
        Err(PyTypeError::new_err(
            "'list' object does not support item assignment",
        ))
    }
}

#[cfg(target_os = "macos")]
const VECT_VECTI: &str = "_vectNSt3__16vectorIiNS_9allocatorIiEEEE";
#[cfg(not(target_os = "macos"))]
const VECT_VECTI: &str = "_vectSt6vectorIiSaIiEE";

/// A vector of int vectors (`std::vector<std::vector<int>>`, e.g. the
/// clusters of a picker). Items come back as `_vecti` copies.
#[pyclass(module = "rdkit.rdBase", sequence, skip_from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct PyVectVectI {
    pub v: Vec<Vec<i32>>,
}

impl PyVectVectI {
    pub fn new_py(py: Python<'_>, v: Vec<Vec<i32>>) -> PyResult<Py<PyAny>> {
        Ok(Py::new(py, PyVectVectI { v })?.into_any())
    }
}

#[pymethods]
impl PyVectVectI {
    fn __len__(&self) -> usize {
        self.v.len()
    }
    fn __getitem__(&self, py: Python<'_>, idx: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(s) = idx.cast::<PySlice>() {
            let ind = s.indices(self.v.len() as isize)?;
            let mut out = Vec::new();
            let mut i = ind.start;
            while (ind.step > 0 && i < ind.stop) || (ind.step < 0 && i > ind.stop) {
                out.push(self.v[i as usize].clone());
                i += ind.step;
            }
            return Self::new_py(py, out);
        }
        let i = norm_index(idx.extract::<isize>()?, self.v.len())?;
        PyVectI::new_py(py, self.v[i].clone())
    }
    fn __setitem__(&mut self, idx: isize, val: Vec<i32>) -> PyResult<()> {
        let i = norm_index(idx, self.v.len())?;
        self.v[i] = val;
        Ok(())
    }
    fn __delitem__(&mut self, idx: isize) -> PyResult<()> {
        let i = norm_index(idx, self.v.len())?;
        self.v.remove(i);
        Ok(())
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        let items = slf
            .v
            .iter()
            .map(|p| PyVectI::new_py(py, p.clone()))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, items)?
            .into_any()
            .try_iter()?
            .into_any()
            .unbind())
    }
    fn __contains__(&self, val: &Bound<'_, PyAny>) -> bool {
        val.extract::<PyRef<'_, PyVectI>>()
            .is_ok_and(|x| self.v.contains(&x.v))
            || val.extract::<Vec<i32>>().is_ok_and(|x| self.v.contains(&x))
    }
    fn append(&mut self, val: Vec<i32>) {
        self.v.push(val);
    }
}

pub fn add_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyVectI>()?;
    m.add_class::<PyVectJ>()?;
    m.add_class::<PyVectD>()?;
    m.add_class::<PyVectStr>()?;
    m.add_class::<PyListVectI>()?;
    // the mangled name upstream exposes for this class
    let cls = m.getattr("PyListVectI")?;
    cls.setattr("__name__", LIST_VECTI)?;
    cls.setattr("__qualname__", LIST_VECTI)?;
    m.add(LIST_VECTI, &cls)?;
    m.delattr("PyListVectI")?;
    m.add_class::<PyVectVectI>()?;
    let cls = m.getattr("PyVectVectI")?;
    cls.setattr("__name__", VECT_VECTI)?;
    cls.setattr("__qualname__", VECT_VECTI)?;
    m.add(VECT_VECTI, &cls)?;
    m.delattr("PyVectVectI")?;
    Ok(())
}
