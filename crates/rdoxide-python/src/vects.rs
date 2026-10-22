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
        #[pyclass(name = $name, module = "rdoxide.rdBase", sequence, skip_from_py_object)]
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

macro_rules! nested_list {
    ($py:ident, $inner:ident, $t:ty) => {
        #[pyclass(module = "rdoxide.rdBase", sequence, skip_from_py_object)]
        #[derive(Debug, Clone, Default)]
        pub struct $py {
            pub v: Vec<Vec<$t>>,
        }

        impl $py {
            pub fn new_py(py: Python<'_>, v: Vec<Vec<$t>>) -> PyResult<Py<PyAny>> {
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
                $inner::new_py(py, self.v[i].clone())
            }
            fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<PyAny>> {
                let py = slf.py();
                let items = slf
                    .v
                    .iter()
                    .map(|p| $inner::new_py(py, p.clone()))
                    .collect::<PyResult<Vec<_>>>()?;
                Ok(PyList::new(py, items)?
                    .into_any()
                    .try_iter()?
                    .into_any()
                    .unbind())
            }
            fn __contains__(&self, val: &Bound<'_, PyAny>) -> bool {
                val.extract::<PyRef<'_, $inner>>()
                    .is_ok_and(|x| self.v.contains(&x.v))
                    || val.extract::<Vec<$t>>().is_ok_and(|x| self.v.contains(&x))
            }
        }
    };
}

/// Adds the `list_indexing_suite` item assignment error of `std::list`.
macro_rules! list_setitem {
    ($py:ident) => {
        #[pymethods]
        impl $py {
            fn __setitem__(&mut self, _idx: isize, _val: &Bound<'_, PyAny>) -> PyResult<()> {
                Err(PyTypeError::new_err(
                    "'list' object does not support item assignment",
                ))
            }
        }
    };
}

macro_rules! nested_vect_mut {
    ($py:ident, $t:ty) => {
        #[pymethods]
        impl $py {
            fn __setitem__(&mut self, idx: isize, val: Vec<$t>) -> PyResult<()> {
                let i = norm_index(idx, self.v.len())?;
                self.v[i] = val;
                Ok(())
            }
            fn __delitem__(&mut self, idx: isize) -> PyResult<()> {
                let i = norm_index(idx, self.v.len())?;
                self.v.remove(i);
                Ok(())
            }
            fn append(&mut self, val: Vec<$t>) {
                self.v.push(val);
            }
        }
    };
}

#[cfg(target_os = "macos")]
const LIST_VECTI: &str = "_listNSt3__16vectorIiNS_9allocatorIiEEEE";
#[cfg(not(target_os = "macos"))]
const LIST_VECTI: &str = "_listSt6vectorIiSaIiEE";
#[cfg(target_os = "macos")]
const LIST_VECTJ: &str = "_listNSt3__16vectorIjNS_9allocatorIjEEEE";
#[cfg(not(target_os = "macos"))]
const LIST_VECTJ: &str = "_listSt6vectorIjSaIjEE";
#[cfg(target_os = "macos")]
const VECT_VECTI: &str = "_vectNSt3__16vectorIiNS_9allocatorIiEEEE";
#[cfg(not(target_os = "macos"))]
const VECT_VECTI: &str = "_vectSt6vectorIiSaIiEE";
#[cfg(target_os = "macos")]
const VECT_VECTJ: &str = "_vectNSt3__16vectorIjNS_9allocatorIjEEEE";
#[cfg(not(target_os = "macos"))]
const VECT_VECTJ: &str = "_vectSt6vectorIjSaIjEE";
#[cfg(target_os = "macos")]
const VECT_VECTD: &str = "_vectNSt3__16vectorIdNS_9allocatorIdEEEE";
#[cfg(not(target_os = "macos"))]
const VECT_VECTD: &str = "_vectSt6vectorIdSaIdEE";

// `std::list<std::vector<int>>` (e.g. `PATH_LIST`) and `std::list<std::vector<unsigned>>`:
// items come back as copies of the flat vector types
nested_list!(PyListVectI, PyVectI, i32);
list_setitem!(PyListVectI);
nested_list!(PyListVectJ, PyVectJ, u32);
list_setitem!(PyListVectJ);
// `std::vector<std::vector<T>>` (e.g. the clusters of a picker)
nested_list!(PyVectVectI, PyVectI, i32);
nested_vect_mut!(PyVectVectI, i32);
nested_list!(PyVectVectJ, PyVectJ, u32);
nested_vect_mut!(PyVectVectJ, u32);
nested_list!(PyVectVectD, PyVectD, f64);
nested_vect_mut!(PyVectVectD, f64);

fn add_named<T: pyo3::PyClass>(
    m: &Bound<'_, PyModule>,
    rust_name: &str,
    name: &str,
) -> PyResult<()> {
    let py = m.py();
    let cls = py.get_type::<T>();
    cls.setattr("__name__", name)?;
    cls.setattr("__qualname__", name)?;
    m.add(name, cls)?;
    let _ = rust_name;
    Ok(())
}

pub fn add_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyVectI>()?;
    m.add_class::<PyVectJ>()?;
    m.add_class::<PyVectD>()?;
    m.add_class::<PyVectStr>()?;
    add_named::<PyListVectI>(m, "PyListVectI", LIST_VECTI)?;
    add_named::<PyListVectJ>(m, "PyListVectJ", LIST_VECTJ)?;
    add_named::<PyVectVectI>(m, "PyVectVectI", VECT_VECTI)?;
    add_named::<PyVectVectJ>(m, "PyVectVectJ", VECT_VECTJ)?;
    add_named::<PyVectVectD>(m, "PyVectVectD", VECT_VECTD)?;
    Ok(())
}
