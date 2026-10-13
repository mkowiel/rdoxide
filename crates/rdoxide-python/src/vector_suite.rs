//! Boost.Python's `vector_indexing_suite` for a `Vec<T>` held by a pyclass:
//! `__len__`, `__getitem__`/`__setitem__`/`__delitem__` with integers and
//! step-less slices, `__contains__`, `__iter__`, `append`, `extend`, with
//! upstream's error messages. The pyclass forwards to [`VectorSuite`].

use pyo3::exceptions::{PyAttributeError, PyIndexError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PySlice, PySliceMethods};

/// A position in the vector: an element or a `[start, stop)` range.
pub enum Index {
    Item(usize),
    Range(usize, usize),
}

/// The element conversions of one vector class.
pub trait VectorSuite {
    type Item: Clone + PartialEq;

    fn items(&self) -> &Vec<Self::Item>;
    fn items_mut(&mut self) -> &mut Vec<Self::Item>;
    /// The Python object for an element.
    fn to_py<'py>(py: Python<'py>, item: &Self::Item) -> PyResult<Bound<'py, PyAny>>;
    /// The element held by `obj`, if it is one.
    fn from_py(obj: &Bound<'_, PyAny>) -> Option<Self::Item>;
    /// A new vector object of the same class.
    fn new_vector<'py>(py: Python<'py>, items: Vec<Self::Item>) -> PyResult<Bound<'py, PyAny>>;

    fn suite_index(&self, key: &Bound<'_, PyAny>) -> PyResult<Index> {
        let len = self.items().len() as isize;
        if let Ok(slice) = key.cast::<PySlice>() {
            if !slice.getattr("step")?.is_none() {
                return Err(PyIndexError::new_err("slice step size not supported."));
            }
            let ind = slice.indices(len)?;
            let (start, stop) = (ind.start.max(0) as usize, ind.stop.max(0) as usize);
            return Ok(Index::Range(start, stop.max(start)));
        }
        let i: isize = key
            .extract()
            .map_err(|_| PyTypeError::new_err("Invalid index type"))?;
        let i = if i < 0 { i + len } else { i };
        if i < 0 || i >= len {
            return Err(PyIndexError::new_err("Index out of range"));
        }
        Ok(Index::Item(i as usize))
    }

    fn suite_len(&self) -> usize {
        self.items().len()
    }

    fn suite_getitem<'py>(
        &self,
        py: Python<'py>,
        key: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        match self.suite_index(key)? {
            Index::Item(i) => Self::to_py(py, &self.items()[i]),
            Index::Range(a, b) => Self::new_vector(py, self.items()[a..b].to_vec()),
        }
    }

    /// The elements of `value`: one element, or an iterable of them.
    fn elements_of(value: &Bound<'_, PyAny>, err: &str) -> PyResult<Vec<Self::Item>> {
        if let Some(item) = Self::from_py(value) {
            return Ok(vec![item]);
        }
        let iter = value.try_iter().map_err(|_| {
            PyAttributeError::new_err(format!(
                "'{}' object has no attribute '__iter__'",
                value
                    .get_type()
                    .name()
                    .map(|n| n.to_string())
                    .unwrap_or_default()
            ))
        })?;
        iter.map(|x| Self::from_py(&x?).ok_or_else(|| PyTypeError::new_err(err.to_string())))
            .collect()
    }

    fn suite_setitem(&mut self, key: &Bound<'_, PyAny>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        match self.suite_index(key)? {
            Index::Item(i) => {
                let item = Self::from_py(value)
                    .ok_or_else(|| PyTypeError::new_err("Invalid assignment"))?;
                self.items_mut()[i] = item;
            }
            Index::Range(a, b) => {
                let items = Self::elements_of(value, "Invalid assignment")?;
                self.items_mut().splice(a..b, items);
            }
        }
        Ok(())
    }

    fn suite_delitem(&mut self, key: &Bound<'_, PyAny>) -> PyResult<()> {
        match self.suite_index(key)? {
            Index::Item(i) => {
                self.items_mut().remove(i);
            }
            Index::Range(a, b) => {
                self.items_mut().drain(a..b);
            }
        }
        Ok(())
    }

    fn suite_contains(&self, value: &Bound<'_, PyAny>) -> bool {
        Self::from_py(value).is_some_and(|v| self.items().contains(&v))
    }

    fn suite_iter<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let objs = self
            .items()
            .iter()
            .map(|x| Self::to_py(py, x))
            .collect::<PyResult<Vec<_>>>()?;
        PyList::new(py, objs)?
            .into_any()
            .try_iter()
            .map(Bound::into_any)
    }

    fn suite_append(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let item = Self::from_py(value)
            .ok_or_else(|| PyTypeError::new_err("Attempting to append an invalid type"))?;
        self.items_mut().push(item);
        Ok(())
    }

    fn suite_extend(&mut self, values: &Bound<'_, PyAny>) -> PyResult<()> {
        let iter = values.try_iter().map_err(|_| {
            PyAttributeError::new_err(format!(
                "'{}' object has no attribute '__iter__'",
                values
                    .get_type()
                    .name()
                    .map(|n| n.to_string())
                    .unwrap_or_default()
            ))
        })?;
        let items = iter
            .map(|x| {
                Self::from_py(&x?).ok_or_else(|| PyTypeError::new_err("Incompatible Data Type"))
            })
            .collect::<PyResult<Vec<_>>>()?;
        self.items_mut().extend(items);
        Ok(())
    }
}

/// The Python methods of a [`VectorSuite`] pyclass.
macro_rules! vector_suite_methods {
    ($cls:ty) => {
        #[pymethods]
        impl $cls {
            fn __len__(&self) -> usize {
                $crate::vector_suite::VectorSuite::suite_len(self)
            }
            fn __getitem__<'py>(
                &self,
                py: Python<'py>,
                key: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                $crate::vector_suite::VectorSuite::suite_getitem(self, py, key)
            }
            fn __setitem__(
                &mut self,
                key: &Bound<'_, PyAny>,
                value: &Bound<'_, PyAny>,
            ) -> PyResult<()> {
                $crate::vector_suite::VectorSuite::suite_setitem(self, key, value)
            }
            fn __delitem__(&mut self, key: &Bound<'_, PyAny>) -> PyResult<()> {
                $crate::vector_suite::VectorSuite::suite_delitem(self, key)
            }
            fn __contains__(&self, value: &Bound<'_, PyAny>) -> bool {
                $crate::vector_suite::VectorSuite::suite_contains(self, value)
            }
            fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
                $crate::vector_suite::VectorSuite::suite_iter(self, py)
            }
            fn append(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
                $crate::vector_suite::VectorSuite::suite_append(self, value)
            }
            fn extend(&mut self, values: &Bound<'_, PyAny>) -> PyResult<()> {
                $crate::vector_suite::VectorSuite::suite_extend(self, values)
            }
        }
    };
}
pub(crate) use vector_suite_methods;

/// A class created only by the module (Boost.Python's `no_init`): calling
/// it raises upstream's RuntimeError.
macro_rules! no_python_init {
    ($($cls:ty),* $(,)?) => {
        $(
            #[pymethods]
            impl $cls {
                #[new]
                #[pyo3(signature = (*_args, **_kwargs))]
                fn __no_init(
                    _args: &Bound<'_, pyo3::types::PyTuple>,
                    _kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
                ) -> PyResult<Self> {
                    Err(pyo3::exceptions::PyRuntimeError::new_err(
                        "This class cannot be instantiated from Python",
                    ))
                }
            }
        )*
    };
}
pub(crate) use no_python_init;
