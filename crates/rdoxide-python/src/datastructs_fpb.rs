//! `FPBReader` and `MultiFPBReader` (upstream `DataStructs/Wrap/wrap_FPB.cpp`).

use std::sync::{Arc, RwLock};

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple};
use rdoxide::datastructs::fpb_reader::{FpbReader, Neighbor};
use rdoxide::datastructs::multi_fpb_reader::{MultiFpbReader, SharedFpbReader};

use crate::datastructs::{PyExplicitBitVect, bytes_arg};
use crate::errors::to_py;

/// Like `IntoPyResult`, except that file errors are `RuntimeError`s: the
/// upstream DataStructs module registers no translator for
/// `BadFileException`.
trait FpbResult<T> {
    fn py(self) -> PyResult<T>;
}

impl<T> FpbResult<T> for rdoxide::Result<T> {
    fn py(self) -> PyResult<T> {
        self.map_err(|e| match e {
            rdoxide::RDKitError::BadFile(msg) => PyRuntimeError::new_err(msg),
            other => to_py(other),
        })
    }
}

fn neighbors<'py>(py: Python<'py>, hits: Vec<Neighbor>) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, hits)
}

/// A class for reading and searching FPB files from Andrew Dalke's chemfp.
/// Note that this functionality is still experimental and the API may
/// change in future releases.
#[pyclass(name = "FPBReader", module = "rdoxide.DataStructs.cDataStructs", dict)]
#[derive(Debug)]
pub struct PyFpbReader {
    reader: SharedFpbReader,
}

impl PyFpbReader {
    fn read<T>(&self, f: impl FnOnce(&FpbReader) -> rdoxide::Result<T>) -> PyResult<T> {
        f(&self.reader.read().unwrap_or_else(|e| e.into_inner())).py()
    }
}

#[pymethods]
impl PyFpbReader {
    #[new]
    #[pyo3(signature = (filename, lazy = false))]
    fn new(filename: &str, lazy: bool) -> PyResult<Self> {
        let reader = FpbReader::open(filename, lazy).py()?;
        Ok(PyFpbReader {
            reader: Arc::new(RwLock::new(reader)),
        })
    }

    /// Read the fingerprints from the file. This can take a while.
    #[pyo3(name = "Init")]
    fn init(&self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| {
            self.reader
                .write()
                .unwrap_or_else(|e| e.into_inner())
                .init()
        })
        .py()
    }

    fn __len__(&self) -> PyResult<usize> {
        self.read(|r| r.len()).map(|n| n as usize)
    }

    fn __getitem__(&self, which: u32) -> PyResult<(PyExplicitBitVect, String)> {
        let (fp, id) = self.read(|r| r.get(which))?;
        Ok((PyExplicitBitVect::wrap(fp), id))
    }

    /// returns the number of bits in a fingerprint
    #[pyo3(name = "GetNumBits")]
    fn get_num_bits(&self) -> PyResult<u32> {
        self.read(|r| r.n_bits())
    }

    /// returns a particular fingerprint as an ExplicitBitVect
    #[pyo3(name = "GetFP")]
    fn get_fp(&self, idx: u32) -> PyResult<PyExplicitBitVect> {
        self.read(|r| r.get_fp(idx)).map(PyExplicitBitVect::wrap)
    }

    /// returns a particular fingerprint as bytes
    #[pyo3(name = "GetBytes")]
    fn get_bytes<'py>(&self, py: Python<'py>, which: u32) -> PyResult<Bound<'py, PyBytes>> {
        let (bytes, n_bits) = self.read(|r| Ok((r.get_bytes(which)?, r.n_bits()?)))?;
        Ok(PyBytes::new(py, &bytes[..(n_bits / 8) as usize]))
    }

    /// returns the id of a particular fingerprint
    #[pyo3(name = "GetId")]
    fn get_id(&self, idx: u32) -> PyResult<String> {
        self.read(|r| r.get_id(idx))
    }

    /// return the tanimoto similarity of a particular fingerprint to the bytes provided
    #[pyo3(name = "GetTanimoto")]
    fn get_tanimoto(&self, which: u32, bytes: &Bound<'_, PyAny>) -> PyResult<f64> {
        let bv = bytes_arg(bytes)?;
        self.read(|r| r.get_tanimoto(which, &bv))
    }

    /// returns tanimoto similarities to and indices of all neighbors above the specified threshold
    #[pyo3(name = "GetTanimotoNeighbors", signature = (bv, threshold = 0.7))]
    fn get_tanimoto_neighbors<'py>(
        &self,
        py: Python<'py>,
        bv: &Bound<'py, PyAny>,
        threshold: f64,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let bv = bytes_arg(bv)?;
        let hits = py.detach(|| self.read(|r| r.get_tanimoto_neighbors(&bv, threshold, true)))?;
        neighbors(py, hits)
    }

    /// return the Tverksy similarity of a particular fingerprint to the bytes provided
    #[pyo3(name = "GetTversky")]
    fn get_tversky(&self, which: u32, bytes: &Bound<'_, PyAny>, ca: f64, cb: f64) -> PyResult<f64> {
        let bv = bytes_arg(bytes)?;
        self.read(|r| r.get_tversky(which, &bv, ca, cb))
    }

    /// returns Tversky similarities to and indices of all neighbors above the specified threshold
    #[pyo3(name = "GetTverskyNeighbors", signature = (bv, ca, cb, threshold = 0.7))]
    fn get_tversky_neighbors<'py>(
        &self,
        py: Python<'py>,
        bv: &Bound<'py, PyAny>,
        ca: f64,
        cb: f64,
        threshold: f64,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let bv = bytes_arg(bv)?;
        let hits =
            py.detach(|| self.read(|r| r.get_tversky_neighbors(&bv, ca, cb, threshold, true)))?;
        neighbors(py, hits)
    }

    /// returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)
    #[pyo3(name = "GetContainingNeighbors")]
    fn get_containing_neighbors<'py>(
        &self,
        py: Python<'py>,
        bv: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let bv = bytes_arg(bv)?;
        let hits = py.detach(|| self.read(|r| r.get_containing_neighbors(&bv)))?;
        PyTuple::new(py, hits)
    }
}

/// A class for reading and searching multiple FPB files from Andrew Dalke's chemfp.
/// Note that this functionality is still experimental and the API may
/// change in future releases.
#[pyclass(
    name = "MultiFPBReader",
    module = "rdoxide.DataStructs.cDataStructs",
    dict
)]
#[derive(Debug)]
pub struct PyMultiFpbReader {
    multi: RwLock<MultiFpbReader>,
    /// The Python objects behind the readers, returned by `GetReader`.
    readers: RwLock<Vec<Py<PyFpbReader>>>,
}

impl PyMultiFpbReader {
    fn multi(&self) -> std::sync::RwLockReadGuard<'_, MultiFpbReader> {
        self.multi.read().unwrap_or_else(|e| e.into_inner())
    }
}

#[pymethods]
impl PyMultiFpbReader {
    #[new]
    #[pyo3(signature = (initOnSearch = false))]
    #[allow(non_snake_case)]
    fn new(initOnSearch: bool) -> Self {
        PyMultiFpbReader {
            multi: RwLock::new(MultiFpbReader::new(initOnSearch)),
            readers: RwLock::new(Vec::new()),
        }
    }

    /// Call Init() on each of our children. This can take a while.
    #[pyo3(name = "Init")]
    fn init(&self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.multi.write().unwrap_or_else(|e| e.into_inner()).init())
            .py()
    }

    fn __len__(&self) -> usize {
        self.multi().len()
    }

    /// returns the number of bits in a fingerprint
    #[pyo3(name = "GetNumBits")]
    fn get_num_bits(&self) -> PyResult<u32> {
        self.multi().n_bits().py()
    }

    /// adds an FPBReader to our set of readers
    #[pyo3(name = "AddReader")]
    fn add_reader(&self, py: Python<'_>, rdr: Py<PyFpbReader>) -> PyResult<usize> {
        let shared = Arc::clone(&rdr.borrow(py).reader);
        let count = self
            .multi
            .write()
            .unwrap_or_else(|e| e.into_inner())
            .add_reader(shared)
            .py()?;
        self.readers
            .write()
            .unwrap_or_else(|e| e.into_inner())
            .push(rdr);
        Ok(count)
    }

    /// returns one of our readers
    #[pyo3(name = "GetReader")]
    fn get_reader(&self, py: Python<'_>, which: usize) -> PyResult<Py<PyFpbReader>> {
        // validates the index with upstream's error
        self.multi().reader(which).py()?;
        let readers = self.readers.read().unwrap_or_else(|e| e.into_inner());
        Ok(readers[which].clone_ref(py))
    }

    /// returns tanimoto similarities to and indices of all neighbors above the specified threshold
    #[pyo3(name = "GetTanimotoNeighbors", signature = (bv, threshold = 0.7, numThreads = 0))]
    #[allow(non_snake_case)]
    fn get_tanimoto_neighbors<'py>(
        &self,
        py: Python<'py>,
        bv: &Bound<'py, PyAny>,
        threshold: f64,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let bv = bytes_arg(bv)?;
        let hits = py
            .detach(|| {
                self.multi()
                    .get_tanimoto_neighbors(&bv, threshold, numThreads)
            })
            .py()?;
        PyTuple::new(py, hits)
    }

    /// returns Tversky similarities to and indices of all neighbors above the specified threshold
    #[pyo3(name = "GetTverskyNeighbors", signature = (bv, ca, cb, threshold = 0.7, numThreads = 0))]
    #[allow(non_snake_case)]
    fn get_tversky_neighbors<'py>(
        &self,
        py: Python<'py>,
        bv: &Bound<'py, PyAny>,
        ca: f64,
        cb: f64,
        threshold: f64,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let bv = bytes_arg(bv)?;
        let hits = py
            .detach(|| {
                self.multi()
                    .get_tversky_neighbors(&bv, ca, cb, threshold, numThreads)
            })
            .py()?;
        PyTuple::new(py, hits)
    }

    /// returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)
    #[pyo3(name = "GetContainingNeighbors", signature = (bv, numThreads = 0))]
    #[allow(non_snake_case)]
    fn get_containing_neighbors<'py>(
        &self,
        py: Python<'py>,
        bv: &Bound<'py, PyAny>,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let bv = bytes_arg(bv)?;
        let hits = py
            .detach(|| self.multi().get_containing_neighbors(&bv, numThreads))
            .py()?;
        PyTuple::new(py, hits)
    }
}

pub fn add_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFpbReader>()?;
    m.add_class::<PyMultiFpbReader>()?;
    Ok(())
}
