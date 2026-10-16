//! `rdoxide.SimDivFilters.rdSimDivPickers` (upstream `SimDivPickers/Wrap`):
//! `MaxMinPicker` and `LeaderPicker`.

use std::sync::Mutex;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::RDKitError;
use rdoxide::datastructs::ExplicitBitVect;
use rdoxide::datastructs::bit_ops::{Metric, similarity};
use rdoxide::ml::cluster::murtagh::Linkage;
use rdoxide::sim_div_pickers::{hierarchical, leader, max_min};

use crate::datastructs::PyExplicitBitVect;
use crate::enums::add_enum_exported;
use crate::errors::IntoPyResult;
use crate::vects::{PyVectI, PyVectVectI};

/// A Python function `f(i, j) -> float` as a pool distance. The first
/// failure of the function is kept so that it can be re-raised after the
/// picker stopped.
struct PyDistance {
    func: Py<PyAny>,
    error: Mutex<Option<PyErr>>,
}

impl PyDistance {
    fn new(func: &Bound<'_, PyAny>) -> Self {
        PyDistance {
            func: func.clone().unbind(),
            error: Mutex::new(None),
        }
    }

    fn dist(&self, i: u32, j: u32) -> rdoxide::Result<f64> {
        Python::attach(|py| self.func.call1(py, (i, j))?.extract::<f64>(py)).map_err(|e| {
            self.error
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .get_or_insert(e);
            RDKitError::Runtime("the distance function failed".into())
        })
    }

    /// The outcome of a picker run that used [`PyDistance::dist`].
    fn finish<T>(self, result: rdoxide::Result<T>) -> PyResult<T> {
        match self.error.into_inner().unwrap_or_else(|e| e.into_inner()) {
            Some(err) => Err(err),
            None => result.py(),
        }
    }
}

/// The `objects` of the bit vector pickers: the first `pool_size` items.
fn bit_vectors<'py>(
    objects: &Bound<'py, PyAny>,
    pool_size: u32,
) -> PyResult<Vec<PyRef<'py, PyExplicitBitVect>>> {
    (0..pool_size)
        .map(|i| -> PyResult<PyRef<'py, PyExplicitBitVect>> {
            Ok(objects
                .get_item(i)?
                .extract::<PyRef<'py, PyExplicitBitVect>>()?)
        })
        .collect()
}

/// `1 - TanimotoSimilarity` between two pool bit vectors.
fn tanimoto_distance(bvs: &[&ExplicitBitVect], i: u32, j: u32) -> rdoxide::Result<f64> {
    Ok(1.0 - similarity(Metric::Tanimoto, bvs[i as usize], bvs[j as usize])?)
}

fn first_picks(first: &Bound<'_, PyAny>) -> PyResult<Vec<i32>> {
    first.try_iter()?.map(|x| x?.extract::<i32>()).collect()
}

fn warn_use_cache(use_cache: &Option<Bound<'_, PyAny>>) {
    if use_cache.as_ref().is_some_and(|c| !c.is_none()) {
        rdoxide::rd_log_warning!("the useCache argument is deprecated and ignored\n");
    }
}

fn picks_vect(py: Python<'_>, picks: Vec<i32>) -> PyResult<Py<PyAny>> {
    PyVectI::new_py(py, picks)
}

fn picks_with_threshold(py: Python<'_>, res: max_min::MaxMinPicks) -> PyResult<Bound<'_, PyTuple>> {
    PyTuple::new(
        py,
        [
            picks_vect(py, res.picks)?,
            res.threshold.into_pyobject(py)?.into_any().unbind(),
        ],
    )
}

/// A class for diversity picking of items using the MaxMin Algorithm
#[pyclass(
    name = "MaxMinPicker",
    module = "rdoxide.SimDivFilters.rdSimDivPickers",
    dict
)]
#[derive(Debug, Default)]
pub struct PyMaxMinPicker;

#[pymethods]
impl PyMaxMinPicker {
    #[new]
    fn new() -> Self {
        PyMaxMinPicker
    }

    /// Pick a subset of items from a pool of items using the MaxMin Algorithm
    /// Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
    ///
    /// ARGUMENTS:
    ///   - distMat: 1D distance matrix (only the lower triangle elements)
    ///   - poolSize: number of items in the pool
    ///   - pickSize: number of items to pick from the pool
    ///   - firstPicks: (optional) the first items to be picked (seeds the list)
    ///   - seed: (optional) seed for the random number generator
    #[pyo3(name = "Pick", signature = (distMat, poolSize, pickSize, firstPicks = None, seed = -1))]
    #[allow(non_snake_case)]
    fn pick(
        &self,
        py: Python<'_>,
        distMat: &Bound<'_, PyAny>,
        poolSize: u32,
        pickSize: u32,
        firstPicks: Option<Bound<'_, PyAny>>,
        seed: i32,
    ) -> PyResult<Py<PyAny>> {
        if pickSize >= poolSize {
            return Err(PyValueError::new_err("pickSize must be less than poolSize"));
        }
        let np = py.import("numpy")?;
        if !distMat.is_instance(&np.getattr("ndarray")?)? {
            return Err(PyValueError::new_err(
                "distance mat argument must be a numpy matrix",
            ));
        }
        let flat = np
            .call_method1("ascontiguousarray", (distMat, "float64"))?
            .call_method0("ravel")?;
        if distMat.getattr("ndim")?.extract::<usize>()? != 1 {
            return Err(PyValueError::new_err(
                "distance mat argument must be one-dimensional",
            ));
        }
        let mat: Vec<f64> = flat.call_method0("tolist")?.extract()?;
        let first = firstPicks.as_ref().map_or(Ok(Vec::new()), first_picks)?;
        let picks = py
            .detach(|| max_min::pick(&mat, poolSize, pickSize, &first, seed))
            .py()?;
        picks_vect(py, picks)
    }

    /// Pick a subset of items from a pool of items using the MaxMin Algorithm
    /// Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
    /// ARGUMENTS:
    ///
    ///   - distFunc: a function that should take two indices and return the
    ///               distance between those two points.
    ///               NOTE: the implementation caches distance values, so the
    ///               client code does not need to do so; indeed, it should not.
    ///   - poolSize: number of items in the pool
    ///   - pickSize: number of items to pick from the pool
    ///   - firstPicks: (optional) the first items to be picked (seeds the list)
    ///   - seed: (optional) seed for the random number generator
    ///   - useCache: IGNORED
    #[pyo3(name = "LazyPick", signature = (distFunc, poolSize, pickSize, firstPicks = None, seed = -1, useCache = None))]
    #[allow(non_snake_case)]
    fn lazy_pick(
        &self,
        py: Python<'_>,
        distFunc: &Bound<'_, PyAny>,
        poolSize: u32,
        pickSize: u32,
        firstPicks: Option<Bound<'_, PyAny>>,
        seed: i32,
        useCache: Option<Bound<'_, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        warn_use_cache(&useCache);
        let first = firstPicks.as_ref().map_or(Ok(Vec::new()), first_picks)?;
        let dist = PyDistance::new(distFunc);
        let res = max_min::lazy_pick(
            |i, j| dist.dist(i, j),
            poolSize,
            pickSize,
            &first,
            seed,
            -1.0,
        );
        picks_vect(py, dist.finish(res)?.picks)
    }

    /// Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
    /// Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
    /// ARGUMENTS:
    ///
    ///   - vectors: a sequence of the bit vectors that should be picked from.
    ///   - poolSize: number of items in the pool
    ///   - pickSize: number of items to pick from the pool
    ///   - firstPicks: (optional) the first items to be picked (seeds the list)
    ///   - seed: (optional) seed for the random number generator
    ///   - useCache: IGNORED.
    #[pyo3(name = "LazyBitVectorPick", signature = (objects, poolSize, pickSize, firstPicks = None, seed = -1, useCache = None))]
    #[allow(non_snake_case)]
    fn lazy_bit_vector_pick(
        &self,
        py: Python<'_>,
        objects: &Bound<'_, PyAny>,
        poolSize: u32,
        pickSize: u32,
        firstPicks: Option<Bound<'_, PyAny>>,
        seed: i32,
        useCache: Option<Bound<'_, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        warn_use_cache(&useCache);
        let first = firstPicks.as_ref().map_or(Ok(Vec::new()), first_picks)?;
        let guards = bit_vectors(objects, poolSize)?;
        let bvs: Vec<&ExplicitBitVect> = guards.iter().map(|g| &g.bv).collect();
        let res = py
            .detach(|| {
                max_min::lazy_pick(
                    |i, j| tanimoto_distance(&bvs, i, j),
                    poolSize,
                    pickSize,
                    &first,
                    seed,
                    -1.0,
                )
            })
            .py()?;
        picks_vect(py, res.picks)
    }

    /// Pick a subset of items from a pool of items using the MaxMin Algorithm
    /// Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
    /// ARGUMENTS:
    ///
    ///   - distFunc: a function that should take two indices and return the
    ///               distance between those two points.
    ///               NOTE: the implementation caches distance values, so the
    ///               client code does not need to do so; indeed, it should not.
    ///   - poolSize: number of items in the pool
    ///   - pickSize: number of items to pick from the pool
    ///   - threshold: stop picking when the distance goes below this value
    ///   - firstPicks: (optional) the first items to be picked (seeds the list)
    ///   - seed: (optional) seed for the random number generator
    #[pyo3(name = "LazyPickWithThreshold", signature = (distFunc, poolSize, pickSize, threshold, firstPicks = None, seed = -1))]
    #[allow(non_snake_case)]
    fn lazy_pick_with_threshold<'py>(
        &self,
        py: Python<'py>,
        distFunc: &Bound<'py, PyAny>,
        poolSize: u32,
        pickSize: u32,
        threshold: f64,
        firstPicks: Option<Bound<'py, PyAny>>,
        seed: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let first = firstPicks.as_ref().map_or(Ok(Vec::new()), first_picks)?;
        let dist = PyDistance::new(distFunc);
        let res = max_min::lazy_pick(
            |i, j| dist.dist(i, j),
            poolSize,
            pickSize,
            &first,
            seed,
            threshold,
        );
        picks_with_threshold(py, dist.finish(res)?)
    }

    /// Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
    /// Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
    /// ARGUMENTS:
    ///
    ///   - vectors: a sequence of the bit vectors that should be picked from.
    ///   - poolSize: number of items in the pool
    ///   - pickSize: number of items to pick from the pool
    ///   - threshold: stop picking when the distance goes below this value
    ///   - firstPicks: (optional) the first items to be picked (seeds the list)
    ///   - seed: (optional) seed for the random number generator
    #[pyo3(name = "LazyBitVectorPickWithThreshold", signature = (objects, poolSize, pickSize, threshold, firstPicks = None, seed = -1))]
    #[allow(non_snake_case)]
    fn lazy_bit_vector_pick_with_threshold<'py>(
        &self,
        py: Python<'py>,
        objects: &Bound<'py, PyAny>,
        poolSize: u32,
        pickSize: u32,
        threshold: f64,
        firstPicks: Option<Bound<'py, PyAny>>,
        seed: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let first = firstPicks.as_ref().map_or(Ok(Vec::new()), first_picks)?;
        let guards = bit_vectors(objects, poolSize)?;
        let bvs: Vec<&ExplicitBitVect> = guards.iter().map(|g| &g.bv).collect();
        let res = py
            .detach(|| {
                max_min::lazy_pick(
                    |i, j| tanimoto_distance(&bvs, i, j),
                    poolSize,
                    pickSize,
                    &first,
                    seed,
                    threshold,
                )
            })
            .py()?;
        picks_with_threshold(py, res)
    }
}

/// A class for diversity picking of items using Roger Sayle's Leader algorithm
/// (analogous to sphere exclusion). The algorithm is currently unpublished, but
/// a description is available in this presentation from the 2019 RDKit UGM:
/// https://github.com/rdkit/UGM_2019/raw/master/Presentations/Sayle_Clustering.pdf
#[pyclass(
    name = "LeaderPicker",
    module = "rdoxide.SimDivFilters.rdSimDivPickers",
    dict
)]
#[derive(Debug, Default)]
pub struct PyLeaderPicker;

#[pymethods]
impl PyLeaderPicker {
    #[new]
    fn new() -> Self {
        PyLeaderPicker
    }

    /// Pick a subset of items from a collection of bit vectors using Tanimoto
    /// distance. The threshold value is a *distance* (i.e. 1-similarity). The
    /// distances of a round are computed on up to numThreads threads.
    #[pyo3(name = "LazyBitVectorPick", signature = (objects, poolSize, threshold, pickSize = 0, firstPicks = None, numThreads = 1))]
    #[allow(non_snake_case)]
    fn lazy_bit_vector_pick(
        &self,
        py: Python<'_>,
        objects: &Bound<'_, PyAny>,
        poolSize: u32,
        threshold: f64,
        pickSize: u32,
        firstPicks: Option<Bound<'_, PyAny>>,
        numThreads: i32,
    ) -> PyResult<Py<PyAny>> {
        let first = firstPicks.as_ref().map_or(Ok(Vec::new()), first_picks)?;
        let guards = bit_vectors(objects, poolSize)?;
        let bvs: Vec<&ExplicitBitVect> = guards.iter().map(|g| &g.bv).collect();
        let picks = py
            .detach(|| {
                leader::lazy_pick(
                    |i, j| tanimoto_distance(&bvs, i, j),
                    poolSize,
                    pickSize,
                    &first,
                    threshold,
                    numThreads,
                )
            })
            .py()?;
        picks_vect(py, picks)
    }

    /// Pick a subset of items from a pool of items using the user-provided
    /// function to determine distances. The numThreads argument is ignored: the
    /// function is called from Python, one call at a time.
    #[pyo3(name = "LazyPick", signature = (distFunc, poolSize, threshold, pickSize = 0, firstPicks = None, numThreads = 1))]
    #[allow(non_snake_case, unused_variables)]
    fn lazy_pick(
        &self,
        py: Python<'_>,
        distFunc: &Bound<'_, PyAny>,
        poolSize: u32,
        threshold: f64,
        pickSize: u32,
        firstPicks: Option<Bound<'_, PyAny>>,
        numThreads: i32,
    ) -> PyResult<Py<PyAny>> {
        let first = firstPicks.as_ref().map_or(Ok(Vec::new()), first_picks)?;
        let dist = PyDistance::new(distFunc);
        let res = leader::lazy_pick(
            |i, j| dist.dist(i, j),
            poolSize,
            pickSize,
            &first,
            threshold,
            1,
        );
        picks_vect(py, dist.finish(res)?)
    }
}

/// The 1-D distance matrix of the hierarchical picker as a private copy
/// (the clustering overwrites its input).
fn copied_distance_matrix(py: Python<'_>, dist_mat: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    let np = py.import("numpy")?;
    if !dist_mat.is_instance(&np.getattr("ndarray")?)? {
        return Err(PyValueError::new_err(
            "distance mat argument must be a numpy matrix",
        ));
    }
    if dist_mat.getattr("ndim")?.extract::<usize>()? != 1 {
        return Err(PyValueError::new_err(
            "distance mat argument must be one-dimensional",
        ));
    }
    np.call_method1("array", (dist_mat, "float64"))?
        .call_method0("tolist")?
        .extract()
}

/// A class for diversity picking of items using Hierarchical Clustering
#[pyclass(
    name = "HierarchicalClusterPicker",
    module = "rdoxide.SimDivFilters.rdSimDivPickers",
    dict
)]
#[derive(Debug)]
pub struct PyHierarchicalClusterPicker {
    method: Linkage,
}

#[pymethods]
impl PyHierarchicalClusterPicker {
    #[new]
    #[pyo3(signature = (clusterMethod))]
    #[allow(non_snake_case)]
    fn new(clusterMethod: i64) -> PyResult<Self> {
        let method = Linkage::from_option(clusterMethod)
            .ok_or_else(|| PyValueError::new_err("clusterMethod must be a ClusterMethod value"))?;
        Ok(PyHierarchicalClusterPicker { method })
    }

    /// Pick a diverse subset of items from a pool of items using hierarchical
    /// clustering
    ///
    /// ARGUMENTS:
    ///   - distMat: 1D distance matrix (only the lower triangle elements)
    ///   - poolSize: number of items in the pool
    ///   - pickSize: number of items to pick from the pool
    #[pyo3(name = "Pick", signature = (distMat, poolSize, pickSize))]
    #[allow(non_snake_case)]
    fn pick(
        &self,
        py: Python<'_>,
        distMat: &Bound<'_, PyAny>,
        poolSize: u32,
        pickSize: u32,
    ) -> PyResult<Py<PyAny>> {
        if pickSize >= poolSize {
            return Err(PyValueError::new_err("pickSize must be less than poolSize"));
        }
        let mat = copied_distance_matrix(py, distMat)?;
        let picks = py
            .detach(|| hierarchical::pick(&mat, poolSize, pickSize, self.method))
            .py()?;
        picks_vect(py, picks)
    }

    /// Return a list of clusters of item from the pool using hierarchical
    /// clustering
    ///
    /// ARGUMENTS:
    ///   - distMat: 1D distance matrix (only the lower triangle elements)
    ///   - poolSize: number of items in the pool
    ///   - pickSize: number of items to pick from the pool
    #[pyo3(name = "Cluster", signature = (distMat, poolSize, pickSize))]
    #[allow(non_snake_case)]
    fn cluster(
        &self,
        py: Python<'_>,
        distMat: &Bound<'_, PyAny>,
        poolSize: u32,
        pickSize: u32,
    ) -> PyResult<Py<PyAny>> {
        let mat = copied_distance_matrix(py, distMat)?;
        let clusters = py
            .detach(|| hierarchical::cluster(&mat, poolSize, pickSize, self.method))
            .py()?;
        PyVectVectI::new_py(py, clusters)
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing the diversity and similarity pickers",
    )?;
    m.add_class::<PyMaxMinPicker>()?;
    m.add_class::<PyLeaderPicker>()?;
    m.add_class::<PyHierarchicalClusterPicker>()?;
    add_enum_exported(
        m,
        "ClusterMethod",
        &[
            ("WARD", Linkage::Ward as i64),
            ("SLINK", Linkage::Slink as i64),
            ("CLINK", Linkage::Clink as i64),
            ("UPGMA", Linkage::Upgma as i64),
            ("MCQUITTY", Linkage::McQuitty as i64),
            ("GOWER", Linkage::Gower as i64),
            ("CENTROID", Linkage::Centroid as i64),
        ],
        None,
    )?;
    Ok(())
}
