//! `rdoxide.DataStructs.cDataStructs`: bit vectors and similarity functions
//! (upstream `DataStructs/Wrap/{DataStructs,wrap_ExplicitBV,wrap_SparseBV,
//! wrap_BitOps,wrap_Utils}.cpp`).

use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyList, PyTuple};
use rdoxide::datastructs::bit_ops::{
    self, BitOpsVect, Metric, all_probe_bits_match, bit_vect_to_binary_text, bit_vect_to_fps_text,
    bit_vect_to_text, fold_fingerprint, similarity, similarity_wrapper, tversky,
};
use rdoxide::datastructs::bit_vect::{base64_encode, bit_vect_from_binary};
use rdoxide::datastructs::{
    BitVect, ExplicitBitVect, SparseBitVect, convert_to_explicit, from_bit_string,
    from_daylight_string,
};

use crate::datastructs_siv as siv;
use crate::errors::IntoPyResult;
use pyo3::IntoPyObjectExt;
use rdoxide::rdgeneral::threads::par_map;

/// Accepts str or bytes like `std::string` arguments of the wrappers.
pub(crate) fn bytes_arg(o: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(s) = o.extract::<String>() {
        return Ok(s.into_bytes());
    }
    o.extract::<Vec<u8>>()
}

/// `get_VectItem` / `set_VectItem` index normalization
fn norm_index(which: i64, n: u32) -> PyResult<u32> {
    if which < 0 {
        if which + (n as i64) < 0 {
            return Err(PyIndexError::new_err(format!("Index Error: {which}")));
        }
        return Ok((which + n as i64) as u32);
    }
    u32::try_from(which).map_err(|_| PyIndexError::new_err(format!("Index Error: {which}")))
}

macro_rules! bitvect_class {
    ($py:ident, $rs:ty, $name:literal, $doc:literal) => {
        #[doc = $doc]
        #[pyclass(name = $name, module = "rdoxide.DataStructs.cDataStructs", dict, skip_from_py_object)]

        #[derive(Debug)]
        pub struct $py {
            pub bv: $rs,
        }

        impl $py {
            pub(crate) fn wrap(bv: $rs) -> Self {
                $py { bv }
            }
        }

        #[pymethods]
        impl $py {
            /// Turns on a particular bit.  Returns the original state of the bit.
            #[pyo3(name = "SetBit", signature = (which))]
            fn set_bit(&mut self, which: u32) -> PyResult<bool> {
                self.bv.set_bit(which).py()
            }
            /// Turns on a set of bits.  The argument should be a tuple or list of bit ids.
            #[pyo3(name = "SetBitsFromList", signature = (onBitList))]
            #[allow(non_snake_case)]
            fn set_bits_from_list(&mut self, onBitList: Vec<i64>) -> PyResult<()> {
                for b in onBitList {
                    self.bv.set_bit(b as u32).py()?;
                }
                Ok(())
            }
            /// Turns off a particular bit.  Returns the original state of the bit.
            #[pyo3(name = "UnSetBit", signature = (which))]
            fn unset_bit(&mut self, which: u32) -> PyResult<bool> {
                self.bv.unset_bit(which).py()
            }
            /// Turns off a set of bits.  The argument should be a tuple or list of bit ids.
            #[pyo3(name = "UnSetBitsFromList", signature = (offBitList))]
            #[allow(non_snake_case)]
            fn unset_bits_from_list(&mut self, offBitList: Vec<i64>) -> PyResult<()> {
                for b in offBitList {
                    self.bv.unset_bit(b as u32).py()?;
                }
                Ok(())
            }
            /// Returns the value of a bit.
            #[pyo3(name = "GetBit", signature = (which))]
            fn get_bit(&self, which: u32) -> PyResult<bool> {
                self.bv.get_bit(which).py()
            }
            /// Returns the number of bits in the vector (the vector's size).
            #[pyo3(name = "GetNumBits")]
            fn get_num_bits(&self) -> u32 {
                self.bv.num_bits()
            }
            fn __len__(&self) -> usize {
                self.bv.num_bits() as usize
            }
            /// Returns the number of on bits.
            #[pyo3(name = "GetNumOnBits")]
            fn get_num_on_bits(&self) -> u32 {
                self.bv.num_on_bits()
            }
            /// Returns the number of off bits.
            #[pyo3(name = "GetNumOffBits")]
            fn get_num_off_bits(&self) -> u32 {
                self.bv.num_off_bits()
            }
            fn __getitem__(&self, which: i64) -> PyResult<i32> {
                let i = norm_index(which, self.bv.num_bits())?;
                Ok(i32::from(self.bv.get_bit(i).py()?))
            }
            fn __setitem__(&mut self, which: i64, val: i64) -> PyResult<()> {
                let i = norm_index(which, self.bv.num_bits())?;
                if val != 0 {
                    self.bv.set_bit(i).py()?;
                } else {
                    self.bv.unset_bit(i).py()?;
                }
                Ok(())
            }
            /// Returns a tuple containing IDs of the on bits.
            #[pyo3(name = "GetOnBits")]
            fn get_on_bits(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                crate::vects::PyVectI::new_py(py, self.bv.on_bits())
            }
            /// Returns an internal binary representation of the vector.
            #[pyo3(name = "ToBinary")]
            fn to_binary<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
                PyBytes::new(py, &self.bv.to_binary())
            }
            /// Initializes the vector from a base64 encoded binary string.
            #[pyo3(name = "FromBase64", signature = (inD))]
            #[allow(non_snake_case)]
            fn init_from_base64(&mut self, inD: &Bound<'_, PyAny>) -> PyResult<()> {
                self.bv = bit_vect_from_binary(&bytes_arg(inD)?, true).py()?;
                Ok(())
            }
            /// Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
            #[pyo3(name = "ToBase64")]
            fn to_base64(&self) -> String {
                base64_encode(&self.bv.to_binary())
            }
            /// Return the Bitvector as a python list (faster than list(vect))
            #[pyo3(name = "ToList")]
            fn to_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
                let mut v = vec![0i32; self.bv.num_bits() as usize];
                for b in self.bv.on_bits() {
                    v[b as usize] = 1;
                }
                PyList::new(py, v)
            }
            /// Returns a string of zeros and ones representing the bit vector.
            #[pyo3(name = "ToBitString")]
            fn to_bit_string(&self) -> String {
                bit_vect_to_text(&self.bv)
            }
            fn __and__(&self, o: PyRef<'_, Self>) -> Self {
                Self::wrap(self.bv.and(&o.bv))
            }
            fn __or__(&self, o: PyRef<'_, Self>) -> Self {
                Self::wrap(self.bv.or(&o.bv))
            }
            fn __xor__(&self, o: PyRef<'_, Self>) -> Self {
                Self::wrap(self.bv.xor(&o.bv))
            }
            fn __invert__(&self) -> Self {
                Self::wrap(self.bv.not())
            }
            fn __eq__(&self, o: &Bound<'_, PyAny>) -> bool {
                o.extract::<PyRef<'_, Self>>().is_ok_and(|o| o.bv == self.bv)
            }
            fn __ne__(&self, o: &Bound<'_, PyAny>) -> bool {
                !self.__eq__(o)
            }
            fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
                let py = slf.py();
                let bin = PyBytes::new(py, &slf.borrow().bv.to_binary());
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
                PyTuple::new(py, [PyBytes::new(py, &self.bv.to_binary())])
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
                slf.getattr("__dict__")?.call_method1("update", (state.get_item(0)?,))?;
                Ok(())
            }
        }
    };
}

bitvect_class!(
    PyExplicitBitVect,
    ExplicitBitVect,
    "ExplicitBitVect",
    "A class to store explicit bit vectors.\n\nThis class is most useful for situations where the size of the vector\nis relatively small (tens of thousands or smaller).\n\nFor larger vectors, use the _SparseBitVect_ class instead.\n\nAs you would expect, _ExplicitBitVects_ support a set of binary operations\nso you can do things like:\n  bv3 = bv1 & bv2  (bitwise and)\n  bv3 = bv1 | bv2  (bitwise or)\n  bv3 = bv1 ^ bv2  (bitwise xor)\n  bv3 = ~bv1       (bitwise negation)\n\nBits can be set and read using either the Set/UnsetBit() and GetBit() methods\nor by indexing (i.e. bv[i] = 1 or if bv[i]).\n"
);
bitvect_class!(
    PySparseBitVect,
    SparseBitVect,
    "SparseBitVect",
    "A class to store sparse bit vectors.\n\nThis class is most useful for situations where the size of the vector\nis large and relatively few bits are set\n\nFor smaller or denser vectors, the _ExplicitBitVect_ class is much faster.\n\nAs you would expect, _SparseBitVects_ support a set of binary operations\nso you can do things like:\n  bv3 = bv1 & bv2  (bitwise and)\n  bv3 = bv1 | bv2  (bitwise or)\n  bv3 = bv1 ^ bv2  (bitwise xor)\n  bv3 = ~bv1       (bitwise negation)\n\nBits can be set and read using either the Set/UnsetBit() and GetBit() methods\nor by indexing (i.e. bv[i] = 1 or if bv[i]).\n"
);

#[pymethods]
impl PyExplicitBitVect {
    /// `ExplicitBitVect(size)`, `ExplicitBitVect(size, bitsSet)` or `ExplicitBitVect(pkl)`
    #[new]
    #[pyo3(signature = (size, bitsSet = None))]
    #[allow(non_snake_case)]
    fn new(size: &Bound<'_, PyAny>, bitsSet: Option<bool>) -> PyResult<Self> {
        if let Ok(n) = size.extract::<u32>() {
            return Ok(Self::wrap(ExplicitBitVect::new(
                n,
                bitsSet.unwrap_or(false),
            )));
        }
        Ok(Self::wrap(
            bit_vect_from_binary(&bytes_arg(size)?, false).py()?,
        ))
    }
    fn __add__(&self, o: PyRef<'_, Self>) -> Self {
        let mut r = self.bv.clone();
        r.append(&o.bv);
        Self::wrap(r)
    }
    fn __iadd__(&mut self, o: PyRef<'_, Self>) {
        self.bv.append(&o.bv);
    }
}

#[pymethods]
impl PySparseBitVect {
    /// `SparseBitVect(size)` or `SparseBitVect(pkl)`
    #[new]
    fn new(size: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(n) = size.extract::<u32>() {
            return Ok(Self::wrap(SparseBitVect::new(n)));
        }
        Ok(Self::wrap(
            bit_vect_from_binary(&bytes_arg(size)?, false).py()?,
        ))
    }
}

/// A pair of vectors of the same kind, or a vector and a pickle.
enum Pair {
    E(ExplicitBitVect, ExplicitBitVect),
    S(SparseBitVect, SparseBitVect),
}

fn extract_pair(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>, allow_pkl: bool) -> PyResult<Pair> {
    if let Ok(x) = a.extract::<PyRef<'_, PyExplicitBitVect>>() {
        if let Ok(y) = b.extract::<PyRef<'_, PyExplicitBitVect>>() {
            return Ok(Pair::E(x.bv.clone(), y.bv.clone()));
        }
        if allow_pkl && let Ok(p) = bytes_arg(b) {
            return Ok(Pair::E(x.bv.clone(), bit_vect_from_binary(&p, false).py()?));
        }
    } else if let Ok(x) = a.extract::<PyRef<'_, PySparseBitVect>>() {
        if let Ok(y) = b.extract::<PyRef<'_, PySparseBitVect>>() {
            return Ok(Pair::S(x.bv.clone(), y.bv.clone()));
        }
        if allow_pkl && let Ok(p) = bytes_arg(b) {
            return Ok(Pair::S(x.bv.clone(), bit_vect_from_binary(&p, false).py()?));
        }
    }
    Err(PyTypeError::new_err(format!(
        "No registered converter was able to produce a C++ rvalue of type ExplicitBitVect from this Python object of type {}",
        b.get_type().name()?
    )))
}

fn wrapped<T: BitOpsVect>(a: &T, b: &T, m: Metric, dist: bool) -> PyResult<f64> {
    similarity_wrapper(a, b, |x, y| similarity(m, x, y), dist).py()
}

/// Metrics that also have a SparseIntVect overload.
fn siv_metric(name: &str) -> Option<siv::SivMetric> {
    match name {
        "DiceSimilarity" => Some(siv::SivMetric::Dice),
        "TanimotoSimilarity" => Some(siv::SivMetric::Tanimoto),
        _ => None,
    }
}

macro_rules! metric_fns {
    ($( $fname:ident, $bulk:ident, $nbr:ident, $nbr_sparse:ident, $pyname:literal, $pybulk:literal, $pynbr:literal, $pynbrs:literal, $metric:expr, $doc:literal; )*) => {
        $(
            #[doc = $doc]
            #[pyfunction]
            #[pyo3(name = $pyname, signature = (bv1, bv2, returnDistance = false, bounds = None))]
            #[allow(non_snake_case)]
            fn $fname(bv1: &Bound<'_, PyAny>, bv2: &Bound<'_, PyAny>, returnDistance: bool, bounds: Option<f64>) -> PyResult<f64> {
                if let Some(m) = siv_metric($pyname) && let Some(r) = siv::try_similarity(m, bv1, bv2, returnDistance, bounds.unwrap_or(0.0))? {
                    return Ok(r);
                }
                if bounds.is_some() {
                    return Err(PyTypeError::new_err(concat!($pyname, "() got an unexpected keyword argument 'bounds'")));
                }
                match extract_pair(bv1, bv2, true)? {
                    Pair::E(a, b) => wrapped(&a, &b, $metric, returnDistance),
                    Pair::S(a, b) => wrapped(&a, &b, $metric, returnDistance),
                }
            }

            #[doc = $doc]
            #[pyfunction]
            #[pyo3(name = $pybulk, signature = (bv1, bvList, returnDistance = false))]
            #[allow(non_snake_case)]
            fn $bulk(bv1: &Bound<'_, PyAny>, bvList: &Bound<'_, PyAny>, returnDistance: bool) -> PyResult<Vec<f64>> {
                if let Some(m) = siv_metric($pyname) && let Some(r) = siv::try_bulk(m, bv1, bvList, returnDistance)? {
                    return Ok(r);
                }
                bulk(bv1, bvList, returnDistance, |p| match p {
                    Pair::E(a, b) => similarity($metric, a, b),
                    Pair::S(a, b) => similarity($metric, a, b),
                })
            }

            #[doc = $doc]
            #[pyfunction]
            #[pyo3(name = $pynbr, signature = (bvqueries, bvList))]
            #[allow(non_snake_case)]
            fn $nbr(bvqueries: &Bound<'_, PyAny>, bvList: &Bound<'_, PyAny>) -> PyResult<Vec<(u32, f64)>> {
                neighbors::<PyExplicitBitVect, ExplicitBitVect>(bvqueries, bvList, |a, b| similarity($metric, a, b), |p| p.bv.clone())
            }

            #[doc = $doc]
            #[pyfunction]
            #[pyo3(name = $pynbrs, signature = (bvqueries, bvList))]
            #[allow(non_snake_case)]
            fn $nbr_sparse(bvqueries: &Bound<'_, PyAny>, bvList: &Bound<'_, PyAny>) -> PyResult<Vec<(u32, f64)>> {
                neighbors::<PySparseBitVect, SparseBitVect>(bvqueries, bvList, |a, b| similarity($metric, a, b), |p| p.bv.clone())
            }
        )*
        fn add_metric_fns(m: &Bound<'_, PyModule>) -> PyResult<()> {
            $(
                m.add_function(wrap_pyfunction!($fname, m)?)?;
                m.add_function(wrap_pyfunction!($bulk, m)?)?;
                m.add_function(wrap_pyfunction!($nbr, m)?)?;
                m.add_function(wrap_pyfunction!($nbr_sparse, m)?)?;
            )*
            Ok(())
        }
    };
}

/// `BulkWrapper`: no folding; each comparison must have equal lengths.
fn bulk(
    bv1: &Bound<'_, PyAny>,
    list: &Bound<'_, PyAny>,
    dist: bool,
    f: impl Fn(&Pair) -> rdoxide::Result<f64> + Sync,
) -> PyResult<Vec<f64>> {
    let mut pairs = Vec::new();
    for item in list.try_iter()? {
        pairs.push(extract_pair(bv1, &item?, false)?);
    }
    // comparisons of short lists are not worth a thread
    let target = if pairs.len() >= PAR_MIN_COMPARISONS {
        0
    } else {
        1
    };
    let scores = bv1.py().detach(|| par_map(&pairs, target, |p| f(p)));
    scores
        .into_iter()
        .map(|s| s.map(|s| if dist { 1.0 - s } else { s }).py())
        .collect()
}

/// Comparisons below which similarity loops stay on the calling thread.
const PAR_MIN_COMPARISONS: usize = 2_000_000;

/// `NeighborWrapper`
fn neighbors<P: pyo3::PyClass, T: Send + Sync>(
    queries: &Bound<'_, PyAny>,
    bvs: &Bound<'_, PyAny>,
    metric: impl Fn(&T, &T) -> rdoxide::Result<f64> + Sync,
    get: impl Fn(&P) -> T,
) -> PyResult<Vec<(u32, f64)>> {
    let pool: Vec<T> = bvs
        .try_iter()?
        .map(|o| o.and_then(|o| Ok(get(&*o.cast::<P>()?.borrow()))))
        .collect::<PyResult<_>>()?;
    let queries: Vec<T> = queries
        .try_iter()?
        .map(|q| q.and_then(|q| Ok(get(&*q.cast::<P>()?.borrow()))))
        .collect::<PyResult<_>>()?;
    let target = if queries.len() * pool.len() >= PAR_MIN_COMPARISONS {
        0
    } else {
        1
    };
    let nearest = |q: &T| -> rdoxide::Result<(u32, f64)> {
        let mut closest = -1.0;
        let mut nbr = 0u32;
        for (j, b) in pool.iter().enumerate() {
            let s = metric(q, b)?;
            if s > closest {
                closest = s;
                nbr = j as u32;
            }
        }
        Ok((nbr, closest))
    };
    let found = bvs.py().detach(|| par_map(&queries, target, nearest));
    found.into_iter().map(|r| r.py()).collect()
}

metric_fns! {
    tanimoto, bulk_tanimoto, tanimoto_nbrs, tanimoto_nbrs_s, "TanimotoSimilarity", "BulkTanimotoSimilarity", "TanimotoSimilarityNeighbors", "TanimotoSimilarityNeighbors_sparse", Metric::Tanimoto, "B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))";
    cosine, bulk_cosine, cosine_nbrs, cosine_nbrs_s, "CosineSimilarity", "BulkCosineSimilarity", "CosineSimilarityNeighbors", "CosineSimilarityNeighbors_sparse", Metric::Cosine, "B(bv1&bv2) / sqrt(B(bv1) * B(bv2))";
    kulczynski, bulk_kulczynski, kulczynski_nbrs, kulczynski_nbrs_s, "KulczynskiSimilarity", "BulkKulczynskiSimilarity", "KulczynskiSimilarityNeighbors", "KulczynskiSimilarityNeighbors_sparse", Metric::Kulczynski, "B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))";
    dice, bulk_dice, dice_nbrs, dice_nbrs_s, "DiceSimilarity", "BulkDiceSimilarity", "DiceSimilarityNeighbors", "DiceSimilarityNeighbors_sparse", Metric::Dice, "2*B(bv1&bv2) / (B(bv1) + B(bv2))";
    sokal, bulk_sokal, sokal_nbrs, sokal_nbrs_s, "SokalSimilarity", "BulkSokalSimilarity", "SokalSimilarityNeighbors", "SokalSimilarityNeighbors_sparse", Metric::Sokal, "B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))";
    mcconnaughey, bulk_mcconnaughey, mcconnaughey_nbrs, mcconnaughey_nbrs_s, "McConnaugheySimilarity", "BulkMcConnaugheySimilarity", "McConnaugheySimilarityNeighbors", "McConnaugheySimilarityNeighbors_sparse", Metric::McConnaughey, "(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))";
    asymmetric, bulk_asymmetric, asymmetric_nbrs, asymmetric_nbrs_s, "AsymmetricSimilarity", "BulkAsymmetricSimilarity", "AsymmetricSimilarityNeighbors", "AsymmetricSimilarityNeighbors_sparse", Metric::Asymmetric, "B(bv1&bv2) / min(B(bv1),B(bv2))";
    braunblanquet, bulk_braunblanquet, braunblanquet_nbrs, braunblanquet_nbrs_s, "BraunBlanquetSimilarity", "BulkBraunBlanquetSimilarity", "BraunBlanquetSimilarityNeighbors", "BraunBlanquetSimilarityNeighbors_sparse", Metric::BraunBlanquet, "B(bv1&bv2) / max(B(bv1),B(bv2))";
    russel, bulk_russel, russel_nbrs, russel_nbrs_s, "RusselSimilarity", "BulkRusselSimilarity", "RusselSimilarityNeighbors", "RusselSimilarityNeighbors_sparse", Metric::Russel, "B(bv1&bv2) / B(bv1)";
    rogotgoldberg, bulk_rogotgoldberg, rogotgoldberg_nbrs, rogotgoldberg_nbrs_s, "RogotGoldbergSimilarity", "BulkRogotGoldbergSimilarity", "RogotGoldbergSimilarityNeighbors", "RogotGoldbergSimilarityNeighbors_sparse", Metric::RogotGoldberg, "B(bv1&bv2) / B(bv1)";
}

/// B(bv1&bv2) / B(bv1|bv2)
#[pyfunction]
#[pyo3(name = "OnBitSimilarity", signature = (v1, v2))]
fn on_bit_similarity(v1: &Bound<'_, PyAny>, v2: &Bound<'_, PyAny>) -> PyResult<f64> {
    match extract_pair(v1, v2, false)? {
        Pair::E(a, b) => similarity(Metric::OnBit, &a, &b).py(),
        Pair::S(a, b) => similarity(Metric::OnBit, &a, &b).py(),
    }
}

/// (B(bv1) - B(bv1^bv2)) / B(bv1)
#[pyfunction]
#[pyo3(name = "AllBitSimilarity", signature = (v1, v2))]
fn all_bit_similarity(v1: &Bound<'_, PyAny>, v2: &Bound<'_, PyAny>) -> PyResult<f64> {
    match extract_pair(v1, v2, false)? {
        Pair::E(a, b) => similarity(Metric::AllBit, &a, &b).py(),
        Pair::S(a, b) => similarity(Metric::AllBit, &a, &b).py(),
    }
}

/// B(bv1&bv2) / B(bv1|bv2)
#[pyfunction]
#[pyo3(name = "BulkOnBitSimilarity", signature = (v1, v2, returnDistance = false))]
#[allow(non_snake_case)]
fn bulk_on_bit(
    v1: &Bound<'_, PyAny>,
    v2: &Bound<'_, PyAny>,
    returnDistance: bool,
) -> PyResult<Vec<f64>> {
    bulk_e(v1, v2, returnDistance, Metric::OnBit)
}

/// (B(bv1) - B(bv1^bv2)) / B(bv1)
#[pyfunction]
#[pyo3(name = "BulkAllBitSimilarity", signature = (v1, v2, returnDistance = false))]
#[allow(non_snake_case)]
fn bulk_all_bit(
    v1: &Bound<'_, PyAny>,
    v2: &Bound<'_, PyAny>,
    returnDistance: bool,
) -> PyResult<Vec<f64>> {
    bulk_e(v1, v2, returnDistance, Metric::AllBit)
}

/// The bulk OnBit/AllBit functions only take ExplicitBitVects upstream.
fn bulk_e(
    v1: &Bound<'_, PyAny>,
    list: &Bound<'_, PyAny>,
    dist: bool,
    m: Metric,
) -> PyResult<Vec<f64>> {
    let a = v1.extract::<PyRef<'_, PyExplicitBitVect>>()?.bv.clone();
    let mut res = Vec::new();
    for item in list.try_iter()? {
        let b = item?.extract::<PyRef<'_, PyExplicitBitVect>>()?.bv.clone();
        let mut s = similarity(m, &a, &b).py()?;
        if dist {
            s = 1.0 - s;
        }
        res.push(s);
    }
    Ok(res)
}

/// B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)
#[pyfunction]
#[pyo3(name = "TverskySimilarity", signature = (bv1, bv2, a, b, returnDistance = false, bounds = None))]
#[allow(non_snake_case)]
fn tversky_similarity(
    bv1: &Bound<'_, PyAny>,
    bv2: &Bound<'_, PyAny>,
    a: f64,
    b: f64,
    returnDistance: bool,
    bounds: Option<f64>,
) -> PyResult<f64> {
    if let Some(r) =
        siv::try_similarity(siv::SivMetric::Tversky(a, b), bv1, bv2, returnDistance, 0.0)?
    {
        return Ok(r);
    }
    if bounds.is_some() {
        return Err(PyTypeError::new_err(
            "TverskySimilarity() got an unexpected keyword argument 'bounds'",
        ));
    }
    match extract_pair(bv1, bv2, true)? {
        Pair::E(x, y) => {
            similarity_wrapper(&x, &y, |p, q| tversky(p, q, a, b), returnDistance).py()
        }
        Pair::S(x, y) => {
            similarity_wrapper(&x, &y, |p, q| tversky(p, q, a, b), returnDistance).py()
        }
    }
}

/// B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)
#[pyfunction]
#[pyo3(name = "BulkTverskySimilarity", signature = (bv1, bvList, a, b, returnDistance = false))]
#[allow(non_snake_case)]
fn bulk_tversky(
    bv1: &Bound<'_, PyAny>,
    bvList: &Bound<'_, PyAny>,
    a: f64,
    b: f64,
    returnDistance: bool,
) -> PyResult<Vec<f64>> {
    if let Some(r) = siv::try_bulk(siv::SivMetric::Tversky(a, b), bv1, bvList, returnDistance)? {
        return Ok(r);
    }
    bulk(bv1, bvList, returnDistance, |p| match p {
        Pair::E(x, y) => tversky(x, y, a, b),
        Pair::S(x, y) => tversky(x, y, a, b),
    })
}

/// Returns a 2-tuple: (B(bv1&bv2) / B(bv1), B(bv1&bv2) / B(bv2))
#[pyfunction]
#[pyo3(name = "OnBitProjSimilarity", signature = (bv1, bv2))]
fn on_bit_proj(
    py: Python<'_>,
    bv1: &Bound<'_, PyAny>,
    bv2: &Bound<'_, PyAny>,
) -> PyResult<Py<PyAny>> {
    let r = match extract_pair(bv1, bv2, false)? {
        Pair::E(a, b) => bit_ops::on_bit_proj_similarity(&a, &b),
        Pair::S(a, b) => bit_ops::on_bit_proj_similarity(&a, &b),
    }
    .py()?;
    crate::vects::PyVectD::new_py(py, r.to_vec())
}

#[pyfunction]
#[pyo3(name = "OffBitProjSimilarity", signature = (bv1, bv2))]
fn off_bit_proj(
    py: Python<'_>,
    bv1: &Bound<'_, PyAny>,
    bv2: &Bound<'_, PyAny>,
) -> PyResult<Py<PyAny>> {
    let r = match extract_pair(bv1, bv2, false)? {
        Pair::E(a, b) => bit_ops::off_bit_proj_similarity(&a, &b),
        Pair::S(a, b) => bit_ops::off_bit_proj_similarity(&a, &b),
    }
    .py()?;
    crate::vects::PyVectD::new_py(py, r.to_vec())
}

/// Returns the total number of bits in common between the two bit vectors
#[pyfunction]
#[pyo3(name = "NumBitsInCommon", signature = (bv1, bv2))]
fn num_bits_in_common(bv1: &Bound<'_, PyAny>, bv2: &Bound<'_, PyAny>) -> PyResult<i32> {
    match extract_pair(bv1, bv2, false)? {
        Pair::E(a, b) => bit_ops::num_bits_in_common(&a, &b).py(),
        Pair::S(a, b) => bit_ops::num_bits_in_common(&a, &b).py(),
    }
}

/// Returns the number of on bits in common between the two bit vectors
#[pyfunction]
#[pyo3(name = "OnBitsInCommon", signature = (bv1, bv2))]
fn on_bits_in_common<'py>(
    py: Python<'py>,
    bv1: &Bound<'py, PyAny>,
    bv2: &Bound<'py, PyAny>,
) -> PyResult<Py<PyAny>> {
    let v = match extract_pair(bv1, bv2, false)? {
        Pair::E(a, b) => bit_ops::on_bits_in_common(&a, &b),
        Pair::S(a, b) => bit_ops::on_bits_in_common(&a, &b),
    }
    .py()?;
    crate::vects::PyVectI::new_py(py, v)
}

/// Returns the number of off bits in common between the two bit vectors
#[pyfunction]
#[pyo3(name = "OffBitsInCommon", signature = (bv1, bv2))]
fn off_bits_in_common<'py>(
    py: Python<'py>,
    bv1: &Bound<'py, PyAny>,
    bv2: &Bound<'py, PyAny>,
) -> PyResult<Py<PyAny>> {
    let v = match extract_pair(bv1, bv2, false)? {
        Pair::E(a, b) => bit_ops::off_bits_in_common(&a, &b),
        Pair::S(a, b) => bit_ops::off_bits_in_common(&a, &b),
    }
    .py()?;
    crate::vects::PyVectI::new_py(py, v)
}

/// Folds the fingerprint by the provided amount. The default, foldFactor=2, returns a fingerprint that is half the size of the original.
#[pyfunction]
#[pyo3(name = "FoldFingerprint", signature = (bv, foldFactor = 2))]
#[allow(non_snake_case)]
fn fold_fp(py: Python<'_>, bv: &Bound<'_, PyAny>, foldFactor: u32) -> PyResult<Py<PyAny>> {
    if let Ok(e) = bv.extract::<PyRef<'_, PyExplicitBitVect>>() {
        let r = PyExplicitBitVect::wrap(fold_fingerprint(&e.bv, foldFactor).py()?);
        return Ok(Py::new(py, r)?.into_any());
    }
    let s = bv.extract::<PyRef<'_, PySparseBitVect>>()?;
    let r = PySparseBitVect::wrap(fold_fingerprint(&s.bv, foldFactor).py()?);
    Ok(Py::new(py, r)?.into_any())
}

/// Returns True if all bits in the first argument match all bits in the
///   vector defined by the pickle in the second argument.
#[pyfunction]
#[pyo3(name = "AllProbeBitsMatch", signature = (probe, r#ref))]
fn all_probe_bits(probe: &Bound<'_, PyAny>, r#ref: &Bound<'_, PyAny>) -> PyResult<bool> {
    match extract_pair(probe, r#ref, true)? {
        Pair::E(a, b) => Ok(a.is_subset_of(&b)),
        Pair::S(a, b) => all_probe_bits_match(&a, &b).py(),
    }
}

fn with_bv<R>(
    bv: &Bound<'_, PyAny>,
    fe: impl Fn(&ExplicitBitVect) -> R,
    fs: impl Fn(&SparseBitVect) -> R,
) -> PyResult<R> {
    if let Ok(e) = bv.extract::<PyRef<'_, PyExplicitBitVect>>() {
        return Ok(fe(&e.bv));
    }
    Ok(fs(&bv.extract::<PyRef<'_, PySparseBitVect>>()?.bv))
}

/// Returns a string of zeros and ones representing the bit vector.
#[pyfunction]
#[pyo3(name = "BitVectToText", signature = (bv1))]
fn bv_to_text(bv1: &Bound<'_, PyAny>) -> PyResult<String> {
    with_bv(bv1, bit_vect_to_text, bit_vect_to_text)
}

/// Returns an FPS string representing the bit vector.
#[pyfunction]
#[pyo3(name = "BitVectToFPSText", signature = (bv1))]
fn bv_to_fps(bv1: &Bound<'_, PyAny>) -> PyResult<String> {
    with_bv(bv1, bit_vect_to_fps_text, bit_vect_to_fps_text)
}

/// Returns a binary string (byte array) representing the bit vector.
#[pyfunction]
#[pyo3(name = "BitVectToBinaryText", signature = (bv))]
fn bv_to_binary<'py>(py: Python<'py>, bv: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyBytes>> {
    let b = with_bv(bv, bit_vect_to_binary_text, bit_vect_to_binary_text)?;
    Ok(PyBytes::new(py, &b))
}

/// Converts a SparseBitVector to an ExplicitBitVector and returns the ExplicitBitVector
#[pyfunction]
#[pyo3(name = "ConvertToExplicit", signature = (sbv))]
fn convert_to_explicit_py(sbv: PyRef<'_, PySparseBitVect>) -> PyExplicitBitVect {
    PyExplicitBitVect::wrap(convert_to_explicit(&sbv.bv))
}

/// Fill a BitVect using an ASCII (Daylight) encoding of a fingerprint.
///
///    **Arguments**
///      - bv: either a _SparseBitVect_ or an _ExplicitBitVect_
///      - txt: a string with the Daylight encoding (this is the text that
///             the Daylight tools put in the FP field of a TDT)
#[pyfunction]
#[pyo3(name = "InitFromDaylightString", signature = (sbv, s))]
fn init_from_daylight_string(sbv: &Bound<'_, PyAny>, s: &str) -> PyResult<()> {
    if let Ok(mut e) = sbv.extract::<PyRefMut<'_, PyExplicitBitVect>>() {
        return from_daylight_string(&mut e.bv, s).py();
    }
    let mut sp = sbv.extract::<PyRefMut<'_, PySparseBitVect>>()?;
    from_daylight_string(&mut sp.bv, s).py()
}

/// Creates an ExplicitBitVect from a bit string (string of 0s and 1s).
#[pyfunction]
#[pyo3(name = "CreateFromBitString", signature = (bits))]
fn create_from_bit_string(bits: &str) -> PyResult<PyExplicitBitVect> {
    let mut bv = ExplicitBitVect::new(bits.len() as u32, false);
    from_bit_string(&mut bv, bits.as_bytes()).py()?;
    Ok(PyExplicitBitVect::wrap(bv))
}

/// Creates an ExplicitBitVect from an FPS string.
#[pyfunction]
#[pyo3(name = "CreateFromFPSText", signature = (fps))]
fn create_from_fps(fps: &str) -> PyResult<PyExplicitBitVect> {
    if fps.len() % 2 != 0 {
        return Err(PyValueError::new_err(
            "input string must have an even number of characters",
        ));
    }
    let mut bv = ExplicitBitVect::new(fps.len() as u32 * 4, false);
    bit_ops::update_bit_vect_from_fps_text(&mut bv, fps.as_bytes()).py()?;
    Ok(PyExplicitBitVect::wrap(bv))
}

/// Creates an ExplicitBitVect from a binary string (byte array).
#[pyfunction]
#[pyo3(name = "CreateFromBinaryText", signature = (fps))]
fn create_from_binary(fps: &Bound<'_, PyAny>) -> PyResult<PyExplicitBitVect> {
    let data = bytes_arg(fps)?;
    let mut bv = ExplicitBitVect::new(data.len() as u32 * 8, false);
    bit_ops::update_bit_vect_from_binary_text(&mut bv, &data).py()?;
    Ok(PyExplicitBitVect::wrap(bv))
}

/// Fill the numpy array `destArray` with the bits of `bv`.
#[pyfunction]
#[pyo3(name = "ConvertToNumpyArray", signature = (bv, destArray))]
#[allow(non_snake_case)]
fn convert_to_numpy(bv: &Bound<'_, PyAny>, destArray: &Bound<'_, PyAny>) -> PyResult<()> {
    let py = bv.py();
    let vals: Vec<Py<PyAny>> = if let Some(v) = crate::datastructs_vv::numpy_values(bv) {
        v
    } else if let Some(v) = siv::get_siv(bv) {
        let mut vals = vec![0i64; v.length() as usize];
        for (&k, &x) in v.nonzero_elements() {
            vals[k as usize] = x as i64;
        }
        vals.into_iter()
            .map(|x| x.into_py_any(py))
            .collect::<PyResult<_>>()?
    } else {
        let bv = bv.extract::<PyRef<'_, PyExplicitBitVect>>()?;
        let mut vals = vec![0i64; bv.bv.num_bits() as usize];
        for b in bv.bv.on_bits() {
            vals[b as usize] = 1;
        }
        vals.into_iter()
            .map(|x| x.into_py_any(py))
            .collect::<PyResult<_>>()?
    };
    let n = vals.len();
    destArray.call_method1("resize", ((n,),)).or_else(|_| {
        destArray.call_method(
            "resize",
            ((n,),),
            Some(&[("refcheck", false)].into_py_dict(destArray.py())?),
        )
    })?;
    for (i, v) in vals.into_iter().enumerate() {
        destArray.set_item(i, v)?;
    }
    Ok(())
}

use pyo3::types::IntoPyDict;

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing an assortment of functionality for basic data structures.\n\nAt the moment the data structures defined are:\n  Bit Vector classes (for storing signatures, fingerprints and the like:\n    - ExplicitBitVect: class for relatively small (10s of thousands of bits) or\n                       dense bit vectors.\n    - SparseBitVect:   class for large, sparse bit vectors\n  DiscreteValueVect:   class for storing vectors of integers\n  SparseIntVect:       class for storing sparse vectors of integers\n",
    )?;
    m.add_class::<PyExplicitBitVect>()?;
    m.add_class::<PySparseBitVect>()?;
    siv::add_classes(m)?;
    crate::datastructs_fpb::add_classes(m)?;
    crate::datastructs_vv::add_classes(m)?;
    add_metric_fns(m)?;
    for f in [
        wrap_pyfunction!(on_bit_similarity, m)?,
        wrap_pyfunction!(all_bit_similarity, m)?,
        wrap_pyfunction!(bulk_on_bit, m)?,
        wrap_pyfunction!(bulk_all_bit, m)?,
        wrap_pyfunction!(tversky_similarity, m)?,
        wrap_pyfunction!(bulk_tversky, m)?,
        wrap_pyfunction!(on_bit_proj, m)?,
        wrap_pyfunction!(off_bit_proj, m)?,
        wrap_pyfunction!(num_bits_in_common, m)?,
        wrap_pyfunction!(on_bits_in_common, m)?,
        wrap_pyfunction!(off_bits_in_common, m)?,
        wrap_pyfunction!(fold_fp, m)?,
        wrap_pyfunction!(all_probe_bits, m)?,
        wrap_pyfunction!(bv_to_text, m)?,
        wrap_pyfunction!(bv_to_fps, m)?,
        wrap_pyfunction!(bv_to_binary, m)?,
        wrap_pyfunction!(convert_to_explicit_py, m)?,
        wrap_pyfunction!(init_from_daylight_string, m)?,
        wrap_pyfunction!(create_from_bit_string, m)?,
        wrap_pyfunction!(create_from_fps, m)?,
        wrap_pyfunction!(create_from_binary, m)?,
        wrap_pyfunction!(convert_to_numpy, m)?,
    ] {
        m.add_function(f)?;
    }
    Ok(())
}
