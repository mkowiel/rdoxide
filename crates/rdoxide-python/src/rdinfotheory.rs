//! `rdoxide.ML.InfoTheory.rdInfoTheory` (upstream `ML/InfoTheory/Wrap/rdInfoTheory.cpp`):
//! `InfoEntropy`, `InfoGain` and `ChiSquare` of numpy count arrays, the
//! `InfoBitRanker` with its `InfoType`, and the `BitCorrMatGenerator`
//! (`Wrap/InfoBitRanker.cpp`, `Wrap/BitCorrMatGenerator.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::datastructs::bit_vect::BitVect;
use rdoxide::ml::info_bit_ranker::{BitCorrMatGenerator, InfoBitRanker, InfoType};
use rdoxide::ml::info_theory::{chi_square, info_entropy, info_entropy_gain};

use crate::datastructs::{PyExplicitBitVect, PySparseBitVect};
use crate::enums::add_enum_exported;
use crate::errors::{IntoPyResult, to_py};
use crate::pyconv::float64_array;

/// The data of a numpy array in one of the types upstream instantiates its
/// templates for (anything else is converted to doubles, as upstream does).
enum Counts {
    F32(Vec<f32>),
    F64(Vec<f64>),
    I32(Vec<i32>),
    I64(Vec<i64>),
}

/// Runs `$body` with `$data` bound to the slice of whichever type `$counts` holds.
macro_rules! with_counts {
    ($counts:expr, |$data:ident| $body:expr) => {
        match &$counts {
            Counts::F32($data) => $body,
            Counts::F64($data) => $body,
            Counts::I32($data) => $body,
            Counts::I64($data) => $body,
        }
    };
}

fn conversion_error() -> PyErr {
    PyValueError::new_err("could not convert argument to a numeric array")
}

/// `nativeOrDoubleType` + `PyArray_ContiguousFromObject`: the contiguous data
/// of `obj` with `min_dim..=max_dim` dimensions, and its shape.
fn load_counts(
    py: Python<'_>,
    obj: &Bound<'_, PyAny>,
    min_dim: usize,
    max_dim: usize,
) -> PyResult<(Counts, Vec<usize>)> {
    let np = py.import("numpy")?;
    if !obj.is_instance(&np.getattr("ndarray")?)? {
        return Err(PyValueError::new_err("Expecting a Numeric array object"));
    }
    let dtype = obj.getattr("dtype")?;
    // NPY_FLOAT, NPY_DOUBLE, NPY_INT, NPY_LONG and NPY_LONGLONG are read as they are
    let native = match dtype.getattr("num")?.extract::<i32>()? {
        11 => Some("f32"),
        12 => Some("f64"),
        5 | 7 | 9 => Some(if dtype.getattr("itemsize")?.extract::<usize>()? == 4 {
            "i32"
        } else {
            "i64"
        }),
        _ => None,
    };
    let target = match native {
        Some(_) => dtype,
        None => {
            let double = np.getattr("float64")?;
            if !np
                .call_method1("can_cast", (&dtype, &double, "safe"))?
                .is_truthy()?
            {
                return Err(conversion_error());
            }
            double
        }
    };
    let arr = np
        .call_method1("ascontiguousarray", (obj, target))
        .map_err(|_| conversion_error())?;
    let shape: Vec<usize> = arr.getattr("shape")?.extract()?;
    if shape.len() < min_dim || shape.len() > max_dim {
        return Err(conversion_error());
    }
    let flat = arr.call_method0("ravel")?.call_method0("tolist")?;
    let counts = match native {
        Some("f32") => Counts::F32(flat.extract::<Vec<f32>>()?),
        Some("i32") => Counts::I32(flat.extract()?),
        Some("i64") => Counts::I64(flat.extract()?),
        _ => Counts::F64(flat.extract()?),
    };
    Ok((counts, shape))
}

/// calculates the informational entropy of the values in an array
///
///   ARGUMENTS:
///
///     - resMat: pointer to a long int array containing the data
///     - dim: long int containing the length of the _tPtr_ array.
///
///   RETURNS:
///
///     a double
#[pyfunction]
#[pyo3(name = "InfoEntropy", signature = (resArr))]
#[allow(non_snake_case)]
fn info_entropy_py(py: Python<'_>, resArr: &Bound<'_, PyAny>) -> PyResult<f64> {
    let (counts, shape) = load_counts(py, resArr, 1, 1)?;
    // we are expecting a 1 dimensional array
    if shape[0] == 0 {
        return Err(to_py(RDKitError::invariant("")));
    }
    Ok(with_counts!(counts, |data| info_entropy(data)))
}

/// Calculates the information gain for a variable
///
///    ARGUMENTS:
///
///      - varMat: a Numeric Array object
///        varMat is a Numeric array with the number of possible occurrences
///          of each result for reach possible value of the given variable.
///
///        So, for a variable which adopts 4 possible values and a result which
///          has 3 possible values, varMat would be 4x3
///
///    RETURNS:
///
///      - a Python float object
///
///    NOTES
///
///      - this is a dropin replacement for _PyInfoGain()_ in entropy.py
#[pyfunction]
#[pyo3(name = "InfoGain", signature = (resArr))]
#[allow(non_snake_case)]
fn info_gain_py(py: Python<'_>, resArr: &Bound<'_, PyAny>) -> PyResult<f64> {
    let (counts, shape) = load_counts(py, resArr, 2, 2)?;
    Ok(with_counts!(counts, |data| info_entropy_gain(
        data, shape[0], shape[1]
    )))
}

/// Calculates the chi squared value for a variable
///
///    ARGUMENTS:
///
///      - varMat: a Numeric Array object
///        varMat is a Numeric array with the number of possible occurrences
///          of each result for reach possible value of the given variable.
///
///        So, for a variable which adopts 4 possible values and a result which
///          has 3 possible values, varMat would be 4x3
///
///    RETURNS:
///
///      - a Python float object
#[pyfunction]
#[pyo3(name = "ChiSquare", signature = (resArr))]
#[allow(non_snake_case)]
fn chi_square_py(py: Python<'_>, resArr: &Bound<'_, PyAny>) -> PyResult<f64> {
    let (counts, shape) = load_counts(py, resArr, 2, 2)?;
    Ok(with_counts!(counts, |data| chi_square(
        data, shape[0], shape[1]
    )))
}

/// Runs `$body` with `$bv` bound to the bit vector in `$obj`, whichever kind
/// it is.
macro_rules! with_bit_vect {
    ($obj:expr, $message:literal, |$bv:ident| $body:expr) => {{
        if let Ok(e) = $obj.cast::<PyExplicitBitVect>() {
            let $bv = &e.borrow().bv;
            $body
        } else if let Ok(sp) = $obj.cast::<PySparseBitVect>() {
            let $bv = &sp.borrow().bv;
            $body
        } else {
            Err(PyValueError::new_err($message))
        }
    }};
}

fn info_type_from(value: i64) -> PyResult<InfoType> {
    match value {
        1 => Ok(InfoType::Entropy),
        2 => Ok(InfoType::BiasEntropy),
        3 => Ok(InfoType::ChiSquare),
        4 => Ok(InfoType::BiasChiSquare),
        _ => Err(PyValueError::new_err("invalid InfoType value")),
    }
}

/// A class to rank the bits from a series of labelled fingerprints
/// A simple demonstration may help clarify what this class does.
/// Here's a small set of vectors:
///
/// >>> for i,bv in enumerate(bvs): print(bv.ToBitString(),acts[i])
/// ...
/// 0001 0
/// 0101 0
/// 0010 1
/// 1110 1
///
/// Default ranker, using infogain:
///
/// >>> ranker = InfoBitRanker(4,2)
/// >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
/// ...
/// >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
/// ...
/// 3 1.000 2 0
/// 2 1.000 0 2
/// 0 0.311 0 1
///
/// Using the biased infogain:
///
/// >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASENTROPY)
/// >>> ranker.SetBiasList((1,))
/// >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
/// ...
/// >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
/// ...
/// 2 1.000 0 2
/// 0 0.311 0 1
/// 1 0.000 1 1
///
/// A chi squared ranker is also available:
///
/// >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.CHISQUARE)
/// >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
/// ...
/// >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
/// ...
/// 3 4.000 2 0
/// 2 4.000 0 2
/// 0 1.333 0 1
///
/// As is a biased chi squared:
///
/// >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASCHISQUARE)
/// >>> ranker.SetBiasList((1,))
/// >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
/// ...
/// >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
/// ...
/// 2 4.000 0 2
/// 0 1.333 0 1
/// 1 0.000 1 1
#[pyclass(name = "InfoBitRanker", module = "rdoxide.ML.InfoTheory.rdInfoTheory")]
#[derive(Debug)]
struct PyInfoBitRanker {
    ranker: InfoBitRanker,
}

#[pymethods]
impl PyInfoBitRanker {
    #[new]
    #[pyo3(signature = (nBits, nClasses, infoType = None))]
    #[allow(non_snake_case)]
    fn new(nBits: u32, nClasses: u32, infoType: Option<i64>) -> PyResult<Self> {
        let info_type = infoType.map_or(Ok(InfoType::Entropy), info_type_from)?;
        Ok(PyInfoBitRanker {
            ranker: InfoBitRanker::new(nBits, nClasses, info_type),
        })
    }

    /// Accumulate the votes for all the bits turned on in a bit vector
    ///
    /// ARGUMENTS:
    ///
    ///   - bv : bit vector either ExplicitBitVect or SparseBitVect operator
    ///   - label : the class label for the bit vector. It is assumed that 0 <= class < nClasses
    #[pyo3(name = "AccumulateVotes", signature = (bitVect, label))]
    #[allow(non_snake_case)]
    fn accumulate_votes(&mut self, bitVect: &Bound<'_, PyAny>, label: u32) -> PyResult<()> {
        with_bit_vect!(
            bitVect,
            "Accumulate Vote can only take a explicitBitVects or SparseBitvects",
            |bv| self.ranker.accumulate_votes(bv, label).py()
        )
    }

    /// Set the classes to which the entropy calculation should be biased
    ///
    /// This list contains a set of class ids used when in the BIASENTROPY mode of ranking bits.
    /// In this mode, a bit must be correlated higher with one of the biased classes than all the
    /// other classes. For example, in a two class problem with actives and inactives, the fraction of
    /// actives that hit the bit has to be greater than the fraction of inactives that hit the bit
    ///
    /// ARGUMENTS:
    ///
    ///   - classList : list of class ids that we want a bias towards
    #[pyo3(name = "SetBiasList", signature = (classList))]
    #[allow(non_snake_case)]
    fn set_bias_list(&mut self, classList: Vec<i32>) -> PyResult<()> {
        self.ranker.set_bias_list(&classList).py()
    }

    /// Set the mask bits for the calculation
    ///
    /// ARGUMENTS:
    ///
    ///   - maskBits : list of mask bits to use
    #[pyo3(name = "SetMaskBits", signature = (maskBits))]
    #[allow(non_snake_case)]
    fn set_mask_bits(&mut self, maskBits: Vec<i32>) -> PyResult<()> {
        self.ranker.set_mask_bits(&maskBits).py()
    }

    /// Returns the top n bits ranked by the information metric
    /// This is actually the function where most of the work of ranking is happening
    ///
    /// ARGUMENTS:
    ///
    ///   - num : the number of top ranked bits that are required
    #[pyo3(name = "GetTopN", signature = (num))]
    fn get_top_n<'py>(&mut self, py: Python<'py>, num: u32) -> PyResult<Bound<'py, PyAny>> {
        let ncols = self.ranker.num_classes() + 2;
        let rows = self.ranker.get_top_n(num).py()?;
        float64_array(py, rows)?.call_method1("reshape", (num, ncols))
    }

    /// Write the bits that have been ranked to a file
    #[pyo3(name = "WriteTopBitsToFile", signature = (fileName))]
    #[allow(non_snake_case)]
    fn write_top_bits_to_file(&self, fileName: &str) -> PyResult<()> {
        self.ranker.write_top_bits_to_file(fileName).py()
    }

    #[pyo3(name = "Tester", signature = (bitVect))]
    #[allow(non_snake_case)]
    fn tester(&self, bitVect: &Bound<'_, PyAny>) {
        if let Ok(sp) = bitVect.cast::<PySparseBitVect>() {
            println!("Num of on bits: {}", sp.borrow().bv.num_on_bits());
        }
    }
}

/// A class to generate a pairwise correlation matrix between a list of bits
/// The mode of operation for this class is something like this
///
///    >>> cmg = BitCorrMatGenerator()
///    >>> cmg.SetBitList(blist)
///    >>> for fp in fpList:
///    >>>    cmg.CollectVotes(fp)
///    >>> corrMat = cmg.GetCorrMatrix()
///
///    The resulting correlation matrix is a one dimensional nummeric array containing the
///    lower triangle elements
#[pyclass(
    name = "BitCorrMatGenerator",
    module = "rdoxide.ML.InfoTheory.rdInfoTheory"
)]
#[derive(Debug)]
struct PyBitCorrMatGenerator {
    generator: BitCorrMatGenerator,
}

#[pymethods]
impl PyBitCorrMatGenerator {
    #[new]
    fn new() -> Self {
        PyBitCorrMatGenerator {
            generator: BitCorrMatGenerator::new(),
        }
    }

    /// Set the list of bits that need to be correllated
    ///
    ///  This may for example be their top ranking ensemble bits
    ///
    /// ARGUMENTS:
    ///
    ///   - bitList : an integer list of bit IDs
    #[pyo3(name = "SetBitList", signature = (bitList))]
    #[allow(non_snake_case)]
    fn set_bit_list(&mut self, bitList: Vec<i32>) {
        self.generator.set_bit_id_list(&bitList);
    }

    /// For each pair of on bits (bi, bj) in fp increase the correlation count for the pair by 1
    ///
    /// ARGUMENTS:
    ///
    ///   - fp : a bit vector to collect the fingerprints from
    #[pyo3(name = "CollectVotes", signature = (bitVect))]
    #[allow(non_snake_case)]
    fn collect_votes(&mut self, bitVect: &Bound<'_, PyAny>) -> PyResult<()> {
        with_bit_vect!(
            bitVect,
            "CollectVote can only take ExplicitBitVects or SparseBitVects",
            |bv| self.generator.collect_votes(bv).py()
        )
    }

    /// Get the correlation matrix following the collection of votes from a bunch of fingerprints
    #[pyo3(name = "GetCorrMatrix")]
    fn get_corr_matrix<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        float64_array(py, self.generator.corr_mat())
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing bunch of functions for information metrics and a ranker to rank bits",
    )?;
    add_enum_exported(
        m,
        "InfoType",
        &[
            ("ENTROPY", 1),
            ("BIASENTROPY", 2),
            ("CHISQUARE", 3),
            ("BIASCHISQUARE", 4),
        ],
        None,
    )?;
    m.add_class::<PyInfoBitRanker>()?;
    m.add_class::<PyBitCorrMatGenerator>()?;
    m.add_function(wrap_pyfunction!(info_entropy_py, m)?)?;
    m.add_function(wrap_pyfunction!(info_gain_py, m)?)?;
    m.add_function(wrap_pyfunction!(chi_square_py, m)?)
}
