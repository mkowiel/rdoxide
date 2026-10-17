//! Reaction enumeration in `rdoxide.Chem.rdChemReactions` (upstream
//! `ChemReactions/Wrap/Enumerate.cpp`).
//!
//! The library and strategy classes keep upstream's hierarchy:
//! `EnumerateLibrary` derives from `EnumerateLibraryBase`, and the four
//! strategies derive from `EnumerationStrategyBase`. `GetEnumerator`,
//! `GetReaction` and `GetReagents` return copies where upstream returns
//! references. Serialization (`Serialize`, `InitFromString`, `GetState`,
//! `SetState`) needs boost archives, which are not provided; those methods
//! raise like an RDKit built without them.

use pyo3::exceptions::{PyStopIteration, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use rdoxide::RDKitError;
use rdoxide::graphmol::chem_reactions::enumerate::{
    Bbs, EnumerateLibrary, EnumerationParams, EnumerationStrategy, enumerate_library_can_serialize,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdchemreactions::{PyChemicalReaction, mol_ref, products_to_py, synced};

fn no_serialization() -> PyErr {
    crate::errors::to_py(RDKitError::precondition(
        "BOOST SERIALIZATION NOT INSTALLED",
    ))
}

/// `ConvertToVect`: a list or tuple of lists or tuples of molecules.
fn convert_bbs(reagents: &Bound<'_, PyAny>) -> PyResult<Bbs> {
    reagents
        .try_iter()?
        .map(|set| {
            set?.try_iter()?
                .map(|m| {
                    let m = m?;
                    mol_ref(&m).map(|m| m.mol.clone()).map_err(|_| {
                        PyValueError::new_err("reaction called with non molecule reactant")
                    })
                })
                .collect()
        })
        .collect()
}

/// Controls some aspects of how the enumeration is performed.
///
/// Options:
///   reagentMaxMatchCount [ default Infinite ]
///     This specifies how many times the reactant template can match a reagent.
///
///   sanePartialProducts [default false]
///     If true, forces all products of the reagent plus the product templates
///      pass chemical sanitization.  Note that if the product template itself
///      does not pass sanitization, then none of the products will.
#[pyclass(
    name = "EnumerationParams",
    module = "rdoxide.Chem.rdChemReactions",
    from_py_object
)]
#[derive(Clone)]
pub struct PyEnumerationParams {
    inner: EnumerationParams,
}

#[pymethods]
impl PyEnumerationParams {
    #[new]
    fn new() -> Self {
        PyEnumerationParams {
            inner: EnumerationParams::default(),
        }
    }

    #[getter(reagentMaxMatchCount)]
    fn reagent_max_match_count(&self) -> i32 {
        self.inner.reagent_max_match_count
    }

    #[setter(reagentMaxMatchCount)]
    fn set_reagent_max_match_count(&mut self, v: i32) {
        self.inner.reagent_max_match_count = v;
    }

    #[getter(sanePartialProducts)]
    fn sane_partial_products(&self) -> bool {
        self.inner.sane_partial_products
    }

    #[setter(sanePartialProducts)]
    fn set_sane_partial_products(&mut self, v: bool) {
        self.inner.sane_partial_products = v;
    }
}

// ------------------------------------------------------------ strategies

/// Base class of the enumeration strategies.
#[pyclass(
    name = "EnumerationStrategyBase",
    module = "rdoxide.Chem.rdChemReactions",
    subclass
)]
pub struct PyEnumerationStrategyBase {
    inner: EnumerationStrategy,
}

/// Builds the Python object of the strategy's own class.
fn strategy_to_py(py: Python<'_>, inner: EnumerationStrategy) -> PyResult<Py<PyAny>> {
    let base = PyEnumerationStrategyBase {
        inner: inner.clone(),
    };
    Ok(match inner.type_name() {
        "CartesianProductStrategy" => Py::new(
            py,
            PyClassInitializer::from(base).add_subclass(PyCartesianProductStrategy),
        )?
        .into_any(),
        "RandomSampleStrategy" => Py::new(
            py,
            PyClassInitializer::from(base).add_subclass(PyRandomSampleStrategy),
        )?
        .into_any(),
        "RandomSampleAllBBsStrategy" => Py::new(
            py,
            PyClassInitializer::from(base).add_subclass(PyRandomSampleAllBBsStrategy),
        )?
        .into_any(),
        _ => Py::new(
            py,
            PyClassInitializer::from(base).add_subclass(PyEvenSamplePairsStrategy),
        )?
        .into_any(),
    })
}

#[pymethods]
impl PyEnumerationStrategyBase {
    fn __bool__(&self) -> bool {
        self.inner.has_next()
    }

    /// Returns the enumeration strategy type as a string.
    #[pyo3(name = "Type")]
    fn type_name(&self) -> &'static str {
        self.inner.type_name()
    }

    /// Skip the next Nth results. note: this may be an expensive operation
    /// depending on the enumeration strategy used. It is recommended to use
    /// the enumerator state to advance to a known position
    #[pyo3(name = "Skip", signature = (skipCount))]
    #[allow(non_snake_case)]
    fn skip(&mut self, py: Python<'_>, skipCount: u64) -> PyResult<bool> {
        let inner = &mut self.inner;
        py.detach(|| inner.skip(skipCount)).py()?;
        Ok(true)
    }

    fn __copy__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        strategy_to_py(py, self.inner.clone())
    }

    /// Returns the total number of results for this enumeration strategy.
    /// Note that some strategies are effectively infinite.
    #[pyo3(name = "GetNumPermutations")]
    fn num_permutations(&self) -> u64 {
        self.inner.num_permutations()
    }

    /// Return the current indices into the arrays of reagents, as returned by
    /// GetReagents().  They do not necessarily refer to the input reagent
    /// sets as they only refer to reagents compatible with the reaction.
    #[pyo3(name = "GetPosition")]
    fn position(&self) -> Vec<u64> {
        self.inner.position().to_vec()
    }

    /// Return the next indices into the arrays of reagents
    fn next(&mut self) -> PyResult<Vec<u64>> {
        Ok(self.inner.advance().py()?.to_vec())
    }

    fn __next__(&mut self) -> PyResult<Vec<u64>> {
        self.next()
    }

    #[pyo3(name = "Initialize", signature = (rxn, ob))]
    fn initialize(
        &mut self,
        rxn: &Bound<'_, PyChemicalReaction>,
        ob: &Bound<'_, PyList>,
    ) -> PyResult<()> {
        let _ = rxn.borrow();
        let bbs = convert_bbs(ob.as_any())?;
        let sizes: Vec<u64> = bbs.iter().map(|b| b.len() as u64).collect();
        self.inner.initialize(&sizes).py()
    }
}

// rustfmt re-indents the attribute inside the macro on every pass
#[rustfmt::skip]
macro_rules! strategy_class {
    ($rust:ident, $py:literal, $ctor:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyclass(name = $py, module = "rdoxide.Chem.rdChemReactions", extends = PyEnumerationStrategyBase)]
        pub struct $rust;

        #[pymethods]
        impl $rust {
            #[new]
            fn new() -> PyClassInitializer<Self> {
                PyClassInitializer::from(PyEnumerationStrategyBase {
                    inner: EnumerationStrategy::$ctor(),
                })
                .add_subclass($rust)
            }
        }
    };
}

strategy_class!(
    PyCartesianProductStrategy,
    "CartesianProductStrategy",
    cartesian_product,
    "CartesianProductStrategy produces a standard walk through all possible\nreagent combinations:\n\n(0,0,0), (1,0,0), (2,0,0) ...\n"
);
strategy_class!(
    PyRandomSampleStrategy,
    "RandomSampleStrategy",
    random_sample,
    "RandomSampleStrategy simply randomly samples from the reagent sets.\nNote that this strategy never halts and can produce duplicates."
);
strategy_class!(
    PyRandomSampleAllBBsStrategy,
    "RandomSampleAllBBsStrategy",
    random_sample_all_bbs,
    "RandomSampleAllBBsStrategy randomly samples from the reagent sets\nwith the constraint that all building blocks are samples as early as possible.\nNote that this strategy never halts and can produce duplicates."
);
strategy_class!(
    PyEvenSamplePairsStrategy,
    "EvenSamplePairsStrategy",
    even_sample_pairs,
    "Randomly sample Pairs evenly from a collection of building blocks\nThis is a good strategy for choosing a relatively small selection\nof building blocks from a larger set.  As the amount of work needed\nto retrieve the next evenly sample building block grows with the\nnumber of samples, this method performs progressively worse as the\nnumber of samples gets larger.\nSee EnumerationStrategyBase for more details.\n"
);

#[pymethods]
impl PyEvenSamplePairsStrategy {
    /// Return the statistics log of the pairs used in the current enumeration.
    #[pyo3(name = "Stats")]
    fn stats(slf: PyRef<'_, Self>) -> String {
        slf.as_super().inner.stats().unwrap_or_default()
    }
}

// --------------------------------------------------------------- library

/// Base class of the reaction enumeration libraries.
#[pyclass(
    name = "EnumerateLibraryBase",
    module = "rdoxide.Chem.rdChemReactions",
    subclass
)]
pub struct PyEnumerateLibraryBase {
    lib: EnumerateLibrary,
}

#[pymethods]
impl PyEnumerateLibraryBase {
    fn __bool__(&self) -> PyResult<bool> {
        self.lib.has_next().py()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Return the next molecule from the enumeration.
    fn next(&mut self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        if !self.lib.has_next().py()? {
            return Err(PyStopIteration::new_err("Enumerations exhausted"));
        }
        let lib = &mut self.lib;
        let products = py.detach(|| lib.next_products()).py()?;
        products_to_py(py, products)
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.next(py)
    }

    /// Return the next smiles string from the enumeration.
    #[pyo3(name = "nextSmiles")]
    fn next_smiles(&mut self, py: Python<'_>) -> PyResult<Vec<Vec<String>>> {
        let lib = &mut self.lib;
        py.detach(|| lib.next_smiles()).py()
    }

    /// Serialize the library to a binary string.
    /// Note that the position in the library is serialized as well.  Care should
    /// be taken when serializing.  See GetState/SetState for position manipulation.
    #[pyo3(name = "Serialize")]
    fn serialize(&self) -> PyResult<()> {
        Err(no_serialization())
    }

    /// Inititialize the library from a binary string
    #[pyo3(name = "InitFromString", signature = (data))]
    fn init_from_string(&mut self, data: &Bound<'_, PyAny>) -> PyResult<()> {
        let _ = data;
        Err(no_serialization())
    }

    /// Returns the current enumeration position into the reagent vectors, as
    /// returned by GetReagents().  They do not necessarily refer to the input
    /// reagent sets as they only refer to reagents compatible with the reaction.
    #[pyo3(name = "GetPosition")]
    fn position(&self) -> PyResult<Vec<u64>> {
        Ok(self.lib.position().py()?.to_vec())
    }

    /// Returns the current enumeration state (position) of the library.
    /// This position can be used to restart the library from a known position
    #[pyo3(name = "GetState")]
    fn state(&self) -> PyResult<()> {
        Err(no_serialization())
    }

    /// Sets the enumeration state (position) of the library.
    #[pyo3(name = "SetState", signature = (state))]
    fn set_state(&mut self, state: &Bound<'_, PyAny>) -> PyResult<()> {
        let _ = state;
        Err(no_serialization())
    }

    /// Returns the current enumeration state (position) of the library to the start.
    #[pyo3(name = "ResetState")]
    fn reset_state(&mut self) -> PyResult<()> {
        self.lib.reset_state().py()
    }

    /// Returns the chemical reaction for this library
    #[pyo3(name = "GetReaction")]
    fn reaction(&self, py: Python<'_>) -> PyResult<Py<PyChemicalReaction>> {
        Py::new(py, PyChemicalReaction::new(self.lib.reaction().clone()))
    }

    /// Returns the enumation strategy for the current library
    #[pyo3(name = "GetEnumerator")]
    fn enumerator(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        strategy_to_py(py, self.lib.enumerator().py()?.clone())
    }
}

/// This class allows easy enumeration of reactions.  Simply provide a reaction
/// and a set of reagents and you are off the races.
///
/// Note that this functionality should be considered beta and that the API may
/// change in a future release.
///
/// EnumerateLibrary follows the python enumerator protocol, for example:
///
/// library = EnumerateLibrary(rxn, bbs)
/// for products in library:
///    ... do something with the product
///
/// It is useful to sanitize reactions before hand:
///
/// SanitizeRxn(rxn)
/// library = EnumerateLibrary(rxn, bbs)
///
/// Each product has the same output as applying a set of reagents to
/// the libraries reaction.
///
/// The returned data structure is as follows:
///
///    [ [products1], [products2],... ]
/// Where products1 are the molecule products for the reactions first product
/// template and products2 are the molecule products for the second product
/// template.  Since each reactant can match more than once, there may be
/// multiple product molecules for each template.
///
/// To help with large libraries, you can supply an enumeration strategy.
/// The default strategy is a CartesianProductStrategy which enumerates
/// everything.  RandomSampleStrategy randomly samples the products but
/// this strategy never terminates, however, python supplies itertools:
///
/// import itertools
/// library = EnumerateLibrary(rxn, bbs, rdChemReactions.RandomSampleStrategy())
/// for result in itertools.islice(library, 1000):
///     # do something with the first 1000 samples
#[pyclass(
    name = "EnumerateLibrary",
    module = "rdoxide.Chem.rdChemReactions",
    extends = PyEnumerateLibraryBase
)]
pub struct PyEnumerateLibrary;

#[pymethods]
impl PyEnumerateLibrary {
    #[new]
    #[pyo3(signature = (rxn = None, reagents = None, enumerator = None, params = None))]
    fn new(
        rxn: Option<&Bound<'_, PyChemicalReaction>>,
        reagents: Option<&Bound<'_, PyAny>>,
        enumerator: Option<&Bound<'_, PyAny>>,
        params: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyClassInitializer<Self>> {
        let Some(rxn) = rxn else {
            return Ok(PyClassInitializer::from(PyEnumerateLibraryBase {
                lib: EnumerateLibrary::empty(),
            })
            .add_subclass(PyEnumerateLibrary));
        };
        let reagents = reagents.ok_or_else(|| PyValueError::new_err("reagents are required"))?;
        // (rxn, reagents, params) and (rxn, reagents, enumerator[, params])
        let (enumerator, params) = match (enumerator, params) {
            (Some(e), None) if e.extract::<PyEnumerationParams>().is_ok() => (None, Some(e)),
            (e, p) => (e, p),
        };
        let strategy = enumerator
            .map(|e| {
                e.extract::<PyRef<'_, PyEnumerationStrategyBase>>()
                    .map(|s| s.inner.clone())
            })
            .transpose()?;
        let params = params
            .map(|p| p.extract::<PyEnumerationParams>().map(|p| p.inner))
            .transpose()?
            .unwrap_or_default();
        let bbs = convert_bbs(reagents)?;
        let rxn = synced(rxn).rxn.clone();
        let lib = EnumerateLibrary::new(&rxn, &bbs, strategy.as_ref(), &params).py()?;
        Ok(PyClassInitializer::from(PyEnumerateLibraryBase { lib })
            .add_subclass(PyEnumerateLibrary))
    }

    /// Return the reagents used in this library.  These are the subset of the
    /// input reagents that are compatible with the reaction so may be smaller
    /// than the input reagent sets.
    #[pyo3(name = "GetReagents")]
    fn reagents(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<PyList>> {
        let sets = slf
            .as_super()
            .lib
            .reagents()
            .iter()
            .map(|set| {
                let mols = set
                    .iter()
                    .map(|m| Py::new(py, PyMol { mol: m.clone() }))
                    .collect::<PyResult<Vec<_>>>()?;
                PyList::new(py, mols)
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, sets)?.unbind())
    }
}

/// Returns True if the EnumerateLibrary is serializable (requires boost serialization
#[pyfunction]
#[pyo3(name = "EnumerateLibraryCanSerialize")]
fn py_can_serialize() -> bool {
    enumerate_library_can_serialize()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyEnumerateLibraryBase>()?;
    m.add_class::<PyEnumerationParams>()?;
    m.add_class::<PyEnumerateLibrary>()?;
    m.add_class::<PyEnumerationStrategyBase>()?;
    m.add_class::<PyCartesianProductStrategy>()?;
    m.add_class::<PyRandomSampleStrategy>()?;
    m.add_class::<PyRandomSampleAllBBsStrategy>()?;
    m.add_class::<PyEvenSamplePairsStrategy>()?;
    m.add_function(wrap_pyfunction!(py_can_serialize, m)?)
}
