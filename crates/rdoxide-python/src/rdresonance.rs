//! `ResonanceMolSupplier` and `ResonanceMolSupplierCallback` of
//! `rdkit.Chem.rdchem` (upstream `Wrap/ResonanceMolSupplier.cpp`).

use std::sync::{Arc, Mutex};

use pyo3::exceptions::{
    PyAttributeError, PyIndexError, PyRuntimeError, PyStopIteration, PyTypeError,
};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::resonance::{ResonanceMolSupplier, ResonanceProgress};
use rdoxide::graphmol::substruct::SubstructMatchParameters;

use crate::errors::to_py;
use crate::mol::PyMol;

const PURE_VIRTUAL: &str = "Pure virtual function called";

/// Create a derived class from this abstract base class and
/// implement the __call__() method.
/// The __call__() method is called at each iteration of the
/// algorithm, and provides a mechanism to monitor or stop
/// its progress.
///
/// To have your callback called, pass an instance of your
/// derived class to ResonanceMolSupplier.SetProgressCallback()
#[pyclass(
    name = "ResonanceMolSupplierCallback",
    module = "rdkit.Chem.rdchem",
    subclass
)]
pub struct PyResonanceMolSupplierCallback {
    progress: ResonanceProgress,
}

#[pymethods]
impl PyResonanceMolSupplierCallback {
    // subclasses pass their own constructor arguments through `__new__`
    #[new]
    #[pyo3(signature = (*_args, **_kwargs))]
    fn new(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
        PyResonanceMolSupplierCallback {
            progress: ResonanceProgress::default(),
        }
    }

    /// Returns the number of individual conjugated groups in the molecule.
    #[pyo3(name = "GetNumConjGrps")]
    fn get_num_conj_grps(&self) -> usize {
        self.progress.groups.len()
    }

    /// Get the number of conjugated groups this molecule has.
    #[pyo3(name = "GetMaxStructures")]
    fn get_max_structures(&self) -> usize {
        self.progress.max_structs
    }

    /// Get the number of resonance structures generated so far for the passed
    /// conjugated group index.
    #[pyo3(name = "GetNumStructures", signature = (conjGrpIdx))]
    #[allow(non_snake_case)]
    fn get_num_structures(&self, conjGrpIdx: usize) -> PyResult<usize> {
        self.group(conjGrpIdx).map(|g| g.0)
    }

    /// Get the number of non-degenrate resonance structures generated so far
    /// for the passed conjugated group index.
    #[pyo3(name = "GetNumDiverseStructures", signature = (conjGrpIdx))]
    #[allow(non_snake_case)]
    fn get_num_diverse_structures(&self, conjGrpIdx: usize) -> PyResult<usize> {
        self.group(conjGrpIdx).map(|g| g.1)
    }

    /// This must be implemented in the derived class. Return True if the
    /// resonance structure generation should continue; False if the resonance
    /// structure generation should stop.
    fn __call__(&self) -> PyResult<bool> {
        Err(PyRuntimeError::new_err(PURE_VIRTUAL))
    }
}

impl PyResonanceMolSupplierCallback {
    fn group(&self, idx: usize) -> PyResult<(usize, usize)> {
        self.progress.groups.get(idx).copied().ok_or_else(|| {
            PyRuntimeError::new_err("Pre-condition Violation\n\tconjGrpIdx out of bounds")
        })
    }
}

/// A class which supplies resonance structures (as mols) from a mol.
///
///   Usage examples:
///
///     1) Lazy evaluation: the resonance structures are not constructed
///        until we ask for them:
///
///        >>> suppl = ResonanceMolSupplier(mol)
///        >>> for resMol in suppl:
///        ...    resMol.GetNumAtoms()
///
///     2) Lazy evaluation 2:
///
///        >>> suppl = ResonanceMolSupplier(mol)
///        >>> resMol1 = next(suppl)
///        >>> resMol2 = next(suppl)
///        >>> suppl.reset()
///        >>> resMol3 = next(suppl)
///        # resMol3 and resMol1 are the same:
///        >>> MolToSmiles(resMol3)==MolToSmiles(resMol1)
///
///     3) Random Access:
///
///        >>> suppl = ResonanceMolSupplier(mol)
///        >>> resMol1 = suppl[0]
///        >>> resMol2 = suppl[1]
///
///        NOTE: this will generate an IndexError if the supplier doesn't have that many
///        molecules.
///
///     4) Random Access 2: looping over all resonance structures
///        >>> suppl = ResonanceMolSupplier(mol)
///        >>> nResMols = len(suppl)
///        >>> for i in range(nResMols):
///        ...   suppl[i].GetNumAtoms()
#[pyclass(name = "ResonanceMolSupplier", module = "rdkit.Chem.rdchem", subclass)]
pub struct PyResonanceMolSupplier {
    inner: ResonanceMolSupplier,
    callback: Option<Py<PyAny>>,
    /// an exception raised by the callback, re-raised once the enumeration
    /// has stopped
    error: Arc<Mutex<Option<PyErr>>>,
}

impl PyResonanceMolSupplier {
    fn enumerate(&mut self, py: Python<'_>) -> PyResult<()> {
        let inner = &mut self.inner;
        let res = py.detach(|| inner.enumerate());
        if let Some(e) = self.error.lock().unwrap_or_else(|p| p.into_inner()).take() {
            return Err(e);
        }
        res.map_err(to_py)
    }

    fn length(&mut self, py: Python<'_>) -> PyResult<u32> {
        self.enumerate(py)?;
        self.inner.length().map_err(to_py)
    }

    fn item(&mut self, py: Python<'_>, idx: u32) -> PyResult<PyMol> {
        self.enumerate(py)?;
        let inner = &mut self.inner;
        py.detach(|| inner.get(idx))
            .map(|mol| PyMol { mol })
            .map_err(to_py)
    }
}

#[pymethods]
impl PyResonanceMolSupplier {
    #[new]
    #[pyo3(signature = (mol, flags = 0, maxStructs = 1000))]
    #[allow(non_snake_case)]
    fn new(mol: &PyMol, flags: u32, maxStructs: u32) -> PyResult<Self> {
        Ok(PyResonanceMolSupplier {
            inner: ResonanceMolSupplier::new(&mol.mol, flags, maxStructs).map_err(to_py)?,
            callback: None,
            error: Arc::default(),
        })
    }

    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyRefMut<'_, Self>> {
        let py = slf.py();
        slf.reset(py)?;
        Ok(slf)
    }

    /// Returns the next resonance structure in the supplier. Raises
    /// _StopIteration_ on end.
    fn __next__(&mut self, py: Python<'_>) -> PyResult<PyMol> {
        self.enumerate(py)?;
        if self.inner.at_end().map_err(to_py)? {
            return Err(PyStopIteration::new_err("End of supplier hit"));
        }
        let inner = &mut self.inner;
        py.detach(|| inner.next_mol())
            .map_err(to_py)?
            .map(|mol| PyMol { mol })
            .ok_or_else(|| PyStopIteration::new_err("End of supplier hit"))
    }

    fn __getitem__(&mut self, py: Python<'_>, idx: i64) -> PyResult<PyMol> {
        let len = i64::from(self.length(py)?);
        let idx = if idx < 0 { idx + len } else { idx };
        if idx < 0 || idx >= len {
            return Err(PyIndexError::new_err("invalid index"));
        }
        self.item(py, idx as u32)
    }

    fn __len__(&mut self, py: Python<'_>) -> PyResult<usize> {
        Ok(self.length(py)? as usize)
    }

    /// Resets our position in the resonance structure supplier to the
    /// beginning.
    fn reset(&mut self, py: Python<'_>) -> PyResult<()> {
        self.enumerate(py)?;
        self.inner.reset().map_err(to_py)
    }

    /// Returns whether or not we have hit the end of the resonance structure
    /// supplier.
    #[pyo3(name = "atEnd")]
    fn at_end(&mut self, py: Python<'_>) -> PyResult<bool> {
        self.enumerate(py)?;
        self.inner.at_end().map_err(to_py)
    }

    /// Returns the number of individual conjugated groups in the molecule.
    #[pyo3(name = "GetNumConjGrps")]
    fn get_num_conj_grps(&self) -> u32 {
        self.inner.get_num_conj_grps()
    }

    /// Given a bond index, it returns the index of the conjugated groupthe bond
    /// belongs to, or -1 if it is not conjugated.
    #[pyo3(name = "GetBondConjGrpIdx", signature = (bi))]
    fn get_bond_conj_grp_idx(&self, bi: u32) -> PyResult<i32> {
        self.inner.get_bond_conj_grp_idx(bi).map_err(to_py)
    }

    /// Given an atom index, it returns the index of the conjugated groupthe atom
    /// belongs to, or -1 if it is not conjugated.
    #[pyo3(name = "GetAtomConjGrpIdx", signature = (ai))]
    fn get_atom_conj_grp_idx(&self, ai: u32) -> PyResult<i32> {
        self.inner.get_atom_conj_grp_idx(ai).map_err(to_py)
    }

    /// Sets the number of threads to be used to enumerate resonance
    /// structures (defaults to 1; 0 selects the number of concurrent
    /// threads supported by the hardware; negative values are added
    /// to the number of concurrent threads supported by the hardware).
    #[pyo3(name = "SetNumThreads", signature = (numThreads))]
    #[allow(non_snake_case)]
    fn set_num_threads(&mut self, numThreads: i32) {
        self.inner.set_num_threads(numThreads);
    }

    /// Pass an instance of a class derived from
    /// ResonanceMolSupplierCallback, which must implement the
    /// __call__() method.
    #[pyo3(name = "SetProgressCallback", signature = (callback))]
    fn set_progress_callback(&mut self, callback: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = callback.py();
        if callback.is_none() {
            self.inner.set_progress_callback(None);
            self.callback = None;
            return Ok(());
        }
        let base = py.get_type::<PyResonanceMolSupplierCallback>();
        if !callback.is_instance(&base)? {
            return Err(PyTypeError::new_err(
                "Expected an instance of a rdchem.ResonanceMolSupplierCallback subclass",
            ));
        }
        let overridden = callback.get_type().getattr("__call__")?;
        if overridden.is(&base.getattr("__call__")?) || !overridden.is_callable() {
            return Err(PyAttributeError::new_err(
                "The __call__ attribute in the rdchem.ResonanceMolSupplierCallback subclass \
                 must exist and be a callable method",
            ));
        }
        let obj: Py<PyAny> = callback.clone().unbind();
        let held = obj.clone_ref(py);
        let error = Arc::clone(&self.error);
        self.inner
            .set_progress_callback(Some(Box::new(move |progress| {
                Python::attach(|py| {
                    let outcome = (|| -> PyResult<bool> {
                        let cb = held.bind(py);
                        cb.cast::<PyResonanceMolSupplierCallback>()?
                            .borrow_mut()
                            .progress = progress.clone();
                        cb.call_method0("__call__")?.is_truthy()
                    })();
                    outcome.unwrap_or_else(|e| {
                        *error.lock().unwrap_or_else(|p| p.into_inner()) = Some(e);
                        false
                    })
                })
            })));
        self.callback = Some(obj);
        Ok(())
    }

    /// Get the ResonanceMolSupplierCallback subclass instance, or None if none
    /// was set.
    #[pyo3(name = "GetProgressCallback")]
    fn get_progress_callback(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.callback.as_ref().map(|c| c.clone_ref(py))
    }

    /// Returns True if the resonance structure generation was canceled.
    #[pyo3(name = "WasCanceled")]
    fn was_canceled(&self) -> bool {
        self.inner.was_canceled()
    }

    /// Ask ResonanceMolSupplier to enumerate resonance structures
    /// (automatically done as soon as any attempt to access them is made).
    #[pyo3(name = "Enumerate")]
    fn py_enumerate(&mut self, py: Python<'_>) -> PyResult<()> {
        self.enumerate(py)
    }

    /// Returns true if resonance structure enumeration has already happened.
    #[pyo3(name = "GetIsEnumerated")]
    fn get_is_enumerated(&self) -> bool {
        self.inner.is_enumerated()
    }

    /// Returns the indices of the molecule's atoms that match a substructure
    /// query, taking into account all resonance structures in
    /// ResonanceMolSupplier.
    ///
    ///   ARGUMENTS:
    ///     - query: a Molecule
    ///
    ///     - useChirality: enables the use of stereochemistry in the matching
    ///
    ///     - useQueryQueryMatches: use query-query matching logic
    ///
    ///   RETURNS: a tuple of integers
    ///
    ///   NOTES:
    ///      - only a single match is returned
    ///      - the ordering of the indices corresponds to the atom ordering
    ///          in the query. For example, the first index is for the atom in
    ///          this molecule that matches the first atom in the query.
    #[pyo3(
        name = "GetSubstructMatch",
        signature = (query, useChirality = false, useQueryQueryMatches = false)
    )]
    #[allow(non_snake_case)]
    fn get_substruct_match<'py>(
        &mut self,
        py: Python<'py>,
        query: &PyMol,
        useChirality: bool,
        useQueryQueryMatches: bool,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let params = SubstructMatchParameters {
            recursion_possible: true,
            use_chirality: useChirality,
            use_query_query_matches: useQueryQueryMatches,
            max_matches: 1,
            ..SubstructMatchParameters::default()
        };
        let matches = self.matches(py, query, params)?;
        let first: Vec<i32> = matches
            .first()
            .map(|m| m.iter().map(|&(_, a)| a).collect())
            .unwrap_or_default();
        PyTuple::new(py, first)
    }

    /// Returns tuples of the indices of the molecule's atoms that match a
    /// substructure query, taking into account all resonance structures in
    /// ResonanceMolSupplier.
    ///
    ///   ARGUMENTS:
    ///     - query: a Molecule.
    ///     - uniquify: (optional) determines whether or not the matches are uniquified.
    ///                 Defaults to 1.
    ///
    ///     - useChirality: enables the use of stereochemistry in the matching
    ///
    ///     - useQueryQueryMatches: use query-query matching logic
    ///
    ///     - maxMatches: The maximum number of matches that will be returned.
    ///                   In high-symmetry cases with medium-sized molecules, it is
    ///                   very easy to end up with a combinatorial explosion in the
    ///                   number of possible matches. This argument prevents that from
    ///                   having unintended consequences
    ///
    ///     - numThreads: The number of threads to be used (defaults to 1; 0 selects the
    ///                   number of concurrent threads supported by the hardware; negative
    ///                   values are added to the number of concurrent threads supported
    ///                   by the hardware).
    ///
    ///   RETURNS: a tuple of tuples of integers
    ///
    ///   NOTE:
    ///      - the ordering of the indices corresponds to the atom ordering
    ///          in the query. For example, the first index is for the atom in
    ///          this molecule that matches the first atom in the query.
    #[pyo3(
        name = "GetSubstructMatches",
        signature = (
            query,
            uniquify = false,
            useChirality = false,
            useQueryQueryMatches = false,
            maxMatches = 1000,
            numThreads = 1
        )
    )]
    #[allow(non_snake_case)]
    fn get_substruct_matches<'py>(
        &mut self,
        py: Python<'py>,
        query: &PyMol,
        uniquify: bool,
        useChirality: bool,
        useQueryQueryMatches: bool,
        maxMatches: u32,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let params = SubstructMatchParameters {
            uniquify,
            recursion_possible: true,
            use_chirality: useChirality,
            use_query_query_matches: useQueryQueryMatches,
            max_matches: maxMatches,
            num_threads: numThreads,
            ..SubstructMatchParameters::default()
        };
        let matches = self.matches(py, query, params)?;
        PyTuple::new(
            py,
            matches
                .iter()
                .map(|m| PyTuple::new(py, m.iter().map(|&(_, a)| a)))
                .collect::<PyResult<Vec<_>>>()?,
        )
    }
}

impl PyResonanceMolSupplier {
    fn matches(
        &mut self,
        py: Python<'_>,
        query: &PyMol,
        params: SubstructMatchParameters,
    ) -> PyResult<Vec<Vec<(i32, i32)>>> {
        self.enumerate(py)?;
        let inner = &mut self.inner;
        py.detach(|| inner.substruct_matches(&query.mol, &params))
            .map_err(to_py)
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyResonanceMolSupplierCallback>()?;
    m.add_class::<PyResonanceMolSupplier>()?;
    Ok(())
}
