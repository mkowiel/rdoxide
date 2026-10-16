//! The tautomer part of `rdoxide.Chem.MolStandardize.rdMolStandardize`
//! (upstream `GraphMol/MolStandardize/nbWrap/Tautomer.cpp`).

use std::sync::{Arc, Mutex};

use pyo3::exceptions::{PyAttributeError, PyIndexError, PyRuntimeError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::RDKitError;
use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::mol_standardize::tautomer::{
    SubstructTerm, TAUTOMER_SCORING_VERSION, Tautomer, TautomerEnumerator,
    TautomerEnumeratorCallback, TautomerEnumeratorResult, TautomerEnumeratorStatus,
    default_tautomer_score_substructs, ensure_symm_sssr, score_hetero_hs, score_rings,
    score_substructs, score_tautomer,
};

use crate::enums::{add_enum, to_py};
use crate::errors::{IntoPyResult, argument_error_in, boost_type_name, to_py as err_to_py};
use crate::mol::PyMol;
use crate::rdmolstandardize::PyCleanupParameters;
use crate::vector_suite::{VectorSuite, no_python_init, vector_suite_methods};

const STATUS: [(&str, i64); 4] = [
    ("Completed", 0),
    ("MaxTautomersReached", 1),
    ("MaxTransformsReached", 2),
    ("Canceled", 3),
];

/// Holds a Python exception raised by a callback until the enumeration
/// returns, so that it reaches the caller unchanged.
#[derive(Debug, Default, Clone)]
struct ErrSlot(Arc<Mutex<Option<PyErr>>>);

impl ErrSlot {
    fn stash(&self, e: PyErr) -> RDKitError {
        *self.0.lock().unwrap_or_else(|p| p.into_inner()) = Some(e);
        RDKitError::Runtime("exception in Python callback".into())
    }

    fn convert(&self, e: RDKitError) -> PyErr {
        self.0
            .lock()
            .unwrap_or_else(|p| p.into_inner())
            .take()
            .unwrap_or_else(|| err_to_py(e))
    }
}

fn mol_obj(py: Python<'_>, mol: &ROMol) -> PyResult<Py<PyMol>> {
    Py::new(py, PyMol { mol: mol.clone() })
}

/// A Python scoring function as a Rust one; raised exceptions (and
/// non-numeric results) are stashed in `slot`.
fn py_scorer(func: Py<PyAny>, slot: ErrSlot) -> impl Fn(&ROMol) -> rdoxide::Result<i32> + Sync {
    let func = Mutex::new(func);
    move |mol: &ROMol| {
        Python::attach(|py| {
            let f = func.lock().unwrap_or_else(|p| p.into_inner()).clone_ref(py);
            let res = mol_obj(py, mol)
                .and_then(|m| f.call1(py, (m,)))
                .and_then(|v| {
                    v.extract::<i32>(py).map_err(|_| {
                        PyTypeError::new_err("scoring function must return a numeric value")
                    })
                });
            res.map_err(|e| slot.stash(e))
        })
    }
}

// ------------------------------------------------------------------ results

/// used to hold the aromatic and kekulized versions of each tautomer
#[pyclass(
    name = "Tautomer",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyTautomer {
    inner: Tautomer,
}

#[pymethods]
impl PyTautomer {
    /// aromatic version of the tautomer
    #[getter]
    fn tautomer(&self) -> PyMol {
        PyMol {
            mol: (*self.inner.tautomer).clone(),
        }
    }

    /// kekulized version of the tautomer
    #[getter]
    fn kekulized(&self) -> PyResult<PyMol> {
        Ok(PyMol {
            mol: (*self.inner.kekulized().py()?).clone(),
        })
    }
}

/// maps SMILES strings to the respective Tautomer objects
#[pyclass(
    name = "SmilesTautomerMap",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PySmilesTautomerMap {
    res: Arc<TautomerEnumeratorResult>,
}

#[pymethods]
impl PySmilesTautomerMap {
    fn keys<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.res.tautomers.keys())
    }

    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let vals = self
            .res
            .tautomers
            .values()
            .map(|t| Py::new(py, PyTautomer { inner: t.clone() }))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, vals)
    }

    fn items<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let items = self
            .res
            .tautomers
            .iter()
            .map(|(k, t)| {
                let t = Py::new(py, PyTautomer { inner: t.clone() })?;
                PyTuple::new(
                    py,
                    [k.into_pyobject(py)?.into_any(), t.into_bound(py).into_any()],
                )
            })
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    fn __len__(&self) -> usize {
        self.res.tautomers.len()
    }
}

/// used to return tautomer enumeration results
#[pyclass(
    name = "TautomerEnumeratorResult",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyTautomerEnumeratorResult {
    res: Arc<TautomerEnumeratorResult>,
}

impl PyTautomerEnumeratorResult {
    fn mols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let mols = self
            .res
            .tautomers
            .values()
            .map(|t| mol_obj(py, &t.tautomer))
            .collect::<PyResult<Vec<_>>>()?;
        PyList::new(py, mols)
    }
}

#[pymethods]
impl PyTautomerEnumeratorResult {
    /// tautomers generated by the enumerator
    #[getter]
    fn tautomers<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        self.mols(py)
    }

    /// SMILES of tautomers generated by the enumerator
    #[getter]
    fn smiles(&self) -> Vec<String> {
        self.res.smiles()
    }

    /// dictionary mapping SMILES strings to the respective Tautomer objects
    #[getter(smilesTautomerMap)]
    fn smiles_tautomer_map(&self) -> PySmilesTautomerMap {
        PySmilesTautomerMap {
            res: self.res.clone(),
        }
    }

    /// whether the enumeration completed or not; see TautomerEnumeratorStatus for possible values
    #[getter]
    fn status(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "TautomerEnumeratorStatus", self.res.status as i64)
    }

    /// tuple of atom indices modified by the transforms
    #[getter(modifiedAtoms)]
    fn modified_atoms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.res.modified_atom_indices())
    }

    /// tuple of bond indices modified by the transforms
    #[getter(modifiedBonds)]
    fn modified_bonds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.res.modified_bond_indices())
    }

    /// tautomers generated by the enumerator
    fn __call__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        self.mols(py)
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.mols(py)?.into_any().try_iter().map(Bound::into_any)
    }

    fn __getitem__(&self, py: Python<'_>, pos: i64) -> PyResult<Py<PyMol>> {
        let sz = self.res.len() as i64;
        let pos = if pos < 0 { pos + sz } else { pos };
        if pos < 0 || pos >= sz {
            return Err(PyIndexError::new_err("index out of bounds"));
        }
        let mol = self.res.at(pos as usize).expect("in range");
        mol_obj(py, mol)
    }

    fn __len__(&self) -> usize {
        self.res.len()
    }
}

// ----------------------------------------------------------------- callback

/// Create a derived class from this abstract base class and
/// implement the __call__() method.
/// The __call__() method is called in the innermost loop of the
/// algorithm, and provides a mechanism to monitor or stop
/// its progress.
///
/// To have your callback called, pass an instance of your
/// derived class to TautomerEnumerator.SetCallback()
#[pyclass(
    name = "TautomerEnumeratorCallback",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    subclass
)]
#[derive(Debug, Default)]
pub struct PyTautomerEnumeratorCallback;

#[pymethods]
impl PyTautomerEnumeratorCallback {
    #[new]
    #[pyo3(signature = (*_args, **_kwargs))]
    fn new(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
        PyTautomerEnumeratorCallback
    }

    /// This must be implemented in the derived class. Return True if the
    /// tautomer enumeration should continue; False if the tautomer
    /// enumeration should stop.
    #[allow(unused_variables)]
    fn __call__(&self, mol: &Bound<'_, PyAny>, res: &Bound<'_, PyAny>) -> PyResult<bool> {
        Err(PyRuntimeError::new_err(
            "TautomerEnumeratorCallback.__call__(): tried to call a pure virtual function",
        ))
    }
}

// --------------------------------------------------------------- enumerator

#[pyclass(
    name = "TautomerEnumerator",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyTautomerEnumerator {
    inner: TautomerEnumerator,
    callback: Option<Py<PyAny>>,
    slot: ErrSlot,
}

impl PyTautomerEnumerator {
    fn wrap(inner: TautomerEnumerator) -> Self {
        PyTautomerEnumerator {
            inner,
            callback: None,
            slot: ErrSlot::default(),
        }
    }

    fn run<T>(&self, py: Python<'_>, f: impl FnOnce() -> rdoxide::Result<T> + Send) -> PyResult<T>
    where
        T: Send,
    {
        py.detach(f).map_err(|e| self.slot.convert(e))
    }
}

/// `TautomerEnumeratorResult` or an iterable of molecules.
enum Candidates {
    Result(Arc<TautomerEnumeratorResult>),
    Mols(Vec<ROMol>),
}

fn candidates(o: &Bound<'_, PyAny>) -> PyResult<Candidates> {
    if let Ok(r) = o.cast::<PyTautomerEnumeratorResult>() {
        return Ok(Candidates::Result(r.borrow().res.clone()));
    }
    let err = || PyTypeError::new_err("the passed object should be an iterable of Mol objects");
    let mut mols = Vec::new();
    for m in o.try_iter().map_err(|_| err())? {
        let m = m?;
        let m = m.cast::<PyMol>().map_err(|_| err())?;
        mols.push(m.borrow().mol.clone());
    }
    Ok(Candidates::Mols(mols))
}

#[pymethods]
impl PyTautomerEnumerator {
    /// `TautomerEnumerator()`, `TautomerEnumerator(params)` or
    /// `TautomerEnumerator(other)`.
    #[new]
    #[pyo3(signature = (params = None))]
    fn new(params: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let Some(p) = params else {
            return Ok(Self::wrap(TautomerEnumerator::default()));
        };
        if let Ok(other) = p.cast::<PyTautomerEnumerator>() {
            let other = other.borrow();
            return Ok(PyTautomerEnumerator {
                inner: other.inner.clone(),
                callback: other.callback.as_ref().map(|c| c.clone_ref(p.py())),
                slot: other.slot.clone(),
            });
        }
        let p = p
            .cast::<PyCleanupParameters>()
            .map_err(|_| PyTypeError::new_err("__init__(): incompatible function arguments"))?;
        Ok(Self::wrap(
            TautomerEnumerator::from_params(&p.borrow().params).py()?,
        ))
    }

    /// Generates the tautomers for a molecule.
    ///
    /// The enumeration rules are inspired by the publication:
    /// M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
    /// https://doi.org/10.1007/s10822-010-9346-4
    ///
    /// Note:
    /// - the definitions used here are that the atoms modified during
    ///   tautomerization are the atoms at the beginning and end of each tautomer
    ///   transform (the H "donor" and H "acceptor" in the transform) and the bonds
    ///   modified during transformation are any bonds whose order is changed during
    ///   the tautomer transform (these are the bonds between the "donor" and the
    ///   "acceptor").
    /// - if any atoms in the molecule have the "_protected" property set, they
    ///   will not be considered for tautomerization. Any transformations
    ///   involving these atoms will be skipped.
    #[pyo3(name = "Enumerate")]
    fn enumerate(
        &self,
        py: Python<'_>,
        mol: PyRef<'_, PyMol>,
    ) -> PyResult<PyTautomerEnumeratorResult> {
        let mol = &mol.mol;
        let res = self.run(py, || self.inner.enumerate(mol))?;
        Ok(PyTautomerEnumeratorResult { res: Arc::new(res) })
    }

    /// Returns the canonical tautomer for a molecule.
    ///
    /// The default scoring scheme is inspired by the publication:
    /// M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
    /// https://doi.org/10.1007/s10822-010-9346-4
    ///
    /// Note that the canonical tautomer is very likely not the most stable tautomer
    /// for any given conditions. The default scoring rules are designed to produce
    /// "reasonable" tautomers, but the primary concern is that the results are
    /// canonical: you always get the same canonical tautomer for a molecule
    /// regardless of what the input tautomer or atom ordering were.
    #[pyo3(name = "Canonicalize", signature = (mol, scoreFunc = None))]
    #[allow(non_snake_case)]
    fn canonicalize(
        &self,
        py: Python<'_>,
        mol: PyRef<'_, PyMol>,
        scoreFunc: Option<Py<PyAny>>,
    ) -> PyResult<PyMol> {
        let mol = &mol.mol;
        let res = match scoreFunc {
            None => self.run(py, || self.inner.canonicalize(mol, None))?,
            Some(f) => {
                let scorer = py_scorer(f, self.slot.clone());
                self.run(py, || self.inner.canonicalize(mol, Some(&scorer)))?
            }
        };
        Ok(PyMol { mol: res })
    }

    /// picks the canonical tautomer from an iterable of molecules
    #[pyo3(name = "PickCanonical", signature = (iterable, scoreFunc = None))]
    #[allow(non_snake_case)]
    fn pick_canonical(
        &self,
        py: Python<'_>,
        iterable: &Bound<'_, PyAny>,
        scoreFunc: Option<Py<PyAny>>,
    ) -> PyResult<PyMol> {
        let cands = candidates(iterable)?;
        let scorer = scoreFunc.map(|f| py_scorer(f, self.slot.clone()));
        let score = scorer
            .as_ref()
            .map(|s| s as &(dyn Fn(&ROMol) -> rdoxide::Result<i32> + Sync));
        let res = self.run(py, || match &cands {
            Candidates::Result(r) => self.inner.pick_canonical(r, score),
            Candidates::Mols(m) => self.inner.pick_canonical_mols(m, score),
        })?;
        Ok(PyMol { mol: res })
    }

    /// returns the score for a tautomer using the default scoring scheme.
    #[staticmethod]
    #[pyo3(name = "ScoreTautomer")]
    fn score_tautomer_py(mut mol: PyRefMut<'_, PyMol>) -> PyResult<i32> {
        // upstream leaves the symmetrized SSSR on the molecule
        ensure_symm_sssr(&mut mol.mol).py()?;
        score_tautomer(&mol.mol).py()
    }

    /// set the maximum number of tautomers to be generated.
    #[pyo3(name = "SetMaxTautomers")]
    #[allow(non_snake_case)]
    fn set_max_tautomers(&mut self, maxTautomers: u32) {
        self.inner.max_tautomers = maxTautomers;
    }

    /// returns the maximum number of tautomers to be generated.
    #[pyo3(name = "GetMaxTautomers")]
    fn get_max_tautomers(&self) -> u32 {
        self.inner.max_tautomers
    }

    /// set the maximum number of transformations to be applied. This limit is
    /// usually hit earlier than the maxTautomers limit and leads to a more linear
    /// scaling of CPU time with increasing number of tautomeric centers (see
    /// Sitzmann et al.).
    #[pyo3(name = "SetMaxTransforms")]
    #[allow(non_snake_case)]
    fn set_max_transforms(&mut self, maxTransforms: u32) {
        self.inner.max_transforms = maxTransforms;
    }

    /// returns the maximum number of transformations to be applied.
    #[pyo3(name = "GetMaxTransforms")]
    fn get_max_transforms(&self) -> u32 {
        self.inner.max_transforms
    }

    /// set to True if you wish stereochemistry information to be removed from
    /// sp3 atoms involved in tautomerism. This means that S-aminoacids will lose
    /// their stereochemistry after going through tautomer enumeration because of
    /// the amido-imidol tautomerism. This defaults to True in RDKit, and to False
    /// in the workflow described by Sitzmann et al.
    #[pyo3(name = "SetRemoveSp3Stereo")]
    #[allow(non_snake_case)]
    fn set_remove_sp3_stereo(&mut self, removeSp3Stereo: bool) {
        self.inner.remove_sp3_stereo = removeSp3Stereo;
    }

    /// returns whether stereochemistry information will be removed from sp3
    /// atoms involved in tautomerism.
    #[pyo3(name = "GetRemoveSp3Stereo")]
    fn get_remove_sp3_stereo(&self) -> bool {
        self.inner.remove_sp3_stereo
    }

    /// set to True if you wish stereochemistry information to be removed from
    /// double bonds involved in tautomerism. This means that enols will lose
    /// their E/Z stereochemistry after going through tautomer enumeration
    /// because of the keto-enolic tautomerism. This defaults to True in the RDKit
    /// and also in the workflow described by Sitzmann et al.
    #[pyo3(name = "SetRemoveBondStereo")]
    #[allow(non_snake_case)]
    fn set_remove_bond_stereo(&mut self, removeBondStereo: bool) {
        self.inner.remove_bond_stereo = removeBondStereo;
    }

    /// returns whether stereochemistry information will be removed from double
    /// bonds involved in tautomerism.
    #[pyo3(name = "GetRemoveBondStereo")]
    fn get_remove_bond_stereo(&self) -> bool {
        self.inner.remove_bond_stereo
    }

    /// set to True if you wish AssignStereochemistry to be called on each
    /// tautomer generated by the Enumerate() method. This defaults to True.
    #[pyo3(name = "SetReassignStereo")]
    #[allow(non_snake_case)]
    fn set_reassign_stereo(&mut self, reassignStereo: bool) {
        self.inner.reassign_stereo = reassignStereo;
    }

    /// returns whether AssignStereochemistry will be called on each tautomer
    /// generated by the Enumerate() method.
    #[pyo3(name = "GetReassignStereo")]
    fn get_reassign_stereo(&self) -> bool {
        self.inner.reassign_stereo
    }

    /// Pass an instance of a class derived from
    /// TautomerEnumeratorCallback, which must implement the
    /// __call__() method.
    #[pyo3(name = "SetCallback")]
    fn set_callback(&mut self, py: Python<'_>, callback: &Bound<'_, PyAny>) -> PyResult<()> {
        if callback.is_none() {
            self.callback = None;
            self.inner.set_callback(None);
            return Ok(());
        }
        if !callback.is_instance_of::<PyTautomerEnumeratorCallback>() {
            return Err(PyTypeError::new_err(
                "Expected an instance of a rdMolStandardize.TautomerEnumeratorCallback subclass",
            ));
        }
        let cls_dict = callback.get_type().getattr("__dict__")?;
        if !cls_dict.contains("__call__")? {
            return Err(PyAttributeError::new_err(
                "TautomerEnumeratorCallback subclass must override __call__",
            ));
        }
        if !cls_dict.get_item("__call__")?.is_callable() {
            return Err(PyAttributeError::new_err(
                "TautomerEnumeratorCallback.__call__ must be callable",
            ));
        }
        let obj = Mutex::new(callback.clone().unbind());
        let slot = self.slot.clone();
        let cb: Arc<TautomerEnumeratorCallback> = Arc::new(move |mol: &ROMol, res| {
            Python::attach(|py| {
                let f = obj.lock().unwrap_or_else(|p| p.into_inner()).clone_ref(py);
                let snapshot = PyTautomerEnumeratorResult {
                    res: Arc::new(res.clone()),
                };
                mol_obj(py, mol)
                    .and_then(|m| f.call1(py, (m, snapshot)))
                    .and_then(|v| v.bind(py).is_truthy())
                    .map_err(|e| slot.stash(e))
            })
        });
        self.inner.set_callback(Some(cb));
        self.callback = Some(callback.clone().unbind());
        let _ = py;
        Ok(())
    }

    /// Get the TautomerEnumeratorCallback subclass instance,
    /// or None if none was set.
    #[pyo3(name = "GetCallback")]
    fn get_callback(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.callback.as_ref().map(|c| c.clone_ref(py))
    }

    #[classattr]
    #[pyo3(name = "tautomerScoreVersion")]
    fn tautomer_score_version() -> &'static str {
        TAUTOMER_SCORING_VERSION
    }
}

/// return a TautomerEnumerator using v1 of the enumeration rules
#[pyfunction]
#[pyo3(name = "GetV1TautomerEnumerator")]
fn get_v1_tautomer_enumerator() -> PyTautomerEnumerator {
    PyTautomerEnumerator::wrap(TautomerEnumerator::v1())
}

/// scores the ring system of the tautomer for canonicalization
/// Aromatic rings score 100, all carbon aromatic rings score 250
#[pyfunction]
#[pyo3(name = "ScoreRings")]
fn score_rings_py(mut mol: PyRefMut<'_, PyMol>) -> PyResult<i32> {
    // upstream leaves the symmetrized SSSR on the molecule
    ensure_symm_sssr(&mut mol.mol).py()?;
    score_rings(&mol.mol).py()
}

/// scores the number of heteroHs of the tautomer for canonicalization
/// This gives a negative penalty to hydrogens attached to S,P, Se and Te
#[pyfunction]
#[pyo3(name = "ScoreHeteroHs")]
fn score_hetero_hs_py(mol: PyRef<'_, PyMol>) -> PyResult<i32> {
    score_hetero_hs(&mol.mol).py()
}

// ---------------------------------------------------------- substruct terms

/// Sets the score of this particular tautomer substructure, higher scores are more preferable
/// Aromatic rings score 100, all carbon aromatic rings score 250
#[pyclass(
    name = "SubstructTerm",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PySubstructTerm {
    inner: SubstructTerm,
}

#[pymethods]
impl PySubstructTerm {
    #[new]
    fn new(name: &str, smarts: &str, score: i32) -> Self {
        PySubstructTerm {
            inner: SubstructTerm::new(name, smarts, score, Vec::new(), ""),
        }
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[getter]
    fn smarts(&self) -> String {
        self.inner.smarts.clone()
    }

    #[getter]
    fn score(&self) -> i32 {
        self.inner.score
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .cast::<PySubstructTerm>()
            .is_ok_and(|o| o.borrow().inner == self.inner)
    }
}

/// A list of SubstructTerm objects (upstream's `vector_indexing_suite`).
#[pyclass(
    name = "SubstructTermVector",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    sequence
)]
#[derive(Debug, Default)]
pub struct PySubstructTermVector {
    terms: Vec<SubstructTerm>,
}

fn term_of(o: &Bound<'_, PyAny>) -> PyResult<SubstructTerm> {
    Ok(o.cast::<PySubstructTerm>()
        .map_err(|_| PyTypeError::new_err("expected a SubstructTerm"))?
        .borrow()
        .inner
        .clone())
}

impl VectorSuite for PySubstructTermVector {
    type Item = SubstructTerm;

    fn items(&self) -> &Vec<SubstructTerm> {
        &self.terms
    }

    fn items_mut(&mut self) -> &mut Vec<SubstructTerm> {
        &mut self.terms
    }

    fn to_py<'py>(py: Python<'py>, item: &SubstructTerm) -> PyResult<Bound<'py, PyAny>> {
        Ok(Bound::new(
            py,
            PySubstructTerm {
                inner: item.clone(),
            },
        )?
        .into_any())
    }

    fn from_py(obj: &Bound<'_, PyAny>) -> Option<SubstructTerm> {
        obj.cast::<PySubstructTerm>()
            .ok()
            .map(|t| t.borrow().inner.clone())
    }

    fn new_vector<'py>(py: Python<'py>, items: Vec<SubstructTerm>) -> PyResult<Bound<'py, PyAny>> {
        Ok(Bound::new(py, PySubstructTermVector { terms: items })?.into_any())
    }
}

#[pymethods]
impl PySubstructTermVector {
    #[new]
    #[pyo3(signature = (*args))]
    fn new(args: &Bound<'_, PyTuple>) -> PyResult<Self> {
        if !args.is_empty() {
            let got: Vec<String> = std::iter::once("SubstructTermVector".to_string())
                .chain(args.iter().map(|a| boost_type_name(&a)))
                .collect();
            return Err(argument_error_in(
                args.py(),
                "SubstructTermVector",
                "__init__",
                &got,
                "__init__(_object*)",
            ));
        }
        Ok(PySubstructTermVector::default())
    }
}

vector_suite_methods!(PySubstructTermVector);

no_python_init!(PyTautomer, PySmilesTautomerMap, PyTautomerEnumeratorResult);

/// scores the tautomer substructures
#[pyfunction]
#[pyo3(name = "ScoreSubstructs", signature = (mol, terms = None))]
fn score_substructs_py(mol: PyRef<'_, PyMol>, terms: Option<&Bound<'_, PyAny>>) -> PyResult<i32> {
    match terms {
        None => score_substructs(&mol.mol, default_tautomer_score_substructs()).py(),
        Some(t) => {
            let terms: Vec<SubstructTerm> = match t.cast::<PySubstructTermVector>() {
                Ok(v) => v.borrow().terms.clone(),
                Err(_) => t
                    .try_iter()?
                    .map(|x| term_of(&x?))
                    .collect::<PyResult<_>>()?,
            };
            score_substructs(&mol.mol, &terms).py()
        }
    }
}

/// Return the default tautomer substructure scoring terms
#[pyfunction]
#[pyo3(name = "GetDefaultTautomerScoreSubstructs")]
fn get_default_tautomer_score_substructs() -> PySubstructTermVector {
    PySubstructTermVector {
        terms: default_tautomer_score_substructs().to_vec(),
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum(m, "TautomerEnumeratorStatus", &STATUS, None)?;
    m.add_class::<PyTautomerEnumeratorCallback>()?;
    m.add_class::<PyTautomer>()?;
    m.add_class::<PySmilesTautomerMap>()?;
    m.add_class::<PyTautomerEnumeratorResult>()?;
    m.add_class::<PyTautomerEnumerator>()?;
    m.add_function(wrap_pyfunction!(get_v1_tautomer_enumerator, m)?)?;
    m.add_function(wrap_pyfunction!(score_rings_py, m)?)?;
    m.add_function(wrap_pyfunction!(score_hetero_hs_py, m)?)?;
    m.add_class::<PySubstructTerm>()?;
    m.add_class::<PySubstructTermVector>()?;
    m.add_function(wrap_pyfunction!(score_substructs_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_default_tautomer_score_substructs, m)?)?;
    let _ = TautomerEnumeratorStatus::Completed;
    Ok(())
}
