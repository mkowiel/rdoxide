//! `rdkit.Chem.rdEnumerateStereoisomers` (upstream
//! `GraphMol/EnumerateStereoisomers/Wrap/rdEnumerateStereoisomers.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::enumerate_stereoisomers::{
    StereoEnumerationOptions, StereoisomerEnumerator,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// EnumerateSteroisomers options.
#[pyclass(
    name = "StereoEnumerationOptions",
    module = "rdkit.Chem.rdEnumerateStereoisomers",
    skip_from_py_object
)]
#[derive(Clone, Default)]
pub struct PyStereoEnumerationOptions {
    pub(crate) opts: StereoEnumerationOptions,
}

#[pymethods]
impl PyStereoEnumerationOptions {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// If true, the process attempts to generate a standard RDKit distance geometry conformation for the stereoisomer.  If this fails, we assume that the stereoisomer is non-physical and don't return it.  NOTE that this is computationally expensive and is just a heuristic that could result in stereoisomers being lost.  Default=False
    #[getter(tryEmbedding)]
    fn try_embedding(&self) -> bool {
        self.opts.try_embedding
    }
    #[setter(tryEmbedding)]
    fn set_try_embedding(&mut self, v: bool) {
        self.opts.try_embedding = v;
    }

    /// If true, stereocenters which have a specified stereochemistry will not be perturbed unless they are part of a relative stereo group.  Default=True.
    #[getter(onlyUnassigned)]
    fn only_unassigned(&self) -> bool {
        self.opts.only_unassigned
    }
    #[setter(onlyUnassigned)]
    fn set_only_unassigned(&mut self, v: bool) {
        self.opts.only_unassigned = v;
    }

    /// If true, only find stereoisomers that differ at the StereoGroups associated with the molecule.  Default=False.
    #[getter(onlyStereoGroups)]
    fn only_stereo_groups(&self) -> bool {
        self.opts.only_stereo_groups
    }
    #[setter(onlyStereoGroups)]
    fn set_only_stereo_groups(&mut self, v: bool) {
        self.opts.only_stereo_groups = v;
    }

    /// If true, only stereoisomers that differ in canonical CXSmiles will be returned.  Default=True.
    #[getter]
    fn unique(&self) -> bool {
        self.opts.unique
    }
    #[setter(unique)]
    fn set_unique(&mut self, v: bool) {
        self.opts.unique = v;
    }

    /// The maximum number of isomers to yield.  If the number of possible isomers is greater than maxIsomers, a random subset will be yielded.  If 0, there is no maximum.  Since every additional stereocenter doubles the number of results (and execution time) it's important to keep an eye on this.
    #[getter(maxIsomers)]
    fn max_isomers(&self) -> u64 {
        self.opts.max_isomers
    }
    #[setter(maxIsomers)]
    fn set_max_isomers(&mut self, v: u64) {
        self.opts.max_isomers = v;
    }

    /// Seed for random number generator.  Default=-1 means no seed.
    #[getter(randomSeed)]
    fn random_seed(&self) -> i32 {
        self.opts.random_seed
    }
    #[setter(randomSeed)]
    fn set_random_seed(&mut self, v: i32) {
        self.opts.random_seed = v;
    }
}

/// Stereoisomer enumerator.
#[pyclass(
    name = "StereoisomerEnumerator",
    module = "rdkit.Chem.rdEnumerateStereoisomers"
)]
pub struct PyStereoisomerEnumerator {
    enumerator: StereoisomerEnumerator,
}

#[pymethods]
impl PyStereoisomerEnumerator {
    #[new]
    #[pyo3(signature = (mol, options = None, verbose = true))]
    fn new(
        mol: PyRef<'_, PyMol>,
        options: Option<PyRef<'_, PyStereoEnumerationOptions>>,
        verbose: bool,
    ) -> PyResult<Self> {
        let _ = verbose; // only controls logging of embedding failures
        let opts = options.map(|o| o.opts.clone()).unwrap_or_default();
        Ok(PyStereoisomerEnumerator {
            enumerator: StereoisomerEnumerator::new(&mol.mol, opts).py()?,
        })
    }

    /// Get next isomer in the sequence, or None if at the end.
    fn next(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyMol>>> {
        let isomer = py.detach(|| self.enumerator.next()).py()?;
        isomer.map(|mol| Py::new(py, PyMol { mol })).transpose()
    }

    /// Get the number of stereoisomers.
    #[pyo3(name = "GetStereoisomerCount")]
    fn stereoisomer_count(&self) -> u64 {
        self.enumerator.stereoisomer_count()
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to enumerate stereoisomers of a molecule.  Chiral centers and double bonds will be enumerated if unassigned, or, if the appropriate option is set, if assigned.  Atropisomers will only be enumerated if assigned.  There is, as yet, no means of finding  unassigned atropisomers.",
    )?;
    m.add_class::<PyStereoEnumerationOptions>()?;
    m.add_class::<PyStereoisomerEnumerator>()
}
