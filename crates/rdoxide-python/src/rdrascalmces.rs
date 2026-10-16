//! `rdkit.Chem.rdRascalMCES`: the RASCAL maximum common edge substructure and
//! clustering (upstream `GraphMol/RascalMCES/Wrap/rdRascalMCES.cpp`).

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;
use rdoxide::graphmol::rascal_mces::{
    RascalClusterOptions, RascalOptions, RascalResult, rascal_butina_cluster, rascal_cluster,
    rascal_mces,
};
use rdoxide::graphmol::romol::ROMol;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Used to return RASCAL MCES results.
#[pyclass(name = "RascalResult", module = "rdkit.Chem.rdRascalMCES")]
pub struct PyRascalResult {
    res: RascalResult,
}

#[pymethods]
impl PyRascalResult {
    /// SMARTS string defining the MCES.
    #[getter(smartsString)]
    fn smarts_string(&self) -> PyResult<String> {
        self.res.smarts().py()
    }

    /// A function returning a list of list of tuples, each inner list
    /// containing the matching bonds in the MCES as tuples of bond indices
    /// from mol1 and mol2
    #[pyo3(name = "bondMatches")]
    fn bond_matches(&self) -> Vec<(i32, i32)> {
        self.res.bond_matches().to_vec()
    }

    /// Likewise for atoms.
    #[pyo3(name = "atomMatches")]
    fn atom_matches(&self) -> Vec<(i32, i32)> {
        self.res.atom_matches().to_vec()
    }

    /// Function that cuts the MCES down to the single largest frag.  This
    /// cannot be undone.
    #[pyo3(name = "largestFragmentOnly")]
    fn largest_fragment_only(&mut self) -> PyResult<()> {
        self.res.largest_frag_only().py()
    }

    /// Johnson similarity between 2 molecules.
    #[getter]
    fn similarity(&self) -> f64 {
        self.res.similarity()
    }

    /// Number of fragments in MCES.
    #[getter(numFragments)]
    fn num_fragments(&self) -> PyResult<i32> {
        self.res.num_frags().py()
    }

    /// Number of atoms in largest fragment.
    #[getter(largestFragmentSize)]
    fn largest_fragment_size(&self) -> PyResult<u32> {
        self.res.largest_frag_size().py()
    }

    /// The tier 1 similarity estimate.
    #[getter(tier1Sim)]
    fn tier1_sim(&self) -> f64 {
        self.res.tier1_sim()
    }

    /// The tier 2 similarity estimate.
    #[getter(tier2Sim)]
    fn tier2_sim(&self) -> f64 {
        self.res.tier2_sim()
    }

    /// Whether it timed out.
    #[getter(timedOut)]
    fn timed_out(&self) -> bool {
        self.res.timed_out()
    }
}

/// RASCAL Options
#[pyclass(
    name = "RascalOptions",
    module = "rdkit.Chem.rdRascalMCES",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyRascalOptions {
    pub(crate) opts: RascalOptions,
}

macro_rules! option_props {
    ($ty:ident, $($py:literal $get:ident $set:ident: $t:ty => $field:ident, $doc:literal;)+) => {
        #[pymethods]
        impl $ty {
            $(
                #[doc = $doc]
                #[getter($py)]
                fn $get(&self) -> $t {
                    self.opts.$field.clone()
                }
                #[setter($py)]
                fn $set(&mut self, v: $t) {
                    self.opts.$field = v;
                }
            )+
        }
    };
}

#[pymethods]
impl PyRascalOptions {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }
}

option_props! {
    PyRascalOptions,
    "similarityThreshold" similarity_threshold set_similarity_threshold: f64 => similarity_threshold,
        "Threshold below which MCES won't be run.  Between 0.0 and 1.0, default=0.7.";
    "singleLargestFrag" single_largest_frag set_single_largest_frag: bool => single_largest_frag,
        "Return the just single largest fragment of the MCES. It is equivalent to running with allBestMCEs=True, finding the result with the largest largestFragmentSize, and calling its largestFragmentOnly method.  This option may not produce the largest possible single fragment that the molecules have in common. If you definitely want that you may be better off using rdFMCS.";
    "completeAromaticRings" complete_aromatic_rings set_complete_aromatic_rings: bool => complete_aromatic_rings,
        "If True (default), partial aromatic rings won't be returned.";
    "ringMatchesRingOnly" ring_matches_ring_only set_ring_matches_ring_only: bool => ring_matches_ring_only,
        "If True (default is False), ring bonds won't match non-ring bonds.";
    "completeSmallestRings" complete_smallest_rings set_complete_smallest_rings: bool => complete_smallest_rings,
        "If True (default is False), only complete rings present in both input molecule's RingInfo will be returned. Implies completeAromaticRings and ringMatchesRingOnly.";
    "exactConnectionsMatch" exact_connections_match set_exact_connections_match: bool => exact_connections_match,
        "If True (default is False), atoms will only match atoms if they have the same number of explicit connections.  E.g. the central atom of C(C)(C) won't match either atom in CC";
    "minFragSize" min_frag_size set_min_frag_size: i32 => min_frag_size,
        "Imposes a minimum on the number of atoms in a fragment that may be part of the MCES.  Default -1 means no minimum.";
    "maxFragSeparation" max_frag_separation set_max_frag_separation: i32 => max_frag_separation,
        "Maximum number of bonds between fragments in the MCES for both to be reported.  Default -1 means no maximum.  If exceeded, the smaller fragment will be removed.";
    "allBestMCESs" all_best_mcess set_all_best_mcess: bool => all_best_mcess,
        "If True, reports all MCESs found of the same maximum size.  Default False means just report the first found.";
    "maxBestMCESs" max_best_mcess set_max_best_mcess: u32 => max_best_mcess,
        "Some pathological cases produce huge numbers of equivalent solutions that can crash the program due to memory depletion.  This caps the number of such solutions to prevent this happening.  Default=10000.";
    "returnEmptyMCES" return_empty_mces set_return_empty_mces: bool => return_empty_mces,
        "If the estimated similarity between the 2 molecules doesn't meet the similarityThreshold, no results are returned.  If you want to know what the estimates were, set this to True, and examine the tier1Sim and tier2Sim properties of the result then returned.";
    "timeout" timeout set_timeout: i32 => timeout,
        "Maximum time (in seconds) to spend on an individual MCESs determination.  Default 60, -1 means no limit.";
    "maxBondMatchPairs" max_bond_match_pairs set_max_bond_match_pairs: u32 => max_bond_match_pairs,
        "Too many matching bond (vertex) pairs can cause the process to run out of memory.  The default of 1000 is fairly safe.  Increase with caution, as memory use increases with the square of this number.  ";
    "equivalentAtoms" equivalent_atoms set_equivalent_atoms: String => equivalent_atoms,
        "SMARTS strings defining atoms that should be considered equivalent. e.g. [F,Cl,Br,I] so all halogens will match each other. Space-separated list allowing more than 1 class of equivalent atoms.";
    "ignoreBondOrders" ignore_bond_orders set_ignore_bond_orders: bool => ignore_bond_orders,
        "If True, will treat all bonds as the same, irrespective of order.  Default=False.";
    "ignoreAtomAromaticity" ignore_atom_aromaticity set_ignore_atom_aromaticity: bool => ignore_atom_aromaticity,
        "If True, matches atoms solely on atomic number.  If False, will treat aromatic and aliphatic atoms as different.  Default=True.";
    "minCliqueSize" min_clique_size set_min_clique_size: u32 => min_clique_size,
        "Normally, the minimum clique size is specified via the similarityThreshold.  Sometimes it's more convenient to specify it directly.  If this is > 0, it will over-ride the similarityThreshold.  Note that this refers to the minimum number of BONDS in the MCES. Default=0.";
}

/// RASCAL Cluster Options.  Most of these pertain to RascalCluster
/// calculations.  Only similarityCutoff is used by RascalButinaCluster.
#[pyclass(
    name = "RascalClusterOptions",
    module = "rdkit.Chem.rdRascalMCES",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyRascalClusterOptions {
    opts: RascalClusterOptions,
}

#[pymethods]
impl PyRascalClusterOptions {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// The weight of matched bonds over matched atoms. Default=2.
    ///
    /// (As in upstream, whose binding of `b` is the member `a`, this reads and
    /// writes the penalty `a`.)
    #[getter(b)]
    #[allow(clippy::misnamed_getters)]
    fn b(&self) -> f64 {
        self.opts.a
    }
    #[setter(b)]
    fn set_b(&mut self, v: f64) {
        self.opts.a = v;
    }
}

option_props! {
    PyRascalClusterOptions,
    "similarityCutoff" similarity_cutoff set_similarity_cutoff: f64 => similarity_cutoff,
        "Similarity cutoff for molecules to be in the same cluster.  Between 0.0 and 1.0, default=0.7.";
    "minFragSize" min_frag_size set_min_frag_size: u32 => min_frag_size,
        "The minimum number of atoms in a fragment for it to be included in the MCES.  Default=3.";
    "maxNumFrags" max_num_frags set_max_num_frags: u32 => max_num_frags,
        "The maximum number of fragments allowed in the MCES for each pair of molecules. Default=2.  So that the MCES isn't a lot of small fragments scattered around the molecules giving an inflated estimate of similarity.";
    "numThreads" num_threads set_num_threads: i32 => num_threads,
        "Number of threads to use during clustering.  Default=-1 means all the hardware threads less one.";
    "a" a set_a: f64 => a,
        "The penalty score for each unconnected component in the MCES. Default=0.05.";
    "minIntraClusterSim" min_intra_cluster_sim set_min_intra_cluster_sim: f64 => min_intra_cluster_sim,
        "Two pairs of molecules are included in the same cluster if the similarity between their MCESs is greater than this.  Default=0.9.";
    "clusterMergeSim" cluster_merge_sim set_cluster_merge_sim: f64 => cluster_merge_sim,
        "Two clusters are merged if the fraction of molecules they have in common is greater than this.  Default=0.6.";
}

/// Find one or more MCESs between the 2 molecules given.  Returns a list of
/// RascalResult objects.
/// - mol1
/// - mol2 The two molecules for which to find the MCES
/// - opts Optional RascalOptions object changing the default run mode.
#[pyfunction]
#[pyo3(name = "FindMCES", signature = (mol1, mol2, opts = None))]
fn py_find_mces<'py>(
    py: Python<'py>,
    mol1: PyRef<'_, PyMol>,
    mol2: PyRef<'_, PyMol>,
    opts: Option<PyRef<'_, PyRascalOptions>>,
) -> PyResult<Bound<'py, PyList>> {
    let opts = opts.map(|o| o.opts.clone()).unwrap_or_default();
    let (m1, m2) = (&mol1.mol, &mol2.mol);
    let results = py.detach(|| rascal_mces(m1, m2, &opts)).py()?;
    PyList::new(py, results.into_iter().map(|res| PyRascalResult { res }))
}

fn extract_mols(mols: &Bound<'_, PyAny>) -> PyResult<Vec<Arc<ROMol>>> {
    let mut out = Vec::new();
    for m in mols.try_iter()? {
        let m = m?;
        if m.is_none() {
            return Err(PyValueError::new_err("molecule is None"));
        }
        out.push(Arc::new(m.extract::<PyRef<'_, PyMol>>()?.mol.clone()));
    }
    Ok(out)
}

fn cluster_opts(opts: Option<PyRef<'_, PyRascalClusterOptions>>) -> RascalClusterOptions {
    opts.map(|o| o.opts.clone()).unwrap_or_default()
}

/// Use the RASCAL MCES similarity metric to do fuzzy clustering.  Returns a
/// list of lists of molecules, each inner list being a cluster.  The last
/// cluster is all the molecules that didn't fit into another cluster (the
/// singletons).
/// - mols List of molecules to be clustered
/// - opts Optional RascalOptions object changing the default run mode.
#[pyfunction]
#[pyo3(name = "RascalCluster", signature = (mols, opts = None))]
fn py_rascal_cluster(
    py: Python<'_>,
    mols: &Bound<'_, PyAny>,
    opts: Option<PyRef<'_, PyRascalClusterOptions>>,
) -> PyResult<Vec<Vec<u32>>> {
    let (mols, opts) = (extract_mols(mols)?, cluster_opts(opts));
    py.detach(|| rascal_cluster(&mols, &opts)).py()
}

/// Use the RASCAL MCES similarity metric to do Butina clustering (Butina
/// JCICS 39 747-750 (1999)).  Returns a list of lists of molecules, each inner
/// list being a cluster.  The last cluster is all the molecules that didn't
/// fit into another cluster (the singletons).
/// - mols List of molecules to be clustered
/// - opts Optional RascalOptions object changing the default run mode.
#[pyfunction]
#[pyo3(name = "RascalButinaCluster", signature = (mols, opts = None))]
fn py_rascal_butina_cluster(
    py: Python<'_>,
    mols: &Bound<'_, PyAny>,
    opts: Option<PyRef<'_, PyRascalClusterOptions>>,
) -> PyResult<Vec<Vec<u32>>> {
    let (mols, opts) = (extract_mols(mols)?, cluster_opts(opts));
    py.detach(|| rascal_butina_cluster(&mols, &opts)).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing implementation of RASCAL Maximum Common Edge Substructure algorithm.",
    )?;
    m.add_class::<PyRascalResult>()?;
    m.add_class::<PyRascalOptions>()?;
    m.add_class::<PyRascalClusterOptions>()?;
    m.add_function(wrap_pyfunction!(py_find_mces, m)?)?;
    m.add_function(wrap_pyfunction!(py_rascal_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(py_rascal_butina_cluster, m)?)?;
    Ok(())
}
