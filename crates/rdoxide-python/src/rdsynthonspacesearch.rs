//! `rdkit.Chem.rdSynthonSpaceSearch` (upstream
//! `SynthonSpaceSearch/Wrap/rdSynthonSpaceSearch.cpp`): searching of
//! synthon-based chemical libraries.

use std::sync::{Arc, Mutex};

use pyo3::prelude::*;
use pyo3::types::PyList;
use rdoxide::RDKitError;
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::substruct::SubstructMatchParameters;
use rdoxide::graphmol::synthon_space_search::searcher::SearchResultCallback;
use rdoxide::graphmol::synthon_space_search::{
    SearchResults, ShapeBuildParams, SynthonSpace, SynthonSpaceSearchParams, UserConfGenerator,
    convert_text_to_db_file, formatted_integer_string,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdenumeratestereoisomers::PyStereoEnumerationOptions;
use crate::rdfingerprintgenerator::PyFPGen;
use crate::rdgaussianshape::{PyShapeInput, PyShapeOverlayOptions};
use crate::rdgeneralizedsubstruct::PyExtendedQueryMol;
use crate::rdrascalmces::PyRascalOptions;
use crate::substruct::PySubstructMatchParameters;

/// Used to return results of SynthonSpace searches.
#[pyclass(
    name = "SubstructureResult",
    module = "rdkit.Chem.rdSynthonSpaceSearch"
)]
pub struct PySearchResults {
    res: SearchResults,
}

#[pymethods]
impl PySearchResults {
    /// A function returning hits from the search
    #[pyo3(name = "GetHitMolecules")]
    fn get_hit_molecules<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty(py);
        for m in self.res.hit_molecules() {
            list.append(Py::new(py, PyMol { mol: m.clone() })?)?;
        }
        Ok(list)
    }

    /// The upper bound on number of results possible.  There may be fewer
    /// than this in practice for several reasons such as duplicate reagent
    /// sets being removed or the final product not matching the query even
    /// though the synthons suggested they would.
    #[pyo3(name = "GetMaxNumResults")]
    fn get_max_num_results(&self) -> u64 {
        self.res.max_num_results()
    }

    /// Returns whether the search timed out or not.
    #[pyo3(name = "GetTimedOut")]
    fn get_timed_out(&self) -> bool {
        self.res.timed_out()
    }

    /// Returns whether the search was cancelled or not.
    #[pyo3(name = "GetCancelled")]
    fn get_cancelled(&self) -> bool {
        self.res.cancelled()
    }

    /// Returns the best hit found in the similarity search, even when none
    /// were under the search threshold.  May be empty if not a similarity
    /// search or nothing came close to a match.
    #[pyo3(name = "GetBestHit")]
    fn get_best_hit(&self) -> Option<PyMol> {
        self.res.best_hit().map(|m| PyMol { mol: m.clone() })
    }
}

fn user_conf_gen(func: Py<PyAny>) -> UserConfGenerator {
    Arc::new(move |smiles: &str, num_confs: u32| {
        Python::attach(|py| {
            let res = func.call1(py, (smiles, num_confs)).ok()?;
            let m = res.extract::<PyRef<'_, PyMol>>(py).ok()?;
            Some(m.mol.clone())
        })
    })
}

macro_rules! simple_props {
    ($ty:ident, $($py:literal $get:ident $set:ident: $t:ty => $field:ident, $doc:literal;)+) => {
        #[pymethods]
        impl $ty {
            $(
                #[doc = $doc]
                #[getter($py)]
                fn $get(&self) -> $t {
                    self.p.$field.clone()
                }
                #[setter($py)]
                fn $set(&mut self, v: $t) {
                    self.p.$field = v;
                }
            )+
        }
    };
}

/// SynthonSpaceSearch parameters.
#[pyclass(
    name = "SynthonSpaceSearchParams",
    module = "rdkit.Chem.rdSynthonSpaceSearch",
    skip_from_py_object
)]
#[derive(Clone, Default)]
pub struct PySynthonSpaceSearchParams {
    pub(crate) p: SynthonSpaceSearchParams,
}

#[pymethods]
impl PySynthonSpaceSearchParams {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// Options for the shape overlays.
    #[getter(shapeOverlayOptions)]
    fn shape_overlay_options(&self) -> PyShapeOverlayOptions {
        PyShapeOverlayOptions {
            opts: self.p.shape_overlay_options.clone(),
        }
    }
    #[setter(shapeOverlayOptions)]
    fn set_shape_overlay_options(&mut self, v: PyRef<'_, PyShapeOverlayOptions>) {
        self.p.shape_overlay_options = v.opts.clone();
    }

    /// Options for stereoisomer enumeration.
    #[getter(stereoEnumOpts)]
    fn stereo_enum_opts(&self) -> PyStereoEnumerationOptions {
        PyStereoEnumerationOptions {
            opts: self.p.stereo_enum_opts.clone(),
        }
    }
    #[setter(stereoEnumOpts)]
    fn set_stereo_enum_opts(&mut self, v: PyRef<'_, PyStereoEnumerationOptions>) {
        self.p.stereo_enum_opts = v.opts.clone();
    }

    /// Add an excluded volume to use in the shape search.  The volume
    /// overlap and mean overlap over clashing atoms will be reported.
    #[getter(excludedVolume)]
    fn excluded_volume(&self) -> Option<PyShapeInput> {
        self.p.excluded_volume.as_ref().map(|s| PyShapeInput {
            shape: (**s).clone(),
        })
    }
    #[setter(excludedVolume)]
    fn set_excluded_volume(&mut self, v: PyRef<'_, PyShapeInput>) {
        self.p.excluded_volume = Some(Arc::new(v.shape.clone()));
    }

    /// Allows you to provide a function that will be called instead of the
    /// default conformer generator to generate conformers for the synthons.
    /// The function should take a SMILES string and the maximum number of
    /// conformers to generated and return a molecule object.
    #[pyo3(name = "setUserConformerGenerator")]
    fn set_user_conformer_generator(&mut self, func: Py<PyAny>) {
        self.p.user_conformer_generator = Some(user_conf_gen(func));
    }
}

simple_props! {
    PySynthonSpaceSearchParams,
    "maxHits" max_hits set_max_hits: i64 => max_hits,
        "The maximum number of hits to return.  Default=1000.Use -1 for no maximum.";
    "maxNumFrags" max_num_frags set_max_num_frags: u64 => max_num_frag_sets,
        "The maximum number of fragments the query can be broken into.  Big molecules will create huge numbers of fragments that may cause excessive memory use.  If the number of fragments hits this number, fragmentation stops and the search results will likely be incomplete.  Default=100000.";
    "hitStart" hit_start set_hit_start: i64 => hit_start,
        "The sequence number of the hit to start from.  So that you can return the next N hits of a search having already obtained N-1.  Default=0";
    "randomSample" random_sample set_random_sample: bool => random_sample,
        "If True, returns a random sample of the hits, up to maxHits in number.  Default=False.";
    "randomSeed" random_seed set_random_seed: i32 => random_seed,
        "If using randomSample, this seeds the random number generator so as to give reproducible results.  Default=-1 means use a random seed.";
    "buildHits" build_hits set_build_hits: bool => build_hits,
        "If false, reports the maximum number of hits that the search could produce, but doesn't return them.";
    "similarityCutoff" similarity_cutoff set_similarity_cutoff: f64 => similarity_cutoff,
        "Similarity cutoff for returning hits by fingerprint similarity.  At present the fp is hard-coded to be Morgan, bits, radius=2.  Default=0.5.";
    "fragSimilarityAdjuster" frag_similarity_adjuster set_frag_similarity_adjuster: f64 => frag_similarity_adjuster,
        "Similarities of fragments are generally low due to low bit densities.  For the fragment matching, reduce the similarity cutoff off by this amount.  Default=0.1.";
    "approxSimilarityAdjuster" approx_similarity_adjuster set_approx_similarity_adjuster: f64 => approx_similarity_adjuster,
        "The fingerprint search uses an approximate similarity method before building a product and doing a final check.  The similarityCutoff is reduced by this value for the approximate check.  A lower value will give faster run times at the risk of missing some hits.  The value you use should have a positive correlation with your FOMO.  The default of 0.1 is appropriate for Morgan fingerprints.  With RDKit fingerprints, 0.05 is adequate, and higher than that has been seen to produce long run times.";
    "minHitHeavyAtoms" min_hit_heavy_atoms set_min_hit_heavy_atoms: u32 => min_hit_heavy_atoms,
        "Minimum number of heavy atoms in a hit.  Default=0.";
    "maxHitHeavyAtoms" max_hit_heavy_atoms set_max_hit_heavy_atoms: i32 => max_hit_heavy_atoms,
        "Maximum number of heavy atoms in a hit.  Default=-1 means no maximum.";
    "minHitMolWt" min_hit_mol_wt set_min_hit_mol_wt: f64 => min_hit_mol_wt,
        "Minimum molecular weight for a hit.  Default=0.0.";
    "maxHitMolWt" max_hit_mol_wt set_max_hit_mol_wt: f64 => max_hit_mol_wt,
        "Maximum molecular weight for a hit.  Default=0.0 mean no maximum.";
    "minHitChiralAtoms" min_hit_chiral_atoms set_min_hit_chiral_atoms: u32 => min_hit_chiral_atoms,
        "Minimum number of chiral atoms in a hit.  Default=0.";
    "maxHitChiralAtoms" max_hit_chiral_atoms set_max_hit_chiral_atoms: i32 => max_hit_chiral_atoms,
        "Maximum number of chiral atoms in a hit.  Default=-1 means no maximum.";
    "numConformers" num_conformers set_num_conformers: u32 => num_conformers,
        "When doing a shape search, the number of conformers to generate for molecules.  Default=100.";
    "confRMSThreshold" conf_rms_threshold set_conf_rms_threshold: f64 => conf_rms_threshold,
        "When doing a shape search, the RMS threshold to use when pruning conformers.  Default=1.0.";
    "bestHit" best_hit set_best_hit: bool => best_hit,
        "If True, when doing a shape search it will return the hit conformer with the best shape match to the query conformer.  If False, it just returns the first hit conformer that exceeds the similarity cutoff.  The latter will be faster but the returned hit conformations are likely to be less relevant.";
    "enumerateUnspecifiedStereo" enumerate_unspecified_stereo set_enumerate_unspecified_stereo: bool => enumerate_unspecified_stereo,
        "When doing a shape search, if there is unspecified stereochemistry in either the query or potential hit, enumerate test all possibilities.  Default=False.";
    "timeOut" time_out set_time_out: u64 => time_out,
        "Time limit for search, in seconds.  Default is 600s, 0 means no timeout.  Requires an integer";
    "toTryChunkSize" to_try_chunk_size set_to_try_chunk_size: u64 => to_try_chunk_size,
        "Process possible hits using the given chunk size";
    "numThreads" num_threads set_num_threads: i32 => num_threads,
        "The number of threads to use for search.  If > 0, will use that number.  If <= 0, will use the number of hardware threads plus this number.  So if the number of hardware threads is 8, and numThreads is -1, it will use 7 threads.  Default=1.";
    "useProgressBar" use_progress_bar set_use_progress_bar: u32 => use_progress_bar,
        "Makes a progress bar of given width.  The number given is the number of '*' characters in a full bar.  There will be about another 35 characters or so depending on the size of the job.  Default=0 means no bar.";
    "maxExcludedVolume" max_excluded_volume set_max_excluded_volume: f64 => max_excluded_volume,
        "Maximum allowed excluded volume for a hit to be accepted.  Default -1.0 means no maximum.";
    "maxMeanExcludedVolume" max_mean_excluded_volume set_max_mean_excluded_volume: f64 => max_mean_excluded_volume,
        "Maximum mean excluded volume for a hit to be accepted.  The mean is the total excluded volume divided by the number of clashing atoms (within 2 CARBON_RAD of an excluded volume atom).  To try and distinguish between a mild clash over the whole hit and a few atoms having a really bad clash.";
    "possibleHitsFile" possible_hits_file set_possible_hits_file: String => possible_hits_file,
        "Name of a file to save the possible hits to. These are the combinations of synthons that might match the query but need building and final checking.  Each line has a space-separated list of the synthons and the hit's name.  The file will be emptied and re-filled if it already exists.";
    "maxPossibleHitsToWrite" max_possible_hits_to_write set_max_possible_hits_to_write: u64 => max_possible_hits_to_write,
        "Maximum number of lines to write to possibleHitsFile.  When dealing with huge synthon spaces it's very easy to fill a disk.  Default=10M.";
    "writePossibleHitsAndStop" write_possible_hits_and_stop set_write_possible_hits_and_stop: bool => write_possible_hits_and_stop,
        "If True, creates the possibleHitsFile and stops without doing the final building and checking.  Default is False.";
}

/// Parameters for building shape objects for SynthonSpaceSearch.
#[pyclass(
    name = "ShapeBuildParams",
    module = "rdkit.Chem.rdSynthonSpaceSearch",
    skip_from_py_object
)]
#[derive(Clone, Default)]
pub struct PyShapeBuildParams {
    pub(crate) p: ShapeBuildParams,
}

#[pymethods]
impl PyShapeBuildParams {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// Options for stereoisomer enumeration.
    #[getter(stereoEnumOpts)]
    fn stereo_enum_opts(&self) -> PyStereoEnumerationOptions {
        PyStereoEnumerationOptions {
            opts: self.p.stereo_enum_opts.clone(),
        }
    }
    #[setter(stereoEnumOpts)]
    fn set_stereo_enum_opts(&mut self, v: PyRef<'_, PyStereoEnumerationOptions>) {
        self.p.stereo_enum_opts = v.opts.clone();
    }

    /// Allows you to provide a function that will be called instead of the
    /// default conformer generator to generate conformers for the synthons.
    /// The function should take a SMILES string and the maximum number of
    /// conformers to generated and return a molecule object.
    #[pyo3(name = "setUserConformerGenerator")]
    fn set_user_conformer_generator(&mut self, func: Py<PyAny>) {
        self.p.user_conformer_generator = Some(user_conf_gen(func));
    }
}

simple_props! {
    PyShapeBuildParams,
    "numConfs" num_confs set_num_confs: u32 => num_confs,
        "Maximum number of conformers per synthon or query.  Default=10";
    "rmsThreshold" rms_threshold set_rms_threshold: f64 => rms_threshold,
        "RMS threshold to use when pruning conformations.  Default=1.0.";
    "shapeSimThreshold" shape_sim_threshold set_shape_sim_threshold: f64 => shape_sim_threshold,
        "When generating shapes, similarity threshold for pruning.  No 2 shapes for each synthon or query will be more similar than this threshold.  Default=1.9.";
    "numThreads" num_threads set_num_threads: i32 => num_threads,
        "The number of threads to use for shape building.  If > 0, will use that number.  If <= 0, will use the number of hardware threads plus this number.Default=1.";
    "randomSeed" random_seed set_random_seed: i32 => random_seed,
        "Seed for random number generator.  Default=-1 means use system random seed.";
    "useProgressBar" use_progress_bar set_use_progress_bar: u32 => use_progress_bar,
        "Makes a progress bar of given width.  The number given is the number of '*' characters in a full bar.  There will be about another 35 characters or so depending on the size of the job.  Default=0 means no bar.";
    "maxSynthonAtoms" max_synthon_atoms set_max_synthon_atoms: u32 => max_synthon_atoms,
        "If >0, sets a maximum number of heavy atoms, excluding dummies, for synthon to have a shape made.  Default=0.";
    "maxEmbedAttempts" max_embed_attempts set_max_embed_attempts: u32 => max_embed_attempts,
        "Maximum number of attempts for embedding a single synthon.  Default=10.";
    "timeOut" time_out set_time_out: u32 => time_out,
        "Maximum time in seconds to spend on each synthon when generating conformers.  Default=600 means no timeout.";
    "interimFile" interim_file set_interim_file: String => interim_file,
        "Interim file to write the SynthonSpace to during shape generation.  In the event of a failure, a restart from this file will be possible.";
    "interimWrites" interim_writes set_interim_writes: u64 => interim_writes,
        "If an interim file has been given, every this many shapes write a new version of the file.  Default=1000.";
}

/// SynthonSpaceSearch object.
#[pyclass(name = "SynthonSpace", module = "rdkit.Chem.rdSynthonSpaceSearch")]
#[derive(Default)]
pub struct PySynthonSpace {
    space: SynthonSpace,
}

fn match_params(sm: Option<PyRef<'_, PySubstructMatchParameters>>) -> SubstructMatchParameters {
    sm.map(|p| p.base_params()).unwrap_or_default()
}

fn search_params(p: Option<PyRef<'_, PySynthonSpaceSearchParams>>) -> SynthonSpaceSearchParams {
    p.map(|p| p.p.clone()).unwrap_or_default()
}

fn finish(res: rdoxide::Result<SearchResults>, what: &str) -> PyResult<PySearchResults> {
    let res = res.py()?;
    if res.cancelled() {
        return Err(pyo3::exceptions::PyKeyboardInterrupt::new_err(format!(
            "{what} cancelled"
        )));
    }
    Ok(PySearchResults { res })
}

#[pymethods]
impl PySynthonSpace {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// Reads text file of the sort used by ChemSpace/Enamine.
    #[pyo3(name = "ReadTextFile", signature = (inFile))]
    #[allow(non_snake_case)]
    fn read_text_file(&mut self, py: Python<'_>, inFile: std::path::PathBuf) -> PyResult<()> {
        let f = inFile.to_string_lossy().into_owned();
        py.detach(|| self.space.read_text_file(&f)).py()
    }

    /// Reads binary database file.  Takes optional number of threads,default=1.
    #[pyo3(name = "ReadDBFile", signature = (inFile, numThreads = 1))]
    #[allow(non_snake_case)]
    fn read_db_file(
        &mut self,
        py: Python<'_>,
        inFile: std::path::PathBuf,
        numThreads: i32,
    ) -> PyResult<()> {
        let f = inFile.to_string_lossy().into_owned();
        py.detach(|| self.space.read_db_file(&f, numThreads)).py()
    }

    /// Writes binary database file.
    #[pyo3(name = "WriteDBFile", signature = (outFile))]
    #[allow(non_snake_case)]
    fn write_db_file(&self, py: Python<'_>, outFile: std::path::PathBuf) -> PyResult<()> {
        let f = outFile.to_string_lossy().into_owned();
        py.detach(|| self.space.write_db_file(&f)).py()
    }

    /// Writes enumerated library to file.
    #[pyo3(name = "WriteEnumeratedFile", signature = (outFile))]
    #[allow(non_snake_case)]
    fn write_enumerated_file(&self, py: Python<'_>, outFile: std::path::PathBuf) -> PyResult<()> {
        let f = outFile.to_string_lossy().into_owned();
        py.detach(|| self.space.write_enumerated_file(&f)).py()
    }

    /// Returns number of reactions in the SynthonSpace.
    #[pyo3(name = "GetNumReactions")]
    fn get_num_reactions(&self) -> usize {
        self.space.num_reactions()
    }

    /// Returns number of products in the SynthonSpace, with multiple
    /// counting of any duplicates.
    #[pyo3(name = "GetNumProducts")]
    fn get_num_products(&self) -> u64 {
        self.space.num_products()
    }

    /// Returns number of synthons in the SynthonSpace.
    #[pyo3(name = "GetNumSynthons")]
    fn get_num_synthons(&self) -> usize {
        self.space.num_synthons()
    }

    /// Returns the number of synthons in the SynthonSpace that have a shape.
    #[pyo3(name = "GetNumSynthonsWithShapes")]
    fn get_num_synthons_with_shapes(&self) -> u64 {
        self.space.num_synthons_with_shapes()
    }

    /// Writes a summary of the SynthonSpace to stdout.
    #[pyo3(name = "Summarise")]
    fn summarise(&self, py: Python<'_>) -> PyResult<()> {
        let text = self.space.summarise();
        py.import("sys")?
            .getattr("stdout")?
            .call_method1("write", (text,))?;
        Ok(())
    }

    /// Returns the information string for the fingerprint generator used to
    /// create this space.
    #[pyo3(name = "GetSynthonFingerprintType")]
    fn get_synthon_fingerprint_type(&self) -> String {
        self.space.fp_type().to_string()
    }

    /// Does a substructure search in the SynthonSpace.  The query may be a
    /// molecule or an ExtendedQueryMol.  With startLine and finishLine the
    /// contents of params.possibleHitsFile, which is assumed to have been
    /// written by an earlier search, between those lines are extracted that
    /// are indeed hits.
    #[pyo3(
        name = "SubstructureSearch",
        signature = (query, substructMatchParams = None, params = None, startLine = None, finishLine = None)
    )]
    #[allow(non_snake_case)]
    fn substructure_search(
        &self,
        py: Python<'_>,
        query: &Bound<'_, PyAny>,
        substructMatchParams: Option<PyRef<'_, PySubstructMatchParameters>>,
        params: Option<PyRef<'_, PySynthonSpaceSearchParams>>,
        startLine: Option<u64>,
        finishLine: Option<u64>,
    ) -> PyResult<PySearchResults> {
        let sm = match_params(substructMatchParams);
        let ps = search_params(params);
        let lines = startLine.zip(finishLine);
        if let Ok(m) = query.extract::<PyRef<'_, PyMol>>() {
            let mol: &ROMol = &m.mol;
            let res = py.detach(|| match lines {
                None => self.space.substructure_search(mol, &sm, &ps),
                Some((s, f)) => self.space.substructure_search_lines(mol, &sm, &ps, s, f),
            });
            return finish(res, "SubstructureSearch");
        }
        let xq = query.extract::<PyRef<'_, PyExtendedQueryMol>>()?;
        let q = &xq.query;
        let res = py.detach(|| match lines {
            None => self.space.substructure_search_xq(q, &sm, &ps),
            Some((s, f)) => self.space.substructure_search_xq_lines(q, &sm, &ps, s, f),
        });
        finish(res, "SubstructureSearch")
    }

    /// Does a substructure search in the SynthonSpace returning results in
    /// the callback.
    #[pyo3(
        name = "SubstructureSearchIncremental",
        signature = (query, callback, substructMatchParams = None, params = None)
    )]
    #[allow(non_snake_case)]
    fn substructure_search_incremental(
        &self,
        py: Python<'_>,
        query: PyRef<'_, PyMol>,
        callback: Py<PyAny>,
        substructMatchParams: Option<PyRef<'_, PySubstructMatchParameters>>,
        params: Option<PyRef<'_, PySynthonSpaceSearchParams>>,
    ) -> PyResult<()> {
        let sm = match_params(substructMatchParams);
        let ps = search_params(params);
        let mol: &ROMol = &query.mol;
        run_incremental(py, &callback, |cb| {
            self.space
                .substructure_search_incremental(mol, cb, &sm, &ps)
        })
    }

    /// Does a fingerprint search in the SynthonSpace using the
    /// FingerprintGenerator passed in.  With startLine and finishLine the
    /// contents of params.possibleHitsFile between those lines are
    /// extracted that are indeed hits.
    #[pyo3(
        name = "FingerprintSearch",
        signature = (query, fingerprintGenerator, params = None, startLine = None, finishLine = None)
    )]
    #[allow(non_snake_case)]
    fn fingerprint_search(
        &mut self,
        py: Python<'_>,
        query: PyRef<'_, PyMol>,
        fingerprintGenerator: PyRef<'_, PyFPGen>,
        params: Option<PyRef<'_, PySynthonSpaceSearchParams>>,
        startLine: Option<u64>,
        finishLine: Option<u64>,
    ) -> PyResult<PySearchResults> {
        let ps = search_params(params);
        let mol: &ROMol = &query.mol;
        let gen_ = &fingerprintGenerator.gen_;
        let space = &mut self.space;
        let res = py.detach(|| match startLine.zip(finishLine) {
            None => space.fingerprint_search(mol, gen_, &ps),
            Some((s, f)) => space.fingerprint_search_lines(mol, gen_, &ps, s, f),
        });
        finish(res, "FingerprintSearch")
    }

    /// Does a fingerprint search in the SynthonSpace using the
    /// FingerprintGenerator passed in, returning results the callback.
    #[pyo3(
        name = "FingerprintSearchIncremental",
        signature = (query, fingerprintGenerator, callback, params = None)
    )]
    #[allow(non_snake_case)]
    fn fingerprint_search_incremental(
        &mut self,
        py: Python<'_>,
        query: PyRef<'_, PyMol>,
        fingerprintGenerator: PyRef<'_, PyFPGen>,
        callback: Py<PyAny>,
        params: Option<PyRef<'_, PySynthonSpaceSearchParams>>,
    ) -> PyResult<()> {
        let ps = search_params(params);
        let mol: &ROMol = &query.mol;
        let gen_ = &fingerprintGenerator.gen_;
        let space = &mut self.space;
        run_incremental(py, &callback, |cb| {
            space.fingerprint_search_incremental(mol, gen_, cb, &ps)
        })
    }

    /// Does a search using the Rascal similarity score.  The similarity
    /// threshold used is provided by rascalOptions, and the one in params is
    /// ignored.
    #[pyo3(
        name = "RascalSearch",
        signature = (query, rascalOptions, params = None, startLine = None, finishLine = None)
    )]
    #[allow(non_snake_case)]
    fn rascal_search(
        &self,
        py: Python<'_>,
        query: PyRef<'_, PyMol>,
        rascalOptions: PyRef<'_, PyRascalOptions>,
        params: Option<PyRef<'_, PySynthonSpaceSearchParams>>,
        startLine: Option<u64>,
        finishLine: Option<u64>,
    ) -> PyResult<PySearchResults> {
        let ps = search_params(params);
        let mol: &ROMol = &query.mol;
        let opts = &rascalOptions.opts;
        let res = py.detach(|| match startLine.zip(finishLine) {
            None => self.space.rascal_search(mol, opts, &ps),
            Some((s, f)) => self.space.rascal_search_lines(mol, opts, &ps, s, f),
        });
        finish(res, "RascalSearch")
    }

    /// Does a search using the Rascal similarity score, returning results
    /// iteratively in the callback.
    #[pyo3(
        name = "RascalSearchIncremental",
        signature = (query, rascalOptions, callback, params = None)
    )]
    #[allow(non_snake_case)]
    fn rascal_search_incremental(
        &self,
        py: Python<'_>,
        query: PyRef<'_, PyMol>,
        rascalOptions: PyRef<'_, PyRascalOptions>,
        callback: Py<PyAny>,
        params: Option<PyRef<'_, PySynthonSpaceSearchParams>>,
    ) -> PyResult<()> {
        let ps = search_params(params);
        let mol: &ROMol = &query.mol;
        let opts = &rascalOptions.opts;
        run_incremental(py, &callback, |cb| {
            self.space.rascal_search_incremental(mol, opts, cb, &ps)
        })
    }

    /// Perform a shape similarity search with the given query molecule
    /// across the synthonspace library.  Duplicate SMILES strings produced
    /// by different reactions will be returned.  Requires a query with at
    /// least 1 3D conformer.  Only the first conformer will be used in the
    /// search.  With startLine and finishLine the contents of
    /// params.possibleHitsFile, which is assumed to have been written by an
    /// earlier search, between those lines are extracted that are indeed
    /// hits.
    #[pyo3(
        name = "ShapeSearch",
        signature = (query, params = None, startLine = None, finishLine = None)
    )]
    #[allow(non_snake_case)]
    fn shape_search(
        &self,
        py: Python<'_>,
        query: PyRef<'_, PyMol>,
        params: Option<PyRef<'_, PySynthonSpaceSearchParams>>,
        startLine: Option<u64>,
        finishLine: Option<u64>,
    ) -> PyResult<PySearchResults> {
        let ps = search_params(params);
        let mol: &ROMol = &query.mol;
        let res = py.detach(|| match startLine.zip(finishLine) {
            None => self.space.shape_search(mol, &ps),
            Some((s, f)) => self.space.shape_search_lines(mol, &ps, s, f),
        });
        finish(res, "ShapeSearch")
    }

    /// Build shapes for the synthons.  The conformations are generated,
    /// pruned with the given threshold, which is passed directly to
    /// EmbedMultipleConfs.
    #[pyo3(name = "BuildSynthonShapes", signature = (params = None))]
    fn build_synthon_shapes(
        &mut self,
        py: Python<'_>,
        params: Option<PyRef<'_, PyShapeBuildParams>>,
    ) -> PyResult<()> {
        let ps = params.map(|p| p.p.clone()).unwrap_or_default();
        let space = &mut self.space;
        py.detach(|| space.build_synthon_shapes(&ps)).py()
    }

    /// Build the synthon fingerprints ready for similarity searching.  This
    /// is done automatically when the first similarity search is done.
    #[pyo3(
        name = "BuildSynthonFingerprints",
        signature = (fingerprintGenerator, progressBarWidth = 0)
    )]
    #[allow(non_snake_case)]
    fn build_synthon_fingerprints(
        &mut self,
        py: Python<'_>,
        fingerprintGenerator: PyRef<'_, PyFPGen>,
        progressBarWidth: u32,
    ) -> PyResult<()> {
        let _ = progressBarWidth;
        let gen_ = &fingerprintGenerator.gen_;
        let space = &mut self.space;
        py.detach(|| space.build_synthon_fingerprints(gen_)).py()
    }
}

/// Runs a search that reports to a Python callback in chunks; a truthy
/// return value from the callback stops the search and an exception is
/// propagated.
fn run_incremental(
    py: Python<'_>,
    callback: &Py<PyAny>,
    search: impl FnOnce(SearchResultCallback<'_>) -> rdoxide::Result<()> + Send,
) -> PyResult<()> {
    let py_err: Mutex<Option<PyErr>> = Mutex::new(None);
    let res = py.detach(|| {
        let mut cb = |mols: &mut Vec<ROMol>| -> rdoxide::Result<bool> {
            Python::attach(|py| {
                let mut call = || -> PyResult<bool> {
                    let list = PyList::empty(py);
                    for m in mols.drain(..) {
                        list.append(Py::new(py, PyMol { mol: m })?)?;
                    }
                    callback.call1(py, (list,))?.is_truthy(py)
                };
                call().map_err(|e| {
                    *py_err.lock().unwrap() = Some(e);
                    RDKitError::Runtime("callback failed".into())
                })
            })
        };
        search(&mut cb)
    });
    if let Some(e) = py_err.into_inner().unwrap() {
        return Err(e);
    }
    res.py()
}

/// Convert the text file into the binary DB file in our format.
#[pyfunction]
#[pyo3(name = "ConvertTextToDBFile", signature = (inFilename, outFilename, fpGen = None, py_shapeParams = None))]
#[allow(non_snake_case)]
fn py_convert_text_to_db_file(
    py: Python<'_>,
    inFilename: &str,
    outFilename: &str,
    fpGen: Option<PyRef<'_, PyFPGen>>,
    py_shapeParams: Option<PyRef<'_, PyShapeBuildParams>>,
) -> PyResult<()> {
    let shape_params = py_shapeParams.map(|p| p.p.clone());
    let fp_gen = fpGen.as_ref().map(|g| &g.gen_);
    py.detach(|| convert_text_to_db_file(inFilename, outFilename, fp_gen, shape_params.as_ref()))
        .py()
}

/// Format an integer with spaces every 3 digits for ease of reading
#[pyfunction]
#[pyo3(name = "FormattedIntegerString", signature = (value))]
fn py_formatted_integer_string(value: i64) -> String {
    formatted_integer_string(value)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing implementation of SynthonSpace search of Synthon-based chemical libraries such as Enamine REAL.  NOTE: This functionality is experimental and the API and/or results may change in future releases.",
    )?;
    m.add_class::<PySearchResults>()?;
    m.add_class::<PySynthonSpaceSearchParams>()?;
    m.add_class::<PyShapeBuildParams>()?;
    m.add_class::<PySynthonSpace>()?;
    m.add_function(wrap_pyfunction!(py_convert_text_to_db_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_formatted_integer_string, m)?)?;
    Ok(())
}
