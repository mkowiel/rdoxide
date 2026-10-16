//! `rdoxide.Chem.rdDistGeom` (upstream `DistGeomHelpers/Wrap/rdDistGeom.cpp`):
//! embedding molecules in 3D with distance geometry and ETKDG.

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::distgeom::{BoundsMatrix, triangle_smooth_bounds};
use rdoxide::graphmol::dist_geom_helpers::{
    CoordMap, EmbedFF, EmbedParameters, InitialEmbeddingMode, TopolBoundsOptions,
    TopolBoundsOutput, embed_molecule, embed_multiple_confs, embed_parameters_to_json,
    init_bounds_mat, set_topol_bounds,
};
use rdoxide::graphmol::forcefield_helpers::crystal_ff::{
    CrystalFFDetails, TorsionAngleRef, TorsionOptions, get_experimental_torsions_with_bonds,
};

use crate::enums::{add_enum_exported, to_py};
use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;
use crate::pyconv::{float64_array, float64_values};
use crate::rdmolfiles_write::Args;

/// Parameters controlling embedding
#[pyclass(
    name = "EmbedParameters",
    module = "rdoxide.Chem.rdDistGeom",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyEmbedParameters {
    params: EmbedParameters,
}

macro_rules! simple_props {
    ($( ($py:ident, $get:ident, $set:ident, $field:ident, $ty:ty) ),* $(,)?) => {
        #[pymethods]
        #[allow(non_snake_case)]
        impl PyEmbedParameters {
            $(
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.params.$field
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.params.$field = value;
                }
            )*
        }
    };
}

simple_props!(
    (
        maxIterations,
        get_max_iterations,
        set_max_iterations,
        max_iterations,
        u32
    ),
    (
        numThreads,
        get_num_threads,
        set_num_threads,
        num_threads,
        i32
    ),
    (timeout, get_timeout, set_timeout, timeout, u32),
    (
        randomSeed,
        get_random_seed,
        set_random_seed,
        random_seed,
        i32
    ),
    (
        clearConfs,
        get_clear_confs,
        set_clear_confs,
        clear_confs,
        bool
    ),
    (
        useRandomCoords,
        get_use_random_coords,
        set_use_random_coords,
        use_random_coords,
        bool
    ),
    (
        boxSizeMult,
        get_box_size_mult,
        set_box_size_mult,
        box_size_mult,
        f64
    ),
    (
        randNegEig,
        get_rand_neg_eig,
        set_rand_neg_eig,
        rand_neg_eig,
        bool
    ),
    (
        numZeroFail,
        get_num_zero_fail,
        set_num_zero_fail,
        num_zero_fail,
        u32
    ),
    (
        optimizerForceTol,
        get_optimizer_force_tol,
        set_optimizer_force_tol,
        optimizer_force_tol,
        f64
    ),
    (
        basinThresh,
        get_basin_thresh,
        set_basin_thresh,
        basin_thresh,
        f64
    ),
    (
        ignoreSmoothingFailures,
        get_ignore_smoothing_failures,
        set_ignore_smoothing_failures,
        ignore_smoothing_failures,
        bool
    ),
    (
        enforceChirality,
        get_enforce_chirality,
        set_enforce_chirality,
        enforce_chirality,
        bool
    ),
    (
        useExpTorsionAnglePrefs,
        get_use_exp_torsion_angle_prefs,
        set_use_exp_torsion_angle_prefs,
        use_exp_torsion_angle_prefs,
        bool
    ),
    (
        useBasicKnowledge,
        get_use_basic_knowledge,
        set_use_basic_knowledge,
        use_basic_knowledge,
        bool
    ),
    (ETversion, get_et_version, set_et_version, et_version, u32),
    (verbose, get_verbose, set_verbose, verbose, bool),
    (
        pruneRmsThresh,
        get_prune_rms_thresh,
        set_prune_rms_thresh,
        prune_rms_thresh,
        f64
    ),
    (
        onlyHeavyAtomsForRMS,
        get_only_heavy_atoms_for_rms,
        set_only_heavy_atoms_for_rms,
        only_heavy_atoms_for_rms,
        bool
    ),
    (
        embedFragmentsSeparately,
        get_embed_fragments_separately,
        set_embed_fragments_separately,
        embed_fragments_separately,
        bool
    ),
    (
        useSmallRingTorsions,
        get_use_small_ring_torsions,
        set_use_small_ring_torsions,
        use_small_ring_torsions,
        bool
    ),
    (
        useMacrocycleTorsions,
        get_use_macrocycle_torsions,
        set_use_macrocycle_torsions,
        use_macrocycle_torsions,
        bool
    ),
    (
        useMacrocycle14config,
        get_use_macrocycle_14_config,
        set_use_macrocycle_14_config,
        use_macrocycle_14_config,
        bool
    ),
    (
        useLegacyImplementation,
        get_use_legacy_implementation,
        set_use_legacy_implementation,
        use_legacy_implementation,
        bool
    ),
    (
        boundsMatForceScaling,
        get_bounds_mat_force_scaling,
        set_bounds_mat_force_scaling,
        bounds_mat_force_scaling,
        f64
    ),
    (
        useSymmetryForPruning,
        get_use_symmetry_for_pruning,
        set_use_symmetry_for_pruning,
        use_symmetry_for_pruning,
        bool
    ),
    (
        forceTransAmides,
        get_force_trans_amides,
        set_force_trans_amides,
        force_trans_amides,
        bool
    ),
    (
        trackFailures,
        get_track_failures,
        set_track_failures,
        track_failures,
        bool
    ),
    (
        enableSequentialRandomSeeds,
        get_enable_sequential_random_seeds,
        set_enable_sequential_random_seeds,
        enable_sequential_random_seeds,
        bool
    ),
    (
        symmetrizeConjugatedTerminalGroupsForPruning,
        get_symmetrize_conjugated_terminal_groups_for_pruning,
        set_symmetrize_conjugated_terminal_groups_for_pruning,
        symmetrize_conjugated_terminal_groups_for_pruning,
        bool
    ),
    (
        onlyInitialEmbedding,
        get_only_initial_embedding,
        set_only_initial_embedding,
        only_initial_embedding,
        bool
    ),
);

fn coord_map_from(dict: &Bound<'_, PyDict>) -> PyResult<CoordMap> {
    let mut map = CoordMap::new();
    for (k, v) in dict.iter() {
        let p = v.extract::<PyRef<'_, PyPoint3D>>()?;
        map.insert(k.extract::<u32>()? as i32, p.p);
    }
    Ok(map)
}

#[pymethods]
impl PyEmbedParameters {
    #[new]
    fn new() -> Self {
        to_py_params(EmbedParameters::default())
    }

    /// Mode for initial embedding: DG_EMBEDDING, INTERNAL_COORDINATE_EMBEDDING, RANDOM_COORDINATE_EMBEDDING
    #[getter(initialEmbeddingMode)]
    fn get_initial_embedding_mode(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(
            py,
            "InitialEmbeddingMode",
            self.params.initial_embedding_mode as i64,
        )
    }
    #[setter(initialEmbeddingMode)]
    fn set_initial_embedding_mode(&mut self, value: i64) -> PyResult<()> {
        self.params.initial_embedding_mode = match value {
            0 => InitialEmbeddingMode::DgEmbedding,
            1 => InitialEmbeddingMode::InternalCoordinateEmbedding,
            2 => InitialEmbeddingMode::RandomCoordinateEmbedding,
            _ => return Err(PyValueError::new_err("bad initial embedding mode")),
        };
        Ok(())
    }

    /// Force Field to use for ideal 1-2 and 1-3 distances.
    #[getter(embedForceField)]
    fn get_embed_force_field(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "EmbedFF", self.params.embed_force_field as i64)
    }
    #[setter(embedForceField)]
    fn set_embed_force_field(&mut self, value: i64) -> PyResult<()> {
        self.params.embed_force_field = match value {
            0 => EmbedFF::Uff,
            1 => EmbedFF::Mmff,
            _ => return Err(PyValueError::new_err("bad force field")),
        };
        Ok(())
    }

    /// set the distance-bounds matrix to be used (no triangle smoothing will be
    /// done on this) from a Numpy array
    #[pyo3(name = "SetBoundsMat")]
    fn set_bounds_mat(
        &mut self,
        py: Python<'_>,
        bounds_mat_arg: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let (values, shape) = float64_values(py, bounds_mat_arg)
            .map_err(|_| PyValueError::new_err("Argument isn't an array"))?;
        if shape.len() != 2 || shape[0] != shape[1] {
            return Err(PyValueError::new_err("The array has to be square"));
        }
        let n = shape[0];
        if n == 0 {
            return Err(PyValueError::new_err(
                "The array has to have a nonzero size",
            ));
        }
        let mut mat = BoundsMatrix::new(n);
        for i in 0..n {
            for j in 0..n {
                mat.set_val(i, j, values[i * n + j]);
            }
        }
        self.params.bounds_mat = Some(Arc::new(mat));
        Ok(())
    }

    /// set the customised pairwise Columb-like interaction to atom pairs.used
    /// during structural minimisation stage
    #[pyo3(name = "SetCPCI")]
    fn set_cpci(&mut self, cpci_dict: &Bound<'_, PyDict>) -> PyResult<()> {
        let mut map = std::collections::BTreeMap::new();
        for (k, v) in cpci_dict.iter() {
            let (a, b): (u32, u32) = k.extract()?;
            map.insert((a, b), v.extract::<f64>()?);
        }
        self.params.cpci = Some(Arc::new(map));
        Ok(())
    }

    /// returns the counts of each failure type
    #[pyo3(name = "GetFailureCounts")]
    fn get_failure_counts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, &self.params.failures)
    }

    /// sets the coordmap to be used
    #[pyo3(name = "SetCoordMap")]
    fn set_coord_map(&mut self, cmap: &Bound<'_, PyDict>) -> PyResult<()> {
        self.params.coord_map = Some(coord_map_from(cmap)?);
        Ok(())
    }

    /// Only the documented attributes can be set (upstream's `safeSetattr`).
    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        crate::pyconv::set_declared_attribute(slf.as_any(), name, value)
    }
}

fn params_object<'py>(obj: &Bound<'py, PyAny>) -> Option<Bound<'py, PyEmbedParameters>> {
    obj.cast::<PyEmbedParameters>().ok().cloned()
}

fn to_py_params(params: EmbedParameters) -> PyEmbedParameters {
    PyEmbedParameters { params }
}

macro_rules! preset {
    ($fn_name:ident, $py:literal, $ctor:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyfunction]
        #[pyo3(name = $py)]
        fn $fn_name() -> PyEmbedParameters {
            to_py_params(EmbedParameters::$ctor())
        }
    };
}

preset!(
    py_etkdg,
    "ETKDG",
    etkdg,
    "Returns an EmbedParameters object for the ETKDG method - version 1."
);
preset!(
    py_etkdg_v2,
    "ETKDGv2",
    etkdg_v2,
    "Returns an EmbedParameters object for the ETKDG method - version 2."
);
preset!(
    py_sr_etkdg_v3,
    "srETKDGv3",
    sr_etkdg_v3,
    "Returns an EmbedParameters object for the ETKDG method - version 3 (small rings)."
);
preset!(
    py_etkdg_v3,
    "ETKDGv3",
    etkdg_v3,
    "Returns an EmbedParameters object for the ETKDG method - version 3 (macrocycles)."
);
preset!(
    py_etdg,
    "ETDG",
    etdg,
    "Returns an EmbedParameters object for the ETDG method."
);
preset!(
    py_etdg_v2,
    "ETDGv2",
    etdg_v2,
    "Returns an EmbedParameters object for the ETDG method - version 2."
);
preset!(
    py_kdg,
    "KDG",
    kdg,
    "Returns an EmbedParameters object for the KDG method."
);
preset!(
    py_dg,
    "DG",
    dg,
    "Returns an EmbedParameters object for plain distance geometry."
);
preset!(
    py_etkdg_v4,
    "ETKDGv4",
    etkdg_v4,
    "Returns an EmbedParameters object for the ETKDG method - version 4."
);

#[allow(unused_assignments)] // the last `take!` advances `pos` for nothing
fn keyword_params(
    a: &Args<'_, '_>,
    first: usize,
    num_confs_pos: Option<usize>,
) -> PyResult<EmbedParameters> {
    // the legacy keyword-argument form of EmbedMolecule / EmbedMultipleConfs
    let mut p = EmbedParameters::default();
    let mut pos = first;
    macro_rules! take {
        ($key:literal, $field:ident, $ty:ty, $default:expr) => {{
            p.$field = match a.get(pos, $key)? {
                Some(v) => v.extract::<$ty>()?,
                None => $default,
            };
            pos += 1;
        }};
    }
    take!("maxAttempts", max_iterations, u32, 0);
    take!("randomSeed", random_seed, i32, -1);
    take!("clearConfs", clear_confs, bool, true);
    take!("useRandomCoords", use_random_coords, bool, false);
    take!("boxSizeMult", box_size_mult, f64, 2.0);
    take!("randNegEig", rand_neg_eig, bool, true);
    take!("numZeroFail", num_zero_fail, u32, 1);
    if num_confs_pos.is_some() {
        take!("pruneRmsThresh", prune_rms_thresh, f64, -1.0);
    }
    if let Some(v) = a.get(pos, "coordMap")? {
        let d = v.cast::<PyDict>()?;
        if !d.is_empty() {
            p.coord_map = Some(coord_map_from(d)?);
        }
    }
    pos += 1;
    take!("forceTol", optimizer_force_tol, f64, 1e-3);
    take!(
        "ignoreSmoothingFailures",
        ignore_smoothing_failures,
        bool,
        false
    );
    take!("enforceChirality", enforce_chirality, bool, true);
    if num_confs_pos.is_some() {
        take!("numThreads", num_threads, i32, 1);
    }
    take!(
        "useExpTorsionAnglePrefs",
        use_exp_torsion_angle_prefs,
        bool,
        true
    );
    take!("useBasicKnowledge", use_basic_knowledge, bool, true);
    take!("printExpTorsionAngles", verbose, bool, false);
    take!("useSmallRingTorsions", use_small_ring_torsions, bool, false);
    take!("useMacrocycleTorsions", use_macrocycle_torsions, bool, true);
    take!("ETversion", et_version, u32, 2);
    take!(
        "useMacrocycle14config",
        use_macrocycle_14_config,
        bool,
        true
    );
    p.only_heavy_atoms_for_rms = false;
    Ok(p)
}

fn mol_arg<'py>(a: &Args<'_, 'py>) -> PyResult<Bound<'py, PyMol>> {
    Ok(a.get(0, "mol")?
        .ok_or_else(|| pyo3::exceptions::PyTypeError::new_err("missing argument 'mol'"))?
        .cast_into::<PyMol>()?)
}

/// Use distance geometry to obtain initial coordinates for a molecule
///
///  ARGUMENTS:
///
///     - mol : the molecule of interest
///     - maxAttempts : maximum number of embedding attempts to use for a single conformation
///     - randomSeed : provide a seed for the random number generator so that the same
///                    coordinates can be obtained for a molecule on multiple runs.
///                    If -1, the RNG will not be seeded.
///     - clearConfs : clear all existing conformations on the molecule
///     - useRandomCoords : start the embedding from random coordinates
///     - boxSizeMult, randNegEig, numZeroFail, coordMap, forceTol,
///       ignoreSmoothingFailures, enforceChirality, useExpTorsionAnglePrefs,
///       useBasicKnowledge, printExpTorsionAngles, useSmallRingTorsions,
///       useMacrocycleTorsions, ETversion, useMacrocycle14config : see the
///       corresponding EmbedParameters attributes
///
///  The second argument can also be an EmbedParameters object (`params`).
///
///  RETURNS:
///
///     ID of the new conformation added to the molecule or -1 if the embedding fails.
#[pyfunction]
#[pyo3(name = "EmbedMolecule", signature = (*args, **kwargs))]
fn py_embed_molecule(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let a = Args {
        name: "EmbedMolecule",
        args,
        kwargs,
    };
    let mol = mol_arg(&a)?;
    let params_obj = a.get(1, "params")?.and_then(|p| params_object(&p));
    let mut mol = mol.borrow_mut();
    let core = &mut mol.mol;
    match params_obj {
        Some(obj) => {
            let mut p = obj.borrow_mut();
            let params = &mut p.params;
            py.detach(|| embed_molecule(core, params)).py()
        }
        None => {
            let mut params = keyword_params(&a, 1, None)?;
            py.detach(|| embed_molecule(core, &mut params)).py()
        }
    }
}

/// Use distance geometry to obtain multiple sets of coordinates for a molecule
///
///  ARGUMENTS:
///
///   - mol : the molecule of interest
///   - numConfs : the number of conformers to generate
///   - maxAttempts, randomSeed, clearConfs, useRandomCoords, boxSizeMult,
///     randNegEig, numZeroFail, pruneRmsThresh, coordMap, forceTol,
///     ignoreSmoothingFailures, enforceChirality, numThreads,
///     useExpTorsionAnglePrefs, useBasicKnowledge, printExpTorsionAngles, ... :
///     see EmbedMolecule
///
///  The third argument can also be an EmbedParameters object (`params`).
///
///  RETURNS:
///
///     the IDs of the new conformations
#[pyfunction]
#[pyo3(name = "EmbedMultipleConfs", signature = (*args, **kwargs))]
fn py_embed_multiple_confs(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<i32>> {
    let a = Args {
        name: "EmbedMultipleConfs",
        args,
        kwargs,
    };
    let mol = mol_arg(&a)?;
    let num_confs = a.get(1, "numConfs")?.map_or(Ok(10u32), |v| v.extract())?;
    let params_obj = a.get(2, "params")?.and_then(|p| params_object(&p));
    let mut mol = mol.borrow_mut();
    let core = &mut mol.mol;
    match params_obj {
        Some(obj) => {
            let mut p = obj.borrow_mut();
            let params = &mut p.params;
            py.detach(|| embed_multiple_confs(core, num_confs, params))
                .py()
        }
        None => {
            let mut params = keyword_params(&a, 2, Some(1))?;
            py.detach(|| embed_multiple_confs(core, num_confs, &mut params))
                .py()
        }
    }
}

/// Returns the distance bounds matrix for a molecule
///
///  ARGUMENTS:
///
///     - mol : the molecule of interest
///     - set15bounds : set bounds for 1-5 atom distances based on topology
///                     (otherwise stop at 1-4s)
///     - scaleVDW : ignored
///     - doTriangleSmoothing : run triangle smoothing on the bounds matrix
///                  before returning it
///     - useMacrocycle14config : use 1-4 distance bound heuristics for macrocycles
///     - forceTransAmides : constrain amide bonds to be trans
///     - set14bounds : set bounds for 1-4 atom distances based on topology
///     - set13bounds : set bounds for 1-3 atom distances based on topology
///
///  The second argument can also be an EmbedParameters object (`embedParams`),
///  followed by doTriangleSmoothing, scaleVDW, set15bounds, set14bounds and
///  set13bounds.
///
///  RETURNS:
///
///     the bounds matrix as a Numeric array with lower bounds in the lower
///     triangle and upper bounds in the upper triangle
#[pyfunction]
#[pyo3(name = "GetMoleculeBoundsMatrix", signature = (*args, **kwargs))]
fn py_get_molecule_bounds_matrix<'py>(
    py: Python<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyAny>> {
    let a = Args {
        name: "GetMoleculeBoundsMatrix",
        args,
        kwargs,
    };
    let mol = mol_arg(&a)?;
    let params_obj = a.get(1, "embedParams")?.and_then(|p| params_object(&p));
    let (params, smooth, options) = match params_obj {
        Some(obj) => {
            let params = obj.borrow().params.clone();
            let options = TopolBoundsOptions {
                scale_vdw: a.flag(3, "scaleVDW", false)?,
                set15_bounds: a.flag(4, "set15bounds", true)?,
                set14_bounds: a.flag(5, "set14bounds", true)?,
                set13_bounds: a.flag(6, "set13bounds", true)?,
                embed_force_field: EmbedFF::Uff,
            };
            (params, a.flag(2, "doTriangleSmoothing", true)?, options)
        }
        None => {
            let params = EmbedParameters {
                use_macrocycle_14_config: a.flag(4, "useMacrocycle14config", false)?,
                force_trans_amides: a.flag(5, "forceTransAmides", true)?,
                ..Default::default()
            };
            let options = TopolBoundsOptions {
                scale_vdw: a.flag(2, "scaleVDW", false)?,
                set15_bounds: a.flag(1, "set15bounds", true)?,
                set14_bounds: a.flag(6, "set14bounds", true)?,
                set13_bounds: a.flag(7, "set13bounds", true)?,
                embed_force_field: EmbedFF::Uff,
            };
            (params, a.flag(3, "doTriangleSmoothing", true)?, options)
        }
    };
    let embed_ff = params.embed_force_field;
    let options = TopolBoundsOptions {
        embed_force_field: embed_ff,
        ..options
    };
    let mol = mol.borrow();
    let n = mol.mol.num_atoms() as usize;
    let mut mmat = BoundsMatrix::new(n);
    init_bounds_mat(&mut mmat, 0.0, 1000.0).py()?;
    set_topol_bounds(
        &mol.mol,
        &mut mmat,
        &params,
        &options,
        TopolBoundsOutput::default(),
    )
    .py()?;
    if smooth {
        triangle_smooth_bounds(&mut mmat, 0.0);
    }
    let values: Vec<f64> = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j)))
        .map(|(i, j)| mmat.val(i, j))
        .collect();
    float64_array(py, &values)?.call_method1("reshape", (n, n))
}

fn torsion_dicts<'py>(
    py: Python<'py>,
    bonds: &[rdoxide::graphmol::forcefield_helpers::crystal_ff::TorsionBond],
) -> PyResult<Bound<'py, PyTuple>> {
    let result = PyList::empty(py);
    for (bond_idx, atoms, angle) in bonds {
        let d = PyDict::new(py);
        d.set_item("bondIndex", bond_idx)?;
        match angle {
            TorsionAngleRef::Cosine(a) => {
                d.set_item("torsionIndex", a.torsion_idx)?;
                d.set_item("smarts", &a.smarts)?;
                d.set_item("V", &a.v)?;
                d.set_item("signs", &a.signs)?;
            }
            TorsionAngleRef::Gaussian(a) => {
                d.set_item("torsionIndex", a.torsion_idx)?;
                d.set_item("smarts", &a.smarts)?;
                d.set_item("positions", &a.positions)?;
                d.set_item("widths", &a.widths)?;
                d.set_item("heights", &a.heights)?;
            }
        }
        d.set_item("atomIndices", atoms.to_vec())?;
        result.append(d)?;
    }
    PyTuple::new(py, result)
}

/// returns information about the bonds corresponding to experimental torsions
///
///  ARGUMENTS: (mol, useExpTorsionAnglePrefs=True, useSmallRingTorsions=False,
///  useMacrocycleTorsions=True, useBasicKnowledge=True, ETversion=2,
///  printExpTorsionAngles=False) or (mol, embedParams)
#[pyfunction]
#[pyo3(name = "GetExperimentalTorsions", signature = (*args, **kwargs))]
fn py_get_experimental_torsions<'py>(
    py: Python<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let a = Args {
        name: "GetExperimentalTorsions",
        args,
        kwargs,
    };
    let mol = mol_arg(&a)?;
    let options = match a.get(1, "embedParams")?.and_then(|p| params_object(&p)) {
        Some(obj) => {
            let p = obj.borrow();
            TorsionOptions {
                use_exp_torsions: p.params.use_exp_torsion_angle_prefs,
                use_small_ring_torsions: p.params.use_small_ring_torsions,
                use_macrocycle_torsions: p.params.use_macrocycle_torsions,
                use_basic_knowledge: p.params.use_basic_knowledge,
                version: p.params.et_version,
                verbose: p.params.verbose,
            }
        }
        None => TorsionOptions {
            use_exp_torsions: a.flag(1, "useExpTorsionAnglePrefs", true)?,
            use_small_ring_torsions: a.flag(2, "useSmallRingTorsions", false)?,
            use_macrocycle_torsions: a.flag(3, "useMacrocycleTorsions", true)?,
            use_basic_knowledge: a.flag(4, "useBasicKnowledge", true)?,
            version: a.int(5, "ETversion", 2)? as u32,
            verbose: a.flag(6, "printExpTorsionAngles", false)?,
        },
    };
    if !matches!(options.version, 1 | 2 | 4) {
        return Err(PyValueError::new_err("ETversion can only be 1, 2 or 4."));
    }
    let mol = mol.borrow();
    let mut details = CrystalFFDetails::default();
    let bonds = get_experimental_torsions_with_bonds(&mol.mol, &mut details, &options).py()?;
    torsion_dicts(py, &bonds)
}

/// Returns the embedding parameters as a JSON string
///
///  ARGUMENTS:
///
///     - embedParameters : the EmbedParameters object
///
///  RETURNS:
///
///    The Params object as json string
#[pyfunction]
#[pyo3(name = "EmbedParametersToJSON", signature = (embedParameters))]
#[allow(non_snake_case)]
fn py_embed_parameters_to_json(embedParameters: &PyEmbedParameters) -> String {
    embed_parameters_to_json(&embedParameters.params)
}

/// The failure causes in the order they are checked.
#[pyfunction]
#[pyo3(name = "OrderedEmbedFailureCauses", signature = (legacyImplementation = true))]
#[allow(non_snake_case)]
fn py_ordered_embed_failure_causes<'py>(
    py: Python<'py>,
    legacyImplementation: bool,
) -> PyResult<Bound<'py, PyList>> {
    let causes: &[(i64, i64)] = if legacyImplementation {
        &[
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 3),
            (4, 4),
            (5, 5),
            (8, 8),
            (10, 10),
            (6, 6),
            (7, 7),
            (9, 9),
            (11, 11),
        ]
    } else {
        &[
            (0, 0),
            (12, 12),
            (13, 13),
            (2, 2),
            (3, 3),
            (6, 6),
            (7, 7),
            (9, 9),
            (14, 14),
            (8, 9),
            (11, 11),
        ]
    };
    let result = PyList::empty(py);
    for &(cause, idx) in causes {
        result.append(PyTuple::new(
            py,
            [
                to_py(py, "EmbedFailureCauses", cause)?.into_bound(py),
                idx.into_pyobject(py)?.into_any(),
            ],
        )?)?;
    }
    Ok(result)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to compute atomic coordinates in 3D using distance geometry",
    )?;
    add_enum_exported(
        m,
        "InitialEmbeddingMode",
        &[
            ("DG_EMBEDDING", 0),
            ("INTERNAL_COORDINATE_EMBEDDING", 1),
            ("RANDOM_COORDINATE_EMBEDDING", 2),
        ],
        None,
    )?;
    add_enum_exported(
        m,
        "EmbedFailureCauses",
        &[
            ("INITIAL_COORDS", 0),
            ("FIRST_MINIMIZATION", 1),
            ("CHECK_TETRAHEDRAL_CENTERS", 2),
            ("CHECK_CHIRAL_CENTERS", 3),
            ("MINIMIZE_FOURTH_DIMENSION", 4),
            ("ETK_MINIMIZATION", 5),
            ("FINAL_CHIRAL_BOUNDS", 6),
            ("FINAL_CENTER_IN_VOLUME", 7),
            ("LINEAR_DOUBLE_BOND", 8),
            ("BAD_DOUBLE_BOND_STEREO", 9),
            ("CHECK_CHIRAL_CENTERS2", 10),
            ("EXCEEDED_TIMEOUT", 11),
            ("MINIMIZATION", 12),
            ("KTERM_VIOLATION", 13),
            ("CLASH", 14),
        ],
        None,
    )?;
    add_enum_exported(m, "EmbedFF", &[("UFF", 0), ("MMFF", 1)], None)?;
    m.add_class::<PyEmbedParameters>()?;
    for f in [
        wrap_pyfunction!(py_embed_molecule, m)?,
        wrap_pyfunction!(py_embed_multiple_confs, m)?,
        wrap_pyfunction!(py_get_molecule_bounds_matrix, m)?,
        wrap_pyfunction!(py_get_experimental_torsions, m)?,
        wrap_pyfunction!(py_ordered_embed_failure_causes, m)?,
        wrap_pyfunction!(py_embed_parameters_to_json, m)?,
        wrap_pyfunction!(py_etkdg, m)?,
        wrap_pyfunction!(py_etkdg_v2, m)?,
        wrap_pyfunction!(py_sr_etkdg_v3, m)?,
        wrap_pyfunction!(py_etkdg_v3, m)?,
        wrap_pyfunction!(py_etdg, m)?,
        wrap_pyfunction!(py_etdg_v2, m)?,
        wrap_pyfunction!(py_kdg, m)?,
        wrap_pyfunction!(py_dg, m)?,
        wrap_pyfunction!(py_etkdg_v4, m)?,
    ] {
        m.add_function(f)?;
    }
    Ok(())
}
