//! `rdkit.Chem.rdGaussianShape`: Gaussian-based shape overlay and scoring
//! (upstream `GraphMol/GaussianShape/Wrap/rdGaussianShape.cpp`).

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::geometry::Transform3D;
use rdoxide::graphmol::gaussian_shape::{
    CustomFeature, OptimMode, ShapeInput, ShapeInputOptions, ShapeOverlayOptions, StartMode,
    align_molecule, align_molecule_to_shape, align_shape, score_molecule,
    score_molecule_all_conformers, score_molecule_to_shape, score_shape,
};

use crate::enums::{add_enum_exported, to_py};
use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;
use crate::rdmif::gather;

const START_MODES: &[(&str, i64)] = &[
    ("ROTATE_0", 0),
    ("ROTATE_180", 1),
    ("ROTATE_180_WIGGLE", 2),
    ("ROTATE_45", 3),
    ("ROTATE_0_FRAGMENT", 4),
    ("ROTATE_180_FRAGMENT", 5),
    ("ROTATE_45_FRAGMENT", 6),
    ("A_LA_PUBCHEM", 7),
];

const OPTIM_MODES: &[(&str, i64)] = &[
    ("SHAPE_ONLY", 0),
    ("SHAPE_PLUS_COLOR_SCORE", 1),
    ("SHAPE_PLUS_COLOR", 2),
];

fn start_mode_value(m: StartMode) -> i64 {
    match m {
        StartMode::Rotate0 => 0,
        StartMode::Rotate180 => 1,
        StartMode::Rotate180Wiggle => 2,
        StartMode::Rotate45 => 3,
        StartMode::Rotate0Fragment => 4,
        StartMode::Rotate180Fragment => 5,
        StartMode::Rotate45Fragment => 6,
        StartMode::ALaPubchem => 7,
    }
}

fn start_mode_from(v: i64) -> PyResult<StartMode> {
    Ok(match v {
        0 => StartMode::Rotate0,
        1 => StartMode::Rotate180,
        2 => StartMode::Rotate180Wiggle,
        3 => StartMode::Rotate45,
        4 => StartMode::Rotate0Fragment,
        5 => StartMode::Rotate180Fragment,
        6 => StartMode::Rotate45Fragment,
        7 => StartMode::ALaPubchem,
        _ => {
            return Err(PyValueError::new_err(format!(
                "invalid StartMode value: {v}"
            )));
        }
    })
}

fn optim_mode_value(m: OptimMode) -> i64 {
    match m {
        OptimMode::ShapeOnly => 0,
        OptimMode::ShapePlusColorScore => 1,
        OptimMode::ShapePlusColor => 2,
    }
}

fn optim_mode_from(v: i64) -> PyResult<OptimMode> {
    Ok(match v {
        0 => OptimMode::ShapeOnly,
        1 => OptimMode::ShapePlusColorScore,
        2 => OptimMode::ShapePlusColor,
        _ => {
            return Err(PyValueError::new_err(format!(
                "invalid OptimMode value: {v}"
            )));
        }
    })
}

/// ShapeInputOptions - options for setting up ShapeInput objects.
#[pyclass(
    name = "ShapeInputOptions",
    module = "rdkit.Chem.rdGaussianShape",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyShapeInputOptions {
    opts: ShapeInputOptions,
}

#[pymethods]
impl PyShapeInputOptions {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// Whether to use color features in overlay.  Default=True.
    #[getter(useColors)]
    fn use_colors(&self) -> bool {
        self.opts.use_colors
    }
    #[setter(useColors)]
    fn set_use_colors(&mut self, v: bool) {
        self.opts.use_colors = v;
    }

    /// Whether to use the same radius, appropriate for Carbon, for all atoms.
    /// There is a slight accuracy penalty but significant speed gain if used.
    /// Default=True.
    #[getter(allCarbonRadii)]
    fn all_carbon_radii(&self) -> bool {
        self.opts.all_carbon_radii
    }
    #[setter(allCarbonRadii)]
    fn set_all_carbon_radii(&mut self, v: bool) {
        self.opts.all_carbon_radii = v;
    }

    /// If not empty, use just these atoms in the molecule to form the
    /// ShapeInput object.
    #[getter(atomSubset)]
    fn atom_subset<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.opts.atom_subset.iter().copied())
    }
    #[setter(atomSubset)]
    fn set_atom_subset(&mut self, v: Vec<u32>) {
        self.opts.atom_subset = v;
    }

    /// Custom features for the shape.  Requires a list of lists of tuples of
    /// int (the feature type), Point3D (the coordinates), float (the radius)
    /// and optionally a list of indices of the atoms that the feature was
    /// derived from.
    #[getter(customFeatures)]
    fn custom_features<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let all = PyList::empty(py);
        for feats in &self.opts.custom_features {
            let list = PyList::empty(py);
            for f in feats {
                let elem = PyList::empty(py);
                elem.append(f.kind as i32)?;
                elem.append(PyPoint3D { p: f.pos })?;
                elem.append(f.rad)?;
                elem.append(f.atoms.clone())?;
                list.append(elem)?;
            }
            all.append(list)?;
        }
        PyTuple::new(py, all)
    }
    #[setter(customFeatures)]
    fn set_custom_features(&mut self, v: &Bound<'_, PyAny>) -> PyResult<()> {
        let mut all = Vec::new();
        for confs in v.try_iter()? {
            let mut feats = Vec::new();
            for feat in confs?.try_iter()? {
                let feat = feat?;
                let kind: u32 = feat.get_item(0)?.extract()?;
                let pos = feat.get_item(1)?.extract::<PyRef<'_, PyPoint3D>>()?.p;
                let rad: f64 = feat.get_item(2)?.extract()?;
                let atoms = if feat.len()? == 4 {
                    feat.get_item(3)?.extract::<Vec<u32>>()?
                } else {
                    Vec::new()
                };
                feats.push(CustomFeature {
                    kind,
                    pos,
                    rad,
                    atoms,
                });
            }
            all.push(feats);
        }
        self.opts.custom_features = all;
        Ok(())
    }

    /// Non-standard radii to use for the atoms specified by their indices in
    /// the molecule.  Not all atoms need have a radius specified.  A list of
    /// tuples of [int, float].
    #[getter(atomRadii)]
    fn atom_radii<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.opts.atom_radii.iter().map(|&(i, r)| (i as i32, r)))
    }
    #[setter(atomRadii)]
    fn set_atom_radii(&mut self, v: Vec<(u32, f64)>) {
        self.opts.atom_radii = v;
    }

    /// If there is more than 1 conformer for the input molecule, prune the
    /// shapes so that none of them are more similar to each other than the
    /// threshold.  Default -1.0 means no pruning.
    #[getter(shapePruneThreshold)]
    fn shape_prune_threshold(&self) -> f64 {
        self.opts.shape_prune_threshold
    }
    #[setter(shapePruneThreshold)]
    fn set_shape_prune_threshold(&mut self, v: f64) {
        self.opts.shape_prune_threshold = v;
    }

    /// If True (the default), the shapes are sorted into descending order of
    /// total volume.
    #[getter(sortShapes)]
    fn sort_shapes(&self) -> bool {
        self.opts.sort_shapes
    }
    #[setter(sortShapes)]
    fn set_sort_shapes(&mut self, v: bool) {
        self.opts.sort_shapes = v;
    }

    /// Whether to include dummy atoms in the shape or not.  Default=True.
    #[getter(includeDummies)]
    fn include_dummies(&self) -> bool {
        self.opts.include_dummies
    }
    #[setter(includeDummies)]
    fn set_include_dummies(&mut self, v: bool) {
        self.opts.include_dummies = v;
    }
}

/// ShapeOverlayOptions - options for controlling the shape overlay process.
#[pyclass(
    name = "ShapeOverlayOptions",
    module = "rdkit.Chem.rdGaussianShape",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyShapeOverlayOptions {
    opts: ShapeOverlayOptions,
}

#[pymethods]
impl PyShapeOverlayOptions {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// Start modes for optimisation.  Default is A_LA_PUBCHEM - as used by
    /// the PubChem code - either ROTATE_180_WIGGLE or ROTATE_45 depending on
    /// the shape of the two molecules.
    #[getter(startMode)]
    fn start_mode(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "StartMode", start_mode_value(self.opts.start_mode))
    }
    #[setter(startMode)]
    fn set_start_mode(&mut self, v: i64) -> PyResult<()> {
        self.opts.start_mode = start_mode_from(v)?;
        Ok(())
    }

    /// Optimisation mode, controlling what parameters are used to drive the
    /// overlay.  Default=SHAPE_PLUS_COLOR_SCORE.
    #[getter(optimMode)]
    fn optim_mode(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "OptimMode", optim_mode_value(self.opts.optim_mode))
    }
    #[setter(optimMode)]
    fn set_optim_mode(&mut self, v: i64) -> PyResult<()> {
        self.opts.optim_mode = optim_mode_from(v)?;
        Ok(())
    }

    /// When doing a Tversky similarity, the alpha value.
    #[getter(simAlpha)]
    fn sim_alpha(&self) -> f64 {
        self.opts.sim_alpha
    }
    #[setter(simAlpha)]
    fn set_sim_alpha(&mut self, v: f64) {
        self.opts.sim_alpha = v;
    }

    /// When doing a Tversky similarity, the beta value.
    #[getter(simBeta)]
    fn sim_beta(&self) -> f64 {
        self.opts.sim_beta
    }
    #[setter(simBeta)]
    fn set_sim_beta(&mut self, v: f64) {
        self.opts.sim_beta = v;
    }

    /// If using colors, the relative weights of the shape and color scores,
    /// as a fraction of 1.  Default=0.5.
    #[getter(optParam)]
    fn opt_param(&self) -> f64 {
        self.opts.opt_param
    }
    #[setter(optParam)]
    fn set_opt_param(&mut self, v: f64) {
        self.opts.opt_param = v;
    }

    /// Maximum number of steps for the shape overlay process. Default=100.
    #[getter(nSteps)]
    fn n_steps(&self) -> u32 {
        self.opts.n_steps
    }
    #[setter(nSteps)]
    fn set_n_steps(&mut self, v: u32) {
        self.opts.n_steps = v;
    }

    /// Whether to normalize the shapes before overlay by putting them into
    /// their canonical orientation.  Default=True.
    #[getter(normalize)]
    fn normalize(&self) -> bool {
        self.opts.normalize
    }
    #[setter(normalize)]
    fn set_normalize(&mut self, v: bool) {
        self.opts.normalize = v;
    }

    /// Whether to use distance cutoff when calculating the shape volumes.
    /// Default=True.
    #[getter(useDistCutoff)]
    fn use_dist_cutoff(&self) -> bool {
        self.opts.use_dist_cutoff
    }
    #[setter(useDistCutoff)]
    fn set_use_dist_cutoff(&mut self, v: bool) {
        self.opts.use_dist_cutoff = v;
    }

    /// If using a distance cutoff, this is the value used.  Default=4.5.
    #[getter(distCutoff)]
    fn dist_cutoff(&self) -> f64 {
        self.opts.dist_cutoff
    }
    #[setter(distCutoff)]
    fn set_dist_cutoff(&mut self, v: f64) {
        self.opts.dist_cutoff = v;
    }

    /// Optimisation stops when the shape Tversky score changes by less than
    /// this amount after an optimisation step.  Default=0.001.
    #[getter(shapeConvergenceCriterion)]
    fn shape_convergence_criterion(&self) -> f64 {
        self.opts.shape_convergence_criterion
    }
    #[setter(shapeConvergenceCriterion)]
    fn set_shape_convergence_criterion(&mut self, v: f64) {
        self.opts.shape_convergence_criterion = v;
    }
}

/// ShapeInput object
#[pyclass(name = "ShapeInput", module = "rdkit.Chem.rdGaussianShape")]
pub struct PyShapeInput {
    shape: ShapeInput,
}

fn input_opts(o: Option<PyRef<'_, PyShapeInputOptions>>) -> ShapeInputOptions {
    o.map(|o| o.opts.clone()).unwrap_or_default()
}

fn overlay_opts(o: Option<PyRef<'_, PyShapeOverlayOptions>>) -> ShapeOverlayOptions {
    o.map(|o| o.opts.clone()).unwrap_or_default()
}

fn matrix(m: &Transform3D) -> Vec<f64> {
    (0..4)
        .flat_map(|i| (0..4).map(move |j| (i, j)))
        .map(|(i, j)| m.get_val(i, j))
        .collect()
}

#[pymethods]
impl PyShapeInput {
    #[new]
    #[pyo3(signature = (mol, confId = -1, shapeOpt = None, overlayOpts = None))]
    #[allow(non_snake_case)]
    fn py_new(
        py: Python<'_>,
        mol: PyRef<'_, PyMol>,
        confId: i32,
        shapeOpt: Option<PyRef<'_, PyShapeInputOptions>>,
        overlayOpts: Option<PyRef<'_, PyShapeOverlayOptions>>,
    ) -> PyResult<Self> {
        let (so, oo) = (input_opts(shapeOpt), overlay_opts(overlayOpts));
        let mol = &mol.mol;
        let shape = py.detach(|| ShapeInput::new(mol, confId, &so, &oo)).py()?;
        Ok(Self { shape })
    }

    /// Get the SMILES string for the molecule that the shape relates to.
    #[getter(GetSmiles)]
    fn smiles(&self) -> String {
        self.shape.smiles().to_string()
    }

    /// Return the number of the active shape.
    #[getter(getActiveShape)]
    fn active_shape(&self) -> usize {
        self.shape.active_shape()
    }

    /// Set the active shape, the one that will be used for overlays etc.
    #[pyo3(name = "setActiveShape")]
    fn set_active_shape(&mut self, shape_num: usize) -> PyResult<()> {
        self.shape.set_active_shape(shape_num).py()
    }

    /// Get the number of atoms defining the shape.
    #[getter(NumAtoms)]
    fn num_atoms(&self) -> usize {
        self.shape.num_atoms()
    }

    /// Get the number of features in the shape.
    #[getter(NumFeatures)]
    fn num_features(&self) -> usize {
        self.shape.num_features()
    }

    /// Get the number of shapes.  There will be a shape for each conformation
    /// of the input molecule, unless shape pruning was carried out in which
    /// case there may be fewer.
    #[getter(NumShapes)]
    fn num_shapes(&self) -> usize {
        self.shape.num_shapes()
    }

    /// Get the volume due to the atoms for the active shape.
    #[getter(ShapeVolume)]
    fn shape_volume(&self) -> f64 {
        self.shape.shape_volume()
    }

    /// Get the volume of the shape's color features for the active shape.
    #[getter(ColorVolume)]
    fn color_volume(&self) -> f64 {
        self.shape.color_volume()
    }

    /// Align the principal axes to the cartesian axes and centre on the
    /// origin.
    #[pyo3(name = "NormalizeCoords")]
    fn normalize_coords(&mut self) -> PyResult<()> {
        self.shape.normalize_coords().py()
    }

    /// Return a molecule with coordinates of the current active shape.  If
    /// includeColors is True, (default is False) the color features will be
    /// added as xenon atoms.  If withBonds is True (the default) a molecule
    /// with bonds will be created, if not then just atoms at the appropriate
    /// positions will be produced.
    #[pyo3(name = "ShapeToMol", signature = (includeColors = false, withBonds = true))]
    #[allow(non_snake_case)]
    fn shape_to_mol(&self, includeColors: bool, withBonds: bool) -> PyResult<PyMol> {
        Ok(PyMol {
            mol: self.shape.shape_to_mol(includeColors, withBonds).py()?,
        })
    }

    /// Find the best similarity score between all shapes in this shape and
    /// the other one. Stops as soon as it gets something above the threshold.
    /// Returns a tuple of the similarity scores ((-1.0, -1.0, -1.0) if there
    /// was nothing above the threshold), the number of the shape for this
    /// object and the shape number of the fitShape that gave the best
    /// similarity and the transformation matrix (as a list of 16 floats) that
    /// will reproduce the best overlay.
    #[pyo3(name = "BestSimilarity", signature = (fitShape, threshold = -1.0, overlayOpts = None))]
    #[allow(non_snake_case)]
    fn best_similarity<'py>(
        &mut self,
        py: Python<'py>,
        fitShape: PyRef<'_, PyShapeInput>,
        threshold: f64,
        overlayOpts: Option<PyRef<'_, PyShapeOverlayOptions>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let oo = overlay_opts(overlayOpts);
        let (shape, fit) = (&mut self.shape, &fitShape.shape);
        let (sim, best) = py
            .detach(|| shape.best_similarity(fit, threshold, &oo))
            .py()?;
        let (this_shape, fit_shape, xform) = match best {
            Some(b) => (b.this_shape as u32, b.fit_shape as u32, b.xform),
            None => (u32::MAX, u32::MAX, Transform3D::identity()),
        };
        PyTuple::new(
            py,
            [
                PyTuple::new(py, sim)?.into_any(),
                this_shape.into_pyobject(py)?.into_any(),
                fit_shape.into_pyobject(py)?.into_any(),
                matrix(&xform).into_pyobject(py)?.into_any(),
            ],
        )
    }

    /// Get the maximum possible similarity score between all shapes in this
    /// shape and all shapes in the fitShape.  The maximum similarity is when
    /// one shape is entirely inside the other.
    #[pyo3(name = "MaxPossibleSimilarity", signature = (fitShape, overlayOpts = None))]
    #[allow(non_snake_case)]
    fn max_possible_similarity(
        &self,
        fitShape: PyRef<'_, PyShapeInput>,
        overlayOpts: Option<PyRef<'_, PyShapeOverlayOptions>>,
    ) -> PyResult<f64> {
        self.shape
            .max_possible_similarity(&fitShape.shape, &overlay_opts(overlayOpts))
            .py()
    }
}

fn scores_tuple<'py>(py: Python<'py>, s: [f64; 3]) -> PyResult<Vec<Bound<'py, PyAny>>> {
    s.iter()
        .map(|v| Ok(v.into_pyobject(py)?.into_any()))
        .collect()
}

/// A keyword/positional argument that may be missing or None.
fn opt<'py, T: for<'a> FromPyObject<'a, 'py>>(v: &Option<Bound<'py, PyAny>>) -> PyResult<Option<T>>
where
    for<'a> <T as FromPyObject<'a, 'py>>::Error: Into<PyErr>,
{
    match v {
        Some(o) if !o.is_none() => o.extract::<T>().map(Some).map_err(Into::into),
        _ => Ok(None),
    }
}

fn int_or(v: &Option<Bound<'_, PyAny>>, default: i32) -> PyResult<i32> {
    Ok(opt::<i32>(v)?.unwrap_or(default))
}

fn borrow_opts<'a>(
    v: &'a Option<Bound<'_, PyAny>>,
) -> PyResult<Option<PyRef<'a, PyShapeInputOptions>>> {
    match v {
        Some(o) if !o.is_none() => Ok(Some(o.extract()?)),
        _ => Ok(None),
    }
}

fn borrow_overlay<'a>(
    v: &'a Option<Bound<'_, PyAny>>,
) -> PyResult<Option<PyRef<'a, PyShapeOverlayOptions>>> {
    match v {
        Some(o) if !o.is_none() => Ok(Some(o.extract()?)),
        _ => Ok(None),
    }
}

/// Whether a call of `AlignMol`/`ScoreMol` is the overload with a reference
/// shape rather than a reference molecule.
fn is_shape_overload(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<bool> {
    if let Ok(first) = args.get_item(0) {
        return Ok(first.is_instance_of::<PyShapeInput>());
    }
    match kwargs {
        Some(k) => Ok(k.contains("refShape")?),
        None => Err(PyTypeError::new_err("missing the reference")),
    }
}

fn required<'a, 'py>(
    v: &'a Option<Bound<'py, PyAny>>,
    name: &str,
) -> PyResult<&'a Bound<'py, PyAny>> {
    v.as_ref()
        .ok_or_else(|| PyTypeError::new_err(format!("missing required argument '{name}'")))
}

/// Aligns a fit molecule onto a reference molecule (or shape).  The fit is
/// modified.
///
/// AlignMol(ref, fit, refOpts=None, fitOpts=None, overlayOpts=None,
/// refConfId=-1, fitConfId=-1) or AlignMol(refShape, fit, fitOpts=None,
/// overlayOpts=None, fitConfId=-1).
///
/// Returns the 3-tuple (combo_score, shape_score, color_score).  The
/// color_score is 0.0 if color features are not used, in which case
/// combo_score and shape_score will be the same.
#[pyfunction]
#[pyo3(name = "AlignMol", signature = (*args, **kwargs))]
fn py_align_mol<'py>(
    py: Python<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let scores = if is_shape_overload(args, kwargs)? {
        let v = gather(
            args,
            kwargs,
            &["refShape", "fit", "fitOpts", "overlayOpts", "fitConfId"],
        )?;
        let ref_shape = required(&v[0], "refShape")?.extract::<PyRef<'_, PyShapeInput>>()?;
        let mut fit = required(&v[1], "fit")?.extract::<PyRefMut<'_, PyMol>>()?;
        let fo = input_opts(borrow_opts(&v[2])?);
        let oo = overlay_opts(borrow_overlay(&v[3])?);
        let fit_conf = int_or(&v[4], -1)?;
        let (r, f) = (&ref_shape.shape, &mut fit.mol);
        py.detach(|| align_molecule_to_shape(r, f, &fo, &oo, fit_conf))
            .py()?
            .0
    } else {
        let v = gather(
            args,
            kwargs,
            &[
                "ref",
                "fit",
                "refOpts",
                "fitOpts",
                "overlayOpts",
                "refConfId",
                "fitConfId",
            ],
        )?;
        let reference = required(&v[0], "ref")?.extract::<PyRef<'_, PyMol>>()?;
        let mut fit = required(&v[1], "fit")?.extract::<PyRefMut<'_, PyMol>>()?;
        let ro = input_opts(borrow_opts(&v[2])?);
        let fo = input_opts(borrow_opts(&v[3])?);
        let oo = overlay_opts(borrow_overlay(&v[4])?);
        let (rc, fc) = (int_or(&v[5], -1)?, int_or(&v[6], -1)?);
        let (r, f) = (&reference.mol, &mut fit.mol);
        py.detach(|| align_molecule(r, f, &ro, &fo, &oo, rc, fc))
            .py()?
            .0
    };
    PyTuple::new(py, scores_tuple(py, scores)?)
}

/// Aligns a fit shape to a reference shape. The fit is modified.
///
/// Returns the 4-tuple (combo_score, shape_score, color_score, matrix); the
/// matrix is a 16-float list giving the transformation matrix that overlays
/// the fit onto the reference.
#[pyfunction]
#[pyo3(name = "AlignShapes", signature = (refShape, fitShape, overlayOpts = None))]
#[allow(non_snake_case)]
fn py_align_shapes<'py>(
    py: Python<'py>,
    refShape: PyRef<'_, PyShapeInput>,
    mut fitShape: PyRefMut<'_, PyShapeInput>,
    overlayOpts: Option<PyRef<'_, PyShapeOverlayOptions>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let oo = overlay_opts(overlayOpts);
    let (r, f) = (&refShape.shape, &mut fitShape.shape);
    let (scores, xform) = py.detach(|| align_shape(r, f, &oo)).py()?;
    let mut items = scores_tuple(py, scores)?;
    items.push(matrix(&xform).into_pyobject(py)?.into_any());
    PyTuple::new(py, items)
}

fn score_tuple<'py>(
    py: Python<'py>,
    scores: [f64; 3],
    vols: (f64, f64),
) -> PyResult<Bound<'py, PyTuple>> {
    let mut items = scores_tuple(py, scores)?;
    items.push(vols.0.into_pyobject(py)?.into_any());
    items.push(vols.1.into_pyobject(py)?.into_any());
    PyTuple::new(py, items)
}

/// Calculates the scores between a reference molecule (or shape) and a fit
/// molecule without overlay.
///
/// ScoreMol(ref, fit, refOpts=None, fitOpts=None, overlayOpts=None,
/// refConfId=-1, fitConfId=-1) or ScoreMol(refShape, fit, fitOpts=None,
/// overlayOpts=None, fitConfId=-1).
///
/// Returns (combo_score, shape_score, color_score, shape_overlap_volume,
/// color_overlap_volume).
#[pyfunction]
#[pyo3(name = "ScoreMol", signature = (*args, **kwargs))]
fn py_score_mol<'py>(
    py: Python<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (scores, vols) = if is_shape_overload(args, kwargs)? {
        let v = gather(
            args,
            kwargs,
            &["refShape", "fit", "fitOpts", "overlayOpts", "fitConfId"],
        )?;
        let ref_shape = required(&v[0], "refShape")?.extract::<PyRef<'_, PyShapeInput>>()?;
        let fit = required(&v[1], "fit")?.extract::<PyRef<'_, PyMol>>()?;
        let fo = input_opts(borrow_opts(&v[2])?);
        let oo = overlay_opts(borrow_overlay(&v[3])?);
        let fit_conf = int_or(&v[4], -1)?;
        let (r, f) = (&ref_shape.shape, &fit.mol);
        py.detach(|| score_molecule_to_shape(r, f, &fo, &oo, fit_conf))
            .py()?
    } else {
        let v = gather(
            args,
            kwargs,
            &[
                "ref",
                "fit",
                "refOpts",
                "fitOpts",
                "overlayOpts",
                "refConfId",
                "fitConfId",
            ],
        )?;
        let reference = required(&v[0], "ref")?.extract::<PyRef<'_, PyMol>>()?;
        let fit = required(&v[1], "fit")?.extract::<PyRef<'_, PyMol>>()?;
        let ro = input_opts(borrow_opts(&v[2])?);
        let fo = input_opts(borrow_opts(&v[3])?);
        let oo = overlay_opts(borrow_overlay(&v[4])?);
        let (rc, fc) = (int_or(&v[5], -1)?, int_or(&v[6], -1)?);
        let (r, f) = (&reference.mol, &fit.mol);
        py.detach(|| score_molecule(r, f, &ro, &fo, &oo, rc, fc))
            .py()?
    };
    score_tuple(py, scores, vols)
}

/// Calculates the scores between a reference shape and a fit shape without
/// overlay.
///
/// Returns (combo_score, shape_score, color_score, shape_overlap_volume,
/// color_overlap_volume).
#[pyfunction]
#[pyo3(name = "ScoreShape", signature = (refShape, fitShape, overlayOpts = None))]
#[allow(non_snake_case)]
fn py_score_shape<'py>(
    py: Python<'py>,
    refShape: PyRef<'_, PyShapeInput>,
    fitShape: PyRef<'_, PyShapeInput>,
    overlayOpts: Option<PyRef<'_, PyShapeOverlayOptions>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let oo = overlay_opts(overlayOpts);
    let (scores, vols) = score_shape(&refShape.shape, &fitShape.shape, &oo);
    score_tuple(py, scores, vols)
}

/// Calculate the scores for the alignment of all conformers of the fit
/// molecule onto the reference.  The molecules themselves are not altered.
///
/// Returns a tuple of: the scores from aligning the fit conformations onto
/// the reference conformations (scores[0][1] is the score of aligning fit
/// conformation 1 onto ref conformation 0); the ID of the ref conformer and
/// of the fit conformer from the best-scoring alignment; and the
/// transformation that gives the best-scoring alignment as 16 floats.
#[pyfunction]
#[pyo3(
    name = "ScoreMoleculeAllConformers",
    signature = (r#ref, fit, refOpts = None, fitOpts = None, overlayOpts = None)
)]
#[allow(non_snake_case)]
fn py_score_molecule_all_conformers<'py>(
    py: Python<'py>,
    r#ref: PyRef<'_, PyMol>,
    fit: PyRef<'_, PyMol>,
    refOpts: Option<PyRef<'_, PyShapeInputOptions>>,
    fitOpts: Option<PyRef<'_, PyShapeInputOptions>>,
    overlayOpts: Option<PyRef<'_, PyShapeOverlayOptions>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (ro, fo, oo) = (
        input_opts(refOpts),
        input_opts(fitOpts),
        overlay_opts(overlayOpts),
    );
    let (r, f) = (&r#ref.mol, &fit.mol);
    let res = py
        .detach(|| score_molecule_all_conformers(r, f, &ro, &fo, &oo))
        .py()?;
    PyTuple::new(
        py,
        [
            res.comb_scores.into_pyobject(py)?.into_any(),
            res.ref_conf_id.into_pyobject(py)?.into_any(),
            res.fit_conf_id.into_pyobject(py)?.into_any(),
            matrix(&res.xform).into_pyobject(py)?.into_any(),
        ],
    )
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing implementation of Gaussian-based shape overlay and scoring.NOTE: This functionality is experimental and the API and/or results may change in future releases.",
    )?;
    add_enum_exported(m, "StartMode", START_MODES, None)?;
    add_enum_exported(m, "OptimMode", OPTIM_MODES, None)?;
    m.add_class::<PyShapeInputOptions>()?;
    m.add_class::<PyShapeOverlayOptions>()?;
    m.add_class::<PyShapeInput>()?;
    m.add_function(wrap_pyfunction!(py_align_mol, m)?)?;
    m.add_function(wrap_pyfunction!(py_align_shapes, m)?)?;
    m.add_function(wrap_pyfunction!(py_score_mol, m)?)?;
    m.add_function(wrap_pyfunction!(py_score_shape, m)?)?;
    m.add_function(wrap_pyfunction!(py_score_molecule_all_conformers, m)?)?;
    Ok(())
}
