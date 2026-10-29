//! `rdoxide.Chem.rdShapeAlign` (upstream `External/pubchem_shape/Wrap/cshapealign.cpp`).

use pyo3::conversion::FromPyObjectOwned;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::external::pubchem_shape::shape::{
    ShapeInput, ShapeInputOptions, align_molecule, align_molecule_to_shape, align_shape,
    prepare_conformer, score_molecule, score_molecule_to_shape, score_shape,
    transform_conformer_by_overlay,
};

use crate::conformer::PyConformer;
use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;

/// Shape Input Options
#[pyclass(
    name = "ShapeInputOptions",
    module = "rdoxide.Chem.rdShapeAlign",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
struct PyShapeInputOptions {
    opts: ShapeInputOptions,
}

#[pymethods]
impl PyShapeInputOptions {
    #[new]
    fn new() -> Self {
        PyShapeInputOptions::default()
    }

    /// Whether to use colors (pharmacophore features) in the score.  Default=True.
    #[getter(useColors)]
    fn use_colors(&self) -> bool {
        self.opts.use_colors
    }
    #[setter(useColors)]
    fn set_use_colors(&mut self, v: bool) {
        self.opts.use_colors = v;
    }

    /// Whether to use dummy atoms in the alignment. Default=False.
    #[getter(includeDummies)]
    fn include_dummies(&self) -> bool {
        self.opts.include_dummies
    }
    #[setter(includeDummies)]
    fn set_include_dummies(&mut self, v: bool) {
        self.opts.include_dummies = v;
    }

    /// If using dummy atoms in the alignment, what radius to use for them.  Default=2.16 (the radius of Xe).
    #[getter(dummyRadius)]
    fn dummy_radius(&self) -> f64 {
        self.opts.dummy_radius
    }
    #[setter(dummyRadius)]
    fn set_dummy_radius(&mut self, v: f64) {
        self.opts.dummy_radius = v;
    }

    /// If not empty, use just these atoms in the molecule to form the ShapeInput object.
    #[getter(atomSubset)]
    fn atom_subset(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        Ok(PyTuple::new(py, self.opts.atom_subset.iter())?.unbind())
    }
    #[setter(atomSubset)]
    fn set_atom_subset(&mut self, v: Vec<u32>) {
        self.opts.atom_subset = v;
    }

    /// Any atoms mentioned here by index should not be used in a color feature.
    #[getter(notColorAtoms)]
    fn not_color_atoms(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        Ok(PyTuple::new(py, self.opts.not_color_atoms.iter())?.unbind())
    }
    #[setter(notColorAtoms)]
    fn set_not_color_atoms(&mut self, v: Vec<u32>) {
        self.opts.not_color_atoms = v;
    }

    /// Non-standard radii to use for the atoms specified by their indices in the molecule.  A list of tuples of [int, float].
    #[getter(atomRadii)]
    fn atom_radii(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        let items: Vec<(i32, f64)> = self
            .opts
            .atom_radii
            .iter()
            .map(|&(i, r)| (i as i32, r))
            .collect();
        Ok(PyTuple::new(py, items)?.unbind())
    }
    #[setter(atomRadii)]
    fn set_atom_radii(&mut self, v: Vec<(u32, f64)>) {
        self.opts.atom_radii = v;
    }

    /// Custom features for the shape.
    #[getter(customFeatures)]
    fn custom_features(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        let mut items = Vec::new();
        for (t, p, r) in &self.opts.custom_features {
            let elem = PyList::empty(py);
            elem.append(*t as i32)?;
            elem.append(Py::new(py, PyPoint3D { p: *p })?)?;
            elem.append(*r)?;
            items.push(PyTuple::new(py, elem.iter())?.unbind());
        }
        Ok(PyTuple::new(py, items)?.unbind())
    }
    #[setter(customFeatures)]
    fn set_custom_features(&mut self, v: &Bound<'_, PyAny>) -> PyResult<()> {
        self.opts.custom_features.clear();
        for item in v.try_iter()? {
            let item = item?;
            let t: u32 = item.get_item(0)?.extract()?;
            let p: PyRef<'_, PyPoint3D> = item.get_item(1)?.extract()?;
            let r: f64 = item.get_item(2)?.extract()?;
            self.opts.custom_features.push((t, p.p, r));
        }
        Ok(())
    }

    /// Whether to normalise the shape by putting into its inertial frame.  Default=True.
    #[getter]
    fn normalize(&self) -> bool {
        self.opts.normalize
    }
    #[setter(normalize)]
    fn set_normalize(&mut self, v: bool) {
        self.opts.normalize = v;
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        const KNOWN: [&str; 8] = [
            "useColors",
            "includeDummies",
            "dummyRadius",
            "atomSubset",
            "notColorAtoms",
            "atomRadii",
            "customFeatures",
            "normalize",
        ];
        if !KNOWN.contains(&name) {
            return Err(pyo3::exceptions::PyAttributeError::new_err(format!(
                "'ShapeInputOptions' object has no attribute '{name}'"
            )));
        }
        slf.call_method1("__setattr_known__", (name, value))?;
        Ok(())
    }

    #[pyo3(name = "__setattr_known__")]
    fn setattr_known(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        // the generic object protocol would call `__setattr__` again, so route through the descriptors directly
        let ty = slf.get_type();
        let attr = ty.getattr(name)?;
        attr.call_method1("__set__", (slf, value))?;
        Ok(())
    }
}

/// the shape of a molecule
#[pyclass(
    name = "ShapeInput",
    module = "rdoxide.Chem.rdShapeAlign",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
struct PyShapeInput {
    shape: ShapeInput,
}

#[pymethods]
impl PyShapeInput {
    #[getter]
    fn coord(&self) -> Vec<f32> {
        self.shape.coord.clone()
    }
    #[setter(coord)]
    fn set_coord(&mut self, v: Vec<f32>) {
        self.shape.coord = v;
    }
    #[getter]
    fn alpha_vector(&self) -> Vec<f64> {
        self.shape.alpha_vector.clone()
    }
    #[setter(alpha_vector)]
    fn set_alpha_vector(&mut self, v: Vec<f64>) {
        self.shape.alpha_vector = v;
    }
    #[getter]
    fn atom_type_vector(&self) -> Vec<u32> {
        self.shape.atom_type_vector.clone()
    }
    #[setter(atom_type_vector)]
    fn set_atom_type_vector(&mut self, v: Vec<u32>) {
        self.shape.atom_type_vector = v;
    }
    #[getter(volumeAtomIndexVector)]
    fn volume_atom_index_vector(&self) -> Vec<u32> {
        self.shape.volume_atom_index_vector.clone()
    }
    #[setter(volumeAtomIndexVector)]
    fn set_volume_atom_index_vector(&mut self, v: Vec<u32>) {
        self.shape.volume_atom_index_vector = v;
    }
    /// Translation of centre of shape coordinates to origin.
    #[getter]
    fn shift(&self) -> Vec<f64> {
        self.shape.shift.clone()
    }
    #[setter(shift)]
    fn set_shift(&mut self, v: Vec<f64>) {
        self.shape.shift = v;
    }
    /// Rotation applied to put the shape into its principal axes frame of reference.
    #[getter(inertialRot)]
    fn inertial_rot(&self) -> Vec<f64> {
        self.shape.inertial_rot.clone()
    }
    #[getter]
    fn sov(&self) -> f64 {
        self.shape.sov
    }
    #[setter(sov)]
    fn set_sov(&mut self, v: f64) {
        self.shape.sov = v;
    }
    #[getter]
    fn sof(&self) -> f64 {
        self.shape.sof
    }
    #[setter(sof)]
    fn set_sof(&mut self, v: f64) {
        self.shape.sof = v;
    }
}

fn get<'py>(
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    pos: usize,
    name: &str,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    if pos < args.len() {
        return Ok(Some(args.get_item(pos)?));
    }
    match kwargs {
        Some(k) => k.get_item(name),
        None => Ok(None),
    }
}

fn arg<'py, T>(
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    pos: usize,
    name: &str,
    default: T,
) -> PyResult<T>
where
    T: FromPyObjectOwned<'py>,
{
    match get(args, kwargs, pos, name)? {
        Some(v) => v.extract::<T>().map_err(Into::into),
        None => Ok(default),
    }
}

fn arg_required<'py>(
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    pos: usize,
    name: &str,
) -> PyResult<Bound<'py, PyAny>> {
    get(args, kwargs, pos, name)?
        .ok_or_else(|| PyTypeError::new_err(format!("missing argument {name}")))
}

/// Aligns a probe molecule to a reference molecule (or shape). The probe is modified. The forms are
/// AlignMol(ref, probe, refConfId=-1, probeConfId=-1, useColors=True, opt_param=1.0, max_preiters=10, max_postiters=30),
/// AlignMol(ref, probe, refShapeOpts, probeShapeOpts, refConfId=-1, probeConfId=-1, opt_param=1.0, max_preiters=10,
/// max_postiters=30) and
/// AlignMol(refShape, probe, probeConfId=-1, useColors=True, opt_param=1.0, max_preiters=10, max_postiters=30,
/// applyRefShift=False). Returns the tuple (shape_score, color_score).
#[pyfunction]
#[pyo3(name = "AlignMol", signature = (reference, probe, *args, **kwargs))]
fn align_mol(
    reference: &Bound<'_, PyAny>,
    mut probe: PyRefMut<'_, PyMol>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(f64, f64)> {
    let mut matrix = [0f32; 12];
    if let Ok(refshape) = reference.extract::<PyRef<'_, PyShapeInput>>() {
        let conf: i32 = arg(args, kwargs, 0, "probeConfId", -1)?;
        let colors: bool = arg(args, kwargs, 1, "useColors", true)?;
        let opt: f64 = arg(args, kwargs, 2, "opt_param", 1.0)?;
        let pre: u32 = arg(args, kwargs, 3, "max_preiters", 10)?;
        let post: u32 = arg(args, kwargs, 4, "max_postiters", 30)?;
        let apply: bool = arg(args, kwargs, 5, "applyRefShift", false)?;
        let opts = ShapeInputOptions {
            use_colors: colors,
            ..Default::default()
        };
        return align_molecule_to_shape(
            &refshape.shape,
            &mut probe.mol,
            &mut matrix,
            &opts,
            conf,
            opt,
            pre,
            post,
            apply,
        )
        .py();
    }
    let refmol: PyRef<'_, PyMol> = reference.extract()?;
    let third_is_opts = args
        .get_item(0)
        .is_ok_and(|a| a.extract::<PyRef<'_, PyShapeInputOptions>>().is_ok())
        || kwargs.is_some_and(|k| k.contains("refShapeOpts").unwrap_or(false));
    if third_is_opts {
        let ro: PyRef<'_, PyShapeInputOptions> =
            arg_required(args, kwargs, 0, "refShapeOpts")?.extract()?;
        let po: PyRef<'_, PyShapeInputOptions> =
            arg_required(args, kwargs, 1, "probeShapeOpts")?.extract()?;
        let rc: i32 = arg(args, kwargs, 2, "refConfId", -1)?;
        let pc: i32 = arg(args, kwargs, 3, "probeConfId", -1)?;
        let opt: f64 = arg(args, kwargs, 4, "opt_param", 1.0)?;
        let pre: u32 = arg(args, kwargs, 5, "max_preiters", 10)?;
        let post: u32 = arg(args, kwargs, 6, "max_postiters", 30)?;
        align_molecule(
            &refmol.mol,
            &mut probe.mol,
            &mut matrix,
            &ro.opts,
            &po.opts,
            rc,
            pc,
            opt,
            pre,
            post,
        )
        .py()
    } else {
        let rc: i32 = arg(args, kwargs, 0, "refConfId", -1)?;
        let pc: i32 = arg(args, kwargs, 1, "probeConfId", -1)?;
        let colors: bool = arg(args, kwargs, 2, "useColors", true)?;
        let opt: f64 = arg(args, kwargs, 3, "opt_param", 1.0)?;
        let pre: u32 = arg(args, kwargs, 4, "max_preiters", 10)?;
        let post: u32 = arg(args, kwargs, 5, "max_postiters", 30)?;
        let o = ShapeInputOptions {
            use_colors: colors,
            ..Default::default()
        };
        align_molecule(
            &refmol.mol,
            &mut probe.mol,
            &mut matrix,
            &o,
            &o,
            rc,
            pc,
            opt,
            pre,
            post,
        )
        .py()
    }
}

/// Aligns a probe shape to a reference shape. Returns (shape_score, color_score, matrix): the matrix is a 12-float list
/// giving the transformation that overlays the probe onto the reference.
#[pyfunction]
#[pyo3(name = "AlignShapes", signature = (refShape, probeShape, opt_param = 1.0, max_preiters = 10, max_postiters = 30))]
#[allow(non_snake_case)]
fn align_shapes(
    refShape: PyRef<'_, PyShapeInput>,
    mut probeShape: PyRefMut<'_, PyShapeInput>,
    opt_param: f64,
    max_preiters: u32,
    max_postiters: u32,
) -> (f64, f64, Vec<f32>) {
    let mut matrix = [0f32; 12];
    let (st, ct) = align_shape(
        &refShape.shape,
        &mut probeShape.shape,
        &mut matrix,
        opt_param,
        max_preiters,
        max_postiters,
    );
    (st, ct, matrix.to_vec())
}

/// Applies the transformation of an overlay (final translation and rotation, the matrix of AlignShapes) to the conformer.
#[pyfunction]
#[pyo3(name = "TransformConformer", signature = (finalTrans, finalRot, matrix, probeShape, probeConformer))]
#[allow(non_snake_case)]
fn transform_conformer_py(
    py: Python<'_>,
    finalTrans: Vec<f64>,
    finalRot: Vec<f64>,
    matrix: Vec<f32>,
    probeShape: PyRef<'_, PyShapeInput>,
    mut probeConformer: PyRefMut<'_, PyConformer>,
) -> PyResult<()> {
    if matrix.len() != 12 {
        return Err(PyValueError::new_err(format!(
            "The transformation matrix must have 12 values.  It had {}.",
            matrix.len()
        )));
    }
    if finalTrans.len() != 3 {
        return Err(PyValueError::new_err(format!(
            "The final translation vector must have 3 values.  It had {}.",
            finalTrans.len()
        )));
    }
    if finalRot.len() != 9 {
        return Err(PyValueError::new_err(format!(
            "The final rotation vector must have 9 values.  It had {}.",
            finalRot.len()
        )));
    }
    probeConformer.with_mut(py, |c| {
        transform_conformer_by_overlay(&finalTrans, &finalRot, &matrix, &probeShape.shape, c)
    })
}

/// the shape of a conformer of the molecule (centred and in its inertial frame unless the options say otherwise)
#[pyfunction]
#[pyo3(name = "PrepareConformer", signature = (mol, confId = -1, opts = None))]
#[allow(non_snake_case)]
fn prepare_conformer_py(
    mol: PyRef<'_, PyMol>,
    confId: i32,
    opts: Option<PyRef<'_, PyShapeInputOptions>>,
) -> PyResult<PyShapeInput> {
    let o = opts.map(|o| o.opts.clone()).unwrap_or_default();
    Ok(PyShapeInput {
        shape: prepare_conformer(&mol.mol, confId, &o).py()?,
    })
}

/// Score a molecule against a molecule or a shape, without any alignment. Returns (shape_score, color_score).
#[pyfunction]
#[pyo3(name = "ScoreMol", signature = (a, b, *args, **kwargs))]
fn score_mol(
    a: &Bound<'_, PyAny>,
    b: PyRef<'_, PyMol>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(f64, f64)> {
    if let Ok(shape) = a.extract::<PyRef<'_, PyShapeInput>>() {
        let o: PyRef<'_, PyShapeInputOptions> =
            arg_required(args, kwargs, 0, "molShapeOpts")?.extract()?;
        let c: i32 = arg(args, kwargs, 1, "molConfId", -1)?;
        return score_molecule_to_shape(&shape.shape, &b.mol, &o.opts, c).py();
    }
    let m1: PyRef<'_, PyMol> = a.extract()?;
    let o1: PyRef<'_, PyShapeInputOptions> =
        arg_required(args, kwargs, 0, "mol1ShapeOpts")?.extract()?;
    let o2: PyRef<'_, PyShapeInputOptions> =
        arg_required(args, kwargs, 1, "mol2ShapeOpts")?.extract()?;
    let c1: i32 = arg(args, kwargs, 2, "mol1ConfId", -1)?;
    let c2: i32 = arg(args, kwargs, 3, "mol2ConfId", -1)?;
    score_molecule(&m1.mol, &b.mol, &o1.opts, &o2.opts, c1, c2).py()
}

/// Score two shapes against each other, without any alignment. Returns (shape_score, color_score).
#[pyfunction]
#[pyo3(name = "ScoreShape", signature = (shape1, shape2, useColors = true))]
#[allow(non_snake_case)]
fn score_shape_py(
    shape1: PyRef<'_, PyShapeInput>,
    shape2: PyRef<'_, PyShapeInput>,
    useColors: bool,
) -> (f64, f64) {
    score_shape(&shape1.shape, &shape2.shape, useColors)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing the PubChem shape alignment functionality.",
    )?;
    m.add_class::<PyShapeInputOptions>()?;
    m.add_class::<PyShapeInput>()?;
    m.add_function(wrap_pyfunction!(align_mol, m)?)?;
    m.add_function(wrap_pyfunction!(align_shapes, m)?)?;
    m.add_function(wrap_pyfunction!(transform_conformer_py, m)?)?;
    m.add_function(wrap_pyfunction!(prepare_conformer_py, m)?)?;
    m.add_function(wrap_pyfunction!(score_mol, m)?)?;
    m.add_function(wrap_pyfunction!(score_shape_py, m)?)?;
    Ok(())
}
