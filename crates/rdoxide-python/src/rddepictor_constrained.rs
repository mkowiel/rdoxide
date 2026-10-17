//! Constrained and 3D-mimicking depictions of `rdoxide.Chem.rdDepictor`
//! (upstream `Depictor/Wrap/rdDepictor.cpp`).

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::depictor::constrained::{
    ConstrainedDepictionParams, generate_depiction_matching_2d_structure,
    generate_depiction_matching_2d_structure_with_map, generate_depiction_matching_3d_structure,
};
use rdoxide::graphmol::substruct::MatchVect;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Parameters controlling constrained depiction
#[pyclass(
    name = "ConstrainedDepictionParams",
    module = "rdoxide.Chem.rdDepictor",
    from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyConstrainedDepictionParams {
    /// if False (default), a DepictException is thrown if the molecule does
    /// not have a substructure match to the reference; if True, an
    /// unconstrained depiction will be generated
    #[pyo3(get, set, name = "acceptFailure")]
    accept_failure: bool,
    /// if True, use RDKit to generate coordinates even if preferCoordGen is
    /// set to True; defaults to False
    #[pyo3(get, set, name = "forceRDKit")]
    force_rdkit: bool,
    /// if True, terminal dummy atoms in the reference are ignored if they
    /// match an implicit hydrogen in the molecule or if they are attached
    /// top a query atom; defaults to False
    #[pyo3(get, set, name = "allowRGroups")]
    allow_rgroups: bool,
    /// if False (default), a part of the molecule is hard-constrained to
    /// have the same coordinates as the reference, and the rest of the
    /// molecule is built around it; if True, coordinates from conformation
    /// existingConfId are preserved (if they exist) or generated without
    /// constraints (if they do not exist), then the conformation is
    /// rigid-body aligned to the reference
    #[pyo3(get, set, name = "alignOnly")]
    align_only: bool,
    /// if True (default), existing wedging information will be updated or
    /// cleared as required; if False, existing molblock wedging information
    /// will always be preserved
    #[pyo3(get, set, name = "adjustMolBlockWedging")]
    adjust_mol_block_wedging: bool,
    /// conformation id whose 2D coordinates should be rigid-body aligned to
    /// the reference (if alignOnly is True), or used to determine whether
    /// existing molblock wedging information can be preserved following the
    /// constrained depiction (if adjustMolBlockWedging is True
    #[pyo3(get, set, name = "existingConfId")]
    existing_conf_id: i32,
    /// use templates to generate coordinates of complex ring systems
    #[pyo3(get, set, name = "useRingTemplates")]
    use_ring_templates: bool,
}

#[pymethods]
impl PyConstrainedDepictionParams {
    #[new]
    fn new() -> Self {
        let d = ConstrainedDepictionParams::default();
        PyConstrainedDepictionParams {
            accept_failure: d.accept_failure,
            force_rdkit: d.force_rdkit,
            allow_rgroups: d.allow_rgroups,
            align_only: d.align_only,
            adjust_mol_block_wedging: d.adjust_mol_block_wedging,
            existing_conf_id: d.existing_conf_id,
            use_ring_templates: d.use_ring_templates,
        }
    }
}

impl From<&PyConstrainedDepictionParams> for ConstrainedDepictionParams {
    fn from(p: &PyConstrainedDepictionParams) -> Self {
        ConstrainedDepictionParams {
            accept_failure: p.accept_failure,
            force_rdkit: p.force_rdkit,
            allow_rgroups: p.allow_rgroups,
            align_only: p.align_only,
            adjust_mol_block_wedging: p.adjust_mol_block_wedging,
            existing_conf_id: p.existing_conf_id,
            use_ring_templates: p.use_ring_templates,
        }
    }
}

/// Binds positional and keyword arguments to `names` the way a Python
/// signature would (`None` for the ones not given).
fn bind_args<'py>(
    func: &str,
    names: &[&str],
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Vec<Option<Bound<'py, PyAny>>>> {
    if args.len() > names.len() {
        return Err(PyTypeError::new_err(format!(
            "{func}() takes at most {} arguments ({} given)",
            names.len() + 2,
            args.len() + 2
        )));
    }
    let mut bound: Vec<Option<Bound<'py, PyAny>>> = vec![None; names.len()];
    for (slot, arg) in bound.iter_mut().zip(args.iter()) {
        *slot = Some(arg);
    }
    if let Some(kw) = kwargs {
        for (k, v) in kw.iter() {
            let key: String = k.extract()?;
            let idx = names.iter().position(|n| *n == key).ok_or_else(|| {
                PyTypeError::new_err(format!(
                    "{func}() got an unexpected keyword argument '{key}'"
                ))
            })?;
            if bound[idx].is_some() {
                return Err(PyTypeError::new_err(format!(
                    "{func}() got multiple values for argument '{key}'"
                )));
            }
            bound[idx] = Some(v);
        }
    }
    Ok(bound)
}

fn int_or(arg: &Option<Bound<'_, PyAny>>, default: i32) -> PyResult<i32> {
    arg.as_ref().map_or(Ok(default), |v| v.extract())
}

fn bool_or(arg: &Option<Bound<'_, PyAny>>, default: bool) -> PyResult<bool> {
    arg.as_ref().map_or(Ok(default), |v| v.is_truthy())
}

fn params_from(arg: &Option<Bound<'_, PyAny>>) -> PyResult<ConstrainedDepictionParams> {
    match arg {
        Some(v) if !v.is_none() => {
            Ok((&*v.extract::<PyRef<'_, PyConstrainedDepictionParams>>()?).into())
        }
        _ => Ok(ConstrainedDepictionParams::default()),
    }
}

fn atom_map_from(arg: &Bound<'_, PyAny>) -> PyResult<MatchVect> {
    let mut res = MatchVect::new();
    for item in arg.try_iter()? {
        let pair: (i32, i32) = item?
            .extract()
            .map_err(|_| PyTypeError::new_err("atomMap must be a sequence of 2-tuples"))?;
        res.push(pair);
    }
    Ok(res)
}

/// Generate a depiction for a molecule where a piece of the molecule is
/// constrained to have the same coordinates as a reference.
///
/// Called as `(mol, reference, confId=-1, refPatt=None, params=None)` or
/// `(mol, reference, confId=-1, refPatt=None, acceptFailure=False,
/// forceRDKit=False, allowRGroups=False)`; the substructure match between
/// the molecule and the reference is found automatically, and a tuple of
/// (refIdx, molIdx) tuples is returned. With an explicit atom mapping,
/// `(mol, reference, atomMap, confId=-1, params=None)` or
/// `(mol, reference, atomMap, confId, forceRDKit)`, nothing is returned.
///
/// The constraint can be hard (default) or soft
/// (`ConstrainedDepictionParams.alignOnly=True`: existing coordinates are
/// kept, or generated without constraints, then rigid-body aligned to the
/// reference). A DepictException is raised when the molecule has no
/// substructure match to the reference and
/// `ConstrainedDepictionParams.acceptFailure` is False.
#[pyfunction]
#[pyo3(name = "GenerateDepictionMatching2DStructure", signature = (mol, reference, *args, **kwargs))]
fn generate_depiction_matching_2d_structure_py<'py>(
    py: Python<'py>,
    mut mol: PyRefMut<'py, PyMol>,
    reference: PyRef<'py, PyMol>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyAny>> {
    const NAME: &str = "GenerateDepictionMatching2DStructure";
    let with_atom_map = kwargs.is_some_and(|k| k.contains("atomMap").unwrap_or(false))
        || args.get_item(0).is_ok_and(|a| a.extract::<i32>().is_err());
    if with_atom_map {
        let b = bind_args(NAME, &["atomMap", "confId", "params"], args, kwargs)?;
        let atom_map = b[0]
            .as_ref()
            .ok_or_else(|| PyTypeError::new_err("atomMap is required"))
            .and_then(atom_map_from)?;
        let conf_id: i32 = int_or(&b[1], -1)?;
        let params = match &b[2] {
            Some(v) if v.is_instance_of::<pyo3::types::PyBool>() => ConstrainedDepictionParams {
                force_rdkit: v.extract()?,
                ..Default::default()
            },
            other => params_from(other)?,
        };
        let (m, r) = (&mut mol.mol, &reference.mol);
        py.detach(|| {
            generate_depiction_matching_2d_structure_with_map(m, r, &atom_map, conf_id, &params)
        })
        .py()?;
        return Ok(py.None().into_bound(py));
    }
    let b = bind_args(
        NAME,
        &[
            "confId",
            "refPatt",
            "acceptFailure",
            "forceRDKit",
            "allowRGroups",
            "params",
        ],
        args,
        kwargs,
    )?;
    let conf_id: i32 = int_or(&b[0], -1)?;
    let pattern: Option<PyRef<'py, PyMol>> = match &b[1] {
        Some(v) if !v.is_none() => Some(v.extract()?),
        _ => None,
    };
    // the fifth positional argument is either `acceptFailure` or `params`
    let params = match (&b[2], &b[5]) {
        (Some(v), None)
            if v.extract::<PyRef<'_, PyConstrainedDepictionParams>>()
                .is_ok() =>
        {
            params_from(&b[2])?
        }
        (_, Some(_)) => params_from(&b[5])?,
        _ => ConstrainedDepictionParams {
            accept_failure: bool_or(&b[2], false)?,
            force_rdkit: bool_or(&b[3], false)?,
            allow_rgroups: bool_or(&b[4], false)?,
            ..Default::default()
        },
    };
    let (m, r) = (&mut mol.mol, &reference.mol);
    let p = pattern.as_ref().map(|p| &p.mol);
    let res = py
        .detach(|| generate_depiction_matching_2d_structure(m, r, conf_id, p, &params))
        .py()?;
    Ok(PyTuple::new(py, res)?.into_any())
}

/// Generate a depiction for a molecule where a piece of the molecule is
/// constrained to have coordinates similar to those of a 3D reference
/// structure.
///
///   ARGUMENTS:
///
///     - mol: the molecule to be aligned, this will come back with a single
///       conformer containing the 2D coordinates.
///     - reference: a molecule with the reference atoms to align to. By
///       default this should be the same as mol, but with 3D coordinates
///     - confId: (optional) the id of the reference conformation to use
///     - referencePattern: (optional) a query molecule to map a subset of
///       the reference onto the mol, so that only some of the atoms are
///       aligned.
///     - acceptFailure: (optional) if True, standard depictions will be
///       generated for molecules that don't match the reference or the
///       referencePattern; if False, throws a DepictException.
///     - forceRDKit: (optional) use RDKit to generate coordinates even if
///       preferCoordGen is set to true
#[pyfunction]
#[pyo3(name = "GenerateDepictionMatching3DStructure", signature = (mol, reference, confId = -1, refPatt = None, acceptFailure = false, forceRDKit = false))]
#[allow(non_snake_case)]
fn generate_depiction_matching_3d_structure_py(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    reference: PyRef<'_, PyMol>,
    confId: i32,
    refPatt: Option<PyRef<'_, PyMol>>,
    acceptFailure: bool,
    forceRDKit: bool,
) -> PyResult<()> {
    let (m, r) = (&mut mol.mol, &reference.mol);
    let p = refPatt.as_ref().map(|p| &p.mol);
    py.detach(|| {
        generate_depiction_matching_3d_structure(m, r, confId, p, acceptFailure, forceRDKit)
    })
    .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyConstrainedDepictionParams>()?;
    m.add_function(wrap_pyfunction!(
        generate_depiction_matching_2d_structure_py,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        generate_depiction_matching_3d_structure_py,
        m
    )?)?;
    Ok(())
}
