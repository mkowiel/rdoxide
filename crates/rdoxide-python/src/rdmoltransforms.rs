//! `rdoxide.Chem.rdMolTransforms` (upstream
//! `MolTransforms/Wrap/rdMolTransforms.cpp`): 3D operations on conformers.

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::geometry::{Point3D, Transform3D};
use rdoxide::graphmol::Conformer;
use rdoxide::graphmol::mol_transforms as mt;
use rdoxide::graphmol::romol::ROMol;

use crate::conformer::{ConfHandle, PyConformer};
use crate::errors::{Check, IntoPyResult, Site, check_failed};
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

const NO_OWNER: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/Conformer.h",
    line: 92,
    expr: "dp_mol",
};

/// Runs `f` on the molecule that owns the conformer and the conformer's id;
/// conformers outside a molecule fail as upstream's `getOwningMol()` does.
fn with_owner<R>(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    f: impl FnOnce(&mut ROMol, i32) -> PyResult<R>,
) -> PyResult<R> {
    match &conf.borrow().handle {
        ConfHandle::InMol { mol, id } => {
            let mut m = mol.bind(py).borrow_mut();
            f(&mut m.mol, *id as i32)
        }
        ConfHandle::Owned(_) => Err(check_failed(NO_OWNER, "no owner")),
    }
}

/// Runs `f` on the conformer's coordinates alone.
fn with_conf<R>(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    f: impl FnOnce(&Conformer) -> PyResult<R>,
) -> PyResult<R> {
    conf.borrow().with(py, f)?
}

fn ndarray<'py>(py: Python<'py>, data: Vec<f64>, shape: &[usize]) -> PyResult<Bound<'py, PyAny>> {
    py.import("numpy")?
        .call_method1("array", (data,))?
        .call_method1("reshape", (shape.to_vec(),))
}

fn center_arg(center: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Point3D>> {
    match center {
        None => Ok(None),
        Some(c) if c.is_none() => Ok(None),
        Some(c) => Ok(Some(c.extract::<PyRef<'_, PyPoint3D>>()?.p)),
    }
}

/// Compute the centroid of the conformation - hydrogens are ignored and no attention
/// is paid to the difference in sizes of the heavy atoms; however,
/// an optional vector of weights can be passed.
#[pyfunction]
#[pyo3(name = "ComputeCentroid", signature = (conf, ignoreHs = None, weights = None))]
#[allow(non_snake_case)]
fn compute_centroid(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    ignoreHs: Option<&Bound<'_, PyAny>>,
    weights: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyPoint3D> {
    if weights.is_some_and(|w| !w.is_none()) {
        return Err(PyTypeError::new_err(
            "ComputeCentroid(): the weights argument only accepts None",
        ));
    }
    let ignore_hs = bool_arg(ignoreHs, true)?;
    let c = conf.clone();
    with_owner(py, &c, |mol, id| {
        let p = mt::compute_centroid(mol, mol.conformer(id).py()?, ignore_hs, None).py()?;
        Ok(PyPoint3D { p })
    })
}

/// Compute the transformation required aligna conformer so that
/// the principal axes align up with the x,y, z axes
/// The conformer itself is left unchanged
///   ARGUMENTS:
///     - conf : the conformer of interest
///     - center : optional center point to compute the principal axes around (defaults to the centroid)
///     - normalizeCovar : optionally normalize the covariance matrix by the number of atoms
#[pyfunction]
#[pyo3(name = "ComputeCanonicalTransform", signature = (conf, center = None, normalizeCovar = None, ignoreHs = None))]
#[allow(non_snake_case)]
fn compute_canonical_transform<'py>(
    py: Python<'py>,
    conf: &Bound<'py, PyConformer>,
    center: Option<&Bound<'py, PyAny>>,
    normalizeCovar: Option<&Bound<'py, PyAny>>,
    ignoreHs: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    let center = center_arg(center)?;
    let normalize = bool_arg(normalizeCovar, false)?;
    let ignore_hs = bool_arg(ignoreHs, true)?;
    let trans = with_owner(py, conf, |mol, id| {
        let (trans, _) = mt::compute_canonical_transform(
            mol,
            mol.conformer(id).py()?,
            center,
            normalize,
            ignore_hs,
        )
        .py()?;
        Ok(trans)
    })?;
    ndarray(py, trans.data.to_vec(), &[4, 4])
}

/// The weights of the principal-axes functions: a sequence of numbers as long
/// as the conformer has atoms.
fn weights_arg(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    weights: Option<&Bound<'_, PyAny>>,
) -> PyResult<Option<Vec<f64>>> {
    let Some(w) = weights.filter(|w| !w.is_none()) else {
        return Ok(None);
    };
    let n = w.len()?;
    if n != conf.borrow().with(py, Conformer::num_atoms)? as usize {
        return Err(PyValueError::new_err(
            "The Python container must have length equal to conf.GetNumAtoms()",
        ));
    }
    (0..n)
        .map(|i| w.get_item(i)?.extract::<f64>())
        .collect::<PyResult<Vec<_>>>()
        .map(Some)
}

type PrincipalAxesFn =
    fn(&mut ROMol, i32, bool, bool, Option<&[f64]>) -> rdoxide::Result<Option<mt::PrincipalAxes>>;

fn principal_axes<'py>(
    py: Python<'py>,
    conf: &Bound<'py, PyConformer>,
    ignore_hs: Option<&Bound<'py, PyAny>>,
    weights: Option<&Bound<'py, PyAny>>,
    compute: PrincipalAxesFn,
) -> PyResult<Bound<'py, PyTuple>> {
    let ignore_hs = bool_arg(ignore_hs, true)?;
    let weights = weights_arg(py, conf, weights)?;
    let res = with_owner(py, conf, |mol, id| {
        compute(mol, id, ignore_hs, true, weights.as_deref()).py()
    })?;
    match res {
        Some(r) => {
            let axes = ndarray(py, r.axes.iter().flatten().copied().collect(), &[3, 3])?;
            let moments = ndarray(py, r.moments.to_vec(), &[3])?;
            PyTuple::new(py, [axes, moments])
        }
        None => PyTuple::new(py, [py.None(), py.None()]),
    }
}

/// Compute principal axes and moments of inertia for a conformer
/// These values are calculated from the inertia tensor:
/// Iij = - sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
/// Iii = sum_{s=1..N} sum_{j!=i} (w_s * r_{sj} * r_{sj})
/// where the coordinates are relative to the center of mass.
///
///   ARGUMENTS:
///     - conf : the conformer of interest
///     - ignoreHs : if True, ignore hydrogen atoms
///     - weights : if present, used to weight the atomic coordinates
///
///   Returns a (principal axes, principal moments) tuple
#[pyfunction]
#[pyo3(name = "ComputePrincipalAxesAndMoments", signature = (conf, ignoreHs = None, weights = None))]
#[allow(non_snake_case)]
fn compute_principal_axes_and_moments<'py>(
    py: Python<'py>,
    conf: &Bound<'py, PyConformer>,
    ignoreHs: Option<&Bound<'py, PyAny>>,
    weights: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    principal_axes(
        py,
        conf,
        ignoreHs,
        weights,
        mt::compute_principal_axes_and_moments,
    )
}

/// Compute principal axes and moments from the gyration matrix of a conformer
/// These values are calculated from the gyration matrix/tensor:
/// Iij = sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
/// Iii = sum_{s=1..N} sum_{t!=s}(w_s * r_{si} * r_{ti})
/// where the coordinates are relative to the center of mass.
///
///   ARGUMENTS:
///     - conf : the conformer of interest
///     - ignoreHs : if True, ignore hydrogen atoms
///     - weights : if present, used to weight the atomic coordinates
///
///   Returns a (principal axes, principal moments) tuple
#[pyfunction]
#[pyo3(
    name = "ComputePrincipalAxesAndMomentsFromGyrationMatrix",
    signature = (conf, ignoreHs = None, weights = None)
)]
#[allow(non_snake_case)]
fn compute_principal_axes_and_moments_from_gyration_matrix<'py>(
    py: Python<'py>,
    conf: &Bound<'py, PyConformer>,
    ignoreHs: Option<&Bound<'py, PyAny>>,
    weights: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    principal_axes(
        py,
        conf,
        ignoreHs,
        weights,
        mt::compute_principal_axes_and_moments_from_gyration_matrix,
    )
}

/// Transform the coordinates of a conformer
#[pyfunction]
#[pyo3(name = "TransformConformer", signature = (conf, trans))]
fn transform_conformer(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    trans: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let np = py.import("numpy")?;
    if !trans.is_instance(&np.getattr("ndarray")?)? {
        return Err(PyValueError::new_err(
            "Expecting a numeric array for transformation",
        ));
    }
    let n_rows: usize = trans.getattr("shape")?.get_item(0)?.extract()?;
    let flat: Vec<f64> = np
        .call_method1("ascontiguousarray", (trans, "float64"))?
        .call_method0("ravel")?
        .extract()?;
    // as upstream, an n x n array overwrites the first n*n entries of the matrix
    let mut transform = Transform3D::identity();
    let count = (n_rows * n_rows).min(16).min(flat.len());
    transform.data[..count].copy_from_slice(&flat[..count]);
    conf.borrow_mut()
        .with_mut(py, |c| mt::transform_conformer(c, &transform))
}

/// Canonicalize the orientation of a conformer so that its principal axes
/// around the specified center point coincide with the x, y, z axes
///
///   ARGUMENTS:
///     - conf : conformer of interest
///     - center : optionally center point about which the principal axes are computed
///                if not specified the centroid of the conformer will be used
///     - normalizeCovar : Optionally normalize the covariance matrix by the number of atoms
#[pyfunction]
#[pyo3(name = "CanonicalizeConformer", signature = (conf, center = None, normalizeCovar = None, ignoreHs = None))]
#[allow(non_snake_case)]
fn canonicalize_conformer(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    center: Option<&Bound<'_, PyAny>>,
    normalizeCovar: Option<&Bound<'_, PyAny>>,
    ignoreHs: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let center = center_arg(center)?;
    let normalize = bool_arg(normalizeCovar, false)?;
    let ignore_hs = bool_arg(ignoreHs, true)?;
    with_owner(py, conf, |mol, id| {
        mt::canonicalize_conformer(mol, id, center, normalize, ignore_hs).py()
    })
}

/// Loop over the conformers in a molecule and canonicalize their orientation
#[pyfunction]
#[pyo3(name = "CanonicalizeMol", signature = (mol, normalizeCovar = None, ignoreHs = None))]
#[allow(non_snake_case)]
fn canonicalize_mol(
    mol: &Bound<'_, PyMol>,
    normalizeCovar: Option<&Bound<'_, PyAny>>,
    ignoreHs: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let normalize = bool_arg(normalizeCovar, false)?;
    let ignore_hs = bool_arg(ignoreHs, true)?;
    mt::canonicalize_mol(&mut mol.borrow_mut().mol, normalize, ignore_hs).py()
}

/// Returns the bond length in angstrom between atoms i, j
#[pyfunction]
#[pyo3(name = "GetBondLength", signature = (conf, iAtomId, jAtomId))]
#[allow(non_snake_case)]
fn get_bond_length(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
) -> PyResult<f64> {
    with_conf(py, conf, |c| mt::get_bond_length(c, iAtomId, jAtomId).py())
}

/// Sets the bond length in angstrom between atoms i, j; all atoms bonded to atom j are moved
#[pyfunction]
#[pyo3(name = "SetBondLength", signature = (conf, iAtomId, jAtomId, value))]
#[allow(non_snake_case)]
fn set_bond_length(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    value: f64,
) -> PyResult<()> {
    with_owner(py, conf, |mol, id| {
        mt::set_bond_length(mol, id, iAtomId, jAtomId, value).py()
    })
}

/// Returns the angle in radians between atoms i, j, k
#[pyfunction]
#[pyo3(name = "GetAngleRad", signature = (conf, iAtomId, jAtomId, kAtomId))]
#[allow(non_snake_case)]
fn get_angle_rad(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
) -> PyResult<f64> {
    with_conf(py, conf, |c| {
        mt::get_angle_rad(c, iAtomId, jAtomId, kAtomId).py()
    })
}

/// Returns the angle in degrees between atoms i, j, k
#[pyfunction]
#[pyo3(name = "GetAngleDeg", signature = (conf, iAtomId, jAtomId, kAtomId))]
#[allow(non_snake_case)]
fn get_angle_deg(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
) -> PyResult<f64> {
    with_conf(py, conf, |c| {
        mt::get_angle_deg(c, iAtomId, jAtomId, kAtomId).py()
    })
}

/// Sets the angle in radians between atoms i, j, k; all atoms bonded to atom k are moved
#[pyfunction]
#[pyo3(name = "SetAngleRad", signature = (conf, iAtomId, jAtomId, kAtomId, value))]
#[allow(non_snake_case)]
fn set_angle_rad(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
    value: f64,
) -> PyResult<()> {
    with_owner(py, conf, |mol, id| {
        mt::set_angle_rad(mol, id, iAtomId, jAtomId, kAtomId, value).py()
    })
}

/// Sets the angle in degrees between atoms i, j, k; all atoms bonded to atom k are moved
#[pyfunction]
#[pyo3(name = "SetAngleDeg", signature = (conf, iAtomId, jAtomId, kAtomId, value))]
#[allow(non_snake_case)]
fn set_angle_deg(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
    value: f64,
) -> PyResult<()> {
    with_owner(py, conf, |mol, id| {
        mt::set_angle_deg(mol, id, iAtomId, jAtomId, kAtomId, value).py()
    })
}

/// Returns the dihedral angle in radians between atoms i, j, k, l
#[pyfunction]
#[pyo3(name = "GetDihedralRad", signature = (conf, iAtomId, jAtomId, kAtomId, lAtomId))]
#[allow(non_snake_case)]
fn get_dihedral_rad(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
    lAtomId: u32,
) -> PyResult<f64> {
    with_conf(py, conf, |c| {
        mt::get_dihedral_rad(c, iAtomId, jAtomId, kAtomId, lAtomId).py()
    })
}

/// Returns the dihedral angle in degrees between atoms i, j, k, l
#[pyfunction]
#[pyo3(name = "GetDihedralDeg", signature = (conf, iAtomId, jAtomId, kAtomId, lAtomId))]
#[allow(non_snake_case)]
fn get_dihedral_deg(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
    lAtomId: u32,
) -> PyResult<f64> {
    with_conf(py, conf, |c| {
        mt::get_dihedral_deg(c, iAtomId, jAtomId, kAtomId, lAtomId).py()
    })
}

/// Sets the dihedral angle in radians between atoms i, j, k, l; all atoms bonded to atom l are moved
#[pyfunction]
#[pyo3(name = "SetDihedralRad", signature = (conf, iAtomId, jAtomId, kAtomId, lAtomId, value))]
#[allow(non_snake_case)]
fn set_dihedral_rad(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
    lAtomId: u32,
    value: f64,
) -> PyResult<()> {
    with_owner(py, conf, |mol, id| {
        mt::set_dihedral_rad(mol, id, iAtomId, jAtomId, kAtomId, lAtomId, value).py()
    })
}

/// Sets the dihedral angle in degrees between atoms i, j, k, l; all atoms bonded to atom l are moved
#[pyfunction]
#[pyo3(name = "SetDihedralDeg", signature = (conf, iAtomId, jAtomId, kAtomId, lAtomId, value))]
#[allow(non_snake_case)]
fn set_dihedral_deg(
    py: Python<'_>,
    conf: &Bound<'_, PyConformer>,
    iAtomId: u32,
    jAtomId: u32,
    kAtomId: u32,
    lAtomId: u32,
    value: f64,
) -> PyResult<()> {
    with_owner(py, conf, |mol, id| {
        mt::set_dihedral_deg(mol, id, iAtomId, jAtomId, kAtomId, lAtomId, value).py()
    })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to perform 3D operations like rotate and translate conformations",
    )?;
    m.add_function(wrap_pyfunction!(compute_centroid, m)?)?;
    m.add_function(wrap_pyfunction!(compute_canonical_transform, m)?)?;
    m.add_function(wrap_pyfunction!(compute_principal_axes_and_moments, m)?)?;
    m.add_function(wrap_pyfunction!(
        compute_principal_axes_and_moments_from_gyration_matrix,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(transform_conformer, m)?)?;
    m.add_function(wrap_pyfunction!(canonicalize_conformer, m)?)?;
    m.add_function(wrap_pyfunction!(canonicalize_mol, m)?)?;
    m.add_function(wrap_pyfunction!(get_bond_length, m)?)?;
    m.add_function(wrap_pyfunction!(set_bond_length, m)?)?;
    m.add_function(wrap_pyfunction!(get_angle_rad, m)?)?;
    m.add_function(wrap_pyfunction!(get_angle_deg, m)?)?;
    m.add_function(wrap_pyfunction!(set_angle_rad, m)?)?;
    m.add_function(wrap_pyfunction!(set_angle_deg, m)?)?;
    m.add_function(wrap_pyfunction!(get_dihedral_rad, m)?)?;
    m.add_function(wrap_pyfunction!(get_dihedral_deg, m)?)?;
    m.add_function(wrap_pyfunction!(set_dihedral_rad, m)?)?;
    m.add_function(wrap_pyfunction!(set_dihedral_deg, m)?)?;
    Ok(())
}
