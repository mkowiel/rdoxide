//! `rdoxide.Chem.rdDepictor` (upstream `Depictor/Wrap/rdDepictor.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rdoxide::graphmol::depictor::depict_utils::{IntPoint2DMap, bond_len, set_bond_len};
use rdoxide::graphmol::depictor::templates::{
    add_ring_system_templates, load_default_ring_system_templates, set_ring_system_templates,
};
use rdoxide::graphmol::depictor::{
    Compute2DCoordParameters, compute_2d_coords, compute_2d_coords_mimic_dist_mat,
    normalize_depiction, straighten_depiction,
};

use crate::errors::IntoPyResult;
use crate::geometry::PyPoint2D;
use crate::mol::PyMol;

/// Compute 2D coordinates for a molecule.
///   The resulting coordinates are stored on each atom of the molecule
///
///   ARGUMENTS:
///
///      mol - the molecule of interest
///      canonOrient - orient the molecule in a canonical way
///      clearConfs - if true, all existing conformations on the molecule
///              will be cleared
///      coordMap - a dictionary mapping atom Ids -> Point2D objects
///                 with starting coordinates for atoms that should
///                 have their positions locked.
///      nFlipsPerSample - number of rotatable bonds that are
///                 flipped at random at a time.
///      nSample - Number of random samplings of rotatable bonds.
///      sampleSeed - seed for the random sampling process.
///      permuteDeg4Nodes - allow permutation of bonds at a degree 4
///                 node during the sampling process
///      bondLength - change the default bond length for depiction
///      forceRDKit - use RDKit to generate coordinates even if
///                   preferCoordGen is set to true
///      useRingTemplates - use templates to generate coordinates of complex
///                   ring systems
///      usePathAngleExpansion - use path angle expansion to resolve collisions
///
///   RETURNS:
///
///      ID of the conformation added to the molecule
#[pyfunction]
#[pyo3(name = "Compute2DCoords", signature = (mol, canonOrient = true, clearConfs = true, coordMap = None, nFlipsPerSample = 0, nSample = 0, sampleSeed = 0, permuteDeg4Nodes = 0, bondLength = -1.0, forceRDKit = false, useRingTemplates = false, usePathAngleExpansion = true))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn compute_2d_coords_py(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    canonOrient: bool,
    clearConfs: bool,
    coordMap: Option<&Bound<'_, PyDict>>,
    nFlipsPerSample: u32,
    nSample: u32,
    sampleSeed: i32,
    // an int, as boost::python also accepts 0/1 for a bool
    permuteDeg4Nodes: i64,
    bondLength: f64,
    forceRDKit: bool,
    useRingTemplates: bool,
    usePathAngleExpansion: bool,
) -> PyResult<u32> {
    let mut cmap = IntPoint2DMap::new();
    if let Some(d) = coordMap {
        for (k, v) in d.iter() {
            let id: u32 = k.extract()?;
            if id >= mol.mol.num_atoms() {
                return Err(PyValueError::new_err("atom index out of range"));
            }
            let p: PyRef<'_, PyPoint2D> = v.extract()?;
            cmap.insert(id, p.p);
        }
    }
    let old = bond_len();
    if bondLength > 0.0 {
        set_bond_len(bondLength);
    }
    let params = Compute2DCoordParameters {
        coord_map: if cmap.is_empty() { None } else { Some(&cmap) },
        canon_orient: canonOrient,
        clear_confs: clearConfs,
        n_flips_per_sample: nFlipsPerSample,
        n_samples: nSample,
        sample_seed: sampleSeed,
        permute_deg4_nodes: permuteDeg4Nodes != 0,
        force_rdkit: forceRDKit,
        use_ring_templates: useRingTemplates,
        use_path_angle_expansion: usePathAngleExpansion,
    };
    let m = &mut mol.mol;
    let res = py.detach(|| compute_2d_coords(m, &params));
    if bondLength > 0.0 {
        set_bond_len(old);
    }
    res.py()
}

/// Compute 2D coordinates for a molecule such that the inter-atom
/// distances mimic those in a distance matrix.
///
/// This is done by randomly sampling (flipping) the rotatable bonds in the
/// molecule and evaluating a cost function which contains two components.
/// The first component is the sum of the inverse of the squared inter-atom
/// distances, this helps in spreading the atoms far from each other. The
/// second component is the sum of squares of the difference in distance
/// between those in distMat and the generated structure. The user can
/// adjust the relative importance of the two components via an adjustable
/// parameter (see below).
///
///   ARGUMENTS:
///
///     - mol : the molecule of interest
///     - distMat : distance matrix that we want the 2D structure to mimic
///     - canonOrient : canonicalize the orientation so that the long axes
///       align with the x-axis etc.
///     - clearConfs : clear all existing conformations on the molecule
///     - weightDistMat : weight assigned in the cost function to mimicking
///       the distance matrix. This must be between (0.0,1.0). (1.0-weightDistMat)
///       is then the weight assigned to improving the density of the 2D
///       structure i.e. try to make it spread out
///     - nFlipsPerSample : number of rotatable bonds that are flipped at
///       random at a time.
///     - nSample : Number of random samplings of rotatable bonds.
///     - sampleSeed : seed for the random sampling process.
///     - permuteDeg4Nodes : allow permutation of bonds at a degree 4
///       node during the sampling process
///     - bondLength : change the default bond length for depiction
///     - forceRDKit : use RDKit to generate coordinates even if
///       preferCoordGen is set to true
///
///   RETURNS:
///
///     ID of the conformation added to the molecule
#[pyfunction]
#[pyo3(name = "Compute2DCoordsMimicDistmat", signature = (mol, distMat, canonOrient = false, clearConfs = true, weightDistMat = 0.5, nFlipsPerSample = 3, nSample = 100, sampleSeed = 100, permuteDeg4Nodes = 1, bondLength = -1.0, forceRDKit = false))]
#[allow(non_snake_case, clippy::too_many_arguments, unused_variables)]
fn compute_2d_coords_mimic_distmat_py(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    distMat: &Bound<'_, PyAny>,
    canonOrient: bool,
    clearConfs: bool,
    weightDistMat: f64,
    nFlipsPerSample: u32,
    nSample: u32,
    sampleSeed: i32,
    // an int, as boost::python also accepts 0/1 for a bool
    permuteDeg4Nodes: i64,
    bondLength: f64,
    forceRDKit: bool,
) -> PyResult<u32> {
    let dmat: Vec<f64> = distMat
        .extract()
        .map_err(|_| PyValueError::new_err("Argument isn't an array"))?;
    let na = mol.mol.num_atoms() as usize;
    if dmat.len() != na * na.saturating_sub(1) / 2 {
        return Err(PyValueError::new_err(
            "The array size does not match the number of atoms in the molecule",
        ));
    }
    let old = bond_len();
    if bondLength > 0.0 {
        set_bond_len(bondLength);
    }
    let m = &mut mol.mol;
    let res = py.detach(|| {
        compute_2d_coords_mimic_dist_mat(
            m,
            &dmat,
            canonOrient,
            clearConfs,
            weightDistMat,
            nFlipsPerSample,
            nSample,
            sampleSeed,
            permuteDeg4Nodes != 0,
        )
    });
    if bondLength > 0.0 {
        set_bond_len(old);
    }
    res.py()
}

/// Rotates the 2D depiction such that the majority of bonds have a
/// 30-degree angle with the X axis.
///
///   ARGUMENTS:
///
///     - mol: the molecule to be rotated
///     - confId: (optional) the id of the reference conformation to use
///     - minimizeRotation: (optional) if False (the default), the
///       molecule is rotated such that the majority of bonds have an angle
///       with the X axis of 30 or 90 degrees. If True, the minimum rotation
///       is applied such that the majority of bonds have an angle with the
///       X axis of 0, 30, 60, or 90 degrees, with the goal of altering the
///       initial orientation as little as possible.
#[pyfunction]
#[pyo3(name = "StraightenDepiction", signature = (mol, confId = -1, minimizeRotation = false))]
#[allow(non_snake_case)]
fn straighten_depiction_py(
    mut mol: PyRefMut<'_, PyMol>,
    confId: i32,
    minimizeRotation: bool,
) -> PyResult<()> {
    straighten_depiction(&mut mol.mol, confId, minimizeRotation).py()
}

/// Normalizes the 2D depiction.
/// If canonicalize is != 0, the depiction is subjected to a canonical
/// transformation such that its main axis is aligned along the X axis
/// (canonicalize >0, the default) or the Y axis (canonicalize <0).
/// If canonicalize is 0, no canonicalization takes place.
/// If scaleFactor is <0.0 (the default) the depiction is scaled such
/// that bond lengths conform to RDKit standards. The applied scaling
/// factor is returned.
///
///   ARGUMENTS:
///
///     - mol: the molecule to be normalized
///     - confId: (optional) the id of the reference conformation to use
///     - canonicalize: (optional) if != 0, a canonical transformation is
///       applied: if >0 (the default), the main molecule axis is aligned
///       to the X axis, if <0 to the Y axis. If 0, no canonical
///       transformation is applied.
///     - scaleFactor: (optional) if >0.0, the scaling factor to apply. The
///       default (-1.0) means that the depiction is automatically scaled
///       such that the average bond length is the standard RDKit bond
///       length.
///
///   RETURNS:
///
///     - the applied scaling factor.
#[pyfunction]
#[pyo3(name = "NormalizeDepiction", signature = (mol, confId = -1, canonicalize = 1, scaleFactor = -1.0))]
#[allow(non_snake_case)]
fn normalize_depiction_py(
    mut mol: PyRefMut<'_, PyMol>,
    confId: i32,
    canonicalize: i32,
    scaleFactor: f64,
) -> PyResult<f64> {
    normalize_depiction(&mut mol.mol, confId, canonicalize, scaleFactor).py()
}

/// Loads the ring system templates from the specified file to be used in 2D
/// coordinate generation. Each template must be a single line in the file
/// represented using CXSMILES, and the structure should be a single ring
/// system. Throws a DepictException if any templates are invalid.
#[pyfunction]
#[pyo3(name = "SetRingSystemTemplates", signature = (templatePath))]
#[allow(non_snake_case)]
fn set_ring_system_templates_py(templatePath: &str) -> PyResult<()> {
    set_ring_system_templates(templatePath).py()
}

/// Adds the ring system templates from the specified file to be used in 2D
/// coordinate generation. If there are duplicates, the most recently added
/// template will be used. Each template must be a single line in the file
/// represented using CXSMILES, and the structure should be a single ring
/// system. Throws a DepictException if any templates are invalid.
#[pyfunction]
#[pyo3(name = "AddRingSystemTemplates", signature = (templatePath))]
#[allow(non_snake_case)]
fn add_ring_system_templates_py(templatePath: &str) -> PyResult<()> {
    add_ring_system_templates(templatePath).py()
}

/// Loads the default ring system templates and removes existing ones, if present.
#[pyfunction]
#[pyo3(name = "LoadDefaultRingSystemTemplates")]
fn load_default_ring_system_templates_py() {
    load_default_ring_system_templates();
}

/// Returns whether RDKit was built with CoordGen support
/// (always False: CoordGen is not part of rdoxide)
#[pyfunction]
#[pyo3(name = "IsCoordGenSupportAvailable")]
fn is_coord_gen_support_available() -> bool {
    false
}

/// Set whether or not the CoordGen library should be preferred to the RDKit
/// depiction library (which has no effect without CoordGen support).
#[pyfunction]
#[pyo3(name = "SetPreferCoordGen", signature = (val))]
fn set_prefer_coord_gen(val: bool) {
    let _ = val;
}

/// Always returns False (CoordGen support not enabled)
#[pyfunction]
#[pyo3(name = "GetPreferCoordGen")]
fn get_prefer_coord_gen() -> bool {
    false
}

/// Context manager to temporarily set CoordGen library preference in RDKit
/// depiction.
#[pyclass(name = "UsingCoordGen", module = "rdoxide.Chem.rdDepictor")]
#[derive(Debug)]
struct PyUsingCoordGen;

#[pymethods]
impl PyUsingCoordGen {
    #[new]
    #[pyo3(signature = (temp_state))]
    fn new(temp_state: bool) -> Self {
        let _ = temp_state;
        PyUsingCoordGen
    }

    fn __enter__(&self) {}

    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        &self,
        exc_type: &Bound<'_, PyAny>,
        exc_val: &Bound<'_, PyAny>,
        traceback: &Bound<'_, PyAny>,
    ) {
        let _ = (exc_type, exc_val, traceback);
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing the functionality to compute 2D coordinates for a molecule",
    )?;
    m.add_function(wrap_pyfunction!(compute_2d_coords_py, m)?)?;
    m.add_function(wrap_pyfunction!(compute_2d_coords_mimic_distmat_py, m)?)?;
    m.add_function(wrap_pyfunction!(straighten_depiction_py, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_depiction_py, m)?)?;
    m.add_function(wrap_pyfunction!(set_ring_system_templates_py, m)?)?;
    m.add_function(wrap_pyfunction!(add_ring_system_templates_py, m)?)?;
    m.add_function(wrap_pyfunction!(load_default_ring_system_templates_py, m)?)?;
    m.add_function(wrap_pyfunction!(is_coord_gen_support_available, m)?)?;
    m.add_function(wrap_pyfunction!(set_prefer_coord_gen, m)?)?;
    m.add_function(wrap_pyfunction!(get_prefer_coord_gen, m)?)?;
    m.add_class::<PyUsingCoordGen>()?;
    Ok(())
}
