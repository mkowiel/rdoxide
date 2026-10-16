//! `rdoxide.Chem.rdDepictor` (upstream `Depictor/Wrap/rdDepictor.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rdoxide::graphmol::depictor::depict_utils::{IntPoint2DMap, bond_len, set_bond_len};
use rdoxide::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};

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
#[pyo3(name = "Compute2DCoords", signature = (mol, canonOrient = true, clearConfs = true, coordMap = None, nFlipsPerSample = 0, nSample = 0, sampleSeed = 0, permuteDeg4Nodes = false, bondLength = -1.0, forceRDKit = false, useRingTemplates = false, usePathAngleExpansion = true))]
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
    permuteDeg4Nodes: bool,
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
        permute_deg4_nodes: permuteDeg4Nodes,
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

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing the functionality to compute 2D coordinates for a molecule",
    )?;
    m.add_function(wrap_pyfunction!(compute_2d_coords_py, m)?)?;
    Ok(())
}
