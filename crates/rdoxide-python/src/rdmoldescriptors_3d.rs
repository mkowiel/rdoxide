//! The 3D shape descriptors of `rdkit.Chem.rdMolDescriptors` (upstream
//! `Descriptors/PMI.cpp`, `PBF.cpp` and their wrappers in
//! `Wrap/rdMolDescriptors.cpp`): CalcPBF, CalcNPR1/2, CalcPMI1/2/3,
//! CalcRadiusOfGyration, CalcInertialShapeFactor, CalcEccentricity,
//! CalcAsphericity and CalcSpherocityIndex.

use pyo3::prelude::*;
use rdoxide::graphmol::descriptors::pbf::{PBF_VERSION, pbf};
use rdoxide::graphmol::descriptors::pmi as pmi_mod;
use rdoxide::graphmol::romol::ROMol;

use crate::errors::{Check, IntoPyResult, Site, check_failed};
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;

/// The failed `PRECONDITION(mol.getNumConformers() >= 1, ...)` of the entry
/// point at `file:line`; its text names the C++ source line.
pub(crate) fn require_conformer(mol: &ROMol, file: &'static str, line: u32) -> PyResult<()> {
    if mol.num_conformers() == 0 {
        return Err(check_failed(
            Site {
                check: Check::Precondition,
                file,
                line,
                expr: "mol.getNumConformers() >= 1",
            },
            "molecule has no conformers",
        ));
    }
    Ok(())
}

const PMI_FILE: &str = "Code/GraphMol/Descriptors/PMI.cpp";

/// Defines `CalcXxx(mol, confId=-1, useAtomicMasses=True, force=True)` for a
/// moment-based descriptor; `$line` is the line of its precondition upstream.
macro_rules! moment_descriptor {
    ($fn_name:ident, $py_name:literal, $core:path, $line:literal) => {
        #[pyfunction]
        #[pyo3(name = $py_name, signature = (mol, confId = -1, useAtomicMasses = None, force = None))]
        #[allow(non_snake_case)]
        fn $fn_name(
            py: Python<'_>,
            mol: &Bound<'_, PyMol>,
            confId: i32,
            useAtomicMasses: Option<&Bound<'_, PyAny>>,
            force: Option<&Bound<'_, PyAny>>,
        ) -> PyResult<f64> {
            let use_masses = bool_arg(useAtomicMasses, true)?;
            let force = bool_arg(force, true)?;
            let mut m = mol.borrow_mut();
            let core = &mut m.mol;
            require_conformer(core, PMI_FILE, $line)?;
            py.detach(|| $core(core, confId, use_masses, force)).py()
        }
    };
}

moment_descriptor!(calc_npr1, "CalcNPR1", pmi_mod::npr1, 110);
moment_descriptor!(calc_npr2, "CalcNPR2", pmi_mod::npr2, 122);
moment_descriptor!(calc_pmi1, "CalcPMI1", pmi_mod::pmi1, 134);
moment_descriptor!(calc_pmi2, "CalcPMI2", pmi_mod::pmi2, 143);
moment_descriptor!(calc_pmi3, "CalcPMI3", pmi_mod::pmi3, 152);
moment_descriptor!(
    calc_radius_of_gyration,
    "CalcRadiusOfGyration",
    pmi_mod::radius_of_gyration,
    163
);
moment_descriptor!(
    calc_inertial_shape_factor,
    "CalcInertialShapeFactor",
    pmi_mod::inertial_shape_factor,
    175
);
moment_descriptor!(
    calc_eccentricity,
    "CalcEccentricity",
    pmi_mod::eccentricity,
    190
);
moment_descriptor!(
    calc_asphericity,
    "CalcAsphericity",
    pmi_mod::asphericity,
    206
);

#[pyfunction]
#[pyo3(name = "CalcSpherocityIndex", signature = (mol, confId = -1, force = None))]
#[allow(non_snake_case)]
fn calc_spherocity_index(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    confId: i32,
    force: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let force = bool_arg(force, true)?;
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    require_conformer(core, PMI_FILE, 224)?;
    py.detach(|| pmi_mod::spherocity_index(core, confId, force))
        .py()
}

/// Returns the PBF (plane of best fit) descriptor (https://doi.org/10.1021/ci300293f)
#[pyfunction]
#[pyo3(name = "CalcPBF", signature = (mol, confId = -1))]
#[allow(non_snake_case)]
fn calc_pbf(py: Python<'_>, mol: &Bound<'_, PyMol>, confId: i32) -> PyResult<f64> {
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    require_conformer(core, "Code/GraphMol/Descriptors/PBF.cpp", 115)?;
    py.detach(|| pbf(core, confId)).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("_CalcPBF_version", PBF_VERSION)?;
    m.add("_CalcNPR1_version", pmi_mod::NPR1_VERSION)?;
    m.add("_CalcNPR2_version", pmi_mod::NPR2_VERSION)?;
    m.add("_CalcPMI1_version", pmi_mod::PMI1_VERSION)?;
    m.add("_CalcPMI2_version", pmi_mod::PMI2_VERSION)?;
    m.add("_CalcPMI3_version", pmi_mod::PMI3_VERSION)?;
    m.add(
        "_CalcRadiusOfGyration_version",
        pmi_mod::RADIUS_OF_GYRATION_VERSION,
    )?;
    m.add(
        "_CalcInertialShapeFactor_version",
        pmi_mod::INERTIAL_SHAPE_FACTOR_VERSION,
    )?;
    m.add("_CalcEccentricity_version", pmi_mod::ECCENTRICITY_VERSION)?;
    m.add("_CalcAsphericity_version", pmi_mod::ASPHERICITY_VERSION)?;
    m.add(
        "_CalcSpherocityIndex_version",
        pmi_mod::SPHEROCITY_INDEX_VERSION,
    )?;
    m.add_function(wrap_pyfunction!(calc_pbf, m)?)?;
    m.add_function(wrap_pyfunction!(calc_npr1, m)?)?;
    m.add_function(wrap_pyfunction!(calc_npr2, m)?)?;
    m.add_function(wrap_pyfunction!(calc_pmi1, m)?)?;
    m.add_function(wrap_pyfunction!(calc_pmi2, m)?)?;
    m.add_function(wrap_pyfunction!(calc_pmi3, m)?)?;
    m.add_function(wrap_pyfunction!(calc_radius_of_gyration, m)?)?;
    m.add_function(wrap_pyfunction!(calc_inertial_shape_factor, m)?)?;
    m.add_function(wrap_pyfunction!(calc_eccentricity, m)?)?;
    m.add_function(wrap_pyfunction!(calc_asphericity, m)?)?;
    m.add_function(wrap_pyfunction!(calc_spherocity_index, m)?)
}
