//! Enhanced-stereo helpers of `rdkit.Chem.rdmolops` (upstream
//! `GraphMol/Wrap/MolOps.cpp` and `ChiralityOps.cpp`): `AddStereoAnnotations`,
//! `SimplifyEnhancedStereo`, `FindMesoCenters`, `CleanupStereoGroups` and
//! `TranslateChiralFlagToStereoGroups`.

use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::chirality::stereo_groups::{
    StereoAnnotationLabels, add_stereo_annotations, find_meso_centers, simplify_enhanced_stereo,
    translate_chiral_flag_to_stereo_groups,
};
use rdoxide::graphmol::mol_ops::cleanup_stereo_groups;
use rdoxide::graphmol::stereo_group::StereoGroupType;

use crate::enums::to_py;
use crate::enums_impl::PyStereoGroupType;
use crate::errors::{IntoPyResult, argument_error, boost_type_name};
use crate::mol::PyMol;

/// add R/S, relative stereo, and E/Z annotations to atoms and bonds
///
///   Arguments:
///    - mol: molecule to modify
///    - absLabel: label for atoms in an ABS stereo group
///    - orLabel: label for atoms in an OR stereo group
///    - andLabel: label for atoms in an AND stereo group
///    - cipLabel: label for chiral atoms that aren't in a stereo group.
///    - bondLabel: label for CIP stereochemistry on bonds
///
///  If any label is empty, the corresponding annotations will not be added.
///
///  The labels can contain the following placeholders:
///    - {id} - the stereo group's index
///    - {cip} - the atom or bond's CIP stereochemistry
///
///  Note that CIP labels will only be added if CIP stereochemistry has been
///  assigned to the molecule.
#[pyfunction]
#[pyo3(
    name = "AddStereoAnnotations",
    signature = (mol, absLabel = "abs ({cip})", orLabel = "or{id}", andLabel = "and{id}", cipLabel = "({cip})", bondLabel = "({cip})")
)]
#[allow(non_snake_case)]
fn py_add_stereo_annotations(
    mut mol: PyRefMut<'_, PyMol>,
    absLabel: &str,
    orLabel: &str,
    andLabel: &str,
    cipLabel: &str,
    bondLabel: &str,
) {
    let labels = StereoAnnotationLabels {
        abs: absLabel.into(),
        or: orLabel.into(),
        and: andLabel.into(),
        cip: cipLabel.into(),
        bond: bondLabel.into(),
    };
    add_stereo_annotations(&mut mol.mol, &labels);
}

/// Simplifies the stereochemical representation of a molecule where all
/// specified stereocenters are in the same StereoGroup
///
///   Arguments:
///    - mol: molecule to modify
///    - removeAffectedStereoGroups: if set then the affected StereoGroups will be removed
///
/// If all specified stereocenters are in the same AND or OR stereogroup, a
/// moleculeNote property will be set on the molecule with the value "AND
/// enantiomer" or "OR enantiomer". CIP labels, if present, are removed.
#[pyfunction]
#[pyo3(name = "SimplifyEnhancedStereo", signature = (mol, removeAffectedStereoGroups = true))]
#[allow(non_snake_case)]
fn py_simplify_enhanced_stereo(mut mol: PyRefMut<'_, PyMol>, removeAffectedStereoGroups: bool) {
    simplify_enhanced_stereo(&mut mol.mol, removeAffectedStereoGroups);
}

/// returns the meso centers in a molecule (if any).
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - includeIsotopes: (optional) toggles whether or not isotopes should be included in the
///       calculation.
///     - includeAtomMaps: (optional) toggles whether or not atom maps should be included in the
///       calculation.
#[pyfunction]
#[pyo3(name = "FindMesoCenters", signature = (mol, includeIsotopes = true, includeAtomMaps = false))]
#[allow(non_snake_case)]
fn py_find_meso_centers<'py>(
    py: Python<'py>,
    mut mol: PyRefMut<'py, PyMol>,
    includeIsotopes: bool,
    includeAtomMaps: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    let pairs = find_meso_centers(&mut mol.mol, includeIsotopes, includeAtomMaps).py()?;
    let lists = pairs
        .into_iter()
        .map(|(i, j)| PyTuple::new(py, [i, j]))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, lists)
}

/// removes atoms without specified chirality from stereo groups
#[pyfunction]
#[pyo3(name = "CleanupStereoGroups", signature = (mol))]
fn py_cleanup_stereo_groups(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    cleanup_stereo_groups(&mut mol.mol).py()
}

/// Generate enhanced stereo groups based on the status of the chiral flag property.
///
///   Arguments:
///    - mol: molecule to be modified
///    - zeroFlagGroupType: how to handle non-grouped stereo centers when the
///         chiral flag is set to zero
///
///   If the chiral flag is set to a value of 1 then all specified tetrahedral
///   chiral centers which are not already in StereoGroups will be added to an
///   ABS StereoGroup.
///
///   If the chiral flag is set to a value of 0 then all specified tetrahedral
///   chiral centers will be added to a StereoGroup of the type zeroFlagGroupType
///
///   If there is no chiral flag set (i.e. the property is not present), the
///   molecule will not be modified.
#[pyfunction]
#[pyo3(
    name = "TranslateChiralFlagToStereoGroups",
    signature = (mol, zeroFlagGroupType = None)
)]
#[allow(non_snake_case)]
fn py_translate_chiral_flag_to_stereo_groups(
    mut mol: PyRefMut<'_, PyMol>,
    zeroFlagGroupType: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let group_type = match zeroFlagGroupType {
        None => StereoGroupType::And,
        Some(t) => {
            let py = t.py();
            let class = to_py(py, "StereoGroupType", 0)?.bind(py).get_type();
            if !t.is_instance(&class)? {
                return Err(argument_error(
                    py,
                    "TranslateChiralFlagToStereoGroups",
                    &["Mol".to_string(), boost_type_name(t)],
                    "TranslateChiralFlagToStereoGroups(RDKit::ROMol {lvalue} mol, RDKit::StereoGroupType zeroFlagGroupType=rdkit.Chem.rdchem.StereoGroupType.STEREO_AND)",
                ));
            }
            PyStereoGroupType::extract(t)?
        }
    };
    translate_chiral_flag_to_stereo_groups(&mut mol.mol, group_type).py()
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_add_stereo_annotations, m)?)?;
    m.add_function(wrap_pyfunction!(py_simplify_enhanced_stereo, m)?)?;
    m.add_function(wrap_pyfunction!(py_find_meso_centers, m)?)?;
    m.add_function(wrap_pyfunction!(py_cleanup_stereo_groups, m)?)?;
    m.add_function(wrap_pyfunction!(
        py_translate_chiral_flag_to_stereo_groups,
        m
    )?)?;
    Ok(())
}
