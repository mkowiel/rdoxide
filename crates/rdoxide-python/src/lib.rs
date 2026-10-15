//! Python bindings for rdoxide, API-compatible with RDKit's Boost.Python
//! wrappers. Built as the extension module `rdkit._rdoxide`; each RDKit
//! extension module (`rdkit.rdBase`, `rdkit.Chem.rdchem`, ...) is a submodule
//! here whose `__name__` is the upstream module name. The pure-Python shims in
//! `python/rdkit` replace themselves in `sys.modules` with these submodules.

use pyo3::prelude::*;
use pyo3::types::PyDict;

mod atom;
mod bond;
mod conformer;
mod datastructs;
mod datastructs_fpb;
mod datastructs_siv;
mod datastructs_vv;
mod debug;
mod enums;
mod enums_impl;
mod errors;
mod geometry;
mod geometry_grids;
mod mol;
mod monomer_info;
mod perception;
mod props;
mod pyconv;
mod rdabbreviations;
mod rdbase;
mod rdchem;
mod rdchem_atom_funcs;
mod rdchem_stereo_info;
mod rdchemicalfeatures;
mod rdchemreactions;
mod rdclustering;
mod rddepictor;
mod rddeterminebonds;
mod rddistgeom;
mod rddistgeom_dg;
mod rdenumeratestereoisomers;
mod rdfiltercatalog;
mod rdfingerprintgenerator;
mod rdfmcs;
mod rdforcefield;
mod rdinfotheory;
mod rdmetricmatrixcalc;
mod rdmhfpfingerprint;
mod rdmmpa;
mod rdmolalign;
mod rdmolchemicalfeatures;
mod rdmoldescriptors;
mod rdmoldescriptors_3d;
mod rdmoldescriptors_autocorr;
mod rdmoldescriptors_bcut;
mod rdmoldescriptors_connectivity;
mod rdmoldescriptors_molsurf;
mod rdmoldescriptors_profiles;
mod rdmoldescriptors_property;
mod rdmoldescriptors_props;
mod rdmoldescriptors_usr;
mod rdmolfiles;
mod rdmolfiles_fwdsdsupplier;
mod rdmolfiles_molblock;
mod rdmolfiles_molwrite;
mod rdmolfiles_mtsupplier;
mod rdmolfiles_pdb;
mod rdmolfiles_sdsupplier;
mod rdmolfiles_sdwriter;
mod rdmolfiles_sequence;
mod rdmolfiles_smarts;
mod rdmolfiles_smisupplier;
mod rdmolfiles_smiwriter;
mod rdmolfiles_write;
mod rdmolfiles_xyz;
mod rdmolhash;
mod rdmolinterchange;
mod rdmolops;
mod rdmolops_adjust_query;
mod rdmolops_fingerprints;
mod rdmolops_fragments;
mod rdmolops_hs;
mod rdmolops_matrices;
mod rdmolops_stereo_groups;
mod rdmolops_stereo_misc;
mod rdmolops_structure;
mod rdmolops_subgraphs;
mod rdmolops_subset;
mod rdmolops_transforms;
mod rdmolops_wedging;
mod rdmolops_zip;
mod rdmolstandardize;
mod rdmolstandardize_pipeline;
mod rdmolstandardize_tautomer;
mod rdmoltransforms;
mod rdpartialcharges;
mod rdqueries;
mod rdreducedgraphs;
mod rdresonance;
mod rdsimdivpickers;
mod rdsubstructlibrary;
mod ringinfo;
mod stereo_group;
mod substance_group;
mod substruct;
mod vector_suite;
mod vects;
mod writer_sink;

/// Create a submodule named `full_name` (e.g. `rdkit.Chem.rdchem`) and attach
/// it to `parent` under `attr`.
fn add_submodule<'py>(
    parent: &Bound<'py, PyModule>,
    attr: &str,
    full_name: &str,
    init: impl FnOnce(&Bound<'py, PyModule>) -> PyResult<()>,
) -> PyResult<()> {
    let py = parent.py();
    let m = PyModule::new(py, attr)?;
    m.setattr("__name__", full_name)?;
    init(&m)?;
    parent.add(attr, &m)?;
    // make `import rdkit.x.y` find it even before the shim runs
    let modules = py
        .import("sys")?
        .getattr("modules")?
        .cast_into::<PyDict>()?;
    modules.set_item(full_name, &m)?;
    Ok(())
}

#[pymodule]
fn _rdoxide(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", rdoxide::rdgeneral::RDOXIDE_VERSION)?;
    add_submodule(m, "rdBase", "rdkit.rdBase", rdbase::init)?;
    add_submodule(m, "rdGeometry", "rdkit.Geometry.rdGeometry", geometry::init)?;
    add_submodule(m, "rdchem", "rdkit.Chem.rdchem", rdchem::init)?;
    add_submodule(m, "rdmolfiles", "rdkit.Chem.rdmolfiles", rdmolfiles::init)?;
    add_submodule(m, "rdmolops", "rdkit.Chem.rdmolops", rdmolops::init)?;
    add_submodule(m, "rdDepictor", "rdkit.Chem.rdDepictor", rddepictor::init)?;
    add_submodule(m, "rdqueries", "rdkit.Chem.rdqueries", rdqueries::init)?;
    add_submodule(m, "rdMolAlign", "rdkit.Chem.rdMolAlign", rdmolalign::init)?;
    add_submodule(m, "rdDistGeom", "rdkit.Chem.rdDistGeom", rddistgeom::init)?;
    add_submodule(m, "rdFMCS", "rdkit.Chem.rdFMCS", rdfmcs::init)?;
    add_submodule(
        m,
        "rdfiltercatalog",
        "rdkit.Chem.rdfiltercatalog",
        rdfiltercatalog::init,
    )?;
    add_submodule(
        m,
        "DistGeom",
        "rdkit.DistanceGeometry.DistGeom",
        rddistgeom_dg::init,
    )?;
    add_submodule(
        m,
        "rdForceField",
        "rdkit.ForceField.rdForceField",
        rdforcefield::init_force_field,
    )?;
    add_submodule(
        m,
        "rdForceFieldHelpers",
        "rdkit.Chem.rdForceFieldHelpers",
        rdforcefield::init_helpers,
    )?;
    add_submodule(
        m,
        "rdChemicalFeatures",
        "rdkit.Chem.rdChemicalFeatures",
        rdchemicalfeatures::init,
    )?;
    add_submodule(
        m,
        "rdChemReactions",
        "rdkit.Chem.rdChemReactions",
        rdchemreactions::init,
    )?;
    add_submodule(
        m,
        "rdMolChemicalFeatures",
        "rdkit.Chem.rdMolChemicalFeatures",
        rdmolchemicalfeatures::init,
    )?;
    add_submodule(
        m,
        "rdMolTransforms",
        "rdkit.Chem.rdMolTransforms",
        rdmoltransforms::init,
    )?;
    add_submodule(
        m,
        "rdMolDescriptors",
        "rdkit.Chem.rdMolDescriptors",
        rdmoldescriptors::init,
    )?;
    add_submodule(
        m,
        "rdFingerprintGenerator",
        "rdkit.Chem.rdFingerprintGenerator",
        rdfingerprintgenerator::init,
    )?;
    add_submodule(
        m,
        "rdInfoTheory",
        "rdkit.ML.InfoTheory.rdInfoTheory",
        rdinfotheory::init,
    )?;
    add_submodule(
        m,
        "rdReducedGraphs",
        "rdkit.Chem.rdReducedGraphs",
        rdreducedgraphs::init,
    )?;
    add_submodule(
        m,
        "rdSubstructLibrary",
        "rdkit.Chem.rdSubstructLibrary",
        rdsubstructlibrary::init,
    )?;
    add_submodule(m, "rdMolHash", "rdkit.Chem.rdMolHash", rdmolhash::init)?;
    add_submodule(
        m,
        "rdMolStandardize",
        "rdkit.Chem.MolStandardize.rdMolStandardize",
        rdmolstandardize::init,
    )?;
    add_submodule(m, "rdMMPA", "rdkit.Chem.rdMMPA", rdmmpa::init)?;
    add_submodule(
        m,
        "rdMolInterchange",
        "rdkit.Chem.rdMolInterchange",
        rdmolinterchange::init,
    )?;
    add_submodule(
        m,
        "rdDetermineBonds",
        "rdkit.Chem.rdDetermineBonds",
        rddeterminebonds::init,
    )?;
    add_submodule(
        m,
        "rdAbbreviations",
        "rdkit.Chem.rdAbbreviations",
        rdabbreviations::init,
    )?;
    add_submodule(
        m,
        "rdMHFPFingerprint",
        "rdkit.Chem.rdMHFPFingerprint",
        rdmhfpfingerprint::init,
    )?;
    add_submodule(
        m,
        "rdEnumerateStereoisomers",
        "rdkit.Chem.rdEnumerateStereoisomers",
        rdenumeratestereoisomers::init,
    )?;
    add_submodule(
        m,
        "rdPartialCharges",
        "rdkit.Chem.rdPartialCharges",
        rdpartialcharges::init,
    )?;
    add_submodule(
        m,
        "rdMetricMatrixCalc",
        "rdkit.DataManip.Metric.rdMetricMatrixCalc",
        rdmetricmatrixcalc::init,
    )?;
    add_submodule(
        m,
        "Clustering",
        "rdkit.ML.Cluster.Clustering",
        rdclustering::init,
    )?;
    add_submodule(
        m,
        "rdSimDivPickers",
        "rdkit.SimDivFilters.rdSimDivPickers",
        rdsimdivpickers::init,
    )?;
    add_submodule(
        m,
        "cDataStructs",
        "rdkit.DataStructs.cDataStructs",
        datastructs::init,
    )?;
    Ok(())
}
