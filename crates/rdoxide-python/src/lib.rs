//! Python bindings for rdoxide, API-compatible with RDKit's Boost.Python
//! wrappers. Built as the extension module `rdoxide._rdoxide`; each RDKit
//! extension module (`rdoxide.rdBase`, `rdoxide.Chem.rdchem`, ...) is a submodule
//! here whose `__name__` is the upstream module name. The pure-Python shims in
//! `python/rdoxide` replace themselves in `sys.modules` with these submodules.

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
mod molbundle;
mod monomer_info;
mod perception;
mod props;
mod pyconv;
mod rdabbreviations;
mod rdalignment;
mod rdbase;
mod rdchem;
mod rdchem_atom_funcs;
mod rdchem_stereo_info;
mod rdchemicalfeatures;
mod rdchemreactions;
mod rdchemreactions_enumerate;
mod rdcip_labeler;
mod rdclustering;
mod rddepictor;
mod rddepictor_constrained;
mod rddeprotect;
mod rddeterminebonds;
mod rddistgeom;
mod rddistgeom_dg;
mod rdenumeratestereoisomers;
mod rdfiltercatalog;
mod rdfingerprintgenerator;
mod rdfmcs;
mod rdforcefield;
mod rdfragcatalog;
mod rdgaussianshape;
mod rdgeneralizedsubstruct;
mod rdinchi;
mod rdinfotheory;
mod rdmetricmatrixcalc;
mod rdmhfpfingerprint;
mod rdmif;
mod rdmmpa;
mod rdmolalign;
mod rdmolalign_o3a;
mod rdmolcatalog;
mod rdmolchemicalfeatures;
mod rdmoldescriptors;
mod rdmoldescriptors_3d;
mod rdmoldescriptors_autocorr;
mod rdmoldescriptors_bcut;
mod rdmoldescriptors_connectivity;
mod rdmoldescriptors_dclv;
mod rdmoldescriptors_molsurf;
mod rdmoldescriptors_profiles;
mod rdmoldescriptors_property;
mod rdmoldescriptors_props;
mod rdmoldescriptors_usr;
mod rdmoldraw2d;
mod rdmolenumerator;
mod rdmolfiles;
mod rdmolfiles_cdxml;
mod rdmolfiles_cml;
mod rdmolfiles_fwdsdsupplier;
mod rdmolfiles_maesupplier;
mod rdmolfiles_maewriter;
mod rdmolfiles_marvin;
mod rdmolfiles_mol2;
mod rdmolfiles_molblock;
mod rdmolfiles_molwrite;
mod rdmolfiles_mtsupplier;
mod rdmolfiles_pdb;
mod rdmolfiles_png;
mod rdmolfiles_proplists;
mod rdmolfiles_scsr;
mod rdmolfiles_sdsupplier;
mod rdmolfiles_sdwriter;
mod rdmolfiles_sequence;
mod rdmolfiles_smarts;
mod rdmolfiles_smisupplier;
mod rdmolfiles_smiwriter;
mod rdmolfiles_svg;
mod rdmolfiles_tdtsupplier;
mod rdmolfiles_tdtwriter;
mod rdmolfiles_tpl;
mod rdmolfiles_write;
mod rdmolfiles_xyz;
mod rdmolhash;
mod rdmolinterchange;
mod rdmolops;
mod rdmolops_adjust_query;
mod rdmolops_attachment_points;
mod rdmolops_canon_stereo_groups;
mod rdmolops_core_match;
mod rdmolops_fingerprints;
mod rdmolops_fragments;
mod rdmolops_hs;
mod rdmolops_matrices;
mod rdmolops_recursive_queries;
mod rdmolops_stereo_groups;
mod rdmolops_stereo_misc;
mod rdmolops_structure;
mod rdmolops_subgraphs;
mod rdmolops_subset;
mod rdmolops_transforms;
mod rdmolops_wedging;
mod rdmolops_zip;
mod rdmolprocessing;
mod rdmolstandardize;
mod rdmolstandardize_pipeline;
mod rdmolstandardize_tautomer;
mod rdmoltransforms;
mod rdpartialcharges;
mod rdqueries;
mod rdrascalmces;
mod rdreducedgraphs;
mod rdresonance;
mod rdrgroupdecomposition;
mod rdscaffoldnetwork;
mod rdshapehelpers;
mod rdsimdivpickers;
mod rdslnparse;
mod rdstructchecker;
mod rdsubstructlibrary;
mod rdsynthonspacesearch;
mod rdtautomerquery;
mod rdtrajectory;
mod ringinfo;
mod stereo_group;
mod substance_group;
mod substruct;
mod vector_suite;
mod vects;
mod writer_sink;

/// Create a submodule named `full_name` (e.g. `rdoxide.Chem.rdchem`) and attach
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
    // make `import rdoxide.x.y` find it even before the shim runs
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
    add_submodule(m, "rdBase", "rdoxide.rdBase", rdbase::init)?;
    add_submodule(
        m,
        "rdGeometry",
        "rdoxide.Geometry.rdGeometry",
        geometry::init,
    )?;
    add_submodule(m, "rdchem", "rdoxide.Chem.rdchem", rdchem::init)?;
    add_submodule(m, "rdmolfiles", "rdoxide.Chem.rdmolfiles", rdmolfiles::init)?;
    add_submodule(m, "rdmolops", "rdoxide.Chem.rdmolops", rdmolops::init)?;
    add_submodule(m, "rdDepictor", "rdoxide.Chem.rdDepictor", rddepictor::init)?;
    add_submodule(
        m,
        "rdMolDraw2D",
        "rdoxide.Chem.Draw.rdMolDraw2D",
        rdmoldraw2d::init,
    )?;
    add_submodule(m, "rdqueries", "rdoxide.Chem.rdqueries", rdqueries::init)?;
    add_submodule(m, "rdMolAlign", "rdoxide.Chem.rdMolAlign", rdmolalign::init)?;
    add_submodule(m, "rdDistGeom", "rdoxide.Chem.rdDistGeom", rddistgeom::init)?;
    add_submodule(
        m,
        "rdTautomerQuery",
        "rdoxide.Chem.rdTautomerQuery",
        rdtautomerquery::init,
    )?;
    add_submodule(
        m,
        "rdtrajectory",
        "rdoxide.Chem.rdtrajectory",
        rdtrajectory::init,
    )?;
    add_submodule(
        m,
        "rdCIPLabeler",
        "rdoxide.Chem.rdCIPLabeler",
        rdcip_labeler::init,
    )?;
    add_submodule(
        m,
        "rdDeprotect",
        "rdoxide.Chem.rdDeprotect",
        rddeprotect::init,
    )?;
    add_submodule(
        m,
        "rdMolEnumerator",
        "rdoxide.Chem.rdMolEnumerator",
        rdmolenumerator::init,
    )?;
    add_submodule(
        m,
        "rdScaffoldNetwork",
        "rdoxide.Chem.Scaffolds.rdScaffoldNetwork",
        rdscaffoldnetwork::init,
    )?;
    add_submodule(
        m,
        "rdGeneralizedSubstruct",
        "rdoxide.Chem.rdGeneralizedSubstruct",
        rdgeneralizedsubstruct::init,
    )?;
    add_submodule(
        m,
        "rdMolCatalog",
        "rdoxide.Chem.rdMolCatalog",
        rdmolcatalog::init,
    )?;
    add_submodule(
        m,
        "rdfragcatalog",
        "rdoxide.Chem.rdfragcatalog",
        rdfragcatalog::init,
    )?;
    add_submodule(
        m,
        "rdMolProcessing",
        "rdoxide.Chem.rdMolProcessing",
        rdmolprocessing::init,
    )?;
    add_submodule(m, "rdFMCS", "rdoxide.Chem.rdFMCS", rdfmcs::init)?;
    add_submodule(
        m,
        "rdRGroupDecomposition",
        "rdoxide.Chem.rdRGroupDecomposition",
        rdrgroupdecomposition::init,
    )?;
    add_submodule(
        m,
        "rdfiltercatalog",
        "rdoxide.Chem.rdfiltercatalog",
        rdfiltercatalog::init,
    )?;
    add_submodule(
        m,
        "DistGeom",
        "rdoxide.DistanceGeometry.DistGeom",
        rddistgeom_dg::init,
    )?;
    add_submodule(
        m,
        "rdForceField",
        "rdoxide.ForceField.rdForceField",
        rdforcefield::init_force_field,
    )?;
    add_submodule(
        m,
        "rdForceFieldHelpers",
        "rdoxide.Chem.rdForceFieldHelpers",
        rdforcefield::init_helpers,
    )?;
    add_submodule(
        m,
        "rdChemicalFeatures",
        "rdoxide.Chem.rdChemicalFeatures",
        rdchemicalfeatures::init,
    )?;
    add_submodule(
        m,
        "rdChemReactions",
        "rdoxide.Chem.rdChemReactions",
        rdchemreactions::init,
    )?;
    add_submodule(
        m,
        "rdMolChemicalFeatures",
        "rdoxide.Chem.rdMolChemicalFeatures",
        rdmolchemicalfeatures::init,
    )?;
    add_submodule(
        m,
        "rdMolTransforms",
        "rdoxide.Chem.rdMolTransforms",
        rdmoltransforms::init,
    )?;
    add_submodule(
        m,
        "rdMolDescriptors",
        "rdoxide.Chem.rdMolDescriptors",
        rdmoldescriptors::init,
    )?;
    add_submodule(
        m,
        "rdFingerprintGenerator",
        "rdoxide.Chem.rdFingerprintGenerator",
        rdfingerprintgenerator::init,
    )?;
    add_submodule(
        m,
        "rdInfoTheory",
        "rdoxide.ML.InfoTheory.rdInfoTheory",
        rdinfotheory::init,
    )?;
    add_submodule(
        m,
        "rdReducedGraphs",
        "rdoxide.Chem.rdReducedGraphs",
        rdreducedgraphs::init,
    )?;
    add_submodule(
        m,
        "rdSubstructLibrary",
        "rdoxide.Chem.rdSubstructLibrary",
        rdsubstructlibrary::init,
    )?;
    add_submodule(m, "rdSLNParse", "rdoxide.Chem.rdSLNParse", rdslnparse::init)?;
    add_submodule(m, "rdinchi", "rdoxide.Chem.rdinchi", rdinchi::init)?;
    add_submodule(
        m,
        "rdStructChecker",
        "rdoxide.Chem.rdStructChecker",
        rdstructchecker::init,
    )?;
    add_submodule(
        m,
        "rdShapeHelpers",
        "rdoxide.Chem.rdShapeHelpers",
        rdshapehelpers::init,
    )?;
    add_submodule(
        m,
        "rdGaussianShape",
        "rdoxide.Chem.rdGaussianShape",
        rdgaussianshape::init,
    )?;
    add_submodule(
        m,
        "rdRascalMCES",
        "rdoxide.Chem.rdRascalMCES",
        rdrascalmces::init,
    )?;
    add_submodule(
        m,
        "rdSynthonSpaceSearch",
        "rdoxide.Chem.rdSynthonSpaceSearch",
        rdsynthonspacesearch::init,
    )?;
    add_submodule(m, "rdMIF", "rdoxide.Chem.rdMIF", rdmif::init)?;
    add_submodule(m, "rdMolHash", "rdoxide.Chem.rdMolHash", rdmolhash::init)?;
    add_submodule(
        m,
        "rdMolStandardize",
        "rdoxide.Chem.MolStandardize.rdMolStandardize",
        rdmolstandardize::init,
    )?;
    add_submodule(m, "rdMMPA", "rdoxide.Chem.rdMMPA", rdmmpa::init)?;
    add_submodule(
        m,
        "rdMolInterchange",
        "rdoxide.Chem.rdMolInterchange",
        rdmolinterchange::init,
    )?;
    add_submodule(
        m,
        "rdDetermineBonds",
        "rdoxide.Chem.rdDetermineBonds",
        rddeterminebonds::init,
    )?;
    add_submodule(
        m,
        "rdAbbreviations",
        "rdoxide.Chem.rdAbbreviations",
        rdabbreviations::init,
    )?;
    add_submodule(
        m,
        "rdMHFPFingerprint",
        "rdoxide.Chem.rdMHFPFingerprint",
        rdmhfpfingerprint::init,
    )?;
    add_submodule(
        m,
        "rdEnumerateStereoisomers",
        "rdoxide.Chem.rdEnumerateStereoisomers",
        rdenumeratestereoisomers::init,
    )?;
    add_submodule(
        m,
        "rdPartialCharges",
        "rdoxide.Chem.rdPartialCharges",
        rdpartialcharges::init,
    )?;
    add_submodule(
        m,
        "rdAlignment",
        "rdoxide.Numerics.rdAlignment",
        rdalignment::init,
    )?;
    add_submodule(
        m,
        "rdMetricMatrixCalc",
        "rdoxide.DataManip.Metric.rdMetricMatrixCalc",
        rdmetricmatrixcalc::init,
    )?;
    add_submodule(
        m,
        "Clustering",
        "rdoxide.ML.Cluster.Clustering",
        rdclustering::init,
    )?;
    add_submodule(
        m,
        "rdSimDivPickers",
        "rdoxide.SimDivFilters.rdSimDivPickers",
        rdsimdivpickers::init,
    )?;
    add_submodule(
        m,
        "cDataStructs",
        "rdoxide.DataStructs.cDataStructs",
        datastructs::init,
    )?;
    Ok(())
}
