//! Molecule standardization (upstream `GraphMol/MolStandardize`, the RDKit
//! port of MolVS): cleanup, normalization, reionization, fragment and charge
//! parents, tautomers and validation.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/MolStandardize.{h,cpp}`.

pub mod catalogs;
pub mod charge;
pub mod fragment;
pub mod metal;
pub mod normalize;
pub mod pipeline;
pub mod tautomer;
pub mod validate;

use std::sync::{Mutex, OnceLock};

use self::charge::Uncharger;
use self::charge::reionizer_from_params;
use self::fragment::LargestFragmentChooser;
use self::fragment::fragment_remover_from_params;
use self::metal::{MetalDisconnector, MetalDisconnectorOptions};
use self::normalize::normalizer_from_params;
use crate::error::{RDKitError, Result};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::mol_ops::structure::remove_stereochemistry;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::rd_log_warning;
use crate::rdgeneral::ptree::PTree;
use crate::rdgeneral::threads::{num_threads_to_use, par_map};

/// `MolStandardize::CleanupParameters`
#[derive(Debug, Clone, PartialEq)]
pub struct CleanupParameters {
    /// `$RDBASE` when the parameters were created (unused by the port).
    pub rdbase: String,
    /// file with the normalization transforms (`normalizations`)
    pub normalizations: String,
    pub acidbase_file: String,
    pub fragment_file: String,
    /// file with the tautomer transforms (`tautomerTransforms`)
    pub tautomer_transforms: String,
    /// maximum number of times the normalizations are restarted
    pub max_restarts: i32,
    /// prefer organic fragments when choosing the fragment parent
    pub prefer_organic: bool,
    /// apply atom-order dependent normalizations in a canonical order
    pub do_canonical: bool,
    pub max_tautomers: i32,
    pub max_transforms: i32,
    pub tautomer_remove_sp3_stereo: bool,
    pub tautomer_remove_bond_stereo: bool,
    pub tautomer_remove_isotopic_hs: bool,
    pub tautomer_reassign_stereo: bool,
    pub largest_fragment_chooser_use_atom_count: bool,
    pub largest_fragment_chooser_count_heavy_atoms_only: bool,
    /// (name, reaction SMARTS) normalizations used instead of a file
    pub normalization_data: Vec<(String, String)>,
    /// (name, SMARTS) fragments used instead of a file
    pub fragment_data: Vec<(String, String)>,
    /// (name, acid SMARTS, base SMARTS) pairs used instead of a file
    pub acidbase_data: Vec<(String, String, String)>,
    /// (name, SMARTS, bonds, charges) transforms used instead of a file
    pub tautomer_transform_data: Vec<(String, String, String, String)>,
}

impl Default for CleanupParameters {
    fn default() -> Self {
        CleanupParameters {
            rdbase: std::env::var("RDBASE").unwrap_or_default(),
            normalizations: String::new(),
            acidbase_file: String::new(),
            fragment_file: String::new(),
            tautomer_transforms: String::new(),
            max_restarts: 200,
            prefer_organic: false,
            do_canonical: true,
            max_tautomers: 1000,
            max_transforms: 1000,
            tautomer_remove_sp3_stereo: true,
            tautomer_remove_bond_stereo: true,
            tautomer_remove_isotopic_hs: true,
            tautomer_reassign_stereo: true,
            largest_fragment_chooser_use_atom_count: true,
            largest_fragment_chooser_count_heavy_atoms_only: false,
            normalization_data: Vec::new(),
            fragment_data: Vec::new(),
            acidbase_data: Vec::new(),
            tautomer_transform_data: Vec::new(),
        }
    }
}

/// `updateCleanupParamsFromJSON`: the keys present in `json` override the
/// corresponding fields; the `*Data` arrays are appended to.
pub fn update_cleanup_params_from_json(params: &mut CleanupParameters, json: &str) -> Result<()> {
    if json.is_empty() {
        return Ok(());
    }
    let pt = PTree::from_json(json)?;
    let get_str = |key: &str, val: &mut String| {
        if let Some(v) = pt.get_str(key) {
            *val = v.to_string();
        }
    };
    get_str("rdbase", &mut params.rdbase);
    get_str("normalizations", &mut params.normalizations);
    get_str("acidbaseFile", &mut params.acidbase_file);
    get_str("fragmentFile", &mut params.fragment_file);
    get_str("tautomerTransforms", &mut params.tautomer_transforms);
    params.max_restarts = pt.get_i32("maxRestarts", params.max_restarts);
    params.prefer_organic = pt.get_bool("preferOrganic", params.prefer_organic);
    params.do_canonical = pt.get_bool("doCanonical", params.do_canonical);
    params.max_tautomers = pt.get_i32("maxTautomers", params.max_tautomers);
    params.max_transforms = pt.get_i32("maxTransforms", params.max_transforms);
    params.tautomer_remove_sp3_stereo =
        pt.get_bool("tautomerRemoveSp3Stereo", params.tautomer_remove_sp3_stereo);
    params.tautomer_remove_bond_stereo = pt.get_bool(
        "tautomerRemoveBondStereo",
        params.tautomer_remove_bond_stereo,
    );
    params.tautomer_remove_isotopic_hs = pt.get_bool(
        "tautomerRemoveIsotopicHs",
        params.tautomer_remove_isotopic_hs,
    );
    params.tautomer_reassign_stereo =
        pt.get_bool("tautomerReassignStereo", params.tautomer_reassign_stereo);

    let field = |entry: &PTree, key: &str| entry.get_str(key).unwrap_or_default().to_string();
    if let Some(entries) = pt.child("normalizationData") {
        for entry in entries.children() {
            let (nm, smarts) = (field(entry, "name"), field(entry, "smarts"));
            if nm.is_empty() || smarts.is_empty() {
                rd_log_warning!(" empty transformation name or SMARTS");
                continue;
            }
            params.normalization_data.push((nm, smarts));
        }
    }
    if let Some(entries) = pt.child("fragmentData") {
        for entry in entries.children() {
            let (nm, smarts) = (field(entry, "name"), field(entry, "smarts"));
            if nm.is_empty() || smarts.is_empty() {
                rd_log_warning!(" empty transformation name or SMARTS");
                continue;
            }
            params.fragment_data.push((nm, smarts));
        }
    }
    if let Some(entries) = pt.child("acidbaseData") {
        for entry in entries.children() {
            let (nm, acid, base) = (
                field(entry, "name"),
                field(entry, "acid"),
                field(entry, "base"),
            );
            if nm.is_empty() || acid.is_empty() || base.is_empty() {
                rd_log_warning!(" empty component in acidbaseData");
                continue;
            }
            params.acidbase_data.push((nm, acid, base));
        }
    }
    if let Some(entries) = pt.child("tautomerTransformData") {
        for entry in entries.children() {
            let (nm, smarts) = (field(entry, "name"), field(entry, "smarts"));
            let (bonds, charges) = (field(entry, "bonds"), field(entry, "charges"));
            if nm.is_empty() || smarts.is_empty() {
                rd_log_warning!(" empty component in tautomerTransformData");
                continue;
            }
            params
                .tautomer_transform_data
                .push((nm, smarts, bonds, charges));
        }
    }
    Ok(())
}

/// `removeFragments`
pub fn remove_fragments(mol: &RWMol, params: &CleanupParameters) -> Result<RWMol> {
    let mut res = mol.clone();
    remove_fragments_in_place(&mut res, params)?;
    Ok(res)
}

/// `removeFragmentsInPlace`
pub fn remove_fragments_in_place(mol: &mut RWMol, params: &CleanupParameters) -> Result<()> {
    fragment_remover_from_params(params, true, false)?.remove_in_place(mol)
}

/// `defaultCleanupParameters`
pub fn default_cleanup_parameters() -> &'static CleanupParameters {
    static DEFAULT: OnceLock<CleanupParameters> = OnceLock::new();
    DEFAULT.get_or_init(CleanupParameters::default)
}

/// `standardizeMultipleMolsInPlace`: applies `func` to every molecule on up
/// to `num_threads` threads; the first error is returned.
pub fn standardize_multiple_mols_in_place<F>(
    mols: &mut [&mut RWMol],
    num_threads: i32,
    func: F,
) -> Result<()>
where
    F: Fn(&mut RWMol) -> Result<()> + Sync,
{
    let threads = num_threads_to_use(num_threads).min(mols.len());
    if threads <= 1 {
        return mols.iter_mut().try_for_each(|m| func(m));
    }
    let slots: Vec<Mutex<&mut &mut RWMol>> = mols.iter_mut().map(Mutex::new).collect();
    par_map(&slots, threads as i32, |slot| {
        let mut guard = slot.lock().unwrap_or_else(|e| e.into_inner());
        func(&mut guard)
    })
    .into_iter()
    .collect()
}

/// `normalize`: the fragment-wise normalization (`Normalizer::normalize`).
pub fn normalize(mol: &RWMol, params: &CleanupParameters) -> Result<RWMol> {
    normalizer_from_params(params)?.normalize(mol)
}

/// `normalizeInPlace`
pub fn normalize_in_place(mol: &mut RWMol, params: &CleanupParameters) -> Result<()> {
    normalizer_from_params(params)?.normalize_in_place(mol)
}

/// `reionize`
pub fn reionize(mol: &RWMol, params: &CleanupParameters) -> Result<RWMol> {
    reionizer_from_params(params)?.reionize(mol)
}

/// `reionizeInPlace`
pub fn reionize_in_place(mol: &mut RWMol, params: &CleanupParameters) -> Result<()> {
    reionizer_from_params(params)?.reionize_in_place(mol)
}

/// `disconnectOrganometallics(RWMol &, mdo)`; `mdo` defaults to
/// [`MetalDisconnectorOptions::organometallics`].
pub fn disconnect_organometallics_in_place(
    mol: &mut RWMol,
    mdo: Option<MetalDisconnectorOptions>,
) -> Result<()> {
    let mdo = mdo.unwrap_or_else(MetalDisconnectorOptions::organometallics);
    MetalDisconnector::new(mdo).disconnect_in_place(mol)
}

/// `disconnectOrganometallics(const ROMol &, mdo)`
pub fn disconnect_organometallics(
    mol: &ROMol,
    mdo: Option<MetalDisconnectorOptions>,
) -> Result<ROMol> {
    let mdo = mdo.unwrap_or_else(MetalDisconnectorOptions::organometallics);
    MetalDisconnector::new(mdo).disconnect(mol)
}

/// `canonicalTautomer`
pub fn canonical_tautomer(mol: &RWMol, params: &CleanupParameters) -> Result<RWMol> {
    tautomer::tautomer_enumerator_from_params(params)?.canonicalize(mol, None)
}

/// `canonicalTautomerInPlace`
pub fn canonical_tautomer_in_place(mol: &mut RWMol, params: &CleanupParameters) -> Result<()> {
    tautomer::tautomer_enumerator_from_params(params)?.canonicalize_in_place(mol, None)
}

// ------------------------------------------------------------- top level

/// `cleanupInPlace`: removeHs, disconnect metals, normalize, reionize,
/// assign stereochemistry.
pub fn cleanup_in_place(mol: &mut RWMol, params: &CleanupParameters) -> Result<()> {
    remove_hs(mol, &RemoveHsParameters::default(), true)?;
    MetalDisconnector::default().disconnect_in_place(mol)?;
    normalize_in_place(mol, params)?;
    reionize_in_place(mol, params)?;
    assign_stereochemistry(mol, true, true, false)
}

/// `cleanup`
pub fn cleanup(mol: &RWMol, params: &CleanupParameters) -> Result<RWMol> {
    let mut res = mol.clone();
    cleanup_in_place(&mut res, params)?;
    Ok(res)
}

/// `tautomerParentInPlace`: the cleaned-up canonical tautomer.
pub fn tautomer_parent_in_place(
    mol: &mut RWMol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    if !skip_standardize {
        cleanup_in_place(mol, params)?;
    }
    canonical_tautomer_in_place(mol, params)?;
    cleanup_in_place(mol, params)
}

/// `fragmentParentInPlace`: the largest fragment.
pub fn fragment_parent_in_place(
    mol: &mut RWMol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    if !skip_standardize {
        cleanup_in_place(mol, params)?;
    }
    LargestFragmentChooser::new(params.prefer_organic).choose_in_place(mol)
}

/// `stereoParentInPlace`: without stereochemistry.
pub fn stereo_parent_in_place(
    mol: &mut RWMol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    if !skip_standardize {
        cleanup_in_place(mol, params)?;
    }
    remove_stereochemistry(mol)
}

/// `isotopeParentInPlace`: without isotopes.
pub fn isotope_parent_in_place(
    mol: &mut RWMol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    if !skip_standardize {
        cleanup_in_place(mol, params)?;
    }
    for i in 0..mol.num_atoms() {
        mol.atom_mut(i).set_isotope(0);
    }
    Ok(())
}

/// `chargeParentInPlace`: the uncharged fragment parent.
pub fn charge_parent_in_place(
    mol: &mut RWMol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    fragment_parent_in_place(mol, params, skip_standardize)?;
    Uncharger::new(params.do_canonical, false, false).uncharge_in_place(mol)?;
    cleanup_in_place(mol, params)
}

/// `superParentInPlace`: the fragment, charge, isotope, stereo and
/// tautomer parent.
pub fn super_parent_in_place(
    mol: &mut RWMol,
    params: &CleanupParameters,
    skip_standardize: bool,
) -> Result<()> {
    if !skip_standardize {
        cleanup_in_place(mol, params)?;
    }
    // chargeParent includes fragmentParent
    charge_parent_in_place(mol, params, true)?;
    isotope_parent_in_place(mol, params, true)?;
    stereo_parent_in_place(mol, params, true)?;
    tautomer_parent_in_place(mol, params, true)?;
    cleanup_in_place(mol, params)
}

macro_rules! copying_parent {
    ($($(#[$doc:meta])* $name:ident => $in_place:ident;)*) => {
        $(
            $(#[$doc])*
            pub fn $name(mol: &RWMol, params: &CleanupParameters, skip_standardize: bool) -> Result<RWMol> {
                let mut res = mol.clone();
                $in_place(&mut res, params, skip_standardize)?;
                Ok(res)
            }
        )*
    };
}

copying_parent! {
    /// `tautomerParent`
    tautomer_parent => tautomer_parent_in_place;
    /// `fragmentParent`
    fragment_parent => fragment_parent_in_place;
    /// `stereoParent`
    stereo_parent => stereo_parent_in_place;
    /// `isotopeParent`
    isotope_parent => isotope_parent_in_place;
    /// `chargeParent`
    charge_parent => charge_parent_in_place;
    /// `superParent`
    super_parent => super_parent_in_place;
}

/// `standardizeSmiles`: the canonical SMILES of the cleaned-up molecule.
pub fn standardize_smiles(smiles: &str) -> Result<String> {
    let ps = SmilesParserParams {
        sanitize: false,
        remove_hs: false,
        ..Default::default()
    };
    let mut mol = mol_from_smiles(smiles, &ps)?.ok_or_else(|| {
        RDKitError::Value(format!(
            "SMILES Parse Error: syntax error for input: {smiles}"
        ))
    })?;
    cleanup_in_place(&mut mol, default_cleanup_parameters())?;
    mol_to_smiles(&mut mol, &SmilesWriteParams::default())
}

/// `enumerateTautomerSmiles`
pub fn enumerate_tautomer_smiles(smiles: &str, params: &CleanupParameters) -> Result<Vec<String>> {
    let ps = SmilesParserParams {
        sanitize: false,
        remove_hs: false,
        ..Default::default()
    };
    let mut mol = mol_from_smiles(smiles, &ps)?.ok_or_else(|| {
        RDKitError::Value(format!(
            "SMILES Parse Error: syntax error for input: {smiles}"
        ))
    })?;
    cleanup_in_place(&mut mol, params)?;
    sanitize_mol(&mut mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
    let te = tautomer::TautomerEnumerator::from_params(params)?;
    Ok(te.enumerate(&mol)?.smiles())
}
