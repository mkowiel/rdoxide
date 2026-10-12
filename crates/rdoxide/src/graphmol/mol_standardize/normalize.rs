//! The `Normalizer`: applies the normalization transforms (reaction SMARTS)
//! that correct functional groups and recombine charges.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/Normalize.{h,cpp}`.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use super::CleanupParameters;
use super::catalogs::TransformCatalogParams;
use crate::error::{RDKitError, Result};
use crate::geometry::Point3D;
use crate::graphmol::chem_reactions::ChemicalReaction;
use crate::graphmol::chem_reactions::runner::{run_reactant_in_place, run_reactants};
use crate::graphmol::chem_transforms::combine_mols;
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, fast_find_rings, symmetrize_sssr};
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::mol_ops::sanitize_flags::{
    SANITIZE_ALL, SANITIZE_CLEANUP, SANITIZE_PROPERTIES,
};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::rd_log_info;
use crate::rdgeneral::common_properties;
use crate::rdgeneral::dict::RDProps;

/// Sanitization of transformed molecules: atoms with too high a valence are
/// let through, but molecules that cannot be kekulized fail.
const SANITIZE_OPS: u32 = SANITIZE_ALL ^ SANITIZE_CLEANUP ^ SANITIZE_PROPERTIES;

/// `sanitizeMol` that swallows `MolSanitizeException`s (logged); `Ok(false)`
/// when sanitization failed.
pub(crate) fn try_sanitize(mol: &mut RWMol, ops: u32) -> Result<bool> {
    match sanitize_mol(mol, ops) {
        Ok(()) => Ok(true),
        Err((_, RDKitError::Sanitize(_))) => {
            rd_log_info!("FAILED sanitizeMol.");
            Ok(false)
        }
        Err((_, e)) => Err(e),
    }
}

fn rule_name(rxn: &ChemicalReaction) -> String {
    rxn.get_prop(common_properties::_Name)
        .map(|v| v.to_rd_string())
        .unwrap_or_default()
}

/// `Normalizer`
#[derive(Debug, Clone)]
pub struct Normalizer {
    params: Arc<TransformCatalogParams>,
    max_restarts: u32,
}

impl Default for Normalizer {
    /// `Normalizer()`: the default transforms, 200 restarts.
    fn default() -> Self {
        rd_log_info!("Initializing Normalizer");
        Normalizer {
            params: Arc::new(TransformCatalogParams::defaults().clone()),
            max_restarts: 200,
        }
    }
}

impl Normalizer {
    /// `Normalizer(normalizeFile, maxRestarts)`; an empty file name selects
    /// the default transforms.
    pub fn from_file(file_name: &str, max_restarts: u32) -> Result<Self> {
        rd_log_info!("Initializing Normalizer");
        let params = if file_name.is_empty() {
            TransformCatalogParams::defaults().clone()
        } else {
            TransformCatalogParams::from_file(file_name)?
        };
        Ok(Self::with_params(params, max_restarts))
    }

    /// `Normalizer(std::istream &, maxRestarts)`
    pub fn from_text(text: &str, max_restarts: u32) -> Result<Self> {
        rd_log_info!("Initializing Normalizer");
        Ok(Self::with_params(
            TransformCatalogParams::from_text(text)?,
            max_restarts,
        ))
    }

    /// `Normalizer(normalizations, maxRestarts)`
    pub fn from_data(data: &[(String, String)], max_restarts: u32) -> Result<Self> {
        rd_log_info!("Initializing Normalizer");
        Ok(Self::with_params(
            TransformCatalogParams::from_data(data)?,
            max_restarts,
        ))
    }

    pub fn with_params(params: TransformCatalogParams, max_restarts: u32) -> Self {
        Normalizer {
            params: Arc::new(params),
            max_restarts,
        }
    }

    pub fn params(&self) -> &TransformCatalogParams {
        &self.params
    }

    /// `normalizeInPlace`: applies the first transform that changes the
    /// molecule (in place), then starts over, until none applies or
    /// `max_restarts` is reached.
    pub fn normalize_in_place(&self, mol: &mut RWMol) -> Result<()> {
        rd_log_info!("Running Normalizer");
        if mol.num_atoms() == 0 {
            return Ok(());
        }
        let transforms = self.params.transformations();
        for t in transforms {
            if t.num_product_templates() != 1
                || t.num_reactant_templates() != 1
                || t.products()[0].num_atoms() > t.reactants()[0].num_atoms()
            {
                return Err(RDKitError::Value(
                    "normalizeInPlace can only be used with transforms which have a single \
                     reactant and single product. The number of atoms in the product cannot be \
                     larger than the number of atoms in the reactant."
                        .into(),
                ));
            }
        }
        // we might want ring info
        if !mol.ring_info().is_symm_sssr() {
            symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
        }
        for _ in 0..self.max_restarts {
            let mut applied = false;
            for t in transforms {
                if run_reactant_in_place(t, mol, false)? {
                    rd_log_info!("Rule applied: {}", rule_name(t));
                    try_sanitize(mol, SANITIZE_OPS)?;
                    applied = true;
                    break;
                }
            }
            if !applied {
                return Ok(());
            }
        }
        rd_log_info!(
            "Gave up normalization after {} restarts.",
            self.max_restarts
        );
        Ok(())
    }

    /// `normalize`: normalizes each fragment separately (with
    /// `runReactants`) and recombines them.
    pub fn normalize(&self, mol: &ROMol) -> Result<ROMol> {
        rd_log_info!("Running Normalizer");
        if mol.num_atoms() == 0 {
            return Ok(mol.clone());
        }
        let transforms = self.params.transformations();
        let (frags, _) = get_mol_frags_as_mols(mol, false, true)?;
        let mut nfrags = Vec::with_capacity(frags.len());
        for mut frag in frags {
            frag.update_property_cache(false)?;
            nfrags.push(self.normalize_fragment(&frag, transforms)?);
        }
        let mut out = nfrags.pop().expect("at least one fragment");
        for nfrag in &nfrags {
            out = combine_mols(&out, nfrag, Point3D::default())?;
        }
        Ok(out)
    }

    fn normalize_fragment(&self, mol: &ROMol, transforms: &[ChemicalReaction]) -> Result<ROMol> {
        let mut nfrag = mol.clone();
        if !nfrag.ring_info().is_find_fast_or_better() {
            fast_find_rings(&mut nfrag)?;
        }
        let mut seen = BTreeSet::new();
        for _ in 0..self.max_restarts {
            let mut applied = false;
            for t in transforms {
                if let Some((smiles, product)) = apply_transform(&nfrag, t)? {
                    if !smiles.is_empty() && seen.insert(smiles) {
                        rd_log_info!("Rule applied: {}", rule_name(t));
                        nfrag = product;
                        applied = true;
                        break;
                    }
                }
            }
            if !applied {
                return Ok(nfrag);
            }
        }
        rd_log_info!(
            "Gave up normalization after {} restarts.",
            self.max_restarts
        );
        Ok(nfrag)
    }
}

/// `applyTransform`: applies `rule` repeatedly (at most 20 times), keeping
/// the product with the alphabetically first SMILES each time. `None` if the
/// rule does not apply; the SMILES is empty if it does not change `mol`.
fn apply_transform(mol: &ROMol, rule: &ChemicalReaction) -> Result<Option<(String, ROMol)>> {
    let mut current: (String, ROMol) = (String::new(), mol.clone());
    for i in 0..20 {
        let mut pdts: BTreeMap<String, ROMol> = BTreeMap::new();
        for pdt in run_reactants(rule, std::slice::from_ref(&current.1), 1000, 1)? {
            let Some(mut tmol) = pdt.into_iter().next() else {
                continue;
            };
            if try_sanitize(&mut tmol, SANITIZE_OPS)? {
                let smiles = mol_to_smiles(&mut tmol, &SmilesWriteParams::default())?;
                pdts.insert(smiles, tmol);
            }
        }
        match pdts.into_iter().next() {
            Some(first) => current = first,
            None if i > 0 => return Ok(Some(current)),
            None => return Ok(None),
        }
    }
    Ok(Some(current))
}

/// `normalizerFromParams`
pub fn normalizer_from_params(params: &CleanupParameters) -> Result<Normalizer> {
    let max_restarts = params.max_restarts as u32;
    if params.normalization_data.is_empty() {
        Normalizer::from_file(&params.normalizations, max_restarts)
    } else {
        Normalizer::from_data(&params.normalization_data, max_restarts)
    }
}
