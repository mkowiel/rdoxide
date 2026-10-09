//! Ported from RDKit `Code/GraphMol/MolChemicalFeatures/` (`MolChemicalFeature`,
//! `MolChemicalFeatureDef`, `MolChemicalFeatureFactory`): feature definitions
//! (SMARTS patterns with a family, a type and atom weights) and the factory
//! that finds the features of a molecule.

pub mod parser;

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::sync::Arc;

use crate::geometry::Point3D;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

const FEATURE_FILE: &str = "Code/GraphMol/MolChemicalFeatures/MolChemicalFeature.cpp";

/// `MolChemicalFeatureDef`
#[derive(Debug, Clone)]
pub struct MolChemicalFeatureDef {
    family: String,
    feature_type: String,
    smarts: String,
    pattern: ROMol,
    weights: Vec<f64>,
}

impl MolChemicalFeatureDef {
    /// The pattern must be valid SMARTS.
    pub fn new(smarts: &str, family: &str, feature_type: &str) -> Result<Self> {
        let pattern = mol_from_smarts(smarts, &SmartsParserParams::default())?
            .ok_or_else(|| RDKitError::Value(format!("invalid SMARTS: {smarts}")))?;
        Ok(MolChemicalFeatureDef {
            family: family.to_string(),
            feature_type: feature_type.to_string(),
            smarts: smarts.to_string(),
            pattern,
            weights: Vec::new(),
        })
    }

    pub fn family(&self) -> &str {
        &self.family
    }
    pub fn feature_type(&self) -> &str {
        &self.feature_type
    }
    pub fn smarts(&self) -> &str {
        &self.smarts
    }
    pub fn pattern(&self) -> &ROMol {
        &self.pattern
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// `setWeights`: one weight per pattern atom.
    pub fn set_weights(&mut self, weights: &[f64]) -> Result<()> {
        self.weights.splice(0..0, weights.iter().copied());
        if self.weights.len() != self.pattern.num_atoms() as usize {
            return Err(RDKitError::Value(
                " pattern->getNumAtoms() != len(feature weight vector)".into(),
            ));
        }
        Ok(())
    }

    /// `normalizeWeights`: scales the weights to sum to one.
    pub fn normalize_weights(&mut self) {
        let total: f64 = self.weights.iter().sum();
        for w in &mut self.weights {
            *w /= total;
        }
    }
}

/// `MolChemicalFeature`: a feature found in a molecule, as the molecule's
/// atoms (in pattern order) plus the definition it matched. The molecule is
/// passed to the position queries instead of being referenced.
#[derive(Debug, Clone)]
pub struct MolChemicalFeature {
    id: i32,
    def: Arc<MolChemicalFeatureDef>,
    atoms: Vec<u32>,
    active_conf: i32,
    locs: BTreeMap<i32, Point3D>,
}

impl MolChemicalFeature {
    pub fn id(&self) -> i32 {
        self.id
    }
    pub fn family(&self) -> &str {
        self.def.family()
    }
    pub fn feature_type(&self) -> &str {
        self.def.feature_type()
    }
    pub fn def(&self) -> &MolChemicalFeatureDef {
        &self.def
    }
    /// The atom indices defining the feature, in pattern order.
    pub fn atoms(&self) -> &[u32] {
        &self.atoms
    }
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
    pub fn active_conformer(&self) -> i32 {
        self.active_conf
    }
    pub fn set_active_conformer(&mut self, conf_id: i32) {
        self.active_conf = conf_id;
    }
    /// `clearCache`: forgets the positions computed so far.
    pub fn clear_cache(&mut self) {
        self.locs.clear();
    }

    /// `getPos()`: the position in the active conformer.
    pub fn pos(&mut self, mol: &ROMol) -> Result<Point3D> {
        self.pos_in(mol, self.active_conf)
    }

    /// `getPos(confId)`: the weighted centre of the feature's atoms in a
    /// conformer (`-1`: the first one). Positions are cached per conformer
    /// until [`clear_cache`](Self::clear_cache).
    pub fn pos_in(&mut self, mol: &ROMol, conf_id: i32) -> Result<Point3D> {
        if mol.num_conformers() == 0 {
            return Err(RDKitError::precondition_at(
                "molecule has no conformers",
                "dp_mol->getNumConformers()",
                FEATURE_FILE,
                36,
            ));
        }
        let conf_id = if conf_id == -1 {
            mol.conformers()[0].id() as i32
        } else {
            conf_id
        };
        if let Some(p) = self.locs.get(&conf_id) {
            return Ok(*p);
        }
        let conf = mol.conformer(conf_id)?;
        let pos = if self.atoms.len() == 1 {
            conf.atom_pos(self.atoms[0])?
        } else {
            let weights = self.def.weights();
            if weights.len() != self.atoms.len() {
                return Err(RDKitError::precondition_at(
                    "weight/atom mismatch",
                    "dp_def->getNumWeights() == this->getNumAtoms()",
                    FEATURE_FILE,
                    57,
                ));
            }
            let mut res = Point3D::new(0.0, 0.0, 0.0);
            for (&atom, &w) in self.atoms.iter().zip(weights) {
                res += conf.atom_pos(atom)? * w;
            }
            res
        };
        self.locs.insert(conf_id, pos);
        Ok(pos)
    }
}

/// `MolChemicalFeatureFactory`
#[derive(Debug, Clone, Default)]
pub struct MolChemicalFeatureFactory {
    defs: Vec<Arc<MolChemicalFeatureDef>>,
}

impl MolChemicalFeatureFactory {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn feature_defs(&self) -> &[Arc<MolChemicalFeatureDef>] {
        &self.defs
    }
    pub fn num_feature_defs(&self) -> usize {
        self.defs.len()
    }
    pub fn add_feature_def(&mut self, def: Arc<MolChemicalFeatureDef>) {
        self.defs.push(def);
    }

    /// `getFeaturesForMol` on one thread.
    pub fn features_for_mol(
        &self,
        mol: &ROMol,
        include_only: &str,
        conf_id: i32,
    ) -> Result<Vec<MolChemicalFeature>> {
        self.features_for_mol_threaded(mol, include_only, conf_id, 1)
    }

    /// `getFeaturesForMol`: the features of `mol` (limited to one family when
    /// `include_only` is not empty), with the matches of the definitions
    /// searched on up to `num_threads` threads (`<= 0`: all but `|n|`
    /// hardware threads). A match whose atoms are contained in an earlier
    /// match of the same family is dropped, so the result does not depend
    /// on the thread count.
    pub fn features_for_mol_threaded(
        &self,
        mol: &ROMol,
        include_only: &str,
        conf_id: i32,
        num_threads: i32,
    ) -> Result<Vec<MolChemicalFeature>> {
        let wanted: Vec<&Arc<MolChemicalFeatureDef>> = self
            .defs
            .iter()
            .filter(|d| include_only.is_empty() || d.family() == include_only)
            .collect();
        let params = SubstructMatchParameters::default();
        let all_matches = par_map(&wanted, num_threads, |d| {
            substruct_match(mol, d.pattern(), &params)
        });

        let mut seen: Vec<(&str, BTreeSet<i32>)> = Vec::new();
        let mut res = Vec::new();
        for (def, matches) in wanted.into_iter().zip(all_matches) {
            for m in matches? {
                let atom_set: BTreeSet<i32> = m.iter().map(|&(_, a)| a).collect();
                let covered = seen
                    .iter()
                    .any(|(family, set)| *family == def.family() && atom_set.is_subset(set));
                if covered {
                    continue;
                }
                seen.push((def.family(), atom_set));
                let mut atoms = vec![0u32; m.len()];
                for &(query_idx, atom_idx) in &m {
                    atoms[query_idx as usize] = atom_idx as u32;
                }
                res.push(MolChemicalFeature {
                    id: res.len() as i32 + 1,
                    def: Arc::clone(def),
                    atoms,
                    active_conf: conf_id,
                    locs: BTreeMap::new(),
                });
            }
        }
        Ok(res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    #[test]
    fn weights_must_match_the_pattern() {
        let mut def = MolChemicalFeatureDef::new("[N,O;H1]", "Donor", "d").unwrap();
        assert!(def.set_weights(&[1.0, 2.0]).is_err());
        let mut def = MolChemicalFeatureDef::new("[N,O;H1]", "Donor", "d").unwrap();
        def.set_weights(&[4.0]).unwrap();
        def.normalize_weights();
        assert_eq!(def.weights(), &[1.0]);
    }

    #[test]
    fn contained_matches_of_a_family_are_dropped() {
        let mol = mol_from_smiles("OCC(=O)N", &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        let mut factory = MolChemicalFeatureFactory::new();
        factory.add_feature_def(Arc::new(
            MolChemicalFeatureDef::new("[O,N;!H0]", "Donor", "any").unwrap(),
        ));
        factory.add_feature_def(Arc::new(
            MolChemicalFeatureDef::new("[O;!H0]", "Donor", "alcohol").unwrap(),
        ));
        let one = factory.features_for_mol(&mol, "", -1).unwrap();
        let four = factory.features_for_mol_threaded(&mol, "", -1, 4).unwrap();
        assert_eq!(one.len(), 2);
        assert_eq!(
            one.iter().map(|f| f.atoms().to_vec()).collect::<Vec<_>>(),
            four.iter().map(|f| f.atoms().to_vec()).collect::<Vec<_>>()
        );
        assert!(
            factory
                .features_for_mol(&mol, "Acceptor", -1)
                .unwrap()
                .is_empty()
        );
    }
}
