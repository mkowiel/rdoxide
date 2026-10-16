//! Port of RDKit `Code/GraphMol/MolBundle.h`: groups of related molecules.
//! Boost-serialization pickles are not supported.

use std::sync::Arc;

use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::rdgeneral::{Dict, RDProps};
use crate::{RDKitError, Result};

/// `MolBundle`; with `fixed_size` it is a `FixedMolSizeMolBundle`, whose
/// molecules all have the same numbers of atoms and bonds.
#[derive(Debug, Clone, Default)]
pub struct MolBundle {
    props: Dict,
    mols: Vec<Arc<ROMol>>,
    fixed_size: bool,
}

impl RDProps for MolBundle {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

impl MolBundle {
    pub fn new() -> Self {
        Self::default()
    }

    /// `FixedMolSizeMolBundle`
    pub fn new_fixed_size() -> Self {
        MolBundle {
            fixed_size: true,
            ..Default::default()
        }
    }

    pub fn mols(&self) -> &[Arc<ROMol>] {
        &self.mols
    }

    /// `addMol`: returns the new size.
    pub fn add_mol(&mut self, mol: Arc<ROMol>) -> Result<usize> {
        if self.fixed_size
            && let Some(first) = self.mols.first()
        {
            if mol.num_atoms() != first.num_atoms() {
                return Err(RDKitError::Value(
                    "all molecules in a bundle must have the same number of atoms".into(),
                ));
            }
            if mol.num_bonds() != first.num_bonds() {
                return Err(RDKitError::Value(
                    "all molecules in a bundle must have the same number of bonds".into(),
                ));
            }
        }
        self.mols.push(mol);
        Ok(self.mols.len())
    }

    pub fn size(&self) -> usize {
        self.mols.len()
    }

    pub fn is_empty(&self) -> bool {
        self.mols.is_empty()
    }

    /// `getMol`
    pub fn mol(&self, idx: usize) -> Result<&Arc<ROMol>> {
        self.mols
            .get(idx)
            .ok_or_else(|| RDKitError::Index(format!("{idx}")))
    }
}

/// `SubstructMatch(bundle, query, params)`: the matches of the first molecule
/// of the bundle that has any.
pub fn substruct_match_bundle(
    bundle: &MolBundle,
    query: &ROMol,
    ps: &SubstructMatchParameters,
) -> Result<Vec<MatchVect>> {
    for mol in bundle.mols() {
        let res = substruct_match(mol, query, ps)?;
        if !res.is_empty() {
            return Ok(res);
        }
    }
    Ok(Vec::new())
}

/// `SubstructMatch(mol, query bundle, params)`: the matches of the first
/// query that has any.
pub fn substruct_match_query_bundle(
    mol: &ROMol,
    query: &MolBundle,
    ps: &SubstructMatchParameters,
) -> Result<Vec<MatchVect>> {
    for q in query.mols() {
        let res = substruct_match(mol, q, ps)?;
        if !res.is_empty() {
            return Ok(res);
        }
    }
    Ok(Vec::new())
}

/// `SubstructMatch(bundle, query bundle, params)`
pub fn substruct_match_bundles(
    bundle: &MolBundle,
    query: &MolBundle,
    ps: &SubstructMatchParameters,
) -> Result<Vec<MatchVect>> {
    for mol in bundle.mols() {
        let res = substruct_match_query_bundle(mol, query, ps)?;
        if !res.is_empty() {
            return Ok(res);
        }
    }
    Ok(Vec::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    fn mol(smi: &str) -> Arc<ROMol> {
        Arc::new(
            mol_from_smiles(smi, &SmilesParserParams::default())
                .unwrap()
                .unwrap(),
        )
    }

    #[test]
    fn bundle_matching() {
        let mut b = MolBundle::new();
        for smi in ["Fc1c(Cl)cccc1", "Fc1cc(Cl)ccc1", "Fc1ccc(Cl)cc1"] {
            b.add_mol(mol(smi)).unwrap();
        }
        let ps = SubstructMatchParameters::default();
        assert_eq!(
            substruct_match_query_bundle(&mol("Fc1c(Cl)cccc1C"), &b, &ps)
                .unwrap()
                .len(),
            1
        );
        assert!(
            substruct_match_query_bundle(&mol("Fc1c(Br)cccc1"), &b, &ps)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            substruct_match_bundle(&b, &mol("ClC"), &ps).unwrap().len(),
            1
        );
        assert!(!substruct_match_bundles(&b, &b, &ps).unwrap().is_empty());
        let mut fixed = MolBundle::new_fixed_size();
        fixed.add_mol(mol("CCC")).unwrap();
        assert!(fixed.add_mol(mol("CCCC")).is_err());
    }
}
