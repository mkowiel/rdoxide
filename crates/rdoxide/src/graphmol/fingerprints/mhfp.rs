//! MinHash fingerprints (MHFP) and the SECFP binary fingerprint.
//!
//! Ported from RDKit `Code/GraphMol/Fingerprints/MHFP.{h,cpp}` (2026.09.1);
//! the method is described in Probst and Reymond, J. Cheminform. 2018, 10, 66.
//!
//! Copyright (C) 2019 Daniel Probst and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::graphmol::kekulize::{DEFAULT_MAX_BACKTRACKS, kekulize};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, get_atom_smiles, mol_to_smiles};
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::graphmol::subgraphs::{find_atom_environment_of_radius_n, path_to_submol};
use crate::rdgeneral::random::Mt19937;
use crate::{RDKitError, Result};

pub const MHFP_FINGERPRINT_VERSION: &str = "1.0.0";

const MAX_HASH: u32 = u32::MAX;
const PRIME: u64 = 2_305_843_009_213_693_951;

/// The Fowler-Noll-Vo hash (FNV-1) of a string.
pub fn fnv_hash(s: &str) -> u32 {
    const FNV_PRIME: u32 = 0x0100_0193;
    s.bytes().fold(0x811C_9DC5u32, |hash, b| {
        (hash ^ u32::from(b)).wrapping_mul(FNV_PRIME)
    })
}

/// Parameters of the shingling shared by all the `encode*` methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShinglingParams {
    /// the maximum radius of the substructures extracted at each atom
    pub radius: u32,
    /// add the SSSR rings to the shingling
    pub rings: bool,
    /// write isomeric SMILES
    pub isomeric: bool,
    /// kekulize the molecule first (fails if that is impossible)
    pub kekulize: bool,
    /// the minimum radius; 0 also adds the atom SMILES
    pub min_radius: u32,
}

impl Default for ShinglingParams {
    fn default() -> Self {
        ShinglingParams {
            radius: 3,
            rings: true,
            isomeric: false,
            kekulize: false,
            min_radius: 1,
        }
    }
}

/// `MHFPEncoder`: a MinHash encoder with a fixed set of permutations.
#[derive(Debug, Clone)]
pub struct MhfpEncoder {
    perms_a: Vec<u32>,
    perms_b: Vec<u32>,
}

impl Default for MhfpEncoder {
    fn default() -> Self {
        Self::new(2048, 42)
    }
}

impl MhfpEncoder {
    /// Hash functions (`n_permutations` of them) are drawn from `seed`;
    /// minhashes made with different settings are not comparable.
    pub fn new(n_permutations: usize, seed: u32) -> Self {
        let mut rng = Mt19937::new(seed);
        // note: the lookups below deliberately scan the whole, zero-filled,
        // vectors as upstream does
        let mut perms_a = vec![0u32; n_permutations];
        let mut perms_b = vec![0u32; n_permutations];
        for i in 0..n_permutations {
            let mut a = rng.uniform_u32(1, MAX_HASH);
            let mut b = rng.uniform_u32(0, MAX_HASH);
            while perms_a.contains(&a) {
                a = rng.uniform_u32(1, MAX_HASH);
            }
            while perms_b.contains(&b) {
                b = rng.uniform_u32(1, MAX_HASH);
            }
            perms_a[i] = a;
            perms_b[i] = b;
        }
        MhfpEncoder { perms_a, perms_b }
    }

    fn min_hash(&self, hashes: impl Iterator<Item = u32>) -> Vec<u32> {
        let mut mh = vec![MAX_HASH; self.perms_a.len()];
        for hash in hashes {
            for ((slot, &a), &b) in mh.iter_mut().zip(&self.perms_a).zip(&self.perms_b) {
                let tmp = a.wrapping_mul(hash).wrapping_add(b);
                *slot = (*slot).min(fast_mod(u64::from(tmp), PRIME) as u32 & MAX_HASH);
            }
        }
        mh
    }

    /// `FromStringArray`: the MinHash of arbitrary strings.
    pub fn from_string_array<S: AsRef<str>>(&self, strings: &[S]) -> Vec<u32> {
        self.min_hash(strings.iter().map(|s| fnv_hash(s.as_ref())))
    }

    /// `FromArray`: the MinHash of unsigned integers.
    pub fn from_array(&self, values: &[u32]) -> Vec<u32> {
        self.min_hash(values.iter().copied())
    }

    /// `CreateShingling`: the circular n-grams (as SMILES) of a molecule.
    pub fn create_shingling(&self, mol: &ROMol, params: &ShinglingParams) -> Result<Vec<String>> {
        let mut tmol = mol.clone();
        if params.kekulize {
            kekulize(&mut tmol, true, false, DEFAULT_MAX_BACKTRACKS)?;
        }
        let mut shingling = Vec::new();
        if params.rings {
            for ring in tmol.ring_info().bond_rings().to_vec() {
                let (submol, _) = path_to_submol(&tmol, &ring, false)?;
                shingling.push(mol_to_smiles(
                    &mut submol.clone(),
                    &SmilesWriteParams::default(),
                )?);
            }
        }
        let mut min_radius = params.min_radius;
        if min_radius == 0 {
            let atom_params = SmilesWriteParams {
                do_isomeric_smiles: true,
                ..SmilesWriteParams::default()
            };
            for idx in 0..tmol.num_atoms() {
                shingling.push(get_atom_smiles(&tmol, idx, &atom_params)?);
            }
            min_radius += 1;
        }
        for index in 0..tmol.num_atoms() {
            for r in min_radius..=params.radius {
                let (path, _) = find_atom_environment_of_radius_n(&tmol, r, index, false, true)?;
                let (mut submol, atom_map) = path_to_submol(&tmol, &path, false)?;
                let Some(&rooted) = atom_map.get(&index) else {
                    continue;
                };
                let smiles_params = SmilesWriteParams {
                    do_isomeric_smiles: params.isomeric,
                    do_kekule: params.kekulize,
                    rooted_at_atom: rooted as i32,
                    ..SmilesWriteParams::default()
                };
                let smiles = mol_to_smiles(&mut submol, &smiles_params)?;
                if !smiles.is_empty() {
                    shingling.push(smiles);
                }
            }
        }
        Ok(shingling)
    }

    /// `CreateShingling(smiles, ...)`
    pub fn create_shingling_from_smiles(
        &self,
        smiles: &str,
        params: &ShinglingParams,
    ) -> Result<Vec<String>> {
        self.create_shingling(&parse_smiles(smiles)?, params)
    }

    /// `Encode(mol, ...)`: the MHFP vector of a molecule.
    pub fn encode(&self, mol: &ROMol, params: &ShinglingParams) -> Result<Vec<u32>> {
        Ok(self.from_string_array(&self.create_shingling(mol, params)?))
    }

    /// `Encode(smiles, ...)`
    pub fn encode_smiles(&self, smiles: &str, params: &ShinglingParams) -> Result<Vec<u32>> {
        Ok(self.from_string_array(&self.create_shingling_from_smiles(smiles, params)?))
    }

    /// `EncodeSECFP(mol, ..., length)`: the shingles hashed and folded into a bit vector.
    pub fn encode_secfp(
        &self,
        mol: &ROMol,
        params: &ShinglingParams,
        length: u32,
    ) -> Result<ExplicitBitVect> {
        fold(&self.create_shingling(mol, params)?, length)
    }

    /// `EncodeSECFP(smiles, ..., length)`
    pub fn encode_secfp_smiles(
        &self,
        smiles: &str,
        params: &ShinglingParams,
        length: u32,
    ) -> Result<ExplicitBitVect> {
        fold(&self.create_shingling_from_smiles(smiles, params)?, length)
    }

    /// `Distance`: the fraction of positions at which two MHFP vectors differ.
    pub fn distance(a: &[u32], b: &[u32]) -> f64 {
        let mismatches = a.iter().zip(b).filter(|(x, y)| x != y).count();
        mismatches as f64 / a.len() as f64
    }
}

fn fast_mod(input: u64, ceil: u64) -> u64 {
    if input >= ceil { input % ceil } else { input }
}

fn parse_smiles(smiles: &str) -> Result<ROMol> {
    mol_from_smiles(smiles, &SmilesParserParams::default())?
        .ok_or_else(|| RDKitError::precondition("could not parse smiles"))
}

fn fold(shingling: &[String], length: u32) -> Result<ExplicitBitVect> {
    let mut ebv = ExplicitBitVect::new(length, false);
    for shingle in shingling {
        ebv.set_bit(fnv_hash(shingle) % length)?;
    }
    Ok(ebv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv_of_empty_string_is_the_offset_basis() {
        assert_eq!(fnv_hash(""), 0x811C_9DC5);
    }

    #[test]
    fn the_same_seed_gives_the_same_permutations() {
        let (a, b) = (MhfpEncoder::new(16, 42), MhfpEncoder::new(16, 42));
        assert_eq!(a.from_array(&[1, 2, 3]), b.from_array(&[1, 2, 3]));
        assert_ne!(
            a.from_array(&[1, 2, 3]),
            MhfpEncoder::new(16, 43).from_array(&[1, 2, 3])
        );
    }

    #[test]
    fn distance_counts_mismatches() {
        assert_eq!(MhfpEncoder::distance(&[1, 2, 3, 4], &[1, 0, 3, 0]), 0.5);
    }
}
