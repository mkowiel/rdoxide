//! Storage back-ends of the [`SubstructLibrary`](super::SubstructLibrary):
//! molecule holders, the pattern-fingerprint holder and the key holder
//! (upstream `SubstructLibrary.h`).

use std::sync::Arc;

use crate::datastructs::ExplicitBitVect;
use crate::datastructs::bit_ops::all_probe_bits_match;
use crate::error::{RDKitError, Result};
use crate::graphmol::fingerprints::pattern::pattern_fingerprint_mol;
use crate::graphmol::pickle::{get_default_pickle_properties, mol_from_pickle, pickle_mol};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::rdgeneral::RDProps;

fn index_error(idx: u32) -> RDKitError {
    RDKitError::Index(idx.to_string())
}

/// The molecule storage of a library: the four upstream holders, which trade
/// memory for the cost of rebuilding a molecule on every access.
#[derive(Debug, Clone)]
pub enum MolHolder {
    /// `MolHolder`: the molecules themselves (fastest, most memory).
    Mols(Vec<Arc<ROMol>>),
    /// `CachedMolHolder`: molecule pickles.
    Binary(Vec<Vec<u8>>),
    /// `CachedSmilesMolHolder`: isomeric SMILES, fully sanitized on access.
    Smiles(Vec<String>),
    /// `CachedTrustedSmilesMolHolder`: SMILES written by this library, parsed
    /// without sanitization (rings are perceived on demand).
    TrustedSmiles(Vec<String>),
}

impl Default for MolHolder {
    fn default() -> Self {
        MolHolder::Mols(Vec::new())
    }
}

impl MolHolder {
    /// Adds a molecule and returns its index.
    pub fn add_mol(&mut self, mol: &ROMol) -> Result<u32> {
        match self {
            MolHolder::Mols(v) => v.push(Arc::new(mol.clone())),
            MolHolder::Binary(v) => v.push(pickle_mol(mol, get_default_pickle_properties())?),
            MolHolder::Smiles(v) | MolHolder::TrustedSmiles(v) => v.push(mol_to_smiles(
                &mut mol.clone(),
                &SmilesWriteParams::default(),
            )?),
        }
        Ok(self.len() - 1)
    }

    /// `addBinary`: adds a molecule pickle without validating it.
    pub fn add_binary(&mut self, pickle: Vec<u8>) -> Result<u32> {
        match self {
            MolHolder::Binary(v) => v.push(pickle),
            _ => return Err(RDKitError::Value("not a binary molecule holder".into())),
        }
        Ok(self.len() - 1)
    }

    /// `addSmiles`: adds a SMILES string without validating it.
    pub fn add_smiles(&mut self, smiles: &str) -> Result<u32> {
        match self {
            MolHolder::Smiles(v) | MolHolder::TrustedSmiles(v) => v.push(smiles.to_string()),
            _ => return Err(RDKitError::Value("not a SMILES molecule holder".into())),
        }
        Ok(self.len() - 1)
    }

    /// The molecule at `idx`; `None` for a stored SMILES that does not parse
    /// (such entries never match). Out of range is an `IndexError`.
    pub fn get_mol(&self, idx: u32) -> Result<Option<Arc<ROMol>>> {
        let i = idx as usize;
        if i >= self.len() as usize {
            return Err(index_error(idx));
        }
        Ok(match self {
            MolHolder::Mols(v) => Some(Arc::clone(&v[i])),
            MolHolder::Binary(v) => Some(Arc::new(mol_from_pickle(
                &v[i],
                get_default_pickle_properties(),
            )?)),
            // a SMILES that fails sanitization is as unusable as one that
            // does not parse
            MolHolder::Smiles(v) => mol_from_smiles(&v[i], &SmilesParserParams::default())
                .ok()
                .flatten()
                .map(Arc::new),
            MolHolder::TrustedSmiles(v) => {
                let params = SmilesParserParams {
                    sanitize: false,
                    ..Default::default()
                };
                match mol_from_smiles(&v[i], &params)? {
                    Some(mut m) => {
                        m.update_property_cache(true)?;
                        Some(Arc::new(m))
                    }
                    None => None,
                }
            }
        })
    }

    pub fn len(&self) -> u32 {
        (match self {
            MolHolder::Mols(v) => v.len(),
            MolHolder::Binary(v) => v.len(),
            MolHolder::Smiles(v) | MolHolder::TrustedSmiles(v) => v.len(),
        }) as u32
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// `PatternHolder` and `TautomerPatternHolder`: one pattern fingerprint per
/// molecule, used to rule out impossible matches.
#[derive(Debug, Clone)]
pub struct PatternHolder {
    num_bits: u32,
    tautomeric: bool,
    fps: Vec<ExplicitBitVect>,
}

impl Default for PatternHolder {
    fn default() -> Self {
        Self::new(Self::DEFAULT_NUM_BITS)
    }
}

impl PatternHolder {
    pub const DEFAULT_NUM_BITS: u32 = 2048;

    pub fn new(num_bits: u32) -> Self {
        PatternHolder {
            num_bits,
            tautomeric: false,
            fps: Vec::new(),
        }
    }

    /// `TautomerPatternHolder`: fingerprints that are insensitive to tautomer
    /// bond changes.
    pub fn new_tautomeric(num_bits: u32) -> Self {
        PatternHolder {
            tautomeric: true,
            ..Self::new(num_bits)
        }
    }

    pub fn num_bits(&self) -> u32 {
        self.num_bits
    }

    pub fn set_num_bits(&mut self, num_bits: u32) {
        self.num_bits = num_bits;
    }

    pub fn is_tautomeric(&self) -> bool {
        self.tautomeric
    }

    /// `makeFingerprint`
    pub fn make_fingerprint(&self, mol: &ROMol) -> Result<ExplicitBitVect> {
        pattern_fingerprint_mol(mol, self.num_bits, self.tautomeric)
    }

    /// Adds the fingerprint of `mol` and returns its index.
    pub fn add_mol(&mut self, mol: &ROMol) -> Result<u32> {
        let fp = self.make_fingerprint(mol)?;
        Ok(self.add_fingerprint(fp))
    }

    /// Adds a fingerprint made elsewhere; it must be compatible with
    /// [`make_fingerprint`](Self::make_fingerprint).
    pub fn add_fingerprint(&mut self, fp: ExplicitBitVect) -> u32 {
        self.fps.push(fp);
        (self.fps.len() - 1) as u32
    }

    /// False if a substructure search can never match the molecule at `idx`.
    pub fn passes_filter(&self, idx: u32, query: &ExplicitBitVect) -> Result<bool> {
        all_probe_bits_match(query, self.get_fingerprint(idx)?)
    }

    pub fn get_fingerprint(&self, idx: u32) -> Result<&ExplicitBitVect> {
        self.fps.get(idx as usize).ok_or_else(|| index_error(idx))
    }

    pub fn fingerprints(&self) -> &[ExplicitBitVect] {
        &self.fps
    }

    pub fn fingerprints_mut(&mut self) -> &mut Vec<ExplicitBitVect> {
        &mut self.fps
    }

    pub fn len(&self) -> u32 {
        self.fps.len() as u32
    }

    pub fn is_empty(&self) -> bool {
        self.fps.is_empty()
    }
}

/// `KeyFromPropHolder`: user keys, read from a molecule property.
#[derive(Debug, Clone)]
pub struct KeyFromPropHolder {
    propname: String,
    keys: Vec<String>,
}

impl Default for KeyFromPropHolder {
    fn default() -> Self {
        Self::new("_Name")
    }
}

impl KeyFromPropHolder {
    pub fn new(propname: &str) -> Self {
        KeyFromPropHolder {
            propname: propname.to_string(),
            keys: Vec::new(),
        }
    }

    pub fn prop_name(&self) -> &str {
        &self.propname
    }

    pub fn set_prop_name(&mut self, name: &str) {
        self.propname = name.to_string();
    }

    /// Adds the key of `mol`: its property value, or `LIBIDX-<index>` when it
    /// has none.
    pub fn add_mol(&mut self, mol: &ROMol) -> u32 {
        let key = mol
            .get_prop(&self.propname)
            .map(|v| v.to_rd_string())
            .unwrap_or_else(|| format!("LIBIDX-{}", self.keys.len()));
        self.add_key(key)
    }

    /// Adds a key; keys must be added in molecule order.
    pub fn add_key(&mut self, key: String) -> u32 {
        self.keys.push(key);
        (self.keys.len() - 1) as u32
    }

    pub fn get_key(&self, idx: u32) -> Result<&str> {
        self.keys
            .get(idx as usize)
            .map(String::as_str)
            .ok_or_else(|| index_error(idx))
    }

    pub fn get_keys(&self, indices: &[u32]) -> Result<Vec<String>> {
        indices
            .iter()
            .map(|&i| self.get_key(i).map(str::to_string))
            .collect()
    }

    pub fn keys(&self) -> &[String] {
        &self.keys
    }

    pub fn keys_mut(&mut self) -> &mut Vec<String> {
        &mut self.keys
    }

    pub fn len(&self) -> u32 {
        self.keys.len() as u32
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }
}
