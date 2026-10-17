//! Ported from RDKit `External/INCHI-API/inchi.h` and `inchi.cpp`: the RDKit
//! side of InChI support (molecule <-> InChI conversion on top of the InChI
//! library, [`crate::external::inchi`]).
//!
//! Copyright (c) 2011-2025 Novartis Institutes for BioMedical Research Inc.
//! and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-ChemReactions-Novartis-BSD-3-Clause.txt, the same Novartis
//! license text).
//!
//! Status: [`inchi_to_inchi_key`] and [`get_inchi_version`] are complete.
//! [`mol_to_inchi`], [`mol_block_to_inchi`], [`inchi_to_mol`] and
//! [`mol_to_inchi_key`] need the InChI generation / structure restoration
//! pipelines, which are not ported yet; they return
//! [`RDKitError::NotImplemented`] (see `docs/INCHI_PORT.md`).

use crate::external::inchi::{self, InchiKeyError};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::{RDKitError, Result};

/// `ExtraInchiReturnValues`: the extra outputs of the InChI API calls.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExtraInchiReturnValues {
    /// `returnCode`: an `inchi_Ret_*` value ([`inchi::InchiRet`]).
    pub return_code: i32,
    /// `messagePtr`
    pub message: String,
    /// `logPtr`
    pub log: String,
    /// `auxInfoPtr` (not used for `InchiToMol`)
    pub aux_info: String,
}

fn not_ported(what: &str) -> RDKitError {
    RDKitError::NotImplemented(format!(
        "{what} is not available in rdoxide yet: the InChI library port \
         (InChI {}) only covers InChIKey generation so far",
        inchi::CURRENT_VER
    ))
}

/// `InchiToMol`: a molecule from a (standard or non-standard) InChI string.
pub fn inchi_to_mol(
    _inchi: &str,
    _rv: &mut ExtraInchiReturnValues,
    _sanitize: bool,
    _remove_hs: bool,
) -> Result<Option<RWMol>> {
    Err(not_ported("InchiToMol"))
}

/// `MolToInchi`: the InChI of a molecule. `options` is the space-delimited
/// InChI option string passed through to the InChI API (`-` and `/`
/// prefixes are both accepted).
pub fn mol_to_inchi(
    _mol: &ROMol,
    _rv: &mut ExtraInchiReturnValues,
    _options: Option<&str>,
) -> Result<String> {
    Err(not_ported("MolToInchi"))
}

/// `MolBlockToInchi`: the InChI of a mol block, parsed by the InChI library's
/// own MOL reader (`MakeINCHIFromMolfileText`).
pub fn mol_block_to_inchi(
    _mol_block: &str,
    _rv: &mut ExtraInchiReturnValues,
    _options: Option<&str>,
) -> Result<String> {
    Err(not_ported("MolBlockToInchi"))
}

/// `InchiToInchiKey`: the InChIKey of an InChI string; an empty string (with
/// an error logged) when the InChI is not valid.
pub fn inchi_to_inchi_key(inchi_str: &str) -> String {
    let error = match inchi::get_inchikey_from_inchi(inchi_str, false, false) {
        Ok(k) => return k.key,
        Err(InchiKeyError::UnknownError) => "Unknown error",
        Err(InchiKeyError::EmptyInput) => "Empty input",
        Err(InchiKeyError::InvalidInchiPrefix) => "Invalid InChI prefix",
        Err(InchiKeyError::NotEnoughMemory) => "Not enough memory",
        Err(InchiKeyError::InvalidInchi) => "Invalid input InChI string",
        Err(InchiKeyError::InvalidStdInchi) => "Invalid standard InChI string",
    };
    crate::rd_log_error!("{error} in generating InChI Key\n");
    String::new()
}

/// `MolToInchiKey`: `InchiToInchiKey(MolToInchi(mol, rv, options))`.
pub fn mol_to_inchi_key(mol: &ROMol, options: Option<&str>) -> Result<String> {
    let mut rv = ExtraInchiReturnValues::default();
    Ok(inchi_to_inchi_key(&mol_to_inchi(mol, &mut rv, options)?))
}

/// `getInchiVersion`: the version of the InChI software being used.
pub fn get_inchi_version() -> &'static str {
    inchi::CURRENT_VER
}
