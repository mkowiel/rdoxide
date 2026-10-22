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
//! Status: [`inchi_to_inchi_key`], [`get_inchi_version`], [`mol_to_inchi`],
//! [`mol_block_to_inchi`] and [`mol_to_inchi_key`] are complete (structures
//! without polymer data). [`inchi_to_mol`] needs the structure restoration
//! pipeline, which is not ported yet; it returns
//! [`RDKitError::NotImplemented`] (see `docs/INCHI_PORT.md`).

use crate::external::inchi::inchi_dll::{InchiInput as InchiLibInput, get_inchi};
use crate::external::inchi::inchi_dll_b::{InchiOutput, make_inchi_from_molfile_text};
use crate::external::inchi::{self, InchiKeyError};
use crate::external::inchi_input::{fix_option_symbol, mol_to_inchi_input};
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
    mol: &ROMol,
    rv: &mut ExtraInchiReturnValues,
    options: Option<&str>,
) -> Result<String> {
    let Some(input) = mol_to_inchi_input(mol, options)? else {
        return Ok(String::new());
    };
    let (code, out) = get_inchi(&InchiLibInput {
        atom: input.atoms,
        stereo_0d: input.stereo_0d,
        options: input.options,
    });
    rv.return_code = code;
    fill_extra(rv, &out);
    Ok(out.sz_inchi.unwrap_or_default())
}

fn fill_extra(rv: &mut ExtraInchiReturnValues, out: &InchiOutput) {
    if let Some(m) = &out.sz_message {
        rv.message = m.clone();
    }
    if let Some(l) = &out.sz_log {
        rv.log = l.clone();
    }
    if let Some(a) = &out.sz_aux_info {
        rv.aux_info = a.clone();
    }
}

/// `MolBlockToInchi`: the InChI of a mol block, parsed by the InChI library's
/// own MOL reader (`MakeINCHIFromMolfileText`).
pub fn mol_block_to_inchi(
    mol_block: &str,
    rv: &mut ExtraInchiReturnValues,
    options: Option<&str>,
) -> Result<String> {
    let opts = options.map(fix_option_symbol);
    let (code, out) =
        make_inchi_from_molfile_text(mol_block.as_bytes(), opts.as_deref().map(str::as_bytes));
    rv.return_code = code;
    fill_extra(rv, &out);
    Ok(out.sz_inchi.unwrap_or_default())
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
