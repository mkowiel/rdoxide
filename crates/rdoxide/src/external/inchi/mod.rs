//! Rust port of the InChI Software library, version 1.07.3 — the version
//! RDKit 2026.09.1 downloads (`External/INCHI-API`, `FindInchi.cmake`:
//! <https://github.com/IUPAC-InChI/InChI/releases/download/v1.07.3/INCHI-1-SRC.zip>).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! The InChI library and programs are free software developed under the
//! auspices of the International Union of Pure and Applied Chemistry (IUPAC).
//! Originally developed at NIST. Modifications and additions by IUPAC and the
//! InChI Trust.
//!
//! Module layout and port plan: `docs/INCHI_PORT.md`. Each submodule names
//! the C sources (`INCHI_BASE/src/*.c`, `INCHI_API/libinchi/src/*.c`) it is
//! translated from. Ported so far: the InChIKey (`ikey_dll.c`,
//! `ikey_base26.c`) the sorting primitives (`ichisort.c`) and the error codes
//! and messages (`ichierr.c`) and the element helpers (`util.c`). The structure -> InChI and InChI -> structure pipelines
//! follow in later stages; the RDKit-facing API is
//! [`crate::external::inchi_api`].

pub mod el_data;
pub mod ichi;
pub mod ichi_bns;
pub mod ichi_io;
pub mod ichican2;
pub mod ichicano;
pub mod ichierr;
pub mod ichimak2;
pub mod ichimake;
pub mod ichimap;
pub mod ichimap4;
pub mod ichinorm;
pub mod ichiparm;
pub mod ichiprt1;
pub mod ichiprt2;
pub mod ichiprt3;
pub mod ichiread;
pub mod ichiring;
pub mod ichirvr1;
pub mod ichirvr2;
pub mod ichirvr3;
pub mod ichirvr4;
pub mod ichirvr5;
pub mod ichirvrs;
pub mod ichisort;
pub mod ichister;
pub mod ichitaut;
pub mod ikey_base26;
pub mod ikey_dll;
pub mod inchi_dll;
pub mod inchi_dll_b;
pub mod incomdef;
pub mod inpdef;
pub mod mol2atom;
pub mod molfmt;
pub mod molfmt3;
pub mod molfmt4;
pub mod readinch;
pub mod runichi;
pub mod runichi2;
pub mod runichi3;
pub mod runichi4;
pub mod sha256;
pub mod strutil;

pub mod util;

pub use ikey_dll::{
    InchiKey, InchiKeyCheck, InchiKeyError, check_inchikey, get_inchikey_from_inchi,
    get_std_inchikey_from_std_inchi,
};

/// `CURRENT_VER` (`bcf_s.h`): the InChI Software version that is ported.
pub const CURRENT_VER: &str = "1.07.3";

/// `INCHI_STRING_PREFIX` (`inchi_api.h`).
pub const INCHI_STRING_PREFIX: &[u8] = b"InChI=";
/// `LEN_INCHI_STRING_PREFIX` (`inchi_api.h`).
pub const LEN_INCHI_STRING_PREFIX: usize = 6;

/// `inchi_Ret_*` return codes of the InChI API (`inchi_api.h`); RDKit passes
/// them through as `ExtraInchiReturnValues::returnCode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum InchiRet {
    /// `inchi_Ret_BREAK`
    Break = -100,
    /// `inchi_Ret_SKIP`: not used in InChI library
    Skip = -2,
    /// `inchi_Ret_EOF`: no structural data has been provided
    Eof = -1,
    /// `inchi_Ret_OKAY`: success; no errors or warnings
    Okay = 0,
    /// `inchi_Ret_WARNING`: success; warning(s) issued
    Warning = 1,
    /// `inchi_Ret_ERROR`: error, no InChI has been created
    Error = 2,
    /// `inchi_Ret_FATAL`: severe error, no InChI has been created
    Fatal = 3,
    /// `inchi_Ret_UNKNOWN`: unknown program error
    Unknown = 4,
    /// `inchi_Ret_BUSY`: previous call to InChI has not returned yet
    Busy = 5,
}
