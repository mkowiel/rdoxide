//! Error types. Each variant corresponds to an upstream RDKit exception class
//! so that the Python layer can raise the same Python exception types.

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RDKitError {
    /// `Invar::Invariant` (PRECONDITION / POSTCONDITION / CHECK_INVARIANT
    /// violations); surfaces in Python as `RuntimeError`.
    Invariant { kind: &'static str, message: String },
    /// `ValueErrorException` -> `ValueError`
    Value(String),
    /// `KeyErrorException` -> `KeyError`
    Key(String),
    /// `IndexErrorException` -> `IndexError`
    Index(String),
    /// `std::runtime_error` and friends (e.g. `MolPicklerException`) ->
    /// `RuntimeError` with the plain message.
    Runtime(String),
    /// `FileParseException` -> `RuntimeError` (most Python readers catch it
    /// and return `None`).
    FileParse(String),
    /// `BadFileException` -> `OSError`
    BadFile(String),
    /// Sanitization failures (`MolSanitizeException` hierarchy).
    Sanitize(SanitizeError),
}

/// `MolSanitizeException` and subclasses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SanitizeError {
    Generic {
        message: String,
    },
    AtomValence {
        atom_idx: u32,
        message: String,
    },
    AtomKekulize {
        atom_idx: u32,
        message: String,
    },
    Kekulize {
        atom_indices: Vec<u32>,
        message: String,
    },
}

impl RDKitError {
    pub fn precondition(message: impl Into<String>) -> Self {
        RDKitError::Invariant {
            kind: "Pre-condition Violation",
            message: message.into(),
        }
    }
    pub fn postcondition(message: impl Into<String>) -> Self {
        RDKitError::Invariant {
            kind: "Post-condition Violation",
            message: message.into(),
        }
    }
    pub fn invariant(message: impl Into<String>) -> Self {
        RDKitError::Invariant {
            kind: "Invariant Violation",
            message: message.into(),
        }
    }
    /// `URANGE_CHECK(x, hi)` failure; logged to the error log like upstream.
    pub fn range(name: &str, x: impl fmt::Display, hi: impl fmt::Display) -> Self {
        let err = RDKitError::Invariant {
            kind: "Range Error",
            message: format!("{name}\n\tFailed Expression: {x} < {hi}"),
        };
        crate::rd_log_error!("\n\n****\n{err}****\n");
        err
    }
}

impl fmt::Display for RDKitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RDKitError::Invariant { kind, message } => {
                write!(
                    f,
                    "{kind}\n\t{message}\n\tRDKIT: {}\n",
                    crate::rdgeneral::RDKIT_VERSION
                )
            }
            RDKitError::Value(m)
            | RDKitError::Key(m)
            | RDKitError::Index(m)
            | RDKitError::Runtime(m)
            | RDKitError::FileParse(m)
            | RDKitError::BadFile(m) => f.write_str(m),
            RDKitError::Sanitize(e) => match e {
                SanitizeError::Generic { message }
                | SanitizeError::AtomValence { message, .. }
                | SanitizeError::AtomKekulize { message, .. }
                | SanitizeError::Kekulize { message, .. } => f.write_str(message),
            },
        }
    }
}

impl std::error::Error for RDKitError {}

pub type Result<T, E = RDKitError> = std::result::Result<T, E>;
