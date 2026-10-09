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
    /// `MolFileUnhandledFeatureException`: a CTAB feature that is not
    /// supported; the mol file readers skip the molecule instead of failing.
    UnhandledFeature(String),
    /// `ConformerException` -> `ValueError("Bad Conformer Id")`
    Conformer(String),
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
    /// `URANGE_CHECK(x, hi)` failure raised at `file:line` of the C++ sources
    /// (which the error text names); logged to the error log like upstream.
    pub fn range_at(
        name: &str,
        x: impl fmt::Display,
        hi: impl fmt::Display,
        file: &str,
        line: u32,
    ) -> Self {
        crate::rd_log_error!(
            "\n\n****\nRange Error\n{name}\nViolation occurred on line {line} in file {file}\nFailed Expression: {x} < {hi}\n****\n"
        );
        RDKitError::Invariant {
            kind: "Range Error",
            message: format!(
                "{name}\n\tViolation occurred on line {line} in file {file}\n\tFailed Expression: {x} < {hi}"
            ),
        }
    }

    /// A failed `CHECK_INVARIANT(expr, message)` of the C++ sources at
    /// `file:line` (the last line of the statement); logged like upstream.
    pub fn invariant_at(message: impl fmt::Display, expr: &str, file: &str, line: u32) -> Self {
        Self::check_at("Invariant Violation", message, expr, file, line)
    }

    /// A failed `PRECONDITION(expr, message)` of the C++ sources at
    /// `file:line` (the last line of the statement); logged like upstream.
    pub fn precondition_at(message: impl fmt::Display, expr: &str, file: &str, line: u32) -> Self {
        Self::check_at("Pre-condition Violation", message, expr, file, line)
    }

    /// A failed `TEST_ASSERT(expr)` of the C++ sources at `file:line` (the
    /// last line of the statement); logged like upstream.
    pub fn test_assert_at(expr: &str, file: &str, line: u32) -> Self {
        Self::check_at("Test Assert", "Expression Failed: ", expr, file, line)
    }

    fn check_at(
        kind: &'static str,
        message: impl fmt::Display,
        expr: &str,
        file: &str,
        line: u32,
    ) -> Self {
        crate::rd_log_error!(
            "\n\n****\n{kind}\n{message}\nViolation occurred on line {line} in file {file}\nFailed Expression: {expr}\n****\n"
        );
        RDKitError::Invariant {
            kind,
            message: format!(
                "{message}\n\tViolation occurred on line {line} in file {file}\n\tFailed Expression: {expr}"
            ),
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
                    "{kind}\n\t{message}\n\tRDKIT: {}\n\tBOOST: \n",
                    crate::rdgeneral::RDKIT_VERSION
                )
            }
            RDKitError::Value(m)
            | RDKitError::Key(m)
            | RDKitError::Index(m)
            | RDKitError::Runtime(m)
            | RDKitError::FileParse(m)
            | RDKitError::UnhandledFeature(m)
            | RDKitError::Conformer(m)
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
