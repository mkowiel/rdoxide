//! Sequence and FASTA reading and writing (upstream `Wrap/rdmolfiles.cpp`):
//! `MolFromSequence`, `MolToSequence`, `MolFromFASTA`, `MolToFASTA`,
//! `MolFromHELM` and `MolToHELM`.

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::helm_parser::helm_to_mol;
use rdoxide::graphmol::file_parsers::helm_writer::mol_to_helm;
use rdoxide::graphmol::file_parsers::sequence_parser::{fasta_to_mol, sequence_to_mol};
use rdoxide::graphmol::file_parsers::sequence_writer::{mol_to_fasta, mol_to_sequence};

use crate::mol::PyMol;

/// The upstream wrappers log `FileParseException`s as warnings and swallow
/// every other error, returning `None`.
fn finish(res: rdoxide::Result<Option<rdoxide::graphmol::RWMol>>) -> Option<PyMol> {
    match res {
        Ok(mol) => mol.map(|mol| PyMol { mol }),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            None
        }
        Err(_) => None,
    }
}

/// `pyObjectToString`: accepts str or bytes.
pub(crate) fn text_arg(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    match obj.extract::<String>() {
        Ok(s) => Ok(s),
        Err(_) => Ok(String::from_utf8_lossy(&obj.extract::<Vec<u8>>()?).into_owned()),
    }
}

/// Construct a molecule from a sequence string (currently supports standard amino acids, DNA and RNA bases).
///
///   ARGUMENTS:
///
///     - text: string containing the sequence
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - flavor: (optional)
///         - 0 Protein, L amino acids (default)
///         - 1 Protein, D amino acids
///         - 2 RNA, no cap
///         - 3 RNA, 5' cap
///         - 4 RNA, 3' cap
///         - 5 RNA, both caps
///         - 6 DNA, no cap
///         - 7 DNA, 5' cap
///         - 8 DNA, 3' cap
///         - 9 DNA, both caps
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromSequence", signature = (text, sanitize = true, flavor = 0))]
fn py_mol_from_sequence(
    py: Python<'_>,
    text: &Bound<'_, PyAny>,
    sanitize: bool,
    flavor: i32,
) -> PyResult<Option<PyMol>> {
    let text = text_arg(text)?;
    Ok(finish(
        py.detach(|| sequence_to_mol(&text, sanitize, flavor)),
    ))
}

/// Returns the sequence string for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///
///   NOTE: the molecule should contain monomer information in AtomMonomerInfo structures
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToSequence", signature = (mol))]
fn py_mol_to_sequence(mol: PyRef<'_, PyMol>) -> String {
    mol_to_sequence(&mol.mol)
}

/// Construct a molecule from a FASTA string (currently supports standard amino acids, DNA and RNA bases).
///
///   ARGUMENTS:
///
///     - text: string containing the FASTA
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
/// - flavor: (optional)
///     - 0 Protein, L amino acids (default)
///     - 1 Protein, D amino acids
///     - 2 RNA, no cap
///     - 3 RNA, 5' cap
///     - 4 RNA, 3' cap
///     - 5 RNA, both caps
///     - 6 DNA, no cap
///     - 7 DNA, 5' cap
///     - 8 DNA, 3' cap
///     - 9 DNA, both caps
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromFASTA", signature = (text, sanitize = true, flavor = 0))]
fn py_mol_from_fasta(
    py: Python<'_>,
    text: &Bound<'_, PyAny>,
    sanitize: bool,
    flavor: i32,
) -> PyResult<Option<PyMol>> {
    let text = text_arg(text)?;
    Ok(finish(py.detach(|| fasta_to_mol(&text, sanitize, flavor))))
}

/// Returns the FASTA string for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///
///   NOTE: the molecule should contain monomer information in AtomMonomerInfo structures
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToFASTA", signature = (mol))]
fn py_mol_to_fasta(mol: PyRef<'_, PyMol>) -> String {
    mol_to_fasta(&mol.mol)
}

/// Construct a molecule from a HELM string (currently supports standard amino acids, DNA and RNA bases).
///
///   ARGUMENTS:
///
///     - text: string containing the HELM
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromHELM", signature = (text, sanitize = true))]
fn py_mol_from_helm(
    py: Python<'_>,
    text: &Bound<'_, PyAny>,
    sanitize: bool,
) -> PyResult<Option<PyMol>> {
    let text = text_arg(text)?;
    Ok(finish(py.detach(|| helm_to_mol(&text, sanitize))))
}

/// Returns the HELM string for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///
///   NOTE: the molecule should contain monomer information in AtomMonomerInfo structures
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToHELM", signature = (mol))]
fn py_mol_to_helm(mol: PyRef<'_, PyMol>) -> String {
    mol_to_helm(&mol.mol)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_from_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_from_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_from_helm, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_helm, m)?)?;
    Ok(())
}
