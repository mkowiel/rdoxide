//! Sequence and FASTA writers.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SequenceWriters.cpp`
//! (2026.09.1).
//!
//! Copyright (C) 2015,2016 Greg Landrum and NextMove Software.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::sequence_parser::AMINO_ACID_CODES;
use crate::graphmol::monomer_info::{AtomMonomerInfo, AtomMonomerType};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// `getOneLetterAACode`: upper case for L amino acids, lower case for the D
/// forms, `X` for anything else.
fn one_letter_aa_code(residue: &str) -> char {
    AMINO_ACID_CODES
        .iter()
        .find_map(|&(c, l, d)| {
            if residue == l {
                Some(c as char)
            } else if residue == d {
                Some(c.to_ascii_lowercase() as char)
            } else {
                None
            }
        })
        .unwrap_or('X')
}

/// `getOneLetterNACode`: the base of `"  A"` (RNA) and `" DA"` (DNA).
fn one_letter_na_code(residue: &str) -> char {
    match residue.as_bytes() {
        [b' ', b' ' | b'D', c @ (b'A' | b'C' | b'G' | b'T' | b'U')] => *c as char,
        _ => 'X',
    }
}

/// `MolToSequence`: one letter codes of the residues, found through their
/// alpha carbons (`" CA "`) and C1' atoms.
///
/// Upstream tracks the chain to emit `.` between chains but only compares
/// against a chain that is never assigned for the first residue, so no
/// separator is ever written; this is reproduced.
pub fn mol_to_sequence(mol: &ROMol) -> String {
    mol.atoms()
        .iter()
        .filter_map(|atom| {
            let info: &AtomMonomerInfo = atom.monomer_info()?;
            if info.monomer_type != AtomMonomerType::PdbResidue {
                return None;
            }
            match info.name.as_str() {
                " CA " => Some(one_letter_aa_code(&info.residue_name)),
                " C1'" => Some(one_letter_na_code(&info.residue_name)),
                _ => None,
            }
        })
        .collect()
}

/// `MolToFASTA`: empty when the molecule has no residues.
pub fn mol_to_fasta(mol: &ROMol) -> String {
    let seq = mol_to_sequence(mol);
    if seq.is_empty() {
        return String::new();
    }
    let name = mol
        .get_prop(cp::_Name)
        .map(|v| v.to_rd_string())
        .unwrap_or_default();
    format!(">{name}\n{seq}\n")
}
