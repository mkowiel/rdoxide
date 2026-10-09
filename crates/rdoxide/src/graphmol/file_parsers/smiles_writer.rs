//! SMILES text files (`RDKit::SmilesWriter`): one molecule per line with an
//! optional header, name column and property columns.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SmilesWriter.cpp` (2026.09.1).
//! This type only produces text; the output sinks (files, Python file
//! objects) live in the bindings.
//!
//! Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::graphmol::smiles_parse::write::mol_to_smiles;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// Formatting state of a `SmilesWriter`.
#[derive(Debug, Clone)]
pub struct SmilesWriter {
    delimiter: String,
    name_header: String,
    include_header: bool,
    isomeric_smiles: bool,
    kekule_smiles: bool,
    props: Vec<String>,
    molid: u32,
}

impl SmilesWriter {
    /// An empty `name_header` leaves the name column out.
    pub fn new(
        delimiter: &str,
        name_header: &str,
        include_header: bool,
        isomeric_smiles: bool,
        kekule_smiles: bool,
    ) -> Self {
        SmilesWriter {
            delimiter: delimiter.to_string(),
            name_header: name_header.to_string(),
            include_header,
            isomeric_smiles,
            kekule_smiles,
            props: Vec::new(),
            molid: 0,
        }
    }

    /// Number of molecules written so far.
    pub fn num_mols(&self) -> u32 {
        self.molid
    }

    /// Sets the property columns; ignored (with an error log) once a molecule
    /// has been written.
    pub fn set_props(&mut self, prop_names: Vec<String>) {
        if self.molid > 0 {
            crate::rd_log_error!("ERROR: Atleast one molecule has already been written");
            crate::rd_log_error!("ERROR: Cannot set properties now - ignoring setProps");
            return;
        }
        self.props = prop_names;
    }

    fn header(&self) -> String {
        let mut cols = vec!["SMILES".to_string()];
        if !self.name_header.is_empty() {
            cols.push(self.name_header.clone());
        }
        cols.extend(self.props.iter().cloned());
        let mut out = cols.join(&self.delimiter);
        // `SMILES` and the name are always followed by a delimiter, even with
        // no property columns after them
        if self.props.is_empty() {
            out.push_str(&self.delimiter);
        }
        out.push('\n');
        out
    }

    /// The text `write` appends for `mol`: the header before the first
    /// molecule (when enabled), then the molecule's line.
    pub fn format_mol(&mut self, mol: &mut ROMol) -> Result<String> {
        let mut out = String::new();
        if self.molid == 0 && self.include_header {
            out.push_str(&self.header());
        }
        let params = SmilesWriteParams {
            do_isomeric_smiles: self.isomeric_smiles,
            do_kekule: self.kekule_smiles,
            ..SmilesWriteParams::default()
        };
        out.push_str(&mol_to_smiles(mol, &params)?);
        if !self.name_header.is_empty() {
            let name = mol
                .get_prop(cp::_Name)
                .map(|v| v.to_rd_string())
                .filter(|n| !n.is_empty())
                .unwrap_or_else(|| self.molid.to_string());
            out.push_str(&self.delimiter);
            out.push_str(&name);
        }
        for p in &self.props {
            out.push_str(&self.delimiter);
            if let Some(v) = mol.get_prop(p) {
                out.push_str(&v.to_rd_string());
            }
        }
        out.push('\n');
        self.molid += 1;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    fn mol(smi: &str) -> ROMol {
        mol_from_smiles(smi, &SmilesParserParams::default())
            .unwrap()
            .unwrap()
    }

    #[test]
    fn header_name_and_properties() {
        let mut w = SmilesWriter::new(" ", "Name", true, true, false);
        w.set_props(vec!["a".into(), "b".into()]);
        let mut m = mol("OCC");
        m.set_prop("_Name", "ethanol".to_string(), false).unwrap();
        m.set_prop("a", 3i32, false).unwrap();
        assert_eq!(
            w.format_mol(&mut m).unwrap(),
            "SMILES Name a b\nCCO ethanol 3 \n"
        );
        // the default name is the running index
        assert_eq!(w.format_mol(&mut mol("C")).unwrap(), "C 1  \n");
        assert_eq!(w.num_mols(), 2);
    }

    #[test]
    fn no_name_and_no_header() {
        let mut w = SmilesWriter::new(",", "", false, true, false);
        assert_eq!(w.format_mol(&mut mol("c1ccccc1")).unwrap(), "c1ccccc1\n");
        let mut w = SmilesWriter::new("\t", "N", true, true, false);
        assert_eq!(w.format_mol(&mut mol("C")).unwrap(), "SMILES\tN\t\nC\t0\n");
    }

    #[test]
    fn props_are_frozen_after_the_first_molecule() {
        let mut w = SmilesWriter::new(" ", "", false, true, false);
        w.format_mol(&mut mol("C")).unwrap();
        w.set_props(vec!["x".into()]);
        assert_eq!(w.format_mol(&mut mol("C")).unwrap(), "C\n");
    }
}
