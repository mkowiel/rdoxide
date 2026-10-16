//! Daylight TDT text (`RDKit::TDTWriter`): a `$SMI<...>` record per molecule
//! with optional name, coordinates and property records.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/TDTWriter.cpp` (2026.09.1).
//! This type only produces text; the output sinks (files, Python file
//! objects) live in the bindings.
//!
//! Copyright (C) 2010 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::rdgeneral::dict::COMPUTED_PROP_NAME;
use crate::rdgeneral::strings::format_g;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};
use crate::{RDKitError, Result};

/// Formatting state of a `TDTWriter`.
#[derive(Debug, Clone)]
pub struct TdtWriter {
    props: Vec<String>,
    molid: u32,
    num_digits: usize,
    write_2d: bool,
    write_names: bool,
}

impl Default for TdtWriter {
    fn default() -> Self {
        TdtWriter {
            props: Vec::new(),
            molid: 0,
            num_digits: 4,
            write_2d: false,
            write_names: true,
        }
    }
}

impl TdtWriter {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of molecules written so far.
    pub fn num_mols(&self) -> u32 {
        self.molid
    }

    pub fn write_2d(&self) -> bool {
        self.write_2d
    }
    pub fn set_write_2d(&mut self, state: bool) {
        self.write_2d = state;
    }
    pub fn write_names(&self) -> bool {
        self.write_names
    }
    pub fn set_write_names(&mut self, state: bool) {
        self.write_names = state;
    }
    pub fn num_digits(&self) -> usize {
        self.num_digits
    }
    pub fn set_num_digits(&mut self, n: usize) {
        self.num_digits = n;
    }

    /// Sets the properties to write; with none all non-computed properties
    /// are written.
    pub fn set_props(&mut self, prop_names: Vec<String>) {
        if self.molid > 0 {
            crate::rd_log_warning!(
                "WARNING: Setting property list after a few molecules have been written\n"
            );
        }
        self.props = prop_names;
    }

    fn property(mol: &ROMol, name: &str) -> String {
        // we need to remove any line breaks in the output, replace them with
        // spaces
        let val = mol
            .get_prop(name)
            .map(RDValue::to_rd_string)
            .unwrap_or_default()
            .replace('\n', " ");
        format!("{name}<{val}>\n")
    }

    /// The text `write` appends for `mol`.
    pub fn format_mol(&mut self, mol: &mut ROMol, conf_id: i32) -> Result<String> {
        let mut out = String::new();
        // start by writing a "|" line unless this is the first line
        if self.molid > 0 {
            out.push_str("|\n");
        }
        out.push_str(&format!(
            "$SMI<{}>\n",
            mol_to_smiles(mol, &SmilesWriteParams::default())?
        ));
        if self.write_names
            && let Some(name) = mol.get_prop(cp::_Name)
        {
            out.push_str(&format!("NAME<{}>\n", name.to_rd_string()));
        }
        // do we need to write coordinates?
        if mol.num_conformers() > 0 {
            // the ordering of the atoms in the output SMILES:
            let ordering: Vec<u32> = match mol.get_prop(cp::_smilesAtomOutputOrder) {
                Some(RDValue::VecUnsignedInt(v)) => v.clone(),
                Some(RDValue::VecInt(v)) => v.iter().map(|&i| i as u32).collect(),
                _ => return Err(RDKitError::Key(cp::_smilesAtomOutputOrder.into())),
            };
            let conf = mol.conformer(conf_id)?;
            out.push_str(if self.write_2d { "2D<" } else { "3D<" });
            let d = self.num_digits;
            let n = ordering.len();
            for (i, &a) in ordering.iter().enumerate() {
                let p = conf.atom_pos(a)?;
                out.push_str(&format!("{},{}", format_g(p.x, d), format_g(p.y, d)));
                if !self.write_2d {
                    out.push_str(&format!(",{}", format_g(p.z, d)));
                }
                if i != n - 1 {
                    out.push(',');
                }
            }
            out.push_str(";>\n");
        }
        // now write the properties
        if !self.props.is_empty() {
            // the properties the user asked for
            for p in &self.props {
                if mol.has_prop(p) {
                    out.push_str(&Self::property(mol, p));
                }
            }
        } else {
            // otherwise all non computed properties
            let computed: Vec<String> = match mol.get_prop(COMPUTED_PROP_NAME) {
                Some(RDValue::VecString(c)) => c.clone(),
                _ => Vec::new(),
            };
            for p in mol.get_prop_list(true, true) {
                if p == COMPUTED_PROP_NAME
                    || p == cp::_Name
                    || p == "_MolFileInfo"
                    || p == "_MolFileComments"
                    || p == cp::_MolFileChiralFlag
                    || computed.contains(&p)
                {
                    continue;
                }
                out.push_str(&Self::property(mol, &p));
            }
        }
        self.molid += 1;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::file_parsers::tdt_mol_supplier::{TdtMolSupplier, TdtMolSupplierParams};
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    #[test]
    fn write_and_read_back() {
        let mut m = mol_from_smiles("CCO", &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        m.set_prop(cp::_Name, "ethanol".to_string(), false).unwrap();
        m.set_prop("foo", "bar\nbaz".to_string(), false).unwrap();
        let mut w = TdtWriter::new();
        let mut text = w.format_mol(&mut m, -1).unwrap();
        text.push_str(&w.format_mol(&mut m, -1).unwrap());
        assert_eq!(
            text,
            "$SMI<CCO>\nNAME<ethanol>\nfoo<bar baz>\n|\n$SMI<CCO>\nNAME<ethanol>\nfoo<bar baz>\n"
        );
        let params = TdtMolSupplierParams {
            name_record: "NAME".into(),
            ..Default::default()
        };
        let mut s = TdtMolSupplier::from_bytes(text.into_bytes(), params);
        assert_eq!(s.length().unwrap(), 2);
        let m2 = s.get(1).unwrap().unwrap();
        assert_eq!(m2.name(), "ethanol");
        assert_eq!(
            m2.get_prop("foo").unwrap().to_rd_string(),
            "bar baz".to_string()
        );
    }
}
