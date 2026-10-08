//! SD file text (`RDKit::SDWriter::getText`, `_MolToSDStream`).
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SDWriter.cpp` (2026.09.1).
//! The output sinks (files, Python file objects) live in the bindings.
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::mol_file_writer::{MolWriterParams, mol_to_mol_block};
use crate::Result;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::dict::COMPUTED_PROP_NAME;

/// `_writePropToStream`
pub fn write_prop(out: &mut String, mol: &ROMol, name: &str, molid: i32) {
    let Some(v) = mol.get_prop(name) else {
        return;
    };
    let pval = v.to_rd_string();
    if name.contains('\n') {
        crate::rd_log_warning!(
            "WARNING: Skipping property {name} because the name includes a newline"
        );
        return;
    }
    if pval.contains("\r\n\r\n") || pval.contains("\n\n") {
        crate::rd_log_warning!(
            "WARNING: Skipping property {name} because the value includes an illegal blank line"
        );
        return;
    }
    out.push_str(&format!(">  <{name}>  "));
    if molid >= 0 {
        out.push_str(&format!("({}) ", molid + 1));
    }
    out.push('\n');
    out.push_str(&pval);
    out.push_str("\n\n");
}

/// `SDWriter::getText` / `_MolToSDStream`
pub fn mol_to_sd_text(
    mol: &ROMol,
    conf_id: i32,
    kekulize: bool,
    force_v3000: bool,
    molid: i32,
    props: Option<&[String]>,
) -> Result<String> {
    let params = MolWriterParams {
        include_stereo: true,
        kekulize,
        force_v3000,
        ..Default::default()
    };
    let mut out = mol_to_mol_block(mol, &params, conf_id)?;
    match props {
        Some(p) if !p.is_empty() => {
            for pi in p {
                if mol.has_prop(pi) {
                    write_prop(&mut out, mol, pi, molid);
                }
            }
        }
        _ => {
            let computed: Vec<String> = match mol.get_prop(COMPUTED_PROP_NAME) {
                Some(crate::rdgeneral::RDValue::VecString(v)) => v.clone(),
                _ => Vec::new(),
            };
            for pi in mol.get_prop_list(true, true) {
                if pi == COMPUTED_PROP_NAME
                    || pi == cp::_Name
                    || pi == "_MolFileInfo"
                    || pi == "_MolFileComments"
                    || pi == cp::_MolFileChiralFlag
                    || computed.contains(&pi)
                {
                    continue;
                }
                write_prop(&mut out, mol, &pi, molid);
            }
        }
    }
    out.push_str("$$$$\n");
    Ok(out)
}
