//! Reading molecules from the `rdkit:mol` metadata of RDKit SVG drawings.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/SVGParser.cpp` (2026.09.1).
//!
//! Copyright (C) 2016 Greg Landrum
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::marvin::xml::{XmlTree, read_xml};
use crate::geometry::Point3D;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::RWMol;
use crate::graphmol::smiles_parse::{atom_from_smiles, bond_from_smiles};
use crate::{RDKitError, Result, rd_log_warning};

/// `ptree::get<T>(path)` on an attribute: a missing attribute or a value that
/// does not convert is an error.
fn required_attr(node: &XmlTree, name: &str) -> Result<String> {
    node.child("<xmlattr>")
        .and_then(|a| a.child(name))
        .map(|n| n.data.clone())
        .ok_or_else(|| RDKitError::Runtime(format!("No such node (<xmlattr>.{name})")))
}

fn attr_f64(node: &XmlTree, name: &str) -> Result<f64> {
    let s = node.attr(name, "0");
    s.trim().parse().map_err(|_| {
        RDKitError::Runtime(format!(
            "conversion of data to type \"double\" failed ({name})"
        ))
    })
}

fn attr_idx(node: &XmlTree, name: &str) -> Result<u32> {
    let s = required_attr(node, name)?;
    let v: i32 = s.trim().parse().map_err(|_| {
        RDKitError::Runtime(format!(
            "conversion of data to type \"int\" failed ({name})"
        ))
    })?;
    Ok((v as u32).wrapping_sub(1))
}

/// `ptreeToMol`
fn ptree_to_mol(mol: &mut RWMol, mol_e: &XmlTree) -> Result<()> {
    let mut pts = Vec::new();
    let mut is_3d = false;
    for at_e in mol_e.children_named("rdkit:atom") {
        let smi = required_attr(at_e, "atom-smiles")?;
        let Some(atom) = atom_from_smiles(&smi) else {
            rd_log_warning!(" Could not convert SMILES '{smi}' to a atom. Ignoring it.");
            continue;
        };
        mol.add_atom(atom, false);
        pts.push(Point3D::new(
            attr_f64(at_e, "x")?,
            attr_f64(at_e, "y")?,
            attr_f64(at_e, "z")?,
        ));
        if at_e.attr("z", "0") != "0" {
            is_3d = true;
        }
    }
    for b_e in mol_e.children_named("rdkit:bond") {
        let smi = required_attr(b_e, "bond-smiles")?;
        let Some(mut bond) = bond_from_smiles(&smi) else {
            rd_log_warning!(" Could not convert SMILES '{smi}' to a bond. Ignoring it.");
            continue;
        };
        bond.set_begin_atom_idx_unchecked(attr_idx(b_e, "begin-atom-idx")?);
        bond.set_end_atom_idx_unchecked(attr_idx(b_e, "end-atom-idx")?);
        mol.add_bond_object(bond)?;
    }
    let mut conf = Conformer::new(mol.num_atoms() as usize);
    for (i, p) in pts.into_iter().enumerate() {
        conf.set_atom_pos(i as u32, p)?;
    }
    conf.set_3d(is_3d);
    mol.add_conformer(conf, true)?;
    Ok(())
}

/// `RDKitSVGToMol`: the last `rdkit:mol` in the SVG metadata, `None` if
/// there is none. Malformed XML is an error.
pub fn rdkit_svg_to_mol(svg: &str, sanitize: bool, remove_hs_flag: bool) -> Result<Option<RWMol>> {
    let tree = read_xml(svg)?;
    let mut res = None;
    if let Some(metadata) = tree.get_child("svg.metadata") {
        for mol_e in metadata.children_named("rdkit:mol") {
            let mut mol = RWMol::new();
            ptree_to_mol(&mut mol, mol_e)?;
            if mol.num_atoms() > 0 {
                if remove_hs_flag {
                    remove_hs(&mut mol, &RemoveHsParameters::default(), sanitize)?;
                } else if sanitize {
                    sanitize_mol(&mut mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
                }
            }
            res = Some(mol);
        }
    }
    Ok(res)
}
