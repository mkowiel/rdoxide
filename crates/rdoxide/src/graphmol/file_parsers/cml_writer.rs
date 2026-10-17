//! CML (Chemical Markup Language) writer.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/CMLWriter.cpp` (2026.09.1).
//! The output reproduces what `boost::property_tree::write_xml` writes with a
//! two-space indent.
//!
//! Copyright (C) 2020 Greg Landrum and other RDKit contributors
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondType};
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::romol::ROMol;
use crate::graphmol::wedge_bonds::{BondWedgingParameters, wedge_mol_bonds};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result, rd_log_info};

/// `boost::property_tree` XML character entity encoding.
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            c => out.push(c),
        }
    }
    out
}

/// A start tag with its attributes in insertion order.
fn open_tag(name: &str, attrs: &[(&str, String)]) -> String {
    let mut s = format!("<{name}");
    for (k, v) in attrs {
        s.push_str(&format!(" {k}=\"{}\"", escape(v)));
    }
    s
}

fn bond_order(t: BondType) -> String {
    match t {
        BondType::Single
        | BondType::DativeOne
        | BondType::Dative
        | BondType::DativeL
        | BondType::DativeR => "S".into(),
        BondType::Double => "D".into(),
        BondType::Triple => "T".into(),
        BondType::Aromatic => "A".into(),
        // RDKit extensions: higher and half-integer orders
        BondType::Quadruple => "4".into(),
        BondType::Quintuple => "5".into(),
        BondType::Hextuple => "6".into(),
        BondType::ThreeCenter => "0.5".into(),
        BondType::OneAndAHalf => "1.5".into(),
        BondType::TwoAndAHalf => "2.5".into(),
        BondType::ThreeAndAHalf => "3.5".into(),
        BondType::FourAndAHalf => "4.5".into(),
        BondType::FiveAndAHalf => "5.5".into(),
        other => {
            rd_log_info!("CMLWriter: Unsupported BondType {}\n", other as i32);
            String::new()
        }
    }
}

/// `MolToCMLBlock(mol, confId, kekulize)`
pub fn mol_to_cml_block(mol: &ROMol, conf_id: i32, kekulize_mol: bool) -> Result<String> {
    let mut rw = mol.clone();
    if kekulize_mol {
        if rw.needs_update_property_cache() {
            rw.update_property_cache(false)?;
        }
        kekulize(&mut rw, true, true, 100)?;
    }
    let conf = if rw.num_conformers() > 0 {
        let conf = rw.conformer(conf_id)?.clone();
        // wedge bonds so that we can use that info
        wedge_mol_bonds(&mut rw, Some(&conf), &BondWedgingParameters::default())?;
        Some(conf)
    } else {
        None
    };

    let mut mol_charge = 0;
    let mut mol_radicals = 0u32;
    let mut atoms = String::new();
    for i in 0..rw.num_atoms() {
        let a = rw.atom(i);
        let element = if a.atomic_num() != 0 {
            a.symbol()
        } else {
            "Du".to_string()
        };
        let charge = a.formal_charge();
        mol_charge += charge;
        let mut attrs = vec![
            ("id", format!("a{i}")),
            ("elementType", element),
            ("formalCharge", charge.to_string()),
            ("hydrogenCount", rw.total_num_hs(i, true)?.to_string()),
        ];
        if a.isotope() != 0 {
            attrs.push(("isotopeNumber", a.isotope().to_string()));
        }
        mol_radicals += a.num_radical_electrons();
        if let Some(conf) = &conf {
            let pos = conf.atom_pos(i)?;
            if conf.is_3d() {
                attrs.push(("x3", format!("{:.6}", pos.x)));
                attrs.push(("y3", format!("{:.6}", pos.y)));
                attrs.push(("z3", format!("{:.6}", pos.z)));
            } else {
                attrs.push(("x2", format!("{:.6}", pos.x)));
                attrs.push(("y2", format!("{:.6}", pos.y)));
            }
        }
        let parity = match a.chiral_tag() {
            ChiralType::TetrahedralCCW => 1,
            ChiralType::TetrahedralCW => -1,
            _ => 0,
        };
        let tag = open_tag("atom", &attrs);
        if parity == 0 {
            atoms.push_str(&format!("      {tag}/>\n"));
        } else {
            let mut nbrs: Vec<u32> = rw.atom_neighbors(i).collect();
            // CML wants exactly 4 atomRefs: an implicit H is the atom itself
            while nbrs.len() < 4 {
                nbrs.push(i);
            }
            let refs = format!("a{} a{} a{} a{}", nbrs[0], nbrs[1], nbrs[2], nbrs[3]);
            let parity_tag = open_tag("atomParity", &[("atomRefs4", refs)]);
            atoms.push_str(&format!(
                "      {tag}>\n        {parity_tag}>{parity}</atomParity>\n      </atom>\n"
            ));
        }
    }

    let mut mol_attrs = vec![
        ("id", format!("m{conf_id}")),
        ("formalCharge", mol_charge.to_string()),
    ];
    if mol_radicals < 2 {
        mol_attrs.push(("spinMultiplicity", (mol_radicals + 1).to_string()));
    } else {
        rd_log_info!(
            "CMLWriter: Unable to determine molecule/@spinMultiplicity ({mol_radicals} radical electrons)\n"
        );
    }

    let mut bonds = String::new();
    let mut bond_id = 0;
    for src in 0..rw.num_atoms() {
        for b in rw.atom_bonds(src) {
            let dst = b.other_atom_idx(src)?;
            if dst < src {
                continue;
            }
            let mut attrs = vec![
                ("atomRefs2", format!("a{src} a{dst}")),
                ("id", format!("b{bond_id}")),
                ("order", bond_order(b.bond_type())),
            ];
            bond_id += 1;
            match b.bond_dir() {
                BondDir::BeginDash => attrs.push(("bondStereo", "H".into())),
                BondDir::BeginWedge => attrs.push(("bondStereo", "W".into())),
                _ => {}
            }
            bonds.push_str(&format!("      {}/>\n", open_tag("bond", &attrs)));
        }
    }

    let mut out = String::from("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    out.push_str(&open_tag(
        "cml",
        &[
            ("xmlns", "http://www.xml-cml.org/schema".into()),
            (
                "xmlns:convention",
                "http://www.xml-cml.org/convention/".into(),
            ),
            ("convention", "convention:molecular".into()),
        ],
    ));
    out.push_str(">\n");
    out.push_str(&format!("  {}>\n", open_tag("molecule", &mol_attrs)));
    if let Some(name) = rw.get_prop(cp::_Name) {
        let name = name.to_rd_string();
        if !name.is_empty() {
            out.push_str(&format!("    <name>{}</name>\n", escape(&name)));
        }
    }
    for (tag, body) in [("atomArray", &atoms), ("bondArray", &bonds)] {
        if body.is_empty() {
            out.push_str(&format!("    <{tag}/>\n"));
        } else {
            out.push_str(&format!("    <{tag}>\n{body}    </{tag}>\n"));
        }
    }
    out.push_str("  </molecule>\n</cml>\n");
    Ok(out)
}

/// `MolToCMLFile`
pub fn mol_to_cml_file(
    mol: &ROMol,
    file_name: &str,
    conf_id: i32,
    kekulize_mol: bool,
) -> Result<()> {
    let block = mol_to_cml_block(mol, conf_id, kekulize_mol)?;
    std::fs::write(file_name, block)
        .map_err(|_| RDKitError::BadFile(format!("Bad output file {file_name}")))
}
