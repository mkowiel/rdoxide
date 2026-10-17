//! BioCAD TPL file parser and writer.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/{TplFileParser,TplFileWriter}.cpp`
//! (2026.09.1).
//!
//! Copyright (C) 2007-2024 Greg Landrum
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::utils::{LineStream, lexical_cast_int, lexical_cast_unsigned, strip_spaces};
use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::BondType;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_ops::{sanitize_flags::SANITIZE_ALL, sanitize_mol};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::strings::format_g;
use crate::{RDKitError, Result};

/// `TPLParserParams`
#[derive(Debug, Clone, Copy)]
pub struct TplParserParams {
    pub sanitize: bool,
    /// skip the first conformer (needed for TPLs written by the CombiCode)
    pub skip_first_conf: bool,
}

impl Default for TplParserParams {
    fn default() -> Self {
        TplParserParams {
            sanitize: true,
            skip_first_conf: false,
        }
    }
}

fn file_parse_error(message: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(message.into())
}

/// `boost::split(..., is_any_of(" \t"), token_compress_on)`: adjacent
/// separators count as one, but a leading or trailing separator still gives
/// an empty first or last token.
fn split_compress(s: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut in_delim = false;
    for (i, c) in s.char_indices() {
        let is_delim = c == ' ' || c == '\t';
        if is_delim && !in_delim {
            tokens.push(&s[start..i]);
            in_delim = true;
        } else if !is_delim && in_delim {
            start = i;
            in_delim = false;
        }
    }
    tokens.push(if in_delim { "" } else { &s[start..] });
    tokens
}

fn cast_double(s: &str) -> Result<f64> {
    String::from_utf8_lossy(strip_spaces(s.as_bytes()))
        .parse::<f64>()
        .map_err(|_| RDKitError::Value("bad lexical cast".into()))
}

fn cast_unsigned(s: &str) -> Result<u32> {
    lexical_cast_unsigned(strip_spaces(s.as_bytes()))
}

fn cast_int(s: &str) -> Result<i32> {
    lexical_cast_int(strip_spaces(s.as_bytes()))
}

fn line_string(stream: &mut LineStream<'_>) -> String {
    String::from_utf8_lossy(&stream.get_line()).into_owned()
}

/// `ParseTPLAtomLine`
fn parse_atom_line(text: &str, line_num: u32, mol: &mut RWMol, conf: &mut Conformer) -> Result<()> {
    let split = split_compress(text);
    if split.len() < 8 {
        return Err(file_parse_error(format!(
            "Atom line {line_num} only has {} tokens. 8 are required.\n",
            split.len()
        )));
    }
    let atom = Atom::from_symbol(split[1])?;
    let atom_id = mol.add_atom(atom, false);
    let charge = cast_int(split[2])?;
    let partial = cast_double(split[3])?;
    // coords in TPL files are in picometers
    let (x, y, z) = (
        cast_double(split[4])? / 100.0,
        cast_double(split[5])? / 100.0,
        cast_double(split[6])? / 100.0,
    );
    let n_bonds = cast_unsigned(split[7])?;
    let atom = mol.atom_mut(atom_id);
    atom.set_formal_charge(charge);
    atom.set_prop("TPLCharge", partial, false)?;
    conf.set_atom_pos(atom_id, Point3D::new(x, y, z))?;
    // the only other information of interest is the stereochemistry, and
    // only for four bonds
    if n_bonds == 4 && split.len() > 8 + n_bonds as usize {
        atom.set_prop(
            "TPLStereoFlag",
            split[8 + n_bonds as usize].to_string(),
            false,
        )?;
    }
    Ok(())
}

/// `ParseTPLBondLine`
fn parse_bond_line(text: &str, line_num: u32, mol: &mut RWMol) -> Result<()> {
    let split = split_compress(text);
    if split.len() < 5 {
        return Err(file_parse_error(format!(
            "Bond line {line_num} only has {} tokens. 5 are required.\n",
            split.len()
        )));
    }
    let order: String = split[1].trim().chars().take(3).collect();
    let bond_type = match order.as_str() {
        "1.5" => BondType::Aromatic,
        "1.0" => BondType::Single,
        "2.0" => BondType::Double,
        "3.0" => BondType::Triple,
        _ => {
            return Err(file_parse_error(format!(
                "Bond line {line_num} has unknown order: {order}\n"
            )));
        }
    };
    let idx1 = cast_unsigned(split[2])?.wrapping_sub(1);
    let idx2 = cast_unsigned(split[3])?.wrapping_sub(1);
    let bond_idx = mol.add_bond(idx1, idx2, bond_type)? - 1;
    let flag2 = split.get(5).copied().unwrap_or("");
    let bond = mol.bond_mut(bond_idx);
    bond.props.set_val("TPLBondDir1", split[4].to_string());
    bond.props.set_val("TPLBondDir2", flag2.to_string());
    Ok(())
}

/// `ParseConfData`
fn parse_conf_data(
    stream: &mut LineStream<'_>,
    line: &mut u32,
    mol: &mut RWMol,
    conf_id: u32,
) -> Result<Conformer> {
    *line += 1;
    let text = line_string(stream);
    if split_compress(&text)[0] != "NAME" {
        return Err(file_parse_error(format!(
            "Did not find NAME tag on line {line} while reading conformer  {conf_id}\n"
        )));
    }
    let name = text.get(4..).unwrap_or("").trim().to_string();
    let prop_name = format!("Conf_{}{}", mol.num_conformers(), cp::_Name);
    mol.set_prop(&prop_name, name, false)?;

    let mut conf = Conformer::new(mol.num_atoms() as usize);
    for i in 0..mol.num_atoms() {
        *line += 1;
        let text = line_string(stream);
        if stream.eof() {
            return Err(file_parse_error(format!(
                "EOF hit while reading conformer  {conf_id}\n"
            )));
        }
        let trimmed = text.trim_matches([' ', '\t', '\r', '\n']);
        let split = split_compress(trimmed);
        if split.len() < 3 {
            return Err(file_parse_error(format!(
                "Did not find enough fields on line {line} while reading conformer  {conf_id}\n"
            )));
        }
        conf.set_atom_pos(
            i,
            Point3D::new(
                cast_double(split[0])? / 100.0,
                cast_double(split[1])? / 100.0,
                cast_double(split[2])? / 100.0,
            ),
        )?;
    }
    Ok(conf)
}

/// `MolFromTPLDataStream`
fn mol_from_tpl_stream(
    stream: &mut LineStream<'_>,
    line: &mut u32,
    params: &TplParserParams,
) -> Result<Option<RWMol>> {
    // format line and comment line
    for _ in 0..2 {
        *line += 1;
        line_string(stream);
        if stream.eof() {
            return Ok(None);
        }
    }
    // optional name line
    *line += 1;
    let mut text = line_string(stream);
    if stream.eof() {
        return Ok(None);
    }
    let mut res = RWMol::default();
    if text.starts_with("NAME") {
        res.set_prop(cp::_Name, text[4..].trim().to_string(), false)?;
        *line += 1;
        text = line_string(stream);
        if stream.eof() {
            return Ok(Some(res));
        }
    }
    if text.starts_with("PROP") {
        *line += 1;
        text = line_string(stream);
        if stream.eof() {
            return Ok(Some(res));
        }
    }

    // the counts line
    let split = split_compress(&text);
    let n_atoms = cast_unsigned(split[0])?;
    let n_bonds = cast_unsigned(split.get(1).copied().unwrap_or(""))?;
    let mut conf = Conformer::new(n_atoms as usize);
    conf.set_id(0);
    for _ in 0..n_atoms {
        *line += 1;
        let text = line_string(stream);
        if stream.eof() {
            return Err(file_parse_error("EOF hit while reading atoms."));
        }
        parse_atom_line(&text, *line, &mut res, &mut conf)?;
    }
    res.add_conformer(conf, true)?;

    for _ in 0..n_bonds {
        *line += 1;
        let text = line_string(stream);
        if stream.eof() {
            return Err(file_parse_error("EOF hit while reading bonds."));
        }
        parse_bond_line(&text, *line, &mut res)?;
    }

    *line += 1;
    let text = line_string(stream);
    if stream.eof() {
        return Ok(Some(res));
    }
    let mut n_confs = 0;
    if text.starts_with("CONFS") {
        let split = split_compress(&text);
        n_confs = cast_unsigned(split.get(1).copied().unwrap_or(""))?;
    }
    for i in 0..n_confs {
        let mut conf = parse_conf_data(stream, line, &mut res, i + 1)?;
        if i > 0 || !params.skip_first_conf {
            conf.set_id(i + 1);
            res.add_conformer(conf, true)?;
        }
        // there should be a blank line
        *line += 1;
        let text = line_string(stream);
        if !stream.eof() && !text.trim_matches([' ', '\t', '\r', '\n']).is_empty() {
            return Err(file_parse_error(
                "Found a non-blank line between conformers.",
            ));
        }
    }
    if params.sanitize {
        sanitize_mol(&mut res, SANITIZE_ALL).map_err(|(_, e)| e)?;
    }
    Ok(Some(res))
}

/// `MolFromTPLDataStream` on a string; `None` if the data ends early.
pub fn mol_from_tpl_block(block: &str, params: &TplParserParams) -> Result<Option<RWMol>> {
    let mut stream = LineStream::new(block.as_bytes());
    let mut line = 0;
    mol_from_tpl_stream(&mut stream, &mut line, params)
}

/// `MolFromTPLFile`
pub fn mol_from_tpl_file(file_name: &str, params: &TplParserParams) -> Result<Option<RWMol>> {
    let data = std::fs::read(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    if data.is_empty() {
        return Ok(None);
    }
    mol_from_tpl_block(&String::from_utf8_lossy(&data), params)
}

// ------------------------------------------------------------------ writer

fn write_atom(mol: &ROMol, atom_id: u32, conf: &Conformer, dest: &mut String, charge_prop: &str) {
    let atom = mol.atom(atom_id);
    dest.push_str(&format!(
        "{} {} {}",
        atom_id + 1,
        atom.symbol(),
        atom.formal_charge()
    ));
    let partial = atom
        .get_prop(charge_prop)
        .map_or_else(|| "0.0".to_string(), |v| v.to_rd_string());
    dest.push_str(&format!(" {partial}"));
    let pos = conf.positions()[atom_id as usize];
    dest.push_str(&format!(
        " {} {} {}",
        format_g(100.0 * pos.x, 6),
        format_g(100.0 * pos.y, 6),
        format_g(100.0 * pos.z, 6)
    ));
    let nbrs: Vec<u32> = mol.atom_neighbors(atom_id).collect();
    dest.push_str(&format!(" {}", nbrs.len()));
    for n in nbrs {
        dest.push_str(&format!(" {}", n + 1));
    }
    // upstream: "FIX: get this right"
    dest.push_str(" U\n");
}

fn write_bond(mol: &ROMol, bond_id: u32, dest: &mut String) {
    let bond = mol.bond(bond_id);
    let label = match bond.bond_type() {
        BondType::Single if bond.is_aromatic() => "1.5",
        BondType::Single => "1.0",
        BondType::Double if bond.is_aromatic() => "1.5",
        BondType::Double => "2.0",
        BondType::Aromatic => "1.5",
        BondType::Triple => "3.0",
        other => {
            crate::rd_log_warning!(
                "TPL files only support single, double, aromatic, and triple bonds."
            );
            crate::rd_log_warning!(
                "Bond of with type {} written as single in output.",
                other as u8
            );
            "1.0"
        }
    };
    dest.push_str(&format!(
        "{} {} {} {} 0 0\n",
        bond_id + 1,
        label,
        bond.begin_atom_idx() + 1,
        bond.end_atom_idx() + 1
    ));
}

/// `MolToTPLText`: an empty string (and an error log) for molecules without
/// conformers.
pub fn mol_to_tpl_text(
    mol: &ROMol,
    partial_charge_prop: &str,
    write_first_conf_twice: bool,
) -> String {
    if mol.num_conformers() == 0 {
        crate::rd_log_error!("Cannot write molecules with no conformers to TPL files\n");
        return String::new();
    }
    let mut res = String::from("BioCAD format, all rights reserved\nOutput from RDKit\n");
    let name = match mol.get_prop(cp::_Name) {
        Some(v) => v.to_rd_string(),
        None => {
            crate::rd_log_warning!("Molecule has no name; arbitrary name assigned.\n");
            "Unnamed molecule".to_string()
        }
    };
    res.push_str(&format!("NAME {name}\nPROP 7 1\n"));
    res.push_str(&format!("{} {}\n", mol.num_atoms(), mol.num_bonds()));
    let confs = mol.conformers();
    for i in 0..mol.num_atoms() {
        write_atom(mol, i, &confs[0], &mut res, partial_charge_prop);
    }
    for i in 0..mol.num_bonds() {
        write_bond(mol, i, &mut res);
    }
    res.push_str(&format!("CONFS {}\n", confs.len() - 1));
    let first = usize::from(!write_first_conf_twice);
    for (k, conf) in confs.iter().enumerate().skip(first) {
        res.push_str(&format!("NAME conformer_{}\n", conf.id()));
        for p in conf.positions() {
            res.push_str(&format!(
                " {} {} {}\n",
                format_g(100.0 * p.x, 6),
                format_g(100.0 * p.y, 6),
                format_g(100.0 * p.z, 6)
            ));
        }
        if k + 1 != confs.len() {
            res.push('\n');
        }
    }
    res
}

/// `MolToTPLFile`
pub fn mol_to_tpl_file(
    mol: &ROMol,
    file_name: &str,
    partial_charge_prop: &str,
    write_first_conf_twice: bool,
) -> Result<()> {
    std::fs::write(
        file_name,
        mol_to_tpl_text(mol, partial_charge_prop, write_first_conf_twice),
    )
    .map_err(|_| RDKitError::BadFile(format!("Bad output file {file_name}")))
}
