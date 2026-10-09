//! XYZ file parser and writer.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/{XYZFileParser,XYZFileWriter}.cpp` (2026.09.1).
//!
//! Copyright (C) 2022 Rachel Walker and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result};

const WHITESPACE: [char; 2] = [' ', '\t'];

fn file_parse_error(message: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(message.into())
}

/// `std::getline` over a string, tracking the stream's `eof` state.
struct Lines<'a> {
    rest: &'a str,
    eof: bool,
}

impl<'a> Lines<'a> {
    fn new(text: &'a str) -> Self {
        Lines {
            rest: text,
            eof: false,
        }
    }

    /// The next line without its terminator (and without a trailing `\r`).
    fn get_line(&mut self) -> &'a str {
        let line = match self.rest.find('\n') {
            Some(pos) => {
                let line = &self.rest[..pos];
                self.rest = &self.rest[pos + 1..];
                line
            }
            None => {
                self.eof = true;
                std::mem::take(&mut self.rest)
            }
        };
        line.strip_suffix('\r').unwrap_or(line)
    }
}

/// `ParseXYZFileAtomLine`: `symbol x y z` (extra trailing fields are ignored
/// only as part of the z field, as upstream).
fn parse_atom_line(line: &str, line_no: usize) -> Result<(Atom, Point3D)> {
    let missing = || file_parse_error(format!("Missing coordinates on line {line_no}\n"));
    // delimiters alternate between the start of a field and its end
    let mut delims = [0usize; 8];
    let mut prev = 0;
    for i in 0..7 {
        let rest = &line[prev..];
        let found = if i % 2 == 0 {
            rest.find(|c| !WHITESPACE.contains(&c))
        } else {
            rest.find(WHITESPACE)
        };
        delims[i] = prev + found.ok_or_else(missing)?;
        prev = delims[i];
    }
    delims[7] = line
        .rfind(|c| !WHITESPACE.contains(&c))
        .map_or(0, |p| p + 1);
    let coordinate = |from: usize, to: usize| -> Result<f64> {
        let text = line.get(from..to).unwrap_or_default();
        text.parse::<f64>().map_err(|_| {
            file_parse_error(format!(
                "Cannot convert '{text}' to double on line {line_no}\n"
            ))
        })
    };
    let pos = Point3D::new(
        coordinate(delims[2], delims[3])?,
        coordinate(delims[4], delims[5])?,
        coordinate(delims[6], delims[7])?,
    );
    let mut symbol = line[delims[0]..delims[1]].to_string();
    if symbol.len() == 2 && symbol.as_bytes()[1].is_ascii_uppercase() {
        symbol = format!("{}{}", &symbol[..1], symbol[1..].to_ascii_lowercase());
    }
    let atomic_num = PeriodicTable::get()
        .get_atomic_number(&symbol)
        .map_err(|e| file_parse_error(e.to_string()))?;
    Ok((Atom::new(atomic_num), pos))
}

/// `MolFromXYZDataStream`
fn mol_from_xyz_data(text: &str) -> Result<RWMol> {
    let mut lines = Lines::new(text);
    let num = lines.get_line();
    let num_atoms: usize = num.trim_matches(WHITESPACE).parse().map_err(|_| {
        file_parse_error(format!(
            "Unable to recognize the number of atoms: cannot convert '{num}' to unsigned int on line 0\n"
        ))
    })?;
    let comment = lines.get_line();
    let mut mol = ROMol::new();
    if num_atoms > 0 {
        let mut conf = Conformer::new(num_atoms);
        if !comment.is_empty() {
            mol.set_prop("_FileComments", comment, false)?;
        }
        for i in 0..num_atoms {
            if lines.eof {
                return Err(file_parse_error("EOF hit while reading atoms"));
            }
            let (atom, pos) = parse_atom_line(lines.get_line(), i + 2)?;
            let idx = mol.add_atom(atom, false);
            conf.set_atom_pos(idx, pos)?;
        }
        mol.add_conformer(conf, false)?;
    }
    while !lines.eof {
        if lines.get_line().contains(|c| !WHITESPACE.contains(&c)) {
            return Err(file_parse_error("More lines than expected\n"));
        }
    }
    Ok(mol)
}

/// `MolFromXYZBlock`: `None` for an empty block.
pub fn mol_from_xyz_block(block: &str) -> Result<Option<RWMol>> {
    if block.is_empty() {
        return Ok(None);
    }
    mol_from_xyz_data(block).map(Some)
}

/// `MolFromXYZFile`: `None` for an empty file.
pub fn mol_from_xyz_file(file_name: &str) -> Result<Option<RWMol>> {
    let data = std::fs::read(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    mol_from_xyz_block(&String::from_utf8_lossy(&data))
}

/// `MolToXYZBlock`: the atoms of conformer `conf_id` with `precision` decimals.
/// An empty string (and an error log) for molecules without conformers.
pub fn mol_to_xyz_block(mol: &ROMol, conf_id: i32, precision: usize) -> Result<String> {
    if mol.num_conformers() == 0 {
        crate::rd_log_error!("Cannot write molecules with no conformers to XYZ block\n");
        return Ok(String::new());
    }
    let conf = mol.conformer(conf_id)?;
    let width = 5 + precision;
    let mut out = format!("{}\n", mol.num_atoms());
    if let Some(name) = mol.get_prop(cp::_Name) {
        let name = name.to_rd_string();
        out.push_str(name.split('\n').next().unwrap_or_default());
    }
    out.push('\n');
    for (i, atom) in mol.atoms().iter().enumerate() {
        let p = conf.atom_pos(i as u32)?;
        out.push_str(&format!(
            "{:<3} {:>width$.precision$} {:>width$.precision$} {:>width$.precision$}\n",
            atom.symbol(),
            p.x,
            p.y,
            p.z,
        ));
    }
    Ok(out)
}

/// `MolToXYZFile`
pub fn mol_to_xyz_file(mol: &ROMol, file_name: &str, conf_id: i32, precision: usize) -> Result<()> {
    let block = mol_to_xyz_block(mol, conf_id, precision)?;
    std::fs::write(file_name, block)
        .map_err(|_| RDKitError::BadFile(format!("Bad output file {file_name}")))
}
