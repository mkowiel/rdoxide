//! PDB file parser.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/PDBParser.cpp` (2026.09.1).
//!
//! Copyright (C) 2013-2024 Greg Landrum and NextMove Software.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! The parser works on the raw bytes of the file and reads fixed columns, as
//! upstream does; reads past the end of the data yield NUL like the C string
//! the original walks over.

use std::collections::HashMap;

use super::proximity_bonds::{
    CTD_IGNORE_H_H_CONTACTS, connect_the_dots, is_blacklisted_pair,
    standard_pdb_residue_bond_orders,
};
use super::utils::{self, is_bad_cast, to_double_field, to_int};
use crate::geometry::Point3D;
use crate::graphmol::atom::ValenceType;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::BondType;
use crate::graphmol::chirality::from_coords::assign_chiral_types_from_3d;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::monomer_info::{AtomMonomerInfo, AtomMonomerType, PdbResidueData};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result};

/// `PDBParserParams`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PdbParserParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    /// bit 1: keep alternate locations and pseudo atoms; bit 2: one molecule
    /// per MODEL; bit 8: assign standard residue bond orders without
    /// proximity bonding
    pub flavor: u32,
    pub proximity_bonding: bool,
}

impl Default for PdbParserParams {
    fn default() -> Self {
        PdbParserParams {
            sanitize: true,
            remove_hs: true,
            flavor: 0,
            proximity_bonding: true,
        }
    }
}

fn file_parse_error(message: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(message.into())
}

/// One record of the file: a view on the whole buffer starting at the
/// record, with the record length.
#[derive(Clone, Copy)]
struct Line<'a> {
    buf: &'a [u8],
    len: usize,
}

impl Line<'_> {
    /// `ptr[i]`; NUL past the end of the data.
    fn c(&self, i: usize) -> u8 {
        self.buf.get(i).copied().unwrap_or(0)
    }

    /// `std::string(ptr + pos, n)`, clamped to the data.
    fn sub(&self, pos: usize, n: usize) -> &[u8] {
        let end = self.buf.len().min(pos + n);
        self.buf.get(pos..end).unwrap_or_default()
    }

    fn text(&self, pos: usize, n: usize) -> String {
        utils::lossy(self.sub(pos, n))
    }

    fn starts_with(&self, prefix: &[u8]) -> bool {
        self.buf.starts_with(prefix)
    }

    fn int(&self, pos: usize, n: usize) -> Result<i32> {
        to_int(self.sub(pos, n), true)
    }

    fn double(&self, pos: usize, n: usize) -> Result<f64> {
        let field = self.sub(pos, n);
        to_double_field(field, 0, field.len())
    }
}

/// Parsing state that spans records.
#[derive(Default)]
struct State {
    /// PDB serial number -> atom index
    atoms: HashMap<i32, u32>,
    /// bond index -> bitmap of the CONECT entries seen for it
    seen: HashMap<u32, u8>,
}

/// `PDBAtomFromSymbol`: `D` and `T` are hydrogen isotopes.
fn atom_from_symbol(symbol: &str) -> Result<Option<Atom>> {
    let isotope = match symbol {
        "D" => Some(2),
        "T" => Some(3),
        _ => None,
    };
    if let Some(isotope) = isotope {
        let mut atom = Atom::new(1);
        atom.set_isotope(isotope);
        return Ok(Some(atom));
    }
    let num = PeriodicTable::get().get_atomic_number(symbol)?;
    Ok((num > 0).then(|| Atom::new(num)))
}

/// The element symbol from columns 77-78, if they hold one.
fn symbol_from_element_columns(l: &Line<'_>) -> String {
    let (c76, c77) = (l.c(76), l.c(77));
    let mut symbol = String::new();
    if l.len >= 78 {
        if c76.is_ascii_uppercase() {
            symbol.push(c76 as char);
            if c77.is_ascii_uppercase() {
                symbol.push(c77.to_ascii_lowercase() as char);
            } else if c77.is_ascii_lowercase() {
                symbol.push(c77 as char);
            }
        } else if c76 == b' ' && c77.is_ascii_uppercase() {
            symbol.push(c77 as char);
        }
    } else if l.len == 77 && c76.is_ascii_uppercase() {
        symbol.push(c76 as char);
    }
    symbol
}

/// The element symbol deduced from the atom name (columns 13-16).
fn symbol_from_atom_name(l: &Line<'_>) -> String {
    let (c12, c13, c14) = (l.c(12), l.c(13), l.c(14));
    let mut symbol = String::new();
    if !c13.is_ascii_uppercase() {
        return symbol;
    }
    if c12 == b' ' {
        symbol.push(c13 as char);
        if c14.is_ascii_lowercase() {
            symbol.push(c14 as char);
        }
    } else if c12.is_ascii_uppercase() {
        symbol.push(c12 as char);
        symbol.push(c13.to_ascii_lowercase() as char);
        if c12 == b'H' && l.c(0) == b'A' {
            // No He, Hf, Hg, Ho or Hs in ATOM records
            symbol = "H".to_string();
        }
    } else if c12.is_ascii_digit() {
        symbol.push(c13 as char);
    }
    symbol
}

/// The formal charge from columns 79-80.
fn charge_from_columns(l: &Line<'_>) -> i32 {
    let (c78, c79) = (l.c(78), l.c(79));
    let digit = |c: u8| (b'1'..=b'9').contains(&c).then(|| i32::from(c - b'0'));
    match c78 {
        b'1'..=b'9' => {
            let d = i32::from(c78 - b'0');
            if c79 == b'-' {
                -d
            } else if matches!(c79, b'+' | b' ' | 0) {
                d
            } else {
                0
            }
        }
        b'+' => digit(c79).unwrap_or(match c79 {
            b'+' => 2,
            b'0' => 0,
            _ => 1,
        }),
        // (a digit after '-' gives a positive charge, as upstream)
        b'-' => digit(c79).unwrap_or(match c79 {
            b'-' => -2,
            b'0' => 0,
            _ => -1,
        }),
        b' ' => digit(c79).unwrap_or(match c79 {
            b'+' => 1,
            b'-' => -1,
            _ => 0,
        }),
        _ => 0,
    }
}

/// `x`, `y`, `z` of an ATOM/HETATM record (missing fields stay 0).
fn coordinates(l: &Line<'_>, error: impl Fn() -> RDKitError) -> Result<Point3D> {
    let mut pos = Point3D::default();
    let read = |at: usize| -> Result<f64> {
        l.double(at, 8)
            .map_err(|e| if is_bad_cast(&e) { error() } else { e })
    };
    pos.x = read(30)?;
    if l.len >= 46 {
        pos.y = read(38)?;
    }
    if l.len >= 54 {
        pos.z = read(46)?;
    }
    Ok(pos)
}

/// `PDBAtomLine`
fn atom_line(mol: &mut RWMol, l: &Line<'_>, flavor: u32, state: &mut State) -> Result<()> {
    if l.len < 16 {
        return Ok(());
    }
    if flavor & 1 == 0 {
        // Ignore alternate locations of atoms.
        if l.len >= 17 && !matches!(l.c(16), b' ' | b'A' | b'1') {
            return Ok(());
        }
        // Ignore XPLOR pseudo atoms
        if l.len >= 54 && l.sub(30, 24) == b"9999.0009999.0009999.000" {
            return Ok(());
        }
        // Ignore NMR pseudo atoms
        if l.c(12) == b' ' && l.c(13) == b'Q' {
            return Ok(());
        }
        // Ignore PDB dummy residues
        if l.len >= 20 && l.sub(18, 3) == b"DUM" {
            return Ok(());
        }
    }

    let serial_text = l.sub(6, 5);
    let serial = to_int(serial_text, true).map_err(|e| {
        if is_bad_cast(&e) {
            file_parse_error(format!(
                "Non-integer PDB serial number {}",
                utils::lossy(serial_text)
            ))
        } else {
            e
        }
    })?;

    let mut atom = None;
    let symbol = symbol_from_element_columns(l);
    if !symbol.is_empty() {
        atom = atom_from_symbol(&symbol)?;
    }
    if atom.is_none() {
        let symbol = symbol_from_atom_name(l);
        if !symbol.is_empty() {
            atom = atom_from_symbol(&symbol)?;
        }
    }
    let Some(mut atom) = atom else {
        return Err(file_parse_error(format!(
            "Cannot determine element for PDB atom #{serial}"
        )));
    };

    if l.len >= 79 {
        let charge = charge_from_columns(l);
        if charge != 0 {
            atom.set_formal_charge(charge);
        }
    }

    let idx = mol.add_atom(atom, true);
    state.atoms.insert(serial, idx);

    if l.len >= 38 {
        let pos = coordinates(l, || {
            file_parse_error(format!("Problem with coordinates for PDB atom #{serial}"))
        })?;
        if mol.num_conformers() == 0 {
            let mut conf = Conformer::new(mol.num_atoms() as usize);
            conf.set_3d(pos.z != 0.0);
            conf.set_id(0);
            mol.add_conformer(conf, false)?;
        } else if pos.z != 0.0 {
            mol.conformer_mut(-1)?.set_3d(true);
        }
        mol.conformer_mut(-1)?.set_atom_pos(idx, pos)?;
    }

    let residue_name = if l.len >= 20 {
        l.text(17, 3)
    } else {
        "UNL".into()
    };
    let single = |len_needed: usize, pos: usize| {
        if l.len >= len_needed {
            l.text(pos, 1)
        } else {
            " ".into()
        }
    };
    let alt_loc = single(17, 16);
    let chain_id = single(22, 21);
    let insertion_code = single(27, 26);
    let bad_field = |what: &str| {
        let message = format!("Problem with {what} for PDB atom #{serial}");
        move |e: RDKitError| {
            if is_bad_cast(&e) {
                file_parse_error(message.clone())
            } else {
                e
            }
        }
    };
    let residue_number = if l.len >= 26 {
        l.int(22, 4).map_err(bad_field("residue number"))?
    } else {
        1
    };
    let occupancy = if l.len >= 60 {
        l.double(54, 6).map_err(bad_field("occupancy"))?
    } else {
        1.0
    };
    let temp_factor = if l.len >= 66 {
        l.double(60, 6).map_err(bad_field("temperature factor"))?
    } else {
        0.0
    };
    let info = AtomMonomerInfo {
        monomer_type: AtomMonomerType::PdbResidue,
        name: l.text(12, 4),
        residue_name,
        residue_number,
        chain_id,
        monomer_class: String::new(),
        pdb: Some(PdbResidueData {
            serial_number: serial as u32,
            alt_loc,
            insertion_code,
            occupancy,
            temp_factor,
            is_hetero_atom: l.c(0) == b'H',
            ..PdbResidueData::default()
        }),
    };
    mol.atom_mut(idx).set_monomer_info(Some(info));
    Ok(())
}

/// Tracks repeated CONECT entries to encode bond orders. `seen` keeps a
/// nibble per direction (low: `src < dst`, high: `src > dst`) with one bit per
/// entry; the n-th repeat in one direction raises the order to n+1 unless the
/// other direction already did. Returns the new state and the new order.
fn update_bond_order(seen: u8, forward: bool) -> (u8, Option<BondType>) {
    let (shift, other) = if forward { (0, 4) } else { (4, 0) };
    let (bit, order) = match ((seen >> shift) & 0x0f).count_ones() {
        1 => (0x02u8, BondType::Double),
        2 => (0x04, BondType::Triple),
        3 => (0x08, BondType::Quadruple),
        _ => return (seen, None),
    };
    let raised_by_other = seen & (bit << other) != 0;
    (seen | (bit << shift), (!raised_by_other).then_some(order))
}

/// `PDBBondLine`
fn bond_line(mol: &mut RWMol, l: &Line<'_>, state: &mut State) -> Result<()> {
    if l.len < 16 {
        return Ok(());
    }
    let src_text = l.sub(6, 5);
    let failure = || {
        file_parse_error(format!(
            "Problem with CONECT record for PDB atom #{}",
            utils::lossy(src_text)
        ))
    };
    let src = to_int(src_text, true).map_err(|e| if is_bad_cast(&e) { failure() } else { e })?;
    let Some(&src_idx) = state.atoms.get(&src) else {
        return Ok(());
    };

    let len = l.len.min(41);
    let mut pos = 11;
    while pos + 5 <= len {
        if l.sub(pos, 5) == b"     " {
            break;
        }
        let dst =
            to_int(l.sub(pos, 5), true).map_err(|e| if is_bad_cast(&e) { failure() } else { e })?;
        pos += 5;
        let dst_idx = match state.atoms.get(&dst) {
            Some(&d) if dst != src => d,
            _ => continue,
        };
        match mol.get_bond_between_atoms(src_idx, dst_idx)? {
            Some(bond) if bond.bond_type() != BondType::Zero => {
                let bond_idx = bond.idx();
                let seen = state.seen.get(&bond_idx).copied().unwrap_or(0);
                let (seen, order) = update_bond_order(seen, src < dst);
                state.seen.insert(bond_idx, seen);
                if let Some(order) = order {
                    mol.bond_mut(bond_idx).set_bond_type(order);
                }
            }
            None => {
                // Bonds in PDB file are explicit; if they are not sanitize
                // friendly, set their order to zero
                let order = if is_blacklisted_pair(mol.atom(src_idx), mol.atom(dst_idx)) {
                    BondType::Zero
                } else {
                    BondType::Single
                };
                mol.add_bond(src_idx, dst_idx, order)?;
                let bond_idx = mol.num_bonds() - 1;
                state
                    .seen
                    .insert(bond_idx, if src < dst { 0x01 } else { 0x10 });
            }
            Some(_) => break,
        }
    }
    Ok(())
}

/// `PDBTitleLine`
fn title_line(mol: &mut RWMol, l: &Line<'_>, mut len: usize) -> Result<()> {
    while len > 0 && l.c(len - 1) == b' ' {
        len -= 1;
    }
    if len > 0 && l.c(len - 1) == b';' {
        len -= 1;
    }
    let title = if len > 21 && l.sub(10, 11) == b" MOLECULE: " {
        l.text(21, len - 21)
    } else if len > 10 {
        l.text(10, len - 10)
    } else {
        String::new()
    };
    if !title.is_empty() {
        mol.set_prop(cp::_Name, title, false)?;
    }
    Ok(())
}

/// `PDBConformerLine`: coordinates of an additional model.
fn conformer_line(
    mol: &mut RWMol,
    l: &Line<'_>,
    conf_id: &mut Option<u32>,
    atom_idx: &mut u32,
) -> Result<()> {
    if l.len < 38 {
        return Ok(());
    }
    let pos = coordinates(l, || {
        file_parse_error("Problem with multi-conformer coordinates")
    })?;
    match conf_id {
        Some(id) if *atom_idx != 0 => {
            if pos.z != 0.0 {
                mol.conformer_mut(*id as i32)?.set_3d(true);
            }
        }
        _ => {
            let mut conf = Conformer::new(mol.num_atoms() as usize);
            let id = mol.num_conformers();
            conf.set_id(id);
            conf.set_3d(pos.z != 0.0);
            mol.add_conformer(conf, false)?;
            *conf_id = Some(id);
        }
    }
    if *atom_idx < mol.num_atoms() {
        if let Some(id) = *conf_id {
            mol.conformer_mut(id as i32)?.set_atom_pos(*atom_idx, pos)?;
        }
        *atom_idx += 1;
    }
    Ok(())
}

/// `StandardPDBChiralAtom`: whether a standard atom name in a recognized
/// amino acid should be chiral (avoids chirality on VAL.CG and LEU.CG).
fn standard_pdb_chiral_atom(residue: &str, name: &str) -> bool {
    let alpha = name.starts_with(" C") && name.as_bytes().get(3) == Some(&b' ');
    match residue.get(..3).unwrap_or_default() {
        "GLY" => false,
        // alpha and beta carbons (" CA " and " CB ")
        "ILE" | "THR" => alpha && matches!(name.as_bytes().get(2), Some(b'A' | b'B')),
        "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "GLN" | "GLU" | "HIS" | "LEU" | "LYS" | "MET"
        | "PHE" | "PRO" | "SER" | "TRP" | "TYR" | "VAL" => {
            alpha && name.as_bytes().get(2) == Some(&b'A')
        }
        _ => false,
    }
}

fn standard_pdb_residue_chirality(mol: &mut RWMol) {
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom(idx);
        if atom.chiral_tag() == ChiralType::Unspecified {
            continue;
        }
        let Some(info) = atom.monomer_info() else {
            continue;
        };
        let hetero = info.pdb.as_ref().is_some_and(|p| p.is_hetero_atom);
        if info.monomer_type == AtomMonomerType::PdbResidue
            && !hetero
            && !standard_pdb_chiral_atom(&info.residue_name, &info.name)
        {
            let atom = mol.atom_mut(idx);
            atom.set_chiral_tag(ChiralType::Unspecified);
            atom.clear_prop(cp::_CIPCode);
        }
    }
}

/// `BasicPDBCleanup`
fn basic_pdb_cleanup(mol: &mut RWMol) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        mol.calc_explicit_valence(idx, false)?;
        let atom = mol.atom(idx);
        // correct four-valent neutral N -> N+ (github #1029)
        if atom.atomic_num() == 7
            && atom.formal_charge() == 0
            && atom.valence(ValenceType::Explicit)? == 4
        {
            mol.atom_mut(idx).set_formal_charge(1);
        }
    }
    Ok(())
}

/// Splits off the next record: `(start of the next record, record length)`.
/// Like upstream, the search for the line end starts at the second
/// character, so a leading line break belongs to the record.
fn next_record(buf: &[u8], start: usize) -> (usize, usize) {
    let mut next = start + 1;
    loop {
        match buf.get(next).copied().unwrap_or(0) {
            b'\r' => {
                let len = next - start;
                next += if buf.get(next + 1) == Some(&b'\n') {
                    2
                } else {
                    1
                };
                return (next, len);
            }
            b'\n' => return (next + 1, next - start),
            0 => return (next, next - start),
            _ => next += 1,
        }
    }
}

/// `parsePdbBlock`
fn parse_pdb_block(text: &str, params: &PdbParserParams) -> Result<Option<RWMol>> {
    // the C string ends at the first NUL
    let bytes = text.as_bytes();
    let bytes = &bytes[..bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len())];

    let mut state = State::default();
    let mut mol: Option<RWMol> = None;
    let mut multi_conformer = false;
    let mut conf_id: Option<u32> = None;
    let mut conf_atom_idx = 0u32;

    let mut start = 0;
    while start < bytes.len() {
        let (next, len) = next_record(bytes, start);
        let rec = Line {
            buf: &bytes[start..],
            len,
        };
        if rec.starts_with(b"ATOM  ") || rec.starts_with(b"HETATM") {
            if multi_conformer {
                if let Some(m) = mol.as_mut() {
                    conformer_line(m, &rec, &mut conf_id, &mut conf_atom_idx)?;
                }
            } else {
                atom_line(
                    mol.get_or_insert_with(RWMol::new),
                    &rec,
                    params.flavor,
                    &mut state,
                )?;
            }
        } else if rec.starts_with(b"CONECT") {
            if let (Some(m), false) = (mol.as_mut(), multi_conformer) {
                bond_line(m, &rec, &mut state)?;
            }
        } else if rec.starts_with(b"COMPND") {
            let m = mol.get_or_insert_with(RWMol::new);
            if len > 10 && (rec.c(9) == b' ' || rec.sub(9, 12) == b"2 MOLECULE: ") {
                title_line(m, &rec, len)?;
            }
        } else if rec.starts_with(b"HEADER") {
            let m = mol.get_or_insert_with(RWMol::new);
            title_line(m, &rec, len.min(50))?;
        } else if rec.starts_with(b"ENDMDL") {
            if mol.is_none() {
                break;
            }
            multi_conformer = true;
            conf_atom_idx = 0;
            conf_id = None;
        }
        start = next;
    }

    let Some(mut mol) = mol else {
        return Ok(None);
    };

    if params.proximity_bonding {
        connect_the_dots(&mut mol, CTD_IGNORE_H_H_CONTACTS)?;
    }
    // flavor & 8 doesn't encode double bonds
    if params.proximity_bonding || params.flavor & 8 != 0 {
        standard_pdb_residue_bond_orders(&mut mol);
    }

    basic_pdb_cleanup(&mut mol)?;

    if params.sanitize {
        if params.remove_hs {
            remove_hs(&mut mol, &RemoveHsParameters::default(), true)?;
        } else {
            sanitize_mol(&mut mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
        }
    } else {
        // we need some properties for the chiral setup
        mol.update_property_cache(false)?;
    }

    // set tetrahedral chirality from 3D co-ordinates
    assign_chiral_types_from_3d(&mut mol, -1, true)?;
    standard_pdb_residue_chirality(&mut mol);
    Ok(Some(mol))
}

/// `v2::FileParsers::MolFromPDBBlock`: `Ok(None)` when the block holds no
/// atoms or records.
pub fn mol_from_pdb_block(text: &str, params: &PdbParserParams) -> Result<Option<RWMol>> {
    parse_pdb_block(text, params)
}

/// `v2::FileParsers::MolFromPDBDataStream`: reads up to the `END` record (or
/// `ENDMDL` with flavor bit 2) and parses what was read.
pub fn mol_from_pdb_data_stream(text: &str, params: &PdbParserParams) -> Result<Option<RWMol>> {
    let mut buffer = String::new();
    for line in text.split('\n') {
        buffer.push_str(line);
        buffer.push('\n');
        let b = line.as_bytes();
        let at = |i: usize| b.get(i).copied().unwrap_or(0);
        if b.starts_with(b"END") && matches!(at(3), b' ' | b'\r' | b'\n' | 0) {
            break;
        }
        if params.flavor & 2 != 0 && b.starts_with(b"ENDMDL") {
            break;
        }
    }
    mol_from_pdb_block(&buffer, params)
}

/// `v2::FileParsers::MolFromPDBFile`
pub fn mol_from_pdb_file(file_name: &str, params: &PdbParserParams) -> Result<Option<RWMol>> {
    let data = std::fs::read(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    mol_from_pdb_data_stream(&String::from_utf8_lossy(&data), params)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::file_parsers::pdb_writer::mol_to_pdb_block;

    const GLY: &str = "\
ATOM      1  N   GLY A   1      -0.966   0.493   1.500  1.00  0.00           N
ATOM      2  CA  GLY A   1       0.257   0.418   0.692  1.00  0.00           C
ATOM      3  C   GLY A   1      -0.094   0.017  -0.716  1.00  0.00           C
ATOM      4  O   GLY A   1      -1.056  -0.682  -0.923  1.00  0.00           O
END
";

    #[test]
    fn parses_atoms_bonds_and_residue_info() {
        let mol = mol_from_pdb_block(GLY, &PdbParserParams::default())
            .unwrap()
            .unwrap();
        assert_eq!(mol.num_atoms(), 4);
        assert_eq!(mol.num_bonds(), 3);
        let info = mol.atom(1).monomer_info().unwrap();
        assert_eq!(
            (info.name.as_str(), info.residue_name.as_str()),
            (" CA ", "GLY")
        );
        assert_eq!(info.pdb.as_ref().unwrap().serial_number, 2);
        // the standard-residue C=O double bond
        assert_eq!(mol.bond(2).bond_type(), BondType::Double);
    }

    #[test]
    fn writer_round_trips_the_atom_lines() {
        let mol = mol_from_pdb_block(GLY, &PdbParserParams::default())
            .unwrap()
            .unwrap();
        let block = mol_to_pdb_block(&mol, -1, 0).unwrap();
        for line in GLY.lines().take(4) {
            assert!(block.contains(line), "missing {line}");
        }
        assert!(block.ends_with("END\n"));
    }

    #[test]
    fn empty_input_gives_no_molecule() {
        assert!(
            mol_from_pdb_block("", &PdbParserParams::default())
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn repeated_conect_entries_raise_the_bond_order() {
        assert_eq!(
            update_bond_order(0x01, true),
            (0x03, Some(BondType::Double))
        );
        assert_eq!(
            update_bond_order(0x03, true),
            (0x07, Some(BondType::Triple))
        );
        // the other direction already raised it
        assert_eq!(update_bond_order(0x21, true), (0x23, None));
        assert_eq!(
            update_bond_order(0x10, false),
            (0x30, Some(BondType::Double))
        );
    }
}
