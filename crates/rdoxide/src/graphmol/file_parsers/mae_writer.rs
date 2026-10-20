//! Maestro (`.mae`) writer.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MaeWriter.cpp` (2026.09.1).
//! The block layout (property names by type in the order bool, real, int,
//! string, then the values; indexed blocks with one row per atom or bond;
//! two spaces of indentation per level) is that of `maeparser`'s writer.
//!
//! Copyright (C) 2023 Schrodinger, LLC and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::mae_mol_supplier::{
    MAE_BOND_DATIVE_MARK, MAE_BOND_PARITY, MAE_RGROUP_LABEL, MAE_STEREO_STATUS, PDB_ATOM_NAME,
    PDB_CHAIN_NAME, PDB_INSERTION_CODE, PDB_OCCUPANCY, PDB_RESIDUE_NAME, PDB_RESIDUE_NUMBER,
    PDB_TFACTOR,
};
use super::mae_reader::constants::*;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::dict::COMPUTED_PROP_NAME;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_warning};

/// A column of an indexed block; `None` is an undefined value (`<>`).
#[derive(Debug, Clone)]
enum Column {
    Bool(Vec<Option<bool>>),
    Int(Vec<Option<i32>>),
    Real(Vec<Option<f64>>),
    Str(Vec<Option<String>>),
}

/// A value set on a property of an indexed block.
enum Value {
    Bool(bool),
    Int(i32),
    Real(f64),
    Str(String),
}

#[derive(Debug, Default)]
struct IndexedBlock {
    size: usize,
    columns: BTreeMap<String, Column>,
}

impl IndexedBlock {
    fn new(size: usize) -> Self {
        IndexedBlock {
            size,
            columns: BTreeMap::new(),
        }
    }

    /// `setPropertyValue`: the type of the value must be the one named by the
    /// first letter of the property.
    fn set(&mut self, name: &str, idx: u32, value: Value) -> Result<()> {
        let (tag, wanted) = match value {
            Value::Bool(_) => (b'b', "a boolean"),
            Value::Int(_) => (b'i', "an integer"),
            Value::Real(_) => (b'r', "a real"),
            Value::Str(_) => (b's', "a string"),
        };
        if name.as_bytes()[0] != tag {
            return Err(RDKitError::Runtime(format!(
                "Property '{name}' is not {wanted} value"
            )));
        }
        let n = self.size;
        let col = self
            .columns
            .entry(name.to_string())
            .or_insert_with(|| match value {
                Value::Bool(_) => Column::Bool(vec![None; n]),
                Value::Int(_) => Column::Int(vec![None; n]),
                Value::Real(_) => Column::Real(vec![None; n]),
                Value::Str(_) => Column::Str(vec![None; n]),
            });
        let i = idx as usize;
        match (col, value) {
            (Column::Bool(c), Value::Bool(v)) => c[i] = Some(v),
            (Column::Int(c), Value::Int(v)) => c[i] = Some(v),
            (Column::Real(c), Value::Real(v)) => c[i] = Some(v),
            (Column::Str(c), Value::Str(v)) => c[i] = Some(v),
            _ => unreachable!("the column type follows the first letter of its name"),
        }
        Ok(())
    }
}

/// The typed properties of the structure block.
#[derive(Debug, Default)]
struct CtBlock {
    bools: BTreeMap<String, bool>,
    reals: BTreeMap<String, f64>,
    ints: BTreeMap<String, i32>,
    strings: BTreeMap<String, String>,
    atoms: IndexedBlock,
    bonds: Option<IndexedBlock>,
}

/// `[birs]_[^_ ]+_.+` (`std::regex_match`; `.` does not match line breaks)
fn is_mmct_property_name(name: &str) -> bool {
    let b = name.as_bytes();
    if b.len() < 5 || !matches!(b[0], b'b' | b'i' | b'r' | b's') || b[1] != b'_' {
        return false;
    }
    let rest = &name[2..];
    let Some(end) = rest.find(['_', ' ']) else {
        return false;
    };
    if end == 0 || !rest[end..].starts_with('_') {
        return false;
    }
    let tail = &rest[end + 1..];
    !tail.is_empty() && !tail.contains(['\n', '\r', '\u{2028}', '\u{2029}'])
}

/// The Maestro name of a property: its own, or prefixed with `<type>_rdkit_`.
fn mae_name(key: &str, prefix: &str) -> String {
    if is_mmct_property_name(key) {
        key.to_string()
    } else {
        format!("{prefix}{key}")
    }
}

/// `copyProperties`: the (non-computed) properties of `origin`, as typed
/// values through the setters.
fn copy_properties(
    origin: &impl RDProps,
    prop_names: &[String],
    mut set: impl FnMut(String, Value) -> Result<()>,
) -> Result<()> {
    for (key, val) in origin.dict().iter() {
        if key == COMPUTED_PROP_NAME {
            continue;
        }
        if !prop_names.is_empty() && !prop_names.iter().any(|p| p == key) {
            continue;
        }
        match val {
            RDValue::Bool(v) => set(mae_name(key, "b_rdkit_"), Value::Bool(*v))?,
            RDValue::Int(_) | RDValue::UnsignedInt(_) => {
                let v = match val {
                    RDValue::Int(v) => *v,
                    RDValue::UnsignedInt(v) => *v as i32,
                    _ => unreachable!(),
                };
                let name = if key == cp::_MolFileRLabel {
                    MAE_RGROUP_LABEL.to_string()
                } else {
                    mae_name(key, "i_rdkit_")
                };
                set(name, Value::Int(v))?;
            }
            RDValue::Double(v) => set(mae_name(key, "r_rdkit_"), Value::Real(*v))?,
            RDValue::Float(v) => set(mae_name(key, "r_rdkit_"), Value::Real(f64::from(*v)))?,
            RDValue::String(v) => set(mae_name(key, "s_rdkit_"), Value::Str(v.clone()))?,
            _ => {
                // vectors and the like are not exported
                rd_log_warning!(
                    "WARNING: the property {key} has an unsupported type, and will not be exported."
                );
            }
        }
    }
    Ok(())
}

fn bond_type_to_order(mol: &ROMol, idx: u32) -> Result<i32> {
    Ok(match mol.bond(idx).bond_type() {
        BondType::Single => 1,
        BondType::Double => 2,
        BondType::Triple => 3,
        BondType::Zero | BondType::Dative => 0,
        _ => {
            return Err(RDKitError::Value(format!(
                "Bond {idx} has a type that is not supported by maeparser."
            )));
        }
    })
}

fn int_prop_is(mol: &ROMol, idx: u32, key: &str, value: i32) -> bool {
    matches!(mol.bond(idx).get_prop(key), Some(RDValue::Int(v)) if *v == value)
}

fn is_double_any_bond(mol: &ROMol, idx: u32) -> bool {
    let b = mol.bond(idx);
    b.bond_type() == BondType::Double
        && (b.stereo() == BondStereo::StereoAny
            || b.bond_dir() == BondDir::EitherDouble
            || int_prop_is(mol, idx, cp::_MolFileBondCfg, 2)
            || int_prop_is(mol, idx, cp::_MolFileBondStereo, 3))
}

/// `copyAtomNumChirality`
fn copy_atom_num_chirality(mol: &ROMol, ct: &mut CtBlock) {
    // tells Schrodinger software that the stereo and chirality are valid
    ct.ints.insert(MAE_STEREO_STATUS.to_string(), 1);
    let mut chiral_atoms = 0;
    for atom in mol.atoms() {
        let label = match atom.chiral_tag() {
            ChiralType::TetrahedralCW => "ANR",
            ChiralType::TetrahedralCCW => "ANS",
            _ => continue,
        };
        chiral_atoms += 1;
        let mut value = format!("{}_{label}", atom.idx() + 1);
        // CIP ranks are unknown: the neighbors are sorted by index
        let mut neighbors: Vec<u32> = mol.atom_neighbors(atom.idx()).collect();
        neighbors.sort_unstable();
        for nb in neighbors {
            let _ = write!(value, "_{}", nb + 1);
        }
        ct.strings
            .insert(format!("{CT_CHIRALITY_PROP_PREFIX}{chiral_atoms}"), value);
    }
}

fn map_mol_properties(mol: &mut ROMol, prop_names: &[String], ct: &mut CtBlock) -> Result<()> {
    // a title is always written
    let name = match mol.get_prop(cp::_Name) {
        Some(RDValue::String(s)) => s.clone(),
        _ => String::new(),
    };
    ct.strings.insert(CT_TITLE.to_string(), name);
    mol.clear_prop(cp::_Name);
    copy_atom_num_chirality(mol, ct);
    copy_properties(&*mol, prop_names, |name, value| {
        match value {
            Value::Bool(v) => ct.bools.insert(name, v).map(|_| ()),
            Value::Int(v) => ct.ints.insert(name, v).map(|_| ()),
            Value::Real(v) => ct.reals.insert(name, v).map(|_| ()),
            Value::Str(v) => ct.strings.insert(name, v).map(|_| ()),
        };
        Ok(())
    })
}

fn map_atoms(mol: &ROMol, prop_names: &[String], conf_id: i32, ct: &mut CtBlock) -> Result<()> {
    let conf = mol.conformer(conf_id)?;
    let block = &mut ct.atoms;
    for atom in mol.atoms() {
        let idx = atom.idx();
        let pos = conf.atom_pos(idx)?;
        block.set(ATOM_X_COORD, idx, Value::Real(pos.x))?;
        block.set(ATOM_Y_COORD, idx, Value::Real(pos.y))?;
        block.set(ATOM_Z_COORD, idx, Value::Real(pos.z))?;
        let mut z = atom.atomic_num() as i32;
        if z == 0 {
            // dummy atom
            z = -2;
        }
        block.set(ATOM_ATOMIC_NUM, idx, Value::Int(z))?;
        block.set(ATOM_FORMAL_CHARGE, idx, Value::Int(atom.formal_charge()))?;
        if let Some(info) = atom.monomer_info() {
            let pdb = info.pdb.clone().unwrap_or_default();
            block.set(PDB_ATOM_NAME, idx, Value::Str(info.name.clone()))?;
            block.set(PDB_RESIDUE_NAME, idx, Value::Str(info.residue_name.clone()))?;
            block.set(PDB_CHAIN_NAME, idx, Value::Str(info.chain_id.clone()))?;
            block.set(PDB_INSERTION_CODE, idx, Value::Str(pdb.insertion_code))?;
            block.set(PDB_RESIDUE_NUMBER, idx, Value::Int(info.residue_number))?;
            block.set(PDB_OCCUPANCY, idx, Value::Real(pdb.occupancy))?;
            block.set(PDB_TFACTOR, idx, Value::Real(pdb.temp_factor))?;
        }
        copy_properties(atom, prop_names, |name, value| block.set(&name, idx, value))?;
    }
    Ok(())
}

fn map_bonds(mol: &ROMol, prop_names: &[String]) -> Result<IndexedBlock> {
    let mut block = IndexedBlock::new(mol.num_bonds() as usize);
    let any_dative = mol
        .bonds()
        .iter()
        .any(|b| b.bond_type() == BondType::Dative);
    for bond in mol.bonds() {
        let idx = bond.idx();
        // atom indices are 1-based, and written in ascending order
        let mut from = bond.end_atom_idx() as i32 + 1;
        let mut to = bond.begin_atom_idx() as i32 + 1;
        if from > to {
            std::mem::swap(&mut from, &mut to);
        }
        block.set(BOND_ATOM_1, idx, Value::Int(from))?;
        block.set(BOND_ATOM_2, idx, Value::Int(to))?;
        block.set(BOND_ORDER, idx, Value::Int(bond_type_to_order(mol, idx)?))?;
        // an unspecified double bond stereo is the only one set; the others
        // follow from the coordinates
        if is_double_any_bond(mol, idx) {
            block.set(MAE_BOND_PARITY, idx, Value::Int(2))?;
        }
        if any_dative {
            block.set(
                MAE_BOND_DATIVE_MARK,
                idx,
                Value::Bool(bond.bond_type() == BondType::Dative),
            )?;
        }
        copy_properties(bond, prop_names, |name, value| block.set(&name, idx, value))?;
    }
    Ok(block)
}

/// `_MolToMaeCtBlock`
fn mol_to_ct_block(mol: &ROMol, conf_id: i32, prop_names: &[String]) -> Result<CtBlock> {
    if mol.num_atoms() == 0 {
        return Err(RDKitError::Value(
            "molecules without atoms cannot be exported to Maestro files.\n".into(),
        ));
    }
    let mut tmp = mol.clone();
    kekulize(&mut tmp, true, true, 100)?;
    if mol.num_conformers() == 0 {
        compute_2d_coords(&mut tmp, &Compute2DCoordParameters::default())?;
    }
    // computed properties are not exported
    RDProps::clear_computed_props(&mut tmp);
    for i in 0..tmp.num_atoms() {
        tmp.atom_mut(i).clear_computed_props();
    }
    for i in 0..tmp.num_bonds() {
        tmp.bond_mut(i).clear_computed_props();
    }
    let mut ct = CtBlock {
        atoms: IndexedBlock::new(tmp.num_atoms() as usize),
        ..Default::default()
    };
    map_mol_properties(&mut tmp, prop_names, &mut ct)?;
    map_atoms(&tmp, prop_names, conf_id, &mut ct)?;
    if tmp.num_bonds() > 0 {
        ct.bonds = Some(map_bonds(&tmp, prop_names)?);
    }
    Ok(ct)
}

// ------------------------------------------------------------------ output

/// A string value: quoted when empty or containing a space, a quote or a
/// backslash.
fn push_string(out: &mut String, s: &str) {
    if s.is_empty() || s.contains([' ', '"', '\\']) {
        out.push('"');
        for c in s.chars() {
            if c == '"' || c == '\\' {
                out.push('\\');
            }
            out.push(c);
        }
        out.push('"');
    } else {
        out.push_str(s);
    }
}

/// `%f`
fn push_real(out: &mut String, v: f64) {
    if v.is_nan() {
        out.push_str(if v.is_sign_negative() { "-nan" } else { "nan" });
    } else if v.is_infinite() {
        out.push_str(if v < 0.0 { "-inf" } else { "inf" });
    } else {
        let _ = write!(out, "{v:.6}");
    }
}

fn write_indexed_block(out: &mut String, name: &str, block: &IndexedBlock) {
    let _ = writeln!(out, "  {name}[{}] {{", block.size);
    out.push_str("    # First column is Index #\n");
    // the order of the columns: bool, real, int, string
    let mut cols: Vec<(&String, &Column)> = block.columns.iter().collect();
    cols.sort_by_key(|(n, c)| {
        let rank = match c {
            Column::Bool(_) => 0,
            Column::Real(_) => 1,
            Column::Int(_) => 2,
            Column::Str(_) => 3,
        };
        (rank, (*n).clone())
    });
    for (n, _) in &cols {
        let _ = writeln!(out, "    {n}");
    }
    out.push_str("    :::\n");
    for row in 0..block.size {
        let _ = write!(out, "    {}", row + 1);
        for (_, col) in &cols {
            out.push(' ');
            match col {
                Column::Bool(c) => match c[row] {
                    Some(v) => out.push(if v { '1' } else { '0' }),
                    None => out.push_str("<>"),
                },
                Column::Int(c) => match c[row] {
                    Some(v) => {
                        let _ = write!(out, "{v}");
                    }
                    None => out.push_str("<>"),
                },
                Column::Real(c) => match c[row] {
                    Some(v) => push_real(out, v),
                    None => out.push_str("<>"),
                },
                Column::Str(c) => match &c[row] {
                    Some(v) => push_string(out, v),
                    None => out.push_str("<>"),
                },
            }
        }
        out.push('\n');
    }
    out.push_str("    :::\n  }\n");
}

fn write_ct_block(ct: &CtBlock) -> String {
    let mut out = format!("{CT_BLOCK} {{\n");
    for n in ct.bools.keys() {
        let _ = writeln!(out, "  {n}");
    }
    for n in ct.reals.keys() {
        let _ = writeln!(out, "  {n}");
    }
    for n in ct.ints.keys() {
        let _ = writeln!(out, "  {n}");
    }
    for n in ct.strings.keys() {
        let _ = writeln!(out, "  {n}");
    }
    out.push_str("  :::\n");
    for v in ct.bools.values() {
        let _ = writeln!(out, "  {}", u8::from(*v));
    }
    for v in ct.reals.values() {
        out.push_str("  ");
        push_real(&mut out, *v);
        out.push('\n');
    }
    for v in ct.ints.values() {
        let _ = writeln!(out, "  {v}");
    }
    for v in ct.strings.values() {
        out.push_str("  ");
        push_string(&mut out, v);
        out.push('\n');
    }
    write_indexed_block(&mut out, ATOM_BLOCK, &ct.atoms);
    if let Some(bonds) = &ct.bonds {
        write_indexed_block(&mut out, BOND_BLOCK, bonds);
    }
    out.push_str("}\n\n");
    out
}

/// The text of the header block that opens a Maestro file.
const FILE_HEADER: &str = " {\n  s_m_m2io_version\n  :::\n  2.0.0\n}\n\n";

/// `MaeWriter::getText`: the text of the structure block of a molecule.
pub fn mol_to_mae_block(mol: &ROMol, conf_id: i32, prop_names: &[String]) -> Result<String> {
    Ok(write_ct_block(&mol_to_ct_block(mol, conf_id, prop_names)?))
}

/// The state of a `MaeWriter`: the text it produces is written by the caller.
#[derive(Debug, Default)]
pub struct MaeWriter {
    props: Vec<String>,
    mol_id: u32,
    header_written: bool,
}

impl MaeWriter {
    pub fn new() -> Self {
        Self::default()
    }

    /// `setProps`
    pub fn set_props(&mut self, prop_names: Vec<String>) {
        if self.mol_id > 0 {
            rd_log_warning!(
                "WARNING: Setting property list after a few molecules have been written\n"
            );
        }
        self.props = prop_names;
    }

    /// `write`: the text to append to the output (the first molecule comes
    /// with the header of the file).
    pub fn write(&mut self, mol: &ROMol, conf_id: i32) -> Result<String> {
        let block = mol_to_mae_block(mol, conf_id, &self.props)?;
        let text = if self.header_written {
            block
        } else {
            self.header_written = true;
            format!("{FILE_HEADER}{block}")
        };
        self.mol_id += 1;
        Ok(text)
    }

    /// `numMols`
    pub fn num_mols(&self) -> u32 {
        self.mol_id
    }
}
