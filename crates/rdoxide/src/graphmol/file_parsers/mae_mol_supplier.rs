//! Molecules from Maestro files.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MaeMolSupplier.cpp` and
//! `MaestroProperties.h` (2026.09.1); the reader is [`super::mae_reader`].
//!
//! Copyright (C) 2018 Pat Lorton and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::mae_reader::constants::*;
use super::mae_reader::{Block, Column, ParsedIndexedBlock, Reader};
use super::utils::to_int;
use crate::geometry::Point3D;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::chirality::from_coords::{
    assign_chiral_types_from_3d, detect_bond_stereochemistry,
};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_interchange::bond_type_from_order;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::{SANITIZE_ALL, SANITIZE_CLEANUP};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::monomer_info::AtomMonomerInfo;
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result, rd_log_warning};

pub const MAE_ATOM_TYPE: &str = "i_m_mmod_type";
pub const MAE_BOND_DATIVE_MARK: &str = "b_sPrivate_dative_bond";
pub const MAE_BOND_PARITY: &str = "i_sd_original_parity";
pub const MAE_ENHANCED_STEREO_STATUS: &str = "i_m_ct_enhanced_stereo_status";
pub const MAE_RGROUP_LABEL: &str = "i_sd__MolFileRLabel";
pub const MAE_STEREO_STATUS: &str = "i_m_ct_stereo_status";
pub const PDB_ATOM_NAME: &str = "s_m_pdb_atom_name";
pub const PDB_CHAIN_NAME: &str = "s_m_chain_name";
pub const PDB_INSERTION_CODE: &str = "s_m_insertion_code";
pub const PDB_OCCUPANCY: &str = "r_m_pdb_occupancy";
pub const PDB_RESIDUE_NAME: &str = "s_m_pdb_residue_name";
pub const PDB_RESIDUE_NUMBER: &str = "i_m_residue_number";
pub const PDB_TFACTOR: &str = "r_m_pdb_tfactor";

/// `MaeMolSupplierParams`
#[derive(Debug, Clone)]
pub struct MaeMolSupplierParams {
    pub sanitize: bool,
    pub remove_hs: bool,
}

impl Default for MaeMolSupplierParams {
    fn default() -> Self {
        MaeMolSupplierParams {
            sanitize: true,
            remove_hs: true,
        }
    }
}

/// `e.what()` of the exceptions caught in `MaeMolSupplier::next`.
fn what(e: &RDKitError) -> String {
    match e {
        RDKitError::Invariant { message, .. } => message.clone(),
        e => e.to_string(),
    }
}

fn file_parse(msg: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(msg.into())
}

/// `maeparser`'s `out_of_range` / `read_exception`, as `FileParseException`.
fn reader_error(msg: String) -> RDKitError {
    file_parse(msg)
}

fn column_at<T: Clone>(col: &Column<T>, i: usize) -> Option<T> {
    col.get(i).cloned().flatten()
}

/// `PDBInfo`: the PDB columns of the atom block.
struct PdbInfo<'a> {
    atom_name: Option<&'a Column<String>>,
    residue_name: Option<&'a Column<String>>,
    chain_id: Option<&'a Column<String>>,
    insertion_code: Option<&'a Column<String>>,
    resnum: Option<&'a Column<i32>>,
    occupancy: Option<&'a Column<f64>>,
    tempfac: Option<&'a Column<f64>>,
}

impl<'a> PdbInfo<'a> {
    fn new(b: &'a ParsedIndexedBlock) -> Self {
        PdbInfo {
            atom_name: b.string_property(PDB_ATOM_NAME),
            residue_name: b.string_property(PDB_RESIDUE_NAME),
            chain_id: b.string_property(PDB_CHAIN_NAME),
            insertion_code: b.string_property(PDB_INSERTION_CODE),
            resnum: b.ints.get(PDB_RESIDUE_NUMBER),
            occupancy: b.reals.get(PDB_OCCUPANCY),
            tempfac: b.reals.get(PDB_TFACTOR),
        }
    }

    fn add_pdb_data(&self, atom: &mut Atom, i: usize) {
        // a PDB atom name is needed to populate the info
        let Some(name) = self.atom_name.and_then(|c| column_at(c, i)) else {
            return;
        };
        let mut info = AtomMonomerInfo::new_pdb(&name);
        if let Some(v) = self.residue_name.and_then(|c| column_at(c, i)) {
            info.residue_name = v;
        }
        if let Some(v) = self.chain_id.and_then(|c| column_at(c, i)) {
            info.chain_id = v;
        }
        if let Some(v) = self.insertion_code.and_then(|c| column_at(c, i)) {
            info.pdb_mut().insertion_code = v;
        }
        if let Some(v) = self.resnum.and_then(|c| column_at(c, i)) {
            info.residue_number = v;
        }
        if let Some(v) = self.occupancy.and_then(|c| column_at(c, i)) {
            info.pdb_mut().occupancy = v;
        }
        if let Some(v) = self.tempfac.and_then(|c| column_at(c, i)) {
            info.pdb_mut().temp_factor = v;
        }
        atom.set_monomer_info(Some(info));
    }
}

/// `strip_prefix_from_mae_property`
fn strip_prefix_from_mae_property(prop_name: &str) -> String {
    let b = prop_name.as_bytes();
    if matches!(b.first(), Some(b'b' | b'i' | b'r' | b's')) {
        let rest = &prop_name[1..];
        if rest.starts_with("_rdk_") {
            return prop_name[6..].to_string();
        } else if rest.starts_with("_rdkit_") {
            return prop_name[8..].to_string();
        }
    }
    prop_name.to_string()
}

fn is_ignored_property(prop: &str) -> bool {
    prop == MAE_ENHANCED_STEREO_STATUS || prop == MAE_STEREO_STATUS
}

/// boost's `char_separator<char>("_")`: empty tokens are dropped.
fn tokens(s: &str) -> Vec<&str> {
    s.split('_').filter(|t| !t.is_empty()).collect()
}

fn to_int_token(t: &str) -> Result<i32> {
    to_int(t.as_bytes(), true)
}

/// `parseChiralityLabel`: whether the label is valid.
fn parse_chirality_label(mol: &mut RWMol, stereo_prop: &str) -> Result<bool> {
    let toks = tokens(stereo_prop);
    let Some(first) = toks.first() else {
        return Ok(false);
    };
    let chiral_idx = to_int_token(first)? - 1;
    if chiral_idx < 0 || chiral_idx >= mol.num_atoms() as i32 {
        return Ok(false);
    }
    let chiral_idx = chiral_idx as u32;
    let Some(rotation) = toks.get(1).and_then(|t| t.chars().last()) else {
        return Ok(false);
    };
    let mut n_swaps: i32 = match rotation {
        'R' => 0,
        'S' => 1,
        '?' => return Ok(true),
        _ => return Ok(false),
    };
    let mut bond_indexes = Vec::new();
    for t in &toks[2..] {
        let nbr = to_int_token(t)? - 1;
        let bond = if nbr >= 0 && (nbr as u32) < mol.num_atoms() {
            mol.bond_idx_between(chiral_idx, nbr as u32)
        } else {
            None
        };
        match bond {
            Some(b) => bond_indexes.push(b as i32),
            None => return Ok(false),
        }
    }
    if bond_indexes.len() != mol.atom_degree(chiral_idx) as usize {
        return Ok(false);
    }
    n_swaps += mol.perturbation_order(chiral_idx, &bond_indexes)?;
    mol.atom_mut(chiral_idx)
        .set_chiral_tag(if n_swaps % 2 == 0 {
            ChiralType::TetrahedralCW
        } else {
            ChiralType::TetrahedralCCW
        });
    Ok(true)
}

/// `parseStereoBondLabel`
fn parse_stereo_bond_label(mol: &mut RWMol, stereo_prop: &str) -> Result<()> {
    let mut stereo = BondStereo::StereoNone;
    let mut atoms = Vec::new();
    for t in tokens(stereo_prop) {
        match t {
            "E" => stereo = BondStereo::StereoTrans,
            "Z" => stereo = BondStereo::StereoCis,
            _ => atoms.push(to_int_token(t)? - 1),
        }
    }
    if stereo == BondStereo::StereoNone {
        return Err(RDKitError::invariant("bad prop value"));
    }
    // allenes and allene-likes are not supported
    if atoms.len() != 4 {
        return Ok(());
    }
    let in_range = |i: i32| i >= 0 && (i as u32) < mol.num_atoms();
    let bond = if in_range(atoms[1]) && in_range(atoms[2]) {
        mol.bond_idx_between(atoms[1] as u32, atoms[2] as u32)
    } else {
        None
    };
    let Some(bond) = bond else {
        return Err(RDKitError::invariant("bad stereo bond"));
    };
    if mol.bond(bond).bond_type() != BondType::Double {
        return Err(RDKitError::invariant("bad stereo bond"));
    }
    let b = mol.bond_mut(bond);
    *b.stereo_atoms_mut() = vec![atoms[0], atoms[3]];
    b.set_stereo(stereo)
}

/// `set_mol_properties`: the structure properties, with the stereochemistry.
fn set_mol_properties(mol: &mut RWMol, ct: &Block) -> Result<()> {
    for (name, value) in &ct.strings {
        if is_ignored_property(name) {
            continue;
        }
        if name == CT_TITLE {
            mol.set_prop(cp::_Name, value.clone(), false)?;
        } else if name.starts_with(CT_CHIRALITY_PROP_PREFIX)
            || name.starts_with(CT_PSEUDOCHIRALITY_PROP_PREFIX)
        {
            // an invalid label is only reported
            if !parse_chirality_label(mol, value)? {
                rd_log_warning!("Ignoring invalid chirality label: '{value}'\n");
            }
        } else if name.starts_with(CT_EZ_PROP_PREFIX) {
            parse_stereo_bond_label(mol, value)?;
        } else {
            mol.set_prop(&strip_prefix_from_mae_property(name), value.clone(), false)?;
        }
    }
    for (name, v) in &ct.reals {
        if !is_ignored_property(name) {
            mol.set_prop(&strip_prefix_from_mae_property(name), *v, false)?;
        }
    }
    for (name, v) in &ct.ints {
        if !is_ignored_property(name) {
            mol.set_prop(&strip_prefix_from_mae_property(name), *v, false)?;
        }
    }
    for (name, v) in &ct.bools {
        if !is_ignored_property(name) {
            mol.set_prop(&strip_prefix_from_mae_property(name), *v, false)?;
        }
    }
    Ok(())
}

/// `set_atom_properties`: the properties that were not used to build the
/// atom; undefined values are skipped.
fn set_atom_properties(atom: &mut Atom, block: &ParsedIndexedBlock, i: usize) -> Result<()> {
    for (name, col) in &block.strings {
        if [
            PDB_ATOM_NAME,
            PDB_RESIDUE_NAME,
            PDB_CHAIN_NAME,
            PDB_INSERTION_CODE,
        ]
        .contains(&name.as_str())
        {
            continue;
        }
        if let Some(v) = column_at(col, i) {
            atom.set_prop(&strip_prefix_from_mae_property(name), v, false)?;
        }
    }
    for (name, col) in &block.reals {
        if [
            ATOM_X_COORD,
            ATOM_Y_COORD,
            ATOM_Z_COORD,
            PDB_OCCUPANCY,
            PDB_TFACTOR,
        ]
        .contains(&name.as_str())
        {
            continue;
        }
        if let Some(v) = column_at(col, i) {
            atom.set_prop(&strip_prefix_from_mae_property(name), v, false)?;
        }
    }
    for (name, col) in &block.ints {
        if name == ATOM_ATOMIC_NUM || name == PDB_RESIDUE_NUMBER {
            continue;
        }
        let Some(v) = column_at(col, i) else { continue };
        if name == ATOM_FORMAL_CHARGE {
            atom.set_formal_charge(v);
        } else if name == MAE_RGROUP_LABEL {
            atom.set_prop(cp::_MolFileRLabel, v, false)?;
        } else {
            atom.set_prop(&strip_prefix_from_mae_property(name), v, false)?;
        }
    }
    for (name, col) in &block.bools {
        if let Some(v) = column_at(col, i) {
            atom.set_prop(&strip_prefix_from_mae_property(name), v, false)?;
        }
    }
    Ok(())
}

fn required<T>(res: std::result::Result<&T, String>) -> Result<&T> {
    res.map_err(file_parse)
}

/// `addAtoms`; the atoms with reserved atomic numbers are returned for removal.
fn add_atoms(block: &ParsedIndexedBlock, mol: &mut RWMol) -> Result<Vec<u32>> {
    let numbers = required(block.int_property(ATOM_ATOMIC_NUM))?;
    let xs = required(block.real_property(ATOM_X_COORD))?;
    let ys = required(block.real_property(ATOM_Y_COORD))?;
    let zs = required(block.real_property(ATOM_Z_COORD))?;
    let size = numbers.len();
    let mut conf = Conformer::new(size);
    conf.set_id(0);
    let pdb = PdbInfo::new(block);
    let undefined = || file_parse("Indexed property value undefined.");
    let mut to_remove = Vec::new();
    let mut nonzero_z = false;
    for i in 0..size {
        let mut remove = false;
        let mut z = numbers[i].ok_or_else(undefined)?;
        if z == 0 || z == -1 || z == -3 {
            rd_log_warning!(
                "WARNING: atom {} in input Maestro file has atomic number '{z}', which is reserved for internal use, and not allowed in inputs. The atom will be ignored.",
                i + 1
            );
            // removing the atom now would disturb the bond parsing
            to_remove.push(i as u32);
            remove = true;
            z = 0;
        } else if z == -2 {
            // dummy atom
            z = 0;
        }
        if z < 0
            || z as u32
                > crate::graphmol::periodic_table::PeriodicTable::get().get_max_atomic_number()
        {
            return Err(file_parse("Atomic number not found"));
        }
        let idx = mol.add_atom(Atom::new(z as u32), true);
        let pos = Point3D::new(
            xs.get(i).copied().flatten().ok_or_else(undefined)?,
            ys.get(i).copied().flatten().ok_or_else(undefined)?,
            zs.get(i).copied().flatten().ok_or_else(undefined)?,
        );
        conf.set_atom_pos(idx, pos)?;
        // an atom that will be removed needs neither PDB data nor properties
        if !remove {
            let atom = mol.atom_mut(idx);
            pdb.add_pdb_data(atom, i);
            set_atom_properties(atom, block, i)?;
            nonzero_z |= pos.z.abs() > 1.0e-4;
        }
    }
    conf.set_3d(nonzero_z);
    mol.add_conformer(conf, false)?;
    Ok(to_remove)
}

/// `addBonds`
fn add_bonds(block: &ParsedIndexedBlock, mol: &mut RWMol) -> Result<()> {
    let from = required(block.int_property(BOND_ATOM_1))?;
    let to = required(block.int_property(BOND_ATOM_2))?;
    let orders = required(block.int_property(BOND_ORDER))?;
    let undefined = || file_parse("Indexed property value undefined.");
    for i in 0..from.len() {
        // Maestro atoms are 1-based
        let a1 = from[i].ok_or_else(undefined)? - 1;
        let a2 = to[i].ok_or_else(undefined)? - 1;
        let order = bond_type_from_order(i64::from(orders[i].ok_or_else(undefined)?))
            .ok_or_else(|| file_parse("bad bond order"))?;
        let (n1, n2) = (a1 as u32, a2 as u32);
        if a1 < 0 || a2 < 0 || n1 >= mol.num_atoms() || n2 >= mol.num_atoms() {
            return Err(RDKitError::range_at(
                "idx",
                if n1 >= mol.num_atoms() { n1 } else { n2 },
                mol.num_atoms(),
                "Code/GraphMol/ROMol.cpp",
                202,
            ));
        }
        if let Some(existing) = mol.bond_idx_between(n1, n2) {
            // Maestro files may double-list some bonds
            if order != mol.bond(existing).bond_type() {
                rd_log_warning!(
                    "WARNING: bond between atoms {a1} and {a2} is defined more than once with different bond orders. The first definition will be honored, and the rest will be ignored."
                );
            }
            continue;
        }
        mol.add_bond(n1, n2, order)?;
    }
    Ok(())
}

/// `build_mol`
fn build_mol(mol: &mut RWMol, ct: &Block, sanitize: bool, do_remove_hs: bool) -> Result<()> {
    let atom_block = ct.indexed_block(ATOM_BLOCK).map_err(file_parse)?;
    let atoms = atom_block.parse().map_err(file_parse)?;
    let to_remove = add_atoms(&atoms, mol)?;
    // the atom block is mandatory, the bond block is not
    if let Ok(bond_block) = ct.indexed_block(BOND_BLOCK) {
        let bonds = bond_block.parse().map_err(file_parse)?;
        add_bonds(&bonds, mol)?;
    }
    // these need the atoms and bonds
    set_mol_properties(mol, ct)?;

    // `replaceExistingTags` is false; as upstream passes it where
    // `detectBondStereochemistry` expects a conformer id, that is conformer 0
    if sanitize {
        if do_remove_hs {
            // bond stereo detection must happen before the removal of Hs,
            // but after the cleanup of e.g. nitro groups
            sanitize_mol(mol, SANITIZE_CLEANUP).map_err(|(_, e)| e)?;
            detect_bond_stereochemistry(mol, 0)?;
            remove_hs(mol, &RemoveHsParameters::default(), true)?;
        } else {
            sanitize_mol(mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
            detect_bond_stereochemistry(mol, 0)?;
        }
    } else {
        // some properties are needed for the chirality
        mol.update_property_cache(false)?;
        detect_bond_stereochemistry(mol, 0)?;
    }

    // chiralities from 3D coordinates do not override the ones of the properties
    if mol.num_conformers() > 0 && mol.conformer(-1)?.is_3d() {
        assign_chiral_types_from_3d(mol, -1, false)?;
    }
    assign_stereochemistry(mol, false, false, false)?;

    if !to_remove.is_empty() {
        mol.begin_batch_edit()?;
        for idx in to_remove {
            mol.remove_atom(idx)?;
        }
        mol.commit_batch_edit()?;
    }
    Ok(())
}

fn index_error(idx: u32) -> RDKitError {
    file_parse(format!(
        "ERROR: Index error (idx = {idx}) :  we do no have enough ct blocks"
    ))
}

/// `v2::FileParsers::MaeMolSupplier`
#[derive(Debug)]
pub struct MaeMolSupplier {
    reader: Reader,
    params: MaeMolSupplierParams,
    next_struct: Option<Block>,
    stored_exc: String,
    position: u32,
    length: u32,
}

impl MaeMolSupplier {
    pub fn from_bytes(data: Vec<u8>, params: MaeMolSupplierParams) -> Result<Self> {
        let mut s = MaeMolSupplier {
            reader: Reader::new(data),
            params,
            next_struct: None,
            stored_exc: String::new(),
            position: 0,
            length: 0,
        };
        s.init()?;
        Ok(s)
    }

    pub fn from_file(path: &str, params: MaeMolSupplierParams) -> Result<Self> {
        let data = std::fs::read(path)
            .map_err(|_| RDKitError::BadFile(format!("Bad input file {path}")))?;
        Self::from_bytes(data, params)
    }

    fn init(&mut self) -> Result<()> {
        self.position = 0;
        self.length = 0;
        self.next_struct = self.reader.next(CT_BLOCK).map_err(reader_error)?;
        Ok(())
    }

    /// `reset`
    pub fn reset(&mut self) -> Result<()> {
        let length = self.length;
        self.reader.rewind();
        self.init()?;
        self.length = length;
        Ok(())
    }

    /// `MaeMolSupplier::next`
    pub fn next_mol(&mut self) -> Result<RWMol> {
        if !self.stored_exc.is_empty() {
            return Err(file_parse(self.stored_exc.clone()));
        } else if self.at_end() {
            return Err(file_parse("All structures read from Maestro file"));
        }
        let mut mol = RWMol::new();
        let ct = self.next_struct.take().expect("not at the end");
        let res = build_mol(&mut mol, &ct, self.params.sanitize, self.params.remove_hs);
        self.next_struct = Some(ct);
        self.move_to_next_block();
        match res {
            Ok(()) => Ok(mol),
            Err(e) => Err(file_parse(what(&e))),
        }
    }

    fn move_to_next_block(&mut self) {
        match self.reader.next(CT_BLOCK) {
            Ok(next) => self.next_struct = next,
            // the structure stays as it was
            Err(e) => self.stored_exc = e,
        }
        self.position += 1;
    }

    /// `atEnd`
    pub fn at_end(&mut self) -> bool {
        if self.next_struct.is_none() {
            self.length = self.position;
            return true;
        }
        false
    }

    /// `length`: scans the whole data once.
    pub fn length(&mut self) -> Result<u32> {
        if self.length == 0 && !self.at_end() {
            let mut reader = self.reader.clone();
            reader.rewind();
            let mut count = 0u32;
            loop {
                match reader.next(CT_BLOCK) {
                    Ok(Some(_)) => count += 1,
                    Ok(None) => break,
                    // the first structure must be readable
                    Err(e) if count == 0 => return Err(file_parse(e)),
                    Err(_) => break,
                }
            }
            self.length = count;
        }
        Ok(self.length)
    }

    /// `moveTo`
    pub fn move_to(&mut self, idx: u32) -> Result<()> {
        if self.length > 0 && idx > self.length {
            return Err(index_error(idx));
        }
        if idx < self.position {
            self.reset()?;
        }
        while idx > self.position {
            self.move_to_next_block();
            if self.at_end() {
                return Err(index_error(idx));
            }
        }
        Ok(())
    }

    /// `operator[]`
    pub fn get(&mut self, idx: u32) -> Result<RWMol> {
        self.move_to(idx)?;
        self.next_mol()
    }
}
