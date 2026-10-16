//! The structure checker (`RDKit::StructureCheck`): checks and fixes the
//! structure of molecules according to a set of "augmented atom" rules.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker` (2026.09.1). This is
//! beta quality code upstream; its quirks are kept so the results match.
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

mod aa_data;
pub mod atom_symbol;
pub mod options;
pub mod pattern;
pub mod recharge;
pub mod stereo;
pub mod strip_small_fragments;
pub mod tautomer;
pub mod utilities;

pub use options::{
    StructCheckerOptions, load_options_from_files, parse_options_json, string_to_augmented_atom,
};

use crate::Result;
use crate::graphmol::ring_info::FindRingType;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};

/// Charge that matches any charge.
pub const ANY_CHARGE: i32 = 8;

/// `RadicalType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RadicalType {
    None = 0,
    Singlet = 1,
    Doublet = 2,
    Triplet = 3,
    Any = 0xFF,
}

/// `AABondType`: MDL CTFile bond types plus extensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AABondType {
    /// Means: remove the bond.
    None = 0,
    Single = 1,
    Double = 2,
    Triple = 3,
    Aromatic = 4,
    SingleDouble = 5,
    SingleAromatic = 6,
    DoubleAromatic = 7,
    Any = 8,
    AllBondTypes = 0xF,
}

/// `AATopology`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AATopology {
    /// Don't care.
    None = 0,
    Ring = 1,
    Chain = 2,
}

/// A neighbour of an augmented atom.
#[derive(Debug, Clone, PartialEq)]
pub struct Ligand {
    pub atom_symbol: String,
    pub charge: i32,
    pub radical: RadicalType,
    /// Substitution count, 0 = don't care.
    pub substitution_count: u32,
    pub bond_type: AABondType,
}

impl Default for Ligand {
    fn default() -> Self {
        Ligand {
            atom_symbol: String::new(),
            charge: ANY_CHARGE,
            radical: RadicalType::Any,
            substitution_count: 0,
            bond_type: AABondType::Any,
        }
    }
}

/// An atom with its neighbours.
#[derive(Debug, Clone, PartialEq)]
pub struct AugmentedAtom {
    pub atom_symbol: String,
    pub short_name: String,
    pub charge: i32,
    pub radical: RadicalType,
    pub topology: AATopology,
    pub ligands: Vec<Ligand>,
}

impl Default for AugmentedAtom {
    fn default() -> Self {
        AugmentedAtom {
            atom_symbol: String::new(),
            short_name: String::new(),
            charge: ANY_CHARGE,
            radical: RadicalType::Any,
            topology: AATopology::None,
            ligands: Vec::new(),
        }
    }
}

/// `IncEntry`: the (unused upstream: the loader is missing) charge and
/// acidity increments of an atom type.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct IncEntry {
    pub atom_symbol: String,
    pub local_inc: f64,
    pub alpha_inc: f64,
    pub beta_inc: f64,
    pub mult_inc: f64,
}

/// `PathEntry`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PathEntry {
    pub path: AugmentedAtom,
    pub cond: f64,
}

/// `StructChecker::StructureFlags`
pub mod structure_flags {
    pub const NO_CHANGE: u32 = 0;
    pub const BAD_MOLECULE: u32 = 0x0001;
    pub const ALIAS_CONVERSION_FAILED: u32 = 0x0002;
    pub const STEREO_ERROR: u32 = 0x0004;
    pub const STEREO_FORCED_BAD: u32 = 0x0008;
    pub const ATOM_CLASH: u32 = 0x0010;
    pub const ATOM_CHECK_FAILED: u32 = 0x0020;
    pub const SIZE_CHECK_FAILED: u32 = 0x0040;
    pub const TRANSFORMED: u32 = 0x0100;
    pub const FRAGMENTS_FOUND: u32 = 0x0200;
    pub const EITHER_WARNING: u32 = 0x0400;
    pub const DUBIOUS_STEREO_REMOVED: u32 = 0x0800;
    pub const RECHARGED: u32 = 0x1000;
    pub const STEREO_TRANSFORMED: u32 = 0x2000;
    pub const TEMPLATE_TRANSFORMED: u32 = 0x4000;
    pub const TAUTOMER_TRANSFORMED: u32 = 0x8000;
    pub const BAD_SET: u32 = BAD_MOLECULE
        | ALIAS_CONVERSION_FAILED
        | STEREO_ERROR
        | STEREO_FORCED_BAD
        | ATOM_CLASH
        | ATOM_CHECK_FAILED
        | SIZE_CHECK_FAILED;
    pub const TRANSFORMED_SET: u32 = TRANSFORMED
        | FRAGMENTS_FOUND
        | EITHER_WARNING
        | DUBIOUS_STEREO_REMOVED
        | STEREO_TRANSFORMED
        | TEMPLATE_TRANSFORMED
        | TAUTOMER_TRANSFORMED
        | RECHARGED;
}
use structure_flags::*;

/// The names of the bits of the flags (the empty one is reserved).
const FLAG_NAMES: [&str; 16] = [
    "BAD_MOLECULE",
    "ALIAS_CONVERSION_FAILED",
    "STEREO_ERROR",
    "STEREO_FORCED_BAD",
    "ATOM_CLASH",
    "ATOM_CHECK_FAILED",
    "SIZE_CHECK_FAILED",
    "",
    "TRANSFORMED",
    "FRAGMENTS_FOUND",
    "EITHER_WARNING",
    "DUBIOUS_STEREO_REMOVED",
    "RECHARGED",
    "STEREO_TRANSFORMED",
    "TEMPLATE_TRANSFORMED",
    "TAUTOMER_TRANSFORMED",
];

/// Checks and fixes the structure of molecules.
#[derive(Debug, Clone, Default)]
pub struct StructChecker {
    options: StructCheckerOptions,
}

impl StructChecker {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_options(options: StructCheckerOptions) -> Self {
        StructChecker { options }
    }

    pub fn options(&self) -> &StructCheckerOptions {
        &self.options
    }

    pub fn set_options(&mut self, options: StructCheckerOptions) {
        self.options = options;
    }

    /// `StructureFlagsToString`: the names of the set flags, comma separated.
    pub fn structure_flags_to_string(flags: u32) -> String {
        let mut names = Vec::new();
        for (bit, name) in FLAG_NAMES.iter().enumerate() {
            if flags & (1 << bit) != 0 {
                names.push(*name);
            }
        }
        names.join(",")
    }

    /// `StringToStructureFlags`: unknown names are ignored.
    pub fn string_to_structure_flags(s: &str) -> u32 {
        let mut f = 0;
        for token in s.split(|c: char| c == ',' || (c as u32) <= (b' ' as u32)) {
            if let Some(bit) = FLAG_NAMES.iter().position(|n| !n.is_empty() && *n == token) {
                f |= 1 << bit;
            }
        }
        f
    }

    /// `checkMolStructure`: checks the molecule and fixes it if needed;
    /// returns the set of flags that describes what has been done.
    pub fn check_mol_structure(&self, mol: &mut ROMol) -> Result<u32> {
        let opts = &self.options;
        let mut flags = NO_CHANGE;

        if opts.max_mol_size != 0
            && (mol.num_atoms() > opts.max_mol_size || mol.num_bonds() > opts.max_mol_size)
        {
            return Ok(SIZE_CHECK_FAILED);
        }
        if mol.num_atoms() == 0 {
            return Ok(SIZE_CHECK_FAILED);
        }
        if !mol.ring_info().is_initialized() {
            mol.ring_info_mut().initialize(FindRingType::OtherOrUnknown);
        }

        if !opts.augmented_atom_pairs.is_empty()
            && pattern::transform_augmented_atoms(mol, &opts.augmented_atom_pairs, opts.verbose)?
        {
            flags |= TRANSFORMED;
            if opts.verbose {
                log_smiles(mol, "")?;
            }
        }

        let mut stereo_result = stereo::dubious_stereochemistry(mol);
        if stereo::fix_dubious_3d_molecule(mol)? & stereo::CONVERTED_TO_2D != 0 {
            stereo_result = 1;
            flags |= DUBIOUS_STEREO_REMOVED;
            if opts.verbose {
                log_smiles(mol, "")?;
            }
        }

        if opts.remove_minor_fragments {
            strip_small_fragments::add_mw_mf(mol, true)?;
            if strip_small_fragments::strip_small_fragments(mol, opts.verbose)? {
                flags |= FRAGMENTS_FOUND;
                if opts.verbose {
                    log_smiles(mol, "Striped SmallFragments")?;
                }
            }
            strip_small_fragments::add_mw_mf(mol, false)?;
        }

        // tautomer standardization
        for i in 0..opts.from_tautomer.len() {
            if opts.verbose {
                crate::rd_log_info!("tautomerizing with rule {i}\n");
            }
            // limit to 3 runs per rule
            for _ in 0..3 {
                let sct = tautomer::StructCheckTautomer::new(opts);
                if !sct.apply_tautomer(mol, i)? {
                    break;
                }
                flags |= TAUTOMER_TRANSFORMED;
                if opts.verbose {
                    let smi = mol_to_smiles(mol, &SmilesWriteParams::default())?;
                    crate::rd_log_info!("molecule: has been tautomerized with rule {i}:\n{smi}\n");
                }
            }
        }

        if stereo_result == stereo::EITHER_BOND_FOUND {
            flags |= EITHER_WARNING;
            stereo::remove_dubious_stereochemistry(mol);
            flags |= DUBIOUS_STEREO_REMOVED;
        } else if stereo_result > stereo::EITHER_BOND_FOUND {
            flags |= STEREO_ERROR;
            if opts.check_stereo {
                flags |= BAD_MOLECULE;
            } else {
                stereo::remove_dubious_stereochemistry(mol);
                flags |= DUBIOUS_STEREO_REMOVED;
            }
        }

        if recharge::total_charge(mol) != opts.desired_charge {
            let mut ch = recharge::ChargeFix::new(opts, mol);
            if ch.recharge_molecule()?.0 {
                flags |= RECHARGED;
                if opts.verbose {
                    log_smiles(mol, "Recharged Molecule:")?;
                }
            }
        }

        let clash_limit = f64::from(opts.collision_limit_percent) / 100.0;
        if opts.check_collisions && stereo::atom_clash(mol, clash_limit)? {
            flags |= ATOM_CLASH;
            if opts.verbose {
                log_smiles(mol, "AtomClash done:")?;
            }
        }

        if !opts.good_atoms.is_empty() && !pattern::check_atoms(mol, &opts.good_atoms, opts.verbose)
        {
            flags |= ATOM_CHECK_FAILED;
        }

        if opts.check_stereo && !stereo::check_stereo(mol)? {
            flags |= STEREO_ERROR;
        }

        // the end: "sanitize" the molecule
        if flags & TRANSFORMED != 0 {
            if mol.ring_info().is_initialized() {
                mol.ring_info_mut().reset(true);
            }
            mol.ring_info_mut().initialize(FindRingType::OtherOrUnknown);
        }
        Ok(flags)
    }
}

fn log_smiles(mol: &mut ROMol, prefix: &str) -> Result<()> {
    let smi = mol_to_smiles(mol, &SmilesWriteParams::default())?;
    crate::rd_log_info!("{prefix}{smi}\n");
    Ok(())
}
