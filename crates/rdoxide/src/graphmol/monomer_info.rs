//! Ported from RDKit `Code/GraphMol/MonomerInfo.{h,cpp}`.
//!
//! Copyright (C) 2013-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Upstream models this as a small class hierarchy (`AtomMonomerInfo` and the
//! derived `AtomPDBResidueInfo`); here the derived class is an optional
//! extension of the base data so that [`AtomMonomerInfo`] stays a plain value.

/// `AtomMonomerInfo::AtomMonomerType`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum AtomMonomerType {
    #[default]
    Unknown = 0,
    PdbResidue,
    Other,
}

impl AtomMonomerType {
    pub const ALL: [AtomMonomerType; 3] = [Self::Unknown, Self::PdbResidue, Self::Other];
    /// Name of the Python enum value (`rdchem.AtomMonomerType.X`).
    pub fn name(self) -> &'static str {
        match self {
            Self::Unknown => "UNKNOWN",
            Self::PdbResidue => "PDBRESIDUE",
            Self::Other => "OTHER",
        }
    }
    pub fn from_u8(v: u8) -> Option<Self> {
        Self::ALL.get(v as usize).copied()
    }
}

/// The PDB-specific fields of `AtomPDBResidueInfo`.
#[derive(Debug, Clone, PartialEq)]
pub struct PdbResidueData {
    pub serial_number: u32,
    pub alt_loc: String,
    pub insertion_code: String,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub is_hetero_atom: bool,
    pub secondary_structure: u32,
    pub segment_number: u32,
}

impl Default for PdbResidueData {
    fn default() -> Self {
        PdbResidueData {
            serial_number: 0,
            alt_loc: String::new(),
            insertion_code: String::new(),
            occupancy: 1.0,
            temp_factor: 0.0,
            is_hetero_atom: false,
            secondary_structure: 0,
            segment_number: 0,
        }
    }
}

/// `AtomMonomerInfo` (and `AtomPDBResidueInfo` when `pdb` is set).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AtomMonomerInfo {
    pub monomer_type: AtomMonomerType,
    pub name: String,
    pub residue_name: String,
    pub residue_number: i32,
    pub chain_id: String,
    pub monomer_class: String,
    pub pdb: Option<PdbResidueData>,
}

impl AtomMonomerInfo {
    /// `AtomMonomerInfo(typ, nm)`
    pub fn new(monomer_type: AtomMonomerType, name: &str) -> Self {
        AtomMonomerInfo {
            monomer_type,
            name: name.to_string(),
            ..Default::default()
        }
    }

    /// `AtomPDBResidueInfo(atomName)`: all PDB fields at their defaults.
    pub fn new_pdb(name: &str) -> Self {
        AtomMonomerInfo {
            monomer_type: AtomMonomerType::PdbResidue,
            name: name.to_string(),
            pdb: Some(PdbResidueData::default()),
            ..Default::default()
        }
    }

    /// Whether this is an `AtomPDBResidueInfo`.
    pub fn is_pdb_residue(&self) -> bool {
        self.pdb.is_some()
    }

    /// The PDB fields, creating defaults if needed (callers have checked
    /// [`Self::is_pdb_residue`]).
    pub fn pdb_mut(&mut self) -> &mut PdbResidueData {
        self.pdb.get_or_insert_with(PdbResidueData::default)
    }
}
