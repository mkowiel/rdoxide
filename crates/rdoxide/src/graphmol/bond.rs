//! Ported from RDKit `Code/GraphMol/Bond.{h,cpp}`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::romol::ROMol;
use crate::rdgeneral::{Dict, RDProps};
use crate::{RDKitError, Result};

macro_rules! rd_enum {
    ($(#[$m:meta])* $name:ident { $($variant:ident = $py:literal),+ $(,)? }) => {
        $(#[$m])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
        #[repr(u8)]
        pub enum $name {
            #[default]
            $($variant),+
        }
        impl $name {
            pub const ALL: &'static [$name] = &[$($name::$variant),+];
            /// Name of the value in the Python enum.
            pub fn name(self) -> &'static str {
                match self { $($name::$variant => $py),+ }
            }
            pub fn from_u8(v: u8) -> Option<Self> {
                Self::ALL.get(v as usize).copied()
            }
        }
    };
}

rd_enum!(
    /// `Bond::BondType`
    BondType {
        Unspecified = "UNSPECIFIED",
        Single = "SINGLE",
        Double = "DOUBLE",
        Triple = "TRIPLE",
        Quadruple = "QUADRUPLE",
        Quintuple = "QUINTUPLE",
        Hextuple = "HEXTUPLE",
        OneAndAHalf = "ONEANDAHALF",
        TwoAndAHalf = "TWOANDAHALF",
        ThreeAndAHalf = "THREEANDAHALF",
        FourAndAHalf = "FOURANDAHALF",
        FiveAndAHalf = "FIVEANDAHALF",
        Aromatic = "AROMATIC",
        Ionic = "IONIC",
        Hydrogen = "HYDROGEN",
        ThreeCenter = "THREECENTER",
        DativeOne = "DATIVEONE",
        Dative = "DATIVE",
        DativeL = "DATIVEL",
        DativeR = "DATIVER",
        Other = "OTHER",
        Zero = "ZERO",
    }
);

rd_enum!(
    /// `Bond::BondDir`
    BondDir {
        None = "NONE",
        BeginWedge = "BEGINWEDGE",
        BeginDash = "BEGINDASH",
        EndDownRight = "ENDDOWNRIGHT",
        EndUpRight = "ENDUPRIGHT",
        EitherDouble = "EITHERDOUBLE",
        Unknown = "UNKNOWN",
    }
);

rd_enum!(
    /// `Bond::BondStereo`
    BondStereo {
        StereoNone = "STEREONONE",
        StereoAny = "STEREOANY",
        StereoZ = "STEREOZ",
        StereoE = "STEREOE",
        StereoCis = "STEREOCIS",
        StereoTrans = "STEREOTRANS",
        StereoAtropCW = "STEREOATROPCW",
        StereoAtropCCW = "STEREOATROPCCW",
    }
);

impl BondType {
    /// `getBondTypeAsDouble` (upstream raises for types without a value).
    pub fn as_double(self) -> Result<f64> {
        use BondType::*;
        Ok(match self {
            Unspecified | Ionic | Zero | Hydrogen => 0.0,
            Single | DativeOne | Dative => 1.0,
            Double => 2.0,
            Triple => 3.0,
            Quadruple => 4.0,
            Quintuple => 5.0,
            Hextuple => 6.0,
            OneAndAHalf | Aromatic => 1.5,
            TwoAndAHalf => 2.5,
            ThreeAndAHalf => 3.5,
            FourAndAHalf => 4.5,
            FiveAndAHalf => 5.5,
            ThreeCenter | DativeL | DativeR | Other => {
                return Err(RDKitError::invariant("Bad bond type"));
            }
        })
    }

    /// `getTwiceBondType`
    pub fn twice(self) -> Result<u8> {
        Ok((self.as_double()? * 2.0) as u8)
    }
}

/// A bond (`RDKit::Bond`).
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub(crate) bond_type: BondType,
    pub(crate) dir: BondDir,
    pub(crate) stereo: BondStereo,
    pub(crate) begin: u32,
    pub(crate) end: u32,
    pub(crate) idx: u32,
    pub(crate) is_aromatic: bool,
    pub(crate) is_conjugated: bool,
    pub(crate) stereo_atoms: Vec<i32>,
    pub(crate) flags: u64,
    pub(crate) has_owner: bool,
    pub(crate) props: Dict,
    /// `QueryBond::dp_query`
    pub(crate) query: Option<Box<crate::graphmol::queries::BondQuery>>,
}

impl Default for Bond {
    fn default() -> Self {
        Bond::new(BondType::Unspecified)
    }
}

impl RDProps for Bond {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

impl Bond {
    pub fn new(bond_type: BondType) -> Self {
        Bond {
            bond_type,
            dir: BondDir::None,
            stereo: BondStereo::StereoNone,
            begin: 0,
            end: 0,
            idx: 0,
            is_aromatic: false,
            is_conjugated: false,
            stereo_atoms: Vec::new(),
            flags: 0,
            has_owner: false,
            props: Dict::new(),
            query: None,
        }
    }

    /// Copy without ownership (`Bond::copy`); keeps indices like upstream.
    pub fn copy_detached(&self) -> Bond {
        let mut b = self.clone();
        b.has_owner = false;
        b
    }

    pub fn bond_type(&self) -> BondType {
        self.bond_type
    }
    pub fn set_bond_type(&mut self, t: BondType) {
        self.bond_type = t;
    }
    pub fn bond_dir(&self) -> BondDir {
        self.dir
    }
    pub fn set_bond_dir(&mut self, d: BondDir) {
        self.dir = d;
    }
    pub fn stereo(&self) -> BondStereo {
        self.stereo
    }
    /// `setStereo`: cis/trans need stereo atoms.
    pub fn set_stereo(&mut self, s: BondStereo) -> Result<()> {
        if matches!(s, BondStereo::StereoCis | BondStereo::StereoTrans)
            && self.stereo_atoms.len() != 2
        {
            return Err(RDKitError::precondition(
                "Stereo atoms should be specified before specifying CIS/TRANS bond stereochemistry",
            ));
        }
        self.stereo = s;
        Ok(())
    }
    pub fn begin_atom_idx(&self) -> u32 {
        self.begin
    }
    pub fn end_atom_idx(&self) -> u32 {
        self.end
    }
    pub fn idx(&self) -> u32 {
        self.idx
    }
    pub fn set_idx(&mut self, i: u32) {
        self.idx = i;
    }
    pub fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }
    pub fn set_is_aromatic(&mut self, v: bool) {
        self.is_aromatic = v;
    }
    pub fn is_conjugated(&self) -> bool {
        self.is_conjugated
    }
    pub fn set_is_conjugated(&mut self, v: bool) {
        self.is_conjugated = v;
    }
    pub fn stereo_atoms(&self) -> &[i32] {
        &self.stereo_atoms
    }
    pub fn stereo_atoms_mut(&mut self) -> &mut Vec<i32> {
        &mut self.stereo_atoms
    }
    pub fn flags(&self) -> u64 {
        self.flags
    }
    pub fn set_flags(&mut self, f: u64) {
        self.flags = f;
    }
    pub fn has_owning_mol(&self) -> bool {
        self.has_owner
    }
    pub fn has_query(&self) -> bool {
        self.query.is_some()
    }
    /// `getQuery()`
    pub fn query(&self) -> Option<&crate::graphmol::queries::BondQuery> {
        self.query.as_deref()
    }
    /// `setQuery()`
    pub fn set_query(&mut self, q: Option<crate::graphmol::queries::BondQuery>) {
        self.query = q.map(Box::new);
    }
    /// `QueryOps::hasComplexBondTypeQuery`
    pub fn has_complex_bond_type_query(&self) -> bool {
        self.query
            .as_deref()
            .is_some_and(crate::graphmol::queries::has_complex_bond_type_query)
    }

    /// For detached bonds only (owned bonds go through `ROMol`).
    pub fn set_begin_atom_idx_unchecked(&mut self, i: u32) {
        self.begin = i;
    }
    pub fn set_end_atom_idx_unchecked(&mut self, i: u32) {
        self.end = i;
    }

    /// `getOtherAtomIdx`
    pub fn other_atom_idx(&self, this_idx: u32) -> Result<u32> {
        if self.begin == this_idx {
            Ok(self.end)
        } else if self.end == this_idx {
            Ok(self.begin)
        } else {
            Err(RDKitError::postcondition("bad index"))
        }
    }

    /// `getBondTypeAsDouble`
    pub fn bond_type_as_double(&self) -> Result<f64> {
        self.bond_type.as_double()
    }

    /// `getValenceContrib(atom)`; 0 when `atom_idx` is not in the bond.
    pub fn valence_contrib(&self, _mol: &ROMol, atom_idx: u32) -> f64 {
        if atom_idx != self.begin && atom_idx != self.end {
            return 0.0;
        }
        // `QueryBond::getValenceContrib`: bonds with a complex bond type query
        // ("single or aromatic", ...) do not count towards the valence
        if let Some(q) = &self.query
            && crate::graphmol::queries::has_complex_bond_type_query(q)
        {
            return 0.0;
        }
        match self.bond_type {
            BondType::Dative | BondType::DativeOne if atom_idx != self.end => 0.0,
            BondType::Hydrogen => 0.0,
            t => t.as_double().unwrap_or(0.0),
        }
    }

    /// `Bond::Match`
    pub fn matches(&self, what: &Bond) -> bool {
        self.bond_type == BondType::Unspecified
            || what.bond_type == BondType::Unspecified
            || self.bond_type == what.bond_type
    }

    /// `invertChirality` (atropisomers only)
    pub fn invert_chirality(&mut self) -> bool {
        match self.stereo {
            BondStereo::StereoAtropCW => {
                self.stereo = BondStereo::StereoAtropCCW;
                true
            }
            BondStereo::StereoAtropCCW => {
                self.stereo = BondStereo::StereoAtropCW;
                true
            }
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bond_type_values() {
        assert_eq!(BondType::Aromatic.as_double().unwrap(), 1.5);
        assert_eq!(BondType::Aromatic.twice().unwrap(), 3);
        assert_eq!(BondType::Zero as u8, 21);
        assert_eq!(BondType::from_u8(12), Some(BondType::Aromatic));
        assert_eq!(BondStereo::StereoAtropCCW.name(), "STEREOATROPCCW");
        assert!(BondType::Other.as_double().is_err());
    }

    #[test]
    fn stereo_requires_atoms() {
        let mut b = Bond::new(BondType::Double);
        assert!(b.set_stereo(BondStereo::StereoCis).is_err());
        b.stereo_atoms = vec![0, 3];
        assert!(b.set_stereo(BondStereo::StereoCis).is_ok());
        assert!(!b.invert_chirality());
    }
}
