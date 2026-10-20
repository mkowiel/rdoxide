//! The CommonChem / RDKit JSON molecule interchange format (`MolJSON`).
//!
//! Ported from RDKit `Code/GraphMol/MolInterchange/{MolInterchange.h,
//! details.h,Parser.cpp,Writer.cpp}` (2026.09.1). Note that, as upstream says,
//! the format and API should be considered beta.
//!
//! Copyright (C) 2018-2022 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

mod parser;
mod writer;

pub use parser::json_to_mols;
pub use writer::{mol_to_json, mols_to_json};

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::stereo_group::StereoGroupType;

/// Version of the plain CommonChem format.
const CURRENT_MOL_JSON_VERSION: i64 = 10;
/// Version of the RDKit flavour of the format.
const CURRENT_RDKIT_JSON_VERSION: i64 = 12;
const CURRENT_RDKIT_REPRESENTATION_VERSION: i64 = 2;
const CURRENT_CHARGE_REPRESENTATION_VERSION: i64 = 10;
const CURRENT_QUERY_REPRESENTATION_VERSION: i64 = 10;

/// `JSONParseParameters`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JsonParseParameters {
    /// set the bond type of bonds flagged aromatic to `Aromatic`
    pub set_aromatic_bonds: bool,
    /// perform reasonable valence checks
    pub strict_valence_check: bool,
    /// extract molecular properties
    pub parse_properties: bool,
    /// extract conformers
    pub parse_conformers: bool,
    /// use the implicit H counts from the JSON (turn off when parsing queries)
    pub use_h_counts: bool,
}

impl Default for JsonParseParameters {
    fn default() -> Self {
        JsonParseParameters {
            set_aromatic_bonds: true,
            strict_valence_check: false,
            parse_properties: true,
            parse_conformers: true,
            use_h_counts: true,
        }
    }
}

/// `JSONWriteParameters`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JsonWriteParameters {
    /// use RDKit extensions to the CommonChem format
    pub use_rdkit_extensions: bool,
}

impl Default for JsonWriteParameters {
    fn default() -> Self {
        JsonWriteParameters {
            use_rdkit_extensions: true,
        }
    }
}

// ------------------------------------------------------------------ lookups

const CHIRAL_TAGS: [(&str, ChiralType); 4] = [
    ("unspecified", ChiralType::Unspecified),
    ("cw", ChiralType::TetrahedralCW),
    ("ccw", ChiralType::TetrahedralCCW),
    ("other", ChiralType::Other),
];

fn chiral_tag_from_name(name: &str) -> Option<ChiralType> {
    CHIRAL_TAGS
        .iter()
        .find(|(n, _)| *n == name)
        .map(|&(_, t)| t)
}

fn chiral_tag_name(tag: ChiralType) -> Option<&'static str> {
    CHIRAL_TAGS.iter().find(|(_, t)| *t == tag).map(|&(n, _)| n)
}

const BOND_ORDERS: [(i64, BondType); 6] = [
    (0, BondType::Zero),
    (1, BondType::Single),
    (2, BondType::Double),
    (3, BondType::Triple),
    (4, BondType::Quadruple),
    (17, BondType::Dative),
];

pub(crate) fn bond_type_from_order(order: i64) -> Option<BondType> {
    BOND_ORDERS
        .iter()
        .find(|(o, _)| *o == order)
        .map(|&(_, t)| t)
}

pub(crate) fn bond_order(bond_type: BondType) -> Option<i64> {
    BOND_ORDERS
        .iter()
        .find(|(_, t)| *t == bond_type)
        .map(|&(o, _)| o)
}

fn bond_stereo_from_name(name: &str) -> Option<BondStereo> {
    Some(match name {
        "unspecified" => BondStereo::StereoNone,
        "cis" => BondStereo::StereoCis,
        "trans" => BondStereo::StereoTrans,
        "either" => BondStereo::StereoAny,
        "atrop_cw" => BondStereo::StereoAtropCW,
        "atrop_ccw" => BondStereo::StereoAtropCCW,
        _ => return None,
    })
}

fn bond_stereo_name(stereo: BondStereo) -> &'static str {
    match stereo {
        BondStereo::StereoNone => "unspecified",
        BondStereo::StereoCis | BondStereo::StereoZ => "cis",
        BondStereo::StereoTrans | BondStereo::StereoE => "trans",
        BondStereo::StereoAny => "either",
        BondStereo::StereoAtropCW => "atrop_cw",
        BondStereo::StereoAtropCCW => "atrop_ccw",
    }
}

fn stereo_group_type_from_name(name: &str) -> Option<StereoGroupType> {
    Some(match name {
        "abs" => StereoGroupType::Absolute,
        "and" => StereoGroupType::And,
        "or" => StereoGroupType::Or,
        _ => return None,
    })
}

fn stereo_group_type_name(group_type: StereoGroupType) -> &'static str {
    match group_type {
        StereoGroupType::Absolute => "abs",
        StereoGroupType::And => "and",
        StereoGroupType::Or => "or",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::romol::ROMol;
    use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    fn mol(smiles: &str) -> ROMol {
        mol_from_smiles(smiles, &SmilesParserParams::default())
            .unwrap()
            .unwrap()
    }

    fn smiles(mol: &mut ROMol) -> String {
        mol_to_smiles(mol, &SmilesWriteParams::default()).unwrap()
    }

    #[test]
    fn round_trips_stereo_and_charges() {
        for smi in ["C[C@H](F)Cl", "F/C=C/F", "c1ccccc1[NH3+]", "[13CH3][O-]"] {
            let m = mol(smi);
            let json = mol_to_json(&m, &JsonWriteParameters::default()).unwrap();
            let mut back = json_to_mols(&json, &JsonParseParameters::default()).unwrap();
            assert_eq!(back.len(), 1);
            assert_eq!(smiles(&mut back[0]), smiles(&mut mol(smi)), "{smi}");
        }
    }

    #[test]
    fn plain_commonchem_has_no_extensions() {
        let m = mol("C[C@H](O)C[C@@H](C)F |o1:1,4|");
        let rdkit = mol_to_json(&m, &JsonWriteParameters::default()).unwrap();
        assert!(rdkit.contains("\"stereoGroups\""));
        let plain = JsonWriteParameters {
            use_rdkit_extensions: false,
        };
        let json = mol_to_json(&m, &plain).unwrap();
        assert!(json.starts_with("{\"commonchem\":{\"version\":10}"));
        assert!(!json.contains("stereoGroups"));
    }

    #[test]
    fn many_molecules_keep_their_order() {
        let mols: Vec<ROMol> = (1..=40).map(|n| mol(&"C".repeat(n))).collect();
        let refs: Vec<_> = mols.iter().map(Some).collect();
        let json = mols_to_json(&refs, &JsonWriteParameters::default()).unwrap();
        let back = json_to_mols(&json, &JsonParseParameters::default()).unwrap();
        let sizes: Vec<_> = back.iter().map(ROMol::num_atoms).collect();
        assert_eq!(sizes, (1..=40).collect::<Vec<u32>>());
    }

    #[test]
    fn rejects_bad_documents() {
        let params = JsonParseParameters::default();
        for bad in [
            "{",
            "[]",
            "{\"a\":1}",
            "{\"commonchem\":{\"version\":11}}",
            "{\"commonchem\":10,\"molecules\":[]}",
        ] {
            assert!(json_to_mols(bad, &params).is_err(), "{bad}");
        }
        let null = [None];
        assert!(mols_to_json(&null, &JsonWriteParameters::default()).is_err());
    }
}
