//! GETAWAY tests; expected values from RDKit 2026.09.1.

use std::path::PathBuf;

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::getaway::getaway;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::dict::RDProps;

fn ethanol() -> RWMol {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/molfiles/ethanol_3d.mol");
    let params = MolFileParserParams {
        remove_hs: false,
        ..MolFileParserParams::default()
    };
    mol_from_mol_file(path.to_str().unwrap(), &params)
        .unwrap()
        .unwrap()
}

#[test]
fn ethanol_getaway() {
    let mut m = ethanol();
    let g = getaway(&mut m, -1, 2, "").unwrap();
    assert_eq!(g.len(), 273);
    assert_eq!(&g[..6], &[4.755, 1.0, 2.958, 26.69, 3.0, 0.81]);
    assert_eq!(&g[269..], &[0.0, 0.0, 0.0, 2.051]);
    assert!(m.has_prop("_3DDistanceMatrix_Conf0"));
    assert!(m.has_prop("AdjacencyMatrix"));
}

#[test]
fn custom_property_gives_one_weighting() {
    let mut m = ethanol();
    let g = getaway(&mut m, -1, 2, "w").unwrap();
    assert_eq!(g.len(), 45);
    // the weight-independent head is the same as for the full descriptor
    assert_eq!(&g[..4], &[4.755, 1.0, 2.958, 26.69]);
}

#[test]
fn errors() {
    let mut flat = mol_from_smiles("CCO", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assert!(getaway(&mut flat, -1, 2, "").is_err());
    assert!(getaway(&mut ethanol(), 4, 2, "").is_err());
}
