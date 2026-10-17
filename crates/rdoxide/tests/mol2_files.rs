//! Tripos Mol2 reading.

use std::path::Path;

use rdoxide::graphmol::file_parsers::mol2::{Mol2ParserParams, mol_from_mol2_file};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};

const DIR: &str = "../../upstream/rdkit/Code/GraphMol/FileParsers/test_data/";

fn have_data() -> bool {
    Path::new(DIR).join("3505.mol2").exists()
}

#[test]
fn pyrazole_pyridine() {
    if !have_data() {
        return;
    }
    let mut m = mol_from_mol2_file(
        &format!("{DIR}pyrazole_pyridine.mol2"),
        &Mol2ParserParams::default(),
    )
    .unwrap()
    .unwrap();
    assert_eq!(m.num_atoms(), 5);
    assert_eq!(
        mol_to_smiles(&mut m, &SmilesWriteParams::default()).unwrap(),
        "c1cn[nH]c1"
    );
}

#[test]
fn cleanup_substructures_toggle() {
    if !have_data() {
        return;
    }
    let file = format!("{DIR}3505.mol2");
    let m = mol_from_mol2_file(&file, &Mol2ParserParams::default())
        .unwrap()
        .unwrap();
    let b = m.bond_idx_between(3, 12).unwrap();
    assert_eq!(
        m.bond(b).bond_type(),
        rdoxide::graphmol::bond::BondType::Single
    );
    assert_eq!(m.atom(12).formal_charge(), 0);

    let params = Mol2ParserParams {
        cleanup_substructures: false,
        ..Mol2ParserParams::default()
    };
    let m = mol_from_mol2_file(&file, &params).unwrap().unwrap();
    let b = m.bond_idx_between(3, 12).unwrap();
    assert_eq!(
        m.bond(b).bond_type(),
        rdoxide::graphmol::bond::BondType::Double
    );
    assert_eq!(m.atom(12).formal_charge(), 1);
}

#[test]
fn missing_file_is_bad_file() {
    let e = mol_from_mol2_file("/no/such.mol2", &Mol2ParserParams::default()).unwrap_err();
    assert!(matches!(e, rdoxide::RDKitError::BadFile(_)));
}
