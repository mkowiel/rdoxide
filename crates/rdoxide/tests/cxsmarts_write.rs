//! CXSMARTS writing; expectations from RDKit 2026.09.1
//! (`SmilesParse/catch_tests.cpp` and `cxsmiles_test.cpp`).

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::smarts_write::{
    mol_fragment_to_cx_smarts, mol_to_cx_smarts, mol_to_smarts,
};
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smarts(s: &str) -> RWMol {
    mol_from_smarts(s, &SmartsParserParams::default())
        .unwrap()
        .unwrap()
}

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn cx(mol: &mut RWMol) -> String {
    mol_to_cx_smarts(mol, &SmilesWriteParams::default()).unwrap()
}

#[test]
fn atom_labels_follow_the_output_order() {
    let mut m = smarts("CCC |$foo;;bar$|");
    assert_eq!(
        mol_to_smarts(&mut m, &SmilesWriteParams::default()).unwrap(),
        "CCC"
    );
    assert_eq!(cx(&mut m), "CCC |$foo;;bar$|");
    // the atoms are written in a different order from the input
    let mut m = smarts("CC1.OC1 |$foo;;;bar$|");
    assert_eq!(cx(&mut m), "CCCO |$foo;;bar;$|");
}

#[test]
fn enhanced_stereo_is_written() {
    let mut m = smarts("[#6][C@]([#8])(F)Cl |&1:1|");
    assert_eq!(cx(&mut m), "[#6][C@]([#8])(F)Cl |&1:1|");
}

#[test]
fn nothing_to_add() {
    let mut m = smarts("C[C@H]([F,Cl,Br])[C@H](C)[C@@H](C)Br");
    assert_eq!(cx(&mut m), "C[C@&H1]([F,Cl,Br])[C@&H1](C)[C@@&H1](C)Br");
    let mut m = smarts("C[C@H]([F,Cl,Br])[C@H](C)[C@@H](C)Br |a:1,o1:4,5|");
    assert_eq!(
        cx(&mut m),
        "C[C@&H1]([F,Cl,Br])[C@&H1](C)[C@@&H1](C)Br |a:1,o1:4,5|"
    );
    let mut empty = RWMol::new();
    assert_eq!(cx(&mut empty), "");
}

#[test]
fn dative_bonds_are_single_bonds_with_an_extension() {
    let mut m = smiles("[NH3]->[Fe]-[NH2]");
    assert_eq!(cx(&mut m), "[#7]-[Fe]-[#7] |C:0.0|");
    let mut m = smiles("[NH3][Fe][NH3]");
    assert_eq!(cx(&mut m), "[#7]-[Fe]-[#7] |C:0.0,2.1|");
}

#[test]
fn fragments() {
    let mut m = smarts("CCC |$foo;bar;baz$|");
    let res =
        mol_fragment_to_cx_smarts(&mut m, &SmilesWriteParams::default(), &[1, 2], None).unwrap();
    assert_eq!(res, "CC |$bar;baz$|");
    assert!(mol_fragment_to_cx_smarts(&mut m, &SmilesWriteParams::default(), &[], None).is_err());
}
