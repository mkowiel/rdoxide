//! `replace_substructs`, `replace_sidechains`, `replace_core` and
//! `murcko_decompose`; expectations from RDKit 2026.09.1.

use rdoxide::graphmol::chem_transforms::{
    murcko_decompose, replace_core, replace_core_with_match, replace_sidechains, replace_substructs,
};
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{ROMol, RWMol};

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn smarts(s: &str) -> RWMol {
    mol_from_smarts(s, &SmartsParserParams::default())
        .unwrap()
        .unwrap()
}

fn can(mol: &mut ROMol) -> String {
    mol_to_smiles(mol, &SmilesWriteParams::default()).unwrap()
}

fn cans(mols: Vec<RWMol>) -> Vec<String> {
    mols.into_iter().map(|mut m| can(&mut m)).collect()
}

#[test]
fn replace_substructs_one_result_per_match() {
    let res = replace_substructs(
        &smiles("COCCOC"),
        &smarts("O[CH3]"),
        &smiles("NC"),
        false,
        0,
        false,
    )
    .unwrap();
    assert_eq!(cans(res), ["CNCCOC", "CNCCOC"]);
    let res = replace_substructs(
        &smiles("CCOC"),
        &smarts("[CH3]O"),
        &smiles("NC"),
        false,
        0,
        false,
    )
    .unwrap();
    assert_eq!(cans(res), ["CC.CN"]);
}

#[test]
fn replace_substructs_all_at_once_and_connection_point() {
    let res = replace_substructs(
        &smiles("COCCOC"),
        &smarts("O[CH3]"),
        &smiles("CN"),
        true,
        1,
        false,
    )
    .unwrap();
    assert_eq!(cans(res), ["CNCCNC"]);
    assert!(
        replace_substructs(&smiles("CCOC"), &smarts("O"), &smiles("N"), false, 5, false).is_err()
    );
}

#[test]
fn replace_substructs_without_match_copies() {
    let res =
        replace_substructs(&smiles("CCCC"), &smarts("N"), &smiles("O"), false, 0, false).unwrap();
    assert_eq!(cans(res), ["CCCC"]);
}

#[test]
fn replace_sidechains_labels_attachment_points() {
    let core = smiles("C1CCC1");
    let mut res = replace_sidechains(&smiles("CCC1CCC1"), &core, false)
        .unwrap()
        .unwrap();
    assert_eq!(can(&mut res), "[1*]C1CCC1");
    let mut res = replace_sidechains(&smiles("C1CC2C1CCC2"), &core, false)
        .unwrap()
        .unwrap();
    assert_eq!(can(&mut res), "[1*]C1CCC1[2*]");
    assert!(
        replace_sidechains(&smiles("CCC1CC1"), &core, false)
            .unwrap()
            .is_none()
    );
}

#[test]
fn replace_core_labels_and_requirements() {
    let core = smiles("C1CCC1");
    let mut res = replace_core(&smiles("CCC1CCC1"), &core, true, false, false, false)
        .unwrap()
        .unwrap();
    assert_eq!(can(&mut res), "[1*]CC");
    let mut res = replace_core(
        &smiles("CCN1CCC1"),
        &smiles("C1CCN1"),
        true,
        true,
        false,
        false,
    )
    .unwrap()
    .unwrap();
    assert_eq!(can(&mut res), "[3*]CC");
    assert!(
        replace_core(&smiles("CCC1CC1"), &core, true, false, false, false)
            .unwrap()
            .is_none()
    );
    assert!(
        replace_core(
            &smiles("CC1CC1"),
            &smiles("C1CC1"),
            true,
            false,
            true,
            false
        )
        .unwrap()
        .is_none()
    );
}

#[test]
fn replace_core_with_explicit_match() {
    let mol = smiles("C1CNCC1");
    let core = smarts("[*]N[*]");
    let m = vec![(0, 1), (1, 2), (2, 3)];
    let mut res = replace_core_with_match(&mol, &core, &m, true, false, false)
        .unwrap()
        .unwrap();
    assert_eq!(can(&mut res), "[1*]CC[2*]");
    let bad = vec![(0, 1), (1, 2), (5, 3)];
    assert!(replace_core_with_match(&mol, &core, &bad, true, false, false).is_err());
    let bad = vec![(0, 1), (1, 2), (2, 30)];
    assert!(replace_core_with_match(&mol, &core, &bad, true, false, false).is_err());
}

#[test]
fn murcko_decompose_keeps_rings_and_linkers() {
    for (input, scaffold) in [
        ("c1ccccc1CCc1ccncc1", "c1ccc(CCc2ccncc2)cc1"),
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "c1ccccc1"),
        ("CCO", ""),
        ("C[C@H](F)C1CCC(=O)CC1", "O=C1CCCCC1"),
        ("c1ccc2[nH]ccc2c1CCN", "c1ccc2[nH]ccc2c1"),
    ] {
        let mut mol = smiles(input);
        let mut res = murcko_decompose(&mut mol).unwrap();
        assert_eq!(can(&mut res), scaffold, "{input}");
    }
}
