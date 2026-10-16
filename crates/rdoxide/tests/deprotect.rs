//! Deprotection; cases from upstream `Deprotect/deprotectTest.cpp` and
//! `Wrap/rough_test.py`.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::deprotect::{
    DeprotectData, deprotect, deprotect_in_place, get_deprotections,
};
use rdoxide::graphmol::smiles_parse::write::mol_to_smiles;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn parse(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn canon(m: &RWMol) -> String {
    mol_to_smiles(&mut m.clone(), &Default::default()).unwrap()
}

const BOC2: &str = "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C";

#[test]
fn every_example_deprotects() {
    let mut n = 0;
    for data in get_deprotections().iter().filter(|d| !d.example.is_empty()) {
        let (start, end) = data.example.split_once(">>").unwrap();
        let m = deprotect(&parse(start), std::slice::from_ref(data)).unwrap();
        assert_eq!(canon(&m), canon(&parse(end)), "{}", data.full_name);
        let mut m = parse(start);
        deprotect_in_place(&mut m, std::slice::from_ref(data)).unwrap();
        assert_eq!(canon(&m), canon(&parse(end)), "in place {}", data.full_name);
        n += 1;
    }
    assert!(n > 0);
}

#[test]
fn defaults_are_annotated() {
    let m = deprotect(&parse(BOC2), get_deprotections()).unwrap();
    assert_eq!(canon(&m), canon(&parse("NCc1ccccc1N")));
    assert_eq!(
        m.get_prop("DEPROTECTIONS"),
        Some(&RDValue::VecString(vec!["Boc".into(), "Boc".into()]))
    );
    assert_eq!(m.get_prop("DEPROTECTION_COUNT"), Some(&RDValue::Int(2)));

    let mut m = parse(BOC2);
    assert!(deprotect_in_place(&mut m, get_deprotections()).unwrap());
    assert_eq!(canon(&m), canon(&parse("NCc1ccccc1N")));
    assert!(!deprotect_in_place(&mut m, get_deprotections()).unwrap());
}

#[test]
fn custom_and_empty_lists() {
    let boc = DeprotectData::new(
        "amine",
        "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]",
        "Boc",
        "tert-butyloxycarbonyl",
        "",
    )
    .unwrap();
    assert!(boc.is_valid());
    let m = deprotect(&parse(BOC2), &[]).unwrap();
    assert_eq!(canon(&m), canon(&parse(BOC2)));
    assert_eq!(m.get_prop("DEPROTECTION_COUNT"), Some(&RDValue::Int(0)));
    let m = deprotect(&parse(BOC2), &[boc]).unwrap();
    assert_eq!(canon(&m), canon(&parse("NCc1ccccc1N")));
}

#[test]
fn two_products_are_invalid() {
    let d = DeprotectData::new("x", "[N:1]>>[N:1].C", "x", "x", "").unwrap();
    assert!(!d.is_valid());
}
