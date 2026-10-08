//! SDMolSupplier tests (cf. upstream testMolSupplier.cpp). Expected values
//! come from RDKit 2026.09.1; the `sdf` differential mode covers the
//! upstream SD files.

use std::path::PathBuf;

use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::file_parsers::sd_mol_supplier::SDMolSupplier;
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn data(f: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data/molfiles")
        .join(f)
        .to_string_lossy()
        .into_owned()
}

fn prop(m: &rdoxide::graphmol::RWMol, k: &str) -> String {
    match m.get_prop(k) {
        Some(RDValue::String(s)) => s.clone(),
        other => panic!("{k}: {other:?}"),
    }
}

#[test]
fn random_access_and_iteration() {
    let mut s =
        SDMolSupplier::from_file(&data("first_200.props.sdf"), MolFileParserParams::default())
            .unwrap();
    assert_eq!(s.length().unwrap(), 200);
    let mut last = s.get(199).unwrap().unwrap();
    assert_eq!(
        mol_to_smiles(&mut last, &SmilesWriteParams::default()).unwrap(),
        "CC(=O)Nc1cc(C)cc(C)c1"
    );
    let m3 = s.get(3).unwrap().unwrap();
    assert_eq!(prop(&m3, "NUM_RINGS"), "1");
    assert!(matches!(s.get(200), Err(RDKitError::FileParse(_))));
    s.reset().unwrap();
    let mut n = 0;
    while !s.at_end().unwrap() {
        assert!(s.next_mol().unwrap().is_some());
        n += 1;
    }
    assert_eq!(n, 200);
    let text = s.get_item_text(1).unwrap();
    assert!(text.starts_with(b"\n     RDKit          2D\n\n 20 23"));
}

#[test]
fn properties() {
    let mut s =
        SDMolSupplier::from_file(&data("BlankPropLines.sdf"), MolFileParserParams::default())
            .unwrap();
    let m = s.next_mol().unwrap().unwrap();
    assert_eq!(prop(&m, "MultiLineProperty1"), "foo\nbar\n \nbaz");
    assert_eq!(prop(&m, "MultiLineProperty2"), "foo\n>  <fooprop>\nbaz\n ");

    // an unterminated data header swallows the rest of the file
    let mut s = SDMolSupplier::from_file(
        &data("esters.prop_name_trunc.sdf"),
        MolFileParserParams::default(),
    )
    .unwrap();
    // length() indexes `$$$$` lines, so it counts all six records
    assert_eq!(s.length().unwrap(), 6);
    let mut n = 0;
    while !s.at_end().unwrap() {
        s.next_mol().unwrap();
        n += 1;
    }
    assert_eq!(n, 2);
}

#[test]
fn bad_input() {
    assert!(matches!(
        SDMolSupplier::from_file("/nonexistent.sdf", MolFileParserParams::default()),
        Err(RDKitError::BadFile(_))
    ));
    let mut s = SDMolSupplier::default();
    assert!(s.next_mol().is_err());
    s.set_data(b"", MolFileParserParams::default());
    assert_eq!(s.length().unwrap(), 0);
    assert!(s.at_end().unwrap());
}
