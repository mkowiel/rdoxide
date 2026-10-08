//! SmilesMolSupplier tests; behaviour checked against RDKit 2026.09.1.

use rdoxide::graphmol::file_parsers::smiles_mol_supplier::{
    SmilesMolSupplier, SmilesMolSupplierParams,
};
use rdoxide::rdgeneral::{RDProps, RDValue};

const TABLE: &[u8] = b"smiles id logP\nCCO e1 0.5\n# comment\n\nc1ccccc1 b1 2.1\nbadsmiles x 1\n";

fn supplier(text: &[u8]) -> SmilesMolSupplier {
    let mut s = SmilesMolSupplier::default();
    s.set_data(text, SmilesMolSupplierParams::default());
    s
}

fn string_prop(m: &impl RDProps, name: &str) -> String {
    match m.get_prop(name) {
        Some(RDValue::String(s)) => s.clone(),
        other => panic!("{name}: {other:?}"),
    }
}

#[test]
fn sequential_reads_set_name_and_properties() {
    let mut s = supplier(TABLE);
    let first = s.next_mol().unwrap().unwrap();
    assert_eq!(first.num_atoms(), 3);
    assert_eq!(string_prop(&first, "_Name"), "e1");
    assert_eq!(string_prop(&first, "logP"), "0.5");
    assert_eq!(s.next_mol().unwrap().unwrap().num_atoms(), 6);
    // unparsable record: a result without a molecule, not an error
    assert!(s.next_mol().unwrap().is_none());
    assert!(s.at_end().unwrap());
}

#[test]
fn random_access_and_length() {
    let mut s = supplier(TABLE);
    assert_eq!(s.length().unwrap(), 3);
    assert_eq!(s.get(1).unwrap().unwrap().num_atoms(), 6);
    assert_eq!(s.get_item_text(2).unwrap(), "badsmiles x 1");
    assert!(s.get(3).is_err());
}

#[test]
fn reset_rewinds() {
    let mut s = supplier(TABLE);
    s.next_mol().unwrap();
    s.next_mol().unwrap();
    s.reset().unwrap();
    assert_eq!(s.next_mol().unwrap().unwrap().num_atoms(), 3);
}

#[test]
fn no_input_is_an_error() {
    let mut s = SmilesMolSupplier::default();
    assert!(s.length().is_err());
    assert!(s.at_end().is_err());
}
