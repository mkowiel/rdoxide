//! Multithreaded suppliers; expectations from RDKit's `testMultithreadedMolSupplier.py`
//! and the reference build.

use std::path::PathBuf;

use rdoxide::graphmol::file_parsers::mol_file_parser::MolFileParserParams;
use rdoxide::graphmol::file_parsers::multithreaded_supplier::{
    MultithreadedMolSupplier, SupplierParams,
};
use rdoxide::graphmol::file_parsers::smiles_mol_supplier::SmilesMolSupplierParams;
use rdoxide::rdgeneral::{RDProps, common_properties as cp};

fn data(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/molfiles");
    path.join(name).to_string_lossy().into_owned()
}

fn threads(n: i32) -> SupplierParams {
    SupplierParams {
        num_writer_threads: n,
        ..Default::default()
    }
}

fn prop(mol: &impl RDProps, key: &str) -> String {
    mol.get_prop(key)
        .map(|v| v.to_rd_string())
        .unwrap_or_default()
}

fn names(supplier: &mut MultithreadedMolSupplier) -> Vec<Option<String>> {
    let mut res = Vec::new();
    while !supplier.at_end() {
        let mol = supplier.next_mol().flatten();
        res.push(mol.map(|m| prop(&m, cp::_Name)));
    }
    res
}

#[test]
fn sd_molecules_come_back_in_file_order_whatever_the_thread_count() {
    let expected = [
        "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
        "220", "229", "256",
    ];
    for n in [1, 2, 4, 0] {
        let mut s = MultithreadedMolSupplier::from_sd_file(
            &data("NCI_aids_few.sdf"),
            MolFileParserParams::default(),
            threads(n),
        )
        .unwrap();
        let got = names(&mut s);
        assert_eq!(
            got.iter().flatten().map(String::as_str).collect::<Vec<_>>(),
            expected
        );
        assert_eq!(s.last_record_id(), 16);
        assert!(s.last_item_text().contains("$$$$"));
    }
}

#[test]
fn a_record_that_cannot_be_parsed_yields_none() {
    // the data item without a label swallows the rest of its record (upstream's loop)
    let mut s = MultithreadedMolSupplier::from_sd_file(
        &data("sdErrors1.sdf"),
        MolFileParserParams::default(),
        threads(1),
    )
    .unwrap();
    assert_eq!(names(&mut s), [None]);
}

#[test]
fn smiles_records_properties_and_the_final_empty_record() {
    let params = SmilesMolSupplierParams {
        delimiter: ",".into(),
        smiles_column: 1,
        name_column: 0,
        title_line: false,
        ..Default::default()
    };
    let mut s = MultithreadedMolSupplier::from_smiles_file(&data("fewSmi.csv"), params, threads(3))
        .unwrap();
    let mut mols = Vec::new();
    while !s.at_end() {
        mols.push(s.next_mol().flatten());
    }
    // the file ends with a newline: the empty last record is a failure
    assert_eq!(mols.len(), 11);
    assert!(mols[10].is_none());
    let first = mols[0].as_ref().unwrap();
    assert_eq!(prop(first, cp::_Name), "1");
    assert_eq!(prop(first, "Column_2"), "34.14");
}

#[test]
fn smiles_without_a_name_column_are_named_after_upstreams_line_number() {
    let params = SmilesMolSupplierParams {
        delimiter: ",".into(),
        name_column: -1,
        ..Default::default()
    };
    let mut s =
        MultithreadedMolSupplier::from_smiles_file(&data("first_200.tpsa.csv"), params, threads(2))
            .unwrap();
    let first = s.next_mol().flatten().unwrap();
    assert_eq!(prop(&first, cp::_Name), "4294967295");
    assert_eq!(prop(&first, "TPSA"), "34.14");
}

#[test]
fn at_end_is_false_before_the_first_request_and_closing_stops_the_threads() {
    let mut s = MultithreadedMolSupplier::from_sd_file(
        &data("NCI_aids_few.sdf"),
        MolFileParserParams::default(),
        threads(2),
    )
    .unwrap();
    assert!(!s.at_end());
    assert!(s.next_mol().is_some());
    s.close();
    assert!(s.at_end());
    assert!(s.next_mol().is_none());
    drop(s); // must not hang
}

#[test]
fn a_missing_file_is_an_error() {
    assert!(
        MultithreadedMolSupplier::from_sd_file(
            "/does/not/exist.sdf",
            MolFileParserParams::default(),
            SupplierParams::default()
        )
        .is_err()
    );
}
