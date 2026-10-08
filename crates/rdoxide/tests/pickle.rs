//! MolPickler tests (cf. upstream testPickler.cpp). Byte-level compatibility
//! with RDKit is checked by the `pickle` differential mode.

use std::path::PathBuf;

use rdoxide::graphmol::pickle::{flags, mol_from_pickle, mol_from_pickle_stream, pickle_mol};
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn smi(s: &str) -> rdoxide::graphmol::RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn canon(m: &mut rdoxide::graphmol::RWMol) -> String {
    mol_to_smiles(m, &SmilesWriteParams::default()).unwrap()
}

#[test]
fn roundtrip() {
    for s in [
        "C[C@H](O)c1ccccc1",
        "F/C=C/Cl",
        "[Na+].[Cl-]",
        "[CH3:1][13CH2]O",
        &"C".repeat(300),
    ] {
        let mut m = smi(s);
        let mut back =
            mol_from_pickle(&pickle_mol(&m, flags::ALL_PROPS).unwrap(), flags::ALL_PROPS).unwrap();
        assert_eq!(canon(&mut m), canon(&mut back), "{s}");
        assert_eq!(m.ring_info().atom_rings(), back.ring_info().atom_rings());
    }
}

#[test]
fn header_and_props() {
    let mut m = smi("CCO");
    m.set_prop("foo", 3i32, false).unwrap();
    let b = pickle_mol(&m, flags::NO_PROPS).unwrap();
    assert_eq!(&b[..4], &0xDEADBEEFu32.to_le_bytes());
    assert_eq!(&b[8..20], &[16, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0]);
    assert!(
        !mol_from_pickle(&b, flags::ALL_PROPS)
            .unwrap()
            .has_prop("foo")
    );
    let b = pickle_mol(&m, flags::MOL_PROPS).unwrap();
    let back = mol_from_pickle(&b, flags::ALL_PROPS).unwrap();
    assert_eq!(back.get_prop("foo"), Some(&RDValue::Int(3)));
    // properties are skipped when not requested on read
    assert!(
        !mol_from_pickle(&b, flags::NO_PROPS)
            .unwrap()
            .has_prop("foo")
    );
}

#[test]
fn queries_survive() {
    let q = mol_from_smarts("[$(c-C(=O)O)]~[c;!R0]", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    let back =
        mol_from_pickle(&pickle_mol(&q, flags::NO_PROPS).unwrap(), flags::ALL_PROPS).unwrap();
    let m = smi("OC(=O)c1ccccc1");
    let ps = SubstructMatchParameters::default();
    assert_eq!(
        substruct_match(&m, &q, &ps).unwrap(),
        substruct_match(&m, &back, &ps).unwrap()
    );
    assert!(!substruct_match(&m, &back, &ps).unwrap().is_empty());
}

#[test]
fn bad_pickles() {
    assert!(mol_from_pickle(b"xx", flags::ALL_PROPS).is_err());
    assert!(mol_from_pickle(&[0u8; 32], flags::ALL_PROPS).is_err());
}

/// upstream testPickler.cpp `test1`: an old (v2) pickle stream versus SMILES
#[test]
fn old_pickle_stream() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../upstream/rdkit/Code/GraphMol/test_data");
    let (Ok(pkl), Ok(smiles)) = (
        std::fs::read(base.join("canonSmiles.v2.pkl")),
        std::fs::read_to_string(base.join("canonSmiles.smi")),
    ) else {
        eprintln!("upstream checkout not found; skipping");
        return;
    };
    let mut pos = 0;
    for (count, line) in smiles.lines().enumerate().take(100) {
        let s = line.split('\t').next().unwrap();
        let mut m1 = mol_from_smiles(s, &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        let mut m2 = mol_from_pickle_stream(&pkl, &mut pos, flags::ALL_PROPS).unwrap();
        assert_eq!(canon(&mut m1), canon(&mut m2), "line {count}");
    }
}
