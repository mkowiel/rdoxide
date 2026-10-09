//! WHIM and EEM tests; expected values from RDKit 2026.09.1.

use std::path::PathBuf;

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::eem::eem;
use rdoxide::graphmol::descriptors::whim::whim;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

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

fn close(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() <= tol, "{a:?} != {b:?}");
    }
}

#[test]
fn ethanol_whim() {
    let w = whim(&ethanol(), -1, 0.001, "").unwrap();
    assert_eq!(w.len(), 114);
    close(&w[..6], &[1.405, 0.599, 0.363, 0.593, 0.253, 0.24], 1e-9);
    close(&w[110..], &[4.137, 3.06, 4.436, 3.394], 1e-9);
    assert_eq!(whim(&ethanol(), -1, 0.001, "w").unwrap().len(), 17);
}

#[test]
fn ethanol_charges() {
    let q = eem(&ethanol(), -1).unwrap();
    assert_eq!(q.len(), 9);
    close(
        &q[..4],
        &[
            -0.5309430092431018,
            -0.013404997902780557,
            -0.7234167614413173,
            0.19630285976213874,
        ],
        1e-12,
    );
    assert!(q.iter().sum::<f64>().abs() < 1e-12);
}

#[test]
fn molecules_without_conformers_are_rejected() {
    let flat = mol_from_smiles("CCO", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assert!(whim(&flat, -1, 0.001, "").is_err());
    assert!(eem(&flat, -1).is_err());
    assert!(whim(&ethanol(), 9, 0.001, "").is_err());
}
