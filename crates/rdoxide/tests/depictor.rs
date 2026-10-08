//! 2D depiction tests. `tests/data/depict_expected.json` holds RDKit
//! 2026.09.1 coordinates (4 decimals) for molecules whose depiction is the
//! same on arm64 and x86-64 builds; the `depict` differential mode checks
//! exact bits against the reference build of the host platform.

use std::path::PathBuf;

use rdoxide::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use rdoxide::graphmol::file_parsers::mol_file_writer::{MolWriterParams, mol_to_mol_block};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn r4(x: f64) -> f64 {
    // Python's round(x, 4) for the magnitudes used here
    let v = (x * 1e4).round() / 1e4;
    if v == 0.0 { 0.0 } else { v }
}

#[test]
fn matches_rdkit() {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/depict_expected.json");
    let cases: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    for c in cases.as_array().unwrap() {
        let smi = c["smiles"].as_str().unwrap();
        let mol = mol_from_smiles(smi, &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        for (key, canon) in [("canon", true), ("nocanon", false)] {
            let mut m = mol.clone();
            let p = Compute2DCoordParameters {
                canon_orient: canon,
                ..Default::default()
            };
            assert_eq!(compute_2d_coords(&mut m, &p).unwrap(), 0);
            let conf = m.conformer(-1).unwrap();
            assert!(!conf.is_3d());
            let exp = c[key].as_array().unwrap();
            for (i, e) in exp.iter().enumerate() {
                let p = conf.atom_pos(i as u32).unwrap();
                let (ex, ey) = (e[0].as_f64().unwrap(), e[1].as_f64().unwrap());
                assert!(
                    (r4(p.x) - ex).abs() < 1e-9 && (r4(p.y) - ey).abs() < 1e-9,
                    "{smi} {key} atom {i}: ({}, {}) vs ({ex}, {ey})",
                    p.x,
                    p.y
                );
            }
        }
        let block = mol_to_mol_block(&mol, &MolWriterParams::default(), -1).unwrap();
        assert_eq!(block, c["block"].as_str().unwrap(), "{smi}");
    }
}
