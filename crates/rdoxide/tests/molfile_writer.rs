//! MolBlock writer tests. `tests/data/molfiles/writer_expected.json` holds
//! inputs and the exact output of RDKit 2026.09.1 (the `molwrite`
//! differential mode covers full corpora).

use std::path::PathBuf;

use rdoxide::graphmol::file_parsers::mol_file_writer::{
    MolWriterParams, mol_to_mol_block, mol_to_v2k_mol_block, mol_to_v3k_mol_block,
};
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_block};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

#[test]
fn matches_rdkit_output() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data/molfiles/writer_expected.json");
    let cases: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let cases = cases.as_array().unwrap();
    assert!(cases.len() >= 19);
    for c in cases {
        let name = c["name"].as_str().unwrap();
        let input = c["input"].as_str().unwrap();
        let mol = mol_from_mol_block(input, &MolFileParserParams::default())
            .unwrap()
            .unwrap();
        let kw = c["kwargs"].as_object().unwrap();
        let params = MolWriterParams {
            include_stereo: kw
                .get("includeStereo")
                .and_then(|v| v.as_bool())
                .unwrap_or(true),
            kekulize: kw.get("kekulize").and_then(|v| v.as_bool()).unwrap_or(true),
            ..Default::default()
        };
        let out = match c["function"].as_str().unwrap() {
            "MolToMolBlock" => mol_to_mol_block(&mol, &params, -1),
            "MolToV3KMolBlock" => mol_to_v3k_mol_block(&mol, &params, -1),
            "MolToV2KMolBlock" => mol_to_v2k_mol_block(&mol, &params, -1),
            f => panic!("{f}"),
        }
        .unwrap();
        assert_eq!(out, c["expected"].as_str().unwrap(), "{name}");
    }
}

#[test]
fn without_coordinates() {
    let mol = mol_from_smiles("c1ccccc1O", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let p = MolWriterParams {
        include_stereo: false,
        ..Default::default()
    };
    let out = mol_to_mol_block(&mol, &p, -1).unwrap();
    assert!(out.starts_with("\n     RDKit          \n\n  7  7  0  0"));
    // 2D coordinate generation is not ported yet
    assert!(mol_to_mol_block(&mol, &MolWriterParams::default(), -1).is_err());
}
