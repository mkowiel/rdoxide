//! InChI tests against `tests/data/inchi_corpus.jsonl`, generated with the
//! official RDKit 2026.09.1 / InChI 1.07.3 by
//! `tools/differential/make_inchi_corpus.py`.
//!
//! InChIKeys are checked now (every InChI in the corpus, plus hand-written
//! valid and invalid InChI strings). `MolToInchi` / `MolBlockToInchi` are not
//! ported yet: `mol_to_inchi_matches_rdkit` asserts that they report
//! `NotImplemented` and turns into the exact-match check once the generation
//! pipeline lands (see docs/INCHI_PORT.md).

use std::path::PathBuf;

use rdoxide::external::inchi::{InchiKeyCheck, check_inchikey, get_std_inchikey_from_std_inchi};
use rdoxide::external::inchi_api::{
    ExtraInchiReturnValues, get_inchi_version, inchi_to_inchi_key, inchi_to_mol,
    mol_block_to_inchi, mol_to_inchi, mol_to_inchi_key,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use serde_json::Value;

fn corpus() -> Vec<Value> {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/inchi_corpus.jsonl");
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

fn s(v: &Value) -> &str {
    v.as_str().unwrap()
}

#[test]
fn version() {
    assert_eq!(get_inchi_version(), "1.07.3");
}

#[test]
fn inchi_key_matches_rdkit() {
    let recs = corpus();
    let mut n = 0;
    for r in &recs {
        let (inchi, want) = (s(&r["inchi"]), s(&r["key"]));
        let got = inchi_to_inchi_key(inchi);
        assert_eq!(got, want, "InchiToInchiKey({inchi:?})");
        if !got.is_empty() {
            let std = inchi.starts_with("InChI=1S/");
            assert_eq!(
                check_inchikey(&got),
                if std {
                    InchiKeyCheck::ValidStandard
                } else if inchi.starts_with("InChI=1B/") {
                    InchiKeyCheck::InvalidLayout
                } else {
                    InchiKeyCheck::ValidNonStandard
                },
                "{got}"
            );
            if std {
                assert_eq!(get_std_inchikey_from_std_inchi(inchi).unwrap(), got);
            }
        }
        n += 1;
    }
    assert!(n > 500, "corpus too small: {n}");
}

/// Exact `MolToInchi`/`MolBlockToInchi` output (InChI, return code, message,
/// AuxInfo, key) once the generation pipeline is ported; until then the API
/// must say that it is not implemented rather than return something wrong.
#[test]
fn mol_to_inchi_matches_rdkit() {
    for r in corpus() {
        let options = r["options"].as_str().filter(|o| !o.is_empty());
        let mut rv = ExtraInchiReturnValues::default();
        let got = match s(&r["kind"]) {
            "mol" => {
                let mol = mol_from_smiles(s(&r["smiles"]), &SmilesParserParams::default())
                    .unwrap()
                    .unwrap();
                let key = mol_to_inchi_key(&mol, options);
                let inchi = mol_to_inchi(&mol, &mut rv, options);
                if let (Ok(inchi), Ok(key)) = (&inchi, &key) {
                    assert_eq!(key, s(&r["key"]));
                    Ok(inchi.clone())
                } else {
                    inchi
                }
            }
            "molblock" => mol_block_to_inchi(s(&r["molblock"]), &mut rv, options),
            _ => continue,
        };
        match got {
            Ok(inchi) => {
                assert_eq!(inchi, s(&r["inchi"]), "{r}");
                assert_eq!(rv.return_code as i64, r["retcode"].as_i64().unwrap(), "{r}");
                assert_eq!(rv.message, s(&r["message"]), "{r}");
                assert_eq!(rv.aux_info, s(&r["aux"]), "{r}");
                // the usage text of an empty molecule carries the build date and platform
                if !s(&r["log"]).starts_with("InChI version 1, Software") {
                    assert_eq!(rv.log, s(&r["log"]), "{r}");
                }
            }
            Err(e) => panic!("{e}: {r}"),
        }
    }
}

#[test]
fn inchi_to_mol_methane() {
    let mut rv = ExtraInchiReturnValues::default();
    let m = inchi_to_mol("InChI=1S/CH4/h1H4", &mut rv, true, true)
        .unwrap()
        .unwrap();
    assert_eq!(m.num_atoms(), 1);
    assert_eq!(rv.return_code, 0);
}
