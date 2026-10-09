//! SubstructLibrary: holders, filtering, thread-count independence, search
//! order, keys and serialization (cases follow upstream `rough_test.py`).

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::SubstructMatchParameters;
use rdoxide::graphmol::substruct_library::{
    KeyFromPropHolder, MolHolder, PatternHolder, SubstructLibrary, shared,
};
use rdoxide::rdgeneral::RDProps;

fn mol(smiles: &str) -> RWMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn holders() -> Vec<MolHolder> {
    vec![
        MolHolder::Mols(Vec::new()),
        MolHolder::Binary(Vec::new()),
        MolHolder::Smiles(Vec::new()),
        MolHolder::TrustedSmiles(Vec::new()),
    ]
}

fn library(holder: MolHolder, patterns: bool, smiles: &[&str]) -> SubstructLibrary {
    let fps = patterns.then(|| shared(PatternHolder::default()));
    let lib = SubstructLibrary::new(shared(holder), fps, None);
    for smi in smiles {
        lib.add_mol(&mol(smi)).unwrap();
    }
    lib
}

const SMILES: [&str; 8] = [
    "c1ccccc1",
    "CCO",
    "c1ccccc1O",
    "CCCC",
    "C1CCCCC1",
    "c1ccncc1",
    "OCCO",
    "CC(=O)O",
];

fn matches(lib: &SubstructLibrary, smarts: &str, threads: i32, max: i32) -> Vec<u32> {
    let q = mol_from_smarts(smarts, &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    lib.get_matches(
        &q,
        0,
        lib.len(),
        &SubstructMatchParameters::default(),
        threads,
        max,
    )
    .unwrap()
}

#[test]
fn every_holder_finds_the_same_matches() {
    for holder in holders() {
        for patterns in [false, true] {
            let lib = library(holder.clone(), patterns, &SMILES);
            assert_eq!(matches(&lib, "c1ccccc1", 1, -1), [0, 2]);
            assert_eq!(matches(&lib, "[OX2H]", 1, -1), [1, 2, 6, 7]);
            assert_eq!(matches(&lib, "[#6]-[#6]-[#8]", 1, -1), [1, 6, 7]);
            assert_eq!(matches(&lib, "[r6]", 1, -1), [0, 2, 4, 5]);
        }
    }
}

#[test]
fn results_do_not_depend_on_the_thread_count() {
    let many: Vec<&str> = SMILES.iter().cycle().take(103).copied().collect();
    let lib = library(MolHolder::default(), true, &many);
    for smarts in ["[OX2H]", "c", "C", "[#7]"] {
        let expected = matches(&lib, smarts, 1, -1);
        for max in [-1, 1, 3, 7, 1000] {
            let single = matches(&lib, smarts, 1, max);
            for threads in [2, 3, 4, 8] {
                let multi = matches(&lib, smarts, threads, max);
                assert_eq!(multi, single, "{smarts} threads={threads} max={max}");
            }
            let want: Vec<u32> = expected
                .iter()
                .copied()
                .take(if max < 0 { usize::MAX } else { max as usize })
                .collect();
            assert_eq!(single, want);
        }
    }
}

#[test]
fn count_and_has_match() {
    let lib = library(MolHolder::default(), true, &SMILES);
    let params = SubstructMatchParameters::default();
    let q = mol_from_smarts("c1ccccc1", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(lib.count_matches(&q, 0, 8, &params, 2).unwrap(), 2);
    assert_eq!(lib.count_matches(&q, 1, 8, &params, 2).unwrap(), 1);
    assert!(lib.has_match(&q, 0, 8, &params, -1).unwrap());
    assert!(!lib.has_match(&q, 3, 8, &params, 1).unwrap());
    // no work for no results
    assert!(matches(&lib, "C", 1, 0).is_empty());
}

#[test]
fn bad_ranges_are_rejected() {
    let lib = library(MolHolder::default(), false, &SMILES);
    let q = mol("C");
    let params = SubstructMatchParameters::default();
    assert!(lib.get_matches(&q, 8, 9, &params, 1, -1).is_err());
    assert!(lib.get_matches(&q, 3, 3, &params, 1, -1).is_err());
    assert!(
        SubstructLibrary::default()
            .get_matches(&q, 0, 1, &params, 1, -1)
            .is_err()
    );
    assert!(lib.get_mol(8).is_err());
}

#[test]
fn search_order_controls_the_result_order() {
    let lib_smiles = ["CCCOC", "CCCCOCC", "CCOC", "COC", "CCCCCOC"];
    let mut lib = library(MolHolder::default(), false, &lib_smiles);
    lib.set_search_order(vec![3, 2, 0, 1, 4]).unwrap();
    assert_eq!(matches(&lib, "COC", 1, 2), [3, 2]);
    assert_eq!(matches(&lib, "COC", 1, -1), [3, 2, 0, 1, 4]);
    // several threads sort by position in the search order
    assert_eq!(matches(&lib, "COC", 3, 2), [3, 2]);
    assert!(lib.set_search_order(vec![5]).is_err());
    lib.set_search_order(Vec::new()).unwrap();
    assert_eq!(matches(&lib, "COC", 1, -1), [0, 1, 2, 3, 4]);
}

#[test]
fn query_with_ring_features_on_trusted_smiles() {
    let lib = library(
        MolHolder::TrustedSmiles(Vec::new()),
        false,
        &["C1CCCC1C", "CCCCC", "c1ccccc1"],
    );
    assert_eq!(matches(&lib, "[C;R]", 1, -1), [0]);
    assert_eq!(matches(&lib, "[R]", 2, -1), [0, 2]);
    assert_eq!(matches(&lib, "[C;!R]", 1, -1), [0, 1]);
}

#[test]
fn unparsable_smiles_never_match() {
    let mut holder = MolHolder::Smiles(Vec::new());
    holder.add_smiles("&not smiles").unwrap();
    holder.add_smiles("CCO").unwrap();
    let lib = SubstructLibrary::new(shared(holder), None, None);
    assert_eq!(matches(&lib, "C", 1, -1), [1]);
    assert!(lib.get_mol(0).unwrap().is_none());
}

#[test]
fn keys_come_from_a_property() {
    let keys = shared(KeyFromPropHolder::new("id"));
    let lib = SubstructLibrary::new(shared(MolHolder::default()), None, Some(keys.clone()));
    for (i, smi) in ["CCO", "CCN", "CCC"].iter().enumerate() {
        let mut m = mol(smi);
        if i != 1 {
            m.set_prop("id", format!("Z{i}"), false).unwrap();
        }
        lib.add_mol(&m).unwrap();
    }
    let keys = keys.read().unwrap();
    assert_eq!(keys.get_keys(&[2, 0, 1]).unwrap(), ["Z2", "Z0", "LIBIDX-1"]);
    assert!(keys.get_key(3).is_err());
}

#[test]
fn add_patterns_fingerprints_an_existing_library() {
    for tautomeric in [false, true] {
        let mut lib = library(MolHolder::Smiles(Vec::new()), false, &SMILES);
        let holder = if tautomeric {
            PatternHolder::new_tautomeric(1024)
        } else {
            PatternHolder::new(1024)
        };
        let expected = holder.make_fingerprint(&mol(SMILES[2])).unwrap();
        lib.add_patterns(Some(shared(holder)), -1).unwrap();
        let fps = lib.fp_holder().unwrap().read().unwrap();
        assert_eq!(fps.len(), lib.len());
        assert_eq!(fps.get_fingerprint(2).unwrap(), &expected);
        drop(fps);
        assert!(lib.add_patterns(None, 1).is_err());
        assert_eq!(matches(&lib, "c1ccccc1", 2, -1), [0, 2]);
    }
}

#[test]
fn serialization_round_trips() {
    for holder in holders() {
        let mut lib = SubstructLibrary::new(
            shared(holder),
            Some(shared(PatternHolder::new_tautomeric(512))),
            Some(shared(KeyFromPropHolder::default())),
        );
        for smi in SMILES {
            lib.add_mol(&mol(smi)).unwrap();
        }
        lib.set_search_order((0..8).rev().collect()).unwrap();
        let copy = SubstructLibrary::from_text(&lib.to_text().unwrap()).unwrap();
        assert_eq!(copy.len(), lib.len());
        assert_eq!(copy.search_order(), lib.search_order());
        assert!(copy.fp_holder().unwrap().read().unwrap().is_tautomeric());
        assert_eq!(
            matches(&copy, "[OX2H]", 2, -1),
            matches(&lib, "[OX2H]", 2, -1)
        );
    }
    assert!(SubstructLibrary::from_bytes(b"garbage").is_err());
}
