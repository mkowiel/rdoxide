//! Reaction enumeration strategies and `EnumerateLibrary`.

use std::collections::HashSet;

use rdoxide::graphmol::chem_reactions::daylight::{
    ReactionSmartsParserParams, reaction_from_smarts,
};
use rdoxide::graphmol::chem_reactions::enumerate::{
    EnumerateLibrary, EnumerationParams, EnumerationStrategy, compute_num_products,
};
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn mol(smiles: &str) -> ROMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn product_counts() {
    assert_eq!(compute_num_products(&[10, 5, 6]), 300);
    assert_eq!(compute_num_products(&[]), 1);
    assert_eq!(compute_num_products(&[u64::MAX, 2]), u64::MAX);
}

#[test]
fn cartesian_product_walks_the_first_index_fastest() {
    let mut s = EnumerationStrategy::cartesian_product();
    s.initialize(&[2, 3]).unwrap();
    let mut seen = Vec::new();
    while s.has_next() {
        seen.push(s.advance().unwrap().to_vec());
    }
    assert_eq!(seen, [[0, 0], [1, 0], [0, 1], [1, 1], [0, 2], [1, 2]]);
    assert_eq!(s.num_permutations(), 6);
}

#[test]
fn random_strategies_never_stop_and_are_reproducible() {
    for make in [
        EnumerationStrategy::random_sample,
        EnumerationStrategy::random_sample_all_bbs,
    ] {
        let run = || {
            let mut s = make();
            s.initialize(&[10, 5, 6]).unwrap();
            (0..100)
                .map(|_| s.advance().unwrap().to_vec())
                .collect::<Vec<_>>()
        };
        let (a, b) = (run(), run());
        assert_eq!(a, b);
        assert!(a.iter().all(|p| p[0] < 10 && p[1] < 5 && p[2] < 6));
    }
}

#[test]
fn all_bbs_strategy_uses_every_building_block_early() {
    let mut s = EnumerationStrategy::random_sample_all_bbs();
    s.initialize(&[10, 5, 6]).unwrap();
    let first: Vec<Vec<u64>> = (0..10).map(|_| s.advance().unwrap().to_vec()).collect();
    for col in 0..3 {
        let used: HashSet<u64> = first.iter().map(|p| p[col]).collect();
        assert_eq!(used.len() as u64, [10, 5, 6][col].min(10));
    }
}

#[test]
fn even_pairs_gives_unique_positions_and_stats() {
    let mut s = EnumerationStrategy::even_sample_pairs();
    s.initialize(&[4, 4, 4]).unwrap();
    let mut seen = HashSet::new();
    for _ in 0..16 {
        assert!(seen.insert(s.advance().unwrap().to_vec()));
    }
    let stats = s.stats().unwrap();
    assert!(stats.contains("PAIRSTAT"));
    assert!(EnumerationStrategy::cartesian_product().stats().is_none());
}

#[test]
fn library_enumerates_all_matching_reagents() {
    let rxn = reaction_from_smarts(
        "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]",
        &ReactionSmartsParserParams::new(),
    )
    .unwrap();
    let acids = ["CC(=O)O", "OC(=O)c1ccccc1", "CCC"].map(mol).to_vec();
    let amines = ["N", "CN"].map(mol).to_vec();
    let mut lib = EnumerateLibrary::new(
        &rxn,
        &vec![acids, amines],
        None,
        &EnumerationParams::default(),
    )
    .unwrap();
    // the alkane does not match the acid template
    assert_eq!(lib.reagents()[0].len(), 2);
    let mut n = 0;
    while lib.has_next().unwrap() {
        assert_eq!(lib.next_smiles().unwrap().len(), 1);
        n += 1;
    }
    assert_eq!(n, 4);
    assert!(lib.next_products().is_err());
    lib.reset_state().unwrap();
    assert!(lib.has_next().unwrap());
}
