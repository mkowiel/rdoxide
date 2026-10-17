//! `substruct_has_match` must agree with the VF2 search on whether a match exists.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use rdoxide::graphmol::substruct_exists::substruct_has_match;

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}
fn sma(s: &str) -> RWMol {
    mol_from_smarts(s, &SmartsParserParams::default())
        .unwrap()
        .unwrap()
}

const MOLS: &[&str] = &[
    "CC(=O)Nc1ccccc1",
    "OC(=O)c1ccccc1C(=O)O",
    "c1ccc2ccccc2c1",
    "c1ccc2c(c1)oc1ccccc12",
    "CCCCCCS",
    "CCN(CC)CCOC(=O)c1ccc(N)cc1",
    "C1CC1C(F)(F)F",
    "[Na+].[O-]C(=O)C",
    "C[C@H](N)C(=O)O",
    "C1=CC=CC=C1",
    "CC(C)(C)c1ccc(O)cc1",
    "ClC(Cl)(Cl)Cl",
    "O=C1CCCCC1",
    "c1ccncc1",
    "CC",
    "C",
];

const QUERIES: &[&str] = &[
    "c1ccccc1",
    "C(=O)[OH]",
    "C(=O)N",
    "[NX3;H2,H1;!$(NC=O)]",
    "[OH]c",
    "c1ccc2ccccc2c1",
    "[#6]~[#6]~[#6]~[#6]~[#6]~[#6]",
    "[$(C=O)]N",
    "CC(=O)Nc1ccccc1",
    "c1ccc2c(c1)oc1ccccc12",
    "*",
    "c[F,Cl,Br,I]",
    "[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#16]",
    "*~*~*~*~*~*~*",
    "[#7]~[#6]~[#6]~[#6]~[#6]~[#6]~[#15]",
    "c1ccccc1-c1ccccc1",
    "[#6]C(=O)O[#6]",
    "C.C",
    "[C;R]",
    "[Cl]C([Cl])[Cl]",
    "[#6][#6]([#9])([#9])[#9]",
    "[C@H]",
    "[Na+]",
    "[!#6]~[#6]~[!#6]",
    "C1CC1",
    "[#8]=[#6]-[#7]",
    "[R2]",
    "[a;r6]",
];

#[test]
fn existence_agrees_with_the_vf2_search() {
    let mut hits = 0;
    for m in MOLS.iter().map(|s| smi(s)) {
        for q in QUERIES.iter().map(|s| sma(s)) {
            let ps = SubstructMatchParameters {
                max_matches: 1,
                ..SubstructMatchParameters::default()
            };
            let expected = !substruct_match(&m, &q, &ps).unwrap().is_empty();
            assert_eq!(substruct_has_match(&m, &q, &ps).unwrap(), expected);
            hits += expected as usize;
        }
    }
    assert!(hits > 30, "the test molecules must give both answers");
}

#[test]
fn smiles_queries_and_odd_parameters_agree() {
    for q in ["c1ccccc1", "CC(=O)O", "C.C", "[Na+]", "O=CN"].map(smi) {
        for m in MOLS.iter().map(|s| smi(s)) {
            for (chir, qq) in [(false, false), (true, false), (false, true)] {
                let ps = SubstructMatchParameters {
                    max_matches: 1,
                    use_chirality: chir,
                    use_query_query_matches: qq,
                    ..SubstructMatchParameters::default()
                };
                assert_eq!(
                    substruct_has_match(&m, &q, &ps).unwrap(),
                    !substruct_match(&m, &q, &ps).unwrap().is_empty()
                );
            }
        }
    }
}
