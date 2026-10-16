//! Tests of `GeneralizedSubstruct` (ported from upstream's `catch_tests.cpp`
//! and the Python wrapper tests).

use rdoxide::graphmol::adjust_query::AdjustQueryParameters;
use rdoxide::graphmol::generalized_substruct::{
    ExtendedQueryMol, create_extended_query_mol, has_substruct_match_xq, substruct_match_xq,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::SubstructMatchParameters;

fn parse_smiles(s: &str) -> Result<rdoxide::graphmol::RWMol, ()> {
    Ok(mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap())
}

fn xq(smi: &str, enumerate: bool, tautomers: bool) -> ExtendedQueryMol {
    let m = parse_smiles(smi).unwrap();
    create_extended_query_mol(
        &m,
        enumerate,
        tautomers,
        false,
        &AdjustQueryParameters::default(),
    )
    .unwrap()
}

#[test]
fn link_nodes_and_tautomers_select_the_variant() {
    assert!(matches!(
        xq("COCC1OC(N)=N1", false, false),
        ExtendedQueryMol::Mol(_)
    ));
    assert!(matches!(
        xq("COCC1OC(N)=N1", false, true),
        ExtendedQueryMol::Tautomer(_)
    ));
    assert!(matches!(
        xq("COCc1n[nH]c(F)c1 |LN:1:1.3|", true, false),
        ExtendedQueryMol::Bundle(_)
    ));
    assert!(matches!(
        xq("COCC1OC(N)=N1 |LN:1:1.3|", true, true),
        ExtendedQueryMol::TautomerBundle(_)
    ));
}

#[test]
fn search_and_round_trip() {
    let q = xq("COCc1n[nH]c(F)c1 |LN:1:1.3|", true, true);
    let target = parse_smiles("COOCc1n[nH]c(F)c1").unwrap();
    let ps = SubstructMatchParameters::default();
    let matches = substruct_match_xq(&target, &q, &ps).unwrap();
    assert_eq!(matches.len(), 1);
    let atoms: Vec<i32> = matches[0].iter().map(|p| p.1).collect();
    assert_eq!(atoms, [0, 1, 2, 3, 4, 5, 9, 6, 7, 8]);

    let q2 = ExtendedQueryMol::from_binary(&q.to_binary().unwrap()).unwrap();
    assert!(has_substruct_match_xq(&target, &q2, &ps).unwrap());
    let q3 = ExtendedQueryMol::from_json(&q.to_json().unwrap()).unwrap();
    assert!(has_substruct_match_xq(&target, &q3, &ps).unwrap());
    assert!(ExtendedQueryMol::from_binary(b"junk").is_err());
    assert!(ExtendedQueryMol::from_json("{}").is_err());
}
