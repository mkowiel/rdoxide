//! Substructure matching tests; expected values from RDKit 2026.09.1 (cf.
//! upstream Substruct/catch_tests.cpp and testSubstructMatch.cpp).

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::{
    FinalMatchCheck, MatchCheck, SubstructMatchParameters, substruct_match, substruct_match_count,
};

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
fn atoms(ms: &[Vec<(i32, i32)>]) -> Vec<Vec<i32>> {
    ms.iter().map(|m| m.iter().map(|p| p.1).collect()).collect()
}

#[test]
fn basic() {
    let m = smi("OC(=O)c1ccccc1C(=O)O");
    let ps = SubstructMatchParameters::default();
    let r = substruct_match(&m, &sma("C(=O)[OH]"), &ps).unwrap();
    assert_eq!(atoms(&r), vec![vec![1, 2, 0], vec![9, 10, 11]]);
    assert!(substruct_match(&m, &sma("[#7]"), &ps).unwrap().is_empty());
    assert_eq!(substruct_match_count(&m, &sma("c"), &ps).unwrap(), 6);
}

#[test]
fn uniquify_and_order() {
    let m = smi("c1ccccc1");
    let q = smi("c1ccccc1");
    let ps = SubstructMatchParameters {
        uniquify: false,
        ..Default::default()
    };
    let r = substruct_match(&m, &q, &ps).unwrap();
    assert_eq!(r.len(), 12);
    assert_eq!(
        atoms(&r[..2]),
        vec![vec![0, 1, 2, 3, 4, 5], vec![0, 5, 4, 3, 2, 1]]
    );
    assert_eq!(
        substruct_match(&m, &q, &SubstructMatchParameters::default())
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn recursion() {
    let m = smi("OC(=O)c1ccccc1N");
    let ps = SubstructMatchParameters::default();
    assert_eq!(
        atoms(&substruct_match(&m, &sma("[$(c-C(=O)O)]"), &ps).unwrap()),
        vec![vec![3]]
    );
    assert_eq!(
        atoms(&substruct_match(&m, &sma("[c;!$(c-[#7]);!$(c-C)]"), &ps).unwrap()),
        vec![vec![4], vec![5], vec![6], vec![7]]
    );
    let nr = SubstructMatchParameters {
        recursion_possible: false,
        ..Default::default()
    };
    assert!(
        substruct_match(&m, &sma("[$(c-C(=O)O)]"), &nr)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn chirality() {
    let m = smi("C[C@H](N)C(=O)O");
    let ch = SubstructMatchParameters {
        use_chirality: true,
        ..Default::default()
    };
    assert!(
        !substruct_match(&m, &smi("C[C@H](N)C=O"), &ch)
            .unwrap()
            .is_empty()
    );
    assert!(
        substruct_match(&m, &smi("C[C@@H](N)C=O"), &ch)
            .unwrap()
            .is_empty()
    );
    assert!(
        !substruct_match(
            &m,
            &smi("C[C@@H](N)C=O"),
            &SubstructMatchParameters::default()
        )
        .unwrap()
        .is_empty()
    );
    let db = smi("F/C=C/F");
    assert!(
        substruct_match(&db, &smi("F/C=C\\F"), &ch)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn extra_atom_check_filters_and_overrides() {
    // the callback sees every matched pair: only atoms 1 and 2 are accepted
    let m = smi("CCOCC");
    let q = smi("CO");
    let first_carbon: MatchCheck = MatchCheck::new(|_, _, _, mi| mi == 1 || mi == 2);
    let ps = SubstructMatchParameters {
        extra_atom_check: Some(first_carbon.clone()),
        ..Default::default()
    };
    assert_eq!(
        atoms(&substruct_match(&m, &q, &ps).unwrap()),
        vec![vec![1, 2]]
    );

    // overriding the default check ignores element labels
    let ps = SubstructMatchParameters {
        extra_atom_check: Some(MatchCheck::new(|_, qi, _, mi| {
            qi == 0 && mi == 0 || qi == 1
        })),
        extra_atom_check_overrides_default: true,
        ..Default::default()
    };
    let r = substruct_match(&m, &q, &ps).unwrap();
    assert_eq!(atoms(&r), vec![vec![0, 1]]);
}

#[test]
fn extra_bond_check_filters() {
    let m = smi("CC=CC");
    let q = smi("CC");
    let ps = SubstructMatchParameters {
        extra_bond_check: Some(MatchCheck::new(|_, _, mol, mb| {
            mol.bond(mb).bond_type_as_double().unwrap() == 1.0
        })),
        ..Default::default()
    };
    let r = substruct_match(&m, &q, &ps).unwrap();
    assert_eq!(atoms(&r), vec![vec![0, 1], vec![2, 3]]);
}

#[test]
fn extra_final_check_validates_matches() {
    let m = smi("CCOCC");
    let q = smi("CCO");
    let all = SubstructMatchParameters::default();
    assert_eq!(substruct_match(&m, &q, &all).unwrap().len(), 2);
    let none = SubstructMatchParameters {
        extra_final_check: Some(FinalMatchCheck::new(|_, _| false)),
        ..Default::default()
    };
    assert!(substruct_match(&m, &q, &none).unwrap().is_empty());
    // the closure sees the matched atoms in query order
    let large = SubstructMatchParameters {
        extra_final_check: Some(FinalMatchCheck::new(|_, aids| aids.iter().sum::<u32>() > 5)),
        ..Default::default()
    };
    assert_eq!(
        atoms(&substruct_match(&m, &q, &large).unwrap()),
        vec![vec![4, 3, 2]]
    );
}
