//! `graphmol::tautomer_query`; cases from upstream `catch_tests.cpp`.

use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use rdoxide::graphmol::tautomer_query::TautomerQuery;

fn mol(smi: &str) -> ROMol {
    mol_from_smiles(smi, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn template_error() {
    // the template must account for bonds modified when tautomers are sanitized
    let query = TautomerQuery::from_mol(&mol("Cc1nc2ccccc2[nH]1"), "").unwrap();
    let target = mol("CN1C2=C(C(=O)Nc3ccccc3)C(=O)CCN2c2ccccc21");
    let params = SubstructMatchParameters::default();
    assert!(
        query
            .tautomers()
            .iter()
            .any(|t| !substruct_match(&target, t, &params).unwrap().is_empty())
    );
    let (matches, tautomers) = query.substruct_of(&target, &params).unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(tautomers.len(), 1);
}

#[test]
fn uniquify() {
    let query = TautomerQuery::from_mol(&mol("O=C1CCCCC1"), "").unwrap();
    let target = mol("O=C1CCCC(CC)C1");
    let mut params = SubstructMatchParameters {
        uniquify: true,
        ..Default::default()
    };
    let (matches, tautomers) = query.substruct_of(&target, &params).unwrap();
    assert_eq!((matches.len(), tautomers.len()), (1, 1));
    params.uniquify = false;
    let (matches, tautomers) = query.substruct_of(&target, &params).unwrap();
    assert_eq!((matches.len(), tautomers.len()), (2, 2));
}

#[test]
fn keto_enol_match() {
    let query = TautomerQuery::from_mol(&mol("CC(=O)C"), "").unwrap();
    assert!(
        query
            .is_substruct_of(&mol("CC(O)=C"), &Default::default())
            .unwrap()
    );
    assert!(
        !query
            .is_substruct_of(&mol("CCCC"), &Default::default())
            .unwrap()
    );
}
