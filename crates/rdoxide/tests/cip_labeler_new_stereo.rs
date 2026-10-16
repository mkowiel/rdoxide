//! CIP labeling of molecules with the new stereo perception; a binary of its
//! own, since the perception is a global setting. Cases from upstream
//! `CIPLabeler/catch_tests.cpp`.

use rdoxide::graphmol::chirality::set_use_legacy_stereo_perception;
use rdoxide::graphmol::cip_labeler::assign_cip_labels;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn cip(smiles: &str) -> Vec<Option<String>> {
    let mut m = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assign_cip_labels(&mut m, 0).unwrap();
    (0..m.num_atoms())
        .map(|i| match m.atom(i).get_prop("_CIPCode") {
            Some(RDValue::String(s)) => Some(s.clone()),
            _ => None,
        })
        .collect()
}

#[test]
fn para_stereochemistry_with_pseudo_asymmetry() {
    set_use_legacy_stereo_perception(false);
    let s = |l: &str| Some(l.to_string());
    let labels = cip("OC(=O)[C@H]1CC[C@@H](CC1)O[C@@H](F)Cl");
    assert_eq!(
        (labels[3].clone(), labels[6].clone(), labels[10].clone()),
        (s("S"), s("r"), s("S"))
    );

    let labels = cip("OC(=O)[C@H]1CC[C@@H](CC1)O[CH](Cl)Cl");
    assert_eq!(
        (labels[3].clone(), labels[6].clone(), labels[10].clone()),
        (s("r"), s("r"), None)
    );
}
