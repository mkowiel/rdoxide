//! Extended reduced graphs and ErG fingerprints; expectations from the reference RDKit.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::reduced_graphs::{
    ERG_ATOM_TYPES_PROP, generate_erg_fingerprint_for_reduced_graph,
    generate_mol_extended_reduced_graph, get_erg_fingerprint,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn mol(smiles: &str) -> RWMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn non_zero(fp: &[f64]) -> Vec<(usize, f64)> {
    fp.iter()
        .enumerate()
        .filter(|(_, v)| **v != 0.0)
        .map(|(i, v)| (i, (v * 1e6).round() / 1e6))
        .collect()
}

#[test]
fn feature_types_of_the_reduced_graph() {
    let rg = generate_mol_extended_reduced_graph(&mol("OCCc1ccccc1"), None).unwrap();
    let types: Vec<Vec<i32>> = rg
        .atoms()
        .iter()
        .map(|a| match a.get_prop(ERG_ATOM_TYPES_PROP) {
            Some(RDValue::VecInt(v)) => v.clone(),
            other => panic!("unexpected {other:?}"),
        })
        .collect();
    // O: donor + acceptor; the ring atoms are gone, the aromatic centroid stays
    assert_eq!(types, [vec![0, 1], vec![], vec![], vec![], vec![5]]);
}

#[test]
fn fingerprints_are_independent_of_the_ring_type() {
    let a = get_erg_fingerprint(&mol("OCCc1ccccc1"), None, 0.3, 1, 15).unwrap();
    assert_eq!(a.len(), 315);
    assert_eq!(
        non_zero(&a),
        [
            (77, 0.3),
            (78, 1.0),
            (79, 0.3),
            (152, 0.3),
            (153, 1.0),
            (154, 0.3)
        ]
    );
    // a cyclohexane ring is aliphatic (hydrophobic feature) instead
    let b = get_erg_fingerprint(&mol("OCCC1CCCCC1"), None, 0.3, 1, 15).unwrap();
    assert_eq!(
        non_zero(&b),
        [
            (62, 0.3),
            (63, 1.0),
            (64, 0.3),
            (137, 0.3),
            (138, 1.0),
            (139, 0.3)
        ]
    );
}

#[test]
fn path_range_and_fuzz_change_the_vector() {
    let m = mol("OCCc1ccccc1");
    let fp = get_erg_fingerprint(&m, None, 0.1, 2, 6).unwrap();
    assert_eq!(fp.len(), 105);
    assert_eq!(
        non_zero(&fp),
        [
            (26, 0.1),
            (27, 1.0),
            (28, 0.1),
            (51, 0.1),
            (52, 1.0),
            (53, 0.1)
        ]
    );
    // the reduced graph can be fingerprinted directly
    let rg = generate_mol_extended_reduced_graph(&m, None).unwrap();
    let from_graph = generate_erg_fingerprint_for_reduced_graph(&rg, None, 0.1, 2, 6).unwrap();
    assert_eq!(fp, from_graph);
}

#[test]
fn bad_path_range_is_an_error() {
    let rg = generate_mol_extended_reduced_graph(&mol("CCO"), None).unwrap();
    assert!(generate_erg_fingerprint_for_reduced_graph(&rg, None, 0.3, 5, 5).is_err());
}
