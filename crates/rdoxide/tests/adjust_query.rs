//! `adjust_query_properties`; expectations from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::adjust_query::{
    AdjustQueryParameters, adjust_flags::*, adjust_query_properties,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn describe(mol: &RWMol, idx: u32) -> String {
    mol.atom(idx)
        .query()
        .map_or_else(String::new, |q| q.description.clone())
}

#[test]
fn defaults_query_dummies_and_ring_atom_degrees() {
    let mol = smiles("*c1ccccc1");
    let adjusted = adjust_query_properties(&mol, &AdjustQueryParameters::default()).unwrap();
    assert_eq!(describe(&adjusted, 0), "AtomNull");
    assert!(adjusted.atom(1).has_query());
    assert!(!mol.atom(1).has_query(), "the input is not modified");
}

#[test]
fn flags_limit_the_adjusted_atoms() {
    let params = AdjustQueryParameters {
        adjust_degree_flags: ADJUST_IGNORERINGS,
        ..Default::default()
    };
    let adjusted = adjust_query_properties(&smiles("C1CC1C"), &params).unwrap();
    assert!(!adjusted.atom(0).has_query());
    assert!(adjusted.atom(3).has_query());
}

#[test]
fn no_adjustments_is_a_no_op() {
    let adjusted = adjust_query_properties(
        &smiles("*c1ccccc1"),
        &AdjustQueryParameters::no_adjustments(),
    )
    .unwrap();
    assert!(adjusted.atoms().iter().all(|a| !a.has_query()));
}

#[test]
fn conjugated_five_rings_get_bond_queries() {
    let params = AdjustQueryParameters {
        adjust_degree: false,
        adjust_conjugated_five_rings: true,
        ..Default::default()
    };
    let adjusted = adjust_query_properties(&smiles("C1=CC=CC1"), &params).unwrap();
    assert!(adjusted.bonds().iter().all(|b| b.has_query()));
}
