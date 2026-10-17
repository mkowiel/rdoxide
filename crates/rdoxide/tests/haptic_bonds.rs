//! Haptic bonds <-> dative bonds.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::bond::BondType;
use rdoxide::graphmol::mol_ops::haptic::{dative_bonds_to_haptic, haptic_bonds_to_dative};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::common_properties as cp;
use rdoxide::rdgeneral::{RDProps, RDValue};

/// A cyclopentadienyl ring and an iron atom, joined by one dative bond from
/// each ring atom.
fn ferrocene_half() -> RWMol {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let mut mol = mol_from_smiles("C1=CC=C[CH-]1.[Fe+2]", &params)
        .unwrap()
        .unwrap();
    for c in 0..5 {
        mol.add_bond(c, 5, BondType::Dative).unwrap();
    }
    mol
}

#[test]
fn dative_bonds_become_one_haptic_bond() {
    let mut mol = ferrocene_half();
    dative_bonds_to_haptic(&mut mol).unwrap();
    assert_eq!(mol.num_atoms(), 7);
    assert_eq!(mol.atom(6).atomic_num(), 0);
    let dative: Vec<_> = mol
        .bonds()
        .iter()
        .filter(|b| b.bond_type() == BondType::Dative)
        .collect();
    assert_eq!(dative.len(), 1);
    // the atoms are listed in breadth-first order from the first one
    assert_eq!(
        dative[0].get_prop(cp::_MolFileBondEndPts),
        Some(&RDValue::String("(5 1 2 5 3 4)".into()))
    );
    assert_eq!(
        dative[0].get_prop(cp::_MolFileBondAttach),
        Some(&RDValue::String("ALL".into()))
    );
}

#[test]
fn round_trip_restores_the_dative_bonds() {
    let mut mol = ferrocene_half();
    dative_bonds_to_haptic(&mut mol).unwrap();
    haptic_bonds_to_dative(&mut mol).unwrap();
    assert_eq!(mol.num_atoms(), 6);
    let dative = mol
        .bonds()
        .iter()
        .filter(|b| b.bond_type() == BondType::Dative)
        .count();
    assert_eq!(dative, 5);
    for b in mol
        .bonds()
        .iter()
        .filter(|b| b.bond_type() == BondType::Dative)
    {
        assert_eq!(b.end_atom_idx(), 5);
    }
}

#[test]
fn single_dative_bond_is_not_haptic() {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let mut mol = mol_from_smiles("N.[Fe+2]", &params).unwrap().unwrap();
    mol.add_bond(0, 1, BondType::Dative).unwrap();
    dative_bonds_to_haptic(&mut mol).unwrap();
    assert_eq!(mol.num_atoms(), 2);
    haptic_bonds_to_dative(&mut mol).unwrap();
    assert_eq!(mol.num_bonds(), 1);
}
