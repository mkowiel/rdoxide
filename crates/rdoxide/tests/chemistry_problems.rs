//! `detect_chemistry_problems` and the atom/bond SMILES parsers.

use rdoxide::SanitizeError;
use rdoxide::graphmol::atom::Atom;
use rdoxide::graphmol::bond::BondType;
use rdoxide::graphmol::mol_ops::detect_chemistry_problems;
use rdoxide::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use rdoxide::graphmol::smiles_parse::{
    SmilesParserParams, atom_from_smiles, bond_from_smiles, mol_from_smiles,
};

#[test]
fn every_problem_is_reported() {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let mol = mol_from_smiles("CFCc1cc1FC", &params).unwrap().unwrap();
    let problems = detect_chemistry_problems(&mol, SANITIZE_ALL);
    assert_eq!(problems.len(), 3);
    assert!(matches!(
        problems[0],
        SanitizeError::AtomValence { atom_idx: 1, .. }
    ));
    assert!(matches!(
        problems[1],
        SanitizeError::AtomValence { atom_idx: 6, .. }
    ));
    match &problems[2] {
        SanitizeError::Kekulize { atom_indices, .. } => assert_eq!(atom_indices, &[3, 4, 5]),
        other => panic!("unexpected {other:?}"),
    }
}

#[test]
fn a_sound_molecule_has_no_problems() {
    let mol = mol_from_smiles("c1ccccc1O", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assert!(detect_chemistry_problems(&mol, SANITIZE_ALL).is_empty());
}

#[test]
fn atoms_and_bonds_from_smiles() {
    let a: Atom = atom_from_smiles("C").unwrap();
    assert_eq!(a.atomic_num(), 6);
    assert_eq!(atom_from_smiles("[13CH3]").unwrap().isotope(), 13);
    assert!(atom_from_smiles("error").is_none());
    assert!(atom_from_smiles("").is_none());
    assert_eq!(bond_from_smiles("=").unwrap().bond_type(), BondType::Double);
    assert_eq!(bond_from_smiles("-").unwrap().bond_type(), BondType::Single);
    assert!(bond_from_smiles("d").is_none());
}
