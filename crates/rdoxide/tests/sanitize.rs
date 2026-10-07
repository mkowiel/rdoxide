//! sanitizeMol tests; expected values from RDKit 2026.09.1 (cf. upstream
//! catch_graphmol.cpp "Sanitization tests" and molopstest.cpp).

use rdoxide::graphmol::mol_ops::{sanitize_flags::*, sanitize_mol};
use rdoxide::graphmol::smiles_parse::{mol_from_smiles, SmilesParserParams};
use rdoxide::graphmol::{BondType, HybridizationType, RWMol};
use rdoxide::{RDKitError, SanitizeError};

fn raw(s: &str) -> RWMol {
    let p = SmilesParserParams { sanitize: false, remove_hs: false, ..Default::default() };
    mol_from_smiles(s, &p).unwrap().unwrap()
}

fn sanitized(s: &str) -> RWMol {
    let mut m = raw(s);
    sanitize_mol(&mut m, SANITIZE_ALL).unwrap();
    m
}

#[test]
fn aromaticity_from_kekule() {
    let m = sanitized("C1=CC=CC=C1");
    assert!(m.atoms().iter().all(|a| a.is_aromatic() && a.hybridization() == HybridizationType::SP2));
    assert!(m.bonds().iter().all(|b| b.bond_type() == BondType::Aromatic && b.is_conjugated()));
    let m = sanitized("c1ccc2[nH]ccc2c1");
    assert!(m.atoms().iter().all(|a| a.is_aromatic()));
}

#[test]
fn cleanup_nitro_and_perchlorate() {
    let m = sanitized("CN(=O)=O");
    let charges: Vec<i32> = m.atoms().iter().map(|a| a.formal_charge()).collect();
    assert_eq!(charges, vec![0, 1, -1, 0]);
    assert_eq!(m.bond(1).bond_type(), BondType::Single);
    let m = sanitized("O=Cl(=O)(=O)[O-]");
    let charges: Vec<i32> = m.atoms().iter().map(|a| a.formal_charge()).collect();
    assert_eq!(charges, vec![-1, 3, -1, -1, -1]);
}

#[test]
fn hybridization() {
    let m = sanitized("CC=O");
    let h: Vec<HybridizationType> = m.atoms().iter().map(|a| a.hybridization()).collect();
    assert_eq!(h, vec![HybridizationType::SP3, HybridizationType::SP2, HybridizationType::SP2]);
}

#[test]
fn kekulize_failure() {
    let mut m = raw("c1cccc1");
    match sanitize_mol(&mut m, SANITIZE_ALL) {
        Err((flag, RDKitError::Sanitize(SanitizeError::Kekulize { atom_indices, message }))) => {
            assert_eq!(flag, SANITIZE_KEKULIZE);
            assert_eq!(atom_indices, vec![0, 1, 2, 3, 4]);
            assert_eq!(message, "Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 4");
        }
        other => panic!("unexpected {other:?}"),
    }
}

#[test]
fn valence_failure() {
    let mut m = raw("C(C)(C)(C)(C)C");
    match sanitize_mol(&mut m, SANITIZE_ALL) {
        Err((flag, RDKitError::Sanitize(SanitizeError::AtomValence { atom_idx, .. }))) => {
            assert_eq!(flag, SANITIZE_PROPERTIES);
            assert_eq!(atom_idx, 0);
        }
        other => panic!("unexpected {other:?}"),
    }
}
