//! SMILES writer tests; expected values from RDKit 2026.09.1 (cf. upstream
//! catch_smiles.cpp / smiTest.cpp).

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::smiles_parse::write::{
    SmilesWriteParams, mol_fragment_to_smiles, mol_to_smiles,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn parse(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn canon(s: &str) -> String {
    mol_to_smiles(&mut parse(s), &SmilesWriteParams::default()).unwrap()
}

#[test]
fn canonical_smiles() {
    assert_eq!(canon("OCC"), "CCO");
    assert_eq!(canon("C1=CC=CC=C1O"), "Oc1ccccc1");
    assert_eq!(
        canon("OC(=O)[C@@H](N)Cc1ccccc1"),
        "N[C@@H](Cc1ccccc1)C(=O)O"
    );
    assert_eq!(canon("Cl\\C=C\\F"), "F/C=C/Cl");
    assert_eq!(canon("[Na+].[Cl-].C1CC1"), "C1CC1.[Cl-].[Na+]");
    assert_eq!(canon("C12C3C4C1C5C2C3C45"), "C12C3C4C1C1C2C3C41");
    assert_eq!(canon("[NH3]->[Pt]"), "[NH3]->[Pt]");
    assert_eq!(canon(""), "");
}

#[test]
fn writer_options() {
    let mut m = parse("OC(=O)[C@@H](N)Cc1ccccc1");
    let p = |f: fn(&mut SmilesWriteParams)| {
        let mut p = SmilesWriteParams::default();
        f(&mut p);
        p
    };
    assert_eq!(
        mol_to_smiles(&mut m, &p(|p| p.do_isomeric_smiles = false)).unwrap(),
        "NC(Cc1ccccc1)C(=O)O"
    );
    assert_eq!(
        mol_to_smiles(&mut m, &p(|p| p.canonical = false)).unwrap(),
        "OC(=O)[C@@H](N)Cc1ccccc1"
    );
    assert_eq!(
        mol_to_smiles(&mut m, &p(|p| p.rooted_at_atom = 6)).unwrap(),
        "c1(C[C@H](N)C(=O)O)ccccc1"
    );
    assert_eq!(
        mol_to_smiles(&mut m, &p(|p| p.all_bonds_explicit = true)).unwrap(),
        "N-[C@@H](-C-c1:c:c:c:c:c:1)-C(=O)-O"
    );
    assert_eq!(
        mol_to_smiles(&mut m.clone(), &p(|p| p.do_kekule = true)).unwrap(),
        "N[C@@H](CC1=CC=CC=C1)C(=O)O"
    );
}

#[test]
fn fragments() {
    let mut m = parse("OC(=O)[C@@H](N)Cc1ccccc1");
    let ps = SmilesWriteParams::default();
    assert_eq!(
        mol_fragment_to_smiles(&mut m, &ps, &[0, 1, 2, 3], None, None, None).unwrap(),
        "CC(=O)O"
    );
    assert_eq!(
        mol_fragment_to_smiles(&mut m, &ps, &[0, 1, 2], Some(&[0]), None, None).unwrap(),
        "O.CO"
    );
    assert!(mol_fragment_to_smiles(&mut m, &ps, &[], None, None, None).is_err());
}
