//! CIP labeling; cases from upstream `CIPLabeler/catch_tests.cpp`.

use rdoxide::graphmol::chirality::legacy::assign_stereochemistry;
use rdoxide::graphmol::cip_labeler::{
    assign_cip_labels, assign_cip_labels_selected, is_max_iterations_exceeded_error,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{BondStereo, RWMol};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn parse(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn atom_cip(m: &RWMol, idx: u32) -> Option<String> {
    match m.atom(idx).get_prop("_CIPCode") {
        Some(RDValue::String(s)) => Some(s.clone()),
        _ => None,
    }
}

fn bond_cip(m: &RWMol, idx: u32) -> Option<String> {
    match m.bond(idx).get_prop("_CIPCode") {
        Some(RDValue::String(s)) => Some(s.clone()),
        _ => None,
    }
}

/// `view_labels`
fn view_labels(m: &RWMol) -> String {
    let mut out = String::new();
    for i in 0..m.num_atoms() {
        if let Some(l) = atom_cip(m, i) {
            out += &format!("{i}{l} ");
        }
    }
    for i in 0..m.num_bonds() {
        if let Some(l) = bond_cip(m, i) {
            let b = m.bond(i);
            out += &format!("{}={}{l} ", b.begin_atom_idx(), b.end_atom_idx());
        }
    }
    out
}

fn check_atoms(smiles: &str, expected: &[(u32, &str)]) {
    let mut m = parse(smiles);
    assign_cip_labels(&mut m, 0).unwrap();
    for &(idx, label) in expected {
        assert_eq!(
            atom_cip(&m, idx).as_deref(),
            Some(label),
            "{smiles} atom {idx}"
        );
    }
}

#[test]
fn tetrahedral_assignment() {
    let mut m = parse("Br[C@H](Cl)F");
    m.atom_mut(1).clear_prop("_CIPCode");
    assert!(!m.has_prop("_CIPComputed"));
    assign_cip_labels(&mut m, 0).unwrap();
    assert!(m.has_prop("_CIPComputed"));
    assert_eq!(atom_cip(&m, 1).as_deref(), Some("S"));
}

#[test]
fn double_bond_assignment() {
    let mut m = parse(r"CC\C(\C(\C)=N\O)=N\O");
    assert_eq!(m.bond(4).stereo(), BondStereo::StereoE);
    assert_eq!(m.bond(6).stereo(), BondStereo::StereoZ);
    assign_cip_labels(&mut m, 0).unwrap();
    assert_eq!(bond_cip(&m, 4).as_deref(), Some("E"));
    assert_eq!(bond_cip(&m, 6).as_deref(), Some("Z"));
}

#[test]
fn phosphine_and_arsine() {
    for (smiles, label) in [
        ("C[P@](C1CCCC1)C1=CC=CC=C1", "R"),
        ("C[As@@](C1CCCC1)C1=CC=CC=C1", "S"),
        ("C[P@H]C1CCCCC1", "R"),
        ("C[P@@H]C1CCCCC1", "S"),
    ] {
        check_atoms(smiles, &[(1, label)]);
    }
}

#[test]
fn para_stereochemistry() {
    check_atoms(
        r"C\C=C/[C@@H](\C=C\O)[C@H](C)[C@H](\C=C/C)\C=C\O",
        &[(3, "R"), (7, "r"), (9, "S")],
    );
    check_atoms(
        r"C\C=C/[C@@H](\C=C\C)[C@H](C)[C@H](\C=C/C)\C=C\C",
        &[(3, "R"), (7, "r"), (9, "S")],
    );
    // the label of atom 12 depends on atom 2, which is not reached in the
    // initial expansion of atom 9
    check_atoms(
        "CC[C@H](C)CCCCC[C@H]1CC[C@@H](C)CC1",
        &[(2, "S"), (9, "s"), (12, "S")],
    );
}

#[test]
fn dummy_atoms() {
    let mut m = parse("*[C@](F)(Cl)Br");
    assign_stereochemistry(&mut m, true, true, false).unwrap();
    assert_eq!(atom_cip(&m, 1).as_deref(), Some("S"));
    m.atom_mut(1).clear_prop("_CIPCode");
    assign_cip_labels(&mut m, 0).unwrap();
    assert_eq!(atom_cip(&m, 1).as_deref(), Some("S"));

    let mut m = parse("*[C@](*)(Cl)Br");
    assign_cip_labels(&mut m, 0).unwrap();
    assert_eq!(atom_cip(&m, 1), None);

    let mut m = parse("C1CC[C@](*)2CCCC[C@H]2C1");
    assign_cip_labels(&mut m, 0).unwrap();
    assert_eq!(atom_cip(&m, 3).as_deref(), Some("s"));
    assert_eq!(atom_cip(&m, 9).as_deref(), Some("s"));
}

#[test]
fn github_5142_does_not_fail() {
    let mut m = parse("*C1C[C@H](CCC)[C@@H](C)[C@H](C)C1");
    assign_cip_labels(&mut m, 0).unwrap();
}

#[test]
fn bad_skip_side() {
    let mut m = parse(r"C1CC[C@H]2C/C(=C3\C[C@H]4CCCC[C@H]4C3)C[C@H]2C1");
    assign_cip_labels(&mut m, 0).unwrap();
    assert_eq!(view_labels(&m), "3S 8R 13S 16R 5=6E ");
    let mut m = parse(r"C1/C(=C2\C[C@H]3CCCC[C@H]3C2)C[C@H]2CCCC[C@@H]12");
    assign_cip_labels(&mut m, 0).unwrap();
    assert_eq!(view_labels(&m), "4R 9S 12R 17S 1=2E ");
}

#[test]
fn iteration_limit_includes_the_preliminary_pass() {
    let mut m = parse("C[C@H](F)Cl");
    let e = assign_cip_labels(&mut m, 1).unwrap_err();
    assert!(is_max_iterations_exceeded_error(&e));
    // the budget of the failed call must not leak into the next one
    let mut m = parse("C[C@H](F)Cl");
    assign_cip_labels(&mut m, 0).unwrap();
    assert_eq!(atom_cip(&m, 1).as_deref(), Some("R"));
}

#[test]
fn selected_atoms_and_bonds() {
    let smiles = r"C[C@H](F)C(Cl)=C/C=C/C";
    let mut m = parse(smiles);
    let n_atoms = m.num_atoms() as usize;
    let n_bonds = m.num_bonds() as usize;
    let mut atoms = vec![false; n_atoms];
    atoms[1] = true;
    assign_cip_labels_selected(&mut m, &atoms, &vec![false; n_bonds], 0).unwrap();
    assert!(atom_cip(&m, 1).is_some());
    assert!(!m.has_prop("_CIPComputed"));
    // the size of the selections is checked
    assert!(assign_cip_labels_selected(&mut m, &[true], &vec![true; n_bonds], 0).is_err());
    assert!(assign_cip_labels_selected(&mut m, &vec![true; n_atoms], &[true], 0).is_err());
}
