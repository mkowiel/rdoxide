//! `needs_hs`, molfile parities and `set_terminal_atom_coords` error paths;
//! expectations from RDKit 2026.09.1.

use rdoxide::graphmol::atom::{ChiralType, set_atom_rlabel};
use rdoxide::graphmol::chirality::from_coords::assign_chiral_types_from_mol_parity;
use rdoxide::graphmol::mol_ops::add_hs::{needs_hs, set_terminal_atom_coords};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{Atom, RWMol};
use rdoxide::rdgeneral::RDProps;
use rdoxide::rdgeneral::common_properties as cp;

fn smiles(s: &str) -> RWMol {
    let params = SmilesParserParams {
        remove_hs: false,
        ..SmilesParserParams::default()
    };
    mol_from_smiles(s, &params).unwrap().unwrap()
}

#[test]
fn needs_hs_follows_implicit_hydrogens() {
    assert!(needs_hs(&smiles("CC")).unwrap());
    assert!(!needs_hs(&smiles("[H]C([H])([H])[H]")).unwrap());
    assert!(!needs_hs(&smiles("[Na+].[Cl-]")).unwrap());
}

#[test]
fn mol_parity_sets_tetrahedral_tags() {
    let mut mol = smiles("C(F)(Cl)Br");
    for (parity, tag) in [
        (1, ChiralType::TetrahedralCW),
        (2, ChiralType::TetrahedralCCW),
        (0, ChiralType::Unspecified),
    ] {
        mol.atom_mut(0)
            .set_prop(cp::molParity, parity, false)
            .unwrap();
        assign_chiral_types_from_mol_parity(&mut mol, true).unwrap();
        assert_eq!(mol.atom(0).chiral_tag(), tag, "parity {parity}");
    }
}

#[test]
fn mol_parity_keeps_existing_tags_on_request() {
    let mut mol = smiles("C(F)(Cl)Br");
    mol.atom_mut(0).set_chiral_tag(ChiralType::TetrahedralCW);
    mol.atom_mut(0).set_prop(cp::molParity, 2, false).unwrap();
    assign_chiral_types_from_mol_parity(&mut mol, false).unwrap();
    assert_eq!(mol.atom(0).chiral_tag(), ChiralType::TetrahedralCW);
    // atoms with fewer than three neighbors are reset
    assign_chiral_types_from_mol_parity(&mut mol, true).unwrap();
    assert_eq!(mol.atom(0).chiral_tag(), ChiralType::TetrahedralCCW);
    mol.atom_mut(1).set_prop(cp::molParity, 1, false).unwrap();
    mol.atom_mut(1).set_chiral_tag(ChiralType::TetrahedralCW);
    assign_chiral_types_from_mol_parity(&mut mol, true).unwrap();
    assert_eq!(mol.atom(1).chiral_tag(), ChiralType::Unspecified);
}

#[test]
fn set_terminal_atom_coords_checks_its_arguments() {
    let mut mol = smiles("C(F)(Cl)Br");
    let msg = |r: rdoxide::Result<()>| r.unwrap_err().to_string();
    assert!(msg(set_terminal_atom_coords(&mut mol, 1, 1)).contains("degenerate atoms"));
    assert!(msg(set_terminal_atom_coords(&mut mol, 100, 1)).contains("Range Error"));
    assert!(msg(set_terminal_atom_coords(&mut mol, 0, 1)).contains("bad atom degree"));
    assert!(msg(set_terminal_atom_coords(&mut mol, 1, 2)).contains("no bond between atoms"));
    assert!(msg(set_terminal_atom_coords(&mut mol, 1, 50)).contains("idx1"));
}

#[test]
fn rlabel_must_be_a_two_digit_number() {
    let mut atom = Atom::new(6);
    set_atom_rlabel(&mut atom, 99).unwrap();
    let err = set_atom_rlabel(&mut atom, 100).unwrap_err().to_string();
    assert!(err.contains("rlabel out of range for MDL files"), "{err}");
    assert!(
        err.contains("line 873 in file Code/GraphMol/Atom.cpp"),
        "{err}"
    );
}
