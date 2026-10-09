//! MolTransforms tests; expected values from RDKit 2026.09.1.

use std::path::PathBuf;

use rdoxide::RDKitError;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};
use rdoxide::graphmol::mol_transforms::*;

fn load(name: &str) -> RWMol {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data/molfiles")
        .join(name);
    let params = MolFileParserParams {
        remove_hs: false,
        ..MolFileParserParams::default()
    };
    mol_from_mol_file(path.to_str().unwrap(), &params)
        .unwrap()
        .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}

#[test]
fn principal_moments_of_ethanol() {
    let mut mol = load("ethanol_3d.mol");
    let res = compute_principal_axes_and_moments(&mut mol, -1, true, true, None)
        .unwrap()
        .unwrap();
    close(res.moments[0], 0.4960649635983768);
    close(res.moments[1], 2.7251972097349584);
    close(res.moments[2], 3.2212621733333333);
    close(res.axes[0][0], 0.9456520915724027);
    // unweighted results are cached on the molecule and reused
    let again = compute_principal_axes_and_moments(&mut mol, -1, true, false, None)
        .unwrap()
        .unwrap();
    assert_eq!(res, again);
}

#[test]
fn canonical_transform_aligns_the_principal_axes() {
    let mol = load("ethanol_3d.mol");
    let (trans, _) =
        compute_canonical_transform(&mol, mol.conformer(-1).unwrap(), None, false, true).unwrap();
    close(trans.get_val(0, 0), 0.9456520915724027);
    close(trans.get_val(0, 3), -0.3043494585647698);
    close(trans.get_val(3, 3), 1.0);
}

#[test]
fn internal_coordinates_round_trip() {
    let mut mol = load("ethanol_3d.mol");
    set_angle_deg(&mut mol, -1, 0, 1, 2, 100.0).unwrap();
    close(
        get_angle_deg(mol.conformer(-1).unwrap(), 0, 1, 2).unwrap(),
        100.0,
    );
    set_dihedral_deg(&mut mol, -1, 0, 1, 2, 8, 60.0).unwrap();
    close(
        get_dihedral_deg(mol.conformer(-1).unwrap(), 0, 1, 2, 8).unwrap(),
        60.0,
    );
    set_bond_length(&mut mol, -1, 0, 1, 2.0).unwrap();
    close(
        get_bond_length(mol.conformer(-1).unwrap(), 0, 1).unwrap(),
        2.0,
    );
}

#[test]
fn ring_bonds_cannot_be_changed() {
    let mut mol = load("cyclohexane_3d.mol");
    assert!(matches!(
        set_bond_length(&mut mol, -1, 0, 1, 2.0),
        Err(RDKitError::Value(m)) if m.contains("must not belong to a ring")
    ));
    assert!(matches!(
        get_bond_length(mol.conformer(-1).unwrap(), 0, 99),
        Err(RDKitError::Invariant { .. })
    ));
}
