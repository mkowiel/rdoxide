//! Depictions constrained to a reference structure.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::depictor::constrained::{
    ConstrainedDepictionParams, generate_depiction_matching_2d_structure,
    generate_depiction_matching_2d_structure_with_map,
};
use rdoxide::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn depicted(s: &str) -> RWMol {
    let mut mol = smiles(s);
    compute_2d_coords(&mut mol, &Compute2DCoordParameters::default()).unwrap();
    mol
}

fn pos(mol: &RWMol, idx: usize) -> (f64, f64) {
    let p = mol.conformer(-1).unwrap().positions()[idx];
    (p.x, p.y)
}

fn assert_core_matches(mol: &RWMol, reference: &RWMol, map: &[(i32, i32)], tol: f64) {
    for &(r, m) in map {
        let (a, b) = (pos(reference, r as usize), pos(mol, m as usize));
        assert!(
            (a.0 - b.0).abs() < tol && (a.1 - b.1).abs() < tol,
            "ref {r} {a:?} != mol {m} {b:?}"
        );
    }
}

#[test]
fn hard_constraint_reproduces_the_reference_coordinates() {
    let reference = depicted("c1ccccc1C");
    let mut mol = smiles("CCOc1ccc(C)cc1");
    let map = generate_depiction_matching_2d_structure(
        &mut mol,
        &reference,
        -1,
        None,
        &ConstrainedDepictionParams::default(),
    );
    // the reference is a substructure of the molecule only through its ring
    // and methyl, which are there
    let map = map.unwrap();
    assert_eq!(map.len(), 7);
    assert_core_matches(&mol, &reference, &map, 1e-4);
}

#[test]
fn align_only_keeps_the_internal_geometry() {
    let reference = depicted("c1ccccc1C");
    let mut mol = depicted("CCOc1ccc(C)cc1");
    let before: Vec<(f64, f64)> = (0..mol.num_atoms() as usize)
        .map(|i| pos(&mol, i))
        .collect();
    let params = ConstrainedDepictionParams {
        align_only: true,
        ..Default::default()
    };
    let map =
        generate_depiction_matching_2d_structure(&mut mol, &reference, -1, None, &params).unwrap();
    assert!(!map.is_empty());
    assert_eq!(mol.num_conformers(), 1);
    // a rigid-body motion preserves distances
    let d = |p: (f64, f64), q: (f64, f64)| ((p.0 - q.0).powi(2) + (p.1 - q.1).powi(2)).sqrt();
    for i in 1..before.len() {
        let (a, b) = (d(before[0], before[i]), d(pos(&mol, 0), pos(&mol, i)));
        assert!((a - b).abs() < 1e-6);
    }
}

#[test]
fn no_match_is_an_error_unless_failure_is_accepted() {
    let reference = depicted("c1ccccc1");
    let mut mol = smiles("CCCC");
    let err = generate_depiction_matching_2d_structure(
        &mut mol,
        &reference,
        -1,
        None,
        &ConstrainedDepictionParams::default(),
    )
    .unwrap_err();
    assert!(
        err.to_string()
            .contains("Substructure match with reference not found.")
    );
    let params = ConstrainedDepictionParams {
        accept_failure: true,
        ..Default::default()
    };
    let map =
        generate_depiction_matching_2d_structure(&mut mol, &reference, -1, None, &params).unwrap();
    assert!(map.is_empty());
    assert_eq!(mol.num_conformers(), 1);
}

#[test]
fn explicit_atom_map_and_pattern() {
    let reference = depicted("c1ccccc1C");
    let mut mol = smiles("Cc1ccccc1");
    // reference atom -> molecule atom for the ring
    let map: Vec<(i32, i32)> = (0..6).map(|i| (i, i + 1)).collect();
    generate_depiction_matching_2d_structure_with_map(
        &mut mol,
        &reference,
        &map,
        -1,
        &ConstrainedDepictionParams::default(),
    )
    .unwrap();
    assert_core_matches(&mol, &reference, &map, 1e-4);

    let pattern = mol_from_smarts("c1ccccc1", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    let mut mol = smiles("OCc1ccccc1");
    let used = generate_depiction_matching_2d_structure(
        &mut mol,
        &reference,
        -1,
        Some(&pattern),
        &ConstrainedDepictionParams::default(),
    )
    .unwrap();
    assert_eq!(used.len(), 6);
    assert_core_matches(&mol, &reference, &used, 1e-4);
}

#[test]
fn out_of_range_atom_map_is_rejected() {
    let reference = depicted("CC");
    let mut mol = smiles("CC");
    let too_long = vec![(0, 0), (1, 1), (0, 1)];
    assert!(
        generate_depiction_matching_2d_structure_with_map(
            &mut mol,
            &reference,
            &too_long,
            -1,
            &ConstrainedDepictionParams::default()
        )
        .is_err()
    );
}
