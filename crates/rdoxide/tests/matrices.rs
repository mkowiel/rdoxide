//! Distance and adjacency matrix tests; expected values from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::matrices::{
    get_adjacency_matrix, get_distance_mat, get_distance_mat_cached, get_distance_mat_weighted,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::dict::RDProps;

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn bond_order_weights_match_floyd_warshall() {
    let m = smi("C=CC#N");
    let d = get_distance_mat_weighted(&m, true, false).unwrap();
    assert_eq!(&d[..4], &[0.0, 0.5, 1.5, 1.8333333333333333]);
    // without bond orders the weighted variant is the plain BFS matrix
    assert_eq!(
        get_distance_mat_weighted(&m, false, false).unwrap(),
        get_distance_mat(&m)
    );
}

#[test]
fn atom_weights_fill_the_diagonal() {
    let m = smi("CO");
    let d = get_distance_mat_weighted(&m, false, true).unwrap();
    assert_eq!(d, vec![1.0, 1.0, 1.0, 0.75]);
}

#[test]
fn adjacency_matrix_with_bond_orders() {
    let m = smi("C=CC#N");
    let a = get_adjacency_matrix(&m, true, 0);
    assert_eq!(&a[4..8], &[2.0, 0.0, 1.0, 0.0]);
    assert_eq!(
        &get_adjacency_matrix(&m, false, 0)[..4],
        &[0.0, 1.0, 0.0, 0.0]
    );
}

#[test]
fn cached_matrix_is_reused_until_forced() {
    let mut m = smi("CCO");
    let first = get_distance_mat_cached(&mut m, false, false, false, "").unwrap();
    assert!(m.has_prop("DistanceMatrix"));
    assert!(m.has_prop("DistanceMatrix_Paths"));
    assert_eq!(
        get_distance_mat_cached(&mut m, false, false, false, "").unwrap(),
        first
    );
    m.clear_computed_props(true);
    assert!(!m.has_prop("DistanceMatrix"));
}
