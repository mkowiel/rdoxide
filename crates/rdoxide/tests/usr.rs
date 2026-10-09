//! USR, USRCAT and Coulomb matrix tests; expected values from RDKit 2026.09.1.

use std::path::PathBuf;

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::coulomb_mat::coulomb_mat;
use rdoxide::graphmol::descriptors::usr::*;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};

fn mol(file: &str) -> RWMol {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data/molfiles")
        .join(file);
    let params = MolFileParserParams {
        remove_hs: false,
        ..MolFileParserParams::default()
    };
    mol_from_mol_file(path.to_str().unwrap(), &params)
        .unwrap()
        .unwrap()
}

fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() <= 1e-12 * y.abs().max(1.0), "{a:?} != {b:?}");
    }
}

#[test]
fn coulomb_matrix_of_ethanol() {
    let m = coulomb_mat(&mol("ethanol_3d.mol"), -1).unwrap();
    assert_eq!(m.len(), 9);
    close(
        &m[0][..3],
        &[36.85810519942594, 23.63585207528015, 20.584422909870824],
    );
    assert_eq!(m[1][2], m[2][1]);
}

#[test]
fn usr_of_ethanol_has_twelve_values() {
    let m = mol("ethanol_3d.mol");
    assert_eq!(usr(&m, -1).unwrap().len(), 12);
    assert_eq!(usrcat(&m, None, -1).unwrap().len(), 60);
    let sel = usrcat(&m, Some(&[vec![0, 1], vec![]]), -1).unwrap();
    assert_eq!(sel.len(), 36);
    assert_eq!(&sel[24..], &[0.0; 12]);
    assert!(usrcat(&m, Some(&[vec![99]]), -1).is_err());
}

#[test]
fn distributions_and_score() {
    let pts = [
        Point3D::new(0.0, 0.0, 0.0),
        Point3D::new(1.0, 0.0, 0.0),
        Point3D::new(0.0, 2.0, 0.0),
        Point3D::new(0.0, 0.0, 3.0),
    ];
    let (dist, refs) = calc_usr_distributions(&pts).unwrap();
    assert_eq!(
        calc_usr_distributions_from_points(&pts, &refs),
        dist.to_vec()
    );
    let d = calc_usr_from_distributions(&dist);
    assert_eq!(d.len(), 12);
    assert_eq!(calc_usr_score(&d, &d, &[1.0]).unwrap(), 1.0);
    assert!(calc_usr_score(&d, &d, &[]).is_err());
    assert!(calc_usr_distributions(&[]).is_err());
}
