//! `ShapeHelpers/testShapeHelpers.cpp`.

use rdoxide::geometry::grid::UniformGrid3D;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};
use rdoxide::graphmol::mol_transforms::canonicalize_mol;
use rdoxide::graphmol::shape_helpers::{ShapeGridOptions, encode_shape, mol};

const DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/rdkit/Code/GraphMol/ShapeHelpers/test_data/"
);

fn load(name: &str) -> Option<rdoxide::graphmol::RWMol> {
    if !std::path::Path::new(DIR).exists() {
        return None;
    }
    mol_from_mol_file(&format!("{DIR}{name}"), &MolFileParserParams::default()).ok()?
}

#[test]
fn encode() {
    let Some(mut m) = load("1oir.mol") else {
        return;
    };
    canonicalize_mol(&mut m, false, true).unwrap();
    let mut grd = UniformGrid3D::new(30.0, 16.0, 10.0, 0.5, 1, None).unwrap();
    encode_shape(
        &m,
        m.conformer(0).unwrap(),
        &mut grd,
        None,
        0.8,
        0.25,
        -1,
        true,
    )
    .unwrap();
    assert_eq!(grd.occupancy_vect().total_val(), 7405);
}

#[test]
fn compare() {
    let Some(mut m) = load("1oir.mol") else {
        return;
    };
    canonicalize_mol(&mut m, false, true).unwrap();
    let mut dup = load("1oir.mol").unwrap();
    canonicalize_mol(&mut dup, false, true).unwrap();
    let o = ShapeGridOptions::default();
    let d = mol::tanimoto_distance(&m, &dup, -1, -1, &o).unwrap();
    assert!(d.abs() < 1e-8);
    let t = mol::tversky_index(&m, &dup, 1.0, 1.0, -1, -1, &o).unwrap();
    assert!((t - 1.0).abs() < 1e-8);

    // GitHub #4364: the grid spacing is used
    let mut m2 = load("1oir_conf.mol").unwrap();
    canonicalize_mol(&mut m2, false, true).unwrap();
    let d1 = mol::tanimoto_distance(&m, &m2, -1, -1, &o).unwrap();
    assert!((d1 - 0.31).abs() < 0.01, "{d1}");
    let o1 = ShapeGridOptions {
        grid_spacing: 1.0,
        ..o
    };
    let d2 = mol::tanimoto_distance(&m, &m2, -1, -1, &o1).unwrap();
    assert!((d2 - d1).abs() > 0.001);
}
