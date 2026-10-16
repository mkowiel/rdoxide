//! `normalizeDepiction` and `straightenDepiction`; the numbers are those of
//! RDKit 2026.09.1's `testNormalizeStraighten`.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::depictor::{normalize_depiction, straighten_depiction};
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_block};

const NORADRENALINE: &str = "
  MJ201100

 12 12  0  0  1  0  0  0  0  0999 V2000
    2.2687    1.0716    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0312    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437   -0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2062    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437   -1.0716    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2687    0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  2  1  0  0  0  0
  3  4  1  6  0  0  0
  3  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  2  0  0  0  0
  9 10  1  0  0  0  0
  9 11  1  0  0  0  0
 11 12  2  0  0  0  0
  5 12  1  0  0  0  0
M  END
";

fn mol_with_copy() -> RWMol {
    let mut mol = mol_from_mol_block(NORADRENALINE, &MolFileParserParams::default())
        .unwrap()
        .unwrap();
    let copy = mol.conformer(0).unwrap().clone();
    mol.add_conformer(copy, true).unwrap();
    mol
}

fn bond_vec(mol: &RWMol, conf: i32, to: u32, from: u32) -> (f64, f64) {
    let c = mol.conformer(conf).unwrap();
    let d = c.positions()[to as usize] - c.positions()[from as usize];
    (d.x, d.y)
}

fn assert_close(got: (f64, f64), want: (f64, f64)) {
    assert!(
        (got.0 - want.0).abs() < 1e-3 && (got.1 - want.1).abs() < 1e-3,
        "{got:?} != {want:?}"
    );
}

#[test]
fn normalize_then_straighten() {
    let mut mol = mol_with_copy();
    let scale = normalize_depiction(&mut mol, 1, 1, -1.0).unwrap();
    assert!((scale - 1.875).abs() < 1e-3);
    assert_close(bond_vec(&mol, 0, 11, 10), (0.825, 0.0));
    assert_close(bond_vec(&mol, 1, 11, 10), (1.513, -0.321));
    straighten_depiction(&mut mol, 1, false).unwrap();
    assert_close(bond_vec(&mol, 1, 11, 10), (1.340, -0.773));
    assert_close(bond_vec(&mol, 1, 11, 4), (0.0, 1.547));
}

#[test]
fn normalize_to_y_axis() {
    let mut mol = mol_with_copy();
    let scale = normalize_depiction(&mut mol, 1, -1, -1.0).unwrap();
    assert!((scale - 1.875).abs() < 1e-3);
    assert_close(bond_vec(&mol, 1, 11, 10), (0.321, 1.513));
    straighten_depiction(&mut mol, 1, false).unwrap();
    assert_close(bond_vec(&mol, 1, 11, 10), (0.0, 1.547));
}

#[test]
fn explicit_scale_factor() {
    let mut mol = mol_with_copy();
    let scale = normalize_depiction(&mut mol, 1, 0, 3.0).unwrap();
    assert!((scale - 3.0).abs() < 1e-3);
    assert_close(bond_vec(&mol, 1, 11, 10), (2.475, 0.0));
    straighten_depiction(&mut mol, 1, false).unwrap();
    assert_close(bond_vec(&mol, 1, 11, 10), (2.143, -1.237));
    assert_close(bond_vec(&mol, 1, 11, 4), (0.0, 2.475));
}

#[test]
fn no_bonds_is_a_no_op() {
    let mut mol = RWMol::new();
    assert_eq!(normalize_depiction(&mut mol, -1, 1, -1.0).unwrap(), -1.0);
    straighten_depiction(&mut mol, -1, false).unwrap();
}
