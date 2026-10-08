//! Gasteiger charge tests; expected values from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::partial_charges::gasteiger::{
    DEFAULT_ITERATIONS, compute_gasteiger_charges,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::common_properties as cp;
use rdoxide::rdgeneral::dict::{RDProps, RDValue};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-11, "{a} != {b}");
}

#[test]
fn acetic_acid() {
    let mut m = smi("CC(=O)O");
    let charges = compute_gasteiger_charges(&mut m, DEFAULT_ITERATIONS, false).unwrap();
    for (c, expected) in charges.iter().zip([
        0.033767593407,
        0.299684553158,
        -0.252820433497,
        -0.481433475172,
    ]) {
        close(*c, expected);
    }
    // the charges are stored as computed properties, with the implicit-H charge
    assert!(matches!(
        m.atom(3).get_prop(cp::_GasteigerHCharge),
        Some(RDValue::Double(h)) if (h - 0.296256680872).abs() < 1e-11
    ));
}

#[test]
fn zero_iterations_leave_zero_charges() {
    let charges = compute_gasteiger_charges(&mut smi("CC(=O)O"), 0, false).unwrap();
    assert_eq!(charges, vec![0.0; 4]);
}

#[test]
fn conjugated_formal_charge_is_shared() {
    let charges = compute_gasteiger_charges(&mut smi("NC(=[NH2+])c1ccccc1"), 0, false).unwrap();
    close(charges[0], 0.5);
    close(charges[2], 0.5);
}

#[test]
fn missing_parameters_can_raise() {
    let err = compute_gasteiger_charges(&mut smi("O=[Se]=O"), 12, true).unwrap_err();
    assert!(err.to_string().contains("Element: Se Mode: sp2"));
}
