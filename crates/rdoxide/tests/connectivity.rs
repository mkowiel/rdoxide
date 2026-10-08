//! Chi / Kappa tests; expected values from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::connectivity::{
    calc_chi_nn, calc_chi_nv, calc_chi0v, calc_chi1n, calc_hall_kier_alpha, calc_kappa1,
    calc_kappa2, calc_kappa3,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
}

#[test]
fn benzene() {
    let mut m = smi("c1ccccc1");
    close(calc_chi0v(&mut m, false).unwrap(), 3.464101615137755);
    close(calc_chi1n(&mut m, false).unwrap(), 2.0000000000000004);
    close(calc_chi_nv(&mut m, 3, false).unwrap(), 0.6666666666666671);
    close(calc_hall_kier_alpha(&m, None).unwrap(), -0.78);
    close(calc_kappa1(&m).unwrap(), 3.4115708812260532);
    close(calc_kappa2(&m).unwrap(), 1.6057694396735218);
    close(calc_kappa3(&m).unwrap(), 0.5823992601400448);
}

#[test]
fn rings_count_their_first_atom_once() {
    // cyclopropane: the 3-atom ring paths close on their first atom
    let mut m = smi("C1CC1");
    close(calc_chi_nv(&mut m, 2, false).unwrap(), 1.060660171779821);
    close(calc_chi_nn(&mut m, 3, false).unwrap(), 0.3535533905932737);
}

#[test]
fn hall_kier_alpha_fills_the_contributions() {
    let m = smi("CC(=O)Nc1ccccc1");
    let mut contribs = vec![0.0; m.num_atoms() as usize];
    let alpha = calc_hall_kier_alpha(&m, Some(&mut contribs)).unwrap();
    close(alpha, contribs.iter().sum());
    assert!(calc_hall_kier_alpha(&m, Some(&mut [0.0])).is_err());
}

#[test]
fn zeroth_order_variants() {
    let mut m = smi("CCO");
    close(calc_chi0v(&mut m, false).unwrap(), 2.1543203766865053);
    // the general functions count the atom paths of length one instead
    close(calc_chi_nv(&mut m, 0, false).unwrap(), 3.0);
    close(calc_chi_nn(&mut m, 0, false).unwrap(), 3.0);
}
