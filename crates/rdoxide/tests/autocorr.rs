//! MQN and 2D autocorrelation tests; expected values from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::autocorr2d::autocorr2d;
use rdoxide::graphmol::descriptors::mqn::calc_mqns;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::RDProps;

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn mqns_of_paracetamol() {
    assert_eq!(
        calc_mqns(&smi("CC(=O)Nc1ccc(O)cc1")).unwrap(),
        [
            8, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 11, 4, 1, 0, 3, 3, 0, 1, 5, 3, 2, 2, 0, 0, 3, 1, 1, 0,
            4, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
        ]
    );
}

#[test]
fn autocorr2d_of_ethanol() {
    let mut m = smi("CCO");
    let v = autocorr2d(&mut m, "").unwrap();
    assert_eq!(v.len(), 192);
    assert_eq!(&v[..3], [1.204, 0.847, 0.0]);
    assert_eq!(&v[144..146], [0.75, 1.5]);
}

#[test]
fn autocorr2d_with_custom_property() {
    let mut m = smi("CCO");
    for i in 0..3 {
        m.atom_mut(i)
            .set_prop("w", f64::from(i) + 1.0, false)
            .unwrap();
    }
    let v = autocorr2d(&mut m, "w").unwrap();
    assert_eq!(v.len(), 32);
    // lag 1 holds the pairs (1,2) and (2,3): ln(1 + 2 + 6)
    assert_eq!(v[0], 2.197);
    assert_eq!(v[7], 0.0);
}
