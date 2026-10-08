//! BCUT2D tests; expected values from RDKit 2026.09.1 (eigenvalues agree to 1e-9).

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::bcut::{bcut2d, bcut2d_with_prop_name, bcut2d_with_props};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() < 1e-9, "{a:?} != {b:?}");
    }
}

#[test]
fn ethanol() {
    close(
        &bcut2d(&smi("CCO")).unwrap(),
        &[
            16.248991329,
            10.908168035,
            1.340293519,
            -1.522358693,
            1.296242638,
            -1.568363547,
            3.835802336,
            0.318405613,
        ],
    );
}

#[test]
fn explicit_properties() {
    let m = smi("CCO");
    let (hi, lo) = bcut2d_with_props(&m, &[1.0, 2.0, 3.0]).unwrap();
    close(&[hi, lo], &[3.7323843332905975, 0.26828233325261325]);
    assert!(bcut2d_with_props(&m, &[1.0]).is_err());
    assert!(bcut2d_with_prop_name(&m, "missing").is_err());
}

#[test]
fn empty_molecule() {
    assert_eq!(bcut2d(&smi("")).unwrap(), vec![0.0; 8]);
    assert_eq!(bcut2d_with_props(&smi(""), &[]).unwrap(), (0.0, 0.0));
}
