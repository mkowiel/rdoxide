//! AUTOCORR3D, RDF and MORSE tests; expected values from RDKit 2026.09.1.

use std::path::PathBuf;

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::autocorr3d::autocorr3d;
use rdoxide::graphmol::descriptors::morse::morse;
use rdoxide::graphmol::descriptors::rdf::rdf;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::dict::RDProps;

fn ethanol() -> RWMol {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/molfiles/ethanol_3d.mol");
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
        assert!((x - y).abs() < 1e-9, "{a:?} != {b:?}");
    }
}

#[test]
fn ethanol_profiles() {
    let mut m = ethanol();
    let a = autocorr3d(&mut m, -1, "").unwrap();
    assert_eq!(a.len(), 80);
    close(&a[..6], &[0.262, 0.737, 0.905, 0.287, 0.0, 0.0]);
    let r = rdf(&mut m, -1, "").unwrap();
    assert_eq!(r.len(), 210);
    close(&r[..6], &[2.392, 1.19, 2.534, 0.103, 4.522, 0.432]);
    let o = morse(&mut m, -1, "").unwrap();
    assert_eq!(o.len(), 224);
    close(&o[..6], &[36.0, 13.636, -1.214, -0.275, -0.53, -1.042]);
    assert!(m.has_prop("_3DDistanceMatrix_Conf0"));
}

#[test]
fn custom_property_sizes() {
    let mut m = ethanol();
    assert_eq!(autocorr3d(&mut m, -1, "w").unwrap().len(), 10);
    assert_eq!(rdf(&mut m, -1, "w").unwrap().len(), 30);
    assert_eq!(morse(&mut m, -1, "w").unwrap().len(), 32);
}

#[test]
fn errors() {
    let mut flat = mol_from_smiles("CCO", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assert!(autocorr3d(&mut flat, -1, "").is_err());
    assert!(rdf(&mut flat, -1, "").is_err());
    assert!(morse(&mut flat, -1, "").is_err());
    let mut m = ethanol();
    assert!(rdf(&mut m, 5, "").is_err());
}
