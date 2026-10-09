//! PMI and PBF tests; expected values from RDKit 2026.09.1.

use std::path::PathBuf;

use rdoxide::RDKitError;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::descriptors::pbf::pbf;
use rdoxide::graphmol::descriptors::pmi::*;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_file};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

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

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn ethanol_descriptors() {
    let mut m = ethanol();
    close(pbf(&mut m, -1).unwrap(), 0.4823414930196825);
    close(npr1(&mut m, -1, true, true).unwrap(), 0.2540480685156133);
    close(npr2(&mut m, -1, true, true).unwrap(), 0.8728777586729661);
    close(pmi1(&mut m, -1, true, true).unwrap(), 15.11773660225969);
    close(pmi3(&mut m, -1, true, true).unwrap(), 59.507386498120866);
    close(
        radius_of_gyration(&mut m, -1, true, true).unwrap(),
        1.1720394242631202,
    );
    close(
        eccentricity(&mut m, -1, true, true).unwrap(),
        0.9671915936790839,
    );
    close(
        asphericity(&mut m, -1, true, true).unwrap(),
        0.422455319681518,
    );
    close(
        spherocity_index(&mut m, -1, true).unwrap(),
        0.4607287308599609,
    );
}

#[test]
fn molecules_without_conformers_are_rejected() {
    let mut m = mol_from_smiles("CCO", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assert!(matches!(
        pmi1(&mut m, -1, true, true),
        Err(RDKitError::Invariant { .. })
    ));
    assert!(matches!(pbf(&mut m, -1), Err(RDKitError::Invariant { .. })));
}

#[test]
fn missing_conformer_ids_are_conformer_errors() {
    let mut m = ethanol();
    assert!(matches!(
        pmi1(&mut m, 5, true, true),
        Err(RDKitError::Conformer(_))
    ));
}
