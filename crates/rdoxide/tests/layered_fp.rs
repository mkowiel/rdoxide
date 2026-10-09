//! Layered fingerprint tests; expected values from the reference RDKit.

use rdoxide::datastructs::{BitVect, ExplicitBitVect};
use rdoxide::graphmol::fingerprints::layered::{
    LayeredFingerprintArgs, SUBSTRUCT_LAYERS, layered_fingerprint_mol,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn run(smiles: &str, args: &LayeredFingerprintArgs, counts: Option<&mut [u32]>) -> Vec<u32> {
    let mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let fp: ExplicitBitVect = layered_fingerprint_mol(&mol, args, counts, None).unwrap();
    (0..fp.num_bits())
        .filter(|&b| fp.get_bit(b).unwrap())
        .collect()
}

#[test]
fn ethanol_bits_and_atom_counts() {
    let args = LayeredFingerprintArgs {
        fp_size: 256,
        ..Default::default()
    };
    let mut counts = [0u32; 3];
    assert_eq!(
        run("CCO", &args, Some(&mut counts)),
        [20, 84, 87, 92, 98, 99, 104, 162, 247, 248]
    );
    assert_eq!(counts, [2, 3, 2]);
}

#[test]
fn substructure_layers_of_cyclopropane() {
    let args = LayeredFingerprintArgs {
        layer_flags: SUBSTRUCT_LAYERS,
        fp_size: 128,
        ..Default::default()
    };
    assert_eq!(
        run("C1CC1", &args, None),
        [20, 34, 54, 97, 99, 104, 112, 118, 119]
    );
}

#[test]
fn bad_arguments_are_errors() {
    let mol = mol_from_smiles("CCO", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let bad = |f: fn(&mut LayeredFingerprintArgs)| {
        let mut args = LayeredFingerprintArgs::default();
        f(&mut args);
        layered_fingerprint_mol(&mol, &args, None, None).is_err()
    };
    assert!(bad(|a| a.min_path = 0));
    assert!(bad(|a| a.max_path = 0));
    assert!(bad(|a| a.fp_size = 0));
    let mut short = [0u32; 1];
    assert!(layered_fingerprint_mol(&mol, &Default::default(), Some(&mut short), None).is_err());
}
