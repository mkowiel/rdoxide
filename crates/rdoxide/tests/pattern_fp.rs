//! Pattern fingerprint tests; expected bits from the reference RDKit.

use rdoxide::datastructs::BitVect;
use rdoxide::graphmol::fingerprints::pattern::pattern_fingerprint_mol;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn on_bits(smiles: &str, fp_size: u32, tautomeric: bool) -> Vec<u32> {
    let mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let fp = pattern_fingerprint_mol(&mol, fp_size, tautomeric).unwrap();
    (0..fp.num_bits())
        .filter(|&b| fp.get_bit(b).unwrap())
        .collect()
}

#[test]
fn ethanol() {
    assert_eq!(
        on_bits("CCO", 512, false),
        [
            9, 37, 54, 110, 173, 212, 217, 247, 265, 266, 271, 362, 429, 434, 455, 510
        ]
    );
}

#[test]
fn tautomeric_fingerprint_adds_bits_for_tautomerizable_bonds() {
    let plain = on_bits("O=C1NC=CC=C1", 512, false);
    let taut = on_bits("O=C1NC=CC=C1", 512, true);
    assert_ne!(plain, taut);
    assert!(taut.contains(&16) && taut.contains(&45));
}

#[test]
fn zero_size_is_an_error() {
    let mol = mol_from_smiles("C", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assert!(pattern_fingerprint_mol(&mol, 0, false).is_err());
}
