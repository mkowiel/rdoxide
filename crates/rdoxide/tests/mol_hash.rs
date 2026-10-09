//! MolHash functions; expectations from upstream `testMolHash.py`.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::mol_hash::{HashFunction, mol_hash, split_molecule, strip, strip_type};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn mol(smiles: &str) -> RWMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn hash(smiles: &str, func: HashFunction) -> String {
    mol_hash(&mut mol(smiles), func, false, 0).unwrap()
}

#[test]
fn all_functions_on_one_molecule() {
    use HashFunction::*;
    let smi = "C1CCCC(O)C1c1ccnc(OC)c1";
    for (func, expected) in [
        (AnonymousGraph, "***1****(*2*****2*)*1"),
        (ElementGraph, "COC1CC(C2CCCCC2O)CCN1"),
        (CanonicalSmiles, "COc1cc(C2CCCCC2O)ccn1"),
        (MurckoScaffold, "c1cc(C2CCCCC2)ccn1"),
        (ExtendedMurcko, "*c1cc(C2CCCCC2*)ccn1"),
        (MolFormula, "C12H17NO2"),
        (AtomBondCounts, "15,16"),
        (DegreeVector, "0,4,9,2"),
        (Mesomer, "CO[C]1[CH][C](C2CCCCC2O)[CH][CH][N]1_0"),
        (Regioisomer, "*O.*O*.C.C1CCCCC1.c1ccncc1"),
        (NetCharge, "0"),
        (SmallWorldIndexBR, "B16R2"),
        (SmallWorldIndexBRL, "B16R2L9"),
        (ArthorSubstructureOrder, "000f001001000c000300005f000000"),
    ] {
        assert_eq!(hash(smi, func), expected, "{func:?}");
    }
}

#[test]
fn tautomer_hashes() {
    use HashFunction::*;
    assert_eq!(hash("CCC=O", HetAtomTautomer), "CC[CH][O]_0_0");
    assert_eq!(hash("CCC=O", HetAtomTautomerv2), "[CH3]-[C]:[C]:[O]_3_0");
    assert_eq!(hash("CC=CO", HetAtomProtomer), "C[CH][CH][O]_1");
    assert_eq!(hash("CC=CO", HetAtomProtomerv2), "[CH3]-[C]:[C]:[O]_3");
}

#[test]
fn formula_charges_and_empty_molecules() {
    use HashFunction::*;
    assert_eq!(hash("[NH4+]", MolFormula), "H4N+");
    assert_eq!(hash("[O-2]", MolFormula), "O-2");
    assert_eq!(hash("*C", MolFormula), "CH3X");
    assert_eq!(hash("", MolFormula), "");
    // unsigned arithmetic upstream
    assert_eq!(hash("C.C", SmallWorldIndexBR), "B0R4294967295");
}

#[test]
fn cxsmiles_hashes() {
    let m = mol(
        "C[C@@H](O)[C@@H](C)[C@@H](C)C[C@H](C1=CN=CN1)C1=CNC=N1 |o1:8,5,&1:1,3,r,c:11,18,t:9,15|",
    );
    let hash = |skip| mol_hash(&mut m.clone(), HashFunction::HetAtomTautomer, true, skip).unwrap();
    assert_eq!(
        hash(0),
        "C[C@H](CC([C]1[CH][N][CH][N]1)[C]1[CH][N][CH][N]1)[C@@H](C)[C@H](C)[O]_3_0 |o1:1,&1:14,16|"
    );
    // without the enhanced stereo block
    assert!(
        !hash(rdoxide::graphmol::smiles_parse::cx_write::cx_fields::ENHANCEDSTEREO).contains('|')
    );
}

#[test]
fn strip_and_split() {
    let mut m = mol("[13CH3:3][C@H](F)Cl.O");
    strip(
        &mut m,
        strip_type::ATOM_STEREO | strip_type::ISOTOPE | strip_type::ATOM_MAP,
    )
    .unwrap();
    assert_eq!(hash_smiles(&mut m), "CC(F)Cl.O");
    let parts = split_molecule(&m).unwrap();
    assert_eq!(parts.len(), 2);
}

fn hash_smiles(m: &mut RWMol) -> String {
    mol_hash(m, HashFunction::CanonicalSmiles, false, 0).unwrap()
}
