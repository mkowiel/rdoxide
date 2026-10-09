//! `molzip` and friends; expectations from RDKit 2026.09.1.

use std::collections::BTreeMap;

use rdoxide::graphmol::chem_transforms::{
    MolzipLabel, MolzipParams, fragment_on_bonds, molzip, molzip_fragments, molzip_row,
    molzip_single,
};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{ROMol, RWMol};

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn can(mol: &mut ROMol) -> String {
    mol_to_smiles(mol, &SmilesWriteParams::default()).unwrap()
}

fn zipped(a: &str, b: &str, params: &MolzipParams) -> String {
    can(&mut molzip(&smiles(a), &smiles(b), params).unwrap().unwrap())
}

#[test]
fn zips_atom_map_numbers_and_keeps_stereo() {
    let p = MolzipParams::default();
    assert_eq!(zipped("C=C[*:1]", "O/C=N/[*:1]", &p), "C=C/N=C/O");
    assert_eq!(zipped("C[C@H](F)[*:1]", "[*:1]Cl", &p), "C[C@H](F)Cl");
    assert_eq!(zipped("[*:1][C@H](F)Cl", "C[*:1]", &p), "C[C@H](F)Cl");
    assert_eq!(zipped("F/C=C/[*:1]", "[*:1]/C=C/Cl", &p), "F/C=C/C=C/Cl");
}

#[test]
fn zips_isotope_labels_and_link_bonds() {
    let iso = MolzipParams {
        label: MolzipLabel::Isotope,
        ..MolzipParams::default()
    };
    assert_eq!(zipped("C=C[1*]", "O/C=N/[1*]", &iso), "C=C/N=C/O");
    let combined = rdoxide::graphmol::chem_transforms::combine_mols(
        &smiles("C[*:1]"),
        &smiles("S[*:2]"),
        Default::default(),
    )
    .unwrap();
    let mut res = molzip(&smiles("[*:1][*:2]"), &combined, &MolzipParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(can(&mut res), "CS");
}

#[test]
fn zips_within_one_molecule() {
    let mol = smiles("C[*:1].N[*:1]");
    let mut res = molzip_single(&mol, &MolzipParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(can(&mut res), "CN");
}

#[test]
fn round_trips_fragment_on_bonds() {
    let p = MolzipParams {
        label: MolzipLabel::FragmentOnBonds,
        ..MolzipParams::default()
    };
    for (input, cut) in [
        ("CCOC(=O)c1ccccc1", vec![1, 3]),
        ("C[C@H](F)Cl", vec![0, 1]),
    ] {
        let mut mol = smiles(input);
        let frags = fragment_on_bonds(&mol, &cut, true, None, None, None).unwrap();
        let mut res = molzip_single(&frags, &p).unwrap().unwrap();
        assert_eq!(can(&mut res), can(&mut mol), "{input}");
    }
}

#[test]
fn zips_decompositions_and_rows() {
    let p = MolzipParams::default();
    let mut pieces = vec![smiles("[*:1]C[*:2]"), smiles("C[*:1]"), smiles("N[*:2]")];
    let mut res = molzip_fragments(&mut pieces, &p).unwrap().unwrap();
    assert_eq!(can(&mut res), "CCN");

    let mut row = BTreeMap::from([
        ("Core".to_string(), smiles("[*:1]C[*:2]")),
        ("R1".to_string(), smiles("C[*:1]")),
        ("R2".to_string(), smiles("N[*:2]")),
    ]);
    let mut res = molzip_row(&mut row, &p).unwrap().unwrap();
    assert_eq!(can(&mut res), "CCN");
    assert!(molzip_row(&mut BTreeMap::from([("R1".to_string(), smiles("C"))]), &p).is_err());
    assert!(molzip_fragments(&mut [], &p).unwrap().is_none());
}

#[test]
fn generated_coordinates_go_to_the_pieces() {
    let p = MolzipParams {
        generate_coordinates: true,
        ..MolzipParams::default()
    };
    let mut pieces = vec![smiles("[*:1]C[*:2]"), smiles("C[*:1]"), smiles("N[*:2]")];
    let res = molzip_fragments(&mut pieces, &p).unwrap().unwrap();
    assert_eq!(res.num_conformers(), 1);
    assert!(pieces.iter().all(|m| m.num_conformers() == 1));
}

#[test]
fn reports_inconsistent_labels() {
    let p = MolzipParams::default();
    let err = molzip(&smiles("C[*:1]"), &smiles("C[*:1].C[*:1]"), &p).unwrap_err();
    assert!(
        err.to_string()
            .contains("bond info already exists for end atom with label:1")
    );
    let err = molzip(&smiles("[*:1][*:2]"), &smiles("C[*:1]"), &p).unwrap_err();
    assert!(err.to_string().contains("Must have an end atom to bond"));
}
