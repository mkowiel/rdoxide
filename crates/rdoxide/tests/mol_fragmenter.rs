//! Bond fragmentation and the BRICS definitions; expectations from RDKit
//! 2026.09.1.

use rdoxide::graphmol::bond::BondType;
use rdoxide::graphmol::chem_transforms::mol_fragmenter::{
    construct_brics_atom_types, construct_brics_bond_types, construct_fragmenter_atom_types,
    construct_fragmenter_bond_types,
};
use rdoxide::graphmol::chem_transforms::{
    fragment_on_bonds, fragment_on_brics_bonds, fragment_on_some_bonds,
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

#[test]
fn fragment_labels_dummies_by_atom_index() {
    let mol = smiles("CCOC(=O)c1ccccc1");
    let mut res = fragment_on_bonds(&mol, &[1, 3], true, None, None, None).unwrap();
    assert_eq!(can(&mut res), "[1*]OC(=[4*])c1ccccc1.[2*]CC.[3*]=O");
}

#[test]
fn fragment_with_labels_types_and_cut_counts() {
    let mol = smiles("CCOCC");
    let mut cuts = vec![0u32; 5];
    let labels = [(7, 8), (9, 10)];
    let mut res = fragment_on_bonds(
        &mol,
        &[1, 2],
        true,
        Some(&labels),
        Some(&[BondType::Single, BondType::Single]),
        Some(&mut cuts),
    )
    .unwrap();
    assert_eq!(cuts, [0, 1, 2, 1, 0]);
    assert_eq!(can(&mut res), "[7*]O[10*].[8*]CC.[9*]CC");
}

#[test]
fn fragment_without_dummies_keeps_aromatic_hs() {
    let mol = smiles("c1cc[nH]c1CC");
    let mut res = fragment_on_bonds(&mol, &[5], false, None, None, None).unwrap();
    assert_eq!(can(&mut res), "C.Cc1ccc[nH]1");
}

#[test]
fn fragment_keeps_tetrahedral_stereo() {
    let mol = smiles("C[C@H](F)Cl");
    let mut res = fragment_on_bonds(&mol, &[0], true, None, None, None).unwrap();
    assert_eq!(can(&mut res), "*[C@H](F)Cl.[1*]C");
}

#[test]
fn fragment_rejects_bad_arguments() {
    let mol = smiles("CCOC");
    assert!(fragment_on_bonds(&mol, &[0, 0], true, None, None, None).is_err());
    assert!(fragment_on_bonds(&mol, &[0, 1], true, Some(&[(1, 2)]), None, None).is_err());
    assert!(fragment_on_bonds(&mol, &[0, 1], true, None, Some(&[BondType::Single]), None).is_err());
    assert!(fragment_on_bonds(&mol, &[9], true, None, None, None).is_err());
}

#[test]
fn fragment_on_some_bonds_enumerates_combinations() {
    let mol = smiles("CCCCO");
    let (mols, cuts) =
        fragment_on_some_bonds(&mol, &[0, 1, 2, 3], 2, true, None, None, true).unwrap();
    let smi: Vec<String> = mols.into_iter().map(|mut m| can(&mut m)).collect();
    assert_eq!(
        smi,
        [
            "*C[2*].[1*]C.[1*]CCO",
            "*CC[3*].[1*]C.[2*]CO",
            "[1*]C[3*].[2*]CC.[2*]CO",
            "*CCC[4*].[1*]C.[3*]O",
            "[1*]CC[4*].[2*]CC.[3*]O",
            "[2*]C[4*].[3*]CCC.[3*]O",
        ]
    );
    assert_eq!(cuts.len(), 6);
    assert_eq!(cuts[0], [1, 2, 1, 0, 0].map(|c| c as u32));
    let (none, _) = fragment_on_some_bonds(&mol, &[0, 1], 3, true, None, None, false).unwrap();
    assert!(none.is_empty());
}

#[test]
fn brics_definitions_parse() {
    let mut environs = std::collections::BTreeMap::new();
    let atom_types = construct_brics_atom_types(Some(&mut environs)).unwrap();
    assert_eq!(atom_types.len(), 15);
    assert_eq!(environs.len(), 15);
    assert_eq!(construct_brics_bond_types().unwrap().len(), 46);
}

#[test]
fn custom_definition_tables_skip_bad_lines() {
    let atoms =
        construct_fragmenter_atom_types("// comment\n1 [C]\n2 [N]\n2 [O]\n3\n", "//", true, None)
            .unwrap();
    assert_eq!(atoms.len(), 2);
    assert_eq!(atoms[&2], "[N]");
    let bonds =
        construct_fragmenter_bond_types("1 2 -;!@\n1 9 -\n", &atoms, "//", true, true).unwrap();
    assert_eq!(bonds.len(), 1);
    assert_eq!((bonds[0].atom1_label, bonds[0].atom2_label), (1, 2));
    assert_eq!(bonds[0].bond_type, BondType::Single);
}

#[test]
fn brics_fragmentation() {
    for (input, expected) in [
        (
            "CCOC(=O)c1ccccc1",
            "[1*]C([6*])=O.[16*]c1ccccc1.[3*]O[3*].[4*]CC",
        ),
        ("C1CCCCC1NC(=O)C", "[1*]C(C)=O.[15*]C1CCCCC1.[5*]N[5*]"),
        ("c1ccccc1", "c1ccccc1"),
    ] {
        let mut res = fragment_on_brics_bonds(&smiles(input)).unwrap();
        assert_eq!(can(&mut res), expected, "{input}");
    }
}
