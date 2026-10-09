//! `copyMolSubset` and `pathToSubmol`; expectations from RDKit 2026.09.1.

use rdoxide::RDKitError;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::bond::BondStereo;
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::subgraphs::path_to_submol;
use rdoxide::graphmol::subset::{
    SubsetMethod, SubsetOptions, copy_mol_subset, copy_mol_subset_atoms_bonds,
};

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn can(mol: &mut RWMol) -> String {
    mol_to_smiles(mol, &SmilesWriteParams::default()).unwrap()
}

#[test]
fn atom_paths_keep_the_bonds_between_them() {
    let mol = smiles("CCOC");
    let (mut sub, info) = copy_mol_subset(&mol, &[0, 1, 2], &SubsetOptions::default()).unwrap();
    assert_eq!(can(&mut sub), "CCO");
    assert_eq!(info.atom_mapping.len(), 3);
    assert_eq!(
        info.bond_mapping.keys().copied().collect::<Vec<_>>(),
        [0, 1]
    );
    // indices outside the molecule are ignored
    let (sub, _) = copy_mol_subset(&mol, &[1, 2, 99], &SubsetOptions::default()).unwrap();
    assert_eq!(sub.num_atoms(), 2);
}

#[test]
fn bond_paths_bring_their_atoms() {
    let mol = smiles("CCOC");
    let options = SubsetOptions {
        method: SubsetMethod::Bonds,
        ..SubsetOptions::default()
    };
    let (mut sub, info) = copy_mol_subset(&mol, &[1, 2], &options).unwrap();
    assert_eq!(can(&mut sub), "COC");
    let atoms: Vec<_> = info.atom_mapping.into_iter().collect();
    assert_eq!(atoms, [(1, 0), (2, 1), (3, 2)]);
}

#[test]
fn explicit_atoms_and_bonds_are_validated() {
    let mol = smiles("CCOC");
    let options = SubsetOptions::default();
    let (mut sub, _) = copy_mol_subset_atoms_bonds(&mol, &[0, 1, 2], &[0, 1], &options).unwrap();
    assert_eq!(can(&mut sub), "CCO");
    let err = copy_mol_subset_atoms_bonds(&mol, &[0, 2], &[0], &options).unwrap_err();
    assert!(matches!(&err, RDKitError::Value(m) if m.contains("not contained in subset atoms")));
    let err = copy_mol_subset_atoms_bonds(&mol, &[0, 10], &[0], &options).unwrap_err();
    assert!(matches!(&err, RDKitError::Index(m) if m == "10"), "{err}");
    // everything selected: a plain copy without a mapping
    let (mut all, info) =
        copy_mol_subset_atoms_bonds(&mol, &[0, 1, 2, 3], &[0, 1, 2], &options).unwrap();
    assert_eq!(can(&mut all), "CCOC");
    assert!(info.atom_mapping.is_empty());
}

#[test]
fn double_bond_stereo_is_translated_or_dropped() {
    let mol = smiles("F/C=C/C(Cl)[C@H](O)N");
    let (mut sub, _) = copy_mol_subset(&mol, &[0, 1, 2, 3], &SubsetOptions::default()).unwrap();
    assert_eq!(can(&mut sub), "C/C=C/F");
    let double = sub.bond(1);
    assert_eq!(double.stereo(), BondStereo::StereoTrans);
    assert_eq!(double.stereo_atoms(), [0, 3]);
    let (mut dropped, _) = copy_mol_subset(&mol, &[0, 1, 2], &SubsetOptions::default()).unwrap();
    assert_eq!(can(&mut dropped), "C=CF");
}

#[test]
fn copy_as_query_makes_query_atoms_and_bonds() {
    let mol = smiles("CCOC");
    let options = SubsetOptions {
        copy_as_query: true,
        ..SubsetOptions::default()
    };
    let (sub, _) = copy_mol_subset(&mol, &[0, 1, 2], &options).unwrap();
    assert!(sub.atom(0).has_query() && sub.bond(0).has_query());
}

#[test]
fn path_to_submol_maps_atoms() {
    let mol = smiles("OC1CC1C[C@H](F)Cl");
    let (mut sub, mapping) = path_to_submol(&mol, &[0, 1], false).unwrap();
    assert_eq!(can(&mut sub), "CCO");
    assert_eq!(
        mapping.into_iter().collect::<Vec<_>>(),
        [(0, 0), (1, 1), (2, 2)]
    );
    let (empty, mapping) = path_to_submol(&mol, &[], false).unwrap();
    assert_eq!(empty.num_atoms(), 0);
    assert!(mapping.is_empty());
    // out-of-range and negative bond indices are ignored
    let (sub, _) = path_to_submol(&mol, &[0, 99, -1], false).unwrap();
    assert_eq!(sub.num_atoms(), 2);
}
