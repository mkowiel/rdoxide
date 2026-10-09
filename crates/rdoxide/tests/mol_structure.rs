//! `renumber_atoms`, `get_formal_charge`, `remove_stereochemistry` and
//! `combine_mols`; expectations from RDKit 2026.09.1.

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::chem_transforms::combine_mols;
use rdoxide::graphmol::mol_ops::{get_formal_charge, remove_stereochemistry, renumber_atoms};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::stereo_group::{StereoGroup, StereoGroupType};
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
fn renumber_keeps_stereo_and_reorders_atoms() {
    let mut m = smiles("O[C@H](C)C[C@@H](N)C");
    let mut r = renumber_atoms(&m, &[6, 5, 4, 3, 2, 1, 0]).unwrap();
    assert_eq!(can(&mut r), can(&mut m));
    let symbols: Vec<_> = (0..r.num_atoms())
        .map(|i| r.atom(i).symbol().to_string())
        .collect();
    assert_eq!(symbols, ["C", "N", "C", "C", "C", "C", "O"]);
}

#[test]
fn renumber_rejects_bad_orders() {
    let m = smiles("CCO");
    assert!(renumber_atoms(&m, &[0, 1]).is_err());
    assert!(renumber_atoms(&m, &[0, 0, 1]).is_err());
    assert!(renumber_atoms(&m, &[0, 1, 3]).is_err());
}

#[test]
fn renumber_remaps_stereo_groups() {
    let mut m = smiles("C[C@](F)(Cl)[C@@](Br)(I)C");
    let group = StereoGroup::new(StereoGroupType::Or, vec![1, 4], vec![], 1).unwrap();
    m.set_stereo_groups(vec![group]).unwrap();
    // atom i of the result is atom order[i] of the original
    let r = renumber_atoms(&m, &[7, 6, 5, 4, 3, 2, 1, 0]).unwrap();
    assert_eq!(r.stereo_groups().len(), 1);
    assert_eq!(r.stereo_groups()[0].atoms(), [6, 3]);
}

#[test]
fn renumber_preserves_ring_info() {
    let m = smiles("[C@@H]1(F)CC[C@H](Cl)CC1");
    let r = renumber_atoms(&m, &[2, 0, 1, 3, 4, 5, 6, 7]).unwrap();
    assert_eq!(r.ring_info().num_rings().unwrap(), 1);
}

#[test]
fn formal_charge_sums_atoms() {
    assert_eq!(get_formal_charge(&smiles("CC(=O)[O-].[NH3+]C")), 0);
    assert_eq!(get_formal_charge(&smiles("C[NH3+]")), 1);
}

#[test]
fn remove_stereochemistry_clears_tags() {
    let mut m = smiles("F/C=C/Cl");
    remove_stereochemistry(&mut m).unwrap();
    assert_eq!(can(&mut m), "FC=CCl");
    let mut m = smiles("C[C@H](F)Cl");
    remove_stereochemistry(&mut m).unwrap();
    assert_eq!(can(&mut m), "CC(F)Cl");
}

#[test]
fn combine_without_conformers() {
    let mut c = combine_mols(&smiles("CC"), &smiles("O"), Point3D::new(1.0, 2.0, 3.0)).unwrap();
    assert_eq!(c.num_atoms(), 3);
    assert_eq!(can(&mut c), "CC.O");
}
