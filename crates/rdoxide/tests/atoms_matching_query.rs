//! `atoms_matching_query` (`Mol.GetAtomsMatchingQuery`); expectations from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::smiles_parse::smarts::atom_from_smarts;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::atoms_matching_query;

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn matching(mol: &RWMol, smarts: &str) -> Vec<u32> {
    let qa = atom_from_smarts(smarts).expect("valid SMARTS");
    atoms_matching_query(mol, &qa).unwrap()
}

#[test]
fn simple_queries() {
    let mol = smiles("CC(=O)Nc1ccccc1");
    assert_eq!(matching(&mol, "[#7]"), [3]);
    assert_eq!(matching(&mol, "a"), [4, 5, 6, 7, 8, 9]);
    assert_eq!(matching(&mol, "[!a]"), [0, 1, 2, 3]);
    assert_eq!(matching(&mol, "[C;H3]"), [0]);
}

#[test]
fn recursive_queries_see_the_whole_molecule() {
    let mol = smiles("CCO");
    assert_eq!(matching(&mol, "[$(C-O)]"), [1]);
    assert_eq!(matching(&mol, "[$(C-C)]"), [0, 1]);
}
