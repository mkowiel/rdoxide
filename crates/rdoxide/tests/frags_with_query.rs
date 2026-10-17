//! `get_mol_frags_with_query`.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::mol_ops::frags::get_mol_frags_with_query;
use rdoxide::graphmol::monomer_info::AtomMonomerInfo;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn residues() -> RWMol {
    // C-C-O-C, atoms 0,1 in residue ALA of chain A and 2,3 in SER of chain B
    let mut mol = mol_from_smiles("CCOC", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    for (idx, (res, chain)) in [("ALA", "A"), ("ALA", "A"), ("SER", "B"), ("SER", "B")]
        .into_iter()
        .enumerate()
    {
        let mut info = AtomMonomerInfo::new_pdb(" CA ");
        info.residue_name = res.into();
        info.chain_id = chain.into();
        mol.atom_mut(idx as u32).set_monomer_info(Some(info));
    }
    mol
}

fn residue(_: &RWMol, a: &rdoxide::graphmol::atom::Atom) -> String {
    a.monomer_info()
        .map(|i| i.residue_name.clone())
        .unwrap_or_default()
}

#[test]
fn atoms_are_grouped_with_their_bonds() {
    let mol = residues();
    let frags = get_mol_frags_with_query(&mol, residue, false, None, false).unwrap();
    assert_eq!(frags.keys().cloned().collect::<Vec<_>>(), ["ALA", "SER"]);
    assert_eq!(frags["ALA"].num_atoms(), 2);
    assert_eq!(frags["ALA"].num_bonds(), 1);
    assert_eq!(frags["SER"].num_atoms(), 2);
    assert_eq!(frags["SER"].num_bonds(), 1);
}

#[test]
fn white_list_and_negation() {
    let mol = residues();
    let list = ["SER".to_string()];
    let only = get_mol_frags_with_query(&mol, residue, false, Some(&list), false).unwrap();
    assert_eq!(only.keys().cloned().collect::<Vec<_>>(), ["SER"]);
    let rest = get_mol_frags_with_query(&mol, residue, false, Some(&list), true).unwrap();
    assert_eq!(rest.keys().cloned().collect::<Vec<_>>(), ["ALA"]);
}

#[test]
fn molecules_without_pdb_info_go_under_the_empty_key() {
    let mol = mol_from_smiles("CC", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let frags = get_mol_frags_with_query(&mol, residue, false, None, false).unwrap();
    assert_eq!(frags.keys().cloned().collect::<Vec<_>>(), [""]);
    assert_eq!(frags[""].num_atoms(), 2);
}
