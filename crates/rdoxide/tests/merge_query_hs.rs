//! `mergeQueryHs` / `hasQueryHs`; expectations from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::mol_ops::query_hs::{has_query_hs, merge_query_hs};
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;

fn smarts(s: &str) -> RWMol {
    mol_from_smarts(s, &SmartsParserParams::default())
        .unwrap()
        .unwrap()
}

fn merged(s: &str, unmapped_only: bool, isotopes: bool) -> RWMol {
    let mut mol = smarts(s);
    merge_query_hs(&mut mol, unmapped_only, isotopes).unwrap();
    mol
}

#[test]
fn hydrogens_become_negated_counts() {
    let mut mol = merged("C([H])[H]", false, false);
    assert_eq!(mol.num_atoms(), 1);
    assert_eq!(
        mol_to_smarts(&mut mol, &SmilesWriteParams::default()).unwrap(),
        "[C&!H0&!H1]"
    );
}

#[test]
fn lone_hydrogens_stay() {
    assert_eq!(merged("[H][H]", false, false).num_atoms(), 2);
    assert_eq!(merged("[H]", false, false).num_atoms(), 1);
}

#[test]
fn maps_and_isotopes_can_block_merging() {
    assert_eq!(merged("[H:1]C", false, false).num_atoms(), 1);
    assert_eq!(merged("[H:1]C", true, false).num_atoms(), 2);
    assert_eq!(merged("[2H]C", false, false).num_atoms(), 2);
    assert_eq!(merged("[2H]C", false, true).num_atoms(), 1);
}

#[test]
fn detects_query_hs() {
    assert_eq!(has_query_hs(&smarts("C[H]")).unwrap(), (true, false));
    assert_eq!(has_query_hs(&smarts("CC")).unwrap(), (false, false));
    assert_eq!(has_query_hs(&smarts("[#1,#6]C")).unwrap(), (true, true));
    assert_eq!(has_query_hs(&smarts("[$(C[H])]C")).unwrap(), (true, false));
}

#[test]
fn smarts_parser_merges_when_asked() {
    let params = SmartsParserParams {
        merge_hs: true,
        ..Default::default()
    };
    let mol = mol_from_smarts("C[H]", &params).unwrap().unwrap();
    assert_eq!(mol.num_atoms(), 1);
}
