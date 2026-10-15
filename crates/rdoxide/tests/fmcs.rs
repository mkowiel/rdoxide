//! Maximum common substructure search (`graphmol::fmcs`).

use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::fmcs::*;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn mols(smiles: &[&str]) -> Vec<ROMol> {
    smiles
        .iter()
        .map(|s| {
            mol_from_smiles(s, &SmilesParserParams::default())
                .unwrap()
                .unwrap()
        })
        .collect()
}

fn mcs(smiles: &[&str], params: &MCSParameters) -> MCSResult {
    let ms = mols(smiles);
    let refs: Vec<&ROMol> = ms.iter().collect();
    find_mcs(&refs, Some(params)).unwrap()
}

#[test]
fn benzene_in_substituted_benzenes() {
    let r = mcs(
        &["c1ccccc1O", "c1ccccc1N", "c1ccccc1C"],
        &MCSParameters::default(),
    );
    assert_eq!((r.num_atoms, r.num_bonds), (6, 6));
    assert_eq!(r.smarts_string, "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1");
    assert!(r.is_completed());
}

#[test]
fn chains_of_different_length() {
    let r = mcs(&["CCCCC", "CCCC", "CCC"], &MCSParameters::default());
    assert_eq!((r.num_atoms, r.num_bonds), (3, 2));
    assert_eq!(r.smarts_string, "[#6]-[#6]-[#6]");
}

#[test]
fn no_common_substructure() {
    let r = mcs(&["CC", "OO"], &MCSParameters::default());
    assert_eq!((r.num_atoms, r.num_bonds), (0, 0));
    assert_eq!(r.smarts_string, "");
}

#[test]
fn single_common_atom() {
    let r = mcs(&["CO", "NO"], &MCSParameters::default());
    assert_eq!((r.num_atoms, r.num_bonds), (1, 0));
    assert_eq!(r.smarts_string, "[#8]");
}

#[test]
fn complete_rings_only() {
    let mut p = MCSParameters::default();
    p.bond_compare_parameters.complete_rings_only = true;
    let r = mcs(&["C1CCCCC1CCC", "C1CCCCC1C"], &p);
    assert_eq!((r.num_atoms, r.num_bonds), (7, 7));
    p.bond_compare_parameters.complete_rings_only = true;
    let r = mcs(&["C1CCC1CCCC", "C1CCCCC1CCCC"], &p);
    // the cyclobutane is no ring of the second molecule
    assert!(r.num_atoms < 8);
}

#[test]
fn threshold_allows_a_missing_molecule() {
    let p = MCSParameters {
        threshold: 0.6,
        ..Default::default()
    };
    let r = mcs(&["c1ccccc1O", "c1ccccc1N", "CCCC"], &p);
    assert_eq!((r.num_atoms, r.num_bonds), (6, 6));
}

#[test]
fn store_all_degenerate() {
    let p = MCSParameters {
        store_all: true,
        ..Default::default()
    };
    let r = mcs(&["CC(C)N", "CC(C)O"], &p);
    assert_eq!((r.num_atoms, r.num_bonds), (3, 2));
    assert!(!r.degenerate_smarts_query_mol_dict.is_empty());
}

#[test]
fn atom_and_bond_comparators() {
    let mut p = MCSParameters::default();
    p.set_atom_typer_from_enum(AtomComparator::Any);
    p.set_bond_typer_from_enum(BondComparator::Any);
    let r = mcs(&["CCO", "CCN"], &p);
    assert_eq!((r.num_atoms, r.num_bonds), (3, 2));
    let r = mcs(&["C=CC", "CCC"], &MCSParameters::default());
    assert_eq!((r.num_atoms, r.num_bonds), (2, 1));
}

#[test]
fn json_parameters() {
    let mut p = MCSParameters::default();
    parse_mcs_parameters_json(
        r#"{"Threshold": 0.5, "RingMatchesRingOnly": true, "AtomCompare": "Any"}"#,
        &mut p,
    )
    .unwrap();
    assert_eq!(p.threshold, 0.5);
    assert!(p.atom_compare_parameters.ring_matches_ring_only);
    assert!(matches!(p.atom_typer, AtomTyper::Any));
}

#[test]
fn rejects_bad_arguments() {
    let ms = mols(&["CC"]);
    let refs: Vec<&ROMol> = ms.iter().collect();
    assert!(find_mcs(&refs, None).is_err());
}
