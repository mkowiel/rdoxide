//! RASCAL MCES; the cases of upstream `RascalMCES/Wrap/testRascalMCES.py`.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::rascal_mces::{RascalOptions, rascal_mces};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn mol(smiles: &str) -> RWMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn simple_aromatic_match() {
    let res = rascal_mces(
        &mol("c1ccccc1Cl"),
        &mol("c1ccccc1F"),
        &RascalOptions::default(),
    )
    .unwrap();
    assert_eq!(res.len(), 1);
    assert_eq!(res[0].smarts().unwrap(), "c1:c:c:c:c:c:1");
    assert_eq!(res[0].bond_matches().len(), 6);
    assert_eq!(res[0].atom_matches().len(), 6);
}

#[test]
fn largest_fragment() {
    let ad1 = mol("CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1");
    let ad2 = mol("N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1");
    let mut res = rascal_mces(&ad1, &ad2, &RascalOptions::default()).unwrap();
    assert_eq!(res.len(), 1);
    assert_eq!(
        res[0].smarts().unwrap(),
        "N(-C)-c1:c:c:c(-CC(=O)-NCCCCCCCCCC):c:c:1.NC12CC3CC(-C1)-CC(-C2)-C3"
    );
    res[0].largest_frag_only().unwrap();
    assert_eq!(
        res[0].smarts().unwrap(),
        "N(-C)-c1:c:c:c(-CC(=O)-NCCCCCCCCCC):c:c:1"
    );
    let opts = RascalOptions {
        single_largest_frag: true,
        ..Default::default()
    };
    let res = rascal_mces(&ad1, &ad2, &opts).unwrap();
    assert_eq!(res.len(), 1);
    assert_eq!(
        res[0].smarts().unwrap(),
        "CNc1:c:c:c(-CC(=O)-NCCCCCCCCCC):c:c:1"
    );
}

#[test]
fn thresholds() {
    let m1 = mol("Oc1cccc2C(=O)C=CC(=O)c12");
    let m2 = mol("O1C(=O)C=Cc2cc(OC)c(O)cc12");
    assert!(
        rascal_mces(&m1, &m2, &RascalOptions::default())
            .unwrap()
            .is_empty()
    );
    let opts = RascalOptions {
        similarity_threshold: 0.5,
        ..Default::default()
    };
    assert_eq!(rascal_mces(&m1, &m2, &opts).unwrap().len(), 1);

    let ad1 = mol("CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1");
    let ad2 = mol("N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1");
    let mut opts = RascalOptions {
        similarity_threshold: 0.95,
        ..Default::default()
    };
    assert!(rascal_mces(&ad1, &ad2, &opts).unwrap().is_empty());
    opts.return_empty_mces = true;
    assert_eq!(rascal_mces(&ad1, &ad2, &opts).unwrap().len(), 1);
}

#[test]
fn max_bond_match_pairs() {
    let mut opts = RascalOptions {
        similarity_threshold: 0.0,
        return_empty_mces: true,
        complete_aromatic_rings: false,
        timeout: -1,
        ..Default::default()
    };
    let m1 = mol("CCCC=CCCCC=CCCCCCCCCCCCCCCCCCCCCC1CNCCC1");
    let m2 = mol("CCCC=CCCCCC=CCCCCCCCCCCCCCCCCCCCC1CNCCC1");
    let res = rascal_mces(&m1, &m2, &opts).unwrap();
    assert_eq!(res[0].bond_matches().len(), 0);
    opts.max_bond_match_pairs = 1200;
    let res = rascal_mces(&m1, &m2, &opts).unwrap();
    assert_eq!(res[0].bond_matches().len(), 34);
}

#[test]
fn exact_connections_equivalent_atoms_and_bonds() {
    let mut opts = RascalOptions {
        similarity_threshold: 0.5,
        all_best_mcess: true,
        ..Default::default()
    };
    let m1 = mol("c1ccccc1C1CCC(C(C)C)C1");
    let m2 = mol("c1ccccc1C(C)C");
    let res = rascal_mces(&m1, &m2, &opts).unwrap();
    assert_eq!(res[0].num_frags().unwrap(), 1);
    assert_eq!(res[0].smarts().unwrap(), "c1:c:c:c:c:c:1-C(-C)-C");
    opts.exact_connections_match = true;
    let res = rascal_mces(&m1, &m2, &opts).unwrap();
    assert_eq!(res[0].num_frags().unwrap(), 2);
    assert_eq!(
        res[0].smarts().unwrap(),
        "[#6&a&D2]1:[#6&a&D2]:[#6&a&D2]:[#6&a&D2]:[#6&a&D2]:[#6&a&D3]:1.[#6&A&D3](-[#6&A&D1])-[#6&A&D1]"
    );

    let opts = RascalOptions {
        similarity_threshold: 0.5,
        equivalent_atoms: "[F,Cl,Br,I]".into(),
        ..Default::default()
    };
    let res = rascal_mces(&mol("c1ccccc1F"), &mol("c1ccccc1Br"), &opts).unwrap();
    assert_eq!(res[0].num_frags().unwrap(), 1);
    assert_eq!(res[0].smarts().unwrap(), "c1:c:c:c:c:c:1-[F,Cl,Br,I]");

    let opts = RascalOptions {
        similarity_threshold: 0.5,
        ignore_bond_orders: true,
        ..Default::default()
    };
    let res = rascal_mces(&mol("CC=CC"), &mol("CCCC"), &opts).unwrap();
    assert_eq!(res[0].smarts().unwrap(), "C~C~C~C");

    let opts = RascalOptions {
        similarity_threshold: 0.1,
        ignore_atom_aromaticity: false,
        ..Default::default()
    };
    let res = rascal_mces(&mol("c1ccccc1NCC"), &mol("C1CCCCC1NCC"), &opts).unwrap();
    assert_eq!(res[0].smarts().unwrap(), "NCC");

    let opts = RascalOptions {
        similarity_threshold: 0.1,
        min_clique_size: 17,
        ..Default::default()
    };
    let res = rascal_mces(
        &mol("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
        &mol("CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"),
        &opts,
    )
    .unwrap();
    assert!(res.is_empty());
}
