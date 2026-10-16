//! Scaffold networks; cases from upstream `ScaffoldNetwork/catch_tests.cpp`.
#![allow(clippy::field_reassign_with_default)]

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::scaffold_network::{
    EdgeType, ScaffoldNetwork, ScaffoldNetworkParams, add_mol_to_network, create_scaffold_network,
    get_brics_network_params, make_scaffold_generic,
};
use rdoxide::graphmol::smiles_parse::write::mol_to_smiles;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn parse(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn count(net: &ScaffoldNetwork, t: EdgeType) -> usize {
    net.edges.iter().filter(|e| e.edge_type == t).count()
}

fn sorted_nodes(net: &ScaffoldNetwork) -> Vec<String> {
    let mut nodes = net.nodes.clone();
    nodes.sort();
    nodes
}

fn penicillins() -> Vec<RWMol> {
    [
        "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12",
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
        "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na]",
        "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12",
    ]
    .map(parse)
    .to_vec()
}

#[test]
fn larger_multi_mol_network() {
    let mut ps = ScaffoldNetworkParams::default();
    ps.include_generic_scaffolds = false;
    ps.include_scaffolds_without_attachments = false;
    let net = create_scaffold_network(&penicillins(), &ps).unwrap();
    assert_eq!(net.nodes.len(), 11);
    assert_eq!(net.counts.len(), net.nodes.len());
    assert_eq!(net.edges.len(), 14);
    assert_eq!(count(&net, EdgeType::Fragment), 10);
    assert_eq!(count(&net, EdgeType::Initialize), 4);
    let expected = [
        "*C1C(=O)N2CCSC12",
        "*c1ccccc1",
        "*c1conc1*",
        "*c1conc1-c1ccccc1",
        "*c1nocc1C(=O)NC1C(=O)N2CCSC12",
        "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na]",
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
        "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H]1SC(C)(C)[C@@H]2C(=O)O",
        "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H]1SC(C)(C)[C@@H]2C(=O)O",
        "O=C(Cc1ccccc1)NC1C(=O)N2CCSC12",
        "O=C(NC1C(=O)N2CCSC12)c1conc1-c1ccccc1",
    ];
    assert_eq!(sorted_nodes(&net), expected);

    ps.include_generic_scaffolds = true;
    let net = create_scaffold_network(&penicillins(), &ps).unwrap();
    assert_eq!(net.nodes.len(), 18);
    assert_eq!(net.edges.len(), 21);
    assert_eq!(count(&net, EdgeType::Generic), 7);
}

#[test]
fn brics_fragmenter() {
    let m = parse("c1ccccc1C(=O)NC1NC(=O)CCC1");
    let mut ps = ScaffoldNetworkParams::default();
    ps.include_scaffolds_without_attachments = false;
    ps.include_generic_scaffolds = false;
    ps.keep_only_first_fragment = false;
    let mut net = ScaffoldNetwork::default();
    add_mol_to_network(&m, &mut net, &ps).unwrap();
    assert_eq!((net.nodes.len(), net.edges.len()), (6, 8));
    assert_eq!(count(&net, EdgeType::Fragment), 8);
    assert_eq!(net.counts.iter().filter(|&&c| c == 1).count(), 3);
    assert_eq!(net.counts.iter().filter(|&&c| c == 2).count(), 3);

    let mut ps = get_brics_network_params().unwrap();
    ps.include_scaffolds_without_attachments = false;
    ps.include_generic_scaffolds = false;
    let mut net = ScaffoldNetwork::default();
    add_mol_to_network(&m, &mut net, &ps).unwrap();
    assert_eq!((net.nodes.len(), net.edges.len()), (10, 20));
    assert_eq!(count(&net, EdgeType::Fragment), 20);
}

#[test]
fn mol_counts() {
    let ms = [
        "C1CC(C1)C1C(C1C1NCCCC1)C1OCCC1",
        "C1CC(C1)C1C(C1C1NCCCC1)C1CCCC1",
    ]
    .map(parse);
    let mut ps = ScaffoldNetworkParams::default();
    ps.include_generic_scaffolds = false;
    ps.include_scaffolds_without_attachments = false;
    let net = create_scaffold_network(&ms, &ps).unwrap();
    assert_eq!((net.nodes.len(), net.edges.len()), (16, 40));
    for (smi, expected) in [
        ("*C1CCC1", 2),
        ("*C1CCCCN1", 2),
        ("*C1CCCO1", 1),
        ("*C1CCCC1", 1),
        ("*C1C(*)C1*", 2),
    ] {
        let idx = net.nodes.iter().position(|n| n == smi).expect(smi);
        assert_eq!(net.mol_counts[idx], expected, "{smi}");
    }
}

#[test]
fn names_and_generic_isotopes() {
    let m = parse("CC(=O)Oc1ccccc1C(=O)O aspirin");
    let mut net = ScaffoldNetwork::default();
    add_mol_to_network(&m, &mut net, &ScaffoldNetworkParams::default()).unwrap();
    assert_eq!(net.nodes[0], "CC(=O)Oc1ccccc1C(=O)O");
    let mut ps = ScaffoldNetworkParams::default();
    ps.include_names = true;
    let mut net = ScaffoldNetwork::default();
    add_mol_to_network(&m, &mut net, &ps).unwrap();
    assert_eq!(net.nodes[0], "CC(=O)Oc1ccccc1C(=O)O aspirin");

    let m = parse("[235U]1CC1");
    let g = make_scaffold_generic(&m, true, false);
    assert_eq!(
        mol_to_smiles(&mut g.clone(), &Default::default()).unwrap(),
        "*1**1"
    );
}

#[test]
fn bad_parameters_are_rejected() {
    assert!(ScaffoldNetworkParams::new(&["[C:1]>>[C:1]"]).is_err());
    assert!(ScaffoldNetworkParams::new(&["[C:1][C:2]>>[C:1].[C:2]"]).is_ok());
    let mut ps = ScaffoldNetworkParams::default();
    ps.include_scaffolds_with_attachments = false;
    ps.include_scaffolds_without_attachments = false;
    assert!(create_scaffold_network(&[parse("CC")], &ps).is_err());
}
