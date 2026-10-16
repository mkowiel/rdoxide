//! Molecule enumeration; cases from upstream `MolEnumerator/Wrap/rough_test.py`.

use std::sync::Arc;

use rdoxide::graphmol::file_parsers::mol_file_parser::{MolFileParserParams, mol_from_mol_block};
use rdoxide::graphmol::mol_enumerator::{
    LinkNodeOp, MolEnumeratorOp, MolEnumeratorParams, PositionVariationOp, RepeatUnitOp, enumerate,
    enumerate_default,
};
use rdoxide::graphmol::smiles_parse::write::mol_to_smiles;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

const LINK_NODE: &str = "one linknode
  Mrv2007 06222005102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 6.9164 14.4947 0 0
M  V30 4 C 9.5836 14.4947 0 0
M  V30 5 C 9.5836 12.9547 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 END BOND
M  V30 LINKNODE 1 4 2 1 2 1 5
M  V30 END CTAB
M  END";

const POSITION_VARIATION: &str = "two position variation bonds
  Mrv2007 06242006032D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.7083 2.415 0 0
M  V30 2 C -3.042 1.645 0 0
M  V30 3 C -3.042 0.105 0 0
M  V30 4 N -1.7083 -0.665 0 0
M  V30 5 C -0.3747 0.105 0 0
M  V30 6 C -0.3747 1.645 0 0
M  V30 7 * -3.042 0.875 0 0
M  V30 8 F -5.0434 0.875 0 0
M  V30 9 * -1.0415 2.03 0 0
M  V30 10 Cl -1.0415 4.34 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(2 2 3) ATTACH=ANY
M  V30 8 1 9 10 ENDPTS=(2 1 6) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
";

fn parse(block: &str) -> rdoxide::graphmol::RWMol {
    mol_from_mol_block(block, &MolFileParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn link_nodes() {
    let m = parse(LINK_NODE);
    let ps = MolEnumeratorParams::with_operation(Arc::new(LinkNodeOp::new()));
    let bundle = enumerate(&m, &[ps]).unwrap();
    assert_eq!(bundle.size(), 4);
    let sizes: Vec<u32> = bundle.mols().iter().map(|m| m.num_atoms()).collect();
    assert_eq!(sizes, vec![6, 8, 10, 12]);
}

#[test]
fn position_variation() {
    let m = parse(POSITION_VARIATION);
    let ps = MolEnumeratorParams::with_operation(Arc::new(PositionVariationOp::new()));
    let bundle = enumerate(&m, &[ps]).unwrap();
    assert_eq!(bundle.size(), 4);
    assert!(
        bundle
            .mols()
            .iter()
            .all(|m| m.num_atoms() == 8 && m.num_bonds() == 8)
    );
}

#[test]
fn nothing_to_enumerate_gives_an_empty_bundle() {
    let m = parse(LINK_NODE);
    let ps = MolEnumeratorParams::with_operation(Arc::new(PositionVariationOp::new()));
    assert!(enumerate(&m, &[ps]).unwrap().is_empty());
    assert!(enumerate(&m, &[MolEnumeratorParams::default()]).is_err());
}

#[test]
fn max_to_enumerate_limits_the_result() {
    let m = parse(LINK_NODE);
    let bundle = enumerate_default(&m, 2).unwrap();
    assert_eq!(bundle.size(), 2);
}

fn smiles(m: &rdoxide::graphmol::ROMol) -> String {
    mol_to_smiles(&mut m.clone(), &Default::default()).unwrap()
}

fn from_cxsmiles(s: &str) -> rdoxide::graphmol::RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn repeat_units_ht_and_hh() {
    let m = from_cxsmiles("*-CO-* |$star_e;;;star_e$,Sg:n:2,1::ht|");
    let op = RepeatUnitOp::from_mol(&m).unwrap();
    assert_eq!(op.variation_counts(), vec![4]);
    for (which, expected) in [(0, "**"), (1, "*CO*"), (3, "*COCOCO*")] {
        assert_eq!(smiles(&op.apply(&[which]).unwrap()), expected);
    }
    let m = from_cxsmiles("*-CO-* |$star_e;;;star_e$,Sg:n:2,1::hh|");
    let op = RepeatUnitOp::from_mol(&m).unwrap();
    assert_eq!(smiles(&op.apply(&[3]).unwrap()), "*COOCCO*");
}

#[test]
fn several_repeat_units() {
    let m = from_cxsmiles("*-CN[13CH2]OC-* |$star_e;;;;;;star_e$,Sg:n:4,5::ht,Sg:n:1,2::ht|");
    let op = RepeatUnitOp::from_mol(&m).unwrap();
    assert_eq!(op.variation_counts(), vec![4, 4]);
    for (which, expected) in [
        ([0, 0], "*[13CH2]*"),
        ([1, 0], "*CO[13CH2]*"),
        ([0, 1], "*CN[13CH2]*"),
        ([1, 1], "*CN[13CH2]OC*"),
        ([3, 1], "*CN[13CH2]OCOCOC*"),
        ([1, 3], "*CNCNCN[13CH2]OC*"),
        ([2, 2], "*CNCN[13CH2]OCOC*"),
    ] {
        assert_eq!(smiles(&op.apply(&which).unwrap()), expected);
    }
}

#[test]
fn repeat_counts_follow_the_labels() {
    for (label, expected) in [
        ("", 4),
        ("n", 4),
        ("hello", 4),
        ("23-20", 4),
        ("-20", 4),
        ("2-4", 3),
        ("20", 21),
        ("20-20", 1),
    ] {
        let m = from_cxsmiles(&format!("FCN(CC)-* |Sg:n:2,5:{label}:ht|"));
        assert_eq!(
            enumerate_default(&m, 0).unwrap().size(),
            expected,
            "{label}"
        );
        assert_eq!(enumerate_default(&m, 1).unwrap().size(), 1, "{label}");
    }
    let m = from_cxsmiles("FCN(CC)-* |Sg:n:2,5:2-3:ht|");
    let bundle = enumerate_default(&m, 0).unwrap();
    assert_eq!(smiles(&bundle.mols()[0]), "*N(CC)N(*)CF");
    assert_eq!(smiles(&bundle.mols()[1]), "*N(CC)N(*)N(*)CF");
}
