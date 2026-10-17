//! `add_recursive_query` and `parse_query_def_text`.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::chem_transforms::recursive_queries::{
    add_recursive_query, parse_query_def_text,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::{SubstructMatchParameters, substruct_match};

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn matches(mol: &RWMol, query: &RWMol) -> bool {
    !substruct_match(mol, query, &SubstructMatchParameters::default())
        .unwrap()
        .is_empty()
}

#[test]
fn recursive_queries_are_anded_to_the_atom() {
    let mut q1 = smiles("CC");
    add_recursive_query(&mut q1, &smiles("CO"), 1, true).unwrap();
    assert!(matches(&smiles("OCC"), &q1));
    let q3 = smiles("CS");
    add_recursive_query(&mut q1, &q3, 1, true).unwrap();
    assert!(!matches(&smiles("OCC"), &q1));
    assert!(matches(&smiles("OC(S)C"), &q1));
}

#[test]
fn atom_index_is_checked() {
    let mut q = smiles("CC");
    assert!(add_recursive_query(&mut q, &smiles("CO"), 5, true).is_err());
}

#[test]
fn query_definitions_are_parsed_from_text() {
    let text = "// a comment\nCarboxylicAcid\tC(=O)[OH]\n\nbad\t[\nAmine\t[NX3;H2]\n";
    let defs = parse_query_def_text(text, true, "\t", "//", 0, 1);
    assert_eq!(
        defs.keys().cloned().collect::<Vec<_>>(),
        ["amine", "carboxylicacid"]
    );
    let defs = parse_query_def_text(text, false, "\t", "//", 0, 1);
    assert!(defs.contains_key("CarboxylicAcid"));
    assert!(matches(&smiles("CC(=O)O"), &defs["CarboxylicAcid"]));
}
