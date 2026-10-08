//! Query model and SMARTS parser tests (cf. upstream testQuery.cpp,
//! catch_queries.cpp, smaTest.cpp).

use rdoxide::graphmol::queries::*;
use rdoxide::graphmol::smiles_parse::smarts::{
    SmartsParserParams, atom_from_smarts, bond_from_smarts, mol_from_smarts,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{BondType, RWMol};

fn mol(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn atom_matches(smarts: &str, m: &RWMol) -> Vec<u32> {
    let q = atom_from_smarts(smarts).unwrap();
    let q = q.query().unwrap();
    (0..m.num_atoms())
        .filter(|&i| q.matches(m, i).unwrap())
        .collect()
}

#[test]
fn atom_query_matching() {
    let m = mol("OCC(=O)c1ccccc1N");
    assert_eq!(atom_matches("[#6]", &m), vec![1, 2, 4, 5, 6, 7, 8, 9]);
    assert_eq!(atom_matches("[c]", &m), vec![4, 5, 6, 7, 8, 9]);
    assert_eq!(atom_matches("[C]", &m), vec![1, 2]);
    assert_eq!(atom_matches("[O,N]", &m), vec![0, 3, 10]);
    assert_eq!(atom_matches("[!#6;!#1]", &m), vec![0, 3, 10]);
    assert_eq!(atom_matches("[CH2]", &m), vec![1]);
    assert_eq!(atom_matches("[D3]", &m), vec![2, 4, 9]);
    assert_eq!(atom_matches("[D{1-1}]", &m), vec![0, 3, 10]);
    assert_eq!(atom_matches("[R]", &m), vec![4, 5, 6, 7, 8, 9]);
    assert_eq!(atom_matches("[r6]", &m), vec![4, 5, 6, 7, 8, 9]);
    assert_eq!(atom_matches("[x2]", &m), vec![4, 5, 6, 7, 8, 9]);
    assert_eq!(atom_matches("[z1]", &m), vec![1, 2, 9]);
    assert_eq!(atom_matches("[*]", &m).len(), 11);
    assert_eq!(atom_matches("[a]", &m), vec![4, 5, 6, 7, 8, 9]);
    assert_eq!(atom_matches("[A]", &m), vec![0, 1, 2, 3, 10]);
    assert_eq!(atom_matches("[^2]", &m), vec![2, 3, 4, 5, 6, 7, 8, 9, 10]);
}

#[test]
fn bond_query_matching() {
    let m = mol("C=Cc1ccccc1");
    let matches = |s: &str| {
        let b = bond_from_smarts(s).unwrap();
        let q = b.query().unwrap().clone();
        (0..m.num_bonds())
            .filter(|&i| q.matches(&m, i).unwrap())
            .collect::<Vec<u32>>()
    };
    assert_eq!(matches("="), vec![0]);
    assert_eq!(matches(":"), vec![2, 3, 4, 5, 6, 7]);
    assert_eq!(matches("@"), vec![2, 3, 4, 5, 6, 7]);
    assert_eq!(matches("!@"), vec![0, 1]);
    assert_eq!(matches("~"), (0..8).collect::<Vec<u32>>());
    assert_eq!(matches("-,="), vec![0, 1]);
}

#[test]
fn smarts_structure() {
    let m = mol_from_smarts("[C;H2]-,=[N,O]", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(m.num_atoms(), 2);
    assert_eq!(
        m.atom(0).query().unwrap().describe(0),
        "AtomAnd\n  AtomType 6 = val\n  AtomHCount 2 = val\n"
    );
    assert_eq!(m.bond(0).query().unwrap().description, "BondOr");
    let m = mol_from_smarts("cC", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(m.bond(0).bond_type(), BondType::Single);
    assert!(
        mol_from_smarts("[C", &SmartsParserParams::default())
            .unwrap()
            .is_none()
    );
}

#[test]
fn expand_query_null_algebra() {
    // AND with a null query keeps the other query
    let q = expand_query(
        make_atom_null_query(),
        make_atom_num_query(6),
        CompositeQueryType::And,
        true,
        "Atom",
    )
    .unwrap();
    assert_eq!(q.description, "AtomAtomicNum");
    // OR with a null query keeps the null query
    let q = expand_query(
        make_atom_null_query(),
        make_atom_num_query(6),
        CompositeQueryType::Or,
        true,
        "Atom",
    )
    .unwrap();
    assert_eq!(q.description, "AtomNull");
    let q = expand_query(
        make_atom_num_query(7),
        make_atom_num_query(6),
        CompositeQueryType::XOr,
        false,
        "Atom",
    )
    .unwrap();
    assert_eq!(
        q.describe(0),
        "AtomXor\n  AtomAtomicNum 6 = val\n  AtomAtomicNum 7 = val\n"
    );
}

#[test]
fn recursive_labels() {
    use rdoxide::graphmol::smiles_parse::smarts::label_recursive_patterns;
    assert_eq!(
        label_recursive_patterns("[$(C=O),$(C#N)]"),
        "[$(C=O)_100,$(C#N)_101]"
    );
    assert_eq!(
        label_recursive_patterns("[$(C=O)]C[$(C=O)]"),
        "[$(C=O)_100]C[$(C=O)_100]"
    );
    assert_eq!(label_recursive_patterns("[$(C=O)_5]"), "[$(C=O)_5]");
}
