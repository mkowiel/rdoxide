//! `SLNParse/test.cpp`, translated.
#![allow(non_snake_case, unused_assignments, unused_mut, unused_variables)]

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::bond::BondType;
use rdoxide::graphmol::mol_ops::{sanitize_flags::SANITIZE_ALL, sanitize_mol};
use rdoxide::graphmol::sln_parse::{sln_query_to_mol, sln_to_mol};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use rdoxide::rdgeneral::{RDProps, RDValue, common_properties as cp};

fn m(o: &Option<RWMol>) -> &RWMol {
    o.as_ref().expect("null molecule")
}

fn mm(o: &mut Option<RWMol>) -> &mut RWMol {
    o.as_mut().expect("null molecule")
}

fn smiles(s: &str, sanitize: bool) -> Option<RWMol> {
    let params = SmilesParserParams {
        sanitize,
        remove_hs: sanitize,
        ..Default::default()
    };
    mol_from_smiles(s, &params).unwrap()
}

fn sm(mol: &Option<RWMol>, patt: &Option<RWMol>, mv: &mut Vec<MatchVect>) -> usize {
    *mv = substruct_match(m(mol), m(patt), &SubstructMatchParameters::default()).unwrap();
    mv.len()
}

fn prop_s(p: &impl RDProps, key: &str) -> String {
    match p.get_prop(key) {
        Some(RDValue::String(s)) => s.clone(),
        Some(v) => v.to_rd_string(),
        None => panic!("missing property {key}"),
    }
}

#[test]
fn test1() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "CH4";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    sln = "CH3CH3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 2);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 0);
    sln = "C[1]H2CH2CH2@1";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 3);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    sln = "C[2]H2CH2CH2@2";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 3);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    sln = "C[200]H2CH2CH2@200";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 3);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    sln = "C[1:foo]H2CH2CH2@1";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 3);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    sln = "C[foo;bar=1;baz=bletch]H4";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 0);
    sln = "CH3CH(CH3)OCH3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 5);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 0);
    assert_eq!(m(&mol).atom_degree(0), 1);
    assert_eq!(m(&mol).atom_degree(1), 3);
    assert_eq!(m(&mol).atom_degree(2), 1);
    assert_eq!(m(&mol).atom_degree(3), 2);
    sln = "CH3CH(CH3)OCH3.Cl";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 6);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 0);
    assert_eq!(m(&mol).atom_degree(0), 1);
    assert_eq!(m(&mol).atom_degree(1), 3);
    assert_eq!(m(&mol).atom_degree(2), 1);
    assert_eq!(m(&mol).atom_degree(3), 2);
    assert_eq!(m(&mol).atom_degree(4), 1);
    assert_eq!(m(&mol).atom_degree(5), 0);
    sln = "H-O-H";
    mol = sln_to_mol(sln, false, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 3);
    assert_eq!(m(&mol).atom(0).atomic_num(), 1);
    assert_eq!(m(&mol).atom(1).atomic_num(), 8);
    assert_eq!(m(&mol).atom(2).atomic_num(), 1);
    sln = "HC(F)(Cl)Br";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 4);
    sln = "CH3C(=O)H";
    mol = sln_to_mol(sln, false, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 7);
}

#[test]
fn test2() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "CH3CH3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 2);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    sln = "CH2-CH2";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 2);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    sln = "CH2=CH2";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 2);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Double
    );
    sln = "CH#CH";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 2);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Triple
    );
    sln = "C[1]H-CH2-CH2-CH2-CH2-CH=@1";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 6);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(1, 2)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 3)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(3, 4)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(4, 5)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(5, 0)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Double
    );
    sln = "C[1]H:CH:CH:CH:CH:CH:@1";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 6);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Aromatic
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(1, 2)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Aromatic
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 3)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Aromatic
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(3, 4)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Aromatic
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(4, 5)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Aromatic
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(5, 0)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Aromatic
    );
}

#[test]
fn test3() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    let mut pval = String::new();
    sln = "C[-]H3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).formal_charge(), -1);
    sln = "C[charge=-]H3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).formal_charge(), -1);
    sln = "C[charge=-1]H3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).formal_charge(), -1);
    sln = "C[CHARGE=-1]H3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).formal_charge(), -1);
    sln = "C[chARgE=-1]H3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).formal_charge(), -1);
    sln = "C[1:-2]H2";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).formal_charge(), -2);
    sln = "C[1:+]H5";
    let sanitize = false;
    mol = sln_to_mol(sln, sanitize, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 6);
    assert_eq!(m(&mol).atom(0).formal_charge(), 1);
    sln = "C[1:+2]H6";
    mol = sln_to_mol(sln, sanitize, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 7);
    assert_eq!(m(&mol).atom(0).formal_charge(), 2);
    sln = "C[1:foo;bar=baz]H4";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert!(m(&mol).atom(0).has_prop("foo"));
    pval = prop_s(m(&mol).atom(0), "foo");
    assert_eq!(pval, "");
    assert!(m(&mol).atom(0).has_prop("bar"));
    pval = prop_s(m(&mol).atom(0), "bar");
    assert_eq!(pval, "baz");
    assert!(!m(&mol).atom(0).has_prop("baz"));
    sln = "H[I=2]-C(-H[I=2])(F)F";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 5);
    assert_eq!(m(&mol).atom(0).isotope(), 2);
    assert_eq!(m(&mol).atom(2).isotope(), 2);
    sln = "C[*]H3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).num_explicit_hs(), 3);
    assert_eq!(m(&mol).atom(0).num_implicit_hs().unwrap(), 0);
    assert!(m(&mol).atom(0).no_implicit());
    assert_eq!(m(&mol).atom(0).num_radical_electrons(), 1);
    sln = "CH[I=2]";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 2);
    sln = "CH4[I=13]";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_none());
}

#[test]
fn test4() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "CH4";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).num_explicit_hs(), 4);
    assert_eq!(m(&mol).atom(0).num_implicit_hs().unwrap(), 0);
    assert!(m(&mol).atom(0).no_implicit());
    sln = "C[-]H3";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    assert_eq!(m(&mol).atom(0).num_explicit_hs(), 3);
    assert_eq!(m(&mol).atom(0).num_implicit_hs().unwrap(), 0);
    assert!(m(&mol).atom(0).no_implicit());
}

#[test]
fn test5() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    let mut pval = String::new();
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "C[charge=+1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(m(&patt).atom(0).num_explicit_hs(), 0);
    assert!(m(&patt).atom(0).no_implicit());
    assert_eq!(m(&patt).atom(0).num_implicit_hs().unwrap(), 0);
    smi = "C[CH2+](C)C";
    let sanitize = false;
    let debugParse = false;
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    smi = "C[CH](C)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) == 0);
    smi = "C[CH3+2](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) == 0);
    smi = "C(=O)OC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    sln = "AnyAny";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    sln = "Any-Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    sln = "Any=Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "Any#Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    sln = "Any:Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    sln = "Any~Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    sln = "C[charge=+1|charge=+2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C[CH2+](C)[CH+2]";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    smi = "C[CH2+](C)[N+2]";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    smi = "C[N+](C)[CH+2]";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "C[charge=+1;HC=2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C[CH2+](CC)[NH+]";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv[0][0].1, 1);
    sln = "Any[charge=+1;HC=2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C[CH2+](CC)[NH+]";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv[0][0].1, 1);
    sln = "Any[charge=+1;!HC=2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C[CH2+](CC)[NH+]";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv[0][0].1, 4);
    smi = "[H]CC[H]";
    mol = smiles(smi, false);
    assert!(mol.is_some());
    sanitize_mol(mm(&mut mol), SANITIZE_ALL).unwrap();
    assert_eq!(m(&mol).num_atoms(), 4);
    sln = "AnyAny";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    sln = "HevHev";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "AnyHev";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    smi = "FC(Cl)(Br)CI";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    sanitize_mol(mm(&mut mol), SANITIZE_ALL).unwrap();
    assert_eq!(m(&mol).num_atoms(), 6);
    sln = "HalC";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 4);
    smi = "CO[2H]";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    sanitize_mol(mm(&mut mol), SANITIZE_ALL).unwrap();
    assert_eq!(m(&mol).num_atoms(), 3);
    sln = "H";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "H[i=1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    sln = "H[i=2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "Het";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
}

#[test]
fn test6() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    let mut pval = String::new();
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "C[rbc=2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert!(!m(&patt).atom(0).has_prop("rbc"));
    smi = "C1CC1";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    smi = "C1C2C1C2";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    smi = "CCC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) == 0);
    sln = "C[rbc=1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C1CC1";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) == 0);
    smi = "CCC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) == 0);
    sln = "C[rbc=3]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C1CC1";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) == 0);
    smi = "C1C2C1C2";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    sln = "C[rBC=2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert!(!m(&patt).atom(0).has_prop("rbc"));
    smi = "C1CC1C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    sln = "C[r]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C1CC1C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    sln = "C[1:rbc=f]CCC@1";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 4);
    smi = "C1CCC1";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    smi = "C1CC2C1CC2";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    sln = "C[1:rbc=f]C[rbc=f]CC@1";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 4);
    smi = "C1CCC1";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) > 0);
    smi = "C1CC2C1CC2";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    sln = "C[1:rbc=f]C[rbc=f]C[rbc=f]C@1";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 4);
    smi = "C1CCC1";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    smi = "C1CC2C1CC2";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert!(sm(&mol, &patt, &mut mv) == 0);
}

#[test]
fn test7() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    let mut pval = String::new();
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "C[charge=+1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C[CH2+](C)C";
    let debugParse = false;
    let sanitize = false;
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "C[charge<1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C[CH2+](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    smi = "C[CH](C)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 4);
    smi = "C[CH3+2](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    sln = "C[charge<=1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    smi = "C[CH2+](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 4);
    smi = "C[CH](C)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 4);
    smi = "C[CH3+2](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    sln = "C[charge>=1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    smi = "C[CH2+](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    smi = "C[CH](C)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    smi = "C[CH3+2](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "C[charge>1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    smi = "C[CH2+](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    smi = "C[CH](C)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    smi = "C[CH3+2](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "C[charge!=0]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    smi = "C[CH2+](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    smi = "C[CH](C)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    smi = "C[CH3+2](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    sln = "C[charge!=1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    smi = "C[CH2+](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    smi = "C[CH](C)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 4);
    smi = "C[CH3+2](C)C";
    mol = smiles(smi, sanitize);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 4);
}

#[test]
fn test8() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    let mut pval = String::new();
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)OC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 4);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[1][0].1, 1);
    assert_eq!(mv[2][0].1, 2);
    assert_eq!(mv[3][0].1, 3);
    sln = "CO[F]C";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 3);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[0][1].1, 2);
    sln = "Any[TBO=4]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv.len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[1][0].1, 3);
    sln = "Any[TAC=2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0][0].1, 2);
    sln = "Any[TAC>2]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv.len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[1][0].1, 3);
    sln = "Any[HC=1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0][0].1, 0);
    sln = "Any[HC=0]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv.len(), 2);
    assert_eq!(mv[0][0].1, 1);
    assert_eq!(mv[1][0].1, 2);
    sln = "Any[HC=F]H3";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0][0].1, 3);
    sln = "Any[HAC>1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv.len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[1][0].1, 2);
    sln = "Any[HAC=1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv.len(), 2);
    assert_eq!(mv[0][0].1, 1);
    assert_eq!(mv[1][0].1, 3);
    sln = "Any[HAC=F]~Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv.len(), 2);
    assert_eq!(mv[0][0].1, 1);
    assert_eq!(mv[1][0].1, 3);
    sln = "AnyAny[TAC=F]Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 3);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0][1].1, 2);
    smi = "CCC=C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    sln = "CC";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    sln = "CC[F]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    sln = "CC[F]H3";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
}

#[test]
fn test9() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    let mut pval = String::new();
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "Any[is=Cl,Br,I]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)Cl";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 2);
    smi = "C(=O)Br";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 2);
    smi = "C(=O)I";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 2);
    smi = "C(=O)C";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    sln = "C[is=C=O]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)OC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 0);
    sln = "Any[is=C,O;charge=-1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)OC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 0);
    smi = "C(=O)[O-]";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 2);
    sln = "Any[Is=C,O;charge=-1]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "[C-](=O)O";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 0);
    sln = "Any[not=N,O]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)OCN";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv.len(), 2);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[1].len(), 1);
    assert_eq!(mv[1][0].1, 3);
}

#[test]
fn test10() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    let mut pval = String::new();
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "Any[is=C[is=Any=O],O]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)OC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 3);
    assert_eq!(mv.len(), 3);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[1].len(), 1);
    assert_eq!(mv[1][0].1, 1);
    assert_eq!(mv[2].len(), 1);
    assert_eq!(mv[2][0].1, 2);
    sln = "Any[not=C[is=Any=O],O]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)OC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 3);
    sln = "C[is=C=O]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    smi = "C(=O)OC";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 0);
    sln = "C[is=O=C*]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 1);
    assert_eq!(mv[0][0].1, 0);
    sln = "C[is=O=C]";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 1);
    assert!(sm(&mol, &patt, &mut mv) == 0);
}

#[test]
fn test12() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "Any-Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    smi = "C=CN";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv.len(), 1);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 1);
    assert_eq!(mv[0][1].1, 2);
    sln = "Any~Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[0][1].1, 1);
    assert_eq!(mv[1].len(), 2);
    assert_eq!(mv[1][0].1, 1);
    assert_eq!(mv[1][1].1, 2);
    sln = "Any-=Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[0][1].1, 1);
    assert_eq!(mv[1].len(), 2);
    assert_eq!(mv[1][0].1, 1);
    assert_eq!(mv[1][1].1, 2);
    sln = "Any=-Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[0][1].1, 1);
    assert_eq!(mv[1].len(), 2);
    assert_eq!(mv[1][0].1, 1);
    assert_eq!(mv[1][1].1, 2);
    sln = "Any-:Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 1);
    assert_eq!(mv[0][1].1, 2);
    sln = "Any-[type=2]Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[0][1].1, 1);
    sln = "Any-[type=2|type=1]Any";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 0);
    assert_eq!(mv[0][1].1, 1);
    assert_eq!(mv[1].len(), 2);
    assert_eq!(mv[1][0].1, 1);
    assert_eq!(mv[1][1].1, 2);
    sln = "O~[r]C";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    smi = "O=CC1COC1";
    mol = smiles(smi, true);
    assert!(mol.is_some());
    assert_eq!(sm(&mol, &patt, &mut mv), 2);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 4);
    assert_eq!(mv[1][0].1, 4);
    sln = "O~[!r]C";
    patt = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(patt.is_some());
    assert_eq!(m(&patt).num_atoms(), 2);
    assert_eq!(sm(&mol, &patt, &mut mv), 1);
    assert_eq!(mv[0].len(), 2);
    assert_eq!(mv[0][0].1, 0);
}

#[test]
fn test13() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "C[1]H2CH2CH2CH2CH2CH2@1";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 6);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(m(&mol).ring_info().atom_rings()[0].len(), 6);
    sln = "C[1]H2CH2(CH2CH2CH2CH2@1)";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 6);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(m(&mol).ring_info().atom_rings()[0].len(), 6);
    sln = "CH2(C[1]H2)CH2(CH2CH2CH2@1)";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 6);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(m(&mol).ring_info().atom_rings()[0].len(), 6);
    sln = "C[1]H2CH2CH2CH2C[2]HCH@1CH2CH2@2";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 8);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 2);
    assert_eq!(m(&mol).ring_info().atom_rings()[0].len(), 4);
    assert_eq!(m(&mol).ring_info().atom_rings()[1].len(), 6);
    sln = "C[1]H2(CH2(CH2(CH2(C[2]H(CH@1(CH2(CH2@2)))))))";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 8);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 2);
    assert_eq!(m(&mol).ring_info().atom_rings()[0].len(), 4);
    assert_eq!(m(&mol).ring_info().atom_rings()[1].len(), 6);
    sln = "C[1](CH2CH2CH(CH2CH2@1)CH2CH2@1)Cl";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 9);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 3);
    assert_eq!(m(&mol).ring_info().atom_rings()[0].len(), 6);
    assert_eq!(m(&mol).ring_info().atom_rings()[1].len(), 6);
    assert_eq!(m(&mol).ring_info().atom_rings()[2].len(), 6);
    sln = "CH2(CH2@1)CH2(C[1]H2)";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_none());
    sln = "CH2(CH2@1)";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_none());
    sln = "C[1]H-CH2-C(@1)=O";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 4);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(1, 2)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 0)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 3)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Double
    );
    sln = "C[1]H-CH2-CH(@1)O";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 4);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(1, 2)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 0)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 3)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    sln = "C[1]H-CH2-C(=@1)O";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 4);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(1, 2)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 0)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Double
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 3)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    sln = "C[1]H-CH2-C(O)=@1";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 4);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(1, 2)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 0)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Double
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 3)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    sln = "C[1]H-CH2-C(O)(=@1)";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 4);
    assert_eq!(m(&mol).ring_info().num_rings().unwrap(), 1);
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(0, 1)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(1, 2)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 0)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Double
    );
    assert_eq!(
        m(&mol)
            .get_bond_between_atoms(2, 3)
            .unwrap()
            .unwrap()
            .bond_type(),
        BondType::Single
    );
}

#[test]
fn test14() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "CH2(C@1H2)CH2(CH2CH2C[1]H2)";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_none());
    sln = "CH2(CH2[1])CH2(CH2CH2CH2@1)";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_none());
}

#[test]
fn test15() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    {
        sln = "CH4<blah>";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert!(m(&mol).has_prop("blah"));
    }
    {
        sln = "CH4<name=methane>";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert!(m(&mol).has_prop(cp::_Name));
        let mut sval = String::new();
        sval = prop_s(m(&mol), cp::_Name);
        assert_eq!(sval, "methane");
    }
    {
        sln = "CH4<blah;foo=\"1\">";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert!(m(&mol).has_prop("blah"));
        assert!(m(&mol).has_prop("foo"));
        let mut sval = String::new();
        sval = prop_s(m(&mol), "foo");
        assert_eq!(sval, "1");
    }
    {
        sln = "CH4<blah;foo=1>";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert!(m(&mol).has_prop("blah"));
        assert!(m(&mol).has_prop("foo"));
        let mut sval = String::new();
        sval = prop_s(m(&mol), "foo");
        assert_eq!(sval, "1");
    }
    {
        sln = "CH4<name=\"methane\";blah;coord2d=(1,0);too.small;test=lots and-lots of special,characters all at once.>";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert!(m(&mol).has_prop("blah"));
        let mut sval = String::new();
        sval = prop_s(m(&mol), cp::_Name);
        assert_eq!(sval, "methane");
        sval = prop_s(m(&mol), "coord2d");
        assert_eq!(sval, "(1,0)");
        assert!(m(&mol).has_prop("too.small"));
        assert!(m(&mol).has_prop("test"));
        sval = prop_s(m(&mol), "test");
        assert_eq!(sval, "lots and-lots of special,characters all at once.");
    }
    {
        sln = "CH4<foo=1><bar=2>";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert!(m(&mol).has_prop("foo"));
        assert!(m(&mol).has_prop("bar"));
        let mut sval = String::new();
        sval = prop_s(m(&mol), "foo");
        assert_eq!(sval, "1");
        sval = prop_s(m(&mol), "bar");
        assert_eq!(sval, "2");
    }
}

#[test]
fn test16() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    {
        sln = "CH3ZCH3{Z:O}";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 3);
        assert_eq!(m(&mol).atom(1).atomic_num(), 8);
    }
    {
        sln = "CH3ZGCH3{Z:O}{G:N}";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 4);
        assert_eq!(m(&mol).atom(1).atomic_num(), 8);
        assert_eq!(m(&mol).atom(2).atomic_num(), 7);
    }
}

#[test]
fn testissue278() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    {
        sln = "Any[IS=C(=O)]";
        mol = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
    {
        sln = "Any[NOT=C(=O)]";
        mol = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
    {
        sln = "Any[IS=C(=O)]";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
    {
        sln = "Any[NOT=C(=O)]";
        mol = sln_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
}

#[test]
fn testissue277() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    {
        sln = "Any[NOT=C,IS=O]";
        mol = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
    {
        sln = "Any[NOT=C;IS=O]";
        mol = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
    {
        sln = "Any[NOT=C&IS=O]";
        mol = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
    {
        sln = "Any[NOT=C,N&IS=O,S]";
        mol = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(mol.is_some());
        assert_eq!(m(&mol).num_atoms(), 1);
        assert!(m(&mol).atom(0).has_query());
    }
    {
        sln = "Hev[!r;NOT=N]";
        patt = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(patt.is_some());
        assert_eq!(m(&patt).num_atoms(), 1);
        assert!(m(&patt).atom(0).has_query());
        smi = "CC1CC1N";
        mol = smiles(smi, true);
        let mut mv: Vec<MatchVect> = Vec::new();
        assert!(mol.is_some());
        assert_eq!(sm(&mol, &patt, &mut mv), 1);
        assert_eq!(mv.len(), 1);
        assert_eq!(mv[0].len(), 1);
        assert_eq!(mv[0][0].1, 0);
    }
    {
        sln = "Hev[!r&NOT=N]";
        patt = sln_query_to_mol(sln, true, 0).unwrap();
        assert!(patt.is_some());
        assert_eq!(m(&patt).num_atoms(), 1);
        assert!(m(&patt).atom(0).has_query());
        smi = "CC1CC1N";
        mol = smiles(smi, true);
        let mut mv: Vec<MatchVect> = Vec::new();
        assert!(mol.is_some());
        assert_eq!(sm(&mol, &patt, &mut mv), 1);
        assert_eq!(mv.len(), 1);
        assert_eq!(mv[0].len(), 1);
        assert_eq!(mv[0][0].1, 0);
    }
    {
        let slns = [
            "HetC(=O)O-[!R]C[8]:Hev(Any[NOT=O,S[TAC=2],C[TAC=4],N[TAC=3]]):Hev:Hev(Any[IS=Hal,C#N,C(F)(F)F,S(=O)=O,C(=O)&NOT=C(=O)OH]):Hev:Hev(Any[NOT=O,S[TAC=2],C[TAC=4],N[TAC=3]]):@8",
            "HetC(=O)O-[!R]C[8]:Hev(Any[IS=Hal,C#N,C(F)(F)F,S(=O)=O,C(=O)&NOT=C(=O)OH]):Hev:Hev(Any[NOT=O,S[TAC=2],C[TAC=4],N[TAC=3]]):Hev:Hev(Any[NOT=O,S[TAC=2],C[TAC=4],N[TAC=3]]):@8",
            "CCH=[!R]C(Any[IS=H,C])Any[IS=C#N,C(=O)&NOT=C(=O)Any[IS=N,O]]",
            "Any[IS=C,O]CH=[!R]C(Any[IS=C(=O),S(=O),C#N,Hal,C(Hal)(Hal)Hal&NOT=C(=O)OH])Any[IS=C(=O),S(=O),C#N&NOT=C(=O)OH]",
            "C[1](C(=O)OC[5]:C:C:C:C:C:@5CH=@1)Any[IS=C(=O),C(=S),S(=O),C#N,Hal,C(Hal)(Hal)Hal&NOT=C(=O)OH]",
            "C[1](=CHOC[5]:C:C:C:C:C:@5C@1=O)Any[IS=C(=O),C(=S),S(=O),C#N,Hal,C(Hal)(Hal)Hal&NOT=C(=O)OH]",
            "C[1]:N:Any(Any[IS=Hal,S(=O)(=O)C]):Any:Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):Any:@1",
            "C[1]:N:Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):Any:Any(Any[IS=Hal,S(=O)(=O)C]):Any:@1",
            "C[1]:N:Any(Any[IS=Hal,S(=O)(=O)C]):Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):Any:Any:@1",
            "C[1]:N:Any(Any[IS=Hal,S(=O)(=O)C]):Any:Any:Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):@1",
            "C[1](Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):N:Any(Any[IS=Hal,S(=O)(=O)C]):Any(Any[NOT=N]):Any(Any[NOT=N]):Any:@1",
            "C[1]:N:Any:Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):Any(Any[IS=Hal,S(=O)(=O)C]):Any:@1",
            "C[1](Any[NOT=N,O]):C(Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):C(Any[IS=Hal,S(=O)(=O)C,C[r](=O)NC]):C(Any[NOT=N,O]):C(Any[NOT=N,O]):C(Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):@1",
            "C[1](Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):C(Any[IS=Hal,S(=O)(=O)C,C[r](=O)NC]):C(Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):C(Any[NOT=N,O]):C(Any[NOT=N,O]):C(Any[NOT=N,O]):@1",
            "N(Any[IS=H,C[TAC=4]&NOT=C[TAC=4]-[R]C[TAC=4]N])(Any[IS=H,C[TAC=4]&NOT=C[TAC=4]-[R]C[TAC=4]N])(Any[IS=H,C[TAC=4]&NOT=C[TAC=4]-[R]C[TAC=4]N])<max=1>",
            "Any[IS=H,C&NOT=C=O]N[!r](Any[IS=H,C&NOT=C=O])C(=O)C<max=2>",
            "Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]",
            "Hev[!r&NOT=C=O,S(=O)(=O),N*S(=O),N*(C=O)]Hev[!r&NOT=C=O,S(=O)(=O),N*S(=O),N*(C=O)]Hev[!r&NOT=C=O,S(=O)(=O),N*S(=O),N*(C=O)]Hev[!r&NOT=C=O,S(=O)(=O),N*S(=O),N*(C=O)]Hev[!r&NOT=C=O,S(=O)(=O)]Any[IS=CH3,OH,NH2,N(CH3)CH3]",
            "EOF",
        ];
        for sln in slns.iter().take_while(|s| **s != "EOF") {
            mol = sln_query_to_mol(sln, true, 0).unwrap();
            assert!(mol.is_some());
        }
    }
}

#[test]
fn test17() {
    let mut mol: Option<RWMol>;
    let mut patt: Option<RWMol>;
    let mut sln = "";
    let mut smi = "";
    let mut mv: Vec<MatchVect> = Vec::new();
    sln = "CH4 \t";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    sln = "CH4\t";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    assert_eq!(m(&mol).num_atoms(), 1);
    sln = "CH4\tfff";
    mol = sln_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_none());
    sln = "C[charge=+1] \t";
    mol = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_some());
    sln = "C[charge=+1] \tfoo";
    mol = sln_query_to_mol(sln, true, 0).unwrap();
    assert!(mol.is_none());
}
