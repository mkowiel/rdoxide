//! Filter catalogs (`graphmol::filter_catalog`).

use std::sync::Arc;

use rdoxide::graphmol::filter_catalog::*;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::ROMol;

fn mol(smi: &str) -> ROMol {
    mol_from_smiles(smi, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn builtin_catalog_sizes() {
    for (flag, n) in [
        (catalogs::PAINS_A, 16),
        (catalogs::PAINS_B, 55),
        (catalogs::PAINS_C, 409),
        (catalogs::PAINS, 480),
        (catalogs::BRENK, 105),
        (catalogs::NIH, 180),
        (catalogs::ZINC, 50),
        (catalogs::CHEMBL, 770),
    ] {
        assert_eq!(FilterCatalog::from_catalogs(flag).unwrap().num_entries(), n, "{flag}");
    }
}

#[test]
fn pains_match_and_properties() {
    let c = FilterCatalog::from_catalogs(catalogs::PAINS_A).unwrap();
    let m = mol("O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2");
    let e = c.first_match(&m).unwrap().expect("a match");
    assert_eq!(e.description(), "hzone_phenol_A(479)");
    assert_eq!(e.prop("Scope").unwrap(), "PAINS filters (family A)");
    assert!(c.has_match(&m).unwrap());
    assert!(!c.has_match(&mol("CCCC")).unwrap());
}

#[test]
fn match_operators() {
    let carbon: SharedMatcher = Arc::new(SmartsMatcher::from_smarts("Carbon", "[#6]", 1, 2));
    let oxygen: SharedMatcher = Arc::new(SmartsMatcher::from_smarts("Oxygen", "[#8]", 1, 1));
    let or = Or::new(carbon.clone(), oxygen.clone());
    assert_eq!(or.name(), "(Carbon Or Oxygen)");
    assert!(or.has_match(&mol("C")).unwrap());
    assert!(!or.has_match(&mol("CCC")).unwrap());
    let not = Not::new(carbon.clone());
    assert!(not.has_match(&mol("CCC")).unwrap());
    let and = And::new(carbon, oxygen);
    assert!(and.has_match(&mol("CO")).unwrap());
    let mut ms = Vec::new();
    assert!(and.get_matches(&mol("CO"), &mut ms).unwrap());
    assert!(!ms.is_empty());
}

#[test]
fn runner_is_order_and_thread_stable() {
    let c = FilterCatalog::from_catalogs(catalogs::BRENK).unwrap();
    let smiles: Vec<String> = ["CC(=O)Cl", "CCCC", "not a smiles", "C=O"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let names = |threads| -> Vec<Vec<String>> {
        run_filter_catalog(&c, &smiles, threads)
            .into_iter()
            .map(|r| r.unwrap().iter().map(|e| e.description()).collect())
            .collect()
    };
    let one = names(1);
    assert_eq!(one, names(4));
    assert_eq!(one[2], vec!["no valid RDKit molecule"]);
    assert!(one[0].contains(&"acid_halide".to_string()));
    assert!(one[1].is_empty() || !one[1].contains(&"acid_halide".to_string()));
}

#[test]
fn functional_group_hierarchy_reports_most_specific() {
    let fc = functional_group_hierarchy();
    let matches = fc.filter_matches(&mol("c1ccccc1Cl")).unwrap();
    let names: Vec<String> = matches.iter().map(|m| m.filter_match.name()).collect();
    assert_eq!(names, ["Halogen.Aromatic", "Halogen.NotFluorine.Aromatic"]);
}
