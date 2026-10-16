//! Fragment catalogs (cases from upstream `FragCatalog/Wrap/rough_test.py`
//! and `test1.cpp`).

use rdoxide::catalogs::CatalogEntry;
use rdoxide::graphmol::frag_catalog::{
    FragCatGenerator, FragCatParams, FragCatalog, FragFPGenerator,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::utils::first_thousand_primes;

const FGROUPS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/rdkit/Code/GraphMol/FragCatalog/test_data/funcGroups.txt"
);

fn parse(smi: &str) -> rdoxide::graphmol::RWMol {
    mol_from_smiles(smi, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn primes() {
    let p = first_thousand_primes();
    assert_eq!(&p[..6], [2, 3, 5, 7, 11, 13]);
    assert_eq!(p[999], 7919);
}

#[test]
fn params_text_round_trip() {
    if !std::path::Path::new(FGROUPS).exists() {
        return; // the upstream checkout is not available
    }
    let params = FragCatParams::new(1, 6, FGROUPS, 1e-8).unwrap();
    assert_eq!(params.num_func_groups(), 15);
    let text = params.serialize().unwrap();
    assert!(text.starts_with("1 6 1e-08\n15\n"));
    let back = FragCatParams::from_text(&text).unwrap();
    assert_eq!(back.num_func_groups(), 15);
    assert_eq!(back.serialize().unwrap(), text);
    assert!(FragCatParams::new(3, 2, FGROUPS, 1e-8).is_err());
    assert!(FragCatParams::new(1, 2, "/no/such/file", 1e-8).is_err());
}

#[test]
fn fragments_fingerprints_and_pickles() {
    if !std::path::Path::new(FGROUPS).exists() {
        return;
    }
    let mut fcat = FragCatalog::new(FragCatParams::new(2, 2, FGROUPS, 1e-8).unwrap());
    let nent = FragCatGenerator
        .add_frags_from_mol(&parse("Cc1ccccc1"), &mut fcat)
        .unwrap();
    assert_eq!(nent, fcat.num_entries());
    assert_eq!(fcat.fp_length(), 2);
    assert_eq!(
        fcat.entry_with_bit_id(0).unwrap().unwrap().description(),
        "ccC"
    );

    let fp = FragFPGenerator
        .fp_for_mol(&parse("Cc1ccccc1"), &fcat)
        .unwrap();
    assert_eq!(rdoxide::datastructs::BitVect::on_bits(&fp), [0, 1]);
    let fp = FragFPGenerator
        .fp_for_mol(&parse("c1ccccc1-c1ccccc1"), &fcat)
        .unwrap();
    assert_eq!(rdoxide::datastructs::BitVect::on_bits(&fp), [1]);

    let fcat2 = FragCatalog::from_pickle(&fcat.serialize().unwrap()).unwrap();
    assert_eq!(fcat2.num_entries(), fcat.num_entries());
    assert_eq!(fcat2.fp_length(), 2);
    let fp2 = FragFPGenerator
        .fp_for_mol(&parse("Cc1ccccc1"), &fcat2)
        .unwrap();
    assert_eq!(rdoxide::datastructs::BitVect::on_bits(&fp2), [0, 1]);
    assert!(FragCatalog::from_pickle(b"junk").is_err());
}
