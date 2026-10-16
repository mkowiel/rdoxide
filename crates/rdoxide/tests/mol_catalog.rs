//! MolCatalog and the hierarchical catalog framework (upstream
//! `MolCatalog/test1.cpp`).

use rdoxide::catalogs::CatalogEntry;
use rdoxide::graphmol::mol_catalog::{MolCatalog, MolCatalogEntry, MolCatalogParams};
use rdoxide::graphmol::smiles_parse::write::mol_to_smiles;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn entry(smi: &str) -> MolCatalogEntry {
    let m = mol_from_smiles(smi, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    let mut e = MolCatalogEntry::new(m);
    e.set_description(smi);
    e
}

fn catalog() -> MolCatalog {
    let mut cat = MolCatalog::new(MolCatalogParams::default());
    for smi in ["C1CCC1OC", "C1CCC1", "C"] {
        cat.add_entry(entry(smi), true);
    }
    cat.add_edge(0, 1).unwrap();
    cat.add_edge(0, 2).unwrap();
    cat.add_edge(0, 1).unwrap(); // duplicates are ignored
    cat.add_edge(1, 2).unwrap();
    cat
}

#[test]
fn bit_ids_and_edges() {
    let mut cat = catalog();
    assert_eq!(cat.num_entries(), 3);
    assert_eq!(cat.fp_length(), 3);
    assert_eq!(cat.id_of_entry_with_bit_id(2).unwrap(), 2);
    assert_eq!(
        cat.entry_with_bit_id(1).unwrap().unwrap().description(),
        "C1CCC1"
    );
    assert_eq!(cat.down_entry_list(0).unwrap(), [1, 2]);
    assert!(cat.add_edge(0, 3).is_err());
    assert!(cat.entry_with_idx(3).is_err());
    assert!(cat.id_of_entry_with_bit_id(3).is_err());
    assert_eq!(cat.entries_of_order(0).unwrap(), [0, 1, 2]);
    assert!(cat.entries_of_order(1).is_err());
}

#[test]
fn pickle_round_trip() {
    let cat = catalog();
    let cat2 = MolCatalog::from_pickle(&cat.serialize().unwrap()).unwrap();
    assert_eq!(cat2.num_entries(), 3);
    assert_eq!(cat2.fp_length(), 3);
    assert_eq!(cat2.down_entry_list(0).unwrap(), [1, 2]);
    assert_eq!(cat2.down_entry_list(1).unwrap(), [2]);
    for i in 0..3 {
        let (a, b) = (
            cat.entry_with_idx(i).unwrap(),
            cat2.entry_with_idx(i).unwrap(),
        );
        assert_eq!(a.description(), b.description());
        assert_eq!(a.bit_id(), b.bit_id());
        let smiles = |e: &MolCatalogEntry| {
            let mut m = e.mol().unwrap().clone();
            mol_to_smiles(&mut m, &Default::default()).unwrap()
        };
        assert_eq!(smiles(a), smiles(b));
    }
    assert!(MolCatalog::from_pickle(&cat.serialize().unwrap()[..20]).is_err());
}
