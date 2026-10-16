//! SynthonSpaceSearch; the substructure cases of upstream
//! `SynthonSpaceSearch/Wrap/testSynthonSpaceSearch.py`.

use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::substruct::SubstructMatchParameters;
use rdoxide::graphmol::synthon_space_search::{SynthonSpace, SynthonSpaceSearchParams};

const DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/rdkit/Code/GraphMol/SynthonSpaceSearch/data/"
);

/// The upstream test data is not present in every checkout.
fn have_data() -> bool {
    std::path::Path::new(DATA).exists()
}

fn space_a() -> SynthonSpace {
    let mut s = SynthonSpace::new();
    s.read_db_file(&format!("{DATA}idorsia_toy_space_a.spc"), 1)
        .unwrap();
    s
}

fn query() -> rdoxide::graphmol::RWMol {
    mol_from_smarts("c1ccccc1C(=O)N1CCCC1", &SmartsParserParams::default())
        .unwrap()
        .unwrap()
}

#[test]
fn substruct_search() {
    if !have_data() {
        return;
    }
    let space = space_a();
    let params = SynthonSpaceSearchParams {
        max_hits: 10,
        ..Default::default()
    };
    let res = space
        .substructure_search(&query(), &SubstructMatchParameters::default(), &params)
        .unwrap();
    assert_eq!(res.hit_molecules().len(), 10);
}

#[test]
fn cutoffs() {
    if !have_data() {
        return;
    }
    let space = space_a();
    let mut params = SynthonSpaceSearchParams {
        max_hits: -1,
        ..Default::default()
    };
    let ps = SubstructMatchParameters::default();
    let res = space.substructure_search(&query(), &ps, &params).unwrap();
    assert_eq!(res.hit_molecules().len(), 220);
    params.max_hit_heavy_atoms = 30;
    params.min_hit_heavy_atoms = 25;
    let res = space.substructure_search(&query(), &ps, &params).unwrap();
    assert_eq!(res.hit_molecules().len(), 141);
    params.max_hit_heavy_atoms = -1;
    params.min_hit_heavy_atoms = 0;
    params.max_hit_mol_wt = 450.0;
    params.min_hit_mol_wt = 350.0;
    let res = space.substructure_search(&query(), &ps, &params).unwrap();
    assert_eq!(res.hit_molecules().len(), 185);
}

#[test]
fn text_file_read() {
    if !have_data() {
        return;
    }
    let mut space = SynthonSpace::new();
    space
        .read_text_file(&format!("{DATA}amide_space.txt"))
        .unwrap();
    assert_eq!(space.num_reactions(), 1);
    assert_eq!(space.num_products(), 12);
    let params = SynthonSpaceSearchParams {
        max_hits: -1,
        ..Default::default()
    };
    let res = space
        .substructure_search(&query(), &SubstructMatchParameters::default(), &params)
        .unwrap();
    assert!(!res.hit_molecules().is_empty());
}

#[test]
fn extended_query() {
    if !have_data() {
        return;
    }
    use rdoxide::graphmol::generalized_substruct::create_extended_query_mol;
    let mut space = SynthonSpace::new();
    space
        .read_text_file(&format!("{DATA}extended_query.csv"))
        .unwrap();
    assert_eq!(space.num_reactions(), 7);
    let m = mol_from_smarts(
        "[#6]-*.c1nc2cccnc2n1 |m:1:3.10|",
        &SmartsParserParams::default(),
    )
    .unwrap()
    .unwrap();
    let xq = create_extended_query_mol(
        &m,
        true,
        true,
        false,
        &rdoxide::graphmol::adjust_query::AdjustQueryParameters::default(),
    )
    .unwrap();
    let res = space
        .substructure_search_xq(
            &xq,
            &SubstructMatchParameters::default(),
            &SynthonSpaceSearchParams::default(),
        )
        .unwrap();
    assert_eq!(res.hit_molecules().len(), 12);
}
