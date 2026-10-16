//! SynthonSpaceSearch; the substructure cases of upstream
//! `SynthonSpaceSearch/Wrap/testSynthonSpaceSearch.py`.

use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::substruct::SubstructMatchParameters;
use rdoxide::graphmol::synthon_space_search::{SynthonSpace, SynthonSpaceSearchParams};
use rdoxide::rdgeneral::RDProps;

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

#[test]
fn shapes_round_trip() {
    if !have_data() {
        return;
    }
    let mut space = SynthonSpace::new();
    space
        .read_db_file(&format!("{DATA}amide_space_shapes.spc"), 1)
        .unwrap();
    assert!(space.num_synthons_with_shapes() > 0);
    let tmp = std::env::temp_dir().join("rdoxide_shapes_rt.spc");
    space.write_db_file(tmp.to_str().unwrap()).unwrap();
    let a = std::fs::read(format!("{DATA}amide_space_shapes.spc")).unwrap();
    let b = std::fs::read(&tmp).unwrap();
    std::fs::remove_file(&tmp).ok();
    // identical apart from the minor version byte of each molecule pickle
    // (the file was written with an older RDKit)
    assert_eq!(a.len(), b.len());
    let diffs = a.iter().zip(&b).filter(|(x, y)| x != y).count();
    assert!(diffs <= 40, "{diffs} bytes differ");
}

const SHAPE_QUERY: &str = "O=C(c1ccccc1)N1CCCC1 |(0.0443291,-1.81486,-1.76886;0.0506321,-0.858174,-0.921491;1.37975,-0.430412,-0.483603;2.18964,-1.35506,0.144714;3.47088,-1.00454,0.585539;3.93803,0.297573,0.388032;3.1267,1.22739,-0.242406;1.85597,0.849751,-0.670531;-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;-2.73583,1.04666,0.696614;-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401)|";

#[test]
fn shape_search_small() {
    if !have_data() {
        return;
    }
    let mut space = SynthonSpace::new();
    space
        .read_db_file(&format!("{DATA}amide_space_shapes.spc"), 1)
        .unwrap();
    let params = SynthonSpaceSearchParams {
        similarity_cutoff: 0.8,
        frag_similarity_adjuster: 0.2,
        approx_similarity_adjuster: 0.2,
        num_conformers: 100,
        conf_rms_threshold: 0.5,
        time_out: 0,
        random_seed: 0xdac,
        best_hit: true,
        ..Default::default()
    };
    let query = rdoxide::graphmol::smiles_parse::mol_from_smiles(
        SHAPE_QUERY,
        &rdoxide::graphmol::smiles_parse::SmilesParserParams::default(),
    )
    .unwrap()
    .unwrap();
    let res = space.shape_search(&query, &params).unwrap();
    let names: Vec<String> = res.hit_molecules().iter().map(|m| m.name()).collect();
    assert_eq!(res.hit_molecules().len(), 3, "{names:?}");
    for m in res.hit_molecules() {
        let sim = m.get_prop("Similarity").unwrap().as_double().unwrap();
        assert!(sim > 0.9, "{sim}");
        assert!(
            ["1-1;2-1;amide-1", "1-2;2-1;amide-1", "1-3;2-1;amide-1"].contains(&m.name().as_str())
        );
    }
}
