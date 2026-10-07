//! removeHs and legacy stereo perception; expected values from RDKit
//! 2026.09.1 (cf. upstream catch_chirality.cpp / molopstest.cpp).

use rdoxide::graphmol::chirality::legacy::assign_stereochemistry;
use rdoxide::graphmol::mol_ops::{RemoveHsParameters, remove_hs};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{BondStereo, ChiralType, RWMol};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn parse(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn raw(s: &str) -> RWMol {
    let p = SmilesParserParams {
        sanitize: false,
        remove_hs: false,
        ..Default::default()
    };
    mol_from_smiles(s, &p).unwrap().unwrap()
}

fn cip(m: &RWMol) -> Vec<Option<String>> {
    m.atoms()
        .iter()
        .map(|a| match a.get_prop("_CIPCode") {
            Some(RDValue::String(s)) => Some(s.clone()),
            _ => None,
        })
        .collect()
}

#[test]
fn cip_codes() {
    let r = Some("R".to_string());
    let s = Some("S".to_string());
    assert_eq!(cip(&parse("C[C@H](O)F")), vec![None, r.clone(), None, None]);
    assert_eq!(
        cip(&parse("C[C@@H](O)F")),
        vec![None, s.clone(), None, None]
    );
    assert_eq!(cip(&parse("F[C@H](Cl)Br")), vec![None, r, None, None]);
    assert_eq!(cip(&parse("[H][C@](F)(Cl)Br")), vec![s, None, None, None]);
}

#[test]
fn spurious_tag_removed() {
    let m = parse("C[C@H](C)O");
    assert_eq!(m.atom(1).chiral_tag(), ChiralType::Unspecified);
    assert_eq!(m.atom(1).num_explicit_hs(), 0);
}

#[test]
fn ring_stereo() {
    let m = parse("C[C@H]1CC[C@@H](C)CC1");
    assert_eq!(
        m.atom(1).get_prop("_ringStereoAtoms"),
        Some(&RDValue::VecInt(vec![5]))
    );
    assert_eq!(
        m.atom(4).get_prop("_ringStereoAtoms"),
        Some(&RDValue::VecInt(vec![2]))
    );
    assert_eq!(m.atom(1).chiral_tag(), ChiralType::TetrahedralCW);
}

#[test]
fn double_bonds() {
    let m = parse("F/C=C/F");
    assert_eq!(m.bond(1).stereo(), BondStereo::StereoE);
    assert_eq!(m.bond(1).stereo_atoms(), &[0, 3]);
    assert_eq!(parse("F/C=C\\F").bond(1).stereo(), BondStereo::StereoZ);
    assert_eq!(
        parse("C1/C=C/CCCC1").bond(1).stereo(),
        BondStereo::StereoNone
    );
    assert_eq!(parse("C1/C=C/CCCCC1").bond(1).stereo(), BondStereo::StereoE);
    assert_eq!(
        parse("F/C(\\Cl)=C/F").bond(2).stereo(),
        BondStereo::StereoNone
    );
}

#[test]
fn stereochem_done_and_force() {
    let mut m = parse("CC(O)F");
    assert!(m.atom(1).has_prop("_ChiralityPossible"));
    m.atom_mut(1).set_chiral_tag(ChiralType::TetrahedralCW);
    assign_stereochemistry(&mut m, false, false, false).unwrap();
    assert!(!m.atom(1).has_prop("_CIPCode"));
    assign_stereochemistry(&mut m, false, true, false).unwrap();
    assert_eq!(cip(&m)[1].as_deref(), Some("S"));
}

#[test]
fn remove_hs_basics() {
    let mut m = raw("[H]C([H])([H])[C@@H](O)F");
    remove_hs(&mut m, &RemoveHsParameters::default(), true).unwrap();
    assert_eq!(m.num_atoms(), 4);
    assert_eq!(m.atom(1).chiral_tag(), ChiralType::TetrahedralCW);

    for s in ["[H][H]", "*[H]", "C[H-]C", "[H]/N=C/C", "[2H]C"] {
        let mut m = raw(s);
        let n = m.num_atoms();
        remove_hs(&mut m, &RemoveHsParameters::default(), false).unwrap();
        assert_eq!(m.num_atoms(), n, "{s}");
    }

    let mut m = raw("[H]C([H])([H])[H]");
    let ps = RemoveHsParameters {
        update_explicit_count: true,
        ..Default::default()
    };
    remove_hs(&mut m, &ps, true).unwrap();
    assert_eq!(m.atom(0).num_explicit_hs(), 4);
}

#[test]
fn remove_and_track_isotopes() {
    let mut m = raw("[2H]C([2H])O");
    let ps = RemoveHsParameters {
        remove_and_track_isotopes: true,
        ..Default::default()
    };
    remove_hs(&mut m, &ps, true).unwrap();
    assert_eq!(m.num_atoms(), 2);
    assert_eq!(
        m.atom(0).get_prop("_isotopicHs"),
        Some(&RDValue::VecUnsignedInt(vec![2, 2]))
    );
}
