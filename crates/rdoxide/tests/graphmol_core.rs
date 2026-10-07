//! Core molecule-graph tests, ported from RDKit `Code/GraphMol/test1.cpp`
//! (testBookmarks, testDegree, testAddAtomWithConf, ...) plus additional
//! checks of RWMol editing semantics. BSD-3-Clause, RDKit contributors.

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use rdoxide::graphmol::{Atom, BondStereo, BondType, Conformer, RWMol, ValenceType};
use rdoxide::rdgeneral::RDProps;

fn chain(nums: &[u32]) -> RWMol {
    let mut m = RWMol::new();
    for &n in nums {
        m.add_atom(Atom::new(n), true);
    }
    for i in 1..nums.len() as u32 {
        m.add_bond(i - 1, i, BondType::Single).unwrap();
    }
    m
}

#[test]
fn bookmarks() {
    let mut m = chain(&[6, 6, 6]);
    m.set_atom_bookmark(1, 666);
    assert_eq!(m.get_atom_with_bookmark(666).unwrap(), 1);
    m.clear_atom_bookmark(666);
    assert!(m.get_atom_with_bookmark(666).is_err());

    m.set_atom_bookmark(1, 666);
    m.set_atom_bookmark(0, 666);
    assert_eq!(m.get_atom_with_bookmark(666).unwrap(), 1);
    m.clear_atom_bookmark_for(666, 1);
    assert_eq!(m.get_atom_with_bookmark(666).unwrap(), 0);
    m.clear_atom_bookmark_for(666, 0);
    assert!(m.get_atom_with_bookmark(666).is_err());
    m.clear_atom_bookmark(777);

    m.set_bond_bookmark(0, 23);
    assert_eq!(m.get_bond_with_bookmark(23).unwrap(), 0);
    m.clear_bond_bookmark(23);
    assert!(m.get_bond_with_bookmark(23).is_err());
}

#[test]
fn degree_and_hs() {
    let mut m = chain(&[6, 6, 6, 6]);
    m.add_bond(0, 2, BondType::Single).unwrap();
    m.update_property_cache(true).unwrap();
    assert_eq!(m.atom_degree(0), 2);
    assert_eq!(m.total_num_hs(0, false).unwrap(), 2);
    assert_eq!(m.total_degree(0).unwrap(), 4);
    assert_eq!(m.atom_degree(2), 3);
    assert_eq!(m.total_num_hs(2, false).unwrap(), 1);
    assert_eq!(m.num_atoms_with_hs().unwrap(), 4 + 2 + 2 + 1 + 3);

    let mut m = chain(&[6, 6, 6, 6]);
    m.add_bond(0, 2, BondType::Single).unwrap();
    let h = m.add_atom(Atom::new(1), true);
    m.add_bond(0, h, BondType::Single).unwrap();
    m.update_property_cache(true).unwrap();
    assert_eq!(m.atom_degree(0), 3);
    assert_eq!(m.total_num_hs(0, false).unwrap(), 1);
    assert_eq!(m.total_num_hs(0, true).unwrap(), 2);
    assert_eq!(m.total_degree(0).unwrap(), 4);
}

#[test]
fn add_atom_with_conf() {
    let mut m = chain(&[6, 6]);
    m.add_conformer(Conformer::new(2), true).unwrap();
    m.add_atom(Atom::new(6), true);
    assert_eq!(m.conformer(-1).unwrap().num_atoms(), m.num_atoms());
    m.add_new_atom(true);
    assert_eq!(m.conformer(-1).unwrap().num_atoms(), 4);
    assert!(m.add_conformer(Conformer::new(2), true).is_err());
}

#[test]
fn valences() {
    // charged nitrogen
    let mut m = chain(&[7]);
    m.atom_mut(0).set_formal_charge(1);
    m.update_property_cache(true).unwrap();
    assert_eq!(m.atom(0).valence(ValenceType::Implicit).unwrap(), 4);
    // pentavalent carbon fails in strict mode, with an AtomValence error
    let mut m = chain(&[6, 6, 6, 6, 6, 6]);
    for i in 2..6 {
        let _ = m.add_bond(1, i, BondType::Single);
    }
    assert!(m.update_property_cache(true).is_err());
    assert!(m.update_property_cache(false).is_ok());
    assert!(rdoxide::graphmol::atom::has_valence_violation(&m, 1));
    assert!(!rdoxide::graphmol::atom::has_valence_violation(&m, 0));
}

#[test]
fn remove_atom_renumbers() {
    // C0-C1-O2-N3, bond stereo atoms referencing the removed atom are cleared
    let mut m = chain(&[6, 6, 8, 7]);
    m.add_conformer(Conformer::new(4), true).unwrap();
    for i in 0..4 {
        m.conformer_mut(-1)
            .unwrap()
            .set_atom_pos(i, Point3D::new(i as f64, 0.0, 0.0))
            .unwrap();
    }
    m.set_atom_bookmark(3, 1);
    m.remove_atom(1).unwrap();
    assert_eq!(m.num_atoms(), 3);
    assert_eq!(m.num_bonds(), 1);
    let b = m.bond(0);
    assert_eq!((b.begin_atom_idx(), b.end_atom_idx(), b.idx()), (1, 2, 0));
    assert_eq!(m.atom(2).atomic_num(), 7);
    assert_eq!(m.atom(2).idx(), 2);
    assert_eq!(m.get_atom_with_bookmark(1).unwrap(), 2);
    let xs: Vec<f64> = m
        .conformer(-1)
        .unwrap()
        .positions()
        .iter()
        .map(|p| p.x)
        .collect();
    assert_eq!(xs, vec![0.0, 2.0, 3.0]);
    assert_eq!(m.atom_neighbors(1).collect::<Vec<_>>(), vec![2]);
}

#[test]
fn remove_bond_renumbers_and_clears_stereo() {
    // C0=C1 with stereo atoms 2 and 3; removing 0-2 clears the stereo
    let mut m = RWMol::new();
    for _ in 0..4 {
        m.add_atom(Atom::new(6), true);
    }
    m.add_bond(0, 1, BondType::Double).unwrap();
    m.add_bond(0, 2, BondType::Single).unwrap();
    m.add_bond(1, 3, BondType::Single).unwrap();
    m.set_bond_stereo_atoms(0, 2, 3).unwrap();
    m.bond_mut(0).set_stereo(BondStereo::StereoCis).unwrap();
    m.remove_bond(0, 2).unwrap();
    assert_eq!(m.num_bonds(), 2);
    assert_eq!(m.bond(0).stereo(), BondStereo::StereoNone);
    assert!(m.bond(0).stereo_atoms().is_empty());
    assert_eq!(m.bond(1).idx(), 1);
    assert_eq!(
        (m.bond(1).begin_atom_idx(), m.bond(1).end_atom_idx()),
        (1, 3)
    );
    assert_eq!(m.bond_idx_between(3, 1), Some(1));
    assert_eq!(m.bond_idx_between(0, 2), None);
}

#[test]
fn neighbor_order_is_insertion_order() {
    let mut m = RWMol::new();
    for _ in 0..4 {
        m.add_atom(Atom::new(6), false);
    }
    m.add_bond(0, 3, BondType::Single).unwrap();
    m.add_bond(0, 1, BondType::Single).unwrap();
    m.add_bond(2, 0, BondType::Single).unwrap();
    assert_eq!(m.atom_neighbors(0).collect::<Vec<_>>(), vec![3, 1, 2]);
    assert_eq!(m.atom_bond_indices(0).collect::<Vec<_>>(), vec![0, 1, 2]);
    // perturbation order from Atom::getPerturbationOrder docs
    assert_eq!(m.perturbation_order(0, &[1, 0, 2]).unwrap(), 1);
}

#[test]
fn batch_edit() {
    let mut m = chain(&[6, 7, 8, 9, 17]);
    m.add_conformer(Conformer::new(5), true).unwrap();
    m.begin_batch_edit().unwrap();
    assert!(m.begin_batch_edit().is_err());
    m.remove_atom(1).unwrap();
    m.remove_atom(3).unwrap();
    // nothing happens until commit
    assert_eq!(m.num_atoms(), 5);
    m.commit_batch_edit().unwrap();
    let nums: Vec<u32> = m.atoms().iter().map(|a| a.atomic_num()).collect();
    assert_eq!(nums, vec![6, 8, 17]);
    assert_eq!(m.num_bonds(), 0);
    assert_eq!(m.conformer(-1).unwrap().num_atoms(), 3);

    let mut m = chain(&[6, 7, 8]);
    m.begin_batch_edit().unwrap();
    m.remove_bond(0, 1).unwrap();
    m.rollback_batch_edit();
    m.commit_batch_edit().unwrap();
    assert_eq!(m.num_bonds(), 2);
}

#[test]
fn insert_mol_and_stereo_groups() {
    let mut a = chain(&[6, 6]);
    a.set_stereo_groups(vec![
        StereoGroup::new(StereoGroupType::Absolute, vec![0], vec![], 0).unwrap(),
    ])
    .unwrap();
    let mut b = chain(&[8, 7]);
    b.set_stereo_groups(vec![
        StereoGroup::new(StereoGroupType::Absolute, vec![1], vec![], 0).unwrap(),
        StereoGroup::new(StereoGroupType::Or, vec![0], vec![], 3).unwrap(),
    ])
    .unwrap();
    a.insert_mol(&b).unwrap();
    assert_eq!(a.num_atoms(), 4);
    assert_eq!(
        (a.bond(1).begin_atom_idx(), a.bond(1).end_atom_idx()),
        (2, 3)
    );
    let g = a.stereo_groups();
    assert_eq!(g.len(), 2);
    assert_eq!(g[0].group_type(), StereoGroupType::Or);
    assert_eq!(g[0].atoms(), &[2]);
    assert_eq!(g[1].group_type(), StereoGroupType::Absolute);
    assert_eq!(g[1].atoms(), &[0, 3]);
    assert!(StereoGroup::new(StereoGroupType::And, vec![1, 1], vec![], 0).is_err());
}

#[test]
fn mol_props_and_copy() {
    let mut m = chain(&[6, 8]);
    m.set_name("ethanol-ish");
    m.set_prop("foo", 3, true).unwrap();
    m.atom_mut(0).set_prop("p", "x", false).unwrap();
    let full = rdoxide::graphmol::ROMol::copy_from(&m, false, -1);
    assert_eq!(full.name(), "ethanol-ish");
    let quick = rdoxide::graphmol::ROMol::copy_from(&m, true, -1);
    assert_eq!(quick.name(), "");
    assert!(quick.atom(0).has_prop("p"));
    m.clear_computed_props(true);
    assert!(!m.has_prop("foo"));
}

#[test]
fn errors() {
    let mut m = chain(&[6, 6]);
    assert!(m.get_atom(5).is_err());
    assert!(m.add_bond(0, 0, BondType::Single).is_err());
    assert!(m.add_bond(0, 1, BondType::Single).is_err());
    assert!(m.get_bond_between_atoms(0, 7).is_err());
    assert!(m.conformer(-1).is_err());
    let a = rdoxide::graphmol::Atom::from_symbol("Xx");
    assert!(a.is_err());
}
