//! `SubstanceGroup` and its bookkeeping inside molecules; expectations from
//! RDKit 2026.09.1.

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::mol_ops::add_hs::{AddHsParameters, RemoveHsParameters, add_hs, remove_hs};
use rdoxide::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::subset::{SubsetOptions, copy_mol_subset};
use rdoxide::graphmol::substance_group::{
    SGroupBondType, SubstanceGroup, add_substance_group, checks,
};
use rdoxide::graphmol::{ROMol, RWMol};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn parse(smiles: &str) -> RWMol {
    mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

/// `CCOCC` with an SRU over atoms 1-2 whose crossing bonds are 0 and 3.
fn sru_mol(smiles: &str) -> RWMol {
    let mut mol = parse(smiles);
    let mut group = SubstanceGroup::new("SRU");
    group.set_atoms(&mol, vec![1, 2]).unwrap();
    group.set_bonds(&mol, vec![0, 3]).unwrap();
    add_substance_group(&mut mol, group);
    mol
}

fn summary(mol: &ROMol) -> Vec<(Vec<u32>, Vec<u32>)> {
    mol.substance_groups()
        .iter()
        .map(|g| (g.atoms().to_vec(), g.bonds().to_vec()))
        .collect()
}

#[test]
fn groups_validate_indices_against_the_molecule() {
    let mol = parse("CCOCC");
    let mut group = SubstanceGroup::new("SRU");
    assert!(group.add_atom_with_idx(&mol, 9).is_err());
    assert!(group.add_bond_with_idx(&mol, 9).is_err());
    assert!(group.set_atoms(&mol, vec![0, 20]).is_err());
    group.add_atom_with_idx(&mol, 1).unwrap();
    assert!(group.add_parent_atom_with_idx(&mol, 0).is_err());
    group.add_parent_atom_with_idx(&mol, 1).unwrap();
    assert_eq!(group.parent_atoms(), [1]);
}

#[test]
fn bond_type_tells_crossing_from_contained_bonds() {
    let mol = sru_mol("CCOCC");
    let group = &mol.substance_groups()[0];
    assert_eq!(group.bond_type(&mol, 0).unwrap(), SGroupBondType::XBond);
    assert!(
        group.bond_type(&mol, 1).is_err(),
        "bond 1 is not in the group"
    );
}

#[test]
fn cstates_need_crossing_bonds() {
    let mut mol = sru_mol("CCOCC");
    let mut group = mol.substance_groups()[0].clone();
    group
        .add_cstate(&mol, 0, Point3D::new(1.0, 2.0, 3.0))
        .unwrap();
    assert_eq!(group.cstates()[0].bond_idx, 0);
    // a group of two atoms and all of their bonds has no crossing bond to spare
    group.set_bonds(&mol, vec![0, 1, 2, 3]).unwrap();
    assert!(
        group
            .add_cstate(&mol, 1, Point3D::new(0.0, 0.0, 0.0))
            .is_err()
    );
    mol.substance_groups_mut()[0] = group;
}

#[test]
fn removing_atoms_and_bonds_drops_or_renumbers_groups() {
    for (removed, expected) in [(0, 0), (4, 0), (5, 1)] {
        let mut mol = sru_mol("CCOCC.CCN");
        mol.remove_atom(removed).unwrap();
        assert_eq!(mol.substance_groups().len(), expected, "atom {removed}");
    }
    let mut mol = sru_mol("CCOCC");
    mol.remove_bond(0, 1).unwrap();
    assert!(mol.substance_groups().is_empty());
    let mut mol = sru_mol("CCOCC");
    mol.remove_bond(1, 2).unwrap();
    assert_eq!(summary(&mol), [(vec![1, 2], vec![0, 2])]);
}

#[test]
fn inserting_a_molecule_shifts_its_groups() {
    let mut mol = sru_mol("CCOCC.CCN");
    let other = sru_mol("CCOCC.CCN");
    mol.insert_mol(&other).unwrap();
    assert_eq!(
        summary(&mol),
        [(vec![1, 2], vec![0, 3]), (vec![9, 10], vec![6, 9])]
    );
}

#[test]
fn fragments_and_subsets_keep_the_groups_they_contain() {
    let mol = sru_mol("CCOCC.CCN");
    let (frags, _) = get_mol_frags_as_mols(&mol, true, true).unwrap();
    assert_eq!(summary(&frags[0]), [(vec![1, 2], vec![0, 3])]);
    assert!(frags[1].substance_groups().is_empty());

    let options = SubsetOptions::default();
    let (whole, _) = copy_mol_subset(&mol, &[0, 1, 2, 3, 4], &options).unwrap();
    assert_eq!(summary(&whole), [(vec![1, 2], vec![0, 3])]);
    let (part, _) = copy_mol_subset(&mol, &[1, 2, 3], &options).unwrap();
    assert!(part.substance_groups().is_empty());
}

#[test]
fn remove_hs_drops_hydrogens_from_groups() {
    let mut mol = parse("CCOCC");
    add_hs(&mut mol, &AddHsParameters::default(), None).unwrap();
    let mut group = SubstanceGroup::new("SRU");
    group.set_atoms(&mol, vec![1, 2, 3, 4, 10]).unwrap();
    group.set_bonds(&mol, vec![0]).unwrap();
    add_substance_group(&mut mol, group);
    remove_hs(&mut mol, &RemoveHsParameters::default(), true).unwrap();
    assert_eq!(summary(&mol), [(vec![1, 2, 3, 4], vec![0])]);
}

#[test]
fn group_properties_and_type_checks() {
    let mut group = SubstanceGroup::new("DAT");
    group.set_prop("FIELDNAME", "fld", false).unwrap();
    group.set_prop("ID", 3u32, false).unwrap();
    assert_eq!(group.get_prop("TYPE"), Some(&RDValue::String("DAT".into())));
    assert!(checks::is_valid_type("SRU") && !checks::is_valid_type("XXX"));
    assert!(checks::is_valid_sub_type("ALT") && checks::is_valid_connect_type("HT"));

    let mut mol = parse("CC");
    assert!(checks::is_substance_group_id_free(&mol, 3));
    add_substance_group(&mut mol, group);
    assert!(!checks::is_substance_group_id_free(&mol, 3));
}

#[test]
fn groups_survive_a_pickle_round_trip() {
    use rdoxide::graphmol::pickle::{flags, mol_from_pickle, pickle_mol};

    let mut mol = sru_mol("CCOCC");
    mol.substance_groups_mut()[0].add_bracket([
        Point3D::new(0.0, 0.0, 0.0),
        Point3D::new(1.0, 1.0, 0.0),
        Point3D::new(0.0, 0.0, 0.0),
    ]);
    mol.substance_groups_mut()[0].add_attach_point(1, -1, "");
    let bytes = pickle_mol(&mol, flags::ALL_PROPS).unwrap();
    let copy = mol_from_pickle(&bytes, flags::ALL_PROPS).unwrap();
    assert_eq!(summary(&copy), summary(&mol));
    assert_eq!(copy.substance_groups()[0].brackets().len(), 1);
    assert_eq!(copy.substance_groups()[0].attach_points().len(), 1);
    assert_eq!(pickle_mol(&copy, flags::ALL_PROPS).unwrap(), bytes);
}
