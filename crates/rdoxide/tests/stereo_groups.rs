//! Enhanced stereo helpers; expectations from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::chirality::legacy::assign_stereochemistry;
use rdoxide::graphmol::chirality::stereo_groups::{
    StereoAnnotationLabels, add_stereo_annotations, find_meso_centers, simplify_enhanced_stereo,
    translate_chiral_flag_to_stereo_groups,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use rdoxide::rdgeneral::RDProps;
use rdoxide::rdgeneral::common_properties as cp;

use StereoGroupType::{Absolute, And, Or};

fn grouped(smiles: &str, groups: &[(StereoGroupType, &[u32], u32)]) -> RWMol {
    let mut mol = mol_from_smiles(smiles, &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assign_stereochemistry(&mut mol, true, true, false).unwrap();
    let groups = groups
        .iter()
        .map(|&(kind, atoms, id)| StereoGroup::new(kind, atoms.to_vec(), vec![], id).unwrap())
        .collect();
    mol.set_stereo_groups(groups).unwrap();
    mol
}

fn atom_note(mol: &RWMol, idx: u32) -> Option<String> {
    mol.atom(idx)
        .get_prop(cp::atomNote)
        .map(|v| v.to_rd_string())
}

fn group_summary(mol: &RWMol) -> Vec<(StereoGroupType, Vec<u32>)> {
    mol.stereo_groups()
        .iter()
        .map(|g| (g.group_type(), g.atoms().to_vec()))
        .collect()
}

#[test]
fn annotations_combine_groups_and_cip_codes() {
    let mut mol = grouped(
        "C[C@H](F)[C@@H](Cl)Br",
        &[(Absolute, &[1], 0), (Or, &[3], 1)],
    );
    add_stereo_annotations(&mut mol, &StereoAnnotationLabels::default());
    assert_eq!(atom_note(&mol, 1).as_deref(), Some("abs (S)"));
    assert_eq!(atom_note(&mol, 3).as_deref(), Some("or1"));

    let mut custom = grouped(
        "C[C@H](F)[C@@H](Cl)Br",
        &[(Absolute, &[1], 0), (Or, &[3], 1)],
    );
    let labels = StereoAnnotationLabels {
        abs: "A{id}{cip}".into(),
        or: "O{id}".into(),
        and: String::new(),
        cip: String::new(),
        bond: "B{cip}".into(),
    };
    add_stereo_annotations(&mut custom, &labels);
    assert_eq!(atom_note(&custom, 1).as_deref(), Some("A0S"));
    assert_eq!(atom_note(&custom, 3).as_deref(), Some("O1"));
}

#[test]
fn ungrouped_centers_and_double_bonds_get_cip_labels() {
    let mut mol = grouped("F/C=C/[C@H](Cl)C", &[]);
    add_stereo_annotations(&mut mol, &StereoAnnotationLabels::default());
    assert_eq!(atom_note(&mol, 3).as_deref(), Some("(R)"));
    let note = mol.bond(1).get_prop(cp::bondNote).map(|v| v.to_rd_string());
    assert_eq!(note.as_deref(), Some("(E)"));
}

#[test]
fn simplifies_a_single_group_covering_all_centers() {
    let mut mol = grouped("C[C@H](F)[C@@H](Cl)Br", &[(And, &[1, 3], 1)]);
    simplify_enhanced_stereo(&mut mol, true);
    assert_eq!(
        mol.get_prop(cp::molNote)
            .map(|v| v.to_rd_string())
            .as_deref(),
        Some("AND enantiomer")
    );
    assert!(mol.stereo_groups().is_empty());
    assert!(!mol.atom(1).has_prop(cp::_CIPCode));

    let mut kept = grouped("C[C@H](F)Cl", &[(Or, &[1], 1)]);
    simplify_enhanced_stereo(&mut kept, false);
    assert_eq!(kept.stereo_groups().len(), 1);
    assert!(kept.has_prop(cp::molNote));

    let mut untouched = grouped("C[C@H](F)Cl", &[(Absolute, &[1], 0)]);
    simplify_enhanced_stereo(&mut untouched, true);
    assert!(!untouched.has_prop(cp::molNote));
    assert_eq!(untouched.stereo_groups().len(), 1);
}

#[test]
fn finds_meso_centers() {
    let parse = |s: &str| {
        mol_from_smiles(s, &SmilesParserParams::default())
            .unwrap()
            .unwrap()
    };
    let mut meso = parse("C[C@H](O)[C@@H](C)O");
    assert_eq!(find_meso_centers(&mut meso, true, false).unwrap(), [(1, 3)]);
    assert!(meso.atom(1).has_prop(cp::_mesoOtherAtom));
    let mut chiral = parse("C[C@H](O)[C@H](C)O");
    assert!(
        find_meso_centers(&mut chiral, true, false)
            .unwrap()
            .is_empty()
    );
    let mut achiral = parse("CC");
    assert!(
        find_meso_centers(&mut achiral, true, false)
            .unwrap()
            .is_empty()
    );
    let mut ring = parse("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O");
    let pairs = find_meso_centers(&mut ring, true, false).unwrap();
    assert_eq!((pairs[0], pairs.len()), ((1, 2), 10));
}

#[test]
fn chiral_flag_becomes_a_stereo_group() {
    let flagged = |flag: i32, groups: &[(StereoGroupType, &[u32], u32)]| {
        let mut mol = grouped("C[C@H](F)[C@@H](Cl)Br", groups);
        mol.set_prop(cp::_MolFileChiralFlag, flag, false).unwrap();
        mol
    };
    let mut mol = flagged(1, &[]);
    translate_chiral_flag_to_stereo_groups(&mut mol, And).unwrap();
    assert_eq!(group_summary(&mol), [(Absolute, vec![1, 3])]);
    assert!(!mol.has_prop(cp::_MolFileChiralFlag));

    let mut zero = flagged(0, &[]);
    translate_chiral_flag_to_stereo_groups(&mut zero, Or).unwrap();
    assert_eq!(group_summary(&zero), [(Or, vec![1, 3])]);

    // existing groups keep their centers; an ABS group takes the rest
    let mut partial = flagged(1, &[(Or, &[1], 1)]);
    translate_chiral_flag_to_stereo_groups(&mut partial, And).unwrap();
    assert_eq!(
        group_summary(&partial),
        [(Or, vec![1]), (Absolute, vec![3])]
    );
    let mut merged = flagged(0, &[(Absolute, &[1], 0)]);
    translate_chiral_flag_to_stereo_groups(&mut merged, And).unwrap();
    assert_eq!(
        group_summary(&merged),
        [(Absolute, vec![1]), (And, vec![3])]
    );

    let mut no_flag = grouped("C[C@H](F)Cl", &[]);
    translate_chiral_flag_to_stereo_groups(&mut no_flag, And).unwrap();
    assert!(no_flag.stereo_groups().is_empty());
}
