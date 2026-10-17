//! Canonical enhanced stereo groups.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::smiles_parse::canonicalize_stereo_groups::{
    StereoGroupAbsOptions, canonicalize_stereo_groups,
};
use rdoxide::graphmol::smiles_parse::cx_write::{
    RestoreBondDirOption, cx_fields, mol_to_cx_smiles,
};
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn parse(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn cx(mol: &RWMol) -> String {
    mol_to_cx_smiles(
        mol,
        &SmilesWriteParams::default(),
        cx_fields::ALL,
        RestoreBondDirOption::Clear,
    )
    .unwrap()
}

#[test]
fn canonicalization_is_idempotent() {
    let mol = parse("C[C@H](O)C[C@@H](F)Cl |o1:1,o2:4|");
    let once = canonicalize_stereo_groups(
        &mol,
        StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist,
        12,
    )
    .unwrap();
    let twice = canonicalize_stereo_groups(
        &once,
        StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist,
        12,
    )
    .unwrap();
    assert_eq!(cx(&once), cx(&twice));
}

#[test]
fn groups_are_never_split_further() {
    // the groups can only stay as they are or merge
    let mol = parse("C[C@H](O)C[C@@H](O)C |o1:1,o2:4|");
    let res = canonicalize_stereo_groups(
        &mol,
        StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist,
        12,
    )
    .unwrap();
    assert!(res.stereo_groups().len() <= 2);
}

#[test]
fn absolute_group_options() {
    let mol = parse("C[C@H](O)Cl");
    let always =
        canonicalize_stereo_groups(&mol, StereoGroupAbsOptions::AlwaysInclude, 12).unwrap();
    assert_eq!(always.stereo_groups().len(), 1);
    let never = canonicalize_stereo_groups(&mol, StereoGroupAbsOptions::NeverInclude, 12).unwrap();
    assert!(never.stereo_groups().is_empty());
}

#[test]
fn too_many_groups_leave_the_molecule_alone() {
    let mol = parse("C[C@H](O)C[C@@H](F)Cl |o1:1,o2:4|");
    let res = canonicalize_stereo_groups(
        &mol,
        StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist,
        1,
    )
    .unwrap();
    assert_eq!(res.stereo_groups().len(), 2);
    assert_eq!(cx(&res), cx(&mol));
}
