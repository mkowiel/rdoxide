//! Expanding and collapsing explicit attachment points; the expectations are
//! from RDKit 2026.09.1 (`v2_file_parsers_catch.cpp`, `catch_graphmol.cpp`).

use rdoxide::graphmol::attachment_points::{
    collapse_attachment_points, expand_attachment_points, is_marked_attachment_point,
};
use rdoxide::graphmol::file_parsers::mol_file_parser::{MolFileParserParams, mol_from_mol_block};
use rdoxide::graphmol::smiles_parse::cx_write::{
    RestoreBondDirOption, cx_fields, mol_to_cx_smiles,
};
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{ROMol, RWMol};
use rdoxide::rdgeneral::common_properties as cp;
use rdoxide::rdgeneral::{RDProps, RDValue};

const MOL_BLOCK: &str = r#"
  Mrv2211 02082417032D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.7083 3.25 0 0
M  V30 2 O -3.3747 4.02 0 0 ATTCHPT=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"#;

fn cx(mol: &ROMol, flags: u32) -> String {
    mol_to_cx_smiles(
        mol,
        &SmilesWriteParams::default(),
        flags,
        RestoreBondDirOption::Clear,
    )
    .unwrap()
}

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn attach_value(mol: &ROMol, idx: u32) -> Option<i32> {
    mol.atom(idx)
        .get_prop(cp::molAttachPoint)
        .and_then(RDValue::as_int)
}

#[test]
fn molfile_attachment_points_are_expanded_on_request() {
    let plain = mol_from_mol_block(MOL_BLOCK, &MolFileParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(plain.num_atoms(), 2);
    let params = MolFileParserParams {
        expand_attachment_points: true,
        ..Default::default()
    };
    let expanded = mol_from_mol_block(MOL_BLOCK, &params).unwrap().unwrap();
    assert_eq!(expanded.num_atoms(), 3);
    assert_eq!(cx(&expanded, cx_fields::ALL_BUT_COORDS), "*OC |$_AP1;;$|");
    // the new atom has coordinates too
    assert_eq!(
        cx(&expanded, cx_fields::ALL),
        "*OC |(-2.50869,4.52002,;-3.3747,4.02,;-4.7083,3.25,),$_AP1;;$|"
    );
}

#[test]
fn expand_handles_both_positions_and_bad_values() {
    let mut m = smiles("CC");
    m.atom_mut(0)
        .set_prop(cp::molAttachPoint, -1, false)
        .unwrap();
    m.atom_mut(1)
        .set_prop(cp::molAttachPoint, 7, false)
        .unwrap();
    expand_attachment_points(&mut m, false, true).unwrap();
    // two dummies for atom 0; the invalid value on atom 1 is left alone
    assert_eq!(m.num_atoms(), 4);
    assert!(attach_value(&m, 0).is_none());
    assert_eq!(attach_value(&m, 1), Some(7));
    assert!(m.atom(2).query().is_none());
    assert!(is_marked_attachment_point(&m, 2));
    assert!(is_marked_attachment_point(&m, 3));
}

#[test]
fn collapse_turns_dummies_back_into_annotations() {
    let mut m = smiles("*CC");
    m.atom_mut(0)
        .set_prop(cp::_fromAttachPoint, 1, false)
        .unwrap();
    collapse_attachment_points(&mut m, true).unwrap();
    assert_eq!(m.num_atoms(), 2);
    assert_eq!(attach_value(&m, 0), Some(1));

    // two on one atom become -1
    let mut m = smiles("*C*");
    m.atom_mut(0)
        .set_prop(cp::_fromAttachPoint, 1, false)
        .unwrap();
    m.atom_mut(2)
        .set_prop(cp::_fromAttachPoint, 2, false)
        .unwrap();
    collapse_attachment_points(&mut m, true).unwrap();
    assert_eq!(m.num_atoms(), 1);
    assert_eq!(attach_value(&m, 0), Some(-1));

    // unmarked dummies are only collapsed when asked for
    let mut m = smiles("*CC");
    collapse_attachment_points(&mut m, true).unwrap();
    assert_eq!(m.num_atoms(), 3);
    collapse_attachment_points(&mut m, false).unwrap();
    assert_eq!(m.num_atoms(), 2);
    assert_eq!(attach_value(&m, 0), Some(1));
}
