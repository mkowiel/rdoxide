//! BioCAD TPL reading and writing.

use std::path::Path;

use rdoxide::geometry::Point3D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::conformer::Conformer;
use rdoxide::graphmol::file_parsers::tpl::{
    TplParserParams, mol_from_tpl_block, mol_from_tpl_file, mol_to_tpl_text,
};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::rdgeneral::RDProps;

const DATA: &str = "../../upstream/rdkit/Code/GraphMol/FileParsers/test_data/cmpd2.tpl";

fn ethanol_with_two_confs() -> RWMol {
    let mut mol = mol_from_smiles("CCO", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    mol.set_prop("_Name", "ethanol".to_string(), false).unwrap();
    for shift in [0.0, 0.5] {
        let mut conf = Conformer::new(3);
        for i in 0..3 {
            conf.set_atom_pos(i, Point3D::new(1.5 * f64::from(i) + shift, 0.25, -0.5))
                .unwrap();
        }
        mol.add_conformer(conf, true).unwrap();
    }
    mol
}

#[test]
fn written_tpl_reads_back() {
    let mol = ethanol_with_two_confs();
    let text = mol_to_tpl_text(&mol, "_GasteigerCharge", false);
    assert!(text.starts_with(
        "BioCAD format, all rights reserved\nOutput from RDKit\nNAME ethanol\nPROP 7 1\n3 2\n"
    ));
    assert!(text.contains("\nCONFS 1\nNAME conformer_1\n"));
    let back = mol_from_tpl_block(&text, &TplParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(back.num_atoms(), 3);
    assert_eq!(back.num_bonds(), 2);
    assert_eq!(back.num_conformers(), 2);
    let p = back.conformer(0).unwrap().positions()[1];
    assert!((p.x - 1.5).abs() < 1e-4 && (p.y - 0.25).abs() < 1e-4 && (p.z + 0.5).abs() < 1e-4);
}

#[test]
fn first_conformer_can_be_written_twice_and_skipped() {
    let mol = ethanol_with_two_confs();
    let text = mol_to_tpl_text(&mol, "_GasteigerCharge", true);
    assert!(text.contains("CONFS 1\nNAME conformer_0\n"));
    let skip = TplParserParams {
        skip_first_conf: true,
        ..Default::default()
    };
    let back = mol_from_tpl_block(&text, &skip).unwrap().unwrap();
    assert_eq!(back.num_conformers(), 1);
}

#[test]
fn molecules_without_conformers_are_not_written() {
    let mol = mol_from_smiles("CC", &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(mol_to_tpl_text(&mol, "_GasteigerCharge", false), "");
}

#[test]
fn truncated_input_gives_none_or_an_error() {
    assert!(
        mol_from_tpl_block("", &TplParserParams::default())
            .unwrap()
            .is_none()
    );
    let bad = "BioCAD\ncomment\nNAME x\nPROP 7 1\n1 0\n";
    assert!(mol_from_tpl_block(bad, &TplParserParams::default()).is_err());
}

#[test]
fn upstream_file_is_read() {
    if !Path::new(DATA).exists() {
        return;
    }
    let mol = mol_from_tpl_file(DATA, &TplParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(mol.num_atoms(), 12);
    assert_eq!(mol.num_conformers(), 2);
}
