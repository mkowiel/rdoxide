//! CXSMILES reading; expectations from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::bond::{BondDir, BondType};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::stereo_group::StereoGroupType;
use rdoxide::rdgeneral::RDProps;
use rdoxide::rdgeneral::common_properties as cp;

fn parse_with(smiles: &str, sanitize: bool, strict: bool) -> Option<RWMol> {
    let params = SmilesParserParams {
        sanitize,
        strict_cx_smiles: strict,
        ..Default::default()
    };
    mol_from_smiles(smiles, &params).ok().flatten()
}

fn parse(smiles: &str) -> RWMol {
    parse_with(smiles, true, true).expect("parses")
}

#[test]
fn coordinates_make_a_conformer() {
    let mol = parse("CCO |(0,0;1.5,0;3,0.5)|");
    assert_eq!(mol.num_conformers(), 1);
    let conf = &mol.conformers()[0];
    assert!(!conf.is_3d());
    assert_eq!(conf.positions()[2].y, 0.5);
    assert!(parse("CC |(0,0,0;1,1,1)|").conformers()[0].is_3d());
    assert!(!parse("CC |(0,0,0;1,1,0)|").conformers()[0].is_3d());
}

#[test]
fn labels_and_properties_are_stored_on_atoms() {
    let mol = parse("C* |$;R1$|");
    let label = mol.atom(1).get_prop(cp::atomLabel).unwrap().to_rd_string();
    assert_eq!(label, "R1");
    let mol = parse("CC |atomProp:0.foo.bar|");
    assert_eq!(mol.atom(0).get_prop("foo").unwrap().to_rd_string(), "bar");
    let mol = parse("CC |$_AV:;a$|");
    assert!(!mol.atom(0).has_prop(cp::molFileValue));
    assert!(mol.atom(1).has_prop(cp::molFileValue));
}

#[test]
fn radicals_follow_the_electron_count_codes() {
    let counts = |smiles: &str| -> Vec<u32> {
        let mol = parse_with(smiles, false, true).unwrap();
        mol.atoms()
            .iter()
            .map(|a| a.num_radical_electrons())
            .collect()
    };
    assert_eq!(counts("C.C.C |^1:0,2|"), [1, 0, 1]);
    assert_eq!(counts("CC |^4:1|"), [0, 2]);
    assert_eq!(counts("CC |^6:0|"), [3, 0]);
    assert!(parse_with("CC |^8:0|", false, true).is_none());
}

#[test]
fn enhanced_stereo_groups_merge_by_type_and_id() {
    let mol = parse("CC[C@H](F)[C@H](Cl)Br |o1:2,o1:4|");
    let groups = mol.stereo_groups();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].group_type(), StereoGroupType::Or);
    assert_eq!(groups[0].atoms(), [2, 4]);
    assert!(parse_with("C[C@H](F)Cl |o:1,1|", true, true).is_none());
}

#[test]
fn bond_extensions() {
    let mol = parse("[Fe]N |C:1.0|");
    assert_eq!(mol.bond(0).bond_type(), BondType::Dative);
    assert_eq!(mol.bond(0).begin_atom_idx(), 1);
    let mol = parse("C1CCCCC~CCCC1 |Z:5|");
    assert_eq!(mol.bond(5).bond_type(), BondType::Zero);
    let mol = parse_with("CC(O)Cl |w:1.0|", false, true).unwrap();
    assert_eq!(mol.bond(0).bond_dir(), BondDir::Unknown);
}

#[test]
fn names_and_failures() {
    let name = |s: &str| {
        parse(s)
            .get_prop(cp::_Name)
            .map(|v| v.to_rd_string())
            .unwrap_or_default()
    };
    assert_eq!(name("CCO ethanol"), "ethanol");
    assert_eq!(name("CCO |$_AV:;;x$| ethanol"), "ethanol");
    // a lenient parse keeps the molecule, and the name is not read
    let mol = parse_with("CC |o:0,0| name", true, false).unwrap();
    assert!(mol.has_prop("_CXSMILES_Data"));
    assert!(!mol.has_prop(cp::_Name));
    assert!(parse_with("CC |o:0,0| name", true, true).is_none());
}

#[test]
fn polymer_and_data_substance_groups() {
    let mol = parse("*-CCO-* |$star_e;;;;star_e$,Sg:n:1,2,3::eu|");
    let sgs = mol.substance_groups();
    assert_eq!(sgs.len(), 1);
    assert_eq!(sgs[0].get_prop("TYPE").unwrap().to_rd_string(), "SRU");
    assert_eq!(sgs[0].atoms(), [1, 2, 3]);
    // the bonds out of the group were found from its end atoms
    assert_eq!(sgs[0].bonds().len(), 2);
    assert_eq!(sgs[0].get_prop("CONNECT").unwrap().to_rd_string(), "EU");

    let mol = parse("*CC |$star_e;;$,SgD:0:querydata:val::::|");
    let sg = &mol.substance_groups()[0];
    assert_eq!(sg.get_prop("TYPE").unwrap().to_rd_string(), "DAT");
    assert_eq!(
        sg.get_prop("FIELDNAME").unwrap().to_rd_string(),
        "querydata"
    );
    assert_eq!(sg.get_prop("index").unwrap().as_uint(), Some(1));

    // a parent/child pair
    let mol = parse("*CC* |$star_e;;;star_e$,Sg:n:1::ht,Sg:n:2::ht,SgH:0:1|");
    assert_eq!(mol.substance_groups().len(), 2);
    assert!(mol.substance_groups()[1].has_prop("PARENT"));
    // groups with unknown types do not parse
    assert!(parse_with("CC |Sg:zzz:0|", true, true).is_none());
}

#[test]
fn link_nodes_become_a_property() {
    let link_nodes = |s: &str| {
        parse(s)
            .get_prop(cp::molFileLinkNodes)
            .map(|v| v.to_rd_string())
            .unwrap_or_default()
    };
    // the outer atoms are implied for an atom with two neighbors
    assert_eq!(link_nodes("OC1CCCC1 |LN:4:1.3|"), "1 3 2 5 4 5 6");
    assert_eq!(
        link_nodes("FC1CCC(O)C1 |LN:1:1.3.2.6,4:1.4.3.6|"),
        "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
    );
    // more than two neighbors need the outer atoms
    assert!(parse_with("FC1CCC(O)C1 |LN:1:1.3|", true, true).is_none());
}

#[test]
fn variable_attachment_points() {
    let mol = parse("CO*.C1=CC=NC=C1 |m:2:3.5.4|");
    let bond = mol.bond(mol.bond_idx_between(1, 2).unwrap());
    assert_eq!(
        bond.get_prop(cp::_MolFileBondEndPts)
            .unwrap()
            .to_rd_string(),
        "(3 4 6 5)"
    );
    assert_eq!(
        bond.get_prop(cp::_MolFileBondAttach)
            .unwrap()
            .to_rd_string(),
        "ANY"
    );
    // the attached atom must have a single bond
    assert!(parse_with("CC(C)C |m:1:3|", true, true).is_none());
}
