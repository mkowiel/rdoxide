//! `StructChecker/testStructChecker.cpp`.

use std::sync::Arc;

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_block, mol_from_mol_file};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::struct_checker::atom_symbol::atom_symbol_match;
use rdoxide::graphmol::struct_checker::pattern::{aa_match, check_atoms};
use rdoxide::graphmol::struct_checker::stereo::check_stereo;
use rdoxide::graphmol::struct_checker::structure_flags::*;
use rdoxide::graphmol::struct_checker::utilities::setup_neighbourhood;
use rdoxide::graphmol::struct_checker::{
    StructChecker, StructCheckerOptions, parse_options_json, string_to_augmented_atom,
};

const DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/rdkit/Code/GraphMol/StructChecker/test/"
);

fn have_data() -> bool {
    std::path::Path::new(DIR).exists()
}

fn smiles(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn molblock(s: &str) -> RWMol {
    mol_from_mol_block(s, &MolFileParserParams::default())
        .unwrap()
        .unwrap()
}

fn to_smiles(m: &mut RWMol) -> String {
    mol_to_smiles(m, &SmilesWriteParams::default()).unwrap()
}

fn load_options(options: &mut StructCheckerOptions) {
    assert!(options.load_good_augmented_atoms(&format!("{DIR}checkfgs.chk")));
    assert!(options.load_acidic_augmented_atoms(&format!("{DIR}checkfgs.aci")));
    assert!(
        options
            .load_augmented_atom_translations(&format!("{DIR}checkfgs.trn"))
            .unwrap()
    );
    assert!(options.load_tautomer_data(&format!("{DIR}tautomer.sdf")));
    assert!(options.load_tautomer_data(&format!("{DIR}tautomer.rdf")));
}

const NITRO: &str = r#"nitro.mol
  ChemDraw08311606582D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.7145   -0.6188    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    0.0000    0.6188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.2062    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    0.7145   -0.6188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0      
  2  3  2  0      
  3  4  1  0      
M  CHG  2   1  -1   3   1
M  END"#;

const SUBSTANCE_310975001: &str = r#"310975001
  -OEChem-06071611182D

 10 10  0     1  0  0  0  0  0999 V2000
    1.1317   -0.3264    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0
    0.3470    0.7535    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
    0.3470   -0.0715    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1317    1.0084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6166    0.3410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2363    1.3369    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3452   -1.1233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0332    1.1233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7618   -1.7067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6166    1.7067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  5  1  0  0  0  0
  1  7  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  6  1  0  0  0  0
  4  5  1  0  0  0  0
  6  8  1  0  0  0  0
  7  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  END
> <EXPECTED>
['stereo_error']

> <GOT>
['atom_check_failed']

$$$$
"#;

#[test]
fn flags() {
    let f = STEREO_ERROR;
    let s = StructChecker::structure_flags_to_string(f);
    assert_eq!(StructChecker::string_to_structure_flags(&s), f);
    let f = STEREO_ERROR | TRANSFORMED;
    let s = StructChecker::structure_flags_to_string(f);
    assert_eq!(s, "STEREO_ERROR,TRANSFORMED");
    assert_eq!(StructChecker::string_to_structure_flags(&s), f);
    let s = StructChecker::structure_flags_to_string(0xFFFF);
    assert_eq!(
        StructChecker::string_to_structure_flags(&s),
        0xFFFF & !0x0080
    );
    let s = " STEREO_ERROR ,\t TRANSFORMED [xXx}";
    assert_eq!(
        StructChecker::string_to_structure_flags(s),
        STEREO_ERROR | TRANSFORMED
    );
}

#[test]
fn options_json() {
    let mut options = StructCheckerOptions::default();
    assert!(parse_options_json("{}", &mut options));
    assert!(!parse_options_json("{...error..}", &mut options));
    assert!(!parse_options_json("", &mut options));
    assert!(parse_options_json(
        "{\"Verbose\": true, \"CheckStereo\": true, \"MaxMolSize\": 10, \"AcidityLimit\": 1.5}",
        &mut options
    ));
    assert!(options.verbose && options.check_stereo);
    assert_eq!(options.max_mol_size, 10);
    assert_eq!(options.acidity_limit, 1.5);
    assert!(!parse_options_json(
        "{\"Verbose\": \"maybe\"}",
        &mut options
    ));
}

#[test]
fn options_default() {
    let chk = StructChecker::new();
    for smi in [
        "COC(=O)C",
        "COC(=O)C(\\C)=C\\C1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC=CC=C",
    ] {
        let mut m = smiles(smi);
        assert_eq!(chk.check_mol_structure(&mut m).unwrap(), 0, "{smi}");
    }
}

#[test]
fn augmented_atom_strings() {
    let aa = string_to_augmented_atom("@C+1.(=O-1:2)(-N,O)").unwrap();
    assert_eq!(aa.atom_symbol, "C");
    assert_eq!(aa.charge, 1);
    assert_eq!(aa.ligands.len(), 2);
    assert_eq!(aa.ligands[0].substitution_count, 2);
    assert_eq!(aa.ligands[1].atom_symbol, "N,O");
    assert!(string_to_augmented_atom("C(-N").is_none());
    assert!(string_to_augmented_atom("1C").is_none());
}

#[test]
fn check_atom_with_default_good_atoms() {
    let chk = StructChecker::new();
    let mut m = smiles("COC(=O)C");
    assert_eq!(
        chk.check_mol_structure(&mut m).unwrap() & ATOM_CHECK_FAILED,
        0
    );
}

#[test]
fn check_match() {
    assert!(atom_symbol_match("C", "C,N,O"));
    assert!(atom_symbol_match("Na", "alk"));
    assert!(atom_symbol_match("Cl", "X"));
    assert!(!atom_symbol_match("H", "A"));
    let mol = smiles("CO");
    let nbp = setup_neighbourhood(&mol);
    // (the upstream test uses "I+1": the charge is that of the whole list)
    let aa = string_to_augmented_atom("C(-N,O,P,S,I+1)").unwrap();
    assert!(!aa_match(&mol, 0, &aa, &[], &nbp, false));
    let aa = string_to_augmented_atom("C(-N,O,P,S,I)").unwrap();
    assert!(aa_match(&mol, 0, &aa, &[], &nbp, true));
    assert!(!check_atoms(&mol, &[aa], false));
}

#[test]
fn check_atom_files() {
    if !have_data() {
        return;
    }
    let mut options = StructCheckerOptions::default();
    load_options(&mut options);
    let chk = StructChecker::with_options(options);
    let names = [
        "Substance_310925001_310950000-003900.sdf",
        "Substance_310925001_310950000-012197.sdf",
        "Substance_310925001_310950000-021440.sdf",
        "Substance_310925001_310950000-021442.sdf",
    ];
    let reference = [ATOM_CHECK_FAILED | TRANSFORMED, 0, 0, 0];
    for (name, r) in names.iter().zip(reference) {
        let mut m = mol_from_mol_file(
            &format!("{DIR}atom_check_failed/{name}"),
            &MolFileParserParams::default(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(chk.check_mol_structure(&mut m).unwrap(), r, "{name}");
    }
}

#[test]
fn stereo() {
    if !have_data() {
        return;
    }
    let mut options = StructCheckerOptions::default();
    assert!(options.load_good_augmented_atoms(&format!("{DIR}checkfgs.chk")));
    assert!(options.load_acidic_augmented_atoms(&format!("{DIR}checkfgs.aci")));
    let chk = StructChecker::with_options(options);
    for smi in [
        "COC(=O)C(\\C)=C\\C1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC=CC=C",
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2",
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1",
    ] {
        let mut m = smiles(smi);
        chk.check_mol_structure(&mut m).unwrap();
    }
    let mut m = molblock(SUBSTANCE_310975001);
    assert!(!check_stereo(&m).unwrap());
    chk.check_mol_structure(&mut m).unwrap();
}

#[test]
fn nitro() {
    if !have_data() {
        return;
    }
    let mut options = StructCheckerOptions::default();
    assert!(options.load_good_augmented_atoms(&format!("{DIR}checkfgs.chk")));
    let chk = StructChecker::with_options(options);
    let mut m = molblock(NITRO);
    assert_eq!(
        chk.check_mol_structure(&mut m).unwrap() & ATOM_CHECK_FAILED,
        0
    );
}

#[test]
fn transform_tautomer() {
    let mut options = StructCheckerOptions {
        remove_minor_fragments: true,
        check_collisions: true,
        collision_limit_percent: 3,
        check_stereo: true,
        max_mol_size: 999,
        ..Default::default()
    };
    options.set_tautomer_data(vec![Arc::new(smiles("C=C"))], vec![Arc::new(smiles("C-C"))]);
    let chk = StructChecker::with_options(options);
    let mut m = smiles("FC=C(F)O");
    let flags = chk.check_mol_structure(&mut m).unwrap();
    assert_eq!(to_smiles(&mut m), "OC(F)CF");
    assert_eq!(flags, TAUTOMER_TRANSFORMED);
}

#[test]
fn transform_atoms() {
    let mut options = StructCheckerOptions {
        collision_limit_percent: 3,
        max_mol_size: 999,
        ..Default::default()
    };
    options.set_augmented_atom_translations(vec![
        (
            string_to_augmented_atom("C(=O)(-C)").unwrap(),
            string_to_augmented_atom("C(-O)(-C)").unwrap(),
        ),
        (
            string_to_augmented_atom("C(=C)(-F)").unwrap(),
            string_to_augmented_atom("C(-C)(=F)").unwrap(),
        ),
    ]);
    let chk = StructChecker::with_options(options);
    for smi in ["C(=O)C", "FC=CO"] {
        let mut m = smiles(smi);
        let flags = chk.check_mol_structure(&mut m).unwrap();
        assert_ne!(flags & TRANSFORMED, 0, "{smi}");
        assert_eq!(flags & FRAGMENTS_FOUND, 0);
    }
}

#[test]
fn augmented_atom_translations_to_atom_list_query() {
    let mut options = StructCheckerOptions {
        collision_limit_percent: 3,
        max_mol_size: 999,
        ..Default::default()
    };
    options.set_augmented_atom_translations(vec![(
        string_to_augmented_atom("C(=O)(-C)").unwrap(),
        string_to_augmented_atom("C(-O,N,S,Ca)(-C)").unwrap(),
    )]);
    let chk = StructChecker::with_options(options);
    let mut m = smiles("C(=O)C");
    let flags = chk.check_mol_structure(&mut m).unwrap();
    assert_ne!(flags & TRANSFORMED, 0);
    // the former oxygen is now a query atom for O, N, S or Ca
    assert!(
        m.atoms()
            .iter()
            .any(|a| a.has_query() && a.atomic_num() == 8)
    );
}

#[test]
fn carboxylic_acids() {
    if !have_data() {
        return;
    }
    let mut options = StructCheckerOptions::default();
    load_options(&mut options);
    options.remove_minor_fragments = true;
    let chk = StructChecker::with_options(options);
    let mut m = mol_from_mol_file(
        &format!("{DIR}Substance_310925001_310950000-000075-000075_Transf_Acidis.sdf"),
        &MolFileParserParams::default(),
    )
    .unwrap()
    .unwrap();
    let flags = chk.check_mol_structure(&mut m).unwrap();
    assert_eq!(flags & ATOM_CHECK_FAILED, 0);
}

#[test]
fn size_checks() {
    let options = StructCheckerOptions {
        max_mol_size: 2,
        ..Default::default()
    };
    let chk = StructChecker::with_options(options);
    let mut m = smiles("CCCC");
    assert_eq!(chk.check_mol_structure(&mut m).unwrap(), SIZE_CHECK_FAILED);
    let mut empty = RWMol::new();
    assert_eq!(
        StructChecker::new()
            .check_mol_structure(&mut empty)
            .unwrap(),
        SIZE_CHECK_FAILED
    );
}
