//! Expanding SCSR mol files; the expectations are those of RDKit 2026.09.1's
//! `macromols_catch.cpp` and `testSCSR.py`.

use rdoxide::graphmol::atom::ChiralType;
use rdoxide::graphmol::chirality::remove_non_explicit_3d_chirality;
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::file_parsers::scsr::{
    MolFromScsrParams, ScsrBaseHbondOptions, ScsrTemplateNames, mol_from_scsr_block,
    mol_from_scsr_file,
};
use rdoxide::graphmol::substance_group::SubstanceGroup;
use rdoxide::rdgeneral::{RDProps, RDValue};

const DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/rdkit/Code/GraphMol/FileParsers/test_data/macromols/"
);

use ScsrBaseHbondOptions::{Auto, UseSapAll};

/// file, expands, hydrogen bond option, atoms, bonds, groups, and the same
/// without leaving groups, then the chirality of some atoms
type Case = (
    &'static str,
    bool,
    ScsrBaseHbondOptions,
    [usize; 6],
    &'static [(u32, ChiralType)],
);

use ChiralType::{TetrahedralCCW as Ccw, TetrahedralCW as Cw};

const CASES: &[Case] = &[
    (
        "DNASlurpErrorImport.mol",
        true,
        Auto,
        [81, 90, 13, 79, 88, 11],
        &[(0, Cw), (22, Cw), (41, Cw), (61, Cw)],
    ),
    (
        "DNASlurpErrorSketch.mol",
        true,
        Auto,
        [84, 93, 14, 82, 91, 12],
        &[(51, Ccw), (60, Ccw), (68, Ccw), (76, Ccw)],
    ),
    (
        "ValenceErrorScsr.mol",
        true,
        Auto,
        [38, 39, 6, 35, 36, 3],
        &[],
    ),
    (
        "ValenceErrorScsr2.mol",
        true,
        Auto,
        [28, 28, 6, 25, 25, 3],
        &[],
    ),
    (
        "RiboseFullname.mol",
        true,
        Auto,
        [45, 49, 8, 43, 47, 6],
        &[],
    ),
    ("Conjugate.mol", true, Auto, [91, 91, 14, 87, 87, 10], &[]),
    (
        "ModifiedPeptide2.mol",
        true,
        Auto,
        [438, 444, 81, 407, 413, 50],
        &[],
    ),
    (
        "DnaBadPairs_NoCh.mol",
        true,
        Auto,
        [84, 94, 14, 80, 90, 10],
        &[],
    ),
    (
        "DnaBadPairs.mol",
        true,
        UseSapAll,
        [84, 94, 14, 80, 90, 10],
        &[],
    ),
    (
        "DnaTest.mol",
        false,
        UseSapAll,
        [254, 300, 14, 250, 296, 10],
        &[],
    ),
    (
        "wobblePairs2.mol",
        true,
        UseSapAll,
        [169, 196, 26, 165, 192, 22],
        &[],
    ),
    (
        "wobblePairs.mol",
        true,
        UseSapAll,
        [169, 196, 26, 165, 192, 22],
        &[],
    ),
    (
        "KellanRNA.mol",
        true,
        UseSapAll,
        [299, 353, 44, 295, 349, 40],
        &[],
    ),
    (
        "DnaTest2.mol",
        true,
        UseSapAll,
        [83, 97, 14, 79, 93, 10],
        &[],
    ),
    (
        "DnaTest3.mol",
        true,
        UseSapAll,
        [165, 194, 26, 161, 190, 22],
        &[],
    ),
    (
        "KellanError.mol",
        true,
        UseSapAll,
        [244, 263, 39, 236, 255, 31],
        &[],
    ),
    (
        "TestRNA2_fixed.mol",
        true,
        UseSapAll,
        [106, 117, 17, 104, 115, 15],
        &[],
    ),
    ("DnaTest.mol", true, Auto, [254, 300, 38, 250, 296, 34], &[]),
    (
        "wobblePairs2.mol",
        true,
        Auto,
        [169, 196, 26, 165, 192, 22],
        &[],
    ),
    (
        "wobblePairs.mol",
        true,
        Auto,
        [169, 196, 26, 165, 192, 22],
        &[],
    ),
    ("DnaBadPairs.mol", true, Auto, [84, 94, 14, 80, 90, 10], &[]),
    ("DnaTest2.mol", true, Auto, [83, 97, 14, 79, 93, 10], &[]),
    (
        "DnaTest3.mol",
        true,
        Auto,
        [165, 194, 26, 161, 190, 22],
        &[],
    ),
    (
        "KellanError.mol",
        true,
        Auto,
        [244, 263, 39, 236, 255, 31],
        &[],
    ),
    (
        "TestRNA2_fixed.mol",
        true,
        Auto,
        [106, 117, 17, 104, 115, 15],
        &[],
    ),
    (
        "TrastuzumabMaxPlus3Register.mol",
        true,
        UseSapAll,
        [7606, 7793, 1451, 7080, 7267, 925],
        &[],
    ),
    (
        "TrastuzumabMaxRegister.mol",
        true,
        UseSapAll,
        [7576, 7763, 1445, 7053, 7240, 922],
        &[],
    ),
    ("Mixed.mol", true, UseSapAll, [51, 54, 8, 51, 54, 8], &[]),
    (
        "CrossLink.mol",
        true,
        UseSapAll,
        [47, 48, 10, 45, 46, 8],
        &[],
    ),
    ("cyclic.mol", true, UseSapAll, [45, 47, 8, 45, 47, 8], &[]),
    ("Triplet.mol", true, UseSapAll, [30, 30, 6, 27, 27, 3], &[]),
    (
        "FromBioviaDoc.mol",
        true,
        UseSapAll,
        [27, 26, 7, 25, 24, 5],
        &[],
    ),
    (
        "testSCSR.mol",
        true,
        UseSapAll,
        [64, 66, 15, 57, 59, 8],
        &[],
    ),
    ("badAtomName.mol", false, UseSapAll, [0, 0, 0, 0, 0, 0], &[]),
    ("badClass.mol", false, UseSapAll, [0, 0, 0, 0, 0, 0], &[]),
    (
        "badClassTemplate.mol",
        false,
        UseSapAll,
        [0, 0, 0, 0, 0, 0],
        &[],
    ),
    (
        "badMissingTemplate.mol",
        false,
        UseSapAll,
        [0, 0, 0, 0, 0, 0],
        &[],
    ),
    (
        "obj3dTest.mol",
        true,
        UseSapAll,
        [27, 26, 7, 25, 24, 5],
        &[],
    ),
    (
        "obj3dTest2.mol",
        true,
        UseSapAll,
        [27, 26, 7, 25, 24, 5],
        &[],
    ),
    (
        "obj3dFoundTwice.mol",
        false,
        UseSapAll,
        [27, 26, 0, 27, 26, 0],
        &[],
    ),
    (
        "SgroupFoundTwice.mol",
        false,
        UseSapAll,
        [0, 0, 0, 0, 0, 0],
        &[],
    ),
];

fn params(sanitize: bool, strict: bool) -> MolFileParserParams {
    MolFileParserParams {
        sanitize,
        remove_hs: false,
        strict_parsing: strict,
        ..Default::default()
    }
}

#[test]
fn macromolecule_files() {
    for &(file, expands, hbonds, counts, chirality) in CASES {
        let path = format!("{DATA}{file}");
        let mut scsr = MolFromScsrParams {
            include_leaving_groups: true,
            scsr_base_hbond_options: hbonds,
            ..Default::default()
        };
        let pp = params(true, true);
        if !expands {
            assert!(mol_from_scsr_file(&path, &pp, &scsr).is_err(), "{file}");
            continue;
        }
        let mut mol = mol_from_scsr_file(&path, &pp, &scsr)
            .unwrap_or_else(|e| panic!("{file}: {e}"))
            .expect(file);
        remove_non_explicit_3d_chirality(&mut mol);
        assert_eq!(
            [
                mol.num_atoms() as usize,
                mol.num_bonds() as usize,
                mol.substance_groups().len()
            ],
            counts[..3],
            "{file}"
        );
        for &(idx, tag) in chirality {
            assert_eq!(mol.atom(idx).chiral_tag(), tag, "{file} atom {idx}");
        }
        // the same without the leaving groups
        scsr.include_leaving_groups = false;
        let query = mol_from_scsr_file(&path, &pp, &scsr)
            .unwrap_or_else(|e| panic!("{file}: {e}"))
            .expect(file);
        assert_eq!(
            [
                query.num_atoms() as usize,
                query.num_bonds() as usize,
                query.substance_groups().len()
            ],
            counts[3..],
            "{file} without leaving groups"
        );
    }
}

fn label(sg: &SubstanceGroup) -> String {
    match sg.get_prop("LABEL") {
        Some(RDValue::String(s)) => s.clone(),
        other => panic!("no label: {other:?}"),
    }
}

#[test]
fn names_with_nested_parentheses() {
    let scsr = MolFromScsrParams {
        scsr_base_hbond_options: Auto,
        ..Default::default()
    };
    let mol = mol_from_scsr_file(&format!("{DATA}Ugly.mol"), &params(true, true), &scsr)
        .unwrap()
        .unwrap();
    assert_eq!(
        [
            mol.num_atoms(),
            mol.num_bonds(),
            mol.substance_groups().len() as u32
        ],
        [140, 151, 11]
    );
    let sgroups = mol.substance_groups();
    assert!(label(&sgroups[1]).starts_with("((cPr)O(2S-Me)Et)NGly"));
    assert!(label(&sgroups[5]).starts_with("Phe(b-Me2)"));
    assert!(label(&sgroups[7]).starts_with("(PhO(2S-Me)Et)NGly"));
}

#[test]
fn three_letter_codes() {
    let path = format!("{DATA}PepTla.mol");
    let mut scsr = MolFromScsrParams {
        scsr_template_names: ScsrTemplateNames::AsEntered,
        ..Default::default()
    };
    let pp = params(false, false);
    let as_entered = mol_from_scsr_file(&path, &pp, &scsr).unwrap().unwrap();
    assert_eq!(
        [
            as_entered.num_atoms(),
            as_entered.num_bonds(),
            as_entered.substance_groups().len() as u32
        ],
        [26, 25, 7]
    );
    assert_eq!(label(&as_entered.substance_groups()[0]), "Ala_4");
    scsr.scsr_template_names = ScsrTemplateNames::UseFirstName;
    let first = mol_from_scsr_file(&path, &pp, &scsr).unwrap().unwrap();
    assert_eq!(label(&first.substance_groups()[0]), "Ala_4");
    scsr.scsr_template_names = ScsrTemplateNames::UseSecondName;
    let second = mol_from_scsr_file(&path, &pp, &scsr).unwrap().unwrap();
    assert_eq!(second.num_atoms(), 26);
    assert_eq!(label(&second.substance_groups()[0]), "A_4");
}

#[test]
fn triplet_from_a_block_and_a_file() {
    let path = format!("{DATA}Triplet.mol");
    let block = std::fs::read_to_string(&path).unwrap();
    let scsr = MolFromScsrParams {
        scsr_base_hbond_options: ScsrBaseHbondOptions::Ignore,
        ..Default::default()
    };
    let pp = params(false, false);
    for mol in [
        mol_from_scsr_block(&block, &pp, &scsr).unwrap().unwrap(),
        mol_from_scsr_file(&path, &pp, &scsr).unwrap().unwrap(),
    ] {
        assert_eq!(mol.num_atoms(), 30);
        assert_eq!(mol.substance_groups().len(), 6);
    }
    let default = mol_from_scsr_block(
        &block,
        &MolFileParserParams::default(),
        &MolFromScsrParams::default(),
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        (default.num_atoms(), default.substance_groups().len()),
        (30, 6)
    );
}

#[test]
fn double_stranded_rna() {
    let block = std::fs::read_to_string(format!("{DATA}DnaTest.mol")).unwrap();
    let pp = params(false, false);
    let mut scsr = MolFromScsrParams {
        scsr_base_hbond_options: Auto,
        ..Default::default()
    };
    let with_pairs = mol_from_scsr_block(&block, &pp, &scsr).unwrap().unwrap();
    assert_eq!(
        [
            with_pairs.num_atoms(),
            with_pairs.num_bonds(),
            with_pairs.substance_groups().len() as u32
        ],
        [254, 300, 38]
    );
    scsr.scsr_base_hbond_options = ScsrBaseHbondOptions::Ignore;
    let without = mol_from_scsr_block(&block, &pp, &scsr).unwrap().unwrap();
    assert_eq!(
        [
            without.num_atoms(),
            without.num_bonds(),
            without.substance_groups().len() as u32
        ],
        [254, 282, 38]
    );
}
