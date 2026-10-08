//! MolBlock reader tests (cf. upstream file_parsers_catch.cpp / test1.cpp).
//! Expected values come from RDKit 2026.09.1; full-corpus comparison is the
//! `molblock` differential mode.

use std::path::PathBuf;

use rdoxide::RDKitError;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::atom::ChiralType;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_block, mol_from_mol_file};
use rdoxide::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use rdoxide::rdgeneral::{RDProps, RDValue};

fn data(f: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../upstream/rdkit/Code/GraphMol/FileParsers/test_data")
        .join(f)
        .to_string_lossy()
        .into_owned()
}

fn read(f: &str) -> RWMol {
    mol_from_mol_file(&data(f), &MolFileParserParams::default())
        .unwrap()
        .unwrap()
}

fn smiles(m: &mut RWMol) -> String {
    mol_to_smiles(m, &SmilesWriteParams::default()).unwrap()
}

#[test]
fn upstream_files() {
    for (f, smi, sma) in [
        (
            "chiral1.mol",
            "C[C@](F)(Cl)Br",
            "[#9]-[#6@](-[#35])(-[#17])-[#6]",
        ),
        ("chiral3.mol", "C[C@H](F)Cl", "[#9]-[#6@H](-[#17])-[#6]"),
        (
            "list-query.mol",
            "*1ccccc1",
            "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6,#7,#15]:1",
        ),
        ("rgroups1.mol", "[1*]C1CC1[2*]", "[#6]1-[#6](-[#6]-1-*)-*"),
        (
            "RingBondQuery.mol",
            "CC1CCC1",
            "[#6]1-[#6]-[#6](-[#6]-1)-&@[#6]",
        ),
        ("bond-query.mol", "C=CC~CC", "[#6]-[#6]-,=[#6]-[#6]=[#6]"),
        (
            "chiral_3h.mol",
            "[3H][C@@H](C)NC(C)=O",
            "[#8]=[#6](-[#7]-[#6@H](-[#6])-[3H])-[#6]",
        ),
    ] {
        let mut m = read(f);
        assert_eq!(smiles(&mut m), smi, "{f}");
        assert_eq!(
            mol_to_smarts(&mut m, &Default::default()).unwrap(),
            sma,
            "{f}"
        );
    }
    let m = read("issue123.mol");
    assert_eq!(m.num_atoms(), 23);
    assert!(m.conformer(-1).unwrap().is_3d());
}

const ALANINE: &str = "
     RDKit          2D

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  6
  2  4  1  0
  4  5  2  0
  4  6  1  0
M  CHG  1   3   1
M  END
";

#[test]
fn wedges_and_charges() {
    let mut m = mol_from_mol_block(ALANINE, &MolFileParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(smiles(&mut m), "C[C@H]([NH3+])C(=O)O");
    assert_eq!(m.atom(1).chiral_tag(), ChiralType::TetrahedralCCW);
    assert!(!m.conformer(-1).unwrap().is_3d());
    assert_eq!(
        m.get_prop("_MolFileInfo"),
        Some(&RDValue::String("     RDKit          2D".into()))
    );
    // CRLF line endings are accepted
    let crlf = ALANINE.replace('\n', "\r\n");
    let mut m2 = mol_from_mol_block(&crlf, &MolFileParserParams::default())
        .unwrap()
        .unwrap();
    assert_eq!(smiles(&mut m2), "C[C@H]([NH3+])C(=O)O");
}

#[test]
fn errors() {
    let p = MolFileParserParams::default();
    // empty input: no molecule
    assert!(mol_from_mol_block("", &p).unwrap().is_none());
    // missing M  END
    let trunc = ALANINE.replace("M  END\n", "");
    assert!(matches!(
        mol_from_mol_block(&trunc, &p),
        Err(RDKitError::FileParse(_))
    ));
    // too few atom lines
    let short: String = ALANINE.lines().take(6).collect::<Vec<_>>().join("\n");
    assert!(matches!(
        mol_from_mol_block(&short, &p),
        Err(RDKitError::FileParse(_))
    ));
    // a bad element symbol only parses in lax mode
    let bad = ALANINE.replacen(" N   0", " Qq  0", 1);
    assert!(mol_from_mol_block(&bad, &p).is_err());
    let lax = MolFileParserParams {
        strict_parsing: false,
        ..Default::default()
    };
    let m = mol_from_mol_block(&bad, &lax).unwrap().unwrap();
    assert_eq!(m.atom(2).atomic_num(), 0);
    assert!(matches!(
        mol_from_mol_file("/nonexistent/file.mol", &p),
        Err(RDKitError::BadFile(_))
    ));
}
