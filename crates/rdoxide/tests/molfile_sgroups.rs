//! SubstanceGroups in mol blocks; expectations from RDKit 2026.09.1.

use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::bond::BondType;
use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_block};
use rdoxide::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;
use rdoxide::rdgeneral::{RDProps, RDValue};

fn v2000(m_lines: &[&str]) -> String {
    let atoms = (0..5).map(|x| {
        format!(
            "{x:10.4}    0.0000    0.0000 {}   0  0  0  0  0  0  0  0  0  0  0  0",
            "CCOCC".as_bytes()[x] as char
        )
    });
    let bonds = (1..5).map(|i| format!("{i:3}{:3}  1  0", i + 1));
    let mut lines = vec![
        "".to_string(),
        "  test".into(),
        "".into(),
        "  5  4  0  0  0  0  0  0  0  0999 V2000".into(),
    ];
    lines.extend(atoms.chain(bonds));
    lines.extend(m_lines.iter().map(|l| l.to_string()));
    lines.push("M  END".into());
    lines.join("\n") + "\n"
}

fn v3000(sgroups: &[&str], bond_types: [u8; 4]) -> String {
    let mut lines: Vec<String> = [
        "",
        "  test",
        "",
        "  0  0  0     0  0            999 V3000",
        "M  V30 BEGIN CTAB",
    ]
    .map(String::from)
    .into();
    lines.push(format!("M  V30 COUNTS 5 4 {} 0 0", sgroups.len()));
    lines.push("M  V30 BEGIN ATOM".into());
    lines.extend((0..5).map(|i| {
        format!(
            "M  V30 {} {} {i} 0 0 0",
            i + 1,
            "CCOCC".as_bytes()[i] as char
        )
    }));
    lines.push("M  V30 END ATOM".into());
    lines.push("M  V30 BEGIN BOND".into());
    lines.extend((0..4).map(|i| format!("M  V30 {} {} {} {}", i + 1, bond_types[i], i + 1, i + 2)));
    lines.push("M  V30 END BOND".into());
    lines.push("M  V30 BEGIN SGROUP".into());
    lines.extend(sgroups.iter().map(|l| l.to_string()));
    lines.extend(["M  V30 END SGROUP", "M  V30 END CTAB", "M  END"].map(String::from));
    lines.join("\n") + "\n"
}

fn try_parse(block: &str, strict: bool) -> rdoxide::Result<Option<RWMol>> {
    let params = MolFileParserParams {
        sanitize: false,
        strict_parsing: strict,
        ..MolFileParserParams::default()
    };
    mol_from_mol_block(block, &params)
}

fn parse(block: &str, strict: bool) -> Option<RWMol> {
    try_parse(block, strict).unwrap()
}

fn prop(group: &rdoxide::graphmol::substance_group::SubstanceGroup, key: &str) -> Option<String> {
    group.get_prop(key).map(RDValue::to_rd_string)
}

#[test]
fn v2000_group_lines_are_combined() {
    let block = v2000(&[
        "M  STY  1   1 SRU",
        "M  SLB  1   1   1",
        "M  SCN  1   1 HT ",
        "M  SAL   1  2   2   3",
        "M  SBL   1  2   1   3",
        "M  SMT   1 n",
        "M  SDI   1  4    1.0000    1.0000    1.0000    2.0000",
        "M  SST  1   1 ALT",
    ]);
    let mol = parse(&block, true).unwrap();
    let group = &mol.substance_groups()[0];
    assert_eq!((group.atoms(), group.bonds()), (&[1, 2][..], &[0, 2][..]));
    assert_eq!(group.brackets().len(), 1);
    for (key, value) in [
        ("CONNECT", "HT"),
        ("LABEL", "n"),
        ("SUBTYPE", "ALT"),
        ("TYPE", "SRU"),
    ] {
        assert_eq!(prop(group, key).as_deref(), Some(value), "{key}");
    }
    assert_eq!(group.get_prop("ID"), Some(&RDValue::UnsignedInt(1)));
}

#[test]
fn v2000_data_lines_are_joined() {
    let block = v2000(&[
        "M  STY  1   1 DAT",
        "M  SAL   1  1   2",
        "M  SDT   1 FLD                           F                    ",
        "M  SCD   1 abc",
        "M  SED   1 def",
    ]);
    let mol = parse(&block, true).unwrap();
    let group = &mol.substance_groups()[0];
    assert_eq!(
        group.get_prop("DATAFIELDS"),
        Some(&RDValue::VecString(vec!["abcdef".into()]))
    );
    assert_eq!(prop(group, "FIELDNAME").as_deref(), Some("FLD"));
}

#[test]
fn unsupported_group_types_skip_the_molecule_or_the_group() {
    let block = v2000(&["M  STY  1   1 XXX"]);
    assert!(parse(&block, true).is_none());
    assert!(parse(&block, false).unwrap().substance_groups().is_empty());
}

#[test]
fn v3000_labels_quotes_and_arrays() {
    let block = v3000(
        &[
            "M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=\"my field\" FIELDINFO=\"a\"\"b\" FIELDDATA=\"hello world\"",
        ],
        [1; 4],
    );
    let mol = parse(&block, true).unwrap();
    let group = &mol.substance_groups()[0];
    assert_eq!(prop(group, "FIELDNAME").as_deref(), Some("my field"));
    assert_eq!(prop(group, "FIELDINFO").as_deref(), Some("a\"b"));
    assert_eq!(
        group.get_prop("DATAFIELDS"),
        Some(&RDValue::VecString(vec!["hello world".into()]))
    );
    assert_eq!(group.atoms(), [1]);

    let block = v3000(
        &["M  V30 1 SRU 0 ATOMS=(2 2 3) XBONDS=(2 1 3) XBHEAD=(1 1) XBCORR=(2 1 3)"],
        [1; 4],
    );
    let mol = parse(&block, true).unwrap();
    let group = &mol.substance_groups()[0];
    assert_eq!(
        group.get_prop("XBCORR"),
        Some(&RDValue::VecUnsignedInt(vec![0, 2]))
    );
}

#[test]
fn v3000_problems_depend_on_strictness() {
    // an unsupported type skips the molecule, the other problems are errors
    assert!(parse(&v3000(&["M  V30 1 XXX 0 ATOMS=(1 1)"], [1; 4]), true).is_none());
    for line in [
        "M  V30 1 SRU 0 ATOMS=(7 1 1 1 1 1 1 1)",
        "M  V30 1 SRU 0 ATOMS=(1 1) SUBTYPE=QQQ",
    ] {
        assert!(try_parse(&v3000(&[line], [1; 4]), true).is_err(), "{line}");
    }
    let lenient = parse(
        &v3000(&["M  V30 1 SRU 0 ATOMS=(1 1) CONNECT=QQ"], [1; 4]),
        false,
    )
    .unwrap();
    assert!(lenient.substance_groups().is_empty());
}

#[test]
fn structure_groups_are_applied_and_removed() {
    let smartsq =
        "M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=QUERY QUERYTYPE=SMARTSQ FIELDDATA=\"[N,O]\"";
    let mut mol = parse(&v3000(&[smartsq], [1; 4]), true).unwrap();
    assert!(mol.substance_groups().is_empty());
    assert_eq!(
        mol_to_smarts(&mut mol, &SmilesWriteParams::default()).unwrap(),
        "[#6]-[N,O]-[#8]-[#6]-[#6]"
    );

    let coordinate =
        "M  V30 1 DAT 0 ATOMS=(2 2 3) FIELDNAME=MRV_COORDINATE_BOND_TYPE FIELDDATA=\"2\"";
    let mol = parse(&v3000(&[coordinate], [1, 0, 1, 1]), true).unwrap();
    assert_eq!(mol.bond(1).bond_type(), BondType::Dative);

    let zch = "M  V30 1 DAT 0 ATOMS=(2 1 2) FIELDNAME=ZCH FIELDDATA=\"1;-1\"";
    let mol = parse(&v3000(&[zch], [1; 4]), true).unwrap();
    assert_eq!(
        (mol.atom(0).formal_charge(), mol.atom(1).formal_charge()),
        (1, -1)
    );

    let other = "M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=other FIELDDATA=\"x\"";
    assert_eq!(
        parse(&v3000(&[other], [1; 4]), true)
            .unwrap()
            .substance_groups()
            .len(),
        1
    );
}
