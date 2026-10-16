"""SubstanceGroups in mol blocks; expectations validated against RDKit 2026.09.1."""

from rdoxide import Chem

def v2000(mlines):
    atoms = [f"{x:10.4f}    0.0000    0.0000 {s}   0  0  0  0  0  0  0  0  0  0  0  0"
             for x, s in zip(range(5), "CCOCC")]
    bonds = [f"{i:3d}{i + 1:3d}  1  0" for i in range(1, 5)]
    head = ["", "  test", "", f"{5:3d}{4:3d}  0  0  0  0  0  0  0  0999 V2000"]
    return "\n".join(head + atoms + bonds + mlines + ["M  END", ""])


def v3000(sgroups, bond_types=(1, 1, 1, 1), atom_symbols="CCOCC"):
    atoms = [f"M  V30 {i + 1} {s} {i} 0 0 0" for i, s in enumerate(atom_symbols)]
    bonds = [f"M  V30 {i + 1} {bond_types[i]} {i + 1} {i + 2}" for i in range(4)]
    lines = (["", "  test", "", "  0  0  0     0  0            999 V3000",
              "M  V30 BEGIN CTAB", f"M  V30 COUNTS 5 4 {len(sgroups)} 0 0", "M  V30 BEGIN ATOM"]
             + atoms + ["M  V30 END ATOM", "M  V30 BEGIN BOND"] + bonds
             + ["M  V30 END BOND", "M  V30 BEGIN SGROUP"] + sgroups
             + ["M  V30 END SGROUP", "M  V30 END CTAB", "M  END", ""])
    return "\n".join(lines)


def parse(block, strict=True):
    return Chem.MolFromMolBlock(block, sanitize=False, strictParsing=strict)


def groups(mol):
    return Chem.GetMolSubstanceGroups(mol)


def test_v2000_polymer_group():
    mol = parse(v2000([
        "M  STY  1   1 SRU", "M  SLB  1   1   1", "M  SCN  1   1 HT ",
        "M  SAL   1  2   2   3", "M  SBL   1  2   1   3", "M  SMT   1 n",
        "M  SDI   1  4    1.0000    1.0000    1.0000    2.0000",
        "M  SDI   1  4    3.0000    1.0000    3.0000    2.0000", "M  SST  1   1 ALT",
    ]))
    (group,) = groups(mol)
    props = group.GetPropsAsDict(includePrivate=True)
    assert list(props.pop("DATAFIELDS")) == []
    assert props == {
        "CONNECT": "HT", "ID": 1, "LABEL": "n", "SUBTYPE": "ALT", "TYPE": "SRU", "index": 1,
    }
    assert (list(group.GetAtoms()), list(group.GetBonds())) == ([1, 2], [0, 2])
    brackets = [tuple(tuple(p) for p in b) for b in group.GetBrackets()]
    assert brackets[0][:2] == ((1.0, 1.0, 0.0), (1.0, 2.0, 0.0))


def test_v2000_data_group_joins_data_lines():
    mol = parse(v2000([
        "M  STY  1   1 DAT", "M  SAL   1  1   2",
        "M  SDT   1 FLD                           F                    ",
        "M  SDD   1     0.0000    0.0000    DR    ALL  1       5  ",
        "M  SCD   1 abc", "M  SED   1 def",
    ]))
    (group,) = groups(mol)
    assert list(group.GetStringVectProp("DATAFIELDS")) == ["abcdef"]
    assert group.GetProp("FIELDNAME") == "FLD" and group.GetProp("FIELDTYPE") == "F"


def test_v2000_multiple_group_with_parent():
    mol = parse(v2000([
        "M  STY  2   1 MUL   2 SRU", "M  SAL   1  2   2   3", "M  SMT   1 3",
        "M  SAL   2  1   3", "M  SPL  1   2   1", "M  SNC  1   2   2",
        "M  SBT  1   1   1", "M  SDS EXP  1   1",
    ]))
    first, second = groups(mol)
    assert first.GetPropsAsDict(includePrivate=True)["MULT"] == 3
    assert first.GetProp("BRKTYP") == "PAREN" and first.GetProp("ESTATE") == "E"
    assert (second.GetUnsignedProp("PARENT"), second.GetUnsignedProp("COMPNO")) == (1, 2)


def test_v2000_problems():
    assert parse(v2000(["M  STY  1   1 XXX"])) is None
    assert len(groups(parse(v2000(["M  STY  1   1 XXX"]), strict=False))) == 0
    assert parse(v2000(["M  CRS   1   1"]), strict=False) is None
    # references to unknown groups only produce warnings
    (group,) = groups(parse(v2000(["M  STY  1   1 SRU", "M  SAL   7  1   1", "M  SMT   1 n"])))
    assert list(group.GetAtoms()) == [] and group.GetProp("LABEL") == "n"


def test_v3000_data_group_quoting():
    mol = parse(v3000([
        'M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME="my field" FIELDINFO="a""b" '
        'FIELDDATA="hello world"'
    ]))
    (group,) = groups(mol)
    assert group.GetProp("FIELDNAME") == "my field"
    assert group.GetProp("FIELDINFO") == 'a"b'
    assert list(group.GetStringVectProp("DATAFIELDS")) == ["hello world"]


def test_v3000_polymer_group():
    mol = parse(v3000([
        "M  V30 1 SRU 1 ATOMS=(2 2 3) XBONDS=(2 1 3) BRKXYZ=(9 1 1 0 1 2 0 0 0 0) "
        "BRKXYZ=(9 3 1 0 3 2 0 0 0 0) CONNECT=HT LABEL=n SUBTYPE=ALT"
    ]))
    (group,) = groups(mol)
    assert list(group.GetAtoms()) == [1, 2] and list(group.GetBonds()) == [0, 2]
    assert len(group.GetBrackets()) == 2
    assert group.GetProp("CONNECT") == "HT" and group.GetUnsignedProp("ID") == 1
    mol = parse(v3000([
        "M  V30 1 SRU 0 ATOMS=(2 2 3) XBONDS=(2 1 3) XBHEAD=(1 1) XBCORR=(2 1 3)"
    ]))
    (group,) = groups(mol)
    assert list(group.GetUnsignedVectProp("XBHEAD")) == [0]
    assert list(group.GetUnsignedVectProp("XBCORR")) == [0, 2]


def test_v3000_attachment_points_and_parents():
    mol = parse(v3000(["M  V30 1 SUP 0 ATOMS=(2 2 3) XBONDS=(1 3) SAP=(3 3 aidx 1)"]))
    (group,) = groups(mol)
    assert [(a.aIdx, a.lvIdx, a.id) for a in group.GetAttachPoints()] == [(2, 2, "1")]
    mol = parse(v3000([
        "M  V30 1 MUL 0 ATOMS=(2 2 3) XBONDS=(2 1 3) MULT=3 PATOMS=(1 2)",
        "M  V30 2 GEN 0 ATOMS=(1 1) PARENT=1 COMPNO=2",
    ]))
    first, second = groups(mol)
    assert list(first.GetParentAtoms()) == [1]
    assert second.GetUnsignedProp("PARENT") == 1


def test_v3000_default_labels_fill_in():
    mol = parse(v3000([
        "M  V30 DEFAULT CONNECT=HT LABEL=z", "M  V30 1 SRU 0 ATOMS=(1 1) LABEL=q",
    ]).replace("COUNTS 5 4 2", "COUNTS 5 4 1"))
    (group,) = groups(mol)
    assert group.GetProp("CONNECT") == "HT" and group.GetProp("LABEL") == "q"


def test_v3000_problems():
    for line in ("M  V30 1 XXX 0 ATOMS=(1 1)", "M  V30 1 SRU 0 ATOMS=(7 1 1 1 1 1 1 1)",
                 "M  V30 1 SRU 0 ATOMS=(1 1) SUBTYPE=QQQ"):
        assert parse(v3000([line])) is None
    lenient = parse(v3000(["M  V30 1 XXX 0 ATOMS=(1 1)"]), strict=False)
    assert [g.GetProp("TYPE") for g in groups(lenient)] == ["XXX"]
    assert len(groups(parse(v3000(["M  V30 1 SRU 0 ATOMS=(1 1) CONNECT=QQ"]), strict=False))) == 0


def test_structure_groups_are_applied_and_removed():
    smartsq = 'M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=QUERY QUERYTYPE=SMARTSQ FIELDDATA="[N,O]"'
    mol = parse(v3000([smartsq]))
    assert len(groups(mol)) == 0
    assert Chem.MolToSmarts(mol) == "[#6]-[N,O]-[#8]-[#6]-[#6]"
    assert mol.GetAtomWithIdx(1).GetProp("MRV SMA") == "[N,O]"

    coordinate = "M  V30 1 DAT 0 ATOMS=(2 2 3) FIELDNAME=MRV_COORDINATE_BOND_TYPE FIELDDATA=\"2\""
    mol = parse(v3000([coordinate], bond_types=(1, 0, 1, 1)))
    assert mol.GetBondWithIdx(1).GetBondType() == Chem.BondType.DATIVE

    mol = parse(v3000(['M  V30 1 DAT 0 ATOMS=(2 1 2) FIELDNAME=ZCH FIELDDATA="1;-1"']))
    assert [a.GetFormalCharge() for a in mol.GetAtoms()][:2] == [1, -1]
    mol = parse(v3000(['M  V30 1 DAT 0 ATOMS=(2 1 2) FIELDNAME=HYD FIELDDATA="2;"']))
    assert mol.GetAtomWithIdx(0).GetNumExplicitHs() == 2

    mol = parse(v3000(['M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=other FIELDDATA="x"']))
    assert len(groups(mol)) == 1
