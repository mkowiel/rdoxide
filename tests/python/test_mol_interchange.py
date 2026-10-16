"""rdoxide.Chem.rdMolInterchange; expected values come from reference RDKit."""

import pytest
from rdoxide import Chem, rdBase
from rdoxide.Chem import rdMolInterchange as R

HEADER = (
    '{"rdkitjson":{"version":12},"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,'
    '"isotope":0,"stereo":"unspecified"},"bond":{"bo":1,"stereo":"unspecified"}},'
)


def versionless(json):
    return json.replace(rdBase.rdkitVersion, "VER")


def test_write_exact_json_without_rings():
    m = Chem.MolFromSmiles("C[C@H](F)C=CCl")
    m.SetProp("_Name", "nm")
    m.SetProp("foo", "bar")
    m.SetDoubleProp("d", 1.0)
    m.SetIntProp("i", 3)
    assert versionless(R.MolToJSON(m)) == (
        HEADER
        + '"molecules":[{"name":"nm","atoms":[{"impHs":3},{"impHs":1,"stereo":"ccw"},{"z":9},'
        '{"impHs":1},{"impHs":1},{"z":17}],"bonds":[{"atoms":[0,1]},{"atoms":[1,2]},'
        '{"atoms":[1,3]},{"bo":2,"atoms":[3,4]},{"atoms":[4,5]}],'
        '"properties":{"foo":"bar","d":1E0,"i":3},"extensions":[{"name":"rdkitRepresentation",'
        '"formatVersion":2,"toolkitVersion":"VER","cipRanks":[0,2,4,1,3,5],'
        '"cipCodes":[[1,"S"]]}]}]}'
    )


def test_double_format_and_string_escapes():
    m = Chem.MolFromSmiles("C")
    m.SetProp("s", 'a"b\\c\n\x01é/')
    m.SetBoolProp("b", True)
    for name, val in [("a", 1e22), ("b2", 0.1), ("c", 123456789.123), ("d", -0.0)]:
        m.SetDoubleProp(name, val)
    json = R.MolToJSON(m)
    assert '"s":"a\\"b\\\\c\\n\\u0001é/"' in json
    assert '"a":1E22' in json and '"b2":1E-1' in json
    assert '"c":1.23456789123E8' in json and '"d":-0E0' in json
    # booleans have no JSON property form and are skipped
    assert '"b":' not in json.split('"properties"')[1].split('"extensions"')[0].replace('"b2"', "")


def test_header_variants_and_parameters():
    m = Chem.MolFromSmiles("C[C@H](O)C[C@@H](C)F |o1:1,4|")
    json = R.MolToJSON(m)
    assert json.startswith('{"rdkitjson":{"version":12}')
    assert '"stereoGroups":[{"type":"or","id":0,"atoms":[1,4]}]' in json
    ps = R.JSONWriteParameters()
    assert ps.useRDKitExtensions
    ps.useRDKitExtensions = False
    plain = R.MolToJSON(m, ps)
    assert plain.startswith('{"commonchem":{"version":10}')
    assert "stereoGroups" not in plain
    with pytest.raises(AttributeError):
        ps.notAParameter = 1
    pp = R.JSONParseParameters()
    assert (pp.setAromaticBonds, pp.strictValenceCheck, pp.parseConformers) == (True, False, True)
    assert (pp.parseProperties, pp.useHCounts) == (True, True)
    with pytest.raises(AttributeError):
        pp.notAParameter = 1


def test_round_trip_keeps_chemistry():
    smis = ["C[C@H](F)Cl", "F/C=C/F", "[NH4+]", "[13CH3][O-]", "[CH2]", "N#N", "C=C=C"]
    ms = [Chem.MolFromSmiles(s) for s in smis]
    back = R.JSONToMols(R.MolsToJSON(ms))
    assert isinstance(back, tuple)
    assert [Chem.MolToSmiles(x) for x in back] == [Chem.MolToSmiles(x) for x in ms]
    assert [Chem.MolToSmiles(x) for x in R.JSONToMols(R.MolToJSON(ms[0]).encode())] == [
        Chem.MolToSmiles(ms[0])
    ]


def test_round_trip_stereo_groups_and_atropisomers():
    for smi in [
        "C[C@H](F)C[C@@H](Cl)Br |a:1,o2:4|",
        "C[C@H](F)C[C@@H](Cl)Br |&1:1,&2:4|",
        "C1CC(C)C=1c1c(CC)cccc1C |wU:4.3,&1:4|",
        "C1CC(C)C=1c1c(CC)cccc1C |wD:4.3|",
    ]:
        m = Chem.MolFromSmiles(smi)
        (n,) = R.JSONToMols(R.MolToJSON(m))
        assert Chem.MolToCXSmiles(n) == Chem.MolToCXSmiles(m)
    m = Chem.MolFromSmiles("C1CC(C)C=1c1c(CC)cccc1C |wU:4.3|")
    assert '"stereo":"atrop_cw"' in R.MolToJSON(m)


def test_conformers_properties_and_charges():
    from rdoxide.Chem import rdPartialCharges

    m = Chem.MolFromSmiles("CCO")
    c = Chem.Conformer(3)
    c.Set3D(True)
    for i in range(3):
        c.SetAtomPosition(i, (i * 0.5, 0.25 * i * i, -1.0 / (i + 1)))
    m.AddConformer(c)
    rdPartialCharges.ComputeGasteigerCharges(m)
    m.SetProp("note", "x")
    m.SetIntProp("n", 4)
    (n,) = R.JSONToMols(R.MolToJSON(m))
    assert n.GetConformer().Is3D()
    for i in range(3):
        assert tuple(n.GetConformer().GetAtomPosition(i)) == pytest.approx(
            tuple(m.GetConformer().GetAtomPosition(i)), abs=1e-12
        )
    assert [a.GetDoubleProp("_GasteigerCharge") for a in n.GetAtoms()] == pytest.approx(
        [a.GetDoubleProp("_GasteigerCharge") for a in m.GetAtoms()], abs=1e-12
    )
    assert n.GetPropsAsDict() == {"note": "x", "n": 4}
    pp = R.JSONParseParameters()
    pp.parseConformers = False
    pp.parseProperties = False
    (n,) = R.JSONToMols(R.MolToJSON(m), pp)
    assert n.GetNumConformers() == 0 and n.GetPropsAsDict() == {}


def test_queries_round_trip():
    for sma in [
        "C[cH]:[#7+]",
        "[C;R2]=,#!@[$(CO),r6]",
        "[$([OH]C=O),$([O-]C=O)]",
        "[NX3;H2,H1;!$(NC=O)]",
        "[#6;R{2-}]",
        "[C;z{1-2}]",
        "[2H]C",
        "C!@C",
    ]:
        q = Chem.MolFromSmarts(sma)
        json = R.MolToJSON(q)
        assert "rdkitQueries" in json
        pp = R.JSONParseParameters()
        pp.useHCounts = False
        (n,) = R.JSONToMols(json, pp)
        assert Chem.MolToSmarts(n) == Chem.MolToSmarts(q)
    assert '"subquery":' in R.MolToJSON(Chem.MolFromSmarts("[$(CO),r6]"))
    json = R.MolToJSON(Chem.MolFromSmarts("C[cH]:[#7+]"))
    assert '{"impHs":3}' not in json and '"descr":"AtomHCount"' in json


def test_query_atom_mol_block():
    mb = """
  Mrv2102 04092105442D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Q -1.4583 0.4583 0 0
M  V30 2 C -0.1247 1.2283 0 0
M  V30 3 C 1.209 0.4583 0 0
M  V30 4 C 2.5427 1.2283 0 0
M  V30 5 C 1.209 -1.0817 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 4
M  V30 3 5 2 3
M  V30 4 1 3 5 TOPO=2
M  V30 END BOND
M  V30 END CTAB
M  END
"""
    json = R.MolToJSON(Chem.MolFromMolBlock(mb))
    assert '"type":"Q"' in json
    pp = R.JSONParseParameters()
    pp.useHCounts = False
    (n,) = R.JSONToMols(json, pp)
    block = Chem.MolToV3KMolBlock(n)
    assert " Q " in block and "TOPO=2" in block


SRU_BLOCK = """
  Mrv

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  STY  1   1 SRU
M  SAL   1  2   2   3
M  SBL   1  2   1   3
M  SCN  1   1 HT 
M  SMT   1 n
M  SDI   1  4    0.5000    1.0000    0.5000   -1.0000
M  SDI   1  4    2.5000    1.0000    2.5000   -1.0000
M  END
"""


def test_substance_groups_round_trip():
    m = Chem.MolFromMolBlock(SRU_BLOCK)
    json = R.MolToJSON(m)
    assert '"substanceGroups":[{"properties":{' in json
    (n,) = R.JSONToMols(json)
    (a,) = Chem.GetMolSubstanceGroups(m)
    (b,) = Chem.GetMolSubstanceGroups(n)
    # vector valued properties (DATAFIELDS) have no JSON form and are skipped
    expected = {k: v for k, v in a.GetPropsAsDict().items() if k != "DATAFIELDS"}
    assert b.GetPropsAsDict() == expected
    assert list(b.GetAtoms()) == list(a.GetAtoms()) == [1, 2]
    assert list(b.GetBonds()) == list(a.GetBonds()) == [0, 2]

    def flat(sg):
        return [c for br in sg.GetBrackets() for p in br for c in p]

    assert flat(b) == pytest.approx(flat(a))


HEAD = (
    '{"commonchem":{"version":10},"defaults":{"atom":{"stereo":"unspecified"},'
    '"bond":{"stereo":"unspecified","bo":1}},"molecules":[%s]}'
)


@pytest.mark.parametrize(
    "text, message",
    [
        ("{", "incomplete JSON"),
        ("", "incomplete JSON"),
        ("[1]", "JSON should be an object"),
        ('{"a":1}', "missing header in JSON"),
        ('{"commonchem":{"version":11}}', "bad version in JSON"),
        ('{"commonchem":1}', "missing version in JSON"),
        ('{"rdkitjson":{"version":13}}', "bad version in JSON"),
        ('{"commonchem":{"version":10},"molecules":1}', "molecules is not an array"),
        (HEAD % '{"atoms":[{"z":6.0}],"bonds":[]}', "value of z is not an int"),
        (HEAD % '{"atoms":[{"z":6,"stereo":"x"}],"bonds":[]}', "bad stereo value for atom"),
        (HEAD % '{"atoms":[{"z":6}]}', "missing bonds in JSON"),
        (HEAD % '{"atoms":[{"z":6},{"z":6}],"bonds":[{"bo":7,"atoms":[0,1]}]}', "bad bond order"),
        (
            HEAD % '{"atoms":[{"z":6},{"z":6}],"bonds":[{"atoms":[0,1],"stereo":"cis"}]}',
            "without stereoAtoms",
        ),
        (HEAD % '{"atoms":[{"z":6}],"bonds":[],"conformers":[{"dim":4,"coords":[]}]}', "dimension"),
        (
            HEAD % '{"atoms":[{"z":6}],"bonds":[],"stereoGroups":[{"type":"foo","atoms":[0]}]}',
            "bad stereogroup type",
        ),
        (
            HEAD
            % '{"atoms":[{"z":6}],"bonds":[],"substanceGroups":[{"properties":{"TYPE":"FOO"}}]}',
            "TYPE 'FOO' not recognized",
        ),
    ],
)
def test_malformed_input_raises_runtime_error(text, message):
    with pytest.raises(RuntimeError, match=message):
        R.JSONToMols(text)


def test_unusual_arguments():
    assert R.MolsToJSON([]) == ""
    with pytest.raises(ValueError, match="null molecule"):
        R.MolsToJSON([Chem.MolFromSmiles("C"), None])
    assert R.JSONToMols('{"commonchem":{"version":10}}') == ()
    many = [Chem.MolFromSmiles("C" * n) for n in range(1, 60)]
    back = R.JSONToMols(R.MolsToJSON(many))
    assert [m.GetNumAtoms() for m in back] == list(range(1, 60))
