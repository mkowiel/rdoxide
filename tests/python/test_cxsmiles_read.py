"""CXSMILES reading; expectations verified against RDKit 2026.09.1."""

import pytest
from rdkit import Chem


def parse(smi, sanitize=True, strict=True):
    p = Chem.SmilesParserParams()
    p.sanitize = sanitize
    p.strictCXSMILES = strict
    return Chem.MolFromSmiles(smi, p)


def test_coordinates_2d_and_3d():
    m = parse("CCO |(0,0;1.5,0;3,0.5)|")
    conf = m.GetConformer()
    assert not conf.Is3D()
    assert tuple(conf.GetAtomPosition(2)) == (3.0, 0.5, 0.0)
    m = parse("CC |(0,0,0;1,1,1)|")
    assert m.GetConformer().Is3D()
    # all-zero z means 2D even when z was given
    assert not parse("CC |(0,0,0;1,1,0)|").GetConformer().Is3D()
    assert m.GetProp("_CXSMILES_Data") == "|(0,0,0;1,1,1)|"


def test_labels_values_and_atom_props():
    m = parse("CCO |$_AV:;a;b$|")
    assert [a.HasProp("molFileValue") for a in m.GetAtoms()] == [False, True, True]
    assert m.GetAtomWithIdx(1).GetProp("molFileValue") == "a"
    m = parse("C* |$;R1$|")
    assert m.GetAtomWithIdx(1).GetProp("atomLabel") == "R1"
    m = parse("CC |atomProp:0.foo.bar:1.baz.qux|")
    assert m.GetAtomWithIdx(0).GetProp("foo") == "bar"
    assert m.GetAtomWithIdx(1).GetProp("baz") == "qux"
    m = parse("C* |$;a&#59;b$|")
    assert m.GetAtomWithIdx(1).GetProp("atomLabel") == "a;b"


def test_special_labels_become_queries():
    m = parse("C* |$;Q_e$|")
    assert m.GetAtomWithIdx(1).HasQuery()
    assert m.GetAtomWithIdx(1).GetSmarts() == "[!#6&!#1]"
    # without an extension a plain * stays a dummy; with one it is the "A" query
    assert not parse("C*").GetAtomWithIdx(1).HasQuery()
    assert parse("C* |(0,0;1,1)|").GetAtomWithIdx(1).HasQuery()


def test_radicals():
    m = parse("C.C.C |^1:0,2|", sanitize=False)
    assert [a.GetNumRadicalElectrons() for a in m.GetAtoms()] == [1, 0, 1]
    # ^2..^4 are divalent, ^5..^7 trivalent
    for code, count in (("^2", 2), ("^4", 2), ("^5", 3), ("^7", 3)):
        radical = parse(f"CC |{code}:0|", sanitize=False)
        assert radical.GetAtomWithIdx(0).GetNumRadicalElectrons() == count
    assert parse("CC |^8:0|") is None


def test_enhanced_stereo_groups():
    m = parse("C[C@H](F)[C@@H](Cl)Br |&1:1,3|")
    (g,) = m.GetStereoGroups()
    assert g.GetGroupType() == Chem.StereoGroupType.STEREO_AND
    assert [a.GetIdx() for a in g.GetAtoms()] == [1, 3]
    m = parse("C[C@H](F)Cl |a:1|")
    assert m.GetStereoGroups()[0].GetGroupType() == Chem.StereoGroupType.STEREO_ABSOLUTE
    # groups with the same type and id are merged
    m = parse("CC[C@H](F)[C@H](Cl)Br |o1:2,o1:4|")
    (g,) = m.GetStereoGroups()
    assert [a.GetIdx() for a in g.GetAtoms()] == [2, 4]
    # a repeated atom invalidates the extension
    assert parse("C[C@H](F)Cl |o:1,1|") is None


def test_double_bond_stereo():
    smi = "C1CCCC/C=C/CCC1"
    assert parse(smi + " |ctu:5|").GetBondWithIdx(5).GetStereo() == Chem.BondStereo.STEREOANY
    assert parse("C1CCCCC=CCCC1 |c:5|").GetBondWithIdx(5).GetStereo() == Chem.BondStereo.STEREOZ
    assert parse("C1CCCCC=CCCC1 |t:5|").GetBondWithIdx(5).GetStereo() == Chem.BondStereo.STEREOE


def test_wedged_bonds():
    m = parse("CC(O)Cl |wU:1.0|")
    assert m.GetBondWithIdx(0).GetIntProp("_MolFileBondCfg") == 1
    m = parse("CC(O)Cl |w:1.0|", sanitize=False)
    assert m.GetBondWithIdx(0).GetBondDir() == Chem.BondDir.UNKNOWN
    assert m.GetBondWithIdx(0).GetIntProp("_MolFileBondCfg") == 2


def test_coordinate_hydrogen_and_zero_order_bonds():
    m = parse("[Fe]N |C:1.0|")
    assert m.GetBondWithIdx(0).GetBondType() == Chem.BondType.DATIVE
    assert m.GetBondWithIdx(0).GetBeginAtomIdx() == 1
    assert parse("OC |H:0.0|").GetBondWithIdx(0).GetBondType() == Chem.BondType.HYDROGEN
    m = parse("C1CCCCC~CCCC1 |Z:5|")
    assert m.GetBondWithIdx(5).GetBondType() == Chem.BondType.ZERO


def test_query_features():
    m = parse("CC(C)C |u:1|")
    assert m.GetAtomWithIdx(1).HasQuery()
    m = parse("CC(C)C |rb:1:2|")
    assert m.GetAtomWithIdx(1).HasQuery()
    m = parse("CC(C)C |s:1:*|")
    assert m.GetAtomWithIdx(1).HasQuery()
    assert parse("CC(C)C |rb:1:5|") is None


def test_names_and_failures():
    assert parse("CCO ethanol").GetProp("_Name") == "ethanol"
    m = parse("CCO |$_AV:;;x$| ethanol")
    assert m.GetProp("_Name") == "ethanol"
    assert m.GetProp("_CXSMILES_Data") == "|$_AV:;;x$|"
    # extra entries beyond the atoms are ignored
    assert parse("CC |(0,0;1,1;2,2)|") is not None
    assert parse("CCO |$a;b;c;d$|") is not None
    # strict parsing refuses a broken extension, lenient parsing keeps the molecule
    assert parse("CC |o:0,0|") is None
    lenient = parse("CC |o:0,0|", strict=False)
    assert lenient is not None
    assert lenient.GetProp("_CXSMILES_Data").startswith("|")
    assert not lenient.HasProp("_Name")


def test_smarts_with_cx_extensions():
    p = Chem.SmartsParserParams()
    q = Chem.MolFromSmarts("[#6][#6] |(0,0;1,1)|", p)
    assert q.GetNumConformers() == 1


def test_any_bond_in_smiles():
    m = parse("C~C", sanitize=False)
    assert m.GetBondWithIdx(0).HasQuery()
