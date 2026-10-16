"""SubstanceGroup API; expectations validated against RDKit 2026.09.1."""

import pytest

from rdoxide import Chem
from rdoxide.Geometry import Point3D


def sru(smiles="CCOCC", atoms=(1, 2), bonds=(0, 3)):
    mol = Chem.MolFromSmiles(smiles)
    group = Chem.CreateMolSubstanceGroup(mol, "SRU")
    for a in atoms:
        group.AddAtomWithIdx(a)
    for b in bonds:
        group.AddBondWithIdx(b)
    return mol, group


def summary(mol):
    return [
        (g.GetProp("TYPE"), list(g.GetAtoms()), list(g.GetBonds()), list(g.GetParentAtoms()))
        for g in Chem.GetMolSubstanceGroups(mol)
    ]


def test_create_and_query():
    mol, group = sru()
    assert group.GetIndexInMol() == 0
    assert list(group.GetAtoms()) == [1, 2]
    assert list(group.GetBonds()) == [0, 3]
    assert list(group.GetParentAtoms()) == []
    assert group.GetProp("TYPE") == "SRU"
    assert Chem.GetMolSubstanceGroupWithIdx(mol, 0).GetIndexInMol() == 0
    with pytest.raises(IndexError):
        Chem.GetMolSubstanceGroupWithIdx(mol, 3)


def test_copies_are_detached():
    mol, _ = sru()
    copies = Chem.GetMolSubstanceGroups(mol)
    assert len(copies) == 1
    with pytest.raises(RuntimeError, match="Unable to find own index"):
        copies[0].GetIndexInMol()
    copies[0].SetProp("note", "x")
    assert not Chem.GetMolSubstanceGroupWithIdx(mol, 0).HasProp("note")
    # list items persist
    assert copies[0].HasProp("note")
    added = Chem.AddMolSubstanceGroup(mol, copies[0])
    assert added.GetIndexInMol() == 1
    assert added.HasProp("note")
    assert len(Chem.GetMolSubstanceGroups(mol)) == 2


def test_index_validation():
    _, group = sru()
    with pytest.raises(ValueError, match="Atom index out of range"):
        group.AddAtomWithIdx(9)
    with pytest.raises(ValueError, match="Bond index out of range"):
        group.AddBondWithIdx(9)
    with pytest.raises(RuntimeError, match="not a member"):
        group.AddParentAtomWithIdx(0)
    group.AddParentAtomWithIdx(1)
    assert list(group.GetParentAtoms()) == [1]
    with pytest.raises(ValueError, match="Atom index 20 is out of range"):
        group.SetAtoms([0, 1, 20])
    with pytest.raises(TypeError):
        group.SetAtoms("x")
    group.SetAtoms([0, 1, 2])
    group.SetBonds((1, 2))
    assert (list(group.GetAtoms()), list(group.GetBonds())) == ([0, 1, 2], [1, 2])


def test_cstates_brackets_attach_points():
    _, group = sru()
    with pytest.raises(RuntimeError, match="Pre-condition Violation"):
        group.AddCState(1, Point3D(1, 2, 3))
    with pytest.raises(RuntimeError, match="Neither beginning nor ending"):
        group.AddCState(3, Point3D(1, 2, 3))
    group.AddCState(0, Point3D(1, 2, 3))
    (state,) = group.GetCStates()
    assert state.bondIdx == 0 and tuple(state.vector) == (1.0, 2.0, 3.0)

    group.AddAttachPoint(1, -1, "1")
    group.AddAttachPoint(2, 3, "2")
    points = [(a.aIdx, a.lvIdx, a.id) for a in group.GetAttachPoints()]
    assert points == [(1, -1, "1"), (2, 3, "2")]

    group.AddBracket([Point3D(0, 0, 0), Point3D(1, 1, 0)])
    group.AddBracket([Point3D(0, 0, 0), Point3D(1, 1, 0), Point3D(2, 2, 2)])
    with pytest.raises(ValueError, match="length of 2 or 3"):
        group.AddBracket([Point3D(0, 0, 0)])
    brackets = [tuple(tuple(p) for p in b) for b in group.GetBrackets()]
    assert brackets[0][2] == (0.0, 0.0, 0.0) and brackets[1][2] == (2.0, 2.0, 2.0)

    group.ClearBrackets()
    group.ClearCStates()
    group.ClearAttachPoints()
    assert (len(group.GetBrackets()), len(group.GetCStates()), len(group.GetAttachPoints())) == (0, 0, 0)


def test_properties():
    _, group = sru()
    group.SetIntProp("ID", 3)
    group.SetUnsignedProp("u", 5)
    group.SetDoubleProp("d", 1.5)
    group.SetBoolProp("b", True)
    assert group.HasProp("ID") is True and group.HasProp("zz") is False
    assert (group.GetIntProp("ID"), group.GetUnsignedProp("u")) == (3, 5)
    assert (group.GetDoubleProp("d"), group.GetBoolProp("b")) == (1.5, True)
    assert group.GetPropsAsDict() == {"TYPE": "SRU", "ID": 3, "u": 5, "d": 1.5, "b": True}
    with pytest.raises(KeyError):
        group.GetProp("zz")
    assert group.GetProp("zz", False, "dflt") == "dflt"
    with pytest.raises(KeyError):
        group.GetUnsignedVectProp("zz")
    group.ClearProp("d")
    assert not group.HasProp("d")


def test_data_group():
    mol = Chem.MolFromSmiles("CC")
    group = Chem.CreateMolDataSubstanceGroup(mol, "fld", "val")
    assert group.GetProp("TYPE") == "DAT"
    assert group.GetProp("FIELDNAME") == "fld"
    assert list(group.GetStringVectProp("DATAFIELDS")) == ["val"]


def test_groups_follow_molecule_edits():
    mol, _ = sru()
    assert summary(Chem.Mol(mol)) == [("SRU", [1, 2], [0, 3], [])]
    assert summary(Chem.RWMol(mol).GetMol()) == [("SRU", [1, 2], [0, 3], [])]
    # renumbering does not carry the groups along
    assert summary(Chem.RenumberAtoms(mol, [4, 3, 2, 1, 0])) == []

    for atom in (0, 4):
        rw = Chem.RWMol(mol)
        rw.RemoveAtom(atom)
        assert summary(rw) == []
    rw = Chem.RWMol(mol)
    rw.RemoveBond(0, 1)
    assert summary(rw) == []
    rw = Chem.RWMol(mol)
    rw.RemoveBond(1, 2)
    assert summary(rw) == [("SRU", [1, 2], [0, 2], [])]


def test_groups_follow_combination_and_fragments():
    mol, _ = sru("CCOCC.CCN")
    both = Chem.CombineMols(mol, mol)
    assert summary(both) == [("SRU", [1, 2], [0, 3], []), ("SRU", [9, 10], [6, 9], [])]
    frags = [summary(f) for f in Chem.GetMolFrags(mol, asMols=True)]
    assert frags == [[("SRU", [1, 2], [0, 3], [])], []]
    rw = Chem.RWMol(mol)
    rw.InsertMol(mol)
    assert len(Chem.GetMolSubstanceGroups(rw)) == 2


def test_remove_hs_drops_removed_atoms_from_groups():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCOCC"))
    group = Chem.CreateMolSubstanceGroup(mol, "SRU")
    for a in (1, 2, 3, 4, 10):
        group.AddAtomWithIdx(a)
    group.AddBondWithIdx(0)
    assert summary(Chem.RemoveHs(mol)) == [("SRU", [1, 2, 3, 4], [0], [])]


def test_clear():
    mol, group = sru()
    Chem.ClearMolSubstanceGroups(mol)
    assert len(Chem.GetMolSubstanceGroups(mol)) == 0
    with pytest.raises(RuntimeError, match="No SubstanceGroups found"):
        group.GetIndexInMol()
