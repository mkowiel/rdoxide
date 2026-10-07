"""Core rdchem object model (Atom, Bond, Mol, RWMol, EditableMol, Conformer,
RingInfo, StereoGroup, enums, props). Cases adapted from upstream
Code/GraphMol/Wrap/rough_test.py and testConformer.py (BSD-3-Clause) that do
not need SMILES parsing; the full upstream files run via
tools/run_upstream_tests.py."""

import copy
import pickle

import pytest

from rdkit import Chem, rdBase
from rdkit.Geometry import Point3D


def make_mol(nums, bonds):
    m = Chem.RWMol()
    for n in nums:
        m.AddAtom(Chem.Atom(n))
    for b, e, t in bonds:
        m.AddBond(b, e, t)
    return m


@pytest.fixture
def ethanol():
    S = Chem.BondType.SINGLE
    m = make_mol([6, 6, 8], [(0, 1, S), (1, 2, S)])
    m.UpdatePropertyCache()
    return m


def test_wrapper_type():
    assert rdBase._wrapperType == "nanobind"


def test_enums():
    bt = Chem.BondType
    assert bt.SINGLE == 1 and isinstance(bt.SINGLE, int)
    assert str(bt.AROMATIC) == "AROMATIC"
    assert repr(bt.DOUBLE) == "rdkit.Chem.rdchem.BondType.DOUBLE"
    assert bt.values[12] is bt.AROMATIC
    assert bt.names["ZERO"] == 21
    assert bt(2) is bt.DOUBLE
    assert pickle.loads(pickle.dumps(bt.TRIPLE)) is bt.TRIPLE
    assert Chem.ChiralType.CHI_TETRAHEDRAL_CCW == 2
    assert Chem.HybridizationType.SP3 == 4
    assert Chem.BondStereo.STEREOE == 3
    assert Chem.BondDir.ENDUPRIGHT == 4
    assert Chem.ValenceType.EXPLICIT == 1
    assert Chem.rdchem.StereoGroupType.STEREO_AND == 2


def test_atom_basics():
    a = Chem.Atom(6)
    assert a.GetSymbol() == "C" and a.GetAtomicNum() == 6
    assert Chem.Atom("Cl").GetAtomicNum() == 17
    assert not a.HasOwningMol()
    assert a.GetDegree() == 0
    a.SetFormalCharge(-1)
    a.SetIsotope(13)
    assert a.GetFormalCharge() == -1 and a.GetIsotope() == 13
    assert abs(a.GetMass() - 13.003354838) < 1e-6
    b = Chem.Atom(a)
    assert b.GetIsotope() == 13
    c = copy.copy(a)
    assert c.GetFormalCharge() == -1
    with pytest.raises(RuntimeError):
        a.GetOwningMol()
    a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    assert a.InvertChirality()
    assert a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    a.SetAtomMapNum(7)
    assert a.GetAtomMapNum() == 7 and a.GetIntProp("molAtomMapNumber") == 7
    a.SetAtomMapNum(0)
    assert not a.HasProp("molAtomMapNumber")
    assert Chem.Atom.NOATOM == 2**32 - 1
    with pytest.raises(RuntimeError):
        Chem.Atom("Xx")


def test_atom_in_mol(ethanol):
    m = ethanol
    a = m.GetAtomWithIdx(1)
    assert a.GetIdx() == 1 and a.GetDegree() == 2
    assert [n.GetIdx() for n in a.GetNeighbors()] == [0, 2]
    assert [b.GetIdx() for b in a.GetBonds()] == [0, 1]
    assert a.GetTotalNumHs() == 2 and a.GetNumImplicitHs() == 2
    assert a.GetValence(Chem.ValenceType.EXPLICIT) == 2
    assert a.GetValence(Chem.ValenceType.IMPLICIT) == 2
    assert a.GetTotalValence() == 4 and a.GetTotalDegree() == 4
    assert not a.HasValenceViolation()
    a.SetFormalCharge(1)
    assert m.GetAtomWithIdx(1).GetFormalCharge() == 1
    assert a.GetOwningMol() is m
    assert m.GetNumAtoms(onlyExplicit=False) == 9
    assert m.GetNumHeavyAtoms() == 3


def test_bonds(ethanol):
    m = ethanol
    b = m.GetBondBetweenAtoms(2, 1)
    assert b.GetIdx() == 1 and b.GetBeginAtomIdx() == 1 and b.GetEndAtomIdx() == 2
    assert b.GetOtherAtomIdx(1) == 2
    assert b.GetOtherAtom(m.GetAtomWithIdx(2)).GetIdx() == 1
    assert b.GetBondType() == Chem.BondType.SINGLE
    assert b.GetBondTypeAsDouble() == 1.0
    assert b.GetValenceContrib(m.GetAtomWithIdx(1)) == 1.0
    assert b.GetValenceContrib(m.GetAtomWithIdx(0)) == 0.0
    b.SetBondType(Chem.BondType.DOUBLE)
    assert m.GetBondWithIdx(1).GetBondType() == Chem.BondType.DOUBLE
    assert m.GetBondBetweenAtoms(0, 2) is None
    with pytest.raises(RuntimeError):
        m.GetBondWithIdx(5)
    with pytest.raises(RuntimeError):
        b.SetStereo(Chem.BondStereo.STEREOCIS)


def test_sequences(ethanol):
    atoms = ethanol.GetAtoms()
    assert len(atoms) == 3
    assert [a.GetSymbol() for a in atoms] == ["C", "C", "O"]
    assert atoms[-1].GetSymbol() == "O"
    with pytest.raises(IndexError):
        atoms[3]
    assert len(ethanol.GetBonds()) == 2
    it = iter(ethanol.GetAtoms())
    next(it)
    ethanol.AddAtom(Chem.Atom(7))
    with pytest.raises(RuntimeError):
        next(it)


def test_props(ethanol):
    m = ethanol
    m.SetProp("a", "1.5")
    m.SetIntProp("i", 3)
    m.SetDoubleProp("d", 0.1)
    m.SetBoolProp("b", True)
    m.SetUnsignedProp("u", 5, computed=True)
    m.SetProp("_priv", "x")
    assert m.GetProp("i") == "3"
    assert m.GetProp("d") == "0.10000000000000001"
    assert m.GetProp("b") == "1"
    assert m.GetProp("i", autoConvert=True) == 3
    assert m.GetProp("a", autoConvert=True) == 1.5
    assert m.GetDoubleProp("a") == 1.5
    with pytest.raises(ValueError):
        m.GetIntProp("a")
    with pytest.raises(ValueError):
        m.GetDoubleProp("i")
    with pytest.raises(KeyError):
        m.GetProp("missing")
    assert m.GetProp("missing", default=None) is None
    assert m.GetIntProp("missing", default=7) == 7
    assert m.GetIntProp("missing", 8) == 8
    assert m.GetProp("missing", False, "dflt") == "dflt"
    assert m.GetProp("i", True, None) == 3
    assert list(m.GetPropNames()) == ["a", "i", "d", "b"]
    assert "u" in m.GetPropNames(includeComputed=True)
    assert "_priv" in m.GetPropNames(includePrivate=True)
    assert m.GetPropsAsDict() == {"a": 1.5, "i": 3, "d": 0.1, "b": True}
    m.ClearComputedProps()
    assert not m.HasProp("u")
    m.ClearProp("a")
    m.ClearProp("does not exist")
    assert not m.HasProp("a")
    m.SetName("ethanol")
    assert m.GetProp("_Name") == "ethanol" and m.GetName() == "ethanol"
    a = m.GetAtomWithIdx(0)
    a.SetProp("foo", "bar")
    assert m.GetAtomWithIdx(0).GetProp("foo") == "bar"
    assert a.GetPropsAsDict() == {"foo": "bar"}


def test_rwmol_editing():
    S = Chem.BondType.SINGLE
    m = make_mol([6, 7, 8, 9], [(0, 1, S), (1, 2, S), (2, 3, S)])
    m.RemoveAtom(1)
    assert m.GetNumAtoms() == 3 and m.GetNumBonds() == 1
    b = m.GetBondWithIdx(0)
    assert (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) == (1, 2)
    m.RemoveBond(1, 2)
    assert m.GetNumBonds() == 0
    m.ReplaceAtom(0, Chem.Atom(17))
    assert m.GetAtomWithIdx(0).GetSymbol() == "Cl"
    m2 = make_mol([6, 6], [(0, 1, S)])
    m.InsertMol(m2)
    assert m.GetNumAtoms() == 5 and m.GetBondWithIdx(0).GetBeginAtomIdx() == 3
    with Chem.RWMol(m) as rw:
        rw.RemoveAtom(0)
        rw.RemoveAtom(1)
        assert rw.GetNumAtoms() == 5
    assert rw.GetNumAtoms() == 3
    mol = rw.GetMol()
    assert type(mol) is Chem.Mol and mol.GetNumAtoms() == 3
    assert type(copy.deepcopy(rw)) is Chem.RWMol


def test_editable_mol():
    em = Chem.EditableMol(Chem.Mol())
    for n in (6, 6, 8):
        em.AddAtom(Chem.Atom(n))
    assert em.AddBond(0, 1, order=Chem.BondType.DOUBLE) == 1
    em.AddBond(1, 2)
    m = em.GetMol()
    assert m.GetNumAtoms() == 3
    assert m.GetBondWithIdx(1).GetBondType() == Chem.BondType.UNSPECIFIED


def test_conformers(ethanol):
    m = ethanol
    conf = Chem.Conformer(3)
    for i in range(3):
        conf.SetAtomPosition(i, Point3D(float(i), 0.0, 0.0))
    conf.SetAtomPosition(2, (2.0, 1.0, 0.0))
    cid = m.AddConformer(conf, assignId=True)
    assert cid == 0 and m.GetNumConformers() == 1
    c = m.GetConformer()
    assert c.HasOwningMol() and c.GetOwningMol() is m
    p = c.GetAtomPosition(2)
    assert (p.x, p.y, p.z) == (2.0, 1.0, 0.0)
    pos = c.GetPositions()
    assert pos.shape == (3, 3)
    c.SetAtomPosition(0, Point3D(5, 5, 5))
    assert m.GetConformer(0).GetAtomPosition(0).x == 5
    with pytest.raises(ValueError):
        m.GetConformer(7)
    assert len(m.GetConformers()) == 1
    m.RemoveAllConformers()
    assert m.GetNumConformers() == 0
    with pytest.raises(RuntimeError):
        m.AddConformer(Chem.Conformer(2))


def test_ring_info_manual(ethanol):
    ri = ethanol.GetRingInfo()
    with pytest.raises(RuntimeError):
        ri.NumRings()
    ri.AddRing((0, 1, 2), (0, 1, 2))
    assert ri.NumRings() == 1
    assert ri.AtomRings() == ((0, 1, 2),)
    assert ri.NumAtomRings(1) == 1 and ri.IsAtomInRingOfSize(1, 3)


def test_stereo_groups(ethanol):
    g = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, ethanol, [1], readId=4)
    ethanol.SetStereoGroups([g])
    groups = ethanol.GetStereoGroups()
    assert len(groups) == 1
    assert groups[0].GetGroupType() == Chem.rdchem.StereoGroupType.STEREO_OR
    assert [a.GetIdx() for a in groups[0].GetAtoms()] == [1]
    assert groups[0].GetReadId() == 4


def test_sanitize_exceptions():
    m = Chem.RWMol()
    for _ in range(6):
        m.AddAtom(Chem.Atom(6))
    for i in range(1, 6):
        m.AddBond(0, i, Chem.BondType.SINGLE)
    with pytest.raises(Chem.AtomValenceException) as ei:
        m.UpdatePropertyCache()
    assert isinstance(ei.value, Chem.AtomSanitizeException)
    assert isinstance(ei.value, ValueError)
    assert ei.value.cause.GetAtomIdx() == 0
    assert ei.value.cause.GetType() == "AtomValenceException"
    m.UpdatePropertyCache(strict=False)


def test_periodic_table_via_rdchem():
    assert Chem.rdchem.GetPeriodicTable().GetElementSymbol(7) == "N"


def test_point3d():
    p = Point3D(1, 2, 3)
    q = Point3D(1, 0, 0)
    assert len(p) == 3 and p[0] == 1 and p[-1] == 3
    assert (p + q).x == 2 and (p - q).x == 0 and (p * 2).z == 6
    assert p.CrossProduct(q).y == 3
    assert abs(q.Distance(Point3D(0, 0, 0)) - 1) < 1e-12
    r = pickle.loads(pickle.dumps(p))
    assert (r.x, r.y, r.z) == (1, 2, 3)
