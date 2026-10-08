"""MolFromMolBlock / MolFromMolFile (V2000). Expected values were checked
against RDKit 2026.09.1."""
import os

import pytest

from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
UPSTREAM_DATA = os.path.join(HERE, "..", "data", "molfiles")

ALANINE = """
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
"""


def test_wedges_and_charges():
    m = Chem.MolFromMolBlock(ALANINE)
    assert Chem.MolToSmiles(m) == "C[C@H]([NH3+])C(=O)O"
    assert m.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    assert not m.GetConformer().Is3D()
    assert m.GetProp("_MolFileInfo") == "     RDKit          2D"
    assert m.GetProp("_MolFileChiralFlag") == "0"
    assert m.GetBondWithIdx(1).GetBondDir() == Chem.BondDir.NONE
    # "  0" fields are skipped, so no molParity etc.
    a = m.GetAtomWithIdx(1)
    assert not a.HasProp("molParity")
    assert a.GetProp("_CIPCode") == "S"


def test_bytes_and_flags():
    m = Chem.MolFromMolBlock(ALANINE.encode())
    assert m.GetNumAtoms() == 6
    m = Chem.MolFromMolBlock(ALANINE, sanitize=False)
    assert m.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    m = Chem.MolFromMolBlock(ALANINE, removeHs=False)
    assert Chem.MolToSmiles(m) == "C[C@H]([NH3+])C(=O)O"


def test_failures():
    assert Chem.MolFromMolBlock("") is None
    assert Chem.MolFromMolBlock(ALANINE.replace("M  END\n", "")) is None
    bad = ALANINE.replace(" N   0", " Qq  0", 1)
    assert Chem.MolFromMolBlock(bad) is None
    m = Chem.MolFromMolBlock(bad, strictParsing=False)
    assert m.GetAtomWithIdx(2).GetAtomicNum() == 0
    assert m.GetAtomWithIdx(2).GetProp("dummyLabel") == "Qq"
    with pytest.raises(OSError):
        Chem.MolFromMolFile("/nonexistent/file.mol")


@pytest.mark.parametrize("fname,smi", [
    ("chiral1.mol", "C[C@](F)(Cl)Br"),
    ("chiral3.mol", "C[C@H](F)Cl"),
    ("list-query.mol", "*1ccccc1"),
    ("rgroups1.mol", "[1*]C1CC1[2*]"),
    ("RingBondQuery.mol", "CC1CCC1"),
    ("bond-query.mol", "C=CC~CC"),
    ("chiral_3h.mol", "[3H][C@@H](C)NC(C)=O"),
])
def test_upstream_files(fname, smi):
    m = Chem.MolFromMolFile(os.path.join(UPSTREAM_DATA, fname))
    assert Chem.MolToSmiles(m) == smi


def test_3d():
    m = Chem.MolFromMolFile(os.path.join(UPSTREAM_DATA, "issue123.mol"))
    assert m.GetNumAtoms() == 23
    assert m.GetConformer().Is3D()
    mh = Chem.MolFromMolFile(os.path.join(UPSTREAM_DATA, "issue123.mol"), removeHs=False)
    assert mh.GetNumAtoms() >= 23


V3K = """
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.299038 -0.750000 0.000000 0
M  V30 2 C 0.000000 0.000000 0.000000 0
M  V30 3 O 0.000000 1.500000 0.000000 0
M  V30 4 C 1.299038 -0.750000 0.000000 0
M  V30 5 N 1.299038 -2.250000 0.000000 0 CHG=1
M  V30 6 Cl 2.598076 -0.000000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 4 5 CFG=1
M  V30 5 1 4 6
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEREL1 ATOMS=(2 2 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
"""


def test_v3000():
    m = Chem.MolFromMolBlock(V3K)
    assert Chem.MolToSmiles(m) == "C[C@H](O)[C@H]([NH3+])Cl"
    groups = m.GetStereoGroups()
    assert len(groups) == 1
    assert groups[0].GetGroupType() == Chem.StereoGroupType.STEREO_OR
    assert [a.GetIdx() for a in groups[0].GetAtoms()] == [1, 3]
    assert groups[0].GetReadId() == 1
    q = Chem.MolFromMolFile(os.path.join(UPSTREAM_DATA, "AtomQuery1.mol"))
    assert Chem.MolToSmarts(q) == "[#6]1(:[#6](:[#6](:[#7]:[#6]:[#6]:1-[#6&h{3-}])-[#6])-[#6])-[#6]"
