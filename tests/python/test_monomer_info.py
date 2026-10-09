"""AtomMonomerInfo / AtomPDBResidueInfo (rdchem) and their pickling."""

import pickle

import pytest
from rdkit import Chem


def test_defaults():
    assert Chem.AtomMonomerInfo().GetMonomerType() == Chem.AtomMonomerType.UNKNOWN
    pdb = Chem.AtomPDBResidueInfo()
    assert pdb.GetMonomerType() == Chem.AtomMonomerType.PDBRESIDUE
    assert pdb.GetOccupancy() == 1.0
    assert pdb.GetSerialNumber() == 0
    assert Chem.AtomPDBResidueInfo("CA").GetSerialNumber() == 1


def test_isinstance_and_enum():
    assert issubclass(Chem.AtomPDBResidueInfo, Chem.AtomMonomerInfo)
    assert set(Chem.AtomMonomerType.names) == {"UNKNOWN", "PDBRESIDUE", "OTHER"}


def test_atom_holds_live_reference():
    a = Chem.Atom(6)
    assert a.GetMonomerInfo() is None
    assert a.GetPDBResidueInfo() is None
    info = Chem.AtomPDBResidueInfo(" CA ", serialNumber=5, chainId="A")
    a.SetMonomerInfo(info)
    assert isinstance(a.GetMonomerInfo(), Chem.AtomPDBResidueInfo)
    a.GetPDBResidueInfo().SetChainId("Z")
    assert a.GetMonomerInfo().GetChainId() == "Z"
    # the atom keeps a copy
    assert info.GetChainId() == "A"


def test_non_pdb_info():
    a = Chem.Atom(6)
    m = Chem.AtomMonomerInfo(Chem.AtomMonomerType.OTHER, "X", "RES", 3, "C", "LGRP")
    a.SetMonomerInfo(m)
    got = a.GetMonomerInfo()
    assert not isinstance(got, Chem.AtomPDBResidueInfo)
    assert (got.GetName(), got.GetResidueName(), got.GetResidueNumber()) == ("X", "RES", 3)
    with pytest.raises(ValueError, match="not a PDB Residue"):
        a.GetPDBResidueInfo()
    with pytest.raises(ValueError, match="not a PDB Residue"):
        a.SetPDBResidueInfo(m)
    a.SetPDBResidueInfo(None)
    assert a.GetMonomerInfo() is None


def test_pickle_roundtrip():
    m = Chem.MolFromSmiles("CC")
    pdb = Chem.AtomPDBResidueInfo(
        " CA ", 12, "B", "ALA", 7, "A", "x", 0.5, 12.25, True, 2, 3, "cls"
    )
    m.GetAtomWithIdx(0).SetMonomerInfo(pdb)
    m.GetAtomWithIdx(1).SetMonomerInfo(
        Chem.AtomMonomerInfo(Chem.AtomMonomerType.OTHER, "X", "RES", 3, "C", "LGRP")
    )
    for m2 in (pickle.loads(pickle.dumps(m)), Chem.Mol(m.ToBinary()), Chem.Mol(m)):
        p = m2.GetAtomWithIdx(0).GetPDBResidueInfo()
        assert (p.GetName(), p.GetSerialNumber(), p.GetAltLoc(), p.GetResidueName()) == (
            " CA ", 12, "B", "ALA")
        assert (p.GetResidueNumber(), p.GetChainId(), p.GetInsertionCode()) == (7, "A", "x")
        assert (p.GetOccupancy(), p.GetTempFactor(), p.GetIsHeteroAtom()) == (0.5, 12.25, True)
        assert (p.GetSecondaryStructure(), p.GetSegmentNumber(), p.GetMonomerClass()) == (2, 3, "cls")
        o = m2.GetAtomWithIdx(1).GetMonomerInfo()
        assert (o.GetName(), o.GetChainId(), o.GetMonomerClass()) == ("X", "C", "LGRP")
