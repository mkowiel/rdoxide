"""Chem.LayeredFingerprint (expectations from the reference RDKit)."""

import pytest
from rdoxide import Chem
from rdoxide.DataStructs import ExplicitBitVect


def test_ethanol_bits_and_atom_counts():
    counts = [0, 0, 0]
    fp = Chem.LayeredFingerprint(Chem.MolFromSmiles("CCO"), fpSize=256, atomCounts=counts)
    assert list(fp.GetOnBits()) == [20, 84, 87, 92, 98, 99, 104, 162, 247, 248]
    assert counts == [2, 3, 2]


def test_substruct_layers():
    assert Chem.LayeredFingerprint_substructLayers == 0x07
    assert Chem._LayeredFingerprint_version == "0.7.0"
    fp = Chem.LayeredFingerprint(
        Chem.MolFromSmiles("C1CC1"), layerFlags=Chem.LayeredFingerprint_substructLayers, fpSize=128
    )
    assert list(fp.GetOnBits()) == [20, 34, 54, 97, 99, 104, 112, 118, 119]


def test_substructure_screen():
    sub = Chem.LayeredFingerprint(Chem.MolFromSmarts("c1ccccc1C(=O)O"), layerFlags=0x07)
    sup = Chem.LayeredFingerprint(Chem.MolFromSmiles("CCc1ccccc1C(=O)O"), layerFlags=0x07)
    assert (sub & sup) == sub


def test_set_only_bits_restricts_the_result():
    m = Chem.MolFromSmiles("CCO")
    full = Chem.LayeredFingerprint(m, fpSize=256)
    mask = ExplicitBitVect(256)
    mask.SetBit(20)
    only = Chem.LayeredFingerprint(m, fpSize=256, setOnlyBits=mask)
    assert list(only.GetOnBits()) == [20]
    assert (full & only) == only


def test_bad_arguments():
    m = Chem.MolFromSmiles("CCO")
    with pytest.raises(ValueError):
        Chem.LayeredFingerprint(m, atomCounts=[0])
    with pytest.raises(RuntimeError):
        Chem.LayeredFingerprint(m, minPath=0)
    with pytest.raises(RuntimeError):
        Chem.LayeredFingerprint(m, fpSize=64, setOnlyBits=ExplicitBitVect(32))
