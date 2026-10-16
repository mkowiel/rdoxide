"""Chem.PatternFingerprint (expectations from the reference RDKit)."""

import pytest
from rdoxide import Chem
from rdoxide.DataStructs import ExplicitBitVect


def test_ethanol_bits():
    fp = Chem.PatternFingerprint(Chem.MolFromSmiles("CCO"), 512)
    assert list(fp.GetOnBits()) == [9, 37, 54, 110, 173, 212, 217, 247, 265, 266, 271, 362, 429, 434, 455, 510]


def test_tautomer_fingerprints_differ():
    m = Chem.MolFromSmiles("O=C1NC=CC=C1")
    plain = Chem.PatternFingerprint(m, 512)
    taut = Chem.PatternFingerprint(m, 512, tautomerFingerprints=True)
    assert plain != taut
    assert taut.GetNumBits() == 512


def test_substructure_screen():
    # a pattern fingerprint of a substructure is a subset of its superstructure's
    sub = Chem.PatternFingerprint(Chem.MolFromSmarts("c1ccccc1C(=O)O"))
    sup = Chem.PatternFingerprint(Chem.MolFromSmiles("CCc1ccccc1C(=O)O"))
    assert (sub & sup) == sub


def test_arguments():
    m = Chem.MolFromSmiles("CCO")
    assert Chem._PatternFingerprint_version == "1.0.0"
    with pytest.raises(ValueError):
        Chem.PatternFingerprint(m, atomCounts=[1])
    with pytest.raises(RuntimeError):
        Chem.PatternFingerprint(m, 512, setOnlyBits=ExplicitBitVect(64))
    assert Chem.PatternFingerprint(m, 512, atomCounts=[0, 0, 0]) == Chem.PatternFingerprint(m, 512)
