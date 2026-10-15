import pytest

from rdkit import Chem
from rdkit.Chem import rdTautomerQuery


def test_query_properties():
    q = rdTautomerQuery.TautomerQuery(Chem.MolFromSmiles("O=C1CCCCC1"))
    assert len(q.GetTautomers()) == 2
    assert len(q.GetModifiedAtoms()) == 3
    assert len(q.GetModifiedBonds()) == 3
    assert q.GetTemplateMolecule().GetNumAtoms() == 7


def test_searches_and_params():
    q = rdTautomerQuery.TautomerQuery(Chem.MolFromSmiles("O=C1CCCCC1"))
    target = Chem.MolFromSmiles("O=C1CCCC(CC)C1")
    assert q.IsSubstructOf(target)
    assert not q.IsSubstructOf(Chem.MolFromSmiles("CCCC"))
    assert len(q.GetSubstructMatch(target)) == 7
    assert len(q.GetSubstructMatches(target)) == 1
    assert len(q.GetSubstructMatches(target, uniquify=False)) == 2
    params = Chem.SubstructMatchParameters()
    params.uniquify = False
    pairs = q.GetSubstructMatchesWithTautomers(target, params)
    assert len(pairs) == 2
    assert all(isinstance(t, Chem.Mol) for _, t in pairs)


def test_fingerprints():
    q = rdTautomerQuery.TautomerQuery(Chem.MolFromSmiles("O=C1CCCCC1"))
    fp = q.PatternFingerprintTemplate(1024)
    assert fp.GetNumBits() == 1024
    assert rdTautomerQuery.PatternFingerprintTautomerTarget(Chem.MolFromSmiles("OC1=CCCCC1")).GetNumBits() == 2048


def test_serialization_is_unavailable():
    q = rdTautomerQuery.TautomerQuery(Chem.MolFromSmiles("CC=O"))
    assert not rdTautomerQuery.TautomerQueryCanSerialize()
    with pytest.raises(RuntimeError):
        q.ToBinary()
    with pytest.raises(RuntimeError):
        rdTautomerQuery.TautomerQuery("pickle")
