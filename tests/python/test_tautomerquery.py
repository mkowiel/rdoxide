import pytest

from rdoxide import Chem
from rdoxide.Chem import rdTautomerQuery


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


def _library(holder_cls):
    from rdoxide.Chem import rdSubstructLibrary
    mols = rdSubstructLibrary.MolHolder()
    lib = rdSubstructLibrary.SubstructLibrary(mols, holder_cls())
    for smi in ["O=C1CCCC(CC)C1", "OC1=CCCC(CC)C1", "CCCCC", "c1ccccc1"]:
        lib.AddMol(Chem.MolFromSmiles(smi))
    return lib


def test_substruct_library_tautomer_search():
    from rdoxide.Chem import rdSubstructLibrary
    query = rdTautomerQuery.TautomerQuery(Chem.MolFromSmiles("O=C1CCCCC1"))
    for holder in (rdSubstructLibrary.PatternHolder, rdSubstructLibrary.TautomerPatternHolder):
        lib = _library(holder)
        assert sorted(lib.GetMatches(query)) == [0, 1]
        assert lib.CountMatches(query) == 2
        assert lib.HasMatch(query)
        assert not lib.HasMatch(rdTautomerQuery.TautomerQuery(Chem.MolFromSmiles("CC(=O)N")))
    # a plain query only finds the keto form
    plain = _library(rdSubstructLibrary.PatternHolder).GetMatches(Chem.MolFromSmiles("O=C1CCCCC1"))
    assert list(plain) == [0]
    with pytest.raises(TypeError):
        _library(rdSubstructLibrary.PatternHolder).GetMatches("O=C")
