"""Legacy fingerprint functions of rdMolDescriptors (expectations from RDKit 2026.09.1)."""

from rdoxide import Chem, RDLogger
from rdoxide.Chem import rdMolDescriptors as D

RDLogger.DisableLog("rdApp.*")


def test_legacy_morgan():
    m = Chem.MolFromSmiles("CCO")
    bi = {}
    fp = D.GetMorganFingerprint(m, 1, bitInfo=bi)
    assert sorted(fp.GetNonzeroElements().items()) == [
        (864662311, 1), (1535166686, 1), (2245384272, 1), (2246728737, 1), (3542456614, 1), (4018048386, 1)]
    assert bi[2246728737] == ((0, 0),) and bi[3542456614] == ((0, 1),)
    assert D.GetFeatureInvariants(Chem.MolFromSmiles("OC(=O)c1ccccc1N")) == [1, 32, 2, 4, 4, 4, 4, 4, 4, 19]


def test_legacy_atom_pairs_and_torsions():
    m = Chem.MolFromSmiles("CCO")
    ap = D.GetAtomPairFingerprint(m)
    assert type(ap).__name__ == "IntSparseIntVect"
    assert sorted(ap.GetNonzeroElements().items()) == [(558113, 1), (1590306, 1), (1590337, 1)]
    assert type(D.GetTopologicalTorsionFingerprint(m)).__name__ == "LongSparseIntVect"
