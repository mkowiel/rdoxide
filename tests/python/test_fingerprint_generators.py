"""rdFingerprintGenerator Morgan generator (expectations from RDKit 2026.09.1)."""

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as G


def test_morgan_basics():
    m = Chem.MolFromSmiles("c1ccccc1O")
    g = G.GetMorganGenerator(radius=2, fpSize=2048)
    fp = g.GetFingerprint(m)
    assert fp.GetNumBits() == 2048
    assert list(fp.GetOnBits()) == [65, 389, 745, 807, 1088, 1199, 1380, 1602, 1750, 1778, 1873]
    cfp = g.GetSparseCountFingerprint(m)
    assert sum(cfp.GetNonzeroElements().values()) == 20
    assert g.GetOptions().radius == 2
    g.GetOptions().radius = 1
    assert g.GetOptions().radius == 1
    assert "radius=1" in g.GetInfoString()


def test_morgan_additional_output():
    m = Chem.MolFromSmiles("CCO")
    g = G.GetMorganGenerator(radius=1, fpSize=64)
    ao = G.AdditionalOutput()
    ao.AllocateBitInfoMap()
    ao.AllocateAtomCounts()
    g.GetFingerprint(m, additionalOutput=ao)
    assert ao.GetAtomCounts() == (2, 2, 2)
    assert ao.GetAtomToBits() is None
    assert sum(len(v) for v in ao.GetBitInfoMap().values()) == 6


def test_morgan_bulk_and_numpy():
    mols = [Chem.MolFromSmiles(s) for s in ("C", "CC", "CCO")]
    g = G.GetMorganGenerator(radius=2)
    fps = g.GetFingerprints(mols)
    assert [fp.GetNumOnBits() for fp in fps] == [g.GetFingerprint(m).GetNumOnBits() for m in mols]
    arr = g.GetFingerprintAsNumPy(mols[2])
    assert arr.dtype.name == "uint8" and int(arr.sum()) == fps[2].GetNumOnBits()


def test_atom_pair_generator():
    m = Chem.MolFromSmiles("CCO")
    g = G.GetAtomPairGenerator(countSimulation=False)
    assert sorted(g.GetSparseCountFingerprint(m).GetNonzeroElements().items()) == [
        (558113, 1), (1590306, 1), (1590337, 1)]
    assert list(g.GetFingerprint(m).GetOnBits()) == [1375, 1423, 1503]
    opts = g.GetOptions()
    assert opts.maxDistance == 30
    opts.maxDistance = 1
    assert len(g.GetSparseCountFingerprint(m).GetNonzeroElements()) == 2


def test_topological_torsion_generator():
    m = Chem.MolFromSmiles("CCCCO")
    g = G.GetTopologicalTorsionGenerator(countSimulation=False)
    assert sorted(g.GetSparseCountFingerprint(m).GetNonzeroElements().items()) == [
        (4437590048, 1), (12893306913, 1)]
    assert g.GetOptions().torsionAtomCount == 4
