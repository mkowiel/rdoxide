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


def test_rdkit_fp_generator():
    m = Chem.MolFromSmiles("c1ccccc1CO")
    g = G.GetRDKitFPGenerator(maxPath=3, fpSize=512)
    assert list(g.GetFingerprint(m).GetOnBits()) == [
        2, 9, 36, 50, 76, 161, 232, 248, 277, 280, 294, 315, 330, 339, 342, 373, 407, 418, 465]
    assert g.GetOptions().maxPath == 3 and g.GetOptions().numBitsPerFeature == 2


def test_legacy_rdk_fingerprint():
    m = Chem.MolFromSmiles("CCO")
    ab, bi = [], {}
    fp = Chem.RDKFingerprint(m, maxPath=2, fpSize=64, atomBits=ab, bitInfo=bi)
    assert list(fp.GetOnBits()) == [0, 28, 31, 44, 50, 59]
    assert ab == [[28, 44, 0, 31], [28, 44, 59, 50, 0, 31], [59, 50, 0, 31]]
    assert bi == {0: [[0, 1]], 28: [[0]], 31: [[0, 1]], 44: [[0]], 50: [[1]], 59: [[1]]}
    assert len(Chem.UnfoldedRDKFingerprintCountBased(m).GetNonzeroElements()) == 3
