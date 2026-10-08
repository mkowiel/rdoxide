"""Batch fingerprint calls with numThreads (results must not depend on the thread count)."""

import pytest

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as G

SMILES = [
    "CC(=O)Nc1ccc(O)cc1",
    "c1ccccc1",
    "CCO",
    "C1CCCCC1C(=O)O",
    "Oc1ccccc1/C=N/NC(=O)c1ccncc1",
    "CN1CCC[C@H]1c1cccnc1",
    "OC(=O)c1ccncc1C(=O)O",
    "CCBr",
]


def _mols(n=200):
    return [Chem.MolFromSmiles(SMILES[i % len(SMILES)]) for i in range(n)]


def _bits(fps):
    return [tuple(fp.GetOnBits()) for fp in fps]


def _counts(fps):
    return [sorted(fp.GetNonzeroElements().items()) for fp in fps]


@pytest.mark.parametrize("threads", [1, 2, 4, 0, -1])
def test_get_fingerprints_is_independent_of_thread_count(threads):
    gen = G.GetMorganGenerator(radius=2)
    mols = _mols()
    expected = _bits([gen.GetFingerprint(m) for m in mols])
    assert _bits(gen.GetFingerprints(mols, numThreads=threads)) == expected


@pytest.mark.parametrize("threads", [1, 3])
def test_other_batch_flavours(threads):
    gen = G.GetRDKitFPGenerator(maxPath=4, fpSize=512)
    mols = _mols(60)
    assert _counts(gen.GetCountFingerprints(mols, numThreads=threads)) == _counts(
        [gen.GetCountFingerprint(m) for m in mols])
    assert _counts(gen.GetSparseCountFingerprints(mols, numThreads=threads)) == _counts(
        [gen.GetSparseCountFingerprint(m) for m in mols])
    assert [tuple(fp.GetOnBits()) for fp in gen.GetSparseFingerprints(mols, numThreads=threads)
            ] == [tuple(gen.GetSparseFingerprint(m).GetOnBits()) for m in mols]


def test_none_entries_stay_in_place():
    gen = G.GetMorganGenerator(radius=1)
    mols = _mols(10)
    mols[3] = None
    mols[9] = None
    out = gen.GetFingerprints(mols, numThreads=4)
    assert len(out) == 10
    assert out[3] is None and out[9] is None
    assert all(fp is not None for i, fp in enumerate(out) if i not in (3, 9))


def test_empty_input():
    assert len(G.GetMorganGenerator().GetFingerprints([], numThreads=4)) == 0
