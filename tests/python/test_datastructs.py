"""rdkit.DataStructs bit vectors (expectations checked against RDKit 2026.09.1)."""

import pickle

import pytest

from rdkit import DataStructs as DS


def _bv(n, bits, cls=DS.ExplicitBitVect):
    bv = cls(n)
    bv.SetBitsFromList(bits)
    return bv


@pytest.mark.parametrize("cls", [DS.ExplicitBitVect, DS.SparseBitVect])
def test_basics(cls):
    bv = _bv(64, [0, 5, 63], cls)
    assert len(bv) == 64 and bv.GetNumOnBits() == 3 and bv.GetNumOffBits() == 61
    assert list(bv.GetOnBits()) == [0, 5, 63]
    assert bv[-1] == 1 and bv[1] == 0
    bv[1] = 1
    assert bv.GetBit(1)
    with pytest.raises(IndexError):
        bv[64]
    assert pickle.loads(pickle.dumps(bv)) == bv
    assert cls(bv.ToBinary()) == bv
    other = cls(64)
    other.FromBase64(bv.ToBase64())
    assert other == bv
    assert (~bv).GetNumOnBits() == 60


def test_similarities():
    a = _bv(128, [1, 3, 5, 7, 9, 100])
    b = _bv(128, [3, 5, 9, 11, 127])
    expected = {
        "Tanimoto": 0.375,
        "Dice": 0.5454545454545454,
        "Cosine": 0.5477225575051661,
        "Sokal": 0.23076923076923078,
        "Russel": 0.0234375,
        "RogotGoldberg": 0.7625231910946196,
        "AllBit": 0.9609375,
        "Kulczynski": 0.55,
        "McConnaughey": 0.1,
        "Asymmetric": 0.6,
        "BraunBlanquet": 0.5,
    }
    for name, fn, _ in DS.similarityFunctions:
        assert fn(a, b) == expected[name], name
    assert DS.TanimotoSimilarity(a, b, returnDistance=True) == 0.625
    assert DS.TanimotoSimilarity(a, b.ToBinary()) == 0.375
    assert DS.BulkTanimotoSimilarity(a, [a, b]) == [1.0, 0.375]
    assert DS.TverskySimilarity(a, b, 0.5, 0.5) == expected["Dice"]
    assert list(DS.OnBitsInCommon(a, b)) == [3, 5, 9]
    assert DS.NumBitsInCommon(a, b) == 123
    assert DS.TanimotoSimilarityNeighbors([a], [b, a]) == [(1, 1.0)]


def test_text_and_fold():
    a = _bv(16, [0, 9])
    assert DS.BitVectToText(a) == "1000000001000000"
    assert a.ToBitString() == "1000000001000000"
    assert DS.BitVectToFPSText(a) == "0102"
    assert DS.CreateFromFPSText("0102") == a
    assert DS.CreateFromBitString("1000000001000000") == a
    assert DS.CreateFromBinaryText(DS.BitVectToBinaryText(a)) == a
    assert list(DS.FoldFingerprint(a).GetOnBits()) == [0, 1]
    assert DS.FingerprintSimilarity(a, DS.FoldFingerprint(a)) == 1.0
    assert DS.ConvertToExplicit(_bv(16, [0, 9], DS.SparseBitVect)) == a
    assert (a + a).GetNumBits() == 32
