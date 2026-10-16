import numpy
import pytest
from rdoxide import DataStructs
from rdoxide.ML.InfoTheory import rdInfoTheory as rdit

BITS = ["0001", "0101", "0010", "1110"]
ACTS = [0, 0, 1, 1]


def _vects():
    out = []
    for s in BITS:
        bv = DataStructs.ExplicitBitVect(4)
        for i, c in enumerate(s):
            if c == "1":
                bv.SetBit(i)
        out.append(bv)
    return out


def _rank(info_type=None, bias=None):
    args = (4, 2) if info_type is None else (4, 2, info_type)
    ranker = rdit.InfoBitRanker(*args)
    if bias is not None:
        ranker.SetBiasList(bias)
    for bv, act in zip(_vects(), ACTS):
        ranker.AccumulateVotes(bv, act)
    return ranker.GetTopN(3)


@pytest.mark.parametrize(
    "info_type, bias, expected",
    [
        (None, None, [[3, 1.0, 2, 0], [2, 1.0, 0, 2], [0, 0.311, 0, 1]]),
        (rdit.InfoType.BIASENTROPY, (1,), [[2, 1.0, 0, 2], [0, 0.311, 0, 1], [1, 0.0, 1, 1]]),
        (rdit.InfoType.CHISQUARE, None, [[3, 4.0, 2, 0], [2, 4.0, 0, 2], [0, 1.333, 0, 1]]),
        (rdit.InfoType.BIASCHISQUARE, (1,), [[2, 4.0, 0, 2], [0, 1.333, 0, 1], [1, 0.0, 1, 1]]),
    ],
)
def test_ranker_docstring_examples(info_type, bias, expected):
    res = _rank(info_type, bias)
    assert res.shape == (3, 4)
    assert numpy.allclose(res, expected, atol=5e-4)


def test_info_type_is_exported():
    assert rdit.ENTROPY == rdit.InfoType.ENTROPY == 1
    assert rdit.InfoType.BIASCHISQUARE == 4


def test_ranker_errors():
    ranker = rdit.InfoBitRanker(4, 2)
    with pytest.raises(ValueError):
        ranker.AccumulateVotes("0101", 0)
    with pytest.raises(ValueError):
        ranker.GetTopN(5)


def test_write_top_bits(tmp_path):
    ranker = rdit.InfoBitRanker(4, 2)
    for bv, act in zip(_vects(), ACTS):
        ranker.AccumulateVotes(bv, act)
    ranker.GetTopN(2)
    out = tmp_path / "bits.txt"
    ranker.WriteTopBitsToFile(str(out))
    lines = out.read_text().splitlines()
    assert lines[0].split() == ["Bit", "InfoContent", "class0", "class1"]
    assert lines[1].split() == ["3", "1", "2", "0"]


def test_corr_mat_generator():
    gen = rdit.BitCorrMatGenerator()
    gen.SetBitList([0, 1, 2])
    for bv in _vects():
        gen.CollectVotes(bv)
    assert list(gen.GetCorrMatrix()) == [1.0, 1.0, 1.0]
    with pytest.raises(ValueError):
        gen.CollectVotes("nope")
