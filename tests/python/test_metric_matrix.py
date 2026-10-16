"""rdMetricMatrixCalc (expectations from the reference RDKit)."""

import numpy as np
import pytest
from rdoxide import DataStructs
from rdoxide.DataManip.Metric import rdMetricMatrixCalc as mmc


def triangle(dtype):
    return np.array([[0, 0], [1, 0], [1, 1]], dtype)


@pytest.mark.parametrize("dtype", ["d", "f", "i"])
def test_euclidean_array(dtype):
    d = mmc.GetEuclideanDistMat(triangle(dtype))
    assert d.tolist() == pytest.approx([1.0, 2**0.5, 1.0])
    assert d.dtype == np.float64


def test_euclidean_sequences_and_memory_order():
    expected = mmc.GetEuclideanDistMat(triangle("d")).tolist()
    assert mmc.GetEuclideanDistMat([[0, 0], (1, 0), [1, 1]]).tolist() == expected
    # Fortran-ordered input is read by value
    assert mmc.GetEuclideanDistMat(np.asfortranarray(triangle("d"))).tolist() == expected


def test_large_matrix_uses_threads_without_changing_the_values():
    rs = np.random.RandomState(0)
    table = rs.rand(700, 8)
    d = mmc.GetEuclideanDistMat(table)
    assert d.shape == (700 * 699 // 2,)
    assert d[699 * 698 // 2 + 5] == pytest.approx(np.linalg.norm(table[699] - table[5]))


def test_tanimoto_matrices():
    a, b, c = (DataStructs.CreateFromBitString(s) for s in ("1100", "0110", "0001"))
    assert mmc.GetTanimotoSimMat([a, b, c]).tolist() == pytest.approx([1 / 3, 0.0, 0.0])
    assert mmc.GetTanimotoDistMat((a, b, c)).tolist() == pytest.approx([2 / 3, 1.0, 1.0])
    sparse = []
    for bits in ((0, 1), (1, 2)):
        bv = DataStructs.SparseBitVect(8)
        for bit in bits:
            bv.SetBit(bit)
        sparse.append(bv)
    assert mmc.GetTanimotoSimMat(sparse).tolist() == pytest.approx([1 / 3])


def test_errors():
    with pytest.raises(ValueError, match="same length"):
        mmc.GetEuclideanDistMat([[1, 2], [3]])
    with pytest.raises(ValueError, match="int, float, or double"):
        mmc.GetEuclideanDistMat(np.zeros((3, 2), bool))
    with pytest.raises(RuntimeError):
        mmc.GetEuclideanDistMat([])
    with pytest.raises(RuntimeError):
        mmc.GetTanimotoDistMat([DataStructs.ExplicitBitVect(8)])
    with pytest.raises(ValueError):
        mmc.GetTanimotoDistMat([1, 2, 3])
