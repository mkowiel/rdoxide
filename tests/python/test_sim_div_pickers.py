"""rdSimDivPickers MaxMin and Leader pickers (expectations from the reference RDKit)."""

import numpy as np
import pytest
from rdoxide import DataStructs
from rdoxide.SimDivFilters import rdSimDivPickers as pickers

POINTS = np.array([0.0, 1.0, 2.5, 4.0, 4.2, 7.0, 9.0, 9.5, 10.0, 15.0])


def line(i, j):
    return abs(POINTS[i] - POINTS[j])


def test_max_min_distance_matrix():
    n = len(POINTS)
    dists = np.array([line(i, j) for i in range(n) for j in range(i)])
    picker = pickers.MaxMinPicker()
    assert list(picker.Pick(dists, n, 4, [0], seed=1)) == [0, 9, 5, 3]
    assert list(picker.Pick(dists, n, 4, seed=7)) == [0, 9, 5, 3]


def test_max_min_lazy_with_threshold():
    picks, last = pickers.MaxMinPicker().LazyPickWithThreshold(line, 10, 10, 3.0, [0], seed=1)
    assert list(picks) == [0, 9, 5] and last == 7.0


def test_max_min_bit_vectors():
    bvs = []
    for i in range(20):
        bv = DataStructs.ExplicitBitVect(32)
        for b in range(i % 7 + 1):
            bv.SetBit((i * 3 + b * 5) % 32)
        bvs.append(bv)
    a = pickers.MaxMinPicker().LazyBitVectorPick(bvs, 20, 5, seed=3)
    assert len(set(a)) == 5
    assert list(a) == list(pickers.MaxMinPicker().LazyBitVectorPick(bvs, 20, 5, seed=3))


def test_leader():
    leader = pickers.LeaderPicker()
    assert list(leader.LazyPick(line, 10, 1.5)) == [0, 2, 4, 5, 6, 9]
    assert list(leader.LazyPick(line, 10, 1.5, 2, [9])) == [9, 0]


def test_errors():
    picker = pickers.MaxMinPicker()
    with pytest.raises(ValueError):
        picker.LazyPick(line, 3, 5)
    with pytest.raises(ValueError):
        picker.Pick([0.1, 0.2, 0.3], 3, 2)  # not a numpy array
    with pytest.raises(ZeroDivisionError):
        picker.LazyPick(lambda i, j: 1 / 0, 5, 2)


def test_hierarchical_picker():
    pts = np.array([0.0, 1.0, 3.0, 7.0, 8.0, 15.0, 16.5])
    dists = np.array([abs(pts[i] - pts[j]) for i in range(7) for j in range(i)])
    ward = pickers.HierarchicalClusterPicker(pickers.ClusterMethod.WARD)
    assert [list(c) for c in ward.Cluster(dists, 7, 3)] == [[0, 1, 2], [3, 4], [5, 6]]
    assert list(ward.Pick(dists, 7, 3)) == [1, 3, 5]
    slink = pickers.HierarchicalClusterPicker(pickers.SLINK)
    assert list(slink.Pick(dists, 7, 3)) == [0, 3, 5]
    assert int(pickers.CENTROID) == 7
    before = dists.copy()
    ward.Pick(dists, 7, 3)
    assert (dists == before).all()  # the picker works on a copy
    with pytest.raises(ValueError):
        ward.Pick(dists, 7, 7)
