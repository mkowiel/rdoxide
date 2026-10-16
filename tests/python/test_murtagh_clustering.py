"""rdoxide.ML.Cluster.Clustering (expectations from the reference RDKit)."""

import numpy as np
import pytest
from rdoxide.ML.Cluster import Clustering


def test_points_single_link():
    x = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 2.0], [5.0, 5.0], [6.0, 5.0]])
    ia, ib, crit = Clustering.MurtaghCluster(x, 5, 2, 2)
    assert ia.tolist() == [1, 4, 1, 1, 0] and ib.tolist() == [2, 5, 3, 4, 0]
    assert crit.tolist() == [1.0, 1.0, 4.0, 34.0, 0.0]
    assert ia.dtype == np.dtype("l") and crit.dtype == np.float64


def test_distances_update_the_callers_array_in_place():
    d = np.array([1.0, 4.0, 3.0, 9.0, 8.0, 2.0, 7.0, 6.0, 5.0, 10.0])
    ia, ib, crit = Clustering.MurtaghDistCluster(d, 5, 1)  # Ward
    assert ia.tolist() == [1, 3, 1, 1, 0] and ib.tolist() == [2, 4, 5, 3, 0]
    assert crit.tolist()[:2] == [0.5, 1.0]
    assert d[0] == 0.5  # upstream's algorithm works in place


def test_argument_errors():
    with pytest.raises(ValueError):
        Clustering.MurtaghCluster([[0.0, 1.0]], 1, 2, 1)  # not an array
    with pytest.raises(ValueError):
        Clustering.MurtaghDistCluster(np.zeros(3), 5, 1)  # too small
