"""Subgraph / path enumeration and rdBase vector types (expectations from RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem, rdBase


def test_paths_and_subgraphs():
    m = Chem.MolFromSmiles("CC(C)CO")
    assert [list(p) for p in Chem.FindAllPathsOfLengthN(m, 3, useBonds=False)] == [
        [0, 1, 2], [0, 1, 3], [1, 3, 4], [2, 1, 3]]
    assert Chem.FindAllSubgraphsOfLengthMToN(m, 1, 2) == (
        [(0,), (1,), (2,), (3,)], [(0, 2), (0, 1), (1, 2), (2, 3)])
    assert len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 4
    assert len(Chem.FindUniqueSubgraphsOfLengthN(m, 2)) == 2
    d = {}
    env = Chem.FindAtomEnvironmentOfRadiusN(m, 2, 0, atomMap=d)
    assert list(env) == [0, 1, 2] and d == {0: 0, 1: 1, 2: 2, 3: 2}


def test_vect_types():
    m = Chem.MolFromSmiles("CC(C)CO")
    v = Chem.FindAtomEnvironmentOfRadiusN(m, 2, 0)
    assert type(v).__name__ == "_vecti" and isinstance(v, rdBase._vecti)
    assert v != [0, 1, 2]  # no value equality, like Boost.Python's indexing suite
    assert v[-1] == 2 and 1 in v and list(v[0:2]) == [0, 1]
    v.append(9)
    v[0] = 5
    del v[1]
    assert list(v) == [5, 2, 9]
    with pytest.raises(IndexError):
        v[10]
    paths = Chem.FindAllPathsOfLengthN(m, 2)
    assert type(paths[0]).__name__ == "_vecti"
