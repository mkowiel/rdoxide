"""Distance and adjacency matrices (expectations from RDKit 2026.09.1)."""

import numpy as np
import pytest

from rdoxide import Chem
from rdoxide.Chem import rdDepictor


def test_topological_distance_matrix():
  m = Chem.MolFromSmiles("CC(=O)O.[Na+]")
  d = Chem.GetDistanceMatrix(m)
  assert d.dtype == np.float64 and d.shape == (5, 5)
  assert d[0].tolist() == [0.0, 1.0, 2.0, 2.0, 1e8]
  assert d[4].tolist() == [1e8, 1e8, 1e8, 1e8, 0.0]


def test_distance_matrix_with_bond_orders_and_atom_weights():
  m = Chem.MolFromSmiles("C=CC#N")
  d = Chem.GetDistanceMatrix(m, useBO=True)
  assert d[0].tolist() == [0.0, 0.5, 1.5, 1.8333333333333333]
  w = Chem.GetDistanceMatrix(m, useAtomWts=True, force=True)
  assert np.diag(w).tolist() == [1.0, 1.0, 1.0, 6.0 / 7.0]
  benzene = Chem.GetDistanceMatrix(Chem.MolFromSmiles("c1ccccc1"), useBO=True)
  assert benzene[0, 1] == pytest.approx(2.0 / 3.0)


def test_adjacency_matrix():
  m = Chem.MolFromSmiles("C=CC#N")
  a = Chem.GetAdjacencyMatrix(m)
  assert a.dtype == np.intc
  assert a.tolist() == [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]]
  bo = Chem.GetAdjacencyMatrix(m, useBO=True)
  assert bo.dtype == np.float64
  assert bo.tolist() == [[0, 2, 0, 0], [2, 0, 1, 0], [0, 1, 0, 3], [0, 0, 3, 0]]


def test_empty_molecule():
  m = Chem.MolFromSmiles("")
  assert Chem.GetDistanceMatrix(m).shape == (0, 0)
  assert Chem.GetAdjacencyMatrix(m).shape == (0, 0)


def test_matrices_are_cached_as_computed_properties():
  m = Chem.MolFromSmiles("CCO")
  Chem.GetDistanceMatrix(m)
  Chem.GetDistanceMatrix(m, useBO=True, prefix="x")
  names = set(m.GetPropNames(includePrivate=True, includeComputed=True))
  assert {"DistanceMatrix", "DistanceMatrix_Paths", "xDistanceMatrixBO"} <= names
  assert "DistanceMatrix" not in m.GetPropsAsDict(includePrivate=True, includeComputed=True)
  # the first call wins until `force` is given
  first = Chem.GetAdjacencyMatrix(m, useBO=True)
  assert Chem.GetAdjacencyMatrix(m, useBO=True, emptyVal=3).tolist() == first.tolist()
  assert Chem.GetAdjacencyMatrix(m, useBO=True, emptyVal=0, force=True).tolist() == first.tolist()
  m.ClearComputedProps()
  assert not m.HasProp("DistanceMatrix")


def test_cache_of_a_different_size_is_recomputed():
  m = Chem.RWMol(Chem.MolFromSmiles("CCO"))
  Chem.GetDistanceMatrix(m)
  m.AddAtom(Chem.Atom(6))
  assert Chem.GetDistanceMatrix(m).shape == (4, 4)


def test_3d_distance_matrix():
  m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
  rdDepictor.Compute2DCoords(m)
  d = Chem.Get3DDistanceMatrix(m)
  assert d.shape == (9, 9)
  assert d[0, 1] == pytest.approx(1.5)
  assert np.allclose(d, d.T) and np.all(np.diag(d) == 0)
  w = Chem.Get3DDistanceMatrix(m, useAtomWts=True)
  assert w[1, 1] == 1.0 and w[2, 2] == 0.75
  with pytest.raises(ValueError, match="Bad Conformer Id"):
    Chem.Get3DDistanceMatrix(Chem.MolFromSmiles("CC"))
