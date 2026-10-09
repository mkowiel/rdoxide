"""Uniform grids, PointND and dihedral helpers (expectations from RDKit 2026.09.1)."""

import math
import pickle

import pytest

from rdkit.Geometry import rdGeometry as geo


def sphere_grid():
    grid = geo.UniformGrid3D(10, 10, 10, 0.5)
    grid.SetSphereOccupancy(geo.Point3D(1, 1, 1), 1.5, 0.25)
    return grid


def test_grid_layout():
    grid = sphere_grid()
    assert (grid.GetNumX(), grid.GetNumY(), grid.GetNumZ(), grid.GetSize()) == (20, 20, 20, 8000)
    assert grid.GetOffset().x == -5.0
    assert grid.GetSpacing() == 0.5
    assert grid.GetGridIndices(grid.GetGridIndex(1, 2, 3)) == (1, 2, 3)


def test_grid_pickle_and_distances():
    grid = sphere_grid()
    other = pickle.loads(pickle.dumps(grid))
    assert grid.CompareParams(other)
    assert geo.TanimotoDistance(grid, other) == 0.0
    assert geo.ComputeGridCentroid(grid, geo.Point3D(1, 1, 1), 1.0)[0] == 279.0
    assert len(geo.FindGridTerminalPoints(grid, 1.0, 0.5)) == 12


def test_real_value_grid():
    grid = geo.UniformRealValueGrid3D(5, 5, 5, 0.5)
    grid.SetVal(3, 2.5)
    copy = pickle.loads(pickle.dumps(grid))
    assert copy.GetVal(3) == 2.5
    assert grid.CompareGrids(copy)
    assert geo.UniformRealValueGrid3D().GetSize() == 1
    assert geo.UniformRealValueGrid3D(grid).CompareGrids(grid)


def test_point_nd():
    p = geo.PointND(3)
    p[0] = 1
    p[-1] = 2
    assert len(p) == 3
    assert p.Length() == pytest.approx(math.sqrt(5))
    assert pickle.loads(pickle.dumps(p))[2] == 2.0
    with pytest.raises(IndexError):
        p[3]
    with pytest.raises(TypeError):
        p.Distance(p)


def test_point2d_angles_and_dihedral():
    a, b = geo.Point2D(1, 0), geo.Point2D(0, 1)
    assert a.AngleTo(b) == pytest.approx(math.pi / 2)
    assert a.SignedAngleTo(b) == pytest.approx(math.pi / 2)
    assert a.DirectionVector(b).x == pytest.approx(-math.sqrt(0.5))
    pts = [geo.Point3D(*c) for c in ((1, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, 1))]
    assert geo.ComputeDihedralAngle(*pts) == pytest.approx(math.pi / 2)
