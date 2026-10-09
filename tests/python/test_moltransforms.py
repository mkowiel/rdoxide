"""rdMolTransforms (expectations from RDKit 2026.09.1, on coordinates read from tests/data/mols3d.json)."""

import json
import math
import os

import numpy as np
import pytest

from rdkit import Chem
from rdkit.Chem import rdMolTransforms as T
from rdkit.Geometry import Point3D

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)


def load(name):
    return Chem.MolFromMolBlock(BLOCKS[name], removeHs=False)


def test_principal_axes_and_moments():
    conf = load("ethanol").GetConformer()
    axes, moments = T.ComputePrincipalAxesAndMoments(conf)
    assert moments == pytest.approx([0.4960649635983768, 2.7251972097349584, 3.2212621733333333], rel=1e-12)
    assert axes[0] == pytest.approx([0.9456520915724027, 0.31544733007111064, -0.07896267254689608], abs=1e-12)
    axes, moments = T.ComputePrincipalAxesAndMomentsFromGyrationMatrix(conf)
    assert moments == pytest.approx([0.0, 0.16535498786612554, 0.9083990699116522], abs=1e-12)
    # weights change the result and must match the number of atoms
    _, weighted = T.ComputePrincipalAxesAndMoments(conf, True, [1.0 + 0.5 * i for i in range(9)])
    assert weighted == pytest.approx([0.7234171409863588, 3.737469789013641, 4.46088693], rel=1e-12)
    with pytest.raises(ValueError, match="length equal to conf.GetNumAtoms"):
        T.ComputePrincipalAxesAndMoments(conf, True, [1.0])


def test_canonical_transform_and_canonicalization():
    mol = load("ethanol")
    conf = mol.GetConformer()
    trans = T.ComputeCanonicalTransform(conf)
    assert trans.shape == (4, 4)
    assert trans[0] == pytest.approx([0.9456520915724027, 0.12673756086503435, 0.2994657115108879, -0.3043494585647698])
    assert trans[3].tolist() == [0.0, 0.0, 0.0, 1.0]
    T.CanonicalizeMol(mol)
    # the centroid of the heavy atoms sits at the origin afterwards
    heavy = mol.GetConformer().GetPositions()[:3]
    assert np.allclose(heavy.mean(axis=0), 0.0, atol=1e-12)


def test_single_atom_and_centroid():
    mol = load("single")
    conf = mol.GetConformer()
    assert np.array_equal(T.ComputeCanonicalTransform(conf), np.eye(4))
    c = T.ComputeCentroid(conf)
    assert (c.x, c.y, c.z) == (0.0, 0.0, 0.0)
    with pytest.raises(TypeError):
        T.ComputeCentroid(conf, True, [1.0])


def test_transform_conformer():
    conf = load("ethanol").GetConformer()
    before = conf.GetPositions()
    shift = np.eye(4)
    shift[:3, 3] = [1.0, -2.0, 3.0]
    T.TransformConformer(conf, shift)
    assert np.allclose(conf.GetPositions(), before + [1.0, -2.0, 3.0])
    with pytest.raises(ValueError, match="numeric array"):
        T.TransformConformer(conf, [[1]])


def test_internal_coordinates():
    conf = load("ethanol").GetConformer()
    assert T.GetAngleDeg(conf, 0, 1, 2) == pytest.approx(106.73896473395976)
    T.SetAngleDeg(conf, 0, 1, 2, 100.0)
    assert T.GetAngleDeg(conf, 0, 1, 2) == pytest.approx(100.0)
    T.SetDihedralDeg(conf, 0, 1, 2, 8, 60.0)
    assert T.GetDihedralDeg(conf, 0, 1, 2, 8) == pytest.approx(60.0)
    assert T.GetDihedralRad(conf, 0, 1, 2, 8) == pytest.approx(math.pi / 3)
    T.SetBondLength(conf, 0, 1, 2.0)
    assert T.GetBondLength(conf, 0, 1) == pytest.approx(2.0)


def test_errors():
    ring = load("cyclohexane").GetConformer()
    with pytest.raises(ValueError, match=r"bond \(i,j\) must not belong to a ring"):
        T.SetBondLength(ring, 0, 1, 2.0)
    with pytest.raises(ValueError, match="must not both belong to a ring"):
        T.SetAngleDeg(ring, 0, 1, 2, 100)
    with pytest.raises(ValueError, match=r"bond \(j,k\) must not belong to a ring"):
        T.SetDihedralDeg(ring, 0, 1, 2, 3, 100)
    conf = load("ethanol").GetConformer()
    with pytest.raises(ValueError, match="must be bonded"):
        T.SetBondLength(conf, 0, 2, 1.0)
    with pytest.raises(RuntimeError, match="jAtomId"):
        T.GetBondLength(conf, 0, 99)
    with pytest.raises(RuntimeError, match="no owner"):
        T.ComputeCentroid(Chem.Conformer(3))
    assert T.GetBondLength(Chem.Conformer(3), 0, 1) == 0.0
