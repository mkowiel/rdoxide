import numpy as np
import pytest

from rdoxide import Chem
from rdoxide.Chem import AllChem, rdDistGeom


def _mol(smi):
    return Chem.AddHs(Chem.MolFromSmiles(smi))


def test_embed_molecule_legacy_args():
    m = _mol("CCO")
    assert rdDistGeom.EmbedMolecule(m, randomSeed=42) == 0
    assert m.GetNumConformers() == 1


def test_embed_with_params_is_seeded():
    m1, m2 = _mol("C[C@H](N)C(=O)O"), _mol("C[C@H](N)C(=O)O")
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = 7
    assert rdDistGeom.EmbedMolecule(m1, ps) == 0
    assert rdDistGeom.EmbedMolecule(m2, ps) == 0
    assert np.allclose(m1.GetConformer().GetPositions(), m2.GetConformer().GetPositions())


def test_multiple_confs_thread_independent():
    res = []
    for nt in (1, 3):
        m = _mol("CCCCO")
        ps = rdDistGeom.ETKDGv2()
        ps.randomSeed = 11
        ps.numThreads = nt
        ids = rdDistGeom.EmbedMultipleConfs(m, 4, ps)
        assert list(ids) == [0, 1, 2, 3]
        res.append([c.GetPositions() for c in m.GetConformers()])
    for a, b in zip(*res):
        assert np.allclose(a, b)


def test_params_attributes():
    ps = rdDistGeom.EmbedParameters()
    ps.maxIterations = 5
    assert ps.maxIterations == 5
    with pytest.raises(AttributeError):
        ps.notAnAttribute = 1


def test_bounds_matrix():
    m = _mol("CC")
    b = rdDistGeom.GetMoleculeBoundsMatrix(m)
    assert b.shape == (m.GetNumAtoms(),) * 2
    assert b[1, 0] <= b[0, 1]


def test_experimental_torsions():
    m = Chem.MolFromSmiles("CCCCO")
    t = rdDistGeom.GetExperimentalTorsions(m)
    assert len(t) > 0 and "atomIndices" in t[0]


def test_distance_geometry_module():
    from rdoxide.DistanceGeometry import DoTriangleSmoothing, EmbedBoundsMatrix

    m = _mol("CCO")
    b = rdDistGeom.GetMoleculeBoundsMatrix(m, doTriangleSmoothing=False)
    assert DoTriangleSmoothing(b)
    pos = EmbedBoundsMatrix(b, randomSeed=3)
    assert pos.shape == (m.GetNumAtoms(), 3)
