"""UFF force fields; the expectations come from RDKit's own testHelpers.py."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import ChemicalForceFields
from rdoxide.ForceField import rdForceField

BENZENE = """
     RDKit          3D

 12 12  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.3900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2038    0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2038   -0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.3900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2038   -0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2038    0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.4700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1391    1.2350    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1391   -1.2350    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.4700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1391   -1.2350    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1391    1.2350    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  2  8  1  0
  3  9  1  0
  4 10  1  0
  5 11  1  0
  6 12  1  0
M  END
"""


def benzene():
    return Chem.MolFromMolBlock(BENZENE, removeHs=False)


def test_optimize_converges():
    m = benzene()
    assert ChemicalForceFields.UFFOptimizeMolecule(m) == 0
    assert ChemicalForceFields.UFFOptimizeMolecule(benzene(), maxIters=1) == 1
    assert ChemicalForceFields.UFFOptimizeMolecule(benzene(), vdwThresh=2.0) == 0
    with pytest.raises(ValueError):
        ChemicalForceFields.UFFOptimizeMolecule(benzene(), confId=1)


def test_optimize_moves_the_atoms():
    m = benzene()
    before = m.GetConformer().GetPositions().copy()
    ChemicalForceFields.UFFOptimizeMolecule(m)
    assert abs(m.GetConformer().GetPositions() - before).max() > 1e-3
    # a C-C aromatic bond ends up around 1.4 Angstrom
    pos = m.GetConformer().GetPositions()
    assert 1.3 < ((pos[0] - pos[1]) ** 2).sum() ** 0.5 < 1.5


def test_force_field_is_linked_to_the_molecule():
    m = benzene()
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    assert ff.NumPoints() == 12 and ff.Dimension() == 3
    e1 = ff.CalcEnergy()
    assert ff.Minimize() == 0
    e2 = ff.CalcEnergy()
    assert e2 < e1
    # the minimization moved the atoms of the molecule
    assert tuple(m.GetConformer().GetPositions().flatten()) == pytest.approx(ff.Positions())
    # and the force field sees coordinate changes of the molecule
    conf = m.GetConformer()
    p = conf.GetAtomPosition(0)
    conf.SetAtomPosition(0, (p.x + 0.5, p.y, p.z))
    assert ff.CalcEnergy() > e2
    assert ff.Minimize(forceTol=1e-8) == 0


def test_gradient_matches_finite_differences():
    m = benzene()
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    pos = list(ff.Positions())
    grad = ff.CalcGrad()
    for i in (0, 4, 13, 30):
        h = 1e-5
        plus = pos[:i] + [pos[i] + h] + pos[i + 1 :]
        minus = pos[:i] + [pos[i] - h] + pos[i + 1 :]
        numeric = (ff.CalcEnergy(plus) - ff.CalcEnergy(minus)) / (2 * h)
        assert grad[i] == pytest.approx(numeric, abs=1e-4)
    with pytest.raises(ValueError):
        ff.CalcEnergy([0.0, 1.0])


def test_fixed_points_do_not_move():
    m = benzene()
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    start = m.GetConformer().GetPositions()[0].copy()
    ff.AddFixedPoint(0)
    ff.Minimize()
    assert abs(m.GetConformer().GetPositions()[0] - start).max() < 1e-12


def test_all_conformers():
    m = benzene()
    conf = Chem.Conformer(m.GetConformer())
    conf.SetId(1)
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x * 1.1, p.y * 1.1, p.z))
    m.AddConformer(conf)
    res = ChemicalForceFields.UFFOptimizeMoleculeConfs(m, numThreads=1)
    res_mt = ChemicalForceFields.UFFOptimizeMoleculeConfs(benzene_two(), numThreads=0)
    assert len(res) == 2 and all(r[0] == 0 for r in res)
    assert res == pytest.approx(res_mt)


def benzene_two():
    m = benzene()
    conf = Chem.Conformer(m.GetConformer())
    conf.SetId(1)
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x * 1.1, p.y * 1.1, p.z))
    m.AddConformer(conf)
    return m


def test_parameters():
    m = benzene()
    assert ChemicalForceFields.UFFHasAllMoleculeParams(m)
    kb, r0 = ChemicalForceFields.GetUFFBondStretchParams(m, 0, 1)
    assert 1.3 < r0 < 1.5 and kb > 0
    assert ChemicalForceFields.GetUFFBondStretchParams(m, 0, 3) is None
    ka, theta0 = ChemicalForceFields.GetUFFAngleBendParams(m, 0, 1, 2)
    assert theta0 == pytest.approx(120.0) and ka > 0
    assert ChemicalForceFields.GetUFFTorsionParams(m, 0, 1, 2, 3) > 0
    assert ChemicalForceFields.GetUFFInversionParams(m, 0, 1, 2, 3) is None
    assert ChemicalForceFields.GetUFFInversionParams(m, 1, 0, 5, 6) is not None
    x, d = ChemicalForceFields.GetUFFVdWParams(m, 0, 3)
    assert x > 0 and d > 0


def test_module_layout():
    assert rdForceField.ForceField is ChemicalForceFields.ForceField


def test_distance_constraint_pulls_atoms():
    m = benzene()
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    ff.UFFAddDistanceConstraint(0, 3, False, 4.0, 4.0, 500.0)
    assert ff.Minimize(maxIts=500) in (0, 1)
    pos = m.GetConformer().GetPositions()
    assert ((pos[0] - pos[3]) ** 2).sum() ** 0.5 > 3.0


def test_constraint_errors():
    ff = ChemicalForceFields.UFFGetMoleculeForceField(benzene())
    with pytest.raises(RuntimeError):
        ff.UFFAddDistanceConstraint(0, 99, False, 1.0, 2.0, 10.0)
    with pytest.raises(RuntimeError):
        ff.UFFAddDistanceConstraint(0, 1, False, 3.0, 2.0, 10.0)
    with pytest.raises(RuntimeError):
        ff.UFFAddAngleConstraint(0, 1, 2, False, 10.0, 200.0, 10.0)


def test_extra_points():
    m = benzene()
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    assert ff.NumPoints() == 12
    ff.AddExtraPoint(5.0, 0.0, 0.0)
    assert ff.NumPoints() == 13
    assert ff.GetExtraPointPos(0) == (5.0, 0.0, 0.0)
    with pytest.raises(IndexError):
        ff.GetExtraPointPos(1)
    # a fixed extra point can anchor a position-like distance constraint
    ff.UFFAddDistanceConstraint(0, 12, False, 1.0, 1.0, 100.0)
    ff.Minimize(maxIts=500)
    assert ff.GetExtraPointPos(0) == (5.0, 0.0, 0.0)
    assert len(ff.Positions()) == 39
