"""USR, USRCAT and CoulombMat (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolDescriptors as D
from rdoxide.Geometry import Point3D

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

CAFFEINE_USR = [2.8355182336516407, 0.9616021997531262, -0.8539666879802394, 2.893547228247012,
                1.1214230847398174, -0.8477425646803405, 4.645130057318242, 2.2479306093268625,
                -0.7235621092777416, 4.317018016080865, 2.03911001094838, -0.6729159953921767]
CAFFEINE_SEL = [2.746080236328757, 0.4800474445976515, -0.4207502312099483, 2.0114331392175457,
                0.4941022361726196, -0.7613023729941901, 1.9261955424016197, 0.6096637781860845,
                -0.7406984074482842, 6.166643657245026, 0.36612819871326374, -0.8791884230686698]


def load(name):
    return Chem.MolFromMolBlock(BLOCKS[name], removeHs=False)


def test_usr():
    assert D.GetUSR(load("caffeine")) == pytest.approx(CAFFEINE_USR, rel=1e-12)


def test_usrcat_default_and_selection():
    m = load("caffeine")
    cat = D.GetUSRCAT(m)
    assert len(cat) == 60
    assert cat[:12] == pytest.approx(CAFFEINE_USR, rel=1e-12)
    # atom selections are 1-based
    sel = D.GetUSRCAT(m, [[1, 2, 3]])
    assert len(sel) == 24
    assert sel[12:] == pytest.approx(CAFFEINE_SEL, rel=1e-12)
    assert D.GetUSRCAT(m, [[]])[12:] == [0.0] * 12


def test_scores():
    c, a = load("caffeine"), load("aspirin")
    assert D.GetUSRScore(D.GetUSR(c), D.GetUSR(a)) == pytest.approx(0.9043939748094802, rel=1e-12)
    cc, ca = D.GetUSRCAT(c), D.GetUSRCAT(a)
    assert D.GetUSRScore(cc, ca) == pytest.approx(0.1704297834069872, rel=1e-12)
    assert D.GetUSRScore(cc, cc) == 1.0
    weighted = D.GetUSRScore(cc, ca, [1, 2, 3, 4, 5])
    assert weighted < D.GetUSRScore(cc, ca)
    with pytest.raises(ValueError, match="number of weights"):
        D.GetUSRScore(cc, ca, [1, 2])
    with pytest.raises(ValueError, match="same length"):
        D.GetUSRScore(cc, D.GetUSR(a))


def test_distribution_helpers():
    pts = [Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(0, 2, 0), Point3D(0, 0, 3), Point3D(1, 1, 1)]
    ref = []
    dist = D.GetUSRDistributions(pts, ref)
    assert len(dist) == 4 and len(ref) == 4
    assert all(len(d) == len(pts) for d in dist)
    assert [list(d) for d in D.GetUSRDistributionsFromPoints(pts, ref)] == [list(d) for d in dist]
    usr = D.GetUSRFromDistributions(dist)
    assert len(usr) == 12
    with pytest.raises(ValueError, match="no coordinates"):
        D.GetUSRDistributions([])
    with pytest.raises(ValueError, match="no points"):
        D.GetUSRDistributionsFromPoints(pts, [])
    with pytest.raises(ValueError, match="distances missing"):
        D.GetUSRFromDistributions([[1.0], []])
    with pytest.raises(RuntimeError, match="3 times more elements"):
        D.GetUSRFromDistributions(dist[:3])
    with pytest.raises(TypeError):
        D.GetUSRDistributions(pts, (1,))


def test_errors():
    flat = Chem.MolFromSmiles("CCO")
    with pytest.raises(ValueError, match="no conformers"):
        D.GetUSR(flat)
    with pytest.raises(ValueError, match="too few atoms"):
        D.GetUSR(load("single"))
    with pytest.raises(ValueError, match="empty atom selections"):
        D.GetUSRCAT(load("caffeine"), [])
    with pytest.raises(OverflowError):
        D.GetUSRCAT(load("caffeine"), [[-1]])


def test_coulomb_matrix():
    m = D.CalcCoulombMat(load("ethanol"))
    assert len(m) == 9 and all(len(r) == 9 for r in m)
    assert list(m[0][:3]) == pytest.approx([36.85810519942594, 23.63585207528015, 20.584422909870824],
                                           rel=1e-12)
    assert m[1][2] == m[2][1]
    assert D._CalcCoulombMat_version == "1.0.0"
    with pytest.raises(RuntimeError, match="molecule has no conformers"):
        D.CalcCoulombMat(Chem.MolFromSmiles("CCO"))
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        D.CalcCoulombMat(load("ethanol"), confId=4)
