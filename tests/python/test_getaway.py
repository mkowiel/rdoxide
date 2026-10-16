"""GETAWAY descriptors (expectations from RDKit 2026.09.1)."""

import json
import math
import os

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolDescriptors as D

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

# (first six, last four, sum) of the 273 values
EXPECTED = {
    "ethanol": ([4.755, 1.0, 2.958, 26.69, 3.0, 0.81], [0.0, 0.0, 0.0, 2.051], 314.413),
    "caffeine": ([44.548, 0.836, 4.112, 8.069, 3.0, 1.524], [0.183, 0.161, 0.0, 0.802], 538.94),
}
# (first six, last three, sum) of the 45 values with the weights 1 + idx / 2
CUSTOM = {
    "ethanol": ([4.755, 1.0, 2.958, 26.69, 40.334, 4.979], [0.0, 0.0, 4.813], 510.572),
    "caffeine": ([44.548, 0.836, 4.112, 8.069, 262.481, 63.407], [4.452, 0.0, 23.097], 4322.901),
}


def load(name):
    return Chem.MolFromMolBlock(BLOCKS[name], removeHs=False)


@pytest.mark.parametrize("name", EXPECTED)
def test_getaway(name):
    head, tail, total = EXPECTED[name]
    v = D.CalcGETAWAY(load(name))
    assert len(v) == 273
    assert v[:6] == pytest.approx(head, abs=1e-9)
    assert v[-4:] == pytest.approx(tail, abs=1e-9)
    assert sum(v) == pytest.approx(total, abs=1e-6)


@pytest.mark.parametrize("name", CUSTOM)
def test_custom_atom_property(name):
    head, tail, total = CUSTOM[name]
    m = load(name)
    for a in m.GetAtoms():
        a.SetDoubleProp("w", 1.0 + 0.5 * a.GetIdx())
    v = D.CalcGETAWAY(m, CustomAtomProperty="w")
    assert len(v) == 45
    assert v[:6] == pytest.approx(head, abs=1e-9)
    assert v[-3:] == pytest.approx(tail, abs=1e-9)
    assert sum(v) == pytest.approx(total, abs=1e-6)


def test_precision_changes_the_information_indices_only():
    m = load("caffeine")
    coarse, fine = D.CalcGETAWAY(m, precision=2), D.CalcGETAWAY(m, precision=3)
    assert coarse[:2] != fine[:2]
    assert coarse[2:] == fine[2:]


def test_single_atom_is_degenerate():
    v = D.CalcGETAWAY(load("single"))
    assert v[0] == 0.0 and math.isnan(v[1]) and math.isnan(v[2])


def test_caches_and_errors():
    m = load("ethanol")
    D.CalcGETAWAY(m)
    assert m.HasProp("_3DDistanceMatrix_Conf0") and m.HasProp("AdjacencyMatrix")
    with pytest.raises(RuntimeError, match="molecule has no conformers"):
        D.CalcGETAWAY(Chem.MolFromSmiles("CCO"))
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        D.CalcGETAWAY(m, confId=5)
    assert D._CalcGETAWAY_version == "1.0.0"
