"""AUTOCORR3D, RDF and MORSE (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

NAMES = ["AUTOCORR3D", "RDF", "MORSE"]
LENGTHS = {"AUTOCORR3D": 80, "RDF": 210, "MORSE": 224}
# (sum of all values, first six values, last three values)
EXPECTED = {
    ("ethanol", "AUTOCORR3D"): (14.906, [0.262, 0.737, 0.905, 0.287, 0.0, 0.0], [0.0, 0.0, 0.0]),
    ("ethanol", "RDF"): (82.799, [2.392, 1.19, 2.534, 0.103, 4.522, 0.432], [0.0, 0.0, 0.0]),
    ("ethanol", "MORSE"): (319.355, [36.0, 13.636, -1.214, -0.275, -0.53, -1.042], [-0.069, 0.27, -0.388]),
    ("caffeine", "AUTOCORR3D"): (34.227, [0.114, 0.339, 0.526, 0.811, 0.763, 0.841], [0.0, 0.0, 0.0]),
    ("caffeine", "RDF"): (903.526, [3.131, 3.891, 3.661, 12.898, 8.053, 10.012], [0.0, 0.0, 0.0]),
    ("caffeine", "MORSE"): (2461.3, [276.0, 20.732, -2.205, -1.384, -3.223, -1.633], [0.43, 0.601, 1.036]),
}


def load(name):
    return Chem.MolFromMolBlock(BLOCKS[name], removeHs=False)


@pytest.mark.parametrize("name,desc", EXPECTED)
def test_values(name, desc):
    total, head, tail = EXPECTED[name, desc]
    v = getattr(D, "Calc" + desc)(load(name))
    assert len(v) == LENGTHS[desc]
    assert v[:6] == pytest.approx(head, abs=1e-9)
    assert v[-3:] == pytest.approx(tail, abs=1e-9)
    assert sum(v) == pytest.approx(total, abs=1e-6)


@pytest.mark.parametrize("desc", NAMES)
def test_custom_atom_property(desc):
    m = load("ethanol")
    ones = getattr(D, "Calc" + desc)(m, CustomAtomProperty="missing")
    assert len(ones) == LENGTHS[desc] // (8 if desc == "AUTOCORR3D" else 7)
    # a missing property weighs every atom 1, which is the unweighted profile
    if desc != "AUTOCORR3D":
        assert ones == getattr(D, "Calc" + desc)(m)[: len(ones)]
    for a in m.GetAtoms():
        a.SetDoubleProp("w", 2.0)
    assert getattr(D, "Calc" + desc)(m, CustomAtomProperty="w") == pytest.approx(
        [4 * x for x in ones], abs=4e-3)


@pytest.mark.parametrize("desc", NAMES)
def test_distance_matrix_is_cached(desc):
    m = load("ethanol")
    assert not m.HasProp("_3DDistanceMatrix_Conf0")
    getattr(D, "Calc" + desc)(m)
    assert m.HasProp("_3DDistanceMatrix_Conf0")


@pytest.mark.parametrize("desc", NAMES)
def test_errors_and_metadata(desc):
    fn = getattr(D, "Calc" + desc)
    with pytest.raises(RuntimeError, match="molecule has no conformers"):
        fn(Chem.MolFromSmiles("CCO"))
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        fn(load("ethanol"), confId=5)
    assert getattr(D, f"_Calc{desc}_version") == "1.0.0"
