"""PMI, PBF and the derived shape descriptors (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdoxide import Chem
from rdoxide.Chem import Descriptors3D
from rdoxide.Chem import rdMolDescriptors as D

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

NAMES = ["PBF", "NPR1", "NPR2", "PMI1", "PMI2", "PMI3", "RadiusOfGyration", "InertialShapeFactor",
         "Eccentricity", "Asphericity", "SpherocityIndex"]
EXPECTED = {
    "ethanol": [0.4823414930196825, 0.2540480685156133, 0.8728777586729661, 15.11773660225969,
                51.942674150965665, 59.507386498120866, 1.1720394242631202, 0.05773865371767985,
                0.9671915936790839, 0.422455319681518, 0.4607287308599609],
    "caffeine": [0.2516143275290736, 0.39626054804626437, 0.6121120609367914, 458.7295241630343,
                 708.6092113697001, 1157.6462164219197, 2.446679998406103, 0.001334363777987929,
                 0.9181381040247019, 0.27843898691979607, 0.06699208320154518],
}


def load(name):
    return Chem.MolFromMolBlock(BLOCKS[name], removeHs=False)


@pytest.mark.parametrize("name", EXPECTED)
def test_descriptors(name):
    m = load(name)
    values = [getattr(D, "Calc" + n)(m) for n in NAMES]
    assert values == pytest.approx(EXPECTED[name], rel=1e-12)


def test_mass_weighting_and_caching():
    m = load("ethanol")
    weighted = D.CalcPMI1(m)
    unweighted = D.CalcPMI1(m, useAtomicMasses=False)
    assert weighted != unweighted
    assert m.HasProp("_PMI1_mass")  # the moments are cached on the molecule
    # without force the cached value is returned even if the coordinates change
    conf = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (2 * p.x, 2 * p.y, 2 * p.z))
    assert D.CalcPMI1(m, force=False) == pytest.approx(weighted, rel=1e-12)
    assert D.CalcPMI1(m) == pytest.approx(4 * weighted, rel=1e-9)


def test_degenerate_cases():
    assert D.CalcPBF(load("water")) == 0.0  # fewer than four atoms
    assert D.CalcNPR1(load("single")) == 0.0
    assert D.CalcEccentricity(load("single")) == 0.0


def test_errors_and_metadata():
    m = Chem.MolFromSmiles("CCO")
    for name in NAMES:
        with pytest.raises(RuntimeError, match="molecule has no conformers"):
            getattr(D, "Calc" + name)(m)
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        D.CalcPMI1(load("ethanol"), confId=5)
    assert all(getattr(D, f"_Calc{n}_version") == "1.0.0" for n in NAMES)
    assert Descriptors3D.NPR1(load("caffeine")) == pytest.approx(0.39626054804626437, rel=1e-12)
