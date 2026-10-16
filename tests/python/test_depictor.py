"""rdDepictor.Compute2DCoords and MolToMolBlock without coordinates.
tests/data/depict_expected.json holds RDKit 2026.09.1 coordinates (rounded
to 4 decimals) for molecules that come out the same on arm64 and x86-64
builds; exact bits differ between platforms (fused multiply-add)."""
import json
import os

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdDepictor
from rdoxide.Geometry import Point2D

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "depict_expected.json")
CASES = json.load(open(DATA))


def coords(m):
    c = m.GetConformer()
    return [[round(c.GetAtomPosition(i).x, 4), round(c.GetAtomPosition(i).y, 4)]
            for i in range(m.GetNumAtoms())]


@pytest.mark.parametrize("case", CASES, ids=[c["smiles"] for c in CASES])
def test_matches_rdkit(case):
    m = Chem.MolFromSmiles(case["smiles"])
    for key, kw in (("canon", {}), ("nocanon", {"canonOrient": False})):
        mm = Chem.Mol(m)
        assert rdDepictor.Compute2DCoords(mm, **kw) == 0
        assert not mm.GetConformer().Is3D()
        assert coords(mm) == case[key]
    assert Chem.MolToMolBlock(m) == case["block"]


def test_options():
    m = Chem.MolFromSmiles("c1ccccc1O")
    rdDepictor.Compute2DCoords(m, canonOrient=False)
    cid = rdDepictor.Compute2DCoords(m, clearConfs=False)
    assert cid == 1 and m.GetNumConformers() == 2
    rdDepictor.Compute2DCoords(m, bondLength=2.0)
    p0, p1 = m.GetConformer().GetAtomPosition(0), m.GetConformer().GetAtomPosition(1)
    assert abs((p0 - p1).Length() - 2.0) < 1e-6
    # the default bond length is restored afterwards
    rdDepictor.Compute2DCoords(m)
    p0, p1 = m.GetConformer().GetAtomPosition(0), m.GetConformer().GetAtomPosition(1)
    assert abs((p0 - p1).Length() - 1.5) < 1e-6
    # a single fixed atom translates the depiction
    rdDepictor.Compute2DCoords(m, coordMap={0: Point2D(10.0, 10.0)})
    p0 = m.GetConformer().GetAtomPosition(0)
    assert (round(p0.x, 6), round(p0.y, 6)) == (10.0, 10.0)
    with pytest.raises(ValueError):
        rdDepictor.Compute2DCoords(m, coordMap={99: Point2D(0.0, 0.0)})
