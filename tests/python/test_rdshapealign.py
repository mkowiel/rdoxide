"""rdShapeAlign against values recorded from the official RDKit wheel (2026.09.1)."""
from pathlib import Path

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdShapeAlign as S

DATA = Path(__file__).resolve().parents[2] / "upstream/rdkit/External/pubchem_shape/test_data/test1.sdf"

pytestmark = pytest.mark.skipif(not DATA.exists(), reason="upstream sources not fetched")


def mols():
    return list(Chem.SDMolSupplier(str(DATA)))


@pytest.mark.parametrize(
    "kw,expected",
    [
        ({}, (0.7728198701627574, 0.30374482395256774)),
        ({"opt_param": 0.5}, (0.7728091089903107, 0.304772109158682)),
    ],
)
def test_align_mol_matches_wheel(kw, expected):
    ref, probe = mols()[:2]
    st, ct = S.AlignMol(ref, probe, **kw)
    assert st == pytest.approx(expected[0], abs=1e-9)
    assert ct == pytest.approx(expected[1], abs=1e-9)


def test_shapes_and_transform_roundtrip():
    ref, probe = mols()[:2]
    a, b = S.PrepareConformer(ref), S.PrepareConformer(probe)
    st, ct, matrix = S.AlignShapes(a, b)
    assert st == pytest.approx(0.7728198701627574, abs=1e-9)
    assert len(matrix) == 12
    with pytest.raises(ValueError):
        S.TransformConformer([0.0] * 3, [0.0] * 9, matrix[:5], b, probe.GetConformer())


def test_options_reject_unknown_attribute():
    o = S.ShapeInputOptions()
    o.useColors = False
    assert o.useColors is False
    with pytest.raises(AttributeError):
        o.notAnOption = 1
