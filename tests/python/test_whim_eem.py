"""WHIM and EEM charges (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolDescriptors as D

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

# (first six, last four, sum) of the 114 values
WHIM = {
    "ethanol": ([1.405, 0.599, 0.363, 0.593, 0.253, 0.24], [4.137, 3.06, 4.436, 3.394], 86.394),
    "caffeine": ([5.258, 3.506, 0.2, 0.587, 0.391, 0.179], [30.83, 23.323, 33.947, 24.594], 424.544),
}
CUSTOM = {
    "ethanol": [1.581, 0.763, 0.396, 0.577, 0.278, 0.24, 0.24, 0.24, 0.661, 0.761, 0.504, 2.74, 2.133,
                0.24, 0.366, 0.642, 5.35],
    "caffeine": [5.95, 4.475, 0.303, 0.555, 0.417, 0.179, 0.179, 0.196, 0.741, 0.824, 0.56, 10.729,
                 29.791, 0.185, 0.458, 0.708, 48.6],
}
EEM = {
    "ethanol": [-0.5309430092431018, -0.013404997902780557, -0.7234167614413173, 0.19630285976213874],
    "caffeine": [-0.3003101718891141, -0.3653604924530798, 0.10713484840414692, -0.46425225211653404],
}


def load(name):
    return Chem.MolFromMolBlock(BLOCKS[name], removeHs=False)


@pytest.mark.parametrize("name", WHIM)
def test_whim(name):
    head, tail, total = WHIM[name]
    v = D.CalcWHIM(load(name))
    assert len(v) == 114
    assert v[:6] == pytest.approx(head, abs=1e-9)
    assert v[-4:] == pytest.approx(tail, abs=1e-9)
    assert sum(v) == pytest.approx(total, abs=1e-6)


@pytest.mark.parametrize("name", CUSTOM)
def test_whim_custom_property(name):
    m = load(name)
    for a in m.GetAtoms():
        a.SetDoubleProp("w", 1.0 + 0.5 * a.GetIdx())
    assert D.CalcWHIM(m, CustomAtomProperty="w") == pytest.approx(CUSTOM[name], abs=1e-9)


def test_whim_threshold_changes_symmetry_terms():
    m = Chem.MolFromMolBlock(BLOCKS["benzene"], removeHs=False)
    assert D.CalcWHIM(m, thresh=0.5) != D.CalcWHIM(m, thresh=1e-6)


@pytest.mark.parametrize("name", EEM)
def test_eem_charges(name):
    m = load(name)
    q = D.CalcEEMcharges(m)
    assert len(q) == m.GetNumAtoms()
    assert q[:4] == pytest.approx(EEM[name], rel=1e-9)
    assert sum(q) == pytest.approx(0.0, abs=1e-9)  # the total charge is conserved
    assert not m.GetAtomWithIdx(0).HasProp("_EEMCharge")  # upstream works on a copy


@pytest.mark.parametrize("fn", ["CalcWHIM", "CalcEEMcharges"])
def test_errors_and_metadata(fn):
    f = getattr(D, fn)
    with pytest.raises(RuntimeError, match="molecule has no conformers"):
        f(Chem.MolFromSmiles("CCO"))
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        f(load("ethanol"), confId=5)
    assert getattr(D, "_Calc%s_version" % ("WHIM" if fn == "CalcWHIM" else "EMMcharges")) == "1.0.0"
