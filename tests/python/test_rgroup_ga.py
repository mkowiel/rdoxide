"""GA matching strategy of RGroupDecomposition: with a fixed seed the rows equal
those of the official RDKit 2026.09.1 wheel (`tests/data/rgroup_ga.json`)."""

import json
import random
from pathlib import Path

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdRGroupDecomposition as R

DATA = json.loads((Path(__file__).parent.parent / "data" / "rgroup_ga.json").read_text())

SUBS = ["Cl", "F", "Br", "C", "CC", "O", "N", "OC", "C(=O)O", "CN", "C#N", "[N+](=O)[O-]", "S", "CCC"]


def library():
    rnd = random.Random(5)
    smis = []
    for i in range(24):
        s = rnd.sample(SUBS, 3)
        smis.append(("c1(%s)c(%s)c(%s)ccc1" if i % 2 else "c1(%s)cc(%s)cc(%s)c1") % tuple(s))
    return smis


def decompose(method, seed):
    params = R.RGroupDecompositionParameters()
    params.matchingStrategy = R.RGroupMatching.GA
    params.scoreMethod = getattr(R.RGroupScore, method)
    params.gaRandomSeed = seed
    params.gaNumberOperationsWithoutImprovement = 2000
    d = R.RGroupDecomposition(Chem.MolFromSmiles("c1ccccc1"), params)
    for s in library():
        d.Add(Chem.MolFromSmiles(s))
    ok = d.Process()
    rows = [{k: Chem.MolToSmiles(v) for k, v in r.items()} for r in d.GetRGroupsAsRows()]
    return ok, rows


@pytest.mark.parametrize("case", DATA, ids=lambda c: f"{c['method']}-{c['seed']}")
def test_rows_equal_official_wheel(case):
    ok, rows = decompose(case["method"], case["seed"])
    assert ok
    assert rows == case["rows"]


def test_small_problems_fall_back_to_exhaustive():
    params = R.RGroupDecompositionParameters()
    params.matchingStrategy = R.RGroupMatching.GA
    d = R.RGroupDecomposition(Chem.MolFromSmiles("c1ccccc1"), params)
    for s in ["c1ccccc1Cl", "c1ccccc1Cl", "Clc1ccccc1F"]:
        d.Add(Chem.MolFromSmiles(s))
    assert d.Process()
    assert len(d.GetRGroupsAsRows()) == 3
