"""rdSynthonSpaceSearch: API surface not covered by the upstream tests."""
import os
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdSynthonSpaceSearch

DATA = Path(__file__).resolve().parents[2] / "upstream" / "rdkit" / "Code" / "GraphMol" / "SynthonSpaceSearch" / "data"
pytestmark = pytest.mark.skipif(not DATA.exists(), reason="upstream test data not checked out")


def test_text_round_trip(tmp_path):
    space = rdSynthonSpaceSearch.SynthonSpace()
    space.ReadTextFile(str(DATA / "amide_space.txt"))
    assert space.GetNumReactions() == 1
    assert space.GetNumProducts() == 12
    assert space.GetNumSynthons() == 7
    out = tmp_path / "amide.spc"
    space.WriteDBFile(str(out))
    space2 = rdSynthonSpaceSearch.SynthonSpace()
    space2.ReadDBFile(str(out))
    assert space2.GetNumProducts() == 12
    q = Chem.MolFromSmarts("c1ccccc1C(=O)N1CCCC1")
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = -1
    r1 = space.SubstructureSearch(q, params=params)
    r2 = space2.SubstructureSearch(q, params=params)
    n1 = sorted(m.GetProp("_Name") for m in r1.GetHitMolecules())
    n2 = sorted(m.GetProp("_Name") for m in r2.GetHitMolecules())
    assert n1 and n1 == n2


def test_params_are_checked():
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    with pytest.raises(AttributeError):
        params.noSuchOption = 1
    params.maxHits = 5
    assert params.maxHits == 5
    assert params.similarityCutoff == 0.5


def test_enumerated_file(tmp_path):
    space = rdSynthonSpaceSearch.SynthonSpace()
    space.ReadTextFile(str(DATA / "amide_space.txt"))
    out = tmp_path / "enum.smi"
    space.WriteEnumeratedFile(str(out))
    assert len(out.read_text().splitlines()) == 12


def test_formatted_integer_string():
    assert rdSynthonSpaceSearch.FormattedIntegerString(1234567) == "1 234 567"
    assert rdSynthonSpaceSearch.FormattedIntegerString(12) == "12"
